# Validação externa com sinais reais do VitalDB.
# Este script é aditivo: não altera src/ nem o script antigo scripts/07_validacao_vitaldb.jl.

using Pkg
const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
Pkg.activate(ROOT)

using CSV
using DataFrames
using Statistics
using Printf
using Plots
using Optim
using PropofolWindkessel

const DATA_DIR = joinpath(ROOT, "dados", "vitaldb", "processados")
const CASES_DIR = joinpath(DATA_DIR, "casos")
const OUT_DIR = joinpath(ROOT, "resultados", "vitaldb")
const PK_MODEL = :schnider
const WINDOW_SECONDS = 1200.0
const BASELINE_SECONDS = 60.0
const MAX_INTERPOLATION_GAP_S = 5.0

struct PreparedCase
    caseid::Int
    partition::String
    coadmin_status::String
    patient::Paciente
    pv::PacienteVirtual
    pkpd::NamedTuple
    tfim_min::Float64
    t_pam::Vector{Float64}
    y_pam::Vector{Float64}
    t_hr::Vector{Float64}
    y_hr::Vector{Float64}
    t_co::Vector{Float64}
    y_co::Vector{Float64}
    t_pas::Vector{Float64}
    y_pas::Vector{Float64}
    t_pad::Vector{Float64}
    y_pad::Vector{Float64}
end

function numeric_column(df::DataFrame, name::Symbol)
    return [ismissing(v) ? NaN : Float64(v) for v in df[!, name]]
end

function valid_series(t::Vector{Float64}, y::Vector{Float64})
    keep = [isfinite(ti) && isfinite(yi) for (ti, yi) in zip(t, y)]
    return t[keep], y[keep]
end

function median_baseline(df::DataFrame, name::Symbol; minimum_n::Int=5)
    t = numeric_column(df, :t_rel_s)
    y = numeric_column(df, name)
    vals = [y[i] for i in eachindex(y) if t[i] < BASELINE_SECONDS && isfinite(y[i])]
    length(vals) >= minimum_n || error("Baseline insuficiente em $name: $(length(vals)) amostras válidas")
    return median(vals)
end

function dose_events_from_observed_rate(df::DataFrame)
    t_s = numeric_column(df, :t_rel_s)
    rate = numeric_column(df, :propofol_rate_mg_min)
    all(isfinite, rate) || error("Taxa de propofol contém lacunas; caso não será simulada")
    length(rate) == length(t_s) || error("Vetores de tempo e taxa têm comprimentos diferentes")
    events = typeof(infusao(0.0, 1 / 60, 1.0))[]
    i = 1
    while i <= length(rate)
        r = rate[i]
        j = i
        while j < length(rate) && isapprox(rate[j + 1], r; atol=1e-8, rtol=1e-6)
            j += 1
        end
        if r > 0.0
            start_min = t_s[i] / 60.0
            duration_min = (j - i + 1) / 60.0
            push!(events, infusao(start_min, duration_min, r))
        end
        i = j + 1
    end
    isempty(events) && error("Nenhum evento positivo de propofol foi encontrado")
    return events
end

function make_observed_series(df::DataFrame, col::Symbol)
    t = numeric_column(df, :t_rel_s)
    y = numeric_column(df, col)
    return valid_series(t, y)
end

function prepare_case(row::DataFrameRow, partition::String)
    cid = Int(row.caseid)
    csvpath = joinpath(CASES_DIR, @sprintf("case_%05d.csv", cid))
    isfile(csvpath) || error("Arquivo de sinais não encontrado: $csvpath")
    df = CSV.read(csvpath, DataFrame; missingstring=["", "NA", "NaN"])
    nrow(df) >= 1261 || error("Janela incompleta")

    pam0 = median_baseline(df, :pam_obs_mmHg; minimum_n=10)
    pas0 = median_baseline(df, :pas_obs_mmHg; minimum_n=10)
    pad0 = median_baseline(df, :pad_obs_mmHg; minimum_n=10)
    hr0  = median_baseline(df, :hr_obs_bpm; minimum_n=10)
    co0  = median_baseline(df, :co_obs_l_min; minimum_n=3)
    pad0 < pam0 < pas0 || error("Baseline observado não satisfaz PAD < PAM < PAS")
    20 <= hr0 <= 250 || error("FC basal fora dos limites de qualidade")
    0.5 <= co0 <= 20 || error("DC basal fora dos limites de qualidade")

    sex_text = uppercase(strip(String(row.sex)))
    sex_text in ("M", "F") || error("Sexo ausente ou não reconhecido")
    p = Paciente(
        nome="VitalDB_case_$cid",
        idade=Float64(row.age),
        peso=Float64(row.weight_kg),
        altura=Float64(row.height_cm),
        sexo=Symbol(sex_text),
        asa=Int(round(Float64(row.asa))),
        pas=pas0,
        pad=pad0,
        pam=pam0,
        fc=hr0,
        dc=co0,
    )

    doses = dose_events_from_observed_rate(df)
    tfim = maximum(numeric_column(df, :t_rel_s)) / 60.0
    tfim >= 21.0 || error("Tempo total menor que 21 minutos")
    # s=1 é declarado explicitamente como paciente de referência: PVI e IVCI não
    # estão garantidos no dataset, então não usamos os valores padrão como se fossem medidos.
    pv = calibrar(p; modelo_pk=PK_MODEL, s=1.0, ajustar_C=true)
    pkpd = simular_propofol(p, doses; modelo=PK_MODEL, tfim=tfim)

    t_pam, y_pam = make_observed_series(df, :pam_obs_mmHg)
    t_hr, y_hr = make_observed_series(df, :hr_obs_bpm)
    t_co, y_co = make_observed_series(df, :co_obs_l_min)
    t_pas, y_pas = make_observed_series(df, :pas_obs_mmHg)
    t_pad, y_pad = make_observed_series(df, :pad_obs_mmHg)
    status = ismissing(row.coadmin_status) ? "unknown" : String(row.coadmin_status)
    return PreparedCase(cid, partition, status, p, pv, pkpd, tfim,
                        t_pam, y_pam, t_hr, y_hr, t_co, y_co, t_pas, y_pas, t_pad, y_pad)
end

function interpolate_observation(t::Vector{Float64}, y::Vector{Float64}, target::Float64;
                                 max_gap::Float64=MAX_INTERPOLATION_GAP_S)
    isempty(t) && return NaN
    idx = searchsortedlast(t, target)
    if idx >= 1 && isapprox(t[idx], target; atol=1e-8)
        return y[idx]
    end
    if idx < 1 || idx >= length(t)
        return NaN
    end
    gap = t[idx + 1] - t[idx]
    (gap > 0 && gap <= max_gap) || return NaN
    fraction = (target - t[idx]) / gap
    return y[idx] + fraction * (y[idx + 1] - y[idx])
end

function aligned_pair(t_pred_s::Vector{Float64}, y_pred::Vector{Float64},
                      t_obs::Vector{Float64}, y_obs::Vector{Float64})
    predicted = Float64[]
    observed = Float64[]
    times = Float64[]
    for i in eachindex(t_pred_s)
        t = t_pred_s[i]
        BASELINE_SECONDS <= t <= BASELINE_SECONDS + WINDOW_SECONDS || continue
        yp = y_pred[i]
        yo = interpolate_observation(t_obs, y_obs, t)
        if isfinite(yp) && isfinite(yo)
            push!(times, t)
            push!(predicted, yp)
            push!(observed, yo)
        end
    end
    return times, predicted, observed
end

function predict_case(c::PreparedCase, effect::EfeitoHemodinamico)
    bat, _ = simular_hemodinamica(c.patient, c.pkpd;
        efeito=effect, θbasal=c.pv.θ, s=1.0, tfim=c.tfim_min)
    t_pred_s = Float64.(bat.t_min .* 60.0)
    return (bat=bat, t_s=t_pred_s)
end

function metric_summary(pred::Vector{Float64}, obs::Vector{Float64})
    keep = [isfinite(a) && isfinite(b) for (a, b) in zip(pred, obs)]
    p, o = pred[keep], obs[keep]
    if isempty(p)
        return (n=0, mae=NaN, rmse=NaN, bias=NaN, correlation=NaN)
    end
    e = p .- o
    r = length(p) >= 2 && std(p) > 0 && std(o) > 0 ? cor(p, o) : NaN
    return (n=length(e), mae=mean(abs.(e)), rmse=sqrt(mean(e .^ 2)), bias=mean(e), correlation=r)
end

function case_metrics(c::PreparedCase, effect::EfeitoHemodinamico, method::String)
    prediction = predict_case(c, effect)
    bat, ts = prediction.bat, prediction.t_s
    _, ppam, opam = aligned_pair(ts, Float64.(bat.PAM), c.t_pam, c.y_pam)
    _, phr, ohr = aligned_pair(ts, Float64.(bat.FC), c.t_hr, c.y_hr)
    _, pco, oco = aligned_pair(ts, Float64.(bat.DC), c.t_co, c.y_co)
    mp, mh, mc = metric_summary(ppam, opam), metric_summary(phr, ohr), metric_summary(pco, oco)
    return (caseid=c.caseid, partition=c.partition, coadmin_status=c.coadmin_status, method=method,
        n_pam=mp.n, pam_mae_mmHg=mp.mae, pam_rmse_mmHg=mp.rmse, pam_bias_mmHg=mp.bias, pam_correlation=mp.correlation,
        n_hr=mh.n, hr_mae_bpm=mh.mae, hr_rmse_bpm=mh.rmse, hr_bias_bpm=mh.bias,
        n_co=mc.n, co_mae_l_min=mc.mae, co_rmse_l_min=mc.rmse, co_bias_l_min=mc.bias)
end

function prediction_dataframe(c::PreparedCase, effect::EfeitoHemodinamico)
    prediction = predict_case(c, effect)
    bat, ts = prediction.bat, prediction.t_s
    keep = findall(t -> BASELINE_SECONDS <= t <= BASELINE_SECONDS + WINDOW_SECONDS, ts)
    t = ts[keep]
    obs_pam = [interpolate_observation(c.t_pam, c.y_pam, ti) for ti in t]
    obs_pas = [interpolate_observation(c.t_pas, c.y_pas, ti) for ti in t]
    obs_pad = [interpolate_observation(c.t_pad, c.y_pad, ti) for ti in t]
    obs_hr  = [interpolate_observation(c.t_hr, c.y_hr, ti) for ti in t]
    obs_co  = [interpolate_observation(c.t_co, c.y_co, ti) for ti in t]
    pred_pam, pred_hr, pred_co = Float64.(bat.PAM[keep]), Float64.(bat.FC[keep]), Float64.(bat.DC[keep])
    return DataFrame(
        caseid=fill(c.caseid, length(t)),
        time_after_induction_s=t .- BASELINE_SECONDS,
        pam_observed_mmHg=obs_pam, pam_predicted_mmHg=pred_pam,
        pam_error_mmHg=pred_pam .- obs_pam,
        pas_observed_mmHg=obs_pas, pad_observed_mmHg=obs_pad,
        hr_observed_bpm=obs_hr, hr_predicted_bpm=pred_hr,
        hr_error_bpm=pred_hr .- obs_hr,
        co_observed_l_min=obs_co, co_predicted_l_min=pred_co,
        co_error_l_min=pred_co .- obs_co,
    )
end

function effect_from_vector(x::AbstractVector{<:Real})
    # x = [Ce50, gamma, Imax_Rp, Imax_VS, Imax_FC]
    return EfeitoHemodinamico(Ce50=Float64(x[1]), γ=Float64(x[2]),
        Imax_Rp=Float64(x[3]), Imax_VS=Float64(x[4]), Imax_FC=Float64(x[5]))
end

function calibration_objective(x, train_cases)
    # Limites são restrições de otimização, não valores medidos de pacientes.
    lower = [0.1, 0.5, 0.0, 0.0, 0.0]
    upper = [10.0, 5.0, 0.95, 0.95, 0.95]
    if any(x .< lower) || any(x .> upper)
        distance = sum(max.(lower .- x, 0.0).^2) + sum(max.(x .- upper, 0.0).^2)
        return 1e6 + 1e6 * distance
    end
    effect = effect_from_vector(x)
    losses = Float64[]
    for c in train_cases
        try
            prediction = predict_case(c, effect)
            bat, ts = prediction.bat, prediction.t_s
            case_losses = Float64[]
            for (pred, tobs, yobs, scale) in (
                (Float64.(bat.PAM), c.t_pam, c.y_pam, max(abs(c.patient.pam), 1.0)),
                (Float64.(bat.FC), c.t_hr, c.y_hr, max(abs(c.patient.fc), 1.0)),
                (Float64.(bat.DC), c.t_co, c.y_co, max(abs(c.patient.dc::Float64), 0.5)),
            )
                _, p, o = aligned_pair(ts, pred, tobs, yobs)
                length(p) >= 30 || return 1e6
                push!(case_losses, mean(((p .- o) ./ scale).^2))
            end
            push!(losses, mean(case_losses))
        catch
            return 1e6
        end
    end
    return isempty(losses) ? 1e6 : mean(losses)
end

function main()
    cohort_path = joinpath(DATA_DIR, "coorte.csv")
    if !isfile(cohort_path)
        error("Não encontrei $cohort_path. Execute primeiro scripts/vitaldb/01_extrair_dados.py")
    end
    mkpath(OUT_DIR)
    pred_dir = joinpath(OUT_DIR, "predicoes")
    fig_dir = joinpath(OUT_DIR, "figuras")
    mkpath(pred_dir); mkpath(fig_dir)

    cohort = CSV.read(cohort_path, DataFrame; missingstring=["", "NA", "NaN"])
    isempty(cohort) && error("coorte.csv está vazio")
    sort!(cohort, :caseid)
    prepared = PreparedCase[]
    rejected = DataFrame(caseid=Int[], reason=String[])
    for (i, row) in enumerate(eachrow(cohort))
        partition = mod(i, 5) == 0 ? "teste" : "treino"
        try
            c = prepare_case(row, partition)
            push!(prepared, c)
            println("Caso $(c.caseid) carregado; partição=$(c.partition); coadministração=$(c.coadmin_status)")
        catch e
            push!(rejected, (Int(row.caseid), sprint(showerror, e)))
            @warn "Caso excluído da validação" caseid=row.caseid erro=sprint(showerror, e)
        end
    end
    CSV.write(joinpath(OUT_DIR, "casos_excluidos_julia.csv"), rejected)
    length(prepared) >= 3 || error("Menos de 3 casos processáveis. Revise os arquivos e os critérios de seleção.")

    default_effect = EfeitoHemodinamico()
    # Ajuste somente no conjunto de treino e somente quando as trilhas disponíveis
    # não registram infusão concomitante de outro fármaco. 'unknown' não é tratado como zero.
    train_clean = [c for c in prepared if c.partition == "treino" && c.coadmin_status == "none_recorded"]
    test_clean = [c for c in prepared if c.partition == "teste" && c.coadmin_status == "none_recorded"]
    calibrated_effect = default_effect
    calibration_status = "não realizada: poucos casos limpos independentes"
    objective_value = NaN
    if length(train_clean) >= 3 && length(test_clean) >= 2
        println("\nCalibrando Ce50, gamma e Imax no treino: $(length(train_clean)) casos; teste limpo separado: $(length(test_clean)) casos.")
        x0 = [3.0, 2.0, 0.25, 0.15, 0.05]
        objective(x) = calibration_objective(x, train_clean)
        result = Optim.optimize(objective, x0, Optim.NelderMead(); iterations=100, show_trace=false, store_trace=false)
        xfit = Optim.minimizer(result)
        objective_value = Optim.minimum(result)
        calibrated_effect = effect_from_vector(xfit)
        calibration_status = string("realizada; convergiu=", Optim.converged(result), "; iterações=", Optim.iterations(result))
        println("Parâmetros ajustados: Ce50=$(xfit[1]), gamma=$(xfit[2]), Imax_Rp=$(xfit[3]), Imax_VS=$(xfit[4]), Imax_FC=$(xfit[5])")
    else
        @warn "Calibração omitida para evitar ajuste sem treino/teste independentes suficientes" treino_limpo=length(train_clean) teste_limpo=length(test_clean)
    end

    params = DataFrame(
        status=[calibration_status], objective=[objective_value],
        pk_model=[string(PK_MODEL)],
        Ce50_default=[default_effect.Ce50], gamma_default=[default_effect.γ],
        Imax_Rp_default=[default_effect.Imax_Rp], Imax_VS_default=[default_effect.Imax_VS], Imax_FC_default=[default_effect.Imax_FC],
        Ce50_calibrated=[calibrated_effect.Ce50], gamma_calibrated=[calibrated_effect.γ],
        Imax_Rp_calibrated=[calibrated_effect.Imax_Rp], Imax_VS_calibrated=[calibrated_effect.Imax_VS], Imax_FC_calibrated=[calibrated_effect.Imax_FC],
        n_cases=[length(prepared)], n_train_clean=[length(train_clean)], n_test_clean=[length(test_clean)],
    )
    CSV.write(joinpath(OUT_DIR, "parametros_calibrados.csv"), params)

    method_calibrated = startswith(calibration_status, "realizada") ? "calibrado" : "sem_calibracao"
    metric_rows = NamedTuple[]
    test_pam_observed = Float64[]
    test_pam_default = Float64[]
    test_pam_calibrated = Float64[]
    for c in prepared
        for (method, effect) in (("padrao", default_effect), (method_calibrated, calibrated_effect))
            push!(metric_rows, case_metrics(c, effect, method))
            pred_df = prediction_dataframe(c, effect)
            CSV.write(joinpath(pred_dir, @sprintf("case_%05d_%s.csv", c.caseid, method)), pred_df)
            if method == "padrao"
                pred_default = pred_df
            else
                pred_calibrated = pred_df
            end
        end

        # Gráfico de comparação por caso. A linha observada é a pressão real;
        # previsões padrão e calibrada são mostradas separadamente.
        plt = plot(pred_default.time_after_induction_s ./ 60,
            pred_default.pam_observed_mmHg; label="PAM observada", linewidth=2,
            xlabel="Tempo após início da infusão (min)", ylabel="PAM (mmHg)",
            title="VitalDB caso $(c.caseid) — $(c.partition)", legend=:best)
        plot!(plt, pred_default.time_after_induction_s ./ 60, pred_default.pam_predicted_mmHg;
            label="Modelo — parâmetros padrão", linewidth=2, linestyle=:dash)
        plot!(plt, pred_calibrated.time_after_induction_s ./ 60, pred_calibrated.pam_predicted_mmHg;
            label=(method_calibrated == "calibrado" ? "Modelo — calibrado no treino" : "Modelo — sem calibração"), linewidth=2)
        savefig(plt, joinpath(fig_dir, @sprintf("case_%05d_pam.png", c.caseid)))

        if c.partition == "teste"
            for i in eachindex(pred_calibrated.pam_observed_mmHg)
                obs = pred_calibrated.pam_observed_mmHg[i]
                dflt = pred_default.pam_predicted_mmHg[i]
                cal = pred_calibrated.pam_predicted_mmHg[i]
                if isfinite(obs) && isfinite(dflt) && isfinite(cal)
                    push!(test_pam_observed, obs)
                    push!(test_pam_default, dflt)
                    push!(test_pam_calibrated, cal)
                end
            end
        end
    end
    metrics_df = DataFrame(metric_rows)
    CSV.write(joinpath(OUT_DIR, "metricas_por_caso.csv"), metrics_df)

    summary_rows = NamedTuple[]
    for method in ("padrao", method_calibrated)
        for stratum in ("all_test", "none_recorded", "active", "unknown")
            sub = filter(r -> r.method == method && r.partition == "teste" &&
                (stratum == "all_test" || r.coadmin_status == stratum), metrics_df)
            for (label, col) in (("PAM_MAE_mmHg", :pam_mae_mmHg), ("PAM_RMSE_mmHg", :pam_rmse_mmHg),
                                 ("PAM_bias_mmHg", :pam_bias_mmHg), ("HR_MAE_bpm", :hr_mae_bpm),
                                 ("CO_MAE_L_min", :co_mae_l_min))
                vals = [Float64(v) for v in sub[!, col] if !ismissing(v) && isfinite(Float64(v))]
                push!(summary_rows, (method=method, coadmin_stratum=stratum, metric=label,
                    mean_across_test_cases=isempty(vals) ? NaN : mean(vals),
                    median_across_test_cases=isempty(vals) ? NaN : median(vals), n_test_cases=length(vals)))
            end
        end
    end
    CSV.write(joinpath(OUT_DIR, "resumo_teste.csv"), DataFrame(summary_rows))

    if !isempty(test_pam_observed)
        lims = (min(minimum(test_pam_observed), minimum(test_pam_default), minimum(test_pam_calibrated)),
                max(maximum(test_pam_observed), maximum(test_pam_default), maximum(test_pam_calibrated)))
        scatter(test_pam_observed, test_pam_default; label="Parâmetros padrão", alpha=0.45,
            xlabel="PAM observada (mmHg)", ylabel="PAM prevista (mmHg)",
            title="Teste independente — modelo padrão vs observado")
        plot!(lims, lims; label="Igualdade", linestyle=:dash, linewidth=2)
        savefig(joinpath(fig_dir, "dispersao_pam_teste_padrao.png"))

        scatter(test_pam_observed, test_pam_calibrated; label="Parâmetros calibrados", alpha=0.45,
            xlabel="PAM observada (mmHg)", ylabel="PAM prevista (mmHg)",
            title=(method_calibrated == "calibrado" ? "Teste independente — modelo calibrado vs observado" : "Teste independente — modelo sem calibração vs observado"))
        plot!(lims, lims; label="Igualdade", linestyle=:dash, linewidth=2)
        savefig(joinpath(fig_dir, method_calibrated == "calibrado" ? "dispersao_pam_teste_calibrado.png" : "dispersao_pam_teste_sem_calibracao.png"))
    end

    println("\nValidação finalizada.")
    println("Casos processados: $(length(prepared)); treino limpo: $(length(train_clean)); teste limpo: $(length(test_clean))")
    println("Saídas: $OUT_DIR")
    println("Arquivos principais: metricas_por_caso.csv, resumo_teste.csv, parametros_calibrados.csv e figuras/*.png")
    println("IMPORTANTE: use as métricas do conjunto de teste para avaliar generalização; não avalie precisão só pelos casos de treino.")
end

main()
