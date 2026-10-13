# =============================================================================
# metricas.jl — métricas de hipotensão pós-indução e relatório
# =============================================================================
#
# Hipotensão pós-indução: queda de PA nos primeiros 20 min após a indução.
# 17–24% dos adultos têm PAM < 55 mmHg; PAM ≤ 55 mmHg sustentada associa-se
# a lesão miocárdica e lesão renal aguda. Por isso o relatório usa dois
# limiares (65 e 55 mmHg) e a área abaixo do limiar (AUT, mmHg·min).
# =============================================================================

"""
    metricas(r::Resultado; janela=20*60, limiares=(65.0, 55.0), alvo_bis=60.0)

Resumo clínico da curva prevista (apenas batimentos dentro da janela de
`janela` segundos após o início do bolus).
"""
function metricas(r::Resultado; janela = 20 * 60.0, limiares = (65.0, 55.0), alvo_bis = 60.0)
    b = r.bat
    t0 = r.t_inducao
    base = findall(t -> t < t0, b.t)
    pam0 = isempty(base) ? b.PAM[1] : mean(b.PAM[base])
    fc0 = isempty(base) ? b.FC[1] : mean(b.FC[base])
    dc0 = isempty(base) ? b.DC[1] : mean(b.DC[base])
    jan = findall(t -> t0 <= t <= t0 + janela, b.t)
    pam = b.PAM[jan]
    tb = b.t[jan]
    dtb = 60 ./ b.FC[jan] ./ 60                        # duração de cada batimento [min]
    inad = argmin(pam)
    nadir = pam[inad]
    tempo_abaixo = Dict(l => sum(dtb[pam .< l]; init = 0.0) for l in limiares)
    aut = Dict(l => sum(max.(l .- pam, 0.0) .* dtb; init = 0.0) for l in limiares)
    # hipnose
    js = findall(t -> t0 <= t <= t0 + janela, r.t)
    bis = r.sinais.BIS[js]
    ib = findfirst(<(alvo_bis), bis)
    t_hipnose = ib === nothing ? NaN : (r.t[js][ib] - t0) / 60
    dt_s = length(r.t) > 1 ? r.t[2] - r.t[1] : 0.0
    tempo_hipnose = count(<(alvo_bis), bis) * dt_s / 60
    lmin = minimum(limiares)
    lmax = maximum(limiares)
    risco = tempo_abaixo[lmin] > 0 ? "ALTO" :
            (tempo_abaixo[lmax] >= 1.0 || nadir < 0.7 * pam0) ? "MODERADO" : "BAIXO"
    return (pam_basal = pam0, fc_basal = fc0, dc_basal = dc0,
            pam_nadir = nadir, t_nadir = (tb[inad] - t0) / 60, queda_pct = 100 * (1 - nadir / pam0),
            pas_nadir = minimum(b.PAS[jan]), fc_no_nadir = b.FC[jan][inad], dc_no_nadir = b.DC[jan][inad],
            tempo_abaixo = tempo_abaixo, aut = aut,
            bis_min = minimum(bis), t_hipnose = t_hipnose, tempo_hipnose = tempo_hipnose,
            ceb_max = maximum(r.sinais.Ceb), ceh_max = maximum(r.sinais.Ceh),
            risco = risco, limiares = limiares, alvo_bis = alvo_bis)
end

"""Imprime um relatório legível das métricas."""
function relatorio(io::IO, mt; titulo = "Previsão hemodinâmica")
    println(io, "─"^72)
    println(io, titulo)
    println(io, "─"^72)
    @printf(io, "PAM basal            : %6.1f mmHg   (FC %.0f bpm, DC %.2f L/min)\n", mt.pam_basal, mt.fc_basal, mt.dc_basal)
    @printf(io, "PAM nadir            : %6.1f mmHg   em %.1f min (queda de %.0f%%)\n", mt.pam_nadir, mt.t_nadir, mt.queda_pct)
    @printf(io, "  no nadir           : FC %.0f bpm, DC %.2f L/min, PAS %.0f mmHg\n", mt.fc_no_nadir, mt.dc_no_nadir, mt.pas_nadir)
    for l in sort(collect(mt.limiares); rev = true)
        @printf(io, "Tempo com PAM < %2.0f   : %6.1f min    AUT = %.0f mmHg·min\n", l, mt.tempo_abaixo[l], mt.aut[l])
    end
    @printf(io, "BIS mínimo           : %6.0f        (BIS < %.0f em %.1f min; por %.1f min)\n",
            mt.bis_min, mt.alvo_bis, mt.t_hipnose, mt.tempo_hipnose)
    @printf(io, "Ce máx. (SNC/vasos)  : %.2f / %.2f µg/ml\n", mt.ceb_max, mt.ceh_max)
    println(io, "Risco previsto       : ", mt.risco)
    println(io, "─"^72)
end
relatorio(mt; kw...) = relatorio(stdout, mt; kw...)
