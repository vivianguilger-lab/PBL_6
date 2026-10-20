# =============================================================================
# plataforma.jl — interface web local para uso clínico
# =============================================================================
#
#   using PropofolWindkessel
#   plataforma()          # abre o navegador em http://localhost:8080
#
# O navegador mostra um formulário (dados pré-indução + dose) e os gráficos;
# todo o cálculo continua sendo feito por este pacote em Julia. O servidor só
# escuta em 127.0.0.1: os dados do paciente não saem do computador.
# =============================================================================

const _HTML = joinpath(@__DIR__, "..", "app", "index.html")
const _CACHE_PV = Dict{Any,PacienteVirtual}()
const _TRAVA = ReentrantLock()

_num(x) = x isa Number ? Float64(x) : parse(Float64, replace(strip(string(x)), "," => "."))
_opcional(d, k) = (v = get(d, k, nothing); (v === nothing || v == "") ? nothing : _num(v))
_r(x, n = 2) = isfinite(x) ? round(x; digits = n) : nothing
_rv(v, n = 2) = [_r(x, n) for x in v]

function _ler_paciente(d)
    pac = d["paciente"]
    sexo = uppercase(string(get(pac, "sexo", "M")))
    sexo in ("M", "F") || throw(ArgumentError("Sexo deve ser M ou F."))
    idade, peso, altura = _num(pac["idade"]), _num(pac["peso"]), _num(pac["altura"])
    pas, pad, fc = _num(pac["pas"]), _num(pac["pad"]), _num(pac["fc"])
    18 <= idade <= 95 || throw(ArgumentError("Idade fora da faixa do modelo (18–95 anos)."))
    35 <= peso <= 200 || throw(ArgumentError("Peso fora da faixa do modelo (35–200 kg)."))
    130 <= altura <= 215 || throw(ArgumentError("Altura fora da faixa do modelo (130–215 cm)."))
    60 <= pas <= 250 && 30 <= pad <= 150 || throw(ArgumentError("Pressão arterial fora da faixa plausível."))
    pas - pad >= 15 || throw(ArgumentError("A PAS deve ser pelo menos 15 mmHg maior que a PAD."))
    35 <= fc <= 160 || throw(ArgumentError("Frequência cardíaca fora da faixa plausível (35–160 bpm)."))
    return Paciente(nome = string(get(pac, "nome", "Paciente")), idade = idade, peso = peso, altura = altura,
                    sexo = Symbol(sexo), pas = pas, pad = pad, fc = fc,
                    pam = _opcional(pac, "pam"), dc = _opcional(pac, "dc"))
end

function _ler_protocolo(d)
    pr = get(d, "protocolo", Dict())
    dose = _num(get(pr, "dose_mgkg", 2.0))
    0 <= dose <= 4 || throw(ArgumentError("Dose fora da faixa do simulador (0–4 mg/kg)."))
    return (dose_mgkg = dose, duracao_bolus = max(_num(get(pr, "duracao_bolus_s", 20.0)), 1.0),
            manutencao_ugkgmin = max(_num(get(pr, "manutencao_ugkgmin", 0.0)), 0.0))
end

"""Paciente virtual calibrado (guardado em cache: recalibrar só se os dados mudarem)."""
function _virtual(p::Paciente)
    chave = (p.idade, p.peso, p.altura, p.sexo, p.pas, p.pad, p.fc, p.pam, p.dc)
    lock(_TRAVA) do
        length(_CACHE_PV) > 200 && empty!(_CACHE_PV)
        get!(() -> calibrar(p), _CACHE_PV, chave)
    end
end

function _resumo_virtual(pv::PacienteVirtual)
    cv, o = pv.cv, pv.obtido
    return Dict("convergiu" => pv.convergiu,
                "pam" => _r(o.PAM, 1), "pas" => _r(o.PAS, 1), "pad" => _r(o.PAD, 1),
                "pas_aorta" => _r(o.PASao, 1), "pad_aorta" => _r(o.PADao, 1),
                "dc" => _r(o.DC), "dc_estimado" => pv.paciente.dc === nothing, "vs" => _r(o.VS, 0),
                "tpr" => _r(cv.Rs, 3), "ac" => _r(cv.Ca, 2), "tau" => _r(cv.Rs * cv.Ca, 2), "zc" => _r(cv.Zc, 3),
                "cs" => _r(1 / cv.Emax, 2), "cd" => _r(1 / cv.Emin, 1),
                "volume" => _r(cv.Vtotal, 0), "v0v" => _r(cv.Vv0, 0), "pms" => _r(pressao_media_enchimento(cv), 1),
                "lbm" => _r(massa_magra(pv.paciente), 1), "bsa" => _r(superficie_corporal(pv.paciente), 2),
                "v1" => _r(pv.pk.V1, 2), "v2" => _r(pv.pk.V2, 1), "cl1" => _r(pv.pk.Cl1, 2),
                "t12_vascular" => _r(log(2) / pv.pk.ke0_hemo / 60, 1),
                "amp" => _r(pv.obs.amp, 2))
end

function _api_simular(d)
    p = _ler_paciente(d)
    pr = _ler_protocolo(d)
    pv = _virtual(p)
    r = simular_inducao(pv; pr..., saveat = 0.01)
    mt = metricas(r)
    b, s, cv = r.bat, r.sinais, pv.cv
    t0 = r.t_inducao
    # sinais contínuos a cada 1 s
    idx = 1:100:length(r.t)
    Vest(v) = cv.Vtotal - cv.Va0 - v - cv.Vd
    i_nadir = argmin(abs.(r.t .- (t0 + 60 * mt.t_nadir)))
    j_nadir = argmin(abs.(b.t .- (t0 + 60 * mt.t_nadir)))
    return Dict(
        "virtual" => _resumo_virtual(pv),
        "dose_mg" => _r(pr.dose_mgkg * p.peso, 0),
        "metricas" => Dict(
            "pam_basal" => _r(mt.pam_basal, 1), "pam_nadir" => _r(mt.pam_nadir, 1),
            "queda_pct" => _r(mt.queda_pct, 1), "t_nadir" => _r(mt.t_nadir, 1),
            "pas_nadir" => _r(mt.pas_nadir, 0), "pad_nadir" => _r(b.PAD[j_nadir], 0),
            "fc_basal" => _r(mt.fc_basal, 0), "fc_nadir" => _r(mt.fc_no_nadir, 0),
            "dc_basal" => _r(mt.dc_basal), "dc_nadir" => _r(mt.dc_no_nadir),
            "min_abaixo_65" => _r(mt.tempo_abaixo[65.0], 1), "min_abaixo_55" => _r(mt.tempo_abaixo[55.0], 1),
            "aut_65" => _r(mt.aut[65.0], 0), "aut_55" => _r(mt.aut[55.0], 0),
            "bis_min" => _r(mt.bis_min, 0), "t_hipnose" => _r(mt.t_hipnose, 1),
            "tempo_hipnose" => _r(mt.tempo_hipnose, 1),
            "ceb_max" => _r(mt.ceb_max), "ceh_max" => _r(mt.ceh_max), "risco" => mt.risco,
            "tpr_nadir_pct" => _r(100 * (s.Rs[i_nadir] / cv.Rs - 1), 0),
            "contr_nadir_pct" => _r(100 * (s.Emax[i_nadir] / cv.Emax - 1), 0),
            "vest_nadir_pct" => _r(100 * (Vest(s.Vv0[i_nadir]) / Vest(cv.Vv0) - 1), 0)),
        "batimentos" => Dict(
            "t" => _rv((b.t .- t0) ./ 60, 3), "pam" => _rv(b.PAM, 1), "pas" => _rv(b.PAS, 1),
            "pad" => _rv(b.PAD, 1), "fc" => _rv(b.FC, 1), "dc" => _rv(b.DC, 2)),
        "sinais" => Dict(
            "t" => _rv((r.t[idx] .- t0) ./ 60, 3), "bis" => _rv(s.BIS[idx], 1),
            "cp" => _rv(s.Cp[idx], 2), "ceb" => _rv(s.Ceb[idx], 2), "ceh" => _rv(s.Ceh[idx], 2),
            "tpr" => _rv(100 .* (s.Rs[idx] ./ cv.Rs .- 1), 1),
            "contr" => _rv(100 .* (s.Emax[idx] ./ cv.Emax .- 1), 1),
            "vest" => _rv(100 .* (Vest.(s.Vv0[idx]) ./ Vest(cv.Vv0) .- 1), 1)))
end

function _api_incerteza(d)
    p = _ler_paciente(d)
    pr = _ler_protocolo(d)
    n = clamp(round(Int, _num(get(d, "n", 50))), 10, 300)
    inc = previsao_com_incerteza(_virtual(p); n, pr...)
    return Dict("t" => _rv(inc.t, 3), "p05" => _rv(inc.p05, 1), "mediana" => _rv(inc.mediana, 1),
                "p95" => _rv(inc.p95, 1), "prob65" => _r(inc.prob_pam65, 3), "prob55" => _r(inc.prob_pam55, 3),
                "prob_hipnose" => _r(inc.prob_hipnose, 3), "nadir_mediano" => _r(inc.nadir_mediano, 0),
                "nadir_p05" => _r(inc.nadir_p05, 0), "nadir_p95" => _r(inc.nadir_p95, 0), "n" => n)
end

"""Interpola o ponto em que `y` cruza `alvo` numa varredura monótona."""
function _cruzamento(x, y, alvo)
    for i in 1:length(x)-1
        if (y[i] - alvo) * (y[i+1] - alvo) <= 0 && y[i] != y[i+1]
            return x[i] + (alvo - y[i]) * (x[i+1] - x[i]) / (y[i+1] - y[i])
        end
    end
    return nothing
end

function _api_doses(d)
    p = _ler_paciente(d)
    pr = _ler_protocolo(d)
    pv = _virtual(p)
    doses = collect(0.5:0.25:3.0)
    tab = varrer_doses(pv, doses; duracao_bolus = pr.duracao_bolus, manutencao_ugkgmin = pr.manutencao_ugkgmin)
    nadir = [t.pam_nadir for t in tab]
    bis = [t.bis_min for t in tab]
    dmin = bis[1] < 60 ? doses[1] : _cruzamento(doses, bis, 60.0)
    dmax = nadir[end] >= 65 ? doses[end] : (nadir[1] < 65 ? nothing : _cruzamento(doses, nadir, 65.0))
    return Dict("doses" => [Dict("dose" => t.dose, "mg" => _r(t.dose * p.peso, 0), "pam_nadir" => _r(t.pam_nadir, 1),
                                 "queda_pct" => _r(t.queda_pct, 0), "min_abaixo_65" => _r(t.tempo_abaixo[65.0], 1),
                                 "min_abaixo_55" => _r(t.tempo_abaixo[55.0], 1), "bis_min" => _r(t.bis_min, 0),
                                 "risco" => t.risco) for t in tab],
                "janela" => Dict("dose_min" => dmin === nothing ? nothing : _r(dmin, 2),
                                 "dose_max" => dmax === nothing ? nothing : _r(dmax, 2),
                                 "existe" => dmin !== nothing && dmax !== nothing && dmin <= dmax))
end

_json(status, x) = HTTP.Response(status, ["Content-Type" => "application/json; charset=utf-8"], JSON.json(x))

function _rotear(req::HTTP.Request)
    try
        alvo = first(split(req.target, '?'))
        if req.method == "GET" && alvo in ("/", "/index.html")
            return HTTP.Response(200, ["Content-Type" => "text/html; charset=utf-8"], read(_HTML, String))
        elseif req.method == "GET" && alvo == "/api/saude"
            return _json(200, Dict("ok" => true, "versao" => string(pkgversion(@__MODULE__))))
        elseif req.method == "POST" && alvo in ("/api/simular", "/api/incerteza", "/api/doses")
            d = JSON.parse(String(req.body))
            f = alvo == "/api/simular" ? _api_simular : alvo == "/api/incerteza" ? _api_incerteza : _api_doses
            return _json(200, f(d))
        end
        return HTTP.Response(404, "não encontrado")
    catch e
        msg = e isa ArgumentError ? e.msg : e isa KeyError ? "Campo obrigatório ausente: $(e.key)" :
              "Erro no cálculo: " * sprint(showerror, e)
        e isa ArgumentError || e isa KeyError || @error "plataforma" exception = (e, catch_backtrace())
        return _json(400, Dict("erro" => msg))
    end
end

function _abrir_navegador(url)
    try
        if Sys.iswindows()
            run(`cmd /c start "" $url`; wait = false)
        elseif Sys.isapple()
            run(`open $url`; wait = false)
        else
            run(`xdg-open $url`; wait = false)
        end
    catch
    end
end

"""
    plataforma(; porta=8080, abrir=true) -> servidor

Abre a plataforma web do simulador no navegador (http://localhost:porta).
O cálculo é feito por este pacote; o servidor só aceita conexões do próprio
computador. Para encerrar: `close(servidor)` ou feche o Julia.
"""
function plataforma(; porta = 8080, abrir = true, aquecer = true)
    if aquecer
        println("Preparando o simulador (na primeira vez leva cerca de 1 minuto)...")
        ex = Paciente(idade = 50, peso = 70, altura = 170, sexo = :M, pas = 125, pad = 75, fc = 70)
        corpo = Dict("paciente" => Dict("idade" => 50, "peso" => 70, "altura" => 170, "sexo" => "M",
                                        "pas" => 125, "pad" => 75, "fc" => 70),
                     "protocolo" => Dict("dose_mgkg" => 1.0), "n" => 10)
        JSON.json(_api_simular(corpo))
        previsao_com_incerteza(_virtual(ex); n = 2, janela = 60.0)
        metricas(simular_inducao(_virtual(ex); duracao = 60.0, saveat = 0.01); janela = 60.0)
    end
    srv = nothing
    for tentativa in 0:9                      # se a porta estiver ocupada, tenta a seguinte
        try
            srv = HTTP.serve!(_rotear, "127.0.0.1", porta + tentativa; verbose = false)
            porta += tentativa
            break
        catch e
            tentativa == 9 && rethrow()
        end
    end
    url = "http://localhost:$porta"
    println("\nPlataforma pronta: $url")
    println("Se o navegador não abrir sozinho, copie o endereço acima no navegador.")
    println("Para encerrar: feche esta janela do Julia.")
    abrir && _abrir_navegador(url)
    return srv
end
