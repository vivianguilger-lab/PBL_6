# =============================================================================
# populacao.jl — coorte virtual (virtual clinical trial) e variabilidade
# =============================================================================
#
# EXPERIMENTAL — NÃO VALIDADO. As variabilidades abaixo e as distribuições da
# coorte são suposições; os resultados não são usados na plataforma nem nos
# resultados do README (pendências em docs/validacao.md). Leitura das variâncias
# de Schnider (ω² ou coeficiente de variação) ainda a confirmar na tabela do artigo.
#
# Variabilidade interindividual (IIV) log-normal: θᵢ = θ·exp(ηᵢ), ηᵢ ~ N(0, ω²)
#   PK  — variâncias de Schnider et al. 1998 (V1, V2, V3, CL, Q2, Q3); ke0: ω = 0,3 (assumido)
#   Hemodinâmica — Ce50 e Imax do EfeitoHemodinamico do grupo: ω = 0,3 (assumido,
#   os próprios valores centrais são ilustrativos no modelo do grupo)
# =============================================================================

const OMEGA_PK = (V1 = sqrt(0.0404), V2 = sqrt(0.01), V3 = sqrt(0.1435),
                  CL = sqrt(0.1005), Q2 = sqrt(0.01), Q3 = sqrt(0.1179), ke0 = 0.3)
const OMEGA_HEMO = (Ce50 = 0.3, Imax = 0.3)

"""Fatores multiplicativos exp(η) da PK, η ~ N(0, ω²)."""
amostrar_eta_pk(rng) = NamedTuple{keys(OMEGA_PK)}(Tuple(exp(w * randn(rng)) for w in values(OMEGA_PK)))

"""Efeito hemodinâmico de um indivíduo sorteado em torno de `ef`."""
function amostrar_efeito(rng, ef::EfeitoHemodinamico = EfeitoHemodinamico())
    fI() = exp(OMEGA_HEMO.Imax * randn(rng))
    return EfeitoHemodinamico(Ce50 = ef.Ce50 * exp(OMEGA_HEMO.Ce50 * randn(rng)), γ = ef.γ,
                              Imax_Rp = min(ef.Imax_Rp * fI(), 0.9), Imax_VS = min(ef.Imax_VS * fI(), 0.9),
                              Imax_FC = min(ef.Imax_FC * fI(), 0.9))
end

_trunc(rng, μ, σ, a, b) = (x = μ + σ * randn(rng); while !(a <= x <= b); x = μ + σ * randn(rng); end; x)

"Um indivíduo da coorte: dados pré-indução + variabilidade + dose do protocolo."
struct Individuo
    paciente::Paciente
    eta_pk::NamedTuple
    efeito::EfeitoHemodinamico
    dose_mgkg::Float64
end

"""
    dose_padrao(idade) -> mg/kg

"Faixa por peso, ajustada por idade": 2,0 mg/kg (< 55 anos), 1,5 (55–69), 1,2 (≥ 70).
"""
dose_padrao(idade) = idade < 55 ? 2.0 : (idade < 70 ? 1.5 : 1.2)

"""
    gerar_coorte(n; rng=Xoshiro(2026))

`n` adultos virtuais com demografia, PA/FC pré-indução, ASA, PVI, IVCI e IIV de
PK e da resposta hemodinâmica. Experimental: as distribuições são suposições, sem
fonte (docs/validacao.md).
"""
function gerar_coorte(n; rng = Xoshiro(2026), variabilidade_dose = 0.1)
    coorte = Individuo[]
    for i in 1:n
        idade = _trunc(rng, 55.0, 16.0, 18.0, 90.0)
        sexo = rand(rng) < 0.5 ? :M : :F
        altura = sexo == :M ? _trunc(rng, 174.0, 7.0, 150.0, 200.0) : _trunc(rng, 162.0, 6.5, 140.0, 190.0)
        peso = _trunc(rng, 27.0, 4.5, 17.0, 40.0) * (altura / 100)^2
        pas = _trunc(rng, 115.0 + 0.45 * (idade - 18), 16.0, 90.0, 210.0)
        pad = clamp(_trunc(rng, 76.0 + 0.15 * (idade - 18) - 0.004 * (idade - 18)^2, 9.0, 45.0, 115.0), pas - 100, pas - 30)
        fc = _trunc(rng, 76.0, 12.0, 48.0, 115.0)
        u = rand(rng)
        asa = idade < 40 ? (u < 0.6 ? 1 : 2) : idade < 65 ? (u < 0.2 ? 1 : u < 0.8 ? 2 : 3) : (u < 0.55 ? 2 : u < 0.95 ? 3 : 4)
        pvi = _trunc(rng, 12.0, 5.0, 3.0, 35.0)
        ivci = _trunc(rng, 25.0, 13.0, 0.0, 80.0)
        # dados sintéticos: a PAM do indivíduo virtual é gerada a partir de PAS e PAD
        p = Paciente(nome = "virtual_$(i)", idade = idade, peso = peso, altura = altura, sexo = sexo, asa = asa,
                     pas = pas, pad = pad, pam = pad + (pas - pad) / 3, fc = fc, pvi = pvi, ivci = ivci)
        dose = dose_padrao(idade) * (1 + variabilidade_dose * (2 * rand(rng) - 1))
        push!(coorte, Individuo(p, amostrar_eta_pk(rng), amostrar_efeito(rng), dose))
    end
    return coorte
end

"Calibra cada paciente virtual com os SEUS dados pré-indução (em paralelo com `julia -t auto`)."
function calibrar_coorte(coorte)
    pvs = Vector{PacienteVirtual}(undef, length(coorte))
    Threads.@threads for i in eachindex(coorte)
        pvs[i] = calibrar(coorte[i].paciente; eta_pk = coorte[i].eta_pk, efeito = coorte[i].efeito)
    end
    return pvs
end

"""
    simular_coorte(coorte, pvs=calibrar_coorte(coorte); manutencao_ugkgmin=100.0, duracao_bolus=30.0,
                   janela=20.0, dose_fixa=nothing)
"""
function simular_coorte(coorte, pvs = calibrar_coorte(coorte); manutencao_ugkgmin = 100.0, duracao_bolus = 30.0,
                        janela = 20.0, dose_fixa = nothing)
    saida = Vector{Any}(undef, length(coorte))
    Threads.@threads for i in eachindex(coorte)
        ind, pv = coorte[i], pvs[i]
        dose = dose_fixa === nothing ? ind.dose_mgkg : dose_fixa
        mt = metricas(simular_inducao(pv; dose_mgkg = dose, duracao_bolus, manutencao_ugkgmin, duracao = janela); janela)
        p = ind.paciente
        saida[i] = (idade = p.idade, sexo = p.sexo, peso = p.peso, asa = p.asa, pvi = p.pvi, ivci = p.ivci,
                    s = pv.s, pam_basal = mt.pam_basal, dose_mgkg = dose, pam_nadir = mt.pam_nadir,
                    queda_pct = mt.queda_pct, t_nadir = mt.t_nadir,
                    min_abaixo_65 = mt.tempo_abaixo[65.0], min_abaixo_55 = mt.tempo_abaixo[55.0],
                    bis_min = mt.bis_min, risco = mt.risco, calibrou = pv.convergiu)
    end
    return [s for s in saida]
end

"""Resumo da coorte: incidências e estatísticas."""
function resumo_coorte(res)
    n = length(res)
    f(c) = 100 * count(c, res) / n
    q(x, p) = quantile([getfield(r, x) for r in res], p)
    return (n = n, inc_pam65 = f(r -> r.min_abaixo_65 > 0), inc_pam55 = f(r -> r.min_abaixo_55 > 0),
            inc_queda30 = f(r -> r.queda_pct > 30),
            queda_mediana = q(:queda_pct, 0.5), queda_p25 = q(:queda_pct, 0.25), queda_p75 = q(:queda_pct, 0.75),
            t_nadir_mediano = q(:t_nadir, 0.5), bis_ok = f(r -> r.bis_min < 60), calibrados = f(r -> r.calibrou))
end
