# =============================================================================
# populacao.jl — coorte virtual (virtual clinical trial)
# =============================================================================
#
# Resumo, Seção 8.7: um ensaio clínico virtual tem três blocos —
#   I.  características dos pacientes (demografia + sinais vitais);
#   II. características da coorte e protocolo de tratamento (dose por peso,
#       ajustada por idade, como na prática clínica);
#   III. o modelo (escolha, identificação e validação).
#
# Variabilidade interindividual (IIV) log-normal: θᵢ = θ·exp(ηᵢ), ηᵢ ~ N(0, ω²)
#   PK  — variâncias de Schnider et al. 1998 (V1, V2, V3, Cl1, Cl2, Cl3)
#   PD  — C50 da TPR: ω² = 0,44 (Su et al. 2023); C50 do BIS: CV 18%
#         (Bouillon 2004); demais ω assumidos (0,3).
# =============================================================================

"""Desvios-padrão (ω) da variabilidade interindividual log-normal."""
const OMEGA_PK = (V1 = sqrt(0.0404), V2 = sqrt(0.01), V3 = sqrt(0.1435),
                  Cl1 = sqrt(0.1005), Cl2 = sqrt(0.01), Cl3 = sqrt(0.1179),
                  ke0_bis = 0.3, ke0_hemo = 0.3)
const OMEGA_PD = (tpr_C50 = sqrt(0.44), inot = 0.3, veno = 0.3, baro = 0.2, bis_C50 = 0.18)

"""Fatores multiplicativos exp(η) da PK, η ~ N(0, ω²)."""
amostrar_eta_pk(rng) = NamedTuple{keys(OMEGA_PK)}(Tuple(exp(w * randn(rng)) for w in values(OMEGA_PK)))

"""Parâmetros PD de um indivíduo sorteado em torno de `pd_ref`."""
amostrar_pd(rng, pd_ref = ParametrosPD()) =
    com(pd_ref;
        tpr_C50 = pd_ref.tpr_C50 * exp(OMEGA_PD.tpr_C50 * randn(rng)),
        inot_Emax35 = pd_ref.inot_Emax35 * exp(OMEGA_PD.inot * randn(rng)),
        veno_Emax = pd_ref.veno_Emax * exp(OMEGA_PD.veno * randn(rng)),
        baro_Emax = min(pd_ref.baro_Emax * exp(OMEGA_PD.baro * randn(rng)), 0.95),
        bis_C50 = pd_ref.bis_C50 * exp(OMEGA_PD.bis_C50 * randn(rng)))

_trunc(rng, μ, σ, a, b) = (x = μ + σ * randn(rng); while !(a <= x <= b); x = μ + σ * randn(rng); end; x)

"""
Um indivíduo da coorte: dados pré-indução + fatores de variabilidade.
"""
struct Individuo
    paciente::Paciente
    eta_pk::NamedTuple
    pd::ParametrosPD
    dose_mgkg::Float64
end

"""
    dose_padrao(idade) -> mg/kg

Regra clínica simplificada usada na coorte — "faixa por peso, ajustada
por idade": 2,0 mg/kg (< 55 anos), 1,5 mg/kg (55–69) e 1,2 mg/kg (≥ 70).
"""
dose_padrao(idade) = idade < 55 ? 2.0 : (idade < 70 ? 1.5 : 1.2)

"""
    gerar_coorte(n; rng=Xoshiro(2026), pd_ref=ParametrosPD())

Gera `n` adultos virtuais com demografia, PA/FC pré-indução, DC e IIV de
PK e PD. As distribuições são suposições plausíveis para adultos
cirúrgicos (documentadas em `docs/modelo_matematico.md`).
"""
function gerar_coorte(n; rng = Xoshiro(2026), pd_ref = ParametrosPD(), variabilidade_dose = 0.1)
    coorte = Individuo[]
    for i in 1:n
        idade = _trunc(rng, 55.0, 16.0, 18.0, 90.0)
        sexo = rand(rng) < 0.5 ? :M : :F
        altura = sexo == :M ? _trunc(rng, 174.0, 7.0, 150.0, 200.0) : _trunc(rng, 162.0, 6.5, 140.0, 190.0)
        imc_i = _trunc(rng, 27.0, 4.5, 17.0, 40.0)
        peso = imc_i * (altura / 100)^2
        pas = _trunc(rng, 115.0 + 0.45 * (idade - 18), 16.0, 90.0, 210.0)
        pad = _trunc(rng, 76.0 + 0.15 * (idade - 18) - 0.004 * (idade - 18)^2, 9.0, 45.0, 115.0)
        pad = clamp(pad, pas - 100, pas - 30)
        fc = _trunc(rng, 76.0, 12.0, 48.0, 115.0)
        bsa = sqrt(altura * peso / 3600)
        dc = indice_cardiaco_estimado(idade) * bsa * exp(0.12 * randn(rng))
        p = Paciente(nome = "virtual_$(i)", idade = idade, peso = peso, altura = altura, sexo = sexo,
                     pas = pas, pad = pad, fc = fc, dc = dc)
        eta = amostrar_eta_pk(rng)
        pd = amostrar_pd(rng, pd_ref)
        dose = dose_padrao(idade) * (1 + variabilidade_dose * (2 * rand(rng) - 1))
        push!(coorte, Individuo(p, eta, pd, dose))
    end
    return coorte
end

"""
    calibrar_coorte(coorte) -> Vector{PacienteVirtual}

Calibra cada paciente virtual com os SEUS dados pré-indução (em paralelo:
rode o Julia com `-t auto`).
"""
function calibrar_coorte(coorte)
    pvs = Vector{PacienteVirtual}(undef, length(coorte))
    Threads.@threads for i in eachindex(coorte)
        pvs[i] = calibrar(coorte[i].paciente; pd = coorte[i].pd, eta_pk = coorte[i].eta_pk)
    end
    return pvs
end

"""
    simular_coorte(coorte, pvs=calibrar_coorte(coorte); manutencao_ugkgmin=0.0,
                   duracao_bolus=20.0, janela=20*60, dose_fixa=nothing)

Simula a indução de cada indivíduo com a dose do protocolo (ou `dose_fixa`
mg/kg para todos). Devolve um vetor de NamedTuple (um por indivíduo).
"""
function simular_coorte(coorte, pvs = calibrar_coorte(coorte); manutencao_ugkgmin = 0.0,
                        duracao_bolus = 20.0, janela = 20 * 60.0, dose_fixa = nothing)
    saida = Vector{Any}(undef, length(coorte))
    Threads.@threads for i in eachindex(coorte)
        ind, pv = coorte[i], pvs[i]
        dose = dose_fixa === nothing ? ind.dose_mgkg : dose_fixa
        r = simular_inducao(pv; dose_mgkg = dose, duracao_bolus, manutencao_ugkgmin,
                            duracao = janela, saveat = 0.01)
        mt = metricas(r; janela)
        saida[i] = (idade = ind.paciente.idade, sexo = ind.paciente.sexo, peso = ind.paciente.peso,
                    pam_basal = mt.pam_basal, dose_mgkg = dose, pam_nadir = mt.pam_nadir,
                    queda_pct = mt.queda_pct, t_nadir = mt.t_nadir,
                    min_abaixo_65 = mt.tempo_abaixo[65.0], min_abaixo_55 = mt.tempo_abaixo[55.0],
                    bis_min = mt.bis_min, risco = mt.risco, calibrou = pv.convergiu,
                    tpr_C50 = ind.pd.tpr_C50)
    end
    return [s for s in saida]
end

"""Resumo da coorte: incidências e estatísticas."""
function resumo_coorte(res)
    n = length(res)
    f(c) = 100 * count(c, res) / n
    q(x, p) = quantile([getfield(r, x) for r in res], p)
    return (n = n,
            inc_pam65 = f(r -> r.min_abaixo_65 > 0),
            inc_pam55 = f(r -> r.min_abaixo_55 > 0),
            inc_queda30 = f(r -> r.queda_pct > 30),
            queda_mediana = q(:queda_pct, 0.5), queda_p25 = q(:queda_pct, 0.25), queda_p75 = q(:queda_pct, 0.75),
            t_nadir_mediano = q(:t_nadir, 0.5),
            bis_ok = f(r -> r.bis_min < 60),
            calibrados = f(r -> r.calibrou))
end
