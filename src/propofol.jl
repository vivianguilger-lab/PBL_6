# =============================================================================
# propofol.jl — PK/PD do propofol sobre o motor do NeoPKPD
# =============================================================================
#
# Segue src/propofol.jl do PropofolPBL do grupo
# (github.com/gabiadamek30-bit/Simula-o-da-resposta-do-paciente-ao-propofol---PBL),
# que usa o NeoPKPD (github.com/shramish2057/NeoPKPD) como motor:
#   - modelo tricompartimental (ThreeCompIVBolus) com bolus e infusões;
#   - compartimento de efeito (ke0) por simulate_biophase;
#   - BIS de Eleveld et al. (2018), calculado com o modelo de Eleveld completo
#     (PK + ke0 + atraso), independente do modelo PK escolhido para a hemodinâmica.
# Tempo em minutos, doses em mg, volumes em L ⇒ concentrações em mg/L = µg/ml.
# =============================================================================

"""
    parametros_pk(pac, modelo=:schnider; eta=NamedTuple())

Parâmetros PK do propofol `(CL, V1, Q2, V2, Q3, V3, ke0)` [L/min, L, 1/min].

- `:schnider` — PK de Schnider et al. 1998 e ke0 = 0,456/min de Schnider et al. 1999.
- `:marsh` — PK de Marsh et al. 1991. O ke0 = 0,26/min não é desse artigo: é o valor
  das bombas Diprifusor, nunca publicado em revista (Absalom et al. 2009).
- `:eleveld` — PK de Eleveld et al. 2018 (ver `parametros_eleveld`).

`eta` multiplica qualquer parâmetro por um fator (variabilidade interindividual,
usada só nas análises experimentais de incerteza e coorte).
"""
function parametros_pk(pac::Paciente, modelo::Symbol = :schnider; eta = NamedTuple())
    if modelo == :schnider
        lbm = massa_magra(pac)
        V1 = 4.27
        V2 = 18.9 - 0.391 * (pac.idade - 53)
        V3 = 238.0
        CL = 1.89 + 0.0456 * (pac.peso - 77) - 0.0681 * (lbm - 59) + 0.0264 * (pac.altura - 177)
        Q2 = 1.29 - 0.024 * (pac.idade - 53)
        Q3 = 0.836
        ke0 = 0.456
        # as equações lineares de Schnider ficam negativas fora da faixa de adultos
        # (ex.: V2 < 0 acima de 101 anos); nesse caso o modelo não se aplica
        min(V2, CL, Q2) > 0 ||
            throw(ArgumentError("paciente fora da faixa de validade do modelo de Schnider (V2, CL ou Q2 ≤ 0)"))
    elseif modelo == :marsh
        V1 = 0.228 * pac.peso
        CL = 0.119 * V1
        Q2 = 0.112 * V1
        V2 = Q2 / 0.055
        Q3 = 0.0419 * V1
        V3 = Q3 / 0.0033
        ke0 = 0.26
    elseif modelo == :eleveld
        e = parametros_eleveld(pac)
        CL, V1, Q2, V2, Q3, V3, ke0 = e.CL, e.V1, e.Q2, e.V2, e.Q3, e.V3, e.ke0
    else
        throw(ArgumentError("modelo PK desconhecido: $(modelo)"))
    end
    f(s) = Float64(get(eta, s, 1.0))
    return (CL = CL * f(:CL), V1 = V1 * f(:V1), Q2 = Q2 * f(:Q2), V2 = V2 * f(:V2),
            Q3 = Q3 * f(:Q3), V3 = V3 * f(:V3), ke0 = ke0 * f(:ke0), modelo = modelo)
end

"Bolus de `dose_mg` em t (min), administrado em `duracao` min (padrão 30 s)."
bolus(t, dose_mg; duracao = 0.5) = DoseEvent(Float64(t), Float64(dose_mg), Float64(duracao))

"Infusão contínua de `taxa_mg_min` de t_inicio a t_inicio + duracao (min)."
infusao(t_inicio, duracao, taxa_mg_min) =
    DoseEvent(Float64(t_inicio), Float64(taxa_mg_min * duracao), Float64(duracao))

"""
    esquema_inducao(pac; bolus_mgkg=2.0, duracao_bolus=0.5, manutencao_mcgkgmin=100.0,
                    inicio_manutencao=1.0, duracao_manutencao=30.0, t_bolus=0.0)

Esquema típico do grupo: bolus de indução (mg/kg) em `t_bolus` seguido de
manutenção (µg/kg/min) a partir de `t_bolus + inicio_manutencao`. Tempos em min.
"""
function esquema_inducao(pac::Paciente; bolus_mgkg = 2.0, duracao_bolus = 0.5, manutencao_mcgkgmin = 100.0,
                         inicio_manutencao = 1.0, duracao_manutencao = 30.0, t_bolus = 0.0)
    doses = DoseEvent[]
    bolus_mgkg > 0 && push!(doses, bolus(t_bolus, bolus_mgkg * pac.peso; duracao = duracao_bolus))
    if manutencao_mcgkgmin > 0 && duracao_manutencao > 0
        push!(doses, infusao(t_bolus + inicio_manutencao, duracao_manutencao, manutencao_mcgkgmin * pac.peso / 1000))
    end
    return doses
end

"""
    parametros_eleveld(pac; arterial=true, opioides=false)

PK do propofol de Eleveld et al. 2018 (Tabela 2) e ke0 do BIS (Tabela 3), na forma
do Python Anesthesia Simulator (pk_models.py, linhas 256–335, commit 4afc03a).
Massa livre de gordura de Al-Sallami; idade pós-menstrual = idade + 40 semanas.
`arterial = true` usa os parâmetros de amostras arteriais (ke0 = 0,146/min a 70 kg);
`opioides = false` porque o contexto de uso é a indução sem opioide.
"""
function parametros_eleveld(pac::Paciente; arterial = true, opioides = false)
    θ = (6.2830780766822, 25.5013145036879, 272.8166615043603, 1.7895836588902, 1.7500983738779,
         1.1085424008536, 0.191307, 42.2760190602615, 9.0548452392807, -0.015633, -0.00285709,
         33.5531248778544, -0.0138166, 68.2767978846832, 2.1002218877899, 1.3042680471360,
         1.4189043652084, 0.6805003109141)
    idade, peso, homem = pac.idade, pac.peso, pac.sexo == :M
    sig(x, c50, γ) = x^γ / (c50^γ + x^γ)
    fcentral(w) = sig(w, θ[12], 1)
    fenv(x) = exp(x * (idade - 35))
    al_sallami(h, w, a, imc) = h ? (0.88 + 0.12 / (1 + (a / 13.4)^-12.7)) * 9270 * w / (6680 + 216 * imc) :
                                   (1.11 - 0.11 / (1 + (a / 7.1)^-1.1)) * 9270 * w / (8780 + 244 * imc)
    pma, pma_ref = (idade + 40 / 52) * 52, (35 + 40 / 52) * 52          # semanas
    fCLmat, fCLmat_ref = sig(pma, θ[8], θ[9]), sig(pma_ref, θ[8], θ[9])
    fQ3mat, fQ3mat_ref = sig(pma, θ[14], 1), sig(pma_ref, θ[14], 1)
    ffm, ffm_ref = al_sallami(homem, peso, idade, imc(pac)), al_sallami(true, 70.0, 35.0, 70 / 1.7^2)
    fop(x) = opioides ? exp(x * idade) : 1.0
    V1 = θ[1] * fcentral(peso) / fcentral(70.0)
    arterial || (V1 *= 1 + θ[17] * (1 - fcentral(peso)))
    V2 = θ[2] * peso / 70 * fenv(θ[10])
    V3 = θ[3] * ffm / ffm_ref * fop(θ[13])
    CL = (homem ? θ[4] : θ[15]) * (peso / 70)^0.75 * fCLmat / fCLmat_ref * fop(θ[11])
    Q2 = θ[5] * (V2 / θ[2])^0.75 * (1 + θ[16] * (1 - fQ3mat))
    arterial || (Q2 *= θ[18])
    Q3 = θ[6] * (V3 / θ[3])^0.75 * fQ3mat / fQ3mat_ref
    ke0 = (arterial ? 0.146 : 1.24) * (peso / 70)^-0.25
    return (CL = CL, V1 = V1, Q2 = Q2, V2 = V2, Q3 = Q3, V3 = V3, ke0 = ke0)
end

"""
    bis_eleveld(Ce, idade)

BIS em função da concentração no sítio efetor do modelo de Eleveld et al. 2018
(Tabela 3): E0 = 93, Ce50 = 3,08·exp(−0,00635·(idade − 35)) µg/ml,
γ = 1,89 para Ce ≤ Ce50 e 1,47 acima.
"""
function bis_eleveld(Ce, idade)
    E0 = 93.0
    Ce50 = 3.08 * exp(-0.00635 * (idade - 35))
    γ = Ce <= Ce50 ? 1.89 : 1.47
    c = max(Float64(Ce), 0.0)
    return E0 - E0 * c^γ / (Ce50^γ + c^γ)
end

"Atraso do BIS [s] de Eleveld et al. 2018 (Tabela 3, θ = 0,0517): 15 + e^(0,0517·idade)."
atraso_bis(idade) = 15 + exp(0.0517 * idade)

# Cp e Ce na grade (min) para os parâmetros `pk`, pelo NeoPKPD
function _cp_ce(pk, doses, grade)
    spec = ModelSpec(ThreeCompIVBolus(), "propofol",
                     ThreeCompIVBolusParams(pk.CL, pk.V1, pk.Q2, pk.V2, pk.Q3, pk.V3),
                     sort(doses; by = d -> d.time))
    solver = SolverSpec(:Tsit5, 1e-8, 1e-10, 10^7)
    res = NeoPKPD.simulate(spec, SimGrid(grade[1], grade[end], grade), solver)
    # max(…, 0) só remove ruído numérico do integrador (valores da ordem de −1e−12)
    Cp = max.(res.observations[:conc], 0.0)
    Ce, _ = simulate_biophase(res.t, Cp, BiophaseEquilibrationParams(pk.ke0, 0.0, 1.0, 1.0), solver)
    return res.t, Cp, max.(Ce, 0.0)
end

"""
    simular_propofol(pac, doses; modelo=:schnider, tfim=60.0, dt=1/60, pk=parametros_pk(pac, modelo))

Simula o propofol com o NeoPKPD. Retorna `(t, Cp, Ce, bis, Ce_bis, pk)`, t em minutos
e concentrações em µg/ml:

- `Cp`, `Ce`: plasma e sítio efetor do modelo PK escolhido (`pk`), que movem a hemodinâmica;
- `bis`: BIS do modelo de Eleveld completo — PK e ke0 de Eleveld (`Ce_bis`), sigmoide
  de Eleveld e atraso de 15 + e^(0,0517·idade) s —, porque o Ce50 de 3,08 µg/ml só vale
  com o PK e o ke0 com que foi estimado.
"""
function simular_propofol(pac::Paciente, doses::Vector{DoseEvent}; modelo = :schnider,
                          tfim = 60.0, dt = 1 / 60, pk = parametros_pk(pac, modelo))
    grade = collect(0.0:dt:tfim)
    if isempty(doses)          # sem fármaco: tudo zero (NeoPKPD exige ao menos uma dose)
        z = zeros(length(grade))
        return (t = grade, Cp = z, Ce = copy(z), bis = fill(bis_eleveld(0.0, pac.idade), length(grade)),
                Ce_bis = copy(z), pk = pk)
    end
    t, Cp, Ce = _cp_ce(pk, doses, grade)
    _, _, Ce_bis = _cp_ce(parametros_pk(pac, :eleveld), doses, grade)
    bis0 = bis_eleveld.(Ce_bis, pac.idade)
    d = atraso_bis(pac.idade) / 60                  # min
    bis = [interp(t, bis0, max(ti - d, t[1])) for ti in t]
    return (t = t, Cp = Cp, Ce = Ce, bis = bis, Ce_bis = Ce_bis, pk = pk)
end

"Matriz K do sistema dA/dt = K·A (1/min), para o teste de verificação contra exp(K·t)."
function matriz_pk(pk)
    k10, k12, k21, k13, k31 = pk.CL / pk.V1, pk.Q2 / pk.V1, pk.Q2 / pk.V2, pk.Q3 / pk.V1, pk.Q3 / pk.V3
    return [-(k10 + k12 + k13) k21 k31; k12 -k21 0.0; k13 0.0 -k31]
end
