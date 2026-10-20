# =============================================================================
# propofol.jl — PK/PD do propofol sobre o motor do NeoPKPD
# =============================================================================
#
# Segue src/propofol.jl do PropofolPBL do grupo
# (github.com/gabiadamek30-bit/Simula-o-da-resposta-do-paciente-ao-propofol---PBL),
# que usa o NeoPKPD (github.com/shramish2057/NeoPKPD) como motor:
#   - modelo tricompartimental (ThreeCompIVBolus) com bolus e infusões;
#   - compartimento de efeito (ke0) por simulate_biophase;
#   - BIS de Eleveld et al. (2018).
# Tempo em minutos, doses em mg, volumes em L ⇒ concentrações em mg/L = µg/ml.
# =============================================================================

"""
    parametros_pk(pac, modelo=:schnider; eta=NamedTuple())

Parâmetros PK do propofol `(CL, V1, Q2, V2, Q3, V3, ke0)` [L/min, L, 1/min].
`modelo = :schnider` (Schnider et al. 1998/1999) ou `:marsh` (Marsh et al. 1991,
ke0 = 0,26/min). `eta` multiplica qualquer parâmetro por um fator (variabilidade
interindividual usada na faixa de incerteza e na coorte virtual).
"""
function parametros_pk(pac::Paciente, modelo::Symbol = :schnider; eta = NamedTuple())
    if modelo == :schnider
        lbm = massa_magra(pac)
        V1 = 4.27
        V2 = max(18.9 - 0.391 * (pac.idade - 53), 1.0)
        V3 = 238.0
        CL = max(1.89 + 0.0456 * (pac.peso - 77) - 0.0681 * (lbm - 59) + 0.0264 * (pac.altura - 177), 0.2)
        Q2 = max(1.29 - 0.024 * (pac.idade - 53), 0.1)
        Q3 = 0.836
        ke0 = 0.456
    elseif modelo == :marsh
        V1 = 0.228 * pac.peso
        CL = 0.119 * V1
        Q2 = 0.112 * V1
        V2 = Q2 / 0.055
        Q3 = 0.0419 * V1
        V3 = Q3 / 0.0033
        ke0 = 0.26
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
    bis_eleveld(Ce, idade)

BIS em função da concentração no sítio efetor (sigmoide de Eleveld et al. 2018):
E0 = 93, Ce50 = 3,08·exp(−0,00635·(idade − 35)) µg/ml, γ = 1,47 abaixo de Ce50
e 1,89 acima. O Ce vem do ke0 do modelo PK escolhido (aproximação do grupo).
"""
function bis_eleveld(Ce, idade)
    E0 = 93.0
    Ce50 = 3.08 * exp(-0.00635 * (idade - 35))
    γ = Ce <= Ce50 ? 1.47 : 1.89
    c = max(Float64(Ce), 0.0)
    return E0 - E0 * c^γ / (Ce50^γ + c^γ)
end

"""
    simular_propofol(pac, doses; modelo=:schnider, tfim=60.0, dt=1/60, pk=parametros_pk(pac, modelo))

Simula a PK/PD do propofol com o NeoPKPD. Retorna `(t, Cp, Ce, bis, pk)` com
t em minutos e concentrações em µg/ml.
"""
function simular_propofol(pac::Paciente, doses::Vector{DoseEvent}; modelo = :schnider,
                          tfim = 60.0, dt = 1 / 60, pk = parametros_pk(pac, modelo))
    grade = collect(0.0:dt:tfim)
    if isempty(doses)          # sem fármaco: tudo zero (NeoPKPD exige ao menos uma dose)
        z = zeros(length(grade))
        return (t = grade, Cp = z, Ce = copy(z), bis = fill(bis_eleveld(0.0, pac.idade), length(grade)), pk = pk)
    end
    spec = ModelSpec(ThreeCompIVBolus(), "propofol_$(pk.modelo)",
                     ThreeCompIVBolusParams(pk.CL, pk.V1, pk.Q2, pk.V2, pk.Q3, pk.V3),
                     sort(doses; by = d -> d.time))
    grid = SimGrid(0.0, Float64(tfim), grade)
    solver = SolverSpec(:Tsit5, 1e-8, 1e-10, 10^7)
    res = NeoPKPD.simulate(spec, grid, solver)
    Cp = max.(res.observations[:conc], 0.0)
    Ce, _ = simulate_biophase(res.t, Cp, BiophaseEquilibrationParams(pk.ke0, 0.0, 1.0, 1.0), solver)
    Ce = max.(Ce, 0.0)
    return (t = res.t, Cp = Cp, Ce = Ce, bis = bis_eleveld.(Ce, pac.idade), pk = pk)
end

"Matriz K do sistema dA/dt = K·A (1/min), para o teste de verificação contra exp(K·t)."
function matriz_pk(pk)
    k10, k12, k21, k13, k31 = pk.CL / pk.V1, pk.Q2 / pk.V1, pk.Q2 / pk.V2, pk.Q3 / pk.V1, pk.Q3 / pk.V3
    return [-(k10 + k12 + k13) k21 k31; k12 -k21 0.0; k13 0.0 -k31]
end
