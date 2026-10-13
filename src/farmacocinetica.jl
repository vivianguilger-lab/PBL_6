# =============================================================================
# farmacocinetica.jl — PK do propofol: modelo tricompartimental de Schnider
#                      + dois compartimentos de efeito (cérebro e vasos)
# =============================================================================
#
#            infusão u(t) [mg/s]
#                  │
#                  ▼
#   ┌────────┐  k12  ┌────────┐            ┌────────┐
#   │  V2    │◄──────│  V1    │──────────►│  V3    │
#   │ rápido │──────►│central │◄──────────│ lento  │
#   └────────┘  k21  └───┬────┘  k13/k31   └────────┘
#                        │ k10 (eliminação)
#                        ▼
#        Cp = A1/V1 ──ke0_bis──► Ce_bis  (efeito hipnótico / SNC)
#                   ──ke0_hemo─► Ce_hemo (efeito hemodinâmico, mais lento)
#
# Referências
#   Schnider TW et al. Anesthesiology 1998;88:1170-82 (PK) e 1999;90:1502-16 (ke0).
#   Kazama T et al. Anesthesiology 1999;90:1517-27: o efeito do propofol sobre a
#   pressão arterial equilibra-se MAIS DEVAGAR que o BIS (t½ke0 ≈ 5,7 min em
#   jovens e ≈ 10,2 min em idosos, contra ≈ 2,3 min para o BIS).
# =============================================================================

"""
Parâmetros farmacocinéticos de um paciente. Volumes em L, constantes em 1/s.
"""
struct ParametrosPK
    V1::Float64
    V2::Float64
    V3::Float64
    Cl1::Float64      # L/min (guardado só para relatório)
    Cl2::Float64
    Cl3::Float64
    k10::Float64
    k12::Float64
    k21::Float64
    k13::Float64
    k31::Float64
    ke0_bis::Float64
    ke0_hemo::Float64
end

"""
    meia_vida_ke0_hemodinamica(idade) -> minutos

t½ do equilíbrio plasma ↔ local de efeito hemodinâmico (Kazama 1999):
5,7 min no grupo 20–39 anos (≈30 anos) e 10,2 min no grupo 70–85 anos
(≈77 anos); interpolação linear entre os grupos e extrapolação constante fora.
"""
meia_vida_ke0_hemodinamica(idade) = 5.7 + (10.2 - 5.7) * clamp((idade - 30) / (77 - 30), 0.0, 1.0)

"""
    schnider(p::Paciente; eta = NamedTuple()) -> ParametrosPK

Modelo de Schnider para o propofol. `eta` permite multiplicar qualquer
parâmetro (V1, V2, V3, Cl1, Cl2, Cl3, ke0_bis, ke0_hemo) por um fator — é assim
que a coorte virtual introduz variabilidade interindividual.

    V1  = 4,27 L
    V2  = 18,9 − 0,391·(idade − 53) L
    V3  = 238 L
    Cl1 = 1,89 + 0,0456·(peso − 77) − 0,0681·(LBM − 59) + 0,0264·(altura − 177) L/min
    Cl2 = 1,29 − 0,024·(idade − 53) L/min
    Cl3 = 0,836 L/min
    ke0 = 0,456 /min
"""
function schnider(p::Paciente; eta = NamedTuple())
    f(s) = get(eta, s, 1.0)
    lbm = massa_magra(p)
    V1 = 4.27 * f(:V1)
    V2 = max(18.9 - 0.391 * (p.idade - 53), 1.0) * f(:V2)
    V3 = 238.0 * f(:V3)
    Cl1 = max(1.89 + 0.0456 * (p.peso - 77) - 0.0681 * (lbm - 59) + 0.0264 * (p.altura - 177), 0.2) * f(:Cl1)
    Cl2 = max(1.29 - 0.024 * (p.idade - 53), 0.1) * f(:Cl2)
    Cl3 = 0.836 * f(:Cl3)
    ke0_bis = 0.456 * f(:ke0_bis)
    ke0_hemo = log(2) / meia_vida_ke0_hemodinamica(p.idade) * f(:ke0_hemo)
    m = 1 / 60  # 1/min → 1/s
    return ParametrosPK(V1, V2, V3, Cl1, Cl2, Cl3,
                        Cl1 / V1 * m, Cl2 / V1 * m, Cl2 / V2 * m, Cl3 / V1 * m, Cl3 / V3 * m,
                        ke0_bis * m, ke0_hemo * m)
end

"""
Matriz do sistema linear dA/dt = K·A + b·u para os 3 compartimentos (1/s).
Usada nos testes de verificação (solução exata via exponencial de matriz).
"""
function matriz_pk(pk::ParametrosPK)
    return [-(pk.k10 + pk.k12 + pk.k13)  pk.k21   pk.k31;
             pk.k12                     -pk.k21   0.0;
             pk.k13                      0.0     -pk.k31]
end

# --------------------------------------------------------------------------
# Esquema de administração
# --------------------------------------------------------------------------

"""
Um trecho de infusão com taxa constante `taxa` [mg/s] entre `inicio` e `fim` [s].
Um bolus é simplesmente um trecho curto (ex.: 20 s).
"""
struct TrechoInfusao
    inicio::Float64
    fim::Float64
    taxa::Float64
end

"""
    protocolo_inducao(p; dose_mgkg=2.0, duracao_bolus=20.0, t_inicio=60.0,
                      manutencao_ugkgmin=0.0, t_final=21*60)

Monta o esquema de administração: bolus de `dose_mgkg` (mg/kg) injetado em
`duracao_bolus` segundos a partir de `t_inicio`, seguido opcionalmente de
infusão de manutenção em µg/kg/min até `t_final`.
"""
function protocolo_inducao(p::Paciente; dose_mgkg = 2.0, duracao_bolus = 20.0, t_inicio = 60.0,
                           manutencao_ugkgmin = 0.0, t_final = 21 * 60.0)
    trechos = TrechoInfusao[]
    dose_mg = dose_mgkg * p.peso
    dose_mg > 0 && push!(trechos, TrechoInfusao(t_inicio, t_inicio + duracao_bolus, dose_mg / duracao_bolus))
    if manutencao_ugkgmin > 0
        taxa = manutencao_ugkgmin * p.peso / 1000 / 60      # µg/kg/min → mg/s
        push!(trechos, TrechoInfusao(t_inicio + duracao_bolus, t_final, taxa))
    end
    return trechos
end

@inline function taxa_infusao(trechos::Vector{TrechoInfusao}, t)
    s = 0.0
    @inbounds for tr in trechos
        (tr.inicio <= t < tr.fim) && (s += tr.taxa)
    end
    return s
end

"""Instantes em que a taxa muda (o integrador deve parar exatamente neles)."""
instantes_infusao(trechos) = sort!(unique!(vcat([tr.inicio for tr in trechos], [tr.fim for tr in trechos])))

"""Massa total administrada [mg]."""
massa_total(trechos) = sum((tr.fim - tr.inicio) * tr.taxa for tr in trechos; init = 0.0)
