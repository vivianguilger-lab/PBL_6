# =============================================================================
# windkessel.jl — circulação em malha fechada (modelo 0-D, parâmetros concentrados)
# =============================================================================
#
#  Analogia eletromecânica (Resumo, Seções 6 e 10–16):
#     pressão ↔ tensão (across)   fluxo ↔ corrente (through)
#     resistência ↔ resistor      complacência ↔ capacitor     válvula ↔ diodo
#     ventrículo ↔ capacitor de capacitância variável C(t) = 1/E(t)
#
#           Rin (diodo)           Rav+Zc (diodo)
#   ┌──────────►|──────── VE ─────────►|────────┬─────────┐ Pao (PA medida)
#   │                    E(t)                   │         │
#   │                                       Ca ═╪═        ▓ TPR = Rs
#   │                                           │ Pc      │
#   │   Pv                                      ┴         │
#   └───────────────────── Cv (veias) ◄────────────────────┘
#                     (V0v = volume não estressado)
#
#  - Árvore arterial: Windkessel de 3 elementos (Zc, Ca = AC, Rs = TPR).
#  - Árvore venosa: complacência Cv = VC com volume não estressado V0v.
#  - Ventrículo: elastância variável E(t) = Emin + (Emax − Emin)·e(t)
#      ⇒ C_D = 1/Emin (diástole) e C_S = 1/Emax (sístole), exatamente o
#        "capacitor de capacitância variável" do Resumo (Seção 14.3).
#    O coração direito e a circulação pulmonar estão agregados ao ventrículo
#    e a Rin (simplificação de um único ventrículo, como no modelo de
#    Bjørdalsbakke et al. 2022 usado no CirculatorySystemModels.jl).
#  - Barorreflexo: malha fechada de controle (Resumo, Seção 9) atuando em
#    FC, TPR, Emax e V0v, com sensibilidade reduzida pelo propofol.
# =============================================================================

"""
Parâmetros mecânicos do sistema cardiovascular. Unidades: mmHg, ml, s.
Os valores padrão representam um adulto típico e são substituídos pela
calibração do paciente virtual (`calibrar`).
"""
Base.@kwdef struct ParametrosCV
    FC::Float64 = 72.0       # bpm
    Emax::Float64 = 2.0      # mmHg/ml   (C_S = 0,5 ml/mmHg, valor do Resumo)
    Emin::Float64 = 0.06     # mmHg/ml   (C_D ≈ 17 ml/mmHg)
    Vd::Float64 = 10.0       # ml        volume morto do ventrículo
    Rin::Float64 = 0.02      # mmHg·s/ml resistência de enchimento (veias centrais + valva AV)
    Rav::Float64 = 0.005     # mmHg·s/ml resistência da valva aórtica aberta
    Zc::Float64 = 0.03       # mmHg·s/ml impedância característica da aorta
    Ca::Float64 = 1.5        # ml/mmHg   complacência arterial (AC)
    Rs::Float64 = 1.0        # mmHg·s/ml resistência periférica total (TPR)
    Va0::Float64 = 700.0     # ml        volume não estressado arterial
    Cv::Float64 = 100.0      # ml/mmHg   complacência venosa (VC, valor do Resumo)
    Vv0::Float64 = 2500.0    # ml        volume não estressado venoso
    Vtotal::Float64 = 4250.0 # ml        volume circulante do modelo
end

"""
Barorreflexo arterial simplificado (inspirado em Ursino 1998).

    r = tanh(−(Ps − Pset)/S)              resposta aferente normalizada
    τᵢ·dxᵢ/dt = g·r − xᵢ                  dinâmica de cada efetor
    FC  = FC₀ ·(1 + A_FC·x_FC)
    TPR = TPR₀·(1 + A_R·x_R)
    Emax= Emax₀·(1 + A_E·x_E)
    V0v = V0v₀·(1 − A_V·x_V)              venoconstrição recruta volume

`g = 1 + efeito_baro(Ce_bis)` é o ganho, reduzido pelo propofol.
"""
Base.@kwdef struct ParametrosBaro
    Pset::Float64 = 93.0     # mmHg  ponto de ajuste (= PAM basal do paciente)
    S::Float64 = 25.0        # mmHg  escala da sigmoide
    tauS::Float64 = 2.0      # s     filtro do sensor (barorreceptor)
    A_FC::Float64 = 0.40
    A_R::Float64 = 0.35
    A_E::Float64 = 0.25
    A_V::Float64 = 0.08
    tauFC::Float64 = 2.0
    tauR::Float64 = 6.0
    tauE::Float64 = 8.0
    tauV::Float64 = 20.0
end

"""
Modelo de observação da pressão arterial (equação de medida).

O modelo 0-D fornece a pressão na raiz da aorta. O manguito mede no braço,
onde a pressão de pulso é AMPLIFICADA pela reflexão de ondas — fenômeno que
um modelo 0-D não reproduz (Resumo, Seção 12.4). A PAM praticamente não muda
ao longo da árvore, então:

    PP_braquial = amp · PP_aórtica
    PAD_braquial = PAM − ff · PP_braquial     (ff = (PAM−PAD)/(PAS−PAD) do manguito; 1/3 por padrão)
    PAS_braquial = PAD_braquial + PP_braquial
"""
Base.@kwdef struct ObservacaoPA
    amp::Float64 = 1.0
    ff::Float64 = 1 / 3
end

"""
    amplificacao_pp(idade)

Razão PP braquial / PP aórtica. A amplificação cai com a idade porque a
aorta enrijece e as ondas refletidas retornam mais cedo (somando-se à
sístole já na aorta). Aproximação exponencial da ordem de grandeza
relatada no estudo Anglo-Cardiff (McEniery et al., JACC 2005):
≈1,7 aos 20 anos, ≈1,3 aos 50 e ≈1,25 após os 65.
"""
amplificacao_pp(idade) = 1.25 + 0.45 * exp(-(max(idade, 18.0) - 20) / 15)

"""
Modelo completo (PK + PD + circulação + barorreflexo) — é o `p` passado
ao integrador.
"""
struct ModeloInducao
    pk::ParametrosPK
    pd::ParametrosPD
    cv::ParametrosCV
    baro::ParametrosBaro
    inotmax::Float64                 # efeito inotrópico máximo para a idade
    infusao::Vector{TrechoInfusao}
    barorreflexo::Bool               # liga/desliga a malha fechada
    farmaco::Bool                    # liga/desliga os efeitos do propofol
    obs::ObservacaoPA                # conversão aorta → manguito (só para saída)
end

ModeloInducao(pk, pd, cv, baro, inotmax, infusao, barorreflexo, farmaco) =
    ModeloInducao(pk, pd, cv, baro, inotmax, infusao, barorreflexo, farmaco, ObservacaoPA())

"""Cópia de um struct trocando alguns campos: `com(cv; Rs = 0.8)`."""
function com(x::T; kw...) where {T}
    vals = map(f -> haskey(kw, f) ? kw[f] : getfield(x, f), fieldnames(T))
    return T(vals...)
end

# Índices do vetor de estado ---------------------------------------------------
const NESTADOS = 16
const I_A1, I_A2, I_A3, I_CEB, I_CEH = 1, 2, 3, 4, 5      # PK [mg] e Ce [µg/ml]
const I_FASE, I_VLV, I_VA, I_VV = 6, 7, 8, 9               # fase [batimentos], volumes [ml]
const I_PS, I_XFC, I_XR, I_XE, I_XV = 10, 11, 12, 13, 14  # barorreflexo
const I_INTP, I_INTQ = 15, 16                              # ∫Pao dt, ∫Qao dt

"""
    ativacao(tc, T)

Função de ativação ventricular e(t) ∈ [0,1] (cossenoidal, à la Shi et al.):
contração em Ts = 0,35·√T, relaxamento em Ts/2, repouso no resto do ciclo.
"""
@inline function ativacao(tc, T)
    Ts = 0.35 * sqrt(T)
    Tr = 0.5 * Ts
    if tc < Ts
        return 0.5 * (1 - cos(pi * tc / Ts))
    elseif tc < Ts + Tr
        return 0.5 * (1 + cos(pi * (tc - Ts) / Tr))
    else
        return 0.0
    end
end

"""
    algebrico(u, m, t)

Todas as variáveis algébricas (pressões, fluxos, parâmetros efetivos) a
partir do estado `u`. Usada pela EDO e pelo pós-processamento — assim o que
se plota é exatamente o que foi integrado.
"""
@inline function algebrico(u, m::ModeloInducao, t)
    cv, b, pk = m.cv, m.baro, m.pk
    Cp = u[I_A1] / pk.V1                          # mg/L = µg/ml
    ef = efeitos(m.pd, m.inotmax, u[I_CEB], u[I_CEH])
    fat = m.farmaco ? 1.0 : 0.0
    FC = cv.FC * (1 + b.A_FC * u[I_XFC])
    Rs = cv.Rs * (1 + fat * ef.tpr) * (1 + b.A_R * u[I_XR])
    Emax = max(cv.Emax * (1 + fat * ef.inot) * (1 + b.A_E * u[I_XE]), cv.Emin)
    Vv0 = cv.Vv0 * (1 + fat * ef.veno) * (1 - b.A_V * u[I_XV])
    T = 60 / FC
    fase = u[I_FASE]
    tc = (fase - floor(fase)) * T
    E = cv.Emin + (Emax - cv.Emin) * ativacao(tc, T)
    Plv = E * (u[I_VLV] - cv.Vd)
    Pc = (u[I_VA] - cv.Va0) / cv.Ca
    Pv = (u[I_VV] - Vv0) / cv.Cv
    Qin = max(Pv - Plv, 0.0) / cv.Rin
    Qao = max(Plv - Pc, 0.0) / (cv.Rav + cv.Zc)
    Pao = Pc + cv.Zc * Qao
    Qs = (Pc - Pv) / Rs
    return (; Cp, ef, FC, Rs, Emax, Vv0, E, Plv, Pc, Pv, Pao, Qin, Qao, Qs)
end

"""Lado direito da EDO (forma in-place exigida pelo OrdinaryDiffEq)."""
function rhs!(du, u, m::ModeloInducao, t)
    pk, b = m.pk, m.baro
    a = algebrico(u, m, t)
    # --- Farmacocinética (Schnider) ---
    inf = taxa_infusao(m.infusao, t)
    A1, A2, A3 = u[I_A1], u[I_A2], u[I_A3]
    du[I_A1] = inf - (pk.k10 + pk.k12 + pk.k13) * A1 + pk.k21 * A2 + pk.k31 * A3
    du[I_A2] = pk.k12 * A1 - pk.k21 * A2
    du[I_A3] = pk.k13 * A1 - pk.k31 * A3
    du[I_CEB] = pk.ke0_bis * (a.Cp - u[I_CEB])
    du[I_CEH] = pk.ke0_hemo * (a.Cp - u[I_CEH])
    # --- Circulação (conservação de volume em cada nó) ---
    du[I_FASE] = a.FC / 60
    du[I_VLV] = a.Qin - a.Qao
    du[I_VA] = a.Qao - a.Qs
    du[I_VV] = a.Qs - a.Qin
    # --- Barorreflexo (malha fechada) ---
    du[I_PS] = (a.Pao - u[I_PS]) / b.tauS
    g = m.barorreflexo ? (1 + (m.farmaco ? a.ef.baro : 0.0)) : 0.0
    r = tanh(-(u[I_PS] - b.Pset) / b.S)
    du[I_XFC] = (g * r - u[I_XFC]) / b.tauFC
    du[I_XR] = (g * r - u[I_XR]) / b.tauR
    du[I_XE] = (g * r - u[I_XE]) / b.tauE
    du[I_XV] = (g * r - u[I_XV]) / b.tauV
    # --- Integrais para médias exatas por batimento ---
    du[I_INTP] = a.Pao
    du[I_INTQ] = a.Qao
    return nothing
end

"""
Estado inicial "ingênuo": volumes distribuídos para uma pressão uniforme
`P0` (≈ Pms) e farmacocinética zerada.
"""
function estado_inicial(cv::ParametrosCV; P0 = 8.0)
    u = zeros(NESTADOS)
    Vlv = cv.Vd + P0 / cv.Emin
    Va = cv.Va0 + cv.Ca * P0
    u[I_VLV] = Vlv
    u[I_VA] = Va
    u[I_VV] = cv.Vtotal - Vlv - Va
    u[I_PS] = P0
    return u
end

"""Volume total contido no circuito (deve ser constante — teste de verificação)."""
volume_total(u) = u[I_VLV] + u[I_VA] + u[I_VV]
