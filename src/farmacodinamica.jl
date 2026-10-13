# =============================================================================
# farmacodinamica.jl — PD do propofol: equação de Hill (Emax sigmoide)
# =============================================================================
#
#                       Ce^γ
#   E(Ce) = Emax · ─────────────        (Hill, 1910)
#                  C50^γ + Ce^γ
#
# Cada efeito hemodinâmico do propofol descrito no problema vira um parâmetro
# do circuito Windkessel modulado por uma curva de Hill:
#
#   efeito clínico                    parâmetro do modelo      fonte dos valores
#   ─────────────────────────────────  ──────────────────────   ───────────────────────
#   vasodilatação arteriolar           TPR (resistor R)  ↓      Su et al. 2023 (BJA)
#   depressão miocárdica               Emax = 1/C_S      ↓      Su et al. 2023 (efeito
#                                      (inotropismo)            no VS, dependente da idade)
#   venodilatação                      V0 venoso         ↑      de Wit et al. 2016 (BJA)
#                                      (↓ volume estressado,     — magnitude ASSUMIDA
#                                       ↓ Pms)
#   atenuação do barorreflexo          ganho do reflexo  ↓      Cullen 1987; Ebert 1992
#                                                               — magnitude ASSUMIDA
#   hipnose (alvo terapêutico)         BIS               ↓      Bouillon et al. 2004
#
# Efeitos vasculares usam Ce_hemo (equilíbrio lento, Kazama 1999); hipnose e
# atenuação do barorreflexo (mediada pelo SNC) usam Ce_bis.
# =============================================================================

"""
    hill(C, C50, γ)

Fração do efeito máximo (0…1) na concentração `C`.
"""
@inline hill(C, C50, γ) = C <= 0 ? 0.0 : (r = (C / C50)^γ; r / (1 + r))

"""
Parâmetros farmacodinâmicos (concentrações em µg/ml).
Os valores padrão e suas fontes estão documentados em `docs/modelo_matematico.md`.
"""
Base.@kwdef struct ParametrosPD
    # Hipnose — BIS = E0 − E0·hill(Ce_bis)          (Bouillon 2004, sem opioide)
    bis_E0::Float64 = 97.4
    bis_C50::Float64 = 4.47
    bis_gamma::Float64 = 1.43
    # Resistência vascular periférica (TPR)          (Su 2023)
    tpr_Emax::Float64 = -0.778
    tpr_C50::Float64 = 3.21
    tpr_gamma::Float64 = 1.83
    # Inotropismo (elastância sistólica máxima)      (Su 2023, efeito no VS)
    #   Emax(idade) = −0,154 · exp(0,0333·(idade − 35))
    inot_Emax35::Float64 = -0.154
    inot_coef_idade::Float64 = 0.0333
    inot_C50::Float64 = 0.44
    inot_gamma::Float64 = 1.0
    # Venodilatação (aumento do volume não estressado venoso) — ASSUMIDO
    veno_Emax::Float64 = 0.25
    veno_C50::Float64 = 2.5
    veno_gamma::Float64 = 1.5
    # Atenuação do barorreflexo (fração do ganho perdida) — ASSUMIDO
    baro_Emax::Float64 = 0.8
    baro_C50::Float64 = 2.5
    baro_gamma::Float64 = 2.0
end

"""Efeito inotrópico máximo para a idade (limitado a −90%)."""
inot_Emax(pd::ParametrosPD, idade) = max(pd.inot_Emax35 * exp(pd.inot_coef_idade * (idade - 35)), -0.9)

"""
    efeitos(pd, idade, Ce_bis, Ce_hemo)

Variações fracionárias dos parâmetros do circuito causadas pelo propofol:
`tpr`, `inot` (Emax do ventrículo), `veno` (V0 venoso) e `baro` (ganho do
barorreflexo), além do `bis`.
"""
@inline function efeitos(pd::ParametrosPD, inotmax, Ce_bis, Ce_hemo)
    return (tpr  = pd.tpr_Emax * hill(Ce_hemo, pd.tpr_C50, pd.tpr_gamma),
            inot = inotmax * hill(Ce_hemo, pd.inot_C50, pd.inot_gamma),
            veno = pd.veno_Emax * hill(Ce_hemo, pd.veno_C50, pd.veno_gamma),
            baro = -pd.baro_Emax * hill(Ce_bis, pd.baro_C50, pd.baro_gamma),
            bis  = pd.bis_E0 * (1 - hill(Ce_bis, pd.bis_C50, pd.bis_gamma)))
end
