# =============================================================================
# guyton.jl — leitura do modelo pelo "ponto de funcionamento" do curso
# =============================================================================
#
# Resumo, Seção 16: fechado o laço, restam duas retas no plano (RAP, Q):
#
#   (A) retorno venoso   VR = (Pms − RAP) / R_VR,   R_VR = Rv + TPR·AC/(AC+VC)
#   (B) débito cardíaco  CO = FC·(C_D·RAP − C_S·Pa)
#
# e o ponto de funcionamento é a interseção CO = VR. Aqui usamos os
# parâmetros EFETIVOS do paciente virtual (antes e no nadir do propofol)
# para mostrar QUAL mecanismo derruba a PAM. É uma leitura estática e
# linearizada (o modelo pulsátil é a referência quantitativa).
# =============================================================================

"""
    pressao_media_enchimento(cv; Vv0=cv.Vv0, Emin=cv.Emin)

Pms = volume estressado / complacência total (Resumo, eq. 8), incluindo o
ventrículo relaxado (C_D = 1/Emin) e seu volume morto:

    Pms = (Vtotal − Va0 − V0v − Vd) / (AC + VC + C_D)
"""
pressao_media_enchimento(cv::ParametrosCV; Vv0 = cv.Vv0, Emin = cv.Emin) =
    (cv.Vtotal - cv.Va0 - Vv0 - cv.Vd) / (cv.Ca + cv.Cv + 1 / Emin)

"""Resistência ao retorno venoso R_VR = Rv + TPR·AC/(AC+VC) (Resumo, eq. 12)."""
resistencia_retorno_venoso(Rv, TPR, AC, VC) = Rv + TPR * AC / (AC + VC)

"""Curva de retorno venoso [L/min] em função de RAP [mmHg]."""
retorno_venoso(RAP, Pms, RVR) = max((Pms - RAP) / RVR, 0.0) * 60 / 1000

"""Curva do débito cardíaco [L/min] (Resumo, eq. 14), com piso em zero."""
debito_cardiaco_curso(RAP, FC, CD, CS, Pa) = max(FC * (CD * RAP - CS * Pa), 0.0) / 1000

"""
    estado_guyton(cv; FC, Rs, Emax, Vv0, Pa)

Parâmetros das duas retas e o ponto de funcionamento (RAP*, Q*) para um
conjunto de parâmetros efetivos.
"""
function estado_guyton(cv::ParametrosCV; FC = cv.FC, Rs = cv.Rs, Emax = cv.Emax, Vv0 = cv.Vv0, Pa = 90.0)
    Pms = pressao_media_enchimento(cv; Vv0)
    RVR = resistencia_retorno_venoso(cv.Rin, Rs, cv.Ca, cv.Cv)
    CD, CS = 1 / cv.Emin, 1 / Emax
    # interseção: (Pms − R)/RVR = (FC/60)(CD·R − CS·Pa)   [ml/s]
    h = FC / 60
    RAP = (Pms / RVR + h * CS * Pa) / (1 / RVR + h * CD)
    Q = (Pms - RAP) / RVR * 60 / 1000
    return (Pms = Pms, RVR = RVR, CD = CD, CS = CS, FC = FC, Pa = Pa, RAP = RAP, Q = Q)
end
