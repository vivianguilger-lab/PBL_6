# =============================================================================
# calibracao.jl — paciente virtual calibrado com os dados pré-indução
# =============================================================================
#
# Abordagem grey-box (Resumo, Seção 5.2): a ESTRUTURA vem da física
# (Windkessel + coração elastância + veias); os NÚMEROS vêm do paciente.
#
# Ideia: as equações do curso dão a primeira estimativa e um ajuste
# iterativo (ponto fixo) corrige o resultado do modelo pulsátil completo.
#
#   TPR  ≈ (PAM − Pv) / DC          Lei de Ohm hemodinâmica (Seção 12.3)
#   AC   ≈ VS / PP                  complacência = ΔV/ΔP (Seção 6.2)
#   V0v  ↔ DC                       o volume estressado define Pms e o retorno
#                                   venoso (Seções 13.4–13.5)
#   τ    = TPR · AC                 constante de tempo diastólica (Seção 12.2)
# =============================================================================

"""
Paciente virtual: parâmetros individualizados + estado de regime permanente.
"""
struct PacienteVirtual
    paciente::Paciente
    pk::ParametrosPK
    pd::ParametrosPD
    cv::ParametrosCV
    baro::ParametrosBaro
    obs::ObservacaoPA
    u0::Vector{Float64}
    alvo::NamedTuple
    obtido::NamedTuple
    iteracoes::Int
    convergiu::Bool
end

# limites fisiológicos usados na identificação
const CA_MIN, CA_MAX = 0.3, 4.0       # ml/mmHg
const ZC_MIN, ZC_MAX = 0.005, 0.12    # mmHg·s/ml

"""Constante de tempo diastólica τ = TPR·AC [s] do paciente virtual."""
tau_diastolico(pv::PacienteVirtual) = pv.cv.Rs * pv.cv.Ca

"""
Roda `nbat` batimentos com o coração em frequência fixa (sem fármaco e sem
barorreflexo) a partir de `u` e devolve (médias dos 3 últimos batimentos,
estado final exatamente no início de um batimento).
"""
function _regime(cv::ParametrosCV, baro::ParametrosBaro, pk, pd, u, nbat)
    m = ModeloInducao(pk, pd, cv, baro, 0.0, TrechoInfusao[], false, false)
    u = copy(u)
    u[I_FASE] = 0.0; u[I_INTP] = 0.0; u[I_INTQ] = 0.0
    T = 60 / cv.FC
    sol = integrar(m, u, (0.0, nbat * T); saveat = T / 400)
    bat = batimentos(sol.t, sol[I_FASE, :], sol[I_INTP, :], sol[I_INTQ, :],
                     [algebrico(x, m, 0.0).Pao for x in sol.u])
    k = length(bat.t)
    ult = max(1, k - 2):k
    med = (PAM = mean(bat.PAM[ult]), PAS = mean(bat.PASao[ult]), PAD = mean(bat.PADao[ult]),
           DC = mean(bat.DC[ult]), VS = mean(bat.VS[ult]))
    return med, sol.u[end]
end

"""
    calibrar(p::Paciente; pd=ParametrosPD(), eta_pk=NamedTuple(),
             cv_ref=ParametrosCV(), baro_ref=ParametrosBaro(),
             tol=0.005, maxit=60, verbose=false) -> PacienteVirtual

Constrói o paciente virtual: ajusta TPR (Rs), complacência arterial (Ca) e
volume não estressado venoso (V0v) para que o modelo reproduza a PAM, a
pressão de pulso AÓRTICA (PP do manguito ÷ amplificação periférica) e o
débito cardíaco do paciente na frequência cardíaca medida. Parâmetros de "tamanho" (elastâncias, volumes) são escalados pela
superfície corporal e pela volemia de Nadler.
"""
function calibrar(p::Paciente; pd = ParametrosPD(), eta_pk = NamedTuple(),
                  cv_ref = ParametrosCV(), baro_ref = ParametrosBaro(),
                  tol = 0.005, maxit = 80, verbose = false)
    pk = schnider(p; eta = eta_pk)
    s = superficie_corporal(p) / 1.9                 # escala alométrica simples
    Vtot = 0.85 * volemia(p)                          # ~85% da volemia na circulação sistêmica (Resumo 13.2)
    amp = amplificacao_pp(p.idade)
    obs = ObservacaoPA(amp = amp, ff = (pam_alvo(p) - p.pad) / (p.pas - p.pad))
    alvo = (PAM = pam_alvo(p), PP = (p.pas - p.pad) / amp, DC = dc_alvo(p), FC = p.fc, PAS = p.pas, PAD = p.pad)
    Q = alvo.DC * 1000 / 60                           # ml/s
    VS = Q * 60 / p.fc
    # --- primeira estimativa pelas equações do curso ---
    Rs = (alvo.PAM - 8.0) / Q
    Ca = clamp(VS / alvo.PP, CA_MIN, CA_MAX)
    Va0 = 0.165 * Vtot
    cv = com(cv_ref; FC = p.fc, Emax = cv_ref.Emax / s, Emin = cv_ref.Emin / s, Vd = cv_ref.Vd * s,
             Vtotal = Vtot, Va0 = Va0, Rs = Rs, Ca = Ca)
    Vv0 = Vtot - (Va0 + Ca * alvo.PAM) - (cv.Vd + 8.0 / cv.Emin) - cv.Cv * 8.0
    cv = com(cv; Vv0 = Vv0)
    u = estado_inicial(cv; P0 = 8.0)
    u[I_VA] = Va0 + Ca * alvo.PAM
    u[I_VV] = Vtot - u[I_VA] - u[I_VLV]
    baro = com(baro_ref; Pset = alvo.PAM)
    med, u = _regime(cv, baro, pk, pd, u, 40)
    it = 0
    convergiu = false
    for k in 1:maxit
        it = k
        eP = med.PAM / alvo.PAM - 1
        ePP = (med.PAS - med.PAD) / alvo.PP - 1
        eQ = med.DC / alvo.DC - 1
        verbose && @printf("  it %2d: PAM %.1f  PP %.1f  DC %.2f | Rs %.3f Ca %.3f V0v %.0f\n",
                           k, med.PAM, med.PAS - med.PAD, med.DC, cv.Rs, cv.Ca, cv.Vv0)
        if max(abs(eP), abs(ePP), abs(eQ)) < tol
            convergiu = true
            break
        end
        # --- correções de ponto fixo ---
        Pv = (u[I_VV] - cv.Vv0) / cv.Cv
        # Rs ≈ (PAM−Pv)/DC: corrige pela razão entre o desejado e o obtido
        Rs_novo = clamp(cv.Rs * ((alvo.PAM - Pv) / (med.PAM - Pv)) / (alvo.DC / med.DC), 0.2, 5.0)
        razao_pp = (med.PAS - med.PAD) / alvo.PP
        Ca_novo = clamp(cv.Ca * razao_pp^1.2, CA_MIN, CA_MAX)
        # Se a complacência saturou no limite fisiológico, a pressão de pulso
        # restante vem da impedância característica (Zc·Q): ajusta-se Zc.
        Zc_novo = cv.Zc
        if (Ca_novo >= CA_MAX && razao_pp > 1) || (Ca_novo <= CA_MIN && razao_pp < 1)
            Zc_novo = clamp(cv.Zc / razao_pp^1.5, ZC_MIN, ZC_MAX)
        end
        # DC ∝ volume estressado (via Pms): Vest_novo ≈ Vest·DC*/DC  ⇒ V0v_novo = V0v − ΔVest
        Vest = cv.Vtotal - cv.Vv0 - cv.Va0 - cv.Vd    # volume estressado total
        Vv0_novo = clamp(cv.Vv0 + 0.7 * eQ * Vest / (1 + eQ), 0.2 * cv.Vtotal, cv.Vtotal - cv.Va0 - 150)
        # ao mudar Ca ou Va0 mantemos a pressão arterial, realocando volume
        Pc = (u[I_VA] - cv.Va0) / cv.Ca
        dVa = (Ca_novo - cv.Ca) * Pc
        u[I_VA] += dVa
        u[I_VV] -= dVa
        cv = com(cv; Rs = Rs_novo, Ca = Ca_novo, Zc = Zc_novo, Vv0 = Vv0_novo)
        med, u = _regime(cv, baro, pk, pd, u, 15)
    end
    med, u = _regime(cv, baro, pk, pd, u, 10)
    baro = com(baro; Pset = med.PAM)
    u0 = copy(u)
    u0[I_FASE] = 0.0; u0[I_INTP] = 0.0; u0[I_INTQ] = 0.0
    u0[I_A1] = u0[I_A2] = u0[I_A3] = u0[I_CEB] = u0[I_CEH] = 0.0
    u0[I_PS] = med.PAM
    u0[I_XFC] = u0[I_XR] = u0[I_XE] = u0[I_XV] = 0.0
    PPb = amp * (med.PAS - med.PAD)
    PADb = med.PAM - obs.ff * PPb
    obtido = (PAM = med.PAM, PAS = PADb + PPb, PAD = PADb, PASao = med.PAS, PADao = med.PAD,
              DC = med.DC, VS = med.VS, FC = cv.FC)
    return PacienteVirtual(p, pk, pd, cv, baro, obs, u0, alvo, obtido, it, convergiu)
end

function Base.show(io::IO, pv::PacienteVirtual)
    cv, a, o = pv.cv, pv.alvo, pv.obtido
    println(io, "Paciente virtual — ", pv.paciente)
    @printf(io, "  calibração (%s em %d iterações):\n", pv.convergiu ? "convergiu" : "NÃO convergiu", pv.iteracoes)
    @printf(io, "    PAS %6.1f / alvo %6.1f mmHg   (aórtica %.1f)\n", o.PAS, a.PAS, o.PASao)
    @printf(io, "    PAD %6.1f / alvo %6.1f mmHg   (aórtica %.1f)\n", o.PAD, a.PAD, o.PADao)
    @printf(io, "    PAM %6.1f / alvo %6.1f mmHg\n", o.PAM, a.PAM)
    @printf(io, "    amplificação da PP aorta→braço assumida: %.2f\n", pv.obs.amp)
    @printf(io, "    DC  %6.2f / alvo %6.2f L/min   (VS %.0f ml)\n", o.DC, a.DC, o.VS)
    @printf(io, "  parâmetros identificados:\n")
    @printf(io, "    TPR = %.3f mmHg·s/ml   AC = %.2f ml/mmHg   τ = TPR·AC = %.2f s   Zc = %.3f mmHg·s/ml\n",
            cv.Rs, cv.Ca, cv.Rs * cv.Ca, cv.Zc)
    @printf(io, "    C_S = 1/Emax = %.2f ml/mmHg   C_D = 1/Emin = %.1f ml/mmHg\n", 1 / cv.Emax, 1 / cv.Emin)
    Pms = pressao_media_enchimento(cv)
    @printf(io, "    volume circulante %.0f ml, não estressado venoso %.0f ml, Pms = %.1f mmHg\n", cv.Vtotal, cv.Vv0, Pms)
    @printf(io, "  PK (Schnider): V1 %.2f L, V2 %.1f L, V3 %.0f L, Cl1 %.2f L/min, t½ke0 hemodinâmico %.1f min",
            pv.pk.V1, pv.pk.V2, pv.pk.V3, pv.pk.Cl1, log(2) / pv.pk.ke0_hemo / 60)
end
