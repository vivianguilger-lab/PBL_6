# =============================================================================
# calibracao.jl — paciente virtual calibrado com os dados pré-indução
# =============================================================================
#
# Base: theta_basal do grupo — Windkessel de 4 elementos com Rp = PAM/DC e
# C, Rc, L de Stergiopulos (1999). Como o monitor fornece também PAS e PAD,
# identificamos a complacência arterial C do próprio paciente para que a
# pressão de pulso do modelo seja a medida (mesma ideia de identificação do
# windkessel-id, aqui com um único parâmetro e uma única medida).
#
#   PAM = Rp · Q̄            (Z(0) = Rp)        → Rp
#   PP  = PAS − PAD  ↓  quando C ↑              → C (bisseção em log C)
# =============================================================================

const C_MIN, C_MAX = 0.2, 6.0   # ml/mmHg

"""
Paciente virtual: biomarcadores + Windkessel basal identificado + PK + efeito.
"""
struct PacienteVirtual
    paciente::Paciente
    θ::Vector{Float64}             # [Rp, C, Rc, L] (mmHg, ml, s)
    pk::NamedTuple
    efeito::EfeitoHemodinamico
    s::Float64                     # sensibilidade hemodinâmica
    alvo::NamedTuple
    obtido::NamedTuple
    convergiu::Bool
end

"Constante de tempo diastólica τ = Rp·C [s]."
tau_diastolico(pv::PacienteVirtual) = pv.θ[1] * pv.θ[2]

function _pp(θ, fc, vs)
    y, _, _ = batimento_basal(θ, fc, vs)
    return maximum(y) - minimum(y), y
end

"""
    calibrar(p; modelo_pk=:schnider, efeito=EfeitoHemodinamico(), eta_pk=NamedTuple(),
             s=sensibilidade_hemodinamica(p), ajustar_C=true) -> PacienteVirtual
"""
function calibrar(p::Paciente; modelo_pk = :schnider, efeito = EfeitoHemodinamico(), eta_pk = NamedTuple(),
                  s = sensibilidade_hemodinamica(p), ajustar_C = true)
    θ = theta_basal(p)
    vs = debito_basal(p) * 1000 / p.fc
    alvo = (PAS = p.pas, PAD = p.pad, PAM = pam_basal(p), PP = p.pas - p.pad, DC = debito_basal(p), FC = p.fc)
    convergiu = true
    if ajustar_C
        ppC(C) = _pp([θ[1], C, θ[3], θ[4]], p.fc, vs)[1]
        if ppC(C_MAX) > alvo.PP          # pulso pequeno demais até para a maior complacência
            θ[2] = C_MAX; convergiu = false
        elseif ppC(C_MIN) < alvo.PP      # pulso grande demais até para a menor complacência
            θ[2] = C_MIN; convergiu = false
        else
            a, b = log(C_MIN), log(C_MAX)
            for _ in 1:60
                c = (a + b) / 2
                ppC(exp(c)) > alvo.PP ? (a = c) : (b = c)
            end
            θ[2] = exp((a + b) / 2)
        end
    end
    _, y = _pp(θ, p.fc, vs)
    obtido = (PAS = maximum(y), PAD = minimum(y), PAM = mean(y), PP = maximum(y) - minimum(y),
              DC = debito_basal(p), VS = vs, FC = p.fc)
    return PacienteVirtual(p, θ, parametros_pk(p, modelo_pk; eta = eta_pk), efeito, s, alvo, obtido, convergiu)
end

function Base.show(io::IO, pv::PacienteVirtual)
    Rp, C, Rc, L = pv.θ
    a, o = pv.alvo, pv.obtido
    println(io, "Paciente virtual — ", pv.paciente)
    @printf(io, "  PA do modelo %.0f/%.0f (PAM %.1f) mmHg | medida %.0f/%.0f (PAM %.1f)%s\n",
            o.PAS, o.PAD, o.PAM, a.PAS, a.PAD, a.PAM, pv.convergiu ? "" : "  ⚠ pressão de pulso não atingida")
    @printf(io, "  Windkessel 4 elementos: Rp = %.3f mmHg·s/ml, C = %.2f ml/mmHg, Rc = %.3f, L = %.4f | τ = Rp·C = %.2f s\n",
            Rp, C, Rc, L, Rp * C)
    @printf(io, "  DC %.2f L/min, VS %.0f ml | sensibilidade hemodinâmica s = %.2f\n", o.DC, o.VS, pv.s)
    pk = pv.pk
    @printf(io, "  PK (%s): CL %.2f L/min, V1 %.2f L, V2 %.1f L, V3 %.0f L, ke0 %.3f /min",
            pk.modelo, pk.CL, pk.V1, pk.V2, pk.V3, pk.ke0)
end

# =============================================================================
# Simulação da indução a partir do paciente virtual
# =============================================================================

"""
Resultado de uma indução simulada: PK/PD na grade de 1 s (`pkpd`), valores
por batimento (`bat`) e, opcionalmente, a forma de onda (`onda`). Tempos em min.
"""
struct Resultado
    pkpd::NamedTuple
    bat::NamedTuple
    onda::Union{Nothing,NamedTuple}
    t_inducao::Float64
    pv::PacienteVirtual
    protocolo::NamedTuple
end

"""
    simular_inducao(pv; dose_mgkg=2.0, duracao_bolus=30.0, manutencao_ugkgmin=100.0,
                    inicio_manutencao=1.0, t_inducao=1.0, duracao=20.0, forma_de_onda=false,
                    pk=pv.pk, efeito=pv.efeito, s=pv.s) -> Resultado

Bolus de `dose_mgkg` injetado em `duracao_bolus` **segundos** no minuto
`t_inducao` (o minuto anterior é a linha de base), seguido de manutenção a partir
de `inicio_manutencao` min após o bolus (esquema do grupo). Simula `duracao` min.
"""
function simular_inducao(pv::PacienteVirtual; dose_mgkg = 2.0, duracao_bolus = 30.0, manutencao_ugkgmin = 100.0,
                         inicio_manutencao = 1.0, t_inducao = 1.0, duracao = 20.0, forma_de_onda = false,
                         pk = pv.pk, efeito = pv.efeito, s = pv.s)
    p = pv.paciente
    tfim = t_inducao + duracao
    doses = esquema_inducao(p; bolus_mgkg = dose_mgkg, duracao_bolus = duracao_bolus / 60,
                            manutencao_mcgkgmin = manutencao_ugkgmin, inicio_manutencao,
                            duracao_manutencao = max(duracao - inicio_manutencao, 0.0), t_bolus = t_inducao)
    pkpd = simular_propofol(p, doses; tfim, pk)
    bat, onda = simular_hemodinamica(p, pkpd; efeito, θbasal = pv.θ, s, tfim, forma_de_onda)
    return Resultado(pkpd, bat, onda, t_inducao, pv,
                     (dose_mgkg = dose_mgkg, duracao_bolus = duracao_bolus, manutencao_ugkgmin = manutencao_ugkgmin))
end
