# =============================================================================
# acoplamento.jl — PK/PD (NeoPKPD) → hemodinâmica (Windkessel do windkessel-id)
# =============================================================================
#
# Segue src/acoplamento.jl do PropofolPBL do grupo:
#
#   Ce(t) ──Emax sigmoide──► Rp, volume sistólico (VS) e FC
#                               │
#   fluxo aórtico (Stergiopulos 1999, escalado por VS e FC)
#                               ▼
#   Windkessel de 4 elementos θ = [Rp, C, Rc, L]  ──►  pressão arterial
#
# A cada batimento o Windkessel é rediscretizado com os parâmetros daquele
# instante e simulado continuando o estado do batimento anterior.
#
# ATENÇÃO (como no README do grupo): os valores padrão de `EfeitoHemodinamico`
# e de `sensibilidade_hemodinamica` são ilustrativos, não ajustados a dados.
# =============================================================================

"""
Efeito do propofol sobre a hemodinâmica. Cada variável X é multiplicada por
`1 − Imax_X · s · Ceᵞ/(Ce50ᵞ + Ceᵞ)` (mínimo 0,05), com `s` a sensibilidade
do paciente.
"""
Base.@kwdef struct EfeitoHemodinamico
    Ce50::Float64 = 3.0       # µg/ml
    γ::Float64 = 2.0
    Imax_Rp::Float64 = 0.25   # vasodilatação arterial
    Imax_VS::Float64 = 0.15   # depressão miocárdica / redução de pré-carga
    Imax_FC::Float64 = 0.05
end

"Fator multiplicativo 1 − Imax·s·H(Ce), limitado a 0,05."
fator_efeito(ef::EfeitoHemodinamico, Imax, s, Ce) =
    max(1 - Imax * s * Ce^ef.γ / (ef.Ce50^ef.γ + Ce^ef.γ), 0.05)

"""
    sensibilidade_hemodinamica(pac; k_idade=0.01, k_asa=0.15, k_pvi=0.03, k_ivci=0.01)

Fator de sensibilidade hemodinâmica (1 = paciente de referência). Aumenta com
idade > 50 anos, ASA ≥ 3, PVI > 15% e IVCI > 43%, fatores associados à hipotensão
pós-indução. Limiares: PVI de Tsuchiya et al. 2010 (PVI > 15 prevê queda da PAM após
propofol); IVCI de Zhang & Critchley 2016 (corte ótimo de 43%). Os coeficientes `k_*`
são do grupo, ilustrativos, ainda sem calibração (docs/validacao.md).
"""
function sensibilidade_hemodinamica(pac::Paciente; k_idade = 0.01, k_asa = 0.15, k_pvi = 0.03, k_ivci = 0.01)
    return 1.0 + k_idade * max(pac.idade - 50, 0) + k_asa * max(pac.asa - 2, 0) +
           k_pvi * max(pac.pvi - 15, 0) + k_ivci * max(pac.ivci - 43, 0)
end

"""
    theta_basal(pac)

Parâmetros basais θ = [Rp, C, Rc, L] (unidades mL/s) do Windkessel de 4
elementos: Rp = PAM/DC (em regime, Z(0) = Rp); C, Rc e L de Stergiopulos (1999).
"""
function theta_basal(pac::Paciente)
    θ, _ = params_stergiopulos1999total()
    θ[1] = pam_basal(pac) / (debito_basal(pac) * 1000 / 60)
    return θ
end

const _FORMA_FLUXO = Ref{Vector{Float64}}()

"""
Forma de onda normalizada do fluxo aórtico (média 1 em τ ∈ [0, 1)), a partir
do batimento de Stergiopulos (1999) distribuído com o windkessel-id.
"""
function forma_fluxo()
    if !isassigned(_FORMA_FLUXO)
        _, _, q = dados_stergiopulos1999()
        _FORMA_FLUXO[] = q ./ (sum(q) / length(q))
    end
    return _FORMA_FLUXO[]
end

"Reamostra `forma` (um período) para `n` pontos e escala para ejetar `vs` ml em `n·h` s."
function fluxo_batimento(forma, n, vs, h)
    m = length(forma)
    q = map(0:n-1) do k
        s = k * m / n
        i = floor(Int, s)
        f = s - i
        (1 - f) * forma[i+1] + f * forma[mod(i + 1, m)+1]
    end
    return q .* (vs / (sum(q) * h))
end

"Interpolação linear de (tx, yx) em t, constante fora do intervalo."
function interp(tx, yx, t)
    t <= tx[1] && return yx[1]
    t >= tx[end] && return yx[end]
    i = searchsortedlast(tx, t)
    f = (t - tx[i]) / (tx[i+1] - tx[i])
    return (1 - f) * yx[i] + f * yx[i+1]
end

"""
    batimento_basal(θ, fc, vs; h=0.005)

Batimento periódico (regime permanente) do Windkessel θ com FC e VS dados.
Devolve `(pressao, fluxo, x0)`.
"""
function batimento_basal(θ, fc, vs; h = 0.005, forma = forma_fluxo())
    n = round(Int, 60 / fc / h)
    u = fluxo_batimento(forma, n, vs, h)
    sys = G(θ, h)
    x0 = solvex0(sys, u)
    y, = lsim_wk(sys, u, x0)
    return y, u, x0
end

"""
    simular_hemodinamica(pac, pkpd; efeito, θbasal=theta_basal(pac), s=sensibilidade_hemodinamica(pac),
                         h=0.005, tfim=pkpd.t[end], forma_de_onda=false)

Pressão arterial batimento a batimento durante a PK/PD `pkpd` (saída de
`simular_propofol`). Devolve `(batimentos, onda)`: `batimentos` tem, por
batimento, `t_min, Ce, BIS, FC, VS, DC, Rp, PAS, PAD, PAM`; `onda` é
`(t_s, pressao, fluxo)` se `forma_de_onda = true`.
"""
function simular_hemodinamica(pac::Paciente, pkpd; efeito = EfeitoHemodinamico(), θbasal = theta_basal(pac),
                              s = sensibilidade_hemodinamica(pac), h = 0.005, tfim = pkpd.t[end],
                              forma = forma_fluxo(), forma_de_onda = false)
    Rp0, C, Rc, L = θbasal
    vs0 = debito_basal(pac) * 1000 / pac.fc     # ml por batimento
    _, _, x = batimento_basal(θbasal, pac.fc, vs0; h, forma)

    col = Dict(k => Float64[] for k in (:t_min, :Ce, :BIS, :FC, :VS, :DC, :Rp, :PAS, :PAD, :PAM))
    ts, ps, qs = Float64[], Float64[], Float64[]
    t = 0.0                                      # s
    while t < tfim * 60
        Ce = interp(pkpd.t, pkpd.Ce, t / 60)
        fc = pac.fc * fator_efeito(efeito, efeito.Imax_FC, s, Ce)
        vs = vs0 * fator_efeito(efeito, efeito.Imax_VS, s, Ce)
        Rp = Rp0 * fator_efeito(efeito, efeito.Imax_Rp, s, Ce)
        n = round(Int, 60 / fc / h)
        u = fluxo_batimento(forma, n, vs, h)
        sys = G([Rp, C, Rc, L], h)
        y, = lsim_wk(sys, u, x)
        x = estado_final(sys, u, x)
        for (k, v) in zip((:t_min, :Ce, :BIS, :FC, :VS, :DC, :Rp, :PAS, :PAD, :PAM),
                          (t / 60, Ce, interp(pkpd.t, pkpd.bis, t / 60), fc, vs, vs * fc / 1000, Rp,
                           maximum(y), minimum(y), mean(y)))
            push!(col[k], v)
        end
        if forma_de_onda
            append!(ts, t .+ h .* (0:n-1)); append!(ps, y); append!(qs, u)
        end
        t += n * h
    end
    bat = NamedTuple{(:t_min, :Ce, :BIS, :FC, :VS, :DC, :Rp, :PAS, :PAD, :PAM)}(
        Tuple(col[k] for k in (:t_min, :Ce, :BIS, :FC, :VS, :DC, :Rp, :PAS, :PAD, :PAM)))
    return bat, forma_de_onda ? (t_s = ts, pressao = ps, fluxo = qs) : nothing
end

"""
    simular_paciente(pac; doses=esquema_inducao(pac), modelo=:schnider, tfim=30.0, kwargs...)

Atalho do grupo: PK/PD (NeoPKPD) + hemodinâmica (Windkessel). Devolve `(pkpd, batimentos, onda)`.
"""
function simular_paciente(pac::Paciente; doses = esquema_inducao(pac), modelo = :schnider, tfim = 30.0, kwargs...)
    pkpd = simular_propofol(pac, doses; modelo, tfim)
    batimentos, onda = simular_hemodinamica(pac, pkpd; tfim, kwargs...)
    return pkpd, batimentos, onda
end
