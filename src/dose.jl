# =============================================================================
# dose.jl — exploração de doses "antes da injeção"
# =============================================================================

"""
    varrer_doses(pv, doses; janela=20.0, kw...) -> Vector de NamedTuple

Simula cada dose (mg/kg) e devolve as métricas. `kw` é repassado a
`simular_inducao` (ex.: `duracao_bolus`, `manutencao_ugkgmin`).
"""
varrer_doses(pv, doses; janela = 20.0, kw...) =
    [merge((dose = d,), metricas(simular_inducao(pv; dose_mgkg = d, duracao = janela, kw...); janela)) for d in doses]

"""
    bissecao(f, a, b; tol=0.01)

Menor `x` em [a, b] com f(x) == true, supondo f monótona (false…false true…true).
"""
function bissecao(f, a, b; tol = 0.01)
    f(a) && return a
    f(b) || return NaN
    while b - a > tol
        c = (a + b) / 2
        f(c) ? (b = c) : (a = c)
    end
    return b
end

"""
    janela_de_dose(pv; limiar_pam=65.0, alvo_bis=60.0, faixa=(0.5, 3.5), kw...)

Janela terapêutica individual (mg/kg): `dose_min` é a menor dose que leva o BIS
abaixo de `alvo_bis`; `dose_max` é a maior dose cujo nadir de PAM fica ≥ `limiar_pam`.
"""
function janela_de_dose(pv; limiar_pam = 65.0, alvo_bis = 60.0, faixa = (0.5, 3.5), janela = 20.0, kw...)
    sim(d) = metricas(simular_inducao(pv; dose_mgkg = d, duracao = janela, kw...); janela, alvo_bis)
    dmin = bissecao(d -> sim(d).bis_min < alvo_bis, faixa...)
    dmax_viola = bissecao(d -> sim(d).pam_nadir < limiar_pam, faixa...)
    dmax = isnan(dmax_viola) ? faixa[2] : dmax_viola
    return (dose_min = dmin, dose_max = dmax, existe = !isnan(dmin) && dmin <= dmax,
            limiar_pam = limiar_pam, alvo_bis = alvo_bis)
end
