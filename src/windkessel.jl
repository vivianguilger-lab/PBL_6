# Identificação de modelos Windkessel (2, 3 e 4 elementos).
# Porta para Julia moderno de oeid.jl do windkessel-id (Pigot et al., licença MIT,
# github.com/hpigot/windkessel-id), feita pelo grupo em src/windkessel.jl do
# PropofolPBL e reproduzida aqui sem alterações nas equações.
#
# Diferenças em relação ao original:
# - Não depende de ControlSystems.jl. Os três modelos têm matriz A diagonal,
#   então a discretização ZOH é feita de forma analítica (exp elemento a
#   elemento), o que funciona direto com ForwardDiff e dispensa expm.jl.
# - `simulate` virou `simulate_wk` para não colidir com `NeoPKPD.simulate`.
#
# Unidades (como em Stergiopulos 1999): pressão em mmHg, fluxo em mL/s,
# Rp e Rc em mmHg·s/mL, C em mL/mmHg, L em mmHg·s²/mL, tempo em s.

"""
Sistema linear em espaço de estados com matriz A diagonal:
    ẋ = diag(a)·x + b·u,   y = c'x + d·u
Se `h > 0`, é a versão discreta (ZOH) com período `h`:
    x[k+1] = diag(a)·x[k] + b·u[k],   y[k] = c'x[k] + d·u[k]
"""
struct WKSistema{T<:Real}
    a::Vector{T}
    b::Vector{T}
    c::Vector{T}
    d::T
    h::Float64      # 0.0 para tempo contínuo
end

function _wk(a, b, c, d, h, discretize)
    T = promote_type(eltype(a), eltype(b), eltype(c), typeof(d))
    a, b, c, d = T.(a), T.(b), T.(c), T(d)
    discretize || return WKSistema(a, b, c, d, 0.0)
    ad = exp.(a .* h)
    bd = (ad .- 1) ./ a .* b
    return WKSistema(ad, bd, c, d, Float64(h))
end

"""
Windkessel paralelo de 4 elementos, θ = [Rp, C, Rc, L], discretizado (ZOH) com período h.
x1 é o volume no capacitor C (forma paralela de Burattini, Westerhof, Stergiopulos).
"""
function G(θ, h; discretize=true)
    Rp, C, Rc, L = θ
    return _wk([-1 / (C * Rp), -Rc / L], [1, Rc], [1 / C, -Rc / L], Rc, h, discretize)
end

"""
Windkessel de 3 elementos, θ = [Rp, C, Rc]. Corresponde ao de 4 elementos com L = Inf.
"""
function G3(θ, h; discretize=true)
    Rp, C, Rc = θ
    return _wk([-1 / (C * Rp)], [1], [1 / C], Rc, h, discretize)
end

"""
Windkessel de 2 elementos, θ = [Rp, C]. Corresponde ao de 4 elementos com L = 0.
"""
function G2(θ, h; discretize=true)
    Rp, C = θ
    return _wk([-1 / (C * Rp)], [1], [1 / C], 0, h, discretize)
end

"""
Simula o sistema discreto `sys` com entrada `u` a partir do estado `x0`.
Retorna `(y, x)`, com `x[k, :]` o estado na amostra k (mesma convenção do `lsim` do ControlSystems).
"""
function lsim_wk(sys::WKSistema, u::AbstractVector, x0=zeros(length(sys.a)))
    sys.h > 0 || throw(ArgumentError("lsim_wk requer um sistema discreto"))
    T = promote_type(eltype(sys.a), eltype(u), eltype(x0))
    n, nx = length(u), length(sys.a)
    y = Vector{T}(undef, n)
    x = Matrix{T}(undef, n, nx)
    xk = T.(x0)
    for k in 1:n
        x[k, :] = xk
        y[k] = dot(sys.c, xk) + sys.d * u[k]
        xk = sys.a .* xk .+ sys.b .* u[k]
    end
    return y, x
end

"Estado após simular `u` inteiro a partir de `x0` (útil para encadear batimentos)."
function estado_final(sys::WKSistema, u::AbstractVector, x0)
    xk = copy(x0)
    for uk in u
        xk = sys.a .* xk .+ sys.b .* uk
    end
    return xk
end

"""
Simula o modelo `G_(θ, h)` com fluxo `u`, a partir de `x0`. Retorna a pressão
(e os estados, se `xout=true`).
"""
function simulate_wk(θ, x0, u, h; G_=G, xout=false)
    y, x = lsim_wk(G_(θ, h), u, x0)
    return xout ? (y, x) : y
end

"""
Resolve x0 explorando a periodicidade da entrada: se `u` é um período, o estado
inicial periódico satisfaz x0 = Aⁿx0 + xₙ⁰, onde xₙ⁰ é o estado após n passos
partindo do repouso. Com A diagonal: x0ᵢ = xₙ⁰ᵢ / (1 - aᵢⁿ).
"""
function solvex0(sys::WKSistema, u::AbstractVector)
    xn0 = estado_final(sys, u, zeros(eltype(sys.a), length(sys.a)))
    return xn0 ./ (1 .- sys.a .^ length(u))
end

"""
Resposta em frequência do modelo contínuo em `w` (rad/s): retorna `(magnitude, fase_graus)`.
Substitui o `bode` do ControlSystems usado nos scripts originais.
"""
function bode_wk(θ, w; G_=G)
    sys = G_(θ, 0.0; discretize=false)
    H = [sum(sys.c .* sys.b ./ (im * wi .- sys.a)) + sys.d for wi in w]
    return abs.(H), rad2deg.(angle.(H))
end

# ------------------------------------------------------------------
# Identificação (erro de saída, Newton com diferenciação automática)
# ------------------------------------------------------------------

struct IdConditions
    θ0::Vector{Float64}     # e.g. [Rp, C, Rc, L]
    x00::Vector{Float64}
    θid::Vector{Int}        # e.g. [1,1,1,0]
    x0id::Bool              # identificar x0
    x0solve::Bool           # resolver x0 pela periodicidade a cada passo
    G_::Function            # modelo, na forma `G_(θ, h; discretize=true)`
    function IdConditions(θ0, x00, θid, x0id, x0solve, G_=G)
        length(θ0) == length(θid) || error("Length of `θ0` != `θid`.")
        new(θ0, x00, θid, x0id, x0solve, G_)
    end
end

"Monta o vetor de otimização."
function p(θ, x0, id::IdConditions)
    pv = θ[findall(Bool.(id.θid))]
    if id.x0id
        pv = vcat(pv, x0)
    end
    return Float64.(pv)
end

"Recupera a parametrização original (θ, x0) a partir do vetor de otimização."
function deal(pv, id::IdConditions)
    pidx = cumsum(id.θid)
    θ = [Bool(id_) ? pv[i] : θ0_ for (i, θ0_, id_) in zip(pidx, id.θ0, id.θid)]
    x0 = id.x0id ? pv[sum(id.θid).+(1:length(id.x00))] : id.x00
    return θ, x0
end

"Custo (soma dos quadrados do erro de saída) a ser minimizado."
function cost(pv, id::IdConditions, u, y, h)
    θ, x0 = deal(pv, id)
    sys = id.G_(θ, h)
    if id.x0solve
        x0 = solvex0(sys, u)
    end
    ym, = lsim_wk(sys, u, x0)
    e = y .- ym
    return dot(e, e)
end

"""
Identifica um subconjunto de {parâmetros θ, estado inicial x0} pelo método de Newton
(IPNewton com θ ≥ 0 se `θpos=true`). Retorna `(phat, res)`.
"""
function identify(θ0, x00, u, y, h; θid=ones(Int, length(θ0)), x0id=true, x0solve=false,
                  G_=G, θpos=true, v=false)
    id = IdConditions(θ0, x00, θid, x0id, x0solve, G_)
    p0 = p(θ0, x00, id)
    f(pv) = cost(pv, id, u, y, h)
    g!(grad, pv) = ForwardDiff.gradient!(grad, f, pv)
    function h!(hess, pv)
        ForwardDiff.hessian!(hess, f, pv)
        hess .= (hess .+ hess') ./ 2
    end

    df = TwiceDifferentiable(f, g!, h!, p0)
    opts = Optim.Options(store_trace=true, extended_trace=true, show_trace=v)
    res = if θpos
        lp = vcat(zeros(sum(id.θid)), fill(-Inf, id.x0id ? length(x00) : 0))
        up = fill(Inf, length(p0))
        dfc = TwiceDifferentiableConstraints(lp, up)
        optimize(df, dfc, p0, IPNewton(), opts)
    else
        optimize(df, p0, Newton(alphaguess=LineSearches.InitialPrevious(),
                                linesearch=LineSearches.BackTracking()), opts)
    end
    return Optim.minimizer(res), res
end

"""
Identifica com Nelder-Mead (busca global, sem restrição de positividade).
"""
function identifyglobal(θ0, x00, u, y, h; θid=ones(Int, length(θ0)), x0id=true,
                        x0solve=false, G_=G, v=false)
    id = IdConditions(θ0, x00, θid, x0id, x0solve, G_)
    f(pv) = cost(pv, id, u, y, h)
    opts = Optim.Options(iterations=10^5, store_trace=true, extended_trace=true, show_trace=v)
    res = optimize(f, p(θ0, x00, id), NelderMead(), opts)
    return Optim.minimizer(res), res
end

"""
Parâmetros de Stergiopulos et al. (1999), batimento tipo A, humano.
"""
function params_stergiopulos1999total()
    Rp = 0.79     # resistência periférica [mmHg·s/mL]
    C = 1.22      # complacência arterial total [mL/mmHg]
    Rc = 0.056    # resistência característica da aorta [mmHg·s/mL]
    L = 0.0051    # inertância arterial total [mmHg·s²/mL]
    return [Rp, C, Rc, L], [8.066, -8.296]
end
