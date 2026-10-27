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

# -----------------------------------------------------------------------------
# Visão geral (para quem não programa)
# -----------------------------------------------------------------------------
# O Windkessel ("câmara de ar", em alemão) representa a árvore arterial como um
# circuito elétrico. Analogia: pressão = tensão (volts), fluxo de sangue =
# corrente (ampères), volume de sangue = carga elétrica.
#   - Rp: resistência periférica (arteríolas). É ela que define a pressão média.
#   - C: complacência arterial (elasticidade da aorta e das grandes artérias).
#     Age como um capacitor: guarda sangue na sístole e o devolve na diástole.
#   - Rc: impedância característica da aorta proximal. Pequena resistência que o
#     coração "sente" no início da ejeção; pesa sobretudo na pressão sistólica.
#   - L: inertância (massa de sangue que precisa ser acelerada). Age como um indutor.
# Entrada do modelo: fluxo aórtico Q(t), em mL/s. Saída: pressão aórtica P(t), em mmHg.
# A pressão é medida em relação a zero (pressão venosa considerada desprezível).
#
# Por que PAM = Rp · fluxo médio: no regime médio (frequência zero) o capacitor C
# não deixa passar fluxo médio e o indutor L "curto-circuita" Rc. Sobra só Rp, como
# na lei de Ohm (V = R·I). Por isso, no paciente, Rp = PAM/DC (acoplamento.jl).
#
# Como este arquivo entra no modelo: acoplamento.jl monta θ = [Rp, C, Rc, L] a cada
# batimento (o propofol reduz Rp) e usa G, solvex0, lsim_wk e estado_final daqui
# para calcular a curva de pressão. A parte de identificação (identify) só é usada
# nos scripts e testes, para ajustar θ a medidas de pressão e fluxo
# (scripts/06_identificacao_windkessel.jl).
#
# De onde vem cada valor (Stergiopulos 1999, windkessel-id, PropofolPBL do grupo):
# ver docs/validacao.md, seção 4.3.

# WKSistema: guarda os coeficientes de um sistema linear (equações na docstring
# logo abaixo). x é o vetor de estados, u a entrada (fluxo, mL/s) e y a saída
# (pressão, mmHg). Como a matriz A é diagonal, basta guardar a sua diagonal (vetor a).
# Sintaxe Julia: `struct` cria um tipo de dado com campos nomeados; `{T<:Real}` diz
# que os números podem ser de qualquer tipo real (Float64 ou os "números duais" que
# o ForwardDiff usa para calcular derivadas durante a identificação).
"""
Sistema linear em espaço de estados com matriz A diagonal:
    ẋ = diag(a)·x + b·u,   y = c'x + d·u
Se `h > 0`, é a versão discreta (ZOH) com período `h`:
    x[k+1] = diag(a)·x[k] + b·u[k],   y[k] = c'x[k] + d·u[k]
"""
struct WKSistema{T<:Real}
    # a: diagonal de A (contínuo, em 1/s) ou de exp(A·h) (discreto, sem unidade).
    # b: quanto a entrada (fluxo) altera cada estado.
    # c, d: como os estados e o fluxo formam a pressão de saída (y = c'x + d·u).
    # h: passo de tempo em s (5 ms no modelo); 0.0 quer dizer tempo contínuo.
    a::Vector{T}
    b::Vector{T}
    c::Vector{T}
    d::T
    h::Float64      # 0.0 para tempo contínuo
end

# _wk: monta um WKSistema e, se pedido, o discretiza por segurador de ordem zero
# (ZOH, "zero-order hold"). O ZOH supõe que o fluxo fica constante dentro de cada
# passo h (5 ms) e calcula o estado no fim do passo de forma exata.
# Para cada estado i (A diagonal), dx/dt = a·x + b·u com u constante dá:
#   x[k+1] = e^(a·h)·x[k] + (e^(a·h) − 1)/a · b · u[k]
# Entradas: a, b, c, d (coeficientes), h (s) e discretize (true ou false).
# O "_" no início do nome indica função interna (não exportada).
function _wk(a, b, c, d, h, discretize)
    # Converte todos os coeficientes para um mesmo tipo numérico T (necessário para
    # o ForwardDiff calcular derivadas). `T.(a)` aplica a conversão a cada elemento:
    # o ponto "." é o broadcasting de Julia (aplicar elemento a elemento).
    T = promote_type(eltype(a), eltype(b), eltype(c), typeof(d))
    a, b, c, d = T.(a), T.(b), T.(c), T(d)
    # `discretize || return ...`: se discretize for falso, devolve o sistema
    # contínuo e termina aqui (atalho de Julia para "se não ..., então ...").
    discretize || return WKSistema(a, b, c, d, 0.0)
    # Fórmula ZOH exata, elemento a elemento: ad = e^(a·h) e bd = (ad − 1)/a · b.
    ad = exp.(a .* h)
    bd = (ad .- 1) ./ a .* b
    return WKSistema(ad, bd, c, d, Float64(h))
end

# G: Windkessel de 4 elementos (o usado no modelo do paciente).
# Entrada: θ = [Rp, C, Rc, L] em unidades de mL/s (ver cabeçalho) e passo h em s.
# Saída: WKSistema discreto. Estados:
#   x1 = volume guardado na complacência C (mL); a pressão em C é x1/C;
#   x2 = L × fluxo no ramo da inertância (mmHg·s); x2/L é esse fluxo.
# Equações (docs/modelo_matematico.md, seção 4):
#   dx1/dt = −x1/(C·Rp) + Q       (C enche com Q e esvazia através de Rp)
#   dx2/dt = −(Rc/L)·x2 + Rc·Q    (Rc e L em paralelo dividem o fluxo)
#   P = x1/C − (Rc/L)·x2 + Rc·Q
# Frequência zero: |Z(0)| = Rp (PAM = Rp·Q̄); frequência alta: |Z(∞)| = Rc.
# `discretize=true` é um argumento nomeado (keyword) com valor padrão.
"""
Windkessel paralelo de 4 elementos, θ = [Rp, C, Rc, L], discretizado (ZOH) com período h.
x1 é o volume no capacitor C (forma paralela de Burattini, Westerhof, Stergiopulos).
"""
function G(θ, h; discretize=true)
    # Desempacota θ: a linha abaixo atribui os 4 valores do vetor de uma vez.
    Rp, C, Rc, L = θ
    # Argumentos de _wk: diagonal a, vetor b, vetor c, termo direto d = Rc e passo h.
    return _wk([-1 / (C * Rp), -Rc / L], [1, Rc], [1 / C, -Rc / L], Rc, h, discretize)
end

# G3: Windkessel de 3 elementos, θ = [Rp, C, Rc] (mL/s). Um único estado, x1 =
# volume em C. Pressão P = x1/C + Rc·Q: Rc fica em série com o par Rp‖C.
# É o caso L → ∞ (o ramo de L não conduz e todo o fluxo passa por Rc).
# Usado só na identificação (scripts/06) para comparar 2, 3 e 4 elementos.
"""
Windkessel de 3 elementos, θ = [Rp, C, Rc]. Corresponde ao de 4 elementos com L = Inf.
"""
function G3(θ, h; discretize=true)
    Rp, C, Rc = θ
    return _wk([-1 / (C * Rp)], [1], [1 / C], Rc, h, discretize)
end

# G2: Windkessel de 2 elementos (o de Frank), θ = [Rp, C] (mL/s). P = x1/C.
# É o caso L = 0 (L curto-circuita Rc). Na diástole (Q = 0) a pressão cai como
# P0·e^(−t/τ), com τ = Rp·C (constante de tempo, em s), como no Resumo da disciplina.
# scripts/01_verificacao.jl confere essa queda contra a solução exata.
"""
Windkessel de 2 elementos, θ = [Rp, C]. Corresponde ao de 4 elementos com L = 0.
"""
function G2(θ, h; discretize=true)
    Rp, C = θ
    return _wk([-1 / (C * Rp)], [1], [1 / C], 0, h, discretize)
end

# lsim_wk: simula o sistema discreto passo a passo (h = 5 ms no modelo).
# Entradas: sys (feito por G, G3 ou G2), u = fluxo em cada amostra (mL/s) e
# x0 = estado inicial (padrão: zeros, artérias "vazias").
# Saídas: y = pressão em cada amostra (mmHg) e x = matriz de estados (uma linha
# por amostra).
"""
Simula o sistema discreto `sys` com entrada `u` a partir do estado `x0`.
Retorna `(y, x)`, com `x[k, :]` o estado na amostra k (mesma convenção do `lsim` do ControlSystems).
"""
function lsim_wk(sys::WKSistema, u::AbstractVector, x0=zeros(length(sys.a)))
    # Só aceita sistema discreto (h > 0); senão, interrompe com mensagem de erro.
    sys.h > 0 || throw(ArgumentError("lsim_wk requer um sistema discreto"))
    # Escolhe o tipo numérico comum (Float64 ou número dual do ForwardDiff).
    T = promote_type(eltype(sys.a), eltype(u), eltype(x0))
    n, nx = length(u), length(sys.a)
    # Reserva espaço para os resultados (`undef` = ainda sem valores).
    y = Vector{T}(undef, n)
    x = Matrix{T}(undef, n, nx)
    xk = T.(x0)
    # Laço no tempo, amostra k = 1, 2, ..., n:
    #   1) guarda o estado atual;
    #   2) pressão = c·x + d·fluxo (dot = produto escalar);
    #   3) avança o estado um passo: x ← a.*x + b.*u ("." = elemento a elemento).
    for k in 1:n
        x[k, :] = xk
        y[k] = dot(sys.c, xk) + sys.d * u[k]
        xk = sys.a .* xk .+ sys.b .* u[k]
    end
    return y, x
end

# estado_final: mesmo laço de lsim_wk, mas devolve só o estado no fim de `u`.
# No modelo, `u` é um batimento: o estado final de um batimento é o inicial do
# próximo. acoplamento.jl encadeia os batimentos assim, enquanto Rp, VS e FC
# mudam com o propofol.
"Estado após simular `u` inteiro a partir de `x0` (útil para encadear batimentos)."
function estado_final(sys::WKSistema, u::AbstractVector, x0)
    xk = copy(x0)
    # Para cada amostra de fluxo uk, avança o estado um passo de h.
    for uk in u
        xk = sys.a .* xk .+ sys.b .* uk
    end
    return xk
end

# simulate_wk: atalho que monta o modelo G_(θ, h) e chama lsim_wk.
# Entradas: θ (mL/s), x0 (estado inicial), u (fluxo, mL/s), h (passo, s) e G_,
# que escolhe o modelo (G, G3 ou G2; padrão G = 4 elementos).
# Saída: a pressão (mmHg); com xout=true, também os estados.
"""
Simula o modelo `G_(θ, h)` com fluxo `u`, a partir de `x0`. Retorna a pressão
(e os estados, se `xout=true`).
"""
function simulate_wk(θ, x0, u, h; G_=G, xout=false)
    # `y, x = ...` recebe de uma vez os dois resultados de lsim_wk.
    y, x = lsim_wk(G_(θ, h), u, x0)
    # `cond ? a : b` é o "se/senão" em uma linha.
    return xout ? (y, x) : y
end

# solvex0: estado inicial periódico (regime estacionário). Responde à pergunta:
# "se este batimento se repetisse para sempre, em que estado estariam as artérias
# no início de cada batimento?". Assim a primeira pressão simulada já é a de
# regime, sem simular dezenas de batimentos até o transitório sumir.
# Ideia: após n passos, x_n = aⁿ·x0 + xₙ⁰ (xₙ⁰ = estado partindo do repouso).
# Impondo x_n = x0 (periodicidade): x0 = xₙ⁰ ./ (1 − aⁿ), elemento a elemento.
"""
Resolve x0 explorando a periodicidade da entrada: se `u` é um período, o estado
inicial periódico satisfaz x0 = Aⁿx0 + xₙ⁰, onde xₙ⁰ é o estado após n passos
partindo do repouso. Com A diagonal: x0ᵢ = xₙ⁰ᵢ / (1 - aᵢⁿ).
"""
function solvex0(sys::WKSistema, u::AbstractVector)
    # Estado ao fim do batimento partindo de zero (repouso).
    xn0 = estado_final(sys, u, zeros(eltype(sys.a), length(sys.a)))
    # ".^" eleva cada elemento à potência n = número de amostras do batimento.
    return xn0 ./ (1 .- sys.a .^ length(u))
end

# bode_wk: resposta em frequência do modelo contínuo, isto é, a impedância de
# entrada da árvore arterial Z(ω) = P/Q para um fluxo senoidal de frequência ω.
# Entradas: θ e um vetor w de frequências (rad/s). Saída: módulo (mmHg·s/mL) e
# fase (graus). Os testes conferem |Z| → Rp quando ω → 0 e |Z| → Rc quando ω → ∞.
"""
Resposta em frequência do modelo contínuo em `w` (rad/s): retorna `(magnitude, fase_graus)`.
Substitui o `bode` do ControlSystems usado nos scripts originais.
"""
function bode_wk(θ, w; G_=G)
    # Monta o sistema contínuo (sem discretizar).
    sys = G_(θ, 0.0; discretize=false)
    # Com A diagonal: Z(jω) = Σ cᵢ·bᵢ/(jω − aᵢ) + d. `im` é a unidade imaginária j;
    # `[... for wi in w]` calcula um valor para cada frequência (compreensão de lista).
    H = [sum(sys.c .* sys.b ./ (im * wi .- sys.a)) + sys.d for wi in w]
    # abs = módulo; angle = fase em radianos; rad2deg converte para graus.
    return abs.(H), rad2deg.(angle.(H))
end

# ------------------------------------------------------------------
# Identificação (erro de saída, Newton com diferenciação automática)
# ------------------------------------------------------------------

# Identificação de parâmetros: dados o fluxo medido u e a pressão medida y,
# procura θ = [Rp, C, Rc, L] (e, se pedido, o estado inicial x0) que deixa a
# pressão simulada o mais próxima possível da medida. "Erro de saída" =
# pressão medida − pressão simulada. O otimizador minimiza a soma dos quadrados
# desse erro (mmHg²): é um ajuste por mínimos quadrados.
# No paciente virtual, só C é ajustada, e por bisseção (calibracao.jl); identify é
# usado em scripts/06_identificacao_windkessel.jl e nos testes.
#
# IdConditions: guarda as "condições" da identificação.
#   θ0: chute inicial de θ;  x00: chute inicial do estado;
#   θid: máscara, 1 = identificar e 0 = manter fixo (ex.: [1,1,1,0] fixa L);
#   x0id: identificar também x0?;  x0solve: recalcular x0 por solvex0 a cada
#   tentativa (estado periódico);  G_: modelo (G, G3 ou G2).
struct IdConditions
    θ0::Vector{Float64}     # e.g. [Rp, C, Rc, L]
    x00::Vector{Float64}
    θid::Vector{Int}        # e.g. [1,1,1,0]
    x0id::Bool              # identificar x0
    x0solve::Bool           # resolver x0 pela periodicidade a cada passo
    G_::Function            # modelo, na forma `G_(θ, h; discretize=true)`
    # Construtor interno: confere que θ0 e θid têm o mesmo tamanho e cria o objeto
    # com `new`. G_ tem valor padrão G (4 elementos).
    function IdConditions(θ0, x00, θid, x0id, x0solve, G_=G)
        length(θ0) == length(θid) || error("Length of `θ0` != `θid`.")
        new(θ0, x00, θid, x0id, x0solve, G_)
    end
end

# p: junta num único vetor só o que vai ser ajustado (os θ marcados com 1 em θid
# e, se x0id, o estado inicial). É nesse vetor que o otimizador mexe.
"Monta o vetor de otimização."
function p(θ, x0, id::IdConditions)
    # findall devolve as posições marcadas com 1; θ[...] pega só esses valores.
    pv = θ[findall(Bool.(id.θid))]
    # vcat concatena vetores (aqui, põe x0 no fim).
    if id.x0id
        pv = vcat(pv, x0)
    end
    return Float64.(pv)
end

# deal: faz o caminho inverso de p. Recebe o vetor de otimização pv e devolve θ
# completo (com os valores fixos de θ0 onde θid = 0) e x0.
"Recupera a parametrização original (θ, x0) a partir do vetor de otimização."
function deal(pv, id::IdConditions)
    # cumsum (soma acumulada) de θid dá, para cada parâmetro identificado, a sua
    # posição em pv. A compreensão percorre (posição, valor fixo, marca) com zip e
    # escolhe pv[i] se a marca for 1, ou o valor fixo de θ0 se for 0.
    pidx = cumsum(id.θid)
    θ = [Bool(id_) ? pv[i] : θ0_ for (i, θ0_, id_) in zip(pidx, id.θ0, id.θid)]
    # Se x0 foi identificado, ele ocupa as últimas posições de pv; senão, usa x00.
    x0 = id.x0id ? pv[sum(id.θid).+(1:length(id.x00))] : id.x00
    return θ, x0
end

# cost: a função que o otimizador minimiza. Simula a pressão com os parâmetros de
# pv e devolve Σ (pressão medida − pressão simulada)², em mmHg².
# u = fluxo medido e y = pressão medida (mesmo número de amostras); h = passo (s).
"Custo (soma dos quadrados do erro de saída) a ser minimizado."
function cost(pv, id::IdConditions, u, y, h)
    θ, x0 = deal(pv, id)
    sys = id.G_(θ, h)
    # Com x0solve, o estado inicial é o periódico (solvex0), não um parâmetro livre.
    if id.x0solve
        x0 = solvex0(sys, u)
    end
    # `ym, = ...` guarda só o primeiro resultado de lsim_wk (a pressão).
    ym, = lsim_wk(sys, u, x0)
    # e = vetor de erros; dot(e, e) = soma dos quadrados dos erros.
    e = y .- ym
    return dot(e, e)
end

# identify: ajuste local pelo método de Newton (biblioteca Optim.jl).
# Newton usa a derivada primeira (gradiente) e a segunda (hessiana) do custo, que o
# ForwardDiff calcula de forma exata por diferenciação automática (sem aproximação
# por diferenças finitas). Com θpos=true usa IPNewton (Newton de pontos interiores),
# que impede parâmetros negativos (R, C e L negativos não têm sentido físico).
# Precisa de um chute inicial θ0 razoável.
# Entradas: θ0, x00, u (fluxo, mL/s), y (pressão, mmHg), h (s) e opções nomeadas
# (depois do ";"); v=true mostra o andamento. Saída: (phat, res) = vetor ótimo
# (use deal para recuperar θ) e o relatório completo do Optim.
"""
Identifica um subconjunto de {parâmetros θ, estado inicial x0} pelo método de Newton
(IPNewton com θ ≥ 0 se `θpos=true`). Retorna `(phat, res)`.
"""
function identify(θ0, x00, u, y, h; θid=ones(Int, length(θ0)), x0id=true, x0solve=false,
                  G_=G, θpos=true, v=false)
    # Monta as condições e o vetor inicial de otimização.
    id = IdConditions(θ0, x00, θid, x0id, x0solve, G_)
    p0 = p(θ0, x00, id)
    # f: custo em função só de pv. g! e h! preenchem o gradiente e a hessiana
    # (o "!" no nome é a convenção de Julia para função que altera o argumento).
    f(pv) = cost(pv, id, u, y, h)
    g!(grad, pv) = ForwardDiff.gradient!(grad, f, pv)
    function h!(hess, pv)
        ForwardDiff.hessian!(hess, f, pv)
        # Força a hessiana a ser simétrica (média com a sua transposta ').
        hess .= (hess .+ hess') ./ 2
    end

    # TwiceDifferentiable reúne custo, gradiente e hessiana para o Optim.
    # store_trace guarda o histórico das iterações.
    df = TwiceDifferentiable(f, g!, h!, p0)
    opts = Optim.Options(store_trace=true, extended_trace=true, show_trace=v)
    # `res = if ... else ... end`: em Julia o próprio if devolve um valor.
    res = if θpos
        # Limites: θ ≥ 0 (zeros) e x0 livre (−Inf), pois o estado x2 pode ser negativo;
        # limite superior +Inf para todos.
        lp = vcat(zeros(sum(id.θid)), fill(-Inf, id.x0id ? length(x00) : 0))
        up = fill(Inf, length(p0))
        dfc = TwiceDifferentiableConstraints(lp, up)
        optimize(df, dfc, p0, IPNewton(), opts)
    else
        # Sem θpos: Newton sem restrição de sinal, com busca em linha por retrocesso.
        optimize(df, p0, Newton(alphaguess=LineSearches.InitialPrevious(),
                                linesearch=LineSearches.BackTracking()), opts)
    end
    # Optim.minimizer(res) = o melhor vetor encontrado.
    return Optim.minimizer(res), res
end

# identifyglobal: mesmo custo, mas com Nelder-Mead (método do simplex, que não usa
# derivadas). Tolera melhor um chute inicial ruim, mas não garante o mínimo global
# (o nome vem do original) e não impõe θ ≥ 0. Até 10^5 iterações.
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

# params_stergiopulos1999total: valores de referência (humano, batimento tipo A)
# de Stergiopulos, Westerhof e Westerhof 1999, Am J Physiol 276:H81, Tabela 1,
# copiados do windkessel-id (oeid.jl) através do PropofolPBL do grupo.
# Veredito em docs/validacao.md, seção 4.3: CONFIRMADO* (pelo código que cita a
# tabela; abrir o artigo ainda é pendência, seção 6.1).
# No modelo: C, Rc e L valem para todo paciente (C é depois recalibrada com
# PAS − PAD em calibracao.jl) e Rp é trocado por PAM/DC do paciente.
# Saída: (θ, x0), com θ = [Rp, C, Rc, L] em unidades de mL/s e x0 = estado inicial.
"""
Parâmetros de Stergiopulos et al. (1999), batimento tipo A, humano.
"""
function params_stergiopulos1999total()
    # Ordem de grandeza: Rp = 0,79 mmHg·s/mL ≈ 1 050 dyn·s·cm⁻⁵, dentro da faixa
    # normal da resistência vascular sistêmica (800–1 200); Rc/Rp ≈ 7%.
    Rp = 0.79     # resistência periférica [mmHg·s/mL]
    C = 1.22      # complacência arterial total [mL/mmHg]
    Rc = 0.056    # resistência característica da aorta [mmHg·s/mL]
    L = 0.0051    # inertância arterial total [mmHg·s²/mL]
    # x0 = [volume em C, L × fluxo em L] é só um chute inicial: o modelo recalcula o
    # estado periódico com solvex0 (docs/validacao.md, seção 3, item 11).
    # 8.066 * 1000 / 60 converte 8,066 de unidades L/min para mL/s.
    # estado inicial do windkessel-id v1.0 estava em L/min; aqui tudo é mL/s, então
    # converte como na v2.0 do windkessel-id (src/WindkesselIdentification.jl#L218)
    return [Rp, C, Rc, L], [8.066 * 1000 / 60, -8.296]
end
