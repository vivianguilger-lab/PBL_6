# =============================================================================
# populacao.jl — coorte virtual (virtual clinical trial) e variabilidade
# =============================================================================
#
# EXPERIMENTAL — NÃO VALIDADO. As variabilidades abaixo e as distribuições da
# coorte são suposições; os resultados não são usados na plataforma nem nos
# resultados do README (pendências em docs/validacao.md). Leitura das variâncias
# de Schnider (ω² ou coeficiente de variação) ainda a confirmar na tabela do artigo.
#
# Variabilidade interindividual (IIV) log-normal: θᵢ = θ·exp(ηᵢ), ηᵢ ~ N(0, ω²)
#   PK  — variâncias de Schnider et al. 1998 (V1, V2, V3, CL, Q2, Q3); ke0: ω = 0,3 (assumido)
#   Hemodinâmica — Ce50 e Imax do EfeitoHemodinamico do grupo: ω = 0,3 (assumido,
#   os próprios valores centrais são ilustrativos no modelo do grupo)
#
# Como ler θᵢ = θ·exp(ηᵢ): cada indivíduo tem o parâmetro típico multiplicado
# por um fator sorteado. Com ω = 0,3 o fator fica tipicamente entre ~0,74 e
# ~1,35 (±1 desvio-padrão). A forma log-normal mantém os parâmetros positivos.
# Este arquivo também gera uma coorte virtual (ensaio clínico virtual, script
# 04_coorte_virtual.jl) com distribuições de idade, IMC, PA, ASA, PVI e IVCI
# INVENTADAS por nós (docs/validacao.md, seção 6.3). Nada aqui está validado.
# =============================================================================

# OMEGA_PK: desvio-padrão ω de η para cada parâmetro da PK. `sqrt` (raiz
# quadrada) porque os números de Schnider foram lidos como variâncias ω²; se
# forem coeficientes de variação, mudam muito (até 25× de diferença,
# docs/validacao.md, seção 6.3). ke0: ω = 0,3 assumido.
# OMEGA_HEMO: ω = 0,3 para Ce50 e Imax, assumido.
# `const` = valor global fixo; `(V1 = ..., ...)` é um NamedTuple.
const OMEGA_PK = (V1 = sqrt(0.0404), V2 = sqrt(0.01), V3 = sqrt(0.1435),
                  CL = sqrt(0.1005), Q2 = sqrt(0.01), Q3 = sqrt(0.1179), ke0 = 0.3)
const OMEGA_HEMO = (Ce50 = 0.3, Imax = 0.3)

# amostrar_eta_pk (EXPERIMENTAL): sorteia os fatores exp(η) da PK, um por
# parâmetro, com o gerador `rng`. `randn` sorteia da normal padrão N(0, 1);
# vezes ω dá N(0, ω²). `NamedTuple{keys(OMEGA_PK)}(...)` devolve os fatores
# com os mesmos nomes (V1, V2, ..., ke0), lidos por `parametros_pk` (propofol.jl).
"""Fatores multiplicativos exp(η) da PK, η ~ N(0, ω²)."""
amostrar_eta_pk(rng) = NamedTuple{keys(OMEGA_PK)}(Tuple(exp(w * randn(rng)) for w in values(OMEGA_PK)))

# amostrar_efeito (EXPERIMENTAL): cria um EfeitoHemodinamico (acoplamento.jl)
# de um indivíduo sorteado em torno do típico `ef`: Ce50 × exp(η) e cada Imax
# × exp(η), com ω = 0,3 assumido; γ não varia. Cada Imax é limitado a 0,9
# (`min`) para não sortear efeitos absurdos. Os valores centrais já são
# ilustrativos (docs/validacao.md, seção 4.4).
"""Efeito hemodinâmico de um indivíduo sorteado em torno de `ef`."""
function amostrar_efeito(rng, ef::EfeitoHemodinamico = EfeitoHemodinamico())
    # `fI()` sorteia um novo fator exp(η) a cada chamada (um para cada Imax).
    fI() = exp(OMEGA_HEMO.Imax * randn(rng))
    return EfeitoHemodinamico(Ce50 = ef.Ce50 * exp(OMEGA_HEMO.Ce50 * randn(rng)), γ = ef.γ,
                              Imax_Rp = min(ef.Imax_Rp * fI(), 0.9), Imax_VS = min(ef.Imax_VS * fI(), 0.9),
                              Imax_FC = min(ef.Imax_FC * fI(), 0.9))
end

# _trunc: sorteia da normal com média μ e desvio σ, truncada em [a, b]: repete
# o sorteio até cair dentro da faixa. Usado só na coorte (valores assumidos).
_trunc(rng, μ, σ, a, b) = (x = μ + σ * randn(rng); while !(a <= x <= b); x = μ + σ * randn(rng); end; x)

# Individuo: um paciente da coorte virtual. `struct` define um tipo de dado com
# campos fixos: dados pré-indução (`Paciente`), fatores da PK sorteados, efeito
# hemodinâmico sorteado e dose do protocolo (mg/kg).
"Um indivíduo da coorte: dados pré-indução + variabilidade + dose do protocolo."
struct Individuo
    paciente::Paciente
    eta_pk::NamedTuple
    efeito::EfeitoHemodinamico
    dose_mgkg::Float64
end

# dose_padrao: dose de indução do protocolo da coorte, por idade (mg/kg). Na
# ordem da bula do Diprivan (2–2,5 mg/kg no adulto saudável; 1–1,5 no idoso;
# docs/validacao.md, seção 4.7), mas os cortes de 55 e 70 anos são escolha nossa.
"""
    dose_padrao(idade) -> mg/kg

"Faixa por peso, ajustada por idade": 2,0 mg/kg (< 55 anos), 1,5 (55–69), 1,2 (≥ 70).
"""
dose_padrao(idade) = idade < 55 ? 2.0 : (idade < 70 ? 1.5 : 1.2)

# -----------------------------------------------------------------------------
# gerar_coorte (EXPERIMENTAL, NÃO VALIDADO): cria `n` indivíduos sintéticos.
# Todas as médias, desvios e faixas abaixo são SUPOSIÇÕES nossas, sem fonte
# (docs/validacao.md, seção 6.3). `variabilidade_dose = 0.1` faz a dose variar
# ±10% em torno da dose padrão. Semente fixa (Xoshiro(2026)): mesma coorte a
# cada execução.
# -----------------------------------------------------------------------------
"""
    gerar_coorte(n; rng=Xoshiro(2026))

`n` adultos virtuais com demografia, PA/FC pré-indução, ASA, PVI, IVCI e IIV de
PK e da resposta hemodinâmica. Experimental: as distribuições são suposições, sem
fonte (docs/validacao.md).
"""
function gerar_coorte(n; rng = Xoshiro(2026), variabilidade_dose = 0.1)
    # `Individuo[]` cria uma lista vazia de Individuo; `push!` acrescenta ao final.
    coorte = Individuo[]
    for i in 1:n
        # Demografia: idade ~ N(55, 16) entre 18 e 90 anos; sexo 50/50; altura (cm)
        # por sexo; peso (kg) = IMC × altura² (IMC ~ N(27; 4,5) entre 17 e 40 kg/m²).
        idade = _trunc(rng, 55.0, 16.0, 18.0, 90.0)
        sexo = rand(rng) < 0.5 ? :M : :F
        altura = sexo == :M ? _trunc(rng, 174.0, 7.0, 150.0, 200.0) : _trunc(rng, 162.0, 6.5, 140.0, 190.0)
        peso = _trunc(rng, 27.0, 4.5, 17.0, 40.0) * (altura / 100)^2
        # PA e FC pré-indução (mmHg, bpm): a PAS sobe com a idade; a PAD segue uma
        # curva suave em idade; `clamp` força a pressão de pulso (PAS − PAD) entre 30 e 100.
        pas = _trunc(rng, 115.0 + 0.45 * (idade - 18), 16.0, 90.0, 210.0)
        pad = clamp(_trunc(rng, 76.0 + 0.15 * (idade - 18) - 0.004 * (idade - 18)^2, 9.0, 45.0, 115.0), pas - 100, pas - 30)
        fc = _trunc(rng, 76.0, 12.0, 48.0, 115.0)
        # ASA sorteado por faixa etária (`u` uniforme entre 0 e 1): < 40 anos 60% ASA 1
        # e 40% ASA 2; 40–64 anos 20/60/20% ASA 1/2/3; >= 65 anos 55/40/5% ASA 2/3/4.
        # PVI (%) e IVCI (%) por normais truncadas.
        u = rand(rng)
        asa = idade < 40 ? (u < 0.6 ? 1 : 2) : idade < 65 ? (u < 0.2 ? 1 : u < 0.8 ? 2 : 3) : (u < 0.55 ? 2 : u < 0.95 ? 3 : 4)
        pvi = _trunc(rng, 12.0, 5.0, 3.0, 35.0)
        ivci = _trunc(rng, 25.0, 13.0, 0.0, 80.0)
        # dados sintéticos: a PAM do indivíduo virtual é gerada a partir de PAS e PAD
        p = Paciente(nome = "virtual_$(i)", idade = idade, peso = peso, altura = altura, sexo = sexo, asa = asa,
                     pas = pas, pad = pad, pam = pad + (pas - pad) / 3, fc = fc, pvi = pvi, ivci = ivci)
        # Dose do protocolo × fator uniforme entre 0,9 e 1,1 (`rand` é uniforme 0–1).
        dose = dose_padrao(idade) * (1 + variabilidade_dose * (2 * rand(rng) - 1))
        push!(coorte, Individuo(p, amostrar_eta_pk(rng), amostrar_efeito(rng), dose))
    end
    return coorte
end

# calibrar_coorte: cria o paciente virtual de cada indivíduo (calibracao.jl), já
# com os fatores da PK e o efeito sorteados. `Vector{...}(undef, n)` reserva n
# posições vazias; `Threads.@threads` reparte as calibrações entre os núcleos
# (cada i preenche só a sua posição, então não há conflito).
"Calibra cada paciente virtual com os SEUS dados pré-indução (em paralelo com `julia -t auto`)."
function calibrar_coorte(coorte)
    pvs = Vector{PacienteVirtual}(undef, length(coorte))
    Threads.@threads for i in eachindex(coorte)
        pvs[i] = calibrar(coorte[i].paciente; eta_pk = coorte[i].eta_pk, efeito = coorte[i].efeito)
    end
    return pvs
end

# -----------------------------------------------------------------------------
# simular_coorte (EXPERIMENTAL): simula a indução de cada indivíduo e guarda
# uma linha por indivíduo (dados, s, PAM basal e mínima, queda %, tempo abaixo
# de 65 e 55 mmHg, BIS mínimo, risco e se a calibração convergiu).
# `pvs = calibrar_coorte(coorte)` é valor padrão: se não for dado, calibra.
# `dose_fixa`: se informada, todos recebem a mesma dose (mg/kg); se `nothing`
# (o "vazio" de Julia), cada um recebe a sua dose do protocolo.
# -----------------------------------------------------------------------------
"""
    simular_coorte(coorte, pvs=calibrar_coorte(coorte); manutencao_ugkgmin=100.0, duracao_bolus=30.0,
                   janela=20.0, dose_fixa=nothing)
"""
function simular_coorte(coorte, pvs = calibrar_coorte(coorte); manutencao_ugkgmin = 100.0, duracao_bolus = 30.0,
                        janela = 20.0, dose_fixa = nothing)
    # `Vector{Any}` aceita qualquer tipo; cada thread preenche a sua posição.
    saida = Vector{Any}(undef, length(coorte))
    Threads.@threads for i in eachindex(coorte)
        ind, pv = coorte[i], pvs[i]
        dose = dose_fixa === nothing ? ind.dose_mgkg : dose_fixa
        mt = metricas(simular_inducao(pv; dose_mgkg = dose, duracao_bolus, manutencao_ugkgmin, duracao = janela); janela)
        p = ind.paciente
        saida[i] = (idade = p.idade, sexo = p.sexo, peso = p.peso, asa = p.asa, pvi = p.pvi, ivci = p.ivci,
                    s = pv.s, pam_basal = mt.pam_basal, dose_mgkg = dose, pam_nadir = mt.pam_nadir,
                    queda_pct = mt.queda_pct, t_nadir = mt.t_nadir,
                    min_abaixo_65 = mt.tempo_abaixo[65.0], min_abaixo_55 = mt.tempo_abaixo[55.0],
                    bis_min = mt.bis_min, risco = mt.risco, calibrou = pv.convergiu)
    end
    # Copia para um vetor com o tipo dos elementos definido (mais fácil de usar depois).
    return [s for s in saida]
end

# resumo_coorte (EXPERIMENTAL): incidências em % da coorte (PAM < 65, PAM < 55,
# queda > 30%, BIS < 60, calibrados), quartis da queda e mediana do tempo do
# nadir. `count(c, res)` conta as linhas que cumprem `c`; `getfield(r, x)` lê o
# campo de nome `x` (ex.: :queda_pct). "Queda > 30%" é só descritivo: saiu da
# regra de risco por não ter fonte.
"""Resumo da coorte: incidências e estatísticas."""
function resumo_coorte(res)
    n = length(res)
    f(c) = 100 * count(c, res) / n
    q(x, p) = quantile([getfield(r, x) for r in res], p)
    return (n = n, inc_pam65 = f(r -> r.min_abaixo_65 > 0), inc_pam55 = f(r -> r.min_abaixo_55 > 0),
            inc_queda30 = f(r -> r.queda_pct > 30),
            queda_mediana = q(:queda_pct, 0.5), queda_p25 = q(:queda_pct, 0.25), queda_p75 = q(:queda_pct, 0.75),
            t_nadir_mediano = q(:t_nadir, 0.5), bis_ok = f(r -> r.bis_min < 60), calibrados = f(r -> r.calibrou))
end
