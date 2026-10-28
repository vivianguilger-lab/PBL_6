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
# e de `sensibilidade_hemodinamica` são ilustrativos, não ajustados a dados.             # PROBLEMA! Valores ilustrativos devem ser substituidos por algo mais justificado
# =============================================================================
#
# GUIA DE LEITURA (comentários acrescentados para o grupo; o código não mudou)
# -----------------------------------------------------------------------------
# Este arquivo é a "ponte" entre a farmacologia e a circulação:
#   1. A PK/PD (src/propofol.jl) entrega Ce(t), a concentração de propofol no
#      sítio efetor, em µg/ml, numa grade de 1 s (tempos em minutos).
#   2. Aqui, o Ce reduz três variáveis: Rp (vasodilatação), VS (depressão do
#      miocárdio e queda da pré-carga) e FC. Ver `fator_efeito`.
#   3. Com Rp, VS e FC de cada batimento, montamos o fluxo que sai do ventrículo
#      e calculamos a pressão com o Windkessel de 4 elementos (src/windkessel.jl).
#   4. De cada batimento extraímos PAS (máximo), PAD (mínimo) e PAM (média).
# Unidades internas da circulação: mmHg, ml e s (como no windkessel-id).
# Rp e Rc em mmHg·s/ml, C em ml/mmHg, L em mmHg·s²/ml, fluxo em ml/s.
# A origem e o veredito de cada número estão em docs/validacao.md (seção 4).

# -----------------------------------------------------------------------------
# Parâmetros do efeito do propofol sobre a circulação (modelo Emax sigmoide)
# -----------------------------------------------------------------------------
# Sintaxe Julia: `struct` cria um "registro" com campos nomeados. `Base.@kwdef`
# permite criá-lo informando, por nome, só os campos que se quer mudar; os demais
# usam o valor padrão (depois do `=`). Ex.: EfeitoHemodinamico(Imax_Rp = 0.30).
# `::Float64` diz que o campo é um número real (com casas decimais).
# Origem: valores do grupo (PropofolPBL, src/acoplamento.jl), ILUSTRATIVOS, sem
# fonte nem calibração com pacientes (ver docs/validacao.md, seção 4.4).

"""
Efeito do propofol sobre a hemodinâmica. Cada variável X é multiplicada por
`1 − Imax_X · s · Ceᵞ/(Ce50ᵞ + Ceᵞ)` (mínimo 0,05), com `s` a sensibilidade
do paciente.
"""
Base.@kwdef struct EfeitoHemodinamico  # como mencionado, todos os valores precisam ser verificados e mudados p/ dados reais, pois são ilustrativos
    # Ce50: concentração no sítio efetor (µg/ml) que produz metade do efeito
    # máximo. 3,0 µg/ml está na ordem do EC50 da resistência periférica de
    # Su 2022 (2,96 µg/ml), mas foi escolhido pelo grupo (ILUSTRATIVO).
    Ce50::Float64 = 3.0       # µg/ml
    # γ (gama): coeficiente de Hill, sem unidade. Quanto maior, mais "em degrau"
    # é a curva concentração-efeito em torno do Ce50. ILUSTRATIVO (o PAS usa 1,83).
    γ::Float64 = 2.0
    # Imax_*: fração máxima de queda de cada variável (0,25 = até 25% a menos,
    # com s = 1). Rp cai pela vasodilatação arterial do propofol.
    Imax_Rp::Float64 = 0.25   # vasodilatação arterial
    # VS cai até 15%: representa ao mesmo tempo a depressão do miocárdio e a
    # redução da pré-carga (venodilatação), pois não há retorno venoso no modelo.
    Imax_VS::Float64 = 0.15   # depressão miocárdica / redução de pré-carga
    # FC cai até 5% (efeito pequeno). Não há barorreflexo no modelo, então a
    # queda de pressão não provoca taquicardia compensatória.
    Imax_FC::Float64 = 0.05
end

# Acho bem importante mudar isso!!
# Além disso, o modelo não possui barorreflexo, porque é uma simplificação estrutural do modelo.
# O código assume que: propofol → PA ↓
# Mas não permite: PA ↓ → barorreflexo → FC ↑
# Portanto, a FC só pode cair pelo efeito definido.

# -----------------------------------------------------------------------------
# fator_efeito: quanto "sobra" de uma variável para um dado Ce
# -----------------------------------------------------------------------------
# Entradas: `ef` (parâmetros acima), `Imax` (o da variável), `s` (sensibilidade
# do paciente, sem unidade) e `Ce` (µg/ml). Saída: número entre 0,05 e 1 que
# multiplica o valor basal (ex.: Rp = Rp0 · fator).
# H(Ce) = Ceᵞ / (Ce50ᵞ + Ceᵞ) é a curva de Hill: vai de 0 a 1 e vale 0,5 no Ce50.
# Exemplo (s = 1, Imax 0,25): Ce = 3 µg/ml → H = 0,5 → Rp cai 12,5% (fator
# 0,875); Ce = 6 µg/ml → H = 0,8 → Rp cai 20% (fator 0,80).
# `max(..., 0.05)` impede valores nulos ou negativos: piso de 5% do basal,
# ILUSTRATIVO, que nunca é atingido com os Imax atuais. `^` é potência.
# Sintaxe Julia: `f(x) = expressão` define uma função curta, sem `function`.
"Fator multiplicativo 1 − Imax·s·H(Ce), limitado a 0,05."
fator_efeito(ef::EfeitoHemodinamico, Imax, s, Ce) =
    max(1 - Imax * s * Ce^ef.γ / (ef.Ce50^ef.γ + Ce^ef.γ), 0.05)  # Também ilustrativo!

# -----------------------------------------------------------------------------
# Fator de sensibilidade s do paciente
# -----------------------------------------------------------------------------
# s multiplica o Imax: s = 1 é o paciente de referência; com s = 1,5 o propofol
# derruba Rp, VS e FC 50% mais. Cada termo só conta acima de um limiar:
#   idade acima de 50 anos:  +0,01 por ano       (k_idade)
#   ASA acima de 2:          +0,15 por classe    (k_asa; ASA III soma 0,15)
#   PVI acima de 15%:        +0,03 por ponto %   (k_pvi; limiar de Tsuchiya 2010)
#   IVCI acima de 43%:       +0,01 por ponto %   (k_ivci; Zhang & Critchley 2016)
# PVI = índice de variabilidade pletismográfica; IVCI = índice de colapsabilidade
# da veia cava inferior. Valores altos sugerem hipovolemia (dependência de
# pré-carga), que piora a hipotensão da indução.
# Limiares com fonte; coeficientes k_* do grupo, ILUSTRATIVOS (docs/validacao.md,
# seção 4.4; o grupo usava 13% e 40%, corrigidos na seção 3, item 6).
# Exemplo (paciente C do README: 81 anos, ASA III, PVI 18%, IVCI 45%):
# s = 1 + 0,31 + 0,15 + 0,09 + 0,02 = 1,57.
# Sintaxe Julia: os argumentos depois do `;` são nomeados e opcionais
# (ex.: sensibilidade_hemodinamica(p; k_asa = 0.2)). `max(x, 0)` zera a parte
# negativa: abaixo do limiar o termo não soma nada. `pac.idade` lê o campo
# `idade` do registro `pac` (struct Paciente, em src/paciente.jl).
"""
    sensibilidade_hemodinamica(pac; k_idade=0.01, k_asa=0.15, k_pvi=0.03, k_ivci=0.01)

Fator de sensibilidade hemodinâmica (1 = paciente de referência). Aumenta com
idade > 50 anos, ASA ≥ 3, PVI > 15% e IVCI > 43%, fatores associados à hipotensão
pós-indução. Limiares: PVI de Tsuchiya et al. 2010 (PVI > 15 prevê queda da PAM após
propofol); IVCI de Zhang & Critchley 2016 (corte ótimo de 43%). Os coeficientes `k_*`
são do grupo, ilustrativos, ainda sem calibração (docs/validacao.md).
"""
function sensibilidade_hemodinamica(pac::Paciente; k_idade = 0.01, k_asa = 0.15, k_pvi = 0.03, k_ivci = 0.01)
    return 1.0 + k_idade * max(pac.idade - 50, 0) + k_asa * max(pac.asa - 2, 0) +        # ILUSTRATIVO TAMBÉM! 
           k_pvi * max(pac.pvi - 15, 0) + k_ivci * max(pac.ivci - 43, 0)
end

# -----------------------------------------------------------------------------
# θ basal: os 4 parâmetros do Windkessel antes do propofol
# -----------------------------------------------------------------------------
# θ (teta) = [Rp, C, Rc, L]:
#   Rp: resistência periférica (arteríolas), mmHg·s/ml;
#   C:  complacência (elasticidade) das grandes artérias, ml/mmHg;
#   Rc: resistência característica da aorta proximal, mmHg·s/ml;
#   L:  inertância (inércia do sangue na aorta), mmHg·s²/ml.
# Rp vem do próprio paciente; C = 1,22, Rc = 0,056 e L = 0,0051 são de
# Stergiopulos 1999 (via código do windkessel-id). Em calibracao.jl a C é
# reajustada para reproduzir a pressão de pulso do paciente.
"""
    theta_basal(pac)

Parâmetros basais θ = [Rp, C, Rc, L] (unidades mL/s) do Windkessel de 4
elementos: Rp = PAM/DC (em regime, Z(0) = Rp); C, Rc e L de Stergiopulos (1999).
"""
function theta_basal(pac::Paciente)
    # Copia os valores de Stergiopulos (src/windkessel.jl). A função devolve dois
    # resultados (θ e um estado inicial); `θ, _ =` guarda o 1º e descarta o 2º.
    θ, _ = params_stergiopulos1999total()
    # Rp = PAM / DC, a "lei de Ohm" da circulação: em regime, a pressão média
    # depende só de Rp e do fluxo médio (PAM = Rp·Q̄). DC vem em L/min e vira
    # ml/s (× 1000 ml/L ÷ 60 s/min). θ[1] é o 1º elemento (em Julia os índices
    # começam em 1). Ex.: PAM 93 mmHg e DC 5 L/min (83,3 ml/s) → Rp ≈ 1,12
    # mmHg·s/ml. Mesmo cálculo do grupo; ver docs/validacao.md.
    θ[1] = pam_basal(pac) / (debito_basal(pac) * 1000 / 60)
    return θ
end

# -----------------------------------------------------------------------------
# Forma de onda do fluxo aórtico (a entrada do Windkessel)
# -----------------------------------------------------------------------------
# O modelo não simula o ventrículo: ele IMPÕE o fluxo que sai para a aorta.
# Usamos o batimento humano publicado por Stergiopulos 1999 (Fig. 4A), que vem
# com o windkessel-id (dados/windkessel-id/stergiopulos1999_data/aflow.csv:
# 166 amostras a cada 5 ms, cerca de 0,83 s). Só a FORMA é usada; a duração e a
# amplitude são ajustadas à FC e ao VS de cada batimento.
# `_FORMA_FLUXO` é uma "caixa" (`Ref`) que guarda a forma depois da primeira
# leitura, para não reler o arquivo a cada batimento. `const` = nome fixo.
const _FORMA_FLUXO = Ref{Vector{Float64}}()

"""
Forma de onda normalizada do fluxo aórtico (média 1 em τ ∈ [0, 1)), a partir
do batimento de Stergiopulos (1999) distribuído com o windkessel-id.
"""
function forma_fluxo()
    # `!isassigned(...)` = "a caixa ainda está vazia?" (`!` significa "não").
    if !isassigned(_FORMA_FLUXO)
        # Lê tempo, pressão e fluxo (ml/s); os `_` descartam tempo e pressão.
        _, _, q = dados_stergiopulos1999()
        # Divide pelo fluxo médio: a forma passa a ter média 1 (sem unidade).
        # O `.` em `./` aplica a divisão a cada elemento do vetor ("broadcasting").
        # `_FORMA_FLUXO[] = ...` guarda o vetor dentro da caixa; `[]` lê de volta.
        _FORMA_FLUXO[] = q ./ (sum(q) / length(q))
    end
    return _FORMA_FLUXO[]
end

# -----------------------------------------------------------------------------
# Fluxo de um batimento com a FC e o VS do momento
# -----------------------------------------------------------------------------
# Entradas: `forma` (m amostras, média 1), `n` (amostras do batimento, n = 60/FC/h),
# `vs` (volume sistólico, ml) e `h` (passo de tempo, s; 0,005 s = 5 ms).
# Saída: vetor com n valores de fluxo aórtico, em ml/s.
# Ideia: "esticar" ou "encolher" o batimento de referência para durar 60/FC
# segundos e escalá-lo para que a área sob a curva (volume ejetado) seja vs.
# Assim, FC menor → batimento mais longo; VS menor → curva mais baixa.
"Reamostra `forma` (um período) para `n` pontos e escala para ejetar `vs` ml em `n·h` s."
function fluxo_batimento(forma, n, vs, h)
    m = length(forma)
    # Para cada amostra k do novo batimento (k = 0 … n−1), acha a posição
    # equivalente `s` no batimento de referência (aqui `s` é só uma posição, não
    # a sensibilidade) e interpola em linha reta entre as amostras vizinhas i e
    # i+1, com peso f. `map(...) do k ... end` aplica o bloco a cada k e junta os
    # resultados num vetor. `floor(Int, s)` arredonda para baixo.
    q = map(0:n-1) do k
        s = k * m / n
        i = floor(Int, s)
        f = s - i
        # Índices em Julia começam em 1, por isso `i+1`. `mod(i + 1, m) + 1`
        # volta à 1ª amostra depois da última (o batimento é periódico).
        (1 - f) * forma[i+1] + f * forma[mod(i + 1, m)+1]
    end
    # Escala: sum(q)·h ≈ integral do fluxo no batimento = volume ejetado (ml).
    # Multiplicar por vs/(sum(q)·h) faz esse volume ser exatamente vs.
    # `.*` multiplica cada elemento do vetor pelo número.
    return q .* (vs / (sum(q) * h))
end

# -----------------------------------------------------------------------------
# Interpolação linear (ler um valor entre dois pontos de uma tabela)
# -----------------------------------------------------------------------------
# `tx` são os tempos da tabela (em ordem crescente) e `yx` os valores; devolve o
# valor em `t` numa reta entre os dois pontos vizinhos. Usada para ler Ce e BIS
# (calculados numa grade de 1 s, em minutos) no instante de cada batimento.
"Interpolação linear de (tx, yx) em t, constante fora do intervalo."
function interp(tx, yx, t)
    # Antes do início ou depois do fim da tabela, repete o primeiro/último valor.
    # `condição && return valor` = "se a condição for verdadeira, devolve o valor".
    t <= tx[1] && return yx[1]
    t >= tx[end] && return yx[end]
    # `searchsortedlast` acha o último i com tx[i] ≤ t (busca rápida em vetor
    # ordenado); f, de 0 a 1, é a fração do caminho entre tx[i] e tx[i+1].
    i = searchsortedlast(tx, t)
    f = (t - tx[i]) / (tx[i+1] - tx[i])
    return (1 - f) * yx[i] + f * yx[i+1]
end

# -----------------------------------------------------------------------------
# Batimento basal em regime permanente
# -----------------------------------------------------------------------------
# Entradas: θ = [Rp, C, Rc, L], `fc` (bpm), `vs` (ml), `h` (passo, s) e a forma
# do fluxo. Saída: `pressao` (mmHg, uma amostra a cada h), `fluxo` (ml/s) e `x0`,
# o estado interno do Windkessel no início do batimento.
# "Regime permanente" = o paciente vinha batendo assim há muito tempo: o
# batimento termina exatamente no mesmo estado em que começou. É a situação
# pré-indução. Passo h = 5 ms, como no windkessel-id (docs/validacao.md).
"""
    batimento_basal(θ, fc, vs; h=0.005)

Batimento periódico (regime permanente) do Windkessel θ com FC e VS dados.
Devolve `(pressao, fluxo, x0)`.
"""
function batimento_basal(θ, fc, vs; h = 0.005, forma = forma_fluxo())
    # Número de amostras de um batimento: duração (60/FC s) dividida pelo passo.
    # Ex.: FC 75 bpm → 0,8 s / 0,005 s = 160 amostras. `round(Int, ...)` arredonda
    # para um número inteiro.
    n = round(Int, 60 / fc / h)
    u = fluxo_batimento(forma, n, vs, h)
    # G monta o Windkessel de 4 elementos discretizado (src/windkessel.jl);
    # solvex0 calcula o estado inicial periódico (início = fim do batimento).
    sys = G(θ, h)
    x0 = solvex0(sys, u)
    # lsim_wk simula o batimento; a vírgula em `y, =` guarda só o primeiro
    # resultado (a pressão) e descarta o segundo (os estados).
    y, = lsim_wk(sys, u, x0)
    return y, u, x0
end

# -----------------------------------------------------------------------------
# Simulação batimento a batimento durante a anestesia (o centro do acoplamento)
# -----------------------------------------------------------------------------
# Entradas: o paciente `pac`, a saída da PK/PD `pkpd` (tempos `t` em min, `Ce` em
# µg/ml, `bis`), os parâmetros do efeito, o θ basal, a sensibilidade `s`, o passo
# h (s) e o tempo final `tfim` (min). Saída: uma tabela com uma linha por
# batimento (tempo em min, Ce, BIS, FC em bpm, VS em ml, DC em L/min, Rp, e PAS,
# PAD e PAM em mmHg) e, se pedido, a onda de pressão e fluxo amostra a amostra.
# Em cada batimento: (1) lê o Ce do instante; (2) calcula FC, VS e Rp com o
# efeito Emax; (3) gera o fluxo; (4) roda o Windkessel a partir do estado em que
# o batimento anterior terminou; (5) guarda PAS, PAD e PAM.
# C, Rc e L não mudam com o propofol neste modelo; só Rp, VS e FC.
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
    # Desempacota θ basal em quatro nomes (Rp0 = Rp antes do propofol).
    Rp0, C, Rc, L = θbasal
    # VS basal (ml/batimento) = DC (L/min) × 1000 (ml/L) ÷ FC (batimentos/min).
    # Ex.: 5 L/min a 75 bpm → 66,7 ml. Igual ao grupo (docs/validacao.md).
    vs0 = debito_basal(pac) * 1000 / pac.fc     # ml por batimento
    # Estado inicial: o batimento basal em regime permanente (pré-indução).
    # `; h, forma` é um atalho de Julia para `; h = h, forma = forma`.
    _, _, x = batimento_basal(θbasal, pac.fc, vs0; h, forma)

    # Colunas da tabela de saída: um dicionário (`Dict`) que liga cada nome
    # (`:t_min`, `:Ce`, ...; o `:` marca um "símbolo", isto é, um rótulo) a um
    # vetor vazio (`Float64[]`), preenchido batimento a batimento.
    # ts, ps, qs guardam a onda completa (tempo em s, pressão em mmHg, fluxo em ml/s).
    col = Dict(k => Float64[] for k in (:t_min, :Ce, :BIS, :FC, :VS, :DC, :Rp, :PAS, :PAD, :PAM))
    ts, ps, qs = Float64[], Float64[], Float64[]
    t = 0.0                                      # s
    # Laço principal: repete enquanto o tempo (s) não chega ao fim (tfim em
    # min × 60). Cada volta do laço é um batimento.
    while t < tfim * 60
        # (1) Ce no início do batimento; pkpd.t está em minutos, por isso t/60.
        Ce = interp(pkpd.t, pkpd.Ce, t / 60)
        # (2) Efeito do propofol: cada variável = valor basal × fator_efeito.
        # FC (bpm), VS (ml) e Rp (mmHg·s/ml) caem com o Ce, e caem mais se s > 1.
        fc = pac.fc * fator_efeito(efeito, efeito.Imax_FC, s, Ce)
        vs = vs0 * fator_efeito(efeito, efeito.Imax_VS, s, Ce)
        Rp = Rp0 * fator_efeito(efeito, efeito.Imax_Rp, s, Ce)
        # (3) Duração do batimento em amostras (FC menor → batimento mais longo)
        # e o fluxo aórtico desse batimento, escalado para ejetar vs ml.
        n = round(Int, 60 / fc / h)
        u = fluxo_batimento(forma, n, vs, h)
        # (4) Windkessel com o Rp do momento (C, Rc e L fixos). A simulação parte
        # do estado `x` deixado pelo batimento anterior, por isso a pressão é
        # contínua de um batimento para o outro; depois `x` é atualizado.
        sys = G([Rp, C, Rc, L], h)
        y, = lsim_wk(sys, u, x)
        x = estado_final(sys, u, x)
        # (5) Guarda uma linha da tabela: tempo (min), Ce, BIS no mesmo instante,
        # FC, VS, DC = VS·FC/1000 (L/min), Rp, PAS = máximo da pressão no
        # batimento, PAD = mínimo e PAM = média das amostras (média no tempo).
        # `zip` percorre juntas as duas listas (nomes e valores); `push!`
        # acrescenta o valor ao fim do vetor da coluna correspondente.
        for (k, v) in zip((:t_min, :Ce, :BIS, :FC, :VS, :DC, :Rp, :PAS, :PAD, :PAM),
                          (t / 60, Ce, interp(pkpd.t, pkpd.bis, t / 60), fc, vs, vs * fc / 1000, Rp,
                           maximum(y), minimum(y), mean(y)))
            push!(col[k], v)
        end
        # Se a onda completa foi pedida, acrescenta os tempos (s) de cada amostra,
        # a pressão e o fluxo. `t .+ h .* (0:n-1)` gera t, t+h, t+2h, ... (o `.`
        # aplica a operação a cada elemento). `append!` junta um vetor ao fim de outro.
        if forma_de_onda
            append!(ts, t .+ h .* (0:n-1)); append!(ps, y); append!(qs, u)
        end
        # Avança o relógio pela duração deste batimento (n amostras × h s).
        t += n * h
    end
    # Converte o dicionário numa NamedTuple: registro de nomes fixos, lido como
    # `bat.PAM`, `bat.t_min` etc. (cada campo é um vetor, um valor por batimento).
    bat = NamedTuple{(:t_min, :Ce, :BIS, :FC, :VS, :DC, :Rp, :PAS, :PAD, :PAM)}(
        Tuple(col[k] for k in (:t_min, :Ce, :BIS, :FC, :VS, :DC, :Rp, :PAS, :PAD, :PAM)))
    # Devolve a tabela e a onda (ou `nothing`, "nada", se a onda não foi pedida).
    # `condição ? a : b` = "se a condição for verdadeira, a; senão, b".
    return bat, forma_de_onda ? (t_s = ts, pressao = ps, fluxo = qs) : nothing
end

# -----------------------------------------------------------------------------
# Atalho do grupo: farmacologia + circulação numa chamada só
# -----------------------------------------------------------------------------
# Roda a PK/PD do propofol (src/propofol.jl; padrão: Schnider, 30 min, esquema
# do grupo com bolus de 2 mg/kg + manutenção de 100 µg/kg/min) e depois a
# hemodinâmica acima, com o θ basal do grupo (C fixa em 1,22 ml/mmHg).
# Para o paciente virtual calibrado (C do próprio paciente), use `calibrar` e
# `simular_inducao` (src/calibracao.jl). `kwargs...` repassa outras opções
# nomeadas (ex.: forma_de_onda = true) para `simular_hemodinamica`.
"""
    simular_paciente(pac; doses=esquema_inducao(pac), modelo=:schnider, tfim=30.0, kwargs...)

Atalho do grupo: PK/PD (NeoPKPD) + hemodinâmica (Windkessel). Devolve `(pkpd, batimentos, onda)`.
"""
function simular_paciente(pac::Paciente; doses = esquema_inducao(pac), modelo = :schnider, tfim = 30.0, kwargs...)
    pkpd = simular_propofol(pac, doses; modelo, tfim)
    batimentos, onda = simular_hemodinamica(pac, pkpd; tfim, kwargs...)
    return pkpd, batimentos, onda
end
