# =============================================================================
# propofol.jl — PK/PD do propofol sobre o motor do NeoPKPD
# =============================================================================
#
# Segue src/propofol.jl do PropofolPBL do grupo
# (github.com/gabiadamek30-bit/Simula-o-da-resposta-do-paciente-ao-propofol---PBL),
# que usa o NeoPKPD (github.com/shramish2057/NeoPKPD) como motor:
#   - modelo tricompartimental (ThreeCompIVBolus) com bolus e infusões;
#   - compartimento de efeito (ke0) por simulate_biophase;
#   - BIS de Eleveld et al. (2018), calculado com o modelo de Eleveld completo
#     (PK + ke0 + atraso), independente do modelo PK escolhido para a hemodinâmica.
# Tempo em minutos, doses em mg, volumes em L ⇒ concentrações em mg/L = µg/ml.
# =============================================================================
#
# Conceitos usados neste arquivo (para quem não programa):
# - Farmacocinética (PK): o que o corpo faz com o fármaco (distribuição e
#   eliminação). Resultado: concentração plasmática Cp(t) [µg/ml].
# - Farmacodinâmica (PD): o que o fármaco faz com o corpo. Aqui, o BIS; o efeito
#   sobre a pressão (resistência periférica, volume sistólico e FC) fica em
#   acoplamento.jl.
# - Modelo tricompartimental: o corpo é visto como 3 "caixas" ligadas:
#     V1 = compartimento central (plasma e órgãos muito perfundidos) [L];
#     V2 = periférico rápido (costuma-se associar a músculo) [L];
#     V3 = periférico lento (costuma-se associar a gordura) [L].
#   CL = clearance de eliminação (metabolismo, sobretudo hepático) [L/min];
#   Q2, Q3 = clearances intercompartimentais (troca V1<->V2 e V1<->V3) [L/min].
#   A dose entra em V1; Cp = quantidade em V1 / V1.
# - Compartimento de efeito (sítio efetor, o cérebro): não tem volume próprio e
#   segue o plasma com atraso: dCe/dt = ke0·(Cp − Ce). ke0 [1/min] é a
#   velocidade de equilíbrio (meia-vida de equilíbrio = ln 2 / ke0; ex.:
#   ke0 = 0,456/min -> ~1,5 min). O efeito clínico acompanha Ce, não Cp.
# - Unidades: mg / L = µg/ml (a mesma usada por Eleveld 2018).
# Fontes e conferência de cada número: docs/validacao.md (seções 4.1 e 4.2).

# -----------------------------------------------------------------------------
# parametros_pk: calcula os 6 parâmetros PK e o ke0 do paciente, conforme o
# modelo escolhido (:schnider, que é o padrão, :marsh ou :eleveld).
# Entrada: Paciente, nome do modelo (Symbol, escrito com dois-pontos) e, se
# quiser, `eta`. Argumentos depois do `;` são dados pelo nome (ex.:
# parametros_pk(p, :marsh; eta = (CL = 1.2,))).
# Saída: NamedTuple (tupla com nomes; lê-se pk.CL, pk.V1, ...) com CL, V1, Q2,
# V2, Q3, V3 [L/min e L], ke0 [1/min] e o nome do modelo.
# -----------------------------------------------------------------------------
"""
    parametros_pk(pac, modelo=:schnider; eta=NamedTuple())

Parâmetros PK do propofol `(CL, V1, Q2, V2, Q3, V3, ke0)` [L/min, L, 1/min].

- `:schnider` — PK de Schnider et al. 1998 e ke0 = 0,456/min de Schnider et al. 1999.
- `:marsh` — PK de Marsh et al. 1991. O ke0 = 0,26/min não é desse artigo: é o valor
  das bombas Diprifusor, nunca publicado em revista (Absalom et al. 2009).
- `:eleveld` — PK de Eleveld et al. 2018 (ver `parametros_eleveld`).

`eta` multiplica qualquer parâmetro por um fator (variabilidade interindividual,
usada só nas análises experimentais de incerteza e coorte).
"""
function parametros_pk(pac::Paciente, modelo::Symbol = :schnider; eta = NamedTuple())
    if modelo == :schnider
        # Schnider et al. 1998 (PK) e Schnider et al. 1999 (ke0): equações
        # lineares centradas no "paciente médio" do estudo (53 anos, 77 kg,
        # massa magra 59 kg, 177 cm). V1 e V3 são fixos; V2 e Q2 diminuem com a
        # idade; CL depende de peso, massa magra (James, paciente.jl) e altura.
        # Conferido com o PropofolPBL do grupo e com o Python Anesthesia
        # Simulator (docs/validacao.md, seção 4.1).
        lbm = massa_magra(pac)
        V1 = 4.27
        V2 = 18.9 - 0.391 * (pac.idade - 53)
        V3 = 238.0
        CL = 1.89 + 0.0456 * (pac.peso - 77) - 0.0681 * (lbm - 59) + 0.0264 * (pac.altura - 177)
        Q2 = 1.29 - 0.024 * (pac.idade - 53)
        Q3 = 0.836
        ke0 = 0.456
        # as equações lineares de Schnider ficam negativas fora da faixa de adultos
        # (ex.: V2 < 0 acima de 101 anos); nesse caso o modelo não se aplica
        min(V2, CL, Q2) > 0 ||
            throw(ArgumentError("paciente fora da faixa de validade do modelo de Schnider (V2, CL ou Q2 ≤ 0)"))
    elseif modelo == :marsh
        # Marsh et al. 1991: tudo proporcional ao peso. O artigo dá V1 =
        # 0,228 L/kg e microconstantes [1/min]: k10 = 0,119, k12 = 0,112,
        # k21 = 0,055, k13 = 0,0419, k31 = 0,0033. Como k10 = CL/V1,
        # k12 = Q2/V1, k21 = Q2/V2, k13 = Q3/V1 e k31 = Q3/V3, daí saem CL, Q2,
        # V2, Q3 e V3 abaixo. O ke0 = 0,26/min é o das bombas Diprifusor, não do
        # artigo de Marsh (Absalom et al. 2009; docs/validacao.md, seção 4.1).
        V1 = 0.228 * pac.peso
        CL = 0.119 * V1
        Q2 = 0.112 * V1
        V2 = Q2 / 0.055
        Q3 = 0.0419 * V1
        V3 = Q3 / 0.0033
        ke0 = 0.26
    elseif modelo == :eleveld
        # Eleveld et al. 2018: calculado em parametros_eleveld (mais abaixo).
        # `e.CL` lê o campo CL; a linha seguinte copia os 7 valores de uma vez.
        e = parametros_eleveld(pac)
        CL, V1, Q2, V2, Q3, V3, ke0 = e.CL, e.V1, e.Q2, e.V2, e.Q3, e.V3, e.ke0
    else
        # Nome de modelo inválido: para com erro. `$(modelo)` põe o valor no texto.
        throw(ArgumentError("modelo PK desconhecido: $(modelo)"))
    end
    # Variabilidade interindividual (opcional): `eta` é uma NamedTuple como
    # (CL = 1.2, ke0 = 0.8), que multiplica CL por 1,2 e ke0 por 0,8.
    # `get(eta, s, 1.0)` devolve o fator do parâmetro `s`, ou 1,0 (sem mudança)
    # se ele não foi dado. Só as análises experimentais de incerteza e coorte
    # (incerteza.jl, populacao.jl) usam isso; no uso normal `eta` é vazio.
    # Saída: NamedTuple com os parâmetros finais e o nome do modelo.
    f(s) = Float64(get(eta, s, 1.0))
    return (CL = CL * f(:CL), V1 = V1 * f(:V1), Q2 = Q2 * f(:Q2), V2 = V2 * f(:V2),
            Q3 = Q3 * f(:Q3), V3 = V3 * f(:V3), ke0 = ke0 * f(:ke0), modelo = modelo)
end

# -----------------------------------------------------------------------------
# Eventos de dose. `DoseEvent(tempo, quantidade, duração)` é o tipo do NeoPKPD:
# tempo de início [min], quantidade total [mg] e duração [min]. Com duração > 0
# o NeoPKPD dá a quantidade a taxa constante nesse intervalo
# (taxa = quantidade / duração). `Float64(...)` só converte para número real.
# bolus: a dose inteira [mg] injetada em 0,5 min = 30 s (padrão), como na
# prática, e não de forma instantânea.
# -----------------------------------------------------------------------------
"Bolus de `dose_mg` em t (min), administrado em `duracao` min (padrão 30 s)."
bolus(t, dose_mg; duracao = 0.5) = DoseEvent(Float64(t), Float64(dose_mg), Float64(duracao))

# infusao: infusão contínua de `taxa_mg_min` [mg/min] durante `duracao` [min],
# convertida em um DoseEvent com quantidade total = taxa × duração [mg].
"Infusão contínua de `taxa_mg_min` de t_inicio a t_inicio + duracao (min)."
infusao(t_inicio, duracao, taxa_mg_min) =
    DoseEvent(Float64(t_inicio), Float64(taxa_mg_min * duracao), Float64(duracao))

# -----------------------------------------------------------------------------
# esquema_inducao: monta a lista de doses do esquema típico do grupo
# (PropofolPBL): bolus de 2,0 mg/kg em 30 s no tempo 0 e manutenção de
# 100 µg/kg/min de 1 a 31 min. Todos esses valores podem ser trocados pelos
# argumentos por nome. Saída: vetor de DoseEvent para simular_propofol.
# -----------------------------------------------------------------------------
"""
    esquema_inducao(pac; bolus_mgkg=2.0, duracao_bolus=0.5, manutencao_mcgkgmin=100.0,
                    inicio_manutencao=1.0, duracao_manutencao=30.0, t_bolus=0.0)

Esquema típico do grupo: bolus de indução (mg/kg) em `t_bolus` seguido de
manutenção (µg/kg/min) a partir de `t_bolus + inicio_manutencao`. Tempos em min.
"""
function esquema_inducao(pac::Paciente; bolus_mgkg = 2.0, duracao_bolus = 0.5, manutencao_mcgkgmin = 100.0,
                         inicio_manutencao = 1.0, duracao_manutencao = 30.0, t_bolus = 0.0)
    # `DoseEvent[]` cria uma lista vazia de doses; `push!` acrescenta uma dose.
    # `condição && ação` executa a ação só se a condição for verdadeira: o bolus
    # só entra se a dose for > 0. mg/kg × kg = mg.
    doses = DoseEvent[]
    bolus_mgkg > 0 && push!(doses, bolus(t_bolus, bolus_mgkg * pac.peso; duracao = duracao_bolus))
    # Manutenção: µg/kg/min × kg / 1000 = mg/min (1 mg = 1000 µg).
    if manutencao_mcgkgmin > 0 && duracao_manutencao > 0
        push!(doses, infusao(t_bolus + inicio_manutencao, duracao_manutencao, manutencao_mcgkgmin * pac.peso / 1000))
    end
    return doses
end

# -----------------------------------------------------------------------------
# parametros_eleveld: PK do propofol de Eleveld et al. 2018 (modelo "universal",
# estimado com dados de recém-nascidos a idosos, incluindo obesos) e ke0 do BIS.
# Entrada: Paciente; `arterial` (amostras de sangue arterial, padrão) e
# `opioides` (presença de opioide, padrão não). Saída: NamedTuple com CL, V1,
# Q2, V2, Q3, V3 [L/min, L] e ke0 [1/min].
# Ideia: cada parâmetro = valor do paciente de referência (homem, 35 anos,
# 70 kg, 170 cm) × fatores de correção por peso, idade, maturação e sexo.
# Conferido com o Python Anesthesia Simulator (PAS) em 12 casos, erro relativo
# < 1e−12 (teste automático em test/runtests.jl; docs/validacao.md, seção 4.1).
# -----------------------------------------------------------------------------
"""
    parametros_eleveld(pac; arterial=true, opioides=false)

PK do propofol de Eleveld et al. 2018 (parâmetros PK) e ke0 do BIS (Tabela 3), na forma
do Python Anesthesia Simulator (pk_models.py, linhas 256–335, commit 4afc03a).
Massa livre de gordura de Al-Sallami; idade pós-menstrual = idade + 40 semanas.
`arterial = true` usa os parâmetros de amostras arteriais (ke0 = 0,146/min a 70 kg);
`opioides = false` porque o contexto de uso é a indução sem opioide.
"""
function parametros_eleveld(pac::Paciente; arterial = true, opioides = false)
    # θ = constantes estimadas por Eleveld 2018 (parâmetros PK), na mesma ordem do
    # artigo e do PAS. θ[1] é o primeiro (Julia conta a partir de 1):
    #  θ1  V1 de referência = 6,28 L
    #  θ2  V2 de referência = 25,5 L
    #  θ3  V3 de referência = 272,8 L
    #  θ4  CL de referência, homem = 1,79 L/min
    #  θ5  Q2 de referência = 1,75 L/min
    #  θ6  Q3 de referência = 1,11 L/min
    #  θ7  erro residual típico (0,191); não é usado aqui (serve ao ajuste
    #      estatístico do artigo)
    #  θ8  maturação do CL: idade pós-menstrual com 50% do CL adulto = 42,3 sem.
    #  θ9  maturação do CL: inclinação da sigmoide (expoente) = 9,05
    #  θ10 V2 menor com a idade: −0,0156 por ano em relação aos 35 anos
    #  θ11 CL menor com a idade: −0,00286 por ano (só aplicado com opioides)
    #  θ12 peso com 50% do V1 máximo = 33,6 kg
    #  θ13 V3 menor com a idade: −0,0138 por ano (só aplicado com opioides)
    #  θ14 maturação do Q3: idade pós-menstrual com 50% do Q3 adulto = 68,3 sem.
    #  θ15 CL de referência, mulher = 2,10 L/min
    #  θ16 Q2 maior enquanto o Q3 ainda não maturou (1,30; só pesa em crianças)
    #  θ17 V1 maior em amostras venosas (crianças) (1,42)
    #  θ18 Q2 em amostras venosas (fator 0,68)
    # Descrições tiradas dos comentários do PAS (pk_models.py, linhas 269–287,
    # commit 4afc03a); valores conferidos em docs/validacao.md, seção 4.1.
    θ = (6.2830780766822, 25.5013145036879, 272.8166615043603, 1.7895836588902, 1.7500983738779,
         1.1085424008536, 0.191307, 42.2760190602615, 9.0548452392807, -0.015633, -0.00285709,
         33.5531248778544, -0.0138166, 68.2767978846832, 2.1002218877899, 1.3042680471360,
         1.4189043652084, 0.6805003109141)
    # Covariáveis do paciente. `a, b, c = x, y, z` atribui várias de uma vez.
    # `homem` vale true (verdadeiro) se o sexo for :M.
    idade, peso, homem = pac.idade, pac.peso, pac.sexo == :M
    # Funções auxiliares de uma linha (só existem dentro desta função):
    # sig: sigmoide de Hill, x^γ / (c50^γ + x^γ); vale 0,5 quando x = c50.
    # fcentral: quanto o V1 "satura" com o peso (50% do máximo em θ12 = 33,6 kg).
    # fenv: efeito do envelhecimento, exp(x·(idade − 35)); vale 1 aos 35 anos.
    sig(x, c50, γ) = x^γ / (c50^γ + x^γ)
    fcentral(w) = sig(w, θ[12], 1)
    fenv(x) = exp(x * (idade - 35))
    # Massa livre de gordura (FFM) [kg] de Al-Sallami et al. 2015: a fórmula de
    # adulto 9270·peso/(6680 + 216·IMC) (homem) ou 9270·peso/(8780 + 244·IMC)
    # (mulher), multiplicada por um fator de idade que corrige crianças e
    # adolescentes (≈ 1 no adulto). h = homem?, w = peso, a = idade.
    # Citação completa ainda não conferida (docs/validacao.md, referência [7]).
    al_sallami(h, w, a, imc) = h ? (0.88 + 0.12 / (1 + (a / 13.4)^-12.7)) * 9270 * w / (6680 + 216 * imc) :
                                   (1.11 - 0.11 / (1 + (a / 7.1)^-1.1)) * 9270 * w / (8780 + 244 * imc)
    # Idade pós-menstrual (PMA) em semanas = idade + 40 semanas de gestação
    # (supõe nascimento a termo). pma_ref = a do paciente de referência (35 anos).
    pma, pma_ref = (idade + 40 / 52) * 52, (35 + 40 / 52) * 52          # semanas
    # Fatores de maturação (sigmoides da PMA): o CL chega a 50% em θ8 semanas e o
    # Q3 em θ14 semanas. Cada um é dividido pelo valor do adulto de referência,
    # então no adulto a razão é ≈ 1; só pesam em crianças.
    fCLmat, fCLmat_ref = sig(pma, θ[8], θ[9]), sig(pma_ref, θ[8], θ[9])
    fQ3mat, fQ3mat_ref = sig(pma, θ[14], 1), sig(pma_ref, θ[14], 1)
    # FFM do paciente e do de referência (homem, 70 kg, 35 anos, 1,70 m).
    ffm, ffm_ref = al_sallami(homem, peso, idade, imc(pac)), al_sallami(true, 70.0, 35.0, 70 / 1.7^2)
    # Efeito dos opioides: exp(x·idade) se houver opioide; senão 1 (sem efeito).
    fop(x) = opioides ? exp(x * idade) : 1.0
    # V1 [L]: cresce com o peso, mas satura (sigmoide com 50% em θ12 = 33,6 kg).
    # Em amostras venosas, V1 aumenta em crianças (θ17). A forma
    # `arterial || (...)` executa o parêntese só quando `arterial` é falso.
    # `x *= y` é o mesmo que `x = x * y`.
    V1 = θ[1] * fcentral(peso) / fcentral(70.0)
    arterial || (V1 *= 1 + θ[17] * (1 - fcentral(peso)))
    # V2 [L]: proporcional ao peso e menor com a idade (θ10 < 0).
    V2 = θ[2] * peso / 70 * fenv(θ[10])
    # V3 [L]: proporcional à massa livre de gordura; menor com a idade só com
    # opioides (θ13).
    V3 = θ[3] * ffm / ffm_ref * fop(θ[13])
    # CL [L/min]: valor de referência por sexo (θ4 homem, θ15 mulher), escalado
    # de forma alométrica por (peso/70)^0,75, corrigido pela maturação e, com
    # opioides, pela idade (θ11).
    CL = (homem ? θ[4] : θ[15]) * (peso / 70)^0.75 * fCLmat / fCLmat_ref * fop(θ[11])
    # Q2 [L/min]: escala alométrica com V2 (expoente 0,75); maior enquanto o Q3
    # não maturou (θ16, só crianças); com amostras venosas, × θ18.
    Q2 = θ[5] * (V2 / θ[2])^0.75 * (1 + θ[16] * (1 - fQ3mat))
    arterial || (Q2 *= θ[18])
    # Q3 [L/min]: escala alométrica com V3 e corrigido pela maturação do Q3.
    Q3 = θ[6] * (V3 / θ[3])^0.75 * fQ3mat / fQ3mat_ref
    # ke0 [1/min] do BIS (Eleveld 2018, Tabela 3): 0,146 (arterial) ou 1,24
    # (venoso) para 70 kg, escalado por (peso/70)^−0,25 (quem pesa menos
    # equilibra mais rápido). É o ke0 usado sempre para o BIS.
    ke0 = (arterial ? 0.146 : 1.24) * (peso / 70)^-0.25
    return (CL = CL, V1 = V1, Q2 = Q2, V2 = V2, Q3 = Q3, V3 = V3, ke0 = ke0)
end

# -----------------------------------------------------------------------------
# bis_eleveld: transforma a concentração no sítio efetor Ce [µg/ml] em BIS
# (índice biespectral do EEG, 0–100; acordado ≈ 93; alvo na anestesia 40–60).
# Entrada: Ce [µg/ml] e idade [anos]. Saída: BIS (de 0 a 93).
# Modelo: sigmoide de Hill (Emax), BIS = E0 − E0·Ce^γ / (Ce50^γ + Ce^γ).
# Quando Ce = Ce50, o BIS cai à metade (46,5). Valores de Eleveld 2018,
# Tabela 3; conferidos em docs/validacao.md, seção 4.2.
# -----------------------------------------------------------------------------
"""
    bis_eleveld(Ce, idade)

BIS em função da concentração no sítio efetor do modelo de Eleveld et al. 2018
(Tabela 3): E0 = 93, Ce50 = 3,08·exp(−0,00635·(idade − 35)) µg/ml,
γ = 1,89 para Ce ≤ Ce50 e 1,47 acima.
"""
function bis_eleveld(Ce, idade)
    # E0 = BIS basal (acordado) = 93; o efeito máximo também é 93 (BIS -> 0).
    E0 = 93.0
    # Ce50 [µg/ml]: 3,08 aos 35 anos e cai ~0,6% por ano de idade, ou seja, o
    # idoso é mais sensível (ex.: 70 anos -> Ce50 ≈ 2,47 µg/ml).
    Ce50 = 3.08 * exp(-0.00635 * (idade - 35))
    # γ (inclinação da curva): 1,89 até o Ce50 e 1,47 acima (Eleveld 2018). No
    # código do grupo estava invertido; corrigido (docs/validacao.md, seção 3).
    γ = Ce <= Ce50 ? 1.89 : 1.47
    # Ce negativo (ruído numérico) vira 0, para não elevar número negativo a
    # potência fracionária.
    c = max(Float64(Ce), 0.0)
    # BIS = E0 × (1 − fração do efeito máximo).
    return E0 - E0 * c^γ / (Ce50^γ + c^γ)
end

# atraso_bis: o monitor mostra o BIS com atraso (processamento do EEG).
# Eleveld 2018: atraso [s] = 15 + e^(0,0517·idade); ex.: 35 anos -> ~21 s,
# 70 anos -> ~52 s. Sai em segundos (convertido para minutos em
# simular_propofol).
"Atraso do BIS [s] de Eleveld et al. 2018 (Tabela 3, θ = 0,0517): 15 + e^(0,0517·idade)."
atraso_bis(idade) = 15 + exp(0.0517 * idade)

# -----------------------------------------------------------------------------
# _cp_ce (o "_" no início do nome indica função de uso interno do pacote).
# Entrada: parâmetros PK `pk` (NamedTuple de parametros_pk), lista de doses e a
# grade de tempos [min]. Saída: (t, Cp, Ce), concentrações em µg/ml.
# Passos: (1) monta o modelo tricompartimental do NeoPKPD e resolve as equações
# diferenciais -> Cp; (2) a partir de Cp, calcula Ce pelo compartimento de
# efeito.
# -----------------------------------------------------------------------------
# Cp e Ce na grade (min) para os parâmetros `pk`, pelo NeoPKPD
function _cp_ce(pk, doses, grade)
    # ModelSpec: tipo de modelo (3 compartimentos, entrada intravenosa), nome do
    # fármaco, parâmetros e doses ordenadas pelo tempo. `d -> d.time` é uma
    # função anônima ("dado d, devolve d.time"), usada como critério de ordem.
    spec = ModelSpec(ThreeCompIVBolus(), "propofol",
                     ThreeCompIVBolusParams(pk.CL, pk.V1, pk.Q2, pk.V2, pk.Q3, pk.V3),
                     sort(doses; by = d -> d.time))
    # Integrador numérico: Tsit5 (método de Runge-Kutta de 5ª ordem), tolerâncias
    # relativa 1e−8 e absoluta 1e−10 (alta precisão), no máximo 10^7 passos.
    solver = SolverSpec(:Tsit5, 1e-8, 1e-10, 10^7)
    # Resolve o modelo de grade[1] a grade[end] e guarda os valores na grade.
    res = NeoPKPD.simulate(spec, SimGrid(grade[1], grade[end], grade), solver)
    # `res.observations[:conc]` = Cp = quantidade em V1 / V1 [mg/L = µg/ml].
    # O ponto em `max.(...)` aplica `max` a cada elemento do vetor
    # ("broadcasting").
    # max(…, 0) só remove ruído numérico do integrador (valores da ordem de −1e−12)
    Cp = max.(res.observations[:conc], 0.0)
    # Compartimento de efeito: resolve dCe/dt = ke0·(Cp − Ce) sobre a curva Cp.
    # Os números (0.0, 1.0, 1.0) são E0, Emax e EC50 de um efeito genérico do
    # NeoPKPD que não usamos; esse efeito é descartado em `_`.
    Ce, _ = simulate_biophase(res.t, Cp, BiophaseEquilibrationParams(pk.ke0, 0.0, 1.0, 1.0), solver)
    # Devolve tempos [min], Cp e Ce [µg/ml] (Ce também sem ruído negativo).
    return res.t, Cp, max.(Ce, 0.0)
end

# -----------------------------------------------------------------------------
# simular_propofol: função principal deste arquivo. Entrada: Paciente, lista de
# doses (DoseEvent), modelo PK, tempo final `tfim` [min] e passo `dt` [min]
# (padrão 1/60 min = 1 s). `pk` pode vir pronto (ex.: com variabilidade `eta`).
# Saída: NamedTuple (t, Cp, Ce, bis, Ce_bis, pk).
# Por que dois cálculos de PK? A hemodinâmica (acoplamento.jl) usa o Ce do
# modelo escolhido (Schnider por padrão). Já o BIS usa SEMPRE a PK, o ke0 e o
# atraso de Eleveld, porque o Ce50 = 3,08 µg/ml do BIS foi estimado junto com
# eles e só vale nesse conjunto (docs/validacao.md, seção 3, item 2).
# -----------------------------------------------------------------------------
"""
    simular_propofol(pac, doses; modelo=:schnider, tfim=60.0, dt=1/60, pk=parametros_pk(pac, modelo))

Simula o propofol com o NeoPKPD. Retorna `(t, Cp, Ce, bis, Ce_bis, pk)`, t em minutos
e concentrações em µg/ml:

- `Cp`, `Ce`: plasma e sítio efetor do modelo PK escolhido (`pk`), que movem a hemodinâmica;
- `bis`: BIS do modelo de Eleveld completo — PK e ke0 de Eleveld (`Ce_bis`), sigmoide
  de Eleveld e atraso de 15 + e^(0,0517·idade) s —, porque o Ce50 de 3,08 µg/ml só vale
  com o PK e o ke0 com que foi estimado.
"""
function simular_propofol(pac::Paciente, doses::Vector{DoseEvent}; modelo = :schnider,
                          tfim = 60.0, dt = 1 / 60, pk = parametros_pk(pac, modelo))
    # Grade de tempos 0, dt, 2dt, ..., tfim [min]; `collect` transforma o
    # intervalo `início:passo:fim` em vetor.
    grade = collect(0.0:dt:tfim)
    # Caso sem nenhuma dose: concentrações zero e BIS constante no basal (93).
    # `zeros(n)` = vetor de n zeros; `fill(v, n)` = vetor com n cópias de v;
    # `copy` evita que Ce e Ce_bis sejam o mesmo vetor, na memória, que Cp.
    if isempty(doses)          # sem fármaco: tudo zero (NeoPKPD exige ao menos uma dose)
        z = zeros(length(grade))
        return (t = grade, Cp = z, Ce = copy(z), bis = fill(bis_eleveld(0.0, pac.idade), length(grade)),
                Ce_bis = copy(z), pk = pk)
    end
    # Cp e Ce do modelo PK escolhido: são eles que movem a hemodinâmica.
    t, Cp, Ce = _cp_ce(pk, doses, grade)
    # Ce do modelo de Eleveld completo (PK + ke0 de Eleveld), só para o BIS.
    # `_` descarta t e Cp, que não são usados aqui.
    _, _, Ce_bis = _cp_ce(parametros_pk(pac, :eleveld), doses, grade)
    # BIS sem atraso em cada instante. `bis_eleveld.(Ce_bis, idade)`: o ponto
    # aplica a função a cada elemento do vetor Ce_bis.
    bis0 = bis_eleveld.(Ce_bis, pac.idade)
    # Atraso do BIS convertido de segundos para minutos.
    d = atraso_bis(pac.idade) / 60                  # min
    # BIS mostrado no tempo t = BIS calculado no tempo t − atraso, por
    # interpolação linear (função `interp`, em acoplamento.jl). Antes de passar
    # o atraso, repete o valor inicial (basal). `[... for ti in t]` é uma
    # "compreensão de lista": calcula um valor para cada tempo ti.
    bis = [interp(t, bis0, max(ti - d, t[1])) for ti in t]
    # Resultado: tempos [min], Cp, Ce, BIS, Ce de Eleveld [µg/ml] e a PK usada.
    return (t = t, Cp = Cp, Ce = Ce, bis = bis, Ce_bis = Ce_bis, pk = pk)
end

# -----------------------------------------------------------------------------
# matriz_pk: escreve o mesmo modelo tricompartimental como sistema linear
# dA/dt = K·A, em que A = quantidades [mg] em V1, V2 e V3. Serve só para o teste
# de verificação (test/runtests.jl): depois de uma dose em V1, a solução exata
# é A(t) = exp(K·t)·A(0), comparada com o resultado do NeoPKPD.
# Microconstantes [1/min]: k10 = CL/V1 (eliminação), k12 = Q2/V1, k21 = Q2/V2,
# k13 = Q3/V1, k31 = Q3/V3 (transferências entre compartimentos).
# -----------------------------------------------------------------------------
"Matriz K do sistema dA/dt = K·A (1/min), para o teste de verificação contra exp(K·t)."
function matriz_pk(pk)
    k10, k12, k21, k13, k31 = pk.CL / pk.V1, pk.Q2 / pk.V1, pk.Q2 / pk.V2, pk.Q3 / pk.V1, pk.Q3 / pk.V3
    # Matriz 3×3 (`;` separa linhas, espaço separa colunas). Linha 1: V1 perde
    # para eliminação, V2 e V3 e recebe de V2 e V3. Linhas 2 e 3: V2 e V3
    # recebem de V1 e devolvem a V1.
    return [-(k10 + k12 + k13) k21 k31; k12 -k21 0.0; k13 0.0 -k31]
end
