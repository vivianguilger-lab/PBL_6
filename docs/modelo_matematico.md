# Modelo matemático

> **Pergunta:** *como prever, por simulação, a resposta hemodinâmica de um paciente ao propofol?*

O núcleo do modelo segue o **PropofolPBL** do grupo
([gabiadamek30-bit/Simula-o-da-resposta-do-paciente-ao-propofol---PBL](https://github.com/gabiadamek30-bit/Simula-o-da-resposta-do-paciente-ao-propofol---PBL)),
que junta dois repositórios:

| origem | uso |
|---|---|
| [shramish2057/NeoPKPD](https://github.com/shramish2057/NeoPKPD) (pacote Julia registrado) | motor PK/PD: modelo tricompartimental com bolus e infusões, compartimento de efeito (`ke0`) |
| [hpigot/windkessel-id](https://github.com/hpigot/windkessel-id) (MIT) | Windkessel de 2/3/4 elementos, identificação a partir de pressão/fluxo e dados aórticos |

Equações, valores padrão e testes do grupo foram mantidos. Sobre essa base, este pacote
acrescenta: o paciente virtual calibrado também com PAS/PAD, as métricas de hipotensão, a
janela de dose, a faixa de incerteza, a coorte virtual e a plataforma web.

## 1. Fluxo

```
Paciente ──► PK (Schnider/Marsh, NeoPKPD) ──► Cp(t) ──► ke0 ──► Ce(t) ──► BIS (Eleveld)
                                                                  │
  PAM, FC, DC basais ──► θ = [Rp, C, Rc, L] ◄─────────────────────┘ (Emax sigmoide × sensibilidade)
  PAS − PAD ──► C (paciente virtual)  │
                                       ▼
          fluxo aórtico (Stergiopulos) ─► Windkessel de 4 elementos ─► pressão arterial batimento a batimento
```

## 2. Biomarcadores

| biomarcador | uso no modelo |
|---|---|
| idade, peso, altura, sexo | PK de Schnider (massa magra de James); idade também no Ce50 do BIS |
| PAM, FC (e DC, se medido) | Windkessel basal: Rp = PAM/DC |
| PAS, PAD | complacência arterial C do paciente virtual (extensão deste pacote) |
| ASA, PVI, IVCI, idade | fator de sensibilidade hemodinâmica ao propofol |
| SpO₂, glicemia | registrados; ainda **não** entram em nenhuma equação |

DC basal, quando não medido: 3,0 L/min/m² × superfície corporal (Mosteller).

## 3. Farmacocinética (NeoPKPD)

Modelo tricompartimental (`ThreeCompIVBolus`) com bolus (padrão: 30 s) e infusão de manutenção:

| | Schnider (1998/1999) | Marsh (1991) |
|---|---|---|
| V1 | 4,27 L | 0,228·peso |
| V2 | 18,9 − 0,391·(idade − 53) L | Q2/0,055 |
| V3 | 238 L | Q3/0,0033 |
| CL | 1,89 + 0,0456·(peso − 77) − 0,0681·(LBM − 59) + 0,0264·(altura − 177) L/min | 0,119·V1 |
| Q2 | 1,29 − 0,024·(idade − 53) L/min | 0,112·V1 |
| Q3 | 0,836 L/min | 0,0419·V1 |
| ke0 | 0,456 /min | 0,26 /min |

Sítio efetor (`simulate_biophase` do NeoPKPD): dCe/dt = ke0·(Cp − Ce).

**BIS (Eleveld 2018):** BIS = 93 − 93·Ceᵞ/(Ce50ᵞ + Ceᵞ), Ce50 = 3,08·e^{−0,00635(idade−35)} µg/ml,
γ = 1,47 abaixo de Ce50 e 1,89 acima.

## 4. Windkessel de 4 elementos (windkessel-id)

Forma paralela (Burattini, Westerhof, Stergiopulos), estados x₁ = volume em C e x₂ (ramo indutivo):

```
ẋ = diag(−1/(C·Rp), −Rc/L)·x + [1, Rc]·Q(t)
P = x₁/C − (Rc/L)·x₂ + Rc·Q(t)
```

Discretização ZOH analítica (matriz A diagonal), passo h = 5 ms. Propriedades usadas na
verificação: |Z(0)| = Rp ⇒ **PAM = Rp·Q̄**; |Z(∞)| = Rc. O Windkessel de 3 elementos é o
caso L → ∞ e o de 2 elementos, L = 0 (com τ = Rp·C, como no Resumo da disciplina).

Parâmetros basais (Stergiopulos 1999, humano): C = 1,22 ml/mmHg, Rc = 0,056 mmHg·s/ml,
L = 0,0051 mmHg·s²/ml; **Rp = PAM / DC** (DC em ml/s).

**Entrada:** o batimento de fluxo aórtico de Stergiopulos (1999), normalizado e esticado para
a duração 60/FC, escalado para ejetar o volume sistólico VS.

## 5. Acoplamento PK/PD → hemodinâmica (grupo)

A cada batimento, com o Ce daquele instante:

```
X = X₀ · max(1 − Imax_X · s · Ceᵞ/(Ce50ᵞ + Ceᵞ), 0,05)      X ∈ {Rp, VS, FC}
```

| parâmetro | valor | significado |
|---|---|---|
| Ce50 | 3,0 µg/ml | concentração de meio-efeito |
| γ | 2,0 | inclinação de Hill |
| Imax_Rp | 0,25 | vasodilatação arterial |
| Imax_VS | 0,15 | depressão miocárdica / redução de pré-carga |
| Imax_FC | 0,05 | redução da FC |

Sensibilidade do paciente:

```
s = 1 + 0,01·max(idade − 50, 0) + 0,15·max(ASA − 2, 0) + 0,03·max(PVI − 13, 0) + 0,01·max(IVCI − 40, 0)
```

> Como diz o README do grupo: esses valores são **ilustrativos** (ordem de grandeza da queda
> de PAM na indução), não ajustados a dados, e precisam ser calibrados.

## 6. Paciente virtual (extensão)

1. Rp = PAM/DC (como no grupo).
2. **C identificado** por bisseção em log C ∈ [0,2; 6] ml/mmHg para que a pressão de pulso do
   batimento periódico seja PAS − PAD (a PP cai quando C aumenta). É a mesma ideia de
   identificação do windkessel-id, com um único parâmetro e uma única medida.
3. Rc e L ficam nos valores de Stergiopulos. Com a linha arterial do paciente (pressão contínua
   + estimativa de fluxo), θ inteiro pode ser identificado com `identify`
   (`scripts/06_identificacao_windkessel.jl`).

## 7. Incerteza e coorte virtual (extensão)

Variabilidade interindividual log-normal: PK com as variâncias de Schnider (V1 0,0404; V2 0,01;
V3 0,1435; CL 0,1005; Q2 0,01; Q3 0,1179) e ω = 0,3 para ke0; resposta hemodinâmica com
ω = 0,3 em Ce50 e nos três Imax (assumidos). A previsão individual mantém a hemodinâmica
basal calibrada e sorteia a farmacologia (200 sorteios) → mediana, faixa de 90% e
P(PAM < 65), P(PAM < 55).

## 8. Métricas

PAM basal (1 min antes do bolus), PAM mínima e instante, queda (%), tempo com PAM < 65 e
< 55 mmHg e área abaixo do limiar nos 20 min após o bolus, BIS mínimo e tempo até BIS < 60.
Risco: **ALTO** (qualquer PAM < 55), **MODERADO** (PAM < 65 por ≥ 1 min ou queda > 30%),
**BAIXO** (demais).

## 9. Classificação (Resumo, Seções 2, 3 e 7)

| dicotomia | este modelo |
|---|---|
| tempo contínuo × discreto | PK contínua (EDOs no NeoPKPD); Windkessel discretizado (ZOH) batimento a batimento |
| linear × não linear | Windkessel linear; acoplamento não linear (Hill) |
| determinístico × estocástico | determinístico por paciente; estocástico na faixa de incerteza e na coorte |
| tempo-invariante × tempo-variante | tempo-variante (Rp, VS e FC mudam com Ce) |
| concentrado × distribuído | concentrado (0-D) |
| black box × grey box | grey box: estrutura física, parâmetros do paciente e da literatura |
| paramétrico × não paramétrico | paramétrico |

## 10. Diferenças em relação ao modelo anterior deste repositório

O modelo anterior (até o commit `d3dd331`) usava um coração de elastância variável em malha
fechada com barorreflexo e dois compartimentos de efeito. Ele foi substituído pelo modelo do
grupo, que é aberto: o fluxo aórtico é imposto (Stergiopulos) e o propofol age diretamente
sobre Rp, VS e FC. Consequências: não há barorreflexo nem retorno venoso explícitos (a queda
de VS representa tanto depressão miocárdica quanto redução de pré-carga), e o efeito
hemodinâmico usa o mesmo Ce do BIS (ke0 único).
