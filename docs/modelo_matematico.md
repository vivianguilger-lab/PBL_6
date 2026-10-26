# Modelo matemático

> **Pergunta:** *como prever, por simulação, a resposta hemodinâmica de um paciente ao propofol?*

O núcleo do modelo segue o **PropofolPBL** do grupo
([gabiadamek30-bit/Simula-o-da-resposta-do-paciente-ao-propofol---PBL](https://github.com/gabiadamek30-bit/Simula-o-da-resposta-do-paciente-ao-propofol---PBL)),
que junta dois repositórios:

| origem | uso |
|---|---|
| [shramish2057/NeoPKPD](https://github.com/shramish2057/NeoPKPD) (pacote Julia registrado) | motor PK/PD: modelo tricompartimental com bolus e infusões, compartimento de efeito (`ke0`) |
| [hpigot/windkessel-id](https://github.com/hpigot/windkessel-id) (MIT) | Windkessel de 2/3/4 elementos, identificação a partir de pressão/fluxo e dados aórticos |

Equações, valores padrão e testes do grupo foram mantidos, com as correções da validação
([validacao.md](validacao.md), seção 3). Sobre essa base, este pacote acrescenta: o BIS pelo
modelo de Eleveld completo, o paciente virtual calibrado também com PAS/PAD, as métricas de
hipotensão, a janela de dose e a plataforma web. Cada valor citado aqui tem arquivo, linha e
fonte em [validacao.md](validacao.md).

## 1. Fluxo

```
Paciente ──► PK de Eleveld ──► Cp ──► ke0 de Eleveld ──► Ce ──► atraso ──► BIS (Eleveld 2018)
        └──► PK escolhido (Schnider/Marsh/Eleveld, NeoPKPD) ──► Cp(t) ──► ke0 ──► Ce(t)
                                                                               │
  PAM, FC, DC basais ──► θ = [Rp, C, Rc, L] ◄──────────────────────────────────┘ (Emax sigmoide × sensibilidade)
  PAS − PAD ──► C (paciente virtual)  │
                                       ▼
          fluxo aórtico (Stergiopulos) ─► Windkessel de 4 elementos ─► pressão arterial batimento a batimento
```

## 2. Biomarcadores

| biomarcador | uso no modelo |
|---|---|
| idade, peso, altura, sexo | PK (Schnider: massa magra de James; Eleveld: massa livre de gordura de Al-Sallami); BIS de Eleveld |
| PAM (do monitor), FC (e DC, se medido) | Windkessel basal: Rp = PAM/DC |
| PAS, PAD | complacência arterial C do paciente virtual (extensão deste pacote) |
| ASA, PVI, IVCI, idade | fator de sensibilidade hemodinâmica ao propofol |

SpO₂ e glicemia existem na estrutura `Paciente` (como no grupo), mas não entram no modelo e
saíram da plataforma.

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
| ke0 | 0,456 /min (Schnider 1999) | 0,26 /min (Diprifusor; não é de Marsh 1991) |

Se as equações lineares de Schnider derem V2, CL ou Q2 ≤ 0 (fora da faixa de adultos), o
programa para com erro; não há pisos.

**Eleveld (2018):** modelo alométrico com idade, peso, altura, sexo e massa livre de gordura
de Al-Sallami (parâmetros e equações em `parametros_eleveld`, iguais aos do Python Anesthesia
Simulator), amostras arteriais, sem opioide; ke0 = 0,146·(peso/70)^−0,25 /min.

Sítio efetor (`simulate_biophase` do NeoPKPD): dCe/dt = ke0·(Cp − Ce).

**BIS (Eleveld 2018, modelo completo):** calculado sempre com o PK e o ke0 de Eleveld, qualquer
que seja o PK escolhido para a hemodinâmica, porque o Ce50 foi estimado com eles:

```
BIS(t) = 93 − 93·Ce(t − d)ᵞ / (Ce50ᵞ + Ce(t − d)ᵞ)
Ce50 = 3,08·e^(−0,00635·(idade − 35)) µg/ml
γ = 1,89 se Ce ≤ Ce50;  1,47 se Ce > Ce50
d = 15 + e^(0,0517·idade) s   (atraso do BIS)
```

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
s = 1 + 0,01·max(idade − 50, 0) + 0,15·max(ASA − 2, 0) + 0,03·max(PVI − 15, 0) + 0,01·max(IVCI − 43, 0)
```

Limiares de PVI (Tsuchiya 2010) e IVCI (Zhang & Critchley 2016) da literatura; o grupo usava
13% e 40%.

> Como diz o README do grupo: Ce50, γ, Imax e os coeficientes de s são **ilustrativos** (ordem
> de grandeza da queda de PAM na indução), não ajustados a dados, e precisam ser calibrados.

## 6. Paciente virtual (extensão)

1. Rp = PAM/DC (como no grupo).
2. **C identificado** por bisseção em log C ∈ [0,2; 6] ml/mmHg para que a pressão de pulso do
   batimento periódico seja PAS − PAD (a PP cai quando C aumenta). É a mesma ideia de
   identificação do windkessel-id, com um único parâmetro e uma única medida.
3. Rc e L ficam nos valores de Stergiopulos. Com a linha arterial do paciente (pressão contínua
   + estimativa de fluxo), θ inteiro pode ser identificado com `identify`
   (`scripts/06_identificacao_windkessel.jl`).

## 7. Incerteza e coorte virtual (experimental, desativadas)

O código (`incerteza.jl`, `populacao.jl`) sorteia a farmacologia com variabilidade log-normal,
mas as variabilidades não estão validadas (a leitura da tabela de Schnider é incerta e as da
resposta hemodinâmica foram assumidas). Por isso não entram na plataforma nem nos resultados
até a validação ([validacao.md](validacao.md), seção 6).

## 8. Métricas

PAM basal (1 min antes do bolus), PAM mínima e instante, queda (%), tempo com PAM < 65 e
< 55 mmHg e área abaixo do limiar nos 20 min após o bolus, BIS mínimo e tempo até BIS < 60.
Risco pela PAM mínima prevista: **ALTO** (< 55 mmHg), **MODERADO** (55–65 mmHg), **BAIXO**
(≥ 65 mmHg). Os dois limiares vêm da literatura (Walsh 2013; Salmasi 2017).

## 9. Classificação (Resumo, Seções 2, 3 e 7)

| dicotomia | este modelo |
|---|---|
| tempo contínuo × discreto | PK contínua (EDOs no NeoPKPD); Windkessel discretizado (ZOH) batimento a batimento |
| linear × não linear | Windkessel linear; acoplamento não linear (Hill) |
| determinístico × estocástico | determinístico (as análises estocásticas estão desativadas até a validação) |
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
hemodinâmico usa o Ce do PK escolhido, cujo ke0 é mais rápido que o observado para a pressão
(Kazama 1999).
