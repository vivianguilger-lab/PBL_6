# Modelo matemático

> **Pergunta:** *como prever, por simulação, a resposta hemodinâmica de um paciente ao propofol?*
>
> **Resposta proposta:** um paciente virtual (modelo 0-D da circulação em malha fechada)
> calibrado com os dados pré-indução do próprio paciente, acoplado à farmacocinética de
> Schnider e a curvas de Hill para cada efeito hemodinâmico do propofol. O modelo prevê a
> curva de PAM para a dose escolhida **antes** da injeção.

## 1. Visão geral (diagrama de blocos)

```
 dose (mg/kg), tempo de injeção
            │  u(t) [mg/s]
            ▼
 ┌─────────────────────┐ Cp  ┌──────────────┐ Ce_bis ──► BIS (hipnose)        ──┐
 │ PK — Schnider       │────►│ compartimento│        ──► ganho do barorreflexo ↓ │
 │ 3 compartimentos    │     │ de efeito    │ Ce_hemo──► TPR ↓ (vasodilatação)   │ curvas de Hill
 └─────────────────────┘     └──────────────┘        ──► Emax ↓ (inotropismo)    │
                                                     ──► V0 venoso ↑ (venodilat.)┘
                                                              │
                                                              ▼
 ┌──────────────────────────────────────────────────────────────────────────────┐
 │ Circulação 0-D em malha fechada: VE (elastância variável) → valva aórtica →  │
 │ Windkessel 3 elementos (Zc, AC, TPR) → veias (VC, V0v) → enchimento → VE      │
 │                ▲                                            │                 │
 │                └──── barorreflexo (FC, TPR, Emax, V0v) ◄────┘ PA sensoriada   │
 └──────────────────────────────────────────────────────────────────────────────┘
                                │
                                ▼
            PAS, PAD, PAM, FC, DC batimento a batimento → métricas de hipotensão
```

Variáveis de estado (16): massas de propofol nos 3 compartimentos, `Ce_bis`, `Ce_hemo`,
fase do ciclo cardíaco, volumes do ventrículo, das artérias e das veias, pressão sensoriada
pelo barorreceptor, 4 ativações dos efetores do barorreflexo e duas integrais auxiliares
(∫Pao dt e ∫Qao dt, para médias exatas por batimento).

## 2. Farmacocinética — Schnider (1998/1999)

Modelo tricompartimental mamilar, massas `A₁, A₂, A₃` (mg), infusão `u(t)` (mg/s):

```
dA₁/dt = u(t) − (k₁₀ + k₁₂ + k₁₃)·A₁ + k₂₁·A₂ + k₃₁·A₃
dA₂/dt = k₁₂·A₁ − k₂₁·A₂
dA₃/dt = k₁₃·A₁ − k₃₁·A₃
Cp     = A₁ / V₁                                  (mg/L = µg/ml)
```

| parâmetro | equação (idade em anos, peso em kg, altura em cm) |
|---|---|
| V₁ | 4,27 L |
| V₂ | 18,9 − 0,391·(idade − 53) L |
| V₃ | 238 L |
| Cl₁ | 1,89 + 0,0456·(peso − 77) − 0,0681·(LBM − 59) + 0,0264·(altura − 177) L/min |
| Cl₂ | 1,29 − 0,024·(idade − 53) L/min |
| Cl₃ | 0,836 L/min |
| LBM (James) | ♂ 1,1·peso − 128·(peso/altura)²; ♀ 1,07·peso − 148·(peso/altura)² |

`k₁₀ = Cl₁/V₁`, `k₁₂ = Cl₂/V₁`, `k₂₁ = Cl₂/V₂`, `k₁₃ = Cl₃/V₁`, `k₃₁ = Cl₃/V₃`.

**Dois compartimentos de efeito.** O propofol age no cérebro (hipnose) mais rápido do que
nos vasos: Kazama et al. (1999) mediram t½ke0 ≈ 2,3 min para o BIS, mas ≈ 5,7 min (20–39 anos)
a ≈ 10,2 min (70–85 anos) para a pressão arterial. Por isso:

```
dCe_bis/dt  = ke0_bis ·(Cp − Ce_bis)     ke0_bis = 0,456 /min (Schnider)
dCe_hemo/dt = ke0_hemo·(Cp − Ce_hemo)    ke0_hemo = ln2 / t½(idade), t½: 5,7 → 10,2 min (Kazama)
```

## 3. Farmacodinâmica — equação de Hill

```
                    Ceᵞ
E(Ce) = Emax · ─────────────
               C50ᵞ + Ceᵞ
```

| efeito clínico (do problema) | parâmetro do circuito | Emax | C50 (µg/ml) | γ | Ce usada | fonte |
|---|---|---|---|---|---|---|
| hipnose | BIS = 97,4·(1 − H) | −100% | 4,47 | 1,43 | Ce_bis | Bouillon 2004 |
| vasodilatação periférica | TPR × (1 + E) | −77,8% | 3,21 | 1,83 | Ce_hemo | Su 2023 |
| depressão do miocárdio | Emax do VE × (1 + E) | −15,4%·e^{0,0333(idade−35)} | 0,44 | 1 | Ce_hemo | Su 2023 (efeito no VS, dependente da idade) |
| venodilatação (↓ Pms) | V0 venoso × (1 + E) | +25% | 2,5 | 1,5 | Ce_hemo | **assumido**, qualitativamente de de Wit 2016 |
| atenuação do barorreflexo | ganho × (1 + E) | −80% | 2,5 | 2 | Ce_bis | **assumido**, qualitativamente de Cullen 1987 e Ebert 1992 |

Os dois parâmetros **assumidos** foram ajustados (calibração populacional) para que o
paciente típico (53 anos, 2 mg/kg) reproduza a resposta clássica descrita na literatura de
anestesiologia: PAM −25 a −40%, DC ≈ −15%, RVS −15 a −25%, FC pouco alterada
(ver `docs/validacao.md`). Eles são os primeiros candidatos a serem identificados com
dados reais.

## 4. Circulação — Windkessel em malha fechada

Analogia eletromecânica do curso: pressão ↔ tensão (*across*), fluxo ↔ corrente
(*through*), resistência ↔ resistor, complacência ↔ capacitor, valva ↔ diodo.

```
Ventrículo (capacitor variável):  P_VE = E(t)·(V_VE − V_d),   E(t) = E_min + (E_max − E_min)·e(t)
                                   C_D = 1/E_min (diástole),   C_S = 1/E_max (sístole)
Ativação (cossenoidal):            e(t) = ½(1 − cos(π t/T_s))           0 ≤ t < T_s
                                          ½(1 + cos(π (t−T_s)/T_r))     T_s ≤ t < T_s + T_r
                                          0                              resto do ciclo
                                   T_s = 0,35·√T, T_r = T_s/2, T = 60/FC
Artérias (capacitor AC):           P_c = (V_a − V_a0)/AC
Veias (capacitor VC):              P_v = (V_v − V_v0)/VC
Enchimento (diodo + R_in):         Q_in = max(P_v − P_VE, 0)/R_in
Ejeção (diodo + R_av + Z_c):       Q_ao = max(P_VE − P_c, 0)/(R_av + Z_c)
Pressão na raiz da aorta:          P_ao = P_c + Z_c·Q_ao
Microcirculação (resistor TPR):    Q_s = (P_c − P_v)/TPR

Conservação de volume (Kirchhoff dos nós):
dV_VE/dt = Q_in − Q_ao
dV_a/dt  = Q_ao − Q_s
dV_v/dt  = Q_s − Q_in
dφ/dt    = FC/60            (fase do ciclo; φ inteiro = início de um batimento)
```

Na diástole a valva aórtica está fechada (`Q_ao = 0`) e o compartimento arterial
descarrega como o **Windkessel de 2 elementos** do curso: `P(t) = P₀·e^{−t/τ}`, `τ = TPR·AC`
— isso é verificado numericamente (`scripts/01_verificacao.jl`).

Com o coração parado todas as pressões convergem para a **pressão média de enchimento
sistêmico** (Resumo, eq. 8), que o modelo reproduz exatamente:

```
Pms = (V_total − V_a0 − V_v0 − V_d) / (AC + VC + C_D)
```

O coração direito e a circulação pulmonar estão agregados ao ventrículo e a `R_in`
(topologia de um ventrículo, como em Bjørdalsbakke et al. 2022 /
`CirculatorySystemModels.jl`).

## 5. Barorreflexo (controle em malha fechada)

```
dP_s/dt = (P_ao − P_s)/τ_s                         sensor (barorreceptor), τ_s = 2 s
r       = tanh(−(P_s − P_set)/S)                   aferência normalizada, S = 25 mmHg
g       = 1 + E_baro(Ce_bis)                       ganho reduzido pelo propofol
dxᵢ/dt  = (g·r − xᵢ)/τᵢ                            efetores i ∈ {FC, TPR, Emax, V0v}

FC   = FC₀ ·(1 + 0,40·x_FC)      τ = 2 s
TPR  = TPR₀·(1 + E_TPR)·(1 + 0,35·x_R)    τ = 6 s
Emax = Emax₀·(1 + E_inot)·(1 + 0,25·x_E)  τ = 8 s
V0v  = V0v₀·(1 + E_veno)·(1 − 0,08·x_V)   τ = 20 s
```

`P_set` é a PAM basal do próprio paciente. Ganho de malha aberta ≈ 1,2 bpm/mmHg para a FC
(faixa fisiológica 0,5–2). Estrutura inspirada em Ursino (1998), simplificada.

| elemento de controle (Resumo, Seção 9.2) | no modelo |
|---|---|
| saída controlada `y` | PAM |
| sensor | barorreceptor (`P_s`) |
| setpoint `r` | PAM basal do paciente |
| controlador | centro vasomotor (tanh + ganho `g`) |
| atuadores | FC, TPR, contratilidade, tônus venoso |
| distúrbio `x` | propofol (vasodilatação, inotropismo negativo, venodilatação) |
| "defeito" causado pelo fármaco | ganho do controlador reduzido |

## 6. Equação de observação (manguito × aorta)

O modelo 0-D fornece a pressão **aórtica**. No braço a pressão de pulso é amplificada pela
reflexão de ondas — fenômeno que um modelo 0-D não reproduz (Resumo, Seção 12.4). A PAM
praticamente não muda ao longo da árvore arterial. Por isso:

```
amp(idade)   = 1,25 + 0,45·exp(−(idade − 20)/15)       (≈1,7 aos 20; ≈1,3 aos 50; Anglo-Cardiff)
PP_braquial  = amp · PP_aórtica
PAD_braquial = PAM − ff·PP_braquial                     ff = (PAM−PAD)/(PAS−PAD) do manguito
PAS_braquial = PAD_braquial + PP_braquial
```

## 7. Calibração do paciente virtual (grey-box)

Entradas: idade, sexo, peso, altura, PAS, PAD, FC (e, se houver, PAM e DC medidos).

1. **Tamanho:** volume circulante = 0,85 × volemia de Nadler; elastâncias e volume morto
   escalados pela superfície corporal (Mosteller).
2. **Alvos:** PAM, PP aórtica = (PAS − PAD)/amp(idade), DC (medido ou índice cardíaco
   estimado × superfície corporal), FC medida.
3. **Primeira estimativa pelas equações do curso:** `TPR ≈ (PAM − Pv)/DC`, `AC ≈ VS/PP`.
4. **Ajuste de ponto fixo** sobre o modelo pulsátil completo (sem fármaco, sem barorreflexo):
   - `TPR ← TPR·[(PAM* − Pv)/(PAM − Pv)]/(DC*/DC)`
   - `AC ← AC·(PP/PP*)^1,2` (se AC satura em 4 ml/mmHg, ajusta-se `Z_c`)
   - `V0v ← V0v + 0,7·e_DC·V_estressado/(1 + e_DC)` (o volume estressado define Pms e o retorno venoso)
   até erro < 0,5% em PAM, PP e DC.
5. O estado de regime permanente vira a condição inicial; `P_set` = PAM obtida.

## 8. Métricas de saída

- PAM basal, **nadir**, tempo até o nadir, queda percentual;
- tempo com PAM < 65 e < 55 mmHg e área abaixo do limiar (AUT, mmHg·min) nos 20 min;
- BIS mínimo, tempo até BIS < 60;
- risco previsto: **ALTO** (qualquer PAM < 55), **MODERADO** (PAM < 65 por ≥ 1 min ou queda > 30%),
  **BAIXO** (demais).

## 9. Classificação do modelo (Resumo, Seções 2, 3 e 7)

| dicotomia | este modelo |
|---|---|
| tempo contínuo × discreto | contínuo (EDOs); as saídas são lidas batimento a batimento |
| linear × não linear | não linear (valvas = diodos, Hill, tanh do barorreflexo) |
| determinístico × estocástico | determinístico por paciente; estocástico na coorte (IIV log-normal) |
| tempo-invariante × tempo-variante | tempo-variante (E(t) do ventrículo, efeitos do fármaco) |
| concentrado × distribuído | concentrado (0-D) |
| black box × grey box | grey box: estrutura física, parâmetros identificados do paciente |
| paramétrico × não paramétrico | paramétrico |

- **Objetivo da modelagem:** preditivo (com componente explicativa — cada parâmetro é um processo fisiológico).
- **Tipologias de uso:** PARAMETER IDENTIFICATION (calibração), SIMULATION ("e se…?" de doses) e FORECAST (curva prevista).

## 10. Parâmetros cardiovasculares de referência (antes da calibração)

| símbolo | valor | unidade | observação |
|---|---|---|---|
| E_max (= 1/C_S) | 2,0 | mmHg/ml | C_S = 0,5 ml/mmHg, valor do Resumo; escala com 1/superfície corporal |
| E_min (= 1/C_D) | 0,06 | mmHg/ml | C_D ≈ 17 ml/mmHg |
| V_d | 10 | ml | volume morto |
| R_in | 0,02 | mmHg·s/ml | enchimento (veias centrais + valva AV) |
| R_av | 0,005 | mmHg·s/ml | valva aórtica aberta |
| Z_c | 0,03 | mmHg·s/ml | impedância característica (calibrada se necessário) |
| AC | 1,5 | ml/mmHg | **calibrada** |
| TPR | 1,0 | mmHg·s/ml | **calibrada** |
| VC | 100 | ml/mmHg | valor do Resumo |
| V_a0 | 16,5% do volume circulante | ml | volume não estressado arterial |
| V_v0 | — | ml | **calibrado** (define volume estressado e Pms) |
