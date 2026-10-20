# Credibilidade do modelo — ASME V&V 40

Organizado pelas quatro fases do V&V 40 (Resumo, Seção 8.5):
**planejamento → avaliação de risco → verificação → validação**.

## 1. Planejamento — Contexto de Uso (COU)

> Ferramenta **educacional e de pesquisa** para explorar, **antes da injeção**, a curva de PAM
> prevista nos primeiros 20 min após um bolus de propofol (com ou sem manutenção) em adultos
> (18–95 anos), sem estímulo cirúrgico, a partir dos biomarcadores pré-indução do paciente
> (idade, sexo, peso, altura, ASA, PAS, PAD, FC, PVI, IVCI e, se houver, DC).
>
> **Fora do COU:** decidir a dose de um paciente real; crianças; obesos graves (a massa magra de
> James usada por Schnider falha); choque, sepse, cardiopatia grave; coindução com opioides;
> manutenção com halogenados.

## 2. Avaliação de risco

| fator | avaliação |
|---|---|
| influência do modelo na decisão | baixa no COU declarado (apoio didático) |
| consequência de uma decisão errada | alta (hipotensão grave → lesão miocárdica e renal) |
| **risco do modelo** | **médio-alto** → exige calibração e validação com dados humanos antes de qualquer uso clínico |

## 3. Verificação — "resolvi as equações corretamente?"

`scripts/01_verificacao.jl` e `test/runtests.jl` (47 testes, inclusive os do PropofolPBL do grupo,
rodados no GitHub Actions):

| teste | referência | resultado |
|---|---|---|
| Windkessel de 2 elementos (G2) descarregando | `P₀·e^{−t/(Rp·C)}`, τ = 2 s | erro 5 × 10⁻¹⁴ |
| NeoPKPD (Schnider, bolus) | `exp(K·t)·A₀` | erro 8 × 10⁻¹⁰ |
| Windkessel de 4 elementos | PAM = Rp·Q̄; estado periódico de `solvex0` | 8 × 10⁻¹⁶; 2 × 10⁻¹⁶ |
| `bode_wk` | \|Z(0)\| = Rp, \|Z(∞)\| = Rc | < 10⁻³ |
| `identify` (dados sintéticos) | recupera θ verdadeiro | < 10⁻³ |
| paciente virtual | PAM e PP do monitor reproduzidas | < 1% |
| linha de base sem fármaco | PAM estável | < 0,5 mmHg |

## 4. Validação — "resolvi as equações corretas?"

### 4.1 O que existe

| evidência | literatura | modelo (parâmetros do grupo) |
|---|---|---|
| identificação em dados humanos (Stergiopulos 1999) | θ publicado [0,79; 1,22; 0,056; 0,0051] | [0,80; 1,25; 0,056; 0,0055], RMSE 2,4 mmHg |
| idade ↑, ASA ↑, PVI/IVCI ↑ ⇒ queda maior | sim (fatores de risco de hipotensão pós-indução) | sim (fator *s*) |
| dose ↑ ⇒ PAM mínima ↓ | sim | sim (teste automatizado) |
| hipnose em < 1 min após bolus | sim | sim (BIS < 60 em ~0,4 min) |
| **coorte virtual: PAM < 55 mmHg** | **17–24%** | **35%** (superestima) |
| queda típica da PAM após 2 mg/kg | 25–40% | mediana 39% na coorte (alto da faixa) |

### 4.2 Calibração populacional sugerida (não aplicada)

Como os parâmetros do efeito hemodinâmico são ilustrativos, testamos uma escala comum nos três
Imax (Rp, VS, FC) na coorte de 500 adultos com dose ajustada por idade:

| Imax × | PAM < 65 | PAM < 55 | queda mediana |
|---|---|---|---|
| 1,0 (grupo) | 75% | 35% | 38% |
| 0,8 | 49% | 14% | 31% |
| 0,7 | 35% | 7% | 27% |
| 0,6 | 21% | 3% | 24% |

Uma escala ≈ 0,85 reproduziria ≈ 20% de PAM < 55. Fica como **proposta para o grupo**: é
calibração com um dado populacional, não validação — a validação precisa de dados independentes.

### 4.3 O que falta

1. **Dados individuais independentes** (ex.: banco aberto VitalDB, com infusão de propofol e
   pressão arterial): calibrar Ce50, γ, Imax e os coeficientes de *s* num subconjunto e validar
   em outro (erro da PAM mínima, RMSE da curva, AUROC para PAM < 65).
2. **Linha arterial pré-indução** do paciente para identificar θ inteiro com `identify`
   (`scripts/06_identificacao_windkessel.jl`).
3. Análise de sensibilidade global (Sobol/Morris) além da local (`scripts/05_sensibilidade.jl`).

## 5. Limitações conhecidas

- Parâmetros do efeito e da sensibilidade ilustrativos (ver 4.2).
- Fluxo aórtico imposto (modelo aberto): sem barorreflexo nem retorno venoso explícitos; a queda
  de VS agrega depressão miocárdica e redução de pré-carga.
- Um único ke0 para BIS e pressão; na prática o efeito sobre a pressão é mais lento (Kazama 1999).
- Schnider com V1 fixo superestima o pico de concentração após bolus; BIS de Eleveld com o Ce do
  ke0 do modelo PK escolhido (aproximação do grupo).
- Forma de onda fixa (Stergiopulos): PAS e PAD individuais não coincidem com o manguito.
- SpO₂ e glicemia ainda não entram no modelo; sem opioides, laringoscopia, fluidos ou vasopressores.
