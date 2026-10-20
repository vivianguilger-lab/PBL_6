# Credibilidade do modelo — ASME V&V 40

Organizado pelas quatro fases do V&V 40 (Resumo, Seção 8.5):
**planejamento → avaliação de risco → verificação → validação**.

## 1. Planejamento — Contexto de Uso (COU)

> Ferramenta **educacional e de pesquisa** para explorar, **antes da injeção**, a curva de
> PAM prevista nos primeiros 20 min após um bolus de propofol na indução de adultos
> (18–90 anos), sem estímulo cirúrgico, a partir dos dados pré-indução do próprio paciente
> (idade, sexo, peso, altura, PAS, PAD, FC e, se houver, DC). Saídas: curva de PAM,
> nadir, tempo abaixo de 65/55 mmHg, BIS e janela de dose.
>
> **Fora do COU:** decidir a dose de um paciente real; crianças; obesos graves (a fórmula
> de massa magra de James usada por Schnider falha); choque, sepse, cardiopatia grave;
> indução com opioides/outros hipnóticos associados; manutenção com halogenados.

## 2. Avaliação de risco

| fator | avaliação |
|---|---|
| influência do modelo na decisão | baixa no COU declarado (apoio didático; a decisão continua com o anestesista) |
| consequência de uma decisão errada | alta (hipotensão grave → lesão miocárdica e renal) |
| **risco do modelo** | **médio-alto** → exige validação com dados humanos independentes antes de qualquer uso clínico |

## 3. Verificação — "resolvi as equações corretamente?"

`scripts/01_verificacao.jl` e `test/runtests.jl` (68 testes, rodados no GitHub Actions a cada *push*):

| teste | solução de referência | erro obtido |
|---|---|---|
| Windkessel de 2 elementos na diástole | `P(t) = P₀·e^{−t/τ}`, τ = TPR·AC = 2 s | 3 × 10⁻⁸ (relativo) |
| PK de Schnider após bolus | `A(t) = exp(K·t)·A₀` (exponencial de matriz) | 3 × 10⁻¹³ (relativo) |
| parada cardíaca | todas as pressões → `Pms = (Vt − V0)/Ct` (eq. 8 do Resumo) | 4 × 10⁻⁷ mmHg |
| conservação de massa (volume de sangue) | volume total constante | 8 × 10⁻¹² ml |
| calibração | PAM, PP aórtica e DC reproduzidos | < 1% |
| linha de base sem fármaco | PAM e FC estáveis por 5 min | < 1,5 mmHg / 1,5 bpm |
| ponto de funcionamento | interseção satisfaz as duas retas do curso | 10⁻¹⁰ |
| verificação cruzada | mesmo circuito montado com `CirculatorySystemModels.jl` (ModelingToolkit), `extras/` | 10⁻⁴ mmHg |

## 4. Validação — "resolvi as equações corretas?"

### 4.1 O que já foi feito (validade de face e plausibilidade)

| evidência | literatura | modelo |
|---|---|---|
| paciente típico (53 a, 2 mg/kg): queda da PAM | −25 a −40% | −35% |
| … débito cardíaco | ≈ −15% | −15% |
| … resistência vascular sistêmica | −15 a −25% | −24% |
| … frequência cardíaca | pouco alterada | +5 a +8 bpm |
| tempo até o nadir após bolus | primeiros minutos | 2,5–3,5 min |
| efeito na PA mais lento que no BIS (Kazama 1999) | t½ 5,7–10,2 × 2,3 min | imposto por construção |
| idade ↑ ⇒ queda maior com a mesma dose/kg | sim | sim (coorte com 2 mg/kg fixo) |
| dose ↑ ⇒ nadir ↓ (monotonicidade) | sim | sim (teste automatizado) |
| propofol ↓ Pms (de Wit 2016) | sim | sim (venodilatação: Pms −10 a −20%) |
| **coorte virtual (500 adultos, dose ajustada por idade): PAM < 55 mmHg** | **17–24%** | **≈25%** |

**Atenção — calibração ≠ validação (Resumo, Seção 5.1):** os dois parâmetros assumidos
(venodilatação e atenuação do barorreflexo) foram ajustados para reproduzir as quatro
primeiras linhas da tabela. A incidência populacional de PAM < 55 **não** foi usada no
ajuste — é uma checagem semi-independente, mas ainda não é validação com dados de pacientes.

### 4.2 O que falta (proposta)

1. **Dados individuais independentes:** o banco aberto **VitalDB** (vitaldb.net) tem milhares de
   cirurgias com registro de infusão de propofol (TCI) e pressão arterial invasiva/não invasiva.
   Plano: calibrar o paciente virtual com os 5 min pré-indução de cada caso, prever a curva dos
   20 min seguintes e comparar (erro do nadir, RMSE da curva, AUROC para PAM < 65 contra um
   modelo só com dados de base). Separar os casos usados para ajustar os parâmetros assumidos
   dos casos usados para validar.
2. **Linha arterial pré-indução:** usar a forma de onda para identificar AC, TPR e Zc
   (ver `hpigot/windkessel-id`) em vez de estimá-los só por PAS/PAD.
3. **Análise de sensibilidade global** (Sobol/Morris) além da local (`scripts/05_sensibilidade.jl`).

## 5. Limitações conhecidas

- Sem opioide de coindução (fentanil/remifentanil potencializam a hipotensão), sem estímulo de
  laringoscopia/intubação (eleva a PA), sem fluidos ou vasopressores — todos entram facilmente
  como novas entradas do modelo.
- Coração de um ventrículo, sem circulação pulmonar explícita; modelo 0-D sem propagação de
  ondas (a amplificação periférica é tratada pela equação de observação).
- Schnider: V1 fixo produz Cp inicial muito alta após bolus; não serve para obesos graves.
- Os C50 de Su et al. foram estimados sobre a concentração plasmática num modelo de
  *turnover*; aqui são aplicados à Ce vascular (equilíbrio de Kazama) — uma reinterpretação.
- O efeito de Su sobre o volume sistólico foi mapeado inteiramente para a contratilidade (Emax).
- O DC pré-indução, quando não medido, é estimado (±10% no DC ⇒ ±1,2 mmHg no nadir).
- A sensibilidade individual ao fármaco é desconhecida antes da injeção: por isso a previsão
  individual é dada com faixa de incerteza (`previsao_com_incerteza`), não como um número único.
