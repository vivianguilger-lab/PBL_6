# =============================================================================
# dose.jl — exploração de doses "antes da injeção"
#
# Pergunta clínica: "que dose usar neste paciente?", respondida no computador,
# antes da injeção. Cada dose candidata é simulada no paciente virtual
# (calibracao.jl) e resumida por `metricas` (metricas.jl). Três ferramentas:
#   - varrer_doses: tabela de métricas para uma lista de doses;
#   - bissecao: busca numérica do ponto em que um critério passa a valer;
#   - janela_de_dose: menor dose que leva o BIS < 60 e maior dose que mantém
#     a PAM mínima >= 65 mmHg.
# Doses em mg/kg (bolus). O efeito do propofol na pressão usa parâmetros
# ilustrativos, ainda sem calibração com pacientes: a janela é exploratória,
# não uma prescrição (docs/validacao.md, seções 4.4 e 6.2).
# =============================================================================

# -----------------------------------------------------------------------------
# varrer_doses: para cada dose `d` da lista `doses` (mg/kg), simula a indução
# por `janela` min e calcula as métricas. Devolve um vetor, uma linha por dose.
# Sintaxe: `f(x) = expressão` define a função numa linha; `[... for d in doses]`
# monta a lista repetindo a conta para cada d; `merge((dose = d,), mt)` junta o
# campo `dose` às métricas num só NamedTuple; `kw...` repassa outras opções
# (ex.: duracao_bolus em s, manutencao_ugkgmin em µg/kg/min).
# -----------------------------------------------------------------------------
"""
    varrer_doses(pv, doses; janela=20.0, kw...) -> Vector de NamedTuple

Simula cada dose (mg/kg) e devolve as métricas. `kw` é repassado a
`simular_inducao` (ex.: `duracao_bolus`, `manutencao_ugkgmin`).
"""
varrer_doses(pv, doses; janela = 20.0, kw...) =
    [merge((dose = d,), metricas(simular_inducao(pv; dose_mgkg = d, duracao = janela, kw...); janela)) for d in doses]

# -----------------------------------------------------------------------------
# bissecao: método da bisseção. Recebe um critério `f` (verdadeiro/falso) e um
# intervalo [a, b]; supõe que o critério é falso nas doses baixas e verdadeiro
# a partir de certo ponto. A cada passo testa o meio do intervalo e descarta a
# metade que não contém a fronteira, até o intervalo ficar menor que `tol`
# (0,01 mg/kg). Devolve o lado verdadeiro (`b`), ou NaN ("não é um número")
# se o critério não vale nem na ponta superior.
# -----------------------------------------------------------------------------
"""
    bissecao(f, a, b; tol=0.01)

Menor `x` em [a, b] com f(x) == true, supondo f monótona (false…false true…true).
"""
function bissecao(f, a, b; tol = 0.01)
    # Casos de borda: já vale na ponta inferior => devolve `a`; não vale nem na
    # superior => não há fronteira no intervalo (NaN).
    # `cond && x` executa x só se cond for verdadeira; `cond || x`, só se for falsa.
    f(a) && return a
    f(b) || return NaN
    # Laço: repete enquanto o intervalo for maior que a tolerância. Se f(c) é
    # verdadeiro, a fronteira está à esquerda de c (b = c); senão, à direita (a = c).
    while b - a > tol
        c = (a + b) / 2
        f(c) ? (b = c) : (a = c)
    end
    return b
end

# -----------------------------------------------------------------------------
# janela_de_dose: janela terapêutica individual do paciente virtual `pv`.
#   dose_min = menor dose com BIS mínimo < 60 nos 20 min (hipnose; faixa-alvo
#              40–60 do NICE DG6);
#   dose_max = maior dose com PAM mínima >= 65 mmHg (Salmasi 2017).
# Busca em `faixa` = 0,5 a 3,5 mg/kg. Se dose_min <= dose_max, a janela existe
# (`existe = true`); senão, no modelo nenhuma dose dá hipnose sem PAM < 65.
# Critério da dose mínima em revisão: com a manutenção de 100 µg/kg/min, até
# 0,5 mg/kg chega a BIS < 60 depois de vários minutos (docs/validacao.md, 6.2).
# -----------------------------------------------------------------------------
"""
    janela_de_dose(pv; limiar_pam=65.0, alvo_bis=60.0, faixa=(0.5, 3.5), kw...)

Janela terapêutica individual (mg/kg): `dose_min` é a menor dose que leva o BIS
abaixo de `alvo_bis`; `dose_max` é a maior dose cujo nadir de PAM fica ≥ `limiar_pam`.
"""
function janela_de_dose(pv; limiar_pam = 65.0, alvo_bis = 60.0, faixa = (0.5, 3.5), janela = 20.0, kw...)
    # `sim(d)`: função auxiliar local que simula a dose d e devolve as métricas.
    sim(d) = metricas(simular_inducao(pv; dose_mgkg = d, duracao = janela, kw...); janela, alvo_bis)
    # Menor dose que leva o BIS mínimo abaixo do alvo. `faixa...` "desempacota"
    # a dupla (0.5, 3.5) nos dois argumentos a e b; `d -> ...` é função anônima.
    dmin = bissecao(d -> sim(d).bis_min < alvo_bis, faixa...)
    # Menor dose em que a PAM mínima já fica abaixo de 65 mmHg (supõe que doses
    # maiores baixam mais a PAM). A dose máxima é essa fronteira (precisão de
    # 0,01 mg/kg); se nenhuma dose da faixa viola (NaN), vale o topo da faixa.
    dmax_viola = bissecao(d -> sim(d).pam_nadir < limiar_pam, faixa...)
    dmax = isnan(dmax_viola) ? faixa[2] : dmax_viola == faixa[1] ? NaN : dmax_viola
    # `existe`: há janela se as duas doses foram achadas e a mínima não passa da máxima (dose_max = NaN se até a menor dose já derruba a PAM).
    return (dose_min = dmin, dose_max = dmax, existe = !isnan(dmin) && !isnan(dmax) && dmin <= dmax,
            limiar_pam = limiar_pam, alvo_bis = alvo_bis)
end
