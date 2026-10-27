# =============================================================================
# incerteza.jl — previsão individual com faixa de incerteza
# =============================================================================
#
# EXPERIMENTAL — NÃO VALIDADO: usa as variabilidades de populacao.jl, que ainda
# não têm validação (docs/validacao.md). Fora da plataforma e dos resultados.
#
# O paciente virtual tem a SUA hemodinâmica basal (calibrada com os dados
# pré-indução), mas a sensibilidade dele ao propofol é desconhecida antes da
# injeção. Mantemos a calibração e sorteamos a PK (Schnider) e a resposta
# hemodinâmica (Ce50 e Imax do modelo do grupo) com a variabilidade da população.
#
# Ideia: o mesmo paciente é simulado n vezes, cada vez com uma farmacologia
# sorteada; a dispersão das curvas de PAM daria uma faixa (percentis 5–95%)
# e a chance de PAM < 65 ou < 55 mmHg. ATENÇÃO: essa faixa e essas chances NÃO
# são confiáveis, porque as variabilidades (OMEGA_PK e OMEGA_HEMO, em
# populacao.jl) foram assumidas ou têm leitura incerta (docs/validacao.md,
# seção 6.3). Por isso a plataforma mostra só a previsão do paciente típico.
# Não use estes números para decisões clínicas.
# =============================================================================

# -----------------------------------------------------------------------------
# previsao_com_incerteza (EXPERIMENTAL, NÃO VALIDADO)
# Entrada: paciente virtual `pv`; `n` = número de sorteios (200); `rng` =
# gerador de números aleatórios com semente fixa (Xoshiro(1)), para o resultado
# se repetir igual; `dt` = passo da grade comum de tempo (0,1 min); `janela`
# em min; `kw...` vai para `simular_inducao` (ex.: dose_mgkg).
# Saída (NamedTuple): t (min após o bolus); mediana, p05 e p95 da PAM [mmHg];
# PAM mínima mediana e seus percentis 5 e 95; fração dos sorteios com PAM < 65,
# com PAM < 55 e com BIS < 60 (prob_*), e n.
# -----------------------------------------------------------------------------
"""
    previsao_com_incerteza(pv; n=200, rng=Xoshiro(1), dt=0.1, janela=20.0, kw...)

Roda `n` simulações do mesmo paciente virtual com PK/PD sorteadas. Devolve a
grade de tempo (min após o bolus), mediana e percentis 5/95 da PAM e as
probabilidades de hipotensão.
"""
function previsao_com_incerteza(pv::PacienteVirtual; n = 200, rng = Xoshiro(1), dt = 0.1, janela = 20.0, kw...)
    # Sorteia antes do laço, em sequência, os n pares (fatores da PK, efeito
    # hemodinâmico), para o resultado não depender da ordem das threads.
    # `_` é uma variável de laço que não é usada.
    amostras = [(amostrar_eta_pk(rng), amostrar_efeito(rng, pv.efeito)) for _ in 1:n]
    # Bolus no minuto 1 (o minuto 0–1 é a linha de base). Grade comum de tempo
    # para comparar as n curvas; `curvas` é uma matriz tempo × sorteio. Vetores de
    # resultado por sorteio: PAM mínima e marcas verdadeiro/falso (`falses(n)`).
    t_inducao = 1.0
    grade = collect(0.0:dt:(t_inducao + janela))
    curvas = zeros(length(grade), n)
    nadir = zeros(n); abaixo65 = falses(n); abaixo55 = falses(n); bis_ok = falses(n)
    # `Threads.@threads` divide os sorteios entre os núcleos do processador
    # (com `julia -t auto`); cada k escreve só na sua coluna/posição, sem conflito.
    Threads.@threads for k in 1:n
        # PK do paciente com os fatores sorteados (`eta`) e simulação com o efeito
        # sorteado (`ef`); a calibração basal (Rp, C) continua a do paciente.
        # Em `f(x; pk, ...)`, um nome sem "= valor" passa a variável de mesmo nome.
        eta, ef = amostras[k]
        pk = parametros_pk(pv.paciente, pv.pk.modelo; eta)
        r = simular_inducao(pv; pk, efeito = ef, t_inducao, duracao = janela, kw...)
        mt = metricas(r; janela)
        # Leva a PAM por batimento para a grade comum e registra os desfechos.
        curvas[:, k] = interpolar(r.bat.t_min, r.bat.PAM, grade)
        nadir[k] = mt.pam_nadir
        abaixo65[k] = mt.tempo_abaixo[65.0] > 0
        abaixo55[k] = mt.tempo_abaixo[55.0] > 0
        bis_ok[k] = mt.bis_min < 60
    end
    # `q(p)`: percentil p da PAM em cada instante, entre os n sorteios (`view`
    # olha a linha i da matriz sem copiá-la). As probabilidades são a média das
    # marcas verdadeiro/falso, ou seja, a fração dos sorteios (não validadas).
    q(p) = [quantile(view(curvas, i, :), p) for i in eachindex(grade)]
    return (t = grade .- t_inducao, mediana = q(0.5), p05 = q(0.05), p95 = q(0.95),
            nadir_mediano = median(nadir), nadir_p05 = quantile(nadir, 0.05), nadir_p95 = quantile(nadir, 0.95),
            prob_pam65 = mean(abaixo65), prob_pam55 = mean(abaixo55), prob_hipnose = mean(bis_ok), n = n)
end

# interpolar: aplica `interp` (acoplamento.jl, interpolação linear) a cada ponto
# de `xq`, levando a curva irregular (um ponto por batimento) para a grade comum.
"""Interpolação linear de (x, y) nos pontos `xq` (extrapola constante)."""
interpolar(x, y, xq) = [interp(x, y, v) for v in xq]
