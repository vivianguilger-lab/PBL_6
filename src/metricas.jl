# =============================================================================
# metricas.jl — métricas de hipotensão pós-indução e relatório
# =============================================================================
#
# Hipotensão pós-indução: queda de PA nos primeiros 20 min após a indução.
# Limiares (docs/validacao.md):
#   - 55 mmHg: mesmo curtos períodos de PAM < 55 associam-se a lesão renal aguda
#     e lesão miocárdica (Walsh et al. 2013);
#   - 65 mmHg: PAM < 65 associa-se a lesão de órgãos (Salmasi et al. 2017; faixa
#     de 60–70 mmHg no consenso de Sessler et al. 2019).
# O relatório dá o tempo e a área abaixo de cada limiar (AUT, mmHg·min).
#
# Como este arquivo se encaixa no modelo: `simular_inducao` (calibracao.jl)
# gera a curva prevista, batimento a batimento. Aqui essa curva vira poucos
# números para o médico: PAM basal, PAM mínima (nadir), queda em %, tempo e
# área abaixo de 65 e 55 mmHg, BIS mínimo, tempo até BIS < 60 e o risco
# (ALTO/MODERADO/BAIXO). Esses números alimentam a plataforma (plataforma.jl),
# a janela de dose (dose.jl) e os scripts. A regra de risco é nossa, mas usa
# só os dois limiares com fonte (docs/validacao.md, seções 3 e 4.6).
# =============================================================================

# -----------------------------------------------------------------------------
# metricas: resume uma indução simulada.
# Entrada: `r`, um `Resultado` de `simular_inducao` (calibracao.jl):
#   r.bat  = valores por batimento (t_min em min, PAM/PAS/PAD em mmHg,
#            FC em bpm, DC em L/min, VS em ml, Rp em mmHg·s/ml);
#   r.pkpd = PK/PD na grade de 1 s (t em min, Cp e Ce em µg/ml, BIS).
# Opções: `janela` = minutos analisados após o bolus (20); `limiares` de PAM
#   em mmHg (65: Salmasi 2017; 55: Walsh 2013); `alvo_bis` = 60 (limite
#   superior da faixa 40–60 do NICE DG6). Fontes em docs/validacao.md, seção 4.6.
# Saída: um NamedTuple, uma "ficha" com campos nomeados (lidos como mt.pam_nadir).
# -----------------------------------------------------------------------------
"""
    metricas(r::Resultado; janela=20.0, limiares=(65.0, 55.0), alvo_bis=60.0)

Resumo clínico da curva prevista (batimentos nos `janela` minutos após o bolus).
"""
function metricas(r::Resultado; janela = 20.0, limiares = (65.0, 55.0), alvo_bis = 60.0)
    # Atalhos: `b` = dados por batimento; `k` = dados da PK/PD; `t0` = instante
    # do bolus (min). A simulação começa 1 min antes do bolus: esse minuto é a
    # linha de base (paciente ainda sem propofol).
    b, k = r.bat, r.pkpd
    t0 = r.t_inducao
    # Linha de base: posições dos batimentos antes do bolus. `findall(cond, v)`
    # devolve as posições de `v` que cumprem a condição; `t -> t < t0` é uma
    # função anônima ("dado t, teste se t < t0"). `media` tira a média nessas
    # posições; se não houver nenhuma, usa o primeiro valor.
    # `cond ? a : b` quer dizer "se cond, então a; senão, b".
    base = findall(t -> t < t0, b.t_min)
    media(v, idx) = isempty(idx) ? v[1] : mean(v[idx])
    pam0, fc0, dc0 = media(b.PAM, base), media(b.FC, base), media(b.DC, base)
    # Janela de análise: batimentos de t0 até t0 + janela (20 min por padrão),
    # o período da hipotensão pós-indução.
    jan = findall(t -> t0 <= t <= t0 + janela, b.t_min)
    pam = b.PAM[jan]
    # Cada batimento dura 60/FC s = 1/FC min. O ponto em `./` é o "broadcasting"
    # de Julia: aplica a operação a cada elemento do vetor.
    dtb = 1 ./ b.FC[jan]                              # duração de cada batimento [min]
    # Nadir: `argmin` dá a posição da menor PAM dentro da janela (`i`); `j` é a
    # mesma posição contada na série completa de batimentos.
    i = argmin(pam)
    j = jan[i]
    # Para cada limiar l (65 e 55 mmHg), guardados num Dict (dicionário:
    # associa chave => valor, aqui limiar => resultado):
    #   tempo_abaixo[l] = soma da duração dos batimentos com PAM < l [min];
    #   aut[l] = área abaixo do limiar (AUT): soma de (l − PAM) × duração, só
    #            quando PAM < l [mmHg·min]; junta profundidade e duração.
    # `init = 0.0` faz a soma valer 0 quando nenhum batimento fica abaixo.
    tempo_abaixo = Dict(l => sum(dtb[pam .< l]; init = 0.0) for l in limiares)
    aut = Dict(l => sum(max.(l .- pam, 0.0) .* dtb; init = 0.0) for l in limiares)
    # BIS na mesma janela (grade da PK/PD, um ponto por segundo). `ib` = primeiro
    # ponto com BIS < alvo (hipnose adequada); `findfirst` devolve `nothing` se
    # isso nunca acontece. `<(alvo_bis)` é a função "x < alvo_bis".
    # `dt` = passo da grade em min (1/60), para converter contagem de pontos em tempo.
    js = findall(t -> t0 <= t <= t0 + janela, k.t)
    bis = k.bis[js]
    ib = findfirst(<(alvo_bis), bis)
    dt = length(k.t) > 1 ? k.t[2] - k.t[1] : 0.0
    # Menor e maior limiar (55 e 65 mmHg com os valores padrão).
    lmin, lmax = minimum(limiares), maximum(limiares)
    # Regra (docs/validacao.md, seções 3 e 4.6): PAM mínima < 55 => ALTO;
    # 55 a 65 => MODERADO; >= 65 => BAIXO. Não usa a queda em % nem o tempo
    # abaixo do limiar: esses critérios foram retirados por não terem fonte.
    # risco pela PAM mínima prevista frente aos dois limiares
    risco = pam[i] < lmin ? "ALTO" : pam[i] < lmax ? "MODERADO" : "BAIXO"
    # Resultado (NamedTuple). Campos principais:
    #   pam_basal, fc_basal, dc_basal = médias do minuto antes do bolus;
    #   pam_nadir [mmHg] e t_nadir [min após o bolus]; queda_pct = queda em % da
    #     PAM basal; PAS, PAD, FC, DC, Rp e VS no instante do nadir;
    #   bis_min; t_hipnose = min do bolus até BIS < alvo (NaN se não atinge);
    #   tempo_hipnose = min com BIS < alvo; cp_max e ce_max [µg/ml] na simulação toda.
    return (pam_basal = pam0, fc_basal = fc0, dc_basal = dc0,
            pam_nadir = pam[i], t_nadir = b.t_min[j] - t0, queda_pct = 100 * (1 - pam[i] / pam0),
            pas_nadir = b.PAS[j], pad_nadir = b.PAD[j], fc_no_nadir = b.FC[j], dc_no_nadir = b.DC[j],
            rp_no_nadir = b.Rp[j], vs_no_nadir = b.VS[j],
            tempo_abaixo = tempo_abaixo, aut = aut,
            bis_min = minimum(bis), t_hipnose = ib === nothing ? NaN : k.t[js][ib] - t0,
            tempo_hipnose = count(<(alvo_bis), bis) * dt,
            cp_max = maximum(k.Cp), ce_max = maximum(k.Ce),
            risco = risco, limiares = limiares, alvo_bis = alvo_bis)
end

# -----------------------------------------------------------------------------
# relatorio: imprime as métricas como texto (tela ou arquivo); não calcula
# nada novo, só formata. Entrada: `io` = onde escrever (`stdout` é a tela) e
# `mt` = saída de `metricas`. `@printf` escreve com formato: %6.1f = número
# com 1 casa decimal ocupando 6 espaços; %.0f = sem casas; %% imprime o
# símbolo %. `"─"^72` repete o traço 72 vezes (linha separadora).
# -----------------------------------------------------------------------------
"""Imprime um relatório legível das métricas."""
function relatorio(io::IO, mt; titulo = "Previsão hemodinâmica")
    println(io, "─"^72)
    println(io, titulo)
    println(io, "─"^72)
    @printf(io, "PAM basal            : %6.1f mmHg   (FC %.0f bpm, DC %.2f L/min)\n", mt.pam_basal, mt.fc_basal, mt.dc_basal)
    @printf(io, "PAM mínima           : %6.1f mmHg   em %.1f min (queda de %.0f%%)\n", mt.pam_nadir, mt.t_nadir, mt.queda_pct)
    @printf(io, "  nesse momento      : PA %.0f/%.0f, FC %.0f bpm, DC %.2f L/min\n", mt.pas_nadir, mt.pad_nadir, mt.fc_no_nadir, mt.dc_no_nadir)
    # Uma linha por limiar, do maior (65) para o menor (55):
    # `sort(...; rev = true)` ordena em ordem decrescente.
    for l in sort(collect(mt.limiares); rev = true)
        @printf(io, "Tempo com PAM < %2.0f   : %6.1f min    AUT = %.0f mmHg·min\n", l, mt.tempo_abaixo[l], mt.aut[l])
    end
    @printf(io, "BIS mínimo           : %6.0f        (BIS < %.0f em %.1f min; por %.1f min)\n",
            mt.bis_min, mt.alvo_bis, mt.t_hipnose, mt.tempo_hipnose)
    @printf(io, "Ce máxima            : %.2f µg/ml\n", mt.ce_max)
    println(io, "Risco previsto       : ", mt.risco)
    println(io, "─"^72)
end
# Versão curta: `relatorio(mt)` escreve direto na tela. `kw...` repassa as
# opções nomeadas (ex.: titulo = "...") para a versão completa acima.
relatorio(mt; kw...) = relatorio(stdout, mt; kw...)
