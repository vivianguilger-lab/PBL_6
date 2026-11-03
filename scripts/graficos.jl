# =============================================================================
# graficos.jl — figuras do projeto (Plots.jl + GR)
# Incluído pelos scripts: include(joinpath(@__DIR__, "graficos.jl"))
# No REPL: include(joinpath(pkgdir(PropofolWindkessel), "scripts", "graficos.jl"))
#
# Este arquivo não é rodado sozinho: ele só define cores, estilo, a pasta de
# saída e as funções que desenham as figuras. Os scripts 01–06 o incluem.
# Biblioteca de gráficos: Plots.jl, com o motor GR.
# =============================================================================
# rodando como script, o GR só grava arquivos; no REPL, abre janelas
# `a || b` só executa b se a for falso. Aqui: fora do REPL, define a variável
# GKSwstype = "100" (modo sem janela), se ela ainda não estiver definida.
isinteractive() || (ENV["GKSwstype"] = get(ENV, "GKSwstype", "100"))
# Plots = gráficos; Printf = @printf e @sprintf (texto com formato)
using Plots
using Printf
# escolhe o motor GR para desenhar
gr()

# Paleta: slots categóricos em ordem fixa + cores de status reservadas
# COR é uma NamedTuple: cada cor é lida pelo nome (COR.s1, COR.critico...).
# s1, s2, s3... = séries de dados; serio (laranja) = 65 mmHg; critico
# (vermelho) = 55 mmHg; bom (verde) = janela de dose; o resto = texto, grade,
# eixos e fundo.
const COR = (s1 = "#2a78d6", s2 = "#eb6834", s3 = "#1baf7a", s4 = "#eda100", s7 = "#6250d6",
             critico = "#d03b3b", serio = "#ec835a", bom = "#0ca30c",
             tinta = "#0b0b0b", tinta2 = "#52514e", muda = "#898781",
             grade = "#e1e0d9", eixo = "#c3c2b7", fundo = "#fcfcfb")

# estilo padrão de todas as figuras: fontes (legenda grande, para não sobrar espaço em
# branco), cores, sem grade (grid = false), margens e 150 dpi
default(fontfamily = "sans-serif", titlefontsize = 13, guidefontsize = 12, tickfontsize = 11,
        legendfontsize = 12, linewidth = 2.5, background_color = COR.fundo,
        foreground_color_axis = COR.eixo, foreground_color_border = COR.eixo,
        foreground_color_text = COR.tinta2, foreground_color_guide = COR.tinta2,
        foreground_color_title = COR.tinta, foreground_color_legend = nothing,
        background_color_legend = nothing, grid = false, framestyle = :axes, titlelocation = :left, dpi = 150,
        left_margin = 3Plots.mm, bottom_margin = 3Plots.mm)

# resultados/ do projeto; se a pasta não aceitar escrita (pacote instalado pelo
# Pkg), usa ~/PBL6_resultados
# _pasta_resultados(): testa se dá para gravar em resultados/ criando e
# apagando um arquivo de teste. try/catch: se der erro, usa o plano B.
function _pasta_resultados()
    d = normpath(joinpath(@__DIR__, "..", "resultados"))
    try
        mkpath(d)
        teste = joinpath(d, ".teste_escrita")
        touch(teste); rm(teste)
        return d
    catch
        # plano B: pasta PBL6_resultados dentro da pasta do usuário
        d = joinpath(homedir(), "PBL6_resultados")
        mkpath(d)
        return d
    end
end
# pasta escolhida, calculada uma vez quando este arquivo é incluído
const DIR_RESULTADOS = _pasta_resultados()
# saida("x.png") devolve o caminho completo de x.png nessa pasta
saida(nome) = joinpath(DIR_RESULTADOS, nome)

# minutos(r, t): tempo contado a partir do início do bolus (r.t_inducao, min)
minutos(r, t) = t .- r.t_inducao
# segundos(r, t): o mesmo tempo em segundos, usado no eixo x de todas as figuras
segundos(r, t) = 60 .* minutos(r, t)
# rótulo do eixo de tempo das figuras da indução e marcas a cada 120 s (2 min)
const EIXO_T = "Tempo após o início do bolus (s)"
const MARCAS_T = 0:120:1200

# limiares!(p, xmax): desenha no gráfico p as linhas tracejadas de 65 e 55
# mmHg, com o rótulo perto de x = xmax (o de 65 acima da linha, o de 55 abaixo,
# para não se encostarem). O `!` indica que a função altera p.
"""Linhas de limiar de PAM (65 = sério, 55 = crítico) com rótulo direto."""
function limiares!(p, xmax; limiares = (65.0, 55.0))
    # zip junta cada limiar com sua cor: 65 → laranja (sério), 55 → vermelho
    for (l, c) in zip(limiares, (COR.serio, COR.critico))
        hline!(p, [l]; color = c, ls = :dash, lw = 1.2, label = "")
        annotate!(p, xmax, l == maximum(limiares) ? l + 3 : l - 3, text(@sprintf("%.0f mmHg", l), 10, COR.tinta2, :right))
    end
    return p
end

# fig_inducao(r; titulo): figura principal da previsão, com 6 painéis.
# r = resultado de simular_inducao. A = Cp e Ce do propofol; B = BIS com a
# faixa-alvo 40–60; C = faixa PAS–PAD e PAM; D = FC; E = débito cardíaco;
# F = variação (%) de Rp, VS e FC em relação ao basal.
# Usada por 02_prever_paciente.jl (<nome>_inducao.png).
"""Painel 3×2 da indução: concentrações, BIS, PA, FC, DC e mecanismos."""
function fig_inducao(r; titulo = "")
    # b = valores por batimento; k = PK/PD (grade de 1 s); pv = paciente virtual
    b, k, pv = r.bat, r.pkpd, r.pv
    # eixos x em segundos após o início do bolus; xl = limites do eixo
    x = segundos(r, k.t)
    xb = segundos(r, b.t_min)
    xl = (minimum(xb), maximum(xb))
    # teto do eixo y do painel A: 1,4 × Ce máximo (no mínimo 6 µg/ml)
    ymax = max(1.4 * maximum(k.Ce), 6.0)
    # A: Cp é cortada no teto (min.), porque no bolus ela sobe muito acima de Ce
    pA = plot(x, min.(k.Cp, ymax); color = COR.s1, label = "Cp (plasma)", title = "A. Propofol — NeoPKPD",
              ylims = (0, ymax), xlims = xl, legend = :topright, xlabel = EIXO_T, ylabel = "Concentração (µg/ml)", xticks = MARCAS_T)
    plot!(pA, x, k.Ce; color = COR.s2, label = "Ce (sítio efetor)")
    # B: BIS de Eleveld, com a faixa-alvo 40–60 sombreada
    pB = plot(x, k.bis; color = COR.s1, label = "", title = "B. BIS (Eleveld)", ylims = (0, 100), xlims = xl,
              xlabel = EIXO_T, ylabel = "BIS (0–100)", xticks = MARCAS_T)
    hspan!(pB, [40, 60]; color = COR.muda, alpha = 0.12, label = "")
    annotate!(pB, xl[2], 50, text("alvo 40–60", 10, COR.tinta2, :right))
    # C: faixa entre PAS e PAD (fillrange) e linha da PAM, com os limiares
    pC = plot(xb, b.PAS; fillrange = b.PAD, fillalpha = 0.15, fillcolor = COR.muda, lw = 0, label = "PAS–PAD",
              title = "C. Pressão arterial — Windkessel 4 elementos", xlims = xl, legend = :topright,
              xlabel = EIXO_T, ylabel = "Pressão arterial (mmHg)", xticks = MARCAS_T,
              ylims = (min(40, minimum(b.PAD) - 5), maximum(b.PAS) + 25))
    plot!(pC, xb, b.PAM; color = COR.s1, label = "PAM")
    limiares!(pC, xl[2])
    # D e E: frequência cardíaca e débito cardíaco
    pD = plot(xb, b.FC; color = COR.s1, label = "", title = "D. Frequência cardíaca", xlims = xl,
              ylims = (minimum(b.FC) - 5, maximum(b.FC) + 5), xlabel = EIXO_T, ylabel = "Frequência cardíaca (bpm)", xticks = MARCAS_T)
    pE = plot(xb, b.DC; color = COR.s1, label = "", title = "E. Débito cardíaco", xlims = xl,
              ylims = (0, maximum(b.DC) * 1.15), xlabel = EIXO_T, ylabel = "Débito cardíaco (L/min)", xticks = MARCAS_T)
    # F: variação percentual de Rp, VS e FC em relação ao valor basal
    vs0, Rp0, fc0 = pv.obtido.VS, pv.θ[1], pv.paciente.fc
    pF = plot(xb, 100 .* (b.Rp ./ Rp0 .- 1); color = COR.s1, label = "Rp (vasodilatação)",
              title = "F. Efeitos do propofol", xlims = xl, legend = :bottomright,
              xlabel = EIXO_T, ylabel = "Variação em relação ao basal (%)", xticks = MARCAS_T)
    plot!(pF, xb, 100 .* (b.VS ./ vs0 .- 1); color = COR.s2, label = "volume sistólico")
    plot!(pF, xb, 100 .* (b.FC ./ fc0 .- 1); color = COR.s3, label = "frequência cardíaca")
    hline!(pF, [0]; color = COR.eixo, lw = 1, label = "")
    # junta os 6 painéis em 3 linhas × 2 colunas
    return plot(pA, pB, pC, pD, pE, pF; layout = (3, 2), size = (1200, 1150), left_margin = 5Plots.mm, plot_title = titulo, plot_titlefontsize = 14)
end

# fig_onda(r, mt): ondas de pressão e de fluxo aórticos. Compara 3 batimentos
# antes do bolus (basal) com 3 batimentos em volta da PAM mínima (nadir).
# Exige simular_inducao(...; forma_de_onda = true).
# Usada por 02_prever_paciente.jl (<nome>_onda.png).
"""Pressão e fluxo aórticos (Windkessel 4 elementos): batimentos basais × no nadir."""
function fig_onda(r, mt)
    o = r.onda
    # t0 = início do bolus e tn = instante do nadir, em segundos
    t0 = 60 * r.t_inducao
    tn = 60 * (r.t_inducao + mt.t_nadir)
    # jb = índices dos 3 batimentos antes do bolus (60/FC = duração de um
    # batimento, em s); jn = de 1,5 batimento antes a 1,5 depois do nadir
    jb = findall(t -> t0 - 3 * 60 / mt.fc_basal <= t < t0, o.t_s)
    jn = findall(t -> tn - 1.5 * 60 / mt.fc_no_nadir <= t < tn + 1.5 * 60 / mt.fc_no_nadir, o.t_s)
    # folga no alto dos dois painéis para a legenda não cobrir as curvas
    ymax = maximum(o.pressao[jb]) + 30
    # A: pressão; os dois trechos começam em 0 s para ficarem sobrepostos
    pA = plot(o.t_s[jb] .- o.t_s[jb[1]], o.pressao[jb]; color = COR.s1, label = "basal", title = "A. Pressão arterial",
              xlabel = "Tempo (s)", ylabel = "Pressão arterial (mmHg)", ylims = (max(0, min(minimum(o.pressao[jb]), minimum(o.pressao[jn])) - 5), ymax), legend = :topright)
    plot!(pA, o.t_s[jn] .- o.t_s[jn[1]], o.pressao[jn]; color = COR.s2, label = @sprintf("nadir (%.0f s)", 60 * mt.t_nadir))
    # B: fluxo aórtico imposto (forma de Stergiopulos escalada por VS e FC)
    pB = plot(o.t_s[jb] .- o.t_s[jb[1]], o.fluxo[jb]; color = COR.s1, label = "basal", title = "B. Fluxo aórtico (Stergiopulos 1999)",
              xlabel = "Tempo (s)", ylabel = "Fluxo aórtico (ml/s)", legend = :topright,
              ylims = (minimum(o.fluxo[jb]) - 20, 1.3 * maximum(o.fluxo[jb])))
    plot!(pB, o.t_s[jn] .- o.t_s[jn[1]], o.fluxo[jn]; color = COR.s2, label = "nadir")
    return plot(pA, pB; layout = (1, 2), size = (1200, 440), bottom_margin = 6Plots.mm, left_margin = 6Plots.mm)
end

# fig_incerteza(inc, r): EXPERIMENTAL, NÃO VALIDADA. Desenha a faixa de 90%
# (percentis 5–95) da PAM calculada por previsao_com_incerteza (src/incerteza.jl).
# Não é usada pelos scripts atuais nem pela plataforma, porque as
# variabilidades ainda não foram validadas (docs/validacao.md, seção 6).
"""Curva de PAM prevista com faixa de incerteza (5–95%)."""
function fig_incerteza(inc, r; titulo = "")
    # inc.t está em minutos após o bolus; ts = o mesmo tempo em segundos
    ts = 60 .* inc.t
    p = plot(ts, inc.p95; fillrange = inc.p05, fillalpha = 0.18, fillcolor = COR.s1, lw = 0,
             label = "faixa 90% (PK/PD sorteadas)", title = titulo, xlabel = EIXO_T,
             ylabel = "PAM (mmHg)", xticks = MARCAS_T, legend = :topright, size = (850, 450),
             ylims = (min(40, minimum(inc.p05) - 5), maximum(inc.p95) + 15))
    plot!(p, ts, inc.mediana; color = COR.s1, label = "mediana")
    plot!(p, segundos(r, r.bat.t_min), r.bat.PAM; color = COR.s2, lw = 1.2, label = "paciente típico (parâmetros médios)")
    limiares!(p, maximum(ts))
    return p
end

# fig_doses(tab, jan): resultado da varredura de doses (03_janela_de_dose.jl).
# tab = métricas por dose; jan = janela de dose. A = PAM mínima; B = BIS
# mínimo, com linha no alvo; C = segundos com PAM < 65 e < 55.
# Se a janela existe, ela aparece como faixa verde nos três painéis.
"""Varredura de dose: nadir de PAM e BIS mínimo × dose."""
function fig_doses(tab, jan; titulo = "")
    # d = doses testadas (mg/kg)
    d = [t.dose for t in tab]
    pA = plot(d, [t.pam_nadir for t in tab]; color = COR.s1, marker = :circle, ms = 4, msw = 0, label = "",
              title = "A. PAM nadir prevista", xlabel = "Dose (mg/kg)", ylabel = "PAM nadir (mmHg)")
    limiares!(pA, maximum(d))
    pB = plot(d, [t.bis_min for t in tab]; color = COR.s1, marker = :circle, ms = 4, msw = 0, label = "",
              title = "B. BIS mínimo previsto", xlabel = "Dose (mg/kg)", ylabel = "BIS mínimo (0–100)", ylims = (0, 100))
    hline!(pB, [jan.alvo_bis]; color = COR.muda, ls = :dash, lw = 1.2, label = "")
    # tempo_abaixo está em minutos; × 60 para segundos
    pC = plot(d, [60 * t.tempo_abaixo[65.0] for t in tab]; color = COR.s1, marker = :circle, ms = 4, msw = 0,
              label = "PAM < 65", title = "C. Tempo em hipotensão", xlabel = "Dose (mg/kg)",
              ylabel = "Tempo com PAM abaixo do limiar (s)\nnos 20 min após o bolus", legend = :topleft)
    plot!(pC, d, [60 * t.tempo_abaixo[55.0] for t in tab]; color = COR.s2, marker = :circle, ms = 4, msw = 0,
          label = "PAM < 55")
    # destaca a janela (vspan! = faixa vertical) e escreve "janela" no painel A
    if jan.existe
        for p in (pA, pB, pC)
            vspan!(p, [jan.dose_min, jan.dose_max]; color = COR.bom, alpha = 0.10, label = "")
        end
        annotate!(pA, (jan.dose_min + jan.dose_max) / 2, maximum(t.pam_nadir for t in tab),
                  text("janela", 11, COR.tinta2, :center))
    end
    # junta os 3 painéis lado a lado
    return plot(pA, pB, pC; layout = (1, 3), size = (1300, 460), plot_title = titulo,
                plot_titlefontsize = 14, bottom_margin = 6Plots.mm, left_margin = 5Plots.mm)
end

