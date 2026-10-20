# =============================================================================
# graficos.jl — figuras do projeto (Plots.jl + GR)
# Incluído pelos scripts: include(joinpath(@__DIR__, "graficos.jl"))
# No REPL: include(joinpath(pkgdir(PropofolWindkessel), "scripts", "graficos.jl"))
# =============================================================================
# rodando como script, o GR só grava arquivos; no REPL, abre janelas
isinteractive() || (ENV["GKSwstype"] = get(ENV, "GKSwstype", "100"))
using Plots
using Printf
gr()

# Paleta: slots categóricos em ordem fixa + cores de status reservadas
const COR = (s1 = "#2a78d6", s2 = "#eb6834", s3 = "#1baf7a", s4 = "#eda100", s7 = "#6250d6",
             critico = "#d03b3b", serio = "#ec835a", bom = "#0ca30c",
             tinta = "#0b0b0b", tinta2 = "#52514e", muda = "#898781",
             grade = "#e1e0d9", eixo = "#c3c2b7", fundo = "#fcfcfb")

default(fontfamily = "sans-serif", titlefontsize = 10, guidefontsize = 9, tickfontsize = 8,
        legendfontsize = 8, linewidth = 2, background_color = COR.fundo,
        foreground_color_axis = COR.eixo, foreground_color_border = COR.eixo,
        foreground_color_text = COR.tinta2, foreground_color_guide = COR.tinta2,
        foreground_color_title = COR.tinta, foreground_color_legend = nothing,
        background_color_legend = nothing, gridcolor = COR.grade, gridalpha = 1.0,
        gridlinewidth = 0.5, framestyle = :axes, titlelocation = :left, dpi = 150,
        left_margin = 3Plots.mm, bottom_margin = 3Plots.mm)

# resultados/ do projeto; se a pasta não aceitar escrita (pacote instalado pelo
# Pkg), usa ~/PBL6_resultados
function _pasta_resultados()
    d = normpath(joinpath(@__DIR__, "..", "resultados"))
    try
        mkpath(d)
        teste = joinpath(d, ".teste_escrita")
        touch(teste); rm(teste)
        return d
    catch
        d = joinpath(homedir(), "PBL6_resultados")
        mkpath(d)
        return d
    end
end
const DIR_RESULTADOS = _pasta_resultados()
saida(nome) = joinpath(DIR_RESULTADOS, nome)

minutos(r, t) = (t .- r.t_inducao) ./ 60

"""Linhas de limiar de PAM (65 = sério, 55 = crítico) com rótulo direto."""
function limiares!(p, xmax; limiares = (65.0, 55.0))
    for (l, c) in zip(limiares, (COR.serio, COR.critico))
        hline!(p, [l]; color = c, ls = :dash, lw = 1.2, label = "")
        annotate!(p, xmax, l + 2.5, text(@sprintf("%.0f mmHg", l), 7, COR.tinta2, :right))
    end
    return p
end

"""Painel 3×2 da indução: concentrações, BIS, PA, FC, DC e mecanismos."""
function fig_inducao(r, pv; titulo = "")
    b = r.bat
    x = minutos(r, r.t)
    xb = minutos(r, b.t)
    xl = (minimum(xb), maximum(xb))
    s = r.sinais
    # A — concentrações
    ymax = max(1.4 * maximum(s.Ceb), 6.0)
    pA = plot(x, min.(s.Cp, ymax); color = COR.s1, label = "Cp (plasma)", title = "A. Propofol (µg/ml)",
              ylims = (0, ymax), xlims = xl, legend = :topright)
    plot!(pA, x, s.Ceb; color = COR.s2, label = "Ce cérebro (BIS)")
    plot!(pA, x, s.Ceh; color = COR.s3, label = "Ce vascular (PA)")
    # B — BIS
    pB = plot(x, s.BIS; color = COR.s1, label = "", title = "B. BIS (hipnose)", ylims = (0, 100), xlims = xl)
    hspan!(pB, [40, 60]; color = COR.muda, alpha = 0.12, label = "")
    annotate!(pB, xl[2], 50, text("alvo 40–60", 7, COR.tinta2, :right))
    # C — pressão arterial
    pC = plot(xb, b.PAS; fillrange = b.PAD, fillalpha = 0.15, fillcolor = COR.s1, lw = 0, label = "PAS–PAD",
              title = "C. Pressão arterial (mmHg)", xlims = xl, legend = :topright,
              ylims = (min(40, minimum(b.PAD) - 5), maximum(b.PAS) + 25))
    plot!(pC, xb, b.PAM; color = COR.s1, label = "PAM")
    limiares!(pC, xl[2])
    # D — FC
    pD = plot(xb, b.FC; color = COR.s1, label = "", title = "D. Frequência cardíaca (bpm)", xlims = xl,
              ylims = (minimum(b.FC) - 8, maximum(b.FC) + 8))
    # E — DC
    pE = plot(xb, b.DC; color = COR.s1, label = "", title = "E. Débito cardíaco (L/min)", xlims = xl,
              ylims = (0, maximum(b.DC) * 1.15), xlabel = "minutos após o início do bolus")
    # F — mecanismos (% do basal)
    cv = pv.cv
    Vest(Vv0) = cv.Vtotal - cv.Va0 - Vv0 - cv.Vd
    pF = plot(x, 100 .* (s.Rs ./ cv.Rs .- 1); color = COR.s1, label = "TPR (vasodilatação)",
              title = "F. Mecanismos (% do basal)", xlims = xl, legend = :bottomright,
              xlabel = "minutos após o início do bolus")
    plot!(pF, x, 100 .* (s.Emax ./ cv.Emax .- 1); color = COR.s2, label = "Emax = 1/C_S (contratilidade)")
    plot!(pF, x, 100 .* (Vest.(s.Vv0) ./ Vest(cv.Vv0) .- 1); color = COR.s3, label = "volume estressado (venodilatação)")
    hline!(pF, [0]; color = COR.eixo, lw = 1, label = "")
    f = plot(pA, pB, pC, pD, pE, pF; layout = (3, 2), size = (1100, 950),
             plot_title = titulo, plot_titlefontsize = 11)
    return f
end

"""Ondas de pressão e laço pressão–volume: basal × nadir."""
function fig_ondas(r, mt)
    s = r.sinais
    t0 = r.t_inducao
    tn = t0 + 60 * mt.t_nadir
    T = 60 / mt.fc_basal
    jb = findall(t -> t0 - 3T <= t < t0, r.t)
    Tn = 60 / mt.fc_no_nadir
    jn = findall(t -> tn - 1.5Tn <= t < tn + 1.5Tn, r.t)
    ymax = maximum(s.Plv[jb]) + 10
    pA = plot(r.t[jb] .- r.t[jb[1]], s.Pao[jb]; color = COR.s1, label = "aorta", title = "A. Basal",
              ylims = (0, ymax), xlabel = "s", ylabel = "mmHg", legend = :topright)
    plot!(pA, r.t[jb] .- r.t[jb[1]], s.Plv[jb]; color = COR.s2, label = "ventrículo")
    pB = plot(r.t[jn] .- r.t[jn[1]], s.Pao[jn]; color = COR.s1, label = "aorta",
              title = @sprintf("B. No nadir (%.1f min)", mt.t_nadir), ylims = (0, ymax), xlabel = "s", legend = :topright)
    plot!(pB, r.t[jn] .- r.t[jn[1]], s.Plv[jn]; color = COR.s2, label = "ventrículo")
    pC = plot(s.Vlv[jb], s.Plv[jb]; color = COR.s1, label = "basal", title = "C. Laço pressão–volume",
              xlabel = "volume do ventrículo (ml)", ylabel = "mmHg", legend = :topleft)
    plot!(pC, s.Vlv[jn], s.Plv[jn]; color = COR.s2, label = "nadir")
    return plot(pA, pB, pC; layout = (1, 3), size = (1150, 380), bottom_margin = 6Plots.mm,
                left_margin = 6Plots.mm)
end

"""Diagrama de Guyton (curso, Seção 16): basal × nadir do propofol."""
function fig_guyton(pv, r, mt)
    cv = pv.cv
    s = r.sinais
    i = argmin(abs.(r.t .- (r.t_inducao + 60 * mt.t_nadir)))
    g0 = estado_guyton(cv; Pa = mt.pam_basal)
    g1 = estado_guyton(cv; FC = s.FC[i], Rs = s.Rs[i], Emax = s.Emax[i], Vv0 = s.Vv0[i], Pa = mt.pam_nadir)
    rap = range(-2, max(g0.Pms, g1.Pms) + 1; length = 200)
    p = plot(; title = "Ponto de funcionamento: retorno venoso × débito cardíaco",
             xlabel = "pressão de enchimento — RAP (mmHg)", ylabel = "fluxo (L/min)", legend = :topright,
             size = (800, 520), ylims = (0, 1.6 * max(g0.Q, g1.Q)))
    for (k, (g, c, nome)) in enumerate(((g0, COR.s1, "basal"), (g1, COR.s2, "propofol (nadir)")))
        plot!(p, rap, retorno_venoso.(rap, g.Pms, g.RVR); color = c, label = "retorno venoso — " * nome)
        plot!(p, rap, debito_cardiaco_curso.(rap, g.FC, g.CD, g.CS, g.Pa); color = c, lw = 1.2,
              label = "curva cardíaca — " * nome)
        scatter!(p, [g.RAP], [g.Q]; color = c, ms = 6, msw = 2, msc = COR.fundo, label = "")
        annotate!(p, g.RAP + 0.4, g.Q, text(@sprintf("%.1f L/min", g.Q), 8, COR.tinta, :left))
        annotate!(p, g.Pms, (0.06 + 0.07 * k) * max(g0.Q, g1.Q), text(@sprintf("Pms %.1f", g.Pms), 7, COR.tinta2, :center))
    end
    return p, g0, g1
end

"""Varredura de dose: nadir de PAM e BIS mínimo × dose."""
function fig_doses(tab, jan; titulo = "")
    d = [t.dose for t in tab]
    pA = plot(d, [t.pam_nadir for t in tab]; color = COR.s1, marker = :circle, ms = 4, msw = 0, label = "",
              title = "A. PAM nadir prevista (mmHg)", xlabel = "dose (mg/kg)")
    limiares!(pA, maximum(d))
    pB = plot(d, [t.bis_min for t in tab]; color = COR.s1, marker = :circle, ms = 4, msw = 0, label = "",
              title = "B. BIS mínimo previsto", xlabel = "dose (mg/kg)", ylims = (0, 100))
    hline!(pB, [jan.alvo_bis]; color = COR.muda, ls = :dash, lw = 1.2, label = "")
    pC = plot(d, [t.tempo_abaixo[65.0] for t in tab]; color = COR.s1, marker = :circle, ms = 4, msw = 0,
              label = "PAM < 65", title = "C. Tempo em hipotensão (min, 20 min)", xlabel = "dose (mg/kg)",
              legend = :topleft)
    plot!(pC, d, [t.tempo_abaixo[55.0] for t in tab]; color = COR.s2, marker = :circle, ms = 4, msw = 0,
          label = "PAM < 55")
    if jan.existe
        for p in (pA, pB, pC)
            vspan!(p, [jan.dose_min, jan.dose_max]; color = COR.bom, alpha = 0.10, label = "")
        end
        annotate!(pA, (jan.dose_min + jan.dose_max) / 2, maximum(t.pam_nadir for t in tab),
                  text("janela", 8, COR.tinta2, :center))
    end
    return plot(pA, pB, pC; layout = (1, 3), size = (1150, 380), plot_title = titulo,
                plot_titlefontsize = 11, bottom_margin = 6Plots.mm)
end

"""Curva de PAM prevista com faixa de incerteza farmacológica (5–95%)."""
function fig_incerteza(inc, r; titulo = "")
    p = plot(inc.t, inc.p95; fillrange = inc.p05, fillalpha = 0.18, fillcolor = COR.s1, lw = 0,
             label = "faixa 90% (PK/PD sorteadas)", title = titulo, xlabel = "minutos após o início do bolus",
             ylabel = "PAM (mmHg)", legend = :topright, size = (850, 450),
             ylims = (min(40, minimum(inc.p05) - 5), maximum(inc.p95) + 15))
    plot!(p, inc.t, inc.mediana; color = COR.s1, label = "mediana")
    plot!(p, minutos(r, r.bat.t), r.bat.PAM; color = COR.s2, lw = 1.2, label = "paciente típico (parâmetros médios)")
    limiares!(p, maximum(inc.t))
    return p
end
