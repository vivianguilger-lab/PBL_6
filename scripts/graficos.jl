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

minutos(r, t) = t .- r.t_inducao

"""Linhas de limiar de PAM (65 = sério, 55 = crítico) com rótulo direto."""
function limiares!(p, xmax; limiares = (65.0, 55.0))
    for (l, c) in zip(limiares, (COR.serio, COR.critico))
        hline!(p, [l]; color = c, ls = :dash, lw = 1.2, label = "")
        annotate!(p, xmax, l + 2.5, text(@sprintf("%.0f mmHg", l), 7, COR.tinta2, :right))
    end
    return p
end

"""Painel 3×2 da indução: concentrações, BIS, PA, FC, DC e mecanismos."""
function fig_inducao(r; titulo = "")
    b, k, pv = r.bat, r.pkpd, r.pv
    x = minutos(r, k.t)
    xb = minutos(r, b.t_min)
    xl = (minimum(xb), maximum(xb))
    ymax = max(1.4 * maximum(k.Ce), 6.0)
    pA = plot(x, min.(k.Cp, ymax); color = COR.s1, label = "Cp (plasma)", title = "A. Propofol (µg/ml) — NeoPKPD",
              ylims = (0, ymax), xlims = xl, legend = :topright)
    plot!(pA, x, k.Ce; color = COR.s2, label = "Ce (sítio efetor)")
    pB = plot(x, k.bis; color = COR.s1, label = "", title = "B. BIS (Eleveld)", ylims = (0, 100), xlims = xl)
    hspan!(pB, [40, 60]; color = COR.muda, alpha = 0.12, label = "")
    annotate!(pB, xl[2], 50, text("alvo 40–60", 7, COR.tinta2, :right))
    pC = plot(xb, b.PAS; fillrange = b.PAD, fillalpha = 0.15, fillcolor = COR.muda, lw = 0, label = "PAS–PAD",
              title = "C. Pressão arterial (mmHg) — Windkessel 4 elementos", xlims = xl, legend = :topright,
              ylims = (min(40, minimum(b.PAD) - 5), maximum(b.PAS) + 25))
    plot!(pC, xb, b.PAM; color = COR.s1, label = "PAM")
    limiares!(pC, xl[2])
    pD = plot(xb, b.FC; color = COR.s1, label = "", title = "D. Frequência cardíaca (bpm)", xlims = xl,
              ylims = (minimum(b.FC) - 5, maximum(b.FC) + 5))
    pE = plot(xb, b.DC; color = COR.s1, label = "", title = "E. Débito cardíaco (L/min)", xlims = xl,
              ylims = (0, maximum(b.DC) * 1.15), xlabel = "minutos após o início do bolus")
    vs0, Rp0, fc0 = pv.obtido.VS, pv.θ[1], pv.paciente.fc
    pF = plot(xb, 100 .* (b.Rp ./ Rp0 .- 1); color = COR.s1, label = "Rp (vasodilatação)",
              title = "F. Efeitos do propofol (% do basal)", xlims = xl, legend = :bottomright,
              xlabel = "minutos após o início do bolus")
    plot!(pF, xb, 100 .* (b.VS ./ vs0 .- 1); color = COR.s2, label = "volume sistólico")
    plot!(pF, xb, 100 .* (b.FC ./ fc0 .- 1); color = COR.s3, label = "frequência cardíaca")
    hline!(pF, [0]; color = COR.eixo, lw = 1, label = "")
    return plot(pA, pB, pC, pD, pE, pF; layout = (3, 2), size = (1100, 950), plot_title = titulo, plot_titlefontsize = 11)
end

"""Pressão e fluxo aórticos (Windkessel 4 elementos): batimentos basais × no nadir."""
function fig_onda(r, mt)
    o = r.onda
    t0 = 60 * r.t_inducao
    tn = 60 * (r.t_inducao + mt.t_nadir)
    jb = findall(t -> t0 - 3 * 60 / mt.fc_basal <= t < t0, o.t_s)
    jn = findall(t -> tn - 1.5 * 60 / mt.fc_no_nadir <= t < tn + 1.5 * 60 / mt.fc_no_nadir, o.t_s)
    ymax = maximum(o.pressao[jb]) + 10
    pA = plot(o.t_s[jb] .- o.t_s[jb[1]], o.pressao[jb]; color = COR.s1, label = "basal", title = "A. Pressão arterial",
              xlabel = "s", ylabel = "mmHg", ylims = (max(0, min(minimum(o.pressao[jb]), minimum(o.pressao[jn])) - 5), ymax), legend = :topright)
    plot!(pA, o.t_s[jn] .- o.t_s[jn[1]], o.pressao[jn]; color = COR.s2, label = @sprintf("nadir (%.1f min)", mt.t_nadir))
    pB = plot(o.t_s[jb] .- o.t_s[jb[1]], o.fluxo[jb]; color = COR.s1, label = "basal", title = "B. Fluxo aórtico (Stergiopulos 1999)",
              xlabel = "s", ylabel = "ml/s", legend = :topright)
    plot!(pB, o.t_s[jn] .- o.t_s[jn[1]], o.fluxo[jn]; color = COR.s2, label = "nadir")
    return plot(pA, pB; layout = (1, 2), size = (1100, 380), bottom_margin = 6Plots.mm, left_margin = 6Plots.mm)
end

"""Curva de PAM prevista com faixa de incerteza (5–95%)."""
function fig_incerteza(inc, r; titulo = "")
    p = plot(inc.t, inc.p95; fillrange = inc.p05, fillalpha = 0.18, fillcolor = COR.s1, lw = 0,
             label = "faixa 90% (PK/PD sorteadas)", title = titulo, xlabel = "minutos após o início do bolus",
             ylabel = "PAM (mmHg)", legend = :topright, size = (850, 450),
             ylims = (min(40, minimum(inc.p05) - 5), maximum(inc.p95) + 15))
    plot!(p, inc.t, inc.mediana; color = COR.s1, label = "mediana")
    plot!(p, minutos(r, r.bat.t_min), r.bat.PAM; color = COR.s2, lw = 1.2, label = "paciente típico (parâmetros médios)")
    limiares!(p, maximum(inc.t))
    return p
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

