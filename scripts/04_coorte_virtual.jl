# =============================================================================
# 04_coorte_virtual.jl — ensaio clínico virtual (Resumo, Seção 8.7)
#
# Gera N adultos virtuais (com ASA, PVI e IVCI, como no modelo do grupo),
# calibra cada um com os próprios dados pré-indução e simula a indução
# (bolus em 30 s + manutenção de 100 µg/kg/min, esquema do grupo) em:
#   1. dose fixa de 2 mg/kg para todos;
#   2. dose ajustada por idade (2,0 / 1,5 / 1,2 mg/kg);
#   3. dose ajustada por idade, só bolus (sem manutenção).
#
# EXPERIMENTAL — NÃO VALIDADO: as distribuições da coorte e as variabilidades
# são suposições (docs/validacao.md, pendências). Fica fora de rodar_tudo.jl e
# os resultados não são usados no README nem na plataforma.
#
# Uso: julia --project=. -t auto scripts/04_coorte_virtual.jl 500
# =============================================================================
using PropofolWindkessel
using Statistics
include(joinpath(@__DIR__, "graficos.jl"))

n = isempty(ARGS) ? 500 : parse(Int, ARGS[1])
println("Coorte virtual: $n adultos | threads: $(Threads.nthreads())")
coorte = gerar_coorte(n)
t = @elapsed pvs = calibrar_coorte(coorte)
@printf("Calibração: %.0f s; %d/%d pacientes calibrados\n", t, count(pv -> pv.convergiu, pvs), n)

cenarios = [("2 mg/kg para todos", (; dose_fixa = 2.0)),
            ("ajustada por idade", NamedTuple()),
            ("ajustada, só bolus", (; manutencao_ugkgmin = 0.0))]
resultados = Dict{String,Any}()
println("\n protocolo                               PAM<65   PAM<55   queda>30%   queda mediana (IIQ)   nadir mediano")
for (nome, kw) in cenarios
    res = filter(r -> r.calibrou, simular_coorte(coorte, pvs; kw...))
    resultados[nome] = res
    s = resumo_coorte(res)
    @printf(" %-38s %5.1f%%   %5.1f%%   %5.1f%%      %4.1f%% (%4.1f–%4.1f)      %.1f min\n", nome, s.inc_pam65,
            s.inc_pam55, s.inc_queda30, s.queda_mediana, s.queda_p25, s.queda_p75, s.t_nadir_mediano)
end
println(" referência clínica (problema):                   17–24%")

faixas = ((18, 40), (40, 55), (55, 70), (70, 91))
rotulos = ["18–39", "40–54", "55–69", "≥70"]
incidencia(res, f, campo) = 100 * mean(getfield(r, campo) > 0 for r in res if f[1] <= r.idade < f[2])
inc55 = Dict{String,Vector{Float64}}()
for nome in ("2 mg/kg para todos", "ajustada por idade")
    res = resultados[nome]
    println("\nPor faixa etária — $nome:")
    inc55[nome] = Float64[]
    for (f, rot) in zip(faixas, rotulos)
        sel = filter(r -> f[1] <= r.idade < f[2], res)
        push!(inc55[nome], incidencia(res, f, :min_abaixo_55))
        @printf("  %-6s n=%3d  dose média %.2f mg/kg  queda mediana %4.1f%%  PAM<65 %5.1f%%  PAM<55 %5.1f%%\n", rot,
                length(sel), mean(r.dose_mgkg for r in sel), median(r.queda_pct for r in sel),
                incidencia(res, f, :min_abaixo_65), inc55[nome][end])
    end
end
res = resultados["ajustada por idade"]
println("\nPor ASA (dose ajustada por idade):")
for a in 1:4
    sel = filter(r -> r.asa == a, res)
    isempty(sel) && continue
    @printf("  ASA %d  n=%3d  sensibilidade média %.2f  queda mediana %4.1f%%  PAM<55 %5.1f%%\n", a, length(sel),
            mean(r.s for r in sel), median(r.queda_pct for r in sel), 100 * mean(r.min_abaixo_55 > 0 for r in sel))
end

# --- figura -------------------------------------------------------------------
fixa = resultados["2 mg/kg para todos"]
pA = histogram([r.pam_nadir for r in res]; bins = 30:2.5:110, color = COR.s1, lw = 0, label = "",
               title = "A. PAM nadir (dose ajustada por idade)", xlabel = "mmHg", ylabel = "pacientes")
for (l, c) in zip((65, 55), (COR.serio, COR.critico))
    vline!(pA, [l]; color = c, ls = :dash, lw = 1.2, label = "")
    annotate!(pA, l + 1, 1, text("$l", 7, COR.tinta2, :left, :bottom))
end
pB = scatter([r.idade for r in fixa], [r.queda_pct for r in fixa]; color = COR.s2, ms = 3, msw = 0, alpha = 0.6,
             label = "2 mg/kg para todos", title = "B. Queda da PAM × idade", xlabel = "idade (anos)",
             ylabel = "queda (%)", legend = :topleft)
scatter!(pB, [r.idade for r in res], [r.queda_pct for r in res]; color = COR.s1, ms = 3, msw = 0, alpha = 0.6,
         label = "ajustada por idade")
x = 1:length(rotulos)
pC = bar(x .- 0.2, inc55["2 mg/kg para todos"]; bar_width = 0.38, color = COR.s2, lw = 0, label = "2 mg/kg para todos",
         title = "C. PAM < 55 mmHg por faixa etária (%)", xticks = (x, rotulos), ylims = (0, 70), legend = :topleft)
bar!(pC, x .+ 0.2, inc55["ajustada por idade"]; bar_width = 0.38, color = COR.s1, lw = 0, label = "ajustada por idade")
hspan!(pC, [17, 24]; color = COR.muda, alpha = 0.18, lw = 0, label = "17–24% (literatura)")
savefig(plot(pA, pB, pC; layout = (1, 3), size = (1200, 400), bottom_margin = 6Plots.mm, left_margin = 5Plots.mm),
        saida("coorte_virtual.png"))

open(saida("coorte_virtual.csv"), "w") do io
    println(io, "idade,sexo,peso,asa,pvi,ivci,sensibilidade,pam_basal,dose_mgkg,pam_nadir,queda_pct,t_nadir,min_abaixo_65,min_abaixo_55,bis_min,risco")
    for r in res
        @printf(io, "%.1f,%s,%.1f,%d,%.1f,%.1f,%.3f,%.1f,%.2f,%.1f,%.1f,%.2f,%.2f,%.2f,%.0f,%s\n", r.idade, r.sexo, r.peso,
                r.asa, r.pvi, r.ivci, r.s, r.pam_basal, r.dose_mgkg, r.pam_nadir, r.queda_pct, r.t_nadir,
                r.min_abaixo_65, r.min_abaixo_55, r.bis_min, r.risco)
    end
end
println("\nFigura e CSV em: ", normpath(DIR_RESULTADOS))
