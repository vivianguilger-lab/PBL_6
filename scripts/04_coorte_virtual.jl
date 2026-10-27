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
#
# Como rodar (opcional): o número no fim do comando é o tamanho da coorte
# (sem número, 500). Com 500 pacientes leva alguns minutos; use -t auto.
# Saída na tela: incidência de PAM < 65 e < 55 por protocolo, por faixa etária
# e por ASA. Saída em resultados/: coorte_virtual.png e coorte_virtual.csv.
# Os sorteios usam semente fixa (Xoshiro(2026), em src/populacao.jl): o mesmo n
# dá sempre o mesmo resultado. Isso NÃO torna o resultado válido: idade, PA,
# ASA, PVI, IVCI e as variabilidades da PK/PD foram supostas, sem fonte.
# Não use estes números como evidência em trabalho ou apresentação.
# =============================================================================
# carrega o pacote, estatística (mean, median) e as funções de figura
using PropofolWindkessel
using Statistics
include(joinpath(@__DIR__, "graficos.jl"))

# n = tamanho da coorte; parse(Int, ...) converte o texto do comando em número
n = isempty(ARGS) ? 500 : parse(Int, ARGS[1])
println("Coorte virtual: $n adultos | threads: $(Threads.nthreads())")
# sorteia n adultos virtuais (dados pré-indução + variabilidade da PK/PD)
coorte = gerar_coorte(n)
# calibra o paciente virtual de cada um; @elapsed mede o tempo gasto (s)
t = @elapsed pvs = calibrar_coorte(coorte)
# quantos calibraram (C dentro de 0,2–6 ml/mmHg); `pv -> ...` é função anônima
@printf("Calibração: %.0f s; %d/%d pacientes calibrados\n", t, count(pv -> pv.convergiu, pvs), n)

# Três protocolos. Cada um é (nome, opções); as opções vão para simular_coorte.
# `(; dose_fixa = 2.0)` é uma NamedTuple (lista de opções com nome).
# NamedTuple() vazia = sem mudança: dose ajustada por idade, com manutenção.
cenarios = [("2 mg/kg para todos", (; dose_fixa = 2.0)),
            ("ajustada por idade", NamedTuple()),
            ("ajustada, só bolus", (; manutencao_ugkgmin = 0.0))]
# dicionário: nome do protocolo → resultados dos pacientes
resultados = Dict{String,Any}()
# Tabela por protocolo: % com PAM < 65, % com PAM < 55, % com queda > 30%,
# queda mediana com intervalo interquartil (IIQ) e tempo mediano até o nadir.
println("\n protocolo                               PAM<65   PAM<55   queda>30%   queda mediana (IIQ)   nadir mediano")
for (nome, kw) in cenarios
    # simula a coorte e descarta quem não calibrou; `kw...` repassa as opções
    res = filter(r -> r.calibrou, simular_coorte(coorte, pvs; kw...))
    resultados[nome] = res
    # incidências e medianas do protocolo
    s = resumo_coorte(res)
    @printf(" %-38s %5.1f%%   %5.1f%%   %5.1f%%      %4.1f%% (%4.1f–%4.1f)      %.1f min\n", nome, s.inc_pam65,
            s.inc_pam55, s.inc_queda30, s.queda_mediana, s.queda_p25, s.queda_p75, s.t_nadir_mediano)
end

# Faixas etárias (início, fim), com o fim excluído, e os rótulos do gráfico
faixas = ((18, 40), (40, 55), (55, 70), (70, 91))
rotulos = ["18–39", "40–54", "55–69", "≥70"]
# Função de uma linha: % dos pacientes da faixa f com `campo` > 0
# (ex.: minutos com PAM < 55 maior que zero = teve hipotensão).
# getfield(r, campo) lê o campo pelo nome (ex.: :min_abaixo_55).
incidencia(res, f, campo) = 100 * mean(getfield(r, campo) > 0 for r in res if f[1] <= r.idade < f[2])
# incidência de PAM < 55 por faixa, guardada para o painel C da figura
inc55 = Dict{String,Vector{Float64}}()
# Tabela por faixa etária para os dois protocolos com manutenção
for nome in ("2 mg/kg para todos", "ajustada por idade")
    res = resultados[nome]
    println("\nPor faixa etária — $nome:")
    inc55[nome] = Float64[]
    for (f, rot) in zip(faixas, rotulos)
        # pacientes desta faixa etária
        sel = filter(r -> f[1] <= r.idade < f[2], res)
        push!(inc55[nome], incidencia(res, f, :min_abaixo_55))
        @printf("  %-6s n=%3d  dose média %.2f mg/kg  queda mediana %4.1f%%  PAM<65 %5.1f%%  PAM<55 %5.1f%%\n", rot,
                length(sel), mean(r.dose_mgkg for r in sel), median(r.queda_pct for r in sel),
                incidencia(res, f, :min_abaixo_65), inc55[nome][end])
    end
end
# Tabela por ASA (protocolo de dose ajustada por idade).
# `isempty(sel) && continue` pula a classe ASA que não tem pacientes.
res = resultados["ajustada por idade"]
println("\nPor ASA (dose ajustada por idade):")
for a in 1:4
    sel = filter(r -> r.asa == a, res)
    isempty(sel) && continue
    @printf("  ASA %d  n=%3d  sensibilidade média %.2f  queda mediana %4.1f%%  PAM<55 %5.1f%%\n", a, length(sel),
            mean(r.s for r in sel), median(r.queda_pct for r in sel), 100 * mean(r.min_abaixo_55 > 0 for r in sel))
end

# --- figura -------------------------------------------------------------------
# A: histograma da PAM mínima (dose ajustada por idade), com linhas em 65 e 55.
# B: queda da PAM × idade nos dois protocolos (um ponto por paciente).
# C: % com PAM < 55 por faixa etária (barras lado a lado). Não há faixa de
#    referência da literatura: a antiga "17–24%" foi retirada por falta de fonte
#    (docs/validacao.md, seção 3, item 9).
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
savefig(plot(pA, pB, pC; layout = (1, 3), size = (1200, 400), bottom_margin = 6Plots.mm, left_margin = 5Plots.mm),
        saida("coorte_virtual.png"))

# CSV com uma linha por paciente (protocolo de dose ajustada por idade).
# Colunas: dados do paciente, sensibilidade s, PAM basal, dose, PAM mínima,
# queda (%), tempo até o nadir, minutos abaixo de 65 e 55, BIS mínimo e risco.
open(saida("coorte_virtual.csv"), "w") do io
    println(io, "idade,sexo,peso,asa,pvi,ivci,sensibilidade,pam_basal,dose_mgkg,pam_nadir,queda_pct,t_nadir,min_abaixo_65,min_abaixo_55,bis_min,risco")
    for r in res
        @printf(io, "%.1f,%s,%.1f,%d,%.1f,%.1f,%.3f,%.1f,%.2f,%.1f,%.1f,%.2f,%.2f,%.2f,%.0f,%s\n", r.idade, r.sexo, r.peso,
                r.asa, r.pvi, r.ivci, r.s, r.pam_basal, r.dose_mgkg, r.pam_nadir, r.queda_pct, r.t_nadir,
                r.min_abaixo_65, r.min_abaixo_55, r.bis_min, r.risco)
    end
end
# mostra a pasta das saídas
println("\nFigura e CSV em: ", normpath(DIR_RESULTADOS))
