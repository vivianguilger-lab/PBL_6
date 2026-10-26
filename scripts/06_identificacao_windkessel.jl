# =============================================================================
# 06_identificacao_windkessel.jl — identificação do Windkessel (windkessel-id)
#
# Como em exemplos/identificacao_invivo.jl do grupo e id_stergiopulos1999total.jl
# do windkessel-id: estima θ do Windkessel de 2, 3 e 4 elementos a partir de
# pressão e fluxo aórticos medidos (Newton com diferenciação automática).
#   A. dados humanos de Stergiopulos et al. (1999)
#   B. dados in vivo (porcos de 65 kg; windkessel-id v2.0, README)
# É o caminho para individualizar o Windkessel com a linha arterial do paciente.
#
# Uso: julia --project=. scripts/06_identificacao_windkessel.jl
# =============================================================================
using PropofolWindkessel, Statistics
include(joinpath(@__DIR__, "graficos.jl"))
h = 0.005

# --- A. Stergiopulos 1999 (mL/s) ------------------------------------------------
t, y, u = dados_stergiopulos1999()
θpub, _ = params_stergiopulos1999total()
resultadosA = Dict{String,Any}()
for (nome, G_, θ0) in (("2 elementos", G2, [mean(y) / mean(u), 1.0]),
                       ("3 elementos", G3, [mean(y) / mean(u), 1.0, 0.05]),
                       ("4 elementos", G,  [mean(y) / mean(u), 1.0, 0.05, 0.005]))
    x00 = solvex0(G_(θ0, h), u)
    phat, _ = identify(θ0, x00, u, y, h; x0id = false, x0solve = true, G_)
    θhat, _ = deal(phat, PropofolWindkessel.IdConditions(θ0, x00, ones(Int, length(θ0)), false, true, G_))
    yhat = simulate_wk(θhat, solvex0(G_(θhat, h), u), u, h; G_)
    rmse = sqrt(mean((y .- yhat) .^ 2))
    resultadosA[nome] = (θ = θhat, y = yhat, rmse = rmse)
    @printf("A. Stergiopulos — %s: θ̂ = %s  RMSE = %.2f mmHg\n", nome, string(round.(θhat; sigdigits = 3)), rmse)
end
println("   θ publicado (4 elementos) = ", round.(θpub; sigdigits = 3))

# --- B. dados in vivo (L/min) ------------------------------------------------------
P, Q = importPQ(joinpath(WINDKESSEL_ID_DIR, "invivo_data", "Ao_baseline1_Ao"), 4)
uB = vcat([q.vb_va for q in Q]...)
yB = vcat([p.vb_va for p in P]...)
θ0 = [mean(yB) / mean(uB), 1.0, 1.0, 1.0]
x00 = solvex0(G(θ0, h), uB)
phat, _ = identify(θ0, x00, uB, yB, h; x0id = false, x0solve = true)
θB, _ = deal(phat, PropofolWindkessel.IdConditions(θ0, x00, [1, 1, 1, 1], false, true))
yBhat = simulate_wk(θB, solvex0(G(θB, h), uB), uB, h)
@printf("B. in vivo — 4 elementos: θ̂ (L/min) = %s  RMSE = %.2f mmHg\n", string(round.(θB; sigdigits = 3)),
        sqrt(mean((yB .- yBhat) .^ 2)))

pA = plot(t, y; color = COR.tinta2, lw = 2.5, label = "medida", title = "A. Stergiopulos 1999 (humano)",
          xlabel = "s", ylabel = "mmHg", legend = :topright)
for (k, nome) in enumerate(("2 elementos", "3 elementos", "4 elementos"))
    plot!(pA, t, resultadosA[nome].y; color = (COR.s1, COR.s2, COR.s3)[k], lw = 1.5,
          label = @sprintf("%s (RMSE %.1f)", nome, resultadosA[nome].rmse))
end
n = min(length(yB), 800)
tB = h .* (0:n-1)
pB = plot(tB, yB[1:n]; color = COR.tinta2, lw = 2.5, label = "medida", title = "B. In vivo (porco)",
          xlabel = "s", ylabel = "mmHg", legend = :topright)
plot!(pB, tB, yBhat[1:n]; color = COR.s3, lw = 1.5, label = "Windkessel 4 elementos")
savefig(plot(pA, pB; layout = (1, 2), size = (1150, 400), bottom_margin = 6Plots.mm, left_margin = 6Plots.mm),
        saida("identificacao_windkessel.png"))
println("Figura: ", normpath(saida("identificacao_windkessel.png")))
