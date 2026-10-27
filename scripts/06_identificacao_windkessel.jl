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
#
# Não precisa de arquivo de entrada: os dados ficam em dados/windkessel-id/.
# Saída: θ̂ estimado e RMSE (raiz do erro quadrático médio, em mmHg) na tela e a
# figura resultados/identificacao_windkessel.png (medida × modelo).
# θ = [Rp, C, Rc, L]: resistência periférica, complacência, resistência
# característica e inertância. O modelo de 2 elementos usa só [Rp, C] e o de 3,
# [Rp, C, Rc]. Em A o fluxo está em ml/s; em B, em L/min (θ̂ em outra unidade).
# =============================================================================
# carrega o pacote, estatística (mean) e as funções de figura
using PropofolWindkessel, Statistics
include(joinpath(@__DIR__, "graficos.jl"))
# passo de tempo dos dados: 0,005 s (5 ms)
h = 0.005

# --- A. Stergiopulos 1999 (mL/s) ------------------------------------------------
# t = tempo (s), y = pressão aórtica medida (mmHg), u = fluxo aórtico (ml/s)
t, y, u = dados_stergiopulos1999()
# θ publicado por Stergiopulos (1999), para comparar no fim
θpub, _ = params_stergiopulos1999total()
# guarda os resultados de cada modelo (2, 3 e 4 elementos)
resultadosA = Dict{String,Any}()
# Para cada modelo: (nome, função que monta o modelo, chute inicial θ0).
# O chute de Rp é pressão média ÷ fluxo médio; os outros são valores simples.
for (nome, G_, θ0) in (("2 elementos", G2, [mean(y) / mean(u), 1.0]),
                       ("3 elementos", G3, [mean(y) / mean(u), 1.0, 0.05]),
                       ("4 elementos", G,  [mean(y) / mean(u), 1.0, 0.05, 0.005]))
    # estado inicial periódico para o chute inicial
    x00 = solvex0(G_(θ0, h), u)
    # identify ajusta θ para a pressão simulada ficar o mais perto possível da
    # medida (mínimos quadrados, método de Newton com derivadas automáticas).
    # x0solve = true recalcula o estado inicial periódico a cada tentativa.
    phat, _ = identify(θ0, x00, u, y, h; x0id = false, x0solve = true, G_)
    # deal separa o vetor ajustado em θ̂ (parâmetros) e estado inicial
    θhat, _ = deal(phat, PropofolWindkessel.IdConditions(θ0, x00, ones(Int, length(θ0)), false, true, G_))
    # pressão simulada com θ̂ e o erro RMSE em mmHg
    yhat = simulate_wk(θhat, solvex0(G_(θhat, h), u), u, h; G_)
    rmse = sqrt(mean((y .- yhat) .^ 2))
    # guarda θ̂, a curva simulada e o RMSE (NamedTuple: campos com nome)
    resultadosA[nome] = (θ = θhat, y = yhat, rmse = rmse)
    @printf("A. Stergiopulos — %s: θ̂ = %s  RMSE = %.2f mmHg\n", nome, string(round.(θhat; sigdigits = 3)), rmse)
end
println("   θ publicado (4 elementos) = ", round.(θpub; sigdigits = 3))

# --- B. dados in vivo (L/min) ------------------------------------------------------
# B. Dados in vivo de porcos: 4 pares de arquivos de pressão (P) e fluxo (Q).
# vcat(...) junta os 4 trechos num único vetor; `...` espalha a lista.
P, Q = importPQ(joinpath(WINDKESSEL_ID_DIR, "invivo_data", "Ao_baseline1_Ao"), 4)
uB = vcat([q.vb_va for q in Q]...)
yB = vcat([p.vb_va for p in P]...)
# chute inicial e identificação do Windkessel de 4 elementos (como em A)
θ0 = [mean(yB) / mean(uB), 1.0, 1.0, 1.0]
x00 = solvex0(G(θ0, h), uB)
phat, _ = identify(θ0, x00, uB, yB, h; x0id = false, x0solve = true)
θB, _ = deal(phat, PropofolWindkessel.IdConditions(θ0, x00, [1, 1, 1, 1], false, true))
yBhat = simulate_wk(θB, solvex0(G(θB, h), uB), uB, h)
@printf("B. in vivo — 4 elementos: θ̂ (L/min) = %s  RMSE = %.2f mmHg\n", string(round.(θB; sigdigits = 3)),
        sqrt(mean((yB .- yBhat) .^ 2)))

# Painel A: pressão medida (linha grossa) e curvas dos modelos de 2, 3 e 4
# elementos. Painel B: dados in vivo e o modelo de 4 elementos.
pA = plot(t, y; color = COR.tinta2, lw = 2.5, label = "medida", title = "A. Stergiopulos 1999 (humano)",
          xlabel = "s", ylabel = "mmHg", legend = :topright)
for (k, nome) in enumerate(("2 elementos", "3 elementos", "4 elementos"))
    plot!(pA, t, resultadosA[nome].y; color = (COR.s1, COR.s2, COR.s3)[k], lw = 1.5,
          label = @sprintf("%s (RMSE %.1f)", nome, resultadosA[nome].rmse))
end
# mostra só os primeiros 800 pontos (4 s) para o painel B ficar legível
n = min(length(yB), 800)
tB = h .* (0:n-1)
pB = plot(tB, yB[1:n]; color = COR.tinta2, lw = 2.5, label = "medida", title = "B. In vivo (porco)",
          xlabel = "s", ylabel = "mmHg", legend = :topright)
plot!(pB, tB, yBhat[1:n]; color = COR.s3, lw = 1.5, label = "Windkessel 4 elementos")
# junta os dois painéis e grava em resultados/
savefig(plot(pA, pB; layout = (1, 2), size = (1150, 400), bottom_margin = 6Plots.mm, left_margin = 6Plots.mm),
        saida("identificacao_windkessel.png"))
println("Figura: ", normpath(saida("identificacao_windkessel.png")))
