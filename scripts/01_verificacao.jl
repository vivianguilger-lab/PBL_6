# =============================================================================
# 01_verificacao.jl — VERIFICAÇÃO: "resolvi as equações corretamente?"
# (ASME V&V 40 — Resumo, Seção 8.5)
#
#   A. Windkessel de 2 elementos (windkessel-id, G2) na diástole: P(t) = P0·exp(−t/RpC)
#   B. NeoPKPD (Schnider, bolus) contra a solução exata A(t) = exp(K·t)·A0
#   C. Windkessel de 4 elementos: PAM = Rp·Q̄ e estado periódico (solvex0)
#
# Uso: julia --project=. scripts/01_verificacao.jl
# =============================================================================
using PropofolWindkessel
using LinearAlgebra, Statistics
include(joinpath(@__DIR__, "graficos.jl"))

# --- A. Windkessel de 2 elementos: descarga sem fluxo de entrada ---------------
Rp, C, h = 1.0, 2.0, 0.005                    # valores do Resumo: τ = Rp·C = 2 s
sys = G2([Rp, C], h)
x0 = [C * 100.0]                              # volume no capacitor ⇒ P0 = 100 mmHg
y, = lsim_wk(sys, zeros(1200), x0)
t = h .* (0:1199)
exata = 100 .* exp.(-t ./ (Rp * C))
errA = maximum(abs.(y .- exata) ./ exata)
@printf("A. Windkessel 2 elementos (τ = Rp·C = %.1f s): erro relativo máximo = %.2e\n", Rp * C, errA)

# --- B. NeoPKPD contra exp(K·t) -----------------------------------------------------
ref = Paciente(idade = 53, peso = 77, altura = 177, sexo = :M, pas = 125, pad = 75, pam = 92, fc = 70)
pk = parametros_pk(ref, :schnider)
D = 2.0 * ref.peso
r = simular_propofol(ref, [bolus(0.0, D; duracao = 0.0)]; tfim = 20.0)
K = matriz_pk(pk)
Cp_exa = [(exp(K * ti) * [D, 0, 0])[1] / pk.V1 for ti in r.t]
sel = r.t .>= 0.05
errB = maximum(abs.(r.Cp[sel] .- Cp_exa[sel]) ./ Cp_exa[sel])
@printf("B. NeoPKPD (Schnider, bolus %.0f mg): erro relativo máximo em Cp = %.2e\n", D, errB)

# --- C. Windkessel 4 elementos (Stergiopulos 1999) ------------------------------------
θ, _ = params_stergiopulos1999total()
tq, pq, q = dados_stergiopulos1999()
sys4 = G(θ, h)
x0 = solvex0(sys4, q)
y4 = simulate_wk(θ, x0, q, h)
errC = abs(mean(y4) - θ[1] * mean(q)) / (θ[1] * mean(q))
errP = norm(estado_final(sys4, q, x0) .- x0) / norm(x0)
@printf("C. Windkessel 4 elementos: |PAM − Rp·Q̄|/PAM = %.2e; periodicidade do estado = %.2e\n", errC, errP)

pA = plot(t, exata; color = COR.s1, label = "exata  P0·e^(−t/τ)", title = "A. Diástole (Windkessel 2 elementos)",
          xlabel = "s", ylabel = "mmHg")
scatter!(pA, t[1:60:end], y[1:60:end]; color = COR.s2, ms = 4, msw = 0, label = "windkessel-id (G2)")
pB = plot(r.t, Cp_exa; color = COR.s1, label = "exata  exp(K·t)·A0", yscale = :log10,
          title = "B. NeoPKPD: Cp após bolus", xlabel = "min", ylabel = "µg/ml")
scatter!(pB, r.t[sel][1:60:end], r.Cp[sel][1:60:end]; color = COR.s2, ms = 4, msw = 0, label = "NeoPKPD")
pC = plot(tq, pq; color = COR.s1, label = "medida (Stergiopulos)", title = "C. Windkessel 4 elementos",
          xlabel = "s", ylabel = "mmHg", legend = :topright)
plot!(pC, tq, y4; color = COR.s2, label = "modelo (θ publicado)")
savefig(plot(pA, pB, pC; layout = (1, 3), size = (1150, 380), bottom_margin = 6Plots.mm, left_margin = 6Plots.mm),
        saida("verificacao.png"))
println("Figura: ", normpath(saida("verificacao.png")))
