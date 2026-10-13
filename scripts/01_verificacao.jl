# =============================================================================
# 01_verificacao.jl — VERIFICAÇÃO: "resolvi as equações corretamente?"
# (ASME V&V 40 — Resumo, Seção 8.5)
#
# Três problemas com solução exata conhecida:
#   A. Windkessel de 2 elementos na diástole: P(t) = P0·exp(−t/RC)
#   B. Farmacocinética de Schnider (linear):  A(t) = exp(K·t)·A0
#   C. Parada cardíaca: todas as pressões convergem para
#      Pms = (Vtotal − V0)/(AC + VC + C_D)     (Resumo, eq. 8)
#
# Uso: julia --project=. scripts/01_verificacao.jl
# =============================================================================
using PropofolWindkessel
using LinearAlgebra
include(joinpath(@__DIR__, "graficos.jl"))

ref = Paciente(nome = "referência", idade = 53, peso = 77, altura = 177, sexo = :M, pas = 125, pad = 75, fc = 70)
pv = calibrar(ref)

# --- A. Windkessel de 2 elementos ---------------------------------------------
cv = ParametrosCV(Emax = 0.06, Emin = 0.06, Cv = 1e9, Vv0 = 0.0, Ca = 2.0, Rs = 1.0)   # valores do Resumo
u0 = zeros(NESTADOS); P0 = 100.0
u0[I_VLV] = cv.Vd; u0[I_VA] = cv.Va0 + cv.Ca * P0
m = ModeloInducao(pv.pk, pv.pd, cv, pv.baro, 0.0, TrechoInfusao[], false, false)
solA = integrar(m, u0, (0.0, 6.0); saveat = 0.25, reltol = 1e-10, abstol = 1e-10)
tau = cv.Rs * cv.Ca
PA_num = [(u[I_VA] - cv.Va0) / cv.Ca for u in solA.u]
PA_exa = P0 .* exp.(-solA.t ./ tau)
errA = maximum(abs.(PA_num .- PA_exa) ./ PA_exa)
@printf("A. Windkessel (τ = TPR·AC = %.1f s): erro relativo máximo = %.2e\n", tau, errA)

# --- B. PK de Schnider ----------------------------------------------------------
D = 2.0 * ref.peso
u0 = copy(pv.u0); u0[I_A1] = D
m = ModeloInducao(pv.pk, pv.pd, pv.cv, pv.baro, 0.0, TrechoInfusao[], false, false)
solB = integrar(m, u0, (0.0, 1200.0); saveat = 10.0, reltol = 1e-8, abstol = 1e-10)
K = matriz_pk(pv.pk)
Cp_num = [u[I_A1] / pv.pk.V1 for u in solB.u]
Cp_exa = [(exp(K * t) * [D, 0, 0])[1] / pv.pk.V1 for t in solB.t]
errB = maximum(abs.(Cp_num .- Cp_exa) ./ Cp_exa)
@printf("B. PK Schnider (bolus %.0f mg): erro relativo máximo em Cp = %.2e\n", D, errB)

# --- C. Parada cardíaca → Pms ---------------------------------------------------
cvC = com(pv.cv; Emax = pv.cv.Emin)
m = ModeloInducao(pv.pk, pv.pd, cvC, pv.baro, 0.0, TrechoInfusao[], false, false)
solC = integrar(m, pv.u0, (0.0, 30.0); saveat = 0.05, reltol = 1e-10, abstol = 1e-10)
alg = [algebrico(u, m, t) for (t, u) in zip(solC.t, solC.u)]
Pms = pressao_media_enchimento(cvC)
errC = maximum(abs.([alg[end].Pc, alg[end].Pv, alg[end].Plv] .- Pms))
@printf("C. Parada cardíaca: Pms teórico = %.3f mmHg; erro máximo após 30 s = %.2e mmHg\n", Pms, errC)
@printf("   volume total: variação máxima = %.2e ml (conservação de massa)\n",
        maximum(abs.(volume_total.(solC.u) .- volume_total(pv.u0))))

# --- figura --------------------------------------------------------------------
pA = plot(solA.t, PA_exa; color = COR.s1, label = "exata  P0·e^(−t/τ)", title = "A. Diástole do Windkessel",
          xlabel = "s", ylabel = "mmHg")
scatter!(pA, solA.t, PA_num; color = COR.s2, ms = 4, msw = 0, label = "numérica (Tsit5)")
pB = plot(solB.t ./ 60, Cp_exa; color = COR.s1, label = "exata  exp(K·t)·A0", yscale = :log10,
          title = "B. Schnider: Cp após bolus", xlabel = "min", ylabel = "µg/ml")
scatter!(pB, solB.t[1:6:end] ./ 60, Cp_num[1:6:end]; color = COR.s2, ms = 4, msw = 0, label = "numérica")
pC = plot(solC.t, [a.Pc for a in alg]; color = COR.s1, label = "arterial", title = "C. Parada cardíaca → Pms",
          xlabel = "s", ylabel = "mmHg", legend = :topright)
plot!(pC, solC.t, [a.Pv for a in alg]; color = COR.s2, label = "venosa")
plot!(pC, solC.t, [a.Plv for a in alg]; color = COR.s3, label = "ventricular")
hline!(pC, [Pms]; color = COR.muda, ls = :dash, lw = 1.2, label = "")
annotate!(pC, 30, Pms + 6, text(@sprintf("Pms = %.1f mmHg (eq. 8)", Pms), 8, COR.tinta2, :right))
savefig(plot(pA, pB, pC; layout = (1, 3), size = (1150, 380), bottom_margin = 6Plots.mm, left_margin = 6Plots.mm),
        saida("verificacao.png"))
println("Figura: ", normpath(saida("verificacao.png")))
