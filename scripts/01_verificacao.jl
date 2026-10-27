# =============================================================================
# 01_verificacao.jl — VERIFICAÇÃO: "resolvi as equações corretamente?"
# (ASME V&V 40 — Resumo, Seção 8.5)
#
#   A. Windkessel de 2 elementos (windkessel-id, G2) na diástole: P(t) = P0·exp(−t/RpC)
#   B. NeoPKPD (Schnider, bolus) contra a solução exata A(t) = exp(K·t)·A0
#   C. Windkessel de 4 elementos: PAM = Rp·Q̄ e estado periódico (solvex0)
#
# Uso: julia --project=. scripts/01_verificacao.jl
#
# Como reproduzir: na pasta do projeto, rode o comando acima. Não precisa de
# arquivo de entrada. Imprime na tela os três erros relativos e grava a figura
# resultados/verificacao.png (painéis A, B e C).
# Erros muito pequenos (da ordem de 1e-10 ou menores) mostram que o código
# resolve as equações corretamente. Isto é VERIFICAÇÃO (a matemática do
# código); não é VALIDAÇÃO (comparação com pacientes reais).
# =============================================================================
# Carrega o pacote do projeto (pasta src/) e duas bibliotecas padrão do Julia:
# LinearAlgebra (norm = tamanho de um vetor) e Statistics (mean = média).
using PropofolWindkessel
using LinearAlgebra, Statistics
# include(...) executa graficos.jl: cores, pasta de saída e funções de figura
include(joinpath(@__DIR__, "graficos.jl"))

# --- A. Windkessel de 2 elementos: descarga sem fluxo de entrada ---------------
# Teste A: sem fluxo entrando (como na diástole), a pressão do Windkessel de 2
# elementos cai como exponencial: P(t) = P0·exp(−t/(Rp·C)). Comparamos a
# simulação com essa fórmula exata. h = passo de tempo (0,005 s = 5 ms).
# `a, b, c = 1, 2, 3` atribui vários valores de uma vez.
Rp, C, h = 1.0, 2.0, 0.005                    # valores do Resumo: τ = Rp·C = 2 s
# G2 monta o modelo discreto do Windkessel de 2 elementos com θ = [Rp, C]
sys = G2([Rp, C], h)
x0 = [C * 100.0]                              # volume no capacitor ⇒ P0 = 100 mmHg
# Simula 1200 passos (6 s) com fluxo zero, partindo do estado x0.
# `y, = ...` guarda só o primeiro resultado (a pressão) e descarta o resto.
y, = lsim_wk(sys, zeros(1200), x0)
# Vetor de tempos. O ponto em `.*`, `./` e `exp.` aplica a conta a cada
# elemento do vetor.
t = h .* (0:1199)
exata = 100 .* exp.(-t ./ (Rp * C))
# maior erro relativo entre a simulação e a fórmula exata
errA = maximum(abs.(y .- exata) ./ exata)
# @printf imprime com formato: %.1f = 1 casa decimal; %.2e = notação científica
@printf("A. Windkessel 2 elementos (τ = Rp·C = %.1f s): erro relativo máximo = %.2e\n", Rp * C, errA)

# --- B. NeoPKPD contra exp(K·t) -----------------------------------------------------
# Teste B: o paciente de referência de Schnider (53 anos, 77 kg, 177 cm) recebe
# um bolus instantâneo de 2 mg/kg. O modelo de 3 compartimentos é linear, então
# a quantidade de propofol tem solução exata A(t) = exp(K·t)·A0. Comparamos a
# concentração plasmática (Cp) do NeoPKPD com essa solução.
ref = Paciente(idade = 53, peso = 77, altura = 177, sexo = :M, pas = 125, pad = 75, pam = 92, fc = 70)
# parâmetros PK de Schnider para este paciente (V1, V2, V3, CL, Q2, Q3, ke0)
pk = parametros_pk(ref, :schnider)
# dose total em mg (2 mg/kg × peso)
D = 2.0 * ref.peso
# simula 20 min de PK no NeoPKPD (bolus de duração zero no tempo 0)
r = simular_propofol(ref, [bolus(0.0, D; duracao = 0.0)]; tfim = 20.0)
# K = matriz de trocas entre os 3 compartimentos (por minuto)
K = matriz_pk(pk)
# Cp exata em cada tempo: exp(K·t) (exponencial de matriz) vezes a dose no
# compartimento central, dividida por V1. `[... for ti in r.t]` monta um vetor.
Cp_exa = [(exp(K * ti) * [D, 0, 0])[1] / pk.V1 for ti in r.t]
# ignora os primeiros 0,05 min, em que Cp muda muito rápido
sel = r.t .>= 0.05
# maior erro relativo em Cp
errB = maximum(abs.(r.Cp[sel] .- Cp_exa[sel]) ./ Cp_exa[sel])
@printf("B. NeoPKPD (Schnider, bolus %.0f mg): erro relativo máximo em Cp = %.2e\n", D, errB)

# --- C. Windkessel 4 elementos (Stergiopulos 1999) ------------------------------------
# Teste C: com os parâmetros publicados por Stergiopulos (1999) e o fluxo
# aórtico medido, a PAM do modelo deve ser exatamente Rp × fluxo médio (Q̄),
# propriedade do Windkessel (|Z(0)| = Rp). Também conferimos que o estado
# inicial calculado por solvex0 se repete no fim do batimento (regime periódico).
# θ = [Rp, C, Rc, L]; `_` recebe um valor que não será usado.
θ, _ = params_stergiopulos1999total()
# tq = tempo (s), pq = pressão medida (mmHg), q = fluxo medido (ml/s)
tq, pq, q = dados_stergiopulos1999()
# monta o Windkessel de 4 elementos e acha o estado inicial periódico
sys4 = G(θ, h)
x0 = solvex0(sys4, q)
# pressão simulada pelo modelo para o mesmo fluxo
y4 = simulate_wk(θ, x0, q, h)
# erro relativo entre a PAM simulada e Rp·Q̄
errC = abs(mean(y4) - θ[1] * mean(q)) / (θ[1] * mean(q))
# quanto o estado no fim do batimento difere do estado inicial
errP = norm(estado_final(sys4, q, x0) .- x0) / norm(x0)
@printf("C. Windkessel 4 elementos: |PAM − Rp·Q̄|/PAM = %.2e; periodicidade do estado = %.2e\n", errC, errP)

# --- figura: um painel por teste (A, B e C) ----------------------------------
# Linha = solução exata ou dado medido; pontos ou 2ª linha = resultado do código.
# Funções terminadas em `!` (scatter!, plot!) desenham sobre um gráfico já feito.
# t[1:60:end] pega 1 ponto a cada 60, para o gráfico não ficar poluído.
pA = plot(t, exata; color = COR.s1, label = "exata  P0·e^(−t/τ)", title = "A. Diástole (Windkessel 2 elementos)",
          xlabel = "s", ylabel = "mmHg")
scatter!(pA, t[1:60:end], y[1:60:end]; color = COR.s2, ms = 4, msw = 0, label = "windkessel-id (G2)")
pB = plot(r.t, Cp_exa; color = COR.s1, label = "exata  exp(K·t)·A0", yscale = :log10,
          title = "B. NeoPKPD: Cp após bolus", xlabel = "min", ylabel = "µg/ml")
scatter!(pB, r.t[sel][1:60:end], r.Cp[sel][1:60:end]; color = COR.s2, ms = 4, msw = 0, label = "NeoPKPD")
pC = plot(tq, pq; color = COR.s1, label = "medida (Stergiopulos)", title = "C. Windkessel 4 elementos",
          xlabel = "s", ylabel = "mmHg", legend = :topright)
plot!(pC, tq, y4; color = COR.s2, label = "modelo (θ publicado)")
# junta os três painéis lado a lado e grava em resultados/verificacao.png
savefig(plot(pA, pB, pC; layout = (1, 3), size = (1150, 380), bottom_margin = 6Plots.mm, left_margin = 6Plots.mm),
        saida("verificacao.png"))
# mostra onde a figura foi gravada
println("Figura: ", normpath(saida("verificacao.png")))
