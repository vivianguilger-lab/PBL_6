# =============================================================================
# extras/windkessel_csm.jl — o MESMO circuito montado com CirculatorySystemModels.jl
#
# CirculatorySystemModels.jl (TS-CUBED, Sheffield Hallam) é a biblioteca Julia de
# componentes 0-D da circulação construída sobre o ModelingToolkit.jl. Em vez de
# escrever as EDOs à mão, ligamos componentes como num simulador de circuitos;
# o ModelingToolkit aplica as leis de Kirchhoff (`connect`) e reduz o sistema.
#
# Aqui o circuito do paciente virtual de referência (sem fármaco, sem
# barorreflexo) é montado com a biblioteca e comparado com a implementação à
# mão do pacote PropofolWindkessel: VERIFICAÇÃO CRUZADA entre duas
# implementações independentes das mesmas equações.
#
# Uso (ambiente próprio, porque o ModelingToolkit é grande):
#   julia --project=extras -e 'using Pkg; Pkg.instantiate()'
#   julia --project=extras extras/windkessel_csm.jl
# =============================================================================
using CirculatorySystemModels
using ModelingToolkit
using OrdinaryDiffEqTsit5
using PropofolWindkessel
using Printf
using Statistics

# --- 1. paciente virtual de referência (pacote à mão) -----------------------------
ref = Paciente(nome = "referência", idade = 53, peso = 77, altura = 177, sexo = :M, pas = 125, pad = 75, fc = 70)
pv = calibrar(ref)
cv = pv.cv
T = 60 / cv.FC
Ts = 0.35 * sqrt(T)                     # mesma ativação cossenoidal (Shi et al.)
Tr = 0.5 * Ts

# --- 2. componentes da biblioteca ---------------------------------------------------
@independent_variables t
@named LV = ShiChamber(V₀ = cv.Vd, p₀ = 0.0, Eₘᵢₙ = cv.Emin, Eₘₐₓ = cv.Emax, τ = T,
                       τₑₛ = Ts, τₑₚ = Ts + Tr, Eshift = 0.0)
# valva aórtica (diodo) + impedância característica Zc num único elemento, como
# no exemplo de Bjørdalsbakke da biblioteca: dois resistores em série com um diodo
# criariam um nó algébrico (sistema DAE).
@named AV = ResistorDiode(R = cv.Rav + cv.Zc)
@named Ca = Compliance(C = cv.Ca, V₀ = cv.Va0)    # complacência arterial (AC)
@named Rs = Resistor(R = cv.Rs)                  # resistência periférica total (TPR)
@named Cv = Compliance(C = cv.Cv, V₀ = cv.Vv0)    # complacência venosa (VC) + V0 venoso
@named MV = ResistorDiode(R = cv.Rin)            # enchimento (valva AV + veias centrais)

# --- 3. ligações (Kirchhoff) -----------------------------------------------------
eqs = [connect(LV.out, AV.in)
       connect(AV.out, Ca.in)
       connect(Ca.out, Rs.in)
       connect(Rs.out, Cv.in)
       connect(Cv.out, MV.in)
       connect(MV.out, LV.in)]
@named _circ = System(eqs, t)
@named circ = compose(_circ, [LV, AV, Ca, Rs, Cv, MV])
sys = mtkcompile(circ)
println("Estados após a redução simbólica: ", unknowns(sys))

# condição inicial: o mesmo estado de regime do paciente virtual
u0 = [LV.V => pv.u0[I_VLV], Ca.V => pv.u0[I_VA], Cv.V => pv.u0[I_VV]]
nbat = 30
prob = ODEProblem(sys, u0, (0.0, nbat * T))
sol = solve(prob, Tsit5(); reltol = 1e-8, abstol = 1e-8, dtmax = 0.005, saveat = T / 400)

# --- 4. mesma simulação com o pacote à mão ------------------------------------------
m = ModeloInducao(pv.pk, pv.pd, cv, pv.baro, 0.0, TrechoInfusao[], false, false)
sol2 = integrar(m, pv.u0, (0.0, nbat * T); saveat = T / 400, reltol = 1e-8, abstol = 1e-8)

# último batimento de cada uma
ult = findall(x -> x >= (nbat - 1) * T, sol.t)
Q_csm = sol[AV.q][ult]
Pao_csm = sol[Ca.p][ult] .+ cv.Zc .* Q_csm      # pressão na raiz da aorta = Pc + Zc·Q
ult2 = findall(x -> x >= (nbat - 1) * T, sol2.t)
alg = [algebrico(u, m, t) for (t, u) in zip(sol2.t[ult2], sol2.u[ult2])]
Pao_mao = [a.Pao for a in alg]

res(nome, a, b) = @printf("  %-22s %8.2f   %8.2f   (dif. %.2e)\n", nome, a, b, abs(a - b))
println("\nÚltimo batimento            CSM.jl    à mão")
res("PAS aórtica (mmHg)", maximum(Pao_csm), maximum(Pao_mao))
res("PAD aórtica (mmHg)", minimum(Pao_csm), minimum(Pao_mao))
res("PAM (mmHg)", mean(Pao_csm), mean(Pao_mao))
res("VS (ml)", sum(Q_csm[1:end-1]) * T / 400, sum(a.Qao for a in alg[1:end-1]) * T / 400)
