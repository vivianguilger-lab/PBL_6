# =============================================================================
# simulacao.jl — integração numérica e pós-processamento batimento a batimento
# =============================================================================

"""
Resultado de uma simulação: sinais amostrados (`t`, `sinais`) e valores por
batimento (`bat`), como um monitor multiparamétrico mostraria.
"""
struct Resultado
    t::Vector{Float64}
    sinais::NamedTuple
    bat::NamedTuple
    t_inducao::Float64
    modelo::ModeloInducao
    u_final::Vector{Float64}
end

"""
    integrar(m, u0, tspan; saveat=0.005)

Integra o modelo com Tsit5 (Runge–Kutta 5(4) adaptativo). As válvulas
(diodos) tornam o lado direito não suave; `dtmax` impede que o passo
"pule" a sístole e `tstops` força paradas nos instantes de mudança da infusão.
"""
function integrar(m::ModeloInducao, u0, tspan; saveat = 0.005, reltol = 1e-6, abstol = 1e-8)
    prob = ODEProblem(rhs!, copy(u0), tspan, m)
    tst = filter(t -> tspan[1] < t < tspan[2], instantes_infusao(m.infusao))
    return solve(prob, Tsit5(); saveat = saveat, tstops = tst, dtmax = 0.02,
                 reltol = reltol, abstol = abstol, maxiters = 10^8)
end

"""
    batimentos(t, U, sinais)

Detecta cada batimento pela passagem da fase por um inteiro e calcula
PAS/PAD aórticas (máx./mín. de Pao), PAM (média temporal exata via ∫Pao dt),
VS (∫Qao dt), FC e DC de cada batimento completo. `PAS`/`PAD` são os
equivalentes de manguito (ver `ObservacaoPA`); `PASao`/`PADao`, os aórticos.
"""
function batimentos(t, fase, intP, intQ, Pao; obs::ObservacaoPA = ObservacaoPA())
    n = length(t)
    kini = ceil(Int, fase[1] + 1e-9)
    kfim = floor(Int, fase[end])
    # instantes (interpolados) em que a fase cruza cada inteiro
    tb = Float64[]; Ib = Float64[]; Qb = Float64[]; ib = Int[]
    i = 1
    for k in kini:kfim
        while i < n && fase[i+1] < k
            i += 1
        end
        i >= n && break
        if fase[i] >= k     # cruzamento exatamente na primeira amostra
            push!(tb, t[i]); push!(Ib, intP[i]); push!(Qb, intQ[i]); push!(ib, i)
            continue
        end
        w = (k - fase[i]) / (fase[i+1] - fase[i])
        push!(tb, t[i] + w * (t[i+1] - t[i]))
        push!(Ib, intP[i] + w * (intP[i+1] - intP[i]))
        push!(Qb, intQ[i] + w * (intQ[i+1] - intQ[i]))
        push!(ib, i)
    end
    nb = length(tb) - 1
    tm = zeros(nb); PAS = zeros(nb); PAD = zeros(nb); PAM = zeros(nb)
    VS = zeros(nb); FC = zeros(nb); DC = zeros(nb)
    for j in 1:nb
        dt = tb[j+1] - tb[j]
        tm[j] = 0.5 * (tb[j] + tb[j+1])
        faixa = view(Pao, ib[j]:min(ib[j+1] + 1, n))
        PAS[j] = maximum(faixa)
        PAD[j] = minimum(faixa)
        PAM[j] = (Ib[j+1] - Ib[j]) / dt
        VS[j] = Qb[j+1] - Qb[j]
        FC[j] = 60 / dt
        DC[j] = VS[j] * FC[j] / 1000
    end
    # equivalentes de manguito (braquial) — equação de observação
    PPb = obs.amp .* (PAS .- PAD)
    PADb = PAM .- obs.ff .* PPb
    PASb = PADb .+ PPb
    return (t = tm, PAS = PASb, PAD = PADb, PAM = PAM, PASao = PAS, PADao = PAD, VS = VS, FC = FC, DC = DC)
end

"""Converte a solução do integrador em um `Resultado`."""
function Resultado(sol, m::ModeloInducao; t_inducao = 0.0)
    t = sol.t
    n = length(t)
    U = reduce(hcat, sol.u)
    Pao = zeros(n); Plv = zeros(n); Pv = zeros(n); Cp = zeros(n); BIS = zeros(n)
    FC = zeros(n); Rs = zeros(n); Emax = zeros(n); Vv0 = zeros(n); DCinst = zeros(n)
    for j in 1:n
        a = algebrico(view(U, :, j), m, t[j])
        Pao[j] = a.Pao; Plv[j] = a.Plv; Pv[j] = a.Pv; Cp[j] = a.Cp
        BIS[j] = a.ef.bis; FC[j] = a.FC; Rs[j] = a.Rs; Emax[j] = a.Emax; Vv0[j] = a.Vv0
    end
    sinais = (Pao = Pao, Plv = Plv, Pv = Pv, Cp = Cp, Ceb = U[I_CEB, :], Ceh = U[I_CEH, :],
              BIS = BIS, FC = FC, Rs = Rs, Emax = Emax, Vv0 = Vv0, Vlv = U[I_VLV, :])
    bat = batimentos(t, U[I_FASE, :], U[I_INTP, :], U[I_INTQ, :], Pao; obs = m.obs)
    return Resultado(collect(t), sinais, bat, t_inducao, m, sol.u[end])
end

"""
    simular_inducao(pv::PacienteVirtual; dose_mgkg=2.0, duracao_bolus=20.0,
                    manutencao_ugkgmin=0.0, t_inducao=60.0, duracao=20*60,
                    barorreflexo=true, saveat=0.005) -> Resultado

Simula 1 min de linha de base + `duracao` segundos após o início do bolus.
É a "previsão da curva de PAM para a dose escolhida, antes da injeção".
"""
function simular_inducao(pv; dose_mgkg = 2.0, duracao_bolus = 20.0, manutencao_ugkgmin = 0.0,
                         t_inducao = 60.0, duracao = 20 * 60.0, barorreflexo = true,
                         saveat = 0.005, pd = pv.pd, pk = pv.pk)
    tf = t_inducao + duracao
    tr = protocolo_inducao(pv.paciente; dose_mgkg, duracao_bolus, t_inicio = t_inducao,
                           manutencao_ugkgmin, t_final = tf)
    m = ModeloInducao(pk, pd, pv.cv, pv.baro, inot_Emax(pd, pv.paciente.idade), tr,
                      barorreflexo, true, pv.obs)
    sol = integrar(m, pv.u0, (0.0, tf); saveat)
    return Resultado(sol, m; t_inducao)
end

"""Média móvel simples (janela em número de amostras) — só para gráficos."""
function media_movel(x, janela)
    n = length(x)
    y = similar(x, Float64)
    s = 0.0
    for i in 1:n
        s += x[i]
        i > janela && (s -= x[i-janela])
        y[i] = s / min(i, janela)
    end
    return y
end
