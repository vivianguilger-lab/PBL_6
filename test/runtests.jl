# =============================================================================
# Testes de VERIFICAÇÃO ("resolvi as equações corretamente?") e de
# plausibilidade fisiológica (validade de face).  Rodar com:
#     julia --project=. -e 'using Pkg; Pkg.test()'
# =============================================================================
using PropofolWindkessel
using LinearAlgebra
using Statistics
using Test

const REF = Paciente(nome = "referência Schnider", idade = 53, peso = 77, altura = 177, sexo = :M,
                     pas = 125, pad = 75, fc = 70)

@testset "PropofolWindkessel" begin

    @testset "Hill" begin
        @test hill(0.0, 3.0, 2.0) == 0.0
        @test hill(3.0, 3.0, 2.0) ≈ 0.5
        @test hill(1e6, 3.0, 2.0) ≈ 1.0 atol = 1e-6
        @test issorted([hill(c, 2.0, 1.5) for c in 0:0.1:10])
    end

    @testset "Schnider — valores de referência" begin
        pk = schnider(REF)
        lbm = massa_magra(REF)
        @test lbm ≈ 1.1 * 77 - 128 * (77 / 177)^2
        @test pk.V1 ≈ 4.27
        @test pk.V2 ≈ 18.9
        @test pk.V3 ≈ 238
        @test pk.Cl1 ≈ 1.89 - 0.0681 * (lbm - 59)
        @test pk.Cl2 ≈ 1.29
        @test pk.ke0_bis * 60 ≈ 0.456
        # Kazama 1999: t½ke0 da PA mais lento que o do BIS e crescente com a idade
        @test meia_vida_ke0_hemodinamica(30) ≈ 5.7
        @test meia_vida_ke0_hemodinamica(77) ≈ 10.2
        @test log(2) / pk.ke0_hemo > log(2) / pk.ke0_bis
    end

    pv = calibrar(REF)

    @testset "PK: solução numérica = solução exata (exponencial de matriz)" begin
        D = 140.0                         # mg em bolus instantâneo
        u0 = copy(pv.u0); u0[I_A1] = D
        m = ModeloInducao(pv.pk, pv.pd, pv.cv, pv.baro, 0.0, TrechoInfusao[], false, false)
        sol = integrar(m, u0, (0.0, 600.0); saveat = 60.0)
        K = matriz_pk(pv.pk)
        for (t, u) in zip(sol.t, sol.u)
            exata = exp(K * t) * [D, 0.0, 0.0]
            @test u[I_A1:I_A3] ≈ exata rtol = 1e-5
        end
    end

    @testset "Windkessel de 2 elementos: decaimento diastólico P0·exp(−t/RC)" begin
        # coração parado (Emax = Emin, VE vazio), veias com complacência enorme (Pv ≈ 0):
        # a árvore arterial é um capacitor descarregando por um resistor.
        cv = ParametrosCV(Emax = 0.06, Emin = 0.06, Cv = 1e9, Vv0 = 0.0, Ca = 1.5, Rs = 1.2)
        u0 = zeros(NESTADOS)
        P0 = 100.0
        u0[I_VLV] = cv.Vd
        u0[I_VA] = cv.Va0 + cv.Ca * P0
        u0[I_VV] = 0.0
        m = ModeloInducao(pv.pk, pv.pd, cv, pv.baro, 0.0, TrechoInfusao[], false, false)
        sol = integrar(m, u0, (0.0, 5.0); saveat = 0.5, reltol = 1e-10, abstol = 1e-10)
        tau = cv.Rs * cv.Ca
        for (t, u) in zip(sol.t, sol.u)
            Pc = (u[I_VA] - cv.Va0) / cv.Ca
            @test Pc ≈ P0 * exp(-t / tau) rtol = 1e-6
        end
    end

    @testset "Parada cardíaca: todas as pressões → Pms (Resumo, eq. 8)" begin
        cv = com(pv.cv; Emax = pv.cv.Emin)
        m = ModeloInducao(pv.pk, pv.pd, cv, pv.baro, 0.0, TrechoInfusao[], false, false)
        sol = integrar(m, pv.u0, (0.0, 120.0); saveat = 120.0, reltol = 1e-10, abstol = 1e-10)
        a = algebrico(sol.u[end], m, 120.0)
        Pms = pressao_media_enchimento(cv)
        @test a.Pc ≈ Pms atol = 1e-2
        @test a.Pv ≈ Pms atol = 1e-2
        @test a.Plv ≈ Pms atol = 1e-2
        @test 4 < Pms < 15                # faixa fisiológica (Resumo: 7–10 mmHg)
    end

    @testset "Calibração do paciente virtual" begin
        @test pv.convergiu
        @test pv.obtido.PAM ≈ pv.alvo.PAM rtol = 0.01
        @test pv.obtido.PASao - pv.obtido.PADao ≈ pv.alvo.PP rtol = 0.01   # PP aórtica
        @test pv.obtido.PAS ≈ REF.pas rtol = 0.01                          # equivalentes de manguito
        @test pv.obtido.PAD ≈ REF.pad rtol = 0.01
        @test pv.obtido.DC ≈ pv.alvo.DC rtol = 0.01
        @test 0.5 < tau_diastolico(pv) < 4.0
    end

    r0 = simular_inducao(pv; dose_mgkg = 0.0, duracao = 300.0)
    r2 = simular_inducao(pv; dose_mgkg = 2.0, duracao = 20 * 60.0)

    @testset "Conservação de volume e estabilidade da linha de base" begin
        @test volume_total(r2.u_final) ≈ volume_total(pv.u0) atol = 1e-4
        @test maximum(abs.(r0.bat.PAM .- pv.obtido.PAM)) < 1.5
        @test maximum(abs.(r0.bat.FC .- pv.cv.FC)) < 1.5
    end

    @testset "Plausibilidade: 2 mg/kg em adulto de referência" begin
        mt = metricas(r2)
        @test 15 < mt.queda_pct < 45            # queda típica de 20–40% da PAM
        @test 1.5 < mt.t_nadir < 15             # nadir nos primeiros minutos
        @test mt.bis_min < 60                   # hipnose atingida
        @test mt.t_hipnose < 2.5                # perda de consciência < ~1–2 min
        @test mt.ceb_max > mt.ceh_max           # efeito vascular mais lento (Kazama)
    end

    @testset "Monotonicidade dose–resposta e efeito da velocidade de injeção" begin
        nadires = [metricas(simular_inducao(pv; dose_mgkg = d, saveat = 0.01)).pam_nadir for d in (1.0, 2.0, 3.0)]
        @test issorted(nadires; rev = true)
        lento = metricas(simular_inducao(pv; dose_mgkg = 2.0, duracao_bolus = 120.0, saveat = 0.01)).pam_nadir
        @test lento > nadires[2]
    end

    @testset "Idade: idoso cai mais com a mesma dose por kg" begin
        jovem = Paciente(idade = 25, peso = 70, altura = 175, sexo = :M, pas = 125, pad = 75, fc = 70)
        idoso = Paciente(idade = 78, peso = 70, altura = 175, sexo = :M, pas = 125, pad = 75, fc = 70)
        qj = metricas(simular_inducao(calibrar(jovem); dose_mgkg = 2.0, saveat = 0.01)).queda_pct
        qi = metricas(simular_inducao(calibrar(idoso); dose_mgkg = 2.0, saveat = 0.01)).queda_pct
        @test qi > qj
    end

    @testset "Previsão com incerteza farmacológica" begin
        inc = previsao_com_incerteza(pv; n = 16, dose_mgkg = 2.0)
        @test all(inc.p05 .<= inc.mediana .<= inc.p95)
        @test 0 <= inc.prob_pam55 <= inc.prob_pam65 <= 1
        @test inc.nadir_p05 <= inc.nadir_mediano <= inc.nadir_p95
    end

    @testset "Ponto de funcionamento (Guyton) satisfaz as duas retas" begin
        g = estado_guyton(pv.cv; Pa = pv.obtido.PAM)
        @test retorno_venoso(g.RAP, g.Pms, g.RVR) ≈ g.Q rtol = 1e-10
        @test debito_cardiaco_curso(g.RAP, g.FC, g.CD, g.CS, g.Pa) ≈ g.Q rtol = 1e-10
    end
end
