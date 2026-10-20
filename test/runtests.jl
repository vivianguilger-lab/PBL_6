# =============================================================================
# Testes de VERIFICAÇÃO ("resolvi as equações corretamente?") e de
# plausibilidade. Inclui os testes do PropofolPBL do grupo (Windkessel do
# windkessel-id, PK/PD do NeoPKPD e acoplamento).
#     julia --project=. -e 'using Pkg; Pkg.test()'
# =============================================================================
using PropofolWindkessel
using LinearAlgebra
using Statistics
using Test

const h = 0.005
const REF = Paciente(nome = "referência Schnider", idade = 53, peso = 77, altura = 177, sexo = :M,
                     pas = 125, pad = 75, fc = 70)

@testset "PropofolWindkessel" begin

    # ---------------------------------------------------------------- grupo
    @testset "Windkessel (windkessel-id)" begin
        θ, _ = params_stergiopulos1999total()
        _, _, q = dados_stergiopulos1999()
        u = repeat(q, 5)

        @testset "solvex0 dá estado periódico" begin
            sys = G(θ, h)
            x0 = solvex0(sys, q)
            @test estado_final(sys, q, x0) ≈ x0 rtol = 1e-8
        end
        @testset "PAM do modelo ≈ Rp · fluxo médio" begin
            y = simulate_wk(θ, solvex0(G(θ, h), q), q, h)
            @test mean(y) ≈ θ[1] * mean(q) rtol = 0.02
        end
        @testset "bode: |Z(0)| = Rp e |Z(∞)| → Rc" begin
            mag, _ = bode_wk(θ, [1e-6, 1e6])
            @test mag[1] ≈ θ[1] rtol = 1e-4
            @test mag[2] ≈ θ[3] rtol = 1e-3
            mag3, _ = bode_wk(θ[1:3], [1e6]; G_ = G3)
            @test mag3[1] ≈ θ[3] rtol = 1e-3
        end
        @testset "identify recupera parâmetros de dados sintéticos" begin
            y = simulate_wk(θ, solvex0(G(θ, h), u), u, h)
            θ0 = θ .* [1.0, 1.4, 0.7, 1.5]
            x00 = solvex0(G(θ0, h), u)
            phat, res = identify(θ0, x00, u, y, h; x0id = false, x0solve = true)
            θhat, _ = deal(phat, IdConditions(θ0, x00, [1, 1, 1, 1], false, true))
            @test θhat ≈ θ rtol = 1e-3
        end
        @testset "importPQ lê os dados in vivo" begin
            P, Q = importPQ(joinpath(WINDKESSEL_ID_DIR, "invivo_data", "Ao_baseline1_Ao"), 4)
            @test length(P) == 4 && length(Q) == 4
            @test 50 < mean(P[1].vb_va) < 150
        end
    end

    @testset "Propofol (NeoPKPD)" begin
        pac = Paciente(idade = 53, peso = 77, altura = 177, sexo = :M, pas = 125, pad = 75, fc = 70)
        @testset "Schnider no paciente de referência" begin
            pk = parametros_pk(pac, :schnider)
            @test pk.V1 == 4.27
            @test pk.V2 ≈ 18.9
            @test pk.Q2 ≈ 1.29
            @test pk.CL ≈ 1.89 - 0.0681 * (massa_magra(pac) - 59) rtol = 1e-12
        end
        @testset "bolus de 2 mg/kg" begin
            r = simular_propofol(pac, [bolus(0.0, 2 * pac.peso)]; tfim = 20.0)
            ipico = argmax(r.Ce)
            @test 1.0 < r.t[ipico] < 4.0          # pico no sítio efetor em ~1,5–2 min (Schnider)
            @test argmax(r.Cp) < ipico            # histerese: Cp atinge o pico antes de Ce
            @test minimum(r.bis) < 60
            @test r.bis[1] ≈ 93.0
        end
        @testset "NeoPKPD = solução exata exp(K·t)·A0 (verificação)" begin
            pk = parametros_pk(pac, :schnider)
            D = 154.0
            r = simular_propofol(pac, [bolus(0.0, D; duracao = 0.0)]; tfim = 10.0)
            K = matriz_pk(pk)
            for t in (1.0, 2.0, 5.0, 10.0)
                i = argmin(abs.(r.t .- t))
                @test r.Cp[i] ≈ (exp(K * r.t[i]) * [D, 0, 0])[1] / pk.V1 rtol = 1e-3
            end
        end
        @testset "BIS de Eleveld" begin
            @test bis_eleveld(3.08, 35) ≈ 93 / 2
            @test bis_eleveld(0.0, 35) == 93.0
            @test bis_eleveld(2.0, 80) < bis_eleveld(2.0, 20)   # idoso mais sensível
        end
    end

    @testset "Acoplamento PK/PD → Windkessel (grupo)" begin
        pac = Paciente(idade = 45, peso = 70, altura = 170, sexo = :F, pas = 120, pad = 75, pam = 90.0, fc = 75.0)
        pkpd, bat, onda = simular_paciente(pac; tfim = 10.0, forma_de_onda = true)
        @test bat.PAM[1] ≈ pam_basal(pac) rtol = 0.03     # basal reproduz a PAM do paciente
        @test minimum(bat.PAM) < 0.85 * pam_basal(pac)     # hipotensão após a indução
        @test all(bat.PAS .> bat.PAM .> bat.PAD)
        @test length(onda.pressao) == length(onda.t_s)
        pac2 = Paciente(idade = 75, peso = 70, altura = 170, sexo = :F, pas = 120, pad = 75, pam = 90.0, fc = 75.0,
                        asa = 3, pvi = 20.0)
        @test sensibilidade_hemodinamica(pac2) > sensibilidade_hemodinamica(pac)
        _, bat2, _ = simular_paciente(pac2; tfim = 10.0)
        @test minimum(bat2.PAM) < minimum(bat.PAM)
    end

    # ---------------------------------------------------------------- extensões
    pv = calibrar(REF)

    @testset "Paciente virtual: PAM e pressão de pulso do monitor" begin
        @test pv.convergiu
        @test pv.obtido.PAM ≈ pam_basal(REF) rtol = 0.01
        @test pv.obtido.PP ≈ REF.pas - REF.pad rtol = 0.01
        @test 0.3 < tau_diastolico(pv) < 5.0
    end

    r0 = simular_inducao(pv; dose_mgkg = 0.0, manutencao_ugkgmin = 0.0, duracao = 5.0)
    r2 = simular_inducao(pv; dose_mgkg = 2.0)

    @testset "Linha de base estável sem fármaco" begin
        @test maximum(abs.(r0.bat.PAM .- pv.obtido.PAM)) < 0.5
        @test all(r0.pkpd.Cp .== 0)
    end

    @testset "Plausibilidade: 2 mg/kg em adulto de referência" begin
        mt = metricas(r2)
        @test 10 < mt.queda_pct < 45
        @test 0.5 < mt.t_nadir < 15
        @test mt.bis_min < 60
        @test mt.t_hipnose < 2.5
    end

    @testset "Monotonicidade dose–resposta" begin
        nadires = [metricas(simular_inducao(pv; dose_mgkg = d, manutencao_ugkgmin = 0.0)).pam_nadir for d in (1.0, 2.0, 3.0)]
        @test issorted(nadires; rev = true)
    end

    @testset "Previsão com incerteza farmacológica" begin
        inc = previsao_com_incerteza(pv; n = 16, dose_mgkg = 2.0)
        @test all(inc.p05 .<= inc.mediana .<= inc.p95)
        @test 0 <= inc.prob_pam55 <= inc.prob_pam65 <= 1
    end

    @testset "Plataforma: API (sem abrir servidor)" begin
        corpo = Dict("paciente" => Dict("idade" => "53", "sexo" => "M", "peso" => "77", "altura" => "177",
                                        "pas" => "125", "pad" => "75", "fc" => "70", "asa" => "2", "pvi" => "12"),
                     "protocolo" => Dict("dose_mgkg" => "2,0", "duracao_bolus_s" => "30", "modelo_pk" => "schnider"))
        r = PropofolWindkessel._api_simular(corpo)
        @test r["metricas"]["pam_nadir"] < r["metricas"]["pam_basal"]
        @test length(r["batimentos"]["t"]) == length(r["batimentos"]["pam"])
        d = PropofolWindkessel._api_doses(corpo)
        @test length(d["doses"]) == 11
        @test issorted([x["pam_nadir"] for x in d["doses"]]; rev = true)
        ruim = Dict("paciente" => merge(corpo["paciente"], Dict("pas" => "70")))
        @test_throws ArgumentError PropofolWindkessel._ler_paciente(ruim)
    end
end
