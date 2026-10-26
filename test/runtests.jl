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
                     pas = 125, pad = 75, pam = 92, fc = 70)

# Parâmetros de Eleveld gerados executando as linhas 256–335 de pk_models.py do Python
# Anesthesia Simulator (commit 4afc03a): (idade, peso, altura, sexo, arterial, opioides,
# (CL, V1, Q2, V2, Q3, V3, ke0)). Script: docs/validacao.md, seção "Como reproduzir".
const ELEVELD_PAS = [
        (62, 72, 160, :F, true, false, (2.1450677945364385, 6.340142612686683, 1.3372361467303682, 17.19828856728441, 0.9437543998441148, 215.67134527591287, 0.14497537529787108)),
        (62, 72, 160, :F, true, true, (1.7968436137579569, 6.340142612686683, 1.3372361467303682, 17.19828856728441, 0.49640596145648325, 91.5720871040533, 0.14497537529787108)),
        (62, 72, 160, :F, false, false, (2.1450677945364385, 9.19980017095141, 0.9099896136155886, 17.19828856728441, 0.9437543998441148, 215.67134527591287, 1.2312977080093161)),
        (62, 72, 160, :F, false, true, (1.7968436137579569, 9.19980017095141, 0.9099896136155886, 17.19828856728441, 0.49640596145648325, 91.5720871040533, 1.2312977080093161)),
        (28, 80, 180, :M, true, false, (1.9780896277183768, 6.548299315703893, 2.2194522694388246, 32.51470318178505, 1.2071294553489746, 309.15317614147835, 0.14120655867953985)),
        (28, 80, 180, :M, true, true, (1.826009572456856, 6.548299315703893, 2.2194522694388246, 32.51470318178505, 0.9031167796633479, 209.9715551212636, 0.14120655867953985)),
        (28, 80, 180, :M, false, false, (1.9780896277183768, 9.293761905023166, 1.510337959412125, 32.51470318178505, 1.2071294553489746, 309.15317614147835, 1.199288580565955)),
        (28, 80, 180, :M, false, true, (1.826009572456856, 9.293761905023166, 1.510337959412125, 32.51470318178505, 0.9031167796633479, 209.9715551212636, 1.199288580565955)),
        (81, 64, 168, :M, true, false, (1.6732601030874201, 6.0978454630928445, 0.9738745465035903, 11.358880277464147, 1.08034392329858, 256.62543770646465, 0.1493077573570792)),
        (81, 64, 168, :M, true, true, (1.327569186919631, 6.0978454630928445, 0.9738745465035903, 11.358880277464147, 0.4666951528638212, 83.8033498442706, 0.1493077573570792)),
        (81, 64, 168, :M, false, false, (1.6732601030874201, 9.073765971885559, 0.6627219316870213, 11.358880277464147, 1.08034392329858, 256.62543770646465, 1.2680932816628645)),
        (81, 64, 168, :M, false, true, (1.327569186919631, 9.073765971885559, 0.6627219316870213, 11.358880277464147, 0.4666951528638212, 83.8033498442706, 1.2680932816628645)),
        (35, 70, 170, :M, true, false, (1.7895836588902, 6.283078076682199, 1.8309210237746485, 25.5013145036879, 1.1085424008536, 272.8166615043603, 0.146)),
        (35, 70, 170, :M, true, true, (1.619285252944501, 6.283078076682199, 1.8309210237746485, 25.5013145036879, 0.7713294071769037, 168.21097582206153, 0.146)),
        (35, 70, 170, :M, false, false, (1.7895836588902, 9.17173087918622, 1.2459423259378106, 25.5013145036879, 1.1085424008536, 272.8166615043603, 1.24)),
        (35, 70, 170, :M, false, true, (1.619285252944501, 9.17173087918622, 1.2459423259378106, 25.5013145036879, 0.7713294071769037, 168.21097582206153, 1.24)),
        (45, 110, 165, :F, true, false, (2.947719565618916, 7.122257314348119, 2.264016593436523, 34.27393400562731, 1.1314861011932773, 277.4829356373051, 0.13040048134380294)),
        (45, 110, 165, :F, true, true, (2.592085574522928, 7.122257314348119, 2.264016593436523, 34.27393400562731, 0.7097955220721459, 149.0098540163805, 0.13040048134380294)),
        (45, 110, 165, :F, false, false, (2.947719565618916, 9.484318304425138, 1.5406639957482353, 34.27393400562731, 1.1314861011932773, 277.4829356373051, 1.107510937440518)),
        (45, 110, 165, :F, false, true, (2.592085574522928, 9.484318304425138, 1.5406639957482353, 34.27393400562731, 0.7097955220721459, 149.0098540163805, 1.107510937440518))
]

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
        pac = Paciente(idade = 53, peso = 77, altura = 177, sexo = :M, pas = 125, pad = 75, pam = 92, fc = 70)
        @testset "Schnider no paciente de referência" begin
            pk = parametros_pk(pac, :schnider)
            @test pk.V1 == 4.27
            @test pk.V2 ≈ 18.9
            @test pk.Q2 ≈ 1.29
            @test pk.CL ≈ 1.89 - 0.0681 * (massa_magra(pac) - 59) rtol = 1e-12
            # fora da faixa das equações lineares de Schnider: erro explícito, sem pisos
            velho = Paciente(idade = 105, peso = 60, altura = 160, sexo = :F, pas = 130, pad = 70, pam = 90, fc = 70)
            @test_throws ArgumentError parametros_pk(velho, :schnider)
        end
        @testset "Eleveld = Python Anesthesia Simulator (12 casos)" begin
            for (idade, peso, altura, sexo, art, op, ref) in ELEVELD_PAS
                p = Paciente(idade = idade, peso = peso, altura = altura, sexo = sexo, pas = 130, pad = 70, pam = 90, fc = 70)
                e = parametros_eleveld(p; arterial = art, opioides = op)
                @test collect((e.CL, e.V1, e.Q2, e.V2, e.Q3, e.V3, e.ke0)) ≈ collect(ref) rtol = 1e-12
            end
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
            # Tabela 3 de Eleveld 2018: γ = 1,89 abaixo do Ce50 e 1,47 acima
            @test bis_eleveld(3.08 / 2, 35) ≈ 93 - 93 * 0.5^1.89 / (1 + 0.5^1.89)
            @test bis_eleveld(3.08 * 2, 35) ≈ 93 - 93 * 2^1.47 / (1 + 2^1.47)
            @test atraso_bis(35) ≈ 15 + exp(0.0517 * 35)
        end
        @testset "BIS atrasado em relação ao Ce de Eleveld" begin
            r = simular_propofol(pac, [bolus(1.0, 2 * pac.peso)]; tfim = 10.0)
            d = atraso_bis(pac.idade) / 60
            i = findfirst(t -> t >= 5.0, r.t)
            @test r.bis[i] ≈ bis_eleveld(PropofolWindkessel.interp(r.t, r.Ce_bis, r.t[i] - d), pac.idade) rtol = 1e-6
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
                        asa = 3, pvi = 20.0, ivci = 50.0)
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

    @testset "Risco pela PAM mínima frente a 65 e 55 mmHg" begin
        for d in (0.5, 1.0, 2.0, 3.0)
            mt = metricas(simular_inducao(pv; dose_mgkg = d))
            esperado = mt.pam_nadir < 55 ? "ALTO" : mt.pam_nadir < 65 ? "MODERADO" : "BAIXO"
            @test mt.risco == esperado
        end
    end

    @testset "Queda relativa da PAM não depende do DC basal estimado" begin
        quedas = map((4.0, 5.5, 7.0)) do dc
            p = Paciente(idade = 53, peso = 77, altura = 177, sexo = :M, pas = 125, pad = 75, pam = 92, fc = 70, dc = dc)
            metricas(simular_inducao(calibrar(p); dose_mgkg = 2.0)).queda_pct
        end
        @test maximum(quedas) - minimum(quedas) < 0.5      # pontos percentuais
    end

    @testset "Monotonicidade dose–resposta" begin
        nadires = [metricas(simular_inducao(pv; dose_mgkg = d, manutencao_ugkgmin = 0.0)).pam_nadir for d in (1.0, 2.0, 3.0)]
        @test issorted(nadires; rev = true)
    end

    @testset "Incerteza (experimental, fora da plataforma)" begin
        inc = previsao_com_incerteza(pv; n = 16, dose_mgkg = 2.0)
        @test all(inc.p05 .<= inc.mediana .<= inc.p95)
        @test 0 <= inc.prob_pam55 <= inc.prob_pam65 <= 1
    end

    @testset "Plataforma: API (sem abrir servidor)" begin
        corpo = Dict("paciente" => Dict("idade" => "53", "sexo" => "M", "peso" => "77", "altura" => "177",
                                        "pas" => "125", "pad" => "75", "pam" => "92", "fc" => "70", "asa" => "2", "pvi" => "12"),
                     "protocolo" => Dict("dose_mgkg" => "2,0", "duracao_bolus_s" => "30", "modelo_pk" => "schnider"))
        r = PropofolWindkessel._api_simular(corpo)
        @test r["metricas"]["pam_nadir"] < r["metricas"]["pam_basal"]
        @test length(r["batimentos"]["t"]) == length(r["batimentos"]["pam"])
        d = PropofolWindkessel._api_doses(corpo)
        @test length(d["doses"]) == 11
        @test issorted([x["pam_nadir"] for x in d["doses"]]; rev = true)
        ruim = Dict("paciente" => merge(corpo["paciente"], Dict("pas" => "70")))
        @test_throws ArgumentError PropofolWindkessel._ler_paciente(ruim)
        sem_pam = Dict("paciente" => filter(kv -> kv[1] != "pam", corpo["paciente"]))
        @test_throws KeyError PropofolWindkessel._ler_paciente(sem_pam)
        e = PropofolWindkessel._api_simular(merge(corpo, Dict("protocolo" => Dict("dose_mgkg" => 2, "modelo_pk" => "eleveld"))))
        @test e["virtual"]["modelo_pk"] == "Eleveld"
    end

    @testset "docs/validacao.md: links de linha apontam para linhas existentes e não vazias" begin
        raiz = pkgdir(PropofolWindkessel)
        doc = read(joinpath(raiz, "docs", "validacao.md"), String)
        links = collect(eachmatch(r"\]\(\.\./([^)#?]+)(?:\?plain=1)?#L(\d+)(?:-L(\d+))?\)", doc))
        @test length(links) > 30
        for m in links
            arq = joinpath(raiz, m.captures[1])
            a = parse(Int, m.captures[2]); b = m.captures[3] === nothing ? a : parse(Int, m.captures[3])
            linhas = isfile(arq) ? readlines(arq) : String[]
            ok = isfile(arq) && 1 <= a <= b <= length(linhas) && !all(isempty ∘ strip, linhas[a:b])
            ok || @info "link quebrado" m.match
            @test ok
        end
    end
end
