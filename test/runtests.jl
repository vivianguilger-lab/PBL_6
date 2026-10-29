# =============================================================================
# Testes de VERIFICAÇÃO ("resolvi as equações corretamente?") e de
# plausibilidade. Inclui os testes do PropofolPBL do grupo (Windkessel do
# windkessel-id, PK/PD do NeoPKPD e acoplamento).
#     julia --project=. -e "using Pkg; Pkg.test()"
# =============================================================================
# -----------------------------------------------------------------------------
# GUIA RÁPIDO PARA QUEM NÃO PROGRAMA
# Este arquivo roda sozinho no GitHub Actions a cada envio de código.
# Cada teste executa uma parte do modelo e compara o resultado com um valor
# que sabemos ser o certo. No fim, o Julia mostra um resumo: quantos testes
# passaram ("Pass") e quantos falharam ("Fail" ou "Error") em cada grupo.
#
# Há dois tipos de teste:
# - VERIFICAÇÃO: "resolvi as equações certo?". Compara o programa com uma
#   resposta exata (fórmula matemática) ou com outro programa independente.
# - PLAUSIBILIDADE: "o resultado faz sentido clínico?". Confere faixas
#   esperadas (ex.: a PAM cai depois do bolus).
# Nenhum dos dois é VALIDAÇÃO com dados de pacientes reais (ainda não feita;
# ver docs/validacao.md, seções 6 e 7).
#
# Sintaxe usada (só o necessário para ler):
#   @testset "nome" begin ... end   agrupa testes sob um nome (aparece no resumo).
#   @test condição                  passa se a condição for verdadeira.
#   a ≈ b                           "a é aproximadamente igual a b". Sem rtol,
#                                   aceita só diferenças mínimas de arredondamento.
#   a ≈ b rtol = 0.01               aceita erro relativo de até 1% (rtol = 0.01);
#                                   rtol = 1e-12 quer dizer 0,0000000001%.
#   @test_throws TipoDeErro f(...)  passa se f(...) PARAR com esse erro
#                                   (usado para conferir que dados absurdos são recusados).
#   ==  igual exatamente;   &&  "e";   .>  .==  .-  comparam/subtraem item a item.
#   θ[1] = primeiro item da lista θ;  1e-3 = 0,001;  # inicia um comentário.
#
# Como rodar (no terminal, dentro da pasta do projeto):
#   julia --project=. -e "using Pkg; Pkg.test()"
# Com aspas duplas funciona no Windows, no macOS e no Linux (a linha 5 acima usa
# aspas simples, que não funcionam no prompt do Windows). Demora alguns minutos.
#
# Se um teste falhar, o resumo mostra o nome do grupo, a linha deste arquivo e
# os valores obtido e esperado. O comentário acima de cada teste explica o que
# a falha significa. Mais contexto: docs/validacao.md (seção 5) e
# docs/modelo_matematico.md.
# -----------------------------------------------------------------------------
# Carrega o nosso pacote (o modelo) e as bibliotecas de álgebra linear,
# estatística (média) e de testes.
using PropofolWindkessel
using LinearAlgebra
using Statistics
using Test

# Passo de tempo do Windkessel: 0,005 s = 5 ms (o mesmo usado no modelo).
const h = 0.005
# Paciente de referência de Schnider (53 anos, 77 kg, 177 cm, homem): é o
# "paciente médio" das equações de Schnider, em que as correções por idade,
# peso e altura valem zero. Usado nos testes do paciente virtual, mais abaixo.
const REF = Paciente(nome = "referência Schnider", idade = 53, peso = 77, altura = 177, sexo = :M,
                     pas = 125, pad = 75, pam = 92, fc = 70)

# TABELA DE REFERÊNCIA DO ELEVELD (usada no teste "Eleveld = Python Anesthesia
# Simulator", mais abaixo). Estes números NÃO foram digitados à mão nem copiados
# do nosso código: foram gerados rodando o código do Eleveld de outro programa,
# independente do nosso (o Python Anesthesia Simulator, PAS). Se o nosso código
# der os mesmos números, as duas implementações concordam.
# Cada linha é um caso: 5 pacientes × amostra arterial/venosa × com/sem opioide.
# true = sim, false = não; :F = feminino, :M = masculino.
# Unidades: CL, Q2 e Q3 em L/min; V1, V2 e V3 em L; ke0 em 1/min.
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

# Grupo principal: reúne todos os testes abaixo. Se qualquer teste interno
# falhar, este grupo também aparece como falho no resumo final.
@testset "PropofolWindkessel" begin

    # ---------------------------------------------------------------- grupo
    # Testes vindos do repositório do grupo (PropofolPBL), mantidos iguais.
    # WINDKESSEL: o "circuito elétrico" da circulação arterial, com 4 parâmetros
    # θ = [Rp, C, Rc, L]: resistência periférica, complacência arterial, resistência
    # característica da aorta e inertância do sangue (valores de Stergiopulos 1999).
    # Se algum destes falhar, a parte de pressão arterial do modelo não é confiável.
    @testset "Windkessel (windkessel-id)" begin
        # θ = parâmetros publicados; q = um batimento de fluxo aórtico medido
        # (Stergiopulos 1999, um ponto a cada 5 ms); u = 5 batimentos seguidos.
        θ, _ = params_stergiopulos1999total()
        _, _, q = dados_stergiopulos1999()
        u = repeat(q, 5)

        # VERIFICAÇÃO. solvex0 calcula o estado inicial do "regime estável": aquele em
        # que cada batimento termina exatamente como começou (paciente estável,
        # batimento após batimento). O teste simula um batimento a partir desse estado
        # e confere que ele volta ao mesmo ponto (erro relativo < 1e-8).
        # Falha = a simulação começaria num estado transitório e a PA basal sairia distorcida.
        @testset "solvex0 dá estado periódico" begin
            sys = G(θ, h)
            x0 = solvex0(sys, q)
            @test estado_final(sys, q, x0) ≈ x0 rtol = 1e-8
        end
        # VERIFICAÇÃO por uma propriedade exata: na média, o Windkessel obedece à
        # "lei de Ohm" da circulação, PAM = Rp × fluxo médio (como PAM = RVP × DC).
        # θ[1] é o Rp. A tolerância de 2% é só uma folga.
        # Falha = erro nas equações ou na simulação; a PAM simulada não seria confiável.
        @testset "PAM do modelo ≈ Rp · fluxo médio" begin
            y = simulate_wk(θ, solvex0(G(θ, h), q), q, h)
            @test mean(y) ≈ θ[1] * mean(q) rtol = 0.02
        end
        # VERIFICAÇÃO por propriedades exatas da impedância Z (a "oposição" do sistema
        # arterial ao fluxo, em cada frequência). Com fluxo constante (frequência zero),
        # |Z| = Rp (θ[1]); em frequência muito alta, |Z| tende a Rc (θ[3]).
        # As frequências 1e-6 e 1e6 fazem o papel de "zero" e "infinito".
        # O mesmo vale para o Windkessel de 3 elementos (G3, que usa só θ[1:3]).
        # Falha = o circuito foi montado errado.
        @testset "bode: |Z(0)| = Rp e |Z(∞)| → Rc" begin
            mag, _ = bode_wk(θ, [1e-6, 1e6])
            @test mag[1] ≈ θ[1] rtol = 1e-4
            @test mag[2] ≈ θ[3] rtol = 1e-3
            mag3, _ = bode_wk(θ[1:3], [1e6]; G_ = G3)
            @test mag3[1] ≈ θ[3] rtol = 1e-3
        end
        # VERIFICAÇÃO da identificação (ajuste dos parâmetros a partir de medidas):
        # 1) gera uma curva de pressão "sintética" com os θ conhecidos;
        # 2) começa o ajuste com θ errados de propósito (C +40%, Rc −30%, L +50%);
        # 3) confere que o ajuste reencontra os θ verdadeiros (erro < 0,1%).
        # Falha = o ajuste não seria confiável para identificar o Windkessel de um paciente.
        @testset "identify recupera parâmetros de dados sintéticos" begin
            y = simulate_wk(θ, solvex0(G(θ, h), u), u, h)
            θ0 = θ .* [1.0, 1.4, 0.7, 1.5]
            x00 = solvex0(G(θ0, h), u)
            phat, res = identify(θ0, x00, u, y, h; x0id = false, x0solve = true)
            θhat, _ = deal(phat, IdConditions(θ0, x00, [1, 1, 1, 1], false, true))
            @test θhat ≈ θ rtol = 1e-3
        end
        # Lê os dados in vivo que vêm com o windkessel-id (4 arquivos de pressão e 4 de
        # fluxo aórticos medidos). Confere que há 4 de cada e que a pressão média está
        # numa faixa fisiológica (50–150 mmHg). Falha = arquivo ausente ou lido errado.
        @testset "importPQ lê os dados in vivo" begin
            P, Q = importPQ(joinpath(WINDKESSEL_ID_DIR, "invivo_data", "Ao_baseline1_Ao"), 4)
            @test length(P) == 4 && length(Q) == 4
            @test 50 < mean(P[1].vb_va) < 150
        end
    end

    # PROPOFOL (NeoPKPD): farmacocinética (concentração no plasma, Cp, e no sítio
    # efetor, Ce) e BIS. Se estes falharem, a cadeia dose → concentração → efeito
    # do modelo está errada.
    @testset "Propofol (NeoPKPD)" begin
        # pac = paciente de referência de Schnider (53 anos, 77 kg, 177 cm, homem).
        pac = Paciente(idade = 53, peso = 77, altura = 177, sexo = :M, pas = 125, pad = 75, pam = 92, fc = 70)
        # VERIFICAÇÃO contra a fórmula publicada de Schnider. No paciente de referência,
        # V1, V2 e Q2 são exatamente as constantes da tabela (4,27 L; 18,9 L; 1,29 L/min)
        # e o CL só tem a correção pela massa magra (fórmula reescrita aqui à mão).
        # Falha = parâmetro de Schnider digitado ou calculado errado.
        @testset "Schnider no paciente de referência" begin
            pk = parametros_pk(pac, :schnider)
            @test pk.V1 == 4.27
            @test pk.V2 ≈ 18.9
            @test pk.Q2 ≈ 1.29
            @test pk.CL ≈ 1.89 - 0.0681 * (massa_magra(pac) - 59) rtol = 1e-12
            # Segunda parte: paciente de 105 anos. Nessa idade, a fórmula de Schnider daria
            # V2 negativo (sem sentido físico). O programa deve PARAR com ArgumentError em
            # vez de inventar um valor. Falha = o modelo aceitaria um caso fora da sua faixa.
            # fora da faixa das equações lineares de Schnider: erro explícito, sem pisos
            velho = Paciente(idade = 105, peso = 60, altura = 160, sexo = :F, pas = 130, pad = 70, pam = 90, fc = 70)
            @test_throws ArgumentError parametros_pk(velho, :schnider)
        end
        # VERIFICAÇÃO independente: para cada linha da tabela ELEVELD_PAS (início do
        # arquivo), calcula os 7 parâmetros de Eleveld com o nosso código e compara com
        # os do Python Anesthesia Simulator (erro relativo < 1e-12: iguais até cerca da
        # 12ª casa significativa). São 20 casos: 5 pacientes × amostra arterial/venosa ×
        # com/sem opioide. Falha = nossa implementação do Eleveld difere da
        # referência (e o BIS, que usa o Eleveld, também ficaria errado).
        @testset "Eleveld = Python Anesthesia Simulator (20 casos)" begin
            for (idade, peso, altura, sexo, art, op, ref) in ELEVELD_PAS
                # PAS, PAD, PAM e FC não entram no PK: são só valores para criar o paciente.
                p = Paciente(idade = idade, peso = peso, altura = altura, sexo = sexo, pas = 130, pad = 70, pam = 90, fc = 70)
                e = parametros_eleveld(p; arterial = art, opioides = op)
                @test collect((e.CL, e.V1, e.Q2, e.V2, e.Q3, e.V3, e.ke0)) ≈ collect(ref) rtol = 1e-12
            end
        end
        # PLAUSIBILIDADE de um bolus de 2 mg/kg (dose habitual de indução no adulto):
        # - o pico de Ce ocorre entre 1 e 4 min (esperado ~1,5–2 min com Schnider);
        # - o Cp atinge o pico antes do Ce (o fármaco chega ao sítio efetor com atraso);
        # - o BIS cai abaixo de 60 (hipnose) e começa em 93 (acordado; E0 de Eleveld).
        # Falha = curva de concentração ou de BIS fora do esperado clinicamente.
        @testset "bolus de 2 mg/kg" begin
            r = simular_propofol(pac, [bolus(0.0, 2 * pac.peso)]; tfim = 20.0)
            ipico = argmax(r.Ce)
            @test 1.0 < r.t[ipico] < 4.0          # pico no sítio efetor em ~1,5–2 min (Schnider)
            @test argmax(r.Cp) < ipico            # histerese: Cp atinge o pico antes de Ce
            @test minimum(r.bis) < 60
            @test r.bis[1] ≈ 93.0
        end
        # VERIFICAÇÃO contra a solução exata. Num bolus instantâneo, o modelo de
        # 3 compartimentos tem resposta matemática exata: A(t) = exp(K·t)·A0, em que
        # A0 = [dose, 0, 0] (toda a dose no plasma) e K (matriz_pk) reúne as constantes
        # de transferência entre compartimentos. Cp = quantidade no plasma / V1.
        # D = 154 mg (2 mg/kg × 77 kg). Conferido em 1, 2, 5 e 10 min, erro < 0,1%.
        # Falha = o NeoPKPD (ou o nosso uso dele) não resolve as equações corretamente.
        @testset "NeoPKPD = solução exata exp(K·t)·A0 (verificação)" begin
            pk = parametros_pk(pac, :schnider)
            D = 154.0
            r = simular_propofol(pac, [bolus(0.0, D; duracao = 0.0)]; tfim = 10.0)
            K = matriz_pk(pk)
            for t in (1.0, 2.0, 5.0, 10.0)
                # i = instante da simulação mais próximo de t.
                i = argmin(abs.(r.t .- t))
                @test r.Cp[i] ≈ (exp(K * r.t[i]) * [D, 0, 0])[1] / pk.V1 rtol = 1e-3
            end
        end
        # VERIFICAÇÃO contra a Tabela 3 de Eleveld 2018 (curva BIS × Ce):
        # - aos 35 anos, Ce50 = 3,08 µg/ml: nesse Ce o BIS é metade de 93 (E0);
        # - sem propofol (Ce = 0), BIS = 93;
        # - com o mesmo Ce, o idoso (80 anos) tem BIS menor (mais sensível);
        # - a inclinação γ muda acima e abaixo do Ce50 (conferido em Ce50/2 e 2·Ce50);
        # - atraso do BIS = 15 + e^(0,0517·idade) segundos.
        # Falha = a curva do BIS do modelo difere da publicada.
        @testset "BIS de Eleveld" begin
            @test bis_eleveld(3.08, 35) ≈ 93 / 2
            @test bis_eleveld(0.0, 35) == 93.0
            @test bis_eleveld(2.0, 80) < bis_eleveld(2.0, 20)   # idoso mais sensível
            # Tabela 3 de Eleveld 2018: γ = 1,89 abaixo do Ce50 e 1,47 acima
            @test bis_eleveld(3.08 / 2, 35) ≈ 93 - 93 * 0.5^1.89 / (1 + 0.5^1.89)
            @test bis_eleveld(3.08 * 2, 35) ≈ 93 - 93 * 2^1.47 / (1 + 2^1.47)
            @test atraso_bis(35) ≈ 15 + exp(0.0517 * 35)
        end
        # VERIFICAÇÃO: o BIS de cada instante é o BIS do Ce de Eleveld (Ce_bis) de
        # d minutos antes, porque o monitor de BIS tem atraso. Confere aos 5 min, após
        # um bolus no minuto 1. Falha = atraso aplicado errado ou não aplicado.
        @testset "BIS atrasado em relação ao Ce de Eleveld" begin
            r = simular_propofol(pac, [bolus(1.0, 2 * pac.peso)]; tfim = 10.0)
            # d = atraso do BIS, convertido de segundos para minutos.
            d = atraso_bis(pac.idade) / 60
            i = findfirst(t -> t >= 5.0, r.t)
            @test r.bis[i] ≈ bis_eleveld(PropofolWindkessel.interp(r.t, r.Ce_bis, r.t[i] - d), pac.idade) rtol = 1e-6
        end
    end

    # ACOPLAMENTO (teste do grupo): o Ce do propofol reduz Rp (vasodilatação), o
    # volume sistólico e a FC, e o Windkessel calcula a pressão batimento a batimento.
    # PLAUSIBILIDADE: antes do fármaco, a PAM simulada reproduz a da paciente (±3%);
    # após a indução, a PAM cai mais de 15%; em todo batimento PAS > PAM > PAD;
    # a onda de pressão tem um valor para cada instante.
    # Falha = o acoplamento não gera a hipotensão esperada ou gera pressões impossíveis.
    @testset "Acoplamento PK/PD → Windkessel (grupo)" begin
        pac = Paciente(idade = 45, peso = 70, altura = 170, sexo = :F, pas = 120, pad = 75, pam = 90.0, fc = 75.0)
        pkpd, bat, onda = simular_paciente(pac; tfim = 10.0, forma_de_onda = true)
        @test bat.PAM[1] ≈ pam_basal(pac) rtol = 0.03     # basal reproduz a PAM do paciente
        @test minimum(bat.PAM) < 0.85 * pam_basal(pac)     # hipotensão após a indução
        @test all(bat.PAS .> bat.PAM .> bat.PAD)
        @test length(onda.pressao) == length(onda.t_s)
        # pac2: mesma paciente, mas com 75 anos, ASA III, PVI 20% e IVCI 50% (sinais de
        # maior risco). O fator de sensibilidade deve ser maior e a PAM mínima, mais baixa.
        # Falha = os fatores de risco não pioram a resposta, como deveriam.
        pac2 = Paciente(idade = 75, peso = 70, altura = 170, sexo = :F, pas = 120, pad = 75, pam = 90.0, fc = 75.0,
                        asa = 3, pvi = 20.0, ivci = 50.0)
        @test sensibilidade_hemodinamica(pac2) > sensibilidade_hemodinamica(pac)
        _, bat2, _ = simular_paciente(pac2; tfim = 10.0)
        @test minimum(bat2.PAM) < minimum(bat.PAM)
    end

    # ---------------------------------------------------------------- extensões
    # Testes das EXTENSÕES deste repositório (paciente virtual, risco, plataforma).
    # pv = paciente virtual do paciente de referência (REF): Windkessel ajustado
    # para reproduzir a PAM e a pressão de pulso (PAS − PAD) dele.
    pv = calibrar(REF)

    # VERIFICAÇÃO da calibração: o ajuste convergiu; a PAM e a pressão de pulso
    # simuladas batem com as do monitor (±1%); e a constante de tempo diastólica
    # (Rp·C, em segundos) fica numa faixa fisiológica (0,3–5 s).
    # Falha = o paciente virtual não representa o paciente real no basal.
    @testset "Paciente virtual: PAM e pressão de pulso do monitor" begin
        @test pv.convergiu
        @test pv.obtido.PAM ≈ pam_basal(REF) rtol = 0.01
        @test pv.obtido.PP ≈ REF.pas - REF.pad rtol = 0.01
        @test 0.3 < tau_diastolico(pv) < 5.0
    end

    # Duas simulações usadas nos testes abaixo: r0 = sem propofol (dose zero, 5 min);
    # r2 = indução padrão de 2 mg/kg (bolus em 30 s + manutenção de 100 µg/kg/min).
    r0 = simular_inducao(pv; dose_mgkg = 0.0, manutencao_ugkgmin = 0.0, duracao = 5.0)
    r2 = simular_inducao(pv; dose_mgkg = 2.0)

    # VERIFICAÇÃO: sem fármaco, nada deve mudar. A PAM fica a menos de 0,5 mmHg da
    # basal e a concentração no plasma é zero o tempo todo.
    # Falha = o modelo "deriva" sozinho e criaria hipotensão sem causa.
    @testset "Linha de base estável sem fármaco" begin
        @test maximum(abs.(r0.bat.PAM .- pv.obtido.PAM)) < 0.5
        @test all(r0.pkpd.Cp .== 0)
    end

    # PLAUSIBILIDADE clínica (não é validação com pacientes): com 2 mg/kg no adulto
    # de referência, a PAM cai entre 10% e 45%; o valor mínimo ocorre entre 0,5 e
    # 15 min após o bolus; o BIS mínimo fica abaixo de 60; e a hipnose (BIS < 60)
    # chega em menos de 2,5 min. Falha = resultado fora da faixa clínica esperada.
    @testset "Plausibilidade: 2 mg/kg em adulto de referência" begin
        mt = metricas(r2)
        @test 10 < mt.queda_pct < 45
        @test 0.5 < mt.t_nadir < 15
        @test mt.bis_min < 60
        @test mt.t_hipnose < 2.5
    end

    # VERIFICAÇÃO da regra de risco: para 4 doses, a classificação depende só da
    # PAM mínima: < 55 mmHg = ALTO; 55–65 = MODERADO; ≥ 65 = BAIXO (Walsh 2013;
    # Salmasi 2017). Falha = a plataforma mostraria o risco errado.
    @testset "Risco pela PAM mínima frente a 65 e 55 mmHg" begin
        for d in (0.5, 1.0, 2.0, 3.0)
            mt = metricas(simular_inducao(pv; dose_mgkg = d))
            # A linha abaixo reescreve a regra. "cond ? a : b" quer dizer "se cond, a; senão, b".
            esperado = mt.pam_nadir < 55 ? "ALTO" : mt.pam_nadir < 65 ? "MODERADO" : "BAIXO"
            @test mt.risco == esperado
        end
    end

    # Quando o débito cardíaco (DC) não é medido, o modelo o estima (índice cardíaco
    # 3,0 L/min/m²). Este teste troca o DC basal (4,0; 5,5 e 7,0 L/min) e confere que
    # a QUEDA percentual da PAM com 2 mg/kg quase não muda (diferença < 0,5 ponto
    # percentual). Ou seja, a estimativa do DC não altera o resultado principal.
    # Falha = a previsão dependeria de um valor que não medimos.
    @testset "Queda relativa da PAM não depende do DC basal estimado" begin
        quedas = map((4.0, 5.5, 7.0)) do dc
            p = Paciente(idade = 53, peso = 77, altura = 177, sexo = :M, pas = 125, pad = 75, pam = 92, fc = 70, dc = dc)
            metricas(simular_inducao(calibrar(p); dose_mgkg = 2.0)).queda_pct
        end
        @test maximum(quedas) - minimum(quedas) < 0.5      # pontos percentuais
    end

    # PLAUSIBILIDADE: mais dose, PAM mínima mais baixa (1, 2 e 3 mg/kg, sem
    # manutenção). issorted(...; rev = true) confere que a lista está em ordem
    # decrescente. Falha = o modelo daria menos hipotensão com mais propofol.
    @testset "Monotonicidade dose–resposta" begin
        nadires = [metricas(simular_inducao(pv; dose_mgkg = d, manutencao_ugkgmin = 0.0)).pam_nadir for d in (1.0, 2.0, 3.0)]
        @test issorted(nadires; rev = true)
    end

    # Faixa de incerteza (EXPERIMENTAL, fora da plataforma até as variabilidades
    # serem validadas; docs/validacao.md, seção 6). Confere só a coerência interna:
    # percentil 5 ≤ mediana ≤ percentil 95, probabilidades entre 0 e 1 e
    # P(PAM < 55) ≤ P(PAM < 65). São só 16 simulações, para o teste ser rápido.
    @testset "Incerteza (experimental, fora da plataforma)" begin
        inc = previsao_com_incerteza(pv; n = 16, dose_mgkg = 2.0)
        @test all(inc.p05 .<= inc.mediana .<= inc.p95)
        @test 0 <= inc.prob_pam55 <= inc.prob_pam65 <= 1
    end

    # PLATAFORMA WEB: testa as funções que respondem ao navegador, sem abrir o
    # servidor. "corpo" imita o formulário preenchido pelo médico (números como
    # texto, com vírgula decimal). Confere: a PAM mínima fica abaixo da basal; a
    # tabela de doses tem 11 doses (0,5 a 3,0 mg/kg, de 0,25 em 0,25) e a PAM
    # mínima cai com a dose; dados inválidos são recusados com erro.
    # Falha = a tela mostraria resultados errados ou aceitaria dados absurdos.
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
        # PAS de 70 mmHg com PAD de 75 é impossível: deve ser recusado (ArgumentError).
        ruim = Dict("paciente" => merge(corpo["paciente"], Dict("pas" => "70")))
        @test_throws ArgumentError PropofolWindkessel._ler_paciente(ruim)
        # Formulário sem a PAM (campo obrigatório): deve dar KeyError (campo ausente).
        sem_pam = Dict("paciente" => filter(kv -> kv[1] != "pam", corpo["paciente"]))
        @test_throws KeyError PropofolWindkessel._ler_paciente(sem_pam)
        # Também deve funcionar quando o médico escolhe o PK de Eleveld.
        e = PropofolWindkessel._api_simular(merge(corpo, Dict("protocolo" => Dict("dose_mgkg" => 2, "modelo_pk" => "eleveld"))))
        @test e["virtual"]["modelo_pk"] == "Eleveld"
    end

    # RASTREABILIDADE: docs/validacao.md tem links como "propofol.jl#L95" que levam à
    # linha exata do código de cada valor. Este teste lê o documento, acha todos os
    # links (mais de 30) e confere que cada linha citada existe e não está vazia.
    # Falha = link quebrado: o código mudou e o documento ficou desatualizado (o
    # resumo mostra o link, com a mensagem "link quebrado"). Corrija o número no .md.
    # Atenção: acrescentar linhas a ESTE arquivo também desloca os links para ele.
    @testset "docs/validacao.md: links de linha apontam para linhas existentes e não vazias" begin
        raiz = pkgdir(PropofolWindkessel)
        doc = join(read.(joinpath.(raiz, "docs", ["validacao.md"]), String))
        # Expressão regular que reconhece links como (../src/arquivo.jl#L10-L12).
        links = collect(eachmatch(r"\]\(\.\./([^)#?]+)(?:\?plain=1)?#L(\d+)(?:-L(\d+))?\)", doc))
        @test length(links) > 30
        for m in links
            arq = joinpath(raiz, m.captures[1])
            a = parse(Int, m.captures[2]); b = m.captures[3] === nothing ? a : parse(Int, m.captures[3])
            linhas = isfile(arq) ? readlines(arq) : String[]
            # ok = o arquivo existe, as linhas estão dentro dele e não são todas em branco.
            ok = isfile(arq) && 1 <= a <= b <= length(linhas) && !all(isempty ∘ strip, linhas[a:b])
            ok || @info "link quebrado" m.match
            @test ok
        end
    end
end
