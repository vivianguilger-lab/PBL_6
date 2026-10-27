# =============================================================================
# calibracao.jl — paciente virtual calibrado com os dados pré-indução
# =============================================================================
#
# Base: theta_basal do grupo — Windkessel de 4 elementos com Rp = PAM/DC e
# C, Rc, L de Stergiopulos (1999). Como o monitor fornece também PAS e PAD,
# identificamos a complacência arterial C do próprio paciente para que a
# pressão de pulso do modelo seja a medida (mesma ideia de identificação do
# windkessel-id, aqui com um único parâmetro e uma única medida).
#
#   PAM = Rp · Q̄            (Z(0) = Rp)        → Rp
#   PP  = PAS − PAD  ↓  quando C ↑              → C (bisseção em log C)
# =============================================================================
#
# Este arquivo cria o "paciente virtual" (o Windkessel ajustado aos dados que o
# monitor mostra antes da indução) e simula a indução nesse paciente.
#   1. Rp = PAM/DC garante que o modelo reproduz a PAM medida.
#   2. C (complacência arterial) é procurada para que a pressão de pulso do
#      modelo (PP = PAS − PAD) seja a medida. Artérias mais rígidas (C menor,
#      como no idoso) → pulso mais largo.
#   3. Rc e L ficam nos valores de Stergiopulos 1999.
# Só PAM e PP são reproduzidas exatamente; PAS e PAD separadas não, porque a
# forma da onda de fluxo é fixa (README, Limitações).
# Unidades: mmHg, ml, s. Rastreabilidade em docs/validacao.md.

# Faixa de busca da complacência C, em ml/mmHg (Stergiopulos: 1,22). A faixa
# 0,2–6 é SUPOSIÇÃO nossa, sem fonte (docs/validacao.md, seção 6).

const C_MIN, C_MAX = 0.2, 6.0   # ml/mmHg 
#VALORES ASSUMIDOS, DEVEM SER MUDADOS PARA VALORES VERIFICADOS!

# -----------------------------------------------------------------------------
# Estrutura PacienteVirtual: tudo o que define o paciente simulado
# -----------------------------------------------------------------------------
# Campos: os biomarcadores originais (`paciente`); θ = [Rp, C, Rc, L] já
# calibrado; os parâmetros de PK (`pk`: CL, V1, V2, V3, Q2, Q3, ke0 do modelo
# escolhido); o efeito do propofol; a sensibilidade s; `alvo` (o que foi medido),
# `obtido` (o que o modelo reproduz) e `convergiu` (true se a pressão de pulso
# medida foi alcançada dentro da faixa de C).
# `NamedTuple` = registro leve com campos nomeados, ex.: (PAS = 148, PAD = 82),
# lido como alvo.PAS. `Bool` = verdadeiro/falso (true/false).
"""
Paciente virtual: biomarcadores + Windkessel basal identificado + PK + efeito.
"""
struct PacienteVirtual
    paciente::Paciente
    θ::Vector{Float64}             # [Rp, C, Rc, L] (mmHg, ml, s)
    pk::NamedTuple
    efeito::EfeitoHemodinamico
    s::Float64                     # sensibilidade hemodinâmica
    alvo::NamedTuple
    obtido::NamedTuple
    convergiu::Bool
end

# τ (tau) = Rp·C, em s: constante de tempo da queda exponencial da pressão na
# diástole. Ex.: Rp 1,1 mmHg·s/ml × C 1,2 ml/mmHg ≈ 1,3 s. θ[1] = Rp, θ[2] = C.
"Constante de tempo diastólica τ = Rp·C [s]."
tau_diastolico(pv::PacienteVirtual) = pv.θ[1] * pv.θ[2] # exemplo ilustrativo

# Pressão de pulso (mmHg) do batimento basal periódico para um θ, FC (bpm) e VS
# (ml). Devolve a PP (máximo − mínimo da pressão) e a curva de pressão inteira.
# O `_` no início do nome indica função interna (não exportada pelo pacote).
function _pp(θ, fc, vs)
    y, _, _ = batimento_basal(θ, fc, vs)
    return maximum(y) - minimum(y), y
end

# -----------------------------------------------------------------------------
# calibrar: monta o paciente virtual
# -----------------------------------------------------------------------------
# Entrada: `p` (struct Paciente com idade, sexo, peso, altura, PAS, PAD, PAM, FC,
# ASA, PVI, IVCI e, se houver, DC). Opções: modelo de PK (:schnider, :marsh ou
# :eleveld), parâmetros do efeito, desvios individuais da PK (`eta_pk`, usados
# só nas análises experimentais de incerteza; vazio = valores típicos), a
# sensibilidade s e se a complacência C deve ser ajustada (`ajustar_C`).
# Saída: um PacienteVirtual (`-> PacienteVirtual` no texto de ajuda indica o
# tipo devolvido).
"""
    calibrar(p; modelo_pk=:schnider, efeito=EfeitoHemodinamico(), eta_pk=NamedTuple(),
             s=sensibilidade_hemodinamica(p), ajustar_C=true) -> PacienteVirtual
"""
function calibrar(p::Paciente; modelo_pk = :schnider, efeito = EfeitoHemodinamico(), eta_pk = NamedTuple(),
                  s = sensibilidade_hemodinamica(p), ajustar_C = true)
    # θ inicial: Rp = PAM/DC do paciente; C, Rc e L de Stergiopulos (acoplamento.jl).
    # VS basal (ml) = DC (L/min) × 1000 / FC (bpm). debito_basal (src/paciente.jl)
    # usa `something(p.dc, ...)`: o DC medido, se existir; senão a estimativa
    # 3,0 L/min/m² × superfície corporal (Mosteller).
    θ = theta_basal(p)
    vs = debito_basal(p) * 1000 / p.fc
    # `alvo`: os valores medidos que o modelo deve reproduzir; PP = PAS − PAD.
    alvo = (PAS = p.pas, PAD = p.pad, PAM = pam_basal(p), PP = p.pas - p.pad, DC = debito_basal(p), FC = p.fc)
    convergiu = true
    # Ajuste da complacência C (só se ajustar_C = true). Sabemos que a PP cai
    # quando C aumenta (artéria mais elástica amortece o pulso), então basta
    # procurar a C cuja PP seja a medida.
    if ajustar_C
        # `ppC(C)` é uma função local: devolve a PP do modelo para uma C dada,
        # mantendo Rp, Rc e L (o `[1]` pega só a PP, sem a curva).
        ppC(C) = _pp([θ[1], C, θ[3], θ[4]], p.fc, vs)[1]
        # Primeiro confere se a resposta está dentro da faixa [C_MIN, C_MAX].
        # Se não estiver, usa o extremo mais próximo e marca convergiu = false
        # (o paciente virtual avisa que a pressão de pulso não foi atingida).
        # se nenhum C entre 0,2 e 6 consegue reproduzir o PP do paciente, o programa simplesmente coloca C no limite mais próximo   VERIFICAR!
        if ppC(C_MAX) > alvo.PP          # pulso pequeno demais até para a maior complacência 
            θ[2] = C_MAX; convergiu = false
        elseif ppC(C_MIN) < alvo.PP      # pulso grande demais até para a menor complacência
            θ[2] = C_MIN; convergiu = false
        else
    # novamente, o mín e o máx são hipóteses que devem, idealmente, serem trocadas por valores justificados   VERIFICAR!
            
            # Bisseção: começa no intervalo [log C_MIN, log C_MAX] e, a cada volta,
            # testa o ponto do meio c. Se a PP ali ainda é maior que a medida, C
            # precisa ser maior (a = c); senão, menor (b = c).
            # Usar log C dá a mesma precisão relativa em toda a faixa (que vai de
            # 0,2 a 6, 30 vezes). 60 voltas dividem o intervalo por 2⁶⁰: o erro
            # final é desprezível. `for _ in 1:60` repete 60 vezes sem usar o contador.
            a, b = log(C_MIN), log(C_MAX)
            for _ in 1:60
                c = (a + b) / 2
                # `cond ? (a = c) : (b = c)` = "se cond, a = c; senão, b = c".
                ppC(exp(c)) > alvo.PP ? (a = c) : (b = c)
            end
            # C final = centro do último intervalo, de volta de log para ml/mmHg.
            θ[2] = exp((a + b) / 2)
        end
    end
    # Confere o resultado: simula o batimento basal com o θ final e guarda o que
    # o modelo obteve (PAS = máximo, PAD = mínimo, PAM = média, em mmHg).
    # A PAM sai igual à medida (PAM = Rp·Q̄; diferença de décimos de mmHg pela
    # discretização); a PP sai igual à medida se convergiu = true.
    _, y = _pp(θ, p.fc, vs)
    obtido = (PAS = maximum(y), PAD = minimum(y), PAM = mean(y), PP = maximum(y) - minimum(y),
              DC = debito_basal(p), VS = vs, FC = p.fc)
    # Monta o paciente virtual. parametros_pk (src/propofol.jl) calcula a PK do
    # modelo escolhido (Schnider 1998/1999, Marsh 1991 ou Eleveld 2018) para
    # este paciente.
    return PacienteVirtual(p, θ, parametros_pk(p, modelo_pk; eta = eta_pk), efeito, s, alvo, obtido, convergiu)
end

# -----------------------------------------------------------------------------
# Como o paciente virtual aparece quando é mostrado no REPL do Julia
# -----------------------------------------------------------------------------
# `Base.show` é a função que o Julia usa para mostrar um objeto na tela; aqui
# ela ganha uma versão para PacienteVirtual. `@printf` imprime com formato:
# `%.0f` = número sem casas decimais, `%.2f` = duas casas, `%s` = texto.
# Mostra: PA do modelo × PA medida (com aviso se a PP não foi atingida), os 4
# parâmetros do Windkessel e τ, DC, VS, s e os parâmetros de PK.
function Base.show(io::IO, pv::PacienteVirtual)
    Rp, C, Rc, L = pv.θ
    a, o = pv.alvo, pv.obtido
    println(io, "Paciente virtual — ", pv.paciente)
    @printf(io, "  PA do modelo %.0f/%.0f (PAM %.1f) mmHg | medida %.0f/%.0f (PAM %.1f)%s\n",
            o.PAS, o.PAD, o.PAM, a.PAS, a.PAD, a.PAM, pv.convergiu ? "" : "  ⚠ pressão de pulso não atingida")
    @printf(io, "  Windkessel 4 elementos: Rp = %.3f mmHg·s/ml, C = %.2f ml/mmHg, Rc = %.3f, L = %.4f | τ = Rp·C = %.2f s\n",
            Rp, C, Rc, L, Rp * C)
    @printf(io, "  DC %.2f L/min, VS %.0f ml | sensibilidade hemodinâmica s = %.2f\n", o.DC, o.VS, pv.s)
    pk = pv.pk
    @printf(io, "  PK (%s): CL %.2f L/min, V1 %.2f L, V2 %.1f L, V3 %.0f L, ke0 %.3f /min",
            pk.modelo, pk.CL, pk.V1, pk.V2, pk.V3, pk.ke0)
end

# =============================================================================
# Simulação da indução a partir do paciente virtual
# =============================================================================

# Estrutura Resultado: guarda tudo de uma indução simulada.
#   pkpd: tempo (min), Cp e Ce (µg/ml) e BIS, numa grade de 1 s;
#   bat: tabela por batimento (t_min, Ce, BIS, FC, VS, DC, Rp, PAS, PAD, PAM);
#   onda: pressão e fluxo amostra a amostra, ou `nothing` se não foi pedida
#         (`Union{Nothing,NamedTuple}` = "ou nada, ou uma NamedTuple");
#   t_inducao: minuto do bolus (o minuto anterior é a linha de base);
#   pv: o paciente virtual usado; protocolo: dose (mg/kg), duração do bolus (s)
#       e manutenção (µg/kg/min). As métricas (src/metricas.jl) leem daqui.
"""
Resultado de uma indução simulada: PK/PD na grade de 1 s (`pkpd`), valores
por batimento (`bat`) e, opcionalmente, a forma de onda (`onda`). Tempos em min.
"""
struct Resultado
    pkpd::NamedTuple
    bat::NamedTuple
    onda::Union{Nothing,NamedTuple}
    t_inducao::Float64
    pv::PacienteVirtual
    protocolo::NamedTuple
end

# -----------------------------------------------------------------------------
# simular_inducao: a previsão para uma dose
# -----------------------------------------------------------------------------
# Linha do tempo (valores padrão): de 0 a 1 min, linha de base sem fármaco; no
# minuto 1, bolus de 2 mg/kg injetado em 30 s; a partir de 1 min após o bolus,
# manutenção de 100 µg/kg/min até o fim; total de 1 + 20 = 21 min (a janela de
# 20 min após o bolus é a da hipotensão pós-indução). Dose e manutenção padrão
# seguem o esquema do grupo e a bula do Diprivan (docs/validacao.md).                    VERIFICAR ESSA ORIGEM!
# Entradas: o paciente virtual e o protocolo; `pk`, `efeito` e `s` podem ser
# trocados (usado nas análises de sensibilidade e incerteza).
# Saída: um Resultado.
"""
    simular_inducao(pv; dose_mgkg=2.0, duracao_bolus=30.0, manutencao_ugkgmin=100.0,
                    inicio_manutencao=1.0, t_inducao=1.0, duracao=20.0, forma_de_onda=false,
                    pk=pv.pk, efeito=pv.efeito, s=pv.s) -> Resultado

Bolus de `dose_mgkg` injetado em `duracao_bolus` **segundos** no minuto
`t_inducao` (o minuto anterior é a linha de base), seguido de manutenção a partir
de `inicio_manutencao` min após o bolus (esquema do grupo). Simula `duracao` min.
"""
function simular_inducao(pv::PacienteVirtual;             # todos esses valores vem de fontes citadas nos documentos de validação, presentes na file docs
                        dose_mgkg = 2.0,                  # VERIFICAR CADA UM E MUDAR SE NECESSÁRIO!
                        duracao_bolus = 30.0, 
                        manutencao_ugkgmin = 100.0,
                        inicio_manutencao = 1.0, 
                        t_inducao = 1.0, 
                        duracao = 20.0, 
                        forma_de_onda = false,
                        pk = pv.pk, 
                        efeito = pv.efeito, 
                        s = pv.s)
    
    p = pv.paciente
    # Tempo final (min) = linha de base + duração da janela.
    tfim = t_inducao + duracao
    # Lista de doses para o NeoPKPD (src/propofol.jl). A duração do bolus vem em
    # segundos e é passada em minutos (÷ 60). A manutenção começa
    # inicio_manutencao min após o bolus e dura até o fim da simulação;
    # `max(..., 0.0)` evita duração negativa.
    doses = esquema_inducao(p; bolus_mgkg = dose_mgkg, duracao_bolus = duracao_bolus / 60,
                            manutencao_mcgkgmin = manutencao_ugkgmin, inicio_manutencao,
                            duracao_manutencao = max(duracao - inicio_manutencao, 0.0), t_bolus = t_inducao)
    # PK/PD: Cp, Ce e BIS ao longo do tempo (grade de 1 s, tempos em min).
    pkpd = simular_propofol(p, doses; tfim, pk)
    # Circulação batimento a batimento com o θ calibrado do paciente virtual.
    bat, onda = simular_hemodinamica(p, pkpd; efeito, θbasal = pv.θ, s, tfim, forma_de_onda)
    # Junta tudo num Resultado, com o protocolo registrado para os relatórios.
    return Resultado(pkpd, bat, onda, t_inducao, pv,
                     (dose_mgkg = dose_mgkg, duracao_bolus = duracao_bolus, manutencao_ugkgmin = manutencao_ugkgmin))
end
