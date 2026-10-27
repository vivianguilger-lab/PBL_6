# =============================================================================
# plataforma.jl — interface web local para uso clínico
# =============================================================================
#
#   using PropofolWindkessel
#   plataforma()          # abre o navegador em http://localhost:8080
#
# O navegador mostra um formulário (biomarcadores pré-indução + dose) e os
# gráficos; todo o cálculo é feito por este pacote em Julia (NeoPKPD +
# Windkessel do windkessel-id). O servidor só escuta em 127.0.0.1: os dados do
# paciente não saem do computador.
#
# Como funciona (para quem não programa): o Julia vira um pequeno servidor web
# local. O navegador abre app/index.html, o médico preenche o formulário e a
# página envia os dados em JSON (texto estruturado) para dois endereços:
#   /api/simular -> calibra o paciente virtual, simula a dose escolhida e
#                   devolve métricas e curvas (PAM, PAS, PAD, FC, DC, BIS, Cp, Ce);
#   /api/doses   -> simula doses de 0,5 a 3,0 mg/kg e devolve a tabela e a
#                   janela de dose.
# Só a previsão do paciente típico aparece: incerteza e coorte (experimentais,
# não validadas) ficam de fora. Os parâmetros do efeito do propofol na pressão
# são ilustrativos (docs/validacao.md, seção 4.4).
# =============================================================================

# Constantes globais (`const`; o `_` no início do nome indica uso interno):
#   _HTML     = caminho da página (app/index.html); `@__DIR__` = pasta deste arquivo;
#   _CACHE_PV = dicionário que guarda pacientes virtuais já calibrados (calibrar
#               leva tempo; com os mesmos dados, reaproveita);
#   _TRAVA    = trava (lock) que impede dois pedidos simultâneos de mexer no
#               cache ao mesmo tempo.
const _HTML = joinpath(@__DIR__, "..", "app", "index.html")
const _CACHE_PV = Dict{Any,PacienteVirtual}()
const _TRAVA = ReentrantLock()

# Funções auxiliares de uma linha:
#   _num: converte o valor do formulário em número; aceita vírgula decimal
#         ("2,5" -> 2.5) e espaços nas pontas.
#   _opcional: campo opcional. `get(d, k, nothing)` lê a chave k do dicionário
#         ou devolve `nothing` (vazio) se ela não existir; campo em branco também
#         vira nothing.
#   _padrao: `something(a, x)` devolve `a` se não for nothing; senão, o padrão `x`.
#   _r / _rv: arredonda um número / cada número de um vetor para n casas; NaN e
#         infinito viram nothing (null no JSON, que não aceita NaN).
#   _faixa: confere a <= x <= b; senão interrompe com a mensagem `msg`
#         (ArgumentError), que aparece para o usuário na tela.
_num(x) = x isa Number ? Float64(x) : parse(Float64, replace(strip(string(x)), "," => "."))
_opcional(d, k) = (v = get(d, k, nothing); (v === nothing || v == "") ? nothing : _num(v))
_padrao(d, k, x) = something(_opcional(d, k), x)
_r(x, n = 2) = isfinite(x) ? round(x; digits = n) : nothing
_rv(v, n = 2) = [_r(x, n) for x in v]
_faixa(x, a, b, msg) = a <= x <= b || throw(ArgumentError(msg))

# -----------------------------------------------------------------------------
# _ler_paciente: lê e valida os biomarcadores do formulário (`d["paciente"]`).
# Obrigatórios: idade (anos), peso (kg), altura (cm), PAS, PAD, PAM (mmHg) e
# FC (bpm). Opcionais com padrão: sexo (M), ASA (1), PVI (10%), IVCI (20%), DC
# (L/min, estimado se vazio). A faixa de idade, peso e altura é a aceita pelo
# modelo (adultos de 18–95 anos, docs/validacao.md, seções 2 e 3); as de PA, FC,
# PVI e IVCI só barram erros de digitação. Devolve um `Paciente` (paciente.jl).
# -----------------------------------------------------------------------------
function _ler_paciente(d)
    pac = d["paciente"]
    # Sexo em maiúsculas; qualquer valor diferente de M ou F é recusado.
    sexo = uppercase(string(get(pac, "sexo", "M")))
    sexo in ("M", "F") || throw(ArgumentError("Sexo deve ser M ou F."))
    # Campos obrigatórios: `pac["idade"]` não tem padrão; se faltar, o Julia gera
    # um KeyError, que `_rotear` transforma em "Campo obrigatório ausente".
    idade, peso, altura = _num(pac["idade"]), _num(pac["peso"]), _num(pac["altura"])
    pas, pad, pam, fc = _num(pac["pas"]), _num(pac["pad"]), _num(pac["pam"]), _num(pac["fc"])
    # ASA vira inteiro (`round(Int, ...)`). PVI e IVCI têm padrão abaixo dos
    # limiares do fator de sensibilidade (15% e 43%): sem efeito se não informados.
    asa = round(Int, _padrao(pac, "asa", 1.0))
    pvi, ivci = _padrao(pac, "pvi", 10.0), _padrao(pac, "ivci", 20.0)
    # Validações. A PAM precisa ficar entre a PAD e a PAS; a pressão de pulso
    # (PAS − PAD) >= 15 mmHg garante um pulso plausível, usado para calibrar a
    # complacência C do paciente virtual.
    _faixa(idade, 18, 95, "Idade fora da faixa do modelo (18–95 anos).")
    _faixa(peso, 35, 200, "Peso fora da faixa do modelo (35–200 kg).")
    _faixa(altura, 130, 215, "Altura fora da faixa do modelo (130–215 cm).")
    (60 <= pas <= 250 && 30 <= pad <= 150) || throw(ArgumentError("Pressão arterial fora da faixa plausível."))
    pas - pad >= 15 || throw(ArgumentError("A PAS deve ser pelo menos 15 mmHg maior que a PAD."))
    pad < pam < pas || throw(ArgumentError("A PAM deve ficar entre a PAD e a PAS."))
    _faixa(fc, 35, 160, "Frequência cardíaca fora da faixa plausível (35–160 bpm).")
    _faixa(asa, 1, 5, "ASA deve estar entre 1 e 5.")
    _faixa(pvi, 0, 60, "PVI fora da faixa plausível (0–60%).")
    _faixa(ivci, 0, 100, "IVCI fora da faixa plausível (0–100%).")
    # Monta o Paciente. `Symbol(sexo)` transforma "M" em :M (forma usada no modelo).
    # DC vazio: o modelo estima 3,0 L/min/m² × superfície corporal (paciente.jl).
    return Paciente(nome = string(get(pac, "nome", "Paciente")), idade = idade, peso = peso, altura = altura,
                    sexo = Symbol(sexo), asa = asa, pas = pas, pad = pad, fc = fc,
                    pam = pam, dc = _opcional(pac, "dc"), pvi = pvi, ivci = ivci)
end

# -----------------------------------------------------------------------------
# _ler_protocolo: lê o esquema de administração (`d["protocolo"]`, opcional).
#   dose_mgkg: bolus em mg/kg (padrão 2,0; aceito de 0 a 4);
#   modelo_pk: Schnider (padrão), Marsh ou Eleveld; é a PK que move a
#              hemodinâmica (o BIS usa sempre Eleveld);
#   duracao_bolus_s: duração da injeção em segundos (padrão 30; mínimo 1);
#   manutencao_ugkgmin: infusão de manutenção em µg/kg/min (padrão 100; >= 0).
# Padrões coerentes com a bula do Diprivan (docs/validacao.md, seção 4.7).
# Devolve dois valores: o protocolo (NamedTuple) e o modelo (Symbol, ex. :schnider).
# -----------------------------------------------------------------------------
function _ler_protocolo(d)
    pr = get(d, "protocolo", Dict())
    dose = _num(get(pr, "dose_mgkg", 2.0))
    _faixa(dose, 0, 4, "Dose fora da faixa do simulador (0–4 mg/kg).")
    modelo = Symbol(lowercase(string(get(pr, "modelo_pk", "schnider"))))
    modelo in (:schnider, :marsh, :eleveld) || throw(ArgumentError("Modelo PK deve ser Schnider, Marsh ou Eleveld."))
    return (dose_mgkg = dose, duracao_bolus = max(_num(get(pr, "duracao_bolus_s", 30.0)), 1.0),
            manutencao_ugkgmin = max(_num(get(pr, "manutencao_ugkgmin", 100.0)), 0.0)), modelo
end

# -----------------------------------------------------------------------------
# _virtual: paciente virtual (calibracao.jl) com cache. A `chave` reúne todos os
# dados que mudam a calibração; se a mesma chave já foi calculada, reaproveita.
# `lock(_TRAVA) do ... end` executa o bloco com a trava fechada: só um pedido
# por vez lê ou escreve no cache. `get!(f, dict, chave)` devolve o valor
# guardado ou, se não houver, calcula com `f()` e guarda. Acima de 200
# pacientes o cache é esvaziado, para não ocupar memória sem limite.
# -----------------------------------------------------------------------------
"""Paciente virtual calibrado (em cache: só recalibra se os dados mudarem)."""
function _virtual(p::Paciente, modelo::Symbol)
    chave = (p.idade, p.peso, p.altura, p.sexo, p.asa, p.pas, p.pad, p.fc, p.pam, p.dc, p.pvi, p.ivci, modelo)
    lock(_TRAVA) do
        length(_CACHE_PV) > 200 && empty!(_CACHE_PV)
        get!(() -> calibrar(p; modelo_pk = modelo), _CACHE_PV, chave)
    end
end

# -----------------------------------------------------------------------------
# _resumo_virtual: dados do paciente virtual que a página mostra, num Dict
# (vira objeto JSON). Inclui PA, DC e VS reproduzidos pelo modelo; θ = [Rp, C,
# Rc, L] (mmHg·s/ml, ml/mmHg, mmHg·s/ml, mmHg·s²/ml); τ = Rp·C (s); a
# sensibilidade s; massa magra, superfície corporal e IMC; a PK escolhida (CL
# em L/min, volumes em L, ke0 em /min) e o BIS de Eleveld: Ce50 =
# 3,08·e^(−0,00635·(idade − 35)) µg/ml, ke0 (/min) e atraso (s).
# Fontes de cada valor em docs/validacao.md, seções 4.1 a 4.3.
# -----------------------------------------------------------------------------
function _resumo_virtual(pv::PacienteVirtual)
    # Desempacota θ em quatro variáveis e cria atalhos para os campos usados.
    Rp, C, Rc, L = pv.θ
    o, p, pk = pv.obtido, pv.paciente, pv.pk
    # `"chave" => valor` forma cada par do Dict; `_r(x, n)` arredonda para n casas.
    # "dc_estimado" é verdadeiro quando o DC não foi digitado (foi estimado).
    return Dict("convergiu" => pv.convergiu, "pam" => _r(o.PAM, 1), "pas" => _r(o.PAS, 1), "pad" => _r(o.PAD, 1),
                "dc" => _r(o.DC), "dc_estimado" => p.dc === nothing, "vs" => _r(o.VS, 0),
                "rp" => _r(Rp, 3), "c" => _r(C, 2), "rc" => _r(Rc, 3), "l" => _r(L, 4), "tau" => _r(Rp * C, 2),
                "s" => _r(pv.s, 2), "lbm" => _r(massa_magra(p), 1), "bsa" => _r(superficie_corporal(p), 2),
                "imc" => _r(imc(p), 1), "modelo_pk" => uppercasefirst(string(pk.modelo)),
                "cl" => _r(pk.CL, 2), "v1" => _r(pk.V1, 2), "v2" => _r(pk.V2, 1), "v3" => _r(pk.V3, 0),
                "ke0" => _r(pk.ke0, 3), "ce50_bis" => _r(3.08 * exp(-0.00635 * (p.idade - 35)), 2),
                "ke0_bis" => _r(parametros_eleveld(p).ke0, 3), "atraso_bis" => _r(atraso_bis(p.idade), 0))
end

# -----------------------------------------------------------------------------
# _api_simular: responde ao endereço /api/simular. Passos: ler paciente e
# protocolo -> paciente virtual (cache) -> `simular_inducao` (20 min após o
# bolus) -> `metricas`. Devolve um Dict com o resumo do paciente virtual, a
# dose total em mg, as métricas e as séries para os gráficos. Tempos em min
# contados a partir do bolus (negativos = linha de base).
# -----------------------------------------------------------------------------
function _api_simular(d)
    p = _ler_paciente(d)
    pr, modelo = _ler_protocolo(d)
    pv = _virtual(p, modelo)
    # `pr...` "espalha" os campos do protocolo como opções nomeadas
    # (dose_mgkg = ..., duracao_bolus = ..., manutencao_ugkgmin = ...).
    r = simular_inducao(pv; pr...)
    mt = metricas(r)
    b, k = r.bat, r.pkpd
    t0 = r.t_inducao
    # Os sinais da PK/PD estão a cada 1 s; para aliviar a página, envia um ponto
    # a cada 5 (`1:5:n` = 1, 6, 11, ...). `vs0` = VS basal, para as variações em %.
    idx = 1:5:length(k.t)                       # 5 s
    vs0 = pv.obtido.VS
    # Resposta. No nadir, variações em % do basal: Rp (vasodilatação), VS
    # (depressão miocárdica / pré-carga) e FC. Em "batimentos", o ponto em `.-`,
    # `.*` e `./` aplica a conta a cada batimento (broadcasting).
    # `mt.tempo_abaixo[65.0]` lê o tempo abaixo de 65 mmHg no Dict das métricas.
    return Dict(
        "virtual" => _resumo_virtual(pv),
        "dose_mg" => _r(pr.dose_mgkg * p.peso, 0),
        "metricas" => Dict(
            "pam_basal" => _r(mt.pam_basal, 1), "pam_nadir" => _r(mt.pam_nadir, 1),
            "queda_pct" => _r(mt.queda_pct, 1), "t_nadir" => _r(mt.t_nadir, 1),
            "pas_nadir" => _r(mt.pas_nadir, 0), "pad_nadir" => _r(mt.pad_nadir, 0),
            "fc_basal" => _r(mt.fc_basal, 0), "fc_nadir" => _r(mt.fc_no_nadir, 0),
            "dc_basal" => _r(mt.dc_basal), "dc_nadir" => _r(mt.dc_no_nadir),
            "min_abaixo_65" => _r(mt.tempo_abaixo[65.0], 1), "min_abaixo_55" => _r(mt.tempo_abaixo[55.0], 1),
            "aut_65" => _r(mt.aut[65.0], 0), "aut_55" => _r(mt.aut[55.0], 0),
            "bis_min" => _r(mt.bis_min, 0), "t_hipnose" => _r(mt.t_hipnose, 1),
            "tempo_hipnose" => _r(mt.tempo_hipnose, 1), "ce_max" => _r(mt.ce_max), "risco" => mt.risco,
            "rp_nadir_pct" => _r(100 * (mt.rp_no_nadir / pv.θ[1] - 1), 0),
            "vs_nadir_pct" => _r(100 * (mt.vs_no_nadir / vs0 - 1), 0),
            "fc_nadir_pct" => _r(100 * (mt.fc_no_nadir / p.fc - 1), 0)),
        "batimentos" => Dict(
            "t" => _rv(b.t_min .- t0, 3), "pam" => _rv(b.PAM, 1), "pas" => _rv(b.PAS, 1),
            "pad" => _rv(b.PAD, 1), "fc" => _rv(b.FC, 1), "dc" => _rv(b.DC, 2),
            "rp" => _rv(100 .* (b.Rp ./ pv.θ[1] .- 1), 1), "vs" => _rv(100 .* (b.VS ./ vs0 .- 1), 1),
            "fcp" => _rv(100 .* (b.FC ./ p.fc .- 1), 1)),
        "sinais" => Dict("t" => _rv(k.t[idx] .- t0, 3), "bis" => _rv(k.bis[idx], 1),
                         "cp" => _rv(k.Cp[idx], 2), "ce" => _rv(k.Ce[idx], 2)))
end

# _cruzamento: numa tabela (x = doses, y = métrica), acha o primeiro trecho em
# que y passa pelo `alvo` (y − alvo muda de sinal) e interpola linearmente o x
# do cruzamento. Devolve `nothing` se y nunca cruza o alvo.
"""Interpola o ponto em que `y` cruza `alvo` numa varredura monótona."""
function _cruzamento(x, y, alvo)
    # Percorre pares de pontos vizinhos (i, i+1); a condição exige troca de sinal
    # e evita divisão por zero quando y[i] == y[i+1].
    for i in 1:length(x)-1
        if (y[i] - alvo) * (y[i+1] - alvo) <= 0 && y[i] != y[i+1]
            return x[i] + (alvo - y[i]) * (x[i+1] - x[i]) / (y[i+1] - y[i])
        end
    end
    return nothing
end

# -----------------------------------------------------------------------------
# _api_doses: responde a /api/doses. Simula 11 doses, 0,5; 0,75; ...; 3,0 mg/kg
# (`0.5:0.25:3.0`), com o protocolo do formulário, e monta:
#   - a tabela (dose, mg totais, PAM mínima, queda %, min abaixo de 65 e 55,
#     BIS mínimo, risco);
#   - a janela de dose, por interpolação entre doses vizinhas (mais rápida que a
#     bisseção de dose.jl): dose mínima com BIS < 60 e máxima com PAM mínima >= 65.
# -----------------------------------------------------------------------------
function _api_doses(d)
    p = _ler_paciente(d)
    pr, modelo = _ler_protocolo(d)
    pv = _virtual(p, modelo)
    doses = collect(0.5:0.25:3.0)
    tab = varrer_doses(pv, doses; duracao_bolus = pr.duracao_bolus, manutencao_ugkgmin = pr.manutencao_ugkgmin)
    # Colunas da tabela: PAM mínima e BIS mínimo de cada dose.
    nadir = [t.pam_nadir for t in tab]
    bis = [t.bis_min for t in tab]
    # Dose mínima: se a menor dose (0,5) já leva o BIS < 60, fica 0,5 (menor valor
    # testado, não necessariamente o limite real; docs/validacao.md, seção 6.2);
    # senão, interpola onde o BIS mínimo cruza 60.
    # Dose máxima: se até 3,0 mg/kg mantém a PAM mínima >= 65, fica 3,0; se já 0,5
    # fica abaixo de 65, não há (nothing); senão, interpola onde ela cruza 65.
    dmin = bis[1] < 60 ? doses[1] : _cruzamento(doses, bis, 60.0)
    dmax = nadir[end] >= 65 ? doses[end] : (nadir[1] < 65 ? nothing : _cruzamento(doses, nadir, 65.0))
    # Uma linha por dose (`[Dict(...) for t in tab]`). A janela "existe" se as
    # duas doses foram achadas e a mínima não passa da máxima.
    return Dict("doses" => [Dict("dose" => t.dose, "mg" => _r(t.dose * p.peso, 0), "pam_nadir" => _r(t.pam_nadir, 1),
                                 "queda_pct" => _r(t.queda_pct, 0), "min_abaixo_65" => _r(t.tempo_abaixo[65.0], 1),
                                 "min_abaixo_55" => _r(t.tempo_abaixo[55.0], 1), "bis_min" => _r(t.bis_min, 0),
                                 "risco" => t.risco) for t in tab],
                "janela" => Dict("dose_min" => dmin === nothing ? nothing : _r(dmin, 2),
                                 "dose_max" => dmax === nothing ? nothing : _r(dmax, 2),
                                 "existe" => dmin !== nothing && dmax !== nothing && dmin <= dmax))
end

# _json: monta a resposta HTTP com o código `status` (200 = ok, 400 = erro no
# pedido) e o conteúdo `x` convertido para texto JSON.
_json(status, x) = HTTP.Response(status, ["Content-Type" => "application/json; charset=utf-8"], JSON.json(x))

# -----------------------------------------------------------------------------
# _rotear: recebe cada pedido do navegador e decide o que responder:
#   GET  / ou /index.html          -> a página (app/index.html);
#   GET  /api/saude                -> confirma que o servidor está no ar e a versão;
#   POST /api/simular ou /api/doses -> lê o JSON, calcula e devolve JSON;
#   qualquer outro endereço        -> 404 (não encontrado).
# `try ... catch` captura erros: dado inválido (ArgumentError) ou campo ausente
# (KeyError) viram mensagem clara; outros erros também são registrados no
# terminal (`@error`). Em todos os casos de erro a página recebe o código 400.
# -----------------------------------------------------------------------------
function _rotear(req::HTTP.Request)
    try
        # Separa o endereço dos parâmetros depois de "?" (ex.: /api/saude?x=1).
        alvo = first(split(req.target, '?'))
        if req.method == "GET" && alvo in ("/", "/index.html")
            return HTTP.Response(200, ["Content-Type" => "text/html; charset=utf-8"], read(_HTML, String))
        elseif req.method == "GET" && alvo == "/api/saude"
            return _json(200, Dict("ok" => true, "versao" => string(pkgversion(@__MODULE__))))
        elseif req.method == "POST" && alvo in ("/api/simular", "/api/doses")
            # Converte o corpo do pedido (texto JSON) em Dict e escolhe a função.
            d = JSON.parse(String(req.body))
            f = alvo == "/api/simular" ? _api_simular : _api_doses
            return _json(200, f(d))
        end
        return HTTP.Response(404, "não encontrado")
    catch e
        # Mensagem de erro para a página, conforme o tipo de erro.
        msg = e isa ArgumentError ? e.msg : e isa KeyError ? "Campo obrigatório ausente: $(e.key)" :
              "Erro no cálculo: " * sprint(showerror, e)
        e isa ArgumentError || e isa KeyError || @error "plataforma" exception = (e, catch_backtrace())
        return _json(400, Dict("erro" => msg))
    end
end

# _abrir_navegador: abre o endereço no navegador padrão com o comando de cada
# sistema operacional: Windows `cmd /c start`, macOS `open`, Linux `xdg-open`.
# Comandos entre crases (`...`) são programas externos; `wait = false` não espera
# o navegador fechar. Se falhar, ignora: o endereço é impresso na tela.
function _abrir_navegador(url)
    try
        if Sys.iswindows()
            run(`cmd /c start "" $url`; wait = false)
        elseif Sys.isapple()
            run(`open $url`; wait = false)
        else
            run(`xdg-open $url`; wait = false)
        end
    catch
    end
end

# -----------------------------------------------------------------------------
# plataforma: função que o usuário chama para abrir o simulador.
#   porta   = 8080 (se ocupada, tenta 8081, ..., 8089);
#   abrir   = abre o navegador sozinho;
#   aquecer = roda uma simulação de teste antes. A primeira execução de cada
#             função em Julia é lenta (compilação); assim a primeira resposta ao
#             médico já sai rápida.
# Devolve o servidor (`srv`); `close(srv)` o encerra.
# -----------------------------------------------------------------------------
"""
    plataforma(; porta=8080, abrir=true) -> servidor

Abre a plataforma web do simulador no navegador (http://localhost:porta).
O cálculo é feito por este pacote; o servidor só aceita conexões do próprio
computador. Para encerrar: `close(servidor)` ou feche o Julia.
"""
function plataforma(; porta = 8080, abrir = true, aquecer = true)
    if aquecer
        println("Preparando o simulador (na primeira vez leva cerca de 1 minuto)...")
        # Paciente fictício (50 anos, 70 kg) só para compilar; o resultado é descartado.
        corpo = Dict("paciente" => Dict("idade" => 50, "peso" => 70, "altura" => 170, "sexo" => "M",
                                        "pas" => 125, "pad" => 75, "pam" => 92, "fc" => 70),
                     "protocolo" => Dict("dose_mgkg" => 1.0))
        JSON.json(_api_simular(corpo))
        JSON.json(_api_doses(corpo))
    end
    # Liga o servidor só em 127.0.0.1 (o próprio computador): nenhum outro
    # computador da rede consegue acessar, e os dados não saem da máquina.
    # `HTTP.serve!` inicia o servidor em segundo plano, com `_rotear` tratando os
    # pedidos. Se a porta estiver em uso, dá erro e o laço tenta a seguinte; na
    # 10ª falha (`tentativa == 9`), `rethrow()` repassa o erro ao usuário.
    srv = nothing
    for tentativa in 0:9                      # se a porta estiver ocupada, tenta a seguinte
        try
            srv = HTTP.serve!(_rotear, "127.0.0.1", porta + tentativa; verbose = false)
            porta += tentativa
            break
        catch e
            tentativa == 9 && rethrow()
        end
    end
    # Endereço final (com a porta que funcionou), avisos e abertura do navegador.
    # `$porta` dentro de uma string insere o valor da variável.
    url = "http://localhost:$porta"
    println("\nPlataforma pronta: $url")
    println("Se o navegador não abrir sozinho, copie o endereço acima no navegador.")
    println("Para encerrar: feche esta janela do Julia.")
    abrir && _abrir_navegador(url)
    return srv
end
