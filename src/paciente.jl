# =============================================================================
# paciente.jl — biomarcadores do paciente (estrutura do PropofolPBL do grupo)
# =============================================================================
#
# O que este arquivo faz:
# - Define o "Paciente": a ficha com os dados que o anestesista tem ANTES da
#   injeção (idade, peso, altura, sexo, ASA, PAS, PAD, PAM, FC, DC opcional,
#   PVI, IVCI). É a entrada de todo o modelo in silico.
# - Calcula grandezas derivadas usadas pelos outros arquivos: IMC, massa magra
#   (James, usada pela PK de Schnider), superfície corporal (Mosteller), PAM e
#   débito cardíaco (DC) basais, que montam o Windkessel basal.
# - Lê os dados de um arquivo de texto TOML (ex.: dados/paciente_exemplo.toml).
#
# Onde se encaixa no modelo:
#   Paciente -> propofol.jl (PK/PD e BIS) -> acoplamento.jl (efeito do propofol
#   sobre o Windkessel de 4 elementos) -> pressão arterial prevista.
# A estrutura (campos e valores padrão) segue o repositório do grupo,
# PropofolPBL (src/paciente.jl). A origem de cada número está em
# docs/validacao.md (seções 4.1 e 4.5).

# -----------------------------------------------------------------------------
# Estrutura Paciente
# Em Julia, `struct` define um "tipo" de dado com campos nomeados (como uma
# ficha com campos fixos). `Base.@kwdef` (logo abaixo do texto de ajuda)
# permite criar o paciente escrevendo os campos pelo nome, em qualquer ordem:
#   Paciente(idade = 62, peso = 72.0, altura = 160.0, sexo = :F,
#            pas = 148.0, pad = 82.0, fc = 78.0, pam = 104.0)
# Campos com `= valor` têm valor padrão e podem ser omitidos; os outros são
# obrigatórios. O texto entre aspas triplas logo abaixo é a "docstring": a
# ajuda que aparece ao digitar `?Paciente` no Julia.
# -----------------------------------------------------------------------------
"""
    Paciente(; nome, idade, peso, altura, sexo, asa, pas, pad, fc, pam, dc,
               spo2, glicemia, pvi, ivci)

Dados que o anestesista tem antes da injeção. Unidades: idade em anos, peso em
kg, altura em cm, pressões em mmHg, FC em bpm, DC em L/min, PVI e IVCI em %.

Uso no modelo:
- PK (Schnider/Marsh/Eleveld): idade, peso, altura, sexo.
- BIS (Eleveld): idade, peso, altura, sexo.
- Hemodinâmica basal (Windkessel): PAM e FC do monitor, superfície corporal (DC
  opcional); PAS e PAD calibram a complacência arterial C (paciente virtual).
  A PAM é obrigatória: os monitores oscilométricos a medem diretamente, como na
  entrada do PropofolPBL do grupo.
- Sensibilidade hemodinâmica ao propofol: idade, ASA, PVI, IVCI.
- `spo2` e `glicemia` não entram no modelo; existem só para manter a mesma
  estrutura do PropofolPBL do grupo.
"""
Base.@kwdef struct Paciente
    # nome: texto livre, só para identificar o paciente nos relatórios. OK
    nome::String = "Paciente"
    # idade [anos], peso [kg] e altura [cm]: covariáveis da PK (Schnider, Marsh,
    # Eleveld) e do BIS de Eleveld (propofol.jl). OK
    idade::Float64 #OK
    peso::Float64 #OK
    altura::Float64
    sexo::Symbol = :F
    # se você fornece o sexo do paciente, tudo certo. Mas se deixar o default, o programa silenciosamente assume que o paciente é feminino VERIFICAR!
    asa::Int = 1
    # pas, pad [mmHg]: pressões sistólica e diastólica pré-indução; calibram a
    # complacência arterial C do Windkessel (paciente virtual, calibracao.jl).
    # fc [bpm]: frequência cardíaca basal.
    # Paciente sem ASA informado -> ASA 1 -> modelo interpreta paciente como saudável VERIFICAR!
    pas::Float64 # OK 
    pad::Float64 # OK
    fc::Float64 # OK
    pam::Float64 # dá pra calcular a partir dos outros dados também (PAM = PAD + (PAS-PAD/3)) VERIFICAR!
    dc::Union{Nothing,Float64} = nothing
    # o modelo permite que seja o DC medido ou (se não ouver valor medido) 3.0 * superficie_corporal(p), ou seja, DC = 3,0 \times SC VERIFICAR!
    spo2::Float64 = 98.0
    glicemia::Float64 = 95.0
    pvi::Float64 = 10.0          # pleth variability index   VALOR ILUSTRATIVO -> ACHAR ALGO PARA MUDAR!!
    ivci::Float64 = 20.0         # índice de colapsabilidade da veia cava inferior   VALOR ILUSTRATIVO -> ACHAR ALGO PARA MUDAR!!
   
    function Paciente(nome, idade, peso, altura, sexo, asa, pas, pad, fc, pam, dc, spo2, glicemia, pvi, ivci)
        sexo in (:M, :F) || throw(ArgumentError("sexo deve ser :M ou :F"))
        pas > pad || throw(ArgumentError("PAS deve ser maior que PAD"))
        # A PAM medida tem de ficar entre a diastólica e a sistólica.
        pad < pam < pas || throw(ArgumentError("a PAM deve ficar entre PAD e PAS"))
        1 <= asa <= 5 || throw(ArgumentError("ASA deve estar entre 1 e 5"))
        new(nome, idade, peso, altura, sexo, asa, pas, pad, fc, pam, dc, spo2, glicemia, pvi, ivci)
    end
end

# Índice de massa corporal (IMC) [kg/m²] = peso / altura², com a altura
# convertida de cm para m. Usado na massa livre de gordura de Al-Sallami
# (modelo de Eleveld, propofol.jl).
# Sintaxe: `f(x) = expressão` define uma função curta em uma linha; `p.peso` lê
# o campo `peso` do paciente `p`; `p::Paciente` exige que `p` seja um Paciente;
imc(p::Paciente) = p.peso / (p.altura / 100)^2 # OK

# Massa magra (lean body mass, LBM) [kg] pela fórmula de James (1976):
#   homem:  1,1·peso − 128·(peso/altura)²
#   mulher: 1,07·peso − 148·(peso/altura)²      (peso em kg, altura em cm)
# É covariável do clearance (CL) no modelo de Schnider (1998). Confere com o
# PropofolPBL do grupo e com o Python Anesthesia Simulator (docs/validacao.md,
# seção 4.1). Limitação: em obesos o termo negativo cresce mais rápido que o
# positivo e a massa magra calculada passa a CAIR com o ganho de peso (ex.:
# mulher de 160 cm acima de ~92 kg), o que não tem sentido fisiológico.
# Sintaxe: `condição ? a : b` vale `a` se a condição for verdadeira, senão `b`.
"Massa magra (fórmula de James), usada por Schnider. Perde validade para IMC muito alto." # fórmula publicada (James 1976, usada por Schnider 1998) 
function massa_magra(p::Paciente)
    # r = peso / altura [kg/cm]; entra ao quadrado na fórmula.
    r = p.peso / p.altura
    return p.sexo == :M ? 1.1 * p.peso - 128 * r^2 : 1.07 * p.peso - 148 * r^2
end

# Superfície corporal (SC) [m²] pela fórmula de Mosteller (1987):
#   SC = raiz quadrada de (altura [cm] × peso [kg] / 3600).
# Serve para estimar o débito cardíaco basal quando ele não foi medido
# (função debito_basal, abaixo). `sqrt` = raiz quadrada.
# Ex.: 160 cm e 72 kg -> SC ≈ 1,79 m².
"Superfície corporal [m²] (Mosteller)."
superficie_corporal(p::Paciente) = sqrt(p.altura * p.peso / 3600) # dado derivado: altura real -> peso real -> Mosteller -> SC derivada     OK

# PAM basal [mmHg]: é o próprio valor do monitor (não é calculada como
# PAD + (PAS − PAD)/3), porque os monitores oscilométricos medem a PAM
# diretamente. É a PAM de partida do Windkessel: a resistência periférica
# basal sai de Rp = PAM / DC (docs/modelo_matematico.md, seção 2).
"PAM basal [mmHg], medida no monitor."
pam_basal(p::Paciente) = p.pam # OK

# Débito cardíaco basal [L/min]. `something(a, b)` devolve `a` se ele existir
# (não for `nothing`) e `b` caso contrário. Ou seja: usa o DC medido, se houver;
# senão estima DC = índice cardíaco × SC, com índice cardíaco de 3,0 L/min/m².
# O 3,0 é escolha do grupo (PropofolPBL), dentro da faixa normal 2,5–4,0
# L/min/m². Ele não muda a queda relativa da PAM; só muda a forma da onda
# (PAS e PAD), via complacência C calibrada (docs/validacao.md, seção 4.5).
# Ex.: SC 1,79 m² -> DC estimado ≈ 5,4 L/min.
"""
Débito cardíaco basal [L/min]: o medido, ou 3,0 L/min/m² × superfície corporal
(como no grupo; índice cardíaco normal 2,5–4,0 L/min/m²). A queda prevista da PAM
não depende desse valor, porque o modelo multiplica Rp, VS e FC por fatores relativos.
"""
debito_basal(p::Paciente) = something(p.dc, 3.0 * superficie_corporal(p))

function Base.show(io::IO, p::Paciente)
    @printf(io, "%s: %d anos, %s, %.0f kg, %.0f cm, ASA %d | PA %.0f/%.0f (PAM %.0f) mmHg | FC %.0f bpm | DC %.2f L/min%s | PVI %.0f%%, IVCI %.0f%%",
            p.nome, round(Int, p.idade), p.sexo == :M ? "masc." : "fem.", p.peso, p.altura, p.asa,
            p.pas, p.pad, pam_basal(p), p.fc, debito_basal(p), p.dc === nothing ? " (estimado)" : "", p.pvi, p.ivci)
end

# -----------------------------------------------------------------------------
# Leitura do paciente a partir de um arquivo TOML: arquivo de texto simples com
# linhas `chave = valor`, organizadas em seções entre colchetes ([paciente],
# [protocolo]). Entrada: caminho do arquivo. Saída: uma dupla (tupla) com o
# Paciente e um dicionário com o protocolo (dose, modelo PK, limiar de PAM...).
# -----------------------------------------------------------------------------
"""
    ler_paciente(arquivo_toml) -> (paciente, protocolo::Dict)

Lê um arquivo TOML com as seções `[paciente]` e `[protocolo]`
(veja `dados/paciente_exemplo.toml`).
"""
function ler_paciente(arquivo::AbstractString)
    # `d` é um dicionário (Dict): uma tabela "nome -> valor" com todo o arquivo.
    # `d["paciente"]` pega a seção [paciente].
    d = TOML.parsefile(arquivo)
    pd = d["paciente"]
    # Monta o Paciente. `pd["idade"]` exige a chave (dá erro se faltar): idade,
    # peso, altura, sexo, PAS, PAD, FC e PAM são obrigatórios.
    # `get(pd, "asa", 1)` lê a chave se ela existir; senão usa o valor padrão
    # (aqui 1). O sexo vem como texto ("f" ou "F") e vira o símbolo :F
    # (`uppercase` põe em maiúscula, `Symbol` converte). DC ausente fica
    # `nothing` (será estimado). O nome padrão é o nome do arquivo (`basename`).
    p = Paciente(nome = get(pd, "nome", basename(arquivo)),
                 idade = pd["idade"], peso = pd["peso"], altura = pd["altura"],
                 sexo = Symbol(uppercase(pd["sexo"])), asa = get(pd, "asa", 1),
                 pas = pd["pas"], pad = pd["pad"], fc = pd["fc"],
                 pam = pd["pam"], dc = get(pd, "dc", nothing),
                 spo2 = get(pd, "spo2", 98.0), glicemia = get(pd, "glicemia", 95.0),
                 pvi = get(pd, "pvi", 10.0), ivci = get(pd, "ivci", 20.0))
    # Devolve o paciente e a seção [protocolo]; se ela não existir, devolve um
    # dicionário vazio (`Dict{String,Any}()`).
    return p, get(d, "protocolo", Dict{String,Any}())
end
