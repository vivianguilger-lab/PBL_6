# =============================================================================
# paciente.jl — biomarcadores do paciente (estrutura do PropofolPBL do grupo)
# =============================================================================

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
    nome::String = "Paciente"
    idade::Float64
    peso::Float64
    altura::Float64
    sexo::Symbol = :F
    asa::Int = 1
    pas::Float64
    pad::Float64
    fc::Float64
    pam::Float64
    dc::Union{Nothing,Float64} = nothing
    spo2::Float64 = 98.0
    glicemia::Float64 = 95.0
    pvi::Float64 = 10.0          # pleth variability index
    ivci::Float64 = 20.0         # índice de colapsabilidade da veia cava inferior
    function Paciente(nome, idade, peso, altura, sexo, asa, pas, pad, fc, pam, dc, spo2, glicemia, pvi, ivci)
        sexo in (:M, :F) || throw(ArgumentError("sexo deve ser :M ou :F"))
        pas > pad || throw(ArgumentError("PAS deve ser maior que PAD"))
        pad < pam < pas || throw(ArgumentError("a PAM deve ficar entre PAD e PAS"))
        1 <= asa <= 5 || throw(ArgumentError("ASA deve estar entre 1 e 5"))
        new(nome, idade, peso, altura, sexo, asa, pas, pad, fc, pam, dc, spo2, glicemia, pvi, ivci)
    end
end

imc(p::Paciente) = p.peso / (p.altura / 100)^2

"Massa magra (fórmula de James), usada por Schnider. Perde validade para IMC muito alto."
function massa_magra(p::Paciente)
    r = p.peso / p.altura
    return p.sexo == :M ? 1.1 * p.peso - 128 * r^2 : 1.07 * p.peso - 148 * r^2
end

"Superfície corporal [m²] (Mosteller)."
superficie_corporal(p::Paciente) = sqrt(p.altura * p.peso / 3600)

"PAM basal [mmHg], medida no monitor."
pam_basal(p::Paciente) = p.pam

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

"""
    ler_paciente(arquivo_toml) -> (paciente, protocolo::Dict)

Lê um arquivo TOML com as seções `[paciente]` e `[protocolo]`
(veja `dados/paciente_exemplo.toml`).
"""
function ler_paciente(arquivo::AbstractString)
    d = TOML.parsefile(arquivo)
    pd = d["paciente"]
    p = Paciente(nome = get(pd, "nome", basename(arquivo)),
                 idade = pd["idade"], peso = pd["peso"], altura = pd["altura"],
                 sexo = Symbol(uppercase(pd["sexo"])), asa = get(pd, "asa", 1),
                 pas = pd["pas"], pad = pd["pad"], fc = pd["fc"],
                 pam = pd["pam"], dc = get(pd, "dc", nothing),
                 spo2 = get(pd, "spo2", 98.0), glicemia = get(pd, "glicemia", 95.0),
                 pvi = get(pd, "pvi", 10.0), ivci = get(pd, "ivci", 20.0))
    return p, get(d, "protocolo", Dict{String,Any}())
end
