# =============================================================================
# paciente.jl — dados pré-indução do paciente e índices antropométricos
# =============================================================================

"""
    Paciente(; nome, idade, peso, altura, sexo, pas, pad, fc, pam=nothing, dc=nothing)

Dados que o anestesista tem **antes** da injeção (monitor + prontuário):

| campo    | unidade | significado                                              |
|:---------|:--------|:---------------------------------------------------------|
| `idade`  | anos    |                                                          |
| `peso`   | kg      |                                                          |
| `altura` | cm      |                                                          |
| `sexo`   | `:M`/`:F` |                                                        |
| `pas`    | mmHg    | pressão arterial sistólica pré-indução                   |
| `pad`    | mmHg    | pressão arterial diastólica pré-indução                  |
| `fc`     | bpm     | frequência cardíaca pré-indução                          |
| `pam`    | mmHg    | PAM medida (opcional; se ausente, PAD + (PAS−PAD)/3)     |
| `dc`     | L/min   | débito cardíaco medido (opcional; se ausente, estimado)  |
"""
Base.@kwdef struct Paciente
    nome::String = "Paciente"
    idade::Float64
    peso::Float64
    altura::Float64
    sexo::Symbol
    pas::Float64
    pad::Float64
    fc::Float64
    pam::Union{Nothing,Float64} = nothing
    dc::Union{Nothing,Float64} = nothing
    function Paciente(nome, idade, peso, altura, sexo, pas, pad, fc, pam, dc)
        sexo in (:M, :F) || throw(ArgumentError("sexo deve ser :M ou :F"))
        pas > pad || throw(ArgumentError("PAS deve ser maior que PAD"))
        new(nome, idade, peso, altura, sexo, pas, pad, fc, pam, dc)
    end
end

"""Massa magra (James, 1976) [kg] — a covariável usada no modelo de Schnider."""
function massa_magra(p::Paciente)
    r = p.peso / p.altura                    # kg/cm
    lbm = p.sexo == :M ? 1.1 * p.peso - 128 * r^2 : 1.07 * p.peso - 148 * r^2
    return lbm
end

"""Índice de massa corporal [kg/m²]."""
imc(p::Paciente) = p.peso / (p.altura / 100)^2

"""Superfície corporal (Mosteller) [m²]."""
superficie_corporal(p::Paciente) = sqrt(p.altura * p.peso / 3600)

"""Volemia total estimada (Nadler, 1962) [ml]."""
function volemia(p::Paciente)
    h = p.altura / 100
    litros = p.sexo == :M ? 0.3669 * h^3 + 0.03219 * p.peso + 0.6041 :
                            0.3561 * h^3 + 0.03308 * p.peso + 0.1833
    return 1000 * litros
end

"""PAM-alvo da calibração [mmHg]: a medida, ou PAD + (PAS − PAD)/3."""
pam_alvo(p::Paciente) = something(p.pam, p.pad + (p.pas - p.pad) / 3)

"""
    indice_cardiaco_estimado(idade)

Índice cardíaco de repouso [L/min/m²] quando o DC não é medido.
Suposição do modelo: ≈3,4 aos 20 anos com queda de ≈0,01 L/min/m² por ano,
limitado a [2,5; 3,6].
"""
indice_cardiaco_estimado(idade) = clamp(3.4 - 0.01 * (idade - 20), 2.5, 3.6)

"""Débito cardíaco-alvo da calibração [L/min]."""
dc_alvo(p::Paciente) = something(p.dc, indice_cardiaco_estimado(p.idade) * superficie_corporal(p))

function Base.show(io::IO, p::Paciente)
    @printf(io, "%s: %d anos, %s, %.0f kg, %.0f cm | PA %.0f/%.0f (PAM %.0f) mmHg | FC %.0f bpm | DC %.2f L/min%s",
            p.nome, round(Int, p.idade), p.sexo == :M ? "masc." : "fem.", p.peso, p.altura,
            p.pas, p.pad, pam_alvo(p), p.fc, dc_alvo(p), p.dc === nothing ? " (estimado)" : "")
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
                 sexo = Symbol(uppercase(pd["sexo"])),
                 pas = pd["pas"], pad = pd["pad"], fc = pd["fc"],
                 pam = get(pd, "pam", nothing), dc = get(pd, "dc", nothing))
    return p, get(d, "protocolo", Dict{String,Any}())
end
