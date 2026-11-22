
# ============================================================
# remifentanil.jl
# Modelo farmacocinetico populacional de Minto (1997)
#
# Etapa 1: parametros individuais.
# Nao implementa efeitos hemodinamicos.
# ============================================================

"""
    parametros_minto(idade, lbm)

Calcula os parametros farmacocineticos populacionais do
remifentanil segundo Minto et al. (1997).

Unidades:
- idade: anos
- lbm: kg de massa corporal magra
- V1, V2, V3: L
- CL, Q2, Q3: L/min
- ke0: min^-1
"""
function parametros_minto(idade::Real, lbm::Real)

    idade = Float64(idade)
    lbm = Float64(lbm)

    all(isfinite, (idade, lbm)) ||
        throw(ArgumentError("Idade e massa magra devem ser finitas"))

    idade > 0 || throw(ArgumentError("Idade invalida"))
    lbm > 0 || throw(ArgumentError("Massa magra invalida"))

    delta_idade = idade - 40.0
    delta_lbm = lbm - 55.0

    V1 = 5.1 - 0.0201 * delta_idade +
         0.072 * delta_lbm

    V2 = 9.82 - 0.0811 * delta_idade +
         0.108 * delta_lbm

    V3 = 5.42

    CL = 2.6 - 0.0162 * delta_idade +
         0.0191 * delta_lbm

    Q2 = 2.05 - 0.0301 * delta_idade

    Q3 = 0.076 - 0.00113 * delta_idade

    ke0 = 0.595 - 0.007 * delta_idade

    parametros = (V1, V2, V3, CL, Q2, Q3, ke0)

    all(x -> isfinite(x) && x > 0, parametros) ||
        throw(ArgumentError(
            "Parametros nao fisicos; verificar aplicabilidade de Minto"
        ))

    return (
        V1=V1,
        V2=V2,
        V3=V3,
        CL=CL,
        Q2=Q2,
        Q3=Q3,
        ke0=ke0
    )
end

# Reutiliza o tipo Paciente e a massa magra do projeto.
function parametros_minto(p::Paciente)
    return parametros_minto(p.idade, massa_magra(p))
end


# ============================================================
# Simulacao PK do remifentanil
# Tres compartimentos + sitio de efeito
# Integracao numerica RK4
# ============================================================

"""
    simular_remifentanil(tempos_s, taxas_ml_h, pk;
                        concentracao_ug_ml=20.0,
                        permitir_zero_inicial=false)

Simula a PK do remifentanil usando o historico de infusao.

Hipoteses desta primeira versao:
- registros em intervalos de 1 segundo, iniciando em t=0;
- taxa constante entre duas amostras consecutivas;
- quantidades e Ce inicialmente iguais a zero;
- valores ausentes internos nao sao permitidos;
- preenchimento inicial com zero somente quando autorizado.

Saidas de Cp e Ce: ng/mL.
"""
function simular_remifentanil(
    tempos_s::AbstractVector,
    taxas_ml_h::AbstractVector,
    pk;
    concentracao_ug_ml::Float64=20.0,
    permitir_zero_inicial::Bool=false
)

    n = length(tempos_s)

    n >= 2 || throw(ArgumentError(
        "Historico deve conter pelo menos dois instantes"
    ))

    length(taxas_ml_h) == n ||
        throw(ArgumentError(
            "Tempos e taxas possuem tamanhos diferentes"
        ))

    t = Float64.(tempos_s)

    all(isfinite, t) ||
        throw(ArgumentError("Tempos invalidos"))

    isapprox(t[1], 0.0; atol=1e-9) ||
        throw(ArgumentError(
            "Historico deve comecar em t_case_s = 0"
        ))

    all(isapprox.(diff(t), 1.0; atol=1e-9)) ||
        throw(ArgumentError(
            "Historico deve ter amostragem de 1 segundo"
        ))

    concentracao_ug_ml > 0 ||
        throw(ArgumentError("Concentracao invalida"))

    taxa = [
        ismissing(v) ? NaN : Float64(v)
        for v in taxas_ml_h
    ]

    any(v -> isfinite(v) && v < 0, taxa) &&
        throw(ArgumentError("Taxa negativa encontrada"))

    # Admite somente valores ausentes no comeco.
    faltantes_iniciais = 0

    while faltantes_iniciais < n &&
          isnan(taxa[faltantes_iniciais + 1])
        faltantes_iniciais += 1
    end

    if faltantes_iniciais > 0
        permitir_zero_inicial ||
            throw(ArgumentError(
                "Taxa inicial desconhecida"
            ))

        faltantes_iniciais <= 5 ||
            throw(ArgumentError(
                "Lacuna inicial superior a 5 segundos"
            ))

        taxa[1:faltantes_iniciais] .= 0.0
    end

    all(isfinite, taxa) ||
        throw(ArgumentError(
            "Historico possui lacunas internas ou valores invalidos"
        ))

    # Estado: A1, A2, A3 [microgramas], Ce [ng/mL].
    # A1, A2, A3, Ce, dose administrada e dose eliminada
    estado = zeros(Float64, 6)

    administrado_ug = zeros(Float64, n)
    eliminado_ug = zeros(Float64, n)
    remanescente_ug = zeros(Float64, n)

    cp = zeros(Float64, n)
    ce = zeros(Float64, n)

    
    function derivada(x, entrada_ug_min)
        A1, A2, A3, Ce, _, _ = x

        C1 = A1 / pk.V1
        C2 = A2 / pk.V2
        C3 = A3 / pk.V3

        eliminacao = pk.CL * C1

        dA1 = entrada_ug_min -
            eliminacao -
            pk.Q2 * (C1 - C2) -
            pk.Q3 * (C1 - C3)

        dA2 = pk.Q2 * (C1 - C2)
        dA3 = pk.Q3 * (C1 - C3)

        dCe = pk.ke0 * (C1 - Ce)

        # Acumuladores de massa, em microgramas.
        dAdministrado = entrada_ug_min
        dEliminado = eliminacao

        return [
            dA1, dA2, dA3, dCe,
            dAdministrado, dEliminado
        ]
    end


    for i in 1:(n - 1)

        # Conversao mL/h -> microgramas/min.
        entrada = taxa[i] * concentracao_ug_ml / 60.0

        # Passo de integracao: 1 segundo = 1/60 minuto.
        h = (t[i + 1] - t[i]) / 60.0

        k1 = derivada(estado, entrada)
        k2 = derivada(estado .+ (h / 2) .* k1, entrada)
        k3 = derivada(estado .+ (h / 2) .* k2, entrada)
        k4 = derivada(estado .+ h .* k3, entrada)

        estado .= estado .+
                  (h / 6) .* (k1 .+ 2 .* k2 .+
                              2 .* k3 .+ k4)

        all(isfinite, estado) ||
            error("Instabilidade numerica em t=$(t[i + 1])")

        minimum(estado) >= -1e-8 ||
            error("Estado PK negativo em t=$(t[i + 1])")

        cp[i + 1] = estado[1] / pk.V1
        ce[i + 1] = estado[4]
        administrado_ug[i + 1] = estado[5]
        eliminado_ug[i + 1] = estado[6]
        remanescente_ug[i + 1] = sum(estado[1:3])
    end

    return (
        t_case_s=t,
        taxa_ml_h=taxa,
        cp_ng_ml=cp,
        ce_ng_ml=ce,
        administrado_ug=administrado_ug,
        eliminado_ug=eliminado_ug,
        remanescente_ug=remanescente_ug
    )
end
