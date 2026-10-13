"""
    PropofolWindkessel

Modelo in silico para predição hemodinâmica da indução anestésica com
propofol: farmacocinética de Schnider + farmacodinâmica de Hill acopladas a
um modelo de Windkessel em malha fechada com coração de elastância variável
e barorreflexo, calibrado com os dados pré-indução de cada paciente.

Pergunta norteadora: *como prever, por simulação, a resposta hemodinâmica
de um paciente ao propofol?*
"""
module PropofolWindkessel

using OrdinaryDiffEqTsit5
using LinearAlgebra
using Printf
using Random
using Statistics
using TOML

include("paciente.jl")
include("farmacocinetica.jl")
include("farmacodinamica.jl")
include("windkessel.jl")
include("simulacao.jl")
include("guyton.jl")
include("calibracao.jl")
include("metricas.jl")
include("dose.jl")
include("populacao.jl")
include("incerteza.jl")

# paciente
export Paciente, ler_paciente, massa_magra, imc, superficie_corporal, volemia, pam_alvo, dc_alvo
# PK / PD
export ParametrosPK, schnider, matriz_pk, meia_vida_ke0_hemodinamica, TrechoInfusao, protocolo_inducao,
       massa_total, ParametrosPD, hill, efeitos, inot_Emax
# circulação
export ParametrosCV, ParametrosBaro, ObservacaoPA, amplificacao_pp, ModeloInducao, com, ativacao, algebrico, rhs!, estado_inicial,
       volume_total, NESTADOS
export I_A1, I_A2, I_A3, I_CEB, I_CEH, I_FASE, I_VLV, I_VA, I_VV, I_PS, I_XFC, I_XR, I_XE, I_XV, I_INTP, I_INTQ
# simulação e análise
export Resultado, integrar, batimentos, simular_inducao, media_movel
export PacienteVirtual, calibrar, tau_diastolico
export metricas, relatorio, varrer_doses, janela_de_dose, bissecao
export previsao_com_incerteza, amostrar_pd, amostrar_eta_pk
export Individuo, gerar_coorte, calibrar_coorte, simular_coorte, resumo_coorte, dose_padrao, OMEGA_PK, OMEGA_PD
export pressao_media_enchimento, resistencia_retorno_venoso, retorno_venoso, debito_cardiaco_curso, estado_guyton

end # module
