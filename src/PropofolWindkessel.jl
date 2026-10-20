"""
    PropofolWindkessel

Modelo in silico para predição hemodinâmica da indução anestésica com propofol,
seguindo o PropofolPBL do grupo
(github.com/gabiadamek30-bit/Simula-o-da-resposta-do-paciente-ao-propofol---PBL):

- PK/PD do propofol com o NeoPKPD (Schnider/Marsh, compartimento de efeito, BIS de Eleveld);
- Windkessel de 4 elementos do windkessel-id (Stergiopulos 1999), alimentado
  pela onda de fluxo aórtico escalada por volume sistólico e FC;
- acoplamento Ce → Rp, VS e FC por Emax sigmoide com sensibilidade do paciente
  (idade, ASA, PVI, IVCI).

Sobre essa base: paciente virtual calibrado com PAS/PAD pré-indução, métricas de
hipotensão, janela de dose, faixa de incerteza, coorte virtual e plataforma web.
"""
module PropofolWindkessel

using NeoPKPD
using Optim, LineSearches, ForwardDiff
using HTTP
using JSON
using LinearAlgebra
using Printf
using Random
using Statistics
using TOML

include("windkessel.jl")
include("dados.jl")
include("paciente.jl")
include("propofol.jl")
include("acoplamento.jl")
include("calibracao.jl")
include("metricas.jl")
include("dose.jl")
include("populacao.jl")
include("incerteza.jl")
include("plataforma.jl")

# Windkessel (windkessel-id)
export WKSistema, G, G2, G3, lsim_wk, estado_final, simulate_wk, solvex0, bode_wk
export IdConditions, identify, identifyglobal, deal, cost, params_stergiopulos1999total
export WINDKESSEL_ID_DIR, importPQ, dados_stergiopulos1999, mLs2Lmin
# Paciente e propofol (NeoPKPD)
export Paciente, ler_paciente, imc, massa_magra, superficie_corporal, debito_basal, pam_basal
export parametros_pk, bolus, infusao, esquema_inducao, bis_eleveld, simular_propofol, matriz_pk
# Acoplamento
export EfeitoHemodinamico, fator_efeito, sensibilidade_hemodinamica, theta_basal, forma_fluxo,
       fluxo_batimento, batimento_basal, simular_hemodinamica, simular_paciente
# Paciente virtual, previsão e análises
export PacienteVirtual, calibrar, tau_diastolico, Resultado, simular_inducao
export metricas, relatorio, varrer_doses, janela_de_dose, bissecao
export previsao_com_incerteza, amostrar_eta_pk, amostrar_efeito, OMEGA_PK, OMEGA_HEMO
export Individuo, gerar_coorte, calibrar_coorte, simular_coorte, resumo_coorte, dose_padrao
export plataforma

end # module
