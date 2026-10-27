# =============================================================================
# PropofolWindkessel.jl — arquivo principal do pacote (o "módulo")
# =============================================================================
# Em Julia, um módulo agrupa todo o código do pacote sob um nome. Quem digita
# `using PropofolWindkessel` carrega este arquivo, que por sua vez carrega os
# demais arquivos de src/ (include) e diz quais funções ficam disponíveis para o
# usuário (export). Não há cálculo aqui, só organização.
# O texto entre aspas triplas logo abaixo é a docstring do módulo: a ajuda que
# aparece ao digitar ?PropofolWindkessel no Julia.
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

# Bibliotecas externas usadas (a lista com versões está em Project.toml):
#   NeoPKPD — motor de farmacocinética/farmacodinâmica (modelo tricompartimental,
#     compartimento de efeito), como no PropofolPBL do grupo;
#   Optim, LineSearches, ForwardDiff — otimização e derivadas automáticas para
#     identificar os parâmetros do Windkessel (windkessel.jl);
#   HTTP, JSON — servidor web local da plataforma para médicos (plataforma.jl);
#   LinearAlgebra, Printf, Random, Statistics, TOML — álgebra linear, texto
#     formatado, sorteios, médias e leitura dos arquivos de paciente (.toml).
using NeoPKPD
using Optim, LineSearches, ForwardDiff
using HTTP
using JSON
using LinearAlgebra
using Printf
using Random
using Statistics
using TOML

# Ordem dos include: cada arquivo pode usar tipos e constantes definidos nos
# anteriores (em Julia, um tipo precisa existir antes de aparecer na definição de
# outro). A ordem segue o fluxo do modelo, do bloco mais básico ao mais completo.
#
# windkessel.jl: Windkessel de 2, 3 e 4 elementos (porta do windkessel-id feita
# pelo grupo): circuito, simulação, estado periódico e identificação de parâmetros.
include("windkessel.jl")
# dados.jl: lê o batimento de Stergiopulos 1999 (a onda de fluxo aórtico que
# alimenta o Windkessel) e os dados in vivo de porcos.
include("dados.jl")
# paciente.jl: estrutura Paciente com os biomarcadores pré-indução (idade, peso,
# PAS, PAD, PAM, FC, ASA, PVI, IVCI...) e grandezas derivadas (massa magra,
# superfície corporal, DC basal).
include("paciente.jl")
# propofol.jl: PK do propofol (Schnider, Marsh ou Eleveld) com o NeoPKPD,
# concentração no sítio efetor (ke0) e BIS de Eleveld 2018.
include("propofol.jl")
# acoplamento.jl: liga a concentração no efetor (Ce) a Rp, VS e FC (Emax sigmoide
# × sensibilidade do paciente; valores ilustrativos do grupo) e calcula a pressão
# batimento a batimento no Windkessel. Usa os quatro arquivos anteriores.
include("acoplamento.jl")
# calibracao.jl: paciente virtual (Rp = PAM/DC; C ajustada para reproduzir
# PAS − PAD) e simulação da indução.
include("calibracao.jl")
# metricas.jl: PAM mínima, tempo e área abaixo de 65 e 55 mmHg, BIS e risco.
include("metricas.jl")
# dose.jl: varredura de doses e janela de dose (usa calibracao.jl e metricas.jl).
include("dose.jl")
# populacao.jl: coorte virtual. EXPERIMENTAL, não validado (fora da plataforma;
# docs/validacao.md, seção 6.3).
include("populacao.jl")
# incerteza.jl: faixa de incerteza da previsão. EXPERIMENTAL, não validado; usa as
# variabilidades de populacao.jl, por isso vem depois dele.
include("incerteza.jl")
# plataforma.jl: servidor web local para médicos. Usa todos os anteriores, por
# isso vem por último.
include("plataforma.jl")

# export: lista os nomes que ficam disponíveis direto depois de
# `using PropofolWindkessel` (sem precisar escrever PropofolWindkessel.nome).
# Funções internas, com "_" no início do nome (como _wk), não são exportadas.
#
# Windkessel (windkessel-id)
# (modelos G2/G3/G, simulação, estado periódico e resposta em frequência)
export WKSistema, G, G2, G3, lsim_wk, estado_final, simulate_wk, solvex0, bode_wk
# (identificação de parâmetros e valores de referência de Stergiopulos 1999)
export IdConditions, identify, identifyglobal, deal, cost, params_stergiopulos1999total
# (leitura dos dados, em dados.jl)
export WINDKESSEL_ID_DIR, importPQ, dados_stergiopulos1999, mLs2Lmin
# Paciente e propofol (NeoPKPD)
export Paciente, ler_paciente, imc, massa_magra, superficie_corporal, debito_basal, pam_basal
export parametros_pk, parametros_eleveld, atraso_bis, bolus, infusao, esquema_inducao, bis_eleveld,
       simular_propofol, matriz_pk
# Acoplamento
export EfeitoHemodinamico, fator_efeito, sensibilidade_hemodinamica, theta_basal, forma_fluxo,
       fluxo_batimento, batimento_basal, simular_hemodinamica, simular_paciente
# Paciente virtual, previsão e análises
export PacienteVirtual, calibrar, tau_diastolico, Resultado, simular_inducao
export metricas, relatorio, varrer_doses, janela_de_dose, bissecao
# (os dois grupos abaixo são EXPERIMENTAIS, não validados: incerteza.jl e populacao.jl)
export previsao_com_incerteza, amostrar_eta_pk, amostrar_efeito, OMEGA_PK, OMEGA_HEMO
export Individuo, gerar_coorte, calibrar_coorte, simular_coorte, resumo_coorte, dose_padrao
# Plataforma web: plataforma() abre a interface em http://localhost:8080.
export plataforma

# Fim do módulo: tudo entre "module PropofolWindkessel" e "end" pertence ao pacote.
end # module
