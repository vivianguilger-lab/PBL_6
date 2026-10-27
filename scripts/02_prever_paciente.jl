# =============================================================================
# 02_prever_paciente.jl — o fluxo clínico proposto, de ponta a ponta
#
#   1. lê os biomarcadores pré-indução (arquivo TOML)
#   2. calibra o paciente virtual (Windkessel de 4 elementos: Rp e C)
#   3. prevê a curva de PAM para a dose escolhida, ANTES da injeção
#   4. relatório, figuras e série batimento a batimento (CSV)
#
# Uso: julia --project=. -t auto scripts/02_prever_paciente.jl dados/paciente_exemplo.toml
#
# Entrada: um arquivo .toml com os dados do paciente e o protocolo (exemplos em
# dados/: paciente_exemplo, paciente_jovem e paciente_idoso). Sem argumento, usa
# dados/paciente_exemplo.toml. Para outro paciente, copie um desses arquivos,
# edite os valores e passe o caminho do novo arquivo no comando.
# `-t auto` deixa o Julia usar todos os núcleos do processador.
# Saída em resultados/ (<nome> = nome do arquivo .toml, sem a extensão):
#   <nome>_inducao.png     painel com Cp/Ce, BIS, PA, FC, DC e efeitos do propofol
#   <nome>_onda.png        ondas de pressão e fluxo: batimento basal e no nadir
#   <nome>_batimentos.csv  uma linha por batimento (abre no Excel)
# No Windows também dá para usar o rodar_previsao.bat (duplo clique).
# Atenção: os parâmetros do efeito do propofol na pressão são ilustrativos (sem
# calibração com dados de pacientes). Resultado didático; não decide dose real.
# =============================================================================
# carrega o pacote do projeto e as funções de figura (graficos.jl)
using PropofolWindkessel
include(joinpath(@__DIR__, "graficos.jl"))

# Arquivo de entrada. ARGS = argumentos digitados no comando; se não houver
# nenhum, usa o paciente de exemplo. `condição ? a : b` é um "se/senão" curto.
arquivo = isempty(ARGS) ? joinpath(@__DIR__, "..", "dados", "paciente_exemplo.toml") : ARGS[1]
# pac = biomarcadores do paciente; prot = seção [protocolo] do .toml
pac, prot = ler_paciente(arquivo)
# nome do arquivo sem pasta e sem ".toml"; vira o prefixo das saídas
nome = splitext(basename(arquivo))[1]
# Protocolo: get(prot, "chave", padrão) pega o valor do .toml ou, se faltar,
# o padrão (2 mg/kg em 30 s + 100 µg/kg/min, PK de Schnider).
# Symbol(...) transforma o texto "schnider" no símbolo :schnider.
dose = get(prot, "dose_mg_kg", 2.0)
duracao_bolus = get(prot, "duracao_bolus_s", 30.0)
manut = get(prot, "manutencao_ug_kg_min", 100.0)
modelo = Symbol(get(prot, "modelo_pk", "schnider"))

# Passo 1: mostra os dados lidos ("\n" = quebra de linha)
println("\n1) Biomarcadores pré-indução\n   ", pac)
# Passo 2: paciente virtual. calibrar() calcula Rp = PAM/DC e ajusta a
# complacência C para que a pressão de pulso do modelo seja PAS − PAD.
# Rc e L ficam nos valores de Stergiopulos (1999).
println("\n2) Paciente virtual")
pv = calibrar(pac; modelo_pk = modelo)
println(pv)

# Passo 3: mostra o protocolo. $(...) dentro do texto insere o valor calculado.
println("\n3) Previsão para $(dose) mg/kg ($(round(dose * pac.peso; digits = 0)) mg) em $(duracao_bolus) s",
        manut > 0 ? " + manutenção $(manut) µg/kg/min" : "")
# Simula a indução: 1 min de linha de base, o bolus e 20 min depois dele.
# forma_de_onda = true guarda também a curva contínua de pressão e fluxo.
# `duracao_bolus` sozinho equivale a duracao_bolus = duracao_bolus.
r = simular_inducao(pv; dose_mgkg = dose, duracao_bolus, manutencao_ugkgmin = manut, forma_de_onda = true)
# Métricas: PAM basal, PAM mínima (nadir) e quando ocorre, queda (%), tempo
# abaixo de 65 e 55 mmHg, BIS mínimo e risco (ALTO < 55, MODERADO < 65, BAIXO).
mt = metricas(r)
# imprime o relatório na tela (ou no log, quando rodado por rodar_tudo.jl)
relatorio(mt; titulo = "Previsão hemodinâmica — $(pac.nome)")

# Passo 4: figuras. O título resume dose, PAM basal → mínima e risco.
# @sprintf monta um texto com formato (%.1f = 1 casa decimal, %s = texto).
titulo = @sprintf("%s | %.1f mg/kg em %.0f s | PAM %.0f → %.0f mmHg (risco %s)",
                  pac.nome, dose, duracao_bolus, mt.pam_basal, mt.pam_nadir, mt.risco)
# grava as duas figuras em resultados/ (funções definidas em graficos.jl)
savefig(fig_inducao(r; titulo), saida("$(nome)_inducao.png"))
savefig(fig_onda(r, mt), saida("$(nome)_onda.png"))

# Grava a série batimento a batimento em CSV (valores separados por vírgula).
# Colunas: t_min (minutos desde o início do bolus; negativo = linha de base),
# Ce (µg/ml), BIS, PAS, PAD e PAM (mmHg), FC (bpm), VS (ml), DC (L/min) e
# Rp (mmHg·s/ml). `open(...) do io ... end` abre o arquivo e o fecha no fim.
b = r.bat
open(saida("$(nome)_batimentos.csv"), "w") do io
    println(io, "t_min,Ce,BIS,PAS,PAD,PAM,FC,VS,DC,Rp")
    # uma linha por batimento; eachindex percorre todas as posições do vetor
    for j in eachindex(b.t_min)
        @printf(io, "%.4f,%.3f,%.1f,%.2f,%.2f,%.2f,%.2f,%.2f,%.3f,%.4f\n", b.t_min[j] - r.t_inducao, b.Ce[j], b.BIS[j],
                b.PAS[j], b.PAD[j], b.PAM[j], b.FC[j], b.VS[j], b.DC[j], b.Rp[j])
    end
end
# mostra a pasta onde ficaram as saídas
println("\nFiguras e CSV em: ", normpath(DIR_RESULTADOS))
