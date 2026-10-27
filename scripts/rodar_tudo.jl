# Roda todas as análises e regenera a pasta resultados/.
# Uso: julia --project=. scripts/rodar_tudo.jl
#
# O que faz: chama, um por um, os scripts 01, 02, 03, 05 e 06, cada um num
# processo Julia separado. Na tela aparece só o nome de cada script; tudo o
# que o script imprime vai para um arquivo de log em resultados/.
# Rode a partir da pasta do projeto (PBL_6). Leva alguns minutos; na primeira
# vez o Julia ainda baixa e compila os pacotes.
# Saída em resultados/: figuras (.png), tabelas (.csv) e logs (log_*.txt).
#
# RAIZ = pasta do projeto. @__DIR__ é a pasta deste arquivo (scripts/) e ".."
# sobe um nível. `const` marca um valor que não muda depois.
const RAIZ = normpath(joinpath(@__DIR__, ".."))
# comando que inicia o mesmo Julia que está rodando agora (mesma versão)
julia = Base.julia_cmd()
# roda(script, args...; threads): roda um script num processo Julia novo.
#   script  = nome do arquivo em scripts/ (ex.: "01_verificacao.jl");
#   args... = argumentos extras, em geral o .toml do paciente ("..." aceita
#             quantos argumentos vierem, inclusive nenhum);
#   threads = quantos núcleos do processador usar ("auto" = todos).
# cada script grava a própria saída em resultados/log_<script>[_<arquivo>].txt
function roda(script, args...; threads = "1")
    # sufixo do log: nome do .toml sem a extensão (ex.: "_paciente_idoso").
    # `condição ? a : b` é um "se/senão" em uma linha.
    sufixo = isempty(args) ? "" : "_" * splitext(basename(string(args[1])))[1]
    log = joinpath(RAIZ, "resultados", "log_" * splitext(script)[1] * sufixo * ".txt")
    # mostra na tela qual script está rodando e onde fica o log
    println("→ ", script, " ", join(args, " "), "   (log: ", relpath(log, RAIZ), ")")
    # monta o comando (entre crases) sem executá-lo ainda; $(...) insere valores.
    # --project=RAIZ faz o Julia usar os pacotes de Project.toml/Manifest.toml.
    cmd = `$julia --project=$RAIZ -t $threads $(joinpath(RAIZ, "scripts", script)) $args`
    # abre o log para escrita ("w") e executa o comando, mandando para ele a
    # saída normal (stdout) e as mensagens de erro (stderr).
    # `io -> ...` é uma função anônima que recebe o arquivo aberto.
    open(io -> run(pipeline(cmd; stdout = io, stderr = io)), log, "w")
end

# 01: verificação numérica (Windkessel e NeoPKPD contra soluções exatas)
roda("01_verificacao.jl")
# 02: previsão completa para os três pacientes de exemplo de dados/.
# "$p.toml" troca $p pelo nome do paciente da volta atual do laço.
for p in ("paciente_exemplo", "paciente_jovem", "paciente_idoso")
    roda("02_prever_paciente.jl", joinpath(RAIZ, "dados", "$p.toml"); threads = "auto")
end
# 03: janela de dose individual para os pacientes A (exemplo) e C (idoso)
for p in ("paciente_exemplo", "paciente_idoso")
    roda("03_janela_de_dose.jl", joinpath(RAIZ, "dados", "$p.toml"))
end
# 04_coorte_virtual.jl é experimental (não validado) e não entra aqui
# Para rodar o 04 à parte, veja o comando no cabeçalho de 04_coorte_virtual.jl.
# 05: sensibilidade local (±10%) para o paciente A
roda("05_sensibilidade.jl", joinpath(RAIZ, "dados", "paciente_exemplo.toml"); threads = "auto")
# 06: identificação dos parâmetros do Windkessel a partir de dados medidos
roda("06_identificacao_windkessel.jl")
