# Roda todas as análises e regenera a pasta resultados/.
# Uso: julia --project=. scripts/rodar_tudo.jl
const RAIZ = normpath(joinpath(@__DIR__, ".."))
julia = Base.julia_cmd()
# cada script grava a própria saída em resultados/log_<script>[_<arquivo>].txt
function roda(script, args...; threads = "1")
    sufixo = isempty(args) ? "" : "_" * splitext(basename(string(args[1])))[1]
    log = joinpath(RAIZ, "resultados", "log_" * splitext(script)[1] * sufixo * ".txt")
    println("→ ", script, " ", join(args, " "), "   (log: ", relpath(log, RAIZ), ")")
    cmd = `$julia --project=$RAIZ -t $threads $(joinpath(RAIZ, "scripts", script)) $args`
    open(io -> run(pipeline(cmd; stdout = io, stderr = io)), log, "w")
end

roda("01_verificacao.jl")
for p in ("paciente_exemplo", "paciente_jovem", "paciente_idoso")
    roda("02_prever_paciente.jl", joinpath(RAIZ, "dados", "$p.toml"); threads = "auto")
end
for p in ("paciente_exemplo", "paciente_idoso")
    roda("03_janela_de_dose.jl", joinpath(RAIZ, "dados", "$p.toml"))
end
# 04_coorte_virtual.jl é experimental (não validado) e não entra aqui
roda("05_sensibilidade.jl", joinpath(RAIZ, "dados", "paciente_exemplo.toml"); threads = "auto")
roda("06_identificacao_windkessel.jl")
