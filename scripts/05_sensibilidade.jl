# =============================================================================
# 05_sensibilidade.jl — análise de sensibilidade local (±10%)
#
# Quais parâmetros do modelo do grupo mais mudam a PAM mínima prevista?
# Os mais influentes são os que mais precisam ser calibrados com dados.
#
# Uso: julia --project=. -t auto scripts/05_sensibilidade.jl dados/paciente_exemplo.toml
#
# Método: muda um parâmetro por vez em −10% e +10%, mantendo os outros, e mede
# quanto muda a PAM mínima (nadir) nos 20 min após o bolus.
# S normalizada = variação relativa da PAM mínima ÷ variação relativa do
# parâmetro (20%). |S| perto de 1: a PAM mínima muda na mesma proporção.
# Entrada: arquivo .toml (padrão: dados/paciente_exemplo.toml).
# Saída: tabela na tela e resultados/<nome>_sensibilidade.png (barras).
# É sensibilidade LOCAL: vale perto dos valores atuais, que são ilustrativos.
# =============================================================================
# carrega o pacote do projeto e as funções de figura (graficos.jl)
using PropofolWindkessel
include(joinpath(@__DIR__, "graficos.jl"))

# lê o paciente e o protocolo (como em 02_prever_paciente.jl)
arquivo = isempty(ARGS) ? joinpath(@__DIR__, "..", "dados", "paciente_exemplo.toml") : ARGS[1]
pac, prot = ler_paciente(arquivo)
nome = splitext(basename(arquivo))[1]
dose = get(prot, "dose_mg_kg", 2.0)
dur = get(prot, "duracao_bolus_s", 30.0)
manut = get(prot, "manutencao_ug_kg_min", 100.0)
# Referência: efeito hemodinâmico padrão do grupo (Ce50, γ, Imax), paciente
# virtual calibrado e a sensibilidade s dele.
ef0 = EfeitoHemodinamico()
pv0 = calibrar(pac)
s0 = pv0.s

# nadir(...): simula a indução e devolve a PAM mínima. Cada opção tem como
# padrão o valor de referência; basta passar a que se quer mudar.
#   ef = efeito hemodinâmico; eta = fatores da PK (ex.: (V1 = 1.1,));
#   s = sensibilidade; dose (mg/kg), dur (s), manut (µg/kg/min); pv = paciente virtual.
# O `;` no início indica que todos os argumentos são passados pelo nome.
function nadir(; ef = ef0, eta = NamedTuple(), s = s0, dose = dose, dur = dur, manut = manut, pv = pv0)
    # parâmetros PK de Schnider multiplicados pelos fatores eta
    pk = parametros_pk(pv.paciente, :schnider; eta)
    # simula e devolve só o campo pam_nadir das métricas
    return metricas(simular_inducao(pv; pk, efeito = ef, s, dose_mgkg = dose, duracao_bolus = dur,
                                    manutencao_ugkgmin = manut)).pam_nadir
end
# comef(...): cópia do efeito de referência trocando só o que for passado.
# Ex.: comef(Ce50 = 3.3) muda o Ce50 e mantém γ e os Imax.
# `...` "espalha" a lista de campos como opções da função.
comef(; kw...) = EfeitoHemodinamico(; (k => getfield(ef0, k) for k in fieldnames(EfeitoHemodinamico))..., kw...)
# comC(f): multiplica a complacência C calibrada pelo fator f.
# ajustar_C = false evita recalibrar C (o ajuste desfaria a mudança).
function comC(f)
    pv = calibrar(pac; ajustar_C = false)
    pv.θ[2] = pv0.θ[2] * f
    return nadir(; pv)
end

# Parâmetros testados: (rótulo, função do fator f). f = 0.9 é −10% e f = 1.1
# é +10%. `f -> ...` é uma função anônima: recebe f e devolve a PAM mínima.
parametros = [
    ("Ce50 hemodinâmico",            f -> nadir(ef = comef(Ce50 = ef0.Ce50 * f))),
    ("γ (Hill) hemodinâmico",        f -> nadir(ef = comef(γ = ef0.γ * f))),
    ("Imax Rp (vasodilatação)",      f -> nadir(ef = comef(Imax_Rp = ef0.Imax_Rp * f))),
    ("Imax VS (volume sistólico)",   f -> nadir(ef = comef(Imax_VS = ef0.Imax_VS * f))),
    ("Imax FC",                      f -> nadir(ef = comef(Imax_FC = ef0.Imax_FC * f))),
    ("sensibilidade s",              f -> nadir(s = s0 * f)),
    ("ke0 (Schnider)",               f -> nadir(eta = (ke0 = f,))),
    ("V1 (Schnider)",                f -> nadir(eta = (V1 = f,))),
    ("CL (Schnider)",                f -> nadir(eta = (CL = f,))),
    ("C (complacência arterial)",    comC),
    ("dose (mg/kg)",                 f -> nadir(dose = dose * f)),
    ("manutenção (µg/kg/min)",       f -> nadir(manut = manut * f)),
    ("tempo de injeção",             f -> nadir(dur = dur * f)),
]

# PAM mínima de referência (todos os parâmetros nos valores originais)
base = nadir()
@printf("%s | %.1f mg/kg em %.0f s + %.0f µg/kg/min | PAM mínima de referência = %.1f mmHg\n\n", pac.nome, dose, dur, manut, base)
# Calcula −10% e +10% de cada parâmetro em paralelo: Threads.@threads divide
# as voltas do laço entre os núcleos. `undef` cria o vetor ainda vazio.
linhas = Vector{Any}(undef, length(parametros))
Threads.@threads for k in eachindex(parametros)
    rot, f = parametros[k]
    # m e p = PAM mínima com o parâmetro em 90% e em 110%
    m, p = f(0.9), f(1.1)
    # diferenças em relação à referência (mmHg) e S normalizada
    linhas[k] = (rot = rot, menos = m - base, mais = p - base, S = ((p - m) / base) / 0.2)
end
# ordena pelo maior efeito em mmHg; `sort!` altera o próprio vetor
sort!(linhas; by = l -> max(abs(l.menos), abs(l.mais)))
# imprime do mais influente para o menos influente (reverse)
println(" parâmetro                              −10%      +10%     S normalizada")
for l in reverse(linhas)
    @printf(" %-36s %+6.2f   %+6.2f    %+5.2f\n", l.rot, l.menos, l.mais, l.S)
end

# Figura: barras horizontais com a variação da PAM mínima para −10% e +10%.
# xm = limite simétrico do eixo x; os parâmetros mais influentes ficam no topo.
y = 1:length(linhas)
xm = 1.15 * max(maximum(max(abs(l.menos), abs(l.mais)) for l in linhas), 0.1)
p = bar(y, [l.menos for l in linhas]; orientation = :h, bar_width = 0.7, color = COR.s1, lw = 0, label = "−10%",
        yticks = (y, [l.rot for l in linhas]), xlabel = "variação da PAM mínima (mmHg)",
        title = @sprintf("Sensibilidade local — %s (PAM mínima %.1f mmHg)", pac.nome, base),
        size = (950, 560), legend = :bottomright, left_margin = 6Plots.mm,
        xlims = (-xm, xm), ylims = (0.4, length(linhas) + 0.6))
bar!(p, y, [l.mais for l in linhas]; orientation = :h, bar_width = 0.7, color = COR.s2, lw = 0, label = "+10%")
vline!(p, [0]; color = COR.eixo, lw = 1, label = "")
# grava a figura em resultados/
savefig(p, saida("$(nome)_sensibilidade.png"))
println("\nFigura em: ", normpath(saida("$(nome)_sensibilidade.png")))
