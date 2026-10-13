# =============================================================================
# 05_sensibilidade.jl — análise de sensibilidade local (±10%)
#
# Credibilidade (Resumo, Seção 8.3): quais parâmetros mais mudam a previsão?
# Os mais influentes são os que mais precisam de dados para serem
# individualizados — e indicam onde o modelo é frágil.
#
# Saída: variação do nadir de PAM (mmHg) para −10% e +10% de cada parâmetro
# e a sensibilidade normalizada S = (ΔPAM/PAM)/(Δp/p).
#
# Uso: julia --project=. -t auto scripts/05_sensibilidade.jl dados/paciente_exemplo.toml
# =============================================================================
using PropofolWindkessel
include(joinpath(@__DIR__, "graficos.jl"))

arquivo = isempty(ARGS) ? joinpath(@__DIR__, "..", "dados", "paciente_exemplo.toml") : ARGS[1]
pac, prot = ler_paciente(arquivo)
nome = splitext(basename(arquivo))[1]
dose = get(prot, "dose_mg_kg", 2.0)
dur = get(prot, "duracao_bolus_s", 20.0)
pd0 = ParametrosPD()
baro0 = ParametrosBaro()

function nadir(; pd = pd0, eta = NamedTuple(), baro = baro0, dose = dose, dur = dur, p = pac)
    pv = calibrar(p; pd, eta_pk = eta, baro_ref = baro)
    return metricas(simular_inducao(pv; dose_mgkg = dose, duracao_bolus = dur, saveat = 0.01)).pam_nadir
end
pacm(; kw...) = Paciente(; nome = pac.nome, idade = pac.idade, peso = pac.peso, altura = pac.altura, sexo = pac.sexo,
                         pas = pac.pas, pad = pac.pad, fc = pac.fc, pam = pac.pam, dc = dc_alvo(pac), kw...)

# (rótulo, função f(fator) → nadir)
parametros = [
    ("C50 TPR (Su 2023)",            f -> nadir(pd = com(pd0; tpr_C50 = pd0.tpr_C50 * f))),
    ("Emax TPR (Su 2023)",           f -> nadir(pd = com(pd0; tpr_Emax = pd0.tpr_Emax * f))),
    ("Emax inotrópico (Su 2023)",    f -> nadir(pd = com(pd0; inot_Emax35 = pd0.inot_Emax35 * f))),
    ("Emax venodilatação (assumido)", f -> nadir(pd = com(pd0; veno_Emax = pd0.veno_Emax * f))),
    ("C50 venodilatação (assumido)", f -> nadir(pd = com(pd0; veno_C50 = pd0.veno_C50 * f))),
    ("atenuação barorreflexo (assumido)", f -> nadir(pd = com(pd0; baro_Emax = min(pd0.baro_Emax * f, 0.99)))),
    ("ke0 vascular (Kazama 1999)",   f -> nadir(eta = (ke0_hemo = f,))),
    ("V1 (Schnider)",                f -> nadir(eta = (V1 = f,))),
    ("Cl1 (Schnider)",               f -> nadir(eta = (Cl1 = f,))),
    ("ganho barorreflexo TPR",       f -> nadir(baro = com(baro0; A_R = baro0.A_R * f))),
    ("ganho barorreflexo FC",        f -> nadir(baro = com(baro0; A_FC = baro0.A_FC * f))),
    ("DC pré-indução (estimado)",    f -> nadir(p = pacm(dc = dc_alvo(pac) * f))),
    ("FC pré-indução",               f -> nadir(p = pacm(fc = pac.fc * f))),
    ("dose (mg/kg)",                 f -> nadir(dose = dose * f)),
    ("tempo de injeção",             f -> nadir(dur = dur * f)),
]

base = nadir()
@printf("%s | %.1f mg/kg em %.0f s | PAM nadir de referência = %.1f mmHg\n\n", pac.nome, dose, dur, base)
linhas = Vector{Any}(undef, length(parametros))
Threads.@threads for k in eachindex(parametros)
    rot, f = parametros[k]
    m, p = f(0.9), f(1.1)
    linhas[k] = (rot = rot, menos = m - base, mais = p - base, S = ((p - m) / base) / 0.2)
end
sort!(linhas; by = l -> max(abs(l.menos), abs(l.mais)))
println(" parâmetro                              −10%      +10%     S normalizada")
for l in reverse(linhas)
    @printf(" %-36s %+6.2f   %+6.2f    %+5.2f\n", l.rot, l.menos, l.mais, l.S)
end

y = 1:length(linhas)
xm = 1.15 * maximum(max(abs(l.menos), abs(l.mais)) for l in linhas)
p = bar(y, [l.menos for l in linhas]; orientation = :h, bar_width = 0.7, color = COR.s1, lw = 0, label = "−10%",
        yticks = (y, [l.rot for l in linhas]), xlabel = "variação do nadir de PAM (mmHg)",
        title = @sprintf("Sensibilidade local — %s (nadir %.1f mmHg)", pac.nome, base),
        size = (950, 560), legend = :bottomright, left_margin = 6Plots.mm,
        xlims = (-xm, xm), ylims = (0.4, length(linhas) + 0.6))
bar!(p, y, [l.mais for l in linhas]; orientation = :h, bar_width = 0.7, color = COR.s2, lw = 0, label = "+10%")
vline!(p, [0]; color = COR.eixo, lw = 1, label = "")
savefig(p, saida("$(nome)_sensibilidade.png"))
println("\nFigura em: ", normpath(saida("$(nome)_sensibilidade.png")))
