# =============================================================================
# 05_sensibilidade.jl — análise de sensibilidade local (±10%)
#
# Quais parâmetros do modelo do grupo mais mudam a PAM mínima prevista?
# Os mais influentes são os que mais precisam ser calibrados com dados.
#
# Uso: julia --project=. -t auto scripts/05_sensibilidade.jl dados/paciente_exemplo.toml
# =============================================================================
using PropofolWindkessel
include(joinpath(@__DIR__, "graficos.jl"))

arquivo = isempty(ARGS) ? joinpath(@__DIR__, "..", "dados", "paciente_exemplo.toml") : ARGS[1]
pac, prot = ler_paciente(arquivo)
nome = splitext(basename(arquivo))[1]
dose = get(prot, "dose_mg_kg", 2.0)
dur = get(prot, "duracao_bolus_s", 30.0)
manut = get(prot, "manutencao_ug_kg_min", 100.0)
ef0 = EfeitoHemodinamico()
pv0 = calibrar(pac)
s0 = pv0.s

function nadir(; ef = ef0, eta = NamedTuple(), s = s0, dose = dose, dur = dur, manut = manut, pv = pv0)
    pk = parametros_pk(pv.paciente, :schnider; eta)
    return metricas(simular_inducao(pv; pk, efeito = ef, s, dose_mgkg = dose, duracao_bolus = dur,
                                    manutencao_ugkgmin = manut)).pam_nadir
end
comef(; kw...) = EfeitoHemodinamico(; (k => getfield(ef0, k) for k in fieldnames(EfeitoHemodinamico))..., kw...)
function comC(f)
    pv = calibrar(pac; ajustar_C = false)
    pv.θ[2] = pv0.θ[2] * f
    return nadir(; pv)
end

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

base = nadir()
@printf("%s | %.1f mg/kg em %.0f s + %.0f µg/kg/min | PAM mínima de referência = %.1f mmHg\n\n", pac.nome, dose, dur, manut, base)
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
xm = 1.15 * max(maximum(max(abs(l.menos), abs(l.mais)) for l in linhas), 0.1)
p = bar(y, [l.menos for l in linhas]; orientation = :h, bar_width = 0.7, color = COR.s1, lw = 0, label = "−10%",
        yticks = (y, [l.rot for l in linhas]), xlabel = "variação da PAM mínima (mmHg)",
        title = @sprintf("Sensibilidade local — %s (PAM mínima %.1f mmHg)", pac.nome, base),
        size = (950, 560), legend = :bottomright, left_margin = 6Plots.mm,
        xlims = (-xm, xm), ylims = (0.4, length(linhas) + 0.6))
bar!(p, y, [l.mais for l in linhas]; orientation = :h, bar_width = 0.7, color = COR.s2, lw = 0, label = "+10%")
vline!(p, [0]; color = COR.eixo, lw = 1, label = "")
savefig(p, saida("$(nome)_sensibilidade.png"))
println("\nFigura em: ", normpath(saida("$(nome)_sensibilidade.png")))
