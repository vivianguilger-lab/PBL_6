# =============================================================================
# 03_janela_de_dose.jl — "testar a dose no computador antes da injeção"
#
# Varre doses de 0,5 a 3,0 mg/kg, calcula a janela individual (hipnose × PAM ≥ 65)
# e compara velocidades de injeção.
#
# Uso: julia --project=. scripts/03_janela_de_dose.jl dados/paciente_idoso.toml
# =============================================================================
using PropofolWindkessel
include(joinpath(@__DIR__, "graficos.jl"))

arquivo = isempty(ARGS) ? joinpath(@__DIR__, "..", "dados", "paciente_exemplo.toml") : ARGS[1]
pac, prot = ler_paciente(arquivo)
nome = splitext(basename(arquivo))[1]
limiar = get(prot, "limiar_pam", 65.0)
alvo_bis = get(prot, "alvo_bis", 60.0)
duracao_bolus = get(prot, "duracao_bolus_s", 30.0)
manut = get(prot, "manutencao_ug_kg_min", 100.0)
dprot = get(prot, "dose_mg_kg", 2.0)

pv = calibrar(pac; modelo_pk = Symbol(get(prot, "modelo_pk", "schnider")))
println(pac)
tab = varrer_doses(pv, 0.5:0.25:3.0; duracao_bolus, manutencao_ugkgmin = manut)
jan = janela_de_dose(pv; limiar_pam = limiar, alvo_bis, duracao_bolus, manutencao_ugkgmin = manut)

println("\n dose    mg    PAM mín.  queda   t mín.   PAM<65  PAM<55  BIS mín  risco")
for t in tab
    @printf(" %.2f  %5.0f   %6.1f    %4.0f%%   %4.1f min  %4.1f    %4.1f    %5.0f    %s\n",
            t.dose, t.dose * pac.peso, t.pam_nadir, t.queda_pct, t.t_nadir, t.tempo_abaixo[65.0],
            t.tempo_abaixo[55.0], t.bis_min, t.risco)
end
println()
if jan.existe
    @printf("Janela individual (bolus em %.0f s): %.2f a %.2f mg/kg  (%.0f a %.0f mg)\n",
            duracao_bolus, jan.dose_min, jan.dose_max, jan.dose_min * pac.peso, jan.dose_max * pac.peso)
else
    @printf("Sem janela: hipnose exige %.2f mg/kg, mas PAM < %.0f já com %.2f mg/kg.\n", jan.dose_min, limiar, jan.dose_max)
end

println("\nVelocidade de injeção ($(dprot) mg/kg):")
pV = plot(; title = "PAM prevista × velocidade de injeção ($(dprot) mg/kg)", xlabel = "minutos após o início",
          ylabel = "PAM (mmHg)", legend = :bottomright, size = (800, 450))
for (k, dur) in enumerate((10.0, 60.0, 120.0))
    r = simular_inducao(pv; dose_mgkg = dprot, duracao_bolus = dur, manutencao_ugkgmin = manut)
    mt = metricas(r)
    @printf("  bolus em %3.0f s: PAM mínima %.1f mmHg (%.1f min), BIS mín %.0f\n", dur, mt.pam_nadir, mt.t_nadir, mt.bis_min)
    plot!(pV, minutos(r, r.bat.t_min), r.bat.PAM; color = (COR.s1, COR.s2, COR.s3)[k], label = "injeção em $(Int(dur)) s")
end
limiares!(pV, 20.0)
savefig(fig_doses(tab, jan; titulo = "$(pac.nome) — janela de dose individual"), saida("$(nome)_doses.png"))
savefig(pV, saida("$(nome)_velocidade.png"))
println("\nFiguras em: ", normpath(DIR_RESULTADOS))
