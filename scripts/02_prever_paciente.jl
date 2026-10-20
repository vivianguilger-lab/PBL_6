# =============================================================================
# 02_prever_paciente.jl — o fluxo clínico proposto, de ponta a ponta
#
#   1. lê os biomarcadores pré-indução (arquivo TOML)
#   2. calibra o paciente virtual (Windkessel de 4 elementos: Rp e C)
#   3. prevê a curva de PAM para a dose escolhida, ANTES da injeção
#   4. relatório, faixa de incerteza, figuras e série batimento a batimento (CSV)
#
# Uso: julia --project=. -t auto scripts/02_prever_paciente.jl dados/paciente_exemplo.toml
# =============================================================================
using PropofolWindkessel
include(joinpath(@__DIR__, "graficos.jl"))

arquivo = isempty(ARGS) ? joinpath(@__DIR__, "..", "dados", "paciente_exemplo.toml") : ARGS[1]
pac, prot = ler_paciente(arquivo)
nome = splitext(basename(arquivo))[1]
dose = get(prot, "dose_mg_kg", 2.0)
duracao_bolus = get(prot, "duracao_bolus_s", 30.0)
manut = get(prot, "manutencao_ug_kg_min", 100.0)
modelo = Symbol(get(prot, "modelo_pk", "schnider"))

println("\n1) Biomarcadores pré-indução\n   ", pac)
println("\n2) Paciente virtual")
pv = calibrar(pac; modelo_pk = modelo)
println(pv)

println("\n3) Previsão para $(dose) mg/kg ($(round(dose * pac.peso; digits = 0)) mg) em $(duracao_bolus) s",
        manut > 0 ? " + manutenção $(manut) µg/kg/min" : "")
r = simular_inducao(pv; dose_mgkg = dose, duracao_bolus, manutencao_ugkgmin = manut, forma_de_onda = true)
mt = metricas(r)
relatorio(mt; titulo = "Previsão hemodinâmica — $(pac.nome)")

inc = previsao_com_incerteza(pv; n = 200, dose_mgkg = dose, duracao_bolus, manutencao_ugkgmin = manut)
@printf("\nIncerteza (%d sorteios de PK e de resposta hemodinâmica, hemodinâmica basal fixa):\n", inc.n)
@printf("   PAM mínima: mediana %.0f mmHg (faixa 90%%: %.0f–%.0f)\n", inc.nadir_mediano, inc.nadir_p05, inc.nadir_p95)
@printf("   P(PAM < 65) = %.0f%%   P(PAM < 55) = %.0f%%   P(BIS < 60) = %.0f%%\n",
        100 * inc.prob_pam65, 100 * inc.prob_pam55, 100 * inc.prob_hipnose)

titulo = @sprintf("%s | %.1f mg/kg em %.0f s | PAM %.0f → %.0f mmHg (risco %s)",
                  pac.nome, dose, duracao_bolus, mt.pam_basal, mt.pam_nadir, mt.risco)
savefig(fig_inducao(r; titulo), saida("$(nome)_inducao.png"))
savefig(fig_onda(r, mt), saida("$(nome)_onda.png"))
savefig(fig_incerteza(inc, r; titulo = @sprintf("%s — %.1f mg/kg: P(PAM<65) = %.0f%%, P(PAM<55) = %.0f%%",
                                                 pac.nome, dose, 100 * inc.prob_pam65, 100 * inc.prob_pam55)),
        saida("$(nome)_incerteza.png"))

b = r.bat
open(saida("$(nome)_batimentos.csv"), "w") do io
    println(io, "t_min,Ce,BIS,PAS,PAD,PAM,FC,VS,DC,Rp")
    for j in eachindex(b.t_min)
        @printf(io, "%.4f,%.3f,%.1f,%.2f,%.2f,%.2f,%.2f,%.2f,%.3f,%.4f\n", b.t_min[j] - r.t_inducao, b.Ce[j], b.BIS[j],
                b.PAS[j], b.PAD[j], b.PAM[j], b.FC[j], b.VS[j], b.DC[j], b.Rp[j])
    end
end
println("\nFiguras e CSV em: ", normpath(DIR_RESULTADOS))
