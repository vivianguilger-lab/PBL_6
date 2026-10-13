# =============================================================================
# 02_prever_paciente.jl — o fluxo clínico proposto, de ponta a ponta
#
#   1. lê os dados pré-indução do paciente (arquivo TOML)
#   2. calibra o paciente virtual (Windkessel + coração + veias)
#   3. prevê a curva de PAM para a dose escolhida, ANTES da injeção
#   4. gera relatório, figuras e a série batimento a batimento (CSV)
#
# Uso:
#   julia --project=. scripts/02_prever_paciente.jl dados/paciente_exemplo.toml
# =============================================================================
using PropofolWindkessel
include(joinpath(@__DIR__, "graficos.jl"))

arquivo = isempty(ARGS) ? joinpath(@__DIR__, "..", "dados", "paciente_exemplo.toml") : ARGS[1]
pac, prot = ler_paciente(arquivo)
nome = splitext(basename(arquivo))[1]

dose = get(prot, "dose_mg_kg", 2.0)
duracao_bolus = get(prot, "duracao_bolus_s", 20.0)
manut = get(prot, "manutencao_ug_kg_min", 0.0)

println("\n1) Dados pré-indução\n   ", pac)

println("\n2) Calibração do paciente virtual")
pv = calibrar(pac)
println(pv)

println("\n3) Previsão para $(dose) mg/kg ($(round(dose * pac.peso; digits = 0)) mg) em $(duracao_bolus) s",
        manut > 0 ? " + manutenção $(manut) µg/kg/min" : "")
r = simular_inducao(pv; dose_mgkg = dose, duracao_bolus, manutencao_ugkgmin = manut)
mt = metricas(r)
relatorio(mt; titulo = "Previsão hemodinâmica — $(pac.nome)")

# --- incerteza: mesma hemodinâmica, farmacologia sorteada --------------------
inc = previsao_com_incerteza(pv; n = 200, dose_mgkg = dose, duracao_bolus, manutencao_ugkgmin = manut)
@printf("\nIncerteza farmacológica (%d sorteios de PK/PD, hemodinâmica do paciente fixa):\n", inc.n)
@printf("   PAM nadir: mediana %.0f mmHg (faixa 90%%: %.0f–%.0f)\n", inc.nadir_mediano, inc.nadir_p05, inc.nadir_p95)
@printf("   P(PAM < 65) = %.0f%%   P(PAM < 55) = %.0f%%   P(BIS < 60) = %.0f%%\n",
        100 * inc.prob_pam65, 100 * inc.prob_pam55, 100 * inc.prob_hipnose)

# --- figuras --------------------------------------------------------------
titulo = @sprintf("%s | %.1f mg/kg em %.0f s | PAM %.0f → %.0f mmHg (risco %s)",
                  pac.nome, dose, duracao_bolus, mt.pam_basal, mt.pam_nadir, mt.risco)
savefig(fig_inducao(r, pv; titulo), saida("$(nome)_inducao.png"))
savefig(fig_ondas(r, mt), saida("$(nome)_ondas.png"))
savefig(fig_incerteza(inc, r; titulo = @sprintf("%s — %.1f mg/kg: P(PAM<65) = %.0f%%, P(PAM<55) = %.0f%%",
                                                 pac.nome, dose, 100 * inc.prob_pam65, 100 * inc.prob_pam55)),
        saida("$(nome)_incerteza.png"))
fg, g0, g1 = fig_guyton(pv, r, mt)
savefig(fg, saida("$(nome)_guyton.png"))

println("\n4) Leitura pelo ponto de funcionamento (Resumo, Seção 16)")
@printf("   Pms:  %.1f → %.1f mmHg   (venodilatação: volume estressado ↓)\n", g0.Pms, g1.Pms)
@printf("   C_S:  %.2f → %.2f ml/mmHg (depressão miocárdica: curva cardíaca achata)\n", g0.CS, g1.CS)
@printf("   R_VR: %.3f → %.3f mmHg·s/ml (TPR entra × AC/(AC+VC) ≈ %.3f)\n", g0.RVR, g1.RVR, pv.cv.Ca / (pv.cv.Ca + pv.cv.Cv))
@printf("   DC no ponto de funcionamento: %.2f → %.2f L/min\n", g0.Q, g1.Q)

# --- série batimento a batimento --------------------------------------------
b = r.bat
open(saida("$(nome)_batimentos.csv"), "w") do io
    println(io, "t_min,PAS,PAD,PAM,FC,DC,VS")
    for j in eachindex(b.t)
        @printf(io, "%.4f,%.2f,%.2f,%.2f,%.2f,%.3f,%.2f\n", (b.t[j] - r.t_inducao) / 60,
                b.PAS[j], b.PAD[j], b.PAM[j], b.FC[j], b.DC[j], b.VS[j])
    end
end
println("\nFiguras e CSV em: ", normpath(DIR_RESULTADOS))
