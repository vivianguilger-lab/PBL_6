# =============================================================================
# 03_janela_de_dose.jl — "testar a dose no computador antes da injeção"
#
# Varre doses de 0,5 a 3,0 mg/kg, calcula a janela individual (hipnose × PAM ≥ 65)
# e compara velocidades de injeção.
#
# Uso: julia --project=. scripts/03_janela_de_dose.jl dados/paciente_idoso.toml
#
# Entrada: arquivo .toml do paciente (sem argumento, usa dados/paciente_exemplo.toml).
# Da seção [protocolo] usa: limiar_pam, alvo_bis, duracao_bolus_s,
# manutencao_ug_kg_min, dose_mg_kg e modelo_pk.
# Saída na tela: tabela dose × PAM mínima × tempo abaixo de 65/55 × BIS mínimo,
# a janela individual e o efeito da velocidade de injeção.
# Saída em resultados/: <nome>_doses.png e <nome>_velocidade.png.
# Janela = da menor dose que leva o BIS abaixo do alvo até a maior dose que
# mantém a PAM mínima ≥ limiar. O critério da dose mínima ainda está em revisão
# (docs/validacao.md, seção 6) e os parâmetros hemodinâmicos são ilustrativos.
# =============================================================================
# carrega o pacote do projeto e as funções de figura (graficos.jl)
using PropofolWindkessel
include(joinpath(@__DIR__, "graficos.jl"))

# lê o paciente e o protocolo (mesma lógica de 02_prever_paciente.jl)
arquivo = isempty(ARGS) ? joinpath(@__DIR__, "..", "dados", "paciente_exemplo.toml") : ARGS[1]
pac, prot = ler_paciente(arquivo)
nome = splitext(basename(arquivo))[1]
# valores do protocolo; get(prot, "chave", padrão) usa o padrão se faltar no .toml
limiar = get(prot, "limiar_pam", 65.0)
alvo_bis = get(prot, "alvo_bis", 60.0)
duracao_bolus = get(prot, "duracao_bolus_s", 30.0)
manut = get(prot, "manutencao_ug_kg_min", 100.0)
dprot = get(prot, "dose_mg_kg", 2.0)

# paciente virtual calibrado (Rp e C) com o modelo PK escolhido
pv = calibrar(pac; modelo_pk = Symbol(get(prot, "modelo_pk", "schnider")))
println(pac)
# Simula as doses de 0,5 a 3,0 mg/kg, de 0,25 em 0,25 (0.5:0.25:3.0).
# Cada item de `tab` guarda as métricas de uma dose.
tab = varrer_doses(pv, 0.5:0.25:3.0; duracao_bolus, manutencao_ugkgmin = manut)
# Procura os limites da janela por bisseção (divide o intervalo ao meio até
# achar o limite). `alvo_bis,` sozinho equivale a alvo_bis = alvo_bis.
jan = janela_de_dose(pv; limiar_pam = limiar, alvo_bis, duracao_bolus, manutencao_ugkgmin = manut)

# Tabela: uma linha por dose. Em @printf, `%%` imprime o próprio símbolo %.
# t.tempo_abaixo[65.0] = minutos com PAM < 65 nos 20 min após o bolus.
println("\n dose    mg    PAM mín.  queda   t mín.   PAM<65  PAM<55  BIS mín  risco")
for t in tab
    @printf(" %.2f  %5.0f   %6.1f    %4.0f%%   %4.1f min  %4.1f    %4.1f    %5.0f    %s\n",
            t.dose, t.dose * pac.peso, t.pam_nadir, t.queda_pct, t.t_nadir, t.tempo_abaixo[65.0],
            t.tempo_abaixo[55.0], t.bis_min, t.risco)
end
println()
# Janela: se existe, mostra a faixa em mg/kg e em mg; se não existe, a dose
# que já dá hipnose baixa a PAM abaixo do limiar.
if jan.existe
    @printf("Janela individual (bolus em %.0f s): %.2f a %.2f mg/kg  (%.0f a %.0f mg)\n",
            duracao_bolus, jan.dose_min, jan.dose_max, jan.dose_min * pac.peso, jan.dose_max * pac.peso)
else
    @printf("Sem janela: hipnose exige %.2f mg/kg, mas PAM < %.0f já com %.2f mg/kg.\n", jan.dose_min, limiar, jan.dose_max)
end

# Velocidade de injeção: a mesma dose do protocolo injetada em 10, 60 e 120 s.
# Cria um gráfico vazio (pV) e acrescenta uma curva de PAM por velocidade.
println("\nVelocidade de injeção ($(dprot) mg/kg):")
pV = plot(; title = "PAM prevista × velocidade de injeção ($(dprot) mg/kg)", xlabel = "minutos após o início",
          ylabel = "PAM (mmHg)", legend = :bottomright, size = (800, 450))
# enumerate dá o número da volta (k = 1, 2, 3) e o valor (dur, em segundos)
for (k, dur) in enumerate((10.0, 60.0, 120.0))
    r = simular_inducao(pv; dose_mgkg = dprot, duracao_bolus = dur, manutencao_ugkgmin = manut)
    mt = metricas(r)
    @printf("  bolus em %3.0f s: PAM mínima %.1f mmHg (%.1f min), BIS mín %.0f\n", dur, mt.pam_nadir, mt.t_nadir, mt.bis_min)
    # minutos(r, ...) conta o tempo a partir do início do bolus
    plot!(pV, minutos(r, r.bat.t_min), r.bat.PAM; color = (COR.s1, COR.s2, COR.s3)[k], label = "injeção em $(Int(dur)) s")
end
# linhas tracejadas de 65 e 55 mmHg (função de graficos.jl)
limiares!(pV, 20.0)
# grava as duas figuras em resultados/
savefig(fig_doses(tab, jan; titulo = "$(pac.nome) — janela de dose individual"), saida("$(nome)_doses.png"))
savefig(pV, saida("$(nome)_velocidade.png"))
println("\nFiguras em: ", normpath(DIR_RESULTADOS))
