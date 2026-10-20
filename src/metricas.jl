# =============================================================================
# metricas.jl — métricas de hipotensão pós-indução e relatório
# =============================================================================
#
# Hipotensão pós-indução: queda de PA nos primeiros 20 min após a indução.
# 17–24% dos adultos têm PAM < 55 mmHg; PAM ≤ 55 mmHg sustentada associa-se a
# lesão miocárdica e lesão renal aguda. Por isso o relatório usa dois limiares
# (65 e 55 mmHg) e a área abaixo do limiar (AUT, mmHg·min).
# =============================================================================

"""
    metricas(r::Resultado; janela=20.0, limiares=(65.0, 55.0), alvo_bis=60.0)

Resumo clínico da curva prevista (batimentos nos `janela` minutos após o bolus).
"""
function metricas(r::Resultado; janela = 20.0, limiares = (65.0, 55.0), alvo_bis = 60.0)
    b, k = r.bat, r.pkpd
    t0 = r.t_inducao
    base = findall(t -> t < t0, b.t_min)
    media(v, idx) = isempty(idx) ? v[1] : mean(v[idx])
    pam0, fc0, dc0 = media(b.PAM, base), media(b.FC, base), media(b.DC, base)
    jan = findall(t -> t0 <= t <= t0 + janela, b.t_min)
    pam = b.PAM[jan]
    dtb = 1 ./ b.FC[jan]                              # duração de cada batimento [min]
    i = argmin(pam)
    j = jan[i]
    tempo_abaixo = Dict(l => sum(dtb[pam .< l]; init = 0.0) for l in limiares)
    aut = Dict(l => sum(max.(l .- pam, 0.0) .* dtb; init = 0.0) for l in limiares)
    js = findall(t -> t0 <= t <= t0 + janela, k.t)
    bis = k.bis[js]
    ib = findfirst(<(alvo_bis), bis)
    dt = length(k.t) > 1 ? k.t[2] - k.t[1] : 0.0
    lmin, lmax = minimum(limiares), maximum(limiares)
    risco = tempo_abaixo[lmin] > 0 ? "ALTO" :
            (tempo_abaixo[lmax] >= 1.0 || pam[i] < 0.7 * pam0) ? "MODERADO" : "BAIXO"
    return (pam_basal = pam0, fc_basal = fc0, dc_basal = dc0,
            pam_nadir = pam[i], t_nadir = b.t_min[j] - t0, queda_pct = 100 * (1 - pam[i] / pam0),
            pas_nadir = b.PAS[j], pad_nadir = b.PAD[j], fc_no_nadir = b.FC[j], dc_no_nadir = b.DC[j],
            rp_no_nadir = b.Rp[j], vs_no_nadir = b.VS[j],
            tempo_abaixo = tempo_abaixo, aut = aut,
            bis_min = minimum(bis), t_hipnose = ib === nothing ? NaN : k.t[js][ib] - t0,
            tempo_hipnose = count(<(alvo_bis), bis) * dt,
            cp_max = maximum(k.Cp), ce_max = maximum(k.Ce),
            risco = risco, limiares = limiares, alvo_bis = alvo_bis)
end

"""Imprime um relatório legível das métricas."""
function relatorio(io::IO, mt; titulo = "Previsão hemodinâmica")
    println(io, "─"^72)
    println(io, titulo)
    println(io, "─"^72)
    @printf(io, "PAM basal            : %6.1f mmHg   (FC %.0f bpm, DC %.2f L/min)\n", mt.pam_basal, mt.fc_basal, mt.dc_basal)
    @printf(io, "PAM mínima           : %6.1f mmHg   em %.1f min (queda de %.0f%%)\n", mt.pam_nadir, mt.t_nadir, mt.queda_pct)
    @printf(io, "  nesse momento      : PA %.0f/%.0f, FC %.0f bpm, DC %.2f L/min\n", mt.pas_nadir, mt.pad_nadir, mt.fc_no_nadir, mt.dc_no_nadir)
    for l in sort(collect(mt.limiares); rev = true)
        @printf(io, "Tempo com PAM < %2.0f   : %6.1f min    AUT = %.0f mmHg·min\n", l, mt.tempo_abaixo[l], mt.aut[l])
    end
    @printf(io, "BIS mínimo           : %6.0f        (BIS < %.0f em %.1f min; por %.1f min)\n",
            mt.bis_min, mt.alvo_bis, mt.t_hipnose, mt.tempo_hipnose)
    @printf(io, "Ce máxima            : %.2f µg/ml\n", mt.ce_max)
    println(io, "Risco previsto       : ", mt.risco)
    println(io, "─"^72)
end
relatorio(mt; kw...) = relatorio(stdout, mt; kw...)
