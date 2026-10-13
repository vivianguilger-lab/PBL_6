# =============================================================================
# incerteza.jl — previsão individual com faixa de incerteza
# =============================================================================
#
# O paciente virtual tem a SUA hemodinâmica (calibrada com os dados
# pré-indução), mas a sensibilidade dele ao propofol é desconhecida antes da
# injeção. Mantemos a calibração e sorteamos a farmacologia (PK de Schnider e
# PD de Su/Bouillon + parâmetros assumidos) com a variabilidade
# interindividual da população:
#
#   curva prevista = mediana;  faixa = percentis 5–95%;
#   risco = fração das simulações com PAM < 65 / < 55 mmHg.
# =============================================================================

"""
    previsao_com_incerteza(pv; n=200, rng=Xoshiro(1), dt=6.0, kw...)

Roda `n` simulações do mesmo paciente virtual, cada uma com PK/PD sorteadas.
`kw` é repassado a `simular_inducao` (dose_mgkg, duracao_bolus, ...).
Devolve a grade de tempo (min), mediana e percentis 5/95 da PAM, e as
probabilidades de hipotensão.
"""
function previsao_com_incerteza(pv::PacienteVirtual; n = 200, rng = Xoshiro(1), dt = 6.0,
                                janela = 20 * 60.0, kw...)
    amostras = [(amostrar_eta_pk(rng), amostrar_pd(rng, pv.pd)) for _ in 1:n]
    t_inducao = 60.0
    grade = collect(0.0:dt:(t_inducao + janela))
    curvas = zeros(length(grade), n)
    nadir = zeros(n); abaixo65 = falses(n); abaixo55 = falses(n); bis_ok = falses(n)
    Threads.@threads for k in 1:n
        eta, pd = amostras[k]
        pk = schnider(pv.paciente; eta)
        r = simular_inducao(pv; pk, pd, t_inducao, duracao = janela, saveat = 0.01, kw...)
        mt = metricas(r; janela)
        curvas[:, k] = interpolar(r.bat.t, r.bat.PAM, grade)
        nadir[k] = mt.pam_nadir
        abaixo65[k] = mt.tempo_abaixo[65.0] > 0
        abaixo55[k] = mt.tempo_abaixo[55.0] > 0
        bis_ok[k] = mt.bis_min < 60
    end
    q(p) = [quantile(view(curvas, i, :), p) for i in eachindex(grade)]
    return (t = (grade .- t_inducao) ./ 60, mediana = q(0.5), p05 = q(0.05), p95 = q(0.95),
            nadir_mediano = median(nadir), nadir_p05 = quantile(nadir, 0.05), nadir_p95 = quantile(nadir, 0.95),
            prob_pam65 = mean(abaixo65), prob_pam55 = mean(abaixo55), prob_hipnose = mean(bis_ok), n = n)
end

"""Interpolação linear de (x, y) nos pontos `xq` (extrapola constante)."""
function interpolar(x, y, xq)
    out = similar(xq, Float64)
    j = 1
    for (i, v) in enumerate(xq)
        if v <= x[1]
            out[i] = y[1]
        elseif v >= x[end]
            out[i] = y[end]
        else
            while x[j+1] < v
                j += 1
            end
            w = (v - x[j]) / (x[j+1] - x[j])
            out[i] = y[j] + w * (y[j+1] - y[j])
        end
    end
    return out
end
