# =============================================================================
# incerteza.jl — previsão individual com faixa de incerteza
# =============================================================================
#
# O paciente virtual tem a SUA hemodinâmica basal (calibrada com os dados
# pré-indução), mas a sensibilidade dele ao propofol é desconhecida antes da
# injeção. Mantemos a calibração e sorteamos a PK (Schnider) e a resposta
# hemodinâmica (Ce50 e Imax do modelo do grupo) com a variabilidade da população.
# =============================================================================

"""
    previsao_com_incerteza(pv; n=200, rng=Xoshiro(1), dt=0.1, janela=20.0, kw...)

Roda `n` simulações do mesmo paciente virtual com PK/PD sorteadas. Devolve a
grade de tempo (min após o bolus), mediana e percentis 5/95 da PAM e as
probabilidades de hipotensão.
"""
function previsao_com_incerteza(pv::PacienteVirtual; n = 200, rng = Xoshiro(1), dt = 0.1, janela = 20.0, kw...)
    amostras = [(amostrar_eta_pk(rng), amostrar_efeito(rng, pv.efeito)) for _ in 1:n]
    t_inducao = 1.0
    grade = collect(0.0:dt:(t_inducao + janela))
    curvas = zeros(length(grade), n)
    nadir = zeros(n); abaixo65 = falses(n); abaixo55 = falses(n); bis_ok = falses(n)
    Threads.@threads for k in 1:n
        eta, ef = amostras[k]
        pk = parametros_pk(pv.paciente, pv.pk.modelo; eta)
        r = simular_inducao(pv; pk, efeito = ef, t_inducao, duracao = janela, kw...)
        mt = metricas(r; janela)
        curvas[:, k] = interpolar(r.bat.t_min, r.bat.PAM, grade)
        nadir[k] = mt.pam_nadir
        abaixo65[k] = mt.tempo_abaixo[65.0] > 0
        abaixo55[k] = mt.tempo_abaixo[55.0] > 0
        bis_ok[k] = mt.bis_min < 60
    end
    q(p) = [quantile(view(curvas, i, :), p) for i in eachindex(grade)]
    return (t = grade .- t_inducao, mediana = q(0.5), p05 = q(0.05), p95 = q(0.95),
            nadir_mediano = median(nadir), nadir_p05 = quantile(nadir, 0.05), nadir_p95 = quantile(nadir, 0.95),
            prob_pam65 = mean(abaixo65), prob_pam55 = mean(abaixo55), prob_hipnose = mean(bis_ok), n = n)
end

"""Interpolação linear de (x, y) nos pontos `xq` (extrapola constante)."""
interpolar(x, y, xq) = [interp(x, y, v) for v in xq]
