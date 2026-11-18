from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / "resultados" / "vitaldb" / "predicoes"
OUTPUT = ROOT / "resultados" / "vitaldb" / "analise_fases"
PHASES = (("0_5min", 0.0, 300.0), ("5_10min", 300.0, 600.0), ("10_20min", 600.0, 1200.0))
THRESHOLD = 65.0
MAX_GAP_S = 5.0


def hypotension_duration(t, y, threshold=THRESHOLD, max_gap=MAX_GAP_S):
    """Tempo abaixo do limiar em trechos adjacentes validos, com cruzamento linear.
    Intervalos sem dados >max_gap nao sao preenchidos.
    """
    total = observed_span = 0.0
    for i in range(len(t) - 1):
        if not (np.isfinite(t[i]) and np.isfinite(t[i+1]) and
                np.isfinite(y[i]) and np.isfinite(y[i+1])):
            continue
        dt = t[i+1] - t[i]
        if not (0 < dt <= max_gap):
            continue
        observed_span += dt
        a, b = y[i], y[i+1]
        if a < threshold and b < threshold:
            total += dt
        elif (a < threshold) != (b < threshold):
            frac = np.clip((threshold - a) / (b - a), 0, 1)
            total += dt * (frac if a < threshold else (1 - frac))
    return total, observed_span


def main():
    paths = sorted(INPUT.glob("case_*_padrao.csv"))
    if not paths:
        raise SystemExit(f"Nao encontrei arquivos case_*_padrao.csv em {INPUT}. Execute a validacao Julia primeiro.")
    OUTPUT.mkdir(parents=True, exist_ok=True)
    rows = []
    required = ("time_after_induction_s", "pam_observed_mmHg", "pam_predicted_mmHg")
    for path in paths:
        d = pd.read_csv(path)
        if any(c not in d.columns for c in required):
            print(f"IGNORADO: {path.name} (colunas ausentes)")
            continue
        cid = int(path.stem.split("_")[1])
        d = d.loc[:, required].apply(pd.to_numeric, errors="coerce")
        d = d.replace([np.inf, -np.inf], np.nan).sort_values(required[0])
        t_all = d[required[0]].to_numpy(dtype=float)
        if pd.Series(t_all).dropna().duplicated().any():
            print(f"AVISO caso {cid}: instantes repetidos; revisar antes de interpretar.")
        obs_all = d[required[1]].to_numpy(dtype=float)
        pred_all = d[required[2]].to_numpy(dtype=float)
        for name, start, end in PHASES:
            # Inclusao no inicio e exclusao no fim, exceto a ultima fase.
            in_phase = (t_all >= start) & (t_all < end if end < 1200 else t_all <= end)
            t, o, p = t_all[in_phase], obs_all[in_phase], pred_all[in_phase]
            pair = np.isfinite(t) & np.isfinite(o) & np.isfinite(p)
            err = p[pair] - o[pair]
            mae = float(np.mean(np.abs(err))) if len(err) else np.nan
            rmse = float(np.sqrt(np.mean(err**2))) if len(err) else np.nan
            bias = float(np.mean(err)) if len(err) else np.nan
            valid_o, valid_p = np.isfinite(o), np.isfinite(p)
            obs_duration, obs_covered = hypotension_duration(t, o)
            pred_duration, pred_covered = hypotension_duration(t, p)
            rows.append({
                "caseid": cid, "fase": name, "inicio_s": start, "fim_s": end,
                "n_pares": int(pair.sum()),
                "mae_mmHg": mae, "rmse_mmHg": rmse, "vies_mmHg": bias,
                "pam_min_observada_mmHg": float(np.min(o[valid_o])) if valid_o.any() else np.nan,
                "pam_min_prevista_mmHg": float(np.min(p[valid_p])) if valid_p.any() else np.nan,
                "tempo_nadir_obs_s": float(t[valid_o][np.argmin(o[valid_o])]) if valid_o.any() else np.nan,
                "tempo_nadir_pred_s": float(t[valid_p][np.argmin(p[valid_p])]) if valid_p.any() else np.nan,
                "hipotensao_obs_s": obs_duration,
                "hipotensao_pred_s": pred_duration,
                "cobertura_obs_pct": 100 * obs_covered / (end-start),
                "cobertura_pred_pct": 100 * pred_covered / (end-start),
            })
    if not rows:
        raise SystemExit("Nenhum caso processavel.")
    metrics = pd.DataFrame(rows)
    metrics.to_csv(OUTPUT / "metricas_por_fase.csv", index=False)
    summary = (metrics.groupby("fase", as_index=False)
               .agg(n_pacientes=("caseid", "nunique"),
                    media_mae_mmHg=("mae_mmHg", "mean"),
                    mediana_mae_mmHg=("mae_mmHg", "median"),
                    media_rmse_mmHg=("rmse_mmHg", "mean"),
                    mediana_rmse_mmHg=("rmse_mmHg", "median"),
                    menor_cobertura_obs_pct=("cobertura_obs_pct", "min")))
    summary["fase"] = pd.Categorical(summary.fase, categories=[x[0] for x in PHASES], ordered=True)
    summary = summary.sort_values("fase")
    summary.to_csv(OUTPUT / "resumo_fases.csv", index=False)
    print("\nRESUMO DAS FASES (metricas agregadas por paciente):")
    print(summary.to_string(index=False, float_format=lambda x: f"{x:.2f}"))
    print("\nATENCAO: os tempos abaixo de 65 mmHg foram estimados por interpolacao linear somente em lacunas <=5 s.")
    print("Os nadirs sao os minimos das amostras validas; nao ha remocao automatica de artefatos.")
    try:
        import matplotlib.pyplot as plt
        ax = (metrics.pivot(index="caseid", columns="fase", values="rmse_mmHg")
              .reindex(columns=[x[0] for x in PHASES]).plot(kind="bar", figsize=(13, 5)))
        ax.set(xlabel="Caso VitalDB", ylabel="RMSE da PAM (mmHg)",
               title="Erro por fase temporal (intervalos sem sobreposicao)")
        
        ax.legend(title="Fase")
        plt.tight_layout()
        plt.savefig(OUTPUT / "rmse_por_fase.png", dpi=180)
        plt.close()
    except ImportError:
        print("Matplotlib ausente: CSVs gerados; grafico omitido.")
    print(f"\nArquivos gerados: {OUTPUT}")


if __name__ == "__main__":
    main()
