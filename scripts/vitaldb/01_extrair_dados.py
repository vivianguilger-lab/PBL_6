#!/usr/bin/env python3
"""Seleciona casos VitalDB com dados observados e baixa somente as trilhas necessárias.

Usa o clinical_data.csv local e a API pública VitalDB para carregar, por caso,
apenas as trilhas selecionadas. Não baixa o ZIP completo nem escreve no código-fonte.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

try:
    import vitaldb
    try:
        from vitaldb import dataset as vitaldb_dataset
        API_DATASET_VERSION = str(getattr(vitaldb_dataset, "DATASET_VERSION", "desconhecida"))
    except Exception:
        API_DATASET_VERSION = "desconhecida"
except ImportError as exc:
    raise SystemExit(
        "Biblioteca vitaldb ausente. Execute: python -m pip install -r scripts/vitaldb/requirements.txt"
    ) from exc

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CLINICAL = ROOT / "dados" / "vitaldb" / "metadados" / "clinical_data.csv"
OUT = ROOT / "dados" / "vitaldb" / "processados"
CASES_OUT = OUT / "casos"

REQUIRED_TRACKS = [
    "Solar8000/ART_MBP",
    "Solar8000/ART_SBP",
    "Solar8000/ART_DBP",
    "Solar8000/HR",
    "Orchestra/PPF20_RATE",
]
CO_TRACKS = ["Vigileo/CO", "EV1000/CO", "Vigilance/CO", "CardioQ/CO"]


def numeric(series: pd.Series) -> pd.Series:
    # Valores censurados como '>89' viram ausentes; não são convertidos em 89.
    return pd.to_numeric(series, errors="coerce")


def normalized_sex(value) -> str | None:
    if pd.isna(value):
        return None
    value = str(value).strip().upper()
    if value in {"M", "MALE"}:
        return "M"
    if value in {"F", "FEMALE"}:
        return "F"
    return None


def get_case_track_map(case_ids: list[int]) -> dict[int, set[str]]:
    """Consulta apenas o índice de trilhas; sinais são carregados depois, por caso."""
    print("Consultando o índice público de trilhas do VitalDB...")
    try:
        tracks_df = vitaldb.get_track_names(caseids=case_ids)
    except Exception as exc:
        raise RuntimeError(
            "Não foi possível consultar o índice de trilhas da API VitalDB. "
            "Verifique a internet e tente novamente."
        ) from exc

    result: dict[int, set[str]] = {}
    if tracks_df is None or tracks_df.empty:
        return result
    for row in tracks_df.itertuples(index=False):
        cid = int(row.caseid)
        names = getattr(row, "tnames", [])
        if isinstance(names, str):
            # Compatibilidade com versões que possam retornar texto em vez de lista.
            names = [part.strip() for part in names.split(",") if part.strip()]
        result[cid] = set(str(name).strip() for name in names)
    return result


def ffill_short(values: np.ndarray, max_gap: int = 5) -> np.ndarray:
    """Mantém o último ajuste da bomba por no máximo max_gap segundos."""
    return pd.Series(values, dtype="float64").ffill(limit=max_gap).to_numpy(dtype=float)


def first_induction(rate_ml_h: np.ndarray, min_baseline_s: int = 60,
                    sustained_s: int = 10, min_positive_s: int = 5) -> int | None:
    """Primeiro início sustentado de propofol após 60 s de taxa explicitamente zero."""
    n = len(rate_ml_h)
    last = n - 1200 - sustained_s
    for i in range(min_baseline_s, max(min_baseline_s, last + 1)):
        baseline = rate_ml_h[i - min_baseline_s:i]
        future = rate_ml_h[i:i + sustained_s]
        if len(future) != sustained_s:
            break
        if not np.isfinite(baseline).all() or not np.isfinite(future).all():
            continue
        if np.max(baseline) <= 0.0 and np.count_nonzero(future > 0.0) >= min_positive_s:
            return i
    return None


def clean_physiology(values: np.ndarray, low: float, high: float) -> np.ndarray:
    out = values.astype(float, copy=True)
    out[~np.isfinite(out)] = np.nan
    out[(out < low) | (out > high)] = np.nan
    return out


def valid_count(values: np.ndarray) -> int:
    return int(np.isfinite(values).sum())


def extract_one(row: pd.Series, track_names: set[str]) -> tuple[dict | None, str]:
    cid = int(row["caseid"])
    co_track = next((name for name in CO_TRACKS if name in track_names), None)
    if co_track is None:
        return None, "sem_trilha_de_debito_cardiaco_medido"

    mandatory = list(REQUIRED_TRACKS) + [co_track]
    if not all(name in track_names for name in mandatory):
        missing = [name for name in mandatory if name not in track_names]
        return None, "trilhas_obrigatorias_ausentes:" + ";".join(missing)

    # Todas as trilhas de taxa disponíveis no registro são consultadas para sinalizar
    # coadministração registrada. Ausência de trilha nunca é interpretada como taxa zero.
    rate_tracks = sorted(
        name for name in track_names
        if name.startswith("Orchestra/") and name.endswith("_RATE")
    )
    requested = list(dict.fromkeys(mandatory + rate_tracks))
    try:
        arr = vitaldb.load_case(cid, requested, interval=1)
    except Exception as exc:
        return None, f"erro_api:{type(exc).__name__}"
    if arr is None:
        return None, "api_sem_dados"
    arr = np.asarray(arr, dtype=float)
    if arr.ndim != 2 or arr.shape[1] != len(requested) or arr.shape[0] < 1261:
        return None, "formato_ou_duracao_insuficiente"
    col = {name: arr[:, i] for i, name in enumerate(requested)}

    rate_raw = col["Orchestra/PPF20_RATE"]
    rate_ml_h = ffill_short(rate_raw, max_gap=5)
    onset = first_induction(rate_ml_h)
    if onset is None:
        return None, "sem_inducao_com_60s_de_baseline_e_20min_de_seguimento"
    start = onset - 60
    end = onset + 1200
    if start < 0 or end >= arr.shape[0]:
        return None, "janela_temporal_incompleta"

    window_rate = rate_ml_h[start:end + 1]
    if not np.isfinite(window_rate).all():
        return None, "lacunas_na_taxa_de_propofol_maiores_que_5s"
    if np.max(window_rate[:60]) > 0.0:
        return None, "baseline_com_propofol_ativo"

    # Sinais fisiológicos: valores implausíveis viram ausentes; não são interpolados.
    pam_all = clean_physiology(col["Solar8000/ART_MBP"], 20, 220)
    pas_all = clean_physiology(col["Solar8000/ART_SBP"], 40, 300)
    pad_all = clean_physiology(col["Solar8000/ART_DBP"], 20, 200)
    hr_all = clean_physiology(col["Solar8000/HR"], 20, 250)
    co_all = clean_physiology(col[co_track], 0.5, 20)

    base = slice(start, onset)
    baselines = {
        "pam": pam_all[base], "pas": pas_all[base], "pad": pad_all[base],
        "hr": hr_all[base], "co": co_all[base],
    }
    minimum_counts = {"pam": 10, "pas": 10, "pad": 10, "hr": 10, "co": 3}
    for name, vals in baselines.items():
        if valid_count(vals) < minimum_counts[name]:
            return None, f"baseline_insuficiente_{name}"
    med = {name: float(np.nanmedian(vals)) for name, vals in baselines.items()}
    if not (med["pas"] > med["pam"] > med["pad"]):
        return None, "relacao_baseline_PAS_PAM_PAD_invalida"
    if med["co"] <= 0 or med["hr"] <= 0:
        return None, "baseline_CO_ou_FC_invalido"

    # Avalia somente as trilhas de taxa de outros fármacos que realmente existem.
    other_rate_tracks = [t for t in rate_tracks if "PPF" not in t.upper() and "PROPOFOL" not in t.upper()]
    active_other: list[str] = []
    unknown_other: list[str] = []
    for name in other_rate_tracks:
        values = ffill_short(col[name], max_gap=5)[start:end + 1]
        coverage = np.isfinite(values).mean()
        if coverage < 0.95:
            unknown_other.append(name)
        elif np.any(values[np.isfinite(values)] > 0.0):
            active_other.append(name)
    if active_other:
        coadmin_status = "active"
    elif not other_rate_tracks or unknown_other:
        coadmin_status = "unknown"
    else:
        coadmin_status = "none_recorded"

    local_t = np.arange(start, end + 1, dtype=int)
    n = len(local_t)
    output = pd.DataFrame({
        "t_case_s": local_t,
        "t_rel_s": local_t - start,
        "pam_obs_mmHg": pam_all[start:end + 1],
        "pas_obs_mmHg": pas_all[start:end + 1],
        "pad_obs_mmHg": pad_all[start:end + 1],
        "hr_obs_bpm": hr_all[start:end + 1],
        "co_obs_l_min": co_all[start:end + 1],
        "propofol_rate_ml_h": window_rate,
        # PPF20_RATE é mL/h de solução a 20 mg/mL: mL/h × 20 mg/mL ÷ 60.
        "propofol_rate_mg_min": window_rate * 20.0 / 60.0,
    })
    # Preserva todas as taxas de outras bombas disponíveis, inclusive remifentanil.
    # Valores não registrados permanecem NaN; nunca são interpretados como zero.
    # A taxa original é mantida em mL/h, sem presumir concentração farmacológica.
    for track in other_rate_tracks:
        colname = "other_rate_ml_h__" + track.replace("/", "__")
        output[colname] = ffill_short(col[track], max_gap=5)[start:end + 1]

    if len(output) != n:
        return None, "erro_interno_tamanho_janela"

    CASES_OUT.mkdir(parents=True, exist_ok=True)
    case_path = CASES_OUT / f"case_{cid:05d}.csv"
    output.to_csv(case_path, index=False, na_rep="")

    meta = {
        "caseid": cid,
        "clinical_metadata_source": "PhysioNet clinical_data.csv local (versão 1.0.0 informada no download)",
        "signal_api_dataset_version": API_DATASET_VERSION,
        "age": float(row["age"]),
        "sex": str(row["sex_norm"]),
        "weight_kg": float(row["weight"]),
        "height_cm": float(row["height"]),
        "asa": int(row["asa_num"]),
        "onset_case_s": int(onset),
        "window_start_case_s": int(start),
        "co_track": co_track,
        "coadmin_status": coadmin_status,
        "active_other_rate_tracks": ";".join(active_other),
        "unknown_other_rate_tracks": ";".join(unknown_other),
        "other_rate_track_count": len(other_rate_tracks),
        "propofol_rate_fill_method": "forward-fill <=5s only; longer gaps rejected",
        "data_file": str(case_path.relative_to(ROOT)).replace("\\", "/"),
        "n_rows": n,
        "baseline_pam_median_mmHg": med["pam"],
        "baseline_pas_median_mmHg": med["pas"],
        "baseline_pad_median_mmHg": med["pad"],
        "baseline_hr_median_bpm": med["hr"],
        "baseline_co_median_l_min": med["co"],
    }
    return meta, "ok"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--clinical-csv", type=Path, default=DEFAULT_CLINICAL)
    parser.add_argument("--max-cases", type=int, default=30,
                        help="máximo de casos extraídos (padrão: 30)")
    parser.add_argument("--scan-limit", type=int, default=300,
                        help="máximo de casos com trilhas compatíveis a inspecionar (padrão: 300)")
    args = parser.parse_args()

    if not args.clinical_csv.is_file():
        print(f"ERRO: não encontrei {args.clinical_csv}", file=sys.stderr)
        return 2

    df = pd.read_csv(args.clinical_csv)
    needed = {"caseid", "age", "sex", "weight", "height", "asa", "caseend"}
    missing = needed - set(df.columns)
    if missing:
        print(f"ERRO: colunas ausentes em clinical_data.csv: {sorted(missing)}", file=sys.stderr)
        return 2

    for colname in ("caseid", "age", "weight", "height", "asa", "caseend"):
        df[colname + "_num"] = numeric(df[colname])
    df["sex_norm"] = df["sex"].map(normalized_sex)
    df["age"] = df["age_num"]
    df["weight"] = df["weight_num"]
    df["height"] = df["height_num"]
    df["asa_num"] = df["asa_num"].round()
    df["caseid"] = df["caseid_num"]

    # Filtros explícitos de qualidade/demografia; valores censurados ou não numéricos
    # são excluídos, não convertidos em números aproximados.
    clinical_ok = (
        df["caseid"].notna() & df["age"].between(18, 89) &
        df["weight"].between(30, 250) & df["height"].between(120, 220) &
        df["asa_num"].between(1, 5) & df["sex_norm"].isin(["M", "F"]) &
        (df["caseend_num"] >= 1260)
    )
    candidates = df.loc[clinical_ok].copy().sort_values("caseid")
    candidates["caseid"] = candidates["caseid"].astype(int)
    print("Metadados clínicos: arquivo local do PhysioNet (v1.0.0, conforme download informado).")
    print(f"Versão dos sinais usada pela biblioteca/API: {API_DATASET_VERSION}")
    if API_DATASET_VERSION not in ("1.0.0", "desconhecida"):
        print("ATENÇÃO: a versão dos sinais da API difere da versão do CSV local. Registre essa diferença no relatório e não descreva os dados como uma única versão sem confirmar a compatibilidade.")
    print(f"Casos no CSV: {len(df)}")
    print(f"Casos com metadados clínicos elegíveis: {len(candidates)}")

    case_ids = candidates["caseid"].tolist()
    track_map = get_case_track_map(case_ids)
    candidates_with_tracks = []
    for _, row in candidates.iterrows():
        names = track_map.get(int(row["caseid"]), set())
        co_track = next((name for name in CO_TRACKS if name in names), None)
        if all(name in names for name in REQUIRED_TRACKS) and co_track:
            candidates_with_tracks.append(row)
    print(f"Casos com PA/FC/propofol e débito cardíaco medido: {len(candidates_with_tracks)}")

    OUT.mkdir(parents=True, exist_ok=True)
    accepted: list[dict] = []
    rejected: list[dict] = []
    scanned = 0
    for row in candidates_with_tracks:
        if scanned >= args.scan_limit or len(accepted) >= args.max_cases:
            break
        scanned += 1
        cid = int(row["caseid"])
        meta, reason = extract_one(row, track_map.get(cid, set()))
        if meta is not None:
            accepted.append(meta)
            print(f"[{len(accepted)}/{args.max_cases}] Caso {cid}: aceito; coadministração={meta['coadmin_status']}")
        else:
            rejected.append({"caseid": cid, "reason": reason})
            print(f"Caso {cid}: rejeitado ({reason})")

    pd.DataFrame(accepted).to_csv(OUT / "coorte.csv", index=False)
    pd.DataFrame(rejected, columns=["caseid", "reason"]).to_csv(OUT / "casos_rejeitados.csv", index=False)
    print("\nConcluído.")
    print(f"Casos aceitos: {len(accepted)}")
    print(f"Relatório da coorte: {OUT / 'coorte.csv'}")
    print(f"Casos rejeitados: {OUT / 'casos_rejeitados.csv'}")
    if not accepted:
        print("Nenhum caso passou todos os critérios. Consulte casos_rejeitados.csv; nenhum valor foi inventado.")
        return 1
    print("Observação: a API foi usada apenas para trilhas de casos selecionados; o ZIP completo não foi baixado.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
