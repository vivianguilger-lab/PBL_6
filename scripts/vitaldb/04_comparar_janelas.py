"""Analise exploratoria de sensibilidade temporal das previsoes VitalDB.
Nao modifica sinais, previsoes ou parametros do modelo.
Execute na raiz do projeto: python scripts/vitaldb/04_comparar_janelas.py
"""
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / 'resultados' / 'vitaldb' / 'predicoes'
OUTPUT = ROOT / 'resultados' / 'vitaldb' / 'analise_janelas'
WINDOWS = [("0_5min", 300), ("0_10min", 600), ("0_20min", 1200)]


def main():
    paths = sorted(INPUT.glob('case_*_padrao.csv'))
    if not paths:
        raise SystemExit(f'Nao encontrei previsoes padrao em: {INPUT}\nExecute primeiro o script Julia de validacao.')
    OUTPUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for path in paths:
        d = pd.read_csv(path)
        required = ['time_after_induction_s', 'pam_observed_mmHg', 'pam_predicted_mmHg']
        absent = [col for col in required if col not in d.columns]
        if absent:
            print(f'Ignorando {path.name}: faltam colunas {absent}')
            continue
        caseid = int(path.name.split('_')[1])
        t = pd.to_numeric(d[required[0]], errors='coerce')
        o = pd.to_numeric(d[required[1]], errors='coerce')
        p = pd.to_numeric(d[required[2]], errors='coerce')
        for label, stop in WINDOWS:
            included = t.notna() & np.isfinite(t) & (t >= 0) & (t <= stop)
            paired = included & np.isfinite(o) & np.isfinite(p)
            n = int(paired.sum())
            if n:
                err = (p[paired] - o[paired]).to_numpy()
                mae = float(np.mean(np.abs(err)))
                rmse = float(np.sqrt(np.mean(err ** 2)))
                bias = float(np.mean(err))
                tmax = float(t[paired].max())
                span_coverage = 100 * tmax / stop
            else:
                mae = rmse = bias = tmax = span_coverage = float('nan')

            rows.append(dict(caseid=caseid, janela=label, limite_s=stop,
                             n_pares=n, ultimo_tempo_valido_s=tmax,
                             cobertura_temporal_percentual=span_coverage,
                             mae_mmHg=mae, rmse_mmHg=rmse, vies_mmHg=bias))
            
    result = pd.DataFrame(rows)
    if result.empty:
        raise SystemExit('Nenhum arquivo com colunas validas encontrado.')
    result.to_csv(OUTPUT / 'metricas_por_janela.csv', index=False)
    summary = (result.groupby('janela', as_index=False)
               .agg(n_pacientes=('caseid', 'nunique'),
                    media_mae_mmHg=('mae_mmHg', 'mean'),
                    mediana_mae_mmHg=('mae_mmHg', 'median'),
                    media_rmse_mmHg=('rmse_mmHg', 'mean'),
                    mediana_rmse_mmHg=('rmse_mmHg', 'median'),
                    menor_cobertura_pct=('cobertura_temporal_percentual', 'min')))
    order = [name for name, _ in WINDOWS]
    summary['janela'] = pd.Categorical(summary['janela'], categories=order, ordered=True)
    summary = summary.sort_values('janela')
    summary.to_csv(OUTPUT / 'resumo_janelas.csv', index=False)
    print('\nResumo por paciente (sem misturar observacoes de pacientes diferentes):')
    print(summary.to_string(index=False, float_format=lambda x: f'{x:.2f}'))
    print(f'\nArquivos salvos em {OUTPUT}')
    try:
        import matplotlib.pyplot as plt
        pivot = result.pivot(index='caseid', columns='janela', values='rmse_mmHg').reindex(columns=order)
        ax = pivot.plot(kind='bar', figsize=(12, 5), width=.8)
        ax.set(ylabel='RMSE da PAM (mmHg)', xlabel='Caso VitalDB', title='Sensibilidade do erro a janela temporal')
        ax.legend(title='Janela')
        plt.tight_layout()
        plt.savefig(OUTPUT / 'rmse_por_caso_e_janela.png', dpi=180)
        plt.close()
        print('Grafico salvo: rmse_por_caso_e_janela.png')

    except ImportError:
            print('Matplotlib nao instalado; tabelas CSV geradas normalmente.')

if __name__ == '__main__':
    main()

