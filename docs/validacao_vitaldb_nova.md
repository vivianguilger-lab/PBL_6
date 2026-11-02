# Protocolo de validação externa do PBL_6 com VitalDB

## Objetivo

Comparar a previsão hemodinâmica do modelo com PAM, FC e débito cardíaco observados em registros intraoperatórios reais, e avaliar se a calibração de parâmetros melhora a previsão em pacientes não usados no ajuste.

## Etapa 1 — Seleção e aquisição

`01_extrair_dados.py` usa `clinical_data.csv` local para os dados clínicos. A consulta do índice de trilhas e a leitura das séries temporais são feitas pela biblioteca `vitaldb` na API pública. O ZIP completo e os arquivos `.vital` locais não são necessários. A versão dos sinais usada pela API é registrada: o CSV local foi baixado como v1.0.0, enquanto a biblioteca pode consultar uma versão mais recente. Para publicação ou análise formal, é necessário confirmar a compatibilidade de versões ou alinhar metadados e sinais à mesma versão.

Critérios mínimos: idade numérica entre 18 e 89 anos, sexo M/F, peso e altura disponíveis, ASA entre 1 e 5, pressão arterial invasiva (PAS/PAD/PAM), FC, taxa `Orchestra/PPF20_RATE` e um track direto de CO (`Vigileo/CO`, `EV1000/CO`, `Vigilance/CO` ou `CardioQ/CO`). Valores censurados ou ausentes não são convertidos em estimativas.

## Etapa 2 — Alinhamento e qualidade

- A indução é definida pelo primeiro início sustentado de taxa de propofol positiva, precedido por 60 s de taxa observada igual a zero.
- A janela contém 60 s de baseline e 1.200 s após o início.
- A taxa da bomba é mantida constante por no máximo 5 s entre registros; lacunas maiores levam à rejeição do caso.
- A taxa PPF20 em mL/h é convertida para mg/min pela concentração de 20 mg/mL.
- Valores fisiológicos fora de faixas de plausibilidade amplas são marcados como ausentes, nunca substituídos por valores inventados.
- Sinais de pressão/FC/CO não são imputados no arquivo de observações. Na comparação com a previsão, só se interpola entre observações válidas separadas por até 5 s.
- O baseline deve conter observações suficientes e satisfazer PAS > PAM > PAD; o DC basal deve ser medido, não estimado por superfície corporal.

## Etapa 3 — Simulação com entrada real

O vetor de infusão real é convertido em segmentos de taxa constante e fornecido ao modelo PK de Schnider. O paciente é construído com idade, sexo, peso, altura e ASA do CSV e com PAS/PAD/PAM/FC/DC calculados como medianas observadas no baseline. Não são usados 120/80 mmHg, 75 bpm ou DC estimado como fallback.

O parâmetro de sensibilidade individual `s` é fixado explicitamente em 1, pois PVI/IVCI não são garantidos. `C` é ajustado ao pulso basal observado e `Rp` decorre da PAM e DC medidos. A forma de fluxo aórtico, `Rc` e `L` continuam sendo hipóteses estruturais do modelo.

## Etapa 4 — Métricas

Para cada paciente, PAM, FC e DC previstos são comparados às observações no período de 20 min após o início da infusão. São calculados MAE, RMSE, viés médio e correlação quando matematicamente definida. A unidade de cada erro é mantida explícita.

## Etapa 5 — Calibração e teste independente

A cada quinto caso da coorte ordenada é atribuído ao teste; os demais ficam no treino. Os parâmetros `Ce50`, `gamma`, `Imax_Rp`, `Imax_VS` e `Imax_FC` são ajustados somente nos casos de treino com status `none_recorded` para outras taxas de infusão. A função objetivo é a média, por caso, dos erros quadráticos relativos de PAM, FC e DC, dando peso semelhante a cada sinal. A calibração só é executada com pelo menos 3 casos de treino limpos e 2 de teste limpos.

O conjunto de teste é usado exclusivamente para avaliação. Reportar tanto os parâmetros padrão quanto os calibrados. Um melhor ajuste no treino não é evidência de generalização; a conclusão depende dos resultados do teste, da quantidade de pacientes, da cobertura dos sinais e da estabilidade dos parâmetros.

## Referências do dataset

- Lee H, Jung CW. *VitalDB, a high-fidelity multi-parameter vital signs database in surgical patients*. Scientific Data. 2022;9:279. https://doi.org/10.1038/s41597-022-01411-5
- Documentação do VitalDB e API Python: https://vitaldb.net/docs/?documentId=API%2FPython_Library.md
- Descrição das trilhas e unidades: https://vitaldb.net/dataset/
