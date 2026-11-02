# Reconstrução da validação VitalDB — PBL_6

Este patch é **aditivo**: contém somente novos scripts e documentação. Não substitui `src/`, `Project.toml`, `.gitignore` ou `scripts/07_validacao_vitaldb.jl`.

## O que a implementação corrige

- Usa os metadados reais de `dados/vitaldb/metadados/clinical_data.csv`.
- Consulta o índice público de trilhas para identificar casos que realmente tenham pressão arterial invasiva, FC, taxa de propofol e débito cardíaco medido.
- Carrega apenas as trilhas dos casos selecionados pela API VitalDB; não baixa o ZIP de 92 GB. O script registra a versão de dataset que a biblioteca usa para os sinais, pois ela pode diferir da versão 1.0.0 dos CSVs locais.
- Não substitui PAS/PAD/PAM/FC/DC ausentes por números padrão.
- Exclui casos sem DC medido, sem baseline de 60 s, sem 20 min de seguimento, com lacunas longas na taxa de propofol ou com baseline incompatível.
- Converte `Orchestra/PPF20_RATE` de mL/h para mg/min (mL/h × 20 mg/mL ÷ 60).
- Preserva lacunas nos sinais fisiológicos; não interpola os CSVs de observação. A comparação permite interpolação linear apenas entre observações separadas por no máximo 5 s.
- Ajusta os parâmetros PD `Ce50`, `gamma`, `Imax_Rp`, `Imax_VS` e `Imax_FC` somente no conjunto de treino e avalia em pacientes separados para teste.
- Mantém `Rc`, `L`, o modelo PK de Schnider e a forma de fluxo de Stergiopulos fixos. `C` é calibrado por paciente usando a pressão de pulso basal e `Rp` usa PAM e DC medidos.
- Define explicitamente `s=1`, porque PVI/IVCI não estão garantidos nos metadados; não trata os valores padrão desses campos como dados reais.

## Pré-requisitos

- O clone Git do `PBL_6` aberto no VS Code.
- Os cinco CSVs do PhysioNet em `dados/vitaldb/metadados/`.
- Python com acesso à internet para a API VitalDB.
- Julia 1.10 e dependências já declaradas em `Project.toml`.

## Execução no Windows / PowerShell

Na raiz do `PBL_6`, crie um ambiente Python separado:

```powershell
py -m venv .venv-vitaldb
.\.venv-vitaldb\Scripts\Activate.ps1
python -m pip install -r scripts\vitaldb\requirements.txt
```

Copie os arquivos deste patch para o clone preservando a estrutura; **não extraia por cima de arquivos existentes**. Depois, confirme:

```powershell
git status --short
```

Para começar com uma coorte pequena:

```powershell
python scripts\vitaldb\01_extrair_dados.py --max-cases 10 --scan-limit 150
```

Se houver poucos casos elegíveis, aumente gradualmente:

```powershell
python scripts\vitaldb\01_extrair_dados.py --max-cases 30 --scan-limit 500
```

Depois, rode a validação em Julia:

```powershell
julia --project=. scripts\vitaldb\02_validar_calibrar_vitaldb.jl
```

## Saídas

Dados processados (locais; não versionar): `dados/vitaldb/processados/`

Resultados: `resultados/vitaldb/`

- `metricas_por_caso.csv`: métricas por paciente e método.
- `resumo_teste.csv`: médias e medianas das métricas nos casos de teste, separando todos os testes de `none_recorded`, `active` e `unknown` para explicitar possível coadministração.
- `parametros_calibrados.csv`: parâmetros padrão e ajustados, status da otimização.
- `predicoes/`: CSVs com observações, previsões e resíduos por batimento.
- `figuras/`: PAM observada vs prevista e dispersão no teste.
- `casos_rejeitados.csv` e `casos_excluidos_julia.csv`: auditoria dos motivos de exclusão.

Acrescente ao `.gitignore` existente (sem apagar regras):

```gitignore
/dados/vitaldb/
/resultados/vitaldb/
```

## Interpretação científica e limites

0. **Versão dos dados:** o CSV clínico local é do PhysioNet v1.0.0. A biblioteca `vitaldb` pode apontar para uma versão de sinais mais recente (por exemplo, 1.0.1). O script imprime e grava essa versão. Para um estudo formal, confirme a compatibilidade das versões antes de juntar metadados e sinais; não esconda essa diferença no relatório.

1. O split é determinístico por ordem de `caseid`: cada quinto caso vai para teste. O teste não entra no ajuste.
2. A calibração é omitida se não houver pelo menos 3 casos de treino e 2 de teste com status `none_recorded` para outras taxas de infusão. `unknown` nunca é tratado como ausência de coadministração.
3. `none_recorded` significa apenas que as outras taxas de infusão disponíveis no registro não indicaram administração concomitante; não prova que nenhum outro medicamento foi administrado por outra via.
4. O modelo continua tendo hipóteses estruturais: Windkessel, forma de fluxo aórtico, PK de Schnider, `Rc`/`L` fixos e efeito PD sigmoide. A validação pode testar essas hipóteses, não eliminá-las.
5. `Ce50` e `gamma` podem não ser identificáveis separadamente dos `Imax` com uma amostra pequena. Se os parâmetros atingirem os limites, a métrica do teste piorar ou houver grande variação entre divisões, isso deve ser reportado; não se deve escolher parâmetros apenas porque deixam a curva visualmente bonita.
6. Este é um pipeline de pesquisa, não uma ferramenta para orientar cuidado clínico.
