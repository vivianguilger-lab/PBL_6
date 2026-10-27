# Guia do código

Este guia explica **como o código está organizado, o que acontece numa simulação e como reproduzir
tudo do zero**. Cada arquivo também tem comentários em cada bloco; este texto é o mapa para lê-los na
ordem certa. De onde vem cada número está em [validacao.md](validacao.md); as equações, em
[modelo_matematico.md](modelo_matematico.md).

## 1. Reproduzir do zero

1. **Instalar o Julia** (versão 1.10 ou mais nova; os testes automáticos usam a 1.11):
   - Windows, no PowerShell: `winget install --name Julia --id 9NJNWW8PVKMN -e -s msstore`
   - macOS/Linux, no terminal: `curl -fsSL https://install.julialang.org | sh`
2. **Baixar o projeto**: o [ZIP do ramo main](https://github.com/vivianguilger-lab/PBL_6/archive/refs/heads/main.zip)
   (extrair) ou `git clone https://github.com/vivianguilger-lab/PBL_6`.
3. **Instalar as dependências**, no terminal, dentro da pasta do projeto:

   ```
   julia --project=. -e "using Pkg; Pkg.instantiate()"
   ```

   A lista de pacotes e as versões aceitas estão em `Project.toml` (seção `[compat]`).
4. **Rodar os testes** (cerca de 1 minuto depois da compilação):

   ```
   julia --project=. -e "using Pkg; Pkg.test()"
   ```

   Todos devem passar. Se algum falhar, o próprio teste diz o que conferir (ver os comentários em
   `test/runtests.jl`).
5. **Regenerar todas as figuras e tabelas** de `resultados/`:

   ```
   julia --project=. scripts/rodar_tudo.jl
   ```

   Cada análise grava um log em `resultados/log_<script>.txt`. Os números devem sair iguais aos do
   README (o modelo é determinístico).
6. **Abrir a plataforma**: `julia --project=. -e "using PropofolWindkessel; wait(plataforma())"`, ou
   duplo clique em `abrir_plataforma.bat` (Windows) / `abrir_plataforma.command` (macOS/Linux).

## 2. Mapa dos arquivos

| arquivo | o que faz | origem |
|---|---|---|
| `src/PropofolWindkessel.jl` | o módulo: carrega os outros arquivos na ordem e lista as funções públicas | nosso |
| `src/windkessel.jl` | Windkessel de 2, 3 e 4 elementos, simulação, estado periódico, resposta em frequência, identificação de parâmetros | porta do windkessel-id feita pelo grupo |
| `src/dados.jl` | lê os dados de Stergiopulos 1999 e os dados in vivo (`dados/windkessel-id/`) | grupo + leitor nosso |
| `src/paciente.jl` | o `Paciente` (dados pré-indução), massa magra, superfície corporal, PAM e DC basais, leitura de TOML | grupo, com PAS/PAD nossos |
| `src/propofol.jl` | PK do propofol (Schnider, Marsh, Eleveld), bolus e infusão, BIS de Eleveld completo | grupo + Eleveld completo nosso |
| `src/acoplamento.jl` | efeito do propofol sobre Rp, VS e FC; fator de sensibilidade; simulação batimento a batimento | grupo |
| `src/calibracao.jl` | paciente virtual (Rp e complacência C do paciente) e a simulação da indução com linha de base | nosso |
| `src/metricas.jl` | PAM mínima, tempo e área abaixo de 65/55 mmHg, BIS, risco; relatório em texto | nosso |
| `src/dose.jl` | varredura de doses e janela de dose | nosso |
| `src/plataforma.jl` | servidor web local da plataforma (API em JSON) | nosso |
| `src/incerteza.jl`, `src/populacao.jl` | faixa de incerteza e coorte virtual — **experimentais, não validados, fora da plataforma** | nosso |
| `app/index.html` | a página da plataforma (formulário, gráficos em SVG, tabela de doses) | nosso |
| `scripts/01…06*.jl` | análises que geram `resultados/` (verificação, previsão, doses, sensibilidade, identificação) | nosso |
| `scripts/graficos.jl` | funções de figura usadas pelos scripts | nosso |
| `test/runtests.jl` | testes automáticos (verificação e plausibilidade) | grupo + nosso |
| `dados/*.toml` | pacientes de exemplo A, B e C | nosso |
| `docs/` | equações, validação e rastreabilidade, este guia | nosso |

## 3. O caminho de uma simulação

O que acontece quando a plataforma (ou o script 02) simula uma indução:

```
Paciente (idade, sexo, peso, altura, PAS, PAD, PAM, FC, ASA, PVI, IVCI)        paciente.jl
   │
   ├─► calibrar(p)                                                              calibracao.jl
   │     ├─ Rp = PAM / DC          (DC medido, ou 3,0 L/min/m² × sup. corporal)
   │     ├─ C ajustada até a pressão de pulso do modelo = PAS − PAD
   │     ├─ parâmetros PK do modelo escolhido (parametros_pk)                   propofol.jl
   │     └─ sensibilidade s (idade, ASA, PVI, IVCI)                             acoplamento.jl
   │
   └─► simular_inducao(pv; dose_mgkg, duracao_bolus, manutencao)                calibracao.jl
         ├─ esquema_inducao: bolus em t = 1 min + infusão de manutenção         propofol.jl
         ├─ simular_propofol:                                                   propofol.jl
         │     ├─ PK escolhido → Cp(t), Ce(t)           (NeoPKPD)   → move a hemodinâmica
         │     └─ PK e ke0 de Eleveld → Ce_bis(t) → sigmoide + atraso → BIS(t)
         ├─ simular_hemodinamica, a cada batimento:                             acoplamento.jl
         │     ├─ lê Ce no instante do batimento
         │     ├─ Rp, VS e FC = basal × (1 − Imax·s·Ceᵞ/(Ce50ᵞ + Ceᵞ))
         │     ├─ fluxo aórtico de Stergiopulos esticado para 60/FC e escalado para VS
         │     └─ Windkessel de 4 elementos → pressão → PAS, PAD, PAM do batimento  windkessel.jl
         └─ Resultado (curvas de PK/PD e batimentos)
   │
   └─► metricas(r): PAM mínima, queda, tempo < 65 e < 55, BIS, risco            metricas.jl
```

A plataforma faz o mesmo em `_api_simular` (`plataforma.jl`) e devolve tudo em JSON para a página
desenhar; `_api_doses` repete a simulação para 11 doses (0,5 a 3,0 mg/kg).

## 4. Unidades

| grandeza | unidade no código |
|---|---|
| tempo da PK/PD | minutos |
| tempo do Windkessel | segundos (passo h = 0,005 s) |
| dose | mg (bolus), mg/min (infusão); a plataforma recebe mg/kg e µg/kg/min |
| volumes PK | L; depurações em L/min; ke0 em 1/min |
| concentrações (Cp, Ce) | µg/ml (= mg/L) |
| pressão | mmHg |
| fluxo do Windkessel | ml/s; débito cardíaco mostrado em L/min |
| Rp, Rc | mmHg·s/ml; C em ml/mmHg; L em mmHg·s²/ml |
| FC | batimentos/min; VS em ml |

## 5. Onde mudar cada coisa

| quero mudar | onde | atenção |
|---|---|---|
| efeito do propofol na pressão (Ce50, γ, Imax) | `EfeitoHemodinamico` em `src/acoplamento.jl`, ou `calibrar(p; efeito = EfeitoHemodinamico(Imax_Rp = 0.2))` sem editar o arquivo | valores ilustrativos, a calibrar com dados |
| fator de sensibilidade | `sensibilidade_hemodinamica` em `src/acoplamento.jl` (coeficientes `k_*` como argumentos) | coeficientes ilustrativos |
| modelo PK | argumento `modelo_pk` de `calibrar` (`:schnider`, `:marsh`, `:eleveld`), ou a lista na plataforma | o BIS usa sempre o Eleveld |
| limiares de hipotensão e alvo de BIS | argumentos `limiares` e `alvo_bis` de `metricas` | fontes em validacao.md |
| protocolo (dose, tempo de injeção, manutenção) | argumentos de `simular_inducao`, ou os campos `[protocolo]` dos arquivos `dados/*.toml` | |
| paciente de exemplo | `dados/*.toml` | a PAM é obrigatória |

Depois de qualquer mudança no modelo: rode os testes, rode `scripts/rodar_tudo.jl` e registre a
alteração na seção 3 de [validacao.md](validacao.md).

## 6. Julia em 10 linhas

| escrita | significado |
|---|---|
| `f(x; a = 1)` | função com argumento `x` e argumento opcional nomeado `a` (padrão 1) |
| `x .* y`, `f.(v)` | o ponto aplica a operação a cada elemento do vetor |
| `(a = 1, b = 2)` | NamedTuple: um pacote de valores com nome, lido como `t.a` |
| `Base.@kwdef struct P … end` | estrutura cujos campos podem ser dados por nome e ter valor padrão |
| `cond ? a : b` | se `cond`, então `a`, senão `b` |
| `something(x, y)` | `x`, ou `y` se `x` for `nothing` (vazio) |
| `get(d, "k", padrão)` | valor da chave `"k"` do dicionário `d`, ou o padrão se não existir |
| `θ[1]` | primeiro elemento (Julia conta a partir de 1) |
| `@testset`, `@test a ≈ b rtol = 1e-6` | grupo de testes; "a é aproximadamente b, com erro relativo até 1e−6" |
| `using Pacote` | carrega um pacote; `--project=.` usa os pacotes listados neste projeto |

## 7. Glossário do modelo

| termo | significado |
|---|---|
| Cp, Ce | concentração de propofol no plasma e no sítio efetor (cérebro, vasos) |
| ke0 | velocidade com que o sítio efetor acompanha o plasma (1/min) |
| V1, V2, V3; CL, Q2, Q3 | volumes dos 3 compartimentos e depurações entre eles (modelo tricompartimental) |
| Ce50, γ, Imax | concentração de meio-efeito, inclinação da curva (Hill) e efeito máximo |
| BIS | índice bispectral (0–100); 40–60 = hipnose adequada |
| Rp | resistência periférica (vasodilatação baixa a Rp) |
| C | complacência arterial (aorta rígida = C baixa = pressão de pulso alta) |
| Rc, L | impedância característica da aorta e inertância do sangue |
| VS, FC, DC | volume sistólico, frequência cardíaca, débito cardíaco = VS × FC |
| PAM | pressão arterial média = Rp × DC (em regime) |
| s | sensibilidade do paciente ao efeito hemodinâmico (1 = referência) |
| paciente virtual | o modelo com Rp e C ajustados aos dados pré-indução de um paciente |
