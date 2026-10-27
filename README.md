# Modelo *in silico* para predição hemodinâmica da indução anestésica com propofol

**PBL 6 — Modelagem e Simulação de Sistemas Biomédicos**

> **Pergunta:** *Como prever, por simulação, a resposta hemodinâmica de um paciente ao propofol?*

**Proposta:** um **paciente virtual** montado com os dados pré-indução do próprio paciente, em
que a PK/PD do propofol (**NeoPKPD**: Schnider/Marsh + BIS de Eleveld) move os parâmetros de um
**Windkessel de 4 elementos** (**windkessel-id**, Stergiopulos 1999). O modelo prevê a curva de
PAM para a dose escolhida **antes da injeção**.

O núcleo do modelo segue o repositório do grupo,
**[PropofolPBL](https://github.com/gabiadamek30-bit/Simula-o-da-resposta-do-paciente-ao-propofol---PBL)**:
mesmas equações, mesmos valores padrão e os mesmos testes. Sobre ele, este repositório
acrescenta calibração com PAS/PAD, métricas de hipotensão, janela de dose, faixa de incerteza,
coorte virtual e uma plataforma web para médicos.

> **Para usar:** leia o **[MANUAL.md](MANUAL.md)** (passo a passo). Resumo: baixe o
> [ZIP](https://github.com/vivianguilger-lab/PBL_6/archive/refs/heads/main.zip),
> extraia e dê duplo clique em `abrir_plataforma.bat` (Windows) ou `abrir_plataforma.command`
> (macOS/Linux) — o Julia é instalado sozinho se faltar.

![Plataforma PBL-6](docs/figuras/plataforma.png)

## O problema

- O propofol causa **vasodilatação periférica**, **depressão do miocárdio** e **atenuação do
  barorreflexo** → hipotensão.
- **Hipotensão pós-indução** = queda de pressão nos primeiros 20 min. Mesmo curtos períodos de
  **PAM < 55 mmHg** associam-se a lesão renal aguda e lesão miocárdica (Walsh 2013), e PAM < 65 mmHg
  a lesão de órgãos (Salmasi 2017).
- A dose segue uma faixa por peso ajustada por idade, sexo, estado de saúde e medicações, mas a
  resposta varia: **a dose certa só se confirma depois da injeção.**

## Como o modelo funciona

```
Paciente ──► PK, ke0 e atraso de Eleveld ──────────────────────────────► BIS (Eleveld 2018)
        └──► PK escolhido (Schnider/Marsh/Eleveld, NeoPKPD) ──► Cp(t) ──► ke0 ──► Ce(t)
                                                                                 │
  PAM, FC, DC basais ──► θ = [Rp, C, Rc, L] ◄────────────────────────────────────┘ Emax sigmoide × sensibilidade
  PAS − PAD ──► C do paciente virtual  │                                          (idade, ASA, PVI, IVCI)
                                        ▼
          fluxo aórtico (Stergiopulos) ─► Windkessel de 4 elementos ─► pressão arterial batimento a batimento
```

| bloco | modelo | origem |
|---|---|---|
| farmacocinética | tricompartimental (Schnider, Marsh ou Eleveld) com bolus + manutenção e compartimento de efeito | [NeoPKPD](https://github.com/shramish2057/NeoPKPD), como no grupo |
| hipnose | BIS de Eleveld (2018) completo: PK, ke0, sigmoide com Ce50 dependente da idade e atraso | grupo + Eleveld completo (validação) |
| circulação | Windkessel de 4 elementos [Rp, C, Rc, L], alimentado pelo fluxo aórtico de Stergiopulos (1999) escalado por VS e FC | [windkessel-id](https://github.com/hpigot/windkessel-id), porta do grupo |
| efeito do propofol | Emax sigmoide sobre Rp (Imax 0,25), VS (0,15) e FC (0,05), Ce50 = 3 µg/ml, γ = 2, × sensibilidade *s* | grupo (valores ilustrativos) |
| sensibilidade | *s* = 1 + 0,01(idade−50)⁺ + 0,15(ASA−2)⁺ + 0,03(PVI−15)⁺ + 0,01(IVCI−43)⁺ | grupo; limiares de PVI e IVCI da literatura |
| paciente virtual | Rp = PAM/DC; **C identificado** para reproduzir a pressão de pulso medida | extensão |

Equações completas: **[docs/modelo_matematico.md](docs/modelo_matematico.md)** ·
**Validação e rastreabilidade** (arquivo e linha de cada valor, fonte, veredito, pendências e
bibliografia): **[docs/validacao.md](docs/validacao.md)** ·
Revisão inicial de repositórios: [docs/revisao_github_literatura.md](docs/revisao_github_literatura.md)

## Instalação e uso

Requer Julia ≥ 1.10 (instalador oficial **juliaup**:
Windows `winget install --name Julia --id 9NJNWW8PVKMN -e -s msstore`;
Linux/macOS `curl -fsSL https://install.julialang.org | sh`).

**Plataforma para médicos (navegador)** — no Julia:

```julia
using Pkg
Pkg.add(url="https://github.com/vivianguilger-lab/PBL_6")
using PropofolWindkessel
plataforma()          # abre http://localhost:8080
```

A primeira instalação demora (o NeoPKPD traz o DifferentialEquations). Quem já tinha instalado:
`using Pkg; Pkg.update()`. Com a pasta do projeto: duplo clique em `abrir_plataforma.bat`
(Windows) ou `abrir_plataforma.command` (macOS/Linux), que instalam o Julia se faltar. Os dados não saem do computador (o servidor só aceita conexões locais).

**No REPL:**

```julia
using PropofolWindkessel
p  = Paciente(idade = 62, sexo = :F, peso = 72, altura = 160, pas = 148, pad = 82, pam = 104, fc = 78,
              asa = 2, pvi = 12, ivci = 25)
pv = calibrar(p)                                    # paciente virtual
r  = simular_inducao(pv; dose_mgkg = 1.5)           # bolus em 30 s + 100 µg/kg/min
relatorio(metricas(r))
janela_de_dose(pv)

# interface do grupo, inalterada:
pkpd, batimentos, onda = simular_paciente(p; tfim = 30.0, forma_de_onda = true)
```

**Scripts** (`julia --project=. -t auto scripts/<script>`; `rodar_tudo.jl` regenera `resultados/`):

| script | o que faz |
|---|---|
| `01_verificacao.jl` | Windkessel e NeoPKPD contra soluções exatas |
| `02_prever_paciente.jl dados/paciente_exemplo.toml` | paciente virtual, previsão, ondas de pressão e fluxo |
| `03_janela_de_dose.jl dados/paciente_idoso.toml` | varredura de doses, janela individual, velocidade de injeção |
| `04_coorte_virtual.jl 500` | ensaio clínico virtual (**experimental, não validado**; fora de `rodar_tudo.jl`) |
| `05_sensibilidade.jl dados/paciente_exemplo.toml` | sensibilidade local (±10%) dos parâmetros do grupo |
| `06_identificacao_windkessel.jl` | identificação do Windkessel 2/3/4 elementos (Stergiopulos e in vivo) |

Testes: `julia --project=. -e "using Pkg; Pkg.test()"` (163 testes, inclusive os do grupo, a conferência do
Eleveld contra o PAS e a dos links de [docs/validacao.md](docs/validacao.md)).

## Resultados (com os parâmetros ilustrativos do grupo)

**Verificação** — Windkessel de 2 elementos × `P₀·e^{−t/RpC}`: erro 5 × 10⁻¹⁴; NeoPKPD ×
`exp(K·t)·A₀`: 8 × 10⁻¹⁰; Windkessel de 4 elementos: PAM = Rp·Q̄ com erro 8 × 10⁻¹⁶.
A identificação recupera θ de Stergiopulos (Rp 0,80; C 1,25; Rc 0,056; L 0,0055 contra
0,79; 1,22; 0,056; 0,0051 publicados; RMSE 2,4 mmHg).

![Verificação](resultados/verificacao.png)

**Três pacientes** (bolus em 30 s + manutenção de 100 µg/kg/min; risco pela PAM mínima frente a 65 e 55 mmHg):

| paciente | dose | PAM basal → mínima | risco | BIS mínimo (BIS < 60 em) | janela de dose |
|---|---|---|---|---|---|
| A — 62 a, ♀, ASA II, PVI 12% | 1,5 mg/kg | 104 → 66 mmHg (2,0 min) | baixo | 45 (2,2 min) | 0,50–1,63 mg/kg |
| B — 28 a, ♂, ASA I | 2,0 mg/kg | 89 → 57 mmHg (1,8 min) | moderado | 46 (1,7 min) | não existe: BIS < 60 só a partir de ~1,4 mg/kg; PAM ≥ 65 só até ~0,9 mg/kg |
| C — 81 a, ♂, ASA III, PVI 18%, IVCI 45% (s = 1,57) | 1,5 mg/kg | 107 → 54 mmHg (2,3 min) | alto | 39 (2,7 min) | 0,50–0,83 mg/kg |

O limite inferior de 0,50 mg/kg é o menor valor testado: com a manutenção de 100 µg/kg/min, até
0,5 mg/kg chega a BIS < 60 em alguns minutos (critério em revisão, [docs/validacao.md](docs/validacao.md), seção 6).

![Paciente C](resultados/paciente_idoso_inducao.png)

A faixa de incerteza e a coorte virtual foram retiradas dos resultados até que as variabilidades
sejam validadas ([docs/validacao.md](docs/validacao.md), seção 6).

**Sensibilidade:** o que mais muda a PAM mínima é o fator de sensibilidade *s*, o Imax da
vasodilatação, o Ce50 hemodinâmico e a dose. A complacência C quase não influi na PAM (PAM = Rp·Q̄).

## Limitações

- Os parâmetros do efeito hemodinâmico e da sensibilidade são **ilustrativos** (como diz o
  grupo) e precisam ser calibrados com dados reais (proposta com o banco aberto VitalDB em
  [docs/validacao.md](docs/validacao.md)). **Ainda não há validação com dados de pacientes.**
- Modelo aberto: o fluxo aórtico é imposto; não há barorreflexo nem retorno venoso explícitos.
- A hemodinâmica usa o Ce do PK escolhido; o efeito na pressão na prática é mais lento (Kazama 1999).
- PAM e pressão de pulso são reproduzidas exatamente; PAS e PAD individualmente, não (a forma
  da onda de Stergiopulos é fixa).
- Sem opioides, laringoscopia, fluidos ou vasopressores.

## Estrutura

```
src/
  PropofolWindkessel.jl  módulo
  windkessel.jl          porta do windkessel-id (grupo): G2/G3/G, solvex0, identify
  dados.jl               dados de Stergiopulos e in vivo (dados/windkessel-id/)
  paciente.jl            biomarcadores (PAS, PAD, PAM, FC, ASA, PVI, IVCI…)
  propofol.jl            PK Schnider/Marsh/Eleveld + BIS de Eleveld completo sobre o NeoPKPD
  acoplamento.jl         Ce → Rp, VS, FC → Windkessel → PA batimento a batimento (grupo)
  calibracao.jl          paciente virtual e simulação da indução
  metricas.jl, dose.jl, plataforma.jl
  incerteza.jl, populacao.jl   experimentais (não validados, fora da plataforma)
app/index.html           interface da plataforma
scripts/                 análises 01–06 e figuras
dados/                   pacientes de exemplo (TOML) e dados do windkessel-id
resultados/              figuras, CSV e logs
docs/                    modelo, validação e rastreabilidade (V&V 40, bibliografia), revisão inicial
test/                    testes (GitHub Actions)
```

O modelo anterior deste repositório (coração de elastância variável em malha fechada com
barorreflexo) está no histórico do git, até o commit `d3dd331`.
