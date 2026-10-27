# Roteiro de validação: ordem de revisão do código

Ordem para revisar o modelo **arquivo por arquivo, linha por linha**, seguindo o caminho dos dados:
primeiro o que entra (dados e paciente), depois a farmacologia, a circulação, o acoplamento entre as
duas e, por fim, o que sai (métricas, dose, plataforma). Cada etapa depende da anterior; vale fechar
uma antes de passar à próxima.

Para cada arquivo: **o que faz**, **as linhas a conferir** (com link), **de onde vem cada dado**,
**vieses possíveis** e **como validar**. Os comentários dentro de cada arquivo explicam o código bloco
a bloco; as fontes e vereditos de cada número estão em [validacao.md](validacao.md).

> "Viés" aqui = erro sistemático: algo que empurra a previsão sempre para o mesmo lado, ou que vale
> para um grupo de pacientes e não para outro. Os itens marcados **(hipótese)** são suspeitas
> razoáveis que ainda precisam ser medidas; não são defeitos confirmados.

| etapa | arquivo | pergunta da etapa |
|---|---|---|
| 1 | `dados/` e `src/dados.jl` | os dados de entrada são os publicados? |
| 2 | `src/paciente.jl` | o paciente é descrito corretamente? |
| 3 | `src/propofol.jl` | a concentração de propofol e o BIS estão certos? |
| 4 | `src/windkessel.jl` | a circulação basal está certa? |
| 5 | `src/acoplamento.jl` | o efeito do propofol na circulação está certo? |
| 6 | `src/calibracao.jl` | o paciente virtual reproduz o paciente real? |
| 7 | `src/metricas.jl` | o resumo e o risco estão certos? |
| 8 | `src/dose.jl` | a janela de dose faz sentido clínico? |
| 9 | `src/plataforma.jl` e `app/index.html` | a tela mostra o que o modelo calcula? |
| 10 | `test/runtests.jl` | os testes cobrem o que importa? |
| 11 | `scripts/` | os resultados do README se reproduzem? |
| — | `src/incerteza.jl`, `src/populacao.jl` | **pular por ora**: experimentais, fora da plataforma |

---

## Etapa 1 — Dados de entrada: `dados/` e `src/dados.jl`

**O que faz.** Guarda e lê a forma de onda de fluxo e pressão aórtica de Stergiopulos 1999 (que alimenta
todo o Windkessel), os dados in vivo de porcos (só para o exemplo de identificação) e os três pacientes de
exemplo (`dados/*.toml`).

| linhas | o que é | de onde vem |
|---|---|---|
| [dados.jl#L41](../src/dados.jl#L41) | lê fluxo (ml/s) e pressão (mmHg), 166 amostras a 5 ms | cópia byte a byte de `stergiopulos1999_data/` do windkessel-id v1.0, digitalizada da Fig. 4A do artigo |
| [dados.jl#L70](../src/dados.jl#L70) | lê os dados in vivo (pressão em mmHg, fluxo em L/min) | windkessel-id; porcos de 65 kg |
| [paciente_exemplo.toml](../dados/paciente_exemplo.toml), [paciente_jovem.toml](../dados/paciente_jovem.toml), [paciente_idoso.toml](../dados/paciente_idoso.toml) | pacientes A, B e C | **fictícios**, criados por nós para ilustrar |

**Vieses possíveis.**
- A forma de onda vem de **um único batimento humano digitalizado de uma figura** (WebPlotDigitizer): há erro
  de digitalização, e a mesma forma é usada para todos os pacientes, de qualquer idade.
- A forma de onda de pressão muda com a idade (reflexão de onda; tipos A, B e C de Murgo). Usar um batimento
  "tipo A" para todos pode distorcer PAS e PAD em idosos **(hipótese)**; a PAM não é afetada, porque depende
  só do fluxo médio.
- Os pacientes de exemplo são inventados: servem para demonstrar, não para validar.

**Como validar.** Abrir a Fig. 4A e a Tabela 1 de Stergiopulos 1999 e comparar com `aflow.csv`/`apressure.csv`
(o script `scripts/06_identificacao_windkessel.jl` refaz o ajuste e deve recuperar os parâmetros publicados).

---

## Etapa 2 — Paciente: `src/paciente.jl`

**O que faz.** Define os dados pré-indução, calcula IMC, massa magra (James), superfície corporal
(Mosteller), PAM basal (a do monitor) e DC basal (medido ou estimado).

| linhas | o que é | de onde vem |
|---|---|---|
| [paciente.jl#L50](../src/paciente.jl#L50) | campos, unidades e valores padrão | estrutura do PropofolPBL do grupo; PAS/PAD nossos |
| [paciente.jl#L93](../src/paciente.jl#L93) | checagens (PAS > PAD, PAD < PAM < PAS, ASA 1–5) | nossas |
| [paciente.jl#L121](../src/paciente.jl#L121) | massa magra de James | James 1976, usada por Schnider 1998 |
| [paciente.jl#L133](../src/paciente.jl#L133) | superfície corporal | Mosteller 1987 |
| [paciente.jl#L140](../src/paciente.jl#L140) | PAM do monitor | medida |
| [paciente.jl#L154](../src/paciente.jl#L154) | DC = medido, ou 3,0 L/min/m² × superfície | escolha do grupo (índice cardíaco normal 2,5–4,0) |

**Vieses possíveis.**
- **Massa magra de James em obesos:** a fórmula passa a *diminuir* quando o peso aumenta muito (IMC acima de
  cerca de 35–40), e isso entra direto na depuração (CL) de Schnider. Em obesos graves o PK de Schnider fica
  enviesado; por isso eles estão fora do contexto de uso.
- **Manguito × aorta:** PAS e PAD do manguito no braço não são a pressão aórtica. A PAS no braço costuma ser
  maior que na aorta (amplificação do pulso), principalmente em jovens. Como a complacência C é calibrada com
  a pressão de pulso do manguito, C tende a sair menor que a real **(hipótese)**. A PAM é parecida no braço e
  na aorta, então a curva de PAM é pouco afetada.
- **PAM oscilométrica:** tem erro de medida próprio do aparelho; o modelo assume que a PAM digitada é a verdadeira.
- **Sexo binário (M/F):** imposto pelas fórmulas de massa magra e de Eleveld.
- **DC estimado:** 3,0 L/min/m² é um valor típico; não muda a queda relativa da PAM (teste automático), só a
  forma da onda.

**Como validar.** Conferir as fórmulas com James 1976 (via método de Schnider 1998) e Mosteller 1987; para o
viés do manguito, comparar a C calibrada com o manguito × com a linha arterial no mesmo paciente (VitalDB).

---

## Etapa 3 — Farmacocinética e BIS: `src/propofol.jl`

**O que faz.** Calcula os parâmetros PK (Schnider, Marsh ou Eleveld), monta bolus e infusão, simula a
concentração no plasma (Cp) e no sítio efetor (Ce) com o NeoPKPD e calcula o BIS pelo modelo de Eleveld completo.

| linhas | o que é | de onde vem |
|---|---|---|
| [propofol.jl#L66-L72](../src/propofol.jl#L66-L72) | PK de Schnider | Schnider 1998 (PK) e 1999 (ke0) |
| [propofol.jl#L75](../src/propofol.jl#L75) | erro fora da faixa das equações | nosso (substituiu pisos sem fonte) |
| [propofol.jl#L84-L90](../src/propofol.jl#L84-L90) | PK de Marsh | Marsh 1991; ke0 0,26 do Diprifusor (Absalom 2009) |
| [propofol.jl#L141](../src/propofol.jl#L141) | bolus em 30 s + manutenção | protocolo do grupo; bula do Diprivan |
| [propofol.jl#L175](../src/propofol.jl#L175) | PK de Eleveld (18 θ, covariáveis) | Eleveld 2018, na forma do PAS (testado em 20 casos) |
| [propofol.jl#L275](../src/propofol.jl#L275) | BIS: E0 93, Ce50 por idade, γ 1,89/1,47 | Eleveld 2018, Tabela 3 |
| [propofol.jl#L296](../src/propofol.jl#L296) | atraso do BIS | Eleveld 2018, Tabela 3 |
| [propofol.jl#L307](../src/propofol.jl#L307) | Cp e Ce pelo NeoPKPD | NeoPKPD (verificado contra a solução exata) |
| [propofol.jl#L353](../src/propofol.jl#L353) | junta tudo: Cp/Ce do PK escolhido + BIS de Eleveld | nosso |

**Vieses possíveis.**
- **Schnider foi estimado em voluntários saudáveis** (número e faixa etária a conferir no artigo), não em
  pacientes cirúrgicos; idosos frágeis, ASA III–IV e obesos estão pouco representados.
- **V1 fixo de 4,27 L no Schnider:** superestima o pico de concentração logo após o bolus (limitação já
  anotada pelo grupo). Isso tende a antecipar e aprofundar o efeito previsto nos primeiros minutos **(hipótese)**.
- **Marsh escala tudo pelo peso total:** em obesos superestima os volumes e a dose.
- **ke0 0,26 do Marsh** não tem publicação em revista.
- **O modelo prevê o paciente médio:** a variabilidade entre pessoas é grande (o próprio Eleveld traz a
  variabilidade do Ce50) e está desativada até ser validada.
- **Eleveld com amostra arterial e sem opioide:** escolha nossa, coerente com o contexto de uso; com opioide o
  modelo muda.

**Como validar.** Abrir as tabelas de Schnider 1998 e de Eleveld 2018 (pendência 6.1 de validacao.md); no
VitalDB, comparar o BIS previsto com o medido nas induções sem remifentanil.

---

## Etapa 4 — Circulação basal: `src/windkessel.jl`

**O que faz.** O Windkessel de 4 elementos (e os de 2 e 3), discretizado a cada 5 ms; simula a pressão a partir
do fluxo, acha o batimento periódico e identifica parâmetros a partir de dados.

| linhas | o que é | de onde vem |
|---|---|---|
| [windkessel.jl#L104](../src/windkessel.jl#L104) | equações do Windkessel de 4 elementos | windkessel-id (Pigot et al.), porta do grupo |
| [windkessel.jl#L144](../src/windkessel.jl#L144) | simulação passo a passo | windkessel-id |
| [windkessel.jl#L206](../src/windkessel.jl#L206) | estado periódico (regime) | windkessel-id v2.0 (correção da v1.0) |
| [windkessel.jl#L321](../src/windkessel.jl#L321) | identificação por otimização | windkessel-id |
| [windkessel.jl#L383](../src/windkessel.jl#L383) | Rp 0,79; C 1,22; Rc 0,056; L 0,0051 | Stergiopulos 1999, Tabela 1 |

**Vieses possíveis.**
- **Modelo aberto:** o fluxo é imposto; a pressão não age de volta sobre o coração. Não há barorreflexo, que na
  vida real compensa parte da queda. Sem essa compensação, o modelo pode **superestimar a queda e a duração da
  hipotensão (hipótese)** — o propofol atenua o barorreflexo, mas não o anula.
- **Parâmetros de um tipo de batimento humano:** Rc e L ficam fixos para todos os pacientes.
- **Linear e concentrado (0-D):** não representa propagação de onda nem diferenças entre artérias.

**Como validar.** Os testes já verificam as propriedades exatas (PAM = Rp·Q̄, impedância em 0 e ∞, identificação
de dados sintéticos); a validação fisiológica precisa de pressão invasiva + fluxo/DC medidos (VitalDB, casos com
monitor de débito).

---

## Etapa 5 — Efeito do propofol na circulação: `src/acoplamento.jl`

**O que faz.** A cada batimento, lê o Ce e reduz Rp, VS e FC por uma curva de Hill multiplicada pela
sensibilidade do paciente; gera o fluxo do batimento e roda o Windkessel. **É o coração da previsão de PAM.**

| linhas | o que é | de onde vem |
|---|---|---|
| [acoplamento.jl#L48](../src/acoplamento.jl#L48) | Ce50 3, γ 2, Imax Rp 0,25 / VS 0,15 / FC 0,05 | **ilustrativos** (grupo), sem fonte |
| [acoplamento.jl#L80](../src/acoplamento.jl#L80) | fator 1 − Imax·s·Hill(Ce), piso 5% | grupo |
| [acoplamento.jl#L112](../src/acoplamento.jl#L112) | s pela idade, ASA, PVI, IVCI | limiares: Tsuchiya 2010, Zhang 2016; coeficientes **ilustrativos** |
| [acoplamento.jl#L134](../src/acoplamento.jl#L134) | Rp = PAM/DC | física do Windkessel |
| [acoplamento.jl#L163](../src/acoplamento.jl#L163) e [acoplamento.jl#L186](../src/acoplamento.jl#L186) | fluxo de Stergiopulos esticado para 60/FC e escalado para VS | grupo |
| [acoplamento.jl#L278](../src/acoplamento.jl#L278) | laço batimento a batimento | grupo |

**Vieses possíveis — os mais importantes do modelo.**
- **Ce50, γ e Imax ilustrativos:** toda a profundidade da queda depende deles. Sem calibração, a curva de PAM
  pode estar sistematicamente funda ou rasa demais.
- **Um único Ce50 para Rp, VS e FC:** Su 2022 encontrou EC50 de ~3 µg/ml para a resistência periférica, mas
  ~0,3 µg/ml para o volume sistólico. Com um Ce50 só, o efeito no VS aparece tarde e fraco demais **(hipótese)**.
- **O efeito usa o Ce do PK escolhido (ke0 de Schnider = 0,456/min):** Kazama 1999 mostra que o efeito na pressão
  é bem mais lento (ke0 ≈ 0,07–0,12/min). O modelo deve prever a **PAM mínima cedo demais** (≈ 2 min) **(hipótese
  forte)**.
- **Fator s linear e aditivo:** idade, ASA, PVI e IVCI somam efeitos independentes; na realidade eles se
  sobrepõem (um idoso ASA III hipovolêmico não é a soma dos três).
- **FC só cai:** não há taquicardia reflexa.

**Como validar.** Este é o principal alvo do VitalDB (validacao.md, seção 7): ajustar Ce50, γ, Imax e ke0 da
pressão em 70% das induções e testar nos outros 30%.

---

## Etapa 6 — Paciente virtual: `src/calibracao.jl`

**O que faz.** Ajusta Rp (= PAM/DC) e a complacência C (por bisseção) para que o modelo reproduza a PAM e a
pressão de pulso do paciente; depois simula a indução com 1 min de linha de base e 20 min de observação.

| linhas | o que é | de onde vem |
|---|---|---|
| [calibracao.jl#L31](../src/calibracao.jl#L31) | faixa da complacência | **suposição** nossa |
| [calibracao.jl#L84](../src/calibracao.jl#L84) | calibração (Rp, C, PK, s) | nosso |
| [calibracao.jl#L206](../src/calibracao.jl#L206) | protocolo e linha do tempo | nosso, com os padrões do grupo |

**Vieses possíveis.**
- **C calibrada com o manguito** (ver etapa 2) **(hipótese)**.
- **Só PAM e pressão de pulso são calibradas:** PAS e PAD durante a queda saem da forma de onda fixa; no paciente C,
  a PA no momento da PAM mínima sai 92/22 mmHg — a PAD de 22 não é plausível (pendência 6.2).
- **Faixa 0,2–6 ml/mmHg sem fonte:** se a pressão de pulso exigir C fora dela, a calibração avisa ("interpretar com
  cautela").

**Como validar.** Com linha arterial pré-indução, comparar PAS, PAD e forma de onda do paciente virtual com as
medidas; checar quantos pacientes reais caem fora da faixa de C.

---

## Etapa 7 — Métricas e risco: `src/metricas.jl`

**O que faz.** Resume a curva: PAM basal e mínima, queda, tempo e área abaixo de 65 e 55 mmHg, BIS mínimo,
tempo até BIS < 60 e o risco.

| linhas | o que é | de onde vem |
|---|---|---|
| [metricas.jl#L38](../src/metricas.jl#L38) | cálculo das métricas | nosso |
| [metricas.jl#L85](../src/metricas.jl#L85) | ALTO < 55, MODERADO < 65, BAIXO ≥ 65 | limiares de Walsh 2013 e Salmasi 2017 |

**Vieses possíveis.**
- **Limiares absolutos:** num hipertenso crônico, 65 mmHg pode já ser uma queda grande; Salmasi 2017 também
  estudou limiares relativos (queda em % do basal). A regra atual não olha a queda relativa.
- **Risco do paciente típico:** o risco vem só da curva média; não considera a variabilidade.
- **Os limiares vêm de estudos de associação** (hipotensão e lesão de órgãos), não de causalidade.

**Como validar.** Discutir com o grupo se a regra deve incluir a queda relativa; no VitalDB, medir quantas
induções com risco ALTO/MODERADO/BAIXO previsto tiveram de fato PAM < 55/65.

---

## Etapa 8 — Janela de dose: `src/dose.jl`

**O que faz.** Varre doses e acha a menor dose que leva o BIS abaixo de 60 e a maior que mantém a PAM mínima ≥ 65.

| linhas | o que é | de onde vem |
|---|---|---|
| [dose.jl#L30](../src/dose.jl#L30) | simula várias doses | nosso |
| [dose.jl#L77](../src/dose.jl#L77) | critérios da janela | nosso; BIS < 60 (NICE DG6) e PAM ≥ 65 |

**Vieses possíveis.**
- **"BIS < 60 em algum momento dos 20 min"** conta a manutenção: até 0,5 mg/kg atinge o critério depois de vários
  minutos, então a dose mínima sai baixa demais (pendência 6.2: definir um tempo máximo para a hipnose).
- Herda todos os vieses das etapas 3 e 5.

---

## Etapa 9 — Plataforma: `src/plataforma.jl` e `app/index.html`

**O que faz.** Recebe o formulário, valida as faixas, calibra o paciente virtual (com cache), roda a simulação e
devolve JSON; a página desenha os gráficos.

| linhas | o que é |
|---|---|
| [plataforma.jl#L61](../src/plataforma.jl#L61) | faixas aceitas (18–95 anos, 35–200 kg, 130–215 cm, PA, FC, PVI, IVCI) |
| [plataforma.jl#L104](../src/plataforma.jl#L104) | dose, tempo de injeção, manutenção, modelo PK |
| [plataforma.jl#L163](../src/plataforma.jl#L163) | o que vai para a tela |
| [plataforma.jl#L229](../src/plataforma.jl#L229) | tabela de doses e janela |

**Vieses possíveis.** As faixas aceitas não garantem que o paciente esteja dentro do contexto de uso (ex.: um
obeso grave passa no filtro de peso). Conferir se cada número da tela é o mesmo do `relatorio` em Julia.

---

## Etapa 10 — Testes: `test/runtests.jl`

**O que faz.** Verifica as equações (soluções exatas, propriedades do Windkessel, Eleveld igual ao PAS) e a
plausibilidade (direção dos efeitos), e confere os links de validacao.md.

**Viés possível.** Verificação ≠ validação: os testes mostram que o código resolve as equações escolhidas, não que
as equações descrevem pacientes reais. A conferência contra o PAS só pega erros que o PAS não tenha também.

---

## Etapa 11 — Reprodução: `scripts/`

Rodar `julia --project=. scripts/rodar_tudo.jl` e conferir se os números dos logs em `resultados/` são os do README
(o modelo é determinístico). Os scripts não têm dados próprios: herdam as etapas anteriores.

---

## Resumo dos vieses, por impacto provável na curva de PAM

| # | viés | etapa | efeito provável | situação |
|---|---|---|---|---|
| 1 | Ce50, γ e Imax ilustrativos | 5 | profundidade da queda incerta | sem fonte; precisa de dados |
| 2 | ke0 da pressão rápido demais | 5 | PAM mínima prevista cedo demais | hipótese forte (Kazama 1999) |
| 3 | sem barorreflexo | 4 | queda e duração possivelmente superestimadas | hipótese |
| 4 | um Ce50 para Rp, VS e FC | 5 | efeito no VS subestimado no início | hipótese (Su 2022) |
| 5 | coeficientes de s aditivos e ilustrativos | 5 | erro maior em pacientes de risco | sem fonte |
| 6 | V1 fixo de Schnider | 3 | pico de Cp alto após o bolus | conhecido |
| 7 | manguito × aorta | 2, 6 | C e PAS/PAD enviesadas; PAM pouco afetada | hipótese |
| 8 | forma de onda única | 1, 6 | PAS/PAD implausíveis na queda | observado (PAD 22) |
| 9 | critério da dose mínima | 8 | dose mínima baixa demais | observado |
| 10 | obesos e populações fora dos estudos | 2, 3 | PK enviesado | fora do contexto de uso |
