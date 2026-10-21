# Manual do PropoSim — passo a passo

**PropoSim** é a plataforma do PBL 6 que monta um **paciente virtual** com os dados
pré-indução e prevê, **antes da injeção**, a curva de pressão arterial média (PAM) nos 20
minutos após um bolus de propofol. Por dentro: PK/PD do propofol (Schnider ou Marsh no
NeoPKPD, BIS de Eleveld) acoplada a um Windkessel de 4 elementos, como no repositório
PropofolPBL do grupo.

> **Uso educacional e de pesquisa.** Os parâmetros do efeito hemodinâmico são ilustrativos e
> o modelo ainda não foi validado com dados de pacientes. Não use para decidir a dose de um
> paciente real.

| | |
|---|---|
| Repositório | <https://github.com/vivianguilger-lab/PBL_6> |
| Download direto (ZIP) | <https://github.com/vivianguilger-lab/PBL_6/archive/refs/heads/claude/windkessel-propofol-hemodynamic-p72qsk.zip> |
| O que precisa | computador com Windows 10/11, macOS ou Linux; internet na primeira vez; ~3 GB livres |

---

## 1. Instalar (só na primeira vez)

Escolha **um** dos dois caminhos.

### Caminho A — duplo clique (recomendado)

**Windows**

1. Baixe o ZIP pelo link acima.
2. Na pasta Downloads, clique com o botão direito no ZIP → **Extrair tudo…** → **Extrair**.
   *Não rode nada de dentro do ZIP sem extrair.*
3. Abra a pasta extraída (`PBL_6-claude-windkessel-propofol-hemodynamic-p72qsk`) e dê
   **duplo clique em `abrir_plataforma.bat`**.
4. Se o Windows mostrar "O Windows protegeu o computador", clique em **Mais informações →
   Executar assim mesmo**.
5. Se o Julia não estiver instalado, o arquivo instala sozinho pela Microsoft Store
   (`winget`). Aceite as permissões que aparecerem.
6. Espere. **Na primeira vez o Julia baixa e compila os pacotes: 10 a 20 minutos.** A janela
   preta mostra o progresso. Quando terminar, o navegador abre a plataforma sozinho.

**macOS**

1. Baixe o ZIP e dê duplo clique para extrair.
2. Na pasta extraída, clique com o botão direito em **`abrir_plataforma.command`** → **Abrir**
   → **Abrir** (só na primeira vez, por causa do Gatekeeper).
3. O Julia é instalado automaticamente se faltar. Espere a compilação (10 a 20 min).

**Linux**: no terminal, dentro da pasta extraída, `./abrir_plataforma.command`.

### Caminho B — pelo próprio Julia

1. Instale o Julia (versão 1.10 ou mais nova):
   - Windows, no PowerShell: `winget install --name Julia --id 9NJNWW8PVKMN -e -s msstore`
   - macOS/Linux, no terminal: `curl -fsSL https://install.julialang.org | sh`
2. Abra o **Julia** (menu Iniciar → Julia). Aparece o prompt verde `julia>`.
3. Cole, uma linha por vez:

   ```julia
   using Pkg
   Pkg.add(url="https://github.com/vivianguilger-lab/PBL_6", rev="claude/windkessel-propofol-hemodynamic-p72qsk")
   ```

4. Espere terminar (10 a 20 min na primeira vez).

---

## 2. Abrir a plataforma (todas as vezes)

- **Caminho A:** duplo clique em `abrir_plataforma.bat` (Windows) ou
  `abrir_plataforma.command` (macOS/Linux). Abre em menos de 1 minuto.
- **Caminho B:** abra o Julia e digite

  ```julia
  using PropofolWindkessel
  plataforma()
  ```

O navegador abre em **http://localhost:8080** (se a porta estiver ocupada, a janela do Julia
mostra o endereço usado, por exemplo 8081). Se o navegador não abrir, copie esse endereço na
barra do navegador.

**Deixe a janela preta / do Julia aberta enquanto usa.** Para encerrar, feche essa janela
(ou aperte `Ctrl+C` nela).

Os dados digitados **não saem do computador**: o servidor só aceita conexões da própria
máquina.

---

## 3. Preencher os dados

A tela tem o formulário à esquerda e os resultados à direita. Para conhecer a plataforma,
clique num dos **exemplos** no topo do formulário (*62 a, hipertensa*, *28 a, ASA I*,
*81 a, aorta rígida*): ele preenche tudo e já simula.

### Dados pré-indução

| campo | unidade | de onde tirar | observação |
|---|---|---|---|
| Idade | anos | prontuário | 18 a 95 |
| Sexo | — | prontuário | entra na massa magra (Schnider) |
| Peso, Altura | kg, cm | prontuário | calculam IMC, massa magra e superfície corporal |

### Monitor

| campo | unidade | de onde tirar | observação |
|---|---|---|---|
| PAS, PAD | mmHg | manguito ou PAI **antes da indução**, paciente calmo | definem a PAM e a pressão de pulso do paciente virtual |
| FC | bpm | monitor | |
| PAM | mmHg | opcional | vazio = calculada (PAD + (PAS − PAD)/3) |
| DC | L/min | opcional (ecocardiograma, monitor de débito) | vazio = estimado (3,0 L/min/m² × superfície corporal) |

### Risco e volemia

| campo | unidade | de onde tirar | efeito no modelo |
|---|---|---|---|
| ASA | I–V | avaliação pré-anestésica | ASA III ou mais aumenta a sensibilidade ao propofol |
| PVI | % | oxímetro com índice de variabilidade pletismográfica | PVI > 13% sugere hipovolemia → queda maior |
| IVCI | % | ultrassom: índice de colapsabilidade da veia cava inferior | IVCI > 40% sugere hipovolemia → queda maior |
| SpO₂, Glicemia | %, mg/dL | monitor, glicemia capilar | registrados no relatório; ainda não entram nas equações |

Idade acima de 50 anos, ASA ≥ III, PVI > 13% e IVCI > 40% aumentam o **fator de
sensibilidade *s*** (aparece em *Paciente virtual*). *s* = 1 é o paciente de referência;
*s* = 1,5 significa efeito 50% maior do propofol sobre a circulação.

### Protocolo

| campo | opções | observação |
|---|---|---|
| Dose de propofol | 0,5 a 3,0 mg/kg (régua) ou digite | o total em mg aparece ao lado |
| Tempo de injeção | 10 s a 2 min | injeção mais lenta → pico menor de concentração |
| Manutenção | µg/kg/min | começa 1 min após o bolus; **0 = só bolus** |
| Modelo farmacocinético | Schnider ou Marsh | Schnider é o padrão |

Clique em **Simular**. O resultado sai em cerca de 1 segundo; a faixa de incerteza e a
comparação de doses terminam alguns segundos depois.

---

## 4. Ler os resultados

### Faixa de risco (topo)

| cor | quando aparece |
|---|---|
| **Baixo** (verde) | a curva típica não passa abaixo de 65 mmHg por 1 min e a PAM cai menos de 30% |
| **Moderado** (laranja) | PAM < 65 mmHg por ≥ 1 min **ou** queda > 30% da PAM basal |
| **Alto** (vermelho) | a curva típica cruza **55 mmHg** (limiar associado a lesão miocárdica e renal) |

### Indicadores

- **PAM: basal → mínima** — valor previsto no pior momento e quando acontece (minutos após o
  início do bolus).
- **Tempo com PAM < 65** — minutos abaixo de 65 mmHg nos 20 minutos simulados.
- **Chance de PAM < 65** e **Chance de PAM < 55** — porcentagem de 100 pacientes virtuais
  sorteados com a variabilidade entre pessoas (PK, ke0, Ce50, Imax) em que a PAM fica abaixo
  do limiar. É a resposta para "e se este paciente for mais sensível que a média?".
- **Hipnose (BIS < 60)** — em quanto tempo o BIS previsto fica abaixo de 60.

### Gráficos

- **Pressão arterial prevista** — PAM (linha), PAS–PAD (faixa clara) e a **faixa de 90%** da
  incerteza. Linhas tracejadas em 65 e 55 mmHg. Passe o mouse sobre a curva para ver os valores.
- **Frequência cardíaca**, **Débito cardíaco** — evolução prevista.
- **Profundidade anestésica (BIS)** — faixa-alvo 40–60 sombreada.
- **Concentração de propofol** — no plasma (Cp) e no sítio efetor (Ce), em µg/ml.

### Comparar doses

Simula o mesmo paciente com 11 doses e mostra, para cada uma, a PAM mínima, a queda, o tempo
abaixo de 65 e de 55 mmHg, o BIS mínimo e o risco. **Clique numa linha da tabela para simular
aquela dose.**

A **faixa de dose sugerida pelo modelo** (em mg/kg e em mg) vai da menor dose que leva o BIS
abaixo de 60 até a maior dose cuja PAM mínima prevista fica ≥ 65 mmHg. Se aparecer **"Nenhuma
dose desta faixa produz hipnose mantendo a PAM ≥ 65 mmHg"**, o modelo prevê hipotensão já na
dose necessária para a hipnose: considere titular a dose, coindução, otimização volêmica ou
vasopressor.

### Por que a pressão cai?

Mostra quanto cada efeito do propofol contribui, em % do basal, no momento de maior queda:

- **Rp** — resistência periférica (vasodilatação), o principal;
- **VS** — volume sistólico (depressão do miocárdio e menor pré-carga);
- **FC** — frequência cardíaca.

### Paciente virtual

Os parâmetros que o modelo ajustou para este paciente: PAM, PAS e PAD reproduzidas, DC e VS,
os quatro elementos do Windkessel (Rp, C, Rc, L), a constante de tempo τ = Rp·C, o fator de
sensibilidade *s*, a massa magra, a superfície corporal e os parâmetros PK (CL, V1, V2, V3,
ke0) e o Ce50 do BIS. Se a calibração não reproduzir a pressão de pulso, aparece
"interpretar com cautela".

---

## 5. Salvar e compartilhar

- **Imprimir / salvar PDF** — gera uma página com os dados, o risco e os gráficos. Na janela
  de impressão escolha "Salvar como PDF".
- **Baixar curva (CSV)** — arquivo `curva_prevista.csv` com a curva prevista (tempo, PAS,
  PAD, PAM, FC, DC), separado por ponto e vírgula, que abre direto no Excel em português.
- **Tema claro/escuro** — botão no canto superior direito.

---

## 6. Exemplo completo

Paciente: mulher, 62 anos, 72 kg, 160 cm, ASA II, PA 148/82 mmHg, FC 78 bpm, PVI 12%,
IVCI 25%. Propofol 1,5 mg/kg em 30 s + 100 µg/kg/min.

1. Clique no exemplo **62 a, hipertensa** (ou digite os dados acima) → **Simular**.
2. Resultado esperado: PAM **104 → 66 mmHg** por volta de 2 min, risco **Moderado**
   (queda > 30%), chance de PAM < 65 ≈ **40%** e de PAM < 55 ≈ **10%**.
3. Em *Comparar doses*, a faixa sugerida fica em torno de **0,50–1,63 mg/kg**.
4. Mude a dose para 2,5 mg/kg (régua) e simule de novo: a PAM mínima cai abaixo de 65 e o
   risco sobe.
5. Mude o tempo de injeção para 2 min: a queda fica menor, porque o pico de Ce é mais baixo.
6. Troque para o exemplo **81 a, aorta rígida** (ASA III, PVI 18%, IVCI 45%, *s* ≈ 1,66): com a
   mesma dose a PAM cai para ≈ 51 mmHg e o risco é **Alto**.

---

## 7. Atualizar para a versão mais nova

- **Caminho A:** baixe o ZIP de novo, extraia e use a pasta nova (pode apagar a antiga).
- **Caminho B:** no Julia, `using Pkg; Pkg.update()`.

---

## 8. Problemas comuns

| o que aparece | causa | solução |
|---|---|---|
| "Este arquivo precisa estar dentro da pasta do projeto" | o `.bat` foi aberto de dentro do ZIP | extraia o ZIP (botão direito → Extrair tudo) e rode de dentro da pasta extraída |
| "Julia não encontrado" e o `winget` falha | Windows sem Microsoft Store ou sem permissão | instale o Julia pelo site <https://julialang.org/downloads/>, feche e abra o `.bat` de novo |
| janela preta parada por muitos minutos | primeira compilação | normal na primeira vez (até 20 min); não feche |
| `ERROR: ... Project.toml not found` ou `Package PropofolWindkessel not found` | Julia aberto fora da pasta do projeto (caminho A) ou pacote não instalado (caminho B) | use o `.bat`, ou faça o passo 3 do caminho B |
| `ParseError` ao colar `cd`, `dir`, `julia --project` no Julia | comando de terminal digitado no prompt `julia>` | no Julia, digite só código Julia (seção 2); comandos de terminal vão no PowerShell |
| `Permission denied` / `Acesso negado` | Julia aberto em `C:\Windows\system32` tentando gravar ali | use o caminho A ou B deste manual (não precisam gravar nessa pasta) |
| navegador não abre | navegador padrão bloqueado | copie o endereço mostrado na janela (ex.: `http://localhost:8080`) |
| "address already in use" | a plataforma já está aberta em outra janela | use a janela que já está aberta, ou feche-a; a plataforma tenta 8081, 8082… sozinha |
| página carrega mas "Simular" dá erro | a janela do Julia foi fechada | abra a plataforma de novo (seção 2) |
| erro de rede/proxy ao instalar | internet institucional bloqueando downloads | instale de outra rede (ex.: celular roteando) |

Se nada resolver, mande ao grupo um print da janela preta inteira.

---

## 9. Para quem vai mexer no modelo

- Pasta do projeto: `src/` (modelo), `app/index.html` (plataforma), `scripts/` (análises e
  figuras), `docs/` (equações, revisão, credibilidade V&V 40), `test/` (testes).
- Equações: [docs/modelo_matematico.md](docs/modelo_matematico.md). Credibilidade e a
  sugestão de calibração (Imax × 0,85): [docs/validacao.md](docs/validacao.md).
- Usar no Julia sem a plataforma:

  ```julia
  using PropofolWindkessel
  p  = Paciente(idade = 62, sexo = :F, peso = 72, altura = 160, pas = 148, pad = 82, fc = 78,
                asa = 2, pvi = 12, ivci = 25)
  pv = calibrar(p)
  r  = simular_inducao(pv; dose_mgkg = 1.5)
  relatorio(metricas(r))
  janela_de_dose(pv)
  previsao_com_incerteza(pv; dose_mgkg = 1.5)
  ```

- Mudar os parâmetros do efeito: `calibrar(p; efeito = EfeitoHemodinamico(...))`
  (`src/acoplamento.jl`).
- Rodar todas as análises (na pasta do projeto, no terminal):
  `julia --project=. -t auto scripts/rodar_tudo.jl` → figuras e tabelas em `resultados/`.
- Testes: `julia --project=. -e "using Pkg; Pkg.test()"`.

---

## 10. Limitações

- Parâmetros do efeito hemodinâmico ilustrativos: na coorte virtual o modelo prevê PAM < 55 em
  35% dos adultos, contra 17–24% na literatura (superestima a hipotensão).
- Sem barorreflexo, retorno venoso explícito, opioides, laringoscopia, fluidos ou
  vasopressores.
- Um único ke0 para BIS e pressão (na prática o efeito na pressão é mais lento).
- PAS e PAD individuais não coincidem exatamente com o manguito (a forma de onda é fixa); a
  PAM e a pressão de pulso, sim.
- Fora do contexto de uso: crianças, obesos graves, choque, sepse, cardiopatia grave.
