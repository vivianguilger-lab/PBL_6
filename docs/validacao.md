# Validação e rastreabilidade do modelo

Para cada valor do modelo, este documento diz **onde ele está no código** (arquivo e linha), **de onde veio**
(repositório do grupo, biblioteca, artigo) e **se a fonte confirma o valor e a ordem de grandeza**. Também
registra o que foi corrigido em cada rodada de validação e o que ainda falta validar.

**Situação hoje:** a farmacocinética, o BIS e o Windkessel seguem modelos publicados (com a fonte primária
ainda a abrir em alguns itens); **o efeito do propofol sobre a pressão (Ce50, γ, Imax) e os coeficientes do
fator de sensibilidade não têm fonte nem calibração**. A curva de PAM prevista depende deles.

## 1. Como ler

**Links.** Os links para o nosso código são relativos: no GitHub abrem o arquivo da mesma versão deste
documento, já na linha. Links para arquivos `.md` usam `?plain=1`; sem isso o GitHub mostra o texto formatado
e ignora o número da linha. Os outros repositórios estão fixados num commit, então a linha não muda se alguém
editar o repositório depois:

| sigla | repositório | versão |
|---|---|---|
| — | este repositório ([vivianguilger-lab/PBL_6](https://github.com/vivianguilger-lab/PBL_6)) | a deste arquivo |
| grupo | [gabiadamek30-bit/Simula-o-da-resposta-do-paciente-ao-propofol---PBL](https://github.com/gabiadamek30-bit/Simula-o-da-resposta-do-paciente-ao-propofol---PBL/tree/b6995fc931af5a542bbc100f1677d93962f5b2b6) (PropofolPBL) | b6995fc |
| PAS | [AnesthesiaSimulation/Python_Anesthesia_Simulator](https://github.com/AnesthesiaSimulation/Python_Anesthesia_Simulator/tree/4afc03ad2997f725205c202e0d0634ff0664e5c7), conferência independente | 4afc03a |
| windkessel-id | [hpigot/windkessel-id](https://github.com/hpigot/windkessel-id) v1.0 (a que o grupo copiou) e v2.0 | d7cdc0f e 153cff9 |
| NeoPKPD | [shramish2057/NeoPKPD](https://github.com/shramish2057/NeoPKPD/tree/269e0c1b68d2e78f1b4c3bfeaeef3fd245c6efb4/packages/core) v0.1.0, pasta `packages/core` | 269e0c1 |

**Vereditos.**

| veredito | significado |
|---|---|
| CONFIRMADO | o valor é o mesmo da fonte, e a fonte foi aberta (código, dados ou texto) |
| CONFIRMADO\* | o valor confere com o resumo indexado do artigo ou com código aberto que reproduz a tabela, mas **o artigo em si não foi aberto** (a rede usada na conferência só alcançava o GitHub). Fecha quando alguém abrir o DOI e conferir a tabela indicada |
| DIFERENTE | a fonte dá outro valor (anotado na tabela) |
| ILUSTRATIVO | valor escolhido pelo grupo, sem fonte; precisa de calibração com dados |
| SUPOSIÇÃO | escolha nossa, sem fonte; fora da plataforma até ser validada |

**Origem da citação** (coluna da [bibliografia](#8-bibliografia)): **G** = citada no código do grupo;
**P** = citada no código do PAS; **W** = citada no windkessel-id; **B** = encontrada na nossa busca de literatura
(não estava em nenhum repositório).

## 2. Contexto de uso e risco (ASME V&V 40)

> Ferramenta **educacional e de pesquisa** para explorar, **antes da injeção**, a curva de PAM prevista nos
> primeiros 20 min após um bolus de propofol (com ou sem manutenção) em adultos de 18 a 95 anos, sem estímulo
> cirúrgico e **sem opioide**, a partir de idade, sexo, peso, altura, ASA, PAS, PAD, PAM, FC, PVI, IVCI e, se
> houver, DC.
>
> **Fora do contexto de uso:** decidir a dose de um paciente real; crianças; obesos graves (a massa magra de
> James usada por Schnider falha); choque, sepse, cardiopatia grave; coindução com opioides; manutenção com
> halogenados.

| fator | avaliação |
|---|---|
| influência do modelo na decisão | baixa no contexto declarado (apoio didático) |
| consequência de uma decisão errada | alta (hipotensão grave → lesão miocárdica e renal) |
| **risco do modelo** | **médio-alto** → exige calibração e validação com dados humanos antes de qualquer uso clínico |

## 3. O que mudou na validação (rodada 1)

Alterações feitas a partir da conferência linha a linha e dos comentários do grupo. Os números "antes → depois"
são dos pacientes de exemplo A (62 a, 1,5 mg/kg), B (28 a, 2 mg/kg) e C (81 a, 1,5 mg/kg).

| # | alteração | onde | por quê | efeito |
|---|---|---|---|---|
| 1 | **γ do BIS de Eleveld corrigido** para 1,89 abaixo do Ce50 e 1,47 acima (estava invertido, também no grupo) | [propofol.jl#L283](../src/propofol.jl#L283); grupo [grupo propofol.jl#L60](https://github.com/gabiadamek30-bit/Simula-o-da-resposta-do-paciente-ao-propofol---PBL/blob/b6995fc931af5a542bbc100f1677d93962f5b2b6/src/propofol.jl#L60) | Eleveld 2018, Tabela 3; PAS [PAS pd_models.py#L226-L228](https://github.com/AnesthesiaSimulation/Python_Anesthesia_Simulator/blob/4afc03ad2997f725205c202e0d0634ff0664e5c7/src/python_anesthesia_simulator/pd_models.py#L226-L228) | sozinho, sobe o BIS mínimo em ~5 pontos (A 14 → 19) |
| 2 | **BIS pelo modelo de Eleveld completo**: PK de Eleveld, ke0 = 0,146·(peso/70)^−0,25 /min e atraso de 15 + e^(0,0517·idade) s | [propofol.jl#L175-L257](../src/propofol.jl#L175-L257), [propofol.jl#L296](../src/propofol.jl#L296), [propofol.jl#L370-L380](../src/propofol.jl#L370-L380) | o Ce50 de 3,08 µg/ml só vale com o PK e o ke0 com que foi estimado; antes o grupo usava o Ce do ke0 de Schnider/Marsh e o próprio grupo chamava isso de aproximação ([grupo README.md#L69](https://github.com/gabiadamek30-bit/Simula-o-da-resposta-do-paciente-ao-propofol---PBL/blob/b6995fc931af5a542bbc100f1677d93962f5b2b6/README.md?plain=1#L69)) | BIS mínimo A 14 → 45, B 12 → 46, C 10 → 39; BIS < 60 em 2,2 / 1,7 / 2,7 min (antes 0,4 min) |
| 3 | **pisos removidos** em V2, CL e Q2 de Schnider; fora da faixa das equações o programa dá erro | [propofol.jl#L75](../src/propofol.jl#L75) | os pisos (1 L, 0,2 L/min, 0,1 L/min) eram nossos, sem fonte; dentro da faixa aceita pela plataforma (18–95 anos, 35–200 kg, 130–215 cm) nunca eram atingidos | nenhum nos exemplos |
| 4 | **citação do ke0 do Marsh** corrigida: 0,26/min é das bombas Diprifusor, não de Marsh 1991 | [propofol.jl#L50-L51](../src/propofol.jl#L50-L51) | Absalom et al. 2009 | só texto |
| 5 | **PAM obrigatória** (a do monitor); saiu a fórmula PAD + (PAS − PAD)/3 | [paciente.jl#L77](../src/paciente.jl#L77), [paciente.jl#L140](../src/paciente.jl#L140) | os monitores oscilométricos medem a PAM diretamente; a fórmula era uma aproximação sem fonte no código; o grupo já recebia a PAM como entrada ([grupo paciente.jl#L21](https://github.com/gabiadamek30-bit/Simula-o-da-resposta-do-paciente-ao-propofol---PBL/blob/b6995fc931af5a542bbc100f1677d93962f5b2b6/src/paciente.jl#L21)) | nenhum quando a PAM digitada = a calculada |
| 6 | **limiares do fator de sensibilidade**: PVI > 15% (era 13%) e IVCI > 43% (era 40%) | [acoplamento.jl#L114](../src/acoplamento.jl#L114) | PVI > 15 prevê queda da PAM após propofol (Tsuchiya 2010); IVCI ≥ 43% prevê hipotensão pós-indução (Zhang & Critchley 2016) | C: s 1,66 → 1,57; PAM mínima 51 → 54 mmHg |
| 7 | **regra de risco** só com os limiares que têm fonte: ALTO se a PAM mínima < 55, MODERADO se < 65, BAIXO se ≥ 65 | [metricas.jl#L85](../src/metricas.jl#L85) | saíram "queda > 30%" e "≥ 1 min abaixo de 65", que não tinham fonte | A: MODERADO → BAIXO (PAM mínima 65,8) |
| 8 | **faixa de incerteza e "chance de PAM < 65/55" fora da plataforma**; coorte virtual fora dos resultados | [incerteza.jl#L5](../src/incerteza.jl#L5), [populacao.jl#L5](../src/populacao.jl#L5), [rodar_tudo.jl#L49](../scripts/rodar_tudo.jl#L49) | dependem de variabilidades sem validação (seção 6) | a plataforma mostra só a previsão do paciente típico |
| 9 | **afirmações sem fonte retiradas**: "17–24% dos adultos com PAM < 55" e "queda típica de 25–40%" | [metricas.jl#L8](../src/metricas.jl#L8) (comentário reescrito) | não achadas na literatura nesses termos | só texto |
| 10 | **SpO₂ e glicemia fora da plataforma** | [paciente.jl#L47](../src/paciente.jl#L47) | não entram em nenhuma equação; ficam na estrutura só para manter a do grupo | — |
| 11 | **x0 do Windkessel** convertido de L/min para mL/s, como na v2.0 do windkessel-id | [windkessel.jl#L395](../src/windkessel.jl#L395); v2.0 [windkessel-id v2.0 WindkesselIdentification.jl#L218](https://github.com/hpigot/windkessel-id/blob/153cff968c05a83d4f6f8e2da914698cf1decde5/src/WindkesselIdentification.jl#L218) | o resto do código está em mL/s | nenhum (o x0 é recalculado por `solvex0`) |
| 12 | "porco de 70 kg" → "porcos de 65 kg" nos dados in vivo | [README.md#L10](../dados/windkessel-id/README.md?plain=1#L10) | correção da v2.0 do windkessel-id ([windkessel-id v2.0 README.md#L78](https://github.com/hpigot/windkessel-id/blob/153cff968c05a83d4f6f8e2da914698cf1decde5/README.md?plain=1#L78)) | só texto |
| 13 | opção **Eleveld** também para a PK que move a hemodinâmica | [propofol.jl#L91](../src/propofol.jl#L91) | modelo publicado, agora implementado | opcional; Schnider continua o padrão |

## 4. Rastreabilidade

### 4.1 Farmacocinética

Schnider e a massa magra conferem nos três códigos e nas reproduções da tabela; Eleveld confere com o PAS em
20 casos (teste automatizado); o ke0 = 0,26 do Marsh não vem do artigo de Marsh.

| valor | nosso código | grupo | conferência independente | fonte primária (onde) | veredito |
|---|---|---|---|---|---|
| Schnider V1 = 4,27 L | [propofol.jl#L66](../src/propofol.jl#L66) | [grupo propofol.jl#L11](https://github.com/gabiadamek30-bit/Simula-o-da-resposta-do-paciente-ao-propofol---PBL/blob/b6995fc931af5a542bbc100f1677d93962f5b2b6/src/propofol.jl#L11) | [PAS pk_models.py#L151](https://github.com/AnesthesiaSimulation/Python_Anesthesia_Simulator/blob/4afc03ad2997f725205c202e0d0634ff0664e5c7/src/python_anesthesia_simulator/pk_models.py#L151) | Schnider 1998, tabela de parâmetros finais [1] | CONFIRMADO\* |
| V2 = 18,9 − 0,391(idade − 53) L | [propofol.jl#L67](../src/propofol.jl#L67) | [grupo propofol.jl#L12](https://github.com/gabiadamek30-bit/Simula-o-da-resposta-do-paciente-ao-propofol---PBL/blob/b6995fc931af5a542bbc100f1677d93962f5b2b6/src/propofol.jl#L12) | [PAS pk_models.py#L152](https://github.com/AnesthesiaSimulation/Python_Anesthesia_Simulator/blob/4afc03ad2997f725205c202e0d0634ff0664e5c7/src/python_anesthesia_simulator/pk_models.py#L152) | idem [1] | CONFIRMADO\* |
| V3 = 238 L | [propofol.jl#L68](../src/propofol.jl#L68) | [grupo propofol.jl#L13](https://github.com/gabiadamek30-bit/Simula-o-da-resposta-do-paciente-ao-propofol---PBL/blob/b6995fc931af5a542bbc100f1677d93962f5b2b6/src/propofol.jl#L13) | [PAS pk_models.py#L153](https://github.com/AnesthesiaSimulation/Python_Anesthesia_Simulator/blob/4afc03ad2997f725205c202e0d0634ff0664e5c7/src/python_anesthesia_simulator/pk_models.py#L153) | idem [1] | CONFIRMADO\* |
| CL = 1,89 + 0,0456(peso − 77) − 0,0681(LBM − 59) + 0,0264(altura − 177) L/min | [propofol.jl#L69](../src/propofol.jl#L69) | [grupo propofol.jl#L14](https://github.com/gabiadamek30-bit/Simula-o-da-resposta-do-paciente-ao-propofol---PBL/blob/b6995fc931af5a542bbc100f1677d93962f5b2b6/src/propofol.jl#L14) | [PAS pk_models.py#L146-L147](https://github.com/AnesthesiaSimulation/Python_Anesthesia_Simulator/blob/4afc03ad2997f725205c202e0d0634ff0664e5c7/src/python_anesthesia_simulator/pk_models.py#L146-L147); [tci](https://raw.githubusercontent.com/jarretrt/tci/master/R/poppk_mods.R) | idem [1] | CONFIRMADO\* |
| Q2 = 1,29 − 0,024(idade − 53) L/min | [propofol.jl#L70](../src/propofol.jl#L70) | [grupo propofol.jl#L15](https://github.com/gabiadamek30-bit/Simula-o-da-resposta-do-paciente-ao-propofol---PBL/blob/b6995fc931af5a542bbc100f1677d93962f5b2b6/src/propofol.jl#L15) | [PAS pk_models.py#L148](https://github.com/AnesthesiaSimulation/Python_Anesthesia_Simulator/blob/4afc03ad2997f725205c202e0d0634ff0664e5c7/src/python_anesthesia_simulator/pk_models.py#L148) | idem [1] | CONFIRMADO\* |
| Q3 = 0,836 L/min | [propofol.jl#L71](../src/propofol.jl#L71) | [grupo propofol.jl#L16](https://github.com/gabiadamek30-bit/Simula-o-da-resposta-do-paciente-ao-propofol---PBL/blob/b6995fc931af5a542bbc100f1677d93962f5b2b6/src/propofol.jl#L16) | [PAS pk_models.py#L149](https://github.com/AnesthesiaSimulation/Python_Anesthesia_Simulator/blob/4afc03ad2997f725205c202e0d0634ff0664e5c7/src/python_anesthesia_simulator/pk_models.py#L149) | idem [1] | CONFIRMADO\* |
| Schnider ke0 = 0,456 /min | [propofol.jl#L72](../src/propofol.jl#L72) | [grupo propofol.jl#L17](https://github.com/gabiadamek30-bit/Simula-o-da-resposta-do-paciente-ao-propofol---PBL/blob/b6995fc931af5a542bbc100f1677d93962f5b2b6/src/propofol.jl#L17) | [PAS pk_models.py#L155](https://github.com/AnesthesiaSimulation/Python_Anesthesia_Simulator/blob/4afc03ad2997f725205c202e0d0634ff0664e5c7/src/python_anesthesia_simulator/pk_models.py#L155) | Schnider 1999, Resultados [2] | CONFIRMADO\* |
| Marsh V1 = 0,228 L/kg | [propofol.jl#L84](../src/propofol.jl#L84) | [grupo propofol.jl#L19](https://github.com/gabiadamek30-bit/Simula-o-da-resposta-do-paciente-ao-propofol---PBL/blob/b6995fc931af5a542bbc100f1677d93962f5b2b6/src/propofol.jl#L19) | [PAS pk_models.py#L180](https://github.com/AnesthesiaSimulation/Python_Anesthesia_Simulator/blob/4afc03ad2997f725205c202e0d0634ff0664e5c7/src/python_anesthesia_simulator/pk_models.py#L180) | Marsh 1991 [3]; 15,9 L para 70 kg em Absalom 2009 [4] | CONFIRMADO\* |
| Marsh k10 0,119; k12 0,112; k21 0,055; k31 0,0033 /min | [propofol.jl#L85-L89](../src/propofol.jl#L85-L89) | [grupo propofol.jl#L20-L24](https://github.com/gabiadamek30-bit/Simula-o-da-resposta-do-paciente-ao-propofol---PBL/blob/b6995fc931af5a542bbc100f1677d93962f5b2b6/src/propofol.jl#L20-L24) | [PAS pk_models.py#L181-L185](https://github.com/AnesthesiaSimulation/Python_Anesthesia_Simulator/blob/4afc03ad2997f725205c202e0d0634ff0664e5c7/src/python_anesthesia_simulator/pk_models.py#L181-L185) (V2 e V3 já em L/kg) | Marsh 1991 [3] | CONFIRMADO\* |
| Marsh k13 = 0,0419 /min | [propofol.jl#L88](../src/propofol.jl#L88) | [grupo propofol.jl#L23](https://github.com/gabiadamek30-bit/Simula-o-da-resposta-do-paciente-ao-propofol---PBL/blob/b6995fc931af5a542bbc100f1677d93962f5b2b6/src/propofol.jl#L23) | PAS usa 0,042: [PAS pk_models.py#L185](https://github.com/AnesthesiaSimulation/Python_Anesthesia_Simulator/blob/4afc03ad2997f725205c202e0d0634ff0664e5c7/src/python_anesthesia_simulator/pk_models.py#L185) | Marsh 1991 [3] | DIFERENTE só no arredondamento (0,3%) |
| Marsh ke0 = 0,26 /min | [propofol.jl#L90](../src/propofol.jl#L90) | [grupo propofol.jl#L25](https://github.com/gabiadamek30-bit/Simula-o-da-resposta-do-paciente-ao-propofol---PBL/blob/b6995fc931af5a542bbc100f1677d93962f5b2b6/src/propofol.jl#L25) | [PAS pk_models.py#L188](https://github.com/AnesthesiaSimulation/Python_Anesthesia_Simulator/blob/4afc03ad2997f725205c202e0d0634ff0664e5c7/src/python_anesthesia_simulator/pk_models.py#L188) ("Marsh modificado" 1,2 em [PAS pk_models.py#L187-L190](https://github.com/AnesthesiaSimulation/Python_Anesthesia_Simulator/blob/4afc03ad2997f725205c202e0d0634ff0664e5c7/src/python_anesthesia_simulator/pk_models.py#L187-L190)) | **não está em Marsh 1991**: valor das bombas Diprifusor, nunca publicado em revista [4]; alternativa: 1,21 /min de Struys 2000 [5] | DIFERENTE na atribuição (corrigida) |
| Eleveld: θ de referência (V1 6,28; V2 25,5; V3 273 L; CL 1,79 (♂) / 2,10 (♀); Q2 1,75; Q3 1,11 L/min) e covariáveis | [propofol.jl#L199-L202](../src/propofol.jl#L199-L202) | não existe | [PAS pk_models.py#L270-L287](https://github.com/AnesthesiaSimulation/Python_Anesthesia_Simulator/blob/4afc03ad2997f725205c202e0d0634ff0664e5c7/src/python_anesthesia_simulator/pk_models.py#L270-L287) | Eleveld 2018, parâmetros PK [6] | CONFIRMADO\* (igual ao PAS em 20 casos, erro relativo < 1e−12: [runtests.jl#L193](../test/runtests.jl#L193)) |
| Eleveld: massa livre de gordura de Al-Sallami | [propofol.jl#L218-L219](../src/propofol.jl#L218-L219) | não existe | [PAS pk_models.py#L294-L298](https://github.com/AnesthesiaSimulation/Python_Anesthesia_Simulator/blob/4afc03ad2997f725205c202e0d0634ff0664e5c7/src/python_anesthesia_simulator/pk_models.py#L294-L298) | Al-Sallami 2015 [7], usada por Eleveld [6] | CONFIRMADO\* (via PAS) |
| Eleveld ke0 = 0,146·(peso/70)^−0,25 /min (amostras arteriais) | [propofol.jl#L256](../src/propofol.jl#L256) | não existe | [PAS pk_models.py#L335](https://github.com/AnesthesiaSimulation/Python_Anesthesia_Simulator/blob/4afc03ad2997f725205c202e0d0634ff0664e5c7/src/python_anesthesia_simulator/pk_models.py#L335) | Eleveld 2018, Tabela 3 [6] | CONFIRMADO\* |
| Massa magra (James): ♂ 1,1P − 128(P/A)²; ♀ 1,07P − 148(P/A)², A em cm | [paciente.jl#L124](../src/paciente.jl#L124) | [grupo paciente.jl#L35-L37](https://github.com/gabiadamek30-bit/Simula-o-da-resposta-do-paciente-ao-propofol---PBL/blob/b6995fc931af5a542bbc100f1677d93962f5b2b6/src/paciente.jl#L35-L37) | [PAS patient.py#L217-L219](https://github.com/AnesthesiaSimulation/Python_Anesthesia_Simulator/blob/4afc03ad2997f725205c202e0d0634ff0664e5c7/src/python_anesthesia_simulator/patient.py#L217-L219) | James 1976 [8] (sem versão online); usada no método de Schnider 1998 [1] | CONFIRMADO\* |

**Ordem de grandeza:** volumes em L e depurações em L/min nos três códigos (o PAS converte para segundos depois).
Para a paciente A, Schnider dá V1 = 4,27 L e Eleveld 6,34 L, e o pico de concentração depois do bolus é bem
maior no Schnider (limitação do V1 fixo, já anotada pelo grupo em [grupo README.md#L68](https://github.com/gabiadamek30-bit/Simula-o-da-resposta-do-paciente-ao-propofol---PBL/blob/b6995fc931af5a542bbc100f1677d93962f5b2b6/README.md?plain=1#L68)).

**Motor de equações (NeoPKPD):** modelo tricompartimental [NeoPKPD specs.jl#L200-L202](https://github.com/shramish2057/NeoPKPD/blob/269e0c1b68d2e78f1b4c3bfeaeef3fd245c6efb4/packages/core/src/specs/specs.jl#L200-L202);
integração [NeoPKPD solve.jl#L407](https://github.com/shramish2057/NeoPKPD/blob/269e0c1b68d2e78f1b4c3bfeaeef3fd245c6efb4/packages/core/src/engine/solve.jl#L407); sítio efetor dCe/dt = ke0(Cp − Ce)
[NeoPKPD biophase_equilibration.jl#L81](https://github.com/shramish2057/NeoPKPD/blob/269e0c1b68d2e78f1b4c3bfeaeef3fd245c6efb4/packages/core/src/pd/biophase_equilibration.jl#L81). O `max(…, 0)` em Cp e Ce
([propofol.jl#L322](../src/propofol.jl#L322)) só tira ruído numérico do integrador (valores da ordem
de −1e−12). Verificação contra a solução exata exp(K·t)·A₀: [runtests.jl#L220](../test/runtests.jl#L220).

### 4.2 BIS (Eleveld 2018)

| valor | nosso código | grupo | conferência independente | fonte primária (onde) | veredito |
|---|---|---|---|---|---|
| E0 = Emax = 93 | [propofol.jl#L277](../src/propofol.jl#L277) | [grupo propofol.jl#L58](https://github.com/gabiadamek30-bit/Simula-o-da-resposta-do-paciente-ao-propofol---PBL/blob/b6995fc931af5a542bbc100f1677d93962f5b2b6/src/propofol.jl#L58) | [PAS pd_models.py#L232](https://github.com/AnesthesiaSimulation/Python_Anesthesia_Simulator/blob/4afc03ad2997f725205c202e0d0634ff0664e5c7/src/python_anesthesia_simulator/pd_models.py#L232) | Eleveld 2018, Tabela 3, θ "Baseline BIS" = 93,0 [6] | CONFIRMADO\* |
| Ce50 = 3,08·e^(−0,00635·(idade − 35)) µg/ml | [propofol.jl#L280](../src/propofol.jl#L280) | [grupo propofol.jl#L59](https://github.com/gabiadamek30-bit/Simula-o-da-resposta-do-paciente-ao-propofol---PBL/blob/b6995fc931af5a542bbc100f1677d93962f5b2b6/src/propofol.jl#L59) | idade de referência [PAS pd_models.py#L215](https://github.com/AnesthesiaSimulation/Python_Anesthesia_Simulator/blob/4afc03ad2997f725205c202e0d0634ff0664e5c7/src/python_anesthesia_simulator/pd_models.py#L215), função [PAS pd_models.py#L218](https://github.com/AnesthesiaSimulation/Python_Anesthesia_Simulator/blob/4afc03ad2997f725205c202e0d0634ff0664e5c7/src/python_anesthesia_simulator/pd_models.py#L218), valor [PAS pd_models.py#L222](https://github.com/AnesthesiaSimulation/Python_Anesthesia_Simulator/blob/4afc03ad2997f725205c202e0d0634ff0664e5c7/src/python_anesthesia_simulator/pd_models.py#L222) | Tabela 3: Ce50 = 3,08 (IC 2,91–3,26); envelhecimento −0,00635 [6] | CONFIRMADO\* |
| γ = 1,89 se Ce ≤ Ce50; 1,47 acima | [propofol.jl#L283](../src/propofol.jl#L283) | **invertido**: [grupo propofol.jl#L60](https://github.com/gabiadamek30-bit/Simula-o-da-resposta-do-paciente-ao-propofol---PBL/blob/b6995fc931af5a542bbc100f1677d93962f5b2b6/src/propofol.jl#L60) | [PAS pd_models.py#L226-L228](https://github.com/AnesthesiaSimulation/Python_Anesthesia_Simulator/blob/4afc03ad2997f725205c202e0d0634ff0664e5c7/src/python_anesthesia_simulator/pd_models.py#L226-L228), uso [PAS pd_models.py#L343](https://github.com/AnesthesiaSimulation/Python_Anesthesia_Simulator/blob/4afc03ad2997f725205c202e0d0634ff0664e5c7/src/python_anesthesia_simulator/pd_models.py#L343); [propofol-pebble-induction](https://raw.githubusercontent.com/ysuzuki1978/propofol-pebble-induction/main/validation/eleveld-models.reference.js) igual | Tabela 3: 1,89 "(Ce < Ce50)" e 1,47 "(Ce > Ce50)" [6] | CONFIRMADO\* (corrigido nesta rodada) |
| atraso do BIS = 15 + e^(0,0517·idade) s | [propofol.jl#L296](../src/propofol.jl#L296) | não existe | [PAS pd_models.py#L219](https://github.com/AnesthesiaSimulation/Python_Anesthesia_Simulator/blob/4afc03ad2997f725205c202e0d0634ff0664e5c7/src/python_anesthesia_simulator/pd_models.py#L219), [PAS pd_models.py#L236](https://github.com/AnesthesiaSimulation/Python_Anesthesia_Simulator/blob/4afc03ad2997f725205c202e0d0634ff0664e5c7/src/python_anesthesia_simulator/pd_models.py#L236) | Tabela 3, θ = 0,0517 [6] | CONFIRMADO\* (implementado nesta rodada) |
| ke0 do BIS | [propofol.jl#L370](../src/propofol.jl#L370) (ke0 de Eleveld) | Ce do ke0 do PK escolhido (aproximação declarada no [grupo README.md#L69](https://github.com/gabiadamek30-bit/Simula-o-da-resposta-do-paciente-ao-propofol---PBL/blob/b6995fc931af5a542bbc100f1677d93962f5b2b6/README.md?plain=1#L69)) | PAS usa o ke0 de Eleveld | Tabela 3 [6] | CONFIRMADO\* |

O PAS traduz exatamente a mesma equação do Ce50: `faging(x) = exp(x·(idade − 35))` e `c50p = 3,08·faging(−0,00635)`
= 3,08·e^(−0,00635·(idade − 35)). A `fdelay` da linha seguinte é outra coisa: o atraso do BIS.

**Ordem de grandeza (checagem clínica):** com o modelo antigo, 1,5–2 mg/kg levavam o BIS a 10–14, o que não é
plausível para uma dose de indução; com o Eleveld completo, o BIS mínimo fica em 39–46 (faixa-alvo 40–60 [19]),
atingido em 1,7–2,7 min.

### 4.3 Windkessel de 4 elementos

| valor | nosso código | grupo | windkessel-id | fonte primária (onde) | veredito |
|---|---|---|---|---|---|
| Rp = 0,79 mmHg·s/ml | [windkessel.jl#L386](../src/windkessel.jl#L386) | [grupo windkessel.jl#L209](https://github.com/gabiadamek30-bit/Simula-o-da-resposta-do-paciente-ao-propofol---PBL/blob/b6995fc931af5a542bbc100f1677d93962f5b2b6/src/windkessel.jl#L209) | v1.0 [windkessel-id v1.0 oeid.jl#L160](https://github.com/hpigot/windkessel-id/blob/d7cdc0f51bc77ac2070c0e0ff7ad053ab2cd283c/oeid.jl#L160); v2.0 [windkessel-id v2.0 WindkesselIdentification.jl#L213](https://github.com/hpigot/windkessel-id/blob/153cff968c05a83d4f6f8e2da914698cf1decde5/src/WindkesselIdentification.jl#L213) | Stergiopulos 1999, Tabela 1, batimento tipo A [9] | CONFIRMADO\* (via código que cita a tabela) |
| C = 1,22 ml/mmHg | [windkessel.jl#L387](../src/windkessel.jl#L387) | [grupo windkessel.jl#L210](https://github.com/gabiadamek30-bit/Simula-o-da-resposta-do-paciente-ao-propofol---PBL/blob/b6995fc931af5a542bbc100f1677d93962f5b2b6/src/windkessel.jl#L210) | [windkessel-id v1.0 oeid.jl#L161](https://github.com/hpigot/windkessel-id/blob/d7cdc0f51bc77ac2070c0e0ff7ad053ab2cd283c/oeid.jl#L161) | idem [9] | CONFIRMADO\* |
| Rc = 0,056 mmHg·s/ml | [windkessel.jl#L388](../src/windkessel.jl#L388) | [grupo windkessel.jl#L211](https://github.com/gabiadamek30-bit/Simula-o-da-resposta-do-paciente-ao-propofol---PBL/blob/b6995fc931af5a542bbc100f1677d93962f5b2b6/src/windkessel.jl#L211) | [windkessel-id v1.0 oeid.jl#L162](https://github.com/hpigot/windkessel-id/blob/d7cdc0f51bc77ac2070c0e0ff7ad053ab2cd283c/oeid.jl#L162) | idem [9] | CONFIRMADO\* |
| L = 0,0051 mmHg·s²/ml | [windkessel.jl#L389](../src/windkessel.jl#L389) | [grupo windkessel.jl#L212](https://github.com/gabiadamek30-bit/Simula-o-da-resposta-do-paciente-ao-propofol---PBL/blob/b6995fc931af5a542bbc100f1677d93962f5b2b6/src/windkessel.jl#L212) | [windkessel-id v1.0 oeid.jl#L163](https://github.com/hpigot/windkessel-id/blob/d7cdc0f51bc77ac2070c0e0ff7ad053ab2cd283c/oeid.jl#L163) — o comentário diz "impedance" (impedância), mas L é a **inertância**; a v2.0 corrigiu para "inertance" ([windkessel-id v2.0 WindkesselIdentification.jl#L216](https://github.com/hpigot/windkessel-id/blob/153cff968c05a83d4f6f8e2da914698cf1decde5/src/WindkesselIdentification.jl#L216)) | idem [9] | CONFIRMADO\* |
| unidades (mmHg, ml, s) | [windkessel.jl#L13](../src/windkessel.jl#L13) | [grupo windkessel.jl#L11](https://github.com/gabiadamek30-bit/Simula-o-da-resposta-do-paciente-ao-propofol---PBL/blob/b6995fc931af5a542bbc100f1677d93962f5b2b6/src/windkessel.jl#L11) | [windkessel-id v1.0 README.md#L17](https://github.com/hpigot/windkessel-id/blob/d7cdc0f51bc77ac2070c0e0ff7ad053ab2cd283c/README.md?plain=1#L17) | — | CONFIRMADO |
| citação dos valores | — | — | [windkessel-id v1.0 id_stergiopulos1999total.jl#L148](https://github.com/hpigot/windkessel-id/blob/d7cdc0f51bc77ac2070c0e0ff7ad053ab2cd283c/id_stergiopulos1999total.jl#L148) | [9] | — |
| θ identificado por nós [0,80; 1,25; 0,056; 0,0055] | [README.md#L116](../README.md?plain=1#L116) | — | v2.0 [windkessel-id v2.0 README.md#L36](https://github.com/hpigot/windkessel-id/blob/153cff968c05a83d4f6f8e2da914698cf1decde5/README.md?plain=1#L36) | — | CONFIRMADO (mesmo ajuste) |
| fluxo e pressão aórticos (166 amostras a 5 ms) | [aflow.csv](../dados/windkessel-id/stergiopulos1999_data/aflow.csv), [apressure.csv](../dados/windkessel-id/stergiopulos1999_data/apressure.csv) | cópia em `external/` | [v1.0 aflow.csv](https://github.com/hpigot/windkessel-id/blob/d7cdc0f51bc77ac2070c0e0ff7ad053ab2cd283c/stergiopulos1999_data/aflow.csv) (idêntico byte a byte) | Fig. 4A de [9] | CONFIRMADO |
| passo h = 0,005 s | [acoplamento.jl#L241](../src/acoplamento.jl#L241) | [grupo acoplamento.jl#L87](https://github.com/gabiadamek30-bit/Simula-o-da-resposta-do-paciente-ao-propofol---PBL/blob/b6995fc931af5a542bbc100f1677d93962f5b2b6/src/acoplamento.jl#L87) | [windkessel-id v1.0 id_stergiopulos1999total.jl#L14](https://github.com/hpigot/windkessel-id/blob/d7cdc0f51bc77ac2070c0e0ff7ad053ab2cd283c/id_stergiopulos1999total.jl#L14) | — | CONFIRMADO |
| fluxo escalado por VS e FC; VS0 = DC·1000/FC | [acoplamento.jl#L285](../src/acoplamento.jl#L285) | [grupo acoplamento.jl#L93](https://github.com/gabiadamek30-bit/Simula-o-da-resposta-do-paciente-ao-propofol---PBL/blob/b6995fc931af5a542bbc100f1677d93962f5b2b6/src/acoplamento.jl#L93) | — | — | igual ao grupo |
| Rp do paciente = PAM/DC | [acoplamento.jl#L143](../src/acoplamento.jl#L143) | [grupo acoplamento.jl#L44](https://github.com/gabiadamek30-bit/Simula-o-da-resposta-do-paciente-ao-propofol---PBL/blob/b6995fc931af5a542bbc100f1677d93962f5b2b6/src/acoplamento.jl#L44) | — | PAM = Rp·Q̄ (\|Z(0)\| = Rp), verificado em [runtests.jl#L124](../test/runtests.jl#L124) | CONFIRMADO |
| C do paciente ajustada entre 0,2 e 6 ml/mmHg para reproduzir PAS − PAD | [calibracao.jl#L31](../src/calibracao.jl#L31) | não existe (C fixo em 1,22) | — | — | SUPOSIÇÃO de faixa; ver seção 6 |

Ordem de grandeza: Rp = 0,79 mmHg·s/ml ≈ 1 050 dyn·s·cm⁻⁵, dentro da faixa fisiológica da resistência vascular
sistêmica (800–1 200); Rc/Rp ≈ 7%.

### 4.4 Efeito hemodinâmico e sensibilidade

**Nenhum destes valores tem fonte publicada.** O grupo escolheu-os para dar a ordem de grandeza da queda de PAM
na indução e escreveu isso no código ([grupo acoplamento.jl#L10-L11](https://github.com/gabiadamek30-bit/Simula-o-da-resposta-do-paciente-ao-propofol---PBL/blob/b6995fc931af5a542bbc100f1677d93962f5b2b6/src/acoplamento.jl#L10-L11)).
A coluna "para comparar" mostra o que existe de mais próximo na literatura.

| valor | nosso código | grupo | para comparar | veredito |
|---|---|---|---|---|
| Ce50 hemodinâmico = 3,0 µg/ml | [acoplamento.jl#L52](../src/acoplamento.jl#L52) | [grupo acoplamento.jl#L18](https://github.com/gabiadamek30-bit/Simula-o-da-resposta-do-paciente-ao-propofol---PBL/blob/b6995fc931af5a542bbc100f1677d93962f5b2b6/src/acoplamento.jl#L18) | Su 2022: EC50 da resistência periférica 2,96 µg/ml [10]; PAS (modelo Su): 3,21 em [PAS pd_models.py#L964](https://github.com/AnesthesiaSimulation/Python_Anesthesia_Simulator/blob/4afc03ad2997f725205c202e0d0634ff0664e5c7/src/python_anesthesia_simulator/pd_models.py#L964) | ILUSTRATIVO, na ordem da resistência periférica |
| γ = 2,0 | [acoplamento.jl#L55](../src/acoplamento.jl#L55) | [grupo acoplamento.jl#L19](https://github.com/gabiadamek30-bit/Simula-o-da-resposta-do-paciente-ao-propofol---PBL/blob/b6995fc931af5a542bbc100f1677d93962f5b2b6/src/acoplamento.jl#L19) | PAS: 1,83 em [PAS pd_models.py#L966](https://github.com/AnesthesiaSimulation/Python_Anesthesia_Simulator/blob/4afc03ad2997f725205c202e0d0634ff0664e5c7/src/python_anesthesia_simulator/pd_models.py#L966) | ILUSTRATIVO |
| Imax Rp = 0,25 | [acoplamento.jl#L58](../src/acoplamento.jl#L58) | [grupo acoplamento.jl#L20](https://github.com/gabiadamek30-bit/Simula-o-da-resposta-do-paciente-ao-propofol---PBL/blob/b6995fc931af5a542bbc100f1677d93962f5b2b6/src/acoplamento.jl#L20) | PAS: Emax −0,778 em [PAS pd_models.py#L965](https://github.com/AnesthesiaSimulation/Python_Anesthesia_Simulator/blob/4afc03ad2997f725205c202e0d0634ff0664e5c7/src/python_anesthesia_simulator/pd_models.py#L965), num modelo com retroalimentação (não comparável direto) | ILUSTRATIVO |
| Imax VS = 0,15 | [acoplamento.jl#L61](../src/acoplamento.jl#L61) | [grupo acoplamento.jl#L21](https://github.com/gabiadamek30-bit/Simula-o-da-resposta-do-paciente-ao-propofol---PBL/blob/b6995fc931af5a542bbc100f1677d93962f5b2b6/src/acoplamento.jl#L21) | PAS: Emax −0,154·e^(0,0333(idade − 35)) em [PAS pd_models.py#L968](https://github.com/AnesthesiaSimulation/Python_Anesthesia_Simulator/blob/4afc03ad2997f725205c202e0d0634ff0664e5c7/src/python_anesthesia_simulator/pd_models.py#L968), mas EC50 do VS de 0,34 µg/ml [10] (0,44 no PAS, [PAS pd_models.py#L967](https://github.com/AnesthesiaSimulation/Python_Anesthesia_Simulator/blob/4afc03ad2997f725205c202e0d0634ff0664e5c7/src/python_anesthesia_simulator/pd_models.py#L967)), não 3 | ILUSTRATIVO; o Ce50 único não serve para o VS |
| Imax FC = 0,05 | [acoplamento.jl#L64](../src/acoplamento.jl#L64) | [grupo acoplamento.jl#L22](https://github.com/gabiadamek30-bit/Simula-o-da-resposta-do-paciente-ao-propofol---PBL/blob/b6995fc931af5a542bbc100f1677d93962f5b2b6/src/acoplamento.jl#L22) | Su 2022: valores da FC não visíveis no resumo [10] | ILUSTRATIVO |
| piso de 5% do valor basal | [acoplamento.jl#L81](../src/acoplamento.jl#L81) | [grupo acoplamento.jl#L90](https://github.com/gabiadamek30-bit/Simula-o-da-resposta-do-paciente-ao-propofol---PBL/blob/b6995fc931af5a542bbc100f1677d93962f5b2b6/src/acoplamento.jl#L90) | — | ILUSTRATIVO (nunca é atingido com os Imax atuais) |
| s = 1 + 0,01(idade − 50)⁺ + 0,15(ASA − 2)⁺ + 0,03(PVI − 15)⁺ + 0,01(IVCI − 43)⁺ | [acoplamento.jl#L113-L114](../src/acoplamento.jl#L113-L114) | [grupo acoplamento.jl#L31-L32](https://github.com/gabiadamek30-bit/Simula-o-da-resposta-do-paciente-ao-propofol---PBL/blob/b6995fc931af5a542bbc100f1677d93962f5b2b6/src/acoplamento.jl#L31-L32) (PVI 13, IVCI 40) | limiares: PVI > 15 [12], IVCI ≥ 43% [13]; idade ≥ 50 e ASA III–V: Reich 2005 [14] (a confirmar) | limiares com fonte; **coeficientes ILUSTRATIVOS** |

O modelo de Su citado pelo PAS ([PAS pd_models.py#L918](https://github.com/AnesthesiaSimulation/Python_Anesthesia_Simulator/blob/4afc03ad2997f725205c202e0d0634ff0664e5c7/src/python_anesthesia_simulator/pd_models.py#L918)) é o de 2023,
de **interação com remifentanil** [11]; os EC50 só do propofol são do Su 2022 [10]. de Wit 2016 [15] mostra queda da
resistência arterial e do volume estressado com débito cardíaco preservado, um argumento a favor de Imax do VS pequeno.

### 4.5 Valores basais do paciente

| valor | nosso código | grupo | fonte (onde) | veredito |
|---|---|---|---|---|
| superfície corporal = √(altura × peso / 3600) m² (Mosteller) | [paciente.jl#L133](../src/paciente.jl#L133) | [grupo paciente.jl#L44](https://github.com/gabiadamek30-bit/Simula-o-da-resposta-do-paciente-ao-propofol---PBL/blob/b6995fc931af5a542bbc100f1677d93962f5b2b6/src/paciente.jl#L44) | Mosteller 1987 [16] | CONFIRMADO\* |
| IMC = peso/(altura/100)² | [paciente.jl#L109](../src/paciente.jl#L109) | [grupo paciente.jl#L30](https://github.com/gabiadamek30-bit/Simula-o-da-resposta-do-paciente-ao-propofol---PBL/blob/b6995fc931af5a542bbc100f1677d93962f5b2b6/src/paciente.jl#L30) | definição | CONFIRMADO |
| PAM do monitor (obrigatória) | [paciente.jl#L140](../src/paciente.jl#L140) | [grupo paciente.jl#L21](https://github.com/gabiadamek30-bit/Simula-o-da-resposta-do-paciente-ao-propofol---PBL/blob/b6995fc931af5a542bbc100f1677d93962f5b2b6/src/paciente.jl#L21) (entrada, padrão 90) | medida | — |
| DC basal = 3,0 L/min/m² × superfície corporal (se não medido) | [paciente.jl#L154](../src/paciente.jl#L154) | [grupo paciente.jl#L47](https://github.com/gabiadamek30-bit/Simula-o-da-resposta-do-paciente-ao-propofol---PBL/blob/b6995fc931af5a542bbc100f1677d93962f5b2b6/src/paciente.jl#L47) | índice cardíaco normal 2,5–4,0 L/min/m² [17] | CONFIRMADO\* para a faixa; o 3,0 é escolha do grupo |

**Por que o 3,0 não pesa no resultado:** o modelo multiplica Rp, VS e FC por fatores relativos, e a PAM média do
batimento é Rp·Q̄. A queda relativa da PAM, portanto, não depende do DC basal: o teste
[runtests.jl#L341](../test/runtests.jl#L341) varia o DC de 4 a 7 L/min e a queda
muda menos de 0,5 ponto percentual. O DC entra na complacência C calibrada e, por ela, só na forma da onda (PAS e PAD).

### 4.6 Limiares e regra de risco

| valor | nosso código | fonte (onde) | veredito |
|---|---|---|---|
| PAM < 55 mmHg | [metricas.jl#L7](../src/metricas.jl#L7) | "even short durations of an intraoperative MAP less than 55 mmHg are associated with AKI and myocardial injury" — Walsh 2013, Conclusões [18] | CONFIRMADO\* |
| PAM < 65 mmHg | [metricas.jl#L9](../src/metricas.jl#L9) | Salmasi 2017 [20]; consenso POQI: PAM abaixo de 60–70 mmHg associa-se a lesão miocárdica, renal e morte [21] | CONFIRMADO\* (faixa 60–70) |
| BIS-alvo < 60 (faixa 40–60) | [metricas.jl#L38](../src/metricas.jl#L38) | "target range of BIS values during general anaesthesia is 40 to 60" — NICE DG6, seção 4 [19] | CONFIRMADO\* |
| risco: ALTO < 55, MODERADO < 65, BAIXO ≥ 65 (PAM mínima) | [metricas.jl#L85](../src/metricas.jl#L85) | os dois limiares acima | regra nossa, só com limiares com fonte |

### 4.7 Protocolo de indução

| valor | nosso código | grupo | fonte | veredito |
|---|---|---|---|---|
| bolus 2 mg/kg em 30 s, manutenção 100 µg/kg/min a partir de 1 min | [propofol.jl#L141](../src/propofol.jl#L141) | [grupo propofol.jl#L42](https://github.com/gabiadamek30-bit/Simula-o-da-resposta-do-paciente-ao-propofol---PBL/blob/b6995fc931af5a542bbc100f1677d93962f5b2b6/src/propofol.jl#L42) | bula do Diprivan: indução 2–2,5 mg/kg no adulto saudável, 1–1,5 mg/kg em idoso ou ASA III–IV; manutenção 100–200 µg/kg/min [22] | CONFIRMADO\* |

## 5. Verificação numérica (testes automáticos)

Rodam no GitHub Actions a cada envio (`julia --project=. -e "using Pkg; Pkg.test()"`):

| teste | o que garante |
|---|---|
| [runtests.jl#L220](../test/runtests.jl#L220) | o integrador do NeoPKPD reproduz exp(K·t)·A₀ |
| [runtests.jl#L193](../test/runtests.jl#L193) | nossos parâmetros de Eleveld = os do PAS em 20 casos (5 pacientes de ambos os sexos × amostra arterial/venosa × com/sem opioide) |
| [runtests.jl#L238](../test/runtests.jl#L238) | E0, Ce50 e os dois γ da Tabela 3 |
| [runtests.jl#L250](../test/runtests.jl#L250) | o BIS é o do Ce de Eleveld deslocado pelo atraso |
| [runtests.jl#L124](../test/runtests.jl#L124) e [runtests.jl#L134](../test/runtests.jl#L134) | propriedades exatas do Windkessel |
| [runtests.jl#L146](../test/runtests.jl#L146) | a identificação recupera θ de dados sintéticos |
| [runtests.jl#L327](../test/runtests.jl#L327) | a regra de risco usa só os dois limiares |
| [runtests.jl#L341](../test/runtests.jl#L341) | a escolha do índice cardíaco 3,0 não muda a queda da PAM |
| [runtests.jl#L400](../test/runtests.jl#L400) | todo link deste documento para o nosso código aponta para uma linha existente e não vazia |

## 6. Pendências: o que precisamos validar juntos

### 6.1 Essencial para o modelo, com fonte a confirmar (abrir o artigo)

| item | o que conferir | onde |
|---|---|---|
| γ e atraso do BIS | Tabela 3: γ 1,89 (Ce < Ce50) e 1,47 (Ce > Ce50); θ do atraso 0,0517 | Eleveld 2018 [6] |
| PK de Schnider | equações e valores; e se a coluna de variabilidade é variância (ω²) ou coeficiente de variação | Schnider 1998 [1] |
| parâmetros do Windkessel | Tabela 1, batimento tipo A (a página do JSIM chama a Fig. 4A de "tipo C": conferir) | Stergiopulos 1999 [9] |
| limiares de PVI e IVCI | PVI > 15 e IVCI ≥ 43% | Tsuchiya 2010 [12]; Zhang & Critchley 2016 [13] |
| idade ≥ 50 e ASA III–V como fatores de risco | se o artigo traz esses cortes | Reich 2005 [14] (citação ainda não aberta nem por resumo) |

### 6.2 Essencial, sem fonte: precisa de dados de pacientes

- **Ce50, γ e Imax do efeito hemodinâmico** (seção 4.4) e **coeficientes de s**. São o coração da previsão de
  PAM e hoje são ilustrativos.
- **Faixa 0,2–6 ml/mmHg da complacência** na calibração.
- **PAS e PAD durante a queda:** a PAM segue Rp·Q̄, mas a forma da onda vem do batimento fixo de
  Stergiopulos. No paciente C, a PA no momento da PAM mínima sai 92/22 mmHg; uma PAD de 22 não é
  plausível. Só a PAM e a pressão de pulso basais estão calibradas; PAS e PAD previstas precisam de
  validação com linha arterial.
- **Critério da dose mínima da janela** (BIS < 60 em algum momento dos 20 min, inclusive pela manutenção): com o
  Eleveld completo, 0,5 mg/kg + 100 µg/kg/min chega a BIS < 60 depois de vários minutos. Decidir se a hipnose
  precisa ocorrer até um tempo máximo após o bolus.

### 6.3 Desativado até validar

| o que | por quê | onde ficou o código |
|---|---|---|
| faixa de 90% e "chance de PAM < 65/55" | variabilidade de Schnider com leitura incerta (ω² × CV, até 25× de diferença); ω = 0,3 de ke0, Ce50 e Imax assumidos | [populacao.jl#L29](../src/populacao.jl#L29), [incerteza.jl#L39](../src/incerteza.jl#L39) |
| coorte virtual | distribuições de idade, IMC, PA, ASA, PVI e IVCI inventadas por nós | [populacao.jl#L93](../src/populacao.jl#L93) |
| fórmula da PAM, "17–24%", "25–40%", queda > 30%, SpO₂, glicemia | sem fonte ou sem papel no modelo | retirados |

## 7. Como validar com dados: proposta

O banco aberto **VitalDB** [23] tem milhares de cirurgias com o registro da bomba de propofol (taxa de infusão,
Cp e Ce calculados), pressão arterial invasiva, BIS, FC e, em parte dos casos, débito cardíaco, além de idade,
sexo, peso, altura e ASA. Plano:

1. Selecionar induções **sem remifentanil** nos primeiros 20 min (o banco tem a trilha da bomba de
   remifentanil), com linha arterial antes da indução.
2. Para cada caso: PAS, PAD, PAM e FC do minuto antes da indução → paciente virtual; histórico de infusão →
   nosso PK.
3. Ajustar Ce50, γ, os três Imax e os coeficientes de s em 70% dos casos; testar nos outros 30%: erro da PAM
   mínima, RMSE da curva, AUROC para PAM < 65.
4. Comparar com os valores de Su 2022 [10] como referência externa.
5. Com os resíduos, estimar a variabilidade entre pacientes e só então reativar a faixa de incerteza.

O BIS pode ser validado no mesmo banco, comparando o BIS previsto pelo Eleveld com o medido.

## 8. Bibliografia

Nenhum artigo pôde ser aberto da rede usada na conferência (só o GitHub abria); o que foi conferido está na
coluna "conferido por".

| # | referência | link | origem | conferido por | usado para |
|---|---|---|---|---|---|
| 1 | Schnider TW et al. The influence of method of administration and covariates on the pharmacokinetics of propofol in adult volunteers. *Anesthesiology* 1998;88:1170–82 | [PMID 9605675](https://pubmed.ncbi.nlm.nih.gov/9605675/) | G | código PAS e tci | PK de Schnider |
| 2 | Schnider TW et al. The influence of age on propofol pharmacodynamics. *Anesthesiology* 1999;90:1502–16 | [DOI](https://doi.org/10.1097/00000542-199906000-00003) | G, P | código PAS; citação em PMID 23040653 | ke0 de Schnider |
| 3 | Marsh B et al. Pharmacokinetic model driven infusion of propofol in children. *Br J Anaesth* 1991;67:41–8 | [DOI](https://doi.org/10.1093/bja/67.1.41) | G, P | código PAS e tci | PK de Marsh |
| 4 | Absalom AR et al. Pharmacokinetic models for propofol — defining and illuminating the devil in the detail. *Br J Anaesth* 2009;103:26–37 | [PMID 19520702](https://pubmed.ncbi.nlm.nih.gov/19520702/) | B | resumo indexado | origem do ke0 0,26 |
| 5 | Struys MMRF et al. *Anesthesiology* 2000;92:399–406 | [DOI](https://doi.org/10.1097/00000542-200002000-00021) | B | citação em revisão | ke0 alternativo 1,21 |
| 6 | Eleveld DJ et al. Pharmacokinetic–pharmacodynamic model for propofol for broad application in anaesthesia and sedation. *Br J Anaesth* 2018;120:942–59 | [DOI](https://doi.org/10.1016/j.bja.2018.01.018) | G, P | código PAS (executado), resumo indexado da Tabela 3 | PK e BIS de Eleveld |
| 7 | Al-Sallami HS et al. *Clin Pharmacokinet* 2015 (fórmula de massa livre de gordura) | a buscar | P (nome da função) | código PAS; **citação completa ainda não conferida** | massa livre de gordura no Eleveld |
| 8 | James WPT. *Research on obesity*. London: HMSO, 1976 | — | G | código PAS e tci | massa magra |
| 9 | Stergiopulos N, Westerhof BE, Westerhof N. Total arterial inertance as the fourth element of the windkessel model. *Am J Physiol* 1999;276:H81–8 | [DOI](https://doi.org/10.1152/ajpheart.1999.276.1.H81) | W | código do windkessel-id | Windkessel de 4 elementos |
| 10 | Su H, Eleveld DJ, Struys MMRF, Colin PJ. *Br J Anaesth* 2022 (modelo hemodinâmico do propofol em voluntários) | [PMID 35249706](https://pubmed.ncbi.nlm.nih.gov/35249706/) | B | resumo indexado | comparação do efeito hemodinâmico |
| 11 | Su H et al. Pharmacodynamic mechanism-based interaction model for the haemodynamic effects of remifentanil and propofol in healthy volunteers. *Br J Anaesth* 2023;131:222–33 | [PMID 37355412](https://pubmed.ncbi.nlm.nih.gov/37355412/) | P | código PAS | modelo hemodinâmico do PAS |
| 12 | Tsuchiya M et al. *Acta Anaesthesiol Scand* 2010;54:596–602 (PVI e hipotensão na indução) | [PDF](https://www.masimo.fr/siteassets/us/documents/pdf/clinical-evidence/pvi/tsuchiya-pleth-variability-index-predicts-hypotension-may-2010.pdf) | B | resumo | limiar de PVI |
| 13 | Zhang J, Critchley LAH. *Anesthesiology* 2016;124:580–9 (ultrassom da veia cava antes da indução) | [PMID 26771910](https://pubmed.ncbi.nlm.nih.gov/26771910/) | B | resumo indexado | limiar de IVCI |
| 14 | Reich DL et al. *Anesth Analg* 2005;101:622–8 (preditores de hipotensão pós-indução) | a buscar | B | **ainda não conferido** | idade e ASA no fator s |
| 15 | de Wit F et al. The effect of propofol on haemodynamics: cardiac output, venous return, mean systemic filling pressure, and vascular resistances. *Br J Anaesth* 2016;116:784–9 | [PMID 27199311](https://pubmed.ncbi.nlm.nih.gov/27199311/) | B | resumo indexado | mecanismo da queda de PA |
| 16 | Mosteller RD. Simplified calculation of body-surface area. *N Engl J Med* 1987;317:1098 | [PMID 3657876](https://pubmed.ncbi.nlm.nih.gov/3657876/) | G | calculadoras e resumo | superfície corporal |
| 17 | StatPearls. Physiology, Cardiac Index | [NBK539905](https://www.ncbi.nlm.nih.gov/books/NBK539905/) | B | resumo indexado | faixa do índice cardíaco |
| 18 | Walsh M et al. Relationship between intraoperative mean arterial pressure and clinical outcomes after noncardiac surgery. *Anesthesiology* 2013;119:507–15 | [PMID 23835589](https://pubmed.ncbi.nlm.nih.gov/23835589/) | B | resumo indexado | limiar de 55 mmHg |
| 19 | NICE. Depth of anaesthesia monitors (DG6), seção 4 | [link](https://www.nice.org.uk/guidance/DG6/chapter/4-the-diagnostic-tests) | B | página indexada | BIS-alvo 40–60 |
| 20 | Salmasi V et al. Relationship between intraoperative hypotension, defined by either reduction from baseline or absolute thresholds, and acute kidney and myocardial injury after noncardiac surgery. *Anesthesiology* 2017;126:47–65 | [PMID 27792044](https://pubmed.ncbi.nlm.nih.gov/27792044/) | B | resumo indexado | limiar de 65 mmHg |
| 21 | Sessler DI et al. (POQI). Consenso sobre pressão arterial intraoperatória. *Br J Anaesth* 2019;122:563–74 | a buscar | B | resumo indexado | limiar de 60–70 mmHg |
| 22 | Diprivan (propofol), bula aprovada pela FDA, 2017 | [PDF](https://www.accessdata.fda.gov/drugsatfda_docs/label/2017/019627s066lbl.pdf) | B | trecho indexado | doses do protocolo |
| 23 | Lee HC et al. VitalDB, a high-fidelity multi-parameter vital signs database in surgical patients. *Sci Data* 2022 | [vitaldb.net](https://vitaldb.net) | B | **ainda não conferido** | proposta de validação |

**Software:** NeoPKPD (MIT), windkessel-id (MIT, Pigot et al.), Python Anesthesia Simulator (GPL-3, Aubouin-Pairault
et al.), [tci](https://github.com/jarretrt/tci) (R) e
[propofol-pebble-induction](https://github.com/ysuzuki1978/propofol-pebble-induction), estes dois usados só como
conferência das tabelas.

**Como reproduzir a conferência com o PAS:** os valores de referência do teste do Eleveld foram gerados
executando as linhas do Eleveld de `pk_models.py` do PAS (commit 4afc03a) num script Python com numpy, para os
casos listados em [runtests.jl#L70](../test/runtests.jl#L70).
