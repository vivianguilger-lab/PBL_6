# Revisão: repositórios no GitHub e literatura usada

Busca feita no GitHub por `windkessel`, `windkessel language:julia`,
`lumped parameter cardiovascular model julia`, `propofol pharmacokinetic`,
`Schnider propofol`, `anesthesia simulator propofol` e `propofol hemodynamic`
(outubro/2026). Os dois repositórios marcados com ★ foram lidos em detalhe e
são a base direta deste projeto.

## 0. Base do modelo: PropofolPBL do grupo

O núcleo do software segue o repositório do grupo
[gabiadamek30-bit/Simula-o-da-resposta-do-paciente-ao-propofol---PBL](https://github.com/gabiadamek30-bit/Simula-o-da-resposta-do-paciente-ao-propofol---PBL)
(pacote Julia `PropofolPBL`), que junta:

| repositório | o que o modelo usa |
|---|---|
| ★ [shramish2057/NeoPKPD](https://github.com/shramish2057/NeoPKPD) (pacote registrado no Julia, MIT) | motor PK/PD: `ThreeCompIVBolus` com bolus e infusões, `simulate_biophase` (ke0) |
| ★ [hpigot/windkessel-id](https://github.com/hpigot/windkessel-id) (MIT) | Windkessel de 2/3/4 elementos (`G2`, `G3`, `G`), `solvex0`, identificação por Newton (`identify`), dados de Stergiopulos (1999) e in vivo |

Do grupo vieram: o `Paciente` com ASA, PVI, IVCI, SpO₂ e glicemia; a PK de Schnider/Marsh no
NeoPKPD; o BIS de Eleveld; a porta do windkessel-id para Julia atual (`src/windkessel.jl`); o
acoplamento Emax sigmoide Ce → Rp, VS, FC com o fator de sensibilidade; o fluxo aórtico de
Stergiopulos escalado por VS e FC; e os testes. Os dados do windkessel-id estão em
`dados/windkessel-id/` com a licença original.

Os repositórios abaixo (seções 1 e 2) foram a pesquisa inicial; o
`CirculatorySystemModels.jl` e o modelo de malha fechada que ele inspirou foram usados na
versão anterior do software (até o commit `d3dd331`).

## 1. Repositórios de modelos Windkessel / circulação

| repositório | linguagem | o que oferece | como usamos |
|---|---|---|---|
| ★ [TS-CUBED/CirculatorySystemModels.jl](https://github.com/TS-CUBED/CirculatorySystemModels.jl) | **Julia** (ModelingToolkit) | Biblioteca acausal de componentes 0-D: `Resistor`, `Compliance` (com volume não estressado V₀), `Elastance`, `WK3`, `WK4`, `WK5`, câmaras de elastância variável (`DHChamber` duplo-Hill, `ShiChamber`), valvas (`ResistorDiode`, `ShiValve`, `MynardValve`). Exemplo completo do modelo de Bjørdalsbakke et al. (2022): um ventrículo + Windkessel + compartimento venoso em malha fechada. | Usado na versão anterior (malha fechada com elastância variável). O modelo atual segue o PropofolPBL do grupo (seção 0). |
| [dcanuto/WK3](https://github.com/dcanuto/WK3) | Julia | Windkessel de 3 elementos em malha aberta com fonte de pressão ventricular. | Referência para o WK3. |
| [hpigot/windkessel-id](https://github.com/hpigot/windkessel-id) | Julia | Identificação de Windkessel de 2, 3 e 4 elementos a partir de dados de pressão e fluxo. | **Usado no modelo atual** (via PropofolPBL do grupo): Windkessel de 4 elementos e identificação. |
| [lawmurray/Windkessel](https://github.com/lawmurray/Windkessel) | LibBi/Matlab | WK3 com estimação bayesiana. | Ideia de quantificar incerteza dos parâmetros. |
| [MiSpitieris/BC-of-Arterial-WK](https://github.com/MiSpitieris/BC-of-Arterial-WK) | R | Código do artigo "Bayesian Calibration of Arterial Windkessel Model". | Idem — calibração com incerteza. |
| [TAMU-ESP/WPINN](https://github.com/TAMU-ESP/WPINN) | Python | Rede neural informada pela física do Windkessel. | Mostra a fronteira *grey-box*/aprendizado de máquina. |
| [VanessaTsai0828/Windkessel-BloodPressure-Simulator](https://github.com/VanessaTsai0828/Windkessel-BloodPressure-Simulator) | Python | Ondas de pressão com Windkessel de 2, 3 e 4 elementos. | Material didático de comparação. |
| [ahmedfawaz879/virtual-patient-population](https://github.com/ahmedfawaz879/virtual-patient-population) | Python | Ensaio *in silico* no estilo ASME V&V 40 com população virtual Windkessel, verificado contra o regime permanente analítico. | Mesma filosofia de verificação usada aqui (testes contra soluções analíticas). |
| [TS-CUBED/haemoFoam](https://github.com/TS-CUBED/haemoFoam), [EManchester/OpenFOAM-v8-Windkessel-code](https://github.com/EManchester/OpenFOAM-v8-Windkessel-code) | C++/OpenFOAM | Windkessel como **condição de contorno** de CFD 3-D. | Ilustra o uso 0-D como contorno de modelos 3-D (Resumo, Seção 11). |

## 2. Repositórios de anestesia / PK-PD do propofol

| repositório | linguagem | o que oferece | como usamos |
|---|---|---|---|
| ★ [AnesthesiaSimulation/Python_Anesthesia_Simulator](https://github.com/AnesthesiaSimulation/Python_Anesthesia_Simulator) (PAS) | Python | Simulador de anestesia venosa total: PK de Schnider, Marsh, Schüttler e Eleveld; BIS (Bouillon, Vanluchene, Eleveld); modelo hemodinâmico mecanístico de Su et al. (TPR, VS, FC com retroalimentação; o PAS cita o artigo de 2023, de interação com remifentanil; os EC50 só do propofol são do Su 2022); variabilidade interindividual. | Conferência independente das equações de **Schnider**, **Marsh** e **Eleveld** (o nosso Eleveld é testado contra o código do PAS); ordem de grandeza do efeito hemodinâmico. A leitura das variâncias de Schnider no PAS difere da nossa (ver [validacao.md](validacao.md)). |
| [victorbezerra-dev/simethesia-simulator](https://github.com/victorbezerra-dev/simethesia-simulator) | C++/Arduino | Simulador PK/PD do BIS em tempo real. | — |
| [danieltorres-sci/propofol-simulator](https://github.com/danieltorres-sci/propofol-simulator) | Python | Modelo de Schnider com concentração no local de efeito. | Conferência dos valores de Schnider. |
| [shashankag14/Propofol-Predictor](https://github.com/shashankag14/Propofol-Predictor) | Python | Predição de concentração (PK + aprendizado de máquina). | — |

**Lacuna encontrada:** além do PropofolPBL do grupo, nenhum repositório acopla PK/PD do
propofol a um Windkessel calibrado por paciente para prever a curva de PAM antes da injeção.

## 3. Por que Julia

- `DifferentialEquations.jl`/`OrdinaryDiffEq.jl`: integradores de alta ordem com passo
  adaptativo, `tstops` para as descontinuidades da infusão — 20 min de simulação pulsátil
  (≈1 500 batimentos) em menos de 1 s, o que torna viável calibrar e simular **centenas de pacientes
  virtuais** em poucos minutos.
- `NeoPKPD.jl`: motor de farmacometria em Julia (PK/PD, populações, estimação).
- `Optim.jl` + `ForwardDiff.jl`: identificação do Windkessel por Newton com derivadas automáticas.
- Velocidade de C com sintaxe próxima da matemática (dá para escrever `τ = TPR·AC`).
- Recomendação: **fazer direto em Julia** — a sintaxe é próxima de Python/Matlab e o
  ecossistema de EDOs é o ponto forte da linguagem; um protótipo em Python só duplicaria o
  trabalho.

## 4. Literatura

A bibliografia conferida, com origem de cada citação e o que ela sustenta no modelo, está em
[validacao.md](validacao.md#8-bibliografia). A lista abaixo é a da pesquisa inicial.

**Farmacocinética e farmacodinâmica**
- Schnider TW et al. The influence of method of administration and covariates on the pharmacokinetics of propofol in adult volunteers. *Anesthesiology* 1998;88:1170–82.
- Schnider TW et al. The influence of age on propofol pharmacodynamics. *Anesthesiology* 1999;90:1502–16.
- Kazama T et al. Comparison of the effect-site k(eO)s of propofol for blood pressure and EEG bispectral index in elderly and younger patients. *Anesthesiology* 1999;90:1517–27.
- Bouillon TW et al. Pharmacodynamic interaction between propofol and remifentanil regarding hypnosis, tolerance of laryngoscopy, bispectral index, and electroencephalographic approximate entropy. *Anesthesiology* 2004;100:1353–72.
- Su H, Koomen JV, Eleveld DJ, Struys MMRF, Colin PJ. Pharmacodynamic mechanism-based interaction model for the haemodynamic effects of remifentanil and propofol in healthy volunteers. *Br J Anaesth* 2023;131(2):222–233 (interação com remifentanil).
- Su H, Eleveld DJ, Struys MMRF, Colin PJ. *Br J Anaesth* 2022 (efeitos hemodinâmicos só do propofol; PMID 35249706).
- Eleveld DJ et al. Pharmacokinetic–pharmacodynamic model for propofol for broad application in anaesthesia and sedation. *Br J Anaesth* 2018;120:942–59.
- Marsh B et al. Pharmacokinetic model driven infusion of propofol in children. *Br J Anaesth* 1991;67:41–8.

**Mecanismos hemodinâmicos do propofol**
- de Wit F et al. The effect of propofol on haemodynamics: cardiac output, venous return, mean systemic filling pressure, and vascular resistances. *Br J Anaesth* 2016;116(6):784–789.
- Cullen PM et al. Effect of propofol anesthesia on baroreflex activity in humans. *Anesth Analg* 1987;66:1115–20.
- Ebert TJ et al. Sympathetic and hemodynamic effects of moderate and deep sedation with propofol in humans. *Anesthesiology* 1992;76:725–33.
- Claeys MA, Gepts E, Camu F. Haemodynamic changes during anaesthesia induced and maintained with propofol. *Br J Anaesth* 1988;60:3–9.

**Circulação e controle**
- Stergiopulos N, Westerhof BE, Westerhof N. Total arterial inertance as the fourth element of the windkessel model. *Am J Physiol* 1999;276:H81–8.
- Bjørdalsbakke NL et al. Parameter estimation for closed-loop lumped parameter models of the systemic circulation using synthetic data. *Math Biosci* 2022;343:108731.
- Shi Y, Lawford P, Hose R. Review of zero-D and 1-D models of blood flow in the cardiovascular system. *Biomed Eng Online* 2011;10:33.
- Ursino M. Interaction between carotid baroregulation and the pulsating heart: a mathematical model. *Am J Physiol* 1998;275:H1733–47.
- McEniery CM et al. Normal vascular aging: differential effects on wave reflection and aortic pulse wave velocity (Anglo-Cardiff Collaborative Trial). *J Am Coll Cardiol* 2005;46:1753–60.

**Antropometria**
- Nadler SB, Hidalgo JU, Bloch T. Prediction of blood volume in normal human adults. *Surgery* 1962;51:224–32.
- James WPT. *Research on obesity*. London: HMSO, 1976 (fórmula de massa magra).
- Mosteller RD. Simplified calculation of body-surface area. *N Engl J Med* 1987;317:1098.

**Software**
- Schenkel T, Saxton H et al. CirculatorySystemModels.jl (Julia/ModelingToolkit).
- Aubouin-Pairault B et al. Python Anesthesia Simulator (PAS).
