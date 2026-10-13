# extras — o circuito com CirculatorySystemModels.jl

`windkessel_csm.jl` monta o mesmo circuito do paciente virtual (ventrículo de
elastância variável de Shi, valvas como diodos, Windkessel de 3 elementos e
compartimento venoso com volume não estressado) usando a biblioteca
[CirculatorySystemModels.jl](https://github.com/TS-CUBED/CirculatorySystemModels.jl),
construída sobre o [ModelingToolkit.jl](https://github.com/SciML/ModelingToolkit.jl).

Em vez de escrever as EDOs, ligam-se componentes (`connect`) como num simulador
de circuitos elétricos — a analogia eletromecânica do curso levada ao pé da
letra. O ModelingToolkit aplica Kirchhoff e reduz o sistema simbolicamente.

O script compara o último batimento das duas implementações (biblioteca × código
à mão do pacote `PropofolWindkessel`): é uma **verificação cruzada** entre duas
implementações independentes das mesmas equações.

Fica num ambiente separado (Julia ≥ 1.11) porque o ModelingToolkit é grande (a primeira
instalação/compilação leva dezenas de minutos):

```
julia --project=extras -e 'using Pkg; Pkg.instantiate()'
julia --project=extras extras/windkessel_csm.jl
```

## Resultado (paciente virtual de referência, último de 30 batimentos)

| grandeza | CirculatorySystemModels.jl | código à mão | diferença |
|---|---|---|---|
| PAS aórtica (mmHg) | 111,39 | 111,39 | 1 × 10⁻⁴ |
| PAD aórtica (mmHg) | 72,76 | 72,76 | 2 × 10⁻⁴ |
| PAM (mmHg) | 91,53 | 91,53 | 8 × 10⁻⁵ |
| volume sistólico (ml) | 85,20 | 85,21 | 7 × 10⁻⁴ |

Testado com ModelingToolkit 11.45 e CirculatorySystemModels 0.6.0 (Julia 1.13).
Detalhe prático: a valva aórtica e a impedância característica foram unidas num
único `ResistorDiode(R = Rav + Zc)` — com dois elementos em série o nó entre eles
vira uma variável algébrica e o sistema passa a ser uma DAE (o `Tsit5` recusa
matriz de massa). A pressão na raiz da aorta é recuperada como `Pc + Zc·Q`.
