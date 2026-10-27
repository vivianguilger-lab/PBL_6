# =============================================================================
# dados.jl — dados de pressão e fluxo aórticos do windkessel-id
# =============================================================================
# Porta de `importPQ` (idanalysis.jl do windkessel-id) e da leitura dos dados de
# Stergiopulos (1999), como em src/dados.jl do PropofolPBL do grupo, sem
# depender de CSV.jl/DataFrames.jl (leitura direta dos arquivos de texto).

# Para quem não programa: este arquivo só LÊ dados de arquivos de texto (CSV);
# não há cálculo do modelo aqui. São dois conjuntos, copiados do windkessel-id
# (licença MIT; ver dados/windkessel-id/README.md e docs/validacao.md, seção 4.3):
#   1) Stergiopulos 1999 (humano): um batimento de pressão e fluxo aórticos.
#      O FLUXO é a forma de onda que alimenta o Windkessel de todo paciente
#      (acoplamento.jl a estica para durar 60/FC e a escala pelo volume sistólico).
#      Como essa forma é fixa, PAS e PAD previstas ainda precisam de validação
#      (docs/validacao.md, seção 6.2).
#   2) in vivo, porcos de 65 kg: pressão e fluxo aórticos medidos, usados só no
#      exemplo de identificação (scripts/06_identificacao_windkessel.jl).
#
# WINDKESSEL_ID_DIR: caminho da pasta dados/windkessel-id. `@__DIR__` é a pasta
# deste arquivo (src/) e joinpath monta o caminho no formato do sistema
# operacional. `const` indica um valor que não muda.
const WINDKESSEL_ID_DIR = joinpath(@__DIR__, "..", "dados", "windkessel-id")

# _numeros: transforma uma linha de texto como "0.005, 66.27" nos números
# [0.005, 66.27]. split corta nas vírgulas, strip tira os espaços e parse converte
# texto em número (Float64). O "." depois de cada função a aplica a cada pedaço
# (broadcasting).
_numeros(linha) = parse.(Float64, strip.(split(linha, ',')))

# dados_stergiopulos1999: lê o batimento tipo A de Stergiopulos et al. 1999
# (Fig. 4A do artigo), digitalizado pelo windkessel-id: 166 amostras a cada 5 ms
# (cerca de 0,83 s de batimento). Arquivos: apressure.csv (tempo em s, pressão em
# mmHg) e aflow.csv (tempo em s, fluxo em mL/s), idênticos aos do windkessel-id
# v1.0 (docs/validacao.md, seção 4.3).
# Saída: tempo (s), pressão (mmHg) e fluxo (mL/s), cada um como vetor.
"""
    dados_stergiopulos1999() -> (t, pressao_mmHg, fluxo_mLs)

Batimento tipo A de Stergiopulos et al. (1999), digitalizado, amostrado a 5 ms.
"""
function dados_stergiopulos1999()
    pasta = joinpath(WINDKESSEL_ID_DIR, "stergiopulos1999_data")
    # ler(f): função local de uma linha. Lê cada linha não vazia do arquivo f,
    # converte em números e empilha: reduce(hcat, ...) junta as linhas como colunas
    # de uma matriz e o apóstrofo (') transpõe, deixando uma amostra por linha
    # (coluna 1 = tempo, coluna 2 = valor).
    ler(f) = reduce(hcat, [_numeros(l) for l in eachline(joinpath(pasta, f)) if !isempty(strip(l))])'
    P = ler("apressure.csv")
    Q = ler("aflow.csv")
    # P[:, 1] = coluna do tempo; P[:, 2] = pressão; Q[:, 2] = fluxo.
    return P[:, 1], P[:, 2], Q[:, 2]
end

# importPQ: lê os registros in vivo do windkessel-id (porcos de 65 kg; a v1.0
# dizia 70 kg, corrigido na v2.0: docs/validacao.md, seção 3, item 12).
# Cada registro tem um arquivo de pressão (P, mmHg) e um de fluxo (Q, L/min), com
# horário "hh:mm:ss.sss" e valor em cada linha. ATENÇÃO: o fluxo está em L/min,
# não em mL/s; os θ identificados com esses dados ficam em "unidades L/min"
# (ver mLs2Lmin, no fim do arquivo).
# Saída: dois vetores (P e Q) com um NamedTuple por arquivo. NamedTuple é um
# pacote de valores com nome: P[1].t_s (tempos, s) e P[1].vb_va (valores).
"""
    importPQ(prefixo, nfiles) -> (P, Q)

Lê `nfiles` pares `<prefixo>P<i>.csv` / `<prefixo>Q<i>.csv` (colunas `timestamp`,
`vb.va`; linhas com `#` são comentários). Devolve vetores de NamedTuple
`(t_s, vb_va)`: tempo em s desde o início do arquivo, pressão em mmHg ou fluxo em L/min.
Ex.: `importPQ(joinpath(WINDKESSEL_ID_DIR, "invivo_data", "Ao_baseline1_Ao"), 4)`.
"""
function importPQ(prefixo::AbstractString, nfiles::Int)
    # ler(f): lê um arquivo e devolve (tempos em s, valores).
    function ler(f)
        t = Float64[]; v = Float64[]
        for l in eachline(f)
            # Pula comentários (#), o cabeçalho "timestamp" e linhas vazias.
            # `cond && continue` = "se cond for verdadeira, vá para a próxima linha".
            (startswith(l, "#") || startswith(l, "timestamp") || isempty(strip(l))) && continue
            # Separa horário e valor; depois, horas, minutos e segundos.
            ts, vs = split(l, ',')
            h, m, s = split(ts, ':')
            # Converte o horário em segundos desde a meia-noite e guarda o valor.
            push!(t, 3600 * parse(Float64, h) + 60 * parse(Float64, m) + parse(Float64, s))
            push!(v, parse(Float64, vs))
        end
        # t .- t[1]: subtrai o primeiro tempo, para o registro começar em 0 s.
        return (t_s = t .- t[1], vb_va = v)
    end
    # Lê os arquivos <prefixo>P1.csv ... <prefixo>Pn.csv e os Q correspondentes.
    # "P$(i).csv" insere o número i no texto; "*" junta textos em Julia.
    P = [ler(prefixo * "P$(i).csv") for i in 1:nfiles]
    Q = [ler(prefixo * "Q$(i).csv") for i in 1:nfiles]
    # Pressão e fluxo do mesmo registro devem ter o mesmo número de amostras;
    # se não tiverem, o programa só avisa (@warn) e continua.
    for i in 1:nfiles
        length(P[i].vb_va) == length(Q[i].vb_va) || @warn "P$(i) e Q$(i) com tamanhos diferentes."
    end
    return P, Q
end

# mLs2Lmin: converte θ das unidades com fluxo em mL/s (as do modelo) para fluxo em
# L/min, mantendo pressão em mmHg e tempo em s. Útil para comparar com θ
# identificado nos dados in vivo (fluxo em L/min). Fatores: Rp, Rc e L × 1000/60;
# C × 60/1000. `[1:length(θ)]` usa só os fatores necessários, então serve para θ
# de 2, 3 ou 4 elementos.
"Converte θ = [Rp, C, Rc, L] de unidades mL/s para unidades L/min (pressão em mmHg)."
mLs2Lmin(θ) = θ .* [1000 / 60, 60 / 1000, 1000 / 60, 1000 / 60][1:length(θ)]
