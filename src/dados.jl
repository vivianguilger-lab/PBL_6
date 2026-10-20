# =============================================================================
# dados.jl — dados de pressão e fluxo aórticos do windkessel-id
# =============================================================================
# Porta de `importPQ` (idanalysis.jl do windkessel-id) e da leitura dos dados de
# Stergiopulos (1999), como em src/dados.jl do PropofolPBL do grupo, sem
# depender de CSV.jl/DataFrames.jl (leitura direta dos arquivos de texto).

const WINDKESSEL_ID_DIR = joinpath(@__DIR__, "..", "dados", "windkessel-id")

_numeros(linha) = parse.(Float64, strip.(split(linha, ',')))

"""
    dados_stergiopulos1999() -> (t, pressao_mmHg, fluxo_mLs)

Batimento tipo A de Stergiopulos et al. (1999), digitalizado, amostrado a 5 ms.
"""
function dados_stergiopulos1999()
    pasta = joinpath(WINDKESSEL_ID_DIR, "stergiopulos1999_data")
    ler(f) = reduce(hcat, [_numeros(l) for l in eachline(joinpath(pasta, f)) if !isempty(strip(l))])'
    P = ler("apressure.csv")
    Q = ler("aflow.csv")
    return P[:, 1], P[:, 2], Q[:, 2]
end

"""
    importPQ(prefixo, nfiles) -> (P, Q)

Lê `nfiles` pares `<prefixo>P<i>.csv` / `<prefixo>Q<i>.csv` (colunas `timestamp`,
`vb.va`; linhas com `#` são comentários). Devolve vetores de NamedTuple
`(t_s, vb_va)`: tempo em s desde o início do arquivo, pressão em mmHg ou fluxo em L/min.
Ex.: `importPQ(joinpath(WINDKESSEL_ID_DIR, "invivo_data", "Ao_baseline1_Ao"), 4)`.
"""
function importPQ(prefixo::AbstractString, nfiles::Int)
    function ler(f)
        t = Float64[]; v = Float64[]
        for l in eachline(f)
            (startswith(l, "#") || startswith(l, "timestamp") || isempty(strip(l))) && continue
            ts, vs = split(l, ',')
            h, m, s = split(ts, ':')
            push!(t, 3600 * parse(Float64, h) + 60 * parse(Float64, m) + parse(Float64, s))
            push!(v, parse(Float64, vs))
        end
        return (t_s = t .- t[1], vb_va = v)
    end
    P = [ler(prefixo * "P$(i).csv") for i in 1:nfiles]
    Q = [ler(prefixo * "Q$(i).csv") for i in 1:nfiles]
    for i in 1:nfiles
        length(P[i].vb_va) == length(Q[i].vb_va) || @warn "P$(i) e Q$(i) com tamanhos diferentes."
    end
    return P, Q
end

"Converte θ = [Rp, C, Rc, L] de unidades mL/s para unidades L/min (pressão em mmHg)."
mLs2Lmin(θ) = θ .* [1000 / 60, 60 / 1000, 1000 / 60, 1000 / 60][1:length(θ)]
