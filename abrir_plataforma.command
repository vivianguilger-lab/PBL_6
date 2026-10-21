#!/usr/bin/env bash
# Duplo clique (macOS) ou ./abrir_plataforma.command (Linux): instala o Julia
# se faltar, prepara o projeto e abre a plataforma PropoSim no navegador.
cd "$(dirname "$0")" || exit 1

if ! command -v julia >/dev/null 2>&1; then
    if [ ! -x "$HOME/.juliaup/bin/julia" ]; then
        echo "Julia não encontrado. Instalando pelo juliaup..."
        curl -fsSL https://install.julialang.org | sh -s -- --yes || {
            echo "Não consegui instalar o Julia. Veja https://julialang.org/downloads/"
            read -r -p "Enter para sair"; exit 1; }
    fi
    export PATH="$HOME/.juliaup/bin:$PATH"
fi

echo
echo "PropoSim — preparando a plataforma"
echo "Na primeira vez o Julia baixa e compila os pacotes (10 a 20 min)."
echo "O navegador abre sozinho. NÃO feche esta janela enquanto usa."
echo
julia --project=. -t auto -e 'using Pkg; Pkg.instantiate(); using PropofolWindkessel; wait(plataforma())'
read -r -p "Plataforma encerrada. Enter para sair"
