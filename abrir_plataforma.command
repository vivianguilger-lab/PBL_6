#!/usr/bin/env bash
# Duplo clique (macOS) ou ./abrir_plataforma.command (Linux): instala o Julia
# se faltar, prepara o projeto e abre a plataforma PBL-6 no navegador.
#
# Para rodar no Linux pelo terminal: entre na pasta do projeto e digite
#   bash abrir_plataforma.command
# No macOS, se aparecer aviso de segurança: botão direito no arquivo > Abrir.
# Na primeira vez, o Julia baixa e compila os pacotes (10 a 20 min).
#
# Vai para a pasta onde está este arquivo (a pasta do projeto). "$0" é o
# caminho do próprio script. Se não conseguir, encerra (exit 1).
cd "$(dirname "$0")" || exit 1

# Procura o Julia: `command -v julia` diz se o comando julia existe.
# Se não existir e também não estiver em ~/.juliaup/bin, instala com o
# juliaup (instalador oficial; precisa de internet). --yes aceita o padrão.
# Se a instalação falhar, mostra o site oficial e encerra.
if ! command -v julia >/dev/null 2>&1; then
    if [ ! -x "$HOME/.juliaup/bin/julia" ]; then
        echo "Julia não encontrado. Instalando pelo juliaup..."
        curl -fsSL https://install.julialang.org | sh -s -- --yes || {
            echo "Não consegui instalar o Julia. Veja https://julialang.org/downloads/"
            read -r -p "Enter para sair"; exit 1; }
    fi
    # PATH = lista de pastas onde o sistema procura programas; aqui a pasta
    # do juliaup é acrescentada só para esta janela
    export PATH="$HOME/.juliaup/bin:$PATH"
fi

# Avisos para o usuário
echo
echo "PBL-6 — preparando a plataforma"
echo "Na primeira vez o Julia baixa e compila os pacotes (10 a 20 min)."
echo "O navegador abre sozinho. NÃO feche esta janela enquanto usa."
echo
# Inicia o Julia:
#   --project=.        usa os pacotes de Project.toml e Manifest.toml;
#   -t auto            usa todos os núcleos do processador;
#   Pkg.instantiate()  baixa e instala os pacotes que faltam;
#   plataforma()       inicia o servidor local e abre http://localhost:8080;
#   wait(...)          mantém tudo rodando até fechar a janela (ou Ctrl+C).
# Os dados digitados não saem do computador: o servidor só aceita conexões locais.
julia --project=. -t auto -e 'using Pkg; Pkg.instantiate(); using PropofolWindkessel; wait(plataforma())'
# quando a plataforma termina, espera Enter antes de fechar a janela
read -r -p "Plataforma encerrada. Enter para sair"
