# Manual PBL-6

## Links e requisitos

| | |
|---|---|
| Repositório | [github.com/vivianguilger-lab/PBL_6](https://github.com/vivianguilger-lab/PBL_6) |
| Download direto (ZIP) | [PBL_6 — ZIP](https://github.com/vivianguilger-lab/PBL_6/archive/refs/heads/main.zip) |
| Manual no repositório | [MANUAL.md](https://github.com/vivianguilger-lab/PBL_6/blob/main/MANUAL.md) |
| O que precisa | Windows 10/11, macOS ou Linux; internet na primeira vez; cerca de 3 GB livres |

## 1. Instalar (só na primeira vez)

Escolha um dos dois caminhos. Nos dois, a primeira vez leva 10 a 20 minutos, porque o Julia
baixa e compila os pacotes.

### Caminho A — duplo clique (recomendado)

**Windows**

1. Baixe o ZIP pelo link acima.
2. Na pasta Downloads, clique com o botão direito no ZIP → **Extrair tudo…** → **Extrair**.
   Não rode nada de dentro do ZIP sem extrair.
3. Abra a pasta extraída e dê **duplo clique em `abrir_plataforma.bat`**.
4. Se aparecer "O Windows protegeu o computador", clique em **Mais informações → Executar
   assim mesmo**.
5. Se o Julia não estiver instalado, o arquivo instala sozinho pela Microsoft Store (winget).
   Aceite as permissões.
6. Espere. A janela preta mostra o progresso; no fim, o navegador abre a plataforma sozinho.

**macOS**

1. Baixe o ZIP e dê duplo clique para extrair.
2. Na pasta extraída, clique com o botão direito em **`abrir_plataforma.command`** → **Abrir**
   → **Abrir** (só na primeira vez, por causa do Gatekeeper).
3. O Julia é instalado automaticamente se faltar. Espere a compilação.

**Linux**: no terminal, dentro da pasta extraída, `./abrir_plataforma.command`.

### Caminho B — pelo próprio Julia (eu fiz assim)

1. Instale o Julia 1.10 ou mais novo. Windows, no PowerShell:
   `winget install --name Julia --id 9NJNWW8PVKMN -e -s msstore`. macOS/Linux, no terminal:
   `curl -fsSL https://install.julialang.org | sh`.
2. Abra o **Julia** (menu Iniciar → Julia). Aparece o prompt verde `julia>`.
3. Cole, uma linha por vez, e espere terminar:

```julia
using Pkg
Pkg.add(url="https://github.com/vivianguilger-lab/PBL_6")
```

## 2. Abrir a plataforma (todas as vezes)

- **Caminho A:** duplo clique em `abrir_plataforma.bat` (Windows) ou
  `abrir_plataforma.command` (macOS/Linux). Abre em menos de 1 minuto.
- **Caminho B (eu fiz assim):** abra o Julia e digite:

```julia
using PropofolWindkessel
plataforma()
```

O navegador abre em **http://localhost:8080**. Se a porta estiver ocupada, a janela do Julia
mostra o endereço usado (por exemplo 8081); se o navegador não abrir, copie esse endereço na
barra do navegador.

**Deixe a janela preta / do Julia aberta enquanto usa.** Para encerrar, feche essa janela ou
aperte Ctrl+C nela. Os dados digitados não saem do computador: o servidor só aceita conexões
da própria máquina.
