@echo off
rem ==========================================================================
rem  Duplo clique: roda a previsao do paciente de exemplo.
rem  Outro paciente: arraste o arquivo .toml (pasta dados) para cima deste .bat.
rem ==========================================================================
rem Linhas que comecam com rem sao comentarios e nao sao executadas.
rem O que faz: roda scripts\02_prever_paciente.jl para um paciente e abre a
rem pasta resultados com as figuras e o CSV da previsao.
rem Precisa do Julia instalado. Se nao tiver, rode antes o
rem abrir_plataforma.bat, que instala o Julia sozinho.
rem cd /d vai para a pasta onde este arquivo esta, a pasta do projeto.
cd /d "%~dp0"

rem Confere se o comando julia existe. Se nao existir, errorlevel fica 1 e o
rem bloco abaixo mostra como instalar e encerra.
where julia >nul 2>nul
if errorlevel 1 (
    echo Julia nao encontrado. Instale pelo juliaup e abra este arquivo de novo:
    echo   winget install --name Julia --id 9NJNWW8PVKMN -e -s msstore
    pause
    exit /b 1
)

rem ARQ = arquivo .toml arrastado para cima deste .bat, o primeiro argumento.
rem Sem arquivo, usa o paciente de exemplo, dados\paciente_exemplo.toml.
set "ARQ=%~1"
if "%ARQ%"=="" set "ARQ=dados\paciente_exemplo.toml"

rem Instala ou confere os pacotes do projeto com Pkg.instantiate.
rem Se der erro, vai para a mensagem de erro no fim.
echo.
echo Pasta do projeto: %CD%
echo Verificando pacotes (na primeira vez demora de 5 a 15 minutos)...
julia --project=. -e "using Pkg; Pkg.instantiate()"
if errorlevel 1 goto erro

rem Roda a previsao com todos os nucleos do processador, opcao -t auto.
rem Saidas em resultados: figuras .png e a tabela batimento a batimento .csv.
echo.
echo Rodando a previsao para %ARQ% ...
julia --project=. -t auto scripts\02_prever_paciente.jl "%ARQ%"
if errorlevel 1 goto erro

rem Tudo certo: abre a pasta resultados no Explorador de Arquivos.
echo.
echo Pronto. As figuras estao na pasta resultados (vou abri-la agora).
start "" "%~dp0resultados"
pause
exit /b 0

rem Mensagem de erro: o print da janela ajuda o grupo a achar o problema.
:erro
echo.
echo Algo deu errado. Tire um print desta janela e envie.
pause
exit /b 1
