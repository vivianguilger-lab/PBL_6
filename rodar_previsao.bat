@echo off
rem ==========================================================================
rem  Duplo clique: roda a previsao do paciente de exemplo.
rem  Outro paciente: arraste o arquivo .toml (pasta dados) para cima deste .bat.
rem ==========================================================================
cd /d "%~dp0"

where julia >nul 2>nul
if errorlevel 1 (
    echo Julia nao encontrado. Instale pelo juliaup e abra este arquivo de novo:
    echo   winget install --name Julia --id 9NJNWW8PVKMN -e -s msstore
    pause
    exit /b 1
)

set "ARQ=%~1"
if "%ARQ%"=="" set "ARQ=dados\paciente_exemplo.toml"

echo.
echo Pasta do projeto: %CD%
echo Verificando pacotes (na primeira vez demora de 5 a 15 minutos)...
julia --project=. -e "using Pkg; Pkg.instantiate()"
if errorlevel 1 goto erro

echo.
echo Rodando a previsao para %ARQ% ...
julia --project=. -t auto scripts\02_prever_paciente.jl "%ARQ%"
if errorlevel 1 goto erro

echo.
echo Pronto. As figuras estao na pasta resultados (vou abri-la agora).
start "" "%~dp0resultados"
pause
exit /b 0

:erro
echo.
echo Algo deu errado. Tire um print desta janela e envie.
pause
exit /b 1
