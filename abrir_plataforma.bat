@echo off
rem Duplo clique: abre a plataforma PropoSim no navegador.
cd /d "%~dp0"
where julia >nul 2>nul
if errorlevel 1 (
    echo Julia nao encontrado. Instale pelo juliaup:
    echo   winget install --name Julia --id 9NJNWW8PVKMN -e -s msstore
    pause
    exit /b 1
)
echo Preparando a plataforma (na primeira vez demora alguns minutos)...
julia --project=. -t auto -e "using Pkg; Pkg.instantiate(); using PropofolWindkessel; wait(plataforma())"
pause
