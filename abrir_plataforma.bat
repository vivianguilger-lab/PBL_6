@echo off
setlocal
title PBL-6 - simulador da inducao com propofol
rem Duplo clique: instala o Julia se faltar, prepara o projeto e abre a
rem plataforma PBL-6 no navegador. Deixe esta janela aberta durante o uso.
cd /d "%~dp0"

if not exist "Project.toml" goto :sem_projeto

set "JULIA=julia"
where julia >nul 2>nul
if not errorlevel 1 goto :rodar
set "JULIA=%LOCALAPPDATA%\Microsoft\WindowsApps\julia.exe"
if exist "%JULIA%" goto :rodar

echo.
echo Julia nao encontrado. Instalando pela Microsoft Store via winget...
echo Aceite as permissoes que o Windows pedir.
echo.
winget install --name Julia --id 9NJNWW8PVKMN -e -s msstore --accept-package-agreements --accept-source-agreements
if exist "%JULIA%" goto :rodar
echo.
echo Nao consegui instalar o Julia automaticamente.
echo Instale pelo site https://julialang.org/downloads/ e rode este arquivo de novo.
pause
exit /b 1

:rodar
echo.
echo ================================================================
echo  PBL-6 - preparando a plataforma
echo  Na primeira vez o Julia baixa e compila os pacotes: 10 a 20 min.
echo  Depois disso abre em menos de 1 minuto.
echo  O navegador abre sozinho. NAO feche esta janela enquanto usa.
echo ================================================================
echo.
"%JULIA%" --project=. -t auto -e "using Pkg; Pkg.instantiate(); using PropofolWindkessel; wait(plataforma())"
echo.
echo A plataforma foi encerrada.
pause
exit /b 0

:sem_projeto
echo.
echo Este arquivo precisa estar dentro da pasta do projeto.
echo Se voce abriu o ZIP, extraia primeiro: botao direito no ZIP, Extrair tudo,
echo e rode o abrir_plataforma.bat de dentro da pasta extraida.
pause
exit /b 1
