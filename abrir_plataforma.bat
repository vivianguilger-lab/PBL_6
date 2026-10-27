@echo off
rem Linhas que comecam com rem sao comentarios e nao sao executadas.
rem A linha @echo off acima esconde os comandos na janela.
rem setlocal: as variaveis criadas aqui valem so dentro deste arquivo.
setlocal
title PBL-6 - simulador da inducao com propofol
rem Duplo clique: instala o Julia se faltar, prepara o projeto e abre a
rem plataforma PBL-6 no navegador. Deixe esta janela aberta durante o uso.
rem Como usar: extraia o ZIP, abra a pasta do projeto e de duplo clique neste
rem arquivo. Nao precisa digitar nada. Os dados nao saem do computador.
rem Na primeira vez o Julia baixa e compila os pacotes: 10 a 20 min.
rem cd /d vai para a pasta onde este arquivo esta, a pasta do projeto.
cd /d "%~dp0"

rem Confere se Project.toml existe aqui. Se nao existir, o arquivo foi aberto
rem de dentro do ZIP ou fora da pasta: pula para a mensagem de ajuda no fim.
if not exist "Project.toml" goto :sem_projeto

rem Procura o Julia. Primeiro no PATH, com o comando where: se achar,
rem errorlevel fica 0 e pula para o trecho :rodar. Se nao achar, tenta o
rem local de instalacao da Microsoft Store, na pasta WindowsApps do usuario.
set "JULIA=julia"
where julia >nul 2>nul
if not errorlevel 1 goto :rodar
set "JULIA=%LOCALAPPDATA%\Microsoft\WindowsApps\julia.exe"
if exist "%JULIA%" goto :rodar

rem Julia nao encontrado: instala pela Microsoft Store com o winget, o
rem instalador de programas do Windows. Precisa de internet.
echo.
echo Julia nao encontrado. Instalando pela Microsoft Store via winget...
echo Aceite as permissoes que o Windows pedir.
echo.
winget install --name Julia --id 9NJNWW8PVKMN -e -s msstore --accept-package-agreements --accept-source-agreements
if exist "%JULIA%" goto :rodar
rem Se a instalacao falhou: mostra o site oficial e encerra com erro.
echo.
echo Nao consegui instalar o Julia automaticamente.
echo Instale pelo site https://julialang.org/downloads/ e rode este arquivo de novo.
pause
exit /b 1

rem Trecho principal: mostra os avisos e inicia a plataforma.
:rodar
echo.
echo ================================================================
echo  PBL-6 - preparando a plataforma
echo  Na primeira vez o Julia baixa e compila os pacotes: 10 a 20 min.
echo  Depois disso abre em menos de 1 minuto.
echo  O navegador abre sozinho. NAO feche esta janela enquanto usa.
echo ================================================================
echo.
rem Inicia o Julia na pasta do projeto com os pacotes do Project.toml.
rem A opcao -t auto usa todos os nucleos do processador.
rem Pkg.instantiate baixa e instala os pacotes que faltam.
rem plataforma abre o servidor local e o navegador em http://localhost:8080
rem e wait mantem tudo rodando ate a janela ser fechada.
"%JULIA%" --project=. -t auto -e "using Pkg; Pkg.instantiate(); using PropofolWindkessel; wait(plataforma())"
rem Chega aqui quando a plataforma e encerrada.
echo.
echo A plataforma foi encerrada.
pause
exit /b 0

rem Mensagem de ajuda quando este arquivo nao esta dentro da pasta do projeto.
:sem_projeto
echo.
echo Este arquivo precisa estar dentro da pasta do projeto.
echo Se voce abriu o ZIP, extraia primeiro: botao direito no ZIP, Extrair tudo,
echo e rode o abrir_plataforma.bat de dentro da pasta extraida.
pause
exit /b 1
