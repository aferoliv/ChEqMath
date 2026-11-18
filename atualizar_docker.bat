@echo off
setlocal
cd /d "%~dp0"

docker info >nul 2>&1
if errorlevel 1 (
  echo Docker nao esta em execucao. Abra o Docker Desktop e tente novamente.
  pause
  exit /b 1
)

:menu
cls
echo ===============================================
echo   Especiacao - Atualizacao Docker e Pagina
echo ===============================================
echo [1] Atualizar tudo (frontend + backend, sem cache)
echo [2] Atualizar so frontend
echo [3] Atualizar so backend
echo [4] Reimportar planilha (reinicia o backend)
echo [5] Ver status dos containers
echo [6] Ver logs do frontend (web)
echo [7] Ver logs do backend (api)
echo [8] Abrir a pagina no navegador
echo [9] Sair
echo.
set /p choice=Escolha uma opcao: 

if "%choice%"=="1" goto all
if "%choice%"=="2" goto web
if "%choice%"=="3" goto api
if "%choice%"=="4" goto reimport
if "%choice%"=="5" goto ps
if "%choice%"=="6" goto logs_web
if "%choice%"=="7" goto logs_api
if "%choice%"=="8" goto open
if "%choice%"=="9" goto end

echo Opcao invalida.
pause
goto menu

:all
echo.
echo Atualizando tudo...
docker compose build --no-cache web api
if errorlevel 1 goto fail
docker compose up -d --force-recreate web api
goto done

:web
echo.
echo Atualizando frontend...
docker compose build web
if errorlevel 1 goto fail
docker compose up -d web
goto done

:api
echo.
echo Atualizando backend...
docker compose build api
if errorlevel 1 goto fail
docker compose up -d api
goto done

:reimport
echo.
echo Reimportando a planilha de data\ (a planilha vai dentro da imagem: rebuild do backend)...
docker compose build api
if errorlevel 1 goto fail
docker compose up -d api
goto done

:ps
echo.
docker compose ps
goto done

:logs_web
echo.
docker compose logs -f web
goto done

:logs_api
echo.
docker compose logs -f api
goto done

:open
start "" http://localhost:8090
goto menu

:fail
echo.
echo Falha no build. Veja a mensagem acima.
pause
goto menu

:done
echo.
echo Se nao refletir no navegador, faca Ctrl+F5.
pause
goto menu

:end
endlocal
exit /b 0
