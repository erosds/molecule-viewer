@echo off
echo Avvio del Visualizzatore Molecolare...

REM Avvio del backend
echo Avvio del backend...
start cmd /k "cd backend && if exist venv\Scripts\activate.bat (call venv\Scripts\activate.bat) && python main.py"

REM Attendi 3 secondi per dare tempo al backend di avviarsi
timeout /t 3 /nobreak > nul

REM Avvio del frontend
echo Avvio del frontend...
start cmd /k "cd frontend && npm start"

echo.
echo Applicazione avviata!
echo Frontend: http://localhost:3000
echo Backend: http://localhost:8001
echo.
echo Entrambi i processi sono stati avviati in finestre separate.
echo Chiudi le finestre di comando per terminare l'applicazione.
echo.