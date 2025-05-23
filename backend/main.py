import os
import shutil
import csv
import io
from fastapi import FastAPI, HTTPException, UploadFile, File, Response
from fastapi.responses import FileResponse
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel
from typing import List, Optional
import uvicorn
import logging
from pathlib import Path
import tempfile

# Importa le funzioni di utilità per le molecole
from molecule_utils import smiles_to_3d, smiles_to_2d_image

# Configurazione logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Creazione cartelle necessarie
CSV_DIR = "public/csv"
MOLECULES_DIR = "public/molecules"
IMAGES_DIR = "public/images"

# Assicurati che le cartelle esistano
os.makedirs(CSV_DIR, exist_ok=True)
os.makedirs(MOLECULES_DIR, exist_ok=True)
os.makedirs(IMAGES_DIR, exist_ok=True)

# Modelli di dati
class SMILESRequest(BaseModel):
    smiles: str

class ModelResponse(BaseModel):
    model_path: str
    success: bool
    message: Optional[str] = None

# Inizializzazione FastAPI
app = FastAPI(title="Molecular Viewer API")

# Configurazione CORS per permettere al frontend di comunicare con l'API
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # In produzione, limitare alle origini specifiche
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Servire file statici
app.mount("/csv", StaticFiles(directory=CSV_DIR), name="csv")
app.mount("/molecules", StaticFiles(directory=MOLECULES_DIR), name="molecules")
app.mount("/images", StaticFiles(directory=IMAGES_DIR), name="images")

# Endpoint per ottenere la lista dei file CSV disponibili
@app.get("/api/csv-files", response_model=List[str])
async def get_csv_files():
    try:
        # Filtra solo i file .csv nella directory
        csv_files = [f for f in os.listdir(CSV_DIR) if f.endswith('.csv')]
        return csv_files
    except Exception as e:
        logger.error(f"Errore nel recupero dei file CSV: {str(e)}")
        raise HTTPException(status_code=500, detail=f"Errore nel recupero dei file: {str(e)}")

# Endpoint per caricare un nuovo file CSV
@app.post("/api/upload-csv")
async def upload_csv_file(file: UploadFile = File(...)):
    try:
        if not file.filename.endswith('.csv'):
            raise HTTPException(status_code=400, detail="Il file deve essere in formato CSV")
        
        # Salva il file nella directory CSV
        file_path = os.path.join(CSV_DIR, file.filename)
        with open(file_path, "wb") as buffer:
            shutil.copyfileobj(file.file, buffer)
        
        return {"filename": file.filename, "success": True}
    except Exception as e:
        logger.error(f"Errore nel caricamento del file: {str(e)}")
        raise HTTPException(status_code=500, detail=f"Errore nel caricamento del file: {str(e)}")

# Endpoint per ottenere le molecole da un file CSV
@app.get("/api/molecules/{csv_file}")
async def get_molecules_from_csv(csv_file: str):
    try:
        file_path = os.path.join(CSV_DIR, csv_file)
        if not os.path.exists(file_path):
            raise HTTPException(status_code=404, detail=f"File CSV non trovato: {csv_file}")
        
        molecules = []
        with open(file_path, 'r') as f:
            csv_reader = csv.reader(f)
            headers = [h.lower() for h in next(csv_reader)]
            
            # Trova l'indice della colonna SMILES
            smiles_index = None
            for i, header in enumerate(headers):
                if 'smiles' in header:
                    smiles_index = i
                    break
            
            if smiles_index is None:
                raise HTTPException(status_code=400, detail=f"Nessuna colonna SMILES trovata nel file {csv_file}")
            
            # Leggi tutte le righe e estrai i dati SMILES
            for i, row in enumerate(csv_reader):
                if len(row) > smiles_index:
                    smiles = row[smiles_index].strip()
                    if smiles:  # Verifica che il valore SMILES non sia vuoto
                        molecules.append({
                            "id": i,
                            "smiles": smiles
                        })
        
        return molecules
    except HTTPException as he:
        raise he
    except Exception as e:
        logger.error(f"Errore nella lettura del file CSV: {str(e)}")
        raise HTTPException(status_code=500, detail=f"Errore nella lettura del file: {str(e)}")

# Endpoint per generare immagine 2D da SMILES
@app.get("/api/molecule-2d/{smiles}")
async def get_molecule_2d_image(smiles: str):
    try:
        # Crea un nome file sicuro basato sullo SMILES
        safe_filename = "".join(c if c.isalnum() else "_" for c in smiles)
        if len(safe_filename) > 50:
            safe_filename = safe_filename[:50]
        
        image_filename = f"{safe_filename}.png"
        image_path = os.path.join(IMAGES_DIR, image_filename)
        
        # Verifica se l'immagine esiste già
        if os.path.exists(image_path):
            return FileResponse(image_path, media_type="image/png")
        
        # Genera l'immagine 2D direttamente come bytes
        image_bytes = smiles_to_2d_image(smiles)
        if not image_bytes:
            raise HTTPException(status_code=500, detail="Impossibile generare l'immagine della molecola")
        
        # Salva l'immagine per richieste future
        with open(image_path, "wb") as f:
            f.write(image_bytes)
        
        # Restituisci l'immagine come risposta
        return Response(content=image_bytes, media_type="image/png")
    except HTTPException as he:
        raise he
    except Exception as e:
        logger.error(f"Errore nella generazione dell'immagine 2D: {str(e)}")
        raise HTTPException(status_code=500, detail=f"Errore interno: {str(e)}")

# Endpoint per generare modello 3D da SMILES
@app.post("/api/generate-3d", response_model=ModelResponse)
async def generate_3d_model(request: SMILESRequest):
    try:
        # Validazione SMILES
        if not request.smiles or len(request.smiles) < 2:
            raise HTTPException(status_code=400, detail="SMILES non valido")
        
        # Crea un nome file sicuro basato sullo SMILES
        safe_filename = "".join(c if c.isalnum() else "_" for c in request.smiles)
        if len(safe_filename) > 50:
            safe_filename = safe_filename[:50]
        
        model_filename = f"{safe_filename}.xyz"  # Cambiato da .pdb a .xyz
        model_path = os.path.join(MOLECULES_DIR, model_filename)
        
        # Controlla se il file esiste già
        if not os.path.exists(model_path):
            # Genera il modello 3D
            success = smiles_to_3d(request.smiles, model_path)
            if not success:
                raise HTTPException(status_code=500, detail="Impossibile generare il modello 3D")
        
        # Restituisci il percorso relativo per l'accesso dal frontend
        relative_path = f"/molecules/{model_filename}"
        
        return ModelResponse(
            model_path=relative_path,
            success=True,
            message="Modello 3D generato con successo"
        )
    except HTTPException as he:
        # Rilancia le eccezioni HTTP
        raise he
    except Exception as e:
        logger.error(f"Errore nella generazione del modello 3D: {str(e)}")
        raise HTTPException(status_code=500, detail=f"Errore interno: {str(e)}")

# Endpoint per salute dell'API
@app.get("/api/health")
async def health_check():
    return {"status": "healthyss"}

# Avvio dell'applicazione
if __name__ == "__main__":
    uvicorn.run("main:app", host="0.0.0.0", port=8001, reload=True)