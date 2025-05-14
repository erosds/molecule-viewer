import os
import shutil
from fastapi import FastAPI, HTTPException, UploadFile, File, Query
from fastapi.responses import FileResponse, JSONResponse
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel
from typing import List, Optional, Dict, Any
import uvicorn
import logging
from pathlib import Path
import tempfile
import csv
import io

# Importa le funzioni di utilità per le molecole
from molecule_utils import smiles_to_3d, smiles_to_svg, get_molecule_properties

# Configurazione logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Creazione cartelle necessarie
CSV_DIR = "public/csv"
MOLECULES_DIR = "public/molecules"
SVG_DIR = "public/svg"

# Assicurati che le cartelle esistano
os.makedirs(CSV_DIR, exist_ok=True)
os.makedirs(MOLECULES_DIR, exist_ok=True)
os.makedirs(SVG_DIR, exist_ok=True)

# Modelli di dati
class SMILESRequest(BaseModel):
    smiles: str

class ModelResponse(BaseModel):
    model_path: str
    svg_path: Optional[str] = None
    properties: Optional[Dict[str, Any]] = None
    success: bool
    message: Optional[str] = None

class MoleculeData(BaseModel):
    id: int
    smiles: str
    name: Optional[str] = None
    formula: Optional[str] = None

# Inizializzazione FastAPI
app = FastAPI(title="Molecular Viewer API")

# Configurazione CORS per permettere al frontend di comunicare con l'API
app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:3000"],  # In produzione, limitare alle origini specifiche
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Servire file statici
app.mount("/csv", StaticFiles(directory=CSV_DIR), name="csv")
app.mount("/molecules", StaticFiles(directory=MOLECULES_DIR), name="molecules")
app.mount("/svg", StaticFiles(directory=SVG_DIR), name="svg")

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

# Endpoint per leggere il contenuto di un file CSV specifico
@app.get("/api/csv/{filename}", response_model=List[MoleculeData])
async def read_csv_file(filename: str, limit: int = Query(100, gt=0, le=1000)):
    try:
        file_path = os.path.join(CSV_DIR, filename)
        if not os.path.exists(file_path):
            raise HTTPException(status_code=404, detail=f"File '{filename}' non trovato")
        
        molecules = []
        with open(file_path, 'r', newline='', encoding='utf-8') as csvfile:
            # Determina il dialetto e intestazioni
            sample = csvfile.read(1024)
            csvfile.seek(0)
            dialect = csv.Sniffer().sniff(sample)
            has_header = csv.Sniffer().has_header(sample)
            
            reader = csv.reader(csvfile, dialect)
            
            # Legge le intestazioni
            headers = next(reader) if has_header else []
            headers_lower = [h.lower() for h in headers]
            
            # Trova gli indici delle colonne rilevanti
            smiles_idx = next((i for i, h in enumerate(headers_lower) if 'smiles' in h), None)
            name_idx = next((i for i, h in enumerate(headers_lower) if 'name' in h or 'nome' in h), None)
            formula_idx = next((i for i, h in enumerate(headers_lower) if 'formula' in h or 'molecular_formula' in h), None)
            
            if smiles_idx is None:
                raise HTTPException(status_code=400, detail="Nessuna colonna SMILES trovata nel file CSV")
            
            # Legge le righe e crea oggetti molecola
            for i, row in enumerate(reader):
                if i >= limit:
                    break
                    
                if len(row) > smiles_idx:
                    smiles = row[smiles_idx].strip()
                    if smiles:  # Ignora le righe con SMILES vuoti
                        molecule = {
                            "id": i,
                            "smiles": smiles
                        }
                        
                        # Aggiunge il nome se disponibile
                        if name_idx is not None and len(row) > name_idx:
                            molecule["name"] = row[name_idx].strip()
                            
                        # Aggiunge la formula se disponibile
                        if formula_idx is not None and len(row) > formula_idx:
                            molecule["formula"] = row[formula_idx].strip()
                            
                        molecules.append(MoleculeData(**molecule))
        
        return molecules
    except HTTPException as he:
        raise he
    except Exception as e:
        logger.error(f"Errore nella lettura del file CSV: {str(e)}")
        raise HTTPException(status_code=500, detail=f"Errore nella lettura del file: {str(e)}")

# Endpoint per generare modello 3D da SMILES
@app.post("/api/generate-3d", response_model=ModelResponse)
async def generate_3d_model(request: SMILESRequest):
    try:
        # Validazione SMILES
        if not request.smiles or len(request.smiles) < 2:
            raise HTTPException(status_code=400, detail="SMILES non valido")
        
        # Crea un nome file sicuro basato sullo SMILES
        import hashlib
        smiles_hash = hashlib.md5(request.smiles.encode()).hexdigest()
        
        model_filename = f"{smiles_hash}.pdb"
        svg_filename = f"{smiles_hash}.svg"
        
        model_path = os.path.join(MOLECULES_DIR, model_filename)
        svg_path = os.path.join(SVG_DIR, svg_filename)
        
        # Controlla se il file esiste già
        if not os.path.exists(model_path):
            # Genera il modello 3D
            success = smiles_to_3d(request.smiles, model_path)
            if not success:
                raise HTTPException(status_code=500, detail="Impossibile generare il modello 3D")
        
        # Genera SVG se non esiste
        if not os.path.exists(svg_path):
            svg_result = smiles_to_svg(request.smiles, svg_path)
            if not svg_result:
                logger.warning(f"Impossibile generare SVG per SMILES: {request.smiles}")
        
        # Calcola proprietà molecolari
        properties = get_molecule_properties(request.smiles)
        
        # Restituisci i percorsi relativi per l'accesso dal frontend
        relative_model_path = f"/molecules/{model_filename}"
        relative_svg_path = f"/svg/{svg_filename}" if os.path.exists(svg_path) else None
        
        return ModelResponse(
            model_path=relative_model_path,
            svg_path=relative_svg_path,
            properties=properties,
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
    return {"status": "healthy"}

# Avvio dell'applicazione
if __name__ == "__main__":
    uvicorn.run("main:app", host="0.0.0.0", port=8000, reload=True)