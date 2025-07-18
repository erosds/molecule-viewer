import os
import shutil
import csv
import io
from fastapi import FastAPI, HTTPException, UploadFile, File, Response
from fastapi.responses import FileResponse
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel
from typing import List, Optional, Dict
import uvicorn
import logging
from pathlib import Path
import tempfile
from datetime import datetime
import zipfile
from fastapi.responses import StreamingResponse
from rdkit import Chem
from typing import Any  # Aggiungi questo import se manca

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

from typing import Dict
from concurrent.futures import ThreadPoolExecutor
import threading

# Aggiungi questo modello dopo gli altri
class ValidationRequest(BaseModel):
    smiles_list: List[str]

class ValidationResponse(BaseModel):
    total_molecules: int
    valid_molecules: int
    invalid_molecules: int
    valid_smiles: List[str]
    invalid_smiles: List[str]
    processing_time: float

class NoveltyAnalysisRequest(BaseModel):
    main_file: str
    reference_file: Optional[str] = None
    
class NoveltyAnalysisResponse(BaseModel):
    total_molecules: int
    valid_molecules: int
    invalid_molecules: int
    unique_molecules: int
    novel_molecules: int
    uniqueness_rate: float
    novelty_rate: float
    processing_time: float
    molecules: List[Dict]  # Lista con info di novelty per ogni molecola
    
class CoordinationFilterRequest(BaseModel):
    csv_file: str
    min_coordination: int = 0
    max_coordination: int = 12
    selected_metals: Optional[List[str]] = None
    include_non_metal_molecules: bool = False  # <-- AGGIUNGI QUESTO


class CoordinationAnalysisResponse(BaseModel):
    total_molecules: int
    molecules_with_metals: int
    molecules_without_metals: int  # <-- AGGIUNGI QUESTO
    filtered_molecules: int
    coordination_distribution: Dict[int, int]
    molecules: List[Dict]
    selected_metals: Optional[List[str]] = None  # Aggiungi questo campo
    
class NonMetalFilterRequest(BaseModel):
    csv_file: str

class NonMetalFilterResponse(BaseModel):
    total_molecules: int
    molecules_without_metals: int
    filtered_molecules: int
    molecules: List[Dict]
    
# Inizializzazione FastAPI
app = FastAPI(title="Molecular Viewer API")
    
@app.post("/api/coordination-stats")
async def get_coordination_statistics_post(request: dict):
    """
    Ottiene statistiche sulla distribuzione del numero di coordinazione in un file CSV (versione POST).
    """
    csv_file = request.get('csv_file')
    if not csv_file:
        raise HTTPException(status_code=400, detail="csv_file è richiesto")
    
    return await get_coordination_statistics(csv_file)

@app.post("/api/filter-non-metal-molecules", response_model=NonMetalFilterResponse)
async def filter_non_metal_molecules(request: NonMetalFilterRequest):
    """
    Filtra le molecole per mostrare solo quelle senza metalli.
    """
    try:
        # Carica le molecole dal file CSV
        file_path = os.path.join(CSV_DIR, request.csv_file)
        if not os.path.exists(file_path):
            raise HTTPException(status_code=404, detail=f"File CSV non trovato: {request.csv_file}")
        
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
                raise HTTPException(status_code=400, detail=f"Nessuna colonna SMILES trovata nel file {request.csv_file}")
            
            # Leggi tutte le righe e estrai i dati SMILES
            for i, row in enumerate(csv_reader):
                if len(row) > smiles_index:
                    smiles = row[smiles_index].strip()
                    if smiles:
                        molecules.append({
                            "id": i,
                            "smiles": smiles
                        })
        
        # Importa la funzione di analisi
        from molecule_utils import analyze_metal_coordination
        
        # Filtra le molecole per trovare solo quelle senza metalli
        non_metal_molecules = []
        molecules_with_metals = 0
        
        for molecule in molecules:
            coord_info = analyze_metal_coordination(molecule['smiles'])
            
            if coord_info['has_metal']:
                molecules_with_metals += 1
            else:
                # Questa molecola non ha metalli, la includiamo nel risultato
                molecule_copy = molecule.copy()
                molecule_copy['coordination_info'] = coord_info
                non_metal_molecules.append(molecule_copy)
        
        # Prepara la risposta
        response_data = {
            "total_molecules": len(molecules),
            "molecules_without_metals": len(non_metal_molecules),
            "filtered_molecules": len(non_metal_molecules),
            "molecules": non_metal_molecules
        }
        
        logger.info(f"Filtro molecole senza metalli completato: {len(non_metal_molecules)}/{len(molecules)} molecole senza metalli")
        
        return NonMetalFilterResponse(**response_data)
        
    except Exception as e:
        logger.error(f"Errore nel filtro molecole senza metalli: {e}")
        raise HTTPException(status_code=500, detail=f"Errore interno del server: {str(e)}")


# Aggiungi anche questa funzione di utilità al file molecule_utils.py

def filter_molecules_without_metals(molecules_data: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """
    Filtra una lista di molecole per restituire solo quelle senza metalli.
    
    Args:
        molecules_data: Lista di dizionari con dati delle molecole
        
    Returns:
        list: Lista filtrata di molecole senza metalli
    """
    if not isinstance(molecules_data, list):
        logger.error("molecules_data deve essere una lista")
        return []
    
    filtered_molecules = []
    
    for molecule in molecules_data:
        if not isinstance(molecule, dict):
            continue
            
        smiles = molecule.get('smiles', '')
        if not smiles:
            continue
            
        coordination_info = analyze_metal_coordination(smiles)
        
        # Includi solo le molecole che NON hanno metalli
        if not coordination_info['has_metal']:
            molecule_copy = molecule.copy()
            molecule_copy['coordination_info'] = coordination_info
            filtered_molecules.append(molecule_copy)
    
    return filtered_molecules


def get_molecules_statistics(molecules_data: List[Dict[str, Any]]) -> Dict[str, Any]:
    """
    Calcola statistiche sui metalli per una lista di molecole.
    
    Args:
        molecules_data: Lista di dizionari con dati delle molecole
        
    Returns:
        dict: Statistiche con conteggi di molecole con/senza metalli
    """
    if not isinstance(molecules_data, list):
        return {
            'total_molecules': 0,
            'molecules_with_metals': 0,
            'molecules_without_metals': 0
        }
    
    total_molecules = len(molecules_data)
    molecules_with_metals = 0
    
    for molecule in molecules_data:
        if not isinstance(molecule, dict):
            continue
            
        smiles = molecule.get('smiles', '')
        if not smiles:
            continue
            
        coordination_info = analyze_metal_coordination(smiles)
        if coordination_info['has_metal']:
            molecules_with_metals += 1
    
    molecules_without_metals = total_molecules - molecules_with_metals
    
    return {
        'total_molecules': total_molecules,
        'molecules_with_metals': molecules_with_metals,
        'molecules_without_metals': molecules_without_metals
    }

# Funzione helper per validazione veloce di una singola molecola
def validate_single_smiles(smiles: str) -> bool:
    """
    Validazione veloce di un SMILES senza generare file 3D.
    Verifica solo se RDKit può creare una molecola valida.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        
        # Conversione da SMILES a molecola RDKit
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False
            
        # Aggiunta degli idrogeni
        mol = Chem.AddHs(mol)
        
        # Prova a generare una conformazione 3D (senza salvarla)
        result = AllChem.EmbedMolecule(mol, AllChem.ETKDG())
        
        # Se result è -1, la conformazione non è stata generata
        return result != -1
        
    except Exception as e:
        logger.error(f"Errore nella validazione SMILES {smiles}: {str(e)}")
        return False

def analyze_novelty_and_uniqueness(main_file: str, reference_file: Optional[str] = None) -> NoveltyAnalysisResponse:
    """
    Analizza unicità e novelty delle molecole secondo la logica dello script Python fornito.
    """
    import time
    start_time = time.time()
    
    try:
        # Leggi il file principale
        main_file_path = os.path.join(CSV_DIR, main_file)
        if not os.path.exists(main_file_path):
            raise HTTPException(status_code=404, detail=f"File principale non trovato: {main_file}")
        
        input_smiles = []
        with open(main_file_path, 'r') as f:
            csv_reader = csv.reader(f)
            headers = [h.lower() for h in next(csv_reader)]
            
            # Trova l'indice della colonna SMILES
            smiles_index = None
            for i, header in enumerate(headers):
                if 'smiles' in header:
                    smiles_index = i
                    break
            
            if smiles_index is None:
                raise ValueError(f"Nessuna colonna SMILES trovata nel file {main_file}")
            
            # Leggi tutti gli SMILES
            for row in csv_reader:
                if len(row) > smiles_index:
                    smiles = row[smiles_index].strip()
                    if smiles:
                        input_smiles.append(smiles)
        
        # Leggi il file di riferimento se specificato
        reference_smiles = []
        if reference_file:
            ref_file_path = os.path.join(CSV_DIR, reference_file)
            if os.path.exists(ref_file_path):
                with open(ref_file_path, 'r') as f:
                    csv_reader = csv.reader(f)
                    headers = [h.lower() for h in next(csv_reader)]
                    
                    # Trova l'indice della colonna SMILES
                    smiles_index = None
                    for i, header in enumerate(headers):
                        if 'smiles' in header:
                            smiles_index = i
                            break
                    
                    if smiles_index is not None:
                        for row in csv_reader:
                            if len(row) > smiles_index:
                                smiles = row[smiles_index].strip()
                                if smiles:
                                    reference_smiles.append(smiles)
        
        logger.info(f"Analisi di {len(input_smiles)} SMILES del file principale")
        if reference_smiles:
            logger.info(f"Confronto con {len(reference_smiles)} SMILES di riferimento")
        
        # Analisi unicità (canonicalizza gli SMILES)
        unique_smiles = set()
        valid_molecules = []
        invalid_count = 0
        
        for i, smi in enumerate(input_smiles):
            mol = Chem.MolFromSmiles(smi)
            if mol is not None:
                canonical_smi = Chem.MolToSmiles(mol)
                unique_smiles.add(canonical_smi)
                valid_molecules.append({
                    'id': i,
                    'original_smiles': smi,
                    'canonical_smiles': canonical_smi,
                    'is_valid': True
                })
            else:
                invalid_count += 1
                valid_molecules.append({
                    'id': i,
                    'original_smiles': smi,
                    'canonical_smiles': None,
                    'is_valid': False
                })
        
        # Calcola uniqueness
        uniqueness_rate = len(unique_smiles) / len(input_smiles) if len(input_smiles) > 0 else 0
        
        # Analisi novelty se abbiamo un riferimento
        reference_unique = set()
        novel_count = 0
        novelty_rate = 0.0
        
        if reference_smiles:
            # Canonicalizza gli SMILES di riferimento
            for smi in reference_smiles:
                mol = Chem.MolFromSmiles(smi)
                if mol is not None:
                    canonical_smi = Chem.MolToSmiles(mol)
                    reference_unique.add(canonical_smi)
            
            # Conta le molecole novel
            for smi in unique_smiles:
                if smi not in reference_unique:
                    novel_count += 1
            
            novelty_rate = novel_count / len(unique_smiles) if len(unique_smiles) > 0 else 0
        
        # Aggiungi informazioni di novelty alle molecole
        molecules_with_novelty = []
        unique_canonical_to_novel = {}
        
        # Prima determina quali molecole canoniche sono novel
        if reference_smiles:
            for canonical_smi in unique_smiles:
                unique_canonical_to_novel[canonical_smi] = canonical_smi not in reference_unique
        
        # Poi aggiungi le informazioni a ogni molecola
        seen_canonical = set()
        for mol_info in valid_molecules:
            if mol_info['is_valid']:
                canonical = mol_info['canonical_smiles']
                is_unique = canonical not in seen_canonical
                is_novel = unique_canonical_to_novel.get(canonical, False) if reference_smiles else True
                
                molecules_with_novelty.append({
                    'id': mol_info['id'],
                    'smiles': mol_info['original_smiles'],
                    'canonical_smiles': canonical,
                    'is_valid': True,
                    'is_unique': is_unique,
                    'is_novel': is_novel
                })
                
                seen_canonical.add(canonical)
            else:
                molecules_with_novelty.append({
                    'id': mol_info['id'],
                    'smiles': mol_info['original_smiles'],
                    'canonical_smiles': None,
                    'is_valid': False,
                    'is_unique': False,
                    'is_novel': False
                })
        
        processing_time = time.time() - start_time
        
        logger.info(f"Analisi completata in {processing_time:.2f}s")
        logger.info(f"Risultati: {len(valid_molecules) - invalid_count} valide, {len(unique_smiles)} uniche, {novel_count} novel")
        
        return NoveltyAnalysisResponse(
            total_molecules=len(input_smiles),
            valid_molecules=len(valid_molecules) - invalid_count,
            invalid_molecules=invalid_count,
            unique_molecules=len(unique_smiles),
            novel_molecules=novel_count,
            uniqueness_rate=uniqueness_rate,
            novelty_rate=novelty_rate,
            processing_time=round(processing_time, 2),
            molecules=molecules_with_novelty
        )
        
    except Exception as e:
        logger.error(f"Errore nell'analisi di novelty: {str(e)}")
        raise HTTPException(status_code=500, detail=f"Errore interno: {str(e)}")

# Modelli di dati
class SMILESRequest(BaseModel):
    smiles: str

class ModelResponse(BaseModel):
    model_path: str
    success: bool
    message: Optional[str] = None

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

# Nuovo endpoint per l'analisi di novelty e uniqueness
@app.post("/api/analyze-novelty", response_model=NoveltyAnalysisResponse)
async def analyze_novelty_endpoint(request: NoveltyAnalysisRequest):
    """
    Analizza unicità e novelty delle molecole confrontandole con un file di riferimento.
    """
    try:
        return analyze_novelty_and_uniqueness(request.main_file, request.reference_file)
    except HTTPException as he:
        raise he
    except Exception as e:
        logger.error(f"Errore nell'analisi di novelty: {str(e)}")
        raise HTTPException(status_code=500, detail=f"Errore interno: {str(e)}")

# Endpoint per ottenere le molecole da un file CSV (deprecato, sostituito da analyze-novelty)
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
    
@app.post("/api/validate-molecules", response_model=ValidationResponse)
async def validate_molecules_bulk(request: ValidationRequest):
    """
    Valida rapidamente una lista di SMILES per verificare
    quali possono generare strutture 3D valide.
    """
    import time
    start_time = time.time()
    
    try:
        if not request.smiles_list:
            raise HTTPException(status_code=400, detail="Lista SMILES vuota")
        
        valid_smiles = []
        invalid_smiles = []
        
        # Utilizziamo ThreadPoolExecutor per parallelizzare le validazioni
        # Limitiamo il numero di thread per evitare sovraccarico
        max_workers = min(4, len(request.smiles_list))
        
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            # Mappa ogni SMILES al suo future
            future_to_smiles = {
                executor.submit(validate_single_smiles, smiles): smiles 
                for smiles in request.smiles_list
            }
            
            # Raccoglie i risultati
            for future in future_to_smiles:
                smiles = future_to_smiles[future]
                try:
                    is_valid = future.result(timeout=10)  # Timeout di 10 secondi per molecola
                    if is_valid:
                        valid_smiles.append(smiles)
                    else:
                        invalid_smiles.append(smiles)
                except Exception as e:
                    logger.error(f"Errore nella validazione di {smiles}: {str(e)}")
                    invalid_smiles.append(smiles)
        
        processing_time = time.time() - start_time
        
        return ValidationResponse(
            total_molecules=len(request.smiles_list),
            valid_molecules=len(valid_smiles),
            invalid_molecules=len(invalid_smiles),
            valid_smiles=valid_smiles,
            invalid_smiles=invalid_smiles,
            processing_time=round(processing_time, 2)
        )
        
    except HTTPException as he:
        raise he
    except Exception as e:
        logger.error(f"Errore nella validazione bulk: {str(e)}")
        raise HTTPException(status_code=500, detail=f"Errore interno: {str(e)}")

@app.post("/api/filter-by-coordination", response_model=CoordinationAnalysisResponse)
async def filter_molecules_by_coordination_endpoint(request: CoordinationFilterRequest):
    """
    Filtra le molecole di un file CSV in base al numero di coordinazione metallica e tipo di metallo.
    """
    try:
        # Carica le molecole dal file CSV
        file_path = os.path.join(CSV_DIR, request.csv_file)
        if not os.path.exists(file_path):
            raise HTTPException(status_code=404, detail=f"File CSV non trovato: {request.csv_file}")
        
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
                raise HTTPException(status_code=400, detail=f"Nessuna colonna SMILES trovata nel file {request.csv_file}")
            
            # Leggi tutte le righe e estrai i dati SMILES
            for i, row in enumerate(csv_reader):
                if len(row) > smiles_index:
                    smiles = row[smiles_index].strip()
                    if smiles:
                        molecules.append({
                            "id": i,
                            "smiles": smiles
                        })
        
        # Importa le funzioni di filtro
        from molecule_utils import filter_molecules_by_coordination, analyze_metal_coordination, filter_molecules_by_metal_and_coordination
        
        # Analizza tutte le molecole per ottenere statistiche
        molecules_with_metals = 0
        coordination_distribution = {}
        
        for molecule in molecules:
            coord_info = analyze_metal_coordination(molecule['smiles'])
            if coord_info['has_metal']:
                molecules_with_metals += 1
                max_coord = coord_info['max_coordination']
                coordination_distribution[max_coord] = coordination_distribution.get(max_coord, 0) + 1
        
        # Calcola anche le molecole senza metalli per le statistiche
        molecules_without_metals = len(molecules) - molecules_with_metals
        
        # Applica il filtro in base ai parametri
        if request.include_non_metal_molecules:
            # Se include_non_metal_molecules è True, filtra solo le molecole SENZA metalli
            from molecule_utils import filter_molecules_without_metals
            filtered_molecules = filter_molecules_without_metals(molecules)
        else:
            # Filtro normale per molecole CON metalli
            if request.selected_metals:
                filtered_molecules = filter_molecules_by_metal_and_coordination(
                    molecules, 
                    request.min_coordination, 
                    request.max_coordination,
                    request.selected_metals
                )
            else:
                filtered_molecules = filter_molecules_by_coordination(
                    molecules, 
                    request.min_coordination, 
                    request.max_coordination
                )
        
        return CoordinationAnalysisResponse(
            total_molecules=len(molecules),
            molecules_with_metals=molecules_with_metals,
            filtered_molecules=len(filtered_molecules),
            coordination_distribution=coordination_distribution,
            molecules_without_metals=molecules_without_metals,  # <-- AGGIUNGI QUESTO
            molecules=filtered_molecules,
            selected_metals=request.selected_metals if request.selected_metals else []
        )
        
    except HTTPException as he:
        raise he
    except Exception as e:
        logger.error(f"Errore nel filtro di coordinazione: {str(e)}")
        raise HTTPException(status_code=500, detail=f"Errore interno: {str(e)}")
    
@app.get("/api/coordination-stats/{csv_file}")
async def get_coordination_statistics(csv_file: str):
    """
    Ottiene statistiche sulla distribuzione del numero di coordinazione in un file CSV.
    """
    try:
        file_path = os.path.join(CSV_DIR, csv_file)
        if not os.path.exists(file_path):
            raise HTTPException(status_code=404, detail=f"File CSV non trovato: {csv_file}")
        
        molecules = []
        with open(file_path, 'r') as f:
            csv_reader = csv.reader(f)
            headers = [h.lower() for h in next(csv_reader)]
            
            smiles_index = None
            for i, header in enumerate(headers):
                if 'smiles' in header:
                    smiles_index = i
                    break
            
            if smiles_index is None:
                raise HTTPException(status_code=400, detail=f"Nessuna colonna SMILES trovata nel file {csv_file}")
            
            for i, row in enumerate(csv_reader):
                if len(row) > smiles_index:
                    smiles = row[smiles_index].strip()
                    if smiles:
                        molecules.append(smiles)
        
        from molecule_utils import analyze_metal_coordination
        
        # Inizializza contatori e strutture dati
        coordination_stats = {}
        metal_elements = set()
        molecules_with_metals = 0
        molecules_without_metals = 0
        molecule_records = []
        
        # Analizza ogni molecola
        for i, smiles in enumerate(molecules):
            coord_info = analyze_metal_coordination(smiles)
            
            if coord_info['has_metal']:
                molecules_with_metals += 1
                max_coord = coord_info['max_coordination']
                coordination_stats[max_coord] = coordination_stats.get(max_coord, 0) + 1
                
                # Raccogli tutti i metalli trovati
                for metal in coord_info['metal_atoms']:
                    metal_elements.add(metal['symbol'])
                
                # Salva info per il frontend
                molecule_records.append({
                    "id": i,
                    "smiles": smiles,
                    "max_coordination": max_coord,
                    "metals": [m['symbol'] for m in coord_info['metal_atoms']]
                })
            else:
                molecules_without_metals += 1

        # Calcola range di coordinazione
        coordination_range = {
            "min": min(coordination_stats.keys()) if coordination_stats else 0,
            "max": max(coordination_stats.keys()) if coordination_stats else 0
        }

        response_data = {
            "total_molecules": len(molecules),
            "molecules_with_metals": molecules_with_metals,
            "molecules_without_metals": molecules_without_metals,
            "coordination_distribution": coordination_stats,
            "metal_elements_found": sorted(list(metal_elements)),
            "coordination_range": coordination_range,
            "molecules": molecule_records
        }
        
        logger.info(f"Statistiche calcolate per {csv_file}: {len(molecules)} totali, {molecules_with_metals} con metalli, {molecules_without_metals} senza metalli")
        
        return response_data
        
    except HTTPException as he:
        raise he
    except Exception as e:
        logger.error(f"Errore nelle statistiche di coordinazione: {str(e)}")
        raise HTTPException(status_code=500, detail=f"Errore interno: {str(e)}")
    
class BatchGenerationRequest(BaseModel):
    smiles_list: List[str]
    
class BatchGenerationResponse(BaseModel):
    total_requested: int
    successfully_generated: int 
    failed_generation: int
    generated_files: List[str]
    failed_smiles: List[str]
    processing_time: float
    batch_folder: str  # Nuovo campo


def generate_xyz_batch(smiles_list: List[str]) -> BatchGenerationResponse:
    """
    Genera file XYZ per una lista di SMILES in parallelo.
    """
    import time
    from concurrent.futures import ThreadPoolExecutor, as_completed
    import hashlib
    
    start_time = time.time()
    
     # Crea sottocartella con timestamp
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    batch_folder = f"molecules_{timestamp}"
    batch_path = os.path.join(MOLECULES_DIR, batch_folder)
    os.makedirs(batch_path, exist_ok=True)
    
    successfully_generated = 0
    failed_generation = 0
    generated_files = []
    failed_smiles = []
    
    def generate_single_xyz(smiles_with_index: tuple) -> tuple:
        """Genera un singolo file XYZ e restituisce (successo, nome_file, smiles)"""
        smiles, index = smiles_with_index
        try:
            # Nome file progressivo
            safe_filename = f"molecule_{index + 1}.xyz"
            output_path = os.path.join(batch_path, safe_filename)
            
            # Se il file esiste già, considera come successo
            if os.path.exists(output_path):
                return True, safe_filename, smiles
            
            # Genera il file XYZ
            success = smiles_to_3d(smiles, output_path)
            if success:
                return True, safe_filename, smiles
            else:
                return False, None, smiles
                
        except Exception as e:
            logger.error(f"Errore nella generazione XYZ per {smiles}: {str(e)}")
            return False, None, smiles
    
    # Usa ThreadPoolExecutor per parallelizzare la generazione
    max_workers = min(4, len(smiles_list))  # Limita il numero di thread
    
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        # Invia tutti i job con indice
        future_to_smiles = {
            executor.submit(generate_single_xyz, (smiles, i)): smiles 
            for i, smiles in enumerate(smiles_list)
        }
        
        # Raccoglie i risultati man mano che completano
        for future in as_completed(future_to_smiles):
            try:
                success, filename, smiles = future.result(timeout=30)  # Timeout di 30 secondi per molecola
                if success:
                    successfully_generated += 1
                    generated_files.append(filename)
                else:
                    failed_generation += 1
                    failed_smiles.append(smiles)
            except Exception as e:
                failed_generation += 1
                failed_smiles.append(future_to_smiles[future])
                logger.error(f"Errore nel future per {future_to_smiles[future]}: {str(e)}")
    
    processing_time = time.time() - start_time
    
    return BatchGenerationResponse(
        total_requested=len(smiles_list),
        successfully_generated=successfully_generated,
        failed_generation=failed_generation,
        generated_files=generated_files,
        failed_smiles=failed_smiles,
        processing_time=round(processing_time, 2),
        batch_folder=batch_folder  # Nuovo campo
    )

# Aggiungi questo endpoint dopo gli altri endpoint

@app.post("/api/generate-xyz-batch", response_model=BatchGenerationResponse)
async def generate_xyz_batch_endpoint(request: BatchGenerationRequest):
    """
    Genera file XYZ per una lista di SMILES validati.
    """
    try:
        if not request.smiles_list:
            raise HTTPException(status_code=400, detail="Lista SMILES vuota")
        
        # Limita il numero massimo di molecole per evitare sovraccarico del server
        if len(request.smiles_list) > 1000:
            raise HTTPException(
                status_code=400, 
                detail="Troppi SMILES richiesti. Massimo 1000 molecole per volta."
            )
        
        logger.info(f"Avvio generazione batch di {len(request.smiles_list)} file XYZ")
        
        result = generate_xyz_batch(request.smiles_list)
        
        logger.info(f"Generazione batch completata: {result.successfully_generated}/{result.total_requested} successi")
        
        return result
        
    except HTTPException as he:
        raise he
    except Exception as e:
        logger.error(f"Errore nella generazione batch XYZ: {str(e)}")
        raise HTTPException(status_code=500, detail=f"Errore interno: {str(e)}")

@app.get("/api/download-xyz/{filename}")
async def download_single_xyz(filename: str):
    """
    Endpoint per scaricare un singolo file XYZ.
    """
    try:
        # Cerca il file nella cartella molecules
        file_path = os.path.join(MOLECULES_DIR, filename)
        
        # Se non esiste nella cartella principale, cerca nelle sottocartelle
        if not os.path.exists(file_path):
            # Cerca in tutte le sottocartelle molecules_*
            for subdir in os.listdir(MOLECULES_DIR):
                subdir_path = os.path.join(MOLECULES_DIR, subdir)
                if os.path.isdir(subdir_path) and subdir.startswith('molecules_'):
                    potential_file = os.path.join(subdir_path, filename)
                    if os.path.exists(potential_file):
                        file_path = potential_file
                        break
        
        if not os.path.exists(file_path):
            raise HTTPException(status_code=404, detail=f"File XYZ non trovato: {filename}")
        
        # Restituisci il file per il download
        return FileResponse(
            file_path,
            media_type="chemical/x-xyz",
            filename=filename,
            headers={"Content-Disposition": f"attachment; filename={filename}"}
        )
        
    except HTTPException as he:
        raise he
    except Exception as e:
        logger.error(f"Errore nel download del file XYZ: {str(e)}")
        raise HTTPException(status_code=500, detail=f"Errore interno: {str(e)}")

@app.get("/api/download-xyz-batch/{batch_folder}")
async def download_xyz_batch_zip(batch_folder: str):
    """
    Endpoint per scaricare tutti i file XYZ di una cartella batch come archivio ZIP.
    """
    try:
        batch_path = os.path.join(MOLECULES_DIR, batch_folder)
        
        if not os.path.exists(batch_path) or not os.path.isdir(batch_path):
            raise HTTPException(status_code=404, detail=f"Cartella batch non trovata: {batch_folder}")
        
        # Trova tutti i file XYZ nella cartella
        xyz_files = [f for f in os.listdir(batch_path) if f.endswith('.xyz')]
        
        if not xyz_files:
            raise HTTPException(status_code=404, detail="Nessun file XYZ trovato nella cartella batch")
        
        # Crea un file ZIP temporaneo
        import tempfile
        with tempfile.NamedTemporaryFile(delete=False, suffix='.zip') as temp_zip:
            temp_zip_path = temp_zip.name
        
        try:
            # Crea il file ZIP
            with zipfile.ZipFile(temp_zip_path, 'w', zipfile.ZIP_DEFLATED) as zip_file:
                for xyz_file in xyz_files:
                    file_path = os.path.join(batch_path, xyz_file)
                    # Aggiungi il file al ZIP mantenendo solo il nome del file
                    zip_file.write(file_path, xyz_file)
            
            # Preparare il nome del file ZIP
            zip_filename = f"{batch_folder}_structures.zip"
            
            # Restituisci il file ZIP come FileResponse
            return FileResponse(
                temp_zip_path,
                media_type="application/zip",
                filename=zip_filename,
                headers={"Content-Disposition": f"attachment; filename={zip_filename}"}
            )
            
        except Exception as e:
            # Pulizia del file temporaneo in caso di errore
            if os.path.exists(temp_zip_path):
                os.unlink(temp_zip_path)
            raise e
        
    except HTTPException as he:
        raise he
    except Exception as e:
        logger.error(f"Errore nel download del batch ZIP: {str(e)}")
        raise HTTPException(status_code=500, detail=f"Errore interno: {str(e)}")
    
@app.get("/api/batch-info/{batch_folder}")
async def get_batch_info(batch_folder: str):
    """
    Endpoint per ottenere informazioni su una cartella batch.
    """
    try:
        batch_path = os.path.join(MOLECULES_DIR, batch_folder)
        
        if not os.path.exists(batch_path) or not os.path.isdir(batch_path):
            raise HTTPException(status_code=404, detail=f"Cartella batch non trovata: {batch_folder}")
        
        # Conta i file XYZ nella cartella
        xyz_files = [f for f in os.listdir(batch_path) if f.endswith('.xyz')]
        
        # Calcola la dimensione totale
        total_size = 0
        for xyz_file in xyz_files:
            file_path = os.path.join(batch_path, xyz_file)
            total_size += os.path.getsize(file_path)
        
        return {
            "batch_folder": batch_folder,
            "file_count": len(xyz_files),
            "total_size_bytes": total_size,
            "total_size_mb": round(total_size / (1024 * 1024), 2),
            "files": xyz_files[:10]  # Mostra solo i primi 10 nomi di file
        }
        
    except HTTPException as he:
        raise he
    except Exception as e:
        logger.error(f"Errore nel recupero info batch: {str(e)}")
        raise HTTPException(status_code=500, detail=f"Errore interno: {str(e)}")

# Endpoint per salute dell'API
@app.get("/api/health")
async def health_check():
    return {"status": "healthy"}

# Avvio dell'applicazione
if __name__ == "__main__":
    uvicorn.run("main:app", host="0.0.0.0", port=8001, reload=True)