import os
import tempfile
import logging
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from PIL import Image
import io

# Configura il logger
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def smiles_to_2d_image(smiles, output_path=None, size=(300, 300)):
    """
    Converte una stringa SMILES in un'immagine 2D e la salva o restituisce come bytes.
    
    Args:
        smiles (str): Stringa SMILES della molecola
        output_path (str, optional): Percorso di output per l'immagine. Se None, l'immagine viene restituita come bytes.
        size (tuple, optional): Dimensioni dell'immagine (larghezza, altezza).
        
    Returns:
        Se output_path è None, restituisce i bytes dell'immagine, altrimenti True se la conversione è avvenuta con successo, False altrimenti.
    """
    try:
        # Conversione da SMILES a molecola RDKit
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            logger.error(f"Impossibile convertire SMILES: {smiles}")
            return False if output_path else None
        
        # Generazione dell'immagine 2D
        img = Draw.MolToImage(mol, size=size)
        
        if output_path:
            # Salva l'immagine su file
            img.save(output_path)
            logger.info(f"Immagine salvata in: {output_path}")
            return True
        else:
            # Converti l'immagine in bytes
            img_bytes = io.BytesIO()
            img.save(img_bytes, format='PNG')
            return img_bytes.getvalue()
            
    except Exception as e:
        logger.error(f"Errore nella generazione dell'immagine 2D: {str(e)}")
        return False if output_path else None

def smiles_to_3d(smiles, output_path):
    """
    Converte una stringa SMILES in un modello 3D e lo salva nel percorso specificato.
    
    Args:
        smiles (str): Stringa SMILES della molecola
        output_path (str): Percorso di output per il file PDB
        
    Returns:
        bool: True se la conversione è avvenuta con successo, False altrimenti
    """
    try:
        logger.info(f"Elaborazione SMILES: {smiles}")
        
        # Conversione da SMILES a molecola RDKit
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            logger.error(f"Impossibile convertire SMILES: {smiles}")
            return False
            
        # Aggiunta degli idrogeni
        mol = Chem.AddHs(mol)
        
        # Generazione della conformazione 3D
        AllChem.EmbedMolecule(mol, AllChem.ETKDG())
        
        # Ottimizzazione della geometria
        AllChem.MMFFOptimizeMolecule(mol)
        
        # Scrittura del file PDB
        Chem.MolToPDBFile(mol, output_path)
        
        logger.info(f"File PDB creato con successo: {output_path}")
        return True
        
    except Exception as e:
        logger.error(f"Errore nella conversione SMILES to 3D: {str(e)}")
        return False