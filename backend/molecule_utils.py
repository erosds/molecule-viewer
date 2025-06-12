import os
import tempfile
import logging
from pathlib import Path
from typing import Dict, List, Optional, Union, Tuple, Any
from rdkit import Chem
from rdkit.Chem import Draw, AllChem, Descriptors
from rdkit.Chem.Draw import rdMolDraw2D
from PIL import Image
import io

# Disabilita i warning di RDKit
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

# Configura il logger
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Costanti
METAL_ATOMS = frozenset({
    # Metalli di transizione comuni
    'Fe', 'Cu', 'Zn', 'Ni', 'Co', 'Mn', 'Cr', 'V', 'Ti', 'Sc',
    'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'Hf', 'Ta', 'W', 'Re', 'Os',
    'Ir', 'Pt', 'Au', 'Hg', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt',
    'Ds', 'Rg', 'Cn',
    # Lantanidi e Attinidi
    'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy',
    'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Ac', 'Th', 'Pa', 'U', 'Np',
    'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr',
    # Altri metalli
    'Al', 'Ga', 'In', 'Sn', 'Pb', 'Bi', 'Mg', 'Ca', 'Sr', 'Ba',
    'Li', 'Na', 'K', 'Rb', 'Cs'
})

DEFAULT_IMAGE_SIZE = (300, 300)
MAX_COORDINATION_NUMBER = 12


def get_metal_atoms() -> frozenset:
    """
    Restituisce un frozenset degli elementi metallici più comuni nei complessi.
    
    Returns:
        frozenset: Set immutabile dei simboli metallici
    """
    return METAL_ATOMS


def _validate_smiles(smiles: str) -> bool:
    """
    Valida se una stringa SMILES è valida.
    
    Args:
        smiles: Stringa SMILES da validare
        
    Returns:
        bool: True se valida, False altrimenti
    """
    return bool(smiles and isinstance(smiles, str) and smiles.strip())


def _create_mol_from_smiles(smiles: str) -> Optional[Chem.Mol]:
    """
    Crea una molecola RDKit da una stringa SMILES.
    
    Args:
        smiles: Stringa SMILES
        
    Returns:
        Molecola RDKit o None se conversione fallita
    """
    if not _validate_smiles(smiles):
        logger.warning(f"SMILES non valido: {smiles}")
        return None
        
    mol = Chem.MolFromSmiles(smiles.strip())
    if mol is None:
        logger.warning(f"Impossibile convertire SMILES: {smiles}")
    return mol


def analyze_metal_coordination(smiles: str) -> Dict[str, Any]:
    """
    Analizza il numero di legami dell'atomo metallico centrale in una molecola SMILES.
    
    Args:
        smiles: Stringa SMILES della molecola
        
    Returns:
        dict: Dizionario con informazioni sui metalli e i loro legami
    """
    default_result = {
        'has_metal': False,
        'metal_atoms': [],
        'max_coordination': 0,
        'min_coordination': 0,
        'total_metals': 0
    }
    
    try:
        mol = _create_mol_from_smiles(smiles)
        if mol is None:
            return default_result
        
        metal_atoms = []
        
        # Analizza ogni atomo nella molecola
        for atom in mol.GetAtoms():
            symbol = atom.GetSymbol()
            if symbol in METAL_ATOMS:
                coordination_number = atom.GetDegree()  # Numero di legami pesanti (non-H)
                
                metal_info = {
                    'symbol': symbol,
                    'atom_idx': atom.GetIdx(),
                    'coordination_number': coordination_number,
                    'formal_charge': atom.GetFormalCharge(),
                    'total_valence': atom.GetTotalValence()
                }
                metal_atoms.append(metal_info)
        
        # Calcola statistiche
        if metal_atoms:
            coordination_numbers = [m['coordination_number'] for m in metal_atoms]
            return {
                'has_metal': True,
                'metal_atoms': metal_atoms,
                'max_coordination': max(coordination_numbers),
                'min_coordination': min(coordination_numbers),
                'total_metals': len(metal_atoms)
            }
        
        return default_result
        
    except Exception as e:
        logger.error(f"Errore nell'analisi della coordinazione metallica per {smiles}: {e}")
        result = default_result.copy()
        result['error'] = str(e)
        return result


def filter_molecules_by_coordination(
    molecules_data: List[Dict[str, Any]], 
    min_coordination: int = 0, 
    max_coordination: int = MAX_COORDINATION_NUMBER
) -> List[Dict[str, Any]]:
    """
    Filtra una lista di molecole in base al numero di coordinazione dell'atomo metallico.
    
    Args:
        molecules_data: Lista di dizionari con dati delle molecole
        min_coordination: Numero minimo di legami
        max_coordination: Numero massimo di legami
        
    Returns:
        list: Lista filtrata di molecole con info aggiuntive sui metalli
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
        
        if (coordination_info['has_metal'] and 
            min_coordination <= coordination_info['max_coordination'] <= max_coordination):
            
            molecule_copy = molecule.copy()
            molecule_copy['coordination_info'] = coordination_info
            filtered_molecules.append(molecule_copy)
    
    return filtered_molecules


def filter_molecules_by_metal_and_coordination(
    molecules_data: List[Dict[str, Any]], 
    min_coordination: int = 0, 
    max_coordination: int = MAX_COORDINATION_NUMBER, 
    selected_metals: Optional[List[str]] = None
) -> List[Dict[str, Any]]:
    """
    Filtra una lista di molecole in base al numero di coordinazione metallica e al tipo di metallo.
    
    Args:
        molecules_data: Lista di dizionari con dati delle molecole
        min_coordination: Numero minimo di legami
        max_coordination: Numero massimo di legami
        selected_metals: Lista di simboli metallici da includere (es. ['Fe', 'Cu', 'Zn'])
        
    Returns:
        list: Lista filtrata di molecole con info aggiuntive sui metalli
    """
    if not isinstance(molecules_data, list):
        logger.error("molecules_data deve essere una lista")
        return []
    
    # Normalizza selected_metals
    if selected_metals:
        selected_metals = set(selected_metals)
    
    filtered_molecules = []
    
    for molecule in molecules_data:
        if not isinstance(molecule, dict):
            continue
            
        smiles = molecule.get('smiles', '')
        if not smiles:
            continue
            
        coordination_info = analyze_metal_coordination(smiles)
        
        if not coordination_info['has_metal']:
            continue
            
        max_coord = coordination_info['max_coordination']
        if not (min_coordination <= max_coord <= max_coordination):
            continue
            
        # Verifica filtro per metalli specifici
        if selected_metals:
            molecule_metals = {metal['symbol'] for metal in coordination_info['metal_atoms']}
            if not molecule_metals.intersection(selected_metals):
                continue
        
        molecule_copy = molecule.copy()
        molecule_copy['coordination_info'] = coordination_info
        filtered_molecules.append(molecule_copy)
    
    return filtered_molecules


def _create_image_with_draw(mol: Chem.Mol, size: Tuple[int, int]) -> Optional[bytes]:
    """
    Crea un'immagine usando Draw.MolToImage.
    
    Args:
        mol: Molecola RDKit
        size: Dimensioni dell'immagine
        
    Returns:
        bytes dell'immagine o None se fallisce
    """
    try:
        img = Draw.MolToImage(mol, size=size)
        img_bytes = io.BytesIO()
        img.save(img_bytes, format='PNG')
        return img_bytes.getvalue()
    except Exception as e:
        logger.warning(f"Draw.MolToImage fallito: {e}")
        return None


def _create_image_with_drawer(mol: Chem.Mol, size: Tuple[int, int]) -> Optional[bytes]:
    """
    Crea un'immagine usando rdMolDraw2D.
    
    Args:
        mol: Molecola RDKit
        size: Dimensioni dell'immagine
        
    Returns:
        bytes dell'immagine o None se fallisce
    """
    try:
        drawer = rdMolDraw2D.MolDraw2DCairo(size[0], size[1])
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        return drawer.GetDrawingText()
    except Exception as e:
        logger.error(f"rdMolDraw2D fallito: {e}")
        return None


def smiles_to_2d_image(
    smiles: str, 
    output_path: Optional[str] = None, 
    size: Tuple[int, int] = DEFAULT_IMAGE_SIZE
) -> Union[bool, bytes, None]:
    """
    Converte una stringa SMILES in un'immagine 2D e la salva o restituisce come bytes.
    
    Args:
        smiles: Stringa SMILES della molecola
        output_path: Percorso di output per l'immagine. Se None, l'immagine viene restituita come bytes
        size: Dimensioni dell'immagine (larghezza, altezza)
        
    Returns:
        Se output_path è None, restituisce i bytes dell'immagine, 
        altrimenti True se la conversione è avvenuta con successo, False altrimenti
    """
    try:
        mol = _create_mol_from_smiles(smiles)
        if mol is None:
            return False if output_path else None
        
        # Prova prima con Draw.MolToImage
        img_bytes = _create_image_with_draw(mol, size)
        
        # Se fallisce, prova con rdMolDraw2D
        if img_bytes is None:
            img_bytes = _create_image_with_drawer(mol, size)
            
        if img_bytes is None:
            logger.error(f"Impossibile generare immagine per SMILES: {smiles}")
            return False if output_path else None
        
        if output_path:
            with open(output_path, 'wb') as f:
                f.write(img_bytes)
            logger.info(f"Immagine salvata in: {output_path}")
            return True
        else:
            return img_bytes
            
    except Exception as e:
        logger.error(f"Errore nella generazione dell'immagine 2D: {e}")
        return False if output_path else None


def _write_xyz_file(mol: Chem.Mol, output_path: str, smiles: str) -> None:
    """
    Scrive un file XYZ da una molecola RDKit.
    
    Args:
        mol: Molecola RDKit con conformazione 3D
        output_path: Percorso del file di output
        smiles: SMILES originale per il commento
    """
    conf = mol.GetConformer()
    atoms = []
    
    for i in range(mol.GetNumAtoms()):
        atom = mol.GetAtomWithIdx(i)
        symbol = atom.GetSymbol()
        pos = conf.GetAtomPosition(i)
        atoms.append((symbol, pos.x, pos.y, pos.z))
    
    with open(output_path, 'w') as f:
        f.write(f"{len(atoms)}\n")
        f.write(f"Molecule generated from SMILES: {smiles}\n")
        for symbol, x, y, z in atoms:
            f.write(f"{symbol} {x:.6f} {y:.6f} {z:.6f}\n")


def smiles_to_3d(smiles: str, output_path: str) -> bool:
    """
    Converte una stringa SMILES in un modello 3D e lo salva nel percorso specificato.
    
    Args:
        smiles: Stringa SMILES della molecola
        output_path: Percorso di output per il file
        
    Returns:
        bool: True se la conversione è avvenuta con successo, False altrimenti
    """
    try:
        logger.info(f"Elaborazione SMILES: {smiles}")
        
        mol = _create_mol_from_smiles(smiles)
        if mol is None:
            return False
            
        # Aggiunta degli idrogeni
        mol = Chem.AddHs(mol)
        
        # Generazione della conformazione 3D
        if AllChem.EmbedMolecule(mol, AllChem.ETKDG()) != 0:
            logger.error(f"Impossibile generare conformazione 3D per: {smiles}")
            return False
        
        # Ottimizzazione della geometria
        AllChem.MMFFOptimizeMolecule(mol)
        
        # Determina il formato del file dall'estensione
        if output_path.lower().endswith('.pdb'):
            Chem.MolToPDBFile(mol, output_path)
            logger.info(f"File PDB creato con successo: {output_path}")
        else:
            _write_xyz_file(mol, output_path, smiles)
            logger.info(f"File XYZ creato con successo: {output_path}")
        
        return True
        
    except Exception as e:
        logger.error(f"Errore nella conversione SMILES to 3D: {e}")
        return False