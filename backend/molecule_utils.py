import os
import tempfile
import logging
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import Draw, AllChem, Descriptors
from PIL import Image
import io

# Configura il logger
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def get_metal_atoms():
    """
    Restituisce un set degli elementi metallici più comuni nei complessi.
    """
    return {
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
    }

def analyze_metal_coordination(smiles: str):
    """
    Analizza il numero di legami dell'atomo metallico centrale in una molecola SMILES.
    
    Args:
        smiles (str): Stringa SMILES della molecola
        
    Returns:
        dict: Dizionario con informazioni sui metalli e i loro legami
        {
            'has_metal': bool,
            'metal_atoms': list,  # Lista di dizionari con info sui metalli
            'max_coordination': int,  # Numero massimo di legami di un metallo
            'min_coordination': int,  # Numero minimo di legami di un metallo
            'total_metals': int
        }
    """
    try:
        # Conversione da SMILES a molecola RDKit
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            logger.warning(f"Impossibile convertire SMILES: {smiles}")
            return {
                'has_metal': False,
                'metal_atoms': [],
                'max_coordination': 0,
                'min_coordination': 0,
                'total_metals': 0
            }
        
        metal_elements = get_metal_atoms()
        metal_atoms = []
        
        # Analizza ogni atomo nella molecola
        for atom in mol.GetAtoms():
            symbol = atom.GetSymbol()
            if symbol in metal_elements:
                # Calcola il numero di legami
                degree = atom.GetDegree()  # Numero di legami pesanti (non-H)
                total_valence = atom.GetTotalValence()  # Include anche idrogeni impliciti
                
                # Per i metalli, di solito il grado è più significativo
                coordination_number = degree
                
                metal_info = {
                    'symbol': symbol,
                    'atom_idx': atom.GetIdx(),
                    'coordination_number': coordination_number,
                    'formal_charge': atom.GetFormalCharge(),
                    'total_valence': total_valence
                }
                metal_atoms.append(metal_info)
        
        # Calcola statistiche
        if metal_atoms:
            coordination_numbers = [m['coordination_number'] for m in metal_atoms]
            result = {
                'has_metal': True,
                'metal_atoms': metal_atoms,
                'max_coordination': max(coordination_numbers),
                'min_coordination': min(coordination_numbers),
                'total_metals': len(metal_atoms)
            }
        else:
            result = {
                'has_metal': False,
                'metal_atoms': [],
                'max_coordination': 0,
                'min_coordination': 0,
                'total_metals': 0
            }
            
        return result
        
    except Exception as e:
        logger.error(f"Errore nell'analisi della coordinazione metallica per {smiles}: {str(e)}")
        return {
            'has_metal': False,
            'metal_atoms': [],
            'max_coordination': 0,
            'min_coordination': 0,
            'total_metals': 0,
            'error': str(e)
        }

def filter_molecules_by_coordination(molecules_data: list, min_coordination: int = 0, max_coordination: int = 12):
    """
    Filtra una lista di molecole in base al numero di coordinazione dell'atomo metallico.
    
    Args:
        molecules_data (list): Lista di dizionari con dati delle molecole
        min_coordination (int): Numero minimo di legami
        max_coordination (int): Numero massimo di legami
        
    Returns:
        list: Lista filtrata di molecole con info aggiuntive sui metalli
    """
    filtered_molecules = []
    
    for molecule in molecules_data:
        smiles = molecule.get('smiles', '')
        if not smiles:
            continue
            
        # Analizza la coordinazione metallica
        coordination_info = analyze_metal_coordination(smiles)
        
        # Verifica se la molecola soddisfa i criteri di filtro
        if coordination_info['has_metal']:
            max_coord = coordination_info['max_coordination']
            if min_coordination <= max_coord <= max_coordination:
                # Aggiungi le informazioni sulla coordinazione alla molecola
                molecule_copy = molecule.copy()
                molecule_copy['coordination_info'] = coordination_info
                filtered_molecules.append(molecule_copy)
    
    return filtered_molecules

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
    
def smiles_to_3d(smiles, output_path):
    """
    Converte una stringa SMILES in un modello 3D e lo salva nel percorso specificato.
    
    Args:
        smiles (str): Stringa SMILES della molecola
        output_path (str): Percorso di output per il file
        
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
        
        # Verifica se il percorso termina con .pdb o .xyz per determinare il formato
        if output_path.lower().endswith('.pdb'):
            # Scrittura del file PDB
            Chem.MolToPDBFile(mol, output_path)
            logger.info(f"File PDB creato con successo: {output_path}")
        else:
            # Creiamo un file XYZ (formato più semplice senza numerazione atomi)
            # Estrai le coordinate e i simboli atomici
            conf = mol.GetConformer()
            atoms = []
            for i in range(mol.GetNumAtoms()):
                atom = mol.GetAtomWithIdx(i)
                symbol = atom.GetSymbol()
                pos = conf.GetAtomPosition(i)
                atoms.append((symbol, pos.x, pos.y, pos.z))
            
            # Scriviamo il file XYZ
            # Scriviamo il file XYZ
            with open(output_path, 'w') as f:
                f.write(f"{len(atoms)}\n")  # Prima riga: numero di atomi
                f.write(f"Molecule generated from SMILES: {smiles}\n")  # Seconda riga: commento
                # Righe successive: Simbolo X Y Z
                for symbol, x, y, z in atoms:
                    f.write(f"{symbol} {x:.6f} {y:.6f} {z:.6f}\n")
            
            logger.info(f"File XYZ creato con successo: {output_path}")
        
        return True
        
    except Exception as e:
        logger.error(f"Errore nella conversione SMILES to 3D: {str(e)}")
        return False