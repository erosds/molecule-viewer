import os
import tempfile
import logging
from pathlib import Path
import math
import random

# Configura il logger
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def smiles_to_3d(smiles, output_path):
    """
    Crea un file PDB fittizio per la visualizzazione.
    
    Questo è un sostituto semplificato che non richiede RDKit.
    In una implementazione reale, usare RDKit per una corretta generazione 3D.
    
    Args:
        smiles (str): Stringa SMILES della molecola
        output_path (str): Percorso di output per il file PDB
        
    Returns:
        bool: True se la conversione è avvenuta con successo, False altrimenti
    """
    try:
        logger.info(f"Elaborazione SMILES: {smiles}")
        
        # Identificazione di alcuni elementi basilari dalla stringa SMILES
        atoms = []
        # Elementi più comuni in chimica organica
        elements = {'C': 0, 'H': 0, 'O': 0, 'N': 0, 'P': 0, 'S': 0}
        
        # Conteggio molto approssimativo degli elementi (non accurato!)
        for char in smiles:
            if char in elements:
                elements[char] += 1
        
        # Aggiungi almeno un atomo per ogni elemento presente
        atom_id = 1
        center_x, center_y, center_z = 0, 0, 0
        
        with open(output_path, "w") as f:
            f.write("HEADER    MOCK MOLECULE FROM SMILES\n")
            f.write(f"TITLE     Generated from SMILES: {smiles}\n")
            
            # Generiamo atomi fittizi per visualizzazione
            radius = 3.0  # Raggio approssimativo della molecola
            
            for element, count in elements.items():
                if count > 0:
                    # Genera posizioni per questo elemento
                    for i in range(count):
                        # Posizione approssimativa basata su un pattern circolare
                        angle = 2 * math.pi * i / max(count, 1)
                        x = center_x + radius * math.cos(angle)
                        y = center_y + radius * math.sin(angle)
                        z = center_z + random.uniform(-1.0, 1.0)
                        
                        # Aggiungi l'atomo al file PDB
                        f.write(f"ATOM  {atom_id:5d}  {element}   UNK     1    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element}  \n")
                        atoms.append((atom_id, element, x, y, z))
                        atom_id += 1
            
            # Aggiungi connessioni semplici tra atomi vicini
            for i, (id1, el1, x1, y1, z1) in enumerate(atoms):
                f.write(f"CONECT{id1:5d}")
                # Connetti a massimo 4 atomi vicini
                connections = 0
                for j, (id2, el2, x2, y2, z2) in enumerate(atoms):
                    if i != j and connections < 4:
                        # Calcola la distanza
                        dist = math.sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)
                        if dist < radius:  # Se abbastanza vicino
                            f.write(f"{id2:5d}")
                            connections += 1
                f.write("\n")
            
            f.write("END\n")
        
        logger.info(f"File PDB fittizio creato con successo: {output_path}")
        return True
        
    except Exception as e:
        logger.error(f"Errore nella generazione del file PDB: {str(e)}")
        return False

def smiles_to_svg(smiles, output_path=None):
    """
    Crea un SVG semplice per la molecola.
    
    Args:
        smiles (str): Stringa SMILES della molecola
        output_path (str, optional): Percorso di output per il file SVG
        
    Returns:
        str: SVG della molecola o percorso del file se output_path è specificato
    """
    try:
        # Crea un SVG semplice
        svg = f'''<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" width="200" height="200" viewBox="0 0 200 200">
  <circle cx="100" cy="100" r="50" fill="#f0f0f0" stroke="#555555" stroke-width="2" />
  <text x="100" y="105" font-family="Arial" font-size="12" text-anchor="middle">{smiles[:20]}</text>
</svg>'''
        
        # Se è specificato un percorso di output, salva il file
        if output_path:
            with open(output_path, 'w') as f:
                f.write(svg)
            return output_path
        else:
            return svg
            
    except Exception as e:
        logger.error(f"Errore nella generazione SVG: {str(e)}")
        return None

def get_molecule_properties(smiles):
    """
    Genera proprietà fittizie di una molecola.
    
    Args:
        smiles (str): Stringa SMILES della molecola
        
    Returns:
        dict: Dizionario con le proprietà molecolari
    """
    properties = {}
    
    try:
        # Calcolo basilare e approssimativo (non accurato!)
        # Conta elementi (molto approssimativo)
        elements = {'C': 0, 'H': 0, 'O': 0, 'N': 0, 'P': 0, 'S': 0}
        formula_parts = []
        
        for char in smiles:
            if char in elements:
                elements[char] += 1
        
        # Formula approssimativa
        for element, count in elements.items():
            if count > 0:
                formula_parts.append(f"{element}{count if count > 1 else ''}")
        
        properties["formula"] = "".join(formula_parts)
        
        # Peso molecolare approssimativo
        weights = {'C': 12.01, 'H': 1.01, 'O': 16.00, 'N': 14.01, 'P': 30.97, 'S': 32.07}
        mol_weight = sum(count * weights[element] for element, count in elements.items() if count > 0)
        
        properties["exact_mass"] = mol_weight
        properties["mol_weight"] = mol_weight
        properties["heavy_atoms"] = sum(count for element, count in elements.items() if element != 'H' and count > 0)
        properties["rings"] = smiles.count('1') # Approssimativo, non accurato!
        properties["aromatic_rings"] = smiles.count('c') // 6  # Molto approssimativo!
        properties["rotatable_bonds"] = len(smiles) // 5  # Totalmente fittizio
        properties["h_donors"] = elements['O'] + elements['N']  # Approssimativo
        properties["h_acceptors"] = elements['O'] + elements['N']  # Approssimativo
        properties["tpsa"] = (elements['O'] + elements['N']) * 20  # Totalmente fittizio
        properties["logp"] = elements['C'] * 0.5 - elements['O'] * 0.5  # Totalmente fittizio
        
    except Exception as e:
        logger.error(f"Errore nel calcolo delle proprietà molecolari: {str(e)}")
        
    return properties