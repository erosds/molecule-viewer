// Aggiorna il file frontend/src/App.js

import React, { useState } from 'react';
import FileSelector from './components/FileSelector';
import MoleculeGrid from './components/MoleculeGrid';
import ValidationButton from './components/ValidationButton';
import CoordinationFilter from './components/CoordinationFilter';
import './App.css';

function App() {
  const [selectedFile, setSelectedFile] = useState(null);
  const [molecules, setMolecules] = useState([]);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);

  const [referenceFile, setReferenceFile] = useState(null);
  const [referenceMolecules, setReferenceMolecules] = useState([]);
  const [newMolecules, setNewMolecules] = useState(0);
  
  // Stato per i risultati della validazione
  const [validationResults, setValidationResults] = useState(null);
  
  // Nuovo stato per il filtro di coordinazione
  const [isFiltered, setIsFiltered] = useState(false);
  const [originalMolecules, setOriginalMolecules] = useState([]);
  const [coordinationStats, setCoordinationStats] = useState(null);
  const [filterInfo, setFilterInfo] = useState(null);

  const handleFileSelect = (file) => {
    setSelectedFile(file);
    setLoading(true);
    setError(null);
    setValidationResults(null);
    setIsFiltered(false); // Reset del filtro quando cambia il file
    setOriginalMolecules([]);
    setFilterInfo(null);

    fetch(`/api/molecules/${file}`)
      .then(response => {
        if (!response.ok) {
          throw new Error('Impossibile caricare le molecole dal file CSV');
        }
        return response.json();
      })
      .then(parsedMolecules => {
        setMolecules(parsedMolecules);
        setOriginalMolecules(parsedMolecules); // Salva le molecole originali
        setLoading(false);

        // Calcola molecole nuove se il file di riferimento è già caricato
        if (referenceMolecules.length > 0) {
          calculateNewMolecules(parsedMolecules, referenceMolecules);
        } else {
          setNewMolecules(parsedMolecules.length);
        }
      })
      .catch(err => {
        setError(`Errore: ${err.message}`);
        setLoading(false);
      });
  };

  const handleReferenceFileSelect = (file) => {
    setReferenceFile(file);

    if (!file) {
      setReferenceMolecules([]);
      return;
    }

    fetch(`/api/molecules/${file}`)
      .then(response => {
        if (!response.ok) {
          throw new Error('Impossibile caricare le molecole dal file di riferimento');
        }
        return response.json();
      })
      .then(parsedMolecules => {
        setReferenceMolecules(parsedMolecules);

        // Calcola molecole nuove se entrambi i set sono disponibili
        if (molecules.length > 0) {
          calculateNewMolecules(molecules, parsedMolecules);
        }
      })
      .catch(err => {
        console.error(`Errore: ${err.message}`);
        setReferenceMolecules([]);
      });
  };

  const calculateNewMolecules = (mainMolecules, refMolecules) => {
    if (!refMolecules.length) {
      setNewMolecules(mainMolecules.length);
      return;
    }

    const referenceSmiles = new Set(refMolecules.map(mol => mol.smiles));
    const newCount = mainMolecules.filter(mol => !referenceSmiles.has(mol.smiles)).length;
    setNewMolecules(newCount);
  };

  // Handler per i risultati della validazione
  const handleValidationComplete = (results) => {
    setValidationResults(results);
  };

  // Handler per il filtro di coordinazione
  const handleFilterApplied = (filteredMolecules, filterResult) => {
    setMolecules(filteredMolecules);
    setIsFiltered(true);
    setFilterInfo(filterResult);
    setValidationResults(null); // Reset validation quando si applica un filtro
    
    // Ricalcola le molecole nuove per il subset filtrato
    if (referenceMolecules.length > 0) {
      calculateNewMolecules(filteredMolecules, referenceMolecules);
    } else {
      setNewMolecules(filteredMolecules.length);
    }
  };

  // Handler per aggiornare le statistiche di coordinazione
  const handleStatsUpdate = (stats) => {
    setCoordinationStats(stats);
  };

  // Funzione per resettare il filtro
  const resetFilter = () => {
    setMolecules(originalMolecules);
    setIsFiltered(false);
    setFilterInfo(null);
    setValidationResults(null);
    
    // Ricalcola le molecole nuove per il set completo
    if (referenceMolecules.length > 0) {
      calculateNewMolecules(originalMolecules, referenceMolecules);
    } else {
      setNewMolecules(originalMolecules.length);
    }
  };

  return (
    <div className="app">
      <header className="app-header">
        <h1>Visualizzatore Molecole</h1>
      </header>
      <main className="app-content">
        {/* Selettore del file di riferimento */}
        <FileSelector 
          onSelectFile={handleReferenceFileSelect} 
          selectedFile={referenceFile} 
          type="reference"
          label="Seleziona il file CSV di molecole di riferimento:"
        />
        
        {/* Selettore del file principale */}
        <FileSelector 
          onSelectFile={handleFileSelect} 
          selectedFile={selectedFile} 
          type="main"
          label="Seleziona un file CSV di molecole da visualizzare:"
        />
        
        {/* Filtro di coordinazione */}
        {selectedFile && originalMolecules.length > 0 && (
          <CoordinationFilter 
            selectedFile={selectedFile}
            onFilterApplied={handleFilterApplied}
            onStatsUpdate={handleStatsUpdate}
          />
        )}
        
        {/* Informazioni sul filtro applicato */}
        {isFiltered && filterInfo && (
          <div className="filter-info">
            <div className="filter-summary">
              <span className="filter-label">
                🔍 Filtro attivo: coordinazione {filterInfo.coordination_distribution ? 
                  `${Math.min(...Object.keys(filterInfo.coordination_distribution).map(Number))}-${Math.max(...Object.keys(filterInfo.coordination_distribution).map(Number))}` : 
                  'personalizzata'}
              </span>
              <span className="filter-stats">
                {filterInfo.filtered_molecules} di {filterInfo.total_molecules} molecole
              </span>
              <button 
                className="reset-filter-button"
                onClick={resetFilter}
                title="Rimuovi filtro e mostra tutte le molecole"
              >
                ✕ Rimuovi filtro
              </button>
            </div>
          </div>
        )}
        
        {loading ? (
          <div className="loading">Caricamento molecole...</div>
        ) : error ? (
          <div className="error">{error}</div>
        ) : molecules.length > 0 ? (
          <>
            {/* Pulsante di validazione delle strutture 3D */}
            <ValidationButton 
              molecules={molecules} 
              onValidationComplete={handleValidationComplete}
            />
            
            <MoleculeGrid 
              molecules={molecules} 
              newMoleculesCount={newMolecules}
              validationResults={validationResults}
            />
          </>
        ) : selectedFile ? (
          <div className="error">
            {isFiltered ? 
              "Nessuna molecola corrisponde ai criteri di filtro selezionati" :
              "Nessuna molecola SMILES riconosciuta nel file selezionato"
            }
          </div>
        ) : (
          <div className="instructions">
            Seleziona un file CSV contenente strutture molecolari SMILES per iniziare
          </div>
        )}
      </main>
    </div>
  );
}

export default App;