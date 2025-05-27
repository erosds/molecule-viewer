// App.js - Layout con Sidebar

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
    setIsFiltered(false);
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
        setOriginalMolecules(parsedMolecules);
        setLoading(false);

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

  const handleValidationComplete = (results) => {
    setValidationResults(results);
  };

  const handleFilterApplied = (filteredMolecules, filterResult) => {
    setMolecules(filteredMolecules);
    setIsFiltered(true);
    setFilterInfo(filterResult);
    setValidationResults(null);

    if (referenceMolecules.length > 0) {
      calculateNewMolecules(filteredMolecules, referenceMolecules);
    } else {
      setNewMolecules(filteredMolecules.length);
    }
  };

  const handleStatsUpdate = (stats) => {
    setCoordinationStats(stats);
  };

  const resetFilter = () => {
    setMolecules(originalMolecules);
    setIsFiltered(false);
    setFilterInfo(null);
    setValidationResults(null);

    if (referenceMolecules.length > 0) {
      calculateNewMolecules(originalMolecules, referenceMolecules);
    } else {
      setNewMolecules(originalMolecules.length);
    }
  };

  return (
    <div className="app">
      {/* Header */}
      <header className="app-header">
        <h1>Visualizzatore Molecole</h1>
      </header>

      {/* Sidebar */}
      <aside className="app-sidebar">
        <div className="sidebar-content">
          {/* Sezione File di Riferimento */}
          <div className="sidebar-section">
            <div className="sidebar-section-title">
              <span className="icon">📁</span>
              File di Riferimento
            </div>
            <FileSelector
              onSelectFile={handleReferenceFileSelect}
              selectedFile={referenceFile}
              type="reference"
              label="Molecole di riferimento:"
            />
          </div>

          {/* Sezione File Principale */}
          <div className="sidebar-section">
            <div className="sidebar-section-title">
              <span className="icon">📊</span>
              File da Visualizzare
            </div>
            <FileSelector
              onSelectFile={handleFileSelect}
              selectedFile={selectedFile}
              type="main"
              label="File CSV molecole:"
            />
          </div>

          {/* Sezione Filtro di Coordinazione */}
          {selectedFile && originalMolecules.length > 0 && (
            <div className="sidebar-section">
              <CoordinationFilter
                selectedFile={selectedFile}
                onFilterApplied={handleFilterApplied}
                onStatsUpdate={handleStatsUpdate}
              />
            </div>
          )}

          {/* Informazioni Filtro Attivo */}
          {isFiltered && filterInfo && (
            <div className="sidebar-section">
              <div className="sidebar-section-title">
                <span className="icon">🔍</span>
                Filtro Attivo
              </div>
              <div className="filter-info">
                <div className="filter-summary">
                  <div className="stat-item filtered">
                    <div className="stat-number">{filterInfo.filtered_molecules}</div>
                    <div className="stat-label">Molecole Filtrate</div>
                  </div>
                  <div className="filter-details">
                    <span className="filter-stats">
                      su {filterInfo.total_molecules} totali
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
              </div>
            </div>
          )}

          {/* Sezione Validazione */}
          {molecules.length > 0 && (
            <div className="sidebar-section">
              <div className="sidebar-section-title">
                <span className="icon">🧪</span>
                Validazione 3D
              </div>
              <ValidationButton
                molecules={molecules}
                onValidationComplete={handleValidationComplete}
              />
            </div>
          )}
        </div>
      </aside>

      {/* Contenuto Principale */}
      <main className="app-content">
        {loading ? (
          <div className="loading">Caricamento molecole...</div>
        ) : error ? (
          <div className="error">{error}</div>
        ) : molecules.length > 0 ? (
          <MoleculeGrid
            molecules={molecules}
            newMoleculesCount={newMolecules}
            validationResults={validationResults}
          />
        ) : selectedFile ? (
          <div className="error">
            {isFiltered ?
              "Nessuna molecola corrisponde ai criteri di filtro selezionati" :
              "Nessuna molecola SMILES riconosciuta nel file selezionato"
            }
          </div>
        ) : (
          <div className="instructions">
            Seleziona un file CSV contenente strutture molecolari SMILES dalla sidebar per iniziare
          </div>
        )}
      </main>
    </div>
  );
}

export default App;