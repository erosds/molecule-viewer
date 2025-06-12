// App.js - Layout con Sidebar e Analisi Novelty/Uniqueness

import React, { useState } from 'react';
import FileSelector from './components/FileSelector';
import MoleculeGrid from './components/MoleculeGrid';
import ValidationButton from './components/ValidationButton';
import CoordinationFilter from './components/CoordinationFilter';
import './App.css';

function App() {
  const [selectedFile, setSelectedFile] = useState(null);
  const [referenceFile, setReferenceFile] = useState(null);
  
  // Stati per l'analisi di novelty
  const [isAnalyzing, setIsAnalyzing] = useState(false);
  const [analysisComplete, setAnalysisComplete] = useState(false);
  const [analysisResults, setAnalysisResults] = useState(null);
  const [analysisError, setAnalysisError] = useState(null);
  
  // Stati per le molecole e filtri
  const [molecules, setMolecules] = useState([]);
  const [originalMolecules, setOriginalMolecules] = useState([]);
  const [validationResults, setValidationResults] = useState(null);
  
  // Nuovi stati per il filtro di coordinazione
  const [isFiltered, setIsFiltered] = useState(false);
  const [coordinationStats, setCoordinationStats] = useState(null);
  const [filterInfo, setFilterInfo] = useState(null);

  const handleFileSelect = async (file) => {
    if (!file) {
      resetAnalysis();
      return;
    }

    setSelectedFile(file);
    
    // Se abbiamo entrambi i file o solo il principale, avvia l'analisi
    if (file) {
      await performNoveltyAnalysis(file, referenceFile);
    }
  };

  const handleReferenceFileSelect = async (file) => {
    setReferenceFile(file);
    
    // Se abbiamo già un file principale, riavvia l'analisi
    if (selectedFile) {
      await performNoveltyAnalysis(selectedFile, file);
    }
  };

  const performNoveltyAnalysis = async (mainFile, refFile) => {
    setIsAnalyzing(true);
    setAnalysisComplete(false);
    setAnalysisError(null);
    setValidationResults(null);
    setIsFiltered(false);
    setFilterInfo(null);

    try {
      const response = await fetch('/api/analyze-novelty', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify({
          main_file: mainFile,
          reference_file: refFile
        }),
      });

      if (!response.ok) {
        const errorData = await response.json();
        throw new Error(errorData.detail || 'Errore nell\'analisi di novelty');
      }

      const results = await response.json();
      setAnalysisResults(results);
      
      // Filtra solo le molecole valide per la visualizzazione
      const validMolecules = results.molecules.filter(mol => mol.is_valid);
      setMolecules(validMolecules);
      setOriginalMolecules(validMolecules);
      
      setAnalysisComplete(true);
      
    } catch (err) {
      setAnalysisError(`Errore: ${err.message}`);
      setAnalysisComplete(false);
    } finally {
      setIsAnalyzing(false);
    }
  };

  const resetAnalysis = () => {
    setIsAnalyzing(false);
    setAnalysisComplete(false);
    setAnalysisResults(null);
    setAnalysisError(null);
    setMolecules([]);
    setOriginalMolecules([]);
    setValidationResults(null);
    setIsFiltered(false);
    setFilterInfo(null);
  };

  const handleValidationComplete = (results) => {
    setValidationResults(results);
  };

  const handleFilterApplied = (filteredMolecules, filterResult) => {
    // Mantieni le proprietà di novelty dalle molecole originali
    const moleculesWithNovelty = filteredMolecules.map(filteredMol => {
      const originalMol = originalMolecules.find(origMol => origMol.smiles === filteredMol.smiles);
      return {
        ...filteredMol,
        is_unique: originalMol?.is_unique || false,
        is_novel: originalMol?.is_novel || false,
        canonical_smiles: originalMol?.canonical_smiles || filteredMol.smiles
      };
    });
    
    setMolecules(moleculesWithNovelty);
    setIsFiltered(true);
    setFilterInfo(filterResult);
    setValidationResults(null);
  };

  const handleStatsUpdate = (stats) => {
    setCoordinationStats(stats);
  };

  const resetFilter = () => {
    setMolecules(originalMolecules);
    setIsFiltered(false);
    setFilterInfo(null);
    setValidationResults(null);
  };

  // Calcola le statistiche per la visualizzazione
  const getDisplayStats = () => {
    if (!analysisResults) return null;
    
    const currentMolecules = isFiltered ? molecules : originalMolecules;
    const uniqueMolecules = currentMolecules.filter(mol => mol.is_unique).length;
    const novelMolecules = currentMolecules.filter(mol => mol.is_novel).length;
    
    return {
      total: currentMolecules.length,
      unique: uniqueMolecules,
      novel: novelMolecules,
      hasReference: !!referenceFile
    };
  };

  return (
    <div className="app">
      {/* Header */}
      <header className="app-header">
        <h1>Visualizzatore Molecole Generate</h1>
      </header>

      {/* Sidebar */}
      <aside className="app-sidebar">
        <div className="sidebar-content">
          {/* Sezione File di Riferimento */}
          <div className="sidebar-section">
            <div className="sidebar-section-title">
              <span className="icon">📁</span>
              File di Riferimento (Opzionale)
            </div>
            <FileSelector
              onSelectFile={handleReferenceFileSelect}
              selectedFile={referenceFile}
              type="reference"
            />
          </div>

          {/* Sezione File Principale */}
          <div className="sidebar-section">
            <div className="sidebar-section-title">
              <span className="icon">📊</span>
              File da Analizzare
            </div>
            <FileSelector
              onSelectFile={handleFileSelect}
              selectedFile={selectedFile}
              type="main"
            />
          </div>

          {/* Statistiche di Novelty */}
          {analysisResults && analysisComplete && (
            <div className="sidebar-section">
              <div className="sidebar-section-title">
                <span className="icon">📈</span>
                Analisi Unicità e Novelty
              </div>
              <div className="novelty-stats">
                <div className="stats-grid">
                  <div className="stat-item total">
                    <div className="stat-number">{analysisResults.total_molecules}</div>
                    <div className="stat-label">Molecole Totali</div>
                  </div>
                  <div className="stat-item valid">
                    <div className="stat-number">{analysisResults.valid_molecules}</div>
                    <div className="stat-label">Valide</div>
                  </div>
                  <div className="stat-item unique">
                    <div className="stat-number">{analysisResults.unique_molecules}</div>
                    <div className="stat-label">Uniche</div>
                  </div>
                  {referenceFile && (
                    <div className="stat-item novel">
                      <div className="stat-number">{analysisResults.novel_molecules}</div>
                      <div className="stat-label">Novel</div>
                    </div>
                  )}
                </div>
                <div className="novelty-rates">
                  <div className="rate-item">
                    <span className="rate-label">Tasso di Unicità:</span>
                    <span className="rate-value">{(analysisResults.uniqueness_rate * 100).toFixed(1)}%</span>
                  </div>
                  {referenceFile && (
                    <div className="rate-item">
                      <span className="rate-label">Tasso di Novelty:</span>
                      <span className="rate-value">{(analysisResults.novelty_rate * 100).toFixed(1)}%</span>
                    </div>
                  )}
                </div>
                <div className="analysis-time">
                  Analisi completata in {analysisResults.processing_time}s
                </div>
              </div>
            </div>
          )}

          {/* Sezione Filtro di Coordinazione */}
          {selectedFile && originalMolecules.length > 0 && analysisComplete && (
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
          {molecules.length > 0 && analysisComplete && (
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
        {isAnalyzing ? (
          <div className="analysis-loading">
            <div className="loading-spinner"></div>
            <div className="loading-text">
              <h3>Calcolo statistiche in corso...</h3>
              <p>Analisi di <b>unicità</b> e <b>novelty</b> delle molecole, potrebbero volerci <b>alcuni minuti</b></p>
              {referenceFile && <p>Confronto con file di riferimento: {referenceFile}</p>}
            </div>
          </div>
        ) : analysisError ? (
          <div className="error">{analysisError}</div>
        ) : analysisComplete && molecules.length > 0 ? (
          <MoleculeGrid
            molecules={molecules}
            analysisResults={analysisResults}
            displayStats={getDisplayStats()}
            validationResults={validationResults}
            isFiltered={isFiltered}
            referenceFile={referenceFile}
          />
        ) : selectedFile && analysisComplete ? (
          <div className="error">
            {isFiltered ?
              "Nessuna molecola corrisponde ai criteri di filtro selezionati" :
              "Nessuna molecola SMILES valida riconosciuta nel file selezionato"
            }
          </div>
        ) : (
          <div className="instructions">
            <h2>Benvenuto nel Visualizzatore Molecole</h2>
            <p>Seleziona un file CSV contenente SMILES dalla sidebar per iniziare l'analisi.</p>
            <div className="instructions-details">
              <h3>Funzionalità disponibili:</h3>
              <ul>
                <li><strong>Analisi di Unicità:</strong> Identifica molecole duplicate usando canonicalizzazione SMILES</li>
                <li><strong>Analisi di Novelty:</strong> Confronta con un file di riferimento per trovare nuove molecole</li>
                <li><strong>Filtri di Coordinazione:</strong> Filtra per metalli e numero di coordinazione</li>
                <li><strong>Validazione 3D:</strong> Verifica quali molecole possono generare strutture 3D</li>
                <li><strong>Visualizzazione Interattiva:</strong> Esplora le strutture molecolari in 2D e 3D</li>
              </ul>
            </div>
          </div>
        )}
      </main>
    </div>
  );
}

export default App;