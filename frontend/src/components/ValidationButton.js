// Aggiornamento di frontend/src/components/ValidationButton.js

import React, { useState } from 'react';
import './ValidationButton.css';

const ValidationButton = ({ molecules, onValidationComplete }) => {
  const [isValidating, setIsValidating] = useState(false);
  const [validationResult, setValidationResult] = useState(null);
  const [error, setError] = useState(null);
  
  // Nuovi stati per la generazione batch
  const [isGeneratingBatch, setIsGeneratingBatch] = useState(false);
  const [batchResult, setBatchResult] = useState(null);
  const [batchError, setBatchError] = useState(null);
  
  // Nuovi stati per il download ZIP
  const [isDownloadingZip, setIsDownloadingZip] = useState(false);
  const [batchInfo, setBatchInfo] = useState(null);

  const handleValidation = async () => {
    if (!molecules || molecules.length === 0) {
      setError('Nessuna molecola da validare');
      return;
    }

    setIsValidating(true);
    setError(null);
    setValidationResult(null);
    // Reset anche i risultati del batch quando si fa una nuova validazione
    setBatchResult(null);
    setBatchError(null);
    setBatchInfo(null);

    try {
      // Estrai solo gli SMILES dalle molecole
      const smilesList = molecules.map(mol => mol.smiles);

      const response = await fetch('/api/validate-molecules', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify({
          smiles_list: smilesList
        }),
      });

      if (!response.ok) {
        const errorData = await response.json();
        throw new Error(errorData.detail || 'Errore nella validazione');
      }

      const result = await response.json();
      setValidationResult(result);
      
      // Comunica i risultati al componente genitore se necessario
      if (onValidationComplete) {
        onValidationComplete(result);
      }

    } catch (err) {
      setError(`Errore: ${err.message}`);
    } finally {
      setIsValidating(false);
    }
  };

  const getSuccessRate = () => {
    if (!validationResult) return 0;
    return Math.round((validationResult.valid_molecules / validationResult.total_molecules) * 100);
  };

  const handleBatchGeneration = async () => {
    if (!validationResult || !validationResult.valid_smiles || validationResult.valid_smiles.length === 0) {
      setBatchError('Nessuna struttura validata da generare');
      return;
    }

    setIsGeneratingBatch(true);
    setBatchError(null);
    setBatchResult(null);

    try {
      const response = await fetch('/api/generate-xyz-batch', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify({
          smiles_list: validationResult.valid_smiles
        }),
      });

      if (!response.ok) {
        const errorData = await response.json();
        throw new Error(errorData.detail || 'Errore nella generazione batch');
      }

      const result = await response.json();
      setBatchResult(result);

      // Carica le informazioni della cartella batch
      if (result.batch_folder) {
        loadBatchInfo(result.batch_folder);
      }

    } catch (err) {
      setBatchError(`Errore: ${err.message}`);
    } finally {
      setIsGeneratingBatch(false);
    }
  };

  const loadBatchInfo = async (batchFolder) => {
    try {
      const response = await fetch(`/api/batch-info/${batchFolder}`);
      if (response.ok) {
        const info = await response.json();
        setBatchInfo(info);
      }
    } catch (err) {
      console.error('Errore nel caricamento info batch:', err);
    }
  };

  const handleDownloadZip = async () => {
    if (!batchResult || !batchResult.batch_folder || isDownloadingZip) return;

    setIsDownloadingZip(true);
    
    try {
      const response = await fetch(`/api/download-xyz-batch/${batchResult.batch_folder}`);
      
      if (!response.ok) {
        throw new Error(`Errore nel download: ${response.status}`);
      }
      
      // Crea un blob dal contenuto della risposta
      const blob = await response.blob();
      
      // Crea un URL temporaneo per il blob
      const url = window.URL.createObjectURL(blob);
      
      // Crea un elemento <a> per scaricare il file
      const link = document.createElement('a');
      link.href = url;
      link.download = `${batchResult.batch_folder}_structures.zip`;
      
      // Simula il click per avviare il download
      document.body.appendChild(link);
      link.click();
      
      // Pulizia
      document.body.removeChild(link);
      window.URL.revokeObjectURL(url);
      
    } catch (err) {
      console.error("Errore durante il download ZIP:", err);
      alert(`Errore durante il download dell'archivio: ${err.message}`);
    } finally {
      setIsDownloadingZip(false);
    }
  };

  return (
    <div className="validation-container">
      <div className="validation-header">
        <button 
          className="validation-button"
          onClick={handleValidation}
          disabled={isValidating || !molecules || molecules.length === 0}
        >
          {isValidating ? (
            <>
              <div className="validation-spinner"></div>
              Validazione in corso...
            </>
          ) : (
            <>
              🧪 Verifica Strutture 3D
            </>
          )}
        </button>
        
        {validationResult && (
          <div className="validation-summary">
            <span className="success-rate">
              {getSuccessRate()}% validabili ({validationResult.valid_molecules}/{validationResult.total_molecules})
            </span>
            <span className="processing-time">
              in {validationResult.processing_time}s
            </span>
          </div>
        )}
      </div>

      {/* Pulsante per la generazione batch - appare solo dopo la validazione */}
      {validationResult && validationResult.valid_molecules > 0 && (
        <div className="batch-generation-section">
          <button 
            className="batch-generation-button"
            onClick={handleBatchGeneration}
            disabled={isGeneratingBatch}
          >
            {isGeneratingBatch ? (
              <>
                <div className="validation-spinner"></div>
                Generazione file XYZ in corso... ({validationResult.valid_molecules} strutture)
              </>
            ) : (
              <>
                📁 Genera File XYZ ({validationResult.valid_molecules} strutture validate)
              </>
            )}
          </button>
          
          {batchResult && (
            <div className="batch-result">
              <div className="batch-stats">
                <span className="batch-success">
                  ✅ {batchResult.successfully_generated} file generati con successo
                </span>
                {batchResult.failed_generation > 0 && (
                  <span className="batch-failed">
                    ❌ {batchResult.failed_generation} generazioni fallite
                  </span>
                )}
                <span className="batch-time">
                  Completato in {batchResult.processing_time}s
                </span>
              </div>
              
              <div className="batch-location">
                📂 File salvati in: <code>backend/public/molecules/{batchResult.batch_folder}/</code>
                
                {/* Informazioni aggiuntive se disponibili */}
                {batchInfo && (
                  <div className="batch-info-details">
                    <small>
                      {batchInfo.file_count} file • {batchInfo.total_size_mb} MB
                    </small>
                  </div>
                )}
              </div>

              {/* Pulsante Download ZIP */}
              <div className="batch-download-section">
                <button 
                  className="download-zip-button"
                  onClick={handleDownloadZip}
                  disabled={isDownloadingZip}
                  title="Scarica tutti i file XYZ come archivio ZIP"
                >
                  {isDownloadingZip ? (
                    <>
                      <div className="validation-spinner"></div>
                      Download in corso...
                    </>
                  ) : (
                    <>
                      🗜️ Scarica Archivio ZIP ({batchResult.successfully_generated} file)
                    </>
                  )}
                </button>
              </div>
            </div>
          )}
          
          {batchError && (
            <div className="batch-error">
              {batchError}
            </div>
          )}
        </div>
      )}

      {error && (
        <div className="validation-error">
          {error}
        </div>
      )}

      {validationResult && (
        <div className="validation-results">
          <div className="validation-stats">
            <div className="stat-item valid">
              <div className="stat-number">{validationResult.valid_molecules}</div>
              <div className="stat-label">Strutture 3D Generabili</div>
            </div>
            
            <div className="stat-item invalid">
              <div className="stat-number">{validationResult.invalid_molecules}</div>
              <div className="stat-label">Non Generabili</div>
            </div>
          </div>

          {validationResult.invalid_molecules > 0 && (
            <div className="validation-details">
              <details>
                <summary>
                  Visualizza molecole non validabili ({validationResult.invalid_molecules})
                </summary>
                <div className="invalid-smiles-list">
                  {validationResult.invalid_smiles.map((smiles, index) => (
                    <div key={index} className="invalid-smiles-item">
                      {smiles.length > 50 ? `${smiles.substring(0, 47)}...` : smiles}
                    </div>
                  ))}
                </div>
              </details>
            </div>
          )}
        </div>
      )}
    </div>
  );
};

export default ValidationButton;