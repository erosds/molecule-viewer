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

    } catch (err) {
      setBatchError(`Errore: ${err.message}`);
    } finally {
      setIsGeneratingBatch(false);
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
  📂            File salvati in: <code>backend/public/molecules/{batchResult.batch_folder}/</code>
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