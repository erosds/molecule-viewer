// Crea il file frontend/src/components/ValidationButton.js

import React, { useState } from 'react';
import './ValidationButton.css';

const ValidationButton = ({ molecules, onValidationComplete }) => {
  const [isValidating, setIsValidating] = useState(false);
  const [validationResult, setValidationResult] = useState(null);
  const [error, setError] = useState(null);

  const handleValidation = async () => {
    if (!molecules || molecules.length === 0) {
      setError('Nessuna molecola da validare');
      return;
    }

    setIsValidating(true);
    setError(null);
    setValidationResult(null);

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
            
            <div className="stat-item total">
              <div className="stat-number">{validationResult.total_molecules}</div>
              <div className="stat-label">Totale</div>
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