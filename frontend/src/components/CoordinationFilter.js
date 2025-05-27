// Aggiorna il file frontend/src/components/CoordinationFilter.js

import React, { useState, useEffect } from 'react';
import './CoordinationFilter.css';

const CoordinationFilter = ({ selectedFile, onFilterApplied, onStatsUpdate }) => {
  const [coordinationStats, setCoordinationStats] = useState(null);
  const [minCoordination, setMinCoordination] = useState(0);
  const [maxCoordination, setMaxCoordination] = useState(12);
  const [selectedMetals, setSelectedMetals] = useState([]); // Nuovo stato per metalli selezionati
  const [isLoading, setIsLoading] = useState(false);
  const [isFiltering, setIsFiltering] = useState(false);
  const [error, setError] = useState(null);

  // Carica le statistiche quando cambia il file selezionato
  useEffect(() => {
    if (selectedFile) {
      loadCoordinationStats();
    } else {
      setCoordinationStats(null);
      setSelectedMetals([]);  // Reset metalli selezionati
    }
  }, [selectedFile]);

  const loadCoordinationStats = async () => {
    setIsLoading(true);
    setError(null);

    try {
      const response = await fetch(`/api/coordination-stats/${selectedFile}`);
      if (!response.ok) {
        throw new Error('Errore nel caricamento delle statistiche');
      }

      const stats = await response.json();
      setCoordinationStats(stats);

      // Imposta i valori di default basati sui dati
      if (stats.coordination_range) {
        setMinCoordination(stats.coordination_range.min);
        setMaxCoordination(stats.coordination_range.max);
      }

      // Reset metalli selezionati quando cambiano le statistiche
      setSelectedMetals([]);

      // Comunica le statistiche al componente padre
      if (onStatsUpdate) {
        onStatsUpdate(stats);
      }

    } catch (err) {
      setError(`Errore: ${err.message}`);
    } finally {
      setIsLoading(false);
    }
  };

  const handleFilter = async () => {
    if (!selectedFile) return;

    setIsFiltering(true);
    setError(null);

    try {
      const requestBody = {
        csv_file: selectedFile,
        min_coordination: minCoordination === '' ? 0 : minCoordination,
        max_coordination: maxCoordination === '' ? 12 : maxCoordination
      };

      // Aggiungi i metalli selezionati solo se ce ne sono
      if (selectedMetals.length > 0) {
        requestBody.selected_metals = selectedMetals;
      }

      const response = await fetch('/api/filter-by-coordination', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify(requestBody),
      });

      if (!response.ok) {
        const errorData = await response.json();
        throw new Error(errorData.detail || 'Errore nel filtro');
      }

      const result = await response.json();

      // Comunica i risultati al componente padre
      if (onFilterApplied) {
        onFilterApplied(result.molecules, result);
      }

    } catch (err) {
      setError(`Errore: ${err.message}`);
    } finally {
      setIsFiltering(false);
    }
  };

  const resetFilter = () => {
    if (coordinationStats && coordinationStats.coordination_range) {
      setMinCoordination(coordinationStats.coordination_range.min);
      setMaxCoordination(coordinationStats.coordination_range.max);
    } else {
      setMinCoordination(0);
      setMaxCoordination(12);
    }
    setSelectedMetals([]);  // Reset anche i metalli selezionati
  };

  // Funzione per gestire la selezione/deselezione dei metalli
  const toggleMetal = (metal) => {
    setSelectedMetals(prev => {
      if (prev.includes(metal)) {
        return prev.filter(m => m !== metal);
      } else {
        return [...prev, metal];
      }
    });
  };

  // Funzione per selezionare tutti i metalli
  const selectAllMetals = () => {
    if (coordinationStats && coordinationStats.metal_elements_found) {
      setSelectedMetals([...coordinationStats.metal_elements_found]);
    }
  };

  // Funzione per deselezionare tutti i metalli
  const deselectAllMetals = () => {
    setSelectedMetals([]);
  };

  if (!selectedFile) {
    return null;
  }

  return (
    <div className="coordination-filter">
      <div className="coordination-filter-header">
        <h3>🧬 Filtro per Metalli e Numero di Coordinazione</h3>
        <button
          className="refresh-stats-button"
          onClick={loadCoordinationStats}
          disabled={isLoading}
        >
          {isLoading ? '⟳' : '↻'}
        </button>
      </div>

      {isLoading ? (
        <div className="coordination-loading">
          Analisi della coordinazione metallica in corso...
        </div>
      ) : error ? (
        <div className="coordination-error">
          {error}
        </div>
      ) : coordinationStats ? (
        <>
          <div className="coordination-stats">
            <div className="stats-grid">
              <div className="stat-item">
                <div className="stat-number">{coordinationStats.total_molecules}</div>
                <div className="stat-label">Molecole Totali</div>
              </div>
              <div className="stat-item metal">
                <div className="stat-number">{coordinationStats.molecules_with_metals}</div>
                <div className="stat-label">Con Metalli</div>
              </div>
            </div>
          </div>

          {/* Sezione per la selezione dei metalli */}
          {coordinationStats.metal_elements_found.length > 0 && (
            <div className="metal-selection-section">
              <div className="metal-selection-header">
                <h4>Seleziona Metalli da Includere:</h4>
                <div className="metal-selection-controls">
                  <button
                    className="metal-control-button select-all"
                    onClick={selectAllMetals}
                    disabled={isFiltering}
                  >
                    ✓ Tutti
                  </button>
                  <button
                    className="metal-control-button deselect-all"
                    onClick={deselectAllMetals}
                    disabled={isFiltering}
                  >
                    ✗ Nessuno
                  </button>
                </div>
              </div>

              <div className="metal-selection-grid">
                {coordinationStats.metal_elements_found.map(metal => (
                  <button
                    key={metal}
                    className={`metal-toggle ${selectedMetals.includes(metal) ? 'selected' : ''}`}
                    onClick={() => toggleMetal(metal)}
                    disabled={isFiltering}
                  >
                    <span className="metal-symbol">{metal}</span>
                    {selectedMetals.includes(metal) && <span className="check-mark">✓</span>}
                  </button>
                ))}
              </div>

              {selectedMetals.length > 0 && (
                <div className="selected-metals-info">
                  <span className="selected-label">
                    Metalli selezionati ({selectedMetals.length}):
                  </span>
                  <div className="selected-metals-list">
                    {selectedMetals.map(metal => (
                      <span key={metal} className="selected-metal-tag">
                        {metal}
                      </span>
                    ))}
                  </div>
                </div>
              )}
            </div>
          )}

          <div className="coordination-filter-controls">
            <div className="filter-inputs">
              <div className="input-group">
                <label htmlFor="min-coord">Coordinazione Min:</label>
                <input
                  id="min-coord"
                  type="number"
                  min="0"
                  max="20"
                  value={minCoordination}
                  onChange={(e) => {
                    const value = e.target.value;
                    if (value === '') {
                      setMinCoordination('');
                    } else {
                      setMinCoordination(parseInt(value) || 0);
                    }
                  }}
                  onBlur={(e) => {
                    if (e.target.value === '') {
                      setMinCoordination(0);
                    }
                  }}
                />
              </div>

              <div className="input-group">
                <label htmlFor="max-coord">Coordinazione Max:</label>
                <input
                  id="max-coord"
                  type="number"
                  min="0"
                  max="20"
                  value={maxCoordination}
                  onChange={(e) => {
                    const value = e.target.value;
                    if (value === '') {
                      setMaxCoordination('');
                    } else {
                      setMaxCoordination(parseInt(value) || 12);
                    }
                  }}
                  onBlur={(e) => {
                    if (e.target.value === '') {
                      setMaxCoordination(12);
                    }
                  }}
                />
              </div>
            </div>

            <div className="filter-buttons">
              <button
                className="filter-button apply"
                onClick={handleFilter}
                disabled={isFiltering || (minCoordination !== '' && maxCoordination !== '' && minCoordination > maxCoordination)}
              >
                {isFiltering ? (
                  <>
                    <div className="filter-spinner"></div>
                    Filtraggio...
                  </>
                ) : (
                  <>
                    🔍 Applica Filtro
                  </>
                )}
              </button>

              <button
                className="filter-button reset"
                onClick={resetFilter}
                disabled={isFiltering}
              >
                🔄 Reset
              </button>
            </div>
          </div>

          {Object.keys(coordinationStats.coordination_distribution).length > 0 && (
            <div className="coordination-distribution">
              <h4>Distribuzione Numero di Coordinazione:</h4>
              <div className="distribution-chart">
                {(() => {
                  const maxCount = Math.max(...Object.values(coordinationStats.coordination_distribution));
                  return Object.entries(coordinationStats.coordination_distribution)
                    .sort(([a], [b]) => parseInt(a) - parseInt(b))
                    .map(([coord, count]) => (
                      <div key={coord} className="distribution-bar">
                        <div className="bar-label">{coord}</div>
                        <div
                          className="bar-fill"
                          style={{
                            width: `${(count / maxCount) * 100}%`, // Era (count / Math.max(...Object.values(...))) * 100
                            backgroundColor: (parseInt(coord) >= minCoordination && parseInt(coord) <= maxCoordination)
                              ? '#007aff' : '#d2d2d7'
                          }}
                        ></div>
                        <div className="bar-count">{count}</div>
                      </div>
                    ));
                })()}
              </div>
            </div>
          )}
        </>
      ) : null}

      {error && (
        <div className="coordination-error">
          {error}
        </div>
      )}
    </div>
  );
};

export default CoordinationFilter;