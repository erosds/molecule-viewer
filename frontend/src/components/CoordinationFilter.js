// Crea il file frontend/src/components/CoordinationFilter.js

import React, { useState, useEffect } from 'react';
import './CoordinationFilter.css';

const CoordinationFilter = ({ selectedFile, onFilterApplied, onStatsUpdate }) => {
  const [coordinationStats, setCoordinationStats] = useState(null);
  const [minCoordination, setMinCoordination] = useState(0);
  const [maxCoordination, setMaxCoordination] = useState(12);
  const [isLoading, setIsLoading] = useState(false);
  const [isFiltering, setIsFiltering] = useState(false);
  const [error, setError] = useState(null);
  
  // Carica le statistiche quando cambia il file selezionato
  useEffect(() => {
    if (selectedFile) {
      loadCoordinationStats();
    } else {
      setCoordinationStats(null);
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
      const response = await fetch('/api/filter-by-coordination', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify({
          csv_file: selectedFile,
          min_coordination: minCoordination,
          max_coordination: maxCoordination
        }),
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
  };

  if (!selectedFile) {
    return null;
  }

  return (
    <div className="coordination-filter">
      <div className="coordination-filter-header">
        <h3>🧬 Filtro per Numero di Coordinazione</h3>
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
              <div className="stat-item range">
                <div className="stat-number">
                  {coordinationStats.coordination_range.min}-{coordinationStats.coordination_range.max}
                </div>
                <div className="stat-label">Range Coordinazione</div>
              </div>
            </div>

            {coordinationStats.metal_elements_found.length > 0 && (
              <div className="metal-elements">
                <span className="elements-label">Metalli trovati:</span>
                <div className="elements-list">
                  {coordinationStats.metal_elements_found.map(element => (
                    <span key={element} className="element-tag">{element}</span>
                  ))}
                </div>
              </div>
            )}
          </div>

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
                  onChange={(e) => setMinCoordination(parseInt(e.target.value) || 0)}
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
                  onChange={(e) => setMaxCoordination(parseInt(e.target.value) || 12)}
                />
              </div>
            </div>

            <div className="filter-buttons">
              <button 
                className="filter-button apply"
                onClick={handleFilter}
                disabled={isFiltering || minCoordination > maxCoordination}
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
                {Object.entries(coordinationStats.coordination_distribution)
                  .sort(([a], [b]) => parseInt(a) - parseInt(b))
                  .map(([coord, count]) => (
                    <div key={coord} className="distribution-bar">
                      <div className="bar-label">{coord}</div>
                      <div 
                        className="bar-fill"
                        style={{
                          width: `${(count / Math.max(...Object.values(coordinationStats.coordination_distribution))) * 100}%`,
                          backgroundColor: (parseInt(coord) >= minCoordination && parseInt(coord) <= maxCoordination) 
                            ? '#007aff' : '#d2d2d7'
                        }}
                      ></div>
                      <div className="bar-count">{count}</div>
                    </div>
                  ))}
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