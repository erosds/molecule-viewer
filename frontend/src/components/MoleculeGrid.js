import React, { useState } from 'react';
import MoleculeViewer from './MoleculeViewer';
import './MoleculeGrid.css';

const MoleculeGrid = ({ molecules, analysisResults, displayStats, validationResults, isFiltered, referenceFile }) => {
  const [currentPage, setCurrentPage] = useState(1);
  const [selectedMolecule, setSelectedMolecule] = useState(null);
  const [viewerOpen, setViewerOpen] = useState(false);
  const [loading3D, setLoading3D] = useState(false);
  const [modelError, setModelError] = useState(null);
  const [activeBadgeFilters, setActiveBadgeFilters] = useState([]);

  // Configurazione per il numero di molecole visualizzate
  const itemsPerRow = 6;
  const rowsPerPage = 6;
  const itemsPerPage = itemsPerRow * rowsPerPage;

  // Funzione per gestire i filtri badge con logica mutuamente esclusiva
  const handleBadgeFilterToggle = (badgeType) => {
    setActiveBadgeFilters(prev => {
      let newFilters = [...prev];
      
      if (prev.includes(badgeType)) {
        // Se il filtro è già attivo, rimuovilo
        newFilters = newFilters.filter(type => type !== badgeType);
      } else {
        // Se il filtro non è attivo, aggiungilo rimuovendo il mutuamente esclusivo
        if (badgeType === 'unique') {
          newFilters = newFilters.filter(type => type !== 'duplicate');
        } else if (badgeType === 'duplicate') {
          newFilters = newFilters.filter(type => type !== 'unique');
        } else if (badgeType === 'novel') {
          newFilters = newFilters.filter(type => type !== 'not-novel');
        } else if (badgeType === 'not-novel') {
          newFilters = newFilters.filter(type => type !== 'novel');
        }
        newFilters.push(badgeType);
      }
      
      return newFilters;
    });
  };

  // Funzione per filtrare le molecole con logica AND tra categorie diverse
  const getFilteredMolecules = () => {
    if (activeBadgeFilters.length === 0) {
      return molecules;
    }
    
    return molecules.filter(molecule => {
      // Dividi i filtri per categoria
      const uniquenessFilter = activeBadgeFilters.find(f => f === 'unique' || f === 'duplicate');
      const noveltyFilter = activeBadgeFilters.find(f => f === 'novel' || f === 'not-novel');
      const validationFilter = activeBadgeFilters.includes('validated');
      
      // Ogni categoria deve essere soddisfatta (AND)
      let passesUniqueness = true;
      let passesNovelty = true;
      let passesValidation = true;
      
      if (uniquenessFilter) {
        if (uniquenessFilter === 'unique') {
          passesUniqueness = molecule.is_unique;
        } else if (uniquenessFilter === 'duplicate') {
          passesUniqueness = !molecule.is_unique;
        }
      }
      
      if (noveltyFilter && referenceFile) {
        if (noveltyFilter === 'novel') {
          passesNovelty = molecule.is_novel;
        } else if (noveltyFilter === 'not-novel') {
          passesNovelty = !molecule.is_novel;
        }
      }
      
      if (validationFilter) {
        passesValidation = isMoleculeValidated(molecule.smiles);
      }
      
      return passesUniqueness && passesNovelty && passesValidation;
    });
  };

  // Ottieni le molecole filtrate
  const displayedMolecules = getFilteredMolecules();

  // Calcolo del numero totale di pagine basato sulle molecole filtrate
  const totalPages = Math.ceil(displayedMolecules.length / itemsPerPage);

  // Molecole per la pagina corrente
  const currentMolecules = displayedMolecules.slice(
    (currentPage - 1) * itemsPerPage,
    currentPage * itemsPerPage
  );

  // Reset alla prima pagina quando cambiano i filtri
  React.useEffect(() => {
    setCurrentPage(1);
  }, [activeBadgeFilters]);

  // Funzione per verificare se una molecola è validata per il 3D
  const isMoleculeValidated = (smiles) => {
    if (!validationResults) return false;
    return validationResults.valid_smiles && validationResults.valid_smiles.includes(smiles);
  };

  const handleMoleculeClick = async (molecule) => {
    setSelectedMolecule(molecule);
    setViewerOpen(true);
    setLoading3D(true);
    setModelError(null);

    try {
      // Richiesta al backend per generare il modello 3D
      const response = await fetch(`/api/generate-3d`, {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify({ smiles: molecule.smiles }),
      });

      if (!response.ok) {
        const errorData = await response.json();
        throw new Error(errorData.detail || 'Errore nella generazione del modello 3D');
      }

      // La risposta contiene il path del file generato
      const data = await response.json();
      molecule.modelPath = data.model_path;

    } catch (err) {
      setModelError(`Errore: ${err.message}`);
    } finally {
      setLoading3D(false);
    }
  };

  const closeViewer = () => {
    setViewerOpen(false);
    setSelectedMolecule(null);
  };

  const nextPage = () => {
    if (currentPage < totalPages) {
      setCurrentPage(currentPage + 1);
    }
  };

  const prevPage = () => {
    if (currentPage > 1) {
      setCurrentPage(currentPage - 1);
    }
  };

  // Funzione per ottenere l'URL dell'immagine 2D della molecola dal nostro backend
  const getMoleculeImageUrl = (smiles) => {
    return `/api/molecule-2d/${encodeURIComponent(smiles)}`;
  };

  // Funzione per ottenere le classi CSS per ogni molecola
  const getMoleculeCardClasses = (molecule) => {
    let classes = ['molecule-card'];
    
    if (isMoleculeValidated(molecule.smiles)) {
      classes.push('validated-3d');
    }
    
    if (!molecule.is_unique) {
      classes.push('duplicate');
    }
    
    if (referenceFile && !molecule.is_novel) {
      classes.push('not-novel');
    }
    
    return classes.join(' ');
  };

  return (
    <div className="molecule-grid-container">
      <div className="molecule-grid-header">
        <div className="grid-title">
          <h2>Risultati dell'Analisi</h2>
          {(isFiltered || activeBadgeFilters.length > 0) && (
            <span className="filtered-indicator">🔍 Filtrate</span>
          )}
        </div>
        
        <div className="molecule-grid-info">
          <div className="info-stats">
            <span className="stat-display total">
              <strong>{displayedMolecules.length}</strong> molecole visualizzate
            </span>
            
            <span className="stat-display unique">
              <strong>{displayStats?.unique || 0}</strong> uniche
            </span>
            
            {referenceFile && (
              <span className="stat-display novel">
                <strong>{displayStats?.novel || 0}</strong> novel
              </span>
            )}
            
            {validationResults && (
              <span className="stat-display validated">
                <strong>{validationResults.valid_count || 0}</strong> validabili 3D
              </span>
            )}
          </div>
          
          <span className="pagination-info">
            Pagina {currentPage} di {totalPages}
          </span>
        </div>

        {/* Sezione Filtri Badge */}
        <div className="molecule-legend">
          <div className="legend-items">
            <div 
              className={`legend-item ${activeBadgeFilters.includes('unique') ? 'active' : ''}`}
              onClick={() => handleBadgeFilterToggle('unique')}
            >
              <div className="legend-badge unique">U</div>
              <span>Unica</span>
            </div>

            <div 
              className={`legend-item ${activeBadgeFilters.includes('duplicate') ? 'active' : ''}`}
              onClick={() => handleBadgeFilterToggle('duplicate')}
            >
              <div className="legend-badge duplicate">D</div>
              <span>Duplicata</span>
            </div>

            {referenceFile && (
              <>
                <div 
                  className={`legend-item ${activeBadgeFilters.includes('novel') ? 'active' : ''}`}
                  onClick={() => handleBadgeFilterToggle('novel')}
                >
                  <div className="legend-badge novel">N</div>
                  <span>New</span>
                </div>

                <div 
                  className={`legend-item ${activeBadgeFilters.includes('not-novel') ? 'active' : ''}`}
                  onClick={() => handleBadgeFilterToggle('not-novel')}
                >
                  <div className="legend-badge not-novel">O</div>
                  <span>Old</span>
                </div>
              </>
            )}

            {validationResults && (
              <div 
                className={`legend-item ${activeBadgeFilters.includes('validated') ? 'active' : ''}`}
                onClick={() => handleBadgeFilterToggle('validated')}
              >
                <div className="legend-badge validated">✓</div>
                <span>Validata 3D</span>
              </div>
            )}
          </div>
          
          {activeBadgeFilters.length > 0 && (
            <button 
              onClick={() => setActiveBadgeFilters([])}
              style={{ 
                marginTop: '1rem', 
                padding: '0.5rem 1rem', 
                background: '#dc3545', 
                color: 'white', 
                border: 'none', 
                borderRadius: '6px', 
                cursor: 'pointer' 
              }}
            >
              Rimuovi tutti i filtri
            </button>
          )}
        </div>
      </div>

      {/* Griglia delle molecole */}
      <div className="molecule-grid">
        {currentMolecules.map((molecule, index) => (
          <div
            key={`${molecule.smiles}-${index}`}
            className={getMoleculeCardClasses(molecule)}
            onClick={() => handleMoleculeClick(molecule)}
          >
            <div className="molecule-image">
              <img
                src={getMoleculeImageUrl(molecule.smiles)}
                alt={`Struttura di ${molecule.smiles}`}
                onError={(e) => {
                  e.target.onerror = null;
                  e.target.src = '/placeholder-molecule.png';
                }}
              />
              
              {/* Badge per le caratteristiche */}
              <div className="molecule-badges">
                {molecule.is_unique && (
                  <div className="molecule-badge unique" title="Molecola unica">
                    U
                  </div>
                )}
                
                {referenceFile && molecule.is_novel && (
                  <div className="molecule-badge novel" title="Molecola novel (non presente nel riferimento)">
                    N
                  </div>
                )}
                
                {isMoleculeValidated(molecule.smiles) && (
                  <div className="molecule-badge validated" title="Molecola validabile per 3D">
                    ✓
                  </div>
                )}
                
                {!molecule.is_unique && (
                  <div className="molecule-badge duplicate" title="Molecola duplicata">
                    D
                  </div>
                )}
                
                {referenceFile && !molecule.is_novel && (
                  <div className="molecule-badge not-novel" title="Molecola conosciuta (presente nel riferimento)">
                    O
                  </div>
                )}
              </div>
            </div>
            
            <div className="molecule-info">
              <p className="molecule-smiles" title={molecule.smiles}>
                {molecule.smiles.length > 25
                  ? `${molecule.smiles.substring(0, 22)}...`
                  : molecule.smiles}
              </p>
              
              <div className="molecule-status">
                {molecule.is_unique ? (
                  <span className="status-unique">Unica</span>
                ) : (
                  <span className="status-duplicate">Duplicata</span>
                )}
                
                {referenceFile && (
                  <span className={molecule.is_novel ? "status-novel" : "status-known"}>
                    {molecule.is_novel ? "Novel" : "Conosciuta"}
                  </span>
                )}
                
                {isMoleculeValidated(molecule.smiles) && (
                  <span className="status-validated">3D OK</span>
                )}
              </div>
            </div>
          </div>
        ))}
      </div>

      <div className="pagination-controls">
        <button
          onClick={prevPage}
          disabled={currentPage === 1}
          className="pagination-button"
        >
          ← Precedente
        </button>
        <span className="page-indicator">
          {currentPage} / {totalPages}
        </span>
        <button
          onClick={nextPage}
          disabled={currentPage === totalPages}
          className="pagination-button"
        >
          Successivo →
        </button>
      </div>

      {viewerOpen && selectedMolecule && (
        <MoleculeViewer
          molecule={selectedMolecule}
          onClose={closeViewer}
          loading={loading3D}
          error={modelError}
        />
      )}
    </div>
  );
};

export default MoleculeGrid;