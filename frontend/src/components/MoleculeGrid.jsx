import React, { useState } from 'react';
import MoleculeViewer from './MoleculeViewer';
import PlaceholderMolecule from './PlaceholderMolecule';
import axios from 'axios';
import './MoleculeGrid.css';

const MoleculeGrid = ({ molecules }) => {
  const [currentPage, setCurrentPage] = useState(1);
  const [selectedMolecule, setSelectedMolecule] = useState(null);
  const [viewerOpen, setViewerOpen] = useState(false);
  const [loading3D, setLoading3D] = useState(false);
  const [modelError, setModelError] = useState(null);
  
  // Configurazione per il numero di molecole visualizzate
  const itemsPerRow = 5;
  const rowsPerPage = 6;
  const itemsPerPage = itemsPerRow * rowsPerPage;
  
  // Calcolo del numero totale di pagine
  const totalPages = Math.ceil(molecules.length / itemsPerPage);
  
  // Molecole per la pagina corrente
  const currentMolecules = molecules.slice(
    (currentPage - 1) * itemsPerPage,
    currentPage * itemsPerPage
  );

  const handleMoleculeClick = async (molecule) => {
    setSelectedMolecule(molecule);
    setViewerOpen(true);
    setLoading3D(true);
    setModelError(null);
    
    try {
      // Richiesta al backend per generare il modello 3D
      const response = await axios.post('/api/generate-3d', { 
        smiles: molecule.smiles 
      });
      
      // La risposta contiene il path del file generato e proprietà aggiuntive
      const { model_path, svg_path, properties } = response.data;
      
      // Aggiorna l'oggetto molecola con i nuovi dati
      setSelectedMolecule({
        ...molecule,
        modelPath: model_path,
        svgPath: svg_path,
        properties: properties
      });
      
    } catch (err) {
      console.error("Errore nella generazione del modello 3D:", err);
      let errorMessage = "Errore nella generazione del modello 3D";
      
      if (err.response && err.response.data) {
        errorMessage = err.response.data.detail || errorMessage;
      }
      
      setModelError(errorMessage);
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
      window.scrollTo({ top: 0, behavior: 'smooth' });
    }
  };

  const prevPage = () => {
    if (currentPage > 1) {
      setCurrentPage(currentPage - 1);
      window.scrollTo({ top: 0, behavior: 'smooth' });
    }
  };

  // Funzione per gestire il caricamento fallito dell'immagine
  const handleImageError = (e, molecule) => {
    e.target.style.display = 'none';
    const container = e.target.parentNode;
    
    // Verifica se c'è già un elemento SVG (PlaceholderMolecule)
    if (!container.querySelector('svg')) {
      const placeholder = document.createElement('div');
      placeholder.className = 'molecule-placeholder';
      container.appendChild(placeholder);
      
      // Utilizziamo React per renderizzare il componente PlaceholderMolecule
      import('react-dom/client').then(({ createRoot }) => {
        const root = createRoot(placeholder);
        root.render(<PlaceholderMolecule />);
      });
    }
  };

  return (
    <div className="molecule-grid-container">
      <div className="molecule-grid-info">
        <p>
          Visualizzazione di {Math.min(itemsPerPage, molecules.length - (currentPage - 1) * itemsPerPage)} molecole 
          (pagina {currentPage} di {totalPages})
        </p>
      </div>
      
      <div className="molecule-grid">
        {currentMolecules.map((molecule) => (
          <div 
            key={molecule.id} 
            className="molecule-card"
            onClick={() => handleMoleculeClick(molecule)}
          >
            <div className="molecule-image">
              {molecule.svgPath ? (
                <img 
                  src={molecule.svgPath} 
                  alt={`Struttura di ${molecule.smiles}`} 
                  onError={(e) => handleImageError(e, molecule)}
                />
              ) : (
                <img 
                  src={`https://cactus.nci.nih.gov/chemical/structure/${encodeURIComponent(molecule.smiles)}/image`} 
                  alt={`Struttura di ${molecule.smiles}`} 
                  onError={(e) => handleImageError(e, molecule)}
                />
              )}
            </div>
            <div className="molecule-info">
              {molecule.name ? (
                <h4 title={molecule.name}>
                  {molecule.name.length > 25 
                    ? `${molecule.name.substring(0, 22)}...` 
                    : molecule.name}
                </h4>
              ) : null}
              <p title={molecule.smiles}>
                {molecule.smiles.length > 20 
                  ? `${molecule.smiles.substring(0, 17)}...` 
                  : molecule.smiles}
              </p>
              {molecule.formula && <p className="molecule-formula">{molecule.formula}</p>}
            </div>
          </div>
        ))}
      </div>
      
      {totalPages > 1 && (
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
      )}

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