import React, { useState } from 'react';
import MoleculeViewer from './MoleculeViewer';
import './MoleculeGrid.css';

const MoleculeGrid = ({ molecules, newMoleculesCount = 0 }) => {
  const [currentPage, setCurrentPage] = useState(1);
  const [selectedMolecule, setSelectedMolecule] = useState(null);
  const [viewerOpen, setViewerOpen] = useState(false);
  const [loading3D, setLoading3D] = useState(false);
  const [modelError, setModelError] = useState(null);

  // Configurazione per il numero di molecole visualizzate
  const itemsPerRow = 5; // Modificato a 5 molecole per riga
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

  return (
    <div className="molecule-grid-container">
      <div className="molecule-grid-info">
        <p>
          <span className='new-molecules'>Totale molecole uniche generate: {molecules.length}</span> |
          <span className='unique-molecules'> Totale molecole nuove: {newMoleculesCount}</span> |
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
              <img
                src={getMoleculeImageUrl(molecule.smiles)}
                alt={`Struttura di ${molecule.smiles}`}
                onError={(e) => {
                  e.target.onerror = null;
                  e.target.src = '/placeholder-molecule.png';
                }}
              />
            </div>
            <div className="molecule-info">
              <p title={molecule.smiles}>
                {molecule.smiles.length > 20
                  ? `${molecule.smiles.substring(0, 17)}...`
                  : molecule.smiles}
              </p>
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