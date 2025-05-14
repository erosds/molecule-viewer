import React, { useEffect, useRef } from 'react';
import './MoleculeViewer.css';

const MoleculeViewer = ({ molecule, onClose, loading, error }) => {
  const viewerRef = useRef(null);

  useEffect(() => {
    if (!loading && !error && molecule.modelPath) {
      // In un'implementazione reale, qui caricheremmo il modello 3D
      // usando librerie come 3Dmol.js o simili
      console.log(`Caricamento del modello: ${molecule.modelPath}`);
      
      // Simulazione del caricamento di un visualizzatore 3D
      // In un'app reale, useremmo codice come:
      /*
      const viewer = $3Dmol.createViewer(viewerRef.current, {
        backgroundColor: 'white',
      });
      
      fetch(molecule.modelPath)
        .then(response => response.text())
        .then(data => {
          viewer.addModel(data, 'pdb');
          viewer.setStyle({}, {stick: {}});
          viewer.zoomTo();
          viewer.render();
        });
      */
    }
  }, [loading, error, molecule]);

  return (
    <div className="molecule-viewer-overlay">
      <div className="molecule-viewer-container">
        <div className="molecule-viewer-header">
          <h3>Visualizzazione 3D</h3>
          <button className="close-button" onClick={onClose}>×</button>
        </div>
        
        <div className="molecule-viewer-content">
          {loading ? (
            <div className="loading-container">
              <div className="loading-spinner"></div>
              <p>Generazione del modello 3D in corso...</p>
            </div>
          ) : error ? (
            <div className="error-container">
              <p>{error}</p>
              <p>Impossibile generare il modello 3D per questa molecola.</p>
            </div>
          ) : (
            <>
              <div className="molecule-details">
                <p><strong>SMILES:</strong> {molecule.smiles}</p>
              </div>
              <div className="viewer-container" ref={viewerRef}>
                {/* Qui verrebbe renderizzato il modello 3D */}
                <div className="placeholder-message">
                  <p>Qui verrà visualizzato il modello 3D.</p>
                  <p className="small">File: {molecule.modelPath || 'Non disponibile'}</p>
                </div>
              </div>
            </>
          )}
        </div>
      </div>
    </div>
  );
};

export default MoleculeViewer;