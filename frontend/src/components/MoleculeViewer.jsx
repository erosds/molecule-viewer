import React, { useEffect, useRef } from 'react';
import './MoleculeViewer.css';

/* global $3Dmol */ // Diciamo a ESLint che $3Dmol è una variabile globale

const MoleculeViewer = ({ molecule, onClose, loading, error }) => {
  const viewerRef = useRef(null);
  const viewer3DRef = useRef(null);

  useEffect(() => {
    // Pulizia del visualizzatore 3D quando il componente viene smontato
    return () => {
      if (viewer3DRef.current) {
        // Rimuovere tutti gli event listener e pulire
        try {
          viewer3DRef.current.clear();
        } catch (err) {
          console.error("Errore durante la pulizia del visualizzatore:", err);
        }
      }
    };
  }, []);

  useEffect(() => {
    if (!loading && !error && molecule.modelPath && viewerRef.current) {
      try {
        // Assicuriamo che il container sia visibile e renderizzabile
        if (!viewerRef.current.offsetWidth) {
          console.warn("Il container del visualizzatore non ha una larghezza visibile");
          return;
        }

        // Verifichiamo che $3Dmol sia disponibile
        if (typeof $3Dmol === 'undefined') {
          console.error("La libreria 3Dmol.js non è stata caricata correttamente");
          if (viewerRef.current) {
            viewerRef.current.innerHTML = `
              <div class="error-message">
                <p>Errore: La libreria 3Dmol.js non è disponibile.</p>
                <p>Verifica che lo script sia stato caricato correttamente.</p>
              </div>
            `;
          }
          return;
        }

        // Inizializziamo il visualizzatore 3D
        const config = {
          backgroundColor: 'white',
          antialias: true,
          id: 'molecule-viewer-' + molecule.id
        };

        // Creiamo il visualizzatore
        viewer3DRef.current = $3Dmol.createViewer(viewerRef.current, config);
        
        // Carichiamo il modello dal percorso
        fetch(molecule.modelPath)
          .then(response => {
            if (!response.ok) {
              throw new Error(`Errore nella richiesta HTTP: ${response.status}`);
            }
            return response.text();
          })
          .then(pdbData => {
            // Aggiungiamo il modello al visualizzatore
            viewer3DRef.current.addModel(pdbData, "pdb");
            
            // Impostiamo lo stile di visualizzazione
            viewer3DRef.current.setStyle({}, {
              stick: {}, // Visualizzazione a bastoncini
              sphere: { scale: 0.3 } // Piccole sfere per gli atomi
            });
            
            // Aggiungiamo etichette agli atomi
            const atoms = viewer3DRef.current.getModel().selectedAtoms({});
            for (let i = 0; i < atoms.length; i++) {
              const atom = atoms[i];
              if (atom.atom !== "H") { // Evitiamo di etichettare gli idrogeni
                viewer3DRef.current.addLabel(atom.atom, {
                  position: { x: atom.x, y: atom.y, z: atom.z },
                  backgroundColor: "#FFFFFF",
                  backgroundOpacity: 0.5,
                  fontColor: "black",
                  fontSize: 12
                });
              }
            }
            
            // Impostiamo la vista
            viewer3DRef.current.zoomTo();
            
            // Rendiamo il modello
            viewer3DRef.current.render();
            
            // Aggiungiamo controlli per la rotazione
            viewer3DRef.current.rotate(30, { x: 1, y: 0, z: 0 });
            viewer3DRef.current.render();
            
            // Animazione iniziale
            let rotate = 0;
            const animationInterval = setInterval(() => {
              if (rotate < 5) {
                viewer3DRef.current.rotate(15, { x: 0, y: 1, z: 0 });
                viewer3DRef.current.render();
                rotate++;
              } else {
                clearInterval(animationInterval);
              }
            }, 100);
          })
          .catch(err => {
            console.error("Errore nel caricamento del modello PDB:", err);
            if (viewerRef.current) {
              viewerRef.current.innerHTML = `
                <div class="error-message">
                  <p>Errore nel caricamento del modello:</p>
                  <p>${err.message}</p>
                </div>
              `;
            }
          });
      } catch (err) {
        console.error("Errore nell'inizializzazione del visualizzatore 3D:", err);
        if (viewerRef.current) {
          viewerRef.current.innerHTML = `
            <div class="error-message">
              <p>Errore nell'inizializzazione del visualizzatore:</p>
              <p>${err.message}</p>
            </div>
          `;
        }
      }
    }
  }, [loading, error, molecule]);

  // Renderizziamo un semplice visualizzatore se la libreria 3D non è disponibile
  const renderFallbackViewer = () => {
    return (
      <div className="fallback-viewer">
        <p>SMILES: {molecule.smiles}</p>
        <div className="molecule-simple-view">
          <svg width="100%" height="100%" viewBox="0 0 200 200">
            <circle cx="100" cy="100" r="50" fill="#f0f0f0" stroke="#aaaaaa" strokeWidth="2" />
            <text x="100" y="105" textAnchor="middle" fontSize="10">Visualizzazione 3D non disponibile</text>
          </svg>
        </div>
      </div>
    );
  };

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
                {molecule.name && <p><strong>Nome:</strong> {molecule.name}</p>}
                {molecule.formula && <p><strong>Formula:</strong> {molecule.formula}</p>}
                {molecule.properties && molecule.properties.mol_weight && 
                  <p><strong>Peso molecolare:</strong> {molecule.properties.mol_weight.toFixed(2)} g/mol</p>
                }
              </div>
              <div className="viewer-container" ref={viewerRef} id={`molecule-viewer-container-${molecule.id}`}>
                {/* Il visualizzatore 3D verrà inserito qui */}
                {typeof $3Dmol === 'undefined' && renderFallbackViewer()}
              </div>
            </>
          )}
        </div>
      </div>
    </div>
  );
};

export default MoleculeViewer;