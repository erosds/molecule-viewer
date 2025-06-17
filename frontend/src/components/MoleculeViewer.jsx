// Aggiornamento di frontend/src/components/MoleculeViewer.jsx

import React, { useEffect, useRef, useState } from 'react';
import './MoleculeViewer.css';

/* global $3Dmol */ // Diciamo a ESLint che $3Dmol è una variabile globale

const MoleculeViewer = ({ molecule, onClose, loading, error }) => {
  const viewerRef = useRef(null);
  const viewer3DRef = useRef(null);
  const [isDownloading, setIsDownloading] = useState(false);
  const [downloadError, setDownloadError] = useState(null);
  const downloadAbortRef = useRef(null);

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
      
      // Cancella eventuali download in corso
      if (downloadAbortRef.current) {
        downloadAbortRef.current.abort();
      }
    };
  }, []);

  // Funzione per il download del file XYZ con cancellazione
  const handleDownloadXYZ = async () => {
    if (!molecule.modelPath || isDownloading) return;

    // Cancella eventuali download precedenti
    if (downloadAbortRef.current) {
      downloadAbortRef.current.abort();
    }

    setIsDownloading(true);
    setDownloadError(null);
    
    // Crea un nuovo AbortController per il download
    const abortController = new AbortController();
    downloadAbortRef.current = abortController;
    
    try {
      // Estrai il nome del file dal percorso
      const filename = molecule.modelPath.split('/').pop();
      
      // Timeout per il download (15 secondi)
      const timeoutId = setTimeout(() => {
        abortController.abort();
      }, 15000);
      
      // Effettua la richiesta di download
      const response = await fetch(`/api/download-xyz/${filename}`, {
        signal: abortController.signal
      });
      
      clearTimeout(timeoutId);
      
      if (!response.ok) {
        throw new Error(`Errore nel download: ${response.status}`);
      }
      
      // Crea un blob dal contenuto della risposta
      const blob = await response.blob();
      
      // Verifica se il download è stato cancellato
      if (abortController.signal.aborted) {
        return;
      }
      
      // Crea un URL temporaneo per il blob
      const url = window.URL.createObjectURL(blob);
      
      // Crea un elemento <a> per scaricare il file
      const link = document.createElement('a');
      link.href = url;
      link.download = filename;
      
      // Simula il click per avviare il download
      document.body.appendChild(link);
      link.click();
      
      // Pulizia
      document.body.removeChild(link);
      window.URL.revokeObjectURL(url);
      
      console.log(`File ${filename} scaricato con successo`);
      
    } catch (err) {
      if (err.name === 'AbortError') {
        console.log("Download cancellato dall'utente");
        setDownloadError("Download cancellato");
      } else {
        console.error("Errore durante il download:", err);
        setDownloadError(`Errore durante il download: ${err.message}`);
      }
    } finally {
      setIsDownloading(false);
      downloadAbortRef.current = null;
      
      // Pulisci l'errore dopo qualche secondo
      setTimeout(() => {
        setDownloadError(null);
      }, 3000);
    }
  };

  // Funzione per cancellare il download in corso
  const cancelDownload = () => {
    if (downloadAbortRef.current) {
      downloadAbortRef.current.abort();
    }
  };

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
          id: 'molecule-viewer-' + molecule.id,
          width: viewerRef.current.clientWidth,
          height: viewerRef.current.clientHeight
        };

        // Creiamo il visualizzatore
        viewer3DRef.current = $3Dmol.createViewer(viewerRef.current, config);

        // Carichiamo il modello dal percorso con timeout
        const loadModelWithTimeout = () => {
          return Promise.race([
            fetch(molecule.modelPath),
            new Promise((_, reject) => 
              setTimeout(() => reject(new Error('Timeout nel caricamento del modello 3D')), 10000)
            )
          ]);
        };

        loadModelWithTimeout()
          .then(response => {
            if (!response.ok) {
              throw new Error(`Errore nella richiesta HTTP: ${response.status}`);
            }
            console.log(`File scaricato con successo da ${molecule.modelPath}`);
            return response.text();
          })
          .then(modelData => {
            // Log per debug
            console.log(`Contenuto del file (primi 100 caratteri): ${modelData.substring(0, 100)}...`);

            // Determiniamo se il file è XYZ o PDB in base al contenuto
            const isXYZ = !modelData.includes("ATOM") && !modelData.includes("HETATM");

            if (isXYZ) {
              // Parsing di un file XYZ
              console.log("Rilevato file in formato XYZ");

              try {
                // In 3Dmol.js, possiamo passare direttamente il contenuto del file XYZ
                // Non abbiamo bisogno di manipolarlo manualmente, basta specificare il formato corretto
                viewer3DRef.current.addModel(modelData, "xyz");
                console.log("Modello XYZ aggiunto con successo");
              } catch (err) {
                console.error("Errore nell'aggiunta del modello XYZ:", err);

                // Fallback: proviamo a fare il parsing manualmente
                const lines = modelData.split('\n');
                const numAtoms = parseInt(lines[0].trim(), 10);
                console.log(`Numero di atomi: ${numAtoms}`);

                // Formatiamo manualmente il file xyz in un formato che 3Dmol.js può leggere
                let formattedXYZ = `${numAtoms}\n${lines[1]}\n`;
                for (let i = 2; i < numAtoms + 2 && i < lines.length; i++) {
                  formattedXYZ += lines[i].trim() + '\n';
                }

                console.log("XYZ formattato manualmente:", formattedXYZ.substring(0, 200) + "...");
                viewer3DRef.current.addModel(formattedXYZ, "xyz");
              }
            } else {
              // Addizione di un modello PDB standard
              console.log("Rilevato file in formato PDB");
              viewer3DRef.current.addModel(modelData, "pdb");
            }

            // Impostiamo lo stile di visualizzazione
            viewer3DRef.current.setStyle({}, {
              stick: { radius: 0.15 }, // Bastoncini sottili come richiesto
              sphere: { scale: 0.25 }  // Piccole sfere per gli atomi come richiesto
            });

            // Aggiungiamo etichette solo agli atomi che non sono carbonio o idrogeno
            const atoms = viewer3DRef.current.getModel().selectedAtoms({});
            for (let i = 0; i < atoms.length; i++) {
              const atom = atoms[i];
              if (atom.atom !== "C" && atom.atom !== "H") { // Evitiamo di etichettare carbonio e idrogeno, come nel codice Python
                viewer3DRef.current.addLabel(atom.atom, {
                  position: { x: atom.x, y: atom.y, z: atom.z },
                  backgroundColor: "#FFFFFF",
                  backgroundOpacity: 0.7,
                  fontColor: "black",
                  fontSize: 12
                });
              }
            }

            // Impostiamo la vista per centrare perfettamente la molecola
            viewer3DRef.current.zoomTo();
            // Aggiungiamo un margine per migliorare la visualizzazione
            viewer3DRef.current.zoom(0.8);

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
            console.error("Errore nel caricamento del modello:", err);
            if (viewerRef.current) {
              const errorMessage = err.message.includes('Timeout') 
                ? 'Timeout nel caricamento del modello 3D'
                : `Errore nel caricamento del modello: ${err.message}`;
              
              viewerRef.current.innerHTML = `
                <div class="error-message">
                  <p>${errorMessage}</p>
                  <p>Prova a ricaricare o contatta il supporto se il problema persiste.</p>
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

  // Funzione migliorata per la chiusura
  const handleClose = () => {
    // Cancella il download se in corso
    cancelDownload();
    
    // Chiama la funzione di chiusura del parent
    onClose();
  };

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
          <div className="molecule-header-info">
            <p><strong>SMILES:</strong> {molecule.smiles}</p>
          </div>
          <div className="molecule-header-actions">
            {/* Pulsante Download XYZ */}
            {!loading && !error && molecule.modelPath && (
              <div className="download-section">
                <button 
                  className="download-xyz-button"
                  onClick={handleDownloadXYZ}
                  disabled={isDownloading}
                  title="Scarica file XYZ"
                >
                  {isDownloading ? (
                    <>
                      <div className="download-spinner"></div>
                      Download...
                      <button 
                        className="cancel-download-btn"
                        onClick={(e) => {
                          e.stopPropagation();
                          cancelDownload();
                        }}
                        title="Annulla download"
                      >
                        ✕
                      </button>
                    </>
                  ) : (
                    <>
                      📁 Scarica XYZ
                    </>
                  )}
                </button>
                
                {/* Messaggio di errore per il download */}
                {downloadError && (
                  <div className="download-error">
                    {downloadError}
                  </div>
                )}
              </div>
            )}
            <button className="close-button" onClick={handleClose}>×</button>
          </div>
        </div>

        <div className="molecule-viewer-content">
          {loading ? (
            <div className="loading-container">
              <div className="loading-spinner"></div>
              <p>Generazione del modello 3D in corso...</p>
              <button 
                className="cancel-generation-btn"
                onClick={handleClose}
                style={{
                  marginTop: '1rem',
                  padding: '0.5rem 1rem',
                  background: '#dc3545',
                  color: 'white',
                  border: 'none',
                  borderRadius: '4px',
                  cursor: 'pointer'
                }}
              >
                Annulla generazione
              </button>
            </div>
          ) : error ? (
            <div className="error-container">
              <p>{error}</p>
              <p>Impossibile generare il modello 3D per questa molecola.</p>
              {error.includes('cancellata') && (
                <p style={{ color: '#666', fontSize: '0.9em' }}>
                  La generazione è stata interrotta dall'utente.
                </p>
              )}
            </div>
          ) : (
            <>
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