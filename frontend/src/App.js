import React, { useState } from 'react';
import FileSelector from './components/FileSelector';
import MoleculeGrid from './components/MoleculeGrid';
import './App.css';

function App() {
  const [selectedFile, setSelectedFile] = useState(null);
  const [molecules, setMolecules] = useState([]);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);

  const handleFileSelect = (file) => {
    setSelectedFile(file);
    setLoading(true);
    setError(null);

    // Utilizziamo il nuovo endpoint per ottenere le molecole dal file CSV
    fetch(`/api/molecules/${file}`)
      .then(response => {
        if (!response.ok) {
          throw new Error('Impossibile caricare le molecole dal file CSV');
        }
        return response.json();
      })
      .then(parsedMolecules => {
        setMolecules(parsedMolecules);
        setLoading(false);
      })
      .catch(err => {
        setError(`Errore: ${err.message}`);
        setLoading(false);
      });
  };

  return (
    <div className="app">
      <header className="app-header">
        <h1>Visualizzatore Molecolare</h1>
      </header>
      <main className="app-content">
        <FileSelector onSelectFile={handleFileSelect} selectedFile={selectedFile} />
        
        {loading ? (
          <div className="loading">Caricamento molecole...</div>
        ) : error ? (
          <div className="error">{error}</div>
        ) : molecules.length > 0 ? (
          <MoleculeGrid molecules={molecules} />
        ) : (
          <div className="instructions">
            Seleziona un file CSV contenente strutture molecolari SMILES per iniziare
          </div>
        )}
      </main>
    </div>
  );
}

export default App;