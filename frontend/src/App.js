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

    // Parsing del file CSV
    fetch(`/csv/${file}`)
      .then(response => {
        if (!response.ok) {
          throw new Error('Impossibile caricare il file CSV');
        }
        return response.text();
      })
      .then(csvText => {
        try {
          const parsedMolecules = parseCSV(csvText);
          setMolecules(parsedMolecules);
        } catch (err) {
          setError(`Errore nel parsing del file CSV: ${err.message}`);
        } finally {
          setLoading(false);
        }
      })
      .catch(err => {
        setError(`Errore: ${err.message}`);
        setLoading(false);
      });
  };

  // Funzione per il parsing del CSV con attenzione alla colonna SMILES
  const parseCSV = (csvText) => {
    const lines = csvText.trim().split('\n');
    if (lines.length <= 1) {
      throw new Error('Il file CSV è vuoto o contiene solo intestazioni');
    }

    const headers = lines[0].split(',').map(header => header.trim().toLowerCase());
    const smilesIndex = headers.findIndex(header => header.includes('smiles'));

    if (smilesIndex === -1) {
      throw new Error('Nessuna colonna SMILES trovata nel file CSV');
    }

    const molecules = [];
    for (let i = 1; i < lines.length; i++) {
      const values = lines[i].split(',');
      if (values.length > smilesIndex) {
        // Aggiungiamo un'ID univoco basato sull'indice
        molecules.push({
          id: i - 1,
          smiles: values[smilesIndex].trim()
        });
      }
    }

    return molecules;
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