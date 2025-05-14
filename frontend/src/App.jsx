import React, { useState } from 'react';
import FileSelector from './components/FileSelector';
import MoleculeGrid from './components/MoleculeGrid';
import './App.css';
import axios from 'axios';

function App() {
  const [selectedFile, setSelectedFile] = useState(null);
  const [molecules, setMolecules] = useState([]);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);

  const handleFileSelect = async (file) => {
    setSelectedFile(file);
    setLoading(true);
    setError(null);

    try {
      // Utilizziamo axios per una migliore gestione degli errori
      const response = await axios.get(`/api/csv/${file}`);
      setMolecules(response.data);
    } catch (err) {
      console.error("Errore nel caricamento del file:", err);
      
      // Gestione migliorata degli errori
      let errorMessage = "Errore nel caricamento del file CSV";
      if (err.response) {
        // Errore proveniente dal server
        errorMessage = `Errore ${err.response.status}: ${err.response.data.detail || err.response.statusText}`;
      } else if (err.request) {
        // Nessuna risposta ricevuta
        errorMessage = "Nessuna risposta dal server. Verificare la connessione di rete.";
      }
      
      setError(errorMessage);
    } finally {
      setLoading(false);
    }
  };

  const handleFileUpload = async (file) => {
    setLoading(true);
    setError(null);

    try {
      const formData = new FormData();
      formData.append('file', file);

      await axios.post('/api/upload-csv', formData, {
        headers: {
          'Content-Type': 'multipart/form-data'
        }
      });

      // Aggiorna la lista dei file dopo il caricamento
      return true;
    } catch (err) {
      console.error("Errore nel caricamento del file:", err);
      
      let errorMessage = "Errore nel caricamento del file CSV";
      if (err.response) {
        errorMessage = `Errore ${err.response.status}: ${err.response.data.detail || err.response.statusText}`;
      } else if (err.request) {
        errorMessage = "Nessuna risposta dal server. Verificare la connessione di rete.";
      }
      
      setError(errorMessage);
      return false;
    } finally {
      setLoading(false);
    }
  };

  return (
    <div className="app">
      <header className="app-header">
        <h1>Visualizzatore Molecolare</h1>
      </header>
      <main className="app-content">
        <FileSelector 
          onSelectFile={handleFileSelect} 
          onUploadFile={handleFileUpload}
          selectedFile={selectedFile} 
        />
        
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