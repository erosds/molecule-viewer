import React, { useState } from 'react';
import FileSelector from './components/FileSelector';
import MoleculeGrid from './components/MoleculeGrid';
import './App.css';

function App() {
  const [selectedFile, setSelectedFile] = useState(null);
  const [molecules, setMolecules] = useState([]);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);

  const [referenceFile, setReferenceFile] = useState(null);
  const [referenceMolecules, setReferenceMolecules] = useState([]);
  const [newMolecules, setNewMolecules] = useState(0); // per tenere traccia del numero di molecole nuove

  const handleFileSelect = (file) => {
    setSelectedFile(file);
    setLoading(true);
    setError(null);

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

        // Calcola molecole nuove se il file di riferimento è già caricato
        if (referenceMolecules.length > 0) {
          calculateNewMolecules(parsedMolecules, referenceMolecules);
        } else {
          setNewMolecules(parsedMolecules.length); // Tutte sono nuove se non c'è riferimento
        }
      })
      .catch(err => {
        setError(`Errore: ${err.message}`);
        setLoading(false);
      });
  };

  const handleReferenceFileSelect = (file) => {
    setReferenceFile(file);

    if (!file) {
      setReferenceMolecules([]);
      return;
    }

    fetch(`/api/molecules/${file}`)
      .then(response => {
        if (!response.ok) {
          throw new Error('Impossibile caricare le molecole dal file di riferimento');
        }
        return response.json();
      })
      .then(parsedMolecules => {
        setReferenceMolecules(parsedMolecules);

        // Calcola molecole nuove se entrambi i set sono disponibili
        if (molecules.length > 0) {
          calculateNewMolecules(molecules, parsedMolecules);
        }
      })
      .catch(err => {
        console.error(`Errore: ${err.message}`);
        setReferenceMolecules([]);
      });
  };

  const calculateNewMolecules = (mainMolecules, refMolecules) => {
    if (!refMolecules.length) {
      setNewMolecules(mainMolecules.length); // Tutte sono nuove se non c'è riferimento
      return;
    }

    // Crea un set degli SMILES di riferimento per una ricerca efficiente
    const referenceSmiles = new Set(refMolecules.map(mol => mol.smiles));

    // Conta quante molecole del set principale non sono nel set di riferimento
    const newCount = mainMolecules.filter(mol => !referenceSmiles.has(mol.smiles)).length;
    setNewMolecules(newCount);
  };


  return (
  <div className="app">
    <header className="app-header">
      <h1>Visualizzatore Molecole</h1>
    </header>
    <main className="app-content">
      {/* Selettore del file di riferimento */}
      <FileSelector 
        onSelectFile={handleReferenceFileSelect} 
        selectedFile={referenceFile} 
        type="reference"
        label="Seleziona il file CSV di molecole di riferimento:"
      />
      
      {/* Selettore del file principale */}
      <FileSelector 
        onSelectFile={handleFileSelect} 
        selectedFile={selectedFile} 
        type="main"
        label="Seleziona un file CSV di molecole da visualizzare:"
      />
      
      {loading ? (
        <div className="loading">Caricamento molecole...</div>
      ) : error ? (
        <div className="error">{error}</div>
      ) : molecules.length > 0 ? (
        <>
          <MoleculeGrid molecules={molecules} newMoleculesCount={newMolecules} />
        </>
      ) : selectedFile ? (
        <div className="error">Nessuna molecola SMILES riconosciuta nel file selezionato</div>
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