import React, { useState, useEffect } from 'react';
import './FileSelector.css';

const FileSelector = ({ onSelectFile, selectedFile, type = "main", label }) => {
  const [availableFiles, setAvailableFiles] = useState([]);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState(null);

  // Al montaggio, recupera l'elenco dei file CSV disponibili
  useEffect(() => {
    setLoading(true);
    fetch('/api/csv-files')
      .then(response => {
        if (!response.ok) {
          throw new Error('Impossibile recuperare la lista dei file');
        }
        return response.json();
      })
      .then(files => {
        setAvailableFiles(files);
        setLoading(false);
      })
      .catch(err => {
        setError(`Errore nel caricamento dei file: ${err.message}`);
        setLoading(false);
      });
  }, []);

  const handleFileChange = (e) => {
    const selectedFileName = e.target.value;
    onSelectFile(selectedFileName);
  };

  return (
    <div className="file-selector">
      <div className="file-selector-container">
        {loading ? (
          <div className="file-selector-loading">Caricamento file disponibili...</div>
        ) : error ? (
          <div className="file-selector-error">{error}</div>
        ) : availableFiles.length === 0 ? (
          <div className="file-selector-error">
            Nessun file CSV disponibile. Aggiungere file nella cartella /backend/public/csv
          </div>
        ) : (
          <select 
            id={`csv-file-${type}`} 
            value={selectedFile || ''} 
            onChange={handleFileChange}
            className="file-selector-dropdown"
          >
            <option value="">-- Seleziona un file --</option>
            {availableFiles.map(file => (
              <option key={file} value={file}>{file}</option>
            ))}
          </select>
        )}
      </div>
    </div>
  );
};

export default FileSelector;