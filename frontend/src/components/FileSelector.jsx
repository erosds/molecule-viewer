import React, { useState, useEffect, useRef } from 'react';
import axios from 'axios';
import './FileSelector.css';

const FileSelector = ({ onSelectFile, onUploadFile, selectedFile }) => {
  const [availableFiles, setAvailableFiles] = useState([]);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState(null);
  const [uploading, setUploading] = useState(false);
  const fileInputRef = useRef(null);

  // Al montaggio, recupera l'elenco dei file CSV disponibili
  useEffect(() => {
    fetchFiles();
  }, []);

  const fetchFiles = async () => {
    setLoading(true);
    try {
      const response = await axios.get('/api/csv-files');
      setAvailableFiles(response.data);
      setError(null);
    } catch (err) {
      console.error("Errore nel recupero dei file:", err);
      setError("Impossibile recuperare la lista dei file");
    } finally {
      setLoading(false);
    }
  };

  const handleFileChange = (e) => {
    const selectedFileName = e.target.value;
    onSelectFile(selectedFileName);
  };

  const handleUploadClick = () => {
    fileInputRef.current.click();
  };

  const handleFileUpload = async (e) => {
    const file = e.target.files[0];
    if (!file) return;

    if (!file.name.endsWith('.csv')) {
      setError("Per favore, seleziona un file CSV");
      return;
    }

    setUploading(true);
    try {
      const success = await onUploadFile(file);
      if (success) {
        await fetchFiles();
        // Seleziona automaticamente il file appena caricato
        onSelectFile(file.name);
      }
    } finally {
      setUploading(false);
      // Reset del valore dell'input file
      e.target.value = null;
    }
  };

  return (
    <div className="file-selector">
      <div className="file-selector-container">
        <div className="file-selector-left">
          <label htmlFor="csv-file">Seleziona un file CSV:</label>
          {loading ? (
            <div className="file-selector-loading">Caricamento file disponibili...</div>
          ) : error ? (
            <div className="file-selector-error">{error}</div>
          ) : (
            <select 
              id="csv-file" 
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

        <div className="file-selector-right">
          {/* Input file nascosto, attivato dal pulsante */}
          <input 
            type="file" 
            ref={fileInputRef} 
            onChange={handleFileUpload} 
            style={{ display: 'none' }} 
            accept=".csv" 
          />
          
          <button 
            className="file-upload-button" 
            onClick={handleUploadClick}
            disabled={uploading}
          >
            {uploading ? 'Caricamento...' : 'Carica nuovo file CSV'}
          </button>
          
          <button 
            className="file-refresh-button"
            onClick={fetchFiles}
            disabled={loading}
          >
            ↻
          </button>
        </div>
      </div>
    </div>
  );
};

export default FileSelector;