# Visualizzatore Molecolare 3D

Un'applicazione web per visualizzare strutture molecolari in 3D a partire da notazioni SMILES.

## Caratteristiche

- Interfaccia utente intuitiva e reattiva basata su React
- Caricamento e selezione di file CSV contenenti strutture SMILES
- Generazione di modelli 3D delle molecole con RDKit
- Visualizzazione 3D interattiva tramite 3Dmol.js
- Generazione di immagini SVG per le molecole
- Visualizzazione di proprietà molecolari

## Struttura del Progetto

```
molecule-viewer/
├── frontend/               # Applicazione React
│   ├── public/             # File statici
│   └── src/                # Codice sorgente
│       ├── components/     # Componenti React
│       └── ...
└── backend/                # API FastAPI
    ├── main.py             # Server principale
    ├── molecule_utils.py   # Utilità per la gestione delle molecole
    └── ...
```

## Requisiti

### Frontend
- Node.js >= 14
- npm >= 7

### Backend
- Python >= 3.8
- FastAPI
- RDKit
- Uvicorn

## Installazione

1. Clonare il repository:
   ```bash
   git clone https://github.com/tuousername/molecule-viewer.git
   cd molecule-viewer
   ```

2. Configurare l'ambiente virtuale per il backend:
   ```bash
   cd backend
   python -m venv venv
   source venv/bin/activate  # Per Linux/Mac
   # oppure
   venv\Scripts\activate  # Per Windows
   pip install -r requirements.txt
   cd ..
   ```

3. Installare le dipendenze del frontend:
   ```bash
   cd frontend
   npm install
   cd ..
   ```

4. Installare le dipendenze globali del progetto:
   ```bash
   npm install
   ```

## Esecuzione

È possibile avviare l'applicazione in due modi:

### Avvio Simultaneo di Frontend e Backend
Utilizzare lo script npm nella directory principale:

```bash
npm start
```

Questo avvierà sia il backend che il frontend contemporaneamente.

### Avvio Separato

1. Avviare il backend:
   ```bash
   cd backend
   source venv/bin/activate  # Per Linux/Mac o venv\Scripts\activate per Windows
   python main.py
   ```

2. In un altro terminale, avviare il frontend:
   ```bash
   cd frontend
   npm start
   ```

L'applicazione sarà disponibile all'indirizzo `http://localhost:3000`.

## Utilizzo

1. Caricare un file CSV contenente strutture molecolari in formato SMILES
   - Il file deve avere almeno una colonna denominata "SMILES" (case-insensitive)
   - Opzionalmente può avere colonne "Nome" e "Formula"

2. Selezionare un file dalla lista
   - L'applicazione visualizzerà le molecole trovate nel file

3. Fare clic su una molecola per visualizzarla in 3D
   - Si aprirà un visualizzatore interattivo che permette di:
     - Ruotare la molecola
     - Ingrandire/rimpicciolire
     - Visualizzare i dettagli della molecola

## Soluzione dei Problemi

### Errore nel caricamento delle dipendenze

Se si riscontrano problemi nell'installazione di RDKit, può essere utile utilizzare Conda:

```bash
conda create -n rdkit-env python=3.9
conda activate rdkit-env
conda install -c conda-forge rdkit
pip install fastapi uvicorn python-multipart
```

### Il visualizzatore 3D non funziona

- Assicurarsi che jQuery e 3Dmol.js siano correttamente caricati
- Verificare la console del browser per errori JavaScript
- Controllare che il backend sia in esecuzione e risponda alle richieste

### Formato del file CSV non riconosciuto

- Verificare che il file CSV contenga una colonna denominata "SMILES" (o una variante come "smiles", "SMILE", ecc.)
- Controllare l'encoding del file (UTF-8 è consigliato)
- Assicurarsi che il delimitatore sia la virgola

## Licenza

MIT