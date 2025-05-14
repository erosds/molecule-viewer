# Visualizzatore Molecolare

Applicazione per visualizzare strutture molecolari 2D e 3D a partire da file CSV contenenti stringhe SMILES.

## Caratteristiche

- Caricamento di file CSV contenenti stringhe SMILES
- Visualizzazione 2D delle molecole generate con RDKit
- Generazione di modelli 3D (file PDB) per le molecole selezionate
- Interfaccia utente reattiva e moderna

## Requisiti

### Backend
- Python 3.8+
- FastAPI
- RDKit
- Uvicorn
- Pillow

### Frontend
- Node.js 14+
- React 18+
- npm o yarn

## Installazione

### Backend

1. Crea un ambiente virtuale Python (consigliato):
   ```
   python -m venv venv
   ```

2. Attiva l'ambiente virtuale:
   - Windows: `venv\Scripts\activate`
   - macOS/Linux: `source venv/bin/activate`

3. Installa le dipendenze:
   ```
   pip install -r requirements.txt
   ```
   
   Nota: RDKit potrebbe richiedere passaggi aggiuntivi per l'installazione. Su alcuni sistemi, potrebbe essere più facile installarlo tramite conda:
   ```
   conda install -c conda-forge rdkit
   ```

4. Crea le cartelle necessarie per i file statici:
   ```
   mkdir -p public/csv public/molecules public/images
   ```

5. Aggiungi alcuni file CSV di esempio nella cartella `public/csv`

### Frontend

1. Naviga nella directory del frontend:
   ```
   cd frontend
   ```

2. Installa le dipendenze:
   ```
   npm install
   ```
   o
   ```
   yarn install
   ```

## Avvio dell'applicazione

### Backend

1. Dalla directory principale, avvia il server backend:
   ```
   cd backend
   python main.py
   ```
   
   Il server sarà disponibile all'indirizzo: http://localhost:8000

### Frontend

1. In un altro terminale, avvia il server di sviluppo React:
   ```
   cd frontend
   npm start
   ```
   o
   ```
   yarn start
   ```
   
   L'applicazione sarà disponibile all'indirizzo: http://localhost:3000

## Struttura del progetto

```
/
├── backend/
│   ├── main.py                # Server FastAPI
│   ├── molecule_utils.py      # Funzioni per la gestione delle molecole
│   ├── requirements.txt       # Dipendenze Python
│   └── public/
│       ├── csv/               # File CSV con dati molecolari
│       ├── molecules/         # Modelli 3D generati (PDB)
│       └── images/            # Immagini 2D generate (PNG)
│
└── frontend/
    ├── public/
    │   └── index.html
    ├── src/
    │   ├── App.js             # Componente principale
    │   ├── components/
    │   │   ├── FileSelector.js   # Selettore file CSV
    │   │   ├── MoleculeGrid.js   # Griglia di molecole
    │   │   └── MoleculeViewer.js # Visualizzatore 3D
    │   └── index.js
    └── package.json
```

## Utilizzo

1. Avvia sia il backend che il frontend come descritto sopra
2. Apri il browser all'indirizzo http://localhost:3000
3. Seleziona un file CSV dal menu a discesa (deve contenere una colonna SMILES)
4. Le molecole verranno visualizzate in una griglia con immagini 2D generate da RDKit
5. Fai clic su una molecola per aprire il visualizzatore 3D

## Note per lo sviluppo

- Il frontend comunica con il backend tramite API REST
- Le immagini 2D sono generate al volo utilizzando RDKit e memorizzate nella cache
- I modelli 3D sono generati solo quando richiesti e memorizzati nella cache
- Il proxy nel package.json del frontend reindirizza automaticamente le chiamate API al backend

## Risoluzione dei problemi

1. **Errore "No module named 'rdkit'"**:
   - Verifica che RDKit sia installato correttamente
   - Prova ad installarlo usando conda: `conda install -c conda-forge rdkit`

2. **Nessun file CSV disponibile**:
   - Aggiungi almeno un file CSV nella cartella `backend/public/csv`
   - Assicurati che il file CSV contenga una colonna con "smiles" nel nome

3. **Le immagini delle molecole non vengono caricate**:
   - Verifica che il backend sia in esecuzione
   - Controlla i log del backend per eventuali errori nella generazione delle immagini
   - Assicurati che i percorsi delle cartelle siano corretti