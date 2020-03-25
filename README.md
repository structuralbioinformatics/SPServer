# Split-statistical Potentials Server (SPServer)

### 2020 - Structural Bioinformatics Group - Universitat Pompeu Fabra

This is the standalone version of **SPServer** (http://sbi.upf.edu/spserver), a web server to assess the quality of protein folds and protein-protein interactions using knowledge-based potentials (also known as statistical potentials).

## Getting Started

### Assessment of protein folds

```
python SPServerFold.py -f <PDB_PATH> -p <POT_TYPE> -o <OUT_PATH> -x <XML_NAME>
  -f: Input file with the PDB proteins. It can be a PDB/CIF file, a list of paths to PDB/CIF files or a TAR.GZ/ZIP.
  -p: Definition of atom-atom minimum distance. It can be CB or MIN.
  -o: Output directory.
  -x: Name of the XML output file (without the complete path and extension).
```

Example:

```
python SPServerFold.py
  -f 1ivo.pdb
  -p CB
  -o ./output
  -x 1ivo
```

### Assessment of protein-protein interactions


```
python SPServerPPI.py -i <INPUT_FILE> -r <CHAIN_RECEPTOR> -l <CHAIN_LIGAND> -s <PPI_SOURCE> -o <OUTPUT_PATH> -j <JOB_ID> -p <POT_TYPE>
  -i: Input file with the structures of the interactions. It can be a PDB/CIF file, a list of paths to PDB/CIF files or a TAR.GZ/ZIP.
  -r: Receptor chain ID. Only necessary if a PDB/CIF file has been introduced.
  -l: Ligand chain ID. Only necessary if a PDB/CIF file has been introduced.
  -s: Type of input. It can be "pdb_separated" (two PDB structures separated) or "pdb_together" (one PDB structure with 2 chains being the proteins that interact).
  -o: Output directory.
  -j: Job ID.
  -p: Definition of atom-atom minimum distance. It can be CB or MIN.
  -c: If the argument is introduced, it calculates the crashes between the interacting proteins.
```

Example:

```
python SPServerPPI.py
  -i  BAX_BID_reference.pdb
  -r  A
  -l  C
  -s  pdb_together
  -o  ./output
  -j  BAX_BID_reference
  -p  CB
  -c  
```

