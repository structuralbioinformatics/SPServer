import gzip, os, ftplib, sys
from BioLib.Structure.Structure import Structure
from BioLib.Structure.Residue import Residue
from BioLib.Structure.Atom import Atom

def read_pdb(pdb_file, one_model=True, from_string=False, merge_chains=False):
    '''
    Parses a PDB file
    @one_model: If the PDB contains more than one model, it returns only the first one
    @from_string: pdb_file is a PDB string instead of a file name. For example: read from a pipe
    @merge_chains: Merge the chains of the PDB into one structure, instead of returning a chain list
    '''
    models = []
    chains = []
    uniprot_ref = {}

    last_structure = '#####'
    base_residue = 0
    base_atom = 0

    if from_string:                pdb_fo = pdb_file.splitlines()
    elif pdb_file.endswith('.gz'): pdb_fo = gzip.open(pdb_file)
    else:                          pdb_fo = open(pdb_file)

    for line in pdb_fo:
        # Read Database references
        if line.startswith('DBREF'):
            if line[26:31].strip() == 'UNP':
                uniprot_ref[line[12]] = line[42:53].strip()
        # New model found
        if line.startswith('ENDMDL'):
            if one_model: 
                break
            else:
                models.append(chains)
                chains = []
                last_structure = '#####'
                base_residue = 0
                base_atom = 0
        # Read structure data
        if line.startswith('ATOM'):
            if len(line[26].strip()) > 0:
                continue
            chain = line[21:22]
            if last_structure != chain:
                if merge_chains and last_structure != '#####':
                    base_residue += structure.get_last_residue().get_num()
                    base_atom += structure.get_last_residue().get_atoms()[-1].get_num()
                    last_structure = chain
                else:
                    structure = Structure(os.path.basename(pdb_file).split('.')[0], chain)
                    chains.append(structure)
                    last_structure = chain  
            residueNum = int(line[22:26].strip())+base_residue
            ridueType = line[17:20].strip()
            atomNum = int(line[6:12].strip())#+base_atom
            atomType = line[13:17].strip()
            atomCoords = (float(line[27:38].strip()), float(line[38:46].strip()), float(line[46:54].strip()))
            if len(structure) == 0 or residueNum != structure.get_last_residue().get_num():
                residue = Residue(residueNum, ridueType)
                structure.add_residue(residue)
            atom = Atom(atomNum, atomType, atomCoords[0], atomCoords[1], atomCoords[2])
            residue.add_atom(atom)
            
    if not from_string:
        pdb_fo.close()

    if len(models) != 0:
        for model in models:
            for chain in model:
                if chain.get_chain() in uniprot_ref:
                    chain.set_uniprot_ref(uniprot_ref[chain.get_chain()])
        return models
    elif len(chains) > 1:
        for chain in chains:
            if chain.get_chain() in uniprot_ref:
                chain.set_uniprot_ref(uniprot_ref[chain.get_chain()])
        return chains
    else:
        if chains[0].get_chain() in uniprot_ref:
            chains[0].set_uniprot_ref(uniprot_ref[chains[0].get_chain()])
        return chains[0] 

def write_pdb(structures, pdb_name, multi_chain=False, multi_model=False):
    '''
    Creates a PDB file
    @structures = Structure objects to print
    @multi_chain = If true, the pdb has more than one chain                                                            
     - Structures must be a list of structure objects 
        --> (structure1, ..., structureN)
    @multi_model = If true, the pdb has more than one model
     - Structures must be a list of structure objects or chains (lists of structure objects) 
        --> (structure1, ..., structureN) or ((structure1, ..., structureN), ..., (structure1, ..., structureN))
    '''
    pdbFilefd = open(pdb_name, 'w')
    countModel = 0
    pdbFilefd.write('HEADER     PDB automatically created\n')
    if multi_model:
        for model in structures:
            pdbFilefd.write('MODEL%s'%countModel.rjust(9))
            if multi_chain:
                for chain in model:
                    pdbFilefd.write(chain.__str__())
            else:
                pdbFilefd.write(model.__str__())
            pdbFilefd.write('ENDMDL\n')
            countModel += 1
    elif multi_chain:
        for chain in structures:
            pdbFilefd.write(chain.__str__())
    else:
        pdbFilefd.write(structures.__str__())
    pdbFilefd.write('END')
    pdbFilefd.close()

def download_pdb(client_directory, pdb_list):
    '''
    Downloads PDBs from internet
    '''
    ftp_server = "ftp.ebi.ac.uk"
    server_directory = "/pub/databases/rcsb/pdb-remediated/data/structures/divided/pdb/"

    ftp = ftplib.FTP(ftp_server)
    ftp.login()
    for pdb in pdb_list:
        pdb = pdb.lower()
        ftp.cwd(server_directory+pdb[1:3])
        pdb_fo = open(os.path.join(client_directory, pdb+'.pdb.gz'), 'wb')
        try:
            ftp.retrbinary("RETR pdb"+pdb+".ent.gz", pdb_fo.write)
        except:
            raise RuntimeError("ERROR: Unable to download %s\n" % pdb)
        pdb_fo.close()
    ftp.quit()




