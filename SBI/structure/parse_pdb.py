from ..data         import aminoacids3to1, nucleic2to1
from .chain         import Chain, ChainOfProtein, ChainOfNucleotide
from .header.Header import PDBHeader
import re


def read_PDB_header(pdbobject):

    PDBH = PDBHeader()
    for line in pdbobject.pdb_file.read():
        line = re.sub('\'', '\\\'', line)   # JIC
        # print line.strip()
        # CONTROL FOR AFTERWARDS
        if line.startswith('HEADER'):       # BASIC INFORMATION
            PDBH.add_header(line)
        elif line.startswith('OBSLTE'):
            PDBH.add_obslte(line)
        elif line.startswith('TITLE'):      # EXPERIMENT TITLE (ML)
            PDBH.add_title(line)
        elif line.startswith('EXPDTA'):     # EXPERIMENT TYPE (ML)
            PDBH.add_experiment(line, 'type')
        elif line.startswith('REMARK   2 RESOLUTION'):
            PDBH.add_experiment(line, 'resolution')
        elif line.startswith('REMARK   3   R VALUE'):
            PDBH.add_experiment(line, 'rfactor')
        elif line.startswith('REMARK   3   FREE R VALUE     '):
            PDBH.add_experiment(line, 'freeR')
        elif line.startswith('SPRSDE'):
            PDBH.add_deprecated(line)
        elif line.startswith('COMPND'):     # MOLECULE PARSING
            PDBH.add_molecule(line, 'COMPND', re.search(r'MOL_ID\:', line))
        elif line.startswith('SOURCE'):
            PDBH.add_molecule(line, 'SOURCE',  re.search(r'MOL_ID\:', line))
        elif line.startswith('KEYWDS'):     # KEYWORDS (ML)
            PDBH.add_keywords(line)
        elif line.startswith('DBREF '):     # REFERENCES TO OTHER DB
            PDBH.add_dbreference(line)
        elif line.startswith('REMARK 800 SITE_IDENTIFIER'):  # SITES
            PDBH.add_site(line, 'IDENTIFIER')
        elif line.startswith('REMARK 800'):
            PDBH.add_site(line, 'REMARK')
        elif line.startswith('SITE'):
            PDBH.add_site(line, 'SITE')
        elif line.startswith('HET   '):     # HETEROATOMS (ML)
            PDBH.add_hetatom(line, 'HET')
        elif line.startswith('HETNAM'):
            PDBH.add_hetatom(line, 'HETNAM')
        elif line.startswith('FORMUL'):
            PDBH.add_hetatom(line, 'FORMUL')
        elif line.startswith('HELIX '):     # SECONDARY STRUCTURE
            PDBH.add_secondary_structure(line, 'HELIX')
        elif line.startswith('SHEET '):
            PDBH.add_secondary_structure(line, 'SHEET')
        elif line.startswith('TURN '):
            PDBH.add_secondary_structure(line, 'TURN')
        elif line.startswith('REMARK 290   SMTRY'):  # MATRICES (ML)
            PDBH.add_simetry_matrix(line)
        elif line.startswith('REMARK 350 BIOMOLECULE'):
            PDBH.add_biomolecule(line)
        elif line.startswith('REMARK 350 APPLY THE FOLLOWING TO CHAINS'):
            PDBH.link_biomolecule(line)
        elif line.startswith('REMARK 350                    AND CHAINS'):
            PDBH.link_biomolecule(line)
        elif line.startswith('REMARK 350   BIOMT'):
            PDBH.add_biomolecule_matrix(line)
        elif line.startswith('SEQRES'):
            PDBH.add_sequence(line)
        elif line.startswith('ATOM') or line.startswith('HETATM'):
            break

    PDBH.process()

    pdbobject._header = PDBH


def read_PDB_file(pdbobject, biomolecule=False):
    """
    Process and load crystal data from a PDB formated file
    """
    read = True

    # Open takes into account the possibility of a gzipped file
    # pdb_fd = pdbobject.pdb_file.read()

    # To force the creation of a new structure
    old_chain = "OLD_CHAIN"

    for line in pdbobject.pdb_file.read():
        # In NMR, this will force to read only one model
        if line.startswith("ENDMDL"):
            old_chain = "OLD_CHAIN"
            pdbobject._NMR = True

        if line.startswith("ATOM") or line.startswith("HETATM"):
            chain = line[21:22].strip()
            if chain != old_chain:
                # The PDB assigned to the structure is the filename WITHOUT THE PATH
                if pdbobject.is_NMR or (not pdbobject.chain_exists(chain) or biomolecule):
                    chain_type = line[17:20].strip()
                    if chain_type in aminoacids3to1:
                        obj_chain = ChainOfProtein(pdb = pdbobject.pdb_file.first_prefix, chain = chain)
                        pdbobject._has_prot = True
                    elif chain_type in nucleic2to1:
                        obj_chain = ChainOfNucleotide(pdb = pdbobject.pdb_file.first_prefix, chain = chain)
                        pdbobject._has_nucl = True
                    else:
                        obj_chain = Chain(pdb = pdbobject.pdb_file.first_prefix, chain = chain)
                    pdbobject.add_chain(obj_chain, NMR = pdbobject.is_NMR)
                else:
                    putative_old_chain = pdbobject.get_chain_by_id(id = chain)
                    obj_chain = putative_old_chain

                old_chain = chain
                read = True

            if read:
                chain_type = line[17:20].strip()
                if not isinstance(obj_chain, ChainOfNucleotide) and not isinstance(obj_chain, ChainOfProtein):
                    if chain_type in nucleic2to1:
                        newobj_chain = ChainOfNucleotide(pdb = pdbobject.pdb_file.first_prefix, chain = obj_chain.chain)
                        pdbobject._has_nucl = True
                        newobj_chain.fuse(chain = obj_chain, lapl = True)
                        del(obj_chain)
                        obj_chain = newobj_chain
                        pdbobject._chains[pdbobject._get_chain_position_by_id(id = obj_chain.chain)] = obj_chain
                    if chain_type in aminoacids3to1:
                        newobj_chain = ChainOfProtein(pdb = pdbobject.pdb_file.first_prefix, chain = obj_chain.chain)
                        pdbobject._has_prot = True
                        newobj_chain.fuse(chain = obj_chain, lapl = True)
                        del(obj_chain)
                        obj_chain = newobj_chain
                        pdbobject._chains[pdbobject._get_chain_position_by_id(id = obj_chain.chain)] = obj_chain

                obj_chain.read_PDB_line(line)

        if line.startswith("TER"):
            obj_chain.is_term
            old_chain = "OLD_CHAIN"

    pdbobject.pdb_file.close()

    pdbobject._chain_id = set([c.chain for c in pdbobject.chains if not c.is_empty])
    pdbobject._chains   = [c for c in pdbobject.chains if not c.is_empty]


def read_PDB_line(thischain, line, keep_version = "A"):
    """
    Given a PDB-formated line, creates an atom to add to a new or a pre-existent residue

    @type  line: String
    @param line: PDB formated line

    @type  keep_version: String
    @param keep_version: Some residues have two versions, codified in front of the residue_type (line[16:17])
                         By default we keep the A version of doubles Aa, but it can be changed through parameters
    """
    isOK = set([" ", keep_version, thischain._residue_version])

    # This skips several versions of the same residue. It always will grab " " and keep_version
    # but in chains that are all B (1EA4_Z) this will help grab them
    if line[16:17] not in isOK:
        if len(thischain) == 0:
            thischain._residue_version = line[16:17]
        else:
            return

    residue_num = int(line[22:26].strip())
    residue_ver = line[26:27]

    # Check the need of creating a new residue
    if thischain.is_empty or (len(thischain._last_appended_list) > 0 and thischain._last_appended_list[-1].identifier != str(residue_num) + residue_ver):
        residue = thischain._new_Residue(number = residue_num, version = residue_ver, Rtype = line[17:20].strip(), mode = line[:6].strip())
        thischain.add_residue(residue)

    # Create the new atom
    x, y, z = [float(line[30+8*i:38+8*i]) for i in range(3)]
    if len(line) >= 55:
        extra = line[54:].rstrip()
    else:
        extra = ''
    atom = thischain._new_Atom(number = line[6:12].strip(), name = line[12:16].strip(), x = x, y = y, z = z,
                               extra  = extra)

    # Add the atom to the last residue (it applies the same if the residue has been created now or before)
    thischain._last_appended_list[-1].add_atom(atom)

# Abbreviations:
# --------------
# (ML) : MultiLine
