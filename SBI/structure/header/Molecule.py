from . import process_COMPND_line
from . import process_SOURCE_line
from . import fix_chain

from SBI.beans.JSONer import JSONer
import re


class Molecule(JSONer):

    '''
    Includes the data from COMPND and SOURCE.
    Each one includes several tokens; those of COMPND:

    TOKEN           VALUE DEFINITION
    -------------------------------------------------------------------------
    MOL_ID          Numbers each component; also used in  SOURCE to associate
                    the information.
    MOLECULE        Name of the macromolecule.
    CHAIN           Comma-separated list of chain  identifier(s).
    FRAGMENT        Specifies a domain or region of the  molecule.
    SYNONYM         Comma-separated list of synonyms for  the MOLECULE.
    EC              The Enzyme Commission number associated  with the molecule.
                    If there is more than one EC number,  they are presented
                    as a comma-separated list.
    ENGINEERED      Indicates that the molecule was  produced using
                    recombinant technology or by purely  chemical synthesis.
    MUTATION        Indicates if there is a mutation.
    OTHER_DETAILS   Additional comments.

    And those of SOURCE:

    TOKEN                     VALUE  DEFINITION
    --------------------------------------------------------------------------
    MOL_ID                    Numbers each molecule. Same as appears in COMPND
    SYNTHETIC                 Indicates a  chemically-synthesized source.
    FRAGMENT                  A domain or  fragment of the molecule may be
                              specified.
    ORGANISM_SCIENTIFIC       Scientific name of the  organism.
    ORGANISM_COMMON           Common name of the  organism.
    ORGANISM_TAXID            NCBI Taxonomy ID number  of the organism.
    STRAIN                    Identifies the  strain.
    VARIANT                   Identifies the  variant.
    CELL_LINE                 The specific line of cells used in
                              the experiment
    ATCC                      American Type  Culture Collection tissue
                              culture  number.
    ORGAN                     Organized group of  tissues that carries on
                              a specialized function.
    TISSUE                    Organized group  of cells with a common
                              function and  structure.
    CELL                      Identifies the  particular cell type.
    ORGANELLE                 Organized structure  within a cell.
    SECRETION                 Identifies the secretion, such as  saliva, urine,
                              or venom,  from which the molecule was isolated.
    CELLULAR_LOCATION         Identifies the location  inside/outside the cell.
    PLASMID                   Identifies the plasmid  containing the gene.
    GENE                      Identifies the  gene.
    EXPRESSION_SYSTEM         Scientific name of the organism in  which the
                              molecule was expressed.
    EXPRESSION_SYSTEM_COMMON  Common name of the organism in which the
                              molecule was  expressed.
    EXPRESSION_SYSTEM_TAXID   NCBI Taxonomy ID of the organism  used as the
                              expression  system.
    EXPRESSION_SYSTEM_STRAIN  Strain of the organism in which  the molecule
                              was  expressed.
    EXPRESSION_SYSTEM_VARIANT Variant of the organism used as the
                              expression  system.
    EXPRESSION_SYSTEM_CELL_LINE  The specific line of cells used as  the
                                 expression  system.
    EXPRESSION_SYSTEM_ATCC_NUMBER  Identifies the ATCC number of the
                                   expression system.
    EXPRESSION_SYSTEM_ORGAN      Specific organ which expressed  the molecule.
    EXPRESSION_SYSTEM_TISSUE     Specific tissue which expressed  the molecule.
    EXPRESSION_SYSTEM_CELL       Specific cell type which  expressed the
                                 molecule.
    EXPRESSION_SYSTEM_ORGANELLE  Specific organelle which expressed
                                 the molecule.
    EXPRESSION_SYSTEM_CELLULAR_LOCATION  Identifies the location inside or
                                         outside the cell  which expressed
                                         the molecule.
    EXPRESSION_SYSTEM_VECTOR_TYPE   Identifies the type of vector used,  i.e.,
                                    plasmid,  virus, or cosmid.
    EXPRESSION_SYSTEM_VECTOR      Identifies the vector used.
    EXPRESSION_SYSTEM_PLASMID     Plasmid used in the recombinant experiment.
    EXPRESSION_SYSTEM_GENE        Name of the gene used in  recombinant
                                  experiment.
    OTHER_DETAILS                 Used to present  information on the
                                  source which
                                  is not  given elsewhere.
    '''
    def __init__(self, pdb):
        self._pdb      = pdb
        self._COMPND   = ''
        self._SOURCE   = ''
        self._chains   = []
        self._name     = ''
        self._ec       = set()
        self._taxid    = set()
        self._processd = False

    #
    # ATTRIBUTES
    #
    @property
    def pdb(self):
        return self._pdb

    @property
    def chains(self):
        return self._chains

    @property
    def name(self):
        return self._name

    @property
    def ec(self):
        return self._ec

    @property
    def taxid(self):
        return self._taxid

    #
    # BOOLEANS
    #
    @property
    def is_processed(self):
        return self._processd

    #
    # READ FUNCTIONS
    #
    def add_line(self, switch, line):
        if switch == 'COMPND':
            self._COMPND += process_COMPND_line(line)
        elif switch == 'SOURCE':
            self._SOURCE += process_SOURCE_line(line)

    #
    # PIVATE FUNCTIONS
    #
    def _parse(self):
        for field in self._SOURCE.split(';'):
            if field.startswith('ORGANISM_TAXID'):
                for tax in re.split(',', field.split(':')[1].strip()):
                    try:  # avoid add tag 'CHEMICALLY SYNTHESIZED' #4F4Z
                        int(tax.strip())
                        self.taxid.add(tax.strip())
                    except:
                        pass
                continue
            if field.strip() == 'SYNTHETIC: YES':
                if len(self.taxid) == 0:
                    self.taxid.add('32630')  # synthetic construct
                continue
        self._processd = True

    def _parse_cmpnd(self):
        # [EXCEPTION]. Found, at least, in pdb: 1TCR
        self._COMPND = re.sub('\\\;', ',', self._COMPND)

        for field in self._COMPND.split(';'):
            fdata = field.split(':')
            if field.startswith('CHAIN:'):
                self._chains = [fix_chain(f.strip()) for f in fdata[1].strip().split(',')]
                # [EXCEPTION]. Found in 2C14 and others
                if '' in self._chains:
                    self._chains.remove('')
                continue
            if field.startswith('MOLECULE:'):
                self._name = ":".join(fdata[1:]).strip()
                continue
            if field.startswith('EC:'):
                self._ec = set([x.strip() for x in fdata[1].split(',')])
                continue

    #
    # OUT FUNCTIONS
    #
    def as_dict(self):
        nobj = {}
        nobj['name']   = self.name
        nobj['chains'] = self.chains
        nobj['ec']     = list(self.ec)
        nobj['taxid']  = list(self.taxid)
        return nobj

    def __repr__(self):
        return self.json()
