from SBI              import SBIglobals
from .                import process_DBREF_line
from SBI.beans.JSONer import JSONer


class DBref(JSONer):

    '''
    {DBRef} includes the reference data to the PDB-CHAIN from:
    Known database-references include:

    Database name                            database (code in columns 27 - 32)
    ---------------------------------------------------------------------------
    BioMagResBank                            BMRB
    BLOCKS                                   BLOCKS
    European Molecular Biology Laboratory    EMBL
    GenBank                                  GB
    Genome Data Base                         GDB
    Nucleic Acid Database                    NDB
    PROSITE                                  PROSIT
    Protein Data Bank                        PDB
    Protein Identification Resource          PIR
    SWISS-PROT                               SWS
    TREMBL                                   TREMBL
    UNIPROT                                  UNP

    '''
    valid_references = {
        'BMRB'      : 'BioMagResBank',
        'BLOCKS'    : 'BLOCKS',
        'EMBL'      : 'European Molecular Biology Laboratory',
        'GB'        : 'GenBank',
        'GDB'       : 'Genome Data Base',
        'NDB'       : 'Nucleic Acid Database',
        'PROSIT'    : 'PROSITE',
        'PDB'       : 'Protein Data Bank',
        'PIR'       : 'Protein Identification Resource',
        'SWS'       : 'SWISS-PROT',
        'TREMBL'    : 'TREMBL',
        'UNP'       : 'UNIPROT'
    }

    def __init__(self, pdb_line):
        '''
        @type  pdb_line:     String
        @param pdb_pdb_line: Line of a PDB. Starts with DBREF

        @raise AttributeError if line does not start with DBREF
        '''
        if not pdb_line.startswith('DBREF '):
            SBIglobals.throw(self, '{0} cannot create DBref'.format(pdb_line))
        data          = self._process_line(pdb_line)

        self._pdb     = data[0]
        self._chain   = data[1]
        self._start   = data[2]
        self._end     = data[3]
        self._db      = data[4]
        self._ref     = data[5]

    #
    # ATTRIBUTES
    #
    @property
    def pdb(self):
        return self._pdb

    @property
    def chain(self):
        return self._chain

    @property
    def start(self):
        return int(self._start)

    @property
    def idxs(self):
        return self._start.index

    @property
    def end(self):
        return int(self._end)

    @property
    def idxe(self):
        return self._end.index

    @property
    def database(self):
        return self.valid_references[self._db]

    @property
    def reference(self):
        return self._ref

    #
    # BOOLEANS
    #
    @property
    def is_reference_to_db(self, db_minicode):
        db_minicode = db_minicode.upper()
        if not db_minicode in self.valid_references:
            line1 = '{0} is not a valid DBref code'.format(db_minicode)
            line2 = 'Available codes: {0}'.format(self.valid_references.keys())
            SBIglobals.throw(self, "\n".join([line1, line2]))

        #SWS and TREMBL are parts of UNP
        if db_minicode == 'UNP' and self._db in ['SWS', 'TREMBL']:
            return True
        return self._db == db_minicode

    #
    # PRIVATE
    #
    def _process_line(self, pdb_line):
        return process_DBREF_line(pdb_line)

    #
    # OUT FUNCTIONS
    #
    def as_dict(self):
        nobj = {}
        nobj['chain']     = self.chain
        nobj['db']        = self.database
        nobj['start']     = self.start
        nobj['idxs']      = self.idxs
        nobj['end']       = self.end
        nobj['idxe']      = self.idxe
        nobj['reference'] = self.reference
        return nobj

    def __repr__(self):
        return repr(self.json())
