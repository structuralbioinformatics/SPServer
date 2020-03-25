import datetime

from . import process_HEADER_line
from . import process_OBSLTE_line
from . import process_TITLE_line
from . import process_SUPERSEEDED_line
from . import process_MOLKEY
from . import process_KEYWRD_line
from . import process_SITE_IDENTIFIER_line
from . import process_REMARK800_line
from . import process_SITE_lines
from . import process_HETEROATOM_lines
from . import process_SECONDARYSTRUCTURE_lines
from . import process_MATRIX_lines
from . import process_REMARK350_lines
from . import process_SEQRES_lines

from .Experiment         import Experiment
from .Molecule           import Molecule
from .DBreference        import DBref
from .Site               import Site
from .HeteroAtom         import Hetero
from .SecondaryStructure import HelixInfo
from .SecondaryStructure import SheetInfo
from .SecondaryStructure import TurnInfo
from .BioMolecule        import BioMolecule
from SBI.beans.JSONer    import JSONer


class PDBHeader(JSONer):
    # TODO: search hetero by name

    def __init__(self):

        self._pdb        = None     # PDB ID
        self._date       = None     # Deposition Date
        self._obslte     = None     # Retirement Date
        self._replaces   = []       # PDB that replaced this
        self._header     = None     # Classification
        self._title      = None     # Experiment Tilte
        self._keywords   = ''       # Words of interest

        self._experiment = None     # Crystalographic Data
        self._deprecated = []       # Previous PDBs substituted by this one

        self._chaindict  = {}       # Dictionary of chains
                                    # Used to access to data by chain
        self._workchain  = None     # Working chain (full private)

        self._molecules  = {}       # Molecules in the PDB
        self._sites      = {}       # Sites of the PDB
        self._hetero     = {}       # Heteroatoms in the PDB
        self._dbrefs     = []       # References to other DB

        self._secstr     = []       # Secondary Structures

        self._symmetryM  = BioMolecule(0)  # Symmetry Matrices
        self._biomolecM  = []              # BioMolecule Matrices

    #
    # ATTRIBUTES
    #
    @property
    def pdb(self):
        return self._pdb if self._pdb is not None else ''

    @property
    def date(self):
        return self._date if self._date is not None else ''

    @property
    def date_as_date(self):
        if self._date is None:
            return None
        return datetime.date(*[int(x) for x in self.date.split('-')])

    @property
    def obsolete_date(self):
        return self._obslte if self._obslte is not None else ''

    @property
    def obsolete_date_as_date(self):
        if self._obslte is None:
            return None
        return datetime.date(*[int(x) for x in self.obsolete_date.split('-')])

    @property
    def replaced(self):
        return self._replaces

    @property
    def header(self):
        return self._header if self._header is not None else ''

    @property
    def title(self):
        return self._title if self._title is not None else ''

    @property
    def keywords(self):
        return self._keywords

    @property
    def experiment(self):
        return self._experiment

    @property
    def xpdta(self):
        return self.experiment.xpdta

    @property
    def resolution(self):
        return self.experiment.resolution

    @property
    def rfactor(self):
        return self.experiment.rfactor

    @property
    def freeR(self):
        return self.experiment.freeR

    @property
    def deprecated(self):
        return self._deprecated

    @property
    def molecules(self):
        return self._molecules

    @property
    def chains(self):
        return self._chaindict.keys()

    @property
    def sites(self):
        return self._sites

    @property
    def hetero(self):
        return self._hetero

    @property
    def dbrefs(self):
        return self._dbrefs

    @property
    def secondary_structures(self):
        return self._secstr

    @property
    def symmetry_matrix(self):
        return self._symmetryM

    @property
    def biomolecules(self):
        return self._biomolecM

    @property
    def sequences(self):
        data = {}
        for chain in self._chaindict:
            if 'SEQ' in self._chaindict[chain]:
                data[chain] = self._chaindict[chain]['SEQ']
        return data

    #
    # BOOLEANS
    #
    @property
    def are_molecules_processed(self):
        if len(self._molecules) > 0:
            k = self._molecules.keys()
            return self._molecules[k[0]].is_processed
        else:
            return True  # To avoid call the function

    #
    # READ FUNCTIONS
    #
    def add_header(self, line):
        data = process_HEADER_line(line)
        self._pdb    = data[2]
        self._header = data[0]
        self._date   = data[1]

    def add_obslte(self, line):
        data = process_OBSLTE_line(line)
        if data[0] != '':
            self._obslte = data[0]
        for x in data[1:]:
            if x.upper() != 'NONE':
                self._replaces.append(x)

    def add_title(self, line):
        data = process_TITLE_line(line)
        if self._title is None:
            self._title = data
        else:
            self._title += ' ' + data

    def add_experiment(self, line, switch):
        if switch == 'type':
            if self._experiment is None:
                self._experiment = Experiment(line)
            else:
                self._experiment.update_experiment(line)
        elif switch == 'resolution':
            self._experiment.add_resolution(line)
        elif switch == 'rfactor':
            self._experiment.add_rfactor(line)
        elif switch == 'freeR':
            self._experiment.add_freer(line)

    def add_deprecated(self, line):
        self._deprecated.extend(process_SUPERSEEDED_line(line))

    def add_molecule(self, line, switch, initkey):
        if initkey:
            self._workchain = process_MOLKEY(line)
            if switch == 'COMPND':
                self._molecules[self._workchain] = Molecule(self.pdb)
            elif switch == 'SOURCE':
                now_molecule = self._molecules[self._workchain]
                now_molecule._parse_cmpnd()
                for chain in now_molecule.chains:
                    self._chaindict.setdefault(chain, {'MOL': now_molecule})
        else:
            if self._workchain is None:  # Deprecated PDBs like 1RS3
                self._workchain = 1
                self._molecules[self._workchain] = Molecule(self.pdb)
            else:
                self._molecules[self._workchain].add_line(switch, line)

    def add_keywords(self, line):
        self._keywords += process_KEYWRD_line(line)

    def add_dbreference(self, line):
        self._dbrefs.append(DBref(line))
        chain = self.dbrefs[-1].chain
        if chain in self._chaindict:
            self._chaindict[chain].setdefault('DBREF', [])
        else:
            self._chaindict.setdefault(chain, {'DBREF': []})
        self._chaindict[chain]['DBREF'].append(self.dbrefs[-1])

    def add_site(self, line, switch):
        if switch == 'IDENTIFIER':
            self._workchain = process_SITE_IDENTIFIER_line(line)
            self._sites[self._workchain] = Site(self._workchain)
        elif switch == 'REMARK':
            data = process_REMARK800_line(line)
            if self._workchain in self._sites:
                self._sites[self._workchain].add_remark(data)
        elif switch == 'SITE':
            data = process_SITE_lines(line)
            self._sites.setdefault(data[0], Site(data[0]))
            self._sites[data[0]].add_residues(data[1:])

    def add_hetatom(self, line, switch):
        data = process_HETEROATOM_lines(line)
        if data[0] != 'HOH':   # SKIP WATERS
            if switch == 'HET':
                self._hetero.setdefault(data[0], Hetero(data[0], data[1]))
            elif switch == 'HETNAM':
                if data[0] in self._hetero:
                    self._hetero[data[0]].add_name(data[1])
            elif switch == 'FORMUL':
                if data[0] in self._hetero:
                    self._hetero[data[0]].add_formula(data[1])

    def add_secondary_structure(self, line, switch):
        data = process_SECONDARYSTRUCTURE_lines(line)
        if switch == 'HELIX':
            self._secstr.append(HelixInfo(data))
        elif switch == 'SHEET':
            self._secstr.append(SheetInfo(data))
        elif switch == 'TURN':
            self._secstr.append(TurnInfo(data))

    def add_simetry_matrix(self, line):
        data = process_MATRIX_lines(line)
        if int(data[0]) == 1:
            self.symmetry_matrix.new_matrix()
        self.symmetry_matrix.update_last_matrix(*data)

    def add_biomolecule(self, lines):
        data = process_REMARK350_lines(lines)
        self._biomolecM.append(BioMolecule(data))

    def link_biomolecule(self, lines):
        data = process_REMARK350_lines(lines)
        data = [x.strip() for x in data.split(',')]
        if len(self._biomolecM) == 0:
            self.add_biomolecule(':1')
        self._biomolecM[-1].chains = data

    def add_biomolecule_matrix(self, lines):
        data = process_MATRIX_lines(lines)
        if int(data[0]) == 1:
            self._biomolecM[-1].new_matrix()
        self._biomolecM[-1].update_last_matrix(*data)

    def add_sequence(self, line):
        data = process_SEQRES_lines(line)
        self._chaindict[data[0]].setdefault('SEQ', '')
        self._chaindict[data[0]]['SEQ'] += data[1]

    #
    # FINAL PROCESSING FUNCTION
    #
    def process(self):
        # Final processing of Molecule data
        if not self.are_molecules_processed:
            for m in self._molecules:
                self._molecules[m]._parse()
        # Processing keywords
        self._keywords = list(set([x.strip() for x in self._keywords.split(',')]))
        if len(self._keywords) == 1 and self._keywords[0] == '':
            self._keywords = []
        # Assigning sites to chain dictionary
        for k, v in self.sites.iteritems():
            for c in v.chains:
                if c not in self._chaindict:
                    if len(self._molecules) > 0:
                        nm = max([int(x) for x in self._molecules.keys()]) + 1
                    else:
                        nm = 1
                    self._molecules[nm] = Molecule(self.pdb)
                    self._chaindict.setdefault(c, {'MOL': self._molecules[nm]})
                    self._chaindict[c]['MOL']._chains.append(c)
                elif 'MOL' not in self._chaindict[c]:
                    self._chaindict[c]['MOL'] = Molecule(self.pdb)
                self._chaindict[c].setdefault('SITES', [])
                self._chaindict[c]['SITES'].append(v)
        # Assigning heteroatoms to chain dictionary
        for k, v in self.hetero.iteritems():
            for c in v.chain:
                if c not in self._chaindict:
                    if len(self._molecules) > 0:
                        nm = max([int(x) for x in self._molecules.keys()]) + 1
                    else:
                        nm = 1
                    self._molecules[nm] = Molecule(self.pdb)
                    self._chaindict.setdefault(c, {'MOL': self._molecules[nm]})
                    self._chaindict[c]['MOL']._chains.append(c)
                self._chaindict[c].setdefault('HETATM', [])
                self._chaindict[c]['HETATM'].append(v)
        # Assigning structures to chain dictionary
        for ss in self.secondary_structures:
            if ss.chain not in self._chaindict:
                if len(self._molecules) > 0:
                        nm = max([int(x) for x in self._molecules.keys()]) + 1
                else:
                    nm = 1
                self._molecules[nm] = Molecule(self.pdb)
                self._chaindict.setdefault(ss.chain, {'MOL': nm})
            self._chaindict[ss.chain].setdefault('SSTRUC', [])
            self._chaindict[ss.chain]['SSTRUC'].append(ss)
        # Assigning all chains to Symmetry Matrix
        self._symmetryM.chains = self._chaindict.keys()

    #
    # FUNCTIONS
    #
    def get_secondary_structures_definition_by(self, chain, start, end):
        ss = []
        for ssdef in self._chaindict[ss.chain]['SSTRUC']:
            if ssdef.defines(chain, start, end):
                ss.append(ssdef)
        return ss

    def get_chainheader(self, chain):
        if chain not in self.chains:
            raise KeyError('No chain {0} in the PDB header'.format(chain))
        return ChainHeader(self._chaindict[chain], self.experiment,
                           self.pdb,               chain)

    def has_hetero(self, heteroID):
        return heteroID in self.hetero.keys()

    def get_hetero_by_name_like(self, hetero_name):
        myheteros = []
        for h in self.hetero:
            if hetero_name.upper() in self.hetero[h].name:
                myheteros.append(self.hetero[h])
        return myheteros

    #
    # OUT FUNCTIONS
    #
    def as_dict(self):
        data = {'chains':    sorted(self.chains),  'date':      self.date,
                'classification': self.header,     'title':     self.title,
                'deprecated':     self.deprecated, 'id':        self.pdb,
                'keywords':       self.keywords,   'molecules': [],
                'experiment':     {},              'dbrefs':    [],
                'sites':          [],              'hetatm':    [],
                'sstruct':        [],              'sequences': self.sequences,
                'symmetry':       self.symmetry_matrix.as_dict(),
                'biomolecules':   [x.as_dict() for x in self.biomolecules],
                'obsolete_date':  self.obsolete_date,
                'replaces':       self.replaced
                }

        for k, x in self.molecules.iteritems():
            data['molecules'].append(x.as_dict())

        if self.experiment is not None:
            data['experiment'] = self.experiment.as_dict()

        if len(self.dbrefs) > 0:
            data['dbrefs']  = [x.as_dict() for x in self.dbrefs]

        if len(self.sites) > 0:
            data['sites']   = [x.as_dict() for k, x in self.sites.iteritems()]

        if len(self.hetero) > 0:
            data['hetatm']  = [x.as_dict() for k, x in self.hetero.iteritems()]

        if len(self.secondary_structures) > 0:
            data['sstruct'] = [x.as_dict() for x in self.secondary_structures]

        return data

    def __repr__(self):
        return self.json(pretty=True)


class ChainHeader(JSONer):
    def __init__(self, datadict, experiment, pdb, chain):
        self._data = datadict
        self._exp  = experiment
        self._pdb  = pdb
        self._chn  = chain

    #
    # ATTRIBUTES
    #
    @property
    def pdb(self):
        return self._pdb if self._pdb is not None else ''

    @property
    def chain(self):
        return self._chn

    @property
    def molecule(self):
        return self._data['MOL']

    @property
    def sequence(self):
        if 'SEQ' in self._data:
            return self._data['SEQ']
        else:
            return ''

    @property
    def sites(self):
        if 'SITES' not in self._data:
            return []
        return self._data['SITES']

    @property
    def dbreferences(self):
        if 'DBREF' not in self._data:
            return []
        return self._data['DBREF']

    @property
    def hetero(self):
        if 'HETATM' not in self._data:
            return []
        return self._data['HETATM']

    @property
    def secondary_structures(self):
        if 'SSTRUC' not in self._data:
            return []
        return self._data['SSTRUC']

    @property
    def xpdta(self):
        return self._exp.xpdta

    @property
    def resolution(self):
        return self._exp.resolution

    @property
    def rfactor(self):
        return self._exp.rfactor

    @property
    def freeR(self):
        return self._exp.freeR

    #
    # FUNCTIONS
    #
    def has_hetero(self, heteroID):
        for h in self.hetero:
            if h.restype == h:
                return True
        return False

    #
    # OUT FUNCTIONS
    #
    def as_dict(self):
        data = {'id':       self.pdb, 'chain':    self.chain,
                'experiment':     {}, 'dbrefs':   [],
                'sites':          [], 'hetatm':   [],
                'secstruct':      [], 'molecule': self.molecule.as_dict()}

        if self._exp is not None:
            data['experiment'] = self._exp.as_dict()

        if len(self.dbreferences) > 0:
            data['dbrefs']     = [x.as_dict() for x in self.dbreferences]

        if len(self.sites) > 0:
            data['sites']      = [x.as_dict() for x in self.sites]

        if len(self.hetero) > 0:
            data['hetatm']     = [x.as_dict() for x in self.hetero]

        if len(self.secondary_structures) > 0:
            for x in self.secondary_structures:
                data['secstruct'].append(x.as_dict())

        return data

    def __repr__(self):
        return self.json(pretty=True)
