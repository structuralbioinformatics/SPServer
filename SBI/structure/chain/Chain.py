"""
Chain

author: jbonet
date:   02/2013

@oliva's lab
"""

import warnings
import copy
import numpy
import re

from ..atom    import Atom
from ..residue import Residue, ResidueOfAminoAcid, ResidueOfNucleotide
from SBI.data  import aminoacids_main1
from SBI       import SBIglobals


class Chain(object):
    """
    A {Chain} represents a collection of {Residue}s
    """
    chaintype   = 'O'
    atomtype    = Atom
    resitype    = Residue
    dictitype   = {}

    def __init__(self, pdb, chain):
        """
        @type  pdb: String
        @param pdb: PDB identifier

        @type  chain: String
        @param chain: Chain identifier
        """

        if chain == "":
            warnings.warn("Possible old format PDB. No-chain is transformed to chain 'Z'")
            chain = "Z"

        # IDENTIFIERS
        self._pdb   = re.sub(r'^PDB', '', pdb.upper())
        self._chain = chain
        self._altID = None

        # LISTS
        self._structure = []
        self._hetero    = []
        self._HOH       = []

        self._chaindict = {}

        # PROPERTIES
        self._atom_length = None

        # RESIDUE VERSION ON READ
        self._residue_version = None

        # TERM
        self._term = False

        # CONTROL
        self._last_appended_list = []
        self._numeric_version    = 0  # pdb: 1CJ0

    #
    # ATTRIBUTES
    #
    @property
    def globalID(self):
        """
        {Chain} global identifier
        @rtype: String
        """
        return "{0}_{1}".format(self._pdb, self._chain)

    @property
    def pdb(self):
        """
        PDB name
        @rtype: String (uppercase)
        """
        return self._pdb

    @property
    def chain(self):
        """
        Chain identifier
        @rtype: String
        """
        return self._chain

    @chain.setter
    def chain(self, value):
        """
        Changes the chain assigned to the {Chain}
        @type value: String
        """
        self._chain = value

    @property
    def waters(self):
        """
        Lists, if any, the waters in the {Chain}
        @rtype: List of {Atom}s with type HOH
        """
        return self._HOH

    @property
    def heteroatoms(self):
        """
        List, if any, of NOT AMINOACID/NUCLEOTIDE heteroatoms
        @rtype: List of {Atom}s
        """
        return self._hetero

    @property
    def full_chain(self):
        return self._structure + self._hetero + self._HOH

    @property
    def alternative_id(self):
        """
        Returns an alternative ID that can be added in some processes
        @rtype: String
        """
        return self._altID

    @alternative_id.setter
    def alternative_id(self, value):
        """
        Sets an alternative ID
        @typr value: String
        """
        self._altID = value

    @property
    def atom_length(self):
        """
        Total number of atoms of the {Chain}
        @rtype: Integer
        """
        if self._atom_length is not None:
            return self._atom_length
        else:
            self._atom_length = sum(map(len, self._all_residues))
            return self._atom_length

    @property
    def last_residue(self):
        """
        Returns the GLOBAL last residue
        @rtype: {Residue}
        """
        return self._all_residues[-1]

    @property
    def first_structure(self):
        return self._structure[0] if len(self._structure) > 0 else self._hetero[0]

    @property
    def last_structure(self):
        return self._structure[-1] if len(self._structure) > 0 else self._hetero[-1]

    @property
    def is_term(self):
        """
        Mark that we have read untill chain TERM, from here everything will count as heteroatoms
        (only for original PDBs)
        """
        self._term = True

    #
    # COMPLEX GETTERS & SETTERS
    #
    def get_residue_by_identifier(self, identifier):
        if len(self._chaindict) == 0:
            for r in self.full_chain:
                self._chaindict[r.identifier.strip()] = r
        try:
            r = self._chaindict[str(identifier.strip())]
        except KeyError as e:
            SBIglobals.alert('error', self, '{0} is not a valid identifier for chain {1}'.format(identifier, self.globalID), e)
        else:
            return r

    def get_residues_by_type(self, restype):
        for r in self.full_chain:
            if r.type == restype:
                yield r

    def add_residue(self, residue):
        """
        Appends a new {Residue} to the list

        @type residue: {Residue}
        """
        if isinstance(residue, ResidueOfAminoAcid) or isinstance(residue, ResidueOfNucleotide):
            self._structure.append(residue)
            self._last_appended_list = self._structure
        elif residue.type == 'HOH' or residue.type == 'DOD':
            self._HOH.append(residue)
            self._last_appended_list = self._HOH
        else:
            self._hetero.append(residue)
            self._last_appended_list = self._hetero

    def join(self, chain):
        self._hetero.extend(chain._hetero)
        self._structure.extend(chain._structure)
        self._HOH.extend(chain._HOH)

    #
    # PRIVATE ATTRIBUTES
    #
    @property
    def _all_residues(self):
        """
        Returns all the 3 types of residues
        @rtype: List of {Residue}s
        """
        return self._structure + self._hetero + self._HOH

    def _all_atoms_coordinates(self, structure = True, hetero = True, water = True, geometric_center = False):
        all_coord = numpy.array(numpy.zeros(3))
        res_array = []
        if structure:
            res_array = res_array + self._structure
        if hetero:
            res_array = res_array + self._hetero
        if water:
            res_array = res_array + self._HOH
        for residue in res_array:
            if not geometric_center:
                all_coord = numpy.vstack((all_coord, residue.all_coordinates))
            else:
                all_coord = numpy.vstack((all_coord, residue.geometric_center))
        all_coord = numpy.delete(all_coord, 0, 0)
        return all_coord

    def _get_structure_array_coordinate(self, identifier):
        if isinstance(identifier, int):
            identifier = str(identifier + ' ')
        for i in range(len(self._structure)):
            if self._structure[i].identifier == identifier:
                return i
    #
    # METHODS
    #

    def dehydrate(self):
        """
        Removes all the waters
        """
        del(self._HOH)
        self._HOH = []

    def extract(self, init = None, end = None, backbone = False):
        """
        returns a section of a Chain
        """
        new_chain = self.__class__(self.pdb, self.chain)
        try:
            init_n, init_i = int(init), ' '
        except:
            init_n, init_i = int(init[:-1]), init[-1]
        try:
            end_n,  end_i  = int(end),  ' '
        except:
            end_n,  end_i  = int(end[:-1]),  end[-1]

        capture = False
        for residue in self._structure:
            if residue.number == init_n and residue.version == init_i:
                capture = True
            if capture:
                new_chain.add_residue(residue)
            if residue.number == end_n and residue.version == end_i:
                capture = False
                break

        nchain = new_chain.duplicate()
        if backbone:
            for residue in nchain.aminoacids:
                residue.remove_side_chain()

        return nchain

    def rotate(self, matrix):
        """
        Rotates each {Residue} according to a given matrix

        @type matrix: numpy.matrix
        """
        for residue in self._all_residues:
            residue.rotate(matrix = matrix)

    def translate(self, vector):
        """
        Translates each {Residue} according to a translational vector

        @type vector: numpy.array
        """
        for residue in self._all_residues:
            residue.translate(vector = vector)

    def reposition(self, matrix = None, vector = None):
        """
        Rotates and Translates each {Residue} according to a matrix and a translational vector

        @type matrix: numpy.matrix

        @type vector: numpy.array
        """
        if matrix is None:
            matrix = numpy.identity(3, float)
        if vector is None:
            vector = numpy.zeros(3, float)
        for residue in self._all_residues:
            residue.reposition(matrix = matrix, vector = vector)

    def translate_onto_origin(self):
        self.translate(vector=(self.geometric_center()) * -1)

    def geometric_center(self, structure = True, hetero = True, water = True, by_residue = False):

        all_atoms = self._all_atoms_coordinates(structure        = structure,
                                                hetero           = hetero,
                                                water            = water,
                                                geometric_center = by_residue)
        if by_residue:
            return all_atoms
        else:
            return numpy.sum(all_atoms, axis=0) / all_atoms[:, 1].size

    def clean(self, initatom = 1, initresidue = 1, nothetero = True, normalize = True):
        try:
            initatom    = int(initatom)
            initresidue = int(initresidue)
        except:
            raise AttributeError
        old_number = -99999

        if nothetero:
            working_array = self._structure
        else:
            working_array = self._all_residues
        for x in range(len(working_array)):
            # Atoms
            if initatom > 99999:
                initatom = 1
            working_array[x].renumerate_atoms(first_atom_number = initatom)
            initatom = working_array[x].last_atom_number + 1
            # Residues
            original_number = working_array[x].number
            if old_number != -99999:
                if old_number != original_number:
                    if original_number > old_number:
                        initresidue += original_number - old_number
                    else:
                        if working_array[x].follows(working_array[x-1]):
                            initresidue += 1
                        else:
                            initresidue += 2
                else:
                    initresidue += 1
            working_array[x].number  = initresidue
            working_array[x].version = ' '
            if normalize and isinstance(working_array[x], ResidueOfAminoAcid):
                if working_array[x].single_letter in aminoacids_main1:
                    working_array[x].normalize()
            old_number = original_number

    def renumerate_residues(self, init = 1, gaps = True):
        """
        Renumerates a chain starting in a given initial number

        @type  init: Integer
        @param init: Starting number of the renumerated chain.
                     1 by default.

        @type  gaps: Boolean
        @param gaps: Defines if the gaps in the chain have to be kept.
                     True by default.

        @raise AttributeError if init is not or cannot be casted to an integer
        """

        try:
            init = int(init)
        except:
            raise AttributeError

        for x in range(len(self)):
            original_number = self._all_residues[x].number
            self._all_residues[x].number = init
            if gaps:
                if x+1 < len(self):
                    init = init + (self._all_residues[x+1].number - original_number)
            else:
                init = init + 1

    def renumerate_atoms(self, init = 1):
        """
        Renumerates de atoms of a chain consecutively from a given value
        @type  init: Integer
        @param init: Starting atom number of the renumerated chain.
                     1 by default.

        @raise AttributeError if init is not or cannot be casted to an integer
        """

        try:
            init = int(init)
        except:
            raise AttributeError

        for x in range(len(self)):
            if init > 99999:
                init = 1
            self._all_residues[x].renumerate_atoms(first_atom_number = init)
            init = self._all_residues[x].last_atom_number + 1

    def duplicate(self, hetero = True, water = False, backbone = False):
        """
        Returns a {Chain} identical to the original but as a new object
        @rtype: {Chain}
        """
        new_chain = self.__class__(pdb = self._pdb, chain = self._chain)

        if not water:
            HOH_chain    = None
        else:
            HOH_chain    = copy.deepcopy(self._HOH)
        if not hetero:
            hetero_chain = None
        else:
            hetero_chain = copy.deepcopy(self._hetero)

        new_chain._add_residues(structure = copy.deepcopy(self._structure),
                                hetero    = hetero_chain,
                                HOH       = HOH_chain)
        if backbone:
            for residue in new_chain.aminoacids:
                residue.remove_side_chain()

        return new_chain

    def fuse(self, chain, lapl = False):
        """
        Merges another chain over self.
        @type chain: {Chain}
        """
        self._add_residues(structure = chain._structure,
                           hetero    = chain.heteroatoms,
                           HOH       = chain.waters)
        if lapl:
            self._last_appended_list = chain._last_appended_list
            self._residue_version = chain._residue_version

    #
    # PRIVATE METHODS
    #
    def _add_residues(self, structure = None, hetero = None, HOH = None):
        """
        Appends a list of {Residues} to the residue list

        @type residues: Array of {Residue}
        """
        if structure is not None:
            self._structure.extend(structure)
        if hetero is not None:
            self._hetero.extend(hetero)
        if HOH is not None:
            self._HOH.extend(HOH)

    def _new_Residue(self, number, version, Rtype, mode):
        """
        Creates a new {Residue}

        @type  number: Integer
        @param number: Residue number

        @type  version: Char
        @param version: Optional char used on pdbs to change count

        @type  type: String
        @param type: Residue type

        @type  mode: String
        @param mode: Residue mode: ATOM or HETATM

        @rtype: {Residue}
        """
        if self.dictitype.has_key(Rtype) and not self._term:
            return self.resitype(number = number, version = version, Rtype = Rtype, mode = mode)
        else:
            return Residue(number = number, version = version, Rtype = Rtype, mode = mode)

    def _new_Atom(self, number, name, x, y, z, extra):
        """
        Creates a new {Atom}

        @type  number: Integer
        @param number: Atom number

        @type  name: String
        @param name: Atom identifier

        @type  x: Decimal
        @param x: x coordinate

        @type  y: Decimal
        @param y: y coordinate

        @type  z: Decimal
        @param z: z coordinate

        @rtype: {Atom}
        """
        if isinstance(self._last_appended_list[-1], self.resitype):
            return self.atomtype(number = number, name = name, x = x, y = y, z = z, extra = extra)
        else:
            return Atom(number = number, name = name, x = x, y = y, z = z, extra = extra)

    def _chain_idx(self):
        """
        Returns indexes (number + version) of the residues separated by x on gaps
        @rtype: String
        """
        data = []
        pos = self._structure[0].number - 1
        for r in self._structure:
            if r.number == pos + 1:
                data.append(r.identifier)
            else:
                # this jump size is due to PDBs badly annotated
                if r.number - pos - 1 < 650:
                    for x in range(r.number - pos - 1):
                        data.append('X')
                data.append(r.identifier)
            pos = r.number

        return ";".join(data)

    #
    # BOOLEANS
    #
    @property
    def has_heteroatoms(self):
        """
        Returns True if there is any NON AMINOACID/NUCLEOTIDE heteroatom. False otherwise
        @rtype: Boolean
        """
        return len(self.heteroatoms) > 0

    @property
    def has_water(self):
        """
        Returns True if there is any HOH heteroatom. False otherwise
        @rtype: Boolean
        """
        return len(self.waters) > 0

    @property
    def is_empty(self):
        """
        Returns True if there isn't any {Residue} in the {Chain}. False otherwise
        @rtype: Boolean
        """
        return not self.has_water and not self.has_heteroatoms and len(self._structure) == 0

    #
    # IO
    #
    def read_PDB_line(self, line, keep_version = "A"):
        from ..parse_pdb import read_PDB_line
        read_PDB_line(self, line, keep_version)

    def PDB_format(self, terminal = True):
        """
        Strings a {Chain} in PDB format
        @rtype: String
        """
        lines = []
        for residue in self._all_residues:
            lines.extend(["%s%s %s%s%s%s%s%s%s%s" % (residue.mode.ljust(6),
                         str(atom.number).rjust(5),
                         str(atom.pretty_name),
                         residue.type.rjust(4),
                         str(self._chain).rjust(2),
                         str(residue.identifier).rjust(5),
                         ("%.3f" % float(atom.x)).rjust(11),
                         ("%.3f" % float(atom.y)).rjust(8),
                         ("%.3f" % float(atom.z)).rjust(8), atom._extra)
                for atom in residue.atoms])
        if terminal:
            lines.append("TER%s%s%s" % (residue.type.rjust(17),
                                        str(self._chain).rjust(2),
                                        str(residue.number).rjust(4)))

        return "\n".join(lines)

    def js_format(self, center=False):
        init = 'var {0}="'.format('pdb_' + self.globalID)
        return init + self.PDB_format().replace('\n', '\\n') + '"'

    #
    # OVERRIDE DEFAULT METHODS
    #
    def __len__(self):
        return len(self._all_residues)

    def __repr__(self):
        repre = []
        repre.append('<{0.__class__.__name__}: [{0.pdb}, {0.chain}]>'.format(self))
        for residue in self._all_residues:
            repre.append(repr(residue))
        return "\n".join(repre)
