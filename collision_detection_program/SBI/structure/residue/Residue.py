"""
Residue

author: jbonet
date:   02/2013

@oliva's lab
"""

import numpy         as np
import scipy.spatial as sp

from SBI import SBIglobals


class Residue(object):
    """
    A {Residue} collects a series of {Atom}
    """

    def __init__(self, number, version, Rtype, mode):
        """
        @type  number: Integer
        @param number: Residue number

        @type  version: Char
        @param version: Optional char used on pdbs to change count

        @type  type: String
        @param type: Residue type

        @type  mode: String
        @param mode: Residue mode: ATOM or HETATM
        """
        self._number    = int(number)
        self._version   = version
        self._type      = Rtype
        self._mode      = mode

        # self._atoms               = []
        self._backbone_coordinates  = None
        self._backbone_atoms        = []
        self._sidechain_coordinates = None
        self._sidechain_atoms       = []

    #
    # ATTRIBUTES
    #
    @property
    def number(self):
        """
        Residue number
        @rtype: Integer
        """
        return int(self._number)

    @number.setter
    def number(self, value):
        """
        @type value: Integer
        """
        self._number = int(value)

    @property
    def version(self):
        """
        Residue version
        @rtype: CHAR
        """
        return self._version

    @version.setter
    def version(self, value):
        self._version = value

    @property
    def identifier(self):
        """
        Combo number + version
        @rtype: String
        """
        return str(self.number) + self.version

    @property
    def type(self):
        """
        Residue type
        @rtype: String
        """
        return self._type

    @property
    def mode(self):
        """
        Residue mode: ATOM or HETATM
        @rtype: String
        """
        return self._mode

    @property
    def atoms(self):
        """
        List of {Atom} objects of the {Residue}
        @rtype: List of {Atom}
        """
        # return self._atoms
        return self._backbone_atoms + self._sidechain_atoms

    @property
    def backbone_atoms(self):
        return self._backbone_atoms

    @property
    def sidechain_atoms(self):
        return self._sidechain_atoms

    @property
    def all_coordinates(self):
        """
        Coordinates of all the {Atom}s in the {Residue}
        @rtype: numpy.matrix
        """
        if self._backbone_coordinates is not None and self._sidechain_coordinates is not None:
            return np.vstack((self._backbone_coordinates, self._sidechain_coordinates))
        else:
            if self._backbone_coordinates is not None:
                return self._backbone_coordinates
            elif self._sidechain_coordinates is not None:
                return self._sidechain_coordinates

    @property
    def first_atom_number(self):
        """
        Returns the number of the first {Atom} of the {Residue}
        @rtype: Integer
        """
        # return self._atoms[0].number
        return self.atoms[0].number

    @property
    def last_atom_number(self):
        """
        Returns the number of the last {Atom} of the {Residue}
        @rtype: Integer
        """
        # return self._atoms[-1].number
        return self.atoms[-1].number

    @property
    def geometric_center(self):
        """
        Returns the geometric central position of the {Residue}
        @rtype numpy.array
        """
        return sum(self.all_coordinates)/self.all_coordinates[:, 1].size

    #
    # ADVANCED GETTERS & SETTERS
    #
    def add_atom(self, atom):
        """
        Adds a new {Atom} to the {Residue}
        This includes filling the all_coordinates attribute

        @type  atom: {Atom}
        @param atom: New {Atom} added to the {Residue}
        """
        # self._atoms.append(atom)
        if atom.is_backbone:
            self._backbone_atoms.append(atom)
        else:
            self._sidechain_atoms.append(atom)
        self._add_to_matrix(atom)

    #
    # METHODS
    #
    def normalize(self):
        raise NotImplementedError

    def renumerate_atoms(self, first_atom_number):
        """
        Given the number for the first atom, renumerates the rest accordingly
        @type  value: Integer
        """
        for atom in self.atoms:
            atom.number = first_atom_number
            first_atom_number += 1

    def rotate(self, matrix = None):
        """
        Rotates a {Residue} according to a given matrix

        @type matrix: numpy.matrix
        """
        if matrix is None:
            matrix = np.identity(3, float)
        self._backbone_coordinates  = None
        self._sidechain_coordinates = None
        for atom in self.atoms:
            atom.rotate(matrix = matrix)
            self._add_to_matrix(atom)

    def translate(self, vector):
        """
        Translates a {Residue} according to a translational vector

        @type vector: numpy.array
        """
        if vector is None:
            vector = np.zeros(3, float)
        self._backbone_coordinates  = None
        self._sidechain_coordinates = None
        for atom in self.atoms:
            atom.translate(vector = vector)
            self._add_to_matrix(atom)

    def reposition(self, matrix = None, vector = None):
        """
        Rotates and translates the {Residue} according to a matrix and translational vector

        @type matrix: numpy.matrix

        @type vector: numpy.array
        """
        if matrix is None:
            matrix = np.identity(3, float)
        if vector is None:
            vector = np.zeros(3, float)
        SBIglobals.alert('deepdebug', self, 'Reposition residue {0.type}:{0.number}'.format(self))

        self._backbone_coordinates  = None
        self._sidechain_coordinates = None
        for atom in self.atoms:
            SBIglobals.alert('deepdebug', self, 'Atom {0.name} {0.is_backbone}'.format(atom))
            atom.rotate(matrix = matrix)
            atom.translate(vector = vector)
            self._add_to_matrix(atom)

    def _add_to_matrix(self, atom):
        def add_to_matrix(matrix, atom):
            if matrix is None:
                matrix = np.array(np.zeros(3))
                matrix = np.vstack((matrix, atom.coordinates))
                matrix = np.delete(matrix, 0, 0)
            else:
                matrix = np.vstack((matrix, atom.coordinates))
            return matrix
        if atom.is_backbone:
            self._backbone_coordinates  = add_to_matrix(self._backbone_coordinates, atom)
        else:
            self._sidechain_coordinates = add_to_matrix(self._sidechain_coordinates, atom)

    def distance(self, residue, dist_type = 'min'):
        """
        Calculates the distance between two residues

        @type  residue: {Residue}
        @param residue: Residue to which we want to find the distance

        @type  dist_type: String
        @param dist_type: Type of distance to evaluate
        @limit dist_type: Accepted: distance.types (min, geometric)

        @rtype: (AT1, AT2, DISTANCE)
        """
        if dist_type not in set(['min', 'geometric', 'backbone']): dist_type = 'min'

        if dist_type.strip().lower() == 'min':
            return self._distance_min(residue = residue)
        if dist_type.strip().lower() == 'geometric':
            return self._distance_geometric_center(residue = residue)
        if dist_type.strip().lower() == 'backbone':
            return self._distance_min_backbone(residue = residue)

    def follows(self, residue):
        number0  = self.number
        version0 = self.version
        number1  = residue.number
        version1 = residue.version
        if number0 == number1 + 1: return True
        if number0 == number1:
            if version0 == ' ': version0 = '@'
            if version1 == ' ': version1 = '@'
            if ord(version0) == ord(version1) + 1: return True
        return False

    def is_followed(self, residue):
        number0  = self.number
        version0 = self.version
        number1  = residue.number
        version1 = residue.version
        if number0 == number1 - 1: return True
        if number0 == number1:
            if version0 == ' ': version0 = '@'
            if version1 == ' ': version1 = '@'
            if ord(version0) == ord(version1) - 1: return True
        return False

    def identifier_distance(self, residue):
        number0  = self.number
        version0 = self.version
        number1  = residue.number
        version1 = residue.version

        if number0 != number1:
            return number1-number0
        else:
            return abs(ord(version1)-ord(version0))

    def _distance_min(self, residue):
        """
        Calculates the minimum distance between two residues

        @type  residue: {Residue}
        @param residue: Residue to which we want to find the distance

        @rtype:
        """
        distances = sp.distance.cdist(self.all_coordinates, residue.all_coordinates)
        index     = np.unravel_index(distances.argmin(), distances.shape)
        return (self.atoms[index[0]], residue.atoms[index[1]], distances.min())

    def _distance_min_backbone(self, residue):
        if self._backbone_coordinates is None or residue._backbone_coordinates is None:
            return (None, None, -1)
        else:
            distances = sp.distance.cdist(self._backbone_coordinates, residue._backbone_coordinates)
            index     = np.unravel_index(distances.argmin(), distances.shape)
            return (self._backbone_atoms[index[0]], residue._backbone_atoms[index[1]], distances.min())

    def _distance_geometric_center(self, residue):
        """
        Calculates the distance between the geometric center of two residues

        @type  residue: {Residue}
        @param residue: Residue to which we want to find the distance

        @rtype:
        """
        return (None, None, sp.distance.euclidean(self.geometric_center, residue.geometric_center))

    #
    # OVERRIDE DEFAULT OBJECT METHODS
    #
    def __len__(self):
        """
        number of atoms in the Residue
        @rtype: Integer
        """
        return len(self.atoms)

    def __hash__(self):
        return '<{0.__class__.__name__}: [{0.type}, {0.number}, {0.version}]>'.format(self).__hash__()

    def __eq__(self, other):
        return (self.type, self.number, self.version) == (other.type, other.number, other.version)

    def __sub__(self, other):
        return np.subtract(self.geometric_center, other.geometric_center)

    def __repr__(self):
        repre = []
        repre.append('<{0.__class__.__name__}: [{0.type}, {0.number}, {0.version}]>'.format(self))
        for atom in self.atoms:
            repre.append("\t" + repr(atom))
        return "\n".join(repre)
