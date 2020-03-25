"""
Atom

author: jbonet
date:   02/2013

@oliva's lab
"""
import numpy as np


class Atom(object):
    """
    An {Atom} is simply a point in space defined by 3 coordinates
    """

    backbone_atoms = set()

    def __init__(self, number, name, x, y, z, extra):
        """
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
        """
        self._number      = int(number)
        self._name        = name
        self._coordinates = np.array([float(x), float(y), float(z)])
        self._extra       = extra

    #
    # ATTRIBUTES
    #
    @property
    def number(self):
        """
        Atom number
        @rtype: Integer
        """
        return int(self._number)

    @number.setter
    def number(self, value):
        self._number = int(value)

    @property
    def name(self):
        """
        Atom name
        @rtype: String
        """
        return self._name

    @property
    def pretty_name(self):
        """
        Returns a version of the name adapted to print in PDB
        The String has a fixed length for that matter
        @rtype: String{4}
        """
        if len(self._name) <= 3:
            return " "  + self._name.ljust(3)
        else:
            return self._name.ljust(4)

    @property
    def x(self):
        """
        x coordinate
        @rtype: Decimal
        """
        return float(self._coordinates[0])

    @property
    def y(self):
        """
        y coordinate
        @rtype: Decimal
        """
        return float(self._coordinates[1])

    @property
    def z(self):
        """
        z coordinate
        @rtype: Decimal
        """
        return float(self._coordinates[2])

    @property
    def coordinates(self):
        """
        @rtype: numpy.array
        """
        return self._coordinates

    @property
    def is_backbone(self):
        """
        Checks if the Atom is part of the backbone of the residue
        @rtype: Boolean
        """
        return self._name in self.backbone_atoms

    @property
    def is_sidechain(self):
        """
        Checks if the Atom is part of the side chain of the residue
        @rtype: Boolean
        """
        return not self.is_backbone if len(self.backbone_atoms) > 0 else False

    #
    # METHODS
    #
    def rotate(self, matrix = None):
        """
        Rotates an atom according to a given matrix
        """
        if matrix is None:
            matrix = np.identity(3)
        self._coordinates = np.dot(matrix, self._coordinates)

    def translate(self, vector = None):
        """
        Translates an atom according to a translational vector
        """
        if vector is None:
            vector = np.zeros(3)
        self._coordinates = self._coordinates + vector

    def distance(self, atom = None):
        """
        Euclidean distance between two atoms
        """
        if atom is None:
            return np.linalg.norm(self._coordinates - np.zeros(3))
        else:
            return np.linalg.norm(self._coordinates - atom.coordinates)

    #
    # OVERRIDE DEFAULT METHODS
    #
    def __repr__(self):
        return '<{0.__class__.__name__}: [{0.name}, {0.number}]:({0.x:.3f}, {0.y:.3f}, {0.z:.3f})>'.format(self)
