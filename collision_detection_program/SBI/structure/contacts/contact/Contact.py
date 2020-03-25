"""
Contact

author: jbonet
date:   03/2013

@oliva's lab
"""
import hashlib


class Contact(object):

    available_distance_types = set()
    atomic_types             = set()
    description              = ''

    def __init__(self, residue1, residue2):

        self._residue1  = residue1
        r1              = residue1.type + residue1.identifier
        self._residue2  = residue2
        r2              = residue2.type + residue2.identifier
        self._md5       = tuple(sorted([hashlib.sha224(r1 + r2).hexdigest(),
                                        hashlib.sha224(r2 + r1).hexdigest()]))
        self._distance  = {}

        self._underthreshold = False

        self._build()

    #
    # ATTRIBUTES
    #
    @property
    def md5s(self):
        return self._md5

    @property
    def description(self):
        return self.description

    #
    # BOOLEAN
    #
    @property
    def is_underthreshold(self):
        return self._underthreshold

    #
    # METHODS
    #
    def reverse(self):
        (self._residue1, self._residue2) = (self._residue2, self._residue1)
        for dist_type in self.available_distance_types:
            self._distance[dist_type] = (self._distance[dist_type][1], self._distance[dist_type][0], self._distance[dist_type][2])

    #
    # PRIVATE METHODS
    #
    def _build(self):
        raise NotImplementedError

    def _check_type(self, requested_type):
        if requested_type not in self.available_distance_types:
            raise AttributeError("Available distance types are {0}.\n".format(repr(self.available_distance_types)))

    #
    # OVERRIDE DEFAULT METHODS
    #
    def __eq__(self, other): return self.md5s == other.md5s
    def __ne__(self, other): return not self.__eq__(other)
    def __hash__(self):      return self.md5s.__hash__()

# class Contact(object):
#     """
#     Stores the basic information of a contact between two residues
#     """
#     def __init__(self, residue1, residue2, atom1 = None, atom2 = None, distance = 0, altdistance1 = -1, altdistance2 = -1):

#         self._residue1  = residue1
#         self._atom1     = atom1
#         self._residue2  = residue2
#         self._atom2     = atom2
#         self._distance  = [float(distance), float(altdistance1), float(altdistance2)]

#     @property
#     def residue1(self):
#         """
#         Data referring to residue1 of the contact
#         @rtype: {Residue}
#         """
#         return self._residue1

#     @property
#     def atom1(self):
#         """
#         Data referring to the atom of residue1 with the assigned interaction
#         @rtype: {Atom}
#         """
#         return self._atom1

#     @property
#     def residue2(self):
#         """
#         Data referring to residue2 of the contact
#         @rtype: {Residue}
#         """
#         return self._residue2

#     @property
#     def atom2(self):
#         """
#         Data referring to the atom of residue2 with the assigned interaction
#         @rtype: {Atom}
#         """
#         return self._atom2

#     @property
#     def distance(self):
#         """
#         Distance assigned to the contact
#         @rtype: Float
#         """
#         return self._distance[0]

#     #
#     # BOOLEAN
#     #
#     def compare_residues(self, contact):
#         """
#         Compares if the residue positions are the same

#         @type contact: {Contact}
#         @rtype: Boolean
#         """
#         if (self.residue1.number == contact.residue1.number) and \
#            (self.residue2.number == contact.residue2.number):
#             return True
#         return False

#     #
#     # FULL CONTACT DATA
#     #
#     def print_aldist(self):
#         if self._distance[1] == -1:
#             return str(self)
#         elif self._distance[2] == -1:
#             return '{0.residue1.type}:{0.residue1.number}\t{0.residue2.type}:{0.residue2.number}\t{0.distance}\t{0._distance[1]}'.format(self)
#         else:
#             return '{0.residue1.type}:{0.residue1.number}\t{0.residue2.type}:{0.residue2.number}\t{0.distance}\t{0._distance[1]}\t{0._distance[2]}'.format(self)

#     #
#     # OVERRIDE DEFAULT METHODS
#     #
#     def __repr__(self):
#         return '<{0.__class__.__name__}: [{0.residue1.type}, {0.residue1.number}, {0.atom1.name}] [{0.residue2.type}, {0.residue2.number}, {0.atom2.name}] {0.distance}>'.format(self)

#     def __str__(self):
#         return '{0.residue1.type}:{0.residue1.number}:{0.atom1.name}\t{0.residue2.type}:{0.residue2.number}:{0.atom2.name}\t{0.distance}'.format(self)
