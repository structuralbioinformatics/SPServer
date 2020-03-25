"""
Interface

author: jbonet
date:   02/2014

@oliva's lab
"""
import copy, hashlib
import scipy.spatial as sp
import numpy         as np

class Interface(object):
    """docstring for Interface"""
    def __init__(self, chain1, chain2, threshold_type, threshold_distance, chain1_centres = None, chain2_centers = None):
        
        self._chain1             = chain1
        c1                       = chain1.globalID
        self._chain1_centers     = chain1_centres
        self._chain2             = chain2
        c2                       = chain2.globalID
        self._chain2_centers     = chain2_centers
        self._md5                = tuple(sorted([hashlib.sha224(c1 + c2).hexdigest(),
                                                 hashlib.sha224(c2 + c1).hexdigest()])).__hash__()
        self._threshold_type     = threshold_type
        self._threshold_distance = threshold_distance

        self._contacts           = set()

        self._filterdistance     = threshold_distance if threshold_type is "geometric" else threshold_distance + 15
        self._filtered           = None

        self._build()

    #
    # ATTRIBUTES
    #
    @property 
    def identifier(self):         return self._md5

    @property
    def threshold_type(self):     return self._threshold_type

    @property
    def threshold_distance(self): return self._threshold_distance

    @property
    def contacts(self):           return self._contacts
    @contacts.setter
    def contacts(self, contact):  self._contacts.add(contact)


    @property
    def is_empty(self):
        return len(self._contacts) == 0
    #
    # METHODS
    #
    @staticmethod
    def test_identifier(chain1, chain2):
        c1 = chain1.globalID
        c2 = chain2.globalID
        return tuple(sorted([hashlib.sha224(c1 + c2).hexdigest(),
                             hashlib.sha224(c2 + c1).hexdigest()])).__hash__()

    def reverse(self):
        (self._chain1, self._chain2) = (self._chain2, self._chain1)
        map(lambda x: x.reverse(), self._contacts)

    #
    # PRIVATE METHODS
    #
    def _build(self):
        if self._chain1_centers is not None and self._chain2_centers is not None:
            distance_matrix = sp.distance.cdist(self._chain1_centers, self._chain2_centers)
            self._filtered  = np.where(distance_matrix <= self._filterdistance)

    def _list_positions(self, chain):
        positions = set()
        for c in self.contacts:
            if chain   == 1: r = c._residue1
            elif chain == 2: r = c._residue2
            positions.add(r)
        return positions

    def _view_interface_from(self, chain):
        positions = {}
        for c in self.contacts:
            if chain   == 1: r = c._residue1
            elif chain == 2: r = c._residue2
            positions.setdefault(r, []).append(c)
        return positions

    #
    # OVERRIDE DEFAULT METHODS
    #
    def __len__(self): return len(self._contacts)
    def __add__(self, other):
        if self.threshold_type != other.threshold_type:
            raise AttributeError("Different threshold type interfaces can not be added\n")
        if not isinstance(other, self.__class__):
            raise AttributeError("Different types of interfaces can not be joined\n")
        if self._chain1.globalID == other._chain1.globalID  and self._chain2.globalID  == other._chain2.globalID:
            self._contacts.update(other._contacts)
            self._threshold_distance = self.threshold_distance if self.threshold_distance > other.threshold_distance else other.threshold_distance
        elif self._chain1.globalID == other._chain2.globalID  and self._chain2.globalID  == other._chain1.globalID:
            t = copy.deepcopy(other)
            t.reverse()
            self.__add__(t)
        else:
            raise AttributeError("Interfaces from different chains can not be added\n")
    
    #
    # TOSTRING
    #
    def toString(self, all_types = False):
        data = []
        # data.append('{0._chain1.chain}\t{0._chain2.chain}\t{0.threshold_type}\t{0.threshold_distance}'.format(self))
        for contact in self.contacts:
            if not all_types:
                data.append(contact.toString(self.threshold_type))
            else:
                data.append(contact.toString("all"))
        return "\n".join(data)
