"""
ContactAH

author: jbonet
date:   02/2014

@oliva's lab
"""

from . import Contact

class ContactAH(Contact):

    available_distance_types = set(["min", "geometric"])
    atomic_types             = set(["min"])
    description              = 'AminoHetero'

    def __init__(self, aminoacid, heteroatom, threshold_type = "min", threshold_distance = 6):

        self._check_type(requested_type = threshold_type)

        self._threshold_type     = threshold_type
        self._threshold_distance = threshold_distance

        super(ContactAH, self).__init__(residue1 = aminoacid, residue2 = heteroatom)

    #
    # ATTRIBUTES
    #
    @property
    def aminoacid(self):          return self._residue1

    @property
    def heteroatom(self):         return self._residue2

    @property
    def min_distance(self):       return self._distance["min"][2]
    @property
    def min_atoms(self):          return self._distance["min"][0:2]

    @property
    def geometric_distance(self):   return self._distance["geometric"][2]

    #
    # OVERWRITE PARENT FUNCTION
    #
    def _build(self):
        for dist_type in self.available_distance_types:
            self._distance.setdefault(dist_type, None)

        self._distance[self._threshold_type] = self.aminoacid.distance(self.heteroatom, dist_type = self._threshold_type)
        if float(self._distance[self._threshold_type][2]) <= self._threshold_distance and \
           float(self._distance[self._threshold_type][2]) >= 0:
            self._underthreshold = True
            for dist_type in self._distance:
                if dist_type != self._threshold_type:
                    self._distance[dist_type] = self.aminoacid.distance(self.heteroatom, dist_type = dist_type)

    #
    # TOSTRING
    #
    def toString(self, requested_type = "cb"):
        if requested_type != "all":
            self._check_type(requested_type = requested_type)
            ddata = self._distance[requested_type]
            if requested_type in self.atomic_types:
                return '{0.aminoacid.type}:{0.aminoacid.number}:{1[0].name}\t{0.heteroatom.type}:{0.heteroatom.number}:{1[1].name}\t{1[2]:.3f}'.format(self, ddata)
            else:
                return '{0.aminoacid.type}:{0.aminoacid.number}\t{0.heteroatom.type}:{0.heteroatom.number}\t{1[2]:.3f}'.format(self, ddata)
        else:
            return '{0.aminoacid.type}:{0.aminoacid.number}\t{0.heteroatom.type}:{0.heteroatom.number}\t{0.min_distance:.3f}\t{0.geometric_distance:.3f}'.format(self)
        