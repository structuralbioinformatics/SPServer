"""
AminoAtom

author: jbonet
date:   02/2013

@oliva's lab
"""

from . import Atom


class AtomOfAminoAcid(Atom):
    """
    An {AtomOfAminoAcid} is simply a point in space defined by 3 coordinates
    WITH specific functions for atoms in amino acids
    """

    backbone_atoms = set(['N', 'CA', 'C'])

    #
    # BOOLEANS
    #
    @property
    def is_Calpha(self): return self._name == "CA"

    @property
    def is_Cbeta(self):  return self._name == "CB"

    @property
    def is_N(self):      return self._name == "N"

    @property
    def is_C(self):      return self._name == "C"

    @property
    def is_O(self):      return self._name == "O"
