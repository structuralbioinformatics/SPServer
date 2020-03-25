"""
NucleotideAtom

author: jbonet
date:   02/2013

@oliva's lab
"""

from . import Atom


class AtomOfNucleotide(Atom):
    """
    An {AtomOfNucleotide} is simply a point in space defined by 3 coordinates
    WITH specific functions for atoms in amino acids
    """

    backbone_atoms = set(["P", "OP1", "OP2", "OP3", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'",
                                                    "O5*", "C5*", "C4*", "O4*", "C3*", "O3*", "C2*", "O2*", "C1*"])
    #
    # BOOLEANS
    #

    @property
    def is_Phosphate(self):     return self._name == "P"

    @property
    def is_PhosphoOxygen(self): return self._name in set(["OP1", "OP2", "OP3"])
