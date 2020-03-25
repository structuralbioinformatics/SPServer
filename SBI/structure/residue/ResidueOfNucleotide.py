"""
Nucleotide

author: jbonet
date:   02/2013

@oliva's lab
"""

from .          import Residue
from SBI.data   import nucleic2to1
import numpy    as np

class ResidueOfNucleotide(Residue):
    """
    A {Nucleotide} collects a series of {NucleotideAtom}s
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
        super(ResidueOfNucleotide, self).__init__(number = number, version = version, Rtype = Rtype, mode = mode)
        

    #
    # ATTRIBUTES
    #
    @property
    def single_letter(self):
        """
        Returns the AminoAcid identifier as a single letter code
        @rtype: String
        """
        return nucleic2to1[self.type]
