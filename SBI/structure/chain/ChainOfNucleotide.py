"""
NucleotideChain

author: jbonet
date:   02/2013

@oliva's lab
"""

from ..atom       import Atom,    AtomOfNucleotide
from ..residue    import Residue, ResidueOfNucleotide
from .            import Chain
from SBI.data import nucleic2to1

class ChainOfNucleotide(Chain):
    """
    A {NucleotideChain} represents a collection of {AminoAcids}s
    """
    chaintype   = 'N'
    atomtype    = AtomOfNucleotide
    resitype    = ResidueOfNucleotide
    dictitype   = nucleic2to1
    
    def __init__(self, pdb, chain):
        """
        @type  pdb: String
        @param pdb: PDB identifier

        @type  chain: String
        @param chain: Chain identifier
        """
        super (ChainOfNucleotide, self).__init__(pdb = pdb, chain = chain)

    #
    # ATTRIBUTES
    #
    @property
    def nucleotides(self):
        """
        Returns the list of {AminoAcid}s
        @rtype: List of {AminoAcid}s
        """
        return self._structure

    @property
    def first_nucleotide(self):
        """
        Returns the first {AminoAcid}
        @rtype: {AminoAcid}
        """
        return self.nucleotides[0]

    @property
    def last_nucleotide(self):
        """
        Returns the last {AminoAcid}
        @rtype: {AminoAcid}
        """
        return self.nucleotides[-1]

    #
    # METHODS
    #
    def nucleotide_sequence(self):
        """
        Returns the exact sequence in the crystal
        @rtype: String
        """
        seq = ''
        for nc in self.nucleotides:
            seq += nc.single_letter
        return seq

    def gapped_nucleotide_sequence(self):
        """
        Returns the sequence of the crystal with 'X' where
        the crystal contains gaps.
        @rtype: String
        """
        seq = ''
        pos = self.nucleotides[0].number - 1
        for nc in self.nucleotides:
            if nc.number == pos + 1:
                seq += nc.single_letter
            else:
                for x in range(nc.number - pos - 1):
                    seq += 'x'
                seq += nc.single_letter
            pos = nc.number
        return seq

    def nucleotide_idx(self):
        """
        Returns indexes (number + version) of the residues separated by x on gaps
        @rtype: String
        """
        return self._chain_idx()

    def full_nucleotide_sequence(self):
        """
        Returns the sequence of the crystal with 'X' where
        the crystal contains gaps AND before the sequence if it does not start at 1.
        @rtype: String
        """
        if self.nucleotides[0].number <= 1: return self.gapped_nucleotide_sequence()

        seq = ''
        for x in range(self.nucleotides[0].number - 1): seq += 'x'
        return seq + self.gapped_nucleotide_sequence
