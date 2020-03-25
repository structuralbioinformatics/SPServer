"""
ProteinChain

author: jbonet
date:   02/2013

@oliva's lab
"""

from ..atom     import AtomOfAminoAcid
from ..residue  import ResidueOfAminoAcid
from .          import Chain
from ..protein  import SShelper
from ..protein  import Sequencer
from SBI.data   import aminoacids3to1

import numpy         as np


class ChainOfProtein(Chain):
    """
    A {ProteinChain} represents a collection of {AminoAcids}s
    """
    chaintype   = 'P'
    atomtype    = AtomOfAminoAcid
    resitype    = ResidueOfAminoAcid
    dictitype   = aminoacids3to1

    def __init__(self, pdb, chain):
        """
        @type  pdb: String
        @param pdb: PDB identifier

        @type  chain: String
        @param chain: Chain identifier
        """
        super(ChainOfProtein, self).__init__(pdb = pdb, chain = chain)
        self._ss_valid     = []
        self._archs        = []
        self._superarchs   = []
        self._torsionsCA   = '--'
        self._psiphi       = ''
        self._gaps         = []

    #
    # ATTRIBUTES
    #
    @property
    def aminoacids(self):
        """
        Returns the list of {AminoAcid}s
        @rtype: List of {AminoAcid}s
        """
        return self._structure

    @property
    def modified_aminoacids(self):
        for residue in self.aminoacids:
            if residue.mode == 'HETATM':
                yield residue

    @property
    def first_aminoacid(self):
        """
        Returns the first {AminoAcid}
        @rtype: {AminoAcid}
        """
        return self.aminoacids[0]

    @property
    def last_aminoacid(self):
        """
        Returns the last {AminoAcid}
        @rtype: {AminoAcid}
        """
        return self.aminoacids[-1]

    @property
    def torsionsCA(self):
        return self._torsionsCA

    @property
    def psiphi(self):
        data = ['', '', '']
        for i in range(0, len(self._psiphi), 3):
            data[0] += self._psiphi[i]
            data[1] += self._psiphi[i+1]
            data[2] += self._psiphi[i+2]
        return "\n".join(data)

    @property
    def secondary_structures(self):
        if not self.has_valid_secondary_structures:
            SShelper.locate_valid_secondarystructures(self)
        return self._ss_valid

    @property
    def archs(self):
        return self._archs

    @property
    def superarchs(self):
        return self._superarchs

    @property
    def protein_sequence(self):
        """
        Returns the exact sequence in the crystal
        @rtype: String
        """
        seq = ''
        for aa in self.aminoacids:
            seq += aa.single_letter
        return seq

    @property
    def gapped_protein_sequence(self):
        """
        Returns the sequence of the crystal with 'X' where
        the crystal contains gaps.
        @rtype: String
        """
        from ..protein import Sequencer
        return Sequencer._sequencer(self, 'seq')

    @property
    def full_protein_sequence(self):
        """
        Returns the sequence of the crystal with 'X' where
        the crystal contains gaps AND before the sequence
        if it does not start at 1.
        @rtype: String
        """
        if self.aminoacids[0].number <= 1:
            return self.gapped_protein_sequence

        seq = ''
        for x in range(self.aminoacids[0].number - 1):
            seq += 'x'
        return seq + self.gapped_protein_sequence

    @property
    def protein_idx(self):
        """
        Returns indexes (number + version) of the residues
        separated by x on gaps
        @rtype: String
        """
        return Sequencer._sequencer(self, 'idx')

    @property
    def gapped_protein_secondary_structure(self):
        return Sequencer._sequencer(self, 'str')

    @property
    def gapped_protein_exposition(self):
        return Sequencer._sequencer(self, 'exp')

    #
    # BOOLEANS
    #
    @property
    def is_only_ca(self):
        counter = 0
        for aa in self.aminoacids:
            if aa.is_only_ca:
                counter += 1

        return counter >= len(self.aminoacids) * 0.8

    @property
    def has_dssp(self):
        """
        Returns True if dssp calculations have been assigned
        to the {AminoAcids} of the {ProteinChain}
        """
        return self.first_aminoacid._dssp is not None

    @property
    def has_valid_secondary_structures(self):
        return len(self._ss_valid) > 0

    @property
    def has_torsionsCA(self):
        return self._torsionsCA != '--'

    @property
    def has_psiphi(self):
        return self._psiphi != ''

    @property
    def has_archs(self):
        return len(self._archs) > 0

    @property
    def has_superarchs(self):
        return len(self._superarchs) > 0

    #
    # METHODS
    #
    def get_aminoacid(self, identifier):
        """
        Returns the {AminoAcid} with the given position identifier
        @rtype: {AminoAcid}
        """
        if isinstance(identifier, int):
            identifier = str(identifier + ' ')
        for aa in self.aminoacids:
            if aa.identifier == identifier:
                return aa

    def revert_heteroAa(self):
        for residue in self.modified_aminoacids:
            residue.normalize()

    def get_backbone_coordinates(self, ca = False):
        all_coord = np.array(np.zeros(3))
        for residue in self.aminoacids:
            if ca is False:
                for atom in residue.backbone:
                    all_coord = np.vstack((all_coord, atom.coordinates))
            else:
                all_coord = np.vstack((all_coord, residue.ca.coordinates))

        all_coord = np.delete(all_coord, 0, 0)
        return all_coord

    def calculate_dssp(self, tmppdb = None, tmpdssp = None, cleanfiles = True):
        from ..protein import SShelper
        if not self.has_dssp:
            SShelper.calculate_dssp(self, tmppdb, tmpdssp, cleanfiles)

    def calculate_archs(self, limit_internal_ss=100,
                        limit_distance=False, allowed_gaps=0):

        if not self.has_dssp:
            self.calculate_dssp()
        if not self.has_valid_secondary_structures:
            self.secondary_structures
        if not self.has_torsionsCA:
            self.calculate_torsionsCA()
        if not self.has_psiphi:
            self.calculate_psiphi()

        SShelper.calculate_archs(self, limit_internal_ss,
                                 limit_distance, allowed_gaps)

    def calculate_torsionsCA(self, accuracy = 0):
        from ..protein import SShelper
        SShelper.calculate_torsionsCA(self, accuracy)

    def calculate_psiphi(self):
        from ..protein import SShelper
        SShelper.calculate_psiphi(self)

    #
    # REPRESENTATIONS
    #
    def archs2json(self, as_string = True):
        data = {'ID'    : self.globalID,
                'SEQ'   : self.gapped_protein_sequence,
                'STR'   : self.gapped_protein_secondary_structure,
                'ARCHS' : [x.json_format(False) for x in self.archs],
                'SUPER' : [x.json_format(False) for x in self.superarchs]}
        return repr(data) if as_string else data
