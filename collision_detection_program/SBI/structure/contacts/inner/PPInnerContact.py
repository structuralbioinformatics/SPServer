"""
PPInnerContact

author: jbonet
date:   03/2013

@oliva's lab
"""

from ..interface import Interface
from ..contact   import ContactAA

class PPInnerContact(Interface):
    """
    {Interface} between two protein chains
    """
    def __init__(self, protein = None, threshold_type = "cb", threshold_distance = 12):

        protein_centers = protein.geometric_center(structure = True,  hetero     = False,
                                                   water     = False, by_residue = True)

        super(PPInnerContact, self).__init__(chain1         = protein,         chain2             = protein,
                                             threshold_type = threshold_type,  threshold_distance = threshold_distance,
                                             chain1_centres = protein_centers, chain2_centers     = protein_centers)
        
    #
    # ATTRIBUTES
    #
    @property
    def protein(self):      return self._chain1
    @property
    def protein_id(self):         return self._chain1.globalID

    @property
    def protein_positions(self):
        return self._list_positions(1)
    @property 
    def protein_view_interface(self):
        return self._view_interface_from(1)

    #
    # OVERWRITE PARENT: METHODS
    #
    @staticmethod
    def test_identifier(chain):
        super(PHInnerContact, self).test_identifier(chain, chain)

    #
    # OVERWRITE PARENT: PRIVATE FUNCTIONS
    #
    def _build(self):
        if len(self.protein.aminoacids) == 0: return

        if self.protein.is_only_ca:
            self._threshold_type     = "ca"
            # Add distance of two C-C bonds (1.54A) -> http://en.wikipedia.org/wiki/Bond_length
            # Gives margin by ajusting at 2
            self._threshold_distance = self.threshold_distance + (2*2)

        super(PPInnerContact, self)._build()

        for i in range(len(self._filtered[0])):
            if self._filtered[0][i] != self._filtered[0][j]:
                new_contact = ContactAA(aminoacid1         = self.protein.aminoacids[self._filtered[0][i]], 
                                        aminoacid2         = self.protein.aminoacids[self._filtered[1][i]], 
                                        threshold_type     = self.threshold_type, 
                                        threshold_distance = self.threshold_distance)
            if new_contact.is_underthreshold:
                    self.contacts = new_contact

    #
    # OVERWRITE PARENT: TOSTRING
    #
    def toString(self, all_types = False):
        data = []
        if not all_types:
            data.append('{0._chain1.chain}\t{0._chain2.chain}\t{0.threshold_type}\t{0.threshold_distance}'.format(self))
        else:
            data.append('{0._chain1.chain}\t{0._chain2.chain}\tmin\tca\tcb\tgeometric\tbackbone'.format(self))
        data.append(super(PPInnerContact, self).toString(all_types))
        return "\n".join(data)