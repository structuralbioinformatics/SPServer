"""
PHInnerContact

author: jbonet
date:   03/2013

@oliva's lab
"""
from ..interface import Interface
from ..contact   import ContactAH

class PHInnerContact(Interface):
    """
    {Interface} between two protein chains
    """
    def __init__(self, protein = None, threshold_type = "min", threshold_distance = 8):

        protein_centers = protein.geometric_center(structure = True,  hetero     = False,
                                                   water     = False, by_residue = True)
        hetero_centers  = protein.geometric_center(structure = False, hetero     = True,
                                                   water     = False, by_residue = True)

        super(PHInnerContact, self).__init__(chain1         = protein,         chain2             = protein,
                                             threshold_type = threshold_type,  threshold_distance = threshold_distance,
                                             chain1_centres = protein_centers, chain2_centers     = hetero_centers)
        
    #
    # ATTRIBUTES
    #
    @property
    def protein(self):       return self._chain1
    @property
    def protein_id(self):    return self._chain1.globalID

    @property
    def protein_positions(self):
        return self._list_positions(1)
    @property 
    def protein_view_innercontact(self):
        return self._view_interface_from(1)

    @property
    def heteroatoms_positions(self):
        return self._list_positions(2)
    @property 
    def heteroatoms_view_innercontact(self):
        return self._view_interface_from(2)

    #
    # OVERWRITE PARENT: METHODS
    #
    @staticmethod
    def test_identifier(chain):
        return Interface.test_identifier(chain, chain)

    #
    # OVERWRITE PARENT: PRIVATE FUNCTIONS
    #
    def _build(self):
        if len(self.protein.aminoacids) == 0 or len(self.protein.heteroatoms) == 0: return

        super(PHInnerContact, self)._build()

        for i in range(len(self._filtered[0])):
            new_contact = ContactAH(aminoacid          = self.protein.aminoacids[self._filtered[0][i]], 
                                    heteroatom         = self.protein.heteroatoms[self._filtered[1][i]], 
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
            data.append('{0._chain1.chain}\t{0._chain2.chain}\tmin\tgeometric'.format(self))
        data.append(super(PHInnerContact, self).toString(all_types))
        return "\n".join(data)