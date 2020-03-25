"""
PHInterface

author: jbonet
date:   03/2013

@oliva's lab
"""
from .         import Interface
from ..contact import ContactAH

class PHInterface(Interface):
    """
    {Interface} between two protein chains
    """
    def __init__(self, protein         = None,  heteroatom          = None, 
                       threshold_type  = "min", threshold_distance = 8,
                       protein_centers = None,  heteroatom_centers  = None):

        if protein_centers is None:
            protein_centers    = protein.geometric_center(structure = True,  hetero     = False,
                                                          water     = False, by_residue = True)
        if heteroatom_centers is None:
            heteroatom_centers = heteroatom.geometric_center(structure = False, hetero     = True,
                                                             water     = False, by_residue = True)

        super(PHInterface, self).__init__(chain1         = protein,         chain2             = heteroatom,
                                          threshold_type = threshold_type,  threshold_distance = threshold_distance,
                                          chain1_centres = protein_centers, chain2_centers     = heteroatom_centers)

    #
    # ATTRIBUTES
    #
    @property
    def protein(self):             return self._chain1
    @property
    def protein_id(self):          return self._chain1.globalID
    @property
    def protein_centers(self):     return self._chain1_centers

    @property
    def heteroatoms(self):         return self._chain2
    @property
    def heteroatoms_id(self):      return self._chain2.globalID
    @property
    def heteroatoms_centers(self): return self._chain1_centers

    @property
    def protein_positions(self):
        return self._list_positions(1)
    @property 
    def protein_view_interface(self):
        return self._view_interface_from(1)

    @property
    def heteroatoms_positions(self):
        return self._list_positions(2)
    @property 
    def heteroatoms_view_interface(self):
        return self._view_interface_from(2)

    #
    # OVERWRITE PARENT: PRIVATE FUNCTIONS
    #
    def _build(self):
        if len(self.protein.aminoacids) == 0 or len(self.heteroatoms.heteroatoms) == 0: return
        super(PHInterface, self)._build()

        for i in range(len(self._filtered[0])):
            new_contact = ContactAH(aminoacid          = self.protein.aminoacids[self._filtered[0][i]], 
                                    heteroatom         = self.heteroatoms.heteroatoms[self._filtered[1][i]], 
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
        data.append(super(PHInterface, self).toString(all_types))
        return "\n".join(data)