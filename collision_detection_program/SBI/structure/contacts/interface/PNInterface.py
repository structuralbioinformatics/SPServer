"""
PNInterface

author: jbonet
date:   03/2013

@oliva's lab
"""
from .         import Interface
from ..contact import ContactAN
from SBI       import SBIglobals

class PNInterface(Interface):
    """
    {Interface} between two protein chains
    """
    def __init__(self, protein         = None,  nucleotide         = None, 
                       threshold_type  = "min", threshold_distance = 8,
                       protein_centers = None,  nucleotide_centers = None):

        if protein_centers is None:
            protein_centers    = protein.geometric_center(structure = True,  hetero     = False,
                                                          water     = False, by_residue = True)
        if nucleotide_centers is None:
            nucleotide_centers = nucleotide.geometric_center(structure = True,  hetero     = False,
                                                             water     = False, by_residue = True)

        super(PNInterface, self).__init__(chain1         = protein,         chain2             = nucleotide,
                                          threshold_type = threshold_type,  threshold_distance = threshold_distance,
                                          chain1_centres = protein_centers, chain2_centers     = nucleotide_centers)


    #
    # ATTRIBUTES
    #
    @property
    def protein(self):            return self._chain1
    @property
    def protein_id(self):         return self._chain1.globalID
    @property
    def protein_centers(self):    return self._chain1_centers

    @property
    def nucleotide(self):         return self._chain2
    @property
    def nucleotide_id(self):      return self._chain2.globalID
    @property
    def nucleotide_centers(self): return self._chain1_centers

    @property
    def protein_positions(self):
        return self._list_positions(1)
    @property 
    def protein_view_interface(self):
        return self._view_interface_from(1)

    @property
    def nucleotide_positions(self):
        return self._list_positions(2)
    @property 
    def nucleotide_view_interface(self):
        return self._view_interface_from(2)

    #
    # OVERWRITE PARENT: PRIVATE FUNCTIONS
    #
    def _build(self):
        if len(self.protein.aminoacids) == 0 or len(self.nucleotide.nucleotides) == 0: return
        super(PNInterface, self)._build()

        for i in range(len(self._filtered[0])):
            SBIglobals.alert('deepdebug', self, 'Analyze AN Contact for {0.type}:{0.number} - {1.type}:{1.number}'.format(self.protein.aminoacids[self._filtered[0][i]],self.nucleotide.nucleotides[self._filtered[1][i]]))
            new_contact = ContactAN(aminoacid          = self.protein.aminoacids[self._filtered[0][i]], 
                                    nucleotide         = self.nucleotide.nucleotides[self._filtered[1][i]], 
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
            data.append('{0._chain1.chain}\t{0._chain2.chain}\tmin\tgeometric\tcbbackbone'.format(self))
        data.append(super(PNInterface, self).toString(all_types))
        return "\n".join(data)


# from ContactList import ContactList

# class PNInterface(ContactList):
#     """
#     {ContactList} between two chains
#     """
#     def __init__(self, protein_chain = None, nucleotide_chain = None, c_type = "min", max_distance = 12, all_distances = False):
#         self._alldistances = all_distances
#         super(PNInterface, self).__init__(chain = protein_chain, secondary_chain = nucleotide_chain,
#                                           c_type = c_type, max_distance = max_distance)
        

#     #
#     # ATTRIBUTES
#     #
#     @property
#     def protein_chain(self):
#         return self.chain
#     @property
#     def nucleotide_chain(self):
#         """
#         Structure of the interactor chain
#         @rtype: {CHAIN}
#         """
#         return self._schain

#     @property
#     def nucleotide_chain_id(self):
#         """
#         Identifier of interactor chain
#         @rtype: String
#         """
#         return self._schain.globalID()

#     @property
#     def protein_chain_positions(self):
#         """
#         Positions in chain1 related to the interaction
#         [OVERWRITTED FROM PARENT]
#         @rtype: set
#         """
#         positions = set()
#         for contact in self._contacts:
#             positions.add(contact.residue1.number)
#         return positions

#     @property
#     def nucleotide_chain_positions(self):
#         """
#         Positions in chain1 related to the interaction
#         @rtype: set
#         """
#         positions = set()
#         for contact in self._contacts:
#             positions.add(contact.residue2.number)
#         return positions

#     #
#     # BOOLEANS
#     #
#     def fulfills_interface(self, interface, strict = False):
#         """
#         Checks if a list of contacts appear in the interface

#         @type  interface: {Interface}
#         @param contact_list: List of {Contact}s to check

#         @type    strict: Boolean
#         @param   strict: Strict checking... all the given contacts must exist
#         @default strict: False

#         @rtype: Boolean
#         """
#         total = 0

#         for c1 in self.contacts:
#             for c2 in interface.contacts:
#                 if c1 == c2:
#                     total += 1

#         if total == len(interface):  return True
#         if not strict and total > 0: return True
#         return False


#     def fulfills_exclude_residues_chain(self, residue_number_list, strict = False):
#         """
#         Checks if a list of residues (number) DOES NOT appear in the chain

#         @type  residue_number_list: List of Integers
#         @param residue_number_list: List of residue numbers to check

#         @type    strict: Boolean
#         @param   strict: Strict checking... EVERY REQUESTED RESIDUE HAS NOT TO BE IN THE INTERFACE
#         @default strict: False

#         @rtype: Boolean
#         """
#         return not self._fulfills_residues(residue_number_list = residue_number_list, chain = True, strict = strict)

#     def fulfills_residues_chain(self, residue_number_list, strict = False):
#         """
#         Checks if a list of residues (number) appear in the chain

#         @type  residue_number_list: List of Integers
#         @param residue_number_list: List of residue numbers to check

#         @type    strict: Boolean
#         @param   strict: Strict checking... EVERY REQUESTED RESIDUE HAS TO BE IN THE INTERFACE
#         @default strict: False

#         @rtype: Boolean
#         """
#         return self._fulfills_residues(residue_number_list = residue_number_list, chain = True, strict = strict)

#     def fulfills_exclude_residues_interactor(self, residue_number_list, strict = False):
#         """
#         Checks if a list of residues (number) DOES NOT appear in the interactor

#         @type  residue_number_list: List of Integers
#         @param residue_number_list: List of residue numbers to check

#         @type    strict: Boolean
#         @param   strict: Strict checking... EVERY REQUESTED RESIDUE HAS NOT TO BE IN THE INTERFACE
#         @default strict: False

#         @rtype: Boolean
#         """
#         return not self._fulfills_residues(residue_number_list = residue_number_list, interactor = True, strict = strict)

#     def fulfills_residues_interactor(self, residue_number_list, strict = False):
#         """
#         Checks if a list of residues (number) appear in the interactor

#         @type  residue_number_list: List of Integers
#         @param residue_number_list: List of residue numbers to check

#         @type    strict: Boolean
#         @param   strict: Strict checking... EVERY REQUESTED RESIDUE HAS TO BE IN THE INTERFACE
#         @default strict: False

#         @rtype: Boolean
#         """
#         return self._fulfills_residues(residue_number_list = residue_number_list, interactor = True, strict = strict)

#     def _fulfills_residues(self, residue_number_list, chain = False, interactor = False, strict = False):
#         """
#         Checks if a list of residues (number) appear in a given chain

#         @type  residue_number_list: List of Integers
#         @param residue_number_list: List of residue numbers to check

#         @type    chain: Boolean
#         @param   chain: Check on chain?
#         @default chain: False

#         @type    interactor: Boolean
#         @param   interactor: Check on interactor?
#         @default interactor: False

#         @type    strict: Boolean
#         @param   strict: Strict checking... EVERY REQUESTED RESIDUE HAS TO BE IN THE INTERFACE
#         @default strict: False

#         @rtype: Boolean

#         @raise: AttributeError if chain and interactor are both False
#         """
#         if residue_number_list is None: return True

#         if not chain and not interactor:
#             raise AttributeError()

#         if chain:
#             set_of_interest = self.chain_positions
#         if interactor:
#             set_of_interest = self.interactor_positions

#         set_of_requested = set([int(x) for x in residue_number_list])

#         intersection = set_of_interest.intersection(set_of_requested)

#         if len(intersection) == len(set_of_requested): return True
#         if not strict and len(intersection) > 0:       return True
#         return False

#     #
#     # OVERWRITE PARENT
#     #
#     def _build(self):
#         """
#         Builds the contacts between the two chains
#         """
#         for r1 in self.protein_chain.aminoacids + self.protein_chain.heteroatoms:
#             for r2 in self.nucleotide_chain.nucleotides + self.nucleotide_chain.heteroatoms:
#                 new_contact = r1.distance(residue = r2, dist_type = self.c_type)
#                 if new_contact.distance <= self.max_distance:
#                     if self._alldistances:
#                         if   self.c_type == "min":  altd = ["bck"]
#                         elif self.c_type == "bck":  altd = ["min"]
#                         new_contact._distance[1] = r1.distance(residue = r2, dist_type = altd[0]).distance
#                     self.add_contact(contact  = new_contact)

#     def __str__(self):
#         if not self._alldistances:
#             return super(PNInterface, self).__str__()
#         else :
#             data = []
#             # data.append('CHAIN: {0.chain_id}'.format(self))
#             for contact in self.contacts:
#                 data.append(contact.print_aldist())
#             return "\n".join(data)
