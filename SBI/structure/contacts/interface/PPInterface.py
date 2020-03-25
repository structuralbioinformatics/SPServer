"""
PPInterface

author: jbonet
date:   03/2013

@oliva's lab
"""

from .         import Interface
from ..contact import ContactAA

class PPInterface(Interface):
    """
    {Interface} between two protein chains
    """
    def __init__(self, protein_chain   = None, protein_interactor = None, 
                       threshold_type  = "cb", threshold_distance = 12,
                       protein_centers = None, interactor_centers = None):

        if protein_centers is None:
            protein_centers    = protein_chain.geometric_center(structure = True,  hetero     = False,
                                                                water     = False, by_residue = True)
        if interactor_centers is None:
            interactor_centers = protein_interactor.geometric_center(structure = True,  hetero     = False,
                                                                     water     = False, by_residue = True)

        super(PPInterface, self).__init__(chain1         = protein_chain,   chain2             = protein_interactor,
                                          threshold_type = threshold_type,  threshold_distance = threshold_distance,
                                          chain1_centres = protein_centers, chain2_centers     = interactor_centers)
        
    #
    # ATTRIBUTES
    #
    @property
    def protein_chain(self):      return self._chain1
    @property
    def protein_id(self):         return self._chain1.globalID
    @property
    def protein_centers(self):    return self._chain1_centers

    @property
    def protein_interactor(self): return self._chain2
    @property
    def interactor_id(self):      return self._chain2.globalID
    @property
    def interactor_centers(self): return self._chain1_centers

    @property
    def protein_positions(self):
        return self._list_positions(1)
    @property 
    def protein_view_interface(self):
        return self._view_interface_from(1)

    @property
    def interactor_positions(self):
        return self._list_positions(2)
    @property 
    def interactor_view_interface(self):
        return self._view_interface_from(2)


    @property
    def is_homodimer(self):
        return self._chain1.protein_sequence == self._chain2.protein_sequence

    #
    # OVERWRITE PARENT: PRIVATE FUNCTIONS
    #
    def _build(self):
        if len(self.protein_chain.aminoacids) == 0 or len(self.protein_interactor.aminoacids) == 0: return

        if (self.protein_chain.is_only_ca or self.protein_interactor.is_only_ca) and self.threshold_type == "cb":
            self._threshold_type     = "ca"
            # Add distance of two C-C bonds (1.54A) -> http://en.wikipedia.org/wiki/Bond_length
            # Gives margin by ajusting at 2
            self._threshold_distance = self.threshold_distance + (2*2)

        super(PPInterface, self)._build()

        for i in range(len(self._filtered[0])):
            new_contact = ContactAA(aminoacid1         = self.protein_chain.aminoacids[self._filtered[0][i]], 
                                    aminoacid2         = self.protein_interactor.aminoacids[self._filtered[1][i]], 
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
        data.append(super(PPInterface, self).toString(all_types))
        return "\n".join(data)


        # for r1 in self.chain.aminoacids + self.chain.heteroatoms:
        #     for r2 in self.interactor.aminoacids + self.chain.heteroatoms:
        #         new_contact = r1.distance(residue = r2, dist_type = self.c_type)
        #         if new_contact.distance <= self.max_distance:
        #             if self._alldistances:
        #                 if   self.c_type == "cb":  altd = ["ca","min"]
        #                 elif self.c_type == "ca":  altd = ["cb","min"]
        #                 elif self.c_type == "min": altd = ["cb","ca"]
        #                 new_contact._distance[1] = r1.distance(residue = r2, dist_type = altd[0]).distance
        #                 new_contact._distance[2] = r1.distance(residue = r2, dist_type = altd[1]).distance
        #             self.add_contact(contact  = new_contact)

    # def __str__(self):
    #     if not self._alldistances:
    #         return super(PPInterface, self).__str__()
    #     else :
    #         data = []
    #         # data.append('CHAIN: {0.chain_id}'.format(self))
    #         for contact in self.contacts:
    #             data.append(contact.print_aldist())
    #         return "\n".join(data)

# from ContactList import ContactList

# class PPInterface(ContactList):
#     """
#     {ContactList} between two chains
#     """
#     def __init__(self, chain = None, interactor = None, c_type = "cb", max_distance = 12, all_distances = False):
#         self._alldistances = all_distances
#         super(PPInterface, self).__init__(chain = chain, secondary_chain = interactor,
#                                           c_type = c_type, max_distance = max_distance)
        

#     #
#     # ATTRIBUTES
#     #
#     @property
#     def interactor(self):
#         """
#         Structure of the interactor chain
#         @rtype: {CHAIN}
#         """
#         return self._schain

#     @property
#     def interactor_id(self):
#         """
#         Identifier of interactor chain
#         @rtype: String
#         """
#         return self._schain.globalID()

#     @property
#     def chain_positions(self):
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
#     def interactor_positions(self):
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
#         for r1 in self.chain.aminoacids + self.chain.heteroatoms:
#             for r2 in self.interactor.aminoacids + self.chain.heteroatoms:
#                 new_contact = r1.distance(residue = r2, dist_type = self.c_type)
#                 if new_contact.distance <= self.max_distance:
#                     if self._alldistances:
#                         if   self.c_type == "cb":  altd = ["ca","min"]
#                         elif self.c_type == "ca":  altd = ["cb","min"]
#                         elif self.c_type == "min": altd = ["cb","ca"]
#                         new_contact._distance[1] = r1.distance(residue = r2, dist_type = altd[0]).distance
#                         new_contact._distance[2] = r1.distance(residue = r2, dist_type = altd[1]).distance
#                     self.add_contact(contact  = new_contact)

#     def __str__(self):
#         if not self._alldistances:
#             return super(PPInterface, self).__str__()
#         else :
#             data = []
#             # data.append('CHAIN: {0.chain_id}'.format(self))
#             for contact in self.contacts:
#                 data.append(contact.print_aldist())
#             return "\n".join(data)
