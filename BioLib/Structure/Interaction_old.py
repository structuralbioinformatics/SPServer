import sys, os

class Interaction(object):
    '''
    This class represent an interaction between two structures
    '''
 
    def __init__(self, structure1, structure2):
        '''
        Constructor
        '''
        self.structure1 = structure1
        self.structure2 = structure2
        self.positive_signatures = []
        self.negative_signatures = []

    # GETTERS #

    def get_structure1(self):
        return self.structure1

    def get_structure2(self):
        return self.structure2
   
    def get_structures(self):
        return self.structure1, self.structure2

    def get_positive_signatures(self):
        return self.positive_signatures

    def get_negative_signatures(self):
        return self.negative_signatures

    # SETTERS #

    def set_positive_signatures(self, positive_signatures):
        self.positive_signatures = positive_signatures

    def set_negative_signatures(self, negative_signatures):
        self.negative_signatures = negative_signatures

    # CLASS METHODS #

    def get_interface_signatures(self, min_coverage = 0.8, c_type = "ca", max_distance = 8, max_pval = 1, signature_type = 'all'):
        '''
        Get the signatures in the interface and the total combination of signatures (used to normalize)
        NOTE: We assume the signatures are ordered from better to worst pValue
        NOTE: We assume structure 1 contains the loops of the first protein signature: (str1 loops, str2loops, pValue)
        '''
        # Check signature type #
        signature_type_accepted = set(['all','contact','no_contact'])
        if signature_type not in signature_type_accepted:
            signature_type = 'all'
        # Init local variables #
        positive_signatures = []
        negative_signatures = []
        total_combination = 0
        # Ger interacting loops of each strcuture #
        interacting_loops_structure1 = self.get_structure1().get_interacting_loops(self.get_structure2(), min_coverage, c_type, max_distance)
        interacting_loops_structure2 = self.get_structure2().get_interacting_loops(self.get_structure1(), min_coverage, c_type, max_distance)
        # Iterate over loops combinations #
        for loop1 in interacting_loops_structure1:
            for loop2 in interacting_loops_structure2:
                if signature_type == 'all' or (signature_type == 'contact' and loop1.get_interacting_coverage(loop2) > 0) or (signature_type == 'no_contact' and loop1.get_interacting_coverage(loop2) == 0):
                    total_combination += 1
                    for positive_signature in self.get_positive_signatures():
                        if loop1.get_id() in positive_signature[0] and loop2.get_id() in positive_signature[1]:
                            if positive_signature[2] <= max_pval:
                                positive_signatures.append(positive_signature)
                                break
                    for negative_signature in self.get_negative_signatures():
                        if loop1.get_id() in negative_signature[0] and loop2.get_id() in negative_signature[1]:
                            if negative_signature[2] <= max_pval:
                                negative_signatures.append(negative_signature)
                                break
        return positive_signatures, negative_signatures, total_combination

    def get_partial_signatures(self, min_coverage = 0.8, c_type = "ca", max_distance = 8, max_pval = 1):
        '''
        Get the signatures in the partial location and the total combination of signatures (used to normalize)
        NOTE: We assume the signatures are ordered from better to worst pValue
        NOTE: We assume structure 1 contains the loops of the first protein signature: (str1 loops, str2loops, pValue)
        '''
        positive_signatures = []
        negative_signatures = []
        total_combination = 0
        interacting_loops_structure1 = self.get_structure1().get_interacting_loops(self.get_structure2(), min_coverage, c_type, max_distance)
        noninteracting_loops_structure1 = self.get_structure1().get_noninteracting_loops(self.get_structure2(), min_coverage, c_type, max_distance)
        interacting_loops_structure2 = self.get_structure2().get_interacting_loops(self.get_structure1(), min_coverage, c_type, max_distance)
        noninteracting_loops_structure2 = self.get_structure2().get_noninteracting_loops(self.get_structure1(), min_coverage, c_type, max_distance)
        for loop1 in noninteracting_loops_structure1:
            for loop2 in interacting_loops_structure2:
                total_combination += 1
                for positive_signature in self.get_positive_signatures():
                    if loop1.get_id() in positive_signature[0] and loop2.get_id() in positive_signature[1]:
                        if positive_signature[2] <= max_pval:
                            positive_signatures.append(positive_signature)
                            break
                for negative_signature in self.get_negative_signatures():
                     if loop1.get_id() in negative_signature[0] and loop2.get_id() in negative_signature[1]:
                        if negative_signature[2] <= max_pval:
                            negative_signatures.append(negative_signature)
                            break
        for loop1 in interacting_loops_structure1:
            for loop2 in noninteracting_loops_structure2:
                total_combination += 1
                for positive_signature in self.get_positive_signatures():
                    if loop1.get_id() in positive_signature[0] and loop2.get_id() in positive_signature[1]:
                        if positive_signature[2] <= max_pval:
                            positive_signatures.append(positive_signature)
                            break
                for negative_signature in self.get_negative_signatures():
                     if loop1.get_id() in negative_signature[0] and loop2.get_id() in negative_signature[1]:
                        if negative_signature[2] <= max_pval:
                            negative_signatures.append(negative_signature)
                            break
        return positive_signatures, negative_signatures, total_combination

    def get_external_signatures(self, min_coverage = 0.8, c_type = "ca", max_distance = 8, max_pval = 1):
        '''
        Get the signatures in the external location and the total combination of signatures (used to normalize)
        NOTE: We assume the signatures are ordered from better to worst pValue
        NOTE: We assume structure 1 contains the loops of the first protein signature: (str1 loops, str2loops, pValue)
        '''
        positive_signatures = []
        negative_signatures = []
        total_combination = 0
        noninteracting_loops_structure1 = self.get_structure1().get_noninteracting_loops(self.get_structure2(), min_coverage, c_type, max_distance)
        noninteracting_loops_structure2 = self.get_structure2().get_noninteracting_loops(self.get_structure1(), min_coverage, c_type, max_distance)
        for loop1 in noninteracting_loops_structure1:
            for loop2 in noninteracting_loops_structure2:
                total_combination += 1
                for positive_signature in self.get_positive_signatures():
                    if loop1.get_id() in positive_signature[0] and loop2.get_id() in positive_signature[1]:
                        if positive_signature[2] <= max_pval:
                            positive_signatures.append(positive_signature)
                            break
                for negative_signature in self.get_negative_signatures():
                     if loop1.get_id() in negative_signature[0] and loop2.get_id() in negative_signature[1]:
                        if negative_signature[2] <= max_pval:
                            negative_signatures.append(negative_signature)
                            break
        return positive_signatures, negative_signatures, total_combination 

    def get_interacting_residues(self, c_type = "ca", max_distance = 8):
        '''
        Get a list of the interacting residues of both structures and distances
        '''
        interacting_residues_structure1 = self.get_structure1().get_interacting_residues(self.get_structure2(), c_type, max_distance)
        interacting_residues_structure2 = self.get_structure2().get_interacting_residues(self.get_structure1(), c_type, max_distance)    
        return interacting_residues_structure1, interacting_residues_structure2

    def __str__(self):
        '''
        Return the interacting residues tab separated
        '''
        pass
