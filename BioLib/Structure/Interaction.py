import sys, os
from collections import Counter

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

    def get_signatures_by_location(self, min_coverage = 0.8, c_type = 'CB', max_distance = 12, max_pval = 1):
        '''
        Get the signatures in the interface and the total combination of signatures
        NOTE: We assume the signatures are ordered from better to worst pValue
        NOTE: We assume structure 1 contains the loops of the first protein signature: (str1 loops, str2loops, pValue)
        NOTE: I: interface
              E: external
              P: partial
        '''
        # Get interacting loops of each strcuture #
        int_loops_structure1 = [i.get_id() for i in self.get_structure1().get_interacting_loops(self.get_structure2(), min_coverage, c_type, max_distance)]
        int_loops_structure2 = [i.get_id() for i in self.get_structure2().get_interacting_loops(self.get_structure1(), min_coverage, c_type, max_distance)]
        nonint_loops_structure1 = [i.get_id() for i in self.get_structure1().get_noninteracting_loops(self.get_structure2(), min_coverage, c_type, max_distance)]
        nonint_loops_structure2 = [i.get_id() for i in self.get_structure2().get_noninteracting_loops(self.get_structure1(), min_coverage, c_type, max_distance)]
        # Init signatures datatype #
        signatures = {'I': ([], [], len(int_loops_structure1)*len(int_loops_structure2)), 
                      'E' : ([], [], len(nonint_loops_structure1)*len(nonint_loops_structure1)),
                      'P'  : ([], [], len(int_loops_structure1)*len(nonint_loops_structure2)+len(int_loops_structure2)*len(nonint_loops_structure1))}
        # Iterate over signatures # 
        for positive_signature in self.get_positive_signatures():
            in_int_structure1 = not Counter(positive_signature[0])-Counter(int_loops_structure1)
            in_int_structure2 = not Counter(positive_signature[1])-Counter(int_loops_structure2)
            in_noint_structure1 = not Counter(positive_signature[0])-Counter(nonint_loops_structure1)
            in_noint_structure2 = not Counter(positive_signature[1])-Counter(nonint_loops_structure2)
            if in_int_structure1 and in_int_structure2 and positive_signature[2] <= max_pval:
                signatures['I'][0].append(positive_signature)
            if in_noint_structure1 and in_noint_structure2 and positive_signature[2] <= max_pval:
                signatures['E'][0].append(positive_signature)
            if in_noint_structure1 and in_int_structure2 and positive_signature[2] <= max_pval:
                signatures['P'][0].append(positive_signature)
            elif in_int_structure1 and in_noint_structure2 and positive_signature[2] <= max_pval:
                signatures['P'][0].append(positive_signature)
        for negative_signature in self.get_negative_signatures():
            in_int_structure1 = not Counter(negative_signature[0])-Counter(int_loops_structure1)
            in_int_structure2 = not Counter(negative_signature[1])-Counter(int_loops_structure2)
            in_noint_structure1 = not Counter(negative_signature[0])-Counter(nonint_loops_structure1)
            in_noint_structure2 = not Counter(negative_signature[1])-Counter(nonint_loops_structure2)
            if in_int_structure1 and in_int_structure2 and negative_signature[2] <= max_pval:
                signatures['I'][1].append(negative_signature)
            if in_noint_structure1 and in_noint_structure2 and negative_signature[2] <= max_pval:
                signatures['E'][1].append(negative_signature)
            if in_noint_structure1 and in_int_structure2 and negative_signature[2] <= max_pval:
                signatures['P'][1].append(negative_signature)
            elif in_int_structure1 and in_noint_structure2 and negative_signature[2] <= max_pval:
                signatures['P'][1].append(negative_signature)
        # Return signatures and total combination
        return signatures
        
    def get_interacting_residues(self, c_type = 'CB', max_distance = 12, uniq=True, warnings=False):
        '''
        Get a list of the interacting residues of both structures and distances
        '''
        return self.get_structure1().get_interacting_residues(self.get_structure2(), c_type, max_distance, uniq, warnings)

    def get_noninteracting_residues(self, c_type = 'CB', max_distance = 12, warnings=False):
        '''
        Get a list of the noninteracting residues of both structures
        '''
        interacting_residues = self.get_structure1().get_interacting_residues(self.get_structure2(), c_type, max_distance, True, warnings)
        return list(set(self.get_structure1().get_residues()).difference(set(interacting_residues[0]))), list(set(self.get_structure2().get_residues()).difference(set(interacting_residues[1])))

    def get_opposite_interacting_residues(self, size_str1 = 1, size_str2 = 1, c_type = 'CB', max_distance = 12, warnings=False):
        '''
        Get a list of residues in the opposite face of the interaction of both structures
        '''
        return self.get_structure1().get_opposite_interacting_residues(self.get_structure2(), size_str1, size_str2, c_type, max_distance, warnings)

    def __str__(self):
        '''
        Return the interacting residues
        '''
        interacting_residues = self.get_structure1().get_interacting_residues(self.get_structure2(), uniq=True)
        interactionPDB = '%s\n' % self.get_structure1().get_id()
        for residue in interacting_residues[0]:
            interactionPDB += residue.__str__()
        interactionPDB = '%s\n' % self.get_structure2().get_id()
        for residue in interacting_residues[1]:
            interactionPDB += residue.__str__()
        return interactionPDB
