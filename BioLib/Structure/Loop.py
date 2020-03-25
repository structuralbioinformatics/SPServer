from BioLib.Tools.BioExceptions import *

class Loop(object):
    '''
    This class represents a loop
    '''
    def __init__(self, loop_id, start, end):
        self.loop_id = loop_id
        self.seq_start = start
        self.struct_start = None
        self.seq_end = end
        self.struct_end = None
        self.residues = []

    # GETTERS #

    def get_id(self):
        return self.loop_id
    def get_start(self):
        if self.struct_start:
            return self.struct_start
        return self.seq_start
    def get_end(self):
        if self.struct_end:
            return self.struct_end
        return self.seq_end
    def get_residues(self):
        return self.residues

    # SETTERS #   

    def set_residues(self, residues):
        self.residues = residues
    def add_residue(self, residue):
        self.residues.append(residue)
    def redefine_position(self, start, end):
        self.struct_start = start
        self.struct_end = end

    # CLASS METHODS #

    def get_interacting_residues(self, structure, c_type = 'CB', max_distance = 12, uniq = True, warnings = False):
        '''
        Get the interacting residues with an structure (that can be a Structure 
        object or a Loop object since they share __iter__ method)
        '''
        max_distance = float(max_distance)
        interacting_residues_self = []
        interacting_residues_structure = []
        interacting_distances = []
        for residue1 in self:
            for residue2 in structure:
                try:
                    distance = residue1.get_residue_distance(residue2, c_type)
                except ResidueDistanceError as e:
                    if warnings: sys.stderr.write(str(e))
                    continue
                if distance <= max_distance:
                    interacting_residues_self.append(residue1)
                    interacting_residues_structure.append(residue2)
                    interacting_distances.append(distance)
        if uniq:
            return list(set(interacting_residues_self)), list(set(interacting_residues_structure))
        else:
            return interacting_residues_self, interacting_residues_structure, interacting_distances

    def get_interacting_coverage(self, structure, c_type = 'CB', max_distance = 12, warnings=False):
        '''
        Get the coverage of a loop interacting with an structure (that can be a Structure
        object or a Loop object since they share get_residues method)
        '''
        total_residues = float(len(self))
        interacting_residues = len(self.get_interacting_residues(structure = structure, c_type = c_type, max_distance = max_distance, uniq = True, warnings = warnings)[0])
        return interacting_residues/total_residues

    def __len__(self):
        '''
        Length of residues
        '''
        return len(self.residues)

    def __iter__(self):
        '''
        Iterates over each residue
        '''
        for residue in self.residues:
            yield residue

    def __str__(self):
        return "%s, %d-%d" % (self.loop_id, self.get_start(), self.get_end())
