import math, os, sys, re
import BioLib.Tools.Dssp as Dssp
from BioLib.Tools.BioExceptions import *

class Structure(object):
    ''' 
    Structure class represents a protein chain within a PDB
    '''

    def __init__(self, pdb, chain):
        '''
        Constructor
        '''

        if chain == '': chain = 'A'

        self.id = pdb + "_" + chain
        self.pdb = pdb
        self.chain = chain
        self.residues = []
        self.loops = []
        self.uniprot_ref = None

    # GETTERS #

    def get_id(self):
        return self.id

    def get_pdb(self):
        return self.pdb

    def get_chain(self):
        return self.chain

    def get_uniprot_ref(self):
        return self.uniprot_ref

    def get_residues(self):
        return self.residues

    def get_residue_by_num(self, num):
        for residue in self.get_residues():
            if residue.get_num() == num:
                return residue
        sys.stderr.write('Residue number %d not found\n' % num)
        return None

    def get_loops(self):
        return self.loops

    def get_loop_by_id(self, subclass_id):
        for loop in self.get_loops():
            if loop.get_id() == subclass_id:
                return loop
        sys.stderr.write('Loop %s not found\n' % subclass_id)
        return None

    def get_number_residues(self):
        return len(self.get_residues())

    def get_number_atoms(self):
        atoms = 0
        for residue in self.get_residues():
            atoms += residue.get_number_atoms()
        return atoms            

    def get_first_residue(self):
        if self.get_number_residues() == 0:
            return None
        else:
            return self.get_residues()[0]

    def get_last_residue(self):
        if self.get_number_residues() == 0:
            return None
        else:
            return self.get_residues()[-1]

    # SETTERS #

    def set_uniprot_ref(self, uniprot_ref):
        self.uniprot_ref = uniprot_ref

    def set_residues(self, residues):
        self.residues = residues  
 
    def add_residue(self,residue):
        self.residues.append(residue)

    def set_loops(self, loops):
        for loop in loops:
            self.add_loop(loop)

    def add_loop(self, loop):
        loop_length = loop.get_end() - loop.get_start() + 1
        for i in xrange(loop_length):
            position = (loop.get_start()-1)+i+self.get_first_residue().get_num()
            if self.get_residue_by_num(position) != None:
                loop.add_residue(self.get_residue_by_num(position))
        loop.redefine_position((loop.get_start()-1)+self.get_first_residue().get_num(), (loop.get_end()-1)+self.get_first_residue().get_num())
        self.loops.append(loop)

    def set_dssp(self):
        '''
        Executes dssp and sets secondary structure (ss) and accesible surface area (acc) to each residue
        '''
        Dssp.execute_dssp(self)

    # CLASS METHODS #

    def clean(self):
        '''
        Cleans the structure removing residues without CA
        WARNING: DNA structures will empty
        NOTE: residues_tmp is needed to avoid remove items of the iterable during the iteration
        '''
        residues_tmp = [x for x in self.residues]
        for residue in residues_tmp:
            if not residue.has_ca():
                self.residues.remove(residue)

    def normalize_residues(self, start=1):
        '''
        Renumerate the residue numbers
        '''
        res_num = start
        for residue in self:
            residue.set_num(res_num)
            res_num += 1
        return res_num

    def normalize_atoms(self, start=1):
        '''
        Renumerate the atom numbers
        '''
        atom_num = start
        for residue in self:
            for atom in residue:
                atom.num = atom_num
                atom_num += 1
        return atom_num
        
    def get_sequence(self):
        '''
        Get the structure sequence
        '''
        seq = ''
        for residue in self:
            seq += residue.get_type_short()
        return seq 

    def get_self_interactions(self, c_type = 'CB', max_distance = 12, gap = 3, warnings = False):
        '''
        Compute the residue interactions intside the structure (folding interactions)
        '''
        interacting_residues = []
        used_residues = []
        max_distance = float(max_distance)
        for residue1 in self:
            used_residues.append(residue1)
            for residue2 in self:
                if residue2.get_num() not in [i for i in xrange(residue1.get_num()-gap, residue1.get_num()+gap+1)]:
                    try:
                        distance = residue1.get_residue_distance(residue2, c_type)
                    except ResidueDistanceError as e:
                        if warnings: sys.stderr.write(str(e))
                        continue
                    if distance <= max_distance and residue2 not in used_residues:
                        interacting_residues.append((residue1, residue2, distance))
        return interacting_residues

    def get_interacting_residues(self, structure, c_type = 'CB', max_distance = 12, uniq=True, warnings = False):
        '''
        Compute the interacting residues between two structure
        @Return List of interacting residues of both proteins. If unique is false, returns interacting
        pairs and their distances.
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

    def get_opposite_interacting_residues(self, structure, size_sef = 1, size_struct = 1, c_type = 'CB', max_distance = 12, warnings=False):
        '''
        Returns a set of residues (for each structure) in the opposite faces of the interaction sites
        Used to compute back to back and partial interfaces
        '''
        interacting_residues = self.get_interacting_residues(structure, c_type=c_type, max_distance=max_distance, uniq=True, warnings=warnings)
        self_mass_center = self.__get_mass_center(interacting_residues[0])
        structure_mass_center = self.__get_mass_center(interacting_residues[1])
        self_residues_distances = []
        structure_residues_distances = []
        for residue in self:
            self_residues_distances.append((residue.get_ca().get_distance_coords(self_mass_center), residue))
        for residue in structure:
            structure_residues_distances.append((residue.get_ca().get_distance_coords(self_mass_center), residue))
        self_residues_distances.sort(reverse=True)
        structure_residues_distances.sort(reverse=True)
        self_residues_distances = [x[1] for x in self_residues_distances]
        structure_residues_distances = [x[1] for x in structure_residues_distances]
        return self_residues_distances[0:size_sef], structure_residues_distances[0:size_struct]

    def __get_mass_center(self, residues):
        '''
        Gets the mass center of a residue iterable (structure, list of residues...)
        '''
        x = 0
        y = 0
        z = 0
        for residue in residues:
            x += residue.get_ca().get_x()
            y += residue.get_ca().get_y()
            z += residue.get_ca().get_z()

        if len(residues) == 0:
            return None
        return x/len(residues), y/len(residues), z/len(residues)

    def residue_interacts(self, residue_num, structure, c_type = 'CB', max_distance = 12, warnings=False):
        '''
        Checks if a certain residue (specify by the number) interacts with other structure
        '''
        residue1 = self.get_residue_by_num(residue_num)
        interacts = False
        for residue2 in structure:
            try:
                distance = residue1.get_residue_distance(residue2, c_type)
            except ResidueDistanceError as e:
                if warnings: sys.stderr.write(str(e))
                continue
            if distance <= max_distance:
                interacts = True
                break
        return interacts

    def get_interacting_loops(self, structure, min_coverage = 0.8, c_type = 'CB', max_distance = 12, warnings=False):
        '''
        Get the interacting loops with other structure
        '''
        interacting_loops = []
        for loop in self.get_loops():
            coverage = loop.get_interacting_coverage(structure, c_type, max_distance, warnings)
            if min_coverage <= coverage:
                interacting_loops.append(loop)
        return interacting_loops

    def get_noninteracting_loops(self, structure, min_coverage = 0.8, c_type = 'CB', max_distance = 12, warnings=False):
        '''
        Get the non-interacting loops with other structure
        '''
        interacting_loops = self.get_interacting_loops(structure, min_coverage, c_type, max_distance, warnings)
        return list(set(self.get_loops()).difference(set(interacting_loops)))

    def relocate(self, translation_vector, rotation_vector, docking='HEX'):
        '''
        Relocate the structure position 
        @translation_vector = Vector used to translate the structure
        @Rotation_vector = Vector used to rotate the structure
        @docking = Docking program used to get the vectors
        '''
        for residue in self:
            try:
                residue.relocate(translation_vector, rotation_vector, docking)
            except RelocateProgramError as e:
                raise e

    def get_RMSD(self, structure):
        if self.get_number_residues() != structure.get_number_residues():
            sys.stderr.write("The RMSD needs to be between the same structures!!\n")
        E = 0
        for residue in xrange(self.get_number_residues()):
            for atom in xrange(self.get_residues()[residue].get_number_atoms()):
                distance = self.get_residues()[residue].get_atoms()[atom].get_distance(atom=structure.get_residues()[residue].get_atoms()[atom])
                E += math.pow(distance, 2)
        E = E /self.get_number_atoms()

        return math.sqrt(E)

    def get_radius(self):
        '''
        Get structure maximum radius (as defined in FTDock)
        '''
        largest = 0
        for residue in self:
            for atom in residue:
                present = math.pow(atom.get_x(), 2) + math.pow(atom.get_y(), 2) + math.pow(atom.get_z(), 2)
                if present > largest:
                    largest = present

        return math.sqrt(largest)

    def translate_onto_origin(self):
        '''
        Translate the structure to its center and return the translation vector
        '''
        atoms_num = 0
        x = 0
        y = 0
        z = 0

        for residue in self:
            for atom in residue:
                x = x + atom.get_x()
                y = y + atom.get_y()
                z = z + atom.get_z()
                atoms_num += 1

        x = x/atoms_num
        y = y/atoms_num
        z = z/atoms_num

        self.relocate((-x,-y,-z), (0,0,0))

    def get_N_linker(self):
        '''
        Get a list of N-terminal linker residues
        '''
        NLinker = []
        for residue in self:
            if residue.get_ss() == 'C':
                NLinker.append(residue)
            else:
               break
        return NLinker
 
    def get_N_linker_length(self):
        '''
        Get the length of the N-terminal linker
        '''
        return len(self.get_N_linker())

    def get_C_linker(self):
        '''
        Get a list of C-terminal linker residues
        '''
        CLinker = []
        for residue in self:
            if residue.get_ss() != 'C':
                CLinker = []
            else:
                CLinker.append(residue)
        return CLinker

    def get_C_linker_length(self):
        '''
        Get the length of the C-terminal linker
        '''
        return len(self.get_C_linker())

    def get_first_ss(self):
        '''
        Get the first residue with secondary structure
        '''
        firstResidue = None
        for residue in self:
            if residue.get_ss() != 'C':
                firstResidue = residue
                break
        return firstResidue

    def get_last_ss(self):
        '''
        Get the last residue with secondary structure
        '''
        lastResidue = None
        for residue in self:
            if residue.get_ss() != 'C':
                lastResidue = residue
        return lastResidue

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
        '''
        Returns the structure in PDB format
        '''
        structurePDB = ''
        for residue in self:
            for atom in residue:
                structurePDB += 'ATOM%s  %s%s%s%s%s%s%s\n' % (str(atom.get_num()).rjust(7),
                                                            atom.get_type().ljust(4),
                                                            residue.get_type().rjust(3),
                                                            self.get_chain().rjust(2),
                                                            str(residue.get_num()).rjust(4),
                                                            ('%.3f' %float(atom.get_x())).rjust(12),
                                                            ('%.3f' %float(atom.get_y())).rjust(8),
                                                            ('%.3f' %float(atom.get_z())).rjust(8))
        structurePDB += 'TER\n' 
        return structurePDB
