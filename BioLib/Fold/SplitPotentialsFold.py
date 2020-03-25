import os, sys, random, math
import numpy
from BioLib.Tools.BioExceptions import ResidueTriadError

class SplitPotentialsFold(object):
    '''
    Split Potentials for protein Folding. 
    Applications: Evaluate models, rank them and find incorrect folded regions.
    '''
    def __init__(self, c_type = 'CB', cutoff = 12):
        '''
        Contructor
        '''
        c_type_accepted = ('CB','MIN')
        if c_type.upper() not in c_type_accepted:
            sys.stderr.write('WARNING: Incorrect contact type %s, using CB!' % c_type)
            c_type = 'CB'

        if c_type == 'CB':
            self.sp_file = os.path.abspath(os.path.join(os.path.dirname(__file__), 'fold3_spf4_cbR12scr1gap4.out'))
        else:
            self.sp_file = os.path.abspath(os.path.join(os.path.dirname(__file__), 'fold3_spf4_minR5scr1gap4.out'))
        if not os.path.isfile(self.sp_file):
            sys.stderr.write('ERROR: Can not find the split potential file %s\n\n' % self.sp_file)
            sys.exit()

        self.c_type = c_type.upper()  # Contact type
        self.cutoff = cutoff          # Standard cutoffs are 12 for CB and 5 for MIN

        self.ppR     = {}             # Pair Potentials of Residues
        self.ppE     = {}             # Pair Potentials of Environment
        self.LppRE   = {}             # Local Pair Potentials of Residues and Environment
        self.ppRE    = {}             # Pair Potentials of Residues and Environment
        self.ppD     = {}             # Standard Pair Potentials of Distances

        self.read_sp_file()           # Read the file and load the matrices
        
    ##################
    # PUBLIC METHODS #
    ##################

    def calculate_global_energies(self, structure, Zscores=False, randoms=100):
        '''
        Compute the split potentials energies for a protein fold
        if Zscores: Computes also the Zscores
        '''
        # Compute global energies
        real_interactions = structure.get_self_interactions(c_type=self.c_type, max_distance=self.cutoff+1, gap=3)
        global_energies = self.calculate_energies(real_interactions)
        if not Zscores:
            return global_energies
        # Generate random fold energies
        rand_energies_lists = {}
        for random in xrange(randoms):
            rand_interactions = self.randomize_residue_interactions(structure.get_residues(), real_interactions)
            rand_energies = self.calculate_energies(rand_interactions)
            for potential_type in rand_energies:
                rand_energies_lists.setdefault(potential_type, []).append(rand_energies[potential_type])
        # Compute Zscores
        global_Zenergies = self.compute_Zscores(global_energies, rand_energies_lists)
        return global_energies, global_Zenergies

    def calculate_residue_energies(self, structure, Zscores=False, randoms=100):
        '''
        Compute the split potential energies for each residue in a protein fold
        if Zscores: Computes also the Zscores
        '''
        # Compute energies for each residue
        real_interactions = structure.get_self_interactions(c_type=self.c_type, max_distance=self.cutoff+1, gap=3)
        residues_energies = []
        residues_interactions = []
        for residue in structure:
            residue_interactions = []
            for real_interaction in real_interactions:
                if residue == real_interaction[0]:
                    residue_interactions.append(real_interaction)
                elif residue == real_interaction[1]:
                    residue_interactions.append((real_interaction[1], real_interaction[0], real_interaction[2]))
            residue_energies = self.calculate_energies(residue_interactions, local_type='R')
            residues_energies.append((residue, residue_energies))
            residues_interactions.append((residue, residue_interactions, residue_energies))
        if not Zscores:
            return residues_energies
        # Generate Zscore energies for each residues
        residues_energies = []
        for residue in residues_interactions:
            rand_energies_lists = {}
            for random in xrange(randoms):
                rand_interactions = self.randomize_residue_interactions(structure.get_residues(), residue[1])
                rand_energies = self.calculate_energies(rand_interactions, local_type='R')
                for potential_type in rand_energies:
                    rand_energies_lists.setdefault(potential_type, []).append(rand_energies[potential_type])
            residues_energies.append((residue[0], residue[2], self.compute_Zscores(residue[2], rand_energies_lists)))
        return residues_energies

    ########################################################
    # METHODS TO COMPUTE SPLIT POTENTIALS FOR INTERACTIONS #
    ########################################################

    def calculate_energies(self, res_interactions, local_type='G'):
        '''
        Compute the split potential energies for a set of residue-residue interactions
        '''
        local_type_accepted = ('G','R') # Global or Residue calculations
        if local_type.upper() not in local_type_accepted:
            sys.stderr.write('WARNING: Incorrect local type %s, using G!' % local_type)
            local_type = 'G'

        # Init Energy values
        energies = {'D-PAIR':   0, 'C-PAIR':  0, 'U-PAIR':  0, 'L-PAIR':  0, 'M-PAIR':  0,
                    'D-3DC':    0, 'C-3DC':   0, 'U-3DC':   0, 'L-3DC':   0, 'M-3DC':   0,
                    'D-LOCAL':  0,
                    'D-S3DC':   0, 'C-S3DC':  0, 'U-S3DC':  0, 'L-S3DC':  0, 'M-S3DC':  0,
                    'D-3D':     0, 'C-3D':    0,
                    'D-COMB':   0, 'C-COMB':  0, 'U-COMB':  0, 'L-COMB':  0, 'M-COMB':  0}

        # Compute energy for each residue interaction
        for res_interaction in res_interactions:
            residue1 = res_interaction[0]
            residue2 = res_interaction[1]
            distance = res_interaction[2]
            try:
                triad1 = residue1.get_triad()
                triad2 = residue2.get_triad()
            except ResidueTriadError as e:
                sys.stdout.write(str(e))
                sys.stdout.flush()
                continue

            # Energy of Pair Potentials of Residues #
            distance_ppR           = self.get_ppR(triad1[0:1], triad2[0:1], distance)
            energies['D-PAIR']    += distance_ppR
            cutoff_ppR             = self.get_ppR(triad1[0:1], triad2[0:1], self.cutoff)
            energies['C-PAIR']    += cutoff_ppR
            maximum_ppR            = self.get_ppR_stat(triad1[0:1], triad2[0:1], stat='upper')
            energies['U-PAIR']    += maximum_ppR
            minimum_ppR            = self.get_ppR_stat(triad1[0:1], triad2[0:1], stat='lower')
            energies['L-PAIR']    += minimum_ppR
            mean_ppR               = self.get_ppR_stat(triad1[0:1], triad2[0:1], stat='mean')
            energies['M-PAIR']    += mean_ppR
            # Energy of Pair Potentials of Environment
            distance_ppE           = self.get_ppE(triad1[2:], triad2[2:], distance)
            energies['D-3DC']     += distance_ppE
            cutoff_ppE             = self.get_ppE(triad1[2:], triad2[2:], self.cutoff)
            energies['C-3DC']     += cutoff_ppE
            maximum_ppE            = self.get_ppE_stat(triad1[2:], triad2[2:], stat='upper')
            energies['U-3DC']     += maximum_ppE
            minimum_ppE            = self.get_ppE_stat(triad1[2:], triad2[2:], stat='lower')
            energies['L-3DC']     += minimum_ppE
            mean_ppE               = self.get_ppE_stat(triad1[2:], triad2[2:], stat='mean')
            energies['M-3DC']     += mean_ppE
            # Energy of Local Pair Potentials of Residues and Environment
            if local_type == 'R':
                distance_LppRE     = self.get_LppRE(triad1)
            else:
                distance_LppRE     = self.get_LppRE(triad1) + self.get_LppRE(triad2)
            energies['D-LOCAL']   += distance_LppRE
            # Energy of Pair Potentials of Residues and Environment
            distance_ppRE          = self.get_ppRE(triad1, triad2, distance)
            energies['D-S3DC']    += distance_ppRE
            cutoff_ppRE            = self.get_ppRE(triad1, triad2, self.cutoff)
            energies['C-S3DC']    += cutoff_ppRE
            maximum_ppRE           = self.get_ppRE_stat(triad1, triad2, stat='upper')
            energies['U-S3DC']    += maximum_ppRE
            minimum_ppRE           = self.get_ppRE_stat(triad1, triad2, stat='lower')
            energies['L-S3DC']    += minimum_ppRE
            mean_ppRE              = self.get_ppRE_stat(triad1, triad2, stat='mean')
            energies['M-S3DC']    += mean_ppRE
            # Energy of Standard Pair Potentials of Distances
            distance_ppD           = self.get_ppD(distance)
            energies['D-3D']      += distance_ppD
            cutoff_ppD             = self.get_ppD(self.cutoff)
            energies['C-3D']      += cutoff_ppD
            # Energy of Combined energy
            distance_comb          = distance_ppRE + distance_ppD + distance_LppRE + distance_ppE
            energies['D-COMB']  += distance_comb
            cutoff_comb            = cutoff_ppRE + cutoff_ppD + distance_LppRE + cutoff_ppE
            energies['C-COMB']  += cutoff_comb
            maximum_comb           = maximum_ppRE + distance_ppD + distance_LppRE + maximum_ppE
            energies['U-COMB']  += maximum_comb
            minimum_comb           = minimum_ppRE + distance_ppD + distance_LppRE + minimum_ppE
            energies['L-COMB']  += minimum_comb
            mean_comb              = mean_ppRE + distance_ppD + distance_LppRE + mean_ppE
            energies['M-COMB']  += mean_comb

        return energies

    ##################################
    # CONTACT RANDOMIZATION  METHODS #
    ##################################

    def randomize_residue_interactions(self, struct_residues, res_interactions):

        rand_interactions = []
        done = {}
        for res_interaction in res_interactions:
            dist = res_interaction[2]
            if done.has_key(res_interaction[0]):
                res1 = done[res_interaction[0]]
            else:
                res1 = random.choice(struct_residues)
                done[res_interaction[0]] = res1
            if done.has_key(res_interaction[1]):
                res2 = done[res_interaction[1]]
            else:
                res2 = random.choice(struct_residues)
                done[res_interaction[1]] = res2
            rand_interactions.append((res1, res2, dist))
        return rand_interactions

    def compute_Zscores(self, real_energies, rand_energies_lists):

        Zenergies = {} 
        for potential_type in real_energies:
            if potential_type == 'D-3D' or potential_type == 'C-3D':
                continue
            energies_mean = numpy.mean(rand_energies_lists[potential_type])
            energies_std  = numpy.std(rand_energies_lists[potential_type])
            Zenergies[potential_type] = (real_energies[potential_type]-energies_mean)/energies_std
        return Zenergies

    ############################
    # SPLIT POTENTIALS GETTERS #
    ############################

    def get_ppR(self, res1, res2, distance):

        distance = int(math.floor(distance))
        if self.ppR.has_key(res1) and self.ppR[res1].has_key(res2):
            return float(self.ppR[res1][res2][distance])
        elif self.ppR.has_key(res2) and self.ppR[res2].has_key(res1):
            return float(self.ppR[res2][res1][distance])
        else:
            return 0
            #raise RuntimeError('ERROR: Matrix entry not found: %s, %s' % (res1, res2))
        
    def get_ppR_stat(self, res1, res2, stat):

        if self.ppR.has_key(res1) and self.ppR[res1].has_key(res2):
            data = numpy.array(self.ppR[res1][res2], dtype='float64')
        elif self.ppR.has_key(res2) and self.ppR[res2].has_key(res1):
            data = numpy.array(self.ppR[res2][res1], dtype='float64')
        else:
            return 0
            #raise RuntimeError('ERROR: Matrix entry not found: %s, %s, %d' % (res1, res2))
        if stat.upper() == 'UPPER':
            return numpy.max(data)
        elif stat.upper() == 'LOWER':
            return numpy.min(data)
        elif stat.upper() == 'MEAN':
            return numpy.mean(data)
        
    def get_ppE(self, res1, res2, distance):

        distance = int(math.floor(distance))
        if self.ppE.has_key(res1) and self.ppE[res1].has_key(res2):
            return -float(self.ppE[res1][res2][distance])
        elif self.ppE.has_key(res2) and self.ppE[res2].has_key(res1):
            return -float(self.ppE[res2][res1][distance])
        else:
            return 0
            #raise RuntimeError('ERROR: Matrix entry not found: %s, %s' % (res1, res2))

    def get_ppE_stat(self, res1, res2, stat):

        if self.ppE.has_key(res1) and self.ppE[res1].has_key(res2):
            data = numpy.array(self.ppE[res1][res2], dtype='float64')
        elif self.ppE.has_key(res2) and self.ppE[res2].has_key(res1):
            data = numpy.array(self.ppE[res2][res1], dtype='float64')
        else:
            return 0
            #raise RuntimeError('ERROR: Matrix entry not found: %s, %s' % (res1, res2))
        # UPPER and LOWER are inverted in this potential
        if stat.upper() == 'UPPER':
            return -numpy.min(data)
        elif stat.upper() == 'LOWER':
            return -numpy.max(data)
        elif stat.upper() == 'MEAN':
            return -numpy.mean(data)
              
    def get_LppRE(self, res):

        return -float(self.LppRE[res])
    
    def get_ppRE(self, res1, res2, distance):
        '''
        '''
        distance = int(math.floor(distance))
        if self.ppRE.has_key(res1) and self.ppRE[res1].has_key(res2):
            return float(self.ppRE[res1][res2][distance])
        elif self.ppRE.has_key(res2) and self.ppRE[res2].has_key(res1):
            return float(self.ppRE[res2][res1][distance])
        else:
            return 0
            #raise RuntimeError('ERROR: Matrix entry not found: %s, %s' % (res1, res2))

    def get_ppRE_stat(self, res1, res2, stat):
        '''
        '''
        if self.ppRE.has_key(res1) and self.ppRE[res1].has_key(res2):
            data = numpy.array(self.ppRE[res1][res2], dtype='float64')
        elif self.ppRE.has_key(res2) and self.ppRE[res2].has_key(res1):
            data = numpy.array(self.ppRE[res2][res1], dtype='float64')
        else:
            return 0
            #raise RuntimeError('ERROR: Matrix entry not found: %s, %s' % (res1, res2))
        if stat.upper() == 'UPPER':
            return numpy.max(data)
        elif stat.upper() == 'LOWER':
            return numpy.min(data)
        elif stat.upper() == 'MEAN':
            return numpy.mean(data)
        
    def get_ppD(self, distance):
        '''
        '''
        distance = str(int(math.floor(distance)))
        return float(self.ppD[distance])

    #########################################################
    # METHODS TO READ AND PROCESS THE SPLIT POTENTIALS FILE #
    #########################################################

    def read_sp_file(self):

        tag = 0
        sp_file_fo = open(self.sp_file, 'r')
        for line in sp_file_fo:
            if line.startswith("PAIR-POTENTIALS OF RESIDUES AND ENVIRONMENTS"):             tag = 4
            elif line.startswith("PAIR-POTENTIALS OF ENVIRONMENTS"):                        tag = 2
            elif line.startswith("LOCAL PAIR-POTENTIALS OF RESIDUES AND ENVIRONMENTS"):     tag = 3
            elif line.startswith("PAIR-POTENTIALS OF RESIDUES"):                            tag = 1
            elif line.startswith("STANDARD PAIR-POTENTIALS OF DISTANCES"):                  tag = 5
            elif line.startswith("END"):                                                    tag = 0
            self.process_sp_file_line(line, tag)
        sp_file_fo.close()

    def process_sp_file_line(self, line, tag):

        if tag == 0:    
            return
        if tag == 1 and line.startswith(" PAIR"): 
            self.assign_multiple_distances(self.ppR, line.split()[1], line.split()[2:])
        if tag == 2 and line.startswith(" PAIR"): 
            self.assign_multiple_distances(self.ppE, line.split()[1], line.split()[2:])
        if tag == 3 and line.startswith(" AA-ENV"): 
            self.assign_single_distance(self.LppRE, line.split()[1], line.split()[2])
        if tag == 4 and line.startswith(" PAIR"): 
            self.assign_multiple_distances(self.ppRE, line.split()[1], line.split()[2:])
        if tag == 5 and line.startswith(" Distance"): 
            self.assign_single_distance(self.ppD, line.split()[1], line.split()[2])

    def assign_multiple_distances(self, matrix, keys, values):
        
        matrix.setdefault(keys.split(':')[0],{}).setdefault(keys.split(':')[1],[]).extend(values)
        
    def assign_single_distance(self, matrix, key, value):
        
        matrix[key] = value

    ##################################################
    # METHODS TO PRINT THE SPLIT POTENTIALS MATRICES #
    ##################################################

    def two_keys_matrix_2_string_array(self, matrix, prev):
        
        lines = []
        for res1 in matrix:
            for res2 in matrix[res1]:
                line = prev + "%s:%s\t" % (res1, res2)
                line = line + "\t".join(matrix[res1][res2])
                lines.append(line)
        return lines
                
    def one_key_matrix_2_string_array(self, matrix, prev):
        
        lines = []
        for k in matrix:
            line = prev + "\t".join((k,matrix[k]))
            lines.append(line)
        return lines

    def ppR_2_string_array(self):
        
        lines = []
        lines.append("PAIR-POTENTIALS OF RESIDUES")
        lines.extend(self.two_keys_matrix_2_string_array(self.ppR, " PAIR "))
        lines.append("END")
        return lines
    
    def string_ppR(self): 

        return "\n".join(self.ppR_2_string_array())

    def ppE_2_string_array(self):
        
        lines = []
        lines.append("PAIR-POTENTIALS OF ENVIRONMENTS")
        lines.extend(self.two_keys_matrix_2_string_array(self.ppE, " PAIR "))
        lines.append("END")
        return lines
    
    def string_ppE(self): 

        return "\n".join(self.ppE_2_string_array())
    
    def LppRE_2_string_array(self):
        
        lines = []
        lines.append("LOCAL PAIR-POTENTIALS OF RESIDUES AND ENVIRONMENTS")
        lines.extend(self.one_key_matrix_2_string_array(self.LppRE, " AA-ENV "))
        lines.append("END")
        return lines
    
    def string_LppRE(self): 

        return "\n".join(self.LppRE_2_string_array())

    def ppRE_2_string_array(self):
        
        lines = []
        lines.append("PAIR-POTENTIALS OF RESIDUES AND ENVIRONMENTS")
        lines.extend(self.two_keys_matrix_2_string_array(self.ppRE, " PAIR "))
        lines.append("END")
        return lines
    
    def string_ppRE(self): 

        return "\n".join(self.ppRE_2_string_array())

    def ppD_2_string_array(self):
        
        lines = []
        lines.append("STANDARD PAIR-POTENTIALS OF DISTANCES")
        lines.extend(self.one_key_matrix_2_string_array(self.ppD, " Distance "))
        lines.append("END")
        return lines
    
    def string_ppD(self): 

        return "\n".join(self.ppD_2_string_array())

    def __str__(self):
        
        lines = []
        lines.extend(self.ppR_2string_array())
        lines.extend(self.ppE_2string_array())
        lines.extend(self.LppRE_2string_array())
        lines.extend(self.ppRE_2string_array())
        lines.extend(self.ppD_2string_array())
        return "\n".join(lines)
        