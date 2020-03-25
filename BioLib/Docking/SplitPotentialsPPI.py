import os, sys, random, math
import numpy
from BioLib.Tools.BioExceptions import ResidueTriadError

class SplitPotentialsPPI(object):
    '''
    Split Potentials for protein-protein interctions. 
    Applications: Evaluate ppi, rank them (docking) and find the hot spots.
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
            self.sp_file = os.path.abspath(os.path.join(os.path.dirname(__file__), 'ppi3_cbR12.out'))
        else:
            self.sp_file = os.path.abspath(os.path.join(os.path.dirname(__file__), 'ppi3_minR5.out'))
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

    def calculate_global_energies(self, interaction, Zscores=False, randoms=100, funnel=False):
        '''
        Compute the split potentials energies for a protein-protein interaction
        if Zscores: Computes also the Zscores
        if funnel: Compute also the external and partial energies (shuffling the interface with non-interacting exposed residues)
        '''
        # Compute global energies
        real_interactions = interaction.get_interacting_residues(c_type=self.c_type, max_distance=self.cutoff+1, uniq=False)
        real_energies = self.calculate_energies(zip(real_interactions[0], real_interactions[1], real_interactions[2]))
        interaction_info = [real_energies]
        if funnel:
            noninterating_exposedR = self.__get_noninteracting_residues(list(set(real_interactions[0])), interaction.get_structure1())
            noninterating_exposedL = self.__get_noninteracting_residues(list(set(real_interactions[1])), interaction.get_structure2())
            if not noninterating_exposedR or not noninterating_exposedL:
                raise RuntimeError('ERROR: Non noninteracting exposed residues (bad ACC assignation?)')
            native_energies_lists = {}
            external_energies_lists = {}
            partial_energies_lists = {}
            native_energies = {}
            external_energies = {}
            partial_energies = {}
            for random in xrange(randoms):
                rand_nonint_R = self.randomize_residue_interface(noninterating_exposedR, real_interactions[0])
                rand_nonint_L = self.randomize_residue_interface(noninterating_exposedL, real_interactions[1])
                rand_int_R = self.randomize_residue_interface(list(set(real_interactions[0])), real_interactions[0])
                rand_int_L = self.randomize_residue_interface(list(set(real_interactions[1])), real_interactions[1])
                rand_native_energies = self.calculate_energies(zip(rand_int_R, rand_int_L, real_interactions[2], real_interactions[0], real_interactions[1]), randomized=True)
                rand_partialR_energies = self.calculate_energies(zip(rand_int_R, rand_nonint_L, real_interactions[2], real_interactions[0], real_interactions[1]), randomized=True)
                rand_partialL_energies = self.calculate_energies(zip(rand_nonint_R, rand_int_L, real_interactions[2], real_interactions[0], real_interactions[1]), randomized=True)
                rand_external_energies = self.calculate_energies(zip(rand_nonint_R, rand_nonint_L, real_interactions[2], real_interactions[0], real_interactions[1]), randomized=True)
                for potential_type in rand_native_energies:
                    native_energies_lists.setdefault(potential_type, []).append(rand_native_energies[potential_type])  
                for potential_type in rand_partialR_energies:
                    partial_energies_lists.setdefault(potential_type, []).append(rand_partialR_energies[potential_type])
                for potential_type in rand_partialL_energies:
                    partial_energies_lists.setdefault(potential_type, []).append(rand_partialL_energies[potential_type]) 
                for potential_type in rand_external_energies:
                    external_energies_lists.setdefault(potential_type, []).append(rand_external_energies[potential_type])
            for potential in native_energies_lists:
                native_energies[potential] = numpy.mean(native_energies_lists[potential])
            for potential in partial_energies_lists:
                partial_energies[potential] = numpy.mean(partial_energies_lists[potential])
            for potential in external_energies_lists:
                external_energies[potential] = numpy.mean(external_energies_lists[potential])
            interaction_info.extend([native_energies, partial_energies, external_energies])
        if not Zscores:
            return interaction_info
        # Generate random fold energies
        rand_energies_lists = {}
        for random in xrange(randoms):
            rand_interface_R = self.randomize_residue_interface(interaction.get_structure1().get_residues(), real_interactions[0])
            rand_interface_L = self.randomize_residue_interface(interaction.get_structure2().get_residues(), real_interactions[1])
            rand_energies = self.calculate_energies(zip(rand_interface_R, rand_interface_L, real_interactions[2], real_interactions[0], real_interactions[1]), randomized=True)
            for potential_type in rand_energies:
                rand_energies_lists.setdefault(potential_type, []).append(rand_energies[potential_type])
        # Compute Zscores
        real_Zenergies = self.compute_Zscores(real_energies, rand_energies_lists)
        interaction_info.append(real_Zenergies)
        if funnel:
            native_Zenergies = self.compute_Zscores(native_energies, rand_energies_lists)
            partial_Zenergies = self.compute_Zscores(partial_energies, rand_energies_lists)
            external_Zenergies = self.compute_Zscores(external_energies, rand_energies_lists)
            interaction_info.extend([native_Zenergies, partial_Zenergies, external_Zenergies])
        return interaction_info

    def calculate_residue_energies(self, interaction, struc, Zscores=False, randoms=100, alanine_scanning=False, funnel=False):
        '''
        Compute the split potential energies for each interacting residue in a protein-protein interaction
        if Zscores: Computes also the Zscores
        if alanine_scanning: Computes also the alanine mutation of each residue
        if funnel: now the residue energy is the mean of the interaction with the whole protein (random noninterating exposed residues)
        '''
        struc_accepted = ('R','L') # Receptor or Ligand residue energies
        if struc.upper() not in struc_accepted:
            sys.stderr.write('WARNING: Incorrect structure type %s, using R (receptor)!' % struc)
        # Compute energies for each residue
        real_interactions = interaction.get_interacting_residues(c_type=self.c_type, max_distance=self.cutoff+1, uniq=False)
        if struc == 'L':
            interacting_residues = list(set(real_interactions[1]))
            structureA = interaction.get_structure2()
            structureB = interaction.get_structure1()
            noninterating_exposedB = self.__get_noninteracting_residues(list(set(real_interactions[0])), structureB)
        else:
            interacting_residues = list(set(real_interactions[0]))
            structureA = interaction.get_structure1()
            structureB = interaction.get_structure2()
            noninterating_exposedB = self.__get_noninteracting_residues(list(set(real_interactions[1])), structureB)
        residues_interactions = []
        for interacting_residue in interacting_residues:
            residue_interactions = []
            for i in xrange(len(real_interactions[2])):
                if struc == 'L':
                    if real_interactions[1][i] == interacting_residue:
                        residue_interactions.append((real_interactions[1][i], real_interactions[0][i], real_interactions[2][i]))
                else:
                    if real_interactions[0][i] == interacting_residue:
                        residue_interactions.append((real_interactions[0][i], real_interactions[1][i], real_interactions[2][i]))
            residue_energies = self.calculate_energies(residue_interactions, local_type='R')
            residue_info = [interacting_residue, residue_interactions, residue_energies]
            if alanine_scanning:
                alanine_energies = self.calculate_energies(residue_interactions, local_type='R', alanine_mutation=True)
                residue_info.append(alanine_energies)
            if funnel:
                rand_funn_energies_lists = {}
                funnel_energies = {}
                if alanine_scanning:
                    rand_fala_energies_lists = {}
                    funnel_ala_energies = {}
                for random in xrange(randoms):
                    interfaceA, interfaceB, distances = zip(*residue_interactions)
                    rand_interfaceB = self.randomize_residue_interface(noninterating_exposedB, interfaceB)
                    rand_funn_energies = self.calculate_energies(zip(interfaceA, rand_interfaceB, distances, interfaceA, interfaceB), local_type='R', randomized=True)
                    for potential_type in rand_funn_energies:
                        rand_funn_energies_lists.setdefault(potential_type, []).append(rand_funn_energies[potential_type])  
                    if alanine_scanning:
                        rand_fala_energies = self.calculate_energies(zip(interfaceA, rand_interfaceB, distances, interfaceA, interfaceB), local_type='R', randomized=True, alanine_mutation=True)
                        for potential_type in rand_fala_energies:
                            rand_fala_energies_lists.setdefault(potential_type, []).append(rand_fala_energies[potential_type]) 
                for potential in rand_funn_energies_lists:
                    funnel_energies[potential] = numpy.mean(rand_funn_energies_lists[potential])
                residue_info.append(funnel_energies)
                if alanine_scanning:
                    for potential in rand_fala_energies_lists:
                        funnel_ala_energies[potential] = numpy.mean(rand_fala_energies_lists[potential])                  
                    residue_info.append(funnel_ala_energies)
            residues_interactions.append(residue_info)
        if not Zscores:
            return residues_interactions
        # Generate Zscore energies for each residues
        residues_energies = []
        for residue in residues_interactions:
            rand_energies_lists = {}
            for random in xrange(randoms):
                interfaceA, interfaceB, distances = zip(*residue[1])
                rand_interfaceA = self.randomize_residue_interface(structureA.get_residues(), interfaceA)
                rand_interfaceB = self.randomize_residue_interface(structureB.get_residues(), interfaceB)
                rand_energies = self.calculate_energies(zip(rand_interfaceA, rand_interfaceB, distances, interfaceA, interfaceB), local_type='R', randomized=True)
                for potential_type in rand_energies:
                    rand_energies_lists.setdefault(potential_type, []).append(rand_energies[potential_type])
            residue_info = [residue[0], residue[2], self.compute_Zscores(residue[2], rand_energies_lists)]
            info_index = 3
            if alanine_scanning:
                residue_info.extend([residue[info_index], self.compute_Zscores(residue[info_index], rand_energies_lists)])
                info_index += 1
            if funnel:
                residue_info.extend([residue[info_index], self.compute_Zscores(residue[info_index], rand_energies_lists)])
                info_index += 1
            if alanine_scanning and funnel:
                residue_info.extend([residue[info_index], self.compute_Zscores(residue[info_index], rand_energies_lists)])
            residue_info.append(len(residue[1]))
            residues_energies.append(residue_info)
        return residues_energies

    def calculate_residue_energies_between_pairs(self, interaction, struc, Zscores=False, randoms=100, alanine_scanning=False, funnel=False):
        '''
        Compute the split potential energies for each interacting pair of residues in a protein-protein interaction
        if Zscores: Computes also the Zscores
        if alanine_scanning: Computes also the alanine mutation of each residue
        if funnel: now the residue energy is the mean of the interaction with the whole protein (random noninterating exposed residues)
        '''
        struc_accepted = ('R','L') # Receptor or Ligand residue energies
        if struc.upper() not in struc_accepted:
            sys.stderr.write('WARNING: Incorrect structure type %s, using R (receptor)!' % struc)
        # Compute energies for each residue
        real_interactions = interaction.get_interacting_residues(c_type=self.c_type, max_distance=self.cutoff+1, uniq=False)
        if struc == 'L':
            interacting_residues = list(set(real_interactions[1]))
            structureA = interaction.get_structure2()
            structureB = interaction.get_structure1()
            noninterating_exposedB = self.__get_noninteracting_residues(list(set(real_interactions[0])), structureB)
        else:
            interacting_residues = list(set(real_interactions[0]))
            structureA = interaction.get_structure1()
            structureB = interaction.get_structure2()
            noninterating_exposedB = self.__get_noninteracting_residues(list(set(real_interactions[1])), structureB)
        residues_interactions = []
        for interacting_residue in interacting_residues:
            residue_interactions = []
            residue_info = []
            for i in xrange(len(real_interactions[2])):
                if struc == 'L':
                    if real_interactions[1][i] == interacting_residue:
                        res1 = real_interactions[1][i]
                        res2 = real_interactions[0][i]
                        distance = real_interactions[2][i]
                        residue_interactions.append((res1, res2, distance))
                        pair_energies = self.calculate_energies([[res1, res2, distance]], local_type='R')
                        residue_info = [interacting_residue, [(res1, res2, distance)], pair_energies]
                    else:
                        continue
                else:
                    if real_interactions[0][i] == interacting_residue:
                        res1 = real_interactions[0][i]
                        res2 = real_interactions[1][i]
                        distance = real_interactions[2][i]
                        residue_interactions.append((res1, res2, distance))
                        pair_energies = self.calculate_energies([[res1, res2, distance]], local_type='R')
                        residue_info = [interacting_residue, [(res1, res2, distance)], pair_energies]
                    else:
                        continue
                if alanine_scanning:
                    alanine_energies = self.calculate_energies(residue_interactions, local_type='R', alanine_mutation=True)
                    residue_info.append(alanine_energies)
                if funnel:
                    rand_funn_energies_lists = {}
                    funnel_energies = {}
                    if alanine_scanning:
                        rand_fala_energies_lists = {}
                        funnel_ala_energies = {}
                    for random in xrange(randoms):
                        interfaceA, interfaceB, distances = zip(*residue_interactions)
                        rand_interfaceB = self.randomize_residue_interface(noninterating_exposedB, interfaceB)
                        rand_funn_energies = self.calculate_energies(zip(interfaceA, rand_interfaceB, distances, interfaceA, interfaceB), local_type='R', randomized=True)
                        for potential_type in rand_funn_energies:
                            rand_funn_energies_lists.setdefault(potential_type, []).append(rand_funn_energies[potential_type])  
                        if alanine_scanning:
                            rand_fala_energies = self.calculate_energies(zip(interfaceA, rand_interfaceB, distances, interfaceA, interfaceB), local_type='R', randomized=True, alanine_mutation=True)
                            for potential_type in rand_fala_energies:
                                rand_fala_energies_lists.setdefault(potential_type, []).append(rand_fala_energies[potential_type]) 
                    for potential in rand_funn_energies_lists:
                        funnel_energies[potential] = numpy.mean(rand_funn_energies_lists[potential])
                    residue_info.append(funnel_energies)
                    if alanine_scanning:
                        for potential in rand_fala_energies_lists:
                            funnel_ala_energies[potential] = numpy.mean(rand_fala_energies_lists[potential])                  
                        residue_info.append(funnel_ala_energies)
                residues_interactions.append(residue_info)
        if not Zscores:
            return residues_interactions
        # Generate Zscore energies for each residues
        residues_energies = []
        for residue in residues_interactions:
            rand_energies_lists = {}
            for random in xrange(randoms):
                interfaceA, interfaceB, distances = zip(*residue[1])
                rand_interfaceA = self.randomize_residue_interface(structureA.get_residues(), interfaceA)
                rand_interfaceB = self.randomize_residue_interface(structureB.get_residues(), interfaceB)
                rand_energies = self.calculate_energies(zip(rand_interfaceA, rand_interfaceB, distances, interfaceA, interfaceB), local_type='R', randomized=True)
                for potential_type in rand_energies:
                    rand_energies_lists.setdefault(potential_type, []).append(rand_energies[potential_type])
            residue_info = [residue[0], residue[2], self.compute_Zscores(residue[2], rand_energies_lists), residue[1][0][1]]
            info_index = 3
            if alanine_scanning:
                residue_info.extend([residue[info_index], self.compute_Zscores(residue[info_index], rand_energies_lists)])
                info_index += 1
            if funnel:
                residue_info.extend([residue[info_index], self.compute_Zscores(residue[info_index], rand_energies_lists)])
                info_index += 1
            if alanine_scanning and funnel:
                residue_info.extend([residue[info_index], self.compute_Zscores(residue[info_index], rand_energies_lists)])
            residue_info.append(len(residue[1]))
            residues_energies.append(residue_info)
        return residues_energies
    ##################
    # PRIVATE METHODS #
    ##################

    def __get_noninteracting_residues(self, int_residues, structure, exposed=True):
        '''
        Return all exposed non-interacting residues
        '''
        out_residues = []
        for residue in structure:
            if exposed:
                try:
                    residue_exposed = residue.is_exposed()
                except:
                    continue
                if not residue_exposed:
                    continue
            if residue not in int_residues:
                out_residues.append(residue)
        return out_residues

    ########################################################
    # METHODS TO COMPUTE SPLIT POTENTIALS FOR INTERACTIONS #
    ########################################################

    def calculate_energies(self, res_interactions, local_type='G', randomized=False, alanine_mutation=False):
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
            distance = res_interaction[2]
            try:
                triad1 = res_interaction[0].get_triad()
                triad2 = res_interaction[1].get_triad()
                if randomized:
                    triad1 = res_interaction[0].get_triad()[:4] + res_interaction[3].get_triad()[4:]
                    triad2 = res_interaction[1].get_triad()[:4] + res_interaction[4].get_triad()[4:]
                if alanine_mutation: # Generally with local_type R
                    triad1 = 'A-n-' + res_interaction[0].get_triad()[4:]
            except ResidueTriadError as e:
                #sys.stdout.write(str(e))
                #sys.stdout.flush()
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

    def randomize_residue_interface(self, struct_residues, res_interface):

        rand_interface = []
        done = {}
        for res in res_interface:
            if done.has_key(res):
                rand_res = done[res]
            else:
                rand_res = random.choice(struct_residues)
                done[res] = rand_res
            rand_interface.append(rand_res)
        return rand_interface

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
        