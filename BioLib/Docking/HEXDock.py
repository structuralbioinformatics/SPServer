import re, copy
import BioLib.Structure.PDB as PDB

class HEXDock():
    '''
    Class to manage the HEX Docking transform file
    '''
    def __init__(self, HEXdockFile):
        '''
        Call the parsing method
        '''
        self.decoys = []
        self.__parse(HEXdockFile)

    def __parse(self, HEXdockFile):
        '''
        '''
        hexdockfo = open(HEXdockFile, 'r')
        for line in hexdockfo:
            if line[0] != '#' and len(line)>1:
                decoy = re.split('\s+', line.strip('\n'))
                self.decoys.append(HEXDecoy(int(decoy[0]), int(decoy[1]), float(decoy[2]), float(decoy[3]), float(decoy[4]), float(decoy[5]), float(decoy[6]), float(decoy[7]), float(decoy[8]), float(decoy[9])))
        hexdockfo.close()
        
    # GETTERS #

    def get_decoys(self):
        return self.decoys

    def get_decoy_by_id(self, cluster):
        for decoy in self.decoys:
            if decoy.get_cluster() == cluster:
                return decoy
        return None

    # METHODS #

    def get_num_decoys(self):
        return len(self.get_decoys())

    def get_structures(self, mobile_structure):
        '''
        ATTENTION: Memory inefficient
        '''
        structure_list = []
        for decoy in self.get_decoys():
            structure_list.append(decoy.get_structure(mobile_structure))
        return structure_list

    def filter_by_spatial_restrictions(self, static_structure, mobile_structure, spatial_restrictions):
        '''
        Filter the decoys by spatial restrictions
        '''
        filtered_decoys = []
        for decoy in self.get_decoys():
            compatible = True
            for spatial_restriction in spatial_restrictions:
                decoy_structure = decoy.get_structure(mobile_structure)
                distance = static_structure.get_residue_by_num(spatial_restriction[0]).get_residue_distance(decoy_structure.get_residue_by_num(spatial_restriction[1]), 'ca')
                if spatial_restriction[3] == 'max':
                    if distance > spatial_restriction[2]:
                        compatible = False
                if spatial_restriction[3] == 'min':
                    if distance < spatial_restriction[2]:
                        compatible = False
            if compatible:
                filtered_decoys.append(decoy)
                print "decoy %d appended" % decoy.get_cluster()
            else:
                print "decoy %d skipped" % decoy.get_cluster()
        self.decoys = filtered_decoys

    def print_structures(self, static_structure, mobile_structure, pdbFile, singlePDB=False):
        '''
        Print all PDB decoys
        '''
        decoyList = []
        for decoy in self.get_decoys():
            if singlePDB:
                decoy_structure = decoy.get_structure(mobile_structure)
                decoyList.append((static_structure, decoy_structure))
            else:
                decoy.print_structure(static_structure, decoy_structure, pdbFile+'_'+str(decoy.get_cluster())+'.pdb')
        if singlePDB:
            PDB.write_pdb(decoyList, pdbFile+'.pdb', multi_chain=True, multi_model=True)

    def __iter__(self):
        '''
        Iterates over each decoy
        '''
        for decoy in self.decoys:
            yield decoy

    def __str__(self):
        '''
        toString method
        '''
        return '\n'.join([decoy.__str__() for decoy in self.decoys])

class HEXDecoy():
    '''
    Object that stores a decoy result from an HEX docking
    '''

    def __init__(self, cluster, solution, RMS, energy, x, y, z, alpha, beta, gamma):
        '''
        '''
        self.cluster = cluster
        self.solution = solution
        self.RMS = RMS
        self.energy = energy
        self.x = x
        self.y = y
        self.z = z
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma

    # GETTERS #

    def get_cluster(self):
        return self.cluster

    def get_solution(self):
        return self.solution

    def get_RMS(self):
        return self.RMS

    def get_energy(self):
        return self.energy

    def get_tranlation_vector(self):
        return (self.x, self.y, self.z)
   
    def get_rotation_vector(self):
        return (self.alpha, self.beta, self.gamma)

    def get_structure(self, mobile_structure):
        '''
        Gets the structure of the given docking result
        '''
        structure = copy.deepcopy(mobile_structure)
        structure.relocate(self.get_tranlation_vector(), self.get_rotation_vector()) 
        return structure

    def print_structure(self, static_structure, mobile_structure, pdbFile):
        decoy_structure = self.get_structure(mobile_structure)
        PDB.write_pdb((static_structure, decoy_structure), pdbFile+'.pdb', multi_chain=True, multi_model=False)

    def __str__(self):
        return '%d %d %f %f %f %f %f %f %f %f' % (self.cluster, self.solution, self.RMS, self.energy, self.x, self.y, self.z, self.alpha, self.beta, self.gamma)
