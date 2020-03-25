import os, sys, subprocess
import BioLib.Structure.PDB as PDB

class ZDock(object):
    '''
    Class to manage a ZDock transform file
    '''
    def __init__(self, ZDock_file):
        '''
        Call the parsing method
        '''
        self.decoys = []
        self.__parse(ZDock_file)

    def __parse(self, ZDock_file):
        '''
        '''
        ZDock_fo = open(ZDock_file, 'r')
        self.n, self.spacing, self.switch_num = ZDock_fo.readline().strip('\n').split('\t')
        self.rec_rand1 = 0.0
        self.rec_rand2 = 0.0
        self.rec_rand3 = 0.0
        if self.switch_num != '':
            self.rec_rand1, self.rec_rand2, self.rec_rand3 = ZDock_fo.readline().strip('\n').split('\t')
        self.lig_rand1, self.lig_rand2, self.lig_rand3 = ZDock_fo.readline().strip('\n').split('\t')
        self.rec, self.r1, self.r2, self.r3 = ZDock_fo.readline().strip('\n').split('\t')
        self.lig, self.l1, self.l2, self.l3 = ZDock_fo.readline().strip('\n').split('\t')
        if self.switch_num == '1':
            temp_name = self.rec;
            self.rec = self.lig;
            self.lig = temp_name;
        num = 1
        for line in ZDock_fo:
            decoy = line.strip('\n').split('\t')
            self.decoys.append(ZDecoy(num, self, float(decoy[0]), float(decoy[1]), float(decoy[2]), int(decoy[3]), int(decoy[4]), int(decoy[5]), float(decoy[6])))
            num += 1
        ZDock_fo.close()

    # GETTERS #

    def get_decoys(self):
        '''
        Return all decoys
        '''
        return self.decoys

    def get_decoy(self, num):
        '''
        Return a decoy by its number
        '''
        for decoy in self:
            if decoy.get_num() == num:
                return decoy
        return None

    # METHODS #

    def get_num_decoys(self):
        return len(self.decoys)

    def print_structures(self, static_structure, mobile_structure, pdb_file, single_pdb=False):
        '''
        Print all PDB decoys
        '''
        decoy_list = []
        for decoy in self:
            if singlePDB:
                decoy_structure = decoy.get_structure(mobile_structure)
                decoy_list.append((static_structure, decoy_structure))
            else:
                decoy.print_structure(static_structure, decoy_structure, pdb_file+'_'+str(decoy.get_num())+'.pdb')
        if singlePDB:
            PDB.write_pdb(decoy_list, pdb_file+'.pdb', multi_chain=True, multi_model=True)

    # MAGIC METHODS #

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

class ZDecoy(object):
    '''
    Represents a decoy result from a ZDock transform file
    '''
    def __init__(self, num, zdock, alpha, beta, gamma, x, y, z, energy):
        '''
        '''
        self.num = num
        self.zdock = zdock
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma
        self.x = x
        self.y = y
        self.z = z
        self.energy = energy

    # GETTERS #

    def get_num(self):
        return self.num

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
        ZDock_create = os.path.abspath(os.path.join(os.path.dirname(__file__), 'ZDock_create'))
        if self.zdock.switch_num != '':
            create_cmd = '%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s' % (
                          ZDock_create, self.zdock.switch_num,
                          self.zdock.rec_rand1, self.zdock.rec_rand2, self.zdock.rec_rand3,
                          self.zdock.lig_rand1, self.zdock.lig_rand2, self.zdock.lig_rand3, 
                          self.zdock.r1, self.zdock.r2, self.zdock.r3,
                          self.zdock.l1, self.zdock.l2, self.zdock.l3,
                          self.alpha, self.beta, self.gamma,
                          self.x, self.y, self.z,
                          self.zdock.n, self.zdock.spacing)
        else:
            create_cmd = '%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s' % (
                          ZDock_create, 
                          self.zdock.lig_rand1, self.zdock.lig_rand2, self.zdock.lig_rand3, 
                          self.zdock.r1, self.zdock.r2, self.zdock.r3,
                          self.zdock.l1, self.zdock.l2, self.zdock.l3,
                          self.alpha, self.beta, self.gamma,
                          self.x, self.y, self.z,
                          self.zdock.n, self.zdock.spacing)
        
        p = subprocess.Popen(create_cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        create_out, create_err = p.communicate(str(mobile_structure))

        structure = PDB.read_pdb(create_out, from_string=True)
        return structure

    def print_structure(self, static_structure, mobile_structure, pdb_file):
        decoy_structure = self.get_structure(mobile_structure)
        PDB.write_pdb((static_structure, decoy_structure), pdb_file+'.pdb', multi_chain=True, multi_model=False)

    def __str__(self):
        return '%d %f %d %d %d %f %f %f' % (self.num, self.energy, self.x, self.y, self.z, self.alpha, self.beta, self.gamma)
