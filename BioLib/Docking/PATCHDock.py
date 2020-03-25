import re, copy
import BioLib.Structure.PDB as PDB

class PATCHDock(object):
    '''
    Class to manage the PACHDock transform file
    '''
    def __init__(self, PACHDock_file):
        '''
        Call the parsing method
        '''
        self.decoys = []
        self.__parse(PACHDock_file)

    def __parse(self, PACHDock_file):
        '''
        '''
        PACHDock_fo = open(PACHDock_file, 'r')
        for decoy in PACHDock_fo:
            if '#' not in decoy and '||' in decoy:
                decoy = decoy.strip('\n').split('||')
                num = int(decoy[0].split('|')[0])
                score = int(decoy[0].split('|')[1])
                pen = float(decoy[0].split('|')[2])
                area = float(decoy[0].split('|')[3])
                as12 = int(decoy[0].split('|')[6])
                ace = float(decoy[0].split('|')[7])
                vectors = re.split('\s+', decoy[1])
                self.decoys.append(PATCHDecoy(num, score, pen, area, as12, ace,
                                              float(vectors[1]), float(vectors[2]), float(vectors[3]), 
                                              float(vectors[4]), float(vectors[5]), float(vectors[6])))
        PACHDock_fo.close()

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

class PATCHDecoy(object):
    '''
    Represents a decoy result from a PACHDock transform file
    '''
    def __init__(self, num, score, pen, area, as12, ace, alpha, beta, gamma, x, y, z):
        '''
        '''
        self.num = num
        self.score = score
        self.pen = pen
        self.area = area
        self.as12 = as12
        self.ace = ace
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma
        self.x = x
        self.y = y
        self.z = z

    # GETTERS #

    def get_num(self):
        return self.num

    def get_score(self):
        return self.score

    def get_pen(self):
        return self.pen

    def get_area(self):
        return self.area

    def get_as12(self):
        return self.as12

    def get_ace(self):
        return self.ace

    def get_tranlation_vector(self):
        return (self.x, self.y, self.z)
   
    def get_rotation_vector(self):
        return (self.alpha, self.beta, self.gamma)

    def get_structure(self, mobile_structure):
        '''
        Gets the structure of the given docking result
        '''
        structure = copy.deepcopy(mobile_structure)
        structure.relocate(self.get_tranlation_vector(), self.get_rotation_vector(), docking='PATCHDOCK') 
        return structure

    def print_structure(self, static_structure, mobile_structure, pdb_file):
        decoy_structure = self.get_structure(mobile_structure)
        PDB.write_pdb((static_structure, decoy_structure), pdb_file+'.pdb', multi_chain=True, multi_model=False)

    def __str__(self):
        return '%d %d %f %f %f %f %f %f' % (self.num, self.score, self.alpha, self.beta, self.gamma, self.x, self.y, self.z)
