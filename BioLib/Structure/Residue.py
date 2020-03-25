import sys
from BioLib.Tools.BioExceptions import *

class Residue(object):
    '''
    Represents an structure residue
    '''
    aa_dic   = {'ALA':"A", 'ARG':"R", 'ASN':"N", 'ASP':"D", 'ASX':"B", 'CYS':"C", 'GLN':"Q", 'GLU':"E",
                'GLX':"Z", 'GLY':"G", 'HIS':"H", 'ILE':"I", 'LEU':"L", 'LYS':"K", 'MET':"M", 'PHE':"F",
                'PRO':"P", 'SER':"S", 'THR':"T", 'TRP':"W", 'TYR':"Y", 'VAL':"V", 'SEC':"U", 'PYL':"O", 
                'XLE':"J", 'HSE':"H", 'HSD':"H", 'MSE':"M", 'DVA':"V" }
    surface  = {'A': 115, 'C': 149, 'D': 170, 'E': 207, 'F': 230, 'G': 86,  'H': 206,
                'I': 187, 'K': 222, 'L': 192, 'M': 210, 'N': 184, 'P': 140, 'Q': 208,
                'R': 263, 'S': 140, 'T': 164, 'V': 161, 'W': 269, 'Y': 257 }
    polarity = {'A': False, 'C': False, 'D': True,  'E': True,  'F': False, 'G': False, 'H': True,
                'I': False, 'K': True,  'L': False, 'M': False, 'N': True,  'P': False, 'Q': True,
                'R': True,  'S': True,  'T': True,  'V': False, 'W': False, 'Y': True }
    charge   = {'A': 0, 'C': 0, 'D':-1, 'E':-1, 'F': 0, 'G': 0, 'H': 0,
                'I': 0, 'K': 1, 'L': 0, 'M': 0, 'N': 0, 'P': 0, 'Q': 0,
                'R': 1, 'S': 0, 'T': 0, 'V': 0, 'W': 0, 'Y': 0 }

    def __init__(self, num, res_type):
        '''
        Constructor
        '''
        self.num = int(num)
        self.type = res_type
        self.atoms = []
        self.ss = None
        self.acc = None
        self.ca = None
        self.cb = None

    # Getters #
    def get_num(self):
        return self.num

    def get_type(self):
        return self.type

    def get_atoms(self):
        return self.atoms

    def get_atom_by_num(self, num):
        for atom in self:
            if atom.get_num() == num:
                return atom

    def get_number_atoms(self):
        return len(self.atoms)

    def get_ss(self):
        return self.ss

    def get_acc(self):
        return self.acc

    def get_ca(self):
        return self.ca

    def get_cb(self, glyCA=True):
        if glyCA and self.is_gly():
            return self.ca
        return self.cb

    def get_charge(self):
        try:
            charge = Residue.charge[self.get_type_short()]
        except ResidueTypeShortError as e:
            raise e
        except KeyError as e:
            raise ResidueChargeError(self.get_num())
        return charge

    def get_triad(self):
        '''
        Needed by SplitPotentials library
        '''
        try:
            # Polarity
            if Residue.polarity[self.get_type_short()]: p = 'p'
            else:                                       p = 'n'
            # Exposition
            if self.is_exposed():                       e = 'E'
            else:                                       e = 'B'
        except ResidueTypeShortError as e:
            raise ResidueTriadError(self.get_num(), str(e))
        except KeyError as e:
            raise ResidueTriadError(self.get_num(), 'Cannot get polarity\n')
        except ResidueExpositionError as e:
            raise ResidueTriadError(self.get_num(), str(e))

        if self.get_ss() == None:
            raise ResidueTriadError(self.get_num(), 'Cannot get secondary structure\n')

        return '%s-%s-%s-%s' % (self.get_type_short(), p, self.get_ss(), e)

    # Setters #
    def set_num(self, num):
        '''
        Modify the residue number
        '''
        self.num = num

    def set_ss(self, ss):
        '''
        Sets secondary structure
        '''
        self.ss = ss

    def set_acc(self, acc):
        '''
        Sets accesible surface area
        '''
        self.acc = acc
   
    # Class methods #
    def add_atom(self, atom):
        self.atoms.append(atom)
        if atom.get_type().startswith('CA'): self.ca = atom
        elif atom.get_type().startswith('CB'): self.cb = atom

    def is_exposed(self, threshold=0.6):
        '''
        Returns if the residue is exposed (more than 60 percent exposed by default)
        Baldo orifginal script: int(10*float(self.acc)/Residue.surface[self.get_type_short()]) > 5
        DSSP must be executed before in order to get the acc (accessibility)
        '''
        try:
            return float(self.acc)/Residue.surface[self.get_type_short()] >= threshold
        except TypeError as e:
            raise ResidueExpositionError(self.get_num(), 'No ACC\n')
        except ResidueTypeShortError as e:
            raise ResidueExpositionError(self.get_num(), str(e))
        except KeyError as e:
            raise ResidueExpositionError(self.get_num(), 'Cannot get surface\n')

    def is_gly(self):
        '''
        Return True if the residue is a Glycine, else return False
        '''
        if self.get_type() == 'GLY':
            return True
        return False

    def is_polar(self):
        '''
        Return True if the residue is Polar, else return False
        '''
        try:
            polar = Residue.polarity[self.get_type_short()]
        except ResidueTypeShortError as e:
            raise e
        except KeyError as e:
            raise ResiduePolarityError(self.get_num())
        return polar

    def has_ca(self):
        '''
        Return True if residue has CA, else return False
        '''
        if self.ca is not None:
            return True
        return False

    def has_cb(self):
        '''
        Return True if residue has CB, else return False
        '''
        if self.cb is not None:
            return True
        return False

    def get_type_short(self):
        '''
        Converts the 3 letter type into 1 letter type (example: 'ALA' -> 'A')
        '''
        try:
            return Residue.aa_dic[self.get_type()]
        except KeyError:
            raise ResidueTypeShortError(self.get_type())
    
    def _get_min_distance(self, residue):
        '''
        Gets the minimum distance between two residues
        '''
        try:
            return min([x.get_distance(y) for x in self for y in residue])
        except:
            raise ResidueDistanceError('MIN', self.get_num(), residue.get_num())

    def _get_ca_distance(self, residue):
        '''
        Gets the distance between the CA of two residues
        '''
        try:
            return self.ca.get_distance(residue.ca)
        except:
            raise ResidueDistanceError('CA', self.get_num(), residue.get_num())

    def _get_cb_distance(self, residue):
        '''
        Gets the distance between the CB of two residues (CA if Glycine)
        '''
        try:
            if self.is_gly() and residue.is_gly():
                return self.ca.get_distance(residue.ca)
            elif not self.is_gly() and residue.is_gly():
                return self.cb.get_distance(residue.ca)
            elif self.is_gly() and not residue.is_gly():
                return self.ca.get_distance(residue.cb)
            else:
                return self.cb.get_distance(residue.cb)
        except:
            raise ResidueDistanceError('CB', self.get_num(), residue.get_num())

    def get_residue_distance(self, residue, c_type):
        '''
        Gets the distance between two residues
        '''
        c_type_accepted = ('CA','CB','MIN')
        if c_type.upper() not in c_type_accepted:
            sys.stderr.write('Residue %d: Incorrect c_type: %s. Using CB...\n' % (self.get_num(), c_type))
        try:
            if c_type.upper() == 'CA':
                distance = self._get_ca_distance(residue)
            elif c_type.upper() == 'MIN':
                distance = self._get_min_distance(residue)
            else:
                distance = self._get_cb_distance(residue)
        except ResidueDistanceError as e:
            raise e
        return distance

    def relocate(self, translation_vector, rotation_vector, docking='HEX'):
        '''
        Relocate the residue position 
        @translation_vector = Vector used to translate the structure
        @Rotation_vector = Vector used to rotate the structure
        @docking = Docking program used to get the vectors
        '''
        for atom in self:
            try:
                atom.relocate(translation_vector, rotation_vector, docking)
            except RelocateProgramError as e:
                raise e

    def __iter__(self):
        '''
        Iterates over each atom
        '''
        for atom in self.atoms:
            yield atom

    def __str__(self):
        '''
        Return the information of the residue
        '''

        residueInfo = ''
        for atom in self:
            residueInfo += 'ATOM%s  %s%s%s%s%s%s\n' % (str(atom.get_num()).rjust(7),
                                                       atom.get_type().ljust(4),
                                                       self.get_type().rjust(3),
                                                       str(self.get_num()).rjust(4),
                                                       ('%.3f' % float(atom.get_x())).rjust(12),
                                                       ('%.3f' % float(atom.get_y())).rjust(8),
                                                       ('%.3f' % float(atom.get_z())).rjust(8))
        return residueInfo
