import os, sys
class ResidueDistanceError(Exception):
    '''
    The distance between two residues cannot be computed because CA or CB atoms are not found
    '''
    def __init__(self, c_type, residue1, residue2):
        '''
        Constructor
        '''
        self.c_type = c_type
        self.residue1 = residue1
        self.residue2 = residue2

    def __str__(self):
        return 'EXCEPTION: Cannot compute %s distance between residue %s and %s!\n' % (self.c_type.upper(), self.residue1, self.residue2)

class ResidueTypeShortError(Exception):
    '''
    Cannot converts the 3 letter type into 1 letter type (example: 'ALA' -> 'A')
    '''
    def __init__(self, residue_type):
        '''
        Constructor
        '''
        self.residue_type = residue_type

    def __str__(self):
        return '%s residue type has not a valid short translation\n' % self.residue_type

class ResidueExpositionError(Exception):
    '''
    The exposition of a residue cannot be computed
    '''
    def __init__(self, residue_num, reason):
        '''
        Constructor
        '''
        self.residue_num = residue_num
        self.reason = reason

    def __str__(self):
        return 'Cannot compute residue exposition for residue %s: %s' % (self.residue_num, self.reason)

class ResidueChargeError(Exception):
    '''
    The charge of a residue cannot be computed
    '''
    def __init__(self, residue_num):
        '''
        Constructor
        '''
        self.residue_num = residue_num

    def __str__(self):
        return 'Cannot compute charge for residue %s\n' % self.residue_num

class ResiduePolarityError(Exception):
    '''
    The Polarity of a residue cannot be computed
    '''
    def __init__(self, residue_num):
        '''
        Constructor
        '''
        self.residue_num = residue_num

    def __str__(self):
        return 'Cannot get polarity for residue %s\n' % self.residue_num

class ResidueTriadError(Exception):
    '''
    Cannot compute the triad for a residue
    '''
    def __init__(self, residue_num, reason):
        '''
        Constructor
        '''
        self.residue_num = residue_num
        self.reason = reason

    def __str__(self):
        return 'Cannot compute triad for residue %s: %s' % (self.residue_num, self.reason) 

class RelocateProgramError(Exception):
    '''
    Incorrect docking program used to relocate structures
    '''
    def __init__(self, docking):
        '''
        Constructor
        '''
        self.docking = docking

    def __str__(self):
        return 'Cannot relocate: Incorrect docking program %s\n' % self.docking