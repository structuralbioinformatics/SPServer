import os, sys, random
import numpy as np

class IC_scores(object):
    '''
    Compute the IC (inter-residue contacs) as described in Bonvin eLife 2015
    '''
    def __init__(self, interaction):
        '''
        Contructor
        '''
        self.interaction = interaction
        interacting_residues = self.interaction.get_interacting_residues(c_type = 'MIN', max_distance = 5.5, uniq=False)
        self.contact_pos1 = interacting_residues[0]
        self.contact_pos2 = interacting_residues[1]
        self.contact_distances = interacting_residues[2]
        self.res_types = {'HIS' : 'charged',
                          'ARG' : 'charged',
                          'LYS' : 'charged',
                          'ASP' : 'charged',
                          'GLU' : 'charged',
                          'SER' : 'polar',
                          'THR' : 'polar',
                          'ASN' : 'polar',
                          'GLN' : 'polar',
                          'CYS' : 'apolar',
                          'GLY' : 'apolar',
                          'PRO' : 'apolar',
                          'ALA' : 'apolar',
                          'VAL' : 'apolar',
                          'ILE' : 'apolar',
                          'MET' : 'apolar',
                          'LEU' : 'apolar',
                          'PHE' : 'apolar',
                          'TYR' : 'apolar',
                          'TRP' : 'apolar'  }

    def calculate(self):
        '''
        Calculate the IC scores
        '''
        ic_scores =  {'Total'           : 0,
                      'charged-charged' : 0, 
                      'charged-polar'   : 0,
                      'charged-apolar'  : 0,
                      'polar-polar'     : 0,
                      'polar-apolar'    : 0,
                      'apolar-apolar'   : 0 }
        for x in xrange(len(self.contact_distances)):
            res_type1 = self.contact_pos1[x].get_type()
            res_type2 = self.contact_pos2[x].get_type()
            if ic_scores.has_key(self.res_types[res_type1]+'-'+self.res_types[res_type2]):
                ic_scores[self.res_types[res_type1]+'-'+self.res_types[res_type2]] += 1
                ic_scores['Total'] += 1
            elif ic_scores.has_key(self.res_types[res_type2]+'-'+self.res_types[res_type1]):
                ic_scores[self.res_types[res_type2]+'-'+self.res_types[res_type1]] += 1
                ic_scores['Total'] += 1
        return ic_scores

class Elec_potential(object):
    '''
    Compute the electrostatic energy using coulomb's law (q1*q2/r**2), asuming each amino acid
    has the same charge.
    The more negative, more energy.
    '''
    def __init__(self, interaction):
        '''
        Contructor
        '''
        self.interaction = interaction
        interacting_residues = self.interaction.get_interacting_residues(c_type = 'CB', max_distance = 40, uniq=False)
        self.contact_pos1 = interacting_residues[0]
        self.contact_pos2 = interacting_residues[1]
        self.contact_distances = interacting_residues[2]

    def calculate(self):
        '''
        Calculate the absolute electrostatic energy
        '''
        Elec_energy = 0
        for x in xrange(len(self.contact_distances)):
            q1 = self.contact_pos1[x].get_charge()
            q2 = self.contact_pos2[x].get_charge()
            r2  = self.contact_distances[x]**2
            Elec_energy += (q1*q2)/r2
        return Elec_energy

    def calculate_normalized(self, randoms = 100):
        '''
        Calculate the normalized electrostatic energy
        '''
        Elec_distribution = []
        for x in xrange(randoms):
            Elec_energy = 0
            for x in xrange(len(self.contact_distances)):
                random_pos1 = random.choice(self.interaction.get_structure1().get_residues())
                random_pos2 = random.choice(self.interaction.get_structure2().get_residues())
                q1 = random_pos1.get_charge()
                q2 = random_pos2.get_charge()
                r2  = self.contact_distances[x]**2
                Elec_energy += (q1*q2)/r2
            Elec_distribution.append(Elec_energy)
        ZElec_energy = (self.calculate()-np.mean(Elec_distribution))/np.std(Elec_distribution)
        return ZElec_energy

class Polarity_potential(object):
    '''
    Compute a a/polar score based on the polarity of the contact residues.
    '''
    def __init__(self, interaction):
        '''
        Contructor
        '''
        self.interaction = interaction
        interacting_residues = self.interaction.get_interacting_residues(c_type = 'CB', max_distance = 12, uniq=False)
        self.contact_pos1 = interacting_residues[0]
        self.contact_pos2 = interacting_residues[1]
        self.contact_distances = interacting_residues[2]

    def calculate(self):
        '''
        Calculate the absolute potential
        '''
        Polar_score = 0
        Apolar_score = 0
        for x in xrange(len(self.contact_distances)):
            if self.contact_pos1[x].is_polar() and self.contact_pos2[x].is_polar():
                Polar_score += 1
            if not self.contact_pos1[x].is_polar() and not self.contact_pos2[x].is_polar():
                Apolar_score += 1
        return Polar_score, Apolar_score

    def calculate_normalized(self, randoms = 100):
        '''
        Calculate the normalized potential
        '''
        Polar_distribution = []
        Apolar_distribution = []
        for x in xrange(randoms):
            Polar_score = 0
            Apolar_score = 0
            for x in xrange(len(self.contact_distances)):
                random_pos1 = random.choice(self.interaction.get_structure1().get_residues())
                random_pos2 = random.choice(self.interaction.get_structure2().get_residues())
                if self.random_pos1.is_polar() and self.random_pos2.is_polar():
                    Polar_score += 1
                if not self.random_pos1.is_polar() and not self.random_pos2.is_polar():
                    Apolar_score += 1
            Polar_distribution.append(Polar_score)
            Apolar_distribution.append(Apolar_score)
        ZPolar_score = (self.calculate()[0]-np.mean(Polar_distribution))/np.std(Polar_distribution)
        ZApolar_score = (self.calculate()[1]-np.mean(Apolar_distribution))/np.std(Apolar_distribution)
        return ZPolar_score, ZApolar_score
