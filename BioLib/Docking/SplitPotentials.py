import os, sys, random, math
import numpy
from BioLib.Tools.BioExceptions import ResidueTriadError

class SplitPotentialsData(object):
    '''
        Contain the SplitPotential Calculation
    '''
    def __init__(self, ori_ppResidues, ori_combined, ori_ppRE, ori_localppRE, ori_ppEnvironment, ori_ppDist,
                 cof_ppResidues, cof_combined, cof_ppRE, cof_localppRE, cof_ppEnvironment, cof_ppDist,
                 max_ppResidues, max_ppEnvironment, max_ppRE, max_combined,
                 min_ppResidues, min_ppEnvironment, min_ppRE, min_combined,
                 men_ppResidues, men_ppEnvironment, men_ppRE, men_combined):

        self.original    = {"D-Score":ori_ppResidues, "D-Energy":ori_combined, "D-Eaa3Denv":ori_ppRE,
                            "D-Elocal":ori_localppRE, "D-E3Denv":ori_ppEnvironment, "D-E3D":ori_ppDist,
                            "C-Score":cof_ppResidues, "C-Energy":cof_combined, "C-Eaa3Denv":cof_ppRE,
                            "C-Elocal":cof_localppRE, "C-E3Denv":cof_ppEnvironment, "C-E3D":cof_ppDist,
                            "M-Score":men_ppResidues, "M-Energy":men_combined, "M-Eaa3Denv":men_ppRE,
                                                      "M-E3Denv":men_ppEnvironment, 
                            "U-Score":max_ppResidues, "U-Energy":max_combined, "U-Eaa3Denv":max_ppRE,
                                                      "U-E3Denv":max_ppEnvironment, 
                            "L-Score":min_ppResidues, "L-Energy":min_combined, "L-Eaa3Denv":min_ppRE,
                                                      "L-E3Denv":min_ppEnvironment}
        
        self.normalized = {"D-ZScore":0,"D-ZEnergy":0,"D-ZEaa3Denv":0,"D-ZElocal":0,"D-ZE3Denv":0,
                           "C-ZScore":0,"C-ZEnergy":0,"C-ZEaa3Denv":0,"C-ZElocal":0,"C-ZE3Denv":0,
                           "M-ZScore":0,"M-ZEnergy":0,"M-ZEaa3Denv":0,              "M-ZE3Denv":0,
                           "U-ZScore":0,"U-ZEnergy":0,"U-ZEaa3Denv":0,              "U-ZE3Denv":0,
                           "L-ZScore":0,"L-ZEnergy":0,"L-ZEaa3Denv":0,              "L-ZE3Denv":0}

        self.original_key    = ["D-Score", "D-Energy", "D-Eaa3Denv", "D-Elocal", "D-E3Denv", "D-E3D",
                                "C-Score", "C-Energy", "C-Eaa3Denv", "C-Elocal", "C-E3Denv", "C-E3D",
                                "M-Score", "M-Energy", "M-Eaa3Denv",             "M-E3Denv",
                                "U-Score", "U-Energy", "U-Eaa3Denv",             "U-E3Denv",
                                "L-Score", "L-Energy", "L-Eaa3Denv",             "L-E3Denv"]
        self.normalized_key  = ["D-ZScore","D-ZEnergy","D-ZEaa3Denv","D-ZElocal","D-ZE3Denv",
                                "C-ZScore","C-ZEnergy","C-ZEaa3Denv","C-ZElocal","C-ZE3Denv",
                                "M-ZScore","M-ZEnergy","M-ZEaa3Denv",            "M-ZE3Denv",
                                "U-ZScore","U-ZEnergy","U-ZEaa3Denv",            "U-ZE3Denv",
                                "L-ZScore","L-ZEnergy","L-ZEaa3Denv",            "L-ZE3Denv"] 
    """
        Getters
    """
    def get_value_by_label(self, label):
        if self.original.has_key(label):
            return self.original[label]
        
    def get_value_by_label_str(self,label):
        return "%.3e" %self.get_value_by_label(label = label)
    
    def get_tabbed_labels(self):
        return "\t".join(self.original_key)
    
    def get_tabbed_values(self):
        outarray = []
        for label in self.original_key:
            outarray.append(self.get_value_by_label_str(label = label))
        return "\t".join(outarray)
    
    def get_normalized_value_by_label(self, label):
        if self.normalized.has_key(label):
            return self.normalized[label]
        
    def get_normalized_value_by_label_str(self,label):
        return "%.3e" %self.get_normalized_value_by_label(label = label)
    
    def get_normalized_tabbed_labels(self):
        return "\t".join(self.normalized_key)
    
    def get_normalized_tabbed_values(self):
        outarray = []
        for label in self.normalized_key:
            outarray.append(self.get_normalized_value_by_label_str(label = label))
        return "\t".join(outarray)
    
    def get_all_tabbed_labels(self):
        return "\t".join((self.get_tabbed_labels(),self.get_normalized_tabbed_labels()))
    
    def get_all_tabbed_values(self):
        return "\t".join((self.get_tabbed_values(),self.get_normalized_tabbed_values()))
    
    """
        Setters
    """
    def set_random_vector(self, vector, label):    

        if not self.normalized.has_key(label):
            return
        
        mean = numpy.mean(vector)
        std = numpy.std(vector)
        
        original_label = label[:2]+label[3:]
        
        self.normalized[label] = (self.original[original_label] - mean) / std

    """
        toString
    """
    def toString(self, header = True, ID = None, headertag = None):
        
        line = []
        
        if ID == None:    
            if header:
                line.append(self.get_all_tabbed_labels())
            line.append(self.get_all_tabbed_values())
        else:
            if header:
                line.append(str(headertag) + "\t" + self.get_all_tabbed_labels())
            line.append(str(ID) + "\t" + self.get_all_tabbed_values()) 
        
        return "\n".join(line)
    

class SplitPotentialsMatrix(object):
    
    """
        Constructor
    """
    def __init__(self, sp_matrix = 'ppi3_cbR12.out', cutoff = 12):
        
        sp_file = os.path.abspath(os.path.join(os.path.dirname(__file__), sp_matrix))
        if not os.path.isfile(sp_file):
            sys.stderr.write("Can not read the potential file %s\n\n" % sp_file)
            sys.exit()
            
        self.file = sp_file    
        self.cutoff = cutoff
        
        self.ppResidues     = {}        # Pair Potentials of Residues
        self.ppEnvironment  = {}        # Pair Potentials of Environment
        self.localppRE      = {}        # Local Pair Potentials of Residues and Environment
        self.ppRE           = {}        # Pair Potentials of Residues and Environment
        self.ppDist         = {}        # Standard Pair Potentials of Distances
#       self.envWeights     = {}        # Solution for Environmental Weights
#       self.ppAccumPair    = {}        # Pair Potentials of Accumulated Pairs
        
        self.realInteraction = None
        
#       self.stats_ppResidues = {}      # Mean, Max and Min for each pair of Pair Potentials of Residues
#       self.stats_ppEnvironment = {}   # Mean, Max and Min for each pair of Pair Potentials of Environment
#       self.stats_ppRE = {}            # Mean, Max and Min for each pair of Pair Potentials of Residues and Environment
        
        tag = 0
        pp_fd = open(self.file)
        for line in pp_fd:
            if line.startswith("PAIR-POTENTIALS OF RESIDUES AND ENVIRONMENTS"):             tag = 4
            elif line.startswith("PAIR-POTENTIALS OF ENVIRONMENTS"):                        tag = 2
            elif line.startswith("LOCAL PAIR-POTENTIALS OF RESIDUES AND ENVIRONMENTS"):     tag = 3
            elif line.startswith("PAIR-POTENTIALS OF RESIDUES"):                            tag = 1
            elif line.startswith("STANDARD PAIR-POTENTIALS OF DISTANCES"):                  tag = 5
#           elif line.startswith("SOLUTION FOR ENVIRONMENTAL WEIGHTS"):                     tag = 6
#           elif line.startswith("PAIR-POTENTIALS OF ACCUMULATED PAIRS"):                   tag = 7
            elif line.startswith("END"):                                                    tag = 0
            
            self.process_file_line(line = line, tag = tag)

        pp_fd.close()

    def add_interaction(self, interaction):
        self.realInteraction = interaction
        real_interacting_residues = interaction.get_interacting_residues(c_type = 'CB', max_distance = 13, uniq=False)
        self.str1_contact_pos = [x.get_num() for x in real_interacting_residues[0]]
        self.str2_contact_pos = [x.get_num() for x in real_interacting_residues[1]]
        self.contact_distances = real_interacting_residues[2]
        
    def calculateSplitPotentials(self, randoms = 100):

        # Calculate Original Split Potentials
        ori_ppResidues,     cof_ppResidues,     max_ppResidues,     min_ppResidues,        men_ppResidues      = 0, 0, 0, 0 ,0
        ori_ppEnvironment,  cof_ppEnvironment,  max_ppEnvironment,  min_ppEnvironment,     men_ppEnvironment   = 0, 0, 0, 0 ,0
        ori_localppRE,      cof_localppRE       = 0, 0
        ori_ppRE,           cof_ppRE,           max_ppRE,           min_ppRE,               men_ppRE           = 0, 0, 0, 0 ,0
        ori_ppDist,         cof_ppDist          = 0, 0
        ori_combined,       cof_combined,       max_combined,       min_combined,           men_combined        = 0, 0, 0, 0 ,0

        for x in xrange(len(self.str1_contact_pos)):
            residue1 = self.realInteraction.get_structure1().get_residue_by_num(self.str1_contact_pos[x])
            residue2 = self.realInteraction.get_structure2().get_residue_by_num(self.str2_contact_pos[x])
            
            try:
                triad1 = residue1.get_triad()
                triad2 = residue2.get_triad()
            except ResidueTriadError as e:
                sys.stdout.write(str(e))
                sys.stdout.flush()
                continue
            
            positional_ppResidues   = float(self.get_ppResidues(k1=triad1[0:1],k2=triad2[0:1], distance=self.contact_distances[x]))
            ori_ppResidues          += positional_ppResidues
            maximum_ppResidues      = float(self.get_ppResidues_stat(k1=triad1[0:1],k2=triad2[0:1], stat='upper'))
            max_ppResidues          += maximum_ppResidues
            minimum_ppResidues      = float(self.get_ppResidues_stat(k1=triad1[0:1],k2=triad2[0:1], stat='lower'))
            min_ppResidues          += minimum_ppResidues
            mean_ppResidues         = float(self.get_ppResidues_stat(k1=triad1[0:1],k2=triad2[0:1], stat='mean'))
            men_ppResidues          += mean_ppResidues
                    
            positional_ppEnvironment = float(self.get_ppEnvironment(k1 = triad1[2:],k2 = triad2[2:], distance=self.contact_distances[x]))
            ori_ppEnvironment       += positional_ppEnvironment
            maximum_ppEnvironment   = float(self.get_ppEnvironment_stat(k1=triad1[2:],k2=triad2[2:], stat='upper'))
            max_ppEnvironment       += maximum_ppEnvironment
            minimum_ppEnvironment   = float(self.get_ppEnvironment_stat(k1=triad1[2:],k2=triad2[2:], stat='lower'))
            min_ppEnvironment       += minimum_ppEnvironment
            mean_ppEnvironment      = float(self.get_ppEnvironment_stat(k1=triad1[2:],k2=triad2[2:], stat='mean'))
            men_ppEnvironment       += mean_ppEnvironment
            
            positional_localppRE    = float(self.get_localppRE(k=triad1)) + float(self.get_localppRE(k=triad2))
            ori_localppRE           += positional_localppRE
            
            positional_ppRE         = float(self.get_ppRE(k1 = triad1, k2 = triad2, distance=self.contact_distances[x]))
            ori_ppRE                += positional_ppRE
            maximum_ppRE            = float(self.get_ppRE_stat(k1=triad1,k2=triad2, stat='upper'))
            max_ppRE                += maximum_ppRE
            minimum_ppRE            = float(self.get_ppRE_stat(k1=triad1,k2=triad2, stat='lower'))
            min_ppRE                += minimum_ppRE
            mean_ppRE               = float(self.get_ppRE_stat(k1=triad1,k2=triad2, stat='mean'))
            men_ppRE                += mean_ppRE
            
            positional_ppDist       = float(self.get_ppDist(k=self.contact_distances[x]))
            ori_ppDist              += positional_ppDist

            positional_combined     = positional_ppRE + positional_ppDist - positional_localppRE - positional_ppEnvironment
            ori_combined            += positional_combined
            maximum_combined        = maximum_ppRE + positional_ppDist - positional_localppRE - maximum_ppEnvironment
            max_combined            += maximum_combined
            minimum_combined        = minimum_ppRE + positional_ppDist - positional_localppRE - minimum_ppEnvironment
            min_combined            += minimum_combined
            mean_combined           = mean_ppRE + positional_ppDist - positional_localppRE - mean_ppEnvironment
            men_combined            += mean_combined
            
            if math.floor(self.contact_distances[x]) <= self.cutoff:
                cutoff_ppResidues   = float(self.get_ppResidues(k1=triad1[0:1],k2=triad2[0:1], distance=self.cutoff))
                cof_ppResidues      += cutoff_ppResidues
                
                cutoff_ppEnvironment = float(self.get_ppEnvironment(k1 = triad1[2:],k2 = triad2[2:], distance=self.cutoff))
                cof_ppEnvironment   += cutoff_ppEnvironment
                
                cutoff_localppRE    = float(self.get_localppRE(k=triad1)) + float(self.get_localppRE(k=triad2))
                cof_localppRE       += cutoff_localppRE
                
                cutoff_ppRE         = float(self.get_ppRE(k1 = triad1, k2 = triad2, distance=self.cutoff))
                cof_ppRE            += cutoff_ppRE
                
                cutoff_ppDist       = float(self.get_ppDist(k=self.cutoff))
                cof_ppDist          += cutoff_ppDist
    
                cutoff_combined     = cutoff_ppRE + cutoff_ppDist - cutoff_localppRE - cutoff_ppEnvironment
                cof_combined        += cutoff_combined
     
        result = SplitPotentialsData(ori_ppResidues = ori_ppResidues, ori_combined = ori_combined, ori_ppRE = ori_ppRE, 
                                     ori_localppRE = ori_localppRE, ori_ppEnvironment = ori_ppEnvironment, ori_ppDist = ori_ppDist,
                                     cof_ppResidues = cof_ppResidues, cof_combined = cof_combined, cof_ppRE = cof_ppRE,
                                     cof_localppRE = cof_localppRE, cof_ppEnvironment = cof_ppEnvironment, cof_ppDist = cof_ppDist,
                                     max_ppResidues = max_ppResidues, max_ppEnvironment = max_ppEnvironment, max_ppRE = max_ppRE, max_combined = max_combined,
                                     min_ppResidues = min_ppResidues, min_ppEnvironment = min_ppEnvironment, min_ppRE = min_ppRE, min_combined = min_combined,
                                     men_ppResidues = men_ppResidues, men_ppEnvironment = men_ppEnvironment, men_ppRE = men_ppRE, men_combined = men_combined)
        
        # If we are not going to do randoms we can return now
        if randoms == 0:
            #return (ori_ppResidues, ori_combined, ori_ppRE, ori_localppRE, ori_ppEnvironment, ori_ppDist)
            return result
        
        # Calculate randoms: Generate zero positions
        rad_ppResidues      = numpy.zeros((randoms+1), float)
        rdC_ppResidues      = numpy.zeros((randoms+1), float)
        rdM_ppResidues      = numpy.zeros((randoms+1), float)
        rdU_ppResidues      = numpy.zeros((randoms+1), float)
        rdL_ppResidues      = numpy.zeros((randoms+1), float)
        rad_ppEnvironment   = numpy.zeros((randoms+1), float)
        rdC_ppEnvironment   = numpy.zeros((randoms+1), float)
        rdM_ppEnvironment   = numpy.zeros((randoms+1), float)
        rdU_ppEnvironment   = numpy.zeros((randoms+1), float)
        rdL_ppEnvironment   = numpy.zeros((randoms+1), float)
        rad_localppRE       = numpy.zeros((randoms+1), float)
        rdC_localppRE       = numpy.zeros((randoms+1), float)
        rad_ppRE            = numpy.zeros((randoms+1), float)
        rdC_ppRE            = numpy.zeros((randoms+1), float)
        rdM_ppRE            = numpy.zeros((randoms+1), float)
        rdU_ppRE            = numpy.zeros((randoms+1), float)
        rdL_ppRE            = numpy.zeros((randoms+1), float)
        rad_ppDist          = numpy.zeros((randoms+1), float)
        rdC_ppDist          = numpy.zeros((randoms+1), float)
        rad_combined        = numpy.zeros((randoms+1), float)
        rdC_combined        = numpy.zeros((randoms+1), float)
        rdM_combined        = numpy.zeros((randoms+1), float)
        rdU_combined        = numpy.zeros((randoms+1), float)
        rdL_combined        = numpy.zeros((randoms+1), float)
        
        for rad in xrange(randoms):
            # Create random selection of positions for str1 and str2
            # We pick positions by array position as we do not care what are we picking and will avoid
            #    errors in case of gaps in the structure.
            # We need to take into account that we can actually have more interactions than Aa in the 
            #    protein, as the same Aa can have several interactions
            rad_str1 = self.randomise_interaface(list=self.str1_contact_pos, 
                                                 range=len(self.realInteraction.get_structure1().get_residues()))
            rad_str2 = self.randomise_interaface(list=self.str2_contact_pos, 
                                                 range=len(self.realInteraction.get_structure2().get_residues()))
            for x in xrange(len(rad_str1)):
                residue1 = self.realInteraction.get_structure1().get_residues()[rad_str1[x]]
                residue2 = self.realInteraction.get_structure2().get_residues()[rad_str2[x]]
                
                # We keep in the triad the exposition and secondary structure
                try:
                    triad1 = residue1.get_triad()[:4] + self.realInteraction.get_structure1().get_residue_by_num(self.str1_contact_pos[x]).get_triad()[4:]
                    triad2 = residue2.get_triad()[:4] + self.realInteraction.get_structure2().get_residue_by_num(self.str2_contact_pos[x]).get_triad()[4:]
                except ResidueTriadError as e:
                    sys.stdout.write(str(e))
                    sys.stdout.flush()
                    continue

                rad_ppResidues[rad]     += float(self.get_ppResidues(k1=triad1[0:1],k2=triad2[0:1],distance=self.contact_distances[x]))
                rdM_ppResidues[rad]     += float(self.get_ppResidues_stat(k1=triad1[0:1],k2=triad2[0:1], stat='mean'))
                rdU_ppResidues[rad]     += float(self.get_ppResidues_stat(k1=triad1[0:1],k2=triad2[0:1], stat='upper'))
                rdL_ppResidues[rad]     += float(self.get_ppResidues_stat(k1=triad1[0:1],k2=triad2[0:1], stat='lower'))
                
                rad_ppEnvironment[rad]  += float(self.get_ppEnvironment(k1 = triad1[2:],k2 = triad2[2:],distance=self.contact_distances[x]))
                rdM_ppEnvironment[rad]  += float(self.get_ppEnvironment_stat(k1=triad1[2:],k2=triad2[2:], stat='mean'))
                rdU_ppEnvironment[rad]  += float(self.get_ppEnvironment_stat(k1=triad1[2:],k2=triad2[2:], stat='upper'))
                rdL_ppEnvironment[rad]  += float(self.get_ppEnvironment_stat(k1=triad1[2:],k2=triad2[2:], stat='lower'))
                
                rad_localppRE[rad]      += float(self.get_localppRE(k=triad1)) + float(self.get_localppRE(k=triad2))
                
                rad_ppRE[rad]           += float(self.get_ppRE(k1 = triad1, k2 = triad2,distance=self.contact_distances[x]))
                rdM_ppRE[rad]           += float(self.get_ppRE_stat(k1=triad1,k2=triad2, stat='mean'))
                rdU_ppRE[rad]           += float(self.get_ppRE_stat(k1=triad1,k2=triad2, stat='upper'))
                rdL_ppRE[rad]           += float(self.get_ppRE_stat(k1=triad1,k2=triad2, stat='lower'))
                
                rad_ppDist[rad]         += float(self.get_ppDist(k=self.contact_distances[x]))
                
                
                if math.floor(self.contact_distances[x]) <= self.cutoff:
                    rdC_ppResidues[rad]     += float(self.get_ppResidues(k1=triad1[0:1],k2=triad2[0:1], distance=self.cutoff))
                
                    rdC_ppEnvironment[rad]  += float(self.get_ppEnvironment(k1 = triad1[2:],k2 = triad2[2:], distance=self.cutoff))
                
                    rdC_localppRE[rad]      += float(self.get_localppRE(k=triad1)) + float(self.get_localppRE(k=triad2))
                
                    rdC_ppRE[rad]           += float(self.get_ppRE(k1 = triad1, k2 = triad2, distance=self.cutoff))
                
                    rdC_ppDist[rad]         += float(self.get_ppDist(k=self.cutoff))
                    
 
            rad_combined[rad] = rad_ppRE[rad] + rad_ppDist[rad] - rad_localppRE[rad] - rad_ppEnvironment[rad]
            rdC_combined[rad] = rdC_ppRE[rad] + rdC_ppDist[rad] - rdC_localppRE[rad] - rdC_ppEnvironment[rad]
            rdM_combined[rad] = rdM_ppRE[rad] + rad_ppDist[rad] - rad_localppRE[rad] - rdM_ppEnvironment[rad]
            rdU_combined[rad] = rdU_ppRE[rad] + rad_ppDist[rad] - rad_localppRE[rad] - rdU_ppEnvironment[rad]
            rdL_combined[rad] = rdL_ppRE[rad] + rad_ppDist[rad] - rad_localppRE[rad] - rdL_ppEnvironment[rad]

        #sys.stderr.write("\n")

        # Add the original as an extra random
        rad_ppResidues[-1]      = ori_ppResidues
        rad_ppEnvironment[-1]   = ori_ppEnvironment
        rad_localppRE[-1]       = ori_localppRE
        rad_ppRE[-1]            = ori_ppRE
        rad_ppDist[-1]          = ori_ppDist
        rad_combined[-1]        = ori_combined
        
        
        result.set_random_vector(vector = rad_ppResidues,     label="D-ZScore")
        result.set_random_vector(vector = rad_ppEnvironment,  label="D-ZE3Denv")
        result.set_random_vector(vector = rad_localppRE,      label="D-ZElocal")
        result.set_random_vector(vector = rad_ppRE,           label="D-ZEaa3Denv")
        result.set_random_vector(vector = rad_combined,       label="D-ZEnergy")
        
        result.set_random_vector(vector = rdC_ppResidues,     label="C-ZScore")
        result.set_random_vector(vector = rdC_ppEnvironment,  label="C-ZE3Denv")
        result.set_random_vector(vector = rdC_localppRE,      label="C-ZElocal")
        result.set_random_vector(vector = rdC_ppRE,           label="C-ZEaa3Denv")
        result.set_random_vector(vector = rdC_combined,       label="C-ZEnergy")
        
        result.set_random_vector(vector = rdM_ppResidues,     label="M-ZScore")
        result.set_random_vector(vector = rdM_ppEnvironment,  label="M-ZE3Denv")
        result.set_random_vector(vector = rdM_ppRE,           label="M-ZEaa3Denv")
        result.set_random_vector(vector = rdM_combined,       label="M-ZEnergy")
        
        result.set_random_vector(vector = rdU_ppResidues,     label="U-ZScore")
        result.set_random_vector(vector = rdU_ppEnvironment,  label="U-ZE3Denv")
        result.set_random_vector(vector = rdU_ppRE,           label="U-ZEaa3Denv")
        result.set_random_vector(vector = rdU_combined,       label="U-ZEnergy")
        
        result.set_random_vector(vector = rdL_ppResidues,     label="L-ZScore")
        result.set_random_vector(vector = rdL_ppEnvironment,  label="L-ZE3Denv")
        result.set_random_vector(vector = rdL_ppRE,           label="L-ZEaa3Denv")
        result.set_random_vector(vector = rdL_combined,       label="L-ZEnergy")
        

        return result

    """
        Getters
    """
    def get_ppResidues(self,k1,k2,distance):
        d = int(math.floor(distance))
        if self.ppResidues.has_key(k1) and self.ppResidues[k1].has_key(k2):
            return self.ppResidues[k1][k2][d]
        elif self.ppResidues.has_key(k2) and self.ppResidues[k2].has_key(k1):
            return self.ppResidues[k2][k1][d]
        else:
            return 0
        
    def get_ppResidues_stat(self,k1,k2,stat):
        data = None
        if self.ppResidues.has_key(k1) and self.ppResidues[k1].has_key(k2):
            data = numpy.array(self.ppResidues[k1][k2], dtype='float64')
        elif self.ppResidues.has_key(k2) and self.ppResidues[k2].has_key(k1):
            data = numpy.array(self.ppResidues[k2][k1], dtype='float64')
        else:
            return 0
        
        if stat == 'upper':
            return numpy.max(data)
        elif stat == 'lower':
            return numpy.min(data)
        elif stat == 'mean':
            return numpy.mean(data)
        
    def get_ppEnvironment(self,k1,k2,distance):
        d = int(math.floor(distance))
        if self.ppEnvironment.has_key(k1) and self.ppEnvironment[k1].has_key(k2):
            return self.ppEnvironment[k1][k2][d]
        elif self.ppEnvironment.has_key(k2) and self.ppEnvironment[k2].has_key(k1):
            return self.ppEnvironment[k2][k1][d]
        else:
            return 0

    def get_ppEnvironment_stat(self,k1,k2,stat):
        data = None
        if self.ppEnvironment.has_key(k1) and self.ppEnvironment[k1].has_key(k2):
            data = numpy.array(self.ppEnvironment[k1][k2], dtype='float64')
        elif self.ppEnvironment.has_key(k2) and self.ppEnvironment[k2].has_key(k1):
            data = numpy.array(self.ppEnvironment[k2][k1], dtype='float64')
        else:
            return 0
        
        # UPPER and LOWER are inverted in this potential
        if stat == 'upper':
            return numpy.min(data)
        elif stat == 'lower':
            return numpy.max(data)
        elif stat == 'mean':
            return numpy.mean(data)  
              
    def get_localppRE(self,k):
        return self.localppRE[k]
    
    def get_ppRE(self,k1,k2,distance):
        d = int(math.floor(distance))
        if self.ppRE.has_key(k1) and self.ppRE[k1].has_key(k2):
            return self.ppRE[k1][k2][d]
        elif self.ppRE.has_key(k2) and self.ppRE[k2].has_key(k1):
            return self.ppRE[k2][k1][d]
        else:
            return 0

    def get_ppRE_stat(self,k1,k2,stat):
        data = None
        if self.ppRE.has_key(k1) and self.ppRE[k1].has_key(k2):
            data = numpy.array(self.ppRE[k1][k2], dtype='float64')
        elif self.ppRE.has_key(k2) and self.ppRE[k2].has_key(k1):
            data = numpy.array(self.ppRE[k2][k1], dtype='float64')
        else:
            return 0
        
        if stat == 'upper':
            return numpy.max(data)
        elif stat == 'lower':
            return numpy.min(data)
        elif stat == 'mean':
            return numpy.mean(data)  
        
    def get_ppDist(self,k):
        d = str(int(math.floor(k)))
        return self.ppDist[d]
    
#    def get_ppAccumPair(self,k1,k2,distance):
#        d = int(math.floor(distance))
#        if self.ppAccumPair.has_key(k1) and self.ppAccumPair[k1].has_key(k2):
#            return self.ppAccumPair[k1][k2][d]
#        else:
#            return self.ppAccumPair[k2][k1][d]
        
    """
        Create a random array substituting a previous one, it creates a random number to substitute each other number
    """
    def randomise_interaface(self,list,range):
        rand = []
        done = {}
        for pos in list:
            if done.has_key(pos):
                rand.append(done[pos])
            else:
                r = random.randint(0,(range-1))
                rand.append(r)
                done[pos] = r
        return rand
        
    """
        IO FUNCTIONS
        ------------
    """        
    def process_file_line(self, line, tag):
        
        if tag == 0:    return
        if tag == 1 and line.startswith(" PAIR"): self.assign_multiple_distances(dic=self.ppResidues, keys= line.split()[1],
                                                                                 values = line.split()[2:])
        if tag == 2 and line.startswith(" PAIR"): self.assign_multiple_distances(dic=self.ppEnvironment, keys= line.split()[1],
                                                                                 values = line.split()[2:])
        if tag == 3 and line.startswith(" AA-ENV"): self.assign_single_distance(dic=self.localppRE, key= line.split()[1],
                                                                                value = line.split()[2])
        if tag == 4 and line.startswith(" PAIR"): self.assign_multiple_distances(dic=self.ppRE, keys= line.split()[1],
                                                                                 values = line.split()[2:])
        if tag == 5 and line.startswith(" Distance"): self.assign_single_distance(dic=self.ppDist, key= line.split()[1],
                                                                                  value = line.split()[2])
#       if tag == 7 and line.startswith(" PAIR"): self.assign_multiple_distances(dic=self.ppAccumPair, keys= line.split()[1],
#                                                                                values = line.split()[2:])
        
    def assign_multiple_distances(self, dic, keys, values):
        
        dic.setdefault(keys.split(':')[0],{}).setdefault(keys.split(':')[1],[]).extend(values)
        
    def assign_single_distance(self, dic, key, value):
        
        dic[key] = value

    def two_keys_dic_multivalue_2string_array(self, dic, prev):
        
        lines = []
        for k1 in dic:
            for k2 in dic[k1]:
                line = prev + "%s:%s\t" %(k1, k2)
                line = line + "\t".join(dic[k1][k2])
                lines.append(line)
        return lines
                
    def one_key_dic_value_2string_array(self, dic, prev):
        
        lines = []
        for k in dic:
            line = prev + "\t".join((k,dic[k]))
            lines.append(line)
        return lines
    
    def ppResidues_2string_array(self):
        
        lines = []
        
        lines.append("PAIR-POTENTIALS OF RESIDUES")
        lines.extend(self.two_keys_dic_multivalue_2string_array(dic=self.ppResidues, prev = " PAIR "))
        lines.append("END")
        
        return lines
    
    def string_ppResidues(self): return "\n".join(self.ppResidues_2string_array())
    
    def ppEnvironment_2string_array(self):
        
        lines = []
        
        lines.append("PAIR-POTENTIALS OF ENVIRONMENTS")
        lines.extend(self.two_keys_dic_multivalue_2string_array(dic=self.ppEnvironment, prev = " PAIR "))
        lines.append("END")
        
        return lines
    
    def string_ppEnvironment(self): return "\n".join(self.ppEnvironment_2string_array())
    
    def localppRE_2string_array(self):
        
        lines = []
        
        lines.append("LOCAL PAIR-POTENTIALS OF RESIDUES AND ENVIRONMENTS")
        lines.extend(self.one_key_dic_value_2string_array(dic=self.localppRE, prev = " AA-ENV "))
        lines.append("END")
        
        return lines
    
    def string_localppRE(self): return "\n".join(self.localppRE_2string_array())
        
    def ppRE_2string_array(self):
        
        lines = []
        
        lines.append("PAIR-POTENTIALS OF RESIDUES AND ENVIRONMENTS")
        lines.extend(self.two_keys_dic_multivalue_2string_array(dic=self.ppRE, prev = " PAIR "))
        lines.append("END")
        
        return lines
    
    def string_ppRE(self): return "\n".join(self.ppRE_2string_array())
    
    def ppDist_2string_array(self):
        
        lines = []
        
        lines.append("STANDARD PAIR-POTENTIALS OF DISTANCES")
        lines.extend(self.one_key_dic_value_2string_array(dic=self.ppDist, prev = " Distance "))
        lines.append("END")
        
        return lines
    
    def string_ppDist(self): return "\n".join(self.ppDist_2string_array())
    
#    def ppAccumPair_2string_arra(self):
#        
#        lines = []
#        
#        lines.append("PAIR-POTENTIALS OF ACCUMULATED PAIRS")
#        lines.extend(self.two_keys_dic_multivalue_2string_array(dic=self.ppAccumPair, prev = " PAIR "))
#        lines.append("END")
#        
#        return lines        
#    
#    def string_ppAccumPair(self): return "\n".join(self.ppAccumPair_2string_arra())
              
    def __str__(self):
        
        lines = []
        
        lines.extend(self.ppResidues_2string_array())
        lines.extend(self.ppEnvironment_2string_array())
        lines.extend(self.localppRE_2string_array())
        lines.extend(self.ppRE_2string_array())
        lines.extend(self.ppDist_2string_array())
#       lines.extend(self.ppAccumPair_2string_arra())
        
        return "\n".join(lines)
        
        
