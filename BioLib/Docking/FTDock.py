import re, copy
import BioLib.Structure.PDB

class FTDock():
    '''
    '''

    static = None
    mobile = None
    gridSize = None
    searchAngleStep = 12
    surfaceThickness = 1.30
    internalDeterrentDalue = -15.00
    electrostatics = 1
    keepPerRotation = 3
    rotations = None
    totalSpan = None
    gridCellSpan = None
 
    decoys = []

    def __init__(self, ftdockFile):
        '''
        Call the parsing methods (parses an ftdock_global.dat file)
        '''
        self.__parse_info(ftdockFile)
        self.__parse_decoys(ftdockFile)

    def __parse_info(self, ftdockFile):
        '''
        '''
        ftdockfo = open(ftdockFile, 'r')
        for line in ftdockfo.readlines():
            static =  re.search('Static molecule\s+::\s+([\w\.]+)', line) 
            mobile = re.search('Mobile molecule\s+::\s+([\w\.]+)', line)
            gridSize = re.search('Global grid size\s+::\s+([\d\.-]+)', line)
            searchAngleStep = re.search('Global search angle step\s+::\s+([\d\.-]+)', line)
            surfaceThickness = re.search('Global surface thickness\s+::\s+([\d\.-]+)', line)
            internalDeterrentDalue = re.search('Global internal deterrent value\s+::\s+([\d\.-]+)', line)
            electrostatics = re.search('Electrostatics\s+::\s+(\w+)', line)
            keepPerRotation = re.search('Global keep per rotation\s+::\s+([\d\.-]+)', line)
            rotations = re.search('Global rotations\s+::\s+([\d\.-]+)', line)
            totalSpan = re.search('Global total span \(angstroms\)\s+::\s+([\d\.-]+)', line)
            gridCellSpan = re.search('Global grid cell span \(angstroms\)\s+::\s+([\d\.-]+)', line)

            if static: self.static = static.group(1)
            elif mobile: self.mobile = mobile.group(1)
            elif gridSize: self.gridSize = gridSize.group(1)
            elif searchAngleStep: self.searchAngleStep = searchAngleStep.group(1)
            elif surfaceThickness: self.surfaceThickness = surfaceThickness.group(1)
            elif internalDeterrentDalue: self.internalDeterrentDalue = internalDeterrentDalue.group(1)
            elif electrostatics:
                if electrostatics.group(1) != 'on': self.electrostatics = 0
            elif keepPerRotation: self.keepPerRotation = keepPerRotation.group(1)
            elif rotations: self.rotations = rotations.group(1)
            elif totalSpan: self.totalSpan = totalSpan.group(1)
            elif gridCellSpan: self.gridCellSpan = gridCellSpan.group(1)

    def __parse_decoys(self, ftdockFile):
        '''
        '''
        ftdockfo = open(ftdockFile, 'r')
        for line in ftdockfo.readlines():
            if line[0:6] == 'G_DATA':
                self.decoys.append(FTDecoy(int(line[8:14]), float(line[27:34]), float(line[40:49]), int(line[55:59]), int(line[59:64]), int(line[64:69]), int(line[75:79]), int(line[79:83]), int(line[83:])))
   
    def get_grid_size(self):
        return self.gridSize

    def get_decoys(self):
        return self.decoys
 
    def get_decoy_by_id(self,id):

        for decoy in self.decoys:
            if decoy.get_id() == id:
                return decoy
        return None

    def get_grid_total_span(self, static, mobile):
        '''
        Calculate the total span of the structures
        '''
        return 1 + ((static.get_radius() + mobile.get_radius()) * 2)

    def get_one_span(self, static, mobile):
        '''
        Calculate the GRID one span
        '''
        return self.get_grid_total_span(static, mobile)/float(self.gridSize) 

    def get_static_structure(self, static):
        '''
        Gets the static structure of the decoys (same in all decoys)
        '''
        staticDecoy = copy.deepcopy(static)
        staticDecoy.translate_onto_origin()
        return staticDecoy

    def get_mobile_structure(self, static, mobile, decoy):
        '''
        Gets the structure of each decoy (ftdock only relocates the mobile structure)
        '''
        mobileDecoy = copy.deepcopy(mobile)
        mobileDecoy.translate_onto_origin()

        translationVector = decoy.get_corrected_tranlation_vector(self.get_one_span(self.get_static_structure(static), mobileDecoy))
        twistVector = decoy.get_rotation_vector()
        mobileDecoy.relocate(translationVector=translationVector, twistVector=twistVector)        

        return mobileDecoy

    def get_ligand_RMSD(self, static, mobile, native, decoy):
        '''
        Gets the ligand RMSD between one decoy and the original (native) mobile structure
        '''
        originalMobile = copy.deepcopy(native)
        staticDecoy = copy.deepcopy(static)
        mobileDecoy = self.get_mobile_structure(static, mobile, decoy)

        centerTraslationVector = staticDecoy.translate_onto_origin()
        originalMobile.relocate(translationVector=centerTraslationVector)

        return originalMobile.get_RMSD(mobileDecoy)

    def get_decoys_ligand_RMSD(self, static, mobile, native):
        '''
        Gets a dictionary of the ligand RMSD between each decoy and the original (native) mobile structure
        '''
        ligandRMSD = {}
        for decoy in self.get_decoys():
            ligandRMSD[decoy.get_id()] = self.get_ligand_RMSD(static, mobile, native, decoy)
            print str(decoy.get_id())+'\t'+str(ligandRMSD[decoy.get_id()])
        return ligandRMSD
    
    def print_RMSD(self, static, mobile, native, RMSDFile):
        '''
        Prints the ligand RMSD of each decoy in a tab file
        '''
        RMSDfo = open(RMSDFile, 'w')
        ligandRMSD = self.get_decoys_ligand_RMSD(static, mobile, native)
        for key, value in sorted(ligandRMSD.iterkeys(), key=lambda (k,v): (v,k)):
            RMSDfo.write('%d\t%.3f\n' % (key, value))
        RMSDfo.close()

    def print_decoys(self, static, mobile, pdbFile, singlePDB=False):
        '''
        Print each decoy in a PDB format
        @singlePDB = Prints each decoy as a model in a single PDB file
        '''
        staticDecoy = self.get_static_structure(static)
        decoyList = []
        for decoy in self.get_decoys():
            mobileDecoy = self.get_mobile_structure(static, mobile, decoy)
            if singlePDB:
                DecoyList.append((staticDecoy, mobileDecoy))
            else:
                PDB.write_pdb((staticDecoy, mobileDecoy), pdbFile+'_'+str(decoy.get_id())+'.pdb' , multiChain=True, multiModel=False)
        if singlePDB:
            PDB.write_pdb(decoyList, pdbFile+'.pdb', multiChain=True, multiModel=True)
    
    def __str__(self):
        '''
        '''
        return "Static molecule: %s\nMobile molecule: %s\nGlobal grid size: %s\nGlobal search angle step: %s\nGlobal surface thickness: %s\nGlobal internal deterrent value: %s\nElectrostatics: %s\nGlobal keep per rotation: %s\nGlobal rotations: %s\nGlobal total span (angstroms): %s\nGlobal grid cell span (angstroms): %s" % (self.static, self.mobile, self.gridSize, self.searchAngleStep, self.surfaceThickness, self.internalDeterrentDalue, self.electrostatics, self.keepPerRotation, self.rotations, self.totalSpan, self.gridCellSpan)


class FTDecoy():
    '''
    '''
    def __init__(self, id, SCscore, ESratio, x, y, z, ztwist, theta, phi):
        '''
        '''
        self.id = id
        self.SCscore = SCscore
        self.ESratio = ESratio
        self.x = x
        self.y = y
        self.z = z
        self.ztwist = ztwist
        self.theta = theta
        self.phi = phi

    def get_id(self):
        return self.id

    def get_SCscore(self):
        return self.SCscore

    def get_ESratio(self):
        return self.ESratio

    def get_tranlation_vector(self):
        return (self.x, self.y, self.z)
   
    def get_corrected_tranlation_vector(self, one_span):
        '''
        Gets the corrected translation vector (needed for ftdock)
        '''
        return (self.x * one_span, self.y * one_span, self.z * one_span)

    def get_rotation_vector(self):
        return (self.ztwist, self.theta, self.phi)

    def __str__(self):
        '''
        '''
        return ('G_DATA\t%d\t%d\t%.3f\t%s%s%s\t%s%s%s') % (self.id, self.SCscore, self.ESratio, str(self.x).rjust(5), str(self.y).rjust(5), str(self.z).rjust(5), str(self.ztwist).rjust(5), str(self.theta).rjust(5), str(self.phi).rjust(5))
