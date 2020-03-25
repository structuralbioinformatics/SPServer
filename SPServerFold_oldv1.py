from BioLib import *
import optparse, os, sys 
import tarfile, zipfile

def parse_options():
    '''
    Parse arguments
    @return options object
    '''
    parser = optparse.OptionParser('SPServerFold.py -f PDB_PATH -p POT_TYPE -o OUT_PATH -x XML_NAME')
    parser.add_option("-f", dest="pdb_path", action="store", help="Path to compress (tar.gz or zip) file with PDBs [MANDATORY]", metavar="PDB_PATH")
    parser.add_option("-p", dest="pot_type", action="store", help="CB or Min atom-atom distance potential", choices = ['CB', 'MIN'], default='CB', metavar="POT_TYPE")
    parser.add_option("-o", dest="out_path", action="store", help="Out path [MANDATORY]", metavar="OUT_PATH")
    parser.add_option("-x", dest="xml_name", action="store", help="XML name [MANDATORY]", metavar="OUT_PATH")
    (options, args) = parser.parse_args()

    if options.pdb_path == None or options.out_path == None and options.xml_name == None:
        parser.error('Missing arguments\n')
    if not os.path.isfile(options.pdb_path):
        parser.error('PDB_PATH is not a file\n')
    if not os.path.isdir(options.out_path):
        parser.error('OUT_PATH is not a directory\n')
    return options

class SPServerFold(object):
    '''
    Analyze a protein folds using the split potentials and prints an XML file
    '''
    def __init__(self, pdb_path, out_path, pot_type, xml_name):
        '''
        Contructor
        '''
        self.pdb_path = pdb_path
        self.out_path = out_path
        self.pot_type = pot_type
        self.xml_name = xml_name
        self.xml_result = '<?xml version="1.0" encoding="utf-8"?>\n<xml>\n'
        self.xml_errors = '<?xml version="1.0" encoding="utf-8"?>\n<xml>\n'

    def run(self):
        '''
        Runs the computations
        '''
        pdb_filename_path = out_path + '/folds' # Future path containing the extracted pdbs
        # Decompress the file with pdb folds
        if self.pdb_path.endswith('.tar.gz'):
            tfile = tarfile.open(self.pdb_path, 'r:gz')
            tfile.extractall(pdb_filename_path)
            #pdb_filename = os.path.splitext(os.path.splitext(os.path.basename(self.pdb_path))[0])[0]
        elif self.pdb_path.endswith('.tgz'):
            tfile = tarfile.open(self.pdb_path, 'r:gz')
            tfile.extractall(pdb_filename_path)
            #pdb_filename = os.path.splitext(os.path.basename(self.pdb_path))[0]  
        elif self.pdb_path.endswith('.zip'):
            zfile = zipfile.ZipFile(self.pdb_path)
            zfile.extractall(pdb_filename_path)
            #pdb_filename = os.path.splitext(os.path.basename(self.pdb_path))[0]
        elif self.pdb_path.endswith('.pdb'):
            pdb = os.path.basename(self.pdb_path)
            os.mkdir(pdb_filename_path)
            os.rename(pdb_path, pdb_filename_path + '/' + pdb)
            #pdb_filename = os.path.splitext(os.path.basename(self.pdb_path))[0]
        else:
            self.print_error(1)
            return None
        sys.stdout.write('[SUCCESS] File decompression done!\n')

        os.system('chmod 777 {}'.format(pdb_filename_path)) # Give permissions to the extraction path

        # Check if the compressed file was a folder containing the PDBs
        # If it is the case, the variable "pdb_filename_path" is reassigned to this folder
        for item in os.listdir(pdb_filename_path):
            item = pdb_filename_path + '/' + item
            if os.path.isdir(item) == True:
                pdb_filename_path = item
                os.system('chmod 777 {}'.format(pdb_filename_path))

        # Iterate over each pdb fold and compute the split potentials
        job_pdb_path = pdb_filename_path
        #job_pdb_path = os.path.join(self.out_path, pdb_filename)
        os.system('chmod 777 {}/*'.format(pdb_filename_path)) # Give permissions to the pdb files inside the extraction path

        restart = True
        while restart:
            restart = False
            for fold_name in os.listdir(job_pdb_path):
                if not fold_name.endswith('.pdb'):
                    self.print_error(2, fold_name)
                    continue
                error_n = self.on_fold(job_pdb_path, fold_name)
                # If on_fold returns 'split', it means that the user has chosen to do the analysis splitting the chains of a PDB complex
                # And the on_fold function has splitted the chains and removed the original pdb
                # So, we restart the loop in order to analyze each one of the new chains extracted from the complex
                if error_n == 'split':
                    restart = True
                    break
                if error_n > 0:
                    self.print_error(error_n, fold_name)

        sys.stdout.write('[SUCCESS] Split Potentials computed!\n')
        # Print the xml result file
        xml_result_fo = open(os.path.join(self.out_path, self.xml_name+'.xml'), 'w')
        xml_result_fo.write(self.xml_result+'</xml>')
        xml_result_fo.close()
        # Print the error xml file
        xml_err_fo = open(os.path.join(self.out_path, self.xml_name+'.err.xml'), 'w')
        xml_err_fo.write(self.xml_errors+'</xml>')
        xml_err_fo.close()
        sys.stdout.write('[SUCCESS] Finished!\n')

    def on_fold(self, job_pdb_path, fold_name):
        '''
        Compute the split potentials for a fold
        '''
        # Get PDB structure (fold structure)
        fold_path = os.path.join(job_pdb_path, fold_name)
        struct = PDB.read_pdb(fold_path)

        if isinstance(struct, list): # if the structure obtained is an instance of a list, it means that it is not a unique structure: it is a COMPLEX of more than one chain
            process_complex = 'merge' # finally we will not let the user choose, so a complex will always be merged
            # (If the user chooses to split the complex, we split it, remove it and return 'split' in order to restart the analysis with the chains)
            if process_complex == 'split':
                for chain in struct:
                    chain_name = chain.get_chain()
                    chain_path = job_pdb_path + '/' + fold_name.split('.')[0] + '_' + chain_name + '.pdb'
                    PDB.write_pdb(chain, chain_path, multi_chain=False, multi_model=False)
                os.remove(fold_path)
                return 'split'
            # If the user chooses to merge the complex, we merge it using read_pdb()
            elif process_complex == 'merge':
                struct = PDB.read_pdb(fold_path, merge_chains=True)

        struct.set_dssp()
        struct.normalize_residues()
        # Compute split potentials
        if self.pot_type == 'CB':
            cutoff = 12
        else: 
            cutoff = 5
        split_potentials = SplitPotentialsFold(c_type=self.pot_type, cutoff=cutoff) # definition of SplitPotentialsFold class (BioLib.Fold)
        global_energies = split_potentials.calculate_global_energies(struct, Zscores=True) # calculation of global energies using the method of the SplitPotentialFold class
        residues_energies = split_potentials.calculate_residue_energies(struct, Zscores=True) # same for residues energies
        self.fold_xml(fold_name, global_energies, residues_energies) # creation of xml results file

    def fold_xml(self, fold_name, global_energies, residues_energies):
        '''
        Append the fold energies to the xml result file
        '''
        self.xml_result += '<protein>\n'
        self.xml_result += '  <id>%s</id>\n' % fold_name
        # Global Energies (and Zenergies)
        self.xml_result += '  <global_energies>\n'
        self.xml_result += '    <Epair>%.2f</Epair>\n'     % global_energies[0]['D-PAIR']
        self.xml_result += '    <Ecomb>%.2f</Ecomb>\n'     % global_energies[0]['D-COMB']
        self.xml_result += '    <Es3dc>%.2f</Es3dc>\n'     % global_energies[0]['D-S3DC']
        self.xml_result += '    <Elocal>%.2f</Elocal>\n'   % global_energies[0]['D-LOCAL']
        self.xml_result += '    <E3dc>%.2f</E3dc>\n'       % global_energies[0]['D-3DC']
        self.xml_result += '    <E3d>%.2f</E3d>\n'         % global_energies[0]['D-3D']
        self.xml_result += '    <ZEpair>%.2f</ZEpair>\n'   % global_energies[1]['D-PAIR']
        self.xml_result += '    <ZEcomb>%.2f</ZEcomb>\n'   % global_energies[1]['D-COMB']
        self.xml_result += '    <ZEs3dc>%.2f</ZEs3dc>\n'   % global_energies[1]['D-S3DC']
        self.xml_result += '    <ZElocal>%.2f</ZElocal>\n' % global_energies[1]['D-LOCAL']
        self.xml_result += '    <ZE3dc>%.2f</ZE3dc>\n'     % global_energies[1]['D-3DC']
        self.xml_result += '    <ZE3d>NA</ZE3d>\n'
        self.xml_result += '  </global_energies>\n'
        # Energies by residue (and Zenergies)
        self.xml_result += '  <residues>\n'
        for residue_energies in residues_energies:
            self.xml_result += '    <residue>\n'
            self.xml_result += '      <number>%d</number>\n'     % residue_energies[0].get_num()
            self.xml_result += '      <type>%s</type>\n'         % residue_energies[0].get_type_short()
            self.xml_result += '      <Epair>%.2f</Epair>\n'     % residue_energies[1]['D-PAIR']
            self.xml_result += '      <Ecomb>%.2f</Ecomb>\n'     % residue_energies[1]['D-COMB']
            self.xml_result += '      <Es3dc>%.2f</Es3dc>\n'     % residue_energies[1]['D-S3DC']
            self.xml_result += '      <Elocal>%.2f</Elocal>\n'   % residue_energies[1]['D-LOCAL']
            self.xml_result += '      <E3dc>%.2f</E3dc>\n'       % residue_energies[1]['D-3DC']
            self.xml_result += '      <E3d>%.2f</E3d>\n'         % residue_energies[1]['D-3D']
            self.xml_result += '      <ZEpair>%.2f</ZEpair>\n'   % residue_energies[2]['D-PAIR']
            self.xml_result += '      <ZEcomb>%.2f</ZEcomb>\n'   % residue_energies[2]['D-COMB']
            self.xml_result += '      <ZEs3dc>%.2f</ZEs3dc>\n'   % residue_energies[2]['D-S3DC']
            self.xml_result += '      <ZElocal>%.2f</ZElocal>\n' % residue_energies[2]['D-LOCAL']
            self.xml_result += '      <ZE3dc>%.2f</ZE3dc>\n'     % residue_energies[2]['D-3DC']
            self.xml_result += '      <ZE3d>NA</ZE3d>\n'
            self.xml_result += '    </residue>\n'
        self.xml_result += '  </residues>\n'
        self.xml_result += '</protein>\n'

    def print_error(self, error_num, fold=''):
        '''
        Prints the error in the error xml file and the stderr channel
        '''
        xml_template = '<error><number>%d</number><fold>%s</fold><description>%s</description></error>\n'
        stderr_template = '[ERROR %d %s] %s\n'
        if error_num == 1:
            description = 'File is not compressed in tar.gz or .zip (wrong suffix)'
        if error_num == 2:
            description = 'Not a PDB (wrong suffix)'
        if error_num == 3:
            description = 'A PDB file contains more than one chain (ppi problem, not a folding problem)'
        self.xml_errors += xml_template % (error_num, fold, description)
        sys.stderr.write(stderr_template % (error_num, fold, description))

#--------------------------#
#          MAIN            #
#--------------------------#

if __name__ == "__main__":

    options = parse_options()

    pdb_path = os.path.abspath(options.pdb_path)
    out_path = os.path.abspath(options.out_path)
   
    spserverfold = SPServerFold(pdb_path, out_path, options.pot_type, options.xml_name) 
    spserverfold.run()
