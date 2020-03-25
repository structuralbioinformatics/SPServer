import sys, os
import argparse
import tarfile
import zipfile
from BioLib import *
import SBI.structure.PDB


def main():

    options = parse_user_arguments()
    calculate_sps_ppi(options)

def parse_user_arguments(*args, **kwds):
    """
    Parses the arguments of the program
    """

    parser = argparse.ArgumentParser(
        description = "Calculate Split-Statistical Potentials",
        epilog      = "@oliva's lab 2018")
    parser.add_argument('-i','--input_file',dest='input_file',action = 'store', help = """ Path to the input file with the structure(s). """)
    parser.add_argument('-r','--chain_receptor',dest='chain_receptor',action = 'store', default = 'A', help = 'Letters corresponding to the receptor chain in the structure (only when ppi_source=together).')
    parser.add_argument('-l','--chain_ligand',dest='chain_ligand',action = 'store', default = 'B', help = 'Letters corresponding to the ligand chain in the structure (only when ppi_source=together).')
    parser.add_argument('-s','--ppi_source',dest='ppi_source',action = 'store', help = 'Type of protein-protein interaction source. It can be pdb_separated (two PDB structures) or pdb_together (one PDB structure).')
    parser.add_argument('-o','--output_dir',dest='output_dir',action = 'store', help = """ Path to the directory where the results of the job will be saved. """)
    parser.add_argument('-j','--job_id',dest='job_id',action = 'store', help = """ Job ID. """)
    parser.add_argument('-p','--pot_type',dest='pot_type',action = 'store', default = 'CB', help = """CB or MIN atom-atom distance potential (default=CB)""")
    parser.add_argument('-c','--crashes',dest='crashes',action = 'store_true', help = 'Flag to use calculate the atoms crashing between the two structures of the PPI.')

    options=parser.parse_args()

    # Example (ppi_source = 'pdb_together'):
    # python /home/quim/PHD/Projects/SPServer/standalone/SPServerPPI.py -i /home/quim/PHD/Projects/SPServer/standalone/examples/ppi_pdb_together/BAX_G108V_BID.pdb -r A -l C -s pdb_together -o /home/quim/PHD/Projects/SPServer/standalone/examples/ppi_pdb_together -j BAX_G108V_BID_new_potentials_CB -p CB -c -cp /home/quim/PHD/Projects/SPServer/standalone/collision_detection_program
    # python /home/quim/PHD/Projects/SPServer/standalone/SPServerPPI.py -i /home/quim/PHD/Projects/SPServer/standalone/examples/ppi_pdb_together/BAX_BID_reference.pdb -r A -l C -s pdb_together -o /home/quim/PHD/Projects/SPServer/standalone/examples/ppi_pdb_together -j BAX_BID_reference_new_potentials_CB -p CB -c -cp /home/quim/PHD/Projects/SPServer/standalone/collision_detection_program
        # Without collision
        # python /home/quim/PHD/Projects/SPServer/standalone/SPServerPPI.py -i /home/quim/PHD/Projects/SPServer/standalone/examples/ppi_pdb_together/BAX_G108V_BID.pdb -r A -l C -s pdb_together -o /home/quim/PHD/Projects/SPServer/standalone/examples/ppi_pdb_together -j BAX_G108V_BID_new_potentials_CB -p CB
        # python /home/quim/PHD/Projects/SPServer/standalone/SPServerPPI.py -i /home/quim/PHD/Projects/SPServer/standalone/examples/ppi_pdb_together/BAX_BID_reference.pdb -r A -l C -s pdb_together -o /home/quim/PHD/Projects/SPServer/standalone/examples/ppi_pdb_together -j BAX_BID_reference_new_potentials_CB -p CB

    # Example (ppi_source = 'pdb_together') with input list:
    # python /home/quim/PHD/Projects/SPServer/standalone/SPServerPPI.py -i /home/quim/PHD/Projects/SPServer/standalone/examples/ppi_pdb_together_list/BAX_BID.list -r A -l C -s pdb_together -o /home/quim/PHD/Projects/SPServer/standalone/examples/ppi_pdb_together_list -j BAX_BID_new_potentials -p CB

    # Example (ppi_source = 'pdb_separated'):
    # python /home/quim/PHD/Projects/SPServer/standalone/SPServerPPI.py -i /home/quim/PHD/Projects/SPServer/standalone/examples/ppi_pdb_separated/histones.tar.gz -s pdb_separated -o /home/quim/PHD/Projects/SPServer/standalone/examples/ppi_pdb_separated -j new_potentials_CB -p CB -c -cp /home/quim/PHD/Projects/SPServer/standalone/collision_detection_program

    return options


#################
#################
# MAIN FUNCTION #
#################
#################

def calculate_sps_ppi(options):
    """
    From a list of PDBs, it calculates the Split-Statistical Potentials
    """
    spserverppi = SPServerPPI(options.input_file, options.output_dir, options.job_id, options.pot_type)
    spserverppi.chain_receptor = options.chain_receptor.upper()
    spserverppi.chain_ligand = options.chain_ligand.upper()
    spserverppi.extract_ppi()
    if options.crashes:
        spserverppi.crashes = True
    if options.ppi_source == 'pdb_separated':
        spserverppi.prepare_ppi_separated()
    elif options.ppi_source == 'pdb_together':
        spserverppi.prepare_ppi_together()
    else:
        print('Incorrect protein-protein interaction type of source! It must be pdb or pdb_together.\n')
        sys.exit(10)

    spserverppi.run()
    return

class SPServerPPI(object):
    '''
    Analyze a protein folds using the split potentials and prints an XML file
    '''
    def __init__(self, input_file, output_dir, job_id, pot_type):
        '''
        Contructor
        '''
        self.input_file = input_file
        self.output_dir = output_dir
        self.job_id = job_id
        self.pot_type = pot_type
        self.chain_receptor = 'A'
        self.chain_ligand = 'B'
        self.structure_to_chain_receptor = {}
        self.structure_to_chain_ligand = {}
        self.crashes = False
        self.crashes_path = os.path.abspath(os.path.join(os.path.dirname(__file__), 'collision_detection_program'))
        self.xml_result = '<?xml version="1.0" encoding="utf-8"?>\n<xml>\n'
        self.xml_errors = '<?xml version="1.0" encoding="utf-8"?>\n<xml>\n'
        self.xml_err_file = os.path.join(output_dir, '{}.err.xml'.format(job_id))
        self.xml_out_file = os.path.join(output_dir, '{}.xml'.format(job_id))
        self.list_file = os.path.join(output_dir, '{}.list'.format(job_id))
        self.extraction_path = os.path.join(output_dir, 'extracted_structures') # Future path containing the extracted pdbs

    def extract_ppi(self):
        '''
        Extracts the input files
        '''
        if self.input_file.endswith('.tar.gz'):
            tfile = tarfile.open(self.input_file, 'r:gz')
            tfile.extractall(self.extraction_path)
        elif self.input_file.endswith('.tgz'):
            tfile = tarfile.open(self.input_file, 'r:gz')
            tfile.extractall(self.extraction_path)
        elif self.input_file.endswith('.zip'):
            zfile = zipfile.ZipFile(self.input_file)
            zfile.extractall(self.extraction_path)
        elif self.input_file.endswith('.pdb'):
            self.move_pdb(self.input_file)
        elif self.input_file.endswith('.cif'):
            self.move_cif_and_convert_to_pdb(self.input_file)
        elif self.input_file.endswith('.list'):
            with open(self.input_file, 'r') as inp_fd:
                for line in inp_fd:
                    structure_path, chain_r, chain_l = line.strip().split('\t')
                    if structure_path.endswith('.pdb'):
                        self.check_pdb_format(structure_path)
                        new_structure_path = self.move_pdb(structure_path)
                    elif structure_path.endswith('.cif'):
                        new_structure_path = self.move_cif_and_convert_to_pdb(structure_path)
                        self.check_pdb_format(new_structure_path)
                    else:
                        # If unknown extension, check if PDB, if so then move to folder
                        self.check_pdb_format(structure_path)
                        new_structure_path = self.move_pdb(structure_path)
                    self.structure_to_chain_receptor[new_structure_path] = chain_r
                    self.structure_to_chain_ligand[new_structure_path] = chain_l
                # Check if at least one of the elements of the list was a pdb or cif! If not return error
                if (len(self.structure_to_chain_receptor) == 0) & (len(self.structure_to_chain_ligand) == 0):
                    self.print_error(8, self.input_file)
                    self.write_output_file(self.xml_errors, self.xml_err_file)
                    self.write_output_file(self.xml_result, self.xml_out_file) # Output an empty results file
                    sys.exit(10)
        else:
            print('Incorrect type of file. Introduce a .tar.gz/.tgz/.zip/.pdb/.cif files!\n')
            self.print_error(1, self.input_file)
            self.write_output_file(self.xml_errors, self.xml_err_file)
            self.write_output_file(self.xml_result, self.xml_out_file) # Output an empty results file
            sys.exit(10)

        os.system('chmod 777 {}'.format(self.extraction_path)) # Give permissions to the extraction path

        # Check if the content extracted is a folder. 
        # If so, the elements inside the folder are moved to the main one and the folder is removed.
        for item in os.listdir(self.extraction_path):
            item = os.path.join(self.extraction_path, item)
            if os.path.isdir(item):
                os.system('chmod 777 {}'.format(self.extraction_path))
                os.system('mv {} {}'.format(os.path.join(item, '*'), self.extraction_path))
                os.system('rmdir {}'.format(item))
        return

    def move_pdb(self, pdb_path):
        file_name = self.get_fold_name_from_fold_path(pdb_path)
        new_file_name = os.path.join(self.extraction_path, file_name)
        create_directory(self.extraction_path)
        #os.rename(structure_path, new_file_name)
        os.system('cp {} {}'.format(pdb_path, new_file_name))
        return new_file_name

    def move_cif_and_convert_to_pdb(self, structure_path):
        import cif2pdb
        file_name = self.get_fold_name_from_fold_path(structure_path)
        new_file_name = os.path.join(self.extraction_path, file_name)
        create_directory(self.extraction_path)
        #os.rename(structure_path, new_file_name)
        #os.system('cp {} {}'.format(structure_path, new_file_name))
        pdb_file = new_file_name+'.pdb'
        try:
            cif2pdb.convert_cif_to_pdb(structure_path, pdb_file=pdb_file, verbose=True)
        except:
            self.print_error(6)
            self.write_output_file(self.xml_errors, self.xml_err_file)
            self.write_output_file(self.xml_result, self.xml_out_file) # Output an empty results file
            sys.exit(10)
        return pdb_file

    def check_pdb_format(self, structure_path):
        file_name = self.get_fold_name_from_fold_path(structure_path)
        try:
            struct = PDB.read_pdb(structure_path)
        except:
            self.print_error(8, file_name)
            self.write_output_file(self.xml_errors, self.xml_err_file)
            self.write_output_file(self.xml_result, self.xml_out_file) # Output an empty results file
            sys.exit(10)
        return

    def prepare_ppi_separated(self):
        '''
        Prepare the structures to be able to run the potentials
        '''
        with open(self.list_file, 'w') as list_fd:

            # Get the receptors and ligands
            receptors = []
            ligands = []
            for structure_file in os.listdir(self.extraction_path):

                if structure_file[-6:] == '.r.pdb':
                    receptors.append(structure_file)
                elif structure_file[-6:] == '.l.pdb':
                    ligands.append(structure_file)
                else:
                    print('Incorrect file!:{}. It must end with .r.pdb or .l.pdb\n'.format(structure_file))
                    self.print_error(2, structure_file)
                    self.write_output_file(self.xml_errors, self.xml_err_file)
                    self.write_output_file(self.xml_result, self.xml_out_file) # Output an empty results file
                    sys.exit(10)

            if len(receptors) == len(ligands):
                pairs = []
                if len(receptors) == 1:
                    receptor = os.path.join(self.extraction_path, receptors[0])
                    ligand = os.path.join(self.extraction_path, ligands[0])
                    pairs.append([receptor, ligand])
                else:
                    for file in receptors:
                        pdb_name = file[:-6]
                        receptor = os.path.join(self.extraction_path, file)
                        ligand = os.path.join(self.extraction_path, '{}.l.pdb'.format(pdb_name))
                        if fileExist(receptor) and fileExist(ligand):
                            pairs.append([receptor, ligand])
                        else:
                            print('The receptor {} has not a corresponding ligand! Please, name the receptor and ligand equally.\n'.format(file))
                            self.print_error(3, pdb_name)
                            self.write_output_file(self.xml_errors, self.xml_err_file)
                            self.write_output_file(self.xml_result, self.xml_out_file) # Output an empty results file
                            sys.exit(10)

                for pair in pairs:
                    # Merge the chains of the structures
                    receptor, ligand = pair
                    chain_r = self.get_merged_chain_from_pdb(receptor, 'A')
                    chain_l = self.get_merged_chain_from_pdb(ligand, 'B')

                    # Print two files, one for the receptor with merged chains and the other for the ligand merged
                    name_a = os.path.join(self.extraction_path, '{}_merged.pdb'.format(receptor[:-4]))
                    name_b = os.path.join(self.extraction_path, '{}_merged.pdb'.format(ligand[:-4]))
                    pdb_chain_a = SBI.structure.PDB()
                    pdb_chain_b = SBI.structure.PDB()
                    pdb_chain_a.add_chain(chain_r)
                    pdb_chain_b.add_chain(chain_l)
                    pdb_chain_a.write(output_file=name_a, format='PDB')
                    pdb_chain_b.write(output_file=name_b, format='PDB')

                    # Add the PDBs of the two separated chains to the list
                    # to run collision program and potentials program
                    name = '{}.pdb'.format(receptor[:-6])
                    structure_name = os.path.basename(name)
                    list_fd.write('{}\t{}\t{}\n'.format( structure_name, name_a, name_b ))
        return

    def prepare_ppi_together(self):
        '''
        Prepare the structures to be able to run the potentials
        '''
        with open(self.list_file, 'w') as list_fd:
            for structure in os.listdir(self.extraction_path):
                structure_file = os.path.join(self.extraction_path, structure)
                struct = SBI.structure.PDB(structure_file)
                chains = list(struct.chain_identifiers)
                # Get the chains identifiers from dictionary in case the input file is a list file
                if structure_file in self.structure_to_chain_receptor and structure_file in self.structure_to_chain_ligand:
                    self.chain_receptor = self.structure_to_chain_receptor[structure_file].upper()
                    self.chain_ligand = self.structure_to_chain_ligand[structure_file].upper()
                # We rename the chains,
                # because the PPI program only works with structures of two chains named A and B
                chain_r = struct.get_chain_by_id(self.chain_receptor)
                chain_l = struct.get_chain_by_id(self.chain_ligand)
                # Check that the chains are in the structure
                if not chain_r:
                    self.print_error(7, ppi=['receptor', self.chain_receptor, ', '.join(chains)])
                    self.write_output_file(self.xml_errors, self.xml_err_file)
                    self.write_output_file(self.xml_result, self.xml_out_file) # Output an empty results file
                    sys.exit(10)
                if not chain_l:
                    self.print_error(7, ppi=['ligand', self.chain_ligand, ', '.join(chains)])
                    self.write_output_file(self.xml_errors, self.xml_err_file)
                    self.write_output_file(self.xml_result, self.xml_out_file) # Output an empty results file
                    sys.exit(10)
                # If chains are not A and B, rename them to A and B
                if not (self.chain_receptor == 'A' and self.chain_ligand == 'B') and not (self.chain_receptor == 'B' and self.chain_ligand == 'A'):
                    chain_r.chain = 'A'
                    chain_l.chain = 'B'
                    pdb = SBI.structure.PDB()
                    pdb.add_chain(chain_r)
                    pdb.add_chain(chain_l)
                    # Rewrite the structure_file using the different chain names
                    pdb.write(output_file=structure_file, format='PDB', force=True)

                # Then we split the structure in order to run the Gaurav program,
                # because it only works when we have the PPI in two separated structures
                # and also the new potentials calculation
                pdb_chain_r = SBI.structure.PDB()
                pdb_chain_r.add_chain(chain_r)
                name_r = '{}.chainR'.format(structure_file)
                pdb_chain_r.write(output_file=name_r, format='PDB')
                pdb_chain_l = SBI.structure.PDB()
                pdb_chain_l.add_chain(chain_l)
                name_l = '{}.chainL'.format(structure_file)
                pdb_chain_l.write(output_file=name_l, format='PDB')

                # We add to the list the two structures (structure_name    chain_r    chain_l)
                structure_name = os.path.basename(structure_file)
                list_fd.write('{}\t{}\t{}\n'.format( structure_name, name_r, name_l ))


    def run(self):
        '''
        Runs the computations
        '''
        # Compute the crashes if required
        if self.crashes:
            self.on_crashes()
        # Compute the potentials for each pair of structures
        with open(self.list_file, 'r') as list_fd:
            for line in list_fd:
                structure_name, chain_r, chain_l = line.strip().split('\t')
                self.on_ppi(chain_r, chain_l, fold_name=structure_name)

        sys.stdout.write('[SUCCESS] Split Potentials computed!\n')
        # Print the xml result file
        self.write_output_file(self.xml_result, self.xml_out_file)
        # Print the error xml file
        self.write_output_file(self.xml_errors, self.xml_err_file)
        sys.stdout.write('[SUCCESS] Finished!\n')

    def on_crashes(self):
        """
        Calculate the crashes between the receptor and the ligand of the PPI
            receptor_path = Complete path to the structure of the receptor
            ligand_path   = Complete path to the structure of the ligand
        """
        # $1 = SP41_path
        # $2 = distance_type
        # $3 = list_file
        # $4 = job_path
        # $5 = job_id
        # $6 = collision_path
        # $7 = collision_list_file
        #python $6/collision_detection_from_list.py $7 $4 > $4/collision.out.log 2> $4/collision.err.log | echo "$!" > $4/pid_collision.txt &
        command = 'python {} {} {} > {} 2> {} | echo "$!" > {} &'.format( os.path.join(self.crashes_path, 'collision_detection_from_list.py'), self.list_file, self.output_dir, os.path.join(self.output_dir, 'collision.out.log'), os.path.join(self.output_dir, 'collision.err.log'), os.path.join(self.output_dir, 'pid_collision.txt') )
        os.system(command)

    def on_ppi(self, receptor_path, ligand_path, fold_name=None):
        '''
        Compute the split potentials for a ppi
            receptor_path = Complete path to the structure of the receptor
            ligand_path   = Complete path to the structure of the ligand
            fold_name     = Optional parameter, name given to the PPI
        '''
        # Get PDB structure (fold structure)
        receptor = PDB.read_pdb(receptor_path)
        ligand = PDB.read_pdb(ligand_path)
        receptor_name = self.get_fold_name_from_fold_path(receptor_path)
        ligand_name = self.get_fold_name_from_fold_path(ligand_path)
        if not fold_name:
            fold_name = '{}-{}'.format(receptor_name, ligand_name)

        if isinstance(receptor, list): # if the structure obtained is an instance of a list, it means that it is not a unique structure: it is a COMPLEX of more than one chain
            receptor = PDB.read_pdb(receptor_path, merge_chains=True) # The complex will always be merged
        if isinstance(ligand, list): # if the structure obtained is an instance of a list, it means that it is not a unique structure: it is a COMPLEX of more than one chain
            ligand = PDB.read_pdb(ligand_path, merge_chains=True) # The complex will always be merged

        receptor.set_dssp()
        receptor.clean()
        receptor.normalize_residues()
        ligand.set_dssp()
        ligand.clean()
        ligand.normalize_residues()

        ppi = Interaction(receptor, ligand)

        # Compute split potentials
        if self.pot_type == 'CB':
            cutoff = 12
        else: 
            cutoff = 5

        split_potentials = SplitPotentialsPPI(c_type=self.pot_type, cutoff=cutoff) # definition of SplitPotentialsPPI class (BioLib.Docking)
        #global_energies = []
        global_energies = split_potentials.calculate_global_energies(ppi, Zscores=True) # calculation of global energies using the method of the SplitPotentialsPPI class
        print('\nGLOBAL ENERGIES')
        #residues_energies = []
        residues_energies = split_potentials.calculate_residue_energies_between_pairs(ppi, 'R', Zscores=True) # same for residues energies
        print('\nRESIDUE ENERGIES')
        self.ppi_xml(fold_name, global_energies, residues_energies, receptor, ligand) # creation of xml results file
        return

    def ppi_xml(self, fold_name, global_energies, residues_energies, receptor, ligand):
        '''
        Append the ppi energies to the xml result file
        '''
        self.xml_result += '<protein>\n'
        self.xml_result += '  <id>%s</id>\n' % fold_name
        self.xml_result += '  <length_a>%s</length_a>\n' % receptor.get_number_residues()
        self.xml_result += '  <length_b>%s</length_b>\n' % ligand.get_number_residues()
        # Global Energies (and Zenergies)
        self.xml_result += '  <global_energies>\n'
        self.xml_result += '    <Epair>%.3e</Epair>\n'     % global_energies[0]['D-PAIR']
        self.xml_result += '    <Ecomb>%.3e</Ecomb>\n'     % global_energies[0]['D-COMB']
        self.xml_result += '    <Es3dc>%.3e</Es3dc>\n'     % global_energies[0]['D-S3DC']
        self.xml_result += '    <Elocal>%.3e</Elocal>\n'   % global_energies[0]['D-LOCAL']
        self.xml_result += '    <E3dc>%.3e</E3dc>\n'       % global_energies[0]['D-3DC']
        self.xml_result += '    <E3d>%.3e</E3d>\n'         % global_energies[0]['D-3D']
        self.xml_result += '    <ZEpair>%.3e</ZEpair>\n'   % global_energies[1]['D-PAIR']
        self.xml_result += '    <ZEcomb>%.3e</ZEcomb>\n'   % global_energies[1]['D-COMB']
        self.xml_result += '    <ZEs3dc>%.3e</ZEs3dc>\n'   % global_energies[1]['D-S3DC']
        self.xml_result += '    <ZElocal>%.3e</ZElocal>\n' % global_energies[1]['D-LOCAL']
        self.xml_result += '    <ZE3dc>%.3e</ZE3dc>\n'     % global_energies[1]['D-3DC']
        self.xml_result += '    <ZE3d>NA</ZE3d>\n'
        self.xml_result += '    <Zene>%.3e</Zene>\n'       % ( float(global_energies[1]['D-S3DC']) + float(global_energies[1]['D-LOCAL']) + float(global_energies[1]['D-3DC']) )
        self.xml_result += '  </global_energies>\n'
        # Energies by residue (and Zenergies)
        self.xml_result += '  <residues>\n'
        for residue_energies in residues_energies:
            self.xml_result += '    <residue>\n'
            self.xml_result += '      <number_x>%d</number_x>\n'     % residue_energies[0].get_num()
            self.xml_result += '      <number_y>%d</number_y>\n'     % residue_energies[3].get_num()
            self.xml_result += '      <type_x>%s</type_x>\n'         % residue_energies[0].get_type_short()
            self.xml_result += '      <type_y>%s</type_y>\n'         % residue_energies[3].get_type_short()
            self.xml_result += '      <Epair>%.3e</Epair>\n'     % residue_energies[1]['D-PAIR']
            self.xml_result += '      <Ecomb>%.3e</Ecomb>\n'     % residue_energies[1]['D-COMB']
            self.xml_result += '      <Es3dc>%.3e</Es3dc>\n'     % residue_energies[1]['D-S3DC']
            self.xml_result += '      <Elocal>%.3e</Elocal>\n'   % residue_energies[1]['D-LOCAL']
            self.xml_result += '      <E3dc>%.3e</E3dc>\n'       % residue_energies[1]['D-3DC']
            self.xml_result += '      <E3d>%.3e</E3d>\n'         % residue_energies[1]['D-3D']
            self.xml_result += '      <ZEpair>%.3e</ZEpair>\n'   % residue_energies[2]['D-PAIR']
            self.xml_result += '      <ZEcomb>%.3e</ZEcomb>\n'   % residue_energies[2]['D-COMB']
            self.xml_result += '      <ZEs3dc>%.3e</ZEs3dc>\n'   % residue_energies[2]['D-S3DC']
            self.xml_result += '      <ZElocal>%.3e</ZElocal>\n' % residue_energies[2]['D-LOCAL']
            self.xml_result += '      <ZE3dc>%.3e</ZE3dc>\n'     % residue_energies[2]['D-3DC']
            self.xml_result += '      <ZE3d>NA</ZE3d>\n'
            self.xml_result += '      <Zene>%.3e</Zene>\n'       % ( float(residue_energies[2]['D-S3DC']) + float(residue_energies[2]['D-LOCAL']) + float(residue_energies[2]['D-3DC']) )
            self.xml_result += '    </residue>\n'
        self.xml_result += '  </residues>\n'
        self.xml_result += '</protein>\n'

    def print_error(self, error_num, ppi=''):
        '''
        Prints the error in the error xml file and the stderr channel
        '''
        xml_template = '<error><number>%d</number><ppi>%s</ppi><description>%s</description></error>\n'
        stderr_template = '[ERROR %d %s] %s\n'
        if error_num == 1:
            description = 'File is not compressed in tar.gz / .tgz / .zip or is not in .pdb / .cif format'
        if error_num == 2:
            description = 'The file {} does not end with .r.pdb (receptor) or .l.pdb (ligand)'.format(ppi)
        if error_num == 3:
            description = 'There must be a file {}.r.pdb and a file {}.l.pdb'.format(ppi, ppi)
        if error_num == 4:
            description = 'The number of receptors and ligands must be the same! Please, every receptor must have its corresponding ligand.'
        if error_num == 5:
            description = 'Incorrect number of chains ({}) in file {}. The PDB must only have 2 chains!'.format(ppi[0], ppi[1])
        if error_num == 6:
            description = 'Could not convert the .cif file to .pdb. Please, provide another file.'
        if error_num == 7:
            description = 'Wrong {} chain {} provided! Chains available in the structure are {}.'.format(ppi[0], ppi[1], ppi[2])
        if error_num == 8:
            description = 'File "{}" is not in any of the following formats: .pdb / .cif'.format(ppi)
        self.xml_errors += xml_template % (error_num, ppi, description)
        sys.stderr.write(stderr_template % (error_num, ppi, description))

    def write_output_file(self, xml_text, xml_file):
        """
        Prints the XML error file.
        """
        xml_fo = open(xml_file, 'w')
        xml_fo.write(xml_text+'</xml>')
        xml_fo.close()
        return

    def get_fold_name_from_fold_path(self, fold_path):
        """
        Get the name of the PDB from a path.
        """
        if '/' in fold_path:
            # Get PDB name
            fields = fold_path.split('/')
            fold_name = fields[len(fields)-1]
        else:
            fold_name = fold_path
        return fold_name

    def get_merged_chain_from_pdb(self, pdb_file, new_chain_id):
        """
        Merges the chains of the PDB file in one chain.
        Renames the chain using the input new_chain_id.
        Returns the chain with its new name
        """
        struct = SBI.structure.PDB(pdb_file)
        struct_merged = struct.fuse_chains(struct.chain_identifiers)
        for chain_id in struct_merged.chain_identifiers:
            c = struct_merged.get_chain_by_id(chain_id)
        c.chain = new_chain_id
        return c

def create_directory(directory):
    """
    Checks if a directory exists and if not, creates it
    """
    try:
        os.stat(directory)
    except:
        os.mkdir(directory)
    return


if  __name__ == "__main__":
    main()

