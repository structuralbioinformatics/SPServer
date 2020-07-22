import sys, os
import pandas as pd


def main():
    """
    From a folder of SPServer results, creates a results file for CASP.
    python /home/quim/PHD/Projects/SPServer/CASP12_analysis/scripts/analyze_casp_subset.py
    python /Users/quim/Dropbox/UPF/PHD/Projects/SPServer/CASP12_analysis/scripts/analyze_casp_subset.py
    """

    scripts_path = os.path.abspath(os.path.dirname(__file__))
    main_path = os.path.join(scripts_path, '..')
    data_dir = os.path.join(main_path, 'data')
    outputs_dir = os.path.join(main_path, 'outputs')
    create_directory(data_dir)
    create_directory(outputs_dir)
    nearnative_file = os.path.join(data_dir, 'NewNearNative.csv')
    wrong_file = os.path.join(data_dir, 'NewDecoy.csv')

    # Get near-native model names
    nearnative_df = pd.read_csv(nearnative_file, sep='\t')
    nearnative_pdbs = set(nearnative_df['PDB'])
    nearnative_structures = set([pdb.split('.pdb')[0] for pdb in nearnative_pdbs if pdb.endswith('.pdb')])

    # Get wrong model names
    wrong_df = pd.read_csv(wrong_file, sep='\t')
    wrong_pdbs = set(wrong_df['PDB'])
    wrong_structures = set([pdb.split('.pdb')[0] for pdb in wrong_pdbs if pdb.endswith('.pdb')])

    # Get target from near-native names
    target_to_nearnatives = get_target_name(nearnative_structures)

    # Get target from wrong names
    target_to_wrongs = get_target_name(wrong_structures)

    # Get all structures from folders
    structure_dirs = [f for f in os.listdir(data_dir) if os.path.isdir(os.path.join(data_dir, f))]
    for type_structure in structure_dirs:
        structure_dir = os.path.join(data_dir, type_structure)
        models = [model for model in os.listdir(structure_dir) if fileExist(os.path.join(structure_dir, model))]
        print('{} has this number of files: {}'.format(type_structure, len(models)))


    targets_with_models = set(target_to_nearnatives.keys()) | set(target_to_wrongs.keys())
    overlap_nearnative_wrong = nearnative_structures & wrong_structures
    nearnative_wrong = nearnative_structures | wrong_structures
    native_nearnative_wrong = targets_with_models | nearnative_structures | wrong_structures

    print('Number of near-native models: {}'.format(len(nearnative_structures)))
    print('Number of wrong models: {}'.format(len(wrong_structures)))
    print('Number of targets with near-native models: {}'.format(len(target_to_nearnatives)))
    print('Number of targets with wrong models: {}'.format(len(target_to_wrongs)))
    print('Number of targets with near-native or wrong models: {}'.format(len(targets_with_models)))
    print('Number of overlapped near-native and wrong models: {}'.format(len(overlap_nearnative_wrong)))
    print('Number of near-native + wrong models: {}'.format(len(nearnative_wrong)))
    print('Number of native + near-native + wrong models: {}'.format(len(native_nearnative_wrong)))

    # Write output
    output_file = os.path.join(outputs_dir, 'casp12_subset_analysis.txt')
    with open(output_file, 'w') as out_fd:
        out_fd.write('target\tnum_nearnative_decoys\tnum_wrong_decoys\n')
        for target in sorted(targets_with_models):
            target_nearnatives = set()
            target_wrongs = set()
            if target in target_to_nearnatives:
                target_nearnatives = target_to_nearnatives[target]
            if target in target_to_wrongs:
                target_wrongs = target_to_wrongs[target]
            out_fd.write('{}\t{}\t{}\n'.format(target, len(target_nearnatives), len(target_wrongs)))

    # Write output
    output_file = os.path.join(outputs_dir, 'casp12_subset.txt')
    with open(output_file, 'w') as out_fd:
        out_fd.write('target\tnear_native_models\twrong_models\n')
        for target in sorted(targets_with_models):
            if target in target_to_nearnatives and target in target_to_wrongs:
                target_nearnatives = set()
                target_wrongs = set()
                for model in target_to_nearnatives[target]:
                    # T0860D1_s236m5.pdb ==> T0860TS236_5-D1
                    model = model.replace('.pdb', '')
                    model = model[:5] + 'TS' + model[9:12] + '_' + model[-1:] + '-' + model[5:7]
                    target_nearnatives.add(model)
                for model in target_to_wrongs[target]:
                    model = model.replace('.pdb', '')
                    model = model[:5] + 'TS' + model[9:12] + '_' + model[-1:] + '-' + model[5:7]
                    target_wrongs.add(model)
                out_fd.write('{}\t{}\t{}\n'.format(target, ', '.join(sorted(target_nearnatives)), ', '.join(sorted(target_wrongs))))

    return


def fileExist(file):
    """
    Checks if a file exists AND is a file
    """
    return os.path.exists(file) and os.path.isfile(file)


def create_directory(directory):
    """
    Checks if a directory exists and if not, creates it
    """
    try:
        os.stat(directory)
    except:
        os.mkdir(directory)
    return


def get_target_name(structures):
    target_to_structures = {}
    for structure in structures:
        target = structure.split('_')[0]
        target_to_structures.setdefault(target, set()).add(structure)
    return target_to_structures

if  __name__ == "__main__":
    main()

