import sys, os
import matplotlib.pyplot as plt
import numpy as np
import optparse
import pandas as pd
import pylab
import random
from scipy.stats import pearsonr
import seaborn as sns



def main():
    """
    From a folder of SPServer results, creates a results file for CASP.
    python /home/quim/PHD/Projects/CASP14/scripts/make_bootstrapping.py -i -d -o
    python /Users/quim/Dropbox/UPF/PHD/Projects/SPServer/SPServer/CASP12_analysis/scripts/make_bootstrapping.py -i /Users/quim/Documents/DATA/SPServer/CASP12_analysis/results -d /Users/quim/Documents/DATA/SPServer/CASP12_analysis/data -o /Users/quim/Documents/DATA/SPServer/CASP12_analysis/results
    """
    options = parse_options()
    make_bootstrapping(options)


def parse_options():
    """
    This function parses the command line arguments and returns an optparse object.
    """

    parser = optparse.OptionParser("make_bootstrapping.py  -i INPUT_DIR -d DATA_DIR -o OUTPUT_DIR")

    # Directory arguments
    parser.add_option("-i", action="store", type="string", dest="input_dir", help="Folder with the SPServer results that we want to analyze", metavar="INPUT_DIR")
    parser.add_option("-d", action="store", type="string", dest="data_dir", help="Folder with the SPServer data (tables with GDT_TS, QCS, TM score)", metavar="DATA_DIR")
    parser.add_option("-o", action="store", type="string", dest="output_dir", help="Output folder where to store the results", metavar="OUTPUT_DIR")

    (options, args) = parser.parse_args()

    return options


def make_bootstrapping(options):
    """
    From a folder of SPServer results, creates a results file for CASP.
    """
    # Get parameters
    input_dir = options.input_dir
    data_dir = options.data_dir
    output_dir = options.output_dir
    create_directory(output_dir)
    num_repetitions = 1000

    # Load global score results
    results_file = os.path.join(output_dir, 'CASP12_global_results.txt')
    results_df = pd.read_csv(results_file, sep='\t', index_col=None)

    # Check number of native, near-native and wrong structures
    pdbs_native = set(results_df.loc[results_df['TypeStructure'] == 'Native', 'PDB'])
    pdbs_nearnative = set(results_df.loc[results_df['TypeStructure'] == 'Near-native', 'PDB'])
    pdbs_wrong = set(results_df.loc[results_df['TypeStructure'] == 'Wrong', 'PDB'])
    print('Number of native: {}'.format(len(pdbs_native)))
    print('Number of near-native: {}'.format(len(pdbs_nearnative)))
    print('Number of wrong: {}'.format(len(pdbs_wrong)))
    print('Number of structures: {}'.format(len(results_df.index)))

    # We get the targets that have at least 1 model of each type
    targets_with_native = set(results_df.loc[results_df['TypeStructure'] == 'Native', 'Target'])
    targets_with_nearnative = set(results_df.loc[results_df['TypeStructure'] == 'Near-native', 'Target'])
    targets_with_wrong = set(results_df.loc[results_df['TypeStructure'] == 'Wrong', 'Target'])
    targets_all_types = sorted(list(targets_with_native & targets_with_nearnative & targets_with_wrong))

    # We get all the structures for these targets
    results_filtered_df = results_df[results_df['Target'].isin(targets_all_types)]

    # Calculate bootstrapping
    output_bootstrap_file = os.path.join(output_dir, 'CASP12_bootstrapping.txt')
    if not fileExist(output_bootstrap_file):

        target_to_nearnatives_done = {}
        target_to_wrongs_done = {}
        columns = ['Repetition', 'Target', 'PDB', 'TypeStructure', 'ZES3DC', 'ZPAIR', 'PROSA', 'DOPE', 'GDT_TS', 'TM score', 'QCS']
        bootstrapping_df = pd.DataFrame(columns=columns)

        for i in xrange(num_repetitions):

            print(i+1)
            structures = []
            scores = []

            for target in sorted(targets_all_types):

                # Get target models
                target_df = results_filtered_df[results_filtered_df['Target'] == target]
                native = list(set(target_df.loc[target_df['TypeStructure'] == 'Native', 'PDB']))[0]
                nearnatives = set(target_df.loc[target_df['TypeStructure'] == 'Near-native', 'PDB'])
                wrongs = set(target_df.loc[target_df['TypeStructure'] == 'Wrong', 'PDB'])

                # Choose models from the ones that have not been selected in previous repetitions
                nearnative, target_to_nearnatives_done = choose_structure(nearnatives, target_to_nearnatives_done, target)
                wrong, target_to_wrongs_done = choose_structure(wrongs, target_to_wrongs_done, target)

                #print(native, nearnative, wrong)

                # Add the structures to the dataframe
                for model in [native, nearnative, wrong]:
                    #model_data = results_filtered_df[results_filtered_df['PDB'] == model][['TypeStructure', 'ZES3DC', 'ZPAIR', 'PROSA', 'DOPE', 'GDT_TS', 'TM score', 'QCS']].tolist()
                    model_data = results_filtered_df.loc[results_filtered_df['PDB'] == model, ['TypeStructure', 'ZES3DC', 'ZPAIR', 'PROSA', 'DOPE', 'GDT_TS', 'TM score', 'QCS']].values.flatten().tolist()
                    results = [i+1, target, model] + model_data
                    df2 = pd.DataFrame([results], columns=columns, index=[i])
                    bootstrapping_df = bootstrapping_df.append(df2) # Add the information to the main data frame

        # Write results 
        bootstrapping_df.to_csv(output_bootstrap_file, sep='\t', index=False)
        bootstrapping_df = pd.read_csv(output_bootstrap_file, sep='\t', index_col=None)
    else:
        # Load results
        bootstrapping_df = pd.read_csv(output_bootstrap_file, sep='\t', index_col=None)


    # Calculate the correlations for each repetition and pair of scoring functions
    output_correlations_file = os.path.join(output_dir, 'CASP12_correlations.txt')
    if not fileExist(output_correlations_file):
        columns = ['Repetition', 'ScoringFunction1', 'ScoringFunction2', 'PearsonCorrelation', 'P-value']
        correlations_df = pd.DataFrame(columns=columns)
        scoring_functions = ['ZES3DC', 'ZPAIR', 'PROSA', 'DOPE', 'GDT_TS', 'TM score', 'QCS']
        for i in xrange(num_repetitions):
            repetition = i+1
            print(repetition)
            repetition_df = bootstrapping_df[bootstrapping_df['Repetition'] == repetition]
            for scoring_function1 in scoring_functions:
                for scoring_function2 in scoring_functions:
                    scores1 = repetition_df[scoring_function1].astype(float).tolist()
                    scores2 = repetition_df[scoring_function2].astype(float).tolist()
                    r_value, p_value = pearsonr(scores1, scores2)
                    results = [repetition, scoring_function1, scoring_function2, r_value, p_value]
                    df2 = pd.DataFrame([results], columns=columns, index=[i])
                    correlations_df = correlations_df.append(df2)
        # Write results
        correlations_df.to_csv(output_correlations_file, sep='\t', index=False)
        correlations_df = pd.read_csv(output_correlations_file, sep='\t', index_col=None)
    else:
        # Load results
        correlations_df = pd.read_csv(output_correlations_file, sep='\t', index_col=None)


    # Calculate mean of correlations and standard deviation
    output_mean_correlations_file = os.path.join(output_dir, 'CASP12_mean_correlations.txt')
    scoring_functions = ['ZES3DC', 'ZPAIR', 'PROSA', 'DOPE', 'GDT_TS', 'TM score', 'QCS']
    if not fileExist(output_mean_correlations_file):
        columns = ['ScoringFunction1', 'ScoringFunction2', 'MeanCorrelation', 'StandardDeviation']
        mean_correlations_df = pd.DataFrame(columns=columns)
        for scoring_function1 in scoring_functions:
            for scoring_function2 in scoring_functions:
                correlations = correlations_df.loc[ (correlations_df['ScoringFunction1'] == scoring_function1) & (correlations_df['ScoringFunction2'] == scoring_function2) , 'PearsonCorrelation'].astype(float).tolist()
                mean = np.mean(correlations)
                standard_deviation = np.std(correlations)
                results = [scoring_function1, scoring_function2, mean, standard_deviation]
                df2 = pd.DataFrame([results], columns=columns)
                mean_correlations_df = mean_correlations_df.append(df2)
        # Write results
        mean_correlations_df.to_csv(output_mean_correlations_file, sep='\t', index=False)
        mean_correlations_df = pd.read_csv(output_mean_correlations_file, sep='\t', index_col=None)
    else:
        # Load results
        mean_correlations_df = pd.read_csv(output_mean_correlations_file, sep='\t', index_col=None)


    # Make table of mean correlations nice
    output_nice_mean_correlations_file = os.path.join(output_dir, 'CASP12_nice_mean_correlations.txt')
    columns = ['ScoringFunction1']
    for scoring_function in scoring_functions:
        columns = columns + [scoring_function+'-MEAN', scoring_function+'-SD']
    nice_mean_correlations_df = pd.DataFrame(columns=columns)
    for scoring_function1 in scoring_functions:
        results = [scoring_function1]
        for scoring_function2 in scoring_functions:
            mean = mean_correlations_df.loc[ (mean_correlations_df['ScoringFunction1'] == scoring_function1) & (correlations_df['ScoringFunction2'] == scoring_function2) , 'MeanCorrelation'].astype(float).tolist()[0]
            sdev = mean_correlations_df.loc[ (mean_correlations_df['ScoringFunction1'] == scoring_function1) & (correlations_df['ScoringFunction2'] == scoring_function2) , 'StandardDeviation'].astype(float).tolist()[0]
            results = results + ['{:.2f}'.format(mean), '{:.2f}'.format(sdev)]
        df2 = pd.DataFrame([results], columns=columns)
        nice_mean_correlations_df = nice_mean_correlations_df.append(df2)
    nice_mean_correlations_df.to_csv(output_nice_mean_correlations_file, sep='\t', index=False)


    # Plot heatmap
    means = mean_correlations_df.pivot("ScoringFunction1", "ScoringFunction2", "MeanCorrelation")
    plot_name = os.path.join(output_dir, 'heatmap_correlations.png')
    fig = pylab.figure(dpi=300)
    cmap = sns.diverging_palette(240, 10, n=100, as_cmap=True)
    cmap.set_under(".5")
    #mask = np.ones_like(means)
    #mask[np.tril_indices_from(mask)] = False
    #ax = sns.heatmap(means, mask=mask, annot=True, cmap=cmap, center=0, cbar_kws={'label': 'Mean Pearson Correlation'}, vmin=-1, vmax=1)
    ax = sns.heatmap(means, annot=True, cmap=cmap, center=0, cbar_kws={'label': 'Mean Pearson Correlation'}, vmin=-1, vmax=1, fmt='.2f')
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
             rotation_mode="anchor")
    plt.yticks(rotation=0)
    plt.xlabel('')
    plt.ylabel('')
    fig.tight_layout()
    pylab.savefig(plot_name, format='png')

    return



#######################
#######################
# SECONDARY FUNCTIONS #
#######################
#######################

def choose_structure(models, target_to_models_done, target):
    """
    Select one model among the ones that have not been selected yet.
    If all of them have been selected, start again from 0.
    """
    target_to_models_done.setdefault(target, set())
    models_todo = models - target_to_models_done[target]
    if len(models_todo) == 0:
        #print('All done: {}'.format(target_to_models_done[target]))
        target_to_models_done[target] = set()
        model_selected = random.choice(list(models))
    else:
        model_selected = random.choice(list(models_todo))
    target_to_models_done[target].add(model_selected)
    return model_selected, target_to_models_done


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


if  __name__ == "__main__":
    main()
