import sys, os
import cPickle
import matplotlib.pyplot as plt
import numpy as np
import optparse
import pandas as pd
import pylab
from scipy.stats import pearsonr
import seaborn as sns
import xml.etree.ElementTree as ET


def main():
    """
    From a folder of SPServer results, creates a results file for CASP.
    python /home/quim/PHD/Projects/CASP14/scripts/plot_residue_scores.py -i -o
    python /Users/quim/Dropbox/UPF/PHD/Projects/SPServer/CASP12_analysis/scripts/plot_residue_scores.py -i /Users/quim/Documents/DATA/SPServer/CASP12_analysis/outputs -o /Users/quim/Documents/DATA/SPServer/CASP12_analysis/results
    """
    options = parse_options()
    plot_residue_scores(options)


def parse_options():
    """
    This function parses the command line arguments and returns an optparse object.
    """

    parser = optparse.OptionParser("plot_residue_scores.py  -i INPUT_DIR -o OUTPUT_DIR")

    # Directory arguments
    parser.add_option("-i", action="store", type="string", dest="input_dir", help="Folder with the SPServer results that we want to analyze", metavar="INPUT_DIR")
    parser.add_option("-o", action="store", type="string", dest="output_dir", help="Output folder where to store the results", metavar="OUTPUT_DIR")

    (options, args) = parser.parse_args()

    return options


def plot_residue_scores(options):
    """
    From a folder of SPServer results, creates a results file for CASP.
    """
    # Get parameters
    input_dir = options.input_dir
    output_dir = options.output_dir
    create_directory(output_dir)
    plots_dir = os.path.join(output_dir, 'plots_residue_scores')
    create_directory(plots_dir)

    # Get directories
    spserver_input_dir = os.path.join(input_dir, 'SPServer_results_residues')
    dope_input_dir = os.path.join(input_dir, 'DOPE')
    prose_input_dir = os.path.join(input_dir, 'PROSA_residues')
    residues_input_dir = os.path.join(input_dir, 'residue_score_results')
    create_directory(residues_input_dir)

    # Load global results to get all structures
    results_file = os.path.join(output_dir, 'CASP12_global_results.txt')
    results_df = pd.read_csv(results_file, sep='\t', index_col=None)
    structures = sorted(set(results_df['PDB']))

    # Parse residue scores for each model
    output_correlations_file = os.path.join(output_dir, 'CASP12_residue_correlations.txt')
    columns = ['ScoringFunction1', 'ScoringFunction2', 'Target', 'PDB', 'TypeStructure', 'PearsonCorrelation', 'P-value']
    correlations_df = pd.DataFrame(columns=columns)
    spserver_scoring_functions = ['ES3DC', 'PAIR']
    alternative_scoring_functions = ['DOPE', 'PROSA']
    scoring_function_to_scores = {}
    scoring_function_to_types = {}
    score1_to_score2_to_rvalues = {}
    scoring_function_to_scores_file = os.path.join(output_dir, 'scoring_function_to_scores_residues.pcl')
    scoring_function_to_types_file = os.path.join(output_dir, 'scoring_function_to_types_residues.pcl')
    score1_to_score2_to_rvalues_file = os.path.join(output_dir, 'score1_to_score2_to_rvalues_residues.pcl')
    if not fileExist(scoring_function_to_scores_file) or not fileExist(score1_to_score2_to_rvalues_file):

        for structure in structures:

            print(structure)
            structure_name = structure.split('.pdb')[0]
            target = results_df.loc[results_df['PDB'] == structure, 'Target'].values.flatten().tolist()[0]
            type_structure = results_df.loc[results_df['PDB'] == structure, 'TypeStructure'].values.flatten().tolist()[0]

            # Required files
            residues_file = os.path.join(residues_input_dir, '{}.residues.txt'.format(structure))
            spserver_file = os.path.join(spserver_input_dir, '{}_SPSresidues.csv'.format(structure)) # e.g. T0860D1_s005m1.pdb_SPSresidues.csv
            dope_file = os.path.join(dope_input_dir, '{}.profile'.format(structure)) # e.g. T0860D1_s005m1.pdb.profile
            prosa_file = os.path.join(prose_input_dir, '{}_PROSARESIDUE.csv'.format(structure_name)) # e.g. T0860D1_PROSARESIDUE.csv

            # If the residue scores for SPServer, DOPE and PROSA have not been parsed, then we parse them
            if not fileExist(residues_file):
                if fileExist(spserver_file) and fileExist(dope_file) and fileExist(prosa_file):

                    # Parse SPServer residue scores
                    spserver_df = pd.read_csv(spserver_file, sep='\t', index_col=0)

                    # Parse DOPE
                    dope_df = read_DOPE_residues(dope_file)

                    # Parse PROSA
                    prosa_df = pd.read_csv(prosa_file, sep='\t', index_col=0)
                    prosa_df = prosa_df.rename(columns={'ResPROSA': 'PROSA'})

                    residues_df = pd.concat([dope_df['Residue'], spserver_df['ES3DC'], spserver_df['PAIR'], dope_df['DOPE'], prosa_df['PROSA']], axis=1)
                    residues_df.to_csv(residues_file, sep='\t', index=False)

                else:
                    if not fileExist(spserver_file):
                        print('FILE NOT FOUND FOR STRUCTURE {}: {}'.format(structure, spserver_file))
                    elif not fileExist(dope_file):
                        print('FILE NOT FOUND FOR STRUCTURE {}: {}'.format(structure, dope_file))
                    elif not fileExist(prosa_file):
                        print('FILE NOT FOUND FOR STRUCTURE {}: {}'.format(structure, prosa_file))
                    continue
            else:
                residues_df = pd.read_csv(residues_file, sep='\t', index_col=None)

            # We calculate the correlations between pairs of scoring functions for the residue scores of the structure
            for scoring_function1 in spserver_scoring_functions:
                for scoring_function2 in alternative_scoring_functions:

                    # Save scores
                    scoring_function_to_scores.setdefault(scoring_function1, [])
                    scoring_function_to_scores[scoring_function1] += residues_df[scoring_function1].tolist()
                    scoring_function_to_scores.setdefault(scoring_function2, [])
                    scoring_function_to_scores[scoring_function2] += residues_df[scoring_function2].tolist()
                    # Save types
                    scoring_function_to_types.setdefault(scoring_function1, [])
                    scoring_function_to_types[scoring_function1] += [type_structure]*len(residues_df[scoring_function1].tolist())
                    scoring_function_to_types.setdefault(scoring_function2, [])
                    scoring_function_to_types[scoring_function2] += [type_structure]*len(residues_df[scoring_function2].tolist())

                    # Calculate correlation between residue scores
                    r_value, p_value = pearsonr(residues_df[scoring_function1], residues_df[scoring_function2])
                    score1_to_score2_to_rvalues.setdefault(scoring_function1, {})
                    score1_to_score2_to_rvalues[scoring_function1].setdefault(scoring_function2, []).append(r_value)
                    results = [scoring_function1, scoring_function2, target, structure, type_structure, r_value, p_value]
                    df2 = pd.DataFrame([results], columns=columns)
                    correlations_df = correlations_df.append(df2)


        # Write results
        correlations_df.to_csv(output_correlations_file, sep='\t', index=False)
        correlations_df = pd.read_csv(output_correlations_file, sep='\t', index_col=None)
        #cPickle.dump(scoring_function_to_scores, open(scoring_function_to_scores_file, 'w'))
        #cPickle.dump(scoring_function_to_types, open(scoring_function_to_types_file, 'w'))
        #cPickle.dump(score1_to_score2_to_rvalues, open(score1_to_score2_to_rvalues_file, 'w'))
    else:
        # Load results
        correlations_df = pd.read_csv(output_correlations_file, sep='\t', index_col=None)
        scoring_function_to_scores = cPickle.load(open(scoring_function_to_scores_file))
        scoring_function_to_types = cPickle.load(open(scoring_function_to_types_file))
        score1_to_score2_to_rvalues = cPickle.load(open(score1_to_score2_to_rvalues_file))


    # Make plots and table of correlations
    output_nice_mean_correlations_file = os.path.join(output_dir, 'CASP12_residue_nice_mean_correlations.txt')
    columns = ['ScoringFunction1']
    for scoring_function in alternative_scoring_functions:
        columns = columns + [scoring_function+'-MEAN', scoring_function+'-SD']
    nice_mean_correlations_df = pd.DataFrame(columns=columns)

    for scoring_function1 in spserver_scoring_functions:
        results = [scoring_function1]
        for scoring_function2 in alternative_scoring_functions:

            print('Plotting: {} vs {}'.format(scoring_function1, scoring_function2))

            # #### Make scatter plot (not available because there are too many residues!) ####
            # output_plot = os.path.join(plots_dir, 'residue_scores_scatter_{}_vs_{}.png'.format(scoring_function1, scoring_function2))
            # fig = pylab.figure(dpi=300)
            # # Calculate correlation of all residues together
            # r_value, p_value = pearsonr(scoring_function_to_scores[scoring_function1], scoring_function_to_scores[scoring_function2])
            # r_value = round(r_value, 2)
            # all_scores_df = pd.DataFrame()
            # all_scores_df[scoring_function1] = scoring_function_to_scores[scoring_function1]
            # all_scores_df[scoring_function2] = scoring_function_to_scores[scoring_function2]
            # all_scores_df['TypeStructure'] = scoring_function_to_types[scoring_function1]
            # all_scores_df['Color'] = 'red'
            # all_scores_df.loc[all_scores_df['TypeStructure'] == 'Native', 'Color'] = '#24c729'
            # all_scores_df.loc[all_scores_df['TypeStructure'] == 'Near-native', 'Color'] = 'blue'
            # all_scores_df.loc[all_scores_df['TypeStructure'] == 'Wrong', 'Color'] = 'red'
            # ax = sns.regplot(all_scores_df[scoring_function1], all_scores_df[scoring_function2], color="black", order=1, marker="+", label="R: {}".format(r_value), scatter_kws={'facecolors':all_scores_df['Color']})
            # # Make a legend
            # # groupby and plot points of one color
            # all_scores_df[all_scores_df['Color'] == 'red'].plot(kind = 'scatter', x=scoring_function1, y=scoring_function2, c = 'red', ax = ax, marker="+", label = 'Wrong', zorder = 3)
            # all_scores_df[all_scores_df['Color'] == 'blue'].plot(kind = 'scatter', x=scoring_function1, y=scoring_function2, c = 'blue', ax = ax, marker="+", label = 'Near-native', zorder = 3)
            # all_scores_df[all_scores_df['Color'] == '#24c729'].plot(kind = 'scatter', x=scoring_function1, y=scoring_function2, c = '#24c729', ax = ax, marker="+", label = 'Native', zorder = 3)
            # handles, labels = ax.get_legend_handles_labels() # How to remove the default title: https://stackoverflow.com/questions/51579215/remove-seaborn-lineplot-legend-title?rq=1
            # ax.legend(handles=plt.plot([],ls="-", color='black') + handles[1:], labels=labels, loc="best")
            # plt.gray()
            # plt.grid(True, color='w', linestyle='-', linewidth=1, zorder=0)
            # plt.gca().patch.set_facecolor('#EAEAF1')
            # plt.xlabel(scoring_function1, fontweight='bold')
            # plt.ylabel(scoring_function2, fontweight='bold')
            # pylab.savefig(output_plot, format='png')
            # plt.clf()

            #### Make histogram plot ####
            output_plot = os.path.join(plots_dir, 'residue_scores_histogram_{}_vs_{}.png'.format(scoring_function1, scoring_function2))
            fig = pylab.figure(dpi=300)
            sns.distplot(score1_to_score2_to_rvalues[scoring_function1][scoring_function2], rug=True)
            plt.gray()
            plt.grid(True, color='w', linestyle='-', linewidth=1)
            plt.gca().patch.set_facecolor('#EAEAF1')
            plt.xlabel('Pearson correlation values', fontweight='bold')
            plt.ylabel('Percentage of observations', fontweight='bold')
            pylab.savefig(output_plot, format='png')
            plt.clf()

            correlations = score1_to_score2_to_rvalues[scoring_function1][scoring_function2]
            mean = np.mean(correlations)
            sdev = np.std(correlations)
            results = results + ['{:.2f}'.format(mean), '{:.2f}'.format(sdev)]
            print(scoring_function1, scoring_function2, mean, sdev)
        df2 = pd.DataFrame([results], columns=columns)
        nice_mean_correlations_df = nice_mean_correlations_df.append(df2)
    nice_mean_correlations_df.to_csv(output_nice_mean_correlations_file, sep='\t', index=False)


    return



#######################
#######################
# SECONDARY FUNCTIONS #
#######################
#######################


def read_DOPE_residues(DOPE_file):
    """
    Reads the DOPE residues
    """
    dope_df = pd.DataFrame()
    residues = []
    scores = []
    with open(DOPE_file, 'r') as input_fd:
        for line in input_fd:
            if line[0] == '#':
                continue
            fields = line.split()
            if len(fields) > 0:
                residue = int(fields[0])
                score = float(fields[-1])
                residues.append(residue)
                scores.append(score)
    dope_df['Residue'] = residues
    dope_df['DOPE'] = scores
    return dope_df


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
