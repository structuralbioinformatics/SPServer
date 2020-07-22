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
    python /home/quim/PHD/Projects/CASP14/scripts/plot_global_scores.py -i -o
    python /Users/quim/Dropbox/UPF/PHD/Projects/SPServer/CASP12_analysis/scripts/plot_global_scores.py -i /Users/quim/Documents/DATA/SPServer/CASP12_analysis/results -o /Users/quim/Documents/DATA/SPServer/CASP12_analysis/results
    """
    options = parse_options()
    plot_global_scores(options)


def parse_options():
    """
    This function parses the command line arguments and returns an optparse object.
    """

    parser = optparse.OptionParser("plot_global_scores.py  -i INPUT_DIR -o OUTPUT_DIR")

    # Directory arguments
    parser.add_option("-i", action="store", type="string", dest="input_dir", help="Folder with the SPServer results that we want to analyze", metavar="INPUT_DIR")
    parser.add_option("-o", action="store", type="string", dest="output_dir", help="Output folder where to store the results", metavar="OUTPUT_DIR")

    (options, args) = parser.parse_args()

    return options


def plot_global_scores(options):
    """
    From a folder of SPServer results, creates a results file for CASP.
    """
    # Get parameters
    input_dir = options.input_dir
    output_dir = options.output_dir
    create_directory(output_dir)
    plots_dir = os.path.join(output_dir, 'plots_global_scores')
    create_directory(plots_dir)


    # Load results of global scores
    results_file = os.path.join(input_dir, 'CASP12_global_results.txt')
    results_df = pd.read_csv(results_file, sep='\t', index_col=None)
    results_df['Color'] = 'red'
    results_df.loc[results_df['TypeStructure'] == 'Native', 'Color'] = '#24c729'
    results_df.loc[results_df['TypeStructure'] == 'Near-native', 'Color'] = 'blue'
    results_df.loc[results_df['TypeStructure'] == 'Wrong', 'Color'] = 'red'
    metric_to_top = {'GDT_TS':100, 'TM score':1, 'QCS':100}
    metric_to_bottom = {'GDT_TS':0, 'TM score':0, 'QCS':0}


    # Plot results in a scatter plot
    for alternative_scoring_function in ['PROSA', 'DOPE', 'GDT_TS', 'TM score', 'QCS']:
        for spserver_scoring_function in ['ZES3DC', 'ZPAIR']:

            output_plot = os.path.join(plots_dir, '{}_vs_{}.png'.format(spserver_scoring_function, alternative_scoring_function))
            if not fileExist(output_plot):

                # Calculate pearson correlation
                scores_spserver = results_df[spserver_scoring_function].astype(float).tolist()
                scores_alternative = results_df[alternative_scoring_function].astype(float).tolist()
                r_value, p_value = pearsonr(scores_spserver, scores_alternative)
                r_value = round(r_value, 2)
                print ('Correlation for energies {} and {}: {}'.format(spserver_scoring_function, alternative_scoring_function, r_value))

                # Make scatter plot
                fig = pylab.figure(dpi=300)
                ax = sns.regplot(results_df[spserver_scoring_function], results_df[alternative_scoring_function], color="black", order=1, marker="+", label="R: {}".format(r_value), scatter_kws={'facecolors':results_df['Color']})
                if alternative_scoring_function in ['GDT_TS', 'TM score', 'QCS']:
                    ax.set_ylim(metric_to_bottom[alternative_scoring_function], metric_to_top[alternative_scoring_function])
                # Make a legend for the scatter plot
                # groupby and plot points of one color
                results_df[results_df['Color'] == 'red'].plot(kind = 'scatter', x=spserver_scoring_function, y=alternative_scoring_function, c = 'red', ax = ax, marker="+", label = 'Wrong', zorder = 3)
                results_df[results_df['Color'] == 'blue'].plot(kind = 'scatter', x=spserver_scoring_function, y=alternative_scoring_function, c = 'blue', ax = ax, marker="+", label = 'Near-native', zorder = 3)
                results_df[results_df['Color'] == '#24c729'].plot(kind = 'scatter', x=spserver_scoring_function, y=alternative_scoring_function, c = '#24c729', ax = ax, marker="+", label = 'Native', zorder = 3)
                handles, labels = ax.get_legend_handles_labels() # How to remove the default title: https://stackoverflow.com/questions/51579215/remove-seaborn-lineplot-legend-title?rq=1
                ax.legend(handles=plt.plot([],ls="-", color='black') + handles[1:], labels=labels, loc="best")
                plt.gray()
                plt.grid(True, color='w', linestyle='-', linewidth=1, zorder=0)
                plt.gca().patch.set_facecolor('#EAEAF1')
                plt.xlabel(spserver_scoring_function, fontweight='bold')
                plt.ylabel(alternative_scoring_function, fontweight='bold')

                pylab.savefig(output_plot, format='png')
                plt.clf()


    # Plot results of score vs length
    sns.set(color_codes=True)
    sns.set_style(style='white')
    for scoring_function in ['ZES3DC', 'ZPAIR', 'PROSA', 'DOPE']:
        score_results_df = results_df[['TypeStructure', scoring_function, 'Length']]
        output_plot = os.path.join(plots_dir, '{}_vs_length.png'.format(scoring_function))
        fig = pylab.figure(dpi=500)
        ax = sns.lmplot(x='Length', 
                        y=scoring_function ,
                        hue='TypeStructure', 
                        hue_order=["Wrong", "Near-native","Native"],
                        data=score_results_df, 
                        x_jitter= 5,
                        legend="full",
                        palette="Set1")
        plt.gray()
        plt.grid(True, color='w', linestyle='-', linewidth=1)
        plt.gca().patch.set_facecolor('#EAEAF1')
        pylab.savefig(output_plot, format='png')
        plt.clf()

        output_plot = os.path.join(plots_dir, '{}_vs_length_histogram.png'.format(scoring_function))
        fig = pylab.figure(dpi=500)
        sns.kdeplot(results_df.loc[results_df['TypeStructure'] == "Native"][scoring_function], shade=True, color="#58D68D", legend=True);
        sns.kdeplot(results_df.loc[results_df['TypeStructure'] == "Near-native"][scoring_function], shade=True, color="#5DADE2", legend=True);
        sns.kdeplot(results_df.loc[results_df['TypeStructure'] == "Wrong"][scoring_function], shade=True, color="#E74C3C", legend=True);
        sns.despine(fig=None, ax=None, top=True, right=True, left=True, bottom=False, offset=None, trim=False)
        pylab.savefig(output_plot, format='png')
        plt.clf()


    return



#######################
#######################
# SECONDARY FUNCTIONS #
#######################
#######################


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
