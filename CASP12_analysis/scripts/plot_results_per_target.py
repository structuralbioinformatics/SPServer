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
    python /home/quim/PHD/Projects/CASP14/scripts/plot_results_per_target.py -i -o
    python /Users/quim/Dropbox/UPF/PHD/Projects/SPServer/CASP12_analysis/scripts/plot_results_per_target.py -i /Users/quim/Dropbox/UPF/PHD/Projects/SPServer/CASP12_analysis/outputs/all_results -o /Users/quim/Dropbox/UPF/PHD/Projects/SPServer/CASP12_analysis/results
    """
    options = parse_options()
    plot_results_per_target(options)


def parse_options():
    """
    This function parses the command line arguments and returns an optparse object.
    """

    parser = optparse.OptionParser("plot_results_per_target.py  -i INPUT_DIR -o OUTPUT_DIR")

    # Directory arguments
    parser.add_option("-i", action="store", type="string", dest="input_dir", help="Folder with the SPServer results that we want to analyze", metavar="INPUT_DIR")
    parser.add_option("-o", action="store", type="string", dest="output_dir", help="Output folder where to store the results", metavar="OUTPUT_DIR")

    (options, args) = parser.parse_args()

    return options


def plot_results_per_target(options):
    """
    From a folder of SPServer results, creates a results file for CASP.
    """
    # Get parameters
    input_dir = options.input_dir
    output_dir = options.output_dir
    create_directory(output_dir)
    scoring_functions = ['ES3DC', 'PAIR', 'ZES3DC', 'ZPAIR']

    # Read results
    native_results_file = os.path.join(input_dir, 'Native.csv')
    nearnative_results_file = os.path.join(input_dir, 'NearNative.csv')
    wrong_results_file = os.path.join(input_dir, 'Decoy.csv')
    native_results_df = pd.read_csv(native_results_file, sep='\t', index_col=0)
    nearnative_results_df = pd.read_csv(nearnative_results_file, sep='\t', index_col=0)
    wrong_results_df = pd.read_csv(wrong_results_file, sep='\t', index_col=0)
    nearnative_results_df['TypeStructure'] = 'Near-native'
    nearnative_results_df['Color'] = 'blue'
    wrong_results_df['TypeStructure'] = 'Wrong'
    wrong_results_df['Color'] = 'red'
    results_df = pd.concat([nearnative_results_df, wrong_results_df]).reset_index(drop=True)
    new_df = results_df['PDB'].str.split('_', n = 1, expand = True)
    results_df['Target'] = new_df[0]
    targets = set(results_df['Target'])
    results_df.rename(columns={ 'SPS_es3dc':'ES3DC',
                                'SPS_Pair':'PAIR',
                                'SPS_zes3dc':'ZES3DC',
                                'SPS_Zpair':'ZPAIR',
                                'GDT':'GDT_TS',
                                'TM':'TM score'}, 
                    inplace=True)


    # Plot results per target
    results_per_target_dir = os.path.join(output_dir, 'results_per_target')
    create_directory(results_per_target_dir)
    results_per_target_file = os.path.join(output_dir, 'results_per_target.txt')
    columns = ['Target', 'ZES3DC-GDT', 'ZES3DC-TM', 'ZES3DC-QCS', 'ZPAIR-GDT', 'ZPAIR-TM', 'ZPAIR-QCS']
    results_per_target_df = pd.DataFrame(columns=columns)
    metric_to_top = {'GDT_TS':100, 'TM score':1, 'QCS':100}
    metric_to_bottom = {'GDT_TS':0, 'TM score':0, 'QCS':0}

    for target in sorted(targets):
        target_df = results_df[results_df['Target'] == target].fillna(0)
        r_values = []
        for spserver_scoring_function in ['ZES3DC', 'ZPAIR']:
            for metric in ['GDT_TS', 'TM score', 'QCS']:
            #for metric in ['GDT_TS']:

                # Calculate pearson correlation
                scores_spserver = target_df[spserver_scoring_function].astype(float).tolist()
                scores_metric = target_df[metric].astype(float).tolist()
                r_value, p_value = pearsonr(scores_spserver, scores_metric)
                r_value = round(r_value, 2)
                r_values.append(r_value)
                print ('Correlation for energies {} and {}: {}'.format(spserver_scoring_function, metric, r_value))

                # Make plot
                output_plot = os.path.join(results_per_target_dir, '{}_{}_vs_{}.png'.format(target, spserver_scoring_function, metric))
                fig = pylab.figure(dpi=300)
                #sns.set_context("paper")
                #ax = sns.lmplot(x=spserver_scoring_function, y=metric, data=target_df, hue="TypeStructure", palette=['blue', 'red'], markers="+")
                ax = sns.regplot(target_df[spserver_scoring_function], target_df[metric], color="black", order=1, marker="+", label="R: {}".format(r_value), scatter_kws={'facecolors':target_df['Color']})
                ax.set_ylim(metric_to_bottom[metric], metric_to_top[metric])

                # Make a legend
                # groupby and plot points of one color
                target_df[target_df['Color'] == 'blue'].plot(kind = 'scatter', x=spserver_scoring_function, y=metric, c = 'blue', ax = ax, marker="+", label = 'Near-native', zorder = 3)
                target_df[target_df['Color'] == 'red'].plot(kind = 'scatter', x=spserver_scoring_function, y=metric, c = 'red', ax = ax, marker="+", label = 'Wrong', zorder = 3)
                handles, labels = ax.get_legend_handles_labels() # How to remove the default title: https://stackoverflow.com/questions/51579215/remove-seaborn-lineplot-legend-title?rq=1
                ax.legend(handles=plt.plot([],ls="-", color='black') + handles[1:], labels=labels, loc="best")
                #plt.title('{} vs {}'.format(spserver_scoring_function, metric), fontweight='bold')
                plt.gray()
                plt.grid(True, color='w', linestyle='-', linewidth=1, zorder=0)
                plt.gca().patch.set_facecolor('#EAEAF1')
                plt.xlabel(spserver_scoring_function, fontweight='bold')
                plt.ylabel(metric, fontweight='bold')
                
                pylab.savefig(output_plot, format='png')
                plt.clf()

        # Save the r values in the table
        results = [target] + r_values
        df2 = pd.DataFrame([results], columns=columns)
        results_per_target_df = results_per_target_df.append(df2)

        # Plot pair plot
        # output_plot = os.path.join(results_per_target_dir, '{}.png'.format(target))
        # fig, axs = pylab.subplots(ncols=2, dpi=300, figsize=(9, 4))
        # plt.rcParams['axes.grid'] = True
        # sns.regplot(target_df['ZES3DC'], target_df['GDT_TS'], color="black", order=1, marker="+", label="R: {}".format(r_value), scatter_kws={'facecolors':target_df['Color']}, ax=axs[0])
        # sns.regplot(target_df['ZPAIR'], target_df['GDT_TS'], color="black", order=1, marker="+", label="R: {}".format(r_value), scatter_kws={'facecolors':target_df['Color']}, ax=axs[1])
        # target_df[target_df['Color'] == 'blue'].plot(kind = 'scatter', x="ZES3DC", y="GDT_TS", c = 'blue', ax = axs[0], marker="+", label = 'Near-native', zorder = 3)
        # target_df[target_df['Color'] == 'red'].plot(kind = 'scatter', x="ZES3DC", y="GDT_TS", c = 'red', ax = axs[0], marker="+", label = 'Wrong', zorder = 3)
        # target_df[target_df['Color'] == 'blue'].plot(kind = 'scatter', x="ZPAIR", y="GDT_TS", c = 'blue', ax = axs[1], marker="+", label = 'Near-native', zorder = 3)
        # target_df[target_df['Color'] == 'red'].plot(kind = 'scatter', x="ZPAIR", y="GDT_TS", c = 'red', ax = axs[1], marker="+", label = 'Wrong', zorder = 3)
        # handles, labels = axs[0].get_legend_handles_labels() # How to remove the default title: https://stackoverflow.com/questions/51579215/remove-seaborn-lineplot-legend-title?rq=1
        # axs[0].legend(handles=plt.plot([],ls="-", color='black') + handles[1:], labels=labels, loc="best")
        # handles, labels = axs[1].get_legend_handles_labels() # How to remove the default title: https://stackoverflow.com/questions/51579215/remove-seaborn-lineplot-legend-title?rq=1
        # axs[1].legend(handles=plt.plot([],ls="-", color='black') + handles[1:], labels=labels, loc="best")
        # #plt.gray()
        # #plt.grid(True, color='w', linestyle='-', linewidth=1, zorder=0)
        # #plt.gca().patch.set_facecolor('#EAEAF1')
        # pylab.savefig(output_plot, format='png')
        # plt.clf()

    # Write results 
    results_per_target_df.to_csv(results_per_target_file, sep='\t', index=False)

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
