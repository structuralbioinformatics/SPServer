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
    python /home/quim/PHD/Projects/CASP14/scripts/parse_casp_results.py -i -d -o
    python /Users/quim/Dropbox/UPF/PHD/Projects/SPServer/CASP12_analysis/scripts/parse_casp_results.py -i /Users/quim/Documents/DATA/SPServer/CASP12_analysis/outputs -d /Users/quim/Documents/DATA/SPServer/CASP12_analysis/data -o /Users/quim/Documents/DATA/SPServer/CASP12_analysis/results
    """
    options = parse_options()
    parse_casp_results(options)


def parse_options():
    """
    This function parses the command line arguments and returns an optparse object.
    """

    parser = optparse.OptionParser("parse_casp_results.py  -i INPUT_DIR -d DATA_DIR -o OUTPUT_DIR")

    # Directory arguments
    parser.add_option("-i", action="store", type="string", dest="input_dir", help="Folder with the SPServer results that we want to analyze", metavar="INPUT_DIR")
    parser.add_option("-d", action="store", type="string", dest="data_dir", help="Folder with the SPServer data (tables with GDT_TS, QCS, TM score)", metavar="DATA_DIR")
    parser.add_option("-o", action="store", type="string", dest="output_dir", help="Output folder where to store the results", metavar="OUTPUT_DIR")

    (options, args) = parser.parse_args()

    return options


def parse_casp_results(options):
    """
    From a folder of SPServer results, creates a results file for CASP.
    """
    # Get parameters
    input_dir = options.input_dir
    data_dir = options.data_dir
    output_dir = options.output_dir
    create_directory(output_dir)
    metrics = ['GDT_TS', 'TM score', 'QCS']
    spserver_scoring_functions = ['ES3DC', 'PAIR', 'ZES3DC', 'ZPAIR']
    spserver_zscores = ['ZES3DC', 'ZPAIR']
    alternative_scoring_functions = ['PROSA', 'DOPE']
    scoring_function_to_name = {'ES3DC':'Es3dc', 'PAIR':'Epair', 'ZES3DC':'ZEs3dc', 'ZPAIR':'ZEpair', 'DOPE':'DOPE', 'PROSA':'PROSA_PAIR', 'GDT_TS':'GDT', 'QCS':'QCS', 'TM score':'TM'}


    # Get SPServer global scores
    spserver_input_dir = os.path.join(input_dir, 'SPServer_results')
    output_df_file = os.path.join(output_dir, 'CASP12_SPServer_global_results.txt')
    if not fileExist(output_df_file):
        target_to_model = {}
        model_to_scores = {}
        model_to_length = {}
        columns = ['Target', 'CASPName', 'PDB', 'Length'] + spserver_scoring_functions
        spserver_results_df = pd.DataFrame(columns=columns)
        index = 1

        # Gather the result files of SPServer
        spserver_files = [os.path.join(spserver_input_dir, f) for f in os.listdir(spserver_input_dir) if fileExist(os.path.join(spserver_input_dir, f))]
        for spserver_file in spserver_files:
            try:
                xml_tree = ET.parse(spserver_file)
                root = xml_tree.getroot()
                for protein in root:
                    # Get the name
                    protein_id = protein.find('id').text
                    model_id = protein_id.split('.pdb')[0]
                    target_name = model_id.split('_')[0]
                    model_name = '{}---{}'.format(target_name, protein_id)
                    if len(model_id) > 7:
                        casp_name = model_id[:5] + 'TS' + model_id[9:12] + '_' + model_id[-1:] + '-' + model_id[5:7] # T0860D1_s236m5 ==> T0860TS236_5-D1
                    else:
                        casp_name = model_id
                    # Get SPServer score
                    global_energies = protein.find('global_energies')
                    scores = []
                    for scoring_function in spserver_scoring_functions:
                        score = float(global_energies.find(scoring_function_to_name[scoring_function]).text)
                        scores.append(score)
                        residues = protein.find('residues')
                        # Get the protein length
                        residue_elements = residues.findall('residue')
                        protein_length = len(residue_elements)
                        # Save results dicts
                        model_to_scores.setdefault(model_name, {})
                        model_to_scores[model_name][scoring_function] = score
                    # Save results dicts
                    target_to_model.setdefault(target_name, []).append(model_name)
                    model_to_length[model_name] = protein_length
                    # Save results to dataframe
                    results = [target_name, casp_name, protein_id, protein_length] + scores
                    #print(target_name, protein_id, model_name, protein_length, scores)
                    df2 = pd.DataFrame([results], columns=columns, index=[index])
                    spserver_results_df = spserver_results_df.append(df2) # Add the information to the main data frame
                    index+=1
            except:
                print('Incorrect results file: {}'.format(spserver_file))
        # Write results 
        spserver_results_df = spserver_results_df.sort_values(by=['PDB'])
        spserver_results_df.to_csv(output_df_file, sep='\t', index=False)
    else:
        # Load results
        spserver_results_df = pd.read_csv(output_df_file, sep='\t', index_col=None)


    # Parse data of GDT_TS, QCS, TM score, PROSA and DOPE from a file created by Ruben
    ruben_results_file = os.path.join(input_dir, 'finaldata.csv')
    output_data_file = os.path.join(output_dir, 'CASP12_global_results.txt')

    if not fileExist(output_data_file) and fileExist(ruben_results_file):
        ruben_results_df = pd.read_csv(ruben_results_file, sep='\t', index_col=0)
        columns = ['Target', 'CASPName', 'PDB', 'Length'] + spserver_scoring_functions + ['PROSA', 'DOPE', 'GDT_TS', 'TM score', 'QCS']
        results_df = pd.DataFrame(columns=columns)
        for index, row in spserver_results_df.iterrows():
            model = row['PDB']
            target = row['Target']
            model_data = ruben_results_df[ruben_results_df['PDB'] == model]
            if len(model_data) > 0:
                if model.split('.pdb')[0] == target:
                    gdt = 100.0
                    tm = 1.0
                    qcs = 100.0
                else:
                    gdt = model_data['GDT'].astype(float).tolist()[0]
                    tm = model_data['TM'].astype(float).tolist()[0]
                    qcs = model_data['QCS'].astype(float).tolist()[0]
                dope = model_data['DOPE'].astype(float).tolist()[0]
                prosa = model_data['PROSA_PAIR'].astype(float).tolist()[0]
                results = list(row) + [prosa, dope, gdt, tm, qcs]
                df2 = pd.DataFrame([results], columns=columns, index=[index])
                results_df = results_df.append(df2) # Add the information to the main data frame
            else:
                print('Model {} not found'.format(model))

        # Remove columns with missing values
        results_df= results_df.dropna()

        # Classify by GDT
        results_df['TypeStructure'] = 'Wrong'
        results_df.loc[results_df['GDT_TS'] < 65, 'TypeStructure'] = 'Wrong'
        results_df.loc[results_df['GDT_TS'] == 100, 'TypeStructure'] = 'Native'
        results_df.loc[((results_df['GDT_TS'] >= 65) & (results_df['GDT_TS'] < 100)), 'TypeStructure'] = 'Near-native'

        # Write results 
        results_df = results_df.sort_values(by=['PDB'])
        results_df.to_csv(output_data_file, sep='\t', index=False)
        results_df = pd.read_csv(output_data_file, sep='\t', index_col=None)
    else:
        # Load results
        results_df = pd.read_csv(output_data_file, sep='\t', index_col=None)

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
