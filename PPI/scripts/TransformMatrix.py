import argparse
import pandas as pd
import math, re, sys, os


def main():

    options = parse_user_arguments()
    transform(options)


def parse_user_arguments(*args, **kwds):
    """Parses the arguments of the program"""

    parser = argparse.ArgumentParser(
        description = "Transform a matrix of individual residues in a matrix of pairs of residues",
        epilog      = "@oliva's lab 2017")
    parser.add_argument('-i','--input',dest='input_file',action = 'store',
                        help = """ Input file with the matrix of individual residues """)
    parser.add_argument('-o','--output',dest='output_file',action = 'store',
                        help = """ Name of the output file """)
    options=parser.parse_args()

    return options


def transform(options):
    """Runs the program"""

    inp_file = options.input_file
    inp_file_fd = open(inp_file, 'r')

    out_file = options.output_file

    energies = False
    prot_a_bool = False
    prot_b_bool = False
    end = 0

    df_a = pd.DataFrame()
    df_b = pd.DataFrame()

    exp = "[+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?" # exponential number regex
    regex = '\s*([0-9]+)\s+('+exp+')\s+('+exp+')\s+('+exp+')\s+('+exp+')\s+('+exp+')\s+('+exp+')\s+('+exp+')\s+('+exp+')'
    energies_regex = re.compile(regex)


    ###### PARSING OF THE SCORES FILE ######

    for line in inp_file_fd:

        line = line.strip()

        if line == 'ENERGIES PER RESIDUE':
            energies = True
            continue

        if energies:

            if line == 'PROTEIN A':
                prot_a_bool = True
                prot_b_bool = False
            elif line == 'PROTEIN B':
                prot_a_bool = False
                prot_b_bool = True

            if line == 'END':
                prot_a_bool = False
                prot_b_bool = False
                end += 1

            # After two ends (of protein A and protein B) we end up the parsing
            if end > 1:
                energies = False
                break

            m = energies_regex.search(line)
            if m:
                residue = int(m.group(1))
                energies = [m.group(2),m.group(3),m.group(4),m.group(5),m.group(6),m.group(7),m.group(8),m.group(9)]
                df2 = pd.DataFrame([energies], index=[residue])

                if prot_a_bool:
                    df_a = df_a.append(df2)
                elif prot_b_bool:
                    df_b = df_b.append(df2)

    inp_file_fd.close()


    ###### CALCULATE A Z-SCORE FOR EACH PAIR OF RESIDUES AND POTENTIAL ######

    df_a_norm = pd.DataFrame()
    df_b_norm = pd.DataFrame()


    # Normalize the scores of the protein A
    for residue_a, energies_a in df_a.iterrows():
        norm_energies = []
        for x in xrange(len(energies_a)):
            energy = float(energies_a[x])
            column = df_a.iloc[:,x]
            column = pd.to_numeric(column)
            mean = column.mean()
            std = column.std()
            if std != 0:
                norm = (energy - mean) / std
            else:
                norm = 0
            norm_energies.append(norm)
        df2 = pd.DataFrame([norm_energies], index=[residue_a])
        df_a_norm = df_a_norm.append(df2)

    # Normalize the scores of the protein B
    for residue_b, energies_b in df_b.iterrows():
        norm_energies = []
        for x in xrange(len(energies_b)):
            energy = float(energies_b[x])
            column = df_a.iloc[:,x]
            column = pd.to_numeric(column)
            mean = column.mean()
            std = column.std()
            if std != 0:
                norm = (energy - mean) / std
            else:
                norm = 0
            norm_energies.append(norm)
        df2 = pd.DataFrame([norm_energies], index=[residue_b])
        df_b_norm = df_b_norm.append(df2)



    df_pairs = pd.DataFrame()
    out_file_fd = open(out_file, 'w')

    # Calculate the z-score for each pair and potential
    # Write the output dataframe
    for residue_a, energies_a in df_a_norm.iterrows():

        for residue_b, energies_b in df_b_norm.iterrows():

            pair = '{}-{}'.format(residue_a, residue_b)
            product = [  float(energies_a[x]) * float(energies_b[x])  for x in xrange(len(energies_a)) ]
            df2 = pd.DataFrame([product], index=[pair])
            df_pairs = df_pairs.append(df2)

            out_file_fd.write('{}'.format(pair))
            for result in product:
                out_file_fd.write('\t{:.2e}'.format(result))
            out_file_fd.write('\n')

    out_file_fd.close()

    return

if  __name__ == "__main__":
    main()