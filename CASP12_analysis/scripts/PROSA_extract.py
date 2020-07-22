import pandas as pd
import os, sys

dataframe = pd.read_csv("Proteins.csv", sep="\t")
globalenergy = list()

for filename in dataframe.PDB:

    residuenumber = list()
    residuePROSA = list()

    PROSAfile = os.path.join("CASP12DCall", filename + ".ana.ana")

    with open(PROSAfile, "r") as PROSAinfo:
        lines = PROSAinfo.readlines()[2:]
        for line in lines:
            residuenumber.append(int(line.split()[0]))
            residuePROSA.append(float(line.split()[1]))
                
    dataframe_residual = pd.DataFrame()
    dataframe_residual["ResNum"] = residuenumber
    dataframe_residual["ResPROSA"] = residuePROSA

    globalenergy.append(sum(residuePROSA))
    dataframe_residual.to_csv(os.path.join("CASP12DCall", filename[:-4] + "_PROSARESIDUE.csv"), sep="\t")
            
dataframe["PROSA_PAIR"] = globalenergy
dataframe.to_csv("./Proteins_PROSA.csv", sep="\t")