import os
import sys
import xmltodict
import pandas as pd

SPS_Pair = list()
SPS_es3dc = list()
SPS_Zpair = list()
SPS_zes3dc = list()

for pdb in os.listdir("./SPServer_results"):
    SPSfile = os.path.join("SPServer_results", pdb)
    OUTfile = os.path.join( "CASP12DCall", pdb[:-4] + "_SPSresidues.csv")
    if SPSfile.endswith(".DS_Store.xml"):
        continue

        #print (SPSfile)

    residuesPAIR = list()
    residuesES3DC = list()
    residuesZPAIR = list()
    residuesZES3DC = list()

    with open(SPSfile) as fd:
        doc = xmltodict.parse(fd.read())
        PDB = doc["xml"]["protein"]["id"]

        for residue in doc["xml"]["protein"]["residues"]["residue"]:
                
            residuesPAIR.append(residue["Epair"])
            residuesES3DC.append(residue["Es3dc"])
            residuesZPAIR.append(residue["ZEpair"])
            residuesZES3DC.append(residue["ZEs3dc"])

    dataframe = pd.DataFrame()
    dataframe["PAIR"] = residuesPAIR
    dataframe["ES3DC"] = residuesES3DC
    dataframe["ZPAIR"] = residuesZPAIR
    dataframe["ZES3DC"] = residuesZES3DC

    dataframe.to_csv(OUTfile, sep="\t")