import os
import sys
import pandas as pd
import glob

native = pd.DataFrame()

PDBs = list()
DOPE = list()

for filename in os.listdir("./CASP12DCall"):
    if filename.endswith(".pdb"):
        try:
            DOPEfile = os.path.join("./CASP12DCall", filename + ".profile")
            with open(DOPEfile, "r") as DOPEenergy:
                PDBs.append(filename)
                line = DOPEenergy.readlines()[3].split()[-1]
                DOPE.append(float(line))
        
        except:
            print ("Error pringao")
            continue

native["PDB"] = PDBs
native["DOPE"] = DOPE
native.to_csv("Proteins.csv", sep="\t")
print(native)