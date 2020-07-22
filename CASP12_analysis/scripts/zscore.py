import os
import sys



directories = ["/home/ruben/PROSA2/CASP12DCall"]
for directory in directories:
        for filename in os.listdir(directory):
                if filename.endswith(".pdb"): 
                        filename = os.path.join(directory,filename)
                        Prosa.passToProsa('delete *')
                        Prosa.passToProsa('init zscore')
                        Prosa.passToProsa('read pdb ' + filename + ' ' + filename)
                        Prosa.passToProsa('zscore *')
