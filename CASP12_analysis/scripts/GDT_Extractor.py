import os
import pandas as pd
import glob

dataframe_csv = "./Protein_SPS.csv"
dataframe = pd.read_csv(dataframe_csv, sep="\t", index_col = 0)

GDT = list()
TM = list()
QCS = list()


for element in dataframe.PDB:
    txt_file = glob.glob("./TXT/" + element[:5] + "-" + element[5:7] + ".txt" )[0]
    #print (element)
    #print (element[:5] + "TS" + element[9:12] + "_" + element[13:14] + "-" + element[5:7])
    print ("Searching for element", element)

    if len(element) == 11:
        print (element, "Hola")
        GDT.append(1)
        TM.append(1)
        QCS.append(1)
    else:
        with open(txt_file, "r") as txt_input:
            for line in txt_input.readlines()[2:]:
                line = line.split()

                try:
                    if line[1] == (element[:5] + "TS" + element[9:12] + "_" + element[13:14] + "-" + element[5:7]):
                            print ("Found", element)
                            GDT.append(line[3])
                            QCS.append(line[-6])
                            TM.append(line[-4])
                            break

                    else:
                        continue

                except:
                    GDT.append("nan")
                    TM.append("nan")
                    QCS.append("nan")
                    break



print (len(dataframe["PDB"]))
print (len(GDT))
dataframe["GDT"] = GDT
dataframe["QCS"] = QCS
dataframe["TM"] = TM

dataframe.to_csv("finaldata.csv", sep="\t")


