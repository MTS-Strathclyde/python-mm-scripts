#!/usr/bin/python2.7
import sys
import os
import pybel
import csv
import operator

"""
To start:
python SaveCSV.py [mol1.out, mol2.out, ...] OR [file.csv]

Assumes, that files are in mopout format
Creates or updates csv file with molecule names and their energy. Looks like:

mol_name1, energy1
mol_name2, energy2
mol_name3, energy3
...
All energies are sorted, so energy1 < energy2 < energy3
"""

def write_dic_to_csv(dic):
    """takes as input dictionary with filenames and molecule energies
        writes it to the csv file"""
    #created csv filename is same as directory
    csv_name = os.getcwd().split("/")[-1] + ".csv"
    f = open(csv_name, "w")
    writer = csv.writer(f)
    for name, energy in sorted(dic.items(), key=operator.itemgetter(1)):
        writer.writerow([name, energy])
    f.close()
    return csv_name
    
def main(argv):
    #get extension
    ext = argv[0].split(".")[-1]
    #If argument is bunch of mopout files and no csv file is supllied
    if ext == "out":
        #Dictionary, where energies will be stored
        ener_dic = {}
        #Loop over all supplied files
        for f in argv:
            molec = pybel.readfile("mopout", f).next()
            ener_dic[f] = molec.energy
        f_name = write_dic_to_csv(ener_dic) #Create CSV
        return f_name
    #Else
    #Update csv file, if it has been supplied as argument
    elif ext == "csv":
        csv_f = argv[0]      #for simplicity
        f = open(csv_f, "r")
        f_temp = open(csv_f +"temp", "w")
        reader = csv.reader(f)
        writer = csv.writer(f_temp)
        #loop over all existing rows
        for row in reader:
            if os.path.exists(row[0]):
                writer.writerow(row)
        #closing and renaming files, which left
        f.close()
        os.unlink(csv_f)
        f_temp.close()
        os.rename(csv_f+"temp", csv_f)
        return csv_f
    #Else
    #If extension is neither out nor csv
    else:
        print "Unknown file extension"
        raise Exception("Unknown filename")

if __name__ == "__main__":
    main(sys.argv[1:])
