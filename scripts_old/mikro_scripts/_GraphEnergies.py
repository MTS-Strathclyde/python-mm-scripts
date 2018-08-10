#!/usr/bin/python
import sys
import pylab
import pybel
import csv

"""Reads bunch of mopout files or csv file and plots found energies"""

def return_energy(pybel_molec):
    """(pybel_molecule)"""
    return pybel_molec.energy

def f_energies(f_names, f_format):
    """(file_name_list)
    iterates over files in list and returns list of energies"""
    energs = []
    for f in f_names:
        molec = pybel.readfile(f_format, f).next()
        energs.append(return_energy(molec))
    return energs

def csv_energies(csv_f, column):
    energs = []  #list of energies
    
    f = open(csv_f, "r")
    reader = csv.reader(f)
    for row in reader:
        energs.append(float(row[column]))
    f.close()
    return energs

if __name__ == '__main__':
    f_names = sys.argv[1:]
    ext = f_names[0].split(".")[-1]
    print "File extension is:",ext
    if ext == "out":
        sorted_e = sorted(f_energies(f_names, "mopout"))
    elif ext == "csv":
        sorted_e = csv_energies(f_names[0], 1)
    elif ext == "log":
        sorted_e = sorted(f_energies(f_names, "g09"))
    pylab.plot(sorted_e, 'bx')
    pylab.title("Energy distribution among files")
    pylab.ylabel("Energy")
    pylab.show()
