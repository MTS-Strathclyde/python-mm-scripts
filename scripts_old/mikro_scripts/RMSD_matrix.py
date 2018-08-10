#!/usr/bin/python
import pybel
import openbabel
import sys
import os
import csv

FORMAT = "g09"


def gen_RMSD_matrix(mol_list):
    """Input: list of molecule filenames
        Output: dictionary, where keys are mol names and their values are lists
        with RMSD values"""
    #Set up OBAling object
    align = openbabel.OBAlign()
    #Convert to OBMol-s
    PyMols = [pybel.readfile(FORMAT, filename).next() for filename in mol_list]
    rows = [[" "] + mol_list]   #all rows
    #loop
    for mol in PyMols:
        row = [mol.OBMol.GetTitle()]
        #setup reference
        align.SetRefMol(mol.OBMol)
        #inner loop
        for mol_in in PyMols:
            #setup target
            align.SetTargetMol(mol_in.OBMol)
            if align.Align():   #if alignment is possible
                rmsd = align.GetRMSD()
                row.append(round(rmsd, 3))    #add to list
        #end of inner loop
        rows.append(row)
    #end of outer loop
    return rows
    
def write_matrix_to_csv(rows):
    """takes as input rows with filenames and RMSD score and
        writes it to the csv file"""
    #created csv filename is same as directory + matrix
    csv_name = os.getcwd().split("/")[-1] + "_matrix.csv"
    #open file
    f = open(csv_name, "w")
    writer = csv.writer(f)
    #write first row, which contains names of molecules
    writer.writerow(rows[0])
    #loop and write all matrix elements
    for i in range(1, len(rows)):
        writer.writerow(rows[i])
    #close file
    f.close()

if __name__ == '__main__':
    args = sys.argv[1:]
    matrix = gen_RMSD_matrix(args)
    #write to file
    write_matrix_to_csv(matrix)
