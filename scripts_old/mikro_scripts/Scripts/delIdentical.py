#!/usr/bin/python2.7
import sys
import pybel
import os
import openbabel

""" python delIdentical.py [minRMS] mol1.ext mol2.ext [mol3 ... ]
    will delete compounds, which have RMSD bigger then minRMS 
    Default RMSD is 0.1
Based on copyUnique"""

def to_OBMol(file_list):
    """Converts given MOPAC output to pybel format. Returns list"""
    return [pybel.readfile("mopout", filename).next().OBMol for filename in file_list]

def removeClose(mol_list, minRMSD):
    """Input: list of pybel molecules, minRMSD
        Checks RMSD between every pair in list.
        If their RMSD is lower, then minRMSD - removes 2nd filename from
        list.
        Prints, how much each molecule had similair conformers.
        Returns processed list"""
    #Set up OBAling object
    align = openbabel.OBAlign()
    #Convert to OBMol-s
    OBMols = to_OBMol(mol_list)
    print "Finished converting"
    #Loop
    i = 0
    total_removed = 0
    while i < len(mol_list):
        referens = OBMols[i]  #reference
        align.SetRefMol(referens)
        j = i + 1
        removed_confs = 0
        while j < len(mol_list):
            target = OBMols[j] #target
            align.SetTargetMol(target)
            #Align and ret rmsd
            if align.Align():
                rmsd = align.GetRMSD()
                if rmsd < minRMSD:
                    os.unlink(mol_list[j])   #delete file
                    mol_list.pop(j)   #remove from both lists
                    OBMols.pop(j)
                    removed_confs += 1
                else:
                    j = j + 1
            else:
                print "Couldn't align"
                raise Exception()
        #end of inner loop
        print "Molecule",mol_list[i],"had",removed_confs,"similair conformers."
        i = i + 1
        total_removed += removed_confs
    #end of outer loop
    print "finished deleting, total number of removed conformers is",total_removed

def is_number(s):
    """checks, wheter given string is float"""
    try:
        float(s)
        return True
    except ValueError:
        return False
        
def main(argv):
    minRMSD = 0.1  #Default RMSD
    #minRMSD update
    if is_number(argv[0]):
        minRMSD = float(argv[0])
        argv = argv[1:] #Remove RMSD from input molecule list
    removeClose(argv, minRMSD)
    
if __name__ == '__main__':
    main(sys.argv[1:])



