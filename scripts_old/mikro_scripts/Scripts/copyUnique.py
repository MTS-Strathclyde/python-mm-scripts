#!/usr/bin/python2.7
import sys
import pybel
import os
import shutil
import openbabel

""" python copyUnique.py [minRMS] mol1.ext mol2.ext [mol3 ... ]
    will copy compounds, which have RMSD bigger then minRMS to another
    folder
    Default RMSD is 0.3"""

def to_OBMol(file_list):
    """Converts given MOPAC output to pybel format. Returns list"""
    return [pybel.readfile("mopout", filename).next().OBMol for filename in file_list]

def removeClose(mol_list, minRMSD):
    """Input: list of molecule filenames, minRMSD
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
    while i < len(mol_list):
        referens = OBMols[i]  #reference
        align.SetRefMol(referens)
        j = i + 1
        removed_confs = 0
        while j < len(mol_list):
            target = OBMols[j] #target
            align.SetTargetMol(target)
            #Align and ret rmsd
            align.Align()
            rmsd = align.GetRMSD()
            if rmsd < minRMSD:
                mol_list.pop(j)   #remove from both lists
                OBMols.pop(j)
                removed_confs += 1
            else:
                j = j + 1
        #end of inner loop
        print "Molecule",mol_list[i],"had",removed_confs,"similair conformers."
        i = i + 1
    #end of outer loop
    return mol_list

def is_number(s):
    """checks, wheter given string is float"""
    try:
        float(s)
        return True
    except ValueError:
        return False

def copy_to_folder(mol_list, minRMSD):
    """copies given molecule names to folder named _present_folder_name_minRMSD"""
    #Get folder name
    folder_name = os.getcwd().split("/")[-1] + str(minRMSD)
    #Make folder
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)
    else:
        folder_name += "c"
        os.makedirs(folder_name)
    #Copy
    for filename in mol_list:
        shutil.copy(filename, os.getcwd() + "/" + folder_name)    
    
if __name__ == '__main__':
    minRMSD = 0.3  #Default RMSD

    args = sys.argv[1:]
    #minRMSD update
    if is_number(args[0]):
        minRMSD = float(args[0])
        args = args[1:] #Remove RMSD from input molecule list
    unique_mols = removeClose(args, minRMSD)
    #copy everything to folder
    print "Copy start!"
    copy_to_folder(unique_mols, minRMSD)



