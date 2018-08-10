import os
import sys
import csv
import pybel
import openbabel

def naive_RMSD(mol1, mol2):
    mol1.OBMol.ToInertialFrame()
    mol2.OBMol.ToInertialFrame()
    at1 = mol1.atoms
    at2 = mol2.atoms
    sq_dist_t = 0
    tot_at = len(mol1.atoms)
    for at_num in range(tot_at):
        dist = at1[at_num].OBAtom.GetDistance(at2[at_num].OBAtom)
        sq_dist_t += dist**2
    return (sq_dist_t/tot_at)**0.5

def to_PyMol(file_list):
    """Converts given MOPAC output to pybel format. Returns list"""
    return [pybel.readfile("mopout", filename).next() for filename in file_list]

def removeClose(mol_list, minRMSD):
    """Input: list of molecule filenames, minRMSD
        Checks RMSD between every pair in list.
        If their RMSD is lower, then minRMSD - removes 2nd filename from
        list.
        Prints, how much each molecule had similair conformers.
        Returns processed list"""
    #Convert to OBMol-s
    PyMols = to_PyMol(mol_list)
    print "finished converting"
    #Loop
    i = 0
    while i < len(mol_list):
        ref = PyMols[i]  #reference
        j = i + 1
        removed_confs = 0
        while j < len(mol_list):
            targ = PyMols[j] #target
            rmsd = naive_RMSD(ref, targ)
            if rmsd < minRMSD:
                os.unlink(mol_list[j])   #delete file
                mol_list.pop(j)   #remove from both lists
                PyMols.pop(j)
                removed_confs += 1
            else:
                j = j + 1
        #end of inner loop
        print "Molecule",mol_list[i],"had",removed_confs,"similair conformers."
        i = i + 1
    #end of outer loop
    print "finished deleting"


def main():
    min_rmsd = 0.2
    f = open(sys.argv[1], "r")
    reader = csv.reader(f)
    i = 0
    mol_buffer = []
    for line in reader:
        mol_buffer.append(line[0])
        i += 1
        if i % 500 == 0:
            removeClose(mol_buffer, min_rmsd)
            mol_buffer = []
    removeClose(mol_buffer, min_rmsd)

if __name__ == '__main__':
    main()
