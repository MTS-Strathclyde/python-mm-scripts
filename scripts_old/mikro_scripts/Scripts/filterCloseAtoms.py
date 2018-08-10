#!/usr/bin/python2.7
import pybel
import sys

#"""(Will check, whether given as an argument molecules have two atoms closer then)
# TO DO!:
# !!!Extend, so it can take many molecules as argument and have adjustable minimal distance!!!


# Checks, wheather distance between any two atoms in molecule is
# bigger then cutoff.
# Shortest distance between two atoms is 0.74 A (H2)
# 0.65 A as critical distance is probably safe enough for any purpose

def isBiggerThenCutoff(pybelMol, cutoff):
    py_atoms = pybelMol.atoms      # Pybel atoms list
    ob_atoms = [atom.OBAtom for atom in py_atoms]   #Conver them to OBAtoms
    total_atoms = len(py_atoms)
    
    #iterate over all pairs
    i = 0
    while i < total_atoms:
        OBatom1 = ob_atoms[i]
        j = i + 1
        #iterate over atoms, which left
        while j < total_atoms:
            dist = OBatom1.GetDistance(ob_atoms[j])
            #check
            if dist < cutoff:
                return False
            j += 1
        #end of inner loop 
        i += 1
    else:
        return True


def test():
    mol = pybel.readstring("smi", "CCCCC")
    mol.make3D()
    assert isBiggerThenCutoff(mol, 1) == True
    assert isBiggerThenCutoff(mol, 3) == False

    mol = pybel.readfile("xyz", "chainSmallerThen1dist.xyz").next()
    assert isBiggerThenCutoff(mol, 1) == False
    assert isBiggerThenCutoff(mol, 0.9) == True

    mol = pybel.readstring("smi", "O")
    assert isBiggerThenCutoff(mol, 0.01) == True
    
    print "Test Finished"
            
if __name__=='__main__':
    #Prints on the screen result, if molecule is given as argument
    # and molecule format is mopout
    file_name = sys.argv[1]
    print file_name
    py_mol = pybel.readfile("mopout", file_name).next()
    print isBiggerThenCutoff(py_mol, 0.65)
    
