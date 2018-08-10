import pybel
import os
import sys
import openbabel as ob

"""Checks, weather folder contatins any MOPAC output file molecules, which RMSD
with given molecule is less then input.
To start:
python Similairs.py min_rmsd reference_molecule
"""

def check(mol, rmsd):
    align = ob.OBAlign()
    align.SetRefMol(mol.OBMol)
    for f in os.listdir("."):
        if f[-4:] == ".out":
            align.SetTargetMol(pybel.readfile("mopout", f).next().OBMol)
            if align.Align():
                if align.GetRMSD() < rmsd:
                    print f, align.GetRMSD()

def main():
    rmsd = float(sys.argv[1])
    mol = pybel.readfile("mopout", sys.argv[2]).next()
    check(mol, rmsd)

if __name__ == '__main__':
    main()
