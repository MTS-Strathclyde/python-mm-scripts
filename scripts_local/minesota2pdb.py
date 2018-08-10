#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  5 20:25:50 2013

@author: max
"""

import argparse
import sys
import os
import openbabel as ob
import pybel


element_dic = {1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N', 8: 'O', 9: 'F', 10: 'Ne', 11: 'Na',
               12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P', 16: 'S', 17: 'Cl', 18: 'Ar', 19: 'K', 20: 'Ca', 21: 'Sc',
               22: 'Ti', 23: 'V', 24: 'Cr', 25: 'Mn', 26: 'Fe', 27: 'Co', 28: 'Ni', 29: 'Cu', 30: 'Zn', 31: 'Ga',
               32: 'Ge', 33: 'As', 34: 'Se', 35: 'Br', 36: 'Kr', 37: 'Rb', 38: 'Sr', 39: 'Y', 40: 'Zr', 41: 'Nb',
               42: 'Mo', 43: 'Tc', 44: 'Ru', 45: 'Rh', 46: 'Pd', 47: 'Ag', 48: 'Cd', 49: 'In', 50: 'Sn', 51: 'Sb',
               52: 'Te', 53: 'I', 54: 'Xe', 55: 'Cs', 56: 'Ba', 57: 'La', 58: 'Ce', 59: 'Pr', 60: 'Nd', 
               61: 'Pm',62: 'Sm', 63: 'Eu', 64: 'Gd', 65: 'Tb', 66: 'Dy', 67: 'Ho', 68: 'Er', 69: 'Tm', 70: 'Yb',
               71:'Lu', 72: 'Hf', 73: 'Ta', 74: 'W', 75: 'Re', 76: 'Os', 77: 'Ir', 78: 'Pt', 79: 'Au', 80: 'Hg',
               81: 'Tl', 82: 'Pb', 83: 'Bi', 84: 'Po', 85: 'At', 86: 'Rn', 87: 'Fr', 88: 'Ra', 89: 'Ac',
               90:'Th', 91: 'Pa', 92: 'U', 93: 'Np', 94: 'Pu', 95: 'Am', 96: 'Cm', 97: 'Bk', 98: 'Cf', 99: 'Es',
               100: 'Fm', 101:'Md', 102: 'No', 103: 'Lr', 104: 'Rf', 105: 'Db', 106: 'Sg', 107: 'Bh', 108: 'Hs',
               109: 'Mt',110: 'Ds', 111: 'Rg', 112: 'Cn', 114: 'Uuq', 116: 'Uuh'}


def process_command_line(argv):
    """
    Processes arguments and returns namespace of them
    """
    parser = argparse.ArgumentParser(description="""Convert minesota xyz
                        to usual pdb. Ignores charged molecules.""")
    #Positional args
    parser.add_argument('files', metavar='file.xyz',
                        help="""Minesota xyz files.""", nargs='+')
    return parser.parse_args(argv)


#def makeopenbabel(atomcoords, atomnos, charge=0, mult=1):
#    """Create an Open Babel molecule.
#    atomcoords - list of list of floats
#    [[.1, .1, .1], [1.2, 2.1, 0]]
#    atomnos - list of ints
#    [1, 1]
#    """
#    obmol = ob.OBMol()
#    for i, atomno in enumerate(atomnos):
#        coords = atomcoords[i]
#        obatom = ob.OBAtom()
#        obatom.SetAtomicNum(atomno)
#        obatom.SetVector(*coords)
#        obmol.AddAtom(obatom)
#    obmol.ConnectTheDots()
#    obmol.PerceiveBondOrders()
#    obmol.SetTotalSpinMultiplicity(mult)
#    obmol.SetTotalCharge(charge)
#    return obmol


#def pymol2rism_pdb(pymol):
def pymol2rism_pdb(atoms, atomcoords):
    """Returns string with molecule in rism pdb format"""
    #mol_txt_lines = pymol.write('pdb').splitlines()
    mol_new_lines = []
    atom_counts = {}
    line_format = '{0[0]}{0[1]:>7}{0[2]:>4}{0[3]:>5}{0[4]:>6}{0[5]: 12.3f}{0[6]: 8.3f}{0[7]: 8.3f}{0[8]:6.2f}{0[9]:6.2f}{0[10]:>12}'
    num = 0
    for atom, atomcoord in zip(atoms, atomcoords):
        #if line.startswith('HETATM') or line.startswith('ATOM'):
        #strings = line.split()
        num += 1
        strings = ['ATOM', num]
        atom_counts[atom] = atom_counts.get(atom, 0) + 1
        strings.append(atom + str(atom_counts[atom]))
        strings.append('MOL')
        strings.append(1)
        strings.append(atomcoord[0])
        strings.append(atomcoord[1])
        strings.append(atomcoord[2])
        strings.append(1.0)
        strings.append(0.0)
        strings.append(atom)
        new_line = line_format.format(strings)
        mol_new_lines.append(new_line)              
    mol_new_lines.extend(['TER', 'END'])
    return '\n'.join(mol_new_lines)


def min2pdb(mol):
    """Convert molecule to pdb"""
    with open(mol, 'rb') as f:
        lines = f.readlines()
    name = os.path.splitext(os.path.split(mol)[1])[0]
    #title = lines[0]
    charge, mult = map(int, lines[2].split())
    #if charge != 0:
    #    print mol, 'is charged and will be skipped'
    #else:
    atoms = []
    atomcoords = []
    for line in lines[3:]:
        strings = line.split()
        atoms.append(element_dic[int(strings[0])])
        atomcoords.append(map(float, strings[1:]))
    #obmol = makeopenbabel(atomcoords, atomnos)
    #pymol = pybel.Molecule(obmol)
    #print pymol.write('pdb')
    #pymol.title = title               
    rism_pdb = pymol2rism_pdb(atoms, atomcoords)
    with open(name + '.pdb', 'wb') as f:
        f.write(rism_pdb)


def main(argv):
    args = process_command_line(argv)
    for mol in args.files:
        min2pdb(mol)

if __name__ == '__main__':
    main(sys.argv[1:])

        