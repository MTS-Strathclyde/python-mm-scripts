# -*- coding: utf-8 -*-
"""
Created on Sun Dec  2 12:29:42 2012

@author: mishin1991

Module containig various functions for molecule handling.
"""

__version__ = 1.0

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
               
def elemental_compos(pymol):
    """returns list of element symbols, which given pybel molecule contains"""
    atom_nums = [atom.atomicnum for atom in pymol.atoms]
    sorted_no_dub_atoms = sorted(list(set(atom_nums)))
    symbol_list = [element_dic[num] for num in sorted_no_dub_atoms]
    return symbol_list
    
def list_to_pymol(babel_format, file_list):
    """converts all files in given list into pymols.
    Makes most sense for single molecule files.
    """
    return [pybel.readfile(babel_format, file_).next() for file_ in file_list]
    
def create_pymol_dic(babel_format, file_list):
    """links each filename in given list with respective pymol and returns
    dictionary"""
    pymol_dic = {}
    for f_name in file_list:
        pymol_dic[f_name] = pybel.readfile(babel_format, f_name).next()
    return pymol_dic
    
def convert_mol(in_babel_format, input_file, out_babel_format):
    """Reads input file and converts it to other babel format.
    Returns string."""
    mol = pybel.readfile(in_babel_format, input_file).next()
    return mol.write(out_babel_format)

def test():
    mol = pybel.readstring('smi', 'CC(O)(Br)CC(N)[Co]')
    mol.addh()
    #elemental compos
    if not elemental_compos(mol) == ['H', 'C', 'N', 'O', 'Co', 'Br']:
        raise Exception

if __name__ == '__main__':
    test()

    