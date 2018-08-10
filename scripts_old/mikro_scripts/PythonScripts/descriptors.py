#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 25 18:53:53 2013

@author: a92549
"""

import sys
import pybel


CHARGES = []
BONDS = [(5,6), (28, 29)]
#CHARGES = [28, 5, 6, 24, 25, 29]
#BONDS = [(28, 5), (5, 6), (6, 24), (24, 25), (25, 29), (29, 28)]


def main(argv):
    print "Name " + str(CHARGES) + str(BONDS)
    for log in argv:
        pymol = pybel.readfile('g09', log).next()
        atoms = pymol.atoms
        charges_string = ''
        for atom_num in CHARGES:
            atom_num = atom_num - 1
            charges_string += str(atoms[atom_num].partialcharge) + ', '
        bond_string = ''
        for bond in BONDS:
            bond = bond[0] - 1, bond[1] - 1
            dist = atoms[bond[0]].OBAtom.GetDistance(atoms[bond[1]].OBAtom)
            bond_string += str(dist) + ', '
        print log + ' ' + charges_string + bond_string
        
        

        
if __name__ == '__main__':
    main(sys.argv[1:])