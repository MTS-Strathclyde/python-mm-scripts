#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  9 15:46:19 2012
@author: mishin1991

Prints molecules with closest dist between two atoms

"""

import operator
import pybel
import sys

def main(args):
    atom_num1 = 16
    atom_num2 = 19
    dist = {}
    processed = 0
    for mol in args:
        pymol = pybel.readfile("mopout", mol).next()
        py_atoms = pymol.atoms
        obAtom1 = py_atoms[atom_num1 - 1].OBAtom
        obAtom2 = py_atoms[atom_num2 - 1].OBAtom
        dist[mol] = obAtom1.GetDistance(obAtom2)
        processed += 1
        if processed%100 == 0:
            print "Processed {0} molecules".format(processed)
    #sort dictionary
    sorted_tuple = sorted(dist.iteritems(), key=operator.itemgetter(1))
    for mol, dist in sorted_tuple[:15]:
        print "Mol {0} has distance {1}".format(mol, dist)
    
if __name__ == '__main__':
    main(sys.argv[1:])
