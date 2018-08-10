#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  9 16:29:08 2012

@author: mishin1991

Checks, wheather given mopout molecules have ring.
Deletes, if they do.
"""

import sys
import os
import pybel

def main(argv):
    deleted_i = 0
    for mol in argv:
        pyMol = pybel.readfile("mopout", mol).next()
        mol_sssr = pyMol.sssr.iterator()
        if list(mol_sssr):
            print "list done"
            os.unlink(mol)
            print "Molecule " + mol + " had ring!"
            deleted_i += 1
        else:
            print mol + " is totally ring-free!"
    print "Total deleted " + str(deleted_i)
    
if __name__ == '__main__':
    main(sys.argv[1:])