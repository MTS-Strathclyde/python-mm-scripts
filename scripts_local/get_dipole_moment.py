#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 21 08:22:37 2015

@author: max
"""

import sys


def get_charges(prmtop_name):
    chrgs = []
    with open(prmtop_name, 'r') as f:
        for line in f:
            if line.startswith('%FLAG CHARGE'):
                next(f)
                next_line = next(f)
                while next_line.startswith(' '):
                    chrgs.extend(next_line.split())
                    next_line = next(f)
    return [float(chg)/18.2223 for chg in chrgs]
    

def get_coords(pdb_name):
    with open(pdb_name) as f:
        pdb = f.readlines()
    pdb = [l for l in pdb if l[0:4] == "ATOM"]
    coord = [[float(line[30+i*8:38+i*8])for i in range(3)] for line in pdb]
    return coord


def calc_moment(charges, coord):
    """
    Actually do moment calculation on centered coordinates.
    """

    # unit conversion
    elementary_q = 1.602176487e-19 # elementary charge in C
    debye_conv = 3.33564e-30       # 1 Debye in C*m
    k = elementary_q*1e-10/debye_conv  

    num_atoms = len(charges)

    # Calculate the moment by sum(q(i)*r(i)) for x, y, and z
    moment = [sum([coord[j][i]*charges[j] for j in range(num_atoms)])
            for i in range(3)]
    # convert to debyes
    moment = [m*k for m in moment]
    
    # we are only interested in magnitude, so no direction
    return sum([m**2 for m in moment])**(1./2)

def main(argv):
    if '-h' in argv:
        print 'usage: get_dipole_moment.py mol.pdb mol.prmtop'
        return 0
    pdb_name = argv[0]
    prmtop_name = argv[1]
    print calc_moment(get_charges(prmtop_name), get_coords(pdb_name))
    
if __name__ == '__main__':
    main(sys.argv[1:])
    
    