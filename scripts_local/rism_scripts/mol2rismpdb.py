#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 15:58:04 2013

@author: max
"""


import argparse
import sys
import os
import pybel
import subprocess

def process_command_line(argv):
    """
    Processes arguments and returns namespace of them
    """
    parser = argparse.ArgumentParser(description="""Converts molecule
                        to RISM pdb.""")
    #Positional args
    parser.add_argument('files', metavar='file.xyz',
                        help="""Files.""", nargs='+')
    #Optional args
    parser.add_argument('-f', '--format', metavar='babel format',
                        help="""Input molecules format as required
                        by openbabel (xyz).""", default='xyz')
    return parser.parse_args(argv)



def pymol2rism_pdb(mol_txt_lines):
    """Returns string with molecule in rism pdb format"""
    mol_new_lines = []
    atom_counts = {}
    line_format = '{0[0]}{0[1]:>7}{0[2]:>4}{0[3]:>5}{0[4]:>6}{0[5]: 12.3f}{0[6]: 8.3f}{0[7]: 8.3f}{0[8]:>6}{0[9]:>6}{0[10]:>12}'
    for line in mol_txt_lines:
        if line.startswith('HETATM') or line.startswith('ATOM'):
            strings = line.split()
            strings[0] = 'ATOM'
            atom_counts[strings[2]] = atom_counts.get(strings[2], 0) + 1
            strings[2] = strings[2] + str(atom_counts[strings[2]])
            strings[3] = 'MOL'
            strings[5] = float(strings[5])
            strings[6] = float(strings[6])
            strings[7] = float(strings[7])
            new_line = line_format.format(strings)
            mol_new_lines.append(new_line)              
    mol_new_lines.extend(['TER', 'END'])
    return '\n'.join(mol_new_lines)


def min2pdb(mol, mol_format):
    """Convert molecule to pdb"""
    name = os.path.splitext(mol)[0]
    pdb_name = '{}.pdb'.format(name)
    pymol = pybel.readfile(mol_format, mol).next()
    pymol = pybel.readstring('xyz', pymol.write('xyz')) #to standartizise output
#    print pymol.write('pdb')
#    input()
    mol_txt_lines = pymol.write('pdb').splitlines()
    #subprocess.call(['babel', '-i{}'.format(mol_format), mol, '-opdb', pdb_name])
#    with open(pdb_name) as f:        
#        mol_txt_lines = f.readlines()
    rism_pdb = pymol2rism_pdb(mol_txt_lines)
    with open(pdb_name, 'w') as f:
        f.write(rism_pdb)


def main(argv):
    args = process_command_line(argv)
    for mol in args.files:
        min2pdb(mol, args.format)

if __name__ == '__main__':
    main(sys.argv[1:])

        