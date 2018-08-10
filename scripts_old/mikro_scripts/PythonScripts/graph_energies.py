#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  6 20:20:41 2012

@author: mishin1991

based on GraphEnergies.py
"""

DEFAULT_INPUT = 'sdf'

import argparse
import pylab
import pybel
import sys
import molutil

def process_command_line(argv):
    """
    Processes arguments and returns namespace of them
    """
    parser = argparse.ArgumentParser(description='Create graph with molecules energies')
    #Positional arguments
    parser.add_argument('files', nargs='+', help='Input molecules')
    #Optional args with value
    parser.add_argument('-f', '--format', help='Format of the input files (g09)',
                        metavar='<babel format>', default='g09')
    return parser.parse_args(argv)


def plot(float_list):
    pylab.plot(float_list, 'bx')
    pylab.title("Energy distribution among files")
    pylab.ylabel("Energy")
    pylab.show()

def get_energies(files, mol_format):
    mol_gen = molutil.list_to_pymol(mol_format, files)
    return sorted([mol.energy for mol in mol_gen])
    
def main(argv):
    args = process_command_line(argv)
    plot(get_energies(args.files, args.format))
    
if __name__ == '__main__':
    main(sys.argv[1:])