#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
"""
Calculate RMSD of 2 molecules.
For usage run:
    calc_rmsd.py -h
"""

import pybel
import openbabel
import argparse
import sys

DEFAULT = 'xyz'

def process_command_line(argv):
    """Processes arguments and returns namespace of them."""
    parser = argparse.ArgumentParser(description="""Calculates RMSD. If both 
    molecules have the same format, specify it only ones. Default format 
    can be changed in script (line 14).""")    
    #Positional args
    parser.add_argument('file1', help='First molecule', metavar='mol.ext')    
    parser.add_argument('file2', help='Second molecule', metavar='mol.ext')
    #Optional args
    parser.add_argument('-f', '--format', help='Format of both molecules',
                        metavar='babel_format')
    parser.add_argument('-g', '--format2', help='Format of second molecule',
                        metavar='babel_format')
    return parser.parse_args(argv)
    
def calculate_rmsd(file1, format1, file2, format2):
    """Returns float"""
    mol1 = pybel.readfile(format1, file1).next().OBMol
    mol2 = pybel.readfile(format2, file2).next().OBMol
    align = openbabel.OBAlign(mol1, mol2, False, False)
    if align.Align():
        return align.GetRMSD()
    else:
        raise Exception("Couldn't align, probably different molecules")
    
def main(args):
    args = process_command_line(args)
    if args.format:
        format1 = args.format
        format2 = args.format
    elif args.format2:
        format2 = args.format2
    else:
        format1 = DEFAULT
        format2 = DEFAULT
    print calculate_rmsd(args.file1, format1, args.file2, format2)
    
if __name__ == '__main__':
    main(sys.argv[1:])
    