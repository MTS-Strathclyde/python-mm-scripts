#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
"""
Superimposes two molecules and displays their RMSD.
Writes their superimposition as molecular file.
For usage run:
    calc_rmsd.py -h
"""

import pybel
import openbabel
import argparse
import sys
import os

DEFAULT_INPUT = 'xyz'
DEFAULT_OUTPUT = 'xyz'

def process_command_line(argv):
    """Processes arguments and returns namespace of them."""
    parser = argparse.ArgumentParser(description="""Superimposes two molecules,
    so that their resulting RMSD is minimal, and saves superimposition in 
    separate file. Default input and output format can be changed in script.
    Also returns RMSD.""")    
    #Positional args
    parser.add_argument('file1', help='First molecule', metavar='mol.ext')    
    parser.add_argument('file2', help='Second molecule', metavar='mol.ext')
    #Optional args
    parser.add_argument('-f', '--format', help='Format of the first molecule. If\
                        second molecules format is not specified, assumes, that\
                        it is the same as firsts.', metavar='babel_format')
    parser.add_argument('-g', '--format2', help='Format of second molecule',
                        metavar='babel_format')
    parser.add_argument('-o', '--oformat', help='Format of output molecule',
                        metavar='babel_format')
    parser.add_argument('-n', '--nooutput', help='Do not write output file',
                        action='store_true')    
    return parser.parse_args(argv)
    
def convert_molecules(args):
    """Converts input molecules to pybel molecules and returns tuple of them"""
    if args.format:
        format1 = args.format
        format2 = args.format
    elif args.format2:
        format2 = args.format2
    else:
        format1 = DEFAULT_INPUT
        format2 = DEFAULT_INPUT
    mol1 = pybel.readfile(format1, args.file1).next()
    mol2 = pybel.readfile(format2, args.file2).next()
    return mol1, mol2

def superimpose(mol1, mol2):
    """Superimposes second molecule and returns rmsd"""
    align = openbabel.OBAlign(mol1.OBMol, mol2.OBMol, False, False)
    if align.Align():
        align.UpdateCoords(mol2.OBMol)
        return align.GetRMSD()
    else:
        raise Exception("Couldn't align, probably different molecules")

def catenate_molecules(mol1, mol2):
    """Returns combined pybel molecule"""
    for pyatom in mol2.atoms:
        mol1.OBMol.AddAtom(pyatom.OBAtom)
    return mol1
    
def name(args):
    """Returns name for output, composed from input molecules names"""
    f1_name = os.path.splitext(args.file1)[0]
    f2_name = os.path.splitext(args.file2)[0]
    return f1_name + '+' + f2_name + '.'    
    
def main(args):
    #processes comand line arguments and superimposes
    args = process_command_line(args)
    mol1, mol2 = convert_molecules(args)
    print superimpose(mol1, mol2)
    #writes output file
    if not args.nooutput:
        out_mol = catenate_molecules(mol1, mol2)
        mol_name = name(args)
        if args.oformat:
            out_mol.write(args.oformat, mol_name + args.oformat)
        else:
            out_mol.write(DEFAULT_OUTPUT, mol_name + DEFAULT_OUTPUT)
    
if __name__ == '__main__':
    main(sys.argv[1:])
    