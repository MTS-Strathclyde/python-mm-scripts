#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 11 07:56:38 2013

@author: mishin1991

Optimized for citronellals with custom basis sets, but can be changed
pretty easily


"""

import argparse
import sys
import gaussian
import fileutil

CUSTOM_BASIS_LOCATION = "/storage/a92549/data/Def2-SVP.json"

def process_command_line(argv):
    """
    Processes arguments and returns namespace of them
    """
    parser = argparse.ArgumentParser(description="""Replaces atoms in gaussian
                            input files. Writes modified
                            files to the same folder.
                            Old atom should be present in input name,
                            because new input name is derrived from old,
                            by the same replacement.
                            The same is true for chk and description.
                            Input will be recognized as TS opt, if it
                            has TS in its name.""")
    #Positional args
    parser.add_argument('input', help='Input files', nargs='+')
    #Optional args
    parser.add_argument('-o', '--old', help="""Old atom. Case-sensitive.""")
    parser.add_argument('-n', '--new', help="""New atom. Case-sensitive.""")
    parser.add_argument('-t', '--times', help="""Times substitution happen""")
    return parser.parse_args(argv)

def atom_replace(txt, args):
    if args.times:
        txt = txt.replace(args.old, args.new, int(args.times))
    else:
        txt = txt.replace(args.old, args.new)
    return txt

def handle_TS(mol_lst, args):
    """ Writes replaced TS input files"""
    inp = gaussian.TS_input()
    for mol in mol_lst:
        with open(mol, 'rb') as f:
            old_input = f.read()
        new_mol_name = atom_replace(mol, args)
        inp.provide_existing_input_file(old_input)
        inp.set_chk(fileutil.get_filename(new_mol_name))
        inp.descript[0] = atom_replace(inp.descript[0], args)
        inp.descript[1] = atom_replace(inp.descript[1], args)   
        inp.geom[0] = atom_replace(inp.geom[0], args)
        inp.geom[1] = atom_replace(inp.geom[1], args)
        inp.add_custom_basis_from_dic(CUSTOM_BASIS_LOCATION)
        with open(new_mol_name, 'wb') as f:
            inp.write(f)
                
def handle_normal(mol_lst, args):
    """ Writes replaced input files"""
    inp = gaussian.Input()
    for mol in mol_lst:
        with open(mol, 'rb') as f:
            old_input = f.read()
        new_mol_name = atom_replace(mol, args)            
        inp.provide_existing_input_file(old_input)
        inp.set_chk(fileutil.get_filename(new_mol_name))
        inp.descript = atom_replace(inp.descript, args)
        inp.geom = atom_replace(inp.geom, args)
        inp.add_custom_basis_from_dic(CUSTOM_BASIS_LOCATION)
        with open(new_mol_name, 'wb') as f:
            inp.write(f)    
    
def main(argv):
    args = process_command_line(argv)
    file_list = args.input
    TS_list = []
    normal_list = []
    for mol in file_list:
        if mol.find('TS') > -1:
            TS_list.append(mol)
        else:
            normal_list.append(mol)
    handle_TS(TS_list, args)
    handle_normal(normal_list, args)
    
if __name__ == '__main__':
    main(sys.argv[1:])