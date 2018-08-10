#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Mon May 13 15:43:03 2013

@author: a92549
"""


import argparse
import shutil
import sys

STRING_BEGIN = """-----------------------------------------------------------------------"""
STRING_END = """======================================================================="""

ATOMS = [5, 6, 24, 25, 28, 29]



def process_command_line(argv):
    """
    Processes arguments and returns namespace of them
    """
    parser = argparse.ArgumentParser(description="""NBO job.""")
    #Positional args
    parser.add_argument('molecules', metavar='mol.chk',
                        help="""Computation to be resumed.""",
                        nargs='+')
    return parser.parse_args(argv)    

    
def main(argv):
    args = process_command_line(argv)
    print 'Name' + str(ATOMS)
    for mol in args.molecules:
        with open(mol, 'rb') as f:
            full_txt = f.read()
        start = full_txt.split(STRING_BEGIN)[1]
        exact = start.split(STRING_END)[0]
        charges_lines = exact.splitlines()
        file_string = mol
        for atom_num in ATOMS:
            #atom_num = atom_num - 1
            charge_line_list = charges_lines[atom_num].split()
            file_string += ';' + charge_line_list[2]
        print file_string
    
    
if __name__ == '__main__':
    main(sys.argv[1:])