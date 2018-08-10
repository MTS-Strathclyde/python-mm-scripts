#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 27 12:06:03 2012
@author: mishin1991

Works only for gaussian input files

Renum.txt of the form
3
2
1

Will make 1st atom in original file atom number 3 in output.
2nd atom - will remain atom number 2
3rd atom will become atom number 1


In other words, original numbering corresponds to row number in txt file
and new numbering is integer on that row.
"""

import argparse
import sys
import os

def process_command_line(argv):
    """
    Processes arguments and returns namespace of them
    """
    parser = argparse.ArgumentParser(description="""Renumbers given molecule
                                    according to specification file.""")
    #Positional args
    parser.add_argument('input', metavar='molec.com',
                        help="""Input molecule.""")
    parser.add_argument('spec', metavar='renum.txt',
                        help="""Text file with new molecule numbering.
                        Details on formating can be found in script.""")
    #Optional args
#    parser.add_argument('-f', '--format', help='Babel format of input file (gau).',
#                        metavar='<format>', default='gau')
    parser.add_argument('-r', '--remove', help="""Will delete given input file
                        and create renumbered file with the same name, as input
                        (false).""",
                        action='store_true')
    parser.add_argument('-o', '--output', help="""Name of outuput file, where
                        renumbered molecule will be written (input filename + renum)""")
    return parser.parse_args(argv)

def parse_input(args):
    """Returns list with new numbers of old atoms (new number being number of element in list)
    and splitted input gau file with job descript, comment, geom + [basis]"""
    spec_file = open(args.spec)
    renum_list = spec_file.read().split()
    input_file = open(args.input)
    input_lst = input_file.read().split('\n\n')
    return renum_list, input_lst
    
def create_new_numbering(renum_list, initial_geom_with_spin):
    """Renumbers according to specification in renumbering list.
    Returns string with new geometry (with spin and mult)."""
    initial_geom_with_spin_lst = initial_geom_with_spin.split('\n')
    spin_mult = initial_geom_with_spin_lst[0]
    initial_geom = initial_geom_with_spin_lst[1:]
    new_geom = ''
    for num in renum_list:
        num = int(num)
        new_geom += initial_geom[num - 1] + '\n'
    return spin_mult + '\n' + new_geom
    
def handle_output(args, new_geom_str):
    if args.output:
        out_f = open(args.output, 'wb')
    else:
        f_name = os.path.splitext(args.input)[0] + '_renum.com'
        out_f = open(f_name, 'wb')
    out_f.write(new_geom_str)
    if args.remove:
        os.unlink(args.input)
    
def main(argv):
    args = process_command_line(argv)
    renum_list, input_lst = parse_input(args)
    new_geom = create_new_numbering(renum_list, input_lst[2])
    input_lst[2] = new_geom
    new_geom_gau = '\n\n'.join(input_lst) + '\n'
    handle_output(args, new_geom_gau)
    
if __name__ == '__main__':
    main(sys.argv[1:])