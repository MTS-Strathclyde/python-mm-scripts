#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  1 16:30:50 2013

@author: a92549
"""

import os
import sys
import shutil
import argparse

DEFAULT_INPUT = """%nproc=4
%chk={1}
# freq=noraman mpwpw91 IOp(3/76=0572004280) int=ultrafine scrf=(smd,solvent={0}) 
guess=tcheck chkbas genchk geom=allcheck

{2}


"""


def process_command_line(argv):
    """
    Processes arguments and returns namespace of them
    """
    parser = argparse.ArgumentParser(description="""Writes MPW1K solvent
                                    calculation frequency jobs.""")
    #Positional args
    parser.add_argument('input', metavar='mol.chk', help='Input files', nargs='+')
    #Optional args
    parser.add_argument('-d', '--descrip', help='description',
                        default='Frequency calculation.')
    parser.add_argument('-s', '--solvent', help='Solvent name recognized by gaussian (benzene)',
                    metavar='<solvent>', default='benzene')
    return parser.parse_args(argv)


def get_filename(path):
    """Retuns filename, remove path and extension."""
    with_ext = os.path.split(path)[1]
    return os.path.splitext(with_ext)[0]


def main(argv):
    args = process_command_line(argv)
    for chk in args.input:
        name = get_filename(chk) + '_freq'
        formated_input = DEFAULT_INPUT.format(args.solvent, name, args.descrip)
        shutil.copy(chk, name + '.chk')
        with open(name + '.com', 'wb') as f:
            f.write(formated_input)

if __name__ == '__main__':
    main(sys.argv[1:])
