#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 17:32:53 2013

@author: max
View molecule"""

import pybel
import sys
import argparse


def process_command_line(argv):
    """
    Processes arguments and returns namespace of them
    """
    parser = argparse.ArgumentParser(description="""View mol.""")
    #Positional args
    parser.add_argument('mol', metavar='mol',
                        help="""Molecules.""", nargs='+')
    #Optional args
    parser.add_argument('-f', '--format', metavar='babel format',
                        help="""Input molecules format as required
                        by openbabel (pdb).""", default='pdb')
    return parser.parse_args(argv)



def main(argv):
    args = process_command_line(argv)
    for mol in args.mol:
        mol = pybel.readfile(args.format, mol).next()
        mol.draw()

if __name__ == '__main__':
    main(sys.argv[1:])
