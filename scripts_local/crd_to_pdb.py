#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 20 10:33:37 2015

@author: max
"""

import argparse
import sys



def process_command_line(argv):
    """Processes arguments and returns namespace of them."""
    parser = argparse.ArgumentParser(description=""" Run a molecular dynamics
    simulation with gromacs.""")
    #Positional args
    parser.add_argument('crd', help='Amber crd or pdb file', metavar='mol.inpcrd')
    parser.add_argument('top', help='Amber prmtop', metavar='mol.prmtop')
    #Optional args
    parser.add_argument('-n', '--name', help="""Name of molecules [MOL]""",
                        default='MOL')
    return parser.parse_args(argv)
