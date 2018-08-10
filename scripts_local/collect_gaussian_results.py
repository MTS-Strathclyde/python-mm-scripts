#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 19 17:26:48 2015

@author: max
"""


from __future__ import print_function

import argparse
import sys
import pandas as pd
import numpy as np


KCAL_IN_HARTREE = 627.509469

def process_command_line(argv):
    """Processes arguments

    Parameters
    ----------
    argv : list
        Command line arguments.

    Returns
    -------
    out : argparse.Namespace
        Namespace of command line arguments.
    """
    parser = argparse.ArgumentParser(description="""Collect results of GAUSSIAN
                        calculation in csv file.""")
    #Positional args
    parser.add_argument('out',
                        help="""GAUSSIAN out files.""", nargs='+')
    #Optional args
    parser.add_argument('-n', '--name',
                        help="""Name of output csv file with the results 
                        [results.csv].""",
                        default='results.csv')
    parser.add_argument('-d', '--decomp',
                        help="""Decompose solvation energy into 
                        electrostatic/cavity contributions""",
                        action='store_true')
    return parser.parse_args(argv)



def parse_gamess_file(fname, decomp=False):
    name = fname[:-4]
    data = [name]
    last_energy = None
    if decomp:
        last_cav  = None
    with open(fname) as f:
        for l in f:
            if l.startswith(' SCF Done:  E'):
                last_energy = float(l.split('=')[1].split()[0])*KCAL_IN_HARTREE
            if decomp:
                if l.startswith(' SMD-CDS'):
                    last_cav = float(l.split()[-1])
    if decomp:
        data.extend([last_energy, last_energy - last_cav, last_cav])
    else:
        data.append(last_energy)
    return data

def main(argv):
    args = process_command_line(argv)
    columns = ['Name', 'Total']
    solv_columns = ['Name', 'Total', 'elec', 'cav']
    data = []
    for f in args.out:
        data.append(parse_gamess_file(f, args.decomp))
    #print(data)
    if args.decomp:
        df = pd.DataFrame(data, columns=solv_columns)
    else:
        df = pd.DataFrame(data, columns=columns)
    df.sort(columns='Name', inplace=True)
    df.to_csv(args.name, index=False)


if __name__ == '__main__':
    main(sys.argv[1:])
    
    