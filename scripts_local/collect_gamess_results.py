#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 30 17:58:57 2015

@author: max
"""


from __future__ import print_function

import argparse
import sys
import os
import pandas as pd
import numpy as np
import glob
from collections import OrderedDict
from chemistry.amber.readparm import AmberParm
from scipy.spatial import cKDTree




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
    parser = argparse.ArgumentParser(description="""Collect results of GAMESS
                        calculation in csv file.""")
    #Positional args
    parser.add_argument('out',
                        help="""GAMESS out files.""", nargs='+')
    #Optional args
    parser.add_argument('-n', '--name',
                        help="""Name of output csv file with the results 
                        [results.csv].""",
                        default='results.csv')
    return parser.parse_args(argv)



def parse_gamess_file(fname):
    name = fname[:-4]
    data = [name]
    with open(fname) as f:
        for l in f:
            if l.startswith(' DELTA INTERNAL ENERGY') and 'KCAL/MOL' in l:
                data.append(l.split()[-2])
            if l.startswith(' ELECTROSTATIC INTERACTION') and 'KCAL/MOL' in l:
                data.append(l.split()[-2])
            if l.startswith(' CDS INTERACTION') and 'KCAL/MOL' in l:
                data.append(l.split()[-2])
            if l.startswith(' FREE ENERGY OF SOLVATION') and 'KCAL/MOL' in l \
                and 'ATM' not in l:
                data.append(l.split()[-2])
    return data

def main(argv):
    args = process_command_line(argv)
    columns = ['Name', 'DE_int', 'E_elec', 'E_cds', 'DG_solv']
    data = []
    for f in args.out:
        data.append(parse_gamess_file(f))
    #print(data)
    df = pd.DataFrame(data, columns=columns)
    df.sort(columns='Name', inplace=True)
    df.to_csv(args.name, index=False)


if __name__ == '__main__':
    main(sys.argv[1:])
