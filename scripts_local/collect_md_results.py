#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 29 18:47:17 2015

@author: max
"""


from __future__ import print_function

import argparse
import sys
import os
import pandas as pd
import numpy as np

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
    parser = argparse.ArgumentParser(description="""Collect free energies
                        computed from MD using pymbar.""")
    #Positional args
    parser.add_argument('out',
                        help="""MD calculations base directories containing
                        xvg folder with all xvg files analyzed with
                        alchemical_analysis script.""", nargs='+')
    #Optional args
    parser.add_argument('-n', '--name',
                        help="""Name of output csv file with the results 
                        [results.csv].""",
                        default='results.csv')
    parser.add_argument( '--inverse',
                        help=""" Inverse sign of written results.""",
                        action='store_true')
    parser.add_argument( '--nodecomp',
                        help=""" Just collect total free energy.""",
                        action='store_true')
                        
    return parser.parse_args(argv)



def parse_results_file(dname, nodecomp):
    if nodecomp:
        n_last_lines = 1
    else:
        n_last_lines = 3
    data = [dname.rstrip('/')]
    fname = os.path.join(dname, 'xvg/results.txt')
    try:
        with open(fname) as f:
            lines = f.readlines()
            for l in lines[-n_last_lines:]:
                numbers = l[15:].split('+-')
                data.append(float(numbers[-2].split()[1]))
                data.append(float(numbers[-1]))
    except IOError:
        print('results were not found in {}'.format(dname))
    return data

def main(argv):
    args = process_command_line(argv)
    if args.nodecomp:
        columns = ['Name', 'tot_MBAR', 'tot_uMBAR']
    else:        
        columns = ['Name', 'elec_MBAR', 'elec_uMBAR', 'vdw_MBAR', 'vdw_uMBAR',
                   'tot_MBAR', 'tot_uMBAR']
    data = []
    for f in args.out:
        data.append(parse_results_file(f, args.nodecomp))
    df = pd.DataFrame(data, columns=columns)
    if args.inverse:
        for c in ['elec_MBAR',  'vdw_MBAR', 'tot_MBAR']:
            if c in df.columns:
                df[c] = -df[c]
    df.sort(columns='Name', inplace=True)
    df.to_csv(args.name, index=False)


if __name__ == '__main__':
    main(sys.argv[1:])
