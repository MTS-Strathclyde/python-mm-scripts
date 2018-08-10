#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  3 15:53:07 2015

@author: max
"""

from __future__ import print_function

import argparse
import sys
import pandas as pd



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
    parser = argparse.ArgumentParser(description="""Collect results of setschenow
                        calculations in csv file.""")
    #Positional args
    parser.add_argument('out',
                        help="""Setschenow out dirs.""", nargs='+')
    #Optional args
    parser.add_argument('-n', '--name',
                        help="""Name of output csv file with the results 
                        [results.csv].""",
                        default='results.csv')
    return parser.parse_args(argv)



def parse_setsch_dir(dirname):
    dirname = dirname.rstrip('/')
    try:
        with open(dirname + '/setschenow.txt') as f:
            line = f.readline()  # skip empty line
            line = f.readline()
        k_setsch = float(line.split()[2])
    except:
        k_setsch = None
    return dirname, k_setsch

def main(argv):
    args = process_command_line(argv)
    columns = ['Name', 'k_setsch']
    data = []
    for d in args.out:
        data.append(parse_setsch_dir(d))
    #print(data)
    df = pd.DataFrame(data, columns=columns)
    df.sort(columns='Name', inplace=True)
    df.to_csv(args.name, index=False)


if __name__ == '__main__':
    main(sys.argv[1:])



