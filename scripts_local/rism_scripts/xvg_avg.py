#!/usr/bin/env python

# -*- coding: utf-8 -*-
"""
Created on Mon Sep  7 13:45:03 2015

@author: max
"""


import argparse
import sys
import os
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
    parser = argparse.ArgumentParser(description=""" Print column averages
                        in xvg files.""")
    #Positional args
    parser.add_argument('xvg',
                        help=""" Xvg files to consider.""", nargs='+')
    #Optional args
    parser.add_argument('-c', '--columns',
                        help="""  Columns to analyze. Partially supports
                        python indexes so arguments like 1:-3 and 1,5,6 will work,
                        but 2:8:2 or 1,
                        Doe[1].""",
                        default='1')
    parser.add_argument( '-s', '--skip',
                        help=""" Skip n first picoseconds [0].""",
                        default=0, type=float)
    parser.add_argument('--nokcal',
                        help=""" Don't convert kJ to kcal""",
                        action='store_true')
    return parser.parse_args(argv)


def parse_cols(matrix, column_idx_string):
    if ',' in column_idx_string:
        cols = column_idx_string.split(',')
        return matrix[:,cols]  #retain dim info
    elif column_idx_string.count(':') == 0:
        return matrix[:,int(column_idx_string),np.newaxis]
    elif column_idx_string.count(':') == 1:
        idxs = map(int, column_idx_string.split(':'))
        return matrix[:,idxs[0]:idxs[1]]
    else:
        raise TypeError('This column specification is unsupported')


def main(argv):
    args = process_command_line(argv)
    for xvg_f in args.xvg:
        with open(xvg_f, 'r') as f:
            lines = f.readlines()     
        for line in lines:
            if 'coul-lambda = ' in line:
                value = line.split('coul-lambda = ')[1]
                coul_lambda = float(value[:-2])
                break
        else:
            coul_lambda = None
        numbers = [l.split() for l in lines if l.strip()[0].isdigit()]
        xvg = np.array(numbers, dtype=float)
        #print xvg.shape
        xvg_e = xvg[xvg[:,0] > args.skip]
        #print xvg_e.shape
        cols = parse_cols(xvg_e, args.columns)
        if coul_lambda:
            print xvg_f, coul_lambda,
        else:
            print xvg_f,
        for col in cols.T:
            rms = np.std(col)/4.184
            print round(np.mean(col)/4.184, 3), round(rms/np.sqrt(len(col)), 3),
        print ''
        

if __name__ == '__main__':
    main(sys.argv[1:])

