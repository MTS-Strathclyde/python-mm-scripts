#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  5 09:05:41 2015

@author: max
"""


from __future__ import print_function, division
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
    parser = argparse.ArgumentParser(description="""Join two or more csvs on 
                        the first column.""")
    #Positional args
    parser.add_argument('csvs',
                        help="""CSV files for joining. Should be
                        more than 2""", nargs='+')
                        
    #Optional args
#    parser.add_argument('-c', '--column',
#                        help="""Number of column to join on [1].""",
#                        default=1, type=int)
    parser.add_argument('-n', '--name',
                        help=""" Combined name [join.csv].""",
                        default='join.csv')
    parser.add_argument('-s', '--suffixes',
                        help=""" Suffixes for all columns. Must be the same
                        number as csv files. [_1, _2, _3, ..].""", nargs='+')
    parser.add_argument('-i', '--ignore',
                        help=""" Ignore duplicate column names in merged
                        csvs. Switches suffixing off.""",
                        action='store_true')
    return parser.parse_args(argv)
    

def main(argv):
    args = process_command_line(argv)
    if len(args.csvs) < 2:
        raise ValueError('Too few csv files.')
    if args.suffixes:
        #print(len(args.suffixes))
        assert len(args.suffixes) == len(args.csvs)
    elif args.ignore:
        args.suffixes = ['' for i in range(len(args.csvs))]
    else:
        args.suffixes = ['_{}'.format(i+1) for i in range(len(args.csvs))]
    csvs = [pd.read_csv(csv) for csv in args.csvs]
    # create first csv
    new_col_names = [name + args.suffixes[0] for name in csvs[0].columns[1:]]
    csvs[0].columns = [csvs[0].columns[0]] + new_col_names
    joined = csvs[0]
    # add all csvs to first one using left join
    for suf, csv in zip(args.suffixes[1:], csvs[1:]):
        new_col_names = [name + suf for name in csv.columns[1:]]
        # first column name is propogated from the first file
        csv.columns = [joined.columns[0]] + new_col_names
        if args.ignore:
            # drop duplicating columns
            for col_name in csv.columns[1:]:
                if col_name in joined.columns:
                    csv.drop(col_name, axis=1, inplace=True)
        joined = joined.merge(csv, how='left', on=joined.columns[0])
    joined.to_csv(args.name, index=False)
    

if __name__ == '__main__':
    main(sys.argv[1:])    
    
