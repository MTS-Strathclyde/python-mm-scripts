#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  1 12:50:10 2015

@author: max
"""


from __future__ import print_function

import argparse
import sys
import os
import pandas as pd
from subprocess import PIPE, Popen
import glob
import shlex
from collections import OrderedDict


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
    parser = argparse.ArgumentParser(description="""Collect avg quantiities
                         from GROMACS edr files.""")
    #Positional args
    parser.add_argument('out',
                        help="""MD calculations base directories containing
                        edr file in folder called MD.""", nargs='+')
    #Optional args
    parser.add_argument('-a', '--averages',
                        help=""" Input that is going to be passed to all
                        edr files. If there are multiple arguments 
                        they should be surrounded by quotes.""")
    parser.add_argument( '--args',
                        help=""" Additional arguments that are going
                        to be passed to gmx energy. Should be surrounded
                        by quotes.""")
    parser.add_argument('-n', '--name',
                        help="""Name of output csv file with the results 
                        [results.csv].""",
                        default='results.csv')
    parser.add_argument('-d', '--del_xvg',
                        help=""" Delete output energy.xvg file.""",
                        action='store_true')
    parser.add_argument('--nokcal',
                        help=""" Don't convert kJ to kcal""",
                        action='store_true')
    return parser.parse_args(argv)


def format_avgs(avg_lines, nokcal):
    row = OrderedDict()
    for l in avg_lines:
        l = l.split()
        prop_name = l[0]
        units = l[-1]
        if not nokcal and units == '(kJ/mol)':
            conv = 1./4.184
        else:
            conv = 1.
        row[prop_name + '_avg'] = float(l[1])*conv
        row[prop_name + '_err'] = float(l[2])*conv
        row[prop_name + '_rmsd'] = float(l[3])*conv
        row[prop_name + '_drift'] = float(l[4])*conv
    return row

def get_avgs(edr_file, averages, args, del_energy, nokcal):
    """Read volume from edr file"""
    number_of_params = len(averages.split())
    #with open(os.devnull, 'w') as fnull:
    p = Popen(['gmx', 'energy', '-f', edr_file] + shlex.split(args), stdout=PIPE, 
              stdin=PIPE)#, stderr=fnull)
    out = p.communicate(input=averages)[0]
    avg_lines = out.splitlines()[-number_of_params:]
    if del_energy:
        os.unlink('energy.xvg')
    return format_avgs(avg_lines, nokcal)


def get_edr(f):
    """ f is a directory containing edr file.
    returns either filename or None"""
    try:
        edr_name = glob.glob(os.path.join(f, 'MD', '*.edr'))[0]
    except IndexError:
        print('Edr files were not found in {}'.format(f))
        edr_name = None
    return edr_name


def main(argv):
    args = process_command_line(argv)
    data = {}
    for f in args.out:
        f = f.rstrip('/')
        name = os.path.split(f)[0]
        edr_file = get_edr(f)
        if edr_file:
            row = get_avgs(edr_file, args.averages, args.args, args.del_xvg,
                           args.nokcal)
        else:
            row = None
        data[name] = row
    df = pd.DataFrame.from_dict(data, orient='index')
    df.sort(inplace=True)
    df.to_csv(args.name, index_label='Name')


if __name__ == '__main__':
    main(sys.argv[1:])

