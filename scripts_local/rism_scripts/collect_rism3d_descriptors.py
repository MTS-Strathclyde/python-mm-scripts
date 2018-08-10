#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 22 16:40:40 2015

@author: max
"""


from __future__ import print_function

import argparse
import sys
import os
import pandas as pd
from collections import OrderedDict

##### my modules ######
from rism3d_descriptors import compute_all_shell_descriptors


dic = OrderedDict([('Name' ,               None,),
           ('Xvv' ,                '\tmm_options:  xvvfile=',),
           ('Closure' ,            '\tmm_options:  closure=',),
           ('rism_exchem' ,        'rism_exchem',),
           ('GF' ,                 'rism_exchGF',),
           ('UV' ,                 'rism_potUV',),
           ('PMV' ,                'rism_volume',),
           ('ISc' ,                'dGhyd(ISc)=',),
           ('Time (s)' ,           '3D-RISM runtime: '),
            ])


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
    parser = argparse.ArgumentParser(description="""Collect results of 3D-RISM
                        calculation in csv file.""")
    #Positional args
    parser.add_argument('directories',
                        help=""" Directories containing rism3d results.""",
                        nargs='+')
    #Optional args
    parser.add_argument('-n', '--name',
                        help="""Name of output csv file with the results 
                        [descriptors.csv].""",
                        default='descriptors.csv')
    parser.add_argument('-r', '--shells_vdw_relative',
                        help=""" Distances in angstroms relative to VdW surface
                        of molecule. VdW surface is defined as a surface of a
                        body obtained by overlaping VdW spheres of molecules.""",
                        default=[],type=float, nargs='+')
    return parser.parse_args(argv)



def is_calc_dir(path):
    """Checks whether folder is the source of calculation."""
    has_log = False
    has_prmtop = False
    files = os.listdir(path)
    for f in files:
        f = os.path.join(path, f)
        if f.endswith('.log'):
            has_log = True
        if f.endswith('.prmtop'):
            has_prmtop = True
    return has_log and has_prmtop


def collect_calc_data(path, args):
    """Analyzes calculation directory and return its results."""
    convg = False
    results_dic = [(k , None) for k in dic.iterkeys()]
    results_dic = OrderedDict(results_dic)
    # First figure out the name
    for f in os.listdir(path):
        f = os.path.join(path, f)
        if f.endswith('.prmtop'):
            full_name = f[:-7]
            name = os.path.split(full_name)[1]
            results_dic['Name'] = name
            # add files to namespace
            args.pdb = os.path.join(path, name + '.pdb')
            args.prmtop = os.path.join(path, name + '.prmtop')            
    # Load all files
    try:
        with open(os.path.join(path, name + '.log')) as f:
            log_lines = f.readlines()
        
    except IOError, e:
        if e.errno == 2:
            print(e)
        else:
            raise e
    try:  #simply append result.txt to log lines
        with open(os.path.join(path, 'results.txt')) as f:
            log_lines.extend(f.readlines())
            convg = True
    except IOError, e:
        if e.errno == 2:
            print(e)
        else:
            raise e
    # deal with easy values
    for title, value in dic.iteritems():  # iterate over titles and startwith strings
        if value:
            for l in log_lines:
                if l.startswith(value): #found the string
                    l = l.replace(value, '')
                    strings = l.split()
                    if title == 'Xvv':
                        xvv = os.path.join(path, strings[0])
                        strings[0] = os.path.split(strings[0])[1]
                    results_dic[title] = strings[0].strip()
    if convg:
        all_d_names, all_d = \
        compute_all_shell_descriptors(path, name, args.shells_vdw_relative, xvv)
        for d_name, d_values in zip(all_d_names, all_d):
            results_dic[d_name] = d_values
    return results_dic


def main(argv):
    args = process_command_line(argv)
    data = []
    if not args.shells_vdw_relative:
        raise ValueError('No distances provided')
    for path in args.directories:
        print('Analyzing files in: {}'.format(path))
        if is_calc_dir(path):
            data.append(collect_calc_data(path, args))
        else:
            print('{} is not a calc dir'.format(path))
    columns = max([i.keys() for i in data], key=len)
    df = pd.DataFrame(data, columns=columns)
    df.sort(columns='Name', inplace=True)
    df.to_csv(args.name, index=False)


if __name__ == '__main__':
    main(sys.argv[1:])
