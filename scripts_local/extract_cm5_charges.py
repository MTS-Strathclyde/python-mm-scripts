#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  3 12:24:48 2015

@author: max
"""

from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division


import sys
import os
import argparse
import subprocess


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
    parser = argparse.ArgumentParser(description="""Extract cm5 charges.""")
    #Positional args
    parser.add_argument('log_file',
                        help=""" Log file with gaussian calculations.""")
    #Optional args
    parser.add_argument('-o', '--opt',
                        help=""" Output file contains multiple geometries due to 
                        geometry optimization.""",
                        action='store_true')
    return parser.parse_args(argv)


def is_float(number):
    try:
        float(number)
        return True
    except ValueError:
        return False


def cut_file(fname):
    """ Removes all, but last geometry"""
    last_natom = None
    with open(fname) as f:
        for i, l in enumerate(f):
            if l.startswith(' NAtoms= '):
                last_natom = i
        f.seek(0) # rewind file
        txt = f.readlines()
    print(last_natom)
    tail_txt = ''.join(txt[last_natom:])
    tail_fname = fname[:-4] + '_tail.log'
    with open(tail_fname, 'w') as f:
        f.write(tail_txt)
    return tail_fname


def main(argv):
    args = process_command_line(argv)
    if args.opt:
        fname = cut_file(args.log_file)
    else:
        fname = args.log_file
    # call cm5pack
    subprocess.call(['cm5pac', fname])
    chg_f = fname + '.charge'
    # extract charges
    charges = []
    with open(chg_f) as f:
        for l in f:
            strings = l.split()
            if len(strings) == 4 and all(map(is_float, strings)):
                charges.append(float(strings[2]))
    out_f = os.path.splitext(chg_f)[0] + '.chg'
    with open(out_f, 'w') as f:
        while charges:
            for i in range(8):
                try:
                    f.write('{: 10.6f}'.format(charges.pop(0)))
                except IndexError:
                    break
            f.write('\n')
            
if __name__ == '__main__':
    main(sys.argv[1:])

