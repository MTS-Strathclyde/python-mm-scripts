#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Created on Fri Feb  8 15:36:20 2013

@author: a92549
"""

import argparse
import os
import sys
import datetime
import subprocess

START_CALCS_LOG = 'log/gaus_strt.log'

FIN_LOG_SCRIPT = "bin/PythonScripts/log_calc_finish_sipsik.py"

EXEC_SCRIPT = """#!/bin/bash
#PBS -l nodes=1:ppn=4
#PBS -q extralong
#PBS -l walltime=1300:00:00
#PBS -l vmem=30gb
#PBS -l mem=30gb
cd {0}
rungauss {1}
{3} {2} {1}"
rm -f {1}
"""


def process_command_line(argv):
    """
    Processes arguments and returns namespace of them
    """
    parser = argparse.ArgumentParser(description="""Script for starting calculations
                                                from sipsik.""")
    #Positional args
    parser.add_argument('input', metavar='<file>.com',
                        help="""Gaussian calculation files.""", nargs='+')
    #Optional args
    parser.add_argument('-c', '--chk_ext', help="""Check point file extension 
                        for this calculation (chk)""", default='chk')
    return parser.parse_args(argv)


def get_datetime():
    """Return string with current datetime in Y-M-D H:Min format.
    """
    return datetime.datetime.now().strftime("%Y-%m-%d %H:%M")


def get_filename(path):
    """Retuns filename, remove path and extension."""
    with_ext = os.path.split(path)[1]
    return os.path.splitext(with_ext)[0]


def write_script_and_start_calculation(name, cwd, args):
    """Write script for PBS system into folder with calculations."""
    exec_script = EXEC_SCRIPT.format(args.folder, name, cwd, FIN_LOG_SCRIPT)
    f = open(args.folder + '/' + name, "wb")
    f.write(exec_script)
    f.close()
    subprocess.call(["qsub", f.name])


def main(argv):
    args = process_command_line(argv)
    cwd = os.getcwd()
    for molec in args.input:
        name = get_filename(molec)
        if molec[-4:] != ".com":
            print "Molecule " + molec + " doesn't have .com extension and will be ignored."
        write_script_and_start_calculation(name, os.getcwd(), args)
        with open(START_CALCS_LOG, 'ab') as f:
            f.write(name + ';' + cwd + ';' + get_datetime() + '\n')

if __name__ == '__main__':
    main(sys.argv[1:])
