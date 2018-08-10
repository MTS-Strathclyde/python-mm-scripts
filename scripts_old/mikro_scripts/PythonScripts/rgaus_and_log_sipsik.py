#!/opt/python27/bin/python
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  8 15:36:20 2013

@author: a92549
"""

import argparse
import os
import sys
import fileutil
import shutil
import datetime
import subprocess

START_CALCS_LOG = '/home2/a92549/log/gaus_strt.log'

FIN_LOG_SCRIPT = "/home2/a92549/bin/PythonScripts/log_calc_finish_sipsik.py"

EXEC_SCRIPT = """#!/bin/bash
#PBS -l nodes=1:ppn=4
#PBS -q default
cd {0}
/usr/local/bin/rungauss.1 {1}
sshpass -p "Max7120" scp {1}.log sipsik.chem.ut.ee:{2}
sshpass -p "Max7120" scp {1}.chk sipsik.chem.ut.ee:{2}
sshpass -p "Max7120" ssh sipsik.chem.ut.ee "{3} {2} {1}"
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
    parser.add_argument('-f', '--folder', help="""Absolute path of folder in
                        /home/*, where log and chk files will be saved.
                        (/home/mishin1991/calculations)""",
                        default="/home/mishin1991/calculations")
    return parser.parse_args(argv)


def get_datetime():
    """Return string with current datetime in Y-M-D H:Min format.
    """
    return datetime.datetime.now().strftime("%Y-%m-%d %H:%M")


def move_to_home(molec, args):
    """Move to directory under /home from directory under /home2."""
    shutil.copy(molec, args.folder)
    #Check point file copying
    mol_name = molec.split('.')[0]
    try:
        shutil.copy(mol_name + '.' + args.chk_ext, args.folder)
    except IOError:
        pass

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
        name = fileutil.get_filename(molec)
        if molec[-4:] == ".com":
            move_to_home(molec, args)
        else:
            print "Molecule " + molec + " doesn't have .com extension and will be ignored."
        write_script_and_start_calculation(name, cwd, args)
        with open(START_CALCS_LOG, 'ab') as f:
            f.write(name + ';' + cwd + ';' + get_datetime() + '\n')

if __name__ == '__main__':
    main(sys.argv[1:])
