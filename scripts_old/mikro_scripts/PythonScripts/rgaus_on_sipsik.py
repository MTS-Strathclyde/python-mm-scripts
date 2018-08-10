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
import subprocess

#CALC_FOLDER = "/home/mishin1991/calculations/"

EXEC_SCRIPT = """#!/bin/bash
#PBS -l nodes=1:ppn=4
#PBS -q default
cd {0}
/usr/local/bin/rungauss.1 {1}
sshpass -p "Max7120" scp {1}.log sipsik.chem.ut.ee:{2}
sshpass -p "Max7120" scp {1}.chk sipsik.chem.ut.ee:{2}
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
    parser.add_argument('-c', '--chk', help="""Check point file for this calculation. If left empty,
                        script will control, wheather file with the same name but with .chk extension
                        exists in the same folder. Should be used, if chkpoint file
                        has different extension.""")
    parser.add_argument('-f', '--folder', help="""Absolute path of folder in
                        /home/*, where log and chk files will be saved.
                        (/home/mishin1991/calculations/)""",
                        default="/home/mishin1991/calculations/")
    return parser.parse_args(argv)


def main(argv):
    args = process_command_line(argv)
    cwd = os.getcwd()
    for molec in args.input:
        name = fileutil.get_filename(molec)
        if molec[-4:] == ".com":
            shutil.copy(molec, args.folder)
            #Check point file copying
            if args.chk:
                shutil.copy(args.chk, args.folder)
            try:
                shutil.copy(molec[:-4] + ".chk", args.folder)
            except IOError:
                pass
        else:
            print "Molecule " + molec + " doesn't have .com extension and will be ignored."
        exec_script = EXEC_SCRIPT.format(args.folder, name, cwd)
        f = open(args.folder + name, "wb")
        f.write(exec_script)
        f.close()
        subprocess.call(["qsub", f.name])

if __name__ == '__main__':
    main(sys.argv[1:])
