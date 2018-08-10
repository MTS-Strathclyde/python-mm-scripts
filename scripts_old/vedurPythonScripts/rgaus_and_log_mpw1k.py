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

START_CALCS_LOG = os.environ['START_LOGS_PATH']
FIN_LOG_SCRIPT = "log_calc_finish.py"
FREQ_CALC_SCRIPT = "auto_freq_calc.py"

EXEC_SCRIPT = """#!/bin/bash
#PBS -l nodes=1:ppn={processors}
#PBS -q {queue_type}
#PBS -l vmem=30gb
#PBS -l mem=30gb
cd {file_directory}
rungauss {file_name}
#{log_calc_script} {file_directory} {file_name}
{freq_calc} {file_name}
rungauss {file_name}_freq
{log_calc_script} {file_directory} {file_name}_freq
rm -f {file_name}
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
    parser.add_argument('-p', '--proc', help="""Number of processors to use (4).""",
                        default='4')
    parser.add_argument('-j', '--length', help="""Queue type (long)""",
                        default='long')
    return parser.parse_args(argv)


def get_datetime():
    """Return string with current datetime in Y-M-D H:Min format.
    """
    return datetime.datetime.now().strftime("%Y-%m-%d %H:%M")


def get_filename(path):
    """Retuns filename, remove path and extension."""
    with_ext = os.path.split(path)[1]
    return os.path.splitext(with_ext)[0]


def write_script_and_start_calculation(name, cwd, proc_num, queue_type):
    """Write script for PBS system into folder with calculations."""
    #exec_script = EXEC_SCRIPT.format(cwd, name, FIN_LOG_SCRIPT, proc_num, queue_type)
    exec_script = EXEC_SCRIPT.format(file_name=name, file_directory=cwd,
                                     log_calc_script=FIN_LOG_SCRIPT,
                                     processors=proc_num, queue_type=queue_type,
                                     freq_calc=FREQ_CALC_SCRIPT)
    f = open(cwd + '/' + name, "wb")
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
        write_script_and_start_calculation(name, cwd, args.proc, args.length)
        with open(START_CALCS_LOG, 'ab') as f:
            f.write(name + '_freq;' + cwd + ';' + get_datetime() + '\n')

if __name__ == '__main__':
    main(sys.argv[1:])
