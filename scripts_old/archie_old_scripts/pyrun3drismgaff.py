#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  1 19:59:53 2013

@author: max

"""

import os
import sys
import subprocess
import argparse
import time
import shutil


SERIAL_SCRIPT = """#$ -V
#$ -cwd
#$ -P fedorov-iles.prj
#$ -q serial-low.q
#$ -j y
#$ -o out.$JOB_ID
#$ -ac runtime="1h"

date
./run3drismgaff.sh {file_name}
sleep 5
{delete_dx}
date
"""


WATER_LIB_PATH = '/users/xpb13212/lib/3drism/{0}'
RUN3DRISM_PATH = '/users/xpb13212/bin/run3drismgaff.sh'
CALCUC_PATH = '/users/xpb13212/bin/calculate_3drismuc.py'
DELETE_DX = 'rm c_*.dx g_*.dx'


#FIN_LOG_SCRIPT = "log_calc_finish.py"


def process_command_line(argv):
    """
    Processes arguments and returns namespace of them
    """
    parser = argparse.ArgumentParser(description="""Script for starting calculations.""")
    #Positional args
    parser.add_argument('input', metavar='<file>.pdb',
                        help="""PDB molecule file.""", nargs='+')
    #Optional args
    parser.add_argument('-w', '--water', metavar='<file>.xvv',
                        help="""Water molecule filename in water lib path (water298.xvv).""",
                        default='water298.xvv')
    parser.add_argument('-k', '--keep_dx', help="""Keep dx files.""",
                        action='store_true')
    return parser.parse_args(argv)


def get_filename(path):
    """Retuns filename, remove path and extension."""
    with_ext = os.path.split(path)[1]
    return os.path.splitext(with_ext)[0]


def write_files_and_run(name, water_name, delet_dx_files):
    """Write script for que system into pwd."""
    os.mkdir(name)
    if delet_dx_files:
        exec_script = SERIAL_SCRIPT.format(file_name=name, delete_dx=DELETE_DX)
    else:
        exec_script = SERIAL_SCRIPT.format(file_name=name, delete_dx='')        
    water_path = WATER_LIB_PATH.format(water_name)
    shutil.copy(water_path, name + '/wat.xvv')
    shutil.copy(CALCUC_PATH, name)
    shutil.copy(RUN3DRISM_PATH, name)
    shutil.copy(name + '.pdb', name)
    if name[0].isdigit():
        qsub_script_name = 'n' + name + '.sh'
    else:
        qsub_script_name = name + '.sh'
    with open(name + '/' + qsub_script_name, "wb") as f:
        f.write(exec_script)
    subprocess.Popen(["qsub", qsub_script_name], cwd=name)
    time.sleep(2)


def main(argv):
    args = process_command_line(argv)
    for molec in args.input:
        name = get_filename(molec)
        if molec[-4:] != ".pdb":
            print "Molecule " + molec + " doesn't have .pdb extension and will be ignored."
        write_files_and_run(name, args.water, not args.keep_dx)


if __name__ == '__main__':
    main(sys.argv[1:])


