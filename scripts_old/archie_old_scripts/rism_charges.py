#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  1 19:59:53 2013

@author: max

"""


import sys
import subprocess
import argparse
import os
import shutil


CHARGE_TYPES = ['b3lyp', 'mp2']


SERIAL_HEADER = """#$ -V
#$ -cwd
#$ -P fedorov-iles.prj
#$ -q serial-low.q
#$ -j y
#$ -o out.$JOB_ID
#$ -ac runtime="1D"

"""

COMPUTE_RISM = """rism_custom_charges.py {name}.pdb {t} am1
"""

#ERROR_CALC_SCRIPT = "~/log/errors.log"



def write_files_and_run(name, exec_script):
    """Write script for que system into pwd."""
    _, name = os.path.split(name)
    script_name = '{}_mp2_and_b3lyp.sh'.format(name)
    if script_name[0].isdigit():
        script_name = 'n' + script_name
    with open(script_name, "wb") as f:
        f.write(exec_script)
    subprocess.call(["qsub", script_name])
#    print 'Calling {} {} {} {}'.format(name, size, first_T, last_T)
    

def main(argv):
    """Assumes that there are mp2 and b3lyp charges gamess output files in the
    format name_(b3lyp/mp2).out in the same directory as pdb files.
    """
    args = process_command_line(argv)
    print args
    for molec in args.input:
        print molec
        if molec.endswith('.pdb'):
            name, ext = os.path.splitext(molec)
            p, no_p_name = os.path.split(name)
            execute_script = SERIAL_HEADER
            if ext == '.pdb':
                for charge in CHARGE_TYPES:
                    calc_dir = no_p_name + '_' + charge
                    os.mkdir(calc_dir)
                    shutil.copy(molec, calc_dir)
                    shutil.copy(name + '_' + charge + '.out', calc_dir)
                    new_name = os.path.join(calc_dir, no_p_name)
                    execute_script += COMPUTE_RISM.format(name=new_name, chg=charge, t=args.temperature)
                write_files_and_run(no_p_name, execute_script)


if __name__ == '__main__':
    main(sys.argv[1:])


