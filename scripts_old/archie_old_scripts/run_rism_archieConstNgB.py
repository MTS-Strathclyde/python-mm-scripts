#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  1 19:59:53 2013

@author: max

T range:
253.15K < T < 383.15K

"""

import sys
import subprocess
import argparse
import csv
import os
import shutil
import time


JOB_LIMIT = 400
JOBS_PER_SCRIPT = 8


SERIAL_HEADER = """#$ -V
#$ -cwd
#$ -P fedorov-iles.prj
#$ -q serial-low.q
#$ -j y
#$ -o {name}.$JOB_ID
#$ -ac runtime="1D"

"""

COMPUTE_RISM = """rism_custom_chargesConstNgB.py {name}.pdb {t} am1
"""

#ERROR_CALC_SCRIPT = "~/log/errors.log"


def process_command_line(argv):
    """
    Processes arguments and returns namespace of them
    """
    parser = argparse.ArgumentParser(description="""Script for starting calculations.""")
    #Positional args
    parser.add_argument('input', metavar='<file>.pdb',
                        help="""PDB molecule file.""", nargs='+')
    #Optional args
    parser.add_argument('-t', '--temperature', 
                        help="""Temperature of calculation. [298.15]""",
                        type=float, defalut=298.15)
#    parser.add_argument('-j', '--jobs_per_script', 
#                        help="""Number of jobs per single script. [10]""",
#                        type=int, defalut=10)                        
#    parser.add_argument('-k', '--keep_dx', help="""Keep dx files.""",
#                        action='store_true')
    return parser.parse_args(argv)


def write_files_and_run(name, size, first_T, last_T, exec_script):
    """Write script for que system into pwd."""
    _, name = os.path.split(name)
    script_name = '{}_{}_{}_{}.sh'.format(name, size, first_T, last_T)    
    if script_name[0].isdigit():
        script_name = 'n' + script_name
    with open(script_name, "wb") as f:
        f.write(exec_script)
    subprocess.call(["qsub", script_name])
#    print 'Calling {} {} {} {}'.format(name, size, first_T, last_T)
    

def run_jobs(job_buffer):
    """Job buffer should have following format:
    [(job_1, t1), (job_2, t2), ...]
    
    Assumes, that we are running jobs from directory, which contains
    needed pdb file.
    """
    script_content = SERIAL_HEADER.format(name=job_buffer[0][0])
    for name, t in job_buffer:
        script_content += COMPUTE_RISM.format(name=name, t=t)
    name = job_buffer[0][0]
    size = len(job_buffer)
    first_T = job_buffer[0][1]
    last_T = job_buffer[-1][1]
    write_files_and_run(name, size, first_T, last_T, script_content)

#  old
#def main(argv):
#    args = process_command_line(argv)
#    for molec in args.input:
#        if molec[-4:] != ".pdb":
#            print "Molecule " + molec + " doesn't have .pdb extension and will be ignored."
#            continue
#        write_files_and_run(molec[-4:], args.t)


def run_csv_f(csv_f):
    """Run rism for molecules and temperatures in csv file.
    
    File should have following structure:
    <name>
    <t1>
    <t2>
    <t3>
    <name2>
    <t4>,
    ....
    
    Name shoudl have 3 digits and a word, for example:
    001tetr
    """
    number_of_jobs_started = 0
    root_dir = os.getcwd()
    job_buffer = []  # contains tuples of format (job_name, temperature)
    with open(csv_f, 'rb') as f:
        rdr = csv.reader(f)
        for row in rdr:
            if '.' not in row[0]:
                name = row[0]
                try:
                    os.mkdir(name)
                except OSError:
                    pass
                shutil.copy(name + '.pdb', name)
                name = os.path.join(name, name)
            else:
                if number_of_jobs_started < JOB_LIMIT:
                    t = str(round(float(row[0]), 2))
                    print 'Starting job {} at t={}'.format(name, t)  
                    if len(job_buffer) < JOBS_PER_SCRIPT:
                        job_buffer.append((name, t))
                    else:
                        run_jobs(job_buffer)
                        job_buffer = [(name, row[0])]
                        number_of_jobs_started += 1
                else:
                    print 'Reached submited job limit!'
                    os.chdir(root_dir)                        
                    with open(csv_f + '.limit', 'wb') as f:
                        wr = csv.writer(f)
                        _, name = os.path.split(name)                        
                        wr.writerow((name,))
                        wr.writerow(row)
                        for row in rdr:
                            wr.writerow(row)
                        break
        else:
            if job_buffer:
                run_jobs(job_buffer)                            


if __name__ == '__main__':
 #   main(sys.argv[1:])
    run_csv_f(sys.argv[1])


