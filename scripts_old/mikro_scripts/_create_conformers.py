#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  6 16:52:42 2012

@author: mishin1991

TODO:
    implement RMSD delete
    decide about output format (right now disabled)

"""

import os
import pybel
import subprocess
import argparse
import shutil
import sys

CONFAB_NAME = 'confabWraper.sh'
DEFAULT_INPUT = 'sdf'
DEFAULT_OUTPUT = 'sdf'

def process_command_line(argv):
    parser = argparse.ArgumentParser(description='Generate conformers')
    #Positional args   
    parser.add_argument('input', help='Input molecule', metavar='<inputmol>')    
    parser.add_argument('output', help='Output molecule', metavar='<outputmol>')
    #Optional arguments
    parser.add_argument('-i', '--informat', help='Format of the input.',
                        metavar='babel_format')
    parser.add_argument('-o', '--outformat', help='Format of the output.',
                        metavar='babel_format')
    return parser.parse_args(argv)

def convert_input_molecule(argv):
    """Converts input into pybel molecule and returns it"""
    if argv.informat:
        return pybel.readfile(argv.informat, argv.input).next()
    else:
        return pybel.readfile(DEFAULT_INPUT, argv.input).next()

def using_confab(pymol):
    """Produce conformers using confab and return output sdf file name"""
    tmp_in_name = pymol.title + '_confab_in_temp.sdf'
    tmp_out_name = pymol.title + '_confab_out_temp.sdf'
    pymol.write('sdf', tmp_in_name)
    #call confab
    subprocess.call([CONFAB_NAME, tmp_in_name, tmp_out_name])
    os.unlink(tmp_in_name)
    return tmp_out_name
    
def using_baloon(pymol):
    """Produce conformers using balloon and return output sdf file name"""
    nGenerations = str(300)
    nconfs = str(200)
    tmp_in_name = pymol.title + '_balloon_in_temp.sdf'
    tmp_out_name = pymol.title + '_balloon_out_temp.sdf'
    pymol.write('sdf', tmp_in_name)
    #call balloon
    subprocess.call(['balloon', '--nGenerations', nGenerations,\
                     '--nconfs', nconfs,\
                     '--addConformerNumberToName',\
                     tmp_in_name, tmp_out_name])
    os.unlink(tmp_in_name)
    return tmp_out_name
    
def concatenate(filename_tuple, outputmol):
    """concateanate all files in tuple and delete them in process
    return concatenated file filename"""
    conc_filename = outputmol
    conc_out_file = open(conc_filename, 'wb')
    for filename in filename_tuple:
        shutil.copyfileobj(open(filename,'rb'), conc_out_file)
        os.unlink(filename)
    return conc_filename
    
def main(argv):
    argv = process_command_line(argv)
    in_pymol = convert_input_molecule(argv)
    out_filename = concatenate((using_baloon(in_pymol), using_confab(in_pymol)),
    argv.output)
    
if __name__ == '__main__':
    main(sys.argv[1:])
