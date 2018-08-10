#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Sun Oct  7 19:27:07 2012

@author: mishin1991
"""

import sys
import argparse
import os
import pybel
import subprocess

def process_keywords(argv):
    parser = argparse.ArgumentParser(description="""Calculate all molecules in
    given file using MOPAC""")    
    #Positional args
    parser.add_argument('file', help='Molecular file', metavar='mol.ext')  
    #Optional args
    parser.add_argument('-f', '--format', help='Babel format, default is sdf-',
                        metavar='<format>', default='sdf')
    parser.add_argument('-o', '--output_format', help='Babel format of output\
                        file. Default is sdf.', metavar='<format>', 
                        default='sdf')
    parser.add_argument('-k', '--kwards', help='MOPAC arguments. \
                        Multiple arguments should be surrounded by double\
                        quotes. Default kwards are: "PM6 GNORM=0.01".',
                        metavar='<kwards>', default='PM6 GNORM=0.01')
    parser.add_argument('-d', '--directory', help='Directory, where MOPAC\
                        output files will be saved. Default name is the same,\
                        as input filename.', metavar='<dir>')
    return parser.parse_args(argv)

def molecule_path(pymol, folder_path):
    """gives unique name and path for molecule"""
    name = pymol.title + '{0}' + '.dat'
    path = os.path.join(folder_path, name)
    if os.path.isfile(path.format('')):
        i = 0
        while True:
            if os.path.isfile(path.format(i)):
                i += 1
            else:
                return path.format(i)
    else:
        return path.format('')

def write_molecule(pymol, folder_path, kwrds):
    mop_string_no_kwrds = pymol.write("mopin").split('\n', 1)[1]
    mopin_string = kwrds + '\n' + mop_string_no_kwrds
    mol_path = molecule_path(pymol, folder_path)
    f = open(mol_path, 'w')
    f.write(mopin_string)
    f.close()
    return mol_path
    
def calculate_molecule(pymol, folder_path, kwrds):
    mol_path = write_molecule(pymol, folder_path, kwrds)
    subprocess.call(['mopac', mol_path])
    try:
        out_file = os.path.splitext(mol_path)[0] + 'out'
        return pybel.readfile('mopout', out_file)
    except IOError:
        print out_file + ' was not found.'
        return False
        
def calculate_molecules(mol_generator, folder_path, kwrds, pybel_output_file):
    for mol in mol_generator:
        output = calculate_molecule(mol, folder_path, kwrds)
        if output:
            pybel_output_file.write(output)
    return True
        
def main(argv):
    args = process_keywords(argv)
    mol_generator = pybel.readfile(args.format, args.file)
    #create output directory
    mol_name = os.path.splitext(args.file)[0]
    if args.directory:
        os.mkdir(args.directory)
        folder_path = os.path.abspath(args.directory)
    else:
        os.mkdir(mol_name)
        folder_path = os.path.abspath(mol_name)
    #create output file
    out_name = mol_name + '_MOPAC_output_' + args.output_format
    pybel_output_file = pybel.Outputfile(args.output_format, out_name)
    #calculate        
    if calculate_molecules(mol_generator, folder_path, args.kwards, pybel_output_file):
        print "Successfully finished"
                               
if __name__ == '__main__':
    main(sys.argv[1:])
