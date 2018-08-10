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
import datetime
import time
import subprocess

def process_keywords(argv):
    """Process command-line arguments"""
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

def molecule_path(folder_path, mol_name, mol_number):
    """gives unique name and path for molecule"""
    name = mol_name + str(mol_number) + '{0}' + '.dat'
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

def write_molecule(pymol, folder_path, kwrds, mol_name, mol_number):
    """Writes molecule to given directory with given keywords
    and returns it path"""
    mop_string_no_kwrds = pymol.write("mopin").split('\n', 1)[1]
    mopin_string = kwrds + '\n' + mop_string_no_kwrds
    mol_path = molecule_path(folder_path, mol_name, mol_number)
    f = open(mol_path, 'w')
    f.write(mopin_string)
    f.close()
    return mol_path

def calculate_molecule(pymol, folder_path, kwrds, mol_name, mol_number):
    """Returns pybel MOPAC output file or False boolean"""
    mol_path = write_molecule(pymol, folder_path, kwrds, mol_name, mol_number)
    print "Starting " + os.path.basename(mol_path)
    subprocess.call(['mopac', mol_path])
    try:
        out_file = os.path.splitext(mol_path)[0] + '.out'
        mol = pybel.readfile('mopout', out_file).next()
        mol.title = mol_name + str(mol_number)
        return mol
    except IOError:
        print out_file + ' was not found.'
        return False

def calculate_molecules(mol_generator, folder_path, kwrds, pybel_output_file,
                        mol_name):
    """Calculate iteratively all molecules in given multimolecular file"""
    mol_number = 1
    for mol in mol_generator:
        output = calculate_molecule(mol, folder_path, kwrds, mol_name,
                                    mol_number)
        if output:
            pybel_output_file.write(output)
        mol_number += 1
    return True

def create_output_directory_and_file(in_filename, directoryname,
                                     output_format):
    """Create directory, where MOPAC calculation results will be placed
    and multimolecular output file, where results of calculation will be
    writen"""
    #create output directory
    mol_name = os.path.splitext(in_filename)[0]
    if directoryname:
        os.mkdir(directoryname)
        folder_path = os.path.abspath(directoryname)
    else:
        os.mkdir(mol_name + '_MOPAC_DIR')
        folder_path = os.path.abspath(mol_name + '_MOPAC_DIR')
    #create output file
    out_name = mol_name + '_MOPAC_output.' + output_format
    pybel_output_file = pybel.Outputfile(output_format, out_name)
    return folder_path, pybel_output_file, mol_name, out_name

def info(args, out_name, folder_path):
    print "Calculating given molecules in MOPAC"
    print "Input file: " + args.file
    print "Output file: " + out_name
    print "Output directory " + folder_path
    print "Calculation type: " + args.kwards
    print datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + " Started."
    print "Job description: "

def main(argv):
    args = process_keywords(argv)
    mol_generator = pybel.readfile(args.format, args.file)
    folder_path, pybel_output_file, mol_name, out_name\
                        = create_output_directory_and_file(args.file,
                        args.directory, args.output_format)
    info(args, out_name, folder_path)
    start_time = time.time()
    #calculate
    if calculate_molecules(mol_generator, folder_path, args.kwards,
                           pybel_output_file, mol_name):
        print "Successfully finished"
    fin_time = time.time()
    print datetime.datetime.now().strftime("%Y-%m-%d %H:%M") + " Job ended."
    print "Total time in seconds: " + str(round(fin_time - start_time, 3))
    return folder_path

if __name__ == '__main__':
    main(sys.argv[1:])
