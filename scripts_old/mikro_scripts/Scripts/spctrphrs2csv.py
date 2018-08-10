#!/usr/bin/python

#Ver: 1.0
#Created: 26.05.2012
#Author: Maksim Mishin

#Calculates spectrophores for given molecules.

# To run:
# ./spectrophores.py molecule(s) [z-matrix babel format]
# in folder, where molecules are located

"""As an argument, will take from input wildcard (IMPORTANT! Wildcard should be
 surounded by quoutes, for example: "*.out" ! ALSO they can only refer to 
molecules in the folder, where script is executed and can not contatin path). All
 filenames, which will satisfy wildcards criteria will be treated as molecules 
subject to fingerprint calculation.

By default, will assume, that molecules of interest are MOPAC output files.
If there will be 2nd argument from command prompt in form of babel recognized
z-Matrix type, will convert files of that kind.

Output will be written to the spec_list.csv file in csv format. If such file
will already exist, output will be written to spec_list1.csv and so on

The output file will have one molecule on each row.
First value will be filename. Others will be calculated spectrophores.
Everything should be separated by commas and different molecules by rowrakes."""

import pybel
import sys
import fnmatch
import os
import csv

def parse_user_input():
    """Depending on the input will return filename wildcard or filename wildcard
    and format specification as well"""
    user_input = sys.argv
    if len(user_input) == 2:
        return user_input[1]
    elif len(user_input) == 3:
        return user_input[1], user_input[2]
    else:
        print "Wrong number of arguments, should be 1 or 2"
        print "Wildcards must be inserted in \" \""
        raise ValueError()

def create_output_file():
    """Returns file object to which output will be written
    At first, it tries to create file named spec_list.csv, but if such file
    exists already, it will try to create spec_list1.csv file, spec_list2.csv
    file and so on..."""
    filename = "spec_list.csv"
    if os.path.isfile(filename): #if such file exists
        i = 0
        while os.path.isfile(filename): #it will moify it's name and check
            i += 1                  #wheather file with modified name exists
            filename = "spec_list" + str(i) + ".csv"
        else:                       #if no - it will return such fiel
            f = open(filename, 'w')
            return f
    else:   #if we are lucky and "spec_list.csv files do not exist
        f = open(filename, 'w')
        return f

def create_molecule(filename, molec_format):
    """creates and returns pybel type molecule object. Molecule
    z-matrix specification should be written as string and be one of the
    babel's supported molecule formats"""
    return pybel.readfile(molec_format, filename).next()

def calculate_spectroph(pybel_molec):
    """as an input takes pybel molecule object. Returns list of spectrophores
        for given molecule"""
    spec = pybel.ob.OBSpectrophore() #Spectrophore class object
    myspec = spec.GetSpectrophore(pybel_molec.OBMol)
    return myspec

def write_spectrophores(files_wild_name, molec_format="mopout"):
    """writes molecules fingerprints to the output file
    As an input will take possibly wildcard filenames and matrices babel
    format"""
    #first create place to write output
    output_f = create_output_file()

    #then create csv wirter object, which will write to that file
    writer = csv.writer(output_f)
    
    #then iterate over molecules and write spectrophores
    for file in os.listdir('.'):
        if fnmatch.fnmatch(file, files_wild_name):
            print file
            pybel_m = create_molecule(file, molec_format)
            spec_list = calculate_spectroph(pybel_m)
            writer.writerow((file,) + spec_list)

def test():
    wild =  parse_user_input()
    print wild
    print fnmatch.fnmatch("ket_7_10_17_14_zeros_n_6.out", wild)

def main():
    usr_input = parse_user_input()
    if isinstance(usr_input, str):
        write_spectrophores(usr_input)
    else:
        write_spectrophores(usr_input[0], usr_input[1])

if __name__ == '__main__':
    main()
