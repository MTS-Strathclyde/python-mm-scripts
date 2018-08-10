#!/usr/bin/python
# -*- coding: utf-8 -*-

""" 
Created on Sat Oct  6 18:30:39 2012

@author: mishin1991

Based on delIdentical.py
"""

import sys
import pybel
import os
import openbabel
import argparse

DEFAULT_RMSD = 0.1

def process_command_line(argv):
    parser = argparse.ArgumentParser(description='Delete similair conformers\
                                    in sdf file and produce new one')
    #Positional arguments
    parser.add_argument('input', help='Input sdf molecule', metavar='<insdf>')
    #Optional arguments
    parser.add_argument('-r', '--rmsd', help='RMSD for pruning (defalult is 0.1)',
                        type=float, metavar='<float>')
    return parser.parse_args(argv)

def to_pyMol(sdf_filename):
    """Converts given MOPAC output to pybel format. Returns list"""
    sdf_mol_gen = pybel.readfile('sdf', sdf_filename)
    return list(sdf_mol_gen)
    
def removeClose(sdf_filename, minRMSD):
    """"""
    align = openbabel.OBAlign()
    pyMols = to_pyMol(sdf_filename)
    print "Finished converting"
    #Loop
    i = 0
    total_removed = 0
    while i < len(pyMols):
        referens = pyMols[i].OBMol  #reference
        align.SetRefMol(referens)
        j = i + 1
        removed_confs = 0
        while j < len(pyMols):
            target = pyMols[j].OBMol #target
            align.SetTargetMol(target)
            #Align and ret rmsd
            if align.Align():
                rmsd = align.GetRMSD()
                if rmsd < minRMSD:
                    pyMols.pop(j)
                    removed_confs += 1
                else:
                    j = j + 1
            else:
                print "Couldn't align"
                raise Exception()
        #end of inner loop
        total_removed += removed_confs
        i += 1
    #end of outer loop
    print "finished deleting, total number of removed conformers is",total_removed
    #Create SDF file
    out_name = os.path.splitext(sdf_filename)[0] + '_pruned' + str(minRMSD) +\
                                '.sdf'
    out_mol = pybel.Outputfile('sdf', out_name)
    for pymol in pyMols:    
        out_mol.write(pymol)
    out_mol.close()
    return out_name
        
def main(argv):
    args = process_command_line(argv)
    if args.rmsd:
        out_name = removeClose(args.input, args.rmsd)
    else:
        out_name = removeClose(args.input, DEFAULT_RMSD)
    
if __name__ == '__main__':
    main(sys.argv[1:])