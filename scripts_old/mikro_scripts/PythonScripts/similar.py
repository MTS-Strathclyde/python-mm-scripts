#!/usr/bin/python

import pybel
import os
import sys
import openbabel as ob
import argparse
import molutil
import operator

""" 
"""

def process_command_line(argv):
    parser = argparse.ArgumentParser(description="""Comapre reference molecule
                                    and target molecules and sort target molecules
                                    by rmsd""")
    #Positional args   
    parser.add_argument('refer', help='Reference molecule') 
    parser.add_argument('compar', help='Compared molecules',
                        nargs=argparse.REMAINDER)    
    #Optional arguments
    parser.add_argument('-f', '--format', help='Format of the reference molecule. (g09)',
                        metavar='babel_format', default='g09')
    parser.add_argument('-cf', '--cmpr_mols_format',
                        help="""Compared molecules format, by default is the
                        same, as input molecules format.""",
                        metavar='babel_format')
    return parser.parse_args(argv)
    
def compare_rmsd_dic(compar_pymols, reference):
    dic = {}
    align = ob.OBAlign()
    align.SetRefMol(reference.OBMol)
    for pymol in compar_pymols:
        align.SetTargetMol(pymol.OBMol)
        if align.Align():
            dic[pymol] = align.GetRMSD()
    return dic

def main(argv):
    args = process_command_line(argv)
    if not args.cmpr_mols_format:
        cmpr_mols_format = args.format
    else:
        cmpr_mols_format = args.cmpr_mols_format
    compar_pymols = molutil.list_to_pymol(cmpr_mols_format, args.compar)
    reference = pybel.readfile(args.format, args.refer).next()
    dic = compare_rmsd_dic(compar_pymols, reference)
    for pymol, similarity in sorted(dic.items(), key=operator.itemgetter(1)):
        print pymol.title, similarity

if __name__ == '__main__':
    main(sys.argv[1:])
