#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  4 16:03:34 2013

@author: a92549
"""

import pybel
import argparse
import sys
import gaussian


#Keyword and basis
DEFAULT_KWRDS = """ {0} int=ultrafine gen"""
SOLVENT_MODEL = ' scrf=(smd,solvent={0})'
BASIS = {'svp' : """/storage/a92549/data/Def2-SVP.json""",
         'tzvp' : """/storage/a92549/data/Def2-TZVP.json""",
         'tzvpp' : """/storage/a92549/data/Def2-TZVPP.json""",
         'tzvpd' : """/storage/a92549/data/Def2-TZVPD.json""",
         'tzvppd' : """/storage/a92549/data/Def2-TZVPPD.json"""}

FUNCTIONAL = {'MPWB1K' : "mpwb95 IOp(3/76=0560004400)",
              "B2GP-PLYP" : "b2plyp iop(3/125=0360003600,3/76=0350006500,3/78=0640006400)"}
              


def process_command_line(argv):
    """
    Processes arguments and returns namespace of them
    """
    parser = argparse.ArgumentParser(description="""Converts given molecule
                                    into gaussian file.""")
    #Positional args
    parser.add_argument('input', metavar='mol.log', nargs='+')
    #Optional args
    parser.add_argument('-f', '--functional', help="""Functional (B2GP-PLYP)""",
                        default='B2GP-PLYP')
    parser.add_argument('--format', help="babel format of input file (g09)""",
                        default='g09')
    parser.add_argument('-b', '--basis', help="""Def2 basis set type - svp,
                        tzvp, tzvpp, ... (tzvp)""", default='tzvp')
    parser.add_argument('-p', '--proc', help="""Processors number (8)""",
                        default=8, type=int)
    parser.add_argument('-mem', '--memory', help="""Memory in MB. (7000)""",
                        default=7000, type=int)
    parser.add_argument('-s', '--solvent', help="""Solvent name recognized by 
                        gaussian. For gas phase calculation set to gas_phase 
                        (benzene)""", metavar='<solvent>', default='benzene')
    return parser.parse_args(argv)

    
def main(argv):
    args = process_command_line(argv)
    gau_input = gaussian.Input()
    gau_input.set_proc(args.proc)
    gau_input.set_memory(args.memory)
    if args.functional in FUNCTIONAL:
        functional = FUNCTIONAL[args.functional]
    else:
        functional = args.functional
    kwards = DEFAULT_KWRDS.format(functional) + SOLVENT_MODEL.format(args.solvent)
    gau_input.set_kwards(kwards)
    for mol in args.input:
        name = mol[:-4] + '_' + args.functional + '_' + args.basis
        pymol = pybel.readfile(args.format, mol).next()
        gau_input.provide_pymol(pymol)
        gau_input.set_chk(name)
        gau_input.add_custom_basis_from_dic(BASIS[args.basis])
        with open(name + '.com', 'wb') as f:
            gau_input.write(f)
        
    

if __name__ == '__main__':
    main(sys.argv[1:])    
    