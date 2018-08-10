#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Created on Thu Nov 29 15:28:04 2012

@author: a92549

Somewhat updated
"""

import pybel
import argparse
import sys
import fileutil
import molutil
import gaussian

__version__ = 1.1

#Keyword and basis
DEFAULT_KWRDS = """opt=(TS,noeigentest,tight,calcfc) freq {0} int=ultrafine tzvp/tzvpfit"""
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
    parser.add_argument('-f', '--functional', help="""Functional (wb97xd)""",
                        default='wb97xd')
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
    parser.add_argument('-suf', '--suffix', help="""Suffix added to filename (noeig)""",
                        default='_vacuum')
    return parser.parse_args(argv)


def create_kwards(pymol, args):
    kwards = DEFAULT_KWRDS.format(args.functional, args.solvent)
    element_list = molutil.elemental_compos(pymol)
    if gaussian.is_ECP(element_list):
        kwards += ' pseudo=read'
    return kwards
    

def handle_mol(mol_name, args):
    """Accepts path to mol."""
    if args.format != 'gau':
        pymol = pybel.readfile(args.format, mol_name).next()
        mol_txt = pymol.write('gau')
    else:
        with open(mol_name, 'rb') as f:
            mol_txt = f.read()
    mol_name_sans_ext = fileutil.get_filename(mol_name)
    new_name = mol_name_sans_ext + args.suffix
    std_input = gaussian.Input()
    std_input.provide_existing_input_file(mol_txt)
    std_input.set_chk(new_name)
    std_input.set_proc(args.proc)
    std_input.set_memory(args.memory)
    std_input.set_charge(1)
    kwards = create_kwards(pymol, args)
    std_input.set_kwards(kwards)
    #std_input.add_custom_basis_from_dic(BASIS[args.basis])
    with open(new_name + '.com', 'wb') as f:
        std_input.write(f)


    
def main(argv):
    args = process_command_line(argv)
    for mol_name in args.input:
        handle_mol(mol_name, args)
        

if __name__ == '__main__':
    main(sys.argv[1:])
