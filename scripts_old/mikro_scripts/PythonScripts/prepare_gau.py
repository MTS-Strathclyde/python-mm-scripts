#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Created on Thu Nov 29 15:28:04 2012

@author: a92549

Somewhat updated
"""

import pybel
import argparse
import os
import sys
import fileutil
import gaussian

__version__ = 1.1


DEFAULT_KWRDS = """opt freq m06l/gen/w06"""
TS_KWRDS = """opt=(calcfc,tight,ts,noeigentest) freq wb97xd/gen/w06
scrf=(smd,solvent={0}) int=ultrafine"""
DEFAULT_BASIS = """/storage/a92549/data/Def2-TZVP.json"""


def process_command_line(argv):
    """
    Processes arguments and returns namespace of them
    """
    parser = argparse.ArgumentParser(description="""Converts given molecules
                                    into gaussian input files. If input is gaussian
                                    file with TS in its name - will treat it
                                    as qst2 TS search input file.""")
    #Positional args
    parser.add_argument('input', metavar='molec.ext',
                        help="""Input molecules.""", nargs='+')
    #Optional args
    parser.add_argument('-f', '--format', help='Babel format (gau).',
                        metavar='<format>', default='gau')
    parser.add_argument('-k', '--kwards', help="""Gaussian arguments. Should
                        be surounded by double quotes.
                        (opt=tight freq wb97xd/gen/w06 int=ultrafine scrf=(smd,solvent={0}) or in TS case
		                opt=(calcfc,tight,ts,noeigentest) freq rwb97xd
                       scrf=(smd,solvent={0}) int=ultrafine gen/w06))""")
    parser.add_argument('-b', '--basis', help="""In case of custom basis
                        set, specify location of json formated dictionary,
                        which contains this basis.
                        (/storage/a92549/data/Def2-TZVP.json)""",
                        default=DEFAULT_BASIS)
    parser.add_argument('-e', '--ECP', help="""Add to kwards pseudo=read""",
                        action='store_true')
    parser.add_argument('-s', '--solvent', help='Solvent name recognized by gaussian (benzene)',
                        metavar='<solvent>', default='benzene')      
    parser.add_argument('-d', '--dir', help="""Output dir. (.)""", default='.')
    parser.add_argument('-n', '--name', help="""Name modifier.""", default='')    
    return parser.parse_args(argv)


def handle_mol(mol_name, args):
    """Accepts path to mol."""
    if args.format != 'gau':
        pymol = pybel.readfile(args.format, mol_name).next()
        mol_txt = pymol.write('gau')
    else:
        with open(mol_name, 'rb') as f:
            mol_txt = f.read()
    mol_name_sans_ext = fileutil.get_filename(mol_name)
    new_name = mol_name_sans_ext + args.name
    std_input = gaussian.Input()
    std_input.provide_existing_input_file(mol_txt)
    std_input.set_chk(new_name)
    if args.kwards:
        kwrds = args.kwards
    else:
        kwrds = DEFAULT_KWRDS.format(args.solvent)
    if args.ECP:
        kwrds += ' pseudo=read'
    std_input.set_kwards(kwrds)
    std_input.add_custom_basis_from_dic(args.basis)
    with open(os.path.realpath(args.dir) + '/' + new_name + '.com', 'wb') as f:
        std_input.write(f)


def handle_TS(mol_name, args):
    """Accepts path to gaussian input TS molecule."""
    with open(mol_name, 'rb') as f:
        input = gaussian.Input()
        input.provide_existing_input_file(f.read())
    new_name = fileutil.add_to_name_but_keep_ext(mol_name, args.name)
    input.set_chk(fileutil.get_filename(new_name))
    if args.kwards:
        kwrds = args.kwards
    else:
        kwrds = TS_KWRDS.format(args.solvent)
    if args.ECP:
        kwrds += ' pseudo=read'
    input.set_kwards(kwrds)
    input.add_custom_basis_from_dic(args.basis)
    with open(os.path.realpath(args.dir) + '/' + new_name, 'wb') as f:
        input.write(f)

    
def main(argv):
    args = process_command_line(argv)
    for mol_name in args.input:
        if 'TS' in mol_name:
            handle_TS(mol_name, args)
        else:
            handle_mol(mol_name, args)
        

if __name__ == '__main__':
    main(sys.argv[1:])
