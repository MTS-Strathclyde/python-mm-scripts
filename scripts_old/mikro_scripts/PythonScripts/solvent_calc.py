#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 29 16:10:02 2013

@author: a92549

Updated for new Rg system
"""

DEFAULT_KWRDS = """opt freq=noraman mpwpw91 IOp(3/76=0572004280) scrf=(smd,solvent={0}) guess=tcheck 
chkbas genchk geom=allcheck"""
TS_KWRDS = """mpwpw91 IOp(3/76=0572004280) opt=qst2 Freq=NoRaman Gen scrf=(smd,solvent={0})"""


import argparse
import sys
import shutil
import os
import gaussian
import fileutil

def process_command_line(argv):
    """
    Processes arguments and returns namespace of them
    """
    parser = argparse.ArgumentParser(description="""Prepares solvent calculation.
                                    Input is gas phase calculation check file.""")
    #Positional args
    parser.add_argument('input', metavar='input.com', nargs='+',
                        help="""Gas phase input files.""")
    #Optional args
    parser.add_argument('-s', '--solvent', help='Solvent name recognized by gaussian (benzene)',
                        metavar='<solvent>', default='benzene')
    parser.add_argument('-o', '--output_dir', help="""Directrory, to which
                        produced files will be copied (.)""", default='.')
    parser.add_argument('-k', '--kwards', help="""Gaussian arguments. Should
                        be surounded by double quotes.
                        (opt freq=noraman mpwpw91 IOp(3/76=0572004280) scrf=(smd,solvent={0}) guess=tcheck 
                        chkbas genchk geom=allcheck for non TS
                        and mpwpw91 IOp(3/76=0572004280) opt=qst2 Freq=NoRaman Gen scrf=(smd,solvent={0}
                        if name containsTS)""")
    parser.add_argument('-d', '--descrip', help='Job description (solvent calculation)',
                        metavar='<descrip>', default='solvent calculation')
    parser.add_argument('-e', '--ECP', help="""Molecule contains elements,
                        which require ECP. Will automatically add to kwards
                        pseudo=read""", action='store_true')
    parser.add_argument('-c', '--chk', help="""Check point file extension. (chk)""",
                        default='chk')                        
    return parser.parse_args(argv)
            

def handle_TS(args, mol, name, output_dir):
    input_ = gaussian.TS_input()
    with open(mol, 'rb') as f:
        input_txt = f.read()
    input_.provide_existing_input_file(input_txt)
    if args.kwards:
        kwards = args.kwards
    else:
        kwards = TS_KWRDS.format(args.solvent)
    if args.ECP:
        kwards += ' pseudo=read'
    input_.set_kwards(kwards)
    descrip = input_.descript[0] + ' ' + args.descrip.format(args.solvent)
    input_.set_description(descrip, 1)
    input_.set_chk(name)
    with open(output_dir + '/' + name + '.com', 'wb') as f:
        input_.write(f)


def handle_normal(args, mol, name, output_dir):
    input_ = gaussian.Input()
    if args.kwards:
        kwards = args.kwards
    else:
        kwards = DEFAULT_KWRDS.format(args.solvent)
    if args.ECP:
        kwards += ' pseudo=read'
    input_.set_kwards(kwards)
    descrip = input_.descript + ' ' + args.descrip.format(args.solvent)
    input_.set_description(descrip)
    input_.set_chk(name)
    shutil.copy2(mol[:-4] + '.' + args.chk, output_dir + '/' + name + '.' + args.chk)
    with open(output_dir + '/' + name + ".com", 'wb') as f:
        input_.write(f)
            
            
def main(argv):
    args = process_command_line(argv)
    output_dir = os.path.realpath(args.output_dir)
    for mol in args.input:
        name = fileutil.get_filename(mol) + '_' + args.solvent
        if 'TS' in mol:
            handle_TS(args, mol, name, output_dir)
        else:
            handle_normal(args, mol, name, output_dir)
            
            
if __name__ == '__main__':
    main(sys.argv[1:])
