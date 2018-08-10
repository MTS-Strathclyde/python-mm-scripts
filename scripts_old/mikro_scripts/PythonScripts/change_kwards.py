#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 16 20:51:51 2013

@author: a92549
"""

import argparse
import fileutil
import sys
import gaussian
import stringutil
import os
import shutil


METHOD_DIC = {"MPW1K" : "mpwpw91 IOp(3/76=0572004280)"}


def process_command_line(argv):
    """
    Processes arguments and returns namespace of them
    """
    parser = argparse.ArgumentParser(description="""Creates new input files from
                                    old.
                                    If filename doesn't contain TS substring,
                                    it will use chk file to extract geometry and
                                    basis.
                                    If filename contains TS - simply substitutes
                                    one method with another and creates new
                                    com files.""")
    #Positional args
    parser.add_argument('input', metavar='mol.com', help='Input files', nargs='+')
    #Optional args
    parser.add_argument('-m', '--mode', help="""Changes script mode. 1 - means 
                        substitute. This is ideal for changing method or some
                        kwards detail. 2 - means add. Thought for adding
                        additional details to kwards. 3 - is for completely
                        replacing old kwards with new ones (1).""",
                        default=1, type=int)
    parser.add_argument('--add', help="""Will add following kwards (int=ultrafine
                        opt=tight).""", default="int=ultrafine opt=tight")
    parser.add_argument('-o', '--old', help="""Old method. Case-insensitive (B3LYP).""", default="B3LYP")
    parser.add_argument('-n', '--new', help="""New method (MPW1K). This string will
                        be added to name in any mode.""", default="MPW1K")
    parser.add_argument('-d', '--dir', help="""Output dir. (.)""", default='.')
    parser.add_argument('-c', '--chk', help="""Check point file extension. (chk)""",
                        default='chk')
    return parser.parse_args(argv)


def replace_method(args, old_kwards):
    """Returns kwards with changed method."""
    #remove /
    old_kward_ls = old_kwards.split('/')
    old_kwards = ' '.join(old_kward_ls)
    #Not std method
    if args.new in METHOD_DIC:
        method = METHOD_DIC[args.new]
    else:
        method = args.new
    return stringutil.case_insensitive_replace(old_kwards, args.old, method)
    

def replace(args, old_kwards):
    """This function will either change method or replace kwards depending
    on chosen mode."""
    if args.mode == 2:
        return old_kwards + " " + args.add
    elif args.mode == 3:
        return args.new
    elif args.mode == 1:
        return replace_method(args, old_kwards)
    else:
        print "Unknown mode."
        raise Exception


def create_new_TS_input(mol, args, name):
    """Returns new input file with right kwards.""" 
    TS_inp = gaussian.TS_input()
    with open(mol, 'rb') as f:
        TS_inp.provide_existing_input_file(f.read())
    TS_inp.kwards = replace(args, TS_inp.kwards)
    TS_inp.set_chk(fileutil.get_filename(name))
    print TS_inp.kwards
    return TS_inp.write()
    
    
def create_new_normal_input(mol, args, name):
    """Returns new input file."""
    inp = gaussian.Input()
    with open(mol, 'rb') as f:
        inp.provide_existing_input_file(f.read())
    inp.kwards = replace(args, inp.kwards)
    if 'chkbas' not in inp.kwards:
        inp.kwards += " chkbas"
    if "geom=allcheck" not in inp.kwards:
        inp.kwards += " geom=allcheck"
    inp.kwards = stringutil.case_insensitive_replace(inp.kwards, ' gen ', ' genchk ')
    inp.set_chk(fileutil.get_filename(name))
    inp.set_geometry('')
    inp.set_charge('')
    inp.set_mult('')
    inp.custom_basis = ''
    print inp.kwards
    return inp.write()


def main(argv):
    args = process_command_line(argv)
    for mol in args.input:
        name = fileutil.add_to_name_but_keep_ext(mol, '_' + args.new)
        if 'TS' in mol:
            new_inp = create_new_TS_input(mol, args, name)
        else:
            new_inp = create_new_normal_input(mol, args, name)
            old_chk_name = mol.split('.')[0] + '.' + args.chk
            new_chk_name = name.split('.')[0] + '.' + args.chk
            shutil.copy(old_chk_name, os.path.realpath(args.dir) + '/' + new_chk_name)
        with open(os.path.realpath(args.dir) + '/' + name, 'wb') as f:
            f.write(new_inp)


if __name__ == '__main__':
    main(sys.argv[1:])
