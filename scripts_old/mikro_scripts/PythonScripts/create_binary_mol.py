#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 17:37:30 2013

@author: a92549
"""

import sys
import pybel
import argparse
import fileutil
import subprocess
import os
import gaussian


DEFAULT_KWRDS = """opt=tight freq=noraman b3lyp int=ultrafine gen scrf=(smd,solvent={0})"""
MOPAC_KWARDS = "PM7"
MOPAC_LOCATION = '/opt/mopac/MOPAC2012.exe'
DEFAULT_BASIS = """/storage/a92549/data/Def2-SVP.json"""

BASE_SET_LOCATION = """/storage/a92549/baka/minimum_set/base_files"""

GZ_MAT_DIC = {0 : "M\n",
              1 : "M\nX 1 1.5",
              2 : "M\nX 1 1.5\nX 1 1.5 2 180.",
              3 : "M\nX 1 1.5\nX 1 1.5 2 120.\nX 1 1.5 2 120. 3 180.",
              4 : """M\nX 1 1.5\nX 1 1.5 2 109.\nX 1 1.5 2 109. 3 120.
              X 1 1.5 2 109. 3 240.""",
              5 : """M\nX 1 1.5\nX 1 1.5 2 90.\nX 1 1.5 2 90. 3 180.
              X 1 1.5 2 120. 3 -90.\nX 1 1.5 5 120. 3 -90.""",
              6 : """M\nX 1 1.5\nX 1 1.5 2 90.\nX 1 1.5 2 90. 3 -90.
              X 1 1.5 3 90. 2 180.\nX 1 1.5 2 90. 3 180.
              X 1 1.5 2 90. 3 90."""}

def process_command_line(argv):
    """
    Processes arguments and returns namespace of them
    """
    parser = argparse.ArgumentParser(description="""Makes binary molecule of
                                    the form A B x 3D using MOPAC.""")
    #Positional args
    parser.add_argument('bin_mol', metavar='MXn',
                        help="""Molecules formula.
                        Elements and number should be separated by space. For
                        example: Be F 2 .""", nargs=3)
    #Optional args
    parser.add_argument('-p', '--proc', help="""Processors number (8)""",
                        default=8, type=int)
    parser.add_argument('-c', '--charge', help="""Charge (0).""",
                        default=0, type=int)
    parser.add_argument('-m', '--mult', help="""Multiplicity (1).""",
                        default=1, type=int)
    parser.add_argument('-s', '--solvent', help='Solvent name recognized by gaussian (benzene)',
                        metavar='<solvent>', default='benzene')                                                
    return parser.parse_args(argv)


def mopopt(mop_string, name):
    """Optimize using mopac."""
    with open(name + '.dat', 'wb') as f:
        f.write(mop_string)    
    subprocess.call([MOPAC_LOCATION, name + '.dat'])
    if fileutil.file_exist(name + '.arc'):
        opt_pymol = pybel.readfile('mopout', name + '.out').next()
        os.unlink(name + '.dat')
        os.unlink(name + '.arc')
        os.unlink(name + '.out')
        return opt_pymol
    else:
        raise ValueError("Molecule coulnd'b optimized wiht MOPAC")


def create_mopin(bin_mol, mult, charge):
    """Creates MOPAC input"""
    try:
        mol_gzmat = '#a\n\nd\n\n0 1\n' + GZ_MAT_DIC[int(bin_mol[2])] + '\n'
    except IndexError:
        raise ValueError("Implement Z Matrix for this number of atoms")
    mol_gzmat = mol_gzmat.replace('M', bin_mol[0]).replace('X', bin_mol[1])
    pymol = pybel.readstring('gzmat', mol_gzmat)
    mop_string = pymol.write('mopin')
    mopac_kwards = MOPAC_KWARDS
    if mult == 2:
        mopac_kwards += ' DOUBLET '
    elif mult > 2:
        raise ValueError("Make new mopac kwards to support this multiplicity!.")
    if charge > 0:
        mopac_kwards += ' CHARGE=' + str(charge)
    mop_lst = mop_string.splitlines()
    mop_lst[0] = mopac_kwards
    return '\n'.join(mop_lst)


def write_mol(pymol, args, name):
    inp = gaussian.Input()
    inp.provide_pymol(pymol)
    inp.set_proc(args.proc)
    inp.set_memory(4000)
    inp.set_chk(name)
    kwards = DEFAULT_KWRDS.format(args.solvent)
    if gaussian.is_ECP(args.bin_mol[:-1]):
        kwards += ' pseudo=read'
    inp.set_kwards(kwards)
    inp.set_description(name)    
    inp.set_charge(args.charge)
    inp.set_mult(args.mult)
    inp.add_custom_basis_from_dic(DEFAULT_BASIS)
    with open(name + '.com', 'wb') as f:
        inp.write(f)


def main(argv):
    args = process_command_line(argv)
    name = ''.join(args.bin_mol) + '_' + args.solvent
    mopin = create_mopin(args.bin_mol, args.mult, args.charge)
    opt_pymol = mopopt(mopin, name)
    write_mol(opt_pymol, args, name)
    
    
if __name__ == '__main__':
    main(sys.argv[1:])