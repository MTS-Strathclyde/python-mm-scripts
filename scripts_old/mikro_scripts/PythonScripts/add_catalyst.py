#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  1 07:43:43 2013

@author: mishin1991

Script is optimized to add catalysts to citronellal and isopulegol and is best
suited for this purpose.

"""

import pybel
import argparse
import sys
import gaussian
import molutil
import fileutil
import subprocess
import os
import mopac


DEFAULT_KWRDS = """opt=tight freq=noraman b3lyp int=ultrafine gen scrf=(smd,solvent={0})"""
TS_KWRDS = """opt=(qst2, tight) freq=noraman b3lyp int=ultrafine gen scrf=(smd,solvent={0})"""
MOPAC_LOCATION = '/opt/mopac/MOPAC2012.exe'
DEFAULT_BASIS = """/storage/a92549/data/Def2-SVP.json"""
MOLEC_NAMES = ['cit', 'iso_iso', 'iso', 'neo_iso_iso', 'neo_iso']
OH_DOWN = ['cit_for_iso_iso_TS.xyz', 'cit_for_neo_iso_TS.xyz', 'neo_iso_min.xyz',
           'iso_iso_min.xyz']



def process_command_line(argv):
    """
    Processes arguments and returns namespace of them
    """
    parser = argparse.ArgumentParser(description="""Converts given molecule
                                    into gaussian file.""")
    #Positional args
    parser.add_argument('catalyst', metavar='cat.dat',
                        help="""Catalyst in mopin format.""")
    parser.add_argument('molecule', metavar='molec.log',
                        help="""Input molecule in g09 format.""", nargs='+')
    #Optional args
    parser.add_argument('-f', '--format', help="""Babel format of input mols (xyz)""",
                        default="xyz")
#    parser.add_argument('-k', '--kwards', help="""Gaussian arguments. Should
#                        be surounded by double quotes.
#                        (B3LYP Opt=Tight Freq Int=UltraFine Gen for non TS
#                        and B3LYP opt=qst2 Freq=NoRaman Gen if name contains
#                        TS)""")
    parser.add_argument('-b', '--basis', help="""In case of custom basis
                        set, specify location of json formated dictionary,
                        which contains this basis.
                        (/storage/a92549/data/Def2-SVP.json)""",
                        default=DEFAULT_BASIS)
    parser.add_argument('-m', '--mult', help="""Multiplicity.""",
                        default=1, type=int)
    parser.add_argument('-n', '--no_opt', help="""In no_opt mode script will
                        not preform any optimization.""", action='store_true')
    parser.add_argument('--localopt', help="""Will use uff, built in pybel
                        for optimization.""", action='store_true')
    parser.add_argument('--no_TS', help="""In no_TS mode script won't write
                        TS input files.""", action='store_true')
    parser.add_argument('-e', '--ECP', help="""Add to kwards pseudo=read""",
                        action='store_true')
    parser.add_argument('-s', '--solvent', help='Solvent name recognized by gaussian (benzene)',
                        metavar='<solvent>', default='benzene')                        
#    parser.add_argument('--OH_down', help="""Special catalyst for compounds, which have
#                        OH facing downward (neo_iso and iso_iso).""")
    return parser.parse_args(argv)


def make_OH_down_catalyst(catalyst_txt):
    """Reads file with catalyst in MOPAC format writen in such way, that it can
    be attached to molecule by catenation."""
    lines = catalyst_txt.splitlines()
    fst_line_list = lines[0].split()
    fst_line_list[5] = str(-float(fst_line_list[5]))
    fst_line = ' '.join(fst_line_list)
    lines[0] = fst_line
    return '\n'.join(lines)


def add_catalyst(mol, input_mol_format, mult, cat_txt):
    """Adds catalyst to molecule and returns molecule in mopin format.
    Catalyst should be supplied as string (because they can be different)"""
    mop_mol = molutil.convert_mol(input_mol_format, mol, 'mopin')
    if mult == 2:
        non_opt_mop_mol = mopac.remove_optimization_and_set_kwards(mop_mol, 'PM7 DOUBLET')
    elif mult == 1:
        non_opt_mop_mol = mopac.remove_optimization_and_set_kwards(mop_mol)
    else:
        print "Make new mopac kwards to support this multiplicity!."
        raise Exception("Make new mopac kwards to support this multiplicity!.")
    return non_opt_mop_mol + cat_txt
    
    
def localopt(mol_with_cat_mopin, name):
    """Optimize using uff built in pybel."""
    pymol = pybel.readstring('mopin', mol_with_cat_mopin)
    pymol.localopt('uff')
    print "Optimized with UFF: " + name
    pymol.write('xyz', name + '.xyz')
    return pymol


def mopopt(mol_with_cat_mopin, name):
    """Optimize using catalyst using mopac."""
    with open(name + '.dat', 'wb') as f:
        f.write(mol_with_cat_mopin)    
    subprocess.call([MOPAC_LOCATION, name + '.dat'])
    if fileutil.file_exist(name + '.arc'):
        print "Optimized with MOPAC: " + name
        return pybel.readfile('mopout', name + '.out').next()
    else:
        print "MOPAC optimization of " + name + " failed."
        return localopt(mol_with_cat_mopin, name)     

    
def pre_opimize(mol_with_cat_mopin, name, locopt):
    """Optimizes only attached catalyst using PM7."""
    if locopt:
        return localopt(mol_with_cat_mopin, name)
    else:
        return mopopt(mol_with_cat_mopin, name)
    
    
def create_molecule_with_cat(mol, input_mol_format, mult, catalyst, name,
                             locopt, no_opt):
    """Creates optimized pymol.
    Catalyst should be supplied as string (because they can be different)"""
    mopin_with_cat = add_catalyst(mol, input_mol_format, mult, catalyst)
    if no_opt:
        return pybel.readstring('mopin', mopin_with_cat)        
    else:
        return pre_opimize(mopin_with_cat, name, locopt)



def write_non_TS_mol(pymol, args, name):
    inp = gaussian.Input()
    inp.provide_pymol(pymol)
    inp.set_proc(8)
    inp.set_memory(4000)
    inp.set_chk(name)
    kwards = DEFAULT_KWRDS.format(args.solvent)
    if args.ECP:
        kwards += ' pseudo=read'    
    inp.set_kwards(kwards)
    inp.set_description(name)    
    inp.set_charge(0)
    inp.set_mult(args.mult)
    inp.add_custom_basis_from_dic(args.basis)
    with open(name + '.com', 'wb') as f:
        inp.write(f)


def write_TS_mol(pymol1, pymol2, args, name):
    inp = gaussian.TS_input()
    inp.provide_pymol(pymol1, 1)
    inp.provide_pymol(pymol2, 2)
    inp.set_proc(8)
    inp.set_memory(4000)
    inp.set_chk(name)
    kwards = TS_KWRDS.format(args.solvent)
    if args.ECP:
        kwards += ' pseudo=read'    
    inp.set_kwards(kwards)
    inp.set_description(name, 1)    
    inp.set_description(name, 2)    
    inp.set_charge('0', 1)
    inp.set_charge('0', 2)
    inp.set_mult(args.mult, 1)
    inp.set_mult(args.mult, 2)
    inp.add_custom_basis_from_dic(args.basis)
    with open(name + '.com', 'wb') as f:
        inp.write(f)


def write_input_files(min_dic, for_dic, args):
    """Function can work only with input from base_files folder."""
    for mol in MOLEC_NAMES:
        if mol in min_dic:
            write_non_TS_mol(min_dic[mol][0], args, min_dic[mol][1])
            if not args.no_TS:
                if mol in for_dic:
                    TS_name = min_dic[mol][1].replace('min', 'TS')
                    write_TS_mol(for_dic[mol][0], min_dic[mol][0], args, TS_name)
                        
                        
def main(argv):
    args = process_command_line(argv)
    min_dic = {}
    for_dic = {}
    with open(args.catalyst, 'rb') as f:
        catalyst_txt = f.read()
    OH_down_catalyst = make_OH_down_catalyst(catalyst_txt)
    for mol in args.molecule:
        name = fileutil.get_filename(mol) + '+' + fileutil.get_filename(os.path.basename(args.catalyst)) + '_' + args.solvent
        if fileutil.remove_path(mol) in OH_DOWN:
            pymol = create_molecule_with_cat(mol, args.format, args.mult, OH_down_catalyst,
                                                 name, args.localopt, args.no_opt)
        else:
            pymol = create_molecule_with_cat(mol, args.format, args.mult, catalyst_txt,
                                                 name, args.localopt, args.no_opt)
        if 'min' in mol:
            mol_type = os.path.basename(mol).split('_min')[0]
            min_dic[mol_type] = pymol, name
        if 'for' in mol:
            mol_type = os.path.basename(mol).split('for_')[1].split('_TS')[0]
            for_dic[mol_type] = pymol, name
    write_input_files(min_dic, for_dic, args)
        
        
if __name__ == '__main__':
    main(sys.argv[1:])
