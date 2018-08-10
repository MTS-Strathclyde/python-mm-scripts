#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 16 11:15:08 2013

@author: a92549

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
from threading import Thread


#Keyword and basis
DEFAULT_KWRDS = """opt=tight freq=noraman {functional} int=ultrafine gen"""
TS_KWRDS = """opt=(qst2, tight) freq=noraman {functional} int=ultrafine gen"""
GAS_PHASE_CALC_KWRD = 'gas_phase'
SOLVENT_MODEL = 'scrf=(smd,solvent={solvent})'
HROTOR = 'freq=HinderedRotor'
TEMPR = 'Temperature={tempr}'
BASIS = {'svp' : """/storage/a92549/data/Def2-SVP.json""",
         'tzvp' : """/storage/a92549/data/Def2-TZVP.json""",
         'tzvpp' : """/storage/a92549/data/Def2-TZVPP.json""",
         'tzvpd' : """/storage/a92549/data/Def2-TZVPD.json""",
         'tzvppd' : """/storage/a92549/data/Def2-TZVPPD.json"""}

#Isopulegols and citronellal geometries
BASE_SET_LOCATION = """/storage/a92549/baka/minimum_set/base_files"""
BASE_FORMAT = 'xyz'
MIN_FILES = {'cit' : ['cit_min.xyz', False]
            #'iso' : ['iso_min.xyz', False],
            #'neo_iso' : ['neo_iso_min.xyz', True],
            #'iso_iso' : ['iso_iso_min.xyz', False],
            #'neo_iso_iso' : ['neo_iso_iso_min.xyz', True]
            }
TS_PAIRS = {'iso' : ['iso_min.xyz', 'cit_for_iso_TS.xyz', False],
            'neo_iso' : ['neo_iso_min.xyz', 'cit_for_neo_iso_TS.xyz', True]
            #'iso_iso' : ['iso_iso_min.xyz', 'cit_for_iso_iso_TS.xyz', False],
            #'neo_iso_iso' : ['neo_iso_iso_min.xyz', 'cit_for_neo_iso_iso_TS.xyz', True]
            }

#Catalyst
GENERAL_CATALYST_M_LANE = """{M} {M_dist} 1 {M_C_angle} 1 {M_H_dihed} 1 28 5 14"""
GENERAL_CATALYST_X_LANE = """{X} {X_dist} 1 {X_O_angle} 1 {X_C_dihed} 1 30 28 5"""

#Mopac
MOPAC_KWARDS = "PM7"
MOPAC_LOCATION = '/opt/mopac/MOPAC2012.exe'
MOPOPT_VERT_ANGLE = 90
MOPIN_HORIS_ANGLE = 45


def process_command_line(argv):
    """
    Processes arguments and returns namespace of them
    """
    parser = argparse.ArgumentParser(description="""Converts given molecule
                                    into gaussian file.""")
    #Positional args
    parser.add_argument('catalyst', metavar='MXn',
                        help="""Catalysts general formula.
                        Elements and number should be separated by space. For
                        example: Be F 2 .""", nargs=3)
    #Optional args
    parser.add_argument('-f', '--functional', help="""Functional (B3LYP)""",
                        default='B3LYP')
    parser.add_argument('-b', '--basis', help="""Def2 basis set type - svp,
                        tzvp or tzvpp (svp)""", default='svp')
    parser.add_argument('-p', '--proc', help="""Processors number (8)""",
                        default=8, type=int)
    parser.add_argument('-mem', '--memory', help="""Memory in MB. (6000)""",
                        default=6000, type=int)
    parser.add_argument('-c', '--charge', help="""Charge (0).""",
                        default=0, type=int)
    parser.add_argument('-m', '--mult', help="""Multiplicity (1).""",
                        default=1, type=int)
    parser.add_argument('-s', '--solvent', help="""Solvent name recognized by 
                        gaussian. For gas phase calculation set to gas_phase 
                        (benzene)""", metavar='<solvent>', default='benzene')
    parser.add_argument('-t', '--tempr', help="""Sets temperature in kwards.""")
    parser.add_argument('-hr', '--hrotor', help="""Sets HinderedRotor option in kwards.""",
                        action='store_true')
    parser.add_argument('-md', '--m_dist', help='Distance between metal and oxygen (2.0)',
                        default='2.0', type=float)
    parser.add_argument('-xd', '--x_dist', help='Distance between X and metal (1.9)',
                        default='1.9', type=float)                        
    parser.add_argument('-mc', '--m_c_angle', help='Angle between metal and carbon (120.0)',
                        default='120.0', type=float)
    parser.add_argument('-xo', '--x_o_angle', help='Angle between X and oxygen (120.0)',
                        default='120.0', type=float)                        
    parser.add_argument('-mh', '--m_h_diheds', help='Dihed between metal and hydrogen (45.0 60.0)',
                        default=['45.0', '60.0'], nargs='+')
    parser.add_argument('-xc', '--x_c_dihed', help='Starting angle for x_c dihedrals (180)',
                        default='180', type=float)
    parser.add_argument('-n', '--nomopac', help="Don't create optimized with mopac TS",
                        action='store_true')
    parser.add_argument('-l', '--localopt', help="Optimize with localopt as well",
                        action='store_true')
    return parser.parse_args(argv)
    
    
def get_kwards(job_type, args):
    """Prepares kwards for given job type."""
    if job_type == 'min':
        kwards = DEFAULT_KWRDS
    else:
        kwards = TS_KWRDS
    kwards = kwards.format(functional=args.functional)
    if args.solvent != GAS_PHASE_CALC_KWRD:
        solvent_model = ' ' + SOLVENT_MODEL.format(solvent=args.solvent)
        kwards += solvent_model
    if gaussian.is_ECP(args.catalyst[:-1]):
        kwards += ' pseudo=read'
    if args.tempr:
        tempr = ' ' + TEMPR.format(tempr=args.tempr)
        kwards += tempr
    if args.hrotor:
        kwards += ' ' + HROTOR
    return kwards


def make_folder(args):
    """Creates folder for minimum set and returns its name"""
    folder_name = ''.join(args.catalyst) + '_' + args.solvent + '_' + args.functional + '_' + args.basis
    if args.tempr:
        folder_name += '_' + args.tempr + 'K'
    os.mkdir(folder_name)
    return folder_name
    

def add_catalyst(mol, input_mol_format, mult, charge, cat_txt):
    """Adds catalyst to molecule and returns molecule in mopin format.
    Catalyst should be supplied as string (because they can be different)
    Optimization parameters (1) are put only for catalyst atoms"""
    mop_mol = molutil.convert_mol(input_mol_format, mol, 'mopin')
    mopac_kwards = MOPAC_KWARDS
    if mult == 2:
        mopac_kwards += ' DOUBLET '
    elif mult > 2:
        raise ValueError("Make new mopac kwards to support this multiplicity!.")
    if charge > 0:
        mopac_kwards += ' CHARGE=' + str(charge)
    non_opt_mop_mol = mopac.remove_optimization_and_set_kwards(mop_mol, mopac_kwards)    
    return non_opt_mop_mol + cat_txt
    
    
def localopt(mol_with_cat_mopin, _):
    """Optimize using uff built in pybel."""
    pymol = pybel.readstring('mopin', mol_with_cat_mopin)
    pymol.localopt('uff')
    return pymol


def mopopt(mol_with_cat_mopin, name):
    """Optimize using catalyst using mopac."""
    with open(name + '.dat', 'wb') as f:
        f.write(mol_with_cat_mopin)    
    subprocess.call([MOPAC_LOCATION, name + '.dat'])
    if fileutil.file_exist(name + '.arc'):
        pymol = pybel.readfile('mopout', name + '.out').next()
        os.unlink(name + '.dat')
        os.unlink(name + '.arc')
        os.unlink(name + '.out')
        return pymol
    else:
        return None


def create_catalyst_txt(args):
    """Creates string with oh down and oh up catalysts."""
    m, x, times = args.catalyst
    m_lane_OH_up = GENERAL_CATALYST_M_LANE.format(M=m, M_dist=args.m_dist, 
                                                  M_C_angle=args.m_c_angle, 
                                                  M_H_dihed='{M_H_dihed}')
    m_lane_OH_down = GENERAL_CATALYST_M_LANE.format(M=m, M_dist=args.m_dist, 
                                                    M_C_angle=args.m_c_angle, 
                                                    M_H_dihed='{M_H_dihed}')                                                  
    x_lanes = '\n'
    angle_increment = 360/float(times)
    for i in range(int(times)):
        x_c_dihed = args.x_c_dihed + angle_increment*i
        x_lane = GENERAL_CATALYST_X_LANE.format(X=x, X_dist=args.x_dist,
                                                X_O_angle=args.x_o_angle,
                                                X_C_dihed=x_c_dihed)
        x_lanes += x_lane + '\n'
    cat_oh_up = m_lane_OH_up + x_lanes
    cat_oh_d = m_lane_OH_down + x_lanes
    return cat_oh_up, cat_oh_d


def write_gau_TS_file(folder_name, pypair, name, ts_kwards, args):
    """Writes gaussian input TS file"""
    inp = gaussian.TS_input()
    inp.provide_pymol(pypair[0], 1)
    inp.provide_pymol(pypair[1], 2)
    inp.set_proc(args.proc)
    inp.set_memory(args.memory)   
    inp.set_chk(name)
    inp.set_kwards(ts_kwards)
    inp.set_description(name, 1)    
    inp.set_description(name, 2)    
    inp.set_charge(args.charge, 1)
    inp.set_charge(args.charge, 2)
    inp.set_mult(args.mult, 1)
    inp.set_mult(args.mult, 2)
    inp.add_custom_basis_from_dic(BASIS[args.basis])
    with open(folder_name + '/' + name + '.com', 'wb') as f:
        inp.write(f)


def write_gau_min_file(folder_name, pymol, name, min_kwards, args):
    """Writes gaussian input min file"""
    inp = gaussian.Input()
    inp.provide_pymol(pymol)
    inp.set_proc(args.proc)
    inp.set_memory(args.memory)   
    inp.set_chk(name)
    inp.set_kwards(min_kwards)
    inp.set_description(name)    
    inp.set_charge(args.charge)
    inp.set_mult(args.mult)
    inp.add_custom_basis_from_dic(BASIS[args.basis])
    with open(folder_name + '/' + name + '.com', 'wb') as f:
        inp.write(f)


def optimize_pair(folder_name, function, pair, name, ts_kwards, args):
    """Optimizes pair of files with given optimization function, which 
    accepts mopin string and molecules name and returns pymols."""
    opt_pair = [function(pair[i], name + str(i)) for i in range(len(pair))]
    try:
        write_gau_TS_file(folder_name, opt_pair, name, ts_kwards, args)
    except AttributeError:
        print name + " wasn't optimized with MOPAC"


def optimize_min(folder_name, function, mol, name, min_kwards, args):
    """Optimizes molecule with given optimization function, which 
    accepts mopin string and molecules name and returns pymols."""
    opt_pair = function(mol, name)
    try:
        write_gau_min_file(folder_name, opt_pair, name, min_kwards, args)
    except AttributeError:
        print name + " wasn't optimized with MOPAC"    


def write_min(folder_name, min_type, mol_mopin, min_kwards, args):
    """Writes minimization file optimized with mopac with metal hydrogen
    dihedral = MOPIN_HORIS_ANGLE."""
    name = min_type + '_min_' + ''.join(args.catalyst) + '_' + args.solvent + '_' + args.functional + '_' + args.basis
    mol_mopin_horisont = mol_mopin.format(M_H_dihed=MOPIN_HORIS_ANGLE)
    if not args.nomopac:
        t = Thread(target=optimize_min, args=(folder_name, mopopt, mol_mopin_horisont, name,  min_kwards, args))
        t.start()
    else:
        pymol = pybel.readstring('mopin', mol_mopin_horisont)
        write_gau_min_file(folder_name, pymol, name, min_kwards, args)
        
    
def write_TSs(folder_name, TS, mols_pair_mopin, ts_kwards, args):
    """Writes TS's input files optimized with mopac, with UFF and for 
    different dihedrals, given in input."""
    base_name = TS + '_TS_' + ''.join(args.catalyst) + '_' + args.solvent + '_' + args.functional + '_' + args.basis
    mols_pair_mopin_vert = [mol.format(M_H_dihed=MOPOPT_VERT_ANGLE) for mol in mols_pair_mopin]
#    mols_pair_mopin_horisont = [mol.format(M_H_dihed=MOPIN_HORIS_ANGLE) for mol in mols_pair_mopin]
    if not args.nomopac:
        t1 = Thread(target=optimize_pair, args=(folder_name, mopopt, 
                                                mols_pair_mopin_vert, 
                                                base_name + '_mopopt_vert', 
                                                ts_kwards, args))
        t1.start()
#        t2 = Thread(target=optimize_pair, args=(folder_name, mopopt, 
#                                                mols_pair_mopin_vert, 
#                                                base_name + '_mopopt_hor', 
#                                                ts_kwards, args))
#        t2.start()
#    if args.localopt:
#        optimize_pair(folder_name, localopt, mols_pair_mopin_horisont, base_name + '_localopt_vert', ts_kwards, args)
#        optimize_pair(folder_name, localopt, mols_pair_mopin_horisont, base_name + '_localopt_hor', ts_kwards, args)
    for angle in args.m_h_diheds:
        mols_pair_mopin_angle = [mol.format(M_H_dihed=angle) for mol in mols_pair_mopin]
        mols_pair_angle = [pybel.readstring('mopin', mol) for mol in mols_pair_mopin_angle]
        angle_name = base_name + '_' + str(int(float(angle)))
        write_gau_TS_file(folder_name, mols_pair_angle, angle_name, ts_kwards, args)

        
def main(argv):
    args = process_command_line(argv)
    min_kwards = get_kwards('min', args)
    ts_kwards = get_kwards('ts', args)
    folder_name = make_folder(args)
    catalyst_OH_up_txt, catalyst_OH_down_txt = create_catalyst_txt(args)
    #minimums
    for minimum, file_spec in MIN_FILES.iteritems():
        if file_spec[1]:
            cat_txt = catalyst_OH_down_txt
        else:
            cat_txt = catalyst_OH_up_txt
        mol_mopin = add_catalyst(BASE_SET_LOCATION + '/' + file_spec[0],
                            BASE_FORMAT, args.mult, args.charge, cat_txt)
        write_min(folder_name, minimum, mol_mopin, min_kwards, args)
    #TSs
    for TS, input_pair in TS_PAIRS.iteritems():
        if input_pair[2]:
            cat_txt = catalyst_OH_down_txt
        else:
            cat_txt = catalyst_OH_up_txt
        mols_pair_mopin = [add_catalyst(BASE_SET_LOCATION + '/' + input_pair[i],
                            BASE_FORMAT, args.mult, args.charge, cat_txt) \
                            for i in range(len(input_pair) - 1)]
        write_TSs(folder_name, TS, mols_pair_mopin, ts_kwards, args)
        
 
if __name__ == '__main__':
    main(sys.argv[1:])
