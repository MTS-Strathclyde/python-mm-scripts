#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 13 11:15:08 2016

@author: max
"""

import numpy as np
import sys
import subprocess
import argparse
import asa
import os
import shlex
import pybel


def process_command_line(argv):
    """Processes arguments

    Parameters
    ----------
    argv : list
        Command line arguments.

    Returns
    -------
    out : argparse.Namespace
        Namespace of command line arguments.
    """
    parser = argparse.ArgumentParser(description=""" Compute ks using areas.""")
    #Positional args
    parser.add_argument('file', metavar='molec.pdb',
                        help="""Input file.""")
    return parser.parse_args(argv)


def run_ffld_server(args, name):
    p, no_p_name = os.path.split(name)
    if p == '':
        p = '.'
    schrod_path = os.environ['SCHRODINGER']
    command = """{}/utilities/ffld_server -ipdb {}.pdb -print_parameters"""\
                .format(schrod_path,name)
    out = subprocess.check_output(shlex.split(command)).splitlines()
    # check if ff parms were generated 
    parm_data_start = None
    for i, l in enumerate(out):
        if l.startswith('OPLSAA FORCE FIELD TYPE ASSIGNED'):
            parm_data_start = i + 4
    if not parm_data_start:
        # try to use mol2 for generation
        print('Failed to assigne parameters for pdb, trying mol2')
        # gen mol2 using babel - it recognizes atoms
        mol = pybel.readfile('pdb', name + '.pdb').next()
        mol_mol2 = mol.write('mol2')
        # in babel mol2 there is a problem atoms called CLx are regarded as
        # carbons, not as clorines
        fixed_mol2 = []
        for l in mol_mol2.splitlines():
            ls = l.split()
            if len(ls) == 9: # atom row
                if ls[1].startswith('CL') or ls[1].startswith('Cl'):
                    ls[5] = 'Cl'
                    l = '{:>7}  {:<7}{:>10}{:>10}{:>10} {:<8}{:<3}{:<8}{:>10}'.format(*ls)
                if ls[1].startswith('BR') or ls[1].startswith('Br'):
                    ls[5] = 'Br'
                    l = '{:>7}  {:<7}{:>10}{:>10}{:>10} {:<8}{:<3}{:<8}{:>10}'.format(*ls)                
            fixed_mol2.append(l)
        fixed_mol2 = '\n'.join(fixed_mol2)
        with open(name + '.mol2', 'w') as f:
            f.write(fixed_mol2)
        command = """{}/utilities/ffld_server -imol2 {}.mol2 -print_parameters"""\
                    .format(schrod_path,name)
        out = subprocess.check_output(shlex.split(command)).splitlines()
        # check again
        for i, l in enumerate(out):
            if l.startswith('OPLSAA FORCE FIELD TYPE ASSIGNED'):
                parm_data_start = i + 4
        if not parm_data_start:
            raise ValueError('Failed to assign oplsaa parameters to {}'.format(name))
    return out, parm_data_start


def get_opls_parameters(args, name):
    out, parm_data_start = run_ffld_server(args, name)
    radii = []
    epss = []
    chgs = []    
    for l in out[parm_data_start:]:
        if not l.startswith('-----'):
            l = l.split()
            radii.append(float(l[5])/2*2**(1./6))   # rmin/2, A
            epss.append(float(l[6]))                # kcal/mol 
            chgs.append(float(l[4]))                # e
        else:
            break
    return radii, epss, chgs


def get_ks(eps):
    # returns Ks for particular eps in (kcal/mol/A^2*M^(-1))
    A = 0.00402488643382
    b = -0.649969253273
    return A*np.exp(b*eps)

    
def get_ks2(eps, chg):
    A, b, d, c1, c2 = [ 0.00499462, -0.75687269,  0.0484793,  -0.00396178, -1.82300813]
    return A*np.exp(b*eps) +  d*chg**2 + c1*np.log(eps)*(chg+c2*chg**2)


def main(argv):
    args = process_command_line(argv)
    k = 8.314/4.184/1000 # kcal/mol/K
    T = 298.15  # K
    name = args.file[:-4]
    radii, epss, chgs = get_opls_parameters(args, name)
    with open(args.file) as f:
        pdb_txt = f.read()
    areas = np.array(asa.main(pdb_txt, 960))   # A^2
    #ks_s = np.array([get_ks(eps) for eps in epss])  # kcal/mol/A^2/M
    ks_s = np.array([get_ks2(eps,chg) for eps,chg in zip(epss,chgs)])  # kcal/mol/A^2/M
    dG_per_conc = np.sum(ks_s*areas)
    print dG_per_conc/k/T*np.log10(np.e)  # 1/M


if __name__ == '__main__':
    main(sys.argv[1:])