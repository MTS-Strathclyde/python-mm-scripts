#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 30 09:31:23 2015

@author: max
"""
import numpy as np
import scipy.stats
import os
import sys
import rism3d_pc
import argparse
import csv


SETSCHENOW_OUT = """
k_conc = {:.4f}

intercept = {:.4f}
r^2 = {:.4f}
p-value = {:.4f}
std_err = {:.4f}
"""


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
    parser = argparse.ArgumentParser(description="""Run a series of 3D-RISM 
            calculations at a range of salt concentration to determine
            setschenows constant.""")
    #Positional args
    parser.add_argument('file', metavar='molec.pdb',
                        help="""Input file. Must be in pdb format
                        acceptable by Antechamber. Must have pdb
                        extension.""")
    parser.add_argument('prmtop',
                        help=""" Existing prmtop file.""")
    parser.add_argument('xvv_dir',
                        help=""" Directory with existing xvv files,
                        thermodynamic outputs and submission scripts.
                        Filenames should have following format:
                        x{conc}_y.format .""")
    #Optional args
    parser.add_argument('--extra_args',
                        help="""Extra arguments for rism3d_pc.py surrounded
                        by quotes""", default='')
    return parser.parse_args(argv)



#def get_total_density(script_name):
#    """ Returns total species density read from the 1drism job script
#    as long as all spcecies concentrations are in M."""
#    M_to_mol_per_angstrom3 = 6.0221413E-4
#    with open(script_name) as f:
#        txt = f.readlines()
#    total_density = 0
#    for l in txt:
#        if 'DENSITY=' in l:
#            dens = re.findall(r'DENSITY=(\d*\.\d*)', l)[0]
#            total_density += float(dens)
#    return total_density*M_to_mol_per_angstrom3


def analyze_xvv_directory(dirname):
    xvv_files = [os.path.join(dirname, f) for f in os.listdir(dirname) if f.endswith('xvv')]
    #print xvv_files
    return xvv_files
    
    
def main(argv):
    args = process_command_line(argv)
    xvv_files = analyze_xvv_directory(args.xvv_dir)
    iscs = []
    for xvv in xvv_files:
        xvv_name = os.path.split(xvv)[1]
        conc = float(xvv_name.split('_')[1])
        isc  = rism3d_pc.main("{} -p {} --xvv {} --dir_name conc_{} {}"\
            .format(args.file, args.prmtop, xvv, conc, args.extra_args).split())
        iscs.append((float(conc),isc))
    iscs = np.array(iscs)
    iscs = iscs[np.argsort(iscs[:,0])]
    # convert from dG to log(H_0/H_i)
    log_ratio = (iscs[:,1] - iscs[0,1])*4184/(298.15*8.314)  
    log_ratio = log_ratio*np.log10(np.e) # convert to log base 10
    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(iscs[:,0], log_ratio)
    with open('setschenow.txt', 'w') as f:
        f.write(SETSCHENOW_OUT.format(slope, intercept, r_value,p_value, std_err))
    with open('setschenow.csv', 'w') as f:
        wr = csv.writer(f)
        wr.writerow(['Conc [M]', 'dG', 'log_ratio'])
        data  = np.c_[iscs, log_ratio]
        for conc, dg, log_r in data:
            wr.writerow([conc, dg, log_r])        
    print 'Constant = {:.4f}'.format(slope)
    
    
if __name__=='__main__':
    main(sys.argv[1:])
    
    
    
    
    