#!/usr/bin/env python

# -*- coding: utf-8 -*-
"""
Created on Wed May 28 12:49:25 2014

@author: max
"""

import os
import shutil
import sys
from cookbook import utility

RISM3D_GRID = '0.300000,0.300000,0.300000'
RISM3D_TOLLERANCE = '0.0000000001'        
TIMEOUT = 45*60

RESULTS = """dGhyd(KH)= {kh} kcal/mol
dGhyd(GF)= {gf} kcal/mol
PMV= {pmv} AA^3
"""


def run_rism(name, T):
        #Run 3DRISM
    log_name = '{}.log'.format(name)
    logfile = open(log_name, 'wb')
    p, no_p_name = os.path.split(name)
    if p == '':
        p = '.'
    print 'Starting 3DRISM for {} at T={}'.format(name, T)
    xvv_name = 'water_{temp}.xvv'.format(temp=T)    
    rism3d_command = ['rism3d.snglpnt',
                     '--pdb', '{}.pdb'.format(no_p_name),
                     '--prmtop', '{}.prmtop'.format(no_p_name),
                     '--closure', 'hnc',
                     '--guv', 'g_{}'.format(no_p_name), #root name for solvent pair  
                                                   #distribution files
                     '--cuv', 'c_{}'.format(no_p_name), #root name for solvent direct  
                                                   #correlation files
                     '--xvv', xvv_name,
                     '--buffer', '30.000000', #distance between solute and the
                                              #edge of solvent box
                     '--grdspc', RISM3D_GRID,
#                     '--polarDecomp',      # decompose solvation FE
#                                          # into polar and non-polar parts
                     '--tolerance', RISM3D_TOLLERANCE]
    stdout, stderr = utility.RunCmd(rism3d_command, TIMEOUT, cwd=p).Run()
    if stdout:
        logfile.write(stdout)
        logfile.flush()
    if stderr:
        logfile.write(stdout)
        logfile.flush()        
    # Parse rism output
    kh, gf, pmv = None, None, None
    with open(log_name, 'rb') as f:
        for line in f:
            if line[0:11] == "rism_exchem":
                kh = float(line.split()[1])
            if line[0:11] == "rism_exchGF":
                gf = float(line.split()[1])
            if line[0:11] == "rism_volume":
                pmv = float(line.split()[1])
    results = RESULTS.format(kh=kh, gf=gf, pmv=pmv)
    logfile.write(results)
    logfile.flush()    
    #Write timestamp and runtime
    with open(p + '/results.txt', 'wb') as f:
        f.write(results)


    

def main(pdb, prmtop, temps):
    """Calculate uc and ngb given folders with temperature and xvv files in them."""
    for t in temps:
        shutil.copy(pdb, t)
        shutil.copy(prmtop, t)
        os.chdir(t)
        name = pdb[:-4]
        run_rism(name, t)
        os.chdir('..')


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], sys.argv[3:])






