#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 30 16:04:49 2014

@author: max

Some variables are used extensively through the script:
name - name of the pdb molecule without extension
T or temp - temperature of the calculation
diel - dielectric constant
conc - concentration of water

Script creates a lot of files; some of them are:
{name}.log - logfile
water_{temp}_script.sh - water succ file generation script
water_{temp}.xvv - water succ file



"""
import sys
import subprocess
import shutil
import os
import glob
import datetime
import time


REQUIRED_EXECUTABLES = ['antechamber', 'parmchk', 'tleap', 'rism3d.snglpnt']

RISM3D_GRID = '0.300000,0.300000,0.300000'

RISM3D_TOLLERANCE = '0.0000000001'

SOLV_SUCEPT_SCRPT = """#!/bin/csh -f

cat > water_{temp}.inp <<EOF
&PARAMETERS
	THEORY='DRISM', CLOSUR='KH',           !Theory
	NR=16384, DR=0.025,                    !Grid
	OUTLST='xCGT', routup=384, toutup=0,   !Output
	NIS=20, DELVV=0.3, TOLVV=1.e-12,       !MDIIS
	KSAVE=-1, KSHOW=1, maxstep=10000,      !Check pointing and iterations
	SMEAR=1, ADBCOR=0.5,                   !Electrostatics
	TEMPER={temp}, DIEps={diel},           !bulk solvent properties
	NSP=1
/
	&SPECIES                               !SPC water
	DENSITY={conc}d0,
	MODEL="$AMBERHOME/dat/rism1d/model/SPC.mdl"
/
EOF

rism1d water_{temp} > water_{temp}.out || goto error

"""

RUNLEAP = """source leaprc.gaff
loadamberprep {name}.prepin
check MOL
loadamberparams {name}.frcmod
SaveAmberParm MOL {name}.prmtop {name}.incrd
SavePdb MOL {name}.pdb
quit
"""

RESULTS = """dGhyd(KH)= {kh} kcal/mol
dGhyd(GF)= {gf} kcal/mol
PMV= {pmv} AA^3
dGhyd(UC)= {ucorr} kcal/mol
"""

def make_uq_dir(name):
    try:
        os.mkdir(name)
        return name
    except OSError:
        i = 1
        while True:
            try:
                dir_name = name + '_{}'.format(i)                
                os.mkdir(dir_name)
                return dir_name
            except OSError:
                i += 1

def water_dielectric_const(T):
    """Return water dielectric constant for temperature 253.15K < T < 383.15K.
    
    Uses correlating equation (eq. 9) for static dielectri constant found in
    the doucment by The International Association for the Properties of 
    Water and Steam from 2011
    (http://www.iapws.org/relguide/LiquidWater.pdf)
    Pressure = 0.1 MPa
    
    >>> round(water_dielectric_const(273.15), 3)
    87.927
    >>> round(water_dielectric_const(298.15), 3)
    78.375
    >>> round(water_dielectric_const(375), 3)
    55.266
    """
    T_star = T/300.0
    coefs = [-43.7527, 299.504, -399.364, 221.327]
    exp_f = [-0.05, -1.47, -2.11, -2.31]
    e = 0
    for i in range(4):
        e += coefs[i]*T_star**(exp_f[i])
    return e
    
    
def water_concentration(T):
    """Return water concentration for temperature range 253.15K < T < 383.15K.
    
    Uses correlating equation (eq. 2) for specific volume found in
    the doucment by The International Association for the Properties of 
    Water and Steam from 2011
    (http://www.iapws.org/relguide/LiquidWater.pdf)
    Pressure = 0.1 MPa    
    
    >>> round(water_concentration(273.15), 3)
    55.498
    >>> round(water_concentration(298.15), 3)
    55.343
    """
    p0 = 10.0**5    # Pa
    R = 8.31464     # J/mol/K
    Tr = 10.0
    Ta = 593.0
    Tb = 232.0
    a = [1.93763157E-2,
         6.74458446E+3,
        -2.22521604E+5,
         1.00231247E+8,
        -1.63552118E+9,
         8.32299658E+9]
    b = [5.78545292E-3,
        -1.53195665E-2,
         3.11337859E-2,
        -4.23546241E-2,
         3.38713507E-2,
        -1.19946761E-2]
    n = [None, 4., 5., 7., 8., 9.]
    m = [1., 2., 3., 4., 5., 6.]
    def alpha(T): 
        return Tr/(Ta - T)
    def beta(T):
        return Tr/(T - Tb)
    coef = a[0] + b[0]*beta(T)**m[0]
    for i in range(1, 6):
        coef += a[i]*alpha(T)**n[i] + b[i]*beta(T)**m[i]
    v0 = R*Tr/p0*coef  # m3/mol
    return 1/(v0*1000)    # mol/L
    

def uc(gf, pmv, t):
    """Return universal correction."""
    gf = float(gf)
    pmv = float(pmv)
    t = float(t)
    density = water_concentration(t)*6.0221413E-4  #molecules / A^3
    return gf - 3.2217*density*pmv + 0.5783


def water_succ_rism_script(T):
    """Create water susceptibility script for given temperature T in K."""
    diel = round(water_dielectric_const(T), 3)
    conc = round(water_concentration(T), 3)
    return SOLV_SUCEPT_SCRPT.format(temp=T, diel=diel, conc=conc)
    

def run_3drism(name, T=298.15):
    """Run 3DRISM calculation at given temperature.
    
    name is the name of a pdb file containing solute without extension.
    Function will automatically generate an xvv solvent file for given
    temperature and then use it to preform a rism calculation.
    
    Returns (kovalenko-hirata, gaussian fluctation hydration energy,
             molecular volume, uc-corrected hydration energy)
    
    """
    p, no_p_name = os.path.split(name)
    log_name = '{}.log'.format(name)
    logfile = open(log_name, 'wb')
    # Wirte timestamp and start measuring runtime
    start_time = time.time()
    logfile.write(str(datetime.datetime.now()))
    logfile.flush()
    #Firstly we use antechamber to recognize atom and bonding types, and 
    #generate topology
    ante_out = subprocess.check_output(['antechamber', 
                     '-i', '{}.pdb'.format(no_p_name),
                     '-fi', 'pdb', 
                     '-o', '{}.prepin'.format(no_p_name), #output file
                     '-fo', 'prepi',   #output format describing each residue
                     '-c', 'bcc',      #charge method  (AM1-BCC)
                     '-s', '2',    #status info ; 2 means verbose
                     '-nc', '0',   #Net molecule charge
                     '-m', '1'],   #Multiplicity
                     cwd=p)
    logfile.write(ante_out)
    #Run parmchk to generate missing gaff force field parameters
    parm_out = subprocess.check_output(['parmchk',
                     '-i', '{}.prepin'.format(no_p_name),
                     '-f', 'prepi',
                     '-o', '{}.frcmod'.format(no_p_name)], #file with missing FF params
                     cwd=p)
    logfile.write(parm_out)
    logfile.flush()
    #Run tleap to generate topology and coordinates for the molecule
    leap_input_name = os.path.join(p, 'runleap.in')
    with open(leap_input_name, 'wb') as f:
        f.write(RUNLEAP.format(name=no_p_name))
    leap_out = subprocess.check_output(['tleap', '-f', 'runleap.in'], cwd=p)
    logfile.write(leap_out)
    logfile.flush()    
    #Generate water susceptibility file
    xvv_script_name_no_p = 'water_{}_script.sh'.format(T)
    xvv_script_name = os.path.join(p, xvv_script_name_no_p)
    succ_srcirpt = water_succ_rism_script(T)
    with open(xvv_script_name, 'wb') as f:
        f.write(succ_srcirpt)
    xvv_out = subprocess.check_output(['bash', xvv_script_name_no_p], cwd=p)
    logfile.write(xvv_out)
    logfile.flush()    
    #Run 3DRISM
    print 'Starting 3DRISM for {} at T={}'.format(name, T)
    xvv_name = 'water_{temp}.xvv'.format(temp=T)    
    rism_out = subprocess.check_output(['rism3d.snglpnt',
                     '--pdb', '{}.pdb'.format(no_p_name),
                     '--prmtop', '{}.prmtop'.format(no_p_name),
                     '--closure', 'kh',
                     '--guv', 'g_{}'.format(no_p_name), #root name for solvent pair  
                                                   #distribution files
                     '--cuv', 'c_{}'.format(no_p_name), #root name for solvent direct  
                                                   #correlation files
                     '--xvv', xvv_name,
                     '--buffer', '30.000000', #distance between solute and the
                                              #edge of solvent box
                     '--grdspc', RISM3D_GRID,
                     '--polarDecomp',      # decompose solvation FE
#                                          # into polar and non-polar parts
                     '--tolerance', RISM3D_TOLLERANCE], cwd=p)
    logfile.write(rism_out)
    logfile.flush()
    # Parse rism output
    kh, gf, pmv, ucorr = None, None, None, None
    with open(log_name, 'rb') as f:
        for line in f:
            if line[0:11] == "rism_exchem":
                kh = float(line.split()[1])
            if line[0:11] == "rism_exchGF":
                gf = float(line.split()[1])
            if line[0:11] == "rism_volume":
                pmv = float(line.split()[1])
    ucorr = uc(gf, pmv, T)
    results = RESULTS.format(kh=kh, gf=gf, pmv=pmv, ucorr=ucorr)
    logfile.write(results)
    logfile.flush()    
    #Write timestamp and runtime
    end_time = time.time()
    logfile.write(str(datetime.datetime.now()) + '\n')
    runtime = end_time - start_time
    logfile.write(str(round(runtime)))
    logfile.flush()
    logfile.close()
    with open(p + '/results.txt', 'wb') as f:
        f.write(results)
    return kh, gf, pmv, ucorr


def prepare_calc_directory(mol_path, T):
    """Create new calcuation directory and copy pdb file there.
    Return path of molecule in directory without extension."""
    pdb_path, name_without_path = os.path.split(mol_path)
    dir_name = os.path.join(pdb_path, str(T))
    #dir_name = make_uq_dir(dir_name) we don't want unique names
    os.mkdir(dir_name)
    print dir_name
    print name_without_path
    name = dir_name + '/' + name_without_path    
    shutil.copy(mol_path, name)
    return name[:-4]
    
    
def delete_dx(name):
    """Delete files with correlation functions."""
    p, _ = os.path.split(name)
    dx_files = glob.glob(os.path.join(p, '*.dx'))
    dx_files.extend(glob.glob(os.path.join(p, '*.cvv*')))
    dx_files.extend(glob.glob(os.path.join(p, '*.gvv*')))
    dx_files.extend(glob.glob(os.path.join(p, '*.xvv*')))
    dx_files.extend(glob.glob(os.path.join(p, '*.sav*')))
    for f in dx_files:
        os.unlink(f)
    

def run_rism(mol_path, T=298.15, delete_dx_files=True):
    """Run RISM calculation and delete dx files."""
    T = round(float(T), 2)
    if 253.15 <= T <= 383.15:
        try:
            name = prepare_calc_directory(mol_path, T)
            #kh, gf, pmv, ucorr = run_3drism(name, T)
        except OSError:
            print 'Skipping job {} T = {}'.format(mol_path, T)
            return 1
        run_3drism(name, T)
        if delete_dx_files:
            delete_dx(name)
    else:
        print 'Temperature outside of scripts range!'
    

def main(argv):
    run_rism(argv[0], argv[1])

if __name__ == '__main__':
    main(sys.argv[1:])



    