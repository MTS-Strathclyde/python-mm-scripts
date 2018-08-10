#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu May 15 15:31:07 2014

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
import numpy as np

import charge_extractor


REQUIRED_EXECUTABLES = ['antechamber', 'parmchk', 'tleap', 'rism3d.snglpnt']

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

ng_10u_c_np={ng_u_cor_with_kh} kcal/mol

o_np_dcfi_v_in={0[0][0]}
o_np_dcfi_v_out={0[1][0]}
o_np_dcfi_times_h_v_in={0[0][1]}
o_np_dcfi_times_h_v_out={0[1][1]}
o_np_gf_v_in={0[0][2]}
o_np_gf_v_out={0[1][2]}
"""

C_O_DX = 'c_{}.O.1.dx'
C_ONP_DX = 'c_{}_0chg.O.1.dx'
H_ONP_DX = 'h_{}_0chg.O.1.dx'
H_O_DX = 'h_{}.O.1.dx'
U_O_DX = 'u_{}.O.1.dx'
C_HNP_DX = 'c_{}_0chg.H1.1.dx'
H_HNP_DX = 'h_{}_0chg.H1.1.dx'


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
    

def gamma_cor(integral, T, gamma=0.38):
    t = float(T)
    density = water_concentration(t)*6.0221413E-4  #molecules / A^3
    R = 1.9872041E-3 #kcal/K/Mol
    return R*t*density/2*integral*(1-gamma)


def uc(gf, pmv, t):
    """Return universal correction."""
    gf = float(gf)
    pmv = float(pmv)
    t = float(t)
    density = water_concentration(t)*6.0221413E-4  #molecules / A^3
    return gf - 3.2217*density*pmv + 0.5783  #kcal/mol


def water_succ_rism_script(T):
    """Create water susceptibility script for given temperature T in K."""
    diel = round(water_dielectric_const(T), 3)
    conc = round(water_concentration(T), 3)
    return SOLV_SUCEPT_SCRPT.format(temp=T, diel=diel, conc=conc)
    
    
def two_dx_file_correction(integrable_dx_m, contstraint_dx_m,
                           spacing3, constraint_value, bigger=True):
    """
    Integrate integrable_dx matrix over the region where
    constraint_dx matrix <= constraint_value. (<= if bigger=False)
    spacing3 = dx*dy*dz"""
    if bigger:
        res = np.where(contstraint_dx_m >= constraint_value, integrable_dx_m
                       , np.NAN).flatten()
    else:
        res = np.where(contstraint_dx_m <= constraint_value, integrable_dx_m
                       , np.NAN).flatten()
    res = res[np.logical_not(np.isnan(res))]
    return sum(res*spacing3)
    

def load_dx(fname):
    """Return dx matrix, origin coordinates tuple and spacing between points tuple."""
    with open(fname, 'rb') as f:
        txt = f.readlines()
    if txt[0].startswith('object 1 class gridpositions counts'):
        size = txt[0][35:].split()
        size = map(float, size)
    else:
        print 'Wrong file format'
        raise ValueError
    dx_m = []
    for line in txt:
        ln = line.split()
        if len(ln) == 3:
            dx_m.extend(ln)
        elif (ln[0] == 'origin'):
            OrX=float(ln[1])
            OrY=float(ln[2])
            OrZ=float(ln[3])
        elif (ln[0]=='delta' and ln[1]!='0'):
            dX = float(ln[1])
        elif (ln[0]=='delta' and ln[2]!='0'):
            dY = float(ln[2])
        elif (ln[0]=='delta' and ln[3]!='0'):
            dZ = float(ln[3])
    dx_m = np.array(dx_m)
    dx_m = dx_m.reshape(size)
    return dx_m.astype('float'), (OrX, OrY, OrZ), (dX, dY, dZ)    


def calculate_ng_corrections(name, T):
    """Calculate 12 Ng corrections:
    """
    T = float(T)
    p, no_p_name = os.path.split(name)
    c_np_dx, _, spacing_l = load_dx(os.path.join(p, C_ONP_DX.format(no_p_name)))
    h_onp_dx = load_dx(os.path.join(p, H_ONP_DX.format(no_p_name)))[0]    
    u_dx = load_dx(os.path.join(p, U_O_DX.format(no_p_name)))[0]
    spacing3 = reduce(lambda x, y: x*y, spacing_l)
    kt_2_kcal = 1.9872041E-3*T  #kt/kcal/mol
    # produces a list of tuples with corrections,
    # [(ng10unp, ng10u), (ng5unp, ng5u), (ng2unp, ng2u)]
    ng_u_cor = two_dx_file_correction(c_np_dx, u_dx, spacing3, 10*kt_2_kcal)
    # [(ng099hnp, ng099h), (ng09hnp, ng09h), (ng075hnp, ng075h)]
    # hydrogen and oxygen corellation function integrals
    # v_in is defined as region where h < 0.9
    # [(-c_v_in, -c*h_v_in, -c-0.5c*h_v_in), (same, but v_out)]
    gf_o = [(two_dx_file_correction(-c_np_dx, u_dx, spacing3, 10*kt_2_kcal, i),
             two_dx_file_correction(-c_np_dx*h_onp_dx, u_dx, spacing3, 10*kt_2_kcal, i),
             two_dx_file_correction(-c_np_dx-0.5*h_onp_dx*c_np_dx, u_dx, spacing3, 10*kt_2_kcal, i))
             for i in [False, True]]
    return ng_u_cor, gf_o
    

class Prmtop_file(object):
    def __init__(self, prmtop_name):
        self.prmtop_name = prmtop_name
        with open(prmtop_name, 'rb') as f:
            self.prm_lines = f.readlines()
        
    def set_charges(self, charge_l=None):
        """If charge_l = None, sets charges to 0"""
        i = 0
        while i < len(self.prm_lines):
            line = self.prm_lines[i]
            if line.startswith('%FLAG CHARGE'):
                if self.prm_lines[i+1].startswith('%FORMAT(5E16.8)'):
                    j = i + 2
                    if not charge_l:
                        while self.prm_lines[j].startswith(' '):
                            chrgs = self.prm_lines[j].split()
                            new_chrgs = '  0.00000000E+00'*len(chrgs)
                            self.prm_lines[j] = new_chrgs + '\n'
                            j += 1
                        break
                    else:
                        while self.prm_lines[j].startswith(' '):
                            chrgs = self.prm_lines[j].split()
                            new_chrgs = ['{: .8E}'.format(charge_l.pop(0)) for chg in chrgs]
                            self.prm_lines[j] = ' ' + ' '.join(new_chrgs) + '\n'
                            j += 1
                        break
                else:
                    raise ValueError('Charge given in unknown format')
            i += 1
    
    def write_prmtop(self, name=None):
        if not name:
            name = self.prmtop_name
        with open(name, 'wb') as f:
            f.writelines(self.prm_lines)


class RISM3D_Singlpnt(object):
    RISM3D_GRID = '0.300000,0.300000,0.300000'
    RISM3D_TOLLERANCE = '0.0000000001'        
    def __init__(self, name, T, logfile, prmtop_name=None, xvv_name=None,
                 start_time=None):
        self.name = name
        self.T = T
        self.p, self.no_p_name = os.path.split(name)
        self.logfile = logfile
        if prmtop_name:
            self.prmtop_name = prmtop_name
        else:
            self.prmtop_name = '{}.prmtop'.format(self.no_p_name)
        if xvv_name:
            self.xvv_name = xvv_name
        else:
            self.xvv_name = 'water_{temp}.xvv'.format(temp=self.T)
        if start_time:
            self.start_time = time.time()
        else:
            self.start_time = start_time
        self.run_flags_list = None
    
    def setup_calclation(self, closure='kh', calculate_h=True, calculate_c=True,
                         calculate_u=True, buffer_distance=30, grdspc=None,
                         tollerance=None, polar_decomp=False):
        """If grdspc and tollerance are set to None method will use default
        values."""
        if not grdspc:
            grdspc = self.RISM3D_GRID
        if not tollerance:
            tollerance = self.RISM3D_TOLLERANCE
        self.run_flags_list = ['rism3d.snglpnt',
             '--pdb', '{}.pdb'.format(self.no_p_name),
             '--prmtop', self.prmtop_name,
             '--closure', closure,
             '--xvv', self.xvv_name,                     
             '--buffer', '{:.6f}'.format(buffer_distance), #distance between solute 
                                                      #and the edge of solvent box
             '--grdspc', grdspc,
             '--tolerance', tollerance]
        if calculate_h:
            self.run_flags_list.extend(['--huv', 
                                         'h_{}'.format(self.no_p_name)])
        if calculate_c:
            self.run_flags_list.extend(['--cuv', 
                                         'c_{}'.format(self.no_p_name)])
        if calculate_u:
            self.run_flags_list.extend(['--uuv', 
                                         'u_{}'.format(self.no_p_name)])
        if polar_decomp:
            self.run_flags_list.extend(['--polarDecomp'])
            
    def run_calculation_and_log(self):
        rism_out = subprocess.check_output(self.run_flags_list, cwd=self.p)
        self.logfile.write(rism_out)
        self.logfile.flush()
        #write timestamp and close
        end_time = time.time()
        self.logfile.write(str(datetime.datetime.now()) + '\n')
        runtime = end_time - self.start_time
        self.logfile.write(str(round(runtime)))
        self.logfile.flush()
        self.logfile.close()
        
        
def run_3drism_0chrg(name, T):
    log_name = '{}_0chg.log'.format(name)
    logfile = open(log_name, 'wb')
    # Wirte timestamp and start measuring runtime
    start_time = time.time()
    logfile.write(str(datetime.datetime.now()))
    logfile.flush()    
    #prepare files
    mol_pdb = '{}.pdb'.format(name)
    mol_0chg_pdb = '{}_0chg.pdb'.format(name)
    shutil.copy(mol_pdb, mol_0chg_pdb)
    prmtop_name = name + '.prmtop'
    name_0chg = name + '_0chg'
    prmtop_0chg_name = name + '_0chg.prmtop'    
    chg0_prmtop = Prmtop_file(prmtop_name)
    chg0_prmtop.set_charges()
    chg0_prmtop.write_prmtop(prmtop_0chg_name)

    _, prmtop_0chg_name = os.path.split(prmtop_0chg_name)

    rism_calc = RISM3D_Singlpnt(name_0chg, T, logfile,
                                prmtop_name=prmtop_0chg_name,
                                start_time=start_time)
    rism_calc.setup_calclation(calculate_u=False)
    rism_calc.run_calculation_and_log()


def prepare_3drism_calc(name, T=298.15):
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
    prmtop_name = '{}.prmtop'.format(name)
    return logfile, prmtop_name, start_time


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


def post_processing(name, T, del_dx_files=True):
    p, no_p_name = os.path.split(name)
    log_name = '{}.log'.format(name)
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
    ng_u_cor, gf_o = calculate_ng_corrections(name, T)
    ng_u_cor_with_kh = kh + gamma_cor(ng_u_cor, T)

    results = RESULTS.format(gf_o, ng_u_cor_with_kh=ng_u_cor_with_kh,
                             kh=kh, gf=gf, pmv=pmv, ucorr=ucorr)
    #Write results
    with open(p + '/results.txt', 'wb') as f:
        f.write(results)
    if del_dx_files:
        delete_dx(name)
    return kh, gf, pmv, ucorr
    

def calculate_mol(mol_path, gamess_charges_f, T=298.15, delete_dx_files=True):
    """Run RISM calculation and delete dx files."""
    T = round(float(T), 2)
    if 253.15 <= T <= 383.15:
        try:
            name = prepare_calc_directory(mol_path, T)
            #kh, gf, pmv, ucorr = run_3drism(name, T)
        except OSError:
            print 'Skipping job {} T = {}'.format(mol_path, T)
            return 1
        charge_l = charge_extractor.extract_charges(gamess_charges_f)  #in electron charge
        charge_l = map(lambda x: x*18.2223, charge_l) #convert to strange amber format
        logfile, prmtop_name, start_time = prepare_3drism_calc(name, T)
        custom_prmtop = Prmtop_file(prmtop_name)
        custom_prmtop.set_charges(charge_l)
        custom_prmtop.write_prmtop()
        rism_calc = RISM3D_Singlpnt(name, T, logfile, start_time=start_time)
        rism_calc.setup_calclation()
        rism_calc.run_calculation_and_log()
        run_3drism_0chrg(name, T)
        post_processing(name, T)
        if delete_dx_files:
            delete_dx(name)
    else:
        print 'Temperature outside of scripts range!'
    

def main(argv):
    calculate_mol(argv[0], argv[1], argv[2])

if __name__ == '__main__':
    main(sys.argv[1:])



    