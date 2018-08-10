# -*- coding: utf-8 -*-
"""
Created on Fri May 23 09:45:00 2014

@author: max
"""

import os
import glob
import time
import datetime
import subprocess
import shutil

import charge
import water


RUNLEAP = """source leaprc.gaff
loadamberprep {name}.prepin
check MOL
loadamberparams {name}.frcmod
SaveAmberParm MOL {name}.prmtop {name}.incrd
SavePdb MOL {name}.pdb
quit
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
    

        
def prepare_3drism_calc(name, charge_model_l, T=298.15, closure='KH',
                        wmodel='SPC'):
    """Run 3DRISM calculation at given temperature.
    
    name is the name of a pdb file containing solute without extension.
    Function will automatically generate an xvv solvent file for given
    temperature and then use it to preform a rism calculation.
    
    Returns (kovalenko-hirata, gaussian fluctation hydration energy,
             molecular volume, uc-corrected hydration energy)
    
    """
    p, no_p_name = os.path.split(name)
    log_name = '{}/pre_{}.log'.format(p, no_p_name)
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
    succ_srcirpt = water.water_succ_rism_script(T, closure, wmodel)
    with open(xvv_script_name, 'wb') as f:
        f.write(succ_srcirpt)
    xvv_out = subprocess.check_output(['bash', xvv_script_name_no_p], cwd=p)
    logfile.write(xvv_out)
    logfile.flush()
    #write timestamp and close
    end_time = time.time()
    logfile.write(str(datetime.datetime.now()) + '\n')
    runtime = end_time - start_time
    logfile.write(str(round(runtime)))
    logfile.flush()
    logfile.close()
    prmtop_names = charge.generate_prmtops(name, charge_model_l)
    return prmtop_names


def prepare_calc_directory(mol_path, T):
    """mol_path is pdb filename *with* extension.
    Create new calcuation directory and copy pdb file there.
    Return path of molecule in directory without extension."""
    pdb_path, name_without_path = os.path.split(mol_path)
    dir_name = os.path.join(pdb_path, str(T))
    #dir_name = make_uq_dir(dir_name) we don't want unique names
    os.mkdir(dir_name)
#    try:
#        os.mkdir(dir_name)
#    except OSError, e:
#        if e.errno == 17: #directory already exists, all is well
#            pass
#        else:
#            raise e
    print dir_name
    print name_without_path
    name = dir_name + '/' + name_without_path    
    shutil.copy(mol_path, name)
    return name[:-4]
    
    