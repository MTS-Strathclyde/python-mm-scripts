#!/usr/bin/env python

# -*- coding: utf-8 -*-
"""
Created on Wed May 28 12:49:25 2014

@author: max
"""

import os
import numpy as np
import csv


from rism3d import rism_singlpnt as rism_singlpnt
from rism3d import water as water

PRMTOP = """%VERSION  VERSION_STAMP = V0001.000  DATE = 04/24/14  17:32:20                  
%FLAG TITLE                                                                     
%FORMAT(20a4)                                                                   
MOL                                                                             
%FLAG POINTERS                                                                  
%FORMAT(10I8)                                                                   
       1       1       1       0       0       0       0       0       0       0
       1       1       0       0       0       1       1       0       1       0
       0       0       0       0       0       0       0       0       3       0
       0
%FLAG ATOM_NAME                                                                 
%FORMAT(20a4)                                                                   
C1
%FLAG CHARGE                                                                    
%FORMAT(5E16.8)                                                                 
 {chg: .8E}
%FLAG ATOMIC_NUMBER                                                             
%FORMAT(10I8)                                                                   
       6
%FLAG MASS                                                                      
%FORMAT(5E16.8)                                                                 
  1.60000000E+01
%FLAG ATOM_TYPE_INDEX                                                           
%FORMAT(10I8)                                                                   
       1
%FLAG NUMBER_EXCLUDED_ATOMS                                                     
%FORMAT(10I8)                                                                   
       1
%FLAG NONBONDED_PARM_INDEX                                                      
%FORMAT(10I8)                                                                   
       1
%FLAG RESIDUE_LABEL                                                             
%FORMAT(20a4)                                                                   
MOL 
%FLAG RESIDUE_POINTER                                                           
%FORMAT(10I8)                                                                   
       1
%FLAG BOND_FORCE_CONSTANT                                                       
%FORMAT(5E16.8)                                                                 
  
%FLAG BOND_EQUIL_VALUE                                                          
%FORMAT(5E16.8)                                                                 
  
%FLAG ANGLE_FORCE_CONSTANT                                                      
%FORMAT(5E16.8)                                                                 
  
%FLAG ANGLE_EQUIL_VALUE                                                         
%FORMAT(5E16.8)                                                                 
  
%FLAG DIHEDRAL_FORCE_CONSTANT                                                   
%FORMAT(5E16.8)                                                                 

%FLAG DIHEDRAL_PERIODICITY                                                      
%FORMAT(5E16.8)                                                                 

%FLAG DIHEDRAL_PHASE                                                            
%FORMAT(5E16.8)                                                                 

%FLAG SCEE_SCALE_FACTOR                                                         
%FORMAT(5E16.8)                                                                 

%FLAG SCNB_SCALE_FACTOR                                                         
%FORMAT(5E16.8)                                                                 

%FLAG SOLTY                                                                     
%FORMAT(5E16.8)                                                                 
  0.00000000E+00
%FLAG LENNARD_JONES_ACOEF                                                       
%FORMAT(5E16.8)                                                                 
 {acoef: .8E}
%FLAG LENNARD_JONES_BCOEF                                                       
%FORMAT(5E16.8)                                                                 
 {bcoef: .8E}
%FLAG BONDS_INC_HYDROGEN                                                        
%FORMAT(10I8)                                                                   
     
%FLAG BONDS_WITHOUT_HYDROGEN                                                    
%FORMAT(10I8)                                                                   

%FLAG ANGLES_INC_HYDROGEN                                                       
%FORMAT(10I8)                                                                   
       
%FLAG ANGLES_WITHOUT_HYDROGEN                                                   
%FORMAT(10I8)                                                                   

%FLAG DIHEDRALS_INC_HYDROGEN                                                    
%FORMAT(10I8)                                                                   

%FLAG DIHEDRALS_WITHOUT_HYDROGEN                                                
%FORMAT(10I8)                                                                   

%FLAG EXCLUDED_ATOMS_LIST                                                       
%FORMAT(10I8)                                                                   
       1
%FLAG HBOND_ACOEF                                                               
%FORMAT(5E16.8)                                                                 

%FLAG HBOND_BCOEF                                                               
%FORMAT(5E16.8)                                                                 

%FLAG HBCUT                                                                     
%FORMAT(5E16.8)                                                                 

%FLAG AMBER_ATOM_TYPE                                                           
%FORMAT(20a4)                                                                   
c 
%FLAG TREE_CHAIN_CLASSIFICATION                                                 
%FORMAT(20a4)                                                                   
M   
%FLAG JOIN_ARRAY                                                                
%FORMAT(10I8)                                                                   
       0
%FLAG IROTAT                                                                    
%FORMAT(10I8)                                                                   
       0
%FLAG RADIUS_SET                                                                
%FORMAT(1a80)                                                                   
                                                   
%FLAG RADII                                                                     
%FORMAT(5E16.8)                                                                 
  0.00000000E+00  
%FLAG SCREEN                                                                    
%FORMAT(5E16.8)                                                                 
  0.00000000E
%FLAG IPOL                                                                      
%FORMAT(1I8)                                                                    
       0
"""

PDB = """ATOM      1 C    MOL     1       0.001   1.000   0.000  1.00  0.00
TER
END
"""

LJ_EPSILONS = [0.0078125,
               0.015625,
               0.03125,
               0.0625,
               0.125,
               0.25,
               0.5,
               1,
               2,
               4]

LJ_SIGMA = [0.6,
            1.4,
            2.2,
            3,
            3.8,
            4.6,
            5.4,
            6.2,
            7]


#Rows correspond to lj_epsilon and columns to lj_sigma
ENERGY_M = """0.417	0.907	1.66	2.71	4.06	5.69	7.58	9.73	12.09
0.412	0.905	1.68	2.73	4.1	5.75	7.59	9.71	12.15
0.386	0.878	1.61	2.68	4.01	5.63	7.42	9.57	11.78
0.332	0.799	1.51	2.49	3.79	5.26	7	8.89	10.85
0.244	0.632	1.29	2.17	3.26	4.57	6.08	7.65	9.47
0.088	0.385	0.88	1.59	2.47	3.48	4.63	5.73	6.98
-0.139	-0.023	0.25	0.67	1.15	1.72	2.16	2.84	3.15
-0.504	-0.64	-0.72	-0.75	-0.85	-0.96	-1.32	-1.82	-2.51
-1.047	-1.58	-2.18	-2.9	-3.84	-5.11	-6.62	-8.58	-11.06
-1.868	-3	-4.36	-6.1	-8.3	-11.07	-14.45	-18.62	-23.81"""


XVV_298_p = "/home/max/Documents/PhD/calculations/LJ_Spheres/lj_with_hnc_closure/water_298.15.xvv"

C_O_DX = 'c_{}.O.1.dx'
C_ONP_DX = 'c_{}.O.1.dx'
H_ONP_DX = 'h_{}.O.1.dx'
H_O_DX = 'h_{}.O.1.dx'
U_O_DX = 'u_{}.O.1.dx'
C_HNP_DX = 'c_{}.H1.1.dx'
H_HNP_DX = 'h_{}.H1.1.dx'

T = 298.15

CSV_HEADER_NAMES = ['LJ_epsilon', 'LJ_sigma', 'MD energy', 'HNC', 'GF',
                    'UC', 'PMV']
                    

def gamma_cor(integral, T, gamma=0.38):
    t = float(T)
    density = water.concentration(t)*6.0221413E-4  #molecules / A^3
    R = 1.9872041E-3 #kcal/K/Mol
    return R*t*density/2*integral*(1-gamma)


def uc(gf, pmv, t):
    """Return universal correction."""
    gf = float(gf)
    pmv = float(pmv)
    t = float(t)
    density = water.concentration(t)*6.0221413E-4  #molecules / A^3
    return gf - 3.2217*density*pmv + 0.5783  #kcal/mol
    
    
def two_dx_file_correction(integrable_dx_m, contstraint_dx_m,
                           spacing3, constraint_value, bigger=True):
    """
    Integrate integrable_dx matrix over the region where
    constraint_dx matrix >= constraint_value. (<= if bigger=False)
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
        if len(ln) <= 3:
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


def calculate_ng_corrections(out_name, T, zero_name):
    """Calculate different nonpolar corrections of correlation functions:
    """
    T = float(T)
    p, no_p_name = os.path.split(out_name)
    _, no_p_zero_name = os.path.split(zero_name)
    # Load dx files
    c_np_dx, _, spacing_l = load_dx(os.path.join(p, C_ONP_DX.format(no_p_zero_name)))
    h_onp_dx = load_dx(os.path.join(p, H_ONP_DX.format(no_p_zero_name)))[0]    
    c_hnp_dx = load_dx(os.path.join(p, C_HNP_DX.format(no_p_zero_name)))[0]
    h_hnp_dx = load_dx(os.path.join(p, H_HNP_DX.format(no_p_zero_name)))[0]    

    u_dx = load_dx(os.path.join(p, U_O_DX.format(no_p_name)))[0]    
    spacing3 = reduce(lambda x, y: x*y, spacing_l)
    kt_2_kcal = 4184/8.3144621/T  #kt/kcal/mol
    # hydrogen and oxygen corellation function integrals
    # v_in is defined as region where u > 10 kcal/mol
    # [(-c_v_in, -c*h_v_in, -c-0.5c*h_v_in), (same, but v_out)] for oxygen
    gf_o = [(two_dx_file_correction(c_np_dx, u_dx, spacing3, 10*kt_2_kcal, i),
             two_dx_file_correction(c_np_dx*h_onp_dx, u_dx, spacing3, 10*kt_2_kcal, i),
             two_dx_file_correction(-c_np_dx-0.5*h_onp_dx*c_np_dx, u_dx, spacing3, 10*kt_2_kcal, i))
             for i in [True, False]]
    # [(-c_v_in, -c*h_v_in, -c-0.5c*h_v_in), (same, but v_out)] for hydrogen
    gf_h = [(two_dx_file_correction(c_hnp_dx, u_dx, spacing3, 10*kt_2_kcal, i),
             two_dx_file_correction(c_hnp_dx*h_hnp_dx, u_dx, spacing3, 10*kt_2_kcal, i),
             two_dx_file_correction(-c_hnp_dx-0.5*c_hnp_dx*h_hnp_dx, u_dx, spacing3, 10*kt_2_kcal, i))
             for i in [True, False]]
    # NgB correction
    ng_u_cor = gf_o[0][0]
    return ng_u_cor, gf_o, gf_h


def gen_energy_matrix():
    """Convert matrix of MD computed solvation energies from text format to 
    numpy array matrix."""
    m = []
    energy_l = ENERGY_M.splitlines()
    for l in energy_l:
        m.append(map(float, l.split()))
    return np.array(m)


#def gen_csv_line(e, s, md_e, log_name, ng_u_cor, gf_o, gf_h):
#    csv_line = [e, s, md_e]
#    with open(log_name, 'rb') as f:
#        for line in f:
#            if line[0:11]=="rism_exchem":
#                kh=float(line.split()[1])
#            if line[0:11]=="rism_exchGF":
#                gf=float(line.split()[1])
#            if line[0:11]=="rism_volume":
#                pmv=float(line.split()[1])
#    uc_value = uc(gf, pmv, T)
#    ngb_value = gamma_cor(ng_u_cor, T) + kh
#    csv_line.extend([kh, gf, ngb_value, uc_value, pmv])
#    csv_line.extend([gf_o[0][0],
#                 gf_o[1][0],
#                 gf_o[0][1],
#                 gf_o[1][1],
#                 gf_o[0][2],
#                 gf_o[1][2],
#                 gf_h[0][0],
#                 gf_h[1][0],
#                 gf_h[0][1],
#                 gf_h[1][1],
#                 gf_h[0][2],
#                 gf_h[1][2]])
#    return csv_line

def gen_csv_line(e, s, md_e, log_name):
    csv_line = [e, s, md_e]
    with open(log_name, 'rb') as f:
        for line in f:
            if line[0:11]=="rism_exchem":
                kh=float(line.split()[1])
            if line[0:11]=="rism_exchGF":
                gf=float(line.split()[1])
            if line[0:11]=="rism_volume":
                pmv=float(line.split()[1])
    uc_value = uc(gf, pmv, T)
    csv_line.extend([kh, gf, uc_value, pmv])
    return csv_line
    

def main():
    dir_name = '.'
#    os.mkdir(dir_name)
#    energy_m = gen_energy_matrix()
#    csv_rows = [CSV_HEADER_NAMES]
    chg = 18.2223
    for k, e in enumerate(LJ_EPSILONS):
        for l, s in enumerate(LJ_SIGMA):
            p = os.path.join(dir_name, '{}_{}'.format(e, s))
            os.mkdir(p)
            pdb_name = os.path.join(p, '{}_{}.pdb'.format(e, s))
            name = pdb_name[:-4]
            with open(pdb_name, 'wb') as f:
                f.write(PDB)
            acoef = 4*e*s**12
            bcoef = 4*e*s**6
            prmtop_name = os.path.join(p, '{}_{}.prmtop'.format(e, s))
            with open(prmtop_name, 'wb') as f:
                f.write(PRMTOP.format(chg=chg,acoef=acoef, bcoef=bcoef))
#            rism_calc = rism_singlpnt.RISM3D_Singlpnt(name, T,
#                                                      prmtop_name, XVV_298_p)
#            rism_calc.setup_calclation(closure='hnc', calculate_u=False)
#            rism_calc.run_calculation_and_log()
#            log_name = rism_calc.log_name
            #ng_u_cor, gf_o, gf_h = calculate_ng_corrections(name, T, name)
 #           md_e = energy_m[k, l]
            #csv_rows.append(gen_csv_line(e, s, md_e, log_name, ng_u_cor, gf_o, gf_h))
#            csv_rows.append(gen_csv_line(e, s, md_e, log_name))
    
    #with open('lj_results.csv', 'wb') as f:
    #    wr = csv.writer(f)
    #    wr.writerows(csv_rows)


if __name__ == '__main__':
    main()






