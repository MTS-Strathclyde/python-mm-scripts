#!/usr/bin/env python

# -*- coding: utf-8 -*-
"""
Created on Wed May 28 12:49:25 2014

@author: max
"""

import os
import numpy as np
import csv
import sys

from rism3d import water as water


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



C_O_DX = 'c_{}.O.1.dx'
H_O_DX = 'h_{}.O.1.dx'
U_O_DX = 'u_{}.O.1.dx'
U_H_DX = 'u_{}.H1.1.dx'
C_H_DX = 'c_{}.H1.1.dx'
H_H_DX = 'h_{}.H1.1.dx'

T = 298.15


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
    # Load dx files
#    c_o_dx, _, spacing_l = load_dx(os.path.join(p, C_O_DX.format(no_p_name)))
   # h_o_dx = load_dx(os.path.join(p, H_O_DX.format(no_p_name)))[0]    
    c_h_dx,_,spacing_l = load_dx(os.path.join(p, C_H_DX.format(no_p_name)))
    
   # h_h_dx = load_dx(os.path.join(p, H_H_DX.format(no_p_name)))[0]    

#    u_o_dx = load_dx(os.path.join(p, U_O_DX.format(no_p_name)))[0]
    u_h_dx = load_dx(os.path.join(p, U_H_DX.format(no_p_name)))[0]    
    spacing3 = reduce(lambda x, y: x*y, spacing_l)
    #kt_2_kcal = 4184/8.3144621/T  #kt/kcal/mol
    
    # the structure is following:
    # V-in
    # kt = -3, -2, -1, 0, 1, 2, 3
    
    
    
#    gf_o = [two_dx_file_correction(c_o_dx, u_dx, spacing3, j)
#             for j in range(-3, 4)]
    gf_h = [two_dx_file_correction(c_h_dx, u_h_dx, spacing3, j)
             for j in range(0, 1)]
#    gf_np_o = [two_dx_file_correction(c_np_o_dx, u_dx, spacing3, j)
#             for j in range(-3, 4)]
#    gf_np_h = [two_dx_file_correction(c_np_h_dx, u_dx, spacing3, j)
#             for j in range(-3, 4)]                 
    return gf_h
#    return  gf_o, gf_h, gf_np_o, gf_np_h




CSV_HEADER_NAMES = ['LJ_epsilon', 'LJ_sigma', 'KH', 'GF',
                    'UC', 'PMV', 
                    'h_dcfi_v_in_0kt',]
                    


def gen_csv_line(e, s, log_name, gf_h):
    csv_line = [e, s]
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
    csv_line.extend(gf_h)
    return csv_line
    

def main(args):
    dir_name = args[0]
    csv_rows = [CSV_HEADER_NAMES]
    
    for k, e in enumerate(LJ_EPSILONS):
        for l, s in enumerate(LJ_SIGMA):
            p = os.path.join(dir_name, '{}_{}'.format(e, s))

            pdb_name = os.path.join(p, '{}_{}.pdb'.format(e, s))
            name = pdb_name[:-4]

            
            gf_h = calculate_ng_corrections(name, T, name)
            csv_rows.append(gen_csv_line(e, s, name + '.log', gf_h))
            print name, "analyzed"
    
    with open('lj_results_h_vin.csv', 'wb') as f:
        wr = csv.writer(f)
        wr.writerows(csv_rows)


if __name__ == '__main__':
    main(sys.argv[1:])






