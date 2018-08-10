# -*- coding: utf-8 -*-
"""
Created on Fri May 23 09:42:01 2014

@author: max
"""

import numpy as np
import os

import water

RESULTS_TXT = """dGhyd(KH)= {kh} kcal/mol
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

h_np_dcfi_v_in={1[0][0]}
h_np_dcfi_v_out={1[1][0]}
h_np_dcfi_times_h_v_in={1[0][1]}
h_np_dcfi_times_h_v_out={1[1][1]}
h_np_gf_v_in={1[0][2]}
h_np_gf_v_out={1[1][2]}


"""

C_O_DX = 'c_{}.O.1.dx'
C_ONP_DX = 'c_{}.O.1.dx'
H_ONP_DX = 'h_{}.O.1.dx'
H_O_DX = 'h_{}.O.1.dx'
U_O_DX = 'u_{}.O.1.dx'
C_HNP_DX = 'c_{}.H1.1.dx'
H_HNP_DX = 'h_{}.H1.1.dx'



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
    

def write_results(out_name, T):
    # TODO 
    # More robust zero name convention   
    zero_name = out_name.split('_')[0] + '_zero'
    p, no_p_name = os.path.split(out_name)
    log_name = '{}.log'.format(out_name)
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
    ng_u_cor, gf_o, gf_h = calculate_ng_corrections(out_name, T, zero_name)
    ng_u_cor_with_kh = kh + gamma_cor(ng_u_cor, T)
    
    results_txt = RESULTS_TXT.format(gf_o, gf_h, ng_u_cor_with_kh=ng_u_cor_with_kh,
                                     kh=kh, gf=gf, pmv=pmv, ucorr=ucorr)
    #Write results
    with open(os.path.join(p, 'results_{}.txt'.format(no_p_name)), 'wb') as f:
        f.write(results_txt)


    
    