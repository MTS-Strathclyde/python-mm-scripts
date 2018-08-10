#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 11 14:19:15 2014

@author: max
"""

import os
import sys
import glob



def collect_result_names(base_dir):
    res_names = []
    for p, dnames, fnames in os.walk(base_dir):
        for fn in fnames:
            if fn.startswith('results') and fn.endswith('.txt'):
                res_names.append((p, fn))
    return res_names

def collect_dbf_names(base_dir):
    res_names = []
    for p, dnames, fnames in os.walk(base_dir):
    
        for fn in fnames:
            if fn.endswith('.dbf'):
                res_names.append((p, fn))
    return res_names



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


def density(t): return water_concentration(t)*6.0221413E-4


def mdft_cor(kh, pmv, temp, compres):
    """Return fixed and free solute correction."""
    w_density = density(temp)
    brackets = 1 - 1e-20/(w_density*1.3806488E-23*temp*compres)
    cor = 1.9872041E-3*temp/2*brackets*w_density*pmv
    fix_cor = -w_density*1.9872041E-3*temp*pmv
    return kh + cor + fix_cor, kh + cor
    

def parse_res_names(res_names):
    """Iterates over results and computes mdft correction."""
    for p, resname in res_names:
        res_path = os.path.join(p, resname)
        with open(res_path, 'rb') as f:
            lines = f.readlines()
        try:
            kh = float(lines[0].split()[1])
            pmv = float(lines[2].split()[1])
        except IndexError, e:
            print res_path
            raise e
        temp = float(p.split('/')[-1])
        water_therm_name = 'water_{}.therm'.format(temp)
        water_therm_p = os.path.join(p, water_therm_name)
        with open(water_therm_p, 'rb') as f:
            w_lines = f.readlines()
        compres = float(w_lines[2].split()[-1])
        fix_kh_cor, kh_cor = mdft_cor(kh, pmv, temp, compres)
        mdft_cor_fname = os.path.join(p, 'mdft_cor.txt')
        with open(mdft_cor_fname, 'wb') as f:
            f.write('kh_fix {}\n'.format(fix_kh_cor))
            f.write('kh {}\n'.format(kh_cor))
        try:
            with open(os.path.join(p, 'results.txt'), 'ab') as f:
                f.write('kh_fix {}\n'.format(fix_kh_cor))
                f.write('kh {}\n'.format(kh_cor))
        except IOError:
            print 'results.txt wasn\'t found in {}'.format(res_path)
        try:
            dbf = glob.glob(os.path.join(p, '*.dbf'))[0]
            with open(dbf, 'ab') as f:
                f.write('pmvc_fix_kh, {}, -, -\n'.format(fix_kh_cor))
                f.write('pmvc_kh, {}, -, -\n'.format(kh_cor))
        except IndexError:
            print 'Warning: *.dbf wasn\'t found in {}'.format(res_path)  
            print 'You can safely ignore this warning if you don\'t need to use *.dbf file.'

def parse_dbf_names(dbf_names):
    """Iterates over results and computes mdft correction."""
    for p, dbf_name in dbf_names:
        dbf_path = os.path.join(p, dbf_name)
        try:
            with open(dbf_path, 'rb') as f:
                lines = f.readlines()
            for line in lines:
                if line.startswith('rism_exchem,'):
                    kh = float(line.split(',')[1])
                if line.startswith('rism_volume,'):
                    pmv = float(line.split(',')[1])                
        except IndexError, e:
            print dbf_path
            raise e
        temp = float(p.split('/')[-1])
        water_therm_name = 'water_{}.therm'.format(temp)
        water_therm_p = os.path.join(p, water_therm_name)
        with open(water_therm_p, 'rb') as f:
            w_lines = f.readlines()
        compres = float(w_lines[2].split()[-1])
        fix_kh_cor, kh_cor = mdft_cor(kh, pmv, temp, compres)
        mdft_cor_fname = os.path.join(p, 'mdft_cor.txt')
        with open(mdft_cor_fname, 'wb') as f:
            f.write('kh_fix {}\n'.format(fix_kh_cor))
            f.write('kh {}\n'.format(kh_cor))
        try:
            with open(os.path.join(p, 'results.txt'), 'ab') as f:
                f.write('kh_fix {}\n'.format(fix_kh_cor))
                f.write('kh {}\n'.format(kh_cor))
        except IOError:
            print 'results.txt wasn\'t found in {}'.format(p)
        try:
            dbf = glob.glob(os.path.join(p, '*.dbf'))[0]
            with open(dbf, 'ab') as f:
                f.write('pmvc_fix_kh, {}, -, -\n'.format(fix_kh_cor))
                f.write('pmvc_kh, {}, -, -\n'.format(kh_cor))
        except IndexError:
            print 'Warning: *.dbf wasn\'t found in {}'.format(p)  
            print 'You can safely ignore this warning if you don\'t need to use *.dbf file.'


def main(argv):
    
    if argv[1] == 'dbf':
        parse_dbf_names(collect_dbf_names(argv[0]))
    else:
        parse_res_names(collect_result_names(argv[0]))


if __name__ == '__main__':
    main(sys.argv[1:])
        
        
