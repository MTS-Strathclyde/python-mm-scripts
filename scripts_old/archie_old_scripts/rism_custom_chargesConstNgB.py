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

import rism3dConstNgB.pre_post_calc as pre_post_calc
import rism3dConstNgB.rism_singlpnt as rism_singlpnt
import rism3dConstNgB.results as results

REQUIRED_EXECUTABLES = ['antechamber', 'parmchk', 'tleap', 'rism3d.snglpnt']


def calculate(name, T, prmtop_name):
    rism_calc = rism_singlpnt.RISM3D_Singlpnt(name, T, prmtop_name)
    rism_calc.setup_calclation()
    rism_calc.run_calculation_and_log()
    out_name = rism_calc.out_name
#    zero_prmtop_name = charge.generate_prmtops(out_name, ['zero'])[0]
#    rism_zero_calc = rism_singlpnt.RISM3D_Singlpnt(name, T, zero_prmtop_name)
#    rism_zero_calc.setup_calclation()
#    rism_zero_calc.run_calculation_and_log()
    return out_name


def rism3d_custom_charge_models(mol_path, T=298.15, charge_model_l=None, 
                                delete_dx_files=True):
    """Run RISM calculation and delete dx files."""
    T = round(float(T), 2)
    charge_model_l.append('zero')
    if 253.15 <= T <= 383.15:
        try:
            name = pre_post_calc.prepare_calc_directory(mol_path, T)
            #kh, gf, pmv, ucorr = run_3drism(name, T)
        except OSError:
            print 'Skipping job {} T = {}'.format(mol_path, T)
            return 1
        
        prmtop_names = \
                    pre_post_calc.prepare_3drism_calc(name, charge_model_l, T)

        out_names = []
        for prmtop_name in prmtop_names:
            out_names.append(calculate(name, T, prmtop_name))
        for out_name in out_names[:-1]:
            results.write_results(out_name, T)
        if delete_dx_files:
            pre_post_calc.delete_dx(name)
    else:
        print 'Temperature outside of scripts range!'
    

def main(argv):
    rism3d_custom_charge_models(argv[0], argv[1], argv[2:])


if __name__ == '__main__':
    main(sys.argv[1:])



    
