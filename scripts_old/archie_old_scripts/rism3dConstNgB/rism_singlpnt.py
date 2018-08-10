# -*- coding: utf-8 -*-
"""
Created on Fri May 23 09:43:52 2014

@author: max
"""

import os
import time
import datetime
import subprocess


class RISM3D_Singlpnt(object):
    RISM3D_GRID = '0.300000,0.300000,0.300000'
    RISM3D_TOLLERANCE = '0.0000000001'        
    def __init__(self, name, T, prmtop_name=None, xvv_name=None):
        self.name = name
        self.T = T
        self.p, self.no_p_name = os.path.split(name)
        self.out_name = prmtop_name[:-7]
        self.no_p_out_name = os.path.split(self.out_name)[1]
        self.log_name = '{}.log'.format(self.out_name)
        self.logfile = open(self.log_name, 'wb')
        # Wirte timestamp and start measuring runtime
        self.start_time = time.time()
        self.logfile.write(str(datetime.datetime.now()))
        self.logfile.flush()
        if prmtop_name:
            self.prmtop_name = self.no_p_out_name + '.prmtop'
        else:
            self.prmtop_name = '{}.prmtop'.format(self.no_p_name)
        if xvv_name:
            self.xvv_name = xvv_name
        else:
            self.xvv_name = 'water_{temp}.xvv'.format(temp=self.T)
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
                                         'h_{}'.format(self.no_p_out_name)])
        if calculate_c:
            self.run_flags_list.extend(['--cuv', 
                                         'c_{}'.format(self.no_p_out_name)])
        if calculate_u:
            self.run_flags_list.extend(['--uuv', 
                                         'u_{}'.format(self.no_p_out_name)])
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
