# -*- coding: utf-8 -*-
"""
Created on Sat Oct  6 16:42:43 2012

@author: mishin1991

DECIDEDE TO ABONDON FOR NOW



"""

import pybel
import os
import subprocess
import sys


def convert_to_tinker(pymol):
    tmp_name_in = pymol.title + '_tmp_to_sdf2tinker'
    tmp_name_out = pymol.title + '_tmp_from_sdf2tinker'
    pymol.write('sdf', tmp_name)
    subprocess.call(['sdf2tinkerxyz', '<', tmp_name_in])
    pass

    
    
def tinker_scan()
    pass