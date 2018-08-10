#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 31 15:22:43 2016

@author: max
"""

import sys

with open(sys.argv[1]) as f:
    txt = f.readlines()

data = []    
for l in txt:
    if l.startswith('rism_exchem'):
        data.append(float(l.split()[1]))
    if l.startswith('rism_exchGF'):
        data.append(float(l.split()[1]))
    if l.startswith('rism_volume'):
        data.append(float(l.split()[1]))

print ','.join(map(str,data)) + ','
    