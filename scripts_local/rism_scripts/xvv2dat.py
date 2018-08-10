#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Created on Mon Dec 7 15:17:42 2015

@author: max
"""

import sys
import numpy as np



def read_xvv(fname):
    data = []
    with open(fname) as f:
        start_xvv = False
        for l in f:
            if l.startswith('%FLAG POINTERS'):
                f.next()
                nr, nsites, _ = f.next().split()
                nr = int(nr)
                nsites = int(nsites)
            if l.startswith('%FLAG THERMO'):
                f.next()
                dr = f.next().split()[-1]
                dr = float(dr)
            if l.startswith('%FLAG XVV'):
                start_xvv = True
                f.next()
                continue
            if start_xvv:
                data.extend(l.split())
            if l.startswith('%FLAG XVV_DT'):
                break
    dk = np.pi/(nr*dr)
    xvv = []
    for i in range(nsites**2):
        xvv.append([])
        for j in range(nr):
            xvv[-1].append(data[j+i*nr])
    k = [i*dk for i in range(nr)]
    xvv = np.array(xvv, dtype=float)
    return k, xvv
                    

if __name__  == '__main__':
    k, xvv = read_xvv(sys.argv[1])
    combined = np.c_[k, xvv.T]
    np.savetxt(sys.argv[1] + '.xvg', combined)

