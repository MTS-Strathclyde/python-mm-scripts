#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri May  6 11:01:35 2016

@author: max
"""

import sys

if __name__ == '__main__':
    quadrupole = 0
    traceless_quadrupole = 0
    with open(sys.argv[1]) as f:
        for l in f:
            if l.startswith(' Dipole moment '):
                mom = f.next().split()
                mom = float(mom[-1])
            if l.startswith(' Quadrupole moment '):
                trace = f.next().split()
                trace = float(trace[1])**2 + float(trace[3])**2 + float(trace[5])**2
                offdiag = f.next().split()
                offdiag = float(offdiag[1])**2 + float(offdiag[3])**2 + float(offdiag[5])**2
                quadrupole = trace + 2*offdiag
            if l.startswith(' Traceless Quadrupole moment'):
                trace = f.next().split()
                trace = float(trace[1])**2 + float(trace[3])**2 + float(trace[5])**2
                offdiag = f.next().split()
                offdiag = float(offdiag[1])**2 + float(offdiag[3])**2 + float(offdiag[5])**2
                traceless_quadrupole = trace + 2*offdiag
    print mom, quadrupole, traceless_quadrupole
                