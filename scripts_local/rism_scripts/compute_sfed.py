#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 13:27:35 2015

@author: max
"""

import numpy as np
import sys


def load_dx(fname):
    """Return dx matrix, origin coordinates tuple and spacing between points tuple."""
    with open(fname, 'rb') as f:
        txt = f.readlines()
    if txt[0].startswith('object 1 class gridpositions counts'):
        size = txt[0][35:].split()
        size = map(int, size)
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


def hnc_fed(g, c):
    """ hnc free energy, g - is pair correlation function
    and c - is direct correlation function """
    h = g - 1
    return .5*(h**2 - c*h) - c


def main():
    kT = 1.9872/1000*298.15  # kcal/mol
    rho = 3.3328e-02         # 1/A^3
    # assume arguments are  g_O.dx, g_H.dx, c_O.dx, c_H.dx
    grids = []
    for f in sys.argv[1:]:
        grid, origin, (dx,dy,dz) = load_dx(f)
        grids.append(grid)
    o_energy = rho*kT*hnc_fed(grids[0], grids[2])    # kcal/mol/A^3
    h_energy = 2*rho*kT*hnc_fed(grids[1], grids[3])  # kcal/mol/A^3
    total = o_energy + h_energy
    # save total
    np.savetxt(total, 'sfed.txt')
    # if we want to load sfed - simply use:
    #sfed = np.loadtxt('sfed.txt')
    
    # check results
    total_fe = np.sum(total)*dx*dy*dz
    print 'Total FE:',total_fe
    
                
if __name__ == '__main__':
    main()

