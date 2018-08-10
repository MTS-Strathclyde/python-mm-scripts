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


def main():
    # assume arguments are  g_O.dx, g_H.dx, c_O.dx, c_H.dx, acr.dx
    grids = []
    for f in sys.argv[1:]:
        grids.append(load_dx(f)[0])
    # compute short contributions and total value
    c_short_o = grids[2] - -2.0066044493929137E+01*grids[4]
    c_short_h = grids[3] -  1.0033022246964569E+01*grids[4]
    o_contrib_short = 1.9872041/1000*298.15*3.3328311138810005E-02*\
                np.sum((.5*(grids[0]-1)**2 - c_short_o- \
                0.5*(grids[0]-1)*c_short_o)*0.5*0.5*0.5)
    h_contrib_short = 2*1.9872041/1000*298.15*3.3328311138810005E-02*\
                np.sum((.5*(grids[1]-1)**2 - c_short_h- \
                0.5*(grids[1]-1)*c_short_h)*0.5*0.5*0.5)
    total = o_contrib_short + h_contrib_short
    
    print 'total, {}, o_short, {}, h_short, {}'.format(total, o_contrib_short,
                                                        h_contrib_short)
    
                
if __name__ == '__main__':
    main()

