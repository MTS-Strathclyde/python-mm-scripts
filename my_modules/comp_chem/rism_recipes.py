# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 13:29:06 2015

@author: max
"""
import numpy as np

def load_dx(fname, return_xyz=False):
    """Return dx matrix, origin coordinates tuple and spacing between points tuple.
    If return_xyz is true returns dx, x, y, and z matrices."""
    with open(fname, 'rb') as f:
        txt = f.readlines()
    if txt[0].startswith('object 1 class gridpositions counts'):
        size = txt[0][35:].split()
        size = map(int, size)
    else:
        print('Wrong file format')
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
    if return_xyz:
        x = np.arange(OrX, OrX + dX*size[0] - dX/2., dX)
        y = np.arange(OrY, OrY + dY*size[1] - dY/2., dY)
        z = np.arange(OrZ, OrZ + dZ*size[2] - dZ/2., dZ)
        x_mat, y_mat, z_mat = np.meshgrid(x, y, z, indexing='ij')
        return dx_m.astype('float'), x_mat, y_mat, z_mat
    else:
        return dx_m.astype('float'), (OrX, OrY, OrZ), (dX, dY, dZ)    
