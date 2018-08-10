#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 13:06:13 2015

@author: max
"""

import numpy as np
import sys
import argparse
import matplotlib.pyplot as plt




def process_command_line(argv):
    """Processes arguments

    Parameters
    ----------
    argv : list
        Command line arguments.

    Returns
    -------
    out : argparse.Namespace
        Namespace of command line arguments.
    """
    parser = argparse.ArgumentParser(description="""Vizualize volumetric distributions.""")
    parser.add_argument('dx',
                        help=""" DX file.""")
    parser.add_argument('dx2',
                        help=""" DX file.""")
#    parser.add_argument('idx',
#                        help=""" idx.""", type=int)
    return parser.parse_args(argv)



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
        x = np.arange(OrX, OrX + dX*size[0], dX)
        y = np.arange(OrY, OrY + dY*size[1], dY)
        z = np.arange(OrZ, OrZ + dZ*size[2], dZ)
        x_mat, y_mat, z_mat = np.meshgrid(x, y, z, indexing='ij')
        return dx_m.astype('float'), x_mat, y_mat, z_mat
    else:
        return dx_m.astype('float'), (OrX, OrY, OrZ), (dX, dY, dZ)    


def potential(dx,rho_q):
    vacuum_p = 8.854187817e-12 #F/m
    conver = 3.694707135311445e-8 # m*e/(F*A) to kcal/mol/e
    coef = 1./(vacuum_p)*conver
    x = np.arange(0,dx*len(rho_q) - dx/2.,dx)
    phi_x = []
    # density times direction towards vapour
    # $$\phi (z) = -\frac{1}{\epsilon_{0}} \int \limits_{- \infty}^{z} \rho(z') (z - z') dz'$$
    for i, x_c in enumerate(x):
        phi_x.append(-np.sum((x_c-x[:i+1])*rho_q[:i+1])*dx*coef )
    return x,phi_x


def main(argv):
    args = process_command_line(argv)
    dx_m, _, (dx,dy,dz) = load_dx(args.dx)
    dx_m2, _, (dx,dy,dz) = load_dx(args.dx2)
    #idx = args.idx
    idx = 0
    go  = np.array([np.mean(dx_m[i,:,:]) for i, _ in enumerate(dx_m[:,idx,idx])])
    gh  = np.array([np.mean(dx_m2[i,:,:]) for i, _ in enumerate(dx_m2[:,idx,idx])])
    rho_q = 3.3328311138810005E-02 * (go*-0.8476 + 2*gh*0.4238)
    rho_q = rho_q - np.mean(rho_q)
    print np.sum(rho_q)
    #rho_q[40] = rho_q[40] - np.sum(rho_q)
    x,phi_x = potential(dx,rho_q)
    print phi_x[len(phi_x)/2]
    mat = np.c_[x,go,gh]
    np.savetxt('gr.txt',mat)
    plt.plot(x,rho_q)
    plt.figure()
    plt.plot(x,go)
    plt.plot(x,gh)
    plt.figure()
    plt.plot(x,phi_x)
    plt.show()
    
    
#    plt.matshow(dx_m[:,:,args.idx])
#    plt.figure()
#    plt.plot(dx_m[:, args.idx, args.idx])
#    plt.show()
    
    
if __name__=='__main__':
    main(sys.argv[1:])
    
    
