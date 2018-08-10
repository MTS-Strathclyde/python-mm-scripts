#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 28 19:08:51 2016

@author: max
"""
from scipy import integrate
from scipy import ndimage
import numpy as np
import scipy as sp
import numpy.ma as ma
import sys
import math
from scipy import interpolate as interp
import argparse
import warnings

warnings.filterwarnings('error')

T=298.15
kB = 1.9872041E-3 # boltzmann const in kcal/mol/K


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
    parser = argparse.ArgumentParser(description=""" Compute short and long contributions.""")
    parser.add_argument('z_r',
                        help=""" z matrix.""")
    parser.add_argument('-m', '--multiplicities',
                        help="""Solvent sites multiplicities.""", nargs='+',
                        type=int)
    parser.add_argument('-g', '--guv',
                        help=""" Spatial distribution functions.""", nargs='+')
    parser.add_argument('-t', '--huv',
                        help=""" Total distribution functions.""", nargs='+')
    parser.add_argument('-u', '--uuv',
                        help=""" Potential energy dist.""", nargs='+')
#    parser.add_argument('-z', '--huv',
#                        help=""" Total correlation functions.""", nargs='+')
    parser.add_argument('-d', '--density',
                        help=""" Solvent number density. [0.03332]""",
                         default=0.03332, type=float)
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
        x = np.arange(OrX, OrX + dX*size[0] - dX/2., dX)
        y = np.arange(OrY, OrY + dY*size[1] - dY/2., dY)
        z = np.arange(OrZ, OrZ + dZ*size[2] - dZ/2., dZ)
        x_mat, y_mat, z_mat = np.meshgrid(x, y, z, indexing='ij')
        return dx_m.astype('float'), x_mat, y_mat, z_mat
    else:
        return dx_m.astype('float'), (OrX, OrY, OrZ), (dX, dY, dZ)    
        


#def interp_k(zk):
#    """ Interpolate for the grid used in 3d-rism """
#    

def F_ideal_int(g_s, density, multiplicities):
    rho_s = [density*g_i for g_i in g_s]
    rho_s_ma = [ma.masked_values(rho_i/density, 0) for rho_i in rho_s]
    log_s = [ma.log(rho_i_ma) for rho_i_ma in rho_s_ma]
    integral = sum([m*(rho_i * ma.filled(log_i, 0) - (rho_i - density)) \
                for m, rho_i, log_i in zip(multiplicities, rho_s, log_s)])
    return np.sum(integral)
    
    
def F_ext_int(g_s,  u_s,  density, multiplicities):
    integrand = sum([m*g_i*u_i*density for m, g_i, u_i in \
                     zip(multiplicities, g_s, u_s)])
    return np.sum(integrand)
    

def slow_convolution(delta_rho, zr, r, x, y, z, voxel):
    zr_interp = interp.interp1d(r, zr, 'nearest')
    integral = np.zeros_like(delta_rho)
    count = 0
    # main iterator
    it = np.nditer([x, y, z], flags=['multi_index'])
    # convolution iterator
    convol_it = np.nditer([delta_rho, x,y,z])
    # compute convolution
    while not it.finished:
        xi, yi, zi = it[0], it[1], it[2]
        convol = 0
        while not convol_it.finished:
            xj, yj, zj = convol_it[1:]
            r = math.sqrt((xi - xj)**2 + (yi - yj)**2 + (zi - zj)**2)
            zr_ij = zr_interp(r)
            convol += zr_ij*convol_it[0]
            count += 1
            if count % 10000 == 0:
                print(count)
            convol_it.iternext()
        #end convo
        integral[it.multi_index] = convol*voxel
        convol_it.reset()
        it.iternext()
    #end main
    return integral


def F_exc_int_new(h1, h2, zr, r, x, y, z, voxel):
    """ updated version"""
    convol = slow_convolution(h2, zr, r, x, y, z, voxel)
    return np.sum(convol*h1)


def F_exc_int(h1, u1, h2, u2, zr, r, x, y, z, voxel):
    """ original ver"""
    dr = r[1] - r[0]
    integral = 0
    count = 0
    mask1 = np.ones(u1.shape, dtype=bool)
    mask2 = np.ones(u2.shape, dtype=bool)
    print 'Total calcs:',np.sum(mask1)*np.sum(mask2)
    for h1i, xi, yi, zi in np.nditer([h1[mask1],x[mask1],y[mask1],z[mask1]]):
        convol = 0
#        if h1i < 0.:
#            print h1i
        for h2j, xj, yj, zj in np.nditer([h2[mask2],x[mask2],y[mask2],z[mask2]]):
            r = math.sqrt((xi - xj)**2 + (yi - yj)**2 + (zi - zj)**2)
            zr_ij = zr[int(round(r/dr))] # nearest interpolation
            convol += zr_ij*h2j
            count += 1
            if count % 1000000 == 0:
                print count
        #end for
        #print convol, xi, yi, zi
        try:
            integral += h1i*convol*voxel
        except Warning:
            print h1i
            print convol
            print voxel
            integral += h1i*convol*voxel
    #end for
    return integral


def F_exc_int_3d(h1, h2, z_3d, r, x, y, z, voxel):
    """ z3d  version"""
    #zr_interp = interp.interp1d(r, zr, 'nearest')
    #count = 0
    convol = ndimage.filters.convolve(h2, z_3d,mode='constant')
    #end for
    #print np.sum(convol)
    return np.sum(h1*convol)*voxel



def load_files(args):
    g_s = []
    u_s = []
    x = None
    if args.huv:
        args.guv = args.huv
    for g_name, u_name in zip(args.guv, args.uuv):
        if x is None:
            g_i, x, y, z = load_dx(g_name, True)
        else:
            g_i, _, (dx, dy, dz) = load_dx(g_name)
        u_i, _, (dx, dy, dz) = load_dx(u_name)
#        g_i[:,:,:] = 1.
#        u_i[:,:,:] = 0.
        if args.huv:
            g_i += 1
        g_s.append(g_i)
        u_s.append(u_i)
    voxel = dx*dy*dz
    

    print g_s[0].shape
    return g_s, u_s, x, y, z, voxel


def make_z_3d(r, z_r, x, y, z):
    zr_interp = interp.interp1d(r, z_r)
    z_r_3d = np.zeros_like(x)
    for z_r_at_r, xi, yi, zi in np.nditer([z_r_3d, x, y, z], 
           op_flags=[['readwrite'], ['readonly'], ['readonly'], ['readonly']]):
        r = math.sqrt(xi**2+yi**2+zi**2)
        z_r_at_r[...]  = zr_interp(r)
    return z_r_3d


def main(argv):
    args = process_command_line(argv)
    g_s, u_s, x, y, z, voxel = load_files(args)
    z_r = np.loadtxt(args.z_r, dtype=float)
    density = args.density
    multiplicities = np.array(args.multiplicities)
    #print np.sum(g_s[1]-g_s[0])
    r = z_r[:,0]
    F_id = F_ideal_int(g_s, density, multiplicities)*kB*T*voxel
    F_ext = F_ext_int(g_s, u_s, density, multiplicities)*voxel*kB*T
    F_double = 0
    interaction_pair = 1
    for i in range(len(multiplicities)):
        for j in range(i+1):
            print i,j
#            z_3d = make_z_3d(r, z_r[:, interaction_pair], x, y, z)
            dbl_integral = F_exc_int(g_s[i] - 1., u_s[i], g_s[j] - 1., u_s[i],
                                     z_r[:,interaction_pair], r, x, y, z, voxel)
#            dbl_integral = F_exc_int_3d(g_s[i] - 1., g_s[j] - 1.,
#                                        z_3d, r, x, y, z, voxel)
            if i != j:
                dbl_integral *= 2
            dbl_integral = dbl_integral*multiplicities[j]*multiplicities[j]*density**2
            print dbl_integral
            F_double += dbl_integral
            interaction_pair += 1
    F_exc = -F_double*voxel*kB*T/2.
    print ''
    print 'F_ideal',F_id
    print 'F_ext',F_ext
    print 'F_exc',F_exc
    print 'Total',F_id + F_ext + F_exc
    
        
        
if __name__ == '__main__':
    main(sys.argv[1:])
