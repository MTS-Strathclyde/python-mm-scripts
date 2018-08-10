#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 28 19:08:51 2016

@author: max
"""
import os
import numpy as np
import numpy.ma as ma
import glob
import sys
import math
from scipy import interpolate as interp
import argparse
import warnings

# my imports
from rism3d_pressure import Xvv

warnings.filterwarnings('error')

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
    parser = argparse.ArgumentParser(description=""" Compute short and long contributions.
    Optimizied only for a single solvent""")
    parser.add_argument('-x', '--xvv',
                        help=""" Susceptibility functions file. If present,
                        densities, z_k matrix and multiplicities will be ignored.""")
    parser.add_argument('-n', '--dir_name',
                        help=""" Calculation directory with dx files. If present,
                        huv and uuv arguments will be ignored.""")
    parser.add_argument('-t', '--huv',
                        help=""" Total distribution functions. Guv will be ignored.""", nargs='+')
    parser.add_argument('-u', '--uuv',
                        help=""" Potential energy dist.""", nargs='+')
    parser.add_argument('--z_k',
                        help=""" z matrix.""")
    parser.add_argument('-m', '--multiplicities',
                        help="""Solvent sites multiplicities.""", nargs='+',
                        type=int)
    parser.add_argument('-d', '--density',
                        help=""" Solvent number density. [0.03332]""",
                         default=3.3328311138810005E-02, type=float)
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
        

def F_ideal_int(g_s, density, multiplicities):
    # we evaluate \rho(i) * ln(\rho(i)/\rho_0) - \delta \rho(i)
    # this can be rewritten as:
    #    \ln(\rho(i)**\rho(i)) - \rho(i) * (\ln(rho_0)  + 1) + \rho_0
    #
    # old method: start
    # rho_s = [density*g_i for g_i in g_s]
    # rho_s_ma = [ma.masked_values(rho_i/density, 0,atol=1.0e-15) for rho_i in rho_s]
    # log_s = [ma.log(rho_i_ma) for rho_i_ma in rho_s_ma]
    # integral = sum([m*(rho_i * ma.filled(log_i, 0) - (rho_i - density)) \
    #            for m, rho_i, log_i in zip(multiplicities, rho_s, log_s)])
    # return np.sum(integral)
    # old method: stop
    #
    log_s = [density*np.log(g_i**g_i) for g_i in g_s]
    lin_s = [-density*(g_i - 1) for g_i in g_s]
    integral = sum([m*(log_i + lin_i) for m, log_i, lin_i in \
                   zip(multiplicities, log_s, lin_s)])
    return np.sum(integral)
    
    
def F_ext_int(g_s,  u_s,  density, multiplicities):
    integrand = sum([m*g_i*u_i*density for m, g_i, u_i in \
                     zip(multiplicities, g_s, u_s)])
    return np.sum(integrand)


def F_exc_int_3d(h1, u1, h2, z_3d, voxel):
    """ zk_3d  version"""
    h2_k = np.fft.fftn(h2)
    h2_k = np.fft.fftshift(h2_k)  # shift 0 freq to the middle
    convol = np.fft.ifftn(np.fft.ifftshift(h2_k*z_3d))
    return np.sum((h1*np.real(convol))[u1<5])
    

def load_files(args):
    g_s = []
    u_s = []
    x = None
    if args.dir_name:
        h_name = os.path.join(args.dir_name, 'h_*.dx')
        u_name = os.path.join(args.dir_name, 'u_*.dx')
        h_dxs = glob.glob(h_name)
        u_dxs = glob.glob(u_name)
        print 'Tot functions:',h_dxs
        print 'U functions:',u_dxs
        assert len(h_dxs) == len(u_dxs)
        assert len(h_dxs) > 0
    else:
        h_dxs = args.huv
        u_dxs = args.uuv
    #print h_dxs
    #print u_dxs
    for h_name, u_name in zip(h_dxs, u_dxs):
        if x is None:
            h_i, x, y, z = load_dx(h_name, True)
        else:
            h_i, _, (dx, dy, dz) = load_dx(h_name)
        u_i, _, (dx, dy, dz) = load_dx(u_name)
        h_i += 1   # make it pair correlation function
        g_s.append(h_i)
        u_s.append(u_i)
    voxel = dx*dy*dz
    #print 'Array size:',g_s[0].shape
    return g_s, u_s, x, y, z, voxel
    

def get_bulk_solvent_parameters(args):
    if args.xvv:
        xvv_inst = Xvv(args.xvv)
        k, zk_all = xvv_inst.compute_zk()
        z_k = k.T
        for m in range(xvv_inst.nsites):
            for n in range(m+1):
                z_k = np.c_[z_k, zk_all[:,m,n].T]
        T = xvv_inst.temperature
        density = xvv_inst.normalized_densities[0]
        multiplicities = xvv_inst.multiplicities
    else:
        T = 298.15
        z_k = np.loadtxt(args.z_k, dtype=float)
        density = args.density
        multiplicities = np.array(args.multiplicities)
        k = z_k[:,0]        
    return density, multiplicities, T, k, z_k
        

def make_z_3d(k, z_k, x, y, z):
    """ Makes 3d cube out of z, centered in the origin. """
    zk_interp = interp.interp1d(k, z_k)
    z_k_3d = np.zeros_like(x)
    dx = x[1,0,0] - x[0,0,0]
    nx = x.shape[0]
    dy = y[0,1,0] - y[0,0,0]
    ny = y.shape[1]
    dz = z[0,0,1] - z[0,0,0]
    nz = z.shape[2]
    #print dx,nx,dy,ny,dz,nz
    for z_k_at_k, xi, yi, zi in np.nditer([z_k_3d, x, y, z], 
           op_flags=[['readwrite'], ['readonly'], ['readonly'], ['readonly']]):
        # using dk = (2*pi)/(N*dr), we get k_i = x_i * (2*pi)/(N*dr**2)
        xki = xi*(np.pi*2)/(nx*(dx)**(2))
        yki = yi*(np.pi*2)/(ny*(dy)**(2))
        zki = zi*(np.pi*2)/(nz*(dz)**(2))
        k = math.sqrt(xki**2+yki**2+zki**2)
        z_k_at_k[...]  = zk_interp(k)
    return z_k_3d
    


def main(argv):
    args = process_command_line(argv)
    g_s, u_s, x, y, z, voxel = load_files(args)
    density, multiplicities, T, k, z_k = get_bulk_solvent_parameters(args)
    print 'Number of dxs',len(g_s)
    print 'Density:',density
    print 'Mult:',multiplicities
    print 'T:',T
    volume = x.shape[0]*y.shape[1]*z.shape[2]*voxel
    F_id = F_ideal_int(g_s, density, multiplicities)*kB*T*voxel
    F_ext = F_ext_int(g_s, u_s, density, multiplicities)*voxel*kB*T
    F_double = 0
    interaction_pair = 1
    for i in range(len(multiplicities)):
        for j in range(i+1):
            print 'Convolving',i,j
            z_3d = make_z_3d(k, z_k[:, interaction_pair], x, y, z)
            dbl_integral = F_exc_int_3d(g_s[i] - 1., u_s[i], g_s[j] - 1.,
                                        z_3d, voxel)
            if i != j:
                dbl_integral *= 2
            #dbl_integral = dbl_integral*multiplicities[j]*multiplicities[j]*density**2
            dbl_integral = dbl_integral*density**2
            print dbl_integral
            F_double += dbl_integral
            interaction_pair += 1
    F_exc = -F_double*kB*T/2.*voxel
    print ''
    print 'F_ideal: {:.8f}'.format(F_id)
    print 'F_ext: {:.8f}'.format(F_ext)
    print 'F_exc: {:.8f}'.format(F_exc)
#    print 'F_mod: {:.8f}'.format(F_mod)
#    print 'F_mod2: {:.8f}'.format(F_mod2)
#    print 'F_mod3: {:.8f}'.format(F_mod3)
    print 'F_total: {:.8f}'.format(F_id + F_ext + F_exc)
#    with open('decomp.txt', 'w') as f:
#        f.write('F_ideal: {:.8f}\n'.format(F_id))
#        f.write('F_ext: {:.8f}\n'.format(F_ext))
#        f.write('F_exc: {:.8f}\n'.format(F_exc))
#        f.write('F_total: {:.8f}\n'.format(F_id + F_ext + F_exc))
        
        
if __name__ == '__main__':
    main(sys.argv[1:])
