#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 28 19:08:51 2016

@author: max
"""
from __future__ import print_function

from scipy import optimize
from scipy import ndimage
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import numpy.ma as ma
import sys
import math
from scipy import interpolate as interp
import argparse

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
    parser = argparse.ArgumentParser(description="""Use sergievskis formula
                        to minimize rism.""")
    parser.add_argument('z_k',
                        help=""" z matrix.""")
    parser.add_argument('-g', '--guv',
                        help=""" Pair correlation functions.""", nargs='+')
    parser.add_argument('-t', '--huv',
                        help=""" Total distribution functions.""", nargs='+')

    parser.add_argument('-u', '--uuv',
                        help=""" Potential energy dist.""", nargs='+')
    parser.add_argument('-d', '--rho0_s',
                        help=""" Solvent sites number densities. 
                        [0.03332, 0.06664]""",
                         default=[0.03332, 0.03332*2], type=float, nargs='+')
    parser.add_argument('-m', '--multiplicities',
                        help="""Solvent sites multiplicities.""", nargs='+',
                        type=int)

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

def F_ideal_int(rho_s, rho0_s):
    rho_s_ma = np.array([ma.masked_values(rho_i/rho0_i,0,atol=1.0e-15) for rho_i, rho0_i \
                           in zip(rho_s, rho0_s)])
    log_s = [ma.log(rho_i_ma) for rho_i_ma in rho_s_ma]
    integral = sum([rho_i * ma.filled(log_i, 0) - (rho_i - rho0_i) \
                        for rho_i, rho0_i, log_i in zip(rho_s, rho0_s, log_s)])
    return np.sum(integral)
    
    
def F_ext_int(rho_s, u_s):
    integrand = sum([rho_i*u_i for rho_i, u_i in zip(rho_s, u_s)])
    return np.sum(integrand)
        
    

def F_exc_int_3d(h1, h2, z_3d):
    """ zk_3d  version"""
    h2_k = np.fft.fftn(h2)
    h2_k = np.fft.fftshift(h2_k)  # shift 0 freq to the middle
    convol = np.fft.ifftn(np.fft.ifftshift(h2_k*z_3d))
    return np.sum(h1*np.real(convol))


def F_cor(rho_s, rho0_s, kBT, voxel, args):
    rho_O  = rho_s[0]
    bar_V = -np.sum(rho_O-rho0_s[0])*voxel/rho0_s[0]
    print('pmv',bar_V)
    return bar_V*0.122324502591


def F_energy(rho_s, rho0_s, u_s, z_3ds,  kBT, voxel, args, shape=False, verbose=False):
    if np.any(shape):
        rho_s = rho_s.reshape(shape)
    verboseprint = print if verbose else lambda *a, **k: None
    F_id = F_ideal_int(rho_s, rho0_s)*kBT*voxel
    F_ext = F_ext_int(rho_s, u_s)*voxel*kBT
    F_double = 0
    interaction_pair = 0
    verboseprint('  Computing double int...')
    for i in range(len(rho0_s)):
        for j in range(i + 1):
            verboseprint('    Interaction',i,j)
            dbl_integral = F_exc_int_3d(rho_s[i] - rho0_s[i],
                                        rho_s[j] - rho0_s[j], z_3ds[i+j])
            if i != j:
                dbl_integral *= 2
            F_double += dbl_integral/args.multiplicities[i]/args.multiplicities[j]
            verboseprint('    Energy: ',dbl_integral/args.multiplicities[i]/args.multiplicities[j])
            verboseprint()
            interaction_pair += 1
    F_exc = -F_double*voxel*kBT/2.
    verboseprint('  F_id',F_id)
    verboseprint('  F_ext',F_ext)
    verboseprint('  F_exc',F_exc)
    verboseprint()
    #print('energy',F_id + F_exc + F_ext - F_cor(rho_s, rho0_s, kBT, voxel, args))
    return F_id + F_exc + F_ext# - F_cor(rho_s, rho0_s, kBT, voxel, args)
    
    
def F_derviative(rho_s, rho0_s, u_s, z_3ds, kBT, voxel, args, shape=None):
    if np.any(shape):
        rho_s = np.reshape(rho_s, shape)
    rho_s_ma = [ma.masked_values(rho_i/rho0_i, 0, atol=1.0e-15) for rho_i, rho0_i in\
                                                           zip(rho_s, rho0_s)]
    log_s = [ma.log(rho_i_ma) for rho_i_ma in rho_s_ma]
    deriv = []
    for i, (log_i, u_i, rho_i, z_3d_i) in\
                                     enumerate(zip(log_s, u_s, rho_s, z_3ds)):
        mean_f_part = kBT*(ma.filled(log_i, -34.53877639491068) + u_i)
        ex_part = 0
        for j, (rho_j, rho0_j, z_3d_ij) in enumerate(zip(rho_s, rho0_s, z_3ds[i,:])):
            delta_rho = rho_j - rho0_j
            delta_rho_k = np.fft.fftn(delta_rho)
            delta_rho_k = np.fft.fftshift(delta_rho_k)  # shift 0 freq to the middle
            convol = np.real(np.fft.ifftn(np.fft.ifftshift(delta_rho_k*z_3d_ij)))
            ex_part += kBT*convol/args.multiplicities[i]
        deriv_i = mean_f_part - ex_part
        #fix discontinuty
        deriv_i = np.where(deriv_i > 500, 500, deriv_i)
        deriv.append(deriv_i)
        print('median',np.median(deriv_i))
        print('avg',np.average(deriv_i))
    return np.array(deriv) #- 0.122324502591
    
    
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
            xj, yj, zj = convol_it[1], convol_it[2], convol_it[3]
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


def F_derviative_slow(rho_s, rho0_s, u_s, zr_s, r, x, y, z, kBT, voxel, shape=None):
    """ Derivative of 3D-RISM functional:
    F'_3drism_i  = kBT*np.log(g_i (r)) + u_i(r) - kBT*sum(z_ij * h_j)
    
    return array of derivative with the shape n_s * grid
    """
    if np.any(shape):
        rho_s = np.reshape(rho_s, shape)
    if len(zr_s.shape) == 1:
        zr_s = zr_s.reshape((-1, 1))
    rho_s_ma = [ma.masked_values(rho_i/rho0_i, 0, atol=1.0e-15) for rho_i, rho0_i in\
                                                           zip(rho_s, rho0_s)]
    log_s = [ma.log(rho_i_ma) for rho_i_ma in rho_s_ma]
    deriv = []
    for i, (u_i, rho_i, log_i) in enumerate(zip(u_s, rho_s, log_s)):
        # we multipy u_i by kbt as well as it is actually u_i*beta
        mean_f_part = kBT*(ma.filled(log_i, -34.53877639491068) + u_i)
        ex_part = 0
        for rho_j, rho0_j, zr in zip(rho_s, rho0_s, zr_s.T):
            delta_rho = rho_j - rho0_j
            #convol = ndimage.filters.convolve(delta_rho, z_3d_ij, mode='constant')
            convol = slow_convolution(delta_rho, zr, r, x, y, z, voxel)
            ex_part += convol
        deriv.append(mean_f_part + ex_part*kBT)
    return np.array(deriv)
    
    
def check_Z(rho_s, rho0_s, u_s, zk_3ds, kBT, voxel, args, shape=None):
    """ Checks validity of Z by computing g(r) from it.
    To do it uses the following relationship:
    g_i(r) = exp(-u_i/kt + \sum_j \rho_j (z_{ij}*h_j) (r))
    where * stands for convolution
    """
    if np.any(shape):
        rho_s = np.reshape(rho_s, shape)
    g_r = []
    for i, (u_i, rho_i) in enumerate(zip(u_s, rho_s)):
        u_part = -u_i
        ex_part = 0
        for j, (rho_j, rho0_j) in enumerate(zip(rho_s, rho0_s)):
            interaction_pair = i + j
            print('Check_Z convolution')
            print('Checking pair',interaction_pair)
            zk_3d = zk_3ds[interaction_pair]
            delta_rho = rho_j - rho0_j
            delta_rho_k = np.fft.fftn(delta_rho)
            delta_rho_k = np.fft.fftshift(delta_rho_k)  # shift 0 freq to the middle
            convol = np.real(np.fft.ifftn(np.fft.ifftshift(delta_rho_k*zk_3d)))
            ex_part += convol/args.multiplicities[j]/args.multiplicities[i]
#            print('Convolution at pair',interaction_pair,i,j)
            print(np.sum(convol))
            print(np.sum(ex_part))
#            raw_input('paused')
#        print('u_r')
#        print(u_part)
#        print('!!!!!!!!!!!!!!!!!!!!!!!!')
#        print('Sum of convolutions')        
#        print(ex_part)
#        print('!!!!!!!!!!!!!!!!!!!!!!!!')
#        print('total')
#        print(u_part+ex_part)
#        raw_input('paused')
        g_r.append(np.exp(u_part + ex_part))
    return g_r
    

def finite_dif_derivative(rho_s, rho0_s, u_s, z_3ds,  kBT, voxel, 
                          shape, args, eps=1.0e-6):
    rho_s  = rho_s.reshape(shape)
    finite_deriv = np.zeros_like(rho_s)
    it = np.nditer(rho_s, op_flags=['readwrite'], flags=['multi_index'])
    while not it.finished:
        rho_i = it[0]
        idx = it.multi_index
        orig_rho = float(rho_i)
        # plus/minus eps
        rho_i[...] = orig_rho + eps
        e_plus = F_energy(rho_s, rho0_s, u_s, z_3ds,  kBT, voxel, args)
        rho_i[...] = orig_rho - eps
        e_minus = F_energy(rho_s, rho0_s, u_s, z_3ds,  kBT, voxel, args)
        # deriv
        finite_deriv[idx] = (e_plus - e_minus)/(2*eps)
        it.iternext()
    return finite_deriv 
           
    
def F_opt(rho_s, rho0_s, u_s, z_3ds, kBT, voxel, shape, args, verbose=False):
    rho_s = np.reshape(rho_s, shape)
    energy = F_energy(rho_s, rho0_s, u_s, z_3ds,  kBT, voxel, args,verbose=False)
    deriv = F_derviative(rho_s, rho0_s, u_s, z_3ds, kBT, voxel, args)
    #deriv = finite_dif_derivative(rho_s, rho0_s, u_s, z_3ds, kBT, voxel,shape)
    return energy, deriv.flatten()

    
def make_z_3d(k, z_k, x, y, z):
    """ Makes 3d cube out of z, centered in the origin. """
    zk_interp = interp.interp1d(k, z_k, kind='nearest')
    z_k_3d = np.zeros_like(x)
    # using dk = (2*pi)/(N*dr), we get k_i = x_i * (2*pi)/(N*dr**2)
    dx = x[1,0,0] - x[0,0,0]
    nx = x.shape[0]
    dx_dk = (np.pi*2)/(nx*(dx)**(2))
    dy = y[0,1,0] - y[0,0,0]
    ny = y.shape[0]
    dy_dk = (np.pi*2)/(ny*(dy)**(2))
    dz = z[0,0,1] - z[0,0,0]
    nz = z.shape[0]
    dz_dk = (np.pi*2)/(nz*(dz)**(2))
    #print dx,nx,dy,ny,dz,nz
    for z_k_at_k, xi, yi, zi in np.nditer([z_k_3d, x, y, z], 
           op_flags=[['readwrite'], ['readonly'], ['readonly'], ['readonly']]):
        xki = xi*dx_dk
        yki = yi*dy_dk
        zki = zi*dz_dk
        k = math.sqrt(xki**2+yki**2+zki**2)
        z_k_at_k[...]  = zk_interp(k)
    return z_k_3d


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
        if args.huv:
            g_i += 1
        #g_i = np.ones_like(g_i)
        g_s.append(g_i)
        u_s.append(u_i)
    voxel = dx*dy*dz
    rho_s = []
    for g_i, rho0_i in zip(g_s, args.rho0_s):
        rho_s.append(g_i*rho0_i)
    rho_s = np.array(rho_s)
    return rho_s, u_s, x, y, z, voxel


def main(argv):
    args = process_command_line(argv)
    kBT = kB*T
    rho_s, u_s, x, y, z, voxel = load_files(args)
    rho0_s = np.array(args.rho0_s)
    z_k = np.loadtxt(args.z_k, dtype=float)
    k = z_k[:,0]
    zk_3ds = np.array([make_z_3d(k, z_k[:,i], x, y, z) for i in range(1,z_k.shape[1])])
    print(zk_3ds.shape)
#    print(np.sum(np.abs(zk_3ds[0] - make_z_3d(k, z_k[:,1], x, y, z))))
#    raw_input('Paused')
    print(F_energy(rho_s, rho0_s, u_s, zk_3ds, kBT, voxel, args, verbose=True))
    raw_input('Paused')
    #prepare rho_s for optimization
    shape = rho_s.shape
    rho_s = rho_s.flatten()
    print('Total number of points:',rho_s.size)
    print()
    print('Max u',np.amax(u_s))
    print()
#    print('Removing max u value')
#    
    print('Check Z')
    computed_g_r = check_Z(rho_s, rho0_s, u_s, zk_3ds, kBT, voxel, args, shape)
    g_r  = np.array([rho_i/rho0_i for rho_i, rho0_i  in zip(rho_s.reshape(shape), rho0_s)])
    for comp_g_r_i, g_r_i in zip(computed_g_r, g_r):        
        plt.figure()
        plt.plot(comp_g_r_i.flatten(), g_r_i.flatten(), 'o')
    plt.show()
    raw_input('Paused')
    print('Slow derivative')
    sderiv = F_derviative(rho_s, rho0_s, u_s, zk_3ds, kBT, voxel, args, shape)
#    print(sderiv)
    print('Sum of derivative',np.sum(sderiv))
    print('Largest derivative',np.amax(sderiv), 'At',np.argmax(sderiv))
    print('abs avg',np.mean(np.abs(sderiv)))
    raw_input('Paused')
    print('Evaluate starting guess...')
    ener, derivative = F_opt(rho_s, rho0_s, u_s, zk_3ds, kBT, voxel, shape, args, verbose=True)
    print('Energy',ener)
    raw_input('press enter to minimize')
    print('Derivative',derivative)
    print('Derivative sum',np.sum(derivative))
    print('Derivative shape',derivative.shape)
    print('')
#    print('Testing finite derivative')
#    f_deriv = finite_dif_derivative(rho_s, rho0_s, u_s, z_3ds,  kBT, voxel, shape, args, eps=1.0e-14)
#    print('Derivative',np.c_[f_deriv.flatten(), derivative, abs(f_deriv.flatten()-derivative)])
#    print('Derivative sum',np.sum(f_deriv))
#    print('Derivative shape',f_deriv.shape)
    
    bounds= []
    for i in range(rho_s.size):
        bounds.append((0, None))
#    cons = ({'type' : 'ineq',
#             'fun'  : })
    print('')
    print('Starting optimization')
#    opt = optimize.minimize(F_energy, rho_s,
#                            args=(rho0_s, u_s, z_3ds, kBT, voxel, shape),
#                            options=dict(disp=True), bounds=bounds)
    opt = optimize.minimize(F_opt, rho_s, jac=True,
                            args=(rho0_s, u_s, zk_3ds, kBT, voxel, shape,args),
                            options=dict(disp=True,ftol=1.0e-15,maxfun=1000,
                                         maxiter=1000),
                            bounds=bounds,)
    dif = np.abs(opt.x.reshape(shape)[0] - rho_s.reshape(shape)[0])
    print(dif[dif > 0.000001])
    print('MAE',np.sum(np.abs(opt.x.reshape(shape)[0] - rho_s.reshape(shape)[0])))
    
    
        
        
if __name__ == '__main__':
    main(sys.argv[1:])
