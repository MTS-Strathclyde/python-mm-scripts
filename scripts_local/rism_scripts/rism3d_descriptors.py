#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 22 14:25:17 2015

@author: max
"""

from __future__ import print_function

import sys
import os
import numpy as np
from chemistry.amber.readparm import AmberParm
from scipy.spatial import cKDTree

G_NAME = 'g_{name}.{aname}.1.dx'
C_NAME = 'c_{name}.{aname}.1.dx'
ACR_NAME = 'a_{name}cr.1.dx'


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
        # substraction is done to avoid weird floating point errors
        x = np.arange(OrX, OrX + dX*size[0] - dX/5., dX)
        y = np.arange(OrY, OrY + dY*size[1] - dY/5., dY)
        z = np.arange(OrZ, OrZ + dZ*size[2] - dZ/5., dZ)
        x_mat, y_mat, z_mat = np.meshgrid(x, y, z, indexing='ij')
        return dx_m.astype('float'), x_mat, y_mat, z_mat
    else:
        return dx_m.astype('float'), size, (dX, dY, dZ)    


def get_rmin2_radii(prmtop):
    """ Return a list of r_min/2 values """
    try:
        parm = AmberParm(prmtop)
        radii = []
        for atom in parm.atom_list:
            nbidx = parm.LJ_types[atom.attype]
            radii.append(float(parm.LJ_radius[nbidx - 1]))
    except IndexError:
        print('Couldnt open {} using AmberParm'.format(prmtop))
        print('Assuming single atom topology')
        with open(prmtop) as f:
            for l in f:
                if 'JONES_ACOEF' in l:
                    f.next()
                    acoef = float(f.next())
                    f.next()
                    f.next()
                    bcoef = float(f.next())
                    radii = [(acoef/bcoef*2)**(1./6)/2]
    return radii


def get_vdw_shell_regions(pdb, prmtop, x, y, z, size, shell_boundaries):
    """ 
    Parameters
    ----------
    x : 3d-matrix
    y : 3d-matrix
    z : 3d-matrix
    size : tuple
        shape of the matrices
    shell_boundaries : list
        list of tuples containing minimum and maximum radiuses of shells
        
    """
    with open(pdb) as f:
        lines = f.readlines()
    atoms = []
    for l in lines:
        if l.startswith('ATOM'):
            l = l.split()
            atoms.append(l[5:8])
    atoms = np.array(atoms, dtype=float)
    radii = get_rmin2_radii(prmtop)
    shell_regions = []
    xyz = np.c_[x.flatten(), y.flatten(), z.flatten()]
    ktree = cKDTree(xyz)
    for s_min, s_max in shell_boundaries:
        mask_outer = []
        mask_inner = []
        for atom, r in zip(atoms, radii):
            outer_r = r + s_max
            inner_r = r + s_min
            # negative radii don't make much sense
            if outer_r < 0:
                outer_r = 0
            if inner_r < 0:
                inner_r = 0
            outer_indxs = ktree.query_ball_point(atom, outer_r)
            inner_indx = ktree.query_ball_point(atom, inner_r)
            mask_outer.extend(outer_indxs)
            mask_inner.extend(inner_indx)
        mask_outer = set(mask_outer)
        mask_inner = set(mask_inner)
        shell_mask = list(mask_outer - mask_inner)
        # create shrot range from indexes
        shell_idxs = np.zeros(xyz.shape[0], dtype='bool')  #all false array
        shell_idxs[shell_mask] = True
        shell_regions.append(shell_idxs.reshape(size))
    return shell_regions

    
def get_grids(p, name, atomname, atomcharge):
    """ Return x, y, z, h, c, [c_short],
    given molpath, molname, atname, atchg."""
    g_path = os.path.join(p, G_NAME.format(name=name, aname=atomname))
    c_path = os.path.join(p, C_NAME.format(name=name, aname=atomname))
    acr_path = os.path.join(p, ACR_NAME.format(name=name))
    g, x, y, z = load_dx(g_path, True)
    h  = g - 1
    c, size, step_vector = load_dx(c_path)
    #voxel = step_vector[0]*step_vector[1]*step_vector[2]
    if atomcharge:
        acr, _, _ = load_dx(acr_path)
        c_short = c - atomcharge*acr    
        return x, y, z, size, h, c, c_short
    else:
        return x, y, z, size, h, c        
    
    
def get_avg_std_max_min(array):
    if len(array) > 0:
        return np.average(array), np.std(array), np.amax(array) - np.amin(array)    
    else:
        raise ValueError('Empty array - sth went wrong')
    
    
def get_descriptors(shell_region, h, c, c_short=None):
    """ For each region compute average, std, max-min. """
    h_sr = h[shell_region]
    c_sr = c[shell_region]
    fe_sr = .5*h_sr**2 - c_sr - .5*h_sr*c_sr
    sr_names = ['h', 'c', 'fe']
    sr_points = [h_sr, c_sr, fe_sr]
    if c_short:
        c_short_sr = c_short[shell_region]
        fe_short_sr = .5*h_sr**2 - c_short_sr - .5*h_sr*c_short_sr
        sr_names.extend(['c_short', 'fe_short'])
        sr_points.extend([c_short_sr, fe_short_sr])
    descriptors = []
    d_names = []
    for sr_name, sr in zip(sr_names, sr_points):
        descriptors.extend(get_avg_std_max_min(sr))
        d_names.extend(('avg_' + sr_name, 'std_' + sr_name, 'minmax_' + sr_name ))
    return d_names, descriptors
    

def get_shell_boundaries(vdw_distances):
    """ return shell boundaries --- list of tuples, where each tuple contains
    beginning and end distance of each shell relative to atom vdw radius."""
    # add -5000 that ensure that first boundary goes from 0 
    vdw_distances = [-5000] + vdw_distances
    shell_boundaries = []    
    for i in range(1, len(vdw_distances)):
        shell_boundaries.append((vdw_distances[i-1], vdw_distances[i]))
    vdw_distances.pop(0)
    # create regions of shell idxs
    return shell_boundaries


def get_solnames_solcharges(xvv):
    """ return names of solvent atoms and chgs of solvent atoms
    charges are in amber prmtop format (e chgs multiplied by 18.2223) """
    solnames = []
    solcharges = []
    with open(xvv) as f:
        for l in f:
            if l.startswith('%FLAG ATOM_NAME'):
                f.next()
                solnames = f.next().split()
            if l.startswith('%FLAG QV'):
                f.next()
                solcharges = map(float, f.next().split())
                break
    return solnames, solcharges


def compute_all_shell_descriptors(p, name, vdw_distances, xvv):
    """ given path to directory containing all files;
    molecule base name (without extension)
    distances from vdw surfaces defining shells,
    xvv file of solvent
    ------
    compute descriptros for each atom in each shell.
    """
    #print(p, name, vdw_distances, xvv)
    pdb = os.path.join(p, name + '.pdb')
    prmtop = os.path.join(p, name + '.prmtop')
    solatoms, solcharges = get_solnames_solcharges(xvv)
    # for each solatom solgrid contains : x, y, z, size, h, c, [c_short]
    # c_short would be used only if atomcharge != 0
    solgrids  = []
    for atomname, atomcharge in zip(solatoms, solcharges):
        solgrids.append(get_grids(p, name, atomname, atomcharge))
    # !! we assume that x==y==z==size for all solatoms to save time here
    # that's why we load x y z size for only one solvent atom
    x, y, z, size = solgrids[0][:4]
    shell_boundaries = get_shell_boundaries(vdw_distances)
    shell_regions = get_vdw_shell_regions(pdb, prmtop, x, y, z, size, 
                                          shell_boundaries)
    # now compute discriptro for each atom and for each shell region
    all_d_names = []
    all_d = []
    for solatom_name, solchg, solgrid in zip(solatoms, solcharges, solgrids):
        for (s_min, s_max), sr in zip(shell_boundaries, shell_regions):
            if solchg:
                h, c, c_short = solgrid[4:]
            else:
                h, c = solgrid[4:]
                c_short = None
            d_names, descriptors = get_descriptors(sr, h, c, c_short)
            # write full names of descriptors
            d_names_full = []
            for d_name in d_names:
                full_name = '{}_{}_{}_to_{}'.format(solatom_name, d_name,
                                                     s_min, s_max)
                d_names_full.append(full_name)
            all_d_names.extend(d_names_full)
            all_d.extend(descriptors)
    # finished calculating descriptors
    return all_d_names, all_d
    

def main(argv):
    print('USAGE:')
    print('path name shell_distance1 [shell_distance2 shell_distance3 ...]  xvv')
    print(' ')
    p = argv[0]
    name = argv[1]
    vdw_d = map(float, argv[2:-1])
    xvv = argv[-1]
    all_d_names, all_d = compute_all_shell_descriptors(p, name, vdw_d, xvv)
    for name, d in zip(all_d_names, all_d):
        print('{:<50} {: .3e}'.format(name, d))


if __name__== '__main__':
    main(sys.argv[1:])
    



