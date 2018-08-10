#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 15:51:33 2015

@author: max
"""

from __future__ import print_function, division

import sys
import numpy as np
from scipy.spatial import cKDTree
import argparse
from chemistry.amber.readparm import AmberParm


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
    parser.add_argument('pdb', help="""PDB file""")
    parser.add_argument('prmtop', help="""prmtop file""")
    #Optional args
    parser.add_argument('-g', '--guv',
                        help=""" Spatial distribution functions.""", nargs='+')
    parser.add_argument('-c', '--cuv',
                        help=""" Direct correlation functions.""", nargs='+')
    parser.add_argument('-r', '--vdw_multiple',
                        help=""" Multiple of vdw considered as radius. [1]""", default=1., type=float)
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
        return dx_m.astype('float'), size, (dX, dY, dZ)    


def get_shrot_region(args, x, y, z, size):
    with open(args.pdb) as f:
        lines = f.readlines()
    atoms = []
    for l in lines:
        if l.startswith('ATOM'):
            l = l.split()
            atoms.append(l[5:8])
    atoms = np.array(atoms, dtype=float)
    parm = AmberParm(args.prmtop)
    radii = []
    for atom in parm.atom_list:
        nbidx = parm.LJ_types[atom.attype]
        radii.append(float(parm.LJ_radius[nbidx - 1]))
    # convert meshgrid to a matrix and find indexes of all points
    # which are inside molecule
    xyz = np.c_[x.flatten(), y.flatten(), z.flatten()]
    ktree = cKDTree(xyz)
    mask = []
    for atom, r in zip(atoms, radii):
        excluded_r = r * args.vdw_multiple
        excl_indxs = ktree.query_ball_point(atom, excluded_r)
        mask.extend(excl_indxs)
    mask = list(set(mask))
    # create shrot range from indexes
    short_range = np.zeros(xyz.shape[0], dtype='bool')  #all false array
    short_range[mask] = True
    return short_range.reshape(size)


def main(argv):
    args = process_command_line(argv)
    fe_prefactor = 1.9872041/1000*298.15*3.3328311138810005E-02
    
    g_o, x, y, z = load_dx(args.guv[0], return_xyz=True)
    g_h, size, (dx,dy,dz) = load_dx(args.guv[1])
    c_o = load_dx(args.cuv[0])[0]
    c_h = load_dx(args.cuv[1])[0] 
    cors = [(g_o, c_o, 'O'), (g_h, c_h,'H')]
    
    short_region = get_shrot_region(args, x, y, z, size)
    long_region = np.invert(short_region)

    short_total = 0
    long_total = 0
    for g, c, name in cors:
        short_g = g[short_region]
        #print(short_g.shape)
       # print(short_g)
        short_c = c[short_region]
        long_g = g[long_region]
        long_c = c[long_region]
        fe_short = np.sum(.5*(short_g-1)**2 - short_c - .5*(short_g - 1)*short_c)*dx*dy*dz*fe_prefactor
        fe_long = np.sum(.5*(long_g-1)**2 - long_c - .5*(long_g - 1)*long_c)*dx*dy*dz*fe_prefactor
        if name == 'H':
            fe_short = fe_short*2
            fe_long = fe_long*2
        short_total += fe_short
        long_total += fe_long
        print('{} contributes {:.4f} in short region and {:.4f} in long region.'.format(name, fe_short, fe_long))
    print()
    print('Short region contributes {:.4f} to total free energy, long region {:.4f}'.format(short_total, long_total))
    print()
    print('calculated total hnc_fe: {}'.format(short_total + long_total))
        
        
    
if __name__ == '__main__':
    main(sys.argv[1:])


