#!/usr/bin/python

# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 10:56:46 2015

@author: max
"""

from mayavi import mlab
import numpy as np
import sys
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
    parser = argparse.ArgumentParser(description="""Vizualize volumetric distributions.""")
    #Optional args
    parser.add_argument('-g', '--guv',
                        help=""" Spatial distribution functions.""", nargs='+')
    parser.add_argument('-c', '--cuv',
                        help=""" Direct correlation functions.""", nargs='+')
    parser.add_argument('-f', '--free_energy',
                        help=""" Compute hnc free energy.""",
                        action='store_true')
    parser.add_argument('--o_contrib',
                        help=""" Consider oxygen radial functions.""",
                        action='store_true')
    parser.add_argument('--h_contrib',
                        help=""" Consider hydrogen radial functions.""",
                        action='store_true')
    parser.add_argument('-p', '--pdb',
                        help=""" pdb molecule.""")
    parser.add_argument('-t', '--prmtop',
                        help=""" Amber topology.""")
    parser.add_argument('--invert',
                        help=""" Plot -vv.""",
                        action='store_true')
    parser.add_argument('-ls', '--lskip_part',
                        help=""" Part of the data at the beginning that is going to
                        be skipped.""", default=.2, type=float)
    parser.add_argument('-rs', '--rskip_part',
                        help=""" Part of the data at the end that is going to
                        be skipped.""", default=.2, type=float)
    parser.add_argument('-ln', '--lskipnum_part',
                        help=""" Minimal value of the data to visualize.""",
                        type=float)
    parser.add_argument('-rn', '--rskipnum_part',
                        help=""" Maximal value of the data to visualize.""",
                        type=float)
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


def paint_atoms(names):
    colours = []
    for name in names:
        if 'Cl' in name:
            colours.append((0, 1, 0)) # green
        elif 'C' in name:
            colours.append((0, 0, 0))  # black
        elif 'O' in name:
            colours.append((1, 0, 0))  # red
        elif 'H' in name:
            colours.append((1, 1, 1))  # white
        elif 'N' in name:
            colours.append((0, 0, 1))  #blue
        elif 'P' in name:
            colours.append((1, 140/255, 0)) #dark orange
        elif 'S' in name:
            colours.append((1, 1, 0)) # yellow            
        elif 'Br' in name:
            colours.append((205/255,92/255,92/255))  # dark red
        elif 'I' in name:
            colours.append((160/255,32/255,240/255)) #pruple
        elif 'F' in name:
            colours.append((0,1,1))  #cyan
        else:
            colours.append((190/255,190/255,190/255))  #gray
    return colours


def visualize_molecule(args):
    with open(args.pdb) as f:
        lines = f.readlines()
    names = []
    atoms = []
    for l in lines:
        if l.startswith('ATOM'):
            l = l.split()
            names.append(l[2])
            atoms.append(l[5:8])
    atoms = np.array(atoms, dtype=float)
    colours = paint_atoms(names)
    #get vdw and plot
    if args.prmtop:
        parm = AmberParm(args.prmtop)
        radii = []
        for atom in parm.atom_list:
            nbidx = parm.LJ_types[atom.attype]
            radii.append(float(parm.LJ_radius[nbidx - 1]))
    else:
        radii = [1. for i in atoms]
    for (x, y, z), r, c in zip(atoms, radii, colours):
        #print('visualizing atom')
        mlab.points3d(x, y, z, scale_factor=r*1.25,
                      resolution=20, color=c)
    # draw bonds
#    mlab.plot3d(atoms[:, 0], atoms[:,1], atoms[:,2],
#            tube_radius=0.4, colormap='Reds')



def get_free_energy(args):
    print('Assuming water solvent')
    g_o, x, y, z = load_dx(args.guv[0], True)
    g_h = load_dx(args.guv[1])[0]
    c_o = load_dx(args.cuv[0])[0]
    c_h = load_dx(args.cuv[1])[0]
    o_fe = .5*(g_o-1)**2 - c_o - .5*(g_o - 1)*c_o
    h_fe = .5*(g_h-1)**2 - c_h - .5*(g_h - 1)*c_o
    if args.o_contrib:
        return o_fe, x, y, z
    if args.h_contrib:
        return h_fe, x, y, z
    return 2*h_fe + o_fe, x, y, z
    

def main(argv):
    args = process_command_line(argv)
    if not args.free_energy:
        vv, x, y, z = load_dx(args.guv[0], True)
    else:
        vv, x, y, z = get_free_energy(args)
        
    if args.pdb:
        visualize_molecule(args)
    if args.invert:
        print('Inverted data!')
        vv = -vv
    src = mlab.pipeline.scalar_field(x, y, z, vv)
    print('Smallest value in volumetric data: {}'.format(vv.min()))
    print('Largest value in volumetric data: {}'.format(vv.max()))
    if args.lskipnum_part != None:
        vmin = args.lskipnum_part
    else:
        vmin = vv.min() + args.lskip_part*vv.ptp()
    print('Darkest colour on map corresponds to: {}'.format(vmin))
    if args.rskipnum_part:
        vmax = args.rskipnum_part
    else:
        vmax=vv.max() - args.rskip_part*vv.ptp()
    print('Warmest colour on map corresponds to: {}'.format(vmax))    
    mlab.pipeline.volume(src, vmin=vmin, vmax=vmax)
    mlab.show()
    


if __name__ == '__main__':
    main(sys.argv[1:])
