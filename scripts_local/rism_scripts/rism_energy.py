#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  2 15:36:06 2015

@author: max
"""


from __future__ import print_function

import argparse
#from chemistry.amber.readparm import AmberParm

import sys
import numpy as np

##### my modules ######
from comp_chem.rism_recipes import load_dx

MIN_RADIUS = 0.002 # A, minimum distance how close partilces can get together
# taken from rism3d_potential



def process_command_line(argv):
    """
    Processes arguments and returns namespace of them
    """
    parser = argparse.ArgumentParser(description=""" Compute lj contribution
                        to total potential energy.""")
    #Positional args
    parser.add_argument('guv', metavar='guv.dx',
                        help="""Guv dx files.""", nargs=2)
    #Optional args
#    parser.add_argument('-u', '--uuv',
#                        help=""" Potential energy file containig potential
#                        for which computation is needed. Sphere options
#                        are going to be ignored""")
    parser.add_argument('-t', '--huv',
                        help=""" Total correlation function instead of g(r).""",
                        action='store_true')
    parser.add_argument('-c', '--coords',
                        help=""" pdb file.""")
    parser.add_argument('-p', '--prmtop',
                        help=""" prmtop file. Doesn't work!""")
    parser.add_argument('--noelec',
                        help=""" Skip electrostatic energy.""",
                        action='store_true')
    parser.add_argument('--nolj',
                        help=""" Skip vdw energy.""",
                        action='store_true')
    parser.add_argument('-e', '--eps',
                        help=""" LJ epsilon of sphere in the box [0].""",
                        type=float, default=0)
    parser.add_argument('-q', '--charge',
                        help=""" Charge in the scpecified centre [1].""",
                        type=float, default=1.)
    parser.add_argument('-s', '--sigma',
                        help=""" LJ sigma of sphere in the box [0].""",
                        type=float, default=0)
    parser.add_argument('-d', '--density',
                        help=""" Density of solvent. Default is 298.15
                        water density. [3.3328311138810005E-02].""",
                        type=float, default=3.3328311138810005E-02)            
#    parser.add_argument('-t', '--temperature',
#                        help=""" Temperature. Only needed if uuv is provided
#                        [298.15].""",
#                        type=float, default=298.15)  
#    parser.add_argument('-r', '--radius',
#                        help=""" Complete energy contribution within given
#                        radius.""",
#                        type=float)  
#    parser.add_argument('-v', '--site',
#                        help=""" Name of cspc/e water site {O, H}. By
#                        default computes both.""")
    parser.add_argument('--centre', 
                        help=""" Centre of the lj sphere [0., 0., 0.].""",
                        default=[0., 0., 0.], type=float, nargs=3)
    return parser.parse_args(argv)
    

def calculate_pot_coefs(sigma, eps, q, site):
    sigma_mol, eps_mol, q_mol = sigma, eps, q
    #cspc/e parmas
    if site == 'O':
        sigma_sol = 3.1658  #A
        eps_sol   = 0.1553    #kcal/mol
        q_sol     = -0.8476   #e
    elif site == 'H':
        sigma_sol = 1.1658  #A
        eps_sol   = 0.01553     #kcal/mol    
        q_sol     = 0.4238    #e
    else:
        raise ValueError('Unknown water site')
    # using lorentz-berthelot rules
    sigma = (sigma_sol + sigma_mol)/2.
    eps = np.sqrt(eps_sol*eps_mol)
    A = 4*eps*sigma**12
    B = 4*eps*sigma**6
    q2 = q_sol * q_mol
    vacuum_p = 8.854187817e-12 #F/m
    conver = 3.694707135311445e-8 # m*e^2/(F*A) to kcal/mol
    coef = 1./(4*np.pi*vacuum_p)*conver
    q2_cor = coef*q2
    return A, B, q2_cor


def get_distance_matrices(x, y, z, centre):
    """ Return ri and r6i"""
    min_radius2 = MIN_RADIUS**2
    r2 = (x - centre[0])**2 + (y - centre[1])**2 + \
         (z - centre[2])**2
    r2[r2 < min_radius2] = min_radius2
    ri  = 1/np.sqrt(r2)
    r6i = 1/r2**3
    return ri, r6i
    
    
def compute_pots(sigma, eps, chg, ri, r6i, args):
    """ given sigma [A], eps [kcal/mol], chg [electrons], ri [1/A], r6i [1/A^6]
    and args returns energy of oxygen and hydrogen (multiplied by 2 for H)."""
    A_o, B_o, q2_o  = calculate_pot_coefs(sigma, eps, chg, 'O')
    A_h, B_h, q2_h  = calculate_pot_coefs(sigma, eps, chg, 'H')    
    u_o, u_h = 0, 0
    if not args.nolj:
        u_o += A_o*r6i**2 - B_o*r6i
        u_h += A_h*r6i**2 - B_h*r6i
    if not args.noelec:
        u_o += q2_o*ri
        u_h += q2_h*ri
    return u_o, u_h*2
    
    
def get_parm_parameters(prmtop):
    """ Return list of atom sigmas [A], epsilons [kcal/mol], charges [electrons]"""
    parm = AmberParm(prmtop)
    charges = []
    sigmas = []
    epss = []
    for i, atom in enumerate(parm.atom_list):
        attyp, _, attchg = atom.attype, atom.atname, float(atom.charge)
        charges.append(attchg)
        nbidx = parm.LJ_types[attyp]
        # convert to sigma
        sigmas.append(float(parm.LJ_radius[nbidx - 1])*2./(2**(1./6)))
        epss.append(float(parm.LJ_depth[nbidx - 1]))
    return sigmas, epss, charges

    
def get_coords_and_parms(pdb, prmtop):
    """ return list of atom centres and list of (sigma, eps, chg)"""
    with open(pdb) as f:
        lines = f.readlines()
    atoms = []
    for l in lines:
        if l.startswith('ATOM'):
            l = l.split()
            atoms.append(l[5:8])
    atoms = np.array(atoms, dtype=float)
    radii, epss, charges = get_parm_parameters(prmtop)
    parms = zip(radii, epss, charges)
    return atoms, parms


def main(argv):
    r""" print 
    kT * \int g(r)*u(r)*\rho*dV
    in kcal/mol
    """
    args = process_command_line(argv)
    k = 1.9872041E-3  # R in kcal/mol/K
    # load g(r)
    g_o, x, y, z  = load_dx(args.guv[0], True)
    dx = x[1, 0, 0] - x[0, 0, 0]
    dy = y[0, 1, 0] - y[0, 0, 0]
    dz = z[0, 0, 1] - z[0, 0, 0]
    g_h, _, _ = load_dx(args.guv[1], False)
    if args.huv:
        g_o = g_o + 1
        g_h = g_h + 1
    # calc pot
#    if args.uuv:
#        print('Using uuv file')
#        u = load_dx(args.uuv)[0]
    if args.coords:
        #print('Using pdb')
        atoms, parms = get_coords_and_parms(args.coords, args.prmtop)
        u_o, u_h = 0, 0
        for centre, (sigma, eps, chg) in zip(atoms, parms):
            ri, r6i = get_distance_matrices(x, y, z, centre)
            u_o_atom, u_h_atom = compute_pots(sigma, eps, chg, ri, r6i, args)
            u_o += u_o_atom
            u_h += u_h_atom
    else:
        #print('Assuming potential is relative to a single point at given centre')
        ri, r6i = get_distance_matrices(x, y, z, args.centre)
        u_o, u_h = compute_pots(args.sigma, args.eps, args.charge, ri, r6i, args)
    voxel = dx*dy*dz
    coef = voxel*args.density
#    if args.uuv:
#        coef  = k*args.temperature*coef
    energy_o = coef*np.sum(g_o*u_o)
    energy_h = coef*np.sum(g_h*u_h)
    print('Total: {:.4f} energy O: {:.4f} energy H: {:.4f}'.format(energy_o + energy_h,
          energy_o, energy_h))
    

if __name__ == '__main__':
    main(sys.argv[1:])

    
    
    
