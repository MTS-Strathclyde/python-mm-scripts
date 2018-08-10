#!/usr/bin/python

# -*- coding: utf-8 -*-
"""
Created on Thu May 21 16:40:52 2015

@author: max
"""

from grid import grid

import matplotlib.pyplot as plt
import numpy as np
import argparse
import sys
import os


def process_command_line(argv):
    """
    Processes arguments and returns namespace of them
    """
    parser = argparse.ArgumentParser(description="""Convert volumetric dx file
                        to 1d, normalized rdf.""")
    #Positional args
    parser.add_argument('files', metavar='file.dx',
                        help="""Volumetric dx file.""", nargs='+')
    #Optional args
    parser.add_argument('-d', '--delta',
                        help="""Step in 1d rdf [0.01].""", type=float, default=0.01)
    parser.add_argument('-l', '--limit',
                        help="""Extent of 1d rdf [15].""", type=float, default=15)
    parser.add_argument('-c', '--centre', 
                        help="""Origin of 1d rdf [0., 0., 0.].""",
                        default=[0., 0., 0.], type=float, nargs=3)
    parser.add_argument('-p', '--plot', 
                        help=""" Plot rdf.""",
                        action='store_true')
    parser.add_argument('-n', '--name', 
                        help=""" output name.""")
    parser.add_argument('--potential',
                        help=""" Should contain 4 elements, 
                        indication whether we compute energy for O or for H
                        site of water, 2 lj paramters for solute, solute chg.
                        [O/H] Lj simga [A], lj eps [kcal/mol, elementary chg.
                        Assumes that solvent is cspc/e water.
                        """, nargs=4)
    return parser.parse_args(argv)


def calculate_pot_coefs(potential):
    #cspc/e parmas
    if potential[0] == 'O':
        sigma_sol = 3.1658  #A
        eps_sol   = 0.1553    #kcal/mol
        q_sol     = -0.8476   #e
    elif potential[0] == 'H':
        sigma_sol = 1.1658  #A
        eps_sol   = 0.046     #kcal/mol    
        q_sol     = 0.4238    #e
    else:
        raise ValueError('Unknown water site')
    sigma_mol, eps_mol, q_mol = map(float, potential[1:])
    # using lorentz-berthelot rules
    sigma = (sigma_sol + sigma_mol)/2.
    eps = np.sqrt(eps_sol*eps_mol)
    A = 4*eps*sigma**12
    B = 4*eps*sigma**6
    q2 = q_sol * q_mol
    return A, B, q2
    

def main(argv):
    args = process_command_line(argv)
    coord = args.centre
    delta = args.delta
    limit = args.limit
    r = np.arange(0, limit, delta)
    
    
    if args.potential:
        pot_coefs = calculate_pot_coefs(args.potential)

    for f in args.files:
        if args.name:
            name = args.name
        else:
            name = f[:-3]
            if args.potential:
                name = name +'_p'
        g = grid.dx2Grid(f)
        if args.potential:
            rdf = g.interpRDF(coord, delta, limit, potential=pot_coefs)
        else:
            rdf = g.interpRDF(coord, delta, limit)
        np.savetxt(name+'.rdf', np.c_[r, rdf])
        if args.plot:
            plt.plot(r, rdf, label=f)
    if args.plot:                
        plt.legend(loc='best')
        #plt.title(os.path.split(os.path.dirname(os.path.realpath(__file__)))[1])
        plt.show()
    
    
if __name__ == '__main__':
    main(sys.argv[1:])

