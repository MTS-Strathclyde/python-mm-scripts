#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 14 09:06:17 2015

@author: max
"""
from __future__ import print_function, division
import argparse
import sys
import numpy as np
from scipy.signal import argrelextrema
import matplotlib.pyplot as plt
from scipy.integrate import simps
import os
from subprocess import Popen, PIPE

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
    parser = argparse.ArgumentParser(description="""Aanlyze xvg file containing
                        Gromacs rdf.""")
    #Positional args
    parser.add_argument('rdf',
                        help="""RDF file in xvg for""")
    #Optional args
    parser.add_argument('-n', '--num_molecules',
                        help="""Number of molecules in box [100].""",
                        default=100, type=int)
    parser.add_argument('-m', '--multiplicity',
                        help="""Multiplicity of atom in molecule [1].""",
                        default=1, type=int)
    parser.add_argument('-e', '--edr_file',
                        help="""GROMACS edr file. Curretly is used to 
                        extract average volume.""")
    parser.add_argument('-v', '--volume',
                        help="""Volume of box in nm [60.0]. Use it if edr is
                        not available.""",
                        default=60.0, type=float)

    return parser.parse_args(argv)


def get_volume(edr_file):
    """Read volume from edr file"""
    with open(os.devnull, 'w') as fnull:
        p = Popen(['gmx', 'energy', '-f', edr_file], stdout=PIPE, 
                  stdin=PIPE, stderr=fnull)
        out = p.communicate(input='Volume')[0]
    vol_line = out.splitlines()[-1]
    return float(vol_line.split()[1])


def load_xvg(fname):
    r = []
    rdf = []
    with open(fname) as f:
        for l in f:
            if not (l.startswith('#') or l.startswith('@')):
                l = l.split()
                r.append(l[0])
                rdf.append(l[1])
    return np.array(r, dtype=float), np.array(rdf, dtype=float)
    
    
def plot_and_pick(r, rdf, loc_min_r, loc_min_value):
    picked_point = []
    def onpick1(event):
        if isinstance(event.artist, plt.Line2D):
            thisline = event.artist
            xdata = thisline.get_xdata()
            ydata = thisline.get_ydata()
            ind = event.ind
            picked_point.append((np.take(xdata, ind), np.take(ydata, ind)))
            print('onpick line:', picked_point[-1])
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.set_title("""If the minimum is found incorrectly pick correct point and close figure. Otherwise
                    just close the figure""")
    ax1.set_xlabel('r [nm]')
    ax1.set_ylabel('g(r)')
    line, = ax1.plot(r, rdf, 'o-', picker=4, label='rdf')
    print(loc_min_r, loc_min_value)
    ax1.plot(loc_min_r, loc_min_value, 'rx', ms=15, label='minimum')
    fig.canvas.mpl_connect('pick_event', onpick1)
    ax1.legend(loc='best')
    plt.show()
    if picked_point:
        return picked_point[-1]
    else:
        return loc_min_r, loc_min_value
        

def find_first_minima_idx(r, rdf):
    """Finds first local minima. The search is relatively unsofisticated.
    Plots the figure to confirm minima value. If the user disagrees with minima
    found asks user to provide estimate.
    """
    # Get first local minima idx and values
    grt_than_one_first_idx = np.where(rdf > 1.0)[0][0]
    cut_rdf = rdf[grt_than_one_first_idx:]
    loc_min_idx = argrelextrema(cut_rdf, np.less, order=3)[0][0]
    loc_min_value = cut_rdf[loc_min_idx]
    loc_min_r = r[np.where(rdf==loc_min_value)[0][0]]
    # Plot the found minimum and confirm it
    #return plot_and_pick(r, rdf, loc_min_r, loc_min_value)
    return loc_min_r, loc_min_value

def main(argv):
    args = process_command_line(argv)
    if args.edr_file:
        volume = get_volume(args.edr_file)
    else:
        volume = args.volume
    r, rdf = load_xvg(args.rdf)
    loc_min_r, loc_min_value = find_first_minima_idx(r, rdf)
    min_idx = np.argmin(abs(r - loc_min_r)) + 1
    dx = r[1] - r[0]
    coord_num = simps(rdf[:min_idx]*4*np.pi*r[:min_idx]**2, dx=dx)
    coord_num = coord_num*args.multiplicity*args.num_molecules/volume
    print()
    print('Coordination number: {}'.format(coord_num))
    print() 
    
if __name__ == '__main__':
    main(sys.argv[1:])

