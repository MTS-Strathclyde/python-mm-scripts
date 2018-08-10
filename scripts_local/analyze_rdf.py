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

K_B = 1.9872041E-3 # boltzmann const in kcal/mol/K
N_A = 6.022141e23 # avogadro's constant


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
    parser = argparse.ArgumentParser(description="""Aanlyze gvv file containing
                        (RISM) rdfs.""")
    #Positional args
    parser.add_argument('rdf',
                        help="""RDF file in gvv format""")
    parser.add_argument('-x', '--xvv',
                        help=""" Xvv file to get site densities""")
    parser.add_argument('-d', '--densities',
                        help=""" List numeric site densities in the same
                        order as they are appear in rism1d input file.
                        Units are: 1/A^3.""",
                        nargs='+', type=float)
    parser.add_argument('--plot',
                        help=""" Plot rdfs and positions of the first
                        minima/maxima""",
                        action='store_true')
    parser.add_argument('-c', '--columns',
                        help="""Specify number(s) of column in vv file to be
                        analyzed. Numbering starts from 1.""", nargs='+',
                        type=int)
    parser.add_argument('--h_of_r',
                        help=""" File contains total correlation functions.
                        Adds 1 to all points.""",
                        action='store_true')
    parser.add_argument('-t', '--threshold',
                        help=""" How high should peak be to be considered a
                        peak [1.0].""",
                        type=float, default=1.0)
                        
    return parser.parse_args(argv)


class Xvv(object):
    """ Wrapper around xvvfile used to compute 3d-rism pressure """
    def __init__(self, fname):
        """ Read xvvfile and set instance attributes 
        
        Parameters
        ----------
        
        fname : string
            Path to a valid xvv file
        """
        self.fname = fname
        self.ngrid = None
        self.nsites = None
        self.nspecies = None
        self.temperature = None
        self.dr = None
        self.atom_names = None
        self.densities = None
        self.xvv_data = None
        self.multiplicities = None
        self.unique_sites_per_species = None
        self.total_sites_per_species = None
        self.species_densities = None
        self.normalized_densities = None
        self._read_xvvfile()
        self._compute_species_properties()

    def _read_xvvfile(self):
        with open(self.fname) as f:
            lines = f.readlines()
        tot_lines = len(lines)
        for i, line in enumerate(lines):
            line = line.split()
            if len(line) <= 1:
                continue
            if line[1] == 'POINTERS':
                data = map(int, lines[i+2].split())
                self.ngrid, self.nsites, self.nspecies = data
            if line[1] == 'MTV':
                self.multiplicities = map(int, lines[i+2].split())
            if line[1] == 'NVSP':
                self.unique_sites_per_species = map(int, lines[i+2].split())
            if line[1] == 'THERMO':
                data = lines[i+2].split()
                self.temperature = float(data[0]) # K
                self.dr = float(data[4]) # Angstrom
            if line[1] == 'ATOM_NAME':
                data = lines[i+2].strip()
                #split into groups of 4
                self.atom_names = [data[i:i+4].strip() for i in range(0, len(data), 4)]
            if line[1] == 'RHOV' and len(line) == 2:
                self.densities = map(float, lines[i+2].split())
                #are there more lines with density?
                counter = 3
                while lines[i+counter].startswith(' '):
                    self.densities.extend(map(float, lines[i+counter].split()))
                    counter += 1
                try:
                    assert len(self.densities) == len(self.atom_names)
                except AssertionError:
                    print('Inconsistent number of densities and atom names')
                    print(self.densities)
                    print(self.atom_names)
                    raise ValueError
            if line[1] == 'XVV' and len(line) == 2:
                self.xvv_data = []
                xvv_ind = i + 2
                while xvv_ind < tot_lines and not lines[xvv_ind].startswith('%'):
                    self.xvv_data.extend(lines[xvv_ind].split())
                    xvv_ind += 1
                break
        assert len(self.xvv_data) == self.ngrid*self.nsites*self.nsites
        self.xvv_data = np.array(self.xvv_data, dtype=float)
        self.xvv_data = np.reshape(self.xvv_data,
                                   (self.ngrid, self.nsites, self.nsites),
                                   order='F')

    def _compute_species_properties(self):
        self.normalized_densities = []
        for density, multiplicity in zip(self.densities, self.multiplicities):
            self.normalized_densities.append(density/multiplicity)
        self.species_densities = []
        self.total_sites_per_species = []
        pointer = 0 
        for sp_sites in self.unique_sites_per_species:
            pointer += sp_sites
            total_sites = sum(self.multiplicities[pointer - sp_sites:pointer])
            self.total_sites_per_species.append(total_sites)
            self.species_densities.append(self.normalized_densities[pointer - 1])
        assert len(self.species_densities) == self.nspecies
    
    def compute_3drism_pressures(self, k=0):
        """ Compute 3drism pressure using loaded xvv file.
        Uses equation 20 from the article by Sergiievskyi et al. 
        (http://dx.doi.org/10.1063/1.4935065). 

        Parameters
        ----------
        k : int
            Which k value to use to compute pressure. The pressure can be pretty
            sensitive to it. It is recommended to experiment with a couple of
            k values or better, plot dependency of pressure on it to see
            which value works best.
            
        Return
        ------
        pressures : tuple of floats
            Tuple containeing two pressures.
            First element is 3D-RISM pressure (used in PC), second element is
            3D-RISM pressure minus ideal gas pressure (used in PC+).
            Both have units of kcal/mol/A^3.
        """
        xvv_k = self.xvv_data[k,:,:]
        density_vec = np.array(self.normalized_densities)
        mult_vec = np.array(self.multiplicities)
        # Z_k from sergievskyi's article
        z_k = mult_vec/density_vec*(np.identity(self.nsites) - np.linalg.inv(xvv_k))
        z_k_sum_densities2 = np.sum(density_vec*z_k*density_vec.T)
        densities_times_sites = [sites*dens for sites, dens in zip(self.total_sites_per_species,
                                                                   self.species_densities)]
        pressure = sum(densities_times_sites) - .5*z_k_sum_densities2
        pressure = pressure*self.temperature*K_B
        ideal_pressure  = sum(self.species_densities)*K_B*self.temperature
        return pressure, pressure - ideal_pressure



def load_gvv(fname):
    with open(fname) as f:
        lines = f.readlines()
    if lines[0].startswith('#'):
        #amber rism
        labels = lines[3][1:].split()[1:]
        data = [l.split() for l in lines[4:]]
    else:
        # assume it is in numpy format
        data = [l.split() for l in lines]
        labels = [str(i+1) for i in range(len(data[0]))]
    data = np.array(data, dtype=float)
    return data, labels
    

def plot_and_pick(r, rdf, loc_max_r, loc_max_value, loc_min_r, loc_min_value):
#    picked_point = []
#    def onpick1(event):
#        if isinstance(event.artist, plt.Line2D):
#            thisline = event.artist
#            xdata = thisline.get_xdata()
#            ydata = thisline.get_ydata()
#            ind = event.ind
#            picked_point.append((np.take(xdata, ind), np.take(ydata, ind)))
#            print('onpick line:', picked_point[-1])
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
#    ax1.set_title("""If the minimum or maximum are found incorrectly,
#pick correct points and close the figure. Otherwise just close the figure""")
    ax1.set_xlim((0, 15))
    ax1.set_title('Check that everything is ok')
    ax1.set_xlabel('r [A]')
    ax1.set_ylabel('g(r)')
    line, = ax1.plot(r, rdf, '-', picker=4, label='rdf')
    #print('maximum:',loc_max_r, loc_max_value)
    #print('minimum:',loc_min_r, loc_min_value)
    ax1.plot(loc_min_r, loc_min_value, 'bo', ms=15, label='minimum')
    ax1.plot(loc_max_r, loc_max_value, 'ro', ms=15, label='maximum')
    #fig.canvas.mpl_connect('pick_event', onpick1)
    ax1.legend(loc='best')
    plt.show()
    #print(picked_point)
#    if picked_point:
#        # which is which?
#        p1, p2 = picked_point[-1], picked_point[-2]
#        if p1[1] >= p2[1]:
#            return p1, p2
#        else:
#            return p2, p1
#    else:
#        return (loc_max_r, loc_max_value), (loc_min_r, loc_min_value)
    return (loc_max_r, loc_max_value), (loc_min_r, loc_min_value)


def find_first_extrema(r, rdf, threshold, plot=False):
    """Finds first local minima and maxima and returns 2 tuples with maximum
    and mimum coordinates. The search is relatively unsofisticated.
    Can plot results if needed.
    """
    # cut the beginning to avoid strange results
    grt_than_threshold_first_idx = np.where(rdf > threshold)[0][0]
    cut_rdf = rdf[grt_than_threshold_first_idx:]
    # find extrema
    loc_max_idx = argrelextrema(cut_rdf, np.greater, order=3)[0][0]
    loc_max_value = cut_rdf[loc_max_idx]
    loc_max_r = r[np.where(rdf==loc_max_value)[0][0]]

    loc_min_idx = argrelextrema(cut_rdf, np.less, order=3)[0][0]
    loc_min_value = cut_rdf[loc_min_idx]
    loc_min_r = r[np.where(rdf==loc_min_value)[0][0]]
    
    if plot:
        return plot_and_pick(r, rdf, loc_max_r, loc_max_value, 
                             loc_min_r, loc_min_value)
    else:
        return (loc_max_r, loc_max_value), (loc_min_r, loc_min_value)
    

def compute_site_site_densities(site_densities):
    ss_densities = []
    # essentiall, create upper diagonal matrix out of site densities vector
    # ordered in columns first order
    for i, _ in enumerate(site_densities):
        j = 0
        while j <= i:
            ss_densities.append(site_densities[j])
            j += 1
    return ss_densities


def main(argv):
    args = process_command_line(argv)
    data, labels = load_gvv(args.rdf)
    r, rdfs = data[:,0], data[:,1:]
    if args.h_of_r:
        rdfs = rdfs + 1
    dx = r[1] - r[0]
    # compute densities
    if args.xvv:
        xvv_obj = Xvv(args.xvv)
        site_densities = xvv_obj.densities
    elif args.densities:
        site_densities = args.densities
    else:
        print("No densities supplied. The coordination numbers will be incorrect \
        and will have to be multiplied by at least site multiplicities to provide \
        any meaningful information")
        site_densities = [1 for i in labels]
    ss_densities = compute_site_site_densities(site_densities)
#    print(site_densities)
#    print (ss_densities)
    # check what we need to analyze
    if args.columns:
        columns_to_analyze = map(lambda x: x - 1, args.columns)
    else:
        columns_to_analyze = [i for i in range(len(labels))]
    for i, (label, rdf, density) in enumerate(zip(labels, rdfs.T, ss_densities)):   
        if i in columns_to_analyze:
            max_p, min_p = find_first_extrema(r, rdf, args.threshold, args.plot)
            min_idx = np.argmin(abs(r - min_p[0])) + 1
            coord_num = simps(rdf[:min_idx]*4*np.pi*r[:min_idx]**2, dx=dx)
            coord_num = coord_num*density
            print('{:28} {:^6} {:^6}'.format(label, 'A', 'g(r)'))
            print('location of the first peak: {p[0]: .3f} {p[1]: .3f}'.format(p=max_p))
            print('location of the first min:  {p[0]: .3f} {p[1]: .3f}'.format(p=min_p))
            print('Coordination number:        {: .3f}'.format(coord_num))
            print() 
    
if __name__ == '__main__':
    main(sys.argv[1:])

