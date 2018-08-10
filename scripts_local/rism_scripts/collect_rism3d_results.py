#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 26 12:01:13 2015

@author: max

TODO!
Built-in join with experimental data


"""



from __future__ import print_function

import argparse
import sys
import os
import pandas as pd
import numpy as np
import glob
from collections import OrderedDict
from scipy.spatial import cKDTree


SHORT = OrderedDict([('Name' ,               None,),
         ('Xvv' ,                'mm_options:  xvvfile=',),
         ('Closure' ,            'mm_options:  closure=',),
         ('rism_exchem' ,        'rism_exchem',),
         ('PMV' ,                'rism_volume',),
         ('ISc' ,                'dGhyd(ISc)='),
         ])

DEFAULT = OrderedDict([('Name' ,               None,),
           ('Xvv' ,                '\tmm_options:  xvvfile=',),
           ('Closure' ,            '\tmm_options:  closure=',),
           ('rism_exchem' ,        'rism_exchem',),
           ('GF' ,                 'rism_exchGF',),
           ('UV' ,                 'rism_potUV',),
           ('PMV' ,                'rism_volume',),
           ('PC+' ,                'dGsolv(PC+)=',),
           ('PC' ,                 'dGsolv(PC)=',),
           ('P+' ,                 'P_minus_ideal_gas_pressure=',),
           ('P' ,                  'P=',),
           ('Time (s)' ,           '3D-RISM runtime: ',),
           ('KBI' ,            'rism_KB',),
           ('DCFI' ,           'rism_DCFI',)
           ])
               
ALL =  OrderedDict([('Name' ,                 None,                     ),
        (   'Xvv' ,                '\tmm_options:  xvvfile=',),
        (  'Closure' ,             '\tmm_options:  closure=',),
        (   'Grdx' ,               '\tmm_options:  grdspcx=',),
        (   'Grdy' ,               '\tmm_options:  grdspcy=',),
        (   'Grdz' ,               '\tmm_options:  grdspcz=',),
        (   'Tolerance' ,          '\tmm_options:  tolerance',),
        (   'Buffer' ,             '\tmm_options:  buffer=',),
        (   'rism_exchem' ,        'rism_exchem',),
        (   'GF' ,                 'rism_exchGF',),
        (   'UV' ,                 'rism_potUV',),
        (   'PMV' ,                'rism_volume',),
        (   'ExNumb' ,             'rism_exNumb',),
        (   'ExChg' ,              'rism_exChrg',),
        (   'UC' ,                 None,),
        (   'ISc' ,                'dGhyd(ISc)=',),
        (   'Time (s)' ,           '3D-RISM runtime: ',),
        (   'ConvgSteps' ,         '|RXRISM converged in')
        ])
           

PC_RE = OrderedDict([('pc', None),
                     ('pcp', None)])
           
           
POLAR = OrderedDict([ ('polar_exchem' ,            'rism_polar',),
          ('apolar_exchem' ,           'rism_apolar',),
          ('polar_GF' ,                'rism_polGF',),
          ('apolar_GF' ,               'rism_apolGF',),
          ('polar_UV' ,                'rism_pol_potUV', ),
          ('apolar_UV' ,               'rism_apol_potUV'),
        ])

FUNCTIONAL = OrderedDict([ #('ideal' ,            'F_ideal:',),
          #('potential' ,               'F_ext:',),
          #('excess' ,                  'F_exc:',),
            ('entrop_vv',  'F_exc:'),
#          ('mod1' ,                  'F_mod:',),
#          ('mod2' ,                  'F_mod2:',),
        ])

INTS = OrderedDict([ ('Ideal_log',       ' Ideal integral (log part): '),
                     ('Ideal_tot',       ' Ideal integral tot:')])



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
    parser = argparse.ArgumentParser(description="""Collect results of 3D-RISM
                        calculation in csv file.""")
    #Positional args
    parser.add_argument('root_directory',
                        help="""Reads all files under this directory recursively
                        trying to identify calculation directories. (containing
                        .log files).""")
    #Optional args
    parser.add_argument('-n', '--name',
                        help="""Name of output csv file with the results 
                        [results.csv].""",
                        default='results.csv')
    parser.add_argument( '--full_name',
                        help="""Name of the molecule will also contain path.""",
                        action='store_true')
    parser.add_argument( '--dir_name',
                        help="""Use name of calculation directory as a molecules name.""",
                        action='store_true')
    parser.add_argument('-o', '--output',
                        help="""Level of details in the output. 0 - least,
                        2 - most verbose. [1]""",
                        default='1',choices=['0', '1', '2'])
    parser.add_argument('--atom_decomp',
                        help="""All thermodynamic parameters will be additionally
                        decomposed by site contributions""",
                        action='store_true')
    parser.add_argument('--polar_decomp',
                        help="""Some thermodynamic parameters will be additionally
                        split into polar/nonpolar contributions""",
                        action='store_true')
    parser.add_argument('--functional_decomp',
                        help=""" Save results of functional decomposition (ideal,
                        potential, excess). """,
                        action='store_true')
    parser.add_argument('--integrals',
                        help=""" Save KB and direct cor. function integrals. """,
                        action='store_true')
    parser.add_argument('--pc_recalc',
                        help=""" Calculate proper pc correction using xvv file. """,
                        action='store_true')
    parser.add_argument('--old_results',
                        help=""" results.txt contains dG(hyd) instead of
                        dG(solv). """,
                        action='store_true')
    parser.add_argument('--new_results',
                        help=""" Results obtaind with Ambertools 16. """,
                        action='store_true')
    parser.add_argument('--extra_output_name',
                        help="""The name of the file containing extra output. [results.txt]""",
                        default='results.txt')
                        
    return parser.parse_args(argv)
    
    
K_B = 1.9872041E-3 # boltzmann const in kcal/mol/K

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
        # Unique sites per species
        self.unique_sites_per_species = None
        # Actual number of sites per species
        self.total_sites_per_species = None
        # Density of each species
        self.species_densities = None
        # Sites from the same species have densities equal to densities of species
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

    def get_k(self):
        """ Return array of k values"""
        dk = np.pi/(self.ngrid*self.dr)
        return np.array([i*dk for i in range(self.ngrid)])


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
        self.mol_densities = []
        for i, sp_sites in enumerate(self.unique_sites_per_species):
            self.mol_densities.append(self.species_densities[i])
            if sp_sites > 1:
                self.mol_densities.extend([0,]*(sp_sites-1))
        assert len(self.mol_densities) == self.nsites


    def _get_compressibility_from_therm(self, therm_p):
        """ Return solvent compressibility.
    
        Parameters
        ----------
        therm_p : string
            path to .therm file
    
        Returns
        -------
        compres : float
            Units: 1/MPa
        """
        with open(therm_p, 'r') as f:
            therm_lines = f.readlines()
        compres = float(therm_lines[2].split()[-1])
        units = therm_lines[2].split()[-2]
        if units == '[10e-4/MPa]':
            # Old version of Ambertools
            return compres*10e-4
        if units == '[1/kPa]':
            # # !! This is ambertools 14, where compressiblity is bugged.
            # Units are shown to be [1/kPa], while in reality compres thea are in [1/MPa]
            # http://archive.ambermd.org/201503/0651.html
            return compres
        if units == '[1/MPa]':
            # This is ambertools 15
            # All is good
            return compres
        else:
            raise ValueError('Unknown compressiblity format, check *.therm file')


    def compute_3drism_pressures(self, k=0):
        """ Compute 3drism pressure using loaded xvv file.
        Uses equation 20 from the article by Sergiievskyi et al. 
        (http://dx.doi.org/10.1063/1.4935065). 

        Parameters
        ----------
        k : int
            Which k value to use to compute pressure. 1 is recommended.
            
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
        #print self.total_sites_per_species
        ideal_pressure  = sum(self.species_densities)*K_B*self.temperature
        #print 'under_pressure',pressure - 2*ideal_pressure
        return pressure - ideal_pressure, pressure
        

    def proper_pressure(self, k=0):
        kT = self.temperature*K_B
        xvv_k = self.xvv_data[k,:,:]
        density_vec = np.array(self.densities)
        mult_vec = np.array(self.multiplicities)
        # C_k0_rho is a vector of full site-site direct correlation functions
        # summed up column wise
        C_k0_rho = np.identity(self.nsites) - np.linalg.inv(xvv_k)
        C_k0_rho = np.sum(C_k0_rho, axis=0)
        # pressure_mult - the multiplier of sum rho_i V_i
        pressure_mult = kT*(1 - C_k0_rho/2)
        ideal_pressure  = np.array(self.mol_densities)*kT
        return density_vec, ideal_pressure, pressure_mult


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


def is_calc_dir(files):
    """Checks whether folder is the source of calculation."""
    has_log = False
    has_prmtop = False
    for f in files:
        if f.endswith('.log'):
            has_log = True
        if f.endswith('.prmtop'):
            has_prmtop = True
    return has_log and has_prmtop
        
        
def prepare_dictionary(args):
    #print(args.output)
    if args.output == '0':
        dic = SHORT
    elif args.output == '1':
        dic = DEFAULT
    elif args.output == '2':
        dic = ALL
    else:
        raise ValueError('Incorrect args.output value')
    if args.pc_recalc:
        dic.update(PC_RE)
    if args.polar_decomp:
        dic.update(POLAR)
    if args.functional_decomp:
        dic.update(FUNCTIONAL)
    if args.integrals:
        dic.update(INTS)
    return dic
    

def collect_calc_data(path, files, dic, args):
    """Analyzes calculation directory and return its results."""
    results_dic = [(k , None) for k in dic.iterkeys()]
    results_dic = OrderedDict(results_dic)
    log_file = None
    # First figure out the name
    for f in files:
        if f.endswith('.prmtop'):
            name = f[:-7]
            if args.full_name:
                results_dic['Name'] = os.path.join(path, name)
            elif args.dir_name:
                results_dic['Name'] = os.path.basename(path)
            else:
                results_dic['Name'] = name
            # add files to namespace
#!HACK            
            args.pdb = os.path.join(path, name + '.pdb')
            args.prmtop = os.path.join(path, name + '.prmtop')            
        if f.endswith('.log'):
            #assume this is the log file, but throw error if more than one
            #are found
            if not log_file:
                log_file = os.path.join(path, f)
                with open(log_file) as logf:
                    log_lines = logf.readlines()
            else:
                #raise ValueError('More than one log file found in {}'.format(path))
                print('WARNING: Multiple log files found in {}'.format(path))
    try:  #simply append result.txt to log lines
        with open(os.path.join(path, args.extra_output_name)) as f:
            log_lines.extend(f.readlines())
    except IOError, e:
        if e.errno == 2:
            print(e)
        else:
            raise e
    if args.functional_decomp:
        # see if decomp.txt is present
        # if not - no biggy
        try:
            with open(os.path.join(path, 'decomp.out')) as f:
                log_lines.extend(f.readlines())
        except IOError, e:
            pass
    # deal with easy values
    for title, value in dic.iteritems():  # iterate over titles and startwith strings
        if value:
            for l in log_lines:
#                if 'xvvfile' in l:
#                    print(l)
                if l.startswith(value): #found the string
                    l = l.replace(value, '')
                    strings = l.split()  # all the numbers and units and stuff
                    #if title == 'Xvv':
                    #    strings[0] = os.path.split(strings[0])[1]
                    if title in ['ExNumb', 'ExChg']:
                        results_dic[title] = ', '.join(strings)
                    # each number gets a separate name
                    if title in ['KBI', 'DCFI'] + INTS.keys():
                        results_dic.pop(title)
                        for i, value in enumerate(strings):
                            name = title + '_s{}'.format(i)
                            results_dic[name] = value.strip()
                    else:
                        results_dic[title] = strings[0].strip()
    if args.pc_recalc:
        # computing proper pc correction
        try:
            xvv_path = os.path.join(path, results_dic['Xvv'])
            xvv_inst = Xvv(xvv_path)            
            kb_vect = np.array([float(results_dic['KBI_s{}'.format(i)]) for i in range(xvv_inst.nsites)])
            density_vec, ideal_pressure, pressure_mult = xvv_inst.proper_pressure()
            pv = np.sum(density_vec*kb_vect*pressure_mult)
            ppv = pv - np.sum(kb_vect*ideal_pressure)
            results_dic['pc'] = float(results_dic['rism_exchem'])  + pv
            results_dic['pcp'] = float(results_dic['rism_exchem'])  + ppv
        except KeyError, e:
            print('Key: ',e,'was not found. Not computing PC correction')
        except IOError, e:
            print('Xvv file: ',e,'was not found. Not computing PC correction')
    if 'UC' in dic:
        try:
            results_dic['UC'] = float(results_dic['GF']) - 3.2217 * 0.0333 * float(results_dic['PMV']) + 0.5783
        except KeyError:
            pass
    # fix xvv name
    results_dic['Xvv'] = os.path.split(results_dic['Xvv'])[1]
    return results_dic


def main(argv):
    args = process_command_line(argv)
    if args.old_results:
        DEFAULT['PC+'] = 'dGhyd(PC+)='
        DEFAULT['PC'] = 'dGhyd(PC)='        
    if args.new_results:
        DEFAULT['rism_exchem'] = 'rism_excessChemicalPotential'
        DEFAULT['UV'] = 'rism_solventPotentialEnergy'
        DEFAULT['PMV'] = 'rism_partialMolarVolume'
        DEFAULT['KBI'] = 'rism_KirkwoodBuff'
        DEFAULT['DCFI'] = 'rism_DCFintegral'
        DEFAULT.pop('GF')
    dic = prepare_dictionary(args)
    data = []        
    for subdir, dirs, files in os.walk(args.root_directory):
        if is_calc_dir(files):
            print('Analyzing files in: {}'.format(subdir))
            data.append(collect_calc_data(subdir, files, dic, args))
    columns = max([i.keys() for i in data], key=len)
    df = pd.DataFrame(data, columns=columns)
    df.sort(columns='Name', inplace=True)
    df.to_csv(args.name, index=False)


if __name__ == '__main__':
    main(sys.argv[1:])
