#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 17 14:04:42 2015

@author: max
"""

from __future__ import print_function

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
    parser = argparse.ArgumentParser(description="""Run 3D-RISM single point
            calculation at given temperature. Uses GAFF force field and
            AM1-BCC charges.""")
    #Positional args
    parser.add_argument('pdb', metavar='molec.pdb',
                        help="""Input file. Must be in pdb format
                        acceptable by Antechamber. Must end with .pdb""")
    parser.add_argument('prmtop', metavar='molec.prmtop',
                        help=""" Input file. Must be in prmtop format
                        acceptable by Antechamber.""")
    #Optional args
#    parser.add_argument('-c', '--molcharge',
#                        help="""Charge of the solute [0]""", default=0,
#                        type=int)
#    parser.add_argument('--multiplicity',
#                        help="""Multiplicity of the solute [1]""", default=1,
#                        type=int)
    return parser.parse_args(argv)


def get_sigma2_radii(prmtop):
    """ Return a list of sigma/2 radii """
    try:
        parm = AmberParm(prmtop)
        radii = []
        for atom in parm.atom_list:
            nbidx = parm.LJ_types[atom.attype]
            radii.append(float(parm.LJ_radius[nbidx - 1])/(2**(1./6)))
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
                    radii = [(acoef/bcoef)**(1./6)/2]
    return radii


def get_charges(prmtop):
    charges = []
    with open(prmtop) as f:
        for line in f:
            if line.startswith('%FLAG CHARGE'):
                f.next() # skip one line
                next_line = next(f)
                while next_line.startswith(' '):
                    l_chrgs = next_line.split()
                    charges.extend(map(lambda x: float(x)/18.2223, l_chrgs))
                    next_line = next(f)
    return charges


def get_atnames_resnames(pdb):
    """ return atom names and resnames from pdb file """
    names = []
    with open(pdb) as f:
        for line in f:
            if line.startswith('ATOM'):
                names.append(line[13:20].split())
    return names


def main(argv):
    args = process_command_line(argv)
    fname = args.pdb[:-4]
    names = get_atnames_resnames(args.pdb)
    charges = get_charges(args.prmtop)
    radii = get_sigma2_radii(args.prmtop)
    # write charges
    format_string = '{:<4}  {:<8} {: .3f}\n'
    with open(fname + '.crg', 'w') as f:
        f.write('atom__resnumbc_charge_\n')
        for (atom, res), chg in zip(names, charges):
            f.write(format_string.format(atom, res, chg))
    # write radii
    with open(fname + '.siz', 'w') as f:
        f.write('atom__resnumbc_radius_\n')
        for (atom, res), r in zip(names, radii):
            f.write(format_string.format(atom, res, r))



if __name__ == '__main__':
    main(sys.argv[1:])

