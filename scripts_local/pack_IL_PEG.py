#!/usr/bin/env python

# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 10:26:33 2015

@author: max
"""

import argparse
import os
import sys
import shutil
import subprocess


PACK_INP = """
tolerance 2.0
filetype pdb
output {name}.pdb


structure {peg}
  number {number_peg_mol}
  inside box 0. 0. 0. {box_size} {box_size} {box_size}
end structure

structure {C}
  number {number_il_mol}
  inside box 0. 0. 0. {box_size} {box_size} {box_size}
end structure

structure {A}
  number {number_il_mol}
  inside box 0. 0. 0. {box_size} {box_size} {box_size}
end structure

"""

BOXBUILDER = """
#!/bin/bash

for i in `seq 1 {N}`;
do
	mkdir $i
	cp *pdb *inp $i
	cd $i
	packmol.exe < *inp
	yes 0 | gmx trjconv -f {name}.pdb -o {name}.gro -s {name}.pdb -box {box_size} {box_size} {box_size}
	cd ..
	echo $i
done

"""


def process_command_line(argv):
    """
    Processes arguments and returns namespace of them
    """
    parser = argparse.ArgumentParser(description="""Prepare N GROMACS gro files
                        with a given molar fraction of PEG.""")
    #Positional args
    parser.add_argument('peg', metavar='peg.pdb',
                        help="""PEG structure in pdb format.""")
    parser.add_argument('C', metavar='C.pdb',
                        help="""Cation of ionic liquid in pdb format.""")
    parser.add_argument('A', metavar='A.pdb',
                        help="""Anion of ionic liquid pdb format.""")       
    parser.add_argument('fraction', metavar='fraction',
                        help="""Mole fraction of PEG in the mixture. Valid
                        numbers are between 0 and 1. 
                        fraction=(num of peg mols)/(num of peg + num of il pairs)""",
                        type=float)
    #Optional args
    parser.add_argument('-pd', '--peg_density',
                        help=""" Number density of peg bulk phase in 1/nm^3
                        [1.0] """,
                        default=1.0, type=float)
    parser.add_argument('-id', '--il_density',
                        help="""Number density of IL bulk phase in 1/nm^3.
                        An anion and cation pair is 1, not 2.
                        [1.0]""",
                        default=1.0, type=float)
    parser.add_argument('-N', '--number', help="""Number of produced packings
                        with different seed.""", default=5, type=int)
    parser.add_argument('-s', '--size', help="""The length of the side of
                        cubic box side in nm [10].""", default=10, type=float)
    return parser.parse_args(argv)


def main(argv):
    args = process_command_line(argv)
    # Prepare dir
    name = '{}_{}_{}_PegFraction_{}_Size_{}'.format(args.peg[:-4], args.C[:-4],
                                                args.A[:-4], args.fraction,
                                                args.size)
    os.mkdir(name)
    shutil.copy(args.peg, name)
    shutil.copy(args.C, name)
    shutil.copy(args.A, name)
    os.chdir(name)
    # Create inp file for packmol
    # Box size should be in A
    # Number of PEG = V_{tot} * [ 1/rho_{peg} + 1/rho_{il}*(x^{-1} - 1) ]^{-1}
    # Number of IL = #PEG*(1/x - 1)
    # here x = molar fraction of PEG
    number_peg_mol = args.size**3 * ( 1/args.peg_density + \
                     1/args.il_density*(args.fraction**(-1) - 1))**(-1)
    number_il_mol = number_peg_mol * (args.fraction**(-1) - 1)
    number_peg_mol = int(round(number_peg_mol))
    number_il_mol = int(round(number_il_mol))
    pack_inp = PACK_INP.format(name=name, peg=args.peg, C=args.C, A=args.A,
                    box_size=args.size*10, number_peg_mol=number_peg_mol,
                    number_il_mol=number_il_mol)
    with open('pack.inp', 'w') as f:
        f.write(pack_inp)
    # Write boxbuilder script and run it
    # Increment box size by 0.1 nm to be sure nothing is left out
    # Box size is in nm
    box_size = args.size + 0.1
    boxbuilder = BOXBUILDER.format(N=args.number, name=name,
                                   box_size=box_size)
    with open('boxbuilder.sh', 'w') as f:
        f.write(boxbuilder)
    subprocess.call(['bash', 'boxbuilder.sh'])


if __name__ == '__main__':
    main(sys.argv[1:])

