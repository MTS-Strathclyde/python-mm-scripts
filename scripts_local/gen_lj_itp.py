#!/usr/bin/env python2

import sys
import argparse

#epsilon0 =  8.854187817e-12  #F/m
#N = 6.02214129e23
#k = 1.9872041E-3  # kcal/mol/K
R = 8.3144598E-3  # kJ/mol/K



Mr = 20.0 # amu. Any mass will do, this will only affect dynamics, not energies
# with reduced units we can vary either T, rho or eps, sigma
# we will do the later to simplify batch submissions
T = 298.15 #K
rho_A = 7.5e-3 # 1/A^3

itp_txt = """[ atomtypes ]
; name      at.num  mass     charge ptype  sigma      epsilon
LJs           8      {mass:.2f}     0.0000  A   {sig:.5e}  {eps:.5e}

[ moleculetype ]
; molname    nrexcl
SOL        1

[ atoms ]
; id  at type     res nr  res name  at name  cg nr  charge    mass
  1   LJs          1       SOL       LJs      1       0.0000   {mass:.2f}
"""


def process_command_line(argv):
    """
    Processes arguments and returns namespace of them
    """
    parser = argparse.ArgumentParser(description="""Create itp of LJ solvent
                         for gromacs. """)
    parser.add_argument('--rhor', help="""Reduced density.""",
                        type=float)
    parser.add_argument('--Tr', help="""Reduced temperature.""",
                        type=float)
    parser.add_argument('--sigma', help=""" LJ sigma [A].""",
                        type=float)
    parser.add_argument('--epsilon', help=""" LJ eps [kcal/mol]""",
                        type=float)
    parser.add_argument('--mass', help=""" Mass in Da [20.0]""",
                        type=float, default=20.0)
    return parser.parse_args(argv)


def main(argv):
    args = process_command_line(argv)
    if args.rhor or args.Tr:
        sigma = (args.rhor/rho_A)**(1./3)  # A
        epsilon = R*T/args.Tr  # kJ/mol
    else:
        sigma = args.sigma # A
        epsilon = args.epsilon*4.184  #kJ/mol
    txt = itp_txt.format(sig=sigma/10., eps=epsilon, mass=args.mass)
    for l in txt.splitlines():
        print l

if __name__ == '__main__':
    main(sys.argv[1:])


