#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 21 16:49:38 2015

@author: max
"""

import sys
import argparse
import os
import subprocess


RISM1D_SCRPT = """#!/bin/csh -f

cat > {name}.inp <<EOF
&PARAMETERS
   THEORY='{rism1d}', 
   CLOSURE='{closure}',
   NR={npoints}, DR={grdsp},                  !Grid Size
   OUTLIST='XGTSC', SELFTEST=1,            !Output
   KSAVE={ksave},
   MDIIS_NVEC={mdiis_nvec}, MDIIS_DEL={mdiis_step}, !MDIIS                
   TOLERANCE={tolerance},
   MAXSTEP={npoints},                         !Iterations
   SMEAR=1, ADBCOR={adbcor},                   !Electrostatics
   DIEps={diel},                          !Solvent description
   NSP=1
/
	&SPECIES                               !solvent
   UNITS='{units}',
	DENSITY={conc}d0,
	MODEL="{path}"
/
EOF

rism1d {name} > {name}.out || goto error

"""

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
    parser = argparse.ArgumentParser(description="""Run 1D-RISM single point
            calculation.""")
    #Positional args
    parser.add_argument('file', metavar='solvent.mdl',
                        help="""Solvent mdl file.""")
    parser.add_argument('conc',
                        help="""Concentration of solvent species. Units
                        are set up with separate option (default = g/cm^3).""",
                        type=float)
    parser.add_argument('diel', 
                        help="""Dielectric constant of solvent.""",
                        type=float)
    #Optional args
    parser.add_argument('--closure',
                        help="""Brdige closure which will be used in both
                        1D-RISM and 3D-RISM simmulations. Either HNC, PSEn or
                        KH (n in PSEn should be an integer). [HNC]""",
                        default="HNC")
    parser.add_argument('--units',
                        help="""Units of the concentration [g/cm^3]""",
                        default="g/cm^3")                        
    parser.add_argument('--rism1d',
                        help="""Type of 1D-RISM theory. Only DRISM has been
                        extensively tested [DRISM]""",
                        default="DRISM")                        
    parser.add_argument('-t', '--temperature',
                        help="""Temperature in K at which simulation will be
                        run [298.15]""", default=298.15, type=float)
    parser.add_argument('--tolerance',
                        help=""" Tollerance of MDIIS algorithm in  1D-RISM
                        [1E-12]""",
                        default=1.0e-12, type=float)
    parser.add_argument('--grdsp',
                        help="""Linear grid spacings [0.025]""",
                        default=0.025, type=float)
    parser.add_argument('--mdiis_step',
                        help="""Step for MDIIS algorithm [0.3]""",
                        default=0.3, type=float)
    parser.add_argument('--mdiis_nvec',
                        help="""Number of vectors for MDIIS algorithm [20]""",
                        default=20, type=int)
    parser.add_argument('--adbcor',
                        help="""Numeric parameter for RISM [0.5]""",
                        default=0.5, type=float)
                        
    parser.add_argument('--ksave',
                        help="""Save every N iteration steps [-1]""",
                        default=-1, type=int)
    parser.add_argument('--nsteps',
                        help="""Maximum number of iterations [10000] .""",
                        default=10000, type=int)
    parser.add_argument('--npoints',
                        help="""Number of points [16384] .""",
                        default=16384, type=int)                        
    return parser.parse_args(argv)



def main(argv):
    args = process_command_line(argv)
    _, name = os.path.split(args.file)
    name = name[:-4] #file should end with .mdl
    script_name = name + '_1drism.sh'
    rism1d_srcirpt = RISM1D_SCRPT.format(name=name, rism1d=args.rism1d,
                                       closure=args.closure,
                                       grdsp=args.grdsp,
                                       tolerance=args.tolerance,
                                       temp=args.temperature,
                                       diel=args.diel,
                                       units=args.units,
                                       conc=args.conc,
                                       path=args.file,
                                       mdiis_step=args.mdiis_step,
                                       mdiis_nvec=args.mdiis_nvec,
                                       adbcor=args.adbcor,
                                       ksave=args.ksave,
                                       nsteps=args.nsteps,
                                       npoints=args.npoints
                                       
                                       )
    with open(script_name, 'wb') as f:
        f.write(rism1d_srcirpt)
    print(subprocess.check_output(['bash', script_name]))


if __name__ == '__main__':
    main(sys.argv[1:])