#!/usr/bin/env python

# -*- coding: utf-8 -*-
"""
Created on Fri Sep  5 13:45:03 2014

@author: max
"""
import sys
import subprocess
import shutil
import argparse
import os
import glob
#import subprocess as sub
#import threading
#import signal


##############################################################################
                              # Machine specific defaults #
##############################################################################

GMX_ARCHIE = ['mpirun', '-np', '1', 'gmx_mpi',]  # For 1 core stuff

JOB_SCRIPT_HEADER_ARCHIE = r"""
export PROCS_ON_EACH_NODE=12

# ************* SGE qsub options ****************
#Export env variables and keep current working directory
#$ -V -cwd
#Select parallel environment and number of parallel queue slots (nodes)
#$ -pe mpi-verbose {nnodes}
#Combine STDOUT/STDERR
#$ -j y
#Specify output file
#$ -o out.$JOB_ID
#Request resource reservation (reserve slots on each scheduler run until enough have been gathered to run the job
#$ -R y
#$ -l h_rt=24:00:00
#$ -P fedorov-iles.prj


# ************** END SGE qsub options ************

# supress dump
export GMX_SUPPRESS_DUMP=Y

export NCORES=`expr $PROCS_ON_EACH_NODE \* {nnodes}`

export OMPI_MCA_btl=openib,self

# Start MDRUN in parallel
"""


GMX_FERRARI = ['gmx',]  # For 1 core stuff

JOB_SCRIPT_HEADER_FERRARI = r"""
# supress dump
export NCORES=`expr 12 \* {nnodes}`
export GMX_SUPPRESS_DUMP=Y
# start jobs
"""



##############################################################################
                                #   MDP Options  #
##############################################################################

# topology
# we want lorentz-bertholet rules so we include amber force field
#INCLUDED_FF = "amber03.ff/ffnonbonded.itp"
# for oplsaa rules include oplaa.ff/forcefield.itp

# Free energy stuff
NSTDHDL = 100 # frequency of writing to free energy file

# Run control
EQUILIBRATION_N_STEPS = 100000 # 200
#MD_N_STEPS = 5000000    # 5 ns / timestep in fs

# Neighbour list update
NSTLIST_EM = 1 # for energy minimization
#NSTLIST = 40 # for everything else


LAMBDA_ALL = ("""
; init_lambda_state        0    1    2    3    4    5    6    7    8    9    10   11   12   13   14   15   16   17   18   19   20   21   22   23   24   25
vdw_lambdas              = 0.00 0.00 0.00 0.00 0.00 0.00 0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50 0.55 0.60 0.65 0.70 0.75 0.80 0.85 0.90 0.95 1.00
coul_lambdas             = 0.00 0.20 0.40 0.60 0.80 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00
""",
25)


LAMBDA_Q = ("""
; init_lambda_state        0    1    2    3    4    5
coul_lambdas             = 0.00 0.20 0.40 0.60 0.80 1.00
""",
5)

LAMBDA_VDW = ("""
; init_lambda_state        0    1    2    3    4    5    6    7    8    9    10   11   12   13   14   15   16   17   18   19   20
vdw_lambdas              = 0.00 0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50 0.55 0.60 0.65 0.70 0.75 0.80 0.85 0.90 0.95 1.00
""",
20)

LAMBDA_ALL_S = ("""
; init_lambda_state        0    1    2    3    4    5    6    7    8    9    10   11   12   13   14   15   16   17   18   19 
vdw_lambdas              = 0.0  0.00 0.0  0.00 0.0  0.05 0.1  0.2  0.3  0.4  0.5  0.6 0.65  0.7  0.75 0.8  0.85 0.9  0.95 1.0
coul_lambdas             = 0.0  0.25 0.5  0.75 1.0  1.00 1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
""",
19)


LAMBDA_Q_S = ("""
; init_lambda_state        0    1    2    3    4
coul_lambdas             = 0.0  0.25 0.5  0.75 1.0
""",
4)

LAMBDA_Q_LONG = ("""
; init_lambda_state        0   1   2   3   4   5   6   7   8   9   10
coul_lambdas             = 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0
""",
10)

LAMBDA_VDW_S = ("""
; init_lambda_state        0    1    2    3    4    5    6    7    8    9    10   11   12   13   14   15
vdw_lambdas              = 0.0  0.05 0.1  0.2  0.3  0.4  0.5  0.6 0.65  0.7  0.75 0.8  0.85 0.9  0.95 1.0
""",
15)

FREE_ENERGY = """ ; Free energy control
free_energy              = yes
init_lambda_state        = {lambda_state}
delta_lambda             = 0
calc_lambda_neighbors    = {calc_lambda_neighbors}        ; -1 for All lambdas
; Vectors of lambda specified here
; Each combination is an index that is retrieved from init_lambda_state for each simulation
{lambda_vec}
; Options for the decoupling
sc-alpha                 = 0.5
sc-coul                  = no       ; linear interpolation of Coulomb (none in this case)
sc-power                 = 1.0
sc-sigma                 = 0.3
couple-moltype           = MOL      ; name of moleculetype to decouple
couple-lambda0           = {interactions_lambda0}    ; which interactions are on
couple-lambda1           = {interactions_lambda1}    ; which interactions are on
couple-intramol          = {couple_intramol}
nstdhdl                  = {nstdhdl} ; how frequently we save free energy data
"""

EM_L_BFGS_RUN_CONTROL = """ ; Run control
integrator               = l-bfgs
nsteps                   = 5000
define                   = -DFLEXIBLE
"""

EM_STEEP_RUN_CONTROL = """; Run control
integrator               = steep
nsteps                   = 5000
"""

MD_RUN_CONTROL = """; MD Run control
integrator               = {integrator}       ; Langevin dynamics
tinit                    = 0
dt                       = {dt}             ; ps timestep
nsteps                   = {nsteps}       ; 100 ps for equilibration, 5 ns for MD
nstcomm                  = 100
"""

EM_CRITERIA = """ Options for energy minimization and EM output
emtol                    = 1000
emstep                   = 0.01
niter                    = 20
nbfgscorr                = 10
; Output control
nstlog                   = 1
nstenergy                = 1
"""

MD_OUTPUT_CONTROL = """ ; Output control
nstxout                  = 0
nstvout                  = 0
nstfout                  = 0
nstlog                   = 500
nstenergy                = 500
nstxout-compressed       = {nstxout}
energygrps          = MOL SOL
"""

NEIGHBOURSEARCH = """ Neighboursearch and short-range nondbonded interactions
cutoff-scheme            = {cutoff}
nstlist                  = {nstlist}
ns_type                  = grid
pbc                      = {pbc}
rlist                    = {rlist}
"""

ELECTROSTATICS_AND_WDW_AND_PME = """; Electrostatics
 coulombtype              = {col_type}
rcoulomb                 = {coulomb_cutoff}
epsilon-rf               = {epsilon_rf}
; van der Waals
vdwtype                  = cutoff
vdw-modifier             = {vdw_modifier}
rvdw-switch              = 1.0
rvdw                     = {vdw_cutoff}
; Apply long range dispersion corrections for Energy and Pressure
DispCorr                  = {dispcor}
; Spacing for the PME/PPPM FFT grid
fourierspacing           = {fourierspacing}  ; should be multiple of box size
; EWALD/PME/PPPM parameters
pme_order                = {pme_order}
ewald_rtol               = 1e-06
epsilon_surface          = 0
"""

TEMPERTURE_COUPLING_ON = """; Temperature coupling
; tcoupl is implicitly handled by the sd integrator
tcoupl                   = {thermostat}
tc_grps                  = system
tau-t                    = 1.0
ref-t                    = {temperature}
"""

PRESURE_COUPLING_ON = """; Presure coupling
Pcoupl                   = {barostat}   ; Parrinello-Rahman for all MD but NVT
tau_p                    = 1.0
compressibility          = 4.5e-05
ref_p                    = {pressure}
"""

TEMPERTURE_COUPLING_OFF = """; Temperature coupling is off during EM
tcoupl                   = no
"""

PRESURE_COUPLING_OFF = """; Presure coupling is off
pcoupl                   = no
"""

VELOCITY_GENERATION = """; Velocity generation
gen-vel                 = {gen_vel}        ; yes for nvt, otherwise no
gen-temp                = {temperature}    ; run temperature
"""

BOND_CONSTRAINT = """ ; Bond constraint
constraints              = {bond_contraint}     ; h-bonds, for all but L-BFGS
; Type of constraint algorithm
constraint-algorithm     = lincs
; Do not constrain the starting configuration
continuation             = {continuing_run}    ; yes or no depending on whether we continue simulation
; Highest order in the expansion of the constraint coupling matrix
lincs-order              = {lincs_order}
lincs-iter               = {lincs_iter}
"""

FREEZE = """ ; Freeze MOL
energygrp-excl      = MOL MOL
freezegrps          = MOL
freezedim           = Y Y Y
"""

EM_MDP = """
{run_control}
{em_control}
{neighboursearch}
{estat_vdw_pme}
{tcoupling}
{pcoupling}
{free_energy}
{vel_gen}
{bond_constraints}
"""

EQ_MDP = """
{run_control}
{output_control}
{neighboursearch}
{estat_vdw_pme}
{tcoupling}
{pcoupling}
{vel_gen}
{bond_constraints}
"""

MD_MDP = """
{run_control}
{output_control}
{neighboursearch}
{estat_vdw_pme}
{tcoupling}
{pcoupling}
{free_energy}
{vel_gen}
{bond_constraints}
"""



##############################################################################
                          #   Solvent stuff #
##############################################################################



# amber-03.ff spce.itp
SPCE_WATER = """[ atomtypes ]
OW_spc       8      15.9994  0.0000  A   3.16557e-01  6.50629e-01
HW_spc       1       1.0080  0.0000  A   0.00000e+00  0.00000e+00

[ moleculetype ]
; molname    nrexcl
SOL        2

[ atoms ]
; id  at type     res nr  res name  at name  cg nr  charge    mass
  1   OW_spc      1       SOL       OW       1      -0.8476   15.99940
  2   HW_spc      1       SOL       HW1      1       0.4238    1.00800
  3   HW_spc      1       SOL       HW2      1       0.4238    1.00800

#ifndef FLEXIBLE

[ settles ]
; OW	funct	doh	dhh
1       1       0.1     0.16330

[ exclusions ]
1	2	3
2	1	3
3	1	2

#else

[ bonds ]
; i     j       funct   length  force.c.
1       2       1       0.1     345000  0.1     345000
1       3       1       0.1     345000  0.1     345000

[ angles ]
; i     j       k       funct   angle   force.c.
2       1       3       1       109.47  383     109.47  383

#endif"""


CSPCE_WATER = """[ atomtypes ]
OW_spc       8      15.9994  0.0000  A   3.16557e-01  6.50629e-01
HW_spc       1       1.0080  0.0000  A   1.16580e-01  6.50629e-02

[ moleculetype ]
; molname    nrexcl
SOL        2

[ atoms ]
; id  at type     res nr  res name  at name  cg nr  charge    mass
  1   OW_spc      1       SOL       OW       1      -0.8476   15.99940
  2   HW_spc      1       SOL       HW1      1       0.4238    1.00800
  3   HW_spc      1       SOL       HW2      1       0.4238    1.00800

#ifndef FLEXIBLE

[ settles ]
; OW  funct doh dhh
1       1       0.1     0.16330

[ exclusions ]
1 2 3
2 1 3
3 1 2

#else

[ bonds ]
; i     j       funct   length  force.c.
1       2       1       0.1     345000  0.1     345000
1       3       1       0.1     345000  0.1     345000

[ angles ]
; i     j       k       funct   angle   force.c.
2       1       3       1       109.47  383     109.47  383

#endif"""


# LJ solute pdb file for packmol. Assume that the atomname is LJs
SOLV_FILENAME = 'solv.pdb'
LJ_PDB = """ATOM      1  LJs SOL     1       0.000   0.000   0.000  1.00  0.00             
TER
END"""
BOX_PDB_FILENAME = 'box.pdb'

PACMOL_INP = """# All atoms from diferent molecules will be at least 2.0 Angstroms apart
# from each other at the solution.
tolerance 2.0

filetype pdb

output {box_pdb}

# Solute is fixed
# The first three number indicate the position where the center of the
# solute will be fixed and the last three number indicate its rotation.
structure {solu_pdb}
  center
  fixed  {solu_pos[0]:.3f} {solu_pos[1]:.3f} {solu_pos[2]:.3f} 0. 0. 0.
end structure

# Solvent
structure {solv_pdb}
  number {nsolv}
  inside box 0. 0. 0. {box_vec[0]:.3f} {box_vec[1]:.3f} {box_vec[2]:.3f}
end structure
"""


##############################################################################
                            #   Globals #
##############################################################################



MDOUT_FOLDER = './mdout'
#FOLDER_NAMES = ['EM1', 'NVT', 'MD']
#FOLDER_NAMES = ['EM1', 'EM2', 'MD']
FOLDER_NAMES = ['EM1', 'NPT', 'MD']

MOL_GRO = 'MOL.gro'
MOL_TOP = 'MOL.top'
#WATER_FOR_BOX = 'spc216.gro'
MDP_FILES_ROOT = './mdp'
MOL_FILES_DIR = './MOL'
XVG_FILES = './xvg'


##############################################################################
                              #   Code #
##############################################################################



def process_command_line(argv):
    """Processes arguments and returns namespace of them."""
    parser = argparse.ArgumentParser(description=""" Run a molecular dynamics
    simulation with gromacs.""")
    #Positional args
    parser.add_argument('crd', help='Amber crd or pdb file', metavar='mol.inpcrd')
    parser.add_argument('top', help='Amber prmtop', metavar='mol.prmtop')
    #Optional args
    # computational resources
    resources = parser.add_argument_group('Computational resources',
                            """Set number of processors/nodes.""")
    resources.add_argument('-p', '--procs', help="""Number of nodes to
                        use during md run. On Archie each node has 12 cores, thus
                        n_procs = PROCS*12.
                        There are no nodes on ferrari, but this number will
                        still be multiplied by 12. If multidir
                        option is selected, run is going to set ncores = n lambda,
                        but this number is going to be passed to que system [2]""",
                        default='2')
    resources.add_argument('--multidir',help=""" Run each lambda on a single cpu instead of 
                        running 1 lambda on multiple CPUs.""", action='store_true')
    resources.add_argument('-cs', '--custom_cores', help=""" When you need to specify
                        exact number of processors for all mdruns. Useful on 
                        ferrari.""",
                        type=int)
    # preprocessing otpions
    preprocess = parser.add_argument_group('Preprocessing options',
                        """Control grompp and solvation parameters. Solvation
                        can be done either via gmx solvate (called by 
                        setting a single box_size in --box_size option) or 
                        packmol (called by setting 3 box sizes in --box_size 
                        option).""")
    preprocess.add_argument('--maxwarn', help=""" Ignore n warnings [0]. Applies
                        to every grompp call, including energy minimization.""",
                        default=0, type=int)
    preprocess.add_argument('-st', '--solvent_top', help=""" Topology of the
                        solvent. Script has predefined topologies for spc/e
                        and cspc/e water models which can be used by 
                        specifying spce or cspce repsectively. Othewise
                        path to existing .itp file should be provided [spce].""",
                        default='spce')
    preprocess.add_argument('-sp', '--solvent_pdb', help=""" pdb structure of
                        solvent to use with packmol. If none is provided
                        will assume we are dealing with spce water.""")
    preprocess.add_argument('-sb', '--solvent_box', help=""" Path to the box
                        with pre-equilibrated solvent. This option is not compatible
                        with pacmol and can be only used with gmx solvate.
                        [spc216.gro].""",
                        default='spc216.gro')
    preprocess.add_argument('--box_size', help="""Box vectors. Units are nm.
                        Specifying all three dimensions results in packmol call.
                        Specifying a single size calls gmx solvate.
                        [3.22]""", default=['3.22'],
                        metavar='nm', nargs='+')
    preprocess.add_argument('--max_solv', help=""" Maximum number of solvent molecules
                        in simulation box. 0 for None. [1024] """, 
                        default='1024') 
    preprocess.add_argument('--box_type', help=""" Type of box to simulate in.
                        Chose from {cubic, dodecahedron} [cubic]""",
                        default='cubic')
    preprocess.add_argument('--pbc', help=""" Specify pbc. Can be either {no,
                        xyz, xy} [xyz]. """, 
                        default='xyz')
    preprocess.add_argument('--no_solvate', help=""" Do not solvate. """, 
                        action='store_true')
    preprocess.add_argument('--nstxout', help=""" Trajectory output frequency [0]. """, 
                        default='0')
#    preprocess.add_argument('--nb_list', help="""Provide a file containing a list
#                        of non-bonded interactions that will be appended to .top
#                        file. """)
    preprocess.add_argument('--acpype', help="""Path to acpype.py file.""",
                        metavar='acpype_path')
    preprocess.add_argument('--submit', help="""Start calculation from the script. """, 
                        action='store_true')
    preprocess.add_argument('--dir_name', help="""Name of the calculation directory.""")
    preprocess.add_argument('--verbose', help="""Be verbose.""", 
                        action='store_true')
    # md dynamics parameters
    dynamics = parser.add_argument_group('MD dynamic options',
                            """Controls dynamics and baro/thermostats.""")
    dynamics.add_argument('--equilibration', help=""" Run pressure 
                        equilibration with berendsen barostat at 1 bar.""",
                        action='store_true')
    dynamics.add_argument('--integrator', help=""" Integrator for production run 
                         [sd]""",
                        default='sd')
    dynamics.add_argument('--pcoupl', help=""" Barostat for production run.
                        Equilibration is always done using berendsen barostat.
                        For NVT run simply use no.
                        {no, berendsen, Parrinello-Rahman}.
                        [berendsen] """, default='berendsen')
    dynamics.add_argument('--pressure', help=""" System pressure in bars [1.0].
                        """, default=1.0, type=float)
    dynamics.add_argument('--tcoupl', help=""" Default integrator is sd, thus
                        thermostat is implicitly handled by it. For non-sd 
                        integrators you can specify:
                        {no, v-rescale, nose-hoover}.
                        [no] """, default='no')
    dynamics.add_argument('-t', '--temperature', help='Simulation temperature K [298.15]',
                        default='298.15')
    dynamics.add_argument('--nsteps', help=""" Number of calculation steps [750000] """,
                        type=int, default=750000)
    dynamics.add_argument('--eqnsteps', help=""" Number of equlibration steps [100000] """,
                        type=int, default=100000)
    dynamics.add_argument('--time_step', help=""" Time step in fs for both
                        equilibration and production [2].""",
                        type=float, default=2.0)
    # bond constraints
    constr = parser.add_argument_group('Bond constraints',
                            """Types and accuracy of constraints.""")
    constr.add_argument('--bond_const', help=""" Bond constraints
                        Chose from {none, h-bonds, all-bonds, h-angles, all-angles}
                        [h-bonds]""",
                        default='h-bonds')
    constr.add_argument('--freeze', help=""" Freeze and exlude MOL from energy 
                        grps. Only done for production run. Usually works
                        only with barostat turned off.""",
                        action='store_true')
    constr.add_argument('--lincs_order', help=""" Decrease for better paralelisation [4].""",
                        default='4')
    constr.add_argument('--lincs_iter', help=""" Increase for better accuracy [1].""",
                        default='1')
    # non-bonded interactions
    nb_inter = parser.add_argument_group('Non-bonded interactions',
                            """Cutoffs, coulomb, vdws, pme.""")
    nb_inter.add_argument('--cutoff', help=""" Cut-off scheme {Verlet, group}.
                        [verlet] """, default='verlet')
    nb_inter.add_argument('--rlist', help=""" Rlist size for both cutoff and
                        coulomb cutoff. Set larger
                        value for group scheme [1.2].""", default='1.2')
    nb_inter.add_argument('--coulomb_cutoff', help=""" Real space cutoff
                         coulomb interactions [1.2].""", default='1.2')
    nb_inter.add_argument('--vdw_cutoff', help=""" Real space cutoff
                         vdw interactions [1.2].""", default='1.2')
    nb_inter.add_argument('--vdw_modifier', help=""" Controls how vdw potential
                        behaves. Sensible options are (potential-switch: 
                        smoothly decay vdw potential between 1.0 nm and
                        vdw_cutoff; Potential-shift-Verlet : shifts potential so
                        it is 0 at cutoff on the whole range 
                        [Potential-shift-Verlet]""", default='Potential-shift-Verlet')
    nb_inter.add_argument('--nstlist', help=""" Nstlist size [40].""",
                        default='40')
    nb_inter.add_argument('--col_type', help=""" Coulmbtype {Cut-off, PME, ..}
                        For full list of options check mdp manual [PME]. """, 
                        default='pme')
    nb_inter.add_argument('--epsilon_rf', help=""" Dielectric constant outside
                        coulomb reaction field. Default value corresponds to
                        SPC-E water [73.5]. """, 
                        default=73.5, type=float)
    nb_inter.add_argument('--fourierspacing', help=""" Spacing between grid points
                        in PME. Increase for parallel [0.15]. """, 
                        default=0.15, type=float)
    nb_inter.add_argument('--pme_order', help=""" Order of interpolation for 
                        particle mesh ewald [6]. """, 
                        default=6, type=int)
    nb_inter.add_argument('--dispcor', help=""" Dispersion correction {no, EnerPres}.
                        [EnerPres] """, default='EnerPres')
    # free energy options
    fe_options = parser.add_argument_group('Free energy options',
                            """Number of lambdas and xvg related options.""")
    fe_options.add_argument('--fe_type', help=""" Which interactions to decouple.
                        Choose from {all, all_s, vdw, vdw_s, q, q_s, q_long, 
                        q_int, vdw_int}.
                        Subscript s indicates that calculation is going to
                        use fewer lambdas. In q_int int stands for integer
                        and means that simulation is going to be run for int
                        number of lambdas equally spread between 0 and 1.
                        [all_s]""", default='all_s')
    fe_options.add_argument('--calc_lambda_neighbors', help=""" Number of
                        neighbour lambdas to calculate [-1].""",
                        default='-1')
    fe_options.add_argument('--subset', help="""Run only for certain lambdas.""",
                        metavar='lambda',  nargs='+')
    fe_options.add_argument('--couple_intramol', help=""" Decouple intermolecular
                        interactions in addition to intramolecular. Yes or no
                         [no]. """, default='no'
                        )
    return parser.parse_args(argv)


def calc_fe_paramters(fe_type):
    """ Check the fe type of calc and return lambda parameters """
    if fe_type == 'all':
        interactions_lambda0 = 'vdw-q'
        interactions_lambda1 = 'none'
        lambda_vec = LAMBDA_ALL[0]
        lambda_states = LAMBDA_ALL[1] + 1  # +1 cuz we start from 0
    elif fe_type == 'all_s':
        interactions_lambda0 = 'vdw-q'
        interactions_lambda1 = 'none'
        lambda_vec = LAMBDA_ALL_S[0]
        lambda_states = LAMBDA_ALL_S[1] + 1  # +1 cuz we start from 0
    elif fe_type == 'q':
        interactions_lambda0 = 'vdw-q'
        interactions_lambda1 = 'vdw'
        lambda_vec = LAMBDA_Q[0]
        lambda_states = LAMBDA_Q[1] + 1
    elif fe_type == 'q_s':
        interactions_lambda0 = 'vdw-q'
        interactions_lambda1 = 'vdw'
        lambda_vec = LAMBDA_Q_S[0]
        lambda_states = LAMBDA_Q_S[1] + 1  # +1 cuz we start from 0
    elif fe_type == 'q_long':
        interactions_lambda0 = 'vdw-q'
        interactions_lambda1 = 'vdw'
        lambda_vec = LAMBDA_Q_LONG[0]
        lambda_states = LAMBDA_Q_LONG[1] + 1  # +1 cuz we start from 0
    elif fe_type.startswith('q_') and fe_type[2:].isdigit():
        interactions_lambda0 = 'vdw-q'
        interactions_lambda1 = 'vdw'
        lambda_states = int(fe_type[2:])
        if lambda_states < 2:
            raise ValueError('Number of states should be bigger than 1.')
        coul_lambdas = [round(float(i)/(lambda_states-1), 3) for i in range(lambda_states)]
        lambda_vec = 'coul_lambdas      = ' + ' '.join(map(str, coul_lambdas))
    elif fe_type == 'vdw':
        interactions_lambda0 = 'vdw'
        interactions_lambda1 = 'none'
        lambda_vec = LAMBDA_VDW_S[0]
        lambda_states = LAMBDA_VDW_S[1] + 1
    elif fe_type == 'vdw_s':
        interactions_lambda0 = 'vdw'
        interactions_lambda1 = 'none'
        lambda_vec = LAMBDA_VDW_S[0]
        lambda_states = LAMBDA_VDW_S[1] + 1
    elif fe_type.startswith('vdw_') and fe_type[4:].isdigit():
        interactions_lambda0 = 'vdw'
        interactions_lambda1 = 'none'
        lambda_states = int(fe_type[4:])
        if lambda_states < 2:
            raise ValueError('Number of states should be bigger than 1.')
        coul_lambdas = [round(float(i)/(lambda_states-1), 3) for i in range(lambda_states)]
        lambda_vec = 'vdw_lambdas       = ' + ' '.join(map(str, coul_lambdas))
    else:
        raise ValueError('Unknown free energy type (fe_type).')
    return interactions_lambda0, interactions_lambda1, lambda_vec, lambda_states


def mdrun(job_list, job_name, args):
    """ Run mdrun. 
    job_list is a list containing tuples with following data:
    lambda_state, ((grompp1 call, mdrun1), (grompp2 call, mdrun2), ...)
    """
    # check where we are
    if args.host  == 'ferrari':
        job_script_header = JOB_SCRIPT_HEADER_FERRARI
        run_cmd = 'bash'
        mdrun_exe = 'mdrun'
        single_run_string = mdrun_exe + \
                       ' -nt {ncores} -deffnm {mdrun_name} -s {mdrun_name}.tpr \n'

    elif args.host == 'archie':
        job_script_header = JOB_SCRIPT_HEADER_ARCHIE
        run_cmd = 'qsub'
        mdrun_exe = 'mdrun_mpi'
        single_run_string = 'mpirun -np {ncores} ' + mdrun_exe + \
                       ' -deffnm {mdrun_name} -s {mdrun_name}.tpr \n'
    p, _ = os.path.split(job_name)
    if not p:
        p = '.'
    if args.verbose:
        single_run_string += '-v'
    run_script_name = '{}_mdrun_script.sh'.format(job_name)
    if run_script_name[0].isdigit():
        # on archie script name cannot start with digit
        run_script_name = 'n' + run_script_name
    run_script_name = os.path.join(p, run_script_name)
    # determine number of procs
    if args.custom_cores:
        ncores = args.custom_cores
    else:
        if args.multidir:
            ncores = len(job_list)
        else:
            ncores = '$NCORES'
    run_strings = ''
    if args.multidir:
        if args.host == 'ferrari':
            raise ValueError('Caution, Not tested on ferrari')
        # grompp calls + names separated by type
        mdrun_series = [[] for i in job_list[0][1]]  
        for _, mdrun_pairs in job_list:
            for i, mdrun_pair in enumerate(mdrun_pairs):
                mdrun_series[i].append(mdrun_pair)
        for job_collection in mdrun_series:
            name = os.path.split(job_collection[0][1])[1]
            base = 'mpirun -np {} {} -deffnm {} -multidir '.format(ncores, mdrun_exe, name)
            folders = ''
            for grompp_call, mdrun_name in job_collection:
                folders += os.path.split(mdrun_name)[0] + ' '
                run_strings += ' '.join(grompp_call) + '\n'
            run_strings += base + folders + '\n\n'
    else:
        for lambda_state, mdrun_names  in job_list:
            for grompp_call, mdrun_name in mdrun_names:
                run_strings += ' '.join(grompp_call) + '\n'
                run_strings += single_run_string.format(ncores=ncores, mdrun_name=mdrun_name)
    run_script_header = job_script_header.format(nnodes=args.procs)
    run_script = run_script_header + '\n' + run_strings
    with open(run_script_name, 'w') as f:
        f.write(run_script)
    if args.submit:
        subprocess.Popen([run_cmd, run_script_name])
    

def convert_to_grom(args):
    """Takes as input amber crd and top files, returns gromacs
    gro and top files."""
    # find acpype first
    if args.acpype:
        acpype = args.acpype
    else:
        if args.host == 'ferrari':
            acpype = '/home/max/opt/acpype/acpype.py'
        elif args.host == 'archie':
            acpype = '/users/xpb13212/opt/acpype/acpype.py'
    # flag b ensures that the name is MOL        
    subprocess.call(['python', acpype, '-p', args.top, '-x', args.crd, '-b', 'MOL'])
    return './MOL_GMX.gro', './MOL_GMX.top'


def create_box_and_solvate(gmx_gro, gmx_top, args):
    """Solvate molecule with 3 point water model. Everything takes place
    insed the job directory."""
    if args.host == 'ferrari':
        singlpnt_gmx = GMX_FERRARI
    elif args.host == 'archie':
        singlpnt_gmx = GMX_ARCHIE
    if len(args.box_size) == 1:
        # Use standard gromacs solvation in which we assume that exact
        # density equilibration will be performed for every lambda
        subprocess.call(singlpnt_gmx + [ 'editconf', '-f', gmx_gro, '-c',
                        '-bt', args.box_type, '-o', 'MOL_BOX.gro', '-box'] + args.box_size)
        if args.no_solvate:
            shutil.move('MOL_BOX.gro', MOL_GRO)
        else:
            subprocess.call(singlpnt_gmx + ['solvate', '-cp', 'MOL_BOX.gro', 
                            '-cs', args.solvent_box,
                            '-o', MOL_GRO, '-p', gmx_top, '-maxsol', args.max_solv,
                            '-box'] + args.box_size )        
        shutil.copy(gmx_top, MOL_TOP)
    elif len(args.box_size) == 3:
        # Specific number of molecules and box geometry required
        # pacmol needed
#        if args.solvent_box:
#            raise ValueError('Solvent box and packmol are not compatible.')
        if not args.pdb:
            raise ValueError('Pacmol only works with pdb files')
        # create a solvent molecule
        if args.solvent_pdb:
            # we don't need to know real name of working directory,
            # so 'job_dir' works
            solv_relpath = os.path.relpath(args.solvent_pdb, 'job_dir')
            shutil.copy(solv_relpath, SOLV_FILENAME)
        else:
            with open(SOLV_FILENAME, 'w') as f:
                f.write(LJ_PDB)
        # convert box vec to angstroms for pacmol
        box_vec = map(lambda x: float(x)*10, args.box_size)
        solu_pos = [i/2. for i in box_vec]
        packmol_script = PACMOL_INP.format(solu_pdb=args.pdb, solv_pdb=SOLV_FILENAME,
                                           box_pdb=BOX_PDB_FILENAME,
                                           box_vec=box_vec, solu_pos=solu_pos,
                                           nsolv=args.max_solv)
        p = subprocess.Popen(['packmol'], stdin=subprocess.PIPE, 
                             stdout=subprocess.PIPE)
        packmol_out = p.communicate(packmol_script)[0]
        with open('packmol.log', 'w') as f:
            f.write(packmol_out)
        p = subprocess.Popen(['babel', '-ipdb', BOX_PDB_FILENAME, '-ogro', MOL_GRO],
                             stderr=subprocess.PIPE)
        p.communicate()
        # fix box dimensions
        with open(MOL_GRO, 'r') as f:
            txt = f.readlines()
        txt[-1] = '   {0[0]:.5f}   {0[0]:.5f}   {0[0]:.5f}\n'.format(
                                                map(float, args.box_size))
        with open(MOL_GRO, 'w') as f:
            f.writelines(txt)
        # fix toplogy
        shutil.copy(gmx_top, MOL_TOP)
        with open(MOL_TOP, 'a') as f:
            f.write('SOL     ' + args.max_solv + '\n')
    return MOL_GRO, MOL_TOP


def clean_gro_and_top_files(gro_file, top_file, amber_crd, amber_top):
    """Creates MOL_FILE_DIR and then copies gro file, top file and
    amber input files there.
    Afterwards removes intermediate files."""
    try_to_make_dir(MOL_FILES_DIR)
    try:
        shutil.copy(gro_file, MOL_FILES_DIR)
        shutil.copy(top_file, MOL_FILES_DIR)
        shutil.move(amber_top, MOL_FILES_DIR)
        shutil.move(amber_crd, MOL_FILES_DIR)
    except shutil.Error:
        pass
    del_cmd = ['rm', 'em.mdp', 'md.mdp']
    del_cmd.extend(glob.glob('*.gro'))
    del_cmd.extend(glob.glob('*.top*'))
    subprocess.call(del_cmd)
    return os.path.join(MOL_FILES_DIR, gro_file), \
           os.path.join(MOL_FILES_DIR, top_file)


def gen_em_mdp_files(lambda_state, args):
    """Create steep andL-BFGS energy minimization mdp files."""
    interactions_lambda0, interactions_lambda1, lambda_vec, _ = \
                                                calc_fe_paramters(args.fe_type)
    # if args.freeze:
    #     BOND_CONST = BOND_CONSTRAINT + FREEZE
    # else:
    #     BOND_CONST = BOND_CONSTRAINT
    em_mdp_files = []
    em_steep_file = EM_MDP.format(run_control=EM_STEEP_RUN_CONTROL,
                                   em_control=EM_CRITERIA,
                                   neighboursearch=NEIGHBOURSEARCH,
                                   estat_vdw_pme=ELECTROSTATICS_AND_WDW_AND_PME,
                                   tcoupling=TEMPERTURE_COUPLING_OFF,
                                   pcoupling=PRESURE_COUPLING_OFF,
                                   free_energy=FREE_ENERGY,
                                   vel_gen=VELOCITY_GENERATION,
                                   bond_constraints=BOND_CONSTRAINT)
#    em_l_bfgs_file = EM_MDP.format(run_control=EM_L_BFGS_RUN_CONTROL,
#                                   em_control=EM_CRITERIA,
#                                   neighboursearch=NEIGHBOURSEARCH,
#                                   estat_vdw_pme=ELECTROSTATICS_AND_WDW_AND_PME,
#                                   tcoupling=TEMPERTURE_COUPLING_OFF,
#                                   pcoupling=PRESURE_COUPLING_OFF,
#                                   free_energy=FREE_ENERGY,
#                                   vel_gen=VELOCITY_GENERATION,
#                                   bond_constraints=BOND_CONSTRAINT)
    em_steep_file = em_steep_file.format(nstlist=NSTLIST_EM,
                          coulomb_cutoff=args.coulomb_cutoff,
                          vdw_cutoff=args.vdw_cutoff,
                          vdw_modifier=args.vdw_modifier,
                          temperature=args.temperature,
                          gen_vel='no',
                          bond_contraint=args.bond_const,
                          lincs_order=args.lincs_order,
                          lincs_iter=args.lincs_iter,
                          continuing_run='no',
                          nstdhdl=NSTDHDL,
                          pbc=args.pbc,
                          dispcor=args.dispcor,
                          epsilon_rf=args.epsilon_rf,
                          fourierspacing=args.fourierspacing,
                          pme_order=args.pme_order,
                          cutoff=args.cutoff,
                          rlist=args.rlist,
                          col_type=args.col_type,
                          couple_intramol=args.couple_intramol,
                          lambda_state=lambda_state,
                          interactions_lambda0 = interactions_lambda0,
                          interactions_lambda1 = interactions_lambda1,
                          calc_lambda_neighbors = args.calc_lambda_neighbors,
                          lambda_vec = lambda_vec 
                          )
    em_mdp_files.append(em_steep_file)
#    em_l_bfgs_file = em_l_bfgs_file.format(nstlist=NSTLIST_EM,
#                          columb_cutoff=COLUMB_CUTOFF,
#                          vdw_cutoff=VDW_CUTOFF,
#                          temperature=args.temperature,
#                          gen_vel='no',
#                          bond_contraint='none',
#                          continuing_run='no',
#                          nstdhdl=NSTDHDL,
#                          dispcor=args.dispcor,
#                          cutoff=args.cutoff,                          
#                          pbc=args.pbc,
#                          col_type=args.col_type,                      
#                          lambda_state=lambda_state,
#                          interactions_lambda0 = interactions_lambda0,
#                          interactions_lambda1 = interactions_lambda1,
#                          lambda_vec = lambda_vec 
#                          )
#    return em_steep_file, em_l_bfgs_file
    return em_mdp_files


def gen_md_mdp_files(lambda_state, nsteps, dt, args):
    """Create mdp files."""
    interactions_lambda0, interactions_lambda1, lambda_vec, _ = \
                                                calc_fe_paramters(args.fe_type)
    if args.freeze:
        BOND_CONST = BOND_CONSTRAINT + FREEZE
    else:
        BOND_CONST = BOND_CONSTRAINT
    md_mdp_files = []
    eq_file = EQ_MDP.format(run_control=MD_RUN_CONTROL,
                                    output_control=MD_OUTPUT_CONTROL,
                                    neighboursearch=NEIGHBOURSEARCH,
                                    estat_vdw_pme=ELECTROSTATICS_AND_WDW_AND_PME,
                                    tcoupling=TEMPERTURE_COUPLING_ON,
                                    pcoupling=PRESURE_COUPLING_ON,
                                    vel_gen=VELOCITY_GENERATION,
                                    bond_constraints=BOND_CONSTRAINT)   # no freeze here
    md_file = MD_MDP.format(run_control=MD_RUN_CONTROL,
                                    output_control=MD_OUTPUT_CONTROL,
                                    neighboursearch=NEIGHBOURSEARCH,
                                    estat_vdw_pme=ELECTROSTATICS_AND_WDW_AND_PME,
                                    tcoupling=TEMPERTURE_COUPLING_ON,
                                    pcoupling=PRESURE_COUPLING_ON,
                                    free_energy=FREE_ENERGY,
                                    vel_gen=VELOCITY_GENERATION,
                                    bond_constraints=BOND_CONST)
#    nvt_file = generic_md_file.format(nsteps=EQUILIBRATION_N_STEPS,
#                                      nstlist=NSTLIST,
#                                      columb_cutoff=COLUMB_CUTOFF,
#                                      vdw_cutoff=VDW_CUTOFF,
#                                      temperature=temperature,
#                                      barostat='no',
#                                      gen_vel='yes',
#                                      bond_contraint='h-bonds',
#                                      continuing_run='no',
#                                      nstdhdl=NSTDHDL,
#                                      lambda_state=lambda_state)
    if args.equilibration:
        npt_file = eq_file.format(nsteps=args.eqnsteps,
                                  dt=dt,
                                  nstlist=args.nstlist,
                                  coulomb_cutoff=args.coulomb_cutoff,
                                  vdw_cutoff=args.vdw_cutoff,
                                  vdw_modifier=args.vdw_modifier,
                                  pbc=args.pbc,
                                  integrator=args.integrator,
                                  # don't write equilibration trajectory
                                  nstxout='0',
                                  col_type=args.col_type,
                                  dispcor=args.dispcor,
                                  cutoff=args.cutoff,
                                  rlist=args.rlist,
                                  epsilon_rf=args.epsilon_rf,
                                  # new run should always start with berendsen
                                  # and this is npt, so barostat is always there
                                  pressure=args.pressure,
                                  barostat='berendsen',      
                                  temperature=args.temperature,
                                  thermostat=args.tcoupl,
                                  gen_vel='yes',
                                  bond_contraint=args.bond_const,
                                  lincs_order=args.lincs_order,
                                  lincs_iter=args.lincs_iter,
                                  continuing_run='no',
                                  fourierspacing=args.fourierspacing,
                                  pme_order=args.pme_order,
                                  )
        md_mdp_files.append(npt_file)
    md_file = md_file.format(nsteps=nsteps,
                              dt=dt,
                              nstlist=args.nstlist,
                              coulomb_cutoff=args.coulomb_cutoff,
                              vdw_cutoff=args.vdw_cutoff,
                              vdw_modifier=args.vdw_modifier,
                              temperature=args.temperature,
                              thermostat=args.tcoupl,
                              pressure=args.pressure,
                              barostat=args.pcoupl,
                              pbc=args.pbc,
                              integrator=args.integrator,
                              nstxout=args.nstxout,
                              col_type=args.col_type,
                              dispcor=args.dispcor,
                              cutoff=args.cutoff,
                              rlist=args.rlist,
                              epsilon_rf=args.epsilon_rf,
                              gen_vel='no',
                              continuing_run='yes',
                              bond_contraint=args.bond_const,
                              lincs_order=args.lincs_order,
                              lincs_iter=args.lincs_iter,
                              fourierspacing=args.fourierspacing,
                              pme_order=args.pme_order,
                              nstdhdl=NSTDHDL,
                              calc_lambda_neighbors = args.calc_lambda_neighbors,
                              lambda_state=lambda_state,
                              interactions_lambda0 = interactions_lambda0,
                              interactions_lambda1 = interactions_lambda1,
                              couple_intramol=args.couple_intramol,
                              lambda_vec = lambda_vec 
                              )
    md_mdp_files.append(md_file)
    #return nvt_file, npt_file, md_file
    return md_mdp_files


def modify_topology(top_file, args):
    """Add a solvent itp (spce by default) to topology file."""
    with open(top_file, 'rb') as f:
        top_text = f.read()
    top_sections = top_text.split('\n\n')
    if args.solvent_top == 'spce':
        solvent_def = SPCE_WATER
    elif args.solvent_top == 'cspce':
        solvent_def = CSPCE_WATER
    else:
        with open(os.path.relpath(args.solvent_top, args.name)) as f:
            solvent_def = f.read() 
    top_sections.insert(3, solvent_def)
    new_top = '\n\n'.join(top_sections)
    with open(top_file, 'wb') as f:
        f.write(new_top)
    return top_file


def generate_mdp_files(lambda_states_to_evavluate_on, nsteps, dt, args):
    """Create a tree of mdp files:
    MDP_FILES_ROOT/
    |---EM1/
        |---em_steep_{lambda_state}.mdp
    |---EM2/
        |---em_l_bfgs_{lambda_state}.mdp
    |---NVT
        |---nvt_{lambda_state}.mdp
    |---NPT
        |---npt_{lambda_state}.mdp
    |---MD
        |---md_{lambda_state}.mdp
    """
    # make root folder for mdps
    try_to_make_dir(MDP_FILES_ROOT)
    em1_dirname, npt_dirname, md_dirname = FOLDER_NAMES
    if args.equilibration:
        folders = FOLDER_NAMES
    else:
        folders = [em1_dirname, md_dirname]
    # create individual dirs and define paths
    mdp_paths = []
    em1_path = os.path.join(MDP_FILES_ROOT, em1_dirname, 'em_steep_{lambda_state}.mdp')
    try_to_make_dir(os.path.join(MDP_FILES_ROOT, em1_dirname))
    mdp_paths.append(em1_path)
#    em2_path = os.path.join(MDP_FILES_ROOT, em2_dirname, 'em_l_bfgs_{lambda_state}.mdp')
#    try_to_make_dir(os.path.join(MDP_FILES_ROOT, em2_dirname))
#    nvt_path = os.path.join(MDP_FILES_ROOT, nvt_dirname, 'nvt_{lambda_state}.mdp')
#    try_to_make_dir(os.path.join(MDP_FILES_ROOT, nvt_dirname))
    if args.equilibration:
        npt_path = os.path.join(MDP_FILES_ROOT, npt_dirname, 'npt_{lambda_state}.mdp')
        try_to_make_dir(os.path.join(MDP_FILES_ROOT, npt_dirname))
        mdp_paths.append(npt_path)
    md_path = os.path.join(MDP_FILES_ROOT, md_dirname, 'md_{lambda_state}.mdp')
    try_to_make_dir(os.path.join(MDP_FILES_ROOT, md_dirname))
    mdp_paths.append(md_path)
    # generate necessary mdp files and write them where needed
    for lambda_state in lambda_states_to_evavluate_on:
        em_mdp_files = gen_em_mdp_files(lambda_state, args)
        md_mdp_files = gen_md_mdp_files(lambda_state, nsteps, dt, args)
        mdp_files = em_mdp_files + md_mdp_files
        for mdp_path, mdp_file in zip(mdp_paths, mdp_files):
            with open(mdp_path.format(lambda_state=lambda_state), 'wb') as f:
                f.write(mdp_file)
    return folders, mdp_paths


def run_md_simulations(job_name, lambda_states_to_evavluate_on,
                       gro_file, top_file, mdp_files, folders,
                       args):
    """Preforms all MD simulation (takes quite a while)"""
    if args.host == 'ferrari':
        singlpnt_gmx = GMX_FERRARI
        #thread_specifier = '-nt'
    elif args.host == 'archie':
        singlpnt_gmx = GMX_ARCHIE
        #thread_specifier = '-ntomp'
#    log_file = open('md_{}.log'.format(job_name), 'wb')
    jobs_list = []
    for lambda_state in lambda_states_to_evavluate_on:
        print 'Starting Lambda state: {}'.format(lambda_state)
#        log_file.write("""\n##########################\n\n\n\n
#                       Lambda state {}:""".format(lambda_state))
        lambda_dir = 'Lambda_{}'.format(lambda_state)
        try_to_make_dir(lambda_dir)
        lambda_job_list = []
        input_gro_file = gro_file  # first input gro file is simply solvated MOL        
        for dirname, mdp_name in zip(folders, mdp_files):
            job_directory = os.path.join(lambda_dir, dirname)
            mdrun_name = os.path.join(job_directory,
                                    '{}_{}'.format(dirname.lower(),lambda_state))
            try_to_make_dir(job_directory)
            mdout_mdp_name = os.path.join(MDOUT_FOLDER, job_name + '_{}'.format(lambda_state))
            grompp_call = singlpnt_gmx + ['grompp',
                           '-f', mdp_name.format(lambda_state=lambda_state),
                           '-po', mdout_mdp_name,
                           '-c', input_gro_file,
                           '-p', top_file,
                           '-o', mdrun_name + '.tpr']
            if args.maxwarn:
                grompp_call.extend(['-maxwarn', '{}'.format(args.maxwarn)])
#            if dirname == 'EM2':
#                grompp_call.extend(['-maxwarn', '1'])
#            if dirname == 'NPT':
#                nvt_trajectory = os.path.join(lambda_dir,
#                                              'NVT',
#                                              'nvt.{}'.format(lambda_state))
#                grompp_call.extend(['-t', nvt_trajectory])
            if dirname == 'MD' and args.equilibration:
                npt_trajectory = os.path.join(lambda_dir,
                                              'NPT',
                                              'npt_{}.cpt'.format(lambda_state))
                grompp_call.extend(['-t', npt_trajectory])
            if dirname.startswith('EM'):
                subprocess.call(grompp_call)
                em_call = singlpnt_gmx + ['mdrun', '-ntomp', '1',
                              '-deffnm', mdrun_name, '-s', mdrun_name + '.tpr']
                if args.verbose:
                    em_call.append('-v')
                subprocess.call(em_call)
            else:
                lambda_job_list.append((grompp_call, mdrun_name))
            # If the run will be successful new coordinates will be written to
            # the output gro file.
            input_gro_file = mdrun_name + '.gro'
        # end for
        jobs_list.append((lambda_state, lambda_job_list))
    mdrun(jobs_list, job_name, args)
#    log_file.close()


def try_to_make_dir(dirname):
    """Tries to make directory, passes in case it exists."""
    try:
        os.mkdir(dirname)
    except OSError, e:
        if e.errno == 17:
            pass # Directory already exists, all is well
        else:
            raise e


def check_input(args):
    """Checks if amber files are here."""
    if not os.path.isfile(args.top) or not os.path.isfile(args.crd):
        raise ValueError('Can not find amber files!')
    

def prepare_dir(args):
    """creatates or cd's to the directory where whole calculation is organized."""
    if args.dir_name:        
        job_name = args.dir_name
    else:
        job_name = args.top[:-7]
    _, _, _, lambda_states = calc_fe_paramters(args.fe_type)
    if args.subset:
        lambda_states_to_evavluate_on = args.subset
    else:
        lambda_states_to_evavluate_on = list(range(lambda_states))
    print 'Starting job {}'.format(job_name)
    print 'Simulation will be run for {} states'.format(lambda_states_to_evavluate_on)
    try_to_make_dir(job_name)
    incrd_f = os.path.join(job_name, args.crd)
    prmtop_f = os.path.join(job_name, args.top)
    shutil.copy(args.crd, incrd_f)
    shutil.copy(args.top, prmtop_f)
    os.chdir(job_name)
    try_to_make_dir(MDOUT_FOLDER)    
    return job_name, lambda_states_to_evavluate_on
    
    
def calculate_free_energy(mbar_path, temperature):
    try_to_make_dir(XVG_FILES)
    xvg_files = glob.glob('./*/*/md*xvg')
    for f in xvg_files:
        shutil.copy(f, XVG_FILES)
    os.chdir(XVG_FILES)
    mbar_log = subprocess.check_output(['python', mbar_path, '-p', 'md',
                                        '-t', temperature, '-u', 'kcal',
                                        '-s', '0'])
    with open('mbar_out.log', 'w') as f:
        f.write(mbar_log)


def main(argv):
    print('cli args: {}').format(argv[1:])
    args = process_command_line(argv)
    check_input(args)
    # Check where we are
    import socket
    host = socket.gethostname()
    if host  == 'ferrari':
        args.host = host
    elif host.startswith('archie'):
        args.host = 'archie'
    else:
        raise ValueError("I don't know on which computer I'm being run. Set appropriate hostname.")
    # prepare directories
    job_name, lambda_states_to_evavluate_on = prepare_dir(args)
    args.name = job_name
    # convert to inpcrd if needed
    if args.crd.endswith('.pdb'):
        data = []
        with open(args.crd) as f:
            txt = f.readlines()
        for l in txt:
            if l.startswith('ATOM'):
                data.extend(l.split()[5:8])
        new_crdname = args.crd[:-3] + 'incrd'
        with open(new_crdname, 'w') as f:
            f.write('MOL\n')
            f.write('{:6d}\n'.format(int(len(data)/3)))
            for i, coord in enumerate(data):
                #print coord
                coord = float(coord)
                f.write('{: 12.7f}'.format(coord))
                if (i+1)%6 == 0:
                    f.write('\n')
        args.pdb = args.crd
        args.crd = new_crdname
    # md paramters
    dt = args.time_step
    nsteps = args.nsteps
    dt = dt/1000   # convert to ps
    # prepare and run calc
    gmx_gro, gmx_top = convert_to_grom(args)
    gro_file, top_file = create_box_and_solvate(gmx_gro, gmx_top, args)
    if not args.no_solvate:
        top_file = modify_topology(top_file, args)
    gro_file, top_file = clean_gro_and_top_files(gro_file, top_file,
                                                 args.crd, args.top)
    folders, mdp_files = generate_mdp_files(lambda_states_to_evavluate_on,
                                   nsteps, dt, args)
    run_md_simulations(job_name, lambda_states_to_evavluate_on,
                       gro_file, top_file, mdp_files, folders,
                       args)
    #calculate_free_energy(args.mbar, args.temperature)
    

if __name__ == '__main__':
    main(sys.argv[1:])

