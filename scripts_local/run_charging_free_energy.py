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
                              # Gromacs on archie #
##############################################################################

GMX = ['gmx',]  # For 1 core stuff

JOB_SCRIPT_HEADER = r"""
# supress dump
export GMX_SUPPRESS_DUMP=Y
# start jobs
"""


##############################################################################
                                #   MDP Options  #
##############################################################################

# Free energy stuff
NSTDHDL = 10 # frequency of writing to free energy file

 
# Neighbour list update
NSTLIST_EM = 1 # for energy minimization
NSTLIST = 40 # for everything else

# Electrostatics
COLUMB_CUTOFF = 1.2
VDW_CUTOFF = 1.2

LAMBDA_Q = ("""
; init_lambda_state        0    1    2
coul_lambdas             = 0.00 0.5 1.00
""",
2)


FREE_ENERGY = """ ; Free energy control
free_energy              = yes
init_lambda_state        = {lambda_state}
delta_lambda             = 0
calc_lambda_neighbors    = -1        ; All lambdas
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
integrator               = sd       ; Langevin dynamics
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
nstxout-compressed       = 0
"""

NEIGHBOURSEARCH = """ Neighboursearch and short-range nondbonded interactions
cutoff-scheme            = {cutoff}
nstlist                  = {nstlist}
ns_type                  = grid
pbc                      = {pbc}
rlist                    = 1.2
"""

ELECTROSTATICS_AND_WDW_AND_PME = """; Electrostatics
 coulombtype              = {coltype}
rcoulomb                 = {columb_cutoff}
epsilon-rf               = {epsilon_rf}
; van der Waals
vdwtype                  = cutoff
vdw-modifier             = potential-switch
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
tc_grps                  = system
tau-t                    = 1.0
ref-t                    = {temperature}
"""

PRESURE_COUPLING_ON = """; Presure coupling
Pcoupl                   = no   ; Parrinello-Rahman for all MD but NVT
tau_p                    = 1.0
compressibility          = 4.5e-05
ref_p                    = 1.0
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
lincs-order              = 12
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
                          #   MDP Options END #
##############################################################################



MDOUT_FOLDER = './mdout'
#FOLDER_NAMES = ['EM1', 'EM2', 'NVT', 'NPT', 'MD']
#FOLDER_NAMES = ['EM1', 'EM2', 'MD']
FOLDER_NAMES = [ 'MD']

MOL_GRO = 'MOL.gro'
MOL_TOP = 'MOL.top'
MDP_FILES_ROOT = './mdp'
MOL_FILES_DIR = './MOL'
XVG_FILES = './xvg'


def process_command_line(argv):
    """Processes arguments and returns namespace of them."""
    parser = argparse.ArgumentParser(description=""" Run a molecular dynamics
    simulation with gromacs.""")
    #Positional args
    parser.add_argument('crd', help='Amber crd or pdb file', metavar='mol.inpcrd')
    parser.add_argument('top', help='Amber prmtop', metavar='mol.prmtop')
    #Optional args
    parser.add_argument('-n', '--name', help=""" Name of calculation """)
    parser.add_argument('-t', '--temperature', help='Simulation temperature K [298.15]',
                        default='298.15')
    parser.add_argument('--nsteps', help=""" Number of calculation steps [500] """,
                        type=int, default=500)
    parser.add_argument('--time_step', help=""" Time step in fs [2].""",
                        type=float, default=2.0)
    parser.add_argument('-p', '--procs', help="""Number of processores to
                        use during md run. [4]""",
                        default='4')
    parser.add_argument('--box_size', help="""The distance from the molecule
                        to the end of the box in nm [1.6] """, default='1.6',
                        metavar='box_size')
    parser.add_argument('--box_type', help=""" Type of box to simulate in.
                        Chose from {cubic, dodecahedron} [cubic]""",
                        default='cubic')
    parser.add_argument('--fe_type', help=""" Which interactions to decouple.
                        Choose from {q or q_int}.
                        In q_int int stands for integer
                        and means that simulation is going to be run for int
                        number of lambdas equally spread between 0 and 1.
                        [q]""", default='q')
    parser.add_argument('--lambda_state', help="""Evaluate at this lambda [1].""",
                        metavar='lambda',  default='1')
    parser.add_argument('--bond_const', help=""" Bond constraints
                        Chose from {none, h-bonds, all-bonds} [all-bonds]""",
                        default='all-bonds')
    parser.add_argument('--couple_intramol', help=""" Decouple intermolecular
                        interactions in addition to intramolecular. Yes or no
                         [no]. """, default='no'
                        )
    parser.add_argument('--pbc', help=""" Specify pbc. Can be either {no,
                        xyz, xy} [xyz]. """, 
                        default='xyz')
    parser.add_argument('--coltype', help=""" Coulmbtype {Cut-off, PME, ..}
                        For full list of options check mdp manual [PME]. """, 
                        default='pme')
    parser.add_argument('--epsilon_rf', help=""" Dielectric constant outside
                        coulomb reaction field. Default value corresponds to
                        SPC-E water [73.5]. """, 
                        default=73.5, type=float)
    parser.add_argument('--fourierspacing', help=""" Spacing between grid points
                        in PME. Increase for parallel [0.15]. """, 
                        default=0.15, type=float)
    parser.add_argument('--pme_order', help=""" Order of interpolation for 
                        particle mesh ewald [6]. """, 
                        default=6, type=int)
    parser.add_argument('--cutoff', help=""" Cut-off scheme {Verlet, group}.
                        [verlet] """, default='verlet')
    parser.add_argument('--dispcor', help=""" Dispersion correction {no, EnerPres}.
                        [EnerPres] """, default='EnerPres')
    parser.add_argument('--acpype', help="""Path to acpype.py file.
                        [ /home/max/opt/acpype/acpype.py ]""",
                        default='/home/max/opt/acpype/acpype.py',
                        metavar='acpype_path')
                        
    return parser.parse_args(argv)


def calc_fe_paramters(fe_type):
    """ Check the fe type of calc and return lambda parameters """
    if fe_type == 'q':
        interactions_lambda0 = 'vdw-q'
        interactions_lambda1 = 'vdw'
        lambda_vec = LAMBDA_Q[0]
    elif fe_type.startswith('q_') and fe_type[2:].isdigit():
        interactions_lambda0 = 'vdw-q'
        interactions_lambda1 = 'vdw'
        lambda_states = int(fe_type[2:])
        if lambda_states < 2:
            raise ValueError('Number of states should be bigger than 1.')
        coul_lambdas = [round(float(i)/(lambda_states-1), 3) for i in range(lambda_states)]
        lambda_vec = 'coul_lambdas      = ' + ' '.join(map(str, coul_lambdas))
    else:
        raise ValueError('Unknown free energy type (fe_type).')
    return interactions_lambda0, interactions_lambda1, lambda_vec
    

def mdrun_ferrari(mdrun_name, job_name, args):
    """ Run mdrun on Ferrari.
    job_list is a list containing tuples with following data:
    mdrun_name, lambda_state
    """
    p, _ = os.path.split(job_name)
    if not p:
        p = '.'
    run_script_name = '{}_mdrun_script.sh'.format(job_name)
    if run_script_name[0].isdigit():
        # on archie script name cannot start with digit
        run_script_name = 'n' + run_script_name
    run_script_name = os.path.join(p, run_script_name)
    run_strings = ''
    single_run_string = 'mdrun -deffnm {mdrun_name} -s {mdrun_name}.tpr -v -nt {procs}\n'
    run_strings += single_run_string.format(mdrun_name=mdrun_name, procs=args.procs)
    run_script_header = JOB_SCRIPT_HEADER
    run_script = run_script_header + '\n' + run_strings
    with open(run_script_name, 'w') as f:
        f.write(run_script)
    subprocess.call(['bash', run_script_name])


def convet_to_grom(crd, top, acpype):
    """Takes as input amber crd and top files, returns gromacs
    gro and top files."""
    # flag b ensures that the name is MOL
    subprocess.call(['python', acpype, '-p', top, '-x', crd, '-b', 'MOL'])
    return './MOL_GMX.gro', './MOL_GMX.top'


def create_box_and_solvate(gmx_gro, gmx_top, args):
    """Solvate molecule with  3 point water mode."""
    subprocess.call(GMX + [ 'editconf', '-f', gmx_gro, '-c', '-d', args.box_size,
                    '-bt', args.box_type, '-o', 'MOL_BOX.gro'])
    shutil.move('MOL_BOX.gro', MOL_GRO)
    shutil.copy(gmx_top, MOL_TOP)
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


def gen_md_mdp_file(nsteps, dt, args):
    """Create NVT, NPT, and MD dp file."""
    interactions_lambda0, interactions_lambda1, lambda_vec = \
                                                calc_fe_paramters(args.fe_type)
    generic_md_file = MD_MDP.format(run_control=MD_RUN_CONTROL,
                                    output_control=MD_OUTPUT_CONTROL,
                                    neighboursearch=NEIGHBOURSEARCH,
                                    estat_vdw_pme=ELECTROSTATICS_AND_WDW_AND_PME,
                                    tcoupling=TEMPERTURE_COUPLING_ON,
                                    pcoupling=PRESURE_COUPLING_OFF,
                                    free_energy=FREE_ENERGY,
                                    vel_gen=VELOCITY_GENERATION,
                                    bond_constraints=BOND_CONSTRAINT)
    md_file = generic_md_file.format(nsteps=nsteps,
                                      dt=dt,
                                      nstlist=NSTLIST,
                                      columb_cutoff=COLUMB_CUTOFF,
                                      vdw_cutoff=VDW_CUTOFF,
                                      temperature=args.temperature,
                                      pbc=args.pbc,
                                      coltype=args.coltype,
                                      dispcor=args.dispcor,
                                      cutoff=args.cutoff,
                                      epsilon_rf=args.epsilon_rf,
                                      gen_vel='yes',
                                      couple_intramol=args.couple_intramol,
                                      bond_contraint=args.bond_const,
                                      continuing_run='no',
                                      fourierspacing=args.fourierspacing,
                                      pme_order=args.pme_order,
                                      nstdhdl=NSTDHDL,
                                      lambda_state=args.lambda_state,
                                      interactions_lambda0 = interactions_lambda0,
                                      interactions_lambda1 = interactions_lambda1,
                                      lambda_vec = lambda_vec 
)
    return md_file


def generate_mdp_file(lambda_state, nsteps, dt, args):
    """Create a tree of mdp file:
    MDP_FILE_ROOT/
    |---MD
        |---md_{lambda_state}.mdp
    """
    mdp_path = os.path.join(MDP_FILES_ROOT, 
                    'md_{lambda_state}.mdp'.format(lambda_state=lambda_state))
    try_to_make_dir(MDP_FILES_ROOT)
    mdp = gen_md_mdp_file(nsteps, dt, args)
    with open(mdp_path, 'wb') as f:
        f.write(mdp)
    return mdp_path


def run_md_simulations(job_name, lambda_state,
                       gro_file, top_file, mdp_file,
                       args):
    """Preforms all MD simulation (takes quite a while)"""
    log_file = open('md_{}.log'.format(job_name), 'wb')

    print 'Starting Lambda state: {}'.format(lambda_state)
    log_file.write("""\n##########################\n\n\n\n
                   Lambda state {}:""".format(lambda_state))
    dirname = 'MD'
    mdrun_name = os.path.join(dirname,
                            '{}.{}'.format(dirname.lower(), lambda_state))
    try_to_make_dir(dirname)
    mdout_mdp_name = os.path.join(MDOUT_FOLDER, job_name + '_{}'.format(lambda_state))
    grompp_call = GMX + ['grompp',
                   '-f', mdp_file,
                   '-po', mdout_mdp_name,
                   '-c', gro_file,
                   '-p', top_file,
                   '-o', mdrun_name + '.tpr']
    grompp_out = subprocess.check_output(grompp_call)
    log_file.write('\n\n\nGrompp output:')
    log_file.write(grompp_out)
    mdrun_ferrari(mdrun_name, job_name, args)
    log_file.close()


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
    if args.name:
        job_name = args.name
    else:
        job_name = args.top[:-7]
    lambda_state = args.lambda_state
    print 'Starting job {}'.format(job_name)
    print 'Simulation will be run for {} state'.format(lambda_state)
    try_to_make_dir(job_name)
    incrd_f = os.path.join(job_name, args.crd)
    prmtop_f = os.path.join(job_name, args.top)
    shutil.copy(args.crd, incrd_f)
    shutil.copy(args.top, prmtop_f)
    os.chdir(job_name)
    try_to_make_dir(MDOUT_FOLDER)    
    return job_name, lambda_state
    

def main(argv):
    args = process_command_line(argv)
    check_input(args)
    job_name, lambda_state = prepare_dir(args)
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
        args.crd = new_crdname
    # md paramters
    dt = args.time_step
    nsteps = args.nsteps
    dt = dt/1000   # convert to ps
    # prepare and run calc
    gmx_gro, gmx_top = convet_to_grom(args.crd, args.top, args.acpype)
    gro_file, top_file = create_box_and_solvate(gmx_gro, gmx_top, args)
    gro_file, top_file = clean_gro_and_top_files(gro_file, top_file,
                                                 args.crd, args.top)
    mdp_file = generate_mdp_file(lambda_state,
                                   nsteps, dt, args)
    run_md_simulations(job_name, lambda_state,
                       gro_file, top_file, mdp_file,
                       args)
    

if __name__ == '__main__':
    main(sys.argv[1:])

