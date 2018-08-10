#!/usr/bin/python2

import sys
import argparse
import math
import os
import pandas as pd

epsilon0 =  8.854187817e-12  #F/m
Na = 6.02214129e23
k = 1.9872041E-3  # kcal/mol/K


YAWS_PATH = '/home/max/Documents/PhD/projects/nonaqueous_solvation/\
new_attempt/_solvent_data/all_data_yaws2014.csv'



sh_txt = """#!/bin/csh -f

cat > {name}.inp <<EOF
&PARAMETERS
        THEORY='DRISM', CLOSUR='$1',           !Theory
        NR=16384, DR=0.025,                    !Grid
        OUTLST='xCGT',
        NIS=20, DELVV=0.3, TOLVV=1e-12,
        KSAVE=-1, KSHOW=1, maxstep=1000,      !Check pointing and iterations
        SMEAR=1, ADBCOR=0.5,                   !Electrostatics
        TEMPER={T}, DIEps=1.0,              !bulk solvent properties
        NSP=1
/
        &SPECIES                               !solvent
        units='1/A^3',
        DENSITY={dens:.5e},
        MODEL="{name}.mdl"
/
EOF

rism1d {name} > {name}.out
"""

mdl_txt = """%VERSION  VERSION_STAMP = V0001.000  DATE = 01/22/15  11:58:10
%FLAG TITLE
%FORMAT(20a4)
{name}
%FLAG POINTERS
%FORMAT(10I8)
       1       1
%FLAG ATMTYP
%FORMAT(10I8)
       1
%FLAG ATMNAME
%FORMAT(20a4)
{name:.4}
%FLAG MASS
%FORMAT(5e16.8)
 {mass: .8E}
%FLAG CHG
%FORMAT(5e16.8)
  0.00000000e+00
%FLAG LJEPSILON
%FORMAT(5e16.8)
 {eps: .8E}
%FLAG LJSIGMA
%FORMAT(5e16.8)
 {sig: .8E}
%FLAG MULTI
%FORMAT(10I8)
       1
%FLAG COORD
%FORMAT(5e16.8)
  0.00000000e+00  0.00000000e+00  0.00000000e+00
"""

# name : (  Mr,  Tc,   Vc,   rho_298.15)
# units   g/mol, K, cm^3/mol, g/cm^3
#solvent_dic = {
#'acetonitrile' : (41.05, 545.5, 173, 0.779),
#'benzene' : (78.11, 562.1, 256, 0.873),
#'carbtet' : (153.82, 556.35, 276, 1.613),
#'chloroform' : (119.38, 536.4, 239, 1.48),
#'cs2' : (76.139, 552.0, 160, 1.256),
#'decane' : (142.29, 617.7, 600, .728),
#'dichloroethane' : (84.93, 561.6, 220, 1.318),
#'dmso' : (78.13, 729.0, 227, 1.095),
#'heptane' : (100.21, 540.2, 428, .682),
#'toluene' : (92.14, 591.75, 316, .865),
#'diethylether': (74.12, 466.70, 280.0, 0.78 )
#}


# name : (eps_mult, v_mult)
model_multipliers = {
'A' : (1.2593, 0.302), # ref: galliero 2005
'B' : (1.2593, 0.3189), # ref Chung 1984
'C' : (1.313,  0.304)  # ref Okumura 2000
}

solvent_synonyms = {
'carbtet' : 'carbon tetrachloride',
'cs2' : 'carbon disulfide',
'dichloroethane' : '1,2-dichloroethane',
'diethylether' : 'diethyl ether',
'dmso' : 'dimethyl sulfoxide',
'decalin' : 'trans-decahydronaphthalene',
'ethylacetate' : 'ethyl acetate',
'isooctane' : '2,2,4-trimethylpentane',
'octanol' : '1-octanol',
'xylene' : 'xylenes',
#'hexadecane' : 'n-hexadecane'
}

def process_command_line(argv):
    """
    Processes arguments and returns namespace of them
    """
    parser = argparse.ArgumentParser(description="""Create LJ solvent.""")
    parser.add_argument('-p', '--path',
                        help="""Output path [./]""", default='./')
    parser.add_argument('--noout',
                        help="""Don't write output""", action='store_true')
    # Realistic solvents
    real_solv = parser.add_argument_group('Real solvents',
                            'Options to create a model of realistic solvents.')
    real_solv.add_argument('-m', '--model',
                        help="""Type of the model. Can be either A, B or C.
                        A: galliero 2005, B: Chung 1984, C: Okumura 2000 [C].""",
                        default='C')
    real_solv.add_argument('-s', '--solvent',
                        help="""Name of the solvent matching one of the names
                        in the following document: {}. The following properties
                        will be read: name, mass, crit temp, crit colume. You 
                        can overwrite any of these properties using
                        corresponding flag""".format(YAWS_PATH))
    real_solv.add_argument('-t', '--temperature',
                        help="""Temperature at which density is computed [298.15]""",
                        type=float, default=298.15)
    real_solv.add_argument('-n', '--name',
                        help="""Name of the model.""")
    real_solv.add_argument('-mr',
                        help="""Molar mass.""", type=float)
    real_solv.add_argument('-tc',
                        help="""Critical temperature. Alternaticely, one can
                        provide epsilon""", type=float)
    real_solv.add_argument('--eps',
                        help=""" Epsilon [kcal/mol]. Specifying this option
                        will force script to ignore critical temperatrue""", type=float)
    real_solv.add_argument('-vc',
                        help="""Critical volume [cm^3/mol]. Alternatively, one can define
                        sigma""", type=float)
    real_solv.add_argument('--sigma',
                        help=""" Sigma [A]. Specifying this option
                        will force script to ignore critical volume""", type=float)
    real_solv.add_argument('-d', '--rho', 
                        help=""" Density [g/cm^3].""", type=float)
    # LJ solvent
    lj_solv = parser.add_argument_group('Lennard-Jones solvent',
                        """Reduced units to define lj solvent. Defining any 
                        of these options will make script ignore realistic 
                        solvation parameters""")
    lj_solv.add_argument('-tr', help="""Reduced temperature.""",
                        type=float)
    lj_solv.add_argument('-dr', help="""Reduced density.""",
                        type=float)
    return parser.parse_args(argv)


def compute_yaws_density(yaws_row, t=298.15):
    """Returns g/cm^3"""
#    if float(yaws_row.T_min) > t or float(yaws_row.T_max) < t:
#        raise ValueError('Extrapolation equation can not be used at t={}'.format(t))
    A = float(yaws_row.A)
    B = float(yaws_row.B)
    C = float(yaws_row.C)
    n = float(yaws_row.n)
    #print A,B,C,n
    density = A*B**(-(1.-t/C)**n)
    return density


def main(argv):
    args = process_command_line(argv)
    T = args.temperature
    if args.tr or args.dr:
        ## LJ-reduced solvent
        name = 'LJs'
        mass = 20.0   # g/mol
        rho_A = 7.5e-3  # 1/A^3
        sigma = (args.dr/rho_A)**(1./3)  # A
        eps = k*T/args.tr  # kCal/mol
    else:
        if args.solvent:
            df = pd.read_csv(YAWS_PATH)
            # get solvent
            solvent = solvent_synonyms.get(args.solvent, args.solvent)
            if not any(df.Name.isin([solvent])):
                raise ValueError('Unknown solvent. Stopped')
            yaws_row = df[df.Name==solvent]
            mass = float(yaws_row.mass)  # g/mol
            Tc = float(yaws_row.crit_t)    # K
            Vc = float(yaws_row.crit_vol)    # cm^3/mol
            rho = compute_yaws_density(yaws_row)   # g/cm^3
            name = args.solvent
        else:
            if not (args.tc and args.vc and args.rho):
                raise ValueError('Neither solvent name or critical properties \
have been passed. Quiting.')
        # check if we need to overwrite defaults
        if args.name:
            name = args.name
        if args.mr:
            mass = args.mr # g/mol
        if args.tc:
            Tc = args.tc # K
        if args.vc:
            Vc = args.vc # A
        if args.rho:
            rho = args.rho # g/cm^3
        rho_A = rho/mass*Na/1.0e24 # 1/A^3
        eps_mult, v_mult  = model_multipliers[args.model.upper()]
        # compute lj parameters
        if args.sigma:
            sigma = args.sigma # A
        else:
            sigma = (v_mult*Vc/Na*1.0e24)**(1./3)  # A
        if args.eps:
            eps = args.eps # kcal/mol
        else:
            eps = k*Tc/eps_mult  # kcal/mol
        tr = k*T/eps
        dr = sigma**3*rho_A
        print 'Tr:    {:.3f} '.format(tr)
        print 'Dr:    {:.3f} '.format(dr)
#        dipole = float(args[4])
#        F = 1 + dipole**4/(12*k*T*eps*sigma**6*(4*math.pi*epsilon0)**2)*1.0e60*N**2*(3.33564e-30)**4/4184.**2
#        eps = eps*F**2
#        sigma = sigma/F**(1./6)
#        print F
    print 'Eps:   {:.3f} [kcal/mol]'.format(eps)
    print 'Sigma: {:.3f} [A]'.format(sigma)
    if not args.noout:
        with open(os.path.join(args.path, name) + '.mdl', 'w') as f:
            f.write(mdl_txt.format(name=name, mass=mass, eps=eps, 
                                  sig=sigma/2*2**(1./6)))
        with open(os.path.join(args.path, name) + '.sh', 'w') as f:
            f.write(sh_txt.format(name=name, T=T,dens=rho_A))
    

if __name__ == '__main__':
    main(sys.argv[1:])


