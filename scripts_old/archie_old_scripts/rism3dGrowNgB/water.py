# -*- coding: utf-8 -*-
"""
Created on Fri May 23 09:40:18 2014

@author: max
"""

SOLV_SUCEPT_SCRPT = """#!/bin/csh -f

cat > water_{temp}.inp <<EOF
&PARAMETERS
	THEORY='DRISM', CLOSUR='KH',           !Theory
	NR=16384, DR=0.025,                    !Grid
	OUTLST='xCGT', routup=384, toutup=0,   !Output
	NIS=20, DELVV=0.3, TOLVV=1.e-12,       !MDIIS
	KSAVE=-1, KSHOW=1, maxstep=10000,      !Check pointing and iterations
	SMEAR=1, ADBCOR=0.5,                   !Electrostatics
	TEMPER={temp}, DIEps={diel},           !bulk solvent properties
	NSP=1
/
	&SPECIES                               !SPC water
	DENSITY={conc}d0,
	MODEL="$AMBERHOME/dat/rism1d/model/SPC.mdl"
/
EOF

rism1d water_{temp} > water_{temp}.out || goto error

"""




def dielectric_const(T):
    """Return water dielectric constant for temperature 253.15K < T < 383.15K.
    
    Uses correlating equation (eq. 9) for static dielectri constant found in
    the doucment by The International Association for the Properties of 
    Water and Steam from 2011
    (http://www.iapws.org/relguide/LiquidWater.pdf)
    Pressure = 0.1 MPa
    
    >>> round(water_dielectric_const(273.15), 3)
    87.927
    >>> round(water_dielectric_const(298.15), 3)
    78.375
    >>> round(water_dielectric_const(375), 3)
    55.266
    """
    T_star = T/300.0
    coefs = [-43.7527, 299.504, -399.364, 221.327]
    exp_f = [-0.05, -1.47, -2.11, -2.31]
    e = 0
    for i in range(4):
        e += coefs[i]*T_star**(exp_f[i])
    return e
    
    
def concentration(T):
    """Return water concentration for temperature range 253.15K < T < 383.15K.
    
    Uses correlating equation (eq. 2) for specific volume found in
    the doucment by The International Association for the Properties of 
    Water and Steam from 2011
    (http://www.iapws.org/relguide/LiquidWater.pdf)
    Pressure = 0.1 MPa    
    
    >>> round(water_concentration(273.15), 3)
    55.498
    >>> round(water_concentration(298.15), 3)
    55.343
    """
    p0 = 10.0**5    # Pa
    R = 8.31464     # J/mol/K
    Tr = 10.0
    Ta = 593.0
    Tb = 232.0
    a = [1.93763157E-2,
         6.74458446E+3,
        -2.22521604E+5,
         1.00231247E+8,
        -1.63552118E+9,
         8.32299658E+9]
    b = [5.78545292E-3,
        -1.53195665E-2,
         3.11337859E-2,
        -4.23546241E-2,
         3.38713507E-2,
        -1.19946761E-2]
    n = [None, 4., 5., 7., 8., 9.]
    m = [1., 2., 3., 4., 5., 6.]
    def alpha(T): 
        return Tr/(Ta - T)
    def beta(T):
        return Tr/(T - Tb)
    coef = a[0] + b[0]*beta(T)**m[0]
    for i in range(1, 6):
        coef += a[i]*alpha(T)**n[i] + b[i]*beta(T)**m[i]
    v0 = R*Tr/p0*coef  # m3/mol
    return 1/(v0*1000)    # mol/L


def water_succ_rism_script(T):
    """Create water susceptibility script for given temperature T in K."""
    diel = round(dielectric_const(T), 3)
    conc = round(concentration(T), 3)
    return SOLV_SUCEPT_SCRPT.format(temp=T, diel=diel, conc=conc)


