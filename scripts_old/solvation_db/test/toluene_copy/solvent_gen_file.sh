#!/bin/csh -f

cat > SPC.inp <<EOF
&PARAMETERS
THEORY='DRISM', CLOSUR='KH',   !Theory
NR=16384, DR=0.025,            !Grid
OUTLST='xCGT', routup=384, toutup=0,  !Output
NIS=20, DELVV=0.3, TOLVV=1.e-12, !MDIIS
KSAVE=-1, KSHOW=1, maxste=10000,  !iter
SMEAR=1, ADBCOR=0.5,              !Electrostatics
!bulk solvent properties
TEMPER=298, DIEps=78.497,
NSP=1
/
&SPECIES
DENSITY=55.296d0,
MODEL="$AMBERHOME/dat/rism1d/model/SPC.mdl"
/
EOF

rism1d SPC > SPC.out || goto error

