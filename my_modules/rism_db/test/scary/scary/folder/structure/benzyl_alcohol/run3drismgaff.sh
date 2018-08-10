#!/bin/bash
#################################################################
#
# Run 3DRISM calculations
#
#
# Usage:
#	bash run3drismgaff.sh [mol]
#
# where:
#	[mol] is the name of a pdb file containing the solute (without the file extension)"
#
#################################################################

#Check input
if [ -z "$1" ];
then
  echo
  echo "ERROR! Correct syntax ./run3drism_tleap.sh [mol]"
  echo " ... where [mol] is the name of a pdb file containing the solute (without the file extension)"
  echo 
  exit
fi

###########################################
#
# Step 1: Assign forcefield parameters 
#         to solute
#
###########################################

antechamber -i $1.pdb -fi pdb -o $1.prepin -fo prepi -c bcc -s 2 -nc 0 -m 1
parmchk -i $1.prepin -f prepi -o $1.frcmod

cat > runleap.in <<EOF
source leaprc.gaff 
loadamberprep $1.prepin
check MOL
loadamberparams $1.frcmod
SaveAmberParm MOL $1.prmtop $1.incrd
SavePdb MOL $1.pdb
quit
EOF

#Run tleap
tleap -f runleap.in

###########################################
#
# Step 2: Solve 3D RISM equations for 
#         solute/solvent system
#	  (using pre-calculated 1D RISM output)
#
###########################################

#Run 3D RISM (Parameters from JPCM 2010)
rism3d.snglpnt --pdb $1.pdb --prmtop $1.prmtop --closure kh --guv g_$1 --cuv c_$1 --xvv wat.xvv --buffer 30.000000 --grdspc 0.300000,0.300000,0.300000 --tolerance 0.0000000001 | tee rism3d.log

###########################################
#
# Step 3: Calculate hydration free energy
#
###########################################

#Parse output file and calculate hydration free energy using 3DRISM/UC functional
python calculate_3drismuc.py


