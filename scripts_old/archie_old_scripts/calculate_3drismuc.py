import os,sys
#########################################################
#
# - parse output of rism3d.snglpnt
# - calculate hydration free energy using 3DRISM/UC
# - write output to results.txt
#
#########################################################

#Read data
for line in open("rism3d.log", "r"):
  if line[0:11]=="rism_exchem":
    kh=line.split()[1] 
  if line[0:11]=="rism_exchGF":
    gf=line.split()[1]
  if line[0:11]=="rism_volume":
    v=line.split()[1]

#Calculate dGhyd using 3DRISM/UC
uc=float(gf) - 3.2217 * 0.0333 * float(v) + 0.5783

#Write output to results.txt
output=["dGhyd(KH)= "+kh+" kcal/mol", "dGhyd(GF)= "+gf+" kcal/mol", "PMV= "+v+" AA^3", "dGhyd(UC)= "+str(uc)+" kcal/mol"]
outfile=open("results.txt", "w")
for line in output:
  outfile.write(line+"\n")
outfile.close()



