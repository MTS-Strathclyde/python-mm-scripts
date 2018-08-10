#!/usr/bin/env python
# Andrey I. Frolov, December 2010, MPI MIS, Leipzig
#
# The script converts *.xvv file format into *.dat.
# *.dx contains 1D solvent susceptibility functions in k-space (output from rism1d in AmberTools)
# *.dat: 1st column - support, 2nd-Xth - susceptibility functions .
# Mind: standard Amber output order - has to be clarified.

# USAGE: xvv2dat.py mol.xvv > mol.dat

import sys

fname=sys.argv[1]
ifile=open(fname,"r")
lines=ifile.readlines()
ifile.close()

start=False

data = []

print len(lines)
for ind in range(len(lines)):
    line=lines[ind]
    ln=line.split()
    if (len(ln)<=1): 
        continue
    if (ln[1] == 'POINTERS'):
        #line2 = lines()
        line2 = lines[ind+2]
        ln=line2.split()
        NGrid = int(ln[0])
        NSites = int(ln[1])
    if (ln[1] == 'ATOM_NAME'):
        #line2 = ifile.readline()
        line2 = lines[ind+2]
        l = len(line2[:-1])
        AtomName = []
        i=-1
        print l
        for j in range(0,l,4):
             print i,j
             i+=1
             AtomName.append(line2[j:j+4])
        print AtomName
        
    if ln[1] == 'XVV' and len(ln) == 2:
        xvv_ind = ind + 2
        while xvv_ind < len(lines) and not lines[xvv_ind].startswith('%'):
            data.extend(lines[xvv_ind].split())
            xvv_ind += 1
        break

print len(data)
#xvv=numpy.zeros([NGrid,NSites*NSites+1])
xvv_ind = 0
#print len(data), NGrid*NSites*NSites 
xvv=[]
print NSites*NSites
for i in range(NSites*NSites):
    xvv.append([])
    for j in range(NGrid):
        xvv[i].append(data[xvv_ind])
        xvv_ind += 1
        
cnt = 0
for j in range(NGrid):
    cnt += 1
    for i in range(NSites*NSites):
        strtmp = xvv[i][j] + ' '
    strtmp = str(cnt) + ' ' + strtmp
    print strtmp
    


### EXAMPLE XVV FILE
        
#%VERSION  VERSION_STAMP = V0000.001  DATE = 07:04:10 17:20:44
#%COMMENT NR,NATV,NSP
#%FORMAT(10I8)
#   16384       4       3
#%COMMENT TEMPER,DIEPS,XAPPA,XIKT,DR,SMEAR
#%FLAG THERMO
#%FORMAT(5E16.8)
#  0.30000000E+03  0.78497000E+02  0.17949689E+00  0.18825661E+01  0.25000000E-01
#  0.10000000E+01
#%FLAG ATOM_NAME
#%FORMAT(20a4)
#O   H1  Cs_cCl_c
#%FLAG MTV
#%FORMAT(10I8)
#       1       2       1       1
#%FLAG RHOV
#%FORMAT(5E16.8)
#  0.33300007E-01  0.66600014E-01  0.18066410E-03  0.18066410E-03
#%FLAG QV
#%FORMAT(5E16.8)
# -0.83400000E+00  0.41700000E+00  0.10000000E+01 -0.10000000E+01
#%FLAG QSPV
#%FORMAT(5E16.8)
#  0.00000000E+00  0.00000000E+00  0.10000000E+01 -0.10000000E+01
#%FLAG EPSV
#%FORMAT(5E16.8)
#  0.63596800E+03  0.63596800E+02  0.17009608E+04  0.14891274E+03
#%FLAG SIGV
#%FORMAT(5E16.8)
#  0.17683000E+01  0.69387933E+00  0.19760000E+01  0.25130000E+01
#%FLAG DELHV0
#%FORMAT(5E16.8)
#  0.28590226E+02  0.28585391E+02 -0.16401054E+04 -0.16395915E+04
#%COMMENT COLUMN MAJOR NR X NAT X NAT
#%FLAG XVV
#%FORMAT(5E16.8)
#  0.90277788E-01  0.90275959E-01  0.90270903E-01  0.90262663E-01  0.90251282E-01
#  0.90236728E-01  0.90218898E-01  0.90197642E-01  0.90172778E-01  0.90144104E-01
#  0.90111420E-01  0.90074540E-01  0.90033302E-01  0.89987576E-01  0.89937277E-01
#  0.89882359E-01  0.89822826E-01  0.89758726E-01  0.89690153E-01  0.89617239E-01
#  0.89540157E-01  0.89459111E-01  0.89374335E-01  0.89286087E-01  0.89194646E-01
#  0.89100306E-01  0.89003373E-01  0.88904161E-01  0.88802990E-01  0.88700184E-01