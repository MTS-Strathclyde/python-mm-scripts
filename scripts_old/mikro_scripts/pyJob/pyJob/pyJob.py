#!/usr/bin/python
'''
Created on Mar 17, 2012

@author: Max
'''

import sys
from Molecule import Molecule
from Conformers import conformers

if __name__ == '__main__':
    #print("filename, min_dist, max_dist")
    filename = sys.argv[1]
    min_dist = float(sys.argv[2])
    max_dist = float(sys.argv[3])
    molecula = Molecule(filename, min_dist, max_dist)
    molecula.edit()
    accepted_filenames_list = (conformers(molecula))
    
    if len(accepted_filenames_list) == 0:
        print "None of the conformers satisfies given distance constraints"
    else:
        print "creating script"
        out = "#!/bin/sh\n"
        for name in accepted_filenames_list:
            out += "mopac " + str(name) + "\n"
        f = open("pyJob.sh", 'w')
        f.write(out)
        f.close()
