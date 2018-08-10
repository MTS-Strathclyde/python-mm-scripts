#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 25 14:39:30 2013

@author: a92549
"""

import sys

def main(argv):
    for log in argv:
        with open(log, 'rb') as f:
            lines = f.readlines()
        for line in lines:
           # print line
            if line.startswith(' Alpha  occ. eigenvalues --'):
                try:
                    HOMO = float(line.split()[-1])
                except ValueError:
                    pass
            elif line.startswith(' Alpha virt. eigenvalues --'):
                LUMO = float(line[28:].split()[0])
                break
        print log
        print 'HOMO ' + str(HOMO)
        print 'LUMO ' + str(LUMO)
        print 'dif. ' + str(LUMO-HOMO)
        
        
if __name__ == '__main__':
    main(sys.argv[1:])