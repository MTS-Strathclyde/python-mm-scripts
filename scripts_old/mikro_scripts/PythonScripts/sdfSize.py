#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  6 19:36:25 2012

@author: mishin1991
"""

import pybel
import sys

def size(filename):
    mol = pybel.readfile('sdf', filename)
    return sum(1 for x in mol)

def main(argv):
    for filename in argv:
        print filename + ' length is: ' + str(size(filename))
        
if __name__ == '__main__':
    main(sys.argv[1:])