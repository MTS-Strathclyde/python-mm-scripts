#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 17:07:05 2013

@author: a92549

Makes json dictionaries from basis set text files using gaussian classes.
"""

from comp_chem import gaussian
import sys


def main(argv):
    bs = gaussian.CustomBasis(argv[0])
    bs.write_json(argv[0] + '.json')
    try:
        ecp = gaussian.ECPBasis(argv[0] + '_pseudo')
        ecp.write_json(argv[0] + '_pseudo.json')
    except IOError:
        print 'Couldn\'t find ECP basis'
    
    
if __name__ == '__main__':
    main(sys.argv[1:])
    