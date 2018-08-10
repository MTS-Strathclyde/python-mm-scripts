#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 17:07:05 2013

@author: a92549

Makes json dictionaries from basis set text files using gaussian classes.
"""

import gaussian
import sys


def main(argv):
    bs = gaussian.CustomBasis(argv[0])
    bs.write_json(argv[0] + '.json')
    ecp = gaussian.ECPBasis(argv[0] + '_pseudo')
    ecp.write_json(argv[0] + '_pseudo.json')
    
    
if __name__ == '__main__':
    main(sys.argv[1:])
    