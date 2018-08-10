#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri May 16 09:39:47 2014

@author: max
"""

import sys

def is_float(number):
    try:
        float(number)
        return True
    except ValueError:
        return False


def extract_charges(name):
    """Extracts lowrin charges from GAMESS output file.
    Returns them as a list of floats."""
    with open(name, 'rb') as lines:
        charges = []
        for line in lines:
            if line.startswith('          TOTAL MULLIKEN AND LOWDIN ATOMIC'):
                lines.next()
                for chg_line in lines:
                    chg_line = chg_line.split()
                    if len(chg_line) == 6 and is_float(chg_line[-1]):
                        charges.append(float(chg_line[-1]))
                    else:
                        break
                break
    return charges


def main(name):
    charges = extract_charges(name)
    for chg in charges:
        print chg

if __name__ == '__main__':
    main(sys.argv[1])