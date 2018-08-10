#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 12 16:38:15 2013

@author: a92549


Fixes lack of / between tzvp and tzvpfit

"""

import sys

def main(argv):
    for com in argv:
        with open(com, 'rb') as f:
            txt = f.read()
        if 'tzvp tzvpfit' in txt:
            parts = txt.split('tzvp tzvpfit',1)
            new_txt = parts[0] + 'tzvp/tzvpfit' + parts[1]
            with open(com, 'wb') as f:
                f.write(new_txt)
        elif 'tzvp\ntzvpfit' in txt:
            parts = txt.split('tzvp\ntzvpfit',1)
            new_txt = parts[0] + 'tzvp/tzvpfit' + parts[1]
            with open(com, 'wb') as f:
                f.write(new_txt)
    
    
if __name__ == '__main__':
    main(sys.argv[1:])