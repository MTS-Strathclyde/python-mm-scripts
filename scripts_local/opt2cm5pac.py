#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 11 13:43:55 2015

@author: max
"""

import sys

def main(argv):
    if ('-h' or '--help') in argv:
        print 'Convert Gaussian opt to format suitable for cm5pac'
        return 0
    for gauopt in argv:
        if not gauopt.endswith('.log'):
            print 'Oputput file should have extension log'
            print 'Skipping {}'.format(gauopt)
            continue
        try:
            lines = []
            natoms_line_number = 0
            last_input_orientation = 0
            with open(gauopt) as f:
                for i, l in enumerate(f):
                    lines.append(l)
                    if 'NAtoms=' in l:
                        natoms_line_number = i
                    if 'Input orientation:' in l:
                        last_input_orientation = i
            pruned_file = []
            if natoms_line_number and last_input_orientation:
                pruned_file.append(lines[natoms_line_number])
                pruned_file.extend(lines[last_input_orientation:])
            with open(gauopt[:-4] + '_tail.log', 'w') as f:
                f.writelines(pruned_file)
        except IOError as e:
            print e
            print 'Skipping {}'.format(gauopt)
            
            
if __name__ == '__main__':
    main(sys.argv[1:])
                
                
            
            