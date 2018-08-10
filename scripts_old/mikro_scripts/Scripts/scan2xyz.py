#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  2 17:01:53 2012

@author: mishin1991

Converts TINKER scan output files to babel xyz files

"""

import sys
import os


def main(args):
    """accepts list of tinke scan output file=s as argument,
    tries to convert each of them into xyz file"""
    for file_ in args:
        #open files
        tink_f = open(file_, 'r')
        xyz_f = open(file_ + '.xyz', 'w')
        try:
            #read and write file head
            num_of_atoms = tink_f.next().split()[0]
            xyz_f.write(num_of_atoms + '\n' + file_ + '\n')
            #read and write file main body
            lines = []
            for line in tink_f:
                line_list = line.split()
                atom = line_list[1][0]
                converted_line = atom + ' ' +line_list[2] + ' ' +  line_list[3] + \
                                 ' ' + line_list[4] + '\n'
                lines.append(converted_line)
            xyz_f.writelines(lines)
            xyz_f.close()
            tink_f.close()
        except IndexError:
            print "couldn't convert file " + file_
    
if __name__ == '__main__':
    main(sys.argv[1:])
            