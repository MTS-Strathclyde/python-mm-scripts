#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 18:33:02 2013

@author: a92549
"""

import sys
import fnmatch
import os


def file_exist(path):
    try:
        with open(path) as f:
            return True
    except IOError:
       return False



def main(argv):
    for name in argv:
        with open(name, 'rb') as f:
            lines = f.readlines()
        if file_exist(name[:-4]):
            for file in os.listdir('.'):
                if fnmatch.fnmatch(file, name[:-4] + '.e*'):
                    print name + " Job run out of resources!"
                    break
            else:
                print name + " is still running"
        elif lines[-1][:19] == ' Normal termination':
            print name + ' has finished'
        elif lines[-9][:42] == ' >>>>>>>>>> Convergence criterion not met.':
            print name + ' didnt converge after n cycles'
        elif lines[-4] == ' Error termination request processed by link 9999.\n':
            print name + ' didnt converge (9999).'
        else:
            print name + ' didnt finish.'

    
    
if __name__ == '__main__':
    main(sys.argv[1:])

