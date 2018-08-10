#!/usr/bin/python

"""Draw mopout molecule"""

import pybel
import sys
#import argparse

def main(argv):
    mol = pybel.readfile("mopout", argv.pop()).next()
    mol.draw()

if __name__ == '__main__':
    main(sys.argv[1:])
