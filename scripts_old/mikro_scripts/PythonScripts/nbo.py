#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Mon May 13 00:17:19 2013

@author: a92549
"""

import argparse
import shutil
import sys

STRING = """%nproc=8
%mem=8GB
%chk={0}_{1}
# Guess=Read Geom=AllCheck ChkBasis pop=nbo


"""



def process_command_line(argv):
    """
    Processes arguments and returns namespace of them
    """
    parser = argparse.ArgumentParser(description="""NBO job.""")
    #Positional args
    parser.add_argument('molecules', metavar='mol.chk',
                        help="""Computation to be resumed.""",
                        nargs='+')
    #Optional args
    parser.add_argument('-c', '--chk_ext',
                        help="""Checkpoint file extenison. (chk)""",
                        default='.chk')
    parser.add_argument('-s', '--suffix',
                        help="""Suffix, which will be added to the molecule
                        name end before extension (nbo)""",
                        default='nbo')

    return parser.parse_args(argv)    

    
def main(argv):
    args = process_command_line(argv)
    for mol in args.molecules:
        no_ext_name = mol[:-4]
        txt = STRING.format(no_ext_name,args.suffix)
        with open(no_ext_name + '_' + args.suffix + '.com' , 'wb') as f:
            f.write(txt)
        shutil.copy(no_ext_name + args.chk_ext, no_ext_name + '_' + args.suffix + args.chk_ext)
    
    
if __name__ == '__main__':
    main(sys.argv[1:])