#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  3 02:52:15 2013

@author: a92549
"""

import sys
import argparse

def process_command_line(argv):
    """
    Processes arguments and returns namespace of them
    """
    parser = argparse.ArgumentParser(description="""Makes supplied line first 
                                    row of all given file.""")
    #Positional args
    parser.add_argument('files', metavar='<f.ext>',
                        help="""All files, to wich line will be added""", nargs='+')
    #Optional args
    parser.add_argument('-l', '--line', help="""Surround by quotes, if it has special chars. (' ')""",
                        default=" ")
    return parser.parse_args(argv)


def main(argv):
    args = process_command_line(argv)
    for f_name in args.files:
        with open(f_name, 'rb') as f:
            f_txt = f.read()
        new_txt = args.line + '\n' + f_txt
        with open(f_name, 'wb') as f:
            f.write(new_txt)
            
if __name__ == '__main__':
    main(sys.argv[1:])
    