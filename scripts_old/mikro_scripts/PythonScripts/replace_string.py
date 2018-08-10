#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 11 01:10:52 2013

@author: mishin1991
"""
import os
import sys
import argparse

def process_command_line(argv):
    """
    Processes arguments and returns namespace of them
    """
    parser = argparse.ArgumentParser(description="""Replaces old substring with
                            new one in all files in given list. Writes modified
                            files to the same folder and adds to their name 
                            new substring.""")
    #Positional args
    parser.add_argument('input', help='Input files', nargs='+')
    #Optional args
    parser.add_argument('-o', '--old', help="""Old string. Case-sensitive.""")
    parser.add_argument('-n', '--new', help="""New string. Case-sensitive.""")
    parser.add_argument('-t', '--times', help="""Times substitution happen""")
    return parser.parse_args(argv)

def main(argv):
    args = process_command_line(argv)
    for file_ in args.input:
        #read
        i = open(file_, 'rb')
        txt = i.read()
        i.close
        #substitute
        if args.times:
            txt = txt.replace(args.old, args.new, int(args.times))
        else:
            txt = txt.replace(args.old, args.new)
        #write
        name, extension = os.path.splitext(file_)
        new_name = name + '_' + args.new + extension
        o = open(new_name, 'wb')
        o.write(txt)    
        o.close()
    
if __name__ == '__main__':
    main(sys.argv[1:])