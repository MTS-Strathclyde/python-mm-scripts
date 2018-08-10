#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 28 19:23:41 2012

@author: a92549
"""

import sipsik
import sys
import argparse


def process_command_line(argv):
    """
    Processes arguments and returns namespace of them
    """
    parser = argparse.ArgumentParser(description="""Runs given command on sipsik.""")
    #Optional args
    parser.add_argument('-c', '--command', help="""Command to run (qstat)""",
                        default="qstat")
    return parser.parse_args(argv)


def main(argv):
    args = process_command_line(argv)
    sips = sipsik.Connection()    
    sips.connect()
    sips.execute_command(args.command)
    sips.disconnect()
    
if __name__ == '__main__':
    main(sys.argv[1:])
