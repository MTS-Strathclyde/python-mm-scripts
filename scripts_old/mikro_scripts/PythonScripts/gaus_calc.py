#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 12 13:25:44 2013

@author: a92549
"""


import sipsik
import os
import sys
import argparse

def process_command_line(argv):
    """
    Processes arguments and returns namespace of them
    """
    parser = argparse.ArgumentParser(description="""Runs gaus script on given
                                    molecules in the same directory (but
                                    relative to home2) on sipsik.""")
    #Positional args
    parser.add_argument('molecules', metavar='molec.com', nargs='+',
                        help="""Input molecules.""")
    #Optional args
    parser.add_argument('-s', '--script', help="""Name of used script
                        (rgaus_and_log_sipsik.py)""",
                        default="rgaus_and_log_sipsik.py")
    return parser.parse_args(argv)


def main(argv):
    args = process_command_line(argv)
    sips = sipsik.Connection()
    sips.connect()
    sips_cwd = sips.change_path(os.getcwd())
    for mol in argv:
        sips.execute_command("cd " + sips_cwd + ";" + args.script + " " + mol)
    sips.disconnect()


if __name__ == '__main__':
    main(sys.argv[1:])
