#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 26 17:09:23 2012

@author: a92549

ver 0.1

a) Generate confs
b) Prune
c) Optimize with mopac
d) Prune

TODO: remove high energy compounds
TODO: implement gaussian file creation and optimization

"""

import sys
import argparse
import generate_conformers as gc
import prune_conformers as pc
import mopac_calculation as mop_calc
import datetime
import time

def process_command_line(argv):
    """
    Processes arguments and returns namespace of them
    """
    parser = argparse.ArgumentParser(description="""Generate conformers
                                    and optimize them with MOPAC.""")
    #Positional args
    parser.add_argument('file', metavar='molec.sdf',
                        help="""Input files. Should be in sdf format.""")
    #Optional args
    return parser.parse_args(argv)

def main(argv):
    args = process_command_line(argv)
    filename = args.file[:-4]
    print ""
    print datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + " Starting conformers."
    print "Input file: " + args.file
    start_time = time.time()
    print ""
    FF_conformer_filename = filename + "_FF_confs.sdf"
    MOPAC_folder = filename + "_MOPAC"
    FF_conformers = gc.main([args.file, FF_conformer_filename])
    print ""
    pruned_FF_conformers = pc.main([FF_conformers, '-o', FF_conformers[:-4] + "_p.sdf"])
    print ""
    mop_calc.main([pruned_FF_conformers, "-d", MOPAC_folder])
    print ""
    pc.main([MOPAC_folder, "-m"])
    end_time = time.time()
    print ""
    print datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + " Finished."
    print "Total time: " + str(round(end_time - start_time, 3))
    
if __name__ == '__main__':
    main(sys.argv[1:])