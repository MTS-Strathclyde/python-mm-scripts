#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 14 10:51:37 2013

@author: a92549
"""

INPUT = """%NProc={0}
%Mem={1}MB
%NoSave
%Chk={2}{3}.chk
# """
IRC_KWARDS = """{func} irc=(rcfc,recalc={recalc},MaxPoints={mp},{dir},tight,stepsize={ss}) tzvp/tzvpfit Geom=AllCheckpoint scrf=checkpoint"""

import argparse
import sys
import shutil

def process_command_line(argv):
    """
    Processes arguments and returns namespace of them
    """
    parser = argparse.ArgumentParser(description="""Prepares IRC calculation
                                    for TS files. Input is TS output file.
                                    """)
    #Positional args
    parser.add_argument('input', metavar='mol.log',
                        help="""Output TS file.""", nargs='+')
    #Optional args
    parser.add_argument('-c', '--chk_ext',
                        help="""Checkpoint file extenison. (chk)""",
                        default='.chk')
    parser.add_argument('-mp', '--max_points',
                                        help='Maximum number of points in direction (70)', default='70')
    parser.add_argument('-ss', '--step_size',
                                        help='Size of single step (5) ', default='5')
    parser.add_argument('-k', '--kwards', help="""Gaussian arguments. Should
                        be surounded by double quotes.
                        ({func} irc=(rcfc,recalc={recalc},MaxPoints={mp},{dir},tight,stepsize={ss}) tzvp/tzvpfit Geom=AllCheckpoint scrf=checkpoint).""",
                        default=IRC_KWARDS)
    parser.add_argument('-p', '--proc', help="""Processors number (8)""",
                        default=8, type=int)
    parser.add_argument('-suf', '--suffix', help="Suffix for produced IRC filename",
                        default='')
    parser.add_argument('-mem', '--memory', help="""Memory in MB. (7000)""",
                        default=7000, type=int)
    parser.add_argument('-f', '--func', help='functional (wb97xd)', default='wb97xd')
    parser.add_argument('-r', '--recalc', help="Recalc fc every n steps (10)",
                        default='10')
    return parser.parse_args(argv)


def main(argv):
    args = process_command_line(argv)
    forward_kwards = args.kwards.format(func=args.func, recalc=args.recalc, dir='forward',ss=args.step_size,mp=args.max_points)
    reverse_kwards = args.kwards.format(func=args.func, recalc=args.recalc, dir='reverse',ss=args.step_size,mp=args.max_points)
    for filename in args.input:
        name = filename[:-4].replace('TS', 'IRC') + args.suffix
        forward_input_params = INPUT.format(args.proc, args.memory, name, 'forward')
        reverse_input_params = INPUT.format(args.proc, args.memory, name, 'reverse')    
        forward_input_string = forward_input_params + forward_kwards + '\n\n'
        reverse_input_string = reverse_input_params + reverse_kwards + '\n\n'
        with open(name + 'forward.com', 'wb') as f:
            f.write(forward_input_string)
        with open(name + 'reverse.com', 'wb') as f:
            f.write(reverse_input_string)
        try:
            shutil.copy(filename[:-4] + args.chk_ext, name + 'forward' + args.chk_ext)
            shutil.copy(filename[:-4] + args.chk_ext, name + 'forward' + args.chk_ext + '.bck')
            shutil.copy(filename[:-4] + args.chk_ext, name + 'reverse' +args.chk_ext)
            shutil.copy(filename[:-4] + args.chk_ext, name + 'reverse' +args.chk_ext + '.bck')            
        except IOError:
            print "Couldnt copy .chk file for " + filename
        
if __name__ == '__main__':
    main(sys.argv[1:])
