#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 30 18:08:06 2013

@author: a92549
"""

import argparse
import shutil
import sys


def process_command_line(argv):
    """
    Processes arguments and returns namespace of them
    """
    parser = argparse.ArgumentParser(description="""Resume job.""")
    #Positional args
    parser.add_argument('molecule', metavar='mol.com',
                        help="""Computation to be resumed.""")
    #Optional args
    parser.add_argument('-c', '--chk_ext',
                        help="""Checkpoint file extenison. (chk)""",
                        default='.chk')
    parser.add_argument('-s', '--suffix',
                        help="""Suffix, which will be added to the molecule
                        name end before extension (resume)""",
                        default='resume')
    parser.add_argument('--no_read_geom',
                        help="""Don't add Geom=AllCheckpoint option""",
                        action='store_true')
    parser.add_argument('--no_chk_basis',
                        help="""Don't add ChkBasis option""",
                        action='store_true')
    parser.add_argument('--no_guess',
                        help="""Don't add Guess=Read option""",
                        action='store_true')
    return parser.parse_args(argv)
    

def handle_kwards(kwards, args):
    """Removes gen and pseudo=read kwards, if present and adds specified options"""
    kward_list = kwards.split()
    for i in range(len(kward_list)):
        if kward_list[i].capitalize() == 'Gen' or \
        kward_list[i].capitalize().startswith('Pseudo'):
            kward_list[i] = ''
    if not args.no_read_geom:
        kward_list.append('Geom=AllCheckpoint')
    if not args.no_chk_basis:
        kward_list.append('ChkBasis')
    if not args.no_guess:
        kward_list.append('Guess=Read')
    return ' '.join(kward_list)
            

def process_header(header, args):
    """Modifies chkpoint filename to mach new one and keywords by adding restart
    and Max number of cycles parameters."""
    specs, kwards = header.split('#')
    #chk modificaton
    spec_lines = specs.splitlines()
    for i in range(len(spec_lines)):
        if spec_lines[i].startswith('%chk='):
            spec_lines[i] = spec_lines[i] + '_' + args.suffix
    new_specs = '\n'.join(spec_lines) + '\n'
    new_kwards = handle_kwards(kwards, args)
    return '#'.join([new_specs, new_kwards])
    
    
def main(argv):
    args = process_command_line(argv)
    with open(args.molecule, 'rb') as f:
        com_txt = f.read()
    header = com_txt.split('\n\n')[0]
    new_header = process_header(header, args) + '\n\n'
    no_ext_name = args.molecule[:-4]
    new_name = no_ext_name + '_' + args.suffix + '.com'
    with open(new_name, 'wb') as f:
        f.write(new_header)
    shutil.copy(no_ext_name + args.chk_ext, no_ext_name + '_' + args.suffix + args.chk_ext)
    
    
if __name__ == '__main__':
    main(sys.argv[1:])