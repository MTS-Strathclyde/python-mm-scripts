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
    parser = argparse.ArgumentParser(description="""Restart unexpectedly ended
                                    calculation.""")
    #Positional args
    parser.add_argument('molecule', metavar='mol.com',
                        help="""Computation to be restarted.""")
    #Optional args
    parser.add_argument('-c', '--chk_ext',
                        help="""Checkpoint file extenison. (chk)""",
                        default='.chk')
    parser.add_argument('-s', '--suffix',
                        help="""Suffix, which will be added to the molecule
                        name end before extension (restart)""",
                        default='restart')
    return parser.parse_args(argv)
    

def handle_opt(opt):
    """Adds restart and max cycle to different opt specification"""
    opt_list = opt.split('=')
    if len(opt_list) == 1:
        return 'opt=(Restart)'
    else:
        if opt_list[1].endswith(')'):    
            new_specs = '(' + opt_list[1].strip('()') + ',Restart)'
            return 'opt=' + new_specs
        elif not opt_list[1].startswith('('):    #single kward case
            return 'opt=(' + opt_list[1] + ',Restart)'
        else:
            return 'opt=' + opt_list[1] + 'Restart,'
    

def process_header(header, suffix):
    """Modifies chkpoint filename to mach new one and keywords by adding restart
    and Max number of cycles parameters."""
    specs, kwards = header.split('#')
    #chk modificaton
    spec_lines = specs.splitlines()
    for i in range(len(spec_lines)):
        if spec_lines[i].startswith('%chk='):
            spec_lines[i] = spec_lines[i] + '_' + suffix
    new_specs = '\n'.join(spec_lines) + '\n'
    #kwards modification
    separate_kwards = kwards.split()
    for i in range(len(separate_kwards)):
        if separate_kwards[i].startswith('opt'):
            new_opt = handle_opt(separate_kwards[i])
            separate_kwards[i] = new_opt
    new_kwards = ' '.join(separate_kwards)
    return '#'.join([new_specs, new_kwards])
    
    
def main(argv):
    args = process_command_line(argv)
    with open(args.molecule, 'rb') as f:
        com_txt = f.read()
    header = com_txt.split('\n\n')[0]
    new_header = process_header(header, args.suffix) + '\n\n'
    no_ext_name = args.molecule[:-4]
    new_name = no_ext_name + '_' + args.suffix + '.com'
    with open(new_name, 'wb') as f:
        f.write(new_header)
    shutil.copy(no_ext_name + args.chk_ext, no_ext_name + '_' + args.suffix + args.chk_ext)
    
    
if __name__ == '__main__':
    main(sys.argv[1:])