#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  3 03:06:12 2013

@author: a92549
"""


B3LYP = "B3LYP freq=NoRaman"
MPW1K = "mpwpw91 IOp(3/76=0572004280)"


import sys
import argparse


def process_command_line(argv):
    """
    Processes arguments and returns namespace of them
    """
    parser = argparse.ArgumentParser(description="""Converts MPW1K to B3LYP""")
    #Positional args
    parser.add_argument('input', metavar='mol.com', help='Input files', nargs='+')
#    #Optional args
    return parser.parse_args(argv)


def replace_kwards(kwards):
    new_txt = None
#    if B3LYP in kwards:
#        kwards = kwards.replace(B3LYP, MPW1K)
#        new_txt = 'MPW1K'
    if MPW1K in kwards:
        kwards = kwards.replace(MPW1K, B3LYP)
        new_txt = 'B3LYP'
    else:
        raise ValueError("Unknown kwards!")
    return kwards, new_txt
    
    

def main(argv):
    args = process_command_line(argv)
    for f_name in args.input:
        with open(f_name, 'rb') as f:
            com_txt = f.read()
        sections = com_txt.split('\n\n')
        header = sections[0]
        specs, kwards = header.split('#')
        new_kwards, new_txt = replace_kwards(kwards)
        sections[0] = '#'.join([specs, new_kwards])
        new_input_txt = '\n\n'.join(sections)
        new_f_name = f_name[:-4] + new_txt + '.com'
        with open(new_f_name, 'wb') as f:
            f.write(new_input_txt)
        

if __name__ == '__main__':
    main(sys.argv[1:])
