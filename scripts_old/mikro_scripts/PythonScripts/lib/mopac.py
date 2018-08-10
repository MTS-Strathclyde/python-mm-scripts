# -*- coding: utf-8 -*-
"""
Created on Fri Feb 15 16:56:42 2013

@author: a92549
"""

def set_kwards(header, kwards='PM7'):
    lines = header.split('\n')
    lines[0] = kwards
    return '\n'.join(lines)

def split_mop(mop):
    """Returns file header and geometry"""
    lines = mop.split('\n')
    header = '\n'.join(lines[:3]) + '\n'
    geom = '\n'.join(lines[3:])
    return header, geom
    
def remove_optimization_from_line(line):
    """Removes all optimized parameters from line"""
    strings = line.split()
    for i in (2, 4, 6):
        strings[i] = '0'
    return ' '.join(strings)
    
def remove_optimization_and_set_kwards(mop, kwards='PM7'):
    """Removes all optimized paramteters from given mopac file.
    And sets kwards."""
    head, geom = split_mop(mop)
    head = set_kwards(head, kwards)
    lines = geom.rstrip().split('\n')
    for i in range(len(lines)):
        non_opt_line = remove_optimization_from_line(lines[i])
        lines[i] = non_opt_line
    new_geom = '\n'.join(lines)
    return head + new_geom + '\n'