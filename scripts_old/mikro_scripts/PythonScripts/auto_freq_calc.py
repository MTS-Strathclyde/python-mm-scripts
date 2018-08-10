#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  3 00:35:46 2013

@author: a92549
"""

import sys


REMOVE_KWARDS_STARTING_WITH = ['opt', 'gen']
ADD_KWARDS = ['freq=noraman', 'genchk', 'chkbas', 'guess=tcheck', 'geom=allcheck']


def is_finished(log_file_txt):
    """ Checks, weather job ended successfully"""
    log_file_lines = log_file_txt.split('\n')
    return log_file_lines[-2][:19] == ' Normal termination'


def change_kwards(kward_txt):
    """Removes everything, which starts with REMOVE_KWARDS_STARTING_WITH
    and adds all from ADD_KWARDS"""
    old_kward_list = kward_txt.split()
    new_kward_list = []
    for kw in old_kward_list:
        for rem_kw in REMOVE_KWARDS_STARTING_WITH:
            if not kw.startswith(rem_kw):
                new_kward_list.append(kw)
    new_kward_list += ADD_KWARDS
    return ' '.join(new_kward_list)


def main(argv):
    with open(argv[0] + '.log', 'rb') as f:
        log_file_txt = f.read()
    if is_finished(log_file_txt):
        with open(argv[0] + '.com', 'rb') as f:
            com_file_txt = f.read()
            header = com_file_txt.split('\n\n')[0]
            specs, kwards = header.split('#')
            new_kwards = change_kwards(kwards)
            new_header = '#'.join([specs, new_kwards])
            freq_calc_txt = new_header + '\n\nFreq calculation.\n\n'        
        with open(argv[0] + '_freq.com', 'wb') as f:
            f.write(freq_calc_txt)


if __name__ == '__main__':
    main(sys.argv[1:]) #passed argument is filename of calculation without extension