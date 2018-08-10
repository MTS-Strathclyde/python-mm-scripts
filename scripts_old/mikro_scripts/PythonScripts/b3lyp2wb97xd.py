#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Sat Apr  6 12:38:17 2013

@author: a92549
"""

import sys

method = 'wb97xd'
basis = 'tzvp/tzvpfit'

def main(argv):
    for qst3 in argv:
        with open(qst3, 'rb') as f:
            txt = f.read()
        parts = txt.split('\n\n')
        header = parts[0]
        header_parts = header.split('#')
        kwards = header_parts[1]
        new_kwards = kwards.replace('b3lyp', method)
        new_kwards = new_kwards.replace('gen', basis)
        header = header_parts[0]
        new_header = header.replace('.chk', '_wb97xd.chk')
        parts[0] = new_header + '#' + new_kwards
        new_txt = '\n\n'.join(parts[:-1]) + '\n\n'
        with open(qst3[:-4] + '_wb97xd.com', 'wb') as f:
            f.write(new_txt)
            
if __name__ == '__main__':
    main(sys.argv[1:])