#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 12 16:38:15 2013

@author: a92549


Creates gas phase com files.
"""

import sys

def main(argv):
    solvent_model = 'scrf=(smd,solvent=benzene)'
    gas_phase = ' '
    for com in argv:
        with open(com, 'rb') as f:
            txt = f.read()
        if solvent_model in txt:
            parts = txt.split(solvent_model,1)
            new_txt = parts[0] + gas_phase + parts[1]
            if '.chk' in new_txt:
                new_txt.replace('.chk', '_gas.chk')
            elif '\n#' in new_txt:
                new_txt.replace('\n#', '_gas\n#')
            else:
                print "Couldnt fix .chk file in " + com
                continue
            with open(com[:-4] + '_gas.com', 'wb') as f:
                f.write(new_txt)            
            

    
    
if __name__ == '__main__':
    main(sys.argv[1:])