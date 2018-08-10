#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 17:15:47 2016

@author: max
"""

print('Ambertools 16 is required!')
print('First argument should be micelle topology')

import sys
from parmed import gromacs, amber

DEFAULTS = """[ defaults ]
; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ
1               2               yes             0.5         0.5"""

def create_water_free_gmx_top(micelle_top):
    name = micelle_top[:-4]
    with open(micelle_top) as f:
        txt = f.read()
    sections = txt.split('\n\n')
    # get rid of water
    sections[1] = DEFAULTS
    sections.pop(3)
    with_water = sections[-1]
    without_water = []
    for l in with_water.splitlines():
        if not 'SOL' in l:
            without_water.append(l)
    without_water = '\n'.join(without_water)
    sections[-1] = without_water
    non_water_name = name + '_mic.top'
    with open(non_water_name, 'w') as f:
        f.write('\n\n'.join(sections))
    return non_water_name
    

def main(micelle_top):
    no_water_top = create_water_free_gmx_top(micelle_top)
    amb_top_name = no_water_top[:-4] + '.prmtop'
    gmx_top = gromacs.GromacsTopologyFile(no_water_top)
    # ugly fix to get rid of rb_torsions, which are curently unsuported in
    # amber topology (but should be supproted soon)
    gmx_top.rb_torsions = gmx_top.torsion_torsions
    amb_prm = amber.AmberParm.from_structure(gmx_top)
    amb_prm.write_parm(amb_top_name)
    
    
if __name__ == '__main__':
    main(sys.argv[1])
    
    