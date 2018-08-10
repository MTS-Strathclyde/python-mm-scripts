#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 18:47:28 2013

@author: max
"""

import csv
import sys
from cclib.parser.gamessparser import GAMESS
from cclib.parser.data import ccData


SMD_attribute_dic = {'smd_internal_energy' : float,
                     'smd_electr_energy' : float,
                     'smd_cds_energy' : float,
                     'smd_slv_free_energy' : float,
                     'smd_std_free_energy' : float}


class SMD_ccDATA(ccData):
    def __init__(self, attributes=None):
        super(SMD_ccDATA, self).__init__(attributes)
        self._attrlist.extend(SMD_attribute_dic.keys())
        self._attrtypes.update(SMD_attribute_dic)


class SMD_Gamess(GAMESS):
    def __init__(self, *args, **kwargs):
        # Call the __init__ method of the superclass
        super(SMD_Gamess, self).__init__(*args, **kwargs)
    
    def extract(self, inputfile, line):
        super(SMD_Gamess, self).extract(inputfile, line)
        
        if line[1:22] == 'DELTA INTERNAL ENERGY' and line.find('A.U.') == -1:
            temp = line.split()
            #Take the next number after =
            #In KCAL/MOL
            self.smd_internal_energy = float(temp[temp.index("=") + 1])

        if line[1:26] == 'ELECTROSTATIC INTERACTION' and line.find('A.U.') == -1:
            temp = line.split()
            #Take the next number after =
            #In KCAL/MOL
            self.smd_electr_energy = float(temp[temp.index("=") + 1])
            
        if line[1:16] == 'CDS INTERACTION' and line.find('A.U.') == -1:
            temp = line.split()
            #Take the next number after =
            #In KCAL/MOL
            self.smd_cds_energy = float(temp[temp.index("=") + 1])
            
        if line[1:25] == 'FREE ENERGY OF SOLVATION' and line.find('1 ATM') == -1:
            temp = line.split()
            #Take the next number after =
            #In KCAL/MOL
            self.smd_slv_free_energy = float(temp[temp.index("=") + 1])            

        if line[1:25] == 'FREE ENERGY OF SOLVATION' and line.find('1 ATM') > 0:
            temp = line.split()
            #Take the next number after =
            #In KCAL/MOL
            self.smd_std_free_energy = float(temp[temp.index("=") + 1])
            

def main(argv):
    with open('results2.csv', 'wb') as f:
        writer = csv.writer(f)
        lines = []
        for f in sys.argv[1:]:
            if f[-4:] == '.out':
                line = []
                line.append(f[:-4])
                parser = SMD_Gamess(f, datatype=SMD_ccDATA)
                data = parser.parse()
                if hasattr(data, 'smd_internal_energy'):
                    line.append(data.smd_internal_energy)
                    line.append(data.smd_electr_energy)
                    line.append(data.smd_cds_energy)
                    line.append(data.smd_slv_free_energy)
                    line.append(data.smd_std_free_energy)
                lines.append(line)
        writer.writerows(lines)


if __name__ == "__main__":
    main(sys.argv[1:])
