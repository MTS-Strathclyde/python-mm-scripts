# -*- coding: utf-8 -*-
"""
Created on Fri May 23 09:43:31 2014

@author: max
"""


class PrmtopFile(object):
    def __init__(self, prmtop_name):
        self.prmtop_name = prmtop_name
        with open(prmtop_name, 'rb') as f:
            self.prm_lines = f.readlines()
        
    def set_charges(self, charge_l=None):
        """If charge_l = None, sets charges to 0"""
        i = 0
        while i < len(self.prm_lines):
            line = self.prm_lines[i]
            if line.startswith('%FLAG CHARGE'):
                if self.prm_lines[i+1].startswith('%FORMAT(5E16.8)'):
                    j = i + 2
                    if not charge_l:
                        while self.prm_lines[j].startswith(' '):
                            chrgs = self.prm_lines[j].split()
                            new_chrgs = '  0.00000000E+00'*len(chrgs)
                            self.prm_lines[j] = new_chrgs + '\n'
                            j += 1
                        break
                    else:
                        while self.prm_lines[j].startswith(' '):
                            chrgs = self.prm_lines[j].split()
                            new_chrgs = ['{: .8E}'.format(charge_l.pop(0)) for chg in chrgs]
                            self.prm_lines[j] = ' ' + ' '.join(new_chrgs) + '\n'
                            j += 1
                        break
                else:
                    raise ValueError('Charge given in unknown format')
            i += 1
    
    def write_prmtop(self, name=None):
        if not name:
            name = self.prmtop_name
        with open(name, 'wb') as f:
            f.writelines(self.prm_lines)
