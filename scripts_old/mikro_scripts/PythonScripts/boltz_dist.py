#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 18 21:29:32 2013

@author: a92549
"""

import gaussian
import math
import sys

KCAL_GAS_CONST = 1.9858/1000


def boltzman(kcal_list, temp=273.15):
    min_e = min(kcal_list)
    boltz_list = [math.exp(-(energy-min_e)/temp/KCAL_GAS_CONST) for energy in kcal_list]
    Z = sum(boltz_list)
    print [i/Z*100 for i in boltz_list]


def main(args):
    en_dic = gaussian.create_energy_dic(args)
    name_list = []
    kcal_lsit = []
    for mol, vals in en_dic.iteritems():
        name_list.append(mol)
        kcal_lsit.append(vals[2])
    print name_list
    boltzman(kcal_lsit)
    

if __name__ == '__main__':
    main(sys.argv[1:])