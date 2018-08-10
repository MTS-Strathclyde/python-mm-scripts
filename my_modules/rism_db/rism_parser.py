# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 14:57:23 2013

@author: max

This script 

"""

import os
import solvation_database as sd


def add_meta_f_to_db(meta_f, p, dbi):
    """Adds a meta file to database."""
    rism_attributes = sd.Water3DRISM.__dict__.keys()
    extra_attributes = sd.Water3DRISMExtra.__dict__.keys()
    with open(os.path.join(p, meta_f), 'rb') as f:
        txt = f.readlines()
    inchi_line = txt[0]
    if inchi_line.startswith('InChI'):
        print inchi_line
        _, inchi = inchi_line.split(', ')
        inchi = inchi.strip()
        dbmol = dbi.get_molecule(inchi)
        rism = sd.Water3DRISM()
        rism_extra = sd.Water3DRISMExtra()
    else:
        raise ValueError('dbf file must start with InChI, <inchi code>')
    for line in txt[1:]:
        if ',' in line:
            line_l = line.split(', ')
            name = line_l[0].strip()
            values = map(lambda x: x.strip(), line_l[1:])
            if len(line_l) == 2:
                if name in rism_attributes:
                    rism.__setattr__(name, values[0])
                elif name in extra_attributes:
                    if name == 'UCorrMult':
                        rism_extra.__setattr__(name, values[0])
                    else:
                        with open(os.path.join(p, values[0]), 'rb') as f:
                            value = f.read()
                        rism_extra.__setattr__(name, value)
            elif len(line_l) == 4:
                rism_therm = sd.ThermodynamicOutput(Property=name)
                if values[0] != '-':
                    rism_therm.TotalValue = values[0]
                if values[1] != '-':
                    rism_therm.OContrib = values[1]
                if values[2] != '-':
                    rism_therm.HContrib = values[2]
                rism.ThermOut.append(rism_therm)
            else:
                print 'Unknown attribute: {}'.format(name)
    rism.Extra = rism_extra
    dbmol.RISMCalcs.append(rism)
    dbi.add_molecule(dbmol)
    print 'Added molecule {}'.format(dbmol)
                

