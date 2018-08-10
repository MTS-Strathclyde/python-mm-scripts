# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 12:09:04 2014

@author: max

Rows in csv should have following format

<INCHI>,<IUPAC>,<SMILES>
<temperature1>, <energy1>, <method1>, <source1>
<temperature2>, <energy2>, <method2>, <source2>
<temperature3>, <energy3>, <method3>, <source3>,
...

"""

import db_interface
import csv
from sqlalchemy.orm.exc import NoResultFound
from solvation_database import Experiment


def read_csv(csv_file_path, session):
    dbi = db_interface.DBInterface(session)
    with open(csv_file_path, 'rb') as f:
        rdr = csv.reader(f)
        mol = None
        for row in rdr:
            if len(row) == 3:
                inchi, iupac, smi = row
                if mol:
                    dbi.session.add(mol)
                try:
                    mol = dbi.get_molecule(inchi)
                except NoResultFound:                    
                    mol = dbi.create_mol(inchi, iupac, smi)                
            else:
                mol.Experiments.append(Experiment(Temperature=row[0],
                                       SolvEnergy=row[1], Method=row[2],
                                       Source=row[3]))
        dbi.session.add(mol)
    session.commit()             





