# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 16:35:51 2014

@author: max
"""

import sys
import csv
import db_interface
import solv_utility
from sqlalchemy.orm.exc import NoResultFound
from db_interface import create_session


def write_row(row, dbi):
    try:
        mol = dbi.get_molecule(row[1])
    except NoResultFound:
        mol = dbi.create_mol(row[1])
        dbi.add_molecule(mol)
        print 'new mol added ' + row[0]
    with open(row[0], 'rb') as f:
        pdb = f.read()
    mol = solv_utility.mol_converter(pdb, 'pdb', 'mol')
    conformation = dbi.bind_conf(mol)
    conformation.Source = row[3]
    conformation.Properties = '{{"filename" : "{0}"}}'.format(row[0])


def read_csv(csv_file_path, session):
    dbi = db_interface.DBInterface(session)
    with open(csv_file_path, 'rb') as f:
        reader =  csv.reader(f)
        reader.next()
        for row in reader:
            write_row(row, dbi)
    session.commit()
    

def main(csv_path):
    ses = create_session()
    read_csv(csv_path, ses)
    ses.close()
    

if __name__ == '__main__':
    main(sys.argv[1])