#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  9 18:14:26 2013

@author: max

Solvation database command line inteface.

Usage:
    slvdb_cli.py add rism <base_folder>
    slvdb_cli.py add exp_data <csv_file>
    slvdb_cli.py create_db
    slvdb_cli.py -h | --help
    slvdb_cli.py --version
    
Options:
    --temp=<K>          Temperature of calculation in K [default: 298].
    --3source=<Name>    Source of 3D structures
    -h --help           Show this screen.
    --version           Show version.
"""

from docopt import docopt
import os


from db_interface import DBInterface
from db_interface import check_connection
from rism_parser import add_meta_f_to_db
from db_interface import create_session
from sqlalchemy.orm.exc import NoResultFound


def get_dbf(filenames):
    return [f for f in filenames if f.endswith('.dbf')]


def parse_rism_folders(base_folder):
    ses = create_session()
    dbi = DBInterface(ses)
    dirs = os.walk(base_folder)
    for p, _, filenames in dirs:
        dbfs = get_dbf(filenames)
        for dbf in dbfs:
            try:
                add_meta_f_to_db(dbf, p, dbi)
            except NoResultFound:
                print 'No results found for {}'.format(dbf)
    ses.commit()
    ses.close()


def main():
    """Should only be used to handle arguments from command line.
    Accepts list of folders and calls parse_rism_folders to handle the rest
    """
    args = docopt(__doc__, version=0.01)
    if args['create_db']:
        check_connection()
    if args['add']:
        if args['rism']:
            parse_rism_folders(args['<base_folder>'])

#        if args['exp_data']:
#        ses = create_session()
#        read_csv(args['<csv_file>'], ses)
#        ses.close()
            

if __name__ == '__main__':
    main()
    
