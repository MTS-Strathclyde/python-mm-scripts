#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  9 18:14:26 2013

@author: max

Solvation database command line inteface.

Usage:
    slvdb_cli.py add [--temp=<K> --3source=<Name>] (rism|smd|smvle) <folder> ...
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
import warnings
from sqlalchemy.exc import IntegrityError

from db_interface import check_connection
from rism_parser import RISM_3D_calculation
from smd_parser import SMD_calculation
from smvle_parser import SMVLE_calculation
from calc_writer import Writer
from db_interface import create_session
from exp_writer import read_csv


conf_attributes = {'Source' : ''}
rism_attributes = {'Temperature' : 298}


def parse_args(args):
    """Find important options in arguments and add them to dictionaries."""
    rism_attributes['Temperature'] = args['--temp']
    conf_attributes['Source'] = args['--3source']


def parse_rism_folders(folders):
    """Accepts list (or tuple) of folders and writes
    3DRISM_calculation data to db.
    """
    ses = create_session()
    wr = Writer(ses)    
    for folder in folders:
        try:
            rism = RISM_3D_calculation(folder)
            wr.write_rism(rism, conf_attributes, rism_attributes)
            ses.commit()
            print folder + ' written'
        except OSError, e:
            warnings.warn(str(e), UserWarning)
        except TypeError, e:
            warnings.warn(str(e), UserWarning)
#        except IntegrityError:
#            warnings.warn('Integrity Error', UserWarning)
            
    ses.close()


def parse_smd_folders(folders):
    """Accepts list (or tuple) of folders and writes
    smd_calculation data to db.
    """
    ses = create_session()
    wr = Writer(ses)    
    for folder in folders:
        try:
            smd = SMD_calculation(folder)
            wr.write_smd(smd, conf_attributes)
            ses.commit()
            print folder + ' written'
        except OSError, e:
            warnings.warn(str(e), UserWarning)
        except TypeError, e:
            warnings.warn(str(e), UserWarning)
#        except IntegrityError:
#            warnings.warn('Integrity Error', UserWarning)
            
    ses.close()


def parse_smvle_folders(folders):
    """Accepts list (or tuple) of folders and writes
    smvle_calculation data to db.
    """
    ses = create_session()
    wr = Writer(ses)    
    for folder in folders:
        print folder
        try:
            smvle = SMVLE_calculation(folder)
            wr.write_smvle(smvle, conf_attributes)
            ses.commit()
            print folder + ' written'
        except OSError, e:
            warnings.warn(str(e), UserWarning)
        except TypeError, e:
            warnings.warn(str(e), UserWarning)
#        except IntegrityError:
#            warnings.warn('Integrity Error', UserWarning)
            
    ses.close()


def main():
    """Should only be used to handle arguments from command line.
    Accepts list of folders and calls parse_rism_folders to handle the rest
    """
    args = docopt(__doc__, version=0.01)
#    print args
    parse_args(args)
    if args['create_db']:
        check_connection()
    if args['add']:
        if args['rism']:
            parse_rism_folders(args['<folder>'])
        if args['smd']:
            parse_smd_folders(args['<folder>'])      
        if args['smvle']:
            parse_smvle_folders(args['<folder>'])        
        if args['exp_data']:
            ses = create_session()
            read_csv(args['<csv_file>'], ses)
            ses.close()
            

if __name__ == '__main__':
    main()
    
