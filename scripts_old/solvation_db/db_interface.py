# -*- coding: utf-8 -*-

from sqlalchemy.orm.exc import NoResultFound
from sqlalchemy.orm import sessionmaker

from solvation_database import db_connect
from solvation_database import create_tables
from solvation_database import Molecule, Conformation
from solv_utility import get_IUPAC_name, mol_2_inchi, inchi_2_smi


def create_session():
    """Connects to a database and creates new session instance."""
    engine = db_connect()
    create_tables(engine)
    Session = sessionmaker(bind=engine)
    return Session()


class DBInterface(object):
    def __init__(self, session):
        self.session = session
         
    def get_molecule(self, inchi):
        """Accepts inchi string.
        If molecule doesn't exist returns false."""
        query = self.session.query(Molecule).filter(\
                Molecule.INCHI == inchi)
        return query.one()

    def create_mol(self, inchi, iupac=None, smiles=None):
        """Will create new molecule row. If iupac isn't provided
        will try to find it in the internet."""
        if not iupac:
            iupac = get_IUPAC_name(inchi)
        if not smiles:
            smiles = inchi_2_smi(inchi)
        return Molecule(INCHI=inchi, IUPACName=iupac, RDK_SMILES=smiles)
            
    def add_molecule(self, dbmol):
        """Adds given molecule to the database.
        Molecule should be instance of class Molecule from solvation_database."""
        print 'Adding molecule {}'.format(dbmol.IUPACName)
        print 'It\'s InChI is {}'.format(dbmol.INCHI)
        raw_input('Press key to continue ')
        self.session.add(dbmol)                  
        
    def bind_conf(self, mol):
        """Returns conformation row corresponding to given coordinates.
        Default tolerance is 0.05.
        If corresponding cofnormation doesn't exist will create a new entry
        and return it.
        If corresponding molecule doesn't exist: will create a new entry, try
        to find IUPAC name for it, add conformation to it and return conformation.
        
        In all cases will add created conformation to session"""
        inchi = mol_2_inchi(mol)
        try:
            dbmol = self.get_molecule(inchi)
            try:
                dbconf = dbmol.get_conformation(mol)
            except NoResultFound:
                dbconf = Conformation(Mol=mol)
                dbmol.add_conformation(dbconf)
        except NoResultFound:
            dbmol = self.create_mol(inchi)
            dbconf = Conformation(Mol=mol)
            dbmol.add_conformation(dbconf)
            print 'Adding molecule {}'.format(dbmol.IUPACName)
            print 'It\'s InChI is {}'.format(dbmol.INCHI)
            raw_input('Press key to continue ')
            
        self.session.add(dbconf)
        return dbconf
        
    def update_dbconf(self, dbconf, conf_attributes):
        for conf_attr in conf_attributes:
            setattr(dbconf, conf_attr, conf_attributes[conf_attr])



def check_connection():
    ses = create_session()
    dbi = DBInterface(ses)
    dbi.get_molecule('CC')
    