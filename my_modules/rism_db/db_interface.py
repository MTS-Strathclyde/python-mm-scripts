# -*- coding: utf-8 -*-

from sqlalchemy.orm import sessionmaker

from solvation_database import db_connect
from solvation_database import create_tables
from solvation_database import Molecule
from solv_utility import get_IUPAC_name, inchi_2_smi


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
                Molecule.InChI == inchi)
        return query.one()

    def create_mol(self, inchi, mol, iupac=None, smiles=None):
        """Will create new molecule row. If iupac isn't provided
        will try to find it in the internet."""
        if not iupac:
            iupac = get_IUPAC_name(inchi)
        if not smiles:
            smiles = inchi_2_smi(inchi)
        return Molecule(InChI=inchi, IUPACName=iupac, SMILES=smiles,
                        Mol=mol)
            
    def add_molecule(self, dbmol):
        """Adds given molecule to the database.
        Molecule should be instance of class Molecule from solvation_database."""
        self.session.add(dbmol)                  
        


def check_connection():
    ses = create_session()
    dbi = DBInterface(ses)
    dbi.get_molecule('CC')
    