# -*- coding: utf-8 -*-
"""
Created on Tue Nov 19 20:13:03 2013

@author: max
"""

from sqlalchemy import create_engine
from sqlalchemy import Column, Integer, UniqueConstraint
from sqlalchemy.engine.url import URL
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import ForeignKey
from sqlalchemy.orm import relationship, backref 
from sqlalchemy.dialects.postgresql import TEXT, FLOAT, DATE
from sqlalchemy.orm.exc import NoResultFound
from razi.orm import ChemColumn
from razi.chemtypes import Molecule as MOL

import settings
from solv_utility import get_RMSD_value, mol_2_inchi
from operator import itemgetter


Base = declarative_base()


def db_connect():
    """
    Preforms database connection using database settings from settings.py.
    Returns sqlalchemy engine instance.
    """
    return create_engine(URL(**settings.DATABASE))


def create_tables(engine):
    Base.metadata.create_all(engine)
    

class Molecule(Base):
    __tablename__ = 'Molecule'
    
    INCHI = Column(TEXT, primary_key=True, nullable=False)
    IUPACName = Column(TEXT)
    Properties = Column(TEXT)
    RDK_SMILES = ChemColumn(MOL)
    
    Conformations = relationship('Conformation', backref=backref('molecule',
                                 cascade="all"))
                                 
    Experiments = relationship('Experiment', backref=backref('molecule', 
                                                            cascade='all'))

    def find_0rmsd_conformation(self):
        f = lambda conformation : getattr(conformation, 'RMSD') == 0
        try:
            return filter(f, self.Conformations)[0]
        except IndexError:
            raise NoResultFound
            
    def find_most_similair_conf(self, Mol):
        """Returns tuple with conformation and rmsd value showing difference
        between submited structure and structure in db."""
        conformers_and_rms_diff = []
        for conf in self.Conformations:
            diff = get_RMSD_value(Mol, conf.Mol)
            conformers_and_rms_diff.append((conf, diff)) 
        if len(conformers_and_rms_diff) > 0:
            return min(conformers_and_rms_diff, key=itemgetter(1))
        else:
            raise NoResultFound
        
    def get_conformation(self, Mol, tolerance=0.05):
        """Finds conformation similair withing given tolerance."""
        try:
            most_sim_conf, rmsd = self.find_most_similair_conf(Mol)
            if rmsd <= tolerance:
                return most_sim_conf
            else:
                raise NoResultFound
        except TypeError:
            raise NoResultFound
            
    def add_conformation(self, dbconf):
        """Method will check whether given conformation belongs to this molecule,
        and set rmsd for given dbcoonf."""
        inchi = mol_2_inchi(dbconf.Mol)
        try:
            assert inchi == str(self.INCHI)
        except AssertionError:
            print inchi
            print str(self.INCHI)
            raise AssertionError
        if len(self.Conformations) == 0:
            dbconf.RMSD = 0
            self.Conformations.append(dbconf)
        else:
            base_conf = self.find_0rmsd_conformation()
            rms = get_RMSD_value(base_conf.Mol, dbconf.Mol)
            dbconf.RMSD = rms
            self.Conformations.append(dbconf)
        
    def __repr__(self):
        return 'Molecule : INCHI = {}'.format(self.INCHI)

        
class Experiment(Base):
    __tablename__ = 'Experiment'
    
    ExpID = Column(Integer, primary_key=True, nullable=False)
    SolvEnergy = Column(FLOAT, nullable=False)
    Temperature = Column(FLOAT, nullable=False)    
    Method = Column(TEXT)    
    Source = Column(TEXT)
    
    Molecule_INCHI = Column(TEXT, ForeignKey('Molecule.INCHI'), nullable=False)
    
    def __repr__(self):
        return 'Experiment : INCHI = {}, Temperature = {}, Energy = {}'.\
                                format(self.Molecule_INCHI, self.Temperature,
                                       self.SolvEnergy)
    

class Conformation(Base):
    __tablename__ = 'Conformation'

    StrID = Column(Integer, primary_key=True, nullable=False)
    Mol = Column(TEXT, nullable=False)
    Source = Column(TEXT)
    RMSD = Column(FLOAT)
    Properties = Column(TEXT)

    Molecule_INCHI = Column(TEXT, ForeignKey('Molecule.INCHI'), nullable=False)

    Calculations = relationship('Calculation', backref=backref('conformation',
                                cascade="all"))
    
    __table_args__ = (UniqueConstraint('RMSD', 'Molecule_INCHI', 
                                       name='UQ_RMSDForGivenMol'),)
    
    def __repr__(self):
        return 'Conformation : INCHI = {}'.format(self.Molecule_INCHI)

    def get_rism_calculations(self):
        """Returns all rism calculations associated with this conformaton."""
        f = lambda x : getattr(x, 'CalculationType') == 'RISMCalculation'
        return filter(f, self.Calculations)        
        
        
class Calculation(Base):
    __tablename__ = 'Calculation'
    
    CalcID = Column(Integer, primary_key=True, nullable=False)
    Temperature = Column(FLOAT, nullable=False)
    SolventDetails = Column(TEXT)
    Date = Column(DATE)
    Software = Column(TEXT, nullable=False)
    Version = Column(TEXT, nullable=False)
    Description = Column(TEXT)
    Server = Column(TEXT, nullable=False)
    Path = Column(TEXT, nullable=False)
    Runtime = Column(FLOAT) # seconds
    ParsedInput = Column(TEXT)
    ParsedOutput = Column(TEXT)
    
    Conformation_StrID = Column(Integer, ForeignKey('Conformation.StrID'), 
                                nullable=False)
    
    CalculationType = Column(TEXT)
    
    __table_args__ = (UniqueConstraint('Server', 'Path', 
                                       name='UQ_Path'),)
        
    __mapper_args__ = {
        'polymorphic_identity': 'Calculation',
        'polymorphic_on': CalculationType
    }



class CMCalculation(Calculation):
    __tablename__ = 'CMCalculation'
    
    CalcID = Column(Integer, ForeignKey('Calculation.CalcID'), 
                    primary_key=True, nullable=False)
    SolvE = Column(FLOAT, nullable=False)
    InputFile = Column(TEXT, nullable=False)
    SolvMethod = Column(TEXT, nullable=False)
    Theory = Column(TEXT, nullable=False)
    Output = Column(TEXT, nullable=False)

    __mapper_args__ = {
        'polymorphic_identity': 'CMCalculation',
    }


class RISMCalculation(Calculation):
    __tablename__ = 'RISMCalculation'
    
    CalcID = Column(Integer, ForeignKey('Calculation.CalcID'), 
                    primary_key=True, nullable=False)
    SolvE = Column(FLOAT, nullable=False)
    InputFile = Column(TEXT, nullable=False)
    ParametersFile = Column(TEXT, nullable=False)
    Topology = Column(TEXT, nullable=False)
    SolvGenFile = Column(TEXT, nullable=False)
    Results = Column(TEXT, nullable=False)
    StdOutput = Column(TEXT)

    __mapper_args__ = {
        'polymorphic_identity': 'RISMCalculation',
    }


class Units(Base):
    __tablename__ = 'Units'
    
    Quantity = Column(TEXT, primary_key=True)
    Unit = Column(TEXT)

    
