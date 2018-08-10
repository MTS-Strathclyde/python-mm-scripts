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
from razi.orm import ChemColumn
from razi.chemtypes import Molecule as MOL

import settings


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
    
    ID = Column(Integer, primary_key=True, nullable=False)    
    InChI = Column(TEXT, nullable=False)
    IUPACName = Column(TEXT)
    SMILES = ChemColumn(MOL)
    XYZ = Column(TEXT, nullable=False)
    GeomSource = Column(TEXT)
    Properties = Column(TEXT)

    RISMCalcs = relationship('Water3DRISM', backref=backref('molecule'))
                                 
    Experiments = relationship('Experiment', backref=backref('molecule'))
                                                            
    __table_args__ = (UniqueConstraint('InChI', name='UQ_InChI'),)
                                                                    
    def __repr__(self):
        return 'Molecule : ID = {}, name = {}'.format(self.ID,
                                                         self.IUPACName)

        
class Experiment(Base):
    __tablename__ = 'Experiment'
    
    ExpID = Column(Integer, primary_key=True, nullable=False)
    Temperature = Column(FLOAT, nullable=False)        
    SolvEnergy = Column(FLOAT, nullable=False)
    Method = Column(TEXT)    
    Source = Column(TEXT)
    
    Molecule_ID = Column(Integer, ForeignKey('Molecule.ID'), nullable=False)
    
    def __repr__(self):
        return 'Experiment : ExpID = {}, Temperature = {}, Energy = {}'.\
                                format(self.ExpID, self.Temperature,
                                       self.SolvEnergy)


class Water3DRISM(Base):
    __tablename__ = 'Water3DRISM'
    
    CalcID = Column(Integer, primary_key=True, nullable=False)
    Temperature = Column(FLOAT, nullable=False)
    UCorr = Column(FLOAT, nullable=False)
    ChargeModel = Column(TEXT, nullable=False)
    Software = Column(TEXT, nullable=False)
    Version = Column(TEXT, nullable=False)    
    ForceField = Column(TEXT, nullable=False)
    WaterModel = Column(TEXT, nullable=False)
    Bridge = Column(TEXT, nullable=False)
    Date = Column(DATE)
    Runtime = Column(FLOAT) # seconds
    Description = Column(TEXT)
    
    Molecule_ID = Column(Integer, ForeignKey('Molecule.ID'), 
                                nullable=False)
    
    Extra = relationship("Water3DRISMExtra", uselist=False, backref="water_rism")
    ThermOut = relationship("ThermodynamicOutput", backref="water_rism")


class ThermodynamicOutput(Base):
    __tablename__ = 'ThermodynamicOutput'
    
    OutID = Column(Integer, primary_key=True, nullable=False)
    Property = Column(TEXT, nullable=False)
    TotalValue = Column(FLOAT)
    OContrib = Column(FLOAT)
    HContrib = Column(FLOAT)

    Water3DRISM_CalcID = Column(Integer, ForeignKey('Water3DRISM.CalcID'), 
                                nullable=False)


class Water3DRISMExtra(Base):
    __tablename__ = 'Water3DRISMExtra'
    
    CalcID = Column(Integer, ForeignKey('Water3DRISM.CalcID'), 
                    primary_key=True, nullable=False)
    UCorrMult = Column(TEXT, nullable=False)
    InputFile = Column(TEXT, nullable=False)
    OutputFile = Column(TEXT, nullable=False)
    Topology = Column(TEXT, nullable=False)
    SolventGenFile = Column(TEXT, nullable=False)


class Units(Base):
    __tablename__ = 'Units'
    
    Quantity = Column(TEXT, primary_key=True)
    Unit = Column(TEXT)

    
