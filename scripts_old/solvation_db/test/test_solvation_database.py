# -*- coding: utf-8 -*-


from sqlalchemy.engine import Engine

from solvation_database import db_connect
from solvation_database import create_tables


def test_db_connect():
    engine = db_connect()
    assert isinstance(engine, Engine)


def test_create_tables():
    engine = db_connect()
    create_tables(engine)
    assert set(engine.table_names()) == set(['Molecule', 'Conformation',
                                     'Calculation', 'RISMCalculation',
                                     'CMCalculation', 'Units'])



