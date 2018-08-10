# -*- coding: utf-8 -*-

from solvation_database import Molecule, Conformation
from solv_utility import get_RMSD_value
from rism_parser import RISM_3D_calculation
from calc_writer import Writer
from db_interface import create_session, DBInterface

import pytest
import pybel



class TestDBWriter:
    ses = create_session()
    def test_write_new_mol_and_new_0_conf(self):
        rism_toluene = RISM_3D_calculation('test/toluene')
        can = 'Cc1ccccc1'
        wr = Writer(self.ses)
        wr.write_rism(rism_toluene)
        dbi = DBInterface(self.ses)
        db_mol = dbi.get_molecule(can)
        db_can = str(db_mol.SMILES)
        db_conf = db_mol.find_0rmsd_conformation()
        pymol = pybel.readstring('mol', str(db_conf.Mol))
        assert can == db_can
#        assert xyz == db_xyz
        self.ses.rollback()

    def test_write_new_mol_and_new_conf_and_conf_props(self):
        rism_toluene = RISM_3D_calculation('test/toluene')
        can = 'Cc1ccccc1'
        wr = Writer(self.ses)
        wr.write_rism(rism_toluene, {'Source' : 'David'})
        dbi = DBInterface(self.ses)
        db_mol = dbi.get_molecule(can)
        db_conf = db_mol.find_0rmsd_conformation()
        assert str(db_conf.Source) == 'David'
        self.ses.rollback()
                
    def test_update_existing_conf(self):
        rism_toluene = RISM_3D_calculation('test/toluene')
        rism_toluene_copy = RISM_3D_calculation('test/toluene_copy')        
        can = 'Cc1ccccc1'
        wr = Writer(self.ses)
        wr.write_rism(rism_toluene)
        dbi = DBInterface(self.ses)
        db_mol = dbi.get_molecule(can)
        db_conf = db_mol.find_0rmsd_conformation()        
        assert str(db_conf.Source) == 'None'
        wr.write_rism(rism_toluene_copy, {'Source' : 'David'})        
        dbi = DBInterface(self.ses)
        db_mol = dbi.get_molecule(can)
        db_conf = db_mol.find_0rmsd_conformation()
        assert str(db_conf.Source) == 'David'
        self.ses.rollback()
        
    def test_write_rism_calculation(self):
        rism_toluene = RISM_3D_calculation('test/toluene')
        can = 'Cc1ccccc1'
        wr = Writer(self.ses)
        wr.write_rism(rism_toluene, {'Source' : 'David'})
        dbi = DBInterface(self.ses)
        db_mol = dbi.get_molecule(can)
        db_conf = db_mol.find_0rmsd_conformation()
        db_calc = db_conf.get_rism_calculations()[0]
        assert db_calc.SolvE == 1.34305025E+001
        assert db_calc.Temperature == '298'
        self.ses.rollback()
    

TestDBWriter.ses.rollback()
TestDBWriter.ses.close()