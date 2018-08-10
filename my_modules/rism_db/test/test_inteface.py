# -*- coding: utf-8 -*-

from db_interface import create_session
from db_interface import DBInterface
from solvation_database import Molecule, Conformation
from solv_utility import mol_2_can, get_RMSD_value
from sqlalchemy.orm.exc import NoResultFound

import pytest

class TestDBInterface:
    ses = create_session()
      
    def test_check_molecule_non_existing(self):
        SMILES='CCCCC'
        dbi = DBInterface(self.ses)
        with pytest.raises(NoResultFound) as e:
            dbi.get_molecule(SMILES)
                
    def test_add_molecule(self):
        smiles = 'CCCCCCC'
        mol = Molecule(SMILES=smiles, IUPACName='Heptane')
        dbi = DBInterface(self.ses)
        dbi.add_molecule(mol)
        data_mol = dbi.get_molecule(smiles)
        assert str(data_mol.SMILES) == str(mol.SMILES)

    def test_find_most_similair_conformation(self):
        molref = """
 OpenBabel12071321472D

  2  1  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
M  END"""
        mol = """
 OpenBabel12071321472D

  2  1  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
M  END
"""
        can = mol_2_can(mol)
        dbi = DBInterface(self.ses)
        db_ethane = Molecule(SMILES='CC')
        dbi.add_molecule(db_ethane)
        db_mol = dbi.get_molecule(can)
        db_mol.add_conformation(Conformation(Mol=molref))
        self.ses.add(db_mol)
        db_conf = db_mol.find_most_similair_conf(mol)[0]
        assert str(db_conf.Mol) == molref
        
    def test_add_0_rmsd_conf_existing_mol(self):
        coord = """
 OpenBabel12071320143D

  5  4  0  0  0  0  0  0  0  0999 V2000
    1.0663    0.0688   -0.0295 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.1585    0.0688   -0.0295 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.7023    0.7365    0.7544 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.7023   -0.9440    0.1568 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.7023    0.4138   -0.9997 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  3  1  0  0  0  0
  1  4  1  0  0  0  0
  1  5  1  0  0  0  0
M  END
"""
        conf = Conformation(Mol=coord, Molecule_SMILES='C')
        can = mol_2_can(coord)
        dbi = DBInterface(self.ses)        
        db_methane = Molecule(SMILES='C')
        dbi.add_molecule(db_methane)
        db_mol = dbi.get_molecule(can)
        db_mol.add_conformation(conf)
        self.ses.add(db_mol)
        db_conf = db_mol.get_conformation(coord, 0.0000001)
        assert db_conf == conf
        #self.ses.commit()

TestDBInterface.ses.rollback()
TestDBInterface.ses.close()