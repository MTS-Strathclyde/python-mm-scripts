# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 14:18:08 2013

@author: max
"""

from rism_parser import listdir_fullpath
from rism_parser import RISM_3D_calculation


import os
import pytest


TEST_DIR = os.path.realpath('./test')


def test_listdir_fullpath():
    test_script_path = os.path.join(TEST_DIR, 'test_rism_parser.py')
    assert test_script_path in listdir_fullpath(TEST_DIR)


class TestRISM_3D_calculation:
    def test_initialization(self):
        with pytest.raises(OSError) as excinfo:
            RISM_3D_calculation('not_a_folder')
        assert excinfo.value.message == 'not_a_folder is not a folder.'

    def test_parse_rism_folder_empty_folder(self):
        folder = 'test/stupid_empty_folder'
        with pytest.raises(TypeError) as excinfo:
            RISM_3D_calculation(folder)
        info = """Couldn't find *.prmtop file in folder {}.""".format(folder)
        assert excinfo.value.message == info

    def test_parse_rism_folder_unambigious(self):
        folder = 'test/scary'
        with pytest.raises(TypeError) as excinfo:
            RISM_3D_calculation(folder)
        info = """Found multiple files of type *.prmtop in folder {}.""".format(folder)
        assert excinfo.value.message == info

    def test_parse_rism_folder_rism_paths(self):
        folder = 'test/toluene'
        calculation = RISM_3D_calculation(folder)
        assert calculation.rism3d_folder == "test/toluene"
        assert calculation.file_path_dic['topology'] == 'test/toluene/toluene.prmtop'
        assert calculation.file_path_dic['input'] == \
            'test/toluene/toluene.pdb'
        assert calculation.file_path_dic['parameters'] == \
            'test/toluene/run3drismgaff.sh'
        assert calculation.file_path_dic['results'] == \
            'test/toluene/results.txt'
        assert calculation.file_path_dic['output'] == \
            'test/toluene/out.425966'
        assert calculation.file_path_dic['solvent'] == \
            'test/toluene/solvent_gen_file.sh'
        assert calculation.runtime == 2790 
            
    def test_parse_rism_folder_input(self):
        folder = 'test/toluene'
        calculation = RISM_3D_calculation(folder)
        assert calculation.file_dic['input'] == \
"""ATOM      1  C1  MOL     1       3.537   1.423   0.000  1.00  0.00
ATOM      2  H1  MOL     1       3.935   2.311  -0.499  1.00  0.00
ATOM      3  H2  MOL     1       3.942   0.539  -0.501  1.00  0.00
ATOM      4  H3  MOL     1       3.909   1.424   1.030  1.00  0.00
ATOM      5  C2  MOL     1       2.028   1.417  -0.027  1.00  0.00
ATOM      6  C3  MOL     1       1.306   2.620  -0.021  1.00  0.00
ATOM      7  H4  MOL     1       1.844   3.565  -0.031  1.00  0.00
ATOM      8  C4  MOL     1      -0.093   2.617  -0.010  1.00  0.00
ATOM      9  H5  MOL     1      -0.635   3.561  -0.009  1.00  0.00
ATOM     10  C5  MOL     1      -0.793   1.407  -0.002  1.00  0.00
ATOM     11  H6  MOL     1      -1.879   1.403   0.003  1.00  0.00
ATOM     12  C6  MOL     1      -0.085   0.203  -0.010  1.00  0.00
ATOM     13  H7  MOL     1      -0.619  -0.744  -0.009  1.00  0.00
ATOM     14  C7  MOL     1       1.314   0.210  -0.021  1.00  0.00
ATOM     15  H8  MOL     1       1.860  -0.731  -0.031  1.00  0.00
TER
END
"""

    def test_parse_rism_folder_files(self):
        folder = 'test/toluene'
        calculation = RISM_3D_calculation(folder)
        results = """dGhyd(KH)= 2.09294808E+001 kcal/mol
dGhyd(GF)= 1.34305025E+001 kcal/mol
PMV= 1.51209629E+002 AA^3
dGhyd(UC)= -2.21336115625 kcal/mol
"""
        assert calculation.file_dic['results'] == results
        
    def test_parse_rism_folder_energy(self):
        folder = 'test/toluene'
        calculation = RISM_3D_calculation(folder)
        energy = 1.34305025E+001
        assert calculation.solvation_energy == energy
    
    def test_parse_rism_folder_mol_conversion(self):
        folder = 'test/toluene'
        calculation = RISM_3D_calculation(folder)
        smi = 'Cc1ccccc1'
        assert calculation.pymol.write('can').strip() == smi




    #assert parse_rism_folder(['toluene']) ==


