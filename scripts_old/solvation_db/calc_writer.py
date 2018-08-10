# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 11:32:06 2013

@author: max
"""

from solvation_database import RISMCalculation, CMCalculation
from db_interface import DBInterface

    
class Writer(object):
    def __init__(self, session):
        """Takes session instance as an argument"""
        self.session = session
        self.dbi = DBInterface(self.session)
        
    def write_rism(self, rism_calc_instance, conf_attributes={},
                                             rism_attributes={}):
        """Save rism calculation.
        Source of conformation,
        Temperature of RISM calculation
        are set thorugh attrubutes dictionaries."""
        mol = rism_calc_instance.pymol.write('mol')
        dbconf = self.dbi.bind_conf(mol)
        rism = RISMCalculation(**rism_attributes)
        rism.Date = rism_calc_instance.date
        rism.InputFile = rism_calc_instance.file_dic['input']
        rism.ParametersFile = rism_calc_instance.file_dic['parameters']
        rism.Path = rism_calc_instance.calc_abs_path
        rism.Results = rism_calc_instance.file_dic['results']
        rism.Runtime = rism_calc_instance.runtime
        rism.Server = rism_calc_instance.server
        rism.Software = rism_calc_instance.software
        rism.SolvGenFile = rism_calc_instance.file_dic['solvent']
        rism.SolvE = rism_calc_instance.solvation_energy
        rism.StdOutput = rism_calc_instance.file_dic['output']
        rism.Temperature = rism_calc_instance.temperature
        rism.Topology = rism_calc_instance.file_dic['topology']
        rism.Version = rism_calc_instance.version
        dbconf.Calculations.append(rism)
        self.dbi.update_dbconf(dbconf, conf_attributes)
        self.session.add(dbconf)
        
    def write_smd(self, smd_calc_instance, conf_attributes={}):
        """Save smd calculation."""
        mol = smd_calc_instance.pymol.write('mol')
        dbconf = self.dbi.bind_conf(mol)
        smd = CMCalculation()
        smd.Date = smd_calc_instance.date
        smd.InputFile = smd_calc_instance.file_dic['input']
        smd.Path = smd_calc_instance.calc_abs_path
        smd.Runtime = smd_calc_instance.runtime
        smd.Server = smd_calc_instance.server
        smd.Software = smd_calc_instance.software
        smd.SolvE = smd_calc_instance.solvation_energy
        smd.SolvMethod = smd_calc_instance.solv_method
        smd.Output = smd_calc_instance.file_dic['output']
        smd.Theory = smd_calc_instance.theory
        smd.Temperature = smd_calc_instance.temperature
        smd.Version = smd_calc_instance.version
        dbconf.Calculations.append(smd)
        self.dbi.update_dbconf(dbconf, conf_attributes)
        self.session.add(dbconf)
        
    def write_smvle(self, smvle_calc_instance, conf_attributes={}):
        """Save smvle calculation."""
        mol = smvle_calc_instance.pymol.write('mol')
        dbconf = self.dbi.bind_conf(mol)
        smvle = CMCalculation()
        smvle.Date = smvle_calc_instance.date
        smvle.InputFile = smvle_calc_instance.file_dic['gas_input'] + '\n\n\n' + \
                        smvle_calc_instance.file_dic['solv_input']
        smvle.Path = smvle_calc_instance.calc_abs_path
        smvle.Runtime = smvle_calc_instance.runtime
        smvle.Server = smvle_calc_instance.server
        smvle.Software = smvle_calc_instance.software
        smvle.SolvE = smvle_calc_instance.solvation_energy
        smvle.SolvMethod = smvle_calc_instance.solv_method
        smvle.Output = smvle_calc_instance.file_dic['gas_output'] + '\n\n\n' + \
                     smvle_calc_instance.file_dic['solv_output']
        smvle.Theory = smvle_calc_instance.theory
        smvle.Temperature = smvle_calc_instance.temperature
        smvle.Version = smvle_calc_instance.version
        dbconf.Calculations.append(smvle)
        self.dbi.update_dbconf(dbconf, conf_attributes)
        self.session.add(dbconf)
                
