# -*- coding: utf-8 -*-
"""
Created on Wed Jan  8 17:31:24 2014

@author: max
"""

import os
import glob
import re
import pybel
from datetime import datetime


KCAL_MOL_IN_HARTREE = 627.509469


class SMVLE_calculation(object):
    """
    Main class for storing and parsing results of smvle calculation.
    
    Correct smvle calculation folder should contain both gas and smvle 
    calculations named *gas and *smvle respectively.
    Inputs should be in gamess input format and have extension inp
    Outputs should have an extension out.
    
    As a result of successful parsing populates instance variable
    file_path_dic with following keys:
    'gas_input', 'solv_input', 'gas_output', 'solv_output'
    
    The same keys are used in file_dic. This dictionary contains content of
    the files themselves.
    """
    job_name = '*gas.inp'  # Extension of one of the files
                                     # which name should be assumed to be
                                     # a job name
    
    named_files = {'{job_name}gas.inp' : 'gas_input',
                   '{job_name}gas.out' : 'gas_output',
                   '{job_name}smvle.inp' : 'solv_input',
                   '{job_name}smvle.out' : 'solv_output'}
          
    input_format = 'inp'
    
    timestring_format = "%a %b %d %H:%M:%S GMT %Y"
    
    # ============= Calculation Defaults ==============================
    
    software = 'GAMESS'
    solvent_details = '{"solvent": "water"}'
    server = 'ferrari.phys.strath.ac.uk'
    solv_method = 'SMVLE'
    
    
    def __init__(self, rism3d_folder):
        if not os.path.isdir(rism3d_folder):
            raise OSError("{} is not a folder.".format(rism3d_folder))
            
        self.file_path_dic = {}
        self.file_dic = {}
        self.solvation_energy = None
        self.extra_output = {}
        self.pymol = None
        self.date = None
        self.runtime = None
        self.temperature = None
        self.theory = None
        self.version = None
        
        self.calc_abs_path = os.path.realpath(rism3d_folder)
        self.rism3d_folder = os.path.relpath(self.calc_abs_path)
        self.parse_folder()
        
    def parse_folder(self):
        """
        Populates self.file_path_dic with file locations.
        """
        self._set_job_name()
        self._find_named_files()
        self._load_files()
        self._parse_input()
        self._check_results()
        try:
            self._parse_results()
        except ValueError:
            print "Couldn't parse results correctly for {}".format(self.job_name)
            raise TypeError
        self._load_molecule()
        self._calculate_runtime()
    
    def _set_job_name(self):
        """Finds name of the calculation."""
        job_path = '{}/{}'.format(self.rism3d_folder, self.job_name)
        job_name_list = glob.glob(job_path)
        if len(job_name_list) < 1:
            self._not_found_error(self.job_name)
        elif len(job_name_list) > 1:
            self._many_files_error(self.job_name)
        else:
            name = os.path.split(job_name_list[0])[1]
            len_of_suffix = len(self.job_name) - 1
            self.job_name =  name[:-len_of_suffix]

    def _find_named_files(self):
        """Finds paths of files with definite names."""
        for name, description in self.named_files.iteritems():
            name = name.format(job_name=self.job_name)
            f_path = '{}/{}'.format(self.rism3d_folder, name)
            if os.path.isfile(f_path):
                self.file_path_dic[description] = f_path
            else:
                self._not_found_error(f_path)
        
    def _load_files(self):
        """Populates file_dic."""
        for description, path in self.file_path_dic.items():
            with open(path, 'rb') as f:
                txt = f.read()
            self.file_dic[description] = txt
            
    def _parse_input(self):
        """Parses solvent generation script."""
        #temperature
        regex = re.compile("TEMP=(\d+\.\d*|\d+)")
        r = regex.search(self.file_dic['gas_input'])
        if r:
            self.temperature = r.groups()[0]
        else:
            self.temperature = 298.15
        #theory
        regex = re.compile('(\$contrl.+\$end|\$basis.+ \$end)')
        temp_theory = regex.findall(self.file_dic['gas_input'])
        contrl = temp_theory[0][:-4][7:].strip()
        basis = temp_theory[1][:-4][6:].strip()
        self.theory = contrl + ' ' + basis

    def _check_results(self):
        """Checks whether calculation terminated successfully"""
        if not 'EXECUTION OF GAMESS TERMINATED NORMALLY' in self.file_dic['gas_output']:
            print self.job_name + " didn't finish"
            raise TypeError('Calculation didn\'t finish')
        if not 'EXECUTION OF GAMESS TERMINATED NORMALLY' in self.file_dic['solv_output']:
            print self.job_name + " didn't finish"
            raise TypeError('Calculation didn\'t finish')
            
            
    def _parse_results(self):
        """Uses modified cclib Gaussian class.
        """
        for line in self.file_dic['gas_output'].splitlines():
            if line.startswith('          *         GAMESS VERSION =  '):
                temp = line.split('=')[1]
                temp = temp.split('*')[0]
                self.version = temp.strip()

            if line.startswith('                       TOTAL ENERGY ='):
                temp = line.split('=')[1]
                gas_energy_h = float(temp.strip())

        for line in self.file_dic['solv_output'].splitlines():
            if line.startswith(' TOTAL SMVLE FREE ENERGY      GSMVLE ='):
                temp = line.split('=')[1]
                smvle_energy_h = float(temp.strip())

        self.solvation_energy = (smvle_energy_h - gas_energy_h)*KCAL_MOL_IN_HARTREE
    
        
    def _load_molecule(self):
        """Creates pymol molecule from input file
        Additional code is needed to supress all warnings openbabel is 
        printing"""
        self.pymol = pybel.readstring(self.input_format, self.file_dic['gas_input'])

        
    def _calculate_runtime(self):
        """Finds runtime and date. 
        Assumes that time strings are printed at the beginning
        and at the end of output file."""
        lines = self.file_dic['solv_output'].splitlines()
        start_time = datetime.strptime(lines[0].strip(), self.timestring_format)
        fin_time = datetime.strptime(lines[-1].strip(), self.timestring_format)
        dif = fin_time - start_time
        self.date = fin_time.strftime('%d %b %Y')
        self.runtime = dif.total_seconds()
        
    def _not_found_error(self, ftype):
        info = "Couldn't find {} file in folder {}.".format(ftype, self.rism3d_folder)
        print info
        raise TypeError(info)

    def _many_files_error(self, ftype):
        info = "Found multiple files of type {} in folder {}.".format(ftype, self.rism3d_folder)
        print info        
        raise TypeError(info)
