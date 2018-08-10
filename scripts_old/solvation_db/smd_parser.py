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


class SMD_calculation(object):
    """
    Main class for storing and parsing results of smd calculation.
    
    As a result of successful parsing populates instance variable
    file_path_dic with following keys:
    'input', 'output'
    
    The same keys are used in file_dic. This dictionary contains content of
    the files themselves.
    """
    job_name = '*.inp'  # Extension of one of the files
                                     # which name should be assumed to be
                                     # a job name
    
    named_files = {'{job_name}.inp' : 'input',
                   '{job_name}.out' : 'output'}
          
    input_format = 'inp'
    
    timestring_format = "%a %b %d %H:%M:%S GMT %Y"
    
    # ============= Calculation Defaults ==============================
    
    software = 'GAMESS'
    solvent_details = '{"solvent": "water"}'
    server = 'ferrari.phys.strath.ac.uk'
    solv_method = 'SMD'
    
    
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
        self.parse_rism_folder()
        
    def parse_rism_folder(self):
        """
        Populates self.file_path_dic with rism calculation file locations.
        In case of raises a warning and returns 1     

        >>> r3d = RISM_3D_calculation('test/toluene')
        >>> len(r3d.file_path_dic)
        6
        >>> r3d.file_path_dic['input']
        'test/toluene/toluene.pdb'
        
        rism folder should contain:
            <input structure name>.pdb
            <input structure name>.prmtop
                topology file
                this file is used to find out what <input structure name> is
            run3drismgaff.sh
                parameters file
            results.txt
            out.*
        """
        self._set_job_name()
        self._find_named_files()
        self._load_files()
        self._parse_input()
        self._check_results()
        self._parse_results()
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
            self.job_name =  '.'.join(name.split('.')[:-1])

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
        r = regex.search(self.file_dic['input'])
        if r:
            self.temperature = r.groups()[0]
        else:
            self.temperature = 298.15
        #theory
        regex = re.compile('(\$contrl.+\$end|\$basis.+ \$end)')
        temp_theory = regex.findall(self.file_dic['input'])
        contrl = temp_theory[0][:-4][7:].strip()
        basis = temp_theory[1][:-4][6:].strip()
        self.theory = contrl + ' ' + basis

    def _check_results(self):
        """Checks whether calculation terminated successfully"""
        if not 'EXECUTION OF GAMESS TERMINATED NORMALLY' in self.file_dic['output']:
            print self.job_name + " didn't finish"
            raise TypeError('Calculation didn\'t finish')
            
    def _parse_results(self):
        """Uses modified cclib Gaussian class.
        """
        for line in self.file_dic['output'].splitlines():
            if line.startswith('          *         GAMESS VERSION =  '):
                temp = line.split('=')[1]
                temp = temp.split('*')[0]
                self.version = temp.strip()

            if line[1:25] == 'FREE ENERGY OF SOLVATION' and line.find('1 ATM') == -1:
                temp = line.split()
                #Take the next number after =
                #In KCAL/MOL
                self.solvation_energy = float(temp[temp.index("=") + 1])            
        
    def _load_molecule(self):
        """Creates pymol molecule from input file
        Additional code is needed to supress all warnings openbabel is 
        printing"""
        self.pymol = pybel.readstring(self.input_format, self.file_dic['input'])

        
    def _calculate_runtime(self):
        """Finds runtime and date. 
        Assumes that time strings are printed at the beginning
        and at the end of output file."""
        lines = self.file_dic['output'].splitlines()
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
