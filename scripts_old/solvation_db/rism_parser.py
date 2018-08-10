# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 14:57:23 2013

@author: max

This script 

"""

import glob
import os
import sys
from datetime import datetime
import solv_utility
import re


def listdir_fullpath(d):
    """Lists fullpaths of all files in directory"""
    return [os.path.join(d, f) for f in os.listdir(d)]


class RISM_3D_calculation(object):
    """
    Main class for storing and parsing results of 3drism calculation.
    
    As a result of successful parsing populates instance variable
    file_path_dic with following keys:
    'input', 'output', 'parameters', 'topology',
    'results', 'solvent'
    
    The same keys are used in file_dic. This dictionary contains content of
    the files themselves.
    
    Throws OSError in case submited argument is not a folder.
    Throws TypeError in case folder in not a RISM folder.
    Throws ValueError in case class can't recognize 
    
    """
    job_name = '*.prmtop'  # Extension of one of the files
                                     # which name should be assumed to be
                                     # a job name
    
    named_files = {'{job_name}.pdb' : 'input',
                   'run3drismgaff.sh' : 'parameters',
                   '{job_name}.prmtop' : 'topology',
                   'solvent_gen_file.sh' : 'solvent',
                   'results.txt' : 'results'}

    varied_name_files = {'out.*' : 'output'}      
      
    solv_energy_key = 'dGhyd(GF)'
    
    input_format = 'pdb'
    
    timestring_format = "%a %b %d %H:%M:%S GMT %Y"
    
    # ============= Calculation Defaults ==============================
    
    software = 'ambertools'
    version = '12.7'
    solvent_details = '{"solvent": "water"}'
    server = 'ferrari.phys.strath.ac.uk'
    
    
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
        self._find_varied_name_files()
        self._load_files()
        self._check_results()
        self._parse_results()
        self._parse_solvent_gen_file()
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
    
    def _find_varied_name_files(self):
        """Find paths of files with unspecified names."""
        for name, description in self.varied_name_files.iteritems():
            f_path = '{}/{}'.format(self.rism3d_folder, name)
            candidates_list = glob.glob(f_path)
            if len(candidates_list) < 1:
                self._not_found_error(f_path)
            elif len(candidates_list) > 1:
                self._many_files_error(f_path)
            else:
                self.file_path_dic[description] = candidates_list[0]
    
    def _load_files(self):
        """Populates file_dic."""
        for description, path in self.file_path_dic.items():
            with open(path, 'rb') as f:
                txt = f.read()
            self.file_dic[description] = txt
            
    def _check_results(self):
        """Checks whether calculation terminated successfully"""
        if not '3D-RISM processing complete.' in self.file_dic['output']:
            print self.job_name + " didn't finish"
            raise TypeError('Calculation didn\'t finish')
            
    def _parse_results(self):
        """Assumes that results has following format:
        key= value unit
        key= value unit
        ...
        Value should be convertable to float.
        """
        for line in self.file_dic['results'].splitlines():
            key, value, _ = line.split()
            key = key[:-1]
            self.extra_output[key] = float(value)
        self.solvation_energy = self.extra_output[self.solv_energy_key]
        
    def _parse_solvent_gen_file(self):
        """Parses solvent generation script."""
        regex = re.compile("TEMPER=(\d+\.\d*|\d+)")
        r = regex.search(self.file_dic['solvent'])
        self.temperature = r.groups()[0]
        
    def _load_molecule(self):
        """Load silently to suppress all pybel warnings"""
        self.pymol = solv_utility.silent_mol_converter(self.file_dic['input'], 
                                                       self.input_format)
        
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


if __name__ == '__main__':
    RISM_3D_calculation(sys.argv[1])

