# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 10:50:29 2013

@author: max
"""

import sys
import csv
import os

from rism_parser import RISM_3D_calculation


class JustGiveMeEnergies(RISM_3D_calculation):
    def __init__(self, rism_folder):
        super(JustGiveMeEnergies, self).__init__(rism_folder)
        self.results = []
        self.parse_energies()
        
    def parse_energies(self):
        with open(self.file_path_dic['results_path'], 'rb') as f:
            results_lines = f.readlines()
        for line in results_lines:
            value = line.split()[1]
            self.results.append(value)
        

def dirty_parse_rism_folder(folders):
    """Accepts list (or tuple) of folders and returns a list of 
    3DRISM_calculation objects.
    """
    if not isinstance(folders, (list, tuple)):
        raise TypeError("""Wrong input type.
                        Should be a list of folders.""")    
    f = open('out.csv', 'ab')
    writer = csv.writer(f)
    for folder in folders:
        if os.path.isdir(folder):
            calculation = JustGiveMeEnergies(folder)
            job_name = os.path.split(calculation.job_name)[1]
            writer.writerow([job_name] + calculation.results)
    f.close()
        
    


def main(argv):
    """Should only be used to handle arguments from command line.
    Accepts list of folders and calls parse_rism_folders to handle the rest
    """
    dirty_parse_rism_folder(argv)
    

if __name__ == '__main__':
    main(sys.argv[1:])
    