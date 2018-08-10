# -*- coding: utf-8 -*-

from collections import OrderedDict
import os
import csv
import glob
import re
import pybel

def mol_converter(txt, informat, outformat=None):
    """Converts molecule using pybel."""
    pymol = pybel.readstring(informat, txt)
    if outformat:
        return pymol.write(outformat)
    else:
        return pymol    


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
                   '{job_name}.prmtop' : 'topology',
                   'results.txt' : 'results'}

    varied_name_files = {'water*.sh' : 'solvent'}      
      
    solv_energy_key = 'dGhyd(GF)'
    
    input_format = 'pdb'
    
    timestring_format = "%a %b %d %H:%M:%S GMT %Y"
    
    # ============= Calculation Defaults ==============================
    
    software = 'ambertools'
    version = '12.7'
    force_field = 'gaff'
    water_model = 'SPC'
    bridge  = 'KH'
    solvent_details = '{"solvent": "water"}'
    server = 'ferrari.phys.strath.ac.uk'
    
    
    def __init__(self, rism3d_folder):
        if not os.path.isdir(rism3d_folder):
            raise OSError("{} is not a folder.".format(rism3d_folder))
            
        self.file_path_dic = {}
        self.file_dic = {}
        self.GF = None
        self.extra_output = {}
        #self.date = None
        #self.runtime = None
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
        #self._check_results()
        self._parse_results()
        self._parse_solvent_gen_file()
        self._load_molecule()
        #self._calculate_runtime()
    
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
        self.GF = self.extra_output[self.solv_energy_key]
        
    def _parse_solvent_gen_file(self):
        """Parses solvent generation script."""
        regex = re.compile("TEMPER=(\d+\.\d*|\d+)")
        r = regex.search(self.file_dic['solvent'])
        self.temperature = r.groups()[0]
                
    def _not_found_error(self, ftype):
        info = "Couldn't find {} file in folder {}.".format(ftype, self.rism3d_folder)
        print info
        raise TypeError(info)

    def _many_files_error(self, ftype):
        info = "Found multiple files of type {} in folder {}.".format(ftype, self.rism3d_folder)
        print info
        raise TypeError(info)



def water_concentration(T):
    """Return water concentration for temperature range 253.15K < T < 383.15K.
    
    Uses correlating equation (eq. 2) for specific volume found in
    the doucment by The International Association for the Properties of 
    Water and Steam from 2011
    (http://www.iapws.org/relguide/LiquidWater.pdf)
    Pressure = 0.1 MPa    
    
    >>> round(water_concentration(273.15), 3)
    55.498
    >>> round(water_concentration(298.15), 3)
    55.343
    """
    p0 = 10.0**5    # Pa
    R = 8.31464     # J/mol/K
    Tr = 10.0
    Ta = 593.0
    Tb = 232.0
    a = [1.93763157E-2,
         6.74458446E+3,
        -2.22521604E+5,
         1.00231247E+8,
        -1.63552118E+9,
         8.32299658E+9]
    b = [5.78545292E-3,
        -1.53195665E-2,
         3.11337859E-2,
        -4.23546241E-2,
         3.38713507E-2,
        -1.19946761E-2]
    n = [None, 4., 5., 7., 8., 9.]
    m = [1., 2., 3., 4., 5., 6.]
    def alpha(T): 
        return Tr/(Ta - T)
    def beta(T):
        return Tr/(T - Tb)
    coef = a[0] + b[0]*beta(T)**m[0]
    for i in range(1, 6):
        coef += a[i]*alpha(T)**n[i] + b[i]*beta(T)**m[i]
    v0 = R*Tr/p0*coef  # m3/mol
    return 1/(v0*1000)    # mol/L
    


def comput_ucorr(gf, pmv, t):
    gf = float(gf)
    pmv = float(pmv)
    t = float(t)
    density = water_concentration(t)*6.0221413E-4
    return gf - 3.2217*density*pmv + 0.5783
    


def write_meta_file(temp, dirname):
    name = dirname.split('/')[-1]
    meta_f = OrderedDict()
    inchi = names[name]
    meta_f['InChI'] = inchi
    rism = RISM_3D_calculation(temp)
    meta_f['Temperature'] = rism.temperature
    meta_f['GF'] = rism.GF
    meta_f['PMV'] = rism.file_dic['results'].splitlines()[2].split()[1]
    ucorr = comput_ucorr(meta_f['GF'], meta_f['PMV'], rism.temperature)
    meta_f['UCorr'] = ucorr
    meta_f['Software'] = rism.software
    meta_f['Version'] = rism.version
    meta_f['ForceField'] = rism.force_field
    meta_f['WaterModel'] = rism.water_model
    meta_f['Bridge'] = rism.bridge
#    meta_f['Date'] = rism.date
#    meta_f['Runtime'] = rism.runtime
    meta_f['UCorrMult'] = '-3.2217,0.5783'
    meta_f['InputFile'] = rism.file_path_dic['input'].split('/')[-1]
    meta_f['OutputFile'] = rism.job_name + '.log'
    meta_f['Topology'] = rism.file_path_dic['topology'].split('/')[-1]
    meta_f['SolventGenFile'] = rism.file_path_dic['solvent'].split('/')[-1]
    meta_f_txt = ''
    for k, v in meta_f.iteritems():
        meta_f_txt += '{}, {}\n'.format(k, v)
    return meta_f_txt
    
    

def load_names():
    data_f = open('/home/max/Documents/PhD/calculations/RISM_calcs/Temp_dep_RISM/meta_files/temperature_no_method_no_source.csv', 'rb')
    names = OrderedDict()
    rdr = csv.reader(data_f)
    for row in rdr:
        if row[0].startswith('In'):
            names[row[1]] = row[0]
    data_f.close()
    return names


def iterate_over_calcs():
    root_dir = '/home/max/Documents/PhD/calculations/RISM_calcs/Temp_dep_RISM/db_files/rism_part/raw_calc_files'
    dirs = [os.path.join(root_dir, d) for d in os.listdir(root_dir) if os.path.isdir(os.path.join(root_dir, d))]
    for dirname in dirs:
        temps = [os.path.join(root_dir, dirname, d) for d in os.listdir(os.path.join(root_dir, dirname)) if os.path.isdir(os.path.join(root_dir, dirname, d))]
        for temp in temps:
            meta_f = write_meta_file(temp, dirname)     
            name = dirname.split('/')[-1]
            t = temp.split('/')[-1]
            with open(os.path.join(temp, name + '_' + t + '.dbf'), 'wb') as f:
                f.write(meta_f)

names = load_names()
iterate_over_calcs()


