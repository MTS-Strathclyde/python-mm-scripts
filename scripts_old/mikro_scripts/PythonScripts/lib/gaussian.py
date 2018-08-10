# -*- coding: utf-8 -*-
"""
Created on Sun Dec  2 12:14:35 2012

@author: mishin1991

Module to prepare and send files for gaussian calculation.
"""

__version__ = 1.3

KCAL_IN_HARTREE = 627.509469

import json
import re
import os
import fileutil


class Input(object):
    """Class to prepare gaussian input file.

        Default value for nproc is 4

        Default value for memory is 256 mb

        Chk file does not have default. If pymol is provided, chk = pymol.title

        kwards should be suplied without #
        kwards can be suplied with constructor or in separate method.

        Default value for description is chk name

        Multiplicity, charge and geometry do not have defaults. They are extracted from
        pymol or charge-mult-geometry string, if either is provided.

        Custom basis set should be set manually

        Methods provide_first_section and provide_existing_input_file can parse
        first section or existing gaussian input file and extract various info
        from them.

    """
    def __init__(self):
        self.proc = '%NProc=' + str(8)
        self.memory = '%Mem=' + str(4000) + 'MB'
        self.chk = None
        self.kwards = None
        self.descript = ''
        self.charge = ''
        self.mult = ''
        self.geom = ''
        self.add_info = None
        self.custom_basis = None

    def set_proc(self, proc):
        """Sets number of processors, default is 4"""
        self.proc = '%NProc=' + str(proc)

    def set_memory(self, memory):
        """Sets size of memory used for calculation in MB.
        Default is 256"""
        self.memory = '%Mem=' + str(memory) + 'MB'

    def set_chk(self, chk):
        """Set calculation chkpoint filename (no extension needed).
        It will be set as description as well, if no other description is
        provided."""
        self.chk = '%chk=' + chk
        if not self.descript:
            self.set_description(chk)

    def set_kwards(self, kwards):
        self.kwards = '# ' + kwards

    def set_description(self, descript):
        self.descript = descript

    def set_charge(self, charge):
        self.charge = str(charge)

    def set_mult(self, mult):
        self.mult = str(mult)

    def set_geometry(self, geom):
        """Set molecules geometry. Provided geometry must be a string."""
        self.geom = geom
        
    def get_element_list(self):
        element_list = []
        if isinstance(self.geom, list):
            geom = self.geom[0]
        else:
            geom = self.geom
        for line in geom.split('\n'):
            line_list = line.split()
            element_list.append(line_list[0])
        element_list = list(set(element_list))
        return element_list

    def add_additional_info(self, txt, overwrite=False):
        """This string will go after geometry and before custom basis set
        specifications. Here it is possible to add second geometry in case of TS
        search, additional info about connectivity, solvent models and simmilair stuff.

        Right now add_info can contain many different sections, so by default
        this method will add new additional info after blank line, if existing one
        is present. If overwrite it will simply set new info, destroying previous.

        """
        if self.add_info:
            self.add_info += '\n\n' + txt
        if overwrite or not self.add_info:
            self.add_info = txt

    def provide_first_section(self, txt):
        """Method will parse first section of input file (before first blank line),
        and try to extract nproc, memory, chk and kwards from there."""
        procs = re.findall('NProc=(\d+)', txt, re.IGNORECASE)
        if len(procs) > 0:
            self.set_proc(procs[0])
        mem = re.findall('%Mem=(\d+)', txt, re.IGNORECASE)
        if len(mem) > 0:
            self.set_memory(mem[0])
        chk = re.findall('%chk=(.+)', txt, re.IGNORECASE)
        if len(chk) > 0:
            self.set_chk(chk[0])
        kwrd_lst = txt.split('#')
#        for i in kwrd_lst:
#            print i
        self.set_kwards(kwrd_lst[1])

    def provide_chrg_mult_geom(self, geom_with_charge_and_mult):
        """Provide geometry with charge and multiplicity attached on top of it.
        It can be easily extracted, if gaussian input is split using \n\n and
        3 element of list is selected."""
        chrg_mult_geom_list = geom_with_charge_and_mult.split('\n')
        chrg_mult_list = chrg_mult_geom_list[0].split()
        if len(chrg_mult_list) == 2:
            self.set_charge(chrg_mult_list[0])
            self.set_mult(chrg_mult_list[1])
            geom_list = chrg_mult_geom_list[1:]
            self.set_geometry('\n'.join(geom_list))
        else: #molecule doesn't have geometry specified
            pass

    def provide_existing_input_file(self, old_input, additional_info=False):
        """Provide existing gaussian input file.
        Input is string!
        Method will try to extract kwards, description, charge, multiplicity and
        geometry. if there is additional
        information after geometry section it will
        be extracted as well.

        """
        gaus_lst = old_input.split('\n\n')
        self.provide_first_section(gaus_lst[0])
        self.set_description(gaus_lst[1])
        self.provide_chrg_mult_geom(gaus_lst[2])
        custom_basis_section_number = 3
        if additional_info:
            custom_basis_section_number += 1
            self.add_additional_info(gaus_lst[3])
        if 'GEN' in self.kwards.upper() and 'GENC' not in self.kwards.upper():
            self.custom_basis = gaus_lst[custom_basis_section_number]
            if 'GENECP' in self.kwards.upper() or 'PSEUDO=READ' in self.kwards.upper():
                self.custom_basis = self.custom_basis + '\n\n' + gaus_lst[custom_basis_section_number + 1]
            

    def provide_pymol(self, pymol):
        """Class will use pymol to extract chkpoint name, multiplicity, charge and geometry."""
        self.set_chk(fileutil.get_filename(pymol.title))
        gaus_str = pymol.write('gau')
        gaus_lst = gaus_str.split('\n\n')
        self.provide_chrg_mult_geom(gaus_lst[2])

    def add_custom_basis_from_dic(self, json_basis_dic_path):
        """Takes as input path to  custom basis set dictionary, saved in
        json format. Reads molecules geometry (self.geom) and adds parameters
        for present elements.
        
        Loads ECP as well.
        There has to be json_basis_dictionary with ECP-s
        in the same directory, starting with the same name but having _pseudo
        suffix (Def2-SVP_pseudo.json)"""
        with open(json_basis_dic_path) as dic:
            basis_dic = json.load(dic)
        ECP_path = fileutil.add_to_name_but_keep_ext(json_basis_dic_path, '_pseudo')
        with open(ECP_path) as dic:
            ECP_basis_dic = json.load(dic)
        element_list = self.get_element_list()
        basis_list = [basis_dic[element] for element in element_list]
        ecp_basis_list = [ECP_basis_dic[element.upper()] for element in element_list if element.upper() in ECP_basis_dic]
        basis_string = '****\n'.join(basis_list) + '****\n\n'
        ecp_basis_list = ''.join(ecp_basis_list) + '\n'
        self.custom_basis = basis_string + ecp_basis_list

    def write(self, f=None):
        """Prepares gausian input file. If no argument is supplied,
        returns it as string. Otherwise writes it into given file
        """
        #Create input string
        input_lst = []
        input_lst.append(self.memory)
        input_lst.append(self.proc)
        input_lst.append(self.chk)
        input_lst.append(self.kwards)
        input_lst.append('')
        input_lst.append(self.descript)
        input_lst.append('')
        input_lst.append(self.charge + ' ' + self.mult)
        input_lst.append(self.geom)
        input_lst.append('')
        if self.add_info:
            input_lst.append(self.add_info)
            input_lst.append('')
        if self.custom_basis:
            input_lst.append(self.custom_basis)
            input_lst.append('')
        gaus_input = '\n'.join(input_lst)
        if gaus_input[-2:] != '\n':
            gaus_input += '\n'

        #write or return
        if f:
            f.write(gaus_input)
            f.close()
        else:
            return gaus_input


class TS_input(Input):
    """Class to create qst2 TS optimization.
    Description, charge, multiplicity and geometry are now tuples, not strings!"""
    def __init__(self):
        super(TS_input, self).__init__()
        self.descript = ['', '']
        self.charge = ['', '']
        self.mult = ['', '']
        self.geom = ['', '']

    def set_description(self, descript, number):
        """Number is 1 and 2 and indicates, weather property belongs to
        1st or 2nd geometry."""
        self.descript[number - 1] = descript

    def set_charge(self, charge, number):
        """Number is 1 and 2 and indicates, weather property belongs to
        1st or 2nd geometry."""
        self.charge[number - 1] = str(charge)

    def set_mult(self, mult, number):
        """Number is 1 and 2 and indicates, weather property belongs to
        1st or 2nd geometry."""
        self.mult[number - 1] = str(mult)

    def set_geometry(self, geom, number):
        """Number is 1 and 2 and indicates, weather property belongs to
        1st or 2nd geometry."""
        self.geom[number - 1] = geom

    def provide_chrg_mult_geom(self, geom_with_charge_and_mult, number):
        """Provide geometry with charge and multiplicity attached on top of it.
        It can be easily extracted, if gaussian input is split using \n\n and
        3 element of list is selected.
        Number indicates, weather provided chrg_mult_geom are from the first
        or second molecule.
        """
        chrg_mult_geom_list = geom_with_charge_and_mult.split('\n')
        chrg_mult_list = chrg_mult_geom_list[0].split()
        self.set_charge(chrg_mult_list[0], number)
        self.set_mult(chrg_mult_list[1], number)
        geom_list = chrg_mult_geom_list[1:]
        self.set_geometry('\n'.join(geom_list), number)

    def provide_existing_input_file(self, old_input,
                                    additional_info=False):
        """
        Extracts info for both geometries.
        """
        gaus_lst = old_input.split('\n\n')
        self.provide_first_section(gaus_lst[0])
        self.set_description(gaus_lst[1], 1)
        self.set_description(gaus_lst[1], 2)
        self.provide_chrg_mult_geom(gaus_lst[2], 1)
        self.provide_chrg_mult_geom(gaus_lst[4], 2)
        custom_basis_section_number = 5
        if additional_info:
            custom_basis_section_number += 1
            self.add_additional_info(gaus_lst[5])
        if 'GEN' in self.kwards.upper() and 'GENC' not in self.kwards.upper():
            self.custom_basis = gaus_lst[custom_basis_section_number]
            if 'GENECP' in self.kwards.upper() or 'PSEUDO=READ' in self.kwards.upper():
                self.custom_basis = self.custom_basis + '\n\n' + gaus_lst[custom_basis_section_number + 1]

    def provide_pymol(self, pymol, number):
        """Class will use pymol to extract multiplicity,
        charge and geometry."""
        gaus_str = pymol.write('gau')
        gaus_lst = gaus_str.split('\n\n')
        self.provide_chrg_mult_geom(gaus_lst[2], number)

    def write(self, f=None):
        """Prepares gausian input file. If no argument is supplied,
        returns it as string. Otherwise writes it into given file
        """
        #Create input string
        input_lst = []
        input_lst.append(self.memory)
        input_lst.append(self.proc)
        input_lst.append(self.chk)
        input_lst.append(self.kwards)
        input_lst.append('')
        input_lst.append(self.descript[0])
        input_lst.append('')
        input_lst.append(self.charge[0] + ' ' + self.mult[0])
        input_lst.append(self.geom[0])
        input_lst.append('')
        input_lst.append(self.descript[1])
        input_lst.append('')
        input_lst.append(self.charge[1] + ' ' + self.mult[1])
        input_lst.append(self.geom[1])
        input_lst.append('')
        if self.add_info:
            input_lst.append(self.add_info)
            input_lst.append('')
        if self.custom_basis:
            input_lst.append(self.custom_basis)
            input_lst.append('')
        gaus_input = '\n'.join(input_lst)
        if gaus_input[-2:] != '\n':
            gaus_input += '\n'
        #write or return
        if f:
            f.write(gaus_input)
            f.close()
        else:
            return gaus_input


class CustomBasis(object):
    """Class for creating and writing custom basis sets.

    Takes as input basis set string, downloaded from https://bse.pnl.gov/
    (ESML Basis Set Exchange), in Gaussian format. Then, parses it into
    dictionary in format element symbol : parameters. For example:
    {'H' :      0 \nS   3   1.00..., }

    Can save created dictionary in json format.
    """
    def __init__(self, basis_text_path):
        """Parses supplied text file path and creates basis dictionary.
        """
        self.basis_text_path  = basis_text_path
        self.basis_dic = {}
        self._parse()
        self.abs_path = ''

    def _parse(self):
        with open(self.basis_text_path, 'rb') as f:
            basis_string = f.read()
        atom_list = basis_string.split('****\n')
        for atom_str in atom_list:
            strings = atom_str.split(' ', 1)
            atom = strings[0]
            self.basis_dic[atom] = atom_str

    def write_json(self, filename, path=''):
        """Writes basis dictionary in json format. By default saves
        it in the working directory.
        """
        self.abs_path = os.path.join(path, filename)
        with open(self.abs_path, 'wb') as f:
            json.dump(self.basis_dic, f)
            
            
class ECPBasis(CustomBasis):
    """Class for creating ECP basis sets"""
    def _parse(self):
        with open(self.basis_text_path, 'rb') as f:
            element = None
            param_str = None
            for line in f:
                line_list = line.split()
                if len(line_list) == 2 and line_list[1].isdigit() and \
                    not line_list[0].isdigit() and len(line_list[0]) < 3:
                        #new element parameters string has started
                        if element:
                            self.basis_dic[element] = param_str
                        element = line_list[0]
                        param_str = line
                else:
                    param_str += line
            self.basis_dic[element] = param_str            


def get_unique_non_TS_and_TS_kwrds_lists(com_names):
    """Return list of unique kwards, used in given com_filenames iterable."""
    non_TS_kwards = []
    TS_kwards = []
    for com in com_names:
        with open(com, 'rb') as f:
            for line in f:
                if line.startswith('#'):
                    if 'TS' in com:
                        TS_kwards.append(line.strip())
                    else:
                        non_TS_kwards.append(line.strip())
    return list(set(non_TS_kwards)), list(set(TS_kwards))


def create_energy_dic(input_log_names):
    """ Returns dictionary with ZPE, enthalpy and gibbs energies of given
    finished log files.
    Dictionary has format: mol_name : (ZPE, H, G)
    Energyies are in kcal/mol"""
    energy_dic = {}
    for log in input_log_names:
        with open(log, 'rb') as f:
            log_txt = f.read()
        ZPE_abs = re.findall('Sum of electronic and zero-point Energies=\s+(-\d+\.\d+)', log_txt)[0]
        H_abs = re.findall("Sum of electronic and thermal Enthalpies=\s+(-\d+\.\d+)", log_txt)[0]
        G_abs = re.findall("Sum of electronic and thermal Free Energies=\s+(-\d+\.\d+)", log_txt)[0]
        energy_dic[log] = (float(ZPE_abs)*KCAL_IN_HARTREE, float(H_abs)*KCAL_IN_HARTREE, float(G_abs)*KCAL_IN_HARTREE)
    return energy_dic


def is_finished(log_file):
    """Returns, weather job finished successfully."""
    with open(log_file, 'rb') as f:
        log_file_txt = f.read()
    log_file_lines = log_file_txt.split('\n')
    return log_file_lines[-2][:19] == ' Normal termination'
        
        
def get_calc_time(file_path, filename):
    """Accepts log file path and name.
    If job ended successfully, returns string with job length,
    otherwise returns False."""
    with open(file_path + '/' + filename, 'rb') as f:
        log_file_txt = f.read()
    log_file_lines = log_file_txt.split('\n')
    if log_file_lines[-2][:19] == ' Normal termination':
        time_list = re.findall('(\d+)', log_file_lines[-4])
        return time_list[0] + " d " + time_list[1] + " h " + time_list[2] + \
                " min"
    else:
        return False


def get_NImag(log_str):
    """Returns int with the number of NImags in terminated log file with
    calculated frequeincies.
    If no NImag's were found, returns -1"""
    NImag_lst = re.findall('NImag=(\d+)', log_str)[0]
    if len(NImag_lst) == 0:
        return -1
    else:
        return int(NImag_lst[0])


def is_ECP(element_list, ECP_path="/storage/a92549/data/Def2-SVP_pseudo.json"):
    """Returns, weather any of given molecules have pseudo potentials."""
    with open(ECP_path) as dic:
        ECP_basis_dic = json.load(dic)
    for element in element_list:
        if element.upper() in ECP_basis_dic:
            return True
    return False











