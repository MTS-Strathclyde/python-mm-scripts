#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 14 16:58:43 2013

@author: a92549

Prepares csv summary about relative energies of citronellals and isopulegols
in folder. Input should be completed log files.
"""

import argparse
import sys
import csv
import gaussian
import os
import simplog
import math

RESULTS_REL_TO = 'cit_min*'
KCAL_GAS_CONST = 1.9858/1000

def process_command_line(argv):
    """
    Processes arguments and returns namespace of them
    """
    parser = argparse.ArgumentParser(description="""Prepares solvent calculation
                                    for TS files. Input is gas phase input file.
                                    """)
    #Positional args
    parser.add_argument('input', metavar='mol.log', nargs='+',
                        help="""Calculated molecules.""")
    #Optional args
    parser.add_argument('-t', '--temp', help='Temperature of thermodynamic output (298.15)',
                        default=298.15, type=float)
    parser.add_argument('-n', '--name', help="""Name of summary (by default is the
                        same as folder name, but with _summary suffix.""")
    return parser.parse_args(argv)


def write_header(csv_wr, completed_log_files, temp):
    for i in range(3):
        csv_wr.writerow([''])
    csv_wr.writerow(['Date created:', simplog.get_datetime()])
    csv_wr.writerow(['Description:'])
    csv_wr.writerow(['Path:', os.getcwd()])
    com_names = (log[:-4] + '.com' for log in completed_log_files)
    non_TS_kwards, TS_kwards = gaussian.get_unique_non_TS_and_TS_kwrds_lists(com_names)
    csv_wr.writerow(['Kwards used in non TS calcs:'] + non_TS_kwards)
    csv_wr.writerow(['Kwards used in TS calcs:'] + TS_kwards)
    csv_wr.writerow('')
    csv_wr.writerow(('All given energy is in Kcal/mol',))
    csv_wr.writerow(('Temperature:', temp))
    csv_wr.writerow('')


def write_energy_comparison(csv_wr, energy_dic):
    """Write energies for every mol in dictionary"""
    def write_energy_line(mol):
        """Write single row with energies."""
        zpe, h, g = energy_dic[mol]
        row = (mol, zpe, zpe - rel_zpe, h, h - rel_h, g, g - rel_g)
        csv_wr.writerow(row)
    csv_wr.writerow(('Compound', 'Abs. ZPE', 'Rel ZPE', 'Abs. H', 'Rel. H', 'Abs. G', 'Rel. G'))
    #Write relative row
    for mol in energy_dic.iterkeys():
        if mol.startswith(RESULTS_REL_TO[:-1]):
            rel_zpe, rel_h, rel_g = energy_dic[mol]
            write_energy_line(mol)
            energy_dic.pop(mol)
            break
    #Sort by gibbs and write values
    get_gibbs = lambda x: x[1][2]
    for mol, _ in sorted(energy_dic.iteritems(), key=get_gibbs):
        write_energy_line(mol)
        

def write_percentages(csv_wr, energy_dic, temp):
    """Works only for cit iso min TS naming style.
    Uses gibbs free energy."""
    def write_percentage(list_):
        min_energy = min([energy for _, energy in list_])
        boltz_list = [math.exp(-(energy-min_energy)/temp/KCAL_GAS_CONST) for _, energy in list_]
        Z = sum(boltz_list)
        csv_wr.writerow([mol_name for mol_name, _ in list_])
        csv_wr.writerow([p_i/Z*100 for p_i in boltz_list])
    min_list = []
    TS_list = []
    for mol in energy_dic:
        if 'min' in mol and 'cit' not in mol:
            min_list.append((mol, energy_dic[mol][2]))
        elif 'TS' in mol:
            TS_list.append((mol, energy_dic[mol][2]))
    csv_wr.writerow('')
    csv_wr.writerow(('Thermodynamic distribution',))
    write_percentage(min_list)
    csv_wr.writerow('')
    csv_wr.writerow(('Kinetic distribution',))
    write_percentage(TS_list)


def main(argv):
    args = process_command_line(argv)
    completed_log_files = []
    for f_name in args.input:
        if gaussian.is_finished(f_name):
            if 'IRC' not in f_name:
                completed_log_files.append(f_name)
    energy_dic = gaussian.create_energy_dic(completed_log_files)
    if args.name:
        csv_name = args.name
    else:
        csv_name = os.getcwd().split('/')[-1] + '_summary.csv'
    with open(csv_name, 'wb') as csv_f:
        csv_wr = csv.writer(csv_f)
        write_header(csv_wr, completed_log_files, args.temp)
        write_energy_comparison(csv_wr, energy_dic)
        write_percentages(csv_wr, energy_dic, args.temp)

if __name__ == '__main__':
    main(sys.argv[1:])
