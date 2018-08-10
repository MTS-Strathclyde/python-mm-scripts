#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  3 18:43:28 2012

@author: a92549

sends file for gaussian calculation through sipsik

Right now adds arguments:
B3LYP Opt=Tight Freq Int=UltraFine Gen;
uses Def2-SVP basis and saves checkpoint file in same folder, with same
filename but .chk extension.

"""

import argparse
import fileutil
import molutil
import gaussian
import os
import sys
import sipsik
import simplog

def process_command_line(argv):
    parser = argparse.ArgumentParser(description="""Send molecules or folder with
                                    MOPAC output to sipsik for gaussian 
                                    calculation.
                                    
~~~!!!DOESNT WORK ANYMORE AFTER Rg98 CHANGE!!!~~~
                                    
                                    
                                    """)
    #Positional args
    parser.add_argument('files', metavar='molec1.ext',
                        nargs='+', help="""Input files or folder, if MOPAC
                        mode is selected.""")
    #Optional arguments
    parser.add_argument('-m', '--mopac', help="""MOPAC mode, as argument will
                        take folder, with output MOPAC files.""",
                        action='store_true')
    parser.add_argument('-f', '--format', help="""Babel type of input
                        molecules (default is mopout).""",
                        default='mopout')
    parser.add_argument('-d', '--out_folder', help="""Output folder name (default =
                        current working directory name + gaussian.)""")
#    parser.add_argument('-k', '--kwards', help="""Gaussian arguments. Should
#                        be surounded by double quotes.
#                        (B3LYP Opt=Tight Freq Int=UltraFine Gen)""",
#                        default="B3LYP Opt=Tight Freq Int=UltraFine Gen")
#    parser.add_argument('-b', '--basis', help="""In case of custom basis
#                        set, specify location of json formated dictionary,
#                        which contains this basis.""")
#                        (/storage/a92549/data/Def2-SVP.json),
#                        default="/storage/a92549/data/Def2-SVP.json")
#    parser.add_argument('-c', '--check', help="""Checkpoint file filename 
#                        (by default is the same, as input filename, but with chk
#                        extension).""")
    return parser.parse_args(argv)
    
def setup_converter(kwards, basis_path):
    """Creates converter object, which is used to prepare gaussian input
    files"""
    converter = gaussian.Input(kwards)
    converter.set_custom_basis(basis_path)
    return converter
    
def setup_sipsik():
    """Creates sipsik connection object and connects to sipsik."""
    sips = sipsik.Connection()
    sips.connect()
    return sips

def setup_io(args, logger):
    """Get files and setup output folder"""
    if args.mopac:
        out_files = fileutil.list_ext('.out', args.files[0])
        pymols_dic = molutil.create_pymol_dic('mopout', out_files)
        logger.input_directory = args.files[0]
    else:
        pymols_dic = molutil.create_pymol_dic(args.format, args.files)
        logger.input_directory = os.getcwd()     #Fails, if started from different dir.
    if args.out_folder:
        out_folder = args.out_folder
    else:
        out_folder = os.path.split(os.getcwd())[1] + "_gaussian"
    out_folder = os.path.abspath(out_folder)
    os.mkdir(out_folder)
    logger.output_directory = out_folder
    return out_folder, pymols_dic
    
def write_mols(out_folder, pymols_dic, kwards, basis_path):
    """Write molecules to output folder
    """
    gaus_input_pathes = []
    converter = setup_converter(kwards, basis_path)
    for f_name, pymol in pymols_dic.iteritems():
        #prepare file to write into
        g_input_name = fileutil.change_ext(f_name, '.com')
        g_input_path = os.path.join(out_folder, os.path.split(g_input_name)[1])
        gaus_input_pathes.append(g_input_path)
        #write there
        converter.set_molecule(pymol)
        g_input_file = open(g_input_path, 'wb')
        converter.write(g_input_file)
    return gaus_input_pathes

    
def main(argv):
    """Proceeds through 3 stages.
    Firstly it collects all specified files, converts them into 
    pybels dictionary and sets up output folder.
    
    Then it writes gaussian input molecules into this folder with right 
    keywords and custom basis.
    
    Finally, it invokes rg98 script on sipsik, thus sending gaussian input
    files for calculation.
    """
    args = process_command_line(argv)
    logger = simplog.Logger()
    kwards = "B3LYP Opt=Tight Freq Int=UltraFine Gen"
    basis_path = '/storage/a92549/data/Def2-SVP.json'
    out_folder, pymols_dic = setup_io(args, logger)
    gaus_input_pathes = write_mols(out_folder, pymols_dic, kwards, basis_path)
    print "Sending files for calculation on sipsik."
    print logger.datetime()
    print logger.io_information()
    sips = setup_sipsik()
    for g_input_path in gaus_input_pathes:
        sip_path = sips.change_path(g_input_path)
        split_sip_path = os.path.split(sip_path)
        print split_sip_path[1]
        sip_command = 'cd ' + split_sip_path[0] + ' && rg98 ' + split_sip_path[1]
        sips.execute_command(sip_command)
    sips.disconnect()
        
if __name__ == '__main__':
    main(sys.argv[1:])
