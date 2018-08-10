#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 25 21:35:09 2012

@author: mishin1991

"""

import sys
import argparse
import os
import openbabel
import shutil
import pybel
import datetime
import time

STRICT_ENERGY_DIF = 0.1

def process_command_line(argv):
    """
    Processes arguments and returns namespace of them
    """
    parser = argparse.ArgumentParser(description="""Prunes molecules. By default
                                    will move pruned files to new location""")
    #Positional args
    parser.add_argument('files', metavar='molec1.ext',
                        nargs='+', help="""Input files. If there is only one,
                        script assumes, that it is multimolecular file and
                        tries to prune it, if script gets more then one file
                        as input, it asumes, that those are files with diffrent
                        conformers and prunes them.""")
    #Optional true/false args
    parser.add_argument('--delete', help="""Will delete non-unique files.""",
                        action='store_true')
    parser.add_argument('-m', '--mopac', help="""MOPAC mode, script will also
                        delete/move .arc and .dat files. As argument will
                        take folder, with MOPAC files.""",
                        action='store_true')
    parser.add_argument('-s', '--strict', help="""Will only remove conformers,
                        which differ not only by rmsd, but by energy as well.
                        Default energy difference value - 0.1 kJ/mol.
                        Works only for multiple files  input.""",
                        action='store_true')
        #Removed for now for simplicity
#    parser.add_argument('-d', '--folder',
#                        help="""In folder mode script takes folder as argument
#                        and prunes molecules, which this folder contains.
#                        Works only for non multimolecular files""",
#                        action='store_true')
    #Optional args with value
    parser.add_argument('-o', '--out_filename', help="""Output filename
                        curtainly only for sdf input.""")
    parser.add_argument('-r', '--rmsd', help="""RMSD value, by which molecules,
                        should differ (0.3).""",
                        type=float, metavar='RMSD', default=0.3)
    parser.add_argument('-f', '--format', help="""Babel type of input
                        molecules (default is sdf, if mopac option is specified,
                        then mopout).""",
                        default='sdf')
    parser.add_argument('-d', '--out_folder', help="""Output folder name (default =
                        current working directory name + r + rmsd value)""")
    return parser.parse_args(argv)

def get_mopac_molecules(folder):
    """parses folder with MOPAC calculation and returns list with out files
    and mopac path"""
    if folder[-1] == '/':
        folder = folder[:-1]
    mopac_path = os.path.join(os.getcwd(), folder)
    all_files = os.listdir(mopac_path)
    out_files = []
    for molec in all_files:
        if molec[-4:] == ".out":
            out_files.append(os.path.join(mopac_path, molec))
    return out_files, mopac_path

def get_sdf_molecules(sdf_file, rmsd):
    """returns list of pybel molecules and name of sdf file"""
    sdf_mol_gen = pybel.readfile('sdf', sdf_file)
    sdf_name = os.path.splitext(sdf_file)[0] + '_r' + str(rmsd) + '.sdf'
    return list(sdf_mol_gen), sdf_name

def convert_to_pybel(files, babel_format):
    """Converts given filelist in specified babel format into pybel
    molecules."""
    pybel_list = [pybel.readfile(babel_format, filename).next() for filename in files]
    return pybel_list

def get_filename(path):
    """Retuns filename, remove path and extension."""
    with_ext = os.path.split(path)[1]
    return os.path.splitext(with_ext)[0]    
    

def find_duplicates_and_unique(pybel_list, min_rmsd, stict=False, energ_dif=0.1):
    """return two lists with unique and duplicate molecules in
    pybel format"""
    align = openbabel.OBAlign()
    i = 0
    duplicate_list = []
    #iterate over all
    while i < len(pybel_list):
        align.SetRefMol(pybel_list[i].OBMol)
        j = i + 1
        #iterate over reference pairs
        print "Molecule " + get_filename(pybel_list[i].title) + " has been set."
        while j < len(pybel_list):
            align.SetTargetMol(pybel_list[j].OBMol)
            if align.Align():
                rmsd = align.GetRMSD()
                if rmsd < min_rmsd:
                    if stict:
                        if abs(pybel_list[i].energy - pybel_list[j].energy) < energ_dif:
                            print pybel_list[i].title + " and " + pybel_list[j].title + " are simmilair"
                            print "Enegy difference: " + str(pybel_list[i].energy - pybel_list[j].energy)
                            print "RMSD score: " + str(rmsd)
                            print "Not copying " + pybel_list[j].title
                    duplicate_list.append(pybel_list.pop(j))   
                else:
                    j = j + 1
            else:
                print "Couldn't align"
                raise Exception()
        i += 1
    print "Total number of duplicates: " + str(len(duplicate_list))
    print "Molecules remaining: " + str(len(pybel_list))
    return pybel_list, duplicate_list

def delete_files(duplicate_list, mopac_path=None):
    """removes dublicats and arc and dat files, if necessary"""
    for pymol in duplicate_list:
        filename = pymol.title
        if mopac_path:
            path_out = os.path.join(mopac_path, filename)
            path_arc = os.path.join(mopac_path, filename[:-4] + ".arc")
            path_dat = os.path.join(mopac_path, filename[:-4] + ".dat")
            os.unlink(path_out)
            os.unlink(path_arc)
            os.unlink(path_dat)
        else:
            os.unlink(filename)
    return True

def copy_files(unique_list, out_folder, min_rmsd, mopac_path=None):
    """copy unique files to a new folder"""
    if not mopac_path:
        if not out_folder:
            out_folder = os.path.split(os.getcwd())[1] + "_r" + str(min_rmsd)
        os.mkdir(out_folder)
        for pymol in unique_list:
            filename = pymol.title
            shutil.copyfile(filename, out_folder + "/" + filename)
    else:
        #mopac option
        if not out_folder:
            out_folder = os.path.split(mopac_path)[1] + "_r" + str(min_rmsd)
        os.mkdir(out_folder)
        for pymol in unique_list:
            out_filename = pymol.title.split("/")[-1]
            arc_filename = out_filename[:-4] + ".arc"
            dat_filename = out_filename[:-4] + ".dat"
            shutil.copyfile(mopac_path + "/" + out_filename, out_folder + "/" + out_filename)
            shutil.copyfile(mopac_path + "/" + dat_filename, out_folder + "/" + dat_filename)
            try:
                shutil.copyfile(mopac_path + "/" + arc_filename, out_folder + "/" + arc_filename)
            except:
                print "Couldn't find " + arc_filename + ". Ommiting copying it."
    return out_folder

def write_mols_to_sdf(unique_list, sdf_filename, min_rmsd, out_filename):
    """Creates new sdf file and writes unique molecules into it"""
    if out_filename:
        out_mol = pybel.Outputfile('sdf', out_filename)
    else:
        out_mol = pybel.Outputfile('sdf', sdf_filename)
    for pymol in unique_list:
        out_mol.write(pymol)
    out_filename = out_mol.filename
    out_mol.close()
    return out_filename

def handle_arguments_and_run(args):
    """main function, which analyzes command line input and preforms
    actions according to it"""
    if args.delete:
        print "Will delete non-unique conformers."
    else:
        print "Will copy unique conformers."
    if len(args.files) > 1:
    #bunch of non-mopac files (for example: xyz) handling
        print "File type: " + args.format
        pybel_list = convert_to_pybel(args.files, args.format)
        print "Total number of molecules: " + str(len(pybel_list))
        if args.strict:
            unique_list, duplicate_list = find_duplicates_and_unique(pybel_list, args.rmsd, True, STRICT_ENERGY_DIF)
        else:
            unique_list, duplicate_list = find_duplicates_and_unique(pybel_list, args.rmsd)
        if args.delete:
            return delete_files(duplicate_list)
        else:
            return copy_files(unique_list, args.out_folder, args.rmsd)
    else:
        if args.mopac:
            #mopac handling
            print "File type: MOPAC output"
            out_files, mopac_path = get_mopac_molecules(args.files.pop())
            pybel_list = convert_to_pybel(out_files, "mopout")
            unique_list, duplicate_list = find_duplicates_and_unique(pybel_list, args.rmsd)
            if args.delete:
                return delete_files(duplicate_list, mopac_path)
            else:
                return copy_files(unique_list, args.out_folder, args.rmsd, mopac_path)
        else:
            #multi-molecule sdf handling
            print "File type: multi-conformer sdf file"
            sdf_file = args.files.pop()
            pybel_list, sdf_name = get_sdf_molecules(sdf_file, args.rmsd)
            unique_list, duplicate_list = find_duplicates_and_unique(pybel_list, args.rmsd)
            if args.delete:
                os.unlink(sdf_file)
            return write_mols_to_sdf(unique_list, sdf_name, args.rmsd, args.out_filename)

def main(argv):
    args = process_command_line(argv)
    print "Pruning conformers."
    print "Minimum rmsd difference: " + str(args.rmsd)
    a = time.time()
    output = handle_arguments_and_run(args)
    b = time.time()
    if output:
        print "Finished pruning successfully."
    else:
        print "Error occured, didn't finish: " + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        raise ValueError()
    print "Total time: " + str(round(b - a, 3)) + " seconds."
    return output

if __name__ == '__main__':
    main(sys.argv[1:])
