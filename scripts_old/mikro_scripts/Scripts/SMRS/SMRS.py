#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  4 10:47:13 2012

@author: mishin1991
                    UNFINISHED DESCRIPTION
    Representative objects.

    This script takes as input filenames of MOPOUT conformers and tries
    to select few representative conformers from given set.

"""

import argparse
import sys
import pybel
import openbabel
import molecule_grouping
import csv
import subprocess
import os
import shutil

RUN_PY_SMRS_PATH = '/home/mishin1991/Tools/SMRS_v1.0/run_py_smrs.m'
FILE_WITH_REP_MOL_IDXS = "repInd_temp.csv"

def process_command_line(argv):
    """
    Processes arguments and returns namespace of them
    """
    parser = argparse.ArgumentParser(description='Run SMRS algorithm.')
    #Positional args
    parser.add_argument('files', metavar='molec1.ext',
                        nargs='+', help='files to be processed')
    #Optional true/false args
    parser.add_argument('-d', '--dihed',
                        help='Use dihedral angles as coordinates',
                        action='store_true')
    parser.add_argument('-n', '--nonH',
                        help='Do not include hydrogen in coords',
                        action='store_true')
    parser.add_argument('-e', '--energy', help='Append energy at the end of \
                        molecules vector', action='store_true')
    parser.add_argument('--delCoordCSV', help='Delete CSV file with molecule \
                        coordinates', action='store_true')
    parser.add_argument('--delCoefCSV', help='Delete C matrix CSV file',
                        action='store_true')
    parser.add_argument('--folder', help='Name of folder, where representative\
                        molecules will be saved', action='store_true')
    parser.add_argument('--sdf', help='If you want to calculate pybel\
                        fingerprints', action='store_true')
    #Optional args with value
    parser.add_argument('--alpha', help='Specify lambda paramter',
                        type=int, metavar='A')
    parser.add_argument('--division', help='Specify type and parametrs\
                        of molecule division into groups. Parametrs should\
                        be separated by comma without spaces', metavar='name')
    parser.add_argument('--pruneStart', help='Specify minimum RMSD by which \
                        molecules should be separated in starting set',
                        type=float, metavar='RMSD')
    parser.add_argument('--pruneFinish', help='Specify minimum RMSD by which \
                        representive molecules should be separated',
                        type=float, metavar='RMSD')
    parser.add_argument('--format', help='Babel type of input molecules')
    return parser.parse_args(argv)

def convert_to_pybel(file_list, mol_format):
    """converts list of filenames into pybel objects"""
    if not mol_format:
        mol_format = "mopout"
    return [pybel.readfile(mol_format, name).next() for name in file_list]

def job_string(args):
    """summary of paramters specified for this run in a string"""
    name = ''
    if args.pruneStart:
        name += str(args.pruneStart)
    if args.dihed:
        name += '_dihed_'
    else:
        name += 'xyz_'
    if args.nonH:
        name += 'nonH_'
    else:
        name += 'H_'
    if args.division:
        name += args.division
    if args.pruneFinish:
        name += str(args.pruneFinish)
    return name

def prune(pybel_list, min_RMSD):
    """
    iterates over all pairs and deletes second molecule in pair if RMSD
    between paired molecules is less then min_RMSD
    """
    #Set up OBAling object
    align = openbabel.OBAlign()
    #Loop
    i = 0
    total_removed = 0
    while i < len(pybel_list):
        referens = pybel_list[i].OBMol  #reference
        align.SetRefMol(referens)
        j = i + 1
        while j < len(pybel_list):
            target = pybel_list[j].OBMol #target
            align.SetTargetMol(target)
            #Align and ret rmsd
            if align.Align():
                rmsd = align.GetRMSD()
                if rmsd < min_RMSD:
                    pybel_list.pop(j)   #remove from both lists
                    total_removed += 1
                else:
                    j = j + 1
            else:
                print "Couldn't align"
                raise Exception()
        #end of inner loop
        i = i + 1
    #end of outer loop
    print "finished deleting, total number of \
            removed conformers is", total_removed
    return pybel_list

def call_octave(coord_filename, c_matrix_filename, alpha=5):
    """calls rup_py_smrs.m octave scritpt"""
    subprocess.call([RUN_PY_SMRS_PATH, str(alpha), coord_filename,
                     c_matrix_filename])

def get_rep_mol_indexes():
    """returns representative molecule indexes (starting from 0) and
    deletes csv file, which contains them"""
    f = open(FILE_WITH_REP_MOL_IDXS, "r")
    rd = csv.reader(f)
    mols = rd.next()
    f.close()
    mol_idxs = [int(i) - 1 for i in mols]
    os.unlink(FILE_WITH_REP_MOL_IDXS)
    return mol_idxs

def analyze_coords(c_matrix_filename):
    row_sums = []
    tot_sum = 0    
    f = file(c_matrix_filename, 'r')
    rdr = csv.reader(f)
    for line in rdr:
        pos_floats = map(abs, map(float, line))
        row_tot = sum(pos_floats)
        row_sums.append(row_tot)
        tot_sum += row_tot
    percent = lambda x : x/tot_sum*100
    weights = map(percent, row_sums)
    return weights
    
def vector(molec, dihed, nonH, energy):
    """
    Creates single molecule vector
    reurns it as tuple
    """
    #Torison
    if dihed:
        pass
    #XYZ
    else:
        coords = ()
        if nonH:
            for atom in molec.atoms:
                coords += atom.coords
        else:
            for atom in molec.atoms:
                if atom.atomicnum > 1:
                    coords += atom.coords
    #Energy
    if energy:
        coords += (molec.energy/10.0,)
    return coords

def coord_file(pybel_group, dihed, nonH, energy, name):
    """Creates coord file for single group"""
    csv_name = "coords" + name + ".csv"
    #Open file, create writer
    f  = open(csv_name, "w")
    wr = csv.writer(f)
    #Generate coords and write them
    for py_molec in pybel_group:
        wr.writerow(vector(py_molec, dihed, nonH, energy))
    f.close()
    return csv_name

def run_smrs(grouped_pybels, dihed, nonH, energy, alpha, delCoordCSV,
             delCoefCSV, name):
    """runs SMRS algorithm
    delete csv-s, if specified
    return list of lists. Each inner list contains filenames of representative
    molecules of particular group
    """
    all_weights = []
    groups_reps = []
    group_num = 0
    for group in grouped_pybels:
        #name manipulations
        group_num += 1
        name += str(group_num)
        coord_filename = coord_file(group, dihed, nonH, energy, name)
        c_matrix_filename = "c_matrix_" + name + ".csv"
        #call octave
        if alpha:
            call_octave(coord_filename, c_matrix_filename, alpha)
        else:
            call_octave(coord_filename, c_matrix_filename)
        #get reps
        group_idxs = get_rep_mol_indexes()
        #get their weights
        all_weights.append(analyze_coords(c_matrix_filename))
        #delete csv-s
        if delCoordCSV:
            os.unlink(coord_filename)
        if delCoefCSV:
            os.unlink(c_matrix_filename)
        #save filenames
        group_reps = [group[idx] for idx in group_idxs]
        groups_reps.append(group_reps)
    return groups_reps, all_weights

def handle_sdf(args):
    #load sdf
    mols = list(pybel.readfile("sdf", args.files[0]))
    #create sdf file for fps
    f = open("fps.csv", "w")
    wr = csv.writer(f)
    #safety
    i = 0
    #write fingerprings
    for mol in mols:
        fp = list(mol.calcfp().fp)
        wr.writerow(fp)
        i += 1
        if i > 500:
            break
    #run SMRS and save
    call_octave("fps.csv", "fps_c_matrix.csv")
    group_idxs = get_rep_mol_indexes()
    for idx in group_idxs:
        mols[idx].write("xyz", str(idx) + ".xyz")
    print "done"

def main(argv):
    """
    main function of SMRS algorithm, parses arguments using argparse
    for details: run in terminal with argument --help"""
    args = process_command_line(argv)
    name = job_string(args)
    #That feel when no torison ;_;
    if args.dihed:
        raise Exception("Dihed is not supported right now")
    #SDFS!
    if args.sdf:
        handle_sdf(args)
    #Conversion, pruning
    pybel_mols = convert_to_pybel(args.files, args.format)
    if args.pruneStart:
        pybel_mols = prune(pybel_mols, args.pruneStart)
    print "Total number of molecules to process is", len(pybel_mols)
    #Division
    if args.division:
        grouped_pybels = molecule_grouping.main(args.division, pybel_mols)
    else:
        grouped_pybels = [pybel_mols]
    #Run algorithm
    groups_reps, weights = run_smrs(grouped_pybels, args.dihed, args.nonH, args.energy,
                           args.alpha, args.delCoordCSV, args.delCoefCSV, name)
    prune_finished = False
    #Pruning representatives
    if args.pruneFinish:
        all_reps = []
        for group in groups_reps:
            all_reps += group
        all_reps = prune(all_reps, args.pruneFinish)
        prune_finished = True
    #Save all groups into one folder
    folder_name = 'rep_' + name
    if args.folder:
        #folder creation
        while True:
            if not os.path.exists(folder_name):
                os.mkdir(folder_name)
                break
            else:
                folder_name = folder_name + 'c'
        #copying
        if prune_finished:
            for mol in all_reps:
                shutil.copy(mol.title, os.getcwd() + "/" + folder_name)
        else:
            for group in groups_reps:
                for mol in group:
                    shutil.copy(mol.title, os.getcwd() + "/" + folder_name)
    print "Coeficient matrix results"
    for i in range(len(grouped_pybels)):
        for j in range(len(grouped_pybels[i])):
            print grouped_pybels[i][j].title, weights[i][j]
    print ""
    print "Rep mols"
    for group in groups_reps:
        for mol in group:
            print mol.title
    return groups_reps

if __name__ == '__main__':
    main(sys.argv[1:])