#!/usr/bin/python
import sys
import os
import pybel
import prune_conformers as pc
import datetime
import argparse
import shutil

""" Copies files with energy value bigger then given into separate folder 
"""

def process_command_line(argv):
    parser = argparse.ArgumentParser(description="""Remove molecules with energy
                                    higher, then given.""")
    #Positional args
    parser.add_argument('files', metavar='molec1.ext',
                        nargs='+', help="""Input files or folder, if MOPAC
                        mode is selected.""")
    #Optional arguments
    parser.add_argument('-m', '--mopac', help="""MOPAC mode, script will also
                        move .arc and .dat files. As argument will
                        take folder, with MOPAC files.""",
                        action='store_true', default=True)
    parser.add_argument('-e', '--energy', help="""Max allowed energy of molecules,
                        which will be copied (0).""",
                        type=float, metavar='energy', default=0)
    parser.add_argument('-f', '--format', help="""Babel type of input
                        molecules (default is mopout).""",
                        default='mopout')
    parser.add_argument('-d', '--out_folder', help="""Output folder name (default =
                        current working directory name + _rem_high + energy.)""")
    return parser.parse_args(argv)

def copy_files(low_energy_pybel, out_folder, max_energy, mopac_path=None):
    """copy unique files to a new folder"""
    if not mopac_path:
        if not out_folder:
            out_folder = os.path.split(os.getcwd())[1] + "_rem_high_" + str(max_energy)
        os.mkdir(out_folder)
        for pymol in low_energy_pybel:
            filename = pymol.title
            shutil.copyfile(filename, out_folder + "/" + filename)
    else:
        #mopac option
        if not out_folder:
            out_folder = os.path.split(mopac_path)[1] + "_r" + str(max_energy)
        os.mkdir(out_folder)
        for pymol in low_energy_pybel:
            out_filename = pymol.title.split("/")[-1]
            arc_filename = out_filename[:-4] + ".arc"
            dat_filename = out_filename[:-4] + ".dat"
            shutil.copyfile(mopac_path + "/" + out_filename, out_folder + "/" + out_filename)
            shutil.copyfile(mopac_path + "/" + arc_filename, out_folder + "/" + arc_filename)
            shutil.copyfile(mopac_path + "/" + dat_filename, out_folder + "/" + dat_filename)
    return out_folder

def decide_energies(pybel_list, max_energy):
    """return two lists with low and high energy molecules in
    pybel format"""
    low_energy_pybel = []
    high_energy_pybel = []
    for mol in pybel_list:
        if mol.energy >= max_energy:
            high_energy_pybel.append(mol)
        else:
            low_energy_pybel.append(mol)
    return low_energy_pybel, high_energy_pybel

def main(argv):
    args = process_command_line(argv)
    if args.mopac:
        out_files, mopac_path = pc.get_mopac_molecules(args.files[0])
        pybel_list = pc.convert_to_pybel(out_files, "mopout")
        low_energy_pybel, high_energy_pybel = decide_energies(pybel_list, args.energy)
        out_folder = copy_files(low_energy_pybel, args.out_folder, args.energy, mopac_path)
    else:
        pybel_list = pc.convert_to_pybel(out_files, args.format)
        low_energy_pybel, high_energy_pybel = decide_energies(pybel_list, args.energy)
        out_folder = copy_files(low_energy_pybel, args.out_folder, args.energy)
    #Info
    if out_folder:
        print datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print "Copying molecules with energy lower or equal to: " + str(args.energy)
        print "Total files: " + str(len(pybel_list))
        print "High energy molecules: "  + str(len(high_energy_pybel))
        print "Copied molecules: " + str(len(low_energy_pybel))
        print "Output folder: " + out_folder
    else:
        print "Finished unsucsessfully."
    
if __name__ == '__main__':
    main(sys.argv[1:])








        
