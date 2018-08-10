#!/usr/bin/python
import sys
import os
import pybel
import csv
import operator
import argparse

"""
ver 1.1
9.12.2012

Creates or updates .csv file with sorted molecules energies.

"""

def process_command_line(argv):
    """
    Processes arguments and returns namespace of them
    """
    parser = argparse.ArgumentParser(description="""Creates or updates
                                    csv file with sorted molecule energies.""")
    #Positional args
    parser.add_argument('files', metavar='molec1.ext',
                        nargs='+', help="""Input.""")
    #Optional args
    parser.add_argument('-f', '--format', help="""Babel type of input
                        molecules (default is mopout).""",
                        default='mopout')
    parser.add_argument('-u', '--update', help="""Update mode, as argument will
                        take existing csv file and update it.""",
                        action='store_true')
    return parser.parse_args(argv)

def write_dic_to_csv(dic):
    """takes as input dictionary with filenames and molecule energies
        writes it to the csv file"""
    #created csv filename is same as directory
    csv_name = os.getcwd().split("/")[-1] + ".csv"
    f = open(csv_name, "w")
    writer = csv.writer(f)
    for name, energy in sorted(dic.items(), key=operator.itemgetter(1)):
        writer.writerow([name, energy])
    f.close()
    return csv_name
    
def main(argv):
    args = process_command_line(argv)
    if not args.update:
        ener_dic = {}
        for f in args.files:
            molec = pybel.readfile(args.format, f).next()
            ener_dic[f] = molec.energy
        f_name = write_dic_to_csv(ener_dic) #Create CSV
        return f_name
    else:
        csv_f = args.files[0]
        f = open(csv_f, "r")
        f_temp = open(csv_f +"temp", "w")
        reader = csv.reader(f)
        writer = csv.writer(f_temp)
        #loop over all existing rows
        for row in reader:
            if os.path.exists(row[0]):
                writer.writerow(row)
        #closing and renaming files, which left
        f.close()
        os.unlink(csv_f)
        f_temp.close()
        os.rename(csv_f+"temp", csv_f)
        return csv_f

if __name__ == "__main__":
    main(sys.argv[1:])
