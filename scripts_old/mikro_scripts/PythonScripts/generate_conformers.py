#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 22 16:56:29 2012

@author: Maksim Mišin

Conformer generating script.

Takes molecule in sdf format and writes conformers into sdf file.
"""

import os
import subprocess
import argparse
import shutil
import sys
import sdfSize
import time

#Conformer generating program names and options
BALLOON = ["balloon", "--nGenerations", "3000", #number of generations
            "--nconfs", "1000", #numer of produced conformers
            "--RMSDtol", "0.3", #RMSD value for pruning
            "--pStereoMutation", "0.0", #don't change stereochemistry
            "--pRingFlipMutation", "4.0", #flip rings wiht this possibility
            "--tournamentSize", "3", #increase diversity
            "--addConformerNumberToName" #number conformers
            ]
MARVIN = ["cxcalc", "-Xmx1024m", #More memory for java
          "conformers", "-m", "500", #number of conformers to produce
            "-d", "0.3", #RMSD value for pruning
            ]

def process_command_line(argv):
    parser = argparse.ArgumentParser(description='Generate conformers')
    #Positional args
    parser.add_argument('input', help='Input filename', metavar='<inputmol>')
    parser.add_argument('output', help='Output filename', metavar='<outputmol>')
    #Optional arguments
    #TODO!
    return parser.parse_args(argv)


class ConformerGenerator:
    """Class, which is used to create conformer generation objects from
    existing programs.
    """
    def __init__(self, in_file, options, redirect=False):
        """Initialization requires input filename, options list
        and redirect boolean.
        If redirect is False, will use following syntax to call program:
            options in_file out_file
        Otherwise:
            options in_file > out_file
        """
        self.in_file = in_file
        self.out_file = in_file + options[0] + "_temp.sdf"
        self.redirect = redirect
        self.options = options

    def generate_conformers(self):
        """Produces conformers and returns output filename"""
        command_line_command = []
        command_line_command += self.options
        if self.redirect:
            out_f = file(self.out_file, 'w')
            command_line_command += [self.in_file]
            subprocess.call(command_line_command, stdout=out_f)
            out_f.close()
        else:
            command_line_command.append(self.in_file)
            command_line_command.append(self.out_file)
            #make subprocess silent
            fnull = open(os.devnull, 'w')
            subprocess.call(command_line_command, stdout=fnull)
            fnull.close()
        return self.out_file


def concatenate(individual_program_outputs_list, output):
    """concateanate all files in tuple and delete them in process
    return concatenated file filename"""
    conc_out_file = open(output, 'wb')
    for filename in individual_program_outputs_list:
        shutil.copyfileobj(open(filename,'rb'), conc_out_file)
        os.unlink(filename)
    conc_out_file.close()
    return output

def run_script(argv):
    """returns produced ouput file filename"""
    marvin = ConformerGenerator(argv.input, MARVIN, True)
    balloon = ConformerGenerator(argv.input, BALLOON)
    individual_program_outputs_list = [marvin.generate_conformers(),
                                       balloon.generate_conformers()]
    return concatenate(individual_program_outputs_list, argv.output)
    
def info(args):
    print "Generating conformers using FF"
    print "Input file: " + args.input
    print "Output file: " + args.output
    print "Output directory " + os.getcwd()
    print "Marvin commands: " + " ".join(MARVIN)
    print "Balloon commands: " + " ".join(BALLOON)
    print "Job description: "
    
def main(argv):
    args = process_command_line(argv)
    #initialize programs and call them
    info(args)
    start_time = time.time()
    output_file = run_script(args)
    finish_time = time.time()
    print "Generation finished"
    print "Total number of conformers produced " + str(sdfSize.size(output_file))
    print "Total time: " + str(round(finish_time - start_time, 3))
    return output_file

if __name__ == '__main__':
    main(sys.argv[1:])
