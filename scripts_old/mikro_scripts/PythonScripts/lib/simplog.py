# -*- coding: utf-8 -*-
"""
Created on Tue Dec  4 14:20:31 2012

@author: a92549

Small logging module.

"""
import datetime
import time
import os
import csv


class Logger(object):
    """Class for returning different logging messages in string format"""
    def __init__(self):
        self._t_start = None
        self._t_finish = None
        self.input_file = None
        self.output_file = None
        self.input_directory = None
        self.output_directory = None

#---------------Time notifications--------------------

    def datetime(self):
        """Return string with current datetime in Y-M-D H:Min format.
        """
        return datetime.datetime.now().strftime("%Y-%m-%d %H:%M")

    def start_stopper(self):
        """ Start stopper.
        """
        self._t_start = time.time()

    def stop_stopper(self):
        """Stop stopper.
        """
        self._t_finish = time.time()

    def get_stopper(self):
        """Return string with elapsed time between stopper start and finish.
        """
        return str(self._t_finish - self._t_start)

#--------------File system notifications--------------------------

    def io_information(self):
        """Returns string with info about input and output.
        """
        info = os.linesep
        if self.input_file:
            info += 'Input file is ' + self.input_file + os.linesep
        if self.input_directory:
            info += 'Input directroy is ' + self.input_directory + os.linesep
        if self.output_file:
            info += 'Output file is ' + self.output_file + os.linesep
        if self.output_directory:
            info += 'Output directroy is ' + self.output_directory + os.linesep
        return info


class CSV_Logger(Logger):
    """Prints logging messages into csv file."""
    def __init__(self, csv_file):
        super(CSV_Logger, self).__init__()
        self.csv_file = csv_file

    def gaussian_line(self):
        """Creates entry for gaussian calculations logging file.
        Format is:
        [self.input_file, calc_type(blank), self.input_directory,self.datetime]
        """
        return [self.input_file, '', self.input_directory, self.datetime()]

    def write(self, line):
        """Appends line to csv file.
        Line must be a list. Each list element goes in separate cell."""
        with open(self.csv_file, 'ab') as f:
            wr = csv.writer(f)
            wr.writerow(line)
            
def get_datetime():
    """Return string with current datetime in Y-M-D H:Min format.
    """
    return datetime.datetime.now().strftime("%Y-%m-%d %H:%M")

