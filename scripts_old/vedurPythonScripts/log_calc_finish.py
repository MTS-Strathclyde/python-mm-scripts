#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Created on Tue Feb 12 14:57:39 2013

@author: a92549
"""

import os
import re
import sys
import csv
import datetime

START_CALCS_LOG = os.environ['START_LOGS_PATH']
FIN_CALCS_LOG = os.environ["FIN_LOGS_PATH"]
CALCS_SUMMRY_LOG = os.environ["SMRY_CSV_PATH"]


def get_NImag(log_str):
    """Returns int with the number of NImags in terminated log file with
    calculated frequeincies.
    If no NImag's were found, returns -1"""
    NImag_lst = re.findall('NImag=(\d+)', log_str)
    if len(NImag_lst) == 0:
        return -1
    else:
        return int(NImag_lst[0])

def NImag_check(log_str, right_number):
    """right_number should be an int.
    Returns empty string if NImag = right_number
    Else, returns string with amount of imaginary frequencies found and warning"""
    NImag = get_NImag(log_str)
    if NImag == right_number:
        return ''
    elif NImag == -1:
        return ''
    else:
        return "Warning, Nimag={0}".format(NImag)


def get_datetime():
    """Return string with current datetime in Y-M-D H:Min format.
    """
    return datetime.datetime.now().strftime("%Y-%m-%d %H:%M")


def get_calc_time(log_file_txt):
    """ If job ended successfully, returns string with job length,
    otherwise returns None."""
    log_file_lines = log_file_txt.split('\n')
    if log_file_lines[-2][:19] == ' Normal termination':
        time_list = re.findall(r'(\d+)', log_file_lines[-4])
        return time_list[0] + " d " + time_list[1] + " h " + time_list[2] + \
                " min"
    else:
        return None

def log_into_fin_log_file(file_path, filename):
    """Write line about finished calculation in FIN_CALC_FILE.
    Returns writen line."""
    log_file = open(FIN_CALCS_LOG, 'ab')
    log_line = '{0};' + get_datetime() + ';{1};{2};{3}'
    time = False
    try:
        with open(file_path + '/' + filename + '.log', 'rb') as f:
            log_file_txt = f.read()
        time = get_calc_time(log_file_txt)
    except IOError:
        log_line = log_line.format(filename, 'Output log file not found', 'F', '\n')
    except IndexError:
        log_line = log_line.format(filename, 'Output file is too short', 'F', '\n')
    if time:
        if 'TS' in filename:
            NImag_str = NImag_check(log_file_txt, 1)
        else:
            NImag_str = NImag_check(log_file_txt, 0)            
        log_line = log_line.format(filename, time, 'x', NImag_str + '\n')
    else:
        log_line = log_line.format(filename, 'Calculation did not finish', 'F', '\n')        
    log_file.write(log_line)
    log_file.close()
    return log_line


def log_into_calc_smry(filename, fin_log_line):
    """Adds information about finished calculation into CALC_SMRY_LOG
    """
    start_log_line = False
    with open(START_CALCS_LOG, 'rb') as f:
        for line in reversed(f.readlines()):
            if line.startswith(filename):
                start_log_line = line
                break
    with open(CALCS_SUMMRY_LOG, 'ab') as f:
        wr = csv.writer(f, dialect=csv.excel)
        fin_log_line_list = fin_log_line.split(';')
        fin_log_line_list[-1] = fin_log_line_list[-1].strip('\n')
        if start_log_line:
            start_log_line_list = start_log_line.split(';')
            start_log_line_list[-1] = start_log_line_list[-1].strip('\n')
            sum_log_line_list = start_log_line_list + fin_log_line_list[1:]
        else:
            sum_log_line_list = ['Start line was not found'] + fin_log_line_list
        wr.writerow(sum_log_line_list)
    
def main(argv):
    reslts_path = argv[0]
    jobname = argv[1]              #No extension!
    fin_log_line = log_into_fin_log_file(reslts_path, jobname)
    log_into_calc_smry(jobname, fin_log_line)


if __name__ == '__main__':
    main(sys.argv[1:])

