#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 19 09:54:27 2013

@author: a92549
"""

import time

FIN_CALCS_LOG = '/storage/a92549/log/gaus_fin.log'
START_CALCS_LOG = '/storage/a92549/log/gaus_strt.log'

def make_job_name_discription_list(f):
    """Creates list of the form:
        [(molecule_name, molecule_related input), ..]"""
    job_list = []
    for line in f:
        strings = line.split(';')
        job_list.append((strings[0], ' '.join(strings)))
    return job_list
    
    
def main():
    with open(START_CALCS_LOG, 'rb') as f:
        strted_job_list = make_job_name_discription_list(f)
    with open(FIN_CALCS_LOG, 'rb') as f:
        finished_job_dict = dict(make_job_name_discription_list(f))
    total = 0
    for strted_mol, descript in strted_job_list:
        try:
            finished_job_dict[strted_mol]
        except KeyError:
            print "Running: {0}".format(descript.strip())
            total += 1
    print "{0} jobs are running.".format(total)

if __name__ == '__main__':
    a = time.time()
    main()    
    b = time.time()
    print b - a