#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 23 11:29:21 2015

@author: max
"""

import argparse
import sys
import matplotlib.pyplot as plt
from itertools import izip_longest
import numpy as np
from scipy.signal import argrelextrema
from scipy.integrate import simps



def process_command_line(argv):
    """Processes arguments

    Parameters
    ----------
    argv : list
        Command line arguments.

    Returns
    -------
    out : argparse.Namespace
        Namespace of command line arguments.
    """
    parser = argparse.ArgumentParser(description="""RISM plotter.""")
    #Positional args
    parser.add_argument('file', metavar='file',
                        help="""Supported graph file.""")
    parser.add_argument('-c', '--columns',
                        help="""Specify number(s) of column in vv file to be
                        plotted. Numbering starts from 1.""", nargs='+',
                        type=int)
    parser.add_argument('--group',
                        help="""Group graphs on one figure.""",
                        action='store_true')
    parser.add_argument('--do_coord',
                        help=""" Find coordination number / density.""",
                        action='store_true')

    return parser.parse_args(argv)    
                

def grouper(iterable, n, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx
    args = [iter(iterable)] * n
    return izip_longest(fillvalue=fillvalue, *args)
                

def find_first_minima_idx(r, rdf):
    """Finds first local minima. The search is relatively unsofisticated.
    Plots the figure to confirm minima value. If the user disagrees with minima
    found asks user to provide estimate.
    """
    # Get first local minima idx and values
    grt_than_one_first_idx = np.where(rdf > 1.0)[0][0]
    cut_rdf = rdf[grt_than_one_first_idx:]
    loc_min_idx = argrelextrema(cut_rdf, np.less, order=3)[0][0]
    loc_min_value = cut_rdf[loc_min_idx]
    loc_min_r = r[np.where(rdf==loc_min_value)[0][0]]
    # Plot the found minimum and confirm it
    #return plot_and_pick(r, rdf, loc_min_r, loc_min_value)
    return loc_min_r, loc_min_value

                
def plot_gav(filename, group, do_coord, columns_to_plot=None):
    with open(filename, 'rb') as f:
        lines = f.readlines()
    # header
    roux=False
    if 'R/ANGSTROM' in lines[4]:
        # q-espresso
        strings = lines[4].split()[1:]
        columns = []
        for i1, i2, i3 in grouper(strings, 3):
            column = ' '.join([i1, i2, i3])
            columns.append(column)
        numbers = lines[5:]
    elif lines[0].startswith('#'):
        #amber rism
        columns = lines[3][1:].split()[1:]
        numbers = lines[4:]
        print(columns)
    elif lines[0].startswith('*'):
        #roux rism
        roux=True
        numbers = lines[5:]
        columns = [str(i+1) for i in range(len(numbers))]
        numbers = [map(float, i.split()) for i in numbers]
        numbers = np.array(numbers).T
    else:
        # assume it is in numpy format
        numbers = lines
        columns = [str(i+1) for i in range(len(numbers[0].split()) - 1)]
    # check what  we need to plot
    if columns_to_plot:
        columns_to_plot = map(lambda x: x - 1, columns_to_plot)
    else:
        columns_to_plot = [i+1 for i in range(len(columns))]
    print(columns_to_plot)
    # get data
    distance = []   # x-axis
    g_of_r = []     # y-axis
    labels = []
    for i, column in enumerate(columns):
        i = i + 1
        if i in columns_to_plot:
            distance.append([])
            g_of_r.append([])
            labels.append(column)
            if roux:
                i=i-1
            for l_num,l in enumerate(numbers):
                if not roux:
                    nums = l.split()
                    distance[-1].append(float(nums[0]))
                else:
                    distance[-1].append(l_num)
                try:
                    g_of_r[-1].append(float(nums[i]))
                except ValueError:
                    print i
                    print l[i]
    # plot
    fig, ax = plt.subplots()
    for i, (x, y, label) in enumerate(zip(distance, g_of_r, labels)):
        if do_coord:
            x, y = np.array(x), np.array(y)
            loc_min_x, loc_min_y = find_first_minima_idx(x, y)
            min_idx = np.argmin(abs(x - loc_min_x)) + 1
            dx = x[1] - x[0]
            coord_num = simps(y[:min_idx]*4*np.pi*x[:min_idx]**2, dx=dx)
            ax.plot(loc_min_x, loc_min_y, 'ro')
            print(label, round(coord_num, 3))
        #print(x)
        #print(y)
        #print(label)
        ax.plot(x, y, label=label)
        if not group:
            plt.title(label)
            # check if we need new plot
            if i + 1 < len(labels):
                fig, ax = plt.subplots()
    if group:
        plt.legend(loc='best')    
    plt.show()
    
                

def main(argv):
    args = process_command_line(argv)
    #print args.file
    #if args.file[-3:] == 'gvv' or 'cvv':
    #print '1D-RISM g(r) file'
    plot_gav(args.file, args.group, args.do_coord, args.columns)
#    else:
#        print('Unsupported file format.')
        
if __name__ == '__main__':
    main(sys.argv[1:])
