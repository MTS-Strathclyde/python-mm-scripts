#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 30 17:33:09 2015

@author: max
"""

import sys
import argparse
import numpy as np


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
    parser = argparse.ArgumentParser(description="""Smooth multiple rdfs.""")
    #Positional args
    parser.add_argument('rdfs',
                        help="""RDF files readable by numpy""", nargs='+')
    #Optional args
    parser.add_argument('-w', '--window_size',
                        help=""" Window in savitzky golay filter. Extent of
                        smoothing. [51]""",
                        default=51, type=int)
    parser.add_argument('-o', '--order',
                        help=""" Order of smoothing polynomial [3].""",
                        default=3, type=int)
    parser.add_argument('-p', '--plot',
                        help=""" Plot results.""", action='store_true')
    parser.add_argument('-n', '--name',
                        help=""" Name of produced rdf [rdf.mat]""",
                        default='rdf.mat')
    parser.add_argument('--to_angs',
                        help=""" Write output coordinates in Angstroms""",
                        action='store_true')
    return parser.parse_args(argv)
    
    
def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    r"""Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
    The Savitzky-Golay filter removes high frequency noise from data.
    It has the advantage of preserving the original shape and
    features of the signal better than other types of filtering
    approaches, such as moving averages techniques.
    Parameters
    ----------
    y : array_like, shape (N,)
        the values of the time history of the signal.
    window_size : int
        the length of the window. Must be an odd integer number.
    order : int
        the order of the polynomial used in the filtering.
        Must be less then `window_size` - 1.
    deriv: int
        the order of the derivative to compute (default = 0 means only smoothing)
    Returns
    -------
    ys : ndarray, shape (N)
        the smoothed signal (or it's n-th derivative).
    Notes
    -----
    The Savitzky-Golay is a type of low-pass filter, particularly
    suited for smoothing noisy data. The main idea behind this
    approach is to make for each point a least-square fit with a
    polynomial of high order over a odd-sized window centered at
    the point.
    Examples
    --------
    t = np.linspace(-4, 4, 500)
    y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
    ysg = savitzky_golay(y, window_size=31, order=4)
    import matplotlib.pyplot as plt
    plt.plot(t, y, label='Noisy signal')
    plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
    plt.plot(t, ysg, 'r', label='Filtered signal')
    plt.legend()
    plt.show()
    References
    ----------
    .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
       Data by Simplified Least Squares Procedures. Analytical
       Chemistry, 1964, 36 (8), pp 1627-1639.
    .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       Cambridge University Press ISBN-13: 9780521880688
    """
    from math import factorial

    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError, msg:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')


def main(argv):
    args = process_command_line(argv)
    mats = [np.loadtxt(f) for f in args.rdfs]
    all_mat = np.vstack(mats)
    all_mat = all_mat[all_mat[:,0].argsort()]
    smooth = savitzky_golay(all_mat[:,1], args.window_size, args.order)
    # rdf < 0 is unphysical
    smooth[smooth<0] = 0 
    if args.to_angs:
        r = all_mat[:,0]*10.
    else:
        r = all_mat[:,0]
    smooth = np.c_[r, smooth]
    np.savetxt(args.name, smooth)    
    if args.plot:
        import matplotlib.pyplot as plt
        plt.plot(smooth[:,0], smooth[:, 1])
        plt.show()
        
if  __name__ == '__main__':
    main(sys.argv[1:])
    
    