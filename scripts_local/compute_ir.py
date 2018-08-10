#!/usr/bin/env python2
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import argparse


def parse_argv(argv):
    parser = argparse.ArgumentParser(description="Compute IR from dipole correlation")
    parser.add_argument('file', help="Dipole.traj")
    parser.add_argument('-t', '--timestep', help="Dipole printing timestep in fs [0.5]",
                        default=0.5, type=float)
    parser.add_argument('--deriv', help="Derivative of dipoles is provided",
                        action='store_true')
    return parser.parse_args(argv)


def read_dipole_traj(fname):
    dipole_vecs = []
    with open(fname) as f:
        for l in f:
            if l.startswith('    X='):
                l = l.split()
                dipole_vecs.append([l[1], l[3], l[5]])
            #elif l.startswith(' DIPOLE [Non Periodic] DERIVATIVE(A.U.)|'):
            elif l.startswith(' DIPOLE [Non Periodic](Debye)|'):
                l = l.split()
                dipole_vecs.append([l[-3], l[-2], l[-1]])
    return np.array(dipole_vecs, dtype=float)


def compute_correlation(dipoles, max_cor_steps):
    N = len(dipoles)
    correlation = np.zeros(max_cor_steps)
    for i in range(N - max_cor_steps):
        dip = dipoles[i, :]
        if i%10000 == 0:
            print 'Computed correlation for first {} dipoles'.format(i)
        correlation += np.dot(dipoles[i:i+max_cor_steps, :], dip)
    # normalize correlation
    return correlation/(N - np.arange(max_cor_steps))*N/correlation[0]


def write_correlation(correlation, fname, plot=False):
    np.savetxt(fname, correlation)
    if plot:
        plt.clf()
        plt.plot(correlation[:,0], correlation[:,1], '-')
        plt.savefig(fname + '.png', dpi=300)


def compute_fourier(correlation, timestep, max_cor_steps, derivative):
    spectrum = []
    N = len(correlation)
    # max frequency = 4000 cm^-1, time is assumed to be in fs (hence 1.0e-15)
    max_freq = N*4000*timestep*1.0e-15*29979245800.0
    for i in range(int(round(max_freq))):
        omega = 2*np.pi*i/N
        omegas = 2*np.pi*i/N*np.arange(-max_cor_steps + 1, max_cor_steps)
        integral = np.sum(np.cos(omegas)*correlation)
        if not derivative:
            integral = integral*omega**2
        spectrum.append((omega/(2*np.pi*timestep*1.0e-15*29979245800.0), integral))
    return np.array(spectrum)


def main(argv):
    args = parse_argv(argv)
    dipoles = read_dipole_traj(args.file)
    N = len(dipoles)
    max_cor_steps = min(N/5, 10000)
    print 'Using at most {} steps'.format(max_cor_steps)
    correlation = compute_correlation(dipoles, max_cor_steps)
    correlation = np.r_[correlation[::-1][:-1], correlation]
    spectrum = compute_fourier(correlation, args.timestep, max_cor_steps, args.deriv)
    fname = os.path.splitext(args.file)[0]
    time_correlation = np.c_[np.arange(-max_cor_steps+1, max_cor_steps)*args.timestep, correlation]
    write_correlation(time_correlation, fname+'.time', True)
    write_correlation(spectrum, fname+'.freq', True)
    

if __name__ == '__main__':
    main(sys.argv[1:])

