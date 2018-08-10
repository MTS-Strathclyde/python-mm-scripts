#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 16 18:25:38 2013

@author: max
"""

import sys
import subprocess
import argparse



input_dic = {'cpcm' : """ $system mwords=125 $end
 $contrl dfttyp=b3lyp runtyp=energy coord=cart units=angs scftyp=rhf $end
 $contrl icharg={0} $end
 $basis gbasis=N31 ngauss=6 ndfunc=1 $end
 $pcm smd=.TRUE. solvnt=water $end
 $guess  guess=huckel $end""",
 
 'gas' : """ $contrl dfttyp=b3lyp $end
 $basis GBASIS=N31 NGAUSS=6 NDFUNC=1 $end
 $contrl icharg={0} $end""",

 'smd' : """ $contrl coord=cart units=angs scftyp=rhf dfttyp=b3lyp runtyp=energy $end
 $contrl icharg={0} $end
 $system mwords=125 $end
 $pcm smd=.TRUE. solvnt=water ief=-3 $end
 $tescav mthall=2 ntsall=240 $end
 $guess  guess=huckel $end
 $basis gbasis=N31 ngauss=6 ndfunc=1 $end""",
 
 'm06_smd' : """ $contrl coord=cart units=angs scftyp=rhf dfttyp=m06-2x runtyp=energy $end
 $contrl icharg={0} $end
 $system mwords=125 $end
 $pcm smd=.TRUE. solvnt=water ief=-3 $end
 $tescav mthall=2 ntsall=240 $end
 $guess  guess=huckel $end
 $basis gbasis=N31 ngauss=6 ndfunc=1 $end"""
}

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
    parser = argparse.ArgumentParser(description=""" Prepare AMBER prmtop file.""")
    #Positional args
    parser.add_argument('file', metavar='molec.pdb',
                        help="""Input file""", nargs='+')
    #Optional args
    parser.add_argument('-t', '--calctype',
                        help=""" Gamess calc type. Can select more than one.
                        [gas, smd, cpcm, m06_smd]""", default=['smd'], nargs='+')
    parser.add_argument('-c', '--molcharge',
                        help="""Charge of the solute [0]""", default=0,
                        type=int)

    return parser.parse_args(argv)


def main(argv):
    args = process_command_line(argv)
    charge = args.molcharge
    for name in args.calctype:
        input_string = input_dic[name]
        with open(name, 'wb') as f:
            f.write(input_string.format(charge))
    for mol in args.file:
        for name in args.calctype:
            subprocess.call(['babel', '-xf', name, '-ipdb', mol, '-ogamin', mol[:-4] + '_' + name + '.inp'])


if __name__ == '__main__':
    main(sys.argv[1:])
