#!/usr/bin/python

"""
Converts mopout file to mopcrt, without specifing partial charges
"""

import sys
import pybel

def xyz_to_mop(xyz_string):
    xyz_coord_string = xyz_string.split('\n', 2)[2]
    xyz_coords  = xyz_coord_string.split('\n')
#    print xyz_coords
    mop_coords = ""
    for line in xyz_coords:
        strs = line.split()
        if len(strs) == 4:
            mop_coords += strs[0] + ' ' + strs[1] + ' +1 ' + strs[2] + ' +1 ' +\
                          strs[3] + ' +1\n'
    return mop_coords


def main(args):
    for mol in args:
        #conversion
        pymol = pybel.readfile("mopout", mol).next()
        mop_coords = xyz_to_mop(pymol.write("xyz"))
        #writing
        out_name = mol[:-4] + ".mop"
        out_file = file(out_name, 'w')
        out_file.write('\n ' + mol + '\n\n' + mop_coords)
    print 'Job finished'

if __name__ == '__main__':
    main(sys.argv[1:])
