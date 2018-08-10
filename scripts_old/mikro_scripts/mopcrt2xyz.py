#!/usr/bin/python

"""
Converts mopcrt to xyz, without specifing partial charges
"""

import sys

def mopcrt_to_xyz(mop_string, mol_name):
    mop_coord_string = mop_string.split('\n', 3)[3]
    mop_coords  = mop_coord_string.split('\n')
 #   print mop_coords
    num_of_atoms = 0
    xyz_coords = "{0}\n{1}\n"
    for line in mop_coords:
        strs = line.split()
#        print line
        if len(strs) > 0:
            num_of_atoms += 1
            xyz_coords += strs[0] + ' ' + strs[1] + ' ' + strs[3] + ' ' + strs[5] + '\n'
    return xyz_coords.format(num_of_atoms, mol_name)


def main(args):
    for mol in args:
        #conversion
        mop_file = file(mol, 'r')
        xyz_coords = mopcrt_to_xyz(mop_file.read(), mol)
        #writing
        out_name = mol[:-4] + ".xyz"
        out_file = file(out_name, 'w')
        out_file.write(xyz_coords)
    print 'Job finished'

if __name__ == '__main__':
    main(sys.argv[1:])
