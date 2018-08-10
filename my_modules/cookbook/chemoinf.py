# -*- coding: utf-8 -*-
import pybel

def xyz_to_rism(xyz):
    """Convert xyz file to rism pdb.
    Doesn't use any 3rd party modules."""
    xyz_lines = xyz.splitlines()
    rism_lines = []
    line_format = '{0[0]}{0[1]:>7}{0[2]:>4}{0[3]:>5}{0[4]:>6}{0[5]: 12.3f}{0[6]: 8.3f}{0[7]: 8.3f}{0[8]:>6}{0[9]:>6}'
    atom_count = 1
    element_counts = {}
    for line in xyz_lines[2:]:
        rism_line = ['ATOM']
        rism_line.append(atom_count)
        atom_count += 1
        symbol, x, y, z = line.split()
        element_counts[symbol] = element_counts.get(symbol, 0) + 1 # add element to frequency dictionary
        rism_line.append(symbol + str(element_counts[symbol]))
        rism_line.append('MOL')
        rism_line.append(1)        
        rism_line.append(float(x))
        rism_line.append(float(y))
        rism_line.append(float(z))
        rism_line.extend(['1.00', '0.00'])        
        new_line = line_format.format(rism_line)
        rism_lines.append(new_line)
    rism_lines.extend(['TER', 'END'])
    return '\n'.join(rism_lines)


def mol_converter(txt, informat, outformat=None):
    """Converts molecule using pybel."""
    pymol = pybel.readstring(informat, txt)
    if outformat:
        return pymol.write(outformat)
    else:
        return pymol    


