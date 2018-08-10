# -*- coding: utf-8 -*-
"""
Created on Wed Sep 12 13:32:22 2012

@author: mishin1991

Division functions

"""
        
def sep_by_en(pybel_mols, incr_bound=0.1):
    """
    Produces list of groups, which contain filenames of molecules divided by energies    
    """
    #prepare vars
    groups = []
    pr_en = -float("inf")
    #Main cycle
    for mol in pybel_mols:
        ener = mol.energy
        if ener - pr_en > incr_bound:
            groups.append([])
            groups[-1].append(mol)
        else:
            groups[-1].append(mol)
        pr_en = ener
    #remove groups with one molecule
    corrected_groups = [groups[0]]
    for i in range(1, len(groups)):
        if len(groups[i]) < 2:
            corrected_groups[-1] += groups[i]
        else:
            corrected_groups.append(groups[i])
    #Info 
    print "Total number of groups is", len(corrected_groups)
    return corrected_groups

def main(func_name, pybel_mols):
    if func_name == 'sep_by_en':
        return sep_by_en(pybel_mols)
    else:
        raise Exception("unknown division method")