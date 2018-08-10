# -*- coding: utf-8 -*-

"""Module with utility functions for solvation DB"""

import os
import pybel
import signal

from cinfony import rdk
from cinfony import webel
from cinfony import opsin

def web_service_down_handler(signum, frame):
    raise IOError("web server didn't respond after 3 sec.")


def get_best_name(names):
    for name in names:
        name_no_signs = name.replace('-', '').replace(' ', '')
        if not name_no_signs.isdigit():
            return name
    else:
        return names[0]


def get_IUPAC_name_from_web(inchi):
    """Will try getting name from web service."""
    webmol = webel.readstring('inchi', inchi)
    try:
        name = webmol.write('iupac').capitalize()
        if not name:
            names = webmol.write('names')
            if len(names) == 1:
                name = webmol.write('names')[0].capitalize()
            else:
                name = get_best_name(names)
    except IndexError:
        name = ''
    except AttributeError:
        name = ''
    return name


def get_IUPAC_name(inchi):
    """
    Will try to find molecule IUPAC name on the internet.
    In case of fail will return empty string.
    >>> print get_IUPAC_name('CCC')
    PROPANE
    If molecule doesn't have iupac name will try to find any name by using
    names. If even that fails will return an empty string.
    """
    signal.signal(signal.SIGALRM, web_service_down_handler)
    signal.alarm(1)   # number of seconds before function is killed
    try:
        name = get_IUPAC_name_from_web(inchi)
    except IOError:
        name = ''
    signal.alarm(0)
    return name
#    return ''    
   
def get_RMSD_value(refmol, probemol):
    """Input is 2 mol files."""
    rdref = rdk.readstring('mol', str(refmol))
    rdprobe = rdk.readstring('mol', str(probemol))
    return rdk.Chem.AllChem.GetBestRMS(rdref.Mol, rdprobe.Mol)


#def PK_inchi_from_pymol(pymol):
#    """Converts pybel molecule into db primary key."""
#    xyz = pymol.write('xyz')
#    xyz_pymol = pybel.readstring('xyz', xyz)
#    return xyz_pymol.write('inchi').strip()


def pymol_2_inchi(pymol):
    return pymol.write('inchi').strip()


def std_inchi(inchi):
    """Converts InChI=1/... to InChI=1S/."""
    head, tail = inchi.split('/', 1)
    return 'InChI=1S/' + tail


def mol_2_inchi(mol):
    pymol = pybel.readstring('mol', mol)
    return pymol.write('inchi').strip()


def inchi_2_smi(inchi):
    """Converts inchi to smiles."""
    pymol = pybel.readstring('inchi', inchi)
    return pymol.write('smi').strip()


def smi_2_inchi(smi):
    """Converts smiles to inchi."""
    pymol = pybel.readstring('smi', smi)
    return pymol.write('inchi').strip()
        

def mol_converter(txt, informat, outformat=None):
    """Converts molecule using pybel."""
    pymol = pybel.readstring(informat, txt)
    if outformat:
        return pymol.write(outformat)
    else:
        return pymol    


def silent_mol_converter(txt, informat, outformat=None):
    """Converts molecule suppressing all pybel output.
    Returns pymol if outformat is None."""
    #open 2 fds
    null_fds = [os.open(os.devnull, os.O_RDWR) for x in xrange(2)]
    # save the current file descriptors to a tuple
    save = os.dup(1), os.dup(2)
    # put /dev/null fds on 1 and 2
    os.dup2(null_fds[0], 1)
    os.dup2(null_fds[1], 2)
    # *** run the function ***
    pymol = pybel.readstring(informat, txt)
    # restore file descriptors so I can print the results
    os.dup2(save[0], 1)
    os.dup2(save[1], 2)
    # close the temporary fds
    os.close(null_fds[0])
    os.close(null_fds[1])
    if outformat:
        return pymol.write(outformat)
    else:
        return pymol
