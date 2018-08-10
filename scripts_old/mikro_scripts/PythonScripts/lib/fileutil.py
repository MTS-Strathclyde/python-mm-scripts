# -*- coding: utf-8 -*-
"""
Created on Mon Dec  3 18:33:09 2012

@author: a92549

This module can be used to store typical file system operations, which
are used in scrypts.

"""

import os

def list_ext(ext, folder='.'):
    """Returns absolute pathes of files, which have given extension.
    By default looks in working directory.
    Returns them as list.
    """
    if folder[-1] == '/':
        folder = folder[:-1]
    folder_path = os.path.join(os.getcwd(), folder)
    all_files = os.listdir(folder_path)
    ext_list = []
    for filename in all_files:
        if os.path.splitext(filename)[1] == ext:
            ext_list.append(os.path.join(folder_path, filename))
    return ext_list            
            
def change_ext(path, new_ext):
    """Accepts filenames as well. Extension should be supplied with dot.
    Example:
        >>> change_extension('f.ext', '.o')
        f.o
    """
    return os.path.splitext(path)[0] + new_ext
    
def get_filename(path):
    """Retuns filename, remove path and extension."""
    with_ext = os.path.split(path)[1]
    return os.path.splitext(with_ext)[0]


def remove_path(path):
    """Returns filename without path"""
    return os.path.split(path)[1]


def remove_extenison(path):
    """Returns path without extension"""
    return os.path.splitext(path)[0]


def add_to_name_but_keep_ext(path, str_to_add):
    """Returns path with appended string between end of old file name and 
    extension."""
    name, ext = os.path.splitext(path)
    return name + str_to_add + ext
    
def file_exist(path):
    try:
        with open(path) as f:
            return True
    except IOError:
       return False
