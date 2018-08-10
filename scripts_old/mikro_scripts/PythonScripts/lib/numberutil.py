# -*- coding: utf-8 -*-
"""
Created on Thu Feb 14 12:33:50 2013

@author: a92549
"""

def is_float(string):
    try:
        float(string)
        return True
    except ValueError:
        return False
        