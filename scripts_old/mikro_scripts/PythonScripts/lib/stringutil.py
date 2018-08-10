# -*- coding: utf-8 -*-
"""
Created on Sat Feb 16 21:55:37 2013

@author: a92549
"""

import re


def case_insensitive_replace(string, old, new):
    case_insensitive_old = re.compile(re.escape(old), re.IGNORECASE)
    return case_insensitive_old.sub(new, string)
    
