#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 22 14:03:03 2016

@author: max
"""

import re
import sys
import fileinput

txt = ''
for line in fileinput.input():
    txt += line

# regexp from http://stackoverflow.com/a/385597

# modified so it ignores integers in the form a
floats = re.findall(r"[+-]? *(?:\d+(?:\.\d*)|\.\d+)(?:[eE][+-]?\d+)?", txt)

# modified so it only looks for the form a.bEc
#floats  = re.findall(r"[+-]? *(?:\d+(?:\.\d*)|\.\d+)(?:[eE][+-]?\d+)", txt)


for f in floats:
    print float(f)
    
