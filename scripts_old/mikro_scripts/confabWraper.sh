#!/bin/bash

#Setups enviromental variables for confab and runs it

export LD_LIBRARY_PATH=~/Tools/confab-install/lib
export BABEL_DATADIR=~/Tools/confab-install/share/openbabel/2.2.99
~/Tools/confab-install/bin/confab $*