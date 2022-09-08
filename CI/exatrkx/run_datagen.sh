#!/bin/bash
source $HOME/Documents/acts_project/acts/build/python/setup.sh

$PREFIX python3 $HOME/Documents/acts_project/acts/CI/exatrkx/datagen.py $@
#rm -f timing.tsv
