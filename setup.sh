#!/bin/bash

# Compiling Cython module: cyadd
python3 setup.py build_ext --inplace

mv cyadd*.so corefuns/cyadd.so

# seting up permissions
#chmod u+x preprocessgwas.sh
chmod u+x scripts/*.sh
chmod u+x scripts/plink
chmod u+x scripts/cassi-run.sh
chmod u+x cassi/cassi

export CURRENTDIR=`pwd`
#export PYTHONPATH=$CURRENTDIR/scripts
export PYTHONPATH=$CURRENTDIR:$CURRENTDIR/scripts:$CURRENTDIR/corefuns:$CURRENTDIR/datatools:$CURRENTDIR/classes
export PATH=$CURRENTDIR/scripts/:$PATH
export PATH=$CURRENTDIR/:$PATH
export PATH=$CURRENTDIR/cassi:$PATH
alias plink=$CURRENTDIR/scripts/plink

