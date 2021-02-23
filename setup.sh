#!/bin/bash

# Compiling Cython module: cyadd
python3 setup.py build_ext --inplace

mv cyadd*.so corefuns/cyadd.so
