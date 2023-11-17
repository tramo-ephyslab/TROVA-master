#!/bin/bash

############################################################
# Install fortran functions to use in python
############################################################
rm *.so
f2py -c -m functions functions.f90 
