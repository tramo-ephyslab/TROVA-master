#!/bin/bash

############################################################
# Install fortran functions to use in python
############################################################
f2py -c -m functions functions.f90 
