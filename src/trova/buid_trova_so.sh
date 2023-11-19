#!/bin/bash

############################################################
# Install fortran functions to use in python
############################################################

FILE=`ls -a *.so`
echo $FILE

if [ -f $FILE ]
then
    rm *.so
    f2py -c --f90flags=-Wno-tabs -m functions functions.f90
else
    f2py -c --f90flags=-Wno-tabs -m functions functions.f90
fi
