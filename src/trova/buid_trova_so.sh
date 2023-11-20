#!/bin/bash

############################################################
# Install fortran functions to use in python
############################################################

FILE=`ls -a *.so` 2>/dev/null
echo $FILE

if [ -f "$FILE" ]
then
    echo "There are compiled fortran functions"
    rm *.so
    f2py -c --f90flags=-Wno-tabs -m functions functions.f90
else
    echo "Compiled fortran functions do not exist"
    echo "------------------------------------"
    echo "!!!!!!!Compiling functions...."
    echo "------------------------------------"
    f2py -c --f90flags=-Wno-tabs -m functions functions.f90
fi
