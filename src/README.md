# TROVA: TRansport Of water VApor

TRansport Of water VApor (TROVA) is a software developed in Python and Fortran for the study of moisture sources and sinks. It has been developed within the LAGRIMA project at the EPhysLab (Environmental Physics Laboratory) at the University of Vigo.

# What is TROVA?

TROVA allows the use of the FLEXible PARTicle global dispersion model and the FLEXPART-WRF regional model at different spatial resolutions. It also include the methodologies of Stohl and James (2005) and Sodemann et al. (2008). It contains two main modules:

1- Developed in Python that is responsible for reading the files, configuring TROVA and generating the outputs of the moisture budget *(Evaporation (E)-Precipitation (P))* for the number of days selected in the simulations.

2- Developed in Fortran that is used in interface with Python so that the calculations of great computational demand are carried out in the shortest possible time.


# What do I need to get and run TROVA?

## Prerequisites and Installation
This section describes the prerequisites required to run TROVA, as well as the steps to install it.

### Prerequisites

To run TROVA, you need
* Python3
* [![Git](https://git-scm.com/)](https://git-scm.com/)
* [![Anaconda3](https://www.anaconda.com/)](https://www.anaconda.com/)
* Linux

The main packages required to run TROVA are:
```
numpy 
mpi4py 
time
struct
datetime
netCDF4 
scipy
functools
pathlib 
gzip
shutil
```
