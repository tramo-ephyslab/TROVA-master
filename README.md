# TROVA: TRansport Of water VApor

TRansport Of water VApor (TROVA) is a software developed in Python and Fortran for the study of moisture sources and sinks. It has been developed within the LAGRIMA project at the EPhysLab (Environmental Physics Laboratory) at the University of Vigo.

# What is TROVA?

TROVA allows the use of the FLEXible PARTicle global dispersion model and the FLEXPART-WRF regional model at different spatial resolutions. They also include the methodologies of Sthol and James (2005) and Sodemann et al. (2008). 

It contains 2 main modules:

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
### Installation
 1- Clone the repository 
 ```
 git clone https://github.com/tramo-ephyslab/TROVA-master.git
 ```
 2- Install python requirements in anaconda environment.
 
 3- Enter the src/ directory and execute the *install.sh* code. This creates the executable of the fortran functions used by python for numerical computation.
 ```
 sh install.sh
 ```
 
### Modification of the input file to run TROVA

A description of each parameter is shown below:
 ```
 path                              Input data path
 path_output                       Outputs path
 mode                              Mode for particle tracking [1 forward, -1 backward]
 mass                              Atmospheric mass of particles
 numP                              Particles number
 type_file                         Parameter to use data from FLEXPART-global (2) or FLEXPART-WRF (1)
 resolution                        Output data resolution 
 numPdX                            Points number for output mesh in direction X
 numPdY                            Points number for output mesh in direction Y
 x_lower_left                      Bottom left corner longitude
 y_lower_left                      Bottom left corner latitude
 step_numbers                      Number of time steps to analyze (for each day are 4 steps)
 year                              Initialization year for particle tracking  
 month                             Initialization month for particle tracking 
 days                              Initialization day for particle tracking   
 hour                              Initialization hour for particle tracking 
 ndays                             Number of days to perform particle tracking
 file_mask                         Mask path 
 name_mascara                      Mask variable name in the mask file
 name_variable_lat                 Latitude variable name in the mask file
 name_variable_lon                 Longitude variable name in the mask file
 x_left_lower_corner               lower left longitude of input data
 y_left_lower_corner               lower left latitude of input data
 x_rigth_upper_corner              upper rigth longitude of input data
 y_rigth_upper_corner              upper rigth latitude of input data
 model                             Model type (e.g. FLEXPART or FLEXPART-WRF)
 type_lon                          Format longitude [1 ([0-360]), 2 (lon [-180-180])
 method                            Method to use [1 (Sthol), 2 (Sodemann)]
 threshold                         Threshold for filtering precipitating particles for both methods [e.g. -0.00025 kg/kg]
 filter_value                      Parameter to activate filter of precipitating particles [0 no activate, 1 activate]
 output_txt                        Save outputs in .txt for each day [0 no, 1 yes]
 output_npy                        Save outputs in .npy for each day [0 no, 1 yes] 
 output_nc                         Save outputs in .nc for each day [0 no, 1 yes] 
 name_target_region                Target region name
 value_mask                        Mask value in file mask
 file_gz                           Parameter to know whether the input data is compressed or not [0 no, 1 yes]
 save_position_part                Save the dq/dt values for each position of each particle [0 no, 1 yes]
```
### Input data

In addition, to run TROVA the following data sets are needed:
- Outputs of traces air parcels and their properties from the global dispersion model FLEXPARTv9 (Piso et al., 2019) or higher or from its version for regional domains FLEXPART-WRFv3.3.2 (Brioude et al., 2013).

- Target regions mask that will be used for tracking the particles.

# How do I run TROVA?

To run TROVA, change into the src directory

**num_CPU:** *CPU numbers to use (preferably divisible by 4).*

**input_file_path:** *Input file path for TROVA run.*

On a Linux computer:

```
cd src
mpirun -np num_CPU python TROVA.py input_file_path
```

On a HPC with Linux:

```
cd src
Create an execution code for example for a queue manager like slurm (See example/run_example.sh)
Then execute: sbatch run_example.sh
```

# Examples
Here, we provide two examples:

1- TROVA is used to determine moisture sources for a tropical cyclone (October 17, 2014, at 18 UTC). In this case, a regional mask is used and the outputs of the FLEXPART dispersion model are used (Masks/mask_AL082014_20141017_18.nc). The configuration file for this case is shown in the Inputs/input_back_TC.cfg directory. The result is displayed in the directory: Figures/E_P_10-day_TC_backward.png.

2- TROVA is used in forward in time to determine the moisture sinks (October 10, 2014, at 00 UTC) associated with the main source of the North Atlantic Ocean (NATL). In this case, a global mask is used and the outputs of the FLEXPART dispersion model are used (Masks/NATL.nc). The configuration file for this case is shown in the directory Inputs/input_forw_NATL.cfg. 

The necessary data to be able to carry out tests with TROVA can be downloaded at the link: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6490365.svg)](https://doi.org/10.5281/zenodo.6490365)

# Important notes

This code is not bug-free. Please report any bugs through 'Issues': https://github.com/tramo-ephyslab/TROVA-master/issues

## Contact and support

José Carlos fernández Alvarez (jose.carlos.fernandez.alvarez@uvigo.es) and Albenis Pérez Alarcón (albenis.perez.alarcon@uvigo.es)
 
# References
[1] Stohl A, James PA. A Lagrangian analysis of the atmospheric branch of the global water cycle: Part II: Earth’s river catchments ocean basins, and moisture transports between them. J. Hydrometeorol. 2005; 6:961–984. https://doi.org/10.1175/JHM470.1.

[2] Sodemann H, Schwierz C, Wernli H. Interannual variability of Greenland winter precipitation sources: Lagrangian moisture diagnostic and North Atlantic Oscillation influence. J. Geophys. Res.-Atmos. 2008; 113:D03107. https://doi.org/10.1029/2007JD008503. 

[3] Piso I et al. The Lagrangian particle dispersion model FLEXPART version 10.3. Geosci. Model Dev. Discuss. 2019; https://doi.org/10.5194/gmd-2018-333.

[4] Brioude J et al. The Lagrangian particle dispersion model FLEXPART-WRF version 3.1. Geosci. Model Dev. 2013; 6:1889–1904. https://doi.org/10.5194/gmd-6-1889-2013.
