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
* [Git] https://git-scm.com/
* [Anaconda3] https://www.anaconda.com/
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
 
 3- Enter the src/ directory and execute the install.sh code. This creates the executable of the fortran functions used by python for numerical computation.
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
 model = FLEXPART-WRF              Model type (e.g. FLEXPART or FLEXPART-WRF)
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

 
# References
[1] Stohl A, James PA. A Lagrangian analysis of the atmospheric branch of the global water cycle: Part II: Earth’s river catchments ocean basins, and moisture transports between them. J. Hydrometeorol. 2005; 6:961–984. https://doi.org/10.1175/JHM470.1.

[2] Sodemann H, Schwierz C, Wernli H. Interannual variability of Greenland winter precipitation sources: Lagrangian moisture diagnostic and North Atlantic Oscillation influence. J. Geophys. Res.-Atmos. 2008; 113:D03107. https://doi.org/10.1029/2007JD008503. 
