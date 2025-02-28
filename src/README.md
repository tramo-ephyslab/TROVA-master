# TROVA: TRansport Of water VApor

TTRansport Of water VApor (TROVAv1.1.1) is a software developed in Python and Fortran
for the study of moisture sources and sinks. It has been developed within the LAGRIMA and 
SETESTRELO projects at the EPhysLab (Environmental Physics Laboratory) at the University of Vigo. 
Subsequently, its development and updating have continued within a collaboration between the University 
of Vigo and the Galician Supercomputing Center. Many investigations use this software to obtain scientific results. 
These can be consulted at the following web address: https://ephyslab.uvigo.es/en/staff/.
**This is an update of the software presented by Fernández-Alvarez et al. (2022)**

```
#*****************************************************************************************
#*                    EPhysLab (Environmental Physics Laboratory), Spain                 *
#*                        Galician Supercomputing Center, Spain                          *
#*                        TRansport Of water VApor (TROVA)                               *
#*                             version 1.1.1 (15-02-2025)                                *
#*                        _____ __    ____                                               *
#*                          |  |  |  /    \ \        //\                                 *
#*                          |  |__| /      \ \      //__\                                *
#*                          |  |  \ \      /  \    //    \                               *
#*                          |  |   \ \____/    \__//      \                              *
#*                                                                                       *
#*                       Edificio Campus da Auga/Edificio CESGA                          *
#*                            University of Vigo/CESGA                                   *
#*                          www.ephyslab.uvigo.es/www.cesga.es                           *
#*      contact: jose.carlos.fernandez.alvarez@uvigo.es (jcfernandez@cesga.es),          * 
#*                         albenis.perez.alarcon@uvigo.es                                *
#*****************************************************************************************
```

# What is TROVA?

TROVA allows the use of the FLEXible PARTicle global dispersion model and the FLEXPART-WRF 
regional model at different spatial resolutions. It also includes the methodologies 
of Stohl and James (2005) and Sodemann et al. (2008). We herein refer to these methodologies
as STHOL2005 and SOD2008 respectively. It contains two main modules:

1- Developed in Python, responsible for reading files, configuring TROVA, and generating 
the outputs of the moisture balance (Evaporation (E)-Precipitation (P)) for the number 
of days selected in the simulations.

2- Developed in Fortran, used in interface with Python to perform computationally 
demanding calculations in the shortest possible time. It also includes a parallel 
implementation using the MPI library to reduce TROVA's processing time.

3- This new version includes the analysis of moisture sources and sinks by vertical layers.

4- This version allows the calculation of the residence time of water vapor in the 
atmosphere for particles in a target region applying the methodology of Läderach and Sodemann (2016).

5- This version has functions that allow the representation of moisture source and sink patterns and the representation in a 2D graph of the residence time values of water vapor in the atmosphere for particles in a target region.

For a more detailed understanding of TROVA you can review the official TROVA API documentation: https://trova-docs.readthedocs.io/en/latest/.


# What do I need to get and run TROVA?

## Prerequisites and Installation
This section describes the prerequisites required to run TROVA, as well as the steps to install it.

### Prerequisites

To run TROVA, you need

- [![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/release/python-380/)

- [![Git](https://img.shields.io/badge/Git-blue.svg)](https://git-scm.com/)

- [![Anaconda 3](https://img.shields.io/badge/Anaconda-3-green.svg)](https://www.anaconda.com/)

- [![Linux](https://img.shields.io/badge/Linux-red.svg)](https://www.linux.org/)

- [![Fortran](https://img.shields.io/badge/Fortran-yellow.svg)](https://fortran-lang.org/)

1- Create a python environment with conda, for example

 ```
conda create -n py38 python=3.8
```

2- The main python packages that must be installed are the following (consider using the proposed options):

- numpy (*conda install numpy*)
- mpi4py (*pip install mpi4py* (or *conda install mpi4py*))
- time (*conda install -c conda-forge time*)
- netCDF4 (*conda install -c conda-forge netcdf4*)
- scipy (*conda install scipy*)
- importlib (*conda install -c conda-forge importlib*)
- cartoy (*conda install -c conda-forge cartopy*)
- setuptools (*pip install setuptools==58.2.0*)

### Installation
 1- Clone the repository 
 ```
 git clone https://github.com/tramo-ephyslab/TROVA.git
 ```
 2- Enter the TROVA/src/ directory and run the *install_trova.sh* code.
 ```
 sh install_trova.sh
 ```
 **_NOTE:_** With this option it is installed in the python environment and can be used directly.
 
## Input file and data to run TROVA

### Input file

To use TROVA you must modified the input file depending on the problem to be solved. This is an example that can be used as a test with the data presented below. A description of each parameter is shown below:
 ```
 #*****************************************************************************************
#*                    EPhysLab (Environmental Physics Laboratory), Spain                 *
#*                        Galician Supercomputing Center, Spain                          *
#*                        TRansport Of water VApor (TROVA)                               *
#*                             version 1.1.1 (15-02-2025)                                *
#*                        _____ __    ____                                               *
#*                          |  |  |  /    \ \        //\                                 *
#*                          |  |__| /      \ \      //__\                                *
#*                          |  |  \ \      /  \    //    \                               *
#*                          |  |   \ \____/    \__//      \                              *
#*                                                                                       *
#*                       Edificio Campus da Auga/Edificio CESGA                          *
#*                            University of Vigo/CESGA                                   *
#*                          www.ephyslab.uvigo.es/www.cesga.es                           *
#*      contact: jose.carlos.fernandez.alvarez@uvigo.es (jcfernandez@cesga.es),          * 
#*                         albenis.perez.alarcon@uvigo.es                                *
#*****************************************************************************************
#------------------------------------------------------------------------------------------
#Path to FLEXPART or FLEXPART-WRF partposit binary files [str]
path_data = "/home/jose/Documentos/TROVArun/Data/"

#Path for TROVA outputs [str]
path_output = "output/"

#Lagrangian tracking mode: ('backward' / 'forward'/ 'wvrt' / 'partposit') [str]
#backward: moisture sources, forward: moisture sinks, wvrt: water vapor residence time, partposit: particle variables over target region
mode = "backward"

#Atmospheric mass [float]
mass = 1.165725e+18

#Total number of atmospheric parcels in model simulation [int]
numP = 2045128

#Type of file: Set 1 for FELXPART-WRF and FLEXPART newler than version 9. Set 2 for FLEXPART older than version 9.  [int]
type_file = 1

#Spatial resolution for TROVA outputs [float]
resolution = 0.25 

#Number of point in x-direction for TROVA outputs [int]
numPdX = 600

#Number of point in y-direction for TROVA outputs [int]
numPdY = 325

#Lower longitude for TROVA output domain [float]
x_lower_left = -110

#Lower latitude for TROVA output domain [float]
y_lower_left = -15

#Time step for parcel tracking (minutes) [int]
dtime = 360

#Total time for parcel tracking (minutes) [int]
totaltime = 1440

#Start date for tracking [int]
year = 2014
month = 10
day = 17
hour = 00
min = 00

#Number of days to perform parcel tracking from start day [int]
ndays = 1

#path to mask fil (netcdf)
file_mask = "Masks/CAN.nc"

#Mask name variable in the mask filee [str]
maskname = "mask"     

#Latitude variable name  in the mask file [str]
maskvar_lat = "lat"

#Longitude variable name in the mask file [str]
maskvar_lon = "lon"

#Mask value for filterirng parcels [int]
mask_value = 1

#Subdomain limits for regional models [float]
#x_left_lower_corner: longitude min, y_left_lower_corner: latitude min, x_right_upper_corner: longitude max, y_right_upper_corner: latitude max
x_left_lower_corner = -100.0
y_left_lower_corner = -15.0
x_right_upper_corner = 39.86
y_right_upper_corner = 57.0

#model type: ['FLEXPART' / 'FLEXPART-WRF'] [str]
model = "FLEXPART-WRF"

#Set method = 1 for Stohl and James (2004,2005). Set method = 2 for Sodemann et al. (2008) [int]
method = 1

#To filter precipitating parcels ["True" / "False"]  [str]
filter_parcels_dqdt = False

#Threshold for filtering precipitating parcels [float]. It is only necessary if filter_parcels_dqdt = True.
dqdt_threshold = -0.0001

#To filter parcels by heigh ["True" / "False"]  [str]
filter_parcels_height = False

#Vertical layer for filtering parcels by height [lower_layer, upper_layer] [meters]. It is only necessary if filter_parcels_height = True.
filter_vertical_layers = [0,25000]

#To compute the moisture uptake in vertical layers ["True" / "False"]  [str]
use_vertical_layers = False

#Vertical layers to compute moisture uptake
vertical_layers = [0, 750, 900, 1500, 2250, 3000, 4000, 6000, 9000, 12000, 15000, 20000]

#File output format. Set 1 to activate output format and 0 to deactivate [int]
output_txt = 0
output_npy = 0
output_nc = 1

#Target region name [str]
name_target_region = "CAN"

#Set file_gz=1 if partposit files are compressed in gz format, else file_gz=0 [int]
file_gz = 0

#To save particle positions for each time step [str]
save_position_part = True

#To save dqdt positions for each dt [str]
save_position_dqdt = True

#Plotting identified parcels within the target region at time t0 (year_month_day_hour_min) [True /  False] [str]
plotting_parcels_t0 = False

#Ploting identified parcels trajectories on a map [True /  False] [str]
plotting_parcels_tracks_on_map = False

#Map limits for plotting [latmin, lonmin, latmax, lonmax, mapcenter, dlat, dlon] [float]
#map center must be 0 or 180. If center=180, provide lonmin and lonmax in 0-360 format
maps_limits = [0, -110, 75, 15, 0, 5, 25]

#Plotting 3D parcels trajectories [True /  False]
plotting_3Dparcels_tracks = False

#Calendar leap/noleap ["True" / "False"]  [str]
noleap = False

#Parameter to limit the particles to the domain limits. Consider only in regional models ["True" / "False"]  [str]
limit_domain = True

#Parameter to activate the calculation of the water vapor residence time together with the calculation of the humidity sources in backward mode ["True" / "False"]  [str]
method_wvrt = False

#Plotting the pattern of moisture sources or sinks [True /  False] [str]
plotting_moisture_sink_source = True

#Color pallete limits for plotting [min, max, step] [float]
limits_plot = [0, 3.2, 0.2]
```
### Input data

In addition, to run TROVA the following data sets are needed:
- Outputs of traces air parcels and their properties from the global dispersion model FLEXPARTv9 (Piso et al., 2019) or higher or from its version for regional domains FLEXPART-WRFv3.3.2 (Brioude et al., 2013).

- Target regions mask that will be used for tracking the particles.


# How do I run TROVA?

To run TROVA, create a file with the following code (could be **run_TROVA.py**)

```
#!/usr/bin/env python3
import trova as tv
import sys

input_file=sys.argv[1]
tv.TROVA_main(input_file)
```
Once this code is created, TROVA can be used in the following way:

On a Linux computer:

```
mpirun -np num_CPU python TROVA.py input_file_path

e.g (mpirun -np 4 python TROVA.py input_back_TC.cfg)
```

On a HPC with Linux (See https://github.com/tramo-ephyslab/TROVA-master/blob/main/run_example_HPC/run_example.sh):

```
#!/bin/bash -l

#SBATCH --mem=120GB
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -t 00:10:00

module --purge
module load cesga/2020
module load miniconda3/4.9.2
conda activate envname

srun -n $SLURM_NTASKS --mpi=pmi2 python  run_TROVA.py input_file_path >> py_${SLURM_JOB_ID}.log

```
**num_CPU:** *CPU numbers to use (preferably divisible by 4).*

**input_file_path:** *Input file path for TROVA run.*

## Illustrative examples for a specific day

Using the data corresponding to the outputs of FLEXPART-WRF forced with ERA5  provided in the following repositories [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14886887.svg)](https://zenodo.org/records/14886887)
and [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14939160.svg)]( https://zenodo.org/records/14939160), the results shown below for October 17, 2014 (moisture source, using as CAN target region), and October 6, 2014 (moisture sink, using as MED target region), can be reproduced. Configuration files and masks are available at https://github.com/tramo-ephyslab/TROVA-master/tree/main/Masks and https://github.com/tramo-ephyslab/TROVA-master/tree/main/Inputs. In these cases, the methodology of Stohl and James (2005) is used. Check that once TROVA is installed you can reproduce these results.

**1- Moisture source pattern (17/10/2014)** 

<p align="center" width="100%">
 <img src="https://github.com/tramo-ephyslab/TROVA-master/blob/develop/Figures/moisture_source_20141017000000.png" width=50% height=50%>
 </p>

**2- Moisture sink pattern (06/10/2014)**

 <p align="center" width="100%">
 <img src="https://github.com/tramo-ephyslab/TROVA-master/blob/develop/Figures/moisture_sink_20141006000000.png" width=50% height=50%>
 </p>

**3- Water vapor residence time (17/10/2014)**

In addition, we show the values ​​of the water vapor residence time applying the methodology of Läderach and Sodemann (2016) for October 17, 2014.

<p align="center" width="100%">
 <img src="https://github.com/tramo-ephyslab/TROVA-master/blob/develop/Figures/WVRT_plot_20141017000000.png" width=50% height=50%>
 </p>

## Climatological analysis

**1- TROVA is used in backward in time to determine moisture sources.**

 In this case, it is the moisture source pattern associated with the Iberian Peninsula for the month of October 2001. In this analysis, TROVA is used with the methodology of Sodemann et al. (2008) and as input data the outputs of the FLEXPART-WRF model forced with the Community Earth System Model 2 (CESM2) climate model. The mask used is represented in red (the Iberian Peninsula itself).

<p align="center" width="100%">
 <img src="https://github.com/tramo-ephyslab/TROVA/blob/main/Figures/Fig1-git.jpg" width=50% height=50%>
 </p>

**2- TROVA is used in forward in time to determine the moisture sinks.**

The following Figure presents the moisture sink pattern associated with the Mediterranean Sea for the month of October 2014. The methodology of Stohl and James (2005) is considered and how it masks the geographical limits of the Mediterranean Sea. The input data for TROVA are the outputs of FLEXPART forced with ERA-Interim.


<p align="center" width="100%">
 <img src="https://github.com/tramo-ephyslab/TROVA/blob/main/Figures/Fig2-git.png" width=50% height=50%>
 </p>

# Important notes

This code is not bug-free. Please report any bugs through 'Issues': https://github.com/tramo-ephyslab/TROVA/issues

## Contact and support

José Carlos Fernández Alvarez (jose.carlos.fernandez.alvarez@uvigo.es, jcfernandez@cesga.es) and Albenis Pérez Alarcón (albenis.perez.alarcon@uvigo.es).

# LICENSE

Copyright 2025 Fernández-Alvarez et al. (2022)

This software is published under the GPLv3 license. This means:

- Anyone can copy, modify and distribute this software.
- You have to include the license and copyright notice with each and every distribution.
- You can use this software privately.
- You can use this software for commercial purposes.
- If you dare build your business solely from this code, you risk open-sourcing the whole code base.
- If you modify it, you have to indicate changes made to the code.
- Any modifications of this code base MUST be distributed with the same license, GPLv3.
- This software is provided without warranty.
- The software author or license can not be held liable for any damages inflicted by the software.

# References

[1] Stohl A, James PA. A Lagrangian analysis of the atmospheric branch of the global water cycle: Part II: Earth’s river catchments ocean basins, and moisture transports between them. J. Hydrometeorol. 2005; 6:961–984. https://doi.org/10.1175/JHM470.1.

[2] Sodemann H, Schwierz C, Wernli H. Interannual variability of Greenland winter precipitation sources: Lagrangian moisture diagnostic and North Atlantic Oscillation influence. J. Geophys. Res.-Atmos. 2008; 113:D03107. https://doi.org/10.1029/2007JD008503. 

[3] Piso I et al. The Lagrangian particle dispersion model FLEXPART version 10.3. Geosci. Model Dev. Discuss. 2019; https://doi.org/10.5194/gmd-2018-333.

[4] Brioude J et al. The Lagrangian particle dispersion model FLEXPART-WRF version 3.1. Geosci. Model Dev. 2013; 6:1889–1904. https://doi.org/10.5194/gmd-6-1889-2013.

[5] Fernández-Alvarez, J. C., Pérez-Alarcón, A., Nieto, R., & Gimeno, L. (2022). TROVA: TRansport of water VApor. SoftwareX, 20, 101228. https://doi.org/10.1016/j.softx.2022.101228.