Steps to run the TROVA software:

1- Have the anaconda3 environment installed with the default packages and also netcdf4 and mpi4py.

2- Enter the src/ directory and execute the install.sh code. This creates the executable of the fortran functions used by python for numerical computation 

3- Modify the configuration input file of the TROVA software.

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
