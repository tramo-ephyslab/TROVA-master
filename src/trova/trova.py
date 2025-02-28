#!/usr/bin/env python

# -*- coding: utf-8 -*-
#
# MAIN SCRIPT OF TROVA
# 
# This file is part of TROVA, 
# originally created by José Carlos Fernández Alvarez, Albenis Pérez Alarcón, Raquel Nieto and Luis Gimeno
# at the EPhysLab, University of Vigo. Currently updated by Jose C. Fernández Alvarez of the Galicia 
# Supercomputing Center. In addition, it has contributions from predoctoral students Yiying Wang and 
# Gleisis Alvarez Socorro from the University of Vigo.
# 
# https://github.com/tramo-ephyslab/TROVA-master/
# 
# TROVA is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation v3.
#
# TROVA is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with TROVA. If not, see <http://www.gnu.org/licenses/>.
#

import numpy as np
import sys
from time import time
import time as pytime
import struct
import datetime
from datetime import datetime, timedelta
from netCDF4 import Dataset
import netCDF4
import os
import scipy.interpolate as interp
import functools
from pathlib import Path
import warnings
import gzip
import shutil
from numpy import dtype
import matplotlib.pylab as plt
from mpl_toolkits.mplot3d import Axes3D
import math
import importlib.machinery
from mpi4py import MPI
from .functions import k_dq as K_dq
#from .tensor_operations import k_dq as K_dq
from .functions import k_dq_layers as K_dq_layers
from .functions import read_binary_file as RBF
from .functions import len_file as lf
from .functions import k_dq_so as K_dq_So
from .functions import k_dq_por as K_dq_POR
from .functions import filter_part as Filter_Part
from .functions import filter_part2 as Filter_Part2
from .functions import filter_part_by_height as Filter_by_Height
from .functions import determined_id as D_id
from .functions import search_row as sRow
#from .functions import kdif as kdif
warnings.filterwarnings("ignore", category=DeprecationWarning)
print = functools.partial(print, flush=True)

def Kdif_python(matrix1, matrix2, paso):
    """
    Computes the difference between two matrices based on a given step.

    Parameters:
    matrix1 (numpy.ndarray): The first input matrix.
    matrix2 (numpy.ndarray): The second input matrix.
    paso (float): The step value, can be either -1 (backward) or 1 (forward).

    Returns:
    numpy.ndarray: The output matrix with computed differences and selected values from the input matrices.
    """
    dx, dy = matrix1.shape
    output = np.full((dx, dy - 1), -999.9, dtype=matrix1.dtype)
    
    valid_mask1 = (matrix1[:, 0] != -999.9)
    valid_mask2 = (matrix2[:, 0] != -999.9)
    
    if paso == -1.:
        valid_indices = np.where(valid_mask1 & valid_mask2)[0]
        output[valid_indices, 2] = matrix2[valid_indices, 3] - matrix1[valid_indices, 3]
        output[valid_indices, 1] = matrix1[valid_indices, 2]
        output[valid_indices, 0] = matrix1[valid_indices, 1]
        output[valid_indices, 3] = matrix1[valid_indices, 4]
    
    elif paso == 1.:
        valid_indices = np.where(valid_mask1 & valid_mask2)[0]
        output[valid_indices, 2] = matrix2[valid_indices, 3] - matrix1[valid_indices, 3]
        output[valid_indices, 1] = matrix2[valid_indices, 2]
        output[valid_indices, 0] = matrix2[valid_indices, 1]
        output[valid_indices, 3] = matrix2[valid_indices, 4]
    
    return output

def search_row_python(matrix, lista):
    """
    Searches for rows in a matrix that match values in a given list.

    This function takes a matrix and a list of values, and returns a new matrix
    where each row corresponds to a row in the input matrix whose first element
    matches a value in the list. If a value from the list is not found in the matrix,
    the corresponding row in the output matrix is filled with -999.9.

    Parameters:
    matrix (numpy.ndarray): The input matrix to search within.
    lista (list or numpy.ndarray): The list of values to search for in the first column of the matrix.

    Returns:
    numpy.ndarray: A new matrix with rows from the input matrix that match the values in the list.
                   Rows that do not match are filled with -999.9.
    """
    matrix = np.asarray(matrix)
    lista = np.asarray(lista)
    
    output = np.full((len(lista), matrix.shape[1]), -999.9, dtype=matrix.dtype)
    matrix_dict = {int(row[0]): row for row in matrix}
    
    for i, val in enumerate(lista):
        if int(val) in matrix_dict:
            output[i, :] = matrix_dict[int(val)]
    return output

def determined_id_python(value_mascara, value_mask):
    """
    Determines the indices of elements in a mask that match a given value.

    This function takes a mask array and a value to match, and returns an array
    where each element corresponds to the index of the matching value in the mask.
    If a value from the mask is not found, the corresponding element in the output
    array is filled with -999.

    Parameters:
    value_mascara (numpy.ndarray): The mask array to search within.
    value_mask (int): The value to search for in the mask array.

    Returns:
    numpy.ndarray: An array with indices of the matching values in the mask.
                   Indices that do not match are filled with -999.
    """
    
    value_mascara = np.asarray(value_mascara)
    vector = np.full(value_mascara.shape, -999, dtype=int)
    matching_indices = np.where(value_mascara == value_mask)[0]
    vector[matching_indices] = matching_indices
    return vector

def check_paths(pfile, path):
    """
    Checks if a given path attribute exists in the provided file object.

    This function attempts to retrieve the value of a specified path attribute
    from a given file object. If the attribute does not exist, it returns an
    empty string.

    Parameters:
    pfile (object): The file object to check for the path attribute.
    path (str): The name of the path attribute to retrieve.

    Returns:
    str: The value of the path attribute if it exists, otherwise an empty string.
    """
    try:
        fpath = getattr(pfile, path)
    except:
        fpath = ""
    
    return fpath

def str2boolean(arg):
    """
    Converts a string representation of truth to a boolean value.

    This function takes a string argument and returns its corresponding boolean value.
    It recognizes several common string representations of true and false values.

    Parameters:
    arg (str): The string to convert to a boolean value. Recognized true values are
               "yes", "true", "t", "y", "1". Recognized false values are "no", "false",
               "f", "n", "0". The comparison is case-insensitive.

    Returns:
    bool: The boolean value corresponding to the input string.

    Raises:
    argparse.ArgumentTypeError: If the input string does not match any recognized true or false values.
    """

    if isinstance(arg, bool):
        return arg
    if arg.lower() in ("yes", "true", "t", "y", "1"):
        return True
    elif arg.lower() in ("no", "false", "f", "n", "0"):
        return False
    else:
        raise argparse.ArgumentTypeError("Boolean value expected.")

def ProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '»', printEnd = "\r"):
    """
    Displays a progress bar in the terminal.

    This function prints a progress bar to the terminal to indicate the progress of a task.
    The progress bar updates with each iteration and shows the percentage of completion.

    Parameters:
    iteration (int): Current iteration (must be between 0 and total).
    total (int): Total number of iterations.
    prefix (str): Prefix string (optional).
    suffix (str): Suffix string (optional).
    decimals (int): Positive number of decimals in percent complete (optional).
    length (int): Character length of the bar (optional).
    fill (str): Bar fill character (optional).
    printEnd (str): End character (e.g. "\r", "\r\n") (optional).

    Returns:
    None
    """

    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end = "\r")
    if iteration == total:
        print()


def get_currentversion():
    """
    Retrieves the current version of the TROVA software.

    This function reads the version information from the VERSION file located
    in the same directory as the script and returns it as a string.

    Returns:
    str: The current version of the TROVA software.
    """

    pathpkg = os.path.dirname(__file__)
    version_file = pathpkg+"/VERSION"

    with open(version_file) as vfile:
        version = vfile.readlines()[0].strip()
    return(version)


def get_lastupdate():
    """
    Retrieves the last update date of the TROVA software.

    This function reads the last update date from the LAST_UPDATE file located
    in the same directory as the script and returns it as a string.

    Returns:
    str: The last update date of the TROVA software.
    """

    pathpkg = os.path.dirname(__file__)
    update_file = pathpkg+"/LAST_UPDATE"

    with open(update_file) as ufile:
        lastupdate = ufile.readlines()[0].strip()
    return(lastupdate)


def plotting_tracks_3d(particle_positions, fname):
    """
    Plots 3D tracks of parcels.

    This function creates a 3D plot of parcel tracks using their positions and saves the plot to a file.

    Parameters:
    particle_positions (numpy.ndarray): Array containing the positions of the parcels.
    fname (str): The filename to save the plot.

    Returns:
    None
    """
    fig = plt.figure(figsize=(15,15))
    ax = fig.add_subplot(111, projection='3d')
    ProgressBar(0, particle_positions.shape[1], prefix=" Plotting 3D parcels' tracks       ------>", suffix='', decimals=1, length=40, printEnd="\r")
    for i in range(0, particle_positions.shape[1]):
        lat = particle_positions[:, i, 2]
        lon = particle_positions[:, i, 1]
        ids = particle_positions[:, i, 0]
        z = particle_positions[:, i, 4]

        if ids[0] != -999.9:
            if all(value >= -180 for value in lon):
                ax.plot3D(lon, lat, z / 1000)

        pytime.sleep(0.0001)
        ProgressBar(i + 1, particle_positions.shape[1], prefix=" Plotting 3D parcels' tracks       ------>", suffix='', decimals=1, length=40, printEnd="\r")
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)

    plt.xlabel("Longitude", fontsize=25, labelpad=25)
    plt.ylabel("Latitude", fontsize=25, labelpad=25)
    ax.set_zlabel("Height (km)", fontsize=25, labelpad=25)

    ax.tick_params(axis='z', labelsize=25)

    plt.savefig(fname, dpi=600)
    plt.close()

def ploting_parcels_tracks_map(particle_positions, maps_limits, paso, lat_masked, lon_masked, mascara, value_mask, fname):
    """
    Plots parcel tracks on a map.

    This function creates a 2D map plot of parcel tracks using their positions and saves the plot to a file.

    Parameters:
    particle_positions (numpy.ndarray): Array containing the positions of the parcels.
    maps_limits (list): List containing the map limits [latmin, lonmin, latmax, lonmax, center, dlat, dlon].
    paso (int): Step value indicating the direction of the plot (-1 for backward, 1 for forward).
    lat_masked (numpy.ndarray): Array containing the masked latitudes.
    lon_masked (numpy.ndarray): Array containing the masked longitudes.
    mascara (numpy.ndarray): Array containing the mask values.
    value_mask (int): The value to use for the mask.
    fname (str): The filename to save the plot.

    Returns:
    None
    """
    plt.figure(figsize=(18,12))
    mapa, crs = create_map(maps_limits)

    ProgressBar(0, particle_positions.shape[1], prefix=" Plotting parcels' tracks on a map   ---->", suffix='', decimals=1, length=40, printEnd="\r")
    for i in range(0, particle_positions.shape[1]):
        lat = particle_positions[:, i, 2]
        lon = particle_positions[:, i, 1]
        ids = particle_positions[:, i, 0]

        if ids[0] != -999.9:
            if all(value >= -180 for value in lon):
                mapa.plot(lon, lat, transform=crs)

                if paso == -1:
                    mapa.plot(lon[0], lat[0], color="k", marker="o", markersize=5, transform=crs)
                if paso == 1:
                    mapa.plot(lon[-1], lat[-1], color="k", marker="o", markersize=5, transform=crs)
        pytime.sleep(0.00001)
        ProgressBar(i + 1, particle_positions.shape[1], prefix=" Plotting parcels' tracks on a map   ---->", suffix='', decimals=1, length=40, printEnd="\r")
    mapa.contour(lon_masked, lat_masked, mascara, value_mask, colors="b", linewidths=4)
    plt.savefig(fname, bbox_inches="tight", dpi=600)
    plt.close()

def create_map(maps_limits):
    """
    Creates a map with specified limits.

    This function creates a map with the given latitude and longitude limits and returns the map and its coordinate reference system (CRS).

    Parameters:
    maps_limits (list): List containing the map limits [latmin, lonmin, latmax, lonmax, center, dlat, dlon].

    Returns:
    tuple: A tuple containing the map and its CRS.
    """
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    import matplotlib.ticker as mticker
    from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
    from cartopy.io.shapereader import Reader
    from cartopy.feature import ShapelyFeature

    latmin = maps_limits[0]
    lonmin = maps_limits[1]
    latmax = maps_limits[2]
    lonmax = maps_limits[3]
    center = maps_limits[4]
    dlat = int(maps_limits[5])
    dlon = int(maps_limits[6])

    crs = ccrs.PlateCarree()
    mapa = plt.subplot(111, projection=ccrs.PlateCarree(center))
    mapa.add_feature(cfeature.COASTLINE.with_scale('10m'), linewidth=0.95)
    mapa.add_feature(cfeature.BORDERS, linestyle="-", linewidth=0.75)
    mapa.add_feature(cfeature.LAKES, alpha=0.5)
    mapa.set_extent([lonmin, lonmax, latmin, latmax], crs=ccrs.PlateCarree())

    gl = mapa.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5, color='black', alpha=1, linestyle='--')
    if int(lonmin) != int(lonmax):
        lons = np.arange(int(math.ceil(lonmin)), int(math.floor(lonmax)) + dlon, dlon)
    else:
        lons = np.arange(lonmin, lonmax + dlon, dlon)

    gl_lon_info = []
    if center == 180:
        for clons in lons:
            if clons < 180:
                gl_lon_info = np.append(gl_lon_info, clons)
            else:
                gl_lon_info = np.append(gl_lon_info, clons - 360)
    elif center == 0:
        gl_lon_info = lons

    gl_loc = [True, False, False, True]
    if float(sys.version[0:3]) >= 3.7:
        gl.left_labels = gl_loc[0]
        gl.right_labels = gl_loc[1]
        gl.top_labels = gl_loc[2]
        gl.bottom_labels = gl_loc[3]
    else:
        gl.ylabels_left = gl_loc[0]
        gl.ylabels_right = gl_loc[1]
        gl.xlabels_top = gl_loc[2]
        gl.xlabels_bottom = gl_loc[3]

    gl.xlocator = mticker.FixedLocator(gl_lon_info)
    if dlat >= 1:
        gl.ylocator = mticker.MultipleLocator(dlat)
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 25, 'color': 'k'}
    gl.ylabel_style = {'size': 25, 'color': 'k'}
    return mapa, crs

def plotting_parcels_within_target_region(particle_positions, maps_limits, paso, lat_masked, lon_masked, mascara, value_mask, fname):
    """
    Plots parcels within the target region on a map.

    This function creates a 2D map plot of parcels within the target region using their positions and saves the plot to a file.

    Parameters:
    particle_positions (numpy.ndarray): Array containing the positions of the parcels.
    maps_limits (list): List containing the map limits [latmin, lonmin, latmax, lonmax, center, dlat, dlon].
    paso (int): Step value indicating the direction of the plot (-1 for backward, 1 for forward).
    lat_masked (numpy.ndarray): Array containing the masked latitudes.
    lon_masked (numpy.ndarray): Array containing the masked longitudes.
    mascara (numpy.ndarray): Array containing the mask values.
    value_mask (int): The value to use for the mask.
    fname (str): The filename to save the plot.

    Returns:
    None
    """
    plt.figure(figsize=(18,12))
    mapa, crs = create_map(maps_limits)

    if paso == -1:
        idx = -1
    if paso == 1:
        idx = 0
    ProgressBar(0, particle_positions.shape[1], prefix=" Plotting parcels within the target region", suffix='', decimals=1, length=40, printEnd="\r")
    for i in range(0, particle_positions.shape[1]):
        lat = particle_positions[idx, i, 2]
        lon = particle_positions[idx, i, 1]
        ids = particle_positions[idx, i, 0]

        if ids != -999.9:
            if lon >= -180:
                mapa.plot(lon, lat, color="r", marker="o", markersize=3, transform=crs)
        pytime.sleep(0.001)
        ProgressBar(i + 1, particle_positions.shape[1], prefix=" Plotting parcels within the target region", suffix='', decimals=1, length=40, printEnd="\r")
    mapa.contour(lon_masked, lat_masked, mascara, value_mask, colors="b", linewidths=4)
    plt.savefig(fname, bbox_inches="tight", dpi=600)
    plt.close()

def plot_moisture_sink_source(lon, lat, data, paso, path_output, folder, limit_plot):
    """
    Plots moisture sources or sinks on a map.

    This function creates a contour plot of moisture sources or sinks based on the input data and saves the plot to a file.

    Parameters:
    lon (numpy.ndarray): Array of longitude values.
    lat (numpy.ndarray): Array of latitude values.
    data (numpy.ndarray): Array of data values to plot.
    paso (int): Step value indicating the direction of the simulation (-1 for backward, 1 for forward).
    path_output (str): Path to save the output plot.
    folder (str): Folder name to save the output plot.
    limit_plot (list): List containing the plot limits [latmin, lonmin, latmax, lonmax].

    Returns:
    None
    """
    
    a, b, c = limit_plot[0], limit_plot[1], limit_plot[2]

    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    import matplotlib.ticker as mticker
    from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
    from cartopy.io.shapereader import Reader
    from cartopy.feature import ShapelyFeature
    import warnings
    warnings.filterwarnings("ignore", category=UserWarning, module="cartopy.mpl.gridliner")

    lon, lat = np.meshgrid(lon, lat)
    fig = plt.figure(figsize=(14, 10))
    ax1 = plt.subplot(111, projection=ccrs.PlateCarree())
    ax1.set_extent([np.min(lon), np.max(lon), np.min(lat), np.max(lat)], crs=ccrs.PlateCarree())
    ax1.add_feature(cfeature.COASTLINE.with_scale('10m'), linewidth=0.9)    
    
    if paso == -1:
        band_a = 0.001
        data[np.abs(data) < band_a] = np.nan
        unit = "E-P >0 (mm/day)"
        levels=np.arange(a, b, c)
        pallete = "Reds"
        name = "moisture_source"

    if paso == 1:
        band_a = 0.001
        data[data>0]=np.nan
        data = np.abs(data)
        data[np.abs(data) < band_a] = np.nan
        unit = "P-E >0 (mm/day)"
        levels=np.arange(a, b, c)
        pallete = "Blues"
        name = "moisture_sink"

    cf = ax1.contourf(lon, lat, data, levels, cmap=plt.get_cmap(pallete), extend='both', transform=ccrs.PlateCarree())
    cb = plt.colorbar(cf, orientation="horizontal", pad=0.06, shrink=0.8)
    cb.set_label(label=unit, size=20)
    cb.ax.tick_params(labelsize=20)
    
    gl = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='gray', alpha=0.5, linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_left = True
    gl.ylabels_right = False
    gl.xlines = True
    
    paso_h = 30
    dlat = 10
    gl.xlocator = mticker.FixedLocator(np.arange(-180, 181, paso_h))
    gl.ylocator = mticker.MultipleLocator(dlat)
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 15, 'color': 'black'}
    gl.ylabel_style = {'size': 15, 'color': 'black'}

    plt.savefig(path_output+folder+"/"+name+"_"+folder + ".png", bbox_inches='tight', dpi=300)
    plt.close()

def generate_fecha_simulation(ndias, cyear, cmonth, cday, chours, cminutes):
    """
    Generates a list of simulation dates.

    This function generates a list of dates for the simulation based on the number of days and the initial date and time components provided.

    Parameters:
    ndias (int): Number of days for the simulation.
    cyear (int or list): Initial year(s) of the simulation.
    cmonth (int or list): Initial month(s) of the simulation.
    cday (int or list): Initial day(s) of the simulation.
    chours (int or list): Initial hour(s) of the simulation.
    cminutes (int or list): Initial minute(s) of the simulation.

    Returns:
    tuple: A tuple containing lists of years, months, days, hours, and minutes for the simulation dates.
    """
    nhour = int(ndias * 24)
    year = []
    mes = []
    dia = []
    hora = []
    mins = []
    array = np.arange(0, nhour, 24)

    if not isinstance(chours, list):
        chours = [chours]
    if not isinstance(cminutes, list):
        cminutes = [cminutes]
    if not isinstance(cyear, list):
        cyear = [cyear]
    if not isinstance(cmonth, list):
        cmonth = [cmonth]
    if not isinstance(cday, list):
        cday = [cday]

    for i in array:
        for yy in cyear:
            yy = str(int(yy)).zfill(4)
            for mm in cmonth:
                mm = str(int(mm)).zfill(2)
                for dd in cday:
                    dd = str(int(dd)).zfill(2)
                    for hh in chours:
                        for mmin in cminutes:
                            fecha = yy + "-" + mm + "-" + dd + " " + str(int(hh)).zfill(2) + ":" + str(int(mmin)).zfill(2) + ":00"
                            a = str(time_calc(fecha, float(i)))
                            var1 = a.split(" ")
                            var11 = var1[0].split("-")
                            var12 = var1[1].split(":")
                            year_ = str(var11[0])
                            mes_ = str(var11[1])
                            dia_ = str(var11[2])
                            hora_ = str(var12[0])
                            minn_ = str(var12[1])
                            year.append(year_)
                            mes.append(mes_)
                            dia.append(dia_)
                            hora.append(hora_)
                            mins.append(minn_)

    return year, mes, dia, hora, mins


def function(latitude, longitude, var, var_layers, use_vlayers, vlayers, method, varpor, filename, path, name_var, unit_var, date_save):
    """
    Creates a NetCDF file with the given data.

    This function creates a NetCDF file with the specified latitude, longitude, variable data, and other attributes.

    Parameters:
    latitude (numpy.ndarray): Array of latitude values.
    longitude (numpy.ndarray): Array of longitude values.
    var (numpy.ndarray): Array of variable data.
    var_layers (numpy.ndarray): Array of variable data for layers.
    use_vlayers (bool): Whether to use vertical layers.
    vlayers (list): List of vertical layers.
    method (int): Method used for processing.
    varpor (numpy.ndarray): Array of variable data for sources contribution.
    filename (str): Name of the output file.
    path (str): Path to save the output file.
    name_var (str): Name of the variable.
    unit_var (str): Unit of the variable.
    date_save (numpy.ndarray): Array of dates for the time dimension.

    Returns:
    None
    """
    ncout = Dataset(path + filename + ".nc", 'w', format='NETCDF4')
    ncout.createDimension('lat', len(latitude))
    ncout.createDimension('lon', len(longitude))
    ncout.createDimension('time', len(var[:, 0, 0]))
    time_var = np.arange(len(var[:, 0, 0]))

    if use_vlayers:
        ncout.createDimension('layers', len(vlayers) - 1)
        layer_vec = []
        for ilayer in range(0, len(vlayers) - 1):
            clayer = str(int(vlayers[ilayer])) + "_" + str(int(vlayers[ilayer + 1]))
            layer_vec = np.append(layer_vec, clayer)

    lat = ncout.createVariable('lat', np.dtype('float64').char, ('lat'))
    lat.standard_name = 'latitude'
    lat.long_name = 'latitude'
    lat.units = 'degrees'
    lat.axis = 'Y'

    lon = ncout.createVariable('lon', np.dtype('float64').char, ('lon'))
    lon.standard_name = 'longitude'
    lon.long_name = 'longitude'
    lon.units = 'degrees'
    lon.axis = 'X'

    time = ncout.createVariable('time', np.dtype('float64').char, ('time'))
    time.standard_name = 'time'
    time.long_name = 'time'
    time.units = 'day'
    time.axis = 't'
    time.calendar = "gregorian"
    time.description = "days since 1900-01-01"
    time.units = "days since 1900-01-01"

    vout = ncout.createVariable(name_var, np.dtype('float').char, ('time', 'lat', 'lon'), zlib=True)
    vout.long_name = name_var
    vout.units = unit_var
    vout.standard_name = name_var
    vout.coordinates = "time, lat, lon"

    if name_var == "E_P":
        vvout = ncout.createVariable("E_P_integrated", np.dtype('float').char, ('lat', 'lon'), zlib=True)
        vvout.long_name = "E_P integrated for ndays considered"
        vvout.units = unit_var
        vvout.standard_name = "E_P_integrate"
        vvout.coordinates = "lat, lon"
    vout.original_name = name_var

    if method == 2:
        voutpor = ncout.createVariable("POR", np.dtype('float').char, ('time', 'lat', 'lon'), zlib=True)
        voutpor.long_name = "Sources Contribution for each parcel"
        voutpor.units = "%"
        voutpor.standard_name = "Sources Contribution for each parcel"
        voutpor.coordinates = "time, lat, lon"
        voutpor.original_name = "Sources Contribution for each parcel"
        voutpor.coordinates = "time, lat, lon"

    if use_vlayers:
        voutlayer = ncout.createVariable(name_var + "_layers", np.dtype('float').char, ('time', 'layers', 'lat', 'lon'), zlib=True)
        voutlayer.long_name = name_var + "_layers"
        voutlayer.units = unit_var
        voutlayer.standard_name = name_var
        voutlayer.coordinates = "time, layers, lat, lon"
        voutlayer.original_name = name_var + "_layers"

        vvoutlayer = ncout.createVariable("E_P_integrated_layers", np.dtype('float').char, ('layers', 'lat', 'lon'), zlib=True)
        vvoutlayer.long_name = "E_P integrated for ndays in vertical layers"
        vvoutlayer.units = unit_var
        vvoutlayer.standard_name = "E_P_integrated_by_layers"
        vvoutlayer.coordinates = "layers, lat, lon"
        vvoutlayer.original_name = "E_P integrated for ndays in vertical layers"

        v_layers = ncout.createVariable('vertical_layers', 'str', 'layers')
        v_layers.long_name = 'vertical layers'
        v_layers.units = 'vertical layers'
        v_layers.axis = 'vertical layers'

    lon[:] = longitude
    lat[:] = latitude
    time[:] = date_save[:]
    vout[:, :, :] = var
    if name_var == "E_P":
        vvout[:, :] = np.sum(var, axis=0)

    if method == 2:
        voutpor[:] = varpor

    if use_vlayers:
        voutlayer[:, :, :, :] = var_layers
        vvoutlayer[:, :, :] = np.sum(var_layers, axis=0)
        v_layers[:] = layer_vec

    # Add global attributes
    ncout.setncatts({
        'Institution': 'Galicia Supercomputing Center (CESGA) and Environmental Physics Laboratory (EPhysLab), Centro de Investigación Mariña, Universidade de Vigo, Spain',
        'Author': 'José Carlos Fernández Alvarez',
        'Documentation': 'https://trova-docs.readthedocs.io/en/latest/',
        'Code origin': 'José Carlos Fernández Alvarez et al. 2022, CESGA-UVIGO, Spain',
        'Application': 'TROVA: TRansport Of water VApor'
    })
    ncout.close()

def write_nc(dates, tensor, vartype, filename="output"):
    """
    Writes data to a NetCDF file.

    This function writes the given tensor data to a NetCDF file with the specified filename and variable type.

    Parameters:
    dates (numpy.ndarray): Array of dates for the time dimension.
    tensor (numpy.ndarray): Tensor data to be written to the NetCDF file.
    vartype (str): Type of variable data (e.g., "partpos" or "dqdt").
    filename (str): Name of the output file (default is "output").

    Returns:
    None
    """
    ncout = Dataset(filename + ".nc", 'w', format='NETCDF4')
    ncout.history = 'Parcels positions'
    if vartype == "partpos":
        ncout.history = "partpos[:,:,0] - parcel ID, partpos[:,:,1] - longitude, partpos[:,:,2] - latitude, partpos[:,:,3] - specific humidity, partpos[:,:,4] - vertical position (m), partpos[:,:,5] - topography high (m), partpos[:,:,6] - density (kg/m3), partpos[:,:,7] - PBL high (m), partpos[:,:,8] - Tropopause high (m), partpos[:,:,9] - temperature (K), partpos[:,:,10] - parcel mass (kg)"
    if vartype == "dqdt":
        ncout.history = "partpos[:,:,0] - longitude, partpos[:,:,1] - latitude, partpos[:,:,2] - dq/dt, partpos[:,:,3] - vertical position (m), partpos[:,:,4] - parcel ID, partpos[:,:,5] - specific humidity at starting tracking point (time t0)"
    ndates = len(dates)
    npart = tensor.shape[1]
    vprop = tensor.shape[2]

    ncout.createDimension('time', ndates)
    ncout.createDimension('npart', npart)
    ncout.createDimension('properties', vprop)

    times = ncout.createVariable('times', np.dtype('float64').char, ('time'))
    times.standard_name = 'times'
    times.long_name = 'times'
    times.units = 'day'
    times.axis = 't'
    times.calendar = "gregorian"
    times.description = "days since 1900-01-01"
    times.units = "days since 1900-01-01"

    parts = ncout.createVariable('parcels', dtype('float32').char, ('npart'))
    parts.standard_name = 'parcels'
    parts.long_name = 'Parcels IDs'
    parts.units = ''
    parts.axis = ''

    vout = ncout.createVariable('partpos', dtype('float32').char, ('time', 'npart', 'properties'), zlib=True)
    vout.long_name = 'Parcels position'
    vout.units = ''
    vout.standard_name = "Parcels position"
    vout.coordinates = "times, npart, properties"
    vout.original_name = "Parcels position"

    times[:] = dates
    parts[:] = tensor[0, :, 0]
    vout[:] = tensor[:, :, :]

    # Add global attributes
    ncout.setncatts({
        'Institution': 'Galicia Supercomputing Center (CESGA) and Environmental Physics Laboratory (EPhysLab), Centro de Investigación Mariña, Universidade de Vigo, Spain',
        'Author': 'José Carlos Fernández Alvarez',
        'Documentation': 'https://trova-docs.readthedocs.io/en/latest/',
        'Code origin': 'José Carlos Fernández Alvarez et al. 2022, CESGA-UVIGO, Spain',
        'Application': 'TROVA: TRansport Of water VApor'
    })
    ncout.close()

def create_directory(path):
    """
    Creates a directory if it does not exist.

    This function checks if a directory exists at the specified path, and if not, it creates the directory.

    Parameters:
    path (str): The path of the directory to create.

    Returns:
    None
    """
    try:
        if not os.path.exists(path):
            os.mkdir(path)
    except OSError:
        pass

def read_binaryFile_fortran(filename, type_file,x_left_lower_corner,y_left_lower_corner,x_right_upper_corner,y_right_upper_corner, limit_domain):
    """
    Reads a binary file using Fortran routines.

    This function reads a binary file based on the specified type and domain limits, and returns the data.

    Parameters:
    filename (str): The name of the binary file to read.
    type_file (int): The type of file (1 for FLEXPART-WRF, 2 for FLEXPART-ERAI and FLEXPART-ERA5).
    x_left_lower_corner (float): X-coordinate of the lower left corner of the domain.
    y_left_lower_corner (float): Y-coordinate of the lower left corner of the domain.
    x_right_upper_corner (float): X-coordinate of the upper right corner of the domain.
    y_right_upper_corner (float): Y-coordinate of the upper right corner of the domain.
    limit_domain (int): Whether to limit the domain (1 for yes, 0 for no).

    Returns:
    numpy.ndarray: The data read from the binary file.
    """

    if type_file==1:
        with open(filename,'rb') as inputfile:
            a=b''.join([line for line in inputfile])
        npart=struct.unpack('iiii', a[0:16])
        npart=npart[2]
        data= RBF(filename,npart,x_left_lower_corner,y_left_lower_corner, x_right_upper_corner,y_right_upper_corner, limit_domain)

    if type_file==2:
        len_a=lf(filename)
        npart=((len_a-12)/60)-1
        data= RBF(filename,npart, x_left_lower_corner,y_left_lower_corner, x_right_upper_corner,y_right_upper_corner, limit_domain)

    if limit_domain == 1:
        ind=np.where(data[:, 0]==-999.)
        data = np.delete(data, ind[0], axis=0)
    else:
        data=data
    return data

def load_mask_grid_NR(filename, name_mascara,name_variable_lon, name_variable_lat):
    """
    Loads a mask grid from a NetCDF file.

    This function loads the latitude, longitude, and mask variables from a NetCDF file.

    Parameters:
    filename (str): The name of the NetCDF file to read.
    name_mascara (str): The name of the mask variable in the NetCDF file.
    name_variable_lon (str): The name of the longitude variable in the NetCDF file.
    name_variable_lat (str): The name of the latitude variable in the NetCDF file.

    Returns:
    tuple: A tuple containing the latitude, longitude, and mask arrays.
    """

    wrfile = Dataset(filename)
    lat  = wrfile.variables[name_variable_lat][:]
    lon  = wrfile.variables[name_variable_lon][:]
    mask  = wrfile.variables[name_mascara][:]

    if len(lon.shape)<2:
        lon, lat=np.meshgrid(lon,lat)

    for i in range(0,lon.shape[0]):
       for j in range(0,lon.shape[1]):
          if lon[i,j]>180:
                lon[i,j]=lon[i,j]-360

    return lat, lon,mask

def funtion_interpol_mascara (lat_mascara, lon_mascara, mascara, data):
    """
    Interpolates a mask onto data points.

    This function interpolates the values of a mask onto the given data points using nearest neighbor interpolation.

    Parameters:
    lat_mascara (numpy.ndarray): Array of latitudes for the mask.
    lon_mascara (numpy.ndarray): Array of longitudes for the mask.
    mascara (numpy.ndarray): Array of mask values.
    data (numpy.ndarray): Array of data points to interpolate the mask onto.

    Returns:
    numpy.ndarray: The interpolated mask values at the data points.
    """

    lat_lon=np.empty((len(lat_mascara), 2))
    lat_lon[:,0]=lon_mascara
    lat_lon[:,1]=lat_mascara

    prsInterpu = interp.NearestNDInterpolator(lat_lon,mascara)
    si = np.empty((data[:,1].size, 2))
    si[:,0] = data[:,1]
    si[:,1] = data[:,2]
    result=prsInterpu(si)

    return result

def plot_point_(lat, lon,mascara):
    """
    Plots points on a map using a mask.

    This function creates a scatter plot of points on a map using the given latitude, longitude, and mask values.

    Parameters:
    lat (numpy.ndarray): Array of latitude values.
    lon (numpy.ndarray): Array of longitude values.
    mascara (numpy.ndarray): Array of mask values.

    Returns:
    None
    """

    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    import matplotlib.ticker as mticker
    from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
    import os
    from cartopy.io.shapereader import Reader
    from cartopy.feature import ShapelyFeature
    import matplotlib.pyplot as plt
    import matplotlib

    fig=plt.figure(figsize=(10,8))
    ax1=plt.subplot(111,projection=ccrs.PlateCarree())

    ax1.add_feature(cfeature.COASTLINE.with_scale('10m'), linewidth=0.9)
    ax1.add_feature(cfeature.STATES, edgecolor="black",zorder=10)
    ax1.set_extent([-180.,180.,-90,90], crs=ccrs.PlateCarree())

    mascara=mascara.astype("bool")
    plt.scatter(lon[mascara], lat[mascara], s=10)
    
    gl = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=1, color='gray', alpha=0.2, linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_left = True
    gl.ylabels_right = False
    gl.xlines = True

    paso_h=10
    lons=np.arange(np.ceil(-180),np.ceil(180),paso_h)
    gl.xlocator = mticker.FixedLocator(lons)
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 15, 'color': 'gray'}
    gl.xlabel_style = {'color': 'black'}
    
    plt.savefig("mask.png")
    plt.close()

def funtion_interpol_mascara_2(lat_mascara, lon_mascara, mascara, data):
    """
    Interpolates a mask onto data points.

    This function interpolates the values of a mask onto the given data points using nearest neighbor interpolation.

    Parameters:
    lat_mascara (numpy.ndarray): Array of latitudes for the mask.
    lon_mascara (numpy.ndarray): Array of longitudes for the mask.
    mascara (numpy.ndarray): Array of mask values.
    data (numpy.ndarray): Array of data points to interpolate the mask onto.

    Returns:
    numpy.ndarray: The interpolated mask values at the data points.
    """

    lat_lon=np.empty((len(lat_mascara), 2))
    lat_lon[:,0]=lon_mascara
    lat_lon[:,1]=lat_mascara
    prsInterpu = interp.NearestNDInterpolator(lat_lon,mascara)
    si = np.empty((data[:,1].size, 2))
    si[:,0] = data[:,0]
    si[:,1] = data[:,1]
    result=prsInterpu(si)
    return result

def determine_id_binary_grid_NR_fortran(data, lat_mascara, lon_mascara, value_mascara, value_mask):
    """
    Determines the indices of elements in a binary grid that match a given value using Fortran routines.

    This function interpolates the mask values onto the data points and determines the indices of elements that match the given value.

    Parameters:
    data (numpy.ndarray): Array of data points.
    lat_mascara (numpy.ndarray): Array of latitudes for the mask.
    lon_mascara (numpy.ndarray): Array of longitudes for the mask.
    value_mascara (numpy.ndarray): Array of mask values.
    value_mask (int): The value to search for in the mask array.

    Returns:
    numpy.ndarray: A submatrix of data points that match the given value.
    """
    
    value_mascara = funtion_interpol_mascara(lat_mascara, lon_mascara, value_mascara, data)
    id_vector = determined_id_python(value_mascara.astype(int), value_mask)

    submatrix=[]
    ind=[]
    for ii in id_vector:
        if ii !=-999:
            submatrix=np.append(submatrix, data[ii,:])
            ind.append(ii)
    submatrix=np.reshape(submatrix,(len(ind), 11))
    return submatrix

def search_row_fortran(lista, matrix):
    """
    Searches for rows in a matrix that match values in a given list using Fortran routines.

    This function takes a matrix and a list of values, and returns a new matrix where each row corresponds to a row in the input matrix whose first element matches a value in the list.

    Parameters:
    lista (list or numpy.ndarray): The list of values to search for in the first column of the matrix.
    matrix (numpy.ndarray): The input matrix to search within.

    Returns:
    numpy.ndarray: A new matrix with rows from the input matrix that match the values in the list.
    """
    matrix_ = search_row_python(matrix, lista) #python 
    return matrix_

def calc_A(resolution, lat, lon):
    """
    Calculates the area of grid cells based on latitude and longitude.

    This function calculates the area of each grid cell defined by the given latitude and longitude arrays and the specified resolution.

    Parameters:
    resolution (float): The resolution of the grid cells.
    lat (numpy.ndarray): Array of latitude values.
    lon (numpy.ndarray): Array of longitude values.

    Returns:
    numpy.ndarray: An array of the same shape as the input latitude and longitude arrays, containing the area of each grid cell.
    """

    rt = 6371000.
    gr = np.pi/180.
    a,b=lat.shape
    area=np.empty((a-1,b-1))
    area[:,:]=0
    for j in range(len(lat[0,:])-1):
        for i in range(len(lat[:,0])-1):
            area[i,j]=np.abs((gr*rt**2)*( np.sin(gr*lat[i,j]) - np.sin(gr*lat[i+1,j])))*np.abs(resolution)
    return area

def grid_point (resolution, numPdX, numPdY, x_lower_left, y_lower_left):
    """
    Generates a grid of points based on the specified resolution and domain limits.

    This function generates a grid of latitude and longitude points based on the specified resolution and the coordinates of the lower left corner of the domain.

    Parameters:
    resolution (float): The resolution of the grid cells.
    numPdX (int): Number of grid points in the X direction.
    numPdY (int): Number of grid points in the Y direction.
    x_lower_left (float): X-coordinate of the lower left corner of the domain.
    y_lower_left (float): Y-coordinate of the lower left corner of the domain.

    Returns:
    tuple: A tuple containing two numpy arrays: the latitude and longitude points of the grid.
    """

    lat_new=[]
    lon_new=[]
    lat_min=y_lower_left
    lon_min=x_lower_left
    lat_new=np.append(lat_new, lat_min)
    lon_new=np.append(lon_new, lon_min)

    for i in range(numPdY):
        lat_min=lat_min+resolution
        lat_new=np.append(lat_new, lat_min)
    for j in range(numPdX):
        lon_min=lon_min+resolution
        lon_new= np.append(lon_new, lon_min)
    lon, lat=np.meshgrid(lon_new, lat_new)
    return lat, lon

def grid_plot_final(lat, lon):
    """
    Generates a grid of points for plotting based on the input latitude and longitude arrays.

    This function generates a grid of latitude and longitude points for plotting, based on the input latitude and longitude arrays.

    Parameters:
    lat (numpy.ndarray): Array of latitude values.
    lon (numpy.ndarray): Array of longitude values.

    Returns:
    tuple: A tuple containing two numpy arrays: the latitude and longitude points for plotting.
    """

    lat_new=[]
    lon_new=[]
    for i in range(len(lat[:,0])-1):
        lat_new= np.append(lat_new, (lat[i+1,0]+ lat[i,0])/2.)
    for j in range(len(lon[0,:])-1):
        lon_new= np.append(lon_new, (lon[0,j+1]+ lon[0,j])/2.)
    lon_plot, lat_plot=np.meshgrid(lon_new, lat_new)
    return lat_plot,lon_plot

def time_calc(init_time,h_diff):
    """
    Calculates a new time based on the initial time and a time difference in hours.

    This function calculates a new time by adding the specified time difference in hours to the initial time.

    Parameters:
    init_time (str): The initial time in the format "YYYY-MM-DD HH:MM:SS".
    h_diff (float): The time difference in hours.

    Returns:
    datetime: The calculated time.
    """

    formatted_time = datetime.strptime(init_time, "%Y-%m-%d %H:%M:%S")
    calculated_time=formatted_time+timedelta(hours=h_diff)
    return calculated_time

def time_calcminutes(init_time,h_diff):
    """
    Calculates a new time based on the initial time and a time difference in minutes.

    This function calculates a new time by adding the specified time difference in minutes to the initial time.

    Parameters:
    init_time (str): The initial time in the format "YYYY-MM-DD HH:MM:SS".
    h_diff (float): The time difference in minutes.

    Returns:
    datetime: The calculated time.
    """

    formatted_time = datetime.strptime(init_time, "%Y-%m-%d %H:%M:%S")
    calculated_time=formatted_time+timedelta(minutes=h_diff)
    return calculated_time

def generate_file(paso, dtime, totaltime, fecha, path, key_gz, noleap):
    """
    Generates a list of file names and dates for the simulation.

    This function generates a list of file names and corresponding dates for the simulation based on the specified parameters.

    Parameters:
    paso (int): Step value indicating the direction of the simulation (-1 for backward, 1 for forward).
    dtime (int): Time step in minutes.
    totaltime (int): Total simulation time in minutes.
    fecha (str): Initial date and time in the format "YYYY-MM-DD HH:MM:SS".
    path (str): Path to save the output files.
    key_gz (int): Whether to use gzip compression (1 for yes, 0 for no).
    noleap (int): Whether to exclude leap years (1 for yes, 0 for no).

    Returns:
    tuple: A tuple containing two lists: the list of file names and the list of corresponding dates.
    """
    
    nhour = int(totaltime) + dtime
    list_fecha = []
    listdates = []
    
    if paso in [-1, -2, -3]: 
        array = np.arange(nhour, 0, -dtime)
        for i in array:
            a = str(time_calcminutes(fecha, float(i) * (-1)))
            var1 = a.split(" ")
            var11 = var1[0].split("-")
            var12 = var1[1].split(":")
            fecha_dia = str(var11[0] + var11[1] + var11[2] + var12[0] + var12[1] + var12[2])
            name = path + "partposit_" + fecha_dia
            if key_gz == 1:
                name = path + "partposit_" + fecha_dia + ".gz"
            list_fecha = np.append(list_fecha, name)
            listdates = np.append(listdates, int(fecha_dia))
        
        fecha_ = fecha.split(" ")
        var11 = fecha_[0].split("-")
        var12 = fecha_[1].split(":")
        fecha_dia = str(var11[0] + var11[1] + var11[2] + var12[0] + var12[1] + var12[2])
        if key_gz == 1:
            name = path + "partposit_" + fecha_dia + ".gz"
        else:
            name = path + "partposit_" + fecha_dia
        list_fecha = np.append(list_fecha, name)
        listdates = np.append(listdates, int(fecha_dia))
        
        if noleap == 1:
            ind = []
            for i in range(len(list_fecha)):
                name = list_fecha[i].split("/")[-1]
                date = name.split("_")[-1]
                if date[4:6] == "02" and date[6:8] == "29":
                    ind.append(i)
            if ind:
                list_fecha = np.delete(list_fecha, ind)
                listdates = np.delete(listdates, ind)

                first_date_str = list_fecha[0].split("_")[-1].replace(".gz", "")
                first_date = datetime.strptime(first_date_str, "%Y%m%d%H%M%S")
        
                for _ in range(len(ind)):
                    first_date -= timedelta(minutes=dtime)
                    while first_date.month == 2 and first_date.day == 29:
                        first_date -= timedelta(days=1)
                    new_date_str = first_date.strftime("%Y%m%d%H%M%S")
                    new_name = path + "partposit_" + new_date_str
                    if key_gz == 1:
                        new_name += ".gz"
                    list_fecha = np.insert(list_fecha, 0, new_name)
                    listdates = np.insert(listdates, 0, int(new_date_str))
            
    if paso == 1:
        array =np.arange(0,nhour+dtime, dtime)
        for i in array:
            a=str(time_calcminutes(fecha,float(i)))
            var1=a.split(" ")
            var11=var1[0].split("-")
            var12=var1[1].split(":")
            fecha_dia=str(var11[0]+var11[1]+var11[2]+var12[0]+var12[1]+var12[2])
            name=path+"partposit_"+fecha_dia
            if key_gz==1:
                name=path+"partposit_"+fecha_dia+".gz"
            else:
                name=path+"partposit_"+fecha_dia
            list_fecha=np.append(list_fecha, name)
            listdates=np.append(listdates, int(fecha_dia))
        
        if noleap == 1:
            ind = []
            for i in range(len(list_fecha)):
                name = list_fecha[i].split("/")[-1]
                date = name.split("_")[-1]
                if date[4:6] == "02" and date[6:8] == "29":
                    ind.append(i)
            if ind:
                list_fecha = np.delete(list_fecha, ind)
                listdates = np.delete(listdates, ind)
            
                last_date_str = list_fecha[-1].split("_")[-1].replace(".gz", "")
                last_date = datetime.strptime(last_date_str, "%Y%m%d%H%M%S")
        
                for _ in range(len(ind)):
                    last_date += timedelta(minutes=dtime)
                    while last_date.month == 2 and last_date.day == 29:
                        last_date += timedelta(days=1)
                    new_date_str = last_date.strftime("%Y%m%d%H%M%S")
                    new_name = path + "partposit_" + new_date_str
                    if key_gz == 1:
                        new_name += ".gz"
                    list_fecha = np.append(list_fecha, new_name)
                    listdates = np.append(listdates, new_date_str)

    return list_fecha, listdates

def read_proccesor(lista_partposi,submatrix, rank, x_left_lower_corner,y_left_lower_corner,
               x_right_upper_corner,y_right_upper_corner, model, key_gz, type_file, limit_domain):
    
    """
    Reads and processes binary files in parallel.

    This function reads binary files in parallel using MPI, processes the data, and returns a tensor of the processed data.

    Parameters:
    lista_partposi (list): List of file paths to read.
    submatrix (numpy.ndarray): Submatrix of data points to process.
    rank (int): Rank of the current MPI process.
    x_left_lower_corner (float): X-coordinate of the lower left corner of the domain.
    y_left_lower_corner (float): Y-coordinate of the lower left corner of the domain.
    x_right_upper_corner (float): X-coordinate of the upper right corner of the domain.
    y_right_upper_corner (float): Y-coordinate of the upper right corner of the domain.
    model (str): Model type (e.g., "FLEXPART").
    key_gz (int): Whether to use gzip compression (1 for yes, 0 for no).
    type_file (int): Type of file (1 for FLEXPART-WRF, 2 for FLEXPART-ERAI and FLEXPART-ERA5).
    limit_domain (int): Whether to limit the domain (1 for yes, 0 for no).

    Returns:
    numpy.ndarray: A tensor of the processed data.
    """

    a1=np.arange(len(lista_partposi))
    dx,dy =submatrix.shape
    tensor_local=np.ones((len(lista_partposi),dx,dy))*(-999.9)
    for i in a1:
        print ("Reading | " + model+" -> ",  lista_partposi[i])
        if key_gz==1:
            desc_gz(lista_partposi[i])
            part_post_i=read_binaryFile_fortran(lista_partposi[i][:-3], type_file,x_left_lower_corner,y_left_lower_corner,
               x_right_upper_corner,y_right_upper_corner, limit_domain)
            cmd_rm= "rm -rf "+lista_partposi[i][:-3]
            os.system(cmd_rm)
        else:
            part_post_i=read_binaryFile_fortran(lista_partposi[i], type_file,x_left_lower_corner,y_left_lower_corner,
               x_right_upper_corner,y_right_upper_corner, limit_domain)

        matrix_i=search_row_fortran(submatrix[:,0],part_post_i)         
        tensor_local[i,:,:]=matrix_i
    return tensor_local

def remove_rows_with_value(tensor, tensor_por, idPart, qIni, ref_index=0, value=-999.9):
    """
    Removes rows from a tensor that contain a specified value.

    This function removes rows from the input tensor, tensor_por, idPart, and qIni arrays where any element in the specified reference index row contains the given value.

    Parameters:
    tensor (numpy.ndarray): The input tensor to filter.
    tensor_por (numpy.ndarray): The tensor containing percentage values to filter.
    idPart (numpy.ndarray): Array of parcel IDs to filter.
    qIni (numpy.ndarray): Array of initial specific humidity values to filter.
    ref_index (int): The reference index of the row to check for the specified value (default is 0).
    value (float): The value to check for in the reference row (default is -999.9).

    Returns:
    tuple: A tuple containing the filtered tensor, tensor_por, idPart, and qIni arrays.
    """
    
    if tensor.shape[0] == 0:
        return tensor

    if ref_index >= tensor.shape[0]:
        raise IndexError(f"Reference index {ref_index} is outside the range of the tensor with {tensor.shape[0]} sub-arrays")

    mask = ~(tensor[ref_index] == value).any(axis=-1)
    filtered_tensor = tensor[:, mask, :]
    filtered_tensor_por = tensor_por[:, mask, :]
    filtered_idPart = idPart[mask]
    filtered_qIni = qIni[mask]

    return filtered_tensor, filtered_tensor_por, filtered_idPart, filtered_qIni

def plot_residence_time(residence_time_particles, residence_time_mean, output_dir, date, rank):
    """
    Plot residence time for all particles and display the mean value in the title.

    Parameters:
    residence_time_particles (numpy.ndarray): Array of residence times for each particle.
    residence_time_mean (float): Mean residence time.
    output_dir (str): Directory where the plot will be saved.
    date (str): Date string for output file naming.

    Returns:
    None
    """
    
    import matplotlib.pyplot as plt
    plt.figure(figsize=(10, 6))
    plt.plot(residence_time_particles, 'o', markersize=2)
    plt.axhline(y=residence_time_mean, color='r', linestyle='--', label=f'Mean: {residence_time_mean:.2f} days')
    plt.xlabel('Particle number')
    plt.ylabel('Residence time (days)')
    plt.legend()
    plot_file = f"{output_dir}WVRT_plot_{date}.png"
    plt.savefig(plot_file, bbox_inches='tight', dpi=300)
    plt.close()
    if rank == 0:
        print("--------------------------------------------------------------------")
        print(f"Plot for the residence time for all particles saved to {plot_file}")
        print("--------------------------------------------------------------------")
   
def compute_residence_time_and_save(dqdt, output_dir, date, dtime, totaltime, folder, rank):
    """
    Compute water vapor residence time from a dq/dt tensor and save results to a NetCDF file.
    
    Parameters:
    dqdt (numpy.ndarray): dq/dt tensor.
    output_dir (str): Directory where the output NetCDF file will be saved.
    date (str): Date string for output file naming.
    dtime (int): Time step interval in minutes.
    totaltime (int): Total time in minutes.
    
    Returns:
    None
    """
   
    dtime = dtime / 60  # Convert dtime to hours
    total_hours = totaltime / 60  # Convert totaltime to hours
    total_days = total_hours / 24  # Convert total_hours to days

    # Compute specific humidity change ∆q = (dq/dt) * ∆t
    delta_q = dqdt * dtime

    # Keep only positive ∆q values (increase in humidity)
    delta_q_positive = np.where(delta_q > 0, delta_q, np.nan)  # Set negative values to NaN

    # Compute total increase in specific humidity for each particle
    total_q_increase = np.nansum(delta_q_positive, axis=0) 

    # Compute contribution ratio f_i = ∆q_i / sum(∆q)
    fi = delta_q_positive / total_q_increase 

    # Generate time steps (hours) based on the correct time order
    # Generate time steps from -total_hours to 0 (steps of dtime hours)
    time_steps = np.arange(-(total_hours-dtime), dtime, dtime)[:, np.newaxis]  # Reshape to (n_steps, 1) for broadcasting
    #print(time_steps)
  
    # Compute residence time for each particle τ_i (in hours)
    residence_time_particles = np.nansum(fi * time_steps, axis=0) / 24 

    # Compute the average residence time for the entire region (convert to days)
    residence_time_mean = np.nanmean(residence_time_particles)   # Convert hours to days

    # Print results
    #print(f"Computed average residence time: {residence_time_mean:.2f} days")

    # Create output filename with timestamp
    output_file = f"{output_dir}WVRT_{folder}.nc"

    # Save results to a NetCDF file
    ds_out = Dataset(output_file, 'w', format='NETCDF4')

    # Create dimensions
    ds_out.createDimension('particle', residence_time_particles.shape[0])
    ds_out.createDimension('time', 1)

    # Create variables
    times = ds_out.createVariable('times', np.dtype('float64').char, ('time',))
    times.standard_name = 'times'
    times.long_name = 'times'
    times.units = 'day'
    times.axis = 't'
    times.calendar = "gregorian"
    times.description = "days since 1900-01-01"
    times.units = "days since 1900-01-01"

    residence_time_particles_var = ds_out.createVariable('residence_time_particles', 'f4', ('particle',))
    residence_time_mean_var = ds_out.createVariable('residence_time_mean', 'f4')

    residence_time_particles_var.long_name = "Water vapor residence time of particles"
    residence_time_particles_var.units = "days"
    residence_time_particles_var.standard_name = "residence_time_particles"
    residence_time_particles_var.negative_values = "Values are negative because they are considered in backward mode in time"
    
    residence_time_mean_var.long_name = "Mean water vapor residence time"
    residence_time_mean_var.units = "days"
    residence_time_mean_var.standard_name = "residence_time_mean"
    residence_time_mean_var.negative_values = "Value is negative because it is considered in backward mode in time"

    # Assign data to variables
    residence_time_particles_var[:] = residence_time_particles
    residence_time_mean_var.assignValue(residence_time_mean) 
    date_ = list(map(lambda date: convert_date_to_ordinal(*decompose_date(date)), date))
    times[:] = date_ 

    # Add global attributes
    ds_out.setncatts({
        'Institution': 'Galicia Supercomputing Center (CESGA) and Environmental Physics Laboratory (EPhysLab), Centro de Investigación Mariña, Universidade de Vigo, Spain',
        'Author': 'José Carlos Fernández Alvarez',
        'Documentation': 'https://trova-docs.readthedocs.io/en/latest/',
        'Code origin': 'José Carlos Fernández Alvarez et al. 2022, CESGA-UVIGO, Spain',
        'Application': 'TROVA: TRansport Of water VApor'
    })

    # Close datasets
    ds_out.close()

    plot_residence_time(np.abs(residence_time_particles), np.abs(residence_time_mean), output_dir, str(int(date[0])), rank)

    return residence_time_mean

def _backward_dq(lista_partposi ,file_mask, name_mascara,name_variable_lon, name_variable_lat,lat_f, lon_f,rank,size, comm, type_file,
                 x_left_lower_corner,y_left_lower_corner, x_right_upper_corner,y_right_upper_corner, model, method,threshold,filter_value,
                 value_mask, key_gz, path_output,use_vertical_layers, vertical_layers, filter_parcels_height, filter_vertical_layers,limit_domain, dates,
                 dtime, totaltime, folder, method_wvrt):
    """
    Processes backward parcel tracking data.

    This function processes backward parcel tracking data, filters the data based on specified criteria, and returns the results.

    Parameters:
    lista_partposi (list): List of file paths to read.
    file_mask (str): Path to the mask file.
    name_mascara (str): Name of the mask variable.
    name_variable_lon (str): Name of the longitude variable in the mask file.
    name_variable_lat (str): Name of the latitude variable in the mask file.
    lat_f (numpy.ndarray): Array of latitude values.
    lon_f (numpy.ndarray): Array of longitude values.
    rank (int): Rank of the current MPI process.
    size (int): Total number of MPI processes.
    comm (MPI.Comm): MPI communicator.
    type_file (int): Type of file (1 for FLEXPART-WRF, 2 for FLEXPART-ERAI and FLEXPART-ERA5).
    x_left_lower_corner (float): X-coordinate of the lower left corner of the domain.
    y_left_lower_corner (float): Y-coordinate of the lower left corner of the domain.
    x_right_upper_corner (float): X-coordinate of the upper right corner of the domain.
    y_right_upper_corner (float): Y-coordinate of the upper right corner of the domain.
    model (str): Model type (e.g., "FLEXPART").
    method (int): Method used for processing.
    threshold (float): Threshold value for filtering.
    filter_value (int): Value used for filtering.
    value_mask (int): Value to search for in the mask array.
    key_gz (int): Whether to use gzip compression (1 for yes, 0 for no).
    path_output (str): Path to save the output files.
    use_vertical_layers (bool): Whether to use vertical layers.
    vertical_layers (list): List of vertical layers.
    filter_parcels_height (bool): Whether to filter parcels by height.
    filter_vertical_layers (list): List of vertical layers for filtering.
    limit_domain (int): Whether to limit the domain (1 for yes, 0 for no).
    dates (list): List of dates for the simulation.
    dtime (int): Time step interval in minutes.
    totaltime (int): Total time in minutes.
    folder (str): Folder name to save the output file within the output directory.
    method_wvrt (int): Method used for calculating water vapor residence time.

    Returns:
    tuple: A tuple containing the processed data and additional information.
    """

    name_file=lista_partposi[-1]
    if rank==0:
        print ("Reading | " + model+" -> ",  name_file)
        name_txt_part=name_file.split("/")
        f=open(path_output+name_txt_part[-1].split("_")[-1].split(".")[0]+"/"+name_txt_part[-1].split("_")[-1].split(".")[0]+".txt", "a")
    
    if key_gz==1:
        desc_gz(name_file)
        part_post=read_binaryFile_fortran(name_file[:-3], type_file, x_left_lower_corner,y_left_lower_corner,
               x_right_upper_corner,y_right_upper_corner, limit_domain)
        cmd_rm= "rm -rf "+name_file[:-3]
        os.system(cmd_rm)
    else:
        part_post=read_binaryFile_fortran(name_file, type_file, x_left_lower_corner,y_left_lower_corner,
               x_right_upper_corner,y_right_upper_corner, limit_domain)
        
    lat_masked, lon_masked, mascara=load_mask_grid_NR(file_mask, name_mascara,name_variable_lon, name_variable_lat)
    submatrix=determine_id_binary_grid_NR_fortran(part_post, lat_masked.flatten(), lon_masked.flatten(), mascara.flatten(), value_mask)
    submatrix=submatrix[np.argsort(submatrix[:, 0])]
        
    if filter_parcels_height:
       submatrix, counter_part_height = Filter_by_Height(submatrix,submatrix,-1,filter_vertical_layers[0],filter_vertical_layers[1], len(submatrix[:, 0]))

    if rank==0:
        print ("Reading | " + model+" -> ",  lista_partposi[-2])

    if key_gz==1:
        desc_gz(lista_partposi[-2])
        part_post_i=read_binaryFile_fortran(lista_partposi[-2][:-3], type_file, x_left_lower_corner,y_left_lower_corner,x_right_upper_corner,y_right_upper_corner, limit_domain)
        cmd_rm= "rm -rf "+lista_partposi[-2][:-3]
        os.system(cmd_rm)
    else:
        part_post_i=read_binaryFile_fortran(lista_partposi[-2], type_file, x_left_lower_corner,y_left_lower_corner,x_right_upper_corner,y_right_upper_corner, limit_domain)
    
    matrix_i=search_row_fortran(submatrix[:,0],part_post_i)
    
    dimX, dimY=matrix_i.shape

    if filter_value!=0:
        tmp_matrix=submatrix[:,3] - matrix_i[:,3]
        aux_matrix=np.copy(submatrix)
        aux_matrix[:,3]=tmp_matrix
        omatrix, counter_part=Filter_Part2(aux_matrix,aux_matrix,-1,threshold,len(aux_matrix[:, 0]))
        submatrix[omatrix==-999.9]=-999.9
        matrix_i[omatrix==-999.9]=-999.9

    tensor_t=np.ones((len(lista_partposi)-1,dimX ,dimY ))*(-999.9)
    tensor_t[-1,:,:]=matrix_i

    tensor_org=np.ones((len(lista_partposi),dimX ,dimY ))*(-999.9)
    tensor_org[-1,:,:]=submatrix
    tensor_org[-2,:,:]=matrix_i
    
    n = len(lista_partposi)-2
    count = n // size
    remainder = n % size

    if rank < remainder:
        start = rank * (count + 1)
        stop = start + count + 1
    else:
        start = rank * count + remainder
        stop = start + count

    local_list=lista_partposi[start:stop]
    local_results = np.empty((len(local_list), dimX, dimY))
    local_results = read_proccesor(local_list, submatrix, rank,x_left_lower_corner,y_left_lower_corner,
               x_right_upper_corner,y_right_upper_corner, model, key_gz, type_file, limit_domain)

    if rank > 0:
        comm.Send(local_results, dest=0, tag=14)
    else:
        i_start=[]
        i_stop=[]
        for i in range(size):
            if i < remainder:
                i_start = np.append(i_start,i * (count + 1))
                i_stop = np.append(i_stop,i_start + count + 1)
            else:
                ii_start=i * count + remainder
                ii_stop=ii_start + count
                i_start = np.append(i_start,ii_start)
                i_stop = np.append(i_stop, ii_stop)
        final_results = np.copy(local_results)
        tensor_t[int(i_start[0]):int(i_stop[0]),:,:]= final_results
        tensor_org[int(i_start[0]):int(i_stop[0]),:,:]= final_results
        for i in range(1, size):
            if i < remainder:
                rank_size = count + 1
            else:
                rank_size = count
            tmp = np.empty((rank_size, final_results.shape[1],final_results.shape[2]), dtype=np.float64)
            comm.Recv(tmp, source=i, tag=14)
            tensor_t[int(i_start[i]):int(i_stop[i]),:,:]=tmp
            tensor_org[int(i_start[i]):int(i_stop[i]),:,:]=tmp
    
    comm.Bcast(tensor_t, root=0)
    comm.Bcast(tensor_org, root=0)
        
    matrix_result=np.ones((len(tensor_t[:,0,0])-1, len(submatrix[:,0]),4))*(-999.9)
    a2=np.arange(len(tensor_t[:,0,0])-1)
     
    
    for i in a2[::-1]:
        matrix = Kdif_python (tensor_t[i,:,:5], tensor_t[i+1,:,:5],-1.)
        matrix_result[i,:,2]=matrix[:,2]
        matrix_result[i,:,1]=matrix[:,1]
        matrix_result[i,:,0]=matrix[:,0]
        matrix_result[i,:,3]=matrix[:,3]

    if method_wvrt==1:
        submatrix_wvrt = matrix_result[:, matrix_result[0, :, 0] != -999.9, :]
        wvrt_mean = compute_residence_time_and_save(submatrix_wvrt[:, :, 2], path_output+folder+"/", dates, dtime, totaltime, folder)

    if method==1:
        matrix_result=matrix_result
    matrix_result_por=np.ones((len(tensor_t[:,0,0])-1, len(submatrix[:,0]),3))*(-999.9)

    if method==2:
        dqdt=matrix_result[:,:, 2]
        position_mask=np.empty((dqdt.shape))
        for i in range(dqdt.shape[0]):
            position_mask[i,:]=funtion_interpol_mascara_2 (lat_masked.flatten(), lon_masked.flatten(), mascara.flatten(), matrix_result[i,:,:2])
        dqdt = dqdt.swapaxes(0, 1)
        position_mask=position_mask.swapaxes(0, 1)
        result=K_dq_So(dqdt, position_mask, threshold,len(dqdt[:,0]), len(dqdt[0,:])).swapaxes(0, 1)
        result=np.expand_dims(result, axis=2)
        matrix_result[:,:,2]=result[:,:,0]
        matrix_result_por[:,:,0]=matrix_result[:,:,0]
        matrix_result_por[:,:,1]=matrix_result[:,:,1]
        for i in range(len(matrix_result_por[:,0,0])):
            matrix_result_por[i,:,2]=(matrix_result[i,:, 2]/tensor_t[-1,:,3])*100
    
    matrix_result, matrix_result_por, Filtered_idPart, Filtered_qIni = remove_rows_with_value(matrix_result, matrix_result_por, tensor_t[-1,:,0], tensor_t[-1,:,3], ref_index=matrix_result.shape[0]-1, value=-999.9)
        
    if rank==0:
        print("")
        print ("   + Number of detected parcels within the target region => ", len(matrix_i[:,0]))
        f.write("%s %d\n"%("NumP: ",len(matrix_i[:,0])))
        if filter_parcels_height:
            print ("   + Number of parcels filter by height: Layer [" + str(filter_vertical_layers[0]) + " - " +str(filter_vertical_layers[1]) +" meters] => ", counter_part_height)
            f.write("%s %d\n"%("by height: ",counter_part_height))
        if filter_value!=0:
            print ("   + Number of precipitating parcels within the target region => ", counter_part)
            f.write("%s %d\n"%("Filtered NumP : ",counter_part))
        if method_wvrt:
            print("   + Mean water vapor residence time => %.2f days" % np.abs(wvrt_mean))
            f.write("%s %.2f\n"%("Mean_WVRT: ",np.abs(wvrt_mean)))
    
    return matrix_result, Filtered_idPart, matrix_result_por, Filtered_qIni, tensor_org, lat_masked, lon_masked, mascara

def _forward_dq(lista_partposi ,file_mask, name_mascara,name_variable_lon, name_variable_lat,lat_f, lon_f,rank,size,comm, type_file,
                x_left_lower_corner,y_left_lower_corner, x_right_upper_corner,y_right_upper_corner, model, value_mask, key_gz,path_output,use_vertical_layers, vertical_layers, filter_parcels_height, filter_vertical_layers, limit_domain):
    """
    Processes forward parcel tracking data.

    This function processes forward parcel tracking data, filters the data based on specified criteria, and returns the results.

    Parameters:
    lista_partposi (list): List of file paths to read.
    file_mask (str): Path to the mask file.
    name_mascara (str): Name of the mask variable.
    name_variable_lon (str): Name of the longitude variable in the mask file.
    name_variable_lat (str): Name of the latitude variable in the mask file.
    lat_f (numpy.ndarray): Array of latitude values.
    lon_f (numpy.ndarray): Array of longitude values.
    rank (int): Rank of the current MPI process.
    size (int): Total number of MPI processes.
    comm (MPI.Comm): MPI communicator.
    type_file (int): Type of file (1 for FLEXPART-WRF, 2 for FLEXPART-ERAI and FLEXPART-ERA5).
    x_left_lower_corner (float): X-coordinate of the lower left corner of the domain.
    y_left_lower_corner (float): Y-coordinate of the lower left corner of the domain.
    x_right_upper_corner (float): X-coordinate of the upper right corner of the domain.
    y_right_upper_corner (float): Y-coordinate of the upper right corner of the domain.
    model (str): Model type (e.g., "FLEXPART").
    value_mask (int): Value to search for in the mask array.
    key_gz (int): Whether to use gzip compression (1 for yes, 0 for no).
    path_output (str): Path to save the output files.
    use_vertical_layers (bool): Whether to use vertical layers.
    vertical_layers (list): List of vertical layers.
    filter_parcels_height (bool): Whether to filter parcels by height.
    filter_vertical_layers (list): List of vertical layers for filtering.
    limit_domain (int): Whether to limit the domain (1 for yes, 0 for no).

    Returns:
    tuple: A tuple containing the processed data and additional information.
    """
    name_file=lista_partposi[0]
    
    if rank==0:
        print ("Reading | " + model+" -> ",  name_file)
        name_txt_part=name_file.split("/")
        f=open(path_output+name_txt_part[-1].split("_")[-1].split(".")[0]+"/"+name_txt_part[-1].split("_")[-1].split(".")[0]+".txt", "a")
    
    if key_gz==1:
        desc_gz(name_file)
        part_post=read_binaryFile_fortran(name_file[:-3], type_file,x_left_lower_corner,y_left_lower_corner, x_right_upper_corner,y_right_upper_corner, limit_domain)
        cmd_rm= "rm -rf "+name_file[:-3]
        os.system(cmd_rm)
    else:
        part_post=read_binaryFile_fortran(name_file, type_file,x_left_lower_corner,y_left_lower_corner, x_right_upper_corner,y_right_upper_corner, limit_domain)

    lat_masked, lon_masked, mascara=load_mask_grid_NR(file_mask, name_mascara,name_variable_lon, name_variable_lat)
    submatrix=determine_id_binary_grid_NR_fortran(part_post, lat_masked.flatten(), lon_masked.flatten(), mascara.flatten(), value_mask)
    submatrix=submatrix[np.argsort(submatrix[:, 0])]

    if filter_parcels_height:
       submatrix, counter_part_height = Filter_by_Height(submatrix,submatrix,-1,filter_vertical_layers[0],filter_vertical_layers[1], len(submatrix[:, 0]))
    if rank==0:
        print ("Reading | " + model+" -> ",  lista_partposi[1])
    if key_gz==1:
        desc_gz(lista_partposi[1])
        part_post_i=read_binaryFile_fortran(lista_partposi[1][:-3], type_file, x_left_lower_corner,y_left_lower_corner, x_right_upper_corner,y_right_upper_corner, limit_domain)
        cmd_rm= "rm -rf "+lista_partposi[1][:-3]
        os.system(cmd_rm)
    else:
        part_post_i=read_binaryFile_fortran(lista_partposi[1], type_file, x_left_lower_corner,y_left_lower_corner, x_right_upper_corner,y_right_upper_corner, limit_domain)
    matrix_i=search_row_fortran(submatrix[:,0],part_post_i)

    dimX, dimY=matrix_i.shape
    tensor_t=np.ones((len(lista_partposi)-1,dimX ,dimY ))*(-999.9)
    tensor_t[0,:,:]=matrix_i
    tensor_org=np.ones((len(lista_partposi),dimX ,dimY ))*(-999.9)
    tensor_org[0,:,:]=submatrix
    tensor_org[1,:,:]=matrix_i
    n = len(lista_partposi)-2
    count = n // size
    remainder = n % size
    if rank < remainder:
        start = rank * (count + 1)
        stop = start + count + 1
    else:
        start = rank * count + remainder
        stop = start + count

    local_list=lista_partposi[start+2:stop+2]
    local_results = np.empty((len(local_list), dimX, dimY))
    local_results= read_proccesor(local_list, submatrix, rank,x_left_lower_corner,y_left_lower_corner,x_right_upper_corner,y_right_upper_corner, model, key_gz, type_file, limit_domain)

    if rank > 0:
        comm.Send(local_results, dest=0, tag=14)
    else:
        i_start=[]
        i_stop=[]
        for i in range(size):
            if i < remainder:
                i_start = np.append(i_start,i * (count + 1))
                i_stop = np.append(i_stop,i_start + count + 1)
            else:
                ii_start=i * count + remainder
                ii_stop=ii_start + count
                i_start = np.append(i_start,ii_start)
                i_stop = np.append(i_stop, ii_stop)

        final_results = np.copy(local_results)
        tensor_t[int(i_start[0])+1:int(i_stop[0])+1,:,:]= final_results
        tensor_org[int(i_start[0])+2:int(i_stop[0])+2,:,:]= final_results
        for i in range(1, size):
            if i < remainder:
                rank_size = count + 1
            else:
                rank_size = count
            tmp = np.empty((rank_size, final_results.shape[1],final_results.shape[2]), dtype=np.float64)
            comm.Recv(tmp, source=i, tag=14)
            tensor_t[int(i_start[i])+1:int(i_stop[i])+1,:,:]=tmp
            tensor_org[int(i_start[i]):int(i_stop[i]),:,:]=tmp
    comm.Bcast(tensor_t, root=0)
    comm.Bcast(tensor_org, root=0)

    matrix_result=np.ones((len(tensor_t[:,0,0])-1, len(submatrix[:,0]),4))*(-999.9)
    a2=np.arange(len(tensor_t[:,0,0])-1)
    for i in a2:
        matrix = Kdif_python(tensor_t[i,:,:5], tensor_t[i+1,:,:5],1.)
        matrix_result[i,:,2]=matrix[:,2]
        matrix_result[i,:,1]=matrix[:,1]
        matrix_result[i,:,0]=matrix[:,0]
        matrix_result[i,:,3]=matrix[:,3]

    if rank==0:
        print("")
        print ("   + Number of detected parcels within the target region => ", len(matrix_i[:,0]))
        f.write("%s %d\n"%("NumP: ",len(matrix_i[:,0])))
        if filter_parcels_height:
             print ("   + Number of parcels filter by height: Layer [" + str(filter_vertical_layers[0]) + " - " +str(filter_vertical_layers[1]) +" meters] => ", counter_part_height)
             f.write("%s %d\n"%("by height: ",counter_part_height))
    return matrix_result, tensor_t[0,:,0], tensor_t[0,:,3],tensor_org,lat_masked, lon_masked, mascara

def _vector_wvrt(lista_partposi ,file_mask, name_mascara,name_variable_lon, name_variable_lat,lat_f, lon_f,rank,size, comm, type_file,
                 x_left_lower_corner,y_left_lower_corner, x_right_upper_corner,y_right_upper_corner, model, method,threshold,filter_value,
                 value_mask, key_gz, path_output,use_vertical_layers, vertical_layers, filter_parcels_height, filter_vertical_layers,limit_domain, dates,
                 dtime, totaltime, folder):
    """
    Processes vector water vapor residence time (WVRT) data.

    This function processes vector water vapor residence time data, filters the data based on specified criteria, and returns the results.

    Parameters:
    lista_partposi (list): List of file paths to read.
    file_mask (str): Path to the mask file.
    name_mascara (str): Name of the mask variable.
    name_variable_lon (str): Name of the longitude variable in the mask file.
    name_variable_lat (str): Name of the latitude variable in the mask file.
    lat_f (numpy.ndarray): Array of latitude values.
    lon_f (numpy.ndarray): Array of longitude values.
    rank (int): Rank of the current MPI process.
    size (int): Total number of MPI processes.
    comm (MPI.Comm): MPI communicator.
    type_file (int): Type of file (1 for FLEXPART-WRF, 2 for FLEXPART-ERAI and FLEXPART-ERA5).
    x_left_lower_corner (float): X-coordinate of the lower left corner of the domain.
    y_left_lower_corner (float): Y-coordinate of the lower left corner of the domain.
    x_right_upper_corner (float): X-coordinate of the upper right corner of the domain.
    y_right_upper_corner (float): Y-coordinate of the upper right corner of the domain.
    model (str): Model type (e.g., "FLEXPART").
    method (int): Method used for processing.
    threshold (float): Threshold value for filtering.
    filter_value (int): Value used for filtering.
    value_mask (int): Value to search for in the mask array.
    key_gz (int): Whether to use gzip compression (1 for yes, 0 for no).
    path_output (str): Path to save the output files.
    use_vertical_layers (bool): Whether to use vertical layers.
    vertical_layers (list): List of vertical layers.
    filter_parcels_height (bool): Whether to filter parcels by height.
    filter_vertical_layers (list): List of vertical layers for filtering.
    limit_domain (int): Whether to limit the domain (1 for yes, 0 for no).
    dates (list): List of dates for the simulation.
    dtime (int): Time step interval in minutes.
    totaltime (int): Total time in minutes.
    folder (str): Folder name to save the output file within the output directory.

    Returns:
    tuple: A tuple containing the processed data and additional information.
    """

    name_file=lista_partposi[-1]
    if rank==0:
        print ("Reading | " + model+" -> ",  name_file)
        name_txt_part=name_file.split("/")
        f=open(path_output+name_txt_part[-1].split("_")[-1].split(".")[0]+"/"+name_txt_part[-1].split("_")[-1].split(".")[0]+".txt", "a")
    
    if key_gz==1:
        desc_gz(name_file)
        part_post=read_binaryFile_fortran(name_file[:-3], type_file, x_left_lower_corner,y_left_lower_corner,
               x_right_upper_corner,y_right_upper_corner, limit_domain)
        cmd_rm= "rm -rf "+name_file[:-3]
        os.system(cmd_rm)
    else:
        part_post=read_binaryFile_fortran(name_file, type_file, x_left_lower_corner,y_left_lower_corner,
               x_right_upper_corner,y_right_upper_corner, limit_domain)
        
    lat_masked, lon_masked, mascara=load_mask_grid_NR(file_mask, name_mascara,name_variable_lon, name_variable_lat)
    submatrix=determine_id_binary_grid_NR_fortran(part_post, lat_masked.flatten(), lon_masked.flatten(), mascara.flatten(), value_mask)
    submatrix=submatrix[np.argsort(submatrix[:, 0])]
    
    if rank==0:
        print ("Reading | " + model+" -> ",  lista_partposi[-2])

    if key_gz==1:
        desc_gz(lista_partposi[-2])
        part_post_i=read_binaryFile_fortran(lista_partposi[-2][:-3], type_file, x_left_lower_corner,y_left_lower_corner,x_right_upper_corner,y_right_upper_corner, limit_domain)
        cmd_rm= "rm -rf "+lista_partposi[-2][:-3]
        os.system(cmd_rm)
    else:
        part_post_i=read_binaryFile_fortran(lista_partposi[-2], type_file, x_left_lower_corner,y_left_lower_corner,x_right_upper_corner,y_right_upper_corner, limit_domain)
    
    matrix_i=search_row_fortran(submatrix[:,0],part_post_i)

    dimX, dimY=matrix_i.shape

    tensor_t=np.ones((len(lista_partposi)-1,dimX ,dimY ))*(-999.9)
    tensor_t[-1,:,:]=matrix_i

    tensor_org=np.ones((len(lista_partposi),dimX ,dimY ))*(-999.9)
    tensor_org[-1,:,:]=submatrix
    tensor_org[-2,:,:]=matrix_i
    
    n = len(lista_partposi)-2
    count = n // size
    remainder = n % size

    if rank < remainder:
        start = rank * (count + 1)
        stop = start + count + 1
    else:
        start = rank * count + remainder
        stop = start + count

    local_list=lista_partposi[start:stop]
    local_results = np.empty((len(local_list), dimX, dimY))
    local_results = read_proccesor(local_list, submatrix, rank,x_left_lower_corner,y_left_lower_corner,
               x_right_upper_corner,y_right_upper_corner, model, key_gz, type_file, limit_domain)

    if rank > 0:
        comm.Send(local_results, dest=0, tag=14)
    else:
        i_start=[]
        i_stop=[]
        for i in range(size):
            if i < remainder:
                i_start = np.append(i_start,i * (count + 1))
                i_stop = np.append(i_stop,i_start + count + 1)
            else:
                ii_start=i * count + remainder
                ii_stop=ii_start + count
                i_start = np.append(i_start,ii_start)
                i_stop = np.append(i_stop, ii_stop)
        final_results = np.copy(local_results)
        tensor_t[int(i_start[0]):int(i_stop[0]),:,:]= final_results
        tensor_org[int(i_start[0]):int(i_stop[0]),:,:]= final_results
        for i in range(1, size):
            if i < remainder:
                rank_size = count + 1
            else:
                rank_size = count
            tmp = np.empty((rank_size, final_results.shape[1],final_results.shape[2]), dtype=np.float64)
            comm.Recv(tmp, source=i, tag=14)
            tensor_t[int(i_start[i]):int(i_stop[i]),:,:]=tmp
            tensor_org[int(i_start[i]):int(i_stop[i]),:,:]=tmp
    
    comm.Bcast(tensor_t, root=0)
    comm.Bcast(tensor_org, root=0)
        
    matrix_result=np.ones((len(tensor_t[:,0,0])-1, len(submatrix[:,0]),4))*(-999.9)
    a2=np.arange(len(tensor_t[:,0,0])-1)
     
    for i in a2[::-1]:
        matrix = Kdif_python (tensor_t[i,:,:5], tensor_t[i+1,:,:5],-1.)
        matrix_result[i,:,2]=matrix[:,2]
        matrix_result[i,:,1]=matrix[:,1]
        matrix_result[i,:,0]=matrix[:,0]
        matrix_result[i,:,3]=matrix[:,3]

    submatrix_wvrt = matrix_result[:, matrix_result[0, :, 0] != -999.9, :]
    wvrt_mean = compute_residence_time_and_save(submatrix_wvrt[:, :, 2], path_output+folder+"/", dates, dtime, totaltime, folder, rank)

    if rank==0:
        print("")
        print("   + Number of detected parcels within the target region => ", len(matrix_i[:,0]))
        print("   + Mean water vapor residence time => %.2f days" % np.abs(wvrt_mean))
        f.write("%s %d %s %.2f\n" % ("NumP: ", len(matrix_i[:, 0]), "Mean_WVRT: ", np.abs(wvrt_mean)))

def _only_partposit_particles(lista_partposi ,file_mask, name_mascara,name_variable_lon, name_variable_lat,lat_f, lon_f,rank,size, comm, type_file,
                 x_left_lower_corner,y_left_lower_corner, x_right_upper_corner,y_right_upper_corner, model, method,threshold,filter_value,
                 value_mask, key_gz, path_output,use_vertical_layers, vertical_layers, filter_parcels_height, filter_vertical_layers,limit_domain): 
    """
    Processes only particle positions data.

    This function processes only particle positions data, filters the data based on specified criteria, and returns the results.

    Parameters:
    lista_partposi (list): List of file paths to read.
    file_mask (str): Path to the mask file.
    name_mascara (str): Name of the mask variable.
    name_variable_lon (str): Name of the longitude variable in the mask file.
    name_variable_lat (str): Name of the latitude variable in the mask file.
    lat_f (numpy.ndarray): Array of latitude values.
    lon_f (numpy.ndarray): Array of longitude values.
    rank (int): Rank of the current MPI process.
    size (int): Total number of MPI processes.
    comm (MPI.Comm): MPI communicator.
    type_file (int): Type of file (1 for FLEXPART-WRF, 2 for FLEXPART-ERAI and FLEXPART-ERA5).
    x_left_lower_corner (float): X-coordinate of the lower left corner of the domain.
    y_left_lower_corner (float): Y-coordinate of the lower left corner of the domain.
    x_right_upper_corner (float): X-coordinate of the upper right corner of the domain.
    y_right_upper_corner (float): Y-coordinate of the upper right corner of the domain.
    model (str): Model type (e.g., "FLEXPART").
    method (int): Method used for processing.
    threshold (float): Threshold value for filtering.
    filter_value (int): Value used for filtering.
    value_mask (int): Value to search for in the mask array.
    key_gz (int): Whether to use gzip compression (1 for yes, 0 for no).
    path_output (str): Path to save the output files.
    use_vertical_layers (bool): Whether to use vertical layers.
    vertical_layers (list): List of vertical layers.
    filter_parcels_height (bool): Whether to filter parcels by height.
    filter_vertical_layers (list): List of vertical layers for filtering.
    limit_domain (int): Whether to limit the domain (1 for yes, 0 for no).

    Returns:
    tuple: A tuple containing the processed data and additional information.
    """

    name_file=lista_partposi[-1]
    if rank==0:
        print ("Reading | " + model+" -> ",  name_file)
        name_txt_part=name_file.split("/")
        f=open(path_output+name_txt_part[-1].split("_")[-1].split(".")[0]+"/"+name_txt_part[-1].split("_")[-1].split(".")[0]+".txt", "a")
    
    if key_gz==1:
        desc_gz(name_file)
        part_post=read_binaryFile_fortran(name_file[:-3], type_file, x_left_lower_corner,y_left_lower_corner,
               x_right_upper_corner,y_right_upper_corner, limit_domain)
        cmd_rm= "rm -rf "+name_file[:-3]
        os.system(cmd_rm)
    else:
        part_post=read_binaryFile_fortran(name_file, type_file, x_left_lower_corner,y_left_lower_corner,
               x_right_upper_corner,y_right_upper_corner, limit_domain)
        
    lat_masked, lon_masked, mascara=load_mask_grid_NR(file_mask, name_mascara,name_variable_lon, name_variable_lat)
    submatrix=determine_id_binary_grid_NR_fortran(part_post, lat_masked.flatten(), lon_masked.flatten(), mascara.flatten(), value_mask)
    submatrix=submatrix[np.argsort(submatrix[:, 0])]
        
    if rank==0:
        print ("Reading | " + model+" -> ",  lista_partposi[-2])

    if key_gz==1:
        desc_gz(lista_partposi[-2])
        part_post_i=read_binaryFile_fortran(lista_partposi[-2][:-3], type_file, x_left_lower_corner,y_left_lower_corner,x_right_upper_corner,y_right_upper_corner, limit_domain)
        cmd_rm= "rm -rf "+lista_partposi[-2][:-3]
        os.system(cmd_rm)
    else:
        part_post_i=read_binaryFile_fortran(lista_partposi[-2], type_file, x_left_lower_corner,y_left_lower_corner,x_right_upper_corner,y_right_upper_corner, limit_domain)
    
    matrix_i=search_row_fortran(submatrix[:,0],part_post_i)
    
    dimX, dimY=matrix_i.shape

    tensor_t=np.ones((len(lista_partposi)-1,dimX ,dimY ))*(-999.9)
    tensor_t[-1,:,:]=matrix_i

    tensor_org=np.ones((len(lista_partposi),dimX ,dimY ))*(-999.9)
    tensor_org[-1,:,:]=submatrix
    tensor_org[-2,:,:]=matrix_i
    
    n = len(lista_partposi)-2
    count = n // size
    remainder = n % size

    if rank < remainder:
        start = rank * (count + 1)
        stop = start + count + 1
    else:
        start = rank * count + remainder
        stop = start + count

    local_list=lista_partposi[start:stop]
    local_results = np.empty((len(local_list), dimX, dimY))
    local_results = read_proccesor(local_list, submatrix, rank,x_left_lower_corner,y_left_lower_corner,
               x_right_upper_corner,y_right_upper_corner, model, key_gz, type_file, limit_domain)

    if rank > 0:
        comm.Send(local_results, dest=0, tag=14)
    else:
        i_start=[]
        i_stop=[]
        for i in range(size):
            if i < remainder:
                i_start = np.append(i_start,i * (count + 1))
                i_stop = np.append(i_stop,i_start + count + 1)
            else:
                ii_start=i * count + remainder
                ii_stop=ii_start + count
                i_start = np.append(i_start,ii_start)
                i_stop = np.append(i_stop, ii_stop)
        final_results = np.copy(local_results)
        tensor_t[int(i_start[0]):int(i_stop[0]),:,:]= final_results
        tensor_org[int(i_start[0]):int(i_stop[0]),:,:]= final_results
        for i in range(1, size):
            if i < remainder:
                rank_size = count + 1
            else:
                rank_size = count
            tmp = np.empty((rank_size, final_results.shape[1],final_results.shape[2]), dtype=np.float64)
            comm.Recv(tmp, source=i, tag=14)
            tensor_t[int(i_start[i]):int(i_stop[i]),:,:]=tmp
            tensor_org[int(i_start[i]):int(i_stop[i]),:,:]=tmp
    
    comm.Bcast(tensor_t, root=0)
    comm.Bcast(tensor_org, root=0)
            
    if rank==0:
        print("")
        print ("   + Number of detected parcels within the target region => ", len(matrix_i[:,0]))
        print ("   + Saved the partposit file in netcdf of variables corresponding to the particles within the target region")
        f.write("%s %d\n"%("NumP: ",len(matrix_i[:,0])))
    return tensor_org

def time_calc_day(init_time,day_diff):
    """
    Calculates a new date based on the initial date and a day difference.

    This function calculates a new date by adding the specified day difference to the initial date.

    Parameters:
    init_time (str): The initial date in the format "YYYY-MM-DD HH:MM:SS".
    day_diff (int): The day difference to add.

    Returns:
    datetime: The calculated date.
    """

    formatted_time = datetime.strptime(init_time, "%Y-%m-%d %H:%M:%S")
    calculated_time=formatted_time+timedelta(days=int(day_diff))
    return calculated_time

def convert_date_to_ordinal(year, month, day, hour, minute, second):
    """
    Converts a date to an ordinal number.

    This function converts the specified date and time components to an ordinal number based on the NetCDF convention.

    Parameters:
    year (int): The year component of the date.
    month (int): The month component of the date.
    day (int): The day component of the date.
    hour (int): The hour component of the date.
    minute (int): The minute component of the date.
    second (int): The second component of the date.

    Returns:
    float: The ordinal number representing the date.
    """

    date=datetime(year, month, day, hour, minute, second)
    date_=netCDF4.date2num(date, "days since 1900-01-01", "gregorian")
    return date_

def decompose_date(value):
    """
    Decomposes an ordinal date value into its components.

    This function decomposes an ordinal date value into its year, month, day, hour, minute, and second components.

    Parameters:
    value (float): The ordinal date value to decompose.

    Returns:
    tuple: A tuple containing the year, month, day, hour, minute, and second components.
    """

    str_value = str(int(value))
    year = int(str_value[:4])
    month = int(str_value[4:6])
    day = int(str_value[6:8])
    hour = int(str_value[8:10])
    minute = int(str_value[10:12])
    second = int(str_value[12:14])
    return year, month, day, hour, minute, second

class InputNotInRangeError(Exception):
    """
    Exception raised for errors in the input parameters.

    This exception is raised when an input parameter is not within the expected range.

    Attributes:
    message (str): Explanation of the error.
    """

    def __init__(self, message):
        self.message = message
        super().__init__(self.message)

def to_check_params(paso, type_file, numPdX, numPdY, method, resolution, file_mask):
    """
    Checks the validity of input parameters.

    This function checks if the input parameters are within the expected range and raises an exception if they are not.

    Parameters:
    paso (int): Step value indicating the direction of the simulation (-1 for backward, 1 for forward, -2 for WVRT, -3 for saving variables for particles in the target region).
    type_file (int): Type of file (1 for FLEXPART-WRF, 2 for FLEXPART-ERAI and FLEXPART-ERA5).
    numPdX (int): Number of grid points in the X direction.
    numPdY (int): Number of grid points in the Y direction.
    method (int): Method used for processing (1 for Stohl and James, 2 for Sodemann).
    resolution (float): The resolution of the grid cells.
    file_mask (str): Path to the mask file.

    Returns:
    None

    Raises:
    InputNotInRangeError: If any input parameter is not within the expected range.
    """

    valid_paso_values = [-1, 1, -2, -3]
    valid_type_file_values = [1, 2]
    valid_method_values = [1, 2]

    if paso not in valid_paso_values:
        raise InputNotInRangeError("Only calculations are allowed for moisture sources (mode = -1, backward), sinks (mode = 1, forward), water vapor residence time (mode = -2, wvrt) or saving variables for particles in the target region (mode =-3)")

    if type_file not in valid_type_file_values:
        raise InputNotInRangeError("Only FLEXPART-WRF (type_file = 1) and FLEXPART model (type_file = 2) are allowed")

    if method not in valid_method_values:
        raise InputNotInRangeError("Only Stohl and James (method = 1) and Sodemann methodology (method = 2) are allowed")

    if numPdX <= 0 or numPdY <= 0:
        raise InputNotInRangeError("NumPdX and numPdY must be greater than zero")

    if resolution <= 0:
        raise InputNotInRangeError("Resolution must be greater than zero")

    my_file = Path(file_mask)
    if not my_file.is_file():
        raise InputNotInRangeError("File mask not found")

def function_proof(lat, lon):
    """
    Checks if latitude and longitude points are within the valid range.

    This function checks if the latitude and longitude points are within the valid range and raises an exception if they are not.

    Parameters:
    lat (numpy.ndarray): Array of latitude values.
    lon (numpy.ndarray): Array of longitude values.

    Returns:
    None

    Raises:
    InputNotInRangeError: If any latitude or longitude point is not within the valid range.
    """

    lat_=lat.flatten()
    lon_=lon.flatten()
    for i in range(len(lat_)):
        if -90 <= lat_[i] <= 90:
            pass
        else:
            raise InputNotInRangeError("A latitude point is not in the range -90 and 90")
        if -180 <= lon_[i]<= 180 :
            pass
        else:
            raise InputNotInRangeError("A longitude point is not in the range -180 and 180")

def desc_gz(name_file):
    """
    Decompresses a gzip file.

    This function decompresses the specified gzip file and saves the decompressed content to a new file.

    Parameters:
    name_file (str): The name of the gzip file to decompress.

    Returns:
    None
    """

    with gzip.open(name_file, 'rb') as f_in:
        with open(name_file[:-3], 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

def is_binary(file_path):
    """
    Checks if a file is binary by looking for non-text characters.
    """
    try:
        with open(file_path, "rb") as f:
            chunk = f.read(1024)  
            if b"\x00" in chunk: 
                return True
    except Exception:
        return False
    return False 

def verify_binary_files(file_list):
    """
    Verifies if a list of files meet the following conditions:
    1. The filename (after the last '/') starts with 'partposit'.
    2. The file is binary (contains non-text characters).
    3. The file size is within ±10% of the first valid file.

    Parameters:
    - file_list (list): List of file paths.

    Returns:
    - valid_files (list): List of files that meet the conditions.
    
    Prints:
    - Errors for files that do not meet the conditions.
    """
    valid_files = []
    errors = []
    reference_size = None  # To store the size of the first valid file

    for file_path in file_list:
        # Check if the file exists
        if not os.path.exists(file_path):
            errors.append(f"Error: The file '{file_path}' does not exist.")
            break

        # Extract only the filename (after the last '/')
        filename = os.path.basename(file_path)

        # Check if the filename starts with 'partposit'
        if not filename.startswith("partposit"):
            errors.append(f"Error: The filename '{filename}' does not start with 'partposit'.")
            break

        # Check if the file is binary
        if not is_binary(file_path):
            errors.append(f"Error: The file '{filename}' is not a valid binary file.")
            break

        # Get file size
        file_size = os.path.getsize(file_path)

        # Set reference size with the first valid file
        if reference_size is None:
            reference_size = file_size

        # Check if the file size is within ±10% of the reference
        lower_bound = reference_size * 0.9
        upper_bound = reference_size * 1.1
        if not (lower_bound <= file_size <= upper_bound):
            errors.append(f"Error: The file '{filename}' has an unusual size ({file_size} bytes).")
            break

        # If all conditions are met, add the file to the valid list
        valid_files.append(file_path)

    # Get the rank of the current process
    rank = MPI.COMM_WORLD.Get_rank()

    # If there were errors, stop the program and print the message (only rank == 0)
    if errors:
        if rank == 0:
            for error in errors:
                print(error)
            print ('                                                            ')
            print ('------------------------------------------------------------')
            print("ERROR: Please check your input data!")
            print ('                                                            ')  
        sys.exit(1)  # Stop the program with an error code
    # If all files passed and rank == 0, print success message
    if rank == 0:
        print ('                                                            ')
        print ('------------------------------------------------------------')
        print("!!!! IMPORTANT: All files passed the checks and are valid!!!!")
        print ('------------------------------------------------------------')
        print ('                                                            ')
    return valid_files

def TROVA_LOGO():
    """
    Prints the TROVA logo.

    This function prints the TROVA logo to the terminal.
    """
    print(" *                        _____ __    ____                                               *")
    print(" *                          |  |  |  /    \ \        //\                                 *")
    print(" *                          |  |__| /      \ \      //__\                                *")
    print(" *                          |  |  \ \      /  \    //    \                               *")
    print(" *                          |  |   \ \____/    \__//      \                              *")
    print(" *                                                                                       *")

def main_process(path, paso, comm, size, rank, resolution, numPdX, numPdY, dtime, totaltime, year, month,
         day, hour, minn, time, path_output,file_mask, name_mascara,name_variable_lon, name_variable_lat,x_lower_left,y_lower_left, 
         type_file,masa,numP, x_left_lower_corner,y_left_lower_corner, x_right_upper_corner,y_right_upper_corner, model, method,
         threshold, filter_value, output_txt, output_npy,output_nc, value_mask, key_gz, save_position_part, use_vertical_layers, 
         vertical_layers,save_position_dqdt, filter_parcels_height,filter_vertical_layers, plotting_parcels_t0, 
         plotting_parcels_tracks_on_map, plotting_3Dparcels_tracks, maps_limits, noleap, limit_domain, method_wvrt, 
         plotting_moisture_sink_source, limit_plot):
    """
    Main processing function for TROVA.

    This function performs the main processing tasks for TROVA, including reading input files, processing data,
    and generating output files and plots.

    Parameters:
    path (str): Path to the input data.
    paso (int): Step value indicating the direction of the simulation (-1 for backward, 1 for forward).
    comm (MPI.Comm): MPI communicator.
    size (int): Total number of MPI processes.
    rank (int): Rank of the current MPI process.
    resolution (float): The resolution of the grid cells.
    numPdX (int): Number of grid points in the X direction.
    numPdY (int): Number of grid points in the Y direction.
    dtime (int): Time step in minutes.
    totaltime (int): Total simulation time in minutes.
    year (str): Initial year of the simulation.
    month (str): Initial month of the simulation.
    day (str): Initial day of the simulation.
    hour (str): Initial hour of the simulation.
    minn (str): Initial minute of the simulation.
    time (float): Start time of the simulation.
    path_output (str): Path to save the output files.
    file_mask (str): Path to the mask file.
    name_mascara (str): Name of the mask variable.
    name_variable_lon (str): Name of the longitude variable in the mask file.
    name_variable_lat (str): Name of the latitude variable in the mask file.
    x_lower_left (float): X-coordinate of the lower left corner of the domain.
    y_lower_left (float): Y-coordinate of the lower left corner of the domain.
    type_file (int): Type of file (1 for FLEXPART-WRF, 2 for FLEXPART-ERAI and FLEXPART-ERA5).
    masa (float): Mass of the parcels.
    numP (int): Number of parcels.
    x_left_lower_corner (float): X-coordinate of the lower left corner of the domain.
    y_left_lower_corner (float): Y-coordinate of the lower left corner of the domain.
    x_right_upper_corner (float): X-coordinate of the upper right corner of the domain.
    y_right_upper_corner (float): Y-coordinate of the upper right corner of the domain.
    model (str): Model type (e.g., "FLEXPART").
    method (int): Method used for processing.
    threshold (float): Threshold value for filtering.
    filter_value (int): Value used for filtering.
    output_txt (int): Whether to output results in TXT format (1 for yes, 0 for no).
    output_npy (int): Whether to output results in NPY format (1 for yes, 0 for no).
    output_nc (int): Whether to output results in NetCDF format (1 for yes, 0 for no).
    value_mask (int): Value to search for in the mask array.
    key_gz (int): Whether to use gzip compression (1 for yes, 0 for no).
    save_position_part (bool): Whether to save parcel positions at each time step.
    use_vertical_layers (bool): Whether to use vertical layers.
    vertical_layers (list): List of vertical layers.
    save_position_dqdt (bool): Whether to save dq/dt at each time step.
    filter_parcels_height (bool): Whether to filter parcels by height.
    filter_vertical_layers (list): List of vertical layers for filtering.
    plotting_parcels_t0 (bool): Whether to plot identified parcels within the target region at time t0.
    plotting_parcels_tracks_on_map (bool): Whether to plot identified parcels' trajectories on a map.
    plotting_3Dparcels_tracks (bool): Whether to plot 3D parcels' trajectories.
    maps_limits (list): List containing the map limits [latmin, lonmin, latmax, lonmax, center, dlat, dlon].
    noleap (bool): Whether to exclude leap years.
    limit_domain (int): Whether to limit the domain (1 for yes, 0 for no).
    method_wvrt (int): Method used for calculating water vapor residence time.
    plotting_moisture_sink_source (bool): Whether to plot moisture sink and source patterns.
    limit_plot (list): List containing the plot limits [latmin, lonmin, latmax, lonmax].

    Returns:
    None
    """

    create_directory(path_output)

    filesperday=int(1440/dtime)

    t=np.arange(0,int((totaltime/dtime))+filesperday, filesperday)

    lat,lon= grid_point (resolution, numPdX, numPdY,x_lower_left,y_lower_left)
    area=calc_A(resolution, lat, lon)
    dimX, dimY=lat.shape
    lat_plot,lon_plot= grid_plot_final(lat, lon)
    function_proof(lat_plot, lon_plot)
    np.savetxt(path_output+"lat_plot.txt", lat_plot)
    np.savetxt(path_output+"lon_plot.txt", lon_plot)
    folder=year+month+day+hour+minn+"00"
    create_directory(path_output+folder)
    fecha= year+"-"+month+"-"+day+" "+hour+":"+minn+":"+"00"
    lista_partposi, listdates=generate_file(paso, dtime, totaltime, fecha, path, key_gz, noleap)
  
    if rank==0:
       if (len(lista_partposi)-2)%size!=0:
           print("TROVA ERROR: The number of processors must exactly divide the number of partposit files to process (totaltime/dtime).\n Based on your configuration file, the recommended number of processors is " + str(int(totaltime/dtime))) 
           raise SystemExit("Bye :)")
    elif (len(lista_partposi)-2)%size!=0:
      raise SystemExit()
    
    verify_binary_files(lista_partposi)
    
    if paso==1:
        matrix_result, id_part, q_ini, partpos, lat_masked, lon_masked, mascara =_forward_dq(lista_partposi, file_mask, name_mascara,name_variable_lon, name_variable_lat,lat, lon,rank,size,comm, type_file,x_left_lower_corner,y_left_lower_corner, 
                                                                            x_right_upper_corner,y_right_upper_corner, model, value_mask, key_gz, path_output,use_vertical_layers, vertical_layers,filter_parcels_height,filter_vertical_layers, limit_domain)
        matrix_save=np.empty((matrix_result.shape[0],matrix_result.shape[1],6))
        matrix_save[:,:,:-2]=matrix_result
        matrix_save[:,:,-2]=id_part
        matrix_save[:,:,-1]=q_ini

        if rank == 0:
            if save_position_dqdt:
                ordinal_dates_dqdt = list(map(lambda date: convert_date_to_ordinal(*decompose_date(date)), listdates[2:]))
                write_nc(ordinal_dates_dqdt, matrix_save, "dqdt", filename=path_output+folder+"/"+folder+"_dqdt_forw")
            if  save_position_part:
                ordinal_dates_part = list(map(lambda date: convert_date_to_ordinal(*decompose_date(date)), listdates))
                write_nc(ordinal_dates_part, partpos, "partpos", filename=path_output+folder+"/"+folder+"_parposit_forw")

    if paso == -1:
        matrix_result, id_part, matrix_result_por,q_ini,partpos, lat_masked, lon_masked, mascara =_backward_dq(lista_partposi, file_mask, name_mascara,name_variable_lon, name_variable_lat,lat, lon,rank,size,comm, type_file, x_left_lower_corner, y_left_lower_corner, 
                                                                                                            x_right_upper_corner, y_right_upper_corner, model, method,threshold,filter_value, value_mask, key_gz, path_output,use_vertical_layers, vertical_layers,filter_parcels_height,
                                                                                                            filter_vertical_layers, limit_domain, listdates[-1:],dtime, totaltime, folder, method_wvrt)

        matrix_save=np.empty((matrix_result.shape[0],matrix_result.shape[1],6))
        matrix_save[:,:,:-2]=matrix_result
        matrix_save[:,:,-2]=id_part
        matrix_save[:,:,-1]=q_ini
         
        if rank==0:
            if save_position_dqdt:
               ordinal_dates_dqdt = list(map(lambda date: convert_date_to_ordinal(*decompose_date(date)), listdates[1:-1]))
               write_nc(ordinal_dates_dqdt, matrix_save, "dqdt", filename=path_output+folder+"/"+folder+"_dqdt_back")
            if save_position_part:
               ordinal_dates_part = list(map(lambda date: convert_date_to_ordinal(*decompose_date(date)), listdates))
               write_nc(ordinal_dates_part, partpos, "partpos",filename=path_output+folder+"/"+folder+"_parposit_back")
    
    if paso == -2:
        _vector_wvrt(lista_partposi, file_mask, name_mascara,name_variable_lon, name_variable_lat,lat, lon,rank,size,comm, type_file, x_left_lower_corner, y_left_lower_corner, 
                                                                                                            x_right_upper_corner, y_right_upper_corner, model, method,threshold,filter_value, value_mask, key_gz, path_output,use_vertical_layers, vertical_layers,filter_parcels_height,
                                                                                                            filter_vertical_layers, limit_domain, listdates[-1:],dtime, totaltime, folder)
    
    if paso == -3:
        partpos =_only_partposit_particles(lista_partposi, file_mask, name_mascara,name_variable_lon, name_variable_lat,lat, lon,rank,size,comm, type_file, x_left_lower_corner, y_left_lower_corner, 
                                                                                                            x_right_upper_corner, y_right_upper_corner, model, method,threshold,filter_value, value_mask, key_gz, path_output,use_vertical_layers, vertical_layers,filter_parcels_height,
                                                                                                            filter_vertical_layers, limit_domain)         
        if rank==0:
            if save_position_part:
               ordinal_dates_part = list(map(lambda date: convert_date_to_ordinal(*decompose_date(date)), listdates))
               write_nc(ordinal_dates_part, partpos, "partpos",filename=path_output+folder+"/"+folder+"_parposit_back")
        
    if rank == 0:
        if paso in [-1, 1]:
            cant_plazo=totaltime/dtime
            ndf=np.arange(1,len(t),1)
            ndb=np.arange(len(t)-1,0,-1)
            density=masa/numP
            array_day=np.empty((len(t)-1,dimX-1, dimY-1))
            array_day_por=np.empty((len(t)-1,dimX-1, dimY-1))
            date_save=[]
            date_save.append(convert_date_to_ordinal(int(year), int(month), int(day), int(hour), int(minn), 0))

            if use_vertical_layers:
                array_day_layers=np.empty((len(t)-1,len(vertical_layers)-1,dimX-1, dimY-1))
                array_day_por_layers=np.empty((len(t)-1,len(vertical_layers)-1,dimX-1, dimY-1))
            else:
                vertical_layers=np.arange(0,2,1)
                array_day_layers=np.empty((len(t)-1,len(vertical_layers)-1,dimX-1, dimY-1))
                array_day_layers[:]=-999.9
                array_day_por_layers=np.empty((len(t)-1,len(vertical_layers)-1,dimX-1, dimY-1))
                array_day_por_layers[:]=-999.9

            for ii in range(len(t)-1):
                final_results=np.array(K_dq(matrix_result[t[ii]:t[ii+1],:,:-1],lon,lat,numPdY,numPdX,len(matrix_result[t[ii]:t[ii+1],0,0]),len(matrix_result[0,:,0])),dtype=np.float64)
                E_P=final_results*(density)/area

                if use_vertical_layers:
                    for ilayer in range(0,len(vertical_layers)-1):
                        temp_res=np.array(K_dq_layers(matrix_result[t[ii]:t[ii+1],:,:],vertical_layers[ilayer],vertical_layers[ilayer+1], lon,lat,numPdY,numPdX,len(matrix_result[t[ii]:t[ii+1],0,0]), len(matrix_result[0,:,0])),dtype=np.float64)
                        vars()["EP"+str(int(vertical_layers[ilayer]))+"_"+str(int(vertical_layers[ilayer+1]))]=temp_res*density/area
                
                if method==2:
                    POR=np.array(K_dq_POR(matrix_result_por[t[ii]:t[ii+1],:,:],lon,lat,numPdY,numPdX,len(matrix_result_por[t[ii]:t[ii+1],0,0]),len(matrix_result_por[0,:,0])),dtype=np.float64)

                if paso ==-1:
                    array_day[int (ndb[ii]-1), :,:]=E_P
                    date_back=time_calc_day(str(year)+"-"+str(month)+"-"+str(day)+" "+str(hour)+":"+str(minn)+":00", ndb[ii]*(-1))
                    date_save.append(int(convert_date_to_ordinal(int(date_back.year), int(date_back.month), int(date_back.day), int(date_back.hour), int(date_back.minute), 0)))
                
                    if use_vertical_layers:
                        for ilayer in range(0,len(vertical_layers)-1):
                            array_day_layers[int (ndb[ii]-1), ilayer, :,:]=vars()["EP"+str(int(vertical_layers[ilayer]))+"_"+str(int(vertical_layers[ilayer+1]))]

                    if output_txt!=0:
                        np.savetxt(path_output+folder+"/day_"+str(ndb[ii])+"_"+folder+".txt", E_P)
                        if use_vertical_layers:
                            for ilayer in range(0,len(vertical_layers)-1):
                                np.savetxt(path_output+folder+"/day_"+str(ndb[ii])+"_"+str(int(vertical_layers[ilayer]))+"_"+str(int(vertical_layers[ilayer+1]))+"_"+folder+".txt", vars()["POR"+str(int(vertical_layers[ilayer]))+"_"+str(int(vertical_layers[ilayer+1]))]) 

                    if method==2:
                        array_day_por[int (ndb[ii]-1), :,:]=POR
                        if output_txt!=0:
                            np.savetxt(path_output+folder+"/day_POR_"+str(ndb[ii])+"_"+folder+".txt", POR)

                if paso ==1:
                    array_day[int(ndf[ii]-1), :,:]=E_P
                    date_forw=time_calc_day(year+"-"+month+"-"+day+" "+hour+":"+minn+":00", ndf[ii])
                    date_save.append(int(convert_date_to_ordinal(int(date_forw.year), int(date_forw.month), int(date_forw.day), int(date_forw.hour), int(date_forw.minute), 0)))
                    if use_vertical_layers:
                        for ilayer in range(0,len(vertical_layers)-1):
                            array_day_layers[int(ndf[ii]-1), ilayer,:,:]=vars()["EP"+str(int(vertical_layers[ilayer]))+"_"+str(int(vertical_layers[ilayer+1]))]

                    if output_txt!=0:
                        np.savetxt(path_output+folder+"/day_"+str(ndf[ii])+"_"+folder+".txt", E_P)
                        if use_vertical_layers:
                            for ilayer in range(0,len(vertical_layers)-1):
                                np.savetxt(path_output+folder+"/day_"+str(ndb[ii])+"_"+str(int(vertical_layers[ilayer]))+"_"+str(int(vertical_layers[ilayer+1]))+"_"+folder+".txt", vars()["EP"+str(int(vertical_layers[ilayer]))+"_"+str(int(vertical_layers[ilayer+1]))])

            if paso ==-1:
                if output_nc!=0:
                    function(lat_plot[:,0], lon_plot[0,:], array_day, array_day_layers, use_vertical_layers, vertical_layers, method, array_day_por, "back_"+folder, path_output+folder+"/", "E_P", "mm/day", date_save[1:])
                if output_npy!=0:
                    np.save(path_output+folder+"/"+"back_"+folder+".npy", array_day)
                    if use_vertical_layers:
                        np.save(path_output+folder+"/"+"back_layers_"+folder+".npy", array_day_layers)
                if method==2:
                    if output_npy!=0:
                        np.save(path_output+folder+"/"+"POR_back_"+folder+".npy", array_day_por)
                if plotting_moisture_sink_source:
                    plot_moisture_sink_source(lon_plot[0,:], lat_plot[:,0], np.sum(array_day, axis=0), paso, path_output, folder, limit_plot)
                
            if paso == 1:
                if output_nc!=0:
                    array_day_por=np.empty_like(array_day)
                    array_day_por[:]=-999.9
                    function(lat_plot[:,0], lon_plot[0,:], array_day, array_day_layers, use_vertical_layers, vertical_layers,  method, array_day_por, "forw_"+folder, path_output+folder+"/", "E_P", "mm/day", date_save[:-1])
                if output_npy!=0:
                    np.save(path_output+folder+"/"+"forw_"+folder+".npy", array_day)
                    if use_vertical_layers:
                        np.save(path_output+folder+"/"+"forw_layers_"+folder+".npy", array_day_layers)
                if plotting_moisture_sink_source:
                    plot_moisture_sink_source(lon_plot[0,:], lat_plot[:,0], np.sum(array_day, axis=0), paso, path_output, folder, limit_plot)
                
            if rank==0:
                print()
                if plotting_3Dparcels_tracks or plotting_parcels_tracks_on_map or plotting_parcels_t0:
                    print ('------------------------------------------------------------------------------------------')
                    print ('                                   TROVA PLOTTING                                         ')
                    print ('------------------------------------------------------------------------------------------')
                    if  plotting_3Dparcels_tracks:
                        plotting_tracks_3d(partpos, path_output+folder+"/"+folder+"_3Dtracks.png")
                    if plotting_parcels_tracks_on_map:
                        ploting_parcels_tracks_map(partpos, maps_limits, paso, lat_masked, lon_masked, mascara, value_mask, path_output+folder+"/"+folder+"_parcels_tracks.png" )
                    if plotting_parcels_t0:
                        plotting_parcels_within_target_region(partpos, maps_limits, paso, lat_masked, lon_masked, mascara, value_mask, path_output+folder+"/"+folder+"_parcels_at_time_t0.png" )
                    print()


def TROVA_main(input_file):
    """
    Main function for TROVA.

    This function reads the input configuration file, initializes parameters, and starts the main processing function.

    Parameters:
    input_file (str): Path to the input configuration file.

    Returns:
    None
    """
    
    loader = importlib.machinery.SourceFileLoader("", input_file)
    content = loader.load_module()
    year_init = check_paths(content, "year_case_init")
    path = check_paths(content, "path_data")
    path_output = check_paths(content, "path_output")
    mode = check_paths(content, "mode")
    masa = check_paths(content, "mass")
    numP = check_paths(content, "numP")
    resolution = check_paths(content,"resolution")
    numPdX = check_paths(content,"numPdX")
    numPdY = check_paths(content,"numPdY")
    x_lower_left = check_paths(content,"x_lower_left")
    y_lower_left = check_paths(content,"y_lower_left")
    dtime = check_paths(content,"dtime")
    totaltime = check_paths(content,"totaltime")
    year = check_paths(content,"year")
    month = check_paths(content,"month")
    day = check_paths(content,"day")
    hour = check_paths(content,"hour")
    minn = check_paths(content,"min")
    ndias = check_paths(content,"ndays")
    file_mask = check_paths(content,"file_mask")
    name_mascara = check_paths(content,"maskname")
    name_variable_lat = check_paths(content,"maskvar_lat")
    name_variable_lon = check_paths(content,"maskvar_lon")
    x_left_lower_corner = check_paths(content,"x_left_lower_corner")
    y_left_lower_corner = check_paths(content,"y_left_lower_corner")
    x_right_upper_corner = check_paths(content,"x_right_upper_corner")
    y_right_upper_corner = check_paths(content,"y_right_upper_corner")
    model = check_paths(content,"model")
    method = check_paths(content,"method")
    threshold = check_paths(content,"dqdt_threshold")
    filter_parcels_dq = check_paths(content,"filter_parcels_dqdt")
    filter_parcels_height = check_paths(content,"filter_parcels_height")
    filter_vertical_layers = check_paths(content,"filter_vertical_layers")
    output_txt = check_paths(content,"output_txt")
    output_npy = check_paths(content,"output_npy")
    output_nc = check_paths(content,"nc")
    name_target_region = check_paths(content,"name_target_region")
    value_mask = check_paths(content,"mask_value")
    key_gz = check_paths(content,"file_gz")
    save_position_part = check_paths(content,"save_position_part")
    save_position_dqdt = check_paths(content,"save_position_dqdt")
    use_vertical_layers = check_paths(content,"use_vertical_layers")
    vertical_layers = check_paths(content,"vertical_layers")
    plotting_parcels_t0 = check_paths(content,"plotting_parcels_t0")
    plotting_parcels_tracks_on_map = check_paths(content,"plotting_parcels_tracks_on_map")
    maps_limits = check_paths(content,"maps_limits")
    plotting_3Dparcels_tracks = check_paths(content,"plotting_3Dparcels_tracks")   
    noleap =  check_paths(content,"noleap")
    limit_domain =  check_paths(content,"limit_domain")
    method_wvrt =  check_paths(content,"method_wvrt")
    plotting_moisture_sink_source = check_paths(content,"plotting_moisture_sink_source") 
    limit_plot = check_paths(content,"limits_plot")

    save_position_part = str2boolean(save_position_part)
    use_vertical_layers = str2boolean(use_vertical_layers)
    save_position_dqdt = str2boolean(save_position_dqdt)
    filter_parcels_dq = str2boolean(filter_parcels_dq)
    filter_parcels_height = str2boolean(filter_parcels_height)

    plotting_parcels_t0 = str2boolean(plotting_parcels_t0)
    plotting_parcels_tracks_on_map = str2boolean(plotting_parcels_tracks_on_map)
    plotting_3Dparcels_tracks =  str2boolean(plotting_3Dparcels_tracks)
    plotting_moisture_sink_source = str2boolean(plotting_moisture_sink_source)
    noleap = str2boolean(noleap)
    limit_domain = str2boolean(limit_domain)
    method_wvrt = str2boolean(method_wvrt)

    if filter_parcels_dq:
       filter_value=1
    else:
       filter_value=0
       
    if limit_domain:
        limit_domain = 1
    else:
        limit_domain = 0

    if mode=="backward":
       paso=-1
    elif mode=="forward":
       paso=1
    elif mode=="wvrt":
       paso=-2
    elif mode=="partposit":
       paso=-3    
    else:
       print ("Parameter mode must be backward, forward, wvrt or partposit")
       raise SystemExit("Bye :)")

    if model=="FLEXPART":
       type_file=2
    elif model=="FLEXPART-WRF":
       type_file=1
    else:
       print ("Model must be FLEXPART or FLEXPART-WRF")
       raise SystemExit("Bye :)")

    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()

    list_year, list_month, list_day, list_hour, list_min=generate_fecha_simulation(ndias, year, month, day, hour, minn)
    
    if noleap:
        ind=[]
        for i in range(len(list_month)):
            if list_month[i]=="02" and list_day[i]=="29":
                ind.append(int(i))
        list_year = np.delete(list_year, ind)
        list_month = np.delete(list_month, ind)
        list_day = np.delete(list_day, ind)
        list_hour = np.delete(list_hour, ind)
        list_min = np.delete(list_min, ind)

    to_check_params(paso,type_file,numPdX , numPdY, method, resolution, file_mask)

    for year, month, day, hour, minn in zip(list_year, list_month, list_day, list_hour, list_min):

            fecha=year+"-"+month+"-"+day+" "+hour+":"+minn+":00"

            if rank==0:
                print ("                                                                                          ")
                print ('--------------------------------TROVA has started-----------------------------------------')
                print (' *****************************************************************************************')
                print (" *                    EPhysLab (Environmental Physics Laboratory), Spain                 *")
                print (" *                      Galician Supercomputing Center, Spain                            *")
                print (" *                        TRansport Of water VApor (TROVA)                               *")
                print (" *                         Version " +str(get_currentversion())+" ("+ str(get_lastupdate())+")                                    *")
                TROVA_LOGO()
                print (" *                       Edificio Campus da Auga/Edificio CESGA                          *")
                print (" *                             University of Vigo/CESGA                                  *")
                print (" *                           www.ephyslab.uvigo.es/www.cesga.es                          *")
                print (" *      contact: jose.carlos.fernandez.alvarez@uvigo.es (jcfernandez@cesga.es),          *")
                print (" *                            albenis.perez.alarcon@uvigo.es                             *")
                print (" *****************************************************************************************")
                print ('------------------------------------------------------------------------------------------')
                print ("                                                                                          ")
                print ("----------------------------------- RUN INFORMATION --------------------------------------\n")
                print ('+ Configuration file ->  '+input_file)
                
                if paso==1 or paso==-1:
                    if method==1:
                        print ("+ You are using methodology of Stohl and James (2005) (DOI:10.1175/1525-7541(2004)005<0656:ALAOTA>2.0.CO;2)")
                    if method==2:
                        print ("+ You are using methodology of Sodemann et al. (2008) (DOI:10.1002/2017JD027543)")
                    print ('                                                                                          ')
                    print ('+ Target region ->  '+name_target_region)
                    print ('+ CPUs for tracking ->   '+ str(size))
                    print("+ Tracking mode -> " + str(mode))
                    print("+ Simulation starts -> " + fecha)
                    print("+ Lagrangian Model -> " + str(model))
                    print("+ Filter precipitating parcels -> " + str(filter_parcels_dq))
                    print("+ Filter parcels by height -> " + str(filter_parcels_height))
                    print("+ Mask file -> " + file_mask)
                    print("+ Save parcels' positions at each time step -> " + str(save_position_part))
                    print("+ Save dqdt at each dt -> " + str(save_position_dqdt))
                    print("+ Plot identified parcels within the target region at time t0 -> " + str(plotting_parcels_t0))
                    print("+ Plot identified parcels trajectories on a map -> " + str(plotting_parcels_tracks_on_map))
                    print("+ Plot 3D parcels trajectories -> " + str(plotting_3Dparcels_tracks))
                    print("+ Plot the pattern of moisture sources or sinks -> " + str(plotting_moisture_sink_source))
                    print("+ Calculate the residence time of water vapor -> " + str(method_wvrt))
                if paso==-2:
                    print ("+ You are using methodology of Läderach and Sodemann (2016) to calculate the residence time of water vapor (DOI: 10.1002/2015GL067449)")
                    print ('+ Target region ->  '+name_target_region)
                    print ('+ CPUs for tracking ->   '+ str(size))
                    print("+ Tracking mode -> " + str(mode))
                    print("+ Simulation starts -> " + fecha)
                    print("+ Lagrangian Model -> " + str(model))
                    print("+ Mask file -> " + file_mask)
                if paso==-3:
                    print ("+ You are saving the variables associated with the particles within the target region")
                    print ('+ Target region ->  '+name_target_region)
                    print ('+ CPUs for tracking ->   '+ str(size))
                    print("+ Tracking mode -> " + str(mode))
                    print("+ Simulation starts -> " + fecha)
                    print("+ Lagrangian Model -> " + str(model))
                    print("+ Mask file -> " + file_mask)

                print ('                                                                                          ')
                print ('------------------------------------------------------------------------------------------')
                print ('                              PROCESSING PARTPOSIT FILES                                  ')
                print ('------------------------------------------------------------------------------------------')

            start_time = time()
            main_process(path,paso, comm, size, rank, resolution, numPdX, numPdY, dtime,totaltime,
               year, month, day, hour, minn,start_time, path_output,file_mask,name_mascara,name_variable_lon, name_variable_lat,x_lower_left,y_lower_left, 
               type_file,masa,numP,x_left_lower_corner,y_left_lower_corner,x_right_upper_corner,y_right_upper_corner, model, method,threshold, filter_value,
               output_txt,output_npy,output_nc, value_mask, key_gz, save_position_part,use_vertical_layers, vertical_layers,save_position_dqdt, filter_parcels_height,
               filter_vertical_layers, plotting_parcels_t0, plotting_parcels_tracks_on_map, plotting_3Dparcels_tracks, maps_limits, noleap, limit_domain, method_wvrt, 
               plotting_moisture_sink_source, limit_plot)
            elapsed_time = time() - start_time
            if rank==0:
                print ('                                                                                          ')
                print ('-------------------------------TROVA has ended--------------------------------------------')
                print ("Congratulations, the run has been successfully completed for " +fecha +". Run time: %.2f seconds." % np.round(elapsed_time, 2))
                print ('------------------------------------------------------------------------------------------')