#!/usr/bin/env python
import numpy as np
import sys
from mpi4py import MPI
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
import imp
import matplotlib.pylab as plt
from mpl_toolkits.mplot3d import Axes3D
import math
from .functions import k_dq as K_dq
from .functions import k_dq_layers as K_dq_layers
from .functions import read_binary_file as RBF
from .functions import determined_id as D_id
from .functions import search_row as sRow
from .functions import len_file as lf
from .functions import kdif as kdif
from .functions import k_dq_so as K_dq_So
from .functions import k_dq_por as K_dq_POR
from .functions import filter_part as Filter_Part
from .functions import filter_part2 as Filter_Part2
from .functions import filter_part_by_height as Filter_by_Height


warnings.filterwarnings("ignore", category=DeprecationWarning)
print = functools.partial(print, flush=True)


def check_paths(pfile, path):
	try:
			fpath = getattr(pfile, path)
	except:
			fpath = ""
	return fpath


def str2boolean(arg):
    if isinstance(arg, bool):
        return arg
    if arg.lower() in ("yes", "true", "t", "y", "1"):
        return True
    elif arg.lower() in ("no", "false", "f", "n", "0"):
        return False
    else:
        raise argparse.ArgumentTypeError("Boolean value expected.")

def ProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '»', printEnd = "\r"):
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end = "\r")
    if iteration == total:
        print()


def get_currentversion():
    pathpkg = os.path.dirname(__file__)
    version_file = pathpkg+"/VERSION"

    with open(version_file) as vfile:
        version = vfile.readlines()[0].strip()
    return(version)


def get_lastupdate():
    pathpkg = os.path.dirname(__file__)
    update_file = pathpkg+"/LAST_UPDATE"

    with open(update_file) as ufile:
        lastupdate = ufile.readlines()[0].strip()
    return(lastupdate)


def plotting_tracks_3d(particle_positions, fname):
	fig = plt.figure(figsize=(15,15))

	ax = fig.add_subplot(111, projection='3d')
	ProgressBar (0, particle_positions.shape[1], prefix = " Plotting 3D parcels' tracks       ------>", suffix = '', decimals = 1, length = 40, printEnd = "\r")
	for i in range(0,particle_positions.shape[1]):
		lat=particle_positions[:,i,2]
		lon=particle_positions[:,i,1]
		ids=particle_positions[:,i,0]
		z=particle_positions[:,i,4]

		if ids[0]!=-999.9:
			if all(value>=-180 for value in lon):
				ax.plot3D(lon, lat, z/1000)

		pytime.sleep(0.0001)
		ProgressBar (i+1, particle_positions.shape[1], prefix = " Plotting 3D parcels' tracks       ------>", suffix = '', decimals = 1, length = 40, printEnd = "\r")
	plt.xticks(fontsize=25)
	plt.yticks(fontsize=25)

	plt.xlabel("Longitude", fontsize=25, labelpad=25)
	plt.ylabel("Latitude", fontsize=25, labelpad=25)
	ax.set_zlabel("Height (km)", fontsize=25, labelpad=25)

	ax.tick_params(axis='z', labelsize=25)

	plt.savefig(fname,dpi=600)
	plt.close()

def ploting_parcels_tracks_map(particle_positions, maps_limits, paso, lat_masked, lon_masked, mascara, value_mask, fname):
	plt.figure(figsize=(18,12))
	mapa,crs=create_map(maps_limits)

	ProgressBar (0, particle_positions.shape[1], prefix = " Plotting parcels' tracks on a map   ---->", suffix = '', decimals = 1, length = 40, printEnd = "\r")
	for i in range(0,particle_positions.shape[1]):
		lat=particle_positions[:,i,2]
		lon=particle_positions[:,i,1]
		ids=particle_positions[:,i,0]

		if ids[0]!=-999.9:
			if all(value>=-180 for value in lon):
				mapa.plot(lon, lat, transform=crs)

				if paso==-1:
					mapa.plot(lon[0], lat[0], color="k", marker="o",markersize=5, transform=crs)
				if paso==1:
					mapa.plot(lon[-1], lat[-1], color="k", marker="o",markersize=5, transform=crs)
		pytime.sleep(0.00001)
		ProgressBar (i+1, particle_positions.shape[1], prefix = " Plotting parcels' tracks on a map   ---->", suffix = '', decimals = 1, length = 40, printEnd = "\r")
	mapa.contour(lon_masked, lat_masked, mascara, value_mask, colors="b", linewidths=4)
	plt.savefig(fname, bbox_inches="tight",dpi=600)
	plt.close()

def create_map(maps_limits):
	import cartopy.crs as ccrs
	import cartopy.feature as cfeature
	import matplotlib.ticker as mticker
	from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
	from cartopy.io.shapereader import Reader
	from cartopy.feature import ShapelyFeature

	latmin=maps_limits[0]
	lonmin=maps_limits[1]
	latmax=maps_limits[2]
	lonmax=maps_limits[3]
	center=maps_limits[4]
	dlat=int(maps_limits[5])
	dlon=int(maps_limits[6])

	crs = ccrs.PlateCarree()
	mapa=plt.subplot(111,projection=ccrs.PlateCarree(center))
	mapa.add_feature(cfeature.COASTLINE.with_scale('10m'), linewidth=0.95)
	mapa.add_feature(cfeature.BORDERS,linestyle="-", linewidth=0.75)
	mapa.add_feature(cfeature.LAKES, alpha=0.5)
	mapa.set_extent([lonmin,lonmax,latmin,latmax], crs=ccrs.PlateCarree())

	gl = mapa.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=0.5, color='black', alpha=1, linestyle='--')
	if int(lonmin)!=int(lonmax):
		lons=np.arange(int(math.ceil(lonmin)),int(math.floor(lonmax))+dlon,dlon)
	else:
		lons=np.arange(lonmin,lonmax+dlon,dlon)

	gl_lon_info=[]
	if center==180:
		for clons in lons:
			if clons<180:
				gl_lon_info=np.append(gl_lon_info,clons)
			else:
				gl_lon_info=np.append(gl_lon_info,clons-360)
	elif center==0:
		gl_lon_info=lons

	gl_loc=[True,False,False,True]
	if float(sys.version[0:3]) >= (3.7):
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
	if dlat>=1:
		gl.ylocator = mticker.MultipleLocator(dlat)
	gl.xformatter = LONGITUDE_FORMATTER
	gl.yformatter = LATITUDE_FORMATTER
	gl.xlabel_style = {'size': 25, 'color': 'k'}
	gl.ylabel_style = {'size': 25, 'color': 'k'}
	return mapa,crs


def plotting_parcels_within_target_region(particle_positions, maps_limits, paso, lat_masked, lon_masked, mascara, value_mask, fname):
	plt.figure(figsize=(18,12))
	mapa,crs=create_map(maps_limits)

	if paso==-1:
		idx=-1
	if paso==1:
		idx=0
	ProgressBar (0, particle_positions.shape[1], prefix = " Plotting parcels within the target region", suffix = '', decimals = 1, length = 40, printEnd = "\r")
	for i in range(0,particle_positions.shape[1]):
		lat=particle_positions[idx,i,2]
		lon=particle_positions[idx,i,1]
		ids=particle_positions[idx,i,0]

		if ids!=-999.9:
			if lon>=-180:
				mapa.plot(lon, lat, color="r", marker="o", markersize=3, transform=crs)
		pytime.sleep(0.001)
		ProgressBar (i+1, particle_positions.shape[1], prefix = " Plotting parcels within the target region", suffix = '', decimals = 1, length = 40, printEnd = "\r")
	mapa.contour(lon_masked, lat_masked, mascara, value_mask, colors="b", linewidths=4)
	plt.savefig(fname,bbox_inches="tight",dpi=600)
	plt.close()

def generate_fecha_simulation(ndias, cyear, cmonth, cday, chours, cminutes ):

	nhour=int(ndias*24)
	year=[]
	mes=[]
	dia=[]
	hora=[]
	mins=[]
	array =np.arange(0,nhour,24)

	if not isinstance(chours, list):
		chours=[chours]
	if not isinstance(cminutes, list):
		cminutes=[cminutes]
	if not isinstance(cyear, list):
		cyear=[cyear]
	if not isinstance(cmonth, list):
		cmonth=[cmonth]
	if not isinstance(cday, list):
		cday=[cday]


	for i in array:

		for yy in cyear:
			yy=str(int(yy)).zfill(4)
			for mm in cmonth:
				mm=str(int(mm)).zfill(2)
				for dd in cday:
					dd=str(int(dd)).zfill(2)
					for hh in chours:
						for mmin in cminutes:
							fecha=yy+"-"+mm+"-"+dd+" "+str(int(hh)).zfill(2)+":"+str(int(mmin)).zfill(2)+":00"
							a=str(time_calc(fecha,float(i)))
							var1=a.split(" ")
							var11=var1[0].split("-")
							var12=var1[1].split(":")
							year_=str(var11[0])
							mes_=str(var11[1])
							dia_=str(var11[2])
							hora_=str(var12[0])
							minn_=str(var12[1])
							year.append(year_)
							mes.append(mes_)
							dia.append(dia_)
							hora.append(hora_)
							mins.append(minn_)

	return year, mes, dia, hora,mins


def function(latitude, longitude, var, var_layers, use_vlayers, vlayers, method, varpor, filename,path, name_var, unit_var,date_save):

    ncout = Dataset(path+filename+".nc", 'w', format='NETCDF4')
    ncout.createDimension('lat', len(latitude))
    ncout.createDimension('lon', len(longitude))
    ncout.createDimension('time', len(var[:,0,0]))
    time_var=np.arange(len(var[:,0,0]))

    if use_vlayers:
       ncout.createDimension('layers', len(vlayers)-1)
       layer_vec=[]
       for ilayer in range(0,len(vlayers)-1):

           clayer=str(int(vlayers[ilayer]))+"_"+str(int(vlayers[ilayer+1]))
           layer_vec=np.append(layer_vec, clayer)

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
    time.calendar="gregorian"
    time.description = "days since 1900-01-01"
    time.units = "days since 1900-01-01"

    vout = ncout.createVariable(name_var, np.dtype('float').char, ('time','lat','lon'),zlib=True)
    vout.long_name = name_var
    vout.units = unit_var
    vout.standard_name = name_var
    vout.coordinates = "time, lat,lon"

    if name_var=="E_P":
        vvout = ncout.createVariable("E_P_integrated", np.dtype('float').char, ('lat','lon'),zlib=True)
        vvout.long_name = "E_P integrated for ndays considered"
        vvout.units = unit_var
        vvout.standard_name = "E_P_integrate"
        vvout.coordinates = "lat,lon"
    vout.original_name = name_var

    if method==2:
         voutpor = ncout.createVariable("POR", np.dtype('float').char, ('time','lat','lon'),zlib=True)
         voutpor.long_name = "Sources Contribution for each parcel"
         voutpor.units = "%"
         voutpor.standard_name = "Sources Contribution for each parcel"
         voutpor.coordinates = "time, lat,lon"
         voutpor.original_name = "Sources Contribution for each parcel"
         voutpor.coordinates = "time, lat,lon"

    if use_vlayers:
       voutlayer = ncout.createVariable(name_var+"_layers", np.dtype('float').char, ('time','layers','lat','lon'),zlib=True)
       voutlayer.long_name = name_var+"_layers"
       voutlayer.units = unit_var
       voutlayer.standard_name = name_var
       voutlayer.coordinates = "time, layers, lat,lon"
       voutlayer.original_name = name_var+"_layers"


       vvoutlayer = ncout.createVariable("E_P_integrated_layers", np.dtype('float').char, ('layers','lat','lon'),zlib=True)
       vvoutlayer.long_name = "E_P integrated for ndays in vertical layers"
       vvoutlayer.units = unit_var
       vvoutlayer.standard_name = "E_P_integrated_by_layers"
       vvoutlayer.coordinates = "layers, lat,lon"
       vvoutlayer.original_name = "E_P integrated for ndays in vertical layers"

       v_layers = ncout.createVariable('vertical_layers', 'str', 'layers')
       v_layers.long_name = 'vertical layers'
       v_layers.units = 'vertical layers'
       v_layers.axis = 'vertical layers'

    lon[:] = longitude
    lat[:]= latitude
    time[:]= date_save[:]
    vout[:,:,:] = var
    if name_var=="E_P":
        vvout[:,:] = np.sum(var, axis=0)

    if method==2:
         voutpor[:]=varpor

    if use_vlayers:
        voutlayer[:,:,:,:]=var_layers
        vvoutlayer[:,:,:]=np.sum(var_layers, axis=0)
        v_layers[:]=layer_vec
    ncout.close()

def write_nc(dates, tensor, vartype,filename="output"):
	ncout = Dataset(filename+".nc", 'w', format='NETCDF4')
	ncout.history='Parcels positions'

	if vartype=="partpos":
		ncout.history= "partpos[:,:,0] - parcel ID, partpos[:,:,1] - longitude, partpos[:,:,2] - latitude, partpos[:,:,3] - specific humidity, partpos[:,:,4] - vertical position (m), partpos[:,:,5] - topography high (m), partpos[:,:,6] - density (kg/m3), partpos[:,:,7] - PBL high (m), partpos[:,:,8] - Tropopause high (m), partpos[:,:,9] - temperature (K), partpos[:,:,10] - parcel mass (kg)"
	if vartype=="dqdt":
		ncout.history= "partpos[:,:,0] - longitude, partpos[:,:,1] - latitude, partpos[:,:,2] - dq/dt, partpos[:,:,3] - vertical position (m), partpos[:,:,4] - parcel ID, partpos[:,:,5] - specific humidity at starting trcking point (time t0)"

	ndates=len(dates)
	npart=tensor.shape[1]
	vprop=tensor.shape[2]

	ncout.createDimension('time', ndates)
	ncout.createDimension('npart', npart)
	ncout.createDimension('properties', vprop)

	times = ncout.createVariable('times', dtype('int').char, ('time'))
	times.standard_name = 'times'
	times.long_name = 'times (YYYYMMDDHHMM)'
	times.units = ''
	times.axis = ''

	parts = ncout.createVariable('parcels', dtype('float32').char, ('npart'))
	parts.standard_name = 'parcels'
	parts.long_name = 'Parcels IDs'
	parts.units = ''
	parts.axis = ''

	vout = ncout.createVariable('partpos', dtype('float32').char, ('time','npart','properties'))
	vout.long_name = 'Parcels position'
	vout.units = ''
	vout.standard_name = "Parcels position";
	vout.coordinates = "times,npart, properties" ;
	vout.original_name = "Parcels position"

	times[:] = dates
	parts[:]= tensor[0,:,0]
	vout[:] = tensor[:,:,:]
	ncout.close()


def create_directory(path):

    try:
        if not os.path.exists(path):
            os.mkdir(path)
    except OSError:
        pass

def read_binaryFile_fortran(filename, type_file,x_left_lower_corner,y_left_lower_corner,x_right_upper_corner,y_right_upper_corner):

    if type_file==1:
        with open(filename,'rb') as inputfile:
            a=b''.join([line for line in inputfile])
        npart=struct.unpack('iiii', a[0:16])
        npart=npart[2]
        data= RBF(filename,npart,x_left_lower_corner,y_left_lower_corner, x_right_upper_corner,y_right_upper_corner)

    if type_file==2:
        len_a=lf(filename)
        npart=((len_a-12)/60)-1
        data= RBF(filename,npart, x_left_lower_corner,y_left_lower_corner, x_right_upper_corner,y_right_upper_corner)
    ind=np.where(data[:, 0]==-999.)
    data=data[:int(ind[0][0]), :]

    return data

def load_mask_grid_NR(filename, name_mascara,name_variable_lon, name_variable_lat):

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
    ax1.set_extent([-170.,60.,-30,85], crs=ccrs.PlateCarree())

    mascara=mascara.astype("bool")
    plt.scatter(lon, lat, s=5)
    plt.scatter(-60,20,s=20,color="r")
    gl = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=1, color='gray', alpha=0.2, linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_left = True
    gl.ylabels_right = False
    gl.xlines = True

    paso_h=10
    lons=np.arange(np.ceil(-110),np.ceil(40),paso_h)
    gl.xlocator = mticker.FixedLocator(lons)
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 15, 'color': 'gray'}
    gl.xlabel_style = {'color': 'black'}

    plt.savefig("mask.png")
    plt.close()

def funtion_interpol_mascara_2(lat_mascara, lon_mascara, mascara, data):

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

    value_mascara=funtion_interpol_mascara (lat_mascara, lon_mascara, value_mascara, data)
    id_vector=np.array(D_id(value_mascara, value_mask, len(value_mascara)),dtype=int)
    submatrix=[]
    ind=[]
    for ii in id_vector:
        if ii !=-999:
            submatrix=np.append(submatrix, data[ii,:])
            ind.append(ii)
    submatrix=np.reshape(submatrix,(len(ind), 11))
    return submatrix

def search_row_fortran(lista, matrix):

    matrix_=np.array(sRow(matrix, lista, len(lista), len(matrix[:,0])), np.float64)

    return matrix_

def calc_A(resolution, lat, lon):

    rt = 6371000.
    gr = np.pi/180.
    a,b=lat.shape
    area=np.empty((a-1,b-1))
    area[:,:]=0
    for j in range(len(lat[0,:])-1):
        for i in range(len(lat[:,0])-1):
            area[i,j]=np.abs((gr*rt**2)*( np.sin(gr*lat[i,j]) - np.sin(gr*lat[i+1,j])))*np.abs(resolution)
    return area

def grid_point (resolution, numPdX, numPdY, x_lower_left,y_lower_left):

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

    lat_new=[]
    lon_new=[]
    for i in range(len(lat[:,0])-1):
        lat_new= np.append(lat_new, (lat[i+1,0]+ lat[i,0])/2.)
    for j in range(len(lon[0,:])-1):
        lon_new= np.append(lon_new, (lon[0,j+1]+ lon[0,j])/2.)
    lon_plot, lat_plot=np.meshgrid(lon_new, lat_new)
    return lat_plot,lon_plot

def time_calc(init_time,h_diff):

    formatted_time = datetime.strptime(init_time, "%Y-%m-%d %H:%M:%S")
    calculated_time=formatted_time+timedelta(hours=h_diff)
    return calculated_time

def time_calcminutes(init_time,h_diff):

    formatted_time = datetime.strptime(init_time, "%Y-%m-%d %H:%M:%S")
    calculated_time=formatted_time+timedelta(minutes=h_diff)
    return calculated_time

def generate_file(paso, dtime, totaltime, fecha, path, key_gz):

    nhour=int(totaltime)+dtime
    list_fecha=[]
    listdates=[]
    if paso == -1:
        array =np.arange(nhour,0, -dtime)
        for i in array:
            a=str(time_calcminutes(fecha,float(i)*(-1)))
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
        fecha_=fecha.split(" ")
        var11=fecha_[0].split("-")
        var12=fecha_[1].split(":")
        fecha_dia=str(var11[0]+var11[1]+var11[2]+var12[0]+var12[1]+var12[2])
        if key_gz==1:
            name=path+"partposit_"+fecha_dia+".gz"
        else:
            name=path+"partposit_"+fecha_dia
        list_fecha=np.append(list_fecha, name)
        listdates=np.append(listdates, int(fecha_dia))
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

    return list_fecha, listdates

def read_proccesor(lista_partposi,submatrix, rank, x_left_lower_corner,y_left_lower_corner,
               x_right_upper_corner,y_right_upper_corner, model, key_gz, type_file):

    a1=np.arange(len(lista_partposi))
    dx,dy =submatrix.shape
    tensor_local=np.ones((len(lista_partposi),dx,dy))*(-999.9)
    for i in a1:
        print ("Reading | " + model+" -> ",  lista_partposi[i])
        if key_gz==1:
            desc_gz(lista_partposi[i])
            part_post_i=read_binaryFile_fortran(lista_partposi[i][:-3], type_file,x_left_lower_corner,y_left_lower_corner,
               x_right_upper_corner,y_right_upper_corner)
            cmd_rm= "rm -rf "+lista_partposi[i][:-3]
            os.system(cmd_rm)
        else:
            part_post_i=read_binaryFile_fortran(lista_partposi[i], type_file,x_left_lower_corner,y_left_lower_corner,
               x_right_upper_corner,y_right_upper_corner)
        matrix_i=search_row_fortran(submatrix[:,0],part_post_i)
        tensor_local[i,:,:]=matrix_i
    return tensor_local

def _backward_dq(lista_partposi ,file_mask, name_mascara,name_variable_lon, name_variable_lat,lat_f, lon_f,rank,size, comm, type_file,
                 x_left_lower_corner,y_left_lower_corner, x_right_upper_corner,y_right_upper_corner, model, method,threshold,filter_value, value_mask, key_gz, path_output,use_vertical_layers, vertical_layers, filter_parcels_height, filter_vertical_layers):

    name_file=lista_partposi[-1]
    if rank==0:
        print ("Reading | " + model+" -> ",  name_file)
        name_txt_part=name_file.split("/")
        f=open(path_output+name_txt_part[-1].split("_")[-1].split(".")[0]+"/"+name_txt_part[-1].split("_")[-1].split(".")[0]+".txt", "a")
    if key_gz==1:
        desc_gz(name_file)
        part_post=read_binaryFile_fortran(name_file[:-3], type_file, x_left_lower_corner,y_left_lower_corner,
               x_right_upper_corner,y_right_upper_corner)
        cmd_rm= "rm -rf "+name_file[:-3]
        os.system(cmd_rm)
    else:
        part_post=read_binaryFile_fortran(name_file, type_file, x_left_lower_corner,y_left_lower_corner,
               x_right_upper_corner,y_right_upper_corner)


    lat_masked, lon_masked, mascara=load_mask_grid_NR(file_mask, name_mascara,name_variable_lon, name_variable_lat)

    submatrix=determine_id_binary_grid_NR_fortran(part_post, lat_masked.flatten(), lon_masked.flatten(), mascara.flatten(), value_mask)
    submatrix=submatrix[np.argsort(submatrix[:, 0])]


    if filter_parcels_height:
       submatrix, counter_part_height = Filter_by_Height(submatrix,submatrix,-1,filter_vertical_layers[0],filter_vertical_layers[1], len(submatrix[:, 0]))

    if rank==0:
        print ("Reading | " + model+" -> ",  lista_partposi[-2])

    if key_gz==1:
        desc_gz(lista_partposi[-2])
        part_post_i=read_binaryFile_fortran(lista_partposi[-2][:-3], type_file, x_left_lower_corner,y_left_lower_corner,x_right_upper_corner,y_right_upper_corner)
        cmd_rm= "rm -rf "+lista_partposi[-2][:-3]
        os.system(cmd_rm)
    else:
        part_post_i=read_binaryFile_fortran(lista_partposi[-2], type_file, x_left_lower_corner,y_left_lower_corner,x_right_upper_corner,y_right_upper_corner)
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
    local_results= read_proccesor(local_list, submatrix, rank,x_left_lower_corner,y_left_lower_corner,
               x_right_upper_corner,y_right_upper_corner, model, key_gz, type_file)

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

    matrix_results=np.ones((len(tensor_org[:,0,0])-1, len(submatrix[:,0]),4))*(-999.9)
    a3=np.arange(len(tensor_org[:,0,0])-1)

    for i in a2[::-1]:

        matrix=kdif(tensor_t[i,:,:5], tensor_t[i+1,:,:5],-1.,len(submatrix[:,0]), 5)
        matrix_result[i,:,2]=matrix[:,2]
        matrix_result[i,:,1]=matrix[:,1]
        matrix_result[i,:,0]=matrix[:,0]
        matrix_result[i,:,3]=matrix[:,3]

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

    return matrix_result, tensor_t[-1,:,0], matrix_result_por, tensor_t[-1,:,3], tensor_org, lat_masked, lon_masked, mascara

def _forward_dq(lista_partposi ,file_mask, name_mascara,name_variable_lon, name_variable_lat,lat_f, lon_f,rank,size,comm, type_file,
                x_left_lower_corner,y_left_lower_corner, x_right_upper_corner,y_right_upper_corner, model, value_mask, key_gz,path_output,use_vertical_layers, vertical_layers, filter_parcels_height, filter_vertical_layers):
    name_file=lista_partposi[0]
    if rank==0:
        print ("Reading | " + model+" -> ",  name_file)
        name_txt_part=name_file.split("/")
        f=open(path_output+name_txt_part[-1].split("_")[-1].split(".")[0]+"/"+name_txt_part[-1].split("_")[-1].split(".")[0]+".txt", "a")
    if key_gz==1:
        desc_gz(name_file)
        part_post=read_binaryFile_fortran(name_file[:-3], type_file,x_left_lower_corner,y_left_lower_corner, x_right_upper_corner,y_right_upper_corner)
        cmd_rm= "rm -rf "+name_file[:-3]
        os.system(cmd_rm)
    else:
        part_post=read_binaryFile_fortran(name_file, type_file,x_left_lower_corner,y_left_lower_corner, x_right_upper_corner,y_right_upper_corner)

    lat_masked, lon_masked , mascara=load_mask_grid_NR(file_mask, name_mascara,name_variable_lon, name_variable_lat)
    submatrix=determine_id_binary_grid_NR_fortran(part_post, lat_masked.flatten(), lon_masked.flatten(), mascara.flatten(), value_mask)
    submatrix=submatrix[np.argsort(submatrix[:, 0])]


    if filter_parcels_height:
       submatrix, counter_part_height = Filter_by_Height(submatrix,submatrix,-1,filter_vertical_layers[0],filter_vertical_layers[1], len(submatrix[:, 0]))
    if rank==0:
        print ("Reading | " + model+" -> ",  lista_partposi[1])
    if key_gz==1:
        desc_gz(lista_partposi[1])
        part_post_i=read_binaryFile_fortran(lista_partposi[1][:-3], type_file, x_left_lower_corner,y_left_lower_corner, x_right_upper_corner,y_right_upper_corner)
        cmd_rm= "rm -rf "+lista_partposi[1][:-3]
        os.system(cmd_rm)
    else:
        part_post_i=read_binaryFile_fortran(lista_partposi[1], type_file, x_left_lower_corner,y_left_lower_corner, x_right_upper_corner,y_right_upper_corner)
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
    local_results= read_proccesor(local_list, submatrix, rank,x_left_lower_corner,y_left_lower_corner,x_right_upper_corner,y_right_upper_corner, model, key_gz, type_file)

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
        matrix=kdif(tensor_t[i,:,:5], tensor_t[i+1,:,:5],1.,len(submatrix[:,0]), 5)
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

def time_calc_day(init_time,day_diff):

    formatted_time = datetime.strptime(init_time, "%Y-%m-%d %H:%M:%S")
    calculated_time=formatted_time+timedelta(days=int(day_diff))
    return calculated_time

def convert_date_to_ordinal(year, month, day, hour, minute, second):
    date=datetime(year, month, day, hour, minute, second)
    date_=netCDF4.date2num(date, "days since 1900-01-01", "gregorian")
    return date_

class InputNotInRangeError(Exception):

    def __init__(self, message):
        self.message = message
        super().__init__(self.message)

def to_check_params(paso,type_file,numPdx , numPdy, method, resolution, cant_plazo, file_mask):

    if paso==1 or paso==-1:
        pass
    else:
        raise InputNotInRangeError("Only forward (mode = 1) and backward mode (mode = - 1) are allowed")

    if type_file==1 or type_file==2:
        pass
    else:
        raise InputNotInRangeError("Only FLEXPART-WRF (type_file = 1) and FLEXPART model (type_file = 2) are allowed")


    if method==1 or method==2:
        pass
    else:
        raise  InputNotInRangeError("Only Sthol and James ( method = 1) and Sodemman methodology (method = 2) are allowed")

    if numPdx <=0  or  numPdy<=0:
        raise InputNotInRangeError(" numPdX  and  numPdY  must be greater than zero ")
    if resolution<=0:
        raise InputNotInRangeError(" resolution must be greater than zero ")
    my_file=Path(file_mask)
    if not my_file.is_file():
        raise InputNotInRangeError("File mask no found")

def function_proof(lat, lon):

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

    with gzip.open(name_file, 'rb') as f_in:
        with open(name_file[:-3], 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)


def TROVA_LOGO():
    print(" *                        _____ __    ____                                               *")
    print(" *                          |  |  |  /    \ \        //\                                 *")
    print(" *                          |  |__| /      \ \      //__\                                *")
    print(" *                          |  |  \ \      /  \    //    \                               *")
    print(" *                          |  |   \ \____/    \__//      \                              *")
    print(" *                                                                                       *")

def main_process(path, paso, comm, size, rank, resolution, numPdX, numPdY, dtime, totaltime, year, month,
         day, hour, minn, time, path_output,file_mask, name_mascara,name_variable_lon, name_variable_lat,x_lower_left,y_lower_left, type_file,masa,numP, x_left_lower_corner,y_left_lower_corner, x_right_upper_corner,y_right_upper_corner, model, method, threshold, filter_value, output_txt, output_npy,output_nc, value_mask, key_gz, save_position_part, use_vertical_layers, vertical_layers,save_position_dqdt, filter_parcels_height,filter_vertical_layers, plotting_parcels_t0, plotting_parcels_tracks_on_map, plotting_3Dparcels_tracks, maps_limits):

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
    lista_partposi, listdates=generate_file(paso, dtime, totaltime, fecha, path, key_gz)

    if rank==0:
       if (len(lista_partposi)-2)%size!=0:
           print("TROVA ERROR: The number of processors must exactly divide the number of partposit files to process (totaltime/dtime).\n Based on your configuration file, the recommended number of processors is " + str(int(totaltime/dtime))) 
           raise SystemExit("Bye :)")
    elif (len(lista_partposi)-2)%size!=0:
       raise SystemExit()

    for i in lista_partposi:
        my_file=Path(i)
        if not my_file.is_file():
            if rank==0:
                print (i +"-> not exits" )
            raise SystemExit("Bye :)")

    if paso==1:
        matrix_result, id_part, q_ini, partpos, lat_masked, lon_masked, mascara =_forward_dq(lista_partposi, file_mask, name_mascara,name_variable_lon, name_variable_lat,lat, lon,rank,size,comm, type_file,x_left_lower_corner,y_left_lower_corner, x_right_upper_corner,y_right_upper_corner, model, value_mask, key_gz, path_output,use_vertical_layers, vertical_layers,filter_parcels_height,filter_vertical_layers)
        matrix_save=np.empty((matrix_result.shape[0],matrix_result.shape[1],6))
        matrix_save[:,:,:-2]=matrix_result
        matrix_save[:,:,-2]=id_part
        matrix_save[:,:,-1]=q_ini


        if rank==0:
            if save_position_dqdt:
               write_nc(listdates[2:], matrix_save, "dqdt", filename=path_output+folder+"/"+folder+"_dqdt_forw")
            if  save_position_part:
                 write_nc(listdates, partpos, "partpos", filename=path_output+folder+"/"+folder+"_parposit_forw")

    if paso ==-1:
        matrix_result, id_part, matrix_result_por,q_ini,partpos, lat_masked, lon_masked, mascara =_backward_dq(lista_partposi, file_mask, name_mascara,name_variable_lon, name_variable_lat,lat, lon,rank,size,comm, type_file, x_left_lower_corner, y_left_lower_corner, x_right_upper_corner, y_right_upper_corner, model, method,threshold,filter_value, value_mask, key_gz, path_output,use_vertical_layers, vertical_layers,filter_parcels_height,filter_vertical_layers)

        matrix_save=np.empty((matrix_result.shape[0],matrix_result.shape[1],6))
        matrix_save[:,:,:-2]=matrix_result
        matrix_save[:,:,-2]=id_part
        matrix_save[:,:,-1]=q_ini

        if rank==0:

            if save_position_dqdt:
               write_nc(listdates[1:-1], matrix_save, "dqdt", filename=path_output+folder+"/"+folder+"_dqdt_back")
            if  save_position_part:

               write_nc(listdates, partpos, "partpos",filename=path_output+folder+"/"+folder+"_parposit_back")

    if rank==0:
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


        ProgressBar (0, len(t)-1, prefix = " -> Processing data", suffix = '', decimals = 1, length = 65, printEnd = "\r")
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

            pytime.sleep(0.0001)
            ProgressBar (ii+1, len(t)-1, prefix = " -> Processing data", suffix = '', decimals = 1, length = 65, printEnd = "\r")
            
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
        if paso == 1:
            if output_nc!=0:
                array_day_por=np.empty_like(array_day)
                array_day_por[:]=-999.9
                function(lat_plot[:,0], lon_plot[0,:], array_day, array_day_layers, use_vertical_layers, vertical_layers,  method, array_day_por, "forw_"+folder, path_output+folder+"/", "E_P", "mm/day", date_save[:-1])
            if output_npy!=0:
                np.save(path_output+folder+"/"+"forw_"+folder+".npy", array_day)
                if use_vertical_layers:
                   np.save(path_output+folder+"/"+"forw_layers_"+folder+".npy", array_day_layers)

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
    
    content = imp.load_source("", input_file)
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

    save_position_part = str2boolean(save_position_part)
    use_vertical_layers = str2boolean(use_vertical_layers)
    save_position_dqdt = str2boolean(save_position_dqdt)
    filter_parcels_dq = str2boolean(filter_parcels_dq)
    filter_parcels_height = str2boolean(filter_parcels_height)

    plotting_parcels_t0 = str2boolean(plotting_parcels_t0)
    plotting_parcels_tracks_on_map = str2boolean(plotting_parcels_tracks_on_map)
    plotting_3Dparcels_tracks =  str2boolean(plotting_3Dparcels_tracks)

    if filter_parcels_dq:
       filter_value=1
    else:
       filter_value=0

    if mode=="backward":
       paso=-1
    elif mode=="forward":
       paso=1
    else:
       print ("Parameter mode must be backward or forward")
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


    list_year, list_month, list_day, list_hour, list_min=generate_fecha_simulation(ndias, year, month, day, hour, minn )

    to_check_params(paso,type_file,numPdX , numPdY, method, resolution, dtime,file_mask)

    for year, month, day, hour, minn in zip(list_year, list_month, list_day, list_hour, list_min):

            fecha=year+"-"+month+"-"+day+" "+hour+":"+minn+":00"

            if rank==0:
                print ("                                                                                          ")
                print ('--------------------------------TROVA has started-----------------------------------------')
                print (' *****************************************************************************************')
                print (" *                    EPhysLab (Environmental Physics Laboratory)                        *")
                print (" *                        TRansport Of water VApor (TROVA)                               *")
                print (" *                             Version " +str(get_currentversion())+" ("+ str(get_lastupdate())+")                                  *")
                TROVA_LOGO()
                print (" *                            Edificio Campus da Auga                                    *")
                print (" *                               University of Vigo                                      *")
                print (" *                                ephyslab.uvigo.es                                      *")
                print (" *  contact: jose.carlos.fernandez.alvarez@uvigo.es, albenis.perez.alarcon@uvigo.es      *")
                print (" *****************************************************************************************")
                print ('------------------------------------------------------------------------------------------')
                print ("                                                                                          ")
                print ("----------------------------------- RUN INFORMATION ------------------------------------\n")
                print ('+ Configuration file ->  '+input_file)
                if method==1:
                    print ("+ You are using methodology of Stohl and James (2005) (DOI:10.1175/1525-7541(2004)005<0656:ALAOTA>2.0.CO;2)")
                if method==2:
                    print ("+ You are using methodology of Sodemann et al. (2008) (DOI:10.1002/2017JD027543)")
                print ('+ Target region ->  '+name_target_region)
                print ('+ CPUs for tracking ->   '+ str(size))
                print("+ Tracking mode -> " + str(mode))
                print("+ Simulation starts -> " + fecha)
                print("+ Lagrangian Model -> " + str(model))
                print("+ Filter precipitating parcels -> " + str(filter_parcels_dq))
                print("+ Filter parcels by height -> " + str(filter_parcels_height))
                print("+ Mask file -> " + file_mask)
                print("")
                print("AUXILIAR TOOLS")
                print("--------------")
                print("+ Save parcels' positions at each time step -> " + str(save_position_part))
                print("+ Save dqdt at each dt -> " + str(save_position_dqdt))
                print("+ Plot identified parcels within the target region at time t0 -> " + str(plotting_parcels_t0))
                print("+ Plot identified parcels trajectories on a map -> " + str(plotting_parcels_tracks_on_map))
                print("+ Plot 3D parcels trajectories -> " + str(plotting_3Dparcels_tracks))
                print ('                                                                                          ')
                print ('------------------------------------------------------------------------------------------')
                print ('                              PROCESSING PARTPOSIT FILES                                  ')
                print ('------------------------------------------------------------------------------------------')

            start_time = time()
            main_process(path,paso, comm, size, rank, resolution, numPdX, numPdY, dtime,totaltime,
               year, month, day, hour, minn,start_time, path_output,file_mask,name_mascara,name_variable_lon, name_variable_lat,x_lower_left,y_lower_left, type_file,masa,numP,x_left_lower_corner,y_left_lower_corner,x_right_upper_corner,y_right_upper_corner, model, method,threshold, filter_value,output_txt,output_npy,output_nc, value_mask, key_gz, save_position_part,use_vertical_layers, vertical_layers,save_position_dqdt, filter_parcels_height,filter_vertical_layers, plotting_parcels_t0, plotting_parcels_tracks_on_map, plotting_3Dparcels_tracks, maps_limits)
            elapsed_time = time() - start_time
            if rank==0:
                print ('                                                                                          ')
                print ('-------------------------------TROVA has ended--------------------------------------------')
                print ("Congratulations, the run has been successfully completed for " +fecha +" \n---> Run time: %.2f seconds." % np.round(elapsed_time, 2))
                print ('------------------------------------------------------------------------------------------')
