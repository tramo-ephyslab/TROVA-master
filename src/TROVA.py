#!/usr/bin/env python
import numpy as np
import sys
from mpi4py import MPI
from time import time
import struct
import datetime
from datetime import datetime, timedelta
from netCDF4 import Dataset
import os
import scipy.interpolate as interp
import functools
from pathlib import Path
import warnings
import gzip
import shutil
import netCDF4
from netCDF4 import date2num
from functions import k_dq as K_dq
from functions import read_binary_file as RBF
from functions import determined_id as D_id
from functions import search_row as sRow
from functions import len_file as lf
from functions import kdif as kdif
from functions import k_dq_so as K_dq_So
from functions import k_dq_por as K_dq_POR
from functions import filter_part as Filter_Part
warnings.filterwarnings("ignore", category=DeprecationWarning)
print = functools.partial(print, flush=True)

def generate_fecha_simulation(ndias, fecha):

    nhour=int(ndias*24)
    year=[]
    mes=[]
    dia=[]
    hora=[]
    array =np.arange(0,nhour,24)
    for i in array:
        a=str(time_calc(fecha,float(i)))
        var1=a.split(" ")
        var11=var1[0].split("-")
        var12=var1[1].split(":")
        year_=str(var11[0])
        mes_=str(var11[1])
        dia_=str(var11[2])
        hora_=str(var12[0])
        year.append(year_)
        mes.append(mes_)
        dia.append(dia_)
        hora.append(hora_)
    return year, mes, dia, hora

def load_parameters(input_file_name):

    label_attr_map = {
        "path_data" : ["path_data", str], 
        "path_output" : ["path_output", str],
        "mode" : ["mode", float],
        "mass" : ["mass", float],
        "numP" : ["numP", float],
        "type_file" : ["type_file", int],
        "resolution" : ["resolution", float],
        "numPdX" : ["numPdX", int],
        "numPdY" : ["numPdY", int],
        "x_lower_left" : ["x_lower_left", float],
        "y_lower_left" : ["y_lower_left", float],
        "step_numbers" : ["step_numbers", int],
        "year" : ["year", str],
        "month" : ["month", str],
        "day" : ["day", str],
        "hour" : ["hour", str],
        "ndays" : ["ndays", int],
        "file_mask" : ["file_mask", str],
        "name_mascara" : ["name_mascara", str],
        "name_variable_lat" : ["name_variable_lat", str],
        "name_variable_lon" : ["name_variable_lon", str],
        "x_left_lower_corner" : ["x_left_lower_corner", float],
        "y_left_lower_corner" : ["y_left_lower_corner", float],
        "x_rigth_upper_corner" : ["x_rigth_upper_corner", float],
        "y_rigth_upper_corner" : ["y_rigth_upper_corner", float],
        "model" : ["model", str],
        "type_lon" : ["type_lon", int],
        "method" : ["method", int],
        "threshold" : ["threshold", float],
        "filter_value" : ["filter_value", int],
        "output_txt" : ["output_txt", int],
        "output_npy" : ["output_npy", int],
        "output_nc" : ["output_nc", int],
        "name_target_region" : ["name_target_region", str], 
        "value_mask" : ["value_mask", int],
        "file_gz" : ["file_gz", int],
        "save_position_part": ["save_position_part", str]
    }

    def Params (input_file_name):

        parametros=[]
        with open(input_file_name, 'r') as input_file:
            for line in input_file:
                row = line.split("=")
                label = row[0]
                data = row[1:]
                label=label.replace(" ", "")
                attr = label_attr_map[label][0]
                datatypes = label_attr_map[label][1:]
                data_=[]
                for i in data:
                    data_.append(i.replace(" ", "").rstrip("\n"))
                values = [(datatypes[i](data_[i])) for i in range(len(data_))]
                parametros.append(values[0])
        return parametros
    params = Params(input_file_name)
    return params

def convert_date_to_ordinal(year, month, day, hour, minute, second):
    date=datetime(year, month, day, hour, minute, second)
    date_=netCDF4.date2num(date, "days since 1900-01-01", "gregorian")
    return date_

def function(latitude, longitude, var, filename,path, name_var, unit_var, date_save):

    ncout = Dataset(path+filename+".nc", 'w', format='NETCDF4')
    ncout.createDimension('lat', len(latitude))
    ncout.createDimension('lon', len(longitude))
    ncout.createDimension('time', len(var[:,0,0]))
    time_var=np.arange(len(var[:,0,0]))

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
    time.axis = 't'
    time.calendar="gregorian"
    time.description = "days since 1900-01-01" ;
    time.units = "days since 1900-01-01" ;


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
    lon[:] = longitude
    lat[:]= latitude
    time[:]= date_save[:]
    vout[:,:,:] = var
    if name_var=="E_P":
        vvout[:,:] = np.sum(var, axis=0)
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

def load_mask_grid_NR(filename, name_mascara,name_variable_lon, name_variable_lat, type_lon):

    wrfile = Dataset(filename)
    lat  = wrfile.variables[name_variable_lat][:]
    lon  = wrfile.variables[name_variable_lon][:]
    mask  = wrfile.variables[name_mascara][:]

    if type_lon==1:
        for i in range(len(lon)):
            if lon[i]>180:
                lon[i]=lon[i]-360
    if type_lon==2:
        lon=lon
    lon, lat=np.meshgrid(lon,lat)
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
    submatrix=np.reshape(submatrix,(len(ind), 5))
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

def time_calc_day(init_time,day_diff):

    formatted_time = datetime.strptime(init_time, "%Y-%m-%d %H:%M:%S")
    calculated_time=formatted_time+timedelta(days=day_diff)
    return calculated_time

def generate_file(paso, cant_plazo, fecha, path, key_gz):

    n=cant_plazo/4
    nhour=int(n*24)+6
    list_fecha=[]

    if paso == -1:
        array =np.arange(nhour,0, -6)
        for i in array:
            a=str(time_calc(fecha,float(i)*(-1)))
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

        fecha_=fecha.split(" ")
        var11=fecha_[0].split("-")
        var12=fecha_[1].split(":")
        fecha_dia=str(var11[0]+var11[1]+var11[2]+var12[0]+var12[1]+var12[2])
        if key_gz==1:
            name=path+"partposit_"+fecha_dia+".gz"
        else:
            name=path+"partposit_"+fecha_dia
        list_fecha=np.append(list_fecha, name)

    if paso == 1:
        array =np.arange(0,nhour+6, 6)
        for i in array:
            a=str(time_calc(fecha,float(i)))
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
    return list_fecha

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

def _backward_dq(lista_partposi ,file_mask, name_mascara,name_variable_lon, name_variable_lat,lat_f, lon_f,rank,size, type_file,
                 x_left_lower_corner,y_left_lower_corner, x_right_upper_corner,y_right_upper_corner, model, type_lon, method,threshold,filter_value, value_mask, key_gz, path_output):

    name_file=lista_partposi[-1]
    if rank==0:
        print ("Reading | " + model+" -> ",  name_file)
        name_txt_part=name_file.split("/")
        f=open(path_output+name_txt_part[-1][10:20]+"/"+name_txt_part[-1][10:20]+".txt", "a")
    if key_gz==1:
        desc_gz(name_file)
        part_post=read_binaryFile_fortran(name_file[:-3], type_file, x_left_lower_corner,y_left_lower_corner,
               x_right_upper_corner,y_right_upper_corner)
        cmd_rm= "rm -rf "+name_file[:-3]
        os.system(cmd_rm)
    else:
        part_post=read_binaryFile_fortran(name_file, type_file, x_left_lower_corner,y_left_lower_corner,
               x_right_upper_corner,y_right_upper_corner)

    lat_masked, lon_masked, mascara=load_mask_grid_NR(file_mask, name_mascara,name_variable_lon, name_variable_lat, type_lon)
    submatrix=determine_id_binary_grid_NR_fortran(part_post, lat_masked.flatten(), lon_masked.flatten(), mascara.flatten(), value_mask)
    submatrix=submatrix[np.argsort(submatrix[:, 0])]

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

    if rank==0:
        print ("NumP => ", len(matrix_i[:,0]))
        f.write("%s %d\n"%("NumP: ",len(matrix_i[:,0])))
    dimX, dimY=matrix_i.shape

    tensor_t=np.ones((len(lista_partposi)-1,dimX ,dimY ))*(-999.9)
    tensor_t[-1,:,:]=matrix_i
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
        for i in range(1, size):
            if i < remainder:
                rank_size = count + 1
            else:
                rank_size = count
            tmp = np.empty((rank_size, final_results.shape[1],final_results.shape[2]), dtype=np.float64)
            comm.Recv(tmp, source=i, tag=14)
            tensor_t[int(i_start[i]):int(i_stop[i]),:,:]=tmp
    comm.Bcast(tensor_t, root=0)
    matrix_result=np.ones((len(tensor_t[:,0,0])-1, len(submatrix[:,0]),4))*(-999.9)
    a2=np.arange(len(tensor_t[:,0,0])-1)

    for i in a2[::-1]:
        matrix=kdif(tensor_t[i,:,:], tensor_t[i+1,:,:],-1.,len(submatrix[:,0]), 5)
        matrix_result[i,:,2]=matrix[:,2]
        matrix_result[i,:,1]=matrix[:,1]
        matrix_result[i,:,0]=matrix[:,0]
        matrix_result[i,:,3]=matrix[:,3]

    if filter_value!=0:
        for i in range(len(matrix_result[:, 0, 0])):
            output_matrix, counter_part=Filter_Part(matrix_result[i, :, :],matrix_result[-1, :, :],-1,threshold,len(matrix_result[-1, :, 2]))
            matrix_result[i,:,:]=output_matrix
        if rank==0:
            print ("Number of filtered particles => ", counter_part)
            f.write("%s %d\n"%("Filtered NumP : ",counter_part))
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
    return matrix_result, tensor_t[-1,:,0], matrix_result_por, tensor_t[-1,:,3]

def _forward_dq(lista_partposi ,file_mask, name_mascara,name_variable_lon, name_variable_lat,lat_f, lon_f,rank,size, type_file,
                x_left_lower_corner,y_left_lower_corner, x_right_upper_corner,y_right_upper_corner, model, type_lon, value_mask, key_gz,path_output):
    name_file=lista_partposi[0]
    if rank==0:
        print ("Reading | " + model+" -> ",  name_file)
        name_txt_part=name_file.split("/")
        f=open(path_output+name_txt_part[-1][10:20]+"/"+name_txt_part[-1][10:20]+".txt", "a")
    if key_gz==1:
        desc_gz(name_file)
        part_post=read_binaryFile_fortran(name_file[:-3], type_file,x_left_lower_corner,y_left_lower_corner, x_right_upper_corner,y_right_upper_corner)
        cmd_rm= "rm -rf "+name_file[:-3]
        os.system(cmd_rm)
    else:
        part_post=read_binaryFile_fortran(name_file, type_file,x_left_lower_corner,y_left_lower_corner, x_right_upper_corner,y_right_upper_corner)

    lat_masked, lon_masked , mascara=load_mask_grid_NR(file_mask, name_mascara,name_variable_lon, name_variable_lat, type_lon)
    submatrix=determine_id_binary_grid_NR_fortran(part_post, lat_masked.flatten(), lon_masked.flatten(), mascara.flatten(), value_mask)
    submatrix=submatrix[np.argsort(submatrix[:, 0])]
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
        for i in range(1, size):
            if i < remainder:
                rank_size = count + 1
            else:
                rank_size = count
            tmp = np.empty((rank_size, final_results.shape[1],final_results.shape[2]), dtype=np.float64)
            comm.Recv(tmp, source=i, tag=14)
            tensor_t[int(i_start[i])+1:int(i_stop[i])+1,:,:]=tmp
    comm.Bcast(tensor_t, root=0)
    if rank==0:
        print ("NumP => ", len(submatrix[:,0]))
        f.write("%s %d\n"%("NumP: ",len(submatrix[:,0])))
    matrix_result=np.ones((len(tensor_t[:,0,0])-1, len(submatrix[:,0]),4))*(-999.9)
    a2=np.arange(len(tensor_t[:,0,0])-1)
    for i in a2:
        matrix=kdif(tensor_t[i,:,:], tensor_t[i+1,:,:],1.,len(submatrix[:,0]), 5)
        matrix_result[i,:,2]=matrix[:,2]
        matrix_result[i,:,1]=matrix[:,1]
        matrix_result[i,:,0]=matrix[:,0]
        matrix_result[i,:,3]=matrix[:,3]
    return matrix_result, tensor_t[0,:,0], tensor_t[0,:,3]

class InputNotInRangeError(Exception):

    def __init__(self, message):
        self.message = message
        super().__init__(self.message)

def to_check_params(paso,type_file,numPdx , numPdy, type_lon, method, resolution, cant_plazo, file_mask):

    if paso==1 or paso==-1:
        pass
    else:
        raise InputNotInRangeError("Only forward (mode = 1) and backward mode (mode = - 1) are allowed")

    if type_file==1 or type_file==2:
        pass
    else:
        raise InputNotInRangeError("Only FLEXPART-WRF (type_file = 1) and FLEXPART model (type_file = 2) are allowed")

    if type_lon==1 or type_lon==2:
        pass
    else:
        raise InputNotInRangeError("Only (type_lon = 1) and (type_lon = 2) are allowed. Review the Readme.txt")

    if method==1 or method==2:
        pass
    else:
        raise  InputNotInRangeError("Only Sthol and James ( method = 1) and Sodemman methodology (method = 2) are allowed")

    if numPdx <=0  or  numPdy<=0:
        raise InputNotInRangeError(" numPdX  and  numPdY  must be greater than zero ")
    if resolution<=0:
        raise InputNotInRangeError(" resolution must be greater than zero ")
    if not cant_plazo%4==0:
        raise InputNotInRangeError("Step_numbers must be divisible by 4")
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

def main_process(path, paso, comm, size, rank, resolution, numPdX, numPdY, cant_plazo, year, month,
         day, hour, time, path_output,file_mask, name_mascara,name_variable_lon, name_variable_lat,x_lower_left,y_lower_left, type_file,masa,numP, x_left_lower_corner,y_left_lower_corner, x_right_upper_corner,y_right_upper_corner, model, type_lon, method, threshold, filter_value, output_txt, output_npy,output_nc, value_mask, key_gz, save_position_part):

    create_directory(path_output)
    t=np.arange(0,cant_plazo+4, 4)
    lat,lon= grid_point (resolution, numPdX, numPdY,x_lower_left,y_lower_left)
    area=calc_A(resolution, lat, lon)
    dimX, dimY=lat.shape
    lat_plot,lon_plot= grid_plot_final(lat, lon)
    function_proof(lat_plot, lon_plot)
    np.savetxt(path_output+"lat_plot.txt", lat_plot)
    np.savetxt(path_output+"lon_plot.txt", lon_plot)
    folder=year+month+day+hour
    create_directory(path_output+folder)
    fecha= year+"-"+month+"-"+day+" "+hour+":"+"00"+":"+"00"
    lista_partposi=generate_file(paso, cant_plazo, fecha, path, key_gz)

    for i in lista_partposi:
        my_file=Path(i)
        if not my_file.is_file():
            if rank==0:
                print (i +"-> not exits" )
            sys.exit()

    if paso==1:
        matrix_result, id_part, q_ini=_forward_dq(lista_partposi, file_mask, name_mascara,name_variable_lon, name_variable_lat,lat, lon,rank,size, type_file,x_left_lower_corner,y_left_lower_corner, x_right_upper_corner,y_right_upper_corner, model, type_lon, value_mask, key_gz, path_output)
        matrix_save=np.empty((matrix_result.shape[0],matrix_result.shape[1],6))
        matrix_save[:,:,:-2]=matrix_result
        matrix_save[:,:,-2]=id_part
        matrix_save[:,:,-1]=q_ini
        if save_position_part=="True":
            np.save(path_output+folder+"/"+"forw_"+folder+"_dq_dt.npy", matrix_save)

    if paso ==-1:
        matrix_result, id_part, matrix_result_por,q_ini=_backward_dq(lista_partposi, file_mask, name_mascara,name_variable_lon, name_variable_lat,lat, lon,rank,size, type_file,x_left_lower_corner,y_left_lower_corner,x_right_upper_corner,y_right_upper_corner, model, type_lon, method,threshold,filter_value, value_mask, key_gz, path_output)
        matrix_save=np.empty((matrix_result.shape[0],matrix_result.shape[1],6))
        matrix_save[:,:,:-2]=matrix_result
        matrix_save[:,:,-2]=id_part
        matrix_save[:,:,-1]=q_ini
        if save_position_part=="True":
            np.save(path_output+folder+"/"+"back_"+folder+"_dq_dt.npy", matrix_save)

    if rank==0:
        ndf=np.arange(1,(cant_plazo/4)+1,1)
        ndb=np.arange(cant_plazo/4,0,-1)
        density=masa/numP
        array_day=np.empty((len(t)-1,dimX-1, dimY-1))
        array_day_por=np.empty((len(t)-1,dimX-1, dimY-1))

        date_save=[]
        date_save.append(convert_date_to_ordinal(int(year), int(month), int(day), int(hour), 0, 0))

        for ii in range(len(t)-1):
            final_results=np.array(K_dq(matrix_result[t[ii]:t[ii+1],:,:-1],lon,lat,numPdY,numPdX,len(matrix_result[t[ii]:t[ii+1],0,0]),len(matrix_result[0,:,0])),dtype=np.float64)
            E_P=final_results*(density)/area
            if method==2:
                POR=np.array(K_dq_POR(matrix_result_por[t[ii]:t[ii+1],:,:],lon,lat,numPdY,numPdX,len(matrix_result_por[t[ii]:t[ii+1],0,0]),len(matrix_result_por[0,:,0])),dtype=np.float64)
            if paso ==-1:
                array_day[int (ndb[ii]-1), :,:]=E_P
                date_back=time_calc_day(year+"-"+month+"-"+day+" "+hour+":00:00", ndb[ii]*(-1))
                
                date_save.append(convert_date_to_ordinal(int(date_back.year), int(date_back.month), int(date_back.day), int(date_back.hour), 0, 0))
                if output_txt!=0:
                    np.savetxt(path_output+folder+"/day_"+str(ndb[ii])+"_"+folder+".txt", E_P)
                if method==2:
                    array_day_por[int (ndb[ii]-1), :,:]=POR
                    if output_txt!=0:
                        np.savetxt(path_output+folder+"/day_POR_"+str(ndb[ii])+"_"+folder+".txt", POR)
            if paso ==1:
                array_day[int(ndf[ii]-1), :,:]=E_P
                date_forw=time_calc_day(year+"-"+month+"-"+day+" "+hour+":00:00", ndf[ii])
                date_save.append(convert_date_to_ordinal(int(date_forw.year), int(date_forw.month), int(date_forw.day), int(date_forw.hour), 0, 0))
                if output_txt!=0:
                    np.savetxt(path_output+folder+"/day_"+str(ndf[ii])+"_"+folder+".txt", E_P)
        if paso ==-1:
            if output_nc!=0:
                function(lat_plot[:,0], lon_plot[0,:], array_day, "back_"+folder, path_output+folder+"/", "E_P", "mm/day", date_save[1:])
            if output_npy!=0:
                np.save(path_output+folder+"/"+"back_"+folder+".npy", array_day)
            if method==2:
                if output_nc!=0:
                    function(lat_plot[:,0], lon_plot[0,:], array_day_por, "POR_back_"+folder, path_output+folder+"/","P_C", "%", date_save[1:])
                if output_npy!=0:
                    np.save(path_output+folder+"/"+"POR_back_"+folder+".npy", array_day_por)
        if paso == 1:
            if output_nc!=0:
                function(lat_plot[:,0], lon_plot[0,:], array_day, "forw_"+folder, path_output+folder+"/", "E_P", "mm/day", date_save[:-1])
            if output_npy!=0:
                np.save(path_output+folder+"/"+"forw_"+folder+".npy", array_day)

if __name__=='__main__':

    input_file=sys.argv[1]
    params=load_parameters(input_file)
    path=params[0]
    path_output=params[1]
    paso =params[2]
    masa=params[3]
    numP=params[4]
    type_file=params[5]
    resolution=params[6]
    numPdX=params[7]
    numPdY=params[8]
    x_lower_left=params[9]
    y_lower_left=params[10]
    cant_plazo=params[11]
    year=params[12]
    month=params[13]
    day=params[14]
    hour=params[15]
    ndias=params[16]
    file_mask=params[17]
    name_mascara=params[18]
    name_variable_lat=params[19]
    name_variable_lon=params[20]
    x_left_lower_corner=params[21]
    y_left_lower_corner=params[22]
    x_right_upper_corner=params[23]
    y_right_upper_corner=params[24]
    model=params[25]
    type_lon=params[26]
    method=params[27]
    threshold=params[28]
    filter_value=params[29]
    output_txt = params[30]
    output_npy = params[31]
    output_nc = params[32]
    name_target_region = params[33]
    value_mask = params[34]
    key_gz = params[35]
    save_position_part=params[36]

    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()

    fecha=year+"-"+month+"-"+day+" "+hour+":00:00"
    list_year, list_month, list_day, list_hour=generate_fecha_simulation(ndias, fecha)
    to_check_params(paso,type_file,numPdX , numPdY, type_lon, method, resolution, cant_plazo,file_mask)

    for year, month, day, hour in zip(list_year, list_month, list_day, list_hour):
            if rank==0:
                print ("                                                                                          ")
                print ('--------------------------------TROVA has started-----------------------------------------')
                print (' *****************************************************************************************')
                print (" *                    EPhysLab (Environmental Physics Laboratory)                        *")
                print (" *                        TRansport Of water VApor (TROVA)                               *")
                print (" *                             version 1.1 (28-09-2022)                                  *")
                TROVA_LOGO()
                print (" *                            Edificio Campus da Auga                                    *")
                print (" *                               University of Vigo                                      *")
                print (" *                                ephyslab.uvigo.es                                      *")
                print (" *  contact: jose.carlos.fernandez.alvarez@uvigo.es, albenis.perez.alarcon@uvigo.es      *")
                print (" *****************************************************************************************")
                print ('------------------------------------------------------------------------------------------')
                print ("                                                                                          ")
                print ("-------------------------------------Information of the run-------------------------------")
                if method==1:
                    print ("You are using methodology of Stohl and James (2005) (https://doi.org/10.1175/1525-7541(2004)005<0656:ALAOTA>2.0.CO;2)")
                if method==2:
                    print ("You are using methodology of Sodemann et al. (2008) (http://dx.doi.org/10.1002/2017JD027543)")
                print ('Target region ->  '+name_target_region)
                print ('CPU used ->   '+ str(size))
                print ('                                                                                          ')
                print ('------------------------------------------------------------------------------------------')
            start_time = time()
            main_process(path,paso, comm, size, rank, resolution, numPdX, numPdY, cant_plazo,
               year, month, day, hour, start_time, path_output,file_mask,name_mascara,name_variable_lon, name_variable_lat,x_lower_left,y_lower_left, type_file,masa,numP,x_left_lower_corner,y_left_lower_corner,x_right_upper_corner,y_right_upper_corner, model, type_lon, method,threshold, filter_value,output_txt,output_npy,output_nc, value_mask, key_gz, save_position_part)
            elapsed_time = time() - start_time
            if rank==0:
                print ('                                                                                          ')
                print ('-------------------------------TROVA has ended--------------------------------------------')
                print ("Congratulations, the run has been successfully completed: ---> Run time: %.2f seconds." % np.round(elapsed_time, 2))
                print ('------------------------------------------------------------------------------------------')
