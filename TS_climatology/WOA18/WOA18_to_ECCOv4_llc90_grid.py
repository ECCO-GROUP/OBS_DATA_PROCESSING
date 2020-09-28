#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul  5 13:26:55 2020

@author: ifenty
"""



import sys

sys.path.append('/home/ifenty/ECCOv4-py')
import ecco_v4_py as ecco
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from pathlib import Path
import importlib
import netCDF4 as nc4
import pyresample as pr

import scipy as scipy
# requires seawater package
# https://pythonhosted.org/seawater
import seawater as sw

import sys
sys.path.append('/home/ifenty/ECCOv4-py/')
import ecco_v4_py as ecco
import json
import glob
import numpy as np
import collections
import numpy as np
import xarray as xr
import xmitgcm as xm
import os
from pathlib import Path
import netCDF4 as nc4

# requires ecco-cloud package
from datetime import datetime

sys.path.append('/home/ifenty/ECCO-GROUP/ECCO-ACCESS/ecco-cloud-utils')
import ecco_cloud_utils as ec

#%%


#%%
def save_to_disk(data_DS,
                 output_filename,
                 binary_fill_value, netcdf_fill_value,
                 netcdf_output_dir, binary_output_dir, binary_output_dtype,
                 model_grid_type, save_binary=True, save_netcdf=True):

    if save_binary:
        print('saving binary')

        # define binary file output filetype
        dt_out = np.dtype(binary_output_dtype)

        # create directory
        binary_output_dir.mkdir(exist_ok=True)

        # define binary output filename
        binary_output_filename = binary_output_dir / output_filename


        for data_var in list(data_DS.data_vars):
            data_DA = data_DS[data_var]
            print(data_DA.name)
            
            # replace nans with the binary fill value (something like -9999)
            tmp_field = np.where(np.isnan(data_DA.values),
                                  binary_fill_value, data_DA.values)
    
            # SAVE FLAT BINARY
            # loop through each record of the year, save binary fields one at a time
            # appending each record as we go
            fd1 = open(str(binary_output_filename), 'wb')
            fd1 = open(str(binary_output_filename), 'ab')
        
            # if we have an llc grid, then we have to reform to compact
            if model_grid_type == 'llc':
                tmp_field = ecco.llc_tiles_to_compact(
                    tmp_field, less_output=True)

            # otherwise assume grid is x,y (2 dimensions)
            #elif model_grid_type == 'latlon':
             
            # make sure we have something to save...
            if len(tmp_field) > 0:
                # if this is the first record, create new binary file
                tmp_field.astype(dt_out).tofile(fd1)
    
            # close the file at the end of the operation
            fd1.close()

    if save_netcdf:
        print('saving netcdf record')

        # create directory
        netcdf_output_dir.mkdir(exist_ok=True)

        # define netcdf output filename
        netcdf_output_filename = netcdf_output_dir / \
            Path(output_filename + '.nc')

        # SAVE NETCDF
        # replace the binary fill value (-9999) with the netcdf fill value
        # which is much more interesting

        for data_var in list(data_DS.data_vars):
            data_DA = data_DS[data_var]
            print(data_DA.name)
            
            # replace nans with the binary fill value (something like -9999)
            data_DA.values = \
                np.where(np.isnan(data_DA.values),
                         netcdf_fill_value, data_DA.values)

            encoding_each = {'zlib': True,
                             'complevel': 1,
                             'fletcher32': True,
                             '_FillValue': netcdf_fill_value}


            encoding_coords = {'zlib': True,
                             'complevel': 1,
                             'fletcher32': True,
                             '_FillValue': None}

        encoding_dv = {var: encoding_each for var in data_DS.data_vars}
        encoding_coords = {var: encoding_coords for var in data_DS.coords}

        encoding= {**encoding_coords, **encoding_dv}
        print(encoding)
        
        # the actual saving (so easy with xarray!)
        data_DS.to_netcdf(netcdf_output_filename,  encoding=encoding)
        data_DS.close()
        
        
#####################
# Define precision of output files, float32 is standard
array_precision = np.float32

# Define fill values for binary and netcdf
if array_precision == np.float32:
    binary_output_dtype = '>f4'
    netcdf_fill_value = nc4.default_fillvals['f4']

elif array_precision == np.float64:
    binary_output_dtype = '>f8'
    netcdf_fill_value = nc4.default_fillvals['f8']

# ECCO always uses -9999 for missing data.
binary_fill_value = -9999


#%%

# directories with the 1995-2004 and 2005-2017 monthly climatology files
source_dir_1995_2004 = Path('/mnt/intraid/ian3/ifenty/data/observations/TS_climatology/WOA_2018/1995-2004/monthly')
source_dir_2005_2017 = Path('/mnt/intraid/ian3/ifenty/data/observations/TS_climatology/WOA_2018/2005-2017/monthly')

source_dir_1995_2004_ann = Path('/mnt/intraid/ian3/ifenty/data/observations/TS_climatology/WOA_2018/1995-2004/annual')
source_dir_2005_2017_ann = Path('/mnt/intraid/ian3/ifenty/data/observations/TS_climatology/WOA_2018/2005-2017/annual')

combined_output_dir_1995_2017 = Path('/mnt/intraid/ian3/ifenty/data/observations/TS_climatology/WOA_2018/1995-2017_combined/')
combined_output_dir_1995_2017.mkdir()

#%% get all filenames (sorted)
t_files_95A4 = np.sort(list(source_dir_1995_2004.glob('*95A4_t*')))
t_files_95A4

s_files_95A4 = np.sort(list(source_dir_1995_2004.glob('*95A4_s*')))
s_files_95A4

t_files_0517 = np.sort(list(source_dir_2005_2017.glob('*A5B7_t*')))
t_files_0517

s_files_0517 = np.sort(list(source_dir_2005_2017.glob('*A5B7_s*')))
s_files_0517

#%%
t_files_95A4_ann = np.sort(list(source_dir_1995_2004_ann.glob('*95A4_t*')))
print(t_files_95A4_ann)

s_files_95A4_ann = np.sort(list(source_dir_1995_2004_ann.glob('*95A4_s*')))
print(s_files_95A4_ann)

t_files_0517_ann = np.sort(list(source_dir_2005_2017_ann.glob('*A5B7_t*')))
print(t_files_0517_ann)

s_files_0517_ann = np.sort(list(source_dir_2005_2017_ann.glob('*A5B7_s*')))
print(s_files_0517_ann)


# load t and s climatologies
tfs_95A4_ann = xr.open_dataset(t_files_95A4_ann[0], decode_times=False)
tfs_0517_ann = xr.open_dataset(t_files_0517_ann[0], decode_times=False)

sfs_95A4_ann = xr.open_dataset(s_files_95A4_ann[0], decode_times=False)
sfs_0517_ann = xr.open_dataset(s_files_0517_ann[0], decode_times=False)

# load vertical levels

clim_Z_bnds_ann = sfs_95A4_ann.depth_bnds.values
clim_Z_centers_ann  = []

nz = len(sfs_95A4_ann.depth)

for k in range(len(sfs_95A4_ann.depth)):
    clim_Z_centers_ann.append(np.mean(clim_Z_bnds_ann[k]))
    
clim_Z_centers_ann = -np.array(clim_Z_centers_ann)

# Merge the two decades of clim

tfs_9517_ann = 0.5*(tfs_95A4_ann.t_an.values[0,:]+ tfs_0517_ann.t_an.values[0,:])
t_clim_1995_2017_ann_DA = xr.DataArray(tfs_9517_ann, dims=['Z','lat','lon'])
t_clim_1995_2017_ann_DA=t_clim_1995_2017_ann_DA.assign_coords({'lat' : tfs_95A4_ann.lat})
t_clim_1995_2017_ann_DA=t_clim_1995_2017_ann_DA.assign_coords({'lon' : tfs_95A4_ann.lon})
t_clim_1995_2017_ann_DA=t_clim_1995_2017_ann_DA.assign_coords({'Z' : clim_Z_centers_ann})


# Merge the two decades of clim

sfs_9517_ann = 0.5*(sfs_95A4_ann.s_an.values[0,:] + sfs_0517_ann.s_an.values[0,:])
s_clim_1995_2017_ann_DA = xr.DataArray(sfs_9517_ann, dims=['Z','lat','lon'])
s_clim_1995_2017_ann_DA=s_clim_1995_2017_ann_DA.assign_coords({'lat' : tfs_95A4_ann.lat})
s_clim_1995_2017_ann_DA=s_clim_1995_2017_ann_DA.assign_coords({'lon' : tfs_95A4_ann.lon})
s_clim_1995_2017_ann_DA=s_clim_1995_2017_ann_DA.assign_coords({'Z' : clim_Z_centers_ann})

#%%

# Load netcdf files, combine into single objects using xr.concat

tfs = []
for tf in t_files_95A4:
    tfs.append(xr.open_dataset(tf, decode_times=False))
    
tfs_95A4 = xr.concat((tfs), dim='time')

tfs = []
for tf in t_files_0517:
    tfs.append(xr.open_dataset(tf, decode_times=False))
   
tfs_0517 = xr.concat((tfs), dim='time')

sfs = []
for sf in s_files_95A4:
    sfs.append(xr.open_dataset(sf, decode_times=False))
    
sfs_95A4 = xr.concat((sfs), dim='time')

sfs = []
for tf in s_files_0517:
    sfs.append(xr.open_dataset(sf, decode_times=False))
   
sfs_0517 = xr.concat((sfs), dim='time')



#%%
# FIND Z COORDS OF ORIGINAL DATA, MIDDLE OF THE DEPTH_BNDS

clim_Z_bnds = sfs_95A4.depth_bnds[0,:].values
clim_Z_centers  = []

for k in range(clim_Z_bnds.shape[0]):
    clim_Z_centers.append(np.mean(clim_Z_bnds[k]))
    
clim_Z_centers = -np.array(clim_Z_centers)


#%%
# calculate the mean of the two time periods
tmp = 0.5 *( tfs_95A4.t_an.values + tfs_0517.t_an.values)

# make a new temperature dataset this time with the records organized with
# the time dimension 'month' corresponding to 1..12 month of year
t_clim_1995_2017_DA = xr.DataArray(tmp, dims=['month','Z','lat','lon'])
t_clim_1995_2017_DA=t_clim_1995_2017_DA.assign_coords({'month' : range(1,13)})
t_clim_1995_2017_DA=t_clim_1995_2017_DA.assign_coords({'lat' : tfs_0517.lat})
t_clim_1995_2017_DA=t_clim_1995_2017_DA.assign_coords({'lon' : tfs_0517.lon})
t_clim_1995_2017_DA=t_clim_1995_2017_DA.assign_coords({'Z' : clim_Z_centers})

# repeat for salinity
tmp = 0.5 *( sfs_95A4.s_an.values + sfs_0517.s_an.values)

s_clim_1995_2017_DA = xr.DataArray(tmp, dims=['month','Z','lat','lon'])
s_clim_1995_2017_DA=s_clim_1995_2017_DA.assign_coords({'month' : range(1,13)})
s_clim_1995_2017_DA=s_clim_1995_2017_DA.assign_coords({'lat' : tfs_0517.lat})
s_clim_1995_2017_DA=s_clim_1995_2017_DA.assign_coords({'lon' : tfs_0517.lon})
s_clim_1995_2017_DA=s_clim_1995_2017_DA.assign_coords({'Z' : clim_Z_centers})

#%%
# CALCULATE PRESSURE from Z

# seawater.eos80.pres(depth, lat)
#    Calculates pressure in dbars from depth in meters.
#    Parameters:	
#
#    depth : array_like
#        depth [meters]
#    lat : array_like
#        latitude in decimal degrees north [-90..+90]
#
#    Returns:	
#    p : array_like
#        pressure [db]


#%%
# ADD PRESSURE TO DATA ARRAY
pressure = np.zeros((len(s_clim_1995_2017_DA.Z.values),\
                     len(s_clim_1995_2017_DA.lat.values)))

for zi, z in enumerate(s_clim_1995_2017_DA.Z.values):
    pressure[zi,:] = sw.eos80.pres(-z, s_clim_1995_2017_DA.lat.values)
    
pressure_3D = (np.tile(pressure.T, (360,1,1))).T;

s_clim_1995_2017_DA = s_clim_1995_2017_DA.assign_coords({'pressure' : (('Z', 'lat','lon'), pressure_3D)})
t_clim_1995_2017_DA = t_clim_1995_2017_DA.assign_coords({'pressure' : (('Z', 'lat','lon'), pressure_3D )})


#%%

pressure_ann = np.zeros((len(s_clim_1995_2017_ann_DA.Z.values),\
                     len(s_clim_1995_2017_ann_DA.lat.values)))

for zi, z in enumerate(s_clim_1995_2017_ann_DA.Z.values):
    pressure_ann[zi,:] = sw.eos80.pres(-z, s_clim_1995_2017_ann_DA.lat.values)
    
pressure_3D_ann = (np.tile(pressure_ann.T, (360,1,1))).T;

s_clim_1995_2017_ann_DA = s_clim_1995_2017_ann_DA.assign_coords({'pressure' : (('Z', 'lat','lon'), \
                                                                               pressure_3D_ann)})
t_clim_1995_2017_ann_DA = t_clim_1995_2017_ann_DA.assign_coords({'pressure' : (('Z', 'lat','lon'),\
                                                                               pressure_3D_ann )})


#%%  
# CONVERT IN SITU TEMPERATURE TO POTENTIAL TEMPERATURE

# seawater.eos80.ptmp(s, t, p, pr=0)
#    Calculates potential temperature as per UNESCO 1983 report.
#    Parameters:	
#
#    s(p) : array_like
#        salinity [psu (PSS-78)]
#    t(p) : array_like
#        temperature [? (ITS-90)]
#    p : array_like
#        pressure [db].
#    pr : array_like
#        reference pressure [db], default = 0
#
#    Returns:	
#
#    pt : array_like
#        potential temperature relative to PR [? (ITS-90)]

pott_clim_1995_2017_DA = t_clim_1995_2017_DA.copy(deep=True)
pott_clim_1995_2017_DA.values = pott_clim_1995_2017_DA.values*np.nan
pott_clim_1995_2017_DA.name = 'Potential Temperature'

# loop through each month, calculate potential temperature
for m in range(12):
    print(m)
    pott_clim_1995_2017_DA.values[m,:] = sw.eos80.ptmp(s_clim_1995_2017_DA.values[m,:,:], \
                                                       t_clim_1995_2017_DA.values[m,:,:], \
                                                       t_clim_1995_2017_DA.pressure.values)
    

#%% annual mean
pott_clim_1995_2017_ann_DA = t_clim_1995_2017_ann_DA.copy(deep=True)
pott_clim_1995_2017_ann_DA.values = pott_clim_1995_2017_ann_DA.values*np.nan
pott_clim_1995_2017_ann_DA.name = 'Potential Temperature'


pott_clim_1995_2017_ann_DA.values = sw.eos80.ptmp(s_clim_1995_2017_ann_DA.values[:,:], \
                                                       t_clim_1995_2017_ann_DA.values[:,:], \
                                                       t_clim_1995_2017_ann_DA.pressure.values)
    
    
    
#%%
# sanity check, potential temperature should be less than in situ temperature
# at depth.
    
cmin=5;
cmax=25
x = t_clim_1995_2017_ann_DA.lon.values
y = t_clim_1995_2017_ann_DA.pressure.values[:,0,0]

plt.figure(1,clear=True, figsize=(10,15));
plt.subplot(311)
#plt.imshow(x,y, t_clim_1995_2017_ann_DA[:,90,:].values,origin='upper',vmin=cmin, vmax=cmax)
t_clim_1995_2017_ann_DA[:,90,:].plot(vmin=0,vmax=10)
#vmin=cmin,vmax=cmax)
#plt.colorbar()
plt.title('in situ temperature')

plt.subplot(312)
pott_clim_1995_2017_ann_DA[:,90,:].plot(vmin=0,vmax=10)
plt.title('potential temperature')


cmin=-.6;cmax=.6
plt.subplot(313)
(pott_clim_1995_2017_ann_DA[:,90,:] - t_clim_1995_2017_ann_DA[:,90,:]).plot(vmin=cmin, vmax=cmax, cmap='bwr')
plt.title('potential - in situ temperature')
plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=.5)
plt.suptitle('T and potential T along equator')



#%% save clim fields on original grid to this point

print(combined_output_dir_1995_2017)
combined_output_dir_1995_2017.mkdir()
#%%mkdir()

pott_clim_1995_2017_ann_DA.to_netcdf(combined_output_dir_1995_2017 / '1995-2017_potential_T_annual.nc')
pott_clim_1995_2017_DA.to_netcdf(combined_output_dir_1995_2017 / '1995-2017_potential_T.nc')


s_clim_1995_2017_ann_DA.to_netcdf(combined_output_dir_1995_2017 / '1995-2017_S_annual.nc')
s_clim_1995_2017_DA.to_netcdf(combined_output_dir_1995_2017 / '1995-2017_S.nc')

#%%
# DEFINE 
model_grid_dir = Path('/home/ifenty/data/grids/ecco_pipline_grids/')
model_grid_filename  = 'ECCO_llc90_demo.nc'
model_grid = xr.open_dataset(model_grid_dir / model_grid_filename)

model_grid_dir_2 = Path('/home/ifenty/ECCOv4-release/Release4/nctiles_grid')
model_grid_2 = xr.open_dataset(model_grid_dir_2 / 'ECCO-GRID.nc')
model_grid['Z'] = model_grid_2.Z.values


#%%
model_grid_id   = model_grid.name
model_grid_type =  model_grid.type

effective_grid_radius = model_grid.effective_grid_radius


#%%
# data names
product_name = 'World Ocean Atlas 2018 : sea_water_salinity 1.00 degree'
#product_source = 'https://www.nodc.noaa.gov/OC5/woa18/'

product_source = 'woa18'

# output parameters
# output folder is specified with data idenifier
mapping_time =  datetime.now().strftime("%Y%m%dT%H%M%S")
netcdf_output_dir = Path(f'/mnt/intraid/ian3/ifenty/data/observations/TS_climatology/WOA_2018/data_output_{product_source}/mapped_to_' + model_grid_id + '/' + mapping_time + '/netcdf')
binary_output_dir = Path(f'/mnt/intraid/ian3/ifenty/data/observations/TS_climatology/WOA_2018/data_output_{product_source}/mapped_to_' + model_grid_id + '/' + mapping_time + '/binary')

# Define precision of output files, float32 is standard
array_precision = np.float32

# Define fill values for binary and netcdf
if array_precision == np.float32:
    binary_output_dtype = '>f4'
    netcdf_fill_value = nc4.default_fillvals['f4']

elif array_precision == np.float64:
    binary_output_dtype = '>f8'
    netcdf_fill_value = nc4.default_fillvals['f8']

# ECCO always uses -9999 for missing data.
binary_fill_value = -9999
           
#%%
# fields to process
data_field_T   = {'name': 'THETA',
                'long_name':'Objectively analyzed mean fields for sea_water_potential_temperature converted from woa18 T and S.',
                'standard_name': 'sea_water_potential_temperature',
                'units':'degrees_celsius',
                'comments':'potential T converted from WOA18 in situ T and S using seawater.eos80'}

# fields to process
data_field_S   = {'name': 'SALT' ,
                'long_name':'Objectively analyzed mean fields for sea_water_salinity from woa18.',
                'standard_name':'sea_water_salinity',
                'units': sfs_95A4.s_an.units}

# setup output attributes
new_data_attr_T = {'original_dataset_title':'World Ocean Atlas 2018 : sea_water_temperature January 1995-2004 1.00 degree + 2005-2017 1.00 degree',
                 'original_dataset_url':'https://www.nodc.noaa.gov/OC5/woa18/',
                 'original_dataset_reference':tfs_95A4.references,
                 'original_dataset_product_id':'woa18',
                 'comments':'potential T converted from woa in situ T and S using EOS80'}


# setup output attributes
new_data_attr_S = {'original_dataset_title':'World Ocean Atlas 2018 : sea_water_salinity_January 1995-2004 1.00 degree + 2005-2017 1.00 degree',
                 'original_dataset_url':'https://www.nodc.noaa.gov/OC5/woa18/',
                 'original_dataset_reference':sfs_95A4.references,
                 'original_dataset_product_id':'woa18'}

# data grid information
data_res = 1 
data_max_lat = 89.5 
area_extent = [-180, 90, 180, -90]
dims = [360, 180]

# Grid projection information
proj_info = {'area_id':'longlat',
             'area_name':'Plate Carree',
             'proj_id':'EPSG:4326',
             'proj4_args':'+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'}    

#%%

   #######################################################
## BEGIN GRID PRODUCT                                ##

source_grid_min_L, source_grid_max_L, source_grid, \
 data_grid_lons, data_grid_lats = ec.generalized_grid_product(product_name,
                                                                data_res,
                                                                data_max_lat,
                                                                area_extent,
                                                                dims,
                                                                proj_info)

plt.figure(num=3,clear=True,figsize=(10,15));
plt.subplot(311);plt.imshow(data_grid_lons, origin='lower');plt.colorbar()
plt.title('lons')
plt.subplot(312);plt.imshow(data_grid_lats, origin='lower');plt.colorbar()
plt.title('lats')
plt.subplot(313);plt.imshow(pott_clim_1995_2017_DA.values[0,0,:], origin='lower');plt.colorbar()
plt.title('potential clim surface first time level')
plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=.5)
plt.suptitle('WOA on WOA grid')

#%%
# make output directories
netcdf_output_dir.mkdir(exist_ok=True, parents=True)
binary_output_dir.mkdir(exist_ok=True, parents=True)

# load the model grid
model_grid = xr.open_dataset(model_grid_dir / model_grid_filename)
model_grid  = model_grid.reset_coords()
  
# Define the 'swath' (in the terminology of the pyresample module)
# as the lats/lon pairs of the model grid
# The routine needs the lats and lons to be one-dimensional vectors.
target_grid  = \
    pr.geometry.SwathDefinition(lons=model_grid.XC.values.ravel(), 
                                lats=model_grid.YC.values.ravel())
   
target_grid_radius = model_grid.effective_grid_radius.values.ravel()

# Compute the mapping between the data and model grid
source_indices_within_target_radius_i,\
num_source_indices_within_target_radius_i,\
nearest_source_index_to_target_index_i = \
    ec.find_mappings_from_source_to_target(source_grid,\
                                           target_grid,\
                                           target_grid_radius, \
                                           source_grid_min_L, \
                                           source_grid_max_L)
    
#%%
# MAP THE CLIMATOLOGY T AND S TO MODEL HORIZONTAL GRID
# NOTE, RESULTING ARRAY STILL HAS VERTICAL DIMENSION OF WOA 18 GRID
    
pott_model_grid_xy = np.zeros(([12, 57] + list(model_grid.XC.shape)))
s_model_grid_xy = np.zeros(([12, 57] + list(model_grid.XC.shape)))

for m in range(12):
    for k in range(57):
        
        # POTENTIAL TEMPERATURE
        orig_data = pott_clim_1995_2017_DA[m,k,:].values
        
        data_model_projection = ec.transform_to_target_grid(source_indices_within_target_radius_i, \
                                                            num_source_indices_within_target_radius_i, \
                                                            nearest_source_index_to_target_index_i, \
                                                            orig_data, model_grid.XC.shape)
            
        pott_model_grid_xy[m,k,:] = data_model_projection
        
        # SALINITY
        orig_data = s_clim_1995_2017_DA[m,k,:].values
        
        data_model_projection = ec.transform_to_target_grid(source_indices_within_target_radius_i, \
                                                            num_source_indices_within_target_radius_i, \
                                                            nearest_source_index_to_target_index_i, \
                                                            orig_data, model_grid.XC.shape)
            
        s_model_grid_xy[m,k,:] = data_model_projection
        
        print(m,k)
    
#%%
# annual mean
      
pott_model_grid_ann_xy = np.zeros(([102] + list(model_grid.XC.shape)))
s_model_grid_ann_xy = np.zeros(([102] + list(model_grid.XC.shape)))
  
# POTENTIAL TEMPERATURE
for k in range(102):
    print(k)
    orig_data = pott_clim_1995_2017_ann_DA[k,:].values
    
    data_model_projection = ec.transform_to_target_grid(source_indices_within_target_radius_i, \
                                                        num_source_indices_within_target_radius_i, \
                                                        nearest_source_index_to_target_index_i, \
                                                        orig_data, model_grid.XC.shape)
    
    pott_model_grid_ann_xy[k,:] = data_model_projection
  
    orig_data = s_clim_1995_2017_ann_DA[k,:].values
    
    data_model_projection = ec.transform_to_target_grid(source_indices_within_target_radius_i, \
                                                        num_source_indices_within_target_radius_i, \
                                                        nearest_source_index_to_target_index_i, \
                                                        orig_data, model_grid.XC.shape)
    
    s_model_grid_ann_xy[k,:] = data_model_projection
        
#%%
if model_grid_type == 'llc':
    plt.close('all')        
    ecco.plot_tiles(pott_model_grid_xy[0,9,:], \
                    layout='latlon', rotate_to_latlon=True, show_tile_labels=False)

    plt.close('all')        
    ecco.plot_tiles(pott_model_grid_xy[0,9,:], \
                    layout='latlon', rotate_to_latlon=True, show_tile_labels=False)

           

#%% save intermediate solution to disk
# WOA mapped to model grid in x and y but not depth
    
    
pott_ann_model_grid_xy_DA = xr.DataArray(pott_model_grid_ann_xy, dims=['Z','tile','j','i'])
pott_ann_model_grid_xy_DA=pott_ann_model_grid_xy_DA.assign_coords({'j' :model_grid.j})
pott_ann_model_grid_xy_DA=pott_ann_model_grid_xy_DA.assign_coords({'i' : model_grid.i})
pott_ann_model_grid_xy_DA=pott_ann_model_grid_xy_DA.assign_coords({'tile' :model_grid.tile})
pott_ann_model_grid_xy_DA=pott_ann_model_grid_xy_DA.assign_coords({'Z' : clim_Z_centers_ann})
pott_ann_model_grid_xy_DA.name = 'annual_clim_potential_T'
 
s_ann_model_grid_xy_DA = xr.DataArray(s_model_grid_ann_xy, dims=['Z','tile','j','i'])
s_ann_model_grid_xy_DA=s_ann_model_grid_xy_DA.assign_coords({'j' :model_grid.j})
s_ann_model_grid_xy_DA=s_ann_model_grid_xy_DA.assign_coords({'i' : model_grid.i})
s_ann_model_grid_xy_DA=s_ann_model_grid_xy_DA.assign_coords({'tile' :model_grid.tile})
s_ann_model_grid_xy_DA=s_ann_model_grid_xy_DA.assign_coords({'Z' : clim_Z_centers_ann})
s_ann_model_grid_xy_DA.name = 'annual_clim_S'


pott_model_grid_xy_DA = xr.DataArray(pott_model_grid_xy, dims=['month','Z','tile','j','i'])
pott_model_grid_xy_DA=pott_model_grid_xy_DA.assign_coords({'j' :model_grid.j})
pott_model_grid_xy_DA=pott_model_grid_xy_DA.assign_coords({'i' : model_grid.i})
pott_model_grid_xy_DA=pott_model_grid_xy_DA.assign_coords({'tile' :model_grid.tile})
pott_model_grid_xy_DA=pott_model_grid_xy_DA.assign_coords({'Z' : clim_Z_centers})
pott_model_grid_xy_DA=pott_model_grid_xy_DA.assign_coords({'month' : range(1,13)})
pott_model_grid_xy_DA.name = 'mon_clim_potential_T'


s_model_grid_xy_DA = xr.DataArray(s_model_grid_xy, dims=['month','Z','tile','j','i'])
s_model_grid_xy_DA=s_model_grid_xy_DA.assign_coords({'j' :model_grid.j})
s_model_grid_xy_DA=s_model_grid_xy_DA.assign_coords({'i' : model_grid.i})
s_model_grid_xy_DA=s_model_grid_xy_DA.assign_coords({'tile' :model_grid.tile})
s_model_grid_xy_DA=s_model_grid_xy_DA.assign_coords({'Z' : clim_Z_centers})
s_model_grid_xy_DA=s_model_grid_xy_DA.assign_coords({'month' : range(1,13)})
s_model_grid_xy_DA.name = "mon_clim_S"

#%%
       
pott_ann_model_grid_xy_DA.to_netcdf(netcdf_output_dir / 'pott_ann_model_grid_xy_DA.nc')
s_ann_model_grid_xy_DA.to_netcdf(netcdf_output_dir / 's_ann_model_grid_xy_DA.nc')

pott_model_grid_xy_DA.to_netcdf(netcdf_output_dir / 'pott_model_grid_xy_DA.nc')
s_model_grid_xy_DA.to_netcdf(netcdf_output_dir / 's_model_grid_xy_DA.nc')



#%%## MASK OUT T AND S VALUES SOUTH OF 40S BECAUSE OF DISCONTINUITIES
tmp = s_ann_model_grid_xy_DA[-1,:]
s_ann_model_grid_xy_DA[-1,:] = np.where(model_grid.YC.values < -40, np.nan, tmp)

tmp = pott_ann_model_grid_xy_DA[-1,:]
pott_ann_model_grid_xy_DA[-1,:] = np.where(model_grid.YC.values < -40, np.nan, tmp)


#%% APPEND ANNUAL CLIM VALUES TO MONTHLY MEAN CLIM STARTING AT THE LEVEL
# WITH 1475-1500 M [MONTHLY] AND 1457-1525M [ANNUAL]
# THAT IS THE LAST LEVEL OF CLIM_Z_BNDS (56) and the 56th level of clim_Z_bnds_ann

s_monthly_plus_ann_model_grid_xy = np.zeros([12, 102] + list(s_model_grid_xy_DA.shape[2:]))
pott_monthly_plus_ann_model_grid_xy = np.zeros([12, 102] + list(pott_model_grid_xy_DA.shape[2:]))

for m in range(12):
    s_monthly_plus_ann_model_grid_xy[m,0:56,:] = s_model_grid_xy_DA[m,0:56,:]
    s_monthly_plus_ann_model_grid_xy[m,56:,:] = s_ann_model_grid_xy_DA[56:,:]
    
    pott_monthly_plus_ann_model_grid_xy[m,0:56,:] = pott_model_grid_xy_DA[m,0:56,:]
    pott_monthly_plus_ann_model_grid_xy[m,56:,:] = pott_ann_model_grid_xy_DA[56:,:]
    
    
#%%
plt.figure(102,clear=True);
plt.subplot(121)
for m in range(12):
     plt.plot(s_monthly_plus_ann_model_grid_xy[m,:,9,30,1],clim_Z_centers_ann,'r.-')

plt.plot(s_model_grid_ann_xy[:,9, 30, 1],clim_Z_centers_ann,'k.')
plt.grid()


plt.subplot(122)
for m in range(12):
     plt.plot(pott_monthly_plus_ann_model_grid_xy[m,:,9,30,1],clim_Z_centers_ann,'r.-')

plt.plot(pott_model_grid_ann_xy[:,9,30,1],clim_Z_centers_ann,'k.')
plt.grid()




#%%
# DO VERTICAL INTERPOLATION TO MODEL GRID
from scipy.interpolate import interp1d

#%%
f = interp1d(clim_Z_centers_ann, s_monthly_plus_ann_model_grid_xy,axis=1,bounds_error=False,\
             fill_value=np.nan)

s_monthly_plus_ann_model_grid_xyz = f(model_grid.Z.values)

print('monthly salinity done...')

f = interp1d(clim_Z_centers_ann, pott_monthly_plus_ann_model_grid_xy,axis=1,bounds_error=False,\
             fill_value=np.nan)

pott_monthly_plus_ann_model_grid_xyz = f(model_grid.Z.values)

print('monthly temperature done...')

#%%
f = interp1d(clim_Z_centers_ann, s_ann_model_grid_xy_DA.values,axis=0,bounds_error=False,\
             fill_value=np.nan)

s_ann_model_grid_xyz = f(model_grid.Z.values)

print('annual salinity done...')

f = interp1d(clim_Z_centers_ann, pott_ann_model_grid_xy_DA.values,axis=0,bounds_error=False,\
             fill_value=np.nan)

pott_ann_model_grid_xyz = f(model_grid.Z.values)

print('annual temperature done...')


#%%

plt.figure(200,clear=True)

x=np.mean(s_monthly_plus_ann_model_grid_xyz, axis=0)
x2=np.mean(s_model_grid_xy, axis=0)

for ii in range(25):

    plt.subplot(5,5,ii+1)

    tile = np.int(random.uniform(0,12))
    i = np.int(random.uniform(0,89))
    j = np.int(random.uniform(0,89))
    
    tmp = s_ann_model_grid_xy_DA[:,tile,j,i]
    
    while np.isnan(tmp[-1]):
        tile = np.int(random.uniform(0,12))
        i = np.int(random.uniform(0,89))
        j = np.int(random.uniform(0,89))
        
        tmp = s_ann_model_grid_xy_DA[:,tile,j,i]

    plt.plot(x[:,tile,j,i],model_grid.Z.values,  'b.-')
    plt.plot(x2[:,tile,j,i], clim_Z_centers,'ro-')
    plt.plot(tmp, clim_Z_centers_ann,'g.-',markersize=0.5)
    plt.title(("%i %i %i ") % (tile, j, i))
    plt.grid()
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=.5, hspace=.5)
    plt.suptitle('S on model + WOA grid')
    
    
#%%
plt.figure(201,clear=True)

x=np.mean(pott_monthly_plus_ann_model_grid_xyz, axis=0)
x2=np.mean(pott_model_grid_xy, axis=0)

for ii in range(25):

    plt.subplot(5,5,ii+1)

    tile = np.int(random.uniform(0,12))
    i = np.int(random.uniform(0,89))
    j = np.int(random.uniform(0,89))
    
    tmp = pott_ann_model_grid_xy_DA[:,tile,j,i]
    
    while np.isnan(tmp[-1]):
        tile = np.int(random.uniform(0,12))
        i = np.int(random.uniform(0,89))
        j = np.int(random.uniform(0,89))
        
        tmp = pott_ann_model_grid_xy_DA[:,tile,j,i]

    plt.plot(x[:,tile,j,i],model_grid.Z.values,  'b.-')
    plt.plot(x2[:,tile,j,i], clim_Z_centers,'ro-')
    plt.plot(tmp, clim_Z_centers_ann,'g.-',markersize=0.5)
    plt.title(("%i %i %i ") % (tile, j, i))
    plt.grid()
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=.5, hspace=.5)
    plt.suptitle('potential T on model + WOA grid')

    
#%%
# MONTHLY S
s_monthly_plus_ann_model_grid_xyz_DA = xr.DataArray(s_monthly_plus_ann_model_grid_xyz, \
                                        dims=['month','k', 'tile','j','i'])
s_monthly_plus_ann_model_grid_xyz_DA = s_monthly_plus_ann_model_grid_xyz_DA.assign_coords({'month' : range(1,13)})
s_monthly_plus_ann_model_grid_xyz_DA = s_monthly_plus_ann_model_grid_xyz_DA.assign_coords({'j' : range(90)})
s_monthly_plus_ann_model_grid_xyz_DA = s_monthly_plus_ann_model_grid_xyz_DA.assign_coords({'i' : range(90)})
s_monthly_plus_ann_model_grid_xyz_DA = s_monthly_plus_ann_model_grid_xyz_DA.assign_coords({'k' : range(50)})
s_monthly_plus_ann_model_grid_xyz_DA = s_monthly_plus_ann_model_grid_xyz_DA.assign_coords({'tile' : range(13)})
s_monthly_plus_ann_model_grid_xyz_DA = s_monthly_plus_ann_model_grid_xyz_DA.assign_coords({'XC' : (('tile', 'j', 'i'), model_grid.XC)})
s_monthly_plus_ann_model_grid_xyz_DA = s_monthly_plus_ann_model_grid_xyz_DA.assign_coords({'YC' : (('tile', 'j', 'i'), model_grid.YC)})
s_monthly_plus_ann_model_grid_xyz_DA = s_monthly_plus_ann_model_grid_xyz_DA.assign_coords({'Z'  : ('k', model_grid.Z)})


s_monthly_plus_ann_model_grid_xyz_DA.name = 'SALT'
s_monthly_plus_ann_model_grid_xyz_DA.attrs = data_field_S

s_monthly_plus_ann_model_grid_xyz_DS = s_monthly_plus_ann_model_grid_xyz_DA.to_dataset()
s_monthly_plus_ann_model_grid_xyz_DS.attrs = new_data_attr_S


# MONTHLY THETA

pott_monthly_plus_ann_model_grid_xyz_DA = xr.DataArray(pott_monthly_plus_ann_model_grid_xyz, \
                                        dims=['month','k', 'tile','j','i'])
pott_monthly_plus_ann_model_grid_xyz_DA = pott_monthly_plus_ann_model_grid_xyz_DA.assign_coords({'month' : range(1,13)})
pott_monthly_plus_ann_model_grid_xyz_DA = pott_monthly_plus_ann_model_grid_xyz_DA.assign_coords({'j' : range(90)})
pott_monthly_plus_ann_model_grid_xyz_DA = pott_monthly_plus_ann_model_grid_xyz_DA.assign_coords({'i' : range(90)})
pott_monthly_plus_ann_model_grid_xyz_DA = pott_monthly_plus_ann_model_grid_xyz_DA.assign_coords({'k' : range(50)})
pott_monthly_plus_ann_model_grid_xyz_DA = pott_monthly_plus_ann_model_grid_xyz_DA.assign_coords({'tile' : range(13)})
pott_monthly_plus_ann_model_grid_xyz_DA = pott_monthly_plus_ann_model_grid_xyz_DA.assign_coords({'XC' : (('tile', 'j', 'i'), model_grid.XC)})
pott_monthly_plus_ann_model_grid_xyz_DA = pott_monthly_plus_ann_model_grid_xyz_DA.assign_coords({'YC' : (('tile', 'j', 'i'), model_grid.YC)})
pott_monthly_plus_ann_model_grid_xyz_DA = pott_monthly_plus_ann_model_grid_xyz_DA.assign_coords({'Z'  : ('k', model_grid.Z)})

pott_monthly_plus_ann_model_grid_xyz_DA.name = 'THETA'
pott_monthly_plus_ann_model_grid_xyz_DA.attrs = data_field_T

pott_monthly_plus_ann_model_grid_xyz_DS = pott_monthly_plus_ann_model_grid_xyz_DA.to_dataset()
pott_monthly_plus_ann_model_grid_xyz_DS.attrs = new_data_attr_T

# ANNUAL S
s_ann_model_grid_xyz_DA = xr.DataArray(s_ann_model_grid_xyz, \
                                        dims=['k', 'tile','j','i'])
s_ann_model_grid_xyz_DA = s_ann_model_grid_xyz_DA.assign_coords({'j' : range(90)})
s_ann_model_grid_xyz_DA = s_ann_model_grid_xyz_DA.assign_coords({'i' : range(90)})
s_ann_model_grid_xyz_DA = s_ann_model_grid_xyz_DA.assign_coords({'k' : range(50)})
s_ann_model_grid_xyz_DA = s_ann_model_grid_xyz_DA.assign_coords({'tile' : range(13)})
s_ann_model_grid_xyz_DA = s_ann_model_grid_xyz_DA.assign_coords({'XC' : (('tile', 'j', 'i'), model_grid.XC)})
s_ann_model_grid_xyz_DA = s_ann_model_grid_xyz_DA.assign_coords({'YC' : (('tile', 'j', 'i'), model_grid.YC)})
s_ann_model_grid_xyz_DA = s_ann_model_grid_xyz_DA.assign_coords({'Z'  : ('k', model_grid.Z)})

s_ann_model_grid_xyz_DA.name = 'SALT'
s_ann_model_grid_xyz_DA.attrs = data_field_S

s_ann_model_grid_xyz_DS = s_ann_model_grid_xyz_DA.to_dataset()
s_ann_model_grid_xyz_DS.attrs = new_data_attr_S

# ANNUAL THETA

pott_ann_model_grid_xyz_DA = xr.DataArray(pott_ann_model_grid_xyz, \
                                        dims=['k', 'tile','j','i'])
pott_ann_model_grid_xyz_DA = pott_ann_model_grid_xyz_DA.assign_coords({'j' : range(90)})
pott_ann_model_grid_xyz_DA = pott_ann_model_grid_xyz_DA.assign_coords({'i' : range(90)})
pott_ann_model_grid_xyz_DA = pott_ann_model_grid_xyz_DA.assign_coords({'k' : range(50)})
pott_ann_model_grid_xyz_DA = pott_ann_model_grid_xyz_DA.assign_coords({'tile' : range(13)})
pott_ann_model_grid_xyz_DA = pott_ann_model_grid_xyz_DA.assign_coords({'XC' : (('tile', 'j', 'i'), model_grid.XC)})
pott_ann_model_grid_xyz_DA = pott_ann_model_grid_xyz_DA.assign_coords({'YC' : (('tile', 'j', 'i'), model_grid.YC)})
pott_ann_model_grid_xyz_DA = pott_ann_model_grid_xyz_DA.assign_coords({'Z'  : ('k', model_grid.Z)})

pott_ann_model_grid_xyz_DA.name = 'THETA'
pott_ann_model_grid_xyz_DA.attrs = data_field_T

pott_ann_model_grid_xyz_DS = pott_ann_model_grid_xyz_DA.to_dataset()
pott_ann_model_grid_xyz_DS.attrs = new_data_attr_T


# ADD MAPPING TIME

s_monthly_plus_ann_model_grid_xyz_DS.attrs['created'] = mapping_time
pott_monthly_plus_ann_model_grid_xyz_DS.attrs['created'] = mapping_time
s_ann_model_grid_xyz_DS.attrs['created'] = mapping_time
pott_ann_model_grid_xyz_DS.attrs['created'] = mapping_time


# SAVE TO DISK
fname= 'MON_CLIM_SALT_WOA18_' + model_grid.name
#s_monthly_plus_ann_model_grid_xyz_DS.to_netcdf(netcdf_output_dir / fname)

save_to_disk(s_monthly_plus_ann_model_grid_xyz_DS, fname,binary_fill_value,\
             netcdf_fill_value, netcdf_output_dir,binary_output_dir,\
             binary_output_dtype, 'llc')


fname= 'MON_CLIM_THETA_WOA18_' + model_grid.name
save_to_disk(pott_monthly_plus_ann_model_grid_xyz_DS, fname,binary_fill_value,\
             netcdf_fill_value, netcdf_output_dir,binary_output_dir,\
             binary_output_dtype, 'llc')

fname= 'ANN_CLIM_SALT_WOA18_' + model_grid.name
#s_ann_model_grid_xyz_DS.to_netcdf(netcdf_output_dir / fname)
save_to_disk(s_ann_model_grid_xyz_DS, fname,binary_fill_value,\
             netcdf_fill_value, netcdf_output_dir,binary_output_dir,\
             binary_output_dtype, 'llc')


fname= 'ANN_CLIM_THETA_WOA18_' + model_grid.name
#pott_ann_model_grid_xyz_DS.to_netcdf(netcdf_output_dir / fname)
save_to_disk(pott_ann_model_grid_xyz_DS, fname,binary_fill_value,\
             netcdf_fill_value, netcdf_output_dir,binary_output_dir,\
             binary_output_dtype, 'llc')


