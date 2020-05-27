#!/usr/bin/env python
# coding: utf-8

import numpy as np
import xarray as xr
import pyresample as pr
import sys
from datetime import datetime
sys.path.append('/home/ifenty/ECCOv4-py')
import ecco_v4_py as ecco
from pathlib import Path
import netCDF4 as nc4

np.warnings.filterwarnings('ignore')

#%%
# ### For Reference: Python binary file i/o specifications
# ```
# endianness:
#     big     : '>'
#     little  : '<'
# 
# precision  
#     float32':'f4',
#     float64':'f8',
# ```


## ROUTINE TO MAKE AN EMPTY DATASET WITH SST AND SST ERROR 

# create new SST and SST error DataArray objects that will hold
# the interpolated data for a specific day (current_day)

# the spatial dimensions are the same as the model grid coordinates (XC, YC)
# and a single time dimension (day) is added
#
# arrays are initially filled with 'fill_value'

#Fill value for f8
#    X.attrs['_FillValue'] = nc4.default_fillvals['f8']




def make_time_bounds_from_ds64(rec_avg_end, AVG_PERIOD):
    """

    Given a datetime64 object (rec_avg_end) representing the 'end' of an 
    1M averaging time period  create a time_bounds array
    with two datetime64 variables, one for the averaging period start, one
    for the averaging period end.  Also find the middle time between the
    two.
    
    Parameters
    ----------    
    rec_avg_end : numpy.datetime64 
        the time at the end of an averaging period


    Returns
    -------
    time_bnds : numpy.array(dtype=numpy.datetime64)
        a datetime64 array with the start and end time of the averaging periods
    
    center_times : numpy.datetime64
        the 'center' of the averaging period
    
    """ 

    rec_year = int(str(rec_avg_end)[:4])
    rec_mon = int(str(rec_avg_end)[5:7])
    rec_day = int(str(rec_avg_end)[8:10])
    
    rec_avg_end_as_dt = datetime(rec_year, rec_mon, 
                                          rec_day)
 
    import dateutil as dateutil
    
    if AVG_PERIOD == 'AVG_MON':
        rec_avg_start =  rec_avg_end_as_dt - \
            dateutil.relativedelta.relativedelta(months=1)    

    elif AVG_PERIOD == 'AVG_DAY':
        rec_avg_start =  rec_avg_end_as_dt - \
           dateutil.relativedelta.relativedelta(days=1)    
        
    rec_avg_start =  np.datetime64(rec_avg_start)
    
    rec_avg_delta = rec_avg_end - rec_avg_start
    rec_avg_middle = rec_avg_start + rec_avg_delta/2
    
    rec_time_bnds = np.array([rec_avg_start, rec_avg_end])
    
    return rec_time_bnds, rec_avg_middle

#%%

def make_empty_record(standard_name, long_name, \
                      units, current_time, model_grid, array_precision):
    
    # make empty data array to hold the interpolated 2D field
    data_DA =xr.DataArray(np.ones(np.shape(model_grid.XC.values),\
                                  dtype=array_precision), \
                          dims=model_grid.XC.dims)*np.nan

    for dim in model_grid.XC.dims:
        data_DA=data_DA.assign_coords({dim : model_grid[dim]})

    data_DA=data_DA.assign_coords(time = np.datetime64(current_time,'D'))
    data_DA=data_DA.expand_dims(dim = 'time', axis=0)
    
    # create XC and YC coordinates
    if 'tile' in model_grid.XC.dims:
        # llc grid has 'tile dimension'
        data_DA=data_DA.assign_coords({'XC': (('tile','j','i'), model_grid.XC)})
        data_DA=data_DA.assign_coords({'YC': (('tile','j','i'), model_grid.YC)})
 
        # some grids only have j i dimensions
    else:
        data_DA=data_DA.assign_coords({'XC': (('j','i'), model_grid.XC)})
        data_DA=data_DA.assign_coords({'YC': (('j','i'), model_grid.YC)})
     
    data_DA.XC.attrs = model_grid.XC.attrs
    data_DA.YC.attrs = model_grid.YC.attrs
 
    data_DA.attrs = []
    data_DA.attrs['interpolated_grid'] = model_grid.title
    data_DA.attrs['long_name'] = long_name
    data_DA.attrs['standard_name'] = standard_name
    data_DA.attrs['units'] = units

    return data_DA
#%%



def transform_to_model_grid(
        data_indices_at_model_index_i,
        num_data_indices_at_model_index_i,
        day, file, \
        sst_field, standard_name, long_name, units,array_precision):
                    


    # This product's data resolution is 0.25 degrees
    # So the search should be <= 30km (110 km / degree * 0.25 degrees = 27 km)

    #interpolation_parameters = 'nearest neighbor, 30km search radius'
    interpolation_parameters = 'Gaussian weighting, 55km search radius, sigma=25km'


    # create empty data array     
    sst_DA =  make_empty_record(standard_name, long_name, \
                                    units, day, model_grid,array_precision)

    
    if isinstance(file, Path) == False:
        print('could not read file!')
    else:
        # read the original data dataset     
        print('reading ', file.name)

    ds = xr.open_dataset(file, decode_times=True)
    
    # extract the original sea ice concentration data
    orig_data = ds[sst_field].values[0,:]
    
    
    if np.sum(~np.isnan(orig_data)) > 0:
        
        
        # convert to Celsius
        if sst_field == 'analysed_sst':
            orig_data -= 273.15
       
    
        # make a 1D version of the original data        
        orig_data_r = orig_data.ravel()
    
        # make an empty array that will hold the data interpolated/averaged to
        # the grid of the model
        data_model_projection = np.zeros(model_grid.XC.shape)*np.nan
        
        # get a 1D version of data_model_projection
        dmp_r = data_model_projection.ravel()
     
        # loop through every model grid point
        for i in range(len(dmp_r)):
            # if the number of *data* points at this model grid point is > 0
            if num_data_indices_at_model_index_i[i] > 0:
                # average those together, put the average in dmp_r
                dmp_r[i] = orig_data_r[data_indices_at_model_index_i[i]].mean()
            
        
        # put the new SST values into the sst_DA array.
        # --where the mapped data are not nan, replace the original values
        # --where they are nan, just leave the original values alone
    
        #print ('tmp values shape ', hemi_i, tmp_values.shape)
        sst_DA.values  = np.where(~np.isnan(data_model_projection), \
                                  data_model_projection, sst_DA.values)
                
    
        ##% MASK SST WHERE THERE IS SEA ICE
        [llo, lla] = np.meshgrid(ds.lon.values, ds.lat.values)
        
        llo_ice = llo[np.where(ds.sea_ice_fraction.values[0,:]> 0)]
        lla_ice = lla[np.where(ds.sea_ice_fraction.values[0,:]> 0)]
    
        sea_ice_data = ds.sea_ice_fraction.values[0,:]
        
        if np.sum(~np.isnan(sea_ice_data)) > 0:
        
            # array of zeros, size of mask
            tmp = ds.mask[0,:].values * 0 + 1
            sea_ice_mask = tmp[np.where(ds.sea_ice_fraction.values[0,:]>0)]
        
            # make a special "swath" for where there is sea ice 
            sea_ice_swath = pr.geometry.SwathDefinition(lons=llo_ice, lats=lla_ice)
                
            # project nonzero sea ice concentration values with 50 km to the model grid
            roi = 50.0e3
            
            # do the actual resampling    
            sea_ice_mask_model_projection = \
                pr.kd_tree.resample_nearest(sea_ice_swath, sea_ice_mask, model_swath,\
                                fill_value=np.nan, \
                                radius_of_influence=roi)
            
            # resape to model grid    
            sea_ice_mask_model_projection = np.reshape(sea_ice_mask_model_projection, \
                                                       sst_DA.values.shape)
                
        #    plt.figure(num=2,clear=True);
        #    ecco.plot_proj_to_latlon_grid(model_grid.XC, model_grid.YC, \
        #           sea_ice_mask_model_projection,
        #           show_colorbar=True,cmin=0,cmax=2)
                
            # mask out places with nonzero sea ice mask
            sst_DA.values = np.where(sea_ice_mask_model_projection == 1, np.nan, \
                                     sst_DA.values)

    
    if 1 == 0:    
        plt.figure(num=3,clear=True);
        ecco.plot_proj_to_latlon_grid(model_grid.XC, model_grid.YC, \
           sst_DA.values[0,:],
           show_colorbar=True)
       
    
    if 'time_bnds' in ds.keys():
        # make time bounds for this record (1 day )
        tb, ct =  make_time_bounds_from_ds64(ds.time_bnds.values[0,-1], 'AVG_DAY') 
    else:
        tb, ct = \
            make_time_bounds_from_ds64(ds.time.values[0] + np.timedelta64(1,'D'), 'AVG_DAY') 
         
    # start time is 30 days earlier
    avg_start_time = sst_DA.time.copy(deep=True)
    avg_start_time.values[0] = tb[0]

    avg_end_time = sst_DA.time.copy(deep=True)
    avg_end_time.values[0] = tb[1]    
    
    avg_center_time = sst_DA.time.copy(deep=True)
    avg_center_time.values[0] = ct    

    # we'll make the center of the averaging time
    sst_DA=sst_DA.assign_coords({'time_start': ('time', avg_start_time)})
    sst_DA=sst_DA.assign_coords({'time_end': ('time', avg_end_time)})
    
    # halfway through the approx 1M averaging period.
    sst_DA.time.values[0] = ct
    sst_DA.time.attrs['long_name'] = 'center time of 1D averaging period'
            
    sst_DA.attrs['original_filename'] = file.name
    sst_DA.attrs['original_field_name'] = sst_field
    sst_DA.attrs['interplation_parameters'] = interpolation_parameters
    sst_DA.attrs['interplation_code'] = 'pyresample'
    sst_DA.attrs['interpolation_date'] = str(np.datetime64(datetime.now(),'D'))

    
    ## Return the new data array     
    return sst_DA



#%%
def get_data_filepaths_for_year(year, data_dir):
      
    # make empty dictionaries,
    daily_filenames = dict()

#   # find all etcdf files in this directory that have the year and nh or sh
    all_netcdf_files_year = np.sort(list(data_dir.glob('*' + str(year) + '*.nc')))
     
    print(all_netcdf_files_year)
    # extract just the filenames from the full pathnames
    all_netcdf_files_basename = [x.name for x in all_netcdf_files_year]
    
    # make an array with all of the days in the year
    dates_in_year = \
        np.arange(str(year) + '-01-01', str(year+1) + '-01-01', dtype='datetime64[D]')
    
    # make empty list that will contain the dates in this year in iso format
    # yyyy-mm-dd
    dates_in_year_iso = []
    
    # loop through every day in the year
    for day in dates_in_year:
        
        # construct date string (yyyymmdd)
        date_str = str(day.tolist().year) + \
            str(day.tolist().month).zfill(2) + \
            str(day.tolist().day).zfill(2)
  
        # construct the date string in iso format
        date_str_iso = str(day.tolist().year) + '-' + \
            str(day.tolist().month).zfill(2) + '-' + \
            str(day.tolist().day).zfill(2)
        
        # add iso format date to dates_in_year_iso
        dates_in_year_iso.append(date_str_iso)
        
        #print(date_str)

        # find the filename that matches this day
        # filenames have yyyymmdd in them, not yyyy-mm-dd
        files_day = []
        
        # careful, only look at the first 8 digits of the filename
        # the string "20010112", jan 12 2001 matched two files
        #
        # 20010112120000-NCEI-L4_GHRSST-SSTblend-AVHRR_OI-GLOB-v02.0-fv02.0.nc
        # ^^^^^^^
        # 20200101120000-NCEI-L4_GHRSST-SSTblend-AVHRR_OI-GLOB-v02.0-fv02.0.nc
        #   ^^^^^^^^
        #
        # because this product includes the Time 12:00 in the filename
        for netcdf_file in all_netcdf_files_basename:
            if str.find(netcdf_file[:8], date_str) >= 0:
                files_day = netcdf_file
        
        # add this filename to the dictionary with the date_str_iso as the
        # key
        daily_filenames[date_str_iso] = files_day
        
    # return the dates in iso format, and the filenames for nh and sh
    return dates_in_year_iso, daily_filenames


#%%    

def save_to_disk(sst_daily_DA_year_merged, sst_mon_DA_year_merged,\
                 output_daily_filename, output_monthly_filename, \
                 binary_fill_value, netcdf_fill_value,
                 netcdf_output_dir, binary_output_dir, 
                 binary_output_dtype, model_grid_type)  :  
    
    # define binary file output filetype    
    dt_out = np.dtype(binary_output_dtype)

    # create directories
    netcdf_output_dir.mkdir(exist_ok=True)
    binary_output_dir.mkdir(exist_ok=True)

    # daily first
    
    for freq in range(2):
        
        if freq == 0:  #daily
            netcdf_output_filename = netcdf_output_dir / Path(output_daily_filename + '.nc')
            binary_output_filename = binary_output_dir / Path(output_daily_filename )

            F = sst_daily_DA_year_merged
        elif freq == 1: # monthly
            netcdf_output_filename = netcdf_output_dir / Path(output_monthly_filename + '.nc')
            binary_output_filename = binary_output_dir / Path(output_monthly_filename )

            F = sst_mon_DA_year_merged
            
        # first save binary version of data
        
        # replace nans with the binary fill value (something like -9999)
        F.values = \
            np.where(np.isnan(F.values), binary_fill_value,\
                     F.values)
    

        # loop through each record of the year, save one at a time
        for i in range(len(F.time)):
            print ('saving record: ', str(i))
            
            # if we have an llc grid, then we have to reform to compact
            if model_grid_type == 'llc':
                tmp_data       = \
                    ecco.llc_tiles_to_compact(F.values[i,:], \
                                              less_output=True)
                
            # otherwise assume grid is x,y (2 dimensions)
            else:
                tmp_data  = F.values[i,:]
                
            # if this is the first record, create new files
            if i == 0:
                fd1 = open(str(binary_output_filename),'wb')
               
            # otherwise append to the existing binary file
            else:
                fd1 = open(str(binary_output_filename),'ab')
    
            # the actual save commands.
            tmp_data.astype(dt_out).tofile(fd1)
            
            # cleanup
            fd1.close()
        
            ### SAVE NETCDF
            # replace the binary fill value (-9999) with the netcdf fill value
            # which is much more interesting
            F.values = \
                np.where(F.values == binary_fill_value, \
                         netcdf_fill_value,\
                         F.values)
                
            # record the new fill value 
            F.attrs['_FillValue'] = netcdf_fill_value
            F.XC.attrs['_FillValue'] = netcdf_fill_value
            F.YC.attrs['_FillValue'] = netcdf_fill_value
        
        
            # the actual saving (so easy with xarray!)    
            F.to_netcdf(netcdf_output_filename)
            F.close()
    
        
#%%
##################################################
if __name__== "__main__":
    
    #######################################################
    ##  BEGIN  RUN-TIME SPECIFIC PARAMETERS              ##

    # model grid file to use for interpolation
    # ------------------------------------------
    #   * model grid must be provided as a netcdf file with XC and YC fields (lon, lat)
    #   * the netcdf file must also have an attribute (metadata) field 'title'
    #   * model grid can have multiple tiles (or facets or faces)

    model_grid_dir = Path('/home/ifenty/ECCOv4-release/Release4/nctiles_grid')
    model_grid_filename  = 'ECCO-GRID.nc'
    model_grid_id   = 'llc90'   
    model_grid_type = 'llc'
    model_grid_search_radius_max = 55000.0 # m
    
    # output parameters 
    # -----------------
    mapping_time =  datetime.now().strftime("%Y%m%dT%H%M%S")
    netcdf_output_dir = Path('/mnt/intraid/ian1/ifenty/data/observations/SST/AVHRR_OI-NCEI-L4_GLOB-v2.0/mapped_to_' + model_grid_id + '/' + mapping_time + '/netcdf')
    binary_output_dir = Path('/mnt/intraid/ian1/ifenty/data/observations/SST/AVHRR_OI-NCEI-L4_GLOB-v2.0/mapped_to_' + model_grid_id + '/' + mapping_time + '/binary')
    
    # Define precision of output files, float32 is standard
    # ------------------------------------------------------
    array_precision = np.float32
    
    
    # Define fill values for binary and netcdf
    # ---------------------------------------------
    if array_precision == np.float32:
        binary_output_dtype = '>f4'
        netcdf_fill_value = nc4.default_fillvals['f4']

    elif array_precision == np.float64:
        binary_output_dtype = '>f8'
        netcdf_fill_value = nc4.default_fillvals['f8']

    # ECCO always uses -9999 for missing data.
    binary_fill_value = -9999
           
    
    # fields to process
    # ---------------------
    #sst_fields = ['analysis_error']
    #sst_fields = ['analysed_sst']
    sst_fields = ['analysed_sst', 'analysis_error']
    

    # years to process
    # --------------------------

    years = np.arange(1992,2021) 

    # Location of original net  (easier if all fields are softlinked to one spot)
    # ---------------------------------------------------------------------------
    data_dir = Path('/home/ifenty/data/observations/SST/AVHRR_OI-NCEI-L4_GLOB-v2.0/all/')

    #%%
    ##  END RUN-TIME SPECIFIC PARAMETERS                 ##
    #######################################################
    
    
    #######################################################
    ##  BEGIN DATA PRODUCT PARAMETERS                    ##

    # these fields were taken from the product.
    


    product_name = 'AVHRR_OI-NCEI-L4-GLOB-v2.0' 
    
    areaExtent = (-180.0, 90.0, 180.0, -90.0)

    #projection = 'init=EPSG:4326'  # Use 'EPSG:3409' with pyproj 2.0+
    epsg, proj, pName = '4326', 'longlat', 'Geographic'
    projDict = {'proj': proj, 'datum': 'WGS84', 'units': 'degree'}
    cols = 1440; rows=720
    #data_grid = geom.AreaDefinition(epsg, pName, proj, projDict, cols, rows, areaExtent)

    projs = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
    data_grid = pr.area_config.get_area_def('longlat', 'Plate Carree', \
                                            'EPSG:4326', projs, 1440,720,\
                                            areaExtent)
    
    
    #%%
    ##  END   DATASET PARAMETERS                         ##
    #######################################################



    #######################################################
    ## BEGIN TRANSFORMATION LOOP                         ## 
   
    # make output directories
    netcdf_output_dir.mkdir(exist_ok=True, parents=True)
    binary_output_dir.mkdir(exist_ok=True, parents=True)
    
    # load the model grid
    model_grid = xr.open_dataset(model_grid_dir / model_grid_filename)
    model_grid  = model_grid.reset_coords()

  
    ##### PRE-COMPUTE THE MAPPING BETWEEN DATA GRID AND MODEL GRID
    
    # Define the 'swath' (in the terminology of the pyresample module)
    # as the lats/lon pairs of the model grid
    # The routine needs the lats and lons to be one-dimensional vectors.
    model_swath  = \
        pr.geometry.SwathDefinition(lons=model_grid.XC.values.ravel(), 
                                    lats=model_grid.YC.values.ravel())
    
    
    # find out the closest 24 SST data grid points to each model grid point
    # for this data product, 24 is plenty for grid 1-degree or finer 
    
    # "model_grid_search_radius_max" is the half of the distance between 
    # the most widely spaced model grid cells.  No need to search for SST
    # data grid points more than halfway to the next model grid cell.
    
    # this get_neighbour_info is quite useful.  Ax[2] is the matrix of 
    # closest data grid points for each model grid point
    # Ax[3] is the actual distance
    # also cool is that Ax[3] is sorted, first column is closest, last column
    # is furtherst.
    Ax = pr.kd_tree.get_neighbour_info(data_grid, model_swath, \
                           radius_of_influence=model_grid_search_radius_max,\
                           neighbours=24)
    
    # 1/2 of the average distance between grid cell centeres in both
    # x and y directions
    model_grid_cell_radius = 0.25*(model_grid.dxC.values.ravel() + \
                                   model_grid.dyC.values.ravel())
    
    MGCR_DA = xr.DataArray(model_grid_cell_radius, \
                           dims=('model_grid_cell_index'))
    DIST_DATA_TO_MODEL_DA = xr.DataArray(Ax[3], \
                                         dims=('model_grid_cell_index','data_neighbor_index'))
    
    # T/F values, whether the closest 24 SST points to the model grid
    # area within the model grid cell's 'radius'
    data_within_search_radius = (DIST_DATA_TO_MODEL_DA < MGCR_DA).values
    
    # these two arrays have all of the mapping information
    # -- the first is a dictionary, which SST data grid points are within
    #    the model_grid_cell radius
    data_indices_at_model_index_i = dict()
    
    # -- the second is just a count of the number of SST data grid points
    #    within the model_grid_cell_radius
    num_data_indices_at_model_index_i = np.zeros((model_grid_cell_radius.shape))
    
    # loop through every model grid cell
    for i in range(len(model_grid_cell_radius)):
        
        # Ax[2][i,:] are the closest 24 SST data grid indices 
        #            for model grid cell i
        # data_within_search_radius[i,:] is the T/F array for which
        #            of the closest 24 SST data grid indices are within
        #            the radius of this model grid cell i
        # -- so we're pulling out just those SST data grid indices
        #    that fall within the model grid cell radius
        data_indices_at_model_index_i[i] = \
            Ax[2][i,:][data_within_search_radius[i,:] == True]
            
        # just count the # 
        num_data_indices_at_model_index_i[i] = len(data_indices_at_model_index_i[i] )
        
        # print progress.  always nice 
        if np.mod(i, 100)==0:
            print(i)


        
    #%%
    # process model fields one at a time.
    for sst_field_i, sst_field in enumerate(sst_fields) :
        # loop through different sst fields that are present in the 
        # netcdf files.   
        # they have different long names
        if sst_field == 'analysed_sst':
            long_name = 'analysed sea surface temperature'
            standard_name = 'sea_surface_temperature'
            units = 'Celsius'

        elif sst_field == 'analysis_error':
            long_name = 'estimated error standard deviation of analysed_sst'
            standard_name = []
            units = 'Celsius'
     
        # loop through different years
        for year in years:
            
            # get dates and filenames for this year
            iso_dates_for_year, files = \
                get_data_filepaths_for_year(year, data_dir)
   
            sst_daily_DA_year  = []
            
            # Process each day of the year 
            for day in iso_dates_for_year:
                print(sst_field, day)
                file = []
                
                if len(files[day]) > 0:
                    file = data_dir / files[day]
                
                # send filename to the transformation routine, return data array
                # that contains the sea ice concentration for this day
                sst_daily_DA = transform_to_model_grid(
                    data_indices_at_model_index_i, \
                    num_data_indices_at_model_index_i,
                    day, file, \
                    sst_field, standard_name, long_name, units, \
                    array_precision)
                
                # append this data array to the list of data arrays for this 
                # year
                sst_daily_DA_year.append(sst_daily_DA)
    
            ## END   TRANSFORMATION LOOP                         ## 
            #######################################################
        
            # merge the data arrays for this year into a mega data array 
            # along the time dimension
            sst_daily_DA_year_merged = xr.concat((sst_daily_DA_year), dim='time')
           
            
            # if everything comes back nans it means there were no files
            # to load for the entire year.  don't bother saving the 
            # netcdf or binary files for this year
            
            if np.sum(~np.isnan(sst_daily_DA_year_merged.values)) == 0:
                print('Empty year not writing to disk', year)

            else:
                #######################################################
                ## BEGIN ADD METADATA AND COORDINATES                ## 
                      
                # update the dataset attributes (subset of metadata fields taken from 
                # original dataset )
                sst_daily_DA_year_merged.attrs['original_dataset_title'] = 'NOAA/NCEI 1/4 Degree Daily Optimum Interpolation Sea Surface Temperature (OISST) Analysis, Version 2 - Final'
                sst_daily_DA_year_merged.attrs['original_dataset_url'] =  'http://doi.org/10.7289/V5SQ8XB5'
                sst_daily_DA_year_merged.attrs['original_dataset_reference'] =  'Reynolds, et al.(2009) What is New in Version 2. Available at http://www.ncdc.noaa.gov/sites/default/files/attachments/Reynolds2009_oisst_daily_v02r00_version2-features.pdf;Daily 1/4 Degree Optimum Interpolation Sea Surface Temperature (OISST) - Climate Algorithm Theoretical Basis Document, NOAA Climate Data Record Program CDRP-ATBD-0303 Rev. 2 (2013). Available at http://www1.ncdc.noaa.gov/pub/data/sds/cdr/CDRs/Sea_Surface_Temperature_Optimum_Interpolation/AlgorithmDescription.pdf.'
                sst_daily_DA_year_merged.attrs['original_dataset_product_id'] = 'NCEI-L4LRblend-GLOB-AVHRR_OI'
                sst_daily_DA_year_merged.attrs['interpolated_grid_id'] = model_grid_id
                sst_daily_DA_year_merged.name = sst_field + '_interpolated_to_' + model_grid_id
       
      
                # list holding the dates for the end of each month this year.
                iso_dates_at_end_of_month = []
                iso_dates_at_start_of_month = []
                
                dt64_dates_at_end_of_month = []
                dt64_dates_at_start_of_month = []
                # pull one record per month
                
                sst_mon_DA_year = []
                for month in range(1,13):
                    # to find the last day of the month, we go up one month, 
                    # and back one day
                    #   if Jan-Nov, then we'll go forward one month to Feb-Dec
    
                    if month < 12:
                        cur_mon_year = np.datetime64(str(year) + '-' + \
                                                     str(month+1).zfill(2) + '-'  + str(1).zfill(2))
                                        # for december we go up one year, and set month to january
                    else:
                        cur_mon_year = np.datetime64(str(year+1) + '-' + \
                                                     str('01')+ '-' + str(1).zfill(2))
                    
                    mon_str = str(year) + '-' + str(month).zfill(2)
                    
                    
                    sst_mon_DA = \
                        sst_daily_DA_year_merged.sel(time= mon_str).mean(axis=0, \
                                                    skipna=False, keep_attrs=True )
                    tb, ct = make_time_bounds_from_ds64(cur_mon_year,'AVG_MON')
                    
                    sst_mon_DA =sst_mon_DA.assign_coords({'time' : ct})
                    sst_mon_DA =sst_mon_DA.expand_dims('time', axis=0)
    
                    avg_start_time = sst_mon_DA.time.copy(deep=True)
                    avg_start_time.values[0] = tb[0]
                
                    avg_end_time = sst_mon_DA.time.copy(deep=True)
                    avg_end_time.values[0] = tb[1]    
                    
                    avg_center_time = sst_mon_DA.time.copy(deep=True)
                    avg_center_time.values[0] = ct    
                
                    # we'll make the center of the averaging time
                    sst_mon_DA=sst_mon_DA.assign_coords({'time_start': ('time', avg_start_time)})
                    sst_mon_DA=sst_mon_DA.assign_coords({'time_end': ('time', avg_end_time)})
                    
                    # halfway through the approx 1M averaging period.
                    sst_mon_DA.time.values[0] = ct
                    sst_mon_DA.time.attrs['long_name'] = 'center time of 1M averaging period'
                    
                    sst_mon_DA_year.append(sst_mon_DA)
                    
                sst_mon_DA_year_merged = xr.concat((sst_mon_DA_year), dim='time')
    
    
        
                # OTHER METADATA SHOULD GO HERE!!!
            
            
                ## END   ADD METADATA AND COORDINATES                ## 
                #######################################################
                 
            
            
                #######################################################
                ## BEGIN SAVE TO DISK                                ## 
            
                output_daily_filename = \
                    product_name + '_' + sst_field + '_' + \
                    model_grid_id + '_DAILY_' + str(year)
            
                output_monthly_filename = \
                    product_name + '_' + sst_field + '_' + \
                    model_grid_id + '_MONTHLY_' + str(year)
                     
                save_to_disk(sst_daily_DA_year_merged, sst_mon_DA_year_merged, \
                             output_daily_filename, output_monthly_filename,
                             binary_fill_value, netcdf_fill_value,\
                             netcdf_output_dir, binary_output_dir, 
                             binary_output_dtype, model_grid_type)
                
                ## END   SAVE TO DISK                                ## 
                #######################################################