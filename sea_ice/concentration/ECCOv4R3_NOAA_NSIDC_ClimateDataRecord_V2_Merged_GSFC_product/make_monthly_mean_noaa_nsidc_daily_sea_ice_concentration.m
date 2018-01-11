

% This script makes monthly mean equivalents of the 
%   NOAA/NSIDC Climate Data Record of Passive Microwave Sea Ice 
% Concentration, Version 2 dataset  
% the dataset contains 4 different sea ice products
% (1) NOAA/NSIDC CDR
% (2) NASA Bootstrap
% (3) NASA team
% (4) merged Bootstram + NASA team

% http://nsidc.org/data/docs/noaa/g02202_ice_conc_cdr/index.html#data_description


% Filename; make_monthly_mean_noaa_nsidc_daily_sea_ice_concentration.m
%  ** former filename : 
% Date Created: 2016-02-02
% Last Modified:  


% notes:
 

% 
% This data set provides a Climate Data Record (CDR) of 
% sea ice concentration from passive microwave data.
% It provides a consistent, daily and monthly time series of sea ice
% concentrations from 09 July 1987 through through most recent processing 
% for both the north and south polar regions on a 25 km x 25 km grid.

% The NOAA/NSIDC CDR is based on the recommendations from the 
% National Research Council (NRC) (2004). It is produced from 
% gridded brightness temperatures from the Defense Meteorological 
% Satellite Program (DMSP) F8, F11, and F13 Special Sensor Microwave
% Imager (SSM/I) passive microwave radiometers and the DMSP F17 Special 
% Sensor Microwave Imager/Sounder (SSMIS) passive microwave radiometer.
% 
% The NOAA/NSIDC CDR sea ice concentrations are an estimate of the 
% fraction of ocean area covered by sea ice that is produced by 
% combining concentration estimates created using two algorithms 
% developed at the NASA Goddard Space Flight Center (GSFC): the 
% NASA Team algorithm (Cavalieri et al., 1984) and the Bootstrap
% algorithm (Comiso, 1986). NSIDC applies the individual algorithms
% to brightness temperature data from Remote Sensing Systems, Inc. (RSS).

% The data are gridded on the NSIDC polar stereographic grid with 25 x 25 km
% grid cells and are available in netCDF file format. Each file includes four
% different sea ice concentration variables: a variable with the primary CDR 
% sea ice concentrations created by NSIDC and three variables with sea ice 
% concentrations created by Goddard.
% 
% The three Goddard-produced sea ice concentrations are the Goddard NASA 
% Team algorithm sea ice concentrations, the Goddard Bootstrap sea ice
% concentrations, and a merged version of the two sea ice concentrations. 
% These Goddard-produced sea ice concentrations are included in the data 
% files for a number of reasons. First, the merged Goddard NASA 
% Team/Bootstrap sea ice concentrations are an ancillary data set that is
% analogous to the NOAA/NSIDC CDR data but that adds late 1978 through mid 
% 1987 data to the record. A different instrument, the Scanning Multichannel
% Microwave Radiometer (SMMR), was the source for the brightness temperatures 
% from this period. Sea ice concentrations from the extended period are not 
% part of the primary NSIDC-produced CDR record because complete 
% documentation of the SMMR brightness temperature processing method is 
% not available. Second, the separate Goddard NASA Team and Bootstrap sea 
% ice concentrations are provided for reference.

% flags special for cdr data

% 1 BT_source_for_CDR 
% 2 NT_source_for_CDR 
% 4 no_ice_allowed_per_climatology 
% 8 grid_cell_near_to_coast 
% 32 concentration_below_fifty_percent 
% -128 melt_start_detected

%%

clear all;
%% 

ice_dir = '/ian1/ifenty/data/observations/seaice/NOAA_NSIDC_CLIM_DATA_REC_OF_PM_SEA_ICE/'
cd(ice_dir)

nh_ref = nc_getall([ice_dir '/example_nh_record.nc']);
nh_lats = nh_ref.latitude.data;
nh_lons = nh_ref.longitude.data;

sh_ref = nc_getall([ice_dir '/example_sh_record.nc']);
sh_lats = sh_ref.latitude.data;
sh_lons = sh_ref.longitude.data;

%%
% seaice_conc flag meanings
% -5 pole_hole 
% -4 lakes 
% -3 coastal
% -2 land_mask
% -1 missing_data
% 

tmp_nh_conc = floor(nh_ref.goddard_merged_seaice_conc.data.*100);
tmp_sh_conc = floor(sh_ref.goddard_merged_seaice_conc.data.*100);

nh_pole_ins = find(tmp_nh_conc == -5);
nh_lake_ins = find(tmp_nh_conc == -4);
nh_coastal_ins = find(tmp_nh_conc == -3);
nh_land_ins = find(tmp_nh_conc == -2);

sh_pole_ins = find(tmp_sh_conc == -5);
sh_lake_ins = find(tmp_sh_conc == -4);
sh_coastal_ins = find(tmp_sh_conc == -3);
sh_land_ins = find(tmp_sh_conc == -2);

tmp_nh = zeros(size(tmp_nh_conc));
tmp_sh = zeros(size(tmp_sh_conc));

tmp_nh(nh_lake_ins) = 1;
tmp_nh(nh_land_ins) = 1;
tmp_nh(nh_coastal_ins) = 1;
tmp_nh(nh_pole_ins) = 1;

tmp_sh(sh_lake_ins) = 1;
tmp_sh(sh_land_ins) = 1;
tmp_sh(sh_coastal_ins) = 1;
tmp_sh(sh_pole_ins) = 1;

tmp_nh(1,:) = 0;
tmp_nh(end,:) = 0;
tmp_nh(:,1) = 0;
tmp_nh(:,end) = 0;

tmp_sh(1,:) = 0;
tmp_sh(end,:) = 0;
tmp_sh(:,1) = 0;
tmp_sh(:,end) = 0;

% set the border to be data

nh_no_data_ins = find(tmp_nh == 1);
nh_data_ins    = find(tmp_nh == 0);

sh_no_data_ins = find(tmp_sh == 1);
sh_data_ins    = find(tmp_sh == 0);

figure(1);clf;
subplot(121);
imagesc(tmp_nh);
title('NH no data mask')
subplot(122);
imagesc(tmp_sh);
title('SH no data mask')

 

%% test out the mapping

% seaice_conc flag meanings
% -5 pole_hole 
% -4 lakes 
% -3 coastal
% -2 land_mask
% -1 missing_data
% 
% do the mapping!
 
%%
makefigs=1

for year = 2012:2015
        close all;
     
        NOAA_NSIDC_MONTHLY_MEAN.year = year;
             
        NOAA_NSIDC_MONTHLY_MEAN.flag_meanings =  ...
            nh_ref.seaice_conc_cdr.flag_meanings;
        
        NOAA_NSIDC_MONTHLY_MEAN.flag_values =  ...
            nh_ref.seaice_conc_cdr.flag_values;
        
        NOAA_NSIDC_MONTHLY_MEAN.processing_date = date;
        
        %cur_dir = [ice_dir '/' num2str(year)]
        cur_dir = [ice_dir '/monthly/north_south'];
        cd(cur_dir);
        
        for mon = 1:12
            
            datestr = [padzero(year,4) padzero(mon,2)]
            
            
            ice_files_n = dir(['*' datestr '*nh*nc']);
            ice_files_s = dir(['*' datestr '*sh*nc']);
        
            %%
            
            for f = 1:length(ice_files_n)
            [year f]
            %%
            clear nh_data; clear sh_data;
             cd(cur_dir);
            
            % YYYYMMDD
            nh_fname = ice_files_n(f).name;
            nh_data_date = ice_files_n(f).name(26:33);
            nh_data_year = nh_data_date(1:4);
            nh_data_mon  = nh_data_date(5:6);
            nh_data_day  = nh_data_date(7:8);
            
            sh_fname = ice_files_s(f).name;
            sh_data_date = ice_files_s(f).name(26:33)
            sh_data_year = sh_data_date(1:4);
            sh_data_mon = sh_data_date(5:6);
            sh_data_day = sh_data_date(7:8);
            
            [nh_data_date ' '  nh_data_year ' ' nh_data_mon  ' ' nh_data_day]
            [sh_data_date ' '  sh_data_year ' ' sh_data_mon  ' ' sh_data_day]
            
               
            NOAA_NSIDC_MONTHLY_MEAN.nh_fname{f} = nh_fname;
            NOAA_NSIDC_MONTHLY_MEAN.sh_fname{f} = sh_fname;
            
            %%            
            nh_data = nc_getall(nh_fname);
            sh_data = nc_getall(sh_fname);            
                      
            %%
            nh_g_merged = ...
                (floor(nh_data.goddard_merged_seaice_conc.data.*100));
            
            nh_g_nt = ...
                 (floor(nh_data.goddard_nt_seaice_conc.data.*100));
                       
            nh_g_bt = ...
                 (floor(nh_data.goddard_bt_seaice_conc.data.*100));
                       
            nh_cdr = ...
                 (floor(nh_data.seaice_conc_cdr.data.*100));
               
             %
            sh_g_merged = ...
                 (floor(sh_data.goddard_merged_seaice_conc.data.*100));
            
            sh_g_bt = ...
                 (floor(sh_data.goddard_bt_seaice_conc.data.*100));
            
            sh_g_nt = ....
                 (floor(sh_data.goddard_nt_seaice_conc.data.*100));
            
            sh_cdr= ....
                 (floor(sh_data.seaice_conc_cdr.data.*100));
            
            %%
            nh_goddard_merged(:,:,f) = nh_g_merged;
            nh_goddard_nt(:,:,f)     = nh_g_nt;
            nh_goddard_bt(:,:,f)     = nh_g_bt;
            nh_goddard_cdr(:,:,f)    = nh_cdr;
            
            sh_goddard_merged(:,:,f) = sh_g_merged;
            sh_goddard_nt(:,:,f)     = sh_g_nt;
            sh_goddard_bt(:,:,f)     = sh_g_bt;
            sh_goddard_cdr(:,:,f)    = sh_cdr;
            end % all files in the month
            
            % do a mean here
            %r=input('here')
            
        end % months
        %%
       % cd(ice_dir)
       % fname = ['NOAA_NSIDC_MONTHLY_MEAN_' num2str(year)]
       % save(fname, 'NOAA_NSIDC_MONTHLY_MEAN','-v7.3');
       
end % year