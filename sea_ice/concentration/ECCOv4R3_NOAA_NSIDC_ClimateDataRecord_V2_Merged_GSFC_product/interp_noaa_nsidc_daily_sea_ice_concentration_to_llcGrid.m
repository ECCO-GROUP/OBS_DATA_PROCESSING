% This script interpolates sea ice concentration observations provided
% in the NOAA/NSIDC Climate Data Record of Passive Microwave Sea Ice
% Concentration, Version 2 dataset to an llc grid.

% the dataset contains 4 different sea ice products
% (1) NOAA/NSIDC CDR
% (2) NASA Bootstrap
% (3) NASA team
% (4) merged Bootstram + NASA team

% http://nsidc.org/data/docs/noaa/g02202_ice_conc_cdr/index.html#data_description

% Filename; interp_noaa_nsidc_daily_sea_ice_concentration_to_llcGrid.m
%  ** former filename :
% Date Created: 2014-09-09
% Last Modified: 2017-01-10

% notes:

% 2017-01-10 , updated to end of 2015;
% 2016-03-16, updated to extend dataset to end of 2014.


% 2014-10-06:
%   updated input directory

% -- 2014-09-15, decided to map all values of the original data to the llc
% grid, including points that are 'coast' and 'land'.  The problem came up
% on the ross sea where our llc90 mask has a smaller ice shelf than the
% original data. projecting the nearest neighbor values into that ice shelf
% area is dangerous because coastal values can be quite small. a bunch of
% small coastal values projected over a larger ice shelf area will cause
% adjustments which would create too much open water over a large area.
% best to just keep those points as 'no data'

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

run_codes=[0 1]
ice_dir = '~/data/observations/seaice/NOAA_NSIDC_CLIM_DATA_REC_OF_PM_SEA_ICE/'

for run_code=run_codes
    
    switch run_code
        case 0
            llcN=90;
            years=2015;
            makefigs=1
            input_dir = '/ian1/ifenty/data/observations/seaice/NOAA_NSIDC_CLIM_DATA_REC_OF_PM_SEA_ICE/daily/'
            output_dir = '/ian1/ifenty/data/observations/seaice/NOAA_NSIDC_CLIM_DATA_REC_OF_PM_SEA_ICE/daily_projected_to_matlab/llc090/'
            fig_dir =  '/ian1/ifenty/data/observations/seaice/NOAA_NSIDC_CLIM_DATA_REC_OF_PM_SEA_ICE/daily_figures/llc090/';
            
        case 1
            llcN=270;
            years=2015;
            makefigs=1
            input_dir = '/ian1/ifenty/data/observations/seaice/NOAA_NSIDC_CLIM_DATA_REC_OF_PM_SEA_ICE/daily/'
            output_dir = '/ian1/ifenty/data/observations/seaice/NOAA_NSIDC_CLIM_DATA_REC_OF_PM_SEA_ICE/daily_projected_to_matlab/llc270/'
            fig_dir =  '/ian1/ifenty/data/observations/seaice/NOAA_NSIDC_CLIM_DATA_REC_OF_PM_SEA_ICE/daily_figures/llc270/';
        otherwise
            ['need to specify a run code']
    end
    
    
    %%
    
    switch llcN
        case 90
            load_llc90_grid
            X_llcN = X_90;
            Y_llcN = Y_90;
            Z_llcN = Z_90;
            hf0_N = hf0_90;
            dry_ins_N_k = dry_ins_90_k;
            wet_ins_N_k = wet_ins_90_k;
            lat_llcN = lat_90;
            lon_llcN = lon_90;
            
        case 270
            load_llc270_grid;
            X_llcN = X_270;
            Y_llcN = Y_270;
            Z_llcN = Z_270;
            hf0_N = hf0_270;
            dry_ins_N_k = dry_ins_270_k;
            wet_ins_N_k = wet_ins_270_k;
            
            lat_llcN = lat_270;
            lon_llcN = lon_270;
    end
    
    
    %%
    
    
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
    set(0,'DefaultTextInterpreter','none');
    subplot(121);
    imagesc(tmp_nh);
    title('NH no data mask')
    subplot(122);
    imagesc(tmp_sh);
    title('SH no data mask')
    
    
    
    %%
    % nh
    [X_nh_icec, Y_nh_icec, Z_nh_icec] = ...
        sph2cart(nh_lons(:)*deg2rad, nh_lats(:)*deg2rad, 1);
    xyz_nh =[X_nh_icec, Y_nh_icec, Z_nh_icec];
    
    [X_nh_icec_wet, Y_nh_icec_wet, Z_nh_icec_wet] = ...
        sph2cart(nh_lons(nh_data_ins)*deg2rad, nh_lats(nh_data_ins)*deg2rad, 1);
    
    xy_nh_wet =[X_nh_icec_wet, Y_nh_icec_wet];
    xyz_nh_wet =[X_nh_icec_wet, Y_nh_icec_wet, Z_nh_icec_wet];
    
    % sh
    [X_sh_icec, Y_sh_icec, Z_sh_icec] = ...
        sph2cart(sh_lons(:)*deg2rad, sh_lats(:)*deg2rad, 1);
    xyz_sh =[X_sh_icec, Y_sh_icec, Z_sh_icec];
    
    [X_sh_icec_wet, Y_sh_icec_wet, Z_sh_icec_wet] = ...
        sph2cart(sh_lons(sh_data_ins)*deg2rad, sh_lats(sh_data_ins)*deg2rad, 1);
    
    xy_sh_wet =[X_sh_icec_wet, Y_sh_icec_wet];
    xyz_sh_wet =[X_sh_icec_wet, Y_sh_icec_wet, Z_sh_icec_wet];
    
    
    %%
    AI_nh = 1:length(tmp_nh_conc(:));
    AI_nh_wet = AI_nh(nh_data_ins);
    
    AI_sh = 1:length(tmp_sh_conc(:));
    AI_sh_wet = AI_sh(sh_data_ins);
    
    %%
    F_nh = TriScatteredInterp(xyz_nh, AI_nh','nearest');
    ['Finished F_nh nearest']
    
    F_nh_wet_nearest = TriScatteredInterp(xyz_nh_wet, AI_nh_wet','nearest');
    ['Finished F_nh wet nearest']
    %%
    F_sh = TriScatteredInterp(xyz_sh, AI_sh','nearest');
    ['Finished F_sh nearest']
    
    F_sh_wet_nearest = TriScatteredInterp(xyz_sh_wet, AI_sh_wet','nearest');
    ['Finished F_sh wet nearest']
    
    %%
    figure(10);clf;set(0,'DefaultTextInterpreter','none');
    hold on;
    subplot(221);hold on;
    plot3(xyz_nh(:,1), xyz_nh(:,2), xyz_nh(:,3),'k.')
    plot3(xyz_nh_wet(:,1), xyz_nh_wet(:,2), xyz_nh_wet(:,3),'r.')
    
    subplot(222);hold on;
    plot3(F_nh.X(:,1), F_nh.X(:,2), F_nh.X(:,3),'k.')
    
    subplot(223);hold on;
    plot3(xyz_sh(:,1), xyz_sh(:,2), xyz_sh(:,3),'k.')
    plot3(xyz_sh_wet(:,1), xyz_sh_wet(:,2), xyz_sh_wet(:,3),'r.')
    
    subplot(224);hold on;
    plot3(F_sh.X(:,1), F_sh.X(:,2), F_sh.X(:,3),'k.')
    
    %%
    llcN_NH_ins = find(lat_llcN >= min(nh_lats(:)));
    
    X_llcN_NH = X_llcN(llcN_NH_ins);
    Y_llcN_NH = Y_llcN(llcN_NH_ins);
    Z_llcN_NH = Z_llcN(llcN_NH_ins);
    
    llcN_SH_ins = find(lat_llcN <= max(sh_lats(:)));
    
    X_llcN_SH = X_llcN(llcN_SH_ins);
    Y_llcN_SH = Y_llcN(llcN_SH_ins);
    Z_llcN_SH = Z_llcN(llcN_SH_ins);
    
    %%
    tmp_i = F_nh_wet_nearest(X_llcN_NH, Y_llcN_NH, Z_llcN_NH);
    tmp = hf0_N.*0;
    tmp(llcN_NH_ins) = tmp_i;
    tmp_nh = tmp;
    
    tmp_i = F_sh_wet_nearest(X_llcN_SH, Y_llcN_SH, Z_llcN_SH);
    tmp = hf0_N.*0;
    tmp(llcN_SH_ins) = tmp_i;
    tmp_sh = tmp;
    
    figure(11);clf;set(0,'DefaultTextInterpreter','none');
    subplot(221);
    quikplot_llc(tmp_nh.*hf0_N)
    subplot(222);
    quikplot_llc(tmp_nh);
    subplot(223);
    quikplot_llc(tmp_sh.*hf0_N)
    subplot(224);
    quikplot_llc(tmp_sh)
    
    
    %%
    % make the nearest-neighbor mapping!
    nh_to_llc_mapping = F_nh(X_llcN_NH, Y_llcN_NH, Z_llcN_NH);
    nh_to_llc_mapping (find(isnan(nh_to_llc_mapping )))=1;
    
    nh_wet_nearest_to_llc_mapping = F_nh_wet_nearest(X_llcN_NH, Y_llcN_NH, Z_llcN_NH);
    nh_wet_nearest_to_llc_mapping (find(isnan(nh_wet_nearest_to_llc_mapping )))=1;
    
    %
    sh_to_llc_mapping = F_sh(X_llcN_SH, Y_llcN_SH, Z_llcN_SH);
    sh_to_llc_mapping (find(isnan(sh_to_llc_mapping )))=1;
    
    sh_wet_nearest_to_llc_mapping = F_sh_wet_nearest(X_llcN_SH, Y_llcN_SH, Z_llcN_SH);
    sh_wet_nearest_to_llc_mapping (find(isnan(sh_wet_nearest_to_llc_mapping )))=1;
    ['Finished mapping']
    %%
    
    %% test out the mapping
    
    % seaice_conc flag meanings
    % -5 pole_hole
    % -4 lakes
    % -3 coastal
    % -2 land_mask
    % -1 missing_data
    %
    % do the mapping!
    %nh
    tmp_nh_conc(1,:) = 0;
    tmp_nh_conc(end,:) = 0;
    tmp_nh_conc(:,1) = 0;
    tmp_nh_conc(:,end) = 0;
    
    tmp1 = tmp_nh_conc(nh_to_llc_mapping);
    tmp2 = tmp_nh_conc(nh_wet_nearest_to_llc_mapping);
    
    g1 = hf0_N.*0;
    g1(llcN_NH_ins)=tmp1;
    g2 = hf0_N.*0;
    g2(llcN_NH_ins)=tmp2;
    
    % sh
    tmp_sh_conc(1,:) = 0;
    tmp_sh_conc(end,:) = 0;
    tmp_sh_conc(:,1) = 0;
    tmp_sh_conc(:,end) = 0;
    
    tmp3 = tmp_sh_conc(sh_to_llc_mapping);
    tmp4 = tmp_sh_conc(sh_wet_nearest_to_llc_mapping);
    
    g3 = hf0_N.*0;
    g3(llcN_SH_ins)=tmp3;
    g4 = hf0_N.*0;
    g4(llcN_SH_ins)=tmp4;
    
    
    % plot the results!
    figure(3);clf;set(0,'DefaultTextInterpreter','none');
    subplot(231);
    quikplot_llc((g1));caxis([-5 5]);
    subplot(232);
    quikplot_llc((g2.*hf0_N));caxis([-5 5]);
    subplot(233);
    %g1(find(g1<0))=0;
    quikplot_llc(hf0_N.*(g2-g1));caxis([-5 5]);
    subplot(234);
    quikplot_llc((g3));caxis([-5 5]);
    subplot(235);
    quikplot_llc((g4.*hf0_N));caxis([-5 5]);
    subplot(236);
    %g3(find(g3<0))=0;
    quikplot_llc(hf0_N.*(g4-g3));caxis([-5 5]);
    
    %%
    g6=hf0_N.*0;
    
    g6(llcN_NH_ins)=tmp1;
    g6(llcN_SH_ins)=tmp3;
    
    g7 = hf0_N.*0;
    g7(llcN_NH_ins)=tmp2;
    g7(llcN_SH_ins)=tmp4;
    
    figure(4);clf;set(0,'DefaultTextInterpreter','none');
    subplot(131);
    quikplot_llc(g6);
    title('g6')
    caxis([-5 5]);
    subplot(132);
    quikplot_llc(g7.*hf0_N);
    title('g7.*mask')
    caxis([-5 5]);
    subplot(133);
    quikplot_llc((g7-g6).*hf0_N);
    caxis([-5 5]);colorbar
    ['x']
    
    %% the combined mask of the hFacC + where the sea ice data has a flag <= -2
    %  --- which means lake, coast, land, or pole hole.
    figure(5);clf;set(0,'DefaultTextInterpreter','none');
    g8 = hf0_N.*0;
    g8(find(g6 <= -2)) = 1;
    g9 = hf0_N.*0;
    g9(dry_ins_N_k{1}) = 2;
    
    subplot(131);
    quikplot_llc(g8)
    subplot(132);
    quikplot_llc(g9);
    subplot(133);
    quikplot_llc(g8+g9)
    title('1=sea ice mask only, 2=hf0 mask only, 3=both ice + hf0');colorbar;
    colormap(jet(4));
    caxis([-.5 3.5])
    
    %%
    
    %%
    
    [X_sh_icec, Y_sh_icec, Z_sh_icec] = ...
        sph2cart(sh_lons*deg2rad, sh_lats*deg2rad, 1);
    
    [X_nh_icec, Y_nh_icec, Z_nh_icec] = ...
        sph2cart(nh_lons*deg2rad, nh_lats*deg2rad, 1);
    
    switch llcN
        case 90
            cd(llc90_grid_dir)
            load F_llc90_ALL_INS_SURF_XYZ_to_INDEX;
            F_map = F_llc90_ALL_INS_SURF_XYZ_to_INDEX;
            clear F_llc90_ALL_INS_SURF_XYZ_to_INDEX;
        case 270
            cd(llc270_grid_dir)
            load F_llc270_ALL_INS_SURF_XYZ_to_INDEX;
            F_map = F_llc270_ALL_INS_SURF_XYZ_to_INDEX;
            clear F_llc270_ALL_INS_SURF_XYZ_to_INDEX;
        otherwise
            ['ABORT!!!']
            break;
            return;
    end
    
    llc_all_to_sh_psg_mapping = F_map(X_sh_icec, Y_sh_icec, Z_sh_icec);
    llc_all_to_nh_psg_mapping = F_map(X_nh_icec, Y_nh_icec, Z_nh_icec);
    
    ['finished making llc to sh nh polar stereographic mapping']
    
    %%
    
    for year = years
        
        close all;
        clear NOAA_NSIDC_DAILY_MAPPED_TO_LLC;
        
        NOAA_NSIDC_DAILY_MAPPED_TO_LLC.year = year;
        NOAA_NSIDC_DAILY_MAPPED_TO_LLC.llcN = llcN;
        
        NOAA_NSIDC_DAILY_MAPPED_TO_LLC.flag_meanings =  ...
            nh_ref.seaice_conc_cdr.flag_meanings;
        
        NOAA_NSIDC_DAILY_MAPPED_TO_LLC.flag_values =  ...
            nh_ref.seaice_conc_cdr.flag_values;
        
        NOAA_NSIDC_DAILY_MAPPED_TO_LLC.processing_date = date;
        
        cd(input_dir)
        cur_dir = [input_dir '/' num2str(year)]
        cd(cur_dir);
        
        ice_files_n = dir('*nh*nc');
        ice_files_s = dir('*sh*nc');
        
        %%
        for f = 1:length(ice_files_n)
            [year f]
            %%
            clear nh_data; clear sh_data;
            clear *mapped;
            cd(cur_dir);
            
            % YYYYMMDD
            nh_fname = ice_files_n(f).name;
            nh_data_date = ice_files_n(f).name(26:33);
            nh_data_year = nh_data_date(1:4);
            nh_data_mon = nh_data_date(5:6);
            nh_data_day = nh_data_date(7:8);
            
            sh_fname = ice_files_s(f).name;
            sh_data_date = ice_files_s(f).name(26:33)
            sh_data_year = sh_data_date(1:4);
            sh_data_mon = sh_data_date(5:6);
            sh_data_day = sh_data_date(7:8);
            
            [nh_data_date ' '  nh_data_year ' ' nh_data_mon  ' ' nh_data_day]
            [sh_data_date ' '  sh_data_year ' ' sh_data_mon  ' ' sh_data_day]
            
            NOAA_NSIDC_DAILY_MAPPED_TO_LLC.yyyy_mm_dd(f,:) = ...
                [str2num(nh_data_year) str2num(nh_data_mon) ...
                str2num(nh_data_day)]
            NOAA_NSIDC_DAILY_MAPPED_TO_LLC.nh_fname{f} = nh_fname;
            NOAA_NSIDC_DAILY_MAPPED_TO_LLC.sh_fname{f} = sh_fname;
            
            %%
            nh_data = nc_getall(nh_fname);
            sh_data = nc_getall(sh_fname);
            
            %%
            nh_g_merged = ...
                zero_frame(floor(nh_data.goddard_merged_seaice_conc.data.*100));
            
            nh_g_nt = ...
                zero_frame(floor(nh_data.goddard_nt_seaice_conc.data.*100));
            
            nh_g_bt = ...
                zero_frame(floor(nh_data.goddard_bt_seaice_conc.data.*100));
            
            nh_cdr = ...
                zero_frame(floor(nh_data.seaice_conc_cdr.data.*100));
            
            sh_g_merged = ...
                zero_frame(floor(sh_data.goddard_merged_seaice_conc.data.*100));
            
            sh_g_bt = ...
                zero_frame(floor(sh_data.goddard_bt_seaice_conc.data.*100));
            
            sh_g_nt = ....
                zero_frame(floor(sh_data.goddard_nt_seaice_conc.data.*100));
            
            sh_cdr= ....
                zero_frame(floor(sh_data.seaice_conc_cdr.data.*100));
            
            
            %% do the nearest-neighbor mapping
            nh_g_merged_mapped = nh_g_merged(nh_to_llc_mapping);
            nh_g_nt_mapped     = nh_g_nt(nh_to_llc_mapping);
            nh_g_bt_mapped     = nh_g_bt(nh_to_llc_mapping);
            nh_cdr_mapped      = nh_cdr(nh_to_llc_mapping);
            
            sh_g_merged_mapped = sh_g_merged(sh_to_llc_mapping);
            sh_g_nt_mapped     = sh_g_nt(sh_to_llc_mapping);
            sh_g_bt_mapped     = sh_g_bt(sh_to_llc_mapping);
            sh_cdr_mapped      = sh_cdr(sh_to_llc_mapping);
            
            %%
            global_g_merged_mapped = hf0_N.*0;
            global_g_merged_mapped(llcN_NH_ins)=nh_g_merged_mapped;
            global_g_merged_mapped(llcN_SH_ins)=sh_g_merged_mapped;
            
            global_nt_mapped = hf0_N.*0;
            global_nt_mapped(llcN_NH_ins)=nh_g_nt_mapped;
            global_nt_mapped(llcN_SH_ins)=sh_g_nt_mapped;
            
            global_bt_mapped = hf0_N.*0;
            global_bt_mapped(llcN_NH_ins)=nh_g_bt_mapped;
            global_bt_mapped(llcN_SH_ins)=sh_g_bt_mapped;
            
            global_cdr_mapped = hf0_N.*0;
            global_cdr_mapped(llcN_NH_ins)=nh_cdr_mapped;
            global_cdr_mapped(llcN_SH_ins)=sh_cdr_mapped;
            
            %% wet nearest mapping
            
            %% do the nearest-neighbor mapping
            nh_g_merged_mapped = nh_g_merged(nh_wet_nearest_to_llc_mapping);
            nh_g_nt_mapped     = nh_g_nt(nh_wet_nearest_to_llc_mapping);
            nh_g_bt_mapped     = nh_g_bt(nh_wet_nearest_to_llc_mapping);
            nh_cdr_mapped      = nh_cdr(nh_wet_nearest_to_llc_mapping);
            
            sh_g_merged_mapped = sh_g_merged(sh_wet_nearest_to_llc_mapping);
            sh_g_nt_mapped     = sh_g_nt(sh_wet_nearest_to_llc_mapping);
            sh_g_bt_mapped     = sh_g_bt(sh_wet_nearest_to_llc_mapping);
            sh_cdr_mapped      = sh_cdr(sh_wet_nearest_to_llc_mapping);
            
            %%
            global_g_merged_wet_mapped = hf0_N.*0;
            global_g_merged_wet_mapped(llcN_NH_ins)=nh_g_merged_mapped;
            global_g_merged_wet_mapped(llcN_SH_ins)=sh_g_merged_mapped;
            
            global_nt_wet_mapped = hf0_N.*0;
            global_nt_wet_mapped(llcN_NH_ins)=nh_g_nt_mapped;
            global_nt_wet_mapped(llcN_SH_ins)=sh_g_nt_mapped;
            
            global_bt_wet_mapped = hf0_N.*0;
            global_bt_wet_mapped(llcN_NH_ins)=nh_g_bt_mapped;
            global_bt_wet_mapped(llcN_SH_ins)=sh_g_bt_mapped;
            
            global_cdr_wet_mapped = hf0_N.*0;
            global_cdr_wet_mapped(llcN_NH_ins)=nh_cdr_mapped;
            global_cdr_wet_mapped(llcN_SH_ins)=sh_cdr_mapped;
            %%
            
            if makefigs
                
                dmax =  100;
                dmin = 0
                m=20;
                land_ins = find(hf0_N ==0);
                
                close all;
                figure(10);clf;set(0,'DefaultTextInterpreter','none');
                subplot(221);
                data = global_g_merged_mapped;
                nan_ins = find(data < 0);
                [pdata, pcax, pc] = ...
                    prep_data_land_nan_minmax(data, land_ins, nan_ins, ...
                    dmax, dmin, m,[.5 .5 .5],[1 1 1]);
                quikplot_llc(pdata);  colormap(pc);    caxis(pcax);
                axis tight; axis off;   thincolorbar;
                title('goddard nt + bt merged')
                
                subplot(222);
                data = global_nt_mapped;
                nan_ins = find(data < 0);
                [pdata, pcax, pc] = ...
                    prep_data_land_nan_minmax(data, land_ins, nan_ins, ...
                    dmax, dmin, m,[.5 .5 .5],[1 1 1]);
                quikplot_llc(pdata);  colormap(pc);    caxis(pcax);
                axis tight; axis off;   thincolorbar;
                title('goddard nt ')
                
                subplot(223);
                data = global_bt_mapped;
                nan_ins = find(data < 0);
                [pdata, pcax, pc] = ...
                    prep_data_land_nan_minmax(data, land_ins, nan_ins, ...
                    dmax, dmin, m,[.5 .5 .5],[1 1 1]);
                quikplot_llc(pdata);  colormap(pc);    caxis(pcax);
                axis tight; axis off;   thincolorbar;
                title('goddard bt ')
                
                subplot(224);
                data = global_cdr_mapped;
                nan_ins = find(data < 0);
                [pdata, pcax, pc] = ...
                    prep_data_land_nan_minmax(data, land_ins, nan_ins, ...
                    dmax, dmin, m,[.5 .5 .5],[1 1 1]);
                quikplot_llc(pdata);  colormap(pc);    caxis(pcax);
                axis tight; axis off;   thincolorbar;
                title('noaa/nsidc cdr')
                
                A=['sea ice concentration mapped to llc' padzero(llcN,2)];
                B=ice_files_n(f).name;
                suptitle({A,B})
                tit = ['noaa_nsidc_to_llc' num2str(llcN) '_global_diffs_' B]
                paperx=12;papery=8;prep_figure_for_print;
                
                cd(fig_dir);
                print(gcf,'-dpng','-r100',[tit '.png']);
                
                %%
                
                close all;
                figure(12);clf;
                set(0,'DefaultTextInterpreter','none');
                
                land_psg = hf0_N(llc_all_to_nh_psg_mapping);
                land_ins = find(land_psg ==0);
                
                for i = 1:4
                    switch i
                        case 1
                            tmp = global_g_merged_mapped;
                            orig_nh = nh_g_merged;
                            tit = 'goddard merged'
                        case 2
                            tmp = global_nt_mapped;
                            orig_nh = nh_g_nt;
                            tit = 'goddard nt'
                        case 3
                            tmp = global_bt_mapped;
                            orig_nh = nh_g_bt;
                            tit = 'goddard bt'
                        case 4
                            tmp = global_cdr_mapped;
                            orig_nh = nh_cdr;
                            tit = 'cdr'
                    end
                    dmin=0;
                    dmax=100;
                    m=100;
                    
                    subplot(3,4,i)
                    tmp(dry_ins_N_k{1})=-1;
                    
                    data = tmp(llc_all_to_nh_psg_mapping);
                    nan_ins = find(isnan(data));
                    
                    [pdata, pcax, pc] = ...
                        prep_data_land_nan_minmax(data, land_ins, nan_ins, ...
                        dmax, dmin, m,[.5 .5 .5],[1 1 1]);
                    
                    imagesc(pdata);  colormap(pc);    caxis(pcax);
                    axis tight; axis off;   thincolorbar;
                    title(tit)
                    
                    subplot(3,4,4+i);
                    data = orig_nh;
                    nan_ins = find(isnan(data));
                    
                    [pdata, pcax, pc] = ...
                        prep_data_land_nan_minmax(data, land_ins, nan_ins, ...
                        dmax, dmin, m,[.5 .5 .5],[1 1 1]);
                    
                    imagesc(pdata);  colormap(pc);    caxis(pcax);
                    axis tight; axis off;   thincolorbar;
                    title(tit)
                    
                    subplot(3,4,8+i);
                    dmax=10;
                    dmin = -10
                    
                    data = tmp(llc_all_to_nh_psg_mapping) - orig_nh;
                    nan_ins = find(isnan(data));
                    
                    [pdata, pcax, pc] = ...
                        prep_data_land_nan_minmax(data, land_ins, nan_ins, ...
                        dmax, dmin, m,[.5 .5 .5],[1 1 1]);
                    
                    imagesc(pdata);  colormap(pc);    caxis(pcax);
                    axis tight; axis off;   thincolorbar;
                    title(tit)
                end
                tightfig;
                A='top row mapped, mid row orig last row new-orig';
                B=ice_files_n(f).name;
                suptitle({A,B})
                
                tit = ['noaa_nsidc_to_llc' num2str(llcN) '_nh_' B]
                paperx=12;papery=12;prep_figure_for_print;
                cd(fig_dir);
                print(gcf,'-dpng','-r100',[tit '.png']);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% sh plotting
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                close all;
                figure(13);clf;set(0,'DefaultTextInterpreter','none');
                land_psg = hf0_N(llc_all_to_sh_psg_mapping);
                land_ins = find(land_psg ==0);
                
                for i = 1:4
                    switch i
                        case 1
                            tmp = global_g_merged_mapped;
                            orig_sh = sh_g_merged;
                            tit = 'goddard merged';
                        case 2
                            tmp = global_nt_mapped;
                            orig_sh = sh_g_nt;
                            tit = 'goddard nt';
                        case 3
                            tmp = global_bt_mapped;
                            orig_sh = sh_g_bt;
                            tit = 'goddard bt';
                        case 4
                            tmp = global_cdr_mapped;
                            orig_sh = sh_cdr;
                            tit = 'cdr';
                    end
                    dmin=0;
                    dmax=100;
                    m=100;
                    subplot(3,4,i)
                    tmp(dry_ins_N_k{1})=-1;
                    
                    data = tmp(llc_all_to_sh_psg_mapping);
                    nan_ins = find(isnan(data));
                    [pdata, pcax, pc] = ...
                        prep_data_land_nan_minmax(data, land_ins, nan_ins, ...
                        dmax, dmin, m,[.5 .5 .5],[1 1 1]);
                    
                    imagesc(pdata);  colormap(pc);    caxis(pcax);
                    axis tight; axis off;   thincolorbar;
                    title(tit)
                    
                    subplot(3,4,4+i);
                    data = orig_sh;
                    nan_ins = find(isnan(data));
                    [pdata, pcax, pc] = ...
                        prep_data_land_nan_minmax(data, land_ins, nan_ins, ...
                        dmax, dmin, m,[.5 .5 .5],[1 1 1]);
                    
                    imagesc(pdata);  colormap(pc);    caxis(pcax);
                    axis tight; axis off;   thincolorbar;
                    title(tit);
                    
                    subplot(3,4,8+i);
                    dmax = 10;
                    dmin = -10;
                    
                    data = tmp(llc_all_to_sh_psg_mapping) - orig_sh;
                    nan_ins = find(isnan(data));
                    [pdata, pcax, pc] = ...
                        prep_data_land_nan_minmax(data, land_ins, nan_ins, ...
                        dmax, dmin, m,[.5 .5 .5],[1 1 1]);
                    
                    imagesc(pdata);  colormap(pc);    caxis(pcax);
                    axis tight; axis off;   thincolorbar;
                    title(tit)
                end % plotting
                
                tightfig;
                A='top row mapped, mid row orig, bot row new-orig';
                B=ice_files_s(f).name;
                suptitle({A,B})
                
                tit = ['noaa_nsidc_to_llc' num2str(llcN) '_sh_' B];
                paperx=12;papery=8;prep_figure_for_print;
                cd(fig_dir)
                print(gcf,'-dpng','-r100',[tit '.png']);
            end % makefigs
            %%
            
            NOAA_NSIDC_DAILY_MAPPED_TO_LLC.goddard_merged(:,:,f) = ...
                global_g_merged_mapped;
            
            NOAA_NSIDC_DAILY_MAPPED_TO_LLC.goddard_merged_wet(:,:,f) = ...
                global_g_merged_wet_mapped;
            
            NOAA_NSIDC_DAILY_MAPPED_TO_LLC.goddard_nt(:,:,f) = ...
                global_nt_mapped;
            
            NOAA_NSIDC_DAILY_MAPPED_TO_LLC.goddard_nt_wet(:,:,f) = ...
                global_nt_wet_mapped;
            
            NOAA_NSIDC_DAILY_MAPPED_TO_LLC.goddard_bt(:,:,f) = ...
                global_bt_mapped;
            
            NOAA_NSIDC_DAILY_MAPPED_TO_LLC.goddard_bt_wet(:,:,f) = ...
                global_bt_wet_mapped;
            
            NOAA_NSIDC_DAILY_MAPPED_TO_LLC.cdr(:,:,f) = ...
                global_cdr_mapped;
            
            NOAA_NSIDC_DAILY_MAPPED_TO_LLC.cdr_wet(:,:,f) = ...
                global_cdr_wet_mapped;
        end %day
        %%
        
        cd(output_dir)
        fname = ['NOAA_NSIDC_DAILY_MAPPED_TO_LLC' num2str(llcN) '_' num2str(year)]
        save(fname, 'NOAA_NSIDC_DAILY_MAPPED_TO_LLC','-v7.3');
        
    end % year
end
