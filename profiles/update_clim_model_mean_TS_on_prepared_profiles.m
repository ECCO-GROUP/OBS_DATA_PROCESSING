function [MITprof] = ...
    update_clim_model_mean_TS_on_prepared_profiles(...
    update_clim_run_codes, write_profile_to_netcdf, make_figs, MITprof)


% This script updates ECCO profile.nc files with climatology or model mean
% fields and sends the updated profile files to be written in netcdf format
% fields can be mapped to prof_{T,S}estim, prof_{T,S}clim, or prof_{T,S}_model_mean

% Filename; update_clim_model_mean_TS_on_prepared_profiles.m
%  ** former filename :update_climatologyTS_on_prepared_profiles.m
% Date Created: 2014-04-29
% Last Modified: 2016-08-16
%

% notes:
%    previous version changed sigma and climatology at the same time
%    7/28/14 : changed so climatology is interpolated in the vertical.
%    fixed how files are loaded in, check for steppiness in vertical
%    interpl
%    2014-08-05: add make_roto_dir
%    2015-02-09: update fDataBase to reflect new argo convenction
%    2015-04-23: added ability to put in clim as clim field  and model mean
%    as model mean field
%    2015-04-29 : cleaned up file and made it clearer what it is each run
%    does
%    2015-05-05 : changed 3d interpolation to 2d
%    2015-05-06 : fixed missing transpose of returned 2d matrix
%     2015-06-17 : added exist file check
%    2016-03-08  : added file list as a field in the update_clim_run_code.
%    2016-03-10  : added some debugging figures to ensure correct mapping.
%    2016-04-08  : changed so that I could just use an MITprof file in
%    memory.
%    2016-08-11  : added check to T and Sclim  and T and Sestim,  don't plot
%    them if them don't exist!
%    2016-08-16  : made the whole thing a function, now you pass the run
%    codes and wehther you want to saave to file.  
%    2018-05-10  : added usage text


                    
% usage
%    update_clim_run_codes   : set of 'run codes'
%    write_profile_to_netcdf : self explanatory boolean
%    make_figs               : self explanatory boolean
%    MITprof                 : (optional) the MITprof to operate on

%%%
set(0,'DefaultTextInterpreter','none');
%%
make_root_dir;

%%---------------------------------------------
%% SET INPUT PARAMETERS
%%---------------------------------------------
fillVal=-9999;checkVal=-9000;


%% climatolgy dir and filenames
%

source_dir = ''
f_source_ST={}

%
%% source field type
% 0 == annual mean ** calc from 12 month climatology
% 1 == monthly
% 2 == monthly anomaly from annual mean
% 3 == annual mean from model output as prof_Tmodel_mean

%% target_field 
% 1  : as prof_Xestim
% 2  : as prof_Tclim/Sclim
% 3  : as prof_Tmodel_mmean/Smodel_mean

%% this parameter determines whether to extend the last good value in
%  the source field to the through land to the deepest cell.

% extend_last_good_mean_value_to_depth
extend_last_good_mean_value_to_depth = 0;

% 0 if MITprof is pass

if nargin == 3
    reload_MITprof = 1;
else
    reload_MITprof = 0;
end

        
deg2rad = pi/180.0;


%%---------------------------------------------
%%  BEGIN THE PROGRAM
%%---------------------------------------------

%%---------------------------------------------
%%  Init profile package
%%---------------------------------------------
MITprof_path;

%%---------------------------------------------
%%  Read in grid for source
%%---------------------------------------------

if make_figs
    close all;
end

for rr = 1:length(update_clim_run_codes)
    
    reload_MITprof = 1;
    update_clim_run_code = update_clim_run_codes{rr}
    
    switch update_clim_run_code
            ['apply monthly clim to mitprof variable']
            rootdir = []
            input_dir = [];
            
            file_suffix_in  = [];
            file_suffix_out = [];
            
            output_dir = input_dir;
            
            % the climatology we want to put in is the monthly climatology
            source_field = 1
            
            % add to prof_Tclim Sclim
            target_field = 2
            
            % the climatology
            source_dir = ['/home/ifenty/data/observations/TS_climatology/ECCO_v4r2/']
            f_source_ST={['/S_monthly_woa09'], ['/T_monthly_woa09']};
            
            reload_MITprof = 0;
            
            extend_last_good_mean_value_to_depth = 0;
            ilists = 1;
        

            %%%%%%%%%%%%%%
            % 2018-05-08 %
            %%%%%%%%%%%%%%

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % ALL
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case  'APPLY_WOA13V2_MONTHLY_CLIM_TO_WOD13_20180508'
            ['apply monthly clim to mitprof variable']
            rootdir = ['/home/ifenty/data/observations/insitu/NODC/NODC_20180508/all/']
            input_dir = [rootdir 'step_00_update_tile_and_prof_points'];
            
            file_suffix_in  = ['step_00.nc'];
            file_suffix_out = ['step_01.nc']
            
            output_dir =  [rootdir 'step_01_monthly_clim']
            
            % the climatology we want to put in is the monthly climatology
            source_field = 1
            
            % add to prof_Tclim Sclim
            target_field = 2
            
            cd(input_dir);
            f = dir(['*nc*']);
            
            for i = 1:length(f)
                fDataIn{i}  =   f(i).name;
                fDataOut{i} =  [f(i).name(1:end-length(file_suffix_in)) file_suffix_out];
            end
            ilists = 1:length(f);
            
            % the climatology
            source_dir = ['/ian3/ifenty/data/observations/TS_climatology/WOA_2013_V2/1995-2014 merged']
            f_source_ST={['/WOA13_v2_TS_clim_merged.mat'],['/WOA13_v2_TS_clim_merged.mat']};

            
            %%%%%%%%%%%%%%
            % 2016-08-01 %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        case  'APPLY_MONTHLY_CLIM_TO_WOD13_CTD_20160801'
            ['apply monthly clim to mitprof variable']
            %%
            rootdir = ['/ian4/ifenty/data/observations/insitu/NODC/NODC_20160801/']
            input_dir = [rootdir 'TMPb_sigma_v2'];
            
            file_suffix_in  = ['step_08_sigma_v3.nc'];
            file_suffix_out = ['step_08_sigma_TS_clim.nc'];
            
            output_dir = input_dir;
            
            % the climatology we want to put in is the monthly climatology
            source_field = 1
            
            % add to prof_Tclim Sclim
            target_field = 2
            
            cd(input_dir);
            f = dir(['*step_08*nc*']);
            
            for i = 1:length(f)
                fDataIn{i}  =   f(i).name;
                fDataOut{i} =  [f(i).name(1:end-length(file_suffix_in)) file_suffix_out];
            end
            
            ilists = 1:length(f);
            
            % the climatology
            source_dir = ['/home/ifenty/data/observations/TS_climatology/ECCO_v4r2/']
            f_source_ST={['/S_monthly_woa09'], ['/T_monthly_woa09']};
        
        case  'APPLY_MONTHLY_CLIM_TO_WOD13_PROFILES'
            %%
            ['apply monthly clim to mitprof variable']
            rootdir = ['/ian4/ifenty/data/observations/insitu/NODC/WOD13_AUG2016_CTD_GLD_APB_XBT_CSVFORMAT/']
            input_dir = [rootdir 'APB_PROC'];
            
            file_suffix_in  = ['step_00.nc'];
            file_suffix_out = ['step_02.nc'];
            
            output_dir = input_dir;
            
            % the climatology we want to put in is the monthly climatology
            source_field = 1
            
            % add to prof_Tclim Sclim
            target_field = 2
            
            cd(input_dir);
            f = dir(['*step_00*']);
            
            for i = 1:length(f)
                fDataIn{i}=   f(i).name;
                fDataOut{i} =[f(i).name(1:end-length(file_suffix_in)) file_suffix_out];
            end
            
            ilists = 1:length(f);
            
            % the climatology
            source_dir = ['/home/ifenty/data/observations/TS_climatology/ECCO_v4r2/']
            f_source_ST={['/S_monthly_woa09'], ['/T_monthly_woa09']};
            
            
        case 'APPLY_MODEL_MEAN_TO_WOD13_PROFILES'
            %%
            ['apply monthly clim to mitprof variable']
            rootdir = ['/ian4/ifenty/data/observations/insitu/NODC/WOD13_AUG2016_CTD_GLD_APB_XBT_CSVFORMAT/']
            input_dir = [rootdir 'APB_PROC'];
            
            file_suffix_in  = ['step_02.nc'];
            file_suffix_out = ['step_03.nc'];
            
            output_dir = input_dir;
            
            source_field = 3; % the 'climatology is actually the model time-mean output
            target_field = 3  % add Tmodel_mean, Smodel_mean;
            
            cd(input_dir);
            f = dir(['*step_02*']);
            
            for i = 1:length(f)
                fDataIn{i}=   f(i).name;
                fDataOut{i} =[f(i).name(1:end-length(file_suffix_in)) file_suffix_out];
            end
            
            ilists = 1:length(f);
            
            % th model mean
            source_dir = ['/ian2/ifenty/data/model_output/ECCOv4/r2_iteration_26/model_time_mean']
            f_source_ST={['/sbarmean.data'],['/tbarmean.data']};

                        
      
            
        case  'ECCO_v5_update_to_Tclim_Sclim'
            rootdir = ['/home/ifenty/data/observations/insitu/ECCO_v5r1/']
            input_dir = [rootdir '/20150616_llc270']
            
            file_suffix_in  = '_llc270_20150616.nc';
            file_suffix_out = '_llc270_20150616b.nc';
            
            output_dir = input_dir;
            
            % the climatology we want to put in is the time-mean of the 12
            % month climatology
            source_field = 1
            
            % add to prof_Testim, prof_Sestim
            target_field = 2
            
            % the climatology
            source_dir = ['/home/ifenty/data/observations/TS_climatology/ECCO_v4r2/']
            f_source_ST={['/S_monthly_woa09'], ['/T_monthly_woa09']};
            
            extend_last_good_mean_value_to_depth = 0;
            
        case  'ECCO_v5_update_to_Testim_Sestim'
            rootdir = ['/home/ifenty/data/observations/insitu/ECCO_v5r1/']
            input_dir = [rootdir '/20150616_llc270']
            
            file_suffix_in  = '_llc270_20150616b.nc';
            file_suffix_out = '_llc270_20150616c.nc';
            
            % the climatology we want to put in is the time-mean of the 12
            % month climatology
            source_field = 1
            target_field =1  % add prof_Testim, prof_Sestim;
            
            % the climatology
            source_dir = ['/home/ifenty/data/observations/TS_climatology/ECCO_v4r2/']
            f_source_ST={['/S_monthly_woa09'], ['/T_monthly_woa09']};
            
            extend_last_good_mean_value_to_depth = 0;
            
            
        case  'ECCO_v4_i011_add_clim_mean_TS'
            rootdir = ['/ian2/ifenty/data/model_output/ECCOv4/'];
            input_dir = [rootdir '/r2_iteration_11/profiles/']
            file_suffix_in  = '_AWPP_llc90_20150206_model.nc';
            file_suffix_out = '_step_1.nc';
            
            output_dir = [rootdir '/r2_iteration_11/profiles/add_model_mean_and_clim']
            
            % the climatology we want to put in is the time-mean of the 12
            % month climatology
            source_field = 0
            target_field =2  % add Tclim, Sclim
            
            % the climatology
            source_dir = ['/home/ifenty/data/observations/TS_climatology/ECCO_v4r2/']
            
            f_source_ST={['/S_monthly_woa09'], ['/T_monthly_woa09']};
            
            extend_last_good_mean_value_to_depth = 0;
            
        case 'ECCO_v4_i011_add_model_mean_TS'
            rootdir = ['/ian2/ifenty/data/model_output/ECCOv4/'];
            input_dir = [rootdir '/r2_iteration_11/profiles/add_model_mean_and_clim']
            file_suffix_in  = '_step_1.nc';
            file_suffix_out = '_step_2.nc';
            
            output_dir = input_dir;
            
            source_field = 3; % the 'climatology is actually the model time-mean output
            target_field = 3  % add Tmodel_mean, Smodel_mean;
            
            % th model mean
            source_dir = [rootdir '/r2_iteration_11/model_time_mean'];
            f_source_ST={['/sbarmean.data'],['/tbarmean.data']};
            
            extend_last_good_mean_value_to_depth = 0;
            
            
        case  'ECCO_v4_i026_add_clim_mean_TS'
            rootdir = ['/ian2/ifenty/data/model_output/ECCOv4/'];
            input_dir = [rootdir '/r2_iteration_26/profiles/']
            
            output_dir = input_dir;
            
            % the climatology we want to put in is the time-mean of the 12
            % month climatology
            source_field = 0
            target_field =2  % add Tclim, Sclim
            
            % the climatology
            source_dir = ['/home/ifenty/data/observations/TS_climatology/ECCO_v4r2/']
            f_source_ST={['/S_monthly_woa09'], ['/T_monthly_woa09']};
            
            extend_last_good_mean_value_to_depth = 0;
            
        case 'ECCO_v4_i026_add_model_mean_TS'
            rootdir = ['/ian2/ifenty/data/model_output/ECCOv4/'];
            input_dir = [rootdir '/r2_iteration_26/profiles/add_model_mean_and_clim']
            file_suffix_in  = '_step_1b.nc';
            file_suffix_out = '_step_2b.nc';
            
            output_dir = input_dir;
            
            source_field = 3; % the 'climatology is actually the model time-mean output
            target_field = 3  % add Tmodel_mean, Smodel_mean;
            
            % th model mean
            source_dir = [rootdir '/r2_iteration_26/model_time_mean'];
            f_source_ST={['/sbarmean.data'],['/tbarmean.data']};
            
            extend_last_good_mean_value_to_depth = 0;
            
            
        case  'ECCO_v4_feb2016_add_clim_mean_TS'
            rootdir = ['/home/ifenty/data/observations/insitu/ECCO_v4r2/llc90_20160308']
            input_dir = rootdir;
            
            file_suffix_in  = '_feb2016_llc90_step_07_20160308.nc'
            file_suffix_out = '_feb2016_llc90_step_08_20160308.nc'
            
            output_dir = input_dir;
            
            % the climatology we want to put in is the time-mean of the 12
            % month climatology
            source_field = 0
            target_field =2  % add Tclim, Sclim
            
            % the climatology
            source_dir = ['/home/ifenty/data/observations/TS_climatology/ECCO_v4r2/']
            f_source_ST={['/S_monthly_woa09'], ['/T_monthly_woa09']};
            
            extend_last_good_mean_value_to_depth = 0;
            
            yrs = 1995:2016
            for i = 1:length(yrs)
                fDataBase{i} = ['argo_' num2str(yrs(i))];
            end
            
            clear fDataIn* fDataOut;
            for i = 1:length(fDataBase)
                fDataIn{i}  = [fDataBase{i} file_suffix_in]
                fDataOut{i} = [fDataBase{i} file_suffix_out]
            end
            
            ilists = 1:length(fDataBase)
            
        case 'ECCO_v4_feb2016_add_model_mean_TS'
            rootdir = ['/ian2/ifenty/data/model_output/ECCOv4/'];
            
            input_dir = ['/home/ifenty/data/observations/insitu/ECCO_v4r2/llc90_20160308']
            
            file_suffix_in  = '_feb2016_llc90_step_08_20160308.nc'
            file_suffix_out = '_feb2016_llc90_step_09_20160308.nc'
            
            output_dir = input_dir;
            
            source_field = 3; % the 'climatology is actually the model time-mean output
            target_field = 3  % add Tmodel_mean, Smodel_mean;
            
            % th model mean
            source_dir = [rootdir '/r2_iteration_26/model_time_mean'];
            f_source_ST={['/sbarmean.data'],['/tbarmean.data']};
            
            extend_last_good_mean_value_to_depth = 0;
            
            yrs = 2014:2016
            for i = 1:length(yrs)
                fDataBase{i} = ['argo_' num2str(yrs(i))];
            end
            
            clear fDataIn* fDataOut;
            for i = 1:length(fDataBase)
                fDataIn{i}  = [fDataBase{i} file_suffix_in]
                fDataOut{i} = [fDataBase{i} file_suffix_out]
            end
            
            ilists = 1:length(fDataBase)
            
    end
    
    
    mkdir(output_dir)
    
    
    % source_field
    %% source field
    % 0 == annual mean ** calc from 12 month climatology
    % 1 == monthly
    % 2 == monthly anomaly from annual mean
    % 3 == annual mean from model output as prof_Tmodel_mean
    
    %% target_field
    % 1  : as prof_Xestim
    % 2  : as prof_Tclim/Sclim
    % 3  : as prof_Tmodel_mmean/Smodel_mean
    
    switch source_field
        case {0,1,2}
            source_name = 'climatology'
        case 3
            source_name = 'model mean'
    end
    
    switch target_field
        case 1
            target_name = 'estim'
        case 2
            target_name = 'climatology'
        case 3
            target_name = 'model mean'
    end
    
    
    %%---------------------------------------------
    %%  Read in climatology fields.
    %%---------------------------------------------
    ['reading clim or model field']
    switch source_field
        
        case {0,1,2} % llc90 grid climatology
            
           % MITprof.prof_Tclim = MITprof.prof_T.*NaN;
           % MITprof.prof_Sclim = MITprof.prof_T.*NaN;
            
            % salt
            load_llc90_grid
            cd(source_dir)
            
            source_S = readbin(f_source_ST{1},[llcN 13*llcN 50 12],1,'real*4',0,'ieee-be');
            % theta
            source_T = readbin(f_source_ST{2},[llcN 13*llcN 50 12],1,'real*4',0,'ieee-be');
            
            % NaN out all of the land values in each month.
            for m = 1:12
                tmp = source_S(:,:,:,m);
                tmp(dry_ins_90) = NaN;
                source_S(:,:,:,m) = tmp;
                
                tmp = source_T(:,:,:,m);
                tmp(dry_ins_90) = NaN;
                source_T(:,:,:,m) = tmp;
            end
            
            
            %% calculate annual mean and monthly anomaly for the climatology
            % annual mean
            
            source_Sm = mean(source_S,4);
            source_Tm = mean(source_T,4);
            
            if source_field == 2
                % calc monthly anomaly (month - annual mean)
                for m = 1:12
                    source_Sa(:,:,:,m) = source_S(:,:,:,m) - source_Sm;
                    source_Ta(:,:,:,m) = source_T(:,:,:,m) - source_Tm;
                end % model mean
            end
            
            %% MODEL OUTPUT
        case {3}
            
            switch llcN
                case 90
                    load_llc90_grid
                    nz=50;
                    hFacC =hFacC_90;
                case 270
                    load_ll270_grid
                    nz=50;
                    hFacC =hFacC_270;
                otherwise
                    ['check your grid homie']
                    nz = -99999
                    
            end
            
            cd(source_dir)
            
            source_Sm = readbin(f_source_ST{1},[llcN 13*llcN nz 1],1,'real*4',0,'ieee-be');
            % theta
            source_Tm = readbin(f_source_ST{2},[llcN 13*llcN nz 1],1,'real*4',0,'ieee-be');
            
            
            [' if the model output isnt llc90 this much change!!']
            
            source_Sm(find(hFacC ==0))=NaN;
            source_Tm(find(hFacC ==0))=NaN;
            
    end
    
    stats(source_Sm)
    stats(source_Tm)
    
    
    %%---------------------------------------------
    %% Prepare the nearest neighbor mapping
    %%---------------------------------------------
    ['preparing the nearest neighbor']
    cd(llc90_grid_dir)
    if exist('F_llc90_ALL_INS_SURF_XYZ_to_INDEX.mat','file')
        ['loading F']
        load('F_llc90_ALL_INS_SURF_XYZ_to_INDEX');
        ['loaded F']
    else
        make_F_llc90_ALL_INS_SURF_XYZ_to_INDEX;
    end
    
    
    % verify that our little trick works in 4 parts of the earth
    for i = 1:4
        switch i
            case 1
                test_lat = 56;
                test_lon = -40;
            case 2
                test_lat = 60;
                test_lon = 10;
            case 3
                test_lat = -60;
                test_lon = -120;
            case 4
                test_lat = -69;
                test_lon = 60;
        end
        [test_x, test_y, test_z] = ...
            sph2cart(test_lon*deg2rad, test_lat*deg2rad, 1);
        test_ind = F_llc90_ALL_INS_SURF_XYZ_to_INDEX(test_x, test_y, test_z);
        
        [X_90(test_ind) Y_90(test_ind) Z_90(test_ind) test_x test_y test_z]
        [lat_90(test_ind) lon_90(test_ind) test_lat test_lon]
    end
    
    %%---------------------------------------------
    %% Read and process the profile files
    %%---------------------------------------------
    
    
    if ~reload_MITprof
        ilists = 1
        
        ['Using  MITprof is already in memory']
        whos MITprof
    end
    
    % Cycle through ilists
    for ilist = ilists
        
        if reload_MITprof
            cd(input_dir)
            
            fDataIn{ilist}
            
            clear fileOut MITprof
            
            if exist(fDataIn{ilist})
                MITprof = MITprof_read(fDataIn{ilist});
            else
                MITprof = []
                ['your stupid file does not exist ' fDataIn{ilist}]
            end
        end
        
        if unique(MITprof.prof_depth == -9999)
            MITprof.prof_depth = argo_profile_depths;
        end
        
        num_profs = length(MITprof.prof_lat)
        num_prof_depths = length(MITprof.prof_depth);
        
        
        clear good ins prof_Testim_orig prof_Sestim_orig;
        
        %%
        % bad data = profX_flag > 0 --> 0 weight, valid climatology value.
        %  no data => profX == -9999, valid value in climatology, 0 in weight
        
        % determine the month for every profile
        tmp = num2str(MITprof.prof_YYYYMMDD);
        prof_month = str2num(tmp(:,5:6));
        
        ['mapping profiles to llc grid']
        [prof_x, prof_y, prof_z] = ...
            sph2cart(MITprof.prof_lon*deg2rad, MITprof.prof_lat*deg2rad, 1);
        
        % map a llc90 grid index to each profile.
        prof_llc90_cell_index = ...
            F_llc90_ALL_INS_SURF_XYZ_to_INDEX(prof_x, prof_y, prof_z);
        %%
        unique_prof_index = unique(prof_llc90_cell_index);
        
        % go through each z level in the profile array
        
        clear tmp*;
        tmp_source_field_T = zeros(num_profs, num_prof_depths);
        tmp_source_field_S = zeros(num_profs, num_prof_depths);
        
        switch  source_field
            case {0,3} % annual means climatology
                %%
                [' annual mean']
                nx = size(source_Tm,1);
                ny = size(source_Tm,2);
                nz = size(source_Tm,3);
                new_nz = length(MITprof.prof_depth);
                
                source_Tm_tmp = reshape(source_Tm, [size(source_Tm,1)*size(source_Tm,2) ...
                    size(source_Tm,3)]);
                source_Tm_here = source_Tm_tmp(prof_llc90_cell_index,:)';
                % rows = depth
                % cols = different profiles
                
                source_Sm_tmp = reshape(source_Sm, [size(source_Sm,1)*size(source_Sm,2) ...
                    size(source_Sm,3)]);
                source_Sm_here = source_Sm_tmp(prof_llc90_cell_index,:)';
                
                % interp cTm, cSm to the new vertical levels
                tmp_source_field_T = interp_2D_to_arbitrary_z_levels(source_Tm_here, ....
                    z_cen_90, MITprof.prof_depth, ...
                    extend_last_good_mean_value_to_depth)';
                %%
                
                ['size source_Sm here']
                size(source_Sm_here)
                ['size z_cen_90']
                size(z_cen_90)
                ['size mitprof.prof_depth']
                size(MITprof.prof_depth)
                
                tmp_source_field_S = interp_2D_to_arbitrary_z_levels(source_Sm_here, ....
                    z_cen_90, MITprof.prof_depth, ...
                    extend_last_good_mean_value_to_depth)';
                
                if make_figs
                    figure(101);clf;
                    nx=size(source_Tm_here,2);
                    subplot(221);mypcolor(1:size(source_Tm_here,2), -z_cen_90, source_Tm_here);
                    cax=caxis;
                    axis([0 nx -5000 0]);
                    title({'source Tm here', size(source_Tm_here)});thincolorbar;
                    
                    subplot(222);
                    mypcolor(1:size(source_Tm_here,2), -MITprof.prof_depth, tmp_source_field_T);
                    title(size(tmp_source_field_T))
                    axis([0 nx -5000 0]);
                    caxis(cax);
                    title({'post-interp source Tm ' , size(tmp_source_field_T)});thincolorbar;
                    
                    subplot(223);mypcolor(1:size(source_Sm_here,2), -z_cen_90, source_Sm_here);
                    caxis([30 36.25]);
                    axis([0 nx -5000 0]);
                    title({'source Sm here', size(source_Sm_here)});thincolorbar;
                    
                    subplot(224);
                    mypcolor(1:size(source_Sm_here,2), -MITprof.prof_depth, tmp_source_field_S);
                    title(size(tmp_source_field_S))
                    axis([0 nx -5000 0]);
                    caxis([30 36.25]);
                    title({'post-interp source Sm' , size(tmp_source_field_S)});thincolorbar;
                end
                
                %%
                ['fin interp']
                %clear cTm_tmp cTm_here cSm_tmp cSm_here;
                
            case 1 % standard monthly value
                %%
                [' std clim']
                nx = size(source_T,1);
                ny = size(source_T,2);
                nz = size(source_T,3);
                new_nz = length(MITprof.prof_depth);
                
                clear cT_new_z cS_new_z source_T_new_* source_S_new*;
                for m = 1:12
                    ['looping through months ' padzero(m,2)]
                    
                    % interp cT to the new vertical levels
                    
                    source_Tm = source_T(:,:,:,m);
                    source_Sm = source_S(:,:,:,m);
                    
                    source_Tm_tmp = reshape(source_Tm, ...
                        [size(source_Tm,1)*size(source_Tm,2) size(source_Tm,3)]);
                    
                    source_Tm_here = source_Tm_tmp(unique_prof_index,:)';
                    
                    % rows = depth
                    % cols = different profiles
                    
                    source_Sm_tmp = reshape(source_Sm, ...
                        [size(source_Sm,1)*size(source_Sm,2) size(source_Sm,3)]);
                 
                    source_Sm_here = source_Sm_tmp(unique_prof_index,:)';
                    
                    
                    % interp   to the new vertical levels
                    %%
                    % must be in x-z  form
                    %size(source_Tm_here')
                    tmp_source_field_T = interp_2D_to_arbitrary_z_levels(source_Tm_here, ....
                        z_cen_90, MITprof.prof_depth, ...
                        extend_last_good_mean_value_to_depth)';
                    %size(tmp_source_field_T)
                    
                    %%
                    tmp_source_field_S = interp_2D_to_arbitrary_z_levels(source_Sm_here, ....
                        z_cen_90, MITprof.prof_depth, ...
                        extend_last_good_mean_value_to_depth)';
                    %size(tmp_source_field_S)
                    %%
                    
                    tmp_source_field_T =tmp_source_field_T';
                    tmp_source_field_S =tmp_source_field_S';
                    %%
                    source_T_new_z(:,:,m) = tmp_source_field_T;
                    source_S_new_z(:,:,m) = tmp_source_field_S;
                end
                clear tmp_source_field_T;
                
                %% pull out the profs that match each unique pci
                tic
                clear tmp_source_field_*;
                tmp_source_field_T = MITprof.prof_T.*NaN;
                tmp_source_field_S = MITprof.prof_S.*NaN;
                
                for uci = 1:length( unique_prof_index)
                    uc = unique_prof_index(uci);
                    % these are teh profiles that have the unique cell
                    % index uc
                    
                    ins = find(ismembc2(prof_llc90_cell_index, uc));
             
                    % loop through these profiles one at a time
                    for ii = 1:length(ins)
                        % this is the profile index
                        ins_i = ins(ii);
                        % this is its month
                        pm = prof_month(ins_i);
                        % this is its clim field
                        tmp_source_field_T(ins_i,:) = ...
                            squeeze(source_T_new_z(uci, :, pm));
                        tmp_source_field_S(ins_i,:) = ...
                            squeeze(source_S_new_z(uci, :, pm));
                    end
                end
                toc
            case 2 % monthly anomaly against annual mean
                ['monthly anomaly']
                nx = size(source_T,1);
                ny = size(source_T,2);
                nz = size(source_T,3);
                new_nz = length(MITprof.prof_depth);
                
                clear source_*a_new_z;
                for m = 1:12
                    % interp cT to the new vertical levels
                    tmp = interp_3D_to_arbitrary_z_levels(source_Ta(:,:,:,m), ....
                        z_top_90, z_cen_90, z_bot_90, MITprof.prof_depth);
                    tmp = reshape(tmp, [nx*ny new_nz]);
                    tmp = tmp(prof_llc90_cell_index,:);
                    
                    source_Ta_new_z(:,:,m) = tmp;
                    
                    tmp  = interp_3D_to_arbitrary_z_levels(source_Sa(:,:,:,m), ....
                        z_top_90, z_cen_90, z_bot_90, MITprof.prof_depth);
                    tmp = reshape(tmp, [nx*ny new_nz]);
                    tmp = tmp(prof_llc90_cell_index,:);
                    source_Sa_new_z(:,:,m) = tmp;
                end
                
                %%
                % this little bit of magic is how we pull the correct months out of
                % the climatology.
                for z = 1:new_nz
                    tmpT =  squeeze(source_Ta_new_z(:,z,:));
                    tmpS =  squeeze(source_Sa_new_z(:,z,:));
                    source_Index = sub2ind(size(tmpT), 1:num_profs,prof_month');
                    
                    tmp_source_field_T(:,z) = tmpT(source_Index);
                    tmp_source_field_S(:,z) = tmpS(source_Index);
                end
                ['done mapping in z']
                
        end % climatology type
        
        %%
        if size(tmp_source_field_T,1) ~= size(MITprof.prof_T,1)
            tmp_source_field_T = tmp_source_field_T';
            tmp_source_field_S = tmp_source_field_S';
        else
            ['no need to rotate clim field']
        end
        
        
        %%
        % replace climatology
        
        tmp_source_field_T(find(isnan(tmp_source_field_T)))= fillVal;
        tmp_source_field_S(find(isnan(tmp_source_field_S)))= fillVal;
        
        
        %%
        
        if make_figs
            figure(100);clf;
            
            subplot(241);
            imagesc(MITprof.prof_T);caxis([-2 25]);title('prof T')
            subplot(242);
            if isfield(MITprof,'prof_Tclim')
                imagesc(MITprof.prof_Tclim);caxis([-2 25]);title('Tclim');
            end
            
            subplot(243);
            if isfield(MITprof,'prof_Testim')
                imagesc(MITprof.prof_Testim);caxis([-2 25]);title('Testim')
            end
            subplot(244);
            
            imagesc(tmp_source_field_T);;caxis([-2 25]);
            title([' new T: ' target_name]);caxis([-2 25]);
            
            subplot(245);
            imagesc(MITprof.prof_S);caxis([33 36]);title('prof S')
            subplot(246);
            if isfield(MITprof,'prof_Sclim')
                imagesc(MITprof.prof_Sclim);caxis([33 36]);title('Sclim');
            end
            subplot(247);
            if isfield(MITprof,'prof_Testim')
                imagesc(MITprof.prof_Sestim);caxis([33 36]);title('Sestim')
            end
            subplot(248);
            imagesc(tmp_source_field_S);
            title([' new S: ' target_name]);caxis([33 36]);
        end
      
        switch target_field
            case 1
                MITprof.prof_Testim = tmp_source_field_T;
                if isfield(MITprof,'prof_S')
                    MITprof.prof_Sestim = tmp_source_field_S;
                end
                
            case 2
                MITprof.prof_Tclim = tmp_source_field_T;
                if isfield(MITprof,'prof_S')
                    MITprof.prof_Sclim = tmp_source_field_S;
                end
            case 3
                MITprof.prof_Tmodel_mean = tmp_source_field_T;
                if isfield(MITprof,'prof_S')
                    MITprof.prof_Smodel_mean = tmp_source_field_S;
                end
        end
        
        %%
        
        MITprof.prof_gci = prof_llc90_cell_index;
        
        %  Write output
        
        if write_profile_to_netcdf
            fileOut=[output_dir '/' fDataOut{ilist}];
            fprintf('%s\n',fileOut);
            
            write_profile_structure_to_netcdf(MITprof, fileOut);
        
        else
            ['I am not saving anything, the rest is up to you']
        end % write
    end % ilist
end % run code.

clear *source* *prof*estim_orig* 


%    %prof_Testim_orig = MITprof.prof_Testim;
        
        %if isfield(MITprof,'prof_Tweight')
        %    good_ins = find(MITprof.prof_Tweight > 0);
        %else
        %    good_ins = find(~isnan(MITprof.prof_T));
        %end
        
        %if length(good_ins) > 0
        %    stats(prof_Testim_orig(good_ins))
        %    stats(MITprof.prof_T(good_ins));
        %end
        
        %if isfield(MITprof,'prof_S')
        %    if isfield(MITprof,'prof_STweight')
        %        good_ins = find(MITprof.prof_Sweight > 0);
        %    else
        %        good_ins = find(~isnan(MITprof.prof_T));
        %    end
            
            %prof_Sestim_orig = MITprof.prof_Sestim;
            
            %if length(good_ins) > 0
            %    stats(prof_Sestim_orig(good_ins))
            %    stats(MITprof.prof_S(good_ins));
            %end
            
        %end
        
%                     %%
%                     tmp_source_field_T(ins',:,:) = squeeze(source_T_new_z(uci,:,:));
%                     %%
%                     tmp_source_field_S(ins',:,:) = squeeze(source_S_new_z(uci,:,:));
%                 end
%                 toc
%                 %%
%                 % this little bit of magic is how we pull the correct
%                 % months out of
%                 % the climatology.
%                 for z = 1:new_nz
%                     tmpT =  squeeze(source_T_new_z(:,z,:));
%                     tmpS =  squeeze(source_S_new_z(:,z,:));
%                     source_Index = sub2ind(size(tmpT), 1:num_profs,prof_month');
%                     
%                     tmp_source_field_T(:,z) = tmpT(source_Index);
%                     tmp_source_field_S(:,z) = tmpS(source_Index);
%                 end
%                 ['done mapping in z']
%                 %%
%                 