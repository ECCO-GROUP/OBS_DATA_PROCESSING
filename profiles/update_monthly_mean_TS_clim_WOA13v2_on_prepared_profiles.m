function [MITprof] = ...
    update_monthly_mean_TS_clim_WOA13v2_on_prepared_profiles(...
    run_code, write_profile_to_netcdf, make_figs)

%%
% This script updates ECCO profile.nc files with climatology from the WOA 
% 2013 v2 and sends the updated profile files to be written in netcdf format
% 
% The WOA climatology is in lat/lon at standard WOA levels
% profiles must also be at the standard WOA levels

% Filename; update_monthly_TS_clim_WOA13v2_on_prepared_profiles.m
%  ** former filename :update_clim_model_mean_TS_on_prepared_profiles.m
% Date Created: 2018-05-16
% Last Modified: 2018-06-05
%


% notes:
%    2018-06-01, got it to work ( I hope)
%    2018-06-05, potential temperature instead of in-situ temperature

                    
% usage
%    run_code                : the 'run code'
%    write_profile_to_netcdf : self explanatory boolean
%    make_figs               : self explanatory boolean
%    MITprof                 : (optional) the MITprof to operate on

%%%
set(0,'DefaultTextInterpreter','none');

make_root_dir;

%%---------------------------------------------
%% SET INPUT PARAMETERS
%%---------------------------------------------
fillVal=-9999;checkVal=-9000;


%% climatolgy dir and filenames
%

TS_clim_dir = ''
f_source_ST={}

%
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


%%%%%%%%%%%%%%
% 2018-05-08 %
%%%%%%%%%%%%%%
switch run_code
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ALL
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case  'APPLY_WOA13V2_MONTHLY_CLIM_TO_ALL_WOD13_20180508'
        ['apply monthly clim to mitprof variable']
        rootdir = ['/home/ifenty/data/observations/insitu/NODC/NODC_20180508/all/']
        input_dir = [rootdir 'step_00_update_tile_and_prof_points'];
        
        %%
        file_suffix_in  = ['step_00.nc'];
        file_suffix_out = ['step_01.nc']
        
        output_dir =  [rootdir 'step_01_monthly_clim/']
        
        cd(input_dir);
        f = dir(['*nc*']);
        
        for i = 1:length(f)
            fDataIn{i}  =   f(i).name;
            fDataOut{i} =  [f(i).name(1:end-length(file_suffix_in)) file_suffix_out];
        end
        ilists = 1:length(f);
        
        % the climatology
        TS_clim_dir = ['/ian3/ifenty/data/observations/TS_climatology/WOA_2013_V2/1995-2014 merged']
        
        %TS_clim_fname=['/WOA13_v2_TS_clim_merged.mat'];
        TS_clim_fname = ['/WOA13_v2_TS_clim_merged_with_potential_T.mat'];
        
        %%
     case  'APPLY_WOA13V2_MONTHLY_CLIM_TO_ALL_ARGO_20180612'
        ['apply monthly clim to mitprof variable']
        rootdir = ['/home/ifenty/data/observations/insitu/ARGO/from_gael_June2018/combined_by_latest']
        input_dir = [rootdir 'step_00_update_tile_and_prof_points'];
        
        %%
        file_suffix_in  = ['step_00.nc'];
        file_suffix_out = ['step_01.nc']
        
        output_dir =  [rootdir 'step_01_monthly_clim/']
        
        cd(input_dir);
        f = dir(['*nc*']);
        
        for i = 1:length(f)
            fDataIn{i}  =   f(i).name;
            fDataOut{i} =  [f(i).name(1:end-length(file_suffix_in)) file_suffix_out];
        end
        ilists = 1:length(f);
        
        % the climatology
        TS_clim_dir = ['/ian3/ifenty/data/observations/TS_climatology/WOA_2013_V2/1995-2014 merged']
        
        %TS_clim_fname=['/WOA13_v2_TS_clim_merged.mat'];
        TS_clim_fname = ['/WOA13_v2_TS_clim_merged_with_potential_T.mat'];
        
        %%
        
end

%%
cd(TS_clim_dir);
load(TS_clim_fname);

[lon_woam, lat_woam] = meshgrid(WOA_2013_v2_clim.lon.data, ...
    WOA_2013_v2_clim.lat.data);

if make_figs
    figure(1);clf;
    
    subplot(311);imagesc(lat_woam);axis xy;colorbar
    subplot(312);imagesc(lon_woam);axis xy;colorbar
    subplot(313);imagesc(squeeze(WOA_2013_v2_clim.potential_T_monthly(1,1,:,:)));colorbar;
    
    axis xy
end
%%

deg2rad=pi/180;
[X_woa, Y_woa, Z_woa] = sph2cart(lon_woam*deg2rad, lat_woam*deg2rad, 1);

AI = 1:length(X_woa(:));
AI = reshape(AI, size(X_woa));

if make_figs
    figure(3);clf;
    plot3(X_woa, Y_woa, Z_woa,'r.');axis equal;
end

% these are the x,y,z coordinates of all points in the climatology
xyz= [X_woa(:) Y_woa(:) Z_woa(:)];

cd(TS_clim_dir)

if exist('Scattered_Interpolant_WOA_xyz_element_number.mat','file')
    'loading scattered interpolant'
    load Scattered_Interpolant_WOA_xyz_element_number
else
    'making scattered interpolant'
    F  = scatteredInterpolant(xyz, AI(:),'nearest');
    save('Scattered_Interpolant_WOA_xyz_element_number.mat', 'F');
end

['now I have the interpolant']

%%
    
% verify that our little trick works in 4 parts of the earth
for i = 1:4
    switch i
        case 1
            test_lat = 56;i
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
    
    test_ind = F(test_x, test_y, test_z);
    
    [X_woa(test_ind) Y_woa(test_ind) Z_woa(test_ind) test_x test_y test_z]
    [lat_woam(test_ind) lon_woam(test_ind) test_lat test_lon]
    
end
%%

% Cycle through ilists
for ilist = ilists
    
    cd(input_dir)
    
    fDataIn{ilist}
    
    clear fileOut MITprof
    
    if exist(fDataIn{ilist})
        MITprof = MITprof_read(fDataIn{ilist});
    else
        MITprof = []
        ['your stupid file does not exist ' fDataIn{ilist}]
    end
    
    num_profs = length(MITprof.prof_lat);
    num_prof_depths = length(MITprof.prof_depth);
    
    
    %%
    % bad data = profX_flag > 0 --> 0 weight, valid climatology value.
    %  no data => profX == -9999, valid value in climatology, 0 in weight
    
    % determine the month for every profile
    tmp = num2str(MITprof.prof_YYYYMMDD);
    prof_month = str2num(tmp(:,5:6));
    
    %%
    ['mapping profiles to x,y,z']
    [prof_x, prof_y, prof_z] = ...
        sph2cart(MITprof.prof_lon*deg2rad, MITprof.prof_lat*deg2rad, 1);
    
    % map a climatology grid index to each profile.
    prof_clim_cell_index = F(prof_x, prof_y, prof_z);
    
    % go through each z level in the profile array
    %%
    clear tmp*;
    
    % set the default climatology value to be fillVal (-9999)
    prof_clim_T = ones(num_profs, num_prof_depths).*fillVal;
    prof_clim_S = ones(num_profs, num_prof_depths).*fillVal;
    
    for k = 1:min(num_prof_depths, 102)
        T_clim_k = squeeze(WOA_2013_v2_clim.potential_T_monthly(:,k,:,:));
        S_clim_k = squeeze(WOA_2013_v2_clim.S_monthly(:,k,:,:));
        
        % get the T and S at each profile point at this depth level
        for m = 1:12
            T_clim_mk = squeeze(T_clim_k(m,:,:));
            S_clim_mk = squeeze(S_clim_k(m,:,:));
            
            profs_in_month = find(prof_month == m);
            prof_clim_cell_index_m = prof_clim_cell_index(profs_in_month);
            
            prof_clim_T_tmp = T_clim_mk(prof_clim_cell_index_m);
            prof_clim_S_tmp = S_clim_mk(prof_clim_cell_index_m);
            
            tmp_T = prof_clim_T(:,k);
            tmp_T(profs_in_month) = prof_clim_T_tmp;
            
            tmp_S = prof_clim_S(:,k);
            tmp_S(profs_in_month) = prof_clim_S_tmp;
            
            prof_clim_T(:,k) = tmp_T;
            prof_clim_S(:,k) = tmp_S;
        end % m
    end % k
    ['finished mapping clim T and S to profiles']
  
    
    %% Fill -9999 for NaNs
    prof_clim_S(find(isnan(prof_clim_S)))= fillVal;
    prof_clim_T(find(isnan(prof_clim_T)))= fillVal;
    
    MITprof.prof_Tclim = prof_clim_T;
    MITprof.prof_Sclim = prof_clim_S;
    
    %%    Write output
    
    if write_profile_to_netcdf
        fileOut=[output_dir '/' fDataOut{ilist}];
        fprintf('%s\n',fileOut);
        cd(output_dir);
        
        %%
        write_profile_structure_to_netcdf(MITprof, fileOut);
        
    else
        ['I am not saving anything, the rest is up to you']
    end % write
end % ilist
    


