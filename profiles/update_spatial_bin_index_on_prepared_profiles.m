function [MITprof] = update_spatial_bin_index_on_prepared_profiles(run_code)


% This script updates each profile with a bin index that is specified from
% some file.

% Filename; update_spatial_bin_index_on_prepared_profiles.m

%  ** former filename :
% Date Created: 2014-05-13
% Last Modified: 2018-06-04
%

% notes:
%   2016-03-09  :  made filenames a part of the run_index.   
%   2016-08-18  :  converted to funciton
%   2018-06-04  :  added some debugging


set(0,'DefaultTextInterpreter','none');
%%
make_root_dir;

%%---------------------------------------------
%% SET INPUT PARAMETERS
%%---------------------------------------------
fillVal=-9999;checkVal=-9000;


%% bin_index directory and filename
% bin_dir = ''
% bin_file
% bin_llcN

%%---------------------------------------------
%%  BEGIN THE PROGRAM
%%---------------------------------------------

%%---------------------------------------------
%%  Init profile package
%%---------------------------------------------
MITprof_path;

%make_argo_profile_depths;

switch run_code
    %%%%%%%%%%%%%%%%%
    % 2018-05-10    %
    %%%%%%%%%%%%%%%%%
    
    % CTD 10242
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'WOD_13_20180509_all_add_10242_geodesic'
        input_dir = '/ian4/ifenty/data/observations/insitu/NODC/NODC_20180508/all/step_01_monthly_clim';
        output_dir = '/ian4/ifenty/data/observations/insitu/NODC/NODC_20180508/all/step_02_geobins_10242/'
        
        file_suffix_in  = 'step_01.nc';
        file_suffix_out = 'step_02.nc';
        
        % the bin data have to be projected to the model grid;
        bin_dir = ['~/data/grids/grid_llc90/sphere_point_distribution/']
        bin_file = ['llc090_sphere_point_n_10242_ids.bin']
        
        bin_llcN = 90;
        load_llc90_grid;
        prof_bin_name = 'prof_bin_id_a'
        yrs = 1990:2018
        
        cd(input_dir)
        f  = dir(['*nc']);
        for i = 1:length(f)
            fDataIn{i}=   f(i).name;
            fDataOut{i} =[f(i).name(1:end-length(file_suffix_in)) file_suffix_out];
        end
        
        % CTD 2562
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'WOD_13_20180509_all_add_2562_geodesic'
        input_dir = '/ian4/ifenty/data/observations/insitu/NODC/NODC_20180508/all/step_02_geobins_10242/'
        output_dir = '/ian4/ifenty/data/observations/insitu/NODC/NODC_20180508/all/step_03_geobins_02562/'
        
        file_suffix_in  = 'step_02.nc';
        file_suffix_out = 'step_03.nc';
        
        % the bin data have to be projected to the model grid;
        bin_dir = ['~/data/grids/grid_llc90/sphere_point_distribution/']
        bin_file = ['llc090_sphere_point_n_02562_ids.bin']
        
        bin_llcN = 90;
        load_llc90_grid;
        prof_bin_name = 'prof_bin_id_b'
        yrs = 1990:2018
        
        cd(input_dir)
        f  = dir(['*nc']);
        for i = 1:length(f)
            fDataIn{i}=   f(i).name;
            fDataOut{i} =[f(i).name(1:end-length(file_suffix_in)) file_suffix_out];
        end
        
    case 'argo_20180610_10242_geodesic'
        input_dir = '/home/ifenty/data/observations/insitu/ARGO/from_gael_June2018/combined_by_latest/step_00_update_tile_and_prof_points'
        output_dir = '/home/ifenty/data/observations/insitu/ARGO/from_gael_June2018/combined_by_latest/step_02_geobins_10242/'
        
        file_suffix_in  = 'step_00.nc';
        file_suffix_out = 'step_02.nc';
        
        % the bin data have to be projected to the model grid;
        bin_dir = ['~/data/grids/grid_llc90/sphere_point_distribution/']
        bin_file = ['llc090_sphere_point_n_10242_ids.bin']
        
        bin_llcN = 90;
        load_llc90_grid;
        prof_bin_name = 'prof_bin_id_a'
        yrs = 1990:2018
        
        cd(input_dir)
        f  = dir(['*nc']);
        for i = 1:length(f)
            fDataIn{i}=   f(i).name;
            fDataOut{i} =[f(i).name(1:end-length(file_suffix_in)) file_suffix_out];
        end
        
        % CTD 2562
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'argo_20180610_02562_geodesic'
        input_dir = '/home/ifenty/data/observations/insitu/ARGO/from_gael_June2018/combined_by_latest/step_02_geobins_10242/'
        output_dir= '/home/ifenty/data/observations/insitu/ARGO/from_gael_June2018/combined_by_latest/step_03_geobins_02562/'

        file_suffix_in  = 'step_02.nc';
        file_suffix_out = 'step_03.nc';
        
        % the bin data have to be projected to the model grid;
        bin_dir = ['~/data/grids/grid_llc90/sphere_point_distribution/']
        bin_file = ['llc090_sphere_point_n_02562_ids.bin']
        
        bin_llcN = 90;
        load_llc90_grid;
        prof_bin_name = 'prof_bin_id_b'
        yrs = 1990:2018
        
        cd(input_dir)
        f  = dir(['*nc']);
        for i = 1:length(f)
            fDataIn{i}=   f(i).name;
            fDataOut{i} =[f(i).name(1:end-length(file_suffix_in)) file_suffix_out];
        end
        
        
end

ilists = 1:length(fDataIn);

mkdir(output_dir)
%%---------------------------------------------
%%  Read in bin fioeld
%%---------------------------------------------
cd(bin_dir)

bin = readbin([bin_dir bin_file],[bin_llcN 13*bin_llcN 1 1],1,'real*4',0,'ieee-be');

stats(bin)

[run_code]

[file_suffix_in file_suffix_out]
[input_dir]
[output_dir]



%%---------------------------------------------
%% Prepare the nearest neighbor mapping
%%---------------------------------------------

switch bin_llcN
    case 90
        load_llc90_grid
        
        cd(llc90_grid_dir)
        if exist('F_llc90_ALL_INS_SURF_XYZ_to_INDEX.mat','file')
            ['loading F']
            load('F_llc90_ALL_INS_SURF_XYZ_to_INDEX');
            ['loaded F']
        else
            make_F_llc90_ALL_INS_SURF_XYZ_to_INDEX;
        end
        F = F_llc90_ALL_INS_SURF_XYZ_to_INDEX;
        X = X_90; Y = Y_90; Z = Z_90;
        lon_llc = lon_90;lat_llc = lat_90;
        
    case 270
        load_llc270_grid
        
        cd(llc270_grid_dir)
        if exist('F_llc270_ALL_INS_SURF_XYZ_to_INDEX.mat','file')
            ['loading F']
            load('F_llc270_ALL_INS_SURF_XYZ_to_INDEX');
            ['loaded F']
        else
            make_F_llc270_ALL_INS_SURF_XYZ_to_INDEX;
        end
        F = F_llc270_ALL_INS_SURF_XYZ_to_INDEX;
        X = X_270; Y = Y_270; Z = Z_270;
        lon_llc = lon_270;lat_llc = lat_270;
end

%%
% verify that our little trick works in 4 parts of the earth'
deg2rad = pi/180;
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
    test_ind = F(test_x, test_y, test_z);
    
    [X(test_ind) Y(test_ind) Z(test_ind) test_x test_y test_z]
    [lat_llc(test_ind) lon_llc(test_ind) test_lat test_lon]
end

%%---------------------------------------------
%% Read and process the profile files
%%---------------------------------------------

for ilist = ilists
    %%
    
    fileOut=[output_dir '/' fDataOut{ilist}]

    if ~exist(fileOut,'file')
        
        cd(input_dir)
        fDataIn{ilist}
        
        clear MITprof
        
        MITprof = MITprof_read(fDataIn{ilist});
        
        num_profs = length(MITprof.prof_lat)
        
        ['mapping profiles to bin llc grid']
        [prof_x, prof_y, prof_z] = ...
            sph2cart(MITprof.prof_lon*deg2rad, MITprof.prof_lat*deg2rad, 1);
        
        %%
        % map a grid index to each profile.
        prof_llcN_cell_index = F(prof_x, prof_y, prof_z);
        
        % go through each z level in the profile array
        if isfield(MITprof, 'prof_gci')
            'gci exists'
            if unique(MITprof.prof_gci - prof_llcN_cell_index) ~= 0
                [' prof_llcN_cell index does not equal prof_gci!!!']
                break;
            else
                'prof_gci and prof_llcN_cell_index are the same'
            end
        end
        
        
        MITprof = setfield(MITprof, prof_bin_name, bin(prof_llcN_cell_index));
        
        
        %  Write output
        
        fprintf('%s\n',fileOut);
        
        write_profile_structure_to_netcdf(MITprof, fileOut);
    else
        [fileOut ' exists, skipping....']
    end
end



% %%
%    case 'feb2016_add_02562_geodesic'
%             input_dir = '/ian4/ifenty/data/observations/insitu/ECCO_v4r2/llc90_20160308'
%             
%             file_suffix_in  = '_feb2016_llc90_step_04_20160308.nc';
%             file_suffix_out = '_feb2016_llc90_step_05_20160308.nc';
%             
%             output_dir = input_dir;
%             
%             % the bin data have to be projected to the model grid;
%             bin_dir = ['~/data/grids/grid_llc90/sphere_point_distribution/']
%             bin_file = ['llc090_sphere_point_n_02562_ids.bin']
%             bin_llcN = 90;
%             
%             load_llc90_grid
%             
%             prof_bin_name = 'prof_bin_id_b'
%             
%             yrs = 1995:2016
%             for i = 1:length(yrs)
%                 fDataBase{i} = ['argo_' num2str(yrs(i))];
%             end
%             
%             clear fDataIn* fDataOut;
%             for i = 1:length(fDataBase)
%                 fDataIn{i}  = [fDataBase{i} file_suffix_in]
%                 fDataOut{i} = [fDataBase{i} file_suffix_out]
%             end
%             
%             ilists = 1:length(fDataBase)
%             
%             
%        
%         
%         case 'feb2016_add_10242_geodesic'
%             input_dir = '/ian4/ifenty/data/observations/insitu/ECCO_v4r2/llc90_20160308'
% 
%             file_suffix_in  = '_feb2016_llc90_step_03_20160308.nc';
%             file_suffix_out = '_feb2016_llc90_step_04_20160308.nc';
%             
%             output_dir = input_dir;
%             
%             % the bin data have to be projected to the model grid;
%             bin_dir = ['~/data/grids/grid_llc90/sphere_point_distribution/']
%             bin_file = ['llc090_sphere_point_n_10242_ids.bin']
%             
%             bin_llcN = 90;
%             
%             load_llc90_grid
%             
%             prof_bin_name = 'prof_bin_id_a'
%             
%             yrs = 1995:2016
%             for i = 1:length(yrs)
%                 fDataBase{i} = ['argo_' num2str(yrs(i))];
%             end
%             
%             clear fDataIn* fDataOut;
%             for i = 1:length(fDataBase)
%                 fDataIn{i}  = [fDataBase{i} file_suffix_in]
%                 fDataOut{i} = [fDataBase{i} file_suffix_out]
%             end
%             
%             ilists = 1:length(fDataBase)
%             
%         case 'feb2016_add_02562_geodesic'
%             input_dir = '/ian4/ifenty/data/observations/insitu/ECCO_v4r2/llc90_20160308'
%             
%             file_suffix_in  = '_feb2016_llc90_step_04_20160308.nc';
%             file_suffix_out = '_feb2016_llc90_step_05_20160308.nc';
%             
%             output_dir = input_dir;
%             
%             % the bin data have to be projected to the model grid;
%             bin_dir = ['~/data/grids/grid_llc90/sphere_point_distribution/']
%             bin_file = ['llc090_sphere_point_n_02562_ids.bin']
%             bin_llcN = 90;
%             
%             load_llc90_grid
%             
%             prof_bin_name = 'prof_bin_id_b'
%             
%             yrs = 1995:2016
%             for i = 1:length(yrs)
%                 fDataBase{i} = ['argo_' num2str(yrs(i))];
%             end
%             
%             clear fDataIn* fDataOut;
%             for i = 1:length(fDataBase)
%                 fDataIn{i}  = [fDataBase{i} file_suffix_in]
%                 fDataOut{i} = [fDataBase{i} file_suffix_out]
%             end
%             
%             ilists = 1:length(fDataBase)
%             