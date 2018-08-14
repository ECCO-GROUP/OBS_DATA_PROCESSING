function update_sigmaTS_on_prepared_profiles(run_code);
% This script updates ECCO profile.nc files with new sigma and climatology 
% fields and sends the updated profile files to be written in netcdf format

% Filename; update_sigmaTS_on_prepared_profiles.m
%  ** former filename : 
% Date Created: 2014-06-?
% Last Modified: 2016-03-09

% notes:
% -- carry over the zero weights
% -- if a point has a cost greater than 'max_cost_cutoff' set the whole
%    profile to zero weight
%  -- if a profile has fewer than min_vaid_pts_cutoff, set the whole
%    profile to zero weight
%
%  2014-08-04 : interp sigma in Z not nearest Z level.
%  2014-08-05 : use a sigma field which is extrapolated to all dry points
%  2015-01-12 : add the option to adjust the uncertainties 
%               \sigma_{new} = \sigma_{old} * 1./sqrt(alpha)
%               where alpha = [ A(x) / max(A)]
%               weights therefore change by factor alpha
%  2015-02-03 : add a new ARGO file [ARGO_nov2014_2011_to_2014.nc]
%  2015-02-06 : remove 'split by year' (bad idea to put it here)
%               and changed fileIn so that it can read the 20 argo files
%  2015-04-28 : make it so we remember past runs.
%             : remove step_05a and replaced with write_profile_to_netcdf
%  2015-09-15 : added profile weight area factor gamma as a field in the
%               netcdf file
%  2016-03-09 : moved file list to run code (as is the new standard)
%  2016-08-16 : changed to function
%  2018-05-08 : removed questionable 'bogus' flag where if the weights were
%               zero it would not replace them.  this is a problem when
%               creating fields for the first time, since they are all
%               zero at that time.


set(0,'DefaultTextInterpreter','none');
%%
hostname = getenv('HOSTNAME')

if length(strfind(hostname,'skylla2'))
    rootdir = '/net/skylla/data17/IGF_output'
elseif length(strfind(hostname,'skylla'))
    rootdir = '/data17/IGF_output'
elseif length(strfind(hostname,'pfe'))
    rootdir = '/nobackup/ifenty';
elseif length(strfind(hostname,'bridge'))
    rootdir = '/nobackup/ifenty';
elseif length(strfind(hostname,'cnidarian'))
    rootdir = '/home/ifenty/data/';
end

%%---------------------------------------------
%% SET INPUT PARAMETERS
%%---------------------------------------------
fillVal=-9999;checkVal=-9000;

% output and figure directories
%output_dir  = [rootdir '/observations/insitu/ECCO_v4r2']

% profile input dir and filenames
%input_dir  = [rootdir '/observations/insitu/ECCO_v4r2/20140803']

% output and figure directories


clear apply_area_factor apply_TS_floor save_format;


%--- save format  ---
%
%   0 = netcdf
%   1 = matlab
        
%----  apply_area_factor  -----
%
%    flag specifying whether to add a factor 1/sqrt(alpha)
%    where alpha = area/max(area) of the grid cell area in which
%    this profile is found.

% --- apply_TS_floor ---
%
%    flag specifiying whether to add the above T + S floor to the sigma
%    0 = no, 1= yes
        
        
make_figs = 0;

switch run_code
    case 'argo_20180610_gamma'
        
        basedir = ['/home/ifenty/data/observations/insitu/ARGO/from_gael_June2018/combined_by_latest/']
        input_dir = [basedir 'step_05_remove_zero_T_S_weight_profs'];
        output_dir = [basedir 'step_06_sigma_gamma']
        
        file_suffix_in  = 'step_05.nc'
        file_suffix_out = 'final_gamma.nc'
        
        cd(input_dir);
        f = dir(['*nc']);
        for i = 1:length(f)
            fDataIn{i}=   f(i).name;
            fDataOut{i} =[f(i).name(1:end-length(file_suffix_in)) file_suffix_out];
        end
        
        ilists = 1:length(f);
        
        apply_area_factor = 1; % use the gamma or not.
        apply_TS_floor = 0;  %0 means that the sigma field has already been capped
        save_format = 0;  %0 = netcdf, 1=matlab.
        
    case 'argo_20180610'
        
        basedir = ['/home/ifenty/data/observations/insitu/ARGO/from_gael_June2018/combined_by_latest/']
        input_dir = [basedir 'step_03_geobins_02562'];
        output_dir = [basedir 'step_04_sigmas']
        
        file_suffix_in  = 'step_03.nc'
        file_suffix_out = 'step_04.nc'
        
        cd(input_dir);
        f = dir(['*nc']);
        for i = 1:length(f)
            fDataIn{i}=   f(i).name;
            fDataOut{i} =[f(i).name(1:end-length(file_suffix_in)) file_suffix_out];
        end
        
        ilists = 1:length(f);
        
        apply_area_factor = 0; % use the gamma or not.
        apply_TS_floor = 0;  %0 means that the sigma field has already been capped
        save_format = 0;  %0 = netcdf, 1=matlab.
        
    case 13
        %%
        input_dir = ['/ian4/ifenty/data/observations/insitu/NODC/NODC_20180508/all/step_03_geobins_02562']
        output_dir = ['/ian4/ifenty/data/observations/insitu/NODC/NODC_20180508/all/step_04_sigmas/']
        
        
        file_suffix_in  = 'step_03.nc'
        file_suffix_out = 'step_04.nc'
        
        
        cd(input_dir);
        f = dir(['*nc']);
        for i = 1:length(f)
            fDataIn{i}=   f(i).name;
            fDataOut{i} =[f(i).name(1:end-length(file_suffix_in)) file_suffix_out];
        end
        
        ilists = 1:length(f);
        
        apply_area_factor = 0; % use the gamma or not.
        apply_TS_floor = 0;  %0 means that the sigma field has already been capped
        save_format = 0;  %0 = netcdf, 1=matlab.
        
    case 14
        %%
        input_dir = ['/ian4/ifenty/data/observations/insitu/NODC/NODC_20180508/all/step_08b_remove_zero_T_S_weight_profs']
        output_dir = ['/ian4/ifenty/data/observations/insitu/NODC/NODC_20180508/all/step_09c_area_weighted_gamma/']
        
        
        file_suffix_in  = 'step_08.nc'
        file_suffix_out = 'final_gamma.nc'
        
        
        cd(input_dir);
        f = dir(['*nc']);
        for i = 1:length(f)
            fDataIn{i}=   f(i).name;
            fDataOut{i} =[f(i).name(1:end-length(file_suffix_in)) file_suffix_out];
        end
        
        ilists = 1:length(f);
        
        apply_area_factor = 1; % use the gamma or not.
        apply_TS_floor = 0;  %0 means that the sigma field has already been capped
        save_format = 0;  %0 = netcdf, 1=matlab.
        
end

%%
fig_dir    = [output_dir '/figures']
mkdir(output_dir);mkdir(fig_dir)


%%
% sigma dir and filenames
sigma_dir  = [rootdir '/uncertainties/insitu/insitu_T_S_sigma_ecco_v6_capped/llc90/']

fsigST={['/Salt_sigma_smoothed_method_02_masked_merged_capped_extrapolated.bin'],...
       ['/Theta_sigma_smoothed_method_02_masked_merged_capped_extrapolated.bin']};

% the floor for the sigma fields
T_floor = 0.05;
S_floor = 0.01;

% sanity check the profile data with these mins + maxes
minT = -3;
minS = 20;

maxT = 40;
maxS = 40;



%%---------------------------------------------
%%  BEGIN THE PROGRAM
%%---------------------------------------------

%%---------------------------------------------
%%  Init profile package
%%---------------------------------------------
MITprof_path;


%%---------------------------------------------
%%  Read in grid for sigma and climatology
%%---------------------------------------------

load_llc90_grid
make_llc90_z_map

hf = hFacC_90;%
hf(find(hf~=1))=nan;

%%---------------------------------------------
%%  Read in sigma fields.
%%---------------------------------------------
cd(sigma_dir)

% salt
sigma_S = readbin(fsigST{1},[llcN 13*llcN 50],1,'real*4',0,'ieee-be');
sigma_S_90_pf = patchface3D(llcN, llcN*13, 50, sigma_S , 2);

% theta
sigma_T = readbin(fsigST{2},[llcN 13*llcN 50],1,'real*4',0,'ieee-be');
sigma_T_90_pf = patchface3D(llcN, llcN*13, 50, sigma_T , 2);

%%---------------------------------------------
%% Prepare the nearest neighbor mapping
%%---------------------------------------------
cd(llc90_grid_dir)
if exist('F_llc90_WET_INS_SURF_XYZ_to_INDEX.mat','file')
     ['loading F']
     load('F_llc90_WET_INS_SURF_XYZ_to_INDEX');
     ['loaded F']
else
     make_F_llc90_WET_INS_SURF_XYZ_to_INDEX;
end



%%
% verify that our little trick works in 4 parts of the earth
deg2rad=pi/180;

for i = 1:4
    switch i
        case 1
            test_lat = 56;
            test_lon = -40;
        case 2
            test_lat = 60;
            test_lon = 10;
        case 3
            test_lat = -63.9420;
            test_lon = -2.0790;
        case 4
            test_lat = -69;
            test_lon = 60;
    end
    [test_x, test_y, test_z] = sph2cart(test_lon*deg2rad, test_lat*deg2rad, 1);
    test_ind = F_llc90_WET_INS_SURF_XYZ_to_INDEX(test_x, test_y, test_z);

    ['original (line 1) vs closest (line 2) x,y,z']
    [X_90(test_ind) Y_90(test_ind) Z_90(test_ind)]
    [test_x test_y test_z]
    
    ['original (line 1) vs closest (line 2) lat lon']
    [test_lat test_lon]
    [lat_90(test_ind) lon_90(test_ind)]
end


%%---------------------------------------------
%% Read and process the profile files
%%---------------------------------------------

for ilist=ilists
    %%
    cd(input_dir)
    
    clear fileOut MITprof
    if exist(fDataIn{ilist})
        MITprof = MITprof_read(fDataIn{ilist});
        
        num_profs = length(MITprof.prof_lat)
        num_prof_depths = length(MITprof.prof_depth);
        
        orig_profTweight = MITprof.prof_Tweight;
        
        if isfield(MITprof, 'prof_S')
            orig_profSweight = MITprof.prof_Sweight;
        end
        ['finished loading MITprof']
        %%
        
        % interp sigmas to the new vertical levels
        %sigma_S_MITprof_z = interp_3D_to_arbitrary_z_levels(sigma_S, ....
        %    z_top_90, z_cen_90, z_bot_90, MITprof.prof_depth);
        
         
        sigma_S_MITprof_z = interp_3D_to_arbitrary_z_levels(sigma_S, ....
            z_cen_90, MITprof.prof_depth, 0);
        %%
        %sigma_T_MITprof_z = interp_3D_to_arbitrary_z_levels(sigma_T, ....
        %    z_top_90, z_cen_90, z_bot_90, MITprof.prof_depth);
        
        sigma_T_MITprof_z = interp_3D_to_arbitrary_z_levels(sigma_T, ....
            z_cen_90, MITprof.prof_depth, 0);
        ['interpolated sigma to new levels']
        
        
        %%
        if make_figs
            close all;
            figure(1);clf;hold on;
            
            x=25;y=150;
            scatter(squeeze(sigma_S_MITprof_z(x,y,:)), -MITprof.prof_depth,'ro');
            scatter(squeeze(sigma_S(x,y,:)), -z_cen_90,'b.');
            grid;
        end
        
        %%
        
        ['mapping profiles to llc grid']
        [prof_x, prof_y, prof_z] = ...
            sph2cart(MITprof.prof_lon*deg2rad, MITprof.prof_lat*deg2rad, 1);
        
        % map a llc90 grid index to each profile.
        prof_llc90_cell_index = F_llc90_WET_INS_SURF_XYZ_to_INDEX(prof_x, prof_y, prof_z);
        
        %%
        sigma_S_MITprof_z_flat = reshape(sigma_S_MITprof_z,[90*1170, num_prof_depths]);
        tmp_sigma_S = sigma_S_MITprof_z_flat(prof_llc90_cell_index,:);
        
        sigma_T_MITprof_z_flat = reshape(sigma_T_MITprof_z,[90*1170, num_prof_depths]);
        tmp_sigma_T = sigma_T_MITprof_z_flat(prof_llc90_cell_index,:);
        
        %%
        
        % check whether to apply a floor to these sigmas -- they may have
        % been applied so don't double count
        if apply_TS_floor
            [' applying a new floor']
            tmp_sigma_T = tmp_sigma_T + T_floor;
            tmp_sigma_S = tmp_sigma_S + S_floor;
        else
            ['not applying a new floor']
        end
        
        % make the weights from the sigmas
        tmp_weight_S = 1./tmp_sigma_S.^2;
        tmp_weight_T = 1./tmp_sigma_T.^2;
        
        [' min sigma T, min sigma S']
        [mynanmin(tmp_sigma_T(:)) mynanmin(tmp_sigma_S(:))]
        [' max weight T, max weight S']
        [mynanmax(tmp_weight_T(:))  mynanmax(tmp_weight_S(:))]
        
        % Set the weights
        MITprof.prof_Tweight = tmp_weight_T;
        
        if isfield(MITprof,'prof_S')
            MITprof.prof_Sweight = tmp_weight_S;
        end
        
        
        %%
        % Set the field that notes whatever floor has been applied to the sigmas;
        MITprof.prof_Terr = MITprof.prof_Terr.*0 + T_floor;
        if isfield(MITprof,'prof_S')
            MITprof.prof_Serr = MITprof.prof_Serr.*0 + S_floor;
        end
        
        %%
        % ensure that bad data are given zero weights.
        clear bogus_ins;
        
        ['processing S']
        if isfield(MITprof,'prof_S')
            clear bogus_ins
            bogus_ins{1} = find(MITprof.prof_Sflag > 0);
            bogus_ins{2} = find(MITprof.prof_S < minS);
            bogus_ins{3} = find(MITprof.prof_S > maxS);
            % an identically zero salinity is probably bogus
            bogus_ins{4} = find(MITprof.prof_S == 0);

            bogus_s = [];
            for i = 1:4
                
                switch i
                    case 1
                        ['prof_Sflag > 0  '   num2str(length(bogus_ins{i}))]
                    case 2
                        ['prof S < min S  '  num2str(length(bogus_ins{i}))]
                    case 3
                        ['prof S > max S  '  num2str(length(bogus_ins{i}))]
                    case 4
                        ['prof S == 0  ' num2str(length(bogus_ins{i}))]
                end
                
                bogus_s = union(bogus_s, bogus_ins{i});
            end
            MITprof.prof_Sweight(bogus_s) = 0;
        end
        %%
        clear bogus_ins;
        bogus_ins{1} = find(MITprof.prof_Tflag > 0);
        bogus_ins{2} = find(MITprof.prof_T < minT);
        bogus_ins{3} = find(MITprof.prof_T > maxT);
        
        % an identically zero temperature is probably bogus.
        bogus_ins{4} = find(MITprof.prof_T == 0);
        
        bogus_t = [];
        for i = 1:4
            switch i
                case 1
                    ['prof_Tflag > 0  '   num2str(length(bogus_ins{i}))]
                case 2
                    ['prof T < min T  '  num2str(length(bogus_ins{i}))]
                case 3
                    ['prof T > max T  '  num2str(length(bogus_ins{i}))]
                case 4
                    ['prof T == 0  ' num2str(length(bogus_ins{i}))]
            end
            bogus_t = union(bogus_t, bogus_ins{i});
        end
        
        MITprof.prof_Tweight(bogus_t) = 0;
        %%
        
        sigma_T_orig = sqrt(1./orig_profTweight);
        sigma_T_new = sqrt(1./MITprof.prof_Tweight);
        ins = find(MITprof.prof_Tweight == 0);
        sigma_T_new(ins) = NaN;
        
        ins = find(orig_profTweight == 0);
        sigma_T_orig(ins) = NaN;
        
        ['processing S']
        if isfield(MITprof,'prof_S')
            sigma_S_orig = sqrt(1./orig_profSweight);
            sigma_S_new = sqrt(1./MITprof.prof_Sweight);
            
            ins = find(orig_profSweight == 0);
            sigma_S_orig(ins) = NaN;
            ins = find(MITprof.prof_Sweight == 0);
            sigma_S_new(ins) = NaN;
        end
        
        %%
        if apply_area_factor == 1;
            
            ['applying area factor -- ']
            % find global maximum grid cell area
            max_area = max(RAC_90(:));
            % find ratio of each grid cell area to the maximum
            alpha = RAC_90_pf./max_area;
            
            pp = MITprof.prof_point;
            
            alpha_pp = alpha(pp);
            
            MITprof.prof_area_gamma = alpha_pp;

            tmpT = MITprof.prof_Tweight;
            
            if isfield(MITprof,'prof_S')
                tmpS = MITprof.prof_Sweight;
            end
            
            for k = 1:length(MITprof.prof_depth)
                [k]
                tmpTk = tmpT(:,k);
                tmpTk = tmpTk.*alpha_pp;
                tmpT(:,k) = tmpTk;
                
                ['processing S']
                if isfield(MITprof,'prof_S')
                    tmpSk = tmpS(:,k);
                    tmpSk = tmpSk.*alpha_pp;
                    tmpS(:,k) = tmpSk;
                end
            end
            
            MITprof.prof_Tweight = tmpT;
            
            if isfield(MITprof,'prof_S')
                MITprof.prof_Sweight = tmpS;
            end
        end
        
        %%
        
        ['writing ' fDataOut{ilist}]
        fileOut=[output_dir '/'  fDataOut{ilist}]
        fprintf('%s\n',fileOut);
        switch save_format
            case 0
                write_profile_structure_to_netcdf(MITprof, fileOut)
                %%
            case 1
                fileOut = [fileOut '.mat']
                save(fileOut,'MITprof')
        end
        
        
    else % file does not exist
        fDataIn{ilist}
        ['DOES NOT EXIST']
    end
end
   


% %%
% case 0
% %input_dir  = [rootdir '/observations/insitu/ECCO_v5r1/20150209']
% %output_dir  = [rootdir '/observations/insitu/ECCO_v5r1/20150209']
% %file_suffix_in  = '_PP_CLIM_llc270_20150209.nc'
% %file_suffix_out = '_SIG_PP_CLIM_llc270_20150209.nc'
% case 1
% input_dir  = ['/ian2/ifenty/data/model_output/ECCOv4/r2_iteration_11/profiles']
% output_dir  = ['/ian2/ifenty/data/model_output/ECCOv4/r2_iteration_11/profiles_orig_sigma']
% file_suffix_in = '_AWPP_llc90_20150206_model.nc';
% file_suffix_out = '_20150428_orig_sigma_model.nc';
% 
% apply_area_factor = 0;
% apply_TS_floor = 0
% save_format = 0;
% case 2
% input_dir  = ['/home/ifenty/data/observations/insitu/ECCO_v4r2/new/']
% output_dir = input_dir;
% file_suffix_in  = '_llc270_01_20150616.nc';
% file_suffix_out = '_llc270_02_20150616.nc';
% 
% apply_area_factor = 0;
% apply_TS_floor = 0
% save_format = 0;
% case 3
% input_dir  = ['/home/ifenty/data/observations/insitu/ECCO_v5r1/20150601_llc270/']
% output_dir = input_dir;
% file_suffix_in  = '_llc270_20150601.nc';
% file_suffix_out = '_llc270_20150616.nc';
% 
% apply_area_factor = 0;
% apply_TS_floor = 0
% save_format = 0;
% 
% ilists = [1:19 21:33]
% case 4
% input_dir  = ['/home/ifenty/data/observations/insitu/ECCO_v5r1/20150616_llc270/']
% output_dir = input_dir;
% file_suffix_in  = '_llc270_05_20150616.nc';
% file_suffix_out = '_llc270_20150616.nc';
% 
% apply_area_factor = 0;
% apply_TS_floor = 0
% save_format = 0;
% 
% ilists = [1:19 21:33]
% 
% case 5
% input_dir  = ['/home/ifenty/data/observations/insitu/ECCO_v4r2/llc90_20150915']
% output_dir  = ['/home/ifenty/data/observations/insitu/ECCO_v4r2/llc90_20150915']
% 
% file_suffix_in  = '_llc90_01_20150915.nc';
% file_suffix_out = '_GAMMA_20150915.nc';
% 
% apply_area_factor = 1;
% apply_TS_floor = 0
% save_format = 0;
% 
% ilists = [33]
% case 6
% input_dir  =['/ian2/ifenty/data/model_output/ECCOv4/r2_iteration_38_fenty_biharmonic/']
% output_dir  = input_dir;
% 
% file_suffix_in  = '_AWPP_llc90_20150206_model.nc';
% file_suffix_out = '_GAMMA_20150602_model.nc';
% 
% apply_area_factor = 1;
% apply_TS_floor = 0
% save_format = 0;
% 
% ilists = [1:19 21:33]
% case 7
% input_dir  =['/ian2/ifenty/data/model_output/ECCOv4/r2_iteration_38_no_biharmonic/']
% output_dir  = input_dir;
% 
% file_suffix_in  = '_AWPP_llc90_20150206_model.nc';
% file_suffix_out = '_GAMMA_20150206_model.nc';
% 
% apply_area_factor = 1;
% apply_TS_floor = 0
% save_format = 0;
% 
% ilists = [1:19 21:32]
% 
% case 8
% input_dir = '/ian4/ifenty/data/observations/insitu/ECCO_v4r2/llc90_20160308'
% 
% file_suffix_in  = '_feb2016_llc90_step_05_20160308.nc';
% file_suffix_out = '_feb2016_llc90_step_06_20160308.nc';
% 
% output_dir = input_dir;
% 
% apply_area_factor = 1; % use the gamma or not.
% apply_TS_floor = 0;  %0 means that the sigma field has already been capped
% save_format = 0;  %0 = netcdf, 1=matlab.
% 
% yrs = 1995:2016
% for i = 1:length(yrs)
% fDataBase{i} = ['argo_' num2str(yrs(i))];
% end
% 
% clear fDataIn* fDataOut;
% for i = 1:length(fDataBase)
% fDataIn{i}  = [fDataBase{i} file_suffix_in]
% fDataOut{i} = [fDataBase{i} file_suffix_out]
% end
% 
% ilists = 1:length(fDataBase)
% 
% case 9  %% llc270   2016-03-11
% ['fix sigmas on the new llc270 argo profiles']
% input_dir = ['/home/ifenty/data/observations/insitu/ECCO_v5r1/llc270_20160308/']
% output_dir = input_dir;
% 
% file_suffix_in  = '_feb2016_llc270_step_00_20160308.nc'
% file_suffix_out = '_feb2016_llc270_step_01_20160308.nc'
% 
% apply_area_factor = 0; % use the gamma or not.
% apply_TS_floor = 0;  %0 means that the sigma field has already been capped
% save_format = 0;  %0 = netcdf, 1=matlab.
% 
% yrs = 1995:2016
% for i = 1:length(yrs)
% fDataBase{i} = ['argo_' num2str(yrs(i))];
% end
% 
% clear fDataIn* fDataOut;
% for i = 1:length(fDataBase)
% fDataIn{i}  = [fDataBase{i} file_suffix_in]
% fDataOut{i} = [fDataBase{i} file_suffix_out]
% end
% 
% ilists = 1:length(fDataBase)
% 
% case 10
% input_dir = '/ian2/ifenty/data/model_output/ECCOv4/r1_iteration_11_updated_20160414/profiles_no_gamma_factor'
% output_dir = '/ian2/ifenty/data/model_output/ECCOv4/r1_iteration_11_updated_20160414/profiles_gamma_factor'
% 
% file_suffix_in  = '.nc';
% file_suffix_out = '_add_prof_gamma_factor.nc';
% 
% apply_area_factor = 1; % use the gamma or not.
% apply_TS_floor = 0;  %0 means that the sigma field has already been capped
% save_format = 0;  %0 = netcdf, 1=matlab.
% 
% cd(input_dir)
% ncf = dir('ices*nc')
% 
% clear fDataIn* fDataOut;
% for i = 1:length(ncf)
% fDataIn{i}  = [ncf(i).name]
% ncf(i).name
% fDataOut{i} = [ncf(i).name(1:16) file_suffix_out]
% end
% 
% ilists = 1:length(ncf)
% %%
% case 11
% input_dir = ['/ian4/ifenty/data/observations/insitu/NODC/WOD13_AUG2016_CTD_GLD_APB_XBT_CSVFORMAT/WOD_PROC_STEP_05']
% output_dir = ['/ian4/ifenty/data/observations/insitu/NODC/WOD13_AUG2016_CTD_GLD_APB_XBT_CSVFORMAT/WOD_PROC_STEP_06']
% 
% file_suffix_in  = 'step_05.nc'
% file_suffix_out = 'step_06.nc'
% 
% cd(input_dir)
% yrs = 1992:2016
% 
% cd(input_dir);
% f = dir(['*nc']);
% for i = 1:length(f)
% fDataIn{i}=   f(i).name;
% fDataOut{i} =[f(i).name(1:end-length(file_suffix_in)) file_suffix_out];
% end
% 
% ilists = 1:length(f);
% 
% apply_area_factor = 1; % use the gamma or not.
% apply_TS_floor = 0;  %0 means that the sigma field has already been capped
% save_format = 0;  %0 = netcdf, 1=matlab.
% case 12
% %%
% input_dir = ['/ian4/ifenty/data/observations/insitu/NODC/NODC_20160801/TMPb']
% output_dir = ['/ian4/ifenty/data/observations/insitu/NODC/NODC_20160801/TMPb_sigma_v2']
% 
% 
% file_suffix_in  = 'step_08.nc'
% file_suffix_out = 'step_08_sigma_v3.nc'
% 
% cd(input_dir)
% 
% cd(input_dir);
% f = dir(['*nc']);
% for i = 1:length(f)
% fDataIn{i}=   f(i).name;
% fDataOut{i} =[f(i).name(1:end-length(file_suffix_in)) file_suffix_out];
% end
% 
% ilists = 1:length(f);
% 
% apply_area_factor = 1; % use the gamma or not.
% apply_TS_floor = 0;  %0 means that the sigma field has already been capped
% save_format = 0;  %0 = netcdf, 1=matlab.
