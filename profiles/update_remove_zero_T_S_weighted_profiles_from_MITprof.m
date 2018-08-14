function update_remove_zero_T_S_weighted_profiles_from_MITprof(run_code)

% This script removes profiles that have all zero T and S weights
% from MITprof structures

% Filename; update_remove_zero_T_S_weighted_profiles_from_MITprof

%  ** former filename :
% Date Created: 2018-06-07
% Last Modified: 2018-06-07

% notes:
%    we also have the ability to truncate the number of depths.  we'll
%    cut at 6300 m for the llc90.

set(0,'DefaultTextInterpreter','none');
%%
make_root_dir;

%%---------------------------------------------
%% SET INPUT PARAMETERS
%%---------------------------------------------
fillVal=-9999;checkVal=-9000;

%%---------------------------------------------
%%  BEGIN THE PROGRAM
%%---------------------------------------------

debug_code=0;

switch run_code
    case 'argo_20180610'
        basedir = ['/home/ifenty/data/observations/insitu/ARGO/from_gael_June2018/combined_by_latest/']
        input_dir = [basedir 'step_04_sigmas'];
        output_dir = [basedir 'step_05_remove_zero_T_S_weight_profs']
        
        file_suffix_in  = 'step_04.nc'
        file_suffix_out = 'step_05.nc';
        
        cd(input_dir);
        f = dir(['*nc']);
        for i = 1:length(f)
            fDataIn{i}=   f(i).name;
            fDataOut{i} =[f(i).name(1:end-length(file_suffix_in)) ...
                file_suffix_out];
        end
     
        
    case 'NODC_20180508'
        rootdir = ['/home/ifenty/data/observations/insitu/NODC/NODC_20180508/']
        input_dir  = [rootdir '/all/step_07_zero_weights_high_cost_vs_clim/']
        output_dir = [rootdir '/all/step_08_remove_zero_T_S_weight_profs/']
        
        file_suffix_in  = 'step_07.nc'
        file_suffix_out = 'step_08.nc';
        
        cd(input_dir);
        f = dir(['*nc']);
        for i = 1:length(f)
            fDataIn{i}=   f(i).name;
            fDataOut{i} =[f(i).name(1:end-length(file_suffix_in)) ...
                file_suffix_out];
        end
        

        
    case 'NODC_20180508_high_cost_vs_clim_do_not_exclude_high_lats'

        rootdir = ['/home/ifenty/data/observations/insitu/NODC/NODC_20180508/']
        input_dir  = [rootdir '/all/step_07b_zero_weights_high_cost_vs_clim/']
        output_dir = [rootdir '/all/step_08b_remove_zero_T_S_weight_profs/']
        
        file_suffix_in  = 'step_07.nc'
        file_suffix_out = 'step_08.nc';
        
        cd(input_dir);
        f = dir(['*nc']);
        for i = 1:length(f)
            fDataIn{i}=   f(i).name;
            fDataOut{i} =[f(i).name(1:end-length(file_suffix_in)) ...
                file_suffix_out];
        end

        
end
ilists = 1:length(f);

%%
mkdir(output_dir)
%%
for ilist = ilists
    %%
    cd(input_dir)
    fDataIn{ilist}
    clear fileOut MITprof
    
    fileOut=[output_dir '/' fDataOut{ilist}]
    
    if 1==1 %~exist(fileOut,'file')
        
        if exist(fDataIn{ilist})
            clear MITprof;
            %%
            clear *_profs;
            
            MITprof = MITprof_read(fDataIn{ilist});
            MITprof_orig = MITprof;
            
            [fDataIn{ilist}]
            
            total_Tweight = sum(MITprof.prof_Tweight,2);
            total_Sweight = sum(MITprof.prof_Sweight,2);
            
            total_TS_weight = total_Tweight + total_Sweight;
            
            all_profs = 1:length(MITprof.prof_YYYYMMDD);
            
            good_profs = find(total_TS_weight > 0);
            bad_profs  = find(total_TS_weight == 0);
            
            nan_profs = find(isnan(total_TS_weight));
            
            ['# good_profs, bad_profs, nan_profs']
            [length(good_profs) length(bad_profs) length(nan_profs)]
            %%
            if length(nan_profs) > 0
                ['you have nans in your weights, this should never happen']
                stop
            end
            
            if length(good_profs) > 0    
                ['there is at least one good prof to write']
                
                if length(bad_profs) > 0
                    ['at least one bad prof will be removed']
                    MITprof_new = ...
                        extract_profile_subset_from_MITprof(...
                        MITprof, good_profs, []);
                else
                    ['all profs are good, not doing anything special']
                    MITprof_new = MITprof;
                end
                
                %%
                %  Write output
                ['writing output']
                fileOut=[output_dir '/' fDataOut{ilist}];
                fprintf('%s\n',fileOut);
                
                write_profile_structure_to_netcdf(MITprof_new,fileOut);
            else
                ['THERE ARE NO GOOD PROFILES TO SAVE']
            end % length good_profs 
        else
            [fDataIn{ilist} ' does not exist']
        end
    else
        [fDataOut{ilist} ' already exists, skipping!']
    end
end % ilist