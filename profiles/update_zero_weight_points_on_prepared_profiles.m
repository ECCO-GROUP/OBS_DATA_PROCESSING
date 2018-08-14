function update_zero_weight_points_on_prepared_profiles(run_code)

% This script zeros out  profile profTweight and profSweight on 
% points that match some criteria

% Filename; update_zero_weight_points_on_prepared_profiles.m

%  ** former filename :
% Date Created: 2014-05-22
% Last Modified: 2018-06-04

% notes:

%   2015-07-07 : fix logic on finding bad data. cleaned up documentation
%   2016-03-09 : move file list to run_code as per new standard.
%   2018-03-27 : check for bad dates and times as well.
%   2018-06-04 : fixed comments

set(0,'DefaultTextInterpreter','none');
%
make_root_dir;

%%---------------------------------------------
%% SET INPUT PARAMETERS
%%---------------------------------------------
fillVal=-9999;checkVal=-9000;


%%---------------------------------------------
%%  BEGIN THE PROGRAM
%%---------------------------------------------

% zero critera code
% 0;  where there is no data
% 1;  where there is no climatology value
% 2;  where T and S outside range
% 3:  where dates are invalid.
% 4:  lat-lons outside range or at 0,0
% 5:  high cost vs. climatology

clear zero_criteria_code;

debug_code=0;
make_figs=0;
                       

% for zero criteria code 5 (high cost vs. clim) cost_vs_clim_filter_method
%
% cost_vs_clim_filter_method :  0 use cost_threshold percentage
%                               1 use cost_threshold as number
%
% exclude_high_latitude_profiles_from_test_2 : test 2 determines whether
%                   even a single profile point has a cost vs. clim that exceeds
%                   some threshold. since the climatologies at high
%                   latitutdes are not reliable, it is not reasonable to
%                   exclude high cost vs. clim profiles there (here defined
%                   to be 60N-90S and 60N-90N)
%
% bad_profs_to_plot  :  the number of suspect profiles to plot
%
% plot_individual_bad_profiles : 0/1 whether to plot individual profiles
%
% plot_map_bad_profiles : 0/1 whether to make a plot of bad profs locations
                        

switch run_code
   
    case 'argo_20180610'
        
        basedir = ['/home/ifenty/data/observations/insitu/ARGO/from_gael_June2018/combined_by_latest/']
        input_dir = [basedir 'step_04_sigmas']
        output_dir = [basedir '/all/step_05_zero_weights/']
        
        fig_dir = [output_dir '/figures'];
        
        file_suffix_in  = 'step_04.nc';
        file_suffix_out = 'step_05.nc';
        
        zero_criteria_code = [0 2 3 4]
        
        prof_Tmin = -1.96;
        prof_Tmax = 60;
        prof_Smin = 0;
        prof_Smax = 60;
        
        
        cd(input_dir);
        f = dir(['*nc']);
        for i = 1:length(f)
            fDataIn{i}=   f(i).name;
            fDataOut{i} =[f(i).name(1:end-length(file_suffix_in)) ...
                file_suffix_out];
        end
        
        ilists = 1:length(f)
        
        
    case 'NODC_20180508'
        
        rootdir = ['/home/ifenty/data/observations/insitu/NODC/NODC_20180508/']
        input_dir  = [rootdir '/all/step_04_sigmas/']
        output_dir = [rootdir '/all/step_05_zero_weights/']
        fig_dir = [output_dir '/figures'];
        
        file_suffix_in  = 'step_04.nc';
        file_suffix_out = 'step_05.nc';
        
        zero_criteria_code = [0 1 2 3 4]
        
        prof_Tmin = -1.96;
        prof_Tmax = 60;
        prof_Smin = 0;
        prof_Smax = 60;
        
        
        cd(input_dir);
        f = dir(['*nc']);
        for i = 1:length(f)
            fDataIn{i}=   f(i).name;
            fDataOut{i} =[f(i).name(1:end-length(file_suffix_in)) ...
                file_suffix_out];
        end
        
        ilists = 1:length(f)
        
     case 'NODC_20180508_high_cost_vs_clim'
         %%
        rootdir = ['/home/ifenty/data/observations/insitu/NODC/NODC_20180508/']
        
        input_dir  = [rootdir '/all/step_06_potential_temperature']
        output_dir = [rootdir '/all/step_07_zero_weights_high_cost_vs_clim/']
        fig_dir = [output_dir '/figures']
        
        file_suffix_in  = 'step_06.nc';
        file_suffix_out = 'step_07.nc';
        
        zero_criteria_code = [5]
        
        num_bad_profs_to_plot=200;
                           
        %cost_vs_clim_filter_method = 0 % 0 use cost_threshold  percentage
        cost_vs_clim_filter_method = 1 % 1 use cost_threshold as number
        exclude_high_latitude_profiles_from_test_2 = 1;
        
        % cost_threshold = 0.995;
        cost_threshold = 16;
        
        plot_individual_bad_profiles = 0
        plot_map_bad_profiles = 1

        cd(input_dir);
        f = dir(['*nc']);
        for i = 1:length(f)
            fDataIn{i}=   f(i).name;
            fDataOut{i} =[f(i).name(1:end-length(file_suffix_in)) ...
                file_suffix_out];
        end
        
        ilists = 1:length(f);
        make_figs = 1;
        save_figs = 1;
        
    case 'NODC_20180508_high_cost_vs_clim_do_not_exclude_high_lats'
         %%
        rootdir = ['/home/ifenty/data/observations/insitu/NODC/NODC_20180508/']
        
        input_dir  = [rootdir '/all/step_06_potential_temperature']
        output_dir = [rootdir '/all/step_07b_zero_weights_high_cost_vs_clim/']
        fig_dir = [output_dir '/figures']
        
        file_suffix_in  = 'step_06.nc';
        file_suffix_out = 'step_07.nc';
        
        zero_criteria_code = [5]
        
        num_bad_profs_to_plot=200;
                           
        %cost_vs_clim_filter_method = 0 % 0 use cost_threshold  percentage
        cost_vs_clim_filter_method = 1 % 1 use cost_threshold as number
        exclude_high_latitude_profiles_from_test_2 = 0;
        
        % cost_threshold = 0.995;
        cost_threshold = 16;
        
        plot_individual_bad_profiles = 0
        plot_map_bad_profiles = 1

        cd(input_dir);
        f = dir(['*nc']);
        for i = 1:length(f)
            fDataIn{i}=   f(i).name;
            fDataOut{i} =[f(i).name(1:end-length(file_suffix_in)) ...
                file_suffix_out];
        end
        
        ilists = 1:length(f);
        make_figs = 1;
        save_figs = 1;
end
%%
mkdir(fig_dir);
mkdir(output_dir)

for ilist = ilists
    %%
    cd(input_dir)
    
    fDataIn{ilist};
    
    clear fileOut MITprof
    
    fileOut=[output_dir '/' fDataOut{ilist}]
    
    if 1==1~exist(fileOut,'file')
        %%
        if exist(fDataIn{ilist})
            clear MITprof;
            cd(input_dir)
            MITprof = MITprof_read(fDataIn{ilist});
            MITprof_orig = MITprof;
            [fDataIn{ilist}]
            
            %%
            % loop through the zero criteria codes
            for i = zero_criteria_code
                
                %% Make a copy of MITprof.
                MITprof_new = MITprof;
                num_profs = length(MITprof.prof_YYYYMMDD);
                %%
                switch i
                    
                    % 0;  zero weights where there are no observations
                    case 0
                        
                        clear ins*;
                        
                        ['check to see if there is no data']
                        if debug_code
                            ['# nonzero T and S weights']
                            [length(find(MITprof_new.prof_Tweight > 0)) ...
                            length(find(MITprof_new.prof_Tweight > 0))]
                        end
                        
                        %% If data is nan or -9999, then
                        %% zero the weight and change data to -9999
                        ['T']
                        
                        % is T nan or fillVal?
                        ins1 = find(isnan(MITprof.prof_T));
                        ins2 = find(MITprof.prof_T <= checkVal);
                        ins3 = union(ins1,ins2);
                        
                        ['# prof T is nan, is -9999']
                        [length(ins1), length(ins2)]
                        
                        % if so, put fill val there (just in case we had
                        % NaN)
                        MITprof_new.prof_T(ins3) = fillVal;
                        
                        ins4 = find(isnan(MITprof.prof_Tweight));
                        ins5 = find(MITprof.prof_Tweight <= 0);
                        ins6 = union(ins4, ins5);
                        
                        ['# prof T weight is nan, is zero']
                        [length(ins4), length(ins5)]
                        
                        % no data or zero weight or nan weight
                        ins7 = union(ins3,ins6);
                        
                        % if weight is nan or < 0, then zero weight
                        MITprof_new.prof_Tweight(ins7) = 0;
                                                
                        % Now salinity
                        if isfield(MITprof,'prof_S')
                            clear ins*
                            ['S']
                            
                            ins1 = find(isnan(MITprof.prof_S));
                            ins2 = find(MITprof.prof_S <= checkVal);
                            ins3 = union(ins1,ins2);
                            ['# prof S is nan, is -9999']
                            [length(ins2), length(ins3)]
                            
                            MITprof_new.prof_S(ins3) = fillVal;
                            
                            ins4 = find(isnan(MITprof.prof_Sweight));
                            ins5 = find(MITprof.prof_Sweight <= 0);
                            ins6 = union(ins4, ins5);
                            
                            ['# prof S weight is nan, is zero']
                            [length(ins4), length(ins5)]
                            
                            % no data or zero weight or nan weight
                            ins7 = union(ins3, ins6);
                            
                            % if weight is nan or < 0, zero weight
                            MITprof_new.prof_Sweight(ins7) = 0;
                            
                            ins8 = find(MITprof_new.prof_Sweight ==0);

                        end

                        MITprof = MITprof_new;
                        clear MITprof_new;
                        
                     
                    case 1 % zero weights where no climatology value
                        
                        ['check to see if there is no clim value']
                        if debug_code
                            ['# nonzero T and S weights']
                            [length(find(MITprof_new.prof_Tweight > 0)) ...
                            length(find(MITprof_new.prof_Tweight > 0))]
                        end
                        
                        ins1 = find(MITprof.prof_Tclim <= checkVal);
                        ins2 = find(isnan(MITprof.prof_Tclim));
                        ins3 = union(ins1,ins2);
                        length(ins3);
                        
                        MITprof_new.prof_Tweight(ins3) = 0;
                        
                        ['# no T climatology value']
                        [length(ins3)]
                        ins3=[];
                        
                        if isfield(MITprof,'prof_S')
                            ['S clim']
                            ins1 = find(MITprof.prof_Sclim <= checkVal);
                            ins2 = find(isnan(MITprof.prof_Sclim));
                            ins3 = union(ins1,ins2);
                            length(ins3);
                            MITprof_new.prof_Sweight(ins3) = 0;
                            
                            ['# no Sclimatology value']
                            [length(ins3)]
                            
                            zero_S_ins = find(MITprof.prof_Sclim == 0)
                            if length(zero_S_ins) > 0
                                ['I found zeros in your climatology S']
                                ['those should be -9999s most likely']
                                ['which means your T values also probably']
                                ['have zeros where data are missing']
                                stop
                            end
                        end
                        
                        MITprof = MITprof_new;
                        clear MITprof_new;
                        
                    case 2 % prof T, S exceeding min max
                        ['check to see if T and S are out of bounds']
                        if debug_code
                            ['# nonzero T and S weights']
                            [length(find(MITprof_new.prof_Tweight > 0)) ...
                            length(find(MITprof_new.prof_Tweight > 0))]
                        end
                        
                        MITprof_new.prof_T( ...
                            find(MITprof.prof_T == fillVal)) = NaN;
                        
                        ins1 = find(MITprof.prof_T < prof_Tmin & ...
                            MITprof.prof_T ~= fillVal);
                        MITprof_new.prof_Tweight(ins1) = 0;
                        
                        ins2 = find(MITprof.prof_T > prof_Tmax & ...
                            MITprof.prof_T ~= fillVal);
                        
                        MITprof_new.prof_Tweight(ins2) = 0;
                        
                        % We're just changing the weight, not the profile
                        % value itself.
                        
                        ['# too cold, # too warm']
                        [length(ins1) length(ins2)]
                        
                        if isfield(MITprof,'prof_S')
                            
                            ins3 = find(MITprof.prof_S < prof_Smin & ...
                                MITprof.prof_S ~= fillVal);
                            MITprof_new.prof_Sweight(ins3) = 0;
                            
                            ins4 = find(MITprof.prof_S > prof_Smax & ...
                                MITprof.prof_S ~= fillVal);
                            MITprof_new.prof_Sweight(ins4) = 0;
                            
                            ['# too fresh, # too salty']
                            [length(ins3) length(ins4)]
                        end
 
                        MITprof = MITprof_new;
                        clear MITprof_new;
                        
                    case 3 % illegal dates

                        clear bad_profs good_profs all_profs all_ins;
                        
                        ['check to see if there are Illegal dates']
                        if debug_code
                            ['# nonzero T and S weights']
                            [length(find(MITprof_new.prof_Tweight > 0)) ...
                            length(find(MITprof_new.prof_Tweight > 0))]
                        end
                        
                        all_years = floor(MITprof.prof_YYYYMMDD/1e4);
                        years = unique(all_years)
                        
                        clear y m d;
                        for i = 1:num_profs
                            tmp = num2str(MITprof.prof_YYYYMMDD(i));
                            y(i)  = str2num(tmp(1:4));
                            m(i)  = str2num(tmp(5:6));
                            d(i)  = str2num(tmp(7:8));
                        end
                        
                        %%
                        % get current date
                        s = date;
                        
                        todays_year = str2num(s(end-3:end));
                        % bad years are pre 1900 and after today's year
                        bad_years_a = find(y < 1900);
                        bad_years_b = find(y > todays_year);
                        bad_mons_a  = find(m < 1);
                        bad_mons_b  = find(m > 12);
                        bad_days_a  = find(d < 1);
                        bad_days_b  = find(d > 31);
                        
                        bad_years = union(bad_years_a, bad_years_b);
                        bad_mons = union(bad_mons_a, bad_mons_b);
                        bad_days = union(bad_days_a, bad_days_b);
                        
                        bad_times_A = find(MITprof.prof_HHMMSS < 0);
                        bad_times_B = find(MITprof.prof_HHMMSS > 240000);
                        
                        bad_profs = union(union(union(union(bad_years, ...
                            bad_mons), ...
                            bad_days), bad_times_A), bad_times_B);
                        
                        all_profs = 1:num_profs;
                        
                        good_profs = setdiff(all_profs, bad_profs);
                        
                        ['date check, all, good, bad']
                        [length(all_profs) length(good_profs) ...
                            length(bad_profs)]
                        
                        %%
                        MITprof_new.prof_Tweight(bad_profs,:) = 0;
                        MITprof_new.prof_Sweight(bad_profs,:) = 0;
        
                        MITprof = MITprof_new;
                        clear MITprof_new;
                        
                    case 4 % zero weights for lat-lon out of bounds or 0
                        ['check to see if lat-lon = 0 or out of bounds']
                        clear bad_profs good_profs all_profs all_ins;
                        
                        ['bad lats lons check']
                        lats = MITprof.prof_lat;
                        lons = MITprof.prof_lon;
                        
                        outside_lat_lim_ins = ...
                            find(lats <= -90 | lats >= 90);
                        outside_lon_lim_ins = ...
                            find(lons <= -180 | lons >= 180);
                        zero_lat_lon_ins = find(lats == 0 & lons == 0);
                        
                        bad_profs = union(union(outside_lat_lim_ins,...
                            outside_lon_lim_ins), zero_lat_lon_ins);
                       
                        all_profs = 1:num_profs;
                        good_profs = setdiff(all_profs, bad_profs);
                                                
                        ['lat/lon check: outside lat lims, outside lon lims, zero lat/lon']
                        [length(outside_lat_lim_ins) ...
                            length(outside_lon_lim_ins) ...
                            length(zero_lat_lon_ins)]
                        
                        ['all profs, good profs, bad_profs']
                        [length(all_profs) length(good_profs) length(bad_profs)]
                        
                        % set the weights of profiles that are out of
                        % bounds to zero.
                       
                        MITprof_new.prof_Tweight(bad_profs,:) = 0;
                        MITprof_new.prof_Sweight(bad_profs,:) = 0;
                        
                        MITprof = MITprof_new;
                        clear MITprof_new;
                        
                    case 5
                        %%
                        %%
                        if plot_individual_bad_profiles
                            prof_fig_dir = [fig_dir '/bad_profs_' fDataOut{ilist}(1:end-3)];
                            mkdir(prof_fig_dir);
                        end
                                
                        clear bad_profs* good_profs*
                        
                        %% PART 1, set weights to zero for profiles with 
                        % large average costs vs. clim.
                        
                        T_cost_vs_clim = ...
                            (MITprof.prof_T - MITprof.prof_Tclim).^2 .*...
                            MITprof.prof_Tweight;
                        
                        S_cost_vs_clim = ...
                            (MITprof.prof_S - MITprof.prof_Sclim).^2 .*...
                            MITprof.prof_Sweight;
                        
                        % calculate the average cost for each profile
                        T_no_data_ins = find(MITprof.prof_T == fillVal);
                        S_no_data_ins = find(MITprof.prof_S == fillVal);
                        
                        T_cost_vs_clim(T_no_data_ins) = NaN;
                        S_cost_vs_clim(S_no_data_ins) = NaN;
                        
                        T_zero_weight_ins = find(MITprof.prof_Tweight == 0);
                        S_zero_weight_ins = find(MITprof.prof_Sweight == 0);
                        
                        T_cost_vs_clim(T_zero_weight_ins) = NaN;
                        S_cost_vs_clim(S_zero_weight_ins) = NaN;
                        
                        % average T and S costs of the profile
                        avg_T_costs = mynanmean(T_cost_vs_clim,2);
                        avg_S_costs = mynanmean(S_cost_vs_clim,2);
                        
                        if cost_vs_clim_filter_method == 0
                            % use statistics method
                            
                            % cost bins.
                            cost_rng = 0:1:1000;
                            % get the number of profiles whose average costs
                            % fall within the cost bins.
                            ins = find(avg_T_costs > 0);
                            [a,b]=hist(avg_T_costs(ins), cost_rng);
                            % do a cumulative prob. of costs in these bins.
                            c=cumsum(a)./sum(a);
                            
                            % those cost values from the cost bins that
                            % are in the bottom 'cost_threshold' 99.5% of all
                            % costs
                            c_thresh = find(c <= cost_threshold);
                            
                            % the bin index of the highest cost with the bottom
                            % 99.5%
                            c_thresh_i = c_thresh(end)+1;
                            
                            % the maximum average T or S cost that we'll
                            % tolerate
                            avg_cost_max = cost_rng(c_thresh_i);
                            AVG_T_C_THRESH = avg_cost_max;
                            
                            % profs that are bad because their average
                            % salinity cost exceeds some cost_threshold
                            bad_T_profs = find(avg_T_costs >= avg_cost_max);
                        
                        elseif cost_vs_clim_filter_method == 1
                            % use cost_threshold for maximum cost
                            bad_T_profs = find(avg_T_costs >= cost_threshold);
                        end
                        
                        ['num profs, bad T profs (avg)']
                        [num_profs length(bad_T_profs)]
                         
                        if cost_vs_clim_filter_method == 0
                            %% do the same for salinity
                            ins = find(avg_S_costs > 0);
                            [a,b]=hist(avg_S_costs(ins), cost_rng);
                            c=cumsum(a)./sum(a);
                            c_thresh = find(c <= cost_threshold);
                            c_thresh_i = c_thresh(end) +1;
                            avg_cost_max = cost_rng(c_thresh_i);
                            
                            bad_S_profs = find(avg_S_costs >= avg_cost_max);
                        elseif cost_vs_clim_filter_method == 1
                            bad_S_profs = find(avg_S_costs >= cost_threshold);
                        end
                        
                        ['num profs, bad S profs (avg)']
                        [num_profs length(bad_S_profs)]
                        
                        
                        %% show where we have bad T and S profiles 
                        if make_figs
                            close all;
                            %% T costs
                            if plot_individual_bad_profiles
                                tmp = min(length(bad_T_profs), ...
                                    num_bad_profs_to_plot);
                                
                                for bi = 1:tmp
                                    b_here = bad_T_profs(bi);
                                    close all;
                                    figure(1);clf;subplot(131);hold on;
                                    plot_compare_prof_vs_clim_cost(MITprof,...
                                        b_here, ['prof ' num2str(b_here)...
                                        '   avg T cost ' ...
                                        num2str(avg_T_costs(b_here))]);
                                    cd(prof_fig_dir)
                                    
                                    titA =  fDataOut{ilist}(1:end-3);
                                    pltstmp(gcf,0,titA);
                                    paperx=10;papery=6;prep_figure_for_print;
                                    fig_name = ['prof_' padzero(b_here,6) ...
                                        '_high_avg_T_cost.png'];
                                   
                                    print(gcf,'-dpng','-r200',fig_name);
                                end
                            end
                            %%
                            if plot_map_bad_profiles
                                figure(1000);clf;
                                
                                show_profile_locations_3panel(...
                                    MITprof.prof_lon(bad_T_profs,:),...
                                    MITprof.prof_lat(bad_T_profs,:));
                                
                                cd(fig_dir);
                                titA =  fDataOut{ilist}(1:end-3);
                                titB =  ['PROFS WITH TOO HIGH AVG T COST n=' ...
                                    num2str(length(bad_T_profs))];
                                suptitle(titB);
                                pltstmp(gcf,0,titA);
                                paperx=6;papery=6;prep_figure_for_print;
                                fig_name = [fDataOut{ilist}(1:end-3) ...
                                    '_high_avg_T_cost.png']
                                print(gcf,'-dpng','-r200',fig_name);
                            end
                            
                            %% S costs
                            if plot_individual_bad_profiles
                                tmp = min(length(bad_S_profs), ...
                                    num_bad_profs_to_plot);
                                
                                for bi = 1:tmp
                                    b_here = bad_S_profs(bi);
                                    close all;
                                    figure(1);clf;subplot(131);hold on;
                                    plot_compare_prof_vs_clim_cost(MITprof,...
                                        b_here, ['prof ' num2str(b_here)...
                                        '   avg S cost ' ...
                                        num2str(avg_S_costs(b_here))]);
                                    cd(prof_fig_dir)
                                    
                                    titA =  fDataOut{ilist}(1:end-3);
                                    pltstmp(gcf,0,titA);
                                    paperx=10;papery=6;prep_figure_for_print;
                                    fig_name = ['prof_' padzero(b_here,6) ...
                                        '_high_avg_S_cost.png'];
                                   
                                    print(gcf,'-dpng','-r200',fig_name);
                                end
                            end
                            %%
                            if plot_map_bad_profiles
                                figure(1001);clf;
                                show_profile_locations_3panel(...
                                    MITprof.prof_lon(bad_S_profs,:),...
                                    MITprof.prof_lat(bad_S_profs,:));
                                cd(fig_dir);
                                titA =  fDataOut{ilist}(1:end-3);
                                titB =  ['PROFS WITH TOO HIGH AVG S COST n=' ...
                                    num2str(length(bad_S_profs))];
                                suptitle(titB);
                                pltstmp(gcf,0,titA);
                                paperx=6;papery=6;prep_figure_for_print;
                                fig_name = [fDataOut{ilist}(1:end-3) ...
                                    '_high_avg_S_cost.png']
                                print(gcf,'-dpng','-r200',fig_name);                                
                            end
                        end
                        
                        %% now update the weights for those bad T and
                        % S profiles
                        MITprof_new.prof_Tweight(bad_T_profs,:) = 0;
                        MITprof_new.prof_Sweight(bad_S_profs,:) = 0;
                        
                        % calculate new costs based on the new weights
                        T_cost_vs_clim = ...
                            (MITprof.prof_T - MITprof.prof_Tclim).^2 .*...
                            MITprof_new.prof_Tweight;
                        
                        S_cost_vs_clim = ...
                            (MITprof.prof_S - MITprof.prof_Sclim).^2 .*...
                            MITprof_new.prof_Sweight;
                      
                        %% PART 2, remove profiles that have points in the 
                        % profile with costs that exceed some cost_threshold
                        % these are profiles that may have sections of 
                        % bad data that for whatever reason weren't flagged
                        % earlier.
                        
                        % there is an option to exclude profiles from
                        % high latitudes because climatologies at high
                        % latitudes are often unreliable
                            
                        clear bad*
                        
                        % calculate the average cost for each profile
                        T_no_data_ins = find(MITprof.prof_T == fillVal);
                        S_no_data_ins = find(MITprof.prof_S == fillVal);
                        
                        T_cost_vs_clim(T_no_data_ins) = NaN;
                        S_cost_vs_clim(S_no_data_ins) = NaN;
                        
                        T_zero_weight_ins = find(MITprof_new.prof_Tweight == 0);
                        S_zero_weight_ins = find(MITprof_new.prof_Sweight == 0);
                        
                        T_cost_vs_clim(T_zero_weight_ins) = NaN;
                        S_cost_vs_clim(S_zero_weight_ins) = NaN;
                        
                        maxTc = max(T_cost_vs_clim,[],2);
                        maxSc = max(S_cost_vs_clim,[],2);

                        %%
                        % find profs that are bad because they have an T cost
                        % somewhere that exceeds some prob. cost_threshold
                        
                        if cost_vs_clim_filter_method == 0
                            % use statistics
                            ins = find(maxTc > 0);
                            [a,b]=hist(maxTc(ins), cost_rng);
                            c=cumsum(a)./sum(a);
                            c_thresh = find(c <= cost_threshold);
                            c_thresh_i = c_thresh(end) +1;
                            % the highest T cost in a profile that can be
                            % tolerated
                            max_T_cost_max = cost_rng(c_thresh_i);
                            bad_T_profs_a = find(maxTc >= max_T_cost_max);
                            
                        elseif cost_vs_clim_filter_method == 1
                            % use a straight number maximum cost
                            bad_T_profs_a = find(maxTc >= cost_threshold);
                        end
                        
                        if cost_vs_clim_filter_method == 0
                            % profs that are bad because they have an S cost
                            % somewhere that exceeds some prob. cost_threshold
                            
                            ins = find(maxSc > 0);
                            [a,b]=hist(maxSc(ins), cost_rng);
                            c=cumsum(a)./sum(a);
                            c_thresh = find(c <= cost_threshold);
                            c_thresh_i = c_thresh(end) +1;
                            % the highest S cost in a profile that can be
                            % tolerated
                            max_S_cost_max = cost_rng(c_thresh_i);
                            
                            %% identify T and S profiles that have some cost
                            % that exceeds the max_T_cost_max or max_S_cost_max
                            bad_S_profs_a = find(maxSc >= max_S_cost_max);

                        elseif cost_vs_clim_filter_method == 1
                            bad_S_profs_a = find(maxSc >= cost_threshold);
                        end
                        
                        if exclude_high_latitude_profiles_from_test_2
                            % findn all profiles that are outside of high
                            % latitudes (-60 to 60)
                            profs_outside_high_lats = ...
                                find(MITprof.prof_lat >= -60 & ...
                                MITprof.prof_lat <=  60);
                            
                            % only remove those profiles that are outside of
                            % high latitudes and have points with costs
                            % exceeding the maximum tolerated values
                            bad_S_profs = intersect(profs_outside_high_lats, ...
                                bad_S_profs_a);
                            
                            bad_T_profs = intersect(profs_outside_high_lats, ...
                                bad_T_profs_a);
                        else
                            bad_S_profs = bad_S_profs_a;
                            bad_T_profs = bad_T_profs_a;
                        end
                        
                        ['PART 2: num bad T and S profs']
                        [length(bad_T_profs) length(bad_S_profs)]
                        
                        %% show where we have bad T and S profiles 
                        if make_figs
                            close all;
                            %% T costs
                            if plot_individual_bad_profiles
                                tmp = min(length(bad_T_profs), ...
                                    num_bad_profs_to_plot);
                                
                                for bi = 1:tmp
                                    %%
                                    b_here = bad_T_profs(bi);
                                    close all;
                                    set(0,'DefaultTextInterpreter','none');
                                    plot_compare_prof_vs_clim_cost(MITprof,...
                                        b_here, ['prof ' num2str(b_here)...
                                        '   max T cost ' ...
                                        num2str(maxTc(b_here))]);
                                    cd(prof_fig_dir)
                                    
                                    titA =  fDataOut{ilist}(1:end-3);
                                    pltstmp(gcf,0,titA);
                                    paperx=10;papery=6;prep_figure_for_print;
                                    fig_name = ['prof_' padzero(b_here,6) ...
                                        '_high_single_value_T_cost.png'];
                                   
                                    print(gcf,'-dpng','-r200',fig_name);
                                end
                            end
                            %%
                            if plot_map_bad_profiles
                                close all;
                                show_profile_locations_3panel(...
                                    MITprof.prof_lon(bad_T_profs,:),...
                                    MITprof.prof_lat(bad_T_profs,:));
                                
                                cd(fig_dir);
                                titA =  fDataOut{ilist}(1:end-3);
                                titB =  ['PROFS WITH TOO HIGH SINGLE VALUE T COST n=' ...
                                    num2str(length(bad_T_profs))];
                                suptitle(titB);
                                pltstmp(gcf,0,titA);
                                paperx=6;papery=6;prep_figure_for_print;
                                fig_name = [fDataOut{ilist}(1:end-3) ...
                                    '_high_single_value_T_cost.png']
                                print(gcf,'-dpng','-r200',fig_name);
                            end
                            
                            %% S costs
                            if plot_individual_bad_profiles
                                tmp = min(length(bad_S_profs), ...
                                    num_bad_profs_to_plot);
                                
                                for bi = 1:tmp
                                    %%
                                    b_here = bad_S_profs(bi);
                                    close all;
                                    set(0,'DefaultTextInterpreter','none');
                                    
                                    plot_compare_prof_vs_clim_cost(MITprof,...
                                        b_here, ['prof ' num2str(b_here)...
                                        '   max S cost ' ...
                                        num2str(maxSc(b_here))]);
                                    cd(prof_fig_dir)
                                    
                                    titA =  fDataOut{ilist}(1:end-3);
                                    pltstmp(gcf,0,titA);
                                    paperx=10;papery=6;prep_figure_for_print;
                                    fig_name = ['prof_' padzero(b_here,6) ...
                                        '_high_single_value_S_cost.png'];
                                   
                                    print(gcf,'-dpng','-r200',fig_name);
                                end
                            end
                            
                            if plot_map_bad_profiles
                                %%
                                close all;
                                show_profile_locations_3panel(...
                                    MITprof.prof_lon(bad_S_profs,:),...
                                    MITprof.prof_lat(bad_S_profs,:));
                                cd(fig_dir);
                                titA =  fDataOut{ilist}(1:end-3);
                                titB =  ['PROFS WITH TOO HIGH SINGLE VALUE S COST n=' ...
                                    num2str(length(bad_S_profs))];
                                suptitle(titB);
                                pltstmp(gcf,0,titA);
                                paperx=6;papery=6;prep_figure_for_print;
                                fig_name = [fDataOut{ilist}(1:end-3) ...
                                    '_high_single_value_S_cost.png']
                                print(gcf,'-dpng','-r200',fig_name);                                
                            end
                        end
                        
                        % set the weights of profiles that are out of
                        % bounds to zero.
                       
                        MITprof_new.prof_Tweight(bad_T_profs,:) = 0;
                        MITprof_new.prof_Sweight(bad_S_profs,:) = 0;
                        
                        MITprof = MITprof_new;
                        clear MITprof_new;
                        
                end
            end
            %%
            
            if debug_code
                ['after processing # nonzero T and S weights']
                [length(find(MITprof.prof_Tweight > 0)) ...
                    length(find(MITprof.prof_Tweight > 0))]
            end
            
            %%
            if length(MITprof.prof_YYYYMMDD) > 0
                %%
                %  Write output
                ['writing output']
                fileOut=[output_dir '/' fDataOut{ilist}];
                fprintf('%s\n',fileOut);
                
                write_profile_structure_to_netcdf(MITprof,fileOut);
            else
                ['THERE ARE NO GOOD PROFILES LEFT']
            end
        else
            [fDataIn{ilist} ' does not exist']
        end
    else
        [fDataOut{ilist} ' already exists, skipping!']
    end
end % ilist
