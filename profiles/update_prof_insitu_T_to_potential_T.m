function update_prof_insitu_T_to_potential_T(run_code)

% This script updates the profile insitu temperatures so that they are in
% potential temperature

% Filename; update_prof_insitu_T_to_potential_T.m
%  ** former filename : 
% Date Created: 2014-08-04
% Last Modified: 2018-06-04

% notes:

% 2018-06-04, made this look like other processing scripts.  run_code
% also added the ability to detect files that had zero salinity data
% and to then use climatology salinity data if it exists.
%

write_profile_to_netcdf = 1
make_figs = 1

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

TS_clim_dir = ''
f_source_ST={}

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
    case  'NODC_20180508'
        ['apply monthly clim to mitprof variable']
        rootdir = ['/home/ifenty/data/observations/insitu/NODC/NODC_20180508/all/']
        input_dir = [rootdir 'step_05_zero_weights'];
        
        
        file_suffix_in  = ['step_05.nc'];
        file_suffix_out = ['step_06.nc']
        
        output_dir =  [rootdir 'step_06_potential_temperature']
        
        fig_dir = [output_dir '/figures']
       
        cd(input_dir);
        f = dir(['*nc*']);
        
        for i = 1:length(f)
            fDataIn{i}  =   f(i).name;
            fDataOut{i} =  [f(i).name(1:end-length(file_suffix_in)) file_suffix_out];
        end
        ilists = 1:length(f);
        %%
end
mkdir(output_dir);
mkdir(fig_dir);


%%---------------------------------------------
%% Read and process the profile files
%%---------------------------------------------

for ilist = ilists
    
    %%
    cd(input_dir)
    
    
    clear fileOut MITprof
    fileOut=[output_dir '/' fDataOut{ilist}]
    %%
    if ~exist(fileOut,'file')
        
        cd(input_dir)
        MITprof = MITprof_read(fDataIn{ilist});
        MITprof_orig = MITprof;
        
        %
        lats = MITprof.prof_lat;
        prof_T = MITprof.prof_T;
        prof_T_weight = MITprof.prof_Tweight;
        
        % salinity
        prof_S = MITprof.prof_S;
        prof_S_weight = MITprof.prof_Sweight;
        
        % Check to see if **all** salinty values are missing
        S_max=max(prof_S');
       
        if max(S_max) == fillVal
            ['all S are missing']
            prof_S = MITprof.prof_Sclim;
            prof_S_weight = ones(size(prof_S));
        end
        
        % to qualify you need to have a valid T, S, Tweight and Sweight
        good_T_and_S_ins = find(prof_T ~= fillVal & prof_S ~= fillVal & ...
            prof_T_weight > 0 & prof_S_weight > 0);
        
        % Make empty arryas.
        prof_T_tmp = ones(size(prof_T)).*NaN;
        prof_S_tmp = ones(size(prof_S)).*NaN;
        
        % set values at the good T and S pairs to be the original T and S
        prof_T_tmp(good_T_and_S_ins) = prof_T(good_T_and_S_ins);
        prof_S_tmp(good_T_and_S_ins) = prof_S(good_T_and_S_ins);
       
        % define an empty ptemp;
        ptemp = ones(size(prof_T)).*fillVal;
        
        % define a new T weight that is all zeros.
        prof_T_weight_new = zeros(size(prof_T));
        
        
        if length(good_T_and_S_ins) > 0
            % set the new T weights to be the same as prof_T where we have 
            % both good T and S data.
            
            prof_T_weight_new(good_T_and_S_ins) = prof_T_weight(good_T_and_S_ins);
                    
            %% Prepare 2D matrix of pres and lats required for sw_ptmp
            depths = MITprof.prof_depth;
            depths_mat = repmat(depths, [1 length(lats)]);
            
            lats_mat = repmat(lats, [1 length(depths)])';
            
            % calculate equivalent pressure from depth
            pres_mat = sw_pres(depths_mat, lats_mat);
            pres_mat = pres_mat';
            
            %%
            % Calc potential temperature w.r.t. to surf [pres = 0]
            ptemp = sw_ptmp(prof_S_tmp, prof_T_tmp, pres_mat, pres_mat.* 0);
            
            %%
            if make_figs
                %%
                cmap = jet(100);
                cmap(1,:) = [1 1 1];
                
                close all;
                figure(1);clf;
                
                subplot(411);
                imagesc(MITprof_orig.prof_T');
                min_T = min(ptemp(good_T_and_S_ins));
                max_T = max(ptemp(good_T_and_S_ins));
                caxis([min_T max_T]);
                colorbar;
                
                subplot(412);
                
                imagesc(prof_T_tmp'); colorbar;
                caxis([min_T max_T]);
                title('theta');
                %%
                subplot(413);
                imagesc(ptemp'); caxis([min_T max_T]);
                title('pot theta');colorbar;
                
                subplot(414);
                delta = ptemp - prof_T_tmp;
                imagesc(delta')
                title('pot theta - theta')
                caxis([-.1 .1]);colorbar;
                
                colormap(cmap)
                cd(fig_dir)
                paperx=10;papery=10;prep_figure_for_print;
                fig_name = [fDataOut{ilist}(1:end-3) '.png']
                print(gcf,'-dpng','-r200',fig_name);
                
            end
          
        else
            ['There is not a single good T and S pair to use here']
            ['not saving anything!']
        end
        
        MITprof.prof_T = ptemp;
        MITprof.prof_Tweight = prof_T_weight_new;
        
        if write_profile_to_netcdf & length(good_T_and_S_ins) > 0
            fileOut=[output_dir '/' fDataOut{ilist}];
            fprintf('%s\n',fileOut);
            write_profile_structure_to_netcdf(MITprof, fileOut);
        else
            ['I am not saving anything, the rest is up to you']
        end % write
    else
        [fileOut ' exists, skipping!'];
    end
end
