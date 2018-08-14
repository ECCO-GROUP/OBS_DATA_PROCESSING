function update_prof_and_tile_points_on_profiles(run_codes)

% This script updates the prof_points and all tile interpolation point
% fields and sends the updated profile files to be written in netcdf format

% Filename; update_prof_and_tile_points_on_profiles.m
%  ** former filename :
% Date Created: 2014-07-30
% Last Modified: 2018-05-10;
%


% notes:
%    this script exists to allow one to use the same MITprofile package
%    netcdf files in different llc configurations easily.
%
%    2014-08-06 : added some error checking and the ability to specify that
%                 you want to interpolate to only wet points automatically,
%                 if desired.

%    2015-02-05 : add split the netcf file by year in output
%    2015-02-09 : removed split by year in output. it doesn't make sense
%    here actually
%    2016-03-08 : added run_code.. dealing with 'argo_feb2016_setX' files
%    2016-03-11 : aadded text to debuggin figure to denote which profiles
%                 are bad. also added run code to convert latest llc90 to llc270
%    2016-08-17 : converted to a function
%    2018-05-10 : removed script as subroutine, changed write to netcdf
%    arguments
%    2018-06-06 : converted back to a function.

set(0,'DefaultTextInterpreter','none');
%%

make_root_dir;

%%---------------------------------------------
%% SET INPUT PARAMETERS
%%---------------------------------------------
% wet_or_all
% 0: only interpolate to nearest wet point
% 1: interpolate to all points.
wet_or_all = [];%1

% specify the llc code
llcN = [];

clear yrs;

for rci = 1:length(run_codes)
    run_code = run_codes{rci}
    
    clear fDataBase;
    clear fDataIn* fDataOut;
    
    switch run_code
        case '20180610_llc90_argo_all'
            root_dir = '/home/ifenty/data/observations/insitu/ARGO/from_gael_June2018/combined_by_latest'
            input_dir= [root_dir];
            output_dir = [root_dir '/step_00_update_tile_and_prof_points/'];
            
            file_suffix_in = '.nc';
            file_suffix_out = '_step_00.nc'

            cd(input_dir);
            ls
            f = dir(['*nc']);
            f
            for i = 1:length(f)
                fDataIn{i}=   f(i).name;
                fDataOut{i} =[f(i).name(1:end-length(file_suffix_in)) file_suffix_out];
            end
            
            % llc model grid
            llcN = 90;
            wet_or_all = 1;
            
            make_figs = 0;
            
        case '20180508_llc90_NODC_all'
            root_dir = '/home/ifenty/data/observations/insitu/NODC/NODC_20180508/all/'
            input_dir= [root_dir '/first_conversion_to_MITprof/'];
            output_dir = [root_dir '/step_00_update_tile_and_prof_points/'];
            
            file_suffix_in = '.nc';
            file_suffix_out = '_step_00.nc'

            cd(input_dir);
            f = dir(['*nc']);
            for i = 1:length(f)
                fDataIn{i}=   f(i).name;
                fDataOut{i} =[f(i).name(1:end-length(file_suffix_in)) file_suffix_out];
            end
            
            % llc model grid
            llcN = 90;
            wet_or_all = 1;
            
            make_figs = 0;
            
    end
    
    ilists = 1:length(f);
    
    mkdir(output_dir)
    
    %%---------------------------------------------
    %%  Init profile package
    %%---------------------------------------------
    MITprof_path;
    
    
    %%---------------------------------------------
    %%  Read in grid for climatology
    %%---------------------------------------------
    
    nx=llcN;ny=nx*13;nz=50;
    
    switch llcN
        case 90
            load_llc90_grid
            ni = 30;
            nj = 30;
            lon_llc = lon_90;
            lat_llc = lat_90;
            llc_grid_dir = llc90_grid_dir;
            mask_llc = blank_90;
            
            if wet_or_all == 0
                'wet points'
                mask_llc(wet_ins_90_k{1}) = 1;
            else
                'all points'
                mask_llc=ones(size(blank_90));
            end
            
        case 270
            load_llc270_grid
            ni = 30;
            nj = 30;
            lon_llc = lon_270;
            lat_llc = lat_270;
            llc_grid_dir = llc270_grid_dir;
            mask_llc = blank_270;
            
            if wet_or_all ==0
                mask_llc(wet_ins_270_k{1}) = 1;
            else
                mask_llc=ones(size(blank_270));
            end
    end
    
    
    %%---------------------------------------------
    %% Read and process the profile files
    %%---------------------------------------------
    % F is the mapping -- should probably remake this every time.
    clear F;
    
    %%
    for ilist = ilists
        close all;
        fDataIn{ilist}
        cd(input_dir)
        
        clear fileOut MITprof
        
        MITprof = MITprof_read(fDataIn{ilist});
        
        if exist('F')
            [MITprof_pp, F] = ...
                get_profpoint_llc_ian(lon_llc, lat_llc, ...
                mask_llc,  nx, MITprof,F);
        else
            [MITprof_pp, F] = ...
                get_profpoint_llc_ian(lon_llc, lat_llc, ...
                mask_llc,  nx, MITprof);
        end
        
        
        MITprof_tp = get_tile_point_llc_ian(lon_llc, lat_llc, ...
            ni, nj, nx, MITprof_pp);
        
        
        %%
        MITprof = MITprof_tp
        clear MITprof_pp MITprof_tp;
        %%
        if make_figs
            figure(ilist);clf;hold on;
            set(0,'DefaultTextInterpreter','none');
            plot(MITprof.prof_lon, MITprof.prof_lat,'r.')
            plot(MITprof.prof_interp_lon, MITprof.prof_interp_lat,'bo')
            grid;
            
            figure;
            set(0,'DefaultTextInterpreter','none');
            projection_code=1;
            plot_proj=0;
            makemmap_general;
            m_plot(MITprof.prof_lon, MITprof.prof_lat,'r.')
            m_plot(MITprof.prof_interp_lon, MITprof.prof_interp_lat,'bo')
            makemmap_general;
            m_grid;
        end
        
        %% check max distances
        
        tmp_prof_lat = MITprof.prof_lat;
        tmp_prof_lon = MITprof.prof_lon;
        
        bad_lats = find(abs(tmp_prof_lat)>90);
        
        tmp_prof_lat(find(tmp_prof_lat>90))=90;
        tmp_prof_lat(find(tmp_prof_lat<-90))=-90;
        
        
        d = vdist(tmp_prof_lat, tmp_prof_lon, ...
            MITprof.prof_interp_lat, MITprof.prof_interp_lon);
        
        d=d./1e3;  % distance now in KM.
        
        switch llcN
            case 270
                dx=sqrt(RAC_270_pf(MITprof.prof_point))./1e3;
            case 90
                dx=sqrt(RAC_90_pf(MITprof.prof_point))./1e3;
        end
        ins = find(d./dx > 1);
        
        % 
        MITprof.prof_flag(bad_lats) = 0;
        MITprof.prof_flag(ins) = 0;
        
        %%
        if length(ins) > 0
            figure(ilist)
            set(0,'DefaultTextInterpreter','none');
            h1=plot(MITprof.prof_lon(ins), MITprof.prof_lat(ins),'go','MarkerSize',10)
            
            set(h1(:), 'MarkerEdgeColor','k','MarkerFaceColor','g');
            h2=plot(MITprof.prof_interp_lon(ins), MITprof.prof_interp_lat(ins),'yo','MarkerSize',10)
            set(h2(:), 'MarkerEdgeColor','k','MarkerFaceColor','y');
            
            title({'points that are too far away', fDataIn{ilist}})
        end
        %%
        
        ['writing ' fDataOut{ilist}]
        fileOut=[output_dir '/'  fDataOut{ilist}];
        fprintf('%s\n',fileOut);
        write_profile_structure_to_netcdf(MITprof, fileOut);
    end
end
% 
% 
%  case '20180508_llc90_wod13_CTD'
%             tmp ='CTD'
%             input_dir = ['/ian4/ifenty/data/observations/insitu/NODC/NODC_20180508/' tmp '/' tmp '_MITPROF/first_conversion_to_MITprof_format/']
%             output_dir = ['/ian4/ifenty/data/observations/insitu/NODC/NODC_20180508/' tmp '/' tmp '_MITPROF//']
%             file_suffix_in  = '.nc'
%             file_suffix_out = 'step_00.nc'
%             
%             yrs = 1992:2018
%             
%             cd(input_dir);
%             f = dir(['*nc']);
%             for i = 1:length(f)
%                 fDataIn{i}=   f(i).name;
%                 fDataOut{i} =[f(i).name(1:end-length(file_suffix_in)) file_suffix_out];
%             end
%             ilists = 1:length(f);
%             
%             % llc to convert to
%             llcN = 90;
%             wet_or_all = 1;
%             
%         case '20180508_llc90_wod13_GLD'
%             tmp ='GLD'
%             input_dir = ['/ian4/ifenty/data/observations/insitu/NODC/NODC_20180508/' tmp '/' tmp '_MITPROF/first_conversion_to_MITprof_format/']
%             output_dir = ['/ian4/ifenty/data/observations/insitu/NODC/NODC_20180508/' tmp '/' tmp '_MITPROF/step_00_update_tile_and_prof_points/']
%             file_suffix_in  = '.nc'
%             file_suffix_out = 'step_00.nc'
%             
%             yrs = 1992:2018
%             
%             cd(input_dir);
%             f = dir(['*nc']);
%             for i = 1:length(f)
%                 fDataIn{i}=   f(i).name;
%                 fDataOut{i} =[f(i).name(1:end-length(file_suffix_in)) file_suffix_out];
%             end
%             ilists = 1:length(f);
%             
%             % llc to convert to
%             llcN = 90;
%             wet_or_all = 1;
%             
%         case '20180508_llc90_wod13_APB'
%             tmp ='APB'
%             input_dir = ['/ian4/ifenty/data/observations/insitu/NODC/NODC_20180508/' tmp '/' tmp '_MITPROF/first_conversion_to_MITprof_format/']
%             output_dir = ['/ian4/ifenty/data/observations/insitu/NODC/NODC_20180508/' tmp '/' tmp '_MITPROF/step_00_update_tile_and_prof_points/']
%             file_suffix_in  = '.nc'
%             file_suffix_out = 'step_00.nc'
%             
%             yrs = 1992:2018
%             
%             cd(input_dir);
%             f = dir(['*nc']);
%             for i = 1:length(f)
%                 fDataIn{i}=   f(i).name;
%                 fDataOut{i} =[f(i).name(1:end-length(file_suffix_in)) file_suffix_out];
%             end
%             ilists = 1:length(f);
%             
%             % llc to convert to
%             llcN = 90;
%             wet_or_all = 1;
%             
%         case '20180508_llc90_wod13_XBT'
%             tmp ='XBT'
%             input_dir = ['/ian4/ifenty/data/observations/insitu/NODC/NODC_20180508/' tmp '/' tmp '_MITPROF/first_conversion_to_MITprof_format/']
%             output_dir = ['/ian4/ifenty/data/observations/insitu/NODC/NODC_20180508/' tmp '/' tmp '_MITPROF/step_00_update_tile_and_prof_points/']
%             file_suffix_in  = '.nc'
%             file_suffix_out = 'step_00.nc'
%             
%             yrs = 1992:2018
%             
%             cd(input_dir);
%             f = dir(['*nc']);
%             for i = 1:length(f)
%                 fDataIn{i}=   f(i).name;
%                 fDataOut{i} =[f(i).name(1:end-length(file_suffix_in)) file_suffix_out];
%             end
%             
%             ilists = 1:length(f);
%             
%             % llc to convert to
%             llcN = 90;
%             wet_or_all = 1;
%             
%         case '20180508_llc90_wod13_MRB'
%             tmp ='MRB'
%             input_dir = ['/ian4/ifenty/data/observations/insitu/NODC/NODC_20180508/' tmp '/' tmp '_MITPROF/first_conversion_to_MITprof_format/']
%             output_dir = ['/ian4/ifenty/data/observations/insitu/NODC/NODC_20180508/' tmp '/' tmp '_MITPROF/step_00_update_tile_and_prof_points/']
%             file_suffix_in  = '.nc'
%             file_suffix_out = 'step_00.nc'
%             
%             cd(input_dir)
%             yrs = 1992:2018
%             
%             cd(input_dir);
%             f = dir(['*nc']);
%             for i = 1:length(f)
%                 fDataIn{i}=   f(i).name;
%                 fDataOut{i} =[f(i).name(1:end-length(file_suffix_in)) file_suffix_out];
%             end
%             
%             ilists = 1:length(f);
%             
%             % llc to convert to
%             llcN = 90;
%             wet_or_all = 1;