

% This script takes monthly-mean sea ice fields derived from the goddard
% merged mean field from the following dataset and converts it to binary
% for easy use.  see below

% the dataset contains 4 different sea ice products
% (1) NOAA/NSIDC CDR
% (2) NASA Bootstrap
% (3) NASA team
% (4) merged Bootstram + NASA team

% http://nsidc.org/data/docs/noaa/g02202_ice_conc_cdr/index.html#data_description

% Filename; save_noaa_nsidc_monthly_sea_ice_concentration_to_llcGrid_binary.m
%  ** former filename :  convert_monthly_mean_matlab_to_binary.

% Date Created: unknown.
% Last Modified: 2016-03-14

% notes
%   2016-03-14  former version seems to have only saved one month per year as a bug.


clear;
close all;

set(0,'DefaultTextInterpreter','none');

run_codes=[0]
ice_dir = '~/data/observations/seaice/NOAA_NSIDC_CLIM_DATA_REC_OF_PM_SEA_ICE/'

for run_code=run_codes
    
    switch run_code
        case 0
            llcN=90;
            years=2012:2015;
            makefigs=1
            input_dir = '/ian1/ifenty/data/observations/seaice/NOAA_NSIDC_CLIM_DATA_REC_OF_PM_SEA_ICE/monthly_projected_to_matlab/llc090/'
            output_dir = '/ian1/ifenty/data/observations/seaice/NOAA_NSIDC_CLIM_DATA_REC_OF_PM_SEA_ICE/monthly_projected_to_binary/llc090/'
            fig_dir =  '/ian1/ifenty/data/observations/seaice/NOAA_NSIDC_CLIM_DATA_REC_OF_PM_SEA_ICE/monthly_figures/llc090/';
            
        case 1
            llcN=270;
            years=2013;
            makefigs=1
            input_dir = '/ian1/ifenty/data/observations/seaice/NOAA_NSIDC_CLIM_DATA_REC_OF_PM_SEA_ICE/monthly_projected_to_matlab/llc270/'
            output_dir = '/ian1/ifenty/data/observations/seaice/NOAA_NSIDC_CLIM_DATA_REC_OF_PM_SEA_ICE/monthly_projected_to_binary/llc270/'
            fig_dir =  '/ian1/ifenty/data/observations/seaice/NOAA_NSIDC_CLIM_DATA_REC_OF_PM_SEA_ICE/monthly_figures/llc270/';
        otherwise
            ['need to specify a run code']
    end
    
    
    
    %%
    
    %%
    for year = years
        
        [year]
        tic
        cd(input_dir)
        fname = ['NOAA_NSIDC_MONTHLY_MAPPED_TO_LLC' num2str(llcN) '_' num2str(year)]
        load(fname)
        
        
        sea_ice_conc = NOAA_NSIDC_MONTHLY_MAPPED_TO_LLC.goddard_merged_mean./100;
        
        % any sea ice concentraiton value less than 1 should be masked
        sea_ice_conc(find(sea_ice_conc < 0))=-9999;
        
        cd(output_dir)
        fid = fopen([fname '.bin'],'w','b');
        fwrite(fid, sea_ice_conc,'single');
        fclose(fid);
        toc
    end
end
