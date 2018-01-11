%%
 
% This script saves sea ice concentration observations projected to the llc
% grid in 'interp_noaa_nsidc_daily_sea_ice_concentration_to_llcGrid' to
% flat binary format for use within the MITgcm.

% Filename; save_noaa_nsidc_daily_sea_ice_concentration_to_llcGrid_binary.m
%  ** former filename : 
% Date Created: 2014-10-06
% Last Modified:  2017-01-19


% notes:

% 2015-01-19: updated to include 2015 data
% 2014-10-06: updated input directory
% 2016-03-14: converted to run_code format.
% 2016-03-17: added capability of making some dummy empty years when no
% data is avail.

clear all;

run_codes=[0]
ice_dir = '~/data/observations/seaice/NOAA_NSIDC_CLIM_DATA_REC_OF_PM_SEA_ICE/'

for run_code=run_codes
    
    switch run_code
        case 0
            llcN=90;
            years=2012:2015
            makefigs=1
            input_dir = '/ian1/ifenty/data/observations/seaice/NOAA_NSIDC_CLIM_DATA_REC_OF_PM_SEA_ICE/daily_projected_to_matlab/llc090/'
            output_dir = '/ian1/ifenty/data/observations/seaice/NOAA_NSIDC_CLIM_DATA_REC_OF_PM_SEA_ICE/daily_projected_to_binary/llc090/'
            
            zero_years = 2016:2017;
        case 1
            llcN=270;
            years=2015
            makefigs=1
            input_dir = '/ian1/ifenty/data/observations/seaice/NOAA_NSIDC_CLIM_DATA_REC_OF_PM_SEA_ICE/daily_projected_to_matlab/llc270/'
            output_dir = '/ian1/ifenty/data/observations/seaice/NOAA_NSIDC_CLIM_DATA_REC_OF_PM_SEA_ICE/daily_projected_to_binary/llc270/'
            
            zero_years=2016:2017;
        otherwise
            ['need to specify a run code']
    end
    
    
    %%
    for year = years
        
        [year]
        tic
        cd(input_dir)
        fname = ['NOAA_NSIDC_DAILY_MAPPED_TO_LLC' num2str(llcN) '_' num2str(year)]
        load(fname)
        
        % only use the goddard-mergeded field.
        sea_ice_conc = NOAA_NSIDC_DAILY_MAPPED_TO_LLC.goddard_merged./100;
        
        % any sea ice value less than 0 should be masked
        sea_ice_conc(find(sea_ice_conc < 0))=-9999;
        
        cd(output_dir)
        fid = fopen([fname '.bin'],'w','b');
        fwrite(fid, sea_ice_conc,'single');
        fclose(fid);
        toc
        
        clear sea_ice_conc;
    end
end


%%

['zero years']
for year = zero_years
    nd=datenum([year 12 31 0 0 0])-datenum([year 1 1 0 0 0])+1
    [year nd]
    
    fname = ['NOAA_NSIDC_DAILY_MAPPED_TO_LLC' num2str(llcN) '_' num2str(year) '.bin']
    fid = fopen(fname,'w','b');
    f = zeros(1, [llcN*llcN*13*nd]);
    fwrite(fid, f - 9999, 'single');
    fclose(fid);
end
['finished zero years']
