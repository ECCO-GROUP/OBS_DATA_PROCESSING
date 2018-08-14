function [MITprof] = ...
    update_gci_on_prepared_profiles(grid_lons, grid_lats, MITprof)


% This script updates ECCO profile.nc files with their 'grid cell index'
% - an index of which grid cell in the domain that the profile is 
% closest to. 

% Filename; update_gci_on_prepared_profiles.m
%  ** former filename :update_climatologyTS_on_prepared_profiles.m
% Date Created: 2018-08-14
% Last Modified: 2018-08-14
%

% notes:
                    
% usage
%    grid_lons, grid_lats    : matrices with the lons and lats of the grid
%    MITprof                 : The MITprof that will get the 'prof_gci' field 


%%---------------------------------------------
%%  BEGIN THE PROGRAM
%%---------------------------------------------

deg2rad = pi/180.0;

['mapping model grid lat/lon coordinates to cartesian (x,y,z) coordinates']
[grid_x, grid_y, grid_z] =  ...
     sph2cart(grid_lons*deg2rad, grid_lats*deg2rad, 1);
    
    
all_indices = 1:length(grid_x(:));
all_indices = all_indices(:);
% these are the x,y,z coordinates of the 'good' cells in llc90
xyz= [grid_x(:), grid_y(:), grid_z(:)];

size(xyz)
size(all_indices)
% the'tri scattered interpolator' will map the grid cell index 
% values (all_indices) to x,y,z points that you pass it.
F = TriScatteredInterp(xyz, all_indices(:),'nearest');

% cov
['mapping profile lat/lon coordinates to cartesian (x,y,z) coordinates']
[prof_x, prof_y, prof_z] = ...
    sph2cart(MITprof.prof_lon*deg2rad, MITprof.prof_lat*deg2rad, 1);
        
% map a llc90 grid index to each profile.
prof_grid_cell_index = F(prof_x, prof_y, prof_z);
       
% finally add this prof_gci field to the MITprof structure 
MITprof.prof_gci = prof_grid_cell_index;

