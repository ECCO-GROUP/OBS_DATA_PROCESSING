function [point_id, F_grid_PF_XYZ_to_INDEX] = ...
    get_nearest_point_llc_ian(lon_llc, lat_llc, mask_llc, ...
                         nx, point_lon, point_lat, F_grid_PF_XYZ_to_INDEX)
%
% This script finds the 'prof_point' for any arbitrary lat/lon points in an
% llc grid.
%
% Filename; get_nearest_point_llc_ian
% Date Created: 2014-09-23
% Last Modified: 
%
% notes:
%    previous version seems to have been messed up for high northern
%    latitude points
%
%    2014-08-03 : allows for only wet point to be selected via the addition
%    of a mask_llc field
%
%     ... 08-05 : removed the saving of F, seems dangerous to reuse this
%                 file and changed to it to specifying F if you like.
%     ... 09-5

% usage: 
% function prof_point = get_profpoint_llc_ian(grid_lon, grid_lat, mask_llc,
%                       nx, MITprof)
%
%   lon_llc, lat_llc  : the XC and YC of the llc grid in compact format
%   mask_llc          : a mask with 1 or 0 denoting whether to use a point
%                       or not in the search
%   nx  :             : the llcNX
%   point_lon, point_lat : the lon/lat of the profile coordinates
%
%   F_grid_PF_XYZ_to_INDEX : the mapping scatterInterp structure. returned
%                    in case you want to use it again in a loop
%
% returns the index of the nearest neighbor between
% the point_lon, and point_lat and the lon_llc, and lat_llc;

deg2rad = pi/180.0;

ny=nx*13;

% if nargin < 7 that means we do not have the mapping structure and we need
% to create it.
if nargin < 7

    % convert the lon,lat points to x,y,z on the sphere
    [X_grid, Y_grid, Z_grid] = ...
        sph2cart(lon_llc*deg2rad, lat_llc*deg2rad, 1);
    
    % define an array with each element being the index of the domain
    AI_grid = 1:length(X_grid(:));
    
    % exclude points that are nan in the mask
    % why are some lat/lon points NaN?
    good_ins = find(mask_llc == 1);
    
    X_grid = X_grid(good_ins);
    Y_grid = Y_grid(good_ins);
    Z_grid = Z_grid(good_ins);
    AI_grid = AI_grid(good_ins)';
    
    % create a mapping that gives you the nearest neighbor point in the 
    % domain given an x,y,z coordinate

    F_grid_PF_XYZ_to_INDEX = ...
        scatteredInterpolant([X_grid Y_grid Z_grid], AI_grid,'nearest');
end

% at this point we have the mapping structure.

% convert the lon/lat coordinates of the points provided in the arguments
% to x,y,z on the sphere.
[point_x, point_y, point_z] = ...
    sph2cart(point_lon*deg2rad, point_lat*deg2rad, 1);

% use the mapping structure to find the index of the nearest neighbor to
% those points
point_id = F_grid_PF_XYZ_to_INDEX(point_x, point_y, point_z);

% sanity check to see if the point ids returned make sense
residual = point_id - floor(point_id);

if length(find(residual ~= 0)) > 0
    point_id = NaN;
    ['could not find unique nearest neighbors']
else
    ['all profile points matched successfully']
end

%%
%% TESTING
if 1 ==0
        
    [pid, F] = get_nearest_point_llc_ian(lon_90, lat_90, hf0_90, 90, -117, 1);
    [lon_90(pid) lat_90(pid)]
end
