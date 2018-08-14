function [MITprof, F_grid_PF_XYZ_to_INDEX,residual] = ...
    get_profpoint_latlon_ian(lon, lat, mask, ...
                         MITprof,F_grid_PF_XYZ_to_INDEX)
%
% This script finds the equivalent 'prof_point' for the mitgcm profile package for a
% latlon grid.
%
% Filename; get_profpoint_latlon_ian
%  ** former filename :get_profpoint_llc_ian + 
%     updates get_profpoint (but is a complete rewrite of
%  that routine)
% Date Created: 2015-02-16
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
%

% usage: 
% function prof_point = get_profpoint_llc_ian(grid_lon, grid_lat, mask_llc,
%                       nx, MITprof)
%
%   lon, lat  : the XC and YC of the latlon grid in 
%   mask          : a mask with 1 or 0 denoting whether to use a point
%                       or not in the search

%   MITprof           : the MITprofile package structure which
%                       contains the prof_lon, prof_lat the lon/lat
%                       profile coordinates
%   F_grid_PF_XYZ_to_INDEX : the mapping scatterInterp structure. returned
%                    in case you want to use it again in a loop
%
% returns the MITprof structure with an updated prof_point field
% indicating the index of the nearest neighbor between
% the MITprof.prof_lon,prof_lat and grid_lon/grid_lat

deg2rad = pi/180.0;

if nargin < 5
    
    [X_grid, Y_grid, Z_grid] = sph2cart(lon*deg2rad, lat*deg2rad, 1);
    
    
    AI_grid = X_grid.*0;
    AI_grid(1:end) = 1:length(AI_grid(:));
    
    % the count of the every grid point in the patchface universe
    % these AI points are in no domain, they are just a list of points
    % a mapping between these x,y,z points and the AI_grid_pf index
    
    % exclude points that are nan - where do these nan points come from? no
    % one knows.
    good_ins = find(mask == 1);
    
    length(good_ins)
    length(AI_grid(:))
    
    figure(66);clf;imagesc(AI_grid);colorbar
    
    % make a subset of X,Y,Z and AI to include only the non-nan, not masked points
    X_grid = X_grid(good_ins);
    Y_grid = Y_grid(good_ins);
    Z_grid = Z_grid(good_ins);
    AI_grid = AI_grid(good_ins);
    
    %%
    % these are the x,y,z coordinates of the 'good' cells in
    xyz= [X_grid Y_grid Z_grid];
    
    F_grid_PF_XYZ_to_INDEX = scatteredInterpolant(xyz, AI_grid,'nearest');
end


point_lon = MITprof.prof_lon;
point_lat = MITprof.prof_lat;

[prof_x, prof_y, prof_z] = ...
    sph2cart(point_lon*deg2rad, point_lat*deg2rad, 1);

MITprof.prof_point = F_grid_PF_XYZ_to_INDEX(prof_x, prof_y, prof_z);

residual = MITprof.prof_point - floor(MITprof.prof_point);

if length(find(residual ~= 0)) > 0
    MITprof = NaN;
    ['could not find unique nearest neighbors']
else
    ['all profile points matched successfully']
end