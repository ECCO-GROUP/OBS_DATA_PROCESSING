function [MITprof, F_grid_PF_XYZ_to_INDEX] = ...
    get_profpoint_llc_ian(lon_llc, lat_llc, mask_llc, ...
                         nx, MITprof,F_grid_PF_XYZ_to_INDEX)
%
% This script finds the 'prof_point' for the mitgcm profile package for a
% llc grid.
%
% Filename; get_profpoint_llc_ian
%  ** former filename : updates get_profpoint (but is a complete rewrite of
%  that routine)
% Date Created: 2014-07-30
% Last Modified: 2014-08-05
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
%   lon_llc, lat_llc  : the XC and YC of the llc grid in compact format
%   mask_llc          : a mask with 1 or 0 denoting whether to use a point
%                       or not in the search
%   nx  :             : the llcNX
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

ny=nx*13;

if nargin < 6
    lon_pf      = patchface3D(nx,ny,1,lon_llc,2);
    lat_pf      = patchface3D(nx,ny,1,lat_llc,2);
    mask_llc_pf = patchface3D(nx,ny,1,mask_llc,2);
    %
    [X_grid, Y_grid, Z_grid] = sph2cart(lon_llc*deg2rad, lat_llc*deg2rad, 1);
    
    % convert X,Y,Z coords to global view
    X_grid_pf =patchface3D(nx,ny,1,X_grid,2);
    Y_grid_pf =patchface3D(nx,ny,1,Y_grid,2);
    Z_grid_pf =patchface3D(nx,ny,1,Z_grid,2);
    
    AI_grid_pf = X_grid_pf.*0;
    AI_grid_pf(1:end) = 1:length(AI_grid_pf(:));
    
    % the count of the every grid point in the patchface universe
    % these AI points are in no domain, they are just a list of points
    % a mapping between these x,y,z points and the AI_grid_pf index
    
    % exclude points that are nan - where do these nan points come from? no
    % one knows.
    good_ins = find(mask_llc_pf == 1);
    
    % make a subset of X,Y,Z and AI to include only the non-nan, not masked points
    X_grid_pf = X_grid_pf(good_ins);
    Y_grid_pf = Y_grid_pf(good_ins);
    Z_grid_pf = Z_grid_pf(good_ins);
    AI_grid_pf = AI_grid_pf(good_ins);
    
    %%
    % these are the x,y,z coordinates of the 'good' cells in
    xyz= [X_grid_pf Y_grid_pf Z_grid_pf];
    
    F_grid_PF_XYZ_to_INDEX = scatteredInterpolant(xyz, AI_grid_pf,'nearest');
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