% This script loads the important variables of the llc270 grid.

% Filename; load_llc270_grid.m
%  ** former filename : 
% Date Created: unknown
% Last Modified: 2014-08-06
% 
% notes: added 'bad_ins_270' which shows all the points that have 0 in lat
% and lon and bathy. those seem to be undefined tiles.

if length(strfind(hostname,'skylla2'))
    rootdir = '/net/skylla/data17/IGF_output'
elseif length(strfind(hostname,'skylla'))
    rootdir = '/data17/IGF_output'
elseif length(strfind(hostname,'pfe'))
    rootdir = '/nobackup/ifenty';
elseif length(strfind(hostname,'bridge'))
    rootdir = '/nobackup/ifenty';
elseif length(strfind(hostname,'cnidarian'))
    rootdir = '/home/ifenty/data/grids';
elseif length(strfind(hostname,'LMC-054074'))
    rootdir = '/Users/ifenty/My Projects/grids';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load llc270 grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
llcN = 270;

llc270_grid_dir = [rootdir '/grid_llc270']
bathy_270_fname = 'bathy_llc270'

hf_270_fname = 'hFacC.data';

%% load grid
cd(llc270_grid_dir)
%tic;bathy_270 = readbin(bathy_270_fname,[llcN 13*llcN 1],1,'real*8',0,'ieee-be');toc;
tic;bathy_270 = readbin(bathy_270_fname,[llcN 13*llcN 1],1,'real*4',0,'ieee-be');toc;

%lon_270 = readbin('XC.data',[llcN 13*llcN 1],1,'real*4',0,'ieee-be');
%lat_270 = readbin('YC.data',[llcN 13*llcN 1],1,'real*4',0,'ieee-be');

[lat_270, lon_270] = load_grid_fields_from_tile_files(pwd, 270);

XG_270 = readbin('XG.data',[llcN 13*llcN 1],1,'real*4',0,'ieee-be');
YG_270 = readbin('YG.data',[llcN 13*llcN 1],1,'real*4',0,'ieee-be');

RAC_270 = readbin('RAC.data',[llcN 13*llcN 1],1,'real*4',0,'ieee-be');

DXG_270 = readbin('DXG.data',[llcN 13*llcN 1],1,'real*4',0,'ieee-be');
DYG_270 = readbin('DYG.data',[llcN 13*llcN 1],1,'real*4',0,'ieee-be');

DXC_270 = readbin('DXC.data',[llcN 13*llcN 1],1,'real*4',0,'ieee-be');
DYC_270 = readbin('DYC.data',[llcN 13*llcN 1],1,'real*4',0,'ieee-be');

Depth_270 = readbin('Depth.data',[llcN 13*llcN 1],1,'real*4',0,'ieee-be');

hFacC_270 = readbin('hFacC.data',[llcN 13*llcN 50],1,'real*4',0,'ieee-be');
hFacW_270 = readbin('hFacW.data',[llcN 13*llcN 50],1,'real*4',0,'ieee-be');
hFacS_270 = readbin('hFacS.data',[llcN 13*llcN 50],1,'real*4',0,'ieee-be');

%%
basin_mask_270 = readbin('basin_masks_eccollc_llc270.bin',...
    [llcN 13*llcN 1],1,'real*4',0,'ieee-be');
deg2rad = pi/180.0;

%% this is how we find bogus points in the compact grid.
bad_ins_270 = find(lat_270 ==0 & lon_270 == 0 & bathy_270 == 0);
AI_270 = 1:length(lon_270(:));
AI_270 = reshape(AI_270, size(lon_270));

%%
good_ins_270 = setdiff(AI_270(:)', bad_ins_270(:));
%%
[X_270, Y_270, Z_270] = sph2cart(lon_270*deg2rad, lat_270*deg2rad, 1);

X_270(bad_ins_270)= NaN;
Y_270(bad_ins_270)= NaN;
Z_270(bad_ins_270)= NaN;
AI_270(bad_ins_270)= NaN;

lon_270(bad_ins_270) = NaN;
lat_270(bad_ins_270) = NaN;

%wet_ins_270 = find(hFacC_270 > 0);
%dry_ins_270 = find(hFacC_270 ==0);

for k = 1:50
    tmp = hFacC_270(:,:,k);
    dry_ins_270_k{k} = find(tmp == 0);
    wet_ins_270_k{k} = find(tmp > 0);
end

hf0_270 = hFacC_270(:,:,1);

tmp = bathy_270(wet_ins_270_k{1});
max(tmp)
min(tmp)

blank_270 = bathy_270.*NaN;

hf0_270_pf = patchface3D(llcN, llcN*13, 1, hFacC_270(:,:,1), 2);
bathy_270_pf = patchface3D(llcN, llcN*13, 1, bathy_270, 2);

RAC_270_pf = patchface3D(llcN, llcN*13, 1, RAC_270 , 2);

[delR_270, z_top_270, z_bot_270, z_cen_270] = make_llc270_cell_centers;

['finished load llc270 grid']
