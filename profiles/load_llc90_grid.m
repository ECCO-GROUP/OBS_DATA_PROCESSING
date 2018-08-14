% This script loads the important variables of the llc90 grid.

% Filename; load_llc90_grid.m
%  ** former filename : 
% Date Created: unknown
% Last Modified: 2014-06-14

%hostname = getenv('HOSTNAME')
clear hostname;
if length(strfind(hostname,'skylla2'))
    grootdir = '/net/skylla/data17/IGF_output'
elseif length(strfind(hostname,'skylla'))
    grootdir = '/data17/IGF_output'
elseif length(strfind(hostname,'pfe'))
    grootdir = '/nobackup/ifenty';
elseif length(strfind(hostname,'bridge'))
    grootdir = '/nobackup/ifenty';
elseif length(strfind(hostname,'cnidarian'))
    grootdir = '/home/ifenty/data/grids';
elseif length(strfind(hostname,'LMC-054074'))
    grootdir = '/Users/ifenty/My Projects/grids';
else
    grootdir = [homedir 'data/grids']
end

%% BATHY
% CODES 0 is ECCOv4 
%       1 is ice shelf cavity
BATHY_CODE=0

switch BATHY_CODE
    case 0
        bathy_dir = [grootdir '/grid_llc90']
        bathy_90_fname = 'bathy_eccollc_90x50_min2pts.bin';

    case 1
        bathy_dir = [grootdir '/grid_llc90/bathy_ice_version_8_no_single_cavities_no_adjacent_bathy_shallower']
        bathy_90_fname = 'BATHY_ICE_SHELF_CAVITY_PLUS_ICE_FRONT_LLC_0090.bin'

end
           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load llc90 grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

llcN = 90;
llc90_grid_dir =[grootdir '/grid_llc90']

cd(bathy_dir);
tic;bathy_90 = readbin(bathy_90_fname,[llcN 13*llcN 1],1,'real*4',0,'ieee-be');toc;

cd([llc90_grid_dir '/no_blank'])

% clear XC YC DXF DYF RAC XG YG DXV DYU RAZ DXC DYC RAW RAS DXG DYG;

lon_90 = readbin('XC.data',[llcN 13*llcN 1],1,'real*4',0,'ieee-be');
lat_90 = readbin('YC.data',[llcN 13*llcN 1],1,'real*4',0,'ieee-be');

XG_90 = readbin('XG.data',[llcN 13*llcN 1],1,'real*4',0,'ieee-be');
YG_90 = readbin('YG.data',[llcN 13*llcN 1],1,'real*4',0,'ieee-be');

RAC_90 = readbin('RAC.data',[llcN 13*llcN 1],1,'real*4',0,'ieee-be');
RAC_90_pf = patchface3D(llcN, llcN*13, 1, RAC_90, 2);

RC_90 = readbin('RC.data',[1 50],1,'real*4',0,'ieee-be');
RF_90 = readbin('RF.data',[1 50],1,'real*4',0,'ieee-be');

DXG_90 = readbin('DXG.data',[llcN 13*llcN 1],1,'real*4',0,'ieee-be');
DYG_90 = readbin('DYG.data',[llcN 13*llcN 1],1,'real*4',0,'ieee-be');

DXC_90 = readbin('DXC.data',[llcN 13*llcN 1],1,'real*4',0,'ieee-be');
DYC_90 = readbin('DYC.data',[llcN 13*llcN 1],1,'real*4',0,'ieee-be');


AI_90 = 1:length(lon_90(:));
AI_90 = reshape(AI_90, size(lon_90));

blank_90 = bathy_90.*NaN;
bathy_90_pf = patchface3D(llcN, llcN*13, 1, bathy_90, 2);

[delR_90, z_top_90, z_bot_90, z_cen_90] = make_llc90_cell_centers;
deg2rad = pi/180.0;

[X_90, Y_90, Z_90] = sph2cart(lon_90*deg2rad, lat_90*deg2rad, 1);

%%
%% no longer up to date.
%hf0_90_pf   = patchface3D(llcN, llcN*13, 1, hFacC_90(:,:,1), 2);

cd([llc90_grid_dir])
Depth_90 = readbin('Depth.data',[llcN 13*llcN 1],1,'real*4',0,'ieee-be');
hFacC_90 = readbin('hFacC.data',[llcN 13*llcN 50],1,'real*4',0,'ieee-be');

wet_ins_90 = find(hFacC_90 > 0);
dry_ins_90 = find(hFacC_90 == 0);


for k = 1:50
   tmp = hFacC_90(:,:,k);
   dry_ins_90_k{k} = find(tmp == 0);
   wet_ins_90_k{k} = find(tmp > 0);
   nan_ins_90_k{k} = find(isnan(tmp));
end

hf0_90 = hFacC_90(:,:,1);