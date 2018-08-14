% This script makes the mapping between all "good" points in the llc90 grid to the
% 'all index' matrix -- note this is only for the surface points
% note the good points are defined as being defined in XC YC

% Filename; make_F_llc90_ALL_INS_SURF_XYZ_to_INDEX.m
%  ** former filename : 
% Date Created: 2014-08-05
% Last Modified: 

% notes:

set(0,'DefaultTextInterpreter','none');
%%
hostname = getenv('HOSTNAME')

if length(strfind(hostname,'skylla2'))
    rootdir = '/net/skylla/data17/IGF_output'
elseif length(strfind(hostname,'skylla'))
    rootdir = '/data17/IGF_output'
elseif length(strfind(hostname,'pfe'))
    rootdir = '/nobackup/ifenty';
elseif length(strfind(hostname,'bridge'))
    rootdir = '/nobackup/ifenty';
elseif length(strfind(hostname,'cnidarian'))
    rootdir = '/home/ifenty/data/';
end

%%---------------------------------------------
%%  BEGIN THE PROGRAM
%%---------------------------------------------

%%---------------------------------------------
%%  Init profile package
%%---------------------------------------------
MITprof_path;

%%---------------------------------------------
%%  Read in grid for sigma and climatology
%%---------------------------------------------
load_llc90_grid

%%---------------------------------------------
%% Prepare the nearest neighbor mapping
%%---------------------------------------------
cd(llc90_grid_dir)

% do the big mapping.
AI_90 = blank_90;
AI_90(1:end) = 1:length(bathy_90(:));

% all index has been shortened somewhat to include 'wet' (not nan)
% values -- llc90 grid isn't defined in center of continents!
AI_90 = AI_90(good_ins_90);

%%
% these are the x,y,z coordinates of the 'good' cells in llc90
xyz= [X_90(good_ins_90) Y_90(good_ins_90) Z_90(good_ins_90)];

% these AI points are in no domain, they are just a list of points
% a mapping between these x,y,z points and the AI_90 index
% you provide an x,y,z and F_llc90_XYZ_to_INDEX provides you with an
% index in the global llc90 file correpsonding with the nearest
% neighbor.
F_llc90_ALL_INS_SURF_XYZ_to_INDEX = TriScatteredInterp(xyz, AI_90,'nearest');

save('F_llc90_ALL_INS_SURF_XYZ_to_INDEX.mat','F_llc90_ALL_INS_SURF_XYZ_to_INDEX');

['new F saved']

