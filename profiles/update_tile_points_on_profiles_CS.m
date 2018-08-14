% This script updates the prof_points and all tile interpolation point
% fields and sends the updated profile files to be written in netcdf format
% works for Cube Spehre

% Filename; update_prof_points_on_profiles_CS.m
%  ** former filename : update_prof_points_on_profiles.m
% Date Created: 2015-02-12
% Last Modified: 
% 


% notes:
%    this script exists to allow one to use the same MITprofile package
%    netcdf files in different llc configurations easily.
%
%    2014-08-06 : added some error checking and the ability to specify that
%                 you want to interpolate to only wet points automatically, 
%                 if desired.


clear all;
set(0,'DefaultTextInterpreter','none');
%%

make_root_dir;

%%---------------------------------------------
%% SET INPUT PARAMETERS
%%---------------------------------------------

% profile input dir and filenames
%input_dir  = [rootdir '/observations/insitu/ECCO_v4r2/20140806']

% output and figure directories
%output_dir  = [rootdir '/observations/insitu/ECCO_v4r2/20140806b']

input_dir   = [rootdir '/observations/insitu/ECCO_v5r1/20150209/SIG_PP_CLIM']

output_dir  = [rootdir '/observations/insitu/CS/20150212']
mkdir(output_dir)

clear fDataIn*;

ilists = [24:33]

file_suffix_in  = '_SIG_PP_CLIM_llc270_20150209.nc'
file_suffix_out = '_cs510_20150212.nc'

% 1 - 20 area ARGO
yrs = 1995:2014
for j = 1:length(yrs)
    fDataBase{j} = ['argo_' num2str(yrs(j))];
end

fDataBase{21} =  ['ctd'     ];
fDataBase{22}  = ['climode' ];
fDataBase{23}  = ['itp'     ];
fDataBase{24}  = ['seals'   ];

fDataBase{25}  = ['beaufortgyremooring'   ];
fDataBase{26}  = ['beringstraitmooring'   ];

fDataBase{27} = ['ctd_Arctic_NordicSeas' ];
fDataBase{28} = ['ctdhilat_nodupices'    ];
fDataBase{29} = ['ctdlowlat' ];

fDataBase{30} = ['framstraitmooring' ];

fDataBase{31} = ['ices19922012hi_pot_theta' ];
fDataBase{32} = ['ices19922012lo_pot_theta'  ];

fDataBase{33} = ['xbt'  ];


clear fDataIn*
for i = 1:length(fDataBase)
    fDataIn{i}  = [fDataBase{i} file_suffix_in]
    fDataOut{i} = [fDataBase{i} file_suffix_out]
end


%%---------------------------------------------
%%  BEGIN THE PROGRAM
%%---------------------------------------------

%%---------------------------------------------
%%  Init profile package
%%---------------------------------------------
MITprof_path;


%%---------------------------------------------
%%  Read in grid for climatology
%%---------------------------------------------

nz=50;
load_cube86_grid

lon_CS;
lat_CS ;

mask_CS = lon_CS.*0+1;
   

%%---------------------------------------------
%% Read and process the profile files
%%---------------------------------------------
% F is the mapping -- should probably remake this every time.
clear F;

%%
for ilist = ilists
    
    fDataIn{ilist}
    cd(input_dir)
    
    clear fileOut MITprof
    
    MITprof = MITprof_read(fDataIn{ilist});
    
    if exist('F')
        [MITprof, F] = ...
            get_profpoint_cs_ian(lon_CS, lat_CS, ...
            mask_CS,  MITprof,F);
    else
        [MITprof, F] = ...
            get_profpoint_cs_ian(lon_CS, lat_CS, ...
            mask_CS,  MITprof);
    end
    
    
    
    %%
    figure(ilist);clf;hold on;
    plot(MITprof.prof_lon, MITprof.prof_lat,'r.')
    plot(MITprof.prof_interp_lon, MITprof.prof_interp_lat,'bo')
    
    figure;
    projection_code=1;
    plot_proj=0;
    makemmap_general;
    m_plot(MITprof.prof_lon, MITprof.prof_lat,'r.')
    m_plot(MITprof.prof_interp_lon, MITprof.prof_interp_lat,'bo')
    makemmap_general;
    
    %% check max distances
    d = vdist(MITprof.prof_lat, MITprof.prof_lon, ...
          MITprof.prof_interp_lat, MITprof.prof_interp_lon);
    
    d=d./1e3;  % distance now in KM.
    
    dx=sqrt(RAC_CS(MITprof.prof_point))./1e3;
    ins = find(d./dx > 1);
    
    figure
    plot(MITprof.prof_lon(ins), MITprof.prof_lat(ins),'g.','MarkerSize',10)
    title(fDataIn{ilist})
    
    
    ['writing ' fDataOut{ilist}]
    fileOut=[output_dir '/'  fDataOut{ilist}]
    fprintf('%s\n',fileOut);
    write_profile_structure_to_netcdf;
end


