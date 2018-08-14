function [] = write_profile_structure_to_netcdf(MITprof, fileOut)

% This script writes ECCO profile.nc files in netcdf format

% Filename; write_profile_structure_to_netcdf.m
%  ** former filename : step05_write_nc.m
% Date Created: 2014-07-01
% Last Modified: 2016-03-08


% notes:
%   hacked version of An's step05a_write_nc.m
%
%   --- from An'
%   07-Jun-2013: augmented info related to tiles. use this from now on instead of .m
%   after ~/matlab_gcmfaces/MITprof/profiles_IO_v2/MITprof_create.m
%   addpath('/data4/atnguyen/matlab_gcmfaces/MITprof/profiles_misc');
%
%   03-Mar-2015 : If the interpolating fields aren't present, don't demand
%   that they are there before writing the MITprof structure to netcdf

%   29-Apr 2015 : cleaned up output

%   14-May 2015 : added bin ids
%   15-Sep 2015 : added a field called prof_area_gamma that reflects any
%                 scalar factor added to the weights that would need to
%                 be removed for apples-to-apples comparison with unscaled
%                 weight runs.
%   08-Mar 2016 : debugging lines


flagT=isfield(MITprof,'prof_T');
flagS=isfield(MITprof,'prof_S');
flagU=isfield(MITprof,'prof_U');
flagV=isfield(MITprof,'prof_V');
flagPT=isfield(MITprof,'prof_ptr');
flagSSH=isfield(MITprof,'prof_ssh');
flagBP=isfield(MITprof,'prof_bp');



%%  IAN ADDED
% area gamma
%['flag Gamma']
flagGAMMA=isfield(MITprof,'prof_area_gamma');

listVarsGAMMA = {};
if flagGAMMA;listVarsGAMMA = ['prof_area_gamma'];end

% prof_gci
flagGridIndex = isfield(MITprof,'prof_gci');

listVarsGCI = {};
if flagGridIndex;listVarsGCI = ['prof_gci'];end

% T/S clim
flagTClim = isfield(MITprof,'prof_Tclim');
flagSClim = isfield(MITprof,'prof_Sclim');

listVarsClim = {};
if flagTClim;listVarsClim = [listVarsClim 'prof_Tclim'];end
if flagSClim;listVarsClim = [listVarsClim 'prof_Sclim'];end

% model mean
flagModTMean = isfield(MITprof,'prof_Tmodel_mean');
flagModSMean = isfield(MITprof,'prof_Smodel_mean');

listVarsModMean = {};
if flagModTMean;listVarsModMean = [listVarsModMean 'prof_Tmodel_mean'];end
if flagModSMean;listVarsModMean = [listVarsModMean 'prof_Smodel_mean'];end

% bin id
flagBinID_a = isfield(MITprof,'prof_bin_id_a');
flagBinID_b = isfield(MITprof,'prof_bin_id_b');
flagBinID_c = isfield(MITprof,'prof_bin_id_c');
flagBinID_d = isfield(MITprof,'prof_bin_id_d');

listVarsBinID = {};
if flagBinID_a;listVarsBinID = [listVarsBinID 'prof_bin_id_a'];end;
if flagBinID_b;listVarsBinID = [listVarsBinID 'prof_bin_id_b'];end;
if flagBinID_c;listVarsBinID = [listVarsBinID 'prof_bin_id_c'];end;
if flagBinID_d;listVarsBinID = [listVarsBinID 'prof_bin_id_d'];end;

%% END IAN ADD

%%
nLev=length(MITprof.prof_depth);
MITprof.prof_depth=reshape(MITprof.prof_depth,length(MITprof.prof_depth),1);
nProf=length(MITprof.prof_point);

iPROF=nProf; iDEPTH=nLev; 
% the number of interpolation points.  this might not always be one
% (2016-03-08)
%iINTERP=1;

% this makes it depend on the profile itself.  if we have more than one
% interpolation point, then iINTERP will be > 1
iINTERP=size(MITprof.prof_interp_i,2)

lTXT=size(MITprof.prof_descr,2); fillval=double(-9999);

%06.Mar.2014: check Tweight and Sweight here, make sure it's zero and not -9999
if(flagT==1);
  MITprof.prof_Tflag=zeros(size(MITprof.prof_T));
  clear ii;ii=find(isnan(MITprof.prof_Tweight)==1);if(length(ii)>0);MITprof.prof_Tweight(ii)=0;end;
  clear ii;ii=find(MITprof.prof_Tweight(:)<0);if(length(ii)>0);MITprof.prof_Tweight(ii)=0;end;
end;
if(flagS==1);
  MITprof.prof_Sflag=zeros(size(MITprof.prof_S));
  clear ii;ii=find(isnan(MITprof.prof_Sweight)==1);if(length(ii)>0);MITprof.prof_Sweight(ii)=0;end;
  clear ii;ii=find(MITprof.prof_Sweight(:)<0);if(length(ii)>0);MITprof.prof_Sweight(ii)=0;end;
end;

%=============list of variables that will actually be in the file==============%
list_vars={'prof_T','prof_Tweight','prof_Testim','prof_Terr','prof_Tflag'};
if(flagS==1);
    list_vars=[list_vars {'prof_S','prof_Sweight','prof_Sestim',...
        'prof_Serr','prof_Sflag'}];
end;

list_vars_plus=[{'prof_depth','prof_descr','prof_date',...
    'prof_YYYYMMDD','prof_HHMMSS',...
    'prof_lon','prof_lat','prof_basin','prof_point'}, ...
    list_vars listVarsClim listVarsModMean listVarsBinID listVarsGCI listVarsGAMMA];

%profiles_process_main_v2/profiles_prep_mygrid.m
list1d={'prof_interp_XC11','prof_interp_YC11','prof_interp_XCNINJ',...
    'prof_interp_YCNINJ'};
list2d={'prof_interp_i','prof_interp_j','prof_interp_lon',...
    'prof_interp_lat','prof_interp_weights'};
listAll={list1d{:},list2d{:}};

%==========masters table of variables, units, names and dimensions=============%

mt_v={'prof_depth'}; mt_u={'me'}; mt_n={'depth'}; mt_d={'iDEPTH'};

mt_v=[mt_v 'prof_date']; mt_u=[mt_u ' ']; mt_n=[mt_n 'Julian day since Jan-1-0000']; mt_d=[mt_d 'iPROF'];
mt_v=[mt_v 'prof_YYYYMMDD']; mt_u=[mt_u ' ']; mt_n=[mt_n 'year (4 digits), month (2 digits), day (2 digits)']; 
mt_d=[mt_d 'iPROF'];
mt_v=[mt_v 'prof_HHMMSS']; mt_u=[mt_u ' ']; mt_n=[mt_n 'hour (2 digits), minute (2 digits), second (2 digits)']; 
mt_d=[mt_d 'iPROF'];
mt_v=[mt_v 'prof_lon']; mt_u=[mt_u ' ']; mt_n=[mt_n 'Longitude (degree East)']; mt_d=[mt_d 'iPROF'];
mt_v=[mt_v 'prof_lat']; mt_u=[mt_u ' ']; mt_n=[mt_n 'Latitude (degree North)']; mt_d=[mt_d 'iPROF'];
mt_v=[mt_v 'prof_basin']; mt_u=[mt_u ' ']; mt_n=[mt_n 'ocean basin index (ecco 4g)']; mt_d=[mt_d 'iPROF'];
mt_v=[mt_v 'prof_point']; mt_u=[mt_u ' ']; mt_n=[mt_n 'grid point index (ecco 4g)']; mt_d=[mt_d 'iPROF'];


%% IAN ADD
if(flagGAMMA ==  1);
    ['flag gamma'];
    mt_v=[mt_v 'prof_area_gamma']; 
    mt_u=[mt_u ' '];
    mt_n=[mt_n 'scaling factor (real number) applied to the T and S weights']; 
    mt_d=[mt_d 'iPROF'];
end

if(flagGridIndex == 1);
    ['flag grid index'];
    mt_v=[mt_v 'prof_gci']; 
    mt_u=[mt_u ' '];
    mt_n=[mt_n 'index (integer) of the nearest model grid cell cell to the profile (index from 1:n)']; 
    mt_d=[mt_d 'iPROF'];
end

if flagTClim
    mt_v = [mt_v 'prof_Tclim'];
    mt_u=[mt_u 'degree C']; mt_n=[mt_n 'potential temperature']; mt_d=[mt_d 'iPROF,iDEPTH'];
end
if flagSClim
    mt_v = [mt_v 'prof_Sclim'];
    mt_u=[mt_u 'S']; mt_n=[mt_n 'salt fool']; mt_d=[mt_d 'iPROF,iDEPTH'];
end
if flagModTMean
    mt_v = [mt_v 'prof_Tmodel_mean'];
    mt_u=[mt_u 'degree C']; mt_n=[mt_n 'potential temperature']; mt_d=[mt_d 'iPROF,iDEPTH'];
end
if flagModSMean
    mt_v = [mt_v 'prof_Smodel_mean'];
    mt_u=[mt_u 'S']; mt_n=[mt_n 'salt fool']; mt_d=[mt_d 'iPROF,iDEPTH'];
end

if(flagBinID_a == 1);
    mt_v=[mt_v 'prof_bin_id_a'];
    mt_u=[mt_u ' '];
    mt_n=[mt_n 'bin index (integer) A']; 
    mt_d=[mt_d 'iPROF'];
end
if(flagBinID_b == 1);
    mt_v=[mt_v 'prof_bin_id_b'] ;
    mt_u=[mt_u ' '];
    mt_n=[mt_n 'bin index (integer) B']; 
    mt_d=[mt_d 'iPROF'];
end
if(flagBinID_c == 1);
    mt_v=[mt_v 'prof_bin_id_c'] ;
    mt_u=[mt_u ' '];
    mt_n=[mt_n 'bin index (integer) C']; 
    mt_d=[mt_d 'iPROF'];
end
if(flagBinID_d == 1);
    mt_v=[mt_v 'prof_bin_id_d'] ;
    mt_u=[mt_u ' '];
    mt_n=[mt_n 'bin index (integer) D']; 
    mt_d=[mt_d 'iPROF'];
end

%% END IAN ADD

%
mt_v=[mt_v 'prof_T']; mt_u=[mt_u 'degree C']; mt_n=[mt_n 'potential temperature']; mt_d=[mt_d 'iPROF,iDEPTH'];
mt_v=[mt_v 'prof_Tweight']; mt_u=[mt_u '(degree C)^-2']; mt_n=[mt_n 'pot. temp. least-square weight']; mt_d=[mt_d 'iPROF,iDEPTH'];
mt_v=[mt_v 'prof_Testim']; mt_u=[mt_u 'degree C']; mt_n=[mt_n 'pot. temp. estimate (e.g. from atlas)']; mt_d=[mt_d 'iPROF,iDEPTH'];
mt_v=[mt_v 'prof_Terr']; mt_u=[mt_u 'degree C']; mt_n=[mt_n 'pot. temp. instrumental error']; mt_d=[mt_d 'iPROF,iDEPTH'];
mt_v=[mt_v 'prof_Tflag']; mt_u=[mt_u ' ']; mt_n=[mt_n 'flag = i > 0 means test i rejected data.']; mt_d=[mt_d 'iPROF,iDEPTH'];
%
if(flagS==1);
mt_v=[mt_v 'prof_S']; mt_u=[mt_u 'psu']; mt_n=[mt_n 'salinity']; mt_d=[mt_d 'iPROF,iDEPTH'];
mt_v=[mt_v 'prof_Sweight']; mt_u=[mt_u '(psu)^-2']; mt_n=[mt_n 'salinity least-square weight']; mt_d=[mt_d 'iPROF,iDEPTH'];
mt_v=[mt_v 'prof_Sestim']; mt_u=[mt_u 'psu']; mt_n=[mt_n 'salinity estimate (e.g. from atlas)']; mt_d=[mt_d 'iPROF,iDEPTH'];
mt_v=[mt_v 'prof_Serr']; mt_u=[mt_u 'psu']; mt_n=[mt_n 'salinity instrumental error']; mt_d=[mt_d 'iPROF,iDEPTH'];
mt_v=[mt_v 'prof_Sflag']; mt_u=[mt_u ' ']; mt_n=[mt_n 'flag = i > 0 means test i rejected data.']; mt_d=[mt_d 'iPROF,iDEPTH'];
end;
%
if(flagU==1);
mt_v=[mt_v 'prof_U']; mt_u=[mt_u 'm/s']; mt_n=[mt_n 'eastward velocity comp.']; mt_d=[mt_d 'iPROF,iDEPTH'];
mt_v=[mt_v 'prof_Uweight']; mt_u=[mt_u '(m/s)^-2']; mt_n=[mt_n 'east. v. least-square weight']; mt_d=[mt_d 'iPROF,iDEPTH'];
mt_v=[mt_v 'prof_V']; mt_u=[mt_u 'm/s']; mt_n=[mt_n 'northward velocity comp.']; mt_d=[mt_d 'iPROF,iDEPTH'];
mt_v=[mt_v 'prof_Vweight']; mt_u=[mt_u '(m/s)^-2']; mt_n=[mt_n 'north. v. least-square weight']; mt_d=[mt_d 'iPROF,iDEPTH'];
end;
if(flagPT==1);
mt_v=[mt_v 'prof_ptr']; mt_u=[mt_u 'X']; mt_n=[mt_n 'passive tracer']; mt_d=[mt_d 'iPROF,iDEPTH'];
mt_v=[mt_v 'prof_ptrweight']; mt_u=[mt_u '(X)^-2']; mt_n=[mt_n 'pass. tracer least-square weight']; mt_d=[mt_d 'iPROF,iDEPTH'];
end;
%
if(flagBP==1);
mt_v=[mt_v 'prof_bp']; mt_u=[mt_u 'cm']; mt_n=[mt_n 'bottom pressure']; mt_d=[mt_d 'iPROF'];
mt_v=[mt_v 'prof_bpweight']; mt_u=[mt_u '(cm)^-2']; mt_n=[mt_n 'bot. pres. least-square weight']; mt_d=[mt_d 'iPROF'];
end;
if(flagSSH==1);
mt_v=[mt_v 'prof_ssh']; mt_u=[mt_u 'cm']; mt_n=[mt_n 'sea surface height']; mt_d=[mt_d 'iPROF'];
mt_v=[mt_v 'prof_sshweight']; mt_u=[mt_u '(cm)^-2']; mt_n=[mt_n 'ssh least-square weight']; mt_d=[mt_d 'iPROF'];
end;



%=============================create the file=================================%
% write the netcdf structure

global useNativeMatlabNetcdf;
if isempty(useNativeMatlabNetcdf); useNativeMatlabNetcdf = ~isempty(which('netcdf.open')); end;
%if(useNativeMatlabNetcdf);ncmodestr='CLOBBER';else;ncmodestr='CLOBBER';end;
ncmodestr='CLOBBER';
ncid=nccreate(fileOut,ncmodestr);

aa=sprintf(['Format: MITprof netcdf. This file was created using \n' ...
    'the matlab toolbox which can be obtained (see README) from \n'...
    'http://mitgcm.org/viewvc/MITgcm/MITgcm_contrib/gael/profilesMatlabProcessing/']);
ncputAtt(ncid,'','description',aa);
ncputAtt(ncid,'','date',date);

ncdefDim(ncid,'iPROF',iPROF);
ncdefDim(ncid,'iDEPTH',iDEPTH);
ncdefDim(ncid,'lTXT',lTXT);
ncdefDim(ncid,'iINTERP',iINTERP);

['dimensions of the file']
[iPROF iDEPTH lTXT iINTERP]

for ii=1:length(list_vars_plus);
    jj=find(strcmp(mt_v, list_vars_plus{ii}));
    if ~isempty(jj);
        
        if strcmp(mt_d{jj},'iPROF,iDEPTH');
            ncdefVar(ncid,mt_v{jj},'double',{'iDEPTH','iPROF'});%note the direction flip
            %ncdefVar(ncid,mt_v{jj},'double',{'iPROF','iDEPTH'});
        else;
            ncdefVar(ncid,mt_v{jj},'double',{mt_d{jj}});
        end;
    
        ncputAtt(ncid,mt_v{jj},'long_name',mt_n{jj});
        ncputAtt(ncid,mt_v{jj},'units',mt_u{jj});
        ncputAtt(ncid,mt_v{jj},'missing_value',fillval);
        ncputAtt(ncid,mt_v{jj},'_FillValue',fillval);
    
    else;
        if strcmp(list_vars_plus{ii},'prof_descr')
            ncdefVar(ncid,'prof_descr','char',{'lTXT','iPROF'});%think this is a mistake
            %ncdefVar(ncid,'prof_descr','char',{'iPROF','lTXT'});
            ncputAtt(ncid,'prof_descr','long_name','profile description');
        else
            warning([list_vars_plus{ii} ' not included -- it is not a MITprof variable']);
        end
    end;
end;

['defining variables'];
[iPROF];
for ii=1:length(list1d); 
    list1d{ii};
    ncdefVar(ncid,list1d{ii},'double',{'iPROF'});end;

[iINTERP iPROF];
for ii=1:length(list2d); 
    list2d{ii};
    ncdefVar(ncid,list2d{ii},'double',{'iINTERP','iPROF'});end;
  
['end defining variables'];
ncclose(ncid);

%===========================write data ======================================%
% write to file:

['writing data'];

nc=ncopen(fileOut,'write');
vars=ncvars(nc);
list_vars=intersect(vars,list_vars_plus);

for ii=1:length(list_vars);
    fprintf('%i ',ii);
    varname=list_vars{ii};
    data=getfield(MITprof,varname);
    
    if(isnumeric(data)>0);spval=ncgetFillVal(nc,varname);end;
    if isnumeric(data)
        if isempty(spval),
            warning(['no FillVal for ' varname ': use of -9999 default value']);
            spval=double(-9999);
        end
        data(isnan(data))=spval;
        MITprof=setfield(MITprof,varname,data);
    end
    %whos data;
    %varname;
    ncputvar(nc,varname,data);
end

for ii=1:length(listAll);
    %['===== ']
    listAll{ii}
    if isfield(MITprof, listAll{ii})
        tmp =  getfield(MITprof,listAll{ii});
        %size(tmp)
        ncputvar(nc,listAll{ii},tmp);
    end
end;

ncclose(nc);

['finished writing']
%%%%%%%%% done %%%%%%%%%

clear mt_* list_vars* nLev fillval ii jj lTXT nProf aa data spval varname vars
