if 1;

%%%%%%%%%%%%%%%%%%%%%%%
%generate delaunay triangulation

%clear all, close all, 

startup_v4_gcmfaces
gcmfaces_global
gcmfaces_bindata

SSHout_9999=-9999*ones(90,1170);


%%%%%%%%%%%%%%%%%%%%%%%
% generate daily maps  of altimeter ssh (cm) on v4 grid

% RADS data is Bigendian, sequential, binary, real*4
% lat, lon, ssh_meters, time_days_since_1985
% missing value is -9999



load FileDescription_topex_jason.mat;
load FileDescription_gfo.mat;
load FileDescription_ers_env.mat;
load FileDescription_c2.mat;
load FileDescription_sa.mat;


end;

dirr='./';
doOutputArray=0;


%   for alt={'tj','g1','en','c2','sa'}
for alt={'sa'}

        if strcmp(alt{1},'tj')
	     yrs=1993;yre=2015;      % complete
        elseif strcmp(alt{1},'en')   %complete
	     yrs=1992;yre=2013;
        elseif strcmp(alt{1},'g1')    %complete
	     yrs=2001;yre=2008;
        elseif strcmp(alt{1},'c2')
	     yrs=2010;yre=2015;
        elseif strcmp(alt{1},'sa')
	     yrs=2013;yre=2015;
        else
        error(['string ' alt{1} ' not recognized'])
	end


	eval(['dte=cell2mat({' alt{1} '.datenum});'])
   for yr=yrs:yre          
	if doOutputArray; fout=[dirr alt{1} '_daily_ssh_v4array_' int2str(yr)]; end;
	foutm=[dirr alt{1} '_daily_ssh_v4_' int2str(yr)];
	if doOutputArray; fid=fopen(fout,'w','b'); end;
	fidm=fopen(foutm,'w','b');
	n=0; 
        fprintf([alt{1} ' -- ' num2str(yr) ' \n']);
   for t=datenum(yr,1,1):datenum(yr,12,31);
%         disp(datestr(t))
	 lon=[]; lat=[]; ssh=[];
	 it = find( dte >= t & dte <(t+1) );


  for i=1:length(it)
	if strcmp(alt{1},'tj')  
	fn=['RADS_v4_2016/topex_jason/' tj(it(i)).name];
	elseif strcmp(alt{1},'en')
	fn=['RADS_v4_2016/ers_env/' en(it(i)).name];
	elseif strcmp(alt{1},'g1')
	fn=['RADS_v4_2016/gfo/' g1(it(i)).name];
	elseif strcmp(alt{1},'sa')
	fn=['RADS_v4_2016/saral/' sa(it(i)).name];
	elseif strcmp(alt{1},'c2')
	fn=['RADS_v4_2016/cryosat/' c2(it(i)).name];
	else
	error(['string ' alt{1} ' not recognized'])
	end
	eval(['lh=' alt{1} '(it(i)).bytes/16;'])
	tmp=readbin(fn,[4 lh]);
	lat=[lat tmp(1,:)];
	lon=[lon tmp(2,:)];  xx=find(lon>180);lon(xx)=lon(xx)-360;
	ssh=[ssh tmp(3,:)];
   end     % i=1:length(it)

        if(length(it)~=0);   %CKCKC
        [SSHmap,NPT]=gcmfaces_bindata(lon,lat,ssh);    % faces format

	in=find(NPT); SSHmap(in)=100*SSHmap(in)./NPT(in);    % m > cm
	in=find(~NPT); SSHmap(in)=-9999;

        SSHarray=convert2array(SSHmap);   %(360,360)
        SSHout=convert2gcmfaces(SSHmap); 
	fwrite(fidm,SSHout,'float32');      % v4 compact
	else   %CKCKC
        fwrite(fidm,SSHout_9999,'float32');   %CKCKC
	end    %CKCKCKCK

	n=n+1;
	end   % t=datenumm...
	end   % yr=yrs:yre
        if doOutputArray; fclose(fid); end; 
        fclose(fidm);
	end   % alt=

	%%%%%%%%%%%%%%%%%%%%%%%
 
