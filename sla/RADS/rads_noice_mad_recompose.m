function []=rads_noice_mad_recompose(choiceData,years,doTesting);

if isempty(who('doTesting')); doTesting=0; end;

%directories
dirData='RADS_v4_2016/'; %.mat files
dirDataOld='RADS_v4_2016/'; % v4 files
topexfile       = 'tj_daily_ssh_v4';
ersfile         = 'en_daily_ssh_v4';
gfofile         = 'g1_daily_ssh_v4';
safile          = 'sa_daily_ssh_v4';
c2file          = 'c2_daily_ssh_v4';
dirIce='input_nsidc_all/';
fileIce='nsidc79_daily';
dirOutput='./';
skipLoad=0;

nday=zeros(23,1);     %CKCKCK 21 for 92-13
for yy=years;
  tmp1=dir([dirDataOld topexfile '_' num2str(yy)]);
  nday(yy-1991)=tmp1.bytes/90/1170/4;
end;

%testing switches
if doTesting; years=[2009:2010]; choiceData='topexfile'; end;

gcmfaces_global;
eval(['fileData=' choiceData ';']);

reject.l0=[-170:20:170]'*ones(1,9);
reject.L0=ones(18,1)*[-80:20:80];
reject.remain=NaN*zeros(18,9,20);
reject.iceRejectPerc=NaN*zeros(18,9,20);
reject.madRejectPerc=NaN*zeros(18,9,20);

for yy=years;
  fprintf(['reading ' fileData '_' num2str(yy) '\n']);
  tmp1=dir([dirDataOld fileData '_' num2str(yy)]); nrec=tmp1.bytes/90/1170/4;
  if yy==1992; rec00=0; else; rec00=sum(nday(1:yy-1992)); end;
  rec0=rec00+1; rec1=rec00+nrec;
  goodData=NaN*zeros(105300,nrec);
for L0=[-80:20:80];
for l0=[-170:20:170];

%select region
XC=convert2gcmfaces(mygrid.XC.*mygrid.mskC(:,:,1)); XC=XC(:);
YC=convert2gcmfaces(mygrid.YC.*mygrid.mskC(:,:,1)); YC=YC(:);
ii=find(XC>=l0-10&XC<=l0+10&YC>=L0-10&YC<=L0+10);
ni=length(ii);

%read region
tile=load([dirData fileData 'l0is' num2str(l0) 'L0is' num2str(L0)]);

%add to accounting at proper location
[i0,j0]=find(reject.l0==l0&reject.L0==L0); y0=yy-1991;
reject.remain(i0,j0,y0)=sum(~isnan(tile.goodData(:)));
reject.iceRejectPerc(i0,j0,y0)=tile.iceRejectPerc;
reject.madRejectPerc(i0,j0,y0)=tile.madRejectPerc;

%add to goodData at proper location
if length(ii)~=length(tile.ii); error('inconsistent size'); end;
if max(abs(ii-tile.ii))~=0; error('inconsistent size'); end;
goodData(ii,:)=tile.goodData(:,rec0:rec1);

end;
end;

goodData=goodData*100;
goodData(isnan(goodData))=-9999;

fid=fopen([dirOutput fileData '_MAD_' num2str(yy)],'w','b');
fwrite(fid,goodData,'float32');
fclose(fid);
end;

