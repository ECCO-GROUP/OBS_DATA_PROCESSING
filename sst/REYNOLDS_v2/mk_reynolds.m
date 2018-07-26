
%clear all;close('all')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yrs=1992;yre=2017;
for nc=2:2
	switch nc
     case 1
%     startup_llc270_gcmfaces
     dirIn='/bay/data1/king/REYNOLDS_v2/'
     dirOut='REYNOLDS_monthly_llc270/'
     fOut='Reynolds_monthly_llc270_'

     case 2
%     startup_v4_gcmfaces
%    cd ../
%    startup_v4r2
%    cd king_gridded

     dirIn='/pike4/owang/DATA/sst/REYNOLDS_v2/'
     dirOut='20180326/REYNOLDS_monthly_v4/'
     fOut='reynolds_oiv2_r1_mar2018_'
%    cd dirIn

     end

%%%%%%%%%%%%%%%%%%%%%%%

% generate monthly Reynolds on llc270 grid

[LON,LAT]=meshgrid(.5:359.5,-89.5:89.5);LON=LON';LAT=LAT';
fmask=fopen([dirIn 'lstags.onedeg.dat'],'r','b');
sstmsk = fread(fmask, [360 180],'float32');

for nyr=yrs:yre      %1992:2012
  fin=fopen([dirIn 'Reynolds_orig_monthly_' num2str(nyr)],'r','b');
  fout=fopen([dirOut fOut num2str(nyr)],'w','b');
  %
  for nmo=1:12 
    [num2str(nyr) ' ' num2str(nmo)]
    sst=fread(fin,[360 180],'float32');
     xx=find(sstmsk==0);sst(xx)=NaN;
%   figure(1);imagesc(sst',[-2 30]);axis xy;colorbar;
%   title([' SST orig' num2str(nmo) ' ' num2str(nyr)])
    SST=gcmfaces_remap_2d(LON,LAT,sst,3);
%    XC = convert2array(mygrid.XC);YC=convert2array(mygrid.YC);
%    idx = find(abs(XC-61.5)<1e-3&abs(YC+58.836)<1e-3)
%    SSTarr = convert2array(SST);
%'find ind SST'
%    SSTarr(idx)
%    XC(idx)
%    YC(idx)
%
%Do not flag here, since finite does not work for type "gcmfaces" 
% flag with -9999
%   SST(~isfinite(SST))=-9999.;

    obsmap=convert2array(SST);
    obsmapcompact=convert2gcmfaces(SST);
%flag here with -9999.
    obsmap(~isfinite(obsmap))=-9999.; 
    obsmapcompact(~isfinite(obsmapcompact))=-9999.; 

    fwrite(fout,obsmapcompact,'float32');
    cc=[-2 30];
    figure(2);
     imagesc(obsmap',cc);axis xy; colorbar
    title([' SST ' num2str(nmo) ' ' num2str(nyr)])
    pause(1)


  end  % for nmo=
  fclose(fin);fclose(fout);
end   % for nyr=
end   %for nc=
