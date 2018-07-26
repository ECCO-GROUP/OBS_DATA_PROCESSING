%   calculate REYNOLDS_v2 annual temp anom

y1=1992;y2=2012;nyrs=y2-y1+1;
nx=360;ny=180;x=.5:359.5;y=-89.5:89.5;

Qm=zeros(nx*ny,nyrs);
nrec=0;
for n=y1:y2
	nrec=nrec+1
	fin=fopen(['Reynolds_orig_monthly_' num2str(n)],'r','b');
	Q=fread(fin,[nx*ny 12],'float32');fclose(fin);
	xx=find(Q==0);Q(xx)=NaN*ones(size(xx));
	Qm(:,nrec)=nanmean(Q');
	d=reshape(Qm(:,nrec),nx,ny);
	imagesc(d');axis xy;colorbar;title([num2str(n)]);pause(1)
end  % for n=y1:y2
QM=nanmean(Qm');
d=reshape(QM,nx,ny);
imagesc(d');axis xy;colorbar;
title(['Reynolds\_v2 ave temp ' num2str(y1) '-' num2str(y2)]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
years=1992:2012;AnnualSST=reshape(Qm,nx,ny,nyrs);meanSST=reshape(QM,nx,ny);lat=y;lon=x';
save ReynoldsAnnualSST years AnnualSST meanSST nx ny lat lon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

orient tall
Qa=zeros(nx*ny,nyrs);
nrec=0;np=0;clf;npage=0;
for n=y1:y2
        nrec=nrec+1;np=np+1;
	Qa(:,nrec)=Qm(:,nrec)-QM';
	d=reshape(Qa(:,nrec),nx,ny);
	subplot(4,1,np);
        imagesc(x,y,d',[-2 2]);axis xy;colorbar;

xx=find(isnan(d));d(xx)=zeros(size(xx));
ylat=y/57.2958;wsum=0;spmn=0;spsd=0;
for j=1:ny; w=cos(ylat(j));
for i=1:nx;
    if(d(i,j)~=0);
          wsum=wsum+w;spmn=spmn+d(i,j)*w;
	          spsd=spsd+d(i,j)*d(i,j)*w;
    end;
end
end
spmn=spmn/wsum;var=spsd/wsum-spmn*spmn;spsd=sqrt(var);
[spmn var spsd]

title(['Reynolds Temp anom ' num2str(n) '  awm=' num2str(spmn) '  ' num2str(var) '  ' num2str(spsd)]);pause(1)
if(mod(np,4)==0);
npage=npage+1
eval(['print -dpsc Reynolds_annual_anom_' num2str(npage) '.ps'])
eval(['print -djpeg Reynolds_annual_anom_' num2str(npage) '.jpg'])
np=0;pause(1);clf
end
end


