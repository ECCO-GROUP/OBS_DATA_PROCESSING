%.........................
		name='sa'
qm=[];ny1=1992;ny2=2015;
for nyr=ny1:ny2    
	fid=fopen(['RADS_TJ_' num2str(nyr)],'r','b');Q=fread(fid,[105300 365],'float32');
	xx=find(isnan(Q));[nyr  size(xx)]
	xx=find(Q==-9999);Q(xx)=NaN*ones(size(xx));
	QM=nanmean(Q);qm=[qm QM];
end
x=ny1+1/365:1/365:ny2+1;
plot(x,qm)
axis([1992 2016 -8 8])
title('sa\_daily\_ssh\_v4');
%.........................
		















