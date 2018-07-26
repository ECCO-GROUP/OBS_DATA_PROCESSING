orient tall
name='tj'
qm=[];ny1=1993;ny2=2013;
for nyr=ny1:ny2    
	fid=fopen(['tj_daily_ssh_v4_r1_noice_' num2str(nyr)],'r','b');Q=fread(fid,[105300 365],'float32');
	xx=find(Q==-9999);Q(xx)=NaN*ones(size(xx));
	QM=nanmean(Q);qm=[qm QM];
end
x=ny1+1/365:1/365:ny2+1;
subplot(511);plot(x,qm)
axis([1992 2015 -6 6])
title('tj\_daily\_ssh\_v4  RADS ice flag, MAD');
%.........................

name='g1'
qm=[];ny1=2001;ny2=2008;
for nyr=ny1:ny2    
	fid=fopen(['g1_daily_ssh_v4_r1_noice_' num2str(nyr)],'r','b');Q=fread(fid,[105300 365],'float32');
	xx=find(Q==-9999);Q(xx)=NaN*ones(size(xx));
	QM=nanmean(Q);qm=[qm QM];
end
x=ny1+1/365:1/365:ny2+1;
subplot(512);plot(x,qm)
axis([1992 2015 -6 6])
title('g1\_daily\_ssh\_v4  RADS ice flag, MAD');
%.........................
		
		name='en'
qm=[];ny1=1992;ny2=2013;
for nyr=ny1:ny2    
	fid=fopen(['en_daily_ssh_v4_r1_noice_' num2str(nyr)],'r','b');Q=fread(fid,[105300 365],'float32');
	xx=find(Q==-9999);Q(xx)=NaN*ones(size(xx));
	QM=nanmean(Q);qm=[qm QM];
end
x=ny1+1/365:1/365:ny2+1;
subplot(513);plot(x,qm)
axis([1992 2015 -6 6])
title('en\_daily\_ssh\_v4  RADS ice flag, MAD');
%.........................
	% not redone with RADS ice flag, MAD combo	
		name='c2'
qm=[];ny1=2010;ny2=2013;
for nyr=ny1:ny2    
	fid=fopen(['c2_daily_ssh_v4_r1_noice_' num2str(nyr)],'r','b');Q=fread(fid,[105300 365],'float32');
	xx=find(Q==-9999);Q(xx)=NaN*ones(size(xx));
	QM=nanmean(Q);qm=[qm QM];
end
x=ny1+1/365:1/365:ny2+1;
subplot(514);plot(x,qm)
axis([1992 2014 -6 6])
title('c2\_daily\_ssh\_v4');
%.........................
	% not redone with RADS ice flag, MAD combo	
		name='j1c'
qm=[];ny1=2012;ny2=2012;
for nyr=ny1:ny2    
	fid=fopen(['c2_daily_ssh_v4_r1_noice_' num2str(nyr)],'r','b');Q=fread(fid,[105300 365],'float32');
	xx=find(Q==-9999);Q(xx)=NaN*ones(size(xx));
	QM=nanmean(Q);qm=[qm QM];
end
x=ny1+1/365:1/365:ny2+1;
subplot(515);plot(x,qm)
axis([1992 2014 -6 6])
title('j1c\_daily\_ssh\_v4');
%.........................
		
print -dpsc RADS_r1_RADSIceFlag.ps
print -djpeg RADS_r1_RADSIceFlag.jpg















