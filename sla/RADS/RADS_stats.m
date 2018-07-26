clear all;close('all');
startup_v4_gcmfaces
sat=['tj' 'g1' 'en' 'c2' 'sa'];
sat='sa'
%.......................................................
name1=[sat '_daily_ssh_v4_']
name2=[sat '_daily_ssh_v4_MAD_']       
ny1=1993;ny2=2015;
tit1=['VAR  RADS ' sat ' 1992-2015 (RADS ice flags)']   
tit2=['VAR  RADS ' sat ' 1992-2015 (RADS ice flags, MAD)']   
dirIn='RADS_v4_2016/';     

QQ=[];
for nyr=ny1:ny2    
	ndays=365;if(mod(nyr,4)==0);ndays=366;end
	fid=fopen([dirIn name1 num2str(nyr)],'r','b');Q=fread(fid,[105300 ndays],'float32');fclose(fid);
	xx=find(Q==-9999);Q(xx)=NaN*ones(size(xx));
	QQ=[QQ Q];
end
QQv=nanvar(QQ');
QQv=reshape(QQv,90,1170);
q=convert2gcmfaces(QQv);
qS=convert2array(q);
qS1=qS;

orient tall
subplot(311)
imagesc(qS',[0 300]);axis xy;colorbar
title([tit1]);

QQ=[];
for nyr=ny1:ny2    
	ndays=365;if(mod(nyr,4)==0);ndays=366;end
	fid=fopen([dirIn name2 num2str(nyr)],'r','b');Q=fread(fid,[105300 ndays],'float32');fclose(fid);
	xx=find(Q==-9999);Q(xx)=NaN*ones(size(xx));
	QQ=[QQ Q];
end
QQv=nanvar(QQ');
QQv=reshape(QQv,90,1170);
q=convert2gcmfaces(QQv);
qS=convert2array(q);

subplot(312)
imagesc(qS',[0 300]);axis xy;colorbar
title([tit2]);

subplot(313)
qS=qS-qS1;
imagesc(qS',[-20 20]);axis xy;colorbar
title('difference');

eval(['print -djpeg ',name2, '92-15_variance.jpg;'])
eval(['print -dpsc ',name2, '92-15_variance.ps;'])
fid=fopen([dirIn name2 '92-15_variance'],'w','b');
fwrite(fid,QQv,'float32');fclose(fid);

