% Generate daily v4 SSH maps in cm
% from January 1, 1993 to  ?,2015

% RADS data is Bigendian, sequential, binary, real*4
% lat, lon, ssh_meters, time_days_since_1985
% missing value is -9999

% generate file description table
%clear all, close all, 


for numOfDataSet= [5]     %[1:5] 

if numOfDataSet==1;

name='topex_jason'
tj=dir('RADS_v4_2016/topex_jason/*');
for i=3:length(tj) 
  fn=['RADS_v4_2016/topex_jason/' tj(i).name];
  lt=tj(i).bytes/16;
  tmp=readbin(fn,[4 lt]);
  tj(i).datenum=datenum(1985,1,1)+mean(tmp(4,:));
  end
  tj(find(cell2mat({tj.bytes})==0))=[];
  save FileDescription_topex_jason   tj

elseif numOfDataSet==2;

name='gfo'
g1=dir('RADS_v4_2016/gfo/*');
for i=3:length(g1), 
  fn=['RADS_v4_2016/gfo/' g1(i).name];
  lt=g1(i).bytes/16;
  tmp=readbin(fn,[4 lt]);
  g1(i).datenum=datenum(1985,1,1)+mean(tmp(4,:));
end
g1(find(cell2mat({g1.bytes})==0))=[];
save FileDescription_gfo g1

elseif numOfDataSet==3;

name='ers_env'
en=dir('RADS_v4_2016/ers_env/*');
for i=3:length(en), 
  fn=['RADS_v4_2016/ers_env/' en(i).name];
  lt=en(i).bytes/16;
  tmp=readbin(fn,[4 lt]);
  en(i).datenum=datenum(1985,1,1)+mean(tmp(4,:));
end
en(find(cell2mat({en.bytes})==0))=[];
save FileDescription_ers_env en

elseif numOfDataSet==4;

name='c2'
c2=dir('RADS_v4_2016/cryosat/*');
for i=3:length(c2), 
  fn=['RADS_v4_2016/cryosat/' c2(i).name];
  lt=c2(i).bytes/16;
  tmp=readbin(fn,[4 lt]);
  c2(i).datenum=datenum(1985,1,1)+mean(tmp(4,:));
end
c2(find(cell2mat({c2.bytes})==0))=[];
save FileDescription_cryosat c2

elseif numOfDataSet==5;

name='sa'
sa=dir('RADS_v4_2016/saral/*');
for i=3:length(sa), 
  fn=['RADS_v4_2016/saral/' sa(i).name];
  lt=sa(i).bytes/16;
  tmp=readbin(fn,[4 lt]);
  sa(i).datenum=datenum(1985,1,1)+mean(tmp(4,:));
end
sa(find(cell2mat({sa.bytes})==0))=[];
save FileDescription_sa sa

elseif numOfDataSet==6;

name='j1c'        % not used 
j1c=dir('RADS_v4_2016/j1c/*bin');
for i=3:length(j1c), 
  fn=['RADS_v4_2016/j1c/' j1c(i).name];
  lt=j1c(i).bytes/16;
  tmp=readbin(fn,[4 lt]);
  j1c(i).datenum=datenum(1985,1,1)+mean(tmp(4,:));
end
j1c(find(cell2mat({j1c.bytes})==0))=[];
save FileDescription_j1c j1c

end;
end;
