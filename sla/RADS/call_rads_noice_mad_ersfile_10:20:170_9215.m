% To run this script in background:  
% matlab -nodesktop -nodisplay < call_rads_noice_mad_ersfile_10:20:170_9215.m >&! call_rads_noice_mad_ersfile_10:20:170_9215.log&

startup_v4_gcmfaces;
which rads_noice_mad
choiceData='ersfile';years=[1992:2015];doTesting=0;l0list=[10:20:170];L0list=[-80:20:80];doNoice=0;
for l0=l0list;
l0
for L0=L0list;
L0
  rads_noice_mad(choiceData,l0,L0,0,0)
end
end
disp('END OF JOB')
exit
