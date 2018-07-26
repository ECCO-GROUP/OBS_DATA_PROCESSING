% call_rads_noice_mad_recompose invokes rads_noice_mad_recompose e.g.: 
%  rads_noice_mad_recompose('topexfile',[1992:2011],0); %production
%  rads_noice_mad_recompose('topexfile',[1992:2011],1); %test version

startup_v4_gcmfaces;

%rundir='/net/nares/raid8/diana/ecco_v4/RADS_noice_mad/' %EDIT THIS
%eval(['cd ' rundir]) 

rads_noice_mad_recompose('topexfile',[1992:2015],0);       

rads_noice_mad_recompose('ersfile',[1992:2015],0);       

rads_noice_mad_recompose('gfofile',[1992:2015],0);       

rads_noice_mad_recompose('safile',[1992:2015],0);       

rads_noice_mad_recompose('c2file',[1992:2015],0);       
