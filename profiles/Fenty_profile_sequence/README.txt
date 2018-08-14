

Example processing sequence for profile files that have already been formatting in the 'MITprof' format.
========================================================================================================

by Ian Fenty, 08/13/2018


This sequence was used to process profiles from NODC in May/June 2018
for use in llc90 or ll270 grids.

NODC profiles were taken to the 'MITprof' format.  Before they can be
used in the profiles package they need some additional fields.  These
include
0. profile and tile points
1. T/S climatology
2. spatial bin index for the 10242 faced geodesic
3. spatial bin index for the 2562  faced geodesic
4. sigmas/weights
5. a mask put on the sigmas/weights for bad / missing data
6. a conversion from in situ to potential T
7. another mask for where profile costs vs. the climatology
   are too high

After that, I have two routines that merge profiles
8. Merge profiles of a single type into several files each with n profile per file (or fewer)
9. Merge profiles of a single type into a single file

A new helper script writes the 'prof_gci' field.  This field isn't 
used by the profile package and as far I can tell is only used for
offline analysis.
'update_prof_gci_on_prepared_profiles.m'



A number of helper routines are needed for these to work.  At a minimum 
you need to have 

1. Gael's  profilesMatlabProcessing libraries in the path.
   * http://wwwcvs.mitgcm.org/viewvc/MITgcm/MITgcm_contrib/gael/profilesMatlabProcessing/
2. matlab seawater package
3. m_plot package


There are inevitably many other routines that are required that I neglected to include
in this directory.  When you find one, let me know and I'll add it.


Note:
These routines were written as functions where the main argument is the 'run_code'. 
the 'run_code' indicates which arguments the routine should use.  Example argument
include the path to the netcdf files, the name of the output files, the llc grid,
etc.   


========================================================
Example sequence that I used to process profiles from 
the NODC repository in May/June/July of 2018
========================================================

After the NODC profiles have been converted into a skeleton
MITprof structure, one may begin.


%% step 0: Add profile points and tile points on profiles
%%========================================================

clear
run_codes = {'20180508_llc90_NODC_all'}
update_prof_and_tile_points_on_profiles(run_codes)


%% step 1: add T/S monthly mean climatology to the profs
%%========================================================

clear
run_code='APPLY_WOA13V2_MONTHLY_CLIM_TO_ALL_WOD13_20180508'
write_profile_to_netcdf = 1;
make_figs = 0;

update_monthly_mean_TS_clim_WOA13v2_on_prepared_profiles(...
    run_code, write_profile_to_netcdf, make_figs)


%% step 2: find which face of the 10242 geodesic surface
%%           each profile is on
%%========================================================
clear
run_code = 'WOD_13_20180509_all_add_10242_geodesic'
update_spatial_bin_index_on_prepared_profiles(run_code)


%% step 3: find which face of the 2562 geodesic surface
%%           each profile is on
%%========================================================

clear
run_code = 'WOD_13_20180509_all_add_2562_geodesic'
update_spatial_bin_index_on_prepared_profiles(run_code)


%% step 4: add sigma/weight fields
%%========================================================

clear
run_code = 13
update_sigmaTS_on_prepared_profiles(run_code)


%% step 5: zero some weights where data is bad/missing
%%========================================================

clear
run_code = 'NODC_20180508'
update_zero_weight_points_on_prepared_profiles(run_code)


%% step 6: convert in situ T to potential T
%%========================================================

clear
run_code = 'NODC_20180508'
update_prof_insitu_T_to_potential_T(run_code)


%% step 7: zero some weights 
%%========================================================

run_code = 'NODC_20180508_high_cost_vs_clim_do_not_exclude_high_lats'
update_zero_weight_points_on_prepared_profiles(run_code)


%% step 8: zero some more weights
%%========================================================

clear
run_code = 'NODC_20180508_high_cost_vs_clim_do_not_exclude_high_lats'
update_remove_zero_T_S_weighted_profiles_from_MITprof(run_code)


%% step 9: combine profiles of different types to 
%%         MITprof files of fixed length (100000 profs 
%%         per tile)
%%========================================================

clear
prof_types = {'APB', 'CTD', 'GLD', 'MRB', 'XBT'}

max_length = 100000
in_dir = '/home/ifenty/data/observations/insitu/NODC/NODC_20180508/all/step_08b_remove_zero_T_S_weight_profs/'
out_dir = '/home/ifenty/data/observations/insitu/NODC/NODC_20180508/all/step_09b_aggregated_into_fixed_length/'
write_output=1;

for pi = 1:length(prof_types)
    p_type = prof_types{pi}
    in_name = p_type
    out_name = [in_name '_WOD13_20180508_step_09']

    aggregate_profiles_into_fixed_length_files(in_name, in_dir, ...
        out_name, out_dir, max_length)

end



%% step 10: aggregate profiles of different types (APB, CTD
%%          etc., to single files
%%========================================================

clear
prof_types = {'APB', 'CTD', 'GLD', 'MRB', 'XBT'}

in_dir = '/home/ifenty/data/observations/insitu/NODC/NODC_20180508/all/step_09b_aggregated_into_fixed_length/'
out_dir = '/home/ifenty/data/observations/insitu/NODC/NODC_20180508/all/step_10b_aggregated_into_single_files_by_type/'
write_output=1;

for pi = 1:length(prof_types)
    p_type = prof_types{pi}
    in_name = p_type
    out_name = 'AGGREGATED'

    aggregate_profiles_into_single_MITprof(in_name, in_dir, out_name, ...
      out_dir, write_output)
end

