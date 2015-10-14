# this file specifies the locations and name of all input directories
# this line is comment (there MUST be a space after #)

> INPUT_DATA_TIME_BOUNDS
ip_start_yr		2000
ip_end_yr		2013

> FORCING_DATA_DIRS

forcing_data_dir	/media/jaideep/WorkData/Fire_G
# var |	  dir
ts		ncep_reanalysis/ts
rh		ncep_reanalysis/rhum
wsp		ncep_reanalysis/wsp
pr		ncep_reanalysis/precip
npp		GPP_modis

> FORCING_VARIABLE_DATA
# name | unit 	|	prefix   |	start_yr | nyrs/file | val-type |	time_interp_mode
ts		K			air.sig995	2000		1			ins			hold
rh		%			rhum.sig995	2000		1			avg			hold
wsp		m/s			wsp.sig995	2000		1			ins			hold
pr		mm/day		pr_wtr.eatm	2000		1			sum			hold	# 0 means missing values were
npp		gC/m2/s		npp			2000		14			avg			hold	# set to zero to avoid loss of good data in lterp

# file name will be taken as "prefix.yyyy.nc" or "prefix.yyyy-yyyy.nc"
# value types: ins (instantaneous), sum, avg (not used as of now)
# time_interp_modes: auto, hold, lter (not used as of now)
#	hold = hold previous value till next value is available
#	lter = interpolate in-between values (using previous and next times)
#	auto = hold for avg variables, lter for ins variables, sum-conservative random for sum variables

> STATIC_INPUT_FILES
# var	|	file (all these files should be in the input folder)
msk			util_data/masks/surta_global_0.5_sl.nc # surta_india_0.2.nc
vegtype		forest_type/MODIS/ftmap_MOD12q1_0.5_9pft_sl.nc # ftmap_iirs_8pft.nc
albedo		albedo/albedo_monthly_cycle_1998.2011_sl.nc # albedo.avg.2004.nc
elev		util_data/elevation/elev.0.5-deg.nc # elev.0.5-deg.nc

lmois		lmois_spin_end.2000-12-31.nc 

> END



