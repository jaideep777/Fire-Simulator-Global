# this file specifies the locations and name of all input directories
# this line is comment (there MUST be a space after #)

> INPUT_DATA_TIME_BOUNDS
ip_start_yr		2000
ip_end_yr		2008

> FORCING_DATA_DIRS

forcing_data_dir	/media/jaideep/WorkData/Fire/ncep20cen/fire_input
# var |	  dir
ts		temp_sfc
rh		rhum_nsfc
wsp		wsp_10m
pr		precip_daily
npp		npp_monthly

> FORCING_VARIABLE_DATA
# name | unit 	|	prefix   |	start_yr | nyrs/file | val-type |	time_interp_mode
ts		K			air.sfc		2000		1			ins			hold
rh		%			rhum.sig995	2000		1			avg			hold
wsp		m/s			wsp.10m		2000		1			ins			hold
pr		mm/day		pr.iitm0	2000		9			sum			hold	# 0 means missing values were
npp		gC/m2/s		npp			1982		25			avg			hold	# set to zero to avoid loss of good data in lterp

# file name will be taken as "prefix.yyyy.nc" or "prefix.yyyy-yyyy.nc"
# value types: ins (instantaneous), sum, avg (not used as of now)
# time_interp_modes: auto, hold, lter (not used as of now)
#	hold = hold previous value till next value is available
#	lter = interpolate in-between values (using previous and next times)
#	auto = hold for avg variables, lter for ins variables, sum-conservative random for sum variables

> STATIC_INPUT_FILES
# var	|	file (all these files should be in the input folder)
msk			surta_india_0.2.nc
lmois		lmois_spin_end.2000-12-31.nc
vegtype		ftmap_iirs_8pft.nc
albedo		albedo.avg.2004.nc
elev		elev.0.5-deg.nc

> END


