# this file specifies all simulation config parameters
# this line is comment (there MUST be a space after #)
# headings followed by ">" must be exactly as given here and in the same order

> TIME
# time is in GMT
# spinup MUST BE ENABLED for now
spinup		on			# (on/off) off implies IC will be supplied from file
spin_date0	2005-1-1	# time assumed to be 0:0:0. Ideal start date is Apr 1
sim_date0	2006-1-1
sim_t0		0:0:0		# GMT
sim_datef	2012-12-31
sim_tf		21:0:0
dt			6			# hours
base_date	2000-1-1	# base time assumed to be 0:0:0
# spinup end date is automatically taken as 1 time step lesser than sim_date0_t0

> MODEL_GRID
lon0	0.0
lonf	30.0
lat0	-30.0
latf	30.0
dlat	0.5
dlon	0.5
# output params for single point (for testing/characterizing) 
# Some points:
# 300.5 -19.5  - 90% BLD forest in Brazil
# 62.0 -27.5  - 50% GR in Brazil
# 76.5 11.5   - Mudumalai (50% BLE + 25% BLD)
# 21.5 -18.5  - 90% SCD in Africa
# 31.0 9.5    - 80% SCD below Sahara
# 
xlon	22.0	# 76 31 50
xlat	-2.0	# 11 35 41
pointOutFile	../output/point_run.txt
SPout_on	1	# 0/1

> OUTPUT_VARIABLES
# var		nc	sp  <- write flags to nc file and to single point output
pr			0	0
rh			0	0
ts			0	0
wsp			0	0
npp			0	0

evap		0	0
ndr			1	0
ps			0	0
dxl			1	0
canbio		1	0
canbio_max	0	0
lmois		1	0
cmois		1	0
fire		1	0

dfire		1	0
dts			1	0
dwsp		1	0
drh			1	0
dlmois		1	0
ddxl		1	0

# for debugging purposes only.
# using this makes the code slightly tedious but immensely easy to change quickly
> VARS_TO_USE
pr			1
rh			1
ts			1
wsp			1
npp			1

evap		0
ndr			1
ps			1
dxl			1
canbio		1
canbio_max	1
lmois		1
cmois		1
fire		1

dfire		1
dts			1
dwsp		1
drh			1
dlmois		1
ddxl		1

> END

