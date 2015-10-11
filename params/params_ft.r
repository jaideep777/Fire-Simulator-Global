# forest type parameters

> nPFT 8

#	X		AGR		NLE		BLE		MD		DD		GR		SC		<-+ PFT
#	0		1		2		3		4		5		6		7		v Attr
fP	0		1		1		1		1		1		1		1		# fP = packing fraction
aL	0.00	0.45	0.25	0.25	0.30	0.30	0.45	0.30	# allocation to leaves
aS	0.00	0.00	0.50	0.50	0.45	0.45	0.00	0.45	# allocation to stem
LL	1		1		2		2		1		1		1		1		# leaf-lifespan
T	3		3		3		3		3		3		5		3		# halflife of dry litter decomposition (months)
LM	0		4		5		7		5		5		0		4		# LAI max	(REF: pft_2002, Bonan) set to 0 for grass because grass biomass is itself fuel
Lm	0		0		.1		.1		.1		0		0		0		# LAI min	(REF: pft_2002, Bonan)
ZM	-1		-1		-1		-1		-1		2		2		2		# Month until which all leaves are shed / 1st leafless month (-1 if no leafless month)
Wcs	0		0		0.5		0.5		0.5		0.5		0.5		0.5		# canopy water holding capacity per leaf layer (kg/m2) (REF: Ogee_2002)
# sources for litter decomposition T1/2:
# sundarpandian 1999, pandey 2007, nelson 1999

rhobL 		10		# litter bulk density 
theta_sL 	0.8		# saturation water content of litter (m3/m3) 


# modified pheno
# REFS: Bhat 1992,	Bhadra 2011, folder = seasonalphenology	
# F = flush, M = mature, S = shedding, Z = leafless, E = both S & F. For GR, S = drying
#	X	AGR	NLE	BLE	MD	DD	GR	SC	<-+ PFT, X = barren
> PHENO
J	X	E	E	E	S	S	S	S
F	X	E	E	E	F	Z	Z	Z
M	X	E	E	E	F	Z	Z	Z
A	X	E	E	E	F	Z	Z	F
M	X	E	E	E	E	F	Z	F
J	X	E	E	M	E	F	F	F
J	X	E	E	M	E	F	F	F
A	X	E	E	M	E	F	F	F
S	X	E	E	E	E	F	M	M
O	X	E	E	E	E	M	M	S
N	X	E	E	E	S	S	M	S
D	X	E	E	E	S	S	M	S

#J	X	F	B	B	S	S	M	S
#F	X	F	B	B	S	Z	S	Z
#M	X	F	B	B	S	Z	Z	Z
#A	X	F	B	B	Z	Z	Z	Z
#M	X	F	B	B	F	F	Z	F
#J	X	F	B	B	F	F	F	F
#J	X	F	B	B	F	F	F	F
#A	X	F	B	B	F	F	F	F
#S	X	F	B	B	F	F	M	F
#O	X	F	B	B	M	M	M	M
#N	X	F	B	B	M	M	M	M
#D	X	F	B	B	M	S	M	S



# above matrix must be in this same order. see example below.
# There cant be comments at end of row
# rows are months from Jan to Dec
#	Phenology stage
#	------------------
#	F = leaf flushing 		= 0
#	M = mature leaf			= 1
#	L = leaf-fall 			= 2
#	Z = leafless/dormant	= 3
#	D = drying out 			= 4
#	X = invalid 			= -1
# 	(REF: pheno_2005, Singh, Kushwaha)
#	-------------------
#	X	AGR	NLE	BLE	MD	DD	GR	SCE	SCD	...PFTs
#	 original pheno from singh et al.
#	0	1	2	3	4	5	6	7	8	
#	X	D	L	L	L	L	M	L	J	9
#	X	D	L	L	L	Z	D	Z	F	11
#	X	D	F	F	L	Z	D	Z	M	14
#	X	D	F	F	Z	Z	D	Z	A	27
#	X	D	F	F	F	F	D	F	M	52
#	X	F	F	F	F	F	F	F	J	159
#	X	F	F	F	F	F	F	F	J	269
#	X	F	M	M	F	F	F	F	A	230
#	X	F	M	M	F	M	M	M	S	154
#	X	F	M	M	M	M	M	M	O	60
#	X	F	M	M	M	M	M	M	N	26
#	X	D	M	M	M	L	M	L	D	10

# Characteristic Carbon fixation rates, which is just a scale free representation for ... 
# ... monthly npp for 8 pfts (found by regression). Rows are months starting Jan
# generated from NPP data 1982-1991 averaged over the years for each month ...
# ... and then running a regression with 0 intercept over all spatial points
# +ve values -> Biomass is fixed from atmospheric CO2, 
# -ve values -> stem biomass is consumed in respiration, initial leafout 
# 				(depending on pheno stage)
#	 X		AGR		NLE		BLE		MD		DD		GR		SC	...PFTs
> CARBON_FIXATION_RATE
J	2.25	17.19	    0	96.59	52.67	25.12	3.70	10.33
F	2.56	16.85	    0	85.53	47.65	14.09	1.62	 4.68
M	3.29	14.74	20.52	98.82	52.04	    0	1.08	    0
A	2.14	9.16	58.13	115.56	47.54	    0	4.70	    0
M	0.16	7.61	87.37	147.88	54.00	    0	13.23	    0
J	   0	11.22	95.76	116.98	55.47	    0	12.91	 7.15
J	   0	20.82	98.67	99.95	69.14	30.42	14.63	21.33
A	   0	27.06	100.17	95.41	70.28	40.02	11.23	22.78
S	   0	29.61	75.67	114.42	79.54	50.73	9.45	17.75
O	0.91	29.21	39.13	128.21	98.33	49.94	8.27	10.41
N	1.52	21.21	6.00	120.56	82.59	32.73	11.32	 9.64
D	1.59	17.24	   0	106.16	62.40	31.31	8.72	14.86

# Original: Above is +ve numbers only
#J	2.25	17.19	-6.46	96.59	52.67	25.12	3.70	10.33
#F	2.56	16.85	-3.71	85.53	47.65	14.09	1.62	4.68
#M	3.29	14.74	20.52	98.82	52.04	-4.71	1.08	-8.90
#A	2.14	9.16	58.13	115.56	47.54	-12.08	4.70	-12.73
#M	0.16	7.61	87.37	147.88	54.00	-16.84	13.23	-8.90
#J	-1.08	11.22	95.76	116.98	55.47	-2.70	12.91	7.15
#J	-2.38	20.82	98.67	99.95	69.14	30.42	14.63	21.33
#A	-2.58	27.06	100.17	95.41	70.28	40.02	11.23	22.78
#S	-2.13	29.61	75.67	114.42	79.54	50.73	9.45	17.75
#O	0.91	29.21	39.13	128.21	98.33	49.94	8.27	10.41
#N	1.52	21.21	6.00	120.56	82.59	32.73	11.32	9.64
#D	1.59	17.24	-5.85	106.16	62.40	31.31	8.72	14.86


> END

