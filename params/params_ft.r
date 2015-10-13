# forest type parameters

> nPFT 9

#	X		AGR		BLE		NLE		BLD		NLD		GR		SCD		SCX		<-+ PFT
#	0		1		2		3		4		5		6		7		8		v Attr
fP	0		1		1		1		1		1		1		1		1		# fP = packing fraction
aL	0.00	0.45	0.25	0.25	0.30	0.30	0.45	0.30	0.30	# allocation to leaves
aS	0.00	0.00	0.50	0.50	0.45	0.45	0.00	0.45	0.45	# allocation to stem
LL	1		1		2		2		1		1		1		1		1		# leaf-lifespan
T	3		3		3		3		3		3		5		3		3		# halflife of dry litter decomposition (months)
LM	0		4		5		7		5		5		0		4		1		# LAI max	(REF: pft_2002, Bonan) set to 0 for grass because grass biomass is itself fuel
Lm	0		0		.1		.1		.1		0		0		0		0		# LAI min	(REF: pft_2002, Bonan)
ZM	-1		-1		-1		-1		3		3		10		2		3		# Month until which all leaves are shed / 1st leafless month (-1 if no leafless month)
Wcs	0		0		0.5		0.5		0.5		0.5		0.5		0.5		0.5		# canopy water holding capacity per leaf layer (kg/m2) (REF: Ogee_2002)
# sources for litter decomposition T1/2:
# sundarpandian 1999, pandey 2007, nelson 1999

rhobL 		10		# litter bulk density 
theta_sL 	0.8		# saturation water content of litter (m3/m3) 


# modified pheno
# REFS: Bhat 1992,	Bhadra 2011, folder = seasonalphenology	
# F = flush, M = mature, S = shedding, Z = leafless, E = both S & F.
#	X	AGR	BLE	NLE	BLD	NLD	GR	SCD SCX	<-+ PFT, X = barren
> PHENO
Jan	X	X	E	E	S	S	Z	S	S
Feb	X	X	E	E	S	S	Z	Z	S
Mar	X	X	E	E	Z	Z	Z	Z	Z
Apr	X	X	E	E	Z	Z	Z	Z	Z
May	X	X	E	E	F	F	F	F	Z
Jun	X	X	E	F	F	F	F	F	F
Jul	X	X	E	F	M	M	F	F	F
Aug	X	X	E	F	M	M	M	F	F
Sep	X	X	E	S	M	M	S	M	F
Oct	X	X	E	S	M	M	Z	M	M
Nov	X	X	E	S	S	S	Z	S	S
Dec	X	X	E	E	S	S	Z	S	S

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
#	------------------------------------
#	F = leaf flushing 		= 0
#	M = mature leaf			= 1
#	L = leaf-fall 			= 2
#	Z = leafless/dormant	= 3
#	D = drying out 			= 4
#	X = invalid 			= -1
# 	(REF: pheno_2005, Singh, Kushwaha)
#	------------------------------------
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

# Characteristic Carbon fixation rates
# tropics (-30, 30) deg lat, averaged over years: 2001, 2004, 2008, 2011, 2013 
#	 X		AGR		BLE		NLE		BLD		NLD		 GR		SCD		SCX	 ...PFTs
> CARBON_FIXATION_RATE_NORTH
J	1.47	20.56	85.31	 0.00	28.62	0.00	 5.24	12.21	3.27
F	1.71	21.37	77.38	 9.89	16.40	0.00	 2.20	 4.73	2.39
M	2.09	18.65	88.15	30.10	22.97	0.00	 2.22	 6.66	2.73
A	1.90	16.41	86.78	42.61	34.60	0.00	 2.79	 9.35	3.73
M	1.72	17.56	85.84	65.65	51.83	0.00	 2.73	15.78	4.18
J	1.68	22.89	78.04	79.58	55.77	0.00	 3.77	27.48	1.31
J	1.87	32.38	82.26	95.15	53.68	0.00	10.58	42.49	1.07
A	1.80	39.83	85.73	95.37	48.74	0.00	18.81	47.93	3.24
S	1.63	39.93	87.13	70.81	57.25	0.00	15.89	58.81	3.56
O	1.54	36.95	90.90	43.18	73.03	0.00	 5.30	53.28	5.56
N	1.25	25.27	85.79	13.00	67.21	0.00	 3.44	27.62	6.78
D	1.25	21.62	82.26	 0.00	46.23	0.00	 4.76	16.80	6.38

> CARBON_FIXATION_RATE_SOUTH
J	19.96	71.13	 85.13	 79.71	73.70	0.00	23.58	48.59	2.22
F	17.92	66.78	 75.17	 47.67	70.76	0.00	25.82	51.57	3.96
M	21.81	65.12	 88.75	102.85	84.19	0.00	24.07	60.28	6.13
A	19.46	61.36	 86.37	 95.62	72.28	0.00	20.32	56.76	4.88
M	17.03	53.58	 87.67	106.57	51.50	0.00	12.91	50.38	4.86
J	14.40	42.94	 85.82	 50.93	34.37	0.00	 6.39	39.21	3.63
J	15.54	41.47	 91.56	100.66	23.88	0.00	 6.64	32.69	4.14
A	17.79	41.53	 91.62	148.91	15.06	0.00	 9.39	22.96	5.36
S	19.23	41.21	 93.90	138.78	15.21	0.00	16.09	10.93	2.96
O	19.46	54.24	100.12	131.79	35.83	0.00	23.48	10.94	0.85
N	17.69	61.71	 91.56	104.65	57.19	0.00	26.03	18.74	0.22
D	18.96	62.58	 88.64	 75.63	67.72	0.00	22.16	35.43	0.87


# Negative numbers have been set to zero in above tables. For original data, see excel sheet.

> END

