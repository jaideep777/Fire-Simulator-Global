using namespace std;
#include "../include/pheno.h"
#include "../include/globals.h"
#include "../include/vars.h"



/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	--> calc_alloc_pft(...)
	
	Allocate the Total NPP generation among different PFTs.
	
	the carbon fixation rates (R1, R2, ..) in the ft_params were calculated
	by running a regression over NPP in each cell for 10 yrs.
	
	Thus, NPP(cell) = f1*R1 + f2*R2 + ... (where R1, R2 are fixation rates 
										   and f1, f2 are PFT fractions)
	This is expected NPP. But observed NPP may be different due to global 
	trends (increasing CO2 etc). 
		Let NPP_obs = beta* NPP_exp

	Then N1 (NPP fixed in 1st PFT) is  N1 = (f1*R1)*beta, and so on,
	so that the above is satisfied.
	
	Input: month, observed NPP
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
//ofstream fout_allocdebug("../output/alloc_debug.txt");
vector <float> calc_alloc_pft(int m, float npp_obs, vector <float> &pft_fracs){
	vector <float> allocs(npft,0);
	m = m-1;	// convert month from 1-12 to 0-11 for indexing
	
	float npp_exp = 0;
	for (int ipft=0; ipft<npft; ++ipft){
		float f = pft_fracs[ipft];
		float Ni = f*rFixC[IX2(ipft,m, npft)];
		allocs[ipft] = Ni;
		npp_exp += Ni;	// expected NPP from regression relation
	}
	float beta = npp_obs/(npp_exp + 1e-12);  // 1e-12 to avoid NaN
	
//	fout_allocdebug << npp_obs*hrsPerMonth << "\t" << npp_exp << " \t";  
//	for (int i=0; i< npft; ++i){
//		fout_allocdebug << allocs[i] << "\t";
//	}
//	fout_allocdebug << "\n";	

	for (int ipft=0; ipft<npft; ++ipft) allocs[ipft] *= beta;

//	for (int i=0; i< npft; ++i){
//		fout_allocdebug << allocs[i] << "\t";
//	}
//	fout_allocdebug << "\n";	

	return allocs;
	
}


/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	--> set_leaf_alloc()
	
	Allocate NPP to leaves, stem and roots.
	
	npp allocation fractions depend on phenology stage as well as pft.
	this function sets the vectors aL and aS (for all PFTs) given the month.
	
	Allocation scheme:
	F (up-to last F month) -> all NPP to leaves
	F (last F month) -> normal growth 
	M, S, Z -> stem and roots only
	E -> normal growth
	
	input: month (1-12), NPP for each PFT
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
int set_leaf_alloc(int m){
	m = m-1;	// convert month range from (1-12) to (0-11) 

	// get alloc fractions for current month
	for (int i=0; i< npft; ++i){
		int ps_now = phenoStages[IX2(i,m, npft)];			// current pheno stage
		int ps_next = phenoStages[IX2(i,(m+1)%12, npft)];	// next month pheno stage
		if (ps_now == psF){	
			// this month is in flushing stage
			if (ps_next == psF)	{ 
				// next month also in flushing, so full alloc to leaves
				aL[i] = 1;	
				aS[i] = 0;
			}
			else {	
				// next month is not in flushing, so go to normal growth mode (allocs from default values)
				aL[i] = aLf[i];	
				aS[i] = aSf[i];
			}
		}
		else if (ps_now == psE){
			// this month is in flushing as well as shedding (evergreen pft)
			aL[i] = aLf[i];
			aS[i] = aSf[i];
		}
		else{ 
			// this month is not in flushing stage, so aL = 0
			aL[i] = 0;
			aS[i] = aSf[i] / (1-aLf[i]);
		}
		
//		// if plant is in flushing stage but npp_allocation is -ve, aL = 1, aS = -1 i.e. stem biomass converts to leaf biomass
//		if (nppAllocs[i] < 0){
//			aS[i] = 1;	// biomass is lost from stem
//			aL[i] = (ps_now == psF)? -1:0;  // and goes to leaves if flusing, otherwise to air
//			// (note reversed sign because they are going to get multiplied with -ve NPP)
//		}
	}
}


/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	--> calc_litterfall_rate()
	
	calculate of litter-fall, given the month and canopy biomass
	1. for evergreen trees, it is set acc to 1st order rate, dL/dt = -kC
	2. for deciduous trees (those with leafless period) it is linear so 
	   as to exhaust all leaves till leafless month   	
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
vector <float> calc_litterfall_rate(double gtime, vector <float> &canbio_now, float delT){

	vector <float> ls_rates(npft,0);	

	int m = gt2month(gtime)-1; // get current month and convert month (1-12) to index (0-11) 
	for (int i=0; i<npft; ++i){
		// shedding stage
		if (phenoStages[IX2(i,m,npft)] == psS){ 
			if (z1Month[i] >= 0){	// leafless month specified, i.e. deciduous tree
				// shed leaves so as to shed all till 1st leafless day
				int zmonth = z1Month[i];
				int yrnow = gt2year(gtime);
				if (zmonth <= m+1) ++yrnow;
				float leafless_start_day = ymd2gday(yrnow, zmonth, 1);
				
				int shed_tsteps = int((leafless_start_day - gtime)*24/delT) + 1;	// floor function
				float shed_hrs = shed_tsteps*delT; 

				ls_rates[i] = canbio_now[i]/(shed_hrs/hrsPerMonth);	// per month

				//if (ls_rates[i] > canbio_now[i]) ls_rates[i] = canbio_now[i]; // just to avoid -ve canbio
			}
			else{	// no leafless month, i.e. evergreen tree
				// shed leaves at 1st order rate
				ls_rates[i] = 1/(leafLs[i]*12)*canbio_now[i];	// per month. 
			}
		}
		// evergreen tree in E phase
		else if (phenoStages[IX2(i,m,npft)] == psE ){
			// shed leaves at 1st order rate
			ls_rates[i] = 1/(leafLs[i]*12)*canbio_now[i];	// per month. 
		}
		// all other stages for all PFTs
		else{
			ls_rates[i] = 0;	// no shedding if tree not in shedding phase
		}
	}

	return ls_rates;
}


char ps2char(int c);
///*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//	--> calc_pheno()

//	Calculate canbio and dxl
//	
//	input: current GT, timestep (in hrs) 
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
int calc_pheno(float gtime, float delT){

	int curr_month = gt2month(gtime); 	// month from 1-12 
	
	// calculate allocation fractions to leaves (aL) and stem (aS) 
	set_leaf_alloc(curr_month);
	

	for (int ilat=0; ilat<mgnlats; ++ilat){
		for (int ilon=0; ilon<mgnlons; ++ilon){

			float N_obs = npp(ilon, ilat, 0); 
			//N_obs *= hrsPerMonth;		// convert from gm/m2/hr to gm/m2/month 
			vector <float> pft_fracs(npft,0);
			for (int i=0; i<npft; ++i) pft_fracs[i] = vegtype(ilon, ilat, i); 

			// calculate fraction of total NPP going to each PFT (units same as N_obs) 
			vector <float> allocs_pft = calc_alloc_pft(curr_month, N_obs, pft_fracs);

			vector <float> canbiof(npft, 0);
			//vector <float> stembiof(npft, 0);
	
			for (int i=0; i<npft; ++i){
				canbiof[i] = 2*allocs_pft[i]*aL[i] *delT;		// gC/m2/hr* hrs * 2 gm/gC
				//stembiof[i] = allocs_pft[i]*aS[i] *delT;
		
				canbio(ilon, ilat, i) += canbiof[i];
				//stembio_cumm[i] += stembiof[i];
			}

			// calculate litter-fall rates (gm/month)
			vector <float> canbio_c(npft,0);
			for (int i=0; i<npft; ++i) canbio_c[i] = canbio(ilon, ilat, i); 
			vector <float> ls_rates = calc_litterfall_rate(gtime, canbio_c, delT);
	
			// accumulate litter, shed canopy
			for (int i=0; i<npft; ++i){
				float bio_shed = ls_rates[i]/hrsPerMonth*delT;	// gm/month * months/hr * hrs
				canbio(ilon, ilat, i) -= bio_shed;
				dxl(ilon, ilat, 0) += bio_shed;
			}

			// decay litter (dry decay rates for prerun)
			float K_dry = 0.693/(Tdecomp[0]*hrsPerMonth);
			float S = lmois(ilon,ilat,0);
			float beta = 1+0.25*(1-cos(S*pi))*(1-cos(S*pi));	// varies between 1-2 as lmois goes from 0-1

			dxl(ilon, ilat, 0) -=  K_dry*beta * dxl(ilon, ilat, 0)*delT;

		}	// ilon loop ends
	}	// ilat loop ends
	
	// update current time in gVars
	canbio.t = dxl.t = gtime;
	
	// output
	if (spout_on){
		int m = curr_month-1;
		sp_fout << gt2string(gtime) << "\t"; 		
		int u=3;
		sp_fout << ps2char(phenoStages[IX2(u,m, npft)]) << "\t"
				<< canbio.getCellValue(xlon,xlat,u) << "\t" << dxl.getCellValue(xlon,xlat,0) << "\t";
		u=5;
		sp_fout << ps2char(phenoStages[IX2(u,m, npft)]) << "\t"
				<< canbio.getCellValue(xlon,xlat,u) << "\t" << dxl.getCellValue(xlon,xlat,0) << "\t";
		//if (canbio_prerun_on) sp_fout << "\n";
		
	}

}






/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	--> calc_pheno()
	SINGLE POINT VERSION (BACKUP)
	Calculate canbio and dxl
	
	input: current GT, timestep (in hrs) 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
//int calc_pheno(float gtime, float delT){

//	int curr_month = gt2month(gtime); 	// month from 1-12 
//	
//	// calculate allocation fractions to leaves (aL) and stem (aS) 
//	set_leaf_alloc(curr_month);
//	
//	float N_obs = npp(i_xlon, i_xlat, 0); 
//	//N_obs *= hrsPerMonth;		// convert from gm/m2/hr to gm/m2/month 
//	vector <float> pft_fracs(npft,0);
//	for (int i=0; i<npft; ++i) pft_fracs[i] = vegtype(i_xlon, i_xlat, i); 

//	// calculate fraction of total NPP going to each PFT (units same as N_obs) 
//	vector <float> allocs_pft = calc_alloc_pft(curr_month, N_obs, pft_fracs);

//	vector <float> canbiof(npft, 0);
//	vector <float> stembiof(npft, 0);
//	
//	for (int i=0; i<npft; ++i){
//		canbiof[i] = allocs_pft[i]*aL[i] *delT;		// gm/m2/hr* hrs
//		stembiof[i] = allocs_pft[i]*aS[i] *delT;
//		
//		canbio_cumm[i] += canbiof[i];
//		stembio_cumm[i] += stembiof[i];
//	}

//	// calculate litter-fall rates (gm/month)
//	vector <float> ls_rates = calc_litterfall_rate(gtime, canbio_cumm);
//	
//	// accumulate litter, shed canopy, decay litter (dry decay rates for prerun)
//	for (int i=0; i<npft; ++i){
//		float bio_shed = ls_rates[i]/hrsPerMonth*delT;	// gm/month * months/hr * hrs
//		canbio_cumm[i] -= bio_shed;
//		if (canbio_cumm[i] < 0) canbio_cumm[i] = 0; 
//		littbio_cumm[i] += bio_shed;

//		littbio_cumm[i] -= 0.693/(Tdecomp[i]*hrsPerMonth) * littbio_cumm[i]*delT;
//	}

//	// output
//	if (spout_on){
//		int m = curr_month-1;
//		sp_fout << gt2string(gtime) << "\t"; 		
//		sp_fout << N_obs*hrsPerMonth << "\t";
//		//for (int i=0; i<npft; ++i) sp_fout << allocs_pft[i] << "\t";
//		int u=3;
//		sp_fout << allocs_pft[u]*hrsPerMonth << "\t";
//		sp_fout << canbiof[u] << "\t" << stembiof[u] << "\t" << ps2char(phenoStages[IX2(u,m, npft)]) << "\t" << rFixC[IX2(u,m, npft)] << "\t"
//				<< canbio_cumm[u] << "\t" << littbio_cumm[u] << "\t";
//		u=5;
//		sp_fout << allocs_pft[u]*hrsPerMonth << "\t";
//		sp_fout << canbiof[u] << "\t" << stembiof[u] << "\t" << ps2char(phenoStages[IX2(u,m, npft)]) << "\t" << rFixC[IX2(u,m, npft)] << "\t"
//				<< canbio_cumm[u] << "\t" << littbio_cumm[u] << "\t";
//		sp_fout << "\n";
//		
//	}

//}




