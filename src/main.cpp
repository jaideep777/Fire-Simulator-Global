#include <iostream>
#include <cmath>
#include <netcdfcpp.h>
#include <fstream>
#include <vector>
using namespace std;

#include <gsm.h>
#include "../include/globals.h"
#include "../include/vars.h"
#include "../include/init.h"
#include "../include/io.h"
#include "../include/prerun.h"


int main(){

	// set NETCDF error behavior to non-fatal
	NcError err(NcError::silent_nonfatal);

	// Turn on/off messages from GSM
	gsm_info_on = false;
	gsm_debug_on = false;
	gsm_warnings_on = true;
	gsm_errors_on = true;

	// init
	int i = init_infisim();
	if (i != 0) {cout << "** ERROR ** : init failed!\n\n"; return 1;};

	// pre-runs
	prerun_canbio_ic();
	sp_fout << "*************** LMOIS PR ********************\n";
	prerun_lmois_ic();

//	// main run
//	sp_fout << "************** MAIN *********************\n";
//	main_run();
//	
//	// destructor will close all the NC files :) ...phew!
//	cout << "> Check details in log file /output/log.txt\n";
//	cout << "> Successfully wrote NC files.\n\n\n";

}

/*


int main(){

	// set NETCDF error behavior to non-fatal
	NcError err(NcError::silent_nonfatal);

//	// init
//	init_infisim();

	// forest type map
	gVar ftmap, npp;
	NcFile_handle ftmap_handle, npp_handle;
	string ftmap_file = "../input/ftmap_iirs.nc", npp_file = "/media/WorkData/Fire/ncep20cen/fire_input/npp_monthly/npp.1982-2006.nc";

	ftmap_handle.open(ftmap_file, "r", glimits_india);
	ftmap_handle.readCoords(ftmap, log_fout);
	ftmap_handle.readVarAtts(ftmap);
	ftmap_handle.readVar(ftmap, 0);
	ftmap.printGrid();
	
	npp_handle.open(npp_file, "r", glimits_india);
	npp_handle.readCoords(npp, log_fout);
	npp_handle.readVarAtts(npp);
	npp_handle.readVar(npp, 0);
	npp.printGrid();
	
	
}

*/


/*
// define regridding method
#define REGRID lterpCube
#define GETVAL getValue


float rhobL = 10; // kg/m3 = gm/m2/mm
float Kdeci = -log(0.05)/(365.25*24*0.25);	  // 2% mass remains after 1 yr
float Kevgr = -log(0.1)/(365.25*24*0.25);	  // 10% mass remains after 1 yr
float ttsL_BL = 0.80, ttsL_NL = 0.4, pi = 3.14159265;  // dxL in mm
float min_viable_fuel_mm = 0; 		// no fire if fuel < 10 mm


// #define L_EVGR_CONDITION (ftmap.values[i] < 20 || ftmap.values[i] > 50)

int calc_pheno(float d, int mode){
	static int month;	// this will hold month number of last month for which calcs were done
	int curr_month = gday2mt(d);

	if (month != curr_month){
		// set current month from date
		month = curr_month;	
		log_fout << "Calculating phenology for date: " << int2str(gday2yr(d)) << "/" << int2str(month) << '\n';
		
		// read npp data
		nppi_handle.readVar_gt(nppi, d, 0);
		REGRID(nppi, npp, npp_indices);

		
		for (int i=0; i<dxL.nlons*dxL.nlats*dxL.nlevs; ++i){
			if ( isEvergreen(ftmap.values[i]) ){
				// in evergreen, all NPP goes into litter!
				dxL.values[i] += npp.values[i]*(24*30.5)/rhobL;
			}
			else{
				// in deciduous lbio accumulates until senescence,
				// and is dumped into dxL in 3 months (NDJ)
				lbio.values[i] += npp.values[i]*(24*30.5);
				if (month == 11 || month == 12){
					float biomass_dropped = lbio.values[i]*0.4;	// drop @ 40%/month.
					lbio.values[i] -= biomass_dropped;
					dxL.values[i] += biomass_dropped/rhobL;
				}
				else if (month == 1){
					dxL.values[i] += lbio.values[i]/rhobL;	// drop all remaining leaves in jan
					lbio.values[i] = 0;
				}
			}
		}
	}

	// decrease dxL via decomposition
	//dxL = dxL - dxL*(Klitter*dt);
	for (int i=0; i < dxL.nlevs*dxL.nlats*dxL.nlons; ++i){
		// K = K0*beta
		float S = lmois.values[i];
		float beta;
		if (S >= 1) beta = 1;
		else beta = 0.25*(1-cos(S*pi))*(1-cos(S*pi));

		if (isEvergreen(ftmap.values[i])) dxL.values[i] -= Kevgr*beta*dxL.values[i]*dt;
		else dxL.values[i] -= Kdeci*beta*dxL.values[i]*dt;
	}
}

int calc_moisture(){
	if (evap.values.size() != mgnlons*mgnlats*mgnlevs) cout << "Error: evap size does not conform.\n";
	//evap.values.resize(mgnlons*mgnlats*mgnlevs);
	for(int i=0; i< evap.values.size(); ++i){
		if( ndr.values[i] == ndr.missing_value || 
			ps.values[i] == ps.missing_value || 
			rh.values[i] == rh.missing_value || 
			ts.values[i] == ts.missing_value || 
			wsp.values[i] == wsp.missing_value ||
			pr.values[i] == pr.missing_value  ){
				evap.values[i] = evap.missing_value;
				lmois.values[i] = lmois.missing_value;
			}
		else{
			float epsilon = 1e-2;
			// EVAP RATE
			float Rn = ndr.values[i]* 0.0864;	// ndr in MJ/day/m2
			float T = ts.values[i] - 273.16;	// ts in degC
			float RH = rh.values[i]/100;	 	// rh (0,1)
			float U = wsp.values[i];			// wsp in m/s
			float Ps = ps.values[i]/1000; 		// ps in kPa
			float ftype = ftmap.values[i];		// forest type code
			float dzL = dxL.values[i];			// dxL in mm (def = 100)
			float fL = lbio.values[i]/(lbio_max.values[i]+epsilon); // fL
			if (fL > 1) fL = 1;

			// amount of Rn actually reaching ground due to shadowing by canopy:
			float Rn0 = Rn*(1 - fA[ftype]*fP[ftype]*fL);	
			
			// float es = 0.6108*exp(17.27*T/(T+237.3));		// kPa, T in degC (tetens)
			// float m = 4098.17*es/(T+237.3)/(T+237.3);		// kPa/degC (tetens)
			float es = 0.13332*exp(21.07-5336/(T+273.16));	// kPa (Merva); 0.1333 converts mmHg to kPa
			float m = 5336*es/(T+273.16)/(T+273.16); 		// kPa/degC (Merva)
			float lv = (2501 - 2.361*T)*1e-3; 			// MJ/Kg
			float y = 0.0016286*Ps/lv; 					// kPa/degC
			float de = es*(1-RH);						// kPa
			float qevap = (m*Rn0 + 6.43*y*(1+0.536*U)*de)/lv/(m+y);  // mm/day
			
			float S = lmois.values[i];
			// EVAP RATE (dependent on lmois)
			float beta;
			if (S >= 1) beta = 1;
			else beta = 0.25*(1-cos(S*pi))*(1-cos(S*pi));
			evap.values[i] = qevap*beta;

			// LITTER MOISTURE
			float qnet_in = pr.values[i] - evap.values[i];
			
			float ttsL = (isNeedleLeaved(ftype))? ttsL_NL:ttsL_BL;

			// if layer saturated, water will drain out
			float qd = (lmois.values[i] >=1 && qnet_in > 0)? qnet_in:0;
			if (dzL > 1) {	// tol 1 mm
				float dSL = ((qnet_in - qd)*dt/24.0)/dzL/ttsL;
				lmois.values[i] += dSL;
			}				
			else lmois.values[i] = 0;
			if (lmois.values[i] < 0) lmois.values[i] = 0;
			if (lmois.values[i] > 1) lmois.values[i] = 1;

			rn.values[i] = Rn0/0.0864;			// output rn in W/m2 only
		}
	}
	lmois.t = ndr.t;
	evap.t = ndr.t;
	rn.t = ndr.t;
}

int calc_fire(double d){
	static int prev_day = 0;
	static int stepnum = -1;
	int curr_day = int(d);

	if (pfire.values.size() != mgnlons*mgnlats*mgnlevs) pfire.values.resize(mgnlons*mgnlats*mgnlevs);
	for(int i=0; i< pfire.values.size(); ++i){
		if( rh.values[i] == rh.missing_value || 
			ts.values[i] == ts.missing_value || 
			wsp.values[i] == wsp.missing_value ||
			lmois.values[i] == lmois.missing_value  ){
				pfire.values[i] = pfire.missing_value;
			}
		else{
			float T = ts.values[i] - 273.16;	// ts in degC
			float RH = rh.values[i]/100;	 	// rh (0,1)
			float U = wsp.values[i];			// wsp in m/s
			float S = lmois.values[i];			// S (0-1)
			float ftype = ftmap.values[i];

			float F = dxL.values[i] - min_viable_fuel_mm;			// fuel content
			if (F < 0) F = 0;
			
//			if (ftype == grass){
//				F += lbio.values[i];
//			}
			
			if (ftype <= 0){
				pfire.values[i] = 0;
			}
			else{
				// p(ts) etc
				float aT = 5, bT  = 30, aRH = 6, bRH = 0.5, aS = 3, bS = 0.3, aU = 1, bU = 1.5;
				float aF = 5, bF = 30;
				float pT = 1/(1+exp(-aT*(T-bT)/bT));
//				float pU = 1/(1+exp(-aU*(U-bU)/bU));
//				float pRH = 1/(1+exp(aRH*(RH-bRH)/bRH));
//				float pS = 1/(1+exp(aS*(S-bS)/bS));
//				float pF = 1/(1+exp(-aF*(F-bF)/bF));
				// p(fire!)
//				pfire.values[i] = pT*pU*pRH*pS*pF;		 // stat Full
//				pfire.values[i] = pT*pow(1-RH,3)*F/100;	 // Karpa
				pfire.values[i] = pT*pow(1-S,2)*F/100;	 // Jaideep \m/
//				float k0 = 1/(5* 20 * (1-0));
//				pfire.values[i] = U * F * (1-S);		 // Fail!

//				// FFWI
//				float m;
//				float H = RH*100;
//				if (H < 10) m = 0.03 + 0.2626*H - 0.00104*H*T;
//				else if (H >= 10 && H < 50) m = 1.76 + 0.1601*H - 0.0266*T;
//				else m = 21.06 - 0.4944*H + 0.005565*H*H - 0.00063*H*T;
//				
//				float eta = 1- 2*(m/30) + 1.5*pow(m/30,2) -0.5*pow(m/30,3);
//				pfire.values[i] = eta*sqrt(1+U*U)*F/100/8;
				
//				pfire.values[i] = (RH*100/20 + (27-T)/10)*F/100/50;	// angstroem index


			}
		}
	}
	pfire.t = ndr.t;
	
	if (curr_day != prev_day){
		// if condition below is just to avoid writing at 0th step
		if (stepnum != -1) {
			double dd = (prev_day - ofire.tbase)*24.0;
			ofire.times.push_back(dd);
			ofire_handle.writeVar(ofire, ofire_NcVar, stepnum);
		}
		prev_day = curr_day; ++stepnum;
		ofire.values = pfire.values;
	}
	else{
		ofire.values = max_vec(ofire.values, pfire.values);
	}
	
}

int getWeatherVars(double d, int mode){
	// read variables from files, mode 0 for "hold"
	rdlwi_handle.readVar_gt(rdlwi, d, 0);
	rdswi_handle.readVar_gt(rdswi, d, 0);
	rulwi_handle.readVar_gt(rulwi, d, 0);
	ruswi_handle.readVar_gt(ruswi, d, 0);
	tsi_handle.readVar_gt(tsi, d, 0);
	psi_handle.readVar_gt(psi, d, 0);
	uwi_handle.readVar_gt(uwi, d, 0);
	vwi_handle.readVar_gt(vwi, d, 0);
	rhi_handle.readVar_gt(rhi, d, 0);
	pri_handle.readVar_gt(pri, d, 0);
	
	// regrid variables onto model grid
	REGRID(rdlwi, rdlw, rdlw_indices);
	REGRID(rdswi, rdsw, rdsw_indices);
	REGRID(rulwi, rulw, rulw_indices);
	REGRID(ruswi, rusw, rusw_indices);
	REGRID(tsi, ts, ts_indices);
	REGRID(psi, ps, ps_indices);
	REGRID(uwi, uw, uw_indices);
	REGRID(vwi, vw, vw_indices);
	REGRID(rhi, rh, rh_indices);
	REGRID(pri, pr, pr_indices);
	//cout << "regriding done\n";

	pr = mask(pr, msk);	// mask out pr.. masking of 1 variable sufficient to mask all.
	
	// calc net downward radiation (W/m2)
	ndr = rdlw + rdsw - rulw - rusw;
	//ndr = rdlw - rulw;	// only longwave contributes to heating... but doesnt work :P
	
	// calc wind speed (m/s)
	wsp = uw*uw + vw*vw;
	wsp.sqrtVar();

}

int writeSinglePointOutput(ofstream &fout){
	fout << rusw.GETVAL(xlon, xlat) << "\t";
	fout << rulw.GETVAL(xlon, xlat) << "\t";
	fout << rdsw.GETVAL(xlon, xlat) << "\t";
	fout << rdlw.GETVAL(xlon, xlat) << "\t";
	fout << uw.GETVAL(xlon, xlat) << "\t";
	fout << vw.GETVAL(xlon, xlat) << "\t";
	fout << ts.GETVAL(xlon, xlat) << "\t";
	fout << ps.GETVAL(xlon, xlat) << "\t";
	fout << rh.GETVAL(xlon, xlat) << "\t";
	fout << pr.GETVAL(xlon, xlat) << "\t";
	fout << npp.GETVAL(xlon, xlat) << "\t";
	//
	fout << ndr.GETVAL(xlon, xlat) << "\t";
	fout << rn.GETVAL(xlon, xlat) << "\t";
	fout << wsp.GETVAL(xlon, xlat) << "\t";
	fout << evap.GETVAL(xlon, xlat) << "\t";
	fout << lmois.GETVAL(xlon, xlat) << "\t";
	fout << lbio.GETVAL(xlon, xlat) << "\t";
	fout << dxL.GETVAL(xlon, xlat) << "\t";
	fout << lbio.GETVAL(xlon, xlat)/lbio_max.GETVAL(xlon, xlat) << "\t";
	fout << pfire.GETVAL(xlon, xlat) << "\t";
	//
	fout << '\n';

}

#define WRITE_REC(x) if (b_##x) x##_handle.writeVar(x, x##_NcVar, istep)

int writeNC(int istep){
	WRITE_REC(rh);
	WRITE_REC(rusw);
	WRITE_REC(rdsw);
	WRITE_REC(rulw);
	WRITE_REC(rdlw);
	WRITE_REC(ps);
	WRITE_REC(ts);
	WRITE_REC(uw);
	WRITE_REC(vw);
	WRITE_REC(pr);
	WRITE_REC(npp);
	//
	WRITE_REC(wsp);
	WRITE_REC(ndr);
	WRITE_REC(rn);
	WRITE_REC(evap);
	WRITE_REC(lmois);
	WRITE_REC(lbio);
	WRITE_REC(dxL);
	WRITE_REC(pfire);
}

int update_ip_files(double d){
	// check if current year has changed.. if so, open next set of files.
	int yr_now = gday2yr(d);
	if (yr_now > curr_yr) {
		curr_yr = yr_now;
		log_fout << "\n~~~~~~~~~~~~~~~~~~~~~~~~ " << curr_yr << " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
		init_ip_vars();
	}

	return 0;
}

int spin(){

	// ~~~~~~~~~~~~~~ DXL SPINUP ~~~~~~~~~~~~~~~~~~
	// spin npp from april to march. sum over all months in 1 year. 
	// all of that carbon will go into litter at start of april. 
	// at spin end, lai(evergreen) = 1, lai(deciduous, grasses) = 0;
	int ix_npp_spin0 = nppi.gday2ix(spin_gday_t0-365.25+15.5); // start 1 yr before lmois spinup
	cout << "npp index0 = " << ix_npp_spin0 << '\n';
	gVar temp; temp.shallowCopy(nppi);
	temp.values.resize(nppi.nlevs*nppi.nlats*nppi.nlons,0);
	for (int i=ix_npp_spin0; i<ix_npp_spin0+12; ++i){
		log_fout << "npp spinup running for ix = " << i << ", date: " << gdatetime(nppi.ix2gday(i)) << '\n';
		nppi_handle.readVar(nppi, i);	// read npp values for month
		nppi = nppi*(24*30.5);		// convert from gm/m2/hr to gm/m2/month
		temp = temp + nppi;	// sum in temp
	}
	REGRID(temp, npp, npp_indices);	// regrid annual sum npp onto model resolution

	// set dxL0
	dxL = npp/rhobL;	// dxL in mm
	// set lbio0 and lbio_max
	lbio_max = npp;		// lbio_max is the npp accumulated over 1 yr (assuming no trends are seen in NPP)
	
	for (int i=0; i<npp.nlons*npp.nlats*npp.nlevs; ++i){
		if (lbio.nlons != npp.nlons || lbio.nlats != npp.nlats){
			cout << "Error: npp and lbio grids dont conform!!\n";
			return 1;
		}
		if (isEvergreen(ftmap.values[i]))
			lbio.values[i] = npp.values[i]; // evergreen
		else lbio.values[i] = 0;			// deciduous
	}

	// ~~~~~~~~~~~~~~ LMOIS SPINUP ~~~~~~~~~~~~~~~~
	int progress_disp_step = nsteps_spin/40;
	if (progress_disp_step == 0) ++progress_disp_step;

	cout << "\n****************************************************************\n\n";
	cout << "Spinup will run for " << nsteps_spin << " steps.\n";
	cout << "progress will be displayed after every " << progress_disp_step << " steps\n";
	if (progress_disp_step == 1) cout << "progress (#steps) >> ";
	else cout << "progress (%) >> "; cout.flush();

	//	double d = gday_t0;
	for (int istep = 0; istep < nsteps_spin; ++istep){
		double d = spin_gday_t0 + istep*(dt/24.0);
		if (d >= gday_t0) {cout << "Error in spin step count! \n"; return 1;}
		
		update_ip_files(d);
				
		fout << gdatetime(d) << "\t";
		getWeatherVars(d,0);	// read weather data from NC files
		calc_pheno(d,0);
		calc_moisture();		// calc evaporation rate (mm/day) and lmois
		writeSinglePointOutput(fout);	// write single point output to file

		if (istep % progress_disp_step == 0) {cout << "."; cout.flush();}
	}	

	fout << '\n';
	write_endstate("spin_end");
	cout << " > 100%\n";
	
	return 0;
}


int init_infisim(){
	log_fout.open("../output/log.txt");	// open log stream
	log_fout << " ******************* THIS IS LOG FILE ****************************\n\n";
	log_fout << "~~~~~~~~~~~~~~~~~~~~~~~~ " << curr_yr << " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";

	cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
	cout << "~                  I N F I S I M                               ~\n";
	cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
	cout << "\n> Reading config parameters... "; cout.flush();
	read_sim_config();
	cout << "DONE.\n> Reading input filenames... "; cout.flush();
	read_ip_params();
	cout << "DONE.\n> Reading forest type params... "; cout.flush();
	read_ft_params();
	cout << "DONE.\n> Initialising model vars/files... "; cout.flush();
	init_modelvars();
	cout << "DONE.\n> Initialising input vars/files... "; cout.flush();
	init_ip_vars();	// open 1st set of files and init input variables
	cout << "DONE.\n> Setting up initial condition (lmois0)... "; cout.flush();
	init_ic();
	cout << "DONE.\n";

	cout << "Opening file to write point values: " << pointOutFile << '\n';
	fout.open(pointOutFile.c_str());

	fout << "lat:\t " << xlat << "\t lon:\t" << xlon << '\n';
	fout << "ftype: \t" << ftmap.getCellValue(xlon, xlat) << '\n';
	fout << "datetime\trusw\trulw\trdsw\trdlw\tuw\tvw\tts\tps\trh\tpr\tnpp\t" <<
			"ndr\trn\twsp\tevap\tlmois\tlbio\tdxL\tLAI\tFire\n";

	return 0;
}
*/
/*
int main(){

	// set NETCDF error behavior to non-fatal
	NcError err(NcError::silent_nonfatal);

	// init
	init_infisim();

	if (lspinup) spin();

	int progress_disp_step = nsteps/40;
	if (progress_disp_step == 0) ++progress_disp_step;

	cout << "\n****************************************************************\n\n";
	cout << "Simulation will run for " << nsteps << " steps.\n";
	cout << "progress will be displayed after every " << progress_disp_step << " steps\n";
	if (progress_disp_step == 1) cout << "progress (#steps) >> ";
	else cout << "progress (%) >> "; cout.flush();

	//	double d = gday_t0;
	for (int istep = 0; istep < nsteps; ++istep){
		double d = gday_t0 + istep*(dt/24.0);
		if (d > gday_tf) {cout << "Error in sim step count! \n"; return 1;}
		
		update_ip_files(d);
				
		// read weather data from NC files
		getWeatherVars(d,0);

		// run model components
		calc_pheno(d,0); 	// calc dxL, fL and canopy biomass
		calc_moisture(); 	// calc evaporation rate (mm/day) and lmois
		calc_fire(d); 		// calc fire probability

		// write single point output to file
		fout << gdatetime(d) << "\t";
		writeSinglePointOutput(fout); 
		
		// Write desired output variables to nc files
		writeNC(istep);	
		
		if (istep % progress_disp_step == 0) {cout << "."; cout.flush();}
	}	

	cout << " > 100%\n";
	cout << "\n****************************************************************\n";
	write_state("last");
	write_state(dxL, "last");
	write_state(lmois, "last");

	ofire.ntimes = ofire.times.size();
	ofire_handle.writeTimeValues(ofire);	// ofire time values were never written!

	// destructor will close all the NC files :) ...phew!
	cout << "> Check variable details in log file /output/log.txt\n";
	cout << "> Successfully wrote NC files.\n\n\n";
}
*/
