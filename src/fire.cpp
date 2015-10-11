using namespace std;
#include "../include/moisture.h"
#include "../include/globals.h"
#include "../include/vars.h"
#include <cmath> 

const float min_viable_fuel = 0.5;	// min mm of fuel required for fire 

ofstream fout_fire("../output/fireIndex.txt");

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	--> calc_fire(...)
	
	calculate various fire indices 
	calculate daily index as maximum throughout the day
		
	Input: all weather variables
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
int calc_fire(double d){ 

	static int prev_day = 0;
	static int stepnum = -1;
	int curr_day = int(d);

	for (int ilat=0; ilat<mgnlats; ++ilat){
		for (int ilon=0; ilon<mgnlons; ++ilon){

			if( ts(ilon,ilat,0) == ts.missing_value || 
				lmois(ilon,ilat,0) == lmois.missing_value ||
				dxl(ilon,ilat,0) == dxl.missing_value  )
			{
				fire(ilon,ilat,0) = fire.missing_value;
			}
			else{
			
				float T = ts(ilon,ilat,0) - 273.16;	// ts in degC
//				float RH = rh.values[i]/100;	 	// rh (0,1)
//				float U = wsp.values[i];			// wsp in m/s
				float S = lmois(ilon,ilat,0);			// S (0-1)
				float forest_frac = 1 - vegtype(ilon, ilat, 0) - vegtype(ilon, ilat, 1);	// exclude X and AGR from forest types

				float F = dxl(ilon, ilat, 0) - min_viable_fuel;			// fuel content
				if (F < 0) F = 0;
			
//				// J-index
//				float gT = pow(T/20, 4);  if (T<0) gT = 0;						// J1
				float gT = pow(T/15,3)-pow(T/21,4);  if (gT < 0) gT = 0;		// J2
//				fire(ilon, ilat, 0) = gT * 25*pow(1-S, 2)*forest_frac/10000;	// J1,J2
//				fire(ilon, ilat, 0) = gT * 25*pow(1-S, 2)*float(forest_frac>0.3)*dff/10000;	// J3
				fire(ilon, ilat, 0) = gT * 25*pow(1-S, 1)*pow(F,0.3)*float(forest_frac>0.3)/10000;	// J4

				// p(ts) etc
//				float aT = 5, bT  = 30, aRH = 6, bRH = 0.5, aS = 3, bS = 0.3, aU = 1, bU = 1.5;
//				float aF = 5, bF = 30;
//				float pT = 1/(1+exp(-aT*(T-bT)/bT));
//				float pU = 1/(1+exp(-aU*(U-bU)/bU));
//				float pRH = 1/(1+exp(aRH*(RH-bRH)/bRH));
//				float pS = 1/(1+exp(aS*(S-bS)/bS));
//				float pF = 1/(1+exp(-aF*(F-bF)/bF));

				// p(fire!)
//				pfire.values[i] = pT*pU*pRH*pS*pF;		 // stat Full
//				pfire.values[i] = pT*pow(1-RH,3)*F/100;	 // Karpa
//				fire(ilon, ilat, 0) = ;	 // Jaideep \m/
//				float k0 = 1/(5* 20 * (1-0));
//				pfire.values[i] = U * F * (1-S);		 // Fail!

//				// FFWI Index
//				float m;
//				float H = RH*100;
//				if (H < 10) m = 0.03 + 0.2626*H - 0.00104*H*T;
//				else if (H >= 10 && H < 50) m = 1.76 + 0.1601*H - 0.0266*T;
//				else m = 21.06 - 0.4944*H + 0.005565*H*H - 0.00063*H*T;
//				float eta = 1- 2*(m/30) + 1.5*pow(m/30,2) -0.5*pow(m/30,3);
//				pfire.values[i] = eta*sqrt(1+U*U)*F/100/8;
			
				// Angstroem index
//				pfire.values[i] = (RH*100/20 + (27-T)/10)*F/100/50;	

				fout_fire << gt2string(d) << '\t' << T << '\t' << S << '\t' << F << '\t' << forest_frac << '\t' << fire(ilon, ilat, 0) << '\n';
				
			}

		}	// ilon loop ends
	}	// ilat loop ends
	fire.t = ts.t;
	dfire.t = ts.t;

	int hr = (d - int(d))*24; //cout <<  "hr = " << hr << '\n';
	if (hr == 6){
		dfire.values = fire.values;
//		cout <<  "writing\n";
		if (stepnum != -1) dfire.ofile_handle->writeVar(dfire, dfire.outNcVar, stepnum);
		++stepnum;
	}
	
//	if (curr_day != prev_day){
//		if (stepnum != -1) { // avoid writing at 0th step
//			dfire.ofile_handle->writeVar(dfire, dfire.outNcVar, stepnum);
//		}
//		prev_day = curr_day; ++stepnum;
//		dfire.values = fire.values;
//	}
//	else{
//		// set max of current and previous steps within a day in dfire
//		for (int i=0; i<dfire.values.size(); ++i){
//			if (dfire[i] != dfire.missing_value ||
//				fire[i] != fire.missing_value )
//			{
//				dfire[i] = fmax(dfire[i], fire[i]);
//			}
//			else{
//				dfire[i] = dfire.missing_value;
//			}
//		}
//	}
	
	if (spout_on){
		sp_fout << fire(i_xlon, i_xlat, 0) << "\t" << dfire(i_xlon, i_xlat, 0) << "\t";
	}
	
}


