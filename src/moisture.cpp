using namespace std;
#include "../include/moisture.h"
#include "../include/globals.h"
#include "../include/vars.h"


//float rhobL = 10; 				// kg/m3 = gm/m2/mm
//float theta_sL = 0.8; 			// saturation water content of litter

//ofstream fout_mois("../output/mois_tease.txt");


int calc_ndr(double gt){ 

	int m = gt2month(gt);

	float N = gt2daynum(gt);
	float hr = (gt - int(gt))*24;
	//cout << "t = " << t << ", N = " << N << "\n";

	float y = 2*pi/365*(N-1 + (hr-12)/24);

	float eqtime = 229.18*(0.000075 + 0.001868*cos(y) - 0.032077*sin(y)
				 - 0.014615*cos(2*y) - 0.040849*sin(2*y));
	
	float decl = 0.006918 - 0.399912*cos(y) + 0.070257*sin(y) - 0.006758*cos(2*y)
			   + 0.000907*sin(2*y) - 0.002697*cos(3*y) + 0.00148*sin(3*y);

	for (int ilat=0; ilat<mgnlats; ++ilat){
		for (int ilon=0; ilon<mgnlons; ++ilon){

			if (pr(ilon,ilat,0) == pr.missing_value){
				ndr(ilon,ilat,0) = ndr.missing_value;
			}
			else
			{

				float lat = mglats[ilat]*pi/180;
				float lon = mglons[ilon];

				float time_offset = -(eqtime - 4*lon); // + 60*5.5;	// 5.5 is timezone
				float tst = hr*60 + time_offset;
				float ha = (tst/4 - 180)*pi/180;
				float sinY = sin(lat)*sin(decl)+ cos(lat)*cos(decl)*cos(ha);

				float SW_d;
				if (sinY > 0) {
					SW_d = 1370*(0.6+0.2*sinY)*sinY;
				}
				else SW_d = 0;
			
				float cld = pr(ilon,ilat,0)/3e-4/86400;
				if (cld < 0) cld = 0; if (cld > 1) cld = 1;

				ndr(ilon, ilat, 0) = SW_d*(1-albedo(ilon,ilat,0)/100)*(1-cld*0.6) - 100; // 100 is lw_in - lw_out
			}

		}
	}


}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	--> calc_moisture(...)
	
	calculate LAI
	calculate radiation reaching the ground
	calculate evaporation rate
	solve bucket model to calculate moisture content
		
	Input: all weather variables
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
int calc_moisture(){ 

	for (int ilat=0; ilat<mgnlats; ++ilat){
		for (int ilon=0; ilon<mgnlons; ++ilon){

			if( ndr(ilon,ilat,0) == ndr.missing_value 	|| 
				ps(ilon,ilat,0) == ps.missing_value 	|| 
				rh(ilon,ilat,0) == rh.missing_value 	|| 
				ts(ilon,ilat,0) == ts.missing_value 	|| 
				wsp(ilon,ilat,0) == wsp.missing_value 	||
				pr(ilon,ilat,0) == pr.missing_value  	||
				lmois(ilon,ilat,0) == lmois.missing_value  ||
				cmois(ilon,ilat,0) == cmois.missing_value  )
			{
				lmois(ilon,ilat,0) = lmois.missing_value;
				cmois(ilon,ilat,0) = cmois.missing_value;
			}
			else{
			
				// get all variables in the right units
				float Rn  = ndr(ilon,ilat,0)*0.0864;	// convert Rn from W/m2 to MJ/day/m2
				float T   = ts(ilon,ilat,0) - 273.16;	// ts in degC
				float RH  = rh(ilon,ilat,0)/100;	 	// rh (0-1)
				float U   = wsp(ilon,ilat,0);			// wsp in m/s
				float Ps  = ps(ilon,ilat,0)/1000; 		// ps in kPa
				float dzL = dxl(ilon,ilat,0)/rhobL;		// dzL in mm (dxl is in gm/m2)
				float pre = pr(ilon,ilat,0);			// pr in mm/day

				// calculate PFT independent values needed for PER
				// float es = 0.6108*exp(17.27*T/(T+237.3));	// kPa, T in degC (tetens)
				// float m = 4098.17*es/(T+237.3)/(T+237.3);	// kPa/degC (tetens)
				float es = 0.13332*exp(21.07-5336/(T+273.16));	// kPa (Merva); 0.1333 converts mmHg to kPa
				float m = 5336*es/(T+273.16)/(T+273.16); 		// kPa/degC (Merva) 
				float lv = (2501 - 2.361*T)*1e-3; 				// MJ/Kg - latent heat of vapourization
				float y = 0.0016286*Ps/lv; 						// kPa/degC -  slope of vapour pressure curve
				float de = es*(1-RH);							// kPa	- vapour pressure deficit


				// calculate LAI and Radiation reaching ground (Rl)
				float Wc_sat = 0;
//				float Wc_sat_vec[] = {0,0,0.5,0.5,0.5,0.5,0.5};
				float Rl = 0, Rc = 0, lai = 0;
				for (int i=0; i<npft; ++i){

					// calculate LAI = canbio/canbiomax * LAImax
					lai = LAImax[i]*canbio(ilon,ilat,i)/(canbio_max(ilon,ilat,i)+1e-6);	// avoid NaN
					
					Wc_sat += Wc_sat_vec[i]*lai*vegtype(ilon,ilat,0);

					// R = Watts reaching canopy of i'th PFT, will be intercepted by canopy
					//   = R(W/m2)* f * acell(m2)
					// then add all the R's and divide by acell to get avg Rn (W/m2) 
					// since acell appears in both num & den, set it to 1 (normalized)
					float R_i = Rn* vegtype(ilon,ilat,i);	

					// Radiation reaching ground from i'th PFT's canopy
					// beer-lambert law REF: ffmodel_2003_pual_swuf (they have used 0.5)
					float Rl_i = R_i *exp(-0.4*lai);	
					Rl += Rl_i;	// add to total.

				}
				// calculate radiation absorbed by canopy
				// canopy for all PFT's is treated alike
				Rc = Rn - Rl;	// radiation intercepted by canopy is that which does not reach litter layer

				// calculate potential evaporation rates from canopy and litter
				float Ep_l = (m*Rl + 6.43*y*(1+0.536*U)*de)/lv/(m+y);  // mm/day - potential evap rate (litter)
				float Ep_c = (m*Rc + 6.43*y*(1+0.536*U)*de)/lv/(m+y);  // mm/day - potential evap rate (canopy)

				// calculate Precipitation reaching ground = all pr where PFT is X, 0 otherwise
				// where forest is present, all pr is intercepted, then drained as per the water content
				float pr_l = vegtype(ilon,ilat,0)*pre;	// pr reaching litter = f_barren * pr
				float pr_c = pre - pr_l;				// pr intercepted by canopy
			
				// Water balance for canopy
				float Wc = cmois(ilon, ilat, 0);		// mm		
				float beta_c = Wc/Wc_sat;				// Wc_sat (mm = kg/m2)		
				Wc += (pr_c - Ep_c*beta_c)/24.0f*dt;	// mm/day * days/hr *dt = mm/hr*dt = mm (in this timestep)
				
				float Dc = 0; // drain from canopy in (mm/day) 
				if (Wc >= Wc_sat){
					Dc = (Wc - Wc_sat)*24.0/dt;
					Wc = Wc_sat;
				}
				if (Wc < 0) Wc = 0;
				
				// water balance for litter (depends on litter thickness)
				float S = lmois(ilon,ilat,0);
				float beta_l = 0.25*(1-cos(S*pi))*(1-cos(S*pi)); 
				float qnet_in = (pr_l + Dc - Ep_l*beta_l)/24.0f;	// mm/hr
				float q_drain = (S >=1 && qnet_in > 0)? qnet_in:0; // when layer becomes saturated, all excess water should drain out
				// solve for litter moisture
				if (dzL > 0.5) {	// use this tolerance to avoid division by zero. If litter layer is too thin, it will saturate instantly
					S += (qnet_in - q_drain)/dzL/theta_sL *dt;
				}				
				else {
					S = 0;
				}
				if (S < 0) S = 0;
				if (S > 1) S = 1;
				
				// update lmois and cmois
				lmois(ilon,ilat,0) = S;
				cmois(ilon,ilat,0) = Wc;
				//evap(ilon, ilat, 0) = Ep_l + Ep_c;

				// SP output
				if (ilon == i_xlon && ilat == i_xlat && spout_on){
					//cout << "Dc = " << Dc << ", Wc = " << Wc << '\n';
					sp_fout << Rn << "\t" << lai << "\t" << Rc << "\t" << Rl << "\t" 
							<< pre << "\t"  << pr_c << "\t"  << pr_l << "\t" 
							<< Ep_c << "\t"  << Ep_l << "\t"  << Ep_c*beta_c << "\t" << Ep_l*beta_l << "\t" 
							<< dzL << "\t" << S << "\t" << Wc << "\t";
					//sp_fout << "\n";					
				}
				
			}

		}	// ilon loop ends
	}	// ilat loop ends

	lmois.t = ndr.t;
	cmois.t = ndr.t;
}



