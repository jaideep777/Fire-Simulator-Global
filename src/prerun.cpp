#include "../include/prerun.h"
#include "../include/io.h"
#include "../include/pheno.h"
#include "../include/moisture.h"
#include "../include/fire.h"
using namespace std;


char ps2char(int c){
	     if (c == psX) return 'X';
	else if (c == psF) return 'F';
	else if (c == psM) return 'M';
	else if (c == psS) return 'S';
	else if (c == psZ) return 'Z';
	else if (c == psE) return 'E';
	else return 'X';
}


inline int printRunHeader(string s, int ns, int ds){
	cout << "\n****************************************************************\n\n";
	cout << s + " will run for " << ns << " steps.\n";
	cout << "progress will be displayed after every " << ds << " steps\n";
	if (ds == 1) cout << "progress (#steps) >> ";
	else cout << "progress (%) >> "; cout.flush();
}

extern map<string, ip_data> ip_data_map;

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	--> prerun_canbio_ic()
	
	Pre-run to set C0, the initial canopy and litter biomass. On yearly basis, 
		dL/dt = Mshed - kL
		dC/dt = Mnpp  - Mshed
	over the years, the canopy is in equilibrium, so
		Mshed = Mnpp, L = Mshed/k
	
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
int prerun_canbio_ic(){

	int nyrs = 5;	// number of years of bio pre-run
	double gt0 = spin_gday_t0 - 365.2524*nyrs; // start 2n yrs before lmois spinup
	if (gt2day(gt0) < 10 || gt2day(gt0) > 20) gt0 += 15.5;	// bring day to centre of month
	int ix_npp_spin0 = nppi.gt2ix(gt0);
	//cout << ix_npp_spin0 << " " << gt2string(gt0) << " " << gt2day(gt0) << endl; 
	
	int nsteps_npp = nyrs*12;
	int dstep_npp = nsteps_npp/40+1;

	printRunHeader("Canopy Biomass pre-run", nsteps_npp, dstep_npp);

	canbio.fill(0);
	canbio_max.fill(0);
	dxl.fill(0);
	lmois.fill(0);

	for (int istep = ix_npp_spin0 ; istep < nsteps_npp + ix_npp_spin0; ++istep){
		double d = nppi.ix2gt(istep);	
		
		// update input files
		update_ip_files(d);

		// read NPP values
		nppi.ifile_handle->readVar(nppi, istep);	// get total NPP value 
		PREPROC_REGRID(nppi, npp, ip_data_map.find(nppi.varname)->second.lterp_indices);

		// calculate canbio and littbio
		calc_pheno(d, hrsPerMonth); 	// delT is 1 month

		for (int i=0; i<canbio.values.size(); ++i){
			canbio_max[i] = fmax(canbio_max[i], canbio[i]);
		}
		canbio_max.t = d;

		if (spout_on) sp_fout << "\n";
		if (istep % dstep_npp == 0) {cout << "."; cout.flush();}
	}	

	write_state(canbio_max);
	canbio_prerun_on = false;	// finished prerun, so set flag to false
	
	cout << " > 100%\n";
	
}




/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	--> prerun_lmois_ic()
	
	Pre-run to stabilize litter moisture. It should start in April 
	with 0 litter moisture.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
int prerun_lmois_ic(){

	int dstep_spin = nsteps_spin/40+1;

	printRunHeader("Litter Moisture pre-run", nsteps_spin, dstep_spin);

	for (int istep = 0; istep < nsteps_spin; ++istep){
		double d = spin_gday_t0 + istep*(dt/24.0);
		if (d >= gday_t0) {cout << "Error in spin step count! \n"; return 1;}
		
		// update input files
		update_ip_files(d);
				
		// read input (forcing) data from NC files
		int k = read_all_ip_vars(d, 0);

		// run model components
		calc_pheno(d, dt);
		calc_ndr(d);
		calc_moisture();		// calc evaporation rate (mm/day) and lmois

		if (spout_on) sp_fout << "\n";
		if (istep % dstep_spin == 0) {cout << "."; cout.flush();}
	}	

	cout << " > 100%\n";
}


/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	--> main_run()
	
	Main Run!
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
int main_run(){
	// simulation header
	printRunHeader("MAIN RUN", nsteps, dstep);

	// mainloop
	for (int istep = 0; istep < nsteps; ++istep){
		double gtnow = gday_t0 + istep*(dt/24.0);
		if (gtnow > gday_tf) {cout << "Error in sim step count! \n"; return 1;}
		
		// update input files (if time has gone beyond the limits of any of the current files)
		update_ip_files(gtnow);
				
		// read input (forcing) data from NC files
		int k = read_all_ip_vars(gtnow, 0);

		// run model components
		calc_pheno(gtnow, dt); 	// calc dxL, and canopy biomass
		calc_ndr(gtnow);
		calc_moisture(); 		// calc lmois
		calc_fire(gtnow); 		// calc fire index

		// Write desired output variables to nc files and to singlePointOut
		write_all_outvars(istep);	

		if (spout_on) sp_fout << "\n";
		if (istep % dstep == 0) {cout << "."; cout.flush();}
	}	

	cout << " > 100%\n";
	cout << "\n****************************************************************\n";

}



