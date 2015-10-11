#include <iostream>
#include <math.h>
#include <netcdfcpp.h>
#include <fstream>
#include <vector>
#include <map>
using namespace std;

#include <gsm.h>
#include "../include/globals.h"
#include "../include/vars.h"
#include "../include/init.h"
#include "../include/io.h"

int update_ip_files(double gtime){
	// check if current year has changed.. if so, open next set of ip files.
	int yr_now = gt2year(gtime);
	if (yr_now > ip_curr_yr) {
		ip_curr_yr = yr_now;
		logdebug << "Updating input files!\n";
		log_fout << "\n~~~~~~~~~~~~~~~~~~~~~~~~ " << ip_curr_yr << " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
		init_all_ip_vars();
		return 1;
	}
	return 0;
}

extern map<string, ip_data> ip_data_map; 

int read_ip_var(gVar &vi, gVar &v, double gtime, int mode){
	// read into input variable vi
	int k = vi.ifile_handle->readVar_gt(vi, gtime, mode);
	// interpolate and store into model variable v
	PREPROC_REGRID(vi, v, ip_data_map.find(vi.varname)->second.lterp_indices );
	return k;
}


int write_outvar(gVar &v, int istep){
	// write output nc-file record 
	if (v.lwrite) v.ofile_handle->writeVar(v, v.outNcVar, istep);
	// write single point output record
	if (v.lwriteSP) sp_fout << v.PREPROC_GETVAL(xlon, xlat) << "\t";
}


int read_all_ip_vars(double gtime, int mode){
	// read variables from files, mode 0 for "hold", 1 for "interpolate"
	int k; // return value will be stored in this variable

	PREPROC_READ_IP_VAR(pr, k);
//	PREPROC_READ_IP_VAR(ps, k);
	PREPROC_READ_IP_VAR(rh, k);
	PREPROC_READ_IP_VAR(ts, k);
//	PREPROC_READ_IP_VAR(ndr, k);
	PREPROC_READ_IP_VAR(wsp, k);
	PREPROC_READ_IP_VAR(npp, k);

	pr = mask(pr, msk);	// mask out pr.. masking of 1 variable sufficient to mask all.
	return k;
}

int write_all_outvars(int istep){
//	if (spout_on){
//		double gtnow = gday_t0 + istep*(dt/24.0);
//		sp_fout << gt2string(gtnow) << "\t"; 
//	}
	
	PREPROC_WRITE_OUTVAR(pr);
	PREPROC_WRITE_OUTVAR(rh);
	PREPROC_WRITE_OUTVAR(ts);
	PREPROC_WRITE_OUTVAR(wsp);
	PREPROC_WRITE_OUTVAR(npp);
	
	PREPROC_WRITE_OUTVAR(ps);
	PREPROC_WRITE_OUTVAR(ndr);
	
	PREPROC_WRITE_OUTVAR(canbio);
	PREPROC_WRITE_OUTVAR(dxl);
	PREPROC_WRITE_OUTVAR(lmois);
	PREPROC_WRITE_OUTVAR(cmois);
	PREPROC_WRITE_OUTVAR(fire);
	PREPROC_WRITE_OUTVAR(evap);
	

//	if (spout_on) sp_fout << "\n";
}



int write_state(gVar &v, string suffix){	// default suffix = ""

	// get current time in gVar's units
	double x = (v.t - v.tbase)*24.0/v.tscale;	// current time in gVar's units

	// back-up gVar's original time vector
	vector <double> times_orig = v.times;
	if (v.ntimes != v.times.size()) cout << "ERROR!! variable ("+ v.varname +") has inconsistent time and ntimes!\n\n";

	// set new time vector which is just the current time
	v.ntimes = 1; v.times = vector <double> (1, x);	 

	// prepare filename
	if (suffix != "") suffix = "_" + suffix;
	string ncoutfile = "../output/" + v.varname + suffix + "." + gday2ymd(v.t) +"_"+ xhrs2hms(v.t - int(v.t))+".nc";

	// create a new file_handle so that gVar's original handle is not disturbed
	NcFile_handle v_handle;
	v_handle.open(ncoutfile, "w", glimits_india);
	
	// write gVar from created handle
	v_handle.writeCoords(v);
	v_handle.writeTimeValues(v);
	NcVar * vVar = v_handle.createVar(v);
	v_handle.writeVar(v, vVar, 0);	// write at tindex=0 of course!
	v_handle.close();

	// restore original gVar's time vector
	v.times = times_orig;	
	v.ntimes = times_orig.size();
	
	return 0;
}

