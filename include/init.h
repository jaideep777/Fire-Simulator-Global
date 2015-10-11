#ifndef INIT_H
#define INIT_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
using namespace std;

#include <gsm.h>

// global variables for using in init.cpp only
//const string attrbegin = ">";
//int l_ip_init_done = false;
enum ipvar_value_type {ivt_avg, ivt_inst, ivt_sum};
enum ipvar_lterp_mode {ilm_hold, ilm_lin};

// this class stores all meta-info of input variables
class ip_data {
public:
	// data supplied from ip-params file
	string 				name;					// variable name
	string 				unit;					// unit
	string 				fname_prefix;			// prefix in filename
	int 				start_yr, nyrs_file;	// start yr in file, #yrs in file
	ipvar_value_type 	vt;						// value type (inst, avg etc)
	ipvar_lterp_mode 	lm;						// interpolation mode (hold, linear-interp etc)
	
	// final form data generated during init
	vector <int> 		lterp_indices;			// interpolation indices (generated during init)
	vector <string>		fnames;					// list of filenames (genrated during init)
	
	ip_data(string _n, string _u, string _fnp, int _sy, int _ny, 
			string _vt, string _lm); 
	void print_vardata(ofstream &fout1);
};

// READ PARAMS FILES
int read_ip_params_file();
int read_sim_config_file();
int read_veg_params_file();

int init_modelvar(string var_name, string unit, int nl, gVar &v, vector <double> times_vec, ostream& lfout);
int init_ip_var(string var_name, gVar &v, gVar &vout, ostream &lfout);
int read_static_var(gVar &v, ostream &lfout);				

#define PREPROC_INIT_IP_VAR(x) if (varUse_map[#x]) init_ip_var( #x, x##i, x, log_fout) 
#define PREPROC_INIT_MODELVAR(x, unit, nl) if (varUse_map[#x]) init_modelvar(#x, unit, nl, x, mgtimes, log_fout)
#define PREPROC_READ_STATIC_VAR(x) read_static_var(x, log_fout)

int init_all_ip_vars();
int init_all_modelvars();
int read_all_static_vars();

int init_infisim();
int init_ic();


#endif

