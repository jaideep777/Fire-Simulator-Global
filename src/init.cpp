#include "../include/init.h"
#include "../include/globals.h"
#include "../include/vars.h"
using namespace std;

// global variables for use in this file only
const string attrbegin = ">";
static int l_ip_init_done = false;
map <string, ip_data> ip_data_map;	// map from var name to ip data
map <string, string> data_dirs;		// map from var name to data dir
map <string, int> writeVar_flags;	// map from var name to output flag
map <string, int> writeVarSP_flags;	// map from var name to SP-output flag

static string params_dir = "../params";
static string params_ft_file = params_dir + "/" + "params_ft.r";
static string params_ip_file = params_dir + "/" + "params_ip.r";
static string sim_config_file = params_dir + "/" + "sim_config.r";


// Class ip_data
ip_data::ip_data(string _n, string _u, string _fnp, int _sy, int _ny, 
				string _vt, string _lm) : 
				name(_n), unit(_u), fname_prefix(_fnp), 
				start_yr(_sy), nyrs_file(_ny) {
		// 
		if (_vt == "avg") vt = ivt_avg;
		else if (_vt == "inst") vt = ivt_inst;			
		else if (_vt == "sum") vt = ivt_sum;			

		if (_lm == "hold") lm = ilm_hold;		
		else if (_lm == "lin") lm = ilm_lin;
		
		fnames.resize(0);
		lterp_indices.resize(0);		
}
		
void ip_data::print_vardata(ofstream &fout1){
		fout1 << name << " " << unit << " " << fname_prefix << " " << start_yr << " " << nyrs_file << "\n"; 
}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	--> read_ip_params_file()

	READ INPUT PARAMS FILE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
int read_ip_params_file(){
	ifstream fin;
	fin.open(params_ip_file.c_str());
	
	string s, u, v, w, y;
	int n, m, l;
	float f;
	
	while (fin >> s && s != attrbegin);	// read until 1st > is reached
	
	fin >> s; 
	if (s != "INPUT_DATA_TIME_BOUNDS") {cout << "time bounds not found!"; return 1;}
	while (fin >> s && s != attrbegin){
		if (s == "") continue;	// skip empty lines
		if (s == "#") {getline(fin,s,'\n'); continue;}	// skip #followed lines (comments)
		fin >> n;
			 if (s == "ip_start_yr") ip_start_yr = ip_curr_yr = n;
		else if (s == "ip_end_yr") 	 ip_end_yr = n;
	}

	fin >> s; 
	if (s != "FORCING_DATA_DIRS") {cout << "ip data dirs not found!"; return 1;}
	while (fin >> s && s != attrbegin){
		if (s == "") continue;
		if (s == "#") {getline(fin,s,'\n'); continue;}	// skip #followed lines (comments)
		fin >> u;
		
		data_dirs[s] = u;
		logdebug << s << ": " << data_dirs[s] << ".\n";	
	}
	
	fin >> s; 
	if (s != "FORCING_VARIABLE_DATA") {cout << "variable data not found!"; return 1;}
	while (fin >> s && s != attrbegin){
		if (s == "") continue;	// skip empty lines
		if (s == "#") {getline(fin,s,'\n'); continue;}	// skip # following stuff (comments)
		fin >> u >> v >> m >> n >> w >> y;
		//ip_start_yrs[s] = m;	// set start yr of input data for variable s

		// ip_var_names.push_back(s);
		ip_data a(s, u, v, m, n, w, y);
		ip_data_map.insert( pair <string, ip_data> (s,a) );
		logdebug << s << ": "; ip_data_map.find(s)->second.print_vardata(log_fout);
	}	
		
	fin >> s; 
	if (s != "STATIC_INPUT_FILES") {cout << "static input files not found!"; return 1;}
	while (fin >> s && s != attrbegin){
		if (s == "") continue;	// skip empty lines
		if (s == "#") {getline(fin,s,'\n'); continue;}	// skip #followed lines (comments)
		fin >> u;
		if (s == "msk") 			msk.ifname = "../input/" + u;
		else if (s == "lmois") 		lmois.ifname = "../input/" + u;	
		else if (s == "vegtype") 	vegtype.ifname = "../input/" + u;
		else if (s == "albedo") 	albedo.ifname = "../input/" + u;
		else if (s == "elev") 		elev.ifname = "../input/" + u;
	}
	
	fin.close();
}


/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	--> read_sim_config_file()

	READ SIMULATION CONFIGURATION PARAMETERS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
int read_sim_config_file(){
	ifstream fin;
	fin.open(sim_config_file.c_str());

	//const string attrbegin = ">";

	string s, u, v, w;
	int n, m, l;
	float f;
	while (fin >> s && s != attrbegin);	// read until 1st > is reached
	
	fin >> s; 
	if (s != "TIME") {cout << "sim time not found!"; return 1;}
	while (fin >> s && s != attrbegin){
		if (s == "") continue;	// skip empty lines
		if (s == "#") {getline(fin,s,'\n'); continue;}	// skip #followed lines (comments)
		fin >> u;
		if		(s == "sim_date0")	sim_date0 = u;
		else if (s == "sim_t0")		sim_t0 = u;
		else if (s == "sim_datef")	sim_datef = u;
		else if (s == "sim_tf")		sim_tf = u;
		else if (s == "dt")			dt = str2float(u);
		else if (s == "base_date")	gday_tb = ymd2gday(u);
		else if (s == "spinup")		lspinup = (u == "on")? true:false;
		else if (s == "spin_date0")	spin_gday_t0 = ymd2gday(u);	// assume time = 0:0:0
	}

	fin >> s; 
	if (s != "MODEL_GRID") {cout << "model grid not found!"; return 1;}
	while (fin >> s && s != attrbegin){
		if (s == "") continue;	// skip empty lines
		if (s == "#") {getline(fin,s,'\n'); continue;}	// skip #followed lines (comments)
		fin >> u;
		if		(s == "lon0")	mglon0 = str2float(u);
		else if (s == "lonf")	mglonf = str2float(u);
		else if (s == "lat0")	mglat0 = str2float(u);
		else if (s == "latf")	mglatf = str2float(u);
		else if (s == "dlat")	mgdlat = str2float(u);
		else if (s == "dlon")	mgdlon = str2float(u);
		else if (s == "xlon")	  xlon = str2float(u);
		else if (s == "xlat")	  xlat = str2float(u);
		else if (s == "pointOutFile") pointOutFile = u;
		else if (s == "SPout_on") spout_on = (u == "1")? true:false;
	}
	
	fin >> s; 
	if (s != "OUTPUT_VARIABLES") {cout << "output variables not found!"; return 1;}
	while (fin >> s && s != attrbegin){
		if (s == "") continue;	// skip empty lines
		if (s == "#") {getline(fin,s,'\n'); continue;}	// skip #followed lines (comments)
		fin >> n >> m;

		writeVar_flags[s] = n;
		writeVarSP_flags[s] = m;
		logdebug << "Write flag for " << s << ": " << writeVar_flags[s] << ", " << writeVarSP_flags[s] << '\n';
	}

	fin >> s; 
	if (s != "VARS_TO_USE") {cout << "vars to use (debugging) not found!"; return 1;}
	while (fin >> s && s != attrbegin){
		if (s == "") continue;	// skip empty lines
		if (s == "#") {getline(fin,s,'\n'); continue;}	// skip #followed lines (comments)
		fin >> n;

		varUse_map[s] = (n==1)? true:false;
		if (n==0) writeVar_flags[s] = 0;	// if variable is not used at all, it cant be output!
		if (n==0) writeVarSP_flags[s] = 0;	// if variable is not used at all, it cant be output!
		logdebug << "Use variable flag for " << s << ": " << varUse_map[s] << '\n';
	}

	// Set Sim Date and Time
	gday_t0 = ymd2gday(sim_date0) + hms2xhrs(sim_t0);
	gday_tf = ymd2gday(sim_datef) + hms2xhrs(sim_tf);	
	tunits_out = "hours since " + gt2string(gday_tb);
	
	// set spinup flag
	if (lspinup) ip_curr_yr = gt2year(spin_gday_t0);
	else ip_curr_yr = gt2year(gday_t0);

	// Create the model grid...
	mglons = createCoord(mglon0, mglonf, mgdlon, mgnlons);
	mglats = createCoord(mglat0, mglatf, mgdlat, mgnlats);
	mglevs.resize(1,1); mgnlevs = 1;

	// ... and time vector
	nsteps = (gday_tf - gday_t0)*24/dt + 1;
	nsteps_spin = (gday_t0 - spin_gday_t0)*24/dt;	// no +1 because this is 1 step less that t0
	mgtimes.resize(nsteps);
	for (int i=0; i<nsteps; ++i) mgtimes[i] = (gday_t0 + i*dt/24.0 - gday_tb)*24.0;
	
	// number of steps after which to show a dot so that 40 dots make up 100%
	dstep = nsteps/40;	
	if (dstep == 0) ++dstep;
	
	// get indices of cell containing SP-output point
	i_xlon = indexC(mglons, xlon);
	i_xlat = indexC(mglats, xlat);

}



/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	--> read_veg_params_file()
	
	READ VEG TYPE PARAMS FROM FILE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
inline int char2phenostage(char c){
	     if (c == 'X') return psX;
	else if (c == 'F') return psF;
	else if (c == 'M') return psM;
	else if (c == 'S') return psS;
	else if (c == 'Z') return psZ;
	else if (c == 'E') return psE;
	else return psX;
}

int read_veg_params_file(){
	ifstream fin;
	fin.open(params_ft_file.c_str());

	const string attrbegin = ">";

	string s, u, v, w;
	int n, m, l;
	float f;
	while (fin >> s && s != attrbegin);	// read until 1st > is reached
	
	fin >> s; 
	if (s != "nPFT") {cout << "number of PFTs not found!"; return 1;}
	fin >> npft;	// read number of pfts into global variable
	
	// resize all veg related vectors to npft values
	aLf.resize(npft); aSf.resize(npft); aRf.resize(npft);	// allocation fractions during flushing
	aL.resize(npft); aS.resize(npft); aR.resize(npft);		// allocation fractions during current phenology state
	LAImax.resize(npft); LAImin.resize(npft);				// min and max LAI for each pft
	phenoStages.resize(npft*12);	// phenology stages
	rFixC.resize(npft*12);			// monthly carbon fixation rates
	aFixC.resize(npft*12);			// monthly carbon fixation fractions 
	leafLs.resize(npft); Tdecomp.resize(npft);
	z1Month.resize(npft, -1);		// 1st leafless month, -1 if not found.
	Wc_sat_vec.resize(npft, 0);
	
	// start reading values into these vectors
	while (fin >> s && s != attrbegin){
		if (s == "") continue;	// skip empty lines
		if (s == "#") {getline(fin,s,'\n'); continue;}	// skip #followed lines (comments) 
		if		(s == "LM")	for (int i=0; i<npft; ++i) fin >> LAImax[i];
		else if (s == "Lm")	for (int i=0; i<npft; ++i) fin >> LAImin[i];
		else if (s == "aL")	for (int i=0; i<npft; ++i) fin >> aLf[i];
		else if (s == "aS")	for (int i=0; i<npft; ++i) fin >> aSf[i];
		else if (s == "LL") for (int i=0; i<npft; ++i) fin >> leafLs[i];
		else if (s == "T")  for (int i=0; i<npft; ++i) fin >> Tdecomp[i];
		else if (s == "ZM") for (int i=0; i<npft; ++i) fin >> z1Month[i]; 
		else if (s == "Wcs") for (int i=0; i<npft; ++i) fin >> Wc_sat_vec[i]; 

		else if (s == "rhobL") 		fin >> rhobL; 
		else if (s == "theta_sL") 	fin >> theta_sL; 
	}

	fin >> s; //cout << "s = " << s << '\n';
	char c; 
	if (s != "PHENO") {cout << "phenology not found!"; return 1;}
	for (int m=0; m<12; ++m){
		fin >> c;	// ignore the first char which is for month
		for (int i=0; i<npft; ++i){
			fin >> c;
			phenoStages[IX2(i,m, npft)] = char2phenostage(c); 
		}
	}
	
	while (fin >> s && s != attrbegin);	// loop till next >
	fin >> s; 
	if (s != "CARBON_FIXATION_RATE") {cout << "C Fixation rates not found!"; return 1;}
	for (int m=0; m<12; ++m){
		fin >> c;	// ignore the first char which is for month
		for (int i=0; i<npft; ++i){
			fin >> rFixC[IX2(i,m, npft)]; 
		}
	}

	// create aFixC vector (same as rFixC but negative values set to Zero)
	for (int m=0; m<12; ++m){
		for (int i=0; i< npft; ++i) {
			aFixC[IX2(i,m, npft)] = (rFixC[IX2(i,m, npft)] > 0)? rFixC[IX2(i,m, npft)] : 0;  
		}
	}
	
	// check if everything is correct
	log_fout << "--------- C fixation rates ---------------\n";
	for (int i=0; i<12; ++i){
		for (int j=0;j<npft;++j) { log_fout << rFixC[npft*i+j] << "\t";}
		log_fout << "\n";
	}
	log_fout << "\n";
	log_fout << "--------- 1st leafless month ---------------\n";
	for (int j=0;j<npft;++j) { log_fout << z1Month[j] << "\t";}
	log_fout << "\n";
	log_fout << "------------------------\n";

}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	--> init_modelvar(...)
	
	create a single gVar to be used in model and set metadata with model grid.
	if it is to be written to output, init the NcFile_handle and write
	metadata to Nc file.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
int init_modelvar(string var_name, string unit, int nl, gVar &v, vector <double> times_vec, ostream& lfout){
	// create levels vector
	vector <float> levs_vec(nl,0); 
	for (int i=0; i<nl; ++i) levs_vec[i] = i; 

	// create model variable 
	v = gVar(var_name, unit, tunits_out);
	v.setCoords(times_vec, levs_vec, mglats, mglons);
	v.values.resize(mgnlons*mgnlats*nl);
	v.printGrid(lfout);

	// set bool values for variables to output
	v.lwrite = (writeVar_flags[var_name] == 1) ? true:false;
	v.lwriteSP = (writeVarSP_flags[var_name] == 1) ? true:false;
	v.lwriteSP &= spout_on;
	
	// create output files if desired	
	if (v.lwrite){	
		lfout << "~~~~~~~~ Will write variable " << var_name << " to nc file.\n";
		v.ofname = "../output/"+var_name+"."+sim_date0+"-"+sim_datef+".nc";
		v.ofile_handle->open(v.ofname, "w", glimits_india);
		v.ofile_handle->writeCoords(v);
		v.ofile_handle->writeTimeValues(v);
		v.outNcVar = v.ofile_handle->createVar(v);
		//setZero(v.values);
	}
}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	--> init_ip_var(...)

	This function is used to initialize all variable data (for a single 
	variable) for the first time as well as to reload files when the current 
	file ends.
	
	When initialising for 1st time, it will create all the filenames
	which are to be read sequentially, then read variable metadata, and 
	finally create the model variables (with lterp indices) corresponding 
	to the input.  
	
	When entering after 1st init, it will open next set of files and read
	metadata.
	
	input data is read into v, interplated into vout, and messages are sent
	to lfout

 	One time init code is included here to avoid retyping fnctions 
 	for each variable separately! Thats quite lazy of me but its ok! :P
	This code:
		 creates filenames from ip_data 
		 creates corresponding modelvars 
		 does not execute during reload of new input files.

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
int init_ip_var(string var_name, gVar &v, gVar &vout, ostream &lfout){
	
	//cout << "INIT for variable = " << var_name << '\n';
	vector <string> &var_files = ip_data_map.find(var_name)->second.fnames;	// a reference to files vector

	// ~~~~~~~~~~~ one time init ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	if (!l_ip_init_done){
		logdebug << "First time init. Creating input filenames...\n\n";
		// set correct filenames for input data
		int m = ip_data_map.find(var_name)->second.start_yr;
		int n = ip_data_map.find(var_name)->second.nyrs_file;
		string u = ip_data_map.find(var_name)->second.fname_prefix;
		string var_dir = data_dirs[var_name];
		string parent_dir = data_dirs["forcing_data_dir"];
		logdebug << "parent dir = " << data_dirs["forcing_data_dir"] << "\n\n";
		
		if (n > 1) {
			string vfile = parent_dir + "/" + var_dir + "/" + u 
					 + "." + int2str(m) + "-" + int2str(m+n-1) + ".nc";
			lfout << "filename created: " << vfile << '\n';
			var_files.push_back(vfile);
			lfout << '\n';
		}
		else{
			for (int i=0; i<(ip_end_yr-m+1); ++i){
				string vfile = parent_dir + "/" + var_dir + "/" + u 
						 + "." + int2str(m+i) + ".nc";
				lfout << "filename created: " << vfile << '\n';
			var_files.push_back(vfile);
			}
			lfout << '\n';
		}
		// initialize model var corresponding to input var
		string var_unit = ip_data_map.find(var_name)->second.unit;
		init_modelvar(var_name, var_unit, 1, vout, mgtimes, lfout);
	}
	// ~~~~~ one time init done! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	

	// for input vars, open files, read metadata, create lterp indices
	if (var_files.size() > 1 || !l_ip_init_done){	// open only if variable has many files or if its first time
		int filenum; 
		if (var_files.size() == 1) filenum = 0;
		else filenum = ip_curr_yr - ip_data_map.find(var_name)->second.start_yr;
		
		lfout << "Loading (" << var_name << ") input file for year " << ip_curr_yr << ", filenum = " << filenum << "...\n";	
		lfout << "File = " << var_files[filenum] << '\n'; cout.flush();
		v.ifname = var_files[filenum];
		int i = v.ifile_handle->open(v.ifname, "r", glimits_india);
		if (i != 0) cout << "NCFILE NOT VALID!!\n";

		v.ifile_handle->readCoords(v, lfout);
		v.ifile_handle->readVarAtts(v);
		v.printGrid(lfout);
		lfout.flush();

		// calculate regriding indices for {mglons, mglats} from v's grid
		ip_data_map.find(var_name)->second.lterp_indices = bilIndices(v.lons, v.lats, mglons, mglats);


		// ip_data is completely filled. We must finally add an entry to the map with the 
		// this same ip-data but mapped to the varname as in file, since the var_name and the name 
		// in the file may be different
		ip_data a(ip_data_map.find(var_name)->second);
		ip_data_map.insert( pair <string, ip_data> (v.varname,a) );
		logdebug << "updated ipdata: " << v.varname << ": "; a.print_vardata(log_fout);

	}

}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	--> read_static_var(...)
	
	read the metadata and values from file for static variable v 
	and interpolate them to model grid.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
int read_static_var(gVar &v, ostream &lfout){
	int i = v.ifile_handle->open(v.ifname, "r", glimits_india);
	if (i != 0) cout << "NCFILE NOT VALID!!\n";

	v.ifile_handle->readCoords(v, lfout);
	v.ifile_handle->readVarAtts(v);
	v.ifile_handle->readVar(v,0);

	v = lterp(v, mglons, mglats);
	v.printGrid(lfout);
	lfout.flush();
}



/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	--> init_all_ip_vars()
	--> init_all_model_vars()
	--> read_all_static_vars()
	
	Call the single variable functions one by one on each variable.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

int init_all_ip_vars(){

	PREPROC_INIT_IP_VAR(pr);
//	PREPROC_INIT_IP_VAR(ps);
	PREPROC_INIT_IP_VAR(rh);
	PREPROC_INIT_IP_VAR(ts);
//	PREPROC_INIT_IP_VAR(ndr);
	PREPROC_INIT_IP_VAR(wsp);
	PREPROC_INIT_IP_VAR(npp);

	l_ip_init_done = true;
};

int init_all_modelvars(){

	log_fout << "Creating model variables... \n\n"; 
	// create all model gVars 
	PREPROC_INIT_MODELVAR(canbio, "gC/m2", npft);
	PREPROC_INIT_MODELVAR(canbio_max, "gC/m2", npft);
	PREPROC_INIT_MODELVAR(dxl, "cm", 1);
	PREPROC_INIT_MODELVAR(lmois, "kg/m2", 1);
	PREPROC_INIT_MODELVAR(cmois, "kg/m2", 1);
	PREPROC_INIT_MODELVAR(fire, "-", 1);
	PREPROC_INIT_MODELVAR(ps, "Pa", 1);
	PREPROC_INIT_MODELVAR(ndr, "W/m2", 1);
	PREPROC_INIT_MODELVAR(evap, "mm/day", 1);
	
	// create daily fire gVar
	int dfnt = mgtimes.size()*dt/24;
	vector <double> t_temp(dfnt);
	for (int i=0; i<dfnt; ++i) t_temp[i] = mgtimes[i*24/dt];

	if (varUse_map["dfire"]) init_modelvar("dfire", "-", 1, dfire, t_temp, log_fout);

};


int read_all_static_vars(){
	log_fout << "Reading static variables... \n\n"; 
	PREPROC_READ_STATIC_VAR(msk);
	PREPROC_READ_STATIC_VAR(vegtype);
	PREPROC_READ_STATIC_VAR(elev);
	PREPROC_READ_STATIC_VAR(albedo);
} 


/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	--> init_infisim()
	
	Call all init commands to do full sim init and show progress.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
int init_infisim(){
	log_fout.open("../output/log.txt");	// open log stream
	log_fout << " ******************* THIS IS LOG FILE ****************************\n\n";
	log_fout << "~~~~~~~~~~~~~~~~~~~~~~~~ " << ip_curr_yr << " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";

	cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
	cout << "~                  I N F I S I M                               ~\n";
	cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
	cout << "\n> Reading config parameters... "; cout.flush();
	read_sim_config_file();
	cout << "DONE.\n> Reading input filenames... "; cout.flush();
	read_ip_params_file();
	cout << "DONE.\n> Reading forest type params... "; cout.flush();
	read_veg_params_file();
	cout << "DONE.\n> Initialising input vars/files... "; cout.flush();
	init_all_ip_vars();	
	cout << "DONE.\n> Initialising model vars/files... "; cout.flush();
	init_all_modelvars();
	cout << "DONE.\n> Reading static input vars... "; cout.flush();
	read_all_static_vars();
	cout << "DONE.\n";

	// check for consistency in vegtype levels and PFTs
	if (npft != vegtype.nlevs) {
		cout << "** ERROR ** : number of PFTs dont match levels in vegtype!\n\n";
		return 1;
	}

	// convert elevation to surface pressure
	for (int ilat=0; ilat<mgnlats; ++ilat){
		for (int ilon=0; ilon<mgnlons; ++ilon){
	
			ps( ilon, ilat, 0) = 101325 - 1200*elev(ilon,ilat,0)/100;

		}
	}


	if (spout_on){
		cout << "> Opening file to write point values: " << pointOutFile << '\n';
		sp_fout.open(pointOutFile.c_str());
	
		sp_fout << "lat:\t " << xlat << "\t lon:\t" << xlon << '\n';
		sp_fout << "vegtype fractions:\n X\t AGR\t NLE\t BLE\t MD\t DD\t GR\t SC\n";
		for (int i=0; i<npft; ++i){
			sp_fout << vegtype.getCellValue(xlon,xlat, i) << '\t'; 
		}
		sp_fout << '\n';

		sp_fout << "datetime\t"
				<< "pheno_phase\t canbio[u]\t litbio[u]\t pheno_phase\t canbio[u]\t litbio[u]\t"
				<< "ndr\t LAI\t Rn\t pre\t Ep\t Ea\t dzL\t S\t"
				<< "fire\t dfire \n";
	}

	
	return 0;
}






