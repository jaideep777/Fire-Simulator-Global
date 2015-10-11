#include <iostream>
#include <math.h>
#include <netcdfcpp.h>
#include <fstream>
#include <vector>
using namespace std;

#include "/home/jaideep/libgsm/include/gsm.h"

int main(){

	debug_on = false;
	
	float lon0 = 66.5, lat0 = 6.5, lonf = 100.5, latf = 38.5, dx = 0.5, dy = 0.5, dt = 1, dlev = 1;
	//	float lon0 = -90, lat0 = 0, lonf = 90, latf = 90, dx = 0.1, dy = 0.1, dt = 1, dlev = 1;
	int nlons=-1, nlats=-1, nlevs=-1, ntimes=3;
//	bool latOrderSN = true;
//	float *lons, *lats, *levs, *times;
//	static const int NC_ERR = 2;
	ifstream fin;
//	vector <float> timevals;
//	string    fname = "../data/sample2d.nc";
//	string    wfname = "../data/testwrite1.nc";
//	string ascfname = "../data/2005_sorted.csv";
	string prip = "/media/WorkData/Fire/ncep20cen/fire_input/pressure_sfc/pres.sfc.2000.nc";
	string prop = "../output/pres_out_2000_masked.nc";
	string maskadd = "../input/surta_india_0.2.nc";

	vector <float> lons = createCoord(lon0, lonf, dx, nlons);
	vector <float> lats = createCoord(lat0, latf, dy, nlats);
	vector <float> levs = createCoord( 0.f, 0.f, 1.0f, nlevs);
	vector <double> times(3,0); times[1]=1.0; times[2]=2.0;
	
	// set error behavior to non-fatal
	NcError err(NcError::silent_nonfatal);

	// open file
	NcFile_handle infile; infile.open(prip, "r", glimits_india);
	NcFile_handle ofile; ofile.open(prop, "w", glimits_india);
	NcFile_handle mfile; mfile.open(maskadd, "r", glimits_india);
//	infile.setMapLimits(6.5,38.5,66.5,100.5);
//	mfile.setMapLimits(6.5,38.5,66.5,100.5);
	
	gVar pres_in; //("pres", "num fires/day", "days since 2004-12-31");
	gVar vmask;

	infile.readCoords(pres_in);
	infile.readVarAtts(pres_in);
	pres_in.ntimes = ntimes; pres_in.times = times;
	pres_in.printGrid();

	mfile.readCoords(vmask);
	mfile.readVarAtts(vmask);
	vmask.printGrid();
	
	mfile.readVar(vmask, 0);
	vmask = lterp(vmask, lons, lats);
	vmask.printGrid();
	
	cout << " a ****************************\n";
	
	gVar pres_out(pres_in);
	pres_out.setCoords(times, levs, lats, lons);
	pres_out.printGrid();

	ofile.writeCoords(pres_out);
	NcVar * pvar = ofile.createVar(pres_out);

	for (int t=0; t<ntimes; ++t){
		infile.readVar(pres_in, t);
		cout << "t = " << t << "\n\n";
		pres_in.printValues();
		gVar y = lterp(pres_in, lons, lats);
		cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
		y.printValues();
		gVar pres_out = mask(y, vmask);
		pres_out.printValues();
		ofile.writeVar(pres_out, pvar, t);	
		cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
	}	
	
	ofile.writeTimeValues(pres_out);
	ofile.close();	
	cout << "\nSuccessfully wrote file.\n";
}

