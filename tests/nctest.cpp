#include <iostream>
#include <math.h>
#include <netcdfcpp.h>
#include <fstream>
#include <vector>
using namespace std;

#include "../include/defs.h"
#include "../include/globals.h"
#include "../include/arrayutils.h"
#include "../include/time.h"
#include "../include/ncio.h"
#include "../include/grid.h"


int main(){
//	float lon0 = 66.5, lat0 =  6.5, lonf = 100.5, latf = 38.5, dx = 0.5, dy = 0.5, dt = 1, dlev = 1;
	float lon0 = 75.0, lat0 = 15.0, lonf =  80.0, latf = 20.0, dx = 0.1, dy = 0.1, dt = 1, dlev = 1;
	int nlons=0, nlats=0, nlevs=1, ntimes=1;
	ifstream fin;
	string prop = "../output/pres_out_2000_masked.nc";

	// create model grid (lat/lon) some variable may have levels, some not
	vector <float> lons = createCoord(lon0, lonf, dx, nlons);
	vector <float> lats = createCoord(lat0, latf, dy, nlats);
//	vector <float> levs = createCoord( 0.f, 0.f, 1.0f, nlevs);
//	vector <float> times = createCoord(1.f, 3.f, 1.0f, ntimes);
	
	// set NETCDF error behavior to non-fatal
	NcError err(NcError::verbose_nonfatal);

	// open file
	NcFile_handle mfile(maskadd, "r");
	mfile.setMapLimits(lat0,latf,lon0,lonf);
	
	gVar vmask;

	mfile.readCoords(vmask, true);
	mfile.readVarAtts(vmask, vmask.ncoords);
	vmask.printGrid();
	mfile.readVar(vmask, vmask.ncoords, 0);
//	vmask = lterp(vmask, lats, lons);
//	vmask.printGrid();

	gVar msk = lterp(vmask, lats, lons);

	NcFile_handle infile(prip, "r");
	NcFile_handle ofile(prop, "w");
	infile.setMapLimits(6.5,38.5,66.5,100.5);

	gVar pres_in; //("pres", "num fires/day", "days since 2004-12-31");

	infile.readCoords(pres_in, true);
	infile.readVarAtts(pres_in, pres_in.ncoords);
	pres_in.ntimes = ntimes; pres_in.times = times;
	pres_in.printGrid();

	cout << " a ****************************\n";
	
	gVar pres_out(pres_in);
	pres_out.setCoords(times, levs, lats, lons);
	pres_out.printGrid();

	ofile.writeCoords(pres_out, true);
	NcVar * pvar = ofile.createVar(pres_out);

	for (int t=0; t<ntimes; ++t){
		infile.readVar(pres_in, pres_in.ncoords, t);
		gVar y = lterp(pres_in, lats, lons);
		gVar x = mask(y, msk);
		ofile.writeVar(x, pvar, t);	
	}	
	
	ofile.writeTimeValues(pres_out);
	ofile.closeFile();	
	cout << "\nSuccessfully wrote file.\n";
	
	cout << " a ****************************\n";
	
}

