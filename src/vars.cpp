#include "../include/vars.h"


//  Static Variables
// mask (msk) 
string mask_file;
gVar msk;	

// forest type (vegtype) 
string vegtype_file;
gVar vegtype;	

// lmois initial condition
string lmois_ic_file;		// this will be read into variable lmois during init

gVar elev;
gVar albedo;


//	Time varying Input variables

//gVar psi, ps;	// Surface Pressure (ps)
gVar pri, pr;	// Precipitation (pr) 
gVar rhi, rh;	// Relative Humidity (rh) 
gVar tsi, ts;	// Surface Temperature (ts) 
//gVar ndri, ndr;	// Net Downward Radiation (ndr) 
gVar wspi, wsp;	// Wind Speed (wsp) 
gVar nppi, npp;	// NPP (npp) 

//vector <int> pr_indices;

// Time varying model-only variables

gVar canbio;		// canopy biomass
gVar canbio_max;	// max canopy biomass (needed for LAI calculation)
gVar lmois;			// litter moisture content (S)
gVar cmois;			// canopy moisture content (kg) 
gVar dxl;			// litter layer thickness
gVar fire;			// fire!
gVar dfire;			// daily fire indices
gVar ndr; 			// net downward radiation
gVar ps;			// surface pressure
gVar evap;			// potential evaporation rate

