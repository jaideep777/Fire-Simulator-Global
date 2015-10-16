#ifndef VARS_H
#define VARS_H

#include <string>
#include <gsm.h>


//  Static Variables
// mask (msk) 
extern string mask_file;
extern gVar msk;	

// forest type (vegtype) 
extern string vegtype_file;
extern gVar vegtype;	

extern gVar elev;
extern gVar albedo;

// lmois initial condition
extern string lmois_ic_file;		// this will be read into variable lmois during init


//	Time varying Input variables

//extern gVar psi, ps;	// Surface Pressure (ps)
extern gVar pri, pr;	// Precipitation (pr) 
extern gVar rhi, rh;	// Relative Humidity (rh) 
extern gVar tsi, ts;	// Surface Temperature (ts) 
//extern gVar ndri, ndr;	// Net Downward Radiation (ndr) 
extern gVar wspi, wsp;	// Wind Speed (wsp) 
extern gVar nppi, npp;	// NPP (npp) 

//extern vector <int> pr_indices;


// Time varying model-only variables

extern gVar canbio;		// canopy biomass
extern gVar canbio_max;	// max canopy biomass for LAI calculation, set during canbio prerun
extern gVar lmois;		// litter moisture (kg/m2 = mm)
extern gVar cmois;		// canopy moisture content (kg/m2 = mm) 
extern gVar dxl;		// litter layer thickness
extern gVar fire;		// fire!
extern gVar dfire;		// daily fire indices
extern gVar ndr; 		// net downward radiation
extern gVar ps;			// surface pressure
extern gVar evap;			// potential evaporation rate

extern gVar dts, dwsp, drh, dlmois, ddxl;	// daily equivalents of weather variables, to store conditions at max fire

#endif

