#include <iostream>
#include <math.h>
#include <netcdfcpp.h>
#include <fstream>
#include <vector>
using namespace std;

#include "/home/jaideep/libgsm/include/gsm.h"

const float pi = 3.14159265;

int main(){

//	debug_on = false;
	
	float lon0 = 66.5, lat0 = 6.5, lonf = 100.5, latf = 38.5, dx = 0.5, dy = 0.5, dt = 1, dlev = 1;
	//	float lon0 = -90, lat0 = 0, lonf = 90, latf = 90, dx = 0.1, dy = 0.1, dt = 1, dlev = 1;
//	int nlons=-1, nlats=-1, nlevs=-1, ntimes=3;
//	bool latOrderSN = true;
//	float *lons, *lats, *levs, *times;
//	static const int NC_ERR = 2;
	ifstream fin;
//	vector <float> timevals;
//	string    fname = "../data/sample2d.nc";
//	string    wfname = "../data/testwrite1.nc";
//	string ascfname = "../data/2005_sorted.csv";
	string ndrip = "/media/WorkData/Fire/ncep20cen/fire_input/rad_ndr_sfc/ndr.sfc.2004.nc";
	string prip  = "/media/WorkData/Fire/ncep20cen/precip/prate.2004.nc";
	string ndrop = "ndr_out_2004.nc";
	string maskadd = "../input/surta_india_0.2.nc";

	// set error behavior to non-fatal
	NcError err(NcError::silent_nonfatal);

	// open file
	NcFile_handle infile;  infile.open(ndrip, "r", glimits_globe);
	NcFile_handle infilep; infilep.open(prip, "r", glimits_globe);
	NcFile_handle ofile;   ofile.open(ndrop, "w", glimits_globe);
	
	NcFile_handle alfile; 
	alfile.open("/media/WorkData/Fire/ncep20cen/albedo/albedo.avg.2004.nc", "r", glimits_globe);
	
	gVar ndr_in; 
	gVar pr_in;
	gVar albedo;

	infile.readCoords(ndr_in);
	infile.readVarAtts(ndr_in);

	infilep.readCoords(pr_in);
	infilep.readVarAtts(pr_in);

	alfile.readCoords(albedo);
	alfile.readVarAtts(albedo);
	alfile.readVar(albedo, 0);
	albedo.printValues();
	
	cout << " a ****************************\n";
	
	gVar ndr_out(ndr_in);
	ndr_out.printGrid();
	ndr_out.values.resize(ndr_out.nlats*ndr_out.nlons);

	ofile.writeCoords(ndr_out);
	ofile.writeTimeValues(ndr_out);
	NcVar * pvar = ofile.createVar(ndr_out);

	for (int t=0; t<ndr_in.ntimes; ++t){
		double gt = ndr_in.ix2gt(t);
		int m = gt2month(gt);

		infile.readVar(ndr_in, t);
		infilep.readVar(pr_in, t);

		float N = gt2daynum(gt);
		float hr = (gt - int(gt))*24;
		cout << "t = " << t << ", N = " << N << "\n";

		float y = 2*pi/365*(N-1 + (hr-12)/24);

		float eqtime = 229.18*(0.000075 + 0.001868*cos(y) - 0.032077*sin(y)
					 - 0.014615*cos(2*y) - 0.040849*sin(2*y));
		
		float decl = 0.006918 - 0.399912*cos(y) + 0.070257*sin(y) - 0.006758*cos(2*y)
				   + 0.000907*sin(2*y) - 0.002697*cos(3*y) + 0.00148*sin(3*y);

		for (int ilat=0; ilat< ndr_in.nlats; ++ilat){
			for (int ilon=0; ilon< ndr_in.nlons; ++ilon){

				float lat = ndr_in.lats[ilat]*pi/180;
				float lon = ndr_in.lons[ilon];
				if (lon <0) lon -= 360;

				float time_offset = -(eqtime - 4*lon); // + 60*5.5;	// 5.5 is timezone

				float tst = hr*60 + time_offset;
				
				float ha = (tst/4 - 180)*pi/180;
				
				float sinY = sin(lat)*sin(decl)+ cos(lat)*cos(decl)*cos(ha);

				if (ilon == 13 && ilat == 1) cout << "ha = " << ha << '\n'; 
				
				//float sinY = u - v*cos(h);
				float SW_d;
				if (sinY > 0) {
					SW_d = 1370*(0.6+0.2*sinY)*sinY;
//					float Q0 = 3600*1360*(1+2*0.01675*cos(2*pi*N/365));
//					SW_d = sinY*(0.25+0.5*8)*Q0;
				}
				else SW_d = 0;
				
				float cld = pr_in(ilon, ilat,0)/3e-4;
				if (cld < 0) cld = 0; if (cld > 1) cld = 1;
				if (t == 2003) cout << "pr = " << pr_in(ilon, ilat, 0) << " cld = " << cld << '\n'; 

				ndr_out(ilon, ilat, 0) = SW_d*(1-albedo(ilon,ilat,0)/100)*(1-cld*0.6) - 100; // 100 is lw_in - lw_out
			}
		}

		ofile.writeVar(ndr_out, pvar, t);
	}	
	
	ofile.close();	
	cout << "\nSuccessfully wrote file.\n";
}


