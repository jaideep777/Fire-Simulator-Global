#include <iostream>
using namespace std;

#include "/home/jaideep/libgsm/include/gsm.h"

/*
int main(){

	int nlons = 5, nlats = 7, ntimes = 1, nlevs = 1;	
	
	vector <float> lons(nlons);
	for (int i=0; i<lons.size(); ++i) lons[i] = 75 + 5*i;

	vector <float> lats(nlats);
	for (int i=0; i<lats.size(); ++i) lats[i] = 5 + 5*i;

	vector <float> data(nlons*nlats);
	for (int i=0; i<data.size(); ++i) data[i] = i;

	printVar(lons, lats, data);

	while (1){
		float x, y;
		cout << "Enter coordinates: ";
		cin >> x >> y;
	
		float f1 = bilinear(x,y,0, lons, lats, data);
		cout << "lterp value = " << f1 << "\n";
		float f2 = cellVal(x,y,0, lons, lats, data);
		cout << "Cell value  = " << f2 << "\n\n";
		
	}

}
*/

int main(){

	debug_on = false;

	int nlons = 5, nlats = 7, ntimes = 1, nlevs = 1;	
	
	vector <float> lons(nlons);
	for (int i=0; i<lons.size(); ++i) lons[i] = 75 + 5*i;

	vector <float> lats(nlats);
	for (int i=0; i<lats.size(); ++i) lats[i] = 5 + 5*i;

	vector <float> data(nlons*nlats);
	for (int i=0; i<data.size(); ++i) data[i] = i;

	vector <float> mdata(nlons*nlats, 1);
	for (int i=0; i<nlons*nlats; i=i+2) mdata[i] = 0;

	vector <float> mdata2(nlons*nlats, 0);
	for (int i=0; i<nlons*nlats; i=i+2) mdata2[i] = 1;

	vector <double> times(1,1);
	vector <float> levs(1,1);

	int xnlons = 15, xnlats = 17;	
	
	vector <float> xlons(xnlons);
	for (int i=0; i<xlons.size(); ++i) xlons[i] = 70 + 2.5*i;

	vector <float> xlats(xnlats);
	for (int i=0; i<xlats.size(); ++i) xlats[i] = 0 + 2.5*i;

	cout << "1 *************************\n";	
	gVar v("var (v)", "K", "hours since 1988-10-18 8:57:00");
	v.setCoords(times, levs, lats, lons);
	v.values = data;
	v.values[13] = v.missing_value;
	cout << "2 *************************\n";	
	v.printGrid();
	v.printValues();

	// lterp test
	gVar t; t.setCoords(times, levs, xlats, xlons);
	gVar r = lterp(v,t.lats, t.lons);
	r.varname = "lterpvar (r)";
	r.printGrid();
	r.printValues();

	// alt lterp test
	gVar r1 = lterp(v, t.lats, t.lons);
	r1.varname = "lterpvar (r1)";
	r1.printGrid();
	r1.printValues();

	// mask
	gVar m; 
	m.shallowCopy(v);
	m.values = mdata;
	m.varname = "mask (m)"; m.varunits = "-";
	m.printGrid();
	m.printValues();
	
	// mask 2
	gVar m1; 
	m1.shallowCopy(v);
	m1.values = mdata2;
	m.varname = "mask2 (m)"; m.varunits = "-";
	
	// mask test
	gVar w = mask(v, m);
	w.varname = "masked_var (w)"; w.varunits = "varunits";
	w.printGrid();
	w.printValues();

	// sum test
	gVar x = v+m;
	x.varname = "sum (x)";
	x.printGrid();
	x.printValues();
	
	// masked sum test
	gVar z = v+w;
	z.varname = "masked sum (z)";
	z.printGrid();
	z.printValues();

	// double fillvalues sum test
	gVar q = mask(v,m) + mask(v,m1);
	q.varname = "double masked sum (q)";
	q.printGrid();
	q.printValues();

}


