#include <iostream>
using namespace std;

#include "../include/arrayutils.h"
//#include "../include/gvar.h"
//#include "../include/ncio.h"

int main(){

	const int N = 59;
	vector <float> a(N);
	for (int i=0; i<N; ++i){
		a[i] = 7.0 + i*0.5;
	}

	printArray(a,N);
//	reverseArray(a,8);
//	printArray(a,8);
	float y; int m, n;

	y = 6.25;
	m = indexLow(y,  a);
	n = indexHi(y,  a);
	cout << "Bounds: (" << a[m] << " < " << y << " < " << a[n] << ")\n";
	
	y = 7;
	m = indexLow(y,  a);
	n = indexHi(y,  a);
	cout << "Bounds: (" << a[m] << " < " << y << " < " << a[n] << ")\n";

	y = 20.0;
	m = indexLow(y,  a);
	n = indexHi(y,  a);
	cout << "Bounds: (" << a[m] << " < " << y << " < " << a[n] << ")\n";

	y = 20.15;
	m = indexLow(y,  a);
	n = indexHi(y,  a);
	cout << "Bounds: (" << a[m] << " < " << y << " < " << a[n] << ")\n";

	y = 20.75;
	m = indexLow(y,  a);
	n = indexHi(y,  a);
	cout << "Bounds: (" << a[m] << " < " << y << " < " << a[n] << ")\n";

	y = 36;
	m = indexLow(y,  a);
	n = indexHi(y,  a);
	cout << "Bounds: (" << a[m] << " < " << y << " < " << a[n] << ")\n";

	y = 36.25;
	m = indexLow(y,  a);
	n = indexHi(y,  a);
	cout << "Bounds: (" << a[m] << " < " << y << " < " << a[n] << ")\n";

//	int i;
//	i = (false)? 1:0;
//	cout << i << "\n\n";
//	
//	vector <float> b = copyArray(a, m, l);
//	printArray(b);
//	vector <float> lons(5);
//	for (int i=0; i<lons.size(); ++i) lons[i] = 75 + 5*i;

//	vector <float> vals(lons.size()*b.size());
//	for (int i=0; i<vals.size(); ++i) vals[i] = i;
//	
//	vector <float> levs(1), times(1);

//	gVar aa("testvar", "testunits", "days since 2007-12-01 12:0:0");
//	reverseArray(b);
//	aa.setCoords(times, levs, b, lons);
//	aa.values = vals;
//	
//	aa.printGrid();
//	aa.printValues();
//	
////	gVar bb(aa);
////	bb.printGrid();
////	bb.printValues();
////	bb = aa + aa;
////	bb.printGrid();
////	bb.printValues();

}

