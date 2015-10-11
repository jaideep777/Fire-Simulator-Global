#include <iostream>
#include "../include/time.h"
using namespace std;

int main(){

	string stime = "15:35:48";
	string sdate = "1988-10-18";
//	float dayfraction = 15.f/24.f + 35.f/24.f/60.f + 48.f/24.f/3600.f;
//	cout << dayfraction << '\n';
//	cout << gtime(dayfraction) << '\n';
//	cout << gdayfrac(s) << '\n';
	
	double tbase = gdays(sdate) + gdayfrac(stime); // + 5.5/24
	
	cout << "gday float : " << tbase << '\n';
	cout << "day fraction was: " << gdayfrac(stime) << '\n';
	cout << "ripped day frac : " << tbase - int(tbase) << '\n';
	cout << "date is : " << gdate(int(tbase)) << '\n';
	cout << "time is : " << gtime(tbase - int(tbase)) << '\n';
	cout << "time was: " << gtime(gdayfrac(stime)) << '\n';
	return 0;

}
