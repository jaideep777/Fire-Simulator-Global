#include <iostream>
#include <vector>
#include <algorithm>
using namespace std;

const float glimits_globe[4] = {0, 360, -90, 90};
const float glimits_india[4] = {66.5, 100.5, 6.5, 38.5};
float glimits_custom[4] = {0,0,0,0};

bool ascComp(float a, float b){ return (a<b); }
bool dscComp(float a, float b){ return (a>b); }

// lower_bound = 1st element !< (>=) val
// upper_bound = last element !> (<=) val

// returns lower bound in array for val. 
// if val is out of range, returns appropriate edge. does not return missing value
int ncIndexLo(vector <float> &v, float val){
	bool asc = (v[1]>v[0])? true:false;
	if (asc){
		if (val >= v[v.size()-1]) return (v.size()-1);
		else if (val <= v[0]) return 0;
		else return (upper_bound(v.begin(), v.end(), val, ascComp) - v.begin() - 1);
	}
	else{
		if (val <= v[v.size()-1]) return (v.size()-1);
		else if (val >= v[0]) return 0;
		else return (lower_bound(v.begin(), v.end(), val, dscComp) - v.begin());
	}
}

// returns upper bound in array for val. 
// if val is out of range, returns appropriate edge. does not return missing value
int ncIndexHi(vector <float> &v, float val){
	bool asc = (v[1]>v[0])? true:false;
	if (asc){
		if (val >= v[v.size()-1]) return (v.size()-1);
		else if (val <= v[0]) return 0;
		else return (lower_bound(v.begin(), v.end(), val, ascComp) - v.begin());
	}
	else{
		if (val <= v[v.size()-1]) return (v.size()-1);
		else if (val >= v[0]) return 0;
		else return (upper_bound(v.begin(), v.end(), val, dscComp) - v.begin() -1);
	}
}

// find index (m) of grid box such that P lies between m and m+1
// if val is out of range, returns missing value (-999)
//    1------2---P--3  P
//             ^G.C    ^outlier
int lindexSW(vector <float> &v, float val){
	bool asc = (v[1]>v[0])? true:false;
	if (asc){
		if (val > v[v.size()-1] || val < v[0]) return -999;	// if val exceeds edges, return -999
		else if (val == v[v.size()-1]) return v.size() -2;	// if val is on right edge, return 1 less
		else return (upper_bound(v.begin(), v.end(), val, ascComp) - v.begin() - 1);
	}
	else{
		// this case must not be used.. nonetheless, it may not cause trouble in interpolation
		cout << "**WARNING**: lindexSW() invoked on descending vector!\n";
		if (val < v[v.size()-1] || val > v[0]) return -999;	// if val exceeds edges, return -999
		else if (val == v[v.size()-1]) return v.size() -2;	// if val is on right edge, return 1 less
		else return (upper_bound(v.begin(), v.end(), val, dscComp) - v.begin() -1);
	}
}

// find index (m) of grid box such that P lies closer to m than m+1 or m-1
// if val is out of range, returns missing value (-999)
// |---1---|---2-P-|---3-P-|     P
//               ^G.C    ^Sp.C   ^outlier
int indexC(vector <float> &v, float val){
	bool asc = (v[1]>v[0])? true:false;
	if (asc){
		if (val >= v[v.size()-1]){
			if (val <= (v[v.size()-1] + (v[v.size()-1]-v[v.size()-2])/2)) return (v.size()-1);
			else return -999;
		}
		if (val <= v[0]){
			if (val >= (v[0] - (v[1]-v[0])/2)) return 0;
			else return -999;
		} 
		else {
			int m = (upper_bound(v.begin(), v.end(), val, ascComp) - v.begin() - 1); // lower bound
			return ((val - v[m]) < (v[m+1]-val))? m:(m+1);
		}
	}
	else{
		// this case must not be used.. 
		cout << "**ERROR**: indexC() invoked on descending vector!\n";
		return -999;
	}
}


float getA(vector<float> &v, int m){return (m < 0)? 9.9e20:v[m];}

int main(){

	const int N = 12;
	vector <float> a(N);
	for (int i=0; i<N; ++i){
//		a[i] = N-i;
		a[i] = i+1;
		cout << a[i] << " ";
	}
	cout << '\n';
	
	float x, y;
	int m, n;

	while (cin >> y){
		cout << "**** Vector code\n";
		m = ncIndexLo(a,y); //int(lower_bound(a.begin(), a.end(), y, ascComp) - a.begin());
		n = ncIndexHi(a,y); //int(upper_bound(a.begin(), a.end(), y, ascComp) - a.begin());
		cout << "indices:         (" << m << ", " << n << ")\n";
		cout << "Bounds:          (" << a[m] << " < " << y << " < " << a[n] << ")\n";
		cout << "**** GridBoxC\n";
		m = lindexSW(a,y); //int(lower_bound(a.begin(), a.end(), y, ascComp) - a.begin());
//		cout << "index:     " << m << '\n';
		cout << "lindexSW : " << getA(a,m) << "\n";
		m = indexC(a,y);
//		cout << "index:     " << m << '\n';
		cout << "lindexC  : " << getA(a,m) << "\n";

	}
}

