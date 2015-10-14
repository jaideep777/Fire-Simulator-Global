#ifndef PHENO_H
#define PHENO_H

#include <vector>

vector <float> calc_alloc_pft(int m, float lat, float npp_obs, vector <float> &pft_fracs);
int set_leaf_alloc(int m, float lat, vector <float> &nppAllocs);
vector <float> calc_litterfall_rate(int m, float lat, vector <float> &canbio_now, float delT);
int calc_pheno(float gtime, float delT);







#endif

