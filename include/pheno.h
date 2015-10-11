#ifndef PHENO_H
#define PHENO_H

#include <vector>

vector <float> calc_alloc_pft(int m, float npp_obs, vector <float> &pft_fracs);
int set_leaf_alloc(int m, vector <float> &nppAllocs);
vector <float> calc_litterfall_rate(int m, vector <float> &canbio_now, float delT);
int calc_pheno(float gtime, float delT);







#endif

