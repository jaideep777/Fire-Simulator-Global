#ifndef _PRERUN_H
#define _PRERUN_H

#include "../include/globals.h"
#include "../include/vars.h"
#include "../include/init.h"
#include "../include/io.h"
#include "../include/pheno.h"

char ps2char(int c);
int printPreRunHeader(string s, int ns, int ds);
int prerun_canbio_ic();
int prerun_lmois_ic();

int main_run();

#endif
