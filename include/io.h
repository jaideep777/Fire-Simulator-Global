#ifndef _IO_H
#define _IO_H

#include <iostream>
#include <math.h>
#include <netcdfcpp.h>
#include <fstream>
#include <vector>
using namespace std;

#include <gsm.h>

#define PREPROC_READ_IP_VAR(x, ret) if (varUse_map[#x]) ret = read_ip_var(x##i, x, gtime, mode)
#define PREPROC_WRITE_OUTVAR(x) if (varUse_map[#x]) write_outvar(x, istep)

// define regridding method (lterpCube = assign interpolated values, cellRegridCube = assign nearest cell values )
#define PREPROC_REGRID lterpCube
// define method of getting value at (x,y) (getValue = interpolated, getCellValue = nearest cell value)
#define PREPROC_GETVAL getValue

int update_ip_files(double gtime);

int read_ip_var(gVar &vi, gVar &v, double gtime, int mode);
int read_all_ip_vars(double gtime, int mode);

int write_outvar(gVar &v, int istep);
int write_all_outvars(int istep);

//int write_state(string suffix);
int write_state(gVar &v, string suffix = "");
//int write_state(gVar &v, string filename = "");

#endif

