/**
 * @short   Compute displacement Headers
 * @file    compute_displ.h
 * @author  Oluwaseun Fadugba
 *
 * This file contains the prototypes compute_displ.cpp and init_time.cpp
 *
 */

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <cmath>
#include <cstdlib>

using namespace std;

int compute_displ (double RR, double theta, double phi, double *h, double *h_der,
                   displacement *displ, int len, Parameters *params);

int compute_displ_P (double R, double theta, double phi, double *h, double *h_der,
                     displacement *displ, int len, Parameters *params);

int compute_displ_SH_SV (double R, double theta, double phi, double *h, double *h_der,
                         displacement *displ, int len, Parameters *params);
