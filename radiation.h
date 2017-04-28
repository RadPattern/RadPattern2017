/**
 * @short   radiation Headers
 * @file    radiation.h
 * @author  Oluwaseun Fadugba, Yixin Zhang and Eric Jambo
 *
 *This function generates P-, SH- and SV-wave radiation patterns for 
 *single couple force, double couple, force dipole and point forces using
 *the values of theta and phi
*/  

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <cstdlib>
#include <cmath>

using namespace std;

int rad_pattern (double theta, double phi, radiation_pattern *radiation, Parameters *params);

int rad_pattern_x (double theta, double phi, radiation_pattern *radiation, Parameters *params);

int rad_pattern_y (double theta, double phi, radiation_pattern *radiation, Parameters *params);
