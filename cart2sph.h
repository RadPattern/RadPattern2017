/**
 * @short   cartesian coordinates to spherical header
 * @file    cart2sph.h
 * @author  Eric Jambo
 *
 * This file contains the function to convert the cartesian coordinates into spherical coordinates
 * using location(x,y)
 *
 */
 
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <cmath>
#include <cstdlib>
#include <vector>

using namespace std;

double cart_2_sph (double xx, double yy, double *R, double *theta, double *phi);
