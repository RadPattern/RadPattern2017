/**
 * @short   Convert cartesian coordinates to spherical
 * @file    cart2sph.h
 *@author  Eric Jambo
 *
 * This file contains the function to convert the cartesian coordinates into spherical coordinates
 * using location(x,y) from the generated mesh grid
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

// defining the value of PI.
#define PI 3.14159265

// Short description: This function converts the cartesian coordinates into spherical
// coordinates using location(x,y)
// ----------------------------------------------------------------------------------

double cart_2_sph (double xx, double yy, double *R, double *theta, double *phi);
