/**
 * @short   Convert cartesian coordinates to spherical
 * @file    cart2sph.cpp
 *@author  Eric Jambo
 *
 * This file contains the function to convert the cartesian coordinates into spherical coordinates
 * using location(x,y) from the generated mesh grid
 *
 * Return             0 on sucess
 *
 */


#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <cstdlib>
#include <cmath>

#include "cart2sph.h"

using namespace std;

// defining the value of PI.

#define PI 3.14159265
const double Degrees = 180.0/PI; // Radians to Degrees

double cart_2_sph(double xx, double yy, double *R, double *theta, double *phi)
{
    double z = 0.0;

    if (xx !=0)
    {
    	R[1] = sqrt((pow(xx,2))+(pow(yy,2))+(pow(z,2))); // distance
        theta[1] = acos(z/R[1])* Degrees;
    	phi[1] = atan2(yy, xx) * Degrees;
        
        return 0;
   }
   
    return false;
}
