/**
 * @short   gaussian Header
 * @file    gaussian.h
 * @author  Eric Jambo
 *
 * This file contains the prototypes and a short description of all the functions used
 * to generate a gaussian function and its derivative using total time and time steps.
 */
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <cmath>
#include <cstdlib>

using namespace std;

void gauss_func (double *h, double *h_der, int len, Parameters *params);

