/**
 * @short   Initializes time array
 * @file    init_time.cpp
 * @author  Oluwaseun
 *
 * This file contains functions compute_displ.cpp and its helper functions,
 * and init_time.cpp to calculate the dispalcement given necessary parameters
 * and initialize the time vector.
 *
 */

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <cstdlib>
#include <cmath>

#include "struct.h"
#include "error_check.h"

using namespace std;

#define PI 3.14159265

/**
 * Author:            Oluwaseun Fadugba
 *
 * Short description: This function initializes the time vector based on total time and time step
 *
 * Return             0 on sucess
 *
 *                    Exit program on failure
 *
 */

int init_time (double *t, double len, Parameters *params)
{
    int i = 0;
    
    for (double tt=0; tt <= (params->total_time - params->time_step); tt+= params->time_step)
    {
        i++;
        t[i]  = tt;
    }

    // Check the length of the time array
    check_t_len(len, i);

    return 0;
}
