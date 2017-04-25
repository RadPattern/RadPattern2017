/**
 * @short   gaussian
 * @file    gaussian.cpp
 * @author  Eric Jambo
 *
 * short description:This file contains the functions used to generate a guassian
 *                   function and its derivative using total time and time steps.
 *
 * Return   0 on sucess
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
#include "gaussian.h"

using namespace std;

// defining the value of PI.

#define PI 3.14159265

void gauss_func (double *h, double *h_der, int len, Parameters *params)
{
    
    double total_time = params->total_time; 
    double time_step = params->time_step;

    double a = 2.2/ total_time;  //Mean of distribution
    // set standard deviation
    double sigma = 0.75;
    double s = 2.0 * sigma * sigma;
    int i = 0;

    // generate gaussian time function and its derivative
    for (double time = 0; time <= total_time; time += time_step) 	
    {
         i++;
         if(time <=10*sigma)
         {  
            h[i] = (1/sqrt(PI * s))*(exp(-(time-a)*(time-a)/s));
            h_der[i] = -(time - a)/((s*s*s)*sqrt(2.0*PI))*exp(-(time -a)*(time-a)/s);
         }
         else{
              h[i]=0.0;
              h_der[i] =0.0;
             }
    }
    return;
}
