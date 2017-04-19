/**
 * @mainpage RadPattern101
 *
 * Welcome to the Radiation Pattern Program (RadPattern101).
 * This program is a term project in Programming Tools course with Dr. Ricardo Taborda at
 * the University of Memphis in Spring 2017.
 *
 * The team members are Oluwaseun Fadugba, Yixin Zhang and Eric Jambo.
 *
 * In this project, we are going to generate the radiation pattern, wave propagation map,
 * and seismic waves at a specific location in the specified grid. This program requires
 * the type of force (point force, force dipole, double couple or single couple) that 
 * radiation pattern should be generated for. The program also requires the type of wave
 * e.g. P-, SV- or SH- waves and the rock properties through they propagate.
 *
 * Another aspect of the program is to generate a time-map of the wavepropagation with 
 * epicentral distance given an input waveform at the source. The default source input
 * waveform in this program is a Gaussian function.
 *
 * lastly, we also can display the waveform at a particular epicentral distance after
 * incorporating the influence of radiation pattern for the particular waveform as input.
 *
 */

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <cstdlib>
#include <cmath>
#include <time.h>

#include "struct.h"
#include "read_in.h"
#include "error_check.h"
#include "compute_displ.h"
#include "write2file.h"
#include "gaussian.h"
#include "radiation.h"
#include "cart2sph.h"

using namespace std;

#define PI 3.14159265

/**
 * @short   Main program
 * @file    main.cpp
 * @author  Oluwaseun Fadugba, Yixin Zhang and Eric Jambo
 * @param   parameters.in
 * @return  0 on success
 *
 */

int main(int argc, char* argv[])
{  
    clock_t t1 = clock(); //beginning time

    // Declare all parameters and files
    Parameters params; displacement displ;
    
    // Parameters will be read from input file, ckecked  for their reasonability, 
    // stored into memory, and then written into login file
    process_parameter(argc, argv, &params);
      
    int len = params.total_time/ params.time_step;

    double *t = new double[len];   init_time (t, &params);
    double *h = new double[len];   double *h_der = new double[len]; 
    double *rad_P = new double[1]; double *rad_SH = new double[1]; double *rad_SV = new double[1];
    double *R = new double[1];     double *theta = new double[1];  double *phi = new double[1]; 

    // This function generates a guassian function and its derivative
    gauss_func (params.total_time, params.time_step, h, h_der, len);
 
    // Now we want to iterate over the grid centers and determine the radiation 
    // pattern and displacement based on the type of force specified

    for (double xx=-(params.length_x / 2); xx<=(params.length_x / 2); xx+= params.length_x/(params.n_x - 1))  
    {   
        for (double yy=-(params.length_y / 2); yy<=(params.length_y / 2); yy+= params.length_y/(params.n_y - 1)) 
        {
            check_grid(xx, yy, &params);

            cart_2_sph (xx, yy, R, theta, phi);     

            compute_displ (R[1], theta[1], phi[1] , h, h_der, &displ, len, &params); 

            rad_patt (theta[1], phi[1], rad_P, rad_SH, rad_SV, len, params.force_type); 

            write_2_file (&displ, rad_P, rad_SH, rad_SV, t, xx, yy, 
                          "outputfilename", len);

        } // closing inner for loop
    } // outer for loop

    cout << "I have completed running all the functions.\n"
    "Good bye.\n"<< endl;
    std::cout << "Execution time : "<<float( clock () - t1 ) / CLOCKS_PER_SEC<< endl;

    return 0;
}
