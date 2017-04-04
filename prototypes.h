/**
 * @short   Prototype Headers
 * @file    prototypes.h
 * @author  Oluwaseun Fadugba, Yixin Zhang and Eric Jambo
 *
 * This file contains the prototypes and a short description of all the functions used in the RadPattern101 program.
 */
 
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <cmath>
#include <cstdlib>

using namespace std;
// defining the value of PI.

#define PI 3.14159265


//
// Reading in the parameters from a file
//
// Author:            Yixin Zhang
//
// Short description: These functions read the input waveforms in the program
//
// -----------------------------------------------------------------------------------------------
//	

int read_string(string str[], int i, string variable_name,
    string &string_value, ofstream &file);

int read_float(string str[], int i, string variable_name,
    float &float_value, ofstream &file);

int read_int(string str[], int i, string variable_name,
             int &int_value, ofstream &file);

int read_wave_file ();


//
// Function:          gaussian_function and its derivative
//
// Author:            Eric Jambo
//
// Short description: These functions generates a guassian function and its derivative using total time and time steps
//
// -----------------------------------------------------------------------------------------------
//


double gauss_func (double t[],double t_der[]);

double der_wavf_func (double t[],double t_der[]);



//
// Function:          Mesh_generator
//
// Author:            Yixin Zhang
//
// Short description: This function gives the x,y coordinates for every point and returns the spherical coordinates for each grid
//
// -----------------------------------------------------------------------------------------------
//

double mesh_gen_o (float area_x, float area_y, int n_x, int n_y);


//
// Function:          cartesian_to_spherical_function
//
// Author:            Eric Jambo
//
// Short description: This function converts the cartesian coordinates into spherical coordinates using location(x,y)
//
// -----------------------------------------------------------------------------------------------
//


double cart_2_sph (double *x,double *y);


//
// Functions:         Radiation_pattern
//
// Author:            Eric Jambo
//
// Short description: This function generates P-, SH- and SV-wave radiation patterns for single couple force, double couple, 
//                    force dipole and point forces using the values of theta and phi
//
// -----------------------------------------------------------------------------------------------
//

double radp_Pw_si_cpl (double theta,double phi);

double radp_SHw_si_cpl (double theta,double phi);

double radp_SVw_si_cpl (double theta,double phi);

double radp_Pw_do_cpl (double theta,double phi);

double radp_SHw_do_cpl (double theta,double phi);

double radp_SVw_do_cpl (double theta,double phi);

double radp_Pw_fo_dipo (double theta,double phi);

double radp_SHw_fo_dipo (double theta,double phi);

double radp_SVw_fo_dipo (double theta,double phi);

double radp_Pw_pt_fo (double theta,double phi);

double radp_SHw_pt_fo (double theta,double phi);

double radp_SVw_pt_fo (double theta,double phi);


//
// Functions:         displacement
//
// Author:            Oluwaseun Fadugba
//
// Short description: This function calculates the P-, SH- and SH-wave Displacements for single force, double couple, force dipole and point forces 
//                    using the values of theta, phi, distance (R), moment (C1), S-wave velocity (beta), density (rho), the input waveform (h) and its derivative, and the time series (time) and its derivative.
//
// -----------------------------------------------------------------------------------------------
//	


int displ_pt_fo_Pw (
    double theta, double R, double moment, double alpha, 
    double rho, double *h, double *t, double *displ_pt_fo_Pw_o, int len);

double displ_pt_fo_Sw (
    double theta, double phi, double R, double moment, double beta, 
    double rho, double *h, double *t, double *displ_pt_fo_Sw_o, int len);

double displ_si_fo_Pw (
    double theta, double phi, double R, double moment, double alpha, 
    double rho, double *h_der, double *t_der, double *displ_si_fo_Pw_o, int len);

double displ_si_fo_SHw (
    double theta, double phi, double R, double moment, double alpha, 
    double rho, double *h_der, double *t_der, double *displ_si_fo_SHw_o, int len);

double displ_si_fo_SVw (
    double theta, double phi, double R, double moment, double beta, 
    double rho, double *h_der, double *t_der, double *displ_si_fo_SVw_o, int len);

double displ_do_cpl_Pw (
    double theta, double phi, double R, double moment, double alpha, 
    double rho, double *h_der, double *t_der, double *displ_do_cpl_Pw_o, int len);

double displ_do_cpl_SHw (
    double theta, double phi, double R, double moment, double alpha, 
    double rho, double *h_der, double *t_der, double *displ_do_cpl_SHw_o, int len);

double displ_do_cpl_SVw (
    double theta, double phi, double R, double moment, double alpha, 
    double rho, double *h_der, double *t_der, double *displ_do_cpl_SVw_o, int len);

double displ_fo_dipo_Pw (
    double theta, double phi, double R, double moment, double alpha, 
    double rho, double *h, double *t,  double *displ_fo_dipo_Pw_o, int len);

double displ_fo_dipo_SHw (
    double theta, double phi, double R, double moment, double alpha, 
    double rho, double *h, double *t,  double *displ_fo_dipo_SHw_o, int len);

double displ_fo_dipo_SVw (
    double theta, double phi, double R, double moment, double alpha, 
    double rho, double *h, double *t,  double *displ_fo_dipo_SVw_o, int len);

		
		
		
//
// Functions:         write_results_to_file
//
// Author:            Oluwaseun Fadugba
//
// Short description: This function writes the all the P-, SH- and SV-wave displacements, input waveform and its derivative and its location (x,y) into a file
//                   
// -----------------------------------------------------------------------------------------------
//			
		

int wr_Pw_2_file (
    double *displ_pt_fo_Pw_o, double *displ_si_cpl_Pw_o, 
    double *displ_do_cpl_Pw_o, double *displ_fo_dipo_Pw_o, 
    double *h, double *h_der, double *t, double *t_der, double xx, 
    double yy, string outputfilename , int len);

int wr_SVw_2_file (
    double *displ_pt_fo_Sw_o, double *displ_si_cpl_SVw_o, 
    double *displ_do_cpl_SVw_o, double *displ_fo_dipo_SVw_o, 
    double *h, double *h_der, double *t, double *t_der, double xx, 
    double yy, string outputfilename, int len);

int wr_SHw_2_file (
    double *displ_pt_fo_Sw_o, double *displ_si_cpl_SHw_o, 
    double *displ_do_cpl_SHw_o, double *displ_fo_dipo_SHw_o, 
    double *h, double *h_der, double *t, double *t_der, double xy, 
    double yy, string outputfilename, int len);

