/**
 * @short   Prototypes
 * @file    prototypes.cpp
 * @author  Oluwaseun Fadugba, Yixin Zhang and Eric Jambo
 *
 * This file contains the prototypes and a short description of all the functions used in the RadPattern101 program.
 */

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <cstdlib>
#include <cmath>
#include "prototypes.h"
using namespace std;

// defining the value of PI.

#define PI 3.14159265


/**
 * Author:            Yixin Zhang
 *
 * Short description: This function float inputs
 *
 * Return             0 on sucess
 
 * Return             EXIT_FAILURE when the input file is missing
 */
 
// Just test, try to add read in program
// made by Yixin Zhang
//

int read_string(string str[], int i, string variable_name,
    string &string_value, ofstream &file)

{
	
    variable_name = str[i];
    file << variable_name << "=\t" << endl;
    cout << variable_name << "=\t" << endl;
    
    string_value = str[i+1];
    file << string_value << endl << "\n";
    cout << string_value << endl << "\n";
    
    return 0;

}

/**
 * Author:            Yixin Zhang
 *
 * Short description: This function float inputs
 *
 * Return             0 on sucess
 *
 */

int read_float(string str[], int i, string variable_name,
    float &float_value, ofstream &file)

{

    char* c;
    
    variable_name = str[i];
    file << variable_name << "=\t" << endl;
    cout << variable_name << "=\t" << endl;
    
    c = const_cast<char*>(str[i+1].c_str());
    sscanf(c,"%e", &float_value);
    file << float_value << endl << "\n";
    cout << float_value << endl << "\n";
    
    return 0;

}

// read variables in as int

int read_int(string str[], int i, string variable_name,
               int &int_value, ofstream &file)

{
    
    char* c;
    
    variable_name = str[i];
    file << variable_name << "=\t" << endl;
    cout << variable_name << "=\t" << endl;
    
    c = const_cast<char*>(str[i+1].c_str());
    sscanf(c,"%d", &int_value);
    file << int_value << endl << "\n";
    cout << int_value << endl << "\n";
    
    return 0;
    
}


//  End of test by Yixin Zhang

/**
 * Author:            Yixin Zhang
 *
 * Short description: This function gives the x,y coordinates for every point and returns the spherical coordinates for each grid
 *
 * Return             0 on sucess
 *
 */
 
double mesh_gen_o (float area_x, float area_y, int n_x, int n_y)

{
	
    cout <<  "Running: mesh_generator\n";
    cout << endl;
    
    // Perform operation.
    
    float X[n_x], Y[n_y];
    
    for (int i=0; i<n_x; i++)
    {
        
        X[i] = i * area_x / (n_x-1);
        cout << X[i] << "  ";
        
    }
    
    cout << "\n";
    
    for (int j=0; j<n_y; j++)
    {
        
        Y[j] = j * area_y / (n_y-1);
        cout << Y[j] << "  ";
        
    }
    
    cout << "\n";
  
    return 0;

}



/**
 * Author:            Yixin Zhang
 *
 * Short description: This function read the input waveforms in the program
 *
 * Return             0 on sucess
 * The distance between \f$(x_1,y_1)\f$ and \f$(x_2,y_2)\f$ is \f$\sqrt{(x_2-x_1)^2+(y_2-y_1)^2}\f$.
 */
 
int read_wave_file ()

{

    cout << "Running: read_wave_file\n";
    cout << endl;

    return 0;
  
}


 /**
 * Author:            Oluwaseun Fadugba
 *
 * Short description: This function calculates the P-wave Displacement for single force using the values of theta, phi, 
 *                    distance (R), moment (C1), P-wave velocity (alpha), density (rho), the input waveform (h) and the time series (time).
 *
 * Return             0 on sucess
 *
 */

int displ_pt_fo_Pw (
    double theta, double R, double moment, double alpha, 
    double rho, double *h, double *t, double *displ_pt_fo_Pw_o, int len )

{
    cout <<  "Running: displacement_point_force_P_wave\n";
    cout << endl;

    // Perform operation.
    
    for (int i=0; i<len; i++)
    {

        displ_pt_fo_Pw_o[i] = (cos(theta * PI / 180.0) * h[i]) / 
                              (4 * PI * rho * pow(alpha,2) * R);
    
    }

    for (int i=0; i<len; i++)
    {

        cout << "displacement_point_force_P_wave_output = " << displ_pt_fo_Pw_o[i] << endl;

    }
    
    
    return 0;

}

/**
 * Author:            Oluwaseun Fadugba
 *
 * Short description: This function calculates the S-wave Displacement for single force using the values of theta, phi,
 *                    distance (R), moment (C1), S-wave velocity (beta), density (rho), the input waveform (h) and the time series (time).
 *
 * Return             0 on sucess
 *
 */

double displ_pt_fo_Sw (
    double theta, double phi, double R, double moment, double beta, 
    double rho, double *h, double *t , double *displ_pt_fo_Sw_o, int len )

{
// Dont forget to apply time shift to the h and h derivatives	
    cout <<  "Running: displacement_point_force_S_wave\n";
    cout << endl;

    // Perform operation.
    
    for (int i=0; i<len; i++)
    {

        displ_pt_fo_Sw_o[i] = (-sin(theta * PI / 180.0) * h[i]) / 
                              (4 * PI * rho * pow(beta,2) * R);
    }

    for (int i=0; i<len; i++)
    {

        cout << "displacement_point_force_S_wave_output = " << displ_pt_fo_Sw_o[i] << endl;

    }

    return 0;

}

/**
 * Author:            Oluwaseun Fadugba
 *
 * Short description: This function calculates the P-wave Displacement for single force using the values of theta, phi,
 *                    distance (R), moment (C1), P-wave velocity (alpha), density (rho), derivative of the input waveform (h) and derivative of the time series (time). 
 *
 * Return             0 on sucess
 *
 */
 
double displ_si_fo_Pw (
    double theta, double phi, double R, double moment, double alpha, 
    double rho, double *h_der, double *t_der , double *displ_si_fo_Pw_o, int len )

{
    cout <<  "Running: displacement_single_force_P_wave\n";
    cout <<  endl;

    // Perform operation.
    
    for (int i=0; i<len; i++)
    {

        displ_si_fo_Pw_o[i] = (-sin(2.0*phi * PI / 180.0)*pow(sin(theta * PI / 180.0),2) * moment * h_der[i]) / 
                              (8.0 * PI * rho * pow(alpha,3) * R);
       
    }
    
    
    for (int i=0; i<len; i++)
    {
        cout << "displacement_single_force_P_wave_output = " << displ_si_fo_Pw_o[i] << endl;
    }


    return 0;

}


/**
 * Author:            Oluwaseun Fadugba
 *
 * Short description: This function calculates the SH-wave Displacement for single force using the values of theta, phi,
 *                    distance (R), moment (C1), distance (R), moment (C1), S-wave velocity (beta), density (rho), derivative 
 *                    of the input waveform (h) and derivative of the time series (time).
 * 
 * Return             0 on sucess
 *
 */

double displ_si_fo_SHw (
    double theta, double phi, double R, double moment, double beta, 
    double rho, double *h_der, double *t_der, double *displ_si_fo_SHw_o, int len  )

{
	
    cout <<  "Running: displacement_single_force_SH_wave\n";
    cout <<  endl;

    for (int i=0; i<len; i++)
    {
        displ_si_fo_SHw_o[i] = (sin(theta * PI / 180.0)*pow(sin(phi * PI / 180.0),2) * moment * h_der[i]) / 
                               (8.0 * PI * rho * pow(beta,3) * R);
        
    }
    
    for (int i=0; i<len; i++)
    {

        cout << "displacement_single_force_SH_wave_output = " << displ_si_fo_SHw_o[i] << endl;

    }

    return 0;

}


/**
 * Author:            Oluwaseun Fadugba
 *
 * Short description: This function calculates the SV-wave Displacement for single force using the values of theta, phi,
 *                    distance (R), moment (C1), distance (R), moment (C1), S-wave velocity (beta), density (rho), derivative of the input waveform (h) and derivative of the time series (time).
 *
 * Return             0 on sucess
 *
 */

double displ_si_fo_SVw (
    double theta, double phi, double R, double moment, double beta, 
    double rho, double *h_der, double *t_der, double *displ_si_fo_SVw_o, int len )

{

    cout <<  "Running: displacement_single_force_SV_wave\n";
    cout <<  endl;

    for (int i=0; i<len; i++)
    {
        displ_si_fo_SVw_o[i] = (-sin(2.0*theta * PI / 180.0)*(sin(2.0*phi * PI / 180.0)) * moment * h_der[i]) / 
                               (16.0 * PI * rho * pow(beta,3) * R);
        
    }
    
    for (int i=0; i<len; i++)
    {

        cout << "displacement_single_force_SV_wave_output = " << displ_si_fo_SVw_o[i] << endl;

    }

    return 0;

}


/**
 * Author:            Oluwaseun Fadugba
 *
 * Short description: This function calculates the P-wave Displacement for Double Couple using the values of theta, phi,
 *                    distance (R), moment (C1), P-wave velocity (alpha), density (rho), derivative of the input 
 *                    waveform (h) and derivative of the time series (time).
 *
 * Return             0 on sucess
 *
 */

double displ_do_cpl_Pw (
    double theta, double phi, double R, double moment, double alpha, 
    double rho, double *h_der, double *t_der, double *displ_do_cpl_Pw_o, int len  )

{

    cout <<  "Running: displacement_double_couple_P_wave\n";
    cout <<  endl;

    for (int i=0; i<len; i++)
    {

        displ_do_cpl_Pw_o[i] = (-sin(2.0*phi * PI / 180.0)*pow(sin(theta * PI / 180.0),2) * moment * h_der[i]) / 
                               (4.0 * PI * rho * pow(alpha,3) * R);
        
    }
    
    for (int i=0; i<len; i++)
    {

        cout << "displacement_double_couple_P_wave_output = " << displ_do_cpl_Pw_o[i] << endl;

    }


    return 0;

}

/**
 * Author:            Oluwaseun Fadugba
 *
 * Short description: This function calculates the SH-wave Displacement for Double Couple using the values of theta, phi,
 *                    distance (R), moment (C1), S-wave velocity (beta), density (rho), derivative of the input 
 *                    waveform (h) and derivative of the time series (time).
 *
 * Return             0 on sucess
 *
 */
 
double displ_do_cpl_SHw (
    double theta, double phi, double R, double moment, double beta, 
    double rho, double *h_der, double *t_der, double *displ_do_cpl_SHw_o, int len   )

{

    cout <<  "Running: displacement_double_couple_SH_wave\n";
    cout << endl;

    for (int i=0; i<len; i++)
    {

        displ_do_cpl_SHw_o[i] = (-sin(theta * PI / 180.0)*(cos(2.0*phi * PI / 180.0)) * moment * h_der[i]) / 
                                (4.0 * PI * rho * pow(beta,3) * R);
        
    }

    for (int i=0; i<len; i++)
    {

        cout << "displacement_double_couple_SH_wave_output = " << displ_do_cpl_SHw_o[i] << endl;

    }

    return 0;

}


/**
 * Author:            Oluwaseun Fadugba
 *
 * Short description: This function calculates the SV-wave Displacement for Double Couple using the values of theta, phi,
 *                    distance (R), moment (C1), S-wave velocity (beta), density (rho), derivative of the input 
 *                    waveform (h) and derivative of the time series (time).
 *
 * Return             0 on sucess
 *
 */

double displ_do_cpl_SVw (
    double theta, double phi, double R, double moment, double beta, 
    double rho, double *h_der, double *t_der, double *displ_do_cpl_SVw_o, int len   )

{

    cout <<  "Running: displacement_double_couple_SV_wave\n";
    cout <<  endl;

    for (int i=0; i<len; i++)
    {

        displ_do_cpl_SVw_o[i] = (-sin(2.0*theta * PI / 180.0)*(sin(2.0*phi * PI / 180.0)) * moment * h_der[i]) / 
                                (8.0 * PI * rho * pow(beta,3) * R);
       
    }

    for (int i=0; i<len; i++)
    {

         cout << "displacement_double_couple_SV_wave_output = " << displ_do_cpl_SVw_o[i] << endl;

    }
    return 0;

}



/**
 * Author:            Oluwaseun Fadugba
 *
 * Short description: This function calculates the P-wave Displacement for Force Dipole using the values of theta, phi,
 *                    distance (R), moment (C1), P-wave velocity (alpha), density (rho), input waveform (h) and time 
 *                    series (time).
 *
 * Return             0 on sucess
 *
 */

double displ_fo_dipo_Pw (
    double theta, double phi, double R, double moment, double alpha, 
    double rho, double *h, double *t, double *displ_fo_dipo_Pw_o, int len   )

{

    cout <<  "Running: displacement_force_dipole_P_wave\n";
    cout << endl;


    for (int i=0; i<len; i++)
    {

        displ_fo_dipo_Pw_o[i] = (pow(sin(theta * PI / 180.0),2)*pow(cos(phi * PI / 180.0),2) * moment * h[i]) / 
                                (4.0 * PI * rho * pow(alpha,3) * R);
        
    }

    for (int i=0; i<len; i++)
    {

        cout << "displacement_force_dipole_P_wave_output = " << displ_fo_dipo_Pw_o[i] << endl;

    }

    return 0;

}


/**
 * Author:            Oluwaseun Fadugba
 *
 * Short description: This function calculates the SH-wave Displacement for Force Dipole using the values of theta, phi,
 *                    distance (R), moment (C1), S-wave velocity, density (rho), input waveform (h) and time series (time). 
 *
 * Return             0 on sucess
 *
 */

double displ_fo_dipo_SHw (
    double theta, double phi, double R, double moment, double beta, 
    double rho, double *h, double *t, double *displ_fo_dipo_SHw_o, int len  )

{

    cout <<  "Running: displacement_force_dipole_SH_wave\n";
    cout << endl;

    for (int i=0; i<len; i++)
    {

        displ_fo_dipo_SHw_o[i] = (-sin(theta * PI / 180.0)*(sin(2.0*phi * PI / 180.0)) *
                                                      moment * h[i]) / 
                                                      (8.0 * PI * rho * pow(beta,3) * R);
    }


    for (int i=0; i<len; i++)
    {

        cout << "displacement_force_dipole_SH_wave_output = " << displ_fo_dipo_SHw_o[i] << endl;

    }

    return 0;

}



/**
 * Author:            Oluwaseun Fadugba
 *
 * Short description: This function calculates the SV-wave Displacement for Force Dipole using the values of theta, phi,
 *                    distance (R), moment (C1), S-wave velocity, density (rho), input waveform (h) and time series (time).
 *
 * Return             0 on sucess
 *
 */

double displ_fo_dipo_SVw (
    double theta, double phi, double R, double moment, double beta, 
    double rho, double *h, double *t, double *displ_fo_dipo_SVw_o, int len)

{

    cout <<  "Running: displacement_force_dipole_SV_wave\n";
    cout << endl;

    for (int i=0; i<len; i++)
    {

        displ_fo_dipo_SVw_o[i] = (cos(theta * PI / 180.0)*(cos(phi * PI / 180.0)) *
                                                      moment * h[i]) / 
                                                      (4.0 * PI * rho * pow(beta,3) * R);
        
    }

    for (int i=0; i<len; i++)
    {

        cout << "displacement_force_dipole_SV_wave_output = " << displ_fo_dipo_SVw_o[i] << endl;

    }


    return 0;

}


/**
 * Author:            Oluwaseun Fadugba
 *
 * Short description: This function writes the all the P-wave displacements, input waveform and its derivative and its
 *                    location (x,y) into a file
 *
 * Return             0 on sucess
 *
 */

int wr_Pw_2_file (
    double *displ_pt_fo_Pw_o, double *displ_si_cpl_Pw_o, 
    double *displ_do_cpl_Pw_o, double *displ_fo_dipo_Pw_o, 
    double *h, double *h_der, double *t, double *t_der, double xx, 
    double yy, string outputfilename, int len)

{

    cout <<  "Running: write_P_waves_to_file\n";
    cout << endl;

    ofstream fout(outputfilename.c_str(), ios::out);
    
    if (!fout)
    {

        cout << "Cannot create output file!" << endl;
        return 1;

    }


    int width = 10;
 
    // Set precision.
    //fout.precision(3);
    //fout.setf(ios::fixed);   fill by zeros if not enough number of digits after the decimal point

    // Start writing.
        
    // writing xx
    fout.width(width);		
    fout << "x  = " << xx << "\n";

    // writing yy
    fout.width(width);		
    fout << "y  = " << yy << "\n";

    // writing displacement_point_force_P_wave_output
       
    fout << "displacement_point_force_P_wave_output" << "\n";

    for (int i=0; i<len; i++)
    {
    
        fout.width(width);		
            
        fout << displ_pt_fo_Pw_o[i] << "\t";
		
    }
     
    fout << endl;

    // writing displacement_single_couple_P_wave_output

    fout << "displacement_single_couple_P_wave_output"<< "\n";
    
    for (int i=0; i<len; i++)
    {
    
        fout.width(width);		
        fout << displ_si_cpl_Pw_o[i] << "\t";
		
    }
   
    fout << endl;

    // writing displacement_double_couple_P_wave_output

    fout << "displacement_double_couple_P_wave_output" << "\n";
    
    for (int i=0; i<len; i++)
    {
    
        fout.width(width);		
        fout << displ_do_cpl_Pw_o[i] << "\t";
		
    }
    
    fout << endl;

    // writing displacement_force_dipole_P_wave_output
    fout << "displacement_force_dipole_P_wave_output" <<  "\n";
    
    for (int i=0; i<len; i++)
    {
    
        fout.width(width);		
        fout << displ_fo_dipo_Pw_o[i] << "\t";
		
    }
    
    fout << endl;

    fout.close();

    return 0;

}


/**
 * Author:            Oluwaseun Fadugba
 *
 * Short description: This function writes the all the SV-wave displacements, input waveform and its derivative and its
 *                    location (x,y) into a file
 *
 * Return             0 on sucess
 *
 */

int wr_SVw_2_file (
    double *displ_pt_fo_Sw_o, double *displ_si_cpl_SVw_o, 
    double *displ_do_cpl_SVw_o, double *displ_fo_dipo_SVw_o, 
    double *h, double *h_der, double *t, double *t_der, double xx, 
    double yy, string outputfilename, int len)
  
{

    ofstream fout(outputfilename.c_str(), std::ios_base::app); // append to the fout file

    if (!fout)
    {

        cout << "Cannot create output file!" << endl;
        return 1;

    }


    int width = 10;
 
    // Set precision.
    //fout.precision(3);
    //fout.setf(ios::fixed);  // fill by zeros if not enough number of digits after the decimal point

    // Start writing.
        
    // writing xx
    fout.width(width);		
    fout << "x  = " << xx << "\n";

    // writing yy
    fout.width(width);		
    fout << "y  = " << yy << "\n";

    // writing displacement_point_force_S_wave_output
       
    fout << "displacement_point_force_S_wave_output" << "\n";

    for (int i=0; i<len; i++)
    {
    
        fout.width(width);		
            
        fout << displ_pt_fo_Sw_o[i] << "\t";
		
    }
     
    fout << endl;

    // writing displacement_single_couple_SV_wave_output

    fout << "displacement_single_couple_SV_wave_output"<< "\n";
    
    for (int i=0; i<len; i++)
    {
    
        fout.width(width);		
        fout << displ_si_cpl_SVw_o[i] << "\t";
		
    }
   
    fout << endl;

    // writing displacement_double_couple_SV_wave_output

    fout << "displacement_double_couple_SV_wave_output" << "\n";
    
    for (int i=0; i<len; i++)
    {
    
        fout.width(width);		
        fout << displ_do_cpl_SVw_o[i] << "\t";
		
    }
    
    fout << endl;

    // writing displacement_force_dipole_SV_wave_output
    fout << "displacement_force_dipole_SV_wave_output" <<  "\n";
    
    for (int i=0; i<len; i++)
    {
    
        fout.width(width);		
        fout << displ_fo_dipo_SVw_o[i] << "\t";
		
    }
    
    fout << endl;

    fout.close();

    cout <<  "Running: write_SV_waves_to_file\n";
    cout << endl;

    return 0;

}


/**
 * Author:            Oluwaseun Fadugba
 *
 * Short description: This function writes the all the SH-wave displacements, input waveform and its derivative and its
 *                    location (x,y) into a file
 *
 * Return             0 on sucess
 *
 */

int wr_SHw_2_file (
    double *displ_pt_fo_Sw_o, double *displ_si_cpl_SHw_o, 
    double *displ_do_cpl_SHw_o, double *displ_fo_dipo_SHw_o, 
    double *h, double *h_der, double *t, double *t_der, double xx, 
    double yy, string outputfilename, int len)

{

    cout <<  "Running: write_SH_waves_to_file\n";
    cout << endl;

    ofstream fout(outputfilename.c_str(), std::ios_base::app);
    
    if (!fout)
    {

        cout << "Cannot create output file!" << endl;
        return 1;

    }


    int width = 10;
 
    // Set precision.
    //fout.precision(3);
    //fout.setf(ios::fixed);  // fill by zeros if not enough number of digits after the decimal point

    // Start writing.
        
    // writing xx
    fout.width(width);		
    fout << "x  = " << xx << "\n";

    // writing yy
    fout.width(width);		
    fout << "y  = " << yy << "\n";

    // writing displacement_point_force_S_wave_output
       
    fout << "displacement_point_force_S_wave_output" << "\n";

    for (int i=0; i<len; i++)
    {
    
        fout.width(width);		
            
        fout << displ_pt_fo_Sw_o[i] << "\t";
		
    }
     
    fout << endl;

    // writing displacement_single_couple_SH_wave_output

    fout << "displacement_single_couple_SH_wave_output"<< "\n";
    
    for (int i=0; i<len; i++)
    {
    
        fout.width(width);		
        fout << displ_si_cpl_SHw_o[i] << "\t";
		
    }
   
    fout << endl;

    // writing displacement_double_couple_SH_wave_output

    fout << "displacement_double_couple_SH_wave_output" << "\n";
    
    for (int i=0; i<len; i++)
    {
    
        fout.width(width);		
        fout << displ_do_cpl_SHw_o[i] << "\t";
		
    }
    
    fout << endl;

    // writing displacement_force_dipole_SH_wave_output
    fout << "displacement_force_dipole_SH_wave_output" <<  "\n";
   
    for (int i=0; i<len; i++)
    {
    
        fout.width(width);		
        fout << displ_fo_dipo_SHw_o[i] << "\t";
		
    }
    
    fout << endl;

    fout.close();

    return 0;

}


/**
 * Author:            Eric Jambo
 *
 * Short description: This function generates a guassian function using total time and time steps
 *
 * Return             0 on sucess
 *
 */

double gauss_func (double t[],double t_der[])

{

    cout << "Running: gaussian_function\n";
    cout << endl;

    return 0;

}


/**
 * Author:            Eric Jambo
 *
 * Short description: This function calculates the derivative of the input waveform using total time and time steps
 *
 * Return             0 on sucess
 *
 */

double der_wavf_func (double t[],double t_der[])

{

    cout << "Running: derivative_waveform_function\n";
    cout << endl;

    return 0;

}



/**
 * Author:            Eric Jambo
 *
 * Short description: This function converts the cartesian coordinates into spherical coordinates using location(x,y)
 *
 * Return             0 on sucess
 *
 */

double cart_2_sph (double *x,double *y)

{

    cout << "Running: cartesian_to_spherical_function\n";
    cout << endl;

    return 0;

}


/**
 * Author:            Eric Jambo
 *
 * Short description: This function generates P-wave radiation pattern for single couple force using the values of theta and phi
 *
 * Return             0 on sucess
 *
 */

double radp_Pw_si_cpl (double theta,double phi)

{

    cout << "Running: radiation_pattern_P_wave_single_couple_force\n";
    cout << endl;

    return 0;

}


/**
 * Author:            Eric Jambo
 *
 * Short description: This function generates SH-wave radiation pattern for single couple force using the values of theta and phi
 *
 * Return             0 on sucess
 *
 */

double radp_SHw_si_cpl (double theta,double phi)

{

    cout << "radiation_pattern_SH_wave_single_couple_force\n";
    cout << endl;

    return 0;

}



/**
 * Author:            Eric Jambo
 *
 * Short description: This function generates SV-wave radiation pattern for single couple force using the values of theta and phi
 *
 * Return             0 on sucess
 *
 */

double radp_SVw_si_cpl (double theta,double phi)

{

    cout << "radiation_pattern_SV_wave_single_couple_force\n";
    cout << endl;

    return 0;

}



/**
 * Author:            Eric Jambo
 *
 * Short description: This function generates P-wave radiation pattern for double couple force using the values of theta and phi
 *
 * Return             0 on sucess
 *
 */

double radp_Pw_do_cpl (double theta,double phi)

{

    cout << "Running: radiation_pattern_P_wave_double_couple_force\n";
    cout << endl;

    return 0;

}


/**
 * Author:            Eric Jambo
 *
 * Short description: This function generates SH-wave radiation pattern for double couple force using the values of theta and phi
 *
 * Return             0 on sucess
 *
 */

double radp_SHw_do_cpl (double theta,double phi)

{

    cout << "radiation_pattern_SH_wave_double_couple_force\n";
    cout << endl;

    return 0;

}


/**
 * Author:            Eric Jambo
 *
 * Short description: This function generates SV-wave radiation pattern for double couple force using the values of theta and phi
 *
 * Return             0 on sucess
 *
 */

double radp_SVw_do_cpl (double theta,double phi)

{

    cout << "Running: radiation_pattern_SV_wave_double_couple_force\n";
    cout << endl;

    return 0;

}


/**
 * Author:            Eric Jambo
 *
 * Short description: This function generates P-wave radiation pattern for force dipole using the values of theta and phi
 *
 * Return             0 on sucess
 *
 */

double radp_Pw_fo_dipo (double theta,double phi)

{

    cout << "Running: radiation_pattern_P_wave_force_dipole\n";
    cout << endl;

    return 0;

}



/**
 * Author:            Eric Jambo
 *
 * Short description: This function generates SH-wave radiation pattern for force dipole using the values of theta and phi
 *
 * Return             0 on sucess
 *
 */

double radp_SHw_fo_dipo (double theta,double phi)

{

    cout << "radiation_pattern_SH_wave_force_dipole\n";
    cout << endl;

    return 0;

}

/**
 * Author:            Eric Jambo
 *
 * Short description: This function generates SV-wave radiation pattern for force dipole using the values of theta and phi
 *
 * Return             0 on sucess
 *
 */

double radp_SVw_fo_dipo (double theta,double phi)

{

    cout << "radiation_pattern_SV_wave_force_dipole\n";
    cout << endl;

    return 0;

}


/**
 * Author:            Eric Jambo
 *
 * Short description: This function generates P-wave radiation pattern for point force using the values of theta and phi
 *
 * Return             0 on sucess
 *
 */

double radp_Pw_pt_fo (double theta,double phi)

{

    cout << "radiation_pattern_P_wave_point_force\n";
    cout << endl;

    return 0;

}


/**
 * Author:            Eric Jambo
 *
 * Short description: This function generates SH-wave radiation pattern for point force using the values of theta and phi
 *
 * Return             0 on sucess
 *
 */

double radp_SHw_pt_fo (double theta,double phi)

{

    cout << "radiation_pattern_SH_wave_point_force\n";
    cout << endl;

    return 0;

}


/**
 * Author:            Eric Jambo
 *
 * Short description: This function generates SV-wave radiation pattern for point force using the values of theta and phi
 *
 * Return             0 on sucess
 *
 */

double radp_SVw_pt_fo (double theta,double phi)

{

    cout << "radiation_pattern_SV_wave_point_force\n";
    cout << endl;

    return 0;

}
