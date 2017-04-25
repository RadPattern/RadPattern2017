/**
 * @short   error checking
 * @file    error_check.cpp
 * @author  Yixin Zhang
 *
 * This file contains the prototypes and a short description of functions used for checking errors.
 */

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <cstdlib>
#include <cmath>

#include "read_in.h"
#include "struct.h"

using namespace std;

/**
 * Author:            Yixin Zhang
 *
 * Short description: This function checks the number of input file(s).
 *
 * Return             0 on sucess
 *
 * Return             1 on fail
 */

int check_file_num(int argc)
{
    if(argc == 2)
    {
        cout << "//***My program is initializing***// \n" << endl;
        return 0;
    }
    else
    {
        cout << "Redunant inputfile(s) or Missing inputfile, please check\n" << endl;
        exit(EXIT_FAILURE);
    }
}

/**
 * Author:            Yixin Zhang
 *
 * Short description: This function checks if the input file can be open.
 *
 * Return             0 on sucess
 *
 * Return             1 on fail
 */

int check_file_open(ifstream &infile)
{
    if (infile.is_open())
    {
        cout << "File successfully opened \n" << endl;
        return 0;
    }
    else
    {
        cout << "Error opening input file\n" << endl;
        exit(EXIT_FAILURE);
    }
}

/**
 * Author:            Yixin Zhang
 *
 * Short description: This function checks the reasonability of readin values.
 *
 * Return             0 on sucess
 *
 * Return             1 on fail
 */

int check_variables(Parameters *params)
{
    if( (params->alpha > 0) && (params->beta > 0) && (params->alpha > params->beta) && (params->time_step >0)
        && (params->total_time >0) && (params->length_x > 0) && (params->length_y > 0) && (params->n_x > 0)
        && (params->n_y > 0) && (params->moment > 0) && (params->rho > 0) )
    {
        cout << "Variables have been checked, parameters are good to use.\n" << endl;
        return 0;
    }
    else
    {
        cout << "Input file contains illegal variable(s), please check.\n" << endl;
        exit(EXIT_FAILURE);
    }
}

/**
 * Author:            Oluwaseun Fadugba
 *
 * Short description: This function checks the reasonability of the grids.
 *
 */

int check_grid(double xx, double yy, Parameters *params)
{
    if (xx < -(params->length_x / 2) || xx > (params->length_x / 2) + 1 ||
        yy < -(params->length_y / 2) || yy > (params->length_y / 2) + 1)
    {
        cout << "---------------------Error! ------------------------------------- \n";
        cout << "----------------Invalid grid centers!---------------------------- \n";
        cout << "----------------------------------------------------------------- \n";
        cout << endl;
        exit(EXIT_FAILURE);
    }  

    return 0;
}

/**
 * Author:            Oluwaseun Fadugba
 *
 * Short description: 
 * This function checks the reasonability of the locations inputs to compute_displ.cpp
 *
 */

int check_loc(double R, double theta, double phi)
{
    if (R < 0.0 || abs(theta) > 91 || abs(phi) > 181 )
    { 
        cout << "--------------------------Error! ------------------------------ \n";
        cout << "Invalid spherical location for displacement field calculations! \n";
        cout << "--------Check the consistency of the domain geometry----------- \n";
        cout << "--------------------------------------------------------------- \n";
        cout << endl;
        exit(EXIT_FAILURE);
    }

    return 0;
}

/**
 * Author:            Oluwaseun Fadugba
 *
 * Short description: 
 * This function checks the length of the time array
 *
 */

int check_t_len(double len, int i)
{
    if (len != i)
    { 
        cout << "--------------------------Error! ------------------------------ \n";
        cout << "-------------Error initializing time array--------------------- \n";
        cout << "--Required total time should be divisible by time step--------- \n";
        cout << "--------------------------------------------------------------- \n";
        cout << endl;
        exit(EXIT_FAILURE);
    }

    return 0;
}

/**
 * Author:            Oluwaseun Fadugba
 *
 * Short description: This function checks if the output files are successfully opened.
 *
 */

int check_outfile(int index, int str)
{
    if (str == 1)
    {
        if (index != 0)
        {
            cout << "--------------------------Error! ------------------------------ \n";
            cout << "---------Cannot create displacement output file!--------------- \n";
            cout << "---Check output filename for spaces or special characters------ \n";
            cout << "--------------------------------------------------------------- \n";
            cout << endl;

            exit(EXIT_FAILURE);
        }
    }

    if (str == 2)
    {
        if (index != 0)
        {
            cout << "--------------------------Error! ------------------------------ \n";
            cout << "---------Cannot create radiation pattern output file!---------- \n";
            cout << "---Check output filename for spaces or special characters------ \n";
            cout << "--------------------------------------------------------------- \n";
            cout << endl;

            exit(EXIT_FAILURE);
        }
    }

    return 0;
}
