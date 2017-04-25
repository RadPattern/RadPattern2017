/**
 * @short   Compute displacements
 * @file    compute_displ.cpp
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

#include "read_in.h"
#include "struct.h"
#include "compute_displ.h"
#include "error_check.h"

using namespace std;

#define PI 3.14159265

/**
 * Author:            Oluwaseun Fadugba
 *
 * Short description: This function determines the displacement for P, SH and SV waves
 *
 * Return             0 on sucess
 *
 *                    Exit program on issues with coordinate
 *
 */

int compute_displ (double R, double theta, double phi, double *h, double *h_der,
                   displacement *displ, int len, Parameters *params)
{
    // Check variables
    check_loc(R, theta, phi);

    // displacement for P and SH-SV waves. Time shift is included.
    compute_displ_P (R, theta, phi , h, h_der, displ, len, params);

    compute_displ_SH_SV (R, theta, phi , h, h_der, displ, len, params);

    return 0;
}

/**
 * Author:            Oluwaseun Fadugba
 *
 * Short description: This function determines the displacement for P-wave after
 * applying necessary time shift to the guassian function
 *
 * Return             0 on sucess
 *
 */

int compute_displ_P (double R, double theta, double phi, double *h, double *h_der,
                     displacement *displ, int len, Parameters *params)
{
    // making parameters local to avoid multiple "params->"
    double moment = params->moment; double alpha = params->alpha;
    double rho = params->rho; double time_step = params->time_step;
    string force_type = params->force_type;

    // determine time shift
    int shift_P = int(round(R / (alpha * time_step)));

    // leading zeros due to the time shift
    for (int j = 0; j< shift_P; j++)
    {
        displ->P[j] = 0.0;
    }

    // Non-zero displacements
    for (int i = 0; i < (len - shift_P); i++)
    {
        if (force_type == "point_force")
        {
            displ->P[i+shift_P] = ((cos(theta) * PI / 180.0) * h[i]) / (4 * PI * rho * pow(alpha, 2) * R);
        }

        if (force_type == "single_couple")
        {
            displ->P[i+shift_P] = (-sin(2.0*phi * PI / 180.0)*pow(sin(theta * PI / 180.0), 2) *
                                   moment * h_der[i]) / (8.0 * PI * rho * pow(alpha, 3) * R);
        }

        if (force_type == "double_couple")
        {
            displ->P[i+shift_P]  = (-sin(2.0*phi * PI / 180.0)*pow(sin(theta * PI / 180.0), 2) *
                                    moment * h_der[i]) / (4.0 * PI * rho * pow(alpha, 3) * R);
        }

        if (force_type == "force_dipole")
        {
            displ->P[i+shift_P]  = (pow(sin(theta * PI / 180.0),2)*pow(cos(phi * PI / 180.0), 2) *
                                    moment * h[i]) / (4.0 * PI * rho * pow(alpha, 3) * R);
        }
    }

    return 0;
}

/**
 * Author:            Oluwaseun Fadugba
 *
 * Short description: This function determines the displacement for SH- and SV- waves after
 * applying necessary time shift to the guassian function
 *
 * Return             0 on sucess
 *
 */

int compute_displ_SH_SV (double R, double theta, double phi, double *h, double *h_der,
                         displacement *displ, int len, Parameters *params)
{
    // making parameters local to avoid multiple "params->"
    double moment = params->moment; double beta = params->beta;
    double rho = params->rho; double time_step = params->time_step;
    string force_type = params->force_type;

    // determine time shift
    int shift_S = int(round(R / (beta * time_step)));

    // leading zeros due to the time shift
    for (int j=0; j<shift_S; j++)
    {
        displ->SH[j] = 0.0;
        displ->SV[j] = 0.0;
    }

    // Non-zero displacements
    for (int i=0; i<len; i++)
    {
        if (force_type == "point_force")
        {
            displ->SH[i+shift_S] = (-sin(theta * PI / 180.0) * h[i]) / (4 * PI * rho * pow(beta, 2) * R);
            displ->SV[i+shift_S] = displ->SH[i+shift_S];
        }
        if (force_type == "single_couple")
        {
            displ->SH[i+shift_S] = (sin(theta * PI / 180.0)*pow(sin(phi * PI / 180.0), 2) *
                                    moment * h_der[i]) / (8.0 * PI * rho * pow(beta, 3) * R);
            
            displ->SV[i+shift_S] = (-sin(2.0*theta * PI / 180.0)*(sin(2.0*phi * PI / 180.0)) *
                                    moment * h_der[i]) / (16.0 * PI * rho * pow(beta, 3) * R);
        }
        if (force_type == "double_couple")
        {
            displ->SH[i+shift_S] = (-sin(theta * PI / 180.0)*(cos(2.0*phi * PI / 180.0)) *
                                    moment * h_der[i]) / (4.0 * PI * rho * pow(beta, 3) * R);
            displ->SV[i+shift_S] = (-sin(2.0*theta * PI / 180.0)*(sin(2.0*phi * PI / 180.0)) *
                                    moment * h_der[i]) / (8.0 * PI * rho * pow(beta, 3) * R);
        }
        if (force_type == "force_dipole")
        {
            displ->SH[i+shift_S] = (-sin(theta * PI / 180.0)*(sin(2.0*phi * PI / 180.0)) *
                                    moment * h[i]) / (8.0 * PI * rho * pow(beta, 3) * R);
            displ->SV[i+shift_S] = (cos(theta * PI / 180.0)*(cos(phi * PI / 180.0)) *
                                    moment * h[i]) / (4.0 * PI * rho * pow(beta, 3) * R);
        }
    }

    return 0;
}
