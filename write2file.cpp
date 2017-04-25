/**
 * @short   Write
 * @file    write2file.cpp
 * @author  Oluwaseun Fadugba
 *
 * This file contains write_2_file function
 */

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <cstdlib>
#include <cmath>

#include <string.h>
#include "write2file.h"
#include "error_check.h"

using namespace std;

/**
 * Author:            Oluwaseun Fadugba
 *
 * Short description: This function writes the all the displacements,
 *                    radiation and its location (x,y) into a file
 *
 * Return             0 on sucess
 *
 *                    Exit program when output file cannot be created
 *
 */

int write_2_file (displacement *displ, radiation_pattern *radiation,
                  double *t, double xx, double yy, string outputfilename, int len)
{
    // write radiation pattern
    write_rad_patt (radiation, xx, yy, outputfilename);

    // write displacement
    write_displ (displ, t, xx, yy, outputfilename, len);

    return 0;
}

/**
 * Author:            Oluwaseun Fadugba
 *
 * Short description: 
 * This function writes the displacements and its location (x,y) into a file
 *
 * Return             0 on sucess
 *
 *                    Exit program when output file cannot be created
 *
 */

int write_displ (displacement *displ, double *t, double xx, double yy, string outputfilename, int len)
{
    char str[80]; strcpy (str,outputfilename.c_str()); strcat (str,"_displ.txt");

    ofstream fout(str, ios::ios_base::app);

    check_outfile(!fout, 1);

    int width = 10; // allowable width of each field

    // writing xx and yy
    fout.width(width);
    fout << "x  = " << xx << "\t" << "y  = " << yy << "\n";

    // writing displacement results
    fout << "         t   "  "\t" << " displ_P "  "\t" << " displ_SH "  "\t"
         << " displ_SV " << "\n";

    for (int i=1; i<=len; i++)
    {
        fout.width(width);
        fout << t[i] << "\t" << displ->P[i] << "\t" << displ->SH[i] << "\t"
             << displ->SV[i] << "\n";
    }

    fout << "\n"; fout << endl; 
    fout.close();
    
    return 0;
}

/**
 * Author:            Oluwaseun Fadugba
 *
 * Short description: This function writes the radiation pattern results and 
 *                    its location (x,y) into a file
 *
 * Return             0 on sucess
 *
 *                    Exit program when output file cannot be created
 *
 */

int write_rad_patt (
    radiation_pattern *radiation, double xx, double yy, string outputfilename)
{
    char str[80]; strcpy (str, outputfilename.c_str()); strcat (str,"_rad_patt.txt");

    ofstream fout(str, ios::ios_base::app);

    check_outfile(!fout, 2);

    int width = 10; // allowable width of each field

    // writing xx and yy
    fout.width(width);
    fout << xx << "\t" << yy << "\t" << radiation->P_x << "\t" << radiation->P_y <<  "\t" 
         << radiation-> SH_x << "\t"<< radiation-> SH_y << "\t"
         << radiation-> SV_x << "\t"<< radiation-> SV_y << "\n";

    fout.close();

    return 0;
}
