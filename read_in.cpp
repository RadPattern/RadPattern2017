/**
 * @short   read in parameters
 * @file    read_in.cpp
 * @author  Yixin Zhang
 *
 * This file contains the prototypes for reading parameters from the input file 
 * and write into memory and login file.
 */

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <cstdlib>
#include <cmath>

#include "struct.h"
#include "read_in.h"
#include "error_check.h"

/**
 * Author:            Yixin Zhang
 *
 * Short description: This function read all parameters in from parameter.in,
 * check the reasonability, and then store them into login file, out to the screen
 * as well as to the program memory.
 *
 * Return             0 on sucess
 *
 * Return             EXIT_FAILURE when parameter.in is wrong.
 */

int process_parameter(int argc, char* argv[], Parameters *params)
{
    //Parameters *params
    string str[100];

    // Check the number of input files, exit program if check failed
    check_file_num(argc);

    // Read parameters as strings from input file
    read_infile(argc, argv, str);

    // Store parameters into memory as their declared types
    store_parameters(str, params);

    // Check the reasonability of variables
    check_variables(params);

    // Write parameters into login file
    out_logfile(params);

    return 0;
}

/**
 * Author:            Yixin Zhang
 *
 * Short description: This function read information from infile and store them as strings.
 *
 * Return             0 on sucess
 *
 * Return             EXIT_FAILURE when parameter.in is wrong.
 */

int read_infile(int argc, char* argv[], string str[])
{
    int i = 0;
    string s;
    ifstream infile;

    infile.open(argv[1]);

    check_file_open(infile);

    // Read variables as strings from infile
    while( infile >> s )
    {
        str[i] = s;
        i++;
    }

    infile.close();

    return 0;
}

/**
 * Author:            Yixin Zhang
 *
 * Short description: This function read parameters in memory as their declared types.
 *
 * Return             0 on sucess
 *
 * Return             EXIT_FAILURE when parameter.in is wrong.
 */

int store_parameters(string str[], Parameters *params)
{
    cout << "Start reading variables to memory\n" << endl;

    int i = 0;

    while(str[i] != "")
    {
        if(str[i] == "model_name") {read_string(str, i, params->model_name);}

        if(str[i] == "force_type") {read_string(str, i, params->force_type);}

        if(str[i] == "alpha") {read_float(str, i, params->alpha);}

        if(str[i] == "waveform") {read_string(str, i, params->waveform);}

        if(str[i] == "beta") {read_float(str, i, params->beta);}

        if(str[i] == "time_step") {read_float(str, i, params->time_step);}

        if(str[i] == "total_time") {read_float(str, i, params->total_time);}

        if(str[i] == "length_x") {read_float(str, i, params->length_x);}

        if(str[i] == "length_y") {read_float(str, i, params->length_y);}

        if(str[i] == "n_x") {read_int(str, i, params->n_x);}

        if(str[i] == "n_y") {read_int(str, i, params->n_y);}

        if(str[i] == "moment") {read_float(str, i, params->moment);}

        if(str[i] == "density") {read_float(str, i, params->rho);}

        if(str[i] == "logfile_path") {read_string(str, i, params->logfile_path);}

        if(str[i] == "outputfile_path") {read_string(str, i, params->outputfile_path);}

        i++;
    }

    cout << "Finish reading parameters.\n" << endl;

    return 0;
}

/**
 * Author:            Yixin Zhang
 *
 * Short description: This function used for reading string variables in memory, print to screen
 *
 * Return             0 on sucess
 *
 */

int read_string(string str[], int i, string &string_value)
{
    string variable_name;

    variable_name = str[i];
    cout << variable_name << "=\t" << endl;

    string_value = str[i+1];
    cout << string_value << endl << "\n";

    return 0;
}

/**
 * Author:            Yixin Zhang
 *
 * Short description: This function used for reading float variables in memory, print to screen
 *
 * Return             0 on sucess
 *
 */

int read_float(string str[], int i, float &float_value)
{
    char* c;
    string variable_name;

    variable_name = str[i];
    cout << variable_name << "=\t" << endl;

    c = const_cast<char*>(str[i+1].c_str());
    sscanf(c, "%e", &float_value);
    cout << float_value << endl << "\n";
    return 0;
}

/**
 * Author:            Yixin Zhang
 *
 * Short description: This function used for reading int variables in memory, print to screen
 *
 * Return             0 on sucess
 *
 */

int read_int(string str[], int i, int &int_value)
{
    char* c;
    string variable_name;

    variable_name = str[i];
    cout << variable_name << "=\t" << endl;

    c = const_cast<char*>(str[i+1].c_str());
    sscanf(c,"%d", &int_value);
    cout << int_value << endl << "\n";

    return 0;
}

/**
 * Author:            Yixin Zhang
 *
 * Short description: This function used for writing values into login file
 *
 * Return             0 on sucess
 *
 */

int out_logfile(Parameters *params)
{
    const char* path = const_cast<char*>(params->logfile_path.c_str());
    ofstream logfile(path);

    logfile << "model_name =\t" << params->model_name << endl;
    logfile << "force_type =\t" << params->force_type << endl;
    logfile << "alpha =\t" << params->alpha << endl;
    logfile << "beta =\t" << params->beta << endl;
    logfile << "time_step =\t" << params->time_step << endl;
    logfile << "total_time =\t" << params->total_time << endl;
    logfile << "waveform =\t" << params->waveform << endl;
    logfile << "length_x =\t" << params->length_x << endl;
    logfile << "length_y =\t" << params->length_y << endl;
    logfile << "n_x =\t" << params->n_x << endl;
    logfile << "n_y =\t" << params->n_y << endl;
    logfile << "moment =\t" << params->moment << endl;
    logfile << "density =\t" << params->rho << endl;
    logfile << "logfile_path =\t" << params->logfile_path << endl;
    logfile << "outputfile_path =\t" << params->outputfile_path << endl;

    cout << "login file has been written.\n" <<endl;

    logfile.close();
    
    return 0;
}
