/**
 * @short   Read input Headers
 * @file    read_in.h
 * @author  Yixin Zhang
 *
 * This file contains the prototypes all functions used for reading parameters in.
 *
 */

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <cmath>
#include <cstdlib>

#include "struct.h"

using namespace std;

// Short description: These functions read the input parameters in the program
// ----------------------------------------------------------------------------

int process_parameter(int argc, char* argv[], Parameters *params);

int read_infile(int argc, char* argv[], string str[]);

int store_parameters(string str[], Parameters *params);

int read_string(string str[], int i, string &string_value);

int read_float(string str[], int i, float &float_value);

int read_int(string str[], int i, int &int_value);

int out_logfile(Parameters *params);
