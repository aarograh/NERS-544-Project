// AUTHORS: Aaron Graham, Mike Jarrett
// PURPOSE: NERS 544 Course Project
// DATE   : April 30, 2015

#include<iostream>
#include<cstdlib>
#include<limits>
#include<cmath>
#include<vector>
#include<algorithm>

//Variables
const double eps = std::numeric_limits<double>::epsilon()*100.0;
const double nudge = 1.0e-7;
const double pi = 3.14159265358979;
const double neut_mass = 939.565378; // MeV
const double kB = 8.6173324E-11; // MeV K^-1
const double nu = 2.45; // neutrons per fission
const double temp = 293; // Kelvin
const double lightspeed = 299792458.0; // m/s

// Functions
double drand(void);
double Watt(void);
bool approxeq(double, double);
bool approxge(double, double);
bool approxle(double, double);
bool softeq(double, double, double);
