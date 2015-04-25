// AUTHORS: Aaron Graham, Mike Jarrett
// PURPOSE: NERS 544 Course Project
// DATE   : April 3, 2015

#include<iostream>
#include<cstdlib>
#include<limits>
#include<cmath>
#include<vector>
#include<algorithm>

//Variables
const double eps = std::numeric_limits<double>::epsilon()*100.0;
const double nudge = 10.0*eps;
const double pi = 3.14159265358979;
const double neut_mass = 939.565378E6; // eV
const double kB = 8.6173324E-5; // eV K^-1

// Functions
double drand(void);
double Watt(void);
bool approxeq(double,double);
bool approxge(double,double);
bool approxle(double,double);
