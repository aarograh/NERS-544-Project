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

// Functions
double drand(void);
double Watt(void);
double calcEntropy(int,double(*)[3]);
bool approxeq(double,double);
bool approxge(double,double);
bool approxle(double,double);
