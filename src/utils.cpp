// AUTHORS: Aaron Graham, Mike Jarrett
// PURPOSE: NERS 544 Course Project
// DATE   : April 30, 2015

#include<cstdlib>
#include<limits>
#include<cmath>
#include "utils.h"

// Random number generator on [0,1]
double drand(void)
{
  return static_cast<double>(rand())/static_cast<double>(RAND_MAX);
}

// Samples the Watt Spectrum
double Watt(void)
{
  double a = 0.988; // MeV
  double b = 2.249; // MeV^-1

  double x1 = drand();
  double x2 = drand();
  double x3 = drand();
  double x4 = drand();
  
  double W = a*(-log(x1)-log(x2)*cos(x3*pi/2)*cos(x3*pi/2));
  return W + a*a*b/4 + (2*x4-1)*sqrt(a*a*b*W); 
}

// These operators are used for ==, >=, and <= for floating point values
bool approxeq(double x1, double x2)
{
  return (x1 > x2 - eps && x1 < x2 + eps);
}

bool approxge(double x1, double x2)
{
  return (x1 > x2 - eps);
}

bool approxle(double x1, double x2)
{
  return (x1 < x2 + eps);
}

bool softeq(double x1, double x2, double tol)
{
  return (x1 > x2 - tol && x1 < x2 + tol);
}
