// AUTHORS: Aaron Graham, Mike Jarrett
// PURPOSE: NERS 544 Course Project
// DATE   : April 3, 2015

#include<cstdlib>
#include<limits>
#include<cmath>
#include "utils.h"
#include "particles.h"

const double eps = std::numeric_limits<double>::epsilon()*100.0;

// Random number generator on [0,1]
double drand(void)
{
  return static_cast<double>(rand())/static_cast<double>(RAND_MAX);
}

double Watt(void)
{
  double a = 0.988; // MeV
  double b = 2.249; // MeV^-1
  double pi = 3.14159235658979;

  double x1 = drand();
  double x2 = drand();
  double x3 = drand();
  double x4 = drand();
  
  double W = a*(-log(x1)-log(x2)*cos(x3*pi/2)*cos(x3*pi/2));
  return W + a*a*b/4 + (2*x4-1)*sqrt(a*a*b*W); 
}

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
