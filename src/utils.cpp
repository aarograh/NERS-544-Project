// AUTHORS: Aaron Graham, Mike Jarrett
// PURPOSE: NERS 544 Course Project
// DATE   : April 3, 2015

#include<cstdlib>
#include<cmath>

// Random number generator on [0,1]
double normRand(void)
{
  return static_cast<double>(rand())/static_cast<double>(RAND_MAX);
}

double Watt(void)
{
  double a = 0.988; // MeV
  double b = 2.249; // MeV^-1
  double pi = 3.14159235658979;

  double x1 = normRand();
  double x2 = normRand();
  double x3 = normRand();
  double x4 = normRand();
  
  double W = a*(-log(x1)-log(x2)*cos(x3*pi/2)*cos(x3*pi/2));
  return W + a*a*b/4 + (2*x4-1)*sqrt(a*a*b*W); 
}
