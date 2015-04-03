// AUTHORS: Aaron Graham, Mike Jarrett
// PURPOSE: NERS 544 Course Project
// DATE   : April 3, 2015

#include<cstdlib>

// Random number generator on [0,1]
double normRand(void)
{
  return static_cast<double>(rand())/static_cast<double>(RAND_MAX);
}
