// AUTHORS: Aaron Graham, Mike Jarrett
// PURPOSE: NERS 544 Course Project
// DATE   : April 3, 2015

#include<cmath>
#include "geometry.h"
#include "materials.h"
#include "particles.h"

particle::particle(double xyz[3], double gamma, double mu, double E_in,
    int matid_in)
{
  matid = matid_in;
  position[0] = xyz[0];
  position[1] = xyz[1];
  position[2] = xyz[1];

  omega[0] = sqrt(1.0 - mu*mu)*cos(gamma);
  omega[1] = sqrt(1.0 - mu*mu)*sin(gamma);
  omega[2] = mu;

  energy = E_in;
  weight = 1.0;

  return;
}
