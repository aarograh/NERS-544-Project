// AUTHORS: Aaron Graham, Mike Jarrett
// PURPOSE: NERS 544 Course Project
// DATE   : April 3, 2015

#include "materials.h"
#include "geometry.h"
#include "particles.h"

double particles::calcEntropy(int batch_size, double xloc[], double yloc[], double zloc[])
{
  double radius = 1.5;
  int nrad = 10;
  int nz = 10;
  int nbins = nz*nrad;
  int particle_mesh[nbins];
  double dz = 100.0/nz;
  double area = radius*radius/nrad; 
  
  // bin all of the particles
  // uniform axial bins, equal-area radial bins
  double p_rad;
  int z_index, r_index;
  for(int i = 0; i < batch_size; i++){
    p_rad = xloc[i]*xloc[i] + yloc[i]*yloc[i]; 
    r_index = (int)(p_rad/area);
    z_index = (int)(zloc[i]/dz);
    particle_mesh[nrad*z_index + nrad] = particle_mesh[nrad*z_index + nrad] + 1;
  }
  double entropy = 0.0;
  for(int i = 0; i < nz; i++)
  {
    for(int j = 0; j < nrad; j++)
    {
      pi = (double)(particle_mesh[nrad*i + j])/(double)(batch_size);
      entropy = entropy + pi*log2(pi); 
    }
  }
  entropy = -entropy;
  return entropy;
  // calculate Shannon entropy
}
