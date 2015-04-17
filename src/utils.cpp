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

double calcEntropy(int batch_size, double xyz[][3])
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
  double p_rad, pn;
  int z_index, r_index;
  for(int i = 0; i < batch_size; i++){
    p_rad = xyz[i][0]*xyz[i][0] + xyz[i][1]*xyz[i][1]; 
    r_index = (int)(p_rad/area);
    z_index = (int)(xyz[i][2]/dz);
    particle_mesh[nrad*z_index + nrad] = particle_mesh[nrad*z_index + nrad] + 1;
  }
  // calculate Shannon entropy
  double entropy = 0.0;
  for(int i = 0; i < nz; i++)
  {
    for(int j = 0; j < nrad; j++)
    {
      pn = (double)(particle_mesh[nrad*i + j])/(double)(batch_size);
      entropy = entropy + pn*log2(pn); 
    }
  }
  entropy = -entropy;
  return entropy;
}
