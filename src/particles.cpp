// AUTHORS: Aaron Graham, Mike Jarrett
// PURPOSE: NERS 544 Course Project
// DATE   : April 3, 2015

#include<cstdlib>
#include<cmath>
#include<vector>
#include "geometry.h"
#include "materials.h"
#include "particles.h"
#include "utils.h"

particle::particle(double xyz[3], double gamma, double mu, double E_in,
    int cellid_in)
{
  cellid = cellid_in;
  position[0] = xyz[0];
  position[1] = xyz[1];
  position[2] = xyz[2];

  omega[0] = sqrt(1.0 - mu*mu)*cos(gamma);
  omega[1] = sqrt(1.0 - mu*mu)*sin(gamma);
  omega[2] = mu;

  energy = E_in;
  weight = 1.0;

  return;
}

int particle::simulate()
{
  int result, surfid;
  double dcoll, dsurf, intersection[3], tmp;
  cell* currentCell;
  surface* intersectSurface;

  while (isAlive)
  {
    currentCell = getPtr_cell(cellid);
    // Get collision distance
    dcoll = 500.0; // Just to test the surface/cell stuff
    // Get closest surface distance
    dsurf = (*currentCell).distToIntersect(position, omega, intersection, surfid);
    // Move particle to surface
    if (dcoll < dsurf)
    {
      intersectSurface = getPtr_surface(surfid);
      switch((*intersectSurface).boundaryType)
      {
        // Particle hit reflecting boundary
        case reflecting:
          (*intersectSurface).reflect(intersection, omega);
          break;
        // Particle hit vacuum boundary and escaped
        case vacuum:
          // Set return value and "kill" particle
          result = -surfid;
          isAlive = false;
          break;
        // Particle hit interior surface
        case interior:
          // TODO: nudge particle into surface on other side
          break;
        default:
          std::cout << "Error in particle::simulate().  Particle encountered " <<
            "unknown boundary type." << std::endl;
          exit(-2);
      }
    }
    // Move particle to collision point and perform calculations
    else
    {
    }
  }

  return 0;
}

double particle::Coordinate(int index)
{
  return position[index];
}

double particle::Direction(int index)
{
  return omega[index];
}

particle* fissionNeutron(particle* neutron)
{
  double tmp[3] = {0.0,0.0,0.0};
  particle* fissNeutron = new particle(tmp,0.0,1.0,Watt(),1);
  return fissNeutron;
}

void makeSource(std::vector<particle*> fissionBank, std::vector<particle*> sourceBank, int batch_size)
{
  double sourceProb, xi;
  // use fission bank to create source bank for next iteration
  if (fissionBank.size() == batch_size) // Fission bank is correct size
  {
    while (!fissionBank.empty())
    {
      sourceBank.push_back(fissionBank.back());
      fissionBank.pop_back();
    }
  }
  else
  {
    if (fissionBank.size() > batch_size) // Fission bank is too large
    {
      //Add neutrons to source bank with probability batch_size/fissionBank.size()
      while(!fissionBank.empty())
      {
        xi = drand(); 
        sourceProb = static_cast<double>((batch_size-sourceBank.size())/fissionBank.size());
        if(xi < sourceProb)
        {
          sourceBank.push_back(fissionBank.back());
        }
        fissionBank.pop_back();
      }
    }
    else if (fissionBank.size() < batch_size) // Fission bank is too small
    {
      //Add neutrons to source bank with probability batch_size/fissionBank.size()
      while (!fissionBank.empty())
      {
        sourceBank.push_back(fissionBank.back());
      }
      while(!fissionBank.empty())
      {
        xi = drand(); 
        sourceProb = static_cast<double>((batch_size-sourceBank.size())/fissionBank.size());
        if(xi < sourceProb)
        {
          sourceBank.push_back(fissionBank.back());
        }
        fissionBank.pop_back();
      }
    }
  }
  return;
}

double calcEntropy(std::vector<particle*> fissionBank)
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
  double x, y;
  int z_index, r_index;
  for(int i = 0; i < fissionBank.size(); i++){
    x = (*fissionBank.at(i)).Coordinate(0);
    y = (*fissionBank.at(i)).Coordinate(1);
    p_rad = x*x + y*y; 
    r_index = (int)(p_rad/area);
    z_index = (int)(*fissionBank.at(i)).Coordinate(2)/dz;
    particle_mesh[nrad*z_index + nrad] = particle_mesh[nrad*z_index + nrad] + 1;
  }
  // calculate Shannon entropy
  double entropy = 0.0;
  for(int i = 0; i < nz; i++)
  {
    for(int j = 0; j < nrad; j++)
    {
      pn = (double)(particle_mesh[nrad*i + j])/(double)(fissionBank.size());
      entropy = entropy + pn*log2(pn); 
    }
  }
  entropy = -entropy;
  return entropy;
}

