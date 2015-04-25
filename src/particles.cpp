// AUTHORS: Aaron Graham, Mike Jarrett
// PURPOSE: NERS 544 Course Project
// DATE   : April 3, 2015

#include "utils.h"
#include "geometry.h"
#include "materials.h"
#include "particles.h"

particle::particle(double xyz[3], double gamma, double mu, double E_in,
    int cellid_in)
{
  cellid = cellid_in;
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

particle* fissionNeutron(particle* neutron)
{
  double tmp[3] = {0.0,0.0,0.0};
  particle* fissNeutron = new particle(tmp,0.0,0.0,0.0,1);
  return fissNeutron;
}
