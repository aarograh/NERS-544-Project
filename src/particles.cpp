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
  isAlive = true;

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
  cell* cellptr;
  surface* surfptr;

std::cout << "x=" << position[0] << " y=" << position[1] << " z=" << position[2];
std::cout << " omegax=" << omega[0] << " omegay=" << omega[1] << " omegaz=" << omega[2] << std::endl;
  while (isAlive)
  {
    // Get pointer to the current cell
    cellptr = getPtr_cell(cellid);
    // Get collision distance
    dcoll = 500.0; // Just to test the surface/cell stuff
    // Get closest surface distance
    dsurf = (*cellptr).distToIntersect(position, omega, intersection, surfid);
std::cout << "dcoll=" << dcoll << " dsurf=" << dsurf << " surfid=" << surfid <<  std::endl;
    // Move particle to surface
    if (dsurf < dcoll)
    {
      // get pointer to the surface that the particle is colliding with
      surfptr = getPtr_surface(surfid);
      switch((*surfptr).boundaryType)
      {
        // Particle hit reflecting boundary
        case reflecting:
          position[0] += omega[0]*dsurf;
          position[1] += omega[1]*dsurf;
          position[2] += omega[2]*dsurf;
          (*surfptr).reflect(intersection, omega);
          position[0] += omega[0]*nudge;
          position[1] += omega[1]*nudge;
          position[2] += omega[2]*nudge;
std::cout << "x=" << position[0] << " y=" << position[1] << " z=" << position[2];
std::cout << " omegax=" << omega[0] << " omegay=" << omega[1] << " omegaz=" << omega[2] << std::endl;
          break;
        // Particle hit vacuum boundary and escaped
        case vacuum:
          // Set return value and "kill" particle
          result = -surfid;
          isAlive = false;
          break;
        // Particle hit interior surface
        case interior:
          position[0] += omega[0]*(dsurf + nudge);
          position[1] += omega[1]*(dsurf + nudge);
          position[2] += omega[2]*(dsurf + nudge);
std::cout << "x=" << position[0] << " y=" << position[1] << " z=" << position[2];
std::cout << " omegax=" << omega[0] << " omegay=" << omega[1] << " omegaz=" << omega[2] << std::endl;
          
          cellid = getCellID(position);
          cellptr = getPtr_cell(cellid);
          break;
        default:
std::cout << surfid << " " << (*surfptr).boundaryType << std::endl;
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

  return result;
}

particle* fissionNeutron(particle* neutron)
{
  double tmp[3] = {0.0,0.0,0.0};
  particle* fissNeutron = new particle(tmp,0.0,0.0,0.0,1);
  return fissNeutron;
}
