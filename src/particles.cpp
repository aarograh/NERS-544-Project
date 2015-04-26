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
  isAlive = true;

  position[0] = xyz[0];
  position[1] = xyz[1];
  position[2] = xyz[2];

  omega[0] = sqrt(1.0 - mu*mu)*cos(gamma);
  omega[1] = sqrt(1.0 - mu*mu)*sin(gamma);
  omega[2] = mu;

  energy = E_in;
  weight = 1.0;

  totalXS = 0.0;
  f235 = 0.0;
  f238 = 0.0;
  fH = 0.0;
  fcap = 0.0;
  fiss_frac = 0.0;
  abs_frac = 0.0;

  return;
}

int particle::simulate()
{
  int result, surfid;
  int isotope;
  const int fuelid = 0; const int modid = 1; // just temporary, need to get these from elsewhere
  double vn, xi;
  double dcoll, dsurf, intersection[3], tmp;
  cell* cellptr;
  surface* surfptr;
  fuel* thisFuel = new fuel(fuelid);
  moderator* thisMod = new moderator(modid);

std::cout << "x=" << position[0] << " y=" << position[1] << " z=" << position[2];
std::cout << " omegax=" << omega[0] << " omegay=" << omega[1] << " omegaz=" << omega[2] << std::endl;
  while (isAlive)
  {
    // Get pointer to the current cell
    cellptr = getPtr_cell(cellid);
    if((*cellptr).id == fuelid)
    {
        (*thisFuel).fuelMacro(energy,&totalXS,&f235,&f238);
        std::cout << "Total XS = " << totalXS << std::endl;
        std::cout << "U235 fraction = " << f235 << std::endl;
        std::cout << "U238 fraction = " << f238 << std::endl;
    }
    else if((*cellptr).id == modid)
    {
      (*thisMod).modMacro(energy,&totalXS,&fH,&fcap); 
    }
    else
    {
      std::cout << "Not fuel or moderator id." << std::endl;
      exit(-3);
    }

    // Get collision distance
    // dcoll = 500.0; // Just to test the surface/cell stuff
    dcoll = -log(drand())/(totalXS);
    // Get closest surface distance
    dsurf = (*cellptr).distToIntersect(position, omega, intersection, surfid);
std::cout << "dcoll=" << dcoll << " dsurf=" << dsurf << " surfid=" << surfid <<  std::endl;
    
    std::cout << "Total XS = " << totalXS << std::endl;
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
    // Move particle to collision point and sample collision
    else
    {
      moveParticle(dcoll);
      switch((*cellptr).id)
      { 
        case fuelid:
          std::cout << "Particle in fuel." << std::endl;
          isotope = (*thisFuel).sample_U(energy,&f235,&f238,&fiss_frac,&abs_frac); 
          std::cout << "Interaction with isotope " << isotope << std::endl;
          xi = drand();
          if(xi > abs_frac) // scatter
          {
            std::cout << "Particle scattered in fuel." << std::endl;
            vn = sqrt(2.0*energy/neut_mass)*lightspeed;
            std::cout << "Incoming particle velocity = " << vn << std::endl;
            elastic(temp,isotope,&vn,omega);
            energy = neut_mass*(vn/lightspeed)*(vn/lightspeed)/2.0; 
            std::cout << "Outgoing particle velocity = " << vn << std::endl;
            std::cout << "Outgoing particle energy = " << energy << std::endl;
          }
          else // absorption
          {
            isAlive = false;
            if(xi > fiss_frac) // capture, maybe score which isotope
            {
              std::cout << "Particle was captured in fuel." << std::endl;
              result = 0;
            }
            else // fission
            {
              std::cout << "Particle fissioned." << std::endl;
              result = static_cast<int>(nu+drand());
            }
          }
          break;
        case modid:
          if(drand() < fH) // interaction with hydrogen
          {
            if(drand() > fcap)
            {  
              std::cout << "Particle scattered off hydrogen in moderator." << std::endl;
              vn = sqrt(2.0*energy/neut_mass)*lightspeed;
              std::cout << "Incoming particle velocity = " << vn << std::endl;
              elastic(temp,1,&vn,omega);
              energy = neut_mass*(vn/lightspeed)*(vn/lightspeed)/2.0; 
              std::cout << "Outgoing particle velocity = " << vn << std::endl;
              std::cout << "Outgoing particle energy = " << energy << std::endl;
            }
            else // capture; score estimator, end history, etc. 
            {
              result = 0;
              std::cout << "Particle was captured in moderator." << std::endl;
              isAlive = false;
            }
          }
          else // interaction with oxygen; all are scatters
          {
            std::cout << "Particle scattered off oxygen in moderator." << std::endl;
            vn = sqrt(2.0*energy/neut_mass)*lightspeed;
            elastic(temp,16,&vn,omega);
            energy = neut_mass*(vn/lightspeed)*(vn/lightspeed)/2.0; 
            std::cout << "Outgoing particle velocity = " << vn << std::endl;
            std::cout << "Outgoing particle energy = " << energy << std::endl;
          }
          break;
        default:
          std::cout << "Not fuel or moderator id." << std::endl;
          exit(-3);
      }
    }
  //  for(int i = 0; i < 500000000; i++)
  //  {
  //  }
  }

  return result;
}

double particle::Coordinate(int index)
{
  return position[index];
}

double particle::Direction(int index)
{
  return omega[index];
}

void particle::moveParticle(double dist)
{
  for(int i = 0; i < 3; i++)
  {
    position[i] += dist*omega[i];
  }
  return;
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

