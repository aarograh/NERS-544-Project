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

particle::particle(double pos_in[3], double gamma, double mu, double E_in,
    int cellid_in)
{
  cellid = cellid_in;
  isAlive = true;

  position[0] = pos_in[0]; 
  position[1] = pos_in[1]; 
  position[2] = pos_in[2]; 

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

  score = 0.0;
  estimatorTL = 0.0;
  squareTL = 0.0;
  estimatorColl = 0.0;
  squareColl = 0.0;
}

fission::fission()
{
  position[0] = 0.0; 
  position[1] = 0.0; 
  position[2] = 0.0; 
  
  cellid = 0;
}

fission::fission(double pos_in[3], int cellid_in)
{
  position[0] = pos_in[0]; 
  position[1] = pos_in[1]; 
  position[2] = pos_in[2]; 

  cellid = cellid_in;
}

int particle::simulate()
{
  int result, surfid;
  int isotope;
  const int fuelid = 0; const int modid = 1; // just temporary, need to get these from elsewhere
  double vn, xi;
  double dcoll, dsurf, intersection[3];
  cell* cellptr;
  surface* surfptr;
  fuel* thisFuel = new fuel(fuelid);
  moderator* thisMod = new moderator(modid);

//std::cout << "x=" << position[0] << " y=" << position[1] << " z=" << position[2];
//std::cout << " omegax=" << omega[0] << " omegay=" << omega[1] << " omegaz=" << omega[2] << std::endl;
  while (isAlive)
  {
    // Get pointer to the current cell
//std::cout << "Currently in cell " << cellid << std::endl;
    cellptr = getPtr_cell(cellid);
    if(cellptr->id == fuelid)
    {
      thisFuel->fuelMacro(energy,&totalXS,&f235,&f238);
//std::cout << "Total XS = " << totalXS << std::endl;
//std::cout << "U235 fraction = " << f235 << std::endl;
//std::cout << "U238 fraction = " << f238 << std::endl;
    }
    else if(cellptr->id == modid)
    {
      thisMod->modMacro(energy,&totalXS,&fH,&fcap); 
//std::cout << "Total XS = " << totalXS << std::endl;
//std::cout << "Hydrogen fraction = " << fH << std::endl;
//std::cout << "Capture fraction = " << fcap << std::endl;
    }
    else
    {
      std::cout << "Not fuel or moderator id." << std::endl;
      exit(-3);
    }
    // get distance to next collision
    //dcoll = 500.0; // Just to test the surface/cell stuff
    dcoll = -log(drand())/(totalXS);
    // Get closest surface distance
    dsurf = cellptr->distToIntersect(position, omega, intersection, surfid);
//std::cout << std::endl;
//std::cout << "dcoll=" << dcoll << " dsurf=" << dsurf << " surfid=" << surfid <<  std::endl;
    
//    std::cout << "Total XS = " << totalXS << std::endl;
//std::cout << "intersection=(" << intersection[0] << "," << intersection[1] << "," << intersection[2] << ")\n";
    // Move particle to surface
    if (dsurf < dcoll)
    {
      // tally the track length estimator for keff
      if(cellptr->id == fuelid)
      {
        score = dsurf*weight*nu*(thisFuel->fissXS(energy));
        estimatorTL = estimatorTL + score;
        squareTL = squareTL + score*score;
        
      }
      // Move particle
      position[0] = intersection[0];
      position[1] = intersection[1];
      position[2] = intersection[2];
      // get pointer to the surface that the particle is colliding with
      surfptr = getPtr_surface(surfid);
      switch(surfptr->boundaryType)
      {
        // Particle hit reflecting boundary
        case reflecting:
          surfptr->reflect(intersection, omega);
          // Nudge the particle a bit to avoid floating point issues
          position[0] += omega[0]*nudge;
          position[1] += omega[1]*nudge;
          position[2] += omega[2]*nudge;
//std::cout << "x=" << position[0] << " y=" << position[1] << " z=" << position[2];
//std::cout << " omegax=" << omega[0] << " omegay=" << omega[1] << " omegaz=" << omega[2] << std::endl;
          break;
        // Particle hit vacuum boundary and escaped
        case vacuum:
          // Set return value and "kill" particle
          result = -surfid;
          isAlive = false;
//std::cout << "x=" << position[0] << " y=" << position[1] << " z=" << position[2];
//std::cout << " omegax=" << omega[0] << " omegay=" << omega[1] << " omegaz=" << omega[2] << std::endl;
          break;
        // Particle hit interior surface
        case interior:
          position[0] += omega[0]*nudge;
          position[1] += omega[1]*nudge;
          position[2] += omega[2]*nudge;
//std::cout << "x=" << position[0] << " y=" << position[1] << " z=" << position[2];
//std::cout << " omegax=" << omega[0] << " omegay=" << omega[1] << " omegaz=" << omega[2] << std::endl;
          
          cellid = getCellID(position);
//std::cout << "cellid returned from getCellID = " << cellid << std::endl;
          cellptr = getPtr_cell(cellid);
//std::cout << "cellid = " << cellid << std::endl;
          break;
        default:
//std::cout << surfid << " " << (*surfptr).boundaryType << std::endl;
          std::cout << "Error in particle::simulate().  Particle encountered " <<
            "unknown boundary type." << std::endl;
          exit(-2);
      }
    }
    // Move particle to collision point and sample collision
    else
    {
      moveParticle(dcoll);
      switch(cellptr->id)
      { 
        case fuelid:
//std::cout << "Particle in fuel." << std::endl;
          isotope = thisFuel->sample_U(energy,&f235,&f238,&abs_frac,&fiss_frac);
          // tally the track length estimator and the collision estimator
          if(cellptr->id == fuelid)
          {
            score = dcoll*weight*nu*fiss_frac*totalXS;
            estimatorTL = estimatorTL + score;
            squareTL = squareTL + score*score;
            score = weight*nu*fiss_frac;
            estimatorColl = estimatorColl + score;
            squareColl = squareColl + score*score;
          }
//std::cout << "Interaction with isotope " << isotope << std::endl;
//std::cout << "Absorption fraction = " << abs_frac << std::endl;
//std::cout << "Fission fraction = " << fiss_frac << std::endl;
          xi = drand();
          if(xi > abs_frac) // scatter
          {
//std::cout << "Particle scattered in fuel." << std::endl;
            vn = sqrt(2.0*energy/neut_mass)*lightspeed;
//std::cout << "Incoming particle velocity = " << vn << std::endl;
            elastic(temp,isotope,&vn,omega);
            energy = neut_mass*(vn/lightspeed)*(vn/lightspeed)/2.0; 
//std::cout << "Outgoing particle velocity = " << vn << std::endl;
//std::cout << "Outgoing particle energy = " << energy << std::endl;
          }
          else // absorption
          {
            isAlive = false;
            if(xi > fiss_frac) // capture, maybe score which isotope
            {
//std::cout << "Particle was captured in fuel." << std::endl;
              result = 0;
            }
            else // fission
            {
//std::cout << "Particle fissioned." << std::endl;
              result = static_cast<int>(nu+drand());
            }
          }
          break;
        case modid:
//std::cout << "Hydrgon XS fraction = " << fH << std::endl;
          if(drand() < fH) // interaction with hydrogen
          {
            if(drand() > fcap)
            {  
//std::cout << "Particle scattered off hydrogen in moderator." << std::endl;
              vn = sqrt(2.0*energy/neut_mass)*lightspeed;
//std::cout << "Incoming particle velocity = " << vn << std::endl;
              elastic(temp,1,&vn,omega);
              energy = neut_mass*(vn/lightspeed)*(vn/lightspeed)/2.0; 
//std::cout << "Outgoing particle velocity = " << vn << std::endl;
//std::cout << "Outgoing particle energy = " << energy << std::endl;
            }
            else // capture; score estimator, end history, etc. 
            {
              result = 0;
//std::cout << "Particle was captured in moderator." << std::endl;
              isAlive = false;
            }
          }
          else // interaction with oxygen; all are scatters
          {
//std::cout << "Particle scattered off oxygen in moderator." << std::endl;
            vn = sqrt(2.0*energy/neut_mass)*lightspeed;
            elastic(temp,16,&vn,omega);
            energy = neut_mass*(vn/lightspeed)*(vn/lightspeed)/2.0; 
//std::cout << "Outgoing particle velocity = " << vn << std::endl;
//std::cout << "Outgoing particle energy = " << energy << std::endl;
          }
          break;
        default:
          std::cout << "Not fuel or moderator id." << std::endl;
          exit(-3);
      }
    }
  }

  delete thisFuel;
  delete thisMod;
  return result;
}

int particle::getID()
{
  return cellid;
}

double particle::getCoord(int index)
{
  return position[index];
}

double fission::getCoord(int index)
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
fission fissionNeutron(particle neutron)
{
  double tmp[3];
  for(int i = 0; i < 3; i++)
  {
    tmp[i] = neutron.getCoord(i);
  }
//  fission fissNeutron = fission(tmp);
//  return fissNeutron;
  return fission(tmp,0);
}

void makeSource(std::vector<fission> &fissionBank, std::vector<particle> &sourceBank, int batch_size)
{
  double sourceProb, xi;
  double xyz[3];
  if (fissionBank.size() > batch_size) // Fission bank is too large
  {
    //Add neutrons to source bank with probability batch_size/fissionBank.size()
    while(!fissionBank.empty())
    {
      xi = drand(); 
      sourceProb = static_cast<double>((batch_size-sourceBank.size())/fissionBank.size());
      if(xi < sourceProb)
      {
        for(int i = 0; i < 3; i++)
        {
          xyz[i] = (fissionBank.back()).getCoord(i);
        }
        sourceBank.push_back(particle(xyz,2*pi*drand(),drand(),Watt(),0));
      }
      fissionBank.pop_back();
    }
  }
  else if (fissionBank.size() < batch_size) // Fission bank is too small
  {
    //Add neutrons to source bank with probability batch_size/fissionBank.size()
    for(int j = 0; j < (int)(batch_size/fissionBank.size()); j++)
    {
      for(int k = 0; k < fissionBank.size(); k++)
      {
        for(int i = 0; i < 3; i++)
        {
          xyz[i] = (fissionBank.back()).getCoord(i);
        }
        sourceBank.push_back(particle(xyz,2*pi*drand(),drand(),Watt(),0));
//std::cout << "Cell ID = " << (*sourceBank.back()).getID() << std::endl;
      }
    }
    while(!fissionBank.empty())
    {
      xi = drand(); 
      sourceProb = static_cast<double>(batch_size-sourceBank.size())/
                   static_cast<double>(fissionBank.size());
      if(xi < sourceProb)
      {
        for(int i = 0; i < 3; i++)
        {
          xyz[i] = (fissionBank.back()).getCoord(i);
        }
        sourceBank.push_back(particle(xyz,2*pi*drand(),drand(),Watt(),0));
//std::cout << "Cell ID = " << (*sourceBank.back()).getID() << std::endl;
      }
      fissionBank.pop_back();
    }
  }
  return;
}

double particle::getTL(void)
{
  return estimatorTL;
}

double particle::getColl(void)
{
  return estimatorColl;
}

double particle::getTLsq(void)
{
  return squareTL;
}

double particle::getCollsq(void)
{
  return squareColl;
}

double particle::getWeight(void)
{
  return weight;
}

double calcEntropy(std::vector<fission> fissionBank)
{
  double radius = 1.5;
  int nrad = 10;
  int nz = 20;
  int nbins = nz*nrad;
  int particle_mesh[nbins];
  double dz = 100.0/nz;
  double area = radius*radius/nrad; 
  for(int i = 0; i < nbins; i++)
  {
    particle_mesh[i] = 0;
  }
  
  // bin all of the particles
  // uniform axial bins, equal-area radial bins
  double p_rad, pn;
  double x, y;
//  double x, y;
  int z_index, r_index;
  for(int i = 0; i < fissionBank.size(); i++){
    x = fissionBank[i].getCoord(0);
    y = fissionBank[i].getCoord(1);
    p_rad = x*x + y*y;
    r_index = (int)(p_rad/area);
    z_index = (int)fissionBank[i].getCoord(2)/dz;
    particle_mesh[nrad*z_index + nrad] = particle_mesh[nrad*z_index + nrad] + 1;
  }
  // calculate Shannon entropy
  double entropy = 0.0;
  for(int i = 0; i < nz; i++)
  {
    for(int j = 0; j < nrad; j++)
    {
      pn = (double)(particle_mesh[nrad*i + j])/(double)(fissionBank.size());
      if(pn > 0.0)
      {
        entropy = entropy + pn*log2(pn); 
      }
    }
  }
  entropy = -entropy;
  return entropy;
}

