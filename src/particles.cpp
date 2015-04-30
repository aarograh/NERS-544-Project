// AUTHORS: Aaron Graham, Mike Jarrett
// PURPOSE: NERS 544 Course Project
// DATE   : April 30, 2015

#include<cstdlib>
#include<cmath>
#include<vector>
#include "geometry.h"
#include "materials.h"
#include "particles.h"
#include "utils.h"

double fuelSpectrum[1001];
double modSpectrum[1001];
double energyGrid[1001];

particle::particle(const double pos_in[3], double gamma, double mu,
  double E_in, int cellid_in)
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
  cutoff = 0.1;
  survival = 0.4;

  totalXS = 0.0;
  f235 = 0.0;
  f238 = 0.0;
  fH = 0.0;
  fcap = 0.0;
  fiss_frac = 0.0;
  abs_frac = 0.0;

  score = 0.0;
  estimatorTL = 0.0;
  //squareTL = 0.0;
  estimatorColl = 0.0;
  //squareColl = 0.0;
}

fission::fission(const particle& neutron, int cellid_in)
{
  position[0] = neutron.position[0]; 
  position[1] = neutron.position[1]; 
  position[2] = neutron.position[2]; 

  cellid = cellid_in;
}

int particle::simulate()
{
  int result, surfid;
  int isotope;
  const int fuelid = 0; const int modid = 1; 
  double vn, xi;
  double dcoll, dsurf, intersection[3];
  cell* cellptr;
  surface* surfptr;
  fuel* thisFuel = new fuel(fuelid);
  moderator* thisMod = new moderator(modid);

  while (isAlive)
  {
    // Get pointer to the current cell
    cellptr = getPtr_cell(cellid);
    if(cellptr->id == fuelid)
    {
      thisFuel->fuelMacro(energy,&totalXS,&f235,&f238,&fiss_frac,&abs_frac);
    }
    else if(cellptr->id == modid)
    {
      thisMod->modMacro(energy,&totalXS,&fH,&fcap); 
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
    
    // Move particle to surface
    if (dsurf < dcoll)
    {
      // tally the track length estimator for keff
      if(cellptr->id == fuelid)
      {
        score = dsurf*weight*nu*fiss_frac*totalXS;
        estimatorTL = estimatorTL + score;
        //squareTL = squareTL + score*score;
        
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
          break;
        // Particle hit vacuum boundary and escaped
        case vacuum:
          // Set return value and "kill" particle
          result = -surfid;
          isAlive = false;
          break;
        // Particle hit interior surface
        case interior:
          position[0] += omega[0]*nudge;
          position[1] += omega[1]*nudge;
          position[2] += omega[2]*nudge;
          
          cellid = getCellID(position);
          cellptr = getPtr_cell(cellid);
          break;
        default:
          std::cout << "Error in particle::simulate().  Particle encountered " 
            << "unknown boundary type." << std::endl;
          exit(-2);
      }
    }
    // Move particle to collision point and sample collision
    else
    {
      position[0] += omega[0]*dcoll;
      position[1] += omega[1]*dcoll;
      position[2] += omega[2]*dcoll;
      switch(cellptr->id)
      { 
        case fuelid:
          isotope = thisFuel->sample_U(energy,&f235,&f238,&abs_frac,
            &fiss_frac);
          // tally the track length estimator and the collision estimator
          score = dcoll*weight*nu*fiss_frac*totalXS;
          //std::cout << "Particle weight = " << weight << std::endl;
	  //std::cout << "absorption fraction = " << abs_frac << std::endl;
	  //std::cout << "fission fraction = " << fiss_frac << std::endl;
          //std::cout << "total XS = " << totalXS << std::endl;
          estimatorTL = estimatorTL + score;
          //squareTL = squareTL + score*score;
          score = weight*nu*fiss_frac;
          //std::cout << "Score added = " << score << std::endl;
          estimatorColl = estimatorColl + score;
          //squareColl = squareColl + score*score;

          xi = drand();
          if(xi > abs_frac) // scatter
          {
            isotope = thisFuel->sample_U(&f235,&f238);
            vn = sqrt(2.0*energy/neut_mass)*lightspeed;
            elastic(temp,isotope,vn,omega);
            energy = neut_mass*(vn/lightspeed)*(vn/lightspeed)/2.0; 
          }
          else // absorption
          {
            isAlive = false;
            if(xi > fiss_frac) // capture, maybe score which isotope
            {
              result = 0;
            }
            else // fission
            {
              result = static_cast<int>(weight*nu+drand());
            }
          }
          break;
        case modid:
          if(drand() < fH) // interaction with hydrogen
          {
            if(drand() > fcap)
            {  
              vn = sqrt(2.0*energy/neut_mass)*lightspeed;
              elastic(temp,1,vn,omega);
              energy = neut_mass*(vn/lightspeed)*(vn/lightspeed)/2.0; 
            }
            else // capture; score estimator, end history, etc. 
            {
              result = 0;
              isAlive = false;
            }
          }
          else // interaction with oxygen; all are scatters
          {
            vn = sqrt(2.0*energy/neut_mass)*lightspeed;
            elastic(temp,16,vn,omega);
            energy = neut_mass*(vn/lightspeed)*(vn/lightspeed)/2.0; 
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

int particle::simulate_implicit()
{
  int result, surfid;
  int isotope;
  const int fuelid = 0; const int modid = 1; 
  double vn, xi;
  double dcoll, dsurf, intersection[3];
  cell* cellptr;
  surface* surfptr;
  fuel* thisFuel = new fuel(fuelid);
  moderator* thisMod = new moderator(modid);

  while (isAlive)
  {
    // Get pointer to the current cell
    cellptr = getPtr_cell(cellid);
    if(cellptr->id == fuelid)
    {
      thisFuel->fuelMacro(energy,&totalXS,&f235,&f238,&fiss_frac,&abs_frac);
    }
    else if(cellptr->id == modid)
    {
      thisMod->modMacro(energy,&totalXS,&fH,&fcap); 
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
    
    // Move particle to surface
    if (dsurf < dcoll)
    {
      spectrumTally(energy,dsurf*weight,cellptr->id);
      // tally the track length estimator for keff
      if(cellptr->id == fuelid)
      {
        score = dsurf*weight*nu*fiss_frac*totalXS;
        estimatorTL = estimatorTL + score;
        //squareTL = squareTL + score*score;
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
          break;
        // Particle hit vacuum boundary and escaped
        case vacuum:
          // Set return value and "kill" particle
          result = -surfid;
          isAlive = false;
          break;
        // Particle hit interior surface
        case interior:
          position[0] += omega[0]*nudge;
          position[1] += omega[1]*nudge;
          position[2] += omega[2]*nudge;
          
          cellid = getCellID(position);
          cellptr = getPtr_cell(cellid);
          break;
        default:
          std::cout << "Error in particle::simulate().  Particle encountered " 
            << "unknown boundary type." << std::endl;
          exit(-2);
      }
    }
    // Move particle to collision point and sample collision
    else
    {
      position[0] += omega[0]*dcoll;
      position[1] += omega[1]*dcoll;
      position[2] += omega[2]*dcoll;
      if(cellptr->id == fuelid)
      { 
        spectrumTally(energy,dcoll*weight,fuelid);
        isotope = thisFuel->sample_U(&f235,&f238);
        // tally the track length estimator and the collision estimator
        score = dcoll*weight*nu*fiss_frac*totalXS;
        estimatorTL = estimatorTL + score;
        //squareTL = squareTL + score*score;
        score = weight*nu*fiss_frac;
        estimatorColl = estimatorColl + score;
        //squareColl = squareColl + score*score;

        xi = drand();

        vn = sqrt(2.0*energy/neut_mass)*lightspeed;
        elastic(temp,isotope,vn,omega);
        energy = neut_mass*(vn/lightspeed)*(vn/lightspeed)/2.0; 
        if(xi < fiss_frac)
        {
//          std::cout << "Particle fissioned." << std::endl;
          result = static_cast<int>(weight*nu+drand());
          isAlive = false;
        }
        weight = weight*(1.0-abs_frac);
//        std::cout << "Implicit capture: absorption fraction = " << abs_frac << std::endl;
        if(weight < cutoff)
        { 
          isAlive = roulette();
        }
      }
      else if(cellptr->id == modid)
      {
        spectrumTally(energy,dcoll*weight,modid);
        if(drand() < fH) // interaction with hydrogen
        {
          isotope = 1;
        }
        else // interaction with oxygen; all are scatters
        {
          isotope = 16;
        }
        vn = sqrt(2.0*energy/neut_mass)*lightspeed;
        elastic(temp,isotope,vn,omega);
        energy = neut_mass*(vn/lightspeed)*(vn/lightspeed)/2.0; 

        weight = weight*(1.0-fH*fcap);
//        std::cout << "Implicit capture fraction = " << fH*fcap << std::endl;
        if(weight < cutoff)
        { 
          isAlive = roulette();
        }
      }
      else
      {
        std::cout << "Not fuel or moderator id." << std::endl;
        exit(-3);
      }
    }
  }

  delete thisFuel;
  delete thisMod;
  return result;
}

bool particle::roulette()
{
  if(drand() < weight/survival)
  { 
    weight = survival;
    return true;
  }
  else
  {
//    std::cout << "The particle was killed in rouletting with a weight of "
//              << weight << std::endl;
    weight = 0.0;
    return false;
  }
}
void makeSource(std::vector<fission> &fissionBank, 
    std::vector<particle> &sourceBank, int batch_size)
{
  double sourceProb, xi;
  fission* fissptr;

  if (fissionBank.size() == batch_size)
  {
    while(!fissionBank.empty())
    {
      sourceBank.push_back(particle((fissionBank.back()).position,2.0*pi*drand(),2.0*drand()-
        1.0,Watt(),0)); 
      fissionBank.pop_back();
    }
  }
  else if (fissionBank.size() > batch_size) // Fission bank is too large
  {
    // Add to source bank with probability batch_size/fissionBank.size()
    while(!fissionBank.empty())
    {
      xi = drand(); 
      sourceProb = static_cast<double>((batch_size-sourceBank.size())/
                   static_cast<double>(fissionBank.size()));
      if(xi < sourceProb)
      {
        fissptr = &fissionBank.back();
        sourceBank.push_back(particle(fissptr->position,2.0*pi*drand(),2.0*drand()-
          1.0,Watt(),0));
      }
      fissionBank.pop_back();
    }
  }
  else if (fissionBank.size() < batch_size) // Fission bank is too small
  {
    // Add to source bank with probability batch_size/fissionBank.size()
    for(int j = 0; j < (int)(batch_size/fissionBank.size()); j++)
    {
      for(int k = 0; k < fissionBank.size(); k++)
      {
        fissptr = &fissionBank[k];
        sourceBank.push_back(particle(fissptr->position,2.0*pi*drand(),2.0*drand()-
          1.0,Watt(),0));
      }
    }
    while(!fissionBank.empty())
    {
      xi = drand(); 
      sourceProb = static_cast<double>(batch_size-sourceBank.size())/
                   static_cast<double>(fissionBank.size());
      if(xi < sourceProb)
      {
        fissptr = &fissionBank.back();
        sourceBank.push_back(particle(fissptr->position,2.0*pi*drand(),2.0*drand()-
          1.0,Watt(),0));
      }
      fissionBank.pop_back();
    }
  }
  return;
}

double calcEntropy(std::vector<fission> fissionBank)
{
  double radius = 1.5;
  int nrad = 10;
  int nz = 100;
  int nbins = nz*nrad;
  int particle_mesh[nbins];
  double dz = 100.0/nz;
  double area = radius*radius/nrad; 
  for(int i = 0; i < nbins; i++)
  {
    particle_mesh[i] = 0;
  }
  
  int index;

  // bin all of the particles
  // uniform axial bins, equal-area radial bins
  double p_rad, pn;
  double x, y;
  int z_index, r_index;
  for(int i = 0; i < fissionBank.size(); i++){
    x = fissionBank[i].position[0];
    y = fissionBank[i].position[1];
    p_rad = x*x + y*y;
    r_index = (int)(p_rad/area);
    if(r_index >= nrad)
    {
      std::cout << "Outside of the fuel region" << std::endl;
      std::cout << "x = " << fissionBank[i].position[0] << std::endl;
      std::cout << "y = " << fissionBank[i].position[1] << std::endl;
      std::cout << "z = " << fissionBank[i].position[2] << std::endl;
    }
    z_index = (int)fissionBank[i].position[2]/dz;
    particle_mesh[nrad*z_index + r_index] = 
      particle_mesh[nrad*z_index + r_index] + 1;
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

void particle::spectrumTally(double energy, double fluxTally, int id)
{
  int g = 0;
  while(energy > energyGrid[g])
  {
    g = g+1;
  }

  if(id == 0)
  {
    fuelSpectrum[g] = fuelSpectrum[g] + fluxTally; 
  }
  else if(id == 1)
  {
    modSpectrum[g] = modSpectrum[g] + fluxTally; 
  }
}
