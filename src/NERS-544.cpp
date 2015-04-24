// AUTHORS: Aaron Graham, Mike Jarrett
// PURPOSE: NERS 544 Course Project
// DATE   : April 3, 2015

#include<iostream>
#include<iomanip>
#include<cstdlib>
#include<cmath>
#include "utils.h"
#include "particles.h"
#include "geometry.h"
#include "materials.h"

using namespace std;

int main()
{
  srand(time(NULL));

  int fuelid, modid;
  double pitch;
  cout << "Enter the pin pitch in cm (must be greater than 3.0):";
  cin >> pitch;

  if (pitch <= 3.0)
  {
    cout << "Error!  Pin pitch must be greater than pin diameter of 3.0 cm!"
      << endl;
    return -1;
  }

  init_materials(fuelid, modid);
  initPinCell(pitch, fuelid, modid);

  int batch_size = 1E5;
  double En;
  double xyz[3];
  double pinrad = 1.5; // pin radius = 1.5 cm
  double r, gamma, mu;
  double pi = 3.14159265358979;
  vector<particle*> sourceBank;

  // sample neutrons for initial source bank
  for(int i = 0; i < batch_size; i++){
    // sample energy from Watt spectrum
    En = Watt();
    // sample a radial location within the fuel cell
    gamma = 2*pi*drand();
    r = pinrad*sqrt(drand());
    xyz[0] = r*cos(gamma);
    xyz[1] = r*sin(gamma);
    // sample an axial location
    xyz[2] = 100.0*drand(); 
    // sample a direction
    gamma = drand();
    mu = 2.0*drand() - 1.0;
    // add particle to the bank
    sourceBank.push_back(new particle(xyz,gamma,mu,En,fuelid));
  }
  
  // outer loop over power iterations
  bool converged = false;
  double ShannonEntropy = 0.0, PrevEntropy = 0.0;
  double tol_entropy = 1E-3;
  const int max_iters = 100, active_iters = 10;
  int k = 0, l = 0, n = 0, fissions;
  particle* neutron;
  vector<particle*> fissionBank;

  while(k < max_iters){
    k = k+1; // total power iterations 
    
    // inner loop over the source bank
    for (n = batch_size - 1; n > 0; n--)
    {
      neutron = sourceBank.at(n);
      fissions = (*neutron).simulate();
      for (int i = 0; i < fissions; i++)
      {
        fissionBank.push_back(fissionNeutron(neutron));
      }
    }

    // use fission bank to create source bank for next iteration

    // Calculate Shannon Entropy
    PrevEntropy = ShannonEntropy;
    //ShannonEntropy = calcEntropy(batch_size,xyz);
    // some convergence check
    if(converged){
      l = l+1; // power iterations with converged source
      if (l >= active_iters) break;
    }
    else{
      // run some convergence check
      if((ShannonEntropy-PrevEntropy) < tol_entropy){
        converged = true;
      }
    }
  }

  return 0;
}
