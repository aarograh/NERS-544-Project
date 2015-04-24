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
  initPinCell(pitch);

  return 0;

  int batch_size = 1E5;
  double En[batch_size];
  double xyz[batch_size][3];
  double pinrad = 1.5; // pin radius = 1.5 cm
  double r, gamma;
  double pi = 3.14159265358979;

  // sample neutrons for initial source bank
  for(int i = 0; i < batch_size; i++){
    // sample energy from Watt spectrum
    En[i] = Watt();
    // sample a radial location within the fuel cell
    gamma = 2*pi*normRand();
    r = pinrad*sqrt(normRand());
    xyz[i][0] = r*cos(gamma);
    xyz[i][1] = r*sin(gamma);
    // sample an axial location
    xyz[i][2] = 100.0*normRand(); 
  }
  
  // outer loop over power iterations
  bool converged = false;
  double ShannonEntropy = 0.0, PrevEntropy = 0.0;
  double tol_entropy = 1E-3;
  int max_iters = 100;
  int k = 0, l = 0, n = 0;
  while(k < max_iters){
    k = k+1; // total power iterations 
    
    // inner loop over the source bank
    while(n < batch_size){
      n = n+1;
    }

    // use fission bank to create source bank for next iteration

    PrevEntropy = ShannonEntropy;
    ShannonEntropy = calcEntropy(batch_size,xyz);
    // some convergence check
    if(converged){
      l = l+1; // power iterations with converged source
    }
    else{
      // run some convergence check
      if((ShannonEntropy-PrevEntropy) < tol_entropy){
        converged = true;
      }
    }
  }
}
