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
  double pitch;
  cout << "Enter the pin pitch in cm (must be greater than 3.0):";
  cin >> pitch;

  if (pitch <= 3.0)
  {
    cout << "Error!  Pin pitch must be greater than pin diameter of 3.0 cm!"
      << endl;
    return -1;
  }

  initPinCell(pitch);

  srand(1);
  cout << normRand() << endl;
  
  material *water = new material();
  material *fuel = new material();

  cout << Watt() << " MeV" << endl;
  return 0;

  int batch_size = 1E5;
  double En[batch_size];
  double xyz[3][batch_size];
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
    xyz[0][i] = r*cos(gamma);
    xyz[1][i] = r*sin(gamma);
    // sample an axial location
    xyz[2][i] = 100.0*normRand(); 
  }
  
  // outer loop over power iterations
  boolean converged = false;
  double ShannonEntropy = 0.0, PrevEntropy = 0.0;
  double tol_entropy = 1E-3;
  int k = 0; l = 0;
  while(k < max_iters){
    k = k+1; // total power iterations 
    
    // inner loop over the source bank
    while(
      
    }

    // use fission bank to create source bank for next iteration


    PrevEntropy = ShannonEntropy;
    ShannonEntropy = calcEntropy(xyz[0][],xyz[1][],xyz[2][]);
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
