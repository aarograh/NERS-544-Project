// AUTHORS: Aaron Graham, Mike Jarrett
// PURPOSE: NERS 544 Course Project
// DATE   : April 3, 2015

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
    exit(-1);
  }

  init_materials(fuelid, modid);
  initPinCell(pitch, fuelid, modid);

  int batch_size = 1E5;
  int initCycle = 5;
  double En;
  double xi;
//  double xyz[3];
  double xyz[3];
  double pinrad = 1.5; // pin radius = 1.5 cm
  double r, gamma, mu;
  double sourceProb;
  vector<particle> sourceBank;
  vector<fission> fissionBank;
  
  // estimators
  int topCurrent = 0;
  int bottomCurrent = 0;
  int topSurf = -1;
  int bottomSurf = -2;
  double tally_TL = 0.0;
  double tally_coll = 0.0;
  double keff_TL, keff_Coll;

//  double fuelVolume = pi*pinrad*pinrad*100.0;

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
    sourceBank.push_back(particle(xyz,gamma,mu,En,fuelid));
  }
  
  // outer loop over power iterations
  //bool converged = false;
  //double tol_entropy = 1E-1;
  const int max_iters = 100, active_iters = 10, inactive_iters = 0;
  double ShannonEntropy[max_iters];
  double totalEntropy = 0.0, meanEntropy = 0.0;
  int k = 0, l = 0, result, ktot = 0;
  particle neutron = sourceBank.back();

  while(k < max_iters){
    k = k+1; // total power iterations 
    
    // inner loop over the source bank
    while(!sourceBank.empty())
    {
      // Get pointer to particle
      neutron = sourceBank.back();
      // Simulate particle
      result = neutron.simulate();
      // Create fission neutrons (if fissions > 0)
      if(result == topSurf)
      { 
        topCurrent = topCurrent + 1;
      }
      else if(result == bottomSurf)
      {
        bottomCurrent = bottomCurrent + 1;
      }
      else if(result > 0)
      { for(int j = 0; j < 3; j++)
        {
          xyz[j] = neutron.getCoord(j);
        }
        for(int i = 0; i < result; i++)
        {
          fissionBank.push_back(fission(xyz,fuelid));
        }
      }

      // get keff tallies for the history
      tally_TL = tally_TL + neutron.getTL();
      tally_coll = tally_coll + neutron.getColl();

      // Delete pointer to neutron in sourcebank;
      //delete sourceBank.back();
      sourceBank.pop_back();
      //delete neutron;
    }

    // Calculate Shannon Entropy
    ShannonEntropy[k] = calcEntropy(fissionBank);
    // some convergence check
    // let 5 cycles go by before starting to calculate the mean
    if(k > inactive_iters)
    {
      totalEntropy = totalEntropy + ShannonEntropy[k];
      meanEntropy = totalEntropy/(double)(k-initCycle);
      l = l+1; // power iterations with converged source
      if (l >= active_iters) break;
    }
//    if(converged){
//    }
//    else{
//      // run some convergence check
//      if(fabs(ShannonEntropy-PrevEntropy) < tol_entropy){
//        converged = true;
//        cout << "The fission source is converged." << endl;
//      }
//    }

    cout << "Source iteration: " << k << endl;
    ktot = ktot + fissionBank.size();
    cout << "rough keff estimate = " << (double)(ktot)/(double)(batch_size*k) << endl;
    cout << "track length keff estimate = " << tally_TL/(double)(batch_size*(k-inactive_iters)) << endl;
    cout << "collision keff estimate = " << tally_coll/(double)(batch_size*(k-inactive_iters)) << endl;
    cout << "Top leakage estimate = " << (double)(topCurrent)/(double)(batch_size*k) << endl;
    cout << "Bottom leakage estimate = " << (double)(topCurrent)/(double)(batch_size*k) << endl;
    cout << "Shannon Entropy: " << ShannonEntropy[k] << endl;
    cout << "Active cycle: " << l << endl;
    cout << "Fission bank has " << fissionBank.size() << " neutrons." << endl;
    cout << "Making source bank from fission bank..." << endl;
    makeSource(fissionBank,sourceBank,batch_size);
    cout << "Source bank size = " << sourceBank.size() << endl;
    for(int i = 0; i < 500000000; i++)
    {
    }
  }
  return 0;
}
