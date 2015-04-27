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

  int batch_size = 1E2;
  double En;
  double xi;
//  double xyz[3];
  posit xyz;
  double pinrad = 1.5; // pin radius = 1.5 cm
  double r, gamma, mu;
  double sourceProb;
  vector<particle> sourceBank;
  vector<fission> fissionBank;
  

  // sample neutrons for initial source bank
  for(int i = 0; i < batch_size; i++){
    // sample energy from Watt spectrum
    En = Watt();
    // sample a radial location within the fuel cell
    gamma = 2*pi*drand();
    r = pinrad*sqrt(drand());
    xyz.x = r*cos(gamma);
    xyz.y = r*sin(gamma);
    // sample an axial location
    xyz.z = 100.0*drand(); 
    // sample a direction
    gamma = drand();
    mu = 2.0*drand() - 1.0;
    // add particle to the bank
    sourceBank.push_back(particle(xyz,gamma,mu,En,fuelid));
  }
  
  // outer loop over power iterations
  bool converged = false;
  double ShannonEntropy = 0.0, PrevEntropy = 0.0;
  double tol_entropy = 1E-3;
  const int max_iters = 100, active_iters = 10;
  int k = 0, l = 0, n = 0, fissions;
  particle neutron = sourceBank.back();

  while(k < max_iters){
    k = k+1; // total power iterations 
    
    // inner loop over the source bank
    for(int i = 0; i < sourceBank.size(); i++)
    {
      cout << "Source cellid = " << sourceBank[i].getID() << endl;
    }
    while(!sourceBank.empty())
    {
//      cout << "Cell id's in source bank: " << endl;
//      for(int i = 0; i < sourceBank.size(); i++)
//      {
//        cout << (*sourceBank.at(i)).getID() << endl;
//      }
      // Get pointer to particle
      neutron = sourceBank.back();
      cout << "Starting history in cell " << (sourceBank.back()).getID() << endl;
      cout << "Starting history in cell " << neutron.getID() << endl;
      // Simulate particle
      fissions = neutron.simulate();
      // Create fission neutrons (if fissions > 0)
      xyz = neutron.getCoord;
      for(int i = 0; i < fissions; i++)
      {
        fissionBank.push_back(fission(xyz,fuelid));
      }
      // Delete pointer to neutron in sourcebank;
      //delete sourceBank.back();
      sourceBank.pop_back();
      //delete neutron;
    }

    // Calculate Shannon Entropy
    PrevEntropy = ShannonEntropy;
    ShannonEntropy = calcEntropy(fissionBank);
    // some convergence check
    if(converged){
      l = l+1; // power iterations with converged source
      if (l >= active_iters) break;
    }
    else{
      // run some convergence check
      if((ShannonEntropy-PrevEntropy) < tol_entropy){
        converged = true;
        cout << "The fission source is converged." << endl;
      }
    }

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
