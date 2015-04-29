// AUTHORS: Aaron Graham, Mike Jarrett
// PURPOSE: NERS 544 Course Project
// DATE   : April 30, 2015

#include <iostream>
#include <fstream>
#include "utils.h"
#include "particles.h"
#include "geometry.h"
#include "materials.h"

using namespace std;

int main()
{
//  srand(time(NULL));
  srand(10);

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

  int batch_size = 1E4;
  double En;
  double xyz[3];
  double pinrad = 1.5; // pin radius = 1.5 cm
  double r, gamma, mu;
  vector<particle> sourceBank;
  vector<fission> fissionBank;
  
  // estimators
  double topCurrent = 0.0, bottomCurrent = 0.0;
  int topSurf = -1;
  int bottomSurf = -2;
  double tally_TL = 0.0, tally_TLsq = 0.0;
  double tally_coll = 0.0, tally_collsq = 0.0, topleaksq = 0.0;
  double bottomleaksq = 0.0;
  double keff_TL, keff_Coll, sigTL, sigColl, active_particles;
  double topleak, bottomleak, sigtop, sigbottom, score = 0.0;
  particle neutron = particle(xyz,gamma,mu,En,fuelid);

  // outer loop over power iterations
  const int max_iters = 100, active_iters = 80, inactive_iters = 20;
  double ShannonEntropy[max_iters];
  double totalEntropy = 0.0, meanEntropy;
  int k = 0, l = 0, result, ktot = 0;

  const int npitch = 10;
  double pitches[npitch]  = {3.25,3.5,3.75,4.0,4.25,4.5,4.75,5.0,5.5,6.0};

  ofstream myfile;
  myfile.open("keff.out");

  for(int p = 0; p < npitch; p++)
  {
    pitch = pitches[p];

    pitch = 4.0;
    initPinCell(pitch, fuelid, modid);
    // zero all estimators
    k = 0;
    l = 0;
    ktot = 0;
    totalEntropy = 0.0;
    tally_TL = 0;
    tally_coll = 0;
    tally_TLsq = 0;
    tally_collsq = 0;
    topleak = 0;
    bottomleak = 0;
    sigtop = 0.0;
    sigbottom = 0.0;

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
  
    while(k < max_iters){
      k = k+1; // total power iterations 
    
      // inner loop over the source bank
      while(!sourceBank.empty())
      {
        // Get pointer to particle
        neutron = sourceBank.back();
        // Simulate particle
//        result = neutron.simulate_implicit();
        result = neutron.simulate();
        // Create fission neutrons (if fissions > 0)
        if(result > 0)
        { 
          for(int i = 0; i < result; i++)
          {
            fissionBank.push_back(fission(neutron,fuelid));
          }
        }

        // get keff tallies for the history
        if(k > inactive_iters)
        {
          tally_TL = tally_TL + neutron.estimatorTL;
          tally_coll = tally_coll + neutron.estimatorColl;
          tally_TLsq = tally_TLsq + neutron.estimatorTL*neutron.estimatorTL;
          tally_collsq = tally_collsq + neutron.estimatorColl*neutron.estimatorColl;

          // Calculate Leakages
          if(result == topSurf)
          { 
            score = neutron.weight;
            topCurrent = topCurrent + score;
            topleaksq = topleaksq + score*score;
          }
          else if(result == bottomSurf)
          {
            score = neutron.weight;
            bottomCurrent = bottomCurrent + score;
            bottomleaksq = bottomleaksq + score*score;
          }
        }
        // Delete pointer to neutron in sourcebank;
        sourceBank.pop_back();
      }

      // Calculate Shannon Entropy
      ShannonEntropy[k-1] = calcEntropy(fissionBank);
      // let a few cycles go by before starting to calculate the mean
      if(k > inactive_iters)
      {
        totalEntropy = totalEntropy + ShannonEntropy[k-1];
        meanEntropy = totalEntropy/(double)(k-inactive_iters);
        l = l+1; // power iterations with converged source
      }

      // Calculations for output
      active_particles = static_cast<double>(batch_size*(k-inactive_iters));
      keff_TL = tally_TL/active_particles;
      sigTL = sqrt((tally_TLsq/active_particles - 
        keff_TL*keff_TL)/active_particles);
      keff_Coll = tally_coll/active_particles;
      cout << "tally_collsq = " << tally_collsq << endl;
      cout << "sigColl^2 = " << tally_collsq/active_particles - keff_Coll*keff_Coll << endl;
      sigColl = sqrt((tally_collsq/active_particles - 
        keff_Coll*keff_Coll)/active_particles);
      topleak = topCurrent/active_particles;
      sigtop = sqrt((topleaksq/active_particles - 
        topleak*topleak)/active_particles);
      bottomleak = bottomCurrent/active_particles;
      sigbottom = sqrt((bottomleaksq/active_particles - 
        bottomleak*bottomleak)/active_particles);
      ktot = ktot + fissionBank.size();

      // Do some output
      cout << "Source iteration: " << k << endl;
      cout << "rough keff estimate = " << (double)(ktot)/(double)(batch_size*k) 
        << endl;
      cout << "track length keff estimate = " << keff_TL << ", uncertainty = " 
        << sigTL << endl;
      cout << "collision keff estimate = " << keff_Coll << ", unceratainty = " 
        << sigColl << endl;
      cout << "Top leakage estimate = " << topleak << ", uncertainty = " << 
        sigtop << endl;
      cout << "Bottom leakage estimate = " << bottomleak << ", uncertainty = " 
        << sigbottom << endl;
      cout << "Shannon Entropy: " << ShannonEntropy[l] << endl;
      cout << "Active cycle: " << l << endl;
      cout << "Fission bank has " << fissionBank.size() << " neutrons." << endl;

      // Make Source Bank
      cout << "Making source bank from fission bank..." << endl;
      makeSource(fissionBank,sourceBank,batch_size);
      cout << "Source bank size = " << sourceBank.size() << endl << endl;
    }
    cout << "The pin pitch was " << pitch << endl;

    // Write to output file
    myfile << "Pin pitch = " << pitch << endl;
    myfile << "Active cycles: " << active_iters << endl;
    myfile << "Inactive cycles: " << inactive_iters << endl;
    myfile << "track length keff estimate = " << keff_TL << ", uncertianty = "
      << sigTL << endl;
    myfile << "collision keff estimate = " << keff_Coll << ", uncertianty = "
      << sigColl << endl;
    myfile << "Top leakage estimate = " << topleak << ", uncertainty = "
      << sigtop << endl;
    myfile << "Bottom leakage estimate = " << bottomleak << ", uncertainty = " 
      << sigbottom << endl;
    myfile << endl;

    // clean out old geometry
    clearGeom();
  }
  myfile.close();
  return 0;
}
