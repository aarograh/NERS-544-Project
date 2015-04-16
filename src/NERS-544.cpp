// AUTHORS: Aaron Graham, Mike Jarrett
// PURPOSE: NERS 544 Course Project
// DATE   : April 3, 2015

#include<iostream>
#include<iomanip>
#include<cstdlib>
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
}
