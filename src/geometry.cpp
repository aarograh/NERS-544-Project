// AUTHORS: Aaron Graham, Mike Jarrett
// PURPOSE: NERS 544 Course Project
// DATE   : April 3, 2015

#include "geometry.h"

std::vector<surface> surfaceList;
std::vector<cell> cellList;

surface::surface()
{
  id = 0;
}

plane::plane()
{
  id = 0;
}

cylinder::cylinder()
{
  id = 0;
}

cell::cell()
{
  id = 0;
  matid = 0;
}

double cell::distToSurf(double x, double y, double z)
{
  return 0.0;
};

void initPinCell(double pitch)
{
  double halfpitch = pitch/2.0;
  double height = 100.0;
  double radius = 1.5;

  // Construct left plane
  surfaceList.push_back(surface());

  return;
}
