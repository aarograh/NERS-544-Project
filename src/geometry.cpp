// AUTHORS: Aaron Graham, Mike Jarrett
// PURPOSE: NERS 544 Course Project
// DATE   : April 3, 2015

#include "geometry.h"

std::vector<surface> surfaceList;
std::vector<cell> cellList;

const int xplane = 1, yplane = 2, zplane = 3;

surface::surface()
{
  id = 0;
}

plane::plane(int surfid, double position, int orientation)
{
  id = surfid;
  switch(orientation)
  {
    case xplane:
      point[0] = position; point[1] = 0.0; point[2] = 0.0;
      norm[0] = 1.0; norm[1] = 0.0; norm[2] = 0.0;
      break;
    case yplane:
      point[0] = 0.0; point[1] = position; point[2] = 0.0;
      norm[0] = 0.0; norm[1] = 1.0; norm[2] = 0.0;
      break;
    case zplane:
      point[0] = 0.0; point[1] = 0.0; point[2] = position;
      norm[0] = 0.0; norm[1] = 0.0; norm[2] = 1.0;
      break;
  }
}

cylinder::cylinder(int surfid, double x, double y, double z, double R)
{
  id = surfid;
  origin[0] = x; origin[1] = y; origin[2] = z;
  radius = R;
}

cell::cell(int cellid, int size, int* surfs, int* sense)
{
  id = cellid;
  matid = 0;
  for (int i = 0; i < size; i++, surfs++, sense++)
  {
    iSurfs.push_back(*surfs);
    senses.push_back(*sense);
  }
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

// Build the fuel pin
  // Construct cylinder for fuel pin
  surfaceList.push_back(cylinder(surfaceList.size()+1,0.0,0.0,0.0,1.5));
  // Construct top plane
  surfaceList.push_back(plane(surfaceList.size()+1,100.0,zplane));
  // Construct bottom plane
  surfaceList.push_back(plane(surfaceList.size()+1,0.0,zplane));
  // Construct the cell
  {
    int isurfs[3] = {0,1,2};
    int sense[3] = {-1,-1,1};
    cellList.push_back(cell(cellList.size()+1,3,isurfs,sense));
  }

// Construct the "box" for the moderator
  // Construct left plane
  surfaceList.push_back(plane(surfaceList.size()+1,-halfpitch,xplane));
  // Construct right plane
  surfaceList.push_back(plane(surfaceList.size()+1,halfpitch,xplane));
  // Construct front plane
  surfaceList.push_back(plane(surfaceList.size()+1,-halfpitch,yplane));
  // Construct back plane
  surfaceList.push_back(plane(surfaceList.size()+1,halfpitch,yplane));
  // Construct the cell
  {
    int isurfs[7] = {0,1,2,3,4,5,6};
    int sense[7] = {1,-1,1,1,-1,1,-1};
    cellList.push_back(cell(cellList.size()+1,7,isurfs,sense));
  }

  return;
}
