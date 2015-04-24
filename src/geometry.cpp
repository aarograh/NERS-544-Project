// AUTHORS: Aaron Graham, Mike Jarrett
// PURPOSE: NERS 544 Course Project
// DATE   : April 3, 2015

#include<algorithm>
#include<iostream>
#include<cstdlib>
#include<cmath>
#include "utils.h"
#include "geometry.h"

std::vector<surface*> surfaceList;
std::vector<cell*> cellList;

const int xplane = 1, yplane = 2, zplane = 3;

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
    default:
      std::cout << "Error when constructing plane!" << std::endl;
      abort();
  }
}

double plane::distToIntersect(double position[3], double direction[3],
    double intersection[3])
{
  double distance = -1.0;
  double prod = norm[0]*direction[0] + norm[1]*direction[1] +
    norm[2]*direction[2];

  if (prod > 0.0 || prod < 0.0) //TODO: Implement an EPSILON
  {
    distance = ((point[0] - position[0])*norm[0] +
      (point[1] - position[1])*norm[1] + (point[2] - position[2])*norm[2])/
      prod;
  }

  return distance;
}

cylinder::cylinder(int surfid, double x, double y, double z, double R)
{
  id = surfid;
  origin[0] = x; origin[1] = y; origin[2] = z;
  radius = R;
}

double cylinder::distToIntersect(double position[3], double direction[3],
    double intersection[3])
{
  double distance = -1.0;
  double incidence = -1.0;
  double x, y, a, b, c, dis, m;
  double xp, xm, yp, ym, d2o, tmp;

  if (approxeq(direction[0],0.0) && approxeq(direction[1],0.0)) return -1.0;

  x = position[0] - origin[0];
  y = position[1] - origin[1];
  m = direction[1]/direction[0];
  d2o = sqrt(x*x + y*y);

  a = 1.0 + m*m;
  b = -2.0*m*m*x + 2.0*y*m;
  c = y*y + m*m*x*x - 2.0*y*x*m - radius*radius;
  dis = b*b - 4.0*a*c;

  if (approxeq(direction[0],0.0))
  {
    yp = sin(acos(x/radius))*radius;
    ym = -yp;
    if ((d2o < radius && direction[1] > 0.0) || 
      (d2o > radius && direction[1] < 0.0 && y > 0.0 && fabs(x) <= radius))
    {
      intersection[1] = yp;
      intersection[0] = x;
      tmp = (intersection[0] - x)*(intersection[0] - x) +
        (intersection[1] - y)*(intersection[1] - y);
      intersection[2] = position[2] + sqrt(tmp)*direction[2]/
        sqrt(direction[0]*direction[0] + direction[1]*direction[1]);
      distance = sqrt(tmp + (intersection[2] - position[2])*
        (intersection[2] - position[2]));
    }
    else if ((d2o < radius && direction[1] < 0.0) ||
      (d2o > radius && direction[1] > 0.0 && y < 0.0 && fabs(x) <= radius))
    {
      intersection[1] = ym;
      intersection[0] = x;
      tmp = (intersection[0] - x)*(intersection[0] - x) +
        (intersection[1] - y)*(intersection[1] - y);
      intersection[2] = position[2] + sqrt(tmp)*direction[2]/
        sqrt(direction[0]*direction[0] + direction[1]*direction[1]);
      distance = sqrt(tmp + (intersection[2] - position[2])*
        (intersection[2] - position[2]));
    }
  }
  else if (dis > 0.0)
  {
    xp = (-b + sqrt(dis))/(2.0*a);
    xm = (-b - sqrt(dis))/(2.0*a);
    if ((x > xm && x < xp && direction[0] > 0.0) ||
      (x > xp && direction[0] < 0.0))
    {
      intersection[0] = xp;
      intersection[1] = m*(xp - x) + y;
      tmp = (intersection[0] - x)*(intersection[0] - x) +
        (intersection[1] - y)*(intersection[1] - y);
      intersection[2] = position[2] + sqrt(tmp)*direction[2]/
        sqrt(direction[0]*direction[0] + direction[1]*direction[1]);
      distance = sqrt(tmp + (intersection[2] - position[2])*
        (intersection[2] - position[2]));
    }
    else if ((x > xm && x < xp && direction[0] < 0.0) ||
      (x < xm && direction[0] > 0.0))
    {
      intersection[0] = xm;
      intersection[1] = m*(xm - x) + y;
      tmp = (intersection[0] - x)*(intersection[0] - x) +
        (intersection[1] - y)*(intersection[1] - y);
      intersection[2] = position[2] + sqrt(tmp)*direction[2]/
        sqrt(direction[0]*direction[0] + direction[1]*direction[1]);
      distance = sqrt(tmp + (intersection[2] - position[2])*
        (intersection[2] - position[2]));
    }
  }
  else if (approxeq(dis,0.0))
  {
    xp = -b/(2.0*a);
    if ((x < xp && direction[0] > 0.0) || (x > xp && direction[0] < 0.0))
    {
      intersection[0] = xp;
      intersection[1] = m*(xp - x) + y;
      tmp = (intersection[0] - x)*(intersection[0] - x) +
        (intersection[1] - y)*(intersection[1] - y);
      intersection[2] = position[2] + sqrt(tmp)*direction[2]/
        sqrt(direction[0]*direction[0] + direction[1]*direction[1]);
      distance = sqrt(tmp + (intersection[2] - position[2])*
        (intersection[2] - position[2]));
    }
  }

  intersection[0] += origin[0];
  intersection[1] += origin[1];

  return distance;
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

double cell::distToIntersect(double position[3], double direction[3],
  double intersection[3])
{
  double distance = -1.0;
  int surfid = 0;

  for (int i = 0; i < iSurfs.size(); i++)
  {
    surfid = iSurfs.at(i);
    double temp = 
      (*surfaceList.at(surfid)).distToIntersect(position,direction,intersection);
    if (temp > 0.0 && (temp < distance || distance < 0.0)) distance = temp;
  }

  return distance;
}

void initPinCell(double pitch, int fuelid, int modid)
{
  double halfpitch = pitch/2.0;
  double height = 100.0;
  double radius = 1.5;

// Build the fuel pin
  // Construct cylinder for fuel pin
  surfaceList.push_back(new cylinder(surfaceList.size()+1,0.0,0.0,0.0,1.5));
  // Construct top plane
  surfaceList.push_back(new plane(surfaceList.size()+1,100.0,zplane));
  // Construct bottom plane
  surfaceList.push_back(new plane(surfaceList.size()+1,0.0,zplane));
  // Construct the cell
  {
    int isurfs[3] = {0,1,2};
    int sense[3] = {-1,-1,1};
    cellList.push_back(new cell(fuelid,3,isurfs,sense));
  }

// Construct the "box" for the moderator
  // Construct left plane
  surfaceList.push_back(new plane(surfaceList.size()+1,-halfpitch,xplane));
  // Construct right plane
  surfaceList.push_back(new plane(surfaceList.size()+1,halfpitch,xplane));
  // Construct front plane
  surfaceList.push_back(new plane(surfaceList.size()+1,-halfpitch,yplane));
  // Construct back plane
  surfaceList.push_back(new plane(surfaceList.size()+1,halfpitch,yplane));
  // Construct the cell
  {
    int isurfs[7] = {0,1,2,3,4,5,6};
    int sense[7] = {1,-1,1,1,-1,1,-1};
    cellList.push_back(new cell(modid,7,isurfs,sense));
  }

  return;
}
