// AUTHORS: Aaron Graham, Mike Jarrett
// PURPOSE: NERS 544 Course Project
// DATE   : April 30, 2015

#include "geometry.h"
#include "utils.h"

std::vector<surface*> surfaceList;
std::vector<cell*> cellList;

// Constructor for plane class
plane::plane(int surfid, double position, int orientation, int bound_in)
{
  // Set surface id
  id = surfid;
  // Set normal vector depending on which axis the plane intersects
  // Only planes which have normal vectors parallel to an axis are
  // supported right now.
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
  switch(bound_in)
  {
    case reflecting:
    case vacuum:
      boundaryType = bound_in;
      break;
    default:
      boundaryType = -1;
  }
}

// Function to calculate distance between a plane and a point with an
// associated direction vector
double plane::distToIntersect(double position[3], double direction[3])
{
  double distance = -1.0;
  // Take dot product of position vector and plane's normal vector
  double prod = norm[0]*direction[0] + norm[1]*direction[1] +
    norm[2]*direction[2];

  if (!approxeq(prod,0.0))
  {
    // Calculate the distance between position and the plane
    distance = ((point[0] - position[0])*norm[0] +
      (point[1] - position[1])*norm[1] + (point[2] - position[2])*norm[2])/
      prod;
  }

  return distance;
}

// Reflection routine for planar surface
void plane::reflect(double* point_in, double* direction_in)
{
  double dot_product = direction_in[0]*norm[0] + direction_in[1]*norm[1] +
    direction_in[2]*norm[2];
  direction_in[0] -= 2.0*dot_product*norm[0];
  direction_in[1] -= 2.0*dot_product*norm[1];
  direction_in[2] -= 2.0*dot_product*norm[2];

  // "Nudge" point into cell
  point_in[0] += direction_in[0]*eps;
  point_in[1] += direction_in[1]*eps;
  point_in[2] += direction_in[2]*eps;

  return;
}

// Routine to return the sense of a plane with respect to some point
int plane::getSense(double* position)
{
  double position_vec[3], dotproduct;

  position_vec[0] = position[0] - point[0];
  position_vec[1] = position[1] - point[1];
  position_vec[2] = position[2] - point[2];

  dotproduct = position_vec[0]*norm[0] + position_vec[1]*norm[1] +
    position_vec[2]*norm[2];
  if (dotproduct > 0)
  {
    return 1;
  }
  else
  {
    return -1;
  }
}

// Constructor for cylindrical surface
// It is assumed that the cylinder's axis begins at (0,0,0) and
// points along the z-axis
cylinder::cylinder(int surfid, double x, double y, double z, double R, 
    int bound_in)
{
  // Set values
  id = surfid;
  origin[0] = x; origin[1] = y; origin[2] = z;
  radius = R;
  switch(bound_in)
  {
    case reflecting:
    case vacuum:
      boundaryType = bound_in;
      break;
    default:
      boundaryType = -1;
  }
}

// Function to calculate distance between a cylindrical surface and a point
// with an association direction vector
double cylinder::distToIntersect(double position[3], double direction[3])
{
  double distance = -1.0;
  double intersection[3];
  double x, y, a, b, c, dis, m;
  double xp, xm, yp, ym, d2o, tmp;

  // If vector is parallel with z-axis, return
  if (approxeq(direction[0],0.0) && approxeq(direction[1],0.0)) return -1.0;

  // Shift system so that cylinder is center at origin (in x-y)
  x = position[0] - origin[0];
  y = position[1] - origin[1];
  // Calculate 2D slope of vector
  m = direction[1]/direction[0];
  // Calculate 2D distance from position to cylinder's origin
  d2o = sqrt(x*x + y*y);

  // Calculate a, b, c, and discriminant for quadratice formula
  a = 1.0 + m*m;
  b = -2.0*m*m*x + 2.0*y*m;
  c = y*y + m*m*x*x - 2.0*y*x*m - radius*radius;
  dis = b*b - 4.0*a*c;

  // If vector is pointing straight down in 2D, slope is 0, so
  // we treat this case specially
  if (softeq(direction[0],0.0,1.0e-4))
  {
    // Calculate 2 possible intersections
    yp = sin(acos(x/radius))*radius;
    ym = -yp;
    // If inside circle and going up, or outside circle and going down
    // to intersect circle:
    if ((d2o < radius && direction[1] > 0.0) || 
      (d2o > radius && direction[1] < 0.0 && y > 0.0 && fabs(x) <= radius))
    {
      // Select top y-value.  x-value is unchanged
      intersection[1] = yp;
      intersection[0] = x;
      tmp = (intersection[0] - x)*(intersection[0] - x) +
        (intersection[1] - y)*(intersection[1] - y);
      // Calculate z-value and distance
      intersection[2] = position[2] + sqrt(tmp)*direction[2]/
        sqrt(direction[0]*direction[0] + direction[1]*direction[1]);
      distance = sqrt(tmp + (intersection[2] - position[2])*
        (intersection[2] - position[2]));
    }
    // If inside circle and going down, or outside circle and going up
    // to intersect circle:
    else if ((d2o < radius && direction[1] < 0.0) ||
      (d2o > radius && direction[1] > 0.0 && y < 0.0 && fabs(x) <= radius))
    {
      // Select bottom y-value.  x-value is unchanged
      intersection[1] = ym;
      intersection[0] = x;
      tmp = (intersection[0] - x)*(intersection[0] - x) +
        (intersection[1] - y)*(intersection[1] - y);
      // Calculate z-value and distance
      intersection[2] = position[2] + sqrt(tmp)*direction[2]/
        sqrt(direction[0]*direction[0] + direction[1]*direction[1]);
      distance = sqrt(tmp + (intersection[2] - position[2])*
        (intersection[2] - position[2]));
    }
  }
  // If discriminant > 0, then line intersects circle twice
  else if (dis > 0.0)
  {
    // Calculate both x values
    xp = (-b + sqrt(dis))/(2.0*a);
    xm = (-b - sqrt(dis))/(2.0*a);
    // If the rightmost x-value is needed
    if ((x > xm && x < xp && direction[0] > 0.0) ||
      (x > xp && direction[0] < 0.0))
    {
      // Set x-value, calculate y-value using point-slop formula
      intersection[0] = xp;
      intersection[1] = m*(xp - x) + y;
      tmp = (intersection[0] - x)*(intersection[0] - x) +
        (intersection[1] - y)*(intersection[1] - y);
      // Calculate z-value and distance
      intersection[2] = position[2] + sqrt(tmp)*direction[2]/
        sqrt(direction[0]*direction[0] + direction[1]*direction[1]);
      distance = sqrt(tmp + (intersection[2] - position[2])*
        (intersection[2] - position[2]));
    }
    // If the leftmost x-value is needed
    else if ((x > xm && x < xp && direction[0] < 0.0) ||
      (x < xm && direction[0] > 0.0))
    {
      // Set x-value, calculate y-value using point-slope formula
      intersection[0] = xm;
      intersection[1] = m*(xm - x) + y;
      tmp = (intersection[0] - x)*(intersection[0] - x) +
        (intersection[1] - y)*(intersection[1] - y);
      // Calculate z-value and distance
      intersection[2] = position[2] + sqrt(tmp)*direction[2]/
        sqrt(direction[0]*direction[0] + direction[1]*direction[1]);
      distance = sqrt(tmp + (intersection[2] - position[2])*
        (intersection[2] - position[2]));
    }
  }
  // If discriminant == 0, there is only 1 intersection
  else if (approxeq(dis,0.0))
  {
    // Calculate x value
    xp = -b/(2.0*a);
    // Ensure that the vector is pointing the correct direction to
    // use this point
    if ((x < xp && direction[0] > 0.0) || (x > xp && direction[0] < 0.0))
    {
      // Set x-value, calculate y-value with point-slope formula
      intersection[0] = xp;
      intersection[1] = m*(xp - x) + y;
      tmp = (intersection[0] - x)*(intersection[0] - x) +
        (intersection[1] - y)*(intersection[1] - y);
      // Calculate z-value and distance
      intersection[2] = position[2] + sqrt(tmp)*direction[2]/
        sqrt(direction[0]*direction[0] + direction[1]*direction[1]);
      distance = sqrt(tmp + (intersection[2] - position[2])*
        (intersection[2] - position[2]));
    }
  }

  // Shift back to global coordinates
  intersection[0] += origin[0];
  intersection[1] += origin[1];

  // Return the distance value, which is -1.0 if there was no intersection
  return distance;
}

//Reflection routine for cylindrical surface
void cylinder::reflect(double* point_in, double* direction_in)
{
  double shifted[3], norm[3];

  // Shift the point of intersection to have a (0,0,0) origin
  shifted[0] = point_in[0] - origin[0];
  shifted[1] = point_in[1] - origin[1];
  shifted[2] = point_in[2] - origin[2];

  // Calculate normal vector at the point
  norm[0] = cos(shifted[0]/radius);
  norm[1] = cos(shifted[1]/radius);
  norm[2] = 0.0; // cylinder assumed to have axis in z-direction

  // Calculate reflection direction
  double dot_product = direction_in[0]*norm[0] + direction_in[1]*norm[1] +
    direction_in[2]*norm[2];
  direction_in[0] -= 2.0*dot_product*norm[0];
  direction_in[1] -= 2.0*dot_product*norm[1];
  direction_in[2] -= 2.0*dot_product*norm[2];

  // "Nudge" point into cell
  point_in[0] += direction_in[0]*eps;
  point_in[1] += direction_in[1]*eps;
  point_in[2] += direction_in[2]*eps;

  return;
}

// Routine to return the sense of a cylinder with respect to some point
int cylinder::getSense(double* position)
{
  double shifted[2], d2o;

  shifted[0] = position[0] - origin[0];
  shifted[1] = position[1] - origin[1];

  d2o = sqrt(shifted[0]*shifted[0] + shifted[1]*shifted[1]);
  if (d2o > radius)
  {
    return 1;
  }
  else
  {
    return -1;
  }
}

// Constructor for cell class
cell::cell(int cellid, int size, int* surfs, int* sense)
{
  id = cellid;
  matid = 0;
  // Loop over surfaces, adding the surface and its sense to vectors
  for (int i = 0; i < size; i++, surfs++, sense++)
  {
    iSurfs.push_back(*surfs);
    senses.push_back(*sense);
  }
}

// Calculate the distance to the nearest surface for a cell
double cell::distToIntersect(double* position, double* direction,
  double* intersection, int& surfIntersect)
{
  double distance = -1.0;
  int surfid = -1;

  // Loop over surfaces
  for (int i = 0; i < iSurfs.size(); i++)
  {
    // Get surface id
    surfid = iSurfs.at(i);
    // Calculate distancce to that surface
    double tmp = 
      surfaceList.at(surfid)->distToIntersect(position,direction);
    // If this distance is better than the best so far, assign it
    if (tmp > 0.0 && (tmp < distance || distance < 0.0)) 
    {
      distance = tmp;
      surfIntersect = surfid;
    }
  }

  intersection[0] = position[0] + direction[0]*distance;
  intersection[1] = position[1] + direction[1]*distance;
  intersection[2] = position[2] + direction[2]*distance;
  return distance;
}

cell* getPtr_cell(int cellid)
{
  if (cellid < 0 || cellid >= cellList.size())
  {
    std::cout << "Error returning cell ptr.  Cell id " <<
      cellid << " is invalid." << std::endl;
    exit(-4);
  }
  return cellList.at(cellid);
}

surface* getPtr_surface(int surfid)
{
  if (surfid < 0 || surfid >= surfaceList.size())
  {
    std::cout << "Error returning surface ptr.  Surface id " <<
      surfid << " is invalid." << std::endl;
    exit(-5);
  }
  return surfaceList.at(surfid);
}

// Function to initialize a pin cell
void initPinCell(double pitch, int fuelid, int modid)
{
  double halfpitch = pitch/2.0;
  double height = 100.0; // Height is hard-coded, but could easily be changed
  double radius = 1.5; // Radius is hard-coded, but could easily be changed

// Build the fuel pin
  // Construct cylinder for fuel pin
  surfaceList.push_back(new cylinder(surfaceList.size(),0.0,0.0,0.0,radius,
    -1));
  // Construct top plane
  surfaceList.push_back(new plane(surfaceList.size(),height,zplane,0));
  // Construct bottom plane
  surfaceList.push_back(new plane(surfaceList.size(),0.0,zplane,0));
  // Construct the cell
  {
    int isurfs[3] = {0,1,2};
    int sense[3] = {-1,-1,1};
    cellList.push_back(new cell(fuelid,3,isurfs,sense));
  }

// Construct the "box" for the moderator
  // Construct left plane
  surfaceList.push_back(new plane(surfaceList.size(),-halfpitch,xplane,1));
  // Construct right plane
  surfaceList.push_back(new plane(surfaceList.size(),halfpitch,xplane,1));
  // Construct front plane
  surfaceList.push_back(new plane(surfaceList.size(),-halfpitch,yplane,1));
  // Construct back plane
  surfaceList.push_back(new plane(surfaceList.size(),halfpitch,yplane,1));
  // Construct the cell
  {
    int isurfs[7] = {0,1,2,3,4,5,6};
    int sense[7] = {1,-1,1,1,-1,1,-1};
    cellList.push_back(new cell(modid,7,isurfs,sense));
  }

  return;
}

void clearGeom()
{
  while(!cellList.empty())
  {
    delete cellList.back();
    cellList.pop_back();
  }
  while(!surfaceList.empty())
  {
    delete surfaceList.back();
    surfaceList.pop_back();
  }
  return;
}

int getCellID(double* position)
{
  int surfid, j;
  int senses[surfaceList.size()];
  cell* cellptr;
  surface* surfptr;
  std::fill_n(senses, surfaceList.size(), 0);

  for (int i = 0; i < cellList.size(); i++)
  {
    cellptr = cellList.at(i);
    for (j = 0; j < cellptr->iSurfs.size(); j++)
    {
      surfid = cellptr->iSurfs[j];
      surfptr = getPtr_surface(surfid);
      // Surface has not been checked yet
      if (senses[surfid] == 0) senses[surfid] = surfptr->getSense(position);
      // position is on wrong side of surface
      if(!(senses[surfid] == cellptr->senses[j])) break;
    }
    // Checked all surfaces without a break, so return this cell
    if (j == cellptr->iSurfs.size()) return cellptr->id;
  }

  return -1;
}
