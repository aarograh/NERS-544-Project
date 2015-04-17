// AUTHORS: Aaron Graham, Mike Jarrett
// PURPOSE: NERS 544 Course Project
// DATE   : April 3, 2015

#include<vector>

class surface{
  public:
    int id;
    surface();
};

class plane : public surface{
  public:
    double point[3];
    double norm[3];
    plane(int, double, int);
};

class cylinder : public surface{
  public:
    double origin[3];
    double radius;
    cylinder(int, double, double, double, double);
};

class cell{
  public:
    int id; // cell id number
    std::vector<int> iSurfs; // id numbers of surfaces which create this cell
    std::vector<int> senses; // senses for each surface in iSurfs. 1 is +, -1 is -
    int matid;
    cell(int, int, int*, int*);
    double distToSurf(double, double, double);
};

void initPinCell(double);
