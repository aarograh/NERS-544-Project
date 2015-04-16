// AUTHORS: Aaron Graham, Mike Jarrett
// PURPOSE: NERS 544 Course Project
// DATE   : April 3, 2015

#include<vector>

class surface{
  public:
    int id;
    surface();
};

class plane : surface{
  public:
    plane();
};

class cylinder : surface{
  public:
    cylinder();
};

class cell{
  public:
    int id; // cell id number
    std::vector<int> iSurfs; // id numbers of surfaces which create this cell
    std::vector<int> surfSens; // senses for each surface in iSurfs. 1 is +, 0 is -
    int matid;
    cell();
    double distToSurf(double, double, double);
};

void initPinCell(double);
