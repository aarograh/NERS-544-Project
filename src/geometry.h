// AUTHORS: Aaron Graham, Mike Jarrett
// PURPOSE: NERS 544 Course Project
// DATE   : April 3, 2015

#include<vector>

const int xplane = 1, yplane = 2, zplane = 3;
const int interior = -1, vacuum = 0, reflecting = 1;

class surface{
  public:
    int boundaryType;
    int id;
    virtual double distToIntersect(double*, double*) = 0;
    virtual void reflect(double*, double*) = 0;
    virtual int getSense(double*) = 0;
};

class plane : public surface{
  public:
    double point[3];
    double norm[3];
    plane(int, double, int, int);
    double distToIntersect(double*, double*);
    void reflect(double*, double*);
    int getSense(double*);
};

class cylinder : public surface{
  public:
    double origin[3];
    double radius;
    cylinder(int, double, double, double, double, int);
    double distToIntersect(double*, double*);
    void reflect(double*, double*);
    int getSense(double*);
};

class cell{
  public:
    int id; // cell id number
    std::vector<int> iSurfs; // id numbers of surfaces which create this cell
    std::vector<int> senses; // senses for each surface in iSurfs. 1 is +, -1 is -
    int matid;
    cell(int, int, int*, int*);
    double distToIntersect(double*, double*, double*, int&);
};

cell* getPtr_cell(int);
surface* getPtr_surface(int);

void initPinCell(double, int, int);
int getCellID(double*);
