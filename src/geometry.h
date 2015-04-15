// AUTHORS: Aaron Graham, Mike Jarrett
// PURPOSE: NERS 544 Course Project
// DATE   : April 3, 2015

//TODO Extend this type for planes, cylinders, or whatever else we need
class surface{
  public:
    int id;
};

class cell{
  public:
    int id; // cell id number
    int *iSurfs; // id numbers of surfaces which create this cell
    int *surfSens; // senses for each surface in iSurfs. 1 is +, 0 is -
    int matid;
    double distToSurf(double, double, double);
};
