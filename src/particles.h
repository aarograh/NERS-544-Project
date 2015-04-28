// AUTHORS: Aaron Graham, Mike Jarrett
// PURPOSE: NERS 544 Course Project
// DATE   : April 3, 2015

#include<cstdlib>
#include<vector>

class particle{
  private:
    bool isAlive;
    int cellid;
    double position[3];
    double omega[3];
    double energy;
    double weight;
    double totalXS,f235,f238,fH,fcap,fiss_frac,abs_frac;
    double estimatorTL,squareTL,estimatorColl,squareColl,score;;
  
  public:
    particle(double[3], double, double, double, int);
    void moveParticle(double);
    int getID();
    double getCoord(int);
    double Direction(int);
    int simulate();
    double getTL();
    double getColl();
    double getTLsq();
    double getCollsq();
    double getWeight();
};
class fission{
  private:
    double position[3];
    int cellid;
  public:
    fission();
    fission(double[3],int);
    double getCoord(int);
};

fission fissionNeutron(particle neutron);
void makeSource(std::vector<fission>&,std::vector<particle>&,int);
double calcEntropy(std::vector<fission> fissionBank);
