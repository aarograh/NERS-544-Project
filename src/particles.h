// AUTHORS: Aaron Graham, Mike Jarrett
// PURPOSE: NERS 544 Course Project
// DATE   : April 30, 2015

#include<cstdlib>
#include<vector>

class particle{
  private:
    bool isAlive;
    int cellid;
    double position[3];
    double omega[3];
    double energy;
    double score;
    double cutoff;
    double survival;
    double totalXS,f235,f238,fH,fcap,fiss_frac,abs_frac;
  
  public:
    particle(const double[3], double, double, double, int);
    particle(double[3], double, double, double, int);
    double getCoord(int);
    int simulate();
    int simulate_implicit();
    bool roulette();
    double weight;
    double estimatorTL,squareTL,estimatorColl,squareColl;
    friend class fission;
};

class fission{
  private:
    double position[3];
    int cellid;
  public:
    fission(const particle&,int);
    friend void makeSource(std::vector<fission>&,std::vector<particle>&,int);
    friend double calcEntropy(std::vector<fission>);
};

void makeSource(std::vector<fission>&,std::vector<particle>&,int);
double calcEntropy(std::vector<fission> fissionBank);
