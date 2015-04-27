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
    double totalXS;
    double f235; 
    double f238; 
    double fH; 
    double fcap; 
    double fiss_frac; 
    double abs_frac; 
  public:
    particle(double[3], double, double, double, int);
    void moveParticle(double);
    int getID(void);
    double Coordinate(int);
    double Direction(int);
    int simulate();
};

particle* fissionNeutron(particle* neutron);
void makeSource(std::vector<particle*>*,std::vector<particle*>*,int);
double calcEntropy(std::vector<particle*> fissionBank);
