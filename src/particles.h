// AUTHORS: Aaron Graham, Mike Jarrett
// PURPOSE: NERS 544 Course Project
// DATE   : April 3, 2015

#include<cstdlib>
#include<vector>

class particle{
  private:
    int cellid;
    double position[3];
    double omega[3];
    double energy;
    double weight;
  public:
    bool isAlive;
    particle(double[3], double, double, double, int);
    double Coordinate(int);
    double Direction(int);
    int simulate();
};

particle* fissionNeutron(particle* neutron);
void makeSource(std::vector<particle*>,std::vector<particle*>,int);
double calcEntropy(std::vector<particle*> fissionBank);
