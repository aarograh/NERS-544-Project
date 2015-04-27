// AUTHORS: Aaron Graham, Mike Jarrett
// PURPOSE: NERS 544 Course Project
// DATE   : April 3, 2015

#ifndef PARTICLES_H
#define PARTICLES_H

#include<cstdlib>
#include<vector>

class particle{
  private:
    bool isAlive;
    int cellid;
    posit position;
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
    particle(posit, double, double, double, int);
    void moveParticle(double);
    int getID(void);
    posit getCoord(void);
    double Direction(int);
    int simulate();
};
class fission: public particle
{
  public:
    fission();
    fission(posit,int);
};

fission fissionNeutron(particle neutron);
void makeSource(std::vector<fission>,std::vector<particle>,int);
double calcEntropy(std::vector<fission> fissionBank);
#endif
