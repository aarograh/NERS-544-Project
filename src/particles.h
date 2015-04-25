// AUTHORS: Aaron Graham, Mike Jarrett
// PURPOSE: NERS 544 Course Project
// DATE   : April 3, 2015

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
    int simulate();
};

particle* fissionNeutron(particle* neutron);
