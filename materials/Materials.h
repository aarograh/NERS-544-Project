// Header ==> function declarations

#include <iostream>
#include <string>

using namespace std;

#ifndef MATERIALS_H
#define MATERIALS_H

class Material
{
    public:
        double micabs_H[3];
        double micabs_U235[3];
        double micabs_U238[3];
        double micfiss_U235[3];

        double micscat_H[3];
        double micscat_O[3];
        double micscat_U235[3];
        double micscat_U238[3];

        // set coefficients for cross section description
        micabs_H[0] = 0.0E+00;
        micabs_H[1] = 8.0E-05;
        micabs_H[2] = 0.0E+00;
 
        micabs_U235[0] = 4.0E-01;
        micabs_U235[1] = 2.5E-03;
        micabs_U235[2] =-1.0E+00;

        micabs_U238[0] = 1.8E+00;
        micabs_U238[1] = 4.0E-04;
        micabs_U238[2] =-1.5E+00;

        micscat_H[0] = 2.0E+01;
        micscat_H[1] = 3.0E-03;
        micscat_H[2] =-1.2E+00;
 
        micscat_O[0] = 4.0E+00;
        micscat_O[1] = 1.5E-04;
        micscat_O[2] =-6.0E-01;

        micscat_U235[0] = 1.5E+01;
        micscat_U235[1] = 1.5E-04;
        micscat_U235[2] =-4.0E-01;

        micscat_U238[0] = 9.0E+00;
        micscat_U238[1] = 1.0E-04;
        micscat_U238[2] =-1.6E-01;

        fueldens_O =    4.7284E-02;
        fueldens_U235 = 9.4567E-04; 
        fueldens_U238 = 2.2696E-02; 

        moddens_H = 6.6911E-02;
        moddens_O = 3.3455E-02;
};
#endif
