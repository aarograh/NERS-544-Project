#include "Materials.h"

Material::Material()
{  
   matid = -5;
   micabs = 0.0;
   mictot = 0.0;
   micfiss = 0.0;
   ndens = 0.0;
   nu = 0;
};
Material::Material(int ID, double abs, double tot, double fiss, double dens, int v)
{  
   matid = ID;
   micabs = abs;
   mictot = tot;
   micfiss = fiss;
   ndens = dens;
   nu = v;
};
Material::~Material()
{
};
