// AUTHORS: Aaron Graham, Mike Jarrett
// PURPOSE: NERS 544 Course Project
// DATE   : April 3, 2015

#include<vector>
#include<cmath>
#include<cstdlib>
#include "materials.h"
#include "utils.h"

std::vector<material*> materialList;

moderator::moderator(int matid)
{
  id = matid;

  // scattering cross sections
  mod_scat[0][0] = 2.0E+01; 
  mod_scat[0][1] = 3.0E-03; 
  mod_scat[0][2] = 1.2E+00; 
  mod_scat[1][0] = 4.0E+00; 
  mod_scat[1][1] = 1.5E-04; 
  mod_scat[1][2] =-6.0E-01; 

  // capture cross sections
  mod_cap[0] = 0.0E+00;
  mod_cap[0] = 8.0E-05;
  mod_cap[0] = 0.0E+00;

  // isotope number densities in atoms/(b*cm)
  moddens[0] = 6.6911E-02;
  moddens[1] = 3.3455E-02;

}
fuel::fuel(int matid)
{
  id = matid;

  // scattering cross sections
  fuel_scat[0][0] = 4.0E+00;
  fuel_scat[0][1] = 1.5E-04;
  fuel_scat[0][2] =-6.0E-01;
  fuel_scat[1][0] = 1.5E+01;
  fuel_scat[1][1] = 1.5E-04;
  fuel_scat[1][2] =-4.0E-01;
  fuel_scat[2][0] = 9.0E+00;
  fuel_scat[2][1] = 1.0E-04;
  fuel_scat[2][2] =-1.6E-01;

  // capture cross sections
  fuel_cap[0][0] = 4.0E-01;
  fuel_cap[0][1] = 2.5E-03;
  fuel_cap[0][2] =-1.0E+00;
  fuel_cap[1][0] = 1.8E+00;
  fuel_cap[1][1] = 4.0E-04;
  fuel_cap[1][2] =-1.5E+00;

  // fission cross sections
  U235_fiss[0] = 8.0E-01; 
  U235_fiss[1] = 6.0E-02; 
  U235_fiss[2] = 0.0E+00; 

  // resonance data
  nres = 3;
  // energies in MeV
  Eres[0] = 6.6E-06;
  Eres[1] = 2.2E-05;
  Eres[2] = 3.8E-05;
  // peak cross section in barns
  U238_res[0] = 7.0E+03;
  U238_res[1] = 6.0E+03;
  U238_res[2] = 6.5E+03;
  // resonance widths in MeV
  rwidth[0] = 4.0E-08;
  rwidth[1] = 3.0E-08;
  rwidth[2] = 1.0E-07;
  
  dres = 5; // practical width of resonance, we should play with this parameter

  // isotope number densities in atoms/(b*cm)
  fueldens[0] = 4.7284E-02;
  fueldens[1] = 9.4567E-04;
  fueldens[2] = 2.2696E-02;
}

void init_materials(int& fuelid, int& modid)
{
  fuelid = 0;
  materialList.push_back(new fuel(fuelid));
  modid = 1;
  materialList.push_back(new moderator(modid));
  return;
}
material* getPtr_material(int matid)
{
  if (matid < 0 || matid >= materialList.size())
  {
    std::cout << "Error returning material ptr.  Material id " <<
      matid << " is invalid." << std::endl;
    exit(-6);
  }
  return materialList.at(matid);
}

void moderator::modMacro(double E, double *totalxs, double *H_frac, double *abs_frac)
{
  double sqrE = sqrt(E);
  macabs_H = moddens[0]*mod_cap[1]/sqrE;
  macscat_H = moddens[0]*(mod_scat[0][0]+mod_scat[0][1]/sqrE)*exp(mod_scat[0][2]*sqrE);
  macscat_O = moddens[1]*(mod_scat[1][0]+mod_scat[1][1]/sqrE)*exp(mod_scat[1][2]*sqrE);

  *totalxs = macscat_H+macscat_O+macabs_H;
  *H_frac = (macabs_H+macscat_H)/(*totalxs);
  *abs_frac = macabs_H/(macabs_H+macscat_H);

   return;
}

void fuel::fuelMacro(double E, double *totalxs, double *frac_U235, double *frac_U238)
{
  double sqrE = sqrt(E);
  // check for proximity to a resonance
  res_xs = 0; 
  for(int j = 0; j < nres; j++){
    if(fabs(E-Eres[j]) > dres*rwidth[j]){
      y = (2/rwidth[j])*(E-Eres[j]);
      res_xs = fueldens[2]*U238_res[j]*sqrt(Eres[j]/E)/(1+y*y);
    }
  }

  macscat_O = fueldens[0]*(fuel_scat[0][0]+fuel_scat[0][1]/sqrE)*exp(fuel_scat[0][2]*sqrE);
  macscat_U235 = fueldens[1]*(fuel_scat[1][0]+fuel_scat[1][1]/sqrE)*exp(fuel_scat[1][2]*sqrE);
  macscat_U238 = fueldens[2]*(fuel_scat[2][0]+fuel_scat[2][1]/sqrE)*exp(fuel_scat[2][2]*sqrE);

  maccap_U235 = fueldens[1]*(fuel_cap[1][0]+fuel_cap[1][1]/sqrE)*exp(fuel_cap[1][2]*sqrE);
  maccap_U238 = fueldens[2]*(fuel_cap[2][0]+fuel_cap[2][1]/sqrE)*exp(fuel_cap[2][2]*sqrE)+res_xs;
  macfiss_U235 = fueldens[1]*(U235_fiss[0]+U235_fiss[1]/sqrE)*exp(U235_fiss[2]*sqrE);

  *totalxs = macscat_O+macscat_U235+macscat_U238+maccap_U235+maccap_U238;
  *frac_U235 = (macscat_U235+maccap_U235)/(*totalxs);
  *frac_U238 = (macscat_U238+maccap_U238)/(*totalxs);

  return;
}

int fuel::sample_U(double E, double *frac_U235, double *frac_U238, double *abs_frac, double *fiss_frac)
{
  double sqrE = sqrt(E);
  double xi = drand();
  if(xi < *frac_U235){
    macscat_U235 = fueldens[1]*(fuel_scat[1][0]+fuel_scat[1][1]/sqrE)*exp(fuel_scat[1][2]*sqrE);
    maccap_U235 = fueldens[1]*(fuel_cap[0][0]+fuel_cap[0][1]/sqrE)*exp(fuel_cap[0][2]*sqrE);
    macfiss_U235 = fueldens[1]*(U235_fiss[0]+U235_fiss[1]/sqrE)*exp(U235_fiss[2]*sqrE);

    *fiss_frac = macfiss_U235/(macscat_U235+maccap_U235+macfiss_U235);
    *abs_frac = (macfiss_U235+maccap_U235)/(macscat_U235+maccap_U235+macfiss_U235);
//    std::cout << "Capture fraction in 235 = " << *abs_frac << std::endl;
//    std::cout << "Fission fraction in 235 = " << *fiss_frac << std::endl;
    return 235;
  }
  else if(xi < (*frac_U235+*frac_U238)){
    // check for proximity to a resonance
    res_xs = 0; 
    for(int j = 0; j < nres; j++){
      if(fabs(E-Eres[j]) > dres*rwidth[j]){
        y = (2/rwidth[j])*(E-Eres[j]);
        res_xs = fueldens[2]*U238_res[j]*sqrt(Eres[j]/E)/(1+y*y);
      }
    }

    macscat_U238 = fueldens[2]*(fuel_scat[2][0]+fuel_scat[2][1]/sqrE)*exp(fuel_scat[2][2]*sqrE);
    maccap_U238 = fueldens[2]*(fuel_cap[1][0]+fuel_cap[1][1]/sqrE)*exp(fuel_cap[1][2]*sqrE)+res_xs;
    
    *fiss_frac = 0;
    *abs_frac = maccap_U238/(maccap_U238+macscat_U238);
//    std::cout << "Capture fraction in 238 = " << *abs_frac << std::endl;
    return 238;
  }
  else{ // capture in oxygen
    *fiss_frac = 0;
    *abs_frac = 0;
    return 16;
  }
}

void elastic(const double temp, int A_in, double *v_n, double d_n[3])
{
  double A = static_cast<double>(A_in);
  double beta = sqrt(neut_mass/(lightspeed*lightspeed)*A/(2*kB*temp)); 
//std::cout << "beta = " << beta << std::endl;

//std::cout << "temperature = " << temp << std::endl;
//std::cout << "isotope mass = " << A << std::endl;
//std::cout << "incoming velocity = " << *v_n << std::endl;
//std::cout << "directions = " << d_n[0] 
//          << " , " << d_n[1] << " , " << d_n[2] << std::endl;

  double x;
  double y = beta*(*v_n);
  double w1 = sqrt(pi)*y/(2 + sqrt(pi)*y);

  double eta = 1.0;
  double f1 = 0.0;
  double w, x1, x2, x3, Vtil, mutil;
  while(eta > f1){ // sample until Vtil, mutil are accepted
    w = drand(); 
    x1 = drand();
    x2 = drand();  
    if(w < w1){  // sample g1(x)
      x3 = drand(); 
      x = sqrt(-log(x1) - log(x2)*cos(x3*pi/2)*cos(x3*pi/2));
    }
    else{ // sample g2(x)
      x = sqrt(-log(x1*x2));
    }

    Vtil = x/beta; 
    mutil = 2*drand() - 1;
//std::cout << "Vtil = " << Vtil << std::endl;
//std::cout << "mutil = " << mutil << std::endl;

    // check for rejection from scaled f1(V,mu) (Lecture Module 8)
    eta = drand();   
    f1 = sqrt((*v_n)*(*v_n) + Vtil*Vtil - 2*(*v_n)*Vtil*mutil)/((*v_n)+Vtil);
//std::cout << "f1 = " << f1 << std::endl;
  }
  
  // sample direction vector for the target nucleus Omega_T-hat 
  double gamma = 2*pi*drand();
  double Tx = mutil*d_n[0] + (d_n[0]*d_n[2]*cos(gamma) - d_n[1]*sin(gamma)*sqrt((1-mutil*mutil)/(1-d_n[2]*d_n[2])));
  double Ty = mutil*d_n[1] + (d_n[1]*d_n[2]*cos(gamma) + d_n[0]*sin(gamma)*sqrt((1-mutil*mutil)/(1-d_n[2]*d_n[2])));
  double Tz = mutil*d_n[2] - cos(gamma)*sqrt((1-mutil*mutil)*(1-d_n[2]*d_n[2])); 

  // center-of-mass velocity u_xyz
  double ux = ((*v_n)*d_n[0] + A*Vtil*Tx)/(1+A);  
  double uy = ((*v_n)*d_n[1] + A*Vtil*Ty)/(1+A);  
  double uz = ((*v_n)*d_n[2] + A*Vtil*Tz)/(1+A);  
  double uu = sqrt(ux*ux + uy*uy + uz*uz); // center-of-mass speed

  // neutron center-of-mass velocity
  double vcx = (*v_n)*d_n[0] - ux;
  double vcy = (*v_n)*d_n[1] - uy;
  double vcz = (*v_n)*d_n[2] - uz;
  double vcn = sqrt(vcx*vcx + vcy*vcy + vcz*vcz);

  // neutron center-of-mass direction vector
  double ncx = vcx/vcn;
  double ncy = vcy/vcn;
  double ncz = vcz/vcn;

  // outgoing neutron center-of-mass direction 
  gamma = 2*pi*drand();
  double muc = 2*drand() - 1;
  double ncxp = muc*ncx + (ncx*ncz*cos(gamma) - ncy*sin(gamma))*sqrt((1-muc*muc)/(1-ncz*ncz)); 
  double ncyp = muc*ncy + (ncy*ncz*cos(gamma) + ncx*sin(gamma))*sqrt((1-muc*muc)/(1-ncz*ncz)); 
  double nczp = muc*ncz - cos(gamma)*sqrt((1-muc*muc)*(1-ncz*ncz)); 

  // finally, outgoing neutron velocity in lab frame is calculated
  double vncx = vcn*ncxp + ux;
  double vncy = vcn*ncyp + uy;
  double vncz = vcn*nczp + uz; 
  *v_n = sqrt(vncx*vncx + vncy*vncy + vncz*vncz);
  d_n[0] = vncx/(*v_n);
  d_n[1] = vncy/(*v_n);
  d_n[2] = vncz/(*v_n);

}
