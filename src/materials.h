// AUTHORS: Aaron Graham, Mike Jarrett
// PURPOSE: NERS 544 Course Project
// DATE   : April 3, 2015

class material{
  // arrays for isotopes, routines to get cross sections and calculate interaction stuff
  int id;
  int nres;
  double moddens[2]; // 0 = H, 1 = O
  double mod_scat[2][3]; // 0 = H, 1 = O
  double mod_cap[3]; // H only 
  double fueldens[3]; // 0 = U235, 1 = U238, 2 = O
  double fuel_scat[3][3]; // 0 = O, 1 = U235, 2 = U238
  double fuel_cap[2][3]; // 0 = U235, 1 = U238
  double U235_fiss[3];
  double U238_res[3];
  double Eres[3];
  double rwidth[3];
  double dres;

  double macabs_H,macscat_H,macscat_O;
  double macscat_U235,macscat_U238,maccap_U235,maccap_U238,macfiss_U235;
  double res_xs, y;
  
  material init_water();
  material init_fuel();

  double mod_macro(double,double*,double*,double*);
  double fuel_macro(double,double*,double*,double*);
  double sample_U(int,double,double*,double*);
  
  public:
    material();
};
