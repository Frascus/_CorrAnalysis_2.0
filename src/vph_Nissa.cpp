#include "../include/vph_Nissa.h"
#include "numerics.h"

using namespace std;

const double M2PiPhys=pow(0.135,2);
const double alpha = 1/137.04;
const double e2 = alpha*4.0*M_PI;
const int Nboots= 2000;
const int NJ=30;
const double qu = 2.0/3.0; //electric charge of u-type quark
const double qd = -1.0/3.0; //electric charge of d-type quark
const string Meson="Ds";
double Lambda_QCD= 0.3; //300 MeV
bool Is_reph=true;
Vfloat xg_t_list;
Vfloat xg_to_spline;
Vfloat xg_to_spline_VMD;
Vfloat a_to_print;
const int shift=0;
int size_mu_nu= Is_reph?2:4;
string ph_type = Is_reph ? "rph" : "vph";
const double fDs_star = 0.2688;
const double fDs_star_err = 0.0066;
bool Enhanced_xg_only_B64 = false;
const string parameterization_to_be_used_in_Br = "single_pole";
double Vcs_FLAG = 0.9741;
double dVcs_FLAG= 0.0065;
const double re = (0.000510998950 /1.96835); // electron/Ds mass
const double rm = (0.105658/1.96835);
const double rt = (1.77686/1.96835);

// FLAG PARAMETERS
const double fDs_FLAG = 0.2499; //http://flag.unibe.ch/2021/Media?action=AttachFile&do=get&target=FLAG_2023_webupdate.pdf Eq. 179
const double fDs_FLAG_err = 0.0005;





void Get_xg_t_list(int num_xg) {

  for(int ixg=1;ixg<num_xg;ixg++) xg_t_list.push_back( ixg*0.1);
  return;
}



void Get_lattice_spacings_to_print() {

  int Nlat=300;
  double sx= 0.08*1.5/(Nlat-1.0); //fm
  for(int a=0;a<Nlat;a++) a_to_print.push_back( sx*a);


}

void Get_xg_to_spline() {

  int Nxgs=300;
  double sx = 0.9/(Nxgs-1.0);
  for(int ixg=0;ixg<Nxgs;ixg++) { xg_to_spline.push_back(0.1+sx*ixg);}
    
  return;
}

void Get_xg_to_spline_VMD() {

  int Nxgs=300;
  double sx = 1.0/(Nxgs-1.0);
  for(int ixg=0;ixg<Nxgs;ixg++) { xg_to_spline_VMD.push_back(0.0+sx*ixg);}
    
  return;
}

void Get_Tmin_Tmax(string W, int &Tmin, int &Tmax, int ixg, string Ens) {

   
 if(Ens == "cA211a.12.48") {

    if(ixg==0) {
      if(W=="A") {Tmin= 16; Tmax=22;}
      else if(W=="V") {Tmin=20;Tmax=32;}
      else if(W=="Ad") { Tmin=25; Tmax=33;}
      else if(W=="Vd") { Tmin=31; Tmax=42;}
      else if(W=="Au") { Tmin=19; Tmax=23;}
      else if(W=="Vu") { Tmin=30; Tmax=36;}
      else crash("time intervals for fits in channel: "+W+" not implemented");

    }
    else if(ixg==1) {
      if(W=="A") {Tmin= 16; Tmax=22;}
      else if(W=="V") {Tmin=20;Tmax=32;}
      else if(W=="Ad") { Tmin=25; Tmax=33;}
      else if(W=="Vd") { Tmin=31; Tmax=42;}
      else if(W=="Au") { Tmin=19; Tmax=23;}
      else if(W=="Vu") { Tmin=30; Tmax=36;}
      else crash("time intervals for fits in channel: "+W+" not implemented");

    }
    else if(ixg==2) {
      if(W=="A") {Tmin= 17; Tmax=23;}
      else if(W=="V") {Tmin=16;Tmax=25;}
      else if(W=="Ad") {Tmin=25; Tmax=33;}
      else if(W=="Vd") { Tmin=25; Tmax=38;}
      else if(W=="Au") { Tmin=19; Tmax=23;}
      else if(W=="Vu") { Tmin=26; Tmax=36;}
      else crash("time intervals for fits in channel: "+W+" not implemented");

    }
    else if(ixg==3) {
      if(W=="A") {Tmin= 17; Tmax=23;}
      else if(W=="V") {Tmin=18;Tmax=28;}
      else if(W=="Ad") {Tmin=25; Tmax=33;}
      else if(W=="Vd") { Tmin=24; Tmax=38;}
      else if(W=="Au") { Tmin=19; Tmax=23;}
      else if(W=="Vu") { Tmin=26; Tmax=36;}
      else crash("time intervals for fits in channel: "+W+" not implemented");
    
    }
    else if(ixg==4) {
      if(W=="A") {Tmin= 17; Tmax=23;}
      else if(W=="V") {Tmin=18;Tmax=28;}
      else if(W=="Ad") {Tmin=25; Tmax=33;}
      else if(W=="Vd") { Tmin=24; Tmax=34;}
      else if(W=="Au") { Tmin=19; Tmax=23;}
      else if(W=="Vu") { Tmin=26; Tmax=36;}
      else crash("time intervals for fits in channel: "+W+" not implemented");
    
    }
    else if(ixg==5) {
      if(W=="A") {Tmin= 16; Tmax=23;}
      else if(W=="V") {Tmin=18;Tmax=26;}
      else if(W=="Ad") {Tmin=25; Tmax=33;}
      else if(W=="Vd") { Tmin=24; Tmax=34;}
      else if(W=="Au") { Tmin=19; Tmax=23;}
      else if(W=="Vu") { Tmin=25; Tmax=36;}
      else crash("time intervals for fits in channel: "+W+" not implemented");
    
    }
    else if(ixg==6) {
      if(W=="A") {Tmin= 17; Tmax=26;}
      else if(W=="V") {Tmin=17;Tmax=26;}
      else if(W=="Ad") { Tmin=25; Tmax=33;}
      else if(W=="Vd") { Tmin=24; Tmax=34;}
      else if(W=="Au") { Tmin=19; Tmax=23;}
      else if(W=="Vu") { Tmin=25; Tmax=36;}
      else crash("time intervals for fits in channel: "+W+" not implemented");
    
    }
    else if(ixg==7) {
      if(W=="A") {Tmin= 16; Tmax=24;}
      else if(W=="V") {Tmin=16;Tmax=27;}
      else if(W=="Ad") { Tmin=23; Tmax=33;}
      else if(W=="Vd") { Tmin=24; Tmax=34;}
      else if(W=="Au") { Tmin=19; Tmax=23;}
      else if(W=="Vu") { Tmin=25; Tmax=36;}
      else crash("time intervals for fits in channel: "+W+" not implemented");
    
    }
    else if(ixg==8) {
      if(W=="A") {Tmin= 17; Tmax=29;}
      else if(W=="V") {Tmin=15;Tmax=26;}
      else if(W=="Ad") { Tmin=19; Tmax=33;}
      else if(W=="Vd") { Tmin=24; Tmax=34;}
      else if(W=="Au") { Tmin=19; Tmax=23;}
      else if(W=="Vu") { Tmin=25; Tmax=36;}
      else crash("time intervals for fits in channel: "+W+" not implemented");
    
    }
    else if(ixg==9) {
      if(W=="A") {Tmin= 15; Tmax=29;}
      else if(W=="V") {Tmin=15;Tmax=26;}
      else if(W=="Ad") { Tmin=14; Tmax=34;}
      else if(W=="Vd") { Tmin=24; Tmax=34;}
      else if(W=="Au") { Tmin=19; Tmax=23;}
      else if(W=="Vu") { Tmin=25; Tmax=36;}
      else crash("time intervals for fits in channel: "+W+" not implemented");
    
    }
    else if(ixg==10) {
      if(W=="A") {Tmin= 14; Tmax=25;}
      else if(W=="V") {Tmin=15;Tmax=27;}
      else if(W=="Ad") { Tmin=14; Tmax=34;}
      else if(W=="Vd") { Tmin=24; Tmax=34;}
      else if(W=="Au") { Tmin=19; Tmax=23;}
      else if(W=="Vu") { Tmin=25; Tmax=36;}
      else crash("time intervals for fits in channel: "+W+" not implemented");

    }
    else crash("ixg: "+to_string(ixg)+" does not have an established fit range");

  }


  

   else if( (Ens == "cB211b.072.64") || (Ens == "cB211b.072.96")) {

    if(ixg==0) {
      if(W=="A") {Tmin= 20; Tmax=30;}
      else if(W=="V") {Tmin=20;Tmax=30;}
      else if(W=="Ad") { Tmin=20; Tmax=32;}
      else if(W=="Vd") { Tmin=20; Tmax=32;}
      else if(W=="Au") { Tmin=20; Tmax=32;}
      else if(W=="Vu") { Tmin=20; Tmax=32;}
      else crash("time intervals for fits in channel: "+W+" not implemented");

    }
    else if(ixg==1) {
      if(W=="A") {Tmin= 18; Tmax=26;}
      else if(W=="V") {Tmin=23;Tmax=38;}
      else if(W=="Ad") { Tmin=36; Tmax=46;}
      else if(W=="Vd") { Tmin=35; Tmax=43;}
      else if(W=="Au") { Tmin=26; Tmax=33;}
      else if(W=="Vu") { Tmin=35; Tmax=41;}
      else crash("time intervals for fits in channel: "+W+" not implemented");

    }
    else if(ixg==2) {
      if(W=="A") {Tmin= 19; Tmax=27;}
      else if(W=="V") {Tmin=22;Tmax=34;}
      else if(W=="Ad") { Tmin=40; Tmax=46;}
      else if(W=="Vd") { Tmin=32; Tmax=41;}
      else if(W=="Au") { Tmin=26; Tmax=33;}
      else if(W=="Vu") { Tmin=36; Tmax=41;}
      else crash("time intervals for fits in channel: "+W+" not implemented");

    }
    else if(ixg==3) {
      if(W=="A") {Tmin= 19; Tmax=26;}
      else if(W=="V") {Tmin=21;Tmax=34;}
      else if(W=="Ad") { Tmin=40; Tmax=46;}
      else if(W=="Vd") { Tmin=29; Tmax=41;}
      else if(W=="Au") { Tmin=26; Tmax=33;}
      else if(W=="Vu") { Tmin=36; Tmax=41;}
      else crash("time intervals for fits in channel: "+W+" not implemented");
    
    }
    else if(ixg==4) {
      if(W=="A") {Tmin= 20; Tmax=27;}
      else if(W=="V") {Tmin=21;Tmax=33;}
      else if(W=="Ad") { Tmin=40; Tmax=46;}
      else if(W=="Vd") { Tmin=29; Tmax=39;}
      else if(W=="Au") { Tmin=26; Tmax=33;}
      else if(W=="Vu") { Tmin=36; Tmax=41;}
      else crash("time intervals for fits in channel: "+W+" not implemented");
    
    }
    else if(ixg==5) {
      if(W=="A") {Tmin= 19; Tmax=27;}
      else if(W=="V") {Tmin=20;Tmax=34;}
      else if(W=="Ad") { Tmin=40; Tmax=46;}
      else if(W=="Vd") { Tmin=29; Tmax=39;}
      else if(W=="Au") { Tmin=26; Tmax=33;}
      else if(W=="Vu") { Tmin=36; Tmax=41;}
      else crash("time intervals for fits in channel: "+W+" not implemented");
    
    }
    else if(ixg==6) {
      if(W=="A") {Tmin= 20; Tmax=28;}
      else if(W=="V") {Tmin=20;Tmax=31;}
      else if(W=="Ad") { Tmin=40; Tmax=46;}
      else if(W=="Vd") { Tmin=29; Tmax=39;}
      else if(W=="Au") { Tmin=26; Tmax=33;}
      else if(W=="Vu") { Tmin=36; Tmax=41;}
      else crash("time intervals for fits in channel: "+W+" not implemented");
    
    }
    else if(ixg==7) {
      if(W=="A") {Tmin= 18; Tmax=27;}
      else if(W=="V") {Tmin=18;Tmax=31;}
      else if(W=="Ad") { Tmin=38; Tmax=46;}
      else if(W=="Vd") { Tmin=29; Tmax=36;}
      else if(W=="Au") { Tmin=26; Tmax=33;}
      else if(W=="Vu") { Tmin=36; Tmax=41;}
      else crash("time intervals for fits in channel: "+W+" not implemented");
    
    }
    else if(ixg==8) {
      if(W=="A") {Tmin= 19; Tmax=33;}
      else if(W=="V") {Tmin=16;Tmax=31;}
      else if(W=="Ad") { Tmin=27; Tmax=45;}
      else if(W=="Vd") { Tmin=27; Tmax=36;}
      else if(W=="Au") { Tmin=26; Tmax=33;}
      else if(W=="Vu") { Tmin=36; Tmax=41;}
      else crash("time intervals for fits in channel: "+W+" not implemented");
    
    }
    else if(ixg==9) {
      if(W=="A") {Tmin= 17; Tmax=32;}
      else if(W=="V") {Tmin=16;Tmax=29;}
      else if(W=="Ad") { Tmin=23; Tmax=38;}
      else if(W=="Vd") { Tmin=28; Tmax=38;}
      else if(W=="Au") { Tmin=26; Tmax=33;}
      else if(W=="Vu") { Tmin=36; Tmax=41;}
      else crash("time intervals for fits in channel: "+W+" not implemented");
    
    }
    else if(ixg==10) {
      if(W=="A") {Tmin= 16; Tmax=29;}
      else if(W=="V") {Tmin=17;Tmax=31;}
      else if(W=="Ad") { Tmin=23; Tmax=38;}
      else if(W=="Vd") { Tmin=38; Tmax=51;}
      else if(W=="Au") { Tmin=26; Tmax=33;}
      else if(W=="Vu") { Tmin=36; Tmax=41;}
      else crash("time intervals for fits in channel: "+W+" not implemented");
      
    }
    else if(ixg==11) {
      if(W=="A") {Tmin= 20; Tmax=30;}
      else if(W=="V") {Tmin=20;Tmax=30;}
      else if(W=="Ad") { Tmin=24; Tmax=35;}
      else if(W=="Vd") { Tmin=23; Tmax=38;}
      else if(W=="Au") { Tmin=28; Tmax=41;}
      else if(W=="Vu") { Tmin=7; Tmax=42;}
      else crash("time intervals for fits in channel: "+W+" not implemented");

    }
    else if(ixg==12) {
      if(W=="A") {Tmin= 20; Tmax=30;}
      else if(W=="V") {Tmin=20;Tmax=30;}
      else if(W=="Ad") { Tmin=20; Tmax=35;}
      else if(W=="Vd") { Tmin=23; Tmax=38;}
      else if(W=="Au") { Tmin=28; Tmax=41;}
      else if(W=="Vu") { Tmin=7; Tmax=42;}
      else crash("time intervals for fits in channel: "+W+" not implemented");

    }
    else crash("ixg: "+to_string(ixg)+" does not have an established fit range");

  }


  //############ ENSEMBLE C80  ####################

  else  if(Ens == "cC211a.06.80") {

    if(ixg==0) {
      if(W=="A") {Tmin= 19; Tmax=29;}
      else if(W=="V") {Tmin=30;Tmax=47;}
      else if(W=="Ad") { Tmin=40; Tmax=54;}
      else if(W=="Vd") { Tmin=47; Tmax=56;}
      else if(W=="Au") { Tmin=26; Tmax=34;}
      else if(W=="Vu") { Tmin=36; Tmax=50;}
      else crash("time intervals for fits in channel: "+W+" not implemented");

    }
    else if(ixg==1) {
      if(W=="A") {Tmin= 22; Tmax=29;}
      else if(W=="V") {Tmin=27;Tmax=44;}
      else if(W=="Ad") { Tmin=40; Tmax=54;}
      else if(W=="Vd") { Tmin=47; Tmax=56;}
      else if(W=="Au") { Tmin=26; Tmax=34;}
      else if(W=="Vu") { Tmin=36; Tmax=50;}
      else crash("time intervals for fits in channel: "+W+" not implemented");

    }
    else if(ixg==2) {
      if(W=="A") {Tmin= 22; Tmax=29;}
      else if(W=="V") {Tmin=27;Tmax=39;}
      else if(W=="Ad") { Tmin=41; Tmax=55;}
      else if(W=="Vd") { Tmin=41; Tmax=54;}
      else if(W=="Au") { Tmin=26; Tmax=34;}
      else if(W=="Vu") { Tmin=38; Tmax=52;}
      else crash("time intervals for fits in channel: "+W+" not implemented");

    }
    else if(ixg==3) {
      if(W=="A") {Tmin= 22; Tmax=29;}
      else if(W=="V") {Tmin=26;Tmax=39;}
      else if(W=="Ad") { Tmin=41; Tmax=55;}
      else if(W=="Vd") { Tmin=41; Tmax=54;}
      else if(W=="Au") { Tmin=27; Tmax=34;}
      else if(W=="Vu") { Tmin=38; Tmax=52;}
      else crash("time intervals for fits in channel: "+W+" not implemented");
    
    }
    else if(ixg==4) {
      if(W=="A") {Tmin= 23; Tmax=30;}
      else if(W=="V") {Tmin=26;Tmax=39;}
      else if(W=="Ad") { Tmin=41; Tmax=55;}
      else if(W=="Vd") { Tmin=36; Tmax=53;}
      else if(W=="Au") { Tmin=26; Tmax=34;}
      else if(W=="Vu") { Tmin=38; Tmax=52;}
      else crash("time intervals for fits in channel: "+W+" not implemented");
    
    }
    else if(ixg==5) {
      if(W=="A") {Tmin= 22; Tmax=29;}
      else if(W=="V") {Tmin=25;Tmax=39;}
      else if(W=="Ad") { Tmin=41; Tmax=55;}
      else if(W=="Vd") { Tmin=36; Tmax=53;}
      else if(W=="Au") { Tmin=26; Tmax=34;}
      else if(W=="Vu") { Tmin=38; Tmax=52;}
      else crash("time intervals for fits in channel: "+W+" not implemented");
    
    }
    else if(ixg==6) {
      if(W=="A") {Tmin= 22; Tmax=30;}
      else if(W=="V") {Tmin=23;Tmax=37;}
      else if(W=="Ad") { Tmin=41; Tmax=55;}
      else if(W=="Vd") { Tmin=33; Tmax=53;}
      else if(W=="Au") { Tmin=26; Tmax=34;}
      else if(W=="Vu") { Tmin=38; Tmax=52;}
      else crash("time intervals for fits in channel: "+W+" not implemented");
    
    }
    else if(ixg==7) {
      if(W=="A") {Tmin= 22; Tmax=31;}
      else if(W=="V") {Tmin=21;Tmax=36;}
      else if(W=="Ad") { Tmin=39; Tmax=55;}
      else if(W=="Vd") { Tmin=33; Tmax=53;}
      else if(W=="Au") { Tmin=26; Tmax=34;}
      else if(W=="Vu") { Tmin=38; Tmax=52;}
      else crash("time intervals for fits in channel: "+W+" not implemented");
    
    }
    else if(ixg==8) {
      if(W=="A") {Tmin= 24; Tmax=39;}
      else if(W=="V") {Tmin=19;Tmax=37;}
      else if(W=="Ad") { Tmin=35; Tmax=53;}
      else if(W=="Vd") { Tmin=30; Tmax=53;}
      else if(W=="Au") { Tmin=26; Tmax=34;}
      else if(W=="Vu") { Tmin=38; Tmax=52;}
      else crash("time intervals for fits in channel: "+W+" not implemented");
    
    }
    else if(ixg==9) {
      if(W=="A") {Tmin= 18; Tmax=37;}
      else if(W=="V") {Tmin=20;Tmax=35;}
      else if(W=="Ad") { Tmin=26; Tmax=46;}
      else if(W=="Vd") { Tmin=30; Tmax=56;}
      else if(W=="Au") { Tmin=26; Tmax=34;}
      else if(W=="Vu") { Tmin=38; Tmax=52;}
      else crash("time intervals for fits in channel: "+W+" not implemented");
    
    }
    else if(ixg==10) {
      if(W=="A") {Tmin= 16; Tmax=31;}
      else if(W=="V") {Tmin=18;Tmax=34;}
      else if(W=="Ad") { Tmin=26; Tmax=46;}
      else if(W=="Vd") { Tmin=34; Tmax=62;}
      else if(W=="Au") { Tmin=26; Tmax=34;}
      else if(W=="Vu") { Tmin=38; Tmax=52;}
      else crash("time intervals for fits in channel: "+W+" not implemented");

    }

    else crash("ixg: "+to_string(ixg)+" does not have an established fit range");

  }

  else if(Ens=="cD211a.054.96") {

   


    if(ixg==0) {
      if(W=="A") {Tmin= 21; Tmax=28;}
      else if(W=="V") {Tmin=28;Tmax=55;}
      else if(W=="Ad") { Tmin=20; Tmax=32;}
      else if(W=="Vd") { Tmin=20; Tmax=32;}
      else if(W=="Au") { Tmin=20; Tmax=32;}
      else if(W=="Vu") { Tmin=20; Tmax=32;}
      else crash("time intervals for fits in channel: "+W+" not implemented");

    }
    else if(ixg==1) {
      if(W=="A") {Tmin= 26; Tmax=34;}
      else if(W=="V") {Tmin=32;Tmax=53;}
      else if(W=="Ad") { Tmin=51; Tmax=55;}
      else if(W=="Vd") { Tmin=61; Tmax=76;}
      else if(W=="Au") { Tmin=29; Tmax=35;}
      else if(W=="Vu") { Tmin=48; Tmax=59;}
      else crash("time intervals for fits in channel: "+W+" not implemented");

    }
    else if(ixg==2) {
      if(W=="A") {Tmin= 26; Tmax=34;}
      else if(W=="V") {Tmin=32;Tmax=47;}
      else if(W=="Ad") { Tmin=51; Tmax=55;}
      else if(W=="Vd") { Tmin=57; Tmax=68;}
      else if(W=="Au") { Tmin=29; Tmax=35;}
      else if(W=="Vu") { Tmin=48; Tmax=59;}
      else crash("time intervals for fits in channel: "+W+" not implemented");

    }
    else if(ixg==3) {
      if(W=="A") {Tmin= 26; Tmax=34;}
      else if(W=="V") {Tmin=32;Tmax=47;}
      else if(W=="Ad") { Tmin=51; Tmax=55;}
      else if(W=="Vd") { Tmin=57; Tmax=68;}
      else if(W=="Au") { Tmin=29; Tmax=35;}
      else if(W=="Vu") { Tmin=48; Tmax=59;}
      else crash("time intervals for fits in channel: "+W+" not implemented");
    
    }
    else if(ixg==4) {
      if(W=="A") {Tmin= 26; Tmax=34;}
      else if(W=="V") {Tmin=30;Tmax=47;}
      else if(W=="Ad") { Tmin=50; Tmax=55;}
      else if(W=="Vd") { Tmin=57; Tmax=68;}
      else if(W=="Au") { Tmin=29; Tmax=37;}
      else if(W=="Vu") { Tmin=48; Tmax=59;}
      else crash("time intervals for fits in channel: "+W+" not implemented");
    
    }
    else if(ixg==5) {
      if(W=="A") {Tmin= 26; Tmax=34;}
      else if(W=="V") {Tmin=28;Tmax=44;}
      else if(W=="Ad") { Tmin=50; Tmax=55;}
      else if(W=="Vd") { Tmin=57; Tmax=68;}
      else if(W=="Au") { Tmin=29; Tmax=37;}
      else if(W=="Vu") { Tmin=48; Tmax=59;}
      else crash("time intervals for fits in channel: "+W+" not implemented");
    
    }
    else if(ixg==6) {
      if(W=="A") {Tmin= 26; Tmax=34;}
      else if(W=="V") {Tmin=27;Tmax=45;}
      else if(W=="Ad") { Tmin=50; Tmax=55;}
      else if(W=="Vd") { Tmin=57; Tmax=68;}
      else if(W=="Au") { Tmin=29; Tmax=37;}
      else if(W=="Vu") { Tmin=48; Tmax=59;}
      else crash("time intervals for fits in channel: "+W+" not implemented");
    
    }
    else if(ixg==7) {
      if(W=="A") {Tmin= 25; Tmax=35;}
      else if(W=="V")  {Tmin=25;Tmax=44;}
      else if(W=="Ad") { Tmin=48; Tmax=54;}
      else if(W=="Vd") { Tmin=56; Tmax=68;}
      else if(W=="Au") { Tmin=29; Tmax=37;}
      else if(W=="Vu") { Tmin=49; Tmax=59;}
      else crash("time intervals for fits in channel: "+W+" not implemented");
    
    }
    else if(ixg==8) {
      if(W=="A") {Tmin= 27; Tmax=47;}
      else if(W=="V") {Tmin=24;Tmax=43;}
      else if(W=="Ad") { Tmin=47; Tmax=54;}
      else if(W=="Vd") { Tmin=52; Tmax=66;}
      else if(W=="Au") { Tmin=29; Tmax=37;}
      else if(W=="Vu") { Tmin=49; Tmax=59;}
      else crash("time intervals for fits in channel: "+W+" not implemented");
    
    }
    else if(ixg==9) {
      if(W=="A") {Tmin= 26; Tmax=47;}
      else if(W=="V") {Tmin=24;Tmax=42;}
      else if(W=="Ad") { Tmin=21; Tmax=50;}
      else if(W=="Vd") { Tmin=48; Tmax=65;}
      else if(W=="Au") { Tmin=29; Tmax=37;}
      else if(W=="Vu") { Tmin=49; Tmax=59;}
      else crash("time intervals for fits in channel: "+W+" not implemented");
    
    }
    else if(ixg==10) {
      if(W=="A") {Tmin= 22; Tmax=40;}
      else if(W=="V") {Tmin=24;Tmax=43;}
      else if(W=="Ad") { Tmin=26; Tmax=46;}
      else if(W=="Vd") { Tmin=48; Tmax=58;}
      else if(W=="Au") { Tmin=29; Tmax=37;}
      else if(W=="Vu") { Tmin=49; Tmax=59;}
      else crash("time intervals for fits in channel: "+W+" not implemented");

    }

    else crash("ixg: "+to_string(ixg)+" does not have an established fit range");


  }

  else crash("Ensemble: "+Ens+" does not have definite time intervals for FV and FA");


  

  return; 
   
  

 
}

void Compute_form_factors_Nissa() {


  //create directories
  boost::filesystem::create_directory("../data/ph_emission");
  boost::filesystem::create_directory("../data/ph_emission/"+ph_type);
  boost::filesystem::create_directory("../data/ph_emission/"+ph_type+"/"+Meson);
  boost::filesystem::create_directory("../data/ph_emission/"+ph_type+"/"+Meson+"/rate");
  boost::filesystem::create_directory("../data/ph_emission/"+ph_type+"/"+Meson+"/C");
  boost::filesystem::create_directory("../data/ph_emission/"+ph_type+"/"+Meson+"/H");
  boost::filesystem::create_directory("../data/ph_emission/"+ph_type+"/"+Meson+"/mass");
  boost::filesystem::create_directory("../data/ph_emission/"+ph_type+"/"+Meson+"/decay_const");
  boost::filesystem::create_directory("../data/ph_emission/"+ph_type+"/"+Meson+"/FF");
  boost::filesystem::create_directory("../data/ph_emission/"+ph_type+"/"+Meson+"/FF/continuum");
  boost::filesystem::create_directory("../data/ph_emission/"+ph_type+"/"+Meson+"/FF/continuum/BMA");
  boost::filesystem::create_directory("../data/ph_emission/"+ph_type+"/"+Meson+"/FF/continuum/BMA/uniform");
  boost::filesystem::create_directory("../data/ph_emission/"+ph_type+"/"+Meson+"/FF/continuum/BMA/AIC");
  boost::filesystem::create_directory("../data/ph_emission/"+ph_type+"/"+Meson+"/FF/continuum/BMA/random");
  boost::filesystem::create_directory("../data/ph_emission/"+ph_type+"/"+Meson+"/FF/per_kin");
  boost::filesystem::create_directory("../data/ph_emission/"+ph_type+"/"+Meson+"/FF_u");
  boost::filesystem::create_directory("../data/ph_emission/"+ph_type+"/"+Meson+"/FF_d");

  
  Get_xg_to_spline();
  Get_xg_to_spline_VMD();
  Get_lattice_spacings_to_print();

  GaussianMersenne GM(65622223);

  if(Enhanced_xg_only_B64) {
    rt_FF trash= Get_form_factors_Nissa(13,0,0,0,1,"", "_only_B64");
    cout<<"MODE: Enhanced_xg_only_B64"<<endl;
    cout<<"done!"<<endl;
    exit(-1);
  }


  vector<rt_FF> FF_ret_list;

  vector<bool> Use_three_finest_list({0,0,0,0,1,1});
  vector<bool> Include_a4_list({0,0,0,1,0,0});
  vector<bool> UseJack_list({1,0,1,1,0,1});
  vector<int>  num_xg_list({8,11,11,11,11,11});
  vector<double> Perform_continuum_extrapolation_list({0,1,1,1,1,1});
  vector<string> Corr_path_list({"_w_B96", "", "", "", "", ""});

  int N= UseJack_list.size();
  for(int i=0;i<N;i++) {
    xg_t_list.clear();
    Get_xg_t_list(num_xg_list[i]);
    string Fit_tag= ( (Include_a4_list[i]==true)?"wa4_":"");
    Fit_tag = ( (Use_three_finest_list[i]==true)?"wtf_":Fit_tag);
    FF_ret_list.push_back(Get_form_factors_Nissa(num_xg_list[i], Perform_continuum_extrapolation_list[i], Use_three_finest_list[i], Include_a4_list[i], UseJack_list[i], Fit_tag, Corr_path_list[i]));
  }


  
  

  //combine measurement with weight w_1,2
  vector<string> contribs({"FA", "FV", "FA_u", "FV_u", "FA_d", "FV_d"});

  
  
  xg_t_list.clear();
  int nxg=FF_ret_list[2].num_xg;
  Get_xg_t_list(nxg);

  distr_t FDs_star_distr(1);
  GaussianMersenne GM_FDs_star(424233);
  for(int ij=0;ij<NJ;ij++) FDs_star_distr.distr.push_back(fDs_star + GM_FDs_star()*fDs_star_err/sqrt(NJ -1.0));



  distr_t_list FV_pars_for_Br(1), FA_pars_for_Br(1);

  vector<boost::math::interpolators::cardinal_cubic_b_spline<double>> FV_spline_for_Br, FA_spline_for_Br;
  
  

  for( int i=0; i<(signed)contribs.size(); i++) {

   

    distr_t_list FF_jack(1), FF_boot(0);
    distr_t_list FF_jack_AIC(1), FF_boot_AIC(0);

    for(int ixg=1; ixg<nxg;ixg++) {

      rt_FF Fit_A_jack = FF_ret_list[2];
      rt_FF Fit_B_jack = FF_ret_list[5];
      rt_FF Fit_A_boot = FF_ret_list[1];
      rt_FF Fit_B_boot = FF_ret_list[4];

      double wA = 1.0;
      double wB = 1.0;
      double sum=wA+wB;
      wA /= sum;
      wB /= sum;

      double wA_AIC= exp(-0.5*(Fit_A_jack.Get_ch2(i)[ixg-1]*Fit_A_jack.Ndof + 2*Fit_A_jack.Npars - 2*Fit_A_jack.Nmeas));
      double wB_AIC =exp(-0.5*(Fit_B_jack.Get_ch2(i)[ixg-1]*Fit_B_jack.Ndof + 2*Fit_B_jack.Npars - 2*Fit_B_jack.Nmeas));
      double sum_AIC= wA_AIC+wB_AIC;
      wA_AIC /= sum_AIC;
      wB_AIC /= sum_AIC;

    
      

      distr_t D_FF_jack= wA*Fit_A_jack.Get_FF(i).distr_list[ixg-1] + wB*Fit_B_jack.Get_FF(i).distr_list[ixg-1];
      distr_t D_FF_boot= wA*Fit_A_boot.Get_FF(i).distr_list[ixg-1] + wB*Fit_B_boot.Get_FF(i).distr_list[ixg-1];

      distr_t D_FF_jack_AIC= wA_AIC*Fit_A_jack.Get_FF(i).distr_list[ixg-1] + wB_AIC*Fit_B_jack.Get_FF(i).distr_list[ixg-1];
      distr_t D_FF_boot_AIC= wA_AIC*Fit_A_boot.Get_FF(i).distr_list[ixg-1] + wB_AIC*Fit_B_boot.Get_FF(i).distr_list[ixg-1];

      

      double syst= wA*pow( Fit_A_jack.Get_FF(i).ave(ixg-1) - D_FF_jack.ave(),2) + wB*pow( Fit_B_jack.Get_FF(i).ave(ixg-1) - D_FF_jack.ave(),2);
      syst = sqrt(syst);
      distr_t syst_jack(1), syst_boot(0);

      double syst_AIC= wA_AIC*pow( Fit_A_jack.Get_FF(i).ave(ixg-1) - D_FF_jack_AIC.ave(),2) + wB_AIC*pow( Fit_B_jack.Get_FF(i).ave(ixg-1) - D_FF_jack_AIC.ave(),2);
      syst_AIC = sqrt(syst_AIC);
      distr_t syst_AIC_jack(1), syst_AIC_boot(0);

      for(int iboot=0;iboot<Nboots;iboot++) { double rn= GM();  syst_boot.distr.push_back( rn*syst); syst_AIC_boot.distr.push_back( rn*syst_AIC);}

      string contrib_new= contribs[i];

      if(contribs[i].size() > 3 ) { 
	if(contribs[i].substr(3,1) == "d") contrib_new= contrib_new.substr(0,3)+"strange";
	if(contribs[i].substr(3,1) == "u") contrib_new= contrib_new.substr(0,3)+"charm";
      }

      ofstream PrintRand("../data/ph_emission/"+ph_type+"/"+Meson+"/FF/continuum/BMA/random/"+contrib_new+"_xg"+to_string_with_precision(ixg*0.1,1)+".rand");
      for(int ijack=0;ijack<NJ;ijack++)  {
	double rn= GM();
	PrintRand<<rn<<endl;
	syst_jack.distr.push_back( rn*syst/sqrt(NJ-1.0)); syst_AIC_jack.distr.push_back( rn*syst_AIC/sqrt(NJ-1.0));
      }
      PrintRand.close();

      double stat_jack = D_FF_jack.err();
      double stat_jack_AIC= D_FF_jack_AIC.err();

      D_FF_jack = D_FF_jack + syst_jack;
      D_FF_boot = D_FF_boot + syst_boot;

      
      D_FF_jack_AIC = D_FF_jack_AIC + syst_AIC_jack;
      D_FF_boot_AIC = D_FF_boot_AIC + syst_AIC_boot;
      
      FF_jack.distr_list.push_back(  D_FF_jack );
      FF_boot.distr_list.push_back(  D_FF_boot );
      FF_jack_AIC.distr_list.push_back( D_FF_jack_AIC );
      FF_boot_AIC.distr_list.push_back( D_FF_boot_AIC );



      cout<<"#########################"<<endl;
      cout<<"Performing BMA for contrib: "<<contribs[i]<<" ixg: "<<0.1*ixg<<endl;
      cout<<"Four lattice spacings (Meas A): "<<endl;
      cout<<"w(U): "<<wA<<" w(AIC): "<<wA_AIC<<endl;
      cout<<"Nmeas: "<<Fit_A_jack.Nmeas<<endl;
      cout<<"Ndof: "<<Fit_A_jack.Ndof<<endl;
      cout<<"Npars: "<<Fit_A_jack.Npars<<endl;
      cout<<"ch2/dof: "<<Fit_A_jack.Get_ch2(i)[ixg-1]<<endl;
      cout<<"Three finest (Meas B): "<<endl;
      cout<<"w(U): "<<wB<<" w(AIC): "<<wB_AIC<<endl;
      cout<<"Nmeas: "<<Fit_B_jack.Nmeas<<endl;
      cout<<"Ndof: "<<Fit_B_jack.Ndof<<endl;
      cout<<"Npars: "<<Fit_B_jack.Npars<<endl;
      cout<<"ch2/dof: "<<Fit_B_jack.Get_ch2(i)[ixg-1]<<endl;
      cout<<"Meas A: "<<Fit_A_jack.Get_FF(i).ave(ixg-1)<<" +- "<<Fit_A_jack.Get_FF(i).err(ixg-1)<<endl;
      cout<<"Meas B: "<<Fit_B_jack.Get_FF(i).ave(ixg-1)<<" +- "<<Fit_B_jack.Get_FF(i).err(ixg-1)<<endl;
      
      cout<<"Combined(U): "<<D_FF_jack.ave()<<" +- "<<D_FF_jack.err()<<endl;
      cout<<"Expected error stat(U): "<< sqrt( wA*pow(Fit_A_jack.Get_FF(i).err(ixg-1),2) + wB*pow(Fit_B_jack.Get_FF(i).err(ixg-1),2))<<endl;
      cout<<"Expected error syst(U): "<< syst<<endl;
      cout<<"stat(U): "<<stat_jack<<endl;
      cout<<"stat+syst(U): "<< sqrt( syst*syst + wA*pow(Fit_A_jack.Get_FF(i).err(ixg-1),2) + wB*pow(Fit_B_jack.Get_FF(i).err(ixg-1),2))<<endl;
      
      cout<<"Combined(AIC): "<<D_FF_jack_AIC.ave()<<" +- "<<D_FF_jack_AIC.err()<<endl;
      cout<<"Expected error stat(AIC): "<< sqrt( wA_AIC*pow(Fit_A_jack.Get_FF(i).err(ixg-1),2) + wB_AIC*pow(Fit_B_jack.Get_FF(i).err(ixg-1),2))<<endl;
      cout<<"Expected error syst(AIC): "<< syst_AIC<<endl;
      cout<<"stat(AIC): "<<stat_jack_AIC<<endl;
      cout<<"stat+syst(AIC): "<< sqrt( syst_AIC*syst_AIC + wA_AIC*pow(Fit_A_jack.Get_FF(i).err(ixg-1),2) + wB_AIC*pow(Fit_B_jack.Get_FF(i).err(ixg-1),2))<<endl;
      cout<<"#########################"<<endl;
      
      
    }


    //print form factors

    Print_To_File({}, {xg_t_list, FF_jack.ave(), FF_jack.err()},  "../data/ph_emission/"+ph_type+"/"+Meson+"/FF/continuum/BMA/uniform/"+contribs[i]+"_jack.dat", "", "");
    Print_To_File({}, {xg_t_list, FF_boot.ave(), FF_boot.err()},  "../data/ph_emission/"+ph_type+"/"+Meson+"/FF/continuum/BMA/uniform/"+contribs[i]+"_boot.dat", "", "");
    
    Print_To_File({}, {xg_t_list, FF_jack_AIC.ave(), FF_jack_AIC.err()},  "../data/ph_emission/"+ph_type+"/"+Meson+"/FF/continuum/BMA/AIC/"+contribs[i]+"_jack.dat", "", "");
    Print_To_File({}, {xg_t_list, FF_boot_AIC.ave(), FF_boot_AIC.err()},  "../data/ph_emission/"+ph_type+"/"+Meson+"/FF/continuum/BMA/AIC/"+contribs[i]+"_boot.dat", "", "");


    //print covariance matrix

    Eigen::MatrixXd Cov_boot(nxg-1, nxg-1);
    Eigen::MatrixXd Cov_jack(nxg-1, nxg-1);
    Eigen::MatrixXd Cov_boot_AIC(nxg-1, nxg-1);
    Eigen::MatrixXd Cov_jack_AIC(nxg-1, nxg-1);

    Eigen::MatrixXd Corr_boot(nxg-1, nxg-1);
    Eigen::MatrixXd Corr_jack(nxg-1, nxg-1);
    Eigen::MatrixXd Corr_boot_AIC(nxg-1, nxg-1);
    Eigen::MatrixXd Corr_jack_AIC(nxg-1, nxg-1);
 
    for(int x_xg=1; x_xg<nxg;x_xg++) {
      for(int y_xg=1; y_xg<nxg;y_xg++) {

      
	Cov_boot(x_xg-1, y_xg-1) = FF_boot.distr_list[x_xg-1]%FF_boot.distr_list[y_xg-1];
	Cov_jack(x_xg-1, y_xg-1) = FF_jack.distr_list[x_xg-1]%FF_jack.distr_list[y_xg-1];
	Cov_boot_AIC(x_xg-1, y_xg-1) = FF_boot_AIC.distr_list[x_xg-1]%FF_boot_AIC.distr_list[y_xg-1];
	Cov_jack_AIC(x_xg-1, y_xg-1) = FF_jack_AIC.distr_list[x_xg-1]%FF_jack_AIC.distr_list[y_xg-1];

      
	Corr_boot(x_xg-1, y_xg-1) = Cov_boot(x_xg-1,y_xg-1)/(FF_boot.err(x_xg-1)*FF_boot.err(y_xg-1));
	Corr_jack(x_xg-1, y_xg-1) = Cov_jack(x_xg-1,y_xg-1)/(FF_jack.err(x_xg-1)*FF_jack.err(y_xg-1));
	Corr_boot_AIC(x_xg-1, y_xg-1) = Cov_boot_AIC(x_xg-1,y_xg-1)/(FF_boot_AIC.err(x_xg-1)*FF_boot_AIC.err(y_xg-1));
	Corr_jack_AIC(x_xg-1, y_xg-1) = Cov_jack_AIC(x_xg-1,y_xg-1)/(FF_jack_AIC.err(x_xg-1)*FF_jack_AIC.err(y_xg-1));
      }
    }

    //Print To File
  
    ofstream Print_Cov_boot("../data/ph_emission/"+ph_type+"/"+Meson+"/FF/continuum/BMA/uniform/"+contribs[i]+"_boot.cov");
    ofstream Print_Cov_jack("../data/ph_emission/"+ph_type+"/"+Meson+"/FF/continuum/BMA/uniform/"+contribs[i]+"_jack.cov");
    ofstream Print_Cov_boot_AIC("../data/ph_emission/"+ph_type+"/"+Meson+"/FF/continuum/BMA/AIC/"+contribs[i]+"_boot.cov");
    ofstream Print_Cov_jack_AIC("../data/ph_emission/"+ph_type+"/"+Meson+"/FF/continuum/BMA/AIC/"+contribs[i]+"_jack.cov");

    ofstream Print_Corr_boot("../data/ph_emission/"+ph_type+"/"+Meson+"/FF/continuum/BMA/uniform/"+contribs[i]+"_boot.corr");
    ofstream Print_Corr_jack("../data/ph_emission/"+ph_type+"/"+Meson+"/FF/continuum/BMA/uniform/"+contribs[i]+"_jack.corr");
    ofstream Print_Corr_boot_AIC("../data/ph_emission/"+ph_type+"/"+Meson+"/FF/continuum/BMA/AIC/"+contribs[i]+"_boot.corr");
    ofstream Print_Corr_jack_AIC("../data/ph_emission/"+ph_type+"/"+Meson+"/FF/continuum/BMA/AIC/"+contribs[i]+"_jack.corr");



    Print_Cov_boot<<Cov_boot<<endl;
    Print_Cov_boot.close();
    Print_Cov_jack<<Cov_jack<<endl;
    Print_Cov_jack.close();
    Print_Cov_boot_AIC<<Cov_boot_AIC<<endl;
    Print_Cov_boot_AIC.close();
    Print_Cov_jack_AIC<<Cov_jack_AIC<<endl;
    Print_Cov_jack_AIC.close();


  
    Print_Corr_boot<<Corr_boot<<endl;
    Print_Corr_boot.close();
    Print_Corr_jack<<Corr_jack<<endl;
    Print_Corr_jack.close();
    Print_Corr_boot_AIC<<Corr_boot_AIC<<endl;
    Print_Corr_boot_AIC.close();
    Print_Corr_jack_AIC<<Corr_jack_AIC<<endl;
    Print_Corr_jack_AIC.close();


    
    //interpolate the form factors

    vector<boost::math::interpolators::cardinal_cubic_b_spline<double>> FF_cont_interpol_jacks;
    vector<boost::math::interpolators::cardinal_cubic_b_spline<double>> FF_AIC_cont_interpol_jacks;

    //interpolate continuum extrapolated form factors
    for(int ijack=0;ijack<NJ;ijack++) {

      Vfloat FF_ij, FF_AIC_ij;
  
      for(int ixg=1;ixg<nxg;ixg++) {
	FF_ij.push_back( FF_jack.distr_list[ixg-1].distr[ijack]);
	FF_AIC_ij.push_back( FF_jack_AIC.distr_list[ixg-1].distr[ijack]);
      }

      FF_cont_interpol_jacks.emplace_back( FF_ij.begin(), FF_ij.end(), 0.1, 0.1);
      FF_AIC_cont_interpol_jacks.emplace_back( FF_AIC_ij.begin(), FF_AIC_ij.end(), 0.1, 0.1);

    }

    //push_back if FA or FV
    if(contribs[i] == "FA") FA_spline_for_Br = FF_AIC_cont_interpol_jacks;
    if(contribs[i] == "FV") FV_spline_for_Br = FF_AIC_cont_interpol_jacks;


    auto FF_cont_interpol_distr= [&FF_cont_interpol_jacks](double xg) -> distr_t {
      distr_t return_distr(1);
      for(int ijack=0; ijack<NJ;ijack++) { return_distr.distr.push_back( FF_cont_interpol_jacks[ijack](xg));}
      return return_distr;
    };

    auto FF_AIC_cont_interpol_distr= [&FF_AIC_cont_interpol_jacks](double xg) -> distr_t {
      distr_t return_distr(1);
      for(int ijack=0; ijack<NJ;ijack++) { return_distr.distr.push_back( FF_AIC_cont_interpol_jacks[ijack](xg));}
      return return_distr;
    };
  
 
    distr_t_list FF_interpol_cont_to_print_distr(1);
    distr_t_list FF_AIC_interpol_cont_to_print_distr(1);
  
  
    for(auto &X: xg_to_spline) {
      FF_interpol_cont_to_print_distr.distr_list.push_back( FF_cont_interpol_distr(X));
      FF_AIC_interpol_cont_to_print_distr.distr_list.push_back( FF_AIC_cont_interpol_distr(X));
    }

    Print_To_File({}, {xg_to_spline, FF_interpol_cont_to_print_distr.ave(), FF_interpol_cont_to_print_distr.err()}, "../data/ph_emission/"+ph_type+"/"+Meson+"/FF/continuum/BMA/uniform/"+contribs[i]+"_interpol.dat", "", "#xg "+contribs[i]+" "+contribs[i]+"_err");

    Print_To_File({}, {xg_to_spline, FF_AIC_interpol_cont_to_print_distr.ave(), FF_AIC_interpol_cont_to_print_distr.err()}, "../data/ph_emission/"+ph_type+"/"+Meson+"/FF/continuum/BMA/AIC/"+contribs[i]+"_interpol.dat", "", "#xg "+contribs[i]+" "+contribs[i]+"_err");


    //#############################################################################


    vector<string> VMD_fit_types({"single_pole", "double_pole", "pure_VMD", "VMD_background", "double_pole_background", "VMD_background_linear", "single_pole_background"});
    vector<bool> cov_matrix_mode({true, false});
    vector<bool> Full_ranges({true,false});

    distr_t_list g_list(1), R1_list(1), R2_list(1), B_list(1), B_xg_list(1);
    vector<double> ch2_list;
    vector<double> Npars_list;
    vector<string> FF_t;
    vector<double> cmod;
    vector<double> full_mode;
    vector<double> det_boot_AIC_list;
    vector<double> det_jack_AIC_list;

  
    for( auto VMD_fit_type:VMD_fit_types) 
      for(auto Add_cov_matrix: cov_matrix_mode )
	for(auto full_range: Full_ranges) {
      
  
  

	  //###########################################################
	  //###########################################################
	  //###########################################################
	  //###########################################################
	  //###########################################################
	  //################                          #################
	  //################                          #################
	  //################      POLE-FITS           #################
	  //################                          #################
	  //################                          #################
	  //###########################################################
	  //###########################################################
	  //###########################################################
	  //###########################################################
	  //###########################################################

	  class ipar_VMD {
	  public:
	    ipar_VMD() : FF(0.0), FF_err(0.0) {}
	    double FF, FF_err, xg;
	  };
  
	  class fpar_VMD {
	  public:
	    fpar_VMD() {}
	    fpar_VMD(const Vfloat &par) {
	      if((signed)par.size() != 6) crash("In class fpar_VMD  class constructor Vfloat par has size != 5");
	      A=par[0];
	      M=par[1];
	      A2=par[2];
	      M2=par[3];
	      B=par[4];
	      Bxg=par[5];
	    }
	    double A, M, A2, M2,B, Bxg;
	  };

  
	  //init bootstrap fit
	  bootstrap_fit<fpar_VMD,ipar_VMD> bf_VMD(NJ);
	  bf_VMD.set_warmup_lev(0); //sets warmup
	  bf_VMD.Set_number_of_measurements(xg_t_list.size());
	  bf_VMD.Set_verbosity(1);
 
  
	  //################################################
	  //define initial condition for params
	  double A_guess, M_guess, M2_guess;
	  if( contribs[i].substr(1,1) == "V") { M_guess=1.073; M2_guess= 1.38 ; }
	  else if( contribs[i].substr(1,1) == "A") { M_guess=1.25; M2_guess= 1.29 ;}
	  else crash("FF: "+contribs[i]+" is not V or A");
	  A_guess= FF_jack_AIC.ave(0)*M_guess*(M_guess-1.0);
	  //###############################################

  
	  bf_VMD.Add_par("A", A_guess, A_guess/10);
	  bf_VMD.Add_par("M", M_guess, M_guess/10);
	  bf_VMD.Add_par("A2", A_guess, A_guess/10);
	  bf_VMD.Add_par("M2", M2_guess, M2_guess/10);
	  bf_VMD.Add_par("B", A_guess, A_guess/10);
	  bf_VMD.Add_par("Bxg", A_guess, A_guess/10);
	  //fit on mean values to get ch2
	  bootstrap_fit<fpar_VMD,ipar_VMD> bf_VMD_ch2(1);
	  bf_VMD_ch2.set_warmup_lev(0); //sets warmup
	  bf_VMD_ch2.Set_number_of_measurements(xg_t_list.size());
	  bf_VMD_ch2.Set_verbosity(1);
	  bf_VMD_ch2.Add_par("A", A_guess, A_guess/10);
	  bf_VMD_ch2.Add_par("M", M_guess, M_guess/10);
	  bf_VMD_ch2.Add_par("A2", A_guess, A_guess/10);
	  bf_VMD_ch2.Add_par("M2", M2_guess, M2_guess/10);
	  bf_VMD_ch2.Add_par("B", A_guess, A_guess/10);
	  bf_VMD_ch2.Add_par("Bxg", A_guess, A_guess/10);
	  int Npars=0;

	  if(Add_cov_matrix) bf_VMD.Add_covariance_matrix(Cov_boot_AIC);
	  if(Add_cov_matrix) bf_VMD_ch2.Add_covariance_matrix(Cov_boot_AIC);


	  double det_boot_AIC=1;
	  double det_jack_AIC=1;
	  if(Add_cov_matrix) {
	    det_boot_AIC= Corr_boot_AIC.determinant();
	    det_jack_AIC= Corr_jack_AIC.determinant();
	  }
  
	  if(VMD_fit_type =="pure_VMD") {
    
	    bf_VMD.Fix_par("M", M_guess);
	    bf_VMD.Fix_par("M2", M2_guess);
	    bf_VMD.Fix_par("A2", 0.0);
	    bf_VMD.Fix_par("B", 0.0);
	    bf_VMD.Fix_par("Bxg",0.0);

	    bf_VMD_ch2.Fix_par("M", M_guess);
	    bf_VMD_ch2.Fix_par("M2", M2_guess);
	    bf_VMD_ch2.Fix_par("A2", 0.0);
	    bf_VMD_ch2.Fix_par("B", 0.0);
	    bf_VMD_ch2.Fix_par("Bxg",0.0);

	    Npars=1;
	  }
	  else if(VMD_fit_type == "single_pole") {

	    bf_VMD.Fix_par("M2", M2_guess);
	    bf_VMD.Fix_par("A2", 0.0);
	    bf_VMD.Fix_par("B", 0.0);
	    bf_VMD.Fix_par("Bxg",0.0);

	    bf_VMD_ch2.Fix_par("M2", M2_guess);
	    bf_VMD_ch2.Fix_par("A2", 0.0);
	    bf_VMD_ch2.Fix_par("B", 0.0);
	    bf_VMD_ch2.Fix_par("Bxg",0.0);

	    Npars=2;
	  }
	  else if(VMD_fit_type == "double_pole") {

	    bf_VMD.Fix_par("M", M_guess);
	    bf_VMD.Fix_par("M2", M2_guess);
	    bf_VMD.Fix_par("B", 0.0);
	    bf_VMD.Fix_par("Bxg",0.0);
   
	    bf_VMD_ch2.Fix_par("M", M_guess);
	    bf_VMD_ch2.Fix_par("M2", M2_guess);
	    bf_VMD_ch2.Fix_par("B", 0.0);
	    bf_VMD_ch2.Fix_par("Bxg",0.0);
    
	    Npars=2;
	  }
	  else if(VMD_fit_type == "VMD_background") {

	    bf_VMD.Fix_par("M", M_guess);
	    bf_VMD.Fix_par("M2", M2_guess);
	    bf_VMD.Fix_par("A2", 0.0);
	    bf_VMD.Fix_par("Bxg",0.0);
    
	    bf_VMD_ch2.Fix_par("M", M_guess);
	    bf_VMD_ch2.Fix_par("M2", M2_guess);
	    bf_VMD_ch2.Fix_par("A2", 0.0);
	    bf_VMD_ch2.Fix_par("Bxg",0.0);
    
	    Npars=2;

	  }
	  else if(VMD_fit_type == "double_pole_background") {

	    bf_VMD.Fix_par("M", M_guess);
	    bf_VMD.Fix_par("M2", M2_guess);
	    bf_VMD.Fix_par("Bxg",0.0);
    
	    bf_VMD_ch2.Fix_par("M", M_guess);
	    bf_VMD_ch2.Fix_par("M2", M2_guess);
	    bf_VMD_ch2.Fix_par("Bxg",0.0);
       
	    Npars=3;

	  }
	  else if(VMD_fit_type == "VMD_background_linear") {

	    bf_VMD.Fix_par("M", M_guess);
	    bf_VMD.Fix_par("M2", M2_guess);
	    bf_VMD.Fix_par("A2",0.0);


	    bf_VMD_ch2.Fix_par("M", M_guess);
	    bf_VMD_ch2.Fix_par("M2", M2_guess);
	    bf_VMD_ch2.Fix_par("A2",0.0);

	    Npars=3;


	  }

	  else if(VMD_fit_type== "single_pole_background") {

	    bf_VMD.Fix_par("M2", M2_guess);
	    bf_VMD.Fix_par("A2", 0.0);
	    bf_VMD.Fix_par("Bxg", 0.0);

	    bf_VMD_ch2.Fix_par("M2", M2_guess);
	    bf_VMD_ch2.Fix_par("A2",0.0);
	    bf_VMD_ch2.Fix_par("Bxg",0.0);

	    Npars=3;

	  }
	  else crash("Fit type: "+VMD_fit_type+" not yet implemented");
  
 
	  //ansatz
	  bf_VMD.ansatz=  [&full_range](const fpar_VMD &p, const ipar_VMD &ip) {

	    if( (full_range==false) && (ip.xg > 0.5)) return 0.0;

	    double E_res= sqrt( p.M*p.M + ip.xg*ip.xg/4.0);
	    double E_res_2= sqrt( p.M2*p.M2 + ip.xg*ip.xg/4.0);
	    return  p.A/( E_res*( E_res - (1.0 - ip.xg/2.0)) )  + p.A2/(E_res_2*( E_res_2 - (1.0 - ip.xg/2.0)) ) + p.B + p.Bxg*ip.xg  ;
	  };
	  bf_VMD.measurement=  [&full_range](const fpar_VMD &p, const ipar_VMD &ip) {
	    if( (full_range==false) && (ip.xg > 0.5)) return 0.0;
	    return ip.FF;
	  };
	  bf_VMD.error=  [ ](const fpar_VMD &p, const ipar_VMD &ip) {

	    return ip.FF_err;
	  };

	  bf_VMD_ch2.ansatz= bf_VMD.ansatz;
	  bf_VMD_ch2.measurement = bf_VMD.measurement;
	  bf_VMD_ch2.error = bf_VMD.error;


	  //start fitting FA
	  //fill the data
	  vector<vector<ipar_VMD>> data_VMD(NJ);
	  vector<vector<ipar_VMD>> data_VMD_ch2(1);
	  //allocate space for output result
	  boot_fit_data<fpar_VMD> Bt_fit_VMD;
	  boot_fit_data<fpar_VMD> Bt_fit_VMD_ch2;
	  for(auto &data_iboot: data_VMD) data_iboot.resize(xg_t_list.size());
	  for(auto &data_iboot: data_VMD_ch2) data_iboot.resize(xg_t_list.size());
	  for(int ijack=0;ijack<NJ;ijack++) {
	    for(int ix=0;ix<(signed)xg_t_list.size();ix++) {
	      data_VMD[ijack][ix].FF = FF_jack_AIC.distr_list[ix].distr[ijack];
	      data_VMD[ijack][ix].FF_err= FF_jack_AIC.err(ix);
	      data_VMD[ijack][ix].xg= xg_t_list[ix];
	      if(ijack==0) {
		data_VMD_ch2[ijack][ix].FF = FF_jack_AIC.ave(ix);
		data_VMD_ch2[ijack][ix].FF_err= FF_jack_AIC.err(ix);
		data_VMD_ch2[ijack][ix].xg= xg_t_list[ix];

	      }
	    }
	  }

	  //append
	  bf_VMD.Append_to_input_par(data_VMD);
	  bf_VMD_ch2.Append_to_input_par(data_VMD_ch2);
	  //fit
	  cout<<"Fitting "<<contribs[i]<<" using VMD ansatz"<<endl;
	  Bt_fit_VMD= bf_VMD.Perform_bootstrap_fit();
	  Bt_fit_VMD_ch2= bf_VMD_ch2.Perform_bootstrap_fit();
	  int num_eff_xg= (full_range==true)?xg_t_list.size():5;
	  double ch2_red_VMD= Bt_fit_VMD_ch2.get_ch2_ave()/( 1.0*num_eff_xg -1.0*Npars);

	  //retrieve params
	  distr_t Ampl_F(1), pole_F(1), Ampl_F2(1), pole_F2(1), Bg_F(1), Bg_xg_F(1);
	  for(int ijack=0;ijack<NJ;ijack++) {
	    Ampl_F.distr.push_back( Bt_fit_VMD.par[ijack].A); pole_F.distr.push_back( Bt_fit_VMD.par[ijack].M);
	    Ampl_F2.distr.push_back( Bt_fit_VMD.par[ijack].A2); pole_F2.distr.push_back( Bt_fit_VMD.par[ijack].M2);
	    Bg_F.distr.push_back( Bt_fit_VMD.par[ijack].B);
	    Bg_xg_F.distr.push_back( Bt_fit_VMD.par[ijack].Bxg);

	  }


	  distr_t_list F_VMD_fit(1);
	  //plot fit function
	  for(auto &X: xg_to_spline_VMD) {
	    ipar_VMD pp_VMD;
	    pp_VMD.xg=X;
	    distr_t F_VMD_xg(1);
	    for(int ijack=0;ijack<NJ;ijack++) {
	      F_VMD_xg.distr.push_back( bf_VMD.ansatz( Bt_fit_VMD.par[ijack], pp_VMD));
	    }

	    F_VMD_fit.distr_list.push_back( F_VMD_xg);
	  }

	  //print
	  distr_t g= -2*Ampl_F/(pole_F*FDs_star_distr); 
	  string header_F= "#g: "+to_string_with_precision(g.ave(),5)+" +- "+to_string_with_precision(g.err(),5)+"  Ampl1: "+to_string_with_precision(Ampl_F.ave(),5)+" +- "+to_string_with_precision(Ampl_F.err(),5)+" M^res1/Mp: "+to_string_with_precision(pole_F.ave(), 5)+" +- "+to_string_with_precision(pole_F.err(), 5)+" Ampl2: "+to_string_with_precision(Ampl_F2.ave(),5)+" +- "+to_string_with_precision(Ampl_F2.err(),5)+" M^res2/Mp: "+to_string_with_precision(pole_F2.ave(), 5)+" +- "+to_string_with_precision(pole_F2.err(), 5)+" B: "+to_string_with_precision(Bg_F.ave(),5)+" +- "+to_string_with_precision(Bg_F.err(),5)+"Bxg: "+to_string_with_precision(Bg_xg_F.ave(),5)+" +- "+to_string_with_precision(Bg_xg_F.err(),5)+"   ch2/dof: "+to_string_with_precision(ch2_red_VMD  ,5)+" dof: "+to_string(xg_t_list.size()-Npars);

	  string red="";
	  if(full_range == false) red="_xg_upto_0.5";
  
	  Print_To_File({},{ xg_to_spline_VMD, F_VMD_fit.ave(), F_VMD_fit.err()} , "../data/ph_emission/"+ph_type+"/"+Meson+"/FF/continuum/BMA/AIC/"+contribs[i]+"_"+VMD_fit_type+"_corr_"+to_string(Add_cov_matrix)+red+".fit" , "", header_F);

	  
	  //push back parameter of the chosen parameterization to be used for Br_fit
	  if(full_range && Add_cov_matrix && (VMD_fit_type== parameterization_to_be_used_in_Br)) {
	  
	  if(contribs[i] == "FA") {

	    FA_pars_for_Br.distr_list.push_back(Ampl_F);
	    FA_pars_for_Br.distr_list.push_back(pole_F);
	    FA_pars_for_Br.distr_list.push_back(Ampl_F2);
	    FA_pars_for_Br.distr_list.push_back(pole_F2);
	    FA_pars_for_Br.distr_list.push_back(Bg_F);
	    FA_pars_for_Br.distr_list.push_back(Bg_xg_F);
	  }

	  else if(contribs[i] == "FV") {

	    FV_pars_for_Br.distr_list.push_back(Ampl_F);
	    FV_pars_for_Br.distr_list.push_back(pole_F);
	    FV_pars_for_Br.distr_list.push_back(Ampl_F2);
	    FV_pars_for_Br.distr_list.push_back(pole_F2);
	    FV_pars_for_Br.distr_list.push_back(Bg_F);
	    FV_pars_for_Br.distr_list.push_back(Bg_xg_F);
	    


	  }
	  }

	  g_list.distr_list.push_back(g);
	  ch2_list.push_back(ch2_red_VMD);
	  R1_list.distr_list.push_back( pole_F);
	  R2_list.distr_list.push_back( pole_F2);
	  B_list.distr_list.push_back(Bg_F);
	  B_xg_list.distr_list.push_back(Bg_xg_F);
	  cmod.push_back(Add_cov_matrix);
	  full_mode.push_back(full_range);
	  FF_t.push_back( VMD_fit_type);
	  Npars_list.push_back(Npars);
	  det_boot_AIC_list.push_back( det_boot_AIC);
	  det_jack_AIC_list.push_back( det_jack_AIC);
  

   
	  //###########################################################
	  //###########################################################
	  //###########################################################
	  //###########################################################
	  //###########################################################
	  //###########################################################




	  //##################################################################
  
	}
  
  

    //Print_info to file
    Print_To_File({FF_t}, { g_list.ave(), g_list.err(), R1_list.ave(), R1_list.err(), B_list.ave(), B_list.err(), B_xg_list.ave(), B_xg_list.err(), ch2_list, Npars_list, cmod, full_mode, det_boot_AIC_list, det_jack_AIC_list}, "../data/ph_emission/"+ph_type+"/"+Meson+"/FF/continuum/BMA/AIC/"+contribs[i]+"_VMD.info", "", "#fit_type  g  R1 B  Bxg   ch2/dof  Npars  Cov_added Full_range? det_corr_boot  det_corr_jack");

    
  }

  //Perform analysis of the branching ratio

  GaussianMersenne G_FLAGS(8762);

  distr_t fDs_FLAG_distr(1);
  distr_t Vcs_FLAG_distr(1);
  for(int ijack=0;ijack<NJ;ijack++) fDs_FLAG_distr.distr.push_back( fDs_FLAG + fDs_FLAG_err*G_FLAGS()/sqrt(NJ-1.0));
  for(int ijack=0;ijack<NJ;ijack++) Vcs_FLAG_distr.distr.push_back( Vcs_FLAG + dVcs_FLAG*G_FLAGS()/sqrt(NJ-1.0));

  Vfloat xg_cut({2*0.01/1.96835});
  Vfloat xg_diff;
  int Nxg_diff=101;
  int Nxg_cuts=41;
  for(int ixg_cut=1;ixg_cut<(Nxg_cuts-1);ixg_cut++) xg_cut.push_back( ixg_cut*1.0/(Nxg_cuts-1)); //from 0.025 up to ~0.9 in xg
  for(int ixg_diff=0;ixg_diff<Nxg_diff;ixg_diff++) xg_diff.push_back( ixg_diff*1.0/(Nxg_diff-1));
 
  distr_t_list Br_TOT(true, xg_cut.size()), Br_PT(true, xg_cut.size() ), Br_INT(true, xg_cut.size()), Br_SD(true, xg_cut.size());
  distr_t_list Br_diff_TOT(true, xg_diff.size()), Br_diff_PT(true, xg_diff.size() ), Br_diff_INT(true, xg_diff.size()), Br_diff_SD(true, xg_diff.size());

  distr_t_list Br_spline_TOT(true, xg_cut.size()), Br_spline_PT(true, xg_cut.size() ), Br_spline_INT(true, xg_cut.size()), Br_spline_SD(true, xg_cut.size());
  distr_t Br_spline_muon_TOT(true), Br_spline_tau_TOT(true);
  distr_t Br_spline_muon_PT(true), Br_spline_tau_PT(true);
  distr_t_list Br_diff_spline_TOT(true, xg_diff.size()), Br_diff_spline_PT(true, xg_diff.size() ), Br_diff_spline_INT(true, xg_diff.size()), Br_diff_spline_SD(true, xg_diff.size());
  
  for(int ixg_cut=0; ixg_cut<(signed)xg_cut.size(); ixg_cut++) {
    for(int ijack=0; ijack<NJ;ijack++) {

      double FA_A= FA_pars_for_Br.distr_list[0].distr[ijack];
      double FA_M= FA_pars_for_Br.distr_list[1].distr[ijack];
      double FA_A2= FA_pars_for_Br.distr_list[2].distr[ijack];
      double FA_M2= FA_pars_for_Br.distr_list[3].distr[ijack];
      double FA_B = FA_pars_for_Br.distr_list[4].distr[ijack];
      double FA_Bxg = FA_pars_for_Br.distr_list[5].distr[ijack];
      
      auto FA_Br = [&FA_A, &FA_M, &FA_A2, &FA_M2, &FA_B, &FA_Bxg](double x) {

	
	double E_res= sqrt( FA_M*FA_M + x*x/4.0);
	double E_res_2 = sqrt( FA_M2*FA_M2 + x*x/4.0);
	return  FA_A/( E_res*( E_res - (1.0 - x/2.0)) )  + FA_A2/(E_res_2*( E_res_2 - (1.0 - x/2.0)) ) + FA_B + FA_Bxg*x  ;
	
      };


      double FV_A= FV_pars_for_Br.distr_list[0].distr[ijack];
      double FV_M= FV_pars_for_Br.distr_list[1].distr[ijack];
      double FV_A2= FV_pars_for_Br.distr_list[2].distr[ijack];
      double FV_M2= FV_pars_for_Br.distr_list[3].distr[ijack];
      double FV_B = FV_pars_for_Br.distr_list[4].distr[ijack];
      double FV_Bxg = FV_pars_for_Br.distr_list[5].distr[ijack];
      
      auto FV_Br = [&FV_A, &FV_M, &FV_A2, &FV_M2, &FV_B, &FV_Bxg](double x) {

	
	double E_res= sqrt( FV_M*FV_M + x*x/4.0);
	double E_res_2= sqrt( FV_M2*FV_M2 + x*x/4.0);
	return  FV_A/( E_res*( E_res - (1.0 - x/2.0)) )  + FV_A2/(E_res_2*( E_res_2 - (1.0 - x/2.0)) ) + FV_B + FV_Bxg*x  ;
	
      };

      auto FA_Br_spline = [&FA_spline_for_Br, &ijack](double x) { return FA_spline_for_Br[ijack](x);};
      auto FV_Br_spline = [&FV_spline_for_Br, &ijack](double x) { return FV_spline_for_Br[ijack](x);};

      double fp=fDs_FLAG_distr.distr[ijack];
      double Vcs_2 = pow(Vcs_FLAG_distr.distr[ijack],2);
      //compute Br
     
      Br_PT.distr_list[ixg_cut].distr.push_back( Vcs_2*Compute_Ds_lnugamma_decay_rate(re,"PT",ijack, NJ,  FV_Br, FA_Br, fp,  xg_cut[ixg_cut], 1-re*re ));
      Br_INT.distr_list[ixg_cut].distr.push_back( Vcs_2*Compute_Ds_lnugamma_decay_rate(re,"INT",ijack, NJ,  FV_Br, FA_Br, fp,  xg_cut[ixg_cut], 1-re*re ));
      Br_SD.distr_list[ixg_cut].distr.push_back( Vcs_2*Compute_Ds_lnugamma_decay_rate(re,"SD",ijack, NJ,  FV_Br, FA_Br, fp,  xg_cut[ixg_cut], 1-re*re ));
      Br_TOT.distr_list[ixg_cut].distr.push_back( Vcs_2*Compute_Ds_lnugamma_decay_rate(re,"TOT",ijack, NJ,  FV_Br, FA_Br, fp,  xg_cut[ixg_cut], 1-re*re ));
     
  
      Br_spline_PT.distr_list[ixg_cut].distr.push_back( Vcs_2*Compute_Ds_lnugamma_decay_rate(re,"PT",ijack, NJ,  FV_Br_spline, FA_Br_spline, fp,  xg_cut[ixg_cut], 1-re*re ));
      Br_spline_INT.distr_list[ixg_cut].distr.push_back( Vcs_2*Compute_Ds_lnugamma_decay_rate(re,"INT",ijack, NJ,  FV_Br_spline, FA_Br_spline, fp,  xg_cut[ixg_cut], 1-re*re ));
      Br_spline_SD.distr_list[ixg_cut].distr.push_back( Vcs_2*Compute_Ds_lnugamma_decay_rate(re,"SD",ijack, NJ,  FV_Br_spline, FA_Br_spline, fp,  xg_cut[ixg_cut], 1-re*re ));
      Br_spline_TOT.distr_list[ixg_cut].distr.push_back( Vcs_2*Compute_Ds_lnugamma_decay_rate(re, "TOT",ijack, NJ,  FV_Br_spline, FA_Br_spline, fp,  xg_cut[ixg_cut], 1-re*re ));

      if(ixg_cut==0) {
	Br_spline_muon_PT.distr.push_back( Vcs_2*Compute_Ds_lnugamma_decay_rate(rm, "PT",ijack, NJ,  FV_Br_spline, FA_Br_spline, fp,  xg_cut[ixg_cut], 1-rm*rm )   );
	Br_spline_muon_TOT.distr.push_back( Vcs_2*Compute_Ds_lnugamma_decay_rate(rm, "TOT",ijack, NJ,  FV_Br_spline, FA_Br_spline, fp,  xg_cut[ixg_cut], 1-rm*rm ));
	Br_spline_tau_TOT.distr.push_back( Vcs_2*Compute_Ds_lnugamma_decay_rate(rt, "TOT",ijack, NJ,  FV_Br_spline, FA_Br_spline, fp,  xg_cut[ixg_cut], 1-rt*rt ));
	Br_spline_tau_PT.distr.push_back( Vcs_2*Compute_Ds_lnugamma_decay_rate(rt, "PT",ijack, NJ,  FV_Br_spline, FA_Br_spline, fp,  xg_cut[ixg_cut] , 1-rt*rt));
      }

  	  
    }
  }

  for(int ixg_d=0;ixg_d<(signed)xg_diff.size();ixg_d++) {
     for(int ijack=0; ijack<NJ;ijack++) {

       double FA_A= FA_pars_for_Br.distr_list[0].distr[ijack];
       double FA_M= FA_pars_for_Br.distr_list[1].distr[ijack];
       double FA_A2= FA_pars_for_Br.distr_list[2].distr[ijack];
       double FA_M2= FA_pars_for_Br.distr_list[3].distr[ijack];
       double FA_B = FA_pars_for_Br.distr_list[4].distr[ijack];
       double FA_Bxg = FA_pars_for_Br.distr_list[5].distr[ijack];
      
       auto FA_Br = [&FA_A, &FA_M, &FA_A2, &FA_M2, &FA_B, &FA_Bxg](double x) {

	 
	 double E_res= sqrt( FA_M*FA_M + x*x/4.0);
	 double E_res_2 = sqrt( FA_M2*FA_M2 + x*x/4.0);
	 return  FA_A/( E_res*( E_res - (1.0 - x/2.0)) )  + FA_A2/(E_res_2*( E_res_2 - (1.0 - x/2.0)) ) + FA_B + FA_Bxg*x  ;
	
       };


       double FV_A= FV_pars_for_Br.distr_list[0].distr[ijack];
       double FV_M= FV_pars_for_Br.distr_list[1].distr[ijack];
       double FV_A2= FV_pars_for_Br.distr_list[2].distr[ijack];
       double FV_M2= FV_pars_for_Br.distr_list[3].distr[ijack];
       double FV_B = FV_pars_for_Br.distr_list[4].distr[ijack];
       double FV_Bxg = FV_pars_for_Br.distr_list[5].distr[ijack];
       
       auto FV_Br = [&FV_A, &FV_M, &FV_A2, &FV_M2, &FV_B, &FV_Bxg](double x) {
	 
	 
	 double E_res= sqrt( FV_M*FV_M + x*x/4.0);
	 double E_res_2= sqrt( FV_M2*FV_M2 + x*x/4.0);
	 return  FV_A/( E_res*( E_res - (1.0 - x/2.0)) )  + FV_A2/(E_res_2*( E_res_2 - (1.0 - x/2.0)) ) + FV_B + FV_Bxg*x  ;
	 
       };
       
       auto FA_Br_spline = [&FA_spline_for_Br, &ijack](double x) { return FA_spline_for_Br[ijack](x);};
       auto FV_Br_spline = [&FV_spline_for_Br, &ijack](double x) { return FV_spline_for_Br[ijack](x);};
       
       double fp=fDs_FLAG_distr.distr[ijack];

       double Vcs_2 = pow(Vcs_FLAG_distr.distr[ijack],2);
       
       Br_diff_PT.distr_list[ixg_d].distr.push_back( Vcs_2*Compute_Ds_lnugamma_differential_decay_rate(re,xg_diff[ixg_d],ijack, NJ,  FV_Br, FA_Br, fp,  "PT" ));
       Br_diff_INT.distr_list[ixg_d].distr.push_back( Vcs_2*Compute_Ds_lnugamma_differential_decay_rate(re,xg_diff[ixg_d],ijack, NJ,  FV_Br, FA_Br, fp,  "INT" ));
       Br_diff_SD.distr_list[ixg_d].distr.push_back( Vcs_2*Compute_Ds_lnugamma_differential_decay_rate(re,xg_diff[ixg_d],ijack, NJ,  FV_Br, FA_Br, fp,  "SD" ));
       Br_diff_TOT.distr_list[ixg_d].distr.push_back( Vcs_2*Compute_Ds_lnugamma_differential_decay_rate(re,xg_diff[ixg_d], ijack, NJ,  FV_Br, FA_Br, fp, "TOT" ));

       Br_diff_spline_PT.distr_list[ixg_d].distr.push_back( Vcs_2*Compute_Ds_lnugamma_differential_decay_rate(re,xg_diff[ixg_d],ijack, NJ,  FV_Br_spline, FA_Br_spline, fp,  "PT" ));
       Br_diff_spline_INT.distr_list[ixg_d].distr.push_back( Vcs_2*Compute_Ds_lnugamma_differential_decay_rate(re, xg_diff[ixg_d],ijack, NJ,  FV_Br_spline, FA_Br_spline, fp,  "INT" ));
       Br_diff_spline_SD.distr_list[ixg_d].distr.push_back( Vcs_2*Compute_Ds_lnugamma_differential_decay_rate(re,xg_diff[ixg_d],ijack, NJ,  FV_Br_spline, FA_Br_spline, fp,  "SD" ));
       Br_diff_spline_TOT.distr_list[ixg_d].distr.push_back( Vcs_2*Compute_Ds_lnugamma_differential_decay_rate(re,xg_diff[ixg_d], ijack, NJ,  FV_Br_spline, FA_Br_spline, fp, "TOT" ));
     }
  }


  Print_To_File({}, { xg_cut, Br_TOT.ave(), Br_TOT.err(), Br_PT.ave(), Br_PT.err(), Br_INT.ave(), Br_INT.err(), Br_SD.ave(), Br_SD.err()}, "../data/ph_emission/"+ph_type+"/"+Meson+"/rate/rate_Ds_enugamma_"+parameterization_to_be_used_in_Br+".dat", "", "#xg_cut TOT   PT    INT     SD");

  Print_To_File({}, { xg_diff, Br_diff_TOT.ave(), Br_diff_TOT.err(), Br_diff_PT.ave(), Br_diff_PT.err(), Br_diff_INT.ave(), Br_diff_INT.err(), Br_diff_SD.ave(), Br_diff_SD.err()}, "../data/ph_emission/"+ph_type+"/"+Meson+"/rate/rate_diff_Ds_enugamma_"+parameterization_to_be_used_in_Br+".dat", "", "#xg_cut TOT   PT    INT     SD");

  Print_To_File({}, { xg_cut, Br_spline_TOT.ave(), Br_spline_TOT.err(), Br_spline_PT.ave(), Br_spline_PT.err(), Br_spline_INT.ave(), Br_spline_INT.err(), Br_spline_SD.ave(), Br_spline_SD.err()}, "../data/ph_emission/"+ph_type+"/"+Meson+"/rate/rate_Ds_enugamma_spline.dat", "", "#xg_cut TOT   PT    INT     SD");

  Print_To_File({}, { xg_diff, Br_diff_spline_TOT.ave(), Br_diff_spline_TOT.err(), Br_diff_spline_PT.ave(), Br_diff_spline_PT.err(), Br_diff_spline_INT.ave(), Br_diff_spline_INT.err(), Br_diff_spline_SD.ave(), Br_diff_spline_SD.err()}, "../data/ph_emission/"+ph_type+"/"+Meson+"/rate/rate_diff_Ds_enugamma_spline.dat", "", "#xg_cut TOT   PT    INT     SD");

  
  cout<<"With DE= 10MeV: "<<endl;
  cout<<"electron: "<<Br_spline_TOT.ave(0)<<" +- "<<Br_spline_TOT.err(0)<<", PT: "<<Br_spline_PT.ave(0)<<" +- "<<Br_spline_PT.err(0)<<endl;
  cout<<"muon: "<<Br_spline_muon_TOT.ave()<<" +- "<<Br_spline_muon_TOT.err()<<", PT: "<<Br_spline_muon_PT.ave()<<" +- "<<Br_spline_muon_PT.err()<<endl;
  cout<<"tau: "<<Br_spline_tau_TOT.ave()<<" +- "<<Br_spline_tau_TOT.err()<<", PT: "<<Br_spline_tau_PT.ave()<<" +- "<<Br_spline_tau_PT.err()<<endl;
   

  
  return;
}


rt_FF Get_form_factors_Nissa(int num_xg, int Perform_continuum_extrapolation, bool Use_three_finest, bool Include_a4, bool UseJack, string Fit_tag, string path_list ) {

  cout<<"CALLING Get_form_factors_Nissa"<<endl;
  
  //return type
  rt_FF rt_fit(UseJack);
  rt_fit.Use_three_finest=Use_three_finest;
  rt_fit.Include_a4 = Include_a4;
  rt_fit.num_xg=num_xg;
  
  string _Fit_tag = "";

  if(Fit_tag.size() > 0) { _Fit_tag= "_"+Fit_tag.substr(0, Fit_tag.size()-1); } 

  const int Njacks=((UseJack==true)?NJ:Nboots);

   
  string ph_type_mes=ph_type+"/"+Meson;
  
  
  //axial
  vector<vector<vector<data_t>>> C_A_Fu_data(size_mu_nu), C_A_Fd_data(size_mu_nu), C_A_Bu_data(size_mu_nu), C_A_Bd_data(size_mu_nu);
  //vector
  vector<vector<vector<data_t>>> C_V_Fu_data(size_mu_nu), C_V_Fd_data(size_mu_nu), C_V_Bu_data(size_mu_nu), C_V_Bd_data(size_mu_nu);



  auto Sort_light_confs = [](string A, string B) {


    //return A<B;
			     
    int conf_length_A= A.length();
    int conf_length_B= B.length();

    int pos_a_slash=-1;
    int pos_b_slash=-1;
    for(int i=0;i<conf_length_A;i++) if(A.substr(i,1)=="/") pos_a_slash=i;
    for(int j=0;j<conf_length_B;j++) if(B.substr(j,1)=="/") pos_b_slash=j;
    
    string A_bis= A.substr(pos_a_slash+1);
    string B_bis= B.substr(pos_b_slash+1);

			     
    string conf_num_A = A_bis.substr(0,4);
    string conf_num_B = B_bis.substr(0,4);

  
    if(A_bis.length() <= 4) return A_bis < B_bis;  // ################## BE CAREFUL, JUST MODIFIED STANDARD ORDERING!!! ############ return A_bis < B_bis;
     
    string rA = A_bis.substr(A_bis.length()-2);
    string rB = B_bis.substr(B_bis.length()-2);

        
    if(rA.substr(0,1) == "r") { 
      int n1 = stoi(rA.substr(1,1));
      int n2 = stoi(rB.substr(1,1));
      if(rA == rB) {
	if(rA=="r0" || rA=="r2") return conf_num_A > conf_num_B;
	else if(rA=="r1" || rA=="r3") return conf_num_A < conf_num_B;
	else crash("stream not recognized");
      }
      else return n1<n2;
    }
    // return A_bis<B_bis;  ################## BE CAREFUL, JUST MODIFIED STANDARD ORDERING!!! ############àà
    return A_bis < B_bis;
   };
  
  
 

  


  //2pts
  data_t data_2pts, data_2pts_SM;
  data_t data_2pts_V1, data_2pts_V2, data_2pts_V3;
  data_t data_2pts_u_V1, data_2pts_u_V2, data_2pts_u_V3;
  data_t data_2pts_d_V1, data_2pts_d_V2, data_2pts_d_V3;

  for(int mu=0;mu<size_mu_nu;mu++) {

    C_A_Fu_data[mu].resize(size_mu_nu);
    C_A_Fd_data[mu].resize(size_mu_nu);
    C_A_Bu_data[mu].resize(size_mu_nu);
    C_A_Bd_data[mu].resize(size_mu_nu);
    C_V_Fu_data[mu].resize(size_mu_nu);
    C_V_Fd_data[mu].resize(size_mu_nu);
    C_V_Bu_data[mu].resize(size_mu_nu);
    C_V_Bd_data[mu].resize(size_mu_nu);

    for(int nu=0;nu<size_mu_nu;nu++) {

      C_A_Fu_data[mu][nu].resize(num_xg);
      C_A_Fd_data[mu][nu].resize(num_xg);
      C_A_Bu_data[mu][nu].resize(num_xg);
      C_A_Bd_data[mu][nu].resize(num_xg);
      C_V_Fu_data[mu][nu].resize(num_xg);
      C_V_Fd_data[mu][nu].resize(num_xg);
      C_V_Bu_data[mu][nu].resize(num_xg);
      C_V_Bd_data[mu][nu].resize(num_xg);


    }

  }


  
  int off_i = (Is_reph?1:0);
  

 
  //Read data

  //2pts function

  data_2pts.Read("../new_vph_gpu_data_w_123"+path_list, "mes_contr_2pts_3", "P5P5", Sort_light_confs);
  data_2pts_SM.Read("../new_vph_gpu_data_w_123"+path_list, "mes_contr_2pts_SM_3", "P5P5", Sort_light_confs);
  data_2pts_V1.Read("../new_vph_gpu_data_w_123"+path_list, "mes_contr_2pts_3", "V1V1", Sort_light_confs);
  data_2pts_V2.Read("../new_vph_gpu_data_w_123"+path_list, "mes_contr_2pts_3", "V2V2", Sort_light_confs);
  data_2pts_V3.Read("../new_vph_gpu_data_w_123"+path_list, "mes_contr_2pts_3", "V3V3", Sort_light_confs);
  data_2pts_d_V1.Read("../new_vph_gpu_data_w_123"+path_list, "mes_contr_2pts_1", "V1V1", Sort_light_confs);
  data_2pts_d_V2.Read("../new_vph_gpu_data_w_123"+path_list, "mes_contr_2pts_1", "V2V2", Sort_light_confs);
  data_2pts_d_V3.Read("../new_vph_gpu_data_w_123"+path_list, "mes_contr_2pts_1", "V3V3", Sort_light_confs);
  data_2pts_u_V1.Read("../new_vph_gpu_data_w_123"+path_list, "mes_contr_2pts_2", "V1V1" , Sort_light_confs);
  data_2pts_u_V2.Read("../new_vph_gpu_data_w_123"+path_list, "mes_contr_2pts_2", "V2V2", Sort_light_confs);
  data_2pts_u_V3.Read("../new_vph_gpu_data_w_123"+path_list, "mes_contr_2pts_2", "V3V3", Sort_light_confs);
  

  //read data
  for(int ixg=0;ixg<num_xg;ixg++) {


    for(int mu=0;mu<size_mu_nu;mu++) {

      for(int nu=0;nu<size_mu_nu;nu++) {

	//vector
	//Fu
	C_V_Fu_data[mu][nu][ixg].Read("../new_vph_gpu_data_w_123"+path_list, "C_mu_"+to_string(mu+off_i)+"_FF_u_ixg_"+to_string(ixg), "V"+to_string(nu+off_i)+"P5", Sort_light_confs );
	//Fd
	C_V_Fd_data[mu][nu][ixg].Read("../new_vph_gpu_data_w_123"+path_list, "C_mu_"+to_string(mu+off_i)+"_FF_d_ixg_"+to_string(ixg), "V"+to_string(nu+off_i)+"P5", Sort_light_confs );
	//Bu
	C_V_Bu_data[mu][nu][ixg].Read("../new_vph_gpu_data_w_123"+path_list, "C_mu_"+to_string(mu+off_i)+"_BB_u_ixg_"+to_string(ixg), "V"+to_string(nu+off_i)+"P5", Sort_light_confs );
	//Bd
	C_V_Bd_data[mu][nu][ixg].Read("../new_vph_gpu_data_w_123"+path_list, "C_mu_"+to_string(mu+off_i)+"_BB_d_ixg_"+to_string(ixg), "V"+to_string(nu+off_i)+"P5", Sort_light_confs );

	//axial
	C_A_Fu_data[mu][nu][ixg].Read("../new_vph_gpu_data_w_123"+path_list, "C_mu_"+to_string(mu+off_i)+"_FF_u_ixg_"+to_string(ixg), "A"+to_string(nu+off_i)+"P5", Sort_light_confs );
	//Fd
	C_A_Fd_data[mu][nu][ixg].Read("../new_vph_gpu_data_w_123"+path_list, "C_mu_"+to_string(mu+off_i)+"_FF_d_ixg_"+to_string(ixg), "A"+to_string(nu+off_i)+"P5", Sort_light_confs );
	//Bu
	C_A_Bu_data[mu][nu][ixg].Read("../new_vph_gpu_data_w_123"+path_list, "C_mu_"+to_string(mu+off_i)+"_BB_u_ixg_"+to_string(ixg), "A"+to_string(nu+off_i)+"P5", Sort_light_confs );
	//Bd
	C_A_Bd_data[mu][nu][ixg].Read("../new_vph_gpu_data_w_123"+path_list, "C_mu_"+to_string(mu+off_i)+"_BB_d_ixg_"+to_string(ixg), "A"+to_string(nu+off_i)+"P5", Sort_light_confs );

      }
    }
  }


  int Nens = data_2pts.size;
  GaussianMersenne GM(543543);

  //vectors to store FV and FA for all ensembles and values of gamma
  vector<distr_t_list> FA_per_ens(Nens), FV_per_ens(Nens), FA_u_per_ens(Nens), FV_u_per_ens(Nens), FA_d_per_ens(Nens), FV_d_per_ens(Nens);
  vector<distr_t_list> FA_per_kin(num_xg-1), FV_per_kin(num_xg-1), FA_u_per_kin(num_xg-1), FV_u_per_kin(num_xg-1), FA_d_per_kin(num_xg-1), FV_d_per_kin(num_xg-1); //num_xg -1 is to exclude xg=0 which is undefined
  distr_t_list a_distr_list(UseJack);
  distr_t_list a_distr_list_red(UseJack);
  Vfloat L_list, L_list_red;
  //M_Ds, F_Ds
  distr_t_list MP_list(UseJack);
  distr_t_list MP_list_red(UseJack);
  distr_t_list FP_list(UseJack);
  distr_t_list MP_ov_FP_list(UseJack);
  //resize
  
  
  //define lambda function to combine FF and BB

  //resample RCs
  distr_t ZA_A(UseJack), ZA_B(UseJack), ZA_C(UseJack), ZA_D(UseJack);
  distr_t ZV_A(UseJack), ZV_B(UseJack), ZV_C(UseJack), ZV_D(UseJack);
  distr_t a_A(UseJack), a_B(UseJack), a_C(UseJack), a_D(UseJack);

  double fmTGeV= 1.0/0.197327;

  LatticeInfo L_info_A, L_info_B, L_info_C, L_info_D;
  L_info_A.LatInfo_new_ens("cA211a.12.48");
  L_info_B.LatInfo_new_ens("cB211b.072.96");
  L_info_C.LatInfo_new_ens("cC211a.06.80");
  L_info_D.LatInfo_new_ens("cD211a.054.96");
  

  for(int ijack=0; ijack<Njacks;ijack++) {

    ZA_A.distr.push_back( L_info_A.Za_WI_strange + GM()*L_info_A.Za_WI_strange_err/((UseJack==true)?sqrt(Njacks -1.0):1.0));
    ZV_A.distr.push_back( L_info_A.Zv_WI_strange + GM()*L_info_A.Zv_WI_strange_err/((UseJack==true)?sqrt(Njacks -1.0):1.0));

    ZA_B.distr.push_back( L_info_B.Za_WI_strange + GM()*L_info_B.Za_WI_strange_err/((UseJack==true)?sqrt(Njacks -1.0):1.0));
    ZV_B.distr.push_back( L_info_B.Zv_WI_strange + GM()*L_info_B.Zv_WI_strange_err/((UseJack==true)?sqrt(Njacks -1.0):1.0));

    ZA_C.distr.push_back( L_info_C.Za_WI_strange + GM()*L_info_C.Za_WI_strange_err/((UseJack==true)?sqrt(Njacks -1.0):1.0));
    ZV_C.distr.push_back( L_info_C.Zv_WI_strange + GM()*L_info_C.Zv_WI_strange_err/((UseJack==true)?sqrt(Njacks -1.0):1.0));

    ZA_D.distr.push_back( L_info_D.Za_WI_strange + GM()*L_info_D.Za_WI_strange_err/((UseJack==true)?sqrt(Njacks -1.0):1.0));
    ZV_D.distr.push_back( L_info_D.Zv_WI_strange + GM()*L_info_D.Zv_WI_strange_err/((UseJack==true)?sqrt(Njacks -1.0):1.0));

    a_A.distr.push_back( L_info_A.a_from_afp*fmTGeV + GM()*L_info_A.a_from_afp_err*fmTGeV/((UseJack==true)?sqrt(Njacks -1.0):1.0));
    a_B.distr.push_back( L_info_B.a_from_afp*fmTGeV + GM()*L_info_B.a_from_afp_err*fmTGeV/((UseJack==true)?sqrt(Njacks -1.0):1.0));
    a_C.distr.push_back( L_info_C.a_from_afp*fmTGeV + GM()*L_info_C.a_from_afp_err*fmTGeV/((UseJack==true)?sqrt(Njacks -1.0):1.0));
    a_D.distr.push_back( L_info_D.a_from_afp*fmTGeV + GM()*L_info_D.a_from_afp_err*fmTGeV/((UseJack==true)?sqrt(Njacks -1.0):1.0));

  }

 

  for(int iens=0;iens<Nens;iens++) {

    boost::filesystem::create_directory("../data/ph_emission/"+ph_type_mes+"/C/"+data_2pts.Tag[iens]);
    boost::filesystem::create_directory("../data/ph_emission/"+ph_type_mes+"/H/"+data_2pts.Tag[iens]);
    boost::filesystem::create_directory("../data/ph_emission/"+ph_type_mes+"/FF/"+data_2pts.Tag[iens]);
    boost::filesystem::create_directory("../data/ph_emission/"+ph_type_mes+"/FF/"+data_2pts.Tag[iens]+"/fit_results");
    boost::filesystem::create_directory("../data/ph_emission/"+ph_type_mes+"/FF_u/"+data_2pts.Tag[iens]);
    boost::filesystem::create_directory("../data/ph_emission/"+ph_type_mes+"/FF_u/"+data_2pts.Tag[iens]+"/fit_results");
    boost::filesystem::create_directory("../data/ph_emission/"+ph_type_mes+"/FF_d/"+data_2pts.Tag[iens]);
    boost::filesystem::create_directory("../data/ph_emission/"+ph_type_mes+"/FF_d/"+data_2pts.Tag[iens]+"/fit_results");
    boost::filesystem::create_directory("../data/ph_emission/"+ph_type_mes+"/mass/"+data_2pts.Tag[iens]);
    boost::filesystem::create_directory("../data/ph_emission/"+ph_type_mes+"/decay_const/"+data_2pts.Tag[iens]);
    
    cout<<"Analyzing ensemble: "<<data_2pts.Tag[iens]<<endl;

    //Lattice info
    LatticeInfo L_info;
    L_info.LatInfo_new_ens(data_2pts.Tag[iens]);
    CorrAnalysis Corr(UseJack, Njacks,Nboots, iens);
    Corr.Nt = data_2pts.nrows[iens];
    Corr.Reflection_sign=1;
    Corr.Perform_Nt_t_average=1;

    
    //read theta values and loop over them
    Vfloat thetas, masses_u, masses_d, virts;

    thetas= Read_From_File("../new_vph_gpu_data_w_123"+path_list+"/"+data_2pts.Tag[iens]+"/pars_list.dat", 1 , 5);
    virts=  Read_From_File("../new_vph_gpu_data_w_123"+path_list+"/"+data_2pts.Tag[iens]+"/pars_list.dat", 2 , 5);
    masses_u= Read_From_File("../new_vph_gpu_data_w_123"+path_list+"/"+data_2pts.Tag[iens]+"/pars_list.dat", 3 , 5);
    masses_d= Read_From_File("../new_vph_gpu_data_w_123"+path_list+"/"+data_2pts.Tag[iens]+"/pars_list.dat", 4 , 5);

    cout<<"pars_list.dat: Read!"<<endl;

    //if((signed)thetas.size() != num_xg) crash("Number of rows in pars_list.dat does not match num_xg"); 

    int num_xg_per_ens= (signed)thetas.size();
       
    //RCs
    distr_t Za, Zv, a_distr;
    if(data_2pts.Tag[iens].substr(1,1)=="A") { Za= ZA_A; Zv=ZV_A;a_distr=a_A;}
    else if(data_2pts.Tag[iens].substr(1,1)=="B") { Za= ZA_B; Zv=ZV_B;a_distr=a_B;}
    else if(data_2pts.Tag[iens].substr(1,1)=="C") { Za= ZA_C; Zv=ZV_C; a_distr=a_C;}
    else if(data_2pts.Tag[iens].substr(1,1)=="D") { Za= ZA_D; Zv=ZV_D;a_distr=a_D;}
    else crash("Ensemble: "+data_2pts.Tag[iens]+" not recognised");


    a_distr_list.distr_list.push_back(a_distr);
    L_list.push_back(L_info.L);
    
    if(data_2pts.Tag[iens] != "cB211b.072.96") {
      a_distr_list_red.distr_list.push_back(a_distr);
      L_list_red.push_back(L_info.L);
    }
   
    

    //read masses
    double mu= masses_u[0];
    double md= masses_d[0];

    
    cout<<"ZA: "<<Za.ave()<<" +- "<<Za.err()<<endl;
    cout<<"ZV: "<<Zv.ave()<<" +- "<<Zv.err()<<endl;
    cout<<"mu: "<<mu<<endl;
    cout<<"md: "<<md<<endl;
    cout<<"Number of kinematic present: "<<num_xg_per_ens<<endl;
    cout<<"Number of kinematic analyzed: "<<num_xg<<endl;
    
    
    
   
    
    
   
   

    //set time interval for eff_mass_fit
    if(data_2pts.Tag[iens].substr(1,1) =="A") {Corr.Tmin=19; Corr.Tmax=35;}
    else if(data_2pts.Tag[iens] =="cB211b.072.64") {Corr.Tmin=25; Corr.Tmax=45;}
    else if(data_2pts.Tag[iens] =="cB211b.072.96") {Corr.Tmin=25; Corr.Tmax=45;}
    else if(data_2pts.Tag[iens].substr(1,1) == "C") {Corr.Tmin=33;Corr.Tmax=51;}
    else if(data_2pts.Tag[iens].substr(1,1) == "D") {Corr.Tmin=35;Corr.Tmax=54;}
    else crash("In fixing [Tmin, Tmax] for MP, Ensemble: "+data_2pts.Tag[iens]+" not recognized");

    distr_t_list pt2_distr= Corr.corr_t(data_2pts.col(0)[iens], "../data/ph_emission/"+ph_type_mes+"/"+"C/"+data_2pts.Tag[iens]+"/corr_2pt.dat");
    distr_t_list eff_mass = Corr.effective_mass_t(pt2_distr, "../data/ph_emission/"+ph_type_mes+"/"+"mass/"+data_2pts.Tag[iens]+"/eff_mass.dat");
    distr_t_list fp_distr= Corr.decay_constant_t( pow( mu+md,2)*pt2_distr, "../data/ph_emission/"+ph_type_mes+"/"+"decay_const/"+data_2pts.Tag[iens]+"/decay_const.dat");
    distr_t M_P=Corr.Fit_distr(eff_mass);
    distr_t F_P=Corr.Fit_distr(fp_distr);
    MP_list.distr_list.push_back(M_P/a_distr);
    FP_list.distr_list.push_back(F_P/a_distr);
    MP_ov_FP_list.distr_list.push_back( M_P/F_P);
    if(data_2pts.Tag[iens] != "cB211b.072.96") {
      MP_list_red.distr_list.push_back(M_P/a_distr);
    }


    
   
    
    //set time interval for eff_mass_fit SM
    if(data_2pts.Tag[iens].substr(1,1) =="A") {Corr.Tmin=24; Corr.Tmax=35;}
    else if(data_2pts.Tag[iens] =="cB211b.072.64") {Corr.Tmin=20; Corr.Tmax=36;}
    else if(data_2pts.Tag[iens] =="cB211b.072.96") {Corr.Tmin=20; Corr.Tmax=36;}
    else if(data_2pts.Tag[iens].substr(1,1) == "C")  {Corr.Tmin=33;Corr.Tmax=51;}
    else if(data_2pts.Tag[iens].substr(1,1) == "D")  {Corr.Tmin=41;Corr.Tmax=54;}
    else crash("In fixing [Tmin, Tmax] for smeared MP, Ensemble: "+data_2pts.Tag[iens]+" not recognized");
    distr_t_list pt2_distr_SM= Corr.corr_t(data_2pts_SM.col(0)[iens], "../data/ph_emission/"+ph_type_mes+"/"+"C/"+data_2pts.Tag[iens]+"/corr_2pt_SM.dat");
    distr_t_list eff_mass_SM = Corr.effective_mass_t(pt2_distr_SM, "../data/ph_emission/"+ph_type_mes+"/"+"mass/"+data_2pts.Tag[iens]+"/eff_mass_SM.dat");
    distr_t M_P_SM = Corr.Fit_distr(eff_mass_SM);


    cout<<"aM_P: "<<M_P.ave()<<" +- "<<M_P.err()<<" -> "<< (M_P/(a_distr)).ave()<<" +- "<<(M_P/a_distr).err()<<" GeV"<<endl;
    cout<<"aF_P: "<<F_P.ave()<<" +- "<<F_P.err()<<" -> "<< (F_P/(a_distr)).ave()<<" +- "<<(F_P/a_distr).err()<<" GeV"<<endl;
    cout<<"MP/FP: "<<(M_P/F_P).ave()<<" +- "<<(M_P/F_P).err()<<endl;
    cout<<"aM_P(SM): "<<M_P_SM.ave()<<" +- "<<M_P_SM.err()<<" -> "<< (M_P_SM/(a_distr)).ave()<<" +- "<<(M_P_SM/a_distr).err()<<" GeV"<<endl;
    cout<<"MP/FP: "<<(M_P/F_P).ave()<<" +- "<<(M_P/F_P).err()<<endl;



    //vector ud meson
    distr_t_list pt2_vector_distr= (1.0/3.0)*Corr.corr_t(summ_master(data_2pts_V1.col(0)[iens], data_2pts_V2.col(0)[iens], data_2pts_V3.col(0)[iens]), "../data/ph_emission/"+ph_type_mes+"/"+"C/"+data_2pts.Tag[iens]+"/corr_2pt_V_ud.dat");
    distr_t_list eff_mass_V = Corr.effective_mass_t(pt2_vector_distr, "../data/ph_emission/"+ph_type_mes+"/"+"mass/"+data_2pts.Tag[iens]+"/eff_mass_V_ud.dat");
    //vector u meson
    distr_t_list uu_vector_distr= (1.0/3.0)*Corr.corr_t(summ_master(data_2pts_u_V1.col(0)[iens], data_2pts_u_V2.col(0)[iens], data_2pts_u_V3.col(0)[iens]), "../data/ph_emission/"+ph_type_mes+"/"+"C/"+data_2pts.Tag[iens]+"/corr_2pt_V_uu.dat");
    distr_t_list eff_mass_V_uu= Corr.effective_mass_t(uu_vector_distr, "../data/ph_emission/"+ph_type_mes+"/"+"mass/"+data_2pts.Tag[iens]+"/eff_mass_V_uu.dat");
    //vector d meson
    distr_t_list dd_vector_distr= (1.0/3.0)*Corr.corr_t(summ_master(data_2pts_d_V1.col(0)[iens], data_2pts_d_V2.col(0)[iens], data_2pts_d_V3.col(0)[iens]), "../data/ph_emission/"+ph_type_mes+"/"+"C/"+data_2pts.Tag[iens]+"/corr_2pt_V_dd.dat");
    distr_t_list eff_mass_V_dd= Corr.effective_mass_t(dd_vector_distr, "../data/ph_emission/"+ph_type_mes+"/"+"mass/"+data_2pts.Tag[iens]+"/eff_mass_V_dd.dat");


    //define meson mass exponential to be removed
    auto EXP_MES_FUNC = [&] (double a, double b, double c) { return (b<c/2)?1.0/(exp(-a*b)):1.0/(exp(-a*(c-b)));};

    distr_t_list EXP_MES= distr_t_list::f_of_distr( EXP_MES_FUNC, M_P, Corr.Nt);

   

    vector<vector<vector<distr_t_list>>> Ax_glb, Vec_glb, Ax_u_glb, Ax_d_glb, Vec_u_glb, Vec_d_glb;
    distr_t_list FV(UseJack), FA(UseJack), FV_u(UseJack), FV_d(UseJack), FA_u(UseJack), FA_d(UseJack);


    distr_t_list xg_list(UseJack);

    Vfloat ax_Tmin, ax_Tmax;
    Vfloat vec_Tmin, vec_Tmax;

    Vfloat ax_u_Tmin, ax_u_Tmax;
    Vfloat vec_u_Tmin, vec_u_Tmax;

    Vfloat ax_d_Tmin, ax_d_Tmax;
    Vfloat vec_d_Tmin, vec_d_Tmax;
    

    for(int ixg=0;ixg<num_xg;ixg++) {

      vector<vector<distr_t_list>> Ax_tens(size_mu_nu);
      vector<vector<distr_t_list>> Ax_tens_u(size_mu_nu);
      vector<vector<distr_t_list>> Ax_tens_d(size_mu_nu);
      vector<vector<distr_t_list>> Vec_tens(size_mu_nu);
      vector<vector<distr_t_list>> Vec_tens_u(size_mu_nu);
      vector<vector<distr_t_list>> Vec_tens_d(size_mu_nu);
      


      //get xg, Eg, kz from thetas

      double theta=thetas[ixg];

      pt3_momenta pt3_mom(0.0, 0.0, thetas[ixg]/2.0, masses_u[ixg], masses_d[ixg], virts[ixg], L_info.L, L_info.T);

      double Eg= pt3_mom.Egamma();
      distr_t xg= pt3_mom.x_gamma(M_P);
      if(xg.ave() > 1e-10 ) xg_list.distr_list.push_back(xg);
      double kz = pt3_mom.k()[2];
    
   

      //define photon exponential to be removed
      Vfloat EXP_PH(Corr.Nt,0.0);
      Vfloat EXP_PH_RED(Corr.Nt, 0.0);
      for(int t=0; t < Corr.Nt;t++) EXP_PH[t] = exp( Eg*abs( Corr.Nt/2 - t));
      for(int t=0; t < Corr.Nt;t++) EXP_PH_RED[t] = exp( -Eg*t);
     


      cout<<"##### Considering kinematic with......"<<endl;
      cout<<"Eg: "<<Eg<<endl;
      cout<<"xg: "<<xg.ave()<<" +- "<<xg.err()<<endl;
      cout<<"thz: "<<theta<<endl;
      cout<<"kz: "<<kz<<endl;

      

    
    

      //define vectors to combine FF and BB
      Vfloat th_FF(Corr.Nt);
      Vfloat th_BB(Corr.Nt);
      for(int t=0;t<Corr.Nt;t++) {
	if(t<Corr.Nt/2) {th_FF[t]=1.0; th_BB[t] = 0.0;}
	else {th_BB[t] = 1.0; th_FF[t] = 0.0;}
      }

      //loop over mu and nu
      for(int mu=0;mu<size_mu_nu;mu++) {
	for(int nu=0;nu<size_mu_nu;nu++) {

	  cout<<"Analyzing mu: "<<mu+off_i<<" nu: "<<nu+off_i<<endl;
	  int Im_Re;
	  double parity;

	  //vector
	
	  Corr.Reflection_sign = -1;
	  Im_Re=1;
	  parity=1.0;

	  Corr.Perform_Nt_t_average = 0;
	  distr_t_list vec_F_u = qu*Corr.corr_t(C_V_Fu_data[mu][nu][ixg].col(Im_Re)[iens],"");
	  distr_t_list vec_B_u = qu*Corr.corr_t(C_V_Bu_data[mu][nu][ixg].col(Im_Re)[iens],"");
	  distr_t_list vec_F_d = qd*Corr.corr_t(C_V_Fd_data[mu][nu][ixg].col(Im_Re)[iens],"");
	  distr_t_list vec_B_d = qd*Corr.corr_t(C_V_Bd_data[mu][nu][ixg].col(Im_Re)[iens],"");
	
	  distr_t_list vec_u = th_FF*vec_F_u + th_BB*vec_B_u;
	  distr_t_list vec_d = th_FF*vec_F_d + th_BB*vec_B_d;
	  distr_t_list vec = vec_u + vec_d;
	  distr_t_list vec_symm= vec;
	  distr_t_list vec_u_symm= vec_u;
	  distr_t_list vec_d_symm= vec_d;
	  //symmetrize vec
	  for(int t=0; t<Corr.Nt;t++) {
	    vec_symm.distr_list[t] = 0.5*(vec.distr_list[t] + Corr.Reflection_sign*vec.distr_list[( Corr.Nt -t)%Corr.Nt]);
	    vec_u_symm.distr_list[t] = 0.5*(vec_u.distr_list[t] + Corr.Reflection_sign*vec_u.distr_list[( Corr.Nt -t)%Corr.Nt]);
	    vec_d_symm.distr_list[t] = 0.5*(vec_d.distr_list[t] + Corr.Reflection_sign*vec_d.distr_list[( Corr.Nt -t)%Corr.Nt]);
	  }
	  Corr.Perform_Nt_t_average=1;
	
	

	  //axial
	  if( (mu==0 || nu==0) && (mu != 0 || nu != 0)) {Im_Re=1; Corr.Reflection_sign=-1; parity=1.0;}
	  else { Im_Re=0; Corr.Reflection_sign=1; parity=1.0;}

	  Corr.Perform_Nt_t_average=0;
	  distr_t_list ax_F_u = qu*Corr.corr_t(C_A_Fu_data[mu][nu][ixg].col(Im_Re)[iens],"");
	  distr_t_list ax_B_u = qu*Corr.corr_t(C_A_Bu_data[mu][nu][ixg].col(Im_Re)[iens],"");
	  distr_t_list ax_F_d = qd*Corr.corr_t(C_A_Fd_data[mu][nu][ixg].col(Im_Re)[iens],"");
	  distr_t_list ax_B_d = qd*Corr.corr_t(C_A_Bd_data[mu][nu][ixg].col(Im_Re)[iens],"");

	
	  distr_t_list ax_u = parity*(ax_F_u*th_FF  + ax_B_u*th_BB);
	  distr_t_list ax_d = parity*(ax_F_d*th_FF  + ax_B_d*th_BB);
	  distr_t_list ax = ax_u -ax_d;
	  distr_t_list ax_symm=ax;
	  distr_t_list ax_u_symm= ax_u;
	  distr_t_list ax_d_symm= ax_d;
	  //symmetrize ax
	  for(int t=0; t<Corr.Nt;t++) {
	    ax_symm.distr_list[t] = 0.5*(ax.distr_list[t] + Corr.Reflection_sign*ax.distr_list[( Corr.Nt -t)%Corr.Nt]);
	    ax_u_symm.distr_list[t] = 0.5*(ax_u.distr_list[t] + Corr.Reflection_sign*ax_u.distr_list[( Corr.Nt -t)%Corr.Nt]);
	    ax_d_symm.distr_list[t] = 0.5*(ax_d.distr_list[t] + Corr.Reflection_sign*ax_d.distr_list[( Corr.Nt -t)%Corr.Nt]);
	  }
	  Corr.Perform_Nt_t_average=1;

	  //restore standard reflection sign
	  Corr.Reflection_sign=1;


	  //get H tensor
	  distr_t_list HA_u= ax_u*EXP_MES*EXP_PH;
	  distr_t_list HA_d= ax_d*EXP_MES*EXP_PH;
	  distr_t_list HV_u= vec_u*EXP_MES*EXP_PH;
	  distr_t_list HV_d= vec_d*EXP_MES*EXP_PH;
	  distr_t_list HA_u_symm= ax_u_symm*EXP_MES*EXP_PH;
	  distr_t_list HA_d_symm= ax_d_symm*EXP_MES*EXP_PH;
	  distr_t_list HV_u_symm= vec_u_symm*EXP_MES*EXP_PH;
	  distr_t_list HV_d_symm= vec_d_symm*EXP_MES*EXP_PH;
	  distr_t_list HA_symm= ax_symm*EXP_MES*EXP_PH;
	  distr_t_list HV_symm= vec_symm*EXP_MES*EXP_PH;



	  //########################   C TENSOR ###########################


	  //print forward and backward contribution to C
	  Print_To_File({}, {ax_F_u.ave(), ax_F_u.err(), ax_F_d.ave(), ax_F_d.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"C/"+data_2pts.Tag[iens]+"/quark_contr_A_wo_T_symm_FF_mu_"+to_string(mu+off_i)+"_nu_"+to_string(nu+off_i)+"_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "#   ax_F_u    ax_F_d");
	  Print_To_File({}, {ax_B_u.ave(), ax_B_u.err(), ax_B_d.ave(), ax_B_d.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"C/"+data_2pts.Tag[iens]+"/quark_contr_A_wo_T_symm_BB_mu_"+to_string(mu+off_i)+"_nu_"+to_string(nu+off_i)+"_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "#   ax_B_u    ax_B_d");
	  Print_To_File({}, {vec_F_u.ave(), vec_F_u.err(), vec_F_d.ave(), vec_F_d.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"C/"+data_2pts.Tag[iens]+"/quark_contr_V_wo_T_symm_FF_mu_"+to_string(mu+off_i)+"_nu_"+to_string(nu+off_i)+"_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "#   vec_F_u    vec_F_d");
	  Print_To_File({}, {vec_B_u.ave(), vec_B_u.err(), vec_B_d.ave(), vec_B_d.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"C/"+data_2pts.Tag[iens]+"/quark_contr_V_wo_T_symm_BB_mu_"+to_string(mu+off_i)+"_nu_"+to_string(nu+off_i)+"_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "#   vec_B_u    vec_B_d");
	

	  //single contributions to C
	  Print_To_File({}, {ax_u.ave() ,ax_u.err(), ax_d.ave(), ax_d.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"C/"+data_2pts.Tag[iens]+"/quark_contr_A_wo_T_symm_mu_"+to_string(mu+off_i)+"_nu_"+to_string(nu+off_i)+"_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "#   ax_u    ax_d");
	  Print_To_File({}, {vec_u.ave() ,vec_u.err(), vec_d.ave(), vec_d.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"C/"+data_2pts.Tag[iens]+"/quark_contr_V_wo_T_symm_mu_"+to_string(mu+off_i)+"_nu_"+to_string(nu+off_i)+"_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "#   vec_u    vec_d");

	  //single contributions to C with symmetrization
	  Print_To_File({}, {ax_u_symm.ave() ,ax_u_symm.err(), ax_d_symm.ave(), ax_d_symm.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"C/"+data_2pts.Tag[iens]+"/quark_contr_A_mu_"+to_string(mu+off_i)+"_nu_"+to_string(nu+off_i)+"_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "#   ax_u    ax_d");
	  Print_To_File({}, {vec_u_symm.ave() ,vec_u_symm.err(), vec_d_symm.ave(), vec_d_symm.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"C/"+data_2pts.Tag[iens]+"/quark_contr_V_mu_"+to_string(mu+off_i)+"_nu_"+to_string(nu+off_i)+"_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "#   vec_u    vec_d");

	  //error scaling of C^munu(t)
	  for(int t=1;t< Corr.Nt/2;t++) {

	    if(ixg==0) {
	      ofstream Print_scaling_A("../data/ph_emission/"+ph_type_mes+"/"+"C/"+data_2pts.Tag[iens]+"/C_scaling_A_mu_"+to_string(mu+off_i)+"_nu_"+to_string(nu+off_i)+"_t_"+to_string(t)+".dat");
	      Print_scaling_A<<xg.ave()<<" "<<(ax_u_symm*EXP_PH).err(t)<<" "<<(ax_d_symm*EXP_PH).err(t)<<endl;
	      Print_scaling_A.close();

	      ofstream Print_scaling_V("../data/ph_emission/"+ph_type_mes+"/"+"C/"+data_2pts.Tag[iens]+"/C_scaling_V_mu_"+to_string(mu+off_i)+"_nu_"+to_string(nu+off_i)+"_t_"+to_string(t)+".dat");
	      Print_scaling_V<<xg.ave()<<" "<<(vec_u_symm*EXP_PH).err(t)<<" "<<(vec_d_symm*EXP_PH).err(t)<<endl;
	      Print_scaling_V.close();
	    }
	    else {
	      ofstream Print_scaling_A("../data/ph_emission/"+ph_type_mes+"/"+"C/"+data_2pts.Tag[iens]+"/C_scaling_A_mu_"+to_string(mu+off_i)+"_nu_"+to_string(nu+off_i)+"_t_"+to_string(t)+".dat", ofstream::app);
	      Print_scaling_A<<xg.ave()<<" "<<(ax_u_symm*EXP_PH).err(t)<<" "<<(ax_d_symm*EXP_PH).err(t)<<endl;
	      Print_scaling_A.close();

	      ofstream Print_scaling_V("../data/ph_emission/"+ph_type_mes+"/"+"C/"+data_2pts.Tag[iens]+"/C_scaling_V_mu_"+to_string(mu+off_i)+"_nu_"+to_string(nu+off_i)+"_t_"+to_string(t)+".dat", ofstream::app);
	      Print_scaling_V<<xg.ave()<<" "<<(vec_u_symm*EXP_PH).err(t)<<" "<<(vec_d_symm*EXP_PH).err(t)<<endl;
	      Print_scaling_V.close();
	    }

	  }

	  //total contribution without symmetrization to C

	  Print_To_File({}, {ax.ave(), ax.err(), vec.ave(), vec.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"C/"+data_2pts.Tag[iens]+"/wo_T_symm_mu_"+to_string(mu+off_i)+"_nu_"+to_string(nu+off_i)+"_xg_"+to_string_with_precision(xg.ave(),2)+".dat" , "", "#  t    A     V");

	  //total contribution to C
	  Print_To_File({}, {ax_symm.ave(), ax_symm.err(), vec_symm.ave(), vec_symm.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"C/"+data_2pts.Tag[iens]+"/mu_"+to_string(mu+off_i)+"_nu_"+to_string(nu+off_i)+"_xg_"+to_string_with_precision(xg.ave(),2)+".dat" , "", "#  t    A     V");


	  //###############################################################


	  //########################   H TENSOR ###########################


	  //single contributions to H without symmetrization
	   Print_To_File({}, {HA_u.ave() ,HA_u.err(), HA_d.ave(), HA_d.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"H/"+data_2pts.Tag[iens]+"/quark_contr_A_wo_symm_mu_"+to_string(mu+off_i)+"_nu_"+to_string(nu+off_i)+"_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "#   HA_u    HA_d");
	  Print_To_File({}, {HV_u.ave() ,HV_u.err(), HV_d.ave(), HV_d.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"H/"+data_2pts.Tag[iens]+"/quark_contr_V_wo_symm_mu_"+to_string(mu+off_i)+"_nu_"+to_string(nu+off_i)+"_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "#   HV_u    HV_d");
	  

	  
	  //single contributions to H with symmetrization
	  Print_To_File({}, {HA_u_symm.ave() ,HA_u_symm.err(), HA_d_symm.ave(), HA_d_symm.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"H/"+data_2pts.Tag[iens]+"/quark_contr_A_mu_"+to_string(mu+off_i)+"_nu_"+to_string(nu+off_i)+"_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "#   HA_u    HA_d");
	  Print_To_File({}, {HV_u_symm.ave() ,HV_u_symm.err(), HV_d_symm.ave(), HV_d_symm.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"H/"+data_2pts.Tag[iens]+"/quark_contr_V_mu_"+to_string(mu+off_i)+"_nu_"+to_string(nu+off_i)+"_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "#   HV_u    HV_d");

	  //total contribution to H
	  Print_To_File({}, {HA_symm.ave(), HA_symm.err(), HV_symm.ave(), HV_symm.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"H/"+data_2pts.Tag[iens]+"/mu_"+to_string(mu+off_i)+"_nu_"+to_string(nu+off_i)+"_xg_"+to_string_with_precision(xg.ave(),2)+".dat" , "", "#  t    HA     HV");

	  //################################################################

	  //push back
	  Ax_tens[mu].push_back(ax_symm);
	  Ax_tens_u[mu].push_back(ax_u_symm);
	  Ax_tens_d[mu].push_back(ax_d_symm);
	  Vec_tens[mu].push_back(vec_symm);
	  Vec_tens_u[mu].push_back(vec_u_symm);
	  Vec_tens_d[mu].push_back(vec_d_symm);
	  
	

	}
      }

      //push_back Ax_tens and Vec_tens
      Ax_glb.push_back(Ax_tens);
      Ax_u_glb.push_back(Ax_tens_u);
      Ax_d_glb.push_back(Ax_tens_d);
      
      Vec_glb.push_back(Vec_tens);
      Vec_u_glb.push_back(Vec_tens_u);
      Vec_d_glb.push_back(Vec_tens_d);

      //Compute FV and FA
      distr_t_list FA0_distr= 0.5*(Ax_glb[0][1-off_i][1-off_i] + Ax_glb[0][2-off_i][2-off_i]);
      distr_t_list FA_distr = (0.5*(Ax_tens[1-off_i][1-off_i] + Ax_tens[2-off_i][2-off_i])*EXP_PH - FA0_distr)*(1.0/Eg)*F_P/FA0_distr;
      distr_t_list FV_distr = -0.5*(Za/Zv)*(F_P/(-1.0*FA0_distr))*( Vec_tens[1-off_i][2-off_i] - Vec_tens[2-off_i][1-off_i])*EXP_PH/kz;
      distr_t_list FV_sub_distr= -0.5*(Za/Zv)*(F_P/(-1.0*FA0_distr))*( (Vec_tens[1-off_i][2-off_i] - Vec_tens[2-off_i][1-off_i])*EXP_PH - (Vec_glb[0][1-off_i][2-off_i] + Vec_glb[0][2-off_i][1-off_i])*1.0)/kz;
      //compute FV and FA (up-component)
      distr_t_list FA0_u_distr= 0.5*(Ax_u_glb[0][1-off_i][1-off_i] + Ax_u_glb[0][2-off_i][2-off_i]);
      //distr_t_list FA0_d_distr= 0.5*(Ax_d_glb[0][1-off_i][1-off_i] + Ax_d_glb[0][2-off_i][2-off_i]);
      distr_t_list FA_u_distr = (0.5*(Ax_tens_u[1-off_i][1-off_i] + Ax_tens_u[2-off_i][2-off_i])*EXP_PH - qu*FA0_distr )*(1.0/Eg)*F_P/FA0_distr;
      //distr_t_list FA_u_distr = (0.5*(Ax_tens_u[1-off_i][1-off_i] + Ax_tens_u[2-off_i][2-off_i])*EXP_PH/(0.5*FA0_u_distr/qu +0.5*FA0_d_distr/qd) - qu )*(1.0/Eg)*F_P;
      distr_t_list FV_u_distr = -0.5*(Za/Zv)*(F_P/(-1.0*FA0_distr))*( Vec_tens_u[1-off_i][2-off_i] - Vec_tens_u[2-off_i][1-off_i])*EXP_PH/kz;
      distr_t_list FV_u_sub_distr= -0.5*(Za/Zv)*(F_P/(-1.0*FA0_distr))*( (Vec_tens_u[1-off_i][2-off_i] - Vec_tens_u[2-off_i][1-off_i])*EXP_PH - (Vec_u_glb[0][1-off_i][2-off_i] + Vec_u_glb[0][2-off_i][1-off_i])*1.0)/kz;

      //compute FV and FA (d-component)
      distr_t_list FA0_d_distr= 0.5*(Ax_d_glb[0][1-off_i][1-off_i] + Ax_d_glb[0][2-off_i][2-off_i]);
      distr_t_list FA_d_distr = (0.5*(Ax_tens_d[1-off_i][1-off_i] + Ax_tens_d[2-off_i][2-off_i])*EXP_PH -qd*FA0_distr)*(1.0/Eg)*F_P/FA0_distr;
      //distr_t_list FA_d_distr = (0.5*(Ax_tens_d[1-off_i][1-off_i] + Ax_tens_d[2-off_i][2-off_i])*EXP_PH/(0.5*FA0_u_distr/qu +0.5*FA0_d_distr/qd) - qd )*(1.0/Eg)*F_P;
      distr_t_list FV_d_distr = -0.5*(Za/Zv)*(F_P/(-1.0*FA0_distr))*( Vec_tens_d[1-off_i][2-off_i] - Vec_tens_d[2-off_i][1-off_i])*EXP_PH/kz;
      distr_t_list FV_d_sub_distr= -0.5*(Za/Zv)*(F_P/(-1.0*FA0_distr))*( (Vec_tens_d[1-off_i][2-off_i] - Vec_tens_d[2-off_i][1-off_i])*EXP_PH - (Vec_d_glb[0][1-off_i][2-off_i] + Vec_d_glb[0][2-off_i][1-off_i])*1.0)/kz;


      //Print FV and FA
      Print_To_File({}, {FA_distr.ave(), FA_distr.err(), FA0_distr.ave(), FA0_distr.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"FF/"+data_2pts.Tag[iens]+"/FA_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "");
      Print_To_File({}, {FV_distr.ave(), FV_distr.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"FF/"+data_2pts.Tag[iens]+"/FV_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "");
      Print_To_File({}, {FV_sub_distr.ave(), FV_sub_distr.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"FF/"+data_2pts.Tag[iens]+"/FV_sub_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "");

      //Print FV and FA (up-component)
      Print_To_File({}, {FA_u_distr.ave(), FA_u_distr.err(), FA0_u_distr.ave(), FA0_u_distr.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"FF_u/"+data_2pts.Tag[iens]+"/FA_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "");
      Print_To_File({}, {FV_u_distr.ave(), FV_u_distr.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"FF_u/"+data_2pts.Tag[iens]+"/FV_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "");
      Print_To_File({}, {FV_u_sub_distr.ave(), FV_u_sub_distr.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"FF_u/"+data_2pts.Tag[iens]+"/FV_sub_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "");
      
      //Print FV and FA (d-component)
      Print_To_File({}, {FA_d_distr.ave(), FA_d_distr.err(), FA0_d_distr.ave(), FA0_d_distr.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"FF_d/"+data_2pts.Tag[iens]+"/FA_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "");
      Print_To_File({}, {(-1.0*FV_d_distr).ave(), (-1.0*FV_d_distr).err()}, "../data/ph_emission/"+ph_type_mes+"/"+"FF_d/"+data_2pts.Tag[iens]+"/FV_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "");
      Print_To_File({}, {(-1.0*FV_d_sub_distr).ave(), (-1.0*FV_d_sub_distr).err()}, "../data/ph_emission/"+ph_type_mes+"/"+"FF_d/"+data_2pts.Tag[iens]+"/FV_sub_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "");

      //fit the form factors

      //set interval depending on xg
      int Tmin_V, Tmax_V, Tmin_A, Tmax_A;
        


      //axial
      Corr.Tmin= Tmin_A; Corr.Tmax= Tmax_A;
      if(xg.ave() > 1e-10) { //exclude xg=0
	    Get_Tmin_Tmax("A", Tmin_A, Tmax_A, ixg, data_2pts.Tag[iens]);
	    Corr.Tmin=Tmin_A+shift; Corr.Tmax=Tmax_A+shift;
	    FA.distr_list.push_back( Corr.Fit_distr(FA_distr));
	    FA_per_kin[ixg-1].distr_list.push_back(Corr.Fit_distr(FA_distr));
	    ax_Tmin.push_back(Tmin_A+shift);
	    ax_Tmax.push_back(Tmax_A+shift);

	    //shift u-contrib by 0.1 fm
	    int shift_u = (int)(0.1/(a_distr.ave()/fmTGeV));
	    if(data_2pts.Tag[iens].substr(1,1)=="D") shift_u +=1;
	    
	    Get_Tmin_Tmax("A", Tmin_A, Tmax_A, ixg, data_2pts.Tag[iens]);

	    if(xg.ave() > 0.75) {
	      shift_u=0;
	      Tmin_A= (int)(1.7/(a_distr.ave()/fmTGeV));
	      Tmax_A= (int)(2.1/(a_distr.ave()/fmTGeV));
	    }
	    
	    Corr.Tmin=Tmin_A+shift+shift_u; Corr.Tmax=Tmax_A+shift+shift_u;
	    FA_u.distr_list.push_back(Corr.Fit_distr(FA_u_distr));
	    FA_u_per_kin[ixg-1].distr_list.push_back(Corr.Fit_distr(FA_u_distr));
	    ax_u_Tmin.push_back(Tmin_A+shift+shift_u);
	    ax_u_Tmax.push_back(Tmax_A+shift+shift_u);

	    //shift s-contrib by 0.1 fm
	    int shift_d= 0;
	    if(xg.ave() < 0.85) {
	      shift_d = (int)(0.1/(a_distr.ave()/fmTGeV));
	    }
	    
	    Get_Tmin_Tmax("A", Tmin_A, Tmax_A, ixg, data_2pts.Tag[iens]);
	    Corr.Tmin=Tmin_A+shift+shift_d; Corr.Tmax=Tmax_A+shift+shift_d;
	    FA_d.distr_list.push_back(Corr.Fit_distr(FA_d_distr));
	    FA_d_per_kin[ixg-1].distr_list.push_back(Corr.Fit_distr(FA_d_distr));
	    ax_d_Tmin.push_back(Tmin_A+shift+shift_d);
	    ax_d_Tmax.push_back(Tmax_A+shift+shift_d);
	  
      }
      //vector
      if(xg.ave() > 1e-10) { //exclude xg=0
	Get_Tmin_Tmax("V", Tmin_V, Tmax_V, ixg, data_2pts.Tag[iens]);
	Corr.Tmin= Tmin_V+shift; Corr.Tmax= Tmax_V+shift;
	FV.distr_list.push_back( Corr.Fit_distr(FV_distr));
	FV_per_kin[ixg-1].distr_list.push_back(Corr.Fit_distr(FV_distr));
	vec_Tmin.push_back(Tmin_V+shift);
	vec_Tmax.push_back(Tmax_V+shift);

	//shift u-contrib by 0.4 fm
	int shift_u=0;
	if(xg.ave()< 0.25) {
	  shift_u = (int)(0.4/(a_distr.ave()/fmTGeV));
	  if(data_2pts.Tag[iens].substr(1,1)=="A") shift_u +=3;
	}
	else if(xg.ave() < 0.45) {	  
	   shift_u = (int)(0.4/(a_distr.ave()/fmTGeV));
	   if(data_2pts.Tag[iens].substr(1,1)=="A") shift_u +=1;
	   else if(data_2pts.Tag[iens].substr(1,1)=="D") shift_u -= 3;	  
	}
	
	
	Get_Tmin_Tmax("V", Tmin_V, Tmax_V, ixg, data_2pts.Tag[iens]);

	if( xg.ave() > 0.75) {
	  shift_u=0;
	  Tmin_V = (int)(1.4/(a_distr.ave()/fmTGeV));
	  Tmax_V = (int)(2.2/(a_distr.ave()/fmTGeV));

	}
	
	Corr.Tmin= Tmin_V+shift+shift_u; Corr.Tmax= Tmax_V+shift+shift_u;	
	FV_u.distr_list.push_back( Corr.Fit_distr(FV_u_distr));
	FV_u_per_kin[ixg-1].distr_list.push_back(Corr.Fit_distr(FV_u_distr));
	vec_u_Tmin.push_back(Tmin_V+shift+shift_u);
	vec_u_Tmax.push_back(Tmax_V+shift+shift_u);

	//shift s-contrib by 0.3 fm
	int shift_d=0;
	if(xg.ave() < 0.55) {
	  shift_d = (int)(0.3/(a_distr.ave()/fmTGeV));
	}
	else if (xg.ave() < 0.85) {
	  shift_d = (int)(0.1/(a_distr.ave()/fmTGeV));
	}

	
	Get_Tmin_Tmax("V", Tmin_V, Tmax_V, ixg, data_2pts.Tag[iens]);
	Corr.Tmin= Tmin_V+shift+shift_d; Corr.Tmax= Tmax_V+shift+shift_d;
	FV_d.distr_list.push_back( Corr.Fit_distr(FV_d_distr));
	FV_d_per_kin[ixg-1].distr_list.push_back(Corr.Fit_distr(FA_u_distr));
	vec_d_Tmin.push_back(Tmin_V+shift+shift_d);
	vec_d_Tmax.push_back(Tmax_V+shift+shift_d);
      }
      
      
    
      
    }

    FA_per_ens[iens] = FA;
    FV_per_ens[iens] = FV;
    FA_u_per_ens[iens] = FA_u;
    FA_d_per_ens[iens] = FA_d;
    FV_u_per_ens[iens] = FV_u;
    FV_d_per_ens[iens] = FV_d;
    


    //Print fitted form factors
    Print_To_File({}, {xg_list.ave(), xg_list.err(), FV.ave(), FV.err(), vec_Tmin, vec_Tmax}, "../data/ph_emission/"+ph_type_mes+"/FF/"+data_2pts.Tag[iens]+"/fit_results/FV.dat", "", "#xg Dxg  FV DFV Tmin Tmax");
    Print_To_File({}, {xg_list.ave(), xg_list.err(), FA.ave(), FA.err(), ax_Tmin, ax_Tmax}, "../data/ph_emission/"+ph_type_mes+"/FF/"+data_2pts.Tag[iens]+"/fit_results/FA.dat", "", "#xg Dxg  FA DFA Tmin Tmax");

    Print_To_File({}, {xg_list.ave(), xg_list.err(), FV_u.ave(), FV_u.err(), vec_u_Tmin, vec_u_Tmax}, "../data/ph_emission/"+ph_type_mes+"/FF_u/"+data_2pts.Tag[iens]+"/fit_results/FV.dat", "", "#xg Dxg  FV DFV");
    Print_To_File({}, {xg_list.ave(), xg_list.err(), FA_u.ave(), FA_u.err(), ax_u_Tmin, ax_u_Tmax}, "../data/ph_emission/"+ph_type_mes+"/FF_u/"+data_2pts.Tag[iens]+"/fit_results/FA.dat", "", "#xg Dxg  FA DFA");

    Print_To_File({}, {xg_list.ave(), xg_list.err(), (-1.0*FV_d).ave(), (-1.0*FV_d).err(), vec_d_Tmin, vec_d_Tmax}, "../data/ph_emission/"+ph_type_mes+"/FF_d/"+data_2pts.Tag[iens]+"/fit_results/FV.dat", "", "#xg Dxg  FV DFV");
    Print_To_File({}, {xg_list.ave(), xg_list.err(), FA_d.ave(), FA_d.err(), ax_d_Tmin, ax_d_Tmax}, "../data/ph_emission/"+ph_type_mes+"/FF_d/"+data_2pts.Tag[iens]+"/fit_results/FA.dat", "", "#xg Dxg  FA DFA");


  }


  


  if(num_xg==1) crash("Exiting...No interpolation nor continuum extrapolation to perform, only xg available is zero");


  //Print all ensembles for given xg
  for(int ixg=1;ixg<num_xg;ixg++) {

    Print_To_File({}, {(a_distr_list/fmTGeV).ave(), FA_per_kin[ixg-1].ave(), FA_per_kin[ixg-1].err(), L_list}, "../data/ph_emission/"+ph_type_mes+"/FF/per_kin/FA_ixg_"+to_string(ixg)+".dat", "", "#a[fm]  FA FA_err L/a");
    Print_To_File({}, {(a_distr_list/fmTGeV).ave(), FV_per_kin[ixg-1].ave(), FV_per_kin[ixg-1].err(), L_list}, "../data/ph_emission/"+ph_type_mes+"/FF/per_kin/FV_ixg_"+to_string(ixg)+".dat", "", "#a[fm]  FA FA_err L/a");
    
  }

  
  //################################################################################################

  //interpolate form factors for each ensemble

  vector<vector<boost::math::interpolators::cardinal_cubic_b_spline<double>>> FA_interpol_jacks(Nens);
  vector<vector<boost::math::interpolators::cardinal_cubic_b_spline<double>>> FV_interpol_jacks(Nens);

  for(int iens=0;iens<Nens;iens++) {

    int num_xg_iens= FA_per_ens[iens].size();

    for(int ijack=0;ijack<Njacks;ijack++) {

      
      Vfloat FA_jacks, FV_jacks;
      for(int ixg=1;ixg<num_xg_iens;ixg++) { FA_jacks.push_back( FA_per_ens[iens].distr_list[ixg-1].distr[ijack]); FV_jacks.push_back( FV_per_ens[iens].distr_list[ixg-1].distr[ijack]);}
      FA_interpol_jacks[iens].emplace_back( FA_jacks.begin(), FA_jacks.end(), 0.1, 0.1);
      FV_interpol_jacks[iens].emplace_back( FV_jacks.begin(), FV_jacks.end(), 0.1, 0.1);
    }
  }

  auto FA_interpol_distr = [&FA_interpol_jacks, &UseJack, &Njacks](double xg, int iens) -> distr_t {

			     distr_t return_distr(UseJack);

			     for(int ijack=0; ijack<Njacks;ijack++) { return_distr.distr.push_back( FA_interpol_jacks[iens][ijack](xg));}

			     return return_distr;

			   };

  auto FV_interpol_distr = [&FV_interpol_jacks, &UseJack, &Njacks](double xg, int iens) -> distr_t {

			     distr_t return_distr(UseJack);

			     for(int ijack=0; ijack<Njacks;ijack++) { return_distr.distr.push_back( FV_interpol_jacks[iens][ijack](xg));}

			     return return_distr;

			     } ;



  for(int iens=0;iens<Nens;iens++) {

    distr_t_list FA_interpol_to_print_distr(UseJack);
    distr_t_list FV_interpol_to_print_distr(UseJack);
    for(auto &X: xg_to_spline) {
      FA_interpol_to_print_distr.distr_list.push_back( FA_interpol_distr(X, iens));
      FV_interpol_to_print_distr.distr_list.push_back( FV_interpol_distr(X, iens));
    }

    Print_To_File({}, {xg_to_spline, FA_interpol_to_print_distr.ave(), FA_interpol_to_print_distr.err()}, "../data/ph_emission/"+ph_type+"/"+Meson+"/FF/"+data_2pts.Tag[iens]+"/FA_interpol.dat", "", "#xg FA FA_err");
    Print_To_File({}, {xg_to_spline, FV_interpol_to_print_distr.ave(), FV_interpol_to_print_distr.err()}, "../data/ph_emission/"+ph_type+"/"+Meson+"/FF/"+data_2pts.Tag[iens]+"/FV_interpol.dat", "", "#xg FV FV_err");


  }


  //Print MP, FP, MP_ov_FP
  Print_To_File({}, {a_distr_list_red.ave(), MP_list_red.ave(), MP_list_red.err()}, "../data/ph_emission/"+ph_type+"/"+Meson+"/mass/masses.list", "", "#a MP MP_err");
  Print_To_File({}, {a_distr_list.ave(), FP_list.ave(), FP_list.err(), MP_ov_FP_list.ave(), MP_ov_FP_list.err()}, "../data/ph_emission/"+ph_type+"/"+Meson+"/decay_const/fP.list", "", "#a FP FP_err MP/FP MP/FP_err  ");


  //#################################################################################################



  //If ensemble B96 is present estimate FSEs
  int iens_B64=-1;
  int iens_B96=-1;
  for(int iens=0;iens<Nens;iens++) {
    if(data_2pts.Tag[iens] == "cB211b.072.64") iens_B64=iens;
    else if(data_2pts.Tag[iens] == "cB211b.072.96") iens_B96=iens;
  }

  if( (iens_B64 > -1)  && (iens_B96 > -1)) { //compute FSE

    ofstream Print_FVE("../data/ph_emission/"+ph_type+"/"+Meson+"/FF/FVE.dat");
    for(int ixg=1;ixg<num_xg;ixg++) {

      double sA_max= sqrt( pow(FA_per_ens[iens_B96].err(ixg-1),2) +pow(FA_per_ens[iens_B64].err(ixg-1),2));
      double FA_corr= fabs(FA_per_ens[iens_B96].ave(ixg-1)-FA_per_ens[iens_B64].ave(ixg-1));
      double sV_max= sqrt( pow(FV_per_ens[iens_B96].err(ixg-1),2) +pow(FV_per_ens[iens_B64].err(ixg-1),2));
      double FV_corr= fabs(FV_per_ens[iens_B96].ave(ixg-1)-FV_per_ens[iens_B64].ave(ixg-1));
      FA_corr= fabs((FA_corr/FA_per_ens[iens_B64].ave(ixg-1))*erf( FA_corr/(sqrt(2.0)*sA_max)));
      FV_corr= fabs((FV_corr/FV_per_ens[iens_B64].ave(ixg-1))*erf( FV_corr/(sqrt(2.0)*sV_max)));
      Print_FVE<<to_string_with_precision(ixg*0.1,2)<<" "<<FA_corr<<" "<<FV_corr<<" "<<(FA_corr*fabs(FA_per_ens[iens_B64].ave(ixg-1))/FA_per_ens[iens_B64].err(ixg-1))<<" "<<(FV_corr*fabs(FV_per_ens[iens_B64].ave(ixg-1))/FV_per_ens[iens_B64].err(ixg-1))<<endl;

      if(ixg==7) {
	 Print_FVE<<"0.80"<<" "<<FA_corr<<" "<<FV_corr<<" "<<(FA_corr*fabs(FA_per_ens[iens_B64].ave(ixg-1))/FA_per_ens[iens_B64].err(ixg-1))<<" "<<(FV_corr*fabs(FV_per_ens[iens_B64].ave(ixg-1))/FV_per_ens[iens_B64].err(ixg-1))<<endl;
	 Print_FVE<<"0.90"<<" "<<FA_corr<<" "<<FV_corr<<" "<<(FA_corr*fabs(FA_per_ens[iens_B64].ave(ixg-1))/FA_per_ens[iens_B64].err(ixg-1))<<" "<<(FV_corr*fabs(FV_per_ens[iens_B64].ave(ixg-1))/FV_per_ens[iens_B64].err(ixg-1))<<endl;
	 Print_FVE<<"1.00"<<" "<<FA_corr<<" "<<FV_corr<<" "<<(FA_corr*fabs(FA_per_ens[iens_B64].ave(ixg-1))/FA_per_ens[iens_B64].err(ixg-1))<<" "<<(FV_corr*fabs(FV_per_ens[iens_B64].ave(ixg-1))/FV_per_ens[iens_B64].err(ixg-1))<<endl;
      }

      
    }
   
    Print_FVE.close();



    //u contribution
    ofstream Print_FVE_u("../data/ph_emission/"+ph_type+"/"+Meson+"/FF/FVE_u.dat");
    for(int ixg=1;ixg<num_xg;ixg++) {

      double sA_max= sqrt( pow(FA_u_per_ens[iens_B96].err(ixg-1),2) +pow(FA_u_per_ens[iens_B64].err(ixg-1),2));
      double FA_corr= fabs(FA_u_per_ens[iens_B96].ave(ixg-1)-FA_u_per_ens[iens_B64].ave(ixg-1));
      double sV_max= sqrt( pow(FV_u_per_ens[iens_B96].err(ixg-1),2) +pow(FV_u_per_ens[iens_B64].err(ixg-1),2));
      double FV_corr= fabs(FV_u_per_ens[iens_B96].ave(ixg-1)-FV_u_per_ens[iens_B64].ave(ixg-1));
      FA_corr= fabs((FA_corr/FA_u_per_ens[iens_B64].ave(ixg-1))*erf( FA_corr/(sqrt(2.0)*sA_max)));
      FV_corr= fabs((FV_corr/FV_u_per_ens[iens_B64].ave(ixg-1))*erf( FV_corr/(sqrt(2.0)*sV_max)));
      Print_FVE_u<<to_string_with_precision(ixg*0.1,2)<<" "<<FA_corr<<" "<<FV_corr<<" "<<(FA_corr*fabs(FA_u_per_ens[iens_B64].ave(ixg-1))/FA_u_per_ens[iens_B64].err(ixg-1))<<" "<<(FV_corr*fabs(FV_u_per_ens[iens_B64].ave(ixg-1))/FV_u_per_ens[iens_B64].err(ixg-1))<<endl;


      if(ixg==7) {
	Print_FVE_u<<"0.80"<<" "<<FA_corr<<" "<<FV_corr<<" "<<(FA_corr*fabs(FA_u_per_ens[iens_B64].ave(ixg-1))/FA_u_per_ens[iens_B64].err(ixg-1))<<" "<<(FV_corr*fabs(FV_u_per_ens[iens_B64].ave(ixg-1))/FV_u_per_ens[iens_B64].err(ixg-1))<<endl;
	Print_FVE_u<<"0.90"<<" "<<FA_corr<<" "<<FV_corr<<" "<<(FA_corr*fabs(FA_u_per_ens[iens_B64].ave(ixg-1))/FA_u_per_ens[iens_B64].err(ixg-1))<<" "<<(FV_corr*fabs(FV_u_per_ens[iens_B64].ave(ixg-1))/FV_u_per_ens[iens_B64].err(ixg-1))<<endl;
	Print_FVE_u<<"1.00"<<" "<<FA_corr<<" "<<FV_corr<<" "<<(FA_corr*fabs(FA_u_per_ens[iens_B64].ave(ixg-1))/FA_u_per_ens[iens_B64].err(ixg-1))<<" "<<(FV_corr*fabs(FV_u_per_ens[iens_B64].ave(ixg-1))/FV_u_per_ens[iens_B64].err(ixg-1))<<endl;
      }

      
    }
    
    Print_FVE_u.close();



    //d contribution
    ofstream Print_FVE_d("../data/ph_emission/"+ph_type+"/"+Meson+"/FF/FVE_d.dat");
    for(int ixg=1;ixg<num_xg;ixg++) {

      double sA_max= sqrt( pow(FA_d_per_ens[iens_B96].err(ixg-1),2) +pow(FA_d_per_ens[iens_B64].err(ixg-1),2));
      double FA_corr= fabs(FA_d_per_ens[iens_B96].ave(ixg-1)-FA_d_per_ens[iens_B64].ave(ixg-1));
      double sV_max= sqrt( pow(FV_d_per_ens[iens_B96].err(ixg-1),2) +pow(FV_d_per_ens[iens_B64].err(ixg-1),2));
      double FV_corr= fabs(FV_d_per_ens[iens_B96].ave(ixg-1)-FV_d_per_ens[iens_B64].ave(ixg-1));
      FA_corr= fabs((FA_corr/FA_d_per_ens[iens_B64].ave(ixg-1))*erf( FA_corr/(sqrt(2.0)*sA_max)));
      FV_corr= fabs((FV_corr/FV_d_per_ens[iens_B64].ave(ixg-1))*erf( FV_corr/(sqrt(2.0)*sV_max)));
      Print_FVE_d<<to_string_with_precision(ixg*0.1,2)<<" "<<FA_corr<<" "<<FV_corr<<" "<<(FA_corr*fabs(FA_d_per_ens[iens_B64].ave(ixg-1))/FA_d_per_ens[iens_B64].err(ixg-1))<<" "<<(FV_corr*fabs(FV_d_per_ens[iens_B64].ave(ixg-1))/FV_d_per_ens[iens_B64].err(ixg-1))<<endl;


      if(ixg==7) {
	Print_FVE_d<<"0.80"<<" "<<FA_corr<<" "<<FV_corr<<" "<<(FA_corr*fabs(FA_d_per_ens[iens_B64].ave(ixg-1))/FA_d_per_ens[iens_B64].err(ixg-1))<<" "<<(FV_corr*fabs(FV_d_per_ens[iens_B64].ave(ixg-1))/FV_d_per_ens[iens_B64].err(ixg-1))<<endl;
	Print_FVE_d<<"0.90"<<" "<<FA_corr<<" "<<FV_corr<<" "<<(FA_corr*fabs(FA_d_per_ens[iens_B64].ave(ixg-1))/FA_d_per_ens[iens_B64].err(ixg-1))<<" "<<(FV_corr*fabs(FV_d_per_ens[iens_B64].ave(ixg-1))/FV_d_per_ens[iens_B64].err(ixg-1))<<endl;
	Print_FVE_d<<"1.00"<<" "<<FA_corr<<" "<<FV_corr<<" "<<(FA_corr*fabs(FA_d_per_ens[iens_B64].ave(ixg-1))/FA_d_per_ens[iens_B64].err(ixg-1))<<" "<<(FV_corr*fabs(FV_d_per_ens[iens_B64].ave(ixg-1))/FV_d_per_ens[iens_B64].err(ixg-1))<<endl;
      }

      
    }
 
    Print_FVE_d.close();



    
  }
   
  //continuum extrapolation

  
 

  if(Perform_continuum_extrapolation) {


   
  
    vector<boost::math::interpolators::cardinal_cubic_b_spline<double>> FA_cont_interpol_jacks;
    vector<boost::math::interpolators::cardinal_cubic_b_spline<double>> FV_cont_interpol_jacks;

    //u contribution
    vector<boost::math::interpolators::cardinal_cubic_b_spline<double>> FA_u_cont_interpol_jacks;
    vector<boost::math::interpolators::cardinal_cubic_b_spline<double>> FV_u_cont_interpol_jacks;

    //d contribution
    vector<boost::math::interpolators::cardinal_cubic_b_spline<double>> FA_d_cont_interpol_jacks;
    vector<boost::math::interpolators::cardinal_cubic_b_spline<double>> FV_d_cont_interpol_jacks;



    //load FVE

    ifstream Read_FVE("../data/ph_emission/"+ph_type+"/"+Meson+"/FF/FVE.dat");
    for(int iens=0;iens<Nens;iens++) {
      
      ofstream Print_F_after_FVE("../data/ph_emission/"+ph_type+"/"+Meson+"/FF/"+data_2pts.Tag[iens]+"/fit_results/FV_w_FVE.dat");
      Print_F_after_FVE<<"#xg FV DFV"<<endl;
      Print_F_after_FVE.close();
      Print_F_after_FVE.open("../data/ph_emission/"+ph_type+"/"+Meson+"/FF/"+data_2pts.Tag[iens]+"/fit_results/FA_w_FVE.dat");
      Print_F_after_FVE<<"#xg FA DFA"<<endl;
      Print_F_after_FVE.close();
    }

    Vfloat GM_FVE_V(Njacks), GM_FVE_A(Njacks);
    double sum_A=0;
    double sum_V=0;
    for(int ijack=0;ijack<Njacks;ijack++) { double rn_A=GM(); double rn_V=GM();  GM_FVE_V[ijack] = GM(); GM_FVE_A[ijack] = GM();}
    ofstream Print_GM("../data/ph_emission/"+ph_type+"/"+Meson+"/FF/FVE_jack.dat");
    for(int ijack=0;ijack<Njacks;ijack++) Print_GM<<GM_FVE_V[ijack]<<" "<<GM_FVE_A[ijack]<<endl;
    Print_GM.close();


    
    GM_FVE_V = {0.735329, -0.299985, -1.35741, -1.66874, -0.262508, -0.671406 , 0.666785, 1.3346, -0.363563, -0.461873, 1.69269, -0.376367, 0.251096, -0.056725, 0.0406796, 1.27152, 0.801355, 0.193717, 0.821933, -1.61115, -0.24308, -0.271379, 1.22812, -1.88447, 1.33931, -0.955272, -0.768544, -1.06358, 2.11891, -0.721769};
    
    GM_FVE_A = {-0.852467, 1.14618, -2.11582, -0.527966, 0.354243, 0.809818, -0.531842, -1.88762, 0.0162903, -2.18775, -0.556247, -0.0706603,  -0.402653, 0.799146, 2.29603, -2.65147, -0.417837, 1.092, 2.25395, 1.84766, -1.36045, 0.751123, 1.75308, 1.44755, -0.226214, 1.41985, -0.312688,  -0.268878, -1.46655, -0.411726};
       
    for(int ixg=1;ixg<num_xg;ixg++) {
      distr_t FVE_V(UseJack), FVE_A;
      double xx, sV, sA, srel_V, srel_A;
      Read_FVE>>xx>>sA>>sV>>srel_A>>srel_V;
      cout<<ixg<<" "<<sA<<" "<<sV<<endl;
      for(int ijack=0;ijack<Njacks;ijack++) { FVE_A.distr.push_back( sA*GM_FVE_A[ijack]/((UseJack==true)?sqrt(Njacks -1.0):1.0));   }
      for(int ijack=0;ijack<Njacks;ijack++) { FVE_V.distr.push_back( sV*GM_FVE_V[ijack]/((UseJack==true)?sqrt(Njacks -1.0):1.0));   }
    for(int iens=0;iens<Nens;iens++) {
      ofstream Print_FA_after_FVE("../data/ph_emission/"+ph_type+"/"+Meson+"/FF/"+data_2pts.Tag[iens]+"/fit_results/FA_w_FVE.dat", ofstream::app);
      ofstream Print_FV_after_FVE("../data/ph_emission/"+ph_type+"/"+Meson+"/FF/"+data_2pts.Tag[iens]+"/fit_results/FV_w_FVE.dat", ofstream::app);
      FA_per_ens[iens].distr_list[ixg-1] = FA_per_ens[iens].distr_list[ixg-1] + FA_per_ens[iens].ave(ixg-1)*FVE_A;
      FV_per_ens[iens].distr_list[ixg-1] = FV_per_ens[iens].distr_list[ixg-1] + FV_per_ens[iens].ave(ixg-1)*FVE_V;
      Print_FA_after_FVE<<ixg*0.1<<" "<<FA_per_ens[iens].ave(ixg-1)<<" "<<FA_per_ens[iens].err(ixg-1)<<"  "<<FVE_A.ave()<<" "<<FVE_A.err()<<" "<<sA<<endl;
      Print_FV_after_FVE<<ixg*0.1<<" "<<FV_per_ens[iens].ave(ixg-1)<<" "<<FV_per_ens[iens].err(ixg-1)<<"  "<<FVE_V.ave()<<" "<<FVE_V.err()<<" "<<sV<<endl;
      Print_FA_after_FVE.close();
      Print_FV_after_FVE.close();
    }
    }
    Read_FVE.close();
    



    //u contribution
    ifstream Read_FVE_u("../data/ph_emission/"+ph_type+"/"+Meson+"/FF/FVE_u.dat");
    for(int ixg=1;ixg<num_xg;ixg++) {
      distr_t FVE_V(UseJack), FVE_A;
      double xx, sV, sA, srel_V, srel_A;
      Read_FVE_u>>xx>>sA>>sV>>srel_A>>srel_V;
      for(int ijack=0;ijack<Njacks;ijack++) { FVE_A.distr.push_back( sA*GM_FVE_A[ijack]/((UseJack==true)?sqrt(Njacks -1.0):1.0));   }
      for(int ijack=0;ijack<Njacks;ijack++) { FVE_V.distr.push_back( sV*GM_FVE_V[ijack]/((UseJack==true)?sqrt(Njacks -1.0):1.0));   }
    for(int iens=0;iens<Nens;iens++) {
      FA_u_per_ens[iens].distr_list[ixg-1] = FA_u_per_ens[iens].distr_list[ixg-1] + FA_u_per_ens[iens].ave(ixg-1)*FVE_A;
      FV_u_per_ens[iens].distr_list[ixg-1] = FV_u_per_ens[iens].distr_list[ixg-1] + FV_u_per_ens[iens].ave(ixg-1)*FVE_V;
      
    }
    }
    Read_FVE_u.close();


    //d contribution
    ifstream Read_FVE_d("../data/ph_emission/"+ph_type+"/"+Meson+"/FF/FVE_d.dat");
    for(int ixg=1;ixg<num_xg;ixg++) {
      distr_t FVE_V(UseJack), FVE_A(UseJack);
      double xx, sV, sA, srel_V, srel_A;
      Read_FVE_d>>xx>>sA>>sV>>srel_A>>srel_V;
      for(int ijack=0;ijack<Njacks;ijack++) { FVE_A.distr.push_back( sA*GM_FVE_A[ijack]/((UseJack==true)?sqrt(Njacks -1.0):1.0));   }
      for(int ijack=0;ijack<Njacks;ijack++) { FVE_V.distr.push_back( sV*GM_FVE_V[ijack]/((UseJack==true)?sqrt(Njacks -1.0):1.0));   }
    for(int iens=0;iens<Nens;iens++) {
      FA_d_per_ens[iens].distr_list[ixg-1] = FA_d_per_ens[iens].distr_list[ixg-1] + FA_d_per_ens[iens].ave(ixg-1)*FVE_A;
      FV_d_per_ens[iens].distr_list[ixg-1] = FV_d_per_ens[iens].distr_list[ixg-1] + FV_d_per_ens[iens].ave(ixg-1)*FVE_V;
    }
    }
    Read_FVE_d.close();
    

  class ipar_FF_Nissa {
  public:
    ipar_FF_Nissa() : FF(0.0), FF_err(0.0) {}
    double FF, FF_err, a;
    int is;
  };
  
  class fpar_FF_Nissa {
  public:
    fpar_FF_Nissa() {}
    fpar_FF_Nissa(const Vfloat &par) {
      if((signed)par.size() != 3) crash("In class fpar_FF_Nissa  class constructor Vfloat par has size != 3");
      F0=par[0];
      D1=par[1];
      D2=par[2];
    }
    double F0, D1,D2;
  };

  
  //init bootstrap fit
  bootstrap_fit<fpar_FF_Nissa,ipar_FF_Nissa> bf_FF(Njacks);
  bf_FF.set_warmup_lev(1); //sets warmup
  bf_FF.Set_number_of_measurements(Nens);
  bf_FF.Set_verbosity(1);
  bf_FF.Add_par("F0", 0.07, 0.001);
  bf_FF.Add_par("D1", 1.0, 0.1);
  bf_FF.Add_par("D2", 1.0, 0.1);
  //fit on mean values to get ch2
  bootstrap_fit<fpar_FF_Nissa,ipar_FF_Nissa> bf_FF_ch2(1);
  bf_FF_ch2.set_warmup_lev(0); //sets warmup
  bf_FF_ch2.Set_number_of_measurements(Nens);
  bf_FF_ch2.Set_verbosity(1);
  bf_FF_ch2.Add_par("F0", 0.07, 0.001);
  bf_FF_ch2.Add_par("D1", 1.0, 0.1);
  bf_FF_ch2.Add_par("D2", 1.0, 0.1);

  //count number of dof
  int dof= Nens-2;
  if(Include_a4==false) { bf_FF.Fix_par("D2", 0.0); bf_FF_ch2.Fix_par("D2",0.0);}
  else dof-- ;

  if( Use_three_finest) dof--;


  //push_back info
  rt_fit.Nmeas=Nens-(Use_three_finest==true);
  rt_fit.Ndof= dof;
  rt_fit.Npars= rt_fit.Nmeas-dof;
  

  //ansatz
  bf_FF.ansatz=  [&Use_three_finest ](const fpar_FF_Nissa &p, const ipar_FF_Nissa &ip) {

    if( (ip.is==0) && Use_three_finest) return 0.0;
		  
    return p.F0 + p.D1*pow(ip.a*Lambda_QCD,2) + p.D2*pow(ip.a*Lambda_QCD,4);

  };

  
  bf_FF.measurement=  [&Use_three_finest ](const fpar_FF_Nissa &p, const ipar_FF_Nissa &ip) {
    
    if( (ip.is==0) && Use_three_finest) return 0.0;
    
    return ip.FF;
       
  };
  bf_FF.error=  [ ](const fpar_FF_Nissa &p, const ipar_FF_Nissa &ip) {

		 return ip.FF_err;
		 };

  bf_FF_ch2.ansatz= bf_FF.ansatz;
  bf_FF_ch2.measurement = bf_FF.measurement;
  bf_FF_ch2.error = bf_FF.error;
  
  //FIT FA for all xg
  Vfloat ch2_FA;
  distr_t_list F0_A_list(UseJack);
  distr_t_list D1_A_list(UseJack);
  distr_t_list D2_A_list(UseJack);
  for(int ixg=1;ixg<num_xg;ixg++) {
  
  //fill the data
  vector<vector<ipar_FF_Nissa>> data(Njacks);
  vector<vector<ipar_FF_Nissa>> data_ch2(1);
  //allocate space for output result
  boot_fit_data<fpar_FF_Nissa> Bt_fit;
  boot_fit_data<fpar_FF_Nissa> Bt_fit_ch2;
  for(auto &data_iboot: data) data_iboot.resize(Nens);
  for(auto &data_iboot: data_ch2) data_iboot.resize(Nens);
  for(int ijack=0;ijack<Njacks;ijack++) {
    for(int iens=0;iens<Nens;iens++) {
      data[ijack][iens].FF= FA_per_ens[iens].distr_list[ixg-1].distr[ijack];
      data[ijack][iens].FF_err= FA_per_ens[iens].err(ixg-1);
      if(data_2pts.Tag[iens] == "cA211a.12.48") { data[ijack][iens].a = a_A.distr[ijack]; data[ijack][iens].is=0; }
      else if(data_2pts.Tag[iens] == "cB211b.072.64") { data[ijack][iens].a = a_B.distr[ijack]; data[ijack][iens].is=1;}
      else if(data_2pts.Tag[iens] == "cB211b.072.96") { data[ijack][iens].a = a_B.distr[ijack]; data[ijack][iens].is=1; }
      else if(data_2pts.Tag[iens] == "cC211a.06.80")  {data[ijack][iens].a = a_C.distr[ijack]; data[ijack][iens].is=2; }
      else if(data_2pts.Tag[iens] == "cD211a.054.96")  {data[ijack][iens].a = a_D.distr[ijack]; data[ijack][iens].is=3; }
      else crash("Ens_tag: "+data_2pts.Tag[iens]+" not recognized");

      //mean values
      if(ijack==0) {
	data_ch2[ijack][iens].FF= FA_per_ens[iens].ave(ixg-1);
	data_ch2[ijack][iens].FF_err= FA_per_ens[iens].err(ixg-1);
	if(data_2pts.Tag[iens] == "cA211a.12.48") { data_ch2[ijack][iens].a = a_A.distr[ijack]; data_ch2[ijack][iens].is=0; }
	else if(data_2pts.Tag[iens] == "cB211b.072.64") { data_ch2[ijack][iens].a = a_B.distr[ijack]; data_ch2[ijack][iens].is=1;}
	else if(data_2pts.Tag[iens] == "cB211b.072.96") { data_ch2[ijack][iens].a = a_B.distr[ijack]; data_ch2[ijack][iens].is=1; }
	else if(data_2pts.Tag[iens] == "cC211a.06.80")  {data_ch2[ijack][iens].a = a_C.distr[ijack]; data_ch2[ijack][iens].is=2; }
	else if(data_2pts.Tag[iens] == "cD211a.054.96")  {data_ch2[ijack][iens].a = a_D.distr[ijack]; data_ch2[ijack][iens].is=3; }
	else crash("Ens_tag: "+data_2pts.Tag[iens]+" not recognized");
      }
    }
  }
    
  //append
  bf_FF.Append_to_input_par(data);
  bf_FF_ch2.Append_to_input_par(data_ch2);
  //fit
  cout<<"Fitting FA, xg: "<<ixg<<endl;
  Bt_fit= bf_FF.Perform_bootstrap_fit();
  Bt_fit_ch2= bf_FF_ch2.Perform_bootstrap_fit();

  
  //retrieve parameters
  distr_t F0(UseJack), D1(UseJack), D2(UseJack);
  for(int ijack=0;ijack<Njacks;ijack++) { F0.distr.push_back( Bt_fit.par[ijack].F0); D1.distr.push_back( Bt_fit.par[ijack].D1); D2.distr.push_back( Bt_fit.par[ijack].D2);}
  //push_back retrieved parameters
  F0_A_list.distr_list.push_back(F0);
  D1_A_list.distr_list.push_back(D1);
  D2_A_list.distr_list.push_back(D2);
  //push_back ch2
  ch2_FA.push_back( Bt_fit_ch2.get_ch2_ave()/dof);

 
  //print fit func
  distr_t_list FA_xg_to_print(UseJack);
  for(auto &a: a_to_print) FA_xg_to_print.distr_list.push_back( F0 + D1*pow(a*fmTGeV*Lambda_QCD,2) + D2*pow(a*fmTGeV*Lambda_QCD,4));
  Print_To_File({}, {a_to_print, FA_xg_to_print.ave(), FA_xg_to_print.err()}, "../data/ph_emission/"+ph_type+"/"+Meson+"/FF/continuum/FA_"+Fit_tag+"xg_"+to_string_with_precision(0.10*ixg,2)+".fit_func", "", "#a[fm] FA FA_err");
 
  }


  
  //FIT FV for all xg
  Vfloat ch2_FV;
  bf_FF.Set_par_val("F0", -0.1, 0.001);
  bf_FF.Set_par_val("D1", 1.0, 0.1);
  if(Include_a4) bf_FF.Set_par_val("D2", 1.0, 0.1);
  bf_FF_ch2.Set_par_val("F0", -0.1, 0.001);
  bf_FF_ch2.Set_par_val("D1", 1.0, 0.1);
  if(Include_a4) bf_FF_ch2.Set_par_val("D2", 1.0, 0.1);
  distr_t_list F0_V_list(UseJack);
  distr_t_list D1_V_list(UseJack);
  distr_t_list D2_V_list(UseJack);
  for(int ixg=1;ixg<num_xg;ixg++) {

    //fill the data
    vector<vector<ipar_FF_Nissa>> data(Njacks);
    vector<vector<ipar_FF_Nissa>> data_ch2(1);
    //allocate space for output result
    boot_fit_data<fpar_FF_Nissa> Bt_fit;
    boot_fit_data<fpar_FF_Nissa> Bt_fit_ch2;
    for(auto &data_iboot: data) data_iboot.resize(Nens);
    for(auto &data_iboot: data_ch2) data_iboot.resize(Nens);
    for(int ijack=0;ijack<Njacks;ijack++) {
      for(int iens=0;iens<Nens;iens++) {
	data[ijack][iens].FF= FV_per_ens[iens].distr_list[ixg-1].distr[ijack];
	data[ijack][iens].FF_err= FV_per_ens[iens].err(ixg-1);
	if(data_2pts.Tag[iens] == "cA211a.12.48") { data[ijack][iens].a = a_A.distr[ijack]; data[ijack][iens].is=0; }
	else if(data_2pts.Tag[iens] == "cB211b.072.64") { data[ijack][iens].a = a_B.distr[ijack]; data[ijack][iens].is=1;}
	else if(data_2pts.Tag[iens] == "cB211b.072.96") { data[ijack][iens].a = a_B.distr[ijack]; data[ijack][iens].is=1; }
	else if(data_2pts.Tag[iens] == "cC211a.06.80")  {data[ijack][iens].a = a_C.distr[ijack]; data[ijack][iens].is=2; }
	else if(data_2pts.Tag[iens] == "cD211a.054.96")  {data[ijack][iens].a = a_D.distr[ijack]; data[ijack][iens].is=3; }
	else crash("Ens_tag: "+data_2pts.Tag[iens]+" not recognized");

	if(ijack==0) {
	  	data_ch2[ijack][iens].FF= FV_per_ens[iens].ave(ixg-1);
		data_ch2[ijack][iens].FF_err= FV_per_ens[iens].err(ixg-1);
		if(data_2pts.Tag[iens] == "cA211a.12.48") { data_ch2[ijack][iens].a = a_A.distr[ijack]; data_ch2[ijack][iens].is=0; }
		else if(data_2pts.Tag[iens] == "cB211b.072.64") { data_ch2[ijack][iens].a = a_B.distr[ijack]; data_ch2[ijack][iens].is=1;}
		else if(data_2pts.Tag[iens] == "cB211b.072.96") { data_ch2[ijack][iens].a = a_B.distr[ijack]; data_ch2[ijack][iens].is=1; }
		else if(data_2pts.Tag[iens] == "cC211a.06.80")  {data_ch2[ijack][iens].a = a_C.distr[ijack]; data_ch2[ijack][iens].is=2; }
		else if(data_2pts.Tag[iens] == "cD211a.054.96")  {data_ch2[ijack][iens].a = a_D.distr[ijack]; data_ch2[ijack][iens].is=3; }
		else crash("Ens_tag: "+data_2pts.Tag[iens]+" not recognized");
	}
      }
    }
    
  //append
  bf_FF.Append_to_input_par(data);
  bf_FF_ch2.Append_to_input_par(data_ch2);
  //fit
  cout<<"Fitting FV, xg: "<<ixg<<endl;
  Bt_fit= bf_FF.Perform_bootstrap_fit();
  Bt_fit_ch2= bf_FF_ch2.Perform_bootstrap_fit();
  //retrieve parameters
  distr_t F0(UseJack), D1(UseJack), D2(UseJack);
  for(int ijack=0;ijack<Njacks;ijack++) { F0.distr.push_back( Bt_fit.par[ijack].F0); D1.distr.push_back( Bt_fit.par[ijack].D1); D2.distr.push_back( Bt_fit.par[ijack].D2);}
  //push_back retrieved parameters
  F0_V_list.distr_list.push_back(F0);
  D1_V_list.distr_list.push_back(D1);
  D2_V_list.distr_list.push_back(D2);
  //push_back ch2
  ch2_FV.push_back( Bt_fit_ch2.get_ch2_ave()/dof);

  //print fit func
  distr_t_list FV_xg_to_print(UseJack);
  for(auto &a: a_to_print) FV_xg_to_print.distr_list.push_back( F0 + D1*pow(a*fmTGeV*Lambda_QCD,2) + D2*pow(a*fmTGeV*Lambda_QCD,4));
  Print_To_File({}, {a_to_print, FV_xg_to_print.ave(), FV_xg_to_print.err()}, "../data/ph_emission/"+ph_type+"/"+Meson+"/FF/continuum/FV_"+Fit_tag+"xg_"+to_string_with_precision(0.10*ixg,2)+".fit_func", "", "#a[fm] FV FV_err");
  }

 

  //Print continuum extrapolated form factors
  Print_To_File({}, {xg_t_list, F0_A_list.ave(), F0_A_list.err(), (D1_A_list/F0_A_list).ave(), (D1_A_list/F0_A_list).err(), (D2_A_list/F0_A_list).ave(), (D2_A_list/F0_A_list).err(), ch2_FA}, "../data/ph_emission/"+ph_type+"/"+Meson+"/FF/continuum/FA_"+Fit_tag+"cont.dat", "", "#xg  F0   D1   D2  ch2/dof");
  Print_To_File({}, {xg_t_list, F0_V_list.ave(), F0_V_list.err(), (D1_V_list/F0_V_list).ave(), (D1_V_list/F0_V_list).err(), (D2_V_list/F0_V_list).ave(), (D2_V_list/F0_V_list).err(), ch2_FV}, "../data/ph_emission/"+ph_type+"/"+Meson+"/FF/continuum/FV_"+Fit_tag+"cont.dat", "", "#xg  F0   D1   D2  ch2/dof");


  //Print covariance matrix
  string Analysis_tag= ( (UseJack==true)?"jack":"boot")+_Fit_tag;
  Eigen::MatrixXd Cov_FA(num_xg-1, num_xg-1);
  Eigen::MatrixXd Cov_FV(num_xg-1, num_xg-1);
  Eigen::MatrixXd Corr_FA(num_xg-1, num_xg-1);
  Eigen::MatrixXd Corr_FV(num_xg-1, num_xg-1);

  for(int x_xg=1; x_xg<num_xg;x_xg++) {
    for(int y_xg=1; y_xg<num_xg;y_xg++) {
      Cov_FA(x_xg-1, y_xg-1) = F0_A_list.distr_list[x_xg-1]%F0_A_list.distr_list[y_xg-1];
      Cov_FV(x_xg-1, y_xg-1) = F0_V_list.distr_list[x_xg-1]%F0_V_list.distr_list[y_xg-1];
      Corr_FA(x_xg-1, y_xg-1) = Cov_FA(x_xg-1,y_xg-1)/(F0_A_list.err(x_xg-1)*F0_A_list.err(y_xg-1));
      Corr_FV(x_xg-1, y_xg-1) = Cov_FV(x_xg-1,y_xg-1)/(F0_V_list.err(x_xg-1)*F0_V_list.err(y_xg-1));
    }
  }

  //Print To File
  ofstream Print_Cov_FA("../data/ph_emission/"+ph_type+"/"+Meson+"/FF/continuum/FA_"+Analysis_tag+".cov");
  ofstream Print_Cov_FV("../data/ph_emission/"+ph_type+"/"+Meson+"/FF/continuum/FV_"+Analysis_tag+".cov");
  ofstream Print_Corr_FA("../data/ph_emission/"+ph_type+"/"+Meson+"/FF/continuum/FA_"+Analysis_tag+".corr");
  ofstream Print_Corr_FV("../data/ph_emission/"+ph_type+"/"+Meson+"/FF/continuum/FV_"+Analysis_tag+".corr");

  Print_Cov_FA<<Cov_FA<<endl;
  Print_Cov_FV<<Cov_FV<<endl;
  Print_Corr_FA<<Corr_FA<<endl;
  Print_Corr_FV<<Corr_FV<<endl;


  Print_Cov_FA.close();
  Print_Cov_FV.close();
  Print_Corr_FA.close();
  Print_Corr_FV.close();


  
      
	



  //############ Up-quark CONTRIBUTION #########################
  //FIT FA(u) for all xg
  Vfloat ch2_FA_u;
  distr_t_list F0_u_A_list(UseJack);
  distr_t_list D1_u_A_list(UseJack);
  distr_t_list D2_u_A_list(UseJack);
  for(int ixg=1;ixg<num_xg;ixg++) {
  
  //fill the data
  vector<vector<ipar_FF_Nissa>> data(Njacks);
  vector<vector<ipar_FF_Nissa>> data_ch2(1);
  //allocate space for output result
  boot_fit_data<fpar_FF_Nissa> Bt_fit;
  boot_fit_data<fpar_FF_Nissa> Bt_fit_ch2;
  for(auto &data_iboot: data) data_iboot.resize(Nens);
  for(auto &data_iboot: data_ch2) data_iboot.resize(Nens);
  for(int ijack=0;ijack<Njacks;ijack++) {
    for(int iens=0;iens<Nens;iens++) {
      data[ijack][iens].FF= FA_u_per_ens[iens].distr_list[ixg-1].distr[ijack];
      data[ijack][iens].FF_err= FA_u_per_ens[iens].err(ixg-1);
      if(data_2pts.Tag[iens] == "cA211a.12.48") { data[ijack][iens].a = a_A.distr[ijack]; data[ijack][iens].is=0; }
      else if(data_2pts.Tag[iens] == "cB211b.072.64") { data[ijack][iens].a = a_B.distr[ijack]; data[ijack][iens].is=1;}
      else if(data_2pts.Tag[iens] == "cB211b.072.96") { data[ijack][iens].a = a_B.distr[ijack]; data[ijack][iens].is=1; }
      else if(data_2pts.Tag[iens] == "cC211a.06.80")  {data[ijack][iens].a = a_C.distr[ijack]; data[ijack][iens].is=2; }
      else if(data_2pts.Tag[iens] == "cD211a.054.96")  {data[ijack][iens].a = a_D.distr[ijack]; data[ijack][iens].is=3; }
      else crash("Ens_tag: "+data_2pts.Tag[iens]+" not recognized");

      //mean values
      if(ijack==0) {
	data_ch2[ijack][iens].FF= FA_u_per_ens[iens].ave(ixg-1);
	data_ch2[ijack][iens].FF_err= FA_u_per_ens[iens].err(ixg-1);
	if(data_2pts.Tag[iens] == "cA211a.12.48") { data_ch2[ijack][iens].a = a_A.distr[ijack]; data_ch2[ijack][iens].is=0; }
	else if(data_2pts.Tag[iens] == "cB211b.072.64") { data_ch2[ijack][iens].a = a_B.distr[ijack]; data_ch2[ijack][iens].is=1;}
	else if(data_2pts.Tag[iens] == "cB211b.072.96") { data_ch2[ijack][iens].a = a_B.distr[ijack]; data_ch2[ijack][iens].is=1; }
	else if(data_2pts.Tag[iens] == "cC211a.06.80")  {data_ch2[ijack][iens].a = a_C.distr[ijack]; data_ch2[ijack][iens].is=2; }
	else if(data_2pts.Tag[iens] == "cD211a.054.96")  {data_ch2[ijack][iens].a = a_D.distr[ijack]; data_ch2[ijack][iens].is=3; }
	else crash("Ens_tag: "+data_2pts.Tag[iens]+" not recognized");
      }
    }
  }
    
  //append
  bf_FF.Append_to_input_par(data);
  bf_FF_ch2.Append_to_input_par(data_ch2);
  //fit
  cout<<"Fitting FA(u), xg: "<<ixg<<endl;
  Bt_fit= bf_FF.Perform_bootstrap_fit();
  Bt_fit_ch2= bf_FF_ch2.Perform_bootstrap_fit();

  
  //retrieve parameters
  distr_t F0(UseJack), D1(UseJack), D2(UseJack);
  for(int ijack=0;ijack<Njacks;ijack++) { F0.distr.push_back( Bt_fit.par[ijack].F0); D1.distr.push_back( Bt_fit.par[ijack].D1); D2.distr.push_back( Bt_fit.par[ijack].D2);}
  //push_back retrieved parameters
  F0_u_A_list.distr_list.push_back(F0);
  D1_u_A_list.distr_list.push_back(D1);
  D2_u_A_list.distr_list.push_back(D2);
  //push_back ch2
  ch2_FA_u.push_back( Bt_fit_ch2.get_ch2_ave()/dof);

 
  //print fit func
  distr_t_list FA_xg_to_print(UseJack);
  for(auto &a: a_to_print) FA_xg_to_print.distr_list.push_back( F0 + D1*pow(a*fmTGeV*Lambda_QCD,2) + D2*pow(a*fmTGeV*Lambda_QCD,4));
  Print_To_File({}, {a_to_print, FA_xg_to_print.ave(), FA_xg_to_print.err()}, "../data/ph_emission/"+ph_type+"/"+Meson+"/FF/continuum/FA_u_"+Fit_tag+"xg_"+to_string_with_precision(0.10*ixg,2)+".fit_func", "", "#a[fm] FA FA_err");
 
  }


  
  //FIT FV for all xg
  Vfloat ch2_FV_u;
  bf_FF.Set_par_val("F0", -0.1, 0.001);
  bf_FF.Set_par_val("D1", 1.0, 0.1);
  if(Include_a4) bf_FF.Set_par_val("D2", 1.0, 0.1);
  bf_FF_ch2.Set_par_val("F0", -0.1, 0.001);
  bf_FF_ch2.Set_par_val("D1", 1.0, 0.1);
  if(Include_a4) bf_FF_ch2.Set_par_val("D2", 1.0, 0.1);
  distr_t_list F0_u_V_list(UseJack);
  distr_t_list D1_u_V_list(UseJack);
  distr_t_list D2_u_V_list(UseJack);
  for(int ixg=1;ixg<num_xg;ixg++) {

    //fill the data
    vector<vector<ipar_FF_Nissa>> data(Njacks);
    vector<vector<ipar_FF_Nissa>> data_ch2(1);
    //allocate space for output result
    boot_fit_data<fpar_FF_Nissa> Bt_fit;
    boot_fit_data<fpar_FF_Nissa> Bt_fit_ch2;
    for(auto &data_iboot: data) data_iboot.resize(Nens);
    for(auto &data_iboot: data_ch2) data_iboot.resize(Nens);
    for(int ijack=0;ijack<Njacks;ijack++) {
      for(int iens=0;iens<Nens;iens++) {
	data[ijack][iens].FF= FV_u_per_ens[iens].distr_list[ixg-1].distr[ijack];
	data[ijack][iens].FF_err= FV_u_per_ens[iens].err(ixg-1);
	if(data_2pts.Tag[iens] == "cA211a.12.48") { data[ijack][iens].a = a_A.distr[ijack]; data[ijack][iens].is=0; }
	else if(data_2pts.Tag[iens] == "cB211b.072.64") { data[ijack][iens].a = a_B.distr[ijack]; data[ijack][iens].is=1;}
	else if(data_2pts.Tag[iens] == "cB211b.072.96") { data[ijack][iens].a = a_B.distr[ijack]; data[ijack][iens].is=1; }
	else if(data_2pts.Tag[iens] == "cC211a.06.80")  {data[ijack][iens].a = a_C.distr[ijack]; data[ijack][iens].is=2; }
	else if(data_2pts.Tag[iens] == "cD211a.054.96")  {data[ijack][iens].a = a_D.distr[ijack]; data[ijack][iens].is=3; }
	else crash("Ens_tag: "+data_2pts.Tag[iens]+" not recognized");

	if(ijack==0) {
	  	data_ch2[ijack][iens].FF= FV_u_per_ens[iens].ave(ixg-1);
		data_ch2[ijack][iens].FF_err= FV_u_per_ens[iens].err(ixg-1);
		if(data_2pts.Tag[iens] == "cA211a.12.48") { data_ch2[ijack][iens].a = a_A.distr[ijack]; data_ch2[ijack][iens].is=0; }
		else if(data_2pts.Tag[iens] == "cB211b.072.64") { data_ch2[ijack][iens].a = a_B.distr[ijack]; data_ch2[ijack][iens].is=1;}
		else if(data_2pts.Tag[iens] == "cB211b.072.96") { data_ch2[ijack][iens].a = a_B.distr[ijack]; data_ch2[ijack][iens].is=1; }
		else if(data_2pts.Tag[iens] == "cC211a.06.80")  {data_ch2[ijack][iens].a = a_C.distr[ijack]; data_ch2[ijack][iens].is=2; }
		else if(data_2pts.Tag[iens] == "cD211a.054.96")  {data_ch2[ijack][iens].a = a_D.distr[ijack]; data_ch2[ijack][iens].is=3; }
		else crash("Ens_tag: "+data_2pts.Tag[iens]+" not recognized");
	}
      }
    }
    
  //append
  bf_FF.Append_to_input_par(data);
  bf_FF_ch2.Append_to_input_par(data_ch2);
  //fit
  cout<<"Fitting FV(u), xg: "<<ixg<<endl;
  Bt_fit= bf_FF.Perform_bootstrap_fit();
  Bt_fit_ch2= bf_FF_ch2.Perform_bootstrap_fit();
  //retrieve parameters
  distr_t F0(UseJack), D1(UseJack), D2(UseJack);
  for(int ijack=0;ijack<Njacks;ijack++) { F0.distr.push_back( Bt_fit.par[ijack].F0); D1.distr.push_back( Bt_fit.par[ijack].D1); D2.distr.push_back( Bt_fit.par[ijack].D2);}
  //push_back retrieved parameters
  F0_u_V_list.distr_list.push_back(F0);
  D1_u_V_list.distr_list.push_back(D1);
  D2_u_V_list.distr_list.push_back(D2);
  //push_back ch2
  ch2_FV_u.push_back( Bt_fit_ch2.get_ch2_ave()/dof);

  //print fit func
  distr_t_list FV_xg_to_print(UseJack);
  for(auto &a: a_to_print) FV_xg_to_print.distr_list.push_back( F0 + D1*pow(a*fmTGeV*Lambda_QCD,2) + D2*pow(a*fmTGeV*Lambda_QCD,4));
  Print_To_File({}, {a_to_print, FV_xg_to_print.ave(), FV_xg_to_print.err()}, "../data/ph_emission/"+ph_type+"/"+Meson+"/FF/continuum/FV_u_"+Fit_tag+"xg_"+to_string_with_precision(0.10*ixg,2)+".fit_func", "", "#a[fm] FV FV_err");
  }

  //Print continuum extrapolated form factors
  Print_To_File({}, {xg_t_list, F0_u_A_list.ave(), F0_u_A_list.err(), (D1_u_A_list/F0_u_A_list).ave(), (D1_u_A_list/F0_u_A_list).err(), (D2_u_A_list/F0_u_A_list).ave(), (D2_u_A_list/F0_u_A_list).err(), ch2_FA_u}, "../data/ph_emission/"+ph_type+"/"+Meson+"/FF/continuum/FA_u_"+Fit_tag+"cont.dat", "", "#xg  F0   D1   D2   ch2/dof");
  Print_To_File({}, {xg_t_list, F0_u_V_list.ave(), F0_u_V_list.err(), (D1_u_V_list/F0_u_V_list).ave(), (D1_u_V_list/F0_u_V_list).err(), (D2_u_V_list/F0_u_V_list).ave(), (D2_u_V_list/F0_u_V_list).err(), ch2_FV_u}, "../data/ph_emission/"+ph_type+"/"+Meson+"/FF/continuum/FV_u_"+Fit_tag+"cont.dat", "", "#xg  F0   D1   D2   ch2/dof");
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////



  Eigen::MatrixXd Cov_FA_u(num_xg-1, num_xg-1);
  Eigen::MatrixXd Cov_FV_u(num_xg-1, num_xg-1);
  Eigen::MatrixXd Corr_FA_u(num_xg-1, num_xg-1);
  Eigen::MatrixXd Corr_FV_u(num_xg-1, num_xg-1);

  for(int x_xg=1; x_xg<num_xg;x_xg++) {
    for(int y_xg=1; y_xg<num_xg;y_xg++) {
      Cov_FA_u(x_xg-1, y_xg-1) = F0_u_A_list.distr_list[x_xg-1]%F0_u_A_list.distr_list[y_xg-1];
      Cov_FV_u(x_xg-1, y_xg-1) = F0_u_V_list.distr_list[x_xg-1]%F0_u_V_list.distr_list[y_xg-1];
      Corr_FA_u(x_xg-1, y_xg-1) = Cov_FA_u(x_xg-1,y_xg-1)/(F0_u_A_list.err(x_xg-1)*F0_u_A_list.err(y_xg-1));
      Corr_FV_u(x_xg-1, y_xg-1) = Cov_FV_u(x_xg-1,y_xg-1)/(F0_u_V_list.err(x_xg-1)*F0_u_V_list.err(y_xg-1));
    }
  }

  //Print To File
  ofstream Print_Cov_FA_u("../data/ph_emission/"+ph_type+"/"+Meson+"/FF/continuum/FA_u_"+Analysis_tag+".cov");
  ofstream Print_Cov_FV_u("../data/ph_emission/"+ph_type+"/"+Meson+"/FF/continuum/FV_u_"+Analysis_tag+".cov");
  ofstream Print_Corr_FA_u("../data/ph_emission/"+ph_type+"/"+Meson+"/FF/continuum/FA_u_"+Analysis_tag+".corr");
  ofstream Print_Corr_FV_u("../data/ph_emission/"+ph_type+"/"+Meson+"/FF/continuum/FV_u_"+Analysis_tag+".corr");

  Print_Cov_FA_u<<Cov_FA_u<<endl;
  Print_Cov_FV_u<<Cov_FV_u<<endl;
  Print_Corr_FA_u<<Corr_FA_u<<endl;
  Print_Corr_FV_u<<Corr_FV_u<<endl;


  Print_Cov_FA_u.close();
  Print_Cov_FV_u.close();
  Print_Corr_FA_u.close();
  Print_Corr_FV_u.close();





  //############ down-quark CONTRIBUTION #########################
  //FIT FA(d) for all xg
  Vfloat ch2_FA_d;
  distr_t_list F0_d_A_list(UseJack);
  distr_t_list D1_d_A_list(UseJack);
  distr_t_list D2_d_A_list(UseJack);
  for(int ixg=1;ixg<num_xg;ixg++) {

    
  //fill the data
  vector<vector<ipar_FF_Nissa>> data(Njacks);
  vector<vector<ipar_FF_Nissa>> data_ch2(1);
  //allocate space for output result
  boot_fit_data<fpar_FF_Nissa> Bt_fit;
  boot_fit_data<fpar_FF_Nissa> Bt_fit_ch2;
  for(auto &data_iboot: data) data_iboot.resize(Nens);
  for(auto &data_iboot: data_ch2) data_iboot.resize(Nens);
  for(int ijack=0;ijack<Njacks;ijack++) {
    for(int iens=0;iens<Nens;iens++) {
      data[ijack][iens].FF= FA_d_per_ens[iens].distr_list[ixg-1].distr[ijack];
      data[ijack][iens].FF_err= FA_d_per_ens[iens].err(ixg-1);
      if(data_2pts.Tag[iens] == "cA211a.12.48") { data[ijack][iens].a = a_A.distr[ijack]; data[ijack][iens].is=0; }
      else if(data_2pts.Tag[iens] == "cB211b.072.64") { data[ijack][iens].a = a_B.distr[ijack]; data[ijack][iens].is=1;}
      else if(data_2pts.Tag[iens] == "cB211b.072.96") { data[ijack][iens].a = a_B.distr[ijack]; data[ijack][iens].is=1; }
      else if(data_2pts.Tag[iens] == "cC211a.06.80")  {data[ijack][iens].a = a_C.distr[ijack]; data[ijack][iens].is=2; }
      else if(data_2pts.Tag[iens] == "cD211a.054.96")  {data[ijack][iens].a = a_D.distr[ijack]; data[ijack][iens].is=3; }
      else crash("Ens_tag: "+data_2pts.Tag[iens]+" not recognized");

      //mean values
      if(ijack==0) {
	data_ch2[ijack][iens].FF= FA_d_per_ens[iens].ave(ixg-1);
	data_ch2[ijack][iens].FF_err= FA_d_per_ens[iens].err(ixg-1);
	if(data_2pts.Tag[iens] == "cA211a.12.48") { data_ch2[ijack][iens].a = a_A.distr[ijack]; data_ch2[ijack][iens].is=0; }
	else if(data_2pts.Tag[iens] == "cB211b.072.64") { data_ch2[ijack][iens].a = a_B.distr[ijack]; data_ch2[ijack][iens].is=1;}
	else if(data_2pts.Tag[iens] == "cB211b.072.96") { data_ch2[ijack][iens].a = a_B.distr[ijack]; data_ch2[ijack][iens].is=1; }
	else if(data_2pts.Tag[iens] == "cC211a.06.80")  {data_ch2[ijack][iens].a = a_C.distr[ijack]; data_ch2[ijack][iens].is=2; }
	else if(data_2pts.Tag[iens] == "cD211a.054.96")  {data_ch2[ijack][iens].a = a_D.distr[ijack]; data_ch2[ijack][iens].is=3; }
	else crash("Ens_tag: "+data_2pts.Tag[iens]+" not recognized");
      }
    }
  }
    
  //append
  bf_FF.Append_to_input_par(data);
  bf_FF_ch2.Append_to_input_par(data_ch2);
  //fit
  cout<<"Fitting FA(u), xg: "<<ixg<<endl;
  Bt_fit= bf_FF.Perform_bootstrap_fit();
  Bt_fit_ch2= bf_FF_ch2.Perform_bootstrap_fit();

  
  //retrieve parameters
  distr_t F0(UseJack), D1(UseJack), D2(UseJack);
  for(int ijack=0;ijack<Njacks;ijack++) { F0.distr.push_back( Bt_fit.par[ijack].F0); D1.distr.push_back( Bt_fit.par[ijack].D1); D2.distr.push_back( Bt_fit.par[ijack].D2);}
  //push_back retrieved parameters
  F0_d_A_list.distr_list.push_back(F0);
  D1_d_A_list.distr_list.push_back(D1);
  D2_d_A_list.distr_list.push_back(D2);
  //push_back ch2
  ch2_FA_d.push_back( Bt_fit_ch2.get_ch2_ave()/dof);

 
  //print fit func
  distr_t_list FA_xg_to_print(UseJack);
  for(auto &a: a_to_print) FA_xg_to_print.distr_list.push_back( F0 + D1*pow(a*fmTGeV*Lambda_QCD,2) + D2*pow(a*fmTGeV*Lambda_QCD,4));
  Print_To_File({}, {a_to_print, FA_xg_to_print.ave(), FA_xg_to_print.err()}, "../data/ph_emission/"+ph_type+"/"+Meson+"/FF/continuum/FA_d_"+Fit_tag+"xg_"+to_string_with_precision(0.10*ixg,2)+".fit_func", "", "#a[fm] FA FA_err");
 
  }


  
  //FIT FV for all xg
  Vfloat ch2_FV_d;
  bf_FF.Set_par_val("F0", -0.1, 0.001);
  bf_FF.Set_par_val("D1", 1.0, 0.1);
  if(Include_a4) bf_FF.Set_par_val("D2", 1.0, 0.1);
  bf_FF_ch2.Set_par_val("F0", -0.1, 0.001);
  bf_FF_ch2.Set_par_val("D1", 1.0, 0.1);
  if(Include_a4) bf_FF_ch2.Set_par_val("D2", 1.0, 0.1);
  distr_t_list F0_d_V_list(UseJack);
  distr_t_list D1_d_V_list(UseJack);
  distr_t_list D2_d_V_list(UseJack);
  for(int ixg=1;ixg<num_xg;ixg++) {

    //fill the data
    vector<vector<ipar_FF_Nissa>> data(Njacks);
    vector<vector<ipar_FF_Nissa>> data_ch2(1);
    //allocate space for output result
    boot_fit_data<fpar_FF_Nissa> Bt_fit;
    boot_fit_data<fpar_FF_Nissa> Bt_fit_ch2;
    for(auto &data_iboot: data) data_iboot.resize(Nens);
    for(auto &data_iboot: data_ch2) data_iboot.resize(Nens);
    for(int ijack=0;ijack<Njacks;ijack++) {
      for(int iens=0;iens<Nens;iens++) {
	data[ijack][iens].FF= FV_d_per_ens[iens].distr_list[ixg-1].distr[ijack];
	data[ijack][iens].FF_err= FV_d_per_ens[iens].err(ixg-1);
	if(data_2pts.Tag[iens] == "cA211a.12.48") { data[ijack][iens].a = a_A.distr[ijack]; data[ijack][iens].is=0; }
	else if(data_2pts.Tag[iens] == "cB211b.072.64") { data[ijack][iens].a = a_B.distr[ijack]; data[ijack][iens].is=1;}
	else if(data_2pts.Tag[iens] == "cB211b.072.96") { data[ijack][iens].a = a_B.distr[ijack]; data[ijack][iens].is=1; }
	else if(data_2pts.Tag[iens] == "cC211a.06.80")  {data[ijack][iens].a = a_C.distr[ijack]; data[ijack][iens].is=2; }
	else if(data_2pts.Tag[iens] == "cD211a.054.96")  {data[ijack][iens].a = a_D.distr[ijack]; data[ijack][iens].is=3; }
	else crash("Ens_tag: "+data_2pts.Tag[iens]+" not recognized");

	if(ijack==0) {
	  	data_ch2[ijack][iens].FF= FV_d_per_ens[iens].ave(ixg-1);
		data_ch2[ijack][iens].FF_err= FV_d_per_ens[iens].err(ixg-1);
		if(data_2pts.Tag[iens] == "cA211a.12.48") { data_ch2[ijack][iens].a = a_A.distr[ijack]; data_ch2[ijack][iens].is=0; }
		else if(data_2pts.Tag[iens] == "cB211b.072.64") { data_ch2[ijack][iens].a = a_B.distr[ijack]; data_ch2[ijack][iens].is=1;}
		else if(data_2pts.Tag[iens] == "cB211b.072.96") { data_ch2[ijack][iens].a = a_B.distr[ijack]; data_ch2[ijack][iens].is=1; }
		else if(data_2pts.Tag[iens] == "cC211a.06.80")  {data_ch2[ijack][iens].a = a_C.distr[ijack]; data_ch2[ijack][iens].is=2; }
		else if(data_2pts.Tag[iens] == "cD211a.054.96")  {data_ch2[ijack][iens].a = a_D.distr[ijack]; data_ch2[ijack][iens].is=3; }
		else crash("Ens_tag: "+data_2pts.Tag[iens]+" not recognized");
	}
      }
    }
    
  //append
  bf_FF.Append_to_input_par(data);
  bf_FF_ch2.Append_to_input_par(data_ch2);
  //fit
  cout<<"Fitting FV(d), xg: "<<ixg<<endl;
  Bt_fit= bf_FF.Perform_bootstrap_fit();
  Bt_fit_ch2= bf_FF_ch2.Perform_bootstrap_fit();
  //retrieve parameters
  distr_t F0(UseJack), D1(UseJack), D2(UseJack);
  for(int ijack=0;ijack<Njacks;ijack++) { F0.distr.push_back( Bt_fit.par[ijack].F0); D1.distr.push_back( Bt_fit.par[ijack].D1); D2.distr.push_back( Bt_fit.par[ijack].D2);}
  //push_back retrieved parameters
  F0_d_V_list.distr_list.push_back(F0);
  D1_d_V_list.distr_list.push_back(D1);
  D2_d_V_list.distr_list.push_back(D2);
  //push_back ch2
  ch2_FV_d.push_back( Bt_fit_ch2.get_ch2_ave()/dof);

  //print fit func
  distr_t_list FV_xg_to_print(UseJack);
  for(auto &a: a_to_print) FV_xg_to_print.distr_list.push_back( F0 + D1*pow(a*fmTGeV*Lambda_QCD,2) + D2*pow(a*fmTGeV*Lambda_QCD,4));
  Print_To_File({}, {a_to_print, FV_xg_to_print.ave(), FV_xg_to_print.err()}, "../data/ph_emission/"+ph_type+"/"+Meson+"/FF/continuum/FV_d_"+Fit_tag+"xg_"+to_string_with_precision(0.10*ixg,2)+".fit_func", "", "#a[fm] FV FV_err");
  }

  //Print continuum extrapolated form factors
  Print_To_File({}, {xg_t_list, F0_d_A_list.ave(), F0_d_A_list.err(), (D1_d_A_list/F0_d_A_list).ave(), (D1_d_A_list/F0_d_A_list).err(), (D2_d_A_list/F0_d_A_list).ave(), (D2_d_A_list/F0_d_A_list).err(), ch2_FA_d}, "../data/ph_emission/"+ph_type+"/"+Meson+"/FF/continuum/FA_d_"+Fit_tag+"cont.dat", "", "#xg  F0   D1   D2   ch2/dof");
  Print_To_File({}, {xg_t_list, F0_d_V_list.ave(), F0_d_V_list.err(), (D1_d_V_list/F0_d_V_list).ave(), (D1_d_V_list/F0_d_V_list).err(), (D2_d_V_list/F0_d_V_list).ave(), (D2_d_V_list/F0_d_V_list).err(), ch2_FV_d}, "../data/ph_emission/"+ph_type+"/"+Meson+"/FF/continuum/FV_d_"+Fit_tag+"cont.dat", "", "#xg  F0   D1   D2   ch2/dof");
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////


  

  Eigen::MatrixXd Cov_FA_d(num_xg-1, num_xg-1);
  Eigen::MatrixXd Cov_FV_d(num_xg-1, num_xg-1);
  Eigen::MatrixXd Corr_FA_d(num_xg-1, num_xg-1);
  Eigen::MatrixXd Corr_FV_d(num_xg-1, num_xg-1);

  for(int x_xg=1; x_xg<num_xg;x_xg++) {
    for(int y_xg=1; y_xg<num_xg;y_xg++) {
      Cov_FA_d(x_xg-1, y_xg-1) = F0_d_A_list.distr_list[x_xg-1]%F0_d_A_list.distr_list[y_xg-1];
      Cov_FV_d(x_xg-1, y_xg-1) = F0_d_V_list.distr_list[x_xg-1]%F0_d_V_list.distr_list[y_xg-1];
      Corr_FA_d(x_xg-1, y_xg-1) = Cov_FA_d(x_xg-1,y_xg-1)/(F0_d_A_list.err(x_xg-1)*F0_d_A_list.err(y_xg-1));
      Corr_FV_d(x_xg-1, y_xg-1) = Cov_FV_d(x_xg-1,y_xg-1)/(F0_d_V_list.err(x_xg-1)*F0_d_V_list.err(y_xg-1));
    }
  }

  //Print To File
  ofstream Print_Cov_FA_d("../data/ph_emission/"+ph_type+"/"+Meson+"/FF/continuum/FA_d_"+Analysis_tag+".cov");
  ofstream Print_Cov_FV_d("../data/ph_emission/"+ph_type+"/"+Meson+"/FF/continuum/FV_d_"+Analysis_tag+".cov");
  ofstream Print_Corr_FA_d("../data/ph_emission/"+ph_type+"/"+Meson+"/FF/continuum/FA_d_"+Analysis_tag+".corr");
  ofstream Print_Corr_FV_d("../data/ph_emission/"+ph_type+"/"+Meson+"/FF/continuum/FV_d_"+Analysis_tag+".corr");

  Print_Cov_FA_d<<Cov_FA_d<<endl;
  Print_Cov_FV_d<<Cov_FV_d<<endl;
  Print_Corr_FA_d<<Corr_FA_d<<endl;
  Print_Corr_FV_d<<Corr_FV_d<<endl;


  Print_Cov_FA_d.close();
  Print_Cov_FV_d.close();
  Print_Corr_FA_d.close();
  Print_Corr_FV_d.close();




  //##########################################

 
  //push_back FIT INFO
  rt_fit.Ch2_FA= ch2_FA;
  rt_fit.Ch2_FV= ch2_FV;
  rt_fit.Ch2_FA_u = ch2_FA_u;
  rt_fit.Ch2_FV_u = ch2_FV_u;
  rt_fit.Ch2_FA_d = ch2_FA_d;
  rt_fit.Ch2_FV_d = ch2_FV_d;
  rt_fit.FA = F0_A_list;
  rt_fit.FV = F0_V_list;
  rt_fit.FA_u = F0_u_A_list;
  rt_fit.FV_u = F0_u_V_list;
  rt_fit.FA_d = F0_d_A_list;
  rt_fit.FV_d = F0_d_V_list;


  //##########################################



  //VMD-like fit for FV and FA
  //###############################################################################################################

  class ipar_VMD {
  public:
    ipar_VMD() : FF(0.0), FF_err(0.0) {}
    double FF, FF_err, xg;
  };
  
  class fpar_VMD {
  public:
    fpar_VMD() {}
    fpar_VMD(const Vfloat &par) {
      if((signed)par.size() != 2) crash("In class fpar_VMD  class constructor Vfloat par has size != 2");
      A=par[0];
      M=par[1];
    }
    double A, M;
  };

  
  //init bootstrap fit
  bootstrap_fit<fpar_VMD,ipar_VMD> bf_VMD(Njacks);
  bf_VMD.set_warmup_lev(0); //sets warmup
  bf_VMD.Set_number_of_measurements(xg_t_list.size());
  bf_VMD.Set_verbosity(1);
  bf_VMD.Add_par("A", 0.07, 0.001);
  bf_VMD.Add_par("M", 1.3, 0.1);
  //fit on mean values to get ch2
  bootstrap_fit<fpar_VMD,ipar_VMD> bf_VMD_ch2(1);
  bf_VMD_ch2.set_warmup_lev(0); //sets warmup
  bf_VMD_ch2.Set_number_of_measurements(xg_t_list.size());
  bf_VMD_ch2.Set_verbosity(1);
  bf_VMD_ch2.Add_par("A", 0.07, 0.001);
  bf_VMD_ch2.Add_par("M", 1.3, 0.1);
  //bf_VMD.Fix_par("M", 1.25);
  //bf_VMD_ch2.Fix_par("M", 1.25);


  //ansatz
  bf_VMD.ansatz=  [ ](const fpar_VMD &p, const ipar_VMD &ip) {

		    double E_res= sqrt( p.M*p.M + ip.xg*ip.xg/4.0);
		    return p.A/( E_res*( E_res - (1.0 - ip.xg/2.0)) );
		 };
  bf_VMD.measurement=  [ ](const fpar_VMD &p, const ipar_VMD &ip) {

		 return ip.FF;
		 };
  bf_VMD.error=  [ ](const fpar_VMD &p, const ipar_VMD &ip) {

		 return ip.FF_err;
		 };

  bf_VMD_ch2.ansatz= bf_VMD.ansatz;
  bf_VMD_ch2.measurement = bf_VMD.measurement;
  bf_VMD_ch2.error = bf_VMD.error;


  //start fitting FA
  //fill the data
  vector<vector<ipar_VMD>> data_VMD(Njacks);
  vector<vector<ipar_VMD>> data_VMD_ch2(1);
  //allocate space for output result
  boot_fit_data<fpar_VMD> Bt_fit_FA_VMD;
  boot_fit_data<fpar_VMD> Bt_fit_FA_VMD_ch2;
  for(auto &data_iboot: data_VMD) data_iboot.resize(xg_t_list.size());
  for(auto &data_iboot: data_VMD_ch2) data_iboot.resize(xg_t_list.size());
  for(int ijack=0;ijack<Njacks;ijack++) {
    for(int ix=0;ix<(signed)xg_t_list.size();ix++) {
      data_VMD[ijack][ix].FF = F0_A_list.distr_list[ix].distr[ijack];
      data_VMD[ijack][ix].FF_err= F0_A_list.err(ix);
      data_VMD[ijack][ix].xg= xg_t_list[ix];
      if(ijack==0) {
	data_VMD_ch2[ijack][ix].FF = F0_A_list.ave(ix);
	data_VMD_ch2[ijack][ix].FF_err= F0_A_list.err(ix);
	data_VMD_ch2[ijack][ix].xg= xg_t_list[ix];

      }
    }
  }

  //append
  bf_VMD.Append_to_input_par(data_VMD);
  bf_VMD_ch2.Append_to_input_par(data_VMD_ch2);
  //fit
  cout<<"Fitting FA using VMD ansatz"<<endl;
  Bt_fit_FA_VMD= bf_VMD.Perform_bootstrap_fit();
  Bt_fit_FA_VMD_ch2= bf_VMD_ch2.Perform_bootstrap_fit();
  double ch2_red_FA_VMD= Bt_fit_FA_VMD_ch2.get_ch2_ave()/( xg_t_list.size() -2.0);

  //retrieve params
  distr_t Ampl_FA(UseJack), pole_FA(UseJack);
  for(int ijack=0;ijack<Njacks;ijack++) { Ampl_FA.distr.push_back( Bt_fit_FA_VMD.par[ijack].A); pole_FA.distr.push_back( Bt_fit_FA_VMD.par[ijack].M);}


  //start fitting FV

  bf_VMD.Set_par_val("A", -0.07, 0.001);
  bf_VMD_ch2.Set_par_val("A", -0.07, 0.001);
  //bf_VMD.Fix_par("M", 1.073);
  //bf_VMD_ch2.Fix_par("M", 1.073);

  //allocate space for output result
  boot_fit_data<fpar_VMD> Bt_fit_FV_VMD;
  boot_fit_data<fpar_VMD> Bt_fit_FV_VMD_ch2;

  for(int ijack=0;ijack<Njacks;ijack++) {
    for(int ix=0;ix<(signed)xg_t_list.size();ix++) {
      data_VMD[ijack][ix].FF = F0_V_list.distr_list[ix].distr[ijack];
      data_VMD[ijack][ix].FF_err= F0_V_list.err(ix);
      data_VMD[ijack][ix].xg= xg_t_list[ix];
      if(ijack==0) {
	data_VMD_ch2[ijack][ix].FF = F0_V_list.ave(ix);
	data_VMD_ch2[ijack][ix].FF_err= F0_V_list.err(ix);
	data_VMD_ch2[ijack][ix].xg= xg_t_list[ix];

      }
    }
  }

  //append
  bf_VMD.Append_to_input_par(data_VMD);
  bf_VMD_ch2.Append_to_input_par(data_VMD_ch2);
  //fit
  cout<<"Fitting FV using VMD ansatz"<<endl;
  Bt_fit_FV_VMD= bf_VMD.Perform_bootstrap_fit();
  Bt_fit_FV_VMD_ch2= bf_VMD_ch2.Perform_bootstrap_fit();

  double ch2_red_FV_VMD= Bt_fit_FV_VMD_ch2.get_ch2_ave()/( xg_t_list.size() -2.0);
 

  //retrieve params
  distr_t Ampl_FV(UseJack), pole_FV(UseJack);
  for(int ijack=0;ijack<Njacks;ijack++) { Ampl_FV.distr.push_back( Bt_fit_FV_VMD.par[ijack].A); pole_FV.distr.push_back( Bt_fit_FV_VMD.par[ijack].M);}

  
  distr_t_list FA_VMD_fit(UseJack), FV_VMD_fit(UseJack);
  //plot fit function
   for(auto &X: xg_to_spline_VMD) {
     ipar_VMD pp_VMD;
     pp_VMD.xg=X;
     distr_t FA_VMD_xg(UseJack), FV_VMD_xg(UseJack);
     for(int ijack=0;ijack<Njacks;ijack++) {
       FA_VMD_xg.distr.push_back( bf_VMD.ansatz( Bt_fit_FA_VMD.par[ijack], pp_VMD));
       FV_VMD_xg.distr.push_back( bf_VMD.ansatz( Bt_fit_FV_VMD.par[ijack], pp_VMD));
     }

     FA_VMD_fit.distr_list.push_back( FA_VMD_xg);
     FV_VMD_fit.distr_list.push_back( FV_VMD_xg);
  
   }

   //print

   string header_FA= "Ampl: "+to_string_with_precision(Ampl_FA.ave(),5)+" +- "+to_string_with_precision(Ampl_FA.err(),5)+" M^res/Mp: "+to_string_with_precision(pole_FA.ave(), 5)+" +- "+to_string_with_precision(pole_FA.err(), 5)+" ch2/dof: "+to_string_with_precision(ch2_red_FA_VMD  ,5)+" dof: "+to_string(xg_t_list.size()-2);
   string header_FV= "Ampl: "+to_string_with_precision(Ampl_FV.ave(),5)+" +- "+to_string_with_precision(Ampl_FV.err(),5)+" M^res/Mp: "+to_string_with_precision(pole_FV.ave(), 5)+" +- "+to_string_with_precision(pole_FV.err(), 5)+" ch2/dof: "+to_string_with_precision(ch2_red_FV_VMD  ,5)+" dof: "+to_string(xg_t_list.size()-2);
   Print_To_File({},{ xg_to_spline_VMD, FA_VMD_fit.ave(), FA_VMD_fit.err()} , "../data/ph_emission/"+ph_type+"/"+Meson+"/FF/continuum/FA_"+Fit_tag+"VMD.fit" , "", header_FA);
   Print_To_File({},{ xg_to_spline_VMD, FV_VMD_fit.ave(), FV_VMD_fit.err()} , "../data/ph_emission/"+ph_type+"/"+Meson+"/FF/continuum/FV_"+Fit_tag+"VMD.fit" , "", header_FV);







  //#################################################################################################################

  
  
  
  
  
  //interpolate continuum extrapolated form factors
  for(int ijack=0;ijack<Njacks;ijack++) {

    Vfloat FV_cont_ij, FA_cont_ij;
    Vfloat FV_u_cont_ij, FA_u_cont_ij;
    Vfloat FV_d_cont_ij, FA_d_cont_ij;
    for(int ixg=1;ixg<num_xg;ixg++) {
      FA_cont_ij.push_back( F0_A_list.distr_list[ixg-1].distr[ijack]);
      FV_cont_ij.push_back( F0_V_list.distr_list[ixg-1].distr[ijack]);

      //u contribution
      FA_u_cont_ij.push_back( F0_u_A_list.distr_list[ixg-1].distr[ijack]);
      FV_u_cont_ij.push_back( F0_u_V_list.distr_list[ixg-1].distr[ijack]);

      //d contribution
      FA_d_cont_ij.push_back( F0_d_A_list.distr_list[ixg-1].distr[ijack]);
      FV_d_cont_ij.push_back( F0_d_V_list.distr_list[ixg-1].distr[ijack]);
    }

    FA_cont_interpol_jacks.emplace_back( FA_cont_ij.begin(), FA_cont_ij.end(), 0.1, 0.1);
    FV_cont_interpol_jacks.emplace_back( FV_cont_ij.begin(), FV_cont_ij.end(), 0.1, 0.1);

    //u contribution
    FA_u_cont_interpol_jacks.emplace_back( FA_u_cont_ij.begin(), FA_u_cont_ij.end(), 0.1, 0.1);
    FV_u_cont_interpol_jacks.emplace_back( FV_u_cont_ij.begin(), FV_u_cont_ij.end(), 0.1, 0.1);

    //d contribution
    FA_d_cont_interpol_jacks.emplace_back( FA_d_cont_ij.begin(), FA_d_cont_ij.end(), 0.1, 0.1);
    FV_d_cont_interpol_jacks.emplace_back( FV_d_cont_ij.begin(), FV_d_cont_ij.end(), 0.1, 0.1);

  }


  auto FA_cont_interpol_distr= [&FA_cont_interpol_jacks, &UseJack, &Njacks](double xg) -> distr_t {
    distr_t return_distr(UseJack);
    for(int ijack=0; ijack<Njacks;ijack++) { return_distr.distr.push_back( FA_cont_interpol_jacks[ijack](xg));}
    return return_distr;
  };
  auto FV_cont_interpol_distr= [&FV_cont_interpol_jacks,  &UseJack, &Njacks](double xg) -> distr_t {
    distr_t return_distr(UseJack);
    for(int ijack=0; ijack<Njacks;ijack++) { return_distr.distr.push_back( FV_cont_interpol_jacks[ijack](xg));}
    return return_distr;
  };


  //u contribution
  auto FA_u_cont_interpol_distr= [&FA_u_cont_interpol_jacks, &UseJack, &Njacks](double xg) -> distr_t {
    distr_t return_distr(UseJack);
    for(int ijack=0; ijack<Njacks;ijack++) { return_distr.distr.push_back( FA_u_cont_interpol_jacks[ijack](xg));}
    return return_distr;
  };
  auto FV_u_cont_interpol_distr= [&FV_u_cont_interpol_jacks, &UseJack, &Njacks](double xg) -> distr_t {
    distr_t return_distr(UseJack);
    for(int ijack=0; ijack<Njacks;ijack++) { return_distr.distr.push_back( FV_u_cont_interpol_jacks[ijack](xg));}
    return return_distr;
  };
  

  //d contribution
  auto FA_d_cont_interpol_distr= [&FA_d_cont_interpol_jacks, &UseJack, &Njacks](double xg) -> distr_t {
    distr_t return_distr(UseJack);
    for(int ijack=0; ijack<Njacks;ijack++) { return_distr.distr.push_back( FA_d_cont_interpol_jacks[ijack](xg));}
    return return_distr;
  };
  auto FV_d_cont_interpol_distr= [&FV_d_cont_interpol_jacks, &UseJack, &Njacks](double xg) -> distr_t {
    distr_t return_distr(UseJack);
    for(int ijack=0; ijack<Njacks;ijack++) { return_distr.distr.push_back( FV_d_cont_interpol_jacks[ijack](xg));}
    return return_distr;
  };

 
  distr_t_list FA_interpol_cont_to_print_distr(UseJack);
  distr_t_list FV_interpol_cont_to_print_distr(UseJack);

  //u contribution
  distr_t_list FA_u_interpol_cont_to_print_distr(UseJack);
  distr_t_list FV_u_interpol_cont_to_print_distr(UseJack);

  //d contribution
  distr_t_list FA_d_interpol_cont_to_print_distr(UseJack);
  distr_t_list FV_d_interpol_cont_to_print_distr(UseJack);
  
  for(auto &X: xg_to_spline) {
    FA_interpol_cont_to_print_distr.distr_list.push_back( FA_cont_interpol_distr(X));
    FV_interpol_cont_to_print_distr.distr_list.push_back( FV_cont_interpol_distr(X));

    //u contribution
    FA_u_interpol_cont_to_print_distr.distr_list.push_back( FA_u_cont_interpol_distr(X));
    FV_u_interpol_cont_to_print_distr.distr_list.push_back( FV_u_cont_interpol_distr(X));

    //d contribution
    FA_d_interpol_cont_to_print_distr.distr_list.push_back( FA_d_cont_interpol_distr(X));
    FV_d_interpol_cont_to_print_distr.distr_list.push_back( FV_d_cont_interpol_distr(X));
  }

  Print_To_File({}, {xg_to_spline, FA_interpol_cont_to_print_distr.ave(), FA_interpol_cont_to_print_distr.err()}, "../data/ph_emission/"+ph_type+"/"+Meson+"/FF/continuum/FA_"+Fit_tag+"interpol.dat", "", "#xg FA FA_err");
  Print_To_File({}, {xg_to_spline, FV_interpol_cont_to_print_distr.ave(), FV_interpol_cont_to_print_distr.err()}, "../data/ph_emission/"+ph_type+"/"+Meson+"/FF/continuum/FV_"+Fit_tag+"interpol.dat", "", "#xg FV FV_err");

  //u contribution
  Print_To_File({}, {xg_to_spline, FA_u_interpol_cont_to_print_distr.ave(), FA_u_interpol_cont_to_print_distr.err()}, "../data/ph_emission/"+ph_type+"/"+Meson+"/FF/continuum/FA_u_"+Fit_tag+"interpol.dat", "", "#xg FA FA_err");
  Print_To_File({}, {xg_to_spline, FV_u_interpol_cont_to_print_distr.ave(), FV_u_interpol_cont_to_print_distr.err()}, "../data/ph_emission/"+ph_type+"/"+Meson+"/FF/continuum/FV_u_"+Fit_tag+"interpol.dat", "", "#xg FV FV_err");

  //d contribution
  Print_To_File({}, {xg_to_spline, FA_d_interpol_cont_to_print_distr.ave(), FA_d_interpol_cont_to_print_distr.err()}, "../data/ph_emission/"+ph_type+"/"+Meson+"/FF/continuum/FA_d_"+Fit_tag+"interpol.dat", "", "#xg FA FA_err");
  Print_To_File({}, {xg_to_spline, FV_d_interpol_cont_to_print_distr.ave(), FV_d_interpol_cont_to_print_distr.err()}, "../data/ph_emission/"+ph_type+"/"+Meson+"/FF/continuum/FV_d_"+Fit_tag+"interpol.dat", "", "#xg FV FV_err");


  //#################################################################################






  //##########################################################################################
  //##########################################################################################

  //perform continuum extrapolation of MP, fP, and MP/fP
  //MP
  bf_FF.Set_par_val("F0", MP_list.ave(0), MP_list.ave(0)/10);
  bf_FF.Set_par_val("D1", 1.0, 0.1);
  if(Include_a4) bf_FF.Set_par_val("D2", 1.0, 0.1);
  bf_FF_ch2.Set_par_val("F0", MP_list.ave(0), MP_list.ave(0)/10);
  bf_FF_ch2.Set_par_val("D1", 1.0, 0.1);
  if(Include_a4) bf_FF_ch2.Set_par_val("D2", 1.0, 0.1);
  double ch2_MP;
  distr_t F0_MP(UseJack);
  distr_t D1_MP(UseJack);
  distr_t D2_MP(UseJack);
  vector<vector<ipar_FF_Nissa>> data_MP(Njacks);
  vector<vector<ipar_FF_Nissa>> data_MP_ch2(1);
  //allocate space for output result
  boot_fit_data<fpar_FF_Nissa> Bt_fit_MP;
  boot_fit_data<fpar_FF_Nissa> Bt_fit_MP_ch2;
  for(auto &data_iboot: data_MP) data_iboot.resize(Nens);
  for(auto &data_iboot: data_MP_ch2) data_iboot.resize(Nens);
  for(int ijack=0;ijack<Njacks;ijack++) {
    for(int iens=0;iens<Nens;iens++) {
      data_MP[ijack][iens].FF= MP_list.distr_list[iens].distr[ijack];
      data_MP[ijack][iens].FF_err= MP_list.err(iens);
      if(data_2pts.Tag[iens] == "cA211a.12.48") { data_MP[ijack][iens].a = a_A.distr[ijack]; data_MP[ijack][iens].is=0; }
      else if(data_2pts.Tag[iens] == "cB211b.072.64") { data_MP[ijack][iens].a = a_B.distr[ijack]; data_MP[ijack][iens].is=1;}
      else if(data_2pts.Tag[iens] == "cB211b.072.96") { data_MP[ijack][iens].a = a_B.distr[ijack]; data_MP[ijack][iens].is=1; }
      else if(data_2pts.Tag[iens] == "cC211a.06.80")  {data_MP[ijack][iens].a = a_C.distr[ijack]; data_MP[ijack][iens].is=2; }
      else if(data_2pts.Tag[iens] == "cD211a.054.96")  {data_MP[ijack][iens].a = a_D.distr[ijack]; data_MP[ijack][iens].is=3; }
      else crash("Ens_tag: "+data_2pts.Tag[iens]+" not recognized");

      if(ijack==0) {
	data_MP_ch2[ijack][iens].FF= MP_list.ave(iens);
	data_MP_ch2[ijack][iens].FF_err= MP_list.err(iens);
	if(data_2pts.Tag[iens] == "cA211a.12.48") { data_MP_ch2[ijack][iens].a = a_A.distr[ijack]; data_MP_ch2[ijack][iens].is=0; }
	else if(data_2pts.Tag[iens] == "cB211b.072.64") { data_MP_ch2[ijack][iens].a = a_B.distr[ijack]; data_MP_ch2[ijack][iens].is=1;}
	else if(data_2pts.Tag[iens] == "cB211b.072.96") { data_MP_ch2[ijack][iens].a = a_B.distr[ijack]; data_MP_ch2[ijack][iens].is=1; }
	else if(data_2pts.Tag[iens] == "cC211a.06.80")  {data_MP_ch2[ijack][iens].a = a_C.distr[ijack]; data_MP_ch2[ijack][iens].is=2; }
	else if(data_2pts.Tag[iens] == "cD211a.054.96")  {data_MP_ch2[ijack][iens].a = a_D.distr[ijack]; data_MP_ch2[ijack][iens].is=3; }
	else crash("Ens_tag: "+data_2pts.Tag[iens]+" not recognized");
	
      }
    }
  }
  //append
  bf_FF.Append_to_input_par(data_MP);
  bf_FF_ch2.Append_to_input_par(data_MP_ch2);
  //fit
  cout<<"Fitting MP"<<endl;
  Bt_fit_MP= bf_FF.Perform_bootstrap_fit();
  Bt_fit_MP_ch2 = bf_FF_ch2.Perform_bootstrap_fit();
  //retrieve parameters
  for(int ijack=0;ijack<Njacks;ijack++) { F0_MP.distr.push_back( Bt_fit_MP.par[ijack].F0 ); D1_MP.distr.push_back( Bt_fit_MP.par[ijack].D1 ); D2_MP.distr.push_back( Bt_fit_MP.par[ijack].D2 );}
  //get ch2
  ch2_MP= Bt_fit_MP_ch2.get_ch2_ave();


  //FP
  double ch2_FP;
  bf_FF.Set_par_val("F0", FP_list.ave(0), FP_list.ave(0)/10);
  bf_FF.Set_par_val("D1", 1.0, 0.1);
  if(Include_a4) bf_FF.Set_par_val("D2", 1.0, 0.1);
  bf_FF_ch2.Set_par_val("F0", FP_list.ave(0), FP_list.ave(0)/10);
  bf_FF_ch2.Set_par_val("D1", 1.0, 0.1);
  if(Include_a4) bf_FF_ch2.Set_par_val("D2", 1.0, 0.1);
  distr_t F0_FP(UseJack);
  distr_t D1_FP(UseJack);
  distr_t D2_FP(UseJack);
  vector<vector<ipar_FF_Nissa>> data_FP(Njacks);
  vector<vector<ipar_FF_Nissa>> data_FP_ch2(1);
  //allocate space for output result
  boot_fit_data<fpar_FF_Nissa> Bt_fit_FP;
  boot_fit_data<fpar_FF_Nissa> Bt_fit_FP_ch2;
  for(auto &data_iboot: data_FP) data_iboot.resize(Nens);
  for(auto &data_iboot: data_FP_ch2) data_iboot.resize(Nens);
  for(int ijack=0;ijack<Njacks;ijack++) {
    for(int iens=0;iens<Nens;iens++) {
      data_FP[ijack][iens].FF= FP_list.distr_list[iens].distr[ijack];
      data_FP[ijack][iens].FF_err= FP_list.err(iens);
      if(data_2pts.Tag[iens] == "cA211a.12.48") { data_FP[ijack][iens].a = a_A.distr[ijack]; data_FP[ijack][iens].is=0; }
      else if(data_2pts.Tag[iens] == "cB211b.072.64") { data_FP[ijack][iens].a = a_B.distr[ijack]; data_FP[ijack][iens].is=1;}
      else if(data_2pts.Tag[iens] == "cB211b.072.96") { data_FP[ijack][iens].a = a_B.distr[ijack]; data_FP[ijack][iens].is=1; }
      else if(data_2pts.Tag[iens] == "cC211a.06.80")  {data_FP[ijack][iens].a = a_C.distr[ijack]; data_FP[ijack][iens].is=2; }
      else if(data_2pts.Tag[iens] == "cD211a.054.96")  {data_FP[ijack][iens].a = a_D.distr[ijack]; data_FP[ijack][iens].is=3; }
      else crash("Ens_tag: "+data_2pts.Tag[iens]+" not recognized");

      if(ijack==0) {
	data_FP_ch2[ijack][iens].FF= FP_list.ave(iens);
	data_FP_ch2[ijack][iens].FF_err= FP_list.err(iens);
	if(data_2pts.Tag[iens] == "cA211a.12.48") { data_FP_ch2[ijack][iens].a = a_A.distr[ijack]; data_FP_ch2[ijack][iens].is=0; }
	else if(data_2pts.Tag[iens] == "cB211b.072.64") { data_FP_ch2[ijack][iens].a = a_B.distr[ijack]; data_FP_ch2[ijack][iens].is=1;}
	else if(data_2pts.Tag[iens] == "cB211b.072.96") { data_FP_ch2[ijack][iens].a = a_B.distr[ijack]; data_FP_ch2[ijack][iens].is=1; }
	else if(data_2pts.Tag[iens] == "cC211a.06.80")  {data_FP_ch2[ijack][iens].a = a_C.distr[ijack]; data_FP_ch2[ijack][iens].is=2; }
	else if(data_2pts.Tag[iens] == "cD211a.054.96")  {data_FP_ch2[ijack][iens].a = a_D.distr[ijack]; data_FP_ch2[ijack][iens].is=3; }
	else crash("Ens_tag: "+data_2pts.Tag[iens]+" not recognized");
	
      }
    }
  }
  //append
  bf_FF.Append_to_input_par(data_FP);
  bf_FF_ch2.Append_to_input_par(data_FP_ch2);
  //fit
  cout<<"Fitting FP"<<endl;
  Bt_fit_FP= bf_FF.Perform_bootstrap_fit();
  Bt_fit_FP_ch2= bf_FF_ch2.Perform_bootstrap_fit();
  //retrieve parameters
  for(int ijack=0;ijack<Njacks;ijack++) { F0_FP.distr.push_back( Bt_fit_FP.par[ijack].F0); D1_FP.distr.push_back( Bt_fit_FP.par[ijack].D1); D2_FP.distr.push_back( Bt_fit_FP.par[ijack].D2);}
  //get ch2
  ch2_FP = Bt_fit_FP_ch2.get_ch2_ave();



  //MP/FP
  double ch2_MP_ov_FP;
  bf_FF.Set_par_val("F0", MP_ov_FP_list.ave(0), MP_ov_FP_list.ave(0)/10);
  bf_FF.Set_par_val("D1", 1.0, 0.1);
  if(Include_a4) bf_FF.Set_par_val("D2", 1.0, 0.1);
  bf_FF_ch2.Set_par_val("F0", MP_ov_FP_list.ave(0), MP_ov_FP_list.ave(0)/10);
  bf_FF_ch2.Set_par_val("D1", 1.0, 0.1);
  if(Include_a4) bf_FF_ch2.Set_par_val("D2", 1.0, 0.1);
  distr_t F0_MP_ov_FP(UseJack);
  distr_t D1_MP_ov_FP(UseJack);
  distr_t D2_MP_ov_FP(UseJack);
  vector<vector<ipar_FF_Nissa>> data_MP_ov_FP(Njacks);
  vector<vector<ipar_FF_Nissa>> data_MP_ov_FP_ch2(1);
  //allocate space for output result
  boot_fit_data<fpar_FF_Nissa> Bt_fit_MP_ov_FP;
  boot_fit_data<fpar_FF_Nissa> Bt_fit_MP_ov_FP_ch2;
  for(auto &data_iboot: data_MP_ov_FP) data_iboot.resize(Nens);
   for(auto &data_iboot: data_MP_ov_FP_ch2) data_iboot.resize(Nens);
  for(int ijack=0;ijack<Njacks;ijack++) {
    for(int iens=0;iens<Nens;iens++) {
      data_MP_ov_FP[ijack][iens].FF= MP_ov_FP_list.distr_list[iens].distr[ijack];
      data_MP_ov_FP[ijack][iens].FF_err= MP_ov_FP_list.err(iens);
      if(data_2pts.Tag[iens] == "cA211a.12.48") { data_MP_ov_FP[ijack][iens].a = a_A.distr[ijack]; data_MP_ov_FP[ijack][iens].is=0; }
      else if(data_2pts.Tag[iens] == "cB211b.072.64") { data_MP_ov_FP[ijack][iens].a = a_B.distr[ijack]; data_MP_ov_FP[ijack][iens].is=1;}
      else if(data_2pts.Tag[iens] == "cB211b.072.96") { data_MP_ov_FP[ijack][iens].a = a_B.distr[ijack]; data_MP_ov_FP[ijack][iens].is=1; }
      else if(data_2pts.Tag[iens] == "cC211a.06.80")  {data_MP_ov_FP[ijack][iens].a = a_C.distr[ijack]; data_MP_ov_FP[ijack][iens].is=2; }
      else if(data_2pts.Tag[iens] == "cD211a.054.96")  {data_MP_ov_FP[ijack][iens].a = a_D.distr[ijack]; data_MP_ov_FP[ijack][iens].is=3; }
      else crash("Ens_tag: "+data_2pts.Tag[iens]+" not recognized");

      if(ijack==0) {
	data_MP_ov_FP_ch2[ijack][iens].FF= MP_ov_FP_list.ave(iens);
	data_MP_ov_FP_ch2[ijack][iens].FF_err= MP_ov_FP_list.err(iens);
	if(data_2pts.Tag[iens] == "cA211a.12.48") { data_MP_ov_FP_ch2[ijack][iens].a = a_A.distr[ijack]; data_MP_ov_FP_ch2[ijack][iens].is=0; }
	else if(data_2pts.Tag[iens] == "cB211b.072.64") { data_MP_ov_FP_ch2[ijack][iens].a = a_B.distr[ijack]; data_MP_ov_FP_ch2[ijack][iens].is=1;}
	else if(data_2pts.Tag[iens] == "cB211b.072.96") { data_MP_ov_FP_ch2[ijack][iens].a = a_B.distr[ijack]; data_MP_ov_FP_ch2[ijack][iens].is=1; }
	else if(data_2pts.Tag[iens] == "cC211a.06.80")  {data_MP_ov_FP_ch2[ijack][iens].a = a_C.distr[ijack]; data_MP_ov_FP_ch2[ijack][iens].is=2; }
	else if(data_2pts.Tag[iens] == "cD211a.054.96")  {data_MP_ov_FP_ch2[ijack][iens].a = a_D.distr[ijack]; data_MP_ov_FP_ch2[ijack][iens].is=3; }
	else crash("Ens_tag: "+data_2pts.Tag[iens]+" not recognized");
      }
    }
  }
  //append
  bf_FF.Append_to_input_par(data_MP_ov_FP);
  bf_FF_ch2.Append_to_input_par(data_MP_ov_FP_ch2);
  //fit
  cout<<"Fitting MP/FP"<<endl;
  Bt_fit_MP_ov_FP= bf_FF.Perform_bootstrap_fit();
  Bt_fit_MP_ov_FP_ch2= bf_FF_ch2.Perform_bootstrap_fit();
  //retrieve parameters
  for(int ijack=0;ijack<Njacks;ijack++) { F0_MP_ov_FP.distr.push_back( Bt_fit_MP_ov_FP.par[ijack].F0); D1_MP_ov_FP.distr.push_back( Bt_fit_MP_ov_FP.par[ijack].D1); D2_MP_ov_FP.distr.push_back(Bt_fit_MP_ov_FP.par[ijack].D2);}
  //get ch2
  ch2_MP_ov_FP = Bt_fit_MP_ov_FP_ch2.get_ch2_ave();


  //print fitting functions for MP, FP, MP_ov_FP
  distr_t_list MP_to_print(UseJack);
  distr_t_list FP_to_print(UseJack);
  distr_t_list MP_ov_FP_to_print(UseJack);
  for(auto &a: a_to_print) {
    MP_to_print.distr_list.push_back( F0_MP + D1_MP*pow(a*fmTGeV*Lambda_QCD,2) + D2_MP*pow(a*fmTGeV*Lambda_QCD,4));
    FP_to_print.distr_list.push_back( F0_FP + D1_FP*pow(a*fmTGeV*Lambda_QCD,2) + D2_FP*pow(a*fmTGeV*Lambda_QCD,4));
    MP_ov_FP_to_print.distr_list.push_back( F0_MP_ov_FP + D1_MP_ov_FP*pow(a*fmTGeV*Lambda_QCD,2) + D2_MP_ov_FP*pow(a*fmTGeV*Lambda_QCD,4));
  }
  Print_To_File({}, {a_to_print, MP_to_print.ave(), MP_to_print.err()}, "../data/ph_emission/"+ph_type+"/"+Meson+"/mass/masses"+_Fit_tag+".fit_func", "", "#a[fm]  MP MP_err");
  Print_To_File({}, {a_to_print, FP_to_print.ave(), FP_to_print.err(), MP_ov_FP_to_print.ave(), MP_ov_FP_to_print.err()}, "../data/ph_emission/"+ph_type+"/"+Meson+"/decay_const/fP"+_Fit_tag+".fit_func", "", "#a[fm]  FP FP_err   MP/FP  MP/FP_err");


  //print ch2 summary
  cout<<"#######   SUMMARY OF Ch2 ##########"<<endl;
  cout<<"### FA ###"<<endl;
  for(int ixg=1;ixg<num_xg;ixg++) { cout<<"ch2(xg: "<<xg_t_list[ixg-1]<<"): "<<ch2_FA[ixg-1]<<endl; }
  cout<<"### FV ###"<<endl;
  for(int ixg=1;ixg<num_xg;ixg++) { cout<<"ch2(xg: "<<xg_t_list[ixg-1]<<"): "<<ch2_FV[ixg-1]<<endl; }
  cout<<"##########"<<endl;
  cout<<"ch2(MP): "<< ch2_MP<<endl;
  cout<<"ch2(FP): "<< ch2_FP<<endl;
  cout<<"ch2(MP/FP): "<<ch2_MP_ov_FP<<endl;
  cout<<"###################################"<<endl;


 
  
  //##########################################################################################
  //##########################################################################################
  
		
  


 

  }

  
 


  return rt_fit;




}
