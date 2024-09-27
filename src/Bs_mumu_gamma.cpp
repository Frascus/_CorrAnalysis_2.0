#include "../include/Bs_mumu_gamma.h"
#include "Corr_analysis.h"
#include "input.h"
#include "numerics.h"
#include "stat.h"
#include "unsupported/Eigen/src/MatrixFunctions/MatrixExponential.h"
using namespace std;


const int Nboots= 800;
const int NJ=25;
const double Qu = -1.0/3.0; //electric charge of d-type quark
const double Qd = -1.0/3.0; //electric charge of d-type quark
const double qqu = 2.0/3.0; //electric charge of u-type quark
const double qqd = -1.0/3.0; //electric charge of d-type quark
string MESON = "Bs";
const double Lambda_QCD= 0.3; //300 MeV
const bool Is_reph=true;
Vfloat Bs_xg_t_list;
Vfloat Bs_xg_to_spline;
Vfloat Bs_xg_to_spline_VMD;
Vfloat Bs_a_to_print;
const double fT_Jpsi= 0.3927; //https://arxiv.org/pdf/2008.02024.pdf  [MSbar 2GeV]
const double fT_Jpsi_err= 0.0027; //https://arxiv.org/pdf/2008.02024.pdf [MSbar 2GeV]
const int shift=0;
const string ph_type = Is_reph ? "rph" : "vph";
const int Nmasses = 5;
const Vfloat ratio_mh({1, 1.0 / 1.4835, 1.0 / 2.02827, 1.0 / 2.53531, 1.0 / 3.04236});
//const Vfloat ratio_mh({ 1.0/1.4835, 1.0/2.02827, 1.0/2.53531, 1.0/3.04236});
const double mc_MS_bar_2_ave=  1.016392574;
const double mc_MS_bar_2_err = 0.02247780935;
const bool Compute_FF = false;
const bool Skip_virtual_diagram = false;
const bool Generate_data_for_mass_spline = false;
const bool Fit_single_reg = false;
const string Reg_to_fit = "3pt";
const double k_erf = 1e10;
const bool plot_fit_func_fBs = false;
const double MU_GLB = 5.367; // GeV
const bool fit_FB_FT = false;
const bool Perform_global_fit = false;
const bool Get_KMN_JPZ_FF = false;
const bool Compute_rate=false;


void rt_FF_Bs::Print(string path) {

  string Tag="UJ_"+to_string(UseJack)+"_TF_"+to_string(Use_three_finest)+"_a4_"+to_string(Include_a4);

  vector<string> contribs_3pt({"FA", "FA_u", "FA_d", "FV",  "FV_u", "FV_d",  "FA_T", "FA_T_u", "FA_T_d", "FV_T", "FV_T_u", "FV_T_d", "FB", "FB_u", "FB_d", "FT", "FT_u", "FT_d"});
  vector<string> contribs_2pt({"mp", "fp", "fp_a4_2pt", "fp_a4_3pt", "phi", "phi_a4_2pt", "phi_a4_3pt",  "mp_ov_fp", "mp_ov_fp_a4_2pt", "mp_ov_fp_a4_3pt"});
  boost::filesystem::create_directory(path+"/jackknives");
  for(auto &c: contribs_3pt) boost::filesystem::create_directory(path+"/jackknives/"+c);
  boost::filesystem::create_directory(path+"/jackknives/2pts");

  //3pts
  for( int c=0; c<(signed)contribs_3pt.size();c++) {
    for(int ixg=1;ixg<num_xg;ixg++) {
      distr_t D= Get_FF(c).distr_list[ixg-1];
      double ch2 = Get_ch2(c)[ixg-1];
      Print_To_File({}, {D.distr}, path+"/jackknives/"+contribs_3pt[c]+"/ixg_"+to_string(ixg)+"_"+Tag+".jack", "", "ch2/dof: "+to_string_with_precision(ch2,12));
    }
  }

  //2pts
  for( int c=0; c<(signed)contribs_2pt.size();c++) {
    distr_t D= Get_FF_2pts(c);
    double ch2 = Get_ch2_2pts(c);
    Print_To_File({}, {D.distr}, path+"/jackknives/2pts/"+contribs_2pt[c]+"_"+Tag+".jack", "", "ch2/dof: "+to_string_with_precision(ch2,12));
  }
  
  return;
}

void rt_FF_Bs::Read(string path, bool UseJ, bool three_finest, bool inc_a4, int nxg) {

  this->Nmeas_K.clear();
  this->Npars_K.clear();
  this->Ndof_K.clear();

  this->UseJack= UseJ;
  this->Use_three_finest= three_finest;
  this->Include_a4 = inc_a4;
  this->num_xg= nxg;
  this->Nmeas = ( (this->Use_three_finest==true)?3:4);
  this->Npars = ( (this->Include_a4==true)?3:2);
  this->Ndof = this->Nmeas - this->Npars;
  this->Nmeas_K.push_back( ( (this->Use_three_finest==true)?6:8));
  this->Nmeas_K.push_back( ( (this->Use_three_finest==true)?6:8));
  this->Nmeas_K.push_back( ( (this->Use_three_finest==true)?6:8));
  this->Npars_K.push_back( 3);
  this->Npars_K.push_back( 4);
  this->Npars_K.push_back( 4);
  for(int k=0; k < (signed)Npars_K.size();k++) this->Ndof_K.push_back( this->Nmeas_K[k] - this->Npars_K[k]);

  

  

  
  

  string Tag="UJ_"+to_string(UseJack)+"_TF_"+to_string(Use_three_finest)+"_a4_"+to_string(Include_a4);

  vector<string> contribs_3pt({"FA", "FA_u", "FA_d", "FV",  "FV_u", "FV_d",  "FA_T", "FA_T_u", "FA_T_d", "FV_T", "FV_T_u", "FV_T_d", "FB", "FB_u", "FB_d", "FT", "FT_u", "FT_d"});
  vector<string> contribs_2pt({"mp", "fp", "fp_a4_2pt", "fp_a4_3pt", "phi", "phi_a4_2pt", "phi_a4_3pt",  "mp_ov_fp", "mp_ov_fp_a4_2pt", "mp_ov_fp_a4_3pt"});
 

  //3pts
  for(int c=0; c<(signed)contribs_3pt.size();c++) {
    distr_t_list FF(UseJack);
    Vfloat ch2;
    for(int ixg=1;ixg<num_xg;ixg++) {
      ifstream Read_ch2(path+"/jackknives/"+contribs_3pt[c]+"/ixg_"+to_string(ixg)+"_"+Tag+".jack");
      string A; double  B;  Read_ch2>>A>>B; Read_ch2.close(); ch2.push_back(B);
      FF.distr_list.emplace_back( UseJack, Read_From_File(path+"/jackknives/"+contribs_3pt[c]+"/ixg_"+to_string(ixg)+"_"+Tag+".jack", 1,  2,1));
    }
    Fill_FF(c, FF, ch2);
  }
  //2pts
  for(int c=0; c<(signed)contribs_2pt.size();c++) {
    ifstream Read_ch2(path+"/jackknives/2pts/"+contribs_2pt[c]+"_"+Tag+".jack");
    string A; double ch2;  Read_ch2>>A>>ch2; Read_ch2.close();
    distr_t FF(UseJack,  Read_From_File(path+"/jackknives/2pts/"+contribs_2pt[c]+"_"+Tag+".jack", 1,  2,1));
    Fill_FF_2pt(c,FF,ch2);
  }

 

  return;
}




pair<double,double> WC(int id, double q2, const Vfloat &G1, const Vfloat &G2, const Vfloat &BS_F, const Vfloat &Gammas_F) {

  //from https://arxiv.org/pdf/1712.07926.pdf
  if( id == 1) return  make_pair(0.147,0);  // make_pair(0.241,0); Dmitri
  if( id == 2) return make_pair(-1.053,0); //make_pair(-1.1,0); Dmitri
  if( id == 7 ) return make_pair(0.33,0); // make_pair(0.312,0);
  if( id == 10) return make_pair(4.262,0);  //  make_pair(4.41,0); Dmitri
  if( id == 9) {

    double RE_C9 =   -4.327; // -4.21;
    double IM_C9 = 0.0;
    double Cbar=  ( -0.147 + (1.053/3.0));   // 0.359/3.0;


    
    

    vector<string> Res({"J_Psi", "Psi_2S", "Psi_3770", "Psi_4040", "Psi_4160", "Psi_4230", "Psi_4415", "Psi_4660"});
    double K= M_PI/180.0;
    Vfloat Bs({0.05961, 8e-3, 9.6e-6, 10.7e-6, 6.9e-6, 3.2e-5, 2.0e-5, 1.0e-5 });
    Vfloat Bs_err({0.00033, 6e-4, 7e-7, 0.16e-5, 3.3e-6, 2.9e-5, 1e-5, 1e-5});
    //Vfloat Deltas({ 0, 0, 133*K, 301*K, 246*K});
    Vfloat Deltas(Res.size(), 0.0);
    Vfloat Gammas({ 0.0926e-3, 2.94e-4, 2.72e-2, 8e-2, 7e-2, 4.8e-2, 6.2e-2, 7.2e-2}); //GeV
    Vfloat DGammas({ 1.7e-6, 8e-6, 1e-3, 1e-2, 1e-2, 8e-3, 2.0e-2, 1.3e-2}); //GeV
    Vfloat Masses({3.0969, 3.68610,  3.7737, 4.039, 4.191, 4.2225, 4.421, 4.630}); //GeV

    for(int i=0;i<(signed)Res.size(); i++) {
      Bs[i] += BS_F[i]*Bs_err[i];
      Gammas[i] += Gammas_F[i]*DGammas[i];
    }
    
    for(int i=0; i < (signed)Res.size() ; i++) {

      double Cos = cos(G2[i]);
      double Sin = sin(G2[i]);
      
      RE_C9 += ((9.0*M_PI)/(pow(alpha_em,2)))*Cbar*(1+G1[i])*(Masses[i]*Bs[i]*Gammas[i]/( pow( q2- pow(Masses[i],2),2) + pow( Masses[i]*Gammas[i],2)))*( Cos*(q2 - pow(Masses[i],2)) + Sin*Masses[i]*Gammas[i]);
      
      IM_C9 +=  ((9.0*M_PI)/(pow(alpha_em,2)))*Cbar*(1+G1[i])*(Masses[i]*Bs[i]*Gammas[i]/( pow( q2- pow(Masses[i],2),2) + pow( Masses[i]*Gammas[i],2)))*( Sin*(q2 - pow(Masses[i],2)) -Cos*Masses[i]*Gammas[i]);
    }
    
    return make_pair(RE_C9, IM_C9);
  }

  crash("In WC, you should not be here!");
  return make_pair(0,0);
}

//####################################################



void Get_Bs_xg_t_list(int num_xg) {

  for(int ixg=1;ixg<num_xg;ixg++) Bs_xg_t_list.push_back( ixg*0.1);
  return;
}



void Get_Bs_lattice_spacings_to_print() {

  int Nlat=300;
  double sx= 0.08*1.5/(Nlat-1.0); //fm
  for(int a=0;a<Nlat;a++) Bs_a_to_print.push_back( sx*a);
}

void Get_Bs_xg_to_spline() {

  int Nxgs=300;
  double sx = 1.0/(Nxgs-1.0);
  for(int ixg=0;ixg<Nxgs;ixg++) { Bs_xg_to_spline.push_back(0.0+sx*ixg);}
    
  return;
}

void Get_Bs_xg_to_spline_VMD() {

  int Nxgs=300;
  double sx = 1.0/(Nxgs-1.0);
  for(int ixg=0;ixg<Nxgs;ixg++) { Bs_xg_to_spline_VMD.push_back(0.0+sx*ixg);}
    
  return;
}


void Get_Bs_Mh_Tmin_Tmax(string obs, string mh,  int &Tmin, int &Tmax, int ixg, string Ens) {

  if(ixg==0) { Tmin=10; Tmax=15; return;}


  bool use_old=true;

  if(mh == "B0s") {

    if( (Ens=="cA211a.12.48") ||  (Ens=="cB211b.072.64") || (Ens=="cB211b.072.96") || (Ens=="cD211a.054.96")) {

      
      
      if(obs=="FAu") {
	if(ixg==1) {  Tmin=16; Tmax=36;
	  if(Ens=="cA211a.12.48") { Tmin += 12; Tmax+=8;}
	  if(Ens=="cB211b.072.64") { Tmin +=9; Tmax+=9;}
	  if(Ens=="cD211a.054.96") { Tmin -=6; Tmax=Tmin+14;}
	}
	else if(ixg==2) {  Tmin=21; Tmax=32;
	  if(Ens=="cA211a.12.48") { Tmin -= 5; Tmax-=7;}
	  if(Ens=="cD211a.054.96") { Tmax-=5;}
	}
	else if(ixg==3) {  Tmin=21; Tmax=30;
	  if(Ens=="cA211a.12.48") { Tmin -= 7; Tmax-=7;}
	}
	else if(ixg==4) {  Tmin=22; Tmax=32;
	  if(Ens=="cA211a.12.48") { Tmin -= 10; Tmax-=8;}
	  if(Ens=="cD211a.054.96") { Tmin -=1; Tmax-=6;}
	}
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FAd") {	
	if(ixg==1) {  Tmin=18; Tmax=31;
	  if(Ens=="cA211a.12.48") { Tmin += 7; Tmax+=7;}
	}
	else if(ixg==2) {  Tmin=19; Tmax=29;	}
	else if(ixg==3) {  Tmin=19; Tmax=29;  }
	else if(ixg==4) {  Tmin=19; Tmax=28;  }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FVu") {
	if(ixg==1) {  Tmin=39; Tmax=50;  }
	else if(ixg==2) {  Tmin=39; Tmax=49;
	  if(Ens=="cA211a.12.48") { Tmin -= 8; Tmax-=8;}
	  if(Ens=="cB211b.072.64" || Ens=="cD211a.054.96") { Tmin -=10; Tmax-=10;}
	}
	else if(ixg==3) {  Tmin=38; Tmax=48;
	  if(Ens=="cA211a.12.48") { Tmin -= 10; Tmax-=10;}
	  if(Ens=="cB211b.072.64" || Ens=="cD211a.054.96") { Tmin -=10; Tmax-=10;}
	}
	else if(ixg==4) {  Tmin=38; Tmax=48;
	   if(Ens=="cA211a.12.48") { Tmin -= 12; Tmax-=12;}
	   if(Ens=="cB211b.072.64") { Tmin -=11; Tmax-=11;}
	   if(Ens=="cD211a.054.96") { Tmin -=12; Tmax-=12;}
	}
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FVd") {
	if(ixg==1) {  Tmin=35; Tmax=47;  }
	else if(ixg==2) {  Tmin=27; Tmax=38;
	  if(Ens=="cA211a.12.48") { Tmin += 7; Tmax+=7;}
	  if(Ens=="cB211b.072.64" || Ens=="cD211a.054.96") { Tmin +=9; Tmax+=9;}
	}
	else if(ixg==3) {  Tmin=29; Tmax=39;
	   if(Ens=="cB211b.072.64") { Tmin +=9; Tmax+=9;}
	   if(Ens=="cD211a.054.96") { Tmin +=10; Tmax+=10;}
	}
	else if(ixg==4) {  Tmin=28; Tmax=40;
	  if(Ens=="cB211b.072.64") { Tmin +=9; Tmax+=9;}
	   if(Ens=="cD211a.054.96") { Tmin +=10; Tmax+=10;}
	}
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FA_Tu") {
	if(ixg==1) {  Tmin=19; Tmax=31;
	  if(Ens=="cA211a.12.48") { Tmin += 3;}
	}
	else if(ixg==2) {  Tmin=21; Tmax=30;
	  if(Ens=="cB211b.072.64") { Tmin -= 4; }
	}
	else if(ixg==3) {  Tmin=20; Tmax=30;
	  if(Ens=="cB211b.072.64") { Tmin -= 4; }
	}
	else if(ixg==4) {  Tmin=16; Tmax=30;
	  if(Ens=="cB211b.072.64") { Tmin -= 5;  Tmax -=5;}
	  if(Ens=="cD211a.054.96") { Tmin += 5; Tmax+=5;}
	}
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FA_Td") {
	if(ixg==1) {  Tmin=23; Tmax=35;  }
	else if(ixg==2) {  Tmin=19; Tmax=35;
	  if(Ens=="cA211a.12.48") { Tmin += 4; Tmax+=4;}
	  if(Ens=="cD211a.054.96") { Tmin+=5; Tmax+=5;} 
	}
	else if(ixg==3) {  Tmin=18; Tmax=31;
	  if(Ens=="cD211a.054.96") { Tmin += 4; Tmax+=4;}
	}
	else if(ixg==4) {  Tmin=18; Tmax=29;
	  if(Ens=="cD211a.054.96") { Tmin += 4; Tmax+=4;}
	}
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FV_Tu") {
	if(ixg==1) {  Tmin=39; Tmax=50;  }
	else if(ixg==2) {  Tmin=41; Tmax=50;	}
	else if(ixg==3) {  Tmin=38; Tmax=46;
	  if(Ens=="cB211b.072.64") { Tmin -= 12; Tmax-=12;}
	}
	else if(ixg==4) {  Tmin=38; Tmax=44;
	  if(Ens=="cB211b.072.64") { Tmin -= 15; Tmax-=11;}
	}
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FV_Td") {  //quite bad
	if(ixg==1) {  Tmin=40; Tmax=51;  }
	else if(ixg==2) {  Tmin=40; Tmax=51;
	  if(Ens=="cA211a.12.48") { Tmin -= 8; Tmax-=8;}
	}
	else if(ixg==3) {  Tmin=38; Tmax=48;
	  if(Ens=="cA211a.12.48") { Tmin -= 9; Tmax-=9;}
	}
	else if(ixg==4) {  Tmin=37; Tmax=45;
	  if(Ens=="cA211a.12.48") { Tmin -= 10; Tmax-=5;}
	}
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
    }

    else if(Ens=="cC211a.06.80") {

      if(obs=="FAu") {
	if(ixg==1) {  Tmin=16; Tmax=36;  }
	else if(ixg==2) {  Tmin=21; Tmax=32;  }
	else if(ixg==3) {  Tmin=21; Tmax=30;  }
	else if(ixg==4) {  Tmin=22; Tmax=32;  }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FAd") {	
	if(ixg==1) {  Tmin=18; Tmax=31;  }
	else if(ixg==2) {  Tmin=14; Tmax=24;  }
	else if(ixg==3) {  Tmin=14; Tmax=27;  }
	else if(ixg==4) {  Tmin=15; Tmax=24;  }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FVu") {
	if(ixg==1) {  Tmin=38; Tmax=55; use_old=false;  }
	else if(ixg==2) {  Tmin=35; Tmax=54; use_old=false;  }
	else if(ixg==3) {  Tmin=38; Tmax=48;  }
	else if(ixg==4) {  Tmin=27; Tmax=37;  }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FVd") {
	if(ixg==1) {  Tmin=35; Tmax=47;  }
	else if(ixg==2) {  Tmin=36; Tmax=47;  }
	else if(ixg==3) {  Tmin=29; Tmax=39;  }
	else if(ixg==4) {  Tmin=28; Tmax=40;  }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FA_Tu") {
	if(ixg==1) {  Tmin=19; Tmax=40; use_old=false;  }
	else if(ixg==2) {  Tmin=25; Tmax=40; use_old=false;  }
	else if(ixg==3) {  Tmin=28; Tmax=40; use_old=false;  }
	else if(ixg==4) {  Tmin=27; Tmax=41;  }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FA_Td") {
	if(ixg==1) {  Tmin=24; Tmax=36; use_old=false; }
	else if(ixg==2) {  Tmin=28; Tmax=39; use_old=false;  }
	else if(ixg==3) {  Tmin=30; Tmax=38; use_old=false;  }
	else if(ixg==4) {  Tmin=30; Tmax=43; use_old=false;  }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FV_Tu") {
	if(ixg==1) {  Tmin=39; Tmax=50;  }
	else if(ixg==2) {  Tmin=41; Tmax=50;  }
	else if(ixg==3) {  Tmin=35; Tmax=54; use_old=false;  }
	else if(ixg==4) {  Tmin=32; Tmax=54; use_old=false;  }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FV_Td") { 
	if(ixg==1) {  Tmin=40; Tmax=57; use_old=false;  }
	else if(ixg==2) {  Tmin=35; Tmax=55; use_old=false;  }
	else if(ixg==3) {  Tmin=35; Tmax=55; use_old=false; }
	else if(ixg==4) {  Tmin=32; Tmax=55; use_old=false; }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
    }
    
    else crash("Ensemble: "+Ens+" not yet implemented");
  }


  else if(mh == "B1s") {

    if( (Ens=="cA211a.12.48") || (Ens=="cB211b.072.64") || (Ens=="cB211b.072.96") || (Ens=="cD211a.054.96") ) {

      if(obs=="FAu") {
       	if(ixg==1) {  Tmin=24; Tmax=34;
	  if(Ens=="cA211a.12.48") {Tmin+=3;}
	  if(Ens=="cB211b.072.64") {Tmin+=4;}
	}
	else if(ixg==2) {  Tmin=21; Tmax=29;
	  if(Ens=="cA211a.12.48") {Tmin-=4; Tmax -=4;}
	  if(Ens=="cD211a.054.96") {Tmin-=8; Tmax-=8;}
	}
	else if(ixg==3) {  Tmin=18; Tmax=26;	}
	else if(ixg==4) {  Tmin=15; Tmax=26;
	  if(Ens=="cD211a.054.96") { Tmin +=3;}
	}
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FAd") {
	if(ixg==1) {  Tmin=25; Tmax=30;
	   if(Ens=="cA211a.12.48") {Tmin-=3;}
	   if(Ens=="cB211b.072.64") {Tmin-=3; Tmax+=2;}
	}
	else if(ixg==2) {  Tmin=17; Tmax=25;
	  if(Ens=="cD211a.054.96") {Tmin+=6; Tmax+=8;}
	}
	else if(ixg==3) {  Tmin=16; Tmax=29;
	  if(Ens=="cD211a.054.96") {Tmin+=8; Tmax+=8;}
	}
	else if(ixg==4) {  Tmin=16; Tmax=29;
	  if(Ens=="cB211b.072.64") {Tmin-=3;}
	  if(Ens=="cD211a.054.96") {Tmin+=9; Tmax+=9;}
	}
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FVu") {
	if(ixg==1) {  Tmin=34; Tmax=45;
	  if(Ens=="cA211a.12.48") {Tmin-=3;}
	}
	else if(ixg==2) {  Tmin=34; Tmax=49;
	  if(Ens=="cA211a.12.48") {Tmin-=7;}
	}
	else if(ixg==3) {  Tmin=33; Tmax=45;
	  if(Ens=="cA211a.12.48") {Tmin-=7; Tmax-=7;}
	}
	else if(ixg==4) {  Tmin=23; Tmax=37;  }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FVd") {
	if(ixg==1) {  Tmin=33; Tmax=44;
	  if(Ens=="cA211a.12.48") { Tmin+=3;}
	  if(Ens=="cD211a.054.96") { Tmax -=6;    }

	}
	else if(ixg==2) {  Tmin=27; Tmax=36;
	  if(Ens=="cA211a.12.48") {Tmin+=3;}
	  if(Ens=="cB211b.072.64") {Tmin+=2; Tmax+=2;}
	  if(Ens=="cD211a.054.96") {Tmin+=2; Tmax+=2;}
	}
	else if(ixg==3) {  Tmin=26; Tmax=37;  }
	else if(ixg==4) {  Tmin=24; Tmax=29;
	  if(Ens=="cA211a.12.48") {Tmin-=7; Tmax+=4;}
	  if(Ens=="cB211b.072.64") {Tmin-=1; Tmax+=8;}
	  if(Ens=="cD211a.054.96") { Tmax+=4;}
	}
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FA_Tu") {
	if(ixg==1) {  Tmin=20; Tmax=28;  }
	else if(ixg==2) {  Tmin=18; Tmax=26;  }
	else if(ixg==3) {  Tmin=18; Tmax=26;  }
	else if(ixg==4) {  Tmin=18; Tmax=30;
	  if(Ens=="cD211a.054.96") {Tmin-=5; Tmax-=7;}
	}
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FA_Td") {
	if(ixg==1) {  Tmin=24; Tmax=35;  }
	else if(ixg==2) {  Tmin=24; Tmax=31;
	    if(Ens=="cA211a.12.48") {Tmin-=2; Tmax+=2;}
	}
	else if(ixg==3) {  Tmin=18; Tmax=35;
	  if(Ens=="cA211a.12.48") {Tmin+=3; Tmax+=3;}
	}
	else if(ixg==4) {  Tmin=17; Tmax=35;
	  if(Ens=="cA211a.12.48") {Tmin+=4; Tmax+=4;}
	  if(Ens=="cD211a.054.96") {Tmin+=3; Tmax+=1;}
	}
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FV_Tu") {
	if(ixg==1) {  Tmin=37; Tmax=50;
	  if(Ens=="cA211a.12.48") {Tmin-=6;}
	  if(Ens=="cD211a.054.96") {Tmin-=5; Tmax-=5;}
	  if(Ens=="cB211b.072.64") {Tmin+=3;}
	}
	else if(ixg==2) {  Tmin=33; Tmax=50;
	  if(Ens=="cA211a.12.48") {Tmin-=3; Tmax-=5;}
	  if(Ens=="cD211a.054.96") {Tmin-=3;}
	  if(Ens=="cB211b.072.64") {Tmin+=2;}
	}
	else if(ixg==3) {  Tmin=25; Tmax=37;
	  if(Ens=="cB211b.072.64") {Tmax-=2; Tmax-=5;}
	  if(Ens=="cD211a.054.96") {Tmin+=5; Tmax+=5;}
	}
	else if(ixg==4) {  Tmin=23; Tmax=33;
	  if(Ens=="cA211a.12.48") {Tmin-=3;}
	}
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FV_Td") {
	if(ixg==1) {  Tmin=34; Tmax=45;
	  if(Ens=="cD211a.054.96") {Tmin-=3; Tmax+=5;}
	}
	else if(ixg==2) {  Tmin=27; Tmax=34;
	  if(Ens=="cB211b.072.64") {Tmin+=5; Tmax+=7;}
	  if(Ens=="cD211a.054.96") {Tmin+=3; Tmax+=7;}
	  if(Ens=="cA211a.12.48")  {Tmin+=4; Tmax+=4;}
	}
	else if(ixg==3) {  Tmin=26; Tmax=33;
	  if(Ens=="cD211a.054.96") {Tmin+=3; Tmax+=8;}
	  if(Ens=="cB211b.072.64") {Tmin+=6; Tmax+=6;}
	}
	else if(ixg==4) {  Tmin=23; Tmax=31;
	  if(Ens=="cA211a.12.48") {Tmin-=7;}
	  if(Ens=="cB211b.072.64") {Tmin+=5; Tmax+=5;}
	  if(Ens=="cD211a.054.96") {Tmin+=9; Tmax+=9;}
	}
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
    }

    else if(Ens=="cC211a.06.80") {
    
      if(obs=="FAu") {
       	if(ixg==1) {  Tmin=24; Tmax=34;  }
	else if(ixg==2) {  Tmin=17; Tmax=29;  }
	else if(ixg==3) {  Tmin=16; Tmax=26;  }
	else if(ixg==4) {  Tmin=15; Tmax=26;  }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FAd") {
	if(ixg==1) {  Tmin=22; Tmax=37; use_old=false;  }
	else if(ixg==2) {  Tmin=18; Tmax=25;  }
	else if(ixg==3) {  Tmin=13; Tmax=29;  }
	else if(ixg==4) {  Tmin=12; Tmax=29;  }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FVu") {
	if(ixg==1) {  Tmin=38; Tmax=51; use_old=false;  }
	else if(ixg==2) {  Tmin=34; Tmax=49;  }
	else if(ixg==3) {  Tmin=33; Tmax=45;  }
	else if(ixg==4) {  Tmin=24; Tmax=33; use_old=false;  }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FVd") {
	if(ixg==1) {  Tmin=32; Tmax=43;  }
	else if(ixg==2) {  Tmin=27; Tmax=36;  }
	else if(ixg==3) {  Tmin=24; Tmax=31;  }
	else if(ixg==4) {  Tmin=22; Tmax=28;  }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FA_Tu") {
	if(ixg==1) {  Tmin=20; Tmax=28;  }
	else if(ixg==2) {  Tmin=18; Tmax=26;  }
	else if(ixg==3) {  Tmin=18; Tmax=26;  }
	else if(ixg==4) {  Tmin=18; Tmax=30;  }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FA_Td") {
	if(ixg==1) {  Tmin=28; Tmax=47; use_old=false;  }
	else if(ixg==2) {  Tmin=26; Tmax=37; use_old=false;  }
	else if(ixg==3) {  Tmin=20; Tmax=32; use_old=false;  }
	else if(ixg==4) {  Tmin=19; Tmax=29; use_old=false;  }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FV_Tu") {
	if(ixg==1) {  Tmin=37; Tmax=50;  }
	else if(ixg==2) {  Tmin=31; Tmax=42;  }
	else if(ixg==3) {  Tmin=30; Tmax=47; use_old=false;  }
	else if(ixg==4) {  Tmin=24; Tmax=39; use_old=false;  }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FV_Td") {
	if(ixg==1) {  Tmin=35; Tmax=45; }
	else if(ixg==2) {  Tmin=32; Tmax=46; use_old=false;}
	else if(ixg==3) {  Tmin=22; Tmax=33;  }
	else if(ixg==4) {  Tmin=17; Tmax=31;  }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
    }
    else crash("Ensemble: "+Ens+" not yet implemented");
  }


  else if(mh == "B2s") {

    if( (Ens=="cA211a.12.48") ||  (Ens=="cB211b.072.64") || (Ens=="cB211b.072.96") || (Ens=="cD211a.054.96")  ) {

      if(obs=="FAu") {
	if(ixg==1) {  Tmin=26; Tmax=34;
	  if(Ens=="cA211a.12.48") {Tmin-=5; Tmax-=5;}
	  if(Ens=="cB211b.072.64") {Tmin-=4;}
	}
	else if(ixg==2) {  Tmin=27; Tmax=34;
	  if(Ens=="cA211a.12.48") {Tmin-=7; Tmax-=7;}
	  if(Ens=="cD211a.054.96") {Tmin-=3;}
	}
	else if(ixg==3) {  Tmin=25; Tmax=34;
	  if(Ens=="cA211a.12.48") {Tmin-=7; Tmax-=7;}
	  if(Ens=="cD211a.054.96") {Tmin-=14; Tmax-=14;}
	}
	else if(ixg==4) {  Tmin=26; Tmax=34;
	  if(Ens=="cA211a.12.48") {Tmin-=7; Tmax-=7;}
	  if(Ens=="cB211b.072.64") {Tmin-=6;}
	  if(Ens=="cD211a.054.96") {Tmin-=11; Tmax-=11;}
	}
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FAd") {
	if(ixg==1) {  Tmin=27; Tmax=34;
	  if(Ens=="cA211a.12.48") {Tmin-=6; Tmax-=6;}
	  if(Ens=="cB211b.072.64") {Tmin-=8;}
	  if(Ens=="cD211a.054.96") {Tmin-=12;}
	}
	else if(ixg==2) {  Tmin=24; Tmax=35;
	  if(Ens=="cA211a.12.48") {Tmin-=8; Tmax-=7;}
	  if(Ens=="cB211b.072.64") {Tmin-=8; Tmax-=6;}
	  if(Ens=="cD211a.054.96") {Tmin-=2; Tmax -=3;}
	}
	else if(ixg==3) {  Tmin=24; Tmax=35;
	  if(Ens=="cA211a.12.48") {Tmin-=12; Tmax-=13;}
	  if(Ens=="cB211b.072.64") {Tmin-=12; Tmax-=6;}
	  if(Ens=="cD211a.054.96") {Tmin-=1; Tmax-=3;}
	}
	else if(ixg==4) {  Tmin=24; Tmax=35;
	   if(Ens=="cA211a.12.48") {Tmin-=14; Tmax-=15;}
	   if(Ens=="cB211b.072.64") {Tmin-=14; Tmax-=16;}
	   if(Ens=="cD211a.054.96") {Tmin-=4;}
	}
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FVu") {
	if(ixg==1) {  Tmin=31; Tmax=39;
	  if(Ens=="cA211a.12.48") {Tmin-=3; Tmax-=4;}
	  if(Ens=="cB211b.072.64") {Tmin-=1; Tmax+=2;}
	}
	else if(ixg==2) {  Tmin=33; Tmax=43;
	  if(Ens=="cA211a.12.48") {Tmin-=6; Tmax-=8;}
	}
	else if(ixg==3) {  Tmin=22; Tmax=32;
	  if(Ens=="cB211b.072.64") {Tmin+=5; Tmax+=5;}
	  if(Ens=="cD211a.054.96") {Tmin+=5; Tmax+=5;}
	}
	else if(ixg==4) {  Tmin=21; Tmax=32;
	  if(Ens=="cA211a.12.48") {Tmin-=3;}
	  if(Ens=="cB211b.072.64") {Tmin-=3;}
	}
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FVd") {
	if(ixg==1) {  Tmin=31; Tmax=45;
	  if(Ens=="cA211a.12.48") {Tmin-=1;}
	  if(Ens=="cD211a.054.96") { Tmax-=7;}
	  if(Ens=="cB211b.072.64") {Tmin+=2;}
	}
	else if(ixg==2) {  Tmin=27; Tmax=36;
	  if(Ens=="cB211b.072.64") {Tmin+=3; Tmax+=3;}
	  if(Ens=="cD211a.054.96") {Tmin-=1; Tmax+=0;}
	}
	else if(ixg==3) {  Tmin=23; Tmax=30;
	  if(Ens=="cA211a.12.48") {Tmin-=4;}
	  if(Ens=="cB211b.072.64") {Tmin+=3; Tmax+=7;}
	}
	else if(ixg==4) {  Tmin=21; Tmax=30;
	  if(Ens=="cA211a.12.48") {Tmin-=8; Tmax-=8;}
	  if(Ens=="cB211b.072.64") {Tmin-=4; Tmax-=8;}
	}
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FA_Tu") {
	if(ixg==1) {  Tmin=18; Tmax=25;  }
	else if(ixg==2) {  Tmin=18; Tmax=24;
	  if(Ens=="cA211a.12.48") {Tmin-=2; Tmax-=2;}
	}
	else if(ixg==3) {  Tmin=20; Tmax=28;
	  if(Ens=="cA211a.12.48") {Tmin-=2; Tmax-=2;}
	}
	else if(ixg==4) {  Tmin=19; Tmax=28;  }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FA_Td") {
	if(ixg==1) {  Tmin=23; Tmax=30;
	  if(Ens=="cA211a.12.48") {Tmin-=2;}
	}
	else if(ixg==2) {  Tmin=20; Tmax=30;  }
	else if(ixg==3) {  Tmin=20; Tmax=30;  }
	else if(ixg==4) {  Tmin=19; Tmax=30;
	  if(Ens=="cA211a.12.48") {Tmin-=1; Tmax-=3;}
	  if(Ens=="cB211b.072.64") { Tmin -= 4; Tmax-=4;   }
	}
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FV_Tu") {
	if(ixg==1) {  Tmin=31; Tmax=41;
	  if(Ens=="cA211a.12.48") {Tmin-=1; Tmax-=1;}
	}
	else if(ixg==2) {  Tmin=31; Tmax=39;
	  if(Ens=="cA211a.12.48") {Tmin-=7; Tmax-=7;}
	  if(Ens=="cB211b.072.64") {Tmin-=3;}
	}
	else if(ixg==3) {  Tmin=24; Tmax=35;
	  if(Ens=="cD211a.054.96") { Tmin += 4; Tmax+=4;   }
	  if(Ens=="cA211a.12.48") {Tmin+=2;}
	}
	else if(ixg==4) {  Tmin=24; Tmax=38;
	  if(Ens=="cA211a.12.48") {Tmin-=4;}
	  if(Ens=="cB211b.072.64") { Tmin -= 3; Tmax-=3;   }
	  if(Ens=="cD211a.054.96") { Tmin += 3; Tmax+=3;   }
	}
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FV_Td") {
	if(ixg==1) {  Tmin=29; Tmax=35;
	  if(Ens=="cA211a.12.48") {Tmin-=0; Tmax-=1;}
	  if(Ens=="cB211b.072.64") { Tmin += 3; Tmax+=3;   }
	  if(Ens=="cD211a.054.96") { Tmin += 3; Tmax+=3;   }
	}
	else if(ixg==2) {  Tmin=27; Tmax=34;
	  if(Ens=="cA211a.12.48") {Tmin-=2; Tmax-=2;}
	}
	else if(ixg==3) {  Tmin=24; Tmax=34;
	  if(Ens=="cA211a.12.48") {Tmin-=5;}
	}
	else if(ixg==4) {  Tmin=20; Tmax=30;
	  if(Ens=="cA211a.12.48") {Tmin-=5;}
	  if(Ens=="cB211b.072.64") { Tmin -=6; Tmax-=6;   }
	}
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
    }

    else if(Ens=="cC211a.06.80") {

      if(obs=="FAu") {
	if(ixg==1) {  Tmin=26; Tmax=34;  }
	else if(ixg==2) {  Tmin=23; Tmax=30;  }
	else if(ixg==3) {  Tmin=15; Tmax=24;  }
	else if(ixg==4) {  Tmin=17; Tmax=25;  }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FAd") {
	if(ixg==1) {  Tmin=18; Tmax=32;  }
	else if(ixg==2) {  Tmin=12; Tmax=33;  }
	else if(ixg==3) {  Tmin=14; Tmax=33;  }
	else if(ixg==4) {  Tmin=13; Tmax=22;  }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FVu") {
	if(ixg==1) {  Tmin=31; Tmax=43;  }
	else if(ixg==2) {  Tmin=33; Tmax=43;  }
	else if(ixg==3) {  Tmin=28; Tmax=37;  }
	else if(ixg==4) {  Tmin=19; Tmax=35;  }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FVd") {
	if(ixg==1) {  Tmin=30; Tmax=39;  }
	else if(ixg==2) {  Tmin=24; Tmax=36;  }
	else if(ixg==3) {  Tmin=24; Tmax=30;  }
	else if(ixg==4) {  Tmin=16; Tmax=30;  }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FA_Tu") {
	if(ixg==1) {  Tmin=16; Tmax=22;  }
	else if(ixg==2) {  Tmin=18; Tmax=23;  }
	else if(ixg==3) {  Tmin=20; Tmax=28;  }
	else if(ixg==4) {  Tmin=19; Tmax=28;  }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FA_Td") {
	if(ixg==1) {  Tmin=29; Tmax=34;  }
	else if(ixg==2) {  Tmin=20; Tmax=30;  }
	else if(ixg==3) {  Tmin=17; Tmax=27;  }
	else if(ixg==4) {  Tmin=12; Tmax=23;  }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FV_Tu") {
	if(ixg==1) {  Tmin=31; Tmax=41;  }
	else if(ixg==2) {  Tmin=29; Tmax=39;  }
	else if(ixg==3) {  Tmin=24; Tmax=38;  }
	else if(ixg==4) {  Tmin=24; Tmax=38;  }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FV_Td") {
	if(ixg==1) {  Tmin=31; Tmax=33;  }
	else if(ixg==2) {  Tmin=25; Tmax=34;  }
	else if(ixg==3) {  Tmin=24; Tmax=34;  }
	else if(ixg==4) {  Tmin=20; Tmax=30;  }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
    }
    
    else crash("Ensemble: "+Ens+" not yet implemented");
  }



  else if(mh == "B3s") {

    if( (Ens=="cA211a.12.48") || (Ens=="cB211b.072.64") || (Ens=="cB211b.072.96") || (Ens=="cD211a.054.96")  ) {

      if(obs=="FAu") {
       if(ixg==1) {  Tmin=21; Tmax=27;
	 if(Ens=="cA211a.12.48") { Tmin-=3; Tmax-=3;}
	 if(Ens=="cD211a.054.96") { Tmin -= 1; Tmax -=4 ;}
       }
       else if(ixg==2) {  Tmin=20; Tmax=25;
	 if(Ens=="cA211a.12.48") { Tmin-=3; Tmax-=3;}
       }
       else if(ixg==3) {  Tmin=21; Tmax=27;
	 if(Ens=="cA211a.12.48") { Tmin-=4; Tmax-=4;}
       }
       else if(ixg==4) {  Tmin=21; Tmax=27;
	 if(Ens=="cA211a.12.48") { Tmin-=3; Tmax-=3;}
	 if(Ens=="cD211a.054.96") { Tmin -= 3; Tmax -=3 ;}
	 
       }
       else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FAd") { //non bellissimi
	if(ixg==1) {  Tmin=17; Tmax=24;
	  if(Ens=="cA211a.12.48") { Tmin-=1; Tmax+=2;}
	}
	else if(ixg==2) {  Tmin=24; Tmax=31;
	  if(Ens=="cA211a.12.48") { Tmin-=12; Tmax-=12;}
	  if(Ens=="cB211b.072.64") { Tmin-=12; Tmax-=13;}
	  if(Ens=="cD211a.054.96") { Tmin-=3;}
	}
	else if(ixg==3) {  Tmin=22; Tmax=33;
	  if(Ens=="cA211a.12.48") { Tmin-=10; Tmax-=15;}
	  if(Ens=="cB211b.072.64") { Tmin-=13; Tmax-=13;}
	  if(Ens=="cD211a.054.96") { Tmin-=0;}
	}
	else if(ixg==4) {  Tmin=22; Tmax=30;
	  if(Ens=="cA211a.12.48") { Tmin-=14; Tmax-=11;}
	  if(Ens=="cD211a.054.96") { Tmin -= 7; Tmax -=4 ;}
	  if(Ens=="cB211b.072.64") { Tmin -= 13; Tmax -=15;}
	}
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FVu") {
	if(ixg==1) {  Tmin=28; Tmax=40;
	  if(Ens=="cA211a.12.48") { Tmin-=3; Tmax-=3;}
	}
	else if(ixg==2) {  Tmin=27; Tmax=32;
	  if(Ens=="cA211a.12.48") { Tmin-=4; Tmax-=4;}
	  if(Ens=="cD211a.054.96") { Tmin -= 1; Tmax +=4 ;}
	}
	else if(ixg==3) {  Tmin=22; Tmax=30;
	  if(Ens=="cA211a.12.48") { Tmin+=2; Tmax+=2;}
	  if(Ens=="cB211b.072.64") { Tmin+=3; Tmax+=3;}
	}
	else if(ixg==4) {  Tmin=16; Tmax=30;  }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FVd") {
	if(ixg==1) {  Tmin=26; Tmax=32;
	  if(Ens=="cA211a.12.48") { Tmin-=1; Tmax-=4;}
	  if(Ens=="cD211a.054.96") { Tmin+=3; Tmax+=3;}
	  if(Ens=="cB211b.072.64") {Tmin +=2; Tmax+=2;}
	}
	else if(ixg==2) {  Tmin=23; Tmax=30;
	  if(Ens=="cA211a.12.48") { Tmin-=5; Tmax-=5;}
	  if(Ens=="cD211a.054.96") {Tmin+=1; Tmax+=1;}
	  if(Ens=="cB211b.072.64") {Tmin+=2;}
	}
	else if(ixg==3) {  Tmin=20; Tmax=25;
	  if(Ens=="cA211a.12.48") { Tmin-=4;}
	  if(Ens=="cD211a.054.96") { Tmin += 3; Tmax +=5 ;}
	}
	else if(ixg==4) {  Tmin=17; Tmax=24;
	  if(Ens=="cA211a.12.48") { Tmin-=5; Tmax-=5;}
	  if(Ens=="cB211b.072.64") { Tmin-=3; Tmax-=5;}
	  if(Ens=="cD211a.054.96") { Tmin += 2; Tmax +=2 ;}
	}
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FA_Tu") {
	if(ixg==1) {  Tmin=17; Tmax=26;
	  if(Ens=="cA211a.12.48") { Tmin-=4; Tmax-=4;}
	}
	else if(ixg==2) {  Tmin=13; Tmax=24;
	  if(Ens=="cA211a.12.48") { Tmin+=3; Tmax-=2;}
	  if(Ens=="cB211b.072.64") { Tmin+=5; Tmax+=5;}
	}
	else if(ixg==3) {  Tmin=13; Tmax=24;
	   if(Ens=="cA211a.12.48") { Tmin+=3; }
	   if(Ens=="cB211b.072.64") { Tmin+=3;}
	   if(Ens=="cD211a.054.96") { Tmin += 3;}
	}
	else if(ixg==4) {  Tmin=16; Tmax=28;  }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FA_Td") {
	if(ixg==1) {  Tmin=20; Tmax=27;
	   if(Ens=="cA211a.12.48") { Tmin-=2; Tmax-=2;}
	}
	else if(ixg==2) {  Tmin=17; Tmax=27;
	  if(Ens=="cA211a.12.48") { Tmin+=4; Tmax+=2;}
	  if(Ens=="cB211b.072.64") { Tmin-=1; Tmax-=4;}
	}
	else if(ixg==3) {  Tmin=16; Tmax=27;
	  if(Ens=="cB211b.072.64") { Tmin-=3; Tmax-=5;}
	}
	else if(ixg==4) {  Tmin=14; Tmax=22;
	  if(Ens=="cB211b.072.64") { Tmin-=4; Tmax-=7;}
	}
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FV_Tu") { //two plataux for small xg
	if(ixg==1) {  Tmin=31; Tmax=38;
	  if(Ens=="cA211a.12.48") { Tmin-=6; Tmax-=6;}
	  if(Ens=="cB211b.072.64") { Tmin-=6; Tmax-=6;}
	}
	else if(ixg==2) {  Tmin=25; Tmax=30;
	  if(Ens=="cD211a.054.96") { Tmin += 5; Tmax +=5 ;}
	}
	else if(ixg==3) {  Tmin=20; Tmax=28;
	  if(Ens=="cB211b.072.64") { Tmin+=2; Tmax+=2;}
	  if(Ens=="cD211a.054.96") { Tmin += 8; Tmax +=8 ;}
	}
	else if(ixg==4) {  Tmin=14; Tmax=25;
	  if(Ens=="cA211a.12.48") { Tmin+=2;}
	  if(Ens=="cB211b.072.64") { Tmin+=6; Tmax+=6;}
	}
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FV_Td") {
	if(ixg==1) {  Tmin=26; Tmax=40;
	  if(Ens=="cA211a.12.48") { Tmin-=1; Tmax-=1;}
	  if(Ens=="cB211b.072.64") { Tmin+=2; Tmax+=2;}
	  if(Ens=="cD211a.054.96") {Tmin +=3;}
	}
	else if(ixg==2) {  Tmin=25; Tmax=34;
	  if(Ens=="cA211a.12.48") { Tmin-=7; Tmax-=9;}
	  if(Ens=="cD211a.054.96") {Tmin +=3; Tmax+=3;}
	}
	else if(ixg==3) {  Tmin=20; Tmax=25;
	  if(Ens=="cA211a.12.48") { Tmin-=4;}
	  if(Ens=="cB211b.072.64") { Tmin-=2; Tmax-=2;}
	  if(Ens=="cD211a.054.96") {Tmin +=6; Tmax+=6;}
	}
	else if(ixg==4) {  Tmin=17; Tmax=23;
	  if(Ens=="cA211a.12.48") { Tmin-=4; Tmax-=4;}
	  if(Ens=="cB211b.072.64") { Tmin-=3;}
	  if(Ens=="cD211a.054.96") {Tmin +=2; Tmax+=2;}
	}
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
    }

    else if(Ens=="cC211a.06.80") {

      if(obs=="FAu") {
       if(ixg==1) {  Tmin=21; Tmax=27;  }
       else if(ixg==2) {  Tmin=20; Tmax=25;  }
       else if(ixg==3) {  Tmin=21; Tmax=27;  }
       else if(ixg==4) {  Tmin=21; Tmax=27;  }
       else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FAd") { //non bellissimi
	if(ixg==1) {  Tmin=17; Tmax=24;  }
	else if(ixg==2) {  Tmin=12; Tmax=19;  }
	else if(ixg==3) {  Tmin=12; Tmax=22;  }
	else if(ixg==4) {  Tmin=12; Tmax=26;  }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FVu") {
	if(ixg==1) {  Tmin=30; Tmax=40;  }
	else if(ixg==2) {  Tmin=26; Tmax=33;  }
	else if(ixg==3) {  Tmin=22; Tmax=30;  }
	else if(ixg==4) {  Tmin=16; Tmax=30;  }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FVd") {
	if(ixg==1) {  Tmin=26; Tmax=32;  }
	else if(ixg==2) {  Tmin=24; Tmax=31;  }
	else if(ixg==3) {  Tmin=20; Tmax=25;  }
	else if(ixg==4) {  Tmin=16; Tmax=25;  }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FA_Tu") {
	if(ixg==1) {  Tmin=17; Tmax=26;  }
	else if(ixg==2) {  Tmin=19; Tmax=24;  }
	else if(ixg==3) {  Tmin=17; Tmax=24;  }
	else if(ixg==4) {  Tmin=20; Tmax=28;  }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FA_Td") {
	if(ixg==1) {  Tmin=20; Tmax=27;  }
	else if(ixg==2) {  Tmin=19; Tmax=31;  }
	else if(ixg==3) {  Tmin=12; Tmax=27;  }
	else if(ixg==4) {  Tmin=14; Tmax=22;  }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FV_Tu") { //two plataux for small xg
	if(ixg==1) {  Tmin=31; Tmax=38;  }
	else if(ixg==2) {  Tmin=26; Tmax=30;  }
	else if(ixg==3) {  Tmin=25; Tmax=33;  }
	else if(ixg==4) {  Tmin=23; Tmax=31;  }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FV_Td") {
	if(ixg==1) {  Tmin=26; Tmax=36;  }
	else if(ixg==2) {  Tmin=21; Tmax=30;  }
	else if(ixg==3) {  Tmin=20; Tmax=27;  }
	else if(ixg==4) {  Tmin=15; Tmax=23;  }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
    }
    
    else crash("Ensemble: "+Ens+" not yet implemented");
  }

  
  else if(mh == "B4s") {

    if( (Ens=="cA211a.12.48") ||  (Ens=="cB211b.072.64") || (Ens=="cB211b.072.96") || (Ens=="cD211a.054.96")  ) {

      if(obs=="FAu") { //not very nice
       if(ixg==1) {  Tmin=18; Tmax=25;
	 if(Ens=="cA211a.12.48") { Tmin-=1;}
	 if(Ens=="cB211b.072.64") {Tmin+=1;}
       }
       else if(ixg==2) {  Tmin=20; Tmax=25;
	 if(Ens=="cA211a.12.48") { Tmin-=3;}
	 if(Ens=="cB211b.072.64") {Tmin-=1; Tmax+=1;}
       }
       else if(ixg==3) {  Tmin=21; Tmax=25;
	 if(Ens=="cA211a.12.48") { Tmin-=4;}
	 if(Ens=="cB211b.072.64") {Tmin-=1; Tmax+=1;}
       }
       else if(ixg==4) {  Tmin=21; Tmax=25;
	 if(Ens=="cA211a.12.48") { Tmin-=5;}
	 if(Ens=="cB211b.072.64") {Tmin-=2; Tmax+=2;}
       }
       else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FAd") { //not very nice
	if(ixg==1) {  Tmin=21; Tmax=24;
	  if(Ens=="cA211a.12.48") { Tmin-=6; Tmax-=5;}
	  if(Ens=="cB211b.072.64") {Tmin-=6; Tmax+=2;}
	  if(Ens=="cD211a.054.96") { Tmin -= 5; Tmax+=0;}
	}
	else if(ixg==2) {  Tmin=21; Tmax=24;
	  if(Ens=="cA211a.12.48") { Tmin-=7; Tmax-=5;}
	  if(Ens=="cB211b.072.64") {Tmin-=10; Tmax-=10;}
	  if(Ens=="cD211a.054.96") { Tmin += 2; Tmax+=6;}
	}
	else if(ixg==3) {  Tmin=21; Tmax=25;
	  if(Ens=="cA211a.12.48") { Tmin-=7; Tmax-=7;}
	  if(Ens=="cB211b.072.64") {Tmin-=13; Tmax-=13;}
	  if(Ens=="cD211a.054.96") { Tmin -= 3; Tmax+=4;}
	}
	else if(ixg==4) {  Tmin=21; Tmax=26;
	  if(Ens=="cA211a.12.48") { Tmin-=7; Tmax-=5;}
	  if(Ens=="cB211b.072.64") {Tmin-=12; Tmax-=12;}
	  if(Ens=="cD211a.054.96") {Tmin -=11; Tmax-=9;}
	}
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FVu") {
	if(ixg==1) {  Tmin=28; Tmax=33;
	  if(Ens=="cA211a.12.48") { Tmin-=6; Tmax-=5;}
	}
	else if(ixg==2) {  Tmin=22; Tmax=28;
	  if(Ens=="cA211a.12.48") { Tmin-=2; Tmax-=2;}
	  if(Ens=="cB211b.072.64") {Tmin+=3; Tmax+=3;}
	  if(Ens=="cD211a.054.96") { Tmin += 5; Tmax+=5;}
	}
	else if(ixg==3) {  Tmin=19; Tmax=24;
	  if(Ens=="cB211b.072.64") {Tmin+=2;}
	}
	else if(ixg==4) {  Tmin=16; Tmax=24;
	  if(Ens=="cB211b.072.64") {Tmin+=2;}
	}
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FVd") {
	if(ixg==1) {  Tmin=22; Tmax=30;
	  if(Ens=="cA211a.12.48") { Tmin-=3; Tmax-=5;}
	  if(Ens=="cD211a.054.96") { Tmin += 8; Tmax+=4;}
	  if(Ens=="cB211b.072.64") { Tmin +=2; }
	}
	else if(ixg==2) {  Tmin=21; Tmax=27;
	  if(Ens=="cB211b.072.64") { Tmin-=1; Tmax-=1;}
	  if(Ens=="cA211a.12.48") { Tmin-=3; Tmax-=3;}
	  if(Ens=="cD211a.054.96") { Tmin += 3;}
	}
	else if(ixg==3) {  Tmin=17; Tmax=24;
	  if(Ens=="cB211b.072.64") {Tmin+=1; Tmax-=1;}
	  if(Ens=="cD211a.054.96") {Tmin+=4; Tmax += 4; }
	}
	else if(ixg==4) {  Tmin=14; Tmax=23;
	  if(Ens=="cA211a.12.48") { Tmin-=2; Tmax-=6;}
	}
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FA_Tu") {
	if(ixg==1) {  Tmin=15; Tmax=24; if(Ens=="cB211b.072.64") Tmin -=4;  }
	else if(ixg==2) {  Tmin=17; Tmax=24;  }
	else if(ixg==3) {  Tmin=17; Tmax=24;  }
	else if(ixg==4) {  Tmin=18; Tmax=28;  }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FA_Td") {
	if(ixg==1) {  Tmin=18; Tmax=27;
	  if(Ens=="cA211a.12.48") { Tmin-=1; Tmax-=5;}
	  if(Ens=="cD211a.054.96") {Tmin+=1; Tmax += 4; }
	}
	else if(ixg==2) {  Tmin=18; Tmax=28;
	  if(Ens=="cA211a.12.48") { Tmin-=2; Tmax=Tmin+4;}
	  if(Ens=="cD211a.054.96") {Tmin-=2; Tmax += 3; }
	  if(Ens=="cB211b.072.64") {Tmax-=4;}
	}
	else if(ixg==3) {  Tmin=15; Tmax=23;
	  if(Ens=="cA211a.12.48") { Tmin-=5; Tmax-=7;}
	  if(Ens=="cB211b.072.64") { Tmin-=3; }
	  if(Ens=="cD211a.054.96") { Tmin -=6; Tmax -=9;}
	}
	else if(ixg==4) {  Tmin=12; Tmax=21;
	  if(Ens=="cD211a.054.96") {Tmin-=1; Tmax -= 4; }
	}
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FV_Tu") { //bad plateaux at  xg=0.1
	if(ixg==1) {  Tmin=26; Tmax=30;
	  if(Ens=="cA211a.12.48") { Tmin-=5; Tmax-=5;}
	  if(Ens=="cB211b.072.64") {Tmin-=1; Tmax += 2; }
	  if(Ens=="cD211a.054.96") {Tmin+=2; Tmax += 2; }
	}
	else if(ixg==2) {  Tmin=21; Tmax=25;
	  if(Ens=="cA211a.12.48") { Tmin-=2;}
	  if(Ens=="cB211b.072.64") {Tmin+=4; Tmax += 4; }
	  if(Ens=="cD211a.054.96") {Tmin+=8; Tmax += 8; }
	}
	else if(ixg==3) {  Tmin=22; Tmax=25;
	  if(Ens=="cA211a.12.48") { Tmin-=4; Tmax-=2;}
	  if(Ens=="cD211a.054.96") {Tmin+=4; Tmax += 4; }
	}
	else if(ixg==4) {  Tmin=20; Tmax=26;
	  if(Ens=="cA211a.12.48") { Tmin-=3; Tmax-=3;}
	}
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FV_Td") {
	if(ixg==1) {  Tmin=26; Tmax=30;
	  if(Ens=="cA211a.12.48") { Tmin-=6; Tmax-=7;}
	  if(Ens=="cD211a.054.96") {Tmin+=3; Tmax += 6; }
	  if(Ens=="cB211b.072.64") {Tmin-=2; Tmax -=1; }
	}
	else if(ixg==2) {  Tmin=24; Tmax=30;
	    if(Ens=="cA211a.12.48") { Tmin-=6; Tmax-=6;}
	    if(Ens=="cD211a.054.96") {Tmin+=3; Tmax += 4; }
	    if(Ens=="cB211b.072.64") {Tmin-=4; Tmax-=4;}
	}
	else if(ixg==3) {  Tmin=16; Tmax=24;
	  if(Ens=="cA211a.12.48") { Tmin-=3; Tmax-=5;}
	  if(Ens=="cD211a.054.96") {Tmin+=6; Tmax += 2; }
	}
	else if(ixg==4) {  Tmin=11; Tmax=17;
	  if(Ens=="cA211a.12.48") { Tmin+=1; Tmax+=1;}
	}
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
    }
    else if(Ens=="cC211a.06.80") {

      if(obs=="FAu") { //not very nice
       if(ixg==1) {  Tmin=18; Tmax=25;  }
       else if(ixg==2) {  Tmin=18; Tmax=25;  }
       else if(ixg==3) {  Tmin=18; Tmax=25;  }
       else if(ixg==4) {  Tmin=21; Tmax=25;  }
       else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FAd") { //not very nice
	if(ixg==1) {  Tmin=15; Tmax=21;  }
	else if(ixg==2) {  Tmin=12; Tmax=17;  }
	else if(ixg==3) {  Tmin=13; Tmax=17;  }
	else if(ixg==4) {  Tmin=11; Tmax=16;  }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FVu") {
	if(ixg==1) {  Tmin=25; Tmax=33;  }
	else if(ixg==2) {  Tmin=28; Tmax=33;  }
	else if(ixg==3) {  Tmin=21; Tmax=24;  }
	else if(ixg==4) {  Tmin=16; Tmax=24;  }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FVd") {
	if(ixg==1) {  Tmin=24; Tmax=30;  }
	else if(ixg==2) {  Tmin=23; Tmax=26;  }
	else if(ixg==3) {  Tmin=18; Tmax=27;  }
	else if(ixg==4) {  Tmin=15; Tmax=23;  }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FA_Tu") {
	if(ixg==1) {  Tmin=15; Tmax=24;  }
	else if(ixg==2) {  Tmin=17; Tmax=24;  }
	else if(ixg==3) {  Tmin=17; Tmax=24;  }
	else if(ixg==4) {  Tmin=18; Tmax=28;  }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FA_Td") {
	if(ixg==1) {  Tmin=18; Tmax=27;  }
	else if(ixg==2) {  Tmin=17; Tmax=29;  }
	else if(ixg==3) {  Tmin=10; Tmax=16;  }
	else if(ixg==4) {  Tmin=12; Tmax=21;  }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FV_Tu") { //bad plateaux at  xg=0.1
	if(ixg==1) {  Tmin=28; Tmax=32;  }
	else if(ixg==2) {  Tmin=29; Tmax=33;  }
	else if(ixg==3) {  Tmin=26; Tmax=29;  }
	else if(ixg==4) {  Tmin=20; Tmax=26;  }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FV_Td") {
	if(ixg==1) {  Tmin=24; Tmax=30;  }
	else if(ixg==2) {  Tmin=20; Tmax=26;  }
	else if(ixg==3) {  Tmin=16; Tmax=24;  }
	else if(ixg==4) {  Tmin=11; Tmax=17;  }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
    }
    else crash("Ensemble: "+Ens+" not yet implemented");
  }

  else crash("mh: "+mh+" not yet implemented");
  

  
  if(use_old && Ens=="cC211a.06.80") { Tmin= (int)(Tmin*0.07957/0.0682); Tmax= (int)( Tmax*0.07957/0.0682);}

  if(use_old && Ens=="cD211a.054.96") { Tmin= (int)(Tmin*0.07957/0.0569); Tmax= (int)( Tmax*0.07957/0.0569);}

  if(use_old && Ens=="cA211a.12.48") { Tmin=(int)(Tmin*0.07957/0.0907593); Tmax= (int)(Tmax*0.07957/0.0907593);}
 

  return;

}











void Compute_Bs_mumu_gamma() {

  
  //test WC
  //testing WC
  for(int ixg=0;ixg<100;ixg++) {
    double xg= 0.42 -ixg*(0.42/99);
    double q2= MBs*MBs*( 1 - xg);
    cout<<"xg: "<<xg<<" q2: "<<q2<<" C9: -4.21: C9_eff[RE]: "<<WC(9, q2).first<<" C9_eff[IM]: "<<WC(9, q2).second<<endl;
  }


  int Npt_test=999;

  auto F_null = [](double x) { return 0.0;};

  for(int i=0;i<Npt_test;i++) {
    double x= (i+1)*0.2/1000;
    cout<<x<<" "<<Compute_Bs_mumugamma_differential_decay_rate(x, F_null, F_null, F_null, F_null, F_null, F_null, 0.2245, x_b, Vtbs, tBs, Vfloat(8,0), Vfloat(8,0), Vfloat(8,0), Vfloat(8,0),  "PT", "NO_CH")<<endl; 
  }
  

  cout<<"TESTING ERROR PROPAGATION ON C9"<<endl;
  //testing error propagation on C9;

  for(int ixg=0;ixg<100;ixg++) {

    double xg =0.42 - ixg*(0.42/99);
    double q2=MBs*MBs*(1-xg);

    distr_t C9_RE_jack(1), C9_IM_jack(1);
    distr_t C9_RE_boot(0), C9_IM_boot(0);

    GaussianMersenne GPAR(854135);
    UniformMersenne  RAN(45252, 0.0, 2*M_PI);
    int NJJ=500;
    for(int ijack=0;ijack<NJJ;ijack++) {
    
      Vfloat Th(8), kappa(8), Bs_F(8), Gammas_F(8); 
     
     
      for(int j=0; j < 8; j++) {
	double unif=RAN();  {Th[j] = unif/sqrt(NJJ-1.0);}
	double rn1=GPAR();   {kappa[j] = 0.75 + 0.75*rn1/sqrt(NJJ-1.0);}
	double rn2=GPAR();  {Bs_F[j] = 1.0*rn2/sqrt(NJJ-1.0) ;}
	double rn3=GPAR();  {Gammas_F[j] = 1.0*rn3/sqrt(NJJ-1.0) ;}
      }
      
      C9_RE_jack.distr.push_back( WC(9, q2, kappa, Th, Bs_F, Gammas_F).first);
      C9_IM_jack.distr.push_back( WC(9, q2, kappa, Th, Bs_F, Gammas_F).second);
      
    }

    int NNboots=10000;

    for(int iboot=0;iboot<NNboots;iboot++) {

      Vfloat Th(8), kappa(8), Bs_F(8), Gammas_F(8); 
      for(int j=0; j < 8; j++) {

	double unif=RAN();  {Th[j] = unif;}
	double rn1=GPAR();  {kappa[j] = 0.75 + 0.75*rn1;}
	double rn2=GPAR();  {Bs_F[j] = 1.0*rn2 ;}
	double rn3=GPAR();   {Gammas_F[j] = 1.0*rn3 ;}
      }
      
      C9_RE_boot.distr.push_back( WC(9, q2, kappa, Th, Bs_F, Gammas_F).first);
      C9_IM_boot.distr.push_back( WC(9, q2, kappa, Th, Bs_F, Gammas_F).second);

    }


    double C9_RE_val= C9_RE_boot.ave(); double C9_RE_err= C9_RE_boot.err();
    double C9_IM_val= C9_IM_boot.ave(); double C9_IM_err= C9_IM_boot.err();

    distr_t C9_RE_resampled(1), C9_IM_resampled(1);
    for(int i=0;i<NJ;i++) { C9_RE_resampled.distr.push_back( C9_RE_val + GPAR()*C9_RE_err/sqrt(NJ-1.0)); C9_IM_resampled.distr.push_back( C9_IM_val + GPAR()*C9_IM_err/sqrt(NJ-1.0)); }

    if(ixg==50) { Print_To_File({}, {C9_RE_boot.distr, C9_IM_boot.distr}, "test_C9_distr", "", "");}

    cout<<"xg: "<<xg<<endl;
    cout<<"C9(RE): "<<C9_RE_jack.ave()<<" +- "<<C9_RE_jack.err()<<" [jack] , "<<C9_RE_boot.ave()<<" +- "<<C9_RE_boot.err()<<" [boot] "<<C9_RE_resampled.ave()<<" +- "<<C9_RE_resampled.err()<<" [jack resampled]"<<endl;
    cout<<"C9(IM): "<<C9_IM_jack.ave()<<" +- "<<C9_IM_jack.err()<<" [jack] , "<<C9_IM_boot.ave()<<" +- "<<C9_IM_boot.err()<<" [boot] "<<C9_IM_resampled.ave()<<" +- "<<C9_IM_resampled.err()<<" [jack resampled]"<<endl;
    cout<<"######"<<endl;
    
  }
  

  



  
  Get_Bs_xg_to_spline();
  Get_Bs_xg_to_spline_VMD();
  Get_Bs_lattice_spacings_to_print();

  GaussianMersenne GM(61111223);

  distr_t mc_MS_distr(1);
  double Lambda_MS_bar=Get_Lambda_MS_bar(4);
  cout<<"mb(mb): "<<m_MS_bar_m(3, 4, Lambda_MS_bar, 4.683 )<<endl;
  for(int i=0;i<NJ;i++) mc_MS_distr.distr.push_back( mc_MS_bar_2_ave+ GM()*mc_MS_bar_2_err/sqrt(NJ-1.0));
  cout<<"Lambda_MS_bar: "<<Lambda_MS_bar<<endl;

  double mcmc= m_MS_bar_m( 3.0,4.0, Lambda_QCD,mc_MS_distr.ave());

  if(Generate_data_for_mass_spline) {

    //generate file with m(MS-bar, 3GeV), m(MS-bar,m), m_MRS
    Vfloat m_MS_3Gev_list(3001); //between 1 and ~30 GeV
    Vfloat m_MS_m_list(3001);
    Vfloat m_MRS_list(3001);
    #pragma omp parallel for schedule(dynamic)
    for(int i=0;i<=3000;i++) {
      cout<<i<<endl<<flush;
      m_MS_3Gev_list[i] = ( 0.5 + i*0.01);
      pair<double, double > mmu_M_MRS=  MS_bar_to_mm_and_MRS_mass(3.0, 4.0, Lambda_QCD, m_MS_3Gev_list[i], mcmc);
      m_MS_m_list[i] = mmu_M_MRS.first;
      m_MRS_list[i] =  mmu_M_MRS.second;
    }

    //print to File
    Print_To_File({}, { m_MS_3Gev_list, m_MS_m_list, m_MRS_list }, "../data/ph_emission/rph/Bs_extr/MS_bar_masses.txt", "", "");
  }


  //load data for mass spline
  Vfloat m_MS_3GeV_list = Read_From_File("../data/ph_emission/rph/Bs_extr/MS_bar_masses.txt", 1, 4);
  Vfloat m_MS_m_list = Read_From_File("../data/ph_emission/rph/Bs_extr/MS_bar_masses.txt", 2, 4);
  Vfloat m_MRS_list= Read_From_File("../data/ph_emission/rph/Bs_extr/MS_bar_masses.txt", 3, 4);

  //derivatives
  Vfloat m_MS_m_list_der_MS_3GeV_list(m_MS_3GeV_list.size());
  Vfloat m_MRS_list_der_MS_3GeV_list(m_MS_3GeV_list.size());

  Vfloat m_MS_3GeV_list_der_m_MS_m_list(m_MS_3GeV_list.size());
  Vfloat m_MRS_list_der_m_MS_m_list(m_MS_3GeV_list.size());

  Vfloat m_MS_3GeV_list_der_MRS_list(m_MS_3GeV_list.size());
  Vfloat m_MS_m_list_der_MRS_list(m_MS_3GeV_list.size());


  //compute derivatives
  for(int i=0;i<(signed)m_MS_3GeV_list.size();i++) {

    if(i==0) { //forward der
      m_MS_m_list_der_MS_3GeV_list[i] = (m_MS_m_list[i+1] - m_MS_m_list[i])/(m_MS_3GeV_list[i+1] - m_MS_3GeV_list[i]);
      m_MRS_list_der_MS_3GeV_list[i] = (m_MRS_list[i+1] - m_MRS_list[i])/(m_MS_3GeV_list[i+1] - m_MS_3GeV_list[i]);

      m_MS_3GeV_list_der_m_MS_m_list[i] = (m_MS_3GeV_list[i+1] - m_MS_3GeV_list[i])/(m_MS_m_list[i+1] - m_MS_m_list[i]);
      m_MRS_list_der_m_MS_m_list[i] = (m_MRS_list[i+1] - m_MRS_list[i])/(m_MS_m_list[i+1] - m_MS_m_list[i]);

      m_MS_3GeV_list_der_MRS_list[i] = (m_MS_3GeV_list[i+1] - m_MS_3GeV_list[i])/(m_MRS_list[i+1] - m_MRS_list[i]);
      m_MS_m_list_der_MRS_list[i] = (m_MS_m_list[i+1] - m_MS_m_list[i])/(m_MRS_list[i+1] - m_MRS_list[i]);


    }
    else if(i == (signed)m_MS_3GeV_list.size() -1 ) {  //backward der
      m_MS_m_list_der_MS_3GeV_list[i] = (m_MS_m_list[i] - m_MS_m_list[i-1])/(m_MS_3GeV_list[i] - m_MS_3GeV_list[i-1]);
      m_MRS_list_der_MS_3GeV_list[i] = (m_MRS_list[i] - m_MRS_list[i-1])/(m_MS_3GeV_list[i] - m_MS_3GeV_list[i-1]);

      m_MS_3GeV_list_der_m_MS_m_list[i] = (m_MS_3GeV_list[i] - m_MS_3GeV_list[i-1])/(m_MS_m_list[i] - m_MS_m_list[i-1]);
      m_MRS_list_der_m_MS_m_list[i] = (m_MRS_list[i] - m_MRS_list[i-1])/(m_MS_m_list[i] - m_MS_m_list[i-1]);

      m_MS_3GeV_list_der_MRS_list[i] = (m_MS_3GeV_list[i] - m_MS_3GeV_list[i-1])/(m_MRS_list[i] - m_MRS_list[i-1]);
      m_MS_m_list_der_MRS_list[i] = (m_MS_m_list[i] - m_MS_m_list[i-1])/(m_MRS_list[i] - m_MRS_list[i-1]);
      

    }
    else {  //symmetric der
      m_MS_m_list_der_MS_3GeV_list[i] = (m_MS_m_list[i+1] - m_MS_m_list[i-1])/(m_MS_3GeV_list[i+1] - m_MS_3GeV_list[i-1]);
      m_MRS_list_der_MS_3GeV_list[i] = (m_MRS_list[i+1] - m_MRS_list[i-1])/(m_MS_3GeV_list[i+1] - m_MS_3GeV_list[i-1]);

      m_MS_3GeV_list_der_m_MS_m_list[i] = (m_MS_3GeV_list[i+1] - m_MS_3GeV_list[i-1])/(m_MS_m_list[i+1] - m_MS_m_list[i-1]);
      m_MRS_list_der_m_MS_m_list[i] = (m_MRS_list[i+1] - m_MRS_list[i-1])/(m_MS_m_list[i+1] - m_MS_m_list[i-1]);

      m_MS_3GeV_list_der_MRS_list[i] = (m_MS_3GeV_list[i+1] - m_MS_3GeV_list[i-1])/(m_MRS_list[i+1] - m_MRS_list[i-1]);
      m_MS_m_list_der_MRS_list[i] = (m_MS_m_list[i+1] - m_MS_m_list[i-1])/(m_MRS_list[i+1] - m_MRS_list[i-1]);
    }

  }

  


  Vfloat m_MS_3GeV_list_a = m_MS_3GeV_list;
  Vfloat m_MS_3GeV_list_b = m_MS_3GeV_list;
  Vfloat m_MS_3GeV_list_c = m_MS_3GeV_list;
  Vfloat m_MS_3GeV_list_d = m_MS_3GeV_list;


  Vfloat m_MS_m_list_a = m_MS_m_list;
  Vfloat m_MS_m_list_b = m_MS_m_list;
  Vfloat m_MS_m_list_c = m_MS_m_list;
  Vfloat m_MS_m_list_d = m_MS_m_list;

  Vfloat m_MRS_list_a = m_MRS_list;
  Vfloat m_MRS_list_b = m_MRS_list;
  Vfloat m_MRS_list_c = m_MRS_list;
  Vfloat m_MRS_list_d = m_MRS_list; 
  
  //interpolate
  //boost
  //cubic hermite interpolator
  auto MS_bar_to_mm = boost::math::interpolators::cubic_hermite( move(m_MS_3GeV_list_a), move(m_MS_m_list_a), move(m_MS_m_list_der_MS_3GeV_list));
  auto MS_bar_to_MRS= boost::math::interpolators::cubic_hermite( move(m_MS_3GeV_list_b), move(m_MRS_list_a), move(m_MRS_list_der_MS_3GeV_list));
  auto mm_to_MS_bar= boost::math::interpolators::cubic_hermite(move(m_MS_m_list_b), move(m_MS_3GeV_list_c), move(m_MS_3GeV_list_der_m_MS_m_list));
  auto mm_to_MRS= boost::math::interpolators::cubic_hermite(move(m_MS_m_list_c), move(m_MRS_list_b), move(m_MRS_list_der_m_MS_m_list));
  auto MRS_to_MS_bar= boost::math::interpolators::cubic_hermite(move(m_MRS_list_c), move(m_MS_3GeV_list_d), move(m_MS_3GeV_list_der_MRS_list));
  auto MRS_to_mm= boost::math::interpolators::cubic_hermite(move(m_MRS_list_d), move(m_MS_m_list_d), move(m_MS_m_list_der_MRS_list));

  /*
  //check approximations
  Vfloat m_MS_bar_check(100);
  for(int i=0;i<100;i++) m_MS_bar_check[i] = 1.0 + i*0.4653/10 ;
  for(auto &m_MS:m_MS_bar_check) {
    cout<<"------------------------------"<<endl;
    cout<<"m(MS-bar): "<<m_MS<<" m(m, EXACT): "<<m_MS_bar_m( 3.0, 4, Lambda_QCD, m_MS )<<" m(m,INT): "<<MS_bar_to_mm(m_MS)<<endl;
    cout<<"m(MS-bar): "<<m_MS<<" MRS(EXACT): "<<MS_bar_to_MRS_mass( 3.0, 4, Lambda_QCD, m_MS, mc_MS_distr.ave() )<<" MRS(INT): "<<MS_bar_to_MRS(m_MS)<<endl;
    cout<<"Inverse relations: "<<endl;
    cout<<"MS-bar -> mm -> MS-bar: "<<mm_to_MS_bar(MS_bar_to_mm(m_MS))<<endl;
    cout<<"MS-bar -> MRS -> MS-bar: "<<MRS_to_MS_bar(MS_bar_to_MRS(m_MS))<<endl;
    cout<<"MS-bar -> MRS -> mm -> MS-bar: "<<mm_to_MS_bar( MRS_to_mm( MS_bar_to_MRS( m_MS)))<<endl;
    cout<<"-------------------------------"<<endl;     
    } */


  
  for(int r=0; r<(signed)ratio_mh.size(); r++)  cout<<"POLE MASS for m[MS-bar, 3GeV]: "<<mc_MS_distr.ave()/ratio_mh[r]<<" GeV, "<< MS_bar_to_pole_mass_bis( 3.0 , 4.0, Lambda_MS_bar, mc_MS_distr.ave()/ratio_mh[r], mc_MS_distr.ave() )<<" [MRS]: "<<MS_bar_to_MRS(mc_MS_distr.ave()/ratio_mh[r])<<" m[MS-bar, m]: "<<MS_bar_to_mm( mc_MS_distr.ave()/ratio_mh[r])<<endl;


  cout<<"m_MRS(b)/mbmb: "<<mm_to_MRS(4.203)/4.203<<endl;
  cout<<"J(mb): "<<J_MRS(3, 4.203)<<endl;
  cout<<"J(mc): "<<J_MRS(3, 1.280)<<" mc(mc): "<<MS_bar_to_mm(mc_MS_distr.ave())<<endl;
  cout<<"L(MS-BAR): "<<"Nf=4, "<<Get_Lambda_MS_bar(4)<<" Nf=3: "<<Get_Lambda_MS_bar(3)<<endl;
  cout<<"alpha(mc(mc)): "<<Get_4l_alpha_s(m_MS_bar_m(3,4,Lambda_MS_bar,mc_MS_distr.ave()), 4, Lambda_MS_bar)<<" [Nf=4], "<<Get_4l_alpha_s(m_MS_bar_m(3,4,Lambda_MS_bar,mc_MS_distr.ave()), 3, Get_Lambda_MS_bar(3))<<endl;
  cout<<"alpha(mb(mb)): "<<Get_4l_alpha_s(m_MS_bar_m(3,4,Lambda_MS_bar,4.683), 4, Lambda_MS_bar)<<" [Nf=4], "<<Get_4l_alpha_s(m_MS_bar_m(3,4,Lambda_MS_bar,4.683), 3, Get_Lambda_MS_bar(3))<<endl;
  cout<<"alpha(5 GeV): "<<Get_4l_alpha_s(5, 4, Lambda_MS_bar)<<" [Nf=4], "<<Get_4l_alpha_s(4, 3, Get_Lambda_MS_bar(3))<<endl;


  cout<<"Comparing evolutor (mu1=4.203, mu2=3): EXACT: "<<MS_bar_mass_evolutor(4.203, 3, 4, Lambda_MS_bar, 0)<<" , APPROXIMATE: "<<MS_bar_mass_evolutor(4.203, 3, 4, Lambda_MS_bar, 1)<<endl;

  cout<<"mb(5 GeV): "<<4.203*MS_bar_mass_evolutor(5, 4.203, 4, Lambda_MS_bar, 0)<<endl;
 

  //checking MRS calculation
  cout<<"Checking calculation of MRS masses"<<endl;
  cout<<"mc(mc, FNAL): 1.273 [GeV] -> mc(MRS, FNAL): 1.392 [GeV] . This work mc(MRS): "<<MS_bar_to_MRS_mass(1.273, 4, Lambda_MS_bar, 1.273, 1.273)<<" [GeV]"<<endl;
  cout<<"mb(mb, FNAL): 4.201 [GeV] -> mc(MRS, FNAL): 4.749 [GeV] . This work mc(MRS): "<<MS_bar_to_MRS_mass(4.201, 4, Lambda_MS_bar, 4.201, 1.273)<<" [GeV]"<<endl;

 
  //determine the form factor corresponding to the emission of the real photon from O7

  vector<rt_07_Bs> TFF_virtual_ret_list;

  if(!Skip_virtual_diagram) {
  
    vector<bool> UseJack_list_07({1,1,1});
    vector<int>  num_xg_list_07({4,4,4});
    vector<string> Corr_path_list_07({ "../Bs_mumu_gamma_data/07_virtual/mh0", "../Bs_mumu_gamma_data/07_virtual/mh1", "../Bs_mumu_gamma_data/07_virtual/mh3"  });
    vector<string> out_tag_list_07({"mh0", "mh1", "mh3"});
    vector<string> Meson_list_07({"B0s", "B1s", "B3s"});
 
  
  int NN= UseJack_list_07.size();
  for(int i=0;i<NN;i++) {
    MESON=Meson_list_07[i];
    //create output directories
    boost::filesystem::create_directory("../data/ph_emission");
    boost::filesystem::create_directory("../data/ph_emission/"+ph_type);
    boost::filesystem::create_directory("../data/ph_emission/"+ph_type+"/"+MESON);
    boost::filesystem::create_directory("../data/ph_emission/"+ph_type+"/"+MESON+"/07_virtual");
    
    string path_out="../data/ph_emission/"+ph_type+"/"+MESON+"/07_virtual";
    
    TFF_virtual_ret_list.push_back(Get_virtual_tensor_FF(num_xg_list_07[i], UseJack_list_07[i], NJ, MESON, Corr_path_list_07[i],  path_out));
    
  }
  
  }

 



      

  vector<rt_FF_Bs> FF_ret_list;
  vector<rt_FF_Bs> FF_ret_list_boot;
 
  /*
  vector<bool> Use_three_finest_list({0,0,0,0,0,1,1,1,1,1,0,0,0,0,0});
  vector<bool> Include_a4_list({0,0,0,0,0,0,0,0,0,0,1,1,1,1,1});
  vector<bool> UseJack_list({1,1,1,1,1,1,1,1,1,1,1,1,1,1,1});
  vector<int>  num_xg_list({5,5,5,5,5,5,5,5,5,5,5,5,5,5,5});
  vector<double> Perform_continuum_extrapolation_list({1,1,1,1,1,1,1,1,1,1,1,1,1,1,1});
  vector<string> Corr_path_list({ "../Bs_mumu_gamma_data/mh0",  "../Bs_mumu_gamma_data/mh1", "../Bs_mumu_gamma_data/mh2", "../Bs_mumu_gamma_data/mh3", "../Bs_mumu_gamma_data/mh4", "../Bs_mumu_gamma_data/mh0",  "../Bs_mumu_gamma_data/mh1", "../Bs_mumu_gamma_data/mh2", "../Bs_mumu_gamma_data/mh3", "../Bs_mumu_gamma_data/mh4",  "../Bs_mumu_gamma_data/mh0",  "../Bs_mumu_gamma_data/mh1", "../Bs_mumu_gamma_data/mh2", "../Bs_mumu_gamma_data/mh3", "../Bs_mumu_gamma_data/mh4"});
  vector<string> out_tag_list({"mh0", "mh1", "mh2", "mh3", "mh4", "mh0", "mh1", "mh2", "mh3", "mh4",  "mh0", "mh1", "mh2", "mh3", "mh4"});
  vector<string> Meson_list({"B0s", "B1s", "B2s", "B3s", "B4s", "B0s", "B1s", "B2s", "B3s", "B4s", "B0s", "B1s", "B2s", "B3s", "B4s"});
  */

  
  vector<bool> Use_three_finest_list({0,0,0,0,0,1,1,1,1,1});
  vector<bool> Include_a4_list({0,0,0,0,0,0,0,0,0,0});
  vector<bool> UseJack_list({1,1,1,1,1,1,1,1,1,1});
  vector<bool> UseJack_list_boot({0,0,0,0,0,0,0,0,0,0});
  vector<int>  num_xg_list({5,5,5,5,5,5,5,5,5,5});
  vector<double> Perform_continuum_extrapolation_list({1,1,1,1,1,1,1,1,1,1});
  vector<string> Corr_path_list({ "../Bs_mumu_gamma_data/mh0","../Bs_mumu_gamma_data/mh1", "../Bs_mumu_gamma_data/mh2", "../Bs_mumu_gamma_data/mh3", "../Bs_mumu_gamma_data/mh4", "../Bs_mumu_gamma_data/mh0",  "../Bs_mumu_gamma_data/mh1", "../Bs_mumu_gamma_data/mh2", "../Bs_mumu_gamma_data/mh3", "../Bs_mumu_gamma_data/mh4"});
  vector<string> out_tag_list({"mh0", "mh1", "mh2", "mh3", "mh4", "mh0", "mh1", "mh2", "mh3", "mh4"});
  vector<string>Meson_list({"B0s", "B1s", "B2s", "B3s", "B4s", "B0s", "B1s", "B2s", "B3s", "B4s"});
  
  
  /*
  vector<bool> Use_three_finest_list({0,0,0,0,1,1,1,1,0,0,0,0});
  vector<bool> Include_a4_list({0,0,0,0,0,0,0,0,1,1,1,1});
  vector<bool> UseJack_list({1,1,1,1,1,1,1,1,1,1,1,1});
  vector<int>  num_xg_list({5,5,5,5,5,5,5,5,5,5,5,5});
  vector<double> Perform_continuum_extrapolation_list({1,1,1,1,1,1,1,1,1,1,1,1});
  vector<string> Corr_path_list({ "../Bs_mumu_gamma_data/mh1",  "../Bs_mumu_gamma_data/mh2", "../Bs_mumu_gamma_data/mh3", "../Bs_mumu_gamma_data/mh4", "../Bs_mumu_gamma_data/mh1",  "../Bs_mumu_gamma_data/mh2", "../Bs_mumu_gamma_data/mh3", "../Bs_mumu_gamma_data/mh4",  "../Bs_mumu_gamma_data/mh1",  "../Bs_mumu_gamma_data/mh2", "../Bs_mumu_gamma_data/mh3", "../Bs_mumu_gamma_data/mh4"});
  vector<string> out_tag_list({"mh1", "mh2", "mh3", "mh4", "mh1", "mh2", "mh3", "mh4",  "mh1", "mh2", "mh3", "mh4"});
  vector<string> Meson_list({"B1s", "B2s", "B3s", "B4s",  "B1s", "B2s", "B3s", "B4s",  "B1s", "B2s", "B3s", "B4s"});
  */

  
  int N= UseJack_list.size();
  for(int i=0;i<N;i++) {
    Bs_xg_t_list.clear();
    
    Get_Bs_xg_t_list(num_xg_list[i]);
    MESON=Meson_list[i];
    //create output directories
    boost::filesystem::create_directory("../data/ph_emission");
    boost::filesystem::create_directory("../data/ph_emission/"+ph_type);
    boost::filesystem::create_directory("../data/ph_emission/"+ph_type+"/"+MESON);
    boost::filesystem::create_directory("../data/ph_emission/"+ph_type+"/"+MESON+"/mass");
    boost::filesystem::create_directory("../data/ph_emission/"+ph_type+"/"+MESON+"/decay_const");
    boost::filesystem::create_directory("../data/ph_emission/"+ph_type+"/"+MESON+"/decay_const/corr_matrix_fit");
    boost::filesystem::create_directory("../data/ph_emission/"+ph_type+"/"+MESON+"/C");
    boost::filesystem::create_directory("../data/ph_emission/"+ph_type+"/"+MESON+"/H");
    boost::filesystem::create_directory("../data/ph_emission/"+ph_type+"/"+MESON+"/C_T");
    boost::filesystem::create_directory("../data/ph_emission/"+ph_type+"/"+MESON+"/H_T");
    boost::filesystem::create_directory("../data/ph_emission/"+ph_type+"/"+MESON+"/FF");
    boost::filesystem::create_directory("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum");
    boost::filesystem::create_directory("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/BMA");
    boost::filesystem::create_directory("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/BMA/uniform");
    boost::filesystem::create_directory("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/BMA/AIC"); 
    boost::filesystem::create_directory("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/per_kin");
    boost::filesystem::create_directory("../data/ph_emission/"+ph_type+"/"+MESON+"/FF_u");
    boost::filesystem::create_directory("../data/ph_emission/"+ph_type+"/"+MESON+"/FF_u/per_kin");
    boost::filesystem::create_directory("../data/ph_emission/"+ph_type+"/"+MESON+"/FF_d");
    boost::filesystem::create_directory("../data/ph_emission/"+ph_type+"/"+MESON+"/FF_d/per_kin");

    
    string Fit_tag= ( (Include_a4_list[i]==true)?"wa4_":"");
    Fit_tag = ( (Use_three_finest_list[i]==true)?"wtf_":Fit_tag);


    //Compute form factors
    if(Compute_FF) {
      FF_ret_list_boot.push_back(Get_Bs_mumu_gamma_form_factors(num_xg_list[i], Perform_continuum_extrapolation_list[i], Use_three_finest_list[i], Include_a4_list[i], UseJack_list_boot[i], Fit_tag, Corr_path_list[i], out_tag_list[i])); 
      FF_ret_list.push_back(Get_Bs_mumu_gamma_form_factors(num_xg_list[i], Perform_continuum_extrapolation_list[i], Use_three_finest_list[i], Include_a4_list[i], UseJack_list[i], Fit_tag, Corr_path_list[i], out_tag_list[i]));
    }
    //or read them from file
    else {
      FF_ret_list.emplace_back("../data/ph_emission/"+ph_type+"/"+MESON, UseJack_list[i], Use_three_finest_list[i], Include_a4_list[i], num_xg_list[i] );
      FF_ret_list_boot.emplace_back("../data/ph_emission/"+ph_type+"/"+MESON, UseJack_list_boot[i], Use_three_finest_list[i], Include_a4_list[i], num_xg_list[i] );
    }

        
    
  }

  
 
  //Fit F_Bs

  
  vector<distr_t_list> FF_Bs_list(6, 0.0*Get_id_jack_distr_list(4, NJ));  // A, V, T_A, T_V, T , V
  
  
  

  
  
  auto C_HQET = [&](double x) {


    

    double L=Lambda_MS_bar;
    double Nf=4.0;
    double g0 = -4.0;
    double g1 = -254.0/9 -56.0*M_PI*M_PI/27 +20.0*Nf/9;
    double B0 = 11.0 -2.0*Nf/3;
    double B1= 102.0 - 38.0*Nf/3.0;
    double alpha = Get_4l_alpha_s(x,Nf,L);


    double Cmm= 1 -2.0*alpha/(3*M_PI) -(4.20 -0.44*(Nf-1)  +0.07)*pow(alpha/M_PI,2);
    //three-loops anomalous dimensions of HQET currents from https://arxiv.org/pdf/hep-ph/0303113.pdf
    double z3= riemann_zeta(3);
    double G0= -1.0;
    double G1= (-3.04328 + 0.138889*(Nf));
    double G2= (-12.941 + 1.55406*(Nf) + 0.0270062*pow(Nf,2));
    double b0= (33.0-2.0*Nf)/(12.0);
    double b1= (153.0 - 19*Nf)/(24.0);
    double b2= (2857.0 -5033.0*Nf/9.0 + 325.0*Nf*Nf/27.0)/(128.0);
    double b3= ( (149753.0/6.0 + 3564*z3) -(1078361.0/162.0 + 6508.0*z3/27.0)*Nf + (50065.0/162.0 + 6472.0*z3/81.0)*Nf*Nf + 1093.0*pow(Nf,3)/729.0)/(256.0);
  
    auto F= [&](double x) { return 0.5*(G0*x + G1*pow(x,2)+ G2*pow(x,3))/(x*( b0*x + b1*pow(x,2) + b2*pow(x,3) + b3*pow(x,4)));};  // F = g/B
    
    gsl_function_pp<decltype(F)> Fgsl(F);
    gsl_integration_workspace * w = gsl_integration_workspace_alloc (10000);
    gsl_function *G = static_cast<gsl_function*>(&Fgsl);
    double val, err;
    gsl_integration_qags(G, 0.5, Get_4l_alpha_s(x, Nf, L)/M_PI , 0.0, 1e-4, 10000, w, &val, &err);
    gsl_integration_workspace_free (w);

    return Cmm*exp(val);
    
     
    //return pow(alpha, g0/(2*B0) )*( 1 + (alpha/(4.0*M_PI))*( -8.0/3 + g1/(2*B0) -g0*B1/(2*B0*B0)));
    
  };

  
  distr_t_list masses(true), fps(true), phis(true);
  distr_t_list masses_boot(false), fps_boot(false), phis_boot(false);
  int FF_SIZE= (signed)FF_ret_list.size();
  for(int m=0; m< Nmasses ; m++) {

    //jackknife
    vector<distr_t> VAL_MASS, VAL_PHI;
    vector<double> CH2_MASS, CH2_PHI;
    vector<int> Ndof_MASS, Ndof_PHI, Nmeas_MASS, Nmeas_PHI;

    int k=0;
    while( k*Nmasses + m < FF_SIZE) {
      int r= k*Nmasses+m;
      VAL_MASS.push_back( FF_ret_list[r].mass);
      VAL_PHI.push_back( FF_ret_list[r].phi[0]);
      VAL_PHI.push_back( FF_ret_list[r].phi[1]);
      VAL_PHI.push_back( FF_ret_list[r].phi[2]);
      CH2_MASS.push_back( FF_ret_list[r].Ch2_mass*FF_ret_list[r].Ndof);
      CH2_PHI.push_back( FF_ret_list[r].Ch2_phi[0]*FF_ret_list[r].Ndof_K[0]);
      CH2_PHI.push_back( FF_ret_list[r].Ch2_phi[1]*FF_ret_list[r].Ndof_K[1]);
      CH2_PHI.push_back( FF_ret_list[r].Ch2_phi[2]*FF_ret_list[r].Ndof_K[2]);
      Ndof_MASS.push_back( FF_ret_list[r].Ndof);
      Ndof_PHI.push_back( FF_ret_list[r].Ndof_K[0]);
      Ndof_PHI.push_back( FF_ret_list[r].Ndof_K[1]);
      Ndof_PHI.push_back( FF_ret_list[r].Ndof_K[2]);
      Nmeas_MASS.push_back( FF_ret_list[r].Nmeas);
      Nmeas_PHI.push_back( FF_ret_list[r].Nmeas_K[0]);
      Nmeas_PHI.push_back( FF_ret_list[r].Nmeas_K[1]);
      Nmeas_PHI.push_back( FF_ret_list[r].Nmeas_K[2]);
      
      k++;
    }

    masses.distr_list.push_back( AIC_lin( VAL_MASS, CH2_MASS, Ndof_MASS, Nmeas_MASS,0));
    phis.distr_list.push_back( AIC_lin( VAL_PHI, CH2_PHI, Ndof_PHI, Nmeas_PHI,0));
    ofstream Print_phi_AIC("../data/ph_emission/"+ph_type+"/Bs_extr/phi_m_"+to_string(m)+".dat");
    Print_phi_AIC<<phis[m].ave()<<" "<<phis[m].err()<<endl;
    Print_phi_AIC.close();
    fps.distr_list.push_back( phis.distr_list[m]/SQRT_D(masses.distr_list[m]) );


    //bootstrap
    vector<distr_t> VAL_MASS_boot, VAL_PHI_boot;
    vector<double> CH2_MASS_boot, CH2_PHI_boot;
    vector<int> Ndof_MASS_boot, Ndof_PHI_boot, Nmeas_MASS_boot, Nmeas_PHI_boot;

    k=0;
    while( k*Nmasses + m < FF_SIZE) {
      int r= k*Nmasses+m;
      VAL_MASS_boot.push_back( FF_ret_list_boot[r].mass);
      VAL_PHI_boot.push_back( FF_ret_list_boot[r].phi[0]);
      VAL_PHI_boot.push_back( FF_ret_list_boot[r].phi[1]);
      VAL_PHI_boot.push_back( FF_ret_list_boot[r].phi[2]);
      CH2_MASS_boot.push_back( FF_ret_list_boot[r].Ch2_mass*FF_ret_list_boot[r].Ndof);
      CH2_PHI_boot.push_back( FF_ret_list_boot[r].Ch2_phi[0]*FF_ret_list_boot[r].Ndof_K[0]);
      CH2_PHI_boot.push_back( FF_ret_list_boot[r].Ch2_phi[1]*FF_ret_list_boot[r].Ndof_K[1]);
      CH2_PHI_boot.push_back( FF_ret_list_boot[r].Ch2_phi[2]*FF_ret_list_boot[r].Ndof_K[2]);
      Ndof_MASS_boot.push_back( FF_ret_list_boot[r].Ndof);
      Ndof_PHI_boot.push_back( FF_ret_list_boot[r].Ndof_K[0]);
      Ndof_PHI_boot.push_back( FF_ret_list_boot[r].Ndof_K[1]);
      Ndof_PHI_boot.push_back( FF_ret_list_boot[r].Ndof_K[2]);
      Nmeas_MASS_boot.push_back( FF_ret_list_boot[r].Nmeas);
      Nmeas_PHI_boot.push_back( FF_ret_list_boot[r].Nmeas_K[0]);
      Nmeas_PHI_boot.push_back( FF_ret_list_boot[r].Nmeas_K[1]);
      Nmeas_PHI_boot.push_back( FF_ret_list_boot[r].Nmeas_K[2]);
      
      k++;
    }

    masses_boot.distr_list.push_back( AIC_lin( VAL_MASS_boot, CH2_MASS_boot, Ndof_MASS_boot, Nmeas_MASS_boot,0));
    phis_boot.distr_list.push_back( AIC_lin( VAL_PHI_boot, CH2_PHI_boot, Ndof_PHI_boot, Nmeas_PHI_boot,0));
    fps_boot.distr_list.push_back( phis_boot.distr_list[m]/SQRT_D(masses_boot.distr_list[m]) );
    
  }

 
  //extrapolate to the Bs-meson mass
  //start from the ratio mh/mc
  cout<<"Determine mb/mc..."<<endl;
  
  class ipar_mh {
  public:
    ipar_mh() : FF(0.0), FF_err(0.0) {}
    double FF, FF_err, mm, r, alpha, mc;
  };
  
  class fpar_mh {
  public:
    fpar_mh() {}
    fpar_mh(const Vfloat &par) {
      if((signed)par.size() != 4) crash("In class fpar_mh  class constructor Vfloat par has size != 4");
      R=par[0];
      A=par[1];
      B=par[2];
      C=par[3];
    }
    double R,A,B,C;
  };

  //init bootstrap fit
  bootstrap_fit<fpar_mh,ipar_mh> bf_mh(NJ);
  bf_mh.set_warmup_lev(1); //sets warmup
  bf_mh.Set_number_of_measurements(Nmasses);
  bf_mh.Set_verbosity(1);
  

  bf_mh.Add_par("R", 1.0, 0.1);
  bf_mh.Add_par("A", 1.0, 0.1);
  bf_mh.Add_par("B", 1.0 , 0.1);
  bf_mh.Add_par("C", 1.0, 0.1);

   //fit on mean values to get ch2
  bootstrap_fit<fpar_mh,ipar_mh> bf_mh_ch2(1);
  bf_mh_ch2.set_warmup_lev(1); //sets warmup
  bf_mh_ch2.Set_number_of_measurements(Nmasses);
  bf_mh_ch2.Set_verbosity(1);
  bf_mh_ch2.Add_par("R", 1.0, 0.1);
  bf_mh_ch2.Add_par("A", 1.0, 0.1);
  bf_mh_ch2.Add_par("B", 1.0, 0.1);
  bf_mh_ch2.Add_par("C", 1.0, 0.1);

  distr_t mc_HQET(1);
  for(int ijack=0;ijack<NJ;ijack++) mc_HQET.distr.push_back( MS_bar_to_MRS(mc_MS_distr.distr[ijack]));

  bool Fix_R=true;
  if(Fix_R) {bf_mh.Fix_par("R", 1.0); bf_mh_ch2.Fix_par("R",1.0);}

  bool Fix_C=false;
  if(Fix_C) {bf_mh.Fix_par("C", 0.0); bf_mh_ch2.Fix_par("C",0.0);}

      
 
  
  //ansatz
  bf_mh.ansatz=  [](const fpar_mh &p, const ipar_mh &ip) {

    double log_mug= pow(ip.alpha,9.0/25)*( 1 + 0.67235*ip.alpha +  1.284*pow(ip.alpha,2));

    return ip.r*ip.mc*( p.R  +  p.A/ip.r + p.B/pow(ip.r,2) + p.C*log_mug/pow(ip.r,2));
  
    	 	  
  };

  
  bf_mh.measurement=  [ ](const fpar_mh &p, const ipar_mh &ip) {
    
    return ip.FF;
  };
  bf_mh.error=  [ ](const fpar_mh &p, const ipar_mh &ip) {
	  
    return ip.FF_err;
  };
	
  bf_mh_ch2.ansatz= bf_mh.ansatz;
  bf_mh_ch2.measurement = bf_mh.measurement;
  bf_mh_ch2.error = bf_mh.error;

  //insert covariance matrix
  //insert covariance matrix
  Eigen::MatrixXd Cov_Matrix_mh(Nmasses,Nmasses);
  Eigen::MatrixXd Corr_Matrix_mh(Nmasses,Nmasses);
  for(int i=0;i<Nmasses;i++)
    for(int j=0;j<Nmasses;j++) {
      Corr_Matrix_mh(i,j) = masses_boot.distr_list[i]%masses_boot.distr_list[j]/(masses_boot.err(i)*masses_boot.err(j));
      Cov_Matrix_mh(i,j) = Corr_Matrix_mh(i,j)*masses.err(i)*masses.err(j);
    }

  //print on screen
  cout<<"Correlation matrix mass fit"<<endl;
  cout<<Corr_Matrix_mh<<endl;
  
  //add cov matrix to bootstrap fit
  bf_mh.Add_covariance_matrix(Cov_Matrix_mh);
  bf_mh_ch2.Add_covariance_matrix(Cov_Matrix_mh);
  
       
  //start fitting
  //fill the data
  vector<vector<ipar_mh>> data_mh(NJ);
  vector<vector<ipar_mh>> data_mh_ch2(1);
  //allocate space for output result
  boot_fit_data<fpar_mh> Bt_fit_mh;
  boot_fit_data<fpar_mh> Bt_fit_mh_ch2;
  distr_t_list ratio_mh_pole(true, Nmasses, NJ);
  ratio_mh_pole.distr_list[0] = Get_id_jack_distr(NJ);
  for(auto &data_iboot: data_mh) data_iboot.resize(Nmasses);
  for(auto &data_iboot: data_mh_ch2) data_iboot.resize(Nmasses);
  for(int ijack=0;ijack<NJ;ijack++) {
    for(int im=0;im<Nmasses;im++) {
      data_mh[ijack][im].FF = (masses.distr_list[im]).distr[ijack];
      data_mh[ijack][im].FF_err= (masses.distr_list[im]).err();
      data_mh[ijack][im].mm= MS_bar_to_mm(mc_MS_distr.distr[ijack]/ratio_mh[im]);
      data_mh[ijack][im].mc= mc_HQET.distr[ijack];
      data_mh[ijack][im].r = MS_bar_to_MRS(mc_MS_distr.distr[ijack]/ratio_mh[im])/mc_HQET.distr[ijack];
      data_mh[ijack][im].alpha= Get_4l_alpha_s(data_mh[ijack][im].mm, 4, Lambda_MS_bar);
      ratio_mh_pole.distr_list[im].distr[ijack]= 1.0/data_mh[ijack][im].r;
      if(ijack==0) {
	data_mh_ch2[ijack][im].FF = (masses.distr_list[im]).ave();
	data_mh_ch2[ijack][im].FF_err= (masses.distr_list[im]).err();
	data_mh_ch2[ijack][im].mm = MS_bar_to_mm(mc_MS_distr.ave()/ratio_mh[im]);
	data_mh_ch2[ijack][im].mc = mc_HQET.ave();
	data_mh_ch2[ijack][im].r = MS_bar_to_MRS(mc_MS_distr.ave()/ratio_mh[im])/mc_HQET.ave();
	data_mh_ch2[ijack][im].alpha= Get_4l_alpha_s(data_mh_ch2[ijack][im].mm, 4, Lambda_MS_bar);
	
      }
    }
  }
  
  //append
  bf_mh.Append_to_input_par(data_mh);
  bf_mh_ch2.Append_to_input_par(data_mh_ch2);
  //fit
  cout<<"Fitting Bs-meson mass"<<endl;
  Bt_fit_mh= bf_mh.Perform_bootstrap_fit();
  Bt_fit_mh_ch2= bf_mh_ch2.Perform_bootstrap_fit();    
  int Npars= 2+ (Fix_R==true) + (Fix_C==true);
  double ch2_red_mh= Bt_fit_mh_ch2.get_ch2_ave()/( Nmasses -Npars);

 
  
  //retrieve params
  distr_t R_mh(1), A_mh(1), B_mh(1), C_mh(1);
  for(int ijack=0;ijack<NJ;ijack++) {
    R_mh.distr.push_back( Bt_fit_mh.par[ijack].R);
    A_mh.distr.push_back( Bt_fit_mh.par[ijack].A);
    B_mh.distr.push_back( Bt_fit_mh.par[ijack].B);
    C_mh.distr.push_back( Bt_fit_mh.par[ijack].C);
  }
  
  distr_t_list F_mh_fit(1);
  distr_t mb_ov_mc(1);
  distr_t mb(1);
  for(int ijack=0;ijack<NJ;ijack++) {
    auto Fmb=[&R_mh, &A_mh, &B_mh,&C_mh,  &ijack, &mc_HQET, &Lambda_MS_bar, &MRS_to_mm](double x) {

      double mm= MRS_to_mm(x);
      double a= Get_4l_alpha_s(mm, 4.0, Lambda_MS_bar);
      double log_mug= pow(a,9.0/25)*( 1 + 0.67235*a + 1.284*pow(a,2));
      //double MBs=5.36692;
      return (x/MBs)*mc_HQET.distr[ijack]*( R_mh.distr[ijack]  + A_mh.distr[ijack]*(1.0/x) + B_mh.distr[ijack]*pow(1.0/x,2) + C_mh.distr[ijack]*log_mug*pow(1.0/x,2)) -1.0 ;
      
    };
    double mb_ov_mc_pole= R_brent(Fmb,3.0, 4.0);
    double mb_pole= mb_ov_mc_pole*mc_HQET.distr[ijack];
    double mb_MS_bar_3_GeV = MRS_to_MS_bar(mb_pole);
    mb_ov_mc.distr.push_back( mb_MS_bar_3_GeV/mc_MS_distr.distr[ijack]);
  }
  mb= mb_ov_mc*mc_MS_distr;
  distr_t mbmb(1);
  for(int ijack=0;ijack<NJ;ijack++) mbmb.distr.push_back( MS_bar_to_mm(mb.distr[ijack]));
  
  cout<<"mb (3GeV): "<<mb.ave()<<" +- "<<mb.err()<<" [GeV] ,  mb/mc: "<<mb_ov_mc.ave()<<" +- "<<mb_ov_mc.err()<<" mb(mb) : "<<mbmb.ave()<<" +- "<<mbmb.err()<<endl;
  
  Vfloat r_to_print;
  for(int r=0;r<300;r++) { r_to_print.push_back( 1.0 + (10-0.9)*r/299 ) ;}
  distr_t_list r_pole_to_print; 
   
  //######### PRINT FIT FUNCTION ##########
  
  //plot fit function
  for(auto &r: r_to_print) {
    ipar_mh pp_mh;
    distr_t F_mh(1), r_pole(1);
    for(int ijack=0;ijack<NJ;ijack++) {
      pp_mh.mm= MS_bar_to_mm(mc_MS_distr.distr[ijack]*r);
      pp_mh.r=  MS_bar_to_MRS(mc_MS_distr.distr[ijack]*r)/mc_HQET.distr[ijack];
      pp_mh.mc = mc_HQET.distr[ijack];
      pp_mh.alpha= Get_4l_alpha_s(pp_mh.mm, 4.0, Lambda_MS_bar);
      F_mh.distr.push_back( bf_mh.ansatz( Bt_fit_mh.par[ijack], pp_mh));
      r_pole.distr.push_back(pp_mh.r);
    }
    
    F_mh_fit.distr_list.push_back(F_mh);
    r_pole_to_print.distr_list.push_back(r_pole);
  }
  //print
  string out_path_mhs="../data/ph_emission/"+ph_type+"/Bs_extr/mb.fit_func";
  
  Print_To_File({},{ r_to_print, r_pole_to_print.ave(), r_pole_to_print.err(), F_mh_fit.ave(), F_mh_fit.err()} , out_path_mhs , "", "ch2/dof: "+to_string_with_precision(ch2_red_mh,5));
  Print_To_File({}, {ratio_mh, ratio_mh_pole.ave(), (1.0/masses).ave(), (1.0/masses).err(), (masses/masses.distr_list[0]).ave(), (masses/masses.distr_list[0]).err()}, "../data/ph_emission/"+ph_type+"/Bs_extr/mb.dat", "", "");
  
  
  
  
  
  //Decay constant fBs
  
  
  //generate LEC lambda

  
  //extrapolate f_Bs
  class ipar_fB {
  public:
    ipar_fB() : FF(0.0), FF_err(0.0) {}
    double FF, FF_err, M, Mpole;
  };
  
  class fpar_fB {
  public:
    fpar_fB() {}
    fpar_fB(const Vfloat &par) {
      if((signed)par.size() != 3) crash("In class fpar_fB  class constructor Vfloat par has size != 3");
      A=par[0];
      B=par[1];
      C=par[2];
    }
    double A,B,C;
  };

  //init bootstrap fit
  bootstrap_fit<fpar_fB,ipar_fB> bf_fB(NJ);
  bf_fB.set_warmup_lev(1); //sets warmup
  bf_fB.Set_number_of_measurements(Nmasses);
  bf_fB.Set_verbosity(1);
  

  
  bf_fB.Add_par("A", 0.4, 0.003);
  bf_fB.Add_par("B", -0.5 , 0.01);
  bf_fB.Add_par("C", 1.0 , 0.01);
 
  //fit on mean values to get ch2
  bootstrap_fit<fpar_fB,ipar_fB> bf_fB_ch2(1);
  bf_fB_ch2.set_warmup_lev(1); //sets warmup
  bf_fB_ch2.Set_number_of_measurements(Nmasses);
  bf_fB_ch2.Set_verbosity(1);
  bf_fB_ch2.Add_par("A", 0.4, 0.003);
  bf_fB_ch2.Add_par("B", -0.5, 0.01); 
  bf_fB_ch2.Add_par("C", 1.0, 0.01);


  bool Fix_C_fp=true;
  if(Fix_C_fp) { bf_fB.Fix_par("C", 0.0); bf_fB_ch2.Fix_par("C",0.0);}
 
  
  //ansatz
  bf_fB.ansatz=  [ &C_HQET ](const fpar_fB &p, const ipar_fB &ip) {
    return C_HQET( ip.M  )*(p.A + p.B/(ip.M/0.6 ) + p.C/(pow(ip.M/0.6,2)));
  };
  
  bf_fB.measurement=  [ ](const fpar_fB &p, const ipar_fB &ip) {
    return ip.FF;
  };
  
  bf_fB.error=  [ ](const fpar_fB &p, const ipar_fB &ip) {
    return ip.FF_err;
  };
	
  bf_fB_ch2.ansatz= bf_fB.ansatz;
  bf_fB_ch2.measurement = bf_fB.measurement;
  bf_fB_ch2.error = bf_fB.error;


  //insert covariance matrix
  Eigen::MatrixXd Cov_Matrix_phi(Nmasses,Nmasses);
  Eigen::MatrixXd Corr_Matrix_phi(Nmasses,Nmasses);
  for(int i=0;i<Nmasses;i++)
    for(int j=0;j<Nmasses;j++) {
      Corr_Matrix_phi(i,j) = phis_boot.distr_list[i]%phis_boot.distr_list[j]/(phis_boot.err(i)*phis_boot.err(j));
      Cov_Matrix_phi(i,j) = Corr_Matrix_phi(i,j)*phis.err(i)*phis.err(j);
    }

  //print on screen
  cout<<"Correlation matrix decay constant fit"<<endl;
  cout<<Corr_Matrix_phi<<endl;
  
  //add cov matrix to bootstrap fit
  bf_fB.Add_covariance_matrix(Cov_Matrix_phi);
  bf_fB_ch2.Add_covariance_matrix(Cov_Matrix_phi);
  
       
  //start fitting
  //fill the data
  vector<vector<ipar_fB>> data_fB(NJ);
  vector<vector<ipar_fB>> data_fB_ch2(1);
  //allocate space for output result
  boot_fit_data<fpar_fB> Bt_fit_fB;
  boot_fit_data<fpar_fB> Bt_fit_fB_ch2;
  for(auto &data_iboot: data_fB) data_iboot.resize(Nmasses);
  for(auto &data_iboot: data_fB_ch2) data_iboot.resize(Nmasses);
  for(int ijack=0;ijack<NJ;ijack++) {
    for(int im=0;im<Nmasses;im++) {
      data_fB[ijack][im].FF = phis.distr_list[im].distr[ijack];
      data_fB[ijack][im].FF_err= phis.err(im);
      data_fB[ijack][im].M= masses.distr_list[im].distr[ijack];
      data_fB[ijack][im].Mpole= MS_bar_to_MRS(mc_MS_distr.distr[ijack]/ratio_mh[im] ); 
      if(ijack==0) {
	data_fB_ch2[ijack][im].FF = phis.ave(im);
	data_fB_ch2[ijack][im].FF_err= phis.err(im);
	data_fB_ch2[ijack][im].M = masses.ave(im);
	data_fB_ch2[ijack][im].Mpole= MS_bar_to_MRS(mc_MS_distr.ave()/ratio_mh[im] ); 
      }
    }
  }
  
  //append
  bf_fB.Append_to_input_par(data_fB);
  bf_fB_ch2.Append_to_input_par(data_fB_ch2);
  //fit
  cout<<"Fitting Bs-meson decay constant"<<endl;
  Bt_fit_fB= bf_fB.Perform_bootstrap_fit();
  Bt_fit_fB_ch2= bf_fB_ch2.Perform_bootstrap_fit();
  Npars= 3 - (Fix_C_fp==true);
  double ch2_red_fB= Bt_fit_fB_ch2.get_ch2_ave()/( Nmasses -Npars);
  
  //retrieve params
  distr_t A_p(1), B_p(1), C_p(1);
  for(int ijack=0;ijack<NJ;ijack++) {
    A_p.distr.push_back( Bt_fit_fB.par[ijack].A);
    B_p.distr.push_back( Bt_fit_fB.par[ijack].B);
    C_p.distr.push_back( Bt_fit_fB.par[ijack].C);
  }

  //print fit function
  distr_t_list F_fB_fit(1);
  distr_t_list C_HQET_fit(1);
  distr_t f_Bs(1);
  Vfloat Inv_Masses_to_print;
  for(int i=0;i<100;i++) { Inv_Masses_to_print.push_back( 0.1 + 0.25*2*i/99  ) ;}

  //### Determine Bs decay constant
  ipar_fB pp_Bs;
  //double MBs=5.36692;
  pp_Bs.M= MBs;
  for(int ijack=0;ijack<NJ;ijack++) {
    pp_Bs.Mpole= MS_bar_to_MRS(mb.distr[ijack] );
    f_Bs.distr.push_back( bf_fB.ansatz( Bt_fit_fB.par[ijack], pp_Bs)/sqrt(MBs));
  }
  //#################################


  cout<<"fBs: "<<f_Bs.ave()<<" +- "<<f_Bs.err()<<" [GeV] "<<endl;

  if(plot_fit_func_fBs) {
  
  //######### PRINT FIT FUNCTION ##########
  //plot fit function
  for(auto &Minv: Inv_Masses_to_print) {
    ipar_fB pp_fB;
    pp_fB.M=1.0/Minv;
    distr_t F_fB_M(1);
    distr_t C_HQET_M(1);
    for(int ijack=0;ijack<NJ;ijack++) {

     
      auto Fmh=[&R_mh, &A_mh, &B_mh, &C_mh, &ijack, &Minv, &mc_HQET, &MRS_to_mm, &Lambda_MS_bar](double x) {
	double Mh = (1.0/Minv);
	double mm= MRS_to_mm(x*mc_HQET.distr[ijack]);
	double a= Get_4l_alpha_s(mm, 4.0, Lambda_MS_bar);
	double log_mug= pow(a,9.0/25)*( 1 + 0.67235*a + 1.284*pow(a,2));
	return (x/Mh)*mc_HQET.distr[ijack]*( R_mh.distr[ijack] + A_mh.distr[ijack]*(1.0/x) + B_mh.distr[ijack]*pow(1.0/x,2) + C_mh.distr[ijack]*log_mug*pow(1.0/x,2)) -1 ; 
      };
      double m_max=2;
      double m_min=0.6;
      if(pp_fB.M > 10) { m_max= 1.2; m_min=0.9;}
      pp_fB.Mpole= mc_HQET.distr[ijack]*R_brent( Fmh, m_min*(1.0/Minv)/mc_HQET.ave(),  m_max*(1.0/Minv)/mc_HQET.ave() ) ;
      C_HQET_M.distr.push_back(C_HQET(pp_fB.M));
      F_fB_M.distr.push_back( bf_fB.ansatz( Bt_fit_fB.par[ijack], pp_fB)/sqrt(pp_fB.M));
    }
    C_HQET_fit.distr_list.push_back( C_HQET_M);
    F_fB_fit.distr_list.push_back(F_fB_M);
  }
  //print
  string out_path_fBs="../data/ph_emission/"+ph_type+"/Bs_extr/fBs.fit_func";

  //C_HEQ on simulated masses
  Vfloat HQET_on_simulated_m, pole_mass_on_simulated_m;
  for(int i=0;i<masses.size();i++) {  pole_mass_on_simulated_m.push_back( MS_bar_to_MRS_mass( 3.0 , 4.0, Lambda_MS_bar , mc_MS_distr.ave()/ratio_mh[i],  mc_MS_distr.ave() )); HQET_on_simulated_m.push_back( C_HQET(pole_mass_on_simulated_m[i]));}  

  Print_To_File({},{ Inv_Masses_to_print, F_fB_fit.ave(), F_fB_fit.err(), C_HQET_fit.ave()} , out_path_fBs , "", "ch2/dof: "+to_string_with_precision(ch2_red_fB,5));
  Print_To_File({}, {ratio_mh, (1.0/masses).ave(), (1.0/masses).err(), fps.ave(), fps.err(), HQET_on_simulated_m, pole_mass_on_simulated_m}, "../data/ph_emission/"+ph_type+"/Bs_extr/fBs.dat", "", "");

  }


  //###################################################################################################################################
  //###################################################################################################################################
  //###################################################################################################################################
  //###################################################################################################################################
  //##############################################                               ######################################################
  //##############################################                               ######################################################
  //##############################################                               ######################################################
  //##############################################                               ######################################################
  //##############################################       PRETTY DIAGRAM          ######################################################
  //##############################################                               ######################################################
  //##############################################                               ######################################################
  //##############################################                               ######################################################
  //##############################################                               ######################################################
  //###################################################################################################################################
  //###################################################################################################################################
  //###################################################################################################################################
  //###################################################################################################################################

  int nxg_to_print_pretty=199;
  distr_t_list pretty_b_MB(1, nxg_to_print_pretty, NJ);
  distr_t_list pretty_b_MB_red(1, nxg_to_print_pretty, NJ);
  distr_t_list pretty_s_MB_RE(1, nxg_to_print_pretty, NJ);
  distr_t_list pretty_s_MB_IM(1, nxg_to_print_pretty, NJ);

  if(!Skip_virtual_diagram) {

  //mass extrapolation

  
  class ipar_pretty_b {
  public:
    ipar_pretty_b() : FF(0.0), FF_err(0.0) {}
    double FF, FF_err, M, xg;
  };
  
  class fpar_pretty_b {
  public:
    fpar_pretty_b() {}
    fpar_pretty_b(const Vfloat &par) {
      if((signed)par.size() != 3) crash("In class fpar_pretty_b  class constructor Vfloat par has size != 3");
      A=par[0];
      B=par[1];
      C=par[2];
    }
    double A,B,C;
  };


  vector<distr_t_list> bar_FT_b;

 
  for(int i=0; i<(signed)TFF_virtual_ret_list.size();i++) bar_FT_b.push_back( TFF_virtual_ret_list[i].Get_FF(0));

  int Njacks_pretty= bar_FT_b[0][0].distr.size();

  cout<<"Njacks pretty: "<<Njacks_pretty<<endl;

  if(Njacks_pretty != NJ) crash("Njacks pretty and NJ are not equal");

   
  //init bootstrap fit
  bootstrap_fit<fpar_pretty_b,ipar_pretty_b> bf_pretty_b(Njacks_pretty);
  bf_pretty_b.set_warmup_lev(1); //sets warmup
  bf_pretty_b.Set_number_of_measurements(12);
  bf_pretty_b.Set_verbosity(1);
  
  
  bf_pretty_b.Add_par("A", 0.2, 0.01);
  bf_pretty_b.Add_par("B", 0.1 , 0.01);
  bf_pretty_b.Add_par("C", -0.7 , 0.1);
 
  //fit on mean values to get ch2
  bootstrap_fit<fpar_pretty_b,ipar_pretty_b> bf_pretty_b_ch2(1);
  bf_pretty_b_ch2.set_warmup_lev(1); //sets warmup
  bf_pretty_b_ch2.Set_number_of_measurements(12);
  bf_pretty_b_ch2.Set_verbosity(1);
  bf_pretty_b_ch2.Add_par("A", 0.2, 0.01);
  bf_pretty_b_ch2.Add_par("B", 0.1, 0.01); 
  bf_pretty_b_ch2.Add_par("C", -0.7, 0.1);



 

  //ansatz
  bf_pretty_b.ansatz=  [ ](const fpar_pretty_b &p, const ipar_pretty_b &ip) {
    return (p.A + p.B*ip.xg)/(ip.M*( 1 + 0.5*ip.xg + p.C/ip.M ))  ;
  };
  
  bf_pretty_b.measurement=  [ ](const fpar_pretty_b &p, const ipar_pretty_b &ip) {
    return ip.FF;
  };
  
  bf_pretty_b.error=  [ ](const fpar_pretty_b &p, const ipar_pretty_b &ip) {
    return ip.FF_err;
  };
	
  bf_pretty_b_ch2.ansatz= bf_pretty_b.ansatz;
  bf_pretty_b_ch2.measurement = bf_pretty_b.measurement;
  bf_pretty_b_ch2.error = bf_pretty_b.error;

  

    

  //insert covariance matrix
  Eigen::MatrixXd Cov_Matrix_pretty_b(12,12);
  Eigen::MatrixXd Corr_Matrix_pretty_b(12,12);
  for(int i=0;i<12;i++)
    for(int j=0;j<12;j++) {
      Cov_Matrix_pretty_b(i,j) = bar_FT_b[(i+3)%3].distr_list[i/3]%bar_FT_b[(j+3)%3].distr_list[j/3];
    }

  for(int i=0;i<12;i++) {
    for(int j=0;j<12;j++) {
      Corr_Matrix_pretty_b(i,j) = Cov_Matrix_pretty_b(i,j)/(sqrt(Cov_Matrix_pretty_b(i,i)*Cov_Matrix_pretty_b(j,j)));
    }
  }

  //print on screen
  cout<<"Correlation matrix pretty b"<<endl;
  cout<<Corr_Matrix_pretty_b<<endl;
  
  //add cov matrix to bootstrap fit
  //bf_pretty_b.Add_covariance_matrix(Cov_Matrix_pretty_b);
  //bf_pretty_b_ch2.Add_covariance_matrix(Cov_Matrix_pretty_b);

  



  //start fitting
  //fill the data
  vector<vector<ipar_pretty_b>> data_pretty_b(Njacks_pretty);
  vector<vector<ipar_pretty_b>> data_pretty_b_ch2(1);
  //allocate space for output result
  boot_fit_data<fpar_pretty_b> Bt_fit_pretty_b;
  boot_fit_data<fpar_pretty_b> Bt_fit_pretty_b_ch2;
  boot_fit_data<fpar_pretty_b> Bt_fit_pretty_b_red;
  boot_fit_data<fpar_pretty_b> Bt_fit_pretty_b_red_ch2;
  for(auto &data_iboot: data_pretty_b) data_iboot.resize(12);
  for(auto &data_iboot: data_pretty_b_ch2) data_iboot.resize(12);
  for(int ijack=0;ijack<Njacks_pretty;ijack++) {
    for(int im=0;im<12;im++) {
      data_pretty_b[ijack][im].FF = bar_FT_b[(im+3)%3][im/3].distr[ijack];
      data_pretty_b[ijack][im].FF_err= bar_FT_b[(im+3)%3].err(im/3);
      data_pretty_b[ijack][im].M= masses.distr_list[(((im+3)%3) > 1)?3:((im+3)%3)].distr[ijack];
      data_pretty_b[ijack][im].xg= Bs_xg_t_list[im/3];
      if(ijack==0) {
	data_pretty_b_ch2[ijack][im].FF =  bar_FT_b[(im+3)%3][im/3].ave();
	data_pretty_b_ch2[ijack][im].FF_err=  bar_FT_b[(im+3)%3][im/3].err();
	data_pretty_b_ch2[ijack][im].M = masses.ave((((im+3)%3) > 1)?3:((im+3)%3));
	data_pretty_b_ch2[ijack][im].xg= Bs_xg_t_list[im/3];
      }
    }
  }
  //append
  bf_pretty_b.Append_to_input_par(data_pretty_b);
  bf_pretty_b_ch2.Append_to_input_par(data_pretty_b_ch2);
  
  //fit
  cout<<"Fitting Bs-meson decay constant"<<endl;
  Bt_fit_pretty_b= bf_pretty_b.Perform_bootstrap_fit();
  Bt_fit_pretty_b_ch2= bf_pretty_b_ch2.Perform_bootstrap_fit();
  bf_pretty_b.Fix_par("B",0.0);
  bf_pretty_b_ch2.Fix_par("B",0.0);
  Bt_fit_pretty_b_red= bf_pretty_b.Perform_bootstrap_fit();
  Bt_fit_pretty_b_red_ch2= bf_pretty_b_ch2.Perform_bootstrap_fit();


  
  Bt_fit_pretty_b_ch2= bf_pretty_b_ch2.Perform_bootstrap_fit();
  double ch2_red_pretty_b= Bt_fit_pretty_b_ch2.get_ch2_ave()/( 12 -bf_pretty_b.Get_number_of_fit_pars());

  cout<<"number of fit pars: "<<bf_pretty_b.Get_number_of_fit_pars()<<endl;
  cout<<"reduced ch2 of pretty b: "<<ch2_red_pretty_b<<endl;
  
  //retrieve params
  distr_t A_pb(1), B_pb(1), C_pb(1);
  distr_t A_pb_red(1), C_pb_red(1);
  for(int ijack=0;ijack<NJ;ijack++) {
    A_pb.distr.push_back( Bt_fit_pretty_b.par[ijack].A);
    B_pb.distr.push_back( Bt_fit_pretty_b.par[ijack].B);
    C_pb.distr.push_back( Bt_fit_pretty_b.par[ijack].C);


    A_pb_red.distr.push_back( Bt_fit_pretty_b_red.par[ijack].A);
    C_pb_red.distr.push_back( Bt_fit_pretty_b_red.par[ijack].C);
  }

  cout<<"Effective pole position: "<<2*5.367 + C_pb.ave()<<" +- "<<C_pb.err()<<" [GeV] "<<endl;
  cout<<"Lambda_T : "<<C_pb.ave()<<" +- "<<C_pb.err()<<" [GeV] "<<endl;
  cout<<"Lambda_T red: "<<C_pb_red.ave()<<" +- "<<C_pb_red.err()<<" [GeV] "<<endl;
  

  Vfloat xg_to_print_pretty;
 
  for(int n=0;n<nxg_to_print_pretty;n++) xg_to_print_pretty.push_back( 0.0025*(n+1));

  //print fit function
  for(int im=0;im<3;im++) {
    distr_t_list pretty_b_to_print(true, nxg_to_print_pretty, Njacks_pretty);
    distr_t_list pretty_b_to_print_red(true, nxg_to_print_pretty, Njacks_pretty);
    for(int n=0;n<nxg_to_print_pretty;n++) {
      ipar_pretty_b ip_pretty_b;
      ip_pretty_b.xg= xg_to_print_pretty[n];
      for(int ijack=0;ijack<Njacks_pretty;ijack++) {
	ip_pretty_b.M= masses.distr_list[ (im==2)?3:im ].distr[ijack];
	pretty_b_to_print.distr_list[n].distr[ijack] = bf_pretty_b.ansatz( Bt_fit_pretty_b.par[ijack], ip_pretty_b);
	pretty_b_to_print_red.distr_list[n].distr[ijack] = bf_pretty_b.ansatz( Bt_fit_pretty_b_red.par[ijack], ip_pretty_b);
      }
    }
    Print_To_File({}, {xg_to_print_pretty, pretty_b_to_print.ave(), pretty_b_to_print.err(), pretty_b_to_print_red.ave(), pretty_b_to_print_red.err()},  "../data/ph_emission/"+ph_type+"/Bs_extr/pretty_b_im_"+to_string(im)+".fit_func", "" , "");
    Print_To_File({}, {Bs_xg_t_list, bar_FT_b[im].ave(), bar_FT_b[im].err()}, "../data/ph_emission/"+ph_type+"/Bs_extr/pretty_b_im_"+to_string(im)+".dat", "" , "");
  }

  //determine pretty b at MBs mass

  
  for(int n=0;n<nxg_to_print_pretty;n++) {
      ipar_pretty_b ip_pretty_b;
      ip_pretty_b.xg= xg_to_print_pretty[n];
      ip_pretty_b.M= MBs;
      for(int ijack=0;ijack<Njacks_pretty;ijack++) {
	pretty_b_MB.distr_list[n].distr[ijack] =  bf_pretty_b.ansatz( Bt_fit_pretty_b.par[ijack], ip_pretty_b);
	pretty_b_MB_red.distr_list[n].distr[ijack] =  bf_pretty_b.ansatz( Bt_fit_pretty_b_red.par[ijack], ip_pretty_b);
      }
  }
  Print_To_File({}, {xg_to_print_pretty, pretty_b_MB.ave(), pretty_b_MB.err(), pretty_b_MB_red.ave(), pretty_b_MB_red.err()},  "../data/ph_emission/"+ph_type+"/Bs_extr/pretty_b_Bs.dat", "" , "");

  vector<distr_t_list> bar_FT_s_RE_xg(4);
  vector<distr_t_list> bar_FT_s_IM_xg(4);

  vector<distr_t_list> bar_FT_s_RE;
  vector<distr_t_list> bar_FT_s_IM;

  for(int i=0; i<(signed)TFF_virtual_ret_list.size();i++) bar_FT_s_RE.push_back( TFF_virtual_ret_list[i].Get_FF(1));
  for(int i=0; i<(signed)TFF_virtual_ret_list.size();i++) bar_FT_s_IM.push_back( TFF_virtual_ret_list[i].Get_FF(2));

  for(int ixg=0;ixg<(signed)Bs_xg_t_list.size();ixg++) {
    for(int i=0; i<(signed)TFF_virtual_ret_list.size();i++) {
      bar_FT_s_RE_xg[ixg].distr_list.push_back( TFF_virtual_ret_list[i].Get_FF(1).distr_list[ixg]);
      bar_FT_s_IM_xg[ixg].distr_list.push_back( TFF_virtual_ret_list[i].Get_FF(2).distr_list[ixg]);
    }
  }

  vector<double> M_07({ masses.ave(0), masses.ave(1), masses.ave(3)});

   for(int ixg=0;ixg<(signed)Bs_xg_t_list.size();ixg++) {
     Print_To_File({}, {M_07, bar_FT_s_RE_xg[ixg].ave(), bar_FT_s_RE_xg[ixg].err()}, "../data/ph_emission/"+ph_type+"/Bs_extr/pretty_s_RE_ixg_"+to_string(ixg)+".dat", "" , "");
     Print_To_File({}, {M_07, bar_FT_s_IM_xg[ixg].ave(), bar_FT_s_IM_xg[ixg].err()}, "../data/ph_emission/"+ph_type+"/Bs_extr/pretty_s_IM_ixg_"+to_string(ixg)+".dat", "" , "");
   }
  
  /*
 
  //pretty s 


   class ipar_pretty_s {
  public:
    ipar_pretty_s() : FF(0.0), FF_err(0.0) {}
     double FF, FF_err, M, xg, y;
     bool Is_RE;
  };
  
  class fpar_pretty_s {
  public:
    fpar_pretty_s() {}
    fpar_pretty_s(const Vfloat &par) {
      if((signed)par.size() != 5) crash("In class fpar_pretty_s  class constructor Vfloat par has size != 3");
      A=par[0];
      B=par[1];
      B2=par[2];
      C=par[3];
      Eth=par[4];
    }
    double A,B,B2,C,Eth;
  };


  

  vector<distr_t_list> bar_FT_s_RE_xg(4);
  vector<distr_t_list> bar_FT_s_IM_xg(4);

  for(int i=0; i<(signed)TFF_virtual_ret_list.size();i++) bar_FT_s_RE.push_back( TFF_virtual_ret_list[i].Get_FF(1));
  for(int i=0; i<(signed)TFF_virtual_ret_list.size();i++) bar_FT_s_IM.push_back( TFF_virtual_ret_list[i].Get_FF(2));

  for(int ixg=0;ixg<(signed)Bs_xg_t_list.size();ixg++) {
    for(int i=0; i<(signed)TFF_virtual_ret_list.size();i++) {
      bar_FT_s_RE_xg[ixg].distr_list.push_back( TFF_virtual_ret_list[i].Get_FF(1).distr_list[ixg]);
      bar_FT_s_IM_xg[ixg].distr_list.push_back( TFF_virtual_ret_list[i].Get_FF(2).distr_list[ixg]);
    }
  }

 

  Vfloat y_eff= TFF_virtual_ret_list[0].y_eff;

  //init bootstrap fit
  bootstrap_fit<fpar_pretty_s,ipar_pretty_s> bf_pretty_s(Njacks_pretty);
  bf_pretty_s.set_warmup_lev(1); //sets warmup
  bf_pretty_s.Set_number_of_measurements(24);
  bf_pretty_s.Set_verbosity(1);
  
  
  bf_pretty_s.Add_par("A", -0.05, 0.003);
  bf_pretty_s.Add_par("B", -0.05 , 0.001);
  bf_pretty_s.Add_par("B2", 1.5 , 0.01);
  bf_pretty_s.Add_par("C", 1.5 , 0.01);
  bf_pretty_s.Add_par("Eth", 1.0, 0.1);
 
  //fit on mean values to get ch2
  bootstrap_fit<fpar_pretty_s,ipar_pretty_s> bf_pretty_s_ch2(1);
  bf_pretty_s_ch2.set_warmup_lev(1); //sets warmup
  bf_pretty_s_ch2.Set_number_of_measurements(24);
  bf_pretty_s_ch2.Set_verbosity(1);
  bf_pretty_s_ch2.Add_par("A", -0.05, 0.003);
  bf_pretty_s_ch2.Add_par("B", -0.05, 0.001); 
  bf_pretty_s_ch2.Add_par("B2", 1.5, 0.01);
  bf_pretty_s_ch2.Add_par("C", 1.5, 0.01);
  bf_pretty_s_ch2.Add_par("Eth", 1.0, 0.1);




    
  bf_pretty_s.ansatz =  [ ](const fpar_pretty_s &p, const ipar_pretty_s &ip) {

    

    double E_tilde= (ip.y*ip.M + (1-ip.y)*MBs)*(1 -0.5*ip.xg);
    double En= sqrt( pow(p.C/ip.M,2) + pow(0.5*ip.xg,2));

    double eff_th= sqrt( pow(p.Eth,2) + pow(0.5*ip.xg*ip.M,2));

    double D= eff_th*(p.A + p.B*ip.xg + p.B2/ip.M)/(ip.M*ip.M*En*(En - E_tilde/ip.M));

    
    if(!ip.Is_RE) {
      return D/E_tilde;
    }
    else {
      return (p.A + p.B*ip.xg + p.B2/ip.M)/(ip.M*ip.M*En*(En - E_tilde/ip.M))    ;
    }
        
  };
  
  bf_pretty_s.measurement=  [ ](const fpar_pretty_s &p, const ipar_pretty_s &ip) {
    return ip.FF;
  };
  
  bf_pretty_s.error=  [ ](const fpar_pretty_s &p, const ipar_pretty_s &ip) {
    return ip.FF_err;
  };
	
  bf_pretty_s_ch2.ansatz= bf_pretty_s.ansatz;
  bf_pretty_s_ch2.measurement = bf_pretty_s.measurement;
  bf_pretty_s_ch2.error = bf_pretty_s.error;


  
  //insert covariance matrix
  Eigen::MatrixXd Cov_Matrix_pretty_s(24,24);
  Eigen::MatrixXd Corr_Matrix_pretty_s(24,24);
  for(int i=0;i<12;i++)
    for(int j=0;j<12;j++) {
      Cov_Matrix_pretty_s(i,j) = bar_FT_s_RE[(i+3)%3].distr_list[i/3]%bar_FT_s_RE[(j+3)%3].distr_list[j/3];
      
    }

  for(int i=0;i<24;i++) {
    for(int j=0;j<24;j++) {
      Corr_Matrix_pretty_s(i,j) = Cov_Matrix_pretty_s(i,j)/(sqrt(Cov_Matrix_pretty_s(i,i)*Cov_Matrix_pretty_s(j,j)));
    }
  }

  //print on screen
  cout<<"Correlation matrix pretty s"<<endl;
  cout<<Corr_Matrix_pretty_s<<endl;
  
  //add cov matrix to bootstrap fit
  //bf_pretty_s.Add_covariance_matrix(Cov_Matrix_pretty_s);
  //bf_pretty_s_ch2.Add_covariance_matrix(Cov_Matrix_pretty_s);



  //start fitting
  //fill the data
  vector<vector<ipar_pretty_s>> data_pretty_s(Njacks_pretty);
  vector<vector<ipar_pretty_s>> data_pretty_s_ch2(1);
  //allocate space for output result
  boot_fit_data<fpar_pretty_s> Bt_fit_pretty_s;
  boot_fit_data<fpar_pretty_s> Bt_fit_pretty_s_ch2;
  for(auto &data_iboot: data_pretty_s) data_iboot.resize(12);
  for(auto &data_iboot: data_pretty_s_ch2) data_iboot.resize(12);
  for(int ijack=0;ijack<Njacks_pretty;ijack++) {
    for(int im=0;im<12;im++) {
      data_pretty_s[ijack][im].FF = bar_FT_s_RE[(im+3)%3][im/3].distr[ijack];
      data_pretty_s[ijack][im].FF_err= bar_FT_s_RE[(im+3)%3].err(im/3);
      data_pretty_s[ijack][im].M= masses.distr_list[(((im+3)%3) > 1)?3:((im+3)%3)].distr[ijack];
      data_pretty_s[ijack][im].xg= Bs_xg_t_list[im/3];
      data_pretty_s[ijack][im].y= y_eff[im/3];
      data_pretty_s[ijack][im].Is_RE= true;

      data_pretty_s[ijack][im+12].FF = bar_FT_s_IM[(im+3)%3][im/3].distr[ijack];
      data_pretty_s[ijack][im+12].FF_err= bar_FT_s_IM[(im+3)%3].err(im/3);
      data_pretty_s[ijack][im+12].M= masses.distr_list[(((im+3)%3) > 1)?3:((im+3)%3)].distr[ijack];
      data_pretty_s[ijack][im+12].xg= Bs_xg_t_list[im/3];
      data_pretty_s[ijack][im+12].y= y_eff[im/3];
      data_pretty_s[ijack][im+12].Is_RE= false;

      
      if(ijack==0) {
	data_pretty_s_ch2[ijack][im].FF =  bar_FT_s_RE[(im+3)%3][im/3].ave();
	data_pretty_s_ch2[ijack][im].FF_err=  bar_FT_s_RE[(im+3)%3][im/3].err();
	data_pretty_s_ch2[ijack][im].M = masses.ave((((im+3)%3) > 1)?3:((im+3)%3));
	data_pretty_s_ch2[ijack][im].xg= Bs_xg_t_list[im/3];
	data_pretty_s_ch2[ijack][im].y = y_eff[im/3];
	data_pretty_s_ch2[ijack][im].Is_RE= true;


	data_pretty_s_ch2[ijack][im+12].FF =  bar_FT_s_IM[(im+3)%3][im/3].ave();
	data_pretty_s_ch2[ijack][im+12].FF_err=  bar_FT_s_IM[(im+3)%3][im/3].err();
	data_pretty_s_ch2[ijack][im+12].M = masses.ave((((im+3)%3) > 1)?3:((im+3)%3));
	data_pretty_s_ch2[ijack][im+12].xg= Bs_xg_t_list[im/3];
	data_pretty_s_ch2[ijack][im+12].y = y_eff[im/3];
	data_pretty_s_ch2[ijack][im+12].Is_RE= false;
      }
    }
  }
  
  //append
  bf_pretty_s.Append_to_input_par(data_pretty_s);
  bf_pretty_s_ch2.Append_to_input_par(data_pretty_s_ch2);
  //fit
  cout<<"Fitting Bs-meson decay constant"<<endl;
  Bt_fit_pretty_s= bf_pretty_s.Perform_bootstrap_fit();
  
  
  Bt_fit_pretty_s_ch2= bf_pretty_s_ch2.Perform_bootstrap_fit();
  double ch2_red_pretty_s= Bt_fit_pretty_s_ch2.get_ch2_ave()/( 24 -bf_pretty_s.Get_number_of_fit_pars());

  cout<<"reduced ch2 of pretty s: "<<ch2_red_pretty_s<<endl;
  
  //retrieve params
  distr_t A_ps_RE(1), B_ps_RE(1), C_ps_RE(1);
  for(int ijack=0;ijack<Njacks_pretty;ijack++) {
    A_ps_RE.distr.push_back( Bt_fit_pretty_s.par[ijack].A);
    B_ps_RE.distr.push_back( Bt_fit_pretty_s.par[ijack].B);
    C_ps_RE.distr.push_back( Bt_fit_pretty_s.par[ijack].C);
  }

  
  Vfloat Inv_Masses_to_print_pretty;
  for(int i=0;i<200;i++) { Inv_Masses_to_print_pretty.push_back(   0.30*2*(i+1)/200  ) ;}
  //print fit function
  for(int ixg=0;ixg<(signed)Bs_xg_t_list.size();ixg++) {
    distr_t_list pretty_s_to_print_RE(true, Inv_Masses_to_print_pretty.size(), Njacks_pretty);
    distr_t_list pretty_s_to_print_IM(true, Inv_Masses_to_print_pretty.size(), Njacks_pretty);
    for(int n=0;n<200;n++) {
      ipar_pretty_s ip_pretty_s;
      ip_pretty_s.xg= Bs_xg_t_list[ixg];
      ip_pretty_s.y = y_eff[ixg];
      ip_pretty_s.M= 1.0/Inv_Masses_to_print_pretty[n];
      for(int ijack=0;ijack<Njacks_pretty;ijack++) {
	ip_pretty_s.Is_RE=true;
	pretty_s_to_print_RE.distr_list[n].distr[ijack] = bf_pretty_s.ansatz( Bt_fit_pretty_s.par[ijack], ip_pretty_s);
	ip_pretty_s.Is_RE=false;
	pretty_s_to_print_IM.distr_list[n].distr[ijack] = bf_pretty_s.ansatz( Bt_fit_pretty_s.par[ijack], ip_pretty_s);
      }
    }
    Print_To_File({}, {Inv_Masses_to_print_pretty, pretty_s_to_print_RE.ave(), pretty_s_to_print_RE.err()},  "../data/ph_emission/"+ph_type+"/Bs_extr/pretty_s_RE_ixg_"+to_string(ixg)+".fit_func", "" , "");
    Print_To_File({}, {Inv_Masses_to_print_pretty, pretty_s_to_print_IM.ave(), pretty_s_to_print_IM.err()},  "../data/ph_emission/"+ph_type+"/Bs_extr/pretty_s_IM_ixg_"+to_string(ixg)+".fit_func", "" , "");

    vector<double> M_07({ masses.ave(0), masses.ave(1), masses.ave(3)});
    
    Print_To_File({}, {M_07, bar_FT_s_RE_xg[ixg].ave(), bar_FT_s_RE_xg[ixg].err()}, "../data/ph_emission/"+ph_type+"/Bs_extr/pretty_s_RE_ixg_"+to_string(ixg)+".dat", "" , "");
    Print_To_File({}, {M_07, bar_FT_s_IM_xg[ixg].ave(), bar_FT_s_IM_xg[ixg].err()}, "../data/ph_emission/"+ph_type+"/Bs_extr/pretty_s_IM_ixg_"+to_string(ixg)+".dat", "" , "");
  }

  //determine pretty b at MBs mass

 
  for(int n=0;n<nxg_to_print_pretty;n++) {
      ipar_pretty_s ip_pretty_s;
      ip_pretty_s.xg= xg_to_print_pretty[n];
      ip_pretty_s.M= MBs;
      ip_pretty_s.y= 1.0;
      for(int ijack=0;ijack<Njacks_pretty;ijack++) {
	pretty_s_MB_RE.distr_list[n].distr[ijack] = bf_pretty_s.ansatz( Bt_fit_pretty_s.par[ijack], ip_pretty_s);
      }
  }
  Print_To_File({}, {xg_to_print_pretty, pretty_s_MB_RE.ave(), pretty_s_MB_RE.err()},  "../data/ph_emission/"+ph_type+"/Bs_extr/pretty_s_RE_Bs.dat", "" , "");
  Print_To_File({}, {xg_to_print_pretty, pretty_s_MB_IM.ave(), pretty_s_MB_IM.err()},  "../data/ph_emission/"+ph_type+"/Bs_extr/pretty_s_IM_Bs.dat", "" , "");



  */

 

  }
	
  
  //###################################################################################################################################
  //###################################################################################################################################
  //###################################################################################################################################
  //###################################################################################################################################
  //###################################################################################################################################
  //###################################################################################################################################
  //###################################################################################################################################
  //#########################################                                    ######################################################
  //#########################################                                    ######################################################
  //#########################################                                    ######################################################
  //#########################################     UNCONTRAINED EXTR.             ######################################################
  //#########################################                                    ######################################################
  //#########################################                                    ######################################################
  //#########################################                                    ######################################################
  //###################################################################################################################################
  //###################################################################################################################################
  //###################################################################################################################################
  //###################################################################################################################################
  //###################################################################################################################################
  //###################################################################################################################################
  //###################################################################################################################################
  
  //print continuum extrapolated results for each contribution as a function of the mass
  vector<string> contribs({"FA", "FA_u", "FA_d", "FV",  "FV_u", "FV_d",  "FA_T", "FA_T_u", "FA_T_d", "FV_T", "FV_T_u", "FV_T_d", "FB", "FB_u", "FB_d", "FT", "FT_u", "FT_d"});
  vector<string> contribs_boot({"FA", "FA_u", "FA_d", "FV",  "FV_u", "FV_d",  "FA_T", "FA_T_u", "FA_T_d", "FV_T", "FV_T_u", "FV_T_d", "FB", "FB_u", "FB_d", "FT", "FT_u", "FT_d"});
  vector<int> FF_id({0,0,0,1,1,1,2,2,2,3,3,3,4,4,4,5,5,5});
  vector<double> FF_sign({1,1,-1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1});
  vector<vector<distr_t_list>> FF_M_list(contribs.size());
  vector<vector<distr_t_list>> FF_M_list_boot(contribs.size());
  vector<vector<distr_t_list>> FF_M_mass_list(contribs.size());
  for(auto &FF_c: FF_M_list)
    for(int im=0;im<Nmasses;im++)  FF_c.emplace_back(1);

  for(auto &FF_c: FF_M_list_boot)
    for(int im=0;im<Nmasses;im++)  FF_c.emplace_back(1);

  for(auto &FF_c: FF_M_mass_list)
    for(int ixg=1;ixg<num_xg_list[0];ixg++) FF_c.emplace_back(1);
 
  
  for(int c=0; c<(signed)contribs.size(); c++) {

    distr_t_list FF_Bs_xg;
    
   
    vector<distr_t_list> FF_M;
    vector<distr_t_list> FF_M_boot;
    for(int ixg=1; ixg<num_xg_list[0]; ixg++) FF_M.emplace_back(1);
     for(int ixg=1; ixg<num_xg_list[0]; ixg++) FF_M_boot.emplace_back(0);
  
    for(int m=0; m<Nmasses; m++) {
      for(int ixg=1; ixg<num_xg_list[0];ixg++) {

	//jackknife
	vector<distr_t> FF_AIC;
	vector<double> Ch2_AIC;
	vector<int> Ndof_AIC;
	vector<int> Nmeas_AIC;
	int k=0;
	while(k*Nmasses+m < FF_SIZE) {
	  int r= k*Nmasses+m;
	  FF_AIC.push_back( FF_ret_list[r].Get_FF(c).distr_list[ixg-1]);
	  Ch2_AIC.push_back( FF_ret_list[r].Get_ch2(c)[ixg-1]*FF_ret_list[r].Ndof);
	  Ndof_AIC.push_back( FF_ret_list[r].Ndof);
	  Nmeas_AIC.push_back( FF_ret_list[r].Nmeas );
	  k++;
	}

	distr_t RES=AIC_lin( FF_AIC, Ch2_AIC, Ndof_AIC, Nmeas_AIC,0);
	FF_M[ixg-1].distr_list.push_back(RES);
	FF_M_list[c][m].distr_list.push_back((c==2)?(-1.0*RES):RES);
	FF_M_mass_list[c][ixg-1].distr_list.push_back( (c==2)?(-1.0*RES):RES);



	//boot
	vector<distr_t> FF_AIC_boot;
	vector<double> Ch2_AIC_boot;
	vector<int> Ndof_AIC_boot;
	vector<int> Nmeas_AIC_boot;
	k=0;
	while(k*Nmasses+m < FF_SIZE) {
	  int r= k*Nmasses+m;
	  FF_AIC_boot.push_back( FF_ret_list_boot[r].Get_FF(c).distr_list[ixg-1]);
	  Ch2_AIC_boot.push_back( FF_ret_list_boot[r].Get_ch2(c)[ixg-1]*FF_ret_list_boot[r].Ndof);
	  Ndof_AIC_boot.push_back( FF_ret_list_boot[r].Ndof);
	  Nmeas_AIC_boot.push_back( FF_ret_list_boot[r].Nmeas );
	  k++;
	}

	distr_t RES_boot=AIC_lin( FF_AIC_boot, Ch2_AIC_boot, Ndof_AIC_boot, Nmeas_AIC_boot,0);
	FF_M_boot[ixg-1].distr_list.push_back(RES_boot);
	FF_M_list_boot[c][m].distr_list.push_back((c==2)?(-1.0*RES_boot):RES_boot);
	
	
      }
    }
    
  
    distr_t_list B_coeff(1);
    
    boost::filesystem::create_directory("../data/ph_emission/"+ph_type+"/"+"Bs_extr");
    boost::filesystem::create_directory("../data/ph_emission/"+ph_type+"/"+"Bs_extr/comb_fits");
    for(int ixg=1;ixg<num_xg_list[0];ixg++) {
     
    

      //perform extrapolation in the heavy quark mass
      if(contribs[c] != "FA" && contribs[c] != "FV" && contribs[c] != "FV_T" && contribs[c] != "FA_T" && contribs[c] != "FB" && contribs[c] != "FT" ) {
	
	 
	class ipar_MEX {
	public:
	  ipar_MEX() : FF(0.0), FF_err(0.0) {}
	  double FF, FF_err, M, xg;
	};
  
	class fpar_MEX {
	public:
	  fpar_MEX() {}
	  fpar_MEX(const Vfloat &par) {
	    if((signed)par.size() != 3) crash("In class fpar_MEX  class constructor Vfloat par has size != 3");
	    A=par[0];
	    B=par[1];
	    C=par[2];
	  }
	  double A,B,C;
	};
	
	
	//init bootstrap fit
	bootstrap_fit<fpar_MEX,ipar_MEX> bf_MEX(NJ);
	bf_MEX.set_warmup_lev(2); //sets warmup
	bf_MEX.Set_number_of_measurements(Nmasses);
	bf_MEX.Set_verbosity(1);
	
	double guess_B= sqrt(2*0.24/Bs_xg_t_list[ixg-1]);
	if( (contribs[c].substr(0,2) == "FA") )  guess_B=  sqrt(2*0.5/Bs_xg_t_list[ixg-1]);
	if( (contribs[c].substr(0,2) == "FB")) guess_B = sqrt( 2*0.5/Bs_xg_t_list[ixg-1]);

	double G= (contribs[c].substr( contribs[c].length()-1, 1) == "u")?masses.ave(2):1.0;
	
	bf_MEX.Add_par("A", (G*FF_M[ixg-1]/fps).ave(3), (G*FF_M[ixg-1]/fps).ave(3)/10);
	bf_MEX.Add_par("B", guess_B , guess_B/10);
	bf_MEX.Add_par("C", 0.1, 0.01);
	//fit on mean values to get ch2
	bootstrap_fit<fpar_MEX,ipar_MEX> bf_MEX_ch2(1);
	bf_MEX_ch2.set_warmup_lev(1); //sets warmup
	bf_MEX_ch2.Set_number_of_measurements(Nmasses);
	bf_MEX_ch2.Set_verbosity(1);
	bf_MEX_ch2.Add_par("A", (G*FF_M[ixg-1]/fps).ave(3), (G*FF_M[ixg-1]/fps).ave(3)/10);
	bf_MEX_ch2.Add_par("B", guess_B, guess_B/10); 
	bf_MEX_ch2.Add_par("C", 0.1, 0.01); 
	bool C_fixed=true;
	if(contribs[c]=="FA_d") {  C_fixed=true;}
	if(contribs[c]=="FA_u") {  C_fixed=true;}
	if(C_fixed) {
	  bf_MEX.Fix_par("C", 0.0);
	  bf_MEX_ch2.Fix_par("C", 0.0);
	}

	bool axial_mode_fit_2=false;
	//ansatz
	
	bf_MEX.ansatz=  [&contribs, &c, &axial_mode_fit_2 ](const fpar_MEX &p, const ipar_MEX &ip) {

	  bool is_axial= (contribs[c].substr(0,2) == "FA") || (contribs[c].substr(0,2) == "FB") ;

	  double x= ip.M*ip.M;
	  
	  if(is_axial) x= ip.M;

	  double func= (p.A/( 1.0 + p.B*p.B/x ))*( 1 + p.C/ip.M );

	  if(is_axial && axial_mode_fit_2 ) func= p.A*( 1 + p.B/ip.M)*(1 + p.C/ip.M);

	  if(contribs[c].substr( contribs[c].length() -1, 1) == "u" ) func /= ip.M;

	  return func;
	  
	  
	};

	
	bf_MEX.measurement=  [ ](const fpar_MEX &p, const ipar_MEX &ip) {
	  
	  return ip.FF;
	};
	
	bf_MEX.error=  [ ](const fpar_MEX &p, const ipar_MEX &ip) {
	  
	  return ip.FF_err;
	};
	
	
	bf_MEX_ch2.ansatz= bf_MEX.ansatz;
	bf_MEX_ch2.measurement = bf_MEX.measurement;
	bf_MEX_ch2.error = bf_MEX.error;

	cout<<"Fitting contrib: "<<contribs[c]<<endl;



	//insert covariance matrix
	Eigen::MatrixXd Cov_Matrix_FF(Nmasses,Nmasses);
	Eigen::MatrixXd Corr_Matrix_FF(Nmasses,Nmasses);
	for(int i=0;i<Nmasses;i++)
	  for(int j=0;j<Nmasses;j++) {
	    Corr_Matrix_FF(i,j) = ((FF_M_boot[ixg-1]/fps_boot).distr_list[i])%((FF_M_boot[ixg-1]/fps_boot).distr_list[j])/( (FF_M_boot[ixg-1]/fps_boot).err(i)*(FF_M_boot[ixg-1]/fps_boot).err(j));
	    Cov_Matrix_FF(i,j) = Corr_Matrix_FF(i,j)*( FF_M[ixg-1]/fps).err(i)*(FF_M[ixg-1]/fps).err(j);
	  }

	//print on screen
	cout<<"Correlation matrix FF: "<<contribs[c]<<endl;
	cout<<Corr_Matrix_FF<<endl;
	
	//add cov matrix to bootstrap fit
	//bf_MEX.Add_covariance_matrix(Cov_Matrix_FF);
	//bf_MEX_ch2.Add_covariance_matrix(Cov_Matrix_FF);

	
      
	//start fitting
	//fill the data
	vector<vector<ipar_MEX>> data_MEX(NJ);
	vector<vector<ipar_MEX>> data_MEX_ch2(1);
	//allocate space for output result
	boot_fit_data<fpar_MEX> Bt_fit_MEX;
	boot_fit_data<fpar_MEX> Bt_fit_MEX_ch2;
	for(auto &data_iboot: data_MEX) data_iboot.resize(Nmasses);
	for(auto &data_iboot: data_MEX_ch2) data_iboot.resize(Nmasses);
	for(int ijack=0;ijack<NJ;ijack++) {
	  for(int im=0;im<Nmasses;im++) {
	    data_MEX[ijack][im].FF = (FF_M[ixg-1].distr_list[im]/fps.distr_list[im]).distr[ijack];
	    data_MEX[ijack][im].FF_err= (FF_M[ixg-1]/fps).err(im);
	    data_MEX[ijack][im].xg= Bs_xg_t_list[ixg-1];
	    data_MEX[ijack][im].M= masses.distr_list[im].distr[ijack];
	    if(ijack==0) {
	      data_MEX_ch2[ijack][im].FF = (FF_M[ixg-1]/fps).ave(im);
	      data_MEX_ch2[ijack][im].FF_err= (FF_M[ixg-1]/fps).err(im);
	      data_MEX_ch2[ijack][im].xg= Bs_xg_t_list[ixg-1];
	      data_MEX_ch2[ijack][im].M = masses.ave(im);	      
	    }
	  }
	}
	
	//append
	bf_MEX.Append_to_input_par(data_MEX);
	bf_MEX_ch2.Append_to_input_par(data_MEX_ch2);
	//fit
	cout<<"Fitting "<<contribs[c]<<" to extrapolate at the Bs-meson mass"<<endl;
	Bt_fit_MEX= bf_MEX.Perform_bootstrap_fit();
	Bt_fit_MEX_ch2= bf_MEX_ch2.Perform_bootstrap_fit();
	int Npars= C_fixed?2:3;
	double ch2_red_MEX= Bt_fit_MEX_ch2.get_ch2_ave()/( Nmasses -Npars);
	
	//retrieve params
	distr_t A(1), B(1), C(1);
	for(int ijack=0;ijack<NJ;ijack++) {
	  A.distr.push_back( Bt_fit_MEX.par[ijack].A);
	  B.distr.push_back( Bt_fit_MEX.par[ijack].B);
	  C.distr.push_back( Bt_fit_MEX.par[ijack].C);
	}

	B_coeff.distr_list.push_back( B*B*Bs_xg_t_list[ixg-1]/2.0);
	
	//print fit function
	distr_t_list F_MEX_fit(1);
	
	//### Determine FF at Bs meson mass
	ipar_MEX pp_Bs;
	pp_Bs.M= MBs;
	distr_t FF_Bs(1);
	for(int ijack=0;ijack<NJ;ijack++) FF_Bs.distr.push_back( bf_MEX.ansatz( Bt_fit_MEX.par[ijack], pp_Bs));
	FF_Bs_xg.distr_list.push_back( FF_Bs);
	//#################################


	 Vfloat Inv_Masses_to_print_MEX;
	 for(int i=0;i<1000;i++) { Inv_Masses_to_print_MEX.push_back(  0.30*2*i/999  ) ;}

	//######### PRINT FIT FUNCTION ##########
	//plot fit function
	for(auto &Minv: Inv_Masses_to_print_MEX) {
	  ipar_MEX pp_MEX;
	  pp_MEX.M=1.0/Minv;
	  distr_t F_MEX_M(1);
	  for(int ijack=0;ijack<NJ;ijack++) {
	    F_MEX_M.distr.push_back( bf_MEX.ansatz( Bt_fit_MEX.par[ijack], pp_MEX));
	  }
	  
	  F_MEX_fit.distr_list.push_back(F_MEX_M);
	}
	//print
	string out_path="../data/ph_emission/"+ph_type+"/Bs_extr/"+contribs[c]+"_ixg_"+to_string(ixg)+".fit_func";
	Print_To_File({},{ Inv_Masses_to_print_MEX, (F_MEX_fit*f_Bs).ave(), (F_MEX_fit*f_Bs).err(),  F_MEX_fit.ave(), F_MEX_fit.err()} , out_path , "", "ch2/dof: "+to_string_with_precision(ch2_red_MEX,5));
	
      }
    }
    //print FF_Bs
    if(contribs[c] != "FA" && contribs[c] != "FV" && contribs[c] != "FV_T" && contribs[c] != "FA_T" && contribs[c] != "FB" && contribs[c] != "FT" ) {
      Print_To_File({}, {Bs_xg_t_list, (FF_Bs_xg*f_Bs).ave(), (FF_Bs_xg*f_Bs).err(), FF_Bs_xg.ave(), FF_Bs_xg.err(),  B_coeff.ave(), B_coeff.err()}, "../data/ph_emission/"+ph_type+"/Bs_extr/"+contribs[c]+"_Bs.dat", "", "");
      FF_Bs_list[FF_id[c]] = FF_Bs_list[FF_id[c]] +  FF_sign[c]*FF_Bs_xg;
    }
  }
  Print_To_File({}, {Bs_xg_t_list,  (FF_Bs_list[0]*f_Bs).ave(), (FF_Bs_list[0]*f_Bs).err(),  FF_Bs_list[0].ave(), FF_Bs_list[0].err()}, "../data/ph_emission/"+ph_type+"/Bs_extr/FA_Bs.dat", "", "");
  Print_To_File({}, {Bs_xg_t_list,  (FF_Bs_list[1]*f_Bs).ave(), (FF_Bs_list[1]*f_Bs).err(),  FF_Bs_list[1].ave(), FF_Bs_list[1].err()}, "../data/ph_emission/"+ph_type+"/Bs_extr/FV_Bs.dat", "", "");
  Print_To_File({}, {Bs_xg_t_list,  (FF_Bs_list[2]*f_Bs).ave(), (FF_Bs_list[2]*f_Bs).err(),  FF_Bs_list[2].ave(), FF_Bs_list[2].err()}, "../data/ph_emission/"+ph_type+"/Bs_extr/FA_T_Bs.dat", "", "");
  Print_To_File({}, {Bs_xg_t_list,  (FF_Bs_list[3]*f_Bs).ave(), (FF_Bs_list[3]*f_Bs).err(),  FF_Bs_list[3].ave(), FF_Bs_list[3].err()}, "../data/ph_emission/"+ph_type+"/Bs_extr/FV_T_Bs.dat", "", "");

  Print_To_File({}, {Bs_xg_t_list,  (FF_Bs_list[4]*f_Bs).ave(), (FF_Bs_list[4]*f_Bs).err(),  FF_Bs_list[4].ave(), FF_Bs_list[4].err()}, "../data/ph_emission/"+ph_type+"/Bs_extr/FB_Bs.dat", "", "");
  Print_To_File({}, {Bs_xg_t_list,  (FF_Bs_list[5]*f_Bs).ave(), (FF_Bs_list[5]*f_Bs).err(),  FF_Bs_list[5].ave(), FF_Bs_list[5].err()}, "../data/ph_emission/"+ph_type+"/Bs_extr/FT_Bs.dat", "", "");


  //print form factors
  for(int c=0; c < (signed)contribs.size();c++) {

    for(int ixg=1;ixg<num_xg_list[0];ixg++) {

      vector<string> contribs({"FA", "FA_u", "FA_d", "FV",  "FV_u", "FV_d",  "FA_T", "FA_T_u", "FA_T_d", "FV_T", "FV_T_u", "FV_T_d", "FB", "FB_u", "FB_d", "FT", "FT_u", "FT_d"});

      if(c%3 == 0) FF_M_mass_list[c][ixg-1] = FF_M_mass_list[c+1][ixg-1] + FF_M_mass_list[c+2][ixg-1];
    
    Print_To_File({}, { ratio_mh, (1.0/masses).ave(), (1.0/masses).err(), fps.ave(), fps.err(), FF_M_mass_list[c][ixg-1].ave(), FF_M_mass_list[c][ixg-1].err(), (FF_M_mass_list[c][ixg-1]/fps).ave(), (FF_M_mass_list[c][ixg-1]/fps).err() }, "../data/ph_emission/"+ph_type+"/Bs_extr/"+contribs[c]+"_ixg_"+to_string(ixg)+".dat", "", "#mc/mh Mh  fh  FFh   FFh/fh");
    }

  }
    




  //###################################################################################################################################
  //###################################################################################################################################
  //###################################################################################################################################
  //###################################################################################################################################
  //###################################################################################################################################
  //###################################################################################################################################
  //###################################################################################################################################
  //#########################################                                    ######################################################
  //#########################################                                    ######################################################
  //#########################################                                    ######################################################
  //#########################################     CONSTRAINED EXTR.              ######################################################
  //#########################################                                    ######################################################
  //#########################################                                    ######################################################
  //#########################################                                    ######################################################
  //###################################################################################################################################
  //###################################################################################################################################
  //###################################################################################################################################
  //###################################################################################################################################
  //###################################################################################################################################
  //###################################################################################################################################
  //###################################################################################################################################


  distr_t_list FA_uniform(1), FV_uniform(1), FTA_uniform(1), FTV_uniform(1), FA_d_uniform(1);

  int num_xg=4;
  Vfloat xg_scan;
  int separate_ud = false;
  int nxgs_to_spline=199;
  double ZT_MB = evolutor_ZT_MS_bar(2.0,5.367,4);
  Vfloat ZT_Mh(Nmasses);
  double ZT_fit_Mb= evolutor_ZT_MS_bar(2.0,5,4); // 1.0;
  for(int im=0;im<Nmasses;im++)   ZT_Mh[im] = evolutor_ZT_MS_bar(2.0, 5,4);   //; //1.0 ; //evolutor_ZT_MS_bar(2.0, masses.ave(im),4);

  printV(ZT_Mh, "ZT_Mh",1);

   
  double evolutor_ZT = evolutor_ZT_MS_bar(2.0,5.0,4);
  for(int ixg=0;ixg<nxgs_to_spline;ixg++)  xg_scan.push_back( 0.0025*(ixg+1));

  Vfloat Inv_Masses_to_print_MEX;
  for(int i=0;i<200;i++) { Inv_Masses_to_print_MEX.push_back(   0.30*2*(i+1)/200  ) ;}
  Vfloat ZT_Mh_to_print(Inv_Masses_to_print_MEX.size());
  for(int i=0;i<(signed)ZT_Mh_to_print.size(); i++) ZT_Mh_to_print[i] = evolutor_ZT_MS_bar(2.0,5,4); //   1.0;  

  
  if(Perform_global_fit) {

  int Nfits=0;
  int multi_M2=0;

  vector<int> Is_M2_fit;
 

  //define the fit poll

  vector<map<string, int>> Fix_fit_pars_list;

  for(int dR= 0; dR <=1 ; dR ++) {
    for(int dC=0; dC<=1; dC ++) {
      for(int fix_den=0; fix_den<1 ; fix_den++) {
	for(int te1=0;te1<2;te1++) {
	  for(int l3=0; l3<=1; l3++) {
	    for(int V2=0;V2<=2;V2++) {
	      for(int A2=0;A2<=2;A2++) {
		for(int TV2=0;TV2<=2; TV2++) {
		  for(int TA2=0; TA2<=1; TA2++) {
		    
		    map<string,int> Fix_fit_pars;
		    
		    Fix_fit_pars.insert({"A1",dR});
		    Fix_fit_pars.insert({"A2", dC});
		    Fix_fit_pars.insert({"TVE1", te1});
		    Fix_fit_pars.insert({"M3", l3});
		    Fix_fit_pars.insert({"M4", 1});
		    Fix_fit_pars.insert({"SV1", (V2+2)%2});
		    Fix_fit_pars.insert({"SV2", V2/2 });
		    Fix_fit_pars.insert({"SA1", (A2+2)%2});
		    Fix_fit_pars.insert({"SA2", A2/2 });
		    Fix_fit_pars.insert({"STV1", (TV2+2)%2});
		    Fix_fit_pars.insert({"STV2", TV2/2 });
		    Fix_fit_pars.insert({"STA2", TA2 });
		    Fix_fit_pars.insert({"fix_T_den", fix_den});
		    
		    
		    //number of always fitted parameters:  A,DV,DA,M1,M2,TVM1,TAE1,LAMBDA
		    
		    int Npars= 8;
		  
		    Npars += dR + dC + l3 + (V2 != 0) + (A2 != 0) + (TV2 != 0) + (TA2 != 0) + (te1 != 0);
		  
		    if(Npars<=14 &&  (TV2 != 1 || TA2 != 1)) {
		      Nfits++;
		      if( V2+A2+TV2+TA2 > 0) { multi_M2++; Is_M2_fit.push_back(1);}
		      else Is_M2_fit.push_back(0);
		      Fix_fit_pars_list.push_back(Fix_fit_pars);
		    }
		    
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
  
  vector<distr_t> CA_list, CV_list, CA_T_list, R_list, RT_list, M3_list, M4_list;
  vector<distr_t> M1_list, M2_list, M1_T_list, M2_T_list, d_T_list;

  cout<<"TOTAL NUMBER OF FITS: "<<Nfits<<endl;
  
  vector<int> separate_ud_list(Nfits,0);

  vector<double> mult_fits;
  for(int ifit=0;ifit<Nfits;ifit++) {
    if(Is_M2_fit[ifit] == 1) mult_fits.push_back(1.0/(double)multi_M2) ;
    else mult_fits.push_back(1.0);
  }

  
  

  
 
  
  
  vector<vector<distr_t>> FV_extr(nxgs_to_spline), FA_extr(nxgs_to_spline), FTV_extr(nxgs_to_spline), FTA_extr(nxgs_to_spline), FV_d_extr(nxgs_to_spline), FA_d_extr(nxgs_to_spline), FTV_d_extr(nxgs_to_spline), FTA_d_extr(nxgs_to_spline), FV_u_extr(nxgs_to_spline), FA_u_extr(nxgs_to_spline), FTV_u_extr(nxgs_to_spline), FTA_u_extr(nxgs_to_spline);

 
  
  for(int ixg=0;ixg<nxgs_to_spline;ixg++) {
    for(int ifit=0;ifit<Nfits;ifit++) {
      FV_extr[ixg].emplace_back(1); FA_extr[ixg].emplace_back(1); FTV_extr[ixg].emplace_back(1); FTA_extr[ixg].emplace_back(1);
      FV_d_extr[ixg].emplace_back(1); FA_d_extr[ixg].emplace_back(1); FTV_d_extr[ixg].emplace_back(1); FTA_d_extr[ixg].emplace_back(1);
      FV_u_extr[ixg].emplace_back(1); FA_u_extr[ixg].emplace_back(1); FTV_u_extr[ixg].emplace_back(1); FTA_u_extr[ixg].emplace_back(1);
    }
  }

 
  vector<vector<vector<distr_t>>> FV_d_sims_m_dep(num_xg), FA_d_sims_m_dep(num_xg), FTV_d_sims_m_dep(num_xg), FTA_d_sims_m_dep(num_xg);
  vector<vector<vector<distr_t>>> FV_u_sims_m_dep(num_xg), FA_u_sims_m_dep(num_xg), FTV_u_sims_m_dep(num_xg), FTA_u_sims_m_dep(num_xg);
  vector<vector<vector<distr_t>>> FV_sims_m_dep(num_xg), FA_sims_m_dep(num_xg), FTV_sims_m_dep(num_xg), FTA_sims_m_dep(num_xg);


  for(int ixg=0;ixg<num_xg;ixg++) {
    FV_d_sims_m_dep[ixg].resize(Inv_Masses_to_print_MEX.size());
    FA_d_sims_m_dep[ixg].resize(Inv_Masses_to_print_MEX.size());
    FTV_d_sims_m_dep[ixg].resize(Inv_Masses_to_print_MEX.size());
    FTA_d_sims_m_dep[ixg].resize(Inv_Masses_to_print_MEX.size());
    
    FV_u_sims_m_dep[ixg].resize(Inv_Masses_to_print_MEX.size());
    FA_u_sims_m_dep[ixg].resize(Inv_Masses_to_print_MEX.size());
    FTV_u_sims_m_dep[ixg].resize(Inv_Masses_to_print_MEX.size());
    FTA_u_sims_m_dep[ixg].resize(Inv_Masses_to_print_MEX.size());

    FV_sims_m_dep[ixg].resize(Inv_Masses_to_print_MEX.size());
    FA_sims_m_dep[ixg].resize(Inv_Masses_to_print_MEX.size());
    FTV_sims_m_dep[ixg].resize(Inv_Masses_to_print_MEX.size());
    FTA_sims_m_dep[ixg].resize(Inv_Masses_to_print_MEX.size());
  }



  //COMPUTE COVARIANCE MATRIX
  cout<<"Computing covariance matrix of global fit"<<endl;

  vector<vector<distr_t_list>> FF_M_list_boot_new= FF_M_list_boot;
  vector<vector<distr_t_list>> FF_M_list_new = FF_M_list;
  //{"FA", "FA_u", "FA_d", "FV",  "FV_u", "FV_d",  "FA_T", "FA_T_u", "FA_T_d", "FV_T", "FV_T_u", "FV_T_d", "FB", "FB_u", "FB_d", "FT", "FT_u", "FT_d"}
  int id_TA=8; int id_TV=11; 
  if(fit_FB_FT) {
    id_TA=14;
    id_TV=17;
  }

  if(!separate_ud) {

    for(int imass=0;imass<Nmasses;imass++) {
      FF_M_list_boot_new[2][imass] = FF_M_list_boot_new[2][imass] + FF_M_list_boot_new[1][imass];
      //FF_M_list_boot_new[1][imass] = 0.0*FF_M_list_boot_new[1][imass];
      FF_M_list_new[2][imass] = FF_M_list_new[2][imass] + FF_M_list_new[1][imass];
      //FF_M_list_new[1][imass] = 0.0*FF_M_list_new[1][imass];

      FF_M_list_boot_new[5][imass] = FF_M_list_boot_new[5][imass] + FF_M_list_boot_new[4][imass];
      //FF_M_list_boot_new[4][imass] = 0.0*FF_M_list_boot_new[4][imass];
      FF_M_list_new[5][imass] = FF_M_list_new[5][imass] + FF_M_list_new[4][imass];
      //FF_M_list_new[4][imass] = 0.0*FF_M_list_new[5][imass];

      FF_M_list_boot_new[id_TA][imass] = FF_M_list_boot_new[id_TA][imass] + FF_M_list_boot_new[id_TA-1][imass];
      //FF_M_list_boot_new[id_TA-1][imass] = 0.0*FF_M_list_boot_new[id_TA-1][imass];
      FF_M_list_new[id_TA][imass] = FF_M_list_new[id_TA][imass] + FF_M_list_new[id_TA-1][imass];
      //FF_M_list_new[id_TA-1][imass] = 0.0*FF_M_list_new[id_TA-1][imass];

      FF_M_list_boot_new[id_TV][imass] = FF_M_list_boot_new[id_TV][imass] + FF_M_list_boot_new[id_TV-1][imass];
      //FF_M_list_boot_new[id_TV-1][imass] = 0.0*FF_M_list_boot_new[id_TV-1][imass];
      FF_M_list_new[id_TV][imass] = FF_M_list_new[id_TV][imass] + FF_M_list_new[id_TV-1][imass];
      //FF_M_list_new[id_TV-1][imass] = 0.0*FF_M_list_new[id_TV-1][imass];

    }

  }

    //insert covariance matrix
  Eigen::MatrixXd Cov_Matrix_comb_FF(Nmasses*8*num_xg,Nmasses*8*num_xg);
  Eigen::MatrixXd Corr_Matrix_comb_FF(Nmasses*8*num_xg,Nmasses*8*num_xg);
  for(int m=0;m< Nmasses*8*num_xg;m++)
    for(int n=0;n< Nmasses*8*num_xg ;n++) { Cov_Matrix_comb_FF(m,n)=0; Corr_Matrix_comb_FF(m,n)=0;}
  
  for(int i=0;i<Nmasses;i++)
    for(int j=0;j<Nmasses;j++) 
      for(int ix=0; ix < num_xg; ix++)
	for(int jx=0; jx<num_xg;jx++) {
	  Corr_Matrix_comb_FF(i*8*num_xg + 8*ix , j*8*num_xg + 8*jx) = ((FF_M_list_boot_new[2][i].distr_list[ix]/fps_boot.distr_list[i])%(FF_M_list_boot_new[2][j].distr_list[jx]/fps_boot.distr_list[j]))/( (FF_M_list_boot_new[2][i].distr_list[ix]/fps_boot.distr_list[i]).err()*(FF_M_list_boot_new[2][j].distr_list[jx]/fps_boot.distr_list[j]).err());
	  Corr_Matrix_comb_FF(i*8*num_xg + 8*ix + 1 , j*8*num_xg + 8*jx + 1) = ((FF_M_list_boot_new[5][i].distr_list[ix]/fps_boot.distr_list[i])%(FF_M_list_boot_new[5][j].distr_list[jx]/fps_boot.distr_list[j]))/( (FF_M_list_boot_new[5][i].distr_list[ix]/fps_boot.distr_list[i]).err()*(FF_M_list_boot_new[5][j].distr_list[jx]/fps_boot.distr_list[j]).err());
	  Corr_Matrix_comb_FF(i*8*num_xg + 8*ix + 2 , j*8*num_xg + 8*jx + 2) =   ((FF_M_list_boot_new[1][i].distr_list[ix]/fps_boot.distr_list[i])%(FF_M_list_boot_new[1][j].distr_list[jx]/fps_boot.distr_list[j]))/( (FF_M_list_boot_new[1][i].distr_list[ix]/fps_boot.distr_list[i]).err()*(FF_M_list_boot_new[1][j].distr_list[jx]/fps_boot.distr_list[j]).err());
	  Corr_Matrix_comb_FF(i*8*num_xg + 8*ix + 3 , j*8*num_xg + 8*jx + 3) =    ((FF_M_list_boot_new[4][i].distr_list[ix]/fps_boot.distr_list[i])%(FF_M_list_boot_new[4][j].distr_list[jx]/fps_boot.distr_list[j]))/( (FF_M_list_boot_new[4][i].distr_list[ix]/fps_boot.distr_list[i]).err()*(FF_M_list_boot_new[4][j].distr_list[jx]/fps_boot.distr_list[j]).err());
 
	  //tensor
	  Corr_Matrix_comb_FF(i*8*num_xg + 8*ix+4 , j*8*num_xg + 8*jx+4) = ((FF_M_list_boot_new[id_TA][i].distr_list[ix]/fps_boot.distr_list[i])%(FF_M_list_boot_new[id_TA][j].distr_list[jx]/fps_boot.distr_list[j]))/( (FF_M_list_boot_new[id_TA][i].distr_list[ix]/fps_boot.distr_list[i]).err()*(FF_M_list_boot_new[id_TA][j].distr_list[jx]/fps_boot.distr_list[j]).err());
	  Corr_Matrix_comb_FF(i*8*num_xg + 8*ix + 5 , j*8*num_xg + 8*jx + 5) = ((FF_M_list_boot_new[id_TV][i].distr_list[ix]/fps_boot.distr_list[i])%(FF_M_list_boot_new[id_TV][j].distr_list[jx]/fps_boot.distr_list[j]))/( (FF_M_list_boot_new[id_TV][i].distr_list[ix]/fps_boot.distr_list[i]).err()*(FF_M_list_boot_new[id_TV][j].distr_list[jx]/fps_boot.distr_list[j]).err());
	  Corr_Matrix_comb_FF(i*8*num_xg + 8*ix + 6 , j*8*num_xg + 8*jx + 6) =   ((FF_M_list_boot_new[id_TA-1][i].distr_list[ix]/fps_boot.distr_list[i])%(FF_M_list_boot_new[id_TA-1][j].distr_list[jx]/fps_boot.distr_list[j]))/( (FF_M_list_boot_new[id_TA-1][i].distr_list[ix]/fps_boot.distr_list[i]).err()*(FF_M_list_boot_new[id_TA-1][j].distr_list[jx]/fps_boot.distr_list[j]).err());
	  Corr_Matrix_comb_FF(i*8*num_xg + 8*ix + 7 , j*8*num_xg + 8*jx + 7) =    ((FF_M_list_boot_new[id_TV-1][i].distr_list[ix]/fps_boot.distr_list[i])%(FF_M_list_boot_new[id_TV-1][j].distr_list[jx]/fps_boot.distr_list[j]))/( (FF_M_list_boot_new[id_TV-1][i].distr_list[ix]/fps_boot.distr_list[i]).err()*(FF_M_list_boot_new[id_TV-1][j].distr_list[jx]/fps_boot.distr_list[j]).err());

	  
	  
	  Cov_Matrix_comb_FF(i*8*num_xg + 8*ix , j*8*num_xg + 8*jx) = Corr_Matrix_comb_FF(i*8*num_xg + 8*ix , j*8*num_xg + 8*jx)*(FF_M_list_new[2][i].distr_list[ix]/fps.distr_list[i]).err()*(FF_M_list_new[2][j].distr_list[jx]/fps.distr_list[j]).err();
	  Cov_Matrix_comb_FF(i*8*num_xg + 8*ix +1 , j*8*num_xg + 8*jx +1) = Corr_Matrix_comb_FF(i*8*num_xg + 8*ix +1 , j*8*num_xg + 8*jx +1)*(FF_M_list_new[5][i].distr_list[ix]/fps.distr_list[i]).err()*(FF_M_list_new[5][j].distr_list[jx]/fps.distr_list[j]).err();
	  Cov_Matrix_comb_FF(i*8*num_xg + 8*ix +2 , j*8*num_xg + 8*jx +2) = Corr_Matrix_comb_FF(i*8*num_xg + 8*ix +2 , j*8*num_xg + 8*jx +2)*(FF_M_list_new[1][i].distr_list[ix]/fps.distr_list[i]).err()*(FF_M_list_new[1][j].distr_list[jx]/fps.distr_list[j]).err();
	  Cov_Matrix_comb_FF(i*8*num_xg + 8*ix +3 , j*8*num_xg + 8*jx +3) = Corr_Matrix_comb_FF(i*8*num_xg + 8*ix +3 , j*8*num_xg + 8*jx +3)*(FF_M_list_new[4][i].distr_list[ix]/fps.distr_list[i]).err()*(FF_M_list_new[4][j].distr_list[jx]/fps.distr_list[j]).err();

	  Cov_Matrix_comb_FF(i*8*num_xg + 8*ix+4 , j*8*num_xg + 8*jx+4) = Corr_Matrix_comb_FF(i*8*num_xg + 8*ix +4  , j*8*num_xg + 8*jx +4)*(FF_M_list_new[id_TA][i].distr_list[ix]/fps.distr_list[i]).err()*(FF_M_list_new[id_TA][j].distr_list[jx]/fps.distr_list[j]).err();
	  Cov_Matrix_comb_FF(i*8*num_xg + 8*ix +5 , j*8*num_xg + 8*jx +5) = Corr_Matrix_comb_FF(i*8*num_xg + 8*ix +5 , j*8*num_xg + 8*jx +5)*(FF_M_list_new[id_TV][i].distr_list[ix]/fps.distr_list[i]).err()*(FF_M_list_new[id_TV][j].distr_list[jx]/fps.distr_list[j]).err();
	  Cov_Matrix_comb_FF(i*8*num_xg + 8*ix +6 , j*8*num_xg + 8*jx +6) = Corr_Matrix_comb_FF(i*8*num_xg + 8*ix +6 , j*8*num_xg + 8*jx +6)*(FF_M_list_new[id_TA-1][i].distr_list[ix]/fps.distr_list[i]).err()*(FF_M_list_new[id_TA-1][j].distr_list[jx]/fps.distr_list[j]).err();
	  Cov_Matrix_comb_FF(i*8*num_xg + 8*ix +7 , j*8*num_xg + 8*jx +7) = Corr_Matrix_comb_FF(i*8*num_xg + 8*ix +7 , j*8*num_xg + 8*jx +7)*(FF_M_list_new[id_TV-1][i].distr_list[ix]/fps.distr_list[i]).err()*(FF_M_list_new[id_TV-1][j].distr_list[jx]/fps.distr_list[j]).err();

	  
	  //Cov_Matrix_comb_FF(i*2*num_xg + 2*ix , j*2*num_xg + 2*jx + 1) = (FF_M_list_new[2][i].distr_list[ix]/fps.distr_list[i])%(FF_M_list_new[5][j].distr_list[jx]/fps.distr_list[j]);
	  //Cov_Matrix_comb_FF(i*2*num_xg + 2*ix + 1 , j*2*num_xg + 2*jx) = (FF_M_list_new[5][i].distr_list[ix]/fps.distr_list[i])%(FF_M_list_new[2][j].distr_list[jx]/fps.distr_list[j]);
	}

  //for(int m=0;m < Nmasses*4*num_xg; m++)
  //for(int n=0; n< Nmasses*4*num_xg;n++) Corr_Matrix_comb_FF(m,n) = Cov_Matrix_comb_FF(m,n)/sqrt( Cov_Matrix_comb_FF(m,m)*Cov_Matrix_comb_FF(n,n));


      
  cout<<"Corr matrix determinant: "<<Corr_Matrix_comb_FF.determinant()<<endl;    
  
  //print on screen
  cout<<"Correlation matrix comb"<<endl;
  cout<<Corr_Matrix_comb_FF<<endl;

 

  

 
  vector<int> Ndofs, Nmeasurements; vector<double> Ch2s;

  for(int ifit=0;ifit<Nfits;ifit++) {

    cout<<"ifit: "<<ifit<<endl;

  

    //index corresponding to FA_d and FV_d is 2 and 5

   
    class ipar_MEX_comb {
    public:
      ipar_MEX_comb() : FF(0.0), FF_err(0.0) {}
      double FF, FF_err, M, xg, imass, ixg, Z;
      string W;
    };
    
    class fpar_MEX_comb {
    public:
      fpar_MEX_comb() {}
    fpar_MEX_comb(const Vfloat &par) {
      if((signed)par.size() != 39) crash("In class fpar_MEX_comb  class constructor Vfloat par has size != 39");
      
      //A is 1/LB
      A=par[0];

      //parameters entering radiative corrections and related to moments of LCDA
      R1=par[1];
      R2=par[2];
      
      //quasi pole parameters DV for vector and DA for axial
      DV=par[3];
      DA=par[4];

      //common parameters in 1/M, 1/Eg corrections to FVd and FAd

      M1=par[5];
      M2=par[6];
     

      //parameters in 1/(M^2) and 1/(M^2*x*x) and 1/(M^2*x) for FVd and FAd
     
      SV3=par[7];
      SA3=par[8];

      //non-common parameters in 1/M and 1/Mx correction to LO for FVu and FAu

      BV1= par[9];
      BV2= par[10];
      BA1= par[11];
      BA2= par[12];

      //parameters for the tensor currents 
     

      //parameters entering radiative corrections
      RTV1=par[13];
      RTV2=par[14];
      RTA1=par[15];
      RTA2=par[16];

      //parameters of the 1/M and 1/Eg expansion
     
      TVM1 = par[17];
      TAE1 = par[18];
      TAM1 = par[19];

      //parameters in 1/(M^2) and 1/(M^2*x*x) and 1/(M^2*x) for FTVd and FTAd
     
      STV3=par[20];
      STA1=par[21];     
      STA3=par[22];

      //parameters entering 1/M and 1/Mx corrections to LO
      BTV1 = par[23];
      BTV2 = par[24];
      BTA1 = par[25];
      BTA2 = par[26];

      
      //prior pars 
      M3=par[27];
      M4=par[28];
      TVE1 = par[29];
      A1=par[30];
      A2=par[31];
      SV1=par[32];
      SV2=par[33];
      SA1=par[34];
      SA2=par[35];
      STV1=par[36];
      STV2=par[37];
      STA2=par[38];
      
      
     
    }
      double A,R1,R2, DV, DA,M1, M2,  SV3,   SA3, BV1, BV2, BA1, BA2;
      double  RTV1, RTV2, RTA1, RTA2, TVE1, TVM1, TAE1, TAM1, STV3, STA1, STA3, BTV1, BTV2, BTA1, BTA2;
      //prior pars
      double M3, M4, A1, A2, SV1, SV2, SA1, SA2, STV1, STV2, STA2;
    };
    
    
    //init bootstrap fit

    bootstrap_fit<fpar_MEX_comb,ipar_MEX_comb> bf_MEX_comb(NJ);
    bf_MEX_comb.set_warmup_lev(0); //sets warmup
    bf_MEX_comb.Set_number_of_measurements(Nmasses*8*num_xg);
    bf_MEX_comb.Set_verbosity(1);
    //fit on mean values to get ch2
    bootstrap_fit<fpar_MEX_comb,ipar_MEX_comb> bf_MEX_comb_ch2(1);
    bf_MEX_comb_ch2.set_warmup_lev(0); //sets warmup
    bf_MEX_comb_ch2.Set_number_of_measurements(Nmasses*8*num_xg);
    bf_MEX_comb_ch2.Set_verbosity(1);
    
    
    bf_MEX_comb.Add_par("A", 1.0/0.5, 0.1/0.5);
    bf_MEX_comb.Add_par("R1", 0.1, 0.01);
    bf_MEX_comb.Add_par("R2", 0.1, 0.01);
    bf_MEX_comb.Add_par("DV", 0.24, 0.1*0.24);
    bf_MEX_comb.Add_par("DA", 0.5, 0.1*0.5);
    bf_MEX_comb.Add_par("M1", 0.6, 0.06 );
    bf_MEX_comb.Add_par("M2", -1.0, 0.05 );
    bf_MEX_comb.Add_par("SV3", 1.0, 0.1);
    bf_MEX_comb.Add_par("SA3", 1.0, 0.1);
    bf_MEX_comb.Add_par("BV1", 1.0, 0.1);
    bf_MEX_comb.Add_par("BV2", 1.0, 0.1);
    bf_MEX_comb.Add_par("BA1", 1.0, 0.1);
    bf_MEX_comb.Add_par("BA2", 1.0, 0.1);
    
    //parameters of the tensor form factors
    bf_MEX_comb.Add_par("RTV1", 0.1, 0.01);
    bf_MEX_comb.Add_par("RTV2", 0.1, 0.01);
    bf_MEX_comb.Add_par("RTA1", 0.1, 0.01);
    bf_MEX_comb.Add_par("RTA2", 0.1, 0.01);  
    bf_MEX_comb.Add_par("TVM1", 1.0, 0.1);
    bf_MEX_comb.Add_par("TAE1", -1.0, 0.1);
    bf_MEX_comb.Add_par("TAM1", 4.0, 0.1);
    bf_MEX_comb.Add_par("STV3", 1.0, 0.1);
    bf_MEX_comb.Add_par("STA1", 1.0, 0.1);
    bf_MEX_comb.Add_par("STA3", 1.0, 0.1);
    bf_MEX_comb.Add_par("BTV1", 1.0, 0.1);
    bf_MEX_comb.Add_par("BTV2", 1.0, 0.1);
    bf_MEX_comb.Add_par("BTA1", 1.0, 0.1);
    bf_MEX_comb.Add_par("BTA2", 1.0, 0.1);
    
    bf_MEX_comb_ch2.Add_par("A", 1.0/0.5, 0.1/0.5);
    bf_MEX_comb_ch2.Add_par("R1", 0.1, 0.01);
    bf_MEX_comb_ch2.Add_par("R2", 0.1, 0.01);
    bf_MEX_comb_ch2.Add_par("DV", 0.24, 0.1*0.24);
    bf_MEX_comb_ch2.Add_par("DA", 0.5, 0.1*0.5);
    bf_MEX_comb_ch2.Add_par("M1", 0.6, 0.06 );
    bf_MEX_comb_ch2.Add_par("M2", -1.0, 0.05 );
    bf_MEX_comb_ch2.Add_par("SV3", 1.0, 0.1);
    bf_MEX_comb_ch2.Add_par("SA3", 1.0, 0.1);
    bf_MEX_comb_ch2.Add_par("BV1", 1.0, 0.1);
    bf_MEX_comb_ch2.Add_par("BV2", 1.0, 0.1);
    bf_MEX_comb_ch2.Add_par("BA1", 1.0, 0.1);
    bf_MEX_comb_ch2.Add_par("BA2", 1.0, 0.1);
  //parameters of the tensor form factors
    bf_MEX_comb_ch2.Add_par("RTV1", 0.1, 0.01);
    bf_MEX_comb_ch2.Add_par("RTV2", 0.1, 0.01);
    bf_MEX_comb_ch2.Add_par("RTA1", 0.1, 0.01);
    bf_MEX_comb_ch2.Add_par("RTA2", 0.1, 0.01);
    bf_MEX_comb_ch2.Add_par("TVM1", 1.0, 0.1);
    bf_MEX_comb_ch2.Add_par("TAE1", -1.0, 0.1);
    bf_MEX_comb_ch2.Add_par("TAM1", 4.0, 0.1);
    bf_MEX_comb_ch2.Add_par("STV3", 1.0, 0.1);
    bf_MEX_comb_ch2.Add_par("STA1", 1.0, 0.1);
    bf_MEX_comb_ch2.Add_par("STA3", 1.0, 0.1);
    bf_MEX_comb_ch2.Add_par("BTV1", 1.0, 0.1);
    bf_MEX_comb_ch2.Add_par("BTV2", 1.0, 0.1);
    bf_MEX_comb_ch2.Add_par("BTA1", 1.0, 0.1);
    bf_MEX_comb_ch2.Add_par("BTA2", 1.0, 0.1);


    bf_MEX_comb.Add_par("M3", -0.2, 0.01); bf_MEX_comb_ch2.Add_par("M3", -0.2, 0.01);
    bf_MEX_comb.Add_par("M4", -1.2, 0.1); bf_MEX_comb_ch2.Add_par("M4", -1.2,0.1);
    bf_MEX_comb.Add_par("TVE1", 0.1, 0.01);   bf_MEX_comb_ch2.Add_par("TVE1", 0.1, 0.01);
    bf_MEX_comb.Add_par("A1",0.1*1.0/0.5, 0.1*0.1/0.5); bf_MEX_comb_ch2.Add_par("A1",0.1*1.0/0.5, 0.1*0.1/0.5);
    bf_MEX_comb.Add_par("A2",0.0, 0.1*0.1); bf_MEX_comb_ch2.Add_par("A2",0.0, 0.1*0.1);

    bf_MEX_comb.Add_par("SV1",0.0, 0.1); bf_MEX_comb_ch2.Add_par("SV1",0.0, 0.1);
    bf_MEX_comb.Add_par("SV2",0.0, 0.1); bf_MEX_comb_ch2.Add_par("SV2",0.0, 0.1);
    bf_MEX_comb.Add_par("SA1",0.0, 0.1); bf_MEX_comb_ch2.Add_par("SA1",0.0, 0.1);
    bf_MEX_comb.Add_par("SA2",0.0, 0.1); bf_MEX_comb_ch2.Add_par("SA2",0.0, 0.1);
    bf_MEX_comb.Add_par("STV1",0.0, 0.1); bf_MEX_comb_ch2.Add_par("STV1",0.0, 0.1);
    bf_MEX_comb.Add_par("STV2",0.0, 0.1); bf_MEX_comb_ch2.Add_par("STV2",0.0, 0.1);
    bf_MEX_comb.Add_par("STA2",0.0, 0.1); bf_MEX_comb_ch2.Add_par("STA2",0.0, 0.1);
    
   
    //add priors
    map<string,int>::iterator it_prior;
    it_prior= Fix_fit_pars_list[ifit].find("M3");
    if(it_prior != Fix_fit_pars_list[ifit].end()) {
      if(it_prior->second == 1) {
	bf_MEX_comb.Add_prior_par("M3", -0.2, 0.01); bf_MEX_comb_ch2.Add_prior_par("M3", -0.2, 0.01);
      }
    }

    it_prior= Fix_fit_pars_list[ifit].find("M4");
    if(it_prior != Fix_fit_pars_list[ifit].end()) {
      if(it_prior->second == 1) {
	bf_MEX_comb.Add_prior_par("M4", -1.2, 0.1); bf_MEX_comb_ch2.Add_prior_par("M4", -1.2,0.1);
      }
    }

    it_prior= Fix_fit_pars_list[ifit].find("TVE1");
    if(it_prior != Fix_fit_pars_list[ifit].end()) {
      if(it_prior->second == 1) {
	bf_MEX_comb.Add_prior_par("TVE1", 0.1, 0.01); bf_MEX_comb_ch2.Add_prior_par("TVE1", 0.1,0.01);
      }
    }

  
    
    it_prior= Fix_fit_pars_list[ifit].find("A1");
    if(it_prior != Fix_fit_pars_list[ifit].end()) {
      if(it_prior->second == 1) {
	bf_MEX_comb.Add_prior_par("A1", 0.1*1.0/0.5, 0.1*0.1/0.5); bf_MEX_comb_ch2.Add_prior_par("A1", 0.1*1.0/0.5, 0.1*0.1/0.5);
      }
    }

    it_prior= Fix_fit_pars_list[ifit].find("A2");
    if(it_prior != Fix_fit_pars_list[ifit].end()) {
      if(it_prior->second == 1) {
	bf_MEX_comb.Add_prior_par("A2", 0.0, 0.1); bf_MEX_comb_ch2.Add_prior_par("A2", 0.0, 0.1);
      }
    }    
    vector<string> M2_pars({"SV1", "SV2", "SA1", "SA2", "STV1", "STV2", "STA2"});
    for( auto &mpar : M2_pars) {
      it_prior= Fix_fit_pars_list[ifit].find(mpar);
      if(it_prior == Fix_fit_pars_list[ifit].end()) crash("Cannot find par: "+mpar+" in fit list");
      if(it_prior != Fix_fit_pars_list[ifit].end()) {
	if(it_prior->second == 1) {
	  cout<<"Adding prior: "<<mpar<<endl;
	  bf_MEX_comb.Add_prior_par(mpar, 0.0, 0.1); bf_MEX_comb_ch2.Add_prior_par(mpar, 0.0, 0.1);
	}
      }
    }
    
  

    
 
  //############### FIXING PARAMETERS ####################
  

  bf_MEX_comb.Fix_par("R1", 0.0); bf_MEX_comb_ch2.Fix_par("R1", 0.0);
  bf_MEX_comb.Fix_par("R2", 0.0); bf_MEX_comb_ch2.Fix_par("R2", 0.0);
  bf_MEX_comb.Fix_par("SV3", 0.0); bf_MEX_comb_ch2.Fix_par("SV3", 0.0);
  bf_MEX_comb.Fix_par("SA3", 0.0);  bf_MEX_comb_ch2.Fix_par("SA3", 0.0);
   
  //fix parameters of tensor form factors
  bf_MEX_comb.Fix_par("RTV1", 0.0); bf_MEX_comb_ch2.Fix_par("RTV1", 0.0);
  bf_MEX_comb.Fix_par("RTV2", 0.0); bf_MEX_comb_ch2.Fix_par("RTV2", 0.0);
  bf_MEX_comb.Fix_par("RTA1", 0.0); bf_MEX_comb_ch2.Fix_par("RTA1", 0.0);
  bf_MEX_comb.Fix_par("RTA2", 0.0); bf_MEX_comb_ch2.Fix_par("RTA2", 0.0);
  bf_MEX_comb.Fix_par("TAM1", 0.0); bf_MEX_comb_ch2.Fix_par("TAM1", 0.0);
  bf_MEX_comb.Fix_par("STA1", 0.0);  bf_MEX_comb_ch2.Fix_par("STA1", 0.0);
  bf_MEX_comb.Fix_par("STV3", 0.0); bf_MEX_comb_ch2.Fix_par("STV3", 0.0);
  bf_MEX_comb.Fix_par("STA3", 0.0);  bf_MEX_comb_ch2.Fix_par("STA3", 0.0);

  //bf_MEX_comb.Fix_par("TVM1", 0.0); bf_MEX_comb_ch2.Fix_par("TVM1", 0.0);
  //bf_MEX_comb.Fix_par("TAE1", 0.0); bf_MEX_comb_ch2.Fix_par("TAE1", 0.0);

  int fix_T_den= Fix_fit_pars_list[ifit].find("fix_T_den")->second;
  
  map<string, int>::iterator it;
  for (it = Fix_fit_pars_list[ifit].begin(); it != Fix_fit_pars_list[ifit].end(); it++)
    {
      if (it->first != "fix_T_den") {
	if(it->second == 0) { bf_MEX_comb.Fix_par(it->first, 0.0); bf_MEX_comb_ch2.Fix_par(it->first, 0.0);}
      }
    }
  
  if(!separate_ud) {
    bf_MEX_comb.Fix_par("BV1", 0.0);  bf_MEX_comb_ch2.Fix_par("BV1", 0.0);
    bf_MEX_comb.Fix_par("BV2", 0.0);  bf_MEX_comb_ch2.Fix_par("BV2", 0.0);
    bf_MEX_comb.Fix_par("BA1", 0.0);  bf_MEX_comb_ch2.Fix_par("BA1", 0.0);
    bf_MEX_comb.Fix_par("BA2", 0.0);  bf_MEX_comb_ch2.Fix_par("BA2", 0.0);

    bf_MEX_comb.Fix_par("BTV1", 0.0) ; bf_MEX_comb_ch2.Fix_par("BTV1", 0.0);
    bf_MEX_comb.Fix_par("BTV2", 0.0) ; bf_MEX_comb_ch2.Fix_par("BTV2", 0.0);
    bf_MEX_comb.Fix_par("BTA1", 0.0) ; bf_MEX_comb_ch2.Fix_par("BTA1", 0.0);
    bf_MEX_comb.Fix_par("BTA2", 0.0) ; bf_MEX_comb_ch2.Fix_par("BTA2", 0.0);
  }
  
   
 
  double L_MSb=Get_Lambda_MS_bar(4);
  double alpha_MU_GLB=Get_4l_alpha_s(MU_GLB,4,L_MSb);
  double zeta_3= riemann_zeta(3);
  double zeta_4 = riemann_zeta(4);
  double zeta_5 = riemann_zeta(5);
  int Fit_mode=0; //options : 0 = use mh1=mh2=mu=MU_GLB , : 1 = use mh1= mh2 =M, mu=MU_GLB,  : 2 = used mh1=x*M, mh2=M, mu=MU_GLB
  int Nf=4;
  double b0= (33.0-2.0*Nf)/(12.0*M_PI);
  double b1= (153.0 - 19*Nf)/(24.0*M_PI*M_PI);
  double b2= (2857.0 -5033.0*Nf/9.0 + 325.0*Nf*Nf/27.0)/(128.0*pow(M_PI,3));
  double b3= ( (149753.0/6.0 + 3564*zeta_3) -(1078361.0/162.0 + 6508.0*zeta_3/27.0)*Nf + (50065.0/162.0 + 6472.0*zeta_3/81.0)*Nf*Nf + 1093.0*pow(Nf,3)/729.0)/(256.0*pow(M_PI,4));
  double b4= (1.0/pow(4*M_PI,5))*( 8157455.0/16 + 621885.0*zeta_3/2.0 - 88209.0*zeta_4/2 -288090*zeta_5 + Nf*( -336460813.0/1944 - 4811164*zeta_3/81 + 33935*zeta_4/6 + 1358995*zeta_5/27) + pow(Nf,2)*( 25960913.0/1944 + 698531.0*zeta_3/81 - 10526*zeta_4/9 - 381760.0*zeta_5/81) + pow(Nf,3)*( -630559.9/5832 - 48722.0*zeta_3/243 + 1618.0*zeta_4/27 + 460.0*zeta_5/9) + pow(Nf,4)*( 1205.0/2916 -152.0*zeta_3/81))   ;

 

  bf_MEX_comb.ansatz =  [&L_MSb, &alpha_MU_GLB, &zeta_3, &zeta_4, &zeta_5, &b0, &b1,&b2, &b3, &b4,  &Fit_mode, &ZT_MB, &separate_ud, &fix_T_den](const fpar_MEX_comb &p, const ipar_MEX_comb &ip) {

    if((separate_ud == 0)  && ip.W.substr( ip.W.length() -1, 1) == "u") return 0.0;

    //if(ip.W.substr(0,2)=="TA") return 0.0;

    //if(ip.W=="Vu" || ip.W == "Au") return 0.0;

    double LB=0.5;
    double M=ip.M;
    double x = ip.xg;
    double Q= 1.0/3.0;

    /*
    
    double t_Mh=log(M*M/(L_MSb*L_MSb));
    double l_Mh=log(t_Mh);
    double t_Eg=log(M*M*x*x/(L_MSb*L_MSb));
    double l_Eg=log(t_Eg);
    //################### COMPUTATION OF RADIATIVE CORRECTIONS ######################
    //double alpha_Mh= alphas_Mh[ip.imass];
    //double alpha_Eg= alphas_Eg[ip.imass][ip.ixg];


    double alpha_Mh = (1.0/(b0*t_Mh))*( 1 - b1*l_Mh/(b0*b0*t_Mh) + (pow(b1,2)*(pow(l_Mh,2) -l_Mh -1.0) + b0*b2)/(pow(b0,4)*pow(t_Mh,2)));
    alpha_Mh = alpha_Mh + (1.0/(b0*t_Mh))*(  (pow(b1,3)*(-2*pow(l_Mh,3)+5*pow(l_Mh,2)+ 4.0*l_Mh -1) - 6.0*b0*b2*b1*l_Mh + b0*b0*b3)/(2.0*pow(b0,6)*pow(t_Mh,3)));
    alpha_Mh = alpha_Mh + (1.0/(b0*t_Mh))*( (18.0*b0*b2*pow(b1,2)*(2*l_Mh*l_Mh -l_Mh -1.0) + pow(b1,4)*(6*pow(l_Mh,4) -26*pow(l_Mh,3) -9*l_Mh*l_Mh + 24*l_Mh + 7.0))/(6.0*pow(b0,8)*pow(t_Mh,4)));
    alpha_Mh = alpha_Mh + (1.0/(b0*t_Mh))*( ( -pow(b0,2)*b3*b1*(12*l_Mh+1) + 2*pow(b0,2)*5*pow(b2,2))/(6*pow(b0,8)*pow(t_Mh,4)));
    alpha_Mh = alpha_Mh + (1.0/(b0*t_Mh))*( 2*pow(b0,2)*b0*b4)/(6*pow(b0,8)*pow(t_Mh,4));

    double alpha_Eg = (1.0/(b0*t_Eg))*( 1 - b1*l_Eg/(b0*b0*t_Eg) + (pow(b1,2)*(pow(l_Eg,2) -l_Eg -1.0) + b0*b2)/(pow(b0,4)*pow(t_Eg,2)));
    alpha_Eg = alpha_Eg + (1.0/(b0*t_Eg))*(  (pow(b1,3)*(-2*pow(l_Eg,3)+5*pow(l_Eg,2)+ 4.0*l_Eg -1) - 6.0*b0*b2*b1*l_Eg + b0*b0*b3)/(2.0*pow(b0,6)*pow(t_Eg,3)));
    alpha_Eg = alpha_Eg + (1.0/(b0*t_Eg))*( (18.0*b0*b2*pow(b1,2)*(2*l_Eg*l_Eg -l_Eg -1.0) + pow(b1,4)*(6*pow(l_Eg,4) -26*pow(l_Eg,3) -9*l_Eg*l_Eg + 24*l_Eg + 7.0))/(6.0*pow(b0,8)*pow(t_Eg,4)));
    alpha_Eg = alpha_Eg + (1.0/(b0*t_Eg))*( ( -pow(b0,2)*b3*b1*(12*l_Eg+1) + 2*pow(b0,2)*5*pow(b2,2))/(6*pow(b0,8)*pow(t_Eg,4)));
    alpha_Eg = alpha_Eg + (1.0/(b0*t_Eg))*( 2*pow(b0,2)*b0*b4)/(6*pow(b0,8)*pow(t_Eg,4));
    

    double alpha_mh1=0;
    double alpha_mh2=0;
    double mh1=0;
    double mh2=0;

    if(Fit_mode ==0) {

      mh1=MU_GLB; mh2=MU_GLB;
      alpha_mh1= alpha_MU_GLB; alpha_mh2=alpha_MU_GLB;
    }
    else if(Fit_mode==1) {
      mh1=M; mh2=M;
      alpha_mh1= alpha_Mh; alpha_mh2=alpha_Mh;
    }
    else if(Fit_mode==2) {
      mh1=M*x; mh2=M;
      alpha_mh1= alpha_Eg; alpha_mh2=alpha_Mh;
    }
    else crash("Fit_mode: "+to_string(Fit_mode)+" not yet implemented");

        
    double CF= 4.0/3;
    double C= 1.0 + alpha_mh1*(CF/(4.0*M_PI))*( -2*pow(log(x*M/mh1),2) + 5*log(x*M/mh1) -(3-2*x)*log(x)/(1-x) -6 -M_PI*M_PI/12.0 - 2*gsl_sf_dilog(1-x));
    double K= 1.0 +alpha_mh2*CF*( 3*log(M/mh2)  -2 )/(4.0*M_PI);
    double J= 1.0 + alpha_MU_GLB*(CF/(4.0*M_PI))*( pow(log(x*M/(MU_GLB*MU_GLB)),2) -2*p.R1*log(x*M/(MU_GLB*MU_GLB)) -1 -M_PI*M_PI/6.0  + p.R2);


  
    //evolutor
    double B0= pow(4*M_PI,1)*b0;
    double B1= pow(4*M_PI,2)*b1;
    double B2= pow(4*M_PI,3)*b2;
      
    double r1 = alpha_MU_GLB/alpha_mh1;
            
    double g0 = -5*CF;
    double g1 = CF*(-1585.0/18 -5*M_PI*M_PI/6.0 + 34*zeta_3 + 4*(125.0/27 + M_PI*M_PI/3.0));
   
    double G0 = 4*CF;
    double G1 = CF*(268.0/3 - 4*M_PI*M_PI -40.0*4.0/9.0);
    double G2 = CF*( 1470.0 -536.0*M_PI*M_PI/3.0 +44.0*pow(M_PI,4)/5.0 +264*zeta_3 + 4*( -1276.0/9.0 +80*M_PI*M_PI/9.0 -208*zeta_3/3.0) -16.0*pow(4.0,2)/27.0);
    
    double U1 =   exp( -(G0/(4*B0*B0))*( (4*M_PI/alpha_mh1)*( log(r1) - 1 + 1.0/r1) - (B1/(2*B0))*log(r1)*log(r1) + ( (G1/G0) - (B1/B0))*(r1-1-log(r1) ) ));
    U1 *= pow(M*x/mh1, -G0*log(r1)/(2*B0))*pow(r1, -g0/(2*B0))*( 1 - (alpha_mh1/(4*M_PI))*(G0/(4*B0*B0))*( (G2/(2*G0))*pow(1-r1,2) + (B2/(2*B0))*(1 - r1*r1 + 2*log(r1)) - (G1/(2*G0))*(B1/B0)*( 3 - 4*r1 + r1*r1 + 2*r1*log(r1)) + (B1*B1/(2*B0*B0))*(1-r1)*(1 - r1 -2*log(r1)) ) + (alpha_mh1/(4*M_PI))*( log(M*x/mh1)*( G1/(2*B0) - G0*B1/(2*B0*B0) ) + g1/(2*B0) - g0*B1/(2*B0*B0))*(1-r1));

    double r2 = alpha_MU_GLB/alpha_mh2;
    double g0_hl = -3*CF;
    double g1_hl = CF*( -127.0/6.0 -14*M_PI*M_PI/9.0 +5.0*4.0/3.0);
    double U2= pow(r2, -g0_hl/(2*B0))*( 1.0 +(alpha_mh2/(4.0*M_PI))*(1-r2)*( g1_hl/(2*B0) -g0_hl*B1/(2*B0*B0) ));

    
    double Rad= (C*J/K)*U1/U2; */
    double Rad=    1 - p.R1*log(M/L_MSb) + -p.R2*x*log(M/L_MSb);
    double Rad_TA= 1 - p.RTA1*log(M/L_MSb) - p.RTA2*x*log(M/L_MSb);  
    double Rad_TV= 1 - p.RTV1*log(M/L_MSb) - p.RTV2*x*log(M/L_MSb);
    //cout<<"C("<<ip.xg<<","<<ip.M<<"): "<<C<<endl;
    //cout<<"K("<<ip.xg<<","<<ip.M<<"): "<<K<<endl;
    //cout<<"J("<<ip.xg<<","<<ip.M<<"): "<<J<<endl;
    //cout<<"R("<<ip.xg<<","<<ip.M<<"): "<<Rad<<endl;
    //################### COMPUTATION OF RADIATIVE CORRECTIONS ######################
   

          
    if (ip.W =="Vd") {
      return (Q/x)*(1.0/(1 + 2*p.DV/(x*M*M)))*( Rad*p.A + (p.M1)/M + (p.M2)/(M*x) + (1.0+p.M3)/(M*x) + (1.0)/(M +p.M4*LB)  + p.SV1*LB/(M*M) + p.SV2*LB/(M*M*x) + p.SV3*LB/(M*M*x*x));   
    }
    else if (ip.W == "Ad") {
      return (Q/x)*(1.0/( 1  +  2*p.DA/(x*M) ) )*( Rad*p.A + (p.M1)/M + p.M2/(M*x) + 2*p.DA*Rad*p.A/(M*x) -(1.0+p.M3)/(M*x) -(1.0)/(M+p.M4*LB) + p.SA1*LB/(M*M) + p.SA2*LB/(M*M*x) + p.SA3*LB/(M*M*x*x));   
    }
    else if(ip.W =="Vu") { //return 0.0; 
      return (Q/(M*x))*(1.0/(1+2*p.DV/(x*M*M)))*(1.0 + p.BV1*LB/M + p.BV2*LB/(M*x) );
    }

    else if(ip.W == "Au") {// return 0.0;
      return (Q/(M*x))*(1.0/(1+2*p.DA/(x*M)))*(1.0 + p.BA1*LB/M + p.BA2*LB/(M*x) );
    }
    else if (ip.W =="TVd") {
      return (1.0/ip.Z)*(Q/x)*((1.0 + (1-fix_T_den)*2*p.DV/(M*M))/(1 + (1-fix_T_den*x)*2*p.DV/(x*M*M)))*( Rad_TV*(p.A1+p.A) + (p.TVM1+1.0)/M + p.TAE1/(M*x) + (1-x)*(1.0+p.TVE1)/(M*x) + p.STV1*LB/(M*M) + p.STV2*(1-x)*LB/(M*M*x) + p.STV3*LB/(M*M*x*x));   
    }
    else if (ip.W == "TAd") {
      return (1.0/ip.Z)*(Q/x)*((1.0+  (1-fix_T_den)*2*((p.DA+p.A2)/M))/( 1  + (1-fix_T_den*x)*2*(p.DA+p.A2)/(x*M) ) )*( Rad_TA*(p.A+p.A1) + (p.TVM1+1.0)/M + p.TAE1/(M*x) -(1-x)*(1.0+p.TVE1-2*(p.DA+p.A2)*Rad_TA*(p.A+p.A1))/(M*x) + p.STV1*LB/(M*M) + (p.STA2)*(1-x)*LB/(M*M*x) + p.STA3*LB/(M*M*x*x));   
    }
    else if(ip.W =="TVu") { //return 0.0; 
      return (1.0/ip.Z)*(Q/(M*x))*(1.0/(1+2*p.DV/(x*M*M)))*(1.0 + p.BTV1*LB/M + p.BTV2*LB/(M*x) );
    }
    else if(ip.W == "TAu") {// return 0.0;
      return (1.0/ip.Z)*(Q/(M*x))*(1.0/(1+2*p.DA/(x*M)))*(1.0 + p.BTA1*LB/M + p.BTA2*LB/(M*x) );
    }
    else crash("channel: "+ip.W+" not recognized");

        
    return 0.0;
  };

	
  bf_MEX_comb.measurement=  [&separate_ud ](const fpar_MEX_comb &p, const ipar_MEX_comb &ip) {

    if((separate_ud == 0)  && ip.W.substr( ip.W.length() -1, 1) == "u") return 0.0;

    //if(ip.W.substr(0,2)=="TA") return 0.0;
    
    //if(ip.W=="Vu" || ip.W == "Au") return 0.0;
    return ip.FF;
    
  };
	
  bf_MEX_comb.error=  [ ](const fpar_MEX_comb &p, const ipar_MEX_comb &ip) {
	  
    return ip.FF_err;
  };

 	
  bf_MEX_comb_ch2.ansatz= bf_MEX_comb.ansatz;
  bf_MEX_comb_ch2.measurement = bf_MEX_comb.measurement;
  bf_MEX_comb_ch2.error = bf_MEX_comb.error;





  
  //add cov matrix to bootstrap fit
  bf_MEX_comb.Add_covariance_matrix(Cov_Matrix_comb_FF);
  bf_MEX_comb_ch2.Add_covariance_matrix(Cov_Matrix_comb_FF);




  
  
  //fill the data
  vector<vector<ipar_MEX_comb>> data_MEX_comb(NJ);
  vector<vector<ipar_MEX_comb>> data_MEX_comb_ch2(1);
  //allocate space for output result
  boot_fit_data<fpar_MEX_comb> Bt_fit_MEX_comb;
  boot_fit_data<fpar_MEX_comb> Bt_fit_MEX_comb_ch2;
  for(auto &data_iboot: data_MEX_comb) data_iboot.resize(Nmasses*8*num_xg);
  for(auto &data_iboot: data_MEX_comb_ch2) data_iboot.resize(Nmasses*8*num_xg);
  for(int ijack=0;ijack<NJ;ijack++) {
    for(int im=0;im<Nmasses;im++) {
      for(int ixg=0;ixg<num_xg;ixg++) {
	data_MEX_comb[ijack][8*im*num_xg+ 8*ixg].FF = (FF_M_list_new[2][im].distr_list[ixg]/fps.distr_list[im]).distr[ijack];
	data_MEX_comb[ijack][8*im*num_xg+ 8*ixg].FF_err= (FF_M_list_new[2][im].distr_list[ixg]/fps.distr_list[im]).err();
	data_MEX_comb[ijack][8*im*num_xg+ 8*ixg].xg= Bs_xg_t_list[ixg];
	data_MEX_comb[ijack][8*im*num_xg+ 8*ixg].M= masses.distr_list[im].distr[ijack];
	data_MEX_comb[ijack][8*im*num_xg+ 8*ixg].W= "Ad";
	data_MEX_comb[ijack][8*im*num_xg +8*ixg].imass = im;
	data_MEX_comb[ijack][8*im*num_xg +8*ixg].ixg = ixg;
	data_MEX_comb[ijack][8*im*num_xg +8*ixg].Z = ZT_Mh[im];

	data_MEX_comb[ijack][8*im*num_xg+8*ixg+1].FF = (FF_M_list_new[5][im].distr_list[ixg]/fps.distr_list[im]).distr[ijack];
	data_MEX_comb[ijack][8*im*num_xg+8*ixg+1].FF_err= (FF_M_list_new[5][im].distr_list[ixg]/fps.distr_list[im]).err();
	data_MEX_comb[ijack][8*im*num_xg+8*ixg+1].xg= Bs_xg_t_list[ixg];
	data_MEX_comb[ijack][8*im*num_xg+8*ixg+1].M= masses.distr_list[im].distr[ijack];
	data_MEX_comb[ijack][8*im*num_xg+ 8*ixg+1].W= "Vd";
	data_MEX_comb[ijack][8*im*num_xg +8*ixg+1].imass = im;
	data_MEX_comb[ijack][8*im*num_xg +8*ixg+1].ixg = ixg;
	data_MEX_comb[ijack][8*im*num_xg +8*ixg+1].Z = ZT_Mh[im];


	data_MEX_comb[ijack][8*im*num_xg+8*ixg+2].FF = (FF_M_list_new[1][im].distr_list[ixg]/fps.distr_list[im]).distr[ijack];
	data_MEX_comb[ijack][8*im*num_xg+8*ixg+2].FF_err= (FF_M_list_new[1][im].distr_list[ixg]/fps.distr_list[im]).err();
	data_MEX_comb[ijack][8*im*num_xg+8*ixg+2].xg= Bs_xg_t_list[ixg];
	data_MEX_comb[ijack][8*im*num_xg+8*ixg+2].M= masses.distr_list[im].distr[ijack];
	data_MEX_comb[ijack][8*im*num_xg+ 8*ixg+2].W= "Au";
	data_MEX_comb[ijack][8*im*num_xg +8*ixg+2].imass = im;
	data_MEX_comb[ijack][8*im*num_xg +8*ixg+2].ixg = ixg;
	data_MEX_comb[ijack][8*im*num_xg +8*ixg+2].Z = ZT_Mh[im];


	data_MEX_comb[ijack][8*im*num_xg+8*ixg+3].FF = (FF_M_list_new[4][im].distr_list[ixg]/fps.distr_list[im]).distr[ijack];
	data_MEX_comb[ijack][8*im*num_xg+8*ixg+3].FF_err= (FF_M_list_new[4][im].distr_list[ixg]/fps.distr_list[im]).err();
	data_MEX_comb[ijack][8*im*num_xg+8*ixg+3].xg= Bs_xg_t_list[ixg];
	data_MEX_comb[ijack][8*im*num_xg+8*ixg+3].M= masses.distr_list[im].distr[ijack];
	data_MEX_comb[ijack][8*im*num_xg+ 8*ixg+3].W= "Vu";
	data_MEX_comb[ijack][8*im*num_xg +8*ixg+3].imass = im;
	data_MEX_comb[ijack][8*im*num_xg +8*ixg+3].ixg = ixg;
	data_MEX_comb[ijack][8*im*num_xg +8*ixg+3].Z = ZT_Mh[im];



	
	data_MEX_comb[ijack][8*im*num_xg+ 8*ixg+4].FF = (FF_M_list_new[id_TA][im].distr_list[ixg]/fps.distr_list[im]).distr[ijack];
	data_MEX_comb[ijack][8*im*num_xg+ 8*ixg+4].FF_err= (FF_M_list_new[id_TA][im].distr_list[ixg]/fps.distr_list[im]).err();
	data_MEX_comb[ijack][8*im*num_xg+ 8*ixg+4].xg= Bs_xg_t_list[ixg];
	data_MEX_comb[ijack][8*im*num_xg+ 8*ixg+4].M= masses.distr_list[im].distr[ijack];
	data_MEX_comb[ijack][8*im*num_xg+ 8*ixg+4].W= "TAd";
	data_MEX_comb[ijack][8*im*num_xg +8*ixg+4].imass = im;
	data_MEX_comb[ijack][8*im*num_xg +8*ixg+4].ixg = ixg;
	data_MEX_comb[ijack][8*im*num_xg +8*ixg+4].Z = ZT_Mh[im];

	data_MEX_comb[ijack][8*im*num_xg+8*ixg+5].FF = (FF_M_list_new[id_TV][im].distr_list[ixg]/fps.distr_list[im]).distr[ijack];
	data_MEX_comb[ijack][8*im*num_xg+8*ixg+5].FF_err= (FF_M_list_new[id_TV][im].distr_list[ixg]/fps.distr_list[im]).err();
	data_MEX_comb[ijack][8*im*num_xg+8*ixg+5].xg= Bs_xg_t_list[ixg];
	data_MEX_comb[ijack][8*im*num_xg+8*ixg+5].M= masses.distr_list[im].distr[ijack];
	data_MEX_comb[ijack][8*im*num_xg+ 8*ixg+5].W= "TVd";
	data_MEX_comb[ijack][8*im*num_xg +8*ixg+5].imass = im;
	data_MEX_comb[ijack][8*im*num_xg +8*ixg+5].ixg = ixg;
	data_MEX_comb[ijack][8*im*num_xg +8*ixg+5].Z = ZT_Mh[im];

	data_MEX_comb[ijack][8*im*num_xg+8*ixg+6].FF = (FF_M_list_new[id_TA-1][im].distr_list[ixg]/fps.distr_list[im]).distr[ijack];
	data_MEX_comb[ijack][8*im*num_xg+8*ixg+6].FF_err= (FF_M_list_new[id_TA-1][im].distr_list[ixg]/fps.distr_list[im]).err();
	data_MEX_comb[ijack][8*im*num_xg+8*ixg+6].xg= Bs_xg_t_list[ixg];
	data_MEX_comb[ijack][8*im*num_xg+8*ixg+6].M= masses.distr_list[im].distr[ijack];
	data_MEX_comb[ijack][8*im*num_xg+ 8*ixg+6].W= "TAu";
	data_MEX_comb[ijack][8*im*num_xg +8*ixg+6].imass = im;
	data_MEX_comb[ijack][8*im*num_xg +8*ixg+6].ixg = ixg;
	data_MEX_comb[ijack][8*im*num_xg +8*ixg+6].Z = ZT_Mh[im];


	data_MEX_comb[ijack][8*im*num_xg+8*ixg+7].FF = (FF_M_list_new[id_TV-1][im].distr_list[ixg]/fps.distr_list[im]).distr[ijack];
	data_MEX_comb[ijack][8*im*num_xg+8*ixg+7].FF_err= (FF_M_list_new[id_TV-1][im].distr_list[ixg]/fps.distr_list[im]).err();
	data_MEX_comb[ijack][8*im*num_xg+8*ixg+7].xg= Bs_xg_t_list[ixg];
	data_MEX_comb[ijack][8*im*num_xg+8*ixg+7].M= masses.distr_list[im].distr[ijack];
	data_MEX_comb[ijack][8*im*num_xg+ 8*ixg+7].W= "TVu";
	data_MEX_comb[ijack][8*im*num_xg +8*ixg+7].imass = im;
	data_MEX_comb[ijack][8*im*num_xg +8*ixg+7].ixg = ixg;
	data_MEX_comb[ijack][8*im*num_xg +8*ixg+7].Z = ZT_Mh[im];

	
      if(ijack==0) {  
	data_MEX_comb_ch2[ijack][8*im*num_xg+ 8*ixg].FF = (FF_M_list_new[2][im].distr_list[ixg]/fps.distr_list[im]).ave();
	data_MEX_comb_ch2[ijack][8*im*num_xg+ 8*ixg].FF_err= (FF_M_list_new[2][im].distr_list[ixg]/fps.distr_list[im]).err();
	data_MEX_comb_ch2[ijack][8*im*num_xg+ 8*ixg].xg= Bs_xg_t_list[ixg];
	data_MEX_comb_ch2[ijack][8*im*num_xg+ 8*ixg].M= masses.distr_list[im].ave();
	data_MEX_comb_ch2[ijack][8*im*num_xg+ 8*ixg].W= "Ad";
	data_MEX_comb_ch2[ijack][8*im*num_xg +8*ixg].imass = im;
	data_MEX_comb_ch2[ijack][8*im*num_xg +8*ixg].ixg = ixg;
	data_MEX_comb_ch2[ijack][8*im*num_xg +8*ixg].Z = ZT_Mh[im];

	data_MEX_comb_ch2[ijack][8*im*num_xg+8*ixg+1].FF = (FF_M_list_new[5][im].distr_list[ixg]/fps.distr_list[im]).ave();
	data_MEX_comb_ch2[ijack][8*im*num_xg+8*ixg+1].FF_err= (FF_M_list_new[5][im].distr_list[ixg]/fps.distr_list[im]).err();
	data_MEX_comb_ch2[ijack][8*im*num_xg+8*ixg+1].xg= Bs_xg_t_list[ixg];
	data_MEX_comb_ch2[ijack][8*im*num_xg+8*ixg+1].M= masses.distr_list[im].ave();
	data_MEX_comb_ch2[ijack][8*im*num_xg+ 8*ixg+1].W= "Vd";
	data_MEX_comb_ch2[ijack][8*im*num_xg +8*ixg+1].imass = im;
	data_MEX_comb_ch2[ijack][8*im*num_xg +8*ixg+1].ixg = ixg;
	data_MEX_comb_ch2[ijack][8*im*num_xg +8*ixg+1].Z = ZT_Mh[im];
		

	data_MEX_comb_ch2[ijack][8*im*num_xg+8*ixg+2].FF = (FF_M_list_new[1][im].distr_list[ixg]/fps.distr_list[im]).ave();
	data_MEX_comb_ch2[ijack][8*im*num_xg+8*ixg+2].FF_err= (FF_M_list_new[1][im].distr_list[ixg]/fps.distr_list[im]).err();
	data_MEX_comb_ch2[ijack][8*im*num_xg+8*ixg+2].xg= Bs_xg_t_list[ixg];
	data_MEX_comb_ch2[ijack][8*im*num_xg+8*ixg+2].M= masses.distr_list[im].ave();
	data_MEX_comb_ch2[ijack][8*im*num_xg+ 8*ixg+2].W= "Au";
	data_MEX_comb_ch2[ijack][8*im*num_xg +8*ixg+2].imass = im;
	data_MEX_comb_ch2[ijack][8*im*num_xg +8*ixg+2].ixg = ixg;
	data_MEX_comb_ch2[ijack][8*im*num_xg +8*ixg+2].Z = ZT_Mh[im];
       

	data_MEX_comb_ch2[ijack][8*im*num_xg+8*ixg+3].FF = (FF_M_list_new[4][im].distr_list[ixg]/fps.distr_list[im]).ave();
	data_MEX_comb_ch2[ijack][8*im*num_xg+8*ixg+3].FF_err= (FF_M_list_new[4][im].distr_list[ixg]/fps.distr_list[im]).err();
	data_MEX_comb_ch2[ijack][8*im*num_xg+8*ixg+3].xg= Bs_xg_t_list[ixg];
	data_MEX_comb_ch2[ijack][8*im*num_xg+8*ixg+3].M= masses.distr_list[im].ave();
	data_MEX_comb_ch2[ijack][8*im*num_xg+ 8*ixg+3].W= "Vu";
	data_MEX_comb_ch2[ijack][8*im*num_xg +8*ixg+3].imass = im;
	data_MEX_comb_ch2[ijack][8*im*num_xg +8*ixg+3].ixg = ixg;
	data_MEX_comb_ch2[ijack][8*im*num_xg +8*ixg+3].Z = ZT_Mh[im];

	
	data_MEX_comb_ch2[ijack][8*im*num_xg+ 8*ixg+4].FF = (FF_M_list_new[id_TA][im].distr_list[ixg]/fps.distr_list[im]).ave();
	data_MEX_comb_ch2[ijack][8*im*num_xg+ 8*ixg+4].FF_err= (FF_M_list_new[id_TA][im].distr_list[ixg]/fps.distr_list[im]).err();
	data_MEX_comb_ch2[ijack][8*im*num_xg+ 8*ixg+4].xg= Bs_xg_t_list[ixg];
	data_MEX_comb_ch2[ijack][8*im*num_xg+ 8*ixg+4].M= masses.distr_list[im].ave();
	data_MEX_comb_ch2[ijack][8*im*num_xg+ 8*ixg+4].W= "TAd";
	data_MEX_comb_ch2[ijack][8*im*num_xg +8*ixg+4].imass = im;
	data_MEX_comb_ch2[ijack][8*im*num_xg +8*ixg+4].ixg = ixg;
	data_MEX_comb_ch2[ijack][8*im*num_xg +8*ixg+4].Z = ZT_Mh[im];

	data_MEX_comb_ch2[ijack][8*im*num_xg+8*ixg+5].FF = (FF_M_list_new[id_TV][im].distr_list[ixg]/fps.distr_list[im]).ave();
	data_MEX_comb_ch2[ijack][8*im*num_xg+8*ixg+5].FF_err= (FF_M_list_new[id_TV][im].distr_list[ixg]/fps.distr_list[im]).err();
	data_MEX_comb_ch2[ijack][8*im*num_xg+8*ixg+5].xg= Bs_xg_t_list[ixg];
	data_MEX_comb_ch2[ijack][8*im*num_xg+8*ixg+5].M= masses.distr_list[im].ave();
	data_MEX_comb_ch2[ijack][8*im*num_xg+ 8*ixg+5].W= "TVd";
	data_MEX_comb_ch2[ijack][8*im*num_xg +8*ixg+5].imass = im;
	data_MEX_comb_ch2[ijack][8*im*num_xg +8*ixg+5].ixg = ixg;
	data_MEX_comb_ch2[ijack][8*im*num_xg +8*ixg+5].Z = ZT_Mh[im];


	data_MEX_comb_ch2[ijack][8*im*num_xg+8*ixg+6].FF = (FF_M_list_new[id_TA-1][im].distr_list[ixg]/fps.distr_list[im]).ave();
	data_MEX_comb_ch2[ijack][8*im*num_xg+8*ixg+6].FF_err= (FF_M_list_new[id_TA-1][im].distr_list[ixg]/fps.distr_list[im]).err();
	data_MEX_comb_ch2[ijack][8*im*num_xg+8*ixg+6].xg= Bs_xg_t_list[ixg];
	data_MEX_comb_ch2[ijack][8*im*num_xg+8*ixg+6].M= masses.distr_list[im].ave();
	data_MEX_comb_ch2[ijack][8*im*num_xg+ 8*ixg+6].W= "TAu";
	data_MEX_comb_ch2[ijack][8*im*num_xg +8*ixg+6].imass = im;
	data_MEX_comb_ch2[ijack][8*im*num_xg +8*ixg+6].ixg = ixg;
	data_MEX_comb_ch2[ijack][8*im*num_xg +8*ixg+6].Z = ZT_Mh[im];


	data_MEX_comb_ch2[ijack][8*im*num_xg+8*ixg+7].FF = (FF_M_list_new[id_TV-1][im].distr_list[ixg]/fps.distr_list[im]).ave();
	data_MEX_comb_ch2[ijack][8*im*num_xg+8*ixg+7].FF_err= (FF_M_list_new[id_TV-1][im].distr_list[ixg]/fps.distr_list[im]).err();
	data_MEX_comb_ch2[ijack][8*im*num_xg+8*ixg+7].xg= Bs_xg_t_list[ixg];
	data_MEX_comb_ch2[ijack][8*im*num_xg+8*ixg+7].M= masses.distr_list[im].ave();
	data_MEX_comb_ch2[ijack][8*im*num_xg+ 8*ixg+7].W= "TVu";
	data_MEX_comb_ch2[ijack][8*im*num_xg +8*ixg+7].imass = im;
	data_MEX_comb_ch2[ijack][8*im*num_xg +8*ixg+7].ixg = ixg;
	data_MEX_comb_ch2[ijack][8*im*num_xg +8*ixg+7].Z = ZT_Mh[im];
		
      }
      }
    }
  }

  cout<<"####################################################"<<endl;

  //for (auto &d: data_MEX_comb_ch2[0]) cout<<"ch: "<<d.W<<"xg: "<<d.xg<<" FF: "<<d.FF<<" +- "<<d.FF_err<<endl;
 
  cout<<"####################################################"<<endl;

  int Nmeas=Nmasses*num_xg*4*((separate_ud==0)?1:2);
  int Npars=bf_MEX_comb.Get_number_of_fit_pars(); //number of parameters

  //cout<<"Nmeas: "<<Nmeas<<endl;
  //cout<<"Npars: "<<Npars<<endl;


  //append
  bf_MEX_comb.Append_to_input_par(data_MEX_comb);
  bf_MEX_comb_ch2.Append_to_input_par(data_MEX_comb_ch2);
  //append prior pars
  it_prior= Fix_fit_pars_list[ifit].find("M3");
  if(it_prior != Fix_fit_pars_list[ifit].end()) {
    if(it_prior->second == 1) {
      for(int ij=0;ij< NJ;ij++) 
	{bf_MEX_comb.Append_to_prior("M3", 0.0, 1); bf_MEX_comb_ch2.Append_to_prior("M3", 0.0,1);} //0.3
    }
  }
  
  it_prior= Fix_fit_pars_list[ifit].find("M4");
  if(it_prior != Fix_fit_pars_list[ifit].end()) {
    if(it_prior->second == 1) {
      for(int ij=0;ij< NJ;ij++) 
	{bf_MEX_comb.Append_to_prior("M4", -1.5, 1); bf_MEX_comb_ch2.Append_to_prior("M4", -1.5,1);}
    }
  }

  it_prior= Fix_fit_pars_list[ifit].find("TVE1");
  if(it_prior != Fix_fit_pars_list[ifit].end()) {
    if(it_prior->second == 1) {
      for(int ij=0;ij< NJ;ij++) 
	{bf_MEX_comb.Append_to_prior("TVE1", 0.0, 1); bf_MEX_comb_ch2.Append_to_prior("TVE1", 0.0,1);}
    }
  }
  
   
  it_prior= Fix_fit_pars_list[ifit].find("A1");
  if(it_prior != Fix_fit_pars_list[ifit].end()) {
    if(it_prior->second == 1) {
      for(int ij=0;ij< NJ;ij++) 
	{bf_MEX_comb.Append_to_prior("A1", 0.0, 0.4); bf_MEX_comb_ch2.Append_to_prior("A1", 0.0,0.4);}
    }
  }
  it_prior= Fix_fit_pars_list[ifit].find("A2");
  if(it_prior != Fix_fit_pars_list[ifit].end()) {
    if(it_prior->second == 1) {
      for(int ij=0;ij< NJ;ij++) 
	{bf_MEX_comb.Append_to_prior("A2", 0.0, 0.4); bf_MEX_comb_ch2.Append_to_prior("A2", 0.0,0.4);}
    }
  }    
  for( auto &mpar : M2_pars) {
    it_prior= Fix_fit_pars_list[ifit].find(mpar);
    if(it_prior != Fix_fit_pars_list[ifit].end()) {
      if(it_prior->second == 1) {
	for(int ij=0;ij< NJ;ij++) {
	  bf_MEX_comb.Append_to_prior(mpar, 0.0, 5.0); bf_MEX_comb_ch2.Append_to_prior(mpar, 0.0, 5.0); }
	}
    }
  }
  

  
  //fit
  Bt_fit_MEX_comb= bf_MEX_comb.Perform_bootstrap_fit();
  Bt_fit_MEX_comb_ch2= bf_MEX_comb_ch2.Perform_bootstrap_fit();

    
  double ch2_red_MEX_comb= Bt_fit_MEX_comb_ch2.get_ch2_ave()/( Nmeas -Npars);

 
  Ndofs.push_back( Nmeas-Npars);
  Nmeasurements.push_back( Nmeas);
  Ch2s.push_back( Bt_fit_MEX_comb_ch2.get_ch2_ave());


 
 	
  //retrieve params
  distr_t A(1),  R1(1), R2(1), DV(1), DA(1), M1(1), M2(1), M3(1), M4(1),  SV1(1), SV2(1), SV3(1), SA1(1), SA2(1), SA3(1), BV1(1), BV2(1), BA1(1), BA2(1);
  distr_t A1(1), A2(1), RTV1(1), RTV2(1), RTA1(1), RTA2(1),  TVM1(1), TVE1(1), TAM1(1), TAE1(1), STV1(1), STV2(1), STV3(1), STA1(1), STA2(1), STA3(1), BTV1(1), BTV2(1), BTA1(1), BTA2(1); 
  for(int ijack=0;ijack<NJ;ijack++) {
    
    A.distr.push_back( Bt_fit_MEX_comb.par[ijack].A);
    R1.distr.push_back( Bt_fit_MEX_comb.par[ijack].R1);
    R2.distr.push_back( Bt_fit_MEX_comb.par[ijack].R2);

    M1.distr.push_back( Bt_fit_MEX_comb.par[ijack].M1);
    M2.distr.push_back( Bt_fit_MEX_comb.par[ijack].M2);
    M3.distr.push_back( Bt_fit_MEX_comb.par[ijack].M3);
    M4.distr.push_back( Bt_fit_MEX_comb.par[ijack].M4);
    
    DV.distr.push_back( Bt_fit_MEX_comb.par[ijack].DV);
    DA.distr.push_back( Bt_fit_MEX_comb.par[ijack].DA);

    SV1.distr.push_back( Bt_fit_MEX_comb.par[ijack].SV1);
    SV2.distr.push_back( Bt_fit_MEX_comb.par[ijack].SV2);
    SV3.distr.push_back( Bt_fit_MEX_comb.par[ijack].SV3);

    SA1.distr.push_back( Bt_fit_MEX_comb.par[ijack].SA1);
    SA2.distr.push_back( Bt_fit_MEX_comb.par[ijack].SA2);
    SA3.distr.push_back( Bt_fit_MEX_comb.par[ijack].SA3);

    BV1.distr.push_back( Bt_fit_MEX_comb.par[ijack].BV1);
    BV2.distr.push_back( Bt_fit_MEX_comb.par[ijack].BV2);

    BA1.distr.push_back( Bt_fit_MEX_comb.par[ijack].BA1);
    BA2.distr.push_back( Bt_fit_MEX_comb.par[ijack].BA2);



    A1.distr.push_back( Bt_fit_MEX_comb.par[ijack].A1);
    A2.distr.push_back( Bt_fit_MEX_comb.par[ijack].A2);
    RTV1.distr.push_back( Bt_fit_MEX_comb.par[ijack].RTV1);
    RTV2.distr.push_back( Bt_fit_MEX_comb.par[ijack].RTV2);
    RTA1.distr.push_back( Bt_fit_MEX_comb.par[ijack].RTA1);
    RTA2.distr.push_back( Bt_fit_MEX_comb.par[ijack].RTA2);

    TVM1.distr.push_back( Bt_fit_MEX_comb.par[ijack].TVM1);
    TVE1.distr.push_back( Bt_fit_MEX_comb.par[ijack].TVE1);
    TAM1.distr.push_back( Bt_fit_MEX_comb.par[ijack].TAM1);
    TAE1.distr.push_back( Bt_fit_MEX_comb.par[ijack].TAE1);

    STV1.distr.push_back( Bt_fit_MEX_comb.par[ijack].STV1);
    STV2.distr.push_back( Bt_fit_MEX_comb.par[ijack].STV2);
    STV3.distr.push_back( Bt_fit_MEX_comb.par[ijack].STV3);

    STA1.distr.push_back( Bt_fit_MEX_comb.par[ijack].STA1);
    STA2.distr.push_back( Bt_fit_MEX_comb.par[ijack].STA2);
    STA3.distr.push_back( Bt_fit_MEX_comb.par[ijack].STA3);

    BTV1.distr.push_back( Bt_fit_MEX_comb.par[ijack].BTV1);
    BTV2.distr.push_back( Bt_fit_MEX_comb.par[ijack].BTV2);

    BTA1.distr.push_back( Bt_fit_MEX_comb.par[ijack].BTA1);
    BTA2.distr.push_back( Bt_fit_MEX_comb.par[ijack].BTA2);
    

    
  }

 
  //push_back important fit parameters
  CA_list.push_back( DA);
  CV_list.push_back( DV);
  CA_T_list.push_back( DA+A2);
  R_list.push_back( A);
  RT_list.push_back( A + A1);
  M3_list.push_back( M3+1.0);
  M4_list.push_back( M4);
  M1_list.push_back(M1);
  M2_list.push_back(M2);
  M1_T_list.push_back(TVM1);
  M2_T_list.push_back(TAE1);
  d_T_list.push_back(TVE1+1.0);
 

  //print fit function
  distr_t_list F_MEX_fit(1);
	
  //### Determine FF at Bs meson mass
  ipar_MEX_comb pp_Bs_comb;
  pp_Bs_comb.M= MBs;
  pp_Bs_comb.Z= ZT_fit_Mb;
 
  //scan over xg
 
  distr_t_list FA_d_Bs_xg(1,xg_scan.size(),NJ);
  distr_t_list FV_d_Bs_xg(1,xg_scan.size(), NJ);
  distr_t_list FA_u_Bs_xg(1,xg_scan.size(),NJ);
  distr_t_list FV_u_Bs_xg(1,xg_scan.size(), NJ);

  distr_t_list FTA_d_Bs_xg(1,xg_scan.size(),NJ);
  distr_t_list FTV_d_Bs_xg(1,xg_scan.size(), NJ);
  distr_t_list FTA_u_Bs_xg(1,xg_scan.size(),NJ);
  distr_t_list FTV_u_Bs_xg(1,xg_scan.size(), NJ);


  for(int ixg=0;ixg<(signed)xg_scan.size();ixg++) {
    pp_Bs_comb.xg= xg_scan[ixg];
    for(int ijack=0;ijack<NJ;ijack++) {
      
      pp_Bs_comb.W="Ad";
      FA_d_Bs_xg.distr_list[ixg].distr[ijack] =  bf_MEX_comb.ansatz( Bt_fit_MEX_comb.par[ijack], pp_Bs_comb);
      pp_Bs_comb.W="Vd";
      FV_d_Bs_xg.distr_list[ixg].distr[ijack] =  bf_MEX_comb.ansatz( Bt_fit_MEX_comb.par[ijack], pp_Bs_comb);
      pp_Bs_comb.W="Au";
      FA_u_Bs_xg.distr_list[ixg].distr[ijack] =  bf_MEX_comb.ansatz( Bt_fit_MEX_comb.par[ijack], pp_Bs_comb);
      pp_Bs_comb.W="Vu";
      FV_u_Bs_xg.distr_list[ixg].distr[ijack] =  bf_MEX_comb.ansatz( Bt_fit_MEX_comb.par[ijack], pp_Bs_comb);

      pp_Bs_comb.W="TAd";
      FTA_d_Bs_xg.distr_list[ixg].distr[ijack] =  bf_MEX_comb.ansatz( Bt_fit_MEX_comb.par[ijack], pp_Bs_comb);
      pp_Bs_comb.W="TVd";
      FTV_d_Bs_xg.distr_list[ixg].distr[ijack] =  bf_MEX_comb.ansatz( Bt_fit_MEX_comb.par[ijack], pp_Bs_comb);
      pp_Bs_comb.W="TAu";
      FTA_u_Bs_xg.distr_list[ixg].distr[ijack] =  bf_MEX_comb.ansatz( Bt_fit_MEX_comb.par[ijack], pp_Bs_comb);
      pp_Bs_comb.W="TVu";
      FTV_u_Bs_xg.distr_list[ixg].distr[ijack] =  bf_MEX_comb.ansatz( Bt_fit_MEX_comb.par[ijack], pp_Bs_comb);      
      
    }

    //push_back the results
    FA_d_extr[ixg][ifit] = FA_d_Bs_xg.distr_list[ixg];
    FA_u_extr[ixg][ifit] = FA_u_Bs_xg.distr_list[ixg];
    FA_extr[ixg][ifit] = FA_d_Bs_xg.distr_list[ixg] + FA_u_Bs_xg.distr_list[ixg];
    
    FV_d_extr[ixg][ifit] = FV_d_Bs_xg.distr_list[ixg];
    FV_u_extr[ixg][ifit] = FV_u_Bs_xg.distr_list[ixg];
    FV_extr[ixg][ifit] = FV_d_Bs_xg.distr_list[ixg] + FV_u_Bs_xg.distr_list[ixg];
    
    FTA_d_extr[ixg][ifit] = FTA_d_Bs_xg.distr_list[ixg];
    FTA_u_extr[ixg][ifit] = FTA_u_Bs_xg.distr_list[ixg];
    FTA_extr[ixg][ifit] = FTA_d_Bs_xg.distr_list[ixg] + FTA_u_Bs_xg.distr_list[ixg];
    
    FTV_d_extr[ixg][ifit] = FTV_d_Bs_xg.distr_list[ixg];
    FTV_u_extr[ixg][ifit] = FTV_u_Bs_xg.distr_list[ixg];
    FTV_extr[ixg][ifit] = FTV_d_Bs_xg.distr_list[ixg] + FTV_u_Bs_xg.distr_list[ixg];


        
  }

  
  
  
  
  
  //#################################
    
  Print_To_File({}, {xg_scan, (FA_d_Bs_xg*f_Bs).ave(), (FA_d_Bs_xg*f_Bs).err(),  FA_d_Bs_xg.ave(), FA_d_Bs_xg.err()},  "../data/ph_emission/"+ph_type+"/Bs_extr/comb_fits/FA_d_comb_Bs_fit_"+to_string(ifit)+".dat", "", "");
  Print_To_File({}, {xg_scan, (FV_d_Bs_xg*f_Bs).ave(), (FV_d_Bs_xg*f_Bs).err(), FV_d_Bs_xg.ave(), FV_d_Bs_xg.err()},  "../data/ph_emission/"+ph_type+"/Bs_extr/comb_fits/FV_d_comb_Bs_fit_"+to_string(ifit)+".dat", "", "");
  Print_To_File({}, {xg_scan, (FA_u_Bs_xg*f_Bs).ave(), (FA_u_Bs_xg*f_Bs).err(),  FA_u_Bs_xg.ave(), FA_u_Bs_xg.err()},  "../data/ph_emission/"+ph_type+"/Bs_extr/comb_fits/FA_u_comb_Bs_fit_"+to_string(ifit)+".dat", "", "");
  Print_To_File({}, {xg_scan, (FV_u_Bs_xg*f_Bs).ave(), (FV_u_Bs_xg*f_Bs).err(), FV_u_Bs_xg.ave(), FV_u_Bs_xg.err()},  "../data/ph_emission/"+ph_type+"/Bs_extr/comb_fits/FV_u_comb_Bs_fit_"+to_string(ifit)+".dat", "", "");
  
  Print_To_File({}, {xg_scan, (FTA_d_Bs_xg*f_Bs).ave(), (FTA_d_Bs_xg*f_Bs).err(),  FTA_d_Bs_xg.ave(), FTA_d_Bs_xg.err()},  "../data/ph_emission/"+ph_type+"/Bs_extr/comb_fits/FA_T_d_comb_Bs_fit_"+to_string(ifit)+".dat", "", "");
  Print_To_File({}, {xg_scan, (FTV_d_Bs_xg*f_Bs).ave(), (FTV_d_Bs_xg*f_Bs).err(), FTV_d_Bs_xg.ave(), FTV_d_Bs_xg.err()},  "../data/ph_emission/"+ph_type+"/Bs_extr/comb_fits/FV_T_d_comb_Bs_fit_"+to_string(ifit)+".dat", "", "");
  Print_To_File({}, {xg_scan, (FTA_u_Bs_xg*f_Bs).ave(), (FTA_u_Bs_xg*f_Bs).err(),  FTA_u_Bs_xg.ave(), FTA_u_Bs_xg.err()},  "../data/ph_emission/"+ph_type+"/Bs_extr/comb_fits/FA_T_u_comb_Bs_fit_"+to_string(ifit)+".dat", "", "");
  Print_To_File({}, {xg_scan, (FTV_u_Bs_xg*f_Bs).ave(), (FTV_u_Bs_xg*f_Bs).err(), FTV_u_Bs_xg.ave(), FTV_u_Bs_xg.err()},  "../data/ph_emission/"+ph_type+"/Bs_extr/comb_fits/FV_T_u_comb_Bs_fit_"+to_string(ifit)+".dat", "", "");
  
  
  Print_To_File({}, {xg_scan, ((FA_d_Bs_xg+FA_u_Bs_xg)*f_Bs).ave(), ((FA_d_Bs_xg+FA_u_Bs_xg)*f_Bs).err(),  (FA_d_Bs_xg+FA_u_Bs_xg).ave(), (FA_d_Bs_xg+FA_u_Bs_xg).err()},  "../data/ph_emission/"+ph_type+"/Bs_extr/comb_fits/FA_comb_Bs_fit_"+to_string(ifit)+".dat", "", "");
  Print_To_File({}, {xg_scan, ((FV_d_Bs_xg+FV_u_Bs_xg)*f_Bs).ave(), ((FV_d_Bs_xg+FV_u_Bs_xg)*f_Bs).err(), (FV_d_Bs_xg+FV_u_Bs_xg).ave(), (FV_d_Bs_xg+FV_u_Bs_xg).err()},  "../data/ph_emission/"+ph_type+"/Bs_extr/comb_fits/FV_comb_Bs_fit_"+to_string(ifit)+".dat", "", "");
  
  Print_To_File({}, {xg_scan, ((FTA_d_Bs_xg+FTA_u_Bs_xg)*f_Bs).ave(), ((FTA_d_Bs_xg+FTA_u_Bs_xg)*f_Bs).err(),  (FTA_d_Bs_xg+FTA_u_Bs_xg).ave(), (FTA_d_Bs_xg+FTA_u_Bs_xg).err()},  "../data/ph_emission/"+ph_type+"/Bs_extr/comb_fits/FA_T_comb_Bs_fit_"+to_string(ifit)+".dat", "", "");
  Print_To_File({}, {xg_scan, ((FTV_d_Bs_xg+FTV_u_Bs_xg)*f_Bs).ave(), ((FTV_d_Bs_xg+FTV_u_Bs_xg)*f_Bs).err(), (FTV_d_Bs_xg+FTV_u_Bs_xg).ave(), (FTV_d_Bs_xg+FTV_u_Bs_xg).err()},  "../data/ph_emission/"+ph_type+"/Bs_extr/comb_fits/FV_T_comb_Bs_fit_"+to_string(ifit)+".dat", "", "");

  
  
  


  //scan over masses and xg

  

  for(int ixg=0;ixg<num_xg;ixg++) {

    distr_t_list FA_d(1,Inv_Masses_to_print_MEX.size(),NJ), FV_d(1, Inv_Masses_to_print_MEX.size(), NJ);
    distr_t_list FA_u(1,Inv_Masses_to_print_MEX.size(),NJ), FV_u(1, Inv_Masses_to_print_MEX.size(), NJ);
    distr_t_list FA(1,Inv_Masses_to_print_MEX.size(),NJ), FV(1, Inv_Masses_to_print_MEX.size(), NJ);

    distr_t_list FTA_d(1,Inv_Masses_to_print_MEX.size(),NJ), FTV_d(1, Inv_Masses_to_print_MEX.size(), NJ);
    distr_t_list FTA_u(1,Inv_Masses_to_print_MEX.size(),NJ), FTV_u(1, Inv_Masses_to_print_MEX.size(), NJ);
    distr_t_list FTA(1,Inv_Masses_to_print_MEX.size(),NJ), FTV(1, Inv_Masses_to_print_MEX.size(), NJ);
    

    

    
    for(int im=0;im<(signed)Inv_Masses_to_print_MEX.size(); im++) {

      ipar_MEX_comb pp;
      pp.M = 1.0/Inv_Masses_to_print_MEX[im];
      pp.Z= ZT_Mh_to_print[im];
      pp.xg= Bs_xg_t_list[ixg];
      pp.ixg= ixg;
      
      for(int ijack=0;ijack<NJ;ijack++) {
	pp.W="Ad";
	FA_d.distr_list[im].distr[ijack] = bf_MEX_comb.ansatz( Bt_fit_MEX_comb.par[ijack], pp);
	pp.W="Vd";
	FV_d.distr_list[im].distr[ijack] = bf_MEX_comb.ansatz( Bt_fit_MEX_comb.par[ijack], pp);
	pp.W="Au";
	FA_u.distr_list[im].distr[ijack] = bf_MEX_comb.ansatz( Bt_fit_MEX_comb.par[ijack], pp);
	pp.W="Vu";
	FV_u.distr_list[im].distr[ijack] = bf_MEX_comb.ansatz( Bt_fit_MEX_comb.par[ijack], pp);

	pp.W="TAd";
	FTA_d.distr_list[im].distr[ijack] = bf_MEX_comb.ansatz( Bt_fit_MEX_comb.par[ijack], pp);
	pp.W="TVd";
	FTV_d.distr_list[im].distr[ijack] = bf_MEX_comb.ansatz( Bt_fit_MEX_comb.par[ijack], pp);
	pp.W="TAu";
	FTA_u.distr_list[im].distr[ijack] = bf_MEX_comb.ansatz( Bt_fit_MEX_comb.par[ijack], pp);
	pp.W="TVu";
	FTV_u.distr_list[im].distr[ijack] = bf_MEX_comb.ansatz( Bt_fit_MEX_comb.par[ijack], pp);

	FA.distr_list[im].distr[ijack] = FA_d.distr_list[im].distr[ijack] + FA_u.distr_list[im].distr[ijack];
	FV.distr_list[im].distr[ijack] = FV_d.distr_list[im].distr[ijack] + FV_u.distr_list[im].distr[ijack];
	FTA.distr_list[im].distr[ijack] = FTA_d.distr_list[im].distr[ijack] + FTA_u.distr_list[im].distr[ijack];
	FTV.distr_list[im].distr[ijack] = FTV_d.distr_list[im].distr[ijack] + FTV_u.distr_list[im].distr[ijack];
	
      }


      FA_d_sims_m_dep[ixg][im].push_back( FA_d.distr_list[im]);
      FA_u_sims_m_dep[ixg][im].push_back( FA_u.distr_list[im]);
      FA_sims_m_dep[ixg][im].push_back(FA.distr_list[im]);

      FV_d_sims_m_dep[ixg][im].push_back(FV_d.distr_list[im]);
      FV_u_sims_m_dep[ixg][im].push_back(FV_u.distr_list[im]);
      FV_sims_m_dep[ixg][im].push_back(FV.distr_list[im]);
      
      FTA_d_sims_m_dep[ixg][im].push_back(FTA_d.distr_list[im]);
      FTA_u_sims_m_dep[ixg][im].push_back(FTA_u.distr_list[im]);
      FTA_sims_m_dep[ixg][im].push_back(FTA.distr_list[im]);
      
      FTV_d_sims_m_dep[ixg][im].push_back(FTV_d.distr_list[im]);
      FTV_u_sims_m_dep[ixg][im].push_back(FTV_u.distr_list[im]);
      FTV_sims_m_dep[ixg][im].push_back(FTV.distr_list[im]);

     
      
    }

  

    
  
    Print_To_File({}, {Inv_Masses_to_print_MEX, FA_d.ave(), FA_d.err()}, "../data/ph_emission/"+ph_type+"/Bs_extr/comb_fits/FA_d_comb_ixg_"+to_string(ixg+1)+"_fit_"+to_string(ifit)+".fit_func", "", "ch2/dof: "+to_string_with_precision(ch2_red_MEX_comb,5));
    Print_To_File({}, {Inv_Masses_to_print_MEX, FV_d.ave(), FV_d.err()}, "../data/ph_emission/"+ph_type+"/Bs_extr/comb_fits/FV_d_comb_ixg_"+to_string(ixg+1)+"_fit_"+to_string(ifit)+".fit_func", "", "ch2/dof: "+to_string_with_precision(ch2_red_MEX_comb,5));
    Print_To_File({}, {Inv_Masses_to_print_MEX, FA_u.ave(), FA_u.err()}, "../data/ph_emission/"+ph_type+"/Bs_extr/comb_fits/FA_u_comb_ixg_"+to_string(ixg+1)+"_fit_"+to_string(ifit)+".fit_func", "", "ch2/dof: "+to_string_with_precision(ch2_red_MEX_comb,5));
    Print_To_File({}, {Inv_Masses_to_print_MEX, FV_u.ave(), FV_u.err()}, "../data/ph_emission/"+ph_type+"/Bs_extr/comb_fits/FV_u_comb_ixg_"+to_string(ixg+1)+"_fit_"+to_string(ifit)+".fit_func", "", "ch2/dof: "+to_string_with_precision(ch2_red_MEX_comb,5));
    Print_To_File({}, {Inv_Masses_to_print_MEX, FA.ave(), FA.err()}, "../data/ph_emission/"+ph_type+"/Bs_extr/comb_fits/FA_comb_ixg_"+to_string(ixg+1)+"_fit_"+to_string(ifit)+".fit_func", "", "ch2/dof: "+to_string_with_precision(ch2_red_MEX_comb,5));
    Print_To_File({}, {Inv_Masses_to_print_MEX, FV.ave(), FV.err()}, "../data/ph_emission/"+ph_type+"/Bs_extr/comb_fits/FV_comb_ixg_"+to_string(ixg+1)+"_fit_"+to_string(ifit)+".fit_func", "", "ch2/dof: "+to_string_with_precision(ch2_red_MEX_comb,5));
    


    Print_To_File({}, {Inv_Masses_to_print_MEX, FTA_d.ave(), FTA_d.err()}, "../data/ph_emission/"+ph_type+"/Bs_extr/comb_fits/FA_T_d_comb_ixg_"+to_string(ixg+1)+"_fit_"+to_string(ifit)+".fit_func", "", "ch2/dof: "+to_string_with_precision(ch2_red_MEX_comb,5));
    Print_To_File({}, {Inv_Masses_to_print_MEX, FTV_d.ave(), FTV_d.err()}, "../data/ph_emission/"+ph_type+"/Bs_extr/comb_fits/FV_T_d_comb_ixg_"+to_string(ixg+1)+"_fit_"+to_string(ifit)+".fit_func", "", "ch2/dof: "+to_string_with_precision(ch2_red_MEX_comb,5));
    Print_To_File({}, {Inv_Masses_to_print_MEX, FTA_u.ave(), FTA_u.err()}, "../data/ph_emission/"+ph_type+"/Bs_extr/comb_fits/FA_T_u_comb_ixg_"+to_string(ixg+1)+"_fit_"+to_string(ifit)+".fit_func", "", "ch2/dof: "+to_string_with_precision(ch2_red_MEX_comb,5));
    Print_To_File({}, {Inv_Masses_to_print_MEX, FTV_u.ave(), FTV_u.err()}, "../data/ph_emission/"+ph_type+"/Bs_extr/comb_fits/FV_T_u_comb_ixg_"+to_string(ixg+1)+"_fit_"+to_string(ifit)+".fit_func", "", "ch2/dof: "+to_string_with_precision(ch2_red_MEX_comb,5));
    Print_To_File({}, {Inv_Masses_to_print_MEX, FTA.ave(), FTA.err()}, "../data/ph_emission/"+ph_type+"/Bs_extr/comb_fits/FA_T_comb_ixg_"+to_string(ixg+1)+"_fit_"+to_string(ifit)+".fit_func", "", "ch2/dof: "+to_string_with_precision(ch2_red_MEX_comb,5));
    Print_To_File({}, {Inv_Masses_to_print_MEX, FTV.ave(), FTV.err()}, "../data/ph_emission/"+ph_type+"/Bs_extr/comb_fits/FV_T_comb_ixg_"+to_string(ixg+1)+"_fit_"+to_string(ifit)+".fit_func", "", "ch2/dof: "+to_string_with_precision(ch2_red_MEX_comb,5));
    
  }


  cout<<"Print fit parameters for quasi-pole"<<endl;
  cout<<"Axial: "<<(DA).ave()<< "+-"<<(DA).err()<<endl;
  cout<<"vector: "<<(DV).ave()<< "+-"<<(DV).err()<<endl;


  cout<<"First inverse moment of B-meson LCDA: LB: "<<(1.0/A).ave()<<" +- "<<(1.0/A).err()<<endl;

  cout<<"R1: "<<(R1).ave()<< "+-"<<(R1).err()<<endl;
  cout<<"R2: "<<(R2).ave()<< "+-"<<(R2).err()<<endl;
  
  cout<<"M1: "<<(M1).ave()<< "+-"<<(M1).err()<<endl;
  cout<<"M2: "<<(M2).ave()<< "+-"<<(M2).err()<<endl;
  cout<<"M3: "<<(M3).ave()<< "+-"<<(M3).err()<<endl;
  cout<<"M4: "<<(M4).ave()<< "+-"<<(M4).err()<<endl;

  cout<<"SV1: "<<(SV1).ave()<< "+-"<<(SV1).err()<<endl;
  cout<<"SV2: "<<(SV2).ave()<< "+-"<<(SV2).err()<<endl;
  cout<<"SV3: "<<(SV3).ave()<< "+-"<<(SV3).err()<<endl;

  cout<<"SA1: "<<(SA1).ave()<< "+-"<<(SA1).err()<<endl;
  cout<<"SA2: "<<(SA2).ave()<< "+-"<<(SA2).err()<<endl;
  cout<<"SA3: "<<(SA3).ave()<< "+-"<<(SA3).err()<<endl;
  
  cout<<"BV1: "<<(BV1).ave()<< "+-"<<(BV1).err()<<endl;
  cout<<"BV2: "<<(BV2).ave()<< "+-"<<(BV2).err()<<endl;
  
  cout<<"BA1: "<<(BA1).ave()<< "+-"<<(BA1).err()<<endl;
  cout<<"BA2: "<<(BA2).ave()<< "+-"<<(BA2).err()<<endl;


  cout<<"A1: "<<(A1).ave()<< "+-"<<(A1).err()<<endl;
  cout<<"A2: "<<(A2).ave()<< "+-"<<(A2).err()<<endl;

  cout<<"RTV1: "<<(RTV1).ave()<< "+-"<<(RTV1).err()<<endl;
  cout<<"RTV2: "<<(RTV2).ave()<< "+-"<<(RTV2).err()<<endl;

  cout<<"RTA1: "<<(RTA1).ave()<< "+-"<<(RTA1).err()<<endl;
  cout<<"RTA2: "<<(RTA2).ave()<< "+-"<<(RTA2).err()<<endl;


  cout<<"TVM1: "<<(TVM1).ave()<< "+-"<<(TVM1).err()<<endl;
  cout<<"TVE1: "<<(TVE1).ave()<< "+-"<<(TVE1).err()<<endl;

  cout<<"TAM1: "<<(TAM1).ave()<< "+-"<<(TAM1).err()<<endl;
  cout<<"TAE1: "<<(TAE1).ave()<< "+-"<<(TAE1).err()<<endl;

  

  cout<<"STV1: "<<(STV1).ave()<< "+-"<<(STV1).err()<<endl;
  cout<<"STV2: "<<(STV2).ave()<< "+-"<<(STV2).err()<<endl;
  cout<<"STV3: "<<(STV3).ave()<< "+-"<<(STV3).err()<<endl;

  cout<<"STA1: "<<(STA1).ave()<< "+-"<<(STA1).err()<<endl;
  cout<<"STA2: "<<(STA2).ave()<< "+-"<<(STA2).err()<<endl;
  cout<<"STA3: "<<(STA3).ave()<< "+-"<<(STA3).err()<<endl;
  
  cout<<"BTV1: "<<(BTV1).ave()<< "+-"<<(BTV1).err()<<endl;
  cout<<"BTV2: "<<(BTV2).ave()<< "+-"<<(BTV2).err()<<endl;
  
  cout<<"BTA1: "<<(BTA1).ave()<< "+-"<<(BTA1).err()<<endl;
  cout<<"BTA2: "<<(BTA2).ave()<< "+-"<<(BTA2).err()<<endl;
  
 
 


  //evaluate the rate

  }

  cout<<"ZT(2GeV)/ZT(MB) : "<<1.0/ZT_MB<<endl;
  cout<<"ZT(5GeV)/ZT(2GeV) NNNLO: "<<evolutor_ZT<<endl;
  cout<<"ZT(5GeV)/ZT(2GeV) NNLO: "<<evolutor_ZT_MS_bar(2.0, 5.0, 3)<<endl;
  cout<<"ZT(1GeV)/ZT(2GeV) NNNLO: "<<evolutor_ZT_MS_bar(2.0, 1.0, 4)<<endl;
  cout<<"ZT(1GeV)/ZT(2GeV) NNLO: "<<evolutor_ZT_MS_bar(2.0, 1.0, 3)<<endl;


  //perform Akaike analysis of the different fits
  distr_t_list FA_AIC(1), FV_AIC(1), FTA_AIC(1), FTV_AIC(1), FA_d_AIC(1), FV_d_AIC(1), FTA_d_AIC(1), FTV_d_AIC(1), FA_u_AIC(1), FV_u_AIC(1), FTA_u_AIC(1), FTV_u_AIC(1);

  //perform BMA uniform analysis of the different fits
  distr_t_list FA_d_uniform(1), FV_d_uniform(1), FTA_d_uniform(1), FTV_d_uniform(1), FA_u_uniform(1), FV_u_uniform(1), FTA_u_uniform(1), FTV_u_uniform(1);

  distr_t AVE_CA(1), AVE_CV(1), AVE_CA_T(1), AVE_R(1), AVE_RT(1), AVE_M3(1), AVE_M4(1);
  distr_t AVE_M1(1), AVE_M2(1), AVE_M1_T(1), AVE_M2_T(1), AVE_d_T(1);
  distr_t AVE_uniform_CA(1), AVE_uniform_CV(1), AVE_uniform_CA_T(1), AVE_uniform_R(1), AVE_uniform_RT(1), AVE_uniform_M3(1), AVE_uniform_M4(1);
  distr_t AVE_uniform_M1(1), AVE_uniform_M2(1), AVE_uniform_M1_T(1), AVE_uniform_M2_T(1), AVE_uniform_d_T(1);

  AVE_CA = AIC(CA_list, Ch2s, Ndofs, Nmeasurements, 0, mult_fits);
  AVE_CV = AIC(CV_list, Ch2s, Ndofs, Nmeasurements, 0, mult_fits);
  AVE_CA_T = AIC(CA_T_list, Ch2s, Ndofs, Nmeasurements, 0, mult_fits);
  AVE_R = AIC(R_list, Ch2s, Ndofs, Nmeasurements, 0, mult_fits);
  AVE_RT = AIC(RT_list, Ch2s, Ndofs, Nmeasurements, 0, mult_fits);
  AVE_M3 = AIC(M3_list, Ch2s, Ndofs, Nmeasurements, 0, mult_fits);
  AVE_M4 = AIC(M4_list, Ch2s, Ndofs, Nmeasurements, 0, mult_fits);
  AVE_M1 = AIC(M1_list, Ch2s, Ndofs, Nmeasurements, 0, mult_fits);
  AVE_M2 = AIC(M2_list, Ch2s, Ndofs, Nmeasurements, 0, mult_fits);
  AVE_M1_T = AIC(M1_T_list, Ch2s, Ndofs, Nmeasurements, 0, mult_fits);
  AVE_M2_T = AIC(M2_T_list, Ch2s, Ndofs, Nmeasurements, 0, mult_fits);
  AVE_d_T = AIC(d_T_list, Ch2s, Ndofs, Nmeasurements, 0, mult_fits);
 


  AVE_uniform_CA = BMA_uniform(CA_list, Ch2s, Ndofs, Nmeasurements, 1.4);
  AVE_uniform_CV = BMA_uniform(CV_list, Ch2s, Ndofs, Nmeasurements, 1.4);
  AVE_uniform_CA_T = BMA_uniform(CA_T_list, Ch2s, Ndofs, Nmeasurements, 1.4);
  AVE_uniform_R = BMA_uniform(R_list, Ch2s, Ndofs, Nmeasurements, 1.4);
  AVE_uniform_RT = BMA_uniform(RT_list, Ch2s, Ndofs, Nmeasurements, 1.4);
  AVE_uniform_M3 = BMA_uniform(M3_list, Ch2s, Ndofs, Nmeasurements, 1.4);
  AVE_uniform_M4 = BMA_uniform(M4_list, Ch2s, Ndofs, Nmeasurements, 1.4);
  AVE_uniform_M1 = BMA_uniform(M1_list, Ch2s, Ndofs, Nmeasurements, 1.4);
  AVE_uniform_M2 = BMA_uniform(M2_list, Ch2s, Ndofs, Nmeasurements, 1.4);
  AVE_uniform_M1_T = BMA_uniform(M1_T_list, Ch2s, Ndofs, Nmeasurements, 1.4);
  AVE_uniform_M2_T = BMA_uniform(M2_T_list, Ch2s, Ndofs, Nmeasurements, 1.4);
  AVE_uniform_d_T = BMA_uniform(d_T_list, Ch2s, Ndofs, Nmeasurements, 1.4);

 
  vector<double> W(Nfits,0.0);
  double wtot=0;
  for(int i=0; i<Nfits;i++) { W[i] = mult_fits[i]*exp(-0.5*( Ch2s[i] + 2*( Nmeasurements[i] - Ndofs[i]) - Nmeasurements[i])); wtot += W[i]; }
  for(int i=0; i<Nfits;i++) { W[i] /= wtot; }
  for(int i=0;i<Nfits;i++) cout<<"ifit: "<<i<<" mult: "<<mult_fits[i]<<" w[ifit]: "<<W[i]<<endl;

  //print averaged pars
  cout<<"################# PRINT AIC AVERAGED PARAMETERS: ################## "<<endl;
  cout<<"DA: "<<AVE_CA.ave()<<" "<<AVE_CA.err()<<endl;
  cout<<"DV: "<<AVE_CV.ave()<<" "<<AVE_CV.err()<<endl;
  cout<<"DA_T: "<<AVE_CA_T.ave()<<" "<<AVE_CA_T.err()<<endl;
  cout<<"R: "<<AVE_R.ave()<<" "<<AVE_R.err()<<endl;
  cout<<"RT: "<<AVE_RT.ave()<<" "<<AVE_RT.err()<<endl;
  cout<<"M3: "<<AVE_M3.ave()<<" "<<AVE_M3.err()<<endl;
  cout<<"M4: "<<AVE_M4.ave()<<" "<<AVE_M4.err()<<endl;
  cout<<"M1: "<<AVE_M1.ave()<<" "<<AVE_M1.err()<<endl;
  cout<<"M2: "<<AVE_M2.ave()<<" "<<AVE_M2.err()<<endl;
  cout<<"M1_T: "<<AVE_M1_T.ave()<<" "<<AVE_M1_T.err()<<endl;
  cout<<"M2_T: "<<AVE_M2_T.ave()<<" "<<AVE_M2_T.err()<<endl;
  cout<<"d_T: "<<AVE_d_T.ave()<<" "<<AVE_d_T.err()<<endl;
  cout<<"################################################################### "<<endl;
  cout<<"################# PRINT UNIFORM BMA AVERAGED PARAMETERS: ################## "<<endl;
  cout<<"DA: "<<AVE_uniform_CA.ave()<<" "<<AVE_uniform_CA.err()<<endl;
  cout<<"DV: "<<AVE_uniform_CV.ave()<<" "<<AVE_uniform_CV.err()<<endl;
  cout<<"DA_T: "<<AVE_uniform_CA_T.ave()<<" "<<AVE_uniform_CA_T.err()<<endl;
  cout<<"R: "<<AVE_uniform_R.ave()<<" "<<AVE_uniform_R.err()<<endl;
  cout<<"RT: "<<AVE_uniform_RT.ave()<<" "<<AVE_uniform_RT.err()<<endl;
  cout<<"M3: "<<AVE_uniform_M3.ave()<<" "<<AVE_uniform_M3.err()<<endl;
  cout<<"M4: "<<AVE_uniform_M4.ave()<<" "<<AVE_uniform_M4.err()<<endl;
  cout<<"M1: "<<AVE_uniform_M1.ave()<<" "<<AVE_uniform_M1.err()<<endl;
  cout<<"M2: "<<AVE_uniform_M2.ave()<<" "<<AVE_uniform_M2.err()<<endl;
  cout<<"M1_T: "<<AVE_uniform_M1_T.ave()<<" "<<AVE_uniform_M1_T.err()<<endl;
  cout<<"M2_T: "<<AVE_uniform_M2_T.ave()<<" "<<AVE_uniform_M2_T.err()<<endl;
  cout<<"d_T: "<<AVE_uniform_d_T.ave()<<" "<<AVE_uniform_d_T.err()<<endl;
  cout<<"################################################################### "<<endl;



  for(int ixg=0; ixg<nxgs_to_spline;ixg++) {

    FA_AIC.distr_list.push_back( AIC(FA_extr[ixg], Ch2s, Ndofs, Nmeasurements,0, mult_fits)) ;
    FV_AIC.distr_list.push_back( AIC(FV_extr[ixg], Ch2s, Ndofs, Nmeasurements,0, mult_fits)) ;
    FTA_AIC.distr_list.push_back( AIC(FTA_extr[ixg], Ch2s, Ndofs, Nmeasurements,0, mult_fits)) ;
    FTV_AIC.distr_list.push_back( AIC(FTV_extr[ixg], Ch2s, Ndofs, Nmeasurements,0, mult_fits)) ;

    
    FA_d_AIC.distr_list.push_back( AIC(FA_d_extr[ixg], Ch2s, Ndofs, Nmeasurements,0, mult_fits)) ;
    FV_d_AIC.distr_list.push_back( AIC(FV_d_extr[ixg], Ch2s, Ndofs, Nmeasurements,0, mult_fits)) ;
    FTA_d_AIC.distr_list.push_back( AIC(FTA_d_extr[ixg], Ch2s, Ndofs, Nmeasurements,0, mult_fits)) ;
    FTV_d_AIC.distr_list.push_back( AIC(FTV_d_extr[ixg], Ch2s, Ndofs, Nmeasurements,0, mult_fits)) ;


    FA_u_AIC.distr_list.push_back( AIC(FA_u_extr[ixg], Ch2s, Ndofs, Nmeasurements,0, mult_fits)) ;
    FV_u_AIC.distr_list.push_back( AIC(FV_u_extr[ixg], Ch2s, Ndofs, Nmeasurements,0, mult_fits)) ;
    FTA_u_AIC.distr_list.push_back( AIC(FTA_u_extr[ixg], Ch2s, Ndofs, Nmeasurements,0, mult_fits)) ;
    FTV_u_AIC.distr_list.push_back( AIC(FTV_u_extr[ixg], Ch2s, Ndofs, Nmeasurements,0, mult_fits)) ;


    FA_uniform.distr_list.push_back( BMA_uniform(FA_extr[ixg], Ch2s, Ndofs, Nmeasurements,1.4)) ;
    FV_uniform.distr_list.push_back( BMA_uniform(FV_extr[ixg], Ch2s, Ndofs, Nmeasurements,1.4)) ;
    FTA_uniform.distr_list.push_back( BMA_uniform(FTA_extr[ixg], Ch2s, Ndofs, Nmeasurements,1.4)) ;
    FTV_uniform.distr_list.push_back( BMA_uniform(FTV_extr[ixg], Ch2s, Ndofs, Nmeasurements,1.4)) ;

    
    FA_d_uniform.distr_list.push_back( BMA_uniform(FA_d_extr[ixg], Ch2s, Ndofs, Nmeasurements,1.4)) ;
    FV_d_uniform.distr_list.push_back( BMA_uniform(FV_d_extr[ixg], Ch2s, Ndofs, Nmeasurements,1.4)) ;
    FTA_d_uniform.distr_list.push_back( BMA_uniform(FTA_d_extr[ixg], Ch2s, Ndofs, Nmeasurements,1.4)) ;
    FTV_d_uniform.distr_list.push_back( BMA_uniform(FTV_d_extr[ixg], Ch2s, Ndofs, Nmeasurements,1.4)) ;


    FA_u_uniform.distr_list.push_back( BMA_uniform(FA_u_extr[ixg], Ch2s, Ndofs, Nmeasurements,1.4)) ;
    FV_u_uniform.distr_list.push_back( BMA_uniform(FV_u_extr[ixg], Ch2s, Ndofs, Nmeasurements,1.4)) ;
    FTA_u_uniform.distr_list.push_back( BMA_uniform(FTA_u_extr[ixg], Ch2s, Ndofs, Nmeasurements,1.4)) ;
    FTV_u_uniform.distr_list.push_back( BMA_uniform(FTV_u_extr[ixg], Ch2s, Ndofs, Nmeasurements,1.4)) ; 

  }

 

  vector<distr_t_list> FA_mass_dep_AIC, FV_mass_dep_AIC, FTA_mass_dep_AIC, FTV_mass_dep_AIC;
  vector<distr_t_list> FA_d_mass_dep_AIC, FV_d_mass_dep_AIC, FTA_d_mass_dep_AIC, FTV_d_mass_dep_AIC;
  vector<distr_t_list> FA_u_mass_dep_AIC, FV_u_mass_dep_AIC, FTA_u_mass_dep_AIC, FTV_u_mass_dep_AIC;

  vector<distr_t_list> FA_mass_dep_uniform, FV_mass_dep_uniform, FTA_mass_dep_uniform, FTV_mass_dep_uniform;
  vector<distr_t_list> FA_d_mass_dep_uniform, FV_d_mass_dep_uniform, FTA_d_mass_dep_uniform, FTV_d_mass_dep_uniform;
  vector<distr_t_list> FA_u_mass_dep_uniform, FV_u_mass_dep_uniform, FTA_u_mass_dep_uniform, FTV_u_mass_dep_uniform;

  for(int ixg=0;ixg<num_xg;ixg++) {

    FA_mass_dep_AIC.emplace_back(1);
    FV_mass_dep_AIC.emplace_back(1);
    FTA_mass_dep_AIC.emplace_back(1);
    FTV_mass_dep_AIC.emplace_back(1);

    FA_d_mass_dep_AIC.emplace_back(1);
    FV_d_mass_dep_AIC.emplace_back(1);
    FTA_d_mass_dep_AIC.emplace_back(1);
    FTV_d_mass_dep_AIC.emplace_back(1);

    FA_u_mass_dep_AIC.emplace_back(1);
    FV_u_mass_dep_AIC.emplace_back(1);
    FTA_u_mass_dep_AIC.emplace_back(1);
    FTV_u_mass_dep_AIC.emplace_back(1);

    FA_mass_dep_uniform.emplace_back(1);
    FV_mass_dep_uniform.emplace_back(1);
    FTA_mass_dep_uniform.emplace_back(1);
    FTV_mass_dep_uniform.emplace_back(1);

    FA_d_mass_dep_uniform.emplace_back(1);
    FV_d_mass_dep_uniform.emplace_back(1);
    FTA_d_mass_dep_uniform.emplace_back(1);
    FTV_d_mass_dep_uniform.emplace_back(1);

    FA_u_mass_dep_uniform.emplace_back(1);
    FV_u_mass_dep_uniform.emplace_back(1);
    FTA_u_mass_dep_uniform.emplace_back(1);
    FTV_u_mass_dep_uniform.emplace_back(1);
    
  }

  for(int ixg=0; ixg<num_xg;ixg++) {

    for(int im=0;im< (signed)Inv_Masses_to_print_MEX.size(); im++) {
    
      FA_mass_dep_AIC[ixg].distr_list.push_back(  AIC( FA_sims_m_dep[ixg][im], Ch2s, Ndofs, Nmeasurements, 0, mult_fits));
      FV_mass_dep_AIC[ixg].distr_list.push_back(  AIC( FV_sims_m_dep[ixg][im], Ch2s, Ndofs, Nmeasurements, 0, mult_fits));
      FTA_mass_dep_AIC[ixg].distr_list.push_back(  AIC( FTA_sims_m_dep[ixg][im], Ch2s, Ndofs, Nmeasurements, 0, mult_fits));
      FTV_mass_dep_AIC[ixg].distr_list.push_back(  AIC( FTV_sims_m_dep[ixg][im], Ch2s, Ndofs, Nmeasurements, 0, mult_fits));

      FA_d_mass_dep_AIC[ixg].distr_list.push_back(  AIC( FA_d_sims_m_dep[ixg][im], Ch2s, Ndofs, Nmeasurements, 0, mult_fits));
      FV_d_mass_dep_AIC[ixg].distr_list.push_back(  AIC( FV_d_sims_m_dep[ixg][im], Ch2s, Ndofs, Nmeasurements, 0, mult_fits));
      FTA_d_mass_dep_AIC[ixg].distr_list.push_back(  AIC( FTA_d_sims_m_dep[ixg][im], Ch2s, Ndofs, Nmeasurements, 0, mult_fits));
      FTV_d_mass_dep_AIC[ixg].distr_list.push_back(  AIC( FTV_d_sims_m_dep[ixg][im], Ch2s, Ndofs, Nmeasurements, 0, mult_fits));

      FA_u_mass_dep_AIC[ixg].distr_list.push_back(  AIC( FA_u_sims_m_dep[ixg][im], Ch2s, Ndofs, Nmeasurements, 0, mult_fits));
      FV_u_mass_dep_AIC[ixg].distr_list.push_back(  AIC( FV_u_sims_m_dep[ixg][im], Ch2s, Ndofs, Nmeasurements, 0, mult_fits));
      FTA_u_mass_dep_AIC[ixg].distr_list.push_back(  AIC( FTA_u_sims_m_dep[ixg][im], Ch2s, Ndofs, Nmeasurements, 0, mult_fits));
      FTV_u_mass_dep_AIC[ixg].distr_list.push_back(  AIC( FTV_u_sims_m_dep[ixg][im], Ch2s, Ndofs, Nmeasurements, 0, mult_fits));
      

      FA_mass_dep_uniform[ixg].distr_list.push_back(  BMA_uniform( FA_sims_m_dep[ixg][im], Ch2s, Ndofs, Nmeasurements, 1.4));
      FV_mass_dep_uniform[ixg].distr_list.push_back(  BMA_uniform( FV_sims_m_dep[ixg][im], Ch2s, Ndofs, Nmeasurements, 1.4));
      FTA_mass_dep_uniform[ixg].distr_list.push_back(  BMA_uniform( FTA_sims_m_dep[ixg][im], Ch2s, Ndofs, Nmeasurements, 1.4));
      FTV_mass_dep_uniform[ixg].distr_list.push_back(  BMA_uniform( FTV_sims_m_dep[ixg][im], Ch2s, Ndofs, Nmeasurements, 1.4));

      FA_d_mass_dep_uniform[ixg].distr_list.push_back(  BMA_uniform( FA_d_sims_m_dep[ixg][im], Ch2s, Ndofs, Nmeasurements, 1.4));
      FV_d_mass_dep_uniform[ixg].distr_list.push_back(  BMA_uniform( FV_d_sims_m_dep[ixg][im], Ch2s, Ndofs, Nmeasurements, 1.4));
      FTA_d_mass_dep_uniform[ixg].distr_list.push_back(  BMA_uniform( FTA_d_sims_m_dep[ixg][im], Ch2s, Ndofs, Nmeasurements, 1.4));
      FTV_d_mass_dep_uniform[ixg].distr_list.push_back(  BMA_uniform( FTV_d_sims_m_dep[ixg][im], Ch2s, Ndofs, Nmeasurements, 1.4));

      FA_u_mass_dep_uniform[ixg].distr_list.push_back(  BMA_uniform( FA_u_sims_m_dep[ixg][im], Ch2s, Ndofs, Nmeasurements, 1.4));
      FV_u_mass_dep_uniform[ixg].distr_list.push_back(  BMA_uniform( FV_u_sims_m_dep[ixg][im], Ch2s, Ndofs, Nmeasurements, 1.4));
      FTA_u_mass_dep_uniform[ixg].distr_list.push_back(  BMA_uniform( FTA_u_sims_m_dep[ixg][im], Ch2s, Ndofs, Nmeasurements, 1.4));
      FTV_u_mass_dep_uniform[ixg].distr_list.push_back(  BMA_uniform( FTV_u_sims_m_dep[ixg][im], Ch2s, Ndofs, Nmeasurements, 1.4));


    }

   
    
    //print the results

    //AIC
    
    Print_To_File({}, {Inv_Masses_to_print_MEX, FA_d_mass_dep_AIC[ixg].ave(), FA_d_mass_dep_AIC[ixg].err()}, "../data/ph_emission/"+ph_type+"/Bs_extr/FA_d_comb_ixg_"+to_string(ixg+1)+".fit_func", "", "");
    Print_To_File({}, {Inv_Masses_to_print_MEX, FV_d_mass_dep_AIC[ixg].ave(), FV_d_mass_dep_AIC[ixg].err()}, "../data/ph_emission/"+ph_type+"/Bs_extr/FV_d_comb_ixg_"+to_string(ixg+1)+".fit_func", "", "");
    Print_To_File({}, {Inv_Masses_to_print_MEX, FA_u_mass_dep_AIC[ixg].ave(), FA_u_mass_dep_AIC[ixg].err()}, "../data/ph_emission/"+ph_type+"/Bs_extr/FA_u_comb_ixg_"+to_string(ixg+1)+".fit_func", "", "");
    Print_To_File({}, {Inv_Masses_to_print_MEX, FV_u_mass_dep_AIC[ixg].ave(), FV_u_mass_dep_AIC[ixg].err()}, "../data/ph_emission/"+ph_type+"/Bs_extr/FV_u_comb_ixg_"+to_string(ixg+1)+".fit_func", "", "");
    //tot
    Print_To_File({}, {Inv_Masses_to_print_MEX, FA_mass_dep_AIC[ixg].ave(), FA_mass_dep_AIC[ixg].err()}, "../data/ph_emission/"+ph_type+"/Bs_extr/FA_comb_ixg_"+to_string(ixg+1)+".fit_func", "", "");
    Print_To_File({}, {Inv_Masses_to_print_MEX, FV_mass_dep_AIC[ixg].ave(), FV_mass_dep_AIC[ixg].err()}, "../data/ph_emission/"+ph_type+"/Bs_extr/FV_comb_ixg_"+to_string(ixg+1)+".fit_func", "", "");
    
    Print_To_File({}, {Inv_Masses_to_print_MEX, FTA_d_mass_dep_AIC[ixg].ave(), FTA_d_mass_dep_AIC[ixg].err()}, "../data/ph_emission/"+ph_type+"/Bs_extr/FA_T_d_comb_ixg_"+to_string(ixg+1)+".fit_func", "", "");
    Print_To_File({}, {Inv_Masses_to_print_MEX, FTV_d_mass_dep_AIC[ixg].ave(), FTV_d_mass_dep_AIC[ixg].err()}, "../data/ph_emission/"+ph_type+"/Bs_extr/FV_T_d_comb_ixg_"+to_string(ixg+1)+".fit_func", "", "");
    Print_To_File({}, {Inv_Masses_to_print_MEX, FTA_u_mass_dep_AIC[ixg].ave(), FTA_u_mass_dep_AIC[ixg].err()}, "../data/ph_emission/"+ph_type+"/Bs_extr/FA_T_u_comb_ixg_"+to_string(ixg+1)+".fit_func", "", "");
    Print_To_File({}, {Inv_Masses_to_print_MEX, FTV_u_mass_dep_AIC[ixg].ave(), FTV_u_mass_dep_AIC[ixg].err()}, "../data/ph_emission/"+ph_type+"/Bs_extr/FV_T_u_comb_ixg_"+to_string(ixg+1)+".fit_func", "", "");
    Print_To_File({}, {Inv_Masses_to_print_MEX, FTA_mass_dep_AIC[ixg].ave(), FTA_mass_dep_AIC[ixg].err()}, "../data/ph_emission/"+ph_type+"/Bs_extr/FA_T_comb_ixg_"+to_string(ixg+1)+".fit_func", "", "");
    Print_To_File({}, {Inv_Masses_to_print_MEX, FTV_mass_dep_AIC[ixg].ave(), FTV_mass_dep_AIC[ixg].err()}, "../data/ph_emission/"+ph_type+"/Bs_extr/FV_T_comb_ixg_"+to_string(ixg+1)+".fit_func", "", "");

    //uniform

    Print_To_File({}, {Inv_Masses_to_print_MEX, FA_d_mass_dep_uniform[ixg].ave(), FA_d_mass_dep_uniform[ixg].err()}, "../data/ph_emission/"+ph_type+"/Bs_extr/FA_d_comb_uniform_ixg_"+to_string(ixg+1)+".fit_func", "", "");
    Print_To_File({}, {Inv_Masses_to_print_MEX, FV_d_mass_dep_uniform[ixg].ave(), FV_d_mass_dep_uniform[ixg].err()}, "../data/ph_emission/"+ph_type+"/Bs_extr/FV_d_comb_uniform_ixg_"+to_string(ixg+1)+".fit_func", "", "");
    Print_To_File({}, {Inv_Masses_to_print_MEX, FA_u_mass_dep_uniform[ixg].ave(), FA_u_mass_dep_uniform[ixg].err()}, "../data/ph_emission/"+ph_type+"/Bs_extr/FA_u_comb_uniform_ixg_"+to_string(ixg+1)+".fit_func", "", "");
    Print_To_File({}, {Inv_Masses_to_print_MEX, FV_u_mass_dep_uniform[ixg].ave(), FV_u_mass_dep_uniform[ixg].err()}, "../data/ph_emission/"+ph_type+"/Bs_extr/FV_u_comb_uniform_ixg_"+to_string(ixg+1)+".fit_func", "", "");
    //tot
    Print_To_File({}, {Inv_Masses_to_print_MEX, FA_mass_dep_uniform[ixg].ave(), FA_mass_dep_uniform[ixg].err()}, "../data/ph_emission/"+ph_type+"/Bs_extr/FA_comb_uniform_ixg_"+to_string(ixg+1)+".fit_func", "", "");
    Print_To_File({}, {Inv_Masses_to_print_MEX, FV_mass_dep_uniform[ixg].ave(), FV_mass_dep_uniform[ixg].err()}, "../data/ph_emission/"+ph_type+"/Bs_extr/FV_comb_uniform_ixg_"+to_string(ixg+1)+".fit_func", "", "");
    
    Print_To_File({}, {Inv_Masses_to_print_MEX, FTA_d_mass_dep_uniform[ixg].ave(), FTA_d_mass_dep_uniform[ixg].err()}, "../data/ph_emission/"+ph_type+"/Bs_extr/FA_T_d_comb_uniform_ixg_"+to_string(ixg+1)+".fit_func", "", "");
    Print_To_File({}, {Inv_Masses_to_print_MEX, FTV_d_mass_dep_uniform[ixg].ave(), FTV_d_mass_dep_uniform[ixg].err()}, "../data/ph_emission/"+ph_type+"/Bs_extr/FV_T_d_comb_uniform_ixg_"+to_string(ixg+1)+".fit_func", "", "");
    Print_To_File({}, {Inv_Masses_to_print_MEX, FTA_u_mass_dep_uniform[ixg].ave(), FTA_u_mass_dep_uniform[ixg].err()}, "../data/ph_emission/"+ph_type+"/Bs_extr/FA_T_u_comb_uniform_ixg_"+to_string(ixg+1)+".fit_func", "", "");
    Print_To_File({}, {Inv_Masses_to_print_MEX, FTV_u_mass_dep_uniform[ixg].ave(), FTV_u_mass_dep_uniform[ixg].err()}, "../data/ph_emission/"+ph_type+"/Bs_extr/FV_T_u_comb_uniform_ixg_"+to_string(ixg+1)+".fit_func", "", "");
    Print_To_File({}, {Inv_Masses_to_print_MEX, FTA_mass_dep_uniform[ixg].ave(), FTA_mass_dep_uniform[ixg].err()}, "../data/ph_emission/"+ph_type+"/Bs_extr/FA_T_comb_uniform_ixg_"+to_string(ixg+1)+".fit_func", "", "");
    Print_To_File({}, {Inv_Masses_to_print_MEX, FTV_mass_dep_uniform[ixg].ave(), FTV_mass_dep_uniform[ixg].err()}, "../data/ph_emission/"+ph_type+"/Bs_extr/FV_T_comb_uniform_ixg_"+to_string(ixg+1)+".fit_func", "", "");
    
  }
 

  // print results at Bs as a function of xg

  Print_To_File({}, {xg_scan, (FA_AIC*f_Bs).ave(), (FA_AIC*f_Bs).err(),  FA_AIC.ave(), FA_AIC.err()},  "../data/ph_emission/"+ph_type+"/Bs_extr/FA_comb_Bs_fit.dat", "", "");
  Print_To_File({}, {xg_scan, (FV_AIC*f_Bs).ave(), (FV_AIC*f_Bs).err(), FV_AIC.ave(), FV_AIC.err()},  "../data/ph_emission/"+ph_type+"/Bs_extr/FV_comb_Bs_fit.dat", "", "");
  Print_To_File({}, {xg_scan, (FTA_AIC*f_Bs).ave(), (FTA_AIC*f_Bs).err(),  FTA_AIC.ave(), FTA_AIC.err()},  "../data/ph_emission/"+ph_type+"/Bs_extr/FA_T_comb_Bs_fit.dat", "", "");
  Print_To_File({}, {xg_scan, (FTV_AIC*f_Bs).ave(), (FTV_AIC*f_Bs).err(), FTV_AIC.ave(), FTV_AIC.err()},  "../data/ph_emission/"+ph_type+"/Bs_extr/FV_T_comb_Bs_fit.dat", "", "");
 
  
  Print_To_File({}, {xg_scan, (FA_d_AIC*f_Bs).ave(), (FA_d_AIC*f_Bs).err(),  FA_d_AIC.ave(), FA_d_AIC.err()},  "../data/ph_emission/"+ph_type+"/Bs_extr/FA_d_comb_Bs_fit.dat", "", "");
  Print_To_File({}, {xg_scan, (FV_d_AIC*f_Bs).ave(), (FV_d_AIC*f_Bs).err(), FV_d_AIC.ave(), FV_d_AIC.err()},  "../data/ph_emission/"+ph_type+"/Bs_extr/FV_d_comb_Bs_fit.dat", "", "");
  Print_To_File({}, {xg_scan, (FA_u_AIC*f_Bs).ave(), (FA_u_AIC*f_Bs).err(),  FA_u_AIC.ave(), FA_u_AIC.err()},  "../data/ph_emission/"+ph_type+"/Bs_extr/FA_u_comb_Bs_fit.dat", "", "");
  Print_To_File({}, {xg_scan, (FV_u_AIC*f_Bs).ave(), (FV_u_AIC*f_Bs).err(), FV_u_AIC.ave(), FV_u_AIC.err()},  "../data/ph_emission/"+ph_type+"/Bs_extr/FV_u_comb_Bs_fit.dat", "", "");
  
  Print_To_File({}, {xg_scan, (FTA_d_AIC*f_Bs).ave(), (FTA_d_AIC*f_Bs).err(),  FTA_d_AIC.ave(), FTA_d_AIC.err()},  "../data/ph_emission/"+ph_type+"/Bs_extr/FA_T_d_comb_Bs_fit.dat", "", "");
  Print_To_File({}, {xg_scan, (FTV_d_AIC*f_Bs).ave(), (FTV_d_AIC*f_Bs).err(), FTV_d_AIC.ave(), FTV_d_AIC.err()},  "../data/ph_emission/"+ph_type+"/Bs_extr/FV_T_d_comb_Bs_fit.dat", "", "");
  Print_To_File({}, {xg_scan, (FTA_u_AIC*f_Bs).ave(), (FTA_u_AIC*f_Bs).err(),  FTA_u_AIC.ave(), FTA_u_AIC.err()},  "../data/ph_emission/"+ph_type+"/Bs_extr/FA_T_u_comb_Bs_fit.dat", "", "");
  Print_To_File({}, {xg_scan, (FTV_u_AIC*f_Bs).ave(), (FTV_u_AIC*f_Bs).err(), FTV_u_AIC.ave(), FTV_u_AIC.err()},  "../data/ph_emission/"+ph_type+"/Bs_extr/FV_T_u_comb_Bs_fit.dat", "", "");


  Print_To_File({}, {xg_scan, (FA_uniform*f_Bs).ave(), (FA_uniform*f_Bs).err(),  FA_uniform.ave(), FA_uniform.err()},  "../data/ph_emission/"+ph_type+"/Bs_extr/FA_comb_uniform_Bs_fit.dat", "", "");
  Print_To_File({}, {xg_scan, (FV_uniform*f_Bs).ave(), (FV_uniform*f_Bs).err(), FV_uniform.ave(), FV_uniform.err()},  "../data/ph_emission/"+ph_type+"/Bs_extr/FV_comb_uniform_Bs_fit.dat", "", "");
  Print_To_File({}, {xg_scan, (FTA_uniform*f_Bs).ave(), (FTA_uniform*f_Bs).err(),  FTA_uniform.ave(), FTA_uniform.err()},  "../data/ph_emission/"+ph_type+"/Bs_extr/FA_T_comb_uniform_Bs_fit.dat", "", "");
  Print_To_File({}, {xg_scan, (FTV_uniform*f_Bs).ave(), (FTV_uniform*f_Bs).err(), FTV_uniform.ave(), FTV_uniform.err()},  "../data/ph_emission/"+ph_type+"/Bs_extr/FV_T_comb_uniform_Bs_fit.dat", "", "");
 
  
  Print_To_File({}, {xg_scan, (FA_d_uniform*f_Bs).ave(), (FA_d_uniform*f_Bs).err(),  FA_d_uniform.ave(), FA_d_uniform.err()},  "../data/ph_emission/"+ph_type+"/Bs_extr/FA_d_comb_uniform_Bs_fit.dat", "", "");
  Print_To_File({}, {xg_scan, (FV_d_uniform*f_Bs).ave(), (FV_d_uniform*f_Bs).err(), FV_d_uniform.ave(), FV_d_uniform.err()},  "../data/ph_emission/"+ph_type+"/Bs_extr/FV_d_comb_uniform_Bs_fit.dat", "", "");
  Print_To_File({}, {xg_scan, (FA_u_uniform*f_Bs).ave(), (FA_u_uniform*f_Bs).err(),  FA_u_uniform.ave(), FA_u_uniform.err()},  "../data/ph_emission/"+ph_type+"/Bs_extr/FA_u_comb_uniform_Bs_fit.dat", "", "");
  Print_To_File({}, {xg_scan, (FV_u_uniform*f_Bs).ave(), (FV_u_uniform*f_Bs).err(), FV_u_uniform.ave(), FV_u_uniform.err()},  "../data/ph_emission/"+ph_type+"/Bs_extr/FV_u_comb_uniform_Bs_fit.dat", "", "");
  
  Print_To_File({}, {xg_scan, (FTA_d_uniform*f_Bs).ave(), (FTA_d_uniform*f_Bs).err(),  FTA_d_uniform.ave(), FTA_d_uniform.err()},  "../data/ph_emission/"+ph_type+"/Bs_extr/FA_T_d_comb_uniform_Bs_fit.dat", "", "");
  Print_To_File({}, {xg_scan, (FTV_d_uniform*f_Bs).ave(), (FTV_d_uniform*f_Bs).err(), FTV_d_uniform.ave(), FTV_d_uniform.err()},  "../data/ph_emission/"+ph_type+"/Bs_extr/FV_T_d_comb_uniform_Bs_fit.dat", "", "");
  Print_To_File({}, {xg_scan, (FTA_u_uniform*f_Bs).ave(), (FTA_u_uniform*f_Bs).err(),  FTA_u_uniform.ave(), FTA_u_uniform.err()},  "../data/ph_emission/"+ph_type+"/Bs_extr/FA_T_u_comb_uniform_Bs_fit.dat", "", "");
  Print_To_File({}, {xg_scan, (FTV_u_uniform*f_Bs).ave(), (FTV_u_uniform*f_Bs).err(), FTV_u_uniform.ave(), FTV_u_uniform.err()},  "../data/ph_emission/"+ph_type+"/Bs_extr/FV_T_u_comb_uniform_Bs_fit.dat", "", "");


  //print jackknives of global fit
  boost::filesystem::create_directory("../data/ph_emission/"+ph_type+"/Bs_extr/jackknife");
  
  for(int ixg=0; ixg<nxgs_to_spline;ixg++) {
    Print_To_File({}, { (FA_uniform).distr_list[ixg].distr, (FV_uniform).distr_list[ixg].distr, (FTA_uniform).distr_list[ixg].distr, (FTV_uniform).distr_list[ixg].distr }, "../data/ph_emission/"+ph_type+"/Bs_extr/jackknife/FF_ixg_"+to_string(ixg)+".dat", "", "#jackks"); 
  }

  }
  else {
    for(int ixg=0;ixg<nxgs_to_spline;ixg++) {
      FA_uniform.distr_list.emplace_back(true, Read_From_File( "../data/ph_emission/"+ph_type+"/Bs_extr/jackknife/FF_ixg_"+to_string(ixg)+".dat", 1,5, 1));
      FV_uniform.distr_list.emplace_back(true, Read_From_File( "../data/ph_emission/"+ph_type+"/Bs_extr/jackknife/FF_ixg_"+to_string(ixg)+".dat", 2,5, 1));
      FTA_uniform.distr_list.emplace_back(true, Read_From_File( "../data/ph_emission/"+ph_type+"/Bs_extr/jackknife/FF_ixg_"+to_string(ixg)+".dat", 3,5, 1));
      FTV_uniform.distr_list.emplace_back(true, Read_From_File( "../data/ph_emission/"+ph_type+"/Bs_extr/jackknife/FF_ixg_"+to_string(ixg)+".dat", 4, 5, 1));
    }
  }


  
  
  if(Compute_rate) {

  //spline the results for the total

  int Nb=2000;



  vector<boost::math::interpolators::cardinal_cubic_b_spline<double>> FA_interpol;
  vector<boost::math::interpolators::cardinal_cubic_b_spline<double>> FV_interpol;

  vector<boost::math::interpolators::cardinal_cubic_b_spline<double>> FTA_interpol;
  vector<boost::math::interpolators::cardinal_cubic_b_spline<double>> FTV_interpol;

  vector<boost::math::interpolators::cardinal_cubic_b_spline<double>> pretty_b_interpol;


  
  
  for(int ijack=0;ijack<NJ;ijack++) {

    Vfloat FA_jacks, FV_jacks, FTA_jacks, FTV_jacks;

    Vfloat pretty_b_jacks;

    for(int ixg=0;ixg<nxgs_to_spline;ixg++) { FA_jacks.push_back( (FA_uniform*f_Bs).distr_list[ixg].distr[ijack]); FV_jacks.push_back( (FV_uniform*f_Bs).distr_list[ixg].distr[ijack]); FTA_jacks.push_back((FTA_uniform*f_Bs).distr_list[ixg].distr[ijack]); FTV_jacks.push_back( (FTV_uniform*f_Bs).distr_list[ixg].distr[ijack]); pretty_b_jacks.push_back( pretty_b_MB.distr_list[ixg].distr[ijack]);  }
 
    FA_interpol.emplace_back( FA_jacks.begin(), FA_jacks.end(), 0.0025, 0.0025);
    
    FV_interpol.emplace_back( FV_jacks.begin(), FV_jacks.end(), 0.0025, 0.0025);
    
    FTV_interpol.emplace_back( FTV_jacks.begin(), FTV_jacks.end(), 0.0025, 0.0025);
    
    FTA_interpol.emplace_back( FTA_jacks.begin(), FTA_jacks.end(), 0.0025, 0.0025);
    
    pretty_b_interpol.emplace_back( pretty_b_jacks.begin(), pretty_b_jacks.end(), 0.0025, 0.0025);
    
  }



  auto FA_boot = [&FA_interpol](double x, int i) {

    double ave=0;
    for(int ijack=0;ijack<NJ;ijack++) ave += FA_interpol[ijack](x)/NJ;
    return ave + sqrt(NJ-1.0)*(FA_interpol[i](x) - ave);
    
  };

  auto FV_boot = [&FV_interpol](double x, int i) {

    double ave=0;
    for(int ijack=0;ijack<NJ;ijack++) ave += FV_interpol[ijack](x)/NJ;
    return ave + sqrt(NJ-1.0)*(FV_interpol[i](x) - ave);
    
  };


   auto FTA_boot = [&FTA_interpol](double x, int i) {

    double ave=0;
    for(int ijack=0;ijack<NJ;ijack++) ave += FTA_interpol[ijack](x)/NJ;
    return ave + sqrt(NJ-1.0)*(FTA_interpol[i](x) - ave);
    
  };

  auto FTV_boot = [&FTV_interpol](double x, int i) {

    double ave=0;
    for(int ijack=0;ijack<NJ;ijack++) ave += FTV_interpol[ijack](x)/NJ;
    return ave + sqrt(NJ-1.0)*(FTV_interpol[i](x) - ave);
    
  };


  auto pretty_b_boot = [&pretty_b_interpol](double x, int i) {

    double ave=0;
    for(int ijack=0;ijack<NJ;ijack++) ave += pretty_b_interpol[ijack](x)/NJ;
    return ave + sqrt(NJ-1.0)*(pretty_b_interpol[i](x) - ave);
    
  };





   
  cout<<"SPLINE GENERATED"<<endl;
  
  

  MESON="Bs_extr";

  
    

   

  //int Nxgs=99;
  int Nxgs=30;
  Nxgs=53;
  //Nxgs=1;
  distr_t_list rate_Bs_SD(true, Nxgs), rate_Bs_INT(true,Nxgs), rate_Bs_diff_SD(true,Nxgs), AFB_Bs(true,Nxgs), rate_Bs_diff_PT(true,Nxgs), rate_Bs_diff_INT(true,Nxgs);


  distr_t_list rate_Bs_SD_FV_only(true,Nxgs);
  
  vector<distr_t_list> rate_Bs_SD_boot, rate_Bs_INT_boot, rate_Bs_diff_SD_boot, AFB_Bs_boot;
  vector<distr_t_list> rate_Bs_SD_FV_only_boot;
  distr_t_list rate_Bs_SD_NOCH(true,Nxgs), rate_Bs_INT_NOCH(true, Nxgs), rate_Bs_diff_SD_NOCH(true,Nxgs);
  vector<distr_t_list> rate_Bs_diff_PT_boot, rate_Bs_diff_INT_boot;
  distr_t_list rate_Bs_diff_INT_NOCH(true, Nxgs);
  vector<double> xg_min_list(Nxgs), xg_max_list(Nxgs), qmin_list(Nxgs);
  double Mb=5.36692;
  // for(int it=0;it<Nxgs;it++) { xg_min_list[it] = 0.0; xg_max_list[it] =  0.0129166666666*(it+1); qmin_list[it] = Mb*sqrt(1-xg_max_list[it]);  }
 

  //interval 4.0 to 5.3

  for(int it=0; it < Nxgs;it++) { qmin_list[it] = 4.0 + 0.025*it; xg_min_list[it] = 0.0; xg_max_list[it] = 1- pow(qmin_list[it]/Mb,2) ; }
  //for(int it=0; it<Nxgs;it++) { qmin_list[it] = 4.2; xg_min_list[it] = 0.0; xg_max_list[it] = 1- pow(qmin_list[it]/Mb,2) ; }

  //resize boot vectors
  for(int it=0;it<Nxgs;it++) {
    rate_Bs_SD_boot.emplace_back(true, Nb);
    rate_Bs_SD_FV_only_boot.emplace_back(true,Nb);
    rate_Bs_INT_boot.emplace_back(true,Nb);
    rate_Bs_diff_SD_boot.emplace_back(true, Nb);
    AFB_Bs_boot.emplace_back(true, Nb);
    rate_Bs_diff_INT_boot.emplace_back(true, Nb);
    rate_Bs_diff_PT_boot.emplace_back(true,Nb);
  }

 
  //##############       generate input parameters        #####################
  
  distr_t a1_distr(false) ;   distr_t mb_distr(false);
  distr_t xb_distr(false), Vtbs_distr(false), tbs_distr(false);

  distr_t a1_distr_std(true) ;   distr_t mb_distr_std(true);
  distr_t xb_distr_std(true), Vtbs_distr_std(true), tbs_distr_std(true);

  distr_t fbs_boot(false);
  
  vector<vector<double>> kappa_distr(Nb), Th_distr(Nb), W_distr(Nb), Br_distr(Nb);
  GaussianMersenne GPP(854135);
  UniformMersenne  RAN(45252, 0.0, 2*M_PI);
  RandomMersenne RM(6114743,NJ-1);
  vector<int> Rbo;
  for(int i=0;i<Nb;i++) Rbo.push_back( RM());


  for(int i=0;i<5000;i++) {
    cout<<"theta: "<<RAN()<<endl;
  }

  distr_t pretty_s_RE(1), pretty_s_IM(1);

  for(int i=0;i<NJ;i++) { pretty_s_RE.distr.push_back( -0.017 + GPP()*0.017/sqrt(NJ-1.0)) ;};
  for(int i=0;i<NJ;i++) { pretty_s_IM.distr.push_back( 0.02 + GPP()*0.02/sqrt(NJ-1.0)) ; };

  for(int ib=0;ib<Nb;ib++) mb_distr.distr.push_back( 4.18 + 0.03*GPP());
  for(int ib=0;ib<Nb;ib++) a1_distr.distr.push_back( -0.13 + GPP()*0.13);
  for(int ib=0;ib<Nb;ib++) xb_distr.distr.push_back( x_b + GPP()*x_b_err);
  for(int ib=0;ib<Nb;ib++) Vtbs_distr.distr.push_back( Vtbs + GPP()*Vtbs_err);
  for(int ib=0;ib<Nb;ib++) tbs_distr.distr.push_back( tBs + GPP()*tBs_err);

  for(int ib=0;ib<NJ;ib++) mb_distr_std.distr.push_back( 4.18 + (1.0/sqrt(NJ-1.0))*0.03*GPP());
  for(int ib=0;ib<NJ;ib++) a1_distr_std.distr.push_back( -0.13 + (1.0/sqrt(NJ-1.0))*GPP()*0.13);
  for(int ib=0;ib<NJ;ib++) xb_distr_std.distr.push_back( x_b + (1.0/sqrt(NJ-1.0))*GPP()*x_b_err);
  for(int ib=0;ib<NJ;ib++) Vtbs_distr_std.distr.push_back( Vtbs + (1.0/sqrt(NJ-1.0))*GPP()*Vtbs_err);
  for(int ib=0;ib<NJ;ib++) tbs_distr_std.distr.push_back( tBs + (1.0/sqrt(NJ-1.0))*GPP()*tBs_err);

  for(int nres=0;nres<8;nres++) {
    for(int ib=0;ib<Nb;ib++) kappa_distr[ib].push_back( 0.75 + 0.75*GPP());
    for(int ib=0;ib<Nb;ib++) Th_distr[ib].push_back( RAN());
    for(int ib=0;ib<Nb;ib++) W_distr[ib].push_back( GPP());
    for(int ib=0;ib<Nb;ib++) Br_distr[ib].push_back( GPP());
  }

  distr_t G_FV_KMN(false), G_FA_KMN(false), G_FTV_KMN(false), G_FTA_KMN(false), G_barF_KMN(false);

  //for KMN
  for(int iboot=0;iboot<Nb;iboot++) {
    G_FV_KMN.distr.push_back( GPP());
    G_FA_KMN.distr.push_back( GPP());
    G_FTV_KMN.distr.push_back( GPP());
    G_FTA_KMN.distr.push_back( GPP());
    G_barF_KMN.distr.push_back( GPP());
    fbs_boot.distr.push_back( f_Bs.ave() + f_Bs.err()*GPP());
  }

  //read JSON file for JPZ data
  // Open the file for reading
  //FILE* fp = fopen("../data/ph_emission/rph/Bs_extr/BtoGam_FF.json", "r");
 
  // Use a FileReadStream to
  // read the data from the file
  //char readBuffer[65536];
  //rapidjson::FileReadStream is(fp, readBuffer, sizeof(readBuffer));

  ifstream ifs("../data/ph_emission/rph/Bs_extr/BtoGam_FF.json");

  if (!ifs.is_open())
    {
      throw runtime_error("..");
    }
  
  rapidjson::IStreamWrapper isw(ifs);
  
  // Parse the JSON data
  // using a Document object
  rapidjson::Document d;
  d.ParseStream(isw);

  // Check if the document is valid 
    if (d.HasParseError()) { 
        cerr << "Error: failed to parse JSON document"
	     << endl;
	cout<<rapidjson::GetParseError_En( d.GetParseError() )<<" at offset: "<<d.GetErrorOffset()<<endl;
        ifs.close(); 
        exit(-1); 
    } 


  ifs.close();
  
  // Close the file
  //fclose(fp);

  cout<<d["Bs"].GetType()<<endl;

  vector<string> FF_names({"Tpara", "Tperp", "Vpara", "Vperp"});

  vector<vector<double>> A_val(4), A_err(4);

  
  
  for(int i=0;i<4;i++) { A_val[i].resize(4); A_err[i].resize(4);}
  
 
 
 
  
 
  for(int i=0;i<4;i++) {
    for(int s=0; s <(signed)FF_names.size(); s++) {
      
      A_val[s][i]  = d["Bs"]["central"][FF_names[s].c_str()][("a"+to_string(i)).c_str()].GetDouble();
      A_err[s][i]  = d["Bs"]["uncertainty"][FF_names[s].c_str()][("a"+to_string(i)).c_str()].GetDouble();
    }
  }

 

  Eigen::VectorXd params_val(16);

  for(int s=0; s<4;s++) {
    for(int i=0;i<4;i++) {
      params_val(4*s+i) = A_val[s][i];
    }
  }



  //Read covariance matrix
  Eigen::MatrixXd Cov_JPZ(16, 16);

  for(int s1=0;s1<(signed)FF_names.size();s1++) {
    for(int s2=0;s2<(signed)FF_names.size(); s2++) {
      for(int i=0;i<4;i++) {
	for(int j=0; j<4; j++) {
	  Cov_JPZ( 4*s1 + i , 4*s2+ j) = d["Bs"]["covariance"][(FF_names[s1]+FF_names[s2]).c_str()][("a"+to_string(i)+"a"+to_string(j)).c_str()].GetDouble();
	}
      }
    }
  }



  /*
  for(int s1=0;s1<(signed)FF_names.size();s1++) {
    for(int s2=0;s2<(signed)FF_names.size(); s2++) {
      for(int i=0;i<4;i++) {
	for(int j=0; j<4; j++) {
	  Cov_JPZ(4*s1+ i , 4*s2+ j) = Cov_JPZ(4*s1+i, 4*s2+j)*A_err[s1][i]*A_err[s2][j];
	}
      }
    }
    }*/

  vector<vector<double>> Tpara_boot, Tperp_boot, Vpara_boot, Vperp_boot;

  
  GaussianMersenne GJPZ(231450765);
  for(int iboot=0;iboot<Nb;iboot++) {
    Vfloat params=Covariate(Cov_JPZ, params_val, GJPZ);
    vector<vector<double>> TT(4) ;
    for(auto &T: TT) T.resize(4);
    
    for(int s=0;s<4;s++) {
      for(int i=0;i<4;i++) {
	TT[s][i] = params[4*s+i];
      }
    }

    Tpara_boot.push_back(TT[0]);
    Tperp_boot.push_back(TT[1]);
    Vpara_boot.push_back(TT[2]);
    Vperp_boot.push_back(TT[3]);
  }

   
  if(Get_KMN_JPZ_FF) {
  
  

    
  
  
  distr_t_list rate_KMN(false, Nxgs), rate_JPZ(false, Nxgs);

  

 

  //determine FF from KMN
  #pragma omp parallel for schedule(dynamic)
  for(int it=0;it< Nxgs;it++) {
    for(int iboot=0;iboot<Nb;iboot++) {


      //##################### KMN ######################

   
      auto FV_KMN = [&iboot, &G_FV_KMN](double x) {
	double A= 0.111;
	double s1=0.144;
	double s2=0.722;
	double MR2= 5.415;
	MR2 *= MR2;
	double q2= (1-x)*MBs*MBs;
	double v1 = A/(  ( 1 - q2/MR2)*(1 - s1*q2/MR2 + s2*pow( q2/MR2,2)));
	double v2 = A/( 1- q2/MR2);
	return v1 + G_FV_KMN.distr[iboot]*0.68*fabs( v2-v1);
	
      };

      auto FA_KMN = [&iboot, &G_FA_KMN](double x) {
	double A= 0.069;
	double s1=-0.031;
	double s2=0.384;
	double MR2= 5.829;
	MR2 *= MR2;
	double q2= (1-x)*MBs*MBs;
	double v1 = A/(  ( 1 - q2/MR2)*(1 - s1*q2/MR2 + s2*pow( q2/MR2,2)));
	double v2 = A/( 1- q2/MR2);
	return v1 + G_FA_KMN.distr[iboot]*0.68*fabs( v2-v1);
	
      };

      auto FTA_KMN = [&iboot, &G_FTA_KMN, &G_barF_KMN](double x) {
	double A= 0.119;
	double s1=-0.063;
	double s2=0.321;
	double MR2= 5.829;
	MR2 *= MR2;
	double q2= (1-x)*MBs*MBs;
	double v1 = A/(  ( 1 - q2/MR2)*(1 - s1*q2/MR2 + s2*pow( q2/MR2,2)));
	double v2 = A/( 1- q2/MR2);
	double fem = (-1.0/3.0)*(0.2205 + G_barF_KMN.distr[iboot]*0.0045);
	double g= -0.27 + G_barF_KMN.distr[iboot]*0.01;
	double Mphi=1.019;
	double Gamma=0.00427;
	double barF= A -2*fem*g*(q2/Mphi)*(q2-Mphi*Mphi)/( pow(q2 - Mphi*Mphi,2) + pow(Mphi*Gamma,2)); 

	
	return v1 + G_FTA_KMN.distr[iboot]*0.68*fabs( v2-v1) + barF;
	
      };

      auto FTV_KMN = [&iboot, &G_FTV_KMN, &G_barF_KMN](double x) {
	double A= 0.119;
	double s1=0.163;
	double s2=0.751;
	double MR2= 5.415;
	MR2 *= MR2;
	double q2= (1-x)*MBs*MBs;
	double v1 = A/(  ( 1 - q2/MR2)*(1 - s1*q2/MR2 + s2*pow( q2/MR2,2)));
	double v2 = A/( 1- q2/MR2);

	double fem = (-1.0/3.0)*(0.2205 + G_barF_KMN.distr[iboot]*0.0045);
	double g= -0.27 + G_barF_KMN.distr[iboot]*0.01;
	double Mphi=1.019;
	double Gamma=0.00427;
	double barF= A -2*fem*g*(q2/Mphi)*(q2-Mphi*Mphi)/( pow(q2 - Mphi*Mphi,2) + pow(Mphi*Gamma,2)); 
	
	return v1 + G_FTV_KMN.distr[iboot]*0.68*fabs( v2-v1) + barF;
	
      };

      auto bar_F_IM_KMN = [&iboot, &G_barF_KMN](double x) {

	double A= 0.119;
	double q2= (1-x)*MBs*MBs;
	double fem = (-1.0/3.0)*(0.2205 + G_barF_KMN.distr[iboot]*0.0045);
	double g= -0.27 + G_barF_KMN.distr[iboot]*0.01;
	double Mphi=1.019;
	double Gamma=0.00427;
	double barF= A +  2*fem*g*(q2/Mphi)*Mphi*Gamma/( pow(q2 - Mphi*Mphi,2) + pow(Mphi*Gamma,2)); 
	return barF;

      };


      //####################  JPZ  ##################################

    
      auto FV_JPZ = [&iboot, &Vperp_boot](double x) {

	double tp= pow(MBs + 0.770,2);
	double tm= pow(MBs - 0.770,2);
	double t0= tp*( 1 -sqrt( 1.0 - (tm/tp)));
	double q2=(1-x)*MBs*MBs;
	double z= (sqrt(tp-q2) -sqrt(tp-t0))/(sqrt(tp-q2) + sqrt(tp-t0));
	double z0= (sqrt(tp-0.0) -sqrt(tp-t0))/(sqrt(tp-0.0) + sqrt(tp-t0));
	double dz= z-z0;
	double res= 0.0;
	for(int i=0;i<4;i++) res += Vperp_boot[iboot][i]*pow(dz,i);
	double Mstar2= pow(5.4154,2);
	return res/(1- q2/Mstar2);
      };

       auto FA_JPZ = [&iboot, &Vpara_boot](double x) {

	 double tp= pow(MBs + 0.770,2);
	 double tm= pow(MBs - 0.770,2);
	 double t0= tp*( 1 -sqrt( 1.0 - (tm/tp)));
	 double q2=(1-x)*MBs*MBs;
	 double z= (sqrt(tp-q2) -sqrt(tp-t0))/(sqrt(tp-q2) + sqrt(tp-t0));
	 double z0= (sqrt(tp-0.0) -sqrt(tp-t0))/(sqrt(tp-0.0) + sqrt(tp-t0));
	 double dz= z-z0;
	 double res= 0.0;
	 for(int i=0;i<4;i++) res += Vpara_boot[iboot][i]*pow(dz,i);
	 double Mstar2= pow(5.8287,2);
	 return res/(1- q2/Mstar2);
       };
       
       
       auto FTV_JPZ = [&iboot, &Tperp_boot](double x) {
	 
	 double tp= pow(MBs + 0.770,2);
	 double tm= pow(MBs - 0.770,2);
	 double t0= tp*( 1 -sqrt( 1.0 - (tm/tp)));
	 double q2=(1-x)*MBs*MBs;
	 double z= (sqrt(tp-q2) -sqrt(tp-t0))/(sqrt(tp-q2) + sqrt(tp-t0));
	 double z0= (sqrt(tp-0.0) -sqrt(tp-t0))/(sqrt(tp-0.0) + sqrt(tp-t0));
	 double dz= z-z0;
	 double res= 0.0;
	 for(int i=0;i<4;i++) res += Tperp_boot[iboot][i]*pow(dz,i);
	 double Mstar2= pow(5.4154,2);
	 return res/(1- q2/Mstar2);
       };
       
       auto FTA_JPZ = [&iboot, &Tpara_boot](double x) {

	 double tp= pow(MBs + 0.770,2);
	 double tm= pow(MBs - 0.770,2);
	 double t0= tp*( 1 -sqrt( 1.0 - (tm/tp)));
	 double q2=(1-x)*MBs*MBs;
	 double z= (sqrt(tp-q2) -sqrt(tp-t0))/(sqrt(tp-q2) + sqrt(tp-t0));
	 double z0= (sqrt(tp-0.0) -sqrt(tp-t0))/(sqrt(tp-0.0) + sqrt(tp-t0));
	 double dz= z-z0;
	 double res= 0.0;
	 for(int i=0;i<4;i++) res += Tpara_boot[iboot][i]*pow(dz,i);
	 double Mstar2= pow(5.8287,2);
	 return res/(1- q2/Mstar2);
       };

       auto bar_F_IM_JPZ = [](double x) {
	 return 0.0;
       };
       
    
      
       
       double a1_val= a1_distr.distr[iboot];
       double mb_val= mb_distr.distr[iboot];
       double xb_val= xb_distr.distr[iboot];
       double Vtbs_val= Vtbs_distr.distr[iboot];
       double tbs_val = tbs_distr.distr[iboot];
       Vfloat kappa_list = kappa_distr[iboot];
       Vfloat Th_list = Th_distr[iboot];
       Vfloat W_list = W_distr[iboot];
       Vfloat Br_list = Br_distr[iboot];
       
      
       
       auto FV_T_eff_KMN = [&FTV_KMN, &a1_val, &fbs_boot, &mb_val, &iboot ](double x) { return FTV_KMN(x) +  0.0*(16.0/(3.0*WC(7,0).first))*(a1_val*fbs_boot.distr[iboot]/mb_val);};

       auto FV_T_eff_JPZ = [&FTV_JPZ, &a1_val, &fbs_boot, &mb_val, &iboot ](double x) { return FTV_JPZ(x) +  0.0*(16.0/(3.0*WC(7,0).first))*(a1_val*fbs_boot.distr[iboot]/mb_val);};

     
      double rate_SD_KMN =  Compute_Bs_mumugamma_decay_rate( FV_KMN, FA_KMN, FV_T_eff_KMN, FTA_KMN, bar_F_IM_KMN, bar_F_IM_KMN , fbs_boot.distr[iboot] , xb_val, Vtbs_val, tbs_val, kappa_list, Th_list, W_list, Br_list, "SD",  "CH",   xg_max_list[it], xg_min_list[it] );

      double rate_INT_KMN = Compute_Bs_mumugamma_decay_rate( FV_KMN, FA_KMN, FV_T_eff_KMN, FTA_KMN, bar_F_IM_KMN, bar_F_IM_KMN , fbs_boot.distr[iboot] , xb_val, Vtbs_val, tbs_val, kappa_list, Th_list, W_list, Br_list, "INT",  "CH",   xg_max_list[it], xg_min_list[it] );

      rate_KMN.distr_list[it].distr.push_back(rate_SD_KMN+ rate_INT_KMN);

    
      double rate_SD_JPZ =  Compute_Bs_mumugamma_decay_rate( FV_JPZ, FA_JPZ, FV_T_eff_JPZ, FTA_JPZ, bar_F_IM_JPZ, bar_F_IM_JPZ , fbs_boot.distr[iboot] , xb_val, Vtbs_val, tbs_val, kappa_list, Th_list, W_list, Br_list, "SD",  "CH",   xg_max_list[it], xg_min_list[it] );

      double rate_INT_JPZ = Compute_Bs_mumugamma_decay_rate( FV_JPZ, FA_JPZ, FV_T_eff_JPZ, FTA_JPZ, bar_F_IM_JPZ, bar_F_IM_JPZ , fbs_boot.distr[iboot] , xb_val, Vtbs_val, tbs_val, kappa_list, Th_list, W_list, Br_list, "INT",  "CH",   xg_max_list[it], xg_min_list[it] );

      rate_JPZ.distr_list[it].distr.push_back(rate_SD_JPZ+rate_INT_JPZ);

        
	 
          

    }

    cout<<"ixg: "<<it<<" completed!"<<endl;

  }

  distr_t_list FV_plot_KMN(false), FA_plot_KMN(false), FTV_plot_KMN(false), FTA_plot_KMN(false);

  distr_t_list FV_plot_JPZ(false), FA_plot_JPZ(false), FTV_plot_JPZ(false), FTA_plot_JPZ(false);

  for(auto &xg: xg_scan) {

    distr_t FV_KMN_D(false), FA_KMN_D(false), FTV_KMN_D(false), FTA_KMN_D(false);

    distr_t FV_JPZ_D(false), FA_JPZ_D(false), FTV_JPZ_D(false), FTA_JPZ_D(false);
    
    for(int iboot=0;iboot<Nb;iboot++) {

       auto FV_KMN = [&iboot, &G_FV_KMN](double x) {
	double A= 0.111;
	double s1=0.144;
	double s2=0.722;
	double MR2= 5.415;
	MR2 *= MR2;
	double q2= (1-x)*MBs*MBs;
	double v1 = A/(  ( 1 - q2/MR2)*(1 - s1*q2/MR2 + s2*pow( q2/MR2,2)));
	double v2 = A/( 1- q2/MR2);
	return v1 + G_FV_KMN.distr[iboot]*0.68*fabs( v2-v1);
	
      };

      auto FA_KMN = [&iboot, &G_FA_KMN](double x) {
	double A= 0.069;
	double s1=-0.031;
	double s2=0.384;
	double MR2= 5.829;
	MR2 *= MR2;
	double q2= (1-x)*MBs*MBs;
	double v1 = A/(  ( 1 - q2/MR2)*(1 - s1*q2/MR2 + s2*pow( q2/MR2,2)));
	double v2 = A/( 1- q2/MR2);
	return v1 + G_FA_KMN.distr[iboot]*0.68*fabs( v2-v1);
	
      };

      auto FTA_KMN = [&iboot, &G_FTA_KMN](double x) {
	double A= 0.119;
	double s1=-0.063;
	double s2=0.321;
	double MR2= 5.829;
	MR2 *= MR2;
	double q2= (1-x)*MBs*MBs;
	double v1 = A/(  ( 1 - q2/MR2)*(1 - s1*q2/MR2 + s2*pow( q2/MR2,2)));
	double v2 = A/( 1- q2/MR2);
	
	return v1 + G_FTA_KMN.distr[iboot]*0.68*fabs( v2-v1);
	
      };

      auto FTV_KMN = [&iboot, &G_FTV_KMN](double x) {
	double A= 0.119;
	double s1=0.163;
	double s2=0.751;
	double MR2= 5.415;
	MR2 *= MR2;
	double q2= (1-x)*MBs*MBs;
	double v1 = A/(  ( 1 - q2/MR2)*(1 - s1*q2/MR2 + s2*pow( q2/MR2,2)));
	double v2 = A/( 1- q2/MR2);

	
	return v1 + G_FTV_KMN.distr[iboot]*0.68*fabs( v2-v1);
	
      };


        auto FV_JPZ = [&iboot, &Vperp_boot](double x) {

	double tp= pow(MBs + 0.770,2);
	double tm= pow(MBs - 0.770,2);
	double t0= tp*( 1 -sqrt( 1 - (tm/tp)));
	double q2=(1-x)*MBs*MBs;
	double z= (sqrt(tp-q2) -sqrt(tp-t0))/(sqrt(tp-q2) + sqrt(tp-t0));
	double z0= (sqrt(tp-0.0) -sqrt(tp-t0))/(sqrt(tp-0.0) + sqrt(tp-t0));
	double dz= z-z0;
	double res= 0.0;
	for(int i=0;i<4;i++) res += Vperp_boot[iboot][i]*pow(dz,i);
	double Mstar2= pow(5.4154,2);
	return res/(1- q2/Mstar2);
      };

       auto FA_JPZ = [&iboot, &Vpara_boot](double x) {

	 double tp= pow(MBs + 0.770,2);
	 double tm= pow(MBs - 0.770,2);
	 double t0= tp*( 1 -sqrt( 1 - (tm/tp)));
	 double q2=(1-x)*MBs*MBs;
	 double z= (sqrt(tp-q2) -sqrt(tp-t0))/(sqrt(tp-q2) + sqrt(tp-t0));
	 double z0= (sqrt(tp-0.0) -sqrt(tp-t0))/(sqrt(tp-0.0) + sqrt(tp-t0));
	 double dz= z-z0;
	 double res= 0.0;
	 for(int i=0;i<4;i++) res += Vpara_boot[iboot][i]*pow(dz,i);
	 double Mstar2= pow(5.8287,2);
	 return res/(1- q2/Mstar2);
       };
       
       
       auto FTV_JPZ = [&iboot, &Tperp_boot](double x) {
	 
	 double tp= pow(MBs + 0.770,2);
	 double tm= pow(MBs - 0.770,2);
	 double t0= tp*( 1 -sqrt( 1 - (tm/tp)));
	 double q2=(1-x)*MBs*MBs;
	 double z= (sqrt(tp-q2) -sqrt(tp-t0))/(sqrt(tp-q2) + sqrt(tp-t0));
	 double z0= (sqrt(tp-0.0) -sqrt(tp-t0))/(sqrt(tp-0.0) + sqrt(tp-t0));
	 double dz= z-z0;
	 double res= 0.0;
	 for(int i=0;i<4;i++) res += Tperp_boot[iboot][i]*pow(dz,i);
	 double Mstar2= pow(5.4154,2);
	 return res/(1- q2/Mstar2);
       };
       
       auto FTA_JPZ = [&iboot, &Tpara_boot](double x) {

	 double tp= pow(MBs + 0.770,2);
	 double tm= pow(MBs - 0.770,2);
	 double t0= tp*( 1 -sqrt( 1 - (tm/tp)));
	 double q2=(1-x)*MBs*MBs;
	 double z= (sqrt(tp-q2) -sqrt(tp-t0))/(sqrt(tp-q2) + sqrt(tp-t0));
	 double z0= (sqrt(tp-0.0) -sqrt(tp-t0))/(sqrt(tp-0.0) + sqrt(tp-t0));
	 double dz= z-z0;
	 double res= 0.0;
	 for(int i=0;i<4;i++) res += Tpara_boot[iboot][i]*pow(dz,i);
	 double Mstar2= pow(5.8287,2);
	 return res/(1- q2/Mstar2);
      };
      

      FV_KMN_D.distr.push_back( FV_KMN(xg));
      FA_KMN_D.distr.push_back( FA_KMN(xg));
      FTV_KMN_D.distr.push_back( FTV_KMN(xg));
      FTA_KMN_D.distr.push_back( FTA_KMN(xg));

      FV_JPZ_D.distr.push_back( FV_JPZ(xg));
      FA_JPZ_D.distr.push_back( FA_JPZ(xg));
      FTV_JPZ_D.distr.push_back( FTV_JPZ(xg));
      FTA_JPZ_D.distr.push_back( FTA_JPZ(xg));
      
    }

    FV_plot_KMN.distr_list.push_back( FV_KMN_D);
    FA_plot_KMN.distr_list.push_back( FA_KMN_D);
    FTV_plot_KMN.distr_list.push_back( FTV_KMN_D);
    FTA_plot_KMN.distr_list.push_back( FTA_KMN_D);

    FV_plot_JPZ.distr_list.push_back( FV_JPZ_D);
    FA_plot_JPZ.distr_list.push_back( FA_JPZ_D);
    FTV_plot_JPZ.distr_list.push_back( FTV_JPZ_D);
    FTA_plot_JPZ.distr_list.push_back( FTA_JPZ_D);
  }

  boost::filesystem::create_directory("../data/ph_emission/rph/Bs_extr/KMN");

  Print_To_File({}, {xg_scan, FV_plot_KMN.ave(), FV_plot_KMN.err()}, "../data/ph_emission/rph/Bs_extr/KMN/FV.dat", "", "");
  Print_To_File({}, {xg_scan, FA_plot_KMN.ave(), FA_plot_KMN.err()}, "../data/ph_emission/rph/Bs_extr/KMN/FA.dat", "", "");
  Print_To_File({}, {xg_scan, FTV_plot_KMN.ave(), FTV_plot_KMN.err()}, "../data/ph_emission/rph/Bs_extr/KMN/FTV.dat", "", "");
  Print_To_File({}, {xg_scan, FTA_plot_KMN.ave(), FTA_plot_KMN.err()}, "../data/ph_emission/rph/Bs_extr/KMN/FTA.dat", "", "");

  Print_To_File({}, {qmin_list, xg_min_list, xg_max_list, rate_KMN.ave(), rate_KMN.err()}, "../data/ph_emission/rph/Bs_extr/KMN/rate_SD_INT.dat", "", "");


  boost::filesystem::create_directory("../data/ph_emission/rph/Bs_extr/JPZ");

  Print_To_File({}, {xg_scan, FV_plot_JPZ.ave(), FV_plot_JPZ.err()}, "../data/ph_emission/rph/Bs_extr/JPZ/FV.dat", "", "");
  Print_To_File({}, {xg_scan, FA_plot_JPZ.ave(), FA_plot_JPZ.err()}, "../data/ph_emission/rph/Bs_extr/JPZ/FA.dat", "", "");
  Print_To_File({}, {xg_scan, FTV_plot_JPZ.ave(), FTV_plot_JPZ.err()}, "../data/ph_emission/rph/Bs_extr/JPZ/FTV.dat", "", "");
  Print_To_File({}, {xg_scan, FTA_plot_JPZ.ave(), FTA_plot_JPZ.err()}, "../data/ph_emission/rph/Bs_extr/JPZ/FTA.dat", "", "");

  Print_To_File({}, {qmin_list, xg_min_list, xg_max_list, rate_JPZ.ave(), rate_JPZ.err()}, "../data/ph_emission/rph/Bs_extr/JPZ/rate_SD_INT.dat", "", "");

 
  }


  //###########################################################################


  distr_t_list rate_SD_PB(false,Nxgs,Nb), rate_INT_PB(false,Nxgs,Nb), rate_diff_SD_PB(false,Nxgs,Nb), rate_diff_INT_PB(false,Nxgs,Nb);
  distr_t_list rate_diff_PT_PB(false,Nxgs,Nb), rate_SD_FV_only_PB(false,Nxgs,Nb);


  cout<<"Starting calculation of the rate! "<<endl;

    #pragma omp parallel for schedule(dynamic)
    for(int it=0; it < Nxgs;it++) {

      for(int jboot=0;jboot<Nb;jboot++) {

	 int JK= jboot;

	 double a1_val= a1_distr.distr[JK];
	 double mb_val= mb_distr.distr[JK];
	 double xb_val= xb_distr.distr[JK];
	 double Vtbs_val= Vtbs_distr.distr[JK];
	 double tbs_val = tbs_distr.distr[JK];
	 Vfloat kappa_list = kappa_distr[JK];
	 Vfloat Th_list = Th_distr[JK];
	 Vfloat W_list = W_distr[JK];
	 Vfloat Br_list = Br_distr[JK];
	 int R= Rbo[jboot];
	 double pr_s_re= pretty_s_RE.ave() +  (pretty_s_RE.distr[R]- pretty_s_RE.ave())*sqrt(NJ-1.0);
	 double pr_s_im= pretty_s_IM.ave() +(pretty_s_IM.distr[R]- pretty_s_IM.ave())*sqrt(NJ-1.0);
	 double f = f_Bs.ave() + (f_Bs.distr[R]-f_Bs.ave())*sqrt(NJ-1.0);
	 auto FA_B = [ &FA_boot, &R](double x) { return FA_boot(x,R);};
	 auto FV_B = [ &FV_boot, &R](double x) { return FV_boot(x,R);};
	 auto FTA_B = [ &FTA_boot, &pretty_b_boot, &pr_s_re, &R, &evolutor_ZT](double x) { return evolutor_ZT*(FTA_boot(x,R) + pretty_b_boot(x,R) + pr_s_re)  ;};
	 auto FTV_B = [ &FTV_boot, &pretty_b_boot, &pr_s_re, &R, &evolutor_ZT](double x) { return evolutor_ZT*(FTV_boot(x,R) + pretty_b_boot(x,R) + pr_s_re)  ;};
	 auto FT_IM_B = [ &pr_s_im, &evolutor_ZT](double x) { return evolutor_ZT*pr_s_im;};
	 auto zero_F = [](double x) { return 0.0;};

	 rate_SD_PB.distr_list[it].distr[JK] =   Compute_Bs_mumugamma_decay_rate( FV_B, FA_B, FTV_B, FTA_B, FT_IM_B, FT_IM_B , f , xb_val, Vtbs_val, tbs_val, kappa_list, Th_list, W_list, Br_list, "SD",  "CH",   xg_max_list[it], xg_min_list[it] ); 
	 rate_INT_PB.distr_list[it].distr[JK] =  Compute_Bs_mumugamma_decay_rate( FV_B, FA_B, FTV_B, FTA_B, FT_IM_B, FT_IM_B , f , xb_val, Vtbs_val, tbs_val, kappa_list, Th_list, W_list, Br_list, "INT",  "CH",   xg_max_list[it], xg_min_list[it] );
	 rate_diff_SD_PB.distr_list[it].distr[JK] = Compute_Bs_mumugamma_differential_decay_rate(xg_max_list[it], FV_B, FA_B, FTV_B, FTA_B, FT_IM_B, FT_IM_B, f ,  xb_val, Vtbs_val, tbs_val, kappa_list, Th_list, W_list, Br_list, "SD", "CH");
	 rate_diff_INT_PB.distr_list[it].distr[JK] = Compute_Bs_mumugamma_differential_decay_rate(xg_max_list[it], FV_B, FA_B, FTV_B, FTA_B, FT_IM_B, FT_IM_B, f ,  xb_val, Vtbs_val, tbs_val, kappa_list, Th_list, W_list, Br_list, "INT", "CH");
	 rate_diff_PT_PB.distr_list[it].distr[JK] = Compute_Bs_mumugamma_differential_decay_rate(xg_max_list[it], FV_B, FA_B, FTV_B, FTA_B, FT_IM_B, FT_IM_B, f ,  xb_val, Vtbs_val, tbs_val, kappa_list, Th_list, W_list, Br_list, "PT", "NO_CH");
	 rate_SD_FV_only_PB.distr_list[it].distr[JK] =  Compute_Bs_mumugamma_decay_rate( FV_B, FA_B, zero_F, zero_F, zero_F, zero_F , f , xb_val, Vtbs_val, tbs_val, kappa_list, Th_list, W_list, Br_list, "SD",  "CH",   xg_max_list[it], xg_min_list[it] );


	 if(jboot==0) {

	    for(int ijack=0;ijack<NJ;ijack++) {


	      auto FV = [&FV_interpol, &ijack](double x)  { return FV_interpol[ijack](x) ;};
	      auto FA = [&FA_interpol, &ijack](double x)  { return FA_interpol[ijack](x) ;};
	      auto FV_T = [&FTV_interpol, &pretty_b_interpol, &pretty_s_RE,  &evolutor_ZT, &ijack](double x) { return evolutor_ZT*(FTV_interpol[ijack](x) + pretty_b_interpol[ijack](x) +pretty_s_RE.distr[ijack] );};
	      auto FA_T = [&FTA_interpol, &pretty_b_interpol, &pretty_s_RE, &evolutor_ZT, &ijack](double x) { return evolutor_ZT*(FTA_interpol[ijack](x) + pretty_b_interpol[ijack](x)  +pretty_s_RE.distr[ijack] );};
	      auto F_T_IM = [&pretty_s_IM, &ijack, &evolutor_ZT](double x) { return evolutor_ZT*pretty_s_IM.distr[ijack];};
	      auto FV_T_eff_std = [&FV_T, &a1_distr_std, &f_Bs, &mb_distr_std, &ijack](double x) { return FV_T(x) +  0.0*(16.0/(3.0*WC(7,0).first))*(a1_distr_std.distr[ijack]*f_Bs.distr[ijack]/mb_distr_std.distr[ijack]);};
	      
	      
	      rate_Bs_SD_NOCH.distr_list[it].distr.push_back(  Compute_Bs_mumugamma_decay_rate( FV, FA, FV_T_eff_std, FA_T, F_T_IM, F_T_IM , f_Bs.distr[ijack] , xb_distr_std.distr[ijack], Vtbs_distr_std.distr[ijack], tbs_distr_std.distr[ijack], kappa_list, Th_list, W_list, Br_list , "SD", "NO_CH",  xg_max_list[it], xg_min_list[it] ));
	  
	      rate_Bs_INT_NOCH.distr_list[it].distr.push_back(  Compute_Bs_mumugamma_decay_rate(FV, FA, FV_T_eff_std, FA_T, F_T_IM, F_T_IM , f_Bs.distr[ijack] , xb_distr_std.distr[ijack], Vtbs_distr_std.distr[ijack], tbs_distr_std.distr[ijack], kappa_list, Th_list, W_list, Br_list , "INT", "NO_CH", xg_max_list[it], xg_min_list[it] ));
	
	      rate_Bs_diff_INT_NOCH.distr_list[it].distr.push_back( Compute_Bs_mumugamma_differential_decay_rate(xg_max_list[it], FV, FA, FV_T_eff_std, FA_T, F_T_IM, F_T_IM, f_Bs.distr[ijack] , xb_distr_std.distr[ijack], Vtbs_distr_std.distr[ijack], tbs_distr_std.distr[ijack], kappa_list, Th_list, W_list, Br_list , "INT", "NO_CH"));
	      
	      rate_Bs_diff_SD_NOCH.distr_list[it].distr.push_back( Compute_Bs_mumugamma_differential_decay_rate( xg_max_list[it], FV, FA, FV_T_eff_std, FA_T, F_T_IM, F_T_IM, f_Bs.distr[ijack] , xb_distr_std.distr[ijack], Vtbs_distr_std.distr[ijack], tbs_distr_std.distr[ijack], kappa_list, Th_list, W_list, Br_list , "SD" , "NO_CH"));

	      
	    }
	 }
      }

      cout<<"Finished calculation using bootstrap resampling of FF for ixg: "<<it<<endl<<flush;
      cout<<"xg: "<<xg_max_list[it]<<" SD+INT: "<<(rate_SD_PB+rate_INT_PB).ave(it)<<" +- "<<(rate_SD_PB+rate_INT_PB).err(it)<<endl;
      cout<<"xg: "<<xg_max_list[it]<<" SD+INT: "<<(rate_diff_SD_PB+rate_diff_INT_PB).ave(it)<<" +- "<<(rate_diff_SD_PB+rate_diff_INT_PB).err(it)<<endl;

     

      

      /*
      
      for(int ijack=0;ijack<NJ;ijack++) {


	auto FV = [&FV_interpol, &ijack](double x)  { return FV_interpol[ijack](x) ;};
	auto FA = [&FA_interpol, &ijack](double x)  { return FA_interpol[ijack](x) ;};
	auto FV_T = [&FTV_interpol, &pretty_b_interpol, &pretty_s_RE,  &evolutor_ZT, &ijack](double x) { return evolutor_ZT*(FTV_interpol[ijack](x) + pretty_b_interpol[ijack](x) +pretty_s_RE.distr[ijack] );};
	auto FA_T = [&FTA_interpol, &pretty_b_interpol, &pretty_s_RE, &evolutor_ZT, &ijack](double x) { return evolutor_ZT*(FTA_interpol[ijack](x) + pretty_b_interpol[ijack](x)  +pretty_s_RE.distr[ijack] );};
	auto F_T_IM = [&pretty_s_IM, &ijack, &evolutor_ZT](double x) { return evolutor_ZT*pretty_s_IM.distr[ijack];};
	
    

	for(int iboot=0;iboot<Nb;iboot++) {

	  int K= iboot;

	  double a1_val= a1_distr.distr[K];
	  double mb_val= mb_distr.distr[K];
	  double xb_val= xb_distr.distr[K];
	  double Vtbs_val= Vtbs_distr.distr[K];
	  double tbs_val = tbs_distr.distr[K];
	  Vfloat kappa_list = kappa_distr[K];

          Vfloat Th_list = Th_distr[K];
	  Vfloat W_list = W_distr[K];
	  Vfloat Br_list = Br_distr[K];




	  
	  auto FV_T_eff = [&FV_T, &a1_val, &f_Bs, &mb_val, &ijack](double x) { return FV_T(x) +  0.0*(16.0/(3.0*WC(7,0).first))*(a1_val*f_Bs.distr[ijack]/mb_val);};
	  auto FV_T_eff_std = [&FV_T, &a1_distr_std, &f_Bs, &mb_distr_std, &ijack](double x) { return FV_T(x) +  0.0*(16.0/(3.0*WC(7,0).first))*(a1_distr_std.distr[ijack]*f_Bs.distr[ijack]/mb_distr_std.distr[ijack]);};
	 

	

	  auto zero_FA = [](double x) { return 0.0;};
	  auto zero_FA_T= [](double x) { return 0.0;};
	  auto zero_F_T_IM= [](double x) { return 0.0;};
	  auto zero_FV_T_eff=  [&a1_val, &f_Bs, &mb_val, &ijack](double x) { return 0.0;}; 

	 
	  
	  rate_Bs_SD_boot[it].distr_list[iboot].distr.push_back(  Compute_Bs_mumugamma_decay_rate( FV, FA, FV_T_eff, FA_T, F_T_IM, F_T_IM , f_Bs.distr[ijack] , xb_val, Vtbs_val, tbs_val, kappa_list, Th_list, W_list, Br_list, "SD",  "CH",   xg_max_list[it], xg_min_list[it] ));
	  
	  rate_Bs_SD_FV_only_boot[it].distr_list[iboot].distr.push_back(  Compute_Bs_mumugamma_decay_rate( FV, FA, zero_FV_T_eff, zero_FA_T, zero_F_T_IM, zero_F_T_IM , f_Bs.distr[ijack] , xb_val, Vtbs_val, tbs_val, kappa_list, Th_list, W_list, Br_list, "SD",  "NO_CH",   xg_max_list[it], xg_min_list[it] ));
	  
	  
	  
	  rate_Bs_INT_boot[it].distr_list[iboot].distr.push_back(  Compute_Bs_mumugamma_decay_rate(FV, FA, FV_T_eff, FA_T, F_T_IM, F_T_IM , f_Bs.distr[ijack] ,  xb_val, Vtbs_val, tbs_val, kappa_list, Th_list, W_list, Br_list, "INT","CH", xg_max_list[it], xg_min_list[it] ));
	  

	  
	  rate_Bs_diff_PT_boot[it].distr_list[iboot].distr.push_back( Compute_Bs_mumugamma_differential_decay_rate(xg_max_list[it], FV, FA, FV_T_eff, FA_T, F_T_IM, F_T_IM, f_Bs.distr[ijack] ,  xb_val, Vtbs_val, tbs_val, kappa_list, Th_list, W_list, Br_list, "PT", "NO_CH"));


	  
	  rate_Bs_diff_INT_boot[it].distr_list[iboot].distr.push_back( Compute_Bs_mumugamma_differential_decay_rate(xg_max_list[it], FV, FA, FV_T_eff, FA_T, F_T_IM, F_T_IM, f_Bs.distr[ijack] , xb_val, Vtbs_val, tbs_val, kappa_list, Th_list, W_list, Br_list , "INT", "CH"));


	  
	  rate_Bs_diff_SD_boot[it].distr_list[iboot].distr.push_back( Compute_Bs_mumugamma_differential_decay_rate(xg_max_list[it], FV, FA, FV_T_eff, FA_T, F_T_IM, F_T_IM, f_Bs.distr[ijack] ,  xb_val, Vtbs_val, tbs_val, kappa_list, Th_list, W_list, Br_list  ,  "SD" , "CH"));


	  
	  AFB_Bs_boot[it].distr_list[iboot].distr.push_back( Compute_AFB( xg_max_list[it], FV, FA, FV_T_eff, FA_T, F_T_IM, F_T_IM, f_Bs.distr[ijack] ,  xb_val, Vtbs_val, tbs_val, kappa_list, Th_list, W_list, Br_list  ,  "SD", "CH" ));


	  
	  


	  if(iboot==0) {
	  

	  
	  rate_Bs_SD_NOCH.distr_list[it].distr.push_back(  Compute_Bs_mumugamma_decay_rate( FV, FA, FV_T_eff_std, FA_T, F_T_IM, F_T_IM , f_Bs.distr[ijack] , xb_distr_std.distr[ijack], Vtbs_distr_std.distr[ijack], tbs_distr_std.distr[ijack], kappa_list, Th_list, W_list, Br_list , "SD", "NO_CH",  xg_max_list[it], xg_min_list[it] ));
	  
	  rate_Bs_INT_NOCH.distr_list[it].distr.push_back(  Compute_Bs_mumugamma_decay_rate(FV, FA, FV_T_eff_std, FA_T, F_T_IM, F_T_IM , f_Bs.distr[ijack] , xb_distr_std.distr[ijack], Vtbs_distr_std.distr[ijack], tbs_distr_std.distr[ijack], kappa_list, Th_list, W_list, Br_list , "INT", "NO_CH", xg_max_list[it], xg_min_list[it] ));
	  
	  rate_Bs_diff_INT_NOCH.distr_list[it].distr.push_back( Compute_Bs_mumugamma_differential_decay_rate(xg_max_list[it], FV, FA, FV_T_eff_std, FA_T, F_T_IM, F_T_IM, f_Bs.distr[ijack] , xb_distr_std.distr[ijack], Vtbs_distr_std.distr[ijack], tbs_distr_std.distr[ijack], kappa_list, Th_list, W_list, Br_list , "INT", "NO_CH"));
	  
	  rate_Bs_diff_SD_NOCH.distr_list[it].distr.push_back( Compute_Bs_mumugamma_differential_decay_rate( xg_max_list[it], FV, FA, FV_T_eff_std, FA_T, F_T_IM, F_T_IM, f_Bs.distr[ijack] , xb_distr_std.distr[ijack], Vtbs_distr_std.distr[ijack], tbs_distr_std.distr[ijack], kappa_list, Th_list, W_list, Br_list , "SD" , "NO_CH"));

	

	  }
	  
	   

	}

	
	
      }

      cout<<"Computed rate for xg: "<<xg_max_list[it]<<endl;

      //combine jacknife and boostrap

      

      
      rate_Bs_SD.distr_list[it] = BMA_Eq_29( rate_Bs_SD_boot[it].distr_list);
      rate_Bs_SD_FV_only.distr_list[it] = BMA_Eq_29( rate_Bs_SD_FV_only_boot[it].distr_list);
      rate_Bs_INT.distr_list[it] = BMA_Eq_29( rate_Bs_INT_boot[it].distr_list);
      rate_Bs_diff_PT.distr_list[it] = BMA_Eq_29( rate_Bs_diff_PT_boot[it].distr_list);
      rate_Bs_diff_INT.distr_list[it] = BMA_Eq_29( rate_Bs_diff_INT_boot[it].distr_list);
      rate_Bs_diff_SD.distr_list[it] = BMA_Eq_29( rate_Bs_diff_SD_boot[it].distr_list);
      AFB_Bs.distr_list[it] = BMA_Eq_29( AFB_Bs_boot[it].distr_list);
      */
    }


    /*
  
    
    //Print table on screen
    cout<<" & " <<"\\multicolumn{7}{c}{$\\sqrt{q^{2}_{\\rm cut}} [\\rm{GeV}]$} \\hline"<<endl;  
    cout<<" & $4.0$ & $4.1$ & $4.2$ & $4.3$ & $4.4$ & $4.5$ & $4.6$ \\hline"<<endl;
    cout<<" $\\mathcal{B}_{\\rm SD+INT}$  ";
    for(int i=0;i<7;i++) cout<<"& $ "<<(rate_Bs_SD + rate_Bs_INT).distr_list[4*i].ave()<<"("<<(rate_Bs_SD + rate_Bs_INT).distr_list[4*i].err()<<")$ ";
    cout<<"\\hline"<<endl;
    cout<<" & $4.7$ & $4.8$ & $4.9$ & $5.0$ & $5.1$ & $5.2$ & $5.3$ \\hline"<<endl;
    cout<<" $\\mathcal{B}_{\\rm SD+INT}$  ";
    for(int i=7;i<14;i++) cout<<"& $ "<<(rate_Bs_SD + rate_Bs_INT).distr_list[4*i].ave()<<"("<<(rate_Bs_SD + rate_Bs_INT).distr_list[4*i].err()<<")$ ";
    cout<<"\\hline"<<endl;

    cout<<"printing actual q_max: "<<endl;
    for(int i=0;i<14;i++) cout<<qmin_list[i*4]<<endl;
    
    //Bs
    Print_To_File( {}, { qmin_list, xg_min_list, xg_max_list, rate_Bs_SD.ave(), rate_Bs_SD.err(), rate_Bs_SD_NOCH.ave(), rate_Bs_SD_NOCH.err(), rate_Bs_SD_FV_only.ave(), rate_Bs_SD_FV_only.err() }, "../data/ph_emission/"+ph_type+"/Bs_extr/rate_Bs_SD.dat", "", "#qmin[GeV] xg_min xg_max  rate ");
    Print_To_File( {}, { qmin_list, xg_min_list, xg_max_list, rate_Bs_INT.ave(), rate_Bs_INT.err(), rate_Bs_INT_NOCH.ave(), rate_Bs_INT_NOCH.err()}, "../data/ph_emission/"+ph_type+"/Bs_extr/rate_Bs_INT.dat", "", "#qmin[GeV] xg_min xg_max  rate");
    Print_To_File( {}, { qmin_list, xg_min_list, xg_max_list, (rate_Bs_INT+rate_Bs_SD).ave(), (rate_Bs_INT+rate_Bs_SD).err(), (rate_Bs_INT_NOCH+rate_Bs_SD_NOCH).ave(), (rate_Bs_INT_NOCH+rate_Bs_SD_NOCH).err() }, "../data/ph_emission/"+ph_type+"/Bs_extr/rate_Bs_SDINT.dat", "", "#qmin[GeV] xg_min xg_max  rate");
    
    
    
    Print_To_File( {}, { qmin_list, xg_max_list, rate_Bs_diff_PT.ave(), rate_Bs_diff_PT.err() }, "../data/ph_emission/"+ph_type+"/Bs_extr/rate_diff_Bs_PT.dat", "", "#qmin[GeV] xg_max  rate ");
    Print_To_File( {}, { qmin_list, xg_max_list, rate_Bs_diff_INT.ave(), rate_Bs_diff_INT.err(), rate_Bs_diff_INT_NOCH.ave(), rate_Bs_diff_INT_NOCH.err() }, "../data/ph_emission/"+ph_type+"/Bs_extr/rate_diff_Bs_INT.dat", "", "#qmin[GeV] xg_max  rate");
    Print_To_File( {}, { qmin_list, xg_max_list, rate_Bs_diff_SD.ave(), rate_Bs_diff_SD.err(), rate_Bs_diff_SD_NOCH.ave(), rate_Bs_diff_SD_NOCH.err()}, "../data/ph_emission/"+ph_type+"/Bs_extr/rate_diff_Bs_SD.dat", "", "#qmin[GeV] xg_max  rate");
    Print_To_File( {}, { qmin_list, xg_max_list, AFB_Bs.ave(), AFB_Bs.err()},  "../data/ph_emission/"+ph_type+"/Bs_extr/AFB_Bs.dat", "", "#q[GeV] xg   AFG    ");

    */

    
    //PRINT PURE BOOTSTRAP RATE
    Print_To_File( {}, { qmin_list, xg_min_list, xg_max_list, rate_SD_PB.ave(), rate_SD_PB.err(), rate_Bs_SD_NOCH.ave(), rate_Bs_SD_NOCH.err(), rate_SD_FV_only_PB.ave(), rate_SD_FV_only_PB.err()  }, "../data/ph_emission/"+ph_type+"/Bs_extr/rate_Bs_SD_PB.dat", "", "#qmin[GeV] xg_min xg_max  rate ");
    Print_To_File( {}, { qmin_list, xg_min_list, xg_max_list, rate_INT_PB.ave(), rate_INT_PB.err(),  rate_Bs_INT_NOCH.ave(), rate_Bs_INT_NOCH.err()}, "../data/ph_emission/"+ph_type+"/Bs_extr/rate_Bs_INT_PB.dat", "", "#qmin[GeV] xg_min xg_max  rate");
    Print_To_File( {}, { qmin_list, xg_min_list, xg_max_list, (rate_INT_PB+rate_SD_PB).ave(), (rate_INT_PB+rate_SD_PB).err(),  (rate_Bs_INT_NOCH+rate_Bs_SD_NOCH).ave(), (rate_Bs_INT_NOCH+rate_Bs_SD_NOCH).err() }, "../data/ph_emission/"+ph_type+"/Bs_extr/rate_Bs_SDINT_PB.dat", "", "#qmin[GeV] xg_min xg_max  rate");
    
    Print_To_File( {}, { qmin_list, xg_max_list, rate_diff_INT_PB.ave(), rate_diff_INT_PB.err(), rate_Bs_diff_INT_NOCH.ave(), rate_Bs_diff_INT_NOCH.err()  }, "../data/ph_emission/"+ph_type+"/Bs_extr/rate_diff_Bs_INT_PB.dat", "", "#qmin[GeV] xg_max  rate");
    Print_To_File( {}, { qmin_list, xg_max_list, rate_diff_SD_PB.ave(), rate_diff_SD_PB.err(),  rate_Bs_diff_SD_NOCH.ave(), rate_Bs_diff_SD_NOCH.err()}, "../data/ph_emission/"+ph_type+"/Bs_extr/rate_diff_Bs_SD_PB.dat", "", "#qmin[GeV] xg_max  rate");
    Print_To_File( {}, { qmin_list, xg_max_list, rate_diff_PT_PB.ave(), rate_diff_PT_PB.err() }, "../data/ph_emission/"+ph_type+"/Bs_extr/rate_diff_Bs_PT_PB.dat", "", "#qmin[GeV] xg_max  rate ");


    //Print table on screen
    cout<<" & " <<"\\multicolumn{7}{c}{$\\sqrt{q^{2}_{\\rm cut}} [\\rm{GeV}]$} \\hline"<<endl;  
    cout<<" & $4.0$ & $4.1$ & $4.2$ & $4.3$ & $4.4$ & $4.5$ & $4.6$ \\hline"<<endl;
    cout<<" $\\mathcal{B}_{\\rm SD+INT}$  ";
    for(int i=0;i<7;i++) cout<<"& $ "<<(rate_SD_PB + rate_INT_PB).distr_list[4*i].ave()<<"("<<(rate_SD_PB + rate_INT_PB).distr_list[4*i].err()<<")$ ";
    cout<<"\\hline"<<endl;
    cout<<" & $4.7$ & $4.8$ & $4.9$ & $5.0$ & $5.1$ & $5.2$ & $5.3$ \\hline"<<endl;
    cout<<" $\\mathcal{B}_{\\rm SD+INT}$  ";
    for(int i=7;i<14;i++) cout<<"& $ "<<(rate_SD_PB + rate_INT_PB).distr_list[4*i].ave()<<"("<<(rate_SD_PB + rate_INT_PB).distr_list[4*i].err()<<")$ ";
    cout<<"\\hline"<<endl;

    cout<<"printing actual q_max: "<<endl;
    for(int i=0;i<14;i++) cout<<qmin_list[i*4]<<endl;
    
     
    cout<<"Done! Bye"<<endl;

  }
    
    return;
}


rt_FF_Bs Get_Bs_mumu_gamma_form_factors(int num_xg, int Perform_continuum_extrapolation, bool Use_three_finest, bool Include_a4, bool UseJack, string Fit_tag, string path_list, string out_tag) {


  cout<<"CALLING Get_Bs_mumu_gamma_form_factors"<<endl;

  int size_mu_nu= Is_reph?2:4;

  int size_mu_nu_T=2;

  //return type
  rt_FF_Bs rt_fit(UseJack);
  rt_fit.Use_three_finest=Use_three_finest;
  rt_fit.Include_a4 = Include_a4;
  rt_fit.num_xg=num_xg;

  
  
  string _Fit_tag = "";

  if(Fit_tag.size() > 0) { _Fit_tag= "_"+Fit_tag.substr(0, Fit_tag.size()-1); } 

  const int Njacks=((UseJack==true)?NJ:Nboots);

   
  string ph_type_mes=ph_type+"/"+MESON;
  
  




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

  
    if(A_bis.length() <= 4) return A_bis < B_bis;
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
    return A_bis<B_bis;
  };
  
  
 
  int off_i = (Is_reph?1:0);
  int off_T = 1;
  


  //2pts
  data_t data_2pts, data_2pts_SM, data_2pts_SMSM;
  data_t data_2pts_V1, data_2pts_V2, data_2pts_V3;
  data_t data_2pts_u_V1, data_2pts_u_V2, data_2pts_u_V3;
  data_t data_2pts_d_V1, data_2pts_d_V2, data_2pts_d_V3;
  data_t data_2pts_u_T1, data_2pts_u_T2, data_2pts_u_T3;
  data_t data_2pts_u_T1_SMSM, data_2pts_u_T2_SMSM, data_2pts_u_T3_SMSM;
  data_t data_2pts_u_B1, data_2pts_u_B2, data_2pts_u_B3;

  
  //Read data

  //2pts function

  data_2pts.Read(path_list, "mes_contr_2pts_SM_3", "P5P5", Sort_light_confs);
  data_2pts_SM.Read(path_list, "mes_contr_2pts_SM_3", "P5P5", Sort_light_confs);
  data_2pts_SMSM.Read(path_list, "mes_contr_2pts_SMSM_3", "P5P5", Sort_light_confs);
  data_2pts_V1.Read(path_list, "mes_contr_2pts_SM_3", "V1V1", Sort_light_confs);
  data_2pts_V2.Read(path_list, "mes_contr_2pts_SM_3", "V2V2", Sort_light_confs);
  data_2pts_V3.Read(path_list, "mes_contr_2pts_SM_3", "V3V3", Sort_light_confs);
  data_2pts_d_V1.Read(path_list, "mes_contr_2pts_SM_1", "V1V1", Sort_light_confs);
  data_2pts_d_V2.Read(path_list, "mes_contr_2pts_SM_1", "V2V2", Sort_light_confs);
  data_2pts_d_V3.Read(path_list, "mes_contr_2pts_SM_1", "V3V3", Sort_light_confs);
  data_2pts_u_V1.Read(path_list, "mes_contr_2pts_SM_2", "V1V1" , Sort_light_confs);
  data_2pts_u_V2.Read(path_list, "mes_contr_2pts_SM_2", "V2V2", Sort_light_confs);
  data_2pts_u_V3.Read(path_list, "mes_contr_2pts_SM_2", "V3V3", Sort_light_confs);
  data_2pts_u_T1.Read(path_list, "mes_contr_2pts_SM_2", "T1T1" , Sort_light_confs);
  data_2pts_u_T2.Read(path_list, "mes_contr_2pts_SM_2", "T2T2", Sort_light_confs);
  data_2pts_u_T3.Read(path_list, "mes_contr_2pts_SM_2", "T3T3", Sort_light_confs);
  data_2pts_u_T1_SMSM.Read(path_list, "mes_contr_2pts_SMSM_2", "T1T1" , Sort_light_confs);
  data_2pts_u_T2_SMSM.Read(path_list, "mes_contr_2pts_SMSM_2", "T2T2", Sort_light_confs);
  data_2pts_u_T3_SMSM.Read(path_list, "mes_contr_2pts_SMSM_2", "T3T3", Sort_light_confs);
  data_2pts_u_B1.Read(path_list, "mes_contr_2pts_SM_2", "B1B1" , Sort_light_confs);
  data_2pts_u_B2.Read(path_list, "mes_contr_2pts_SM_2", "B2B2", Sort_light_confs);
  data_2pts_u_B3.Read(path_list, "mes_contr_2pts_SM_2", "B3B3", Sort_light_confs);



  //axial
  vector<vector<vector<data_t>>> C_A_Fu_data(size_mu_nu), C_A_Fd_data(size_mu_nu), C_A_Bu_data(size_mu_nu), C_A_Bd_data(size_mu_nu);
  //vector
  vector<vector<vector<data_t>>> C_V_Fu_data(size_mu_nu), C_V_Fd_data(size_mu_nu), C_V_Bu_data(size_mu_nu), C_V_Bd_data(size_mu_nu);
  //TK
  vector<vector<vector<data_t>>> C_TK_Fu_data(size_mu_nu), C_TK_Fd_data(size_mu_nu), C_TK_Bu_data(size_mu_nu), C_TK_Bd_data(size_mu_nu);
  //BK
  vector<vector<vector<data_t>>> C_BK_Fu_data(size_mu_nu), C_BK_Fd_data(size_mu_nu), C_BK_Bu_data(size_mu_nu), C_BK_Bd_data(size_mu_nu);

  for(int mu=0;mu<size_mu_nu;mu++) {

    C_A_Fu_data[mu].resize(size_mu_nu);
    C_A_Fd_data[mu].resize(size_mu_nu);
    C_A_Bu_data[mu].resize(size_mu_nu);
    C_A_Bd_data[mu].resize(size_mu_nu);
    C_V_Fu_data[mu].resize(size_mu_nu);
    C_V_Fd_data[mu].resize(size_mu_nu);
    C_V_Bu_data[mu].resize(size_mu_nu);
    C_V_Bd_data[mu].resize(size_mu_nu);

    C_BK_Fu_data[mu].resize(size_mu_nu_T);
    C_BK_Fd_data[mu].resize(size_mu_nu_T);
    C_BK_Bu_data[mu].resize(size_mu_nu_T);
    C_BK_Bd_data[mu].resize(size_mu_nu_T);
    C_TK_Fu_data[mu].resize(size_mu_nu_T);
    C_TK_Fd_data[mu].resize(size_mu_nu_T);
    C_TK_Bu_data[mu].resize(size_mu_nu_T);
    C_TK_Bd_data[mu].resize(size_mu_nu_T);

    
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

    for(int nu=0;nu<size_mu_nu_T;nu++) {

      C_BK_Fu_data[mu][nu].resize(num_xg);
      C_BK_Fd_data[mu][nu].resize(num_xg);
      C_BK_Bu_data[mu][nu].resize(num_xg);
      C_BK_Bd_data[mu][nu].resize(num_xg);
      C_TK_Fu_data[mu][nu].resize(num_xg);
      C_TK_Fd_data[mu][nu].resize(num_xg);
      C_TK_Bu_data[mu][nu].resize(num_xg);
      C_TK_Bd_data[mu][nu].resize(num_xg);

    }

  }

  //read data
  for(int ixg=0;ixg<num_xg;ixg++) {
    for(int mu=0;mu<size_mu_nu;mu++) {
      for(int nu=0;nu<size_mu_nu;nu++) {

	//vector
	//Fu
	C_V_Fu_data[mu][nu][ixg].Read(path_list, "C_mu_"+to_string(mu+off_i)+"_FF_u_ixg_"+to_string(ixg), "V"+to_string(nu+off_i)+"P5", Sort_light_confs );
	//Fd
	C_V_Fd_data[mu][nu][ixg].Read(path_list, "C_mu_"+to_string(mu+off_i)+"_FF_d_ixg_"+to_string(ixg), "V"+to_string(nu+off_i)+"P5", Sort_light_confs );
	//Bu
	C_V_Bu_data[mu][nu][ixg].Read(path_list, "C_mu_"+to_string(mu+off_i)+"_BB_u_ixg_"+to_string(ixg), "V"+to_string(nu+off_i)+"P5", Sort_light_confs );
	//Bd
	C_V_Bd_data[mu][nu][ixg].Read(path_list, "C_mu_"+to_string(mu+off_i)+"_BB_d_ixg_"+to_string(ixg), "V"+to_string(nu+off_i)+"P5", Sort_light_confs );

	//axial
	//Fu
	C_A_Fu_data[mu][nu][ixg].Read(path_list, "C_mu_"+to_string(mu+off_i)+"_FF_u_ixg_"+to_string(ixg), "A"+to_string(nu+off_i)+"P5", Sort_light_confs );
	//Fd
	C_A_Fd_data[mu][nu][ixg].Read(path_list, "C_mu_"+to_string(mu+off_i)+"_FF_d_ixg_"+to_string(ixg), "A"+to_string(nu+off_i)+"P5", Sort_light_confs );
	//Bu
	C_A_Bu_data[mu][nu][ixg].Read(path_list, "C_mu_"+to_string(mu+off_i)+"_BB_u_ixg_"+to_string(ixg), "A"+to_string(nu+off_i)+"P5", Sort_light_confs );
	//Bd
	C_A_Bd_data[mu][nu][ixg].Read(path_list, "C_mu_"+to_string(mu+off_i)+"_BB_d_ixg_"+to_string(ixg), "A"+to_string(nu+off_i)+"P5", Sort_light_confs );
      }
      for(int nu=0;nu<size_mu_nu_T;nu++) {

	//tensor electric part
	//Fu
	C_TK_Fu_data[mu][nu][ixg].Read(path_list, "C_mu_"+to_string(mu+off_i)+"_FF_u_ixg_"+to_string(ixg), "T"+to_string(nu+off_T)+"P5", Sort_light_confs );
	//Fd
	C_TK_Fd_data[mu][nu][ixg].Read(path_list, "C_mu_"+to_string(mu+off_i)+"_FF_d_ixg_"+to_string(ixg), "T"+to_string(nu+off_T)+"P5", Sort_light_confs );
	//Bu
	C_TK_Bu_data[mu][nu][ixg].Read(path_list, "C_mu_"+to_string(mu+off_i)+"_BB_u_ixg_"+to_string(ixg), "T"+to_string(nu+off_T)+"P5", Sort_light_confs );
	//Bd
	C_TK_Bd_data[mu][nu][ixg].Read(path_list, "C_mu_"+to_string(mu+off_i)+"_BB_d_ixg_"+to_string(ixg), "T"+to_string(nu+off_T)+"P5", Sort_light_confs );
	//tensor magnetic part
	//Fu
	C_BK_Fu_data[mu][nu][ixg].Read(path_list, "C_mu_"+to_string(mu+off_i)+"_FF_u_ixg_"+to_string(ixg), "B"+to_string(nu+off_T)+"P5", Sort_light_confs );
	//Fd
	C_BK_Fd_data[mu][nu][ixg].Read(path_list, "C_mu_"+to_string(mu+off_i)+"_FF_d_ixg_"+to_string(ixg), "B"+to_string(nu+off_T)+"P5", Sort_light_confs );
	//Bu
	C_BK_Bu_data[mu][nu][ixg].Read(path_list, "C_mu_"+to_string(mu+off_i)+"_BB_u_ixg_"+to_string(ixg), "B"+to_string(nu+off_T)+"P5", Sort_light_confs );
	//Bd
	C_BK_Bd_data[mu][nu][ixg].Read(path_list, "C_mu_"+to_string(mu+off_i)+"_BB_d_ixg_"+to_string(ixg), "B"+to_string(nu+off_T)+"P5", Sort_light_confs );
      }
    }
  }

  //#######################################################################################################################################



 

  int Nens = data_2pts.size;
  GaussianMersenne GM(543543);
  GaussianMersenne GM_T(5432);

  //vectors to store FV and FA for all ensembles and values of gamma
  vector<distr_t_list> FA_per_ens(Nens), FV_per_ens(Nens), FA_u_per_ens(Nens), FV_u_per_ens(Nens), FA_d_per_ens(Nens), FV_d_per_ens(Nens);
  vector<distr_t_list> FA_per_kin(num_xg-1), FV_per_kin(num_xg-1), FA_u_per_kin(num_xg-1), FV_u_per_kin(num_xg-1), FA_d_per_kin(num_xg-1), FV_d_per_kin(num_xg-1); //num_xg -1 is to exclude xg=0 which is undefined
  vector<distr_t_list> FA_T_per_ens(Nens), FV_T_per_ens(Nens), FA_T_u_per_ens(Nens), FV_T_u_per_ens(Nens), FA_T_d_per_ens(Nens), FV_T_d_per_ens(Nens);
  vector<distr_t_list> FA_T_per_kin(num_xg-1), FV_T_per_kin(num_xg-1), FA_T_u_per_kin(num_xg-1), FV_T_u_per_kin(num_xg-1), FA_T_d_per_kin(num_xg-1), FV_T_d_per_kin(num_xg-1);
  //TK and BK component
  vector<distr_t_list> FB_per_ens(Nens), FT_per_ens(Nens), FB_u_per_ens(Nens), FT_u_per_ens(Nens), FB_d_per_ens(Nens), FT_d_per_ens(Nens);
  vector<distr_t_list> FB_per_kin(num_xg-1), FT_per_kin(num_xg-1), FB_u_per_kin(num_xg-1), FT_u_per_kin(num_xg-1), FB_d_per_kin(num_xg-1), FT_d_per_kin(num_xg-1);
  distr_t_list a_distr_list(UseJack);
  distr_t_list a_distr_list_red(UseJack);
  Vfloat L_list, L_list_red;
  //M_Ds, F_Ds
  distr_t_list MP_list(UseJack);
  distr_t_list MP_list_red(UseJack);
  distr_t_list FP_list(UseJack);
  distr_t_list FP_3pt_list(UseJack);
  distr_t_list FP_diml_list(UseJack);
  distr_t_list FP_3pt_diml_list(UseJack);
  distr_t_list ZT_list(UseJack);
  distr_t_list MP_ov_FP_list(UseJack);
  distr_t_list MP_ov_FP_3pt_list(UseJack);
  //resize
 

  Vfloat Tmin_fp, Tmax_fp;
  Vfloat Tmin_fp_3pt, Tmax_fp_3pt;

  
  
  
  
  //define lambda function to combine FF and BB

  //resample RCs
  distr_t ZA_A(UseJack), ZA_B(UseJack), ZA_C(UseJack), ZA_D(UseJack);
  distr_t ZV_A(UseJack), ZV_B(UseJack), ZV_C(UseJack), ZV_D(UseJack);
  distr_t ZT_A(UseJack), ZT_B(UseJack), ZT_C(UseJack), ZT_D(UseJack);
  distr_t a_A(UseJack), a_B(UseJack), a_C(UseJack), a_D(UseJack);

  //resample Jpsi tensor decay constant
  distr_t fT_Jpsi_distr(UseJack);

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

    fT_Jpsi_distr.distr.push_back( fT_Jpsi + fT_Jpsi_err*GM_T()/((UseJack==true)?sqrt(Njacks-1.0):1.0));

    
    ZT_A.distr.push_back( L_info_A.ZT_RI2 + GM_T()*L_info_A.ZT_RI2_err/((UseJack==true)?sqrt(Njacks -1.0):1.0));
    ZT_B.distr.push_back( L_info_B.ZT_RI2 + GM_T()*L_info_B.ZT_RI2_err/((UseJack==true)?sqrt(Njacks -1.0):1.0));
    ZT_C.distr.push_back( L_info_C.ZT_RI2 + GM_T()*L_info_C.ZT_RI2_err/((UseJack==true)?sqrt(Njacks -1.0):1.0));
    ZT_D.distr.push_back( L_info_D.ZT_RI2 + GM_T()*L_info_D.ZT_RI2_err/((UseJack==true)?sqrt(Njacks -1.0):1.0));
    

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

    boost::filesystem::create_directory("../data/ph_emission/"+ph_type_mes+"/C_T/"+data_2pts.Tag[iens]);
    boost::filesystem::create_directory("../data/ph_emission/"+ph_type_mes+"/H_T/"+data_2pts.Tag[iens]);
      
    
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

    thetas= Read_From_File(path_list+"/"+data_2pts.Tag[iens]+"/pars_list.dat", 1 , 5);
    virts=  Read_From_File(path_list+"/"+data_2pts.Tag[iens]+"/pars_list.dat", 2 , 5);
    masses_u= Read_From_File(path_list+"/"+data_2pts.Tag[iens]+"/pars_list.dat", 3 , 5);
    masses_d= Read_From_File(path_list+"/"+data_2pts.Tag[iens]+"/pars_list.dat", 4 , 5);

    cout<<"pars_list.dat: Read!"<<endl;

    //if((signed)thetas.size() != num_xg) crash("Number of rows in pars_list.dat does not match num_xg"); 

    int num_xg_per_ens= (signed)thetas.size();
       
    //RCs
    distr_t Za, Zv, Z_T, a_distr;
    if(data_2pts.Tag[iens].substr(1,1)=="A") {      Za= ZA_A; Zv=ZV_A; Z_T=ZT_A; a_distr=a_A;}
    else if(data_2pts.Tag[iens].substr(1,1)=="B") { Za= ZA_B; Zv=ZV_B; Z_T=ZT_B; a_distr=a_B;}
    else if(data_2pts.Tag[iens].substr(1,1)=="C") { Za= ZA_C; Zv=ZV_C; Z_T=ZT_C; a_distr=a_C;}
    else if(data_2pts.Tag[iens].substr(1,1)=="D") { Za= ZA_D; Zv=ZV_D; Z_T=ZT_D; a_distr=a_D;}
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
    cout<<"ZT: "<<Z_T.ave()<<" +- "<<Z_T.err()<<endl;
    cout<<"mu: "<<mu<<endl;
    cout<<"md: "<<md<<endl;
    cout<<"Number of kinematic present: "<<num_xg_per_ens<<endl;
    cout<<"Number of kinematic analyzed: "<<num_xg<<endl;
    
    
   
    
    //set time interval for eff_mass_fit SM
    if(data_2pts.Tag[iens].substr(1,1) =="A") {
      if(MESON == "B1s") { Corr.Tmin=14; Corr.Tmax= 24;}
      else if(MESON == "B2s" ) { Corr.Tmin=13; Corr.Tmax= 24;}
      else if(MESON=="B3s") { Corr.Tmin=14; Corr.Tmax=20;}
      else if(MESON=="B4s") { Corr.Tmin=14; Corr.Tmax=19;}
      else { Corr.Tmin=14; Corr.Tmax=27;}
      
    }
    else if(data_2pts.Tag[iens] =="cB211b.072.64") {
      if(MESON == "B1s") {Corr.Tmin=24; Corr.Tmax=38;}
      else if(MESON == "B2s") {Corr.Tmin=22; Corr.Tmax=34;}
      else if(MESON=="B3s") {Corr.Tmin=17; Corr.Tmax=23;}
      else if(MESON=="B4s") { Corr.Tmin=16; Corr.Tmax=22;}
      else {Corr.Tmin=23; Corr.Tmax=33;}
    }
    else if(data_2pts.Tag[iens] =="cB211b.072.96") {Corr.Tmin=20; Corr.Tmax=36;}
    
    else if(data_2pts.Tag[iens].substr(1,1) == "C")  {
      if(MESON == "B1s") {Corr.Tmin=27; Corr.Tmax=36;}
      else if(MESON == "B2s") {Corr.Tmin=23; Corr.Tmax=36;}
      else if(MESON=="B3s") {Corr.Tmin=18; Corr.Tmax=31;}
      else if(MESON=="B4s") { Corr.Tmin=18; Corr.Tmax=30;}
      else {Corr.Tmin=29; Corr.Tmax=42;}
    }
    else if(data_2pts.Tag[iens].substr(1,1) == "D")  {
      if(MESON == "B1s") {Corr.Tmin=32; Corr.Tmax=53;}
      else if(MESON == "B2s") { Corr.Tmin=28; Corr.Tmax= 37;}
      else if(MESON=="B3s") {Corr.Tmin=25; Corr.Tmax=37;}
      else if(MESON=="B4s") { Corr.Tmin=27; Corr.Tmax=37;}
      else {Corr.Tmin=38; Corr.Tmax=52;}
    }
    else crash("In fixing [Tmin, Tmax] for smeared MP, Ensemble: "+data_2pts.Tag[iens]+" not recognized");
    distr_t_list pt2_distr_SM= Corr.corr_t(data_2pts_SM.col(0)[iens], "../data/ph_emission/"+ph_type_mes+"/"+"C/"+data_2pts.Tag[iens]+"/corr_2pt_SM.dat");
    distr_t_list eff_mass_SM = Corr.effective_mass_t(pt2_distr_SM, "../data/ph_emission/"+ph_type_mes+"/"+"mass/"+data_2pts.Tag[iens]+"/eff_mass_SM.dat");
    distr_t_list pt2_distr_SMSM= Corr.corr_t(data_2pts_SMSM.col(0)[iens], "../data/ph_emission/"+ph_type_mes+"/"+"C/"+data_2pts.Tag[iens]+"/corr_2pt_SMSM.dat");
    distr_t_list eff_mass_SMSM = Corr.effective_mass_t(pt2_distr_SMSM, "../data/ph_emission/"+ph_type_mes+"/"+"mass/"+data_2pts.Tag[iens]+"/eff_mass_SMSM.dat");
    distr_t mel_SMSM= Corr.Fit_distr(Corr.mel_ov_mass_t(pt2_distr_SMSM, ""))/2.0;
    distr_t_list mel_SMSM_distr = SQRT_DL( pt2_distr_SMSM/((EXPT_DL( -1.0*eff_mass_SMSM) + EXP_DL(-1.0*eff_mass_SM*Corr.Nt)*EXPT_DL(eff_mass_SMSM) )/(2.0*eff_mass_SMSM)))/(2.0*eff_mass_SMSM);  


    distr_t_list pt2_T_distr= -1*(1.0/3.0)*Corr.corr_t(summ_master(data_2pts_u_T1.col(0)[iens], data_2pts_u_T2.col(0)[iens], data_2pts_u_T3.col(0)[iens]), "../data/ph_emission/"+ph_type_mes+"/"+"C_T/"+data_2pts.Tag[iens]+"/corr_T_2pt.dat");
    distr_t_list pt2_T_distr_SMSM= -1*(1.0/3.0)*Corr.corr_t(summ_master(data_2pts_u_T1_SMSM.col(0)[iens], data_2pts_u_T2_SMSM.col(0)[iens], data_2pts_u_T3_SMSM.col(0)[iens]), "");
    distr_t_list pt2_B_distr= -1*(1.0/3.0)*Corr.corr_t(summ_master(data_2pts_u_B1.col(0)[iens], data_2pts_u_B2.col(0)[iens], data_2pts_u_B3.col(0)[iens]), "../data/ph_emission/"+ph_type_mes+"/"+"C_T/"+data_2pts.Tag[iens]+"/corr_B_2pt.dat");
    distr_t_list eff_mass_T= Corr.effective_mass_t(pt2_T_distr,"../data/ph_emission/"+ph_type_mes+"/"+"mass/"+data_2pts.Tag[iens]+"/eff_mass_T.dat");
    distr_t_list eff_mass_B= Corr.effective_mass_t(pt2_B_distr,"../data/ph_emission/"+ph_type_mes+"/"+"mass/"+data_2pts.Tag[iens]+"/eff_mass_B.dat"); 
    distr_t_list ft_distr = Corr.residue_t( pt2_T_distr, "")/(Corr.Fit_distr(eff_mass_T)*Corr.matrix_element_t(pt2_T_distr_SMSM, ""));


    int Dx;

    if(MESON == "B0s") { Dx= (int)( 0.40*fmTGeV/a_distr.ave() );}
    if(MESON == "B1s") { Dx= (int)( 0.35*fmTGeV/a_distr.ave() );}
    if(MESON == "B2s") { Dx= (int)( 0.30*fmTGeV/a_distr.ave() );}
    if(MESON == "B3s") { Dx= (int)( 0.27*fmTGeV/a_distr.ave() );}
    if(MESON == "B4s") { Dx= (int)( 0.25*fmTGeV/a_distr.ave() );}


    int Dx_2pt;
    if(MESON == "B0s") { Dx_2pt= (int)( 0.50*fmTGeV/a_distr.ave() );}
    if(MESON == "B1s") { Dx_2pt= (int)( 0.42*fmTGeV/a_distr.ave() );}
    if(MESON == "B2s") { Dx_2pt= (int)( 0.34*fmTGeV/a_distr.ave() );}
    if(MESON == "B3s") { Dx_2pt= (int)( 0.28*fmTGeV/a_distr.ave() );}
    if(MESON == "B4s") { Dx_2pt= (int)( 0.25*fmTGeV/a_distr.ave() );}

    cout<<"Dx_2pt: "<<Dx_2pt<<endl;
    Dx_2pt=Dx;
   
    
    distr_t M_P_SM = Corr.Fit_distr(eff_mass_SM);
    distr_t M_P_SMSM = Corr.Fit_distr(eff_mass_SMSM);
    distr_t M_P= 0.5*(M_P_SM+M_P_SMSM);
    auto SINH= [](double x) { return sinh(x);};
    distr_t_list FP_SM_distr_list = (mu + md )*Corr.residue_t( pt2_distr_SM, "")/(M_P_SM*distr_t::f_of_distr(SINH, M_P_SM)*Corr.matrix_element_t(pt2_distr_SMSM, ""));
    distr_t_list FP_SM_II_distr_list = (mu+md)*(pt2_distr_SM/( 0.5*(EXPT_DL( -1.0*eff_mass_SM) + EXP_DL(-1.0*eff_mass_SM*Corr.Nt)*EXPT_DL(eff_mass_SM) )*SINH_DL(eff_mass_SM) ))/SQRT_DL(  pt2_distr_SMSM/( (EXPT_DL(-1.0*eff_mass_SMSM) + EXP_DL(-1.0*eff_mass_SMSM*Corr.Nt)*EXPT_DL(eff_mass_SMSM))/(2.0*eff_mass_SMSM)));
    Print_To_File({}, {FP_SM_distr_list.ave(), FP_SM_distr_list.err(), FP_SM_II_distr_list.ave(), FP_SM_II_distr_list.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"decay_const/"+data_2pts.Tag[iens]+"/decay_const_SM.dat.t", "", "");
    distr_t FP_SM= Corr.Fit_distr( FP_SM_distr_list  );

    //set time interval for FP (second method)
    int Tmin_old =Corr.Tmin; int Tmax_old= Corr.Tmax;

    //set time interval for eff_mass_fit SM
    if(data_2pts.Tag[iens].substr(1,1) =="A") {
      if(MESON == "B1s") { Corr.Tmin=13; Corr.Tmax= 13+11;}
      else if(MESON == "B2s" ) { Corr.Tmin=13; Corr.Tmax= 32;}
      else if(MESON=="B3s") { Corr.Tmin=14; Corr.Tmax=21;}
      else if(MESON=="B4s") { Corr.Tmin=14; Corr.Tmax=21;}
      else { Corr.Tmin=13; Corr.Tmax=13+16;}
      
    }
    else if(data_2pts.Tag[iens] =="cB211b.072.64") {
      if(MESON == "B1s") {Corr.Tmin=22; Corr.Tmax=22+17;}
      else if(MESON == "B2s") {Corr.Tmin=14; Corr.Tmax=27;}
      else if(MESON=="B3s") {Corr.Tmin=13; Corr.Tmax=13+7;}
      else if(MESON=="B4s") { Corr.Tmin=13; Corr.Tmax=21;}
      else {Corr.Tmin=22; Corr.Tmax=22+13;}
    }
    else if(data_2pts.Tag[iens] =="cB211b.072.96") {Corr.Tmin=20; Corr.Tmax=36;}
    
    else if(data_2pts.Tag[iens].substr(1,1) == "C")  {
      if(MESON == "B1s") {Corr.Tmin=21; Corr.Tmax=21+18;}
      else if(MESON == "B2s") {Corr.Tmin=20; Corr.Tmax=20+16;}
      else if(MESON=="B3s") {Corr.Tmin=18; Corr.Tmax=18+10;}
      else if(MESON=="B4s") { Corr.Tmin=18; Corr.Tmax=18+12;}
      else {Corr.Tmin=26; Corr.Tmax=26+18;}
    }
    else if(data_2pts.Tag[iens].substr(1,1) == "D")  {
      if(MESON == "B1s") {Corr.Tmin=26; Corr.Tmax=26+20;}
      else if(MESON == "B2s") { Corr.Tmin=25; Corr.Tmax= 25+17;}
      else if(MESON=="B3s") {Corr.Tmin=25; Corr.Tmax=25+14;}
      else if(MESON=="B4s") { Corr.Tmin=27; Corr.Tmax=27+11;}
      else {Corr.Tmin=31; Corr.Tmax=31+23;}
    }
    else crash("In fixing [Tmin, Tmax] for smeared FP_II, Ensemble: "+data_2pts.Tag[iens]+" not recognized");
    
    //push back fit interval
    Tmin_fp.push_back( Corr.Tmin);
    Tmax_fp.push_back( Corr.Tmax);
    
    distr_t FP_SM_II = Corr.Fit_distr( FP_SM_II_distr_list);
    Corr.Tmin +=Dx_2pt; Corr.Tmax += Dx_2pt;
    
    distr_t FP_SM_II_syst= Corr.Fit_distr( FP_SM_II_distr_list);
    Corr.Tmin=Tmin_old; Corr.Tmax=Tmax_old;
    
    distr_t F_P= FP_SM;
    distr_t F_P_II= FP_SM_II +  fabs( FP_SM_II_syst.ave() - FP_SM_II.ave())*(FP_SM_II - FP_SM_II.ave())/(FP_SM_II.err()); 
    distr_t F_T=Corr.Fit_distr(ft_distr);

    distr_t F_P_fin= F_P_II;
    
    MP_list.distr_list.push_back(M_P/a_distr);
    FP_list.distr_list.push_back(F_P_fin/a_distr);
    FP_diml_list.distr_list.push_back( F_P_II);
    ZT_list.distr_list.push_back(Z_T);
    MP_ov_FP_list.distr_list.push_back( M_P/F_P_fin);
    if(data_2pts.Tag[iens] != "cB211b.072.96") {
      MP_list_red.distr_list.push_back(M_P/a_distr);
    }

    distr_t ZZT = fT_Jpsi_distr*a_distr/F_T;

    //BUG
    distr_t RF= Get_id_distr(Njacks,UseJack); // (a_distr/FP_SM);

    cout<<"aM_P: "<<M_P.ave()<<" +- "<<M_P.err()<<" -> "<< (M_P/(a_distr)).ave()<<" +- "<<(M_P/a_distr).err()<<" GeV"<<endl;
    cout<<"aF_P: "<<F_P.ave()<<" +- "<<F_P.err()<<" -> "<< (F_P/(a_distr)).ave()<<" +- "<<(F_P/a_distr).err()<<" GeV"<<endl;
    cout<<"aF_P_II: "<<F_P_II.ave()<<" +- "<<F_P_II.err()<<" -> "<< (F_P_II/(a_distr)).ave()<<" +- "<<(F_P_II/a_distr).err()<<" GeV"<<endl;
    cout<<"MP/FP: "<<(M_P/F_P).ave()<<" +- "<<(M_P/F_P).err()<<endl;
    cout<<"aM_P(SM): "<<M_P_SM.ave()<<" +- "<<M_P_SM.err()<<" -> "<< (M_P_SM/(a_distr)).ave()<<" +- "<<(M_P_SM/a_distr).err()<<" GeV"<<endl;
    cout<<"aM_P(SMSM): "<<M_P_SMSM.ave()<<" +- "<<M_P_SMSM.err()<<" -> "<<(M_P_SMSM/a_distr).ave()<<" +- "<<(M_P_SMSM/a_distr).err()<<" GeV"<<endl;
    cout<<"MP/FP: "<<(M_P/F_P).ave()<<" +- "<<(M_P/F_P).err()<<endl;
    cout<<"Ensemble: "<<data_2pts.Tag[iens]<<" ZT(ETMC): "<<Z_T.ave()<<" +- "<<Z_T.err()<<" from F_T: "<<ZZT.ave()<<" +- "<<ZZT.err()<<endl;


    //vector ud meson
    distr_t_list pt2_vector_distr= (1.0/3.0)*Corr.corr_t(summ_master(data_2pts_V1.col(0)[iens], data_2pts_V2.col(0)[iens], data_2pts_V3.col(0)[iens]), "../data/ph_emission/"+ph_type_mes+"/"+"C/"+data_2pts.Tag[iens]+"/corr_2pt_V_ud.dat");
    distr_t_list eff_mass_V = Corr.effective_mass_t(pt2_vector_distr, "../data/ph_emission/"+ph_type_mes+"/"+"mass/"+data_2pts.Tag[iens]+"/eff_mass_V_ud.dat");
    //vector u meson
    distr_t_list uu_vector_distr= (1.0/3.0)*Corr.corr_t(summ_master(data_2pts_u_V1.col(0)[iens], data_2pts_u_V2.col(0)[iens], data_2pts_u_V3.col(0)[iens]), "../data/ph_emission/"+ph_type_mes+"/"+"C/"+data_2pts.Tag[iens]+"/corr_2pt_V_uu.dat");
    distr_t_list eff_mass_V_uu= Corr.effective_mass_t(uu_vector_distr, "../data/ph_emission/"+ph_type_mes+"/"+"mass/"+data_2pts.Tag[iens]+"/eff_mass_V_uu.dat");
    //vector d meson
    distr_t_list dd_vector_distr= (1.0/3.0)*Corr.corr_t(summ_master(data_2pts_d_V1.col(0)[iens], data_2pts_d_V2.col(0)[iens], data_2pts_d_V3.col(0)[iens]), "../data/ph_emission/"+ph_type_mes+"/"+"C/"+data_2pts.Tag[iens]+"/corr_2pt_V_dd.dat");
    distr_t_list eff_mass_V_dd= Corr.effective_mass_t(dd_vector_distr, "../data/ph_emission/"+ph_type_mes+"/"+"mass/"+data_2pts.Tag[iens]+"/eff_mass_V_dd.dat");

    //#######################
    
    //#######################


    //define meson mass exponential to be removed
    auto EXP_MES_FUNC = [&] (double a, double b, double c) { return (b<c/2)?1.0/(exp(-a*b)):1.0/(exp(-a*(c-b)));};

    distr_t_list EXP_MES= distr_t_list::f_of_distr( EXP_MES_FUNC, M_P, Corr.Nt);

    //for FV and FA
    vector<vector<vector<distr_t_list>>> Ax_glb, Vec_glb, Ax_u_glb, Ax_d_glb, Vec_u_glb, Vec_d_glb, Ax_Ds_glb, Ax_Ds_u_glb, Ax_Ds_d_glb;
   
    distr_t_list FV(UseJack), FA(UseJack), FV_u(UseJack), FV_d(UseJack), FA_u(UseJack), FA_d(UseJack);
    //for FV_T and FA_T
    vector<vector<vector<distr_t_list>>> Ax_T_glb, Vec_T_glb, Ax_T_u_glb, Ax_T_d_glb, Vec_T_u_glb, Vec_T_d_glb;
    distr_t_list FV_T(UseJack), FA_T(UseJack), FV_T_u(UseJack), FV_T_d(UseJack), FA_T_u(UseJack), FA_T_d(UseJack);
    distr_t_list FT(UseJack), FB(UseJack), FT_u(UseJack), FT_d(UseJack), FB_u(UseJack), FB_d(UseJack);

    distr_t_list xg_list(UseJack);

    //for FV and FA
    Vfloat ax_Tmin, ax_Tmax;
    Vfloat vec_Tmin, vec_Tmax;
    Vfloat ax_u_Tmin, ax_u_Tmax;
    Vfloat vec_u_Tmin, vec_u_Tmax;
    Vfloat ax_d_Tmin, ax_d_Tmax;
    Vfloat vec_d_Tmin, vec_d_Tmax;


    //for FV_T and FA_T
    Vfloat ax_T_Tmin, ax_T_Tmax;
    Vfloat vec_T_Tmin, vec_T_Tmax;
    Vfloat ax_T_u_Tmin, ax_T_u_Tmax;
    Vfloat vec_T_u_Tmin, vec_T_u_Tmax;
    Vfloat ax_T_d_Tmin, ax_T_d_Tmax;
    Vfloat vec_T_d_Tmin, vec_T_d_Tmax;
    

    for(int ixg=0;ixg<num_xg;ixg++) {

      //for FV and FA
      vector<vector<distr_t_list>> Ax_tens(size_mu_nu);
      vector<vector<distr_t_list>> Ax_tens_u(size_mu_nu);
      vector<vector<distr_t_list>> Ax_tens_d(size_mu_nu);
      vector<vector<distr_t_list>> Vec_tens(size_mu_nu);
      vector<vector<distr_t_list>> Vec_tens_u(size_mu_nu);
      vector<vector<distr_t_list>> Vec_tens_d(size_mu_nu);
      //for FV_T and FA_T
      vector<vector<distr_t_list>> Ax_T_tens(size_mu_nu);
      vector<vector<distr_t_list>> Ax_T_tens_u(size_mu_nu);
      vector<vector<distr_t_list>> Ax_T_tens_d(size_mu_nu);
      vector<vector<distr_t_list>> Vec_T_tens(size_mu_nu);
      vector<vector<distr_t_list>> Vec_T_tens_u(size_mu_nu);
      vector<vector<distr_t_list>> Vec_T_tens_d(size_mu_nu);
      //for FA of the Ds meson
      vector<vector<distr_t_list>> Ax_Ds_tens(size_mu_nu);
      vector<vector<distr_t_list>> Ax_Ds_tens_u(size_mu_nu);
      vector<vector<distr_t_list>> Ax_Ds_tens_d(size_mu_nu);


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

	  cout<<"Analyzing FA and FV for  mu: "<<mu+off_i<<" nu: "<<nu+off_i<<endl;
	  int Im_Re;
	  double parity;

	  //vector
	
	  Corr.Reflection_sign = -1;
	  Im_Re=1;
	  parity=1.0;

	  Corr.Perform_Nt_t_average = 0;
	  distr_t_list vec_F_u = Qu*Corr.corr_t(C_V_Fu_data[mu][nu][ixg].col(Im_Re)[iens],"");
	  distr_t_list vec_B_u = Qu*Corr.corr_t(C_V_Bu_data[mu][nu][ixg].col(Im_Re)[iens],"");
	  distr_t_list vec_F_d = Qd*Corr.corr_t(C_V_Fd_data[mu][nu][ixg].col(Im_Re)[iens],"");
	  distr_t_list vec_B_d = Qd*Corr.corr_t(C_V_Bd_data[mu][nu][ixg].col(Im_Re)[iens],"");
	
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
	  distr_t_list ax_F_u = Qu*Corr.corr_t(C_A_Fu_data[mu][nu][ixg].col(Im_Re)[iens],"");
	  distr_t_list ax_B_u = Qu*Corr.corr_t(C_A_Bu_data[mu][nu][ixg].col(Im_Re)[iens],"");
	  distr_t_list ax_F_d = Qd*Corr.corr_t(C_A_Fd_data[mu][nu][ixg].col(Im_Re)[iens],"");
	  distr_t_list ax_B_d = Qd*Corr.corr_t(C_A_Bd_data[mu][nu][ixg].col(Im_Re)[iens],"");

	  //For Ds meson
	  distr_t_list ax_Ds_F_u = qqu*Corr.corr_t(C_A_Fu_data[mu][nu][ixg].col(Im_Re)[iens],"");
	  distr_t_list ax_Ds_B_u = qqu*Corr.corr_t(C_A_Bu_data[mu][nu][ixg].col(Im_Re)[iens],"");
	  distr_t_list ax_Ds_F_d = qqd*Corr.corr_t(C_A_Fd_data[mu][nu][ixg].col(Im_Re)[iens],"");
	  distr_t_list ax_Ds_B_d = qqd*Corr.corr_t(C_A_Bd_data[mu][nu][ixg].col(Im_Re)[iens],"");

	
	  distr_t_list ax_u = parity*(ax_F_u*th_FF  + ax_B_u*th_BB);
	  distr_t_list ax_d = parity*(ax_F_d*th_FF  + ax_B_d*th_BB);
	  distr_t_list ax = ax_u -ax_d;
	  distr_t_list ax_symm=ax;
	  distr_t_list ax_u_symm= ax_u;
	  distr_t_list ax_d_symm= ax_d;


	  //for Ds meson
	  distr_t_list ax_Ds_u = parity*(ax_Ds_F_u*th_FF  + ax_Ds_B_u*th_BB);
	  distr_t_list ax_Ds_d = parity*(ax_Ds_F_d*th_FF  + ax_Ds_B_d*th_BB);
	  distr_t_list ax_Ds = ax_Ds_u -ax_Ds_d;
	  distr_t_list ax_Ds_symm=ax_Ds;
	  distr_t_list ax_Ds_u_symm= ax_Ds_u;
	  distr_t_list ax_Ds_d_symm= ax_Ds_d;
	  
	  //symmetrize ax
	  for(int t=0; t<Corr.Nt;t++) {
	    ax_symm.distr_list[t] = 0.5*(ax.distr_list[t] + Corr.Reflection_sign*ax.distr_list[( Corr.Nt -t)%Corr.Nt]);
	    ax_u_symm.distr_list[t] = 0.5*(ax_u.distr_list[t] + Corr.Reflection_sign*ax_u.distr_list[( Corr.Nt -t)%Corr.Nt]);
	    ax_d_symm.distr_list[t] = 0.5*(ax_d.distr_list[t] + Corr.Reflection_sign*ax_d.distr_list[( Corr.Nt -t)%Corr.Nt]);
	    //for Ds meson
	    ax_Ds_symm.distr_list[t] = 0.5*(ax_Ds.distr_list[t] + Corr.Reflection_sign*ax_Ds.distr_list[( Corr.Nt -t)%Corr.Nt]);
	    ax_Ds_u_symm.distr_list[t] = 0.5*(ax_Ds_u.distr_list[t] + Corr.Reflection_sign*ax_Ds_u.distr_list[( Corr.Nt -t)%Corr.Nt]);
	    ax_Ds_d_symm.distr_list[t] = 0.5*(ax_Ds_d.distr_list[t] + Corr.Reflection_sign*ax_Ds_d.distr_list[( Corr.Nt -t)%Corr.Nt]);
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
	  //Ds meson
	  Ax_Ds_tens[mu].push_back(ax_Ds_symm);
	  Ax_Ds_tens_u[mu].push_back(ax_Ds_u_symm);
	  Ax_Ds_tens_d[mu].push_back(ax_Ds_d_symm);
	

	}


	for(int nu=0; nu<size_mu_nu_T; nu++) {


	  cout<<"Analyzing FA_T and FV_T for  mu: "<<mu+off_i<<" nu: "<<nu+off_T<<endl;
	  int Im_Re;
	  double parity;

	  //tensor- electric part
	  Corr.Reflection_sign= -1;
	  if( (mu+off_i) != (nu+off_T)) Corr.Reflection_sign = 1;
	  Im_Re=1;
	  parity=1.0;

	  Corr.Perform_Nt_t_average = 0;
	  distr_t_list vec_T_F_u = Qu*Corr.corr_t(C_TK_Fu_data[mu][nu][ixg].col(Im_Re)[iens],"");
	  distr_t_list vec_T_B_u = Qu*Corr.corr_t(C_TK_Bu_data[mu][nu][ixg].col(Im_Re)[iens],"");
	  distr_t_list vec_T_F_d = Qd*Corr.corr_t(C_TK_Fd_data[mu][nu][ixg].col(Im_Re)[iens],"");
	  distr_t_list vec_T_B_d = Qd*Corr.corr_t(C_TK_Bd_data[mu][nu][ixg].col(Im_Re)[iens],"");
	
	  distr_t_list vec_T_u = th_FF*vec_T_F_u + th_BB*vec_T_B_u;
	  distr_t_list vec_T_d = th_FF*vec_T_F_d + th_BB*vec_T_B_d;
	  distr_t_list vec_T = vec_T_u + vec_T_d;
	  distr_t_list vec_T_symm= vec_T;
	  distr_t_list vec_T_u_symm= vec_T_u;
	  distr_t_list vec_T_d_symm= vec_T_d;

	  //symmetrize vec
	  for(int t=0; t<Corr.Nt;t++) {
	    vec_T_symm.distr_list[t] = 0.5*(vec_T.distr_list[t] + Corr.Reflection_sign*vec_T.distr_list[( Corr.Nt -t)%Corr.Nt]);
	    vec_T_u_symm.distr_list[t] = 0.5*(vec_T_u.distr_list[t] + Corr.Reflection_sign*vec_T_u.distr_list[( Corr.Nt -t)%Corr.Nt]);
	    vec_T_d_symm.distr_list[t] = 0.5*(vec_T_d.distr_list[t] + Corr.Reflection_sign*vec_T_d.distr_list[( Corr.Nt -t)%Corr.Nt]);
	  }
	  Corr.Perform_Nt_t_average=1;


	  //tensor - magnetic part
	  if(  (  ((mu+off_i) ==0) || ((nu+off_T)==0)) && ( ((mu+off_i) != 0) || ( (nu+off_T) != 0))) {Im_Re=1; Corr.Reflection_sign=-1; parity=1.0;}
	  else {
	    Im_Re=0;
	    Corr.Reflection_sign=1;
	    if(  (mu+off_i) == (nu + off_T)) Corr.Reflection_sign=-1; 
	    parity=1.0;
	  }


	  Corr.Perform_Nt_t_average=0;
	  distr_t_list ax_T_F_u = Qu*Corr.corr_t(C_BK_Fu_data[mu][nu][ixg].col(Im_Re)[iens],"");
	  distr_t_list ax_T_B_u = Qu*Corr.corr_t(C_BK_Bu_data[mu][nu][ixg].col(Im_Re)[iens],"");
	  distr_t_list ax_T_F_d = Qd*Corr.corr_t(C_BK_Fd_data[mu][nu][ixg].col(Im_Re)[iens],"");
	  distr_t_list ax_T_B_d = Qd*Corr.corr_t(C_BK_Bd_data[mu][nu][ixg].col(Im_Re)[iens],"");

	
	  distr_t_list ax_T_u = parity*(ax_T_F_u*th_FF  + ax_T_B_u*th_BB);
	  distr_t_list ax_T_d = parity*(ax_T_F_d*th_FF  + ax_T_B_d*th_BB);
	  distr_t_list ax_T = ax_T_u +ax_T_d;
	  distr_t_list ax_T_symm=ax_T;

	  distr_t_list ax_T_u_symm= ax_T_u;
	  distr_t_list ax_T_d_symm= ax_T_d;
	  //symmetrize ax
	  for(int t=0; t<Corr.Nt;t++) {
	    ax_T_symm.distr_list[t] = 0.5*(ax_T.distr_list[t] + Corr.Reflection_sign*ax_T.distr_list[( Corr.Nt -t)%Corr.Nt]);
	    ax_T_u_symm.distr_list[t] = 0.5*(ax_T_u.distr_list[t] + Corr.Reflection_sign*ax_T_u.distr_list[( Corr.Nt -t)%Corr.Nt]);
	    ax_T_d_symm.distr_list[t] = 0.5*(ax_T_d.distr_list[t] + Corr.Reflection_sign*ax_T_d.distr_list[( Corr.Nt -t)%Corr.Nt]);
	  }
	  Corr.Perform_Nt_t_average=1;

	  //restore standard reflection sign
	  Corr.Reflection_sign=1;
	  
	  //get H tensor
	  distr_t_list H_T_A_u= ax_T_u*EXP_MES*EXP_PH;
	  distr_t_list H_T_A_d= ax_T_d*EXP_MES*EXP_PH;
	  distr_t_list H_T_V_u= vec_T_u*EXP_MES*EXP_PH;
	  distr_t_list H_T_V_d= vec_T_d*EXP_MES*EXP_PH;
	  distr_t_list H_T_A_u_symm= ax_T_u_symm*EXP_MES*EXP_PH;
	  distr_t_list H_T_A_d_symm= ax_T_d_symm*EXP_MES*EXP_PH;
	  distr_t_list H_T_V_u_symm= vec_T_u_symm*EXP_MES*EXP_PH;
	  distr_t_list H_T_V_d_symm= vec_T_d_symm*EXP_MES*EXP_PH;
	  distr_t_list H_T_A_symm= ax_T_symm*EXP_MES*EXP_PH;
	  distr_t_list H_T_V_symm= vec_T_symm*EXP_MES*EXP_PH;



	  //print forward and backward contribution to C
	  Print_To_File({}, {ax_T_F_u.ave(), ax_T_F_u.err(), ax_T_F_d.ave(), ax_T_F_d.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"C_T/"+data_2pts.Tag[iens]+"/quark_contr_A_wo_T_symm_FF_mu_"+to_string(mu+off_i)+"_nu_"+to_string(nu+off_T)+"_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "#   ax_F_u    ax_F_d");
	  Print_To_File({}, {ax_T_B_u.ave(), ax_T_B_u.err(), ax_T_B_d.ave(), ax_T_B_d.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"C_T/"+data_2pts.Tag[iens]+"/quark_contr_A_wo_T_symm_BB_mu_"+to_string(mu+off_i)+"_nu_"+to_string(nu+off_T)+"_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "#   ax_B_u    ax_B_d");
	  Print_To_File({}, {vec_T_F_u.ave(), vec_T_F_u.err(), vec_T_F_d.ave(), vec_T_F_d.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"C_T/"+data_2pts.Tag[iens]+"/quark_contr_V_wo_T_symm_FF_mu_"+to_string(mu+off_i)+"_nu_"+to_string(nu+off_T)+"_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "#   vec_F_u    vec_F_d");
	  Print_To_File({}, {vec_T_B_u.ave(), vec_T_B_u.err(), vec_T_B_d.ave(), vec_T_B_d.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"C_T/"+data_2pts.Tag[iens]+"/quark_contr_V_wo_T_symm_BB_mu_"+to_string(mu+off_i)+"_nu_"+to_string(nu+off_T)+"_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "#   vec_B_u    vec_B_d");
	

	  //single contributions to C
	  Print_To_File({}, {ax_T_u.ave() ,ax_T_u.err(), ax_T_d.ave(), ax_T_d.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"C_T/"+data_2pts.Tag[iens]+"/quark_contr_A_wo_T_symm_mu_"+to_string(mu+off_i)+"_nu_"+to_string(nu+off_T)+"_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "#   ax_u    ax_d");
	  Print_To_File({}, {vec_T_u.ave() ,vec_T_u.err(), vec_T_d.ave(), vec_T_d.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"C_T/"+data_2pts.Tag[iens]+"/quark_contr_V_wo_T_symm_mu_"+to_string(mu+off_i)+"_nu_"+to_string(nu+off_T)+"_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "#   vec_u    vec_d");

	  //single contributions to C with symmetrization
	  Print_To_File({}, {ax_T_u_symm.ave() ,ax_T_u_symm.err(), ax_T_d_symm.ave(), ax_T_d_symm.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"C_T/"+data_2pts.Tag[iens]+"/quark_contr_A_mu_"+to_string(mu+off_i)+"_nu_"+to_string(nu+off_T)+"_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "#   ax_u    ax_d");
	  Print_To_File({}, {vec_T_u_symm.ave() ,vec_T_u_symm.err(), vec_T_d_symm.ave(), vec_T_d_symm.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"C_T/"+data_2pts.Tag[iens]+"/quark_contr_V_mu_"+to_string(mu+off_i)+"_nu_"+to_string(nu+off_T)+"_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "#   vec_u    vec_d");


	  //total contribution without symmetrization to C

	  Print_To_File({}, {ax_T.ave(), ax_T.err(), vec_T.ave(), vec_T.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"C_T/"+data_2pts.Tag[iens]+"/wo_T_symm_mu_"+to_string(mu+off_i)+"_nu_"+to_string(nu+off_T)+"_xg_"+to_string_with_precision(xg.ave(),2)+".dat" , "", "#  t    A     V");

	  //total contribution to C
	  Print_To_File({}, {ax_T_symm.ave(), ax_T_symm.err(), vec_T_symm.ave(), vec_T_symm.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"C_T/"+data_2pts.Tag[iens]+"/mu_"+to_string(mu+off_i)+"_nu_"+to_string(nu+off_T)+"_xg_"+to_string_with_precision(xg.ave(),2)+".dat" , "", "#  t    A     V");


	  //###############################################################


	  //########################   H TENSOR ###########################


	  //single contributions to H without symmetrization
	   Print_To_File({}, {H_T_A_u.ave() ,H_T_A_u.err(), H_T_A_d.ave(), H_T_A_d.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"H_T/"+data_2pts.Tag[iens]+"/quark_contr_A_wo_symm_mu_"+to_string(mu+off_i)+"_nu_"+to_string(nu+off_T)+"_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "#   HA_u    HA_d");
	  Print_To_File({}, {H_T_V_u.ave() ,H_T_V_u.err(), H_T_V_d.ave(), H_T_V_d.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"H_T/"+data_2pts.Tag[iens]+"/quark_contr_V_wo_symm_mu_"+to_string(mu+off_i)+"_nu_"+to_string(nu+off_T)+"_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "#   HV_u    HV_d");
	  

	  
	  //single contributions to H with symmetrization
	  Print_To_File({}, {H_T_A_u_symm.ave() ,H_T_A_u_symm.err(), H_T_A_d_symm.ave(), H_T_A_d_symm.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"H_T/"+data_2pts.Tag[iens]+"/quark_contr_A_mu_"+to_string(mu+off_i)+"_nu_"+to_string(nu+off_T)+"_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "#   HA_u    HA_d");
	  Print_To_File({}, {H_T_V_u_symm.ave() ,H_T_V_u_symm.err(), H_T_V_d_symm.ave(), H_T_V_d_symm.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"H_T/"+data_2pts.Tag[iens]+"/quark_contr_V_mu_"+to_string(mu+off_i)+"_nu_"+to_string(nu+off_T)+"_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "#   HV_u    HV_d");

	  //total contribution to H
	  Print_To_File({}, {H_T_A_symm.ave(), H_T_A_symm.err(), H_T_V_symm.ave(), H_T_V_symm.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"H_T/"+data_2pts.Tag[iens]+"/mu_"+to_string(mu+off_i)+"_nu_"+to_string(nu+off_T)+"_xg_"+to_string_with_precision(xg.ave(),2)+".dat" , "", "#  t    HA     HV");


	  
	  //push back
	  Ax_T_tens[mu].push_back(ax_T_symm);
	  Ax_T_tens_u[mu].push_back(ax_T_u_symm);
	  Ax_T_tens_d[mu].push_back(ax_T_d_symm);
	  Vec_T_tens[mu].push_back(vec_T_symm);
	  Vec_T_tens_u[mu].push_back(vec_T_u_symm);
	  Vec_T_tens_d[mu].push_back(vec_T_d_symm);
	  

	}
      }


      //######################################################################################################################
      //####################     STANDARD     FV     AND    FA     COMPUTATION     ###########################

      //######

      //push_back Ax_tens and Vec_tens
      Ax_glb.push_back(Ax_tens);
      Ax_u_glb.push_back(Ax_tens_u);
      Ax_d_glb.push_back(Ax_tens_d);
      
      Vec_glb.push_back(Vec_tens);
      Vec_u_glb.push_back(Vec_tens_u);
      Vec_d_glb.push_back(Vec_tens_d);

     
      //For Ds
      Ax_Ds_glb.push_back(Ax_Ds_tens);
      Ax_Ds_u_glb.push_back(Ax_Ds_tens_u);
      Ax_Ds_d_glb.push_back(Ax_Ds_tens_d);

      distr_t_list FA0_distr= 0.5*(Ax_Ds_glb[0][1-off_i][1-off_i] + Ax_Ds_glb[0][2-off_i][2-off_i]);


    
          

      if(ixg==0) {

	distr_t_list F_3pt_rest_u = (0.5/qqu)*(Ax_Ds_u_glb[0][1-off_i][1-off_i] + Ax_Ds_u_glb[0][2-off_i][2-off_i]);
	distr_t_list F_3pt_rest_d = (0.5/qqd)*(Ax_Ds_d_glb[0][1-off_i][1-off_i] + Ax_Ds_d_glb[0][2-off_i][2-off_i]);
	distr_t_list F_3pt_rest= F_3pt_rest_u*qqu -F_3pt_rest_d*qqd;
	distr_t_list F_3pt_rest_opt(UseJack);
	for(int t=0; t < Corr.Nt;t++) {
	  double a1= 1.0/pow(F_3pt_rest_u.err(t),2);
	  double a2= 1.0/pow(F_3pt_rest_d.err(t),2);
	  double r= a1+a2;
	  a1 /= r;
	  a2 /= r;
	  F_3pt_rest_opt.distr_list.push_back( F_3pt_rest_u.distr_list[t]*a1 + F_3pt_rest_d.distr_list[t]*a2);
	}
	
	Print_To_File({}, {F_3pt_rest.ave(), F_3pt_rest.err(), F_3pt_rest_opt.ave(), F_3pt_rest_opt.err(), F_3pt_rest_u.ave(), F_3pt_rest_u.err(), F_3pt_rest_d.ave(), F_3pt_rest_d.err()},  "../data/ph_emission/"+ph_type_mes+"/"+"decay_const/"+data_2pts.Tag[iens]+"/FA_3pt_rest.dat.t", "", "");
	
	distr_t_list FP_3pt_distr_u = -1.0*Zv*F_3pt_rest_u*EXPT_DL(eff_mass_SMSM)/(mel_SMSM_distr);
	distr_t_list FP_3pt_distr_std= -1.0*Zv*F_3pt_rest*EXPT_DL(eff_mass_SMSM)/(mel_SMSM_distr);
	distr_t_list FP_3pt_distr_d= -1.0*Zv*F_3pt_rest_d*EXPT_DL(eff_mass_SMSM)/(mel_SMSM_distr);
	distr_t_list FP_3pt_distr_opt= -1.0*Zv*F_3pt_rest_opt*EXPT_DL(eff_mass_SMSM)/(mel_SMSM_distr);
	
	distr_t_list FP_3pt_distr= FP_3pt_distr_u;
	
	int Tmin_old =Corr.Tmin; int Tmax_old= Corr.Tmax;
	
	//set time interval for eff_mass_fit SM
	if(data_2pts.Tag[iens].substr(1,1) =="A") {
	  if(MESON == "B1s") { Corr.Tmin=10; Corr.Tmax= 21;}
	  else if(MESON == "B2s" ) { Corr.Tmin=10; Corr.Tmax= 20;}
	  else if(MESON=="B3s") { Corr.Tmin=10; Corr.Tmax=19;}
	  else if(MESON=="B4s") { Corr.Tmin=10; Corr.Tmax=19;}
	  else { Corr.Tmin=10; Corr.Tmax=26;}
	  
	}
	else if(data_2pts.Tag[iens] =="cB211b.072.64") {
	  if(MESON == "B1s") {Corr.Tmin=23; Corr.Tmax=33;}
	  else if(MESON == "B2s") {Corr.Tmin=18; Corr.Tmax=23;}
	  else if(MESON=="B3s") {Corr.Tmin=18; Corr.Tmax=23;}
	  else if(MESON=="B4s") { Corr.Tmin=17; Corr.Tmax=26;}
	  else {Corr.Tmin=18; Corr.Tmax=30;}
	}
	else if(data_2pts.Tag[iens] =="cB211b.072.96") {Corr.Tmin=20; Corr.Tmax=36;}
	
	else if(data_2pts.Tag[iens].substr(1,1) == "C")  {
	  if(MESON == "B1s") {Corr.Tmin=20; Corr.Tmax=31;}
	  else if(MESON == "B2s") {Corr.Tmin=17; Corr.Tmax=36;}
	  else if(MESON=="B3s") {Corr.Tmin=17; Corr.Tmax=31;}
	  else if(MESON=="B4s") { Corr.Tmin=17; Corr.Tmax=31;}
	  else {Corr.Tmin=20; Corr.Tmax=38;}
	}
	else if(data_2pts.Tag[iens].substr(1,1) == "D")  {
	  if(MESON == "B1s") {Corr.Tmin=21; Corr.Tmax=39;}
	  else if(MESON == "B2s") { Corr.Tmin=25; Corr.Tmax= 45;}
	  else if(MESON=="B3s") {Corr.Tmin=25; Corr.Tmax=44;}
	  else if(MESON=="B4s") { Corr.Tmin=25; Corr.Tmax=39;}
	  else {Corr.Tmin=30; Corr.Tmax=47;}
	}
	else crash("In fixing [Tmin, Tmax] for smeared FP_III, Ensemble: "+data_2pts.Tag[iens]+" not recognized");

	bool mode_c=true;
	if(mode_c) {
	//set time interval for eff_mass_fit SM
	if(data_2pts.Tag[iens].substr(1,1) =="A") {
	  if(MESON == "B1s") { Corr.Tmin=16; Corr.Tmax= 24;}
	  else if(MESON == "B2s" ) { Corr.Tmin=15; Corr.Tmax= 24;}
	  else if(MESON=="B3s") { Corr.Tmin=13; Corr.Tmax=19;}
	  else if(MESON=="B4s") { Corr.Tmin=13; Corr.Tmax=19;}
	  else { Corr.Tmin=16; Corr.Tmax=26;}
	  
	}
	else if(data_2pts.Tag[iens] =="cB211b.072.64") {
	  if(MESON == "B1s") {Corr.Tmin=18; Corr.Tmax=27;}
	  else if(MESON == "B2s") {Corr.Tmin=16; Corr.Tmax=28;}
	  else if(MESON=="B3s") {Corr.Tmin=16; Corr.Tmax=16+9;}
	  else if(MESON=="B4s") { Corr.Tmin=15; Corr.Tmax=15+7;}
	  else {Corr.Tmin=22; Corr.Tmax=31;}
	}
	else if(data_2pts.Tag[iens] =="cB211b.072.96") {Corr.Tmin=20; Corr.Tmax=36;}
	
	else if(data_2pts.Tag[iens].substr(1,1) == "C")  {
	  if(MESON == "B1s") {Corr.Tmin=20; Corr.Tmax=31;}
	  else if(MESON == "B2s") {Corr.Tmin=18; Corr.Tmax=28;}
	  else if(MESON=="B3s") {Corr.Tmin=19; Corr.Tmax=19+10;}
	  else if(MESON=="B4s") { Corr.Tmin=18; Corr.Tmax=27;}
	  else {Corr.Tmin=14; Corr.Tmax=14+17;}
	}
	else if(data_2pts.Tag[iens].substr(1,1) == "D")  {
	  if(MESON == "B1s") {Corr.Tmin=21; Corr.Tmax=21+16;}
	  else if(MESON == "B2s") { Corr.Tmin=16; Corr.Tmax= 16+12;}
	  else if(MESON=="B3s") {Corr.Tmin=16; Corr.Tmax=16+12;}
	  else if(MESON=="B4s") { Corr.Tmin=16; Corr.Tmax=27;}
	  else {Corr.Tmin=25; Corr.Tmax=25+17;}
	}
	else crash("In fixing [Tmin, Tmax] for smeared FP_III mode_c, Ensemble: "+data_2pts.Tag[iens]+" not recognized");
	}
	

	
	
	//push back fit interval
	Tmin_fp_3pt.push_back( Corr.Tmin);
	Tmax_fp_3pt.push_back( Corr.Tmax);
	
	distr_t FP_3pt = Corr.Fit_distr( FP_3pt_distr);
	Corr.Tmin += Dx_2pt; Corr.Tmax += Dx_2pt;
	distr_t FP_3pt_syst= Corr.Fit_distr(FP_3pt_distr);
	FP_3pt= FP_3pt + fabs( FP_3pt.ave() - FP_3pt_syst.ave())*( FP_3pt - FP_3pt.ave())/FP_3pt.err();
	
	Corr.Tmin=Tmin_old; Corr.Tmax=Tmax_old;
	
	
	FP_3pt_list.distr_list.push_back( FP_3pt/a_distr);
	FP_3pt_diml_list.distr_list.push_back( FP_3pt);
	MP_ov_FP_3pt_list.distr_list.push_back( M_P/FP_3pt);
	
	Print_To_File({}, {FP_3pt_distr_u.ave(), FP_3pt_distr_u.err(),  FP_3pt_distr_d.ave(), FP_3pt_distr_d.err(), FP_3pt_distr_std.ave(), FP_3pt_distr_std.err(), FP_3pt_distr_opt.ave(), FP_3pt_distr_opt.err(),},  "../data/ph_emission/"+ph_type_mes+"/"+"decay_const/"+data_2pts.Tag[iens]+"/decay_const_3pt.dat.t", "", "# u  d   std  opt");
      }  
      
      //define exponential mass factor for each form factor
      //RIGHT ONE IS sign_kz = -1; //relative sign in tensor FF contributions
      double sign_kz= -1;
      vector<distr_t_list> EXP_MASS_CH;
      auto F_EXP = [](double m, double t) { return exp(m*t) ;};
      vector<string> ff_tags({"Au", "Vu", "Ad", "Vd", "ATu", "VTu", "ATd", "VTd"});
      //######
         
      EXP_MASS_CH.push_back( distr_t_list::f_of_distr_list(F_EXP, Corr.effective_mass_t( (Ax_tens_u[1-off_i][1-off_i] + Ax_tens_u[1-off_i][1-off_i])*EXP_PH  -Qu*FA0_distr, "")));
      EXP_MASS_CH.push_back( distr_t_list::f_of_distr_list(F_EXP, Corr.effective_mass_t( (Vec_tens_u[1-off_i][2-off_i] - Vec_tens_u[2-off_i][1-off_i])*EXP_PH, "")));
      
      EXP_MASS_CH.push_back( distr_t_list::f_of_distr_list(F_EXP, Corr.effective_mass_t( (Ax_tens_d[1-off_i][1-off_i] + Ax_tens_d[1-off_i][1-off_i])*EXP_PH  -Qd*FA0_distr, "")));
      EXP_MASS_CH.push_back( distr_t_list::f_of_distr_list(F_EXP, Corr.effective_mass_t( (Vec_tens_d[1-off_i][2-off_i] - Vec_tens_d[2-off_i][1-off_i])*EXP_PH, "")));

      EXP_MASS_CH.push_back( distr_t_list::f_of_distr_list(F_EXP, Corr.effective_mass_t( (Ax_T_tens_u[1-off_i][1-off_T]*(1-xg/2) - sign_kz*Vec_T_tens_u[1-off_i][2-off_T]*(xg/2) + Ax_T_tens_u[2-off_i][2-off_T]*(1-xg/2) +sign_kz*Vec_T_tens_u[2-off_i][1-off_T]*(xg/2))*EXP_PH, "")));
      EXP_MASS_CH.push_back( distr_t_list::f_of_distr_list(F_EXP, Corr.effective_mass_t(  ((Vec_T_tens_u[1-off_i][2-off_T]*(1-xg/2) -sign_kz*Ax_T_tens_u[1-off_i][1-off_T]*(xg/2)) -(Vec_T_tens_u[2-off_i][1-off_i]*(1-xg/2) + sign_kz*Ax_T_tens_u[2-off_i][2-off_T]*(xg/2)))*EXP_PH, "") ));

      EXP_MASS_CH.push_back( distr_t_list::f_of_distr_list(F_EXP, Corr.effective_mass_t( (Ax_T_tens_d[1-off_i][1-off_T]*(1-xg/2) - sign_kz*Vec_T_tens_d[1-off_i][2-off_T]*(xg/2) + Ax_T_tens_d[2-off_i][2-off_T]*(1-xg/2) +sign_kz*Vec_T_tens_d[2-off_i][1-off_T]*(xg/2))*EXP_PH, "")));
      EXP_MASS_CH.push_back( distr_t_list::f_of_distr_list(F_EXP, Corr.effective_mass_t(  ((Vec_T_tens_d[1-off_i][2-off_T]*(1-xg/2) -sign_kz*Ax_T_tens_d[1-off_i][1-off_T]*(xg/2)) -(Vec_T_tens_d[2-off_i][1-off_i]*(1-xg/2) + sign_kz*Ax_T_tens_d[2-off_i][2-off_T]*(xg/2)))*EXP_PH, "") ));

          
    
       
      

      //Compute FV and FA   
      distr_t_list FA_distr = -1.0*RF*(0.5*(Ax_tens[1-off_i][1-off_i] + Ax_tens[2-off_i][2-off_i])*EXP_PH - 0.0*FA0_distr)*(1.0/Eg)*(FP_SM/(-1.0*FA0_distr));
      distr_t_list FA_new_distr= -1.0*RF*Zv*(0.5*(Ax_tens[1-off_i][1-off_i] + Ax_tens[2-off_i][2-off_i])*EXP_PH)*(1.0/(Eg*mel_SMSM))*EXP_MES; 
      distr_t_list FV_distr = 0.5*RF*(Za/Zv)*(FP_SM/(-1.0*FA0_distr))*( Vec_tens[1-off_i][2-off_i] - Vec_tens[2-off_i][1-off_i])*EXP_PH/kz;
      distr_t_list FV_new_distr = (0.5*RF*Za*( Vec_tens[1-off_i][2-off_i] - Vec_tens[2-off_i][1-off_i])*EXP_PH/(kz*mel_SMSM))*EXP_MES;
      distr_t_list FV_sub_distr= 0.5*RF*(Za/Zv)*(FP_SM/(-1.0*FA0_distr))*( (Vec_tens[1-off_i][2-off_i] - Vec_tens[2-off_i][1-off_i])*EXP_PH - (Vec_glb[0][1-off_i][2-off_i] + Vec_glb[0][2-off_i][1-off_i])*1.0)/kz;
      //compute FV and FA (up-component)
      distr_t_list FA0_u_distr= 0.5*RF*(Ax_u_glb[0][1-off_i][1-off_i] + Ax_u_glb[0][2-off_i][2-off_i]);
      distr_t_list FA_u_distr = -1.0*RF*(0.5*(Ax_tens_u[1-off_i][1-off_i] + Ax_tens_u[2-off_i][2-off_i])*EXP_PH - Qu*FA0_distr)*(1.0/Eg)*(FP_SM/(-1.0*FA0_distr));
      distr_t_list FA_u_new_distr= -1.0*RF*Zv*(0.5*(Ax_tens_u[1-off_i][1-off_i] + Ax_tens_u[2-off_i][2-off_i])*EXP_PH - Qu*FA0_distr)*(1.0/(Eg*mel_SMSM))*EXP_MES;
      distr_t_list FA_u_var= -1.0*RF*Zv*(0.5*(Ax_tens_u[1-off_i][1-off_i] + Ax_tens_u[2-off_i][2-off_i])*EXP_PH - Qu*FA0_distr)*(1.0/(Eg*mel_SMSM))*EXP_MASS_CH[0];
      distr_t_list FV_u_distr = 0.5*RF*(Za/Zv)*(FP_SM/(-1.0*FA0_distr))*( Vec_tens_u[1-off_i][2-off_i] - Vec_tens_u[2-off_i][1-off_i])*EXP_PH/kz;
      distr_t_list FV_u_new_distr = (0.5*RF*Za*( Vec_tens_u[1-off_i][2-off_i] - Vec_tens_u[2-off_i][1-off_i])*EXP_PH/(kz*mel_SMSM))*EXP_MES;
      distr_t_list FV_u_var = (0.5*RF*Za*( Vec_tens_u[1-off_i][2-off_i] - Vec_tens_u[2-off_i][1-off_i])*EXP_PH/(kz*mel_SMSM))*EXP_MASS_CH[1];
      distr_t_list FV_u_sub_distr= 0.5*RF*(Za/Zv)*(FP_SM/(-1.0*FA0_distr))*( (Vec_tens_u[1-off_i][2-off_i] - Vec_tens_u[2-off_i][1-off_i])*EXP_PH - (Vec_u_glb[0][1-off_i][2-off_i] + Vec_u_glb[0][2-off_i][1-off_i])*1.0)/kz;
      //compute FV and FA (d-component)
      distr_t_list FA0_d_distr= 0.5*RF*(Ax_d_glb[0][1-off_i][1-off_i] + Ax_d_glb[0][2-off_i][2-off_i]);
      distr_t_list FA_d_distr = -1.0*RF*(0.5*(Ax_tens_d[1-off_i][1-off_i] + Ax_tens_d[2-off_i][2-off_i])*EXP_PH - Qd*FA0_distr)*(1.0/Eg)*(FP_SM/(-1.0*FA0_distr));
      distr_t_list FA_d_new_distr= -1.0*RF*Zv*(0.5*(Ax_tens_d[1-off_i][1-off_i] + Ax_tens_d[2-off_i][2-off_i])*EXP_PH - Qd*FA0_distr)*(1.0/(Eg*mel_SMSM))*EXP_MES;
      distr_t_list FA_d_var= -1.0*RF*Zv*(0.5*(Ax_tens_d[1-off_i][1-off_i] + Ax_tens_d[2-off_i][2-off_i])*EXP_PH - Qd*FA0_distr)*(1.0/(Eg*mel_SMSM))*EXP_MASS_CH[2];
      distr_t_list FV_d_distr = 0.5*RF*(Za/Zv)*(FP_SM/(-1.0*FA0_distr))*( Vec_tens_d[1-off_i][2-off_i] - Vec_tens_d[2-off_i][1-off_i])*EXP_PH/kz;
      distr_t_list FV_d_new_distr=   (0.5*RF*Za*( Vec_tens_d[1-off_i][2-off_i] - Vec_tens_d[2-off_i][1-off_i])*EXP_PH/(kz*mel_SMSM))*EXP_MES;
      distr_t_list FV_d_var=   (0.5*RF*Za*( Vec_tens_d[1-off_i][2-off_i] - Vec_tens_d[2-off_i][1-off_i])*EXP_PH/(kz*mel_SMSM))*EXP_MASS_CH[3];
      distr_t_list FV_d_sub_distr= 0.5*RF*(Za/Zv)*(FP_SM/(-1.0*FA0_distr))*( (Vec_tens_d[1-off_i][2-off_i] - Vec_tens_d[2-off_i][1-off_i])*EXP_PH - (Vec_d_glb[0][1-off_i][2-off_i] + Vec_d_glb[0][2-off_i][1-off_i])*1.0)/kz;

      
      
      //Print FV and FA
      Print_To_File({}, {FA_distr.ave(), FA_distr.err(), FA0_distr.ave(), FA0_distr.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"FF/"+data_2pts.Tag[iens]+"/FA_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "");
      Print_To_File({}, {FV_distr.ave(), FV_distr.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"FF/"+data_2pts.Tag[iens]+"/FV_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "");
      Print_To_File({}, {FV_sub_distr.ave(), FV_sub_distr.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"FF/"+data_2pts.Tag[iens]+"/FV_sub_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "");

      //Print FV and FA (up-component)
      Print_To_File({}, {FA_u_distr.ave(), FA_u_distr.err(), FA_u_new_distr.ave(), FA_u_new_distr.err(),  FA_u_var.ave(), FA_u_var.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"FF_u/"+data_2pts.Tag[iens]+"/FA_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "");
      Print_To_File({}, {FV_u_distr.ave(), FV_u_distr.err(), FV_u_new_distr.ave(), FV_u_new_distr.err(),  FV_u_var.ave(), FV_u_var.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"FF_u/"+data_2pts.Tag[iens]+"/FV_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "");
      Print_To_File({}, {FV_u_sub_distr.ave(), FV_u_sub_distr.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"FF_u/"+data_2pts.Tag[iens]+"/FV_sub_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "");
      
      //Print FV and FA (d-component)
      Print_To_File({}, {FA_d_distr.ave(), FA_d_distr.err(), FA_d_new_distr.ave(), FA_d_new_distr.err(),  FA_d_var.ave(), FA_d_var.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"FF_d/"+data_2pts.Tag[iens]+"/FA_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "");
      Print_To_File({}, {(FV_d_distr).ave(), (FV_d_distr).err(), FV_d_new_distr.ave(), FV_d_new_distr.err(),  FV_d_var.ave(), FV_d_var.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"FF_d/"+data_2pts.Tag[iens]+"/FV_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "");
      Print_To_File({}, {(FV_d_sub_distr).ave(), (FV_d_sub_distr).err()}, "../data/ph_emission/"+ph_type_mes+"/"+"FF_d/"+data_2pts.Tag[iens]+"/FV_sub_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "");

      //fit the form factors


      
      //set interval depending on xg
      int Tmin_V, Tmax_V, Tmin_A, Tmax_A;
      if(xg.ave() > 1e-10) {
      //set time intervals for each heavy mass
      //FA_u
      Get_Bs_Mh_Tmin_Tmax("FAu",MESON, Tmin_A, Tmax_A, ixg, data_2pts.Tag[iens]);
      Corr.Tmin = Tmin_A; Corr.Tmax= Tmax_A;
      distr_t FA_u_fit= Corr.Fit_distr(FA_u_new_distr);
      Corr.Tmin += Dx; Corr.Tmax += Dx;
      distr_t FA_u_syst_fit= Corr.Fit_distr(FA_u_new_distr);
      Corr.Tmin -= Dx; Corr.Tmax -= Dx;
      FA_u_fit = FA_u_fit + fabs(FA_u_fit.ave()-FA_u_syst_fit.ave())*erf(k_erf*fabs((FA_u_fit.ave()-FA_u_syst_fit.ave()))/(sqrt(2.0)*FA_u_syst_fit.err()))*(FA_u_fit - FA_u_fit.ave())/(FA_u_fit.err());
      FA_u.distr_list.push_back(FA_u_fit);
      FA_u_per_kin[ixg-1].distr_list.push_back(FA_u_fit);
      ax_u_Tmin.push_back( Corr.Tmin); ax_u_Tmax.push_back( Corr.Tmax);
      //FA_d
      Get_Bs_Mh_Tmin_Tmax("FAd",MESON, Tmin_A, Tmax_A, ixg, data_2pts.Tag[iens]);
      Corr.Tmin = Tmin_A; Corr.Tmax= Tmax_A;
      distr_t FA_d_fit= Corr.Fit_distr(FA_d_new_distr);
      if(MESON=="B3s" || MESON=="B4s") Dx+= (int)( 0.1*fmTGeV/a_distr.ave());
      Corr.Tmin += Dx; Corr.Tmax += Dx;
      distr_t FA_d_syst_fit= Corr.Fit_distr(FA_d_new_distr);
      Corr.Tmin -= Dx; Corr.Tmax -= Dx;
      if(MESON=="B3s" || MESON=="B4s") Dx-= (int)( 0.1*fmTGeV/a_distr.ave());
      FA_d_fit = FA_d_fit + fabs(FA_d_fit.ave()-FA_d_syst_fit.ave())*erf(k_erf*fabs((FA_d_fit.ave()-FA_d_syst_fit.ave()))/(sqrt(2.0)*FA_d_syst_fit.err()))*(FA_d_fit - FA_d_fit.ave())/(FA_d_fit.err());
      FA_d.distr_list.push_back(FA_d_fit);
      FA_d_per_kin[ixg-1].distr_list.push_back(FA_d_fit);
      ax_d_Tmin.push_back( Corr.Tmin); ax_d_Tmax.push_back( Corr.Tmax);
      ax_Tmin.push_back( Corr.Tmin); ax_Tmax.push_back( Corr.Tmax);
      FA.distr_list.push_back( FA_u_fit+ FA_d_fit);
      FA_per_kin[ixg-1].distr_list.push_back( FA_u_fit + FA_d_fit);

      //FV_u
      Get_Bs_Mh_Tmin_Tmax("FVu",MESON, Tmin_V, Tmax_V, ixg, data_2pts.Tag[iens]);
      Corr.Tmin = Tmin_V; Corr.Tmax= Tmax_V;
      distr_t FV_u_fit= Corr.Fit_distr(FV_u_new_distr);
      Corr.Tmin += Dx; Corr.Tmax += Dx;
      distr_t FV_u_syst_fit= Corr.Fit_distr(FV_u_new_distr);
      Corr.Tmin -= Dx; Corr.Tmax -= Dx;
      FV_u_fit = FV_u_fit + fabs(FV_u_fit.ave()-FV_u_syst_fit.ave())*erf(k_erf*fabs((FV_u_fit.ave()-FV_u_syst_fit.ave()))/(sqrt(2.0)*FV_u_syst_fit.err()))*(FV_u_fit - FV_u_fit.ave())/(FV_u_fit.err());
      FV_u.distr_list.push_back(FV_u_fit);
      FV_u_per_kin[ixg-1].distr_list.push_back(FV_u_fit);
      vec_u_Tmin.push_back( Corr.Tmin); vec_u_Tmax.push_back( Corr.Tmax);
      //FV_d
      Get_Bs_Mh_Tmin_Tmax("FVd",MESON, Tmin_V, Tmax_V, ixg, data_2pts.Tag[iens]);
      Corr.Tmin = Tmin_V; Corr.Tmax= Tmax_V;
      distr_t FV_d_fit= Corr.Fit_distr(FV_d_new_distr);
      Corr.Tmin += Dx; Corr.Tmax += Dx;
      distr_t FV_d_syst_fit= Corr.Fit_distr(FV_d_new_distr);
      Corr.Tmin -= Dx; Corr.Tmax -= Dx;
      FV_d_fit = FV_d_fit + fabs(FV_d_fit.ave()-FV_d_syst_fit.ave())*erf(k_erf*fabs((FV_d_fit.ave()-FV_d_syst_fit.ave()))/(sqrt(2.0)*FV_d_syst_fit.err()))*(FV_d_fit - FV_d_fit.ave())/(FV_d_fit.err());
      FV_d.distr_list.push_back(FV_d_fit);
      FV_d_per_kin[ixg-1].distr_list.push_back(FV_d_fit);
      vec_d_Tmin.push_back( Corr.Tmin); vec_d_Tmax.push_back( Corr.Tmax);
      vec_Tmin.push_back( Corr.Tmin); vec_Tmax.push_back( Corr.Tmax);
      FV.distr_list.push_back( FV_u_fit+ FV_d_fit);
      FV_per_kin[ixg-1].distr_list.push_back( FV_u_fit + FV_d_fit);
      }





      //######################################################################################################################



      //######################################################################################################################
      //####################     FV_T     AND    FA_T     COMPUTATION     ###########################

      //push_back Ax_tens and Vec_tens
      Ax_T_glb.push_back(Ax_T_tens);
      Ax_T_u_glb.push_back(Ax_T_tens_u);
      Ax_T_d_glb.push_back(Ax_T_tens_d);
      
      Vec_T_glb.push_back(Vec_T_tens);
      Vec_T_u_glb.push_back(Vec_T_tens_u);
      Vec_T_d_glb.push_back(Vec_T_tens_d);


      distr_t_list TA0_u_distr_11= Ax_T_u_glb[0][1-off_i][1-off_T];
      distr_t_list TA0_u_distr_22= Ax_T_u_glb[0][2-off_i][2-off_T];
      distr_t_list TA0_d_distr_11= Ax_T_d_glb[0][1-off_i][1-off_T];
      distr_t_list TA0_d_distr_22= Ax_T_d_glb[0][2-off_i][2-off_T];

      //improved estimator (zero-momentum-subtracted FB)

      distr_t_list Ax_T_tens_u_impr_11 = Ax_T_tens_u[1-off_i][1-off_T] -  0.0*1.0/((1.0/TA0_u_distr_11)*EXP_PH);
      distr_t_list Ax_T_tens_u_impr_22 = Ax_T_tens_u[2-off_i][2-off_T] -  0.0*1.0/((1.0/TA0_u_distr_22)*EXP_PH);

      distr_t_list Ax_T_tens_d_impr_11 = Ax_T_tens_d[1-off_i][1-off_T] - 0.0*1.0/((1.0/TA0_d_distr_11)*EXP_PH) ;
      distr_t_list Ax_T_tens_d_impr_22 = Ax_T_tens_d[2-off_i][2-off_T] - 0.0*1.0/((1.0/TA0_d_distr_22)*EXP_PH);
    
      //Compute FV and FA
      distr_t_list FA_T_distr =0.5*RF*(Z_T/Zv)*(FP_SM/(-1.0*FA0_distr))*( Ax_T_tens[1-off_i][1-off_T]*(1-xg/2) - sign_kz*Vec_T_tens[1-off_i][2-off_T]*(xg/2) + Ax_T_tens[2-off_i][2-off_T]*(1-xg/2) +sign_kz*Vec_T_tens[2-off_i][1-off_T]*(xg/2))*(1.0/Eg)*EXP_PH;
      distr_t_list FV_T_distr = 0.5*RF*(Z_T/Zv)*(FP_SM/(-1.0*FA0_distr))*( (Vec_T_tens[1-off_i][2-off_T]*(1-xg/2) -sign_kz*Ax_T_tens[1-off_i][1-off_T]*(xg/2)) -(Vec_T_tens[2-off_i][1-off_i]*(1-xg/2) + sign_kz*Ax_T_tens[2-off_i][2-off_T]*(xg/2)))*(1.0/Eg)*EXP_PH;
      //improved estimator
      distr_t_list FA_T_new_distr = (0.5*RF*(Z_T)*( Ax_T_tens[1-off_i][1-off_T]*(1-xg/2) - sign_kz*Vec_T_tens[1-off_i][2-off_T]*(xg/2) + Ax_T_tens[2-off_i][2-off_T]*(1-xg/2) +sign_kz*Vec_T_tens[2-off_i][1-off_T]*(xg/2))*(1.0/(mel_SMSM*Eg))*EXP_PH)*EXP_MES;
      distr_t_list FV_T_new_distr = (0.5*RF*(Z_T)*( (Vec_T_tens[1-off_i][2-off_T]*(1-xg/2) -sign_kz*Ax_T_tens[1-off_i][1-off_T]*(xg/2)) -(Vec_T_tens[2-off_i][1-off_i]*(1-xg/2) + sign_kz*Ax_T_tens[2-off_i][2-off_T]*(xg/2)))*(1.0/(mel_SMSM*Eg))*EXP_PH)*EXP_MES;
      
      //compute FV and FA (up-component)
      distr_t_list FA_T_u_distr =0.5*RF*(Z_T/Zv)*(FP_SM/(-1.0*FA0_distr))*( Ax_T_tens_u_impr_11*(1-xg/2) - sign_kz*Vec_T_tens_u[1-off_i][2-off_T]*(xg/2) + Ax_T_tens_u_impr_22*(1-xg/2) +sign_kz*Vec_T_tens_u[2-off_i][1-off_T]*(xg/2))*(1.0/Eg)*EXP_PH;
      distr_t_list FV_T_u_distr = 0.5*RF*(Z_T/Zv)*(FP_SM/(-1.0*FA0_distr))*( (Vec_T_tens_u[1-off_i][2-off_T]*(1-xg/2) -sign_kz*Ax_T_tens_u_impr_11*(xg/2)) -(Vec_T_tens_u[2-off_i][1-off_i]*(1-xg/2) + sign_kz*Ax_T_tens_u_impr_22*(xg/2)))*(1.0/Eg)*EXP_PH;
   

      
      distr_t_list FA_T_u_new_distr =(0.5*RF*(Z_T)*( Ax_T_tens_u_impr_11*(1-xg/2) - sign_kz*Vec_T_tens_u[1-off_i][2-off_T]*(xg/2) + Ax_T_tens_u_impr_22*(1-xg/2) +sign_kz*Vec_T_tens_u[2-off_i][1-off_T]*(xg/2))*(1.0/(mel_SMSM*Eg))*EXP_PH)*EXP_MES;
      distr_t_list FV_T_u_new_distr = (0.5*RF*(Z_T)*( (Vec_T_tens_u[1-off_i][2-off_T]*(1-xg/2) -sign_kz*Ax_T_tens_u_impr_11*(xg/2)) -(Vec_T_tens_u[2-off_i][1-off_i]*(1-xg/2) + sign_kz*Ax_T_tens_u_impr_22*(xg/2)))*(1.0/(mel_SMSM*Eg))*EXP_PH)*EXP_MES;
      //another variation
      distr_t_list FA_T_u_var =(0.5*RF*(Z_T)*( Ax_T_tens_u_impr_11*(1-xg/2) - sign_kz*Vec_T_tens_u[1-off_i][2-off_T]*(xg/2) + Ax_T_tens_u_impr_22*(1-xg/2) +sign_kz*Vec_T_tens_u[2-off_i][1-off_T]*(xg/2))*(1.0/(mel_SMSM*Eg))*EXP_PH)*EXP_MASS_CH[4];
      distr_t_list FV_T_u_var = (0.5*RF*(Z_T)*( (Vec_T_tens_u[1-off_i][2-off_T]*(1-xg/2) -sign_kz*Ax_T_tens_u_impr_11*(xg/2)) -(Vec_T_tens_u[2-off_i][1-off_i]*(1-xg/2) + sign_kz*Ax_T_tens_u_impr_22*(xg/2)))*(1.0/(mel_SMSM*Eg))*EXP_PH)*EXP_MASS_CH[5];
      

      //compute FV and FA (d-component)
      distr_t_list FA_T_d_distr =0.5*RF*(Z_T/Zv)*(FP_SM/(-1.0*FA0_distr))*( Ax_T_tens_d_impr_11*(1-xg/2) - sign_kz*Vec_T_tens_d[1-off_i][2-off_T]*(xg/2) + Ax_T_tens_d_impr_22*(1-xg/2) +sign_kz*Vec_T_tens_d[2-off_i][1-off_T]*(xg/2))*(1.0/Eg)*EXP_PH;
      distr_t_list FV_T_d_distr = 0.5*RF*(Z_T/Zv)*(FP_SM/(-1.0*FA0_distr))*( (Vec_T_tens_d[1-off_i][2-off_T]*(1-xg/2) -sign_kz*Ax_T_tens_d_impr_11*(xg/2)) -(Vec_T_tens_d[2-off_i][1-off_i]*(1-xg/2) + sign_kz*Ax_T_tens_d_impr_22*(xg/2)))*(1.0/Eg)*EXP_PH;
      //improved estimator
      distr_t_list FA_T_d_new_distr = (0.5*RF*(Z_T)*( Ax_T_tens_d_impr_11*(1-xg/2) - sign_kz*Vec_T_tens_d[1-off_i][2-off_T]*(xg/2) + Ax_T_tens_d_impr_22*(1-xg/2) +sign_kz*Vec_T_tens_d[2-off_i][1-off_T]*(xg/2))*(1.0/(mel_SMSM*Eg))*EXP_PH)*EXP_MES;
      distr_t_list FV_T_d_new_distr = (0.5*RF*(Z_T)*( (Vec_T_tens_d[1-off_i][2-off_T]*(1-xg/2) -sign_kz*Ax_T_tens_d_impr_11*(xg/2)) -(Vec_T_tens_d[2-off_i][1-off_i]*(1-xg/2) + sign_kz*Ax_T_tens_d_impr_22*(xg/2)))*(1.0/(mel_SMSM*Eg))*EXP_PH)*EXP_MES;
      //another variation
      distr_t_list FA_T_d_var =(0.5*RF*(Z_T)*( Ax_T_tens_d_impr_11*(1-xg/2) - sign_kz*Vec_T_tens_d[1-off_i][2-off_T]*(xg/2) + Ax_T_tens_d_impr_22*(1-xg/2) +sign_kz*Vec_T_tens_d[2-off_i][1-off_T]*(xg/2))*(1.0/(mel_SMSM*Eg))*EXP_PH)*EXP_MASS_CH[6];
      distr_t_list FV_T_d_var = (0.5*RF*(Z_T)*( (Vec_T_tens_d[1-off_i][2-off_T]*(1-xg/2) -sign_kz*Ax_T_tens_d_impr_11*(xg/2)) -(Vec_T_tens_d[2-off_i][1-off_i]*(1-xg/2) + sign_kz*Ax_T_tens_d_impr_22*(xg/2)))*(1.0/(mel_SMSM*Eg))*EXP_PH)*EXP_MASS_CH[7];
   
      //Print FV and FA
      Print_To_File({}, {FA_T_distr.ave(), FA_T_distr.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"FF/"+data_2pts.Tag[iens]+"/FA_T_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "");
      Print_To_File({}, {FV_T_distr.ave(), FV_T_distr.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"FF/"+data_2pts.Tag[iens]+"/FV_T_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "");
      //Print_To_File({}, {FV_sub_distr.ave(), FV_sub_distr.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"FF/"+data_2pts.Tag[iens]+"/FV_T_sub_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "");

      //Print FV and FA (up-component)
      Print_To_File({}, {FA_T_u_distr.ave(), FA_T_u_distr.err(), FA_T_u_new_distr.ave(), FA_T_u_new_distr.err(), FA_T_u_var.ave(), FA_T_u_var.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"FF_u/"+data_2pts.Tag[iens]+"/FA_T_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "");
      Print_To_File({}, {FV_T_u_distr.ave(), FV_T_u_distr.err(), FV_T_u_new_distr.ave(), FV_T_u_new_distr.err(), FV_T_u_var.ave(), FV_T_u_var.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"FF_u/"+data_2pts.Tag[iens]+"/FV_T_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "");
      //Print_To_File({}, {FV_T_u_sub_distr.ave(), FV_T_u_sub_distr.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"FF_u/"+data_2pts.Tag[iens]+"/FV_T_sub_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "");
      
      //Print FV and FA (d-component)
      Print_To_File({}, {FA_T_d_distr.ave(), FA_T_d_distr.err(), FA_T_d_new_distr.ave(), FA_T_d_new_distr.err(), FA_T_d_var.ave(), FA_T_d_var.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"FF_d/"+data_2pts.Tag[iens]+"/FA_T_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "");
      Print_To_File({}, {(FV_T_d_distr).ave(), (FV_T_d_distr).err(), FV_T_d_new_distr.ave(), FV_T_d_new_distr.err(), FV_T_d_var.ave(), FV_T_d_var.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"FF_d/"+data_2pts.Tag[iens]+"/FV_T_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "");
      // Print_To_File({}, {(FV_T_d_sub_distr).ave(), (FV_T_d_sub_distr).err()}, "../data/ph_emission/"+ph_type_mes+"/"+"FF_d/"+data_2pts.Tag[iens]+"/FV_T_sub_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "");

      //fit the form factors

      //set interval depending on xg
      int Tmin_T_V, Tmax_T_V, Tmin_T_A, Tmax_T_A;
      
      if(xg.ave() > 1e-10) {
      //set time intervals for each heavy mass
      //FA_T_u
      Get_Bs_Mh_Tmin_Tmax("FA_Tu",MESON, Tmin_T_A, Tmax_T_A, ixg, data_2pts.Tag[iens]);
      Corr.Tmin = Tmin_T_A; Corr.Tmax= Tmax_T_A;
      distr_t FA_T_u_fit= Corr.Fit_distr(FA_T_u_new_distr);
      Corr.Tmin += Dx; Corr.Tmax += Dx;
      distr_t FA_T_u_syst_fit= Corr.Fit_distr(FA_T_u_new_distr);
      Corr.Tmin -= Dx; Corr.Tmax -= Dx;
      FA_T_u_fit = FA_T_u_fit + fabs(FA_T_u_fit.ave()-FA_T_u_syst_fit.ave())*erf(k_erf*fabs((FA_T_u_fit.ave()-FA_T_u_syst_fit.ave()))/(sqrt(2.0)*FA_T_u_syst_fit.err()))*(FA_T_u_fit - FA_T_u_fit.ave())/(FA_T_u_fit.err());
      FA_T_u.distr_list.push_back(FA_T_u_fit);
      FA_T_u_per_kin[ixg-1].distr_list.push_back(FA_T_u_fit);
      ax_T_u_Tmin.push_back( Corr.Tmin); ax_T_u_Tmax.push_back( Corr.Tmax);
      //FA_T_d
      Get_Bs_Mh_Tmin_Tmax("FA_Td", MESON, Tmin_T_A, Tmax_T_A, ixg, data_2pts.Tag[iens]);
      Corr.Tmin = Tmin_T_A; Corr.Tmax= Tmax_T_A;
      distr_t FA_T_d_fit= Corr.Fit_distr(FA_T_d_new_distr);
      Corr.Tmin += Dx; Corr.Tmax += Dx;
      distr_t FA_T_d_syst_fit= Corr.Fit_distr(FA_T_d_new_distr);
      Corr.Tmin -= Dx; Corr.Tmax -= Dx;
      FA_T_d_fit = FA_T_d_fit + fabs(FA_T_d_fit.ave()-FA_T_d_syst_fit.ave())*erf(k_erf*fabs((FA_T_d_fit.ave()-FA_T_d_syst_fit.ave()))/(sqrt(2.0)*FA_T_d_syst_fit.err()))*(FA_T_d_fit - FA_T_d_fit.ave())/(FA_T_d_fit.err());
      FA_T_d.distr_list.push_back(FA_T_d_fit);
      FA_T_d_per_kin[ixg-1].distr_list.push_back(FA_T_d_fit);
      ax_T_d_Tmin.push_back( Corr.Tmin); ax_T_d_Tmax.push_back( Corr.Tmax);
      ax_T_Tmin.push_back( Corr.Tmin); ax_T_Tmax.push_back( Corr.Tmax);
      FA_T.distr_list.push_back( FA_T_u_fit+ FA_T_d_fit);
      FA_T_per_kin[ixg-1].distr_list.push_back( FA_T_u_fit + FA_T_d_fit);

      //FV_T_u
      Get_Bs_Mh_Tmin_Tmax("FV_Tu", MESON, Tmin_T_V, Tmax_T_V, ixg, data_2pts.Tag[iens]);
      Corr.Tmin = Tmin_T_V; Corr.Tmax= Tmax_T_V;
      distr_t FV_T_u_fit= Corr.Fit_distr(FV_T_u_new_distr);
      Corr.Tmin += Dx; Corr.Tmax += Dx;
      distr_t FV_T_u_syst_fit= Corr.Fit_distr(FV_T_u_new_distr);
      Corr.Tmin -= Dx; Corr.Tmax -= Dx;
      FV_T_u_fit = FV_T_u_fit + fabs(FV_T_u_fit.ave()-FV_T_u_syst_fit.ave())*erf(k_erf*fabs((FV_T_u_fit.ave()-FV_T_u_syst_fit.ave()))/(sqrt(2.0)*FV_T_u_syst_fit.err()))*(FV_T_u_fit - FV_T_u_fit.ave())/(FV_T_u_fit.err());
      FV_T_u.distr_list.push_back(FV_T_u_fit);
      FV_T_u_per_kin[ixg-1].distr_list.push_back(FV_T_u_fit);
      vec_T_u_Tmin.push_back( Corr.Tmin); vec_T_u_Tmax.push_back( Corr.Tmax);
      //FV_T_d
      Get_Bs_Mh_Tmin_Tmax("FV_Td", MESON, Tmin_T_V, Tmax_T_V, ixg, data_2pts.Tag[iens]);
      Corr.Tmin = Tmin_T_V; Corr.Tmax= Tmax_T_V;
      distr_t FV_T_d_fit= Corr.Fit_distr(FV_T_d_new_distr);
      Corr.Tmin += Dx; Corr.Tmax += Dx;
      distr_t FV_T_d_syst_fit= Corr.Fit_distr(FV_T_d_new_distr);
      Corr.Tmin -= Dx; Corr.Tmax -= Dx;
      FV_T_d_fit = FV_T_d_fit + fabs(FV_T_d_fit.ave()-FV_T_d_syst_fit.ave())*erf(k_erf*fabs((FV_T_d_fit.ave()-FV_T_d_syst_fit.ave()))/(sqrt(2.0)*FV_T_d_syst_fit.err()))*(FV_T_d_fit - FV_T_d_fit.ave())/(FV_T_d_fit.err());
      FV_T_d.distr_list.push_back(FV_T_d_fit);
      FV_T_d_per_kin[ixg-1].distr_list.push_back(FV_T_d_fit);
      vec_T_d_Tmin.push_back( Corr.Tmin); vec_T_d_Tmax.push_back( Corr.Tmax);
      vec_T_Tmin.push_back( Corr.Tmin); vec_T_Tmax.push_back( Corr.Tmax);
      FV_T.distr_list.push_back( FV_T_u_fit+ FV_T_d_fit);
      FV_T_per_kin[ixg-1].distr_list.push_back( FV_T_u_fit + FV_T_d_fit);


      //FT and FB
      FT_u.distr_list.push_back( (1-xg/2.0)*FV_T_u_fit -sign_kz*FA_T_u_fit*(xg/2.0));
      FB_u.distr_list.push_back( (1-xg/2.0)*FA_T_u_fit -sign_kz*FV_T_u_fit*(xg/2.0) );
      FT_d.distr_list.push_back( (1-xg/2.0)*FV_T_d_fit -sign_kz*FA_T_d_fit*(xg/2.0));
      FB_d.distr_list.push_back( (1-xg/2.0)*FA_T_d_fit -sign_kz*FV_T_d_fit*(xg/2.0));
      FT.distr_list.push_back( (1-xg/2.0)*(FV_T_u_fit+FV_T_d_fit) - sign_kz*(FA_T_u_fit+FA_T_d_fit)*(xg/2.0));
      FB.distr_list.push_back( (1-xg/2.0)*(FA_T_u_fit+FA_T_d_fit) - sign_kz*(FV_T_u_fit+FV_T_d_fit)*(xg/2.0));

      FT_u_per_kin[ixg-1].distr_list.push_back( (1-xg/2.0)*FV_T_u_fit -sign_kz*FA_T_u_fit*(xg/2.0));
      FB_u_per_kin[ixg-1].distr_list.push_back( (1-xg/2.0)*FA_T_u_fit -sign_kz*FV_T_u_fit*(xg/2.0) );
      FT_d_per_kin[ixg-1].distr_list.push_back( (1-xg/2.0)*FV_T_d_fit -sign_kz*FA_T_d_fit*(xg/2.0));
      FB_d_per_kin[ixg-1].distr_list.push_back( (1-xg/2.0)*FA_T_d_fit -sign_kz*FV_T_d_fit*(xg/2.0));
      FT_per_kin[ixg-1].distr_list.push_back( (1-xg/2.0)*(FV_T_u_fit+FV_T_d_fit) - sign_kz*(FA_T_u_fit+FA_T_d_fit)*(xg/2.0));
      FB_per_kin[ixg-1].distr_list.push_back( (1-xg/2.0)*(FA_T_u_fit+FA_T_d_fit) - sign_kz*(FV_T_u_fit+FV_T_d_fit)*(xg/2.0));

      
      
         
      }

      




      //######################################################################################################################
      
      
    
      
    }

    //for FV and FA
    FA_per_ens[iens] = FA;
    FV_per_ens[iens] = FV;
    FA_u_per_ens[iens] = FA_u;
    FA_d_per_ens[iens] = FA_d;
    FV_u_per_ens[iens] = FV_u;
    FV_d_per_ens[iens] = FV_d;

    //for FV_T and FA_T
    FA_T_per_ens[iens] = FA_T;
    FV_T_per_ens[iens] = FV_T;
    FA_T_u_per_ens[iens] = FA_T_u;
    FA_T_d_per_ens[iens] = FA_T_d;
    FV_T_u_per_ens[iens] = FV_T_u;
    FV_T_d_per_ens[iens] = FV_T_d;

    //for FT and FB
    FB_per_ens[iens] = FB;
    FT_per_ens[iens] = FT;
    FB_u_per_ens[iens] = FB_u;
    FB_d_per_ens[iens] = FB_d;
    FT_u_per_ens[iens] = FT_u;
    FT_d_per_ens[iens] = FT_d;


    //Print fitted form factors
    //FV and FA
    Print_To_File({}, {xg_list.ave(), xg_list.err(), FV.ave(), FV.err(), vec_Tmin, vec_Tmax}, "../data/ph_emission/"+ph_type_mes+"/FF/"+data_2pts.Tag[iens]+"/fit_results/FV.dat", "", "#xg Dxg  FV DFV Tmin Tmax");
    Print_To_File({}, {xg_list.ave(), xg_list.err(), FA.ave(), FA.err(), ax_Tmin, ax_Tmax}, "../data/ph_emission/"+ph_type_mes+"/FF/"+data_2pts.Tag[iens]+"/fit_results/FA.dat", "", "#xg Dxg  FA DFA Tmin Tmax");

    Print_To_File({}, {xg_list.ave(), xg_list.err(), FV_u.ave(), FV_u.err(), vec_u_Tmin, vec_u_Tmax}, "../data/ph_emission/"+ph_type_mes+"/FF_u/"+data_2pts.Tag[iens]+"/fit_results/FV.dat", "", "#xg Dxg  FV DFV");
    Print_To_File({}, {xg_list.ave(), xg_list.err(), FA_u.ave(), FA_u.err(), ax_u_Tmin, ax_u_Tmax}, "../data/ph_emission/"+ph_type_mes+"/FF_u/"+data_2pts.Tag[iens]+"/fit_results/FA.dat", "", "#xg Dxg  FA DFA");

    Print_To_File({}, {xg_list.ave(), xg_list.err(), (FV_d).ave(), (FV_d).err(), vec_d_Tmin, vec_d_Tmax}, "../data/ph_emission/"+ph_type_mes+"/FF_d/"+data_2pts.Tag[iens]+"/fit_results/FV.dat", "", "#xg Dxg  FV DFV");
    Print_To_File({}, {xg_list.ave(), xg_list.err(), FA_d.ave(), FA_d.err(), ax_d_Tmin, ax_d_Tmax}, "../data/ph_emission/"+ph_type_mes+"/FF_d/"+data_2pts.Tag[iens]+"/fit_results/FA.dat", "", "#xg Dxg  FA DFA");


    //FV_T and FA_T
    Print_To_File({}, {xg_list.ave(), xg_list.err(), FV_T.ave(), FV_T.err(), vec_T_Tmin, vec_T_Tmax}, "../data/ph_emission/"+ph_type_mes+"/FF/"+data_2pts.Tag[iens]+"/fit_results/FV_T.dat", "", "#xg Dxg  FV DFV Tmin Tmax");
    Print_To_File({}, {xg_list.ave(), xg_list.err(), FA_T.ave(), FA_T.err(), ax_T_Tmin, ax_T_Tmax}, "../data/ph_emission/"+ph_type_mes+"/FF/"+data_2pts.Tag[iens]+"/fit_results/FA_T.dat", "", "#xg Dxg  FA DFA Tmin Tmax");

    Print_To_File({}, {xg_list.ave(), xg_list.err(), FV_T_u.ave(), FV_T_u.err(), vec_T_u_Tmin, vec_T_u_Tmax}, "../data/ph_emission/"+ph_type_mes+"/FF_u/"+data_2pts.Tag[iens]+"/fit_results/FV_T.dat", "", "#xg Dxg  FV DFV Tmin Tmax");
    Print_To_File({}, {xg_list.ave(), xg_list.err(), FA_T_u.ave(), FA_T_u.err(), ax_T_u_Tmin, ax_T_u_Tmax}, "../data/ph_emission/"+ph_type_mes+"/FF_u/"+data_2pts.Tag[iens]+"/fit_results/FA_T.dat", "", "#xg Dxg  FA DFA Tmin Tmax");

    Print_To_File({}, {xg_list.ave(), xg_list.err(), (FV_T_d).ave(), (FV_T_d).err(), vec_T_d_Tmin, vec_T_d_Tmax}, "../data/ph_emission/"+ph_type_mes+"/FF_d/"+data_2pts.Tag[iens]+"/fit_results/FV_T.dat", "", "#xg Dxg  FV DFV Tmin Tmax");
    Print_To_File({}, {xg_list.ave(), xg_list.err(), FA_T_d.ave(), FA_T_d.err(), ax_T_d_Tmin, ax_T_d_Tmax}, "../data/ph_emission/"+ph_type_mes+"/FF_d/"+data_2pts.Tag[iens]+"/fit_results/FA_T.dat", "", "#xg Dxg  FA DFA Tmin Tmax");


    //FT and FB
    Print_To_File({}, {xg_list.ave(), xg_list.err(), FT.ave(), FT.err()}, "../data/ph_emission/"+ph_type_mes+"/FF/"+data_2pts.Tag[iens]+"/fit_results/FT.dat", "", "#xg Dxg  FV DFV");
    Print_To_File({}, {xg_list.ave(), xg_list.err(), FB.ave(), FB.err()}, "../data/ph_emission/"+ph_type_mes+"/FF/"+data_2pts.Tag[iens]+"/fit_results/FB.dat", "", "#xg Dxg  FA DFA");
    Print_To_File({}, {xg_list.ave(), xg_list.err(), FT_u.ave(), FT_u.err()}, "../data/ph_emission/"+ph_type_mes+"/FF_u/"+data_2pts.Tag[iens]+"/fit_results/FT.dat", "", "#xg Dxg  FV DFV");
    Print_To_File({}, {xg_list.ave(), xg_list.err(), FB_u.ave(), FB_u.err()}, "../data/ph_emission/"+ph_type_mes+"/FF_u/"+data_2pts.Tag[iens]+"/fit_results/FB.dat", "", "#xg Dxg  FA DFA");
    Print_To_File({}, {xg_list.ave(), xg_list.err(), (FT_d).ave(), (FT_d).err()}, "../data/ph_emission/"+ph_type_mes+"/FF_d/"+data_2pts.Tag[iens]+"/fit_results/FT.dat", "", "#xg Dxg  FV DFV");
    Print_To_File({}, {xg_list.ave(), xg_list.err(), FB_d.ave(), FB_d.err()}, "../data/ph_emission/"+ph_type_mes+"/FF_d/"+data_2pts.Tag[iens]+"/fit_results/FB.dat", "", "#xg Dxg  FA DFA");


  }


  


  if(num_xg==1) crash("Exiting...No interpolation nor continuum extrapolation to perform, only xg available is zero");


  //Print all ensembles for given xg
  for(int ixg=1;ixg<num_xg;ixg++) {

    //for FA and FV
    Print_To_File({}, {(a_distr_list/fmTGeV).ave(), FA_per_kin[ixg-1].ave(), FA_per_kin[ixg-1].err(), L_list}, "../data/ph_emission/"+ph_type_mes+"/FF/per_kin/FA_ixg_"+to_string(ixg)+".dat", "", "#a[fm]  FA FA_err L/a");
    Print_To_File({}, {(a_distr_list/fmTGeV).ave(), FV_per_kin[ixg-1].ave(), FV_per_kin[ixg-1].err(), L_list}, "../data/ph_emission/"+ph_type_mes+"/FF/per_kin/FV_ixg_"+to_string(ixg)+".dat", "", "#a[fm]  FA FA_err L/a");

    //for FA_T and FV_T
    Print_To_File({}, {(a_distr_list/fmTGeV).ave(), FA_T_per_kin[ixg-1].ave(), FA_T_per_kin[ixg-1].err(), L_list}, "../data/ph_emission/"+ph_type_mes+"/FF/per_kin/FA_T_ixg_"+to_string(ixg)+".dat", "", "#a[fm]  FA FA_err L/a");
    Print_To_File({}, {(a_distr_list/fmTGeV).ave(), FV_T_per_kin[ixg-1].ave(), FV_T_per_kin[ixg-1].err(), L_list}, "../data/ph_emission/"+ph_type_mes+"/FF/per_kin/FV_T_ixg_"+to_string(ixg)+".dat", "", "#a[fm]  FA FA_err L/a");


    //for FB and FT
    Print_To_File({}, {(a_distr_list/fmTGeV).ave(), FB_per_kin[ixg-1].ave(), FB_per_kin[ixg-1].err(), L_list}, "../data/ph_emission/"+ph_type_mes+"/FF/per_kin/FB_ixg_"+to_string(ixg)+".dat", "", "#a[fm]  FA FA_err L/a");
    Print_To_File({}, {(a_distr_list/fmTGeV).ave(), FT_per_kin[ixg-1].ave(), FT_per_kin[ixg-1].err(), L_list}, "../data/ph_emission/"+ph_type_mes+"/FF/per_kin/FT_ixg_"+to_string(ixg)+".dat", "", "#a[fm]  FA FA_err L/a");


    //########### u contribution ############
    //for FA and FV
    Print_To_File({}, {(a_distr_list/fmTGeV).ave(), FA_u_per_kin[ixg-1].ave(), FA_u_per_kin[ixg-1].err(), L_list}, "../data/ph_emission/"+ph_type_mes+"/FF_u/per_kin/FA_ixg_"+to_string(ixg)+".dat", "", "#a[fm]  FA FA_err L/a");
    Print_To_File({}, {(a_distr_list/fmTGeV).ave(), FV_u_per_kin[ixg-1].ave(), FV_u_per_kin[ixg-1].err(), L_list}, "../data/ph_emission/"+ph_type_mes+"/FF_u/per_kin/FV_ixg_"+to_string(ixg)+".dat", "", "#a[fm]  FA FA_err L/a");

    //for FA_T and FV_T
    Print_To_File({}, {(a_distr_list/fmTGeV).ave(), FA_T_u_per_kin[ixg-1].ave(), FA_T_u_per_kin[ixg-1].err(), L_list}, "../data/ph_emission/"+ph_type_mes+"/FF_u/per_kin/FA_T_ixg_"+to_string(ixg)+".dat", "", "#a[fm]  FA FA_err L/a");
    Print_To_File({}, {(a_distr_list/fmTGeV).ave(), FV_T_u_per_kin[ixg-1].ave(), FV_T_u_per_kin[ixg-1].err(), L_list}, "../data/ph_emission/"+ph_type_mes+"/FF_u/per_kin/FV_T_ixg_"+to_string(ixg)+".dat", "", "#a[fm]  FA FA_err L/a");

    //for FB and FT
    Print_To_File({}, {(a_distr_list/fmTGeV).ave(), FB_u_per_kin[ixg-1].ave(), FB_u_per_kin[ixg-1].err(), L_list}, "../data/ph_emission/"+ph_type_mes+"/FF_u/per_kin/FB_ixg_"+to_string(ixg)+".dat", "", "#a[fm]  FA FA_err L/a");
    Print_To_File({}, {(a_distr_list/fmTGeV).ave(), FT_u_per_kin[ixg-1].ave(), FT_u_per_kin[ixg-1].err(), L_list}, "../data/ph_emission/"+ph_type_mes+"/FF_u/per_kin/FT_ixg_"+to_string(ixg)+".dat", "", "#a[fm]  FA FA_err L/a");


    //#######################################



    //########### d contribution ############
     //for FA and FV
    Print_To_File({}, {(a_distr_list/fmTGeV).ave(), FA_d_per_kin[ixg-1].ave(), FA_d_per_kin[ixg-1].err(), L_list}, "../data/ph_emission/"+ph_type_mes+"/FF_d/per_kin/FA_ixg_"+to_string(ixg)+".dat", "", "#a[fm]  FA FA_err L/a");
    Print_To_File({}, {(a_distr_list/fmTGeV).ave(), FV_d_per_kin[ixg-1].ave(), FV_d_per_kin[ixg-1].err(), L_list}, "../data/ph_emission/"+ph_type_mes+"/FF_d/per_kin/FV_ixg_"+to_string(ixg)+".dat", "", "#a[fm]  FA FA_err L/a");

    //for FA_T and FV_T
    Print_To_File({}, {(a_distr_list/fmTGeV).ave(), FA_T_d_per_kin[ixg-1].ave(), FA_T_d_per_kin[ixg-1].err(), L_list}, "../data/ph_emission/"+ph_type_mes+"/FF_d/per_kin/FA_T_ixg_"+to_string(ixg)+".dat", "", "#a[fm]  FA FA_err L/a");
    Print_To_File({}, {(a_distr_list/fmTGeV).ave(), FV_T_d_per_kin[ixg-1].ave(), FV_T_d_per_kin[ixg-1].err(), L_list}, "../data/ph_emission/"+ph_type_mes+"/FF_d/per_kin/FV_T_ixg_"+to_string(ixg)+".dat", "", "#a[fm]  FA FA_err L/a");

    //for FB and FT
    Print_To_File({}, {(a_distr_list/fmTGeV).ave(), FB_d_per_kin[ixg-1].ave(), FB_d_per_kin[ixg-1].err(), L_list}, "../data/ph_emission/"+ph_type_mes+"/FF_d/per_kin/FB_ixg_"+to_string(ixg)+".dat", "", "#a[fm]  FA FA_err L/a");
    Print_To_File({}, {(a_distr_list/fmTGeV).ave(), FT_d_per_kin[ixg-1].ave(), FT_d_per_kin[ixg-1].err(), L_list}, "../data/ph_emission/"+ph_type_mes+"/FF_d/per_kin/FT_ixg_"+to_string(ixg)+".dat", "", "#a[fm]  FA FA_err L/a");

    //#######################################



    
    
  }

  //################################################################################################

  //Print MP, FP, MP_ov_FP, phi
  Print_To_File({}, {a_distr_list_red.ave(), MP_list_red.ave(), MP_list_red.err()}, "../data/ph_emission/"+ph_type+"/"+MESON+"/mass/masses.list", "", "#a MP MP_err");
  Print_To_File({}, {a_distr_list.ave(), FP_list.ave(), FP_list.err(), FP_3pt_list.ave(), FP_3pt_list.err(),  FP_diml_list.ave(), FP_diml_list.err(), FP_3pt_diml_list.ave(), FP_3pt_diml_list.err(),  Tmin_fp, Tmax_fp, Tmin_fp_3pt, Tmax_fp_3pt}, "../data/ph_emission/"+ph_type+"/"+MESON+"/decay_const/fP.list", "", "#a FP FP_3pt FP_bare FP_bare_3pt Tmin Tmax  ");
  Print_To_File({}, {a_distr_list.ave(),  MP_ov_FP_list.ave(), MP_ov_FP_list.err(), Tmin_fp, Tmax_fp}, "../data/ph_emission/"+ph_type+"/"+MESON+"/decay_const/mP_ov_fP.list", "", "#a FP FP_err MP/FP MP/FP_err Tmin Tmax ");
  Print_To_File({}, {a_distr_list.ave(), (FP_list*SQRT_DL(MP_list)).ave(), (FP_list*SQRT_DL(MP_list)).err(),  (FP_3pt_list*SQRT_DL(MP_list)).ave(), (FP_3pt_list*SQRT_DL(MP_list)).err(),  Tmin_fp, Tmax_fp, Tmin_fp_3pt, Tmax_fp_3pt}, "../data/ph_emission/"+ph_type+"/"+MESON+"/decay_const/phi.list", "", "#a phi phi_list Tmin Tmax");


  //#################################################################################################   
  //continuum extrapolation


  
  if(Perform_continuum_extrapolation) {


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
    bf_FF.set_warmup_lev(0); //sets warmup
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



  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //###########################################                                                                 ##########################################
  //###########################################                                                                 ##########################################
  //###########################################                                                                 ##########################################
  //###########################################                                                                 ##########################################
  //###########################################                                                                 ##########################################
  //###########################################                                                                 ##########################################
  //###########################################                           FITTING FV AND FA                     ##########################################
  //###########################################                     INCLUDING SINGLE CONTRIBUTIONS              ##########################################
  //###########################################                        FROM STRANGE AND "CHARM"                 ##########################################
  //###########################################                                                                 ##########################################
  //###########################################                                                                 ##########################################
  //###########################################                                                                 ##########################################
  //###########################################                                                                 ##########################################
  //###########################################                                                                 ##########################################
  //###########################################                                                                 ##########################################
  //###########################################                                                                 ##########################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################




  
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //###########################################                                                                 ##########################################
  //###########################################                                                                 ##########################################
  //###########################################                   FIT FA for all xg                             ##########################################
  //###########################################                                                                 ##########################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################

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
  for(auto &a: Bs_a_to_print) FA_xg_to_print.distr_list.push_back( F0 + D1*pow(a*fmTGeV*Lambda_QCD,2) + D2*pow(a*fmTGeV*Lambda_QCD,4));
  Print_To_File({}, {Bs_a_to_print, FA_xg_to_print.ave(), FA_xg_to_print.err()}, "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FA_"+Fit_tag+"xg_"+to_string_with_precision(0.10*ixg,2)+".fit_func", "", "#a[fm] FA FA_err");
 
  }


  
    
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //###########################################                                                                 ##########################################
  //###########################################                                                                 ##########################################
  //###########################################                   FIT FV for all xg                             ##########################################
  //###########################################                                                                 ##########################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################

  
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
  for(auto &a: Bs_a_to_print) FV_xg_to_print.distr_list.push_back( F0 + D1*pow(a*fmTGeV*Lambda_QCD,2) + D2*pow(a*fmTGeV*Lambda_QCD,4));
  Print_To_File({}, {Bs_a_to_print, FV_xg_to_print.ave(), FV_xg_to_print.err()}, "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FV_"+Fit_tag+"xg_"+to_string_with_precision(0.10*ixg,2)+".fit_func", "", "#a[fm] FV FV_err");
  }

 

  //Print continuum extrapolated form factors
  Print_To_File({}, {Bs_xg_t_list, F0_A_list.ave(), F0_A_list.err(), (D1_A_list/F0_A_list).ave(), (D1_A_list/F0_A_list).err(), (D2_A_list/F0_A_list).ave(), (D2_A_list/F0_A_list).err(), ch2_FA}, "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FA_"+Fit_tag+"cont.dat", "", "#xg  F0   D1   D2  ch2/dof");
  Print_To_File({}, {Bs_xg_t_list, F0_V_list.ave(), F0_V_list.err(), (D1_V_list/F0_V_list).ave(), (D1_V_list/F0_V_list).err(), (D2_V_list/F0_V_list).ave(), (D2_V_list/F0_V_list).err(), ch2_FV}, "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FV_"+Fit_tag+"cont.dat", "", "#xg  F0   D1   D2  ch2/dof");


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
  ofstream Print_Cov_FA("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FA_"+Analysis_tag+".cov");
  ofstream Print_Cov_FV("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FV_"+Analysis_tag+".cov");
  ofstream Print_Corr_FA("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FA_"+Analysis_tag+".corr");
  ofstream Print_Corr_FV("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FV_"+Analysis_tag+".corr");

  Print_Cov_FA<<Cov_FA<<endl;
  Print_Cov_FV<<Cov_FV<<endl;
  Print_Corr_FA<<Corr_FA<<endl;
  Print_Corr_FV<<Corr_FV<<endl;


  Print_Cov_FA.close();
  Print_Cov_FV.close();
  Print_Corr_FA.close();
  Print_Corr_FV.close();


  
      
	

  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################

                                                  //############ Up-quark CONTRIBUTION ##################



  
    
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //###########################################                                                                 ##########################################
  //###########################################                                                                 ##########################################
  //###########################################                   FIT FA(u) for all xg                          ##########################################
  //###########################################                                                                 ##########################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################


  
  
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
  for(auto &a: Bs_a_to_print) FA_xg_to_print.distr_list.push_back( F0 + D1*pow(a*fmTGeV*Lambda_QCD,2) + D2*pow(a*fmTGeV*Lambda_QCD,4));
  Print_To_File({}, {Bs_a_to_print, FA_xg_to_print.ave(), FA_xg_to_print.err()}, "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FA_u_"+Fit_tag+"xg_"+to_string_with_precision(0.10*ixg,2)+".fit_func", "", "#a[fm] FA FA_err");
 
  }


  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //###########################################                                                                 ##########################################
  //###########################################                                                                 ##########################################
  //###########################################                   FIT FV(u) for all xg                          ##########################################
  //###########################################                                                                 ##########################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
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
  for(auto &a: Bs_a_to_print) FV_xg_to_print.distr_list.push_back( F0 + D1*pow(a*fmTGeV*Lambda_QCD,2) + D2*pow(a*fmTGeV*Lambda_QCD,4));
  Print_To_File({}, {Bs_a_to_print, FV_xg_to_print.ave(), FV_xg_to_print.err()}, "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FV_u_"+Fit_tag+"xg_"+to_string_with_precision(0.10*ixg,2)+".fit_func", "", "#a[fm] FV FV_err");
  }

  //Print continuum extrapolated form factors
  Print_To_File({}, {Bs_xg_t_list, F0_u_A_list.ave(), F0_u_A_list.err(), (D1_u_A_list/F0_u_A_list).ave(), (D1_u_A_list/F0_u_A_list).err(), (D2_u_A_list/F0_u_A_list).ave(), (D2_u_A_list/F0_u_A_list).err(), ch2_FA_u}, "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FA_u_"+Fit_tag+"cont.dat", "", "#xg  F0   D1   D2   ch2/dof");
  Print_To_File({}, {Bs_xg_t_list, F0_u_V_list.ave(), F0_u_V_list.err(), (D1_u_V_list/F0_u_V_list).ave(), (D1_u_V_list/F0_u_V_list).err(), (D2_u_V_list/F0_u_V_list).ave(), (D2_u_V_list/F0_u_V_list).err(), ch2_FV_u}, "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FV_u_"+Fit_tag+"cont.dat", "", "#xg  F0   D1   D2   ch2/dof");
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
  ofstream Print_Cov_FA_u("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FA_u_"+Analysis_tag+".cov");
  ofstream Print_Cov_FV_u("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FV_u_"+Analysis_tag+".cov");
  ofstream Print_Corr_FA_u("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FA_u_"+Analysis_tag+".corr");
  ofstream Print_Corr_FV_u("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FV_u_"+Analysis_tag+".corr");

  Print_Cov_FA_u<<Cov_FA_u<<endl;
  Print_Cov_FV_u<<Cov_FV_u<<endl;
  Print_Corr_FA_u<<Corr_FA_u<<endl;
  Print_Corr_FV_u<<Corr_FV_u<<endl;


  Print_Cov_FA_u.close();
  Print_Cov_FV_u.close();
  Print_Corr_FA_u.close();
  Print_Corr_FV_u.close();





  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################

                                                  //############ down-quark CONTRIBUTION ##################



  
    
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //###########################################                                                                 ##########################################
  //###########################################                                                                 ##########################################
  //###########################################                   FIT FA(d) for all xg                          ##########################################
  //###########################################                                                                 ##########################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
 
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
  for(auto &a: Bs_a_to_print) FA_xg_to_print.distr_list.push_back( F0 + D1*pow(a*fmTGeV*Lambda_QCD,2) + D2*pow(a*fmTGeV*Lambda_QCD,4));
  Print_To_File({}, {Bs_a_to_print, FA_xg_to_print.ave(), FA_xg_to_print.err()}, "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FA_d_"+Fit_tag+"xg_"+to_string_with_precision(0.10*ixg,2)+".fit_func", "", "#a[fm] FA FA_err");
 
  }


  
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //###########################################                                                                 ##########################################
  //###########################################                                                                 ##########################################
  //###########################################                   FIT FV(d) for all xg                          ##########################################
  //###########################################                                                                 ##########################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
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
  for(auto &a: Bs_a_to_print) FV_xg_to_print.distr_list.push_back( F0 + D1*pow(a*fmTGeV*Lambda_QCD,2) + D2*pow(a*fmTGeV*Lambda_QCD,4));
  Print_To_File({}, {Bs_a_to_print, FV_xg_to_print.ave(), FV_xg_to_print.err()}, "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FV_d_"+Fit_tag+"xg_"+to_string_with_precision(0.10*ixg,2)+".fit_func", "", "#a[fm] FV FV_err");
  }

  //Print continuum extrapolated form factors
  Print_To_File({}, {Bs_xg_t_list, F0_d_A_list.ave(), F0_d_A_list.err(), (D1_d_A_list/F0_d_A_list).ave(), (D1_d_A_list/F0_d_A_list).err(), (D2_d_A_list/F0_d_A_list).ave(), (D2_d_A_list/F0_d_A_list).err(), ch2_FA_d}, "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FA_d_"+Fit_tag+"cont.dat", "", "#xg  F0   D1   D2   ch2/dof");
  Print_To_File({}, {Bs_xg_t_list, F0_d_V_list.ave(), F0_d_V_list.err(), (D1_d_V_list/F0_d_V_list).ave(), (D1_d_V_list/F0_d_V_list).err(), (D2_d_V_list/F0_d_V_list).ave(), (D2_d_V_list/F0_d_V_list).err(), ch2_FV_d}, "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FV_d_"+Fit_tag+"cont.dat", "", "#xg  F0   D1   D2   ch2/dof");
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
  ofstream Print_Cov_FA_d("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FA_d_"+Analysis_tag+".cov");
  ofstream Print_Cov_FV_d("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FV_d_"+Analysis_tag+".cov");
  ofstream Print_Corr_FA_d("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FA_d_"+Analysis_tag+".corr");
  ofstream Print_Corr_FV_d("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FV_d_"+Analysis_tag+".corr");

  Print_Cov_FA_d<<Cov_FA_d<<endl;
  Print_Cov_FV_d<<Cov_FV_d<<endl;
  Print_Corr_FA_d<<Corr_FA_d<<endl;
  Print_Corr_FV_d<<Corr_FV_d<<endl;


  Print_Cov_FA_d.close();
  Print_Cov_FV_d.close();
  Print_Corr_FA_d.close();
  Print_Corr_FV_d.close();




  //###################################################################################################################################################################


  

  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //###########################################                                                                 ##########################################
  //###########################################                                                                 ##########################################
  //###########################################                                                                 ##########################################
  //###########################################                                                                 ##########################################
  //###########################################                                                                 ##########################################
  //###########################################                                                                 ##########################################
  //###########################################                           FITTING FV_T AND FA_T                  ##########################################
  //###########################################                     INCLUDING SINGLE CONTRIBUTIONS              ##########################################
  //###########################################                        FROM STRANGE AND "CHARM"                 ##########################################
  //###########################################                                                                 ##########################################
  //###########################################                                                                 ##########################################
  //###########################################                                                                 ##########################################
  //###########################################                                                                 ##########################################
  //###########################################                                                                 ##########################################
  //###########################################                                                                 ##########################################
  //###########################################                                                                 ##########################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################




  
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //###########################################                                                                 ##########################################
  //###########################################                                                                 ##########################################
  //###########################################                   FIT FA_T for all xg                           ##########################################
  //###########################################                                                                 ##########################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################

  Vfloat ch2_FA_T;
  distr_t_list F0_A_T_list(UseJack);
  distr_t_list D1_A_T_list(UseJack);
  distr_t_list D2_A_T_list(UseJack);
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
      data[ijack][iens].FF= FA_T_per_ens[iens].distr_list[ixg-1].distr[ijack];
      data[ijack][iens].FF_err= FA_T_per_ens[iens].err(ixg-1);
      if(data_2pts.Tag[iens] == "cA211a.12.48") { data[ijack][iens].a = a_A.distr[ijack]; data[ijack][iens].is=0; }
      else if(data_2pts.Tag[iens] == "cB211b.072.64") { data[ijack][iens].a = a_B.distr[ijack]; data[ijack][iens].is=1;}
      else if(data_2pts.Tag[iens] == "cB211b.072.96") { data[ijack][iens].a = a_B.distr[ijack]; data[ijack][iens].is=1; }
      else if(data_2pts.Tag[iens] == "cC211a.06.80")  {data[ijack][iens].a = a_C.distr[ijack]; data[ijack][iens].is=2; }
      else if(data_2pts.Tag[iens] == "cD211a.054.96")  {data[ijack][iens].a = a_D.distr[ijack]; data[ijack][iens].is=3; }
      else crash("Ens_tag: "+data_2pts.Tag[iens]+" not recognized");

      //mean values
      if(ijack==0) {
	data_ch2[ijack][iens].FF= FA_T_per_ens[iens].ave(ixg-1);
	data_ch2[ijack][iens].FF_err= FA_T_per_ens[iens].err(ixg-1);
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
  cout<<"Fitting FA_T, xg: "<<ixg<<endl;
  Bt_fit= bf_FF.Perform_bootstrap_fit();
  Bt_fit_ch2= bf_FF_ch2.Perform_bootstrap_fit();

  
  //retrieve parameters
  distr_t F0(UseJack), D1(UseJack), D2(UseJack);
  for(int ijack=0;ijack<Njacks;ijack++) { F0.distr.push_back( Bt_fit.par[ijack].F0); D1.distr.push_back( Bt_fit.par[ijack].D1); D2.distr.push_back( Bt_fit.par[ijack].D2);}
  //push_back retrieved parameters
  F0_A_T_list.distr_list.push_back(F0);
  D1_A_T_list.distr_list.push_back(D1);
  D2_A_T_list.distr_list.push_back(D2);
  //push_back ch2
  ch2_FA_T.push_back( Bt_fit_ch2.get_ch2_ave()/dof);

 
  //print fit func
  distr_t_list FA_T_xg_to_print(UseJack);
  for(auto &a: Bs_a_to_print) FA_T_xg_to_print.distr_list.push_back( F0 + D1*pow(a*fmTGeV*Lambda_QCD,2) + D2*pow(a*fmTGeV*Lambda_QCD,4));
  Print_To_File({}, {Bs_a_to_print, FA_T_xg_to_print.ave(), FA_T_xg_to_print.err()}, "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FA_T_"+Fit_tag+"xg_"+to_string_with_precision(0.10*ixg,2)+".fit_func", "", "#a[fm] FA FA_err");
 
  }


  
    
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //###########################################                                                                 ##########################################
  //###########################################                                                                 ##########################################
  //###########################################                   FIT FV_T for all xg                           ##########################################
  //###########################################                                                                 ##########################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################

  
  Vfloat ch2_FV_T;
  bf_FF.Set_par_val("F0", 0.1, 0.001);
  bf_FF.Set_par_val("D1", 1.0, 0.1);
  if(Include_a4) bf_FF.Set_par_val("D2", 1.0, 0.1);
  bf_FF_ch2.Set_par_val("F0", 0.1, 0.001);
  bf_FF_ch2.Set_par_val("D1", 1.0, 0.1);
  if(Include_a4) bf_FF_ch2.Set_par_val("D2", 1.0, 0.1);
  distr_t_list F0_V_T_list(UseJack);
  distr_t_list D1_V_T_list(UseJack);
  distr_t_list D2_V_T_list(UseJack);
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
	data[ijack][iens].FF= FV_T_per_ens[iens].distr_list[ixg-1].distr[ijack];
	data[ijack][iens].FF_err= FV_T_per_ens[iens].err(ixg-1);
	if(data_2pts.Tag[iens] == "cA211a.12.48") { data[ijack][iens].a = a_A.distr[ijack]; data[ijack][iens].is=0; }
	else if(data_2pts.Tag[iens] == "cB211b.072.64") { data[ijack][iens].a = a_B.distr[ijack]; data[ijack][iens].is=1;}
	else if(data_2pts.Tag[iens] == "cB211b.072.96") { data[ijack][iens].a = a_B.distr[ijack]; data[ijack][iens].is=1; }
	else if(data_2pts.Tag[iens] == "cC211a.06.80")  {data[ijack][iens].a = a_C.distr[ijack]; data[ijack][iens].is=2; }
	else if(data_2pts.Tag[iens] == "cD211a.054.96")  {data[ijack][iens].a = a_D.distr[ijack]; data[ijack][iens].is=3; }
	else crash("Ens_tag: "+data_2pts.Tag[iens]+" not recognized");

	if(ijack==0) {
	  	data_ch2[ijack][iens].FF= FV_T_per_ens[iens].ave(ixg-1);
		data_ch2[ijack][iens].FF_err= FV_T_per_ens[iens].err(ixg-1);
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
  cout<<"Fitting FV_T, xg: "<<ixg<<endl;
  Bt_fit= bf_FF.Perform_bootstrap_fit();
  Bt_fit_ch2= bf_FF_ch2.Perform_bootstrap_fit();
  //retrieve parameters
  distr_t F0(UseJack), D1(UseJack), D2(UseJack);
  for(int ijack=0;ijack<Njacks;ijack++) { F0.distr.push_back( Bt_fit.par[ijack].F0); D1.distr.push_back( Bt_fit.par[ijack].D1); D2.distr.push_back( Bt_fit.par[ijack].D2);}
  //push_back retrieved parameters
  F0_V_T_list.distr_list.push_back(F0);
  D1_V_T_list.distr_list.push_back(D1);
  D2_V_T_list.distr_list.push_back(D2);
  //push_back ch2
  ch2_FV_T.push_back( Bt_fit_ch2.get_ch2_ave()/dof);

  //print fit func
  distr_t_list FV_T_xg_to_print(UseJack);
  for(auto &a: Bs_a_to_print) FV_T_xg_to_print.distr_list.push_back( F0 + D1*pow(a*fmTGeV*Lambda_QCD,2) + D2*pow(a*fmTGeV*Lambda_QCD,4));
  Print_To_File({}, {Bs_a_to_print, FV_T_xg_to_print.ave(), FV_T_xg_to_print.err()}, "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FV_T_"+Fit_tag+"xg_"+to_string_with_precision(0.10*ixg,2)+".fit_func", "", "#a[fm] FV FV_err");
  }

 

  //Print continuum extrapolated form factors
  Print_To_File({}, {Bs_xg_t_list, F0_A_T_list.ave(), F0_A_T_list.err(), (D1_A_T_list/F0_A_T_list).ave(), (D1_A_T_list/F0_A_T_list).err(), (D2_A_T_list/F0_A_T_list).ave(), (D2_A_T_list/F0_A_T_list).err(), ch2_FA_T}, "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FA_T_"+Fit_tag+"cont.dat", "", "#xg  F0   D1   D2  ch2/dof");
  Print_To_File({}, {Bs_xg_t_list, F0_V_T_list.ave(), F0_V_T_list.err(), (D1_V_T_list/F0_V_T_list).ave(), (D1_V_T_list/F0_V_T_list).err(), (D2_V_T_list/F0_V_T_list).ave(), (D2_V_T_list/F0_V_T_list).err(), ch2_FV_T}, "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FV_T_"+Fit_tag+"cont.dat", "", "#xg  F0   D1   D2  ch2/dof");


  //Print covariance matrix
  Analysis_tag= ( (UseJack==true)?"jack":"boot")+_Fit_tag;
  Eigen::MatrixXd Cov_FA_T(num_xg-1, num_xg-1);
  Eigen::MatrixXd Cov_FV_T(num_xg-1, num_xg-1);
  Eigen::MatrixXd Corr_FA_T(num_xg-1, num_xg-1);
  Eigen::MatrixXd Corr_FV_T(num_xg-1, num_xg-1);

  for(int x_xg=1; x_xg<num_xg;x_xg++) {
    for(int y_xg=1; y_xg<num_xg;y_xg++) {
      Cov_FA_T(x_xg-1, y_xg-1) = F0_A_T_list.distr_list[x_xg-1]%F0_A_T_list.distr_list[y_xg-1];
      Cov_FV_T(x_xg-1, y_xg-1) = F0_V_T_list.distr_list[x_xg-1]%F0_V_T_list.distr_list[y_xg-1];
      Corr_FA_T(x_xg-1, y_xg-1) = Cov_FA_T(x_xg-1,y_xg-1)/(F0_A_T_list.err(x_xg-1)*F0_A_T_list.err(y_xg-1));
      Corr_FV_T(x_xg-1, y_xg-1) = Cov_FV_T(x_xg-1,y_xg-1)/(F0_V_T_list.err(x_xg-1)*F0_V_T_list.err(y_xg-1));
    }
  }

  //Print To File
  ofstream Print_Cov_FA_T("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FA_T_"+Analysis_tag+".cov");
  ofstream Print_Cov_FV_T("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FV_T_"+Analysis_tag+".cov");
  ofstream Print_Corr_FA_T("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FA_T_"+Analysis_tag+".corr");
  ofstream Print_Corr_FV_T("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FV_T_"+Analysis_tag+".corr");

  Print_Cov_FA_T<<Cov_FA_T<<endl;
  Print_Cov_FV_T<<Cov_FV_T<<endl;
  Print_Corr_FA_T<<Corr_FA_T<<endl;
  Print_Corr_FV_T<<Corr_FV_T<<endl;


  Print_Cov_FA_T.close();
  Print_Cov_FV_T.close();
  Print_Corr_FA_T.close();
  Print_Corr_FV_T.close();


  
      
	

  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################

                                                  //############ Up-quark CONTRIBUTION ##################



  
    
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //###########################################                                                                 ##########################################
  //###########################################                                                                 ##########################################
  //###########################################                   FIT FA_T(u) for all xg                        ##########################################
  //###########################################                                                                 ##########################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################


  
  
  Vfloat ch2_FA_T_u;
  distr_t_list F0_u_A_T_list(UseJack);
  distr_t_list D1_u_A_T_list(UseJack);
  distr_t_list D2_u_A_T_list(UseJack);
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
      data[ijack][iens].FF= FA_T_u_per_ens[iens].distr_list[ixg-1].distr[ijack];
      data[ijack][iens].FF_err= FA_T_u_per_ens[iens].err(ixg-1);
      if(data_2pts.Tag[iens] == "cA211a.12.48") { data[ijack][iens].a = a_A.distr[ijack]; data[ijack][iens].is=0; }
      else if(data_2pts.Tag[iens] == "cB211b.072.64") { data[ijack][iens].a = a_B.distr[ijack]; data[ijack][iens].is=1;}
      else if(data_2pts.Tag[iens] == "cB211b.072.96") { data[ijack][iens].a = a_B.distr[ijack]; data[ijack][iens].is=1; }
      else if(data_2pts.Tag[iens] == "cC211a.06.80")  {data[ijack][iens].a = a_C.distr[ijack]; data[ijack][iens].is=2; }
      else if(data_2pts.Tag[iens] == "cD211a.054.96")  {data[ijack][iens].a = a_D.distr[ijack]; data[ijack][iens].is=3; }
      else crash("Ens_tag: "+data_2pts.Tag[iens]+" not recognized");

      //mean values
      if(ijack==0) {
	data_ch2[ijack][iens].FF= FA_T_u_per_ens[iens].ave(ixg-1);
	data_ch2[ijack][iens].FF_err= FA_T_u_per_ens[iens].err(ixg-1);
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
  cout<<"Fitting FA_T(u), xg: "<<ixg<<endl;
  Bt_fit= bf_FF.Perform_bootstrap_fit();
  Bt_fit_ch2= bf_FF_ch2.Perform_bootstrap_fit();

  
  //retrieve parameters
  distr_t F0(UseJack), D1(UseJack), D2(UseJack);
  for(int ijack=0;ijack<Njacks;ijack++) { F0.distr.push_back( Bt_fit.par[ijack].F0); D1.distr.push_back( Bt_fit.par[ijack].D1); D2.distr.push_back( Bt_fit.par[ijack].D2);}
  //push_back retrieved parameters
  F0_u_A_T_list.distr_list.push_back(F0);
  D1_u_A_T_list.distr_list.push_back(D1);
  D2_u_A_T_list.distr_list.push_back(D2);
  //push_back ch2
  ch2_FA_T_u.push_back( Bt_fit_ch2.get_ch2_ave()/dof);

 
  //print fit func
  distr_t_list FA_T_xg_to_print(UseJack);
  for(auto &a: Bs_a_to_print) FA_T_xg_to_print.distr_list.push_back( F0 + D1*pow(a*fmTGeV*Lambda_QCD,2) + D2*pow(a*fmTGeV*Lambda_QCD,4));
  Print_To_File({}, {Bs_a_to_print, FA_T_xg_to_print.ave(), FA_T_xg_to_print.err()}, "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FA_T_u_"+Fit_tag+"xg_"+to_string_with_precision(0.10*ixg,2)+".fit_func", "", "#a[fm] FA FA_err");
 
  }


  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //###########################################                                                                 ##########################################
  //###########################################                                                                 ##########################################
  //###########################################                   FIT FV_T(u) for all xg                        ##########################################
  //###########################################                                                                 ##########################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  Vfloat ch2_FV_T_u;
  bf_FF.Set_par_val("F0", 0.1, 0.001);
  bf_FF.Set_par_val("D1", 1.0, 0.1);
  if(Include_a4) bf_FF.Set_par_val("D2", 1.0, 0.1);
  bf_FF_ch2.Set_par_val("F0", 0.1, 0.001);
  bf_FF_ch2.Set_par_val("D1", 1.0, 0.1);
  if(Include_a4) bf_FF_ch2.Set_par_val("D2", 1.0, 0.1);
  distr_t_list F0_u_V_T_list(UseJack);
  distr_t_list D1_u_V_T_list(UseJack);
  distr_t_list D2_u_V_T_list(UseJack);
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
	data[ijack][iens].FF= FV_T_u_per_ens[iens].distr_list[ixg-1].distr[ijack];
	data[ijack][iens].FF_err= FV_T_u_per_ens[iens].err(ixg-1);
	if(data_2pts.Tag[iens] == "cA211a.12.48") { data[ijack][iens].a = a_A.distr[ijack]; data[ijack][iens].is=0; }
	else if(data_2pts.Tag[iens] == "cB211b.072.64") { data[ijack][iens].a = a_B.distr[ijack]; data[ijack][iens].is=1;}
	else if(data_2pts.Tag[iens] == "cB211b.072.96") { data[ijack][iens].a = a_B.distr[ijack]; data[ijack][iens].is=1; }
	else if(data_2pts.Tag[iens] == "cC211a.06.80")  {data[ijack][iens].a = a_C.distr[ijack]; data[ijack][iens].is=2; }
	else if(data_2pts.Tag[iens] == "cD211a.054.96")  {data[ijack][iens].a = a_D.distr[ijack]; data[ijack][iens].is=3; }
	else crash("Ens_tag: "+data_2pts.Tag[iens]+" not recognized");

	if(ijack==0) {
	  	data_ch2[ijack][iens].FF= FV_T_u_per_ens[iens].ave(ixg-1);
		data_ch2[ijack][iens].FF_err= FV_T_u_per_ens[iens].err(ixg-1);
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
  cout<<"Fitting FV_T(u), xg: "<<ixg<<endl;
  Bt_fit= bf_FF.Perform_bootstrap_fit();
  Bt_fit_ch2= bf_FF_ch2.Perform_bootstrap_fit();
  //retrieve parameters
  distr_t F0(UseJack), D1(UseJack), D2(UseJack);
  for(int ijack=0;ijack<Njacks;ijack++) { F0.distr.push_back( Bt_fit.par[ijack].F0); D1.distr.push_back( Bt_fit.par[ijack].D1); D2.distr.push_back( Bt_fit.par[ijack].D2);}
  //push_back retrieved parameters
  F0_u_V_T_list.distr_list.push_back(F0);
  D1_u_V_T_list.distr_list.push_back(D1);
  D2_u_V_T_list.distr_list.push_back(D2);
  //push_back ch2
  ch2_FV_T_u.push_back( Bt_fit_ch2.get_ch2_ave()/dof);

  //print fit func
  distr_t_list FV_T_xg_to_print(UseJack);
  for(auto &a: Bs_a_to_print) FV_T_xg_to_print.distr_list.push_back( F0 + D1*pow(a*fmTGeV*Lambda_QCD,2) + D2*pow(a*fmTGeV*Lambda_QCD,4));
  Print_To_File({}, {Bs_a_to_print, FV_T_xg_to_print.ave(), FV_T_xg_to_print.err()}, "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FV_T_u_"+Fit_tag+"xg_"+to_string_with_precision(0.10*ixg,2)+".fit_func", "", "#a[fm] FV FV_err");
  }

  //Print continuum extrapolated form factors
  Print_To_File({}, {Bs_xg_t_list, F0_u_A_T_list.ave(), F0_u_A_T_list.err(), (D1_u_A_T_list/F0_u_A_T_list).ave(), (D1_u_A_T_list/F0_u_A_T_list).err(), (D2_u_A_T_list/F0_u_A_T_list).ave(), (D2_u_A_T_list/F0_u_A_T_list).err(), ch2_FA_T_u}, "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FA_T_u_"+Fit_tag+"cont.dat", "", "#xg  F0   D1   D2   ch2/dof");
  Print_To_File({}, {Bs_xg_t_list, F0_u_V_T_list.ave(), F0_u_V_T_list.err(), (D1_u_V_T_list/F0_u_V_T_list).ave(), (D1_u_V_T_list/F0_u_V_T_list).err(), (D2_u_V_T_list/F0_u_V_T_list).ave(), (D2_u_V_T_list/F0_u_V_T_list).err(), ch2_FV_T_u}, "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FV_T_u_"+Fit_tag+"cont.dat", "", "#xg  F0   D1   D2   ch2/dof");
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////



  Eigen::MatrixXd Cov_FA_T_u(num_xg-1, num_xg-1);
  Eigen::MatrixXd Cov_FV_T_u(num_xg-1, num_xg-1);
  Eigen::MatrixXd Corr_FA_T_u(num_xg-1, num_xg-1);
  Eigen::MatrixXd Corr_FV_T_u(num_xg-1, num_xg-1);

  for(int x_xg=1; x_xg<num_xg;x_xg++) {
    for(int y_xg=1; y_xg<num_xg;y_xg++) {
      Cov_FA_T_u(x_xg-1, y_xg-1) = F0_u_A_T_list.distr_list[x_xg-1]%F0_u_A_T_list.distr_list[y_xg-1];
      Cov_FV_T_u(x_xg-1, y_xg-1) = F0_u_V_T_list.distr_list[x_xg-1]%F0_u_V_T_list.distr_list[y_xg-1];
      Corr_FA_T_u(x_xg-1, y_xg-1) = Cov_FA_T_u(x_xg-1,y_xg-1)/(F0_u_A_T_list.err(x_xg-1)*F0_u_A_T_list.err(y_xg-1));
      Corr_FV_T_u(x_xg-1, y_xg-1) = Cov_FV_T_u(x_xg-1,y_xg-1)/(F0_u_V_T_list.err(x_xg-1)*F0_u_V_T_list.err(y_xg-1));
    }
  }

  //Print To File
  ofstream Print_Cov_FA_T_u("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FA_T_u_"+Analysis_tag+".cov");
  ofstream Print_Cov_FV_T_u("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FV_T_u_"+Analysis_tag+".cov");
  ofstream Print_Corr_FA_T_u("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FA_T_u_"+Analysis_tag+".corr");
  ofstream Print_Corr_FV_T_u("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FV_T_u_"+Analysis_tag+".corr");

  Print_Cov_FA_T_u<<Cov_FA_T_u<<endl;
  Print_Cov_FV_T_u<<Cov_FV_T_u<<endl;
  Print_Corr_FA_T_u<<Corr_FA_T_u<<endl;
  Print_Corr_FV_T_u<<Corr_FV_T_u<<endl;


  Print_Cov_FA_T_u.close();
  Print_Cov_FV_T_u.close();
  Print_Corr_FA_T_u.close();
  Print_Corr_FV_T_u.close();





  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################

                                                  //############ down-quark CONTRIBUTION ##################



  
    
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //###########################################                                                                 ##########################################
  //###########################################                                                                 ##########################################
  //###########################################                   FIT FA_T(d) for all xg                        ##########################################
  //###########################################                                                                 ##########################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
 
  Vfloat ch2_FA_T_d;
  distr_t_list F0_d_A_T_list(UseJack);
  distr_t_list D1_d_A_T_list(UseJack);
  distr_t_list D2_d_A_T_list(UseJack);
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
      data[ijack][iens].FF= FA_T_d_per_ens[iens].distr_list[ixg-1].distr[ijack];
      data[ijack][iens].FF_err= FA_T_d_per_ens[iens].err(ixg-1);
      if(data_2pts.Tag[iens] == "cA211a.12.48") { data[ijack][iens].a = a_A.distr[ijack]; data[ijack][iens].is=0; }
      else if(data_2pts.Tag[iens] == "cB211b.072.64") { data[ijack][iens].a = a_B.distr[ijack]; data[ijack][iens].is=1;}
      else if(data_2pts.Tag[iens] == "cB211b.072.96") { data[ijack][iens].a = a_B.distr[ijack]; data[ijack][iens].is=1; }
      else if(data_2pts.Tag[iens] == "cC211a.06.80")  {data[ijack][iens].a = a_C.distr[ijack]; data[ijack][iens].is=2; }
      else if(data_2pts.Tag[iens] == "cD211a.054.96")  {data[ijack][iens].a = a_D.distr[ijack]; data[ijack][iens].is=3; }
      else crash("Ens_tag: "+data_2pts.Tag[iens]+" not recognized");

      //mean values
      if(ijack==0) {
	data_ch2[ijack][iens].FF= FA_T_d_per_ens[iens].ave(ixg-1);
	data_ch2[ijack][iens].FF_err= FA_T_d_per_ens[iens].err(ixg-1);
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
  cout<<"Fitting FA_T(u), xg: "<<ixg<<endl;
  Bt_fit= bf_FF.Perform_bootstrap_fit();
  Bt_fit_ch2= bf_FF_ch2.Perform_bootstrap_fit();

  
  //retrieve parameters
  distr_t F0(UseJack), D1(UseJack), D2(UseJack);
  for(int ijack=0;ijack<Njacks;ijack++) { F0.distr.push_back( Bt_fit.par[ijack].F0); D1.distr.push_back( Bt_fit.par[ijack].D1); D2.distr.push_back( Bt_fit.par[ijack].D2);}
  //push_back retrieved parameters
  F0_d_A_T_list.distr_list.push_back(F0);
  D1_d_A_T_list.distr_list.push_back(D1);
  D2_d_A_T_list.distr_list.push_back(D2);
  //push_back ch2
  ch2_FA_T_d.push_back( Bt_fit_ch2.get_ch2_ave()/dof);

 
  //print fit func
  distr_t_list FA_T_xg_to_print(UseJack);
  for(auto &a: Bs_a_to_print) FA_T_xg_to_print.distr_list.push_back( F0 + D1*pow(a*fmTGeV*Lambda_QCD,2) + D2*pow(a*fmTGeV*Lambda_QCD,4));
  Print_To_File({}, {Bs_a_to_print, FA_T_xg_to_print.ave(), FA_T_xg_to_print.err()}, "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FA_T_d_"+Fit_tag+"xg_"+to_string_with_precision(0.10*ixg,2)+".fit_func", "", "#a[fm] FA FA_err");
 
  }


  
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //###########################################                                                                 ##########################################
  //###########################################                                                                 ##########################################
  //###########################################                   FIT FV_T(d) for all xg                        ##########################################
  //###########################################                                                                 ##########################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  Vfloat ch2_FV_T_d;
  bf_FF.Set_par_val("F0", -0.1, 0.001);
  bf_FF.Set_par_val("D1", 1.0, 0.1);
  if(Include_a4) bf_FF.Set_par_val("D2", 1.0, 0.1);
  bf_FF_ch2.Set_par_val("F0", -0.1, 0.001);
  bf_FF_ch2.Set_par_val("D1", 1.0, 0.1);
  if(Include_a4) bf_FF_ch2.Set_par_val("D2", 1.0, 0.1);
  distr_t_list F0_d_V_T_list(UseJack);
  distr_t_list D1_d_V_T_list(UseJack);
  distr_t_list D2_d_V_T_list(UseJack);
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
	data[ijack][iens].FF= FV_T_d_per_ens[iens].distr_list[ixg-1].distr[ijack];
	data[ijack][iens].FF_err= FV_T_d_per_ens[iens].err(ixg-1);
	if(data_2pts.Tag[iens] == "cA211a.12.48") { data[ijack][iens].a = a_A.distr[ijack]; data[ijack][iens].is=0; }
	else if(data_2pts.Tag[iens] == "cB211b.072.64") { data[ijack][iens].a = a_B.distr[ijack]; data[ijack][iens].is=1;}
	else if(data_2pts.Tag[iens] == "cB211b.072.96") { data[ijack][iens].a = a_B.distr[ijack]; data[ijack][iens].is=1; }
	else if(data_2pts.Tag[iens] == "cC211a.06.80")  {data[ijack][iens].a = a_C.distr[ijack]; data[ijack][iens].is=2; }
	else if(data_2pts.Tag[iens] == "cD211a.054.96")  {data[ijack][iens].a = a_D.distr[ijack]; data[ijack][iens].is=3; }
	else crash("Ens_tag: "+data_2pts.Tag[iens]+" not recognized");

	if(ijack==0) {
	  	data_ch2[ijack][iens].FF= FV_T_d_per_ens[iens].ave(ixg-1);
		data_ch2[ijack][iens].FF_err= FV_T_d_per_ens[iens].err(ixg-1);
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
  cout<<"Fitting FV_T(d), xg: "<<ixg<<endl;
  Bt_fit= bf_FF.Perform_bootstrap_fit();
  Bt_fit_ch2= bf_FF_ch2.Perform_bootstrap_fit();
  //retrieve parameters
  distr_t F0(UseJack), D1(UseJack), D2(UseJack);
  for(int ijack=0;ijack<Njacks;ijack++) { F0.distr.push_back( Bt_fit.par[ijack].F0); D1.distr.push_back( Bt_fit.par[ijack].D1); D2.distr.push_back( Bt_fit.par[ijack].D2);}
  //push_back retrieved parameters
  F0_d_V_T_list.distr_list.push_back(F0);
  D1_d_V_T_list.distr_list.push_back(D1);
  D2_d_V_T_list.distr_list.push_back(D2);
  //push_back ch2
  ch2_FV_T_d.push_back( Bt_fit_ch2.get_ch2_ave()/dof);

  //print fit func
  distr_t_list FV_T_xg_to_print(UseJack);
  for(auto &a: Bs_a_to_print) FV_T_xg_to_print.distr_list.push_back( F0 + D1*pow(a*fmTGeV*Lambda_QCD,2) + D2*pow(a*fmTGeV*Lambda_QCD,4));
  Print_To_File({}, {Bs_a_to_print, FV_T_xg_to_print.ave(), FV_T_xg_to_print.err()}, "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FV_T_d_"+Fit_tag+"xg_"+to_string_with_precision(0.10*ixg,2)+".fit_func", "", "#a[fm] FV FV_err");
  }

  //Print continuum extrapolated form factors
  Print_To_File({}, {Bs_xg_t_list, F0_d_A_T_list.ave(), F0_d_A_T_list.err(), (D1_d_A_T_list/F0_d_A_T_list).ave(), (D1_d_A_T_list/F0_d_A_T_list).err(), (D2_d_A_T_list/F0_d_A_T_list).ave(), (D2_d_A_T_list/F0_d_A_T_list).err(), ch2_FA_T_d}, "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FA_T_d_"+Fit_tag+"cont.dat", "", "#xg  F0   D1   D2   ch2/dof");
  Print_To_File({}, {Bs_xg_t_list, F0_d_V_T_list.ave(), F0_d_V_T_list.err(), (D1_d_V_T_list/F0_d_V_T_list).ave(), (D1_d_V_T_list/F0_d_V_T_list).err(), (D2_d_V_T_list/F0_d_V_T_list).ave(), (D2_d_V_T_list/F0_d_V_T_list).err(), ch2_FV_T_d}, "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FV_T_d_"+Fit_tag+"cont.dat", "", "#xg  F0   D1   D2   ch2/dof");
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////


  

  Eigen::MatrixXd Cov_FA_T_d(num_xg-1, num_xg-1);
  Eigen::MatrixXd Cov_FV_T_d(num_xg-1, num_xg-1);
  Eigen::MatrixXd Corr_FA_T_d(num_xg-1, num_xg-1);
  Eigen::MatrixXd Corr_FV_T_d(num_xg-1, num_xg-1);

  for(int x_xg=1; x_xg<num_xg;x_xg++) {
    for(int y_xg=1; y_xg<num_xg;y_xg++) {
      Cov_FA_T_d(x_xg-1, y_xg-1) = F0_d_A_T_list.distr_list[x_xg-1]%F0_d_A_T_list.distr_list[y_xg-1];
      Cov_FV_T_d(x_xg-1, y_xg-1) = F0_d_V_T_list.distr_list[x_xg-1]%F0_d_V_T_list.distr_list[y_xg-1];
      Corr_FA_T_d(x_xg-1, y_xg-1) = Cov_FA_T_d(x_xg-1,y_xg-1)/(F0_d_A_T_list.err(x_xg-1)*F0_d_A_T_list.err(y_xg-1));
      Corr_FV_T_d(x_xg-1, y_xg-1) = Cov_FV_T_d(x_xg-1,y_xg-1)/(F0_d_V_T_list.err(x_xg-1)*F0_d_V_T_list.err(y_xg-1));
    }
  }

  //Print To File
  ofstream Print_Cov_FA_T_d("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FA_T_d_"+Analysis_tag+".cov");
  ofstream Print_Cov_FV_T_d("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FV_T_d_"+Analysis_tag+".cov");
  ofstream Print_Corr_FA_T_d("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FA_T_d_"+Analysis_tag+".corr");
  ofstream Print_Corr_FV_T_d("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FV_T_d_"+Analysis_tag+".corr");

  Print_Cov_FA_T_d<<Cov_FA_T_d<<endl;
  Print_Cov_FV_T_d<<Cov_FV_T_d<<endl;
  Print_Corr_FA_T_d<<Corr_FA_T_d<<endl;
  Print_Corr_FV_T_d<<Corr_FV_T_d<<endl;


  Print_Cov_FA_T_d.close();
  Print_Cov_FV_T_d.close();
  Print_Corr_FA_T_d.close();
  Print_Corr_FV_T_d.close();




  //###################################################################################################################################################################

  



  //here
  //###################################################################################################################################################################


  

  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //###########################################                                                                 ##########################################
  //###########################################                                                                 ##########################################
  //###########################################                                                                 ##########################################
  //###########################################                                                                 ##########################################
  //###########################################                                                                 ##########################################
  //###########################################                                                                 ##########################################
  //###########################################                           FITTING FT AND FB                     ##########################################
  //###########################################                     INCLUDING SINGLE CONTRIBUTIONS              ##########################################
  //###########################################                        FROM STRANGE AND "CHARM"                 ##########################################
  //###########################################                                                                 ##########################################
  //###########################################                                                                 ##########################################
  //###########################################                                                                 ##########################################
  //###########################################                                                                 ##########################################
  //###########################################                                                                 ##########################################
  //###########################################                                                                 ##########################################
  //###########################################                                                                 ##########################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################




  
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //###########################################                                                                 ##########################################
  //###########################################                                                                 ##########################################
  //###########################################                   FIT FB for all xg                           ##########################################
  //###########################################                                                                 ##########################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################

  Vfloat ch2_FB;
  distr_t_list F0_B_list(UseJack);
  distr_t_list D1_B_list(UseJack);
  distr_t_list D2_B_list(UseJack);
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
      data[ijack][iens].FF= FB_per_ens[iens].distr_list[ixg-1].distr[ijack];
      data[ijack][iens].FF_err= FB_per_ens[iens].err(ixg-1);
      if(data_2pts.Tag[iens] == "cA211a.12.48") { data[ijack][iens].a = a_A.distr[ijack]; data[ijack][iens].is=0; }
      else if(data_2pts.Tag[iens] == "cB211b.072.64") { data[ijack][iens].a = a_B.distr[ijack]; data[ijack][iens].is=1;}
      else if(data_2pts.Tag[iens] == "cB211b.072.96") { data[ijack][iens].a = a_B.distr[ijack]; data[ijack][iens].is=1; }
      else if(data_2pts.Tag[iens] == "cC211a.06.80")  {data[ijack][iens].a = a_C.distr[ijack]; data[ijack][iens].is=2; }
      else if(data_2pts.Tag[iens] == "cD211a.054.96")  {data[ijack][iens].a = a_D.distr[ijack]; data[ijack][iens].is=3; }
      else crash("Ens_tag: "+data_2pts.Tag[iens]+" not recognized");

      //mean values
      if(ijack==0) {
	data_ch2[ijack][iens].FF= FB_per_ens[iens].ave(ixg-1);
	data_ch2[ijack][iens].FF_err= FB_per_ens[iens].err(ixg-1);
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
  cout<<"Fitting FB, xg: "<<ixg<<endl;
  Bt_fit= bf_FF.Perform_bootstrap_fit();
  Bt_fit_ch2= bf_FF_ch2.Perform_bootstrap_fit();

  
  //retrieve parameters
  distr_t F0(UseJack), D1(UseJack), D2(UseJack);
  for(int ijack=0;ijack<Njacks;ijack++) { F0.distr.push_back( Bt_fit.par[ijack].F0); D1.distr.push_back( Bt_fit.par[ijack].D1); D2.distr.push_back( Bt_fit.par[ijack].D2);}
  //push_back retrieved parameters
  F0_B_list.distr_list.push_back(F0);
  D1_B_list.distr_list.push_back(D1);
  D2_B_list.distr_list.push_back(D2);
  //push_back ch2
  ch2_FB.push_back( Bt_fit_ch2.get_ch2_ave()/dof);

 
  //print fit func
  distr_t_list FB_xg_to_print(UseJack);
  for(auto &a: Bs_a_to_print) FB_xg_to_print.distr_list.push_back( F0 + D1*pow(a*fmTGeV*Lambda_QCD,2) + D2*pow(a*fmTGeV*Lambda_QCD,4));
  Print_To_File({}, {Bs_a_to_print, FB_xg_to_print.ave(), FB_xg_to_print.err()}, "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FB_"+Fit_tag+"xg_"+to_string_with_precision(0.10*ixg,2)+".fit_func", "", "#a[fm] FA FA_err");
 
  }


  
    
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //###########################################                                                                 ##########################################
  //###########################################                                                                 ##########################################
  //###########################################                   FIT FT for all xg                           ##########################################
  //###########################################                                                                 ##########################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################

  
  Vfloat ch2_FT;
  bf_FF.Set_par_val("F0", 0.1, 0.001);
  bf_FF.Set_par_val("D1", 1.0, 0.1);
  if(Include_a4) bf_FF.Set_par_val("D2", 1.0, 0.1);
  bf_FF_ch2.Set_par_val("F0", 0.1, 0.001);
  bf_FF_ch2.Set_par_val("D1", 1.0, 0.1);
  if(Include_a4) bf_FF_ch2.Set_par_val("D2", 1.0, 0.1);
  distr_t_list F0_T_list(UseJack);
  distr_t_list D1_T_list(UseJack);
  distr_t_list D2_T_list(UseJack);
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
	data[ijack][iens].FF= FT_per_ens[iens].distr_list[ixg-1].distr[ijack];
	data[ijack][iens].FF_err= FT_per_ens[iens].err(ixg-1);
	if(data_2pts.Tag[iens] == "cA211a.12.48") { data[ijack][iens].a = a_A.distr[ijack]; data[ijack][iens].is=0; }
	else if(data_2pts.Tag[iens] == "cB211b.072.64") { data[ijack][iens].a = a_B.distr[ijack]; data[ijack][iens].is=1;}
	else if(data_2pts.Tag[iens] == "cB211b.072.96") { data[ijack][iens].a = a_B.distr[ijack]; data[ijack][iens].is=1; }
	else if(data_2pts.Tag[iens] == "cC211a.06.80")  {data[ijack][iens].a = a_C.distr[ijack]; data[ijack][iens].is=2; }
	else if(data_2pts.Tag[iens] == "cD211a.054.96")  {data[ijack][iens].a = a_D.distr[ijack]; data[ijack][iens].is=3; }
	else crash("Ens_tag: "+data_2pts.Tag[iens]+" not recognized");

	if(ijack==0) {
	  	data_ch2[ijack][iens].FF= FT_per_ens[iens].ave(ixg-1);
		data_ch2[ijack][iens].FF_err= FT_per_ens[iens].err(ixg-1);
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
  cout<<"Fitting FT, xg: "<<ixg<<endl;
  Bt_fit= bf_FF.Perform_bootstrap_fit();
  Bt_fit_ch2= bf_FF_ch2.Perform_bootstrap_fit();
  //retrieve parameters
  distr_t F0(UseJack), D1(UseJack), D2(UseJack);
  for(int ijack=0;ijack<Njacks;ijack++) { F0.distr.push_back( Bt_fit.par[ijack].F0); D1.distr.push_back( Bt_fit.par[ijack].D1); D2.distr.push_back( Bt_fit.par[ijack].D2);}
  //push_back retrieved parameters
  F0_T_list.distr_list.push_back(F0);
  D1_T_list.distr_list.push_back(D1);
  D2_T_list.distr_list.push_back(D2);
  //push_back ch2
  ch2_FT.push_back( Bt_fit_ch2.get_ch2_ave()/dof);

  //print fit func
  distr_t_list FT_xg_to_print(UseJack);
  for(auto &a: Bs_a_to_print) FT_xg_to_print.distr_list.push_back( F0 + D1*pow(a*fmTGeV*Lambda_QCD,2) + D2*pow(a*fmTGeV*Lambda_QCD,4));
  Print_To_File({}, {Bs_a_to_print, FT_xg_to_print.ave(), FT_xg_to_print.err()}, "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FT_"+Fit_tag+"xg_"+to_string_with_precision(0.10*ixg,2)+".fit_func", "", "#a[fm] FV FV_err");
  }

 

  //Print continuum extrapolated form factors
  Print_To_File({}, {Bs_xg_t_list, F0_B_list.ave(), F0_B_list.err(), (D1_B_list/F0_B_list).ave(), (D1_B_list/F0_B_list).err(), (D2_B_list/F0_B_list).ave(), (D2_B_list/F0_B_list).err(), ch2_FB}, "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FB_"+Fit_tag+"cont.dat", "", "#xg  F0   D1   D2  ch2/dof");
  Print_To_File({}, {Bs_xg_t_list, F0_T_list.ave(), F0_T_list.err(), (D1_T_list/F0_T_list).ave(), (D1_T_list/F0_T_list).err(), (D2_T_list/F0_T_list).ave(), (D2_T_list/F0_T_list).err(), ch2_FT}, "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FT_"+Fit_tag+"cont.dat", "", "#xg  F0   D1   D2  ch2/dof");


  //Print covariance matrix
  Analysis_tag= ( (UseJack==true)?"jack":"boot")+_Fit_tag;
  Eigen::MatrixXd Cov_FB(num_xg-1, num_xg-1);
  Eigen::MatrixXd Cov_FT(num_xg-1, num_xg-1);
  Eigen::MatrixXd Corr_FB(num_xg-1, num_xg-1);
  Eigen::MatrixXd Corr_FT(num_xg-1, num_xg-1);

  for(int x_xg=1; x_xg<num_xg;x_xg++) {
    for(int y_xg=1; y_xg<num_xg;y_xg++) {
      Cov_FB(x_xg-1, y_xg-1) = F0_B_list.distr_list[x_xg-1]%F0_B_list.distr_list[y_xg-1];
      Cov_FT(x_xg-1, y_xg-1) = F0_T_list.distr_list[x_xg-1]%F0_T_list.distr_list[y_xg-1];
      Corr_FB(x_xg-1, y_xg-1) = Cov_FB(x_xg-1,y_xg-1)/(F0_B_list.err(x_xg-1)*F0_B_list.err(y_xg-1));
      Corr_FT(x_xg-1, y_xg-1) = Cov_FT(x_xg-1,y_xg-1)/(F0_T_list.err(x_xg-1)*F0_T_list.err(y_xg-1));
    }
  }

  //Print To File
  ofstream Print_Cov_FB("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FB_"+Analysis_tag+".cov");
  ofstream Print_Cov_FT("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FT_"+Analysis_tag+".cov");
  ofstream Print_Corr_FB("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FB_"+Analysis_tag+".corr");
  ofstream Print_Corr_FT("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FT_"+Analysis_tag+".corr");

  Print_Cov_FB<<Cov_FB<<endl;
  Print_Cov_FT<<Cov_FT<<endl;
  Print_Corr_FB<<Corr_FB<<endl;
  Print_Corr_FT<<Corr_FT<<endl;


  Print_Cov_FB.close();
  Print_Cov_FT.close();
  Print_Corr_FB.close();
  Print_Corr_FT.close();


  
      
	

  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################

                                                  //############ Up-quark CONTRIBUTION ##################



  
    
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //###########################################                                                                 ##########################################
  //###########################################                                                                 ##########################################
  //###########################################                   FIT FB(u) for all xg                        ##########################################
  //###########################################                                                                 ##########################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################


  
  
  Vfloat ch2_FB_u;
  distr_t_list F0_u_B_list(UseJack);
  distr_t_list D1_u_B_list(UseJack);
  distr_t_list D2_u_B_list(UseJack);
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
      data[ijack][iens].FF= FB_u_per_ens[iens].distr_list[ixg-1].distr[ijack];
      data[ijack][iens].FF_err= FB_u_per_ens[iens].err(ixg-1);
      if(data_2pts.Tag[iens] == "cA211a.12.48") { data[ijack][iens].a = a_A.distr[ijack]; data[ijack][iens].is=0; }
      else if(data_2pts.Tag[iens] == "cB211b.072.64") { data[ijack][iens].a = a_B.distr[ijack]; data[ijack][iens].is=1;}
      else if(data_2pts.Tag[iens] == "cB211b.072.96") { data[ijack][iens].a = a_B.distr[ijack]; data[ijack][iens].is=1; }
      else if(data_2pts.Tag[iens] == "cC211a.06.80")  {data[ijack][iens].a = a_C.distr[ijack]; data[ijack][iens].is=2; }
      else if(data_2pts.Tag[iens] == "cD211a.054.96")  {data[ijack][iens].a = a_D.distr[ijack]; data[ijack][iens].is=3; }
      else crash("Ens_tag: "+data_2pts.Tag[iens]+" not recognized");

      //mean values
      if(ijack==0) {
	data_ch2[ijack][iens].FF= FB_u_per_ens[iens].ave(ixg-1);
	data_ch2[ijack][iens].FF_err= FB_u_per_ens[iens].err(ixg-1);
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
  cout<<"Fitting FB(u), xg: "<<ixg<<endl;
  Bt_fit= bf_FF.Perform_bootstrap_fit();
  Bt_fit_ch2= bf_FF_ch2.Perform_bootstrap_fit();

  
  //retrieve parameters
  distr_t F0(UseJack), D1(UseJack), D2(UseJack);
  for(int ijack=0;ijack<Njacks;ijack++) { F0.distr.push_back( Bt_fit.par[ijack].F0); D1.distr.push_back( Bt_fit.par[ijack].D1); D2.distr.push_back( Bt_fit.par[ijack].D2);}
  //push_back retrieved parameters
  F0_u_B_list.distr_list.push_back(F0);
  D1_u_B_list.distr_list.push_back(D1);
  D2_u_B_list.distr_list.push_back(D2);
  //push_back ch2
  ch2_FB_u.push_back( Bt_fit_ch2.get_ch2_ave()/dof);

 
  //print fit func
  distr_t_list FB_xg_to_print(UseJack);
  for(auto &a: Bs_a_to_print) FB_xg_to_print.distr_list.push_back( F0 + D1*pow(a*fmTGeV*Lambda_QCD,2) + D2*pow(a*fmTGeV*Lambda_QCD,4));
  Print_To_File({}, {Bs_a_to_print, FB_xg_to_print.ave(), FB_xg_to_print.err()}, "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FB_u_"+Fit_tag+"xg_"+to_string_with_precision(0.10*ixg,2)+".fit_func", "", "#a[fm] FA FA_err");
 
  }


  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //###########################################                                                                 ##########################################
  //###########################################                                                                 ##########################################
  //###########################################                   FIT FT(u) for all xg                        ##########################################
  //###########################################                                                                 ##########################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  Vfloat ch2_FT_u;
  bf_FF.Set_par_val("F0", 0.1, 0.001);
  bf_FF.Set_par_val("D1", 1.0, 0.1);
  if(Include_a4) bf_FF.Set_par_val("D2", 1.0, 0.1);
  bf_FF_ch2.Set_par_val("F0", 0.1, 0.001);
  bf_FF_ch2.Set_par_val("D1", 1.0, 0.1);
  if(Include_a4) bf_FF_ch2.Set_par_val("D2", 1.0, 0.1);
  distr_t_list F0_u_T_list(UseJack);
  distr_t_list D1_u_T_list(UseJack);
  distr_t_list D2_u_T_list(UseJack);
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
	data[ijack][iens].FF= FT_u_per_ens[iens].distr_list[ixg-1].distr[ijack];
	data[ijack][iens].FF_err= FT_u_per_ens[iens].err(ixg-1);
	if(data_2pts.Tag[iens] == "cA211a.12.48") { data[ijack][iens].a = a_A.distr[ijack]; data[ijack][iens].is=0; }
	else if(data_2pts.Tag[iens] == "cB211b.072.64") { data[ijack][iens].a = a_B.distr[ijack]; data[ijack][iens].is=1;}
	else if(data_2pts.Tag[iens] == "cB211b.072.96") { data[ijack][iens].a = a_B.distr[ijack]; data[ijack][iens].is=1; }
	else if(data_2pts.Tag[iens] == "cC211a.06.80")  {data[ijack][iens].a = a_C.distr[ijack]; data[ijack][iens].is=2; }
	else if(data_2pts.Tag[iens] == "cD211a.054.96")  {data[ijack][iens].a = a_D.distr[ijack]; data[ijack][iens].is=3; }
	else crash("Ens_tag: "+data_2pts.Tag[iens]+" not recognized");

	if(ijack==0) {
	  	data_ch2[ijack][iens].FF= FT_u_per_ens[iens].ave(ixg-1);
		data_ch2[ijack][iens].FF_err= FT_u_per_ens[iens].err(ixg-1);
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
  cout<<"Fitting FT(u), xg: "<<ixg<<endl;
  Bt_fit= bf_FF.Perform_bootstrap_fit();
  Bt_fit_ch2= bf_FF_ch2.Perform_bootstrap_fit();
  //retrieve parameters
  distr_t F0(UseJack), D1(UseJack), D2(UseJack);
  for(int ijack=0;ijack<Njacks;ijack++) { F0.distr.push_back( Bt_fit.par[ijack].F0); D1.distr.push_back( Bt_fit.par[ijack].D1); D2.distr.push_back( Bt_fit.par[ijack].D2);}
  //push_back retrieved parameters
  F0_u_T_list.distr_list.push_back(F0);
  D1_u_T_list.distr_list.push_back(D1);
  D2_u_T_list.distr_list.push_back(D2);
  //push_back ch2
  ch2_FT_u.push_back( Bt_fit_ch2.get_ch2_ave()/dof);

  //print fit func
  distr_t_list FT_xg_to_print(UseJack);
  for(auto &a: Bs_a_to_print) FT_xg_to_print.distr_list.push_back( F0 + D1*pow(a*fmTGeV*Lambda_QCD,2) + D2*pow(a*fmTGeV*Lambda_QCD,4));
  Print_To_File({}, {Bs_a_to_print, FT_xg_to_print.ave(), FT_xg_to_print.err()}, "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FT_u_"+Fit_tag+"xg_"+to_string_with_precision(0.10*ixg,2)+".fit_func", "", "#a[fm] FV FV_err");
  }

  //Print continuum extrapolated form factors
  Print_To_File({}, {Bs_xg_t_list, F0_u_B_list.ave(), F0_u_B_list.err(), (D1_u_B_list/F0_u_B_list).ave(), (D1_u_B_list/F0_u_B_list).err(), (D2_u_B_list/F0_u_B_list).ave(), (D2_u_B_list/F0_u_B_list).err(), ch2_FB_u}, "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FB_u_"+Fit_tag+"cont.dat", "", "#xg  F0   D1   D2   ch2/dof");
  Print_To_File({}, {Bs_xg_t_list, F0_u_T_list.ave(), F0_u_T_list.err(), (D1_u_T_list/F0_u_T_list).ave(), (D1_u_T_list/F0_u_T_list).err(), (D2_u_T_list/F0_u_T_list).ave(), (D2_u_T_list/F0_u_T_list).err(), ch2_FT_u}, "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FT_u_"+Fit_tag+"cont.dat", "", "#xg  F0   D1   D2   ch2/dof");
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////



  Eigen::MatrixXd Cov_FB_u(num_xg-1, num_xg-1);
  Eigen::MatrixXd Cov_FT_u(num_xg-1, num_xg-1);
  Eigen::MatrixXd Corr_FB_u(num_xg-1, num_xg-1);
  Eigen::MatrixXd Corr_FT_u(num_xg-1, num_xg-1);

  for(int x_xg=1; x_xg<num_xg;x_xg++) {
    for(int y_xg=1; y_xg<num_xg;y_xg++) {
      Cov_FB_u(x_xg-1, y_xg-1) = F0_u_B_list.distr_list[x_xg-1]%F0_u_B_list.distr_list[y_xg-1];
      Cov_FT_u(x_xg-1, y_xg-1) = F0_u_T_list.distr_list[x_xg-1]%F0_u_T_list.distr_list[y_xg-1];
      Corr_FB_u(x_xg-1, y_xg-1) = Cov_FB_u(x_xg-1,y_xg-1)/(F0_u_B_list.err(x_xg-1)*F0_u_B_list.err(y_xg-1));
      Corr_FT_u(x_xg-1, y_xg-1) = Cov_FT_u(x_xg-1,y_xg-1)/(F0_u_T_list.err(x_xg-1)*F0_u_T_list.err(y_xg-1));
    }
  }

  //Print To File
  ofstream Print_Cov_FB_u("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FB_u_"+Analysis_tag+".cov");
  ofstream Print_Cov_FT_u("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FT_u_"+Analysis_tag+".cov");
  ofstream Print_Corr_FB_u("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FB_u_"+Analysis_tag+".corr");
  ofstream Print_Corr_FT_u("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FT_u_"+Analysis_tag+".corr");

  Print_Cov_FB_u<<Cov_FB_u<<endl;
  Print_Cov_FT_u<<Cov_FT_u<<endl;
  Print_Corr_FB_u<<Corr_FB_u<<endl;
  Print_Corr_FT_u<<Corr_FT_u<<endl;


  Print_Cov_FB_u.close();
  Print_Cov_FT_u.close();
  Print_Corr_FB_u.close();
  Print_Corr_FT_u.close();





  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################

                                                  //############ down-quark CONTRIBUTION ##################



  
    
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //###########################################                                                                 ##########################################
  //###########################################                                                                 ##########################################
  //###########################################                   FIT FB(d) for all xg                        ##########################################
  //###########################################                                                                 ##########################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
 
  Vfloat ch2_FB_d;
  distr_t_list F0_d_B_list(UseJack);
  distr_t_list D1_d_B_list(UseJack);
  distr_t_list D2_d_B_list(UseJack);
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
      data[ijack][iens].FF= FB_d_per_ens[iens].distr_list[ixg-1].distr[ijack];
      data[ijack][iens].FF_err= FB_d_per_ens[iens].err(ixg-1);
      if(data_2pts.Tag[iens] == "cA211a.12.48") { data[ijack][iens].a = a_A.distr[ijack]; data[ijack][iens].is=0; }
      else if(data_2pts.Tag[iens] == "cB211b.072.64") { data[ijack][iens].a = a_B.distr[ijack]; data[ijack][iens].is=1;}
      else if(data_2pts.Tag[iens] == "cB211b.072.96") { data[ijack][iens].a = a_B.distr[ijack]; data[ijack][iens].is=1; }
      else if(data_2pts.Tag[iens] == "cC211a.06.80")  {data[ijack][iens].a = a_C.distr[ijack]; data[ijack][iens].is=2; }
      else if(data_2pts.Tag[iens] == "cD211a.054.96")  {data[ijack][iens].a = a_D.distr[ijack]; data[ijack][iens].is=3; }
      else crash("Ens_tag: "+data_2pts.Tag[iens]+" not recognized");

      //mean values
      if(ijack==0) {
	data_ch2[ijack][iens].FF= FB_d_per_ens[iens].ave(ixg-1);
	data_ch2[ijack][iens].FF_err= FB_d_per_ens[iens].err(ixg-1);
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
  cout<<"Fitting FB(u), xg: "<<ixg<<endl;
  Bt_fit= bf_FF.Perform_bootstrap_fit();
  Bt_fit_ch2= bf_FF_ch2.Perform_bootstrap_fit();

  
  //retrieve parameters
  distr_t F0(UseJack), D1(UseJack), D2(UseJack);
  for(int ijack=0;ijack<Njacks;ijack++) { F0.distr.push_back( Bt_fit.par[ijack].F0); D1.distr.push_back( Bt_fit.par[ijack].D1); D2.distr.push_back( Bt_fit.par[ijack].D2);}
  //push_back retrieved parameters
  F0_d_B_list.distr_list.push_back(F0);
  D1_d_B_list.distr_list.push_back(D1);
  D2_d_B_list.distr_list.push_back(D2);
  //push_back ch2
  ch2_FB_d.push_back( Bt_fit_ch2.get_ch2_ave()/dof);

 
  //print fit func
  distr_t_list FB_xg_to_print(UseJack);
  for(auto &a: Bs_a_to_print) FB_xg_to_print.distr_list.push_back( F0 + D1*pow(a*fmTGeV*Lambda_QCD,2) + D2*pow(a*fmTGeV*Lambda_QCD,4));
  Print_To_File({}, {Bs_a_to_print, FB_xg_to_print.ave(), FB_xg_to_print.err()}, "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FB_d_"+Fit_tag+"xg_"+to_string_with_precision(0.10*ixg,2)+".fit_func", "", "#a[fm] FA FA_err");
 
  }


  
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //###########################################                                                                 ##########################################
  //###########################################                                                                 ##########################################
  //###########################################                   FIT FT(d) for all xg                        ##########################################
  //###########################################                                                                 ##########################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  Vfloat ch2_FT_d;
  bf_FF.Set_par_val("F0", -0.1, 0.001);
  bf_FF.Set_par_val("D1", 1.0, 0.1);
  if(Include_a4) bf_FF.Set_par_val("D2", 1.0, 0.1);
  bf_FF_ch2.Set_par_val("F0", -0.1, 0.001);
  bf_FF_ch2.Set_par_val("D1", 1.0, 0.1);
  if(Include_a4) bf_FF_ch2.Set_par_val("D2", 1.0, 0.1);
  distr_t_list F0_d_T_list(UseJack);
  distr_t_list D1_d_T_list(UseJack);
  distr_t_list D2_d_T_list(UseJack);
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
	data[ijack][iens].FF= FT_d_per_ens[iens].distr_list[ixg-1].distr[ijack];
	data[ijack][iens].FF_err= FT_d_per_ens[iens].err(ixg-1);
	if(data_2pts.Tag[iens] == "cA211a.12.48") { data[ijack][iens].a = a_A.distr[ijack]; data[ijack][iens].is=0; }
	else if(data_2pts.Tag[iens] == "cB211b.072.64") { data[ijack][iens].a = a_B.distr[ijack]; data[ijack][iens].is=1;}
	else if(data_2pts.Tag[iens] == "cB211b.072.96") { data[ijack][iens].a = a_B.distr[ijack]; data[ijack][iens].is=1; }
	else if(data_2pts.Tag[iens] == "cC211a.06.80")  {data[ijack][iens].a = a_C.distr[ijack]; data[ijack][iens].is=2; }
	else if(data_2pts.Tag[iens] == "cD211a.054.96")  {data[ijack][iens].a = a_D.distr[ijack]; data[ijack][iens].is=3; }
	else crash("Ens_tag: "+data_2pts.Tag[iens]+" not recognized");

	if(ijack==0) {
	  	data_ch2[ijack][iens].FF= FT_d_per_ens[iens].ave(ixg-1);
		data_ch2[ijack][iens].FF_err= FT_d_per_ens[iens].err(ixg-1);
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
  cout<<"Fitting FT(d), xg: "<<ixg<<endl;
  Bt_fit= bf_FF.Perform_bootstrap_fit();
  Bt_fit_ch2= bf_FF_ch2.Perform_bootstrap_fit();
  //retrieve parameters
  distr_t F0(UseJack), D1(UseJack), D2(UseJack);
  for(int ijack=0;ijack<Njacks;ijack++) { F0.distr.push_back( Bt_fit.par[ijack].F0); D1.distr.push_back( Bt_fit.par[ijack].D1); D2.distr.push_back( Bt_fit.par[ijack].D2);}
  //push_back retrieved parameters
  F0_d_T_list.distr_list.push_back(F0);
  D1_d_T_list.distr_list.push_back(D1);
  D2_d_T_list.distr_list.push_back(D2);
  //push_back ch2
  ch2_FT_d.push_back( Bt_fit_ch2.get_ch2_ave()/dof);

  //print fit func
  distr_t_list FT_xg_to_print(UseJack);
  for(auto &a: Bs_a_to_print) FT_xg_to_print.distr_list.push_back( F0 + D1*pow(a*fmTGeV*Lambda_QCD,2) + D2*pow(a*fmTGeV*Lambda_QCD,4));
  Print_To_File({}, {Bs_a_to_print, FT_xg_to_print.ave(), FT_xg_to_print.err()}, "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FT_d_"+Fit_tag+"xg_"+to_string_with_precision(0.10*ixg,2)+".fit_func", "", "#a[fm] FV FV_err");
  }

  //Print continuum extrapolated form factors
  Print_To_File({}, {Bs_xg_t_list, F0_d_B_list.ave(), F0_d_B_list.err(), (D1_d_B_list/F0_d_B_list).ave(), (D1_d_B_list/F0_d_B_list).err(), (D2_d_B_list/F0_d_B_list).ave(), (D2_d_B_list/F0_d_B_list).err(), ch2_FB_d}, "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FB_d_"+Fit_tag+"cont.dat", "", "#xg  F0   D1   D2   ch2/dof");
  Print_To_File({}, {Bs_xg_t_list, F0_d_T_list.ave(), F0_d_T_list.err(), (D1_d_T_list/F0_d_T_list).ave(), (D1_d_T_list/F0_d_T_list).err(), (D2_d_T_list/F0_d_T_list).ave(), (D2_d_T_list/F0_d_T_list).err(), ch2_FT_d}, "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FT_d_"+Fit_tag+"cont.dat", "", "#xg  F0   D1   D2   ch2/dof");
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////


  

  Eigen::MatrixXd Cov_FB_d(num_xg-1, num_xg-1);
  Eigen::MatrixXd Cov_FT_d(num_xg-1, num_xg-1);
  Eigen::MatrixXd Corr_FB_d(num_xg-1, num_xg-1);
  Eigen::MatrixXd Corr_FT_d(num_xg-1, num_xg-1);

  for(int x_xg=1; x_xg<num_xg;x_xg++) {
    for(int y_xg=1; y_xg<num_xg;y_xg++) {
      Cov_FB_d(x_xg-1, y_xg-1) = F0_d_B_list.distr_list[x_xg-1]%F0_d_B_list.distr_list[y_xg-1];
      Cov_FT_d(x_xg-1, y_xg-1) = F0_d_T_list.distr_list[x_xg-1]%F0_d_T_list.distr_list[y_xg-1];
      Corr_FB_d(x_xg-1, y_xg-1) = Cov_FB_d(x_xg-1,y_xg-1)/(F0_d_B_list.err(x_xg-1)*F0_d_B_list.err(y_xg-1));
      Corr_FT_d(x_xg-1, y_xg-1) = Cov_FT_d(x_xg-1,y_xg-1)/(F0_d_T_list.err(x_xg-1)*F0_d_T_list.err(y_xg-1));
    }
  }

  //Print To File
  ofstream Print_Cov_FB_d("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FB_d_"+Analysis_tag+".cov");
  ofstream Print_Cov_FT_d("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FT_d_"+Analysis_tag+".cov");
  ofstream Print_Corr_FB_d("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FB_d_"+Analysis_tag+".corr");
  ofstream Print_Corr_FT_d("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FT_d_"+Analysis_tag+".corr");

  Print_Cov_FB_d<<Cov_FB_d<<endl;
  Print_Cov_FT_d<<Cov_FT_d<<endl;
  Print_Corr_FB_d<<Corr_FB_d<<endl;
  Print_Corr_FT_d<<Corr_FT_d<<endl;


  Print_Cov_FB_d.close();
  Print_Cov_FT_d.close();
  Print_Corr_FB_d.close();
  Print_Corr_FT_d.close();




  //###################################################################################################################################################################





  

 
  //push_back FIT INFO
  rt_fit.Ch2_FA= ch2_FA;
  rt_fit.Ch2_FV= ch2_FV;
  rt_fit.FA = F0_A_list;
  rt_fit.FV = F0_V_list;
  rt_fit.Ch2_FA_T= ch2_FA_T;
  rt_fit.Ch2_FV_T= ch2_FV_T;
  rt_fit.FA_T = F0_A_T_list;
  rt_fit.FV_T = F0_V_T_list;
  rt_fit.Ch2_FB= ch2_FB;
  rt_fit.Ch2_FT= ch2_FT;
  rt_fit.FT = F0_T_list;
  rt_fit.FB = F0_B_list;

  rt_fit.Ch2_FA_u= ch2_FA_u;
  rt_fit.Ch2_FV_u= ch2_FV_u;
  rt_fit.FA_u = F0_u_A_list;
  rt_fit.FV_u = F0_u_V_list;
  rt_fit.Ch2_FA_T_u= ch2_FA_T_u;
  rt_fit.Ch2_FV_T_u= ch2_FV_T_u;
  rt_fit.FA_T_u = F0_u_A_T_list;
  rt_fit.FV_T_u = F0_u_V_T_list;

  rt_fit.Ch2_FB_u= ch2_FB_u;
  rt_fit.Ch2_FT_u= ch2_FT_u;
  rt_fit.FT_u = F0_u_T_list;
  rt_fit.FB_u = F0_u_B_list;
  
 

  rt_fit.Ch2_FA_d= ch2_FA_d;
  rt_fit.Ch2_FV_d= ch2_FV_d;
  rt_fit.FA_d = F0_d_A_list;
  rt_fit.FV_d = F0_d_V_list;
  rt_fit.Ch2_FA_T_d= ch2_FA_T_d;
  rt_fit.Ch2_FV_T_d= ch2_FV_T_d;
  rt_fit.FA_T_d = F0_d_A_T_list;
  rt_fit.FV_T_d = F0_d_V_T_list;

  rt_fit.Ch2_FB_d= ch2_FB_d;
  rt_fit.Ch2_FT_d= ch2_FT_d;
  rt_fit.FT_d = F0_d_T_list;
  rt_fit.FB_d = F0_d_B_list;

    



  //##########################################################################################
  //##########################################################################################

  //perform continuum extrapolation of MP, fP, , MP/fP and phi
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

  rt_fit.mass = F0_MP;
  rt_fit.Ch2_mass= ch2_MP/dof;

  //begin constrained fits
  //we perform three different fits a^2, a^2 + a^4_tm and a^2 +a^4_OS

  rt_fit.Ndof_K.resize(3);
  rt_fit.Npars_K.resize(3);
  rt_fit.Ndof_K.resize(3);


  vector<string> Tags_comb({"", "_a4_2pt", "_a4_3pt"});
  vector<string> Tags_comb_({"", "a4_2pt_", "a4_3pt_"});
  int c_comb=0;
  string _Fit_tag_std= _Fit_tag;
  string Fit_tag_std=Fit_tag;
  for( auto &tg : Tags_comb) {

    _Fit_tag= _Fit_tag_std+tg;
    Fit_tag= Fit_tag_std+Tags_comb_[c_comb];
    
    class ipar_FF_K {
    public:
      ipar_FF_K() : FF(0.0), FF_err(0.0) {}
      double FF, FF_err, a;
      int is;
      bool Is_2pt;
    };
    
    class fpar_FF_K {
    public:
      fpar_FF_K() {}
      fpar_FF_K(const Vfloat &par) {
	if((signed)par.size() != 5) crash("In class fpar_FF_K  class constructor Vfloat par has size != 5");
	F0=par[0];
	D1=par[1];
	D2=par[2];
	C1=par[3];
	C2=par[4];
      }
      double F0, D1,D2, C1, C2;
    };

  
    //init bootstrap fit
    bootstrap_fit<fpar_FF_K,ipar_FF_K> bf_FF_K(Njacks);
    bf_FF_K.set_warmup_lev(2); //sets warmup
    bf_FF_K.Set_number_of_measurements(2*Nens);
    bf_FF_K.Set_verbosity(1);
    bf_FF_K.Add_par("F0", 0.07, 0.001);
    bf_FF_K.Add_par("D1", -1.0, 0.01);
    bf_FF_K.Add_par("D2", 1.0, 0.01);
    bf_FF_K.Add_par("C1", 1.0, 0.01);
    bf_FF_K.Add_par("C2", 1.0, 0.1);
    //fit on mean values to get ch2
    bootstrap_fit<fpar_FF_K,ipar_FF_K> bf_FF_K_ch2(1);
    bf_FF_K_ch2.set_warmup_lev(2); //sets warmup
    bf_FF_K_ch2.Set_number_of_measurements(2*Nens);
    bf_FF_K_ch2.Set_verbosity(1);
    bf_FF_K_ch2.Add_par("F0", 0.07, 0.001);
    bf_FF_K_ch2.Add_par("D1", -1.0, 0.001);
    bf_FF_K_ch2.Add_par("D2", 1.0, 0.01);
    bf_FF_K_ch2.Add_par("C1", 1.0, 0.01);
    bf_FF_K_ch2.Add_par("C2", 1.0, 0.01);

    if(Fit_single_reg) crash("Only combined fits allowed for FP, phi, FP/MP");
       

    if(Fit_single_reg && Reg_to_fit=="3pt") {
      bf_FF_K.Fix_par("D1", 0.0);
      bf_FF_K.Fix_par("D2", 0.0);
      bf_FF_K_ch2.Fix_par("D1", 0.0);
      bf_FF_K_ch2.Fix_par("D2", 0.0);
    }

    if(Fit_single_reg && Reg_to_fit=="tm") {
      bf_FF_K.Fix_par("C1", 0.0);
      bf_FF_K.Fix_par("C2", 0.0);
      bf_FF_K_ch2.Fix_par("C1", 0.0);
      bf_FF_K_ch2.Fix_par("C2", 0.0);
    }
    
    //count number of dof
    int dof_K= 2*Nens-3;

    if(tg=="") {
      bf_FF_K.Fix_par("D2",0.0); bf_FF_K_ch2.Fix_par("D2",0.0);
      bf_FF_K.Fix_par("C2",0.0); bf_FF_K_ch2.Fix_par("C2",0.0);
    }
    else if(tg=="_a4_2pt") {
      dof_K --;
      bf_FF_K.Fix_par("C2",0.0); bf_FF_K_ch2.Fix_par("C2", 0.0);
    }
    else if(tg=="_a4_3pt") {
      dof_K--;
      bf_FF_K.Fix_par("D2", 0.0); bf_FF_K_ch2.Fix_par("D2", 0.0);
    }
    
    
    //if(Include_a4==false) { bf_FF_K.Fix_par("D2", 0.0); bf_FF_K_ch2.Fix_par("D2",0.0); bf_FF_K.Fix_par("C2", 0.0); bf_FF_K_ch2.Fix_par("C2",0.0);}
    if( Use_three_finest) dof_K-= 2;

    if(Fit_single_reg) dof_K=dof;


  //push_back info
    rt_fit.Nmeas_K[c_comb]=(Fit_single_reg)?rt_fit.Nmeas:(2*Nens-2*(Use_three_finest==true));
    rt_fit.Ndof_K[c_comb]= dof_K;
    rt_fit.Npars_K[c_comb]= rt_fit.Nmeas_K[c_comb]-dof_K;
  

  //ansatz
  bf_FF_K.ansatz=  [&Use_three_finest ](const fpar_FF_K &p, const ipar_FF_K &ip) {

    if( (ip.is==0) && Use_three_finest) return 0.0;

    double D1= p.D1;
    double D2= p.D2;
    
    if(ip.Is_2pt == false) {D1=p.C1; D2=p.C2;}

    
    if(Fit_single_reg) {
      if(ip.Is_2pt && Reg_to_fit=="3pt") return 0.0;
      if(!ip.Is_2pt && Reg_to_fit=="tm") return 0.0;
    }
		  
    return p.F0 + D1*pow(ip.a*Lambda_QCD,2) + D2*pow(ip.a*Lambda_QCD,4);

  };

  
  bf_FF_K.measurement=  [&Use_three_finest ](const fpar_FF_K &p, const ipar_FF_K &ip) {
    
    if( (ip.is==0) && Use_three_finest) return 0.0;

    if(Fit_single_reg) {
      if(ip.Is_2pt && Reg_to_fit=="3pt") return 0.0;
      if(!ip.Is_2pt && Reg_to_fit=="tm") return 0.0;
    }
    
    return ip.FF;
       
  };
  bf_FF_K.error=  [ ](const fpar_FF_K &p, const ipar_FF_K &ip) {

		 return ip.FF_err;
  };

  bf_FF_K_ch2.ansatz= bf_FF_K.ansatz;
  bf_FF_K_ch2.measurement = bf_FF_K.measurement;
  bf_FF_K_ch2.error = bf_FF_K.error;

  


 

  //FP
  double ch2_FP;
  bf_FF_K.Set_par_val("F0", FP_list.ave(0), FP_list.ave(0)/10);
  bf_FF_K.Set_par_val("D1", 1.0/FP_list.ave(0), 0.03/FP_list.ave(0));
  bf_FF_K.Set_par_val("C1", 1.0/FP_list.ave(0), 0.03/FP_list.ave(0));
  if(Include_a4) { bf_FF_K.Set_par_val("D2", 1.0/pow(FP_list.ave(0),2), 0.03/pow(FP_list.ave(0),2));  bf_FF_K.Set_par_val("C2", 1.0/pow(FP_list.ave(0),2), 0.03/pow(FP_list.ave(0),2)); };
  bf_FF_K_ch2.Set_par_val("F0", FP_list.ave(0), FP_list.ave(0)/10);
  bf_FF_K_ch2.Set_par_val("D1", 1.0/FP_list.ave(0), 0.03/FP_list.ave(0));
  bf_FF_K_ch2.Set_par_val("C1", 1.0/FP_list.ave(0), 0.03/FP_list.ave(0));
  if(Include_a4) { bf_FF_K_ch2.Set_par_val("D2", 1.0/pow(FP_list.ave(0),2), 0.03/pow(FP_list.ave(0),2));  bf_FF_K_ch2.Set_par_val("C2", 1.0/pow(FP_list.ave(0),2), 0.03/pow(FP_list.ave(0),2));}
  distr_t F0_FP(UseJack);
  distr_t D1_FP(UseJack);
  distr_t D2_FP(UseJack);
  distr_t C1_FP(UseJack);
  distr_t C2_FP(UseJack);
  vector<vector<ipar_FF_K>> data_FP(Njacks);
  vector<vector<ipar_FF_K>> data_FP_ch2(1);
  //allocate space for output result
  boot_fit_data<fpar_FF_K> Bt_fit_FP;
  boot_fit_data<fpar_FF_K> Bt_fit_FP_ch2;

  //Add covariance matrix
  int Nmeas= 2*Nens;
  
  Eigen::MatrixXd Cov_Matrix_FP(Nmeas,Nmeas);
  Eigen::MatrixXd Corr_Matrix_FP(Nmeas,Nmeas);
  for(int i=0;i<Nmeas;i++) for(int j=0;j<Nmeas;j++) {Cov_Matrix_FP(i,j)=0; Corr_Matrix_FP(i,j)=0;}
   for(int iens=0; iens<Nens;iens++) {
     Corr_Matrix_FP(iens,iens) = 1;
     Corr_Matrix_FP(iens+Nens,iens+Nens) = 1;
     Cov_Matrix_FP(iens,iens) = pow(FP_list.err(iens),2);
     Cov_Matrix_FP(iens+Nens,iens+Nens) = pow(FP_3pt_list.err(iens),2);
     Cov_Matrix_FP(iens, iens+Nens) = FP_list.distr_list[iens]%FP_3pt_list.distr_list[iens];
     Cov_Matrix_FP(iens+Nens, iens) = Cov_Matrix_FP(iens,iens+Nens);
     Corr_Matrix_FP(iens, iens+Nens) = Cov_Matrix_FP(iens, iens+Nens)/(FP_list.err(iens)*FP_3pt_list.err(iens));
     Corr_Matrix_FP(iens+Nens, iens) = Corr_Matrix_FP(iens,iens+Nens);
   }

   //add cov matrix to bootstrap fit
   bf_FF_K.Add_covariance_matrix(Cov_Matrix_FP);
   bf_FF_K_ch2.Add_covariance_matrix(Cov_Matrix_FP);

   //print correlation matrix
   //print covariance matrix
   ofstream Print_Cov_FP("../data/ph_emission/"+ph_type+"/"+MESON+"/decay_const/corr_matrix_fit/FP_"+Fit_tag+".cov");
   ofstream Print_Corr_FP("../data/ph_emission/"+ph_type+"/"+MESON+"/decay_const/corr_matrix_fit/FP_"+Fit_tag+".corr");
   

   Print_Cov_FP<<Cov_Matrix_FP<<endl;  Print_Corr_FP<<Corr_Matrix_FP<<endl;
   Print_Cov_FP.close();               Print_Corr_FP.close();

   
  for(auto &data_iboot: data_FP) data_iboot.resize(2*Nens);
  for(auto &data_iboot: data_FP_ch2) data_iboot.resize(2*Nens);
  for(int ijack=0;ijack<Njacks;ijack++) {
    for(int iens=0;iens<Nens;iens++) {
      data_FP[ijack][iens].FF= FP_list.distr_list[iens].distr[ijack];
      data_FP[ijack][iens].FF_err= FP_list.err(iens);
      data_FP[ijack][iens].Is_2pt = true;
      data_FP[ijack][iens+Nens].FF= FP_3pt_list.distr_list[iens].distr[ijack];
      data_FP[ijack][iens+Nens].FF_err= FP_3pt_list.err(iens);
      data_FP[ijack][iens+Nens].Is_2pt=false;
      if(data_2pts.Tag[iens] == "cA211a.12.48") { data_FP[ijack][iens].a = a_A.distr[ijack]; data_FP[ijack][iens].is=0; }
      else if(data_2pts.Tag[iens] == "cB211b.072.64") { data_FP[ijack][iens].a = a_B.distr[ijack]; data_FP[ijack][iens].is=1;}
      else if(data_2pts.Tag[iens] == "cB211b.072.96") { data_FP[ijack][iens].a = a_B.distr[ijack]; data_FP[ijack][iens].is=1; }
      else if(data_2pts.Tag[iens] == "cC211a.06.80")  {data_FP[ijack][iens].a = a_C.distr[ijack]; data_FP[ijack][iens].is=2; }
      else if(data_2pts.Tag[iens] == "cD211a.054.96")  {data_FP[ijack][iens].a = a_D.distr[ijack]; data_FP[ijack][iens].is=3; }
      else crash("Ens_tag: "+data_2pts.Tag[iens]+" not recognized");

      data_FP[ijack][iens+Nens].is= data_FP[ijack][iens].is;
      data_FP[ijack][iens+Nens].a = data_FP[ijack][iens].a;
      
      if(ijack==0) {
	data_FP_ch2[ijack][iens].FF= FP_list.ave(iens);
	data_FP_ch2[ijack][iens].FF_err= FP_list.err(iens);
	data_FP_ch2[ijack][iens].Is_2pt = true;
	data_FP_ch2[ijack][iens+Nens].FF= FP_3pt_list.ave(iens);
	data_FP_ch2[ijack][iens+Nens].FF_err= FP_3pt_list.err(iens);
	data_FP_ch2[ijack][iens+Nens].Is_2pt=false;
	if(data_2pts.Tag[iens] == "cA211a.12.48") { data_FP_ch2[ijack][iens].a = a_A.distr[ijack]; data_FP_ch2[ijack][iens].is=0; }
	else if(data_2pts.Tag[iens] == "cB211b.072.64") { data_FP_ch2[ijack][iens].a = a_B.distr[ijack]; data_FP_ch2[ijack][iens].is=1;}
	else if(data_2pts.Tag[iens] == "cB211b.072.96") { data_FP_ch2[ijack][iens].a = a_B.distr[ijack]; data_FP_ch2[ijack][iens].is=1; }
	else if(data_2pts.Tag[iens] == "cC211a.06.80")  {data_FP_ch2[ijack][iens].a = a_C.distr[ijack]; data_FP_ch2[ijack][iens].is=2; }
	else if(data_2pts.Tag[iens] == "cD211a.054.96")  {data_FP_ch2[ijack][iens].a = a_D.distr[ijack]; data_FP_ch2[ijack][iens].is=3; }
	else crash("Ens_tag: "+data_2pts.Tag[iens]+" not recognized");

	data_FP_ch2[ijack][iens+Nens].is= data_FP_ch2[ijack][iens].is;
	data_FP_ch2[ijack][iens+Nens].a = data_FP_ch2[ijack][iens].a;
      
	
      }
    }
  }
  //append
  bf_FF_K.Append_to_input_par(data_FP);
  bf_FF_K_ch2.Append_to_input_par(data_FP_ch2);
  //fit
  cout<<"Fitting FP"<<endl;
  Bt_fit_FP= bf_FF_K.Perform_bootstrap_fit();
  Bt_fit_FP_ch2= bf_FF_K_ch2.Perform_bootstrap_fit();
  //retrieve parameters
  for(int ijack=0;ijack<Njacks;ijack++) {
    F0_FP.distr.push_back( Bt_fit_FP.par[ijack].F0);
    D1_FP.distr.push_back( Bt_fit_FP.par[ijack].D1);
    D2_FP.distr.push_back( Bt_fit_FP.par[ijack].D2);
    C1_FP.distr.push_back( Bt_fit_FP.par[ijack].C1);
    C2_FP.distr.push_back( Bt_fit_FP.par[ijack].C2);

  }
  //get ch2
  ch2_FP = Bt_fit_FP_ch2.get_ch2_ave();


  


  //phi
  double ch2_phi;
  bf_FF_K.Set_par_val("F0", FP_list.ave(0)*sqrt(MP_list.ave(0)), FP_list.ave(0)*sqrt(MP_list.ave(0))/10);
  bf_FF_K.Set_par_val("D1", 1.0, 0.1);
  bf_FF_K.Set_par_val("C1", 1.0, 0.1);
  if(Include_a4) { bf_FF_K.Set_par_val("D2", 1.0, 0.1); bf_FF_K.Set_par_val("C2", 1.0, 0.1); }
  bf_FF_K_ch2.Set_par_val("F0", FP_list.ave(0)*sqrt(MP_list.ave(0)), FP_list.ave(0)*sqrt(MP_list.ave(0))/10);
  bf_FF_K_ch2.Set_par_val("D1", 1.0, 0.1);
  bf_FF_K_ch2.Set_par_val("C1", 1.0, 0.1);
  if(Include_a4) { bf_FF_K_ch2.Set_par_val("D2", 1.0, 0.1); bf_FF_K_ch2.Set_par_val("C2", 1.0, 0.1); }
  distr_t F0_phi(UseJack);
  distr_t D1_phi(UseJack);
  distr_t D2_phi(UseJack);
  distr_t C1_phi(UseJack);
  distr_t C2_phi(UseJack);
  vector<vector<ipar_FF_K>> data_phi(Njacks);
  vector<vector<ipar_FF_K>> data_phi_ch2(1);
  //allocate space for output result
  boot_fit_data<fpar_FF_K> Bt_fit_phi;
  boot_fit_data<fpar_FF_K> Bt_fit_phi_ch2;


  Eigen::MatrixXd Cov_Matrix_phi(Nmeas,Nmeas);
  Eigen::MatrixXd Corr_Matrix_phi(Nmeas,Nmeas);
  for(int i=0;i<Nmeas;i++) for(int j=0;j<Nmeas;j++) {Cov_Matrix_phi(i,j)=0; Corr_Matrix_phi(i,j)=0;}
   for(int iens=0; iens<Nens;iens++) {
     Corr_Matrix_phi(iens,iens) = 1;
     Corr_Matrix_phi(iens+Nens,iens+Nens) = 1;
     Cov_Matrix_phi(iens,iens) = pow((FP_list.distr_list[iens]*SQRT_D(MP_list.distr_list[iens])).err(),2);
     Cov_Matrix_phi(iens+Nens,iens+Nens) = pow((FP_3pt_list.distr_list[iens]*SQRT_D(MP_list.distr_list[iens])).err(),2);
     Cov_Matrix_phi(iens, iens+Nens) = (FP_list.distr_list[iens]*SQRT_D(MP_list.distr_list[iens]))%(FP_3pt_list.distr_list[iens]*SQRT_D(MP_list.distr_list[iens]));
     Cov_Matrix_phi(iens+Nens, iens) = Cov_Matrix_phi(iens,iens+Nens);
     Corr_Matrix_phi(iens, iens+Nens) = Cov_Matrix_phi(iens, iens+Nens)/(sqrt( Cov_Matrix_phi(iens,iens)*Cov_Matrix_phi(iens+Nens,iens+Nens) ));
     Corr_Matrix_phi(iens+Nens, iens) = Corr_Matrix_phi(iens,iens+Nens);
   }

   //add cov matrix to bootstrap fit
   bf_FF_K.Add_covariance_matrix(Cov_Matrix_phi);
   bf_FF_K_ch2.Add_covariance_matrix(Cov_Matrix_phi);

   //print correlation matrix
   //print covariance matrix
   ofstream Print_Cov_phi("../data/ph_emission/"+ph_type+"/"+MESON+"/decay_const/corr_matrix_fit/phi_"+Fit_tag+".cov");
   ofstream Print_Corr_phi("../data/ph_emission/"+ph_type+"/"+MESON+"/decay_const/corr_matrix_fit/phi_"+Fit_tag+".corr");
   

   Print_Cov_phi<<Cov_Matrix_phi<<endl;  Print_Corr_phi<<Corr_Matrix_phi<<endl;
   Print_Cov_phi.close();               Print_Corr_phi.close();


  
  for(auto &data_iboot: data_phi) data_iboot.resize(2*Nens);
  for(auto &data_iboot: data_phi_ch2) data_iboot.resize(2*Nens);
  for(int ijack=0;ijack<Njacks;ijack++) {
    for(int iens=0;iens<Nens;iens++) {
      data_phi[ijack][iens].FF= FP_list.distr_list[iens].distr[ijack]*sqrt(MP_list.distr_list[iens].distr[ijack]);
      data_phi[ijack][iens].FF_err= (FP_list.distr_list[iens]*SQRT_D(MP_list.distr_list[iens])).err();
      data_phi[ijack][iens].Is_2pt=true;
      data_phi[ijack][iens+Nens].FF= FP_3pt_list.distr_list[iens].distr[ijack]*sqrt(MP_list.distr_list[iens].distr[ijack]);
      data_phi[ijack][iens+Nens].FF_err= (FP_3pt_list.distr_list[iens]*SQRT_D(MP_list.distr_list[iens])).err();
      data_phi[ijack][iens+Nens].Is_2pt=false;
      if(data_2pts.Tag[iens] == "cA211a.12.48") { data_phi[ijack][iens].a = a_A.distr[ijack]; data_phi[ijack][iens].is=0; }
      else if(data_2pts.Tag[iens] == "cB211b.072.64") { data_phi[ijack][iens].a = a_B.distr[ijack]; data_phi[ijack][iens].is=1;}
      else if(data_2pts.Tag[iens] == "cB211b.072.96") { data_phi[ijack][iens].a = a_B.distr[ijack]; data_phi[ijack][iens].is=1; }
      else if(data_2pts.Tag[iens] == "cC211a.06.80")  {data_phi[ijack][iens].a = a_C.distr[ijack]; data_phi[ijack][iens].is=2; }
      else if(data_2pts.Tag[iens] == "cD211a.054.96")  {data_phi[ijack][iens].a = a_D.distr[ijack]; data_phi[ijack][iens].is=3; }
      else crash("Ens_tag: "+data_2pts.Tag[iens]+" not recognized");

      data_phi[ijack][iens+Nens].is= data_phi[ijack][iens].is;
      data_phi[ijack][iens+Nens].a= data_phi[ijack][iens].a;

      if(ijack==0) {
	data_phi_ch2[ijack][iens].FF= (FP_list.distr_list[iens]*SQRT_D(MP_list.distr_list[iens])).ave();
	data_phi_ch2[ijack][iens].FF_err= (FP_list.distr_list[iens]*SQRT_D(MP_list.distr_list[iens])).err();
	data_phi_ch2[ijack][iens].Is_2pt=true;
	data_phi_ch2[ijack][iens+Nens].FF= (FP_3pt_list.distr_list[iens]*SQRT_D(MP_list.distr_list[iens])).ave();
	data_phi_ch2[ijack][iens+Nens].FF_err= (FP_3pt_list.distr_list[iens]*SQRT_D(MP_list.distr_list[iens])).err();
	data_phi_ch2[ijack][iens+Nens].Is_2pt=false;
	if(data_2pts.Tag[iens] == "cA211a.12.48") { data_phi_ch2[ijack][iens].a = a_A.distr[ijack]; data_phi_ch2[ijack][iens].is=0; }
	else if(data_2pts.Tag[iens] == "cB211b.072.64") { data_phi_ch2[ijack][iens].a = a_B.distr[ijack]; data_phi_ch2[ijack][iens].is=1;}
	else if(data_2pts.Tag[iens] == "cB211b.072.96") { data_phi_ch2[ijack][iens].a = a_B.distr[ijack]; data_phi_ch2[ijack][iens].is=1; }
	else if(data_2pts.Tag[iens] == "cC211a.06.80")  {data_phi_ch2[ijack][iens].a = a_C.distr[ijack]; data_phi_ch2[ijack][iens].is=2; }
	else if(data_2pts.Tag[iens] == "cD211a.054.96")  {data_phi_ch2[ijack][iens].a = a_D.distr[ijack]; data_phi_ch2[ijack][iens].is=3; }
	else crash("Ens_tag: "+data_2pts.Tag[iens]+" not recognized");

	data_phi_ch2[ijack][iens+Nens].is= data_phi_ch2[ijack][iens].is;
	data_phi_ch2[ijack][iens+Nens].a= data_phi_ch2[ijack][iens].a;
	
      }
    }
  }
  //append
  bf_FF_K.Append_to_input_par(data_phi);
  bf_FF_K_ch2.Append_to_input_par(data_phi_ch2);
  //fit
  cout<<"Fitting phi"<<endl;
  Bt_fit_phi= bf_FF_K.Perform_bootstrap_fit();
  Bt_fit_phi_ch2= bf_FF_K_ch2.Perform_bootstrap_fit();
  //retrieve parameters
  for(int ijack=0;ijack<Njacks;ijack++) {
    F0_phi.distr.push_back( Bt_fit_phi.par[ijack].F0);
    D1_phi.distr.push_back( Bt_fit_phi.par[ijack].D1);
    D2_phi.distr.push_back( Bt_fit_phi.par[ijack].D2);
    C1_phi.distr.push_back( Bt_fit_phi.par[ijack].C1);
    C2_phi.distr.push_back( Bt_fit_phi.par[ijack].C2);
  }
  //get ch2
  ch2_phi = Bt_fit_phi_ch2.get_ch2_ave();
  

  
  


  //MP/FP
  double ch2_MP_ov_FP;
  bf_FF_K.Set_par_val("F0", MP_ov_FP_list.ave(0), MP_ov_FP_list.ave(0)/10);
  bf_FF_K.Set_par_val("D1", 1.0, 0.1);
  bf_FF_K.Set_par_val("C1", 1.0, 0.1);
  if(Include_a4) { bf_FF_K.Set_par_val("D2", 1.0, 0.1); bf_FF_K.Set_par_val("C2", 1.0, 0.1); }
  bf_FF_K_ch2.Set_par_val("F0", MP_ov_FP_list.ave(0), MP_ov_FP_list.ave(0)/10);
  bf_FF_K_ch2.Set_par_val("D1", 1.0, 0.1);
  bf_FF_K_ch2.Set_par_val("C1", 1.0, 0.1);
  if(Include_a4) { bf_FF_K_ch2.Set_par_val("D2", 1.0, 0.1); bf_FF_K_ch2.Set_par_val("C2", 1.0, 0.1); }
  distr_t F0_MP_ov_FP(UseJack);
  distr_t D1_MP_ov_FP(UseJack);
  distr_t D2_MP_ov_FP(UseJack);
  distr_t C1_MP_ov_FP(UseJack);
  distr_t C2_MP_ov_FP(UseJack);
  vector<vector<ipar_FF_K>> data_MP_ov_FP(Njacks);
  vector<vector<ipar_FF_K>> data_MP_ov_FP_ch2(1);
  //allocate space for output result
  boot_fit_data<fpar_FF_K> Bt_fit_MP_ov_FP;
  boot_fit_data<fpar_FF_K> Bt_fit_MP_ov_FP_ch2;
  for(auto &data_iboot: data_MP_ov_FP) data_iboot.resize(2*Nens);
  for(auto &data_iboot: data_MP_ov_FP_ch2) data_iboot.resize(2*Nens);

  Eigen::MatrixXd Cov_Matrix_MP_ov_FP(Nmeas,Nmeas);
  Eigen::MatrixXd Corr_Matrix_MP_ov_FP(Nmeas,Nmeas);
  for(int i=0;i<Nmeas;i++) for(int j=0;j<Nmeas;j++) {Cov_Matrix_MP_ov_FP(i,j)=0; Corr_Matrix_MP_ov_FP(i,j)=0;}
   for(int iens=0; iens<Nens;iens++) {
     Corr_Matrix_MP_ov_FP(iens,iens) = 1;
     Corr_Matrix_MP_ov_FP(iens+Nens,iens+Nens) = 1;
     Cov_Matrix_MP_ov_FP(iens,iens) = pow(MP_ov_FP_list.err(iens),2);
     Cov_Matrix_MP_ov_FP(iens+Nens,iens+Nens) = pow(MP_ov_FP_3pt_list.err(iens),2);
     Cov_Matrix_MP_ov_FP(iens, iens+Nens) = MP_ov_FP_list.distr_list[iens]%MP_ov_FP_3pt_list.distr_list[iens];
     Cov_Matrix_MP_ov_FP(iens+Nens, iens) = Cov_Matrix_MP_ov_FP(iens,iens+Nens);
     Corr_Matrix_MP_ov_FP(iens, iens+Nens) = Cov_Matrix_MP_ov_FP(iens, iens+Nens)/(MP_ov_FP_list.err(iens)*MP_ov_FP_3pt_list.err(iens));
     Corr_Matrix_MP_ov_FP(iens+Nens, iens) = Corr_Matrix_MP_ov_FP(iens,iens+Nens);
   }

   //add cov matrix to bootstrap fit
   bf_FF_K.Add_covariance_matrix(Cov_Matrix_MP_ov_FP);
   bf_FF_K_ch2.Add_covariance_matrix(Cov_Matrix_MP_ov_FP);

   //print correlation matrix
   //print covariance matrix
   ofstream Print_Cov_MP_ov_FP("../data/ph_emission/"+ph_type+"/"+MESON+"/decay_const/corr_matrix_fit/mP_ov_fP_"+Fit_tag+".cov");
   ofstream Print_Corr_MP_ov_FP("../data/ph_emission/"+ph_type+"/"+MESON+"/decay_const/corr_matrix_fit/mP_ov_fP_"+Fit_tag+".corr");
   

   Print_Cov_MP_ov_FP<<Cov_Matrix_MP_ov_FP<<endl;  Print_Corr_MP_ov_FP<<Corr_Matrix_MP_ov_FP<<endl;
   Print_Cov_MP_ov_FP.close();                     Print_Corr_MP_ov_FP.close();

  
   for(int ijack=0;ijack<Njacks;ijack++) {
     for(int iens=0;iens<Nens;iens++) {
       data_MP_ov_FP[ijack][iens].FF= MP_ov_FP_list.distr_list[iens].distr[ijack];
       data_MP_ov_FP[ijack][iens].FF_err= MP_ov_FP_list.err(iens);
       data_MP_ov_FP[ijack][iens].Is_2pt=true;
       data_MP_ov_FP[ijack][iens+Nens].FF= MP_ov_FP_3pt_list.distr_list[iens].distr[ijack];
       data_MP_ov_FP[ijack][iens+Nens].FF_err= MP_ov_FP_3pt_list.err(iens);
       data_MP_ov_FP[ijack][iens+Nens].Is_2pt=false;
       if(data_2pts.Tag[iens] == "cA211a.12.48") { data_MP_ov_FP[ijack][iens].a = a_A.distr[ijack]; data_MP_ov_FP[ijack][iens].is=0; }
       else if(data_2pts.Tag[iens] == "cB211b.072.64") { data_MP_ov_FP[ijack][iens].a = a_B.distr[ijack]; data_MP_ov_FP[ijack][iens].is=1;}
       else if(data_2pts.Tag[iens] == "cB211b.072.96") { data_MP_ov_FP[ijack][iens].a = a_B.distr[ijack]; data_MP_ov_FP[ijack][iens].is=1; }
       else if(data_2pts.Tag[iens] == "cC211a.06.80")  {data_MP_ov_FP[ijack][iens].a = a_C.distr[ijack]; data_MP_ov_FP[ijack][iens].is=2; }
       else if(data_2pts.Tag[iens] == "cD211a.054.96")  {data_MP_ov_FP[ijack][iens].a = a_D.distr[ijack]; data_MP_ov_FP[ijack][iens].is=3; }
       else crash("Ens_tag: "+data_2pts.Tag[iens]+" not recognized");

       data_MP_ov_FP[ijack][iens+Nens].is= data_MP_ov_FP[ijack][iens].is;
       data_MP_ov_FP[ijack][iens+Nens].a= data_MP_ov_FP[ijack][iens].a;

      if(ijack==0) {
	data_MP_ov_FP_ch2[ijack][iens].FF= MP_ov_FP_list.ave(iens);
	data_MP_ov_FP_ch2[ijack][iens].FF_err= MP_ov_FP_list.err(iens);
	data_MP_ov_FP_ch2[ijack][iens].Is_2pt=true;
	data_MP_ov_FP_ch2[ijack][iens+Nens].FF= MP_ov_FP_3pt_list.ave(iens);
	data_MP_ov_FP_ch2[ijack][iens+Nens].FF_err= MP_ov_FP_3pt_list.err(iens);
	data_MP_ov_FP_ch2[ijack][iens+Nens].Is_2pt=false;
	if(data_2pts.Tag[iens] == "cA211a.12.48") { data_MP_ov_FP_ch2[ijack][iens].a = a_A.distr[ijack]; data_MP_ov_FP_ch2[ijack][iens].is=0; }
	else if(data_2pts.Tag[iens] == "cB211b.072.64") { data_MP_ov_FP_ch2[ijack][iens].a = a_B.distr[ijack]; data_MP_ov_FP_ch2[ijack][iens].is=1;}
	else if(data_2pts.Tag[iens] == "cB211b.072.96") { data_MP_ov_FP_ch2[ijack][iens].a = a_B.distr[ijack]; data_MP_ov_FP_ch2[ijack][iens].is=1; }
	else if(data_2pts.Tag[iens] == "cC211a.06.80")  {data_MP_ov_FP_ch2[ijack][iens].a = a_C.distr[ijack]; data_MP_ov_FP_ch2[ijack][iens].is=2; }
	else if(data_2pts.Tag[iens] == "cD211a.054.96")  {data_MP_ov_FP_ch2[ijack][iens].a = a_D.distr[ijack]; data_MP_ov_FP_ch2[ijack][iens].is=3; }
	else crash("Ens_tag: "+data_2pts.Tag[iens]+" not recognized");

	data_MP_ov_FP_ch2[ijack][iens+Nens].is= data_MP_ov_FP_ch2[ijack][iens].is;
	data_MP_ov_FP_ch2[ijack][iens+Nens].a= data_MP_ov_FP_ch2[ijack][iens].a;
      }
    }
  }
  //append
  bf_FF_K.Append_to_input_par(data_MP_ov_FP);
  bf_FF_K_ch2.Append_to_input_par(data_MP_ov_FP_ch2);
  //fit
  cout<<"Fitting MP/FP"<<endl;
  Bt_fit_MP_ov_FP= bf_FF_K.Perform_bootstrap_fit();
  Bt_fit_MP_ov_FP_ch2= bf_FF_K_ch2.Perform_bootstrap_fit();
  //retrieve parameters
  for(int ijack=0;ijack<Njacks;ijack++) {
    F0_MP_ov_FP.distr.push_back( Bt_fit_MP_ov_FP.par[ijack].F0);
    D1_MP_ov_FP.distr.push_back( Bt_fit_MP_ov_FP.par[ijack].D1);
    D2_MP_ov_FP.distr.push_back(Bt_fit_MP_ov_FP.par[ijack].D2);
    C1_MP_ov_FP.distr.push_back( Bt_fit_MP_ov_FP.par[ijack].C1);
    C2_MP_ov_FP.distr.push_back(Bt_fit_MP_ov_FP.par[ijack].C2);
  }
  //get ch2
  ch2_MP_ov_FP = Bt_fit_MP_ov_FP_ch2.get_ch2_ave();


  //print fitting functions for MP, FP, MP_ov_FP
  distr_t_list MP_to_print(UseJack);
  distr_t_list FP_to_print(UseJack);
  distr_t_list phi_to_print(UseJack);
  distr_t_list MP_ov_FP_to_print(UseJack);
  distr_t_list FP_3pt_to_print(UseJack);
  distr_t_list phi_3pt_to_print(UseJack);
  distr_t_list MP_ov_FP_3pt_to_print(UseJack);
  for(auto &a: Bs_a_to_print) {
    MP_to_print.distr_list.push_back( F0_MP + D1_MP*pow(a*fmTGeV*Lambda_QCD,2) + D2_MP*pow(a*fmTGeV*Lambda_QCD,4));
    FP_to_print.distr_list.push_back( F0_FP + D1_FP*pow(a*fmTGeV*Lambda_QCD,2) + D2_FP*pow(a*fmTGeV*Lambda_QCD,4));
    phi_to_print.distr_list.push_back( F0_phi + D1_phi*pow(a*fmTGeV*Lambda_QCD,2) + D2_phi*pow(a*fmTGeV*Lambda_QCD,4));
    MP_ov_FP_to_print.distr_list.push_back( F0_MP_ov_FP + D1_MP_ov_FP*pow(a*fmTGeV*Lambda_QCD,2) + D2_MP_ov_FP*pow(a*fmTGeV*Lambda_QCD,4));
    FP_3pt_to_print.distr_list.push_back( F0_FP + C1_FP*pow(a*fmTGeV*Lambda_QCD,2) + C2_FP*pow(a*fmTGeV*Lambda_QCD,4));
    phi_3pt_to_print.distr_list.push_back( F0_phi + C1_phi*pow(a*fmTGeV*Lambda_QCD,2) + C2_phi*pow(a*fmTGeV*Lambda_QCD,4));
    MP_ov_FP_3pt_to_print.distr_list.push_back( F0_MP_ov_FP + C1_MP_ov_FP*pow(a*fmTGeV*Lambda_QCD,2) + C2_MP_ov_FP*pow(a*fmTGeV*Lambda_QCD,4));
  }
  Print_To_File({}, {Bs_a_to_print, MP_to_print.ave(), MP_to_print.err()}, "../data/ph_emission/"+ph_type+"/"+MESON+"/mass/masses"+_Fit_tag+".fit_func", "", "#a[fm]  MP MP_err");
  Print_To_File({}, {Bs_a_to_print, FP_to_print.ave(), FP_to_print.err(), FP_3pt_to_print.ave(), FP_3pt_to_print.err()} , "../data/ph_emission/"+ph_type+"/"+MESON+"/decay_const/fP"+_Fit_tag+".fit_func", "", "#a[fm]  FP FP_err   MP/FP  MP/FP_err");
  Print_To_File({}, {Bs_a_to_print, MP_ov_FP_to_print.ave(), MP_ov_FP_to_print.err(), MP_ov_FP_3pt_to_print.ave(), MP_ov_FP_3pt_to_print.err()} , "../data/ph_emission/"+ph_type+"/"+MESON+"/decay_const/mP_ov_fP"+_Fit_tag+".fit_func", "", "#a[fm]  FP FP_err   MP/FP  MP/FP_err");
  Print_To_File({}, {Bs_a_to_print, phi_to_print.ave(), phi_to_print.err(), phi_3pt_to_print.ave(), phi_3pt_to_print.err()}, "../data/ph_emission/"+ph_type+"/"+MESON+"/decay_const/phi"+_Fit_tag+".fit_func", "", "#a[fm]  phi phi_err");

  //push back the result for the combined fits
 
  rt_fit.fp[c_comb] = F0_FP;
  rt_fit.Ch2_fp[c_comb]= ch2_FP/dof_K;
  rt_fit.phi[c_comb] = F0_phi;
  rt_fit.mp_ov_fp[c_comb] = F0_MP_ov_FP;
  rt_fit.Ch2_phi[c_comb] = ch2_phi/dof_K;
  rt_fit.Ch2_mp_ov_fp[c_comb]= ch2_MP_ov_FP/dof_K;


  cout<<"#### Tag: "<<tg<<"#####"<<endl;
  cout<<"ch2(MP): "<< ch2_MP/dof<<endl;
  cout<<"ch2(FP): "<< ch2_FP/dof_K<<endl;
  cout<<"ch2(phi): "<< ch2_phi/dof_K<<endl;
  cout<<"ch2(MP/FP): "<<ch2_MP_ov_FP/dof_K<<endl;
  cout<<"###################################"<<endl;
  c_comb++; 
  }

  //print ch2 summa
  
  ofstream print_ch2("ch2_"+MESON+"_"+Analysis_tag+".txt");

  
  print_ch2<<"#######   SUMMARY OF Ch2 ##########"<<endl;
  print_ch2<<"### FA u ###"<<endl;
  for(int ixg=1;ixg<num_xg;ixg++) { print_ch2<<"ch2(xg: "<<Bs_xg_t_list[ixg-1]<<"): "<<ch2_FA_u[ixg-1]<<endl; }
  print_ch2<<"### FA d ###"<<endl;
  for(int ixg=1;ixg<num_xg;ixg++) { print_ch2<<"ch2(xg: "<<Bs_xg_t_list[ixg-1]<<"): "<<ch2_FA_d[ixg-1]<<endl; }
  print_ch2<<"### FV u ###"<<endl;
  for(int ixg=1;ixg<num_xg;ixg++) { print_ch2<<"ch2(xg: "<<Bs_xg_t_list[ixg-1]<<"): "<<ch2_FV_u[ixg-1]<<endl; }
  print_ch2<<"### FV d ###"<<endl;
  for(int ixg=1;ixg<num_xg;ixg++) { print_ch2<<"ch2(xg: "<<Bs_xg_t_list[ixg-1]<<"): "<<ch2_FV_d[ixg-1]<<endl; }
  print_ch2<<"### FAT u ###"<<endl;
  for(int ixg=1;ixg<num_xg;ixg++) { print_ch2<<"ch2(xg: "<<Bs_xg_t_list[ixg-1]<<"): "<<ch2_FA_T_u[ixg-1]<<endl; }
  print_ch2<<"### FAT d ###"<<endl;
  for(int ixg=1;ixg<num_xg;ixg++) { print_ch2<<"ch2(xg: "<<Bs_xg_t_list[ixg-1]<<"): "<<ch2_FA_T_d[ixg-1]<<endl; }
  print_ch2<<"### FVT u ###"<<endl;
  for(int ixg=1;ixg<num_xg;ixg++) { print_ch2<<"ch2(xg: "<<Bs_xg_t_list[ixg-1]<<"): "<<ch2_FV_T_u[ixg-1]<<endl; }
  print_ch2<<"### FVT d ###"<<endl;
  for(int ixg=1;ixg<num_xg;ixg++) { print_ch2<<"ch2(xg: "<<Bs_xg_t_list[ixg-1]<<"): "<<ch2_FV_T_d[ixg-1]<<endl; }
  print_ch2<<"### FB u ###"<<endl;
  for(int ixg=1;ixg<num_xg;ixg++) { print_ch2<<"ch2(xg: "<<Bs_xg_t_list[ixg-1]<<"): "<<ch2_FB_u[ixg-1]<<endl; }
  print_ch2<<"### FB d ###"<<endl;
  for(int ixg=1;ixg<num_xg;ixg++) { print_ch2<<"ch2(xg: "<<Bs_xg_t_list[ixg-1]<<"): "<<ch2_FB_d[ixg-1]<<endl; }
  print_ch2<<"### FT u ###"<<endl;
  for(int ixg=1;ixg<num_xg;ixg++) { print_ch2<<"ch2(xg: "<<Bs_xg_t_list[ixg-1]<<"): "<<ch2_FT_u[ixg-1]<<endl; }
  print_ch2<<"### FT d ###"<<endl;
  for(int ixg=1;ixg<num_xg;ixg++) { print_ch2<<"ch2(xg: "<<Bs_xg_t_list[ixg-1]<<"): "<<ch2_FT_d[ixg-1]<<endl; }

  print_ch2.close();
 



  //####################################################################
  //####### PRINT JACKKNIFE INFO ###########

 



  //########################################
  



  
  
  //##########################################################################################
  //##########################################################################################
  
 

  }

  
  rt_fit.Print("../data/ph_emission/"+ph_type+"/"+MESON);
 
  return rt_fit;




}
