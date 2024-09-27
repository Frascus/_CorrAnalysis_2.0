#ifndef __Bs_mumu_gamma__
#define __Bs_mumu_gamma__


#include "numerics.h"
#include "random.h"
#include "Corr_analysis.h"
#include "stat.h"
#include "Bootstrap_fit.h"
#include "LatInfo.h"
#include "pt3_momenta.h"
#include "T_min.h"
#include "input.h"
#include "Num_integrate_l4_decay_rate.h"
#include "virtual_ff_fit.h"
#include "virtual_FF_t_interval_list.h"
#include "ChPT_form_factors.h"
#include "virtual_07.h"
#include "rapidjson/document.h"
#include "rapidjson/stringbuffer.h"
#include "rapidjson/filereadstream.h"
#include "rapidjson/istreamwrapper.h"
#include "rapidjson/writer.h"
#include "rapidjson/ostreamwrapper.h"
#include "rapidjson/error/en.h"
#include "rapidjson/error/error.h"

using namespace std;


//######### LIST OF INPUT PARAMETERS TO COMPUTE DECAY RATE      ################
const double MBs= 5.36692; //GeV
const double tBs = 1.521e-12; // seconds
const double tBs_err= 0.005e-12; //second                              //  
const double hbar =6.582119569e-25; //GeV * seconds
const double Gamma_Bs = hbar / tBs; // total width of the Bs in GeV                                    //
const double FBs = 0.2303;          //
const double FBs_err= 0.0013;   //
const double mb = 4.18;         // GeV mb(mb) in \bar{MS}
const double mb_err = 0.03;                        //  
//const double mc=1.27;//GeV
const double Eg_min= 0.5*MBs*(1- pow(MBs/MBs,2)); //GeV
const double Eg_max= 0.5*MBs*(1- pow(4.2/MBs,2)); //Eg_max
const double xg_min= 2*Eg_min/MBs; 
const double xg_max= 2*Eg_max/MBs; 
const double m_mu= 0.10565837; //GeV
const double x_mu= m_mu/MBs;
const double x_b = mb / MBs;
const double x_b_err = mb_err/MBs;
const double alpha_em = 1/137.04;
const double e2 = alpha_em*4.0*M_PI;
const double GF= 1.1663787*1e-5; //[GeV^-2]
//CKM matrix element
const double Vus=0.2243; const double Vus_err=0.0008;
const double Vcd=0.221;  const double Vcd_err=0.004;
const double Vcs=0.975;  const double Vcs_err=0.006;
const double Vcb=40.8e-3; const double Vcb_err=1.4e-3;
const double Vub=3.82e-3; const double Vub_err=0.20e-3;
const double Vtd=8.6e-3;  const double Vtd_err=0.2e-3;
const double Vts=41.5e-3; const double Vts_err=0.9e-3;
const double Vtb=1.014;   const double Vtb_err=0.029;
const double Vtbs = Vts*Vtb; const double Vtbs_err= sqrt( pow( Vtb*Vts_err,2) + pow(Vts*Vtb_err,2));
//#################################################################################




class rt_FF_Bs {

public:
  //class constructors
  rt_FF_Bs() : UseJack(1), fp(3), phi(3), mp_ov_fp(3), Ch2_fp(3), Ch2_phi(3), Ch2_mp_ov_fp(3), FA(1), FA_u(1), FA_d(1), FV(1), FV_u(1), FV_d(1), FA_T(1), FA_T_u(1), FA_T_d(1), FV_T(1), FV_T_u(1), FV_T_d(1), FB(1), FB_u(1), FB_d(1), FT(1), FT_u(1), FT_d(1),  Nmeas(0), Npars(0), Ndof(0), Use_three_finest(0), Include_a4(0), num_xg(0), Nmeas_K(3), Npars_K(3), Ndof_K(3)  { };
  
  rt_FF_Bs(bool x) : UseJack(x), fp(3), phi(3), mp_ov_fp(3), Ch2_fp(3), Ch2_phi(3), Ch2_mp_ov_fp(3), FA(x), FA_u(x), FA_d(x), FV(x), FV_u(x), FV_d(x), FA_T(x), FA_T_u(x), FA_T_d(x), FV_T(x), FV_T_u(x), FV_T_d(x), FB(1), FB_u(1), FB_d(1), FT(1), FT_u(1), FT_d(1),  Nmeas(0), Npars(0), Ndof(0), Use_three_finest(0), Include_a4(0), num_xg(0) , Nmeas_K(3), Npars_K(3), Ndof_K(3) {    };

  rt_FF_Bs(string P, bool UJ,  bool three_f, bool Inc_a4, int nxg) : fp(3), phi(3), mp_ov_fp(3), Ch2_fp(3), Ch2_phi(3), Ch2_mp_ov_fp(3),  Nmeas_K(3), Npars_K(3), Ndof_K(3)  { Read(P, UJ, three_f, Inc_a4, nxg);}
  
  
  distr_t_list Get_FF(int i) { vector<distr_t_list> A({FA, FA_u, FA_d, FV,FV_u, FV_d, FA_T, FA_T_u, FA_T_d, FV_T, FV_T_u, FV_T_d, FB, FB_u, FB_d, FT, FT_u, FT_d}); return A[i];}
  Vfloat   Get_ch2(int i) { VVfloat A({Ch2_FA, Ch2_FA_u, Ch2_FA_d, Ch2_FV, Ch2_FV_u, Ch2_FV_d,  Ch2_FA_T, Ch2_FA_T_u, Ch2_FA_T_d,  Ch2_FV_T, Ch2_FV_T_u, Ch2_FV_T_d, Ch2_FB, Ch2_FB_u, Ch2_FB_d, Ch2_FT, Ch2_FT_u, Ch2_FT_d}); return A[i];}

  void Print(string path);
  void Read(string path, bool UseJ, bool three_finest, bool inc_a4, int n_xg);
  distr_t Get_FF_2pts(int i) { vector<distr_t> A({mass,fp[0],fp[1], fp[2], phi[0], phi[1], phi[2], mp_ov_fp[0], mp_ov_fp[1], mp_ov_fp[2]}); return A[i];}
  double Get_ch2_2pts(int i) { Vfloat A({Ch2_mass, Ch2_fp[0], Ch2_fp[1], Ch2_fp[2], Ch2_phi[0], Ch2_phi[1], Ch2_phi[2],  Ch2_mp_ov_fp[0], Ch2_mp_ov_fp[1], Ch2_mp_ov_fp[2]}); return A[i];}
  void Fill_FF(int i, const distr_t_list &A, Vfloat &ch2) {
    if(i==0) {FA=A; Ch2_FA=ch2;} 
    else if(i==1) {FA_u=A; Ch2_FA_u=ch2;}
    else if(i==2) {FA_d=A; Ch2_FA_d=ch2;}
    else if(i==3)  {FV=A; Ch2_FV=ch2;}
    else if(i==4) {FV_u=A; Ch2_FV_u=ch2;}
    else if(i==5) {FV_d=A; Ch2_FV_d=ch2;}
    else if(i==6) {FA_T=A; Ch2_FA_T=ch2;}
    else if(i==7) {FA_T_u=A; Ch2_FA_T_u=ch2;}
    else if(i==8) {FA_T_d=A; Ch2_FA_T_d=ch2;}
    else if(i==9) {FV_T=A; Ch2_FV_T=ch2;}
    else if(i==10) {FV_T_u=A; Ch2_FV_T_u=ch2;}
    else if(i==11) {FV_T_d=A; Ch2_FV_T_d=ch2;}
    else if(i==12) {FB=A; Ch2_FB=ch2;}
    else if(i==13) {FB_u=A; Ch2_FB_u=ch2;}
    else if(i==14) {FB_d=A; Ch2_FB_d=ch2;}
    else if(i==15) {FT=A; Ch2_FT=ch2;}
    else if(i==16) {FT_u=A; Ch2_FT_u=ch2;}
    else if(i==17) {FT_d=A; Ch2_FT_d=ch2;}
    else crash("In rt_FF_Bs::Fill_FF(i,A) i must be smaller than 18");
  }

   void Fill_FF_2pt(int i, const distr_t &A, double ch2) {
    if(i==0) {mass=A; Ch2_mass=ch2;} 
    else if(i==1) {fp[0]=A; Ch2_fp[0]=ch2;}
    else if(i==2) {fp[1]=A; Ch2_fp[1]=ch2;}
    else if(i==3) {fp[2]=A; Ch2_fp[2]=ch2;}
    else if(i==4) {phi[0]=A; Ch2_phi[0]=ch2;}
    else if(i==5) {phi[1]=A; Ch2_phi[1]=ch2;}
    else if(i==6) {phi[2]=A; Ch2_phi[2]=ch2;}
    else if(i==7)  {mp_ov_fp[0]=A; Ch2_mp_ov_fp[0]=ch2;}
    else if(i==8)  {mp_ov_fp[1]=A; Ch2_mp_ov_fp[1]=ch2;}
    else if(i==9)  {mp_ov_fp[2]=A; Ch2_mp_ov_fp[2]=ch2;}
    else crash("In rt_FF_Bs::Fill_FF_2pts(i,A) i must be smaller than 10");
  }

  bool UseJack;
  distr_t mass;
  vector<distr_t> fp;
  vector<distr_t> phi;
  vector<distr_t> mp_ov_fp;
  double Ch2_mass;
  Vfloat Ch2_fp;
  Vfloat Ch2_phi;
  Vfloat Ch2_mp_ov_fp;
  distr_t_list FA, FA_u, FA_d, FV, FV_u, FV_d,  FA_T, FA_T_u, FA_T_d, FV_T, FV_T_u, FV_T_d, FB, FB_u, FB_d, FT, FT_u, FT_d;
  Vfloat Ch2_FA, Ch2_FA_u, Ch2_FA_d,  Ch2_FV, Ch2_FV_u, Ch2_FV_d,  Ch2_FA_T, Ch2_FA_T_u, Ch2_FA_T_d,  Ch2_FV_T, Ch2_FV_T_u, Ch2_FV_T_d, Ch2_FB, Ch2_FB_u, Ch2_FB_d, Ch2_FT, Ch2_FT_u, Ch2_FT_d;
  int Nmeas;
  int Npars;
  int Ndof;
  bool Use_three_finest;
  bool Include_a4;
  int num_xg;
  Vfloat Nmeas_K;
  Vfloat Npars_K;
  Vfloat Ndof_K;

  
  
  
};





pair<double, double> WC(int id, double q2, const Vfloat &G1=Vfloat({0,0,0,0,0,0,0,0}), const Vfloat &G2 = Vfloat({0,0,0,0,0,0,0,0}), const Vfloat &BS_F= Vfloat({0,0,0,0,0,0,0,0}), const Vfloat &Gammas_F=Vfloat({0,0,0,0,0,0,0,0}));
void Get_Bs_xg_t_list(int num_xg);
void Get_Bs_xg_to_spline();
void Get_Bs_xg_to_spline_VMD();
void Get_Bs_lattice_spacings_to_print();
void Get_Bs_Tmin_Tmax(string W, int &Tmin, int &Tmax, int ixg, string Ens);
void Get_Bs_Mh_Tmin_Tmax(string obs, string mh,  int &Tmin, int &Tmax, int ixg, string Ens);
void Compute_Bs_mumu_gamma(); 
rt_FF_Bs Get_Bs_mumu_gamma_form_factors(int num_xg, int Perform_continuum_extrapolation, bool Use_three_finest, bool Include_a4, bool UseJack, string Fit_tag, string path_list, string out_tag );




template <typename T1, typename T2, typename T3, typename T4, typename T5,
          typename T6>
double Compute_AFB( double xg, T1&& FV, T2&& FA, T3&& FV_RE_T, T4&& FA_RE_T, T5&& FV_IM_T ,   T6&& FA_IM_T,  double FBS, double xb_val, double Vtbs_val, double tBs_val, Vfloat kappa_list, const Vfloat &Th_list, const Vfloat &W_list, const  Vfloat &Br_list, string MODE, string CH) {

  double diff_br= Compute_Bs_mumugamma_differential_decay_rate(xg, FV, FA, FV_RE_T, FA_RE_T, FV_IM_T , FA_IM_T, FBS, xb_val, Vtbs_val, tBs_val, kappa_list, Th_list, W_list, Br_list,  MODE, CH); 

  double val_FW, err_FW, val_BW, err_BW;
  double prec=5e-5;

  auto FUNC_DOUBLE_DIFF_RATE = [&](double costh) { return Compute_Bs_mumugamma_double_differential_decay_rate( xg, costh, FV, FA, FV_RE_T, FA_RE_T, FV_IM_T, FA_IM_T, FBS, xb_val, Vtbs_val, tBs_val, kappa_list, Th_list, W_list, Br_list,  MODE, CH);};

  gsl_function_pp<decltype(FUNC_DOUBLE_DIFF_RATE)> integrand(FUNC_DOUBLE_DIFF_RATE);
		      gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
		      gsl_function *G = static_cast<gsl_function*>(&integrand);
		      gsl_integration_qags(G, 0, 1.0,  0.0, prec, 1000, w, &val_FW, &err_FW);
		      gsl_integration_qags(G, -1.0, 0,  0.0, prec, 1000, w, &val_BW, &err_BW);
		      gsl_integration_workspace_free(w);
		      if(fabs(err_FW/val_FW) > 2*prec) crash("In Compute_AFB, FW,  cannot reach target precision: "+to_string_with_precision(prec,5));
		      if(fabs(err_BW/val_BW) > 2*prec) crash("In Compute_AFB, BW,  cannot reach target precision: "+to_string_with_precision(prec,5));
		      return (val_FW-val_BW)/diff_br;

  
  

};



template<typename T1, typename T2, typename T3, typename T4, typename T5 , typename T6> 
double Compute_Bs_mumugamma_decay_rate( T1&& FV, T2&& FA, T3&& FV_RE_T, T4&& FA_RE_T, T5&& FV_IM_T ,   T6&& FA_IM_T,  double FBS, double xb_val, double Vtbs_val, double tBs_val, const  Vfloat &kappa_list, const  Vfloat &Th_list, const Vfloat &W_list, const Vfloat &Br_list,  string MODE, string CH,  double xmax=xg_max,  double xmin=xg_min) {

  auto FUNC_DIFF_RATE = [&](double xg) { return Compute_Bs_mumugamma_differential_decay_rate( xg, FV, FA, FV_RE_T, FA_RE_T, FV_IM_T, FA_IM_T, FBS,  xb_val, Vtbs_val, tBs_val, kappa_list, Th_list, W_list, Br_list,  MODE, CH);};

  
  double val, err;
  double prec=5e-5;
  gsl_function_pp<decltype(FUNC_DIFF_RATE)> integrand(FUNC_DIFF_RATE);
		      gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
		      gsl_function *G = static_cast<gsl_function*>(&integrand);
		      gsl_integration_qags(G, xmin, xmax,  0.0, prec, 1000, w, &val, &err);
		      gsl_integration_workspace_free(w);
		      if(fabs(err/val) > 2*prec) crash("In Compute_differential_decay_rate, cannot reach target precision: "+to_string_with_precision(prec,5));

		      return val;


};



template<typename T1, typename T2, typename T3, typename T4, typename T5 , typename T6 > 
double Compute_Bs_mumugamma_differential_decay_rate( double xg,  T1&& FV, T2&& FA, T3&& FV_RE_T, T4&& FA_RE_T, T5&& FV_IM_T , T6&& FA_IM_T,  double FBS, double xb_val, double Vtbs_val, double tBs_val, const Vfloat &kappa_list, const Vfloat &Th_list, const Vfloat &W_list, const Vfloat &Br_list, string MODE, string CH) {

  //integrate at fixed xg

  //integration variable is costheta

  auto FUNC_DOUBLE_DIFF_RATE = [&](double costh) { return Compute_Bs_mumugamma_double_differential_decay_rate(xg, costh, FV, FA, FV_RE_T, FA_RE_T, FV_IM_T, FA_IM_T, FBS,  xb_val, Vtbs_val, tBs_val, kappa_list, Th_list, W_list, Br_list,  MODE, CH);};

  
  double val, err;
  double prec=1e-4;
  gsl_function_pp<decltype(FUNC_DOUBLE_DIFF_RATE)> integrand(FUNC_DOUBLE_DIFF_RATE);
		      gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
		      gsl_function *G = static_cast<gsl_function*>(&integrand);
		      gsl_integration_qags(G, -1.0, 1.0,  0.0, prec, 1000, w, &val, &err);
		      gsl_integration_workspace_free(w);
		      if(fabs(err/val) > 2*prec) crash("In Compute_differential_decay_rate, cannot reach target precision: "+to_string_with_precision(prec,5));

		      return val;
		      
};



template<typename T1, typename T2, typename T3, typename T4, typename T5 , typename T6  > 
double Compute_Bs_mumugamma_double_differential_decay_rate(  double xg, double costh, T1&& FV, T2&& FA, T3&& FV_RE_T, T4&& FA_RE_T, T5&& FV_IM_T , T6&& FA_IM_T, double FBS, double xb_val, double Vtbs_val, double tBs_val, const Vfloat &kappa_list, const Vfloat &Th_list, const Vfloat &W_list, const Vfloat &Br_list,   string MODE, string CH) {


 

  double ix_b = xb_val;
  double iVtbs = Vtbs_val;
  double iFBs = FBS;
  double itBs = tBs_val;
  Vfloat kappa=kappa_list;
  if(CH=="NO_CH") { kappa.clear() ; kappa.resize(8,-1) ;};
    
  //define Mandelstam variables in terms of xg and th
  // s = 1-s 
  double s=1-xg;
  double q2= s*pow(MBs,2);
  double csi= costh*xg*sqrt( 1 - 4*pow(x_mu,2)/s);
  double t= 0.5*(1+2*pow(x_mu,2) - s - csi);
  double u = t+csi;
  //get Wilson coefficients
  double C9_R = WC(9,q2, kappa, Th_list, Br_list,  W_list).first;
  double C9_I = WC(9,q2, kappa, Th_list, Br_list, W_list).second;
  double C9_MOD = sqrt(pow(C9_R,2) + pow(C9_I,2));
  double C7 = WC(7,q2).first;
  double C10 = WC(10,q2).first;
  //define prefactor in decay rate
  double K = pow(iVtbs,2)*pow(GF,2)*(pow(alpha_em,3)*pow(MBs,5)/(pow(2.0,10)*pow(M_PI,4)));
  
  //define F functions
  auto F1 = [&]() { return  (pow(C9_MOD,2) + pow(C10,2))*pow(FV(xg),2) + pow(2*ix_b/s,2)*pow(C7,2)*( pow(FV_RE_T(xg),2) + pow(FV_IM_T(xg),2)) + (4.0*ix_b/s)*FV(xg)*C7*( C9_R*FV_RE_T(xg) +C9_I*FV_IM_T(xg));};

  
  auto F2 = [&]() { return  (pow(C9_MOD,2) + pow(C10,2))*pow(FA(xg),2) + pow(2*ix_b/s,2)*pow(C7,2)*( pow(FA_RE_T(xg),2) + pow(FA_IM_T(xg),2)) + (4.0*ix_b/s)*FA(xg)*C7*( C9_R*FA_RE_T(xg) +C9_I*FA_IM_T(xg));};

  
  //define B functions
  auto B0 = [&]() { return (s+4*pow(x_mu,2))*(F1()+F2()) -8*pow(x_mu,2)*pow(C10,2)*( pow(FV(xg),2) + pow(FA(xg),2)); };
  
  auto B1 = [&]() { return 8*( s*FV(xg)*FA(xg)*C10*C9_R + ix_b*C7*C10*(FV(xg)*FA_RE_T(xg) + FA(xg)*FV_RE_T(xg)));};
  
  auto B2 = [&]() { return s*(F1()+F2());};


  
  //define G function
  auto G1 = [&]() { return K*( pow(xg,2)*B0() + xg*csi*B1() + pow(csi,2)*B2()); };
  auto G2 = [&]() { return K*pow(8*iFBs/MBs,2)*pow(C10*x_mu,2)*( (s+xg*xg/2)/( (u-x_mu*x_mu)*(t-x_mu*x_mu)) - pow( (xg*x_mu)/( (u-x_mu*x_mu)*(t-x_mu*x_mu)),2));};
  auto G12 = [&]() { return -K*(16*iFBs/MBs)*(pow(xg*x_mu,2)/( (u-x_mu*x_mu)*(t-x_mu*x_mu )))*( 2*(xg*ix_b/s)*C10*C7*FV_RE_T(xg) + csi*FA(xg)*pow(C10,2) + xg*FV(xg)*C10*C9_R);};

  //double Ag= 2.0*sqrt(1.0 -xg)/(xg*sqrt( 1 - xg- 4*pow(x_mu,2)));
  //return (1.0/Ag)*(G1()+0.0*G2()+0.0*G12())*itBs/hbar;

  double Ag= 0.5*xg*sqrt( 1 - 4*pow(x_mu,2)/(1-xg));
  if(MODE=="SD") return Ag*G1()*itBs/hbar;
  else if(MODE=="INT") return Ag*G12()*itBs/hbar;
  else if(MODE=="PT") return Ag*G2()*itBs/hbar;
  else  return Ag*(G1()+G2()+G12())*itBs/hbar;
  

};


#endif
