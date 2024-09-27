#ifndef __vph_Nissa__
#define __vph_Nissa__


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

using namespace std;


//#LIST OF INPUT PARAMETER TO BE USED IN THE CALCULATION OF THE DECAY RATE

const double MDs= 1.96835; //GeV 




class rt_FF {

public:
  rt_FF() : UseJack(1), FA(1), FV(1), FA_u(1), FV_u(1), FA_d(1), FV_d(1), Nmeas(0), Npars(0), Ndof(0), Use_three_finest(0), Include_a4(0), num_xg(0) { };
  rt_FF(bool x) : UseJack(x), FA(x), FV(x), FA_u(x), FV_u(x), FA_d(x), FV_d(x), Nmeas(0), Npars(0), Ndof(0), Use_three_finest(0), Include_a4(0), num_xg(0) {    };
  distr_t_list Get_FF(int i) { vector<distr_t_list> A({FA,FV,FA_u,FV_u,FA_d,FV_d}); return A[i];}
  Vfloat   Get_ch2(int i) { VVfloat A({Ch2_FA, Ch2_FV, Ch2_FA_u, Ch2_FV_u, Ch2_FA_d, Ch2_FV_d}); return A[i];}

  bool UseJack;
  distr_t_list FA, FV, FA_u, FV_u, FA_d, FV_d;
  Vfloat Ch2_FA, Ch2_FV, Ch2_FA_u, Ch2_FV_u, Ch2_FA_d, Ch2_FV_d;
  int Nmeas;
  int Npars;
  int Ndof;
  bool Use_three_finest;
  bool Include_a4;
  int num_xg;
  
  
};

void Get_xg_t_list(int num_xg);
void Get_xg_to_spline();
void Get_xg_to_spline_VMD();
void Get_lattice_spacings_to_print();
void Get_Tmin_Tmax(string W, int &Tmin, int &Tmax, int ixg, string Ens);
void Compute_form_factors_Nissa();
rt_FF Get_form_factors_Nissa(int num_xg, int Perform_continuum_extrapolation, bool Use_three_finest, bool Include_a4, bool UseJack, string Fit_tag, string path_list );


template<typename T1, typename T2> 
double Compute_Ds_lnugamma_decay_rate(double rl, string MODE,int jk, int Nj,  T1&& FV, T2&& FA, double fp,  double xmin,  double xmax=1.0   ) {

  auto FUNC_DIFF_RATE = [&](double xg) { return Compute_Ds_lnugamma_differential_decay_rate(rl, xg,jk, Nj, FV, FA, fp, MODE);};

  
  double val, err;
  double prec=1e-7;
  gsl_function_pp<decltype(FUNC_DIFF_RATE)> integrand(FUNC_DIFF_RATE);
		      gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
		      gsl_function *G = static_cast<gsl_function*>(&integrand);
		      gsl_integration_qags(G, xmin, xmax,  0.0, prec, 1000, w, &val, &err);
		      gsl_integration_workspace_free(w);
		      if(fabs(err/val) > 5*prec) crash("In Compute_differential_decay_rate, cannot reach target precision: "+to_string_with_precision(prec,5));

		      return val;
};


template<typename T1, typename T2> 
double Compute_Ds_lnugamma_differential_decay_rate(double rl, double xg, int jk, int Nj,  T1&& FV, T2&& FA, double fp, string MODE) {

  if(xg>= 1.0 - rl*rl) return 0.0;

  double ran=0;
  GaussianMersenne G_Ds(43539998);
  for(int i=0;i<Nj;i++) {double  rn= G_Ds(); if(i==jk) ran=rn/sqrt(Nj-1.0);} 
  

  double GF=  1.1663787*1e-5; //GeV^-2
  double hbar = 6.582119569e-25 ; //GeV * s
  double tDs= 5.04e-13 + ran*0.04e-13 ; //s
  double Gamma_Ds= hbar/tDs ;
  double alpha_em= 1/137.04 ;
 
  //define prefactor in decay rate
  double K = pow(GF,2)*(pow(MDs,3))*alpha_em*rl*rl*pow(1-(rl*rl),2)/(M_PI*M_PI*32*Gamma_Ds);

  double F_pt= -(1.0/xg)*(  ( (pow(2-xg,2)/(1.0-xg)) -4*rl*rl)*(1-xg-rl*rl) -(2*(1-rl*rl)*(1+rl*rl-xg)+xg*xg)*log( (1-xg)/(rl*rl) ) )*(2.0/pow(1- (rl*rl),2));
  
  double F_SD= pow(xg,3)*((2+rl*rl-2*xg)*pow(1-xg-rl*rl,2)/(6.0*pow(1-xg,2)))*pow(MDs,2)/(2.0*rl*rl*pow(1-rl*rl,2));
  
  double F_plus= (xg/2.0)*( pow(rl,4)/(1.0-xg) -1.0 + xg + 2*rl*rl*log( (1-xg)/(rl*rl)))*(-2*MDs/pow(1-rl*rl,2));

  double F_minus = -F_plus + xg*xg*( rl*rl/(1.0-xg) -1.0 + log( (1-xg)/(rl*rl)))*(-2*MDs/pow(1-rl*rl,2));

  double PT= F_pt*fp*fp;

  double INT = F_plus*(FV(xg) +FA(xg))*fp + F_minus*(FV(xg) - FA(xg))*fp;

  double SD = F_SD*( pow(FV(xg)+FA(xg),2) + pow(FV(xg)-FA(xg),2) );
			
  if(MODE=="PT") return K*PT;
  else if(MODE=="INT") return K*INT ;
  else if(MODE=="SD") return K*SD;
  
  return K*(PT+INT+SD);
    

};


#endif
