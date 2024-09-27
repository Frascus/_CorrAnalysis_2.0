#ifndef __Spectral__
#define __Spectral__

#include "numerics.h"
#include "highPrec.h"
#include "random.h"
#include "stat.h"
#include "Bootstrap_fit.h"
#include "Corr_analysis.h"
#include "LatInfo.h"
#include "input.h"
using namespace std;

double Get_exact_gauss(const double &E, const double &m , const double &s, const double &E0);
PrecFloat Get_exact_gauss(const PrecFloat &E,const PrecFloat &m,const PrecFloat &s,const PrecFloat &E0);

PrecFloat BaseFunc(const PrecFloat& E, int t, int T);
PrecFloat aE0(const PrecFloat &E0,const PrecFloat &t, int n);
PrecFloat aE0_std(const PrecFloat &E0,const PrecFloat &t, int n);
PrecFloat aE0_std_Emax(const PrecFloat &E0,const PrecFloat &t, int n);
void Get_Atr(PrecMatr& Atr, const PrecFloat &E0, int T, int tmin, int tmax, const function<PrecFloat( PrecFloat )> &Atr_gen_NORM) ;
void Get_Atr_std(PrecMatr& Atr, const PrecFloat &E0, int T, int tmin, int tmax) ;
void Get_Atr_std_Emax(PrecMatr& Atr, const PrecFloat &E0, int T, int tmin, int tmax) ;
void Get_Rt(PrecVect& Rt, const PrecFloat &E0,   int T, int tmin, int tmax);
void Get_ft(PrecVect& ft, const PrecFloat &E0, const PrecFloat &m, const PrecFloat &s, int jack_id, int  T, int tmin, int tmax, string SMEARING_FUNC,  const function<PrecFloat(const PrecFloat&, const PrecFloat&,const PrecFloat&,const PrecFloat&, int)> &f, const function<PrecFloat(const PrecFloat &, const PrecFloat &, const PrecFloat &, const PrecFloat &, int)> &F_NORM) ;
void Get_ft_std(PrecVect& ft, const PrecFloat &E0, const PrecFloat &m, const PrecFloat &s, int jack_id, int  T, int tmin, int tmax, string SMEARING_FUNC,  const function<PrecFloat(const PrecFloat&, const PrecFloat&,const PrecFloat&,const PrecFloat&, int)> &f) ;
void Get_ft_std_Emax(PrecVect& ft, const PrecFloat &E0, const PrecFloat &m, const PrecFloat &s, int jack_id, int  T, int tmin, int tmax, string SMEARING_FUNC,  const function<PrecFloat(const PrecFloat&, const PrecFloat&,const PrecFloat&,const PrecFloat&, int)> &f) ;
void Get_bt(PrecVect& bt,const PrecFloat &E,  int T, int tmin, int tmax);
PrecFloat Get_norm_constraint(PrecFloat &m, PrecFloat &s, PrecFloat &E0, int jack_id, string SMEARING_FUNC, const function<PrecFloat(const PrecFloat&,const PrecFloat&,const PrecFloat&,const PrecFloat&, int)> &f);
PrecFloat Get_M2(PrecFloat &m, PrecFloat &s, PrecFloat &E0, int jack_id, const function<PrecFloat(const PrecFloat&,const PrecFloat&,const PrecFloat&,const PrecFloat&, int )> &f, const function<PrecFloat(const PrecFloat &, const PrecFloat &, const PrecFloat &, const PrecFloat &, int)> &F_NORM);
PrecFloat Get_M2_std_norm(PrecFloat &m, PrecFloat &s, PrecFloat &E0, int jack_id, const function<PrecFloat(const PrecFloat&,const PrecFloat&,const PrecFloat&,const PrecFloat&, int )> &f);
PrecFloat Get_M2_std_norm_Emax(PrecFloat &m, PrecFloat &s, PrecFloat &E0, int jack_id, const function<PrecFloat(const PrecFloat&,const PrecFloat&,const PrecFloat&,const PrecFloat&, int )> &f);
void Get_Rt_up_to_N(PrecFloat &E0, int T, int tmin, int tmax, vector<PrecVect> &Rt_n);
void Get_G_matrix(PrecMatr &G,const PrecMatr &Atr_inv, vector<PrecVect> &Rt_n);
void Get_M_N(PrecFloat &m, PrecFloat &s, PrecFloat &E0, int jack_id,  const function<PrecFloat(const PrecFloat&, const PrecFloat&,const PrecFloat&,const PrecFloat&, int)> &f, PrecVect &M_n);
void Get_M_tilde_N( const PrecVect &ft, const PrecMatr &Atr_inv, vector<PrecVect> &Rt_n,  PrecVect &M_tilde_n);
void Compute_covariance_matrix(PrecMatr &B,const PrecMatr &Atr, const distr_t_list &corr, int tmin, int tmax, PrecFloat m, PrecFloat s, string spec_type, string MODE, Vfloat &covariance);
void Get_optimal_lambda(const PrecMatr &Atr,const PrecMatr &Atr_std_norm, const PrecMatr &Atr_std_norm_Emax, const PrecMatr &B,const PrecVect &ft,vector<PrecVect>& ft_jack,  const PrecVect &ft_std_norm, const PrecVect &ft_std_norm_Emax, const PrecFloat & M2, const PrecFloat &M2_std , const PrecFloat &M2_std_Emax, const double &mean, const double &sigma, const double &Estart,  double& lambda_opt, vector<PrecVect> Rt_n, const PrecVect &M_n ,vector<PrecVect>& M_n_jack, const distr_t_list & corr,int T, int tmin, int tmax, const double mult,   string MODE, string curr_type, string SMEARING_FUNC, string CORR_NAME, double Ag_ov_A0_tg, bool JackOnKer, const distr_t &Prefact, const distr_t &offset,  string analysis_name, const function<PrecFloat(const PrecFloat&, const PrecFloat&,const PrecFloat&,const PrecFloat&, int)> &f , const function<double(const function<double(double)>&)> &syst_func, bool Use_guess_density, const function<double(double)> &guess_density);
distr_t Get_Laplace_transfo(double mean, double sigma, double Estart, int T, int tmax, int prec, string SMEARING_FUNC, const function<PrecFloat(const PrecFloat&, const PrecFloat&,const PrecFloat&,const PrecFloat&, int)> &f, const distr_t_list& corr, double &syst,const double mult, double& pull, string MODE, string reg_type, string CORR_NAME, double Ag_ov_A0_target, bool JackOnKer, const distr_t &Prefact,const  double &offset,  string analysis_name, Vfloat &covariance, const function<double(const function<double(double)>&)> &syst_func, bool Use_guess_density, const function<double(double)> &guess_density,  bool Int_up_to_Max=false, double Max_Erg=4.0, double b=2.99, bool ONLY_FW=0, bool GENERALIZED_NORM=0, const function<PrecFloat(const PrecFloat &, const PrecFloat &, const PrecFloat &, const PrecFloat &, int)> F_NORM = [](const PrecFloat E, const PrecFloat m, const PrecFloat s, const PrecFloat E0, int jack_id) { return 1.0;}  , const function<PrecFloat( PrecFloat )> Atr_gen_NORM= [](PrecFloat t) {return 0.0;} );

distr_t Get_Laplace_transfo_tmin(double mean, double sigma, double Estart, int T, int tmin, int tmax, int prec, string SMEARING_FUNC, const function<PrecFloat(const PrecFloat&, const PrecFloat&,const PrecFloat&,const PrecFloat&, int)> &f, const distr_t_list& corr, double &syst,const double mult, double& pull, string MODE, string reg_type, string CORR_NAME, double Ag_ov_A0_target, bool JackOnKer, const distr_t &Prefact,const  distr_t &offset,  string analysis_name, Vfloat &covariance, const function<double(const function<double(double)>&)> &syst_func, bool Use_guess_density, const function<double(double)> &guess_density,  bool Int_up_to_Max=false, double Max_Erg=4.0, double b=2.99, bool ONLY_FW=0, bool GENERALIZED_NORM=0, const function<PrecFloat(const PrecFloat &, const PrecFloat &, const PrecFloat &, const PrecFloat &, int)> F_NORM = [](const PrecFloat E, const PrecFloat m, const PrecFloat s, const PrecFloat E0, int jack_id) { return 1.0;}  , const function<PrecFloat( PrecFloat )> Atr_gen_NORM= [](PrecFloat t) {return 0.0;} );


#endif
