#ifndef __Num_integrate_l4_decay_rate__
#define __Num_integrate_l4_decay_rate__



#include "numerics.h"
#include "random.h"
#include "Corr_analysis.h"
#include "stat.h"
#include "Bootstrap_fit.h"
#include "T_min.h"
#include "input.h"
#include "kernel.h"
#include "ChPT_form_factors.h"


using namespace std;




class Decay_Rate_Integration_Result {

 public:
  Decay_Rate_Integration_Result() {}
  Decay_Rate_Integration_Result(int UseJack) : Jack_Distr_Int_Quad_mumu(UseJack), Jack_Distr_Int_Quad_ee(UseJack), Jack_Distr_Int_MonteCarlo_mumu(UseJack), Jack_Distr_Int_MonteCarlo_ee(UseJack), Jack_Distr_Int_MonteCarlo_mumumu(UseJack), Jack_Distr_Int_MonteCarlo_eee(UseJack) {}
  

  double Int_Quad_val_mumu, Int_Quad_val_ee, Int_MonteCarlo_val_mumu, Int_MonteCarlo_val_ee;

  double Int_Quad_err_mumu, Int_Quad_err_ee, Int_MonteCarlo_err_mumu, Int_MonteCarlo_err_ee;

  double Int_MonteCarlo_val_mumumu, Int_MonteCarlo_val_eee;

  double Int_MonteCarlo_err_mumumu, Int_MonteCarlo_err_eee;

  //int Nsteps_MonteCarlo_mumu, Nsteps_MonteCarlo_ee;

  //int Nfunc_eval_Quad_mumu, Nfunc_eval_Quad_elel;

  double eps_rel_mumu, eps_rel_ee;

  Vfloat Stat_err_MonteCarlo_mumu, Stat_err_MonteCarlo_ee, Stat_err_MonteCarlo_mumumu, Stat_err_MonteCarlo_eee;

  bool Exit_status;

  distr_t Jack_Distr_Int_Quad_mumu, Jack_Distr_Int_Quad_ee, Jack_Distr_Int_MonteCarlo_mumu, Jack_Distr_Int_MonteCarlo_ee;

  distr_t Jack_Distr_Int_MonteCarlo_mumumu, Jack_Distr_Int_MonteCarlo_eee;

};

void display_results (char *title, double result, double error);

double MonteCarlo_integration_extended_phase_space(const function<double(double,double)> &H1, const function<double(double,double)> &H2, const function<double(double,double)> &FA, const function<double(double,double)> &FV, double mk, double fk,string channel, double xk_inf,  bool same_lepton, int MySeed, string mode);

Decay_Rate_Integration_Result Num_Integrate_Decay_Rate(vector<function<double(double, double)>> &H1, vector<function<double(double, double)>> &H2  , vector<function<double(double, double)>> &FA, vector<function<double(double, double)>> &FV, distr_t &m_distr, distr_t & fp_distr, bool UseJack, string Meson, bool Print);






#endif
