#ifndef __corranalysis__
#define __corranalysis__

#include "numerics.h"
#include "stat.h"
#include "T_min.h"


using namespace std;

class CorrAnalysis {

 public:
  CorrAnalysis() {Perform_Nt_t_average=true; seed = seed_list[0]; UseJack=1; Tmin=0; Tmax=0; Nt=0; Reflection_sign=1;}
  CorrAnalysis(bool UseJack, int Njacks, int nboots) : UseJack(UseJack), Njacks(Njacks), Nboots(nboots) { Perform_Nt_t_average=true; seed=seed_list[0]; Tmin=0; Tmax=0; Nt=0; Reflection_sign=1;}
  CorrAnalysis(bool UseJack, int Njacks, int nboots, int iseed) : UseJack(UseJack), Njacks(Njacks), Nboots(nboots) {
    if(iseed >= (signed)seed_list.size()) crash(to_string((signed)seed_list.size())+" available seeds!, asked for nr. "+to_string(iseed));
    Perform_Nt_t_average=true; seed=seed_list[iseed]; Tmin=0; Tmax=0; Nt=0; Reflection_sign=1;
  } 
  Pfloat Fit_(const distr_t_list& data_t_distr);
  distr_t Fit_distr(const distr_t_list& data_t_distr);
  distr_t_list effective_mass_t(const VVfloat &corr_A, string Obs); //returns jackknife or bootstrap distribution of the effective mass m(t) "acts on data"
  distr_t_list effective_mass_t(const distr_t_list& corr_A_distr, string Obs); //returns jackknife or bootstrap distribution of the effective mass m(t) "acts on distribution"
  distr_t_list effective_slope_t(const VVfloat &corr_A,const VVfloat &corr_B, string Obs); //returns jackknife or bootstrap distribution of the effective slope dm(t) "acts on  data"
  distr_t_list effective_slope_t(const distr_t_list& corr_A_distr,const distr_t_list& corr_B_distr, string Obs);
  distr_t_list effective_slope2_t(const VVfloat &corr_A,const VVfloat &corr_B, string Obs); //returns jackknife or bootstrap distribution of the effective slope squared dm2(t) "acts on  data"
  //distr_t_list effective_slope2_t(const distr_t_list& corr_A_distr, const distr_t_list& corr_B_distr, string Obs); //function overload 
  distr_t_list corr_t(const VVfloat &corr, string Obs); //returns jackknife or bootstrap distribution of a correlator
  distr_t_list residue_t(const VVfloat &corr_A, string Obs); //returns jackknife or boostrap distribution of the residue of the pole "acts on data"
  distr_t_list residue_t(const distr_t_list &corr_A, string Obs); //returns jackknife or boostrap distribution of the residue of the pole "acts on distributions"
  distr_t_list residue_t_wo_m_fit(const VVfloat &corr_A, string Obs); //returns jackknife or boostrap distribution of the residue of the pole wo fitting the effective mass "acts on data" 
  distr_t_list residue_t_wo_m_fit(const distr_t_list &corr_A, string Obs); //returns jackknife or boostrap distribution of the residue of the pole wo fitting the effective mass "acts on distributions"
  distr_t_list decay_constant_t(const VVfloat &corr, string Obs); //returns jackknife or boostrap distribution of a decay constant "acts on data"
  distr_t_list decay_constant_t(const distr_t_list &corr_A, string Obs); //returns jackknife or boostrap distribution of a decay constant "acts on distribution"

  distr_t_list matrix_element_t(const VVfloat &corr_A, string Obs); //returns jackknife or boostrap distribution of the matrix element "acts on data"
  distr_t_list matrix_element_t(const distr_t_list &corr_A, string Obs); //returns jackknife or boostrap distribution of the matrix element "acts on distributions"

  distr_t_list mel_ov_mass_t(const VVfloat &corr_A, string Obs); //returns jackknife or boostrap distribution of the matrix element over the mass "acts on data"
  distr_t_list mel_ov_mass_t(const distr_t_list &corr_A, string Obs); //returns jackknife or boostrap distribution of the matrix element over the mass "acts on distributions"
  
 
  distr_t effective_slope_t_tanh_fit(const distr_t_list& corr_A_distr,const distr_t_list& corr_B_distr, string OBS); //effective slope from a + dm*(T/2-t)*tanh((T/2-t)*meff)
  distr_t_list effective_slope_t_sub_t2(const distr_t_list& corr_A_distr,const distr_t_list& corr_B_distr,const distr_t_list& corr_C_distr,string OBS);
  distr_t_list residue_correction_t(const distr_t_list& corr_A_distr, const distr_t_list& corr_B_distr, string Obs); //returns jackknife or bootstrap distribution of the leading order correction to the residue
  distr_t_list residue_correction_t(const VVfloat& corr_A, const VVfloat& corr_B, string Obs); //function overload
  Vfloat ASymm(const VVfloat& data, int t );
 


  bool UseJack;
  int Njacks;
  int Nboots;
  int Tmin;
  int Tmax;
  int Nt;
  bool Perform_Nt_t_average;
  int Reflection_sign;
  int seed;

private:
  const vector<int> seed_list = {988453453, 8455721, 3314219, 22619374, 44639404, 138121, 23862, 280182, 501740, 40285, 58317, 99138, 2875, 210774, 372194, 37161, 34823, 4743272, 22369, 92720, 328210, 40491, 3291031, 728293, 6942058, 132929, 121231, 3191983 };
  


};
  

























#endif
