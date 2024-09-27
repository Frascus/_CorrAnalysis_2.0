#ifndef __gm2_fits__
#define __gm2_fits__

#include "g_minus_2_utilities.h"
#include "Meson_mass_extrapolation.h"
#include "Bootstrap_fit.h"
#include "numerics.h"
#include "stat.h"
#include "Corr_analysis.h"
#include "LatInfo.h"

using namespace std;



void Perform_Akaike_fits(const distr_t_list &meas_tm,const distr_t_list &meas_OS, const distr_t &a_A, const distr_t &a_B, const distr_t &a_C, const distr_t &a_D, Vfloat &L_list,const distr_t_list &a_distr_list, const distr_t_list &Mpi_fit, const distr_t_list &fp_fit, vector<string> &Ens_Tag, bool UseJack, int Njacks, int Nboots,  string W_type, string channel,vector<string> &Inc_a2_list, vector<string> &Inc_FSEs_list, vector<string> &Inc_a4_list, vector<string> &Inc_mass_extr_list, vector<string> &Inc_single_fit_list, VPfloat &Inc_log_corr_list, bool allow_a4_and_log, bool allow_only_finest, int w0s_mult, bool Is_log_fitted, int vol_mult, int mass_mult, double cont_guess, LL_functions &LL, double tmin_SD, distr_t &return_distr, distr_t &return4_distr );


void Perform_Akaike_fits_PI_Q2(double Q2, const distr_t_list &meas_tm,const distr_t_list &meas_OS, const distr_t &a_A, const distr_t &a_B, const distr_t &a_C, const distr_t &a_D, Vfloat &L_list,const distr_t_list &a_distr_list, const distr_t_list &Mpi_fit, const distr_t_list &fp_fit, vector<string> &Ens_Tag, bool UseJack, int Njacks, int Nboots,  string W_type, string channel,vector<string> &Inc_a2_list, vector<string> &Inc_FSEs_list, vector<string> &Inc_a4_list, vector<string> &Inc_mass_extr_list, vector<string> &Inc_single_fit_list, VPfloat &Inc_log_corr_list, bool allow_a4_and_log, bool allow_only_finest, int w0s_mult, bool Is_log_fitted, int vol_mult, int mass_mult, double cont_guess, LL_functions &LL, distr_t &return_distr );







#endif
