#ifndef __scale_setting__
#define __scale_setting__

#include "Meson_mass_extrapolation.h"
#include "Bootstrap_fit.h"
#include "numerics.h"
#include "stat.h"
#include "Corr_analysis.h"
#include "LatInfo.h"

using namespace std;


void Determine_scale_from_w0_fp(const distr_t_list &Mpi, const distr_t_list &fpi, const distr_t_list &w0, const vector<double> &L_list, const vector<string>  &Ensemble_tags,distr_t &a_distr_A, distr_t &a_distr_B, distr_t &a_distr_C, distr_t &a_distr_D, bool UseJack, bool Use_three_finest);

void Determine_scale_from_fp(const distr_t_list &Mpi_scale_setting_phys_point_ens, const distr_t_list &fpi_scale_setting_phys_point_ens, const distr_t_list &Mpi_scale_setting_Bens, const distr_t_list &fpi_scale_setting_Bens, const distr_t_list &Mpi_scale_setting_Aens, const distr_t_list &fpi_scale_setting_Aens,  const vector<double> &L_phys_point, const vector<double> &L_B_ens, const vector<double> &L_A_ens, const vector<string>  &Ensemble_phys_point_tag_list, const vector<string>  &Ensemble_B_tag_list, const vector<string>  &Ensemble_A_tag_list, distr_t &a_from_fp_A, distr_t &a_from_fp_B, distr_t &a_from_fp_C, distr_t &a_from_fp_D,bool UseJack,bool Use_three_finest_in_scale_setting_fp);



void Determine_scale_from_fp_FLAG(const distr_t_list &Mpi_scale_setting_phys_point_ens, const distr_t_list &fpi_scale_setting_phys_point_ens, const distr_t_list &Mpi_scale_setting_Bens, const distr_t_list &fpi_scale_setting_Bens, const distr_t_list &Mpi_scale_setting_Aens, const distr_t_list &fpi_scale_setting_Aens,  const vector<double> &L_phys_point, const vector<double> &L_B_ens, const vector<double> &L_A_ens, const vector<string>  &Ensemble_phys_point_tag_list, const vector<string>  &Ensemble_B_tag_list, const vector<string>  &Ensemble_A_tag_list, distr_t &a_from_fp_A, distr_t &a_from_fp_B, distr_t &a_from_fp_C, distr_t &a_from_fp_D, distr_t &a_from_fp_E, bool UseJack,bool Use_three_finest_in_scale_setting_fp);

#endif
