#ifndef __3pts_meson_gamma_W__
#define __3pts_meson_gamma_W__




#include "numerics.h"
#include "random.h"
#include "Corr_analysis.h"
#include "stat.h"
#include "Bootstrap_fit.h"
#include "header_file_virph.h"
#include "LatInfo.h"
#include "pt3_momenta.h"
#include "T_min.h"
#include "input.h"
#include "VMD_Regge_charmed.h"
#include "Num_integrate_l4_decay_rate.h"
#include "virtual_ff_fit.h"
#include "virtual_FF_t_interval_list.h"
#include "ChPT_form_factors.h"








void Get_Tmin_Tmax(string corr_type, string Ens_tag, CorrAnalysis &corr, double xg, string W);

distr_t_list V_ave_unpolarized(vector<vector<distr_t_list>>& distr_mom_k, vector<vector<distr_t_list>>& distr_mom_0, distr_t& Meson_mass, pt3_momenta& mom, int twall); 

distr_t_list A_ave_unpolarized(vector<vector<distr_t_list>>& distr_mom_k);


distr_t_list H_2(const distr_t_list& H30,const distr_t_list& H11,const distr_t_list& H33,const distr_t_list& H30_0,const distr_t_list& H11_0,const distr_t_list& H33_0,pt3_momenta& Mom,const distr_t& m); //valid for p=0, k= kz
distr_t_list H_1(const distr_t_list& H30,const distr_t_list& H11,const distr_t_list& H33,const distr_t_list& H30_0,const distr_t_list& H11_0,const distr_t_list& H33_0,pt3_momenta& Mom,const distr_t& m); //valid for p=0, k=kz
distr_t_list FA_off(const distr_t_list& H30,const distr_t_list& H11,const distr_t_list& H33,const distr_t_list& H30_0,const distr_t_list& H11_0,const distr_t_list& H33_0,pt3_momenta& Mom,const distr_t& m); //valid for p=0, k=kz

distr_t_list H_1_impr(const distr_t_list & H30, const distr_t_list& H03, const distr_t_list& H11,const distr_t_list& H33,const distr_t_list& H11_0,const distr_t_list& H33_0,pt3_momenta& Mom,const distr_t& m); //valid for p=0, k= kz
distr_t_list H_2_impr(const distr_t_list & H30, const distr_t_list& H03, const distr_t_list& H11,const distr_t_list& H33,const distr_t_list& H11_0,const distr_t_list& H33_0,pt3_momenta& Mom,const distr_t& m); //valid for p=0, k= k
distr_t_list FA_off_impr(const distr_t_list & H30, const distr_t_list& H03, const distr_t_list& H11,const distr_t_list& H33,const distr_t_list& H11_0,const distr_t_list& H33_0,pt3_momenta& Mom,const distr_t& m); //valid for p=0, k= k

distr_t_list H_1_mixed_diag(const distr_t_list & H30, const distr_t_list& H03, const distr_t_list& H11,const distr_t_list& H33,const distr_t_list& H11_0,const distr_t_list& H33_0, const distr_t_list& H33_kz0 ,const distr_t_list& H11_kz0, pt3_momenta& Mom,const distr_t& m); //valid for p=0, k= kz
distr_t_list H_2_mixed_diag(const distr_t_list & H30, const distr_t_list& H03, const distr_t_list& H11,const distr_t_list& H33,const distr_t_list& H11_0,const distr_t_list& H33_0, const distr_t_list& H33_kz0, const distr_t_list& H11_kz0,   pt3_momenta& Mom,const distr_t& m); //valid for p=0, k= k
distr_t_list FA_off_mixed_diag(const distr_t_list & H30, const distr_t_list& H03, const distr_t_list& H11,const distr_t_list& H33,const distr_t_list& H11_0,const distr_t_list& H33_0, const distr_t_list& H33_kz0, const distr_t_list& H11_kz0, pt3_momenta& Mom,const distr_t& m); //valid for p=0, k= k


void Compute_form_factors();



void Plot_form_factors(string W, distr_t_list& F, distr_t &F_fit, int Tmin, int Tmax, int Nt, string Ens_tag, double xg, double offsh, int smearing);

void Fit_contaminations(string W, distr_t_list& F, distr_t &F_fit, string Ens_tag, int Nt, double xg, double offsh, int smearing);

void Fit_contaminations_from_derivative(string W, distr_t_list& F, distr_t &F_fit, string Ens_tag, int Nt, double xg, double offsh, int smearing);


#endif
