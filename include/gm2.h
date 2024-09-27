#ifndef __gm2__
#define __gm2__


#include "g_minus_2_utilities.h"
#include "gm2_fits.h"
#include "Meson_mass_extrapolation.h"
#include "Bootstrap_fit.h"
#include "scale_setting.h"




using namespace std;


distr_t PI_q2( const distr_t_list &V,const distr_t &a, double Q2, int Tmax);
distr_t PI_q2( const distr_t_list &V,const distr_t &a, const distr_t &Q2, int Tmax);
distr_t PI_q2_fixed_t(const distr_t &Vt, const distr_t &a, const distr_t &Q2, int t);
void Add_ens_val_PI_q2( vector<distr_t_list> &PI_master, const distr_t_list &V, const distr_t &a, int Tmax);
void Add_ens_val_PI_q2( vector<distr_t_list> &PI_master, const distr_t_list &PI_per_ens);
void Get_PI_q2( distr_t_list &PI_per_ens, const distr_t_list &V, const distr_t &a, int Tmax);
void Get_PI_q2(distr_t_list &PI_per_ens, const Vfloat &V, const distr_t &a, int Tmax);
void Bounding_amu_W(distr_t &amu_W, const distr_t_list &V, const distr_t &a, string path,const distr_t &Z, distr_t lowest_mass);
void Bounding_PI_q2(distr_t_list &PI_per_ens, const distr_t_list &V, const distr_t &a, string path, Vint &Tdatas_opt, distr_t &lowest_mass);
void Bounding_PI_q2_disco(distr_t_list &PI_per_ens, const distr_t_list &V, const distr_t_list &Conn_guess, const distr_t &a, string path, Vint &Tdatas_opt, distr_t m_rho_GS);
void Get_amu_W_eps( vector<distr_t_list> &eps_win_list, const distr_t_list &V, const distr_t &a);
double Get_GS_artifacts_win(double M, double L, LL_functions &LL, bool mode_mixed);
void Gm2() ;





#endif
