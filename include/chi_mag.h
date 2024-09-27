#ifndef __chi_mag__
#define __chi_mag__


#include "numerics.h"
#include "random.h"
#include "Corr_analysis.h"
#include "stat.h"
#include "Bootstrap_fit.h"
#include "LatInfo.h"
#include "T_min.h"
#include "input.h"
#include "binary_io.h"
#include "Meson_mass_extrapolation.h"

using namespace std;

double Ker_sub(double x);
double susc_pert(double p, double m);
distr_t Ker_sub_distr(const distr_t &x);
void Generate_free_corr_data_VT();
void Compute_free_VT_corr(double am, int Tmax, int L = -1);
double Get_tree_lev_der(double x);
void Compute_magnetic_susc();
void Get_magnetic_susc(bool Include_sea_quark_mass_der);



#endif
