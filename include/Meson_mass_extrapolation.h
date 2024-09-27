#ifndef __Meson_mass_extrapolation__
#define __Meson_mass_extrapolation__

#include "Bootstrap_fit.h"
#include "numerics.h"
#include "stat.h"
#include "Corr_analysis.h"

using namespace std;

distr_t Obs_extrapolation_meson_mass( vector<distr_t> &Meas, vector<distr_t> &Op, double Op_phys, string Print_dir, string Tag, bool UseJack, string EXTRAPOLATION_MODE);
distr_t Obs_extrapolation_meson_mass( vector<distr_t> &Meas, vector<distr_t> &Op,const distr_t &Op_phys, string Print_dir, string Tag, bool UseJack, string EXTRAPOLATION_MODE);

















#endif
