#ifndef __tau_decay_strange__
#define __tau_decay_strange__

#include "numerics.h"
#include "highPrec.h"
#include "Spectral.h"
#include "random.h"
#include "stat.h"
#include "Bootstrap_fit.h"
#include "Corr_analysis.h"
#include "LatInfo.h"
#include "input.h"
#include "Meson_mass_extrapolation.h"
#include "g_minus_2_utilities.h"
#include "binary_io.h"



using namespace std;

double Customized_plateaux_tau_spectre_strange( double alpha, double Emax, string channel, string reg, double s, string Ens );
void tau_decay_analysis_strange();
distr_t Compute_tau_decay_width_strange(bool Is_Emax_Finite, double Emax, double beta);
void get_sigma_list_strange();


#endif
