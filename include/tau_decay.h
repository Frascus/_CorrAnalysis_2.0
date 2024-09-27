#ifndef __tau_decay__
#define __tau_decay__

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

double Customized_plateaux_tau_spectre( double alpha, double Emax, string channel, string reg, double s, string Ens );
void tau_decay_analysis();
void Compute_tau_decay_width(bool Is_Emax_Finite, double Emax, double beta,LL_functions &LL);
void get_sigma_list();


#endif
