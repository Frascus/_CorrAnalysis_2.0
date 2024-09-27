#ifndef __HVP__
#define __HVP__


#include "g_minus_2_utilities.h"
#include "gm2_fits.h"
#include "Meson_mass_extrapolation.h"
#include "Bootstrap_fit.h"
#include "scale_setting.h"
#include "gm2.h"
#include "stat.h"
#include "binary_io.h"




void Bounding_HVP(distr_t &HVP, int &Tcut_opt, const distr_t_list &V, const distr_t &a, string path,  distr_t lowest_mass);
void HVP();


#endif
