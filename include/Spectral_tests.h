#ifndef __spectral_tests__
#define __spectral_tests__

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


using namespace std;

void Spectral_tests();
void Get_exp_smeared_R_ratio(const Vfloat &Ergs_GeV_list, double sigma);


#endif
