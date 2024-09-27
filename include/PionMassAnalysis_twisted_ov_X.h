#ifndef __PionMassAnalysis_twisted_ov_X__
#define __PionMassAnalysis_twisted_ov_X__


#include "LatInfo.h"
#include "random.h"
#include "Bootstrap_fit.h"
#include "input.h"
#include "stat.h"
#include "Corr_analysis.h"
using namespace std;

void Get_plateaux(string Ensemble_tag, string OBS, int &Tmin, int &Tmax);
void Pion_mass_analysis_twisted_ov_X(string CURRENT_TYPE, bool IncludeDisconnected) ;





#endif
