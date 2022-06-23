#ifndef __init__
#define __init__


#include "3pts_mes_gamma_W.h"
#include "PionMassAnalysis.h"
#include "Axion_l7.h"
#include "PionMassAnalysis_twisted.h"
#include "PionMassAnalysis_twisted_adim.h"
#include "PionMassAnalysis_twisted_ov_X.h"
#include "PionMassAnalysis_ov_X.h"
#include "gm2.h"
#include "Spectral_tests.h"
#include "Spectral_test_Vittorio.h"
#include "vph_Nissa.h"
#include "tau_decay.h"

using namespace std;

class MasterClass_analysis {

 public:
  MasterClass_analysis(string Input);


 private:
  void Analysis_manager();
  string Analysis_Mode;
  string Meson_to_analyze;
  bool IncludeDisconnected;
  string CURRENT_TYPE;


} ;




#endif
