#ifndef __virtual_07__
#define __virtual_07__


#include "numerics.h"
#include "random.h"
#include "Corr_analysis.h"
#include "stat.h"
#include "Bootstrap_fit.h"
#include "LatInfo.h"
#include "pt3_momenta.h"
#include "input.h"
#include "highPrec.h"
#include "Spectral.h"
#include "kernel.h"



class rt_07_Bs {

public:
  rt_07_Bs() : UseJack(1),  F_T_u(1), F_T_d_RE(1), F_T_d_IM(1),  num_xg(0), y_eff(0) { MESON="" ; };
  rt_07_Bs(bool x) : UseJack(x),  F_T_u(x), F_T_d_RE(x), F_T_d_IM(x),  num_xg(0), y_eff(0) { MESON="";   };
  distr_t_list Get_FF(int i) { vector<distr_t_list> A({ F_T_u, F_T_d_RE, F_T_d_IM}); return A[i];}
  
  

  bool UseJack;
  distr_t_list F_T_u, F_T_d_RE, F_T_d_IM;
  int num_xg;
  Vfloat y_eff;
  string MESON;
};


rt_07_Bs Get_virtual_tensor_FF(int num_xg, bool UseJack, int Njacks, string MESON,  string Corr_path,  string path_out);


#endif
