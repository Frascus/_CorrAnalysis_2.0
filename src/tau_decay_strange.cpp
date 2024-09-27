#include "../include/tau_decay_strange.h"

// Use (void) to silence unused warnings.
#define assertm(exp, msg) assert(((void)msg, exp))

const double MPiPhys=0.135;
const double alpha = 1.0/137.04;
const bool UseJack=1;
const int Njacks=100;
const int Nboots=5000;
const double ln2_10=3.32192809489;
const double fm_to_inv_Gev= 1.0/0.197327;
const int prec = 128;
const double Nc=3;
bool tau_strange_verbosity_lev=1;
const double GF= 1.1663787*1e-5; //[GeV^-2]
//CKM matrix elements
const double m_tau = 1.77686;
const double m_kappa = 0.494600000; //old is 0.4942
const double m_kappa_err = 0.0004e-8;
const double m_etas = 0.68989;
const double m_etas_err= 0.00050;
const double E0_l = 0.9*m_kappa;
const double E0_sp = 0.9*m_kappa; // 0.9*(m_kappa+MPiPhys); //0.9 * (m_kappa + MPiPhys); //E0_l
double E0_A_sp = 0.9*m_kappa; //0.9*(m_kappa+ 2*MPiPhys); // 0.9*(m_kappa+2*MPiPhys);
const double Rs_HFLAV = 0.163260;
const double D_Rs_HFLAV= 0.0027;
Vfloat sigma_list_strange;
const double C_V = 2*M_PI/(pow(m_tau,3));
const double GAMMA_FACT= 12*M_PI; //12*M_PI*pow(Vud*GF,2);
const string MODE="TANT";
bool Use_t_up_to_T_half_strange=false;
const int sm_func_mode= 0;
const string SM_TYPE_0= "KL_"+to_string(sm_func_mode);
const string SM_TYPE_1= "KT_"+to_string(sm_func_mode);
VVfloat covariance_fake_strange;
const double QCD_scale= 0.3*fm_to_inv_Gev;
bool Skip_spectral_density_analysis_strange=false;
const bool Perform_continuum_extrapolation=false;
bool Use_Customized_plateaux_strange=true;
using namespace std;


double Customized_plateaux_tau_spectre_strange( double alpha, double Emax, string channel, string reg, double s, string Ens ) {

  //if(channel=="T") channel="Aii";
  //if(channel=="L") channel="A0";
 
  double Ra0=-1;
  int alpha_m= (int)(alpha+1);

  if(alpha_m < 3) crash("Customized_plateaux_tau_spectre called with alpha = "+to_string_with_precision(alpha,2));
  //if( (alpha_m != 3) && (alpha_m != 4) && (alpha_m != 5) ) crash("Customized_plateaux_tau_spectre called with (int)alpha: "+to_string(alpha_m));

  if(Ens=="cZ211a.077.64") Ens = "cB211b.072.64";

  //if(Ens=="cE211a.044.112") Ens="cD211a.054.96";
  //if(Ens=="cB211b.072.96") Ens="cB211b.072.64";

  if(alpha_m > 3) {

    if(s< 0.125) {
  
      if( reg=="tm") {
	if(channel=="Aii") {
	  if(Ens == "cB211b.072.64") {  Ra0 =3e6;  }
	  else if(Ens == "cB211b.072.96") { Ra0= 4e8;   }
	  else if(Ens == "cC211a.06.80") { Ra0= 6e6;   } 
	  else if(Ens == "cC211a.06.112") { Ra0= 6e6;   }
	  else if(Ens == "cD211a.054.96") {  Ra0=8e6;  }
	  else if(Ens == "cE211a.044.112") {  Ra0=4e6;  }
	  else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
	}
	else if(channel=="Vii") {
	  if(Ens == "cB211b.072.64") { Ra0=8e5;   }
	  else if(Ens == "cB211b.072.96") { Ra0= 2.5e7;    } 
	  else if(Ens == "cC211a.06.80") {  Ra0=3e7;  } 
	  else if(Ens == "cC211a.06.112") { Ra0= 2e6;   }
	  else if(Ens == "cD211a.054.96") {  Ra0= 2e7;   }
	  else if(Ens == "cE211a.044.112") {  Ra0=3e7;  }
	  else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
	}
	else if(channel=="A0") {
	  if(Ens == "cB211b.072.64") { Ra0=1e7;   }
	  else if(Ens == "cB211b.072.96") { Ra0=1e7;   }
	  else if(Ens == "cC211a.06.80") {  Ra0=1e6;  }
	  else if(Ens == "cC211a.06.112") { Ra0= 1e6;   }
	  else if(Ens == "cD211a.054.96") { Ra0=5e7;   }
	  else if(Ens == "cE211a.044.112") {  Ra0=3.5e8;  }
	  else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
	}
	else if(channel=="V0") {
	  if(Ens == "cB211b.072.64") { Ra0=1e6;   }
	  else if(Ens == "cB211b.072.96") { Ra0=1e6;   }
	  else if(Ens == "cC211a.06.80") {  Ra0=3e6;  }
	  else if(Ens == "cC211a.06.112") { Ra0= 4e7;   }
	  else if(Ens == "cD211a.054.96") { Ra0=2.5e7;   }
	  else if(Ens == "cE211a.044.112") {  Ra0=1e7;  }
	  else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
	}
	else if(channel=="T") {
	  if(Ens == "cB211b.072.64") { Ra0=2e4; if( s > 0.11) Ra0=9e3;   }
	  else if(Ens == "cB211b.072.96") { Ra0=2e5; if( s >0.11) Ra0=1.6e5;   }
	  else if(Ens == "cC211a.06.80") { Ra0=6e5; if (s > 0.11) Ra0 = 3e4;   }
	  else if(Ens == "cC211a.06.112") { Ra0= 4e5;  if(s > 0.11) Ra0=6e4;  }
	  else if(Ens == "cD211a.054.96") { Ra0=7.5e4; if(s > 0.11) Ra0=3e4;  }
	  else if(Ens == "cE211a.044.112") {  Ra0=9e5; if( s > 0.11) Ra0=1e5;  }
	  else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
	}
	else if(channel=="L") {
	  if(Ens == "cB211b.072.64") { Ra0=2e6;   }
	  else if(Ens == "cB211b.072.96") { Ra0=2e6;   }
	  else if(Ens == "cC211a.06.80") { Ra0=5e5;   }
	  else if(Ens == "cC211a.06.112") { Ra0= 5e5;   }
	  else if(Ens == "cD211a.054.96") { Ra0=5e6;  }
	  else if(Ens == "cE211a.044.112") {  Ra0=6e5;  }
	  else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
	}
	else crash("In Customized_plateaux_tau_spectre, channel: "+channel+" not yet implemented");
	
      }
      else if(reg=="OS") {
	
	if(channel=="Aii") {
	  if(Ens == "cB211b.072.64") { Ra0=1.7e6;   }
	  else if(Ens == "cB211b.072.96") { Ra0=1.5e7;   }
	  else if(Ens == "cC211a.06.80") { Ra0=1.2e7;   } 
	  else if(Ens == "cC211a.06.112") { Ra0= 6e6;   }
	  else if(Ens == "cD211a.054.96") {  Ra0=5.5e7;  }
	  else if(Ens == "cE211a.044.112") {  Ra0=8e6;  }
	  else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
	}
	else if(channel=="Vii") {
	  if(Ens == "cB211b.072.64") { Ra0=2.5e6; } // Ra0= 6e5;   }
	  else if(Ens == "cB211b.072.96") {  Ra0= 9e7;    }
	  else if(Ens == "cC211a.06.80") {   Ra0= 2e7;   }
	  else if(Ens == "cC211a.06.112") { Ra0= 2e7;   }
	  else if(Ens == "cD211a.054.96") {  Ra0=2e7;  }
	  else if(Ens == "cE211a.044.112") {  Ra0=2e7;  }
	  else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
	}
	else if(channel=="A0") {
	  if(Ens == "cB211b.072.64") { Ra0=2e5;   }
	  else if(Ens == "cB211b.072.96") { Ra0=1.6e6;   }
	  else if(Ens == "cC211a.06.80") { Ra0=1e6;   }
	  else if(Ens == "cC211a.06.112") { Ra0= 1e7;   }
	  else if(Ens == "cD211a.054.96") { Ra0=3e6;   }
	  else if(Ens == "cE211a.044.112") {  Ra0=1e8;  }
	  else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
	}
	else if(channel=="V0") {
	  if(Ens == "cB211b.072.64") { Ra0=5e7;   }
	  else if(Ens == "cB211b.072.96") { Ra0=2e8;   }
	  else if(Ens == "cC211a.06.80") { Ra0=1e9;   }
	  else if(Ens == "cC211a.06.112") { Ra0= 7e8;   }
	  else if(Ens == "cD211a.054.96") { Ra0=5e9;  }
	  else if(Ens == "cE211a.044.112") {  Ra0=5e3;  }
	  else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
	}
	else if(channel=="T") {
	  if(Ens == "cB211b.072.64") { Ra0=1.3e4;  if(s > 0.11) Ra0= 6e3; }
	  else if(Ens == "cB211b.072.96") { Ra0=4e5;   }
	  else if(Ens == "cC211a.06.80") { Ra0=2e5;  if(s > 0.11) Ra0=2e4; }
	  else if(Ens == "cC211a.06.112") { Ra0= 2e5; if(s > 0.11) Ra0=3e4;  }
	  else if(Ens == "cD211a.054.96") { Ra0=3.5e4;  }
	  else if(Ens == "cE211a.044.112") {  Ra0=2e5; if( s > 0.11) Ra0=5e4;  }
	  else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
	}
	else if(channel=="L") {
	  if(Ens == "cB211b.072.64") { Ra0=1e5;   }
	  else if(Ens == "cB211b.072.96") { Ra0=3e5;   }
	  else if(Ens == "cC211a.06.80") { Ra0=1e6;   }
	  else if(Ens == "cC211a.06.112") { Ra0= 5e6;   }
	  else if(Ens == "cD211a.054.96") { Ra0=3e6;  }
	  else if(Ens == "cE211a.044.112") {  Ra0=2e7;  }
	  else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
	}
	else crash("In Customized_plateaux_tau_spectre, channel: "+channel+" not yet implemented");
	
      }
      
      else crash("In Customized_plateaux_tau_spectre, reg: "+reg+" not yet implemented");
    }
    else { // s > 0.12


      if( reg=="tm") {
	if(channel=="Aii") {
	  if(Ens == "cB211b.072.64") {  Ra0 =3e6;  }
	  else if(Ens == "cB211b.072.96") { Ra0= 4e8;   }
	  else if(Ens == "cC211a.06.80") { Ra0= 6e6;   } 
	  else if(Ens == "cC211a.06.112") { Ra0= 6e6;   }
	  else if(Ens == "cD211a.054.96") {  Ra0=8e6;  }
	  else if(Ens == "cE211a.044.112") {  Ra0=4e6;  }
	  else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
	}
	else if(channel=="Vii") {
	  if(Ens == "cB211b.072.64") { Ra0=8e5;   }
	  else if(Ens == "cB211b.072.96") { Ra0= 2.5e7;    } 
	  else if(Ens == "cC211a.06.80") {  Ra0=3e7;  } 
	  else if(Ens == "cC211a.06.112") { Ra0= 2e6;   }
	  else if(Ens == "cD211a.054.96") {  Ra0= 2e7;   }
	  else if(Ens == "cE211a.044.112") {  Ra0=3e7;  }
	  else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
	}
	else if(channel=="A0") {
	  if(Ens == "cB211b.072.64") { Ra0=1e7;   }
	  else if(Ens == "cB211b.072.96") { Ra0=1e7;   }
	  else if(Ens == "cC211a.06.80") {  Ra0=1e6;  }
	  else if(Ens == "cC211a.06.112") { Ra0= 1e6;   }
	  else if(Ens == "cD211a.054.96") { Ra0=5e7;   }
	  else if(Ens == "cE211a.044.112") {  Ra0=3.5e8;  }
	  else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
	}
	else if(channel=="V0") {
	  if(Ens == "cB211b.072.64") { Ra0=1e6;   }
	  else if(Ens == "cB211b.072.96") { Ra0=1e6;   }
	  else if(Ens == "cC211a.06.80") {  Ra0=3e6;  }
	  else if(Ens == "cC211a.06.112") { Ra0= 4e7;   }
	  else if(Ens == "cD211a.054.96") { Ra0=2.5e7;   }
	  else if(Ens == "cE211a.044.112") {  Ra0=1e7;  }
	  else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
	}
	else if(channel=="T") {
	  if(Ens == "cB211b.072.64") { Ra0=2e4;   }
	  else if(Ens == "cB211b.072.96") { Ra0=2e5;   }
	  else if(Ens == "cC211a.06.80") { Ra0=3e4;   }
	  else if(Ens == "cC211a.06.112") { Ra0= 2e4;   }
	  else if(Ens == "cD211a.054.96") { Ra0=3e4;  }
	  else if(Ens == "cE211a.044.112") {  Ra0=2e5;  }
	  else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
	}
	else if(channel=="L") {
	  if(Ens == "cB211b.072.64") { Ra0=2e6;   }
	  else if(Ens == "cB211b.072.96") { Ra0=2e6;   }
	  else if(Ens == "cC211a.06.80") { Ra0=5e5;   }
	  else if(Ens == "cC211a.06.112") { Ra0= 5e5;   }
	  else if(Ens == "cD211a.054.96") { Ra0=5e6;  }
	  else if(Ens == "cE211a.044.112") {  Ra0=6e5;  }
	  else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
	}
	else crash("In Customized_plateaux_tau_spectre, channel: "+channel+" not yet implemented");
	
      }
      else if(reg=="OS") {
	
	if(channel=="Aii") {
	  if(Ens == "cB211b.072.64") { Ra0=1.7e6;   }
	  else if(Ens == "cB211b.072.96") { Ra0=1.5e7;   }
	  else if(Ens == "cC211a.06.80") { Ra0=1.2e7;   } 
	  else if(Ens == "cC211a.06.112") { Ra0= 6e6;   }
	  else if(Ens == "cD211a.054.96") {  Ra0=5.5e7;  }
	  else if(Ens == "cE211a.044.112") {  Ra0=8e6;  }
	  else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
	}
	else if(channel=="Vii") {
	  if(Ens == "cB211b.072.64") { Ra0=2.5e6; } // Ra0= 6e5;   }
	  else if(Ens == "cB211b.072.96") {  Ra0= 9e7;    }
	  else if(Ens == "cC211a.06.80") {   Ra0= 2e7;   }
	  else if(Ens == "cC211a.06.112") { Ra0= 2e7;   }
	  else if(Ens == "cD211a.054.96") {  Ra0=2e7;  }
	  else if(Ens == "cE211a.044.112") {  Ra0=2e7;  }
	  else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
	}
	else if(channel=="A0") {
	  if(Ens == "cB211b.072.64") { Ra0=2e5;   }
	  else if(Ens == "cB211b.072.96") { Ra0=1.6e6;   }
	  else if(Ens == "cC211a.06.80") { Ra0=1e6;   }
	  else if(Ens == "cC211a.06.112") { Ra0= 1e7;   }
	  else if(Ens == "cD211a.054.96") { Ra0=3e6;   }
	  else if(Ens == "cE211a.044.112") {  Ra0=1e8;  }
	  else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
	}
	else if(channel=="V0") {
	  if(Ens == "cB211b.072.64") { Ra0=5e7;   }
	  else if(Ens == "cB211b.072.96") { Ra0=2e8;   }
	  else if(Ens == "cC211a.06.80") { Ra0=1e9;   }
	  else if(Ens == "cC211a.06.112") { Ra0= 7e8;   }
	  else if(Ens == "cD211a.054.96") { Ra0=5e9;  }
	  else if(Ens == "cE211a.044.112") {  Ra0=5e3;  }
	  else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
	}
	else if(channel=="T") {
	  if(Ens == "cB211b.072.64") { Ra0=2e3;   }
	  else if(Ens == "cB211b.072.96") { Ra0=2e5;   }
	  else if(Ens == "cC211a.06.80") { Ra0=2e4;   }
	  else if(Ens == "cC211a.06.112") { Ra0= 2e4;   }
	  else if(Ens == "cD211a.054.96") { Ra0=3.5e4;  }
	  else if(Ens == "cE211a.044.112") {  Ra0=2e5;  }
	  else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
	}
	else if(channel=="L") {
	  if(Ens == "cB211b.072.64") { Ra0=1e5;   }
	  else if(Ens == "cB211b.072.96") { Ra0=3e5;   }
	  else if(Ens == "cC211a.06.80") { Ra0=1e6;   }
	  else if(Ens == "cC211a.06.112") { Ra0= 5e6;   }
	  else if(Ens == "cD211a.054.96") { Ra0=3e6;  }
	  else if(Ens == "cE211a.044.112") {  Ra0=2e7;  }
	  else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
	}
	else crash("In Customized_plateaux_tau_spectre, channel: "+channel+" not yet implemented");
	
      }
      
      else crash("In Customized_plateaux_tau_spectre, reg: "+reg+" not yet implemented");
      
      
    }
    
  }
  else {  // alpha=3

    if(s < 0.125) {

     if( reg=="tm") {
       if(channel=="Aii") {
	 if(Ens == "cB211b.072.64") {  Ra0 =3e6;  }
	 else if(Ens == "cB211b.072.96") { Ra0= 4e8;   }
	 else if(Ens == "cC211a.06.80") { Ra0= 6e6;   } 
	 else if(Ens == "cC211a.06.112") { Ra0= 6e6;   }
	 else if(Ens == "cD211a.054.96") {  Ra0=8e6;  }
	 else if(Ens == "cE211a.044.112") {  Ra0=4e6;  }
	 else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
       }
       else if(channel=="Vii") {
	 if(Ens == "cB211b.072.64") { Ra0=8e5;   }
	 else if(Ens == "cB211b.072.96") { Ra0= 2.5e7;    } 
	 else if(Ens == "cC211a.06.80") {  Ra0=3e7;  } 
	 else if(Ens == "cC211a.06.112") { Ra0= 2e6;   }
	 else if(Ens == "cD211a.054.96") {  Ra0= 2e7;   }
	 else if(Ens == "cE211a.044.112") {  Ra0=3e7;  }
	 else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
       }
       else if(channel=="A0") {
	 if(Ens == "cB211b.072.64") { Ra0=1e7;   }
	 else if(Ens == "cB211b.072.96") { Ra0=1e7;   }
	 else if(Ens == "cC211a.06.80") {  Ra0=1e6;  }
	 else if(Ens == "cC211a.06.112") { Ra0= 1e6;   }
	 else if(Ens == "cD211a.054.96") { Ra0=5e7;   }
	 else if(Ens == "cE211a.044.112") {  Ra0=3.5e8;  }
	 else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
       }
       else if(channel=="V0") {
	 if(Ens == "cB211b.072.64") { Ra0=1e6;   }
	 else if(Ens == "cB211b.072.96") { Ra0=1e6;   }
	 else if(Ens == "cC211a.06.80") {  Ra0=3e6;  }
	 else if(Ens == "cC211a.06.112") { Ra0= 4e7;   }
	 else if(Ens == "cD211a.054.96") { Ra0=2.5e7;   }
	 else if(Ens == "cE211a.044.112") {  Ra0=1e7;  }
	 else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
       }
       else if(channel=="T") {
	 if(Ens == "cB211b.072.64") { Ra0=8e3; if(s>0.11) Ra0=4.8e3;   }
	 else if(Ens == "cB211b.072.96") { Ra0=2e5;   }
	 else if(Ens == "cC211a.06.80") { Ra0=2e4; if(s>0.11) Ra0=1e4;  }
	 else if(Ens == "cC211a.06.112") { Ra0= 2e4; if(s>0.11) Ra0=1e4;  }
	 else if(Ens == "cD211a.054.96") { Ra0=4e4;  }
	 else if(Ens == "cE211a.044.112") {  Ra0=1.5e5; if(s > 0.11) Ra0=1e5;  }
	 else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
       }
       else if(channel=="L") {
	 if(Ens == "cB211b.072.64") { Ra0=8e5;   }
	 else if(Ens == "cB211b.072.96") { Ra0=8e5;   }
	 else if(Ens == "cC211a.06.80") { Ra0=1.1e5;   }
	 else if(Ens == "cC211a.06.112") { Ra0= 2.5e5;   }
	 else if(Ens == "cD211a.054.96") { Ra0=2e6;  }
	 else if(Ens == "cE211a.044.112") {  Ra0=6e5;  }
	 else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
       }
       else crash("In Customized_plateaux_tau_spectre, channel: "+channel+" not yet implemented");
       
     }
     else if(reg=="OS") {
       
       if(channel=="Aii") {
	 if(Ens == "cB211b.072.64") { Ra0=1.7e6;   }
	 else if(Ens == "cB211b.072.96") { Ra0=1.5e7;   }
	 else if(Ens == "cC211a.06.80") { Ra0=1.2e7;   }  //increase E_th
	 else if(Ens == "cC211a.06.112") { Ra0= 6e6;   }
	 else if(Ens == "cD211a.054.96") {  Ra0=5.5e7;  }
	 else if(Ens == "cE211a.044.112") {  Ra0=8e6;  }
	 else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
       }
       else if(channel=="Vii") {
	 if(Ens == "cB211b.072.64") { Ra0=2.5e6; } // Ra0= 6e5;   }
	 else if(Ens == "cB211b.072.96") {  Ra0= 9e7;    }
	 else if(Ens == "cC211a.06.80") {   Ra0= 2e7;   }
	 else if(Ens == "cC211a.06.112") { Ra0= 2e7;   }
	 else if(Ens == "cD211a.054.96") {  Ra0=2e7;  }
	 else if(Ens == "cE211a.044.112") {  Ra0=2e7;  }
	 else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
       }
       else if(channel=="A0") {
	 if(Ens == "cB211b.072.64") { Ra0=2e5;   }
	 else if(Ens == "cB211b.072.96") { Ra0=1.6e6;   }
	 else if(Ens == "cC211a.06.80") { Ra0=1e6;   }
	 else if(Ens == "cC211a.06.112") { Ra0= 1e7;   }
	 else if(Ens == "cD211a.054.96") { Ra0=3e6;   }
	 else if(Ens == "cE211a.044.112") {  Ra0=1e8;  }
	 else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
    }
       else if(channel=="V0") {
	 if(Ens == "cB211b.072.64") { Ra0=5e7;   }
	 else if(Ens == "cB211b.072.96") { Ra0=2e8;   }
	 else if(Ens == "cC211a.06.80") { Ra0=1e9;   }
	 else if(Ens == "cC211a.06.112") { Ra0= 7e8;   }
	 else if(Ens == "cD211a.054.96") { Ra0=5e9;  }
	 else if(Ens == "cE211a.044.112") {  Ra0=5e3;  }
	 else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
    }
       else if(channel=="T") {
	 if(Ens == "cB211b.072.64") { Ra0=4.5e3;   }
	 else if(Ens == "cB211b.072.96") { Ra0=2e5;   }
	 else if(Ens == "cC211a.06.80") { Ra0=2e4; if( s > 0.11) Ra0= 7e3;  }
	 else if(Ens == "cC211a.06.112") { Ra0= 2e4; if( s > 0.11) Ra0=7e3;   }
	 else if(Ens == "cD211a.054.96") { Ra0=3.5e4;  }
	 else if(Ens == "cE211a.044.112") {  Ra0=4e4;  }
	 else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
       }
       else if(channel=="L") {
	 if(Ens == "cB211b.072.64") { Ra0=1e5;   }
	 else if(Ens == "cB211b.072.96") { Ra0=3e5;   }
	 else if(Ens == "cC211a.06.80") { Ra0=1e6;   }
	 else if(Ens == "cC211a.06.112") { Ra0= 1.5e6;   }
	 else if(Ens == "cD211a.054.96") { Ra0=3e6;  }
	 else if(Ens == "cE211a.044.112") {  Ra0=2e7;  }
	 else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
       }
       else crash("In Customized_plateaux_tau_spectre, channel: "+channel+" not yet implemented");
       
     }
     
     else crash("In Customized_plateaux_tau_spectre, reg: "+reg+" not yet implemented");

    }

    else { // sigma > 0.12

      
      if( reg=="tm") {
	if(channel=="Aii") {
	  if(Ens == "cB211b.072.64") {  Ra0 =3e6;  }
	  else if(Ens == "cB211b.072.96") { Ra0= 4e8;   }
	  else if(Ens == "cC211a.06.80") { Ra0= 6e6;   } 
	  else if(Ens == "cC211a.06.112") { Ra0= 6e6;   }
	  else if(Ens == "cD211a.054.96") {  Ra0=8e6;  }
	  else if(Ens == "cE211a.044.112") {  Ra0=4e6;  }
	  else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
	}
	else if(channel=="Vii") {
	 if(Ens == "cB211b.072.64") { Ra0=8e5;   }
	 else if(Ens == "cB211b.072.96") { Ra0= 2.5e7;    } 
	 else if(Ens == "cC211a.06.80") {  Ra0=3e7;  } 
	 else if(Ens == "cC211a.06.112") { Ra0= 2e6;   }
	 else if(Ens == "cD211a.054.96") {  Ra0= 2e7;   }
	 else if(Ens == "cE211a.044.112") {  Ra0=3e7;  }
	 else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
	}
	else if(channel=="A0") {
	  if(Ens == "cB211b.072.64") { Ra0=1e7;   }
	  else if(Ens == "cB211b.072.96") { Ra0=1e7;   }
	  else if(Ens == "cC211a.06.80") {  Ra0=1e6;  }
	  else if(Ens == "cC211a.06.112") { Ra0= 1e6;   }
	  else if(Ens == "cD211a.054.96") { Ra0=5e7;   }
	  else if(Ens == "cE211a.044.112") {  Ra0=3.5e8;  }
	  else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
	}
	else if(channel=="V0") {
	 if(Ens == "cB211b.072.64") { Ra0=1e6;   }
	 else if(Ens == "cB211b.072.96") { Ra0=1e6;   }
	 else if(Ens == "cC211a.06.80") {  Ra0=3e6;  }
	 else if(Ens == "cC211a.06.112") { Ra0= 4e7;   }
	 else if(Ens == "cD211a.054.96") { Ra0=2.5e7;   }
	 else if(Ens == "cE211a.044.112") {  Ra0=1e7;  }
	 else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
	}
       else if(channel=="T") {
	 if(Ens == "cB211b.072.64") { Ra0=1e3;   }
	 else if(Ens == "cB211b.072.96") { Ra0=4e4;   }
	 else if(Ens == "cC211a.06.80") { Ra0=2e3;   }
	 else if(Ens == "cC211a.06.112") { Ra0= 2e3;   }
	 else if(Ens == "cD211a.054.96") { Ra0=9e3;  }
	 else if(Ens == "cE211a.044.112") {  Ra0=1.5e4;  }
	 else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
       }
       else if(channel=="L") {
	 if(Ens == "cB211b.072.64") { Ra0=8e5;   }
	 else if(Ens == "cB211b.072.96") { Ra0=8e5;   }
	 else if(Ens == "cC211a.06.80") { Ra0=1e5;   }
	 else if(Ens == "cC211a.06.112") { Ra0= 1e5;   }
	 else if(Ens == "cD211a.054.96") { Ra0=2e6;  }
	 else if(Ens == "cE211a.044.112") {  Ra0=1.5e5;  }
	 else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
       }
       else crash("In Customized_plateaux_tau_spectre, channel: "+channel+" not yet implemented");
	
      }
      else if(reg=="OS") {
	
	if(channel=="Aii") {
	  if(Ens == "cB211b.072.64") { Ra0=1.7e6;   }
	  else if(Ens == "cB211b.072.96") { Ra0=1.5e7;   }
	  else if(Ens == "cC211a.06.80") { Ra0=1.2e7;   }  //increase E_th
	  else if(Ens == "cC211a.06.112") { Ra0= 6e6;   }
	  else if(Ens == "cD211a.054.96") {  Ra0=5.5e7;  }
	  else if(Ens == "cE211a.044.112") {  Ra0=8e6;  }
	  else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
	}
	else if(channel=="Vii") {
	  if(Ens == "cB211b.072.64") { Ra0=2.5e6; } // Ra0= 6e5;   }
	  else if(Ens == "cB211b.072.96") {  Ra0= 9e7;    }
	  else if(Ens == "cC211a.06.80") {   Ra0= 2e7;   }
	  else if(Ens == "cC211a.06.112") { Ra0= 2e7;   }
	  else if(Ens == "cD211a.054.96") {  Ra0=2e7;  }
	  else if(Ens == "cE211a.044.112") {  Ra0=2e7;  }
	  else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
	}
	else if(channel=="A0") {
	  if(Ens == "cB211b.072.64") { Ra0=2e5;   }
	  else if(Ens == "cB211b.072.96") { Ra0=1.6e6;   }
	  else if(Ens == "cC211a.06.80") { Ra0=1e6;   }
	  else if(Ens == "cC211a.06.112") { Ra0= 1e7;   }
	  else if(Ens == "cD211a.054.96") { Ra0=3e6;   }
	  else if(Ens == "cE211a.044.112") {  Ra0=1e8;  }
	  else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
	}
	else if(channel=="V0") {
	  if(Ens == "cB211b.072.64") { Ra0=5e7;   }
	  else if(Ens == "cB211b.072.96") { Ra0=2e8;   }
	  else if(Ens == "cC211a.06.80") { Ra0=1e9;   }
	  else if(Ens == "cC211a.06.112") { Ra0= 7e8;   }
	  else if(Ens == "cD211a.054.96") { Ra0=5e9;  }
	  else if(Ens == "cE211a.044.112") {  Ra0=5e3;  }
	  else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
	}
	else if(channel=="T") {
	  if(Ens == "cB211b.072.64") { Ra0=1e3;   }
	  else if(Ens == "cB211b.072.96") { Ra0=4e4;   }
	  else if(Ens == "cC211a.06.80") { Ra0=3e3;   }
	  else if(Ens == "cC211a.06.112") { Ra0= 3e3;   }
	  else if(Ens == "cD211a.054.96") { Ra0=5e3;  }
	  else if(Ens == "cE211a.044.112") {  Ra0=1e4;  }
	  else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
	}
	else if(channel=="L") {
	  if(Ens == "cB211b.072.64") { Ra0=1e5;   }
	  else if(Ens == "cB211b.072.96") { Ra0=3e5;   }
	  else if(Ens == "cC211a.06.80") { Ra0=5e5;   }
	  else if(Ens == "cC211a.06.112") { Ra0= 5e5;   }
	  else if(Ens == "cD211a.054.96") { Ra0=3e6;  }
	  else if(Ens == "cE211a.044.112") {  Ra0=2e7;  }
	  else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
	}
	else crash("In Customized_plateaux_tau_spectre, channel: "+channel+" not yet implemented");
	
     }
      
      else crash("In Customized_plateaux_tau_spectre, reg: "+reg+" not yet implemented");
      
      
      

    }
     
     
     
     
  }
  



  assertm(Ra0 > 0, "Assertion Ra0 > 0 failed");
  return Ra0;


}



void get_sigma_list_strange() {

  bool test=false;
  if(test) { sigma_list_strange.push_back(0.05); return;}


  double s_max= 0.2;
  double s_min= 0.004420;

  //sigma_list_strange= { 0.004, 0.02, 0.04, 0.08, 0.12, 0.16};
  sigma_list_strange= { 0.01, 0.02, 0.04, 0.08, 0.12, 0.14, 0.16};
  return;

 
  return;

}



void tau_decay_analysis_strange() {



   //Analyze new strange run
  //###################################

  bool Get_ASCII= false;

    if(Get_ASCII) {
    //read binary files
    boost::filesystem::create_directory("../tau_decay_strange");
    

    vector<string> Ens_T1({"C.06.80", "C.06.112", "B.72.64", "B.72.96" , "D.54.96", "E.44.112"});
    vector<string> Ens_TT1({"cC211a.06.80", "cC211a.06.112", "cB211b.072.64", "cB211b.072.96" , "cD211a.054.96", "cE211a.044.112"});

    for( int it=0; it<(signed)Ens_T1.size(); it++) {

      vector<string> channels({"mix_l_l",  "mix_l_s1", "mix_l_s2", "mix_s1_s1", "mix_s2_s2"});

      for(auto &channel : channels) {
	boost::filesystem::create_directory("../tau_decay_strange/"+channel);
	boost::filesystem::create_directory("../tau_decay_strange/"+channel+"/"+Ens_TT1[it]);
      }
      //read binary
      vector<string> Corr_tags({"TM_VKVK", "TM_AKAK", "TM_A0A0", "TM_V0V0", "TM_P5P5", "TM_S0S0", "OS_VKVK", "OS_AKAK", "OS_A0A0", "OS_V0V0", "OS_P5P5", "OS_S0S0"});

          
      for(int id=0; id<(signed)Corr_tags.size(); id++) {



	for( auto &channel: channels) {

	FILE *stream = fopen( ("../gm2_tau_rep_bin/"+Ens_T1[it]+"/"+channel+"_"+Corr_tags[id]).c_str(), "rb");
        size_t Nconfs, T, Nhits;
	bin_read(Nconfs, stream);
	bin_read(Nhits, stream);
	bin_read(T, stream);

	cout<<"Nconfs: "<<Nconfs<<endl;
	cout<<"T: "<<T<<" "<<T/2+1<<endl;
	cout<<"Nhits: "<<Nhits<<endl;
	for(size_t iconf=0;iconf<Nconfs;iconf++) {
	  vector<double> C(T/2+1);
	  for(size_t t=0;t<T/2+1;t++) bin_read(C[t], stream);
	  boost::filesystem::create_directory("../tau_decay_strange/"+channel+"/"+Ens_TT1[it]+"/"+to_string(iconf));
	  ofstream PrintCorr("../tau_decay_strange/"+channel+"/"+Ens_TT1[it]+"/"+to_string(iconf)+"/mes_contr_"+channel+"_"+Corr_tags[id]);
	  PrintCorr.precision(16);
	  PrintCorr<<"# "<<Corr_tags[id].substr(3,4)<<endl;
	  for(size_t t=0;t<(T/2+1);t++) PrintCorr<<C[t]<<endl;
	  if(Corr_tags[id].substr(3,4) == "VKTK" || Corr_tags[id].substr(3,4) == "TKVK") { for(size_t t=T/2+1; t<T;t++) PrintCorr<<-1*C[T-t]<<endl;   }
	  else  {for(size_t t=T/2+1; t<T;t++) PrintCorr<<C[T-t]<<endl;}
	  PrintCorr.close();

	}

	fclose(stream);

	}
	
      }
    }
    }

    bool Get_ASCII_dm=false;

    if(Get_ASCII_dm) {

      cout<<"Analyzing light-quark mass corrections"<<endl;

    //mass correction
     boost::filesystem::create_directory("../tau_decay_strange/light_mass_correction");
     boost::filesystem::create_directory("../tau_decay_strange/light_mass_correction/mix_l1_s");
     boost::filesystem::create_directory("../tau_decay_strange/light_mass_correction/mix_l2_s");


     vector<string> Ens_T1 = {"C.06.80", "C.06.112", "B.72.64", "B.72.96" , "D.54.96", "E.44.112"};
     vector<string> Ens_TT1 = {"cC211a.06.80", "cC211a.06.112", "cB211b.072.64", "cB211b.072.96" , "cD211a.054.96", "cE211a.044.112"};


     for( int it=0; it<(signed)Ens_T1.size(); it++) {

       vector<string> channels({"mix_l1_s", "mix_l2_s"});

      for(auto &channel : channels) {
	boost::filesystem::create_directory("../tau_decay_strange/light_mass_correction/"+channel);
	boost::filesystem::create_directory("../tau_decay_strange/light_mass_correction/"+channel+"/"+Ens_TT1[it]);
      }
      //read binary
      vector<string> Corr_tags({"TM_VKVK", "TM_AKAK", "TM_A0A0", "TM_V0V0", "TM_P5P5", "TM_S0S0", "OS_VKVK", "OS_AKAK", "OS_A0A0", "OS_V0V0", "OS_P5P5", "OS_S0S0"});

          
      for(int id=0; id<(signed)Corr_tags.size(); id++) {



	for( auto &channel: channels) {

	FILE *stream = fopen( ("../tau_decay_strange_bin_mu_corr/"+Ens_T1[it]+"/"+channel+"_"+Corr_tags[id]).c_str(), "rb");
        size_t Nconfs, T, Nhits;
	bin_read(Nconfs, stream);
	bin_read(Nhits, stream);
	bin_read(T, stream);

	cout<<"Nconfs: "<<Nconfs<<endl;
	cout<<"T: "<<T<<" "<<T/2+1<<endl;
	cout<<"Nhits: "<<Nhits<<endl;
	for(size_t iconf=0;iconf<Nconfs;iconf++) {
	  vector<double> C(T/2+1);
	  for(size_t t=0;t<T/2+1;t++) bin_read(C[t], stream);
	  boost::filesystem::create_directory("../tau_decay_strange/light_mass_correction/"+channel+"/"+Ens_TT1[it]+"/"+to_string(iconf));
	  ofstream PrintCorr("../tau_decay_strange/light_mass_correction/"+channel+"/"+Ens_TT1[it]+"/"+to_string(iconf)+"/mes_contr_"+channel+"_"+Corr_tags[id]);
	  PrintCorr.precision(16);
	  PrintCorr<<"# "<<Corr_tags[id].substr(3,4)<<endl;
	  for(size_t t=0;t<(T/2+1);t++) PrintCorr<<C[t]<<endl;
	  if(Corr_tags[id].substr(3,4) == "VKTK" || Corr_tags[id].substr(3,4) == "TKVK") { for(size_t t=T/2+1; t<T;t++) PrintCorr<<-1*C[T-t]<<endl;   }
	  else  {for(size_t t=T/2+1; t<T;t++) PrintCorr<<C[T-t]<<endl;}
	  PrintCorr.close();

	}

	fclose(stream);

	}
	
      }
     }

   

    
    }



  get_sigma_list_strange();


  //Vfloat betas({ 3.99, 2.99, 4.99, 5.99, 2.99, 1.99 });
  //Vfloat Emax_list({4.0, 4.0 , 4.0, 4.0, 5.0, 4.0});
  //vector<bool> Is_Emax_Finite({1,1,1,1,1,0});

  //Vfloat betas({ 3.99, 2.99, 4.99, 5.99, 2.99});
  //Vfloat Emax_list({4.0, 4.0 , 4.0, 4.0, 5.0});
  //vector<bool> Is_Emax_Finite({1,1,1,1,1});
  
  Vfloat betas({ 3.99 , 4.99 , 3.99, 2.99 });
  Vfloat Emax_list({4.0, 4.0 , 5.0, 4.0});
  vector<bool> Is_Emax_Finite({1,1,1, 1});


 

  

  int rank, size;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  
  
  
   

  int N= (signed)betas.size();

  if(N%size != 0) crash("MPI called with -np= "+to_string(size)+". np does not divide vector size N="+to_string(N));
  
  cout<<"################# DETERMINATION OF BRANCHING RATIOS FOR SEMI-INCLUSIVE TAU-DECAY#################"<<endl;
  cout<<"Rank: "<<rank<<" pid: "<<getpid()<<" core id: "<<"("<<sched_getcpu()<<")"<<endl;
  cout<<"INVERSE LAPLACE RECONSTRUCTION CALLED FOR:"<<endl;
  for(int i=rank*N/size;i<(rank+1)*N/size;i++) {
    string alpha_Emax_Tag= "{"+to_string_with_precision(betas[i],2)+","+((Is_Emax_Finite[i]==0)?"inf":to_string_with_precision(Emax_list[i],1))+"}";
    cout<<"{alpha,Emax} = "<<alpha_Emax_Tag<<endl;
  }
  cout<<"##########################################"<<endl;

  vector<distr_t> FIN_RES(betas.size());;
  
  for(int i=rank*N/size; i<(rank+1)*N/size;i++) {FIN_RES[i] = Compute_tau_decay_width_strange(Is_Emax_Finite[i], Emax_list[i], betas[i]);}


  double syst=0;
  distr_t FAVE(UseJack);
  for(int i=0;i<(signed)betas.size();i++) {
    if(betas[i] == 3.99 && Emax_list[i] == 4.0) FAVE=FIN_RES[i];
  }
  
  for(int i=0;i<(signed)betas.size();i++) {

    if(betas[i] != 3.99 || Emax_list[i] != 4.0) {

      if(FIN_RES[i].err() > 0) {

	double ss= fabs( FIN_RES[i].ave() - FAVE.ave())*erf(  fabs( FIN_RES[i].ave() - FAVE.ave())/(sqrt(2.0)*FIN_RES[i].err()));
	if ( ss > syst) syst=ss;
      }
    }
  }

  cout.precision(4);
  cout<<"FINAL RESULT: "<<endl;
  cout<<FAVE.ave()<<" +- "<<FAVE.err()<<" +- "<<syst<<endl;


  cout<<"Final determination of Vus"<<endl;

  	
  distr_t R_exp(UseJack);
  double R_exp_ave= 0.1633;
  double R_exp_err=0.0027;
  double sqrt_R_exp_ave= sqrt(R_exp_ave);
  double sqrt_R_exp_err= 0.5*R_exp_err/sqrt_R_exp_ave;


  if(FAVE.err() > 0) {
  

    distr_t FF= FAVE.ave() + (FAVE-FAVE.ave())*sqrt( pow(FAVE.err(),2) + syst*syst)/FAVE.err();
	
	
    distr_t sqrt_Rinv_lat= SQRT_D(1.0/(FF));

    double sqrt_Rinv_lat_bis_err= (0.5/pow(FF.ave(), 3.0/2.0))*0.022;
  
    cout<<"Vus: "<<sqrt_Rinv_lat.ave()*sqrt_R_exp_ave<<"("<<sqrt_Rinv_lat.err()*sqrt_R_exp_ave<<")_lat ("<<sqrt_Rinv_lat.ave()*sqrt_R_exp_err<<")_exp"<<" ("<<sqrt( pow( sqrt_Rinv_lat.err()*sqrt_R_exp_ave,2) + pow(sqrt_Rinv_lat.ave()*sqrt_R_exp_err,2))<<")"<<endl;
    cout<<"Vus_bis"<<sqrt_Rinv_lat.ave()*sqrt_R_exp_ave<<"("<<sqrt_Rinv_lat_bis_err*sqrt_R_exp_ave<<")_lat ("<<sqrt_Rinv_lat.ave()*sqrt_R_exp_err<<")_exp"<<" ("<<sqrt( pow( sqrt_Rinv_lat.err()*sqrt_R_exp_ave,2) + pow(sqrt_Rinv_lat.ave()*sqrt_R_exp_err,2))<<")"<<endl;
  }
}


distr_t Compute_tau_decay_width_strange(bool Is_Emax_Finite, double Emax, double beta) {


  distr_t RET(UseJack);
  RET= Get_id_jack_distr(Njacks);

  PrecFloat::setDefaultPrecision(prec);
  cout<<"max possible exponent: "<<PrecFloat::getEmax_max()<<endl;
  cout<<"current max exponent: "<<PrecFloat::getEmax()<<endl;

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int _hostname_len;
  char _hostname[MPI_MAX_PROCESSOR_NAME];
  MPI_Get_processor_name(_hostname, &_hostname_len);

   
  string Tag_reco_type="Beta_"+to_string_with_precision(beta,2);
  Tag_reco_type+="_Emax_"+((Is_Emax_Finite==0)?"inf":to_string_with_precision(Emax,1));
  string alpha_Emax_Tag= "{"+to_string_with_precision(beta,2)+","+((Is_Emax_Finite==0)?"inf":to_string_with_precision(Emax,1))+"}";


  cout<<"STARTING COMPUTATION OF: {alpha,Emax} : "<<alpha_Emax_Tag<<endl;
  cout<<"Creating output directories...";
  
  cout.precision(5);



 
 
    


  //create directories

  boost::filesystem::create_directory("../data/tau_decay");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type);
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/continuum");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/continuum/fit_func");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/continuum/fit_data");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/continuum/cov");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/continuum/corr");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/continuum/AIC");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/continuum/fit_data/A0A0");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/continuum/fit_data/A0A0/tm");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/continuum/fit_data/A0A0/OS");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/continuum/fit_data/V0V0");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/continuum/fit_data/V0V0/tm");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/continuum/fit_data/V0V0/OS");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/continuum/fit_data/AkAk");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/continuum/fit_data/AkAk/tm");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/continuum/fit_data/AkAk/OS");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/continuum/fit_data/VkVk");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/continuum/fit_data/VkVk/tm");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/continuum/fit_data/VkVk/OS");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/continuum/fit_data/tot");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/continuum/fit_data/tot/tm");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/continuum/fit_data/tot/OS");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/continuum/fit_data/tot_TL");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/continuum/fit_data/tot_TL/tm");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/continuum/fit_data/tot_TL/OS");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/continuum/fit_data/VA");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/continuum/fit_data/VA/tm");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/continuum/fit_data/VA/OS");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/continuum/fit_data/VMA");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/continuum/fit_data/VMA/tm");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/continuum/fit_data/VMA/OS");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/continuum/fit_data/AX");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/continuum/fit_data/AX/tm");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/continuum/fit_data/AX/OS");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/continuum/fit_data/T");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/continuum/fit_data/T/tm");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/continuum/fit_data/T/OS");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/continuum/fit_data/L");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/continuum/fit_data/L/tm");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/continuum/fit_data/L/OS");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/jackknife");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/jackknife/tm");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/jackknife/tm/A0A0");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/jackknife/tm/V0V0");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/jackknife/tm/AkAk");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/jackknife/tm/VkVk");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/jackknife/tm/T");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/jackknife/tm/L");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/jackknife/OS");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/jackknife/OS/A0A0");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/jackknife/OS/V0V0");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/jackknife/OS/AkAk");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/jackknife/OS/VkVk");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/jackknife/OS/T");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/jackknife/OS/L");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/Br");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/corr");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/mass");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/covariance");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/FSE");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/FSE/tm");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/FSE/OS");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/summary");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/pull");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/FK");
  

  cout<<"done!"<<endl;

 

    
 

    



  //Read data


  //Custom sorting for V_light to account for the two replica r0 and r1


   auto Sort_light_confs = [](string A, string B) {


			     //return A<B;
			     
			    int conf_length_A= A.length();
			    int conf_length_B= B.length();

			    int pos_a_slash=-1;
			    int pos_b_slash=-1;
			    for(int i=0;i<conf_length_A;i++) if(A.substr(i,1)=="/") pos_a_slash=i;
			    for(int j=0;j<conf_length_B;j++) if(B.substr(j,1)=="/") pos_b_slash=j;

			    string A_bis= A.substr(pos_a_slash+1);
			    string B_bis= B.substr(pos_b_slash+1);

			    //A_bis=A;
			    //B_bis=B;

			     
			    string conf_num_A = A_bis.substr(0,4);
			    string conf_num_B = B_bis.substr(0,4);
							       
		      
			    string rA = A_bis.substr(A_bis.length()-5);
			    string rB = B_bis.substr(B_bis.length()-5);
			    if(rA.substr(0,1) == "r") { 
			      int n1 = stoi(rA.substr(1,1));
			      int n2 = stoi(rB.substr(1,1));
			      if(rA == rB) {
			      if(rA=="r0.h5" || rA=="r2.h5") return conf_num_A > conf_num_B;
			      else if(rA=="r1.h5" || rA=="r3.h5") return conf_num_A < conf_num_B;
			      else crash("stream not recognized");
			      }
			      else return n1<n2;
			    }
			    return A_bis<B_bis;
			  };

   auto Sort_easy = [](string A, string B) {

      int conf_length_A= A.length();
      int conf_length_B= B.length();
      
      int pos_a_slash=-1;
      int pos_b_slash=-1;
      for(int i=0;i<conf_length_A;i++) if(A.substr(i,1)=="/") pos_a_slash=i;
      for(int j=0;j<conf_length_B;j++) if(B.substr(j,1)=="/") pos_b_slash=j;
      
      string A_bis= A.substr(pos_a_slash+1);
      string B_bis= B.substr(pos_b_slash+1);

      return atoi( A_bis.c_str()) < atoi( B_bis.c_str());
      
  };
   

  
  //light strange mass

  data_t ll_data_tm_P5P5, ss_data_tm_P5P5, ls_data_tm_P5P5;
  data_t ll_data_tm_AKAK, ss_data_tm_AKAK, ls_data_tm_AKAK;
  data_t ll_data_tm_A0A0, ss_data_tm_A0A0, ls_data_tm_A0A0;
  data_t ll_data_tm_V0V0, ss_data_tm_V0V0, ls_data_tm_V0V0;
  data_t ll_data_tm_VKVK, ss_data_tm_VKVK, ls_data_tm_VKVK;
  data_t ll_data_tm_S0S0, ss_data_tm_S0S0, ls_data_tm_S0S0;
  

  data_t ll_data_OS_P5P5, ss_data_OS_P5P5, ls_data_OS_P5P5;
  data_t ll_data_OS_AKAK, ss_data_OS_AKAK, ls_data_OS_AKAK;
  data_t ll_data_OS_A0A0, ss_data_OS_A0A0, ls_data_OS_A0A0;
  data_t ll_data_OS_V0V0, ss_data_OS_V0V0, ls_data_OS_V0V0;
  data_t ll_data_OS_VKVK, ss_data_OS_VKVK, ls_data_OS_VKVK;
  data_t ll_data_OS_S0S0, ss_data_OS_S0S0, ls_data_OS_S0S0;
  
  ll_data_tm_P5P5.Read("../tau_decay_strange/mix_l_l", "mes_contr_mix_l_l_TM_P5P5", "P5P5", Sort_easy);
  ss_data_tm_P5P5.Read("../tau_decay_strange/mix_s1_s1", "mes_contr_mix_s1_s1_TM_P5P5", "P5P5", Sort_easy);
  ls_data_tm_P5P5.Read("../tau_decay_strange/mix_l_s1", "mes_contr_mix_l_s1_TM_P5P5", "P5P5", Sort_easy);
  ll_data_tm_AKAK.Read("../tau_decay_strange/mix_l_l", "mes_contr_mix_l_l_TM_AKAK", "AKAK", Sort_easy);
  ss_data_tm_AKAK.Read("../tau_decay_strange/mix_s1_s1", "mes_contr_mix_s1_s1_TM_AKAK", "AKAK", Sort_easy);
  ls_data_tm_AKAK.Read("../tau_decay_strange/mix_l_s1", "mes_contr_mix_l_s1_TM_AKAK", "AKAK", Sort_easy);

  ll_data_tm_A0A0.Read("../tau_decay_strange/mix_l_l", "mes_contr_mix_l_l_TM_A0A0", "A0A0", Sort_easy);
  ss_data_tm_A0A0.Read("../tau_decay_strange/mix_s1_s1", "mes_contr_mix_s1_s1_TM_A0A0", "A0A0", Sort_easy);
  ls_data_tm_A0A0.Read("../tau_decay_strange/mix_l_s1", "mes_contr_mix_l_s1_TM_A0A0", "A0A0", Sort_easy);

  ll_data_tm_V0V0.Read("../tau_decay_strange/mix_l_l", "mes_contr_mix_l_l_TM_V0V0", "V0V0", Sort_easy);
  ss_data_tm_V0V0.Read("../tau_decay_strange/mix_s1_s1", "mes_contr_mix_s1_s1_TM_V0V0", "V0V0", Sort_easy);
  ls_data_tm_V0V0.Read("../tau_decay_strange/mix_l_s1", "mes_contr_mix_l_s1_TM_V0V0", "V0V0", Sort_easy);
  
  ll_data_tm_VKVK.Read("../tau_decay_strange/mix_l_l", "mes_contr_mix_l_l_TM_VKVK", "VKVK", Sort_easy);
  ss_data_tm_VKVK.Read("../tau_decay_strange/mix_s1_s1", "mes_contr_mix_s1_s1_TM_VKVK", "VKVK", Sort_easy);
  ls_data_tm_VKVK.Read("../tau_decay_strange/mix_l_s1", "mes_contr_mix_l_s1_TM_VKVK", "VKVK", Sort_easy);
  ll_data_tm_S0S0.Read("../tau_decay_strange/mix_l_l", "mes_contr_mix_l_l_TM_S0S0", "S0S0", Sort_easy);
  ss_data_tm_S0S0.Read("../tau_decay_strange/mix_s1_s1", "mes_contr_mix_s1_s1_TM_S0S0", "S0S0", Sort_easy);
  ls_data_tm_S0S0.Read("../tau_decay_strange/mix_l_s1", "mes_contr_mix_l_s1_TM_S0S0", "S0S0", Sort_easy);


  ll_data_OS_P5P5.Read("../tau_decay_strange/mix_l_l", "mes_contr_mix_l_l_OS_P5P5", "P5P5", Sort_easy);
  ss_data_OS_P5P5.Read("../tau_decay_strange/mix_s1_s1", "mes_contr_mix_s1_s1_OS_P5P5", "P5P5", Sort_easy);
  ls_data_OS_P5P5.Read("../tau_decay_strange/mix_l_s1", "mes_contr_mix_l_s1_OS_P5P5", "P5P5", Sort_easy);
  ll_data_OS_AKAK.Read("../tau_decay_strange/mix_l_l", "mes_contr_mix_l_l_OS_AKAK", "AKAK", Sort_easy);
  ss_data_OS_AKAK.Read("../tau_decay_strange/mix_s1_s1", "mes_contr_mix_s1_s1_OS_AKAK", "AKAK", Sort_easy);
  ls_data_OS_AKAK.Read("../tau_decay_strange/mix_l_s1", "mes_contr_mix_l_s1_OS_AKAK", "AKAK", Sort_easy);

  ll_data_OS_A0A0.Read("../tau_decay_strange/mix_l_l", "mes_contr_mix_l_l_OS_A0A0", "A0A0", Sort_easy);
  ss_data_OS_A0A0.Read("../tau_decay_strange/mix_s1_s1", "mes_contr_mix_s1_s1_OS_A0A0", "A0A0", Sort_easy);
  ls_data_OS_A0A0.Read("../tau_decay_strange/mix_l_s1", "mes_contr_mix_l_s1_OS_A0A0", "A0A0", Sort_easy);

  ll_data_OS_V0V0.Read("../tau_decay_strange/mix_l_l", "mes_contr_mix_l_l_OS_V0V0", "V0V0", Sort_easy);
  ss_data_OS_V0V0.Read("../tau_decay_strange/mix_s1_s1", "mes_contr_mix_s1_s1_OS_V0V0", "V0V0", Sort_easy);
  ls_data_OS_V0V0.Read("../tau_decay_strange/mix_l_s1", "mes_contr_mix_l_s1_OS_V0V0", "V0V0", Sort_easy);
  
  ll_data_OS_VKVK.Read("../tau_decay_strange/mix_l_l", "mes_contr_mix_l_l_OS_VKVK", "VKVK", Sort_easy);
  ss_data_OS_VKVK.Read("../tau_decay_strange/mix_s1_s1", "mes_contr_mix_s1_s1_OS_VKVK", "VKVK", Sort_easy);
  ls_data_OS_VKVK.Read("../tau_decay_strange/mix_l_s1", "mes_contr_mix_l_s1_OS_VKVK", "VKVK", Sort_easy);
  ll_data_OS_S0S0.Read("../tau_decay_strange/mix_l_l", "mes_contr_mix_l_l_OS_S0S0", "S0S0", Sort_easy);
  ss_data_OS_S0S0.Read("../tau_decay_strange/mix_s1_s1", "mes_contr_mix_s1_s1_OS_S0S0", "S0S0", Sort_easy);
  ls_data_OS_S0S0.Read("../tau_decay_strange/mix_l_s1", "mes_contr_mix_l_s1_OS_S0S0", "S0S0", Sort_easy);

  

  //heavier strange mass

  data_t  ss_H_data_tm_P5P5, ls_H_data_tm_P5P5;
  data_t  ss_H_data_tm_AKAK, ls_H_data_tm_AKAK;
  data_t  ss_H_data_tm_A0A0, ls_H_data_tm_A0A0;
  data_t  ss_H_data_tm_V0V0, ls_H_data_tm_V0V0;
  data_t  ss_H_data_tm_VKVK, ls_H_data_tm_VKVK;
  data_t  ss_H_data_tm_S0S0, ls_H_data_tm_S0S0;
  

  data_t ss_H_data_OS_P5P5, ls_H_data_OS_P5P5;
  data_t ss_H_data_OS_AKAK, ls_H_data_OS_AKAK;
  data_t ss_H_data_OS_A0A0, ls_H_data_OS_A0A0;
  data_t ss_H_data_OS_V0V0, ls_H_data_OS_V0V0;
  data_t ss_H_data_OS_VKVK, ls_H_data_OS_VKVK;
  data_t ss_H_data_OS_S0S0, ls_H_data_OS_S0S0;


  ss_H_data_tm_P5P5.Read("../tau_decay_strange/mix_s2_s2", "mes_contr_mix_s2_s2_TM_P5P5", "P5P5", Sort_easy);
  ls_H_data_tm_P5P5.Read("../tau_decay_strange/mix_l_s2", "mes_contr_mix_l_s2_TM_P5P5", "P5P5", Sort_easy);
  ss_H_data_tm_AKAK.Read("../tau_decay_strange/mix_s2_s2", "mes_contr_mix_s2_s2_TM_AKAK", "AKAK", Sort_easy);
  ls_H_data_tm_AKAK.Read("../tau_decay_strange/mix_l_s2", "mes_contr_mix_l_s2_TM_AKAK", "AKAK", Sort_easy);

  ss_H_data_tm_A0A0.Read("../tau_decay_strange/mix_s2_s2", "mes_contr_mix_s2_s2_TM_A0A0", "A0A0", Sort_easy);
  ls_H_data_tm_A0A0.Read("../tau_decay_strange/mix_l_s2", "mes_contr_mix_l_s2_TM_A0A0", "A0A0", Sort_easy);
  
  ss_H_data_tm_V0V0.Read("../tau_decay_strange/mix_s2_s2", "mes_contr_mix_s2_s2_TM_V0V0", "V0V0", Sort_easy);
  ls_H_data_tm_V0V0.Read("../tau_decay_strange/mix_l_s2", "mes_contr_mix_l_s2_TM_V0V0", "V0V0", Sort_easy);
  
  ss_H_data_tm_VKVK.Read("../tau_decay_strange/mix_s2_s2", "mes_contr_mix_s2_s2_TM_VKVK", "VKVK", Sort_easy);
  ls_H_data_tm_VKVK.Read("../tau_decay_strange/mix_l_s2", "mes_contr_mix_l_s2_TM_VKVK", "VKVK", Sort_easy);
  ss_H_data_tm_S0S0.Read("../tau_decay_strange/mix_s2_s2", "mes_contr_mix_s2_s2_TM_S0S0", "S0S0", Sort_easy);
  ls_H_data_tm_S0S0.Read("../tau_decay_strange/mix_l_s2", "mes_contr_mix_l_s2_TM_S0S0", "S0S0", Sort_easy);


  ss_H_data_OS_P5P5.Read("../tau_decay_strange/mix_s2_s2", "mes_contr_mix_s2_s2_OS_P5P5", "P5P5", Sort_easy);
  ls_H_data_OS_P5P5.Read("../tau_decay_strange/mix_l_s2", "mes_contr_mix_l_s2_OS_P5P5", "P5P5", Sort_easy);
  ss_H_data_OS_AKAK.Read("../tau_decay_strange/mix_s2_s2", "mes_contr_mix_s2_s2_OS_AKAK", "AKAK", Sort_easy);
  ls_H_data_OS_AKAK.Read("../tau_decay_strange/mix_l_s2", "mes_contr_mix_l_s2_OS_AKAK", "AKAK", Sort_easy);

  ss_H_data_OS_A0A0.Read("../tau_decay_strange/mix_s2_s2", "mes_contr_mix_s2_s2_OS_A0A0", "A0A0", Sort_easy);
  ls_H_data_OS_A0A0.Read("../tau_decay_strange/mix_l_s2", "mes_contr_mix_l_s2_OS_A0A0", "A0A0", Sort_easy);

  
  ss_H_data_OS_V0V0.Read("../tau_decay_strange/mix_s2_s2", "mes_contr_mix_s2_s2_OS_V0V0", "V0V0", Sort_easy);
  ls_H_data_OS_V0V0.Read("../tau_decay_strange/mix_l_s2", "mes_contr_mix_l_s2_OS_V0V0", "V0V0", Sort_easy);
  
  ss_H_data_OS_VKVK.Read("../tau_decay_strange/mix_s2_s2", "mes_contr_mix_s2_s2_OS_VKVK", "VKVK", Sort_easy);
  ls_H_data_OS_VKVK.Read("../tau_decay_strange/mix_l_s2", "mes_contr_mix_l_s2_OS_VKVK", "VKVK", Sort_easy);
  ss_H_data_OS_S0S0.Read("../tau_decay_strange/mix_s2_s2", "mes_contr_mix_s2_s2_OS_S0S0", "S0S0", Sort_easy);
  ls_H_data_OS_S0S0.Read("../tau_decay_strange/mix_l_s2", "mes_contr_mix_l_s2_OS_S0S0", "S0S0", Sort_easy);




  //light-quark mass correction

  data_t ls_ph_data_tm_VKVK, ls_ph_data_tm_V0V0, ls_ph_data_tm_AKAK, ls_ph_data_tm_A0A0;
  data_t ls_ph_data_OS_VKVK, ls_ph_data_OS_V0V0, ls_ph_data_OS_AKAK, ls_ph_data_OS_A0A0;
  data_t ls_ph_data_tm_P5P5;


  //light 
  //tm
  ls_ph_data_tm_VKVK.Read("../tau_decay_strange/light_mass_correction/mix_l2_s", "mes_contr_mix_l2_s_TM_VKVK", "VKVK", Sort_easy);
  ls_ph_data_tm_V0V0.Read("../tau_decay_strange/light_mass_correction/mix_l2_s", "mes_contr_mix_l2_s_TM_V0V0", "V0V0", Sort_easy);
  ls_ph_data_tm_AKAK.Read("../tau_decay_strange/light_mass_correction/mix_l2_s", "mes_contr_mix_l2_s_TM_AKAK", "AKAK", Sort_easy);
  ls_ph_data_tm_A0A0.Read("../tau_decay_strange/light_mass_correction/mix_l2_s", "mes_contr_mix_l2_s_TM_A0A0", "A0A0", Sort_easy);
  ls_ph_data_tm_P5P5.Read("../tau_decay_strange/light_mass_correction/mix_l2_s", "mes_contr_mix_l2_s_TM_P5P5", "P5P5", Sort_easy);
  //OS
  ls_ph_data_OS_VKVK.Read("../tau_decay_strange/light_mass_correction/mix_l2_s", "mes_contr_mix_l2_s_OS_VKVK", "VKVK", Sort_easy);
  ls_ph_data_OS_V0V0.Read("../tau_decay_strange/light_mass_correction/mix_l2_s", "mes_contr_mix_l2_s_OS_V0V0", "V0V0", Sort_easy);
  ls_ph_data_OS_AKAK.Read("../tau_decay_strange/light_mass_correction/mix_l2_s", "mes_contr_mix_l2_s_OS_AKAK", "AKAK", Sort_easy);
  ls_ph_data_OS_A0A0.Read("../tau_decay_strange/light_mass_correction/mix_l2_s", "mes_contr_mix_l2_s_OS_A0A0", "A0A0", Sort_easy);

  


  data_t ls_uni_data_tm_VKVK, ls_uni_data_tm_V0V0, ls_uni_data_tm_AKAK, ls_uni_data_tm_A0A0;
  data_t ls_uni_data_OS_VKVK, ls_uni_data_OS_V0V0, ls_uni_data_OS_AKAK, ls_uni_data_OS_A0A0;
  data_t ls_uni_data_tm_P5P5;



  //light

  //tm
  ls_uni_data_tm_VKVK.Read("../tau_decay_strange/light_mass_correction/mix_l1_s", "mes_contr_mix_l1_s_TM_VKVK", "VKVK", Sort_easy);
  ls_uni_data_tm_V0V0.Read("../tau_decay_strange/light_mass_correction/mix_l1_s", "mes_contr_mix_l1_s_TM_V0V0", "V0V0", Sort_easy);
  ls_uni_data_tm_AKAK.Read("../tau_decay_strange/light_mass_correction/mix_l1_s", "mes_contr_mix_l1_s_TM_AKAK", "AKAK", Sort_easy);
  ls_uni_data_tm_A0A0.Read("../tau_decay_strange/light_mass_correction/mix_l1_s", "mes_contr_mix_l1_s_TM_A0A0", "A0A0", Sort_easy);
  ls_uni_data_tm_P5P5.Read("../tau_decay_strange/light_mass_correction/mix_l1_s", "mes_contr_mix_l1_s_TM_P5P5", "P5P5", Sort_easy);
  //OS
  ls_uni_data_OS_VKVK.Read("../tau_decay_strange/light_mass_correction/mix_l1_s", "mes_contr_mix_l1_s_OS_VKVK", "VKVK", Sort_easy);
  ls_uni_data_OS_V0V0.Read("../tau_decay_strange/light_mass_correction/mix_l1_s", "mes_contr_mix_l1_s_OS_V0V0", "V0V0", Sort_easy);
  ls_uni_data_OS_AKAK.Read("../tau_decay_strange/light_mass_correction/mix_l1_s", "mes_contr_mix_l1_s_OS_AKAK", "AKAK", Sort_easy);
  ls_uni_data_OS_A0A0.Read("../tau_decay_strange/light_mass_correction/mix_l1_s", "mes_contr_mix_l1_s_OS_A0A0", "A0A0", Sort_easy);

    


  int Nens = ll_data_tm_P5P5.size;

  int Nens_mcorr= ls_ph_data_OS_A0A0.size;

  cout<<"Nens_mcorr: "<<Nens_mcorr<<endl;

  boost::filesystem::create_directory("../data/tau_decay_strange");
  
  



  //############################################################################################
  //generate fake jack_distr for lattice spacing a_A a_B, a_C, a_D and RENORMALIZATION CONSTANT
  GaussianMersenne GM(36551294);
  LatticeInfo a_info;
  distr_t a_A(UseJack), a_B(UseJack), a_C(UseJack), a_D(UseJack), a_Z(UseJack), a_E(UseJack);
  distr_t ZV_A(UseJack), ZV_B(UseJack), ZV_C(UseJack), ZV_D(UseJack), ZV_Z(UseJack), ZV_E(UseJack);
  distr_t ZA_A(UseJack), ZA_B(UseJack), ZA_C(UseJack), ZA_D(UseJack), ZA_Z(UseJack) , ZA_E(UseJack);
  distr_t ZV_C112(UseJack), ZA_C112(UseJack), ZV_B96(UseJack), ZA_B96(UseJack);
  
  double a_A_ave, a_A_err, a_B_ave, a_B_err, a_C_ave, a_C_err, a_D_ave, a_D_err, a_Z_ave, a_Z_err, a_E_ave, a_E_err;
  double ZV_A_ave, ZV_A_err, ZV_B_ave, ZV_B_err, ZV_C_ave, ZV_C_err, ZV_D_ave, ZV_D_err, ZV_Z_ave, ZV_Z_err, ZV_E_ave, ZV_E_err ;
  double ZA_A_ave, ZA_A_err, ZA_B_ave, ZA_B_err, ZA_C_ave, ZA_C_err, ZA_D_ave, ZA_D_err, ZA_Z_ave, ZA_Z_err, ZA_E_ave, ZA_E_err;

  double ZA_C112_ave, ZA_C112_err, ZV_C112_ave, ZV_C112_err;
  double ZA_B96_ave, ZA_B96_err, ZV_B96_ave, ZV_B96_err;

  distr_t ZV_A_boot(false), ZV_B_boot(false), ZV_C_boot(false), ZV_D_boot(false), ZV_E_boot(false);
  distr_t ZA_A_boot(false), ZA_B_boot(false), ZA_C_boot(false), ZA_D_boot(false), ZA_E_boot(false);
  
  a_info.LatInfo_new_ens("cA211a.53.24");
  a_A_ave= a_info.a_from_afp_FLAG;
  a_A_err= a_info.a_from_afp_FLAG_err;
  ZA_A_ave = a_info.Za_WI_FLAG;
  ZA_A_err = a_info.Za_WI_FLAG_err;
  ZV_A_ave = a_info.Zv_WI_FLAG;
  ZV_A_err = a_info.Zv_WI_FLAG_err;
  a_info.LatInfo_new_ens("cB211b.072.64");
  a_B_ave= a_info.a_from_afp_FLAG;
  a_B_err= a_info.a_from_afp_FLAG_err;
  ZA_B_ave = a_info.Za_WI_FLAG;
  ZA_B_err = a_info.Za_WI_FLAG_err;
  ZV_B_ave = a_info.Zv_WI_FLAG;
  ZV_B_err = a_info.Zv_WI_FLAG_err;
  a_info.LatInfo_new_ens("cB211b.072.96");
  ZA_B96_ave = a_info.Za_WI_FLAG;
  ZA_B96_err = a_info.Za_WI_FLAG_err;
  ZV_B96_ave = a_info.Zv_WI_FLAG;
  ZV_B96_err = a_info.Zv_WI_FLAG_err;
  a_info.LatInfo_new_ens("cC211a.06.80");
  a_C_ave= a_info.a_from_afp_FLAG;
  a_C_err= a_info.a_from_afp_FLAG_err;
  ZA_C_ave = a_info.Za_WI_FLAG;
  ZA_C_err = a_info.Za_WI_FLAG_err;
  ZV_C_ave = a_info.Zv_WI_FLAG;
  ZV_C_err = a_info.Zv_WI_FLAG_err;
  a_info.LatInfo_new_ens("cC211a.06.112");
  ZA_C112_ave = a_info.Za_WI_FLAG;
  ZA_C112_err = a_info.Za_WI_FLAG_err;
  ZV_C112_ave = a_info.Zv_WI_FLAG;
  ZV_C112_err = a_info.Zv_WI_FLAG_err;
  a_info.LatInfo_new_ens("cD211a.054.96");
  a_D_ave= a_info.a_from_afp_FLAG;
  a_D_err= a_info.a_from_afp_FLAG_err;
  ZA_D_ave = a_info.Za_WI_FLAG;
  ZA_D_err = a_info.Za_WI_FLAG_err;
  ZV_D_ave = a_info.Zv_WI_FLAG;
  ZV_D_err = a_info.Zv_WI_FLAG_err;
  a_info.LatInfo_new_ens("cZ211a.077.64");
  a_Z_ave= a_info.a_from_afp_FLAG;
  a_Z_err= a_info.a_from_afp_FLAG_err;
  ZA_Z_ave = a_info.Za_WI_FLAG;
  ZA_Z_err = a_info.Za_WI_FLAG_err;
  ZV_Z_ave = a_info.Zv_WI_FLAG;
  ZV_Z_err = a_info.Zv_WI_FLAG_err;

  a_info.LatInfo_new_ens("cE211a.044.112");
  a_E_ave= a_info.a_from_afp_FLAG;
  a_E_err= a_info.a_from_afp_FLAG_err;
  ZA_E_ave = a_info.Za_WI_FLAG;
  ZA_E_err = a_info.Za_WI_FLAG_err;
  ZV_E_ave = a_info.Zv_WI_FLAG;
  ZV_E_err = a_info.Zv_WI_FLAG_err;
  
  
  if(UseJack)  { for(int ijack=0;ijack<Njacks;ijack++) {
      a_A.distr.push_back( fm_to_inv_Gev*( a_A_ave + GM()*a_A_err*(1.0/sqrt(Njacks-1.0))));
      a_B.distr.push_back( fm_to_inv_Gev*( a_B_ave + GM()*a_B_err*(1.0/sqrt(Njacks-1.0))));
      a_C.distr.push_back( fm_to_inv_Gev*( a_C_ave + GM()*a_C_err*(1.0/sqrt(Njacks-1.0))));
      a_D.distr.push_back( fm_to_inv_Gev*( a_D_ave + GM()*a_D_err*(1.0/sqrt(Njacks-1.0))));
      a_Z.distr.push_back( fm_to_inv_Gev*( a_Z_ave + GM()*a_Z_err*(1.0/sqrt(Njacks-1.0))));
      a_E.distr.push_back( fm_to_inv_Gev*( a_E_ave + GM()*a_E_err*(1.0/sqrt(Njacks-1.0))));

      
      ZA_A.distr.push_back(  ZA_A_ave + GM()*ZA_A_err*(1.0/sqrt(Njacks -1.0)));
      ZV_A.distr.push_back(  ZV_A_ave + GM()*ZV_A_err*(1.0/sqrt(Njacks -1.0)));
      ZA_B.distr.push_back(  ZA_B_ave + GM()*ZA_B_err*(1.0/sqrt(Njacks -1.0)));
      ZV_B.distr.push_back(  ZV_B_ave + GM()*ZV_B_err*(1.0/sqrt(Njacks -1.0)));
      ZA_C.distr.push_back(  ZA_C_ave + GM()*ZA_C_err*(1.0/sqrt(Njacks -1.0)));
      ZV_C.distr.push_back(  ZV_C_ave + GM()*ZV_C_err*(1.0/sqrt(Njacks -1.0)));
      ZA_D.distr.push_back(  ZA_D_ave + GM()*ZA_D_err*(1.0/sqrt(Njacks -1.0)));
      ZV_D.distr.push_back(  ZV_D_ave + GM()*ZV_D_err*(1.0/sqrt(Njacks -1.0)));
      ZA_Z.distr.push_back(  ZA_Z_ave + GM()*ZA_Z_err*(1.0/sqrt(Njacks -1.0)));
      ZV_Z.distr.push_back(  ZV_Z_ave + GM()*ZV_Z_err*(1.0/sqrt(Njacks -1.0)));
      ZA_E.distr.push_back(  ZA_E_ave + GM()*ZA_E_err*(1.0/sqrt(Njacks -1.0)));
      ZV_E.distr.push_back(  ZV_E_ave + GM()*ZV_E_err*(1.0/sqrt(Njacks -1.0)));

      ZA_B96.distr.push_back(  ZA_B96_ave + GM()*ZA_B96_err*(1.0/sqrt(Njacks -1.0)));
      ZV_B96.distr.push_back(  ZV_B96_ave + GM()*ZV_B96_err*(1.0/sqrt(Njacks -1.0)));
      ZA_C112.distr.push_back(  ZA_C112_ave + GM()*ZA_C112_err*(1.0/sqrt(Njacks -1.0)));
      ZV_C112.distr.push_back(  ZV_C112_ave + GM()*ZV_C112_err*(1.0/sqrt(Njacks -1.0)));
      
    }
  }
  else {
    for (int iboot=0; iboot<Nboots;iboot++) {
      a_A.distr.push_back( fm_to_inv_Gev*( a_A_ave + GM()*a_A_err));
      a_B.distr.push_back( fm_to_inv_Gev*( a_B_ave + GM()*a_B_err));
      a_C.distr.push_back( fm_to_inv_Gev*( a_C_ave + GM()*a_C_err));
      a_D.distr.push_back( fm_to_inv_Gev*( a_D_ave + GM()*a_D_err));
      a_E.distr.push_back( fm_to_inv_Gev*( a_E_ave + GM()*a_E_err));
      ZA_A.distr.push_back(  ZA_A_ave + GM()*ZA_A_err);
      ZV_A.distr.push_back(  ZV_A_ave + GM()*ZV_A_err);
      ZA_B.distr.push_back(  ZA_B_ave + GM()*ZA_B_err);
      ZV_B.distr.push_back(  ZV_B_ave + GM()*ZV_B_err);
      ZA_C.distr.push_back(  ZA_C_ave + GM()*ZA_C_err);
      ZV_C.distr.push_back(  ZV_C_ave + GM()*ZV_C_err);
      ZA_D.distr.push_back(  ZA_D_ave + GM()*ZA_D_err);
      ZV_D.distr.push_back(  ZV_D_ave + GM()*ZV_D_err);
      ZA_Z.distr.push_back(  ZA_Z_ave + GM()*ZA_Z_err);
      ZV_Z.distr.push_back(  ZV_Z_ave + GM()*ZV_Z_err);
      ZA_E.distr.push_back(  ZA_E_ave + GM()*ZA_E_err);
      ZV_E.distr.push_back(  ZV_E_ave + GM()*ZV_E_err);
      
    }
  }

  for(int iboot=0;iboot<Nboots;iboot++) {
    ZA_A_boot.distr.push_back( ZA_A_ave + GM()*ZA_A_err);
    ZV_A_boot.distr.push_back(  ZV_A_ave + GM()*ZV_A_err);
    ZA_B_boot.distr.push_back(  ZA_B_ave + GM()*ZA_B_err);
    ZV_B_boot.distr.push_back(  ZV_B_ave + GM()*ZV_B_err);
    ZA_C_boot.distr.push_back(  ZA_C_ave + GM()*ZA_C_err);
    ZV_C_boot.distr.push_back(  ZV_C_ave + GM()*ZV_C_err);
    ZA_D_boot.distr.push_back(  ZA_D_ave + GM()*ZA_D_err);
    ZV_D_boot.distr.push_back(  ZV_D_ave + GM()*ZV_D_err);
    ZA_E_boot.distr.push_back(  ZA_E_ave + GM()*ZA_E_err);
    ZV_E_boot.distr.push_back(  ZV_E_ave + GM()*ZV_E_err);
    


  }


  //############################################################################################




 
  //resize vector with systematic errors
  vector<distr_t_list> Br_tau_tm, Br_tau_OS;
  vector<distr_t_list> Br_A0_tau_tm, Br_V0_tau_tm, Br_Aii_tau_tm, Br_Vii_tau_tm;
  vector<distr_t_list> Br_A0_tau_OS, Br_V0_tau_OS, Br_Aii_tau_OS, Br_Vii_tau_OS;

  vector<distr_t_list> Br_T_tau_tm, Br_L_tau_tm, Br_T_tau_OS, Br_L_tau_OS;

  VVfloat pull_per_ens_tm_A0(Nens);
  VVfloat pull_per_ens_OS_A0(Nens);
  VVfloat pull_per_ens_tm_V0(Nens);
  VVfloat pull_per_ens_OS_V0(Nens);

  VVfloat pull_per_ens_tm_Aii(Nens);
  VVfloat pull_per_ens_OS_Aii(Nens);
  VVfloat pull_per_ens_tm_Vii(Nens);
  VVfloat pull_per_ens_OS_Vii(Nens);

  VVfloat pull_per_ens_tm_T(Nens);
  VVfloat pull_per_ens_OS_T(Nens);
  VVfloat pull_per_ens_tm_L(Nens);
  VVfloat pull_per_ens_OS_L(Nens);
   
  VVfloat syst_per_ens_tm_A0(Nens);
  VVfloat syst_per_ens_tm_V0(Nens);
  VVfloat syst_per_ens_tm_Ak(Nens);
  VVfloat syst_per_ens_tm_Vk(Nens);
  VVfloat syst_per_ens_OS_A0(Nens);
  VVfloat syst_per_ens_OS_V0(Nens);
  VVfloat syst_per_ens_OS_Ak(Nens);
  VVfloat syst_per_ens_OS_Vk(Nens);
  VVfloat syst_per_ens_tm(Nens);
  VVfloat syst_per_ens_OS(Nens);

  VVfloat syst_per_ens_tm_T(Nens);
  VVfloat syst_per_ens_tm_L(Nens);

  VVfloat syst_per_ens_OS_T(Nens);
  VVfloat syst_per_ens_OS_L(Nens);


  for(int iens=0; iens<Nens; iens++) {
    Br_tau_tm.emplace_back( UseJack, sigma_list_strange.size());
    Br_A0_tau_tm.emplace_back( UseJack, sigma_list_strange.size());
    Br_V0_tau_tm.emplace_back( UseJack, sigma_list_strange.size());
    Br_Aii_tau_tm.emplace_back( UseJack, sigma_list_strange.size());
    Br_Vii_tau_tm.emplace_back( UseJack, sigma_list_strange.size());

    Br_T_tau_tm.emplace_back( UseJack, sigma_list_strange.size());
    Br_L_tau_tm.emplace_back( UseJack, sigma_list_strange.size());
    
    
   
    
    Br_tau_OS.emplace_back( UseJack, sigma_list_strange.size());
    Br_A0_tau_OS.emplace_back( UseJack, sigma_list_strange.size());
    Br_V0_tau_OS.emplace_back( UseJack, sigma_list_strange.size());
    Br_Aii_tau_OS.emplace_back( UseJack, sigma_list_strange.size());
    Br_Vii_tau_OS.emplace_back( UseJack, sigma_list_strange.size());

    Br_T_tau_OS.emplace_back( UseJack, sigma_list_strange.size());
    Br_L_tau_OS.emplace_back( UseJack, sigma_list_strange.size());

    syst_per_ens_tm_A0[iens].resize(sigma_list_strange.size());
    syst_per_ens_tm_V0[iens].resize(sigma_list_strange.size());
    syst_per_ens_tm_Ak[iens].resize(sigma_list_strange.size());
    syst_per_ens_tm_Vk[iens].resize(sigma_list_strange.size());
    syst_per_ens_tm[iens].resize(sigma_list_strange.size());

    syst_per_ens_tm_T[iens].resize(sigma_list_strange.size());
    syst_per_ens_tm_L[iens].resize(sigma_list_strange.size());

    syst_per_ens_OS_A0[iens].resize(sigma_list_strange.size());
    syst_per_ens_OS_V0[iens].resize(sigma_list_strange.size());
    syst_per_ens_OS_Ak[iens].resize(sigma_list_strange.size());
    syst_per_ens_OS_Vk[iens].resize(sigma_list_strange.size());
    syst_per_ens_OS[iens].resize(sigma_list_strange.size());

    syst_per_ens_OS_T[iens].resize(sigma_list_strange.size());
    syst_per_ens_OS_L[iens].resize(sigma_list_strange.size());

    pull_per_ens_tm_A0[iens].resize(sigma_list_strange.size());
    pull_per_ens_OS_A0[iens].resize(sigma_list_strange.size());

    pull_per_ens_tm_V0[iens].resize(sigma_list_strange.size());
    pull_per_ens_OS_V0[iens].resize(sigma_list_strange.size());
    
    pull_per_ens_tm_Aii[iens].resize(sigma_list_strange.size());
    pull_per_ens_OS_Aii[iens].resize(sigma_list_strange.size());
    
    pull_per_ens_tm_Vii[iens].resize(sigma_list_strange.size());
    pull_per_ens_OS_Vii[iens].resize(sigma_list_strange.size());
    
    pull_per_ens_tm_T[iens].resize(sigma_list_strange.size());
    pull_per_ens_OS_T[iens].resize(sigma_list_strange.size());
    
    pull_per_ens_tm_L[iens].resize(sigma_list_strange.size());
    pull_per_ens_OS_L[iens].resize(sigma_list_strange.size());

    
    

  }

 
  distr_t Mk_iso(UseJack);
  for(int ijack=0;ijack<Njacks;ijack++) { Mk_iso.distr.push_back( m_kappa + GM()*m_kappa_err/sqrt(Njacks-1.0));}


  //pion masses in lattice units
  double amp_B_ave= 0.0565313; double amp_B_err= 1.438e-05;
  double amp_C_ave= 0.0472193; double amp_C_err= 3.45183e-05;
  double amp_D_ave= 0.0406214; double amp_D_err= 2.94047e-05;
  double amp_E_ave= 0.0338185; double amp_E_err= 3.18799e-05;

  distr_t aMp_B(UseJack), aMp_C(UseJack), aMp_D(UseJack), aMp_E(UseJack);

  for(int ijack=0;ijack<Njacks;ijack++ ) {

    aMp_B.distr.push_back( amp_B_ave + GM()*amp_B_err/sqrt(Njacks-1.0));
    aMp_C.distr.push_back( amp_C_ave + GM()*amp_C_err/sqrt(Njacks-1.0));
    aMp_D.distr.push_back( amp_D_ave + GM()*amp_D_err/sqrt(Njacks-1.0));
    aMp_E.distr.push_back( amp_E_ave + GM()*amp_E_err/sqrt(Njacks-1.0));

  }
  
 
  
  distr_t_list FK_list(UseJack); distr_t_list a_distr_list(UseJack);

  //loop over the ensembles
  for(int iens=0; iens<Nens;iens++) {


    boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/masses");
    boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/masses/"+ll_data_tm_P5P5.Tag[iens]);
    
     CorrAnalysis Corr(UseJack, Njacks,Nboots);
     Corr.Nt = ls_data_tm_VKVK.nrows[iens];

   
     //effective masses
     //tm
     //pseudoscalar
     distr_t_list M_pi_tm = Corr.effective_mass_t( ll_data_tm_P5P5.col(0)[iens], "../data/tau_decay/"+Tag_reco_type+"/strange/masses/"+ll_data_tm_P5P5.Tag[iens]+"/eff_mass_Mpi_tm");
     distr_t_list M_etas_tm = Corr.effective_mass_t( ss_data_tm_P5P5.col(0)[iens], "../data/tau_decay/"+Tag_reco_type+"/strange/masses/"+ll_data_tm_P5P5.Tag[iens]+"/eff_mass_etas_tm");
     distr_t_list M_etas_H_tm = Corr.effective_mass_t( ss_H_data_tm_P5P5.col(0)[iens], "../data/tau_decay/"+Tag_reco_type+"/strange/masses/"+ll_data_tm_P5P5.Tag[iens]+"/eff_mass_etas_H_tm");
     distr_t_list M_K_tm = Corr.effective_mass_t( ls_data_tm_P5P5.col(0)[iens], "../data/tau_decay/"+Tag_reco_type+"/strange/masses/"+ll_data_tm_P5P5.Tag[iens]+"/eff_mass_K_tm");
     distr_t_list M_K_H_tm = Corr.effective_mass_t( ls_H_data_tm_P5P5.col(0)[iens], "../data/tau_decay/"+Tag_reco_type+"/strange/masses/"+ll_data_tm_P5P5.Tag[iens]+"/eff_mass_K_H_tm");
     
     //vector
     distr_t_list M_rho_tm = Corr.effective_mass_t( ll_data_tm_VKVK.col(0)[iens], "../data/tau_decay/"+Tag_reco_type+"/strange/masses/"+ll_data_tm_P5P5.Tag[iens]+"/eff_mass_Mrho_tm");
     distr_t_list M_phi_tm = Corr.effective_mass_t( ss_data_tm_VKVK.col(0)[iens], "../data/tau_decay/"+Tag_reco_type+"/strange/masses/"+ll_data_tm_P5P5.Tag[iens]+"/eff_mass_phi_tm");
     distr_t_list M_phi_H_tm = Corr.effective_mass_t( ss_H_data_tm_VKVK.col(0)[iens], "../data/tau_decay/"+Tag_reco_type+"/strange/masses/"+ll_data_tm_P5P5.Tag[iens]+"/eff_mass_phi_H_tm");
     distr_t_list M_Kstar_tm = Corr.effective_mass_t( ls_data_tm_VKVK.col(0)[iens], "../data/tau_decay/"+Tag_reco_type+"/strange/masses/"+ll_data_tm_P5P5.Tag[iens]+"/eff_mass_K*_tm");
     distr_t_list M_Kstar_H_tm = Corr.effective_mass_t( ls_H_data_tm_VKVK.col(0)[iens], "../data/tau_decay/"+Tag_reco_type+"/strange/masses/"+ll_data_tm_P5P5.Tag[iens]+"/eff_mass_K*_H_tm");
     //axial-vector
     distr_t_list M_a1_tm = Corr.effective_mass_t( ll_data_tm_AKAK.col(0)[iens], "../data/tau_decay/"+Tag_reco_type+"/strange/masses/"+ll_data_tm_P5P5.Tag[iens]+"/eff_mass_a1_tm");
     distr_t_list M_f1_tm = Corr.effective_mass_t( ss_data_tm_AKAK.col(0)[iens], "../data/tau_decay/"+Tag_reco_type+"/strange/masses/"+ll_data_tm_P5P5.Tag[iens]+"/eff_mass_f1_tm");
     distr_t_list M_f1_H_tm = Corr.effective_mass_t( ss_H_data_tm_AKAK.col(0)[iens], "../data/tau_decay/"+Tag_reco_type+"/strange/masses/"+ll_data_tm_P5P5.Tag[iens]+"/eff_mass_f1_H_tm");
     distr_t_list M_K1_tm = Corr.effective_mass_t( ls_data_tm_AKAK.col(0)[iens], "../data/tau_decay/"+Tag_reco_type+"/strange/masses/"+ll_data_tm_P5P5.Tag[iens]+"/eff_mass_K1_tm");
     distr_t_list M_K1_H_tm = Corr.effective_mass_t( ls_H_data_tm_AKAK.col(0)[iens], "../data/tau_decay/"+Tag_reco_type+"/strange/masses/"+ll_data_tm_P5P5.Tag[iens]+"/eff_mass_K1_H_tm");
     //scalar
     distr_t_list M_f0_tm = Corr.effective_mass_t( ll_data_tm_S0S0.col(0)[iens], "../data/tau_decay/"+Tag_reco_type+"/strange/masses/"+ll_data_tm_P5P5.Tag[iens]+"/eff_mass_f0_tm");
     distr_t_list M_fs_tm = Corr.effective_mass_t( ss_data_tm_S0S0.col(0)[iens], "../data/tau_decay/"+Tag_reco_type+"/strange/masses/"+ll_data_tm_P5P5.Tag[iens]+"/eff_mass_fs_tm");
     distr_t_list M_fs_H_tm = Corr.effective_mass_t( ss_H_data_tm_S0S0.col(0)[iens], "../data/tau_decay/"+Tag_reco_type+"/strange/masses/"+ll_data_tm_P5P5.Tag[iens]+"/eff_mass_fs_H_tm");
     distr_t_list M_Ks0_tm = Corr.effective_mass_t( ls_data_tm_S0S0.col(0)[iens], "../data/tau_decay/"+Tag_reco_type+"/strange/masses/"+ll_data_tm_P5P5.Tag[iens]+"/eff_mass_Ks0_tm");
     distr_t_list M_Ks0_H_tm = Corr.effective_mass_t( ls_H_data_tm_S0S0.col(0)[iens], "../data/tau_decay/"+Tag_reco_type+"/strange/masses/"+ll_data_tm_P5P5.Tag[iens]+"/eff_mass_Ks0_H_tm");

     //OS
     //pseudoscalar
     distr_t_list M_pi_OS = Corr.effective_mass_t( ll_data_OS_P5P5.col(0)[iens], "../data/tau_decay/"+Tag_reco_type+"/strange/masses/"+ll_data_OS_P5P5.Tag[iens]+"/eff_mass_Mpi_OS");
     distr_t_list M_etas_OS = Corr.effective_mass_t( ss_data_OS_P5P5.col(0)[iens], "../data/tau_decay/"+Tag_reco_type+"/strange/masses/"+ll_data_OS_P5P5.Tag[iens]+"/eff_mass_etas_OS");
     distr_t_list M_etas_H_OS = Corr.effective_mass_t( ss_H_data_OS_P5P5.col(0)[iens], "../data/tau_decay/"+Tag_reco_type+"/strange/masses/"+ll_data_OS_P5P5.Tag[iens]+"/eff_mass_etas_H_OS");
     distr_t_list M_K_OS = Corr.effective_mass_t( ls_data_OS_P5P5.col(0)[iens], "../data/tau_decay/"+Tag_reco_type+"/strange/masses/"+ll_data_OS_P5P5.Tag[iens]+"/eff_mass_K_OS");
     distr_t_list M_K_H_OS = Corr.effective_mass_t( ls_H_data_OS_P5P5.col(0)[iens], "../data/tau_decay/"+Tag_reco_type+"/strange/masses/"+ll_data_OS_P5P5.Tag[iens]+"/eff_mass_K_H_OS");
     
     //vector
     distr_t_list M_rho_OS = Corr.effective_mass_t( ll_data_OS_VKVK.col(0)[iens], "../data/tau_decay/"+Tag_reco_type+"/strange/masses/"+ll_data_OS_P5P5.Tag[iens]+"/eff_mass_Mrho_OS");
     distr_t_list M_phi_OS = Corr.effective_mass_t( ss_data_OS_VKVK.col(0)[iens], "../data/tau_decay/"+Tag_reco_type+"/strange/masses/"+ll_data_OS_P5P5.Tag[iens]+"/eff_mass_phi_OS");
     distr_t_list M_phi_H_OS = Corr.effective_mass_t( ss_H_data_OS_VKVK.col(0)[iens], "../data/tau_decay/"+Tag_reco_type+"/strange/masses/"+ll_data_OS_P5P5.Tag[iens]+"/eff_mass_phi_H_OS");
     distr_t_list M_Kstar_OS = Corr.effective_mass_t( ls_data_OS_VKVK.col(0)[iens], "../data/tau_decay/"+Tag_reco_type+"/strange/masses/"+ll_data_OS_P5P5.Tag[iens]+"/eff_mass_K*_OS");
     distr_t_list M_Kstar_H_OS = Corr.effective_mass_t( ls_H_data_OS_VKVK.col(0)[iens], "../data/tau_decay/"+Tag_reco_type+"/strange/masses/"+ll_data_OS_P5P5.Tag[iens]+"/eff_mass_K*_H_OS");
     //axial-vector
     distr_t_list M_a1_OS = Corr.effective_mass_t( ll_data_OS_AKAK.col(0)[iens], "../data/tau_decay/"+Tag_reco_type+"/strange/masses/"+ll_data_OS_P5P5.Tag[iens]+"/eff_mass_a1_OS");
     distr_t_list M_f1_OS = Corr.effective_mass_t( ss_data_OS_AKAK.col(0)[iens], "../data/tau_decay/"+Tag_reco_type+"/strange/masses/"+ll_data_OS_P5P5.Tag[iens]+"/eff_mass_f1_OS");
     distr_t_list M_f1_H_OS = Corr.effective_mass_t( ss_H_data_OS_AKAK.col(0)[iens], "../data/tau_decay/"+Tag_reco_type+"/strange/masses/"+ll_data_OS_P5P5.Tag[iens]+"/eff_mass_f1_H_OS");
     distr_t_list M_K1_OS = Corr.effective_mass_t( ls_data_OS_AKAK.col(0)[iens], "../data/tau_decay/"+Tag_reco_type+"/strange/masses/"+ll_data_OS_P5P5.Tag[iens]+"/eff_mass_K1_OS");
     distr_t_list M_K1_H_OS = Corr.effective_mass_t( ls_H_data_OS_AKAK.col(0)[iens], "../data/tau_decay/"+Tag_reco_type+"/strange/masses/"+ll_data_OS_P5P5.Tag[iens]+"/eff_mass_K1_H_OS");
     //scalar
     distr_t_list M_f0_OS = Corr.effective_mass_t( ll_data_OS_S0S0.col(0)[iens], "../data/tau_decay/"+Tag_reco_type+"/strange/masses/"+ll_data_OS_P5P5.Tag[iens]+"/eff_mass_f0_OS");
     distr_t_list M_fs_OS = Corr.effective_mass_t( ss_data_OS_S0S0.col(0)[iens], "../data/tau_decay/"+Tag_reco_type+"/strange/masses/"+ll_data_OS_P5P5.Tag[iens]+"/eff_mass_fs_OS");
     distr_t_list M_fs_H_OS = Corr.effective_mass_t( ss_H_data_OS_S0S0.col(0)[iens], "../data/tau_decay/"+Tag_reco_type+"/strange/masses/"+ll_data_OS_P5P5.Tag[iens]+"/eff_mass_fs_H_OS");
     distr_t_list M_Ks0_OS = Corr.effective_mass_t( ls_data_OS_S0S0.col(0)[iens], "../data/tau_decay/"+Tag_reco_type+"/strange/masses/"+ll_data_OS_P5P5.Tag[iens]+"/eff_mass_Ks0_OS");
     distr_t_list M_Ks0_H_OS = Corr.effective_mass_t( ls_H_data_OS_S0S0.col(0)[iens], "../data/tau_decay/"+Tag_reco_type+"/strange/masses/"+ll_data_OS_P5P5.Tag[iens]+"/eff_mass_Ks0_H_OS");


   
    

     
   
     
     LatticeInfo L_info;
     
     L_info.LatInfo_new_ens(ls_data_tm_VKVK.Tag[iens]);
     
     double aml= L_info.ml;
     double ams1= L_info.ms_L_new;
     double ams2= L_info.ms_M_new;
     
     
       
     //CorrAnalysis Corr_block_1(0, ls_data_tm_VKVK.Nconfs[iens],Nboots, iens);
     CorrAnalysis Corr_block_1(1, ls_data_tm_VKVK.Nconfs[iens], Nboots);
     Corr_block_1.Nt= ls_data_tm_VKVK.nrows[iens];
     int T = Corr.Nt;
     
     
   
     
     //get lattice spacing
     distr_t a_distr(UseJack);
     distr_t Zv(UseJack), Za(UseJack);
     if(ls_data_tm_VKVK.Tag[iens].substr(1,1)=="B") {a_distr=a_B;
       if(ls_data_tm_VKVK.Tag[iens] == "cB211b.072.64") {  Zv = ZV_B; Za = ZA_B; }
       else if( ls_data_tm_VKVK.Tag[iens] == "cB211b.072.96")  { Zv= ZV_B96; Za =ZA_B96; }
       else crash("Ens: "+ls_data_tm_VKVK.Tag[iens]+" not found");
     }
     else if(ls_data_tm_VKVK.Tag[iens].substr(1,1)=="C") {a_distr=a_C;
       if(ls_data_tm_VKVK.Tag[iens] == "cC211a.06.80") {  Zv = ZV_C; Za = ZA_C; }
       else if ( ls_data_tm_VKVK.Tag[iens] == "cC211a.06.112")  { Zv= ZV_C112; Za =ZA_C112; }
       else crash("Ens: "+ls_data_tm_VKVK.Tag[iens]+" not found");
     }
     else if(ls_data_tm_VKVK.Tag[iens].substr(1,1)=="D") {a_distr=a_D; Zv = ZV_D; Za = ZA_D; }
     else if(ls_data_tm_VKVK.Tag[iens].substr(1,1)=="E") {a_distr=a_E; Zv = ZV_E; Za = ZA_E; }
     else crash("lattice spacing distribution for Ens: "+ls_data_tm_VKVK.Tag[iens]+" not found");


     a_distr_list.distr_list.push_back( a_distr/fm_to_inv_Gev);


     cout<<"RC-used Analyzing Ensemble: "<<ls_data_tm_VKVK.Tag[iens]<<endl;
     cout<<"Zv: "<<Zv.ave()<<" +- "<<Zv.err()<<endl;
     cout<<"Za: "<<Za.ave()<<" +- "<<Za.err()<<endl;
  

   
    
 
    //############# LIGHTER MASS ################//
    
    //tm
    distr_t_list Vk_tm_distr, Ak_tm_distr, A0_tm_distr, V0_tm_distr;
    //tm block1
    distr_t_list Vk_tm_block_1_distr, Ak_tm_block_1_distr, A0_tm_block_1_distr, V0_tm_block_1_distr;
    //OS
    distr_t_list Vk_OS_distr, Ak_OS_distr, A0_OS_distr, V0_OS_distr;
    //OS block1
    distr_t_list Vk_OS_block_1_distr,  Ak_OS_block_1_distr, A0_OS_block_1_distr, V0_OS_block_1_distr;
    distr_t_list P5_tm_distr;


    //tm
    distr_t_list T_tm_distr, L_tm_distr;
    //tm block1
    distr_t_list L_tm_block_1_distr, T_tm_block_1_distr;
    //OS
    distr_t_list T_OS_distr, L_OS_distr;
    //OS block1
    distr_t_list T_OS_block_1_distr,  L_OS_block_1_distr;

   
    //light-tm sector
    Vk_tm_distr = Corr.corr_t(ls_data_tm_VKVK.col(0)[iens], "../data/tau_decay/"+Tag_reco_type+"/strange/corr/Vk_tm_"+ls_data_tm_VKVK.Tag[iens]+".dat");
    Ak_tm_distr = Corr.corr_t(ls_data_tm_AKAK.col(0)[iens], "../data/tau_decay/"+Tag_reco_type+"/strange/corr/Ak_tm_"+ls_data_tm_VKVK.Tag[iens]+".dat");
    A0_tm_distr = Corr.corr_t(ls_data_tm_A0A0.col(0)[iens], "../data/tau_decay/"+Tag_reco_type+"/strange/corr/A0_tm_"+ls_data_tm_VKVK.Tag[iens]+".dat");
    V0_tm_distr = Corr.corr_t(ls_data_tm_V0V0.col(0)[iens], "../data/tau_decay/"+Tag_reco_type+"/strange/corr/V0_tm_"+ls_data_tm_VKVK.Tag[iens]+".dat");
    P5_tm_distr = Corr.corr_t(ls_data_tm_P5P5.col(0)[iens], "../data/tau_decay/"+Tag_reco_type+"/strange/corr/P5_tm_"+ls_data_tm_VKVK.Tag[iens]+".dat");

  

    //light-OS sector
    Vk_OS_distr = Corr.corr_t(ls_data_OS_VKVK.col(0)[iens], "../data/tau_decay/"+Tag_reco_type+"/strange/corr/Vk_OS_"+ls_data_tm_VKVK.Tag[iens]+".dat");
    Ak_OS_distr = Corr.corr_t(ls_data_OS_AKAK.col(0)[iens], "../data/tau_decay/"+Tag_reco_type+"/strange/corr/Ak_OS_"+ls_data_tm_VKVK.Tag[iens]+".dat");
    A0_OS_distr = Corr.corr_t(ls_data_OS_A0A0.col(0)[iens], "../data/tau_decay/"+Tag_reco_type+"/strange/corr/A0_OS_"+ls_data_tm_VKVK.Tag[iens]+".dat");
    V0_OS_distr = Corr.corr_t(ls_data_OS_V0V0.col(0)[iens], "../data/tau_decay/"+Tag_reco_type+"/strange/corr/V0_OS_"+ls_data_tm_VKVK.Tag[iens]+".dat");


    distr_t_list MK_OS_A0_distr = Corr.effective_mass_t(A0_OS_distr, "../data/tau_decay/"+Tag_reco_type+"/strange/masses/"+ll_data_OS_P5P5.Tag[iens]+"/eff_mass_K_OS_A0");


    
    
    //light-tm sector
    Vk_tm_block_1_distr = Corr_block_1.corr_t(ls_data_tm_VKVK.col(0)[iens], "");
    Ak_tm_block_1_distr = Corr_block_1.corr_t(ls_data_tm_AKAK.col(0)[iens], "");
    A0_tm_block_1_distr = Corr_block_1.corr_t(ls_data_tm_A0A0.col(0)[iens], "");
    V0_tm_block_1_distr = Corr_block_1.corr_t(ls_data_tm_V0V0.col(0)[iens], "");

  

    //light-OS sector
    Vk_OS_block_1_distr = Corr_block_1.corr_t(ls_data_OS_VKVK.col(0)[iens], "");
    Ak_OS_block_1_distr = Corr_block_1.corr_t(ls_data_OS_AKAK.col(0)[iens], "");
    A0_OS_block_1_distr = Corr_block_1.corr_t(ls_data_OS_A0A0.col(0)[iens], "");
    V0_OS_block_1_distr = Corr_block_1.corr_t(ls_data_OS_V0V0.col(0)[iens], "");

  

    

    //############# LIGHTER MASS ################//



    //############# HEAVIER MASS ################//
    
    //tm
    distr_t_list Vk_H_tm_distr, Ak_H_tm_distr, A0_H_tm_distr, V0_H_tm_distr;
    //tm block1
    distr_t_list Vk_H_tm_block_1_distr, Ak_H_tm_block_1_distr, A0_H_tm_block_1_distr, V0_H_tm_block_1_distr;
    //OS
    distr_t_list Vk_H_OS_distr, Ak_H_OS_distr, A0_H_OS_distr, V0_H_OS_distr;
    //OS block1
    distr_t_list Vk_H_OS_block_1_distr,  Ak_H_OS_block_1_distr, A0_H_OS_block_1_distr, V0_H_OS_block_1_distr;
    distr_t_list P5_H_tm_distr;


    //tm
    distr_t_list T_H_tm_distr, L_H_tm_distr;
    //tm block1
    distr_t_list L_H_tm_block_1_distr, T_H_tm_block_1_distr;
    //OS
    distr_t_list T_H_OS_distr, L_H_OS_distr;
    //OS block1
    distr_t_list T_H_OS_block_1_distr,  L_H_OS_block_1_distr;
   
    //light-tm sector
    Vk_H_tm_distr = Corr.corr_t(ls_H_data_tm_VKVK.col(0)[iens], "../data/tau_decay/"+Tag_reco_type+"/strange/corr/Vk_H_tm_"+ls_H_data_tm_VKVK.Tag[iens]+".dat");
    Ak_H_tm_distr = Corr.corr_t(ls_H_data_tm_AKAK.col(0)[iens], "../data/tau_decay/"+Tag_reco_type+"/strange/corr/Ak_H_tm_"+ls_H_data_tm_VKVK.Tag[iens]+".dat");
    A0_H_tm_distr = Corr.corr_t(ls_H_data_tm_A0A0.col(0)[iens], "../data/tau_decay/"+Tag_reco_type+"/strange/corr/A0_H_tm_"+ls_H_data_tm_VKVK.Tag[iens]+".dat");
    V0_H_tm_distr = Corr.corr_t(ls_H_data_tm_V0V0.col(0)[iens], "../data/tau_decay/"+Tag_reco_type+"/strange/corr/V0_H_tm_"+ls_H_data_tm_VKVK.Tag[iens]+".dat");
    P5_H_tm_distr = Corr.corr_t(ls_H_data_tm_P5P5.col(0)[iens], "../data/tau_decay/"+Tag_reco_type+"/strange/corr/P5_H_tm_"+ls_H_data_tm_VKVK.Tag[iens]+".dat");
    

    
   

    //light-OS sector
    Vk_H_OS_distr = Corr.corr_t(ls_H_data_OS_VKVK.col(0)[iens], "../data/tau_decay/"+Tag_reco_type+"/strange/corr/Vk_H_OS_"+ls_H_data_tm_VKVK.Tag[iens]+".dat");
    Ak_H_OS_distr = Corr.corr_t(ls_H_data_OS_AKAK.col(0)[iens], "../data/tau_decay/"+Tag_reco_type+"/strange/corr/Ak_H_OS_"+ls_H_data_tm_VKVK.Tag[iens]+".dat");
    A0_H_OS_distr = Corr.corr_t(ls_H_data_OS_A0A0.col(0)[iens], "../data/tau_decay/"+Tag_reco_type+"/strange/corr/A0_H_OS_"+ls_H_data_tm_VKVK.Tag[iens]+".dat");
    V0_H_OS_distr = Corr.corr_t(ls_H_data_OS_V0V0.col(0)[iens], "../data/tau_decay/"+Tag_reco_type+"/strange/corr/V0_H_OS_"+ls_H_data_tm_VKVK.Tag[iens]+".dat");

  

    //analyze data with Njacks=Nconfs
    //light-tm sector
    Vk_H_tm_block_1_distr = Corr_block_1.corr_t(ls_H_data_tm_VKVK.col(0)[iens], "");
    Ak_H_tm_block_1_distr = Corr_block_1.corr_t(ls_H_data_tm_AKAK.col(0)[iens], "");
    A0_H_tm_block_1_distr = Corr_block_1.corr_t(ls_H_data_tm_A0A0.col(0)[iens], "");
    V0_H_tm_block_1_distr = Corr_block_1.corr_t(ls_H_data_tm_V0V0.col(0)[iens], "");

   

    //light-OS sector
    Vk_H_OS_block_1_distr = Corr_block_1.corr_t(ls_H_data_OS_VKVK.col(0)[iens], "");
    Ak_H_OS_block_1_distr = Corr_block_1.corr_t(ls_H_data_OS_AKAK.col(0)[iens], "");
    A0_H_OS_block_1_distr = Corr_block_1.corr_t(ls_H_data_OS_A0A0.col(0)[iens], "");
    V0_H_OS_block_1_distr = Corr_block_1.corr_t(ls_H_data_OS_V0V0.col(0)[iens], "");

  

    //############# HEAVIER MASS ################//


     //light-quark mass corrections

     distr_t_list V0V0_tm_dm_distr(UseJack), VKVK_tm_dm_distr(UseJack),  A0A0_tm_dm_distr(UseJack), AKAK_tm_dm_distr(UseJack);

     distr_t_list V0V0_OS_dm_distr(UseJack), VKVK_OS_dm_distr(UseJack),  A0A0_OS_dm_distr(UseJack), AKAK_OS_dm_distr(UseJack);

     distr_t_list P5P5_tm_L1(UseJack), P5P5_tm_L2(UseJack);


     //light- quark mass corrections bootstrap

     distr_t_list V0V0_tm_block_1_dm_distr(false), VKVK_tm_block_1_dm_distr(false),  A0A0_tm_block_1_dm_distr(false), AKAK_tm_block_1_dm_distr(false);

     distr_t_list V0V0_OS_block_1_dm_distr(false), VKVK_OS_block_1_dm_distr(false),  A0A0_OS_block_1_dm_distr(false), AKAK_OS_block_1_dm_distr(false);

   
     bool read_dm=false;
     int iens_dm=-1;
     for(int j=0;j<Nens_mcorr;j++) {
       if( ls_ph_data_OS_A0A0.Tag[j] == ll_data_OS_P5P5.Tag[iens]) { read_dm=true; iens_dm=j;}
     }

     P5P5_tm_L1 = Corr.corr_t( ls_ph_data_tm_P5P5.col(0)[iens_dm], "../data/tau_decay/"+Tag_reco_type+"/strange/corr/P5P5_L1_tm_"+ls_H_data_tm_V0V0.Tag[iens]+".dat");
     P5P5_tm_L2 = Corr.corr_t( ls_uni_data_tm_P5P5.col(0)[iens_dm], "../data/tau_decay/"+Tag_reco_type+"/strange/corr/P5P5_L2_tm_"+ls_H_data_tm_V0V0.Tag[iens]+".dat");

     

     if(!read_dm) crash("Cannot find dm-corrections for ensemble: "+ll_data_OS_P5P5.Tag[iens]);

     if(read_dm) {

       cout<<"Appling dm corrections..."<<endl;
       
       V0V0_tm_dm_distr = Corr.corr_t( summ_master( ls_ph_data_tm_V0V0.col(0)[iens_dm],  Multiply_Vvector_by_scalar(ls_uni_data_tm_V0V0.col(0)[iens_dm], -1.0))  , "../data/tau_decay/"+Tag_reco_type+"/strange/corr/V0_dm_tm_"+ls_H_data_tm_V0V0.Tag[iens]+".dat");
       VKVK_tm_dm_distr = Corr.corr_t( summ_master( ls_ph_data_tm_VKVK.col(0)[iens_dm],  Multiply_Vvector_by_scalar(ls_uni_data_tm_VKVK.col(0)[iens_dm], -1.0))  , "../data/tau_decay/"+Tag_reco_type+"/strange/corr/VK_dm_tm_"+ls_H_data_tm_VKVK.Tag[iens]+".dat");
       A0A0_tm_dm_distr = Corr.corr_t( summ_master( ls_ph_data_tm_A0A0.col(0)[iens_dm],  Multiply_Vvector_by_scalar(ls_uni_data_tm_A0A0.col(0)[iens_dm], -1.0))  , "../data/tau_decay/"+Tag_reco_type+"/strange/corr/A0_dm_tm_"+ls_H_data_tm_A0A0.Tag[iens]+".dat");
       AKAK_tm_dm_distr = Corr.corr_t( summ_master( ls_ph_data_tm_AKAK.col(0)[iens_dm],  Multiply_Vvector_by_scalar(ls_uni_data_tm_AKAK.col(0)[iens_dm], -1.0))  , "../data/tau_decay/"+Tag_reco_type+"/strange/corr/AK_dm_tm_"+ls_H_data_tm_AKAK.Tag[iens]+".dat");


       V0V0_OS_dm_distr = Corr.corr_t( summ_master( ls_ph_data_OS_V0V0.col(0)[iens_dm],  Multiply_Vvector_by_scalar(ls_uni_data_OS_V0V0.col(0)[iens_dm], -1.0))  , "../data/tau_decay/"+Tag_reco_type+"/strange/corr/V0_dm_OS_"+ls_H_data_OS_V0V0.Tag[iens]+".dat");
       VKVK_OS_dm_distr = Corr.corr_t( summ_master( ls_ph_data_OS_VKVK.col(0)[iens_dm],  Multiply_Vvector_by_scalar(ls_uni_data_OS_VKVK.col(0)[iens_dm], -1.0))  , "../data/tau_decay/"+Tag_reco_type+"/strange/corr/VK_dm_OS_"+ls_H_data_OS_VKVK.Tag[iens]+".dat");
       A0A0_OS_dm_distr = Corr.corr_t( summ_master( ls_ph_data_OS_A0A0.col(0)[iens_dm],  Multiply_Vvector_by_scalar(ls_uni_data_OS_A0A0.col(0)[iens_dm], -1.0))  , "../data/tau_decay/"+Tag_reco_type+"/strange/corr/A0_dm_OS_"+ls_H_data_OS_A0A0.Tag[iens]+".dat");
       AKAK_OS_dm_distr = Corr.corr_t( summ_master( ls_ph_data_OS_AKAK.col(0)[iens_dm],  Multiply_Vvector_by_scalar(ls_uni_data_OS_AKAK.col(0)[iens_dm], -1.0))  , "../data/tau_decay/"+Tag_reco_type+"/strange/corr/AK_dm_OS_"+ls_H_data_OS_AKAK.Tag[iens]+".dat");


       V0V0_tm_block_1_dm_distr = Corr_block_1.corr_t( summ_master( ls_ph_data_tm_V0V0.col(0)[iens_dm],  Multiply_Vvector_by_scalar(ls_uni_data_tm_V0V0.col(0)[iens_dm], -1.0)), "");
       VKVK_tm_block_1_dm_distr = Corr_block_1.corr_t( summ_master( ls_ph_data_tm_VKVK.col(0)[iens_dm],  Multiply_Vvector_by_scalar(ls_uni_data_tm_VKVK.col(0)[iens_dm], -1.0)), "");
       A0A0_tm_block_1_dm_distr = Corr_block_1.corr_t( summ_master( ls_ph_data_tm_A0A0.col(0)[iens_dm],  Multiply_Vvector_by_scalar(ls_uni_data_tm_A0A0.col(0)[iens_dm], -1.0)), "");
       AKAK_tm_block_1_dm_distr = Corr_block_1.corr_t( summ_master( ls_ph_data_tm_AKAK.col(0)[iens_dm],  Multiply_Vvector_by_scalar(ls_uni_data_tm_AKAK.col(0)[iens_dm], -1.0)), "");


       V0V0_OS_block_1_dm_distr = Corr_block_1.corr_t( summ_master( ls_ph_data_OS_V0V0.col(0)[iens_dm],  Multiply_Vvector_by_scalar(ls_uni_data_OS_V0V0.col(0)[iens_dm], -1.0)), "");
       VKVK_OS_block_1_dm_distr = Corr_block_1.corr_t( summ_master( ls_ph_data_OS_VKVK.col(0)[iens_dm],  Multiply_Vvector_by_scalar(ls_uni_data_OS_VKVK.col(0)[iens_dm], -1.0)), "");
       A0A0_OS_block_1_dm_distr = Corr_block_1.corr_t( summ_master( ls_ph_data_OS_A0A0.col(0)[iens_dm],  Multiply_Vvector_by_scalar(ls_uni_data_OS_A0A0.col(0)[iens_dm], -1.0)), "");
       AKAK_OS_block_1_dm_distr = Corr_block_1.corr_t( summ_master( ls_ph_data_OS_AKAK.col(0)[iens_dm],  Multiply_Vvector_by_scalar(ls_uni_data_OS_AKAK.col(0)[iens_dm], -1.0)), "");


    
       
       V0_tm_distr = V0_tm_distr + V0V0_tm_dm_distr;
       Vk_tm_distr = Vk_tm_distr + VKVK_tm_dm_distr;
       A0_tm_distr = A0_tm_distr + A0A0_tm_dm_distr;
       Ak_tm_distr = Ak_tm_distr + AKAK_tm_dm_distr;

       V0_OS_distr = V0_OS_distr + V0V0_OS_dm_distr;
       Vk_OS_distr = Vk_OS_distr + VKVK_OS_dm_distr;
       A0_OS_distr = A0_OS_distr + A0A0_OS_dm_distr;
       Ak_OS_distr = Ak_OS_distr + AKAK_OS_dm_distr;


       V0_tm_block_1_distr = V0_tm_block_1_distr + V0V0_tm_block_1_dm_distr;
       Vk_tm_block_1_distr = Vk_tm_block_1_distr + VKVK_tm_block_1_dm_distr;
       A0_tm_block_1_distr = A0_tm_block_1_distr + A0A0_tm_block_1_dm_distr;
       Ak_tm_block_1_distr = Ak_tm_block_1_distr + AKAK_tm_block_1_dm_distr;

       V0_OS_block_1_distr = V0_OS_block_1_distr + V0V0_OS_block_1_dm_distr;
       Vk_OS_block_1_distr = Vk_OS_block_1_distr + VKVK_OS_block_1_dm_distr;
       A0_OS_block_1_distr = A0_OS_block_1_distr + A0A0_OS_block_1_dm_distr;
       Ak_OS_block_1_distr = Ak_OS_block_1_distr + AKAK_OS_block_1_dm_distr;


       V0_H_tm_distr = V0_H_tm_distr + V0V0_tm_dm_distr;
       Vk_H_tm_distr = Vk_H_tm_distr + VKVK_tm_dm_distr;
       A0_H_tm_distr = A0_H_tm_distr + A0A0_tm_dm_distr;
       Ak_H_tm_distr = Ak_H_tm_distr + AKAK_tm_dm_distr;

       V0_H_OS_distr = V0_H_OS_distr + V0V0_OS_dm_distr;
       Vk_H_OS_distr = Vk_H_OS_distr + VKVK_OS_dm_distr;
       A0_H_OS_distr = A0_H_OS_distr + A0A0_OS_dm_distr;
       Ak_H_OS_distr = Ak_H_OS_distr + AKAK_OS_dm_distr;


       V0_H_tm_block_1_distr = V0_H_tm_block_1_distr + V0V0_tm_block_1_dm_distr;
       Vk_H_tm_block_1_distr = Vk_H_tm_block_1_distr + VKVK_tm_block_1_dm_distr;
       A0_H_tm_block_1_distr = A0_H_tm_block_1_distr + A0A0_tm_block_1_dm_distr;
       Ak_H_tm_block_1_distr = Ak_H_tm_block_1_distr + AKAK_tm_block_1_dm_distr;

       V0_H_OS_block_1_distr = V0_H_OS_block_1_distr + V0V0_OS_block_1_dm_distr;
       Vk_H_OS_block_1_distr = Vk_H_OS_block_1_distr + VKVK_OS_block_1_dm_distr;
       A0_H_OS_block_1_distr = A0_H_OS_block_1_distr + A0A0_OS_block_1_dm_distr;
       Ak_H_OS_block_1_distr = Ak_H_OS_block_1_distr + AKAK_OS_block_1_dm_distr;
       

     }
    
     //compute T and L correlators
     T_tm_distr= Za.ave()*Za.ave()*( Vk_tm_distr + (Zv*Zv/(Za*Za))*Ak_tm_distr);
     L_tm_distr= Zv.ave()*Zv.ave()*( A0_tm_distr + (Za*Za/(Zv*Zv))*V0_tm_distr);
     T_OS_distr= Zv.ave()*Zv.ave()*( Vk_OS_distr + (Za*Za/(Zv*Zv))*Ak_OS_distr);
     L_OS_distr= Za.ave()*Za.ave()*( A0_OS_distr + (Zv*Zv/(Za*Za))*V0_OS_distr);
     T_tm_block_1_distr = Za.ave()*Za.ave()*( Vk_tm_block_1_distr + (Zv*Zv/(Za*Za)).ave()*Ak_tm_block_1_distr);
     L_tm_block_1_distr = Zv.ave()*Zv.ave()*( A0_tm_block_1_distr + (Za*Za/(Zv*Zv)).ave()*V0_tm_block_1_distr);
     T_OS_block_1_distr = Zv.ave()*Zv.ave()*( Vk_OS_block_1_distr + (Za*Za/(Zv*Zv)).ave()*Ak_OS_block_1_distr);
     L_OS_block_1_distr = Za.ave()*Za.ave()*( A0_OS_block_1_distr + (Zv*Zv/(Za*Za)).ave()*V0_OS_block_1_distr);

     
     T_H_tm_distr= Za.ave()*Za.ave()*( Vk_H_tm_distr + (Zv*Zv/(Za*Za)).ave()*Ak_H_tm_distr);
     L_H_tm_distr= Zv.ave()*Zv.ave()*( A0_H_tm_distr + (Za*Za/(Zv*Zv)).ave()*V0_H_tm_distr);
     T_H_OS_distr= Zv.ave()*Zv.ave()*( Vk_H_OS_distr + (Za*Za/(Zv*Zv)).ave()*Ak_H_OS_distr);
     L_H_OS_distr= Za.ave()*Za.ave()*( A0_H_OS_distr + (Zv*Zv/(Za*Za)).ave()*V0_H_OS_distr);
     T_H_tm_block_1_distr = Za.ave()*Za.ave()*( Vk_H_tm_block_1_distr + (Zv*Zv/(Za*Za)).ave()*Ak_H_tm_block_1_distr);
     L_H_tm_block_1_distr = Zv.ave()*Zv.ave()*( A0_H_tm_block_1_distr + (Za*Za/(Zv*Zv)).ave()*V0_H_tm_block_1_distr);
     T_H_OS_block_1_distr = Zv.ave()*Zv.ave()*( Vk_H_OS_block_1_distr + (Za*Za/(Zv*Zv)).ave()*Ak_H_OS_block_1_distr);
     L_H_OS_block_1_distr = Za.ave()*Za.ave()*( A0_H_OS_block_1_distr + (Zv*Zv/(Za*Za)).ave()*V0_H_OS_block_1_distr);
     
     
    //set plateaux for J^P = 1^- and J^P = 1^+

    int Tmin_1_minus_tm=0; int Tmax_1_minus_tm=0;
    int Tmin_1_plus_tm=0;  int Tmax_1_plus_tm=0;
    int Tmin_1_minus_OS=0; int Tmax_1_minus_OS=0;
    int Tmin_1_plus_OS=0;  int Tmax_1_plus_OS=0;
    
    if(ll_data_OS_P5P5.Tag[iens] == "cB211b.072.64" || ll_data_OS_P5P5.Tag[iens] == "cB211b.072.96" ) {
      
      Tmin_1_minus_tm=21;   Tmax_1_minus_tm=26;
      Tmin_1_plus_tm=18;   Tmax_1_plus_tm=21;
      
      Tmin_1_minus_OS=21;   Tmax_1_minus_OS=26;
      Tmin_1_plus_OS=18;   Tmax_1_plus_OS=21;
    }
    else if(ll_data_OS_P5P5.Tag[iens] == "cC211a.06.80") {
       Tmin_1_minus_tm=24;   Tmax_1_minus_tm=32;
       Tmin_1_plus_tm=16;   Tmax_1_plus_tm=20;
       
       Tmin_1_minus_OS=24;   Tmax_1_minus_OS=32;
       Tmin_1_plus_OS=16;   Tmax_1_plus_OS=20;
       
     }
    else if(ll_data_OS_P5P5.Tag[iens] == "cC211a.06.112") {
       
       Tmin_1_minus_tm=24;   Tmax_1_minus_tm=29;
       Tmin_1_plus_tm=16;   Tmax_1_plus_tm=20;
       
       Tmin_1_minus_OS=24;   Tmax_1_minus_OS=29;
       Tmin_1_plus_OS=16;   Tmax_1_plus_OS=20;
       
    }
    else if(ll_data_OS_P5P5.Tag[iens] == "cD211a.054.96") {
       
      Tmin_1_minus_tm=26;   Tmax_1_minus_tm=33;
       Tmin_1_plus_tm=24;   Tmax_1_plus_tm=30;
       
       Tmin_1_minus_OS=26;   Tmax_1_minus_OS=33;
       Tmin_1_plus_OS=24;   Tmax_1_plus_OS=30;
	
    }
    else if(ll_data_OS_P5P5.Tag[iens] == "cE211a.044.112") {
       
      Tmin_1_minus_tm=30;   Tmax_1_minus_tm=38;
       Tmin_1_plus_tm=28;   Tmax_1_plus_tm=36;
       
       Tmin_1_minus_OS=30;   Tmax_1_minus_OS=38;
       Tmin_1_plus_OS=26;   Tmax_1_plus_OS=36;
	
    }
    
    else crash("Cannot find ensemble: "+ll_data_OS_P5P5.Tag[iens]);
    
    Corr.Tmin=Tmin_1_minus_tm; Corr.Tmax=Tmax_1_minus_tm;
    distr_t m_Kstar_tm= Corr.Fit_distr(M_Kstar_tm);
    distr_t A_Kstar_tm= Corr.Fit_distr( Corr.residue_t(Vk_tm_distr, "")/(2.0*m_Kstar_tm));
    distr_t_list Vii_tm_gs= (A_Kstar_tm).ave()*(EXPT_D(-1.0*m_Kstar_tm.ave()*Get_id_jack_distr(Njacks), Corr.Nt) + exp(-m_Kstar_tm.ave()*Corr.Nt)*EXPT_D( m_Kstar_tm.ave()*Get_id_jack_distr(Njacks), Corr.Nt) ); 
    Vii_tm_gs.distr_list[0] = 0.0*Vii_tm_gs.distr_list[0];
    
    Corr.Tmin=Tmin_1_plus_tm; Corr.Tmax=Tmax_1_plus_tm;
    distr_t m_K1_tm= Corr.Fit_distr(M_K1_tm);
    distr_t A_K1_tm= Corr.Fit_distr( Corr.residue_t(Ak_tm_distr, "")/(2.0*m_K1_tm));
    distr_t_list Aii_tm_gs= (A_K1_tm).ave()*(EXPT_D(-1.0*m_K1_tm.ave()*Get_id_jack_distr(Njacks), Corr.Nt) + exp(-m_K1_tm.ave()*Corr.Nt)*EXPT_D( m_K1_tm.ave()*Get_id_jack_distr(Njacks), Corr.Nt) );
    Aii_tm_gs.distr_list[0] = 0.0*Aii_tm_gs.distr_list[0];

    Corr.Tmin=Tmin_1_minus_OS; Corr.Tmax=Tmax_1_minus_OS;
    distr_t m_Kstar_OS= Corr.Fit_distr(M_Kstar_OS);
    distr_t A_Kstar_OS= Corr.Fit_distr( Corr.residue_t(Vk_OS_distr, "")/(2.0*m_Kstar_OS));
    distr_t_list Vii_OS_gs= (A_Kstar_OS).ave()*( EXPT_D(-1.0*m_Kstar_OS.ave()*Get_id_jack_distr(Njacks), Corr.Nt) + exp(-m_Kstar_OS.ave()*Corr.Nt)*EXPT_D( m_Kstar_OS.ave()*Get_id_jack_distr(Njacks), Corr.Nt) );  
    Vii_OS_gs.distr_list[0] = 0.0*Vii_OS_gs.distr_list[0];

    Corr.Tmin=Tmin_1_plus_OS; Corr.Tmax=Tmax_1_plus_OS;
    distr_t m_K1_OS= Corr.Fit_distr(M_K1_OS);
    distr_t A_K1_OS= Corr.Fit_distr( Corr.residue_t(Ak_OS_distr, "")/(2.0*m_K1_OS));
    distr_t_list Aii_OS_gs= (A_K1_OS).ave()*(EXPT_D(-1.0*m_K1_OS.ave()*Get_id_jack_distr(Njacks), Corr.Nt) + exp(-m_K1_OS.ave()*Corr.Nt)*EXPT_D( m_K1_OS.ave()*Get_id_jack_distr(Njacks), Corr.Nt) );
    Aii_OS_gs.distr_list[0] = 0.0*Aii_OS_gs.distr_list[0];

    
    


    int Tmin_P5=0;
    int Tmax_P5=0;
    double amL1=0;
    double amSS1=0;
    distr_t aMP(UseJack);
    if( ls_data_tm_VKVK.Tag[iens] =="cB211b.072.96")     { Tmin_P5=38; Tmax_P5=59; amL1= 0.0006675; amSS1=0.018254 ; aMP=aMp_B;}
    else if(ls_data_tm_VKVK.Tag[iens] =="cB211b.072.64") { Tmin_P5=45; Tmax_P5=59; amL1= 0.0006675; amSS1=0.018254 ; aMP=aMp_B;}
    else if(ls_data_tm_VKVK.Tag[iens].substr(1,1)=="C")  { Tmin_P5=36; Tmax_P5=63; amL1= 0.000585; amSS1=0.016067; aMP=aMp_C;}
    else if(ls_data_tm_VKVK.Tag[iens].substr(1,1)=="D")  { Tmin_P5=37; Tmax_P5=69; amL1= 0.0004964; amSS1=0.013557; aMP=aMp_D;}
    else if(ls_data_tm_VKVK.Tag[iens].substr(1,1)=="E")  { Tmin_P5=60; Tmax_P5=85; amL1= 0.000431; amSS1=0.011759; aMP=aMp_E;}
    else crash("Cannot recognize the ensemble: "+ls_data_tm_VKVK.Tag[iens]+" in assigning Tmin_P5,Tmax_P5 for ensemble: ");

    Corr.Tmin = Tmin_P5; Corr.Tmax= Tmax_P5;

    distr_t MK1= Corr.Fit_distr( M_K_tm )/a_distr;
    distr_t FK1= Corr.Fit_distr( Corr.decay_constant_t( pow(ams1+aml, 2)*P5_tm_distr, "../data/tau_decay/"+Tag_reco_type+"/strange/FK/FK1_"+ls_data_tm_VKVK.Tag[iens]))/a_distr;

    distr_t MK2= Corr.Fit_distr( M_K_H_tm )/a_distr;
    distr_t FK2= Corr.Fit_distr( Corr.decay_constant_t( pow(ams2+aml, 2)*P5_H_tm_distr,  "../data/tau_decay/"+Tag_reco_type+"/strange/FK/FK2_"+ls_data_tm_VKVK.Tag[iens]))/a_distr;

    distr_t FK_L1 = Corr.Fit_distr( Corr.decay_constant_t( pow(amSS1+amL1,2)*P5P5_tm_L1,  "../data/tau_decay/"+Tag_reco_type+"/strange/FK/FKL1_"+ls_data_tm_VKVK.Tag[iens]))/a_distr;
    distr_t FK_L2 = Corr.Fit_distr( Corr.decay_constant_t( pow(amSS1+aml,2)*P5P5_tm_L2,  "../data/tau_decay/"+Tag_reco_type+"/strange/FK/FKL2_"+ls_data_tm_VKVK.Tag[iens]))/a_distr;

    //determine physical strange quark mass

    vector<distr_t> MMK2({MK1*MK1, MK2*MK2});
    vector<distr_t> FFK({ FK1, FK2});
    vector<distr_t> MMS({ Get_id_jack_distr(Njacks)*ams1, Get_id_jack_distr(Njacks)*ams2});

    distr_t Mk_iso_corr= SQRT_D(  Mk_iso*Mk_iso + 0.5*( POW_D(aMP/a_distr,2)    - pow(0.135,2)));

    distr_t ams_phys = Obs_extrapolation_meson_mass(MMS, MMK2, Mk_iso*Mk_iso ,  "../data/tau_decay/"+Tag_reco_type+"/strange"  , "ams_extrapolation_"+ls_data_tm_VKVK.Tag[iens]+".dat",  UseJack, "SPLINE" );
    
    distr_t FK =  Obs_extrapolation_meson_mass(FFK, MMK2, Mk_iso_corr*Mk_iso_corr ,  "../data/tau_decay/"+Tag_reco_type+"/strange"  , "FK_extrapolation_"+ls_data_tm_VKVK.Tag[iens]+".dat",  UseJack, "SPLINE" ) + FK_L1 - FK_L2;

    FK_list.distr_list.push_back(FK);
    
       
    if(!Skip_spectral_density_analysis_strange) {

      double resc_GeV = C_V*GAMMA_FACT/(pow(a_distr.ave(),3));
      distr_t resc_GeV_distr= resc_GeV*Get_id_distr(Njacks,UseJack);
       
     
    //print covariance matrix

    //########### LIGHTER MASS #############

    Vfloat cov_A0_tm, cov_V0_tm,  cov_Ak_tm, cov_Vk_tm, cov_A0_OS, cov_V0_OS,  cov_Ak_OS, cov_Vk_OS, TT, RR;
    Vfloat corr_m_A0_tm, corr_m_V0_tm,  corr_m_Ak_tm, corr_m_Vk_tm, corr_m_A0_OS, corr_m_V0_OS, corr_m_Ak_OS, corr_m_Vk_OS;

    Vfloat cov_T_tm, cov_L_tm, cov_T_OS, cov_L_OS;
    Vfloat corr_m_T_tm, corr_m_L_tm, corr_m_T_OS, corr_m_L_OS;
    
    for(int tt=0;tt<Corr.Nt;tt++)
      for(int rr=0;rr<Corr.Nt;rr++) {
	TT.push_back(tt);
	RR.push_back(rr);

	double err_resc_A0_tm= A0_tm_distr.err(tt)*A0_tm_distr.err(rr)/(A0_tm_block_1_distr.err(tt)*A0_tm_block_1_distr.err(rr));
	double err_resc_V0_tm= V0_tm_distr.err(tt)*V0_tm_distr.err(rr)/(V0_tm_block_1_distr.err(tt)*V0_tm_block_1_distr.err(rr));
	double err_resc_Ak_tm= Ak_tm_distr.err(tt)*Ak_tm_distr.err(rr)/(Ak_tm_block_1_distr.err(tt)*Ak_tm_block_1_distr.err(rr));
	double err_resc_Vk_tm= Vk_tm_distr.err(tt)*Vk_tm_distr.err(rr)/(Vk_tm_block_1_distr.err(tt)*Vk_tm_block_1_distr.err(rr));


	double err_resc_T_tm= T_tm_distr.err(tt)*T_tm_distr.err(rr)/(T_tm_block_1_distr.err(tt)*T_tm_block_1_distr.err(rr));
	double err_resc_L_tm= L_tm_distr.err(tt)*L_tm_distr.err(rr)/(L_tm_block_1_distr.err(tt)*L_tm_block_1_distr.err(rr));

	

	double err_resc_A0_OS= A0_OS_distr.err(tt)*A0_OS_distr.err(rr)/(A0_OS_block_1_distr.err(tt)*A0_OS_block_1_distr.err(rr));
	double err_resc_V0_OS= V0_OS_distr.err(tt)*V0_OS_distr.err(rr)/(V0_OS_block_1_distr.err(tt)*V0_OS_block_1_distr.err(rr));
	double err_resc_Ak_OS= Ak_OS_distr.err(tt)*Ak_OS_distr.err(rr)/(Ak_OS_block_1_distr.err(tt)*Ak_OS_block_1_distr.err(rr));
	double err_resc_Vk_OS= Vk_OS_distr.err(tt)*Vk_OS_distr.err(rr)/(Vk_OS_block_1_distr.err(tt)*Vk_OS_block_1_distr.err(rr));

	double err_resc_T_OS= T_OS_distr.err(tt)*T_OS_distr.err(rr)/(T_OS_block_1_distr.err(tt)*T_OS_block_1_distr.err(rr));
	double err_resc_L_OS= L_OS_distr.err(tt)*L_OS_distr.err(rr)/(L_OS_block_1_distr.err(tt)*L_OS_block_1_distr.err(rr));
	
	
	cov_A0_tm.push_back( (A0_tm_block_1_distr.distr_list[tt]%A0_tm_block_1_distr.distr_list[rr])*err_resc_A0_tm);
	cov_V0_tm.push_back( (V0_tm_block_1_distr.distr_list[tt]%V0_tm_block_1_distr.distr_list[rr])*err_resc_V0_tm);
	cov_Ak_tm.push_back( (Ak_tm_block_1_distr.distr_list[tt]%Ak_tm_block_1_distr.distr_list[rr])*err_resc_Ak_tm);
	cov_Vk_tm.push_back( (Vk_tm_block_1_distr.distr_list[tt]%Vk_tm_block_1_distr.distr_list[rr])*err_resc_Vk_tm);
	cov_A0_OS.push_back( (A0_OS_block_1_distr.distr_list[tt]%A0_OS_block_1_distr.distr_list[rr])*err_resc_A0_OS);
	cov_V0_OS.push_back( (V0_OS_block_1_distr.distr_list[tt]%V0_OS_block_1_distr.distr_list[rr])*err_resc_V0_OS);
	cov_Ak_OS.push_back( (Ak_OS_block_1_distr.distr_list[tt]%Ak_OS_block_1_distr.distr_list[rr])*err_resc_Ak_OS);
	cov_Vk_OS.push_back( (Vk_OS_block_1_distr.distr_list[tt]%Vk_OS_block_1_distr.distr_list[rr])*err_resc_Vk_OS);


	cov_T_tm.push_back( (T_tm_block_1_distr.distr_list[tt]%T_tm_block_1_distr.distr_list[rr])*err_resc_T_tm);
	cov_L_tm.push_back( (L_tm_block_1_distr.distr_list[tt]%L_tm_block_1_distr.distr_list[rr])*err_resc_L_tm);
	cov_T_OS.push_back( (T_OS_block_1_distr.distr_list[tt]%T_OS_block_1_distr.distr_list[rr])*err_resc_T_OS);
	cov_L_OS.push_back( (L_OS_block_1_distr.distr_list[tt]%L_OS_block_1_distr.distr_list[rr])*err_resc_L_OS);


	corr_m_A0_tm.push_back( (A0_tm_block_1_distr.distr_list[tt]%A0_tm_block_1_distr.distr_list[rr])/(A0_tm_block_1_distr.err(tt)*A0_tm_block_1_distr.err(rr)));
	corr_m_V0_tm.push_back( (V0_tm_block_1_distr.distr_list[tt]%V0_tm_block_1_distr.distr_list[rr])/(V0_tm_block_1_distr.err(tt)*V0_tm_block_1_distr.err(rr)));
	corr_m_Ak_tm.push_back( (Ak_tm_block_1_distr.distr_list[tt]%Ak_tm_block_1_distr.distr_list[rr])/( Ak_tm_block_1_distr.err(tt)*Ak_tm_block_1_distr.err(rr)));
	corr_m_Vk_tm.push_back( (Vk_tm_block_1_distr.distr_list[tt]%Vk_tm_block_1_distr.distr_list[rr])/( Vk_tm_block_1_distr.err(tt)*Vk_tm_block_1_distr.err(rr)));
	corr_m_A0_OS.push_back( (A0_OS_block_1_distr.distr_list[tt]%A0_OS_block_1_distr.distr_list[rr])/( A0_OS_block_1_distr.err(tt)*A0_OS_block_1_distr.err(rr)));
	corr_m_V0_OS.push_back( (V0_OS_block_1_distr.distr_list[tt]%V0_OS_block_1_distr.distr_list[rr])/( V0_OS_block_1_distr.err(tt)*V0_OS_block_1_distr.err(rr)));
	corr_m_Ak_OS.push_back( (Ak_OS_block_1_distr.distr_list[tt]%Ak_OS_block_1_distr.distr_list[rr])/( Ak_OS_block_1_distr.err(tt)*Ak_OS_block_1_distr.err(rr)));
	corr_m_Vk_OS.push_back( (Vk_OS_block_1_distr.distr_list[tt]%Vk_OS_block_1_distr.distr_list[rr])/( Vk_OS_block_1_distr.err(tt)*Vk_OS_block_1_distr.err(rr)));


	corr_m_T_tm.push_back( (T_tm_block_1_distr.distr_list[tt]%T_tm_block_1_distr.distr_list[rr])/(T_tm_block_1_distr.err(tt)*T_tm_block_1_distr.err(rr)));
	corr_m_L_tm.push_back( (L_tm_block_1_distr.distr_list[tt]%L_tm_block_1_distr.distr_list[rr])/(L_tm_block_1_distr.err(tt)*L_tm_block_1_distr.err(rr)));
	corr_m_T_OS.push_back( (T_OS_block_1_distr.distr_list[tt]%T_OS_block_1_distr.distr_list[rr])/(T_OS_block_1_distr.err(tt)*T_OS_block_1_distr.err(rr)));
	corr_m_L_OS.push_back( (L_OS_block_1_distr.distr_list[tt]%L_OS_block_1_distr.distr_list[rr])/(L_OS_block_1_distr.err(tt)*L_OS_block_1_distr.err(rr)));
	
      }

    Print_To_File({}, {TT,RR,cov_A0_tm, corr_m_A0_tm}, "../data/tau_decay/"+Tag_reco_type+"/strange/covariance/A0_tm_"+ls_data_tm_VKVK.Tag[iens]+".dat", "", "");
    Print_To_File({}, {TT,RR,cov_V0_tm, corr_m_V0_tm}, "../data/tau_decay/"+Tag_reco_type+"/strange/covariance/V0_tm_"+ls_data_tm_VKVK.Tag[iens]+".dat", "", "");
    Print_To_File({}, {TT,RR,cov_Ak_tm, corr_m_Ak_tm}, "../data/tau_decay/"+Tag_reco_type+"/strange/covariance/Aii_tm_"+ls_data_tm_VKVK.Tag[iens]+".dat", "", "");
    Print_To_File({}, {TT,RR,cov_Vk_tm, corr_m_Vk_tm}, "../data/tau_decay/"+Tag_reco_type+"/strange/covariance/Vii_tm_"+ls_data_tm_VKVK.Tag[iens]+".dat", "", "");
    Print_To_File({}, {TT,RR,cov_A0_OS, corr_m_A0_OS}, "../data/tau_decay/"+Tag_reco_type+"/strange/covariance/A0_OS_"+ls_data_tm_VKVK.Tag[iens]+".dat", "", "");
    Print_To_File({}, {TT,RR,cov_V0_OS, corr_m_V0_OS}, "../data/tau_decay/"+Tag_reco_type+"/strange/covariance/V0_OS_"+ls_data_tm_VKVK.Tag[iens]+".dat", "", "");
    Print_To_File({}, {TT,RR,cov_Ak_OS, corr_m_Ak_OS}, "../data/tau_decay/"+Tag_reco_type+"/strange/covariance/Aii_OS_"+ls_data_tm_VKVK.Tag[iens]+".dat", "", "");
    Print_To_File({}, {TT,RR,cov_Vk_OS, corr_m_Vk_OS}, "../data/tau_decay/"+Tag_reco_type+"/strange/covariance/Vii_OS_"+ls_data_tm_VKVK.Tag[iens]+".dat", "", "");

    Print_To_File({}, {TT,RR,cov_T_tm, corr_m_T_tm}, "../data/tau_decay/"+Tag_reco_type+"/strange/covariance/T_tm_"+ls_data_tm_VKVK.Tag[iens]+".dat", "", "");
    Print_To_File({}, {TT,RR,cov_L_tm, corr_m_L_tm}, "../data/tau_decay/"+Tag_reco_type+"/strange/covariance/L_tm_"+ls_data_tm_VKVK.Tag[iens]+".dat", "", "");

    Print_To_File({}, {TT,RR,cov_T_OS, corr_m_T_OS}, "../data/tau_decay/"+Tag_reco_type+"/strange/covariance/T_OS_"+ls_data_OS_VKVK.Tag[iens]+".dat", "", "");
    Print_To_File({}, {TT,RR,cov_L_OS, corr_m_L_OS}, "../data/tau_decay/"+Tag_reco_type+"/strange/covariance/L_OS_"+ls_data_OS_VKVK.Tag[iens]+".dat", "", "");


    
    //########### HEAVIER MASS #############

    Vfloat cov_A0_H_tm, cov_V0_H_tm,  cov_Ak_H_tm, cov_Vk_H_tm, cov_A0_H_OS, cov_V0_H_OS,  cov_Ak_H_OS, cov_Vk_H_OS;
    Vfloat corr_m_A0_H_tm, corr_m_V0_H_tm,  corr_m_Ak_H_tm, corr_m_Vk_H_tm, corr_m_A0_H_OS, corr_m_V0_H_OS, corr_m_Ak_H_OS, corr_m_Vk_H_OS;

    Vfloat cov_T_H_tm, cov_L_H_tm, cov_T_H_OS, cov_L_H_OS;
    Vfloat corr_m_T_H_tm, corr_m_L_H_tm, corr_m_T_H_OS, corr_m_L_H_OS;
    
    for(int tt=0;tt<Corr.Nt;tt++)
      for(int rr=0;rr<Corr.Nt; rr++) {


	double err_resc_A0_H_tm= A0_H_tm_distr.err(tt)*A0_H_tm_distr.err(rr)/(A0_H_tm_block_1_distr.err(tt)*A0_H_tm_block_1_distr.err(rr));
	double err_resc_V0_H_tm= V0_H_tm_distr.err(tt)*V0_H_tm_distr.err(rr)/(V0_H_tm_block_1_distr.err(tt)*V0_H_tm_block_1_distr.err(rr));
	double err_resc_Ak_H_tm= Ak_H_tm_distr.err(tt)*Ak_H_tm_distr.err(rr)/(Ak_H_tm_block_1_distr.err(tt)*Ak_H_tm_block_1_distr.err(rr));
	double err_resc_Vk_H_tm= Vk_H_tm_distr.err(tt)*Vk_H_tm_distr.err(rr)/(Vk_H_tm_block_1_distr.err(tt)*Vk_H_tm_block_1_distr.err(rr));

	double err_resc_T_H_tm= T_H_tm_distr.err(tt)*T_H_tm_distr.err(rr)/(T_H_tm_block_1_distr.err(tt)*T_H_tm_block_1_distr.err(rr));
	double err_resc_L_H_tm= L_H_tm_distr.err(tt)*L_H_tm_distr.err(rr)/(L_H_tm_block_1_distr.err(tt)*L_H_tm_block_1_distr.err(rr));

	

	double err_resc_A0_H_OS= A0_H_OS_distr.err(tt)*A0_H_OS_distr.err(rr)/(A0_H_OS_block_1_distr.err(tt)*A0_H_OS_block_1_distr.err(rr));
	double err_resc_V0_H_OS= V0_H_OS_distr.err(tt)*V0_H_OS_distr.err(rr)/(V0_H_OS_block_1_distr.err(tt)*V0_H_OS_block_1_distr.err(rr));
	double err_resc_Ak_H_OS= Ak_H_OS_distr.err(tt)*Ak_H_OS_distr.err(rr)/(Ak_H_OS_block_1_distr.err(tt)*Ak_H_OS_block_1_distr.err(rr));
	double err_resc_Vk_H_OS= Vk_H_OS_distr.err(tt)*Vk_H_OS_distr.err(rr)/(Vk_H_OS_block_1_distr.err(tt)*Vk_H_OS_block_1_distr.err(rr));
	
	double err_resc_T_H_OS= T_H_OS_distr.err(tt)*T_H_OS_distr.err(rr)/(T_H_OS_block_1_distr.err(tt)*T_H_OS_block_1_distr.err(rr));
	double err_resc_L_H_OS= L_H_OS_distr.err(tt)*L_H_OS_distr.err(rr)/(L_H_OS_block_1_distr.err(tt)*L_H_OS_block_1_distr.err(rr));

	cov_A0_H_tm.push_back( (A0_H_tm_block_1_distr.distr_list[tt]%A0_H_tm_block_1_distr.distr_list[rr])*err_resc_A0_H_tm);
	cov_V0_H_tm.push_back( (V0_H_tm_block_1_distr.distr_list[tt]%V0_H_tm_block_1_distr.distr_list[rr])*err_resc_V0_H_tm);
	cov_Ak_H_tm.push_back( (Ak_H_tm_block_1_distr.distr_list[tt]%Ak_H_tm_block_1_distr.distr_list[rr])*err_resc_Ak_H_tm);
	cov_Vk_H_tm.push_back( (Vk_H_tm_block_1_distr.distr_list[tt]%Vk_H_tm_block_1_distr.distr_list[rr])*err_resc_Vk_H_tm);
	cov_A0_H_OS.push_back( (A0_H_OS_block_1_distr.distr_list[tt]%A0_H_OS_block_1_distr.distr_list[rr])*err_resc_A0_H_OS);
	cov_V0_H_OS.push_back( (V0_H_OS_block_1_distr.distr_list[tt]%V0_H_OS_block_1_distr.distr_list[rr])*err_resc_V0_H_OS);
	cov_Ak_H_OS.push_back( (Ak_H_OS_block_1_distr.distr_list[tt]%Ak_H_OS_block_1_distr.distr_list[rr])*err_resc_Ak_H_OS);
	cov_Vk_H_OS.push_back( (Vk_H_OS_block_1_distr.distr_list[tt]%Vk_H_OS_block_1_distr.distr_list[rr])*err_resc_Vk_H_OS);

	cov_T_H_tm.push_back( (T_H_tm_block_1_distr.distr_list[tt]%T_H_tm_block_1_distr.distr_list[rr])*err_resc_T_H_tm);
	cov_L_H_tm.push_back( (L_H_tm_block_1_distr.distr_list[tt]%L_H_tm_block_1_distr.distr_list[rr])*err_resc_L_H_tm);
	cov_T_H_OS.push_back( (T_H_OS_block_1_distr.distr_list[tt]%T_H_OS_block_1_distr.distr_list[rr])*err_resc_T_H_OS);
	cov_L_H_OS.push_back( (L_H_OS_block_1_distr.distr_list[tt]%L_H_OS_block_1_distr.distr_list[rr])*err_resc_L_H_OS);
	

	corr_m_A0_H_tm.push_back( (A0_H_tm_block_1_distr.distr_list[tt]%A0_H_tm_block_1_distr.distr_list[rr])/(A0_H_tm_block_1_distr.err(tt)*A0_H_tm_block_1_distr.err(rr)));
	corr_m_V0_H_tm.push_back( (V0_H_tm_block_1_distr.distr_list[tt]%V0_H_tm_block_1_distr.distr_list[rr])/(V0_H_tm_block_1_distr.err(tt)*V0_H_tm_block_1_distr.err(rr)));
	corr_m_Ak_H_tm.push_back( (Ak_H_tm_block_1_distr.distr_list[tt]%Ak_H_tm_block_1_distr.distr_list[rr])/( Ak_H_tm_block_1_distr.err(tt)*Ak_H_tm_block_1_distr.err(rr)));
	corr_m_Vk_H_tm.push_back( (Vk_H_tm_block_1_distr.distr_list[tt]%Vk_H_tm_block_1_distr.distr_list[rr])/( Vk_H_tm_block_1_distr.err(tt)*Vk_H_tm_block_1_distr.err(rr)));
	corr_m_A0_H_OS.push_back( (A0_H_OS_block_1_distr.distr_list[tt]%A0_H_OS_block_1_distr.distr_list[rr])/( A0_H_OS_block_1_distr.err(tt)*A0_H_OS_block_1_distr.err(rr)));
	corr_m_V0_H_OS.push_back( (V0_H_OS_block_1_distr.distr_list[tt]%V0_H_OS_block_1_distr.distr_list[rr])/( V0_H_OS_block_1_distr.err(tt)*V0_H_OS_block_1_distr.err(rr)));
	corr_m_Ak_H_OS.push_back( (Ak_H_OS_block_1_distr.distr_list[tt]%Ak_H_OS_block_1_distr.distr_list[rr])/( Ak_H_OS_block_1_distr.err(tt)*Ak_H_OS_block_1_distr.err(rr)));
	corr_m_Vk_H_OS.push_back( (Vk_H_OS_block_1_distr.distr_list[tt]%Vk_H_OS_block_1_distr.distr_list[rr])/( Vk_H_OS_block_1_distr.err(tt)*Vk_H_OS_block_1_distr.err(rr)));

	corr_m_T_H_tm.push_back( (T_H_tm_block_1_distr.distr_list[tt]%T_H_tm_block_1_distr.distr_list[rr])/(T_H_tm_block_1_distr.err(tt)*T_H_tm_block_1_distr.err(rr)));
	corr_m_L_H_tm.push_back( (L_H_tm_block_1_distr.distr_list[tt]%L_H_tm_block_1_distr.distr_list[rr])/(L_H_tm_block_1_distr.err(tt)*L_H_tm_block_1_distr.err(rr)));
	corr_m_T_H_OS.push_back( (T_H_OS_block_1_distr.distr_list[tt]%T_H_OS_block_1_distr.distr_list[rr])/(T_H_OS_block_1_distr.err(tt)*T_H_OS_block_1_distr.err(rr)));
	corr_m_L_H_OS.push_back( (L_H_OS_block_1_distr.distr_list[tt]%L_H_OS_block_1_distr.distr_list[rr])/(L_H_OS_block_1_distr.err(tt)*L_H_OS_block_1_distr.err(rr)));
	
      }


    //set to zero correlation and covariance if < 0.1
   

    Print_To_File({}, {TT,RR,cov_A0_H_tm, corr_m_A0_H_tm}, "../data/tau_decay/"+Tag_reco_type+"/strange/covariance/A0_H_tm_"+ls_H_data_tm_VKVK.Tag[iens]+".dat", "", "");
    Print_To_File({}, {TT,RR,cov_V0_H_tm, corr_m_V0_H_tm}, "../data/tau_decay/"+Tag_reco_type+"/strange/covariance/V0_H_tm_"+ls_H_data_tm_VKVK.Tag[iens]+".dat", "", "");
    Print_To_File({}, {TT,RR,cov_Ak_H_tm, corr_m_Ak_H_tm}, "../data/tau_decay/"+Tag_reco_type+"/strange/covariance/Aii_H_tm_"+ls_H_data_tm_VKVK.Tag[iens]+".dat", "", "");
    Print_To_File({}, {TT,RR,cov_Vk_H_tm, corr_m_Vk_H_tm}, "../data/tau_decay/"+Tag_reco_type+"/strange/covariance/Vii_H_tm_"+ls_H_data_tm_VKVK.Tag[iens]+".dat", "", "");
    Print_To_File({}, {TT,RR,cov_A0_H_OS, corr_m_A0_H_OS}, "../data/tau_decay/"+Tag_reco_type+"/strange/covariance/A0_H_OS_"+ls_H_data_tm_VKVK.Tag[iens]+".dat", "", "");
    Print_To_File({}, {TT,RR,cov_V0_H_OS, corr_m_V0_H_OS}, "../data/tau_decay/"+Tag_reco_type+"/strange/covariance/V0_H_OS_"+ls_H_data_tm_VKVK.Tag[iens]+".dat", "", "");
    Print_To_File({}, {TT,RR,cov_Ak_H_OS, corr_m_Ak_H_OS}, "../data/tau_decay/"+Tag_reco_type+"/strange/covariance/Aii_H_OS_"+ls_H_data_tm_VKVK.Tag[iens]+".dat", "", "");
    Print_To_File({}, {TT,RR,cov_Vk_H_OS, corr_m_Vk_H_OS}, "../data/tau_decay/"+Tag_reco_type+"/strange/covariance/Vii_H_OS_"+ls_H_data_tm_VKVK.Tag[iens]+".dat", "", "");

    Print_To_File({}, {TT,RR,cov_T_H_tm, corr_m_T_H_tm}, "../data/tau_decay/"+Tag_reco_type+"/strange/covariance/T_H_tm_"+ls_H_data_tm_VKVK.Tag[iens]+".dat", "", "");
    Print_To_File({}, {TT,RR,cov_L_H_tm, corr_m_L_H_tm}, "../data/tau_decay/"+Tag_reco_type+"/strange/covariance/L_H_tm_"+ls_H_data_tm_VKVK.Tag[iens]+".dat", "", "");
    Print_To_File({}, {TT,RR,cov_T_H_OS, corr_m_T_H_OS}, "../data/tau_decay/"+Tag_reco_type+"/strange/covariance/T_H_OS_"+ls_H_data_OS_VKVK.Tag[iens]+".dat", "", "");
    Print_To_File({}, {TT,RR,cov_L_H_OS, corr_m_L_H_OS}, "../data/tau_decay/"+Tag_reco_type+"/strange/covariance/L_H_OS_"+ls_H_data_OS_VKVK.Tag[iens]+".dat", "", "");
       

    distr_t_list A0_tm, V0_tm,  Aii_tm, A0_OS, V0_OS, Aii_OS, Vii_tm, Vii_OS;
    distr_t_list A0_H_tm, V0_H_tm,  Aii_H_tm, A0_H_OS, V0_H_OS, Aii_H_OS, Vii_H_tm, Vii_H_OS;


    distr_t_list T_tm , T_OS, L_tm, L_OS;
    distr_t_list T_H_tm, T_H_OS, L_H_tm, L_H_OS;
    
    
    //######### DEFINE 0th and ii component of C^munu ###########
    //lighter
    //tm
    A0_tm = A0_tm_distr;
    V0_tm = V0_tm_distr;
    Aii_tm =Ak_tm_distr;
    Vii_tm = Vk_tm_distr;
    T_tm = T_tm_distr;
    L_tm = L_tm_distr;
    //OS
    A0_OS = A0_OS_distr;
    V0_OS = V0_OS_distr;
    Aii_OS = Ak_OS_distr;
    Vii_OS = Vk_OS_distr;
    T_OS= T_OS_distr;
    L_OS= L_OS_distr;
    //heavier
    //tm
    A0_H_tm = A0_H_tm_distr;
    V0_H_tm = V0_H_tm_distr;
    Aii_H_tm =Ak_H_tm_distr;
    Vii_H_tm = Vk_H_tm_distr;
    T_H_tm= T_H_tm_distr;
    L_H_tm= L_H_tm_distr;
    //OS
    A0_H_OS = A0_H_OS_distr;
    V0_H_OS = V0_H_OS_distr;
    Aii_H_OS = Ak_H_OS_distr;
    Vii_H_OS = Vk_H_OS_distr;
    T_H_OS= T_H_OS_distr;
    L_H_OS= L_H_OS_distr;
    //###########################################################


    //SET [tmin-tmax]
    //################   LIGHTER   ##############################

    bool Found_error_less_x_percent=false;
    double x=10;
    //tm
    int tmax_tm_A0=1;
    while(!Found_error_less_x_percent && tmax_tm_A0 < Corr.Nt/2 -1 ) {
   
      if( (A0_tm.distr_list[tmax_tm_A0]).err()/fabs( (A0_tm.distr_list[tmax_tm_A0]).ave()) <  0.01*x) tmax_tm_A0++;
      else Found_error_less_x_percent=true;
    }

    Found_error_less_x_percent=false;
    
    int tmax_tm_V0=1;
    while(!Found_error_less_x_percent && tmax_tm_V0 < Corr.Nt/2 -1 ) {
   
      if( (V0_tm.distr_list[tmax_tm_V0]).err()/fabs( (V0_tm.distr_list[tmax_tm_V0]).ave()) <  0.01*x) tmax_tm_V0++;
      else Found_error_less_x_percent=true;
    }

    Found_error_less_x_percent=false;
    

    int tmax_tm_1_Aii=1;
    while(!Found_error_less_x_percent && tmax_tm_1_Aii < Corr.Nt/2 -1 ) {
   
      if( (Aii_tm.distr_list[tmax_tm_1_Aii]).err()/fabs( (Aii_tm.distr_list[tmax_tm_1_Aii]).ave()) <  0.01*x) tmax_tm_1_Aii++;
      else Found_error_less_x_percent=true;
    }

    Found_error_less_x_percent=false;

    int tmax_tm_1_Vii=1;
    while(!Found_error_less_x_percent && tmax_tm_1_Vii < Corr.Nt/2 -1 ) {
   
      if( (Vii_tm.distr_list[tmax_tm_1_Vii]).err()/fabs( (Vii_tm.distr_list[tmax_tm_1_Vii]).ave()) <  0.01*x) tmax_tm_1_Vii++;
      else Found_error_less_x_percent=true;
    }

    Found_error_less_x_percent=false;

    int tmax_tm_T=1;
    while(!Found_error_less_x_percent && tmax_tm_T < Corr.Nt/2 -1 ) {
   
      if( (T_tm.distr_list[tmax_tm_T]).err()/fabs( (T_tm.distr_list[tmax_tm_T]).ave()) <  0.01*x) tmax_tm_T++;
      else Found_error_less_x_percent=true;
    }

    Found_error_less_x_percent=false;


    int tmax_tm_L=1;
    while(!Found_error_less_x_percent && tmax_tm_L < Corr.Nt/2 -1 ) {
   
      if( (L_tm.distr_list[tmax_tm_L]).err()/fabs( (L_tm.distr_list[tmax_tm_L]).ave()) <  0.01*x) tmax_tm_L++;
      else Found_error_less_x_percent=true;
    }

    Found_error_less_x_percent=false;

    //#########################################################################
    //OS

    int tmax_OS_A0=1;
      while(!Found_error_less_x_percent && tmax_OS_A0 < Corr.Nt/2 -1 ) {
   
      if( (A0_OS.distr_list[tmax_OS_A0]).err()/fabs( (A0_OS.distr_list[tmax_OS_A0]).ave()) <  0.01*x) tmax_OS_A0++;
      else Found_error_less_x_percent=true;
    }

    Found_error_less_x_percent=false;


    int tmax_OS_V0=1;
    while(!Found_error_less_x_percent && tmax_OS_V0 < Corr.Nt/2 -1 ) {
   
      if( (V0_OS.distr_list[tmax_OS_V0]).err()/fabs( (V0_OS.distr_list[tmax_OS_V0]).ave()) <  0.01*x) tmax_OS_V0++;
      else Found_error_less_x_percent=true;
    }

    Found_error_less_x_percent=false;

    int tmax_OS_1_Aii=1;
    while(!Found_error_less_x_percent && tmax_OS_1_Aii < Corr.Nt/2 -1) {
   
      if( (Aii_OS.distr_list[tmax_OS_1_Aii]).err()/fabs( (Aii_OS.distr_list[tmax_OS_1_Aii]).ave()) <  0.01*x) tmax_OS_1_Aii++;
      else Found_error_less_x_percent=true;
    }

    Found_error_less_x_percent=false;

    int tmax_OS_1_Vii=1;
    while(!Found_error_less_x_percent && tmax_OS_1_Vii < Corr.Nt/2 -1) {
   
      if( (Vii_OS.distr_list[tmax_OS_1_Vii]).err()/fabs( (Vii_OS.distr_list[tmax_OS_1_Vii]).ave()) <  0.01*x) tmax_OS_1_Vii++;
      else Found_error_less_x_percent=true;
    }

    Found_error_less_x_percent=false;
    
    int tmax_OS_T=1;
    while(!Found_error_less_x_percent && tmax_OS_T < Corr.Nt/2 -1 ) {
   
      if( (T_OS.distr_list[tmax_OS_T]).err()/fabs( (T_OS.distr_list[tmax_OS_T]).ave()) <  0.01*x) tmax_OS_T++;
      else Found_error_less_x_percent=true;
    }

    Found_error_less_x_percent=false;


    int tmax_OS_L=1;
    while(!Found_error_less_x_percent && tmax_OS_L < Corr.Nt/2 -1 ) {
   
      if( (L_OS.distr_list[tmax_OS_L]).err()/fabs( (L_OS.distr_list[tmax_OS_L]).ave()) <  0.01*x) tmax_OS_L++;
      else Found_error_less_x_percent=true;
    }

    Found_error_less_x_percent=false;


    

    //################   HEAVIER   ##############################

    Found_error_less_x_percent=false;
    //tm
    int tmax_H_tm_A0=1;
    while(!Found_error_less_x_percent && tmax_H_tm_A0 < Corr.Nt/2 -1 ) {
   
      if( (A0_H_tm.distr_list[tmax_H_tm_A0]).err()/fabs( (A0_H_tm.distr_list[tmax_H_tm_A0]).ave()) <  0.01*x) tmax_H_tm_A0++;
      else Found_error_less_x_percent=true;
    }

    Found_error_less_x_percent=false;
    
    int tmax_H_tm_V0=1;
    while(!Found_error_less_x_percent && tmax_H_tm_V0 < Corr.Nt/2 -1 ) {
   
      if( (V0_H_tm.distr_list[tmax_H_tm_V0]).err()/fabs( (V0_H_tm.distr_list[tmax_H_tm_V0]).ave()) <  0.01*x) tmax_H_tm_V0++;
      else Found_error_less_x_percent=true;
    }

    Found_error_less_x_percent=false;
    

    int tmax_H_tm_1_Aii=1;
    while(!Found_error_less_x_percent && tmax_H_tm_1_Aii < Corr.Nt/2 -1 ) {
   
      if( (Aii_H_tm.distr_list[tmax_H_tm_1_Aii]).err()/fabs( (Aii_H_tm.distr_list[tmax_H_tm_1_Aii]).ave()) <  0.01*x) tmax_H_tm_1_Aii++;
      else Found_error_less_x_percent=true;
    }

    Found_error_less_x_percent=false;

    int tmax_H_tm_1_Vii=1;
    while(!Found_error_less_x_percent && tmax_H_tm_1_Vii < Corr.Nt/2 -1 ) {
   
      if( (Vii_H_tm.distr_list[tmax_H_tm_1_Vii]).err()/fabs( (Vii_H_tm.distr_list[tmax_H_tm_1_Vii]).ave()) <  0.01*x) tmax_H_tm_1_Vii++;
      else Found_error_less_x_percent=true;
    }

    Found_error_less_x_percent=false;


    int tmax_H_tm_T=1;
    while(!Found_error_less_x_percent && tmax_H_tm_T < Corr.Nt/2 -1 ) {
   
      if( (T_H_tm.distr_list[tmax_H_tm_T]).err()/fabs( (T_H_tm.distr_list[tmax_H_tm_T]).ave()) <  0.01*x) tmax_H_tm_T++;
      else Found_error_less_x_percent=true;
    }

    Found_error_less_x_percent=false;


    int tmax_H_tm_L=1;
    while(!Found_error_less_x_percent && tmax_H_tm_L < Corr.Nt/2 -1 ) {
   
      if( (L_H_tm.distr_list[tmax_H_tm_L]).err()/fabs( (L_H_tm.distr_list[tmax_H_tm_L]).ave()) <  0.01*x) tmax_H_tm_L++;
      else Found_error_less_x_percent=true;
    }

    Found_error_less_x_percent=false;


    

    //OS

    int tmax_H_OS_A0=1;
    while(!Found_error_less_x_percent && tmax_H_OS_A0 < Corr.Nt/2 -1 ) {
   
      if( (A0_H_OS.distr_list[tmax_H_OS_A0]).err()/fabs( (A0_H_OS.distr_list[tmax_H_OS_A0]).ave()) <  0.01*x) tmax_H_OS_A0++;
      else Found_error_less_x_percent=true;
    }

    Found_error_less_x_percent=false;


    int tmax_H_OS_V0=1;
    while(!Found_error_less_x_percent && tmax_H_OS_V0 < Corr.Nt/2 -1 ) {
   
      if( (V0_H_OS.distr_list[tmax_H_OS_V0]).err()/fabs( (V0_H_OS.distr_list[tmax_H_OS_V0]).ave()) <  0.01*x) tmax_H_OS_V0++;
      else Found_error_less_x_percent=true;
    }

    Found_error_less_x_percent=false;

    int tmax_H_OS_1_Aii=1;
    while(!Found_error_less_x_percent && tmax_H_OS_1_Aii < Corr.Nt/2 -1) {
   
      if( (Aii_H_OS.distr_list[tmax_H_OS_1_Aii]).err()/fabs( (Aii_H_OS.distr_list[tmax_H_OS_1_Aii]).ave()) <  0.01*x) tmax_H_OS_1_Aii++;
      else Found_error_less_x_percent=true;
    }

    Found_error_less_x_percent=false;

    int tmax_H_OS_1_Vii=1;
    while(!Found_error_less_x_percent && tmax_H_OS_1_Vii < Corr.Nt/2 -1) {
   
      if( (Vii_H_OS.distr_list[tmax_H_OS_1_Vii]).err()/fabs( (Vii_H_OS.distr_list[tmax_H_OS_1_Vii]).ave()) <  0.01*x) tmax_H_OS_1_Vii++;
      else Found_error_less_x_percent=true;
    }

    Found_error_less_x_percent=false;
    
    int tmax_H_OS_T=1;
    while(!Found_error_less_x_percent && tmax_H_OS_T < Corr.Nt/2 -1 ) {
   
      if( (T_H_OS.distr_list[tmax_H_OS_T]).err()/fabs( (T_H_OS.distr_list[tmax_H_OS_T]).ave()) <  0.01*x) tmax_H_OS_T++;
      else Found_error_less_x_percent=true;
    }

    Found_error_less_x_percent=false;


    int tmax_H_OS_L=1;
    while(!Found_error_less_x_percent && tmax_H_OS_L < Corr.Nt/2 -1 ) {
   
      if( (L_H_OS.distr_list[tmax_H_OS_L]).err()/fabs( (L_H_OS.distr_list[tmax_H_OS_L]).ave()) <  0.01*x) tmax_H_OS_L++;
      else Found_error_less_x_percent=true;
    }

    Found_error_less_x_percent=false;
    
    


   
      
    const auto K0 = [&a_distr](const PrecFloat &E, const PrecFloat &m, const PrecFloat &s, const PrecFloat &E0, int ijack) -> PrecFloat {

      
		      
		      PrecFloat X;
		      PrecFloat X_ave = E/(m_tau*a_distr.ave());
		      PrecFloat XX= E/(m_tau*a_distr.ave());
		      if(X_ave < E0) return 0.0;
		      
		      //if(X_ave > 1e4) return 0.0;

		  

                      if(ijack == -1)  X=X_ave;
		      else 	X= E/(m_tau*a_distr.distr[ijack]);
		      
		      
		      PrecFloat sm_theta;

		      if(sm_func_mode==0) sm_theta= 1/(1+ exp(-(1-X)/s));
		      else if(sm_func_mode==1) sm_theta= 1/(1+ exp(-sinh((1-X)/s)));
		      else if(sm_func_mode==2) sm_theta= (1+erf((1-X)/s))/2;
		      else crash("sm_func_mode: "+to_string(sm_func_mode)+" not yet implemented");
						 
		      return (1/XX)*pow(( 1 -pow(X,2)),2)*sm_theta;
		   
		 };

    auto K0_dub = [&a_distr](double &E, double &s) -> double {

      
		      
      double X = E/(m_tau*a_distr.ave());
      double XX= E/(m_tau*a_distr.ave());
		     
          
      double sm_theta;
      
      if(sm_func_mode==0) sm_theta= 1/(1+ exp(-(1-X)/s));
      else if(sm_func_mode==1) sm_theta= 1/(1+ exp(-sinh((1-X)/s)));
      else if(sm_func_mode==2) sm_theta= (1+erf((1-X)/s))/2;
      else crash("sm_func_mode: "+to_string(sm_func_mode)+" not yet implemented");
      
      return (1/XX)*pow(( 1 -pow(X,2)),2)*sm_theta;
		   
    };

    const auto K1 = [&a_distr](const PrecFloat &E, const PrecFloat &m, const PrecFloat &s, const PrecFloat &E0, int ijack) -> PrecFloat {

		      PrecFloat X;
		      PrecFloat X_ave = E/(m_tau*a_distr.ave());
		      PrecFloat XX= E/(m_tau*a_distr.ave());
		      if( X_ave < E0) return 0.0;

		      //if(X_ave > 1e4) return 0.0;

		  

		      if(ijack==-1) {
			X=X_ave;
		      }
		      else X= E/(m_tau*a_distr.distr[ijack]);

		      
		      PrecFloat sm_theta;
		      if(sm_func_mode==0) sm_theta= 1/(1+ exp(-(1-X)/s));
		      else if(sm_func_mode==1) sm_theta= 1/(1+ exp(-sinh((1-X)/s)));
		      else if(sm_func_mode==2) sm_theta= (1+erf((1-X)/s))/2;
		      else crash("sm_func_mode: "+to_string(sm_func_mode)+" not yet implemented");
		      
		      return (1 + 2*pow(X,2))*(1/(XX))*pow(( 1 -pow(X,2)),2)*sm_theta;
		   
		    };

    auto K1_dub = [&a_distr](double E, double s) -> double {

      
		      
      double X = E/(m_tau*a_distr.ave());
      double XX= E/(m_tau*a_distr.ave());
      
      
      double sm_theta;

           
      if(sm_func_mode==0) sm_theta= 1/(1+ exp(-(1-X)/s));
      else if(sm_func_mode==1) sm_theta= 1/(1+ exp(-sinh((1-X)/s)));
      else if(sm_func_mode==2) sm_theta= (1+erf((1-X)/s))/2;
      else crash("sm_func_mode: "+to_string(sm_func_mode)+" not yet implemented");

      return (1 + 2*pow(X,2))*(1/(XX))*pow(( 1 -pow(X,2)),2)*sm_theta;
     		   
    };



    const auto K0_shifted = [&a_distr](const PrecFloat &E, const PrecFloat &m, const PrecFloat &s, const PrecFloat &E0, int ijack) -> PrecFloat {

      
		      
		      PrecFloat X;
		      PrecFloat X_ave = E/(m_tau*(a_distr.ave() + a_distr.err()));
		      if(X_ave < E0) return 0.0;
		      
		      //if(X_ave > 1e4) return 0.0;

		      PrecFloat FF= pow(a_distr.ave()/(a_distr.ave()+a_distr.err()),3);


                      if(ijack == -1)  X=X_ave;
		      else 	X= E/(m_tau*a_distr.distr[ijack]);
		      
		      
		      PrecFloat sm_theta;

		      if(sm_func_mode==0) sm_theta= 1/(1+ exp(-(1-X)/s));
		      else if(sm_func_mode==1) sm_theta= 1/(1+ exp(-sinh((1-X)/s)));
		      else if(sm_func_mode==2) sm_theta= (1+erf((1-X)/s))/2;
		      else crash("sm_func_mode: "+to_string(sm_func_mode)+" not yet implemented");
						 
		      return FF*(1/X)*pow(( 1 -pow(X,2)),2)*sm_theta;
		   
		 };

    const auto K1_shifted = [&a_distr](const PrecFloat &E, const PrecFloat &m, const PrecFloat &s, const PrecFloat &E0, int ijack) -> PrecFloat {

		      PrecFloat X;
		      PrecFloat X_ave = E/(m_tau*(a_distr.ave() + a_distr.err()));
		      if( X_ave < E0) return 0.0;

		      //if(X_ave > 1e4) return 0.0;

		      PrecFloat FF= pow(a_distr.ave()/(a_distr.ave()+a_distr.err()),3);


		      if(ijack==-1) {
			X=X_ave;
		      }
		      else X= E/(m_tau*a_distr.distr[ijack]);

		      
		      PrecFloat sm_theta;
		      if(sm_func_mode==0) sm_theta= 1/(1+ exp(-(1-X)/s));
		      else if(sm_func_mode==1) sm_theta= 1/(1+ exp(-sinh((1-X)/s)));
		      else if(sm_func_mode==2) sm_theta= (1+erf((1-X)/s))/2;
		      else crash("sm_func_mode: "+to_string(sm_func_mode)+" not yet implemented");
		      
		      return FF*(1 + 2*pow(X,2))*(1/(X))*pow(( 1 -pow(X,2)),2)*sm_theta;
		   
		    };
    


    if(Use_t_up_to_T_half_strange) {
      tmax_tm_1_Aii = tmax_tm_1_Vii= tmax_tm_A0 = tmax_OS_A0 = tmax_tm_V0 = tmax_OS_V0 = tmax_OS_1_Aii = tmax_OS_1_Vii= Corr.Nt/2 -1;
      tmax_H_tm_1_Aii = tmax_H_tm_1_Vii= tmax_H_tm_A0 = tmax_H_OS_A0 = tmax_H_tm_V0 = tmax_H_OS_V0 = tmax_H_OS_1_Aii = tmax_H_OS_1_Vii= Corr.Nt/2 -1;

      tmax_tm_T = tmax_tm_L = tmax_OS_T = tmax_OS_L = Corr.Nt/2 -1;
      tmax_H_tm_T = tmax_H_tm_L = tmax_H_OS_T = tmax_H_OS_L = Corr.Nt/2 -1;
    }


    /*
    if( ls_H_data_tm_VKVK.Tag[iens] == "cB211b.072.96") {

      tmax_tm_A0 = 64;
      tmax_OS_A0 = 64;
      tmax_H_tm_A0 = 64;
      tmax_H_OS_A0 = 64;


      tmax_tm_1_Vii=tmax_OS_1_Vii=tmax_H_tm_1_Vii=tmax_H_OS_1_Vii=30;

    }

    if( ls_H_data_tm_VKVK.Tag[iens] == "cE211a.044.112") {

      tmax_tm_A0= 111;
      tmax_OS_A0= 111;
      tmax_H_tm_A0=111;
      tmax_H_OS_A0=111;

      } */
    
  
    cout<<"sigma list : {"<<sigma_list_strange[0];
    for(int is=1;is<(signed)sigma_list_strange.size();is++) { cout<<","<<sigma_list_strange[is]<<flush;}
    cout<<"}"<<endl<<flush;
    cout<<"Looping over sigma"<<flush;

    //loop over sigma
    vector<tuple<int,double, double, double, double>> thread_times_tm(sigma_list_strange.size()), thread_times_OS(sigma_list_strange.size());
    vector<tuple<int,double, double, double, double>> thread_times_H_tm(sigma_list_strange.size()), thread_times_H_OS(sigma_list_strange.size());
    
    
    
    #pragma omp parallel for schedule(dynamic)
    for(int is=0; is < (signed)sigma_list_strange.size(); is++) {

      double s= sigma_list_strange[is];

      
      //#####    LIGHTER  ##########
      distr_t Br_sigma_A0_tm;
      distr_t Br_sigma_V0_tm;
      distr_t Br_sigma_Aii_tm;
      distr_t Br_sigma_Vii_tm;
      distr_t Br_sigma_A0_OS;
      distr_t Br_sigma_V0_OS;
      distr_t Br_sigma_Aii_OS;
      distr_t Br_sigma_Vii_OS;

      distr_t Br_sigma_T_tm;
      distr_t Br_sigma_L_tm;
      distr_t Br_sigma_T_OS;
      distr_t Br_sigma_L_OS;

      distr_t Br_s_sigma_A0_tm;
      distr_t Br_s_sigma_V0_tm;
      distr_t Br_s_sigma_Aii_tm;
      distr_t Br_s_sigma_Vii_tm;
      distr_t Br_s_sigma_A0_OS;
      distr_t Br_s_sigma_V0_OS;
      distr_t Br_s_sigma_Aii_OS;
      distr_t Br_s_sigma_Vii_OS;

      distr_t Br_s_sigma_T_tm;
      distr_t Br_s_sigma_L_tm;
      distr_t Br_s_sigma_T_OS;
      distr_t Br_s_sigma_L_OS;
   
      //int tmax= T/2 -4;
      double lA0_tm, lAii_tm, lVii_tm, lV0_tm;
      double lA0_OS, lAii_OS, lVii_OS, lV0_OS;
      double lT_tm, lL_tm, lT_OS, lL_OS;

      double slA0_tm, slAii_tm, slVii_tm, slV0_tm;
      double slA0_OS, slAii_OS, slVii_OS, slV0_OS;
      double slT_tm, slL_tm, slT_OS, slL_OS;
      
      double syst_A0_tm, syst_Aii_tm, syst_Vii_tm, syst_V0_tm;
      double syst_A0_OS, syst_Aii_OS, syst_Vii_OS, syst_V0_OS;
      double syst_T_tm, syst_L_tm, syst_T_OS, syst_L_OS;
      

      double syst_s_A0_tm, syst_s_Aii_tm, syst_s_Vii_tm, syst_s_V0_tm;
      double syst_s_A0_OS, syst_s_Aii_OS, syst_s_Vii_OS, syst_s_V0_OS;
      double syst_s_T_tm, syst_s_L_tm, syst_s_T_OS, syst_s_L_OS;
      

      //#####    HEAVIER  ##########
      distr_t Br_sigma_A0_H_tm;
      distr_t Br_sigma_V0_H_tm;
      distr_t Br_sigma_Aii_H_tm;
      distr_t Br_sigma_Vii_H_tm;
      distr_t Br_sigma_A0_H_OS;
      distr_t Br_sigma_V0_H_OS;
      distr_t Br_sigma_Aii_H_OS;
      distr_t Br_sigma_Vii_H_OS;

      distr_t Br_sigma_T_H_tm;
      distr_t Br_sigma_L_H_tm;
      distr_t Br_sigma_T_H_OS;
      distr_t Br_sigma_L_H_OS;

      distr_t Br_s_sigma_A0_H_tm;
      distr_t Br_s_sigma_V0_H_tm;
      distr_t Br_s_sigma_Aii_H_tm;
      distr_t Br_s_sigma_Vii_H_tm;
      distr_t Br_s_sigma_A0_H_OS;
      distr_t Br_s_sigma_V0_H_OS;
      distr_t Br_s_sigma_Aii_H_OS;
      distr_t Br_s_sigma_Vii_H_OS;

      distr_t Br_s_sigma_T_H_tm;
      distr_t Br_s_sigma_L_H_tm;
      distr_t Br_s_sigma_T_H_OS;
      distr_t Br_s_sigma_L_H_OS;
   
      //int tmax= T/2 -4;
      double lA0_H_tm, lAii_H_tm, lVii_H_tm, lV0_H_tm;
      double lA0_H_OS, lAii_H_OS, lVii_H_OS, lV0_H_OS;
      double lT_H_tm, lL_H_tm, lT_H_OS, lL_H_OS;

      double slA0_H_tm, slAii_H_tm, slVii_H_tm, slV0_H_tm;
      double slA0_H_OS, slAii_H_OS, slVii_H_OS, slV0_H_OS;
      double slT_H_tm, slL_H_tm, slT_H_OS, slL_H_OS;
      
      double syst_A0_H_tm, syst_Aii_H_tm, syst_Vii_H_tm, syst_V0_H_tm;
      double syst_A0_H_OS, syst_Aii_H_OS, syst_Vii_H_OS, syst_V0_H_OS;
      double syst_T_H_tm, syst_L_H_tm, syst_T_H_OS, syst_L_H_OS;

      double syst_s_A0_H_tm, syst_s_Aii_H_tm, syst_s_Vii_H_tm, syst_s_V0_H_tm;
      double syst_s_A0_H_OS, syst_s_Aii_H_OS, syst_s_Vii_H_OS, syst_s_V0_H_OS;
      double syst_s_T_H_tm, syst_s_L_H_tm, syst_s_T_H_OS, syst_s_L_H_OS;
      
      
      double mult=1e4;
      if( (beta > 2) && (s < 0.15) ) mult=1e5;


      //###########################################################//
      //###########################################################//
      //###########################################################//
      //###########################################################//
      //##############                              ###############//
      //##############                              ###############//
      //##############                              ###############//
      //##############         LIGHTER              ###############//
      //##############                              ###############//
      //##############                              ###############//
      //##############                              ###############//


      
      auto start = chrono::system_clock::now();
      if( (beta > 2) && Use_Customized_plateaux_strange) mult=  Customized_plateaux_tau_spectre_strange( beta, Emax,  "Aii", "tm" , s, ls_data_tm_VKVK.Tag[iens] );
      Br_sigma_Aii_tm = Get_Laplace_transfo(  0.0,  s, E0_A_sp*a_distr.ave(),  T, tmax_tm_1_Aii, prec, SM_TYPE_1,K1, Aii_tm, syst_Aii_tm, mult, lAii_tm, MODE, "tm", "Aii_strange_"+ls_data_tm_VKVK.Tag[iens], 1e-3,0, resc_GeV*Zv*Zv, 0.0, "tau_decay", cov_Ak_tm, fake_func,0, fake_func_d ,  Is_Emax_Finite, Emax, beta);
      Br_s_sigma_Aii_tm = Get_Laplace_transfo(  0.0,  s, E0_A_sp*a_distr.ave(),  T, tmax_tm_1_Aii, prec, SM_TYPE_1,K1_shifted, Aii_tm, syst_s_Aii_tm, mult, slAii_tm, MODE, "tm", "Aii_s_strange_"+ls_data_tm_VKVK.Tag[iens], 1e-3,0, resc_GeV*Zv*Zv, 0.0, "tau_decay", cov_Ak_tm, fake_func,0, fake_func_d ,  Is_Emax_Finite, Emax, beta);
      //distr_t preco_Aii_tm = Get_Laplace_transfo(  0.0,  s, E0_A_sp*a_distr.ave(),  T, tmax_tm_1_Aii, prec, SM_TYPE_1,K1, Aii_tm - Aii_tm_gs.ave(), syst_Aii_tm, mult, lAii_tm, MODE, "tm", "preco_Aii_strange_"+ls_data_tm_VKVK.Tag[iens], 1e-3,0, resc_GeV*Zv*Zv, A_K1_tm.ave()*K1_dub(m_K1_tm.ave(), s), "tau_decay", cov_Ak_tm, fake_func,0, fake_func_d ,  Is_Emax_Finite, Emax, beta);
      syst_s_Aii_tm = fabs( Br_s_sigma_Aii_tm.ave() - Br_sigma_Aii_tm.ave());
      Br_sigma_Aii_tm= Br_sigma_Aii_tm.ave() + (Br_sigma_Aii_tm-Br_sigma_Aii_tm.ave())*sqrt( pow(Br_sigma_Aii_tm.err(),2) + pow(syst_s_Aii_tm,2))/Br_sigma_Aii_tm.err();
      auto end = chrono::system_clock::now();
      cout<<"node: "<<_hostname<<", rank: "<<rank<<", thread_id: "<<omp_get_thread_num()<<" core-id: "<<sched_getcpu()<<endl<<flush;
      chrono::duration<double> elapsed_seconds = end-start;
      double time_Aii_tm= elapsed_seconds.count();
      if(tau_strange_verbosity_lev) {
	cout<<endl<<flush;
	cout<<"Elapsed time[Aii_tm, sigma: "<<s<<", Ens: "<<ls_data_tm_VKVK.Tag[iens]<<", rank: "<<rank<<", #thread="<<omp_get_thread_num()<<"] : "<<time_Aii_tm<<" s"<<endl<<flush;
      }
      else cout<<"."<<flush;

      //double E0_A_sp_old= E0_l;
      //E0_A_sp = (ls_data_tm_VKVK.Tag[iens] == "cC211a.06.80")?0.9*(m_kappa +  MPiPhys  ):E0_A_sp_old; 
      start = chrono::system_clock::now();
      if( (beta > 2) && Use_Customized_plateaux_strange) mult=  Customized_plateaux_tau_spectre_strange( beta, Emax, "Aii", "OS" , s, ls_data_tm_VKVK.Tag[iens] );
      Br_sigma_Aii_OS = Get_Laplace_transfo(  0.0,  s, E0_A_sp*a_distr.ave(),  T, tmax_OS_1_Aii, prec, SM_TYPE_1,K1, Aii_OS, syst_Aii_OS, mult, lAii_OS, MODE, "OS", "Aii_strange_"+ls_data_tm_VKVK.Tag[iens], 1e-3,0,resc_GeV*Za*Za, 0.0, "tau_decay", cov_Ak_OS, fake_func,0, fake_func_d ,  Is_Emax_Finite, Emax, beta);
      Br_s_sigma_Aii_OS = Get_Laplace_transfo(  0.0,  s, E0_A_sp*a_distr.ave(),  T, tmax_OS_1_Aii, prec, SM_TYPE_1,K1_shifted, Aii_OS, syst_s_Aii_OS, mult, slAii_OS, MODE, "OS", "Aii_s_strange_"+ls_data_tm_VKVK.Tag[iens], 1e-3,0,resc_GeV*Za*Za, 0.0, "tau_decay", cov_Ak_OS, fake_func,0, fake_func_d ,  Is_Emax_Finite, Emax, beta);
      //distr_t preco_Aii_OS = Get_Laplace_transfo(  0.0,  s, E0_A_sp*a_distr.ave(),  T, tmax_OS_1_Aii, prec, SM_TYPE_1,K1, Aii_OS - Aii_OS_gs.ave(), syst_Aii_OS, mult, lAii_OS, MODE, "OS", "preco_Aii_strange_"+ls_data_OS_VKVK.Tag[iens], 1e-3,0, resc_GeV*Za*Za, A_K1_OS.ave()*K1_dub(m_K1_OS.ave(), s), "tau_decay", cov_Ak_OS, fake_func,0, fake_func_d ,  Is_Emax_Finite, Emax, beta);
      syst_s_Aii_OS = fabs( Br_s_sigma_Aii_OS.ave() - Br_sigma_Aii_OS.ave());
      Br_sigma_Aii_OS= Br_sigma_Aii_OS.ave() + (Br_sigma_Aii_OS-Br_sigma_Aii_OS.ave())*sqrt( pow(Br_sigma_Aii_OS.err(),2) + pow(syst_s_Aii_OS,2))/Br_sigma_Aii_OS.err();

      //E0_A_sp= E0_A_sp_old;
     
      end = chrono::system_clock::now();
      elapsed_seconds = end-start;
      double time_Aii_OS= elapsed_seconds.count();
      if(tau_strange_verbosity_lev) {
	cout<<endl<<flush;
	cout<<"Elapsed time[Aii_OS, sigma: "<<s<<", Ens: "<<ls_data_tm_VKVK.Tag[iens]<<", rank: "<<rank<<", #thread="<<omp_get_thread_num()<<"] : "<<time_Aii_OS<<" s"<<endl<<flush;
      }
      else cout<<"."<<flush;
      

     
      start = chrono::system_clock::now();
      if( (beta > 2) && Use_Customized_plateaux_strange) mult=  Customized_plateaux_tau_spectre_strange( beta, Emax, "Vii", "tm" , s, ls_data_tm_VKVK.Tag[iens] );
      Br_sigma_Vii_tm = Get_Laplace_transfo(  0.0,  s, E0_sp*a_distr.ave(),  T, tmax_tm_1_Vii, prec, SM_TYPE_1,K1, Vii_tm, syst_Vii_tm, mult, lVii_tm, MODE, "tm", "Vii_strange_"+ls_data_tm_VKVK.Tag[iens], 1e-3,0, resc_GeV*Za*Za, 0.0, "tau_decay", cov_Vk_tm, fake_func,0, fake_func_d , Is_Emax_Finite, Emax, beta);
      Br_s_sigma_Vii_tm = Get_Laplace_transfo(  0.0,  s, E0_sp*a_distr.ave(),  T, tmax_tm_1_Vii, prec, SM_TYPE_1,K1_shifted, Vii_tm, syst_s_Vii_tm, mult, slVii_tm, MODE, "tm", "Vii_s_strange_"+ls_data_tm_VKVK.Tag[iens], 1e-3,0, resc_GeV*Za*Za, 0.0, "tau_decay", cov_Vk_tm, fake_func,0, fake_func_d , Is_Emax_Finite, Emax, beta);
      //distr_t preco_Vii_tm = Get_Laplace_transfo(  0.0,  s, E0_sp*a_distr.ave(),  T, tmax_tm_1_Vii, prec, SM_TYPE_1,K1, Vii_tm - Vii_tm_gs.ave(), syst_Vii_tm, mult, lVii_tm, MODE, "tm", "preco_Vii_strange_"+ls_data_tm_VKVK.Tag[iens], 1e-3,0, resc_GeV*Za*Za, A_Kstar_tm.ave()*K1_dub(m_Kstar_tm.ave(), s), "tau_decay", cov_Ak_tm, fake_func,0, fake_func_d ,  Is_Emax_Finite, Emax, beta);
      syst_s_Vii_tm = fabs( Br_s_sigma_Vii_tm.ave() - Br_sigma_Vii_tm.ave());
      Br_sigma_Vii_tm= Br_sigma_Vii_tm.ave() + (Br_sigma_Vii_tm-Br_sigma_Vii_tm.ave())*sqrt( pow(Br_sigma_Vii_tm.err(),2) + pow(syst_s_Vii_tm,2))/Br_sigma_Vii_tm.err();
      end = chrono::system_clock::now();
      elapsed_seconds = end-start;
      double time_Vii_tm= elapsed_seconds.count();
      if(tau_strange_verbosity_lev) {
	cout<<endl<<flush;
	cout<<"Elapsed time[Vii_tm, sigma: "<<s<<", Ens: "<<ls_data_tm_VKVK.Tag[iens]<<", rank: "<<rank<<", #thread="<<omp_get_thread_num()<<"] : "<<time_Vii_tm<<" s"<<endl<<flush;
      }
      else cout<<"."<<flush;

     
      start = chrono::system_clock::now();
      if( (beta > 2) && Use_Customized_plateaux_strange) mult=  Customized_plateaux_tau_spectre_strange( beta, Emax,  "Vii", "OS" , s, ls_data_tm_VKVK.Tag[iens] );
      Br_sigma_Vii_OS = Get_Laplace_transfo(  0.0,  s, E0_sp*a_distr.ave(),  T, tmax_OS_1_Vii, prec, SM_TYPE_1,K1, Vii_OS, syst_Vii_OS, mult, lVii_OS, MODE, "OS", "Vii_strange_"+ls_data_tm_VKVK.Tag[iens],1e-3,0, resc_GeV*Zv*Zv, 0.0, "tau_decay", cov_Vk_OS, fake_func,0, fake_func_d  ,  Is_Emax_Finite, Emax, beta);
      Br_s_sigma_Vii_OS = Get_Laplace_transfo(  0.0,  s, E0_sp*a_distr.ave(),  T, tmax_OS_1_Vii, prec, SM_TYPE_1,K1_shifted, Vii_OS, syst_s_Vii_OS, mult, slVii_OS, MODE, "OS", "Vii_s_strange_"+ls_data_tm_VKVK.Tag[iens],1e-3,0, resc_GeV*Zv*Zv, 0.0, "tau_decay", cov_Vk_OS, fake_func,0, fake_func_d ,  Is_Emax_Finite, Emax, beta);
      //distr_t preco_Vii_OS = Get_Laplace_transfo(  0.0,  s, E0_sp*a_distr.ave(),  T, tmax_OS_1_Vii, prec, SM_TYPE_1,K1, Vii_OS - Vii_OS_gs.ave(), syst_Vii_OS, mult, lVii_OS, MODE, "OS", "preco_Vii_strange_"+ls_data_OS_VKVK.Tag[iens], 1e-3,0, resc_GeV*Zv*Zv, A_Kstar_OS.ave()*K1_dub(m_Kstar_OS.ave(), s), "tau_decay", cov_Ak_OS, fake_func,0, fake_func_d ,  Is_Emax_Finite, Emax, beta);
      syst_s_Vii_OS = fabs( Br_s_sigma_Vii_OS.ave() - Br_sigma_Vii_OS.ave());
      Br_sigma_Vii_OS= Br_sigma_Vii_OS.ave() + (Br_sigma_Vii_OS-Br_sigma_Vii_OS.ave())*sqrt( pow(Br_sigma_Vii_OS.err(),2) + pow(syst_s_Vii_OS,2))/Br_sigma_Vii_OS.err();
      end = chrono::system_clock::now();
      elapsed_seconds = end-start;
      double time_Vii_OS= elapsed_seconds.count();
      if(tau_strange_verbosity_lev) {
	cout<<endl<<flush;
	cout<<"Elapsed time[Vii_OS, sigma: "<<s<<", Ens: "<<ls_data_tm_VKVK.Tag[iens]<<", rank: "<<rank<<", #thread="<<omp_get_thread_num()<<"] : "<<time_Vii_OS<<" s"<<endl<<flush;
      }
      else cout<<"."<<flush;



      start = chrono::system_clock::now();
      if( (beta > 2) && Use_Customized_plateaux_strange) mult=  Customized_plateaux_tau_spectre_strange( beta, Emax, "A0", "tm" , s, ls_data_tm_VKVK.Tag[iens] );
      Br_sigma_A0_tm = Get_Laplace_transfo(  0.0,  s, E0_l*a_distr.ave(),  T, tmax_tm_A0, prec, SM_TYPE_0,K0, A0_tm, syst_A0_tm, mult, lA0_tm, MODE, "tm", "A0_strange_"+ls_data_tm_VKVK.Tag[iens], 1e-3,0, resc_GeV*Zv*Zv, 0.0, "tau_decay", cov_A0_tm, fake_func,0, fake_func_d ,  Is_Emax_Finite, Emax, beta );
      Br_s_sigma_A0_tm = Get_Laplace_transfo(  0.0,  s, E0_l*a_distr.ave(),  T, tmax_tm_A0, prec, SM_TYPE_0,K0_shifted, A0_tm, syst_s_A0_tm, mult, slA0_tm, MODE, "tm", "A0_s_strange_"+ls_data_tm_VKVK.Tag[iens], 1e-3,0, resc_GeV*Zv*Zv, 0.0, "tau_decay", cov_A0_tm, fake_func,0, fake_func_d ,  Is_Emax_Finite, Emax, beta );
      syst_s_A0_tm = fabs( Br_s_sigma_A0_tm.ave() - Br_sigma_A0_tm.ave());
      Br_sigma_A0_tm= Br_sigma_A0_tm.ave() + (Br_sigma_A0_tm-Br_sigma_A0_tm.ave())*sqrt( pow(Br_sigma_A0_tm.err(),2) + pow(syst_s_A0_tm,2))/Br_sigma_A0_tm.err();
      end = chrono::system_clock::now();
      elapsed_seconds = end-start;
      double time_A0_tm= elapsed_seconds.count();
      if(tau_strange_verbosity_lev) {
	cout<<endl<<flush;
	cout<<"Elapsed time[A0_tm, sigma: "<<s<<", Ens: "<<ls_data_tm_VKVK.Tag[iens]<<", rank: "<<rank<<", #thread="<<omp_get_thread_num()<<"] : "<<time_A0_tm<<" s"<<endl<<flush;
      }
      else cout<<"."<<flush;
      

      start = chrono::system_clock::now();
      if( (beta > 2) && Use_Customized_plateaux_strange) mult=  Customized_plateaux_tau_spectre_strange( beta, Emax, "A0", "OS" , s, ls_data_tm_VKVK.Tag[iens] );
      Br_sigma_A0_OS = Get_Laplace_transfo(  0.0,  s, E0_l*a_distr.ave(),  T, tmax_OS_A0, prec, SM_TYPE_0,K0, A0_OS, syst_A0_OS, mult, lA0_OS, MODE, "OS", "A0_strange_"+ls_data_tm_VKVK.Tag[iens], 1e-3,0, resc_GeV*Za*Za, 0.0, "tau_decay", cov_A0_OS, fake_func,0, fake_func_d ,  Is_Emax_Finite, Emax, beta );
      Br_s_sigma_A0_OS = Get_Laplace_transfo(  0.0,  s, E0_l*a_distr.ave(),  T, tmax_OS_A0, prec, SM_TYPE_0,K0_shifted, A0_OS, syst_s_A0_OS, mult, slA0_OS, MODE, "OS", "A0_s_strange_"+ls_data_tm_VKVK.Tag[iens], 1e-3,0, resc_GeV*Za*Za, 0.0, "tau_decay", cov_A0_OS, fake_func,0, fake_func_d ,  Is_Emax_Finite, Emax, beta );
      syst_s_A0_OS = fabs( Br_s_sigma_A0_OS.ave() - Br_sigma_A0_OS.ave());
      Br_sigma_A0_OS= Br_sigma_A0_OS.ave() + (Br_sigma_A0_OS-Br_sigma_A0_OS.ave())*sqrt( pow(Br_sigma_A0_OS.err(),2) + pow(syst_s_A0_OS,2))/Br_sigma_A0_OS.err();
      end = chrono::system_clock::now();
      elapsed_seconds = end-start;
      double time_A0_OS= elapsed_seconds.count();
      if(tau_strange_verbosity_lev) {
	cout<<endl<<flush;
	cout<<"Elapsed time[A0_OS, sigma: "<<s<<", Ens: "<<ls_data_tm_VKVK.Tag[iens]<<", rank: "<<rank<<", #thread="<<omp_get_thread_num()<<"] : "<<time_A0_OS<<" s"<<endl<<flush;
      }
      else cout<<"."<<flush;



      start = chrono::system_clock::now();
      if( (beta > 2) && Use_Customized_plateaux_strange) mult=  Customized_plateaux_tau_spectre_strange( beta, Emax, "V0", "tm" , s, ls_data_tm_VKVK.Tag[iens] );
      Br_sigma_V0_tm = Get_Laplace_transfo(  0.0,  s, E0_l*a_distr.ave(),  T, tmax_tm_V0, prec, SM_TYPE_0,K0, V0_tm, syst_V0_tm, mult, lV0_tm, MODE, "tm", "V0_strange_"+ls_data_tm_VKVK.Tag[iens], 1e-3,0, resc_GeV*Za*Za, 0.0, "tau_decay", cov_V0_tm, fake_func,0, fake_func_d ,  Is_Emax_Finite, Emax, beta );
      Br_s_sigma_V0_tm = Get_Laplace_transfo(  0.0,  s, E0_l*a_distr.ave(),  T, tmax_tm_V0, prec, SM_TYPE_0,K0_shifted, V0_tm, syst_s_V0_tm, mult, slV0_tm, MODE, "tm", "V0_s_strange_"+ls_data_tm_VKVK.Tag[iens], 1e-3,0, resc_GeV*Za*Za, 0.0, "tau_decay", cov_V0_tm, fake_func,0, fake_func_d ,  Is_Emax_Finite, Emax, beta );

      syst_s_V0_tm = fabs( Br_s_sigma_V0_tm.ave() - Br_sigma_V0_tm.ave());
      Br_sigma_V0_tm= Br_sigma_V0_tm.ave() + (Br_sigma_V0_tm-Br_sigma_V0_tm.ave())*sqrt( pow(Br_sigma_V0_tm.err(),2) + pow(syst_s_V0_tm,2))/Br_sigma_V0_tm.err();
      end = chrono::system_clock::now();
      elapsed_seconds = end-start;
      double time_V0_tm= elapsed_seconds.count();
      if(tau_strange_verbosity_lev) {
	cout<<endl<<flush;
	cout<<"Elapsed time[V0_tm, sigma: "<<s<<", Ens: "<<ls_data_tm_VKVK.Tag[iens]<<", rank: "<<rank<<", #thread="<<omp_get_thread_num()<<"] : "<<time_V0_tm<<" s"<<endl<<flush;
      }
      else cout<<"."<<flush;
      

      start = chrono::system_clock::now();
      if( (beta > 2) && Use_Customized_plateaux_strange) mult=  Customized_plateaux_tau_spectre_strange( beta, Emax, "V0", "OS" , s, ls_data_tm_VKVK.Tag[iens] );
      Br_sigma_V0_OS = Get_Laplace_transfo(  0.0,  s, E0_l*a_distr.ave(),  T, tmax_OS_V0, prec, SM_TYPE_0,K0, V0_OS, syst_V0_OS, mult, lV0_OS, MODE, "OS", "V0_strange_"+ls_data_tm_VKVK.Tag[iens], 1e-3,0, resc_GeV*Zv*Zv, 0.0, "tau_decay", cov_V0_OS, fake_func,0, fake_func_d ,  Is_Emax_Finite, Emax, beta );
      Br_s_sigma_V0_OS = Get_Laplace_transfo(  0.0,  s, E0_l*a_distr.ave(),  T, tmax_OS_V0, prec, SM_TYPE_0,K0_shifted, V0_OS, syst_s_V0_OS, mult, slV0_OS, MODE, "OS", "V0_s_strange_"+ls_data_tm_VKVK.Tag[iens], 1e-3,0, resc_GeV*Zv*Zv, 0.0, "tau_decay", cov_V0_OS, fake_func,0, fake_func_d ,  Is_Emax_Finite, Emax, beta );
      syst_s_V0_OS = fabs( Br_s_sigma_V0_OS.ave() - Br_sigma_V0_OS.ave());
      Br_sigma_V0_OS= Br_sigma_V0_OS.ave() + (Br_sigma_V0_OS-Br_sigma_V0_OS.ave())*sqrt( pow(Br_sigma_V0_OS.err(),2) + pow(syst_s_V0_OS,2))/Br_sigma_V0_OS.err();
      end = chrono::system_clock::now();
      elapsed_seconds = end-start;
      double time_V0_OS= elapsed_seconds.count();
      if(tau_strange_verbosity_lev) {
	cout<<endl<<flush;
	cout<<"Elapsed time[V0_OS, sigma: "<<s<<", Ens: "<<ls_data_tm_VKVK.Tag[iens]<<", rank: "<<rank<<", #thread="<<omp_get_thread_num()<<"] : "<<time_V0_OS<<" s"<<endl<<flush;
      }
      else cout<<"."<<flush;


      
      
      start = chrono::system_clock::now();

      
      
      if( (beta > 2) && Use_Customized_plateaux_strange) mult=  Customized_plateaux_tau_spectre_strange( beta, Emax, "T", "tm" , s, ls_data_tm_VKVK.Tag[iens] );
      
      
      Br_sigma_T_tm = Get_Laplace_transfo(  0.0,  s, E0_sp*a_distr.ave(),  T, tmax_tm_T, prec, SM_TYPE_1 ,K1, T_tm, syst_T_tm, mult, lT_tm, MODE, "tm", "T_strange_"+ls_data_tm_VKVK.Tag[iens], 1e-3,0, resc_GeV*Za*Za/(Za.ave()*Za.ave()), 0.0, "tau_decay", cov_T_tm, fake_func,0, fake_func_d ,  Is_Emax_Finite, Emax, beta );
      Br_s_sigma_T_tm = Get_Laplace_transfo(  0.0,  s, E0_sp*a_distr.ave(),  T, tmax_tm_T, prec, SM_TYPE_1, K1_shifted, T_tm, syst_s_T_tm, mult, slT_tm, MODE, "tm", "T_s_strange_"+ls_data_tm_VKVK.Tag[iens], 1e-3,0, resc_GeV*Za*Za/(Za.ave()*Za.ave()), 0.0, "tau_decay", cov_T_tm, fake_func,0, fake_func_d ,  Is_Emax_Finite, Emax, beta );

      
      syst_s_T_tm = fabs( Br_s_sigma_T_tm.ave() - Br_sigma_T_tm.ave());
      Br_sigma_T_tm= Br_sigma_T_tm.ave() + (Br_sigma_T_tm-Br_sigma_T_tm.ave())*sqrt( pow(Br_sigma_T_tm.err(),2) + pow(syst_s_T_tm,2))/Br_sigma_T_tm.err();
      end = chrono::system_clock::now();
      elapsed_seconds = end-start;
      double time_T_tm= elapsed_seconds.count();
      if(tau_strange_verbosity_lev) {
	cout<<endl<<flush;
	cout<<"Elapsed time[T_tm, sigma: "<<s<<", Ens: "<<ls_data_tm_VKVK.Tag[iens]<<", rank: "<<rank<<", #thread="<<omp_get_thread_num()<<"] : "<<time_T_tm<<" s"<<endl<<flush;
      }
      else cout<<"."<<flush;


      start = chrono::system_clock::now();
      if( (beta > 2) && Use_Customized_plateaux_strange) mult=  Customized_plateaux_tau_spectre_strange( beta, Emax, "T", "OS" , s, ls_data_tm_VKVK.Tag[iens] );
      Br_sigma_T_OS = Get_Laplace_transfo(  0.0,  s, E0_sp*a_distr.ave(),  T, tmax_OS_T, prec, SM_TYPE_1,K1, T_OS, syst_T_OS, mult, lT_OS, MODE, "OS", "T_strange_"+ls_data_tm_VKVK.Tag[iens], 1e-3,0, resc_GeV*Zv*Zv/(Zv.ave()*Zv.ave()), 0.0, "tau_decay", cov_T_OS, fake_func,0, fake_func_d ,  Is_Emax_Finite, Emax, beta );

            
      Br_s_sigma_T_OS = Get_Laplace_transfo(  0.0,  s, E0_sp*a_distr.ave(),  T, tmax_OS_T, prec, SM_TYPE_1,K1_shifted, T_OS, syst_s_T_OS, mult, slT_OS, MODE, "OS", "T_s_strange_"+ls_data_tm_VKVK.Tag[iens], 1e-3,0, resc_GeV*Zv*Zv/(Zv.ave()*Zv.ave()), 0.0, "tau_decay", cov_T_OS, fake_func,0, fake_func_d ,  Is_Emax_Finite, Emax, beta );

      syst_s_T_OS = fabs( Br_s_sigma_T_OS.ave() - Br_sigma_T_OS.ave());
      Br_sigma_T_OS= Br_sigma_T_OS.ave() + (Br_sigma_T_OS-Br_sigma_T_OS.ave())*sqrt( pow(Br_sigma_T_OS.err(),2) + pow(syst_s_T_OS,2))/Br_sigma_T_OS.err();
      end = chrono::system_clock::now();
      elapsed_seconds = end-start;
      double time_T_OS= elapsed_seconds.count();
      if(tau_strange_verbosity_lev) {
	cout<<endl<<flush;
	cout<<"Elapsed time[T_OS, sigma: "<<s<<", Ens: "<<ls_data_tm_VKVK.Tag[iens]<<", rank: "<<rank<<", #thread="<<omp_get_thread_num()<<"] : "<<time_T_OS<<" s"<<endl<<flush;
      }
      else cout<<"."<<flush;



      start = chrono::system_clock::now();
      if( (beta > 2) && Use_Customized_plateaux_strange) mult=  Customized_plateaux_tau_spectre_strange( beta, Emax, "L", "tm" , s, ls_data_tm_VKVK.Tag[iens] );
      Br_sigma_L_tm = Get_Laplace_transfo(  0.0,  s, E0_sp*a_distr.ave(),  T, tmax_tm_L, prec, SM_TYPE_0,K0, L_tm, syst_L_tm, mult, lL_tm, MODE, "tm", "L_strange_"+ls_data_tm_VKVK.Tag[iens], 1e-3,0, resc_GeV*Zv*Zv/(Zv.ave()*Zv.ave()), 0.0, "tau_decay", cov_L_tm, fake_func,0, fake_func_d ,  Is_Emax_Finite, Emax, beta );
      Br_s_sigma_L_tm = Get_Laplace_transfo(  0.0,  s, E0_sp*a_distr.ave(),  T, tmax_tm_L, prec, SM_TYPE_0,K0_shifted, L_tm, syst_s_L_tm, mult, slL_tm, MODE, "tm", "L_s_strange_"+ls_data_tm_VKVK.Tag[iens], 1e-3,0, resc_GeV*Zv*Zv/(Zv.ave()*Zv.ave()), 0.0, "tau_decay", cov_L_tm, fake_func,0, fake_func_d ,  Is_Emax_Finite, Emax, beta );

      syst_s_L_tm = fabs( Br_s_sigma_L_tm.ave() - Br_sigma_L_tm.ave());
      Br_sigma_L_tm= Br_sigma_L_tm.ave() + (Br_sigma_L_tm-Br_sigma_L_tm.ave())*sqrt( pow(Br_sigma_L_tm.err(),2) + pow(syst_s_L_tm,2))/Br_sigma_L_tm.err();
      end = chrono::system_clock::now();
      elapsed_seconds = end-start;
      double time_L_tm= elapsed_seconds.count();
      if(tau_strange_verbosity_lev) {
	cout<<endl<<flush;
	cout<<"Elapsed time[L_tm, sigma: "<<s<<", Ens: "<<ls_data_tm_VKVK.Tag[iens]<<", rank: "<<rank<<", #thread="<<omp_get_thread_num()<<"] : "<<time_L_tm<<" s"<<endl<<flush;
      }
      else cout<<"."<<flush;


      start = chrono::system_clock::now();
      if( (beta > 2) && Use_Customized_plateaux_strange) mult=  Customized_plateaux_tau_spectre_strange( beta, Emax, "L", "OS" , s, ls_data_tm_VKVK.Tag[iens] );
      Br_sigma_L_OS = Get_Laplace_transfo(  0.0,  s, E0_sp*a_distr.ave(),  T, tmax_OS_L, prec, SM_TYPE_0,K0, L_OS, syst_L_OS, mult, lL_OS, MODE, "OS", "L_strange_"+ls_data_tm_VKVK.Tag[iens], 1e-3,0, resc_GeV*Za*Za/(Za.ave()*Za.ave()), 0.0, "tau_decay", cov_L_OS, fake_func,0, fake_func_d ,  Is_Emax_Finite, Emax, beta );
      Br_s_sigma_L_OS = Get_Laplace_transfo(  0.0,  s, E0_sp*a_distr.ave(),  T, tmax_OS_L, prec, SM_TYPE_0,K0_shifted, L_OS, syst_s_L_OS, mult, slL_OS, MODE, "OS", "L_s_strange_"+ls_data_tm_VKVK.Tag[iens], 1e-3,0, resc_GeV*Za*Za/(Za.ave()*Za.ave()), 0.0, "tau_decay", cov_L_OS, fake_func,0, fake_func_d ,  Is_Emax_Finite, Emax, beta );

      syst_s_L_OS = fabs( Br_s_sigma_L_OS.ave() - Br_sigma_L_OS.ave());
      Br_sigma_L_OS= Br_sigma_L_OS.ave() + (Br_sigma_L_OS-Br_sigma_L_OS.ave())*sqrt( pow(Br_sigma_L_OS.err(),2) + pow(syst_s_L_OS,2))/Br_sigma_L_OS.err();
      end = chrono::system_clock::now();
      elapsed_seconds = end-start;
      double time_L_OS= elapsed_seconds.count();
      if(tau_strange_verbosity_lev) {
	cout<<endl<<flush;
	cout<<"Elapsed time[L_OS, sigma: "<<s<<", Ens: "<<ls_data_tm_VKVK.Tag[iens]<<", rank: "<<rank<<", #thread="<<omp_get_thread_num()<<"] : "<<time_L_OS<<" s"<<endl<<flush;
      }
      else cout<<"."<<flush;
      
      

    



      

      thread_times_tm[is] = make_tuple(omp_get_thread_num(), time_A0_tm, time_V0_tm,  time_Aii_tm, time_Vii_tm);
      thread_times_OS[is] = make_tuple(omp_get_thread_num(), time_A0_OS, time_V0_OS,  time_Aii_OS, time_Vii_OS);





      //###########################################################//
      //###########################################################//
      //###########################################################//
      //###########################################################//
      //##############                              ###############//
      //##############                              ###############//
      //##############                              ###############//
      //##############         HEAVIER              ###############//
      //##############                              ###############//
      //##############                              ###############//
      //##############                              ###############//
  
      start = chrono::system_clock::now();
      if( (beta > 2) && Use_Customized_plateaux_strange) mult=  Customized_plateaux_tau_spectre_strange( beta, Emax,  "Aii", "tm" , s, ls_H_data_tm_VKVK.Tag[iens] );
      Br_sigma_Aii_H_tm = Get_Laplace_transfo(  0.0,  s, E0_A_sp*a_distr.ave(),  T, tmax_H_tm_1_Aii, prec, SM_TYPE_1,K1, Aii_H_tm, syst_Aii_H_tm, mult, lAii_H_tm, MODE, "tm", "Aii_strange_H_"+ls_H_data_tm_VKVK.Tag[iens], 1e-3,0, resc_GeV*Zv*Zv, 0.0, "tau_decay", cov_Ak_H_tm, fake_func,0, fake_func_d ,  Is_Emax_Finite, Emax, beta);
      Br_s_sigma_Aii_H_tm = Get_Laplace_transfo(  0.0,  s, E0_A_sp*a_distr.ave(),  T, tmax_H_tm_1_Aii, prec, SM_TYPE_1,K1_shifted, Aii_H_tm, syst_s_Aii_H_tm, mult, slAii_H_tm, MODE, "tm", "Aii_s_strange_H_"+ls_H_data_tm_VKVK.Tag[iens], 1e-3,0, resc_GeV*Zv*Zv, 0.0, "tau_decay", cov_Ak_H_tm, fake_func,0, fake_func_d ,  Is_Emax_Finite, Emax, beta);
      syst_s_Aii_H_tm = fabs( Br_s_sigma_Aii_H_tm.ave() - Br_sigma_Aii_H_tm.ave());
      Br_sigma_Aii_H_tm= Br_sigma_Aii_H_tm.ave() + (Br_sigma_Aii_H_tm-Br_sigma_Aii_H_tm.ave())*sqrt( pow(Br_sigma_Aii_H_tm.err(),2) + pow(syst_s_Aii_H_tm,2))/Br_sigma_Aii_H_tm.err();
      end = chrono::system_clock::now();
      cout<<"node: "<<_hostname<<", rank: "<<rank<<", thread_id: "<<omp_get_thread_num()<<" core-id: "<<sched_getcpu()<<endl<<flush;
      elapsed_seconds = end-start;
      double time_Aii_H_tm= elapsed_seconds.count();
      if(tau_strange_verbosity_lev) {
	cout<<endl<<flush;
	cout<<"Elapsed time[Aii_H_tm, sigma: "<<s<<", Ens: "<<ls_H_data_tm_VKVK.Tag[iens]<<", rank: "<<rank<<", #thread="<<omp_get_thread_num()<<"] : "<<time_Aii_H_tm<<" s"<<endl<<flush;
      }
      else cout<<"."<<flush;

      //E0_A_sp_old= E0_l;
      //E0_A_sp = (ls_data_tm_VKVK.Tag[iens] == "cC211a.06.80")?0.9*(m_kappa +  MPiPhys  ):E0_A_sp_old; 
     
      start = chrono::system_clock::now();
      if( (beta > 2) && Use_Customized_plateaux_strange) mult=  Customized_plateaux_tau_spectre_strange( beta, Emax, "Aii", "OS" , s, ls_H_data_tm_VKVK.Tag[iens] );
      Br_sigma_Aii_H_OS = Get_Laplace_transfo(  0.0,  s, E0_A_sp*a_distr.ave(),  T, tmax_H_OS_1_Aii, prec, SM_TYPE_1,K1, Aii_H_OS, syst_Aii_H_OS, mult, lAii_H_OS, MODE, "OS", "Aii_strange_H_"+ls_H_data_tm_VKVK.Tag[iens], 1e-3,0,resc_GeV*Za*Za, 0.0, "tau_decay", cov_Ak_H_OS, fake_func,0, fake_func_d ,  Is_Emax_Finite, Emax, beta);
      Br_s_sigma_Aii_H_OS = Get_Laplace_transfo(  0.0,  s, E0_A_sp*a_distr.ave(),  T, tmax_H_OS_1_Aii, prec, SM_TYPE_1,K1_shifted, Aii_H_OS, syst_s_Aii_H_OS, mult, slAii_H_OS, MODE, "OS", "Aii_s_strange_H_"+ls_H_data_tm_VKVK.Tag[iens], 1e-3,0,resc_GeV*Za*Za, 0.0, "tau_decay", cov_Ak_H_OS, fake_func,0, fake_func_d ,  Is_Emax_Finite, Emax, beta);
      syst_s_Aii_H_OS = fabs( Br_s_sigma_Aii_H_OS.ave() - Br_sigma_Aii_H_OS.ave());
      Br_sigma_Aii_H_OS= Br_sigma_Aii_H_OS.ave() + (Br_sigma_Aii_H_OS-Br_sigma_Aii_H_OS.ave())*sqrt( pow(Br_sigma_Aii_H_OS.err(),2) + pow(syst_s_Aii_H_OS,2))/Br_sigma_Aii_H_OS.err();
      end = chrono::system_clock::now();
      elapsed_seconds = end-start;
      double time_Aii_H_OS= elapsed_seconds.count();
      if(tau_strange_verbosity_lev) {
	cout<<endl<<flush;
	cout<<"Elapsed time[Aii_H_OS, sigma: "<<s<<", Ens: "<<ls_H_data_tm_VKVK.Tag[iens]<<", rank: "<<rank<<", #thread="<<omp_get_thread_num()<<"] : "<<time_Aii_H_OS<<" s"<<endl<<flush;
      }
      else cout<<"."<<flush;
      //E0_A_sp = E0_A_sp_old;

     
      start = chrono::system_clock::now();
      if( (beta > 2) && Use_Customized_plateaux_strange) mult=  Customized_plateaux_tau_spectre_strange( beta, Emax, "Vii", "tm" , s, ls_H_data_tm_VKVK.Tag[iens] );
      Br_sigma_Vii_H_tm = Get_Laplace_transfo(  0.0,  s, E0_sp*a_distr.ave(),  T, tmax_H_tm_1_Vii, prec, SM_TYPE_1,K1, Vii_H_tm, syst_Vii_H_tm, mult, lVii_H_tm, MODE, "tm", "Vii_strange_H_"+ls_H_data_tm_VKVK.Tag[iens], 1e-3,0, resc_GeV*Za*Za, 0.0, "tau_decay", cov_Vk_H_tm, fake_func,0, fake_func_d , Is_Emax_Finite, Emax, beta);
      Br_s_sigma_Vii_H_tm = Get_Laplace_transfo(  0.0,  s, E0_sp*a_distr.ave(),  T, tmax_H_tm_1_Vii, prec, SM_TYPE_1,K1_shifted, Vii_H_tm, syst_s_Vii_H_tm, mult, slVii_H_tm, MODE, "tm", "Vii_s_strange_H_"+ls_H_data_tm_VKVK.Tag[iens], 1e-3,0, resc_GeV*Za*Za, 0.0, "tau_decay", cov_Vk_H_tm, fake_func,0, fake_func_d , Is_Emax_Finite, Emax, beta);
      syst_s_Vii_H_tm = fabs( Br_s_sigma_Vii_H_tm.ave() - Br_sigma_Vii_H_tm.ave());
      Br_sigma_Vii_H_tm= Br_sigma_Vii_H_tm.ave() + (Br_sigma_Vii_H_tm-Br_sigma_Vii_H_tm.ave())*sqrt( pow(Br_sigma_Vii_H_tm.err(),2) + pow(syst_s_Vii_H_tm,2))/Br_sigma_Vii_H_tm.err();
      end = chrono::system_clock::now();
      elapsed_seconds = end-start;
      double time_Vii_H_tm= elapsed_seconds.count();
      if(tau_strange_verbosity_lev) {
	cout<<endl<<flush;
	cout<<"Elapsed time[Vii_H_tm, sigma: "<<s<<", Ens: "<<ls_H_data_tm_VKVK.Tag[iens]<<", rank: "<<rank<<", #thread="<<omp_get_thread_num()<<"] : "<<time_Vii_H_tm<<" s"<<endl<<flush;
      }
      else cout<<"."<<flush;

     
      start = chrono::system_clock::now();
      if( (beta > 2) && Use_Customized_plateaux_strange) mult=  Customized_plateaux_tau_spectre_strange( beta, Emax,  "Vii", "OS" , s, ls_H_data_tm_VKVK.Tag[iens] );
      Br_sigma_Vii_H_OS = Get_Laplace_transfo(  0.0,  s, E0_sp*a_distr.ave(),  T, tmax_H_OS_1_Vii, prec, SM_TYPE_1,K1, Vii_H_OS, syst_Vii_H_OS, mult, lVii_H_OS, MODE, "OS", "Vii_strange_H_"+ls_H_data_tm_VKVK.Tag[iens],1e-3,0, resc_GeV*Zv*Zv, 0.0, "tau_decay", cov_Vk_H_OS, fake_func,0, fake_func_d  ,  Is_Emax_Finite, Emax, beta);
      Br_s_sigma_Vii_H_OS = Get_Laplace_transfo(  0.0,  s, E0_sp*a_distr.ave(),  T, tmax_H_OS_1_Vii, prec, SM_TYPE_1,K1_shifted, Vii_H_OS, syst_s_Vii_H_OS, mult, slVii_H_OS, MODE, "OS", "Vii_s_strange_H_"+ls_H_data_tm_VKVK.Tag[iens],1e-3,0, resc_GeV*Zv*Zv, 0.0, "tau_decay", cov_Vk_H_OS, fake_func,0, fake_func_d ,  Is_Emax_Finite, Emax, beta);
      syst_s_Vii_H_OS = fabs( Br_s_sigma_Vii_H_OS.ave() - Br_sigma_Vii_H_OS.ave());
      Br_sigma_Vii_H_OS= Br_sigma_Vii_H_OS.ave() + (Br_sigma_Vii_H_OS-Br_sigma_Vii_H_OS.ave())*sqrt( pow(Br_sigma_Vii_H_OS.err(),2) + pow(syst_s_Vii_H_OS,2))/Br_sigma_Vii_H_OS.err();
      end = chrono::system_clock::now();
      elapsed_seconds = end-start;
      double time_Vii_H_OS= elapsed_seconds.count();
      if(tau_strange_verbosity_lev) {
	cout<<endl<<flush;
	cout<<"Elapsed time[Vii_H_OS, sigma: "<<s<<", Ens: "<<ls_H_data_tm_VKVK.Tag[iens]<<", rank: "<<rank<<", #thread="<<omp_get_thread_num()<<"] : "<<time_Vii_H_OS<<" s"<<endl<<flush;
      }
      else cout<<"."<<flush;



      start = chrono::system_clock::now();
      if( (beta > 2) && Use_Customized_plateaux_strange) mult=  Customized_plateaux_tau_spectre_strange( beta, Emax, "A0", "tm" , s, ls_H_data_tm_VKVK.Tag[iens] );
      Br_sigma_A0_H_tm = Get_Laplace_transfo(  0.0,  s, E0_l*a_distr.ave(),  T, tmax_H_tm_A0, prec, SM_TYPE_0,K0, A0_H_tm, syst_A0_H_tm, mult, lA0_H_tm, MODE, "tm", "A0_strange_H_"+ls_H_data_tm_VKVK.Tag[iens], 1e-3,0, resc_GeV*Zv*Zv, 0.0, "tau_decay", cov_A0_H_tm, fake_func,0, fake_func_d ,  Is_Emax_Finite, Emax, beta );
      Br_s_sigma_A0_H_tm = Get_Laplace_transfo(  0.0,  s, E0_l*a_distr.ave(),  T, tmax_H_tm_A0, prec, SM_TYPE_0,K0_shifted, A0_H_tm, syst_s_A0_H_tm, mult, slA0_H_tm, MODE, "tm", "A0_s_strange_H_"+ls_H_data_tm_VKVK.Tag[iens], 1e-3,0, resc_GeV*Zv*Zv, 0.0, "tau_decay", cov_A0_H_tm, fake_func,0, fake_func_d ,  Is_Emax_Finite, Emax, beta );
      syst_s_A0_H_tm = fabs( Br_s_sigma_A0_H_tm.ave() - Br_sigma_A0_H_tm.ave());
      Br_sigma_A0_H_tm= Br_sigma_A0_H_tm.ave() + (Br_sigma_A0_H_tm-Br_sigma_A0_H_tm.ave())*sqrt( pow(Br_sigma_A0_H_tm.err(),2) + pow(syst_s_A0_H_tm,2))/Br_sigma_A0_H_tm.err();
      end = chrono::system_clock::now();
      elapsed_seconds = end-start;
      double time_A0_H_tm= elapsed_seconds.count();
      if(tau_strange_verbosity_lev) {
	cout<<endl<<flush;
	cout<<"Elapsed time[A0_H_tm, sigma: "<<s<<", Ens: "<<ls_H_data_tm_VKVK.Tag[iens]<<", rank: "<<rank<<", #thread="<<omp_get_thread_num()<<"] : "<<time_A0_H_tm<<" s"<<endl<<flush;
      }
      else cout<<"."<<flush;
      

      start = chrono::system_clock::now();
      if( (beta > 2) && Use_Customized_plateaux_strange) mult=  Customized_plateaux_tau_spectre_strange( beta, Emax, "A0", "OS" , s, ls_H_data_tm_VKVK.Tag[iens] );
      Br_sigma_A0_H_OS = Get_Laplace_transfo(  0.0,  s, E0_l*a_distr.ave(),  T, tmax_H_OS_A0, prec, SM_TYPE_0,K0, A0_H_OS, syst_A0_H_OS, mult, lA0_H_OS, MODE, "OS", "A0_strange_H_"+ls_H_data_tm_VKVK.Tag[iens], 1e-3,0, resc_GeV*Za*Za, 0.0, "tau_decay", cov_A0_H_OS, fake_func,0, fake_func_d ,  Is_Emax_Finite, Emax, beta );
      Br_s_sigma_A0_H_OS = Get_Laplace_transfo(  0.0,  s, E0_l*a_distr.ave(),  T, tmax_H_OS_A0, prec, SM_TYPE_0,K0_shifted, A0_H_OS, syst_s_A0_H_OS, mult, slA0_H_OS, MODE, "OS", "A0_s_strange_H_"+ls_H_data_tm_VKVK.Tag[iens], 1e-3,0, resc_GeV*Za*Za, 0.0, "tau_decay", cov_A0_H_OS, fake_func,0, fake_func_d ,  Is_Emax_Finite, Emax, beta );
      syst_s_A0_H_OS = fabs( Br_s_sigma_A0_H_OS.ave() - Br_sigma_A0_H_OS.ave());
      Br_sigma_A0_H_OS= Br_sigma_A0_H_OS.ave() + (Br_sigma_A0_H_OS-Br_sigma_A0_H_OS.ave())*sqrt( pow(Br_sigma_A0_H_OS.err(),2) + pow(syst_s_A0_H_OS,2))/Br_sigma_A0_H_OS.err();
      end = chrono::system_clock::now();
      elapsed_seconds = end-start;
      double time_A0_H_OS= elapsed_seconds.count();
      if(tau_strange_verbosity_lev) {
	cout<<endl<<flush;
	cout<<"Elapsed time[A0_H_OS, sigma: "<<s<<", Ens: "<<ls_H_data_tm_VKVK.Tag[iens]<<", rank: "<<rank<<", #thread="<<omp_get_thread_num()<<"] : "<<time_A0_H_OS<<" s"<<endl<<flush;
      }
      else cout<<"."<<flush;



      start = chrono::system_clock::now();
      if( (beta > 2) && Use_Customized_plateaux_strange) mult=  Customized_plateaux_tau_spectre_strange( beta, Emax, "V0", "tm" , s, ls_H_data_tm_VKVK.Tag[iens] );
      Br_sigma_V0_H_tm = Get_Laplace_transfo(  0.0,  s, E0_l*a_distr.ave(),  T, tmax_H_tm_V0, prec, SM_TYPE_0,K0, V0_H_tm, syst_V0_H_tm, mult, lV0_H_tm, MODE, "tm", "V0_strange_H_"+ls_H_data_tm_VKVK.Tag[iens], 1e-3,0, resc_GeV*Za*Za, 0.0, "tau_decay", cov_V0_H_tm, fake_func,0, fake_func_d ,  Is_Emax_Finite, Emax, beta );
      Br_s_sigma_V0_H_tm = Get_Laplace_transfo(  0.0,  s, E0_l*a_distr.ave(),  T, tmax_H_tm_V0, prec, SM_TYPE_0,K0_shifted, V0_H_tm, syst_s_V0_H_tm, mult, slV0_H_tm, MODE, "tm", "V0_s_strange_H_"+ls_H_data_tm_VKVK.Tag[iens], 1e-3,0, resc_GeV*Za*Za, 0.0, "tau_decay", cov_V0_H_tm, fake_func,0, fake_func_d ,  Is_Emax_Finite, Emax, beta );
      syst_s_V0_H_tm = fabs( Br_s_sigma_V0_H_tm.ave() - Br_sigma_V0_H_tm.ave());
      Br_sigma_V0_H_tm= Br_sigma_V0_H_tm.ave() + (Br_sigma_V0_H_tm-Br_sigma_V0_H_tm.ave())*sqrt( pow(Br_sigma_V0_H_tm.err(),2) + pow(syst_s_V0_H_tm,2))/Br_sigma_V0_H_tm.err();
      end = chrono::system_clock::now();
      elapsed_seconds = end-start;
      double time_V0_H_tm= elapsed_seconds.count();
      if(tau_strange_verbosity_lev) {
	cout<<endl<<flush;
	cout<<"Elapsed time[V0_H_tm, sigma: "<<s<<", Ens: "<<ls_H_data_tm_VKVK.Tag[iens]<<", rank: "<<rank<<", #thread="<<omp_get_thread_num()<<"] : "<<time_V0_H_tm<<" s"<<endl<<flush;
      }
      else cout<<"."<<flush;
      

      start = chrono::system_clock::now();
      if( (beta > 2) && Use_Customized_plateaux_strange) mult=  Customized_plateaux_tau_spectre_strange( beta, Emax, "V0", "OS" , s, ls_H_data_tm_VKVK.Tag[iens] );
      Br_sigma_V0_H_OS = Get_Laplace_transfo(  0.0,  s, E0_l*a_distr.ave(),  T, tmax_H_OS_V0, prec, SM_TYPE_0,K0, V0_H_OS, syst_V0_H_OS, mult, lV0_H_OS, MODE, "OS", "V0_strange_H_"+ls_H_data_tm_VKVK.Tag[iens], 1e-3,0, resc_GeV*Zv*Zv, 0.0, "tau_decay", cov_V0_H_OS, fake_func,0, fake_func_d ,  Is_Emax_Finite, Emax, beta );
      Br_s_sigma_V0_H_OS = Get_Laplace_transfo(  0.0,  s, E0_l*a_distr.ave(),  T, tmax_H_OS_V0, prec, SM_TYPE_0,K0_shifted, V0_H_OS, syst_s_V0_H_OS, mult, slV0_H_OS, MODE, "OS", "V0_s_strange_H_"+ls_H_data_tm_VKVK.Tag[iens], 1e-3,0, resc_GeV*Zv*Zv, 0.0, "tau_decay", cov_V0_H_OS, fake_func,0, fake_func_d ,  Is_Emax_Finite, Emax, beta );
      syst_s_V0_H_OS = fabs( Br_s_sigma_V0_H_OS.ave() - Br_sigma_V0_H_OS.ave());
      Br_sigma_V0_H_OS= Br_sigma_V0_H_OS.ave() + (Br_sigma_V0_H_OS-Br_sigma_V0_H_OS.ave())*sqrt( pow(Br_sigma_V0_H_OS.err(),2) + pow(syst_s_V0_H_OS,2))/Br_sigma_V0_H_OS.err();
      end = chrono::system_clock::now();
      elapsed_seconds = end-start;
      double time_V0_H_OS= elapsed_seconds.count();
      if(tau_strange_verbosity_lev) {
	cout<<endl<<flush;
	cout<<"Elapsed time[V0_H_OS, sigma: "<<s<<", Ens: "<<ls_H_data_tm_VKVK.Tag[iens]<<", rank: "<<rank<<", #thread="<<omp_get_thread_num()<<"] : "<<time_V0_H_OS<<" s"<<endl<<flush;
      }
      else cout<<"."<<flush;


      

      start = chrono::system_clock::now();
      if( (beta > 2) && Use_Customized_plateaux_strange) mult=  Customized_plateaux_tau_spectre_strange( beta, Emax, "T", "tm" , s, ls_data_tm_VKVK.Tag[iens] );
      Br_sigma_T_H_tm = Get_Laplace_transfo(  0.0,  s, E0_sp*a_distr.ave(),  T, tmax_H_tm_T, prec, SM_TYPE_1,K1, T_H_tm, syst_T_H_tm, mult, lT_H_tm, MODE, "tm", "T_strange_H_"+ls_data_tm_VKVK.Tag[iens], 1e-3,0, resc_GeV*Za*Za/(Za.ave()*Za.ave()), 0.0, "tau_decay", cov_T_H_tm, fake_func,0, fake_func_d ,  Is_Emax_Finite, Emax, beta );
      Br_s_sigma_T_H_tm = Get_Laplace_transfo(  0.0,  s, E0_sp*a_distr.ave(),  T, tmax_H_tm_T, prec, SM_TYPE_1,K1_shifted, T_H_tm, syst_s_T_H_tm, mult, slT_H_tm, MODE, "tm", "T_s_strange_H_"+ls_data_tm_VKVK.Tag[iens], 1e-3,0, resc_GeV*Za*Za/(Za.ave()*Za.ave()), 0.0, "tau_decay", cov_T_H_tm, fake_func,0, fake_func_d ,  Is_Emax_Finite, Emax, beta );

      syst_s_T_H_tm = fabs( Br_s_sigma_T_H_tm.ave() - Br_sigma_T_H_tm.ave());
      Br_sigma_T_H_tm= Br_sigma_T_H_tm.ave() + (Br_sigma_T_H_tm-Br_sigma_T_H_tm.ave())*sqrt( pow(Br_sigma_T_H_tm.err(),2) + pow(syst_s_T_H_tm,2))/Br_sigma_T_H_tm.err();
      end = chrono::system_clock::now();
      elapsed_seconds = end-start;
      double time_T_H_tm= elapsed_seconds.count();
      if(tau_strange_verbosity_lev) {
	cout<<endl<<flush;
	cout<<"Elapsed time[T_H_tm, sigma: "<<s<<", Ens: "<<ls_data_tm_VKVK.Tag[iens]<<", rank: "<<rank<<", #thread="<<omp_get_thread_num()<<"] : "<<time_T_H_tm<<" s"<<endl<<flush;
      }
      else cout<<"."<<flush;


      start = chrono::system_clock::now();
      if( (beta > 2) && Use_Customized_plateaux_strange) mult=  Customized_plateaux_tau_spectre_strange( beta, Emax, "T", "OS" , s, ls_data_tm_VKVK.Tag[iens] );
      Br_sigma_T_H_OS = Get_Laplace_transfo(  0.0,  s, E0_sp*a_distr.ave(),  T, tmax_H_OS_T, prec, SM_TYPE_1,K1, T_H_OS, syst_T_H_OS, mult, lT_H_OS, MODE, "OS", "T_strange_H_"+ls_data_tm_VKVK.Tag[iens], 1e-3,0, resc_GeV*Zv*Zv/(Zv.ave()*Zv.ave()), 0.0, "tau_decay", cov_T_H_OS, fake_func,0, fake_func_d ,  Is_Emax_Finite, Emax, beta );
      Br_s_sigma_T_H_OS = Get_Laplace_transfo(  0.0,  s, E0_sp*a_distr.ave(),  T, tmax_H_OS_T, prec, SM_TYPE_1,K1_shifted, T_H_OS, syst_s_T_H_OS, mult, slT_H_OS, MODE, "OS", "T_s_strange_H_"+ls_data_tm_VKVK.Tag[iens], 1e-3,0, resc_GeV*Zv*Zv/(Zv.ave()*Zv.ave()), 0.0, "tau_decay", cov_T_H_OS, fake_func,0, fake_func_d ,  Is_Emax_Finite, Emax, beta );

      syst_s_T_H_OS = fabs( Br_s_sigma_T_H_OS.ave() - Br_sigma_T_H_OS.ave());
      Br_sigma_T_H_OS= Br_sigma_T_H_OS.ave() + (Br_sigma_T_H_OS-Br_sigma_T_H_OS.ave())*sqrt( pow(Br_sigma_T_H_OS.err(),2) + pow(syst_s_T_H_OS,2))/Br_sigma_T_H_OS.err();
      end = chrono::system_clock::now();
      elapsed_seconds = end-start;
      double time_T_H_OS= elapsed_seconds.count();
      if(tau_strange_verbosity_lev) {
	cout<<endl<<flush;
	cout<<"Elapsed time[T_H_OS, sigma: "<<s<<", Ens: "<<ls_data_tm_VKVK.Tag[iens]<<", rank: "<<rank<<", #thread="<<omp_get_thread_num()<<"] : "<<time_T_H_OS<<" s"<<endl<<flush;
      }
      else cout<<"."<<flush;



      start = chrono::system_clock::now();
      if( (beta > 2) && Use_Customized_plateaux_strange) mult=  Customized_plateaux_tau_spectre_strange( beta, Emax, "L", "tm" , s, ls_data_tm_VKVK.Tag[iens] );
      Br_sigma_L_H_tm = Get_Laplace_transfo(  0.0,  s, E0_sp*a_distr.ave(),  T, tmax_H_tm_L, prec, SM_TYPE_0,K0, L_H_tm, syst_L_H_tm, mult, lL_H_tm, MODE, "tm", "L_strange_H_"+ls_data_tm_VKVK.Tag[iens], 1e-3,0, resc_GeV*Zv*Zv/(Zv.ave()*Zv.ave()), 0.0, "tau_decay", cov_L_H_tm, fake_func,0, fake_func_d ,  Is_Emax_Finite, Emax, beta );
      Br_s_sigma_L_H_tm = Get_Laplace_transfo(  0.0,  s, E0_sp*a_distr.ave(),  T, tmax_H_tm_L, prec, SM_TYPE_0,K0_shifted, L_H_tm, syst_s_L_H_tm, mult, slL_H_tm, MODE, "tm", "L_s_strange_H_"+ls_data_tm_VKVK.Tag[iens], 1e-3,0, resc_GeV*Zv*Zv/(Zv.ave()*Zv.ave()), 0.0, "tau_decay", cov_L_H_tm, fake_func,0, fake_func_d ,  Is_Emax_Finite, Emax, beta );

      syst_s_L_H_tm = fabs( Br_s_sigma_L_H_tm.ave() - Br_sigma_L_H_tm.ave());
      Br_sigma_L_H_tm= Br_sigma_L_H_tm.ave() + (Br_sigma_L_H_tm-Br_sigma_L_H_tm.ave())*sqrt( pow(Br_sigma_L_H_tm.err(),2) + pow(syst_s_L_H_tm,2))/Br_sigma_L_H_tm.err();
      end = chrono::system_clock::now();
      elapsed_seconds = end-start;
      double time_L_H_tm= elapsed_seconds.count();
      if(tau_strange_verbosity_lev) {
	cout<<endl<<flush;
	cout<<"Elapsed time[L_H_tm, sigma: "<<s<<", Ens: "<<ls_data_tm_VKVK.Tag[iens]<<", rank: "<<rank<<", #thread="<<omp_get_thread_num()<<"] : "<<time_L_H_tm<<" s"<<endl<<flush;
      }
      else cout<<"."<<flush;


      start = chrono::system_clock::now();
      if( (beta > 2) && Use_Customized_plateaux_strange) mult=  Customized_plateaux_tau_spectre_strange( beta, Emax, "L", "OS" , s, ls_data_tm_VKVK.Tag[iens] );
      Br_sigma_L_H_OS = Get_Laplace_transfo(  0.0,  s, E0_sp*a_distr.ave(),  T, tmax_H_OS_L, prec, SM_TYPE_0,K0, L_H_OS, syst_L_H_OS, mult, lL_H_OS, MODE, "OS", "L_strange_H_"+ls_data_tm_VKVK.Tag[iens], 1e-3,0, resc_GeV*Za*Za/(Za.ave()*Za.ave()), 0.0, "tau_decay", cov_L_H_OS, fake_func,0, fake_func_d ,  Is_Emax_Finite, Emax, beta );
      Br_s_sigma_L_H_OS = Get_Laplace_transfo(  0.0,  s, E0_sp*a_distr.ave(),  T, tmax_H_OS_L, prec, SM_TYPE_0,K0_shifted, L_H_OS, syst_s_L_H_OS, mult, slL_H_OS, MODE, "OS", "L_s_strange_H_"+ls_data_tm_VKVK.Tag[iens], 1e-3,0, resc_GeV*Za*Za/(Za.ave()*Za.ave()), 0.0, "tau_decay", cov_L_H_OS, fake_func,0, fake_func_d ,  Is_Emax_Finite, Emax, beta );

      syst_s_L_H_OS = fabs( Br_s_sigma_L_H_OS.ave() - Br_sigma_L_H_OS.ave());
      Br_sigma_L_H_OS= Br_sigma_L_H_OS.ave() + (Br_sigma_L_H_OS-Br_sigma_L_H_OS.ave())*sqrt( pow(Br_sigma_L_H_OS.err(),2) + pow(syst_s_L_H_OS,2))/Br_sigma_L_H_OS.err();
      end = chrono::system_clock::now();
      elapsed_seconds = end-start;
      double time_L_H_OS= elapsed_seconds.count();
      if(tau_strange_verbosity_lev) {
	cout<<endl<<flush;
	cout<<"Elapsed time[L_H_OS, sigma: "<<s<<", Ens: "<<ls_data_tm_VKVK.Tag[iens]<<", rank: "<<rank<<", #thread="<<omp_get_thread_num()<<"] : "<<time_L_H_OS<<" s"<<endl<<flush;
      }
      else cout<<"."<<flush;


      

      thread_times_H_tm[is] = make_tuple(omp_get_thread_num(), time_A0_H_tm, time_V0_H_tm,  time_Aii_H_tm, time_Vii_H_tm);
      thread_times_H_OS[is] = make_tuple(omp_get_thread_num(), time_A0_H_OS, time_V0_H_OS,  time_Aii_H_OS, time_Vii_H_OS);


      //########### EXTRAPOLATE TO THE PHYSICAL STRANGE POINT ################
      
      //tm
      vector<distr_t> BR_LIST_Vii_tm({ Br_sigma_Vii_tm, Br_sigma_Vii_H_tm});
      vector<distr_t> BR_LIST_Aii_tm({ Br_sigma_Aii_tm, Br_sigma_Aii_H_tm});
      vector<distr_t> BR_LIST_A0_tm({ Br_sigma_A0_tm, Br_sigma_A0_H_tm});
      vector<distr_t> BR_LIST_V0_tm({ Br_sigma_V0_tm, Br_sigma_V0_H_tm});

      vector<distr_t> BR_LIST_T_tm({ Br_sigma_T_tm, Br_sigma_T_H_tm});
      vector<distr_t> BR_LIST_L_tm({ Br_sigma_L_tm, Br_sigma_L_H_tm});
     
     
      //OS
      vector<distr_t> BR_LIST_Vii_OS({ Br_sigma_Vii_OS, Br_sigma_Vii_H_OS});
      vector<distr_t> BR_LIST_Aii_OS({ Br_sigma_Aii_OS, Br_sigma_Aii_H_OS});
      vector<distr_t> BR_LIST_A0_OS({ Br_sigma_A0_OS, Br_sigma_A0_H_OS});
      vector<distr_t> BR_LIST_V0_OS({ Br_sigma_V0_OS, Br_sigma_V0_H_OS});

      vector<distr_t> BR_LIST_T_OS({ Br_sigma_T_OS, Br_sigma_T_H_OS});
      vector<distr_t> BR_LIST_L_OS({ Br_sigma_L_OS, Br_sigma_L_H_OS});

      //tm
      distr_t Br_sigma_Vii_tm_extr=  Obs_extrapolation_meson_mass(BR_LIST_Vii_tm, MMK2, Mk_iso*Mk_iso ,  "../data/tau_decay/"+Tag_reco_type+"/strange"  , "Vii_tm_extrapolation_"+ls_data_tm_VKVK.Tag[iens]+"_sigma_"+to_string_with_precision(sigma_list_strange[is],3)+".dat",  UseJack, "SPLINE" );
      distr_t Br_sigma_Aii_tm_extr=  Obs_extrapolation_meson_mass(BR_LIST_Aii_tm, MMK2, Mk_iso*Mk_iso ,  "../data/tau_decay/"+Tag_reco_type+"/strange"  , "Aii_tm_extrapolation_"+ls_data_tm_VKVK.Tag[iens]+"_sigma_"+to_string_with_precision(sigma_list_strange[is],3)+".dat",  UseJack, "SPLINE" );
      distr_t Br_sigma_A0_tm_extr=  Obs_extrapolation_meson_mass(BR_LIST_A0_tm, MMK2, Mk_iso*Mk_iso ,  "../data/tau_decay/"+Tag_reco_type+"/strange"  , "A0_tm_extrapolation_"+ls_data_tm_VKVK.Tag[iens]+"_sigma_"+to_string_with_precision(sigma_list_strange[is],3)+".dat",  UseJack, "SPLINE" );
      distr_t Br_sigma_V0_tm_extr=  Obs_extrapolation_meson_mass(BR_LIST_V0_tm, MMK2, Mk_iso*Mk_iso ,  "../data/tau_decay/"+Tag_reco_type+"/strange"  , "V0_tm_extrapolation_"+ls_data_tm_VKVK.Tag[iens]+"_sigma_"+to_string_with_precision(sigma_list_strange[is],3)+".dat",  UseJack, "SPLINE" );
      //OS
      distr_t Br_sigma_Vii_OS_extr=  Obs_extrapolation_meson_mass(BR_LIST_Vii_OS, MMK2, Mk_iso*Mk_iso ,  "../data/tau_decay/"+Tag_reco_type+"/strange"  , "Vii_OS_extrapolation_"+ls_data_OS_VKVK.Tag[iens]+"_sigma_"+to_string_with_precision(sigma_list_strange[is],3)+".dat",  UseJack, "SPLINE" );
      distr_t Br_sigma_Aii_OS_extr=  Obs_extrapolation_meson_mass(BR_LIST_Aii_OS, MMK2, Mk_iso*Mk_iso ,  "../data/tau_decay/"+Tag_reco_type+"/strange"  , "Aii_OS_extrapolation_"+ls_data_OS_VKVK.Tag[iens]+"_sigma_"+to_string_with_precision(sigma_list_strange[is],3)+".dat",  UseJack, "SPLINE" );
      distr_t Br_sigma_A0_OS_extr=  Obs_extrapolation_meson_mass(BR_LIST_A0_OS, MMK2, Mk_iso*Mk_iso ,  "../data/tau_decay/"+Tag_reco_type+"/strange"  , "A0_OS_extrapolation_"+ls_data_OS_VKVK.Tag[iens]+"_sigma_"+to_string_with_precision(sigma_list_strange[is],3)+".dat",  UseJack, "SPLINE" );
      distr_t Br_sigma_V0_OS_extr=  Obs_extrapolation_meson_mass(BR_LIST_V0_OS, MMK2, Mk_iso*Mk_iso ,  "../data/tau_decay/"+Tag_reco_type+"/strange"  , "V0_OS_extrapolation_"+ls_data_OS_VKVK.Tag[iens]+"_sigma_"+to_string_with_precision(sigma_list_strange[is],3)+".dat",  UseJack, "SPLINE" );

      distr_t Br_sigma_T_tm_extr=  Obs_extrapolation_meson_mass(BR_LIST_T_tm, MMK2, Mk_iso*Mk_iso ,  "../data/tau_decay/"+Tag_reco_type+"/strange"  , "T_tm_extrapolation_"+ls_data_OS_VKVK.Tag[iens]+"_sigma_"+to_string_with_precision(sigma_list_strange[is],3)+".dat",  UseJack, "SPLINE" );
      distr_t Br_sigma_L_tm_extr=  Obs_extrapolation_meson_mass(BR_LIST_L_tm, MMK2, Mk_iso*Mk_iso ,  "../data/tau_decay/"+Tag_reco_type+"/strange"  , "L_tm_extrapolation_"+ls_data_OS_VKVK.Tag[iens]+"_sigma_"+to_string_with_precision(sigma_list_strange[is],3)+".dat",  UseJack, "SPLINE" );

      distr_t Br_sigma_T_OS_extr=  Obs_extrapolation_meson_mass(BR_LIST_T_OS, MMK2, Mk_iso*Mk_iso ,  "../data/tau_decay/"+Tag_reco_type+"/strange"  , "T_OS_extrapolation_"+ls_data_OS_VKVK.Tag[iens]+"_sigma_"+to_string_with_precision(sigma_list_strange[is],3)+".dat",  UseJack, "SPLINE" );
      distr_t Br_sigma_L_OS_extr=  Obs_extrapolation_meson_mass(BR_LIST_L_OS, MMK2, Mk_iso*Mk_iso ,  "../data/tau_decay/"+Tag_reco_type+"/strange"  , "L_OS_extrapolation_"+ls_data_OS_VKVK.Tag[iens]+"_sigma_"+to_string_with_precision(sigma_list_strange[is],3)+".dat",  UseJack, "SPLINE" );
      


      
      pull_per_ens_tm_A0[iens][is] = max(lA0_tm, lA0_H_tm);
      pull_per_ens_OS_A0[iens][is] = max(lA0_OS, lA0_H_OS);

      pull_per_ens_tm_V0[iens][is] = max(lV0_tm, lV0_H_tm);
      pull_per_ens_OS_V0[iens][is] = max(lV0_OS, lV0_H_OS);

      pull_per_ens_tm_Aii[iens][is] = max(lAii_tm, lAii_H_tm);
      pull_per_ens_OS_Aii[iens][is] = max(lAii_OS, lAii_H_OS);

      pull_per_ens_tm_Vii[iens][is] = max(lVii_tm, lVii_H_tm);
      pull_per_ens_OS_Vii[iens][is] = max(lVii_OS, lVii_H_OS);

      pull_per_ens_tm_T[iens][is] = max(lT_tm, lT_H_tm);
      pull_per_ens_OS_T[iens][is] = max(lT_OS, lT_H_OS);
      
      pull_per_ens_tm_L[iens][is] = max(lL_tm, lL_H_tm);
      pull_per_ens_OS_L[iens][is] = max(lL_OS, lL_H_OS);
           
      syst_per_ens_tm_A0[iens][is] = Br_sigma_A0_tm_extr.ave()*max(    syst_A0_tm/Br_sigma_A0_tm.ave(),  syst_A0_H_tm/Br_sigma_A0_H_tm.ave() ); // syst_A0_tm;
      syst_per_ens_tm_V0[iens][is] =  Br_sigma_V0_tm_extr.ave()*max(   syst_V0_tm/Br_sigma_V0_tm.ave(), syst_V0_H_tm/Br_sigma_V0_H_tm.ave() ); // syst_V0_tm;
      syst_per_ens_tm_Ak[iens][is]=  Br_sigma_Aii_tm_extr.ave()*max(   syst_Aii_tm/Br_sigma_Aii_tm.ave(), syst_Aii_H_tm/Br_sigma_Aii_H_tm.ave() ); // syst_Aii_tm;
      syst_per_ens_tm_Vk[iens][is]=   Br_sigma_Vii_tm_extr.ave()*max(  syst_Vii_tm/Br_sigma_Vii_tm.ave(),  syst_Vii_H_tm/Br_sigma_Vii_H_tm.ave() ); // syst_Vii_tm;
      syst_per_ens_tm[iens][is]= sqrt( pow(syst_per_ens_tm_A0[iens][is],2)+ pow(syst_per_ens_tm_Ak[iens][is],2)+ pow(syst_per_ens_tm_Vk[iens][is],2) + pow(syst_per_ens_tm_V0[iens][is],2));

      syst_per_ens_tm_T[iens][is]=   Br_sigma_T_tm_extr.ave()*max(  syst_T_tm/Br_sigma_T_tm.ave(),  syst_T_H_tm/Br_sigma_T_H_tm.ave() ); // syst_T_tm;
      syst_per_ens_tm_L[iens][is]=   Br_sigma_L_tm_extr.ave()*max(  syst_L_tm/Br_sigma_L_tm.ave(),  syst_L_H_tm/Br_sigma_L_H_tm.ave() ); // syst_L_tm;
      
              
      syst_per_ens_OS_A0[iens][is] = Br_sigma_A0_OS_extr.ave()*max(  syst_A0_OS/Br_sigma_A0_OS.ave(),  syst_A0_H_OS/Br_sigma_A0_H_OS.ave() ); // syst_A0_OS;
      syst_per_ens_OS_V0[iens][is] =  Br_sigma_V0_OS_extr.ave()*max( syst_V0_OS/Br_sigma_V0_OS.ave(),  syst_V0_H_OS/Br_sigma_V0_H_OS.ave() ); // syst_V0_OS;
      syst_per_ens_OS_Ak[iens][is]=  Br_sigma_Aii_OS_extr.ave()*max( syst_Aii_OS/Br_sigma_Aii_OS.ave(), syst_Aii_H_OS/Br_sigma_Aii_H_OS.ave() ); // syst_Aii_OS;
      syst_per_ens_OS_Vk[iens][is]=   Br_sigma_Vii_OS_extr.ave()*max(  syst_Vii_OS/Br_sigma_Vii_OS.ave(),  syst_Vii_H_OS/Br_sigma_Vii_H_OS.ave() ); // syst_Vii_OS;
      syst_per_ens_OS[iens][is]=  sqrt( pow(syst_per_ens_OS_A0[iens][is],2)+ pow(syst_per_ens_OS_Ak[iens][is],2)+ pow(syst_per_ens_OS_Vk[iens][is],2) + pow(syst_per_ens_OS_V0[iens][is],2));

      syst_per_ens_OS_T[iens][is]=   Br_sigma_T_OS_extr.ave()*max(  syst_T_OS/Br_sigma_T_OS.ave(),  syst_T_H_OS/Br_sigma_T_H_OS.ave() ); // syst_T_OS;
      syst_per_ens_OS_L[iens][is]=   Br_sigma_L_OS_extr.ave()*max(  syst_L_OS/Br_sigma_L_OS.ave(),  syst_L_H_OS/Br_sigma_L_H_OS.ave() ); // syst_L_OS;

     
      
      distr_t Br_sigma_tm_extr = Br_sigma_Aii_tm_extr + Br_sigma_Vii_tm_extr + Br_sigma_A0_tm_extr + Br_sigma_V0_tm_extr;
      distr_t Br_sigma_OS_extr = Br_sigma_Aii_OS_extr + Br_sigma_Vii_OS_extr + Br_sigma_A0_OS_extr + Br_sigma_V0_OS_extr;
      Br_tau_tm[iens].distr_list[is] = Br_sigma_tm_extr;
      Br_Aii_tau_tm[iens].distr_list[is] = Br_sigma_Aii_tm_extr;
      Br_Vii_tau_tm[iens].distr_list[is] = Br_sigma_Vii_tm_extr;
      Br_A0_tau_tm[iens].distr_list[is] = Br_sigma_A0_tm_extr;
      Br_V0_tau_tm[iens].distr_list[is] = Br_sigma_V0_tm_extr;
      Br_tau_OS[iens].distr_list[is] = Br_sigma_OS_extr;
      Br_Aii_tau_OS[iens].distr_list[is] = Br_sigma_Aii_OS_extr;
      Br_Vii_tau_OS[iens].distr_list[is] = Br_sigma_Vii_OS_extr;
      Br_A0_tau_OS[iens].distr_list[is] = Br_sigma_A0_OS_extr;
      Br_V0_tau_OS[iens].distr_list[is] = Br_sigma_V0_OS_extr;

      Br_T_tau_tm[iens].distr_list[is] = Br_sigma_T_tm_extr;
      Br_L_tau_tm[iens].distr_list[is] = Br_sigma_L_tm_extr;

      Br_T_tau_OS[iens].distr_list[is] = Br_sigma_T_OS_extr;
      Br_L_tau_OS[iens].distr_list[is] = Br_sigma_L_OS_extr;


      cout<<"Ensemble: "<<ls_data_tm_VKVK.Tag[iens]<<", sigma: "<<s<<" completed!"<<endl<<flush;
     
   
     
    }
    
    cout<<endl;
    cout<<"Finished ensemble: "<<ls_data_tm_VKVK.Tag[iens]<<"########################"<<endl<<flush;
    cout<<"Summary of performances: "<<endl<<flush;
    cout<<"sigma #thread  A0    V0    Aii     Vii"<<endl<<flush;
    cout<<"- - - - - - - - - - - - - - - - - - - - - "<<endl<<flush;
    for(int is=0; is < (signed)sigma_list_strange.size();is++) {
      cout<<sigma_list_strange[is]<<", "<<get<0>(thread_times_tm[is])<<": "<<get<1>(thread_times_tm[is])<<" s, "<<get<2>(thread_times_tm[is])<<" s, "<<get<3>(thread_times_tm[is])<<" s, "<<get<4>(thread_times_tm[is])<<" s"<<endl<<flush;
      cout<<sigma_list_strange[is]<<", "<<get<0>(thread_times_OS[is])<<": "<<get<1>(thread_times_OS[is])<<" s, "<<get<2>(thread_times_OS[is])<<" s, "<<get<3>(thread_times_OS[is])<<" s, "<<get<4>(thread_times_OS[is])<<" s"<<endl<<flush;
      cout<<"- - - - - - - - - - - - - - - - - - - - "<<endl<<flush;
    }
    
    }

  }
   
    
  
  
  //print FK
  Print_To_File( {}, {a_distr_list.ave(), FK_list.ave(), FK_list.err() }, "../data/tau_decay/"+Tag_reco_type+"/strange/FK/FK.list", "", "");  
  
  if(!Skip_spectral_density_analysis_strange) {


  //Print to File
  for(int iens=0; iens<Nens;iens++) {

    
    
    Print_To_File({}, {sigma_list_strange, Br_tau_tm[iens].ave(), Br_tau_tm[iens].err(), syst_per_ens_tm[iens] , Br_tau_OS[iens].ave(), Br_tau_OS[iens].err(), syst_per_ens_OS[iens]}, "../data/tau_decay/"+Tag_reco_type+"/strange/Br/br_"+MODE+"_sm_func_mode_"+to_string(sm_func_mode)+"_"+ls_data_tm_VKVK.Tag[iens]+".dat", "", "#sigma Br[tm] Br[OS],  Rs(HFLAV): "+to_string_with_precision(Rs_HFLAV,4)+" +- "+to_string_with_precision(D_Rs_HFLAV,4));
    Print_To_File({}, {sigma_list_strange, Br_A0_tau_tm[iens].ave(), Br_A0_tau_tm[iens].err(), syst_per_ens_tm_A0[iens], Br_V0_tau_tm[iens].ave(), Br_V0_tau_tm[iens].err(), syst_per_ens_tm_V0[iens], Br_Aii_tau_tm[iens].ave(), Br_Aii_tau_tm[iens].err(), syst_per_ens_tm_Ak[iens], Br_Vii_tau_tm[iens].ave(), Br_Vii_tau_tm[iens].err(), syst_per_ens_tm_Vk[iens]}, "../data/tau_decay/"+Tag_reco_type+"/strange/Br/br_contrib_tm_"+MODE+"_sm_func_mode_"+to_string(sm_func_mode)+"_"+ls_data_tm_VKVK.Tag[iens]+".dat", "", "#sigma Br_A0 Br_V0 Br_Aii Br_Vii");
    Print_To_File({}, {sigma_list_strange, Br_A0_tau_OS[iens].ave(), Br_A0_tau_OS[iens].err(), syst_per_ens_OS_A0[iens], Br_V0_tau_OS[iens].ave(), Br_V0_tau_OS[iens].err(), syst_per_ens_OS_V0[iens], Br_Aii_tau_OS[iens].ave(), Br_Aii_tau_OS[iens].err(), syst_per_ens_OS_Ak[iens], Br_Vii_tau_OS[iens].ave(), Br_Vii_tau_OS[iens].err(), syst_per_ens_OS_Vk[iens]}, "../data/tau_decay/"+Tag_reco_type+"/strange/Br/br_contrib_OS_"+MODE+"_sm_func_mode_"+to_string(sm_func_mode)+"_"+ls_data_tm_VKVK.Tag[iens]+".dat", "", "#sigma Br_A0 Br_V0 Br_Aii Br_Vii");

    Print_To_File({}, {sigma_list_strange, Br_T_tau_tm[iens].ave(), Br_T_tau_tm[iens].err(), syst_per_ens_tm_T[iens], Br_L_tau_tm[iens].ave(), Br_L_tau_tm[iens].err(), syst_per_ens_tm_L[iens]}, "../data/tau_decay/"+Tag_reco_type+"/strange/Br/br_TL_contrib_tm_"+MODE+"_sm_func_mode_"+to_string(sm_func_mode)+"_"+ls_data_tm_VKVK.Tag[iens]+".dat", "", "#sigma Br_T Br_L");

    Print_To_File({}, {sigma_list_strange, Br_T_tau_OS[iens].ave(), Br_T_tau_OS[iens].err(), syst_per_ens_OS_T[iens], Br_L_tau_OS[iens].ave(), Br_L_tau_OS[iens].err(), syst_per_ens_OS_L[iens]}, "../data/tau_decay/"+Tag_reco_type+"/strange/Br/br_TL_contrib_OS_"+MODE+"_sm_func_mode_"+to_string(sm_func_mode)+"_"+ls_data_tm_VKVK.Tag[iens]+".dat", "", "#sigma Br_T Br_L");


    //pull

    Print_To_File( {}, {sigma_list_strange, pull_per_ens_tm_A0[iens], pull_per_ens_OS_A0[iens] }, "../data/tau_decay/"+Tag_reco_type+"/strange/pull/pull_"+ls_data_tm_VKVK.Tag[iens]+".A0", "", "#sigma tm OS" );
    Print_To_File( {}, {sigma_list_strange, pull_per_ens_tm_V0[iens], pull_per_ens_OS_V0[iens] }, "../data/tau_decay/"+Tag_reco_type+"/strange/pull/pull_"+ls_data_tm_VKVK.Tag[iens]+".V0", "", "#sigma tm OS" );

    Print_To_File( {}, {sigma_list_strange, pull_per_ens_tm_Aii[iens], pull_per_ens_OS_Aii[iens] }, "../data/tau_decay/"+Tag_reco_type+"/strange/pull/pull_"+ls_data_tm_VKVK.Tag[iens]+".Aii", "", "#sigma tm OS" );
    Print_To_File( {}, {sigma_list_strange, pull_per_ens_tm_Vii[iens], pull_per_ens_OS_Vii[iens] }, "../data/tau_decay/"+Tag_reco_type+"/strange/pull/pull_"+ls_data_tm_VKVK.Tag[iens]+".Vii", "", "#sigma tm OS" );
    
    Print_To_File( {}, {sigma_list_strange, pull_per_ens_tm_T[iens], pull_per_ens_OS_T[iens] }, "../data/tau_decay/"+Tag_reco_type+"/strange/pull/pull_"+ls_data_tm_VKVK.Tag[iens]+".T", "", "#sigma tm OS" );
    Print_To_File( {}, {sigma_list_strange, pull_per_ens_tm_L[iens], pull_per_ens_OS_L[iens] }, "../data/tau_decay/"+Tag_reco_type+"/strange/pull/pull_"+ls_data_tm_VKVK.Tag[iens]+".L", "", "#sigma tm OS" );
   
    
  }

  cout<<"Output per ensemble printed!"<<endl;

  //Print all ens for each sigma
  for(int is=0; is<(signed)sigma_list_strange.size();is++) {

    Vfloat syst_tm, syst_OS;
    distr_t_list Br_tau_tm_fixed_sigma(UseJack), Br_tau_OS_fixed_sigma(UseJack);
    for(int iens=0;iens<Nens;iens++) {
      Br_tau_tm_fixed_sigma.distr_list.push_back( Br_tau_tm[iens].distr_list[is]);
      Br_tau_OS_fixed_sigma.distr_list.push_back( Br_tau_OS[iens].distr_list[is]);
      syst_tm.push_back( syst_per_ens_tm[iens][is]);
      syst_OS.push_back( syst_per_ens_OS[iens][is]);
    }

    Print_To_File(ls_data_tm_VKVK.Tag,{ Br_tau_tm_fixed_sigma.ave(), Br_tau_tm_fixed_sigma.err(), syst_tm, Br_tau_OS_fixed_sigma.ave(), Br_tau_OS_fixed_sigma.err(), syst_OS},"../data/tau_decay/"+Tag_reco_type+"/strange/Br/br_sigma_"+to_string_with_precision(sigma_list_strange[is],3)+"_sm_func_mode_"+to_string(sm_func_mode)+".dat", "", "#Ens tm  OS,  Rs(HFLAV): "+to_string_with_precision(Rs_HFLAV,4)+" +- "+to_string_with_precision(D_Rs_HFLAV,4));

    
  }

  cout<<"output per sigma printed!"<<endl;

 



  //STORE JACKKNIFE DISTRIBUTIONS
  //light
  for(int i_ens=0; i_ens<Nens;i_ens++) {
    boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/jackknife/tm/A0A0/"+ls_data_tm_VKVK.Tag[i_ens]);
    boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/jackknife/tm/V0V0/"+ls_data_tm_VKVK.Tag[i_ens]);
    boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/jackknife/tm/AkAk/"+ls_data_tm_VKVK.Tag[i_ens]);
    boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/jackknife/tm/VkVk/"+ls_data_tm_VKVK.Tag[i_ens]);

    boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/jackknife/tm/T/"+ls_data_tm_VKVK.Tag[i_ens]);
    boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/jackknife/tm/L/"+ls_data_tm_VKVK.Tag[i_ens]);

    
    boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/jackknife/OS/A0A0/"+ls_data_tm_VKVK.Tag[i_ens]);
    boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/jackknife/OS/V0V0/"+ls_data_tm_VKVK.Tag[i_ens]);
    boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/jackknife/OS/AkAk/"+ls_data_tm_VKVK.Tag[i_ens]);
    boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/jackknife/OS/VkVk/"+ls_data_tm_VKVK.Tag[i_ens]);

    boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/jackknife/OS/T/"+ls_data_tm_VKVK.Tag[i_ens]);
    boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/jackknife/OS/L/"+ls_data_tm_VKVK.Tag[i_ens]);

    boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/jackknife/tm/A0A0/"+ls_data_tm_VKVK.Tag[i_ens]+"/sm_func_mode_"+to_string(sm_func_mode));
    boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/jackknife/tm/V0V0/"+ls_data_tm_VKVK.Tag[i_ens]+"/sm_func_mode_"+to_string(sm_func_mode));
    boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/jackknife/tm/AkAk/"+ls_data_tm_VKVK.Tag[i_ens]+"/sm_func_mode_"+to_string(sm_func_mode));
    boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/jackknife/tm/VkVk/"+ls_data_tm_VKVK.Tag[i_ens]+"/sm_func_mode_"+to_string(sm_func_mode));
    boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/jackknife/OS/A0A0/"+ls_data_tm_VKVK.Tag[i_ens]+"/sm_func_mode_"+to_string(sm_func_mode));
    boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/jackknife/OS/V0V0/"+ls_data_tm_VKVK.Tag[i_ens]+"/sm_func_mode_"+to_string(sm_func_mode));
    boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/jackknife/OS/AkAk/"+ls_data_tm_VKVK.Tag[i_ens]+"/sm_func_mode_"+to_string(sm_func_mode));
    boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/jackknife/OS/VkVk/"+ls_data_tm_VKVK.Tag[i_ens]+"/sm_func_mode_"+to_string(sm_func_mode));


    boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/jackknife/tm/T/"+ls_data_tm_VKVK.Tag[i_ens]+"/sm_func_mode_"+to_string(sm_func_mode));
    boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/jackknife/tm/L/"+ls_data_tm_VKVK.Tag[i_ens]+"/sm_func_mode_"+to_string(sm_func_mode));

    boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/jackknife/OS/T/"+ls_data_tm_VKVK.Tag[i_ens]+"/sm_func_mode_"+to_string(sm_func_mode));
    boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/jackknife/OS/L/"+ls_data_tm_VKVK.Tag[i_ens]+"/sm_func_mode_"+to_string(sm_func_mode));

  
    for(int is=0; is<(signed)sigma_list_strange.size();is++) {
      //print jackknife distribution for tm and OS
      ofstream Print_tm_A0A0("../data/tau_decay/"+Tag_reco_type+"/strange/jackknife/tm/A0A0/"+ls_data_tm_VKVK.Tag[i_ens]+"/sm_func_mode_"+to_string(sm_func_mode)+"/sigma_"+to_string_with_precision(sigma_list_strange[is],3)+".jack");
      ofstream Print_tm_V0V0("../data/tau_decay/"+Tag_reco_type+"/strange/jackknife/tm/V0V0/"+ls_data_tm_VKVK.Tag[i_ens]+"/sm_func_mode_"+to_string(sm_func_mode)+"/sigma_"+to_string_with_precision(sigma_list_strange[is],3)+".jack");
      ofstream Print_tm_AkAk("../data/tau_decay/"+Tag_reco_type+"/strange/jackknife/tm/AkAk/"+ls_data_tm_VKVK.Tag[i_ens]+"/sm_func_mode_"+to_string(sm_func_mode)+"/sigma_"+to_string_with_precision(sigma_list_strange[is],3)+".jack");
      ofstream Print_tm_VkVk("../data/tau_decay/"+Tag_reco_type+"/strange/jackknife/tm/VkVk/"+ls_data_tm_VKVK.Tag[i_ens]+"/sm_func_mode_"+to_string(sm_func_mode)+"/sigma_"+to_string_with_precision(sigma_list_strange[is],3)+".jack");

      ofstream Print_tm_T("../data/tau_decay/"+Tag_reco_type+"/strange/jackknife/tm/T/"+ls_data_tm_VKVK.Tag[i_ens]+"/sm_func_mode_"+to_string(sm_func_mode)+"/sigma_"+to_string_with_precision(sigma_list_strange[is],3)+".jack");
      ofstream Print_tm_L("../data/tau_decay/"+Tag_reco_type+"/strange/jackknife/tm/L/"+ls_data_tm_VKVK.Tag[i_ens]+"/sm_func_mode_"+to_string(sm_func_mode)+"/sigma_"+to_string_with_precision(sigma_list_strange[is],3)+".jack");

      ofstream Print_OS_A0A0("../data/tau_decay/"+Tag_reco_type+"/strange/jackknife/OS/A0A0/"+ls_data_tm_VKVK.Tag[i_ens]+"/sm_func_mode_"+to_string(sm_func_mode)+"/sigma_"+to_string_with_precision(sigma_list_strange[is],3)+".jack");
      ofstream Print_OS_V0V0("../data/tau_decay/"+Tag_reco_type+"/strange/jackknife/OS/V0V0/"+ls_data_tm_VKVK.Tag[i_ens]+"/sm_func_mode_"+to_string(sm_func_mode)+"/sigma_"+to_string_with_precision(sigma_list_strange[is],3)+".jack");
      ofstream Print_OS_AkAk("../data/tau_decay/"+Tag_reco_type+"/strange/jackknife/OS/AkAk/"+ls_data_tm_VKVK.Tag[i_ens]+"/sm_func_mode_"+to_string(sm_func_mode)+"/sigma_"+to_string_with_precision(sigma_list_strange[is],3)+".jack");
      ofstream Print_OS_VkVk("../data/tau_decay/"+Tag_reco_type+"/strange/jackknife/OS/VkVk/"+ls_data_tm_VKVK.Tag[i_ens]+"/sm_func_mode_"+to_string(sm_func_mode)+"/sigma_"+to_string_with_precision(sigma_list_strange[is],3)+".jack");

      ofstream Print_OS_T("../data/tau_decay/"+Tag_reco_type+"/strange/jackknife/OS/T/"+ls_data_tm_VKVK.Tag[i_ens]+"/sm_func_mode_"+to_string(sm_func_mode)+"/sigma_"+to_string_with_precision(sigma_list_strange[is],3)+".jack");
      ofstream Print_OS_L("../data/tau_decay/"+Tag_reco_type+"/strange/jackknife/OS/L/"+ls_data_tm_VKVK.Tag[i_ens]+"/sm_func_mode_"+to_string(sm_func_mode)+"/sigma_"+to_string_with_precision(sigma_list_strange[is],3)+".jack");
     
      for(int ijack=0; ijack<Njacks;ijack++) {
	Print_tm_A0A0<<Br_A0_tau_tm[i_ens].distr_list[is].distr[ijack]<<endl;
	Print_tm_V0V0<<Br_V0_tau_tm[i_ens].distr_list[is].distr[ijack]<<endl;
	Print_tm_AkAk<<Br_Aii_tau_tm[i_ens].distr_list[is].distr[ijack]<<endl;
	Print_tm_VkVk<<Br_Vii_tau_tm[i_ens].distr_list[is].distr[ijack]<<endl;

	Print_tm_T<<Br_T_tau_tm[i_ens].distr_list[is].distr[ijack]<<endl;
	Print_tm_L<<Br_L_tau_tm[i_ens].distr_list[is].distr[ijack]<<endl;
	
	Print_OS_A0A0<<Br_A0_tau_OS[i_ens].distr_list[is].distr[ijack]<<endl;
	Print_OS_V0V0<<Br_V0_tau_OS[i_ens].distr_list[is].distr[ijack]<<endl;
	Print_OS_AkAk<<Br_Aii_tau_OS[i_ens].distr_list[is].distr[ijack]<<endl;
	Print_OS_VkVk<<Br_Vii_tau_OS[i_ens].distr_list[is].distr[ijack]<<endl;

	Print_OS_T<<Br_T_tau_OS[i_ens].distr_list[is].distr[ijack]<<endl;
	Print_OS_L<<Br_L_tau_OS[i_ens].distr_list[is].distr[ijack]<<endl;

	
      }
      Print_tm_A0A0.close();
      Print_tm_V0V0.close();
      Print_tm_AkAk.close();
      Print_tm_VkVk.close();

      Print_tm_T.close();
      Print_tm_L.close();
      
      Print_OS_A0A0.close();
      Print_OS_V0V0.close();
      Print_OS_AkAk.close();
      Print_OS_VkVk.close();

      Print_OS_T.close();
      Print_OS_L.close();

      
    }

   
  }
  cout<<"Jackknives printed!"<<endl;
  }


  



  if(Perform_continuum_extrapolation) {

    cout<<"Performing continuum limit extrapolation"<<endl;


    
    //vector<string> Contribs({"A0A0", "V0V0", "AkAk", "VkVk", "tot", "VA", "VMA", "AX", "T", "L", "tot_TL"});
    vector<string> Contribs({"T", "L", "tot_TL"});
    //vector<string> Fit_types({"tm", "OS", "comb"});
    vector<string> Fit_types({"comb"});
    //vector<string> poly_types({"const", "linear"});
    vector<string> poly_types({"const", "linear", "tm_linear", "OS_linear", "quad", "tm_quad", "OS_quad"});
    
    distr_t_list test(UseJack);
    map< tuple<string,string,string> , distr_t_list> res_map;
    map< tuple<string, string,string>, vector<double>> ch2_map;
    map< tuple<string, string,string>, vector<int>> Ndof_map;
    map< tuple<string, string,string>, vector<int>> Nmeas_map;
     
     
    for(auto &contr: Contribs) {
      for(auto &ftype: Fit_types) {
	for(auto &ptype: poly_types) {

	  tuple<string,string,string> Keey= {contr,ftype,ptype};
	  res_map.emplace(Keey, distr_t_list(UseJack) );
	  ch2_map.emplace(Keey,vector<double>{});
	  Ndof_map.emplace(Keey, vector<int>{});
	  Nmeas_map.emplace(Keey, vector<int>{});
	}
      }
    }

    
    
    int Nlat=300;
    Vfloat a_to_print;
    double sx= 0.08*1.5/(Nlat-1.0); //fm
    for(int pp=0;pp<Nlat;pp++) a_to_print.push_back( sx*pp);


    //order of the limit:    first L -> infty,   then a -> 0, then  sigma -> 0

    //outer loop is over sigma

    vector<distr_t_list> A0A0_tm_all_s(sigma_list_strange.size());
    vector<distr_t_list> V0V0_tm_all_s(sigma_list_strange.size());
    vector<distr_t_list> AkAk_tm_all_s(sigma_list_strange.size());
    vector<distr_t_list> VkVk_tm_all_s(sigma_list_strange.size());
    vector<distr_t_list> T_tm_all_s(sigma_list_strange.size());
    vector<distr_t_list> L_tm_all_s(sigma_list_strange.size());

    vector<distr_t_list> A0A0_OS_all_s(sigma_list_strange.size());
    vector<distr_t_list> V0V0_OS_all_s(sigma_list_strange.size());
    vector<distr_t_list> AkAk_OS_all_s(sigma_list_strange.size());
    vector<distr_t_list> VkVk_OS_all_s(sigma_list_strange.size());
    vector<distr_t_list> T_OS_all_s(sigma_list_strange.size());
    vector<distr_t_list> L_OS_all_s(sigma_list_strange.size());

    distr_t_list tot_TL_Ens_E_tm(UseJack), tot_TL_Ens_E_OS(UseJack);

    //load systematic errors

     
    VVfloat syst_A0A0_tm(Nens), syst_V0V0_tm(Nens),  syst_AkAk_tm(Nens), syst_VkVk_tm(Nens), syst_T_tm(Nens), syst_L_tm(Nens);
    VVfloat syst_A0A0_OS(Nens), syst_V0V0_OS(Nens),  syst_AkAk_OS(Nens), syst_VkVk_OS(Nens), syst_T_OS(Nens), syst_L_OS(Nens);


    VVfloat stat_A0A0_tm(Nens) , stat_V0V0_tm(Nens), stat_AkAk_tm(Nens), stat_VkVk_tm(Nens), stat_T_tm(Nens), stat_L_tm(Nens);
    VVfloat stat_A0A0_OS(Nens) , stat_V0V0_OS(Nens), stat_AkAk_OS(Nens), stat_VkVk_OS(Nens), stat_T_OS(Nens), stat_L_OS(Nens);

    
    for(int iens=0;iens<Nens;iens++) {
      //A0A0
      syst_A0A0_tm[iens] = Read_From_File("../data/tau_decay/"+Tag_reco_type+"/strange/Br/br_contrib_tm_"+MODE+"_sm_func_mode_"+to_string(sm_func_mode)+"_"+ls_data_tm_VKVK.Tag[iens]+".dat", 4,  14,1);
      syst_A0A0_OS[iens] = Read_From_File("../data/tau_decay/"+Tag_reco_type+"/strange/Br/br_contrib_OS_"+MODE+"_sm_func_mode_"+to_string(sm_func_mode)+"_"+ls_data_tm_VKVK.Tag[iens]+".dat", 4,  14,1);

      //V0V0
      syst_V0V0_tm[iens] = Read_From_File("../data/tau_decay/"+Tag_reco_type+"/strange/Br/br_contrib_tm_"+MODE+"_sm_func_mode_"+to_string(sm_func_mode)+"_"+ls_data_tm_VKVK.Tag[iens]+".dat", 7,  14,1);
      syst_V0V0_OS[iens] = Read_From_File("../data/tau_decay/"+Tag_reco_type+"/strange/Br/br_contrib_OS_"+MODE+"_sm_func_mode_"+to_string(sm_func_mode)+"_"+ls_data_tm_VKVK.Tag[iens]+".dat", 7,  14,1);
      
      //AkAk
      syst_AkAk_tm[iens] = Read_From_File("../data/tau_decay/"+Tag_reco_type+"/strange/Br/br_contrib_tm_"+MODE+"_sm_func_mode_"+to_string(sm_func_mode)+"_"+ls_data_tm_VKVK.Tag[iens]+".dat", 10,  14,1);
      syst_AkAk_OS[iens] = Read_From_File("../data/tau_decay/"+Tag_reco_type+"/strange/Br/br_contrib_OS_"+MODE+"_sm_func_mode_"+to_string(sm_func_mode)+"_"+ls_data_tm_VKVK.Tag[iens]+".dat", 10,  14,1);
      //VkVk
      syst_VkVk_tm[iens] = Read_From_File("../data/tau_decay/"+Tag_reco_type+"/strange/Br/br_contrib_tm_"+MODE+"_sm_func_mode_"+to_string(sm_func_mode)+"_"+ls_data_tm_VKVK.Tag[iens]+".dat", 13,  14,1);
      syst_VkVk_OS[iens] = Read_From_File("../data/tau_decay/"+Tag_reco_type+"/strange/Br/br_contrib_OS_"+MODE+"_sm_func_mode_"+to_string(sm_func_mode)+"_"+ls_data_tm_VKVK.Tag[iens]+".dat", 13,  14,1);
      //T
      syst_T_tm[iens] = Read_From_File("../data/tau_decay/"+Tag_reco_type+"/strange/Br/br_TL_contrib_tm_"+MODE+"_sm_func_mode_"+to_string(sm_func_mode)+"_"+ls_data_tm_VKVK.Tag[iens]+".dat", 4,  8,1);
      syst_T_OS[iens] = Read_From_File("../data/tau_decay/"+Tag_reco_type+"/strange/Br/br_TL_contrib_OS_"+MODE+"_sm_func_mode_"+to_string(sm_func_mode)+"_"+ls_data_tm_VKVK.Tag[iens]+".dat", 4,  8,1);

      //L
      syst_L_tm[iens] = Read_From_File("../data/tau_decay/"+Tag_reco_type+"/strange/Br/br_TL_contrib_tm_"+MODE+"_sm_func_mode_"+to_string(sm_func_mode)+"_"+ls_data_tm_VKVK.Tag[iens]+".dat", 7,  8,1);
      syst_L_OS[iens] = Read_From_File("../data/tau_decay/"+Tag_reco_type+"/strange/Br/br_TL_contrib_OS_"+MODE+"_sm_func_mode_"+to_string(sm_func_mode)+"_"+ls_data_tm_VKVK.Tag[iens]+".dat", 7,  8,1);
      
    }


    //to print FSEs corrections and syst

    Vfloat FSE_C_A0A0_tm, FSE_C_V0V0_tm, FSE_C_AkAk_tm, FSE_C_VkVk_tm, FSE_C_T_tm, FSE_C_L_tm;
    Vfloat FSE_C_A0A0_OS, FSE_C_V0V0_OS, FSE_C_AkAk_OS, FSE_C_VkVk_OS, FSE_C_T_OS, FSE_C_L_OS;

    Vfloat FSE_C_A0A0, FSE_C_V0V0, FSE_C_AkAk, FSE_C_VkVk, FSE_C_T, FSE_C_L;


    Vfloat pl_T_tm_distr, pl_T_OS_distr, pl_L_tm_distr, pl_L_OS_distr;
    Vfloat pl_BIS_T_tm_distr, pl_BIS_T_OS_distr, pl_BIS_L_tm_distr, pl_BIS_L_OS_distr;

   
  
    Vfloat FSE_B_A0A0_tm, FSE_B_V0V0_tm, FSE_B_AkAk_tm, FSE_B_VkVk_tm, FSE_B_T_tm, FSE_B_L_tm;
    Vfloat FSE_B_A0A0_OS, FSE_B_V0V0_OS, FSE_B_AkAk_OS, FSE_B_VkVk_OS, FSE_B_T_OS, FSE_B_L_OS;


    Vfloat FSE_err_C_A0A0_tm, FSE_err_C_V0V0_tm, FSE_err_C_AkAk_tm, FSE_err_C_VkVk_tm, FSE_err_C_T_tm, FSE_err_C_L_tm;
    Vfloat FSE_err_C_A0A0_OS, FSE_err_C_V0V0_OS, FSE_err_C_AkAk_OS, FSE_err_C_VkVk_OS, FSE_err_C_T_OS, FSE_err_C_L_OS;

    Vfloat FSE_err_B_A0A0_tm, FSE_err_B_V0V0_tm, FSE_err_B_AkAk_tm, FSE_err_B_VkVk_tm, FSE_err_B_T_tm, FSE_err_B_L_tm;
    Vfloat FSE_err_B_A0A0_OS, FSE_err_B_V0V0_OS, FSE_err_B_AkAk_OS, FSE_err_B_VkVk_OS, FSE_err_B_T_OS, FSE_err_B_L_OS;

     

    for(int is=0; is < (signed)sigma_list_strange.size() ; is++) { //loop over sigma

      //load data
      distr_t_list  A0A0_tm(UseJack), V0V0_tm(UseJack),  AkAk_tm(UseJack), VkVk_tm(UseJack), T_tm(UseJack), L_tm(UseJack);
      distr_t_list  A0A0_OS(UseJack), V0V0_OS(UseJack),  AkAk_OS(UseJack), VkVk_OS(UseJack), T_OS(UseJack), L_OS(UseJack);


      GaussianMersenne GM_sigma(7254324);
            

      for(int iens=0;iens<Nens;iens++) {


	//load A0A0 jack distribution
	A0A0_tm.distr_list.emplace_back(UseJack, Read_From_File("../data/tau_decay/"+Tag_reco_type+"/strange/jackknife/tm/A0A0/"+ls_data_tm_VKVK.Tag[iens]+"/sm_func_mode_"+to_string(sm_func_mode)+"/sigma_"+to_string_with_precision(sigma_list_strange[is],3)+".jack", 0,1));
	A0A0_OS.distr_list.emplace_back(UseJack, Read_From_File("../data/tau_decay/"+Tag_reco_type+"/strange/jackknife/OS/A0A0/"+ls_data_tm_VKVK.Tag[iens]+"/sm_func_mode_"+to_string(sm_func_mode)+"/sigma_"+to_string_with_precision(sigma_list_strange[is],3)+".jack", 0,1));


	//load V0V0 jack distribution
	V0V0_tm.distr_list.emplace_back(UseJack, Read_From_File("../data/tau_decay/"+Tag_reco_type+"/strange/jackknife/tm/V0V0/"+ls_data_tm_VKVK.Tag[iens]+"/sm_func_mode_"+to_string(sm_func_mode)+"/sigma_"+to_string_with_precision(sigma_list_strange[is],3)+".jack", 0,1));
	V0V0_OS.distr_list.emplace_back(UseJack, Read_From_File("../data/tau_decay/"+Tag_reco_type+"/strange/jackknife/OS/V0V0/"+ls_data_tm_VKVK.Tag[iens]+"/sm_func_mode_"+to_string(sm_func_mode)+"/sigma_"+to_string_with_precision(sigma_list_strange[is],3)+".jack", 0,1));

	
	//load AkAk jack distribution
	AkAk_tm.distr_list.emplace_back(UseJack, Read_From_File("../data/tau_decay/"+Tag_reco_type+"/strange/jackknife/tm/AkAk/"+ls_data_tm_VKVK.Tag[iens]+"/sm_func_mode_"+to_string(sm_func_mode)+"/sigma_"+to_string_with_precision(sigma_list_strange[is],3)+".jack", 0,1));
	AkAk_OS.distr_list.emplace_back(UseJack, Read_From_File("../data/tau_decay/"+Tag_reco_type+"/strange/jackknife/OS/AkAk/"+ls_data_tm_VKVK.Tag[iens]+"/sm_func_mode_"+to_string(sm_func_mode)+"/sigma_"+to_string_with_precision(sigma_list_strange[is],3)+".jack", 0,1));

	
	//load VkVk jack distribution
	VkVk_tm.distr_list.emplace_back(UseJack, Read_From_File("../data/tau_decay/"+Tag_reco_type+"/strange/jackknife/tm/VkVk/"+ls_data_tm_VKVK.Tag[iens]+"/sm_func_mode_"+to_string(sm_func_mode)+"/sigma_"+to_string_with_precision(sigma_list_strange[is],3)+".jack", 0,1));
	VkVk_OS.distr_list.emplace_back(UseJack, Read_From_File("../data/tau_decay/"+Tag_reco_type+"/strange/jackknife/OS/VkVk/"+ls_data_tm_VKVK.Tag[iens]+"/sm_func_mode_"+to_string(sm_func_mode)+"/sigma_"+to_string_with_precision(sigma_list_strange[is],3)+".jack", 0,1));

	//load T jack distribution
	T_tm.distr_list.emplace_back(UseJack, Read_From_File("../data/tau_decay/"+Tag_reco_type+"/strange/jackknife/tm/T/"+ls_data_tm_VKVK.Tag[iens]+"/sm_func_mode_"+to_string(sm_func_mode)+"/sigma_"+to_string_with_precision(sigma_list_strange[is],3)+".jack", 0,1));
	T_OS.distr_list.emplace_back(UseJack, Read_From_File("../data/tau_decay/"+Tag_reco_type+"/strange/jackknife/OS/T/"+ls_data_tm_VKVK.Tag[iens]+"/sm_func_mode_"+to_string(sm_func_mode)+"/sigma_"+to_string_with_precision(sigma_list_strange[is],3)+".jack", 0,1));

	//load L jack distribution
	L_tm.distr_list.emplace_back(UseJack, Read_From_File("../data/tau_decay/"+Tag_reco_type+"/strange/jackknife/tm/L/"+ls_data_tm_VKVK.Tag[iens]+"/sm_func_mode_"+to_string(sm_func_mode)+"/sigma_"+to_string_with_precision(sigma_list_strange[is],3)+".jack", 0,1));
	L_OS.distr_list.emplace_back(UseJack, Read_From_File("../data/tau_decay/"+Tag_reco_type+"/strange/jackknife/OS/L/"+ls_data_tm_VKVK.Tag[iens]+"/sm_func_mode_"+to_string(sm_func_mode)+"/sigma_"+to_string_with_precision(sigma_list_strange[is],3)+".jack", 0,1));


	//push-back statistical errors

	stat_A0A0_tm[iens].push_back( A0A0_tm.err(iens));
	stat_A0A0_OS[iens].push_back( A0A0_OS.err(iens));

	stat_V0V0_tm[iens].push_back( V0V0_tm.err(iens));
	stat_V0V0_OS[iens].push_back( V0V0_OS.err(iens));

	stat_AkAk_tm[iens].push_back( AkAk_tm.err(iens));
	stat_AkAk_OS[iens].push_back( AkAk_OS.err(iens));

	stat_VkVk_tm[iens].push_back( VkVk_tm.err(iens));
	stat_VkVk_OS[iens].push_back( VkVk_OS.err(iens));

	stat_T_tm[iens].push_back( T_tm.err(iens));
	stat_T_OS[iens].push_back( T_OS.err(iens));

	stat_L_tm[iens].push_back( L_tm.err(iens));
	stat_L_OS[iens].push_back( L_OS.err(iens));

	//add systematic error

	//A0A0
	A0A0_tm.distr_list[iens] = A0A0_tm.ave(iens) + (A0A0_tm.distr_list[iens]- A0A0_tm.ave(iens))*( sqrt( pow(A0A0_tm.err(iens),2) + pow(syst_A0A0_tm[iens][is],2)))/A0A0_tm.err(iens);
	A0A0_OS.distr_list[iens] = A0A0_OS.ave(iens) + (A0A0_OS.distr_list[iens]- A0A0_OS.ave(iens))*( sqrt( pow(A0A0_OS.err(iens),2) + pow(syst_A0A0_OS[iens][is],2)))/A0A0_OS.err(iens);
	//V0V0
	V0V0_tm.distr_list[iens] = V0V0_tm.ave(iens) + (V0V0_tm.distr_list[iens]- V0V0_tm.ave(iens))*( sqrt( pow(V0V0_tm.err(iens),2) + pow(syst_V0V0_tm[iens][is],2)))/V0V0_tm.err(iens);
	V0V0_OS.distr_list[iens] = V0V0_OS.ave(iens) + (V0V0_OS.distr_list[iens]- V0V0_OS.ave(iens))*( sqrt( pow(V0V0_OS.err(iens),2) + pow(syst_V0V0_OS[iens][is],2)))/V0V0_OS.err(iens);
	
	//AkAk
	AkAk_tm.distr_list[iens] = AkAk_tm.ave(iens) + (AkAk_tm.distr_list[iens]- AkAk_tm.ave(iens))*( sqrt( pow(AkAk_tm.err(iens),2) + pow(syst_AkAk_tm[iens][is],2)))/AkAk_tm.err(iens);
	AkAk_OS.distr_list[iens] = AkAk_OS.ave(iens) + (AkAk_OS.distr_list[iens]- AkAk_OS.ave(iens))*( sqrt( pow(AkAk_OS.err(iens),2) + pow(syst_AkAk_OS[iens][is],2)))/AkAk_OS.err(iens);
	//VkVk
	VkVk_tm.distr_list[iens] = VkVk_tm.ave(iens) + (VkVk_tm.distr_list[iens]- VkVk_tm.ave(iens))*( sqrt( pow(VkVk_tm.err(iens),2) + pow(syst_VkVk_tm[iens][is],2)))/VkVk_tm.err(iens);
	VkVk_OS.distr_list[iens] = VkVk_OS.ave(iens) + (VkVk_OS.distr_list[iens]- VkVk_OS.ave(iens))*( sqrt( pow(VkVk_OS.err(iens),2) + pow(syst_VkVk_OS[iens][is],2)))/VkVk_OS.err(iens);

	//T
	T_tm.distr_list[iens] = T_tm.ave(iens) + (T_tm.distr_list[iens]- T_tm.ave(iens))*( sqrt( pow(T_tm.err(iens),2) + pow(syst_T_tm[iens][is],2)))/T_tm.err(iens);
	T_OS.distr_list[iens] = T_OS.ave(iens) + (T_OS.distr_list[iens]- T_OS.ave(iens))*( sqrt( pow(T_OS.err(iens),2) + pow(syst_T_OS[iens][is],2)))/T_OS.err(iens);

	//L
	L_tm.distr_list[iens] = L_tm.ave(iens) + (L_tm.distr_list[iens]- L_tm.ave(iens))*( sqrt( pow(L_tm.err(iens),2) + pow(syst_L_tm[iens][is],2)))/L_tm.err(iens);
	L_OS.distr_list[iens] = L_OS.ave(iens) + (L_OS.distr_list[iens]- L_OS.ave(iens))*( sqrt( pow(L_OS.err(iens),2) + pow(syst_L_OS[iens][is],2)))/L_OS.err(iens);
	
      }

   
      //push back to all_s
      A0A0_tm_all_s[is] = A0A0_tm;
      A0A0_OS_all_s[is] = A0A0_OS;
      V0V0_tm_all_s[is] = V0V0_tm;
      V0V0_OS_all_s[is] = V0V0_OS;
      AkAk_tm_all_s[is] = AkAk_tm;
      AkAk_OS_all_s[is] = AkAk_OS;
      VkVk_tm_all_s[is] = VkVk_tm;
      VkVk_OS_all_s[is] = VkVk_OS;
      T_tm_all_s[is] = T_tm;
      T_OS_all_s[is] = T_OS;
      L_tm_all_s[is] = L_tm;
      L_OS_all_s[is] = L_OS;

      distr_t_list a_distr_list(UseJack);


      //estimate FSEs at fixed sigma from difference between B96 and B64

      distr_t_list V0V0_tm_red(UseJack), V0V0_OS_red(UseJack);
      distr_t_list A0A0_tm_red(UseJack), A0A0_OS_red(UseJack);
      distr_t_list AkAk_tm_red(UseJack), AkAk_OS_red(UseJack);
      distr_t_list VkVk_tm_red(UseJack), VkVk_OS_red(UseJack);
      distr_t_list T_tm_red(UseJack), T_OS_red(UseJack);
      distr_t_list L_tm_red(UseJack), L_OS_red(UseJack);

      distr_t_list a_distr_list_red(UseJack);
      vector<string> Tag_ens_red;

      
       
      //###################################################
      //############## COMPUTE FSE ########################
      double corr_fact_A0A0_FSE, corr_fact_A0A0_FSE_tm, corr_fact_A0A0_FSE_OS;
      double corr_fact_V0V0_FSE, corr_fact_V0V0_FSE_tm, corr_fact_V0V0_FSE_OS;
      double corr_fact_AkAk_FSE, corr_fact_AkAk_FSE_tm, corr_fact_AkAk_FSE_OS;
      double corr_fact_VkVk_FSE, corr_fact_VkVk_FSE_tm, corr_fact_VkVk_FSE_OS;
      double corr_fact_T_FSE, corr_fact_T_FSE_tm, corr_fact_T_FSE_OS;
      double corr_fact_L_FSE, corr_fact_L_FSE_tm, corr_fact_L_FSE_OS;

      double corr_fact_err_A0A0_FSE, corr_fact_err_A0A0_FSE_tm, corr_fact_err_A0A0_FSE_OS;
      double corr_fact_err_V0V0_FSE, corr_fact_err_V0V0_FSE_tm, corr_fact_err_V0V0_FSE_OS;
      double corr_fact_err_AkAk_FSE, corr_fact_err_AkAk_FSE_tm, corr_fact_err_AkAk_FSE_OS;
      double corr_fact_err_VkVk_FSE, corr_fact_err_VkVk_FSE_tm, corr_fact_err_VkVk_FSE_OS;
      double corr_fact_err_T_FSE, corr_fact_err_T_FSE_tm, corr_fact_err_T_FSE_OS;
      double corr_fact_err_L_FSE, corr_fact_err_L_FSE_tm, corr_fact_err_L_FSE_OS;

      double pl_T_tm, pl_T_OS, pl_L_tm, pl_L_OS;
      
      int id_C80=-1;
      int id_C112=-1;
      //############# COMPUTE FSE FROM B64-B96 ###########
      double corr_fact_BIS_A0A0_FSE, corr_fact_BIS_A0A0_FSE_tm, corr_fact_BIS_A0A0_FSE_OS;
      double corr_fact_BIS_V0V0_FSE, corr_fact_BIS_V0V0_FSE_tm, corr_fact_BIS_V0V0_FSE_OS;
      double corr_fact_BIS_AkAk_FSE, corr_fact_BIS_AkAk_FSE_tm, corr_fact_BIS_AkAk_FSE_OS;
      double corr_fact_BIS_VkVk_FSE, corr_fact_BIS_VkVk_FSE_tm, corr_fact_BIS_VkVk_FSE_OS;
      double corr_fact_BIS_T_FSE, corr_fact_BIS_T_FSE_tm, corr_fact_BIS_T_FSE_OS;
      double corr_fact_BIS_L_FSE, corr_fact_BIS_L_FSE_tm, corr_fact_BIS_L_FSE_OS;

      double corr_fact_BIS_err_A0A0_FSE, corr_fact_BIS_err_A0A0_FSE_tm, corr_fact_BIS_err_A0A0_FSE_OS;
      double corr_fact_BIS_err_V0V0_FSE, corr_fact_BIS_err_V0V0_FSE_tm, corr_fact_BIS_err_V0V0_FSE_OS;
      double corr_fact_BIS_err_AkAk_FSE, corr_fact_BIS_err_AkAk_FSE_tm, corr_fact_BIS_err_AkAk_FSE_OS;
      double corr_fact_BIS_err_VkVk_FSE, corr_fact_BIS_err_VkVk_FSE_tm, corr_fact_BIS_err_VkVk_FSE_OS;
      double corr_fact_BIS_err_T_FSE, corr_fact_BIS_err_T_FSE_tm, corr_fact_BIS_err_T_FSE_OS;
      double corr_fact_BIS_err_L_FSE, corr_fact_BIS_err_L_FSE_tm, corr_fact_BIS_err_L_FSE_OS;

      double pl_BIS_T_tm, pl_BIS_T_OS, pl_BIS_L_tm, pl_BIS_L_OS;
      
      int id_B64=-1;
      int id_B96=-1;
      //##################################################

      for(int iens=0;iens<Nens;iens++) {
	if(ls_data_tm_VKVK.Tag[iens] == "cC211a.06.80") { id_C80=iens;}
	else if(ls_data_tm_VKVK.Tag[iens] == "cC211a.06.112") { id_C112=iens;}
      }

      if( (id_C80==-1) || (id_C112==-1)) crash("Cannot find id_C80 and/or id_C112");

      if(id_C80 == id_C112) crash("Error: id_C80 == id_C112");


     for(int iens=0;iens<Nens;iens++) {
	if(ls_data_tm_VKVK.Tag[iens] == "cB211b.072.64") { id_B64=iens;}
	else if(ls_data_tm_VKVK.Tag[iens] == "cB211b.072.96") { id_B96=iens;}
      }

      if( (id_B64==-1) || (id_B96==-1)) crash("Cannot find id_B64 and/or id_B96");

      if(id_B64 == id_B96) crash("Error: id_B64 == id_B96");
      

      //A0A0
      corr_fact_A0A0_FSE_tm = ((A0A0_tm.ave(id_C112) - A0A0_tm.ave(id_C80))/(A0A0_tm.ave(id_C80)))*fabs(erf( (A0A0_tm.ave(id_C112)-A0A0_tm.ave(id_C80))/(sqrt(2.0*( pow( A0A0_tm.err(id_C112)  ,2)  + pow( A0A0_tm.err(id_C80) ,2)   )))));
      corr_fact_A0A0_FSE_OS = ((A0A0_OS.ave(id_C112) - A0A0_OS.ave(id_C80))/(A0A0_OS.ave(id_C80)))*fabs(erf( (A0A0_OS.ave(id_C112)-A0A0_OS.ave(id_C80))/(sqrt(2.0*( pow( A0A0_OS.err(id_C112)  ,2)  + pow( A0A0_OS.err(id_C80) ,2)   )))));

      corr_fact_err_A0A0_FSE_tm = fabs((( (A0A0_tm.distr_list[id_C112] - A0A0_tm.distr_list[id_C80]).err()  )/(A0A0_tm.ave(id_C80))))*fabs(erf( (A0A0_tm.ave(id_C112)-A0A0_tm.ave(id_C80))/(sqrt(2.0*( pow( A0A0_tm.err(id_C112)  ,2)  + pow( A0A0_tm.err(id_C80) ,2)   )))));
      corr_fact_err_A0A0_FSE_OS = fabs((( (A0A0_OS.distr_list[id_C112] - A0A0_OS.distr_list[id_C80]).err() )/(A0A0_OS.ave(id_C80))))*fabs(erf( (A0A0_OS.ave(id_C112)-A0A0_OS.ave(id_C80))/(sqrt(2.0*( pow( A0A0_OS.err(id_C112)  ,2)  + pow( A0A0_OS.err(id_C80) ,2)   )))));
      
      corr_fact_A0A0_FSE = max(fabs(corr_fact_A0A0_FSE_tm), fabs(corr_fact_A0A0_FSE_OS));
    
      //V0V0
      corr_fact_V0V0_FSE_tm = ((V0V0_tm.ave(id_C112) - V0V0_tm.ave(id_C80))/(V0V0_tm.ave(id_C80)))*fabs(erf( (V0V0_tm.ave(id_C112)-V0V0_tm.ave(id_C80))/(sqrt(2.0*( pow( V0V0_tm.err(id_C112)  ,2)  + pow( V0V0_tm.err(id_C80) ,2)   )))));
      corr_fact_V0V0_FSE_OS = ((V0V0_OS.ave(id_C112) - V0V0_OS.ave(id_C80))/(V0V0_OS.ave(id_C80)))*fabs(erf( (V0V0_OS.ave(id_C112)-V0V0_OS.ave(id_C80))/(sqrt(2.0*( pow( V0V0_OS.err(id_C112)  ,2)  + pow( V0V0_OS.err(id_C80) ,2)   )))));

      corr_fact_err_V0V0_FSE_tm = fabs((( (V0V0_tm.distr_list[id_C112] - V0V0_tm.distr_list[id_C80]).err()  )/(V0V0_tm.ave(id_C80))))*fabs(erf( (V0V0_tm.ave(id_C112)-V0V0_tm.ave(id_C80))/(sqrt(2.0*( pow( V0V0_tm.err(id_C112)  ,2)  + pow( V0V0_tm.err(id_C80) ,2)   )))));
      corr_fact_err_V0V0_FSE_OS = fabs(( (V0V0_OS.distr_list[id_C112] - V0V0_OS.distr_list[id_C80]).err() )/(V0V0_OS.ave(id_C80)))*fabs(erf( (V0V0_OS.ave(id_C112)-V0V0_OS.ave(id_C80))/(sqrt(2.0*( pow( V0V0_OS.err(id_C112)  ,2)  + pow( V0V0_OS.err(id_C80) ,2)   )))));
      
      corr_fact_V0V0_FSE = max( fabs(corr_fact_V0V0_FSE_tm), fabs(corr_fact_V0V0_FSE_OS));
           
      

      //AkAk
      corr_fact_AkAk_FSE_tm = ((AkAk_tm.ave(id_C112) - AkAk_tm.ave(id_C80))/(AkAk_tm.ave(id_C80)))*fabs(erf( (AkAk_tm.ave(id_C112)-AkAk_tm.ave(id_C80))/(sqrt(2.0*( pow( AkAk_tm.err(id_C112)  ,2)  + pow( AkAk_tm.err(id_C80) ,2)   )))));
      corr_fact_AkAk_FSE_OS = ((AkAk_OS.ave(id_C112) - AkAk_OS.ave(id_C80))/(AkAk_OS.ave(id_C80)))*fabs(erf( (AkAk_OS.ave(id_C112)-AkAk_OS.ave(id_C80))/(sqrt(2.0*( pow( AkAk_OS.err(id_C112)  ,2)  + pow( AkAk_OS.err(id_C80) ,2)   )))));

      corr_fact_err_AkAk_FSE_tm = fabs(( (AkAk_tm.distr_list[id_C112] - AkAk_tm.distr_list[id_C80]).err()  )/(AkAk_tm.ave(id_C80)))*fabs(erf( (AkAk_tm.ave(id_C112)-AkAk_tm.ave(id_C80))/(sqrt(2.0*( pow( AkAk_tm.err(id_C112)  ,2)  + pow( AkAk_tm.err(id_C80) ,2)   )))));
      corr_fact_err_AkAk_FSE_OS = fabs(( (AkAk_OS.distr_list[id_C112] - AkAk_OS.distr_list[id_C80]).err()  )/(AkAk_OS.ave(id_C80)))*fabs(erf( (AkAk_OS.ave(id_C112)-AkAk_OS.ave(id_C80))/(sqrt(2.0*( pow( AkAk_OS.err(id_C112)  ,2)  + pow( AkAk_OS.err(id_C80) ,2)   )))));
      
      corr_fact_AkAk_FSE = max(fabs(corr_fact_AkAk_FSE_tm), fabs(corr_fact_AkAk_FSE_OS));
   
      

      //VkVk
      corr_fact_VkVk_FSE_tm = ((VkVk_tm.ave(id_C112) - VkVk_tm.ave(id_C80))/(VkVk_tm.ave(id_C80)))*fabs(erf( (VkVk_tm.ave(id_C112)-VkVk_tm.ave(id_C80))/(sqrt(2.0*( pow( VkVk_tm.err(id_C112)  ,2)  + pow( VkVk_tm.err(id_C80) ,2)   )))));
      corr_fact_VkVk_FSE_OS = ((VkVk_OS.ave(id_C112) - VkVk_OS.ave(id_C80))/(VkVk_OS.ave(id_C80)))*fabs(erf( (VkVk_OS.ave(id_C112)-VkVk_OS.ave(id_C80))/(sqrt(2.0*( pow( VkVk_OS.err(id_C112)  ,2)  + pow( VkVk_OS.err(id_C80) ,2)   )))));

      corr_fact_err_VkVk_FSE_tm = fabs(( (VkVk_tm.distr_list[id_C112] - VkVk_tm.distr_list[id_C80]).err()  )/(VkVk_tm.ave(id_C80)))*fabs(erf( (VkVk_tm.ave(id_C112)-VkVk_tm.ave(id_C80))/(sqrt(2.0*( pow( VkVk_tm.err(id_C112)  ,2)  + pow( VkVk_tm.err(id_C80) ,2)   )))));
      corr_fact_err_VkVk_FSE_OS = fabs(( (VkVk_OS.distr_list[id_C112] - VkVk_OS.distr_list[id_C80]).err() )/(VkVk_OS.ave(id_C80)))*fabs(erf( (VkVk_OS.ave(id_C112)-VkVk_OS.ave(id_C80))/(sqrt(2.0*( pow( VkVk_OS.err(id_C112)  ,2)  + pow( VkVk_OS.err(id_C80) ,2)   )))));
      
      corr_fact_VkVk_FSE = max( fabs(corr_fact_VkVk_FSE_tm), fabs(corr_fact_VkVk_FSE_OS));

      //T
      corr_fact_T_FSE_tm = ((T_tm.ave(id_C112) - T_tm.ave(id_C80))/(T_tm.ave(id_C80)))*fabs(erf( (T_tm.ave(id_C112)-T_tm.ave(id_C80))/(sqrt(2.0*( pow( T_tm.err(id_C112)  ,2)  + pow( T_tm.err(id_C80) ,2)   )))));
      corr_fact_T_FSE_OS = ((T_OS.ave(id_C112) - T_OS.ave(id_C80))/(T_OS.ave(id_C80)))*fabs(erf( (T_OS.ave(id_C112)-T_OS.ave(id_C80))/(sqrt(2.0*( pow( T_OS.err(id_C112)  ,2)  + pow( T_OS.err(id_C80) ,2)   )))));

      corr_fact_err_T_FSE_tm = fabs(( (T_tm.distr_list[id_C112] - T_tm.distr_list[id_C80]).err()  )/(T_tm.ave(id_C80)))*fabs(erf( (T_tm.ave(id_C112)-T_tm.ave(id_C80))/(sqrt(2.0*( pow( T_tm.err(id_C112)  ,2)  + pow( T_tm.err(id_C80) ,2)   )))));
      corr_fact_err_T_FSE_OS = fabs(( (T_OS.distr_list[id_C112] - T_OS.distr_list[id_C80]).err() )/(T_OS.ave(id_C80)))*fabs(erf( (T_OS.ave(id_C112)-T_OS.ave(id_C80))/(sqrt(2.0*( pow( T_OS.err(id_C112)  ,2)  + pow( T_OS.err(id_C80) ,2)   )))));

      pl_T_tm=  (T_tm.ave(id_C112)-T_tm.ave(id_C80))/(sqrt(( pow( T_tm.err(id_C112)  ,2)  + pow( T_tm.err(id_C80) ,2)   )));
      pl_T_OS=  (T_OS.ave(id_C112)-T_OS.ave(id_C80))/(sqrt(( pow( T_OS.err(id_C112)  ,2)  + pow( T_OS.err(id_C80) ,2)   )));
      
      corr_fact_T_FSE = max( fabs(corr_fact_T_FSE_tm), fabs(corr_fact_T_FSE_OS));

      //L
      corr_fact_L_FSE_tm = ((L_tm.ave(id_C112) - L_tm.ave(id_C80))/(L_tm.ave(id_C80)))*fabs(erf( (L_tm.ave(id_C112)-L_tm.ave(id_C80))/(sqrt(2.0*( pow( L_tm.err(id_C112)  ,2)  + pow( L_tm.err(id_C80) ,2)   )))));
      corr_fact_L_FSE_OS = ((L_OS.ave(id_C112) - L_OS.ave(id_C80))/(L_OS.ave(id_C80)))*fabs(erf( (L_OS.ave(id_C112)-L_OS.ave(id_C80))/(sqrt(2.0*( pow( L_OS.err(id_C112)  ,2)  + pow( L_OS.err(id_C80) ,2)   )))));

      pl_L_tm = (L_tm.ave(id_C112)-L_tm.ave(id_C80))/(sqrt(( pow( L_tm.err(id_C112)  ,2)  + pow( L_tm.err(id_C80) ,2))));
      pl_L_OS = (L_OS.ave(id_C112)-L_OS.ave(id_C80))/(sqrt(( pow( L_OS.err(id_C112)  ,2)  + pow( L_OS.err(id_C80) ,2))));
      
      corr_fact_err_L_FSE_tm = fabs((( (L_tm.distr_list[id_C112] - L_tm.distr_list[id_C80]).err()  )/(L_tm.ave(id_C80)))*erf( (L_tm.ave(id_C112)-L_tm.ave(id_C80))/(sqrt(2.0*( pow( L_tm.err(id_C112)  ,2)  + pow( L_tm.err(id_C80) ,2)   )))));
      corr_fact_err_L_FSE_OS = fabs((( (L_OS.distr_list[id_C112] - L_OS.distr_list[id_C80]).err() )/(L_OS.ave(id_C80)))*erf( (L_OS.ave(id_C112)-L_OS.ave(id_C80))/(sqrt(2.0*( pow( L_OS.err(id_C112)  ,2)  + pow( L_OS.err(id_C80) ,2)   )))));
      
      corr_fact_L_FSE = max(fabs(corr_fact_L_FSE_tm), fabs(corr_fact_L_FSE_OS));


      //generate jackknife distribution for FSE corrections

      distr_t distr_syst_FSE_VkVk(UseJack);
      distr_t distr_syst_FSE_AkAk(UseJack);
      distr_t distr_syst_FSE_A0A0(UseJack);
      distr_t distr_syst_FSE_V0V0(UseJack);
      distr_t distr_syst_FSE_T(UseJack);
      distr_t distr_syst_FSE_L(UseJack);

      for(int ijack=0; ijack<Njacks;ijack++) {

	distr_syst_FSE_VkVk.distr.push_back( 1.0 + GM_sigma()*corr_fact_VkVk_FSE/sqrt(Njacks-1.0));
	distr_syst_FSE_AkAk.distr.push_back( 1.0 + GM_sigma()*corr_fact_AkAk_FSE/sqrt(Njacks-1.0));
	distr_syst_FSE_A0A0.distr.push_back( 1.0 + GM_sigma()*corr_fact_A0A0_FSE/sqrt(Njacks-1.0));
	distr_syst_FSE_V0V0.distr.push_back( 1.0 + GM_sigma()*corr_fact_V0V0_FSE/sqrt(Njacks-1.0));
	distr_syst_FSE_T.distr.push_back( 1.0 + GM_sigma()*corr_fact_T_FSE/sqrt(Njacks-1.0));
	distr_syst_FSE_L.distr.push_back( 1.0 + GM_sigma()*corr_fact_L_FSE/sqrt(Njacks-1.0));
	

      }


      //################### COMPUTATION OF FSEs CORRECTION FROM B64-B96 SPREAD##########################

      //A0A0
      corr_fact_BIS_A0A0_FSE_tm = ((A0A0_tm.ave(id_B96) - A0A0_tm.ave(id_B64))/(A0A0_tm.ave(id_B64)))*fabs(erf( (A0A0_tm.ave(id_B96)-A0A0_tm.ave(id_B64))/(sqrt(2.0*( pow( A0A0_tm.err(id_B96)  ,2)  + pow( A0A0_tm.err(id_B64) ,2)   )))));
      corr_fact_BIS_A0A0_FSE_OS = ((A0A0_OS.ave(id_B96) - A0A0_OS.ave(id_B64))/(A0A0_OS.ave(id_B64)))*fabs(erf( (A0A0_OS.ave(id_B96)-A0A0_OS.ave(id_B64))/(sqrt(2.0*( pow( A0A0_OS.err(id_B96)  ,2)  + pow( A0A0_OS.err(id_B64) ,2)   )))));

      corr_fact_BIS_err_A0A0_FSE_tm = fabs((( (A0A0_tm.distr_list[id_B96] - A0A0_tm.distr_list[id_B64]).err())/(A0A0_tm.ave(id_B64)))*erf( (A0A0_tm.ave(id_B96)-A0A0_tm.ave(id_B64))/(sqrt(2.0*( pow( A0A0_tm.err(id_B96)  ,2)  + pow( A0A0_tm.err(id_B64) ,2)   )))));
      corr_fact_BIS_err_A0A0_FSE_OS = fabs((( (A0A0_OS.distr_list[id_B96] - A0A0_OS.distr_list[id_B64]).err())/(A0A0_OS.ave(id_B64)))*erf( (A0A0_OS.ave(id_B96)-A0A0_OS.ave(id_B64))/(sqrt(2.0*( pow( A0A0_OS.err(id_B96)  ,2)  + pow( A0A0_OS.err(id_B64) ,2)   )))));
      
      corr_fact_BIS_A0A0_FSE = max( fabs(corr_fact_BIS_A0A0_FSE_tm), fabs(corr_fact_BIS_A0A0_FSE_OS));
    
      //V0V0
      corr_fact_BIS_V0V0_FSE_tm = ((V0V0_tm.ave(id_B96) - V0V0_tm.ave(id_B64))/(V0V0_tm.ave(id_B64)))*fabs(erf( (V0V0_tm.ave(id_B96)-V0V0_tm.ave(id_B64))/(sqrt(2.0*( pow( V0V0_tm.err(id_B96)  ,2)  + pow( V0V0_tm.err(id_B64) ,2)   )))));
      corr_fact_BIS_V0V0_FSE_OS = ((V0V0_OS.ave(id_B96) - V0V0_OS.ave(id_B64))/(V0V0_OS.ave(id_B64)))*fabs(erf( (V0V0_OS.ave(id_B96)-V0V0_OS.ave(id_B64))/(sqrt(2.0*( pow( V0V0_OS.err(id_B96)  ,2)  + pow( V0V0_OS.err(id_B64) ,2)   )))));

      corr_fact_BIS_err_V0V0_FSE_tm = fabs((( (V0V0_tm.distr_list[id_B96] - V0V0_tm.distr_list[id_B64]).err())/(V0V0_tm.ave(id_B64)))*erf( (V0V0_tm.ave(id_B96)-V0V0_tm.ave(id_B64))/(sqrt(2.0*( pow( V0V0_tm.err(id_B96)  ,2)  + pow( V0V0_tm.err(id_B64) ,2)   )))));
      corr_fact_BIS_err_V0V0_FSE_OS = fabs((( (V0V0_OS.distr_list[id_B96] - V0V0_OS.distr_list[id_B64]).err())/(V0V0_OS.ave(id_B64)))*erf( (V0V0_OS.ave(id_B96)-V0V0_OS.ave(id_B64))/(sqrt(2.0*( pow( V0V0_OS.err(id_B96)  ,2)  + pow( V0V0_OS.err(id_B64) ,2)   )))));
      
      corr_fact_BIS_V0V0_FSE = max(fabs(corr_fact_BIS_V0V0_FSE_tm), fabs(corr_fact_BIS_V0V0_FSE_OS));
           
      

      //AkAk
      corr_fact_BIS_AkAk_FSE_tm = ((AkAk_tm.ave(id_B96) - AkAk_tm.ave(id_B64))/(AkAk_tm.ave(id_B64)))*fabs(erf( (AkAk_tm.ave(id_B96)-AkAk_tm.ave(id_B64))/(sqrt(2.0*( pow( AkAk_tm.err(id_B96)  ,2)  + pow( AkAk_tm.err(id_B64) ,2)   )))));
      corr_fact_BIS_AkAk_FSE_OS = ((AkAk_OS.ave(id_B96) - AkAk_OS.ave(id_B64))/(AkAk_OS.ave(id_B64)))*fabs(erf( (AkAk_OS.ave(id_B96)-AkAk_OS.ave(id_B64))/(sqrt(2.0*( pow( AkAk_OS.err(id_B96)  ,2)  + pow( AkAk_OS.err(id_B64) ,2)   )))));

      corr_fact_BIS_err_AkAk_FSE_tm = fabs((( (AkAk_tm.distr_list[id_B96] - AkAk_tm.distr_list[id_B64]).err())/(AkAk_tm.ave(id_B64)))*erf( (AkAk_tm.ave(id_B96)-AkAk_tm.ave(id_B64))/(sqrt(2.0*( pow( AkAk_tm.err(id_B96)  ,2)  + pow( AkAk_tm.err(id_B64) ,2)   )))));
      corr_fact_BIS_err_AkAk_FSE_OS = fabs((( (AkAk_OS.distr_list[id_B96] - AkAk_OS.distr_list[id_B64]).err())/(AkAk_OS.ave(id_B64)))*erf( (AkAk_OS.ave(id_B96)-AkAk_OS.ave(id_B64))/(sqrt(2.0*( pow( AkAk_OS.err(id_B96)  ,2)  + pow( AkAk_OS.err(id_B64) ,2)   )))));

      
      corr_fact_BIS_AkAk_FSE = max(fabs(corr_fact_BIS_AkAk_FSE_tm), fabs(corr_fact_BIS_AkAk_FSE_OS));
   
      

      //VkVk
      corr_fact_BIS_VkVk_FSE_tm = ((VkVk_tm.ave(id_B96) - VkVk_tm.ave(id_B64))/(VkVk_tm.ave(id_B64)))*fabs(erf( (VkVk_tm.ave(id_B96)-VkVk_tm.ave(id_B64))/(sqrt(2.0*( pow( VkVk_tm.err(id_B96)  ,2)  + pow( VkVk_tm.err(id_B64) ,2)   )))));
      corr_fact_BIS_VkVk_FSE_OS = ((VkVk_OS.ave(id_B96) - VkVk_OS.ave(id_B64))/(VkVk_OS.ave(id_B64)))*fabs(erf( (VkVk_OS.ave(id_B96)-VkVk_OS.ave(id_B64))/(sqrt(2.0*( pow( VkVk_OS.err(id_B96)  ,2)  + pow( VkVk_OS.err(id_B64) ,2)   )))));

      corr_fact_BIS_err_VkVk_FSE_tm = fabs((( (VkVk_tm.distr_list[id_B96] - VkVk_tm.distr_list[id_B64]).err())/(VkVk_tm.ave(id_B64)))*erf( (VkVk_tm.ave(id_B96)-VkVk_tm.ave(id_B64))/(sqrt(2.0*( pow( VkVk_tm.err(id_B96)  ,2)  + pow( VkVk_tm.err(id_B64) ,2)   )))));
      corr_fact_BIS_err_VkVk_FSE_OS = fabs((( (VkVk_OS.distr_list[id_B96] - VkVk_OS.distr_list[id_B64]).err())/(VkVk_OS.ave(id_B64)))*erf( (VkVk_OS.ave(id_B96)-VkVk_OS.ave(id_B64))/(sqrt(2.0*( pow( VkVk_OS.err(id_B96)  ,2)  + pow( VkVk_OS.err(id_B64) ,2)   )))));
      
      corr_fact_BIS_VkVk_FSE = max(fabs(corr_fact_BIS_VkVk_FSE_tm),fabs(corr_fact_BIS_VkVk_FSE_OS));

      //T
      corr_fact_BIS_T_FSE_tm = ((T_tm.ave(id_B96) - T_tm.ave(id_B64))/(T_tm.ave(id_B64)))*fabs(erf( (T_tm.ave(id_B96)-T_tm.ave(id_B64))/(sqrt(2.0*( pow( T_tm.err(id_B96)  ,2)  + pow( T_tm.err(id_B64) ,2)   )))));
      corr_fact_BIS_T_FSE_OS = ((T_OS.ave(id_B96) - T_OS.ave(id_B64))/(T_OS.ave(id_B64)))*fabs(erf( (T_OS.ave(id_B96)-T_OS.ave(id_B64))/(sqrt(2.0*( pow( T_OS.err(id_B96)  ,2)  + pow( T_OS.err(id_B64) ,2)   )))));

      pl_BIS_T_tm =  (T_tm.ave(id_B96)-T_tm.ave(id_B64))/(sqrt(( pow( T_tm.err(id_B96)  ,2)  + pow( T_tm.err(id_B64) ,2)   )));
      pl_BIS_T_OS=   (T_OS.ave(id_B96)-T_OS.ave(id_B64))/(sqrt(( pow( T_OS.err(id_B96)  ,2)  + pow( T_OS.err(id_B64) ,2)   )));

      corr_fact_BIS_err_T_FSE_tm = fabs((( (T_tm.distr_list[id_B96] - T_tm.distr_list[id_B64]).err())/(T_tm.ave(id_B64)))*erf( (T_tm.ave(id_B96)-T_tm.ave(id_B64))/(sqrt(2.0*( pow( T_tm.err(id_B96)  ,2)  + pow( T_tm.err(id_B64) ,2)   )))));
      corr_fact_BIS_err_T_FSE_OS = fabs((( (T_OS.distr_list[id_B96] - T_OS.distr_list[id_B64]).err())/(T_OS.ave(id_B64)))*erf( (T_OS.ave(id_B96)-T_OS.ave(id_B64))/(sqrt(2.0*( pow( T_OS.err(id_B96)  ,2)  + pow( T_OS.err(id_B64) ,2)   )))));
      
      corr_fact_BIS_T_FSE = max(fabs(corr_fact_BIS_T_FSE_tm), fabs(corr_fact_BIS_T_FSE_OS));

      //L
      corr_fact_BIS_L_FSE_tm = ((L_tm.ave(id_B96) - L_tm.ave(id_B64))/(L_tm.ave(id_B64)))*fabs(erf( (L_tm.ave(id_B96)-L_tm.ave(id_B64))/(sqrt(2.0*( pow( L_tm.err(id_B96)  ,2)  + pow( L_tm.err(id_B64) ,2)   )))));
      corr_fact_BIS_L_FSE_OS = ((L_OS.ave(id_B96) - L_OS.ave(id_B64))/(L_OS.ave(id_B64)))*fabs(erf( (L_OS.ave(id_B96)-L_OS.ave(id_B64))/(sqrt(2.0*( pow( L_OS.err(id_B96)  ,2)  + pow( L_OS.err(id_B64) ,2)   )))));

      pl_BIS_L_tm= (L_tm.ave(id_B96)-L_tm.ave(id_B64))/(sqrt(( pow( L_tm.err(id_B96)  ,2)  + pow( L_tm.err(id_B64) ,2)   )));
      pl_BIS_L_OS= (L_OS.ave(id_B96)-L_OS.ave(id_B64))/(sqrt(( pow( L_OS.err(id_B96)  ,2)  + pow( L_OS.err(id_B64) ,2)   )));

      corr_fact_BIS_err_L_FSE_tm = fabs((( (L_tm.distr_list[id_B96] - L_tm.distr_list[id_B64]).err())/(L_tm.ave(id_B64)))*erf( (L_tm.ave(id_B96)-L_tm.ave(id_B64))/(sqrt(2.0*( pow( L_tm.err(id_B96)  ,2)  + pow( L_tm.err(id_B64) ,2)   )))));
      corr_fact_BIS_err_L_FSE_OS = fabs((( (L_OS.distr_list[id_B96] - L_OS.distr_list[id_B64]).err())/(L_OS.ave(id_B64)))*erf( (L_OS.ave(id_B96)-L_OS.ave(id_B64))/(sqrt(2.0*( pow( L_OS.err(id_B96)  ,2)  + pow( L_OS.err(id_B64) ,2)   )))));
      
      corr_fact_BIS_L_FSE = max(fabs(corr_fact_BIS_L_FSE_tm), fabs(corr_fact_BIS_L_FSE_OS));


      //#################################################################################################



      //PUSH_BACK FSEs corrections
      FSE_C_A0A0_tm.push_back( corr_fact_A0A0_FSE_tm);
      FSE_C_V0V0_tm.push_back( corr_fact_V0V0_FSE_tm);
      FSE_C_AkAk_tm.push_back( corr_fact_AkAk_FSE_tm);
      FSE_C_VkVk_tm.push_back( corr_fact_VkVk_FSE_tm);
      FSE_C_T_tm.push_back( corr_fact_T_FSE_tm);
      FSE_C_L_tm.push_back( corr_fact_L_FSE_tm);

      FSE_C_A0A0_OS.push_back( corr_fact_A0A0_FSE_OS);
      FSE_C_V0V0_OS.push_back( corr_fact_V0V0_FSE_OS);
      FSE_C_AkAk_OS.push_back( corr_fact_AkAk_FSE_OS);
      FSE_C_VkVk_OS.push_back( corr_fact_VkVk_FSE_OS);
      FSE_C_T_OS.push_back( corr_fact_T_FSE_OS);
      FSE_C_L_OS.push_back( corr_fact_L_FSE_OS);

      FSE_B_A0A0_tm.push_back( corr_fact_BIS_A0A0_FSE_tm);
      FSE_B_V0V0_tm.push_back( corr_fact_BIS_V0V0_FSE_tm);
      FSE_B_AkAk_tm.push_back( corr_fact_BIS_AkAk_FSE_tm);
      FSE_B_VkVk_tm.push_back( corr_fact_BIS_VkVk_FSE_tm);
      FSE_B_T_tm.push_back( corr_fact_BIS_T_FSE_tm);
      FSE_B_L_tm.push_back( corr_fact_BIS_L_FSE_tm);

      FSE_B_A0A0_OS.push_back( corr_fact_BIS_A0A0_FSE_OS);
      FSE_B_V0V0_OS.push_back( corr_fact_BIS_V0V0_FSE_OS);
      FSE_B_AkAk_OS.push_back( corr_fact_BIS_AkAk_FSE_OS);
      FSE_B_VkVk_OS.push_back( corr_fact_BIS_VkVk_FSE_OS);
      FSE_B_T_OS.push_back( corr_fact_BIS_T_FSE_OS);
      FSE_B_L_OS.push_back( corr_fact_BIS_L_FSE_OS);




      FSE_err_C_A0A0_tm.push_back( corr_fact_err_A0A0_FSE_tm);
      FSE_err_C_V0V0_tm.push_back( corr_fact_err_V0V0_FSE_tm);
      FSE_err_C_AkAk_tm.push_back( corr_fact_err_AkAk_FSE_tm);
      FSE_err_C_VkVk_tm.push_back( corr_fact_err_VkVk_FSE_tm);
      FSE_err_C_T_tm.push_back( corr_fact_err_T_FSE_tm);
      FSE_err_C_L_tm.push_back( corr_fact_err_L_FSE_tm);

      FSE_err_C_A0A0_OS.push_back( corr_fact_err_A0A0_FSE_OS);
      FSE_err_C_V0V0_OS.push_back( corr_fact_err_V0V0_FSE_OS);
      FSE_err_C_AkAk_OS.push_back( corr_fact_err_AkAk_FSE_OS);
      FSE_err_C_VkVk_OS.push_back( corr_fact_err_VkVk_FSE_OS);
      FSE_err_C_T_OS.push_back( corr_fact_err_T_FSE_OS);
      FSE_err_C_L_OS.push_back( corr_fact_err_L_FSE_OS);

      FSE_err_B_A0A0_tm.push_back( corr_fact_BIS_err_A0A0_FSE_tm);
      FSE_err_B_V0V0_tm.push_back( corr_fact_BIS_err_V0V0_FSE_tm);
      FSE_err_B_AkAk_tm.push_back( corr_fact_BIS_err_AkAk_FSE_tm);
      FSE_err_B_VkVk_tm.push_back( corr_fact_BIS_err_VkVk_FSE_tm);
      FSE_err_B_T_tm.push_back( corr_fact_BIS_err_T_FSE_tm);
      FSE_err_B_L_tm.push_back( corr_fact_BIS_err_L_FSE_tm);

      FSE_err_B_A0A0_OS.push_back( corr_fact_BIS_err_A0A0_FSE_OS);
      FSE_err_B_V0V0_OS.push_back( corr_fact_BIS_err_V0V0_FSE_OS);
      FSE_err_B_AkAk_OS.push_back( corr_fact_BIS_err_AkAk_FSE_OS);
      FSE_err_B_VkVk_OS.push_back( corr_fact_BIS_err_VkVk_FSE_OS);
      FSE_err_B_T_OS.push_back( corr_fact_BIS_err_T_FSE_OS);
      FSE_err_B_L_OS.push_back( corr_fact_BIS_err_L_FSE_OS);


      FSE_C_A0A0.push_back(  corr_fact_A0A0_FSE);
      FSE_C_V0V0.push_back(  corr_fact_V0V0_FSE);
      FSE_C_AkAk.push_back(  corr_fact_AkAk_FSE);
      FSE_C_VkVk.push_back(  corr_fact_VkVk_FSE);
      FSE_C_T.push_back(  corr_fact_T_FSE);
      FSE_C_L.push_back(  corr_fact_L_FSE);



      pl_T_tm_distr.push_back( pl_T_tm);
      pl_T_OS_distr.push_back ( pl_T_OS);
      pl_L_tm_distr.push_back( pl_L_tm);
      pl_L_OS_distr.push_back( pl_L_OS);

      pl_BIS_T_tm_distr.push_back( pl_BIS_T_tm);
      pl_BIS_T_OS_distr.push_back ( pl_BIS_T_OS);
      pl_BIS_L_tm_distr.push_back( pl_BIS_L_tm);
      pl_BIS_L_OS_distr.push_back( pl_BIS_L_OS);
      
      

      
      cout<<"corr_facts (tm, OS, max): "<<corr_fact_VkVk_FSE_tm<<" "<<corr_fact_VkVk_FSE_OS<<" "<<corr_fact_VkVk_FSE<<endl;
         
      for(int iens=0;iens<Nens;iens++) {

	if(ls_data_tm_VKVK.Tag[iens] != "cC211a.06.112") {
	  
	  if(ls_data_tm_VKVK.Tag[iens] != "cB211b.072.96") {


	    /*
	    
	    //A0A0	 
	    A0A0_tm_red.distr_list.push_back( A0A0_tm.ave(iens) + (A0A0_tm.distr_list[iens] - A0A0_tm.ave(iens))*(sqrt( pow(A0A0_tm.err(iens),2) + pow(corr_fact_A0A0_FSE*A0A0_tm.ave(iens),2))/A0A0_tm.err(iens)));
	    
	    A0A0_OS_red.distr_list.push_back( A0A0_OS.ave(iens) + (A0A0_OS.distr_list[iens] - A0A0_OS.ave(iens))*(sqrt( pow(A0A0_OS.err(iens),2) + pow(corr_fact_A0A0_FSE*A0A0_OS.ave(iens),2))/A0A0_OS.err(iens)));
	    //V0V0
	    V0V0_tm_red.distr_list.push_back( V0V0_tm.ave(iens) + (V0V0_tm.distr_list[iens] - V0V0_tm.ave(iens))*(sqrt( pow(V0V0_tm.err(iens),2) + pow(corr_fact_V0V0_FSE*V0V0_tm.ave(iens),2))/V0V0_tm.err(iens)));
	    V0V0_OS_red.distr_list.push_back( V0V0_OS.ave(iens) + (V0V0_OS.distr_list[iens] - V0V0_OS.ave(iens))*(sqrt( pow(V0V0_OS.err(iens),2) + pow(corr_fact_V0V0_FSE*V0V0_OS.ave(iens),2))/V0V0_OS.err(iens)));
	    //AkAk
	    AkAk_tm_red.distr_list.push_back( AkAk_tm.ave(iens) + (AkAk_tm.distr_list[iens] - AkAk_tm.ave(iens))*(sqrt( pow(AkAk_tm.err(iens),2) + pow(corr_fact_AkAk_FSE*AkAk_tm.ave(iens),2))/AkAk_tm.err(iens)));
	    AkAk_OS_red.distr_list.push_back( AkAk_OS.ave(iens) + (AkAk_OS.distr_list[iens] - AkAk_OS.ave(iens))*(sqrt( pow(AkAk_OS.err(iens),2) + pow(corr_fact_AkAk_FSE*AkAk_OS.ave(iens),2))/AkAk_OS.err(iens)));
	    //VkVk
	    VkVk_tm_red.distr_list.push_back( VkVk_tm.ave(iens) + (VkVk_tm.distr_list[iens] - VkVk_tm.ave(iens))*(sqrt( pow(VkVk_tm.err(iens),2) + pow(corr_fact_VkVk_FSE*VkVk_tm.ave(iens),2))/VkVk_tm.err(iens)));
	    VkVk_OS_red.distr_list.push_back( VkVk_OS.ave(iens) + (VkVk_OS.distr_list[iens] - VkVk_OS.ave(iens))*(sqrt( pow(VkVk_OS.err(iens),2) + pow(corr_fact_VkVk_FSE*VkVk_OS.ave(iens),2))/VkVk_OS.err(iens)));
	    //T
	    T_tm_red.distr_list.push_back( T_tm.ave(iens) + (T_tm.distr_list[iens] - T_tm.ave(iens))*(sqrt( pow(T_tm.err(iens),2) + pow(corr_fact_T_FSE*T_tm.ave(iens),2))/T_tm.err(iens)));
	    T_OS_red.distr_list.push_back( T_OS.ave(iens) + (T_OS.distr_list[iens] - T_OS.ave(iens))*(sqrt( pow(T_OS.err(iens),2) + pow(corr_fact_T_FSE*T_OS.ave(iens),2))/T_OS.err(iens)));
	    
	    //L
	    L_tm_red.distr_list.push_back( L_tm.ave(iens) + (L_tm.distr_list[iens] - L_tm.ave(iens))*(sqrt( pow(L_tm.err(iens),2) + pow(corr_fact_L_FSE*L_tm.ave(iens),2))/L_tm.err(iens)));
	    L_OS_red.distr_list.push_back( L_OS.ave(iens) + (L_OS.distr_list[iens] - L_OS.ave(iens))*(sqrt( pow(L_OS.err(iens),2) + pow(corr_fact_L_FSE*L_OS.ave(iens),2))/L_OS.err(iens)));
	    */

	    //A0A0
	    A0A0_tm_red.distr_list.push_back( A0A0_tm.distr_list[iens]*distr_syst_FSE_A0A0);
	    A0A0_OS_red.distr_list.push_back( A0A0_OS.distr_list[iens]*distr_syst_FSE_A0A0);

	    //V0V0
	    V0V0_tm_red.distr_list.push_back( V0V0_tm.distr_list[iens]*distr_syst_FSE_V0V0);
	    V0V0_OS_red.distr_list.push_back( V0V0_OS.distr_list[iens]*distr_syst_FSE_V0V0);

	    //AkAk
	    AkAk_tm_red.distr_list.push_back( AkAk_tm.distr_list[iens]*distr_syst_FSE_AkAk);
	    AkAk_OS_red.distr_list.push_back( AkAk_OS.distr_list[iens]*distr_syst_FSE_AkAk);

	    //VkVk
	    VkVk_tm_red.distr_list.push_back( VkVk_tm.distr_list[iens]*distr_syst_FSE_VkVk);
	    VkVk_OS_red.distr_list.push_back( VkVk_OS.distr_list[iens]*distr_syst_FSE_VkVk);

	    //T
	    T_tm_red.distr_list.push_back( T_tm.distr_list[iens]*distr_syst_FSE_T);
	    T_OS_red.distr_list.push_back( T_OS.distr_list[iens]*distr_syst_FSE_T);

	    //L
	    L_tm_red.distr_list.push_back( L_tm.distr_list[iens]*distr_syst_FSE_L);
	    L_OS_red.distr_list.push_back( L_OS.distr_list[iens]*distr_syst_FSE_L);

	    
	  }
	  else {
	    A0A0_tm_red.distr_list.push_back( A0A0_tm.distr_list[iens]);
	    A0A0_OS_red.distr_list.push_back( A0A0_OS.distr_list[iens]);
	    
	    V0V0_tm_red.distr_list.push_back( V0V0_tm.distr_list[iens]);
	    V0V0_OS_red.distr_list.push_back( V0V0_OS.distr_list[iens]);
	    
	    AkAk_tm_red.distr_list.push_back( AkAk_tm.distr_list[iens]);
	    AkAk_OS_red.distr_list.push_back( AkAk_OS.distr_list[iens]);

	    VkVk_tm_red.distr_list.push_back( VkVk_tm.distr_list[iens]);
	    VkVk_OS_red.distr_list.push_back( VkVk_OS.distr_list[iens]);
	    
	    T_tm_red.distr_list.push_back( T_tm.distr_list[iens]);
	    T_OS_red.distr_list.push_back( T_OS.distr_list[iens]);

	    L_tm_red.distr_list.push_back( L_tm.distr_list[iens]);
	    L_OS_red.distr_list.push_back( L_OS.distr_list[iens]);
	    
	    
	  }
	  
	  
          if(ls_data_tm_VKVK.Tag[iens].substr(1,1) == "B") a_distr_list_red.distr_list.push_back( a_B/fm_to_inv_Gev); //lattice spacing is in fm
	  
	  else if(ls_data_tm_VKVK.Tag[iens].substr(1,1) == "C") a_distr_list_red.distr_list.push_back( a_C/fm_to_inv_Gev); //lattice spacing is in fm
	  
	  else if(ls_data_tm_VKVK.Tag[iens].substr(1,1) == "D") a_distr_list_red.distr_list.push_back( a_D/fm_to_inv_Gev); //lattice spacing is in fm
	  
	  else if(ls_data_tm_VKVK.Tag[iens].substr(1,1) == "E") a_distr_list_red.distr_list.push_back( a_E/fm_to_inv_Gev); //lattice spacing is in fm

	  else crash("While building a_distr_list_red cannot recognize ensemble: "+ls_data_tm_VKVK.Tag[iens]);
	  
	  Tag_ens_red.push_back( ls_data_tm_VKVK.Tag[iens]);
	  
	}

	if(ls_data_tm_VKVK.Tag[iens].substr(1,1) == "B") a_distr_list.distr_list.push_back(a_B/fm_to_inv_Gev);  //lattice spacing is in fm
	else if(ls_data_tm_VKVK.Tag[iens].substr(1,1) == "C") a_distr_list.distr_list.push_back( a_C/fm_to_inv_Gev); //lattice spacing is in fm
	else if(ls_data_tm_VKVK.Tag[iens].substr(1,1) == "D") a_distr_list.distr_list.push_back( a_D/fm_to_inv_Gev); //lattice spacing is in fm
	else if(ls_data_tm_VKVK.Tag[iens].substr(1,1) == "E") a_distr_list.distr_list.push_back( a_E/fm_to_inv_Gev); //lattice spacing is in fm
	else crash("When building a_distr_list cannot recognize ensemble: "+ls_data_tm_VKVK.Tag[iens]);


	if(ls_data_tm_VKVK.Tag[iens].substr(1,1) == "E") { tot_TL_Ens_E_tm.distr_list.push_back( (T_tm+L_tm).distr_list[iens]); tot_TL_Ens_E_OS.distr_list.push_back( (T_OS+L_OS).distr_list[iens]); }
      }
      //###################################################

      
      cout<<"T tm size: "<<T_tm_red.size()<<endl<<flush;
      cout<<"a red size: "<<a_distr_list.size()<<endl<<flush;
      
  
      
      //now we are ready to perform the continuum limit extrapolation

      int Nens_eff=Nens-1;

     

      class ipar_TAU {
	
      public:
	ipar_TAU()  {}
	
	double Br, Br_err, a; //lattice spacing a is in fm
	bool Is_tm;
      };
      
      
      class fpar_TAU {

      public:
	fpar_TAU() {}
	fpar_TAU(const Vfloat &par) {
	  if((signed)par.size() != 5) crash("In class fpar_TAU, class constructor Vfloat par has size != 3");
	  D=par[0];
	  D2_tm=par[1];
	  D4_tm=par[2];
	  D2_OS=par[3];
	  D4_OS=par[4];
	}

	double D,D2_tm, D4_tm, D2_OS, D4_OS;
      };

       for( auto &contr: Contribs) {
	for( auto &fit_type: Fit_types) {
	  for( auto &poly_type: poly_types) {


	    if( (poly_type != "quad" && poly_type.substr(3,4) != "quad") || sigma_list_strange[is] > 0.13 ) {
	    
	    if( (fit_type != "comb") && (poly_type.substr(0,2)=="tm" || poly_type.substr(0,2) == "OS")) crash("Cannot use tm/OS_linear with fit type: "+fit_type);
	    cout<<"###########################################################"<<endl;
	    cout<<"Performing continuum limit extrapolation for sigma: "<<sigma_list_strange[is]<<endl;
	    cout<<"Contribution: "<<contr<<endl;
	    cout<<"Fit type: "<<fit_type<<endl;
	    cout<<"polynomial: "<<poly_type<<endl;

	    int Nmeas= ((fit_type=="comb")?(2*Nens_eff):Nens_eff);
	  
	    
	    bootstrap_fit<fpar_TAU,ipar_TAU> bf_TAU(Njacks);
	    bootstrap_fit<fpar_TAU,ipar_TAU> bf_TAU_ch2(1);
	    //bf_TAU.Disable_correlated_fit();
	    //bf_TAU_ch2.Disable_correlated_fit();
	    bf_TAU.Set_number_of_measurements(Nmeas);
	    bf_TAU.Set_verbosity(1);
	    //ch2
	    bf_TAU_ch2.Set_number_of_measurements(Nmeas);
	    bf_TAU_ch2.Set_verbosity(1);

	    //bf_TAU.set_warmup_lev(1);
	    //bf_TAU_ch2.set_warmup_lev(1);

	    //add fit parameters
	    bf_TAU.Add_par("D", 3.0, 0.1);
	    bf_TAU.Add_par("D2_tm", 2, 0.1);
	    bf_TAU.Add_par("D4_tm", 1, 0.1);
	    bf_TAU.Add_par("D2_OS", 2, 0.1);
	    bf_TAU.Add_par("D4_OS", 1, 0.1);
	    //ch2
	    bf_TAU_ch2.Add_par("D", 3.0, 0.1);
	    bf_TAU_ch2.Add_par("D2_tm", 2, 0.1);
	    bf_TAU_ch2.Add_par("D4_tm", 1, 0.1);
	    bf_TAU_ch2.Add_par("D2_OS", 2, 0.1);
	    bf_TAU_ch2.Add_par("D4_OS", 1, 0.1);

	    //fix parameters depending on fit type
	    if(poly_type=="const") {
	      bf_TAU.Fix_par("D2_tm", 0);
	      bf_TAU.Fix_par("D2_OS", 0);
	      //ch2
	      bf_TAU_ch2.Fix_par("D2_tm", 0);
	      bf_TAU_ch2.Fix_par("D2_OS", 0);

	      bf_TAU.Fix_par("D4_tm",0); bf_TAU_ch2.Fix_par("D4_tm",0);
	      bf_TAU.Fix_par("D4_OS",0); bf_TAU_ch2.Fix_par("D4_OS",0);
	    }
	    else if(poly_type=="linear") {
	      bf_TAU.Fix_par("D4_tm",0); bf_TAU_ch2.Fix_par("D4_tm",0);
	      bf_TAU.Fix_par("D4_OS",0); bf_TAU_ch2.Fix_par("D4_OS",0);
	      
	      if(fit_type=="OS") { bf_TAU.Fix_par("D2_tm", 0); bf_TAU_ch2.Fix_par("D2_tm", 0); }
	      else if(fit_type=="tm") { bf_TAU.Fix_par("D2_OS",0); bf_TAU_ch2.Fix_par("D2_OS",0); }
	    }
	    else if(poly_type=="tm_linear") {
	      bf_TAU.Fix_par("D2_OS",0); bf_TAU_ch2.Fix_par("D2_OS",0);

	      bf_TAU.Fix_par("D4_tm",0); bf_TAU_ch2.Fix_par("D4_tm",0);
	      bf_TAU.Fix_par("D4_OS",0); bf_TAU_ch2.Fix_par("D4_OS",0);
	    }
	    else if(poly_type=="OS_linear") {
	       bf_TAU.Fix_par("D2_tm",0); bf_TAU_ch2.Fix_par("D2_tm",0);

	       bf_TAU.Fix_par("D4_tm",0); bf_TAU_ch2.Fix_par("D4_tm",0);
	       bf_TAU.Fix_par("D4_OS",0); bf_TAU_ch2.Fix_par("D4_OS",0);
	    }
	    else if(poly_type=="tm_quad") {
	      bf_TAU.Fix_par("D4_OS",0); bf_TAU_ch2.Fix_par("D4_OS",0);
	    }
	    else if(poly_type=="OS_quad") {
	      bf_TAU.Fix_par("D4_tm",0); bf_TAU_ch2.Fix_par("D4_tm",0);
	    }
	    else if(poly_type=="quad") {
	      //fix nothing
	    }
	    else crash("poly_type: "+poly_type+" not yet implemented");



	    //Depending on fit considered, determine Nmeas, Npars and Ndof
	   
	    int Npars= bf_TAU.Get_number_of_fit_pars();
	    int Ndof= Nmeas-Npars;
	    cout<<"poly type: "<<poly_type<<endl;
	    cout<<"Nmeas: "<<Nmeas<<endl;
	    cout<<"Npars: "<<Npars<<endl;
	    cout<<"Ndof: "<<Ndof<<endl;
	    cout<<"Nens: "<<Nens_eff<<endl;


	    //ansatz
	    bf_TAU.ansatz=  [ ](const fpar_TAU &p, const ipar_TAU &ip) {
	      double D2=0.0;
	      double D4=0.0;
	      if( ip.Is_tm==true ) { D2=p.D2_tm; D4=p.D4_tm;}
	      else {D2=p.D2_OS; D4=p.D4_OS;}
	      
	      return p.D + D2*pow(ip.a*QCD_scale,2) + D4*pow(ip.a*QCD_scale,4);
	    };
	    //meas
	    bf_TAU.measurement=  [ ](const fpar_TAU &p, const ipar_TAU &ip) {
	      return ip.Br;
	    };
	    //err
	    bf_TAU.error=  [ ](const fpar_TAU &p, const ipar_TAU &ip) {
	      return ip.Br_err;
	    };
	    //ch2
	    bf_TAU_ch2.ansatz= bf_TAU.ansatz;
	    bf_TAU_ch2.measurement= bf_TAU.measurement;
	    bf_TAU_ch2.error= bf_TAU.error;


	 

	    //fill the data
	    int off_OS = ((fit_type=="comb")?(Nens_eff):0);
	    vector<vector<ipar_TAU>> data(Njacks);
	    vector<vector<ipar_TAU>> data_ch2(1);
	    //allocate space for output result
	    boot_fit_data<fpar_TAU> Bt_fit;
	    boot_fit_data<fpar_TAU> Bt_fit_ch2;
	    for(auto &data_iboot: data) data_iboot.resize(Nmeas);
	    for(auto &data_iboot: data_ch2) data_iboot.resize(Nmeas);

		    
	    //if fit type is "comb" insert covariance matrix
	    if(fit_type=="comb") {

	      Eigen::MatrixXd Cov_Matrix(Nmeas,Nmeas);
	      Eigen::MatrixXd Corr_Matrix(Nmeas,Nmeas);
	      for(int i=0;i<Nmeas;i++) for(int j=0;j<Nmeas;j++) {Cov_Matrix(i,j)=0; Corr_Matrix(i,j)=0;}


	      //compute cov matrix between tm and OS
	      for(int iens=0; iens<Nens_eff;iens++) {

		Corr_Matrix(iens,iens) = 1;
		Corr_Matrix(iens+off_OS,iens+off_OS) = 1;
		
		if(contr=="A0A0") {
		  Cov_Matrix(iens,iens) = pow(A0A0_tm_red.err(iens),2); Cov_Matrix(iens+off_OS,iens+off_OS) = pow(A0A0_OS_red.err(iens),2);
		  Cov_Matrix(iens, off_OS+iens) = (A0A0_tm_red.distr_list[iens]%A0A0_OS_red.distr_list[iens]);
		}
		else if(contr=="V0V0") {
		  Cov_Matrix(iens,iens) = pow(V0V0_tm_red.err(iens),2); Cov_Matrix(iens+off_OS,iens+off_OS) = pow(V0V0_OS_red.err(iens),2);
		  Cov_Matrix(iens, off_OS+iens) = (V0V0_tm_red.distr_list[iens]%V0V0_OS_red.distr_list[iens]);
		}
		else if(contr=="AkAk") {
		  Cov_Matrix(iens,iens) = pow(AkAk_tm_red.err(iens),2); Cov_Matrix(iens+off_OS,iens+off_OS) = pow(AkAk_OS_red.err(iens),2); 
		  Cov_Matrix(iens, off_OS+iens) = (AkAk_tm_red.distr_list[iens]%AkAk_OS_red.distr_list[iens]);
		}

		else if(contr=="VkVk") {
		  Cov_Matrix(iens,iens) = pow(VkVk_tm_red.err(iens),2); Cov_Matrix(iens+off_OS,iens+off_OS) = pow(VkVk_OS_red.err(iens),2);
		  Cov_Matrix(iens, off_OS+iens) = (VkVk_tm_red.distr_list[iens]%VkVk_OS_red.distr_list[iens]);
		}

		else if(contr=="VA") {
		  Cov_Matrix(iens,iens) = pow(((VkVk_tm_red+V0V0_tm_red)/(AkAk_tm_red+A0A0_tm_red)).err(iens),2); Cov_Matrix(iens+off_OS,iens+off_OS) = pow(((VkVk_OS_red+V0V0_OS_red)/(AkAk_OS_red+A0A0_OS_red)).err(iens),2);
		  Cov_Matrix(iens, off_OS+iens) = (((VkVk_tm_red+V0V0_tm_red)/(AkAk_tm_red+A0A0_tm_red)).distr_list[iens]%((VkVk_OS_red+V0V0_OS_red)/(AkAk_OS_red+A0A0_OS_red)).distr_list[iens]);
		}

		else if(contr=="VMA") {
		  Cov_Matrix(iens,iens) = pow((( VkVk_tm_red+V0V0_tm_red-AkAk_tm_red-A0A0_tm_red )/(VkVk_tm_red+AkAk_tm_red+A0A0_tm_red+V0V0_tm_red)).err(iens),2); Cov_Matrix(iens+off_OS,iens+off_OS) = pow( (( VkVk_OS_red+V0V0_OS_red-AkAk_OS_red-A0A0_OS_red )/(VkVk_OS_red+AkAk_OS_red+A0A0_OS_red+V0V0_OS_red)).err(iens),2);
		  Cov_Matrix(iens, off_OS+iens) = ((( VkVk_tm_red+V0V0_tm_red-AkAk_tm_red-A0A0_tm_red )/(VkVk_tm_red+AkAk_tm_red+A0A0_tm_red+V0V0_tm_red)).distr_list[iens]%(( VkVk_OS_red+V0V0_OS_red-AkAk_OS_red-A0A0_OS_red )/(VkVk_OS_red+AkAk_OS_red+A0A0_OS_red+V0V0_OS_red)).distr_list[iens]);
		}

		else if(contr=="AX") {
		  Cov_Matrix(iens,iens) = pow( (AkAk_tm_red+A0A0_tm_red).err(iens),2); Cov_Matrix(iens+off_OS,iens+off_OS) = pow( (AkAk_OS_red+ A0A0_OS_red).err(iens),2); 
		  Cov_Matrix(iens, off_OS+iens) = ( (AkAk_tm_red+A0A0_tm_red).distr_list[iens]% (AkAk_OS_red+A0A0_OS_red).distr_list[iens]);
		}

		else if(contr=="T") {
		  Cov_Matrix(iens,iens) = pow( (T_tm_red).err(iens),2); Cov_Matrix(iens+off_OS,iens+off_OS) = pow( (T_OS_red).err(iens),2); 
		  Cov_Matrix(iens, off_OS+iens) = ( (T_tm_red).distr_list[iens]% (T_OS_red).distr_list[iens]);
		}

		else if(contr=="L") {
		  Cov_Matrix(iens,iens) = pow( (L_tm_red).err(iens),2); Cov_Matrix(iens+off_OS,iens+off_OS) = pow( (L_OS_red).err(iens),2); 
		  Cov_Matrix(iens, off_OS+iens) = ( (L_tm_red).distr_list[iens]% (L_OS_red).distr_list[iens]);
		}
		

		else if(contr=="tot") {
		  Cov_Matrix(iens,iens) = pow((VkVk_tm_red+AkAk_tm_red+A0A0_tm_red+ V0V0_tm_red).err(iens),2); Cov_Matrix(iens+off_OS,iens+off_OS) = pow( (VkVk_OS_red+AkAk_OS_red+A0A0_OS_red+V0V0_OS_red).err(iens),2);
		  Cov_Matrix(iens, off_OS+iens) = ((VkVk_tm_red+AkAk_tm_red+A0A0_tm_red+V0V0_tm_red).distr_list[iens]%(VkVk_OS_red+AkAk_OS_red+A0A0_OS_red+V0V0_OS_red).distr_list[iens]);
		}

		else if(contr=="tot_TL") {
		  Cov_Matrix(iens,iens) = pow((T_tm_red+ L_tm_red).err(iens),2); Cov_Matrix(iens+off_OS,iens+off_OS) = pow( (T_OS_red+L_OS_red).err(iens),2);
		  Cov_Matrix(iens, off_OS+iens) = ((T_tm_red+L_tm_red).distr_list[iens]%(T_OS_red+L_OS_red).distr_list[iens]);
		}
		else crash("contr: "+contr+" not recognized");

		Corr_Matrix(iens, off_OS+iens) = Cov_Matrix(iens,off_OS+iens)/sqrt( Cov_Matrix(iens,iens)*Cov_Matrix(off_OS+iens,off_OS+iens));

		
		//symmetrize
		Cov_Matrix(off_OS+iens, iens) = Cov_Matrix(iens, off_OS+iens);
		Corr_Matrix(off_OS+iens, iens) = Corr_Matrix(iens, off_OS+iens);
		

	      }

	 	    	      
	      //add cov matrix to bootstrap fit
	      bf_TAU.Add_covariance_matrix(Cov_Matrix);
	      bf_TAU_ch2.Add_covariance_matrix(Cov_Matrix);

	      //print covariance matrix
	      ofstream Print_Cov("../data/tau_decay/"+Tag_reco_type+"/strange/continuum/cov/"+contr+"_sigma_"+to_string_with_precision(sigma_list_strange[is],3)+".cov");
	      ofstream Print_Corr("../data/tau_decay/"+Tag_reco_type+"/strange/continuum/corr/"+contr+"_sigma_"+to_string_with_precision(sigma_list_strange[is],3)+".corr");

	      Print_Cov<<Cov_Matrix<<endl;  Print_Corr<<Corr_Matrix<<endl;
	      Print_Cov.close();            Print_Corr.close();
	    	    	      
	    }


	    
	    for(int ijack=0;ijack<Njacks;ijack++) {
	      for(int iens=0;iens<Nens_eff;iens++) {

		if(contr=="A0A0") {

		
		  if((fit_type=="tm") || (fit_type=="comb")) {
		   
		    data[ijack][iens].Br = A0A0_tm_red.distr_list[iens].distr[ijack];
		    data[ijack][iens].Br_err= A0A0_tm_red.err(iens);
		    data[ijack][iens].Is_tm = true;
		    data[ijack][iens].a = a_distr_list_red.distr_list[iens].distr[ijack];
		  }

		  if( (fit_type=="OS") || (fit_type=="comb") ) {
		    data[ijack][iens+off_OS].Br = A0A0_OS_red.distr_list[iens].distr[ijack];
		    data[ijack][iens+off_OS].Br_err= A0A0_OS_red.err(iens);
		    data[ijack][iens+off_OS].Is_tm = false;
		    data[ijack][iens+off_OS].a = a_distr_list_red.distr_list[iens].distr[ijack];
		  }

		}

		else if(contr=="V0V0") {

		
		  if((fit_type=="tm") || (fit_type=="comb")) {
		   
		    data[ijack][iens].Br = V0V0_tm_red.distr_list[iens].distr[ijack];
		    data[ijack][iens].Br_err= V0V0_tm_red.err(iens);
		    data[ijack][iens].Is_tm = true;
		    data[ijack][iens].a = a_distr_list_red.distr_list[iens].distr[ijack];
		  }

		  if( (fit_type=="OS") || (fit_type=="comb") ) {
		    data[ijack][iens+off_OS].Br = V0V0_OS_red.distr_list[iens].distr[ijack];
		    data[ijack][iens+off_OS].Br_err= V0V0_OS_red.err(iens);
		    data[ijack][iens+off_OS].Is_tm = false;
		    data[ijack][iens+off_OS].a = a_distr_list_red.distr_list[iens].distr[ijack];
		  }

		}
		
		else if(contr=="AkAk") {
		   if((fit_type=="tm") || (fit_type=="comb")) {
		    data[ijack][iens].Br = AkAk_tm_red.distr_list[iens].distr[ijack];
		    data[ijack][iens].Br_err= AkAk_tm_red.err(iens);
		    data[ijack][iens].Is_tm = true;
		    data[ijack][iens].a = a_distr_list_red.distr_list[iens].distr[ijack];
		   }

		   if((fit_type=="OS") || (fit_type=="comb")) {
		    data[ijack][iens+off_OS].Br = AkAk_OS_red.distr_list[iens].distr[ijack];
		    data[ijack][iens+off_OS].Br_err= AkAk_OS_red.err(iens);
		    data[ijack][iens+off_OS].Is_tm = false;
		    data[ijack][iens+off_OS].a = a_distr_list_red.distr_list[iens].distr[ijack];
		   }
		   
		}

		else if(contr=="AX") {
		   if((fit_type=="tm") || (fit_type=="comb")) {
		     data[ijack][iens].Br = (AkAk_tm_red+A0A0_tm_red).distr_list[iens].distr[ijack];
		     data[ijack][iens].Br_err= (AkAk_tm_red+A0A0_tm_red).err(iens);
		     data[ijack][iens].Is_tm = true;
		     data[ijack][iens].a = a_distr_list_red.distr_list[iens].distr[ijack];
		   }

		   if((fit_type=="OS") || (fit_type=="comb")) {
		    data[ijack][iens+off_OS].Br = (AkAk_OS_red+A0A0_OS_red).distr_list[iens].distr[ijack];
		    data[ijack][iens+off_OS].Br_err= (AkAk_OS_red+A0A0_OS_red).err(iens);
		    data[ijack][iens+off_OS].Is_tm = false;
		    data[ijack][iens+off_OS].a = a_distr_list_red.distr_list[iens].distr[ijack];
		   }
		   
		}

		else if(contr=="VkVk") {
		   if((fit_type=="tm") || (fit_type=="comb")) {
		    data[ijack][iens].Br = VkVk_tm_red.distr_list[iens].distr[ijack];
		    data[ijack][iens].Br_err= VkVk_tm_red.err(iens);
		    data[ijack][iens].Is_tm = true;
		    data[ijack][iens].a = a_distr_list_red.distr_list[iens].distr[ijack];
		  }
		   
		   if((fit_type=="OS") || (fit_type=="comb")) {
		    data[ijack][iens+off_OS].Br = VkVk_OS_red.distr_list[iens].distr[ijack];
		    data[ijack][iens+off_OS].Br_err= VkVk_OS_red.err(iens);
		    data[ijack][iens+off_OS].Is_tm = false;
		    data[ijack][iens+off_OS].a = a_distr_list_red.distr_list[iens].distr[ijack];
		   }
		   
		}

		else if(contr=="VA") {
		  if((fit_type=="tm") || (fit_type=="comb")) {
		    data[ijack][iens].Br = ( (VkVk_tm_red+V0V0_tm_red)/(AkAk_tm_red+A0A0_tm_red)).distr_list[iens].distr[ijack];
		    data[ijack][iens].Br_err= ((VkVk_tm_red+V0V0_tm_red)/(AkAk_tm_red+A0A0_tm_red)).err(iens);
		    data[ijack][iens].Is_tm = true;
		    data[ijack][iens].a = a_distr_list_red.distr_list[iens].distr[ijack];
		  }
		  
		  if((fit_type=="OS") || (fit_type=="comb")) {
		    data[ijack][iens+off_OS].Br = ((VkVk_OS_red+V0V0_OS_red)/(AkAk_OS_red+A0A0_OS_red)).distr_list[iens].distr[ijack];
		    data[ijack][iens+off_OS].Br_err= ((VkVk_OS_red+V0V0_OS_red)/(AkAk_OS_red+A0A0_OS_red)).err(iens);
		    data[ijack][iens+off_OS].Is_tm = false;
		    data[ijack][iens+off_OS].a = a_distr_list_red.distr_list[iens].distr[ijack];
		  }
		  
		}

		else if(contr=="VMA") {
		  if((fit_type=="tm") || (fit_type=="comb")) {
		    data[ijack][iens].Br = ( ( VkVk_tm_red+V0V0_tm_red-AkAk_tm_red-A0A0_tm_red )/(VkVk_tm_red+AkAk_tm_red+A0A0_tm_red+V0V0_tm_red)).distr_list[iens].distr[ijack];
		    data[ijack][iens].Br_err= ( ( VkVk_tm_red+V0V0_tm_red-AkAk_tm_red-A0A0_tm_red )/(VkVk_tm_red+AkAk_tm_red+A0A0_tm_red+V0V0_tm_red)).err(iens);
		    data[ijack][iens].Is_tm = true;
		    data[ijack][iens].a = a_distr_list_red.distr_list[iens].distr[ijack];
		  }

		  if((fit_type=="OS") || (fit_type=="comb")) {
		    data[ijack][iens+off_OS].Br = ( ( VkVk_OS_red+V0V0_OS_red-AkAk_OS_red-A0A0_OS_red )/(VkVk_OS_red+AkAk_OS_red+A0A0_OS_red+V0V0_OS_red)).distr_list[iens].distr[ijack];
		    data[ijack][iens+off_OS].Br_err=  ( ( VkVk_OS_red+V0V0_OS_red-AkAk_OS_red-A0A0_OS_red )/(VkVk_OS_red+AkAk_OS_red+A0A0_OS_red+V0V0_OS_red)).err(iens);
		    data[ijack][iens+off_OS].Is_tm = false;
		    data[ijack][iens+off_OS].a = a_distr_list_red.distr_list[iens].distr[ijack];
		  }
		 
		}

		else if(contr=="T") {
		   if((fit_type=="tm") || (fit_type=="comb")) {
		     data[ijack][iens].Br = (T_tm_red).distr_list[iens].distr[ijack];
		     data[ijack][iens].Br_err= (T_tm_red).err(iens);
		     data[ijack][iens].Is_tm = true;
		     data[ijack][iens].a = a_distr_list_red.distr_list[iens].distr[ijack];
		   }

		   if((fit_type=="OS") || (fit_type=="comb")) {
		    data[ijack][iens+off_OS].Br = (T_OS_red).distr_list[iens].distr[ijack];
		    data[ijack][iens+off_OS].Br_err= (T_OS_red).err(iens);
		    data[ijack][iens+off_OS].Is_tm = false;
		    data[ijack][iens+off_OS].a = a_distr_list_red.distr_list[iens].distr[ijack];
		   }
		   
		}

		else if(contr=="L") {
		   if((fit_type=="tm") || (fit_type=="comb")) {
		     data[ijack][iens].Br = (L_tm_red).distr_list[iens].distr[ijack];
		     data[ijack][iens].Br_err= (L_tm_red).err(iens);
		     data[ijack][iens].Is_tm = true;
		     data[ijack][iens].a = a_distr_list_red.distr_list[iens].distr[ijack];
		   }

		   if((fit_type=="OS") || (fit_type=="comb")) {
		    data[ijack][iens+off_OS].Br = (L_OS_red).distr_list[iens].distr[ijack];
		    data[ijack][iens+off_OS].Br_err= (L_OS_red).err(iens);
		    data[ijack][iens+off_OS].Is_tm = false;
		    data[ijack][iens+off_OS].a = a_distr_list_red.distr_list[iens].distr[ijack];
		   }
		   
		}

		

				
		else if(contr=="tot") {
		  if((fit_type=="tm") || (fit_type=="comb")) {
		    data[ijack][iens].Br = (VkVk_tm_red+AkAk_tm_red+A0A0_tm_red+V0V0_tm_red).distr_list[iens].distr[ijack];
		    data[ijack][iens].Br_err= (VkVk_tm_red+AkAk_tm_red+A0A0_tm_red+V0V0_tm_red).err(iens);
		    data[ijack][iens].Is_tm = true;
		    data[ijack][iens].a = a_distr_list_red.distr_list[iens].distr[ijack];
		  }

		  if((fit_type=="OS") || (fit_type=="comb")) {
		    data[ijack][iens+off_OS].Br = (VkVk_OS_red+AkAk_OS_red+A0A0_OS_red+V0V0_OS_red).distr_list[iens].distr[ijack];
		    data[ijack][iens+off_OS].Br_err= (VkVk_OS_red+AkAk_OS_red+A0A0_OS_red+V0V0_OS_red).err(iens);
		    data[ijack][iens+off_OS].Is_tm = false;
		    data[ijack][iens+off_OS].a = a_distr_list_red.distr_list[iens].distr[ijack];
		  }
		 
		}


		else if(contr=="tot_TL") {
		  if((fit_type=="tm") || (fit_type=="comb")) {
		    data[ijack][iens].Br = (T_tm_red+L_tm_red).distr_list[iens].distr[ijack];
		    data[ijack][iens].Br_err= (T_tm_red+L_tm_red).err(iens);
		    data[ijack][iens].Is_tm = true;
		    data[ijack][iens].a = a_distr_list_red.distr_list[iens].distr[ijack];
		  }

		  if((fit_type=="OS") || (fit_type=="comb")) {
		    data[ijack][iens+off_OS].Br = (T_OS_red+L_OS_red).distr_list[iens].distr[ijack];
		    data[ijack][iens+off_OS].Br_err= (T_OS_red+L_OS_red).err(iens);
		    data[ijack][iens+off_OS].Is_tm = false;
		    data[ijack][iens+off_OS].a = a_distr_list_red.distr_list[iens].distr[ijack];
		  }
		 
		}

		
		else crash("contribution: "+contr+" not recognized");
		
		//mean values
		if(ijack==0) {

		  if(contr=="A0A0") {

		   		    
		    if((fit_type=="tm") || (fit_type=="comb")) {
		      data_ch2[ijack][iens].Br = A0A0_tm_red.ave(iens);
		      data_ch2[ijack][iens].Br_err= A0A0_tm_red.err(iens);
		      data_ch2[ijack][iens].Is_tm = true;
		      data_ch2[ijack][iens].a = a_distr_list_red.ave(iens);
		    }
		    
		    if((fit_type=="OS") || (fit_type=="comb")) {
		      data_ch2[ijack][iens+off_OS].Br = A0A0_OS_red.ave(iens);
		      data_ch2[ijack][iens+off_OS].Br_err= A0A0_OS_red.err(iens);
		      data_ch2[ijack][iens+off_OS].Is_tm = false;
		      data_ch2[ijack][iens+off_OS].a = a_distr_list_red.ave(iens);
		    }
		    
		  
		  }

		  else if(contr=="V0V0") {
		     
		   		    
		    if((fit_type=="tm") || (fit_type=="comb")) {
		      data_ch2[ijack][iens].Br = V0V0_tm_red.ave(iens);
		      data_ch2[ijack][iens].Br_err= V0V0_tm_red.err(iens);
		      data_ch2[ijack][iens].Is_tm = true;
		      data_ch2[ijack][iens].a = a_distr_list_red.ave(iens);
		    }
		    
		    if((fit_type=="OS") || (fit_type=="comb")) {
		      data_ch2[ijack][iens+off_OS].Br = V0V0_OS_red.ave(iens);
		      data_ch2[ijack][iens+off_OS].Br_err= V0V0_OS_red.err(iens);
		      data_ch2[ijack][iens+off_OS].Is_tm = false;
		      data_ch2[ijack][iens+off_OS].a = a_distr_list_red.ave(iens);
		    }

		  }
		  
		  
		  else if(contr=="AkAk") {
		    if((fit_type=="tm") || (fit_type=="comb")) {
		      data_ch2[ijack][iens].Br = AkAk_tm_red.ave(iens);
		      data_ch2[ijack][iens].Br_err= AkAk_tm_red.err(iens);
		      data_ch2[ijack][iens].Is_tm = true;
		      data_ch2[ijack][iens].a = a_distr_list_red.ave(iens);
		    }

		    if((fit_type=="OS") || (fit_type=="comb")) {
		      data_ch2[ijack][iens+off_OS].Br = AkAk_OS_red.ave(iens);
		      data_ch2[ijack][iens+off_OS].Br_err= AkAk_OS_red.err(iens);
		      data_ch2[ijack][iens+off_OS].Is_tm = false;
		      data_ch2[ijack][iens+off_OS].a = a_distr_list_red.ave(iens);
		    }
		    
		  }

		  else if(contr=="AX") {
		    if((fit_type=="tm") || (fit_type=="comb")) {
		      data_ch2[ijack][iens].Br = (AkAk_tm_red+A0A0_tm_red).ave(iens);
		      data_ch2[ijack][iens].Br_err= (AkAk_tm_red+A0A0_tm_red).err(iens);
		      data_ch2[ijack][iens].Is_tm = true;
		      data_ch2[ijack][iens].a = a_distr_list_red.ave(iens);
		    }

		    if((fit_type=="OS") || (fit_type=="comb")) {
		      data_ch2[ijack][iens+off_OS].Br = (AkAk_OS_red+A0A0_OS_red).ave(iens);
		      data_ch2[ijack][iens+off_OS].Br_err= (AkAk_OS_red+A0A0_OS_red).err(iens);
		      data_ch2[ijack][iens+off_OS].Is_tm = false;
		      data_ch2[ijack][iens+off_OS].a = a_distr_list_red.ave(iens);
		    }
		    
		  }

		  else if(contr=="VkVk") {
		    if((fit_type=="tm") || (fit_type=="comb")) {
		      data_ch2[ijack][iens].Br = VkVk_tm_red.ave(iens);
		      data_ch2[ijack][iens].Br_err= VkVk_tm_red.err(iens);
		      data_ch2[ijack][iens].Is_tm = true;
		      data_ch2[ijack][iens].a = a_distr_list_red.ave(iens);
		    }

		    if((fit_type=="OS") || (fit_type=="comb")) {
		      data_ch2[ijack][iens+off_OS].Br = VkVk_OS_red.ave(iens);
		      data_ch2[ijack][iens+off_OS].Br_err= VkVk_OS_red.err(iens);
		      data_ch2[ijack][iens+off_OS].Is_tm = false;
		      data_ch2[ijack][iens+off_OS].a = a_distr_list_red.ave(iens);
		    }
		    
		  }

		  else if(contr=="VA") {
		    if((fit_type=="tm") || (fit_type=="comb")) {
		      data_ch2[ijack][iens].Br =  ((VkVk_tm_red+V0V0_tm_red)/(AkAk_tm_red+A0A0_tm_red)).ave(iens);
		      data_ch2[ijack][iens].Br_err=  ((VkVk_tm_red+V0V0_tm_red)/(AkAk_tm_red+A0A0_tm_red)).err(iens);
		      data_ch2[ijack][iens].Is_tm = true;
		      data_ch2[ijack][iens].a = a_distr_list_red.ave(iens);
		    }

		    if((fit_type=="OS") || (fit_type=="comb")) {
		      data_ch2[ijack][iens+off_OS].Br =  ((VkVk_OS_red+V0V0_OS_red)/(AkAk_OS_red+A0A0_OS_red)).ave(iens);
		      data_ch2[ijack][iens+off_OS].Br_err=  ((VkVk_OS_red+V0V0_OS_red)/(AkAk_OS_red+A0A0_OS_red)).err(iens);
		      data_ch2[ijack][iens+off_OS].Is_tm = false;
		      data_ch2[ijack][iens+off_OS].a = a_distr_list_red.ave(iens);
		    }
		    
		  }

		  else if(contr=="VMA") {
		    if((fit_type=="tm") || (fit_type=="comb")) {
		      data_ch2[ijack][iens].Br = ((VkVk_tm_red+V0V0_tm_red-AkAk_tm_red-A0A0_tm_red)/(VkVk_tm_red+AkAk_tm_red+A0A0_tm_red+V0V0_tm_red)).ave(iens);
		      data_ch2[ijack][iens].Br_err=  ((VkVk_tm_red+V0V0_tm_red-AkAk_tm_red-A0A0_tm_red)/(VkVk_tm_red+AkAk_tm_red+A0A0_tm_red+V0V0_tm_red)).err(iens);
		      data_ch2[ijack][iens].Is_tm = true;
		      data_ch2[ijack][iens].a = a_distr_list_red.ave(iens);
		    }
		    
		    if((fit_type=="OS") || (fit_type=="comb")) {
		      data_ch2[ijack][iens+off_OS].Br =  ((VkVk_OS_red+V0V0_OS_red-AkAk_OS_red-A0A0_OS_red)/(VkVk_OS_red+AkAk_OS_red+A0A0_OS_red+V0V0_OS_red)).ave(iens);
		      data_ch2[ijack][iens+off_OS].Br_err=  ((VkVk_OS_red+V0V0_OS_red-AkAk_OS_red-A0A0_OS_red)/(VkVk_OS_red+AkAk_OS_red+A0A0_OS_red+V0V0_OS_red)).err(iens);
		      data_ch2[ijack][iens+off_OS].Is_tm = false;
		      data_ch2[ijack][iens+off_OS].a = a_distr_list_red.ave(iens);
		    }
		    
		  }

		  else if(contr=="T") {
		    if((fit_type=="tm") || (fit_type=="comb")) {
		      data_ch2[ijack][iens].Br = (T_tm_red).ave(iens);
		      data_ch2[ijack][iens].Br_err= (T_tm_red).err(iens);
		      data_ch2[ijack][iens].Is_tm = true;
		      data_ch2[ijack][iens].a = a_distr_list_red.ave(iens);
		    }

		    if((fit_type=="OS") || (fit_type=="comb")) {
		      data_ch2[ijack][iens+off_OS].Br = (T_OS_red).ave(iens);
		      data_ch2[ijack][iens+off_OS].Br_err= (T_OS_red).err(iens);
		      data_ch2[ijack][iens+off_OS].Is_tm = false;
		      data_ch2[ijack][iens+off_OS].a = a_distr_list_red.ave(iens);
		    }
		    
		  }


		  else if(contr=="L") {
		    if((fit_type=="tm") || (fit_type=="comb")) {
		      data_ch2[ijack][iens].Br = (L_tm_red).ave(iens);
		      data_ch2[ijack][iens].Br_err= (L_tm_red).err(iens);
		      data_ch2[ijack][iens].Is_tm = true;
		      data_ch2[ijack][iens].a = a_distr_list_red.ave(iens);
		    }

		    if((fit_type=="OS") || (fit_type=="comb")) {
		      data_ch2[ijack][iens+off_OS].Br = (L_OS_red).ave(iens);
		      data_ch2[ijack][iens+off_OS].Br_err= (L_OS_red).err(iens);
		      data_ch2[ijack][iens+off_OS].Is_tm = false;
		      data_ch2[ijack][iens+off_OS].a = a_distr_list_red.ave(iens);
		    }
		    
		  }
		  
		  else if(contr=="tot") {
		    if((fit_type=="tm") || (fit_type=="comb")) {
		      data_ch2[ijack][iens].Br = (VkVk_tm_red+AkAk_tm_red+A0A0_tm_red+V0V0_tm_red).ave(iens);
		      data_ch2[ijack][iens].Br_err= (VkVk_tm_red+AkAk_tm_red+A0A0_tm_red+V0V0_tm_red).err(iens);
		      data_ch2[ijack][iens].Is_tm = true;
		      data_ch2[ijack][iens].a = a_distr_list_red.ave(iens);
		    }
		    
		    if((fit_type=="OS") || (fit_type=="comb")) {
		      data_ch2[ijack][iens+off_OS].Br = (VkVk_OS_red+AkAk_OS_red+A0A0_OS_red+V0V0_OS_red).ave(iens);
		      data_ch2[ijack][iens+off_OS].Br_err= (VkVk_OS_red+AkAk_OS_red+A0A0_OS_red+V0V0_OS_red).err(iens);
		      data_ch2[ijack][iens+off_OS].Is_tm = false;
		      data_ch2[ijack][iens+off_OS].a = a_distr_list_red.ave(iens);
		    }
		    
		  }

		  else if(contr=="tot_TL") {
		    if((fit_type=="tm") || (fit_type=="comb")) {
		      data_ch2[ijack][iens].Br = (T_tm_red+L_tm_red).ave(iens);
		      data_ch2[ijack][iens].Br_err= (T_tm_red+L_tm_red).err(iens);
		      data_ch2[ijack][iens].Is_tm = true;
		      data_ch2[ijack][iens].a = a_distr_list_red.ave(iens);
		    }
		    
		    if((fit_type=="OS") || (fit_type=="comb")) {
		      data_ch2[ijack][iens+off_OS].Br = (T_OS_red+L_OS_red).ave(iens);
		      data_ch2[ijack][iens+off_OS].Br_err= (T_OS_red+L_OS_red).err(iens);
		      data_ch2[ijack][iens+off_OS].Is_tm = false;
		      data_ch2[ijack][iens+off_OS].a = a_distr_list_red.ave(iens);
		    }
		    
		  }

		  
		  else crash("contribution: "+contr+" not recognized");
			
		}
	      }
	    }

	   
	    //append
	    bf_TAU.Append_to_input_par(data);

	    bf_TAU_ch2.Append_to_input_par(data_ch2);
	    //fit
	    cout<<"Fitting...."<<endl;
	    Bt_fit= bf_TAU.Perform_bootstrap_fit();
	    Bt_fit_ch2= bf_TAU_ch2.Perform_bootstrap_fit();

	   
	    //retrieve parameters
	    distr_t D(UseJack), D2_tm(UseJack), D2_OS(UseJack), D4_tm(UseJack), D4_OS(UseJack);
	    for(int ijack=0;ijack<Njacks;ijack++) {
	      D.distr.push_back( Bt_fit.par[ijack].D);
	      D2_tm.distr.push_back( Bt_fit.par[ijack].D2_tm);
	      D2_OS.distr.push_back( Bt_fit.par[ijack].D2_OS);
	      D4_tm.distr.push_back( Bt_fit.par[ijack].D4_tm);
	      D4_OS.distr.push_back( Bt_fit.par[ijack].D4_OS);
	    }
	    //reduced ch2
	    double ch2= Bt_fit_ch2.get_ch2_ave()/Ndof;


	    //print fit function
	    distr_t_list FF_tm_to_print(UseJack), FF_OS_to_print(UseJack);

	    
	    for(auto &a: a_to_print) {
	      FF_tm_to_print.distr_list.push_back( D + D2_tm*pow(a*QCD_scale,2) + D4_tm*pow(a*QCD_scale,4));
	      FF_OS_to_print.distr_list.push_back( D + D2_OS*pow(a*QCD_scale,2) + D4_OS*pow(a*QCD_scale,4));
	    }
	    string Fit_tag= "../data/tau_decay/"+Tag_reco_type+"/strange/continuum/fit_func/sigma_"+to_string_with_precision(sigma_list_strange[is],3)+"_contr_"+contr+"_ftype_"+fit_type+"_"+poly_type+".dat";
	    Print_To_File({}, {a_to_print, FF_tm_to_print.ave(), FF_tm_to_print.err(), FF_OS_to_print.ave(), FF_OS_to_print.err()},Fit_tag, "", "#a[fm] tm OS, ch2/dof: "+to_string_with_precision(ch2,4));
	    

	    //push back information on ch2  and on fit result

	    res_map.find({contr, fit_type, poly_type})->second.distr_list.push_back(D);
	    ch2_map.find({contr, fit_type, poly_type})->second.push_back(ch2);
	    Nmeas_map.find({contr, fit_type, poly_type})->second.push_back(Nmeas);
	    Ndof_map.find({contr, fit_type, poly_type})->second.push_back(Ndof);
	    


	    }
	    else {

	      //push back fictitious value for a fit that has not been performed

	      res_map.find({contr, fit_type, poly_type})->second.distr_list.push_back(Get_id_jack_distr(Njacks)*0.0);
	      ch2_map.find({contr, fit_type, poly_type})->second.push_back(-1.0);
	      Nmeas_map.find({contr, fit_type, poly_type})->second.push_back(-1.0);
	      Ndof_map.find({contr, fit_type, poly_type})->second.push_back(-1.0);
	    }
	    
	  }
	}
       }

       //print data
       //all
       //tm
       Print_To_File({ls_data_tm_VKVK.Tag}, { a_distr_list.ave(), A0A0_tm.ave(), A0A0_tm.err() } ,  "../data/tau_decay/"+Tag_reco_type+"/strange/continuum/fit_data/A0A0/tm/sigma_"+to_string_with_precision(sigma_list_strange[is],3)+".dat"  , "", "");
       Print_To_File({ls_data_tm_VKVK.Tag}, { a_distr_list.ave(), V0V0_tm.ave(), V0V0_tm.err() } ,  "../data/tau_decay/"+Tag_reco_type+"/strange/continuum/fit_data/V0V0/tm/sigma_"+to_string_with_precision(sigma_list_strange[is],3)+".dat"  , "", "");
       Print_To_File({ls_data_tm_VKVK.Tag}, { a_distr_list.ave(), AkAk_tm.ave(), AkAk_tm.err() } ,  "../data/tau_decay/"+Tag_reco_type+"/strange/continuum/fit_data/AkAk/tm/sigma_"+to_string_with_precision(sigma_list_strange[is],3)+".dat"  , "", "");
       Print_To_File({ls_data_tm_VKVK.Tag}, { a_distr_list.ave(), (AkAk_tm+A0A0_tm).ave(), (AkAk_tm+A0A0_tm).err() } ,  "../data/tau_decay/"+Tag_reco_type+"/strange/continuum/fit_data/AX/tm/sigma_"+to_string_with_precision(sigma_list_strange[is],3)+".dat"  , "", "");
       Print_To_File({ls_data_tm_VKVK.Tag}, { a_distr_list.ave(), (T_tm).ave(), (T_tm).err() } ,  "../data/tau_decay/"+Tag_reco_type+"/strange/continuum/fit_data/T/tm/sigma_"+to_string_with_precision(sigma_list_strange[is],3)+".dat"  , "", "");
       Print_To_File({ls_data_tm_VKVK.Tag}, { a_distr_list.ave(), (L_tm).ave(), (L_tm).err() } ,  "../data/tau_decay/"+Tag_reco_type+"/strange/continuum/fit_data/L/tm/sigma_"+to_string_with_precision(sigma_list_strange[is],3)+".dat"  , "", "");
       
       Print_To_File({ls_data_tm_VKVK.Tag}, { a_distr_list.ave(), VkVk_tm.ave(), VkVk_tm.err() } ,  "../data/tau_decay/"+Tag_reco_type+"/strange/continuum/fit_data/VkVk/tm/sigma_"+to_string_with_precision(sigma_list_strange[is],3)+".dat"  , "", "");
       Print_To_File({ls_data_tm_VKVK.Tag}, { a_distr_list.ave(), (VkVk_tm+AkAk_tm+A0A0_tm+V0V0_tm).ave(), (VkVk_tm+AkAk_tm+A0A0_tm+V0V0_tm).err() } ,  "../data/tau_decay/"+Tag_reco_type+"/strange/continuum/fit_data/tot/tm/sigma_"+to_string_with_precision(sigma_list_strange[is],3)+".dat"  , "", "");
        Print_To_File({ls_data_tm_VKVK.Tag}, { a_distr_list.ave(), (T_tm+L_tm).ave(), (T_tm+L_tm).err() } ,  "../data/tau_decay/"+Tag_reco_type+"/strange/continuum/fit_data/tot_TL/tm/sigma_"+to_string_with_precision(sigma_list_strange[is],3)+".dat"  , "", "");
       
       Print_To_File({ls_data_tm_VKVK.Tag}, { a_distr_list.ave(), (( (VkVk_tm+V0V0_tm)/(AkAk_tm+A0A0_tm))).ave(), (((VkVk_tm+V0V0_tm)/(AkAk_tm+A0A0_tm))).err() } ,  "../data/tau_decay/"+Tag_reco_type+"/strange/continuum/fit_data/VA/tm/sigma_"+to_string_with_precision(sigma_list_strange[is],3)+".dat"  , "", "");
       Print_To_File({ls_data_tm_VKVK.Tag}, { a_distr_list.ave(), (((VkVk_tm+V0V0_tm-AkAk_tm-A0A0_tm)/(VkVk_tm+AkAk_tm+A0A0_tm+V0V0_tm))).ave(), (((VkVk_tm+V0V0_tm-AkAk_tm-A0A0_tm)/(VkVk_tm+AkAk_tm+A0A0_tm+V0V0_tm))).err() } ,  "../data/tau_decay/"+Tag_reco_type+"/strange/continuum/fit_data/VMA/tm/sigma_"+to_string_with_precision(sigma_list_strange[is],3)+".dat"  , "", "");

       
       //OS
       Print_To_File({ls_data_tm_VKVK.Tag}, { a_distr_list.ave(), A0A0_OS.ave(), A0A0_OS.err() } ,  "../data/tau_decay/"+Tag_reco_type+"/strange/continuum/fit_data/A0A0/OS/sigma_"+to_string_with_precision(sigma_list_strange[is],3)+".dat"  , "", "");
       Print_To_File({ls_data_tm_VKVK.Tag}, { a_distr_list.ave(), V0V0_OS.ave(), V0V0_OS.err() } ,  "../data/tau_decay/"+Tag_reco_type+"/strange/continuum/fit_data/V0V0/OS/sigma_"+to_string_with_precision(sigma_list_strange[is],3)+".dat"  , "", "");
       Print_To_File({ls_data_tm_VKVK.Tag}, { a_distr_list.ave(), AkAk_OS.ave(), AkAk_OS.err() } ,  "../data/tau_decay/"+Tag_reco_type+"/strange/continuum/fit_data/AkAk/OS/sigma_"+to_string_with_precision(sigma_list_strange[is],3)+".dat"  , "", "");
       Print_To_File({ls_data_tm_VKVK.Tag}, { a_distr_list.ave(), (AkAk_OS+A0A0_OS).ave(), (AkAk_OS+A0A0_OS).err() } ,  "../data/tau_decay/"+Tag_reco_type+"/strange/continuum/fit_data/AX/OS/sigma_"+to_string_with_precision(sigma_list_strange[is],3)+".dat"  , "", "");
       Print_To_File({ls_data_tm_VKVK.Tag}, { a_distr_list.ave(), (T_OS).ave(), (T_OS).err() } ,  "../data/tau_decay/"+Tag_reco_type+"/strange/continuum/fit_data/T/OS/sigma_"+to_string_with_precision(sigma_list_strange[is],3)+".dat"  , "", "");
       Print_To_File({ls_data_tm_VKVK.Tag}, { a_distr_list.ave(), (L_OS).ave(), (L_OS).err() } ,  "../data/tau_decay/"+Tag_reco_type+"/strange/continuum/fit_data/L/OS/sigma_"+to_string_with_precision(sigma_list_strange[is],3)+".dat"  , "", "");
       
       Print_To_File({ls_data_tm_VKVK.Tag}, { a_distr_list.ave(), VkVk_OS.ave(), VkVk_OS.err() } ,  "../data/tau_decay/"+Tag_reco_type+"/strange/continuum/fit_data/VkVk/OS/sigma_"+to_string_with_precision(sigma_list_strange[is],3)+".dat"  , "", "");
       Print_To_File({ls_data_tm_VKVK.Tag}, { a_distr_list.ave(), (VkVk_OS+AkAk_OS+A0A0_OS+V0V0_OS).ave(), (VkVk_OS+AkAk_OS+A0A0_OS+V0V0_OS).err() } ,  "../data/tau_decay/"+Tag_reco_type+"/strange/continuum/fit_data/tot/OS/sigma_"+to_string_with_precision(sigma_list_strange[is],3)+".dat"  , "", "");
        Print_To_File({ls_data_tm_VKVK.Tag}, { a_distr_list.ave(), (T_OS+L_OS).ave(), (T_OS+L_OS).err() } ,  "../data/tau_decay/"+Tag_reco_type+"/strange/continuum/fit_data/tot_TL/OS/sigma_"+to_string_with_precision(sigma_list_strange[is],3)+".dat"  , "", "");
       
       Print_To_File({ls_data_tm_VKVK.Tag}, { a_distr_list.ave(), (((VkVk_OS+V0V0_OS)/(AkAk_OS+A0A0_OS))).ave(), (((VkVk_OS+V0V0_OS)/(AkAk_OS+A0A0_OS))).err() } ,  "../data/tau_decay/"+Tag_reco_type+"/strange/continuum/fit_data/VA/OS/sigma_"+to_string_with_precision(sigma_list_strange[is],3)+".dat"  , "", "");
       Print_To_File({ls_data_tm_VKVK.Tag}, { a_distr_list.ave(), (((VkVk_OS+V0V0_OS-AkAk_OS-A0A0_OS)/(VkVk_OS+V0V0_OS+AkAk_OS+A0A0_OS))).ave(), (((VkVk_OS+V0V0_OS-AkAk_OS-A0A0_OS)/(VkVk_OS+V0V0_OS+AkAk_OS+A0A0_OS))).err() } ,  "../data/tau_decay/"+Tag_reco_type+"/strange/continuum/fit_data/VMA/OS/sigma_"+to_string_with_precision(sigma_list_strange[is],3)+".dat"  , "", "");

       
       //red
       Print_To_File({Tag_ens_red}, { a_distr_list_red.ave(), A0A0_tm_red.ave(), A0A0_tm_red.err() } ,  "../data/tau_decay/"+Tag_reco_type+"/strange/continuum/fit_data/A0A0/tm/FSE_corr_sigma_"+to_string_with_precision(sigma_list_strange[is],3)+".dat"  , "", "");
       Print_To_File({Tag_ens_red}, { a_distr_list_red.ave(), V0V0_tm_red.ave(), V0V0_tm_red.err() } ,  "../data/tau_decay/"+Tag_reco_type+"/strange/continuum/fit_data/V0V0/tm/FSE_corr_sigma_"+to_string_with_precision(sigma_list_strange[is],3)+".dat"  , "", "");
       Print_To_File({Tag_ens_red}, { a_distr_list_red.ave(), AkAk_tm_red.ave(), AkAk_tm_red.err() } ,  "../data/tau_decay/"+Tag_reco_type+"/strange/continuum/fit_data/AkAk/tm/FSE_corr_sigma_"+to_string_with_precision(sigma_list_strange[is],3)+".dat"  , "", "");
       Print_To_File({Tag_ens_red}, { a_distr_list_red.ave(), (AkAk_tm_red+A0A0_tm_red).ave(), (AkAk_tm_red+A0A0_tm_red).err() } ,  "../data/tau_decay/"+Tag_reco_type+"/strange/continuum/fit_data/AX/tm/FSE_corr_sigma_"+to_string_with_precision(sigma_list_strange[is],3)+".dat"  , "", "");
       Print_To_File({Tag_ens_red}, { a_distr_list_red.ave(), (T_tm_red).ave(), (T_tm_red).err() } ,  "../data/tau_decay/"+Tag_reco_type+"/strange/continuum/fit_data/T/tm/FSE_corr_sigma_"+to_string_with_precision(sigma_list_strange[is],3)+".dat"  , "", "");
        Print_To_File({Tag_ens_red}, { a_distr_list_red.ave(), (L_tm_red).ave(), (L_tm_red).err() } ,  "../data/tau_decay/"+Tag_reco_type+"/strange/continuum/fit_data/L/tm/FSE_corr_sigma_"+to_string_with_precision(sigma_list_strange[is],3)+".dat"  , "", "");
       Print_To_File({Tag_ens_red}, { a_distr_list_red.ave(), VkVk_tm_red.ave(), VkVk_tm_red.err() } ,  "../data/tau_decay/"+Tag_reco_type+"/strange/continuum/fit_data/VkVk/tm/FSE_corr_sigma_"+to_string_with_precision(sigma_list_strange[is],3)+".dat"  , "", "");
       Print_To_File({Tag_ens_red}, { a_distr_list_red.ave(), (VkVk_tm_red+AkAk_tm_red+A0A0_tm_red+V0V0_tm_red).ave(), (VkVk_tm_red+AkAk_tm_red+A0A0_tm_red+V0V0_tm_red).err() } ,  "../data/tau_decay/"+Tag_reco_type+"/strange/continuum/fit_data/tot/tm/FSE_corr_sigma_"+to_string_with_precision(sigma_list_strange[is],3)+".dat"  , "", "");
        Print_To_File({Tag_ens_red}, { a_distr_list_red.ave(), (T_tm_red+L_tm_red).ave(), (T_tm_red+L_tm_red).err() } ,  "../data/tau_decay/"+Tag_reco_type+"/strange/continuum/fit_data/tot_TL/tm/FSE_corr_sigma_"+to_string_with_precision(sigma_list_strange[is],3)+".dat"  , "", "");
       
       Print_To_File({Tag_ens_red}, { a_distr_list_red.ave(), (((VkVk_tm_red+V0V0_tm_red)/(AkAk_tm_red+A0A0_tm_red))).ave(), (( (VkVk_tm_red+V0V0_tm_red)/(AkAk_tm_red+A0A0_tm_red))).err() } ,  "../data/tau_decay/"+Tag_reco_type+"/strange/continuum/fit_data/VA/tm/FSE_corr_sigma_"+to_string_with_precision(sigma_list_strange[is],3)+".dat"  , "", "");
       Print_To_File({Tag_ens_red}, { a_distr_list_red.ave(), (((VkVk_tm_red +V0V0_tm_red -AkAk_tm_red - A0A0_tm_red )/(VkVk_tm_red+ AkAk_tm_red+A0A0_tm_red+V0V0_tm_red))).ave(), (((VkVk_tm_red+V0V0_tm_red -AkAk_tm_red - A0A0_tm_red )/(VkVk_tm_red+V0V0_tm_red+ AkAk_tm_red+A0A0_tm_red))).err() } ,  "../data/tau_decay/"+Tag_reco_type+"/strange/continuum/fit_data/VMA/tm/FSE_corr_sigma_"+to_string_with_precision(sigma_list_strange[is],3)+".dat"  , "", "");
       //OS
       Print_To_File({Tag_ens_red}, { a_distr_list_red.ave(), A0A0_OS_red.ave(), A0A0_OS_red.err() } ,  "../data/tau_decay/"+Tag_reco_type+"/strange/continuum/fit_data/A0A0/OS/FSE_corr_sigma_"+to_string_with_precision(sigma_list_strange[is],3)+".dat"  , "", "");
       Print_To_File({Tag_ens_red}, { a_distr_list_red.ave(), V0V0_OS_red.ave(), V0V0_OS_red.err() } ,  "../data/tau_decay/"+Tag_reco_type+"/strange/continuum/fit_data/V0V0/OS/FSE_corr_sigma_"+to_string_with_precision(sigma_list_strange[is],3)+".dat"  , "", "");
       Print_To_File({Tag_ens_red}, { a_distr_list_red.ave(), AkAk_OS_red.ave(), AkAk_OS_red.err() } ,  "../data/tau_decay/"+Tag_reco_type+"/strange/continuum/fit_data/AkAk/OS/FSE_corr_sigma_"+to_string_with_precision(sigma_list_strange[is],3)+".dat"  , "", "");
       Print_To_File({Tag_ens_red}, { a_distr_list_red.ave(), (AkAk_OS_red+A0A0_OS_red).ave(), (AkAk_OS_red+A0A0_OS_red).err() } ,  "../data/tau_decay/"+Tag_reco_type+"/strange/continuum/fit_data/AX/OS/FSE_corr_sigma_"+to_string_with_precision(sigma_list_strange[is],3)+".dat"  , "", "");
       Print_To_File({Tag_ens_red}, { a_distr_list_red.ave(), (T_OS_red).ave(), (T_OS_red).err() } ,  "../data/tau_decay/"+Tag_reco_type+"/strange/continuum/fit_data/T/OS/FSE_corr_sigma_"+to_string_with_precision(sigma_list_strange[is],3)+".dat"  , "", "");
       Print_To_File({Tag_ens_red}, { a_distr_list_red.ave(), (L_OS_red).ave(), (L_OS_red).err() } ,  "../data/tau_decay/"+Tag_reco_type+"/strange/continuum/fit_data/L/OS/FSE_corr_sigma_"+to_string_with_precision(sigma_list_strange[is],3)+".dat"  , "", "");
       Print_To_File({Tag_ens_red}, { a_distr_list_red.ave(), VkVk_OS_red.ave(), VkVk_OS_red.err() } ,  "../data/tau_decay/"+Tag_reco_type+"/strange/continuum/fit_data/VkVk/OS/FSE_corr_sigma_"+to_string_with_precision(sigma_list_strange[is],3)+".dat"  , "", "");
       Print_To_File({Tag_ens_red}, { a_distr_list_red.ave(), (VkVk_OS_red+AkAk_OS_red+A0A0_OS_red+V0V0_OS_red).ave(), (VkVk_OS_red+AkAk_OS_red+A0A0_OS_red+V0V0_OS_red).err() } ,  "../data/tau_decay/"+Tag_reco_type+"/strange/continuum/fit_data/tot/OS/FSE_corr_sigma_"+to_string_with_precision(sigma_list_strange[is],3)+".dat"  , "", "");
       Print_To_File({Tag_ens_red}, { a_distr_list_red.ave(), (T_OS_red+L_OS_red).ave(), (T_OS_red+L_OS_red).err() } ,  "../data/tau_decay/"+Tag_reco_type+"/strange/continuum/fit_data/tot_TL/OS/FSE_corr_sigma_"+to_string_with_precision(sigma_list_strange[is],3)+".dat"  , "", "");
       
       Print_To_File({Tag_ens_red}, { a_distr_list_red.ave(), (((VkVk_OS_red+V0V0_OS_red)/(AkAk_OS_red+A0A0_OS_red))).ave(), (((VkVk_OS_red+V0V0_OS_red)/(AkAk_OS_red+A0A0_OS_red))).err() } ,  "../data/tau_decay/"+Tag_reco_type+"/strange/continuum/fit_data/VA/OS/FSE_corr_sigma_"+to_string_with_precision(sigma_list_strange[is],3)+".dat"  , "", "");
       Print_To_File({Tag_ens_red}, { a_distr_list_red.ave(), (((VkVk_OS_red + V0V0_OS_red -AkAk_OS_red - A0A0_OS_red )/(VkVk_OS_red+ AkAk_OS_red+A0A0_OS_red+V0V0_OS_red))).ave(), (((VkVk_OS_red +V0V0_OS_red -AkAk_OS_red - A0A0_OS_red )/(VkVk_OS_red+ V0V0_OS_red+ AkAk_OS_red+A0A0_OS_red))).err() } ,  "../data/tau_decay/"+Tag_reco_type+"/strange/continuum/fit_data/VMA/OS/FSE_corr_sigma_"+to_string_with_precision(sigma_list_strange[is],3)+".dat"  , "", "");


       
      
       
    }
    

    //here I should print the result
    for (auto const& [tag,Br] : res_map) {
      string out_tag= get<0>(tag)+"_"+get<1>(tag)+"_"+get<2>(tag);
      Vfloat ch2 = ch2_map.find(tag)->second;
      Print_To_File({}, { sigma_list_strange, Br.ave(), Br.err(), ch2 }, "../data/tau_decay/"+Tag_reco_type+"/strange/continuum/Extr_"+out_tag+".dat", "", "#sigma ave err ch2/dof ");
    }




    //###############################################
    //###############################################
    //###############################################
    //##############                  ###############
    //##############  Akaike analysis ###############
    //##############                  ###############
    //###############################################
    //###############################################
    //###############################################

    Vfloat s_extr_away_TL;

    vector<distr_t_list> Br_finals;
    VVfloat Br_final_systs(Contribs.size());
    for(int c=0; c<(signed)Contribs.size(); c++)  { Br_finals.emplace_back(UseJack, sigma_list_strange.size()); Br_final_systs[c].resize(sigma_list_strange.size(), 0);}

 
    for(int is=0; is < (signed)sigma_list_strange.size();is++) {

      //create directory for AIC output
      boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/continuum/AIC/sigma_"+to_string_with_precision(sigma_list_strange[is],3));
       
      int c=0;
      //loop over contributions
      for(auto &contr: Contribs ) {

	map< tuple<string,string> , double> AIC_map;
	Vfloat weight_list, ch2_list, Nmeas_list, Npars_list;
	vector<string> ftpt_list;
	
	double w_tot=0;
	
	for(auto &ftype: Fit_types)
	  for(auto &ptype: poly_types) {

	    if( ch2_map.find({contr,ftype,ptype})->second[is]  >= 0) {
	    
	    double ch2_i = ch2_map.find({contr,ftype,ptype})->second[is];
	    double Nmeas_i = Nmeas_map.find({contr, ftype,ptype})->second[is];
	    double Ndof_i = Ndof_map.find({contr, ftype, ptype})->second[is];
	    double Npars_i = Nmeas_i - Ndof_i;
	    double w= exp(-0.5*( ch2_i*Ndof_i +2*Npars_i -Nmeas_i));
	    
	    ftpt_list.push_back( ftype+"_"+ptype);
	    weight_list.push_back( w);
	    ch2_list.push_back( ch2_i);
	    Nmeas_list.push_back( Nmeas_i);
	    Npars_list.push_back( Npars_i);
	    
	    AIC_map.insert( { {ftype, ptype}, w });
	    w_tot += w;
	    }
		
	  }
	
	//normalize AIC weights
	for (auto &[tag,w] : AIC_map) w/= w_tot;
	
	distr_t_list Res_partial(UseJack);  
	distr_t Br(UseJack, UseJack?Njacks:Nboots); //constructor automatically sets Br to zero
	for (auto &[tag,w] : AIC_map) {
	  Br = Br + w*res_map.find({contr, get<0>(tag), get<1>(tag)})->second.distr_list[is] ;
	  Res_partial.distr_list.push_back(res_map.find({contr, get<0>(tag), get<1>(tag)})->second.distr_list[is]);
	}
	
	double syst=0;
	double global_ave= Br.ave();
	
	for (auto &[tag,w] : AIC_map) syst += w*pow( res_map.find({contr, get<0>(tag), get<1>(tag)})->second.ave(is) - global_ave,2);
	
	Br_finals[c].distr_list[is]= Br;
	Br_final_systs[c][is] = sqrt(syst);

	if(Contribs[c] == "tot_TL") {
	  s_extr_away_TL.push_back( min( fabs(Br.ave() - tot_TL_Ens_E_tm.distr_list[is].ave() )/sqrt( pow(Br.err(),2) + syst), fabs( Br.ave() -tot_TL_Ens_E_OS.distr_list[is].ave())/sqrt( pow(Br.err(),2) + syst)));
	}
	
	
	  

	
	//print details on AIC
	Print_To_File({ftpt_list}, {weight_list, ch2_list, Nmeas_list, Npars_list, Res_partial.ave(), Res_partial.err()  } , "../data/tau_decay/"+Tag_reco_type+"/strange/continuum/AIC/sigma_"+to_string_with_precision(sigma_list_strange[is],3)+"/"+contr+".out", "", "#ftpt   w(AIC)   ch2/dof  Nmeas  Npars");
	cout<<"generating histograms for "<<Tag_reco_type<<endl<<flush;
	cout<<"Contrib: "<<Contribs[c]<<endl;

	//generate histograms
	int Nboots_histo=5000;
	int k=50;
	Vfloat x_list, s_list;
	double x_min= 0;
	double x_max= 2*Br.ave();
	if( Contribs[c] == "VMA") { x_min= -300*fabs(Br.ave()+ Br.err()); x_max = 300*fabs(Br.ave()+ Br.err()) ;}
	if( Contribs[c] == "V0V0") { x_min = -100*Br.ave(); x_max= 100*Br.ave();}
	

	double sw=1.0201;
	if(Contribs[c] == "VA" || Contribs[c] == "VMA") sw=1.0;
	double s= sqrt( pow(Br.err(),2) + pow( syst, 2))/4;
	double Nsteps= (int)((x_max - x_min)/s);
	cout<<"size of histogram vector: "<<Nsteps<<endl;
	for(int i=0; i< Nsteps;i++) { x_list.push_back( sw*(x_min + (i+0.5)*s)); s_list.push_back(s);}
	Vfloat hist(Nsteps,0);
	GaussianMersenne r(433295);
	for(int iboot=0; iboot<Nboots_histo;iboot++) {
	  for(auto &[tag,w] : AIC_map) {
	    distr_t res= res_map.find({contr, get<0>(tag), get<1>(tag)})->second.distr_list[is];
	    double x = res.ave() + r()*res.err();
	    if( x < x_min  ) crash("x < xmin, x: "+to_string_with_precision(x,5)+" x_min: "+to_string_with_precision(x_min,5));
	    if( x > x_max  ) crash("x > xmax, x: "+to_string_with_precision(x,5)+" x_max: "+to_string_with_precision(x_max,5));
	    hist[(int)( (x-x_min)/s)] += w/Nboots_histo;
	  }
	}

	//print histogram to File
	Print_To_File({}, {x_list, hist, s_list}, "../data/tau_decay/"+Tag_reco_type+"/strange/continuum/AIC/sigma_"+to_string_with_precision(sigma_list_strange[is],3)+"/"+contr+"_histo.out", "", "# x histo[x] ");

	c++;
      }
    }

 
    

    //print results
    for(int c=0;c<(signed)Contribs.size();c++) {
      string contr= Contribs[c];
      if(contr != "tot_TL") {
      Print_To_File({}, { sigma_list_strange, Br_finals[c].ave(), Br_finals[c].err(), Br_final_systs[c]}, "../data/tau_decay/"+Tag_reco_type+"/strange/continuum/Extr_AIC_"+contr+".dat", "", "#sigma ave err_stat  err_syst ");
      }
      else {
	Print_To_File({}, { sigma_list_strange, Br_finals[c].ave(), Br_finals[c].err(), Br_final_systs[c], s_extr_away_TL}, "../data/tau_decay/"+Tag_reco_type+"/strange/continuum/Extr_AIC_"+contr+".dat", "", "#sigma ave err_stat  err_syst s_extr_away_TL ");
      }
      
    }


    //print FSEs from C80-C112 and B64-B96
    Print_To_File({}, {sigma_list_strange, FSE_C_A0A0_tm, FSE_err_C_A0A0_tm, FSE_B_A0A0_tm, FSE_err_B_A0A0_tm}, "../data/tau_decay/"+Tag_reco_type+"/strange/FSE/tm/A0A0.dat", "", "#sigma C   B ");
    Print_To_File({}, {sigma_list_strange, FSE_C_V0V0_tm, FSE_err_C_V0V0_tm,  FSE_B_V0V0_tm, FSE_err_B_V0V0_tm}, "../data/tau_decay/"+Tag_reco_type+"/strange/FSE/tm/V0V0.dat", "", "#sigma C   B ");
    Print_To_File({}, {sigma_list_strange, FSE_C_AkAk_tm, FSE_err_C_AkAk_tm,  FSE_B_AkAk_tm, FSE_err_B_AkAk_tm}, "../data/tau_decay/"+Tag_reco_type+"/strange/FSE/tm/AkAk.dat", "", "#sigma C   B ");
    Print_To_File({}, {sigma_list_strange, FSE_C_VkVk_tm, FSE_err_C_VkVk_tm,  FSE_B_VkVk_tm, FSE_err_B_VkVk_tm}, "../data/tau_decay/"+Tag_reco_type+"/strange/FSE/tm/VkVk.dat", "", "#sigma C   B ");
    Print_To_File({}, {sigma_list_strange, FSE_C_T_tm, FSE_err_C_T_tm,  FSE_B_T_tm, FSE_err_B_T_tm}, "../data/tau_decay/"+Tag_reco_type+"/strange/FSE/tm/T.dat", "", "#sigma C   B ");
    Print_To_File({}, {sigma_list_strange, FSE_C_L_tm, FSE_err_C_L_tm, FSE_B_L_tm, FSE_err_B_L_tm}, "../data/tau_decay/"+Tag_reco_type+"/strange/FSE/tm/L.dat", "", "#sigma C   B ");

    Print_To_File({}, {sigma_list_strange, FSE_C_A0A0_OS, FSE_err_C_A0A0_OS, FSE_B_A0A0_OS, FSE_err_B_A0A0_OS}, "../data/tau_decay/"+Tag_reco_type+"/strange/FSE/OS/A0A0.dat", "", "#sigma C   B ");
    Print_To_File({}, {sigma_list_strange, FSE_C_V0V0_OS, FSE_err_C_V0V0_OS,  FSE_B_V0V0_OS, FSE_err_B_V0V0_OS}, "../data/tau_decay/"+Tag_reco_type+"/strange/FSE/OS/V0V0.dat", "", "#sigma C   B ");
    Print_To_File({}, {sigma_list_strange, FSE_C_AkAk_OS, FSE_err_C_AkAk_OS,  FSE_B_AkAk_OS, FSE_err_B_AkAk_OS}, "../data/tau_decay/"+Tag_reco_type+"/strange/FSE/OS/AkAk.dat", "", "#sigma C   B ");
    Print_To_File({}, {sigma_list_strange, FSE_C_VkVk_OS, FSE_err_C_VkVk_OS,  FSE_B_VkVk_OS, FSE_err_B_VkVk_OS}, "../data/tau_decay/"+Tag_reco_type+"/strange/FSE/OS/VkVk.dat", "", "#sigma C   B ");
    Print_To_File({}, {sigma_list_strange, FSE_C_T_OS, FSE_err_C_T_OS,  FSE_B_T_OS, FSE_err_B_T_OS}, "../data/tau_decay/"+Tag_reco_type+"/strange/FSE/OS/T.dat", "", "#sigma C   B ");
    Print_To_File({}, {sigma_list_strange, FSE_C_L_OS, FSE_err_C_L_OS, FSE_B_L_OS, FSE_err_B_L_OS}, "../data/tau_decay/"+Tag_reco_type+"/strange/FSE/OS/L.dat", "", "#sigma C   B ");


    Print_To_File({}, {sigma_list_strange, pl_T_tm_distr, pl_BIS_T_tm_distr}, "../data/tau_decay/"+Tag_reco_type+"/strange/FSE/tm/pl_T.dat", "", "#sigma C   B ");
    Print_To_File({}, {sigma_list_strange, pl_L_tm_distr, pl_BIS_L_tm_distr}, "../data/tau_decay/"+Tag_reco_type+"/strange/FSE/tm/pl_L.dat", "", "#sigma C   B ");
    Print_To_File({}, {sigma_list_strange, pl_T_OS_distr, pl_BIS_T_OS_distr}, "../data/tau_decay/"+Tag_reco_type+"/strange/FSE/OS/pl_T.dat", "", "#sigma C   B ");
    Print_To_File({}, {sigma_list_strange, pl_L_OS_distr, pl_BIS_L_OS_distr}, "../data/tau_decay/"+Tag_reco_type+"/strange/FSE/OS/pl_L.dat", "", "#sigma C   B ");
    


    //print summary plots

    //get average statistical error
    Vfloat stat_ave_A0A0(sigma_list_strange.size(),0), stat_ave_V0V0(sigma_list_strange.size(),0), stat_ave_AkAk(sigma_list_strange.size(),0), stat_ave_VkVk(sigma_list_strange.size(),0), stat_ave_T(sigma_list_strange.size(),0), stat_ave_L(sigma_list_strange.size(),0);
    Vfloat syst_ave_A0A0(sigma_list_strange.size(),0), syst_ave_V0V0(sigma_list_strange.size(),0), syst_ave_AkAk(sigma_list_strange.size(),0), syst_ave_VkVk(sigma_list_strange.size(),0), syst_ave_T(sigma_list_strange.size(),0), syst_ave_L(sigma_list_strange.size(),0);

    
    int Nen_s=6;
    for(int is=0;is<(signed)sigma_list_strange.size(); is++ ) {
      for(int i=0; i<Nen_s; i++) {
	
	stat_ave_A0A0[is] += (stat_A0A0_tm[i][is]/A0A0_tm_all_s[is].ave(i))/(2.0*Nens) + (stat_A0A0_OS[i][is]/A0A0_OS_all_s[is].ave(i))/(2.0*Nens);
	stat_ave_V0V0[is] += (stat_V0V0_tm[i][is]/V0V0_tm_all_s[is].ave(i))/(2.0*Nens) + (stat_V0V0_OS[i][is]/V0V0_OS_all_s[is].ave(i))/(2.0*Nens);
	stat_ave_VkVk[is] += (stat_VkVk_tm[i][is]/VkVk_tm_all_s[is].ave(i))/(2.0*Nens) + (stat_VkVk_OS[i][is]/VkVk_OS_all_s[is].ave(i))/(2.0*Nens);
	stat_ave_AkAk[is] += (stat_AkAk_tm[i][is]/AkAk_tm_all_s[is].ave(i))/(2.0*Nens) + (stat_AkAk_OS[i][is]/AkAk_OS_all_s[is].ave(i))/(2.0*Nens);
	stat_ave_T[is] += (stat_T_tm[i][is]/T_tm_all_s[is].ave(i))/(2.0*Nens) + (stat_T_OS[i][is]/T_OS_all_s[is].ave(i))/(2.0*Nens);
	stat_ave_L[is] += (stat_L_tm[i][is]/L_tm_all_s[is].ave(i))/(2.0*Nens) +  (stat_L_OS[i][is]/L_tm_all_s[is].ave(i))/(2.0*Nens);

	syst_ave_A0A0[is] += (syst_A0A0_tm[i][is]/A0A0_tm_all_s[is].ave(i))/(2.0*Nens) + (syst_A0A0_OS[i][is]/A0A0_OS_all_s[is].ave(i))/(2.0*Nens);
	syst_ave_V0V0[is] += (syst_V0V0_tm[i][is]/V0V0_tm_all_s[is].ave(i))/(2.0*Nens) + (syst_V0V0_OS[i][is]/V0V0_OS_all_s[is].ave(i))/(2.0*Nens);
	syst_ave_VkVk[is] += (syst_VkVk_tm[i][is]/VkVk_tm_all_s[is].ave(i))/(2.0*Nens) + (syst_VkVk_OS[i][is]/VkVk_OS_all_s[is].ave(i))/(2.0*Nens);
	syst_ave_AkAk[is] += (syst_AkAk_tm[i][is]/AkAk_tm_all_s[is].ave(i))/(2.0*Nens) + (syst_AkAk_OS[i][is]/AkAk_OS_all_s[is].ave(i))/(2.0*Nens);
	syst_ave_T[is] += (syst_T_tm[i][is]/T_tm_all_s[is].ave(i))/(2.0*Nens) + (syst_T_OS[i][is]/T_OS_all_s[is].ave(i))/(2.0*Nens);
	syst_ave_L[is] += (syst_L_tm[i][is]/L_tm_all_s[is].ave(i))/(2.0*Nens) +  (syst_L_OS[i][is]/L_tm_all_s[is].ave(i))/(2.0*Nens);
      }
    }

    Print_To_File({}, {sigma_list_strange, FSE_C_A0A0, stat_ave_A0A0, syst_ave_A0A0}, "../data/tau_decay/"+Tag_reco_type+"/strange/summary/A0A0.dat", "", "#sigma FSE stat syst");
    Print_To_File({}, {sigma_list_strange, FSE_C_V0V0, stat_ave_V0V0, syst_ave_V0V0}, "../data/tau_decay/"+Tag_reco_type+"/strange/summary/V0V0.dat", "", "#sigma FSE stat syst");
    Print_To_File({}, {sigma_list_strange, FSE_C_AkAk, stat_ave_AkAk, syst_ave_AkAk}, "../data/tau_decay/"+Tag_reco_type+"/strange/summary/AkAk.dat", "", "#sigma FSE stat syst");
    Print_To_File({}, {sigma_list_strange, FSE_C_VkVk, stat_ave_VkVk, syst_ave_VkVk}, "../data/tau_decay/"+Tag_reco_type+"/strange/summary/VkVk.dat", "", "#sigma FSE stat syst");
    Print_To_File({}, {sigma_list_strange, FSE_C_T, stat_ave_T, syst_ave_T}, "../data/tau_decay/"+Tag_reco_type+"/strange/summary/T.dat", "", "#sigma FSE stat syst");
    Print_To_File({}, {sigma_list_strange, FSE_C_L, stat_ave_L, syst_ave_L}, "../data/tau_decay/"+Tag_reco_type+"/strange/summary/L.dat", "", "#sigma FSE stat syst");
    
    //###############################################
    //###############################################
    //###############################################


    //perform sigma to zero extrapolation of the various contributions
    boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/strange/sigma_extr");

    vector<bool> Fix_D6_list({true, false});
    vector<int> sigma_to_exclude_list({2,0});
    vector<vector<distr_t>> Final_extr_data(Fix_D6_list.size());

  
    for(int c=0; c<(signed)Contribs.size();c++) {

      

      class ipar_SIGMA {
	
      public:
	ipar_SIGMA()  {}
	
	double Br, Br_err;
	double sigma;
      };
      
      
      class fpar_SIGMA {

      public:
	fpar_SIGMA() {}
	fpar_SIGMA(const Vfloat &par) {
	  if((signed)par.size() != 3) crash("In class fpar_SIGMA, class constructor Vfloat par has size != 3");
	  D=par[0];
	  D4=par[1];
	  D6=par[2];
	}

	double D,D4, D6;
      };


     

      for(int dd=0;dd<(signed)Fix_D6_list.size() ; dd++ ) {
	
	bool Fix_D6= Fix_D6_list[dd];
	int sigma_to_exclude= sigma_to_exclude_list[dd];
	
	bootstrap_fit<fpar_SIGMA,ipar_SIGMA> bf_SIGMA(Njacks);
	bootstrap_fit<fpar_SIGMA,ipar_SIGMA> bf_SIGMA_ch2(1);
	bf_SIGMA.Set_number_of_measurements(sigma_list_strange.size() - sigma_to_exclude);
	bf_SIGMA.Set_verbosity(1);
	//ch2
	bf_SIGMA_ch2.Set_number_of_measurements(sigma_list_strange.size() - sigma_to_exclude);
	bf_SIGMA_ch2.Set_verbosity(1);

	//add fit parameters
	bf_SIGMA.Add_par("D", 2.0, 0.1);
	bf_SIGMA.Add_par("D4", 2, 0.1);
	bf_SIGMA.Add_par("D6", 2, 0.1);
	//ch2
	bf_SIGMA_ch2.Add_par("D", 2.0, 0.1);
	bf_SIGMA_ch2.Add_par("D4", 2, 0.1);
	bf_SIGMA_ch2.Add_par("D6", 2, 0.1);
	
	
	if(Fix_D6) {
	  bf_SIGMA.Fix_par("D6",0.0);
	  bf_SIGMA_ch2.Fix_par("D6",0.0);
	}

	//ansatz
	bf_SIGMA.ansatz=  [ ](const fpar_SIGMA &p, const ipar_SIGMA &ip) {
       
	  return p.D + p.D4*pow(ip.sigma,4) + p.D6*pow(ip.sigma,6);
	};
	//meas
	bf_SIGMA.measurement=  [ ](const fpar_SIGMA &p, const ipar_SIGMA &ip) {
	  return ip.Br;
	};
	//err
	bf_SIGMA.error=  [ ](const fpar_SIGMA &p, const ipar_SIGMA &ip) {
	  return ip.Br_err;
	};
	//ch2
	bf_SIGMA_ch2.ansatz= bf_SIGMA.ansatz;
	bf_SIGMA_ch2.measurement= bf_SIGMA.measurement;
	bf_SIGMA_ch2.error= bf_SIGMA.error;
	
	
	//fill the data
	vector<vector<ipar_SIGMA>> data(Njacks);
	vector<vector<ipar_SIGMA>> data_ch2(1);
	//allocate space for output result
	boot_fit_data<fpar_SIGMA> Bt_fit;
	boot_fit_data<fpar_SIGMA> Bt_fit_ch2;
	for(auto &data_iboot: data) data_iboot.resize(sigma_list_strange.size() - sigma_to_exclude);
	for(auto &data_iboot: data_ch2) data_iboot.resize(sigma_list_strange.size() - sigma_to_exclude);
	
	GaussianMersenne GS(224223); //GS(15431); //;
	
	for(int ijack=0;ijack<Njacks;ijack++) {
	  
	  
	  for(int is=0;is<(signed)sigma_list_strange.size() - sigma_to_exclude;is++) {
	    
	    
	    
	    data[ijack][is].Br = Br_finals[c].ave(is) + (Br_finals[c].distr_list[is].distr[ijack] - Br_finals[c].ave(is))*sqrt( pow(Br_finals[c].err(is),2) + pow(Br_final_systs[c][is],2))/Br_finals[c].err(is);
	    data[ijack][is].Br_err = sqrt( pow(Br_finals[c].err(is),2) + pow(Br_final_systs[c][is],2));
	    data[ijack][is].sigma= sigma_list_strange[is];
	    
	    
	    if(ijack==0) {
	      
	    data_ch2[ijack][is].Br= Br_finals[c].ave(is);
	    data_ch2[ijack][is].Br_err = sqrt( pow(Br_finals[c].err(is),2) + pow(Br_final_systs[c][is],2));
	    data_ch2[ijack][is].sigma = sigma_list_strange[is];
	    
	    }
	    
	  }
	}
	
	
      
	
	//append
	bf_SIGMA.Append_to_input_par(data);
	
	bf_SIGMA_ch2.Append_to_input_par(data_ch2);
	//fit
	cout<<"Fitting...."<<endl;
	Bt_fit= bf_SIGMA.Perform_bootstrap_fit();
	Bt_fit_ch2= bf_SIGMA_ch2.Perform_bootstrap_fit();
	
	
	//retrieve parameters
	distr_t D(UseJack), D4(UseJack), D6(UseJack);
	for(int ijack=0;ijack<Njacks;ijack++) { D.distr.push_back( Bt_fit.par[ijack].D); D4.distr.push_back( Bt_fit.par[ijack].D4); D6.distr.push_back( Bt_fit.par[ijack].D6);}
	//reduced ch2
	int Ndof= sigma_list_strange.size() - sigma_to_exclude - ( (Fix_D6==true)?2:3 );
	double ch2= Bt_fit_ch2.get_ch2_ave()/Ndof;

	Final_extr_data[dd].push_back(D);
	
      
	//print fit function
	distr_t_list R_at_sigma_to_print(UseJack);
	
	Vfloat sigma_to_print;
	for(int iss=0;iss<1000;iss++) { sigma_to_print.push_back( iss*0.0002);}
	
	for(auto &ss: sigma_to_print) {
	  R_at_sigma_to_print.distr_list.push_back( D + D4*pow(ss,4) + D6*pow(ss,6));
	}
	string Fit_tag= "../data/tau_decay/"+Tag_reco_type+"/strange/sigma_extr/contr_"+Contribs[c]+"_ifit_"+to_string(dd);
	Print_To_File({}, {sigma_to_print, R_at_sigma_to_print.ave(), R_at_sigma_to_print.err()},Fit_tag, "", "#sigma Br Br_err "+to_string_with_precision(ch2,4));
	
	
	
      }
    }
      
    
	double Sew= 1.0201;

	
	/*
      
	cout<<"A0 V0  TA    TV   T   A   tot"<<endl;
	cout<<" ~ $"<<to_string_with_precision(Sew*Final_extr_data[0].ave(),3)<<"~("<<to_string_with_precision(1000*Sew*Final_extr_data[0].err(),0)<<")$ ~ &";
	cout<<" ~ $"<<to_string_with_precision(Sew*Final_extr_data[1].ave(),3)<<"~("<<to_string_with_precision(1000*Sew*Final_extr_data[1].err(),0)<<")$ ~ &";
	cout<<" ~ $"<<to_string_with_precision(Sew*Final_extr_data[2].ave(),3)<<"~("<<to_string_with_precision(1000*Sew*Final_extr_data[2].err(),0)<<")$ ~ &";
	cout<<" ~ $"<<to_string_with_precision(Sew*Final_extr_data[3].ave(),3)<<"~("<<to_string_with_precision(1000*Sew*Final_extr_data[3].err(),0)<<")$ ~ &";
	cout<<" ~ $"<<to_string_with_precision(Sew*Final_extr_data[8].ave(),3)<<"~("<<to_string_with_precision(1000*Sew*Final_extr_data[8].err(),0)<<")$ ~ &";
	cout<<" ~ $"<<to_string_with_precision(Sew*Final_extr_data[7].ave(),3)<<"~("<<to_string_with_precision(1000*Sew*Final_extr_data[7].err(),0)<<")$ ~ &";
	cout<<" ~ $"<<to_string_with_precision(Sew*Final_extr_data[4].ave(),3)<<"~("<<to_string_with_precision(1000*Sew*Final_extr_data[4].err(),0)<<")$ ~ \\"<<endl;
      
	*/
	
	distr_t FF= Final_extr_data[0][2].ave() + (Final_extr_data[0][2] - Final_extr_data[0][2].ave())*sqrt( pow(Final_extr_data[0][2].err(),2) + pow( Final_extr_data[0][2].ave()-Final_extr_data[1][2].ave(),2))/Final_extr_data[0][2].err();
	cout<<"R(a,r): "<<Tag_reco_type<<" : "<<Final_extr_data[0][2].ave()<<" +- "<<Final_extr_data[0][2].err()<<" +- "<< Final_extr_data[0][2].ave()-Final_extr_data[1][2].ave()<<endl;

	

	
         RET=Sew*FF;
      
  }
    
    
    
    

      



  


  


  cout<<"Finished calculation of: "<<Tag_reco_type<<endl;
  cout<<"Bye"<<endl;


  
    
 
  


  return RET;

}
