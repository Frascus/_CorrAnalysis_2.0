#include "../include/PionMassAnalysis_twisted_ov_X.h"

using namespace std;
namespace plt = matplotlibcpp;

constexpr double kappa=2.837297;
const double MPiPhys= 0.134977;
const double alpha = 1.0/137.04;
const double e2 = alpha*4.0*M_PI;
const double fpi_phys = 0.1304;
const double r0 = pow(0.659/0.1975,2);
const double f0 = 0.121;
const bool Enable_f0 = 0;
const double csi_phys = pow(MPiPhys,2)/pow(4.0*M_PI*(Enable_f0?f0:fpi_phys),2);
const double xi_phys = pow(fpi_phys*pow(MPiPhys,4),1.0/5);
const int Nbranches = 8;
const int nboots= 100;
const bool Use_JB_distribution= true;
const bool UseJack=1;
const int Njacks=15;
const double Rerr= Use_JB_distribution?sqrt((double)Njacks-1.0):1.0;
const int Nboots=200;
const bool verbose=1;
const bool CDH_correct_FVE=1;
const bool GL_correct_FVE =CDH_correct_FVE?0:0;
const bool Use_fp=1;
const bool use_tailored_interval=true;
const bool Use_A2_Dm_prior=true;
const bool Use_Fa_prior=false;
const double X_phys_val= (1e+3*(Use_fp?pow(fpi_phys,2)/(2.0*MPiPhys):pow(MPiPhys,2)/(2*MPiPhys)));
const double l1ph= -0.4; //-0.4
const double l2ph= 4.3; //4.3
const double l3ph= 3.2; //3.2
const double l4ph= 4.4; //4.4
const double s0= 2.0-M_PI/2.0;
const double s1 = M_PI/4.0 - 0.5;
const double s2 = 0.5 - M_PI/8.0;
const double s3 = 3.0*M_PI/16.0 - 0.5;
bool Use_M3=false;
bool Use_Landau_Gauge=false;
string LorF = Use_Landau_Gauge?"/Landau_gauge":"";

class Y_t_M {

public:
  Y_t_M() {}
  double L, Mpi, Dm_ov_X_val,  Mpi_err, fp, csi, Dm_ov_X_err, f0;
  int ibeta;
};



class X_t_M {

public:
  X_t_M(const Vfloat &par) : ainv(3), Za(3) {
    if(par.size() != 16) {
      cout<<"In class X_t_M, invalid call to constructor"<<endl;
      crash("par_size: +"+to_string(par.size()));
    }
    this->chir=par[0];
    this->log=par[1];
    this->A_1=par[2];
    this->D=par[3];
    this->log_a=par[4];
    this->F_m = par[5];
    this->F_4 = par[6];
    this->F_a=par[7];
    this->Dm=par[8];
    this->A_2=par[9];
    for(int ibeta=0;ibeta<3;ibeta++) {
      this->ainv[ibeta] = par[10+ibeta];
    }
    for(int ibeta=0; ibeta<3;ibeta++) {
      this->Za[ibeta] = par[13+ibeta];
   }
  }
  X_t_M() : ainv(3), Za(3) {}
  vector<double> ainv;
  vector<double> Za;
  double D, chir, log, log_a, A_1, F_a, F_m, F_4, Dm, A_2;
};


void Get_plateaux(string Ensemble_tag, string OBS, int &Tmin, int &Tmax) {

  if(Ensemble_tag=="A30.32_64") {
    if(OBS == "fp") { Tmin=14; Tmax=24;     }
    else if (OBS=="Mpi") { Tmin=19; Tmax=30;        }
    else if (OBS == "Dm") { Tmin=15; Tmax=24;        }
    else if (OBS == "Dm_exch") {   Tmin=15; Tmax=27;               }
    else if (OBS == "Dm_hand") {   Tmin=15; Tmax=27;                }
    else if (OBS == "Dm_exch_untwist") {Tmin=15; Tmax=28;}
    else if (OBS == "Dm_exch_cons") { Tmin= 14; Tmax=29;              }
    else if (OBS == "Dm_exch_Land") {   Tmin=15; Tmax=25;               }
    else if (OBS == "Dm_hand_Land") {   Tmin=15; Tmax=27;                }
    else if (OBS == "Mpi_Z") { Tmin= 13; Tmax=27;}  
    else if (OBS == "Mpi_OS_Z") { Tmin= 15; Tmax=25;}
    else if (OBS == "Zp_ov_Zs") { Tmin= 14; Tmax=30;}
    else if (OBS == "RV") { Tmin=16; Tmax=30;}
    else if (OBS == "RA") { Tmin=9; Tmax=19;}
    else if (OBS == "RV") { Tmin=8; Tmax=31;}
    else if (OBS == "RA") { Tmin=12; Tmax=20;}
    else if (OBS == "RA_bis") {Tmin=13; Tmax=24;}
    else crash("OBS :"+OBS+" not yet implemented");
  }
  else if(Ensemble_tag=="A40.24_48") {
    if(OBS == "fp") { Tmin=12; Tmax=20;         }
    else if (OBS=="Mpi") {  Tmin=12; Tmax=20;       }
    else if (OBS == "Dm") {  Tmin=12; Tmax=19;       }
    else if (OBS == "Dm_exch") {   Tmin=12; Tmax=20;               }
    else if (OBS == "Dm_hand") {   Tmin=10; Tmax=22;                }
    else if (OBS == "Dm_exch_untwist") {Tmin=6; Tmax=19;}
    else if (OBS == "Dm_exch_cons") { Tmin= 13; Tmax=17;              }
    else if (OBS == "Dm_exch_Land") {   Tmin=13; Tmax=20;               }
    else if (OBS == "Dm_hand_Land") {   Tmin=10; Tmax=22;                }
    else if (OBS == "Mpi_Z") { Tmin= 16; Tmax=22;}  
    else if (OBS == "Mpi_OS_Z") { Tmin= 11; Tmax=18;}
    else if (OBS == "Zp_ov_Zs") { Tmin= 11; Tmax=18;}
    else if (OBS == "RV") { Tmin=10; Tmax=22;}
    else if (OBS == "RA") { Tmin=7; Tmax=18;}
    else if (OBS == "RA_bis") {Tmin=5; Tmax=13;}
    else crash("OBS :"+OBS+" not yet implemented");
  }
  else if(Ensemble_tag=="A40.32_64") {
    if(OBS == "fp") {   Tmin=14; Tmax=28;       }
    else if (OBS=="Mpi") { Tmin=14; Tmax=26;        } //20-30
    else if (OBS == "Dm") {  Tmin=14; Tmax=27;       } 
    else if (OBS == "Dm_exch") {   Tmin=14; Tmax=27;               } //14-30
    else if (OBS == "Dm_hand") {   Tmin=14; Tmax=27;                } //14-27
    else if (OBS == "Dm_exch_untwist") {Tmin=16; Tmax=28;}
    else if (OBS == "Dm_exch_cons") { Tmin= 13; Tmax=26;              }
    else if (OBS == "Dm_exch_Land") {   Tmin=14; Tmax=28;               } //14-30
    else if (OBS == "Dm_hand_Land") {   Tmin=14; Tmax=27;                } //14-27
    else if (OBS == "Mpi_Z") { Tmin= 15; Tmax=24;}  
    else if (OBS == "Mpi_OS_Z") { Tmin= 17; Tmax=27;}
    else if (OBS == "Zp_ov_Zs") { Tmin= 16; Tmax=29;}
    else if (OBS == "RV") { Tmin=10; Tmax=28;}
    else if (OBS == "RA") { Tmin=10; Tmax=25;}
    else if (OBS == "RA_bis") {Tmin=5; Tmax=12;}

    else crash("OBS :"+OBS+" not yet implemented");
  }
   else if(Ensemble_tag=="A40.40_80") {
     if(OBS == "fp") {  Tmin=18; Tmax=30;        }
     else if (OBS=="Mpi") {  Tmin=18; Tmax=30;       }
     else if (OBS == "Dm") {  Tmin=14; Tmax=30;       }
     else if (OBS == "Dm_exch") {   Tmin=21; Tmax=35;               } //14-30
     else if (OBS == "Dm_hand") {   Tmin=11; Tmax=35;                } //11-32
     else if (OBS == "Dm_exch_untwist") {Tmin=21; Tmax=29;}
     else if (OBS == "Dm_exch_cons") { Tmin= 15; Tmax=30;              }
     else if (OBS == "Dm_exch_Land") {   Tmin=15; Tmax=32;               } //14-30
     else if (OBS == "Dm_hand_Land") {   Tmin=21; Tmax=36;                } //11-32
     else if (OBS == "Mpi_Z") { Tmin= 21; Tmax=37;}  
     else if (OBS == "Mpi_OS_Z") { Tmin= 15; Tmax=24;}
     else if (OBS == "Zp_ov_Zs") { Tmin= 15; Tmax=28;}
     else if (OBS == "RV") { Tmin=11; Tmax=38;}
     else if (OBS == "RA") { Tmin=13; Tmax=23;}
     else if (OBS == "RA_bis") {Tmin=6; Tmax=12;}
    else crash("OBS :"+OBS+" not yet implemented");
  }
   else if(Ensemble_tag=="A40.48_96") {
     if(OBS == "fp") {   Tmin=15; Tmax=35;       }
     else if (OBS=="Mpi") {  Tmin=16; Tmax=39;       }
     else if (OBS == "Dm") {  Tmin=15; Tmax=42;       }                //15-42
     else if (OBS == "Dm_exch") {   Tmin=17; Tmax=40;               }  //17-40
     else if (OBS == "Dm_hand") {   Tmin=14; Tmax=40;                } //17-42
     else if (OBS == "Dm_exch_untwist") {Tmin=22; Tmax=40;}
     else if (OBS == "Dm_exch_cons") { Tmin= 16; Tmax=40;              }
     else if (OBS == "Dm_exch_Land") {   Tmin=20; Tmax=40;               }  //17-40
     else if (OBS == "Dm_hand_Land") {   Tmin=24; Tmax=39;                } //17-42
     else if (OBS == "Mpi_Z") { Tmin= 25; Tmax=40;}  
     else if (OBS == "Mpi_OS_Z") { Tmin= 12; Tmax=20;}
     else if (OBS == "Mpi_OS_Z") { Tmin= 12; Tmax=17;} //not so good
     else if (OBS == "Zp_ov_Zs") { Tmin= 10; Tmax=21;}
     else if (OBS == "RV") { Tmin=10; Tmax=43;}
     else if (OBS == "RA") { Tmin=11; Tmax=23;}
     else if (OBS == "RA_bis") {Tmin=6; Tmax=11;}
    else crash("OBS :"+OBS+" not yet implemented");
  }
   else if(Ensemble_tag=="A50.32_64") {
     if(OBS == "fp") {  Tmin=15; Tmax=30;        }
     else if (OBS=="Mpi") {  Tmin=15;Tmax=30;       }
     else if (OBS == "Dm") { Tmin=16; Tmax=28;        }
     else if (OBS == "Dm_exch") {   Tmin=17; Tmax=26;               }
     else if (OBS == "Dm_hand") {   Tmin=16; Tmax=30;                }
     else if (OBS == "Dm_exch_untwist") {Tmin=14; Tmax=29;}
     else if (OBS == "Dm_exch_cons") { Tmin= 15; Tmax=27;              }
     else if (OBS == "Dm_exch_Land") {   Tmin=12; Tmax=28;               }
     else if (OBS == "Dm_hand_Land") {   Tmin=16; Tmax=30;                }
     else if (OBS == "Mpi_Z") { Tmin= 16; Tmax=28;}  
     else if (OBS == "Mpi_OS_Z") { Tmin= 12; Tmax=17;} //not so good
     else if (OBS == "Zp_ov_Zs") { Tmin= 9; Tmax=18;}
     else if (OBS == "RV") { Tmin=10; Tmax=31;}
     else if (OBS == "RA") { Tmin=4; Tmax=11;}
     else if (OBS == "RA_bis") {Tmin=13; Tmax=22;}
    else crash("OBS :"+OBS+" not yet implemented");
  }
   else if(Ensemble_tag=="A60.24_48") {
     if(OBS == "fp") {  Tmin=13; Tmax=22;        }
     else if (OBS=="Mpi") { Tmin=12; Tmax=21;         }
     else if (OBS == "Dm") {  Tmin=11;Tmax=22;       }
     else if (OBS == "Dm_exch") {   Tmin=12; Tmax=22;               }
     else if (OBS == "Dm_hand") {   Tmin=13; Tmax=21;                }
     else if (OBS == "Dm_exch_untwist") {Tmin=15; Tmax=19;}
     else if (OBS == "Dm_exch_cons") { Tmin= 12; Tmax=22;              }
     else if (OBS == "Dm_exch_Land") {   Tmin=12; Tmax=18;               }
     else if (OBS == "Dm_hand_Land") {   Tmin=13; Tmax=21;                }
     else if (OBS == "Mpi_Z") { Tmin= 16; Tmax=22;}  //not so good
     else if (OBS == "Mpi_OS_Z") { Tmin= 14; Tmax=19;} //not so good
     else if (OBS == "Zp_ov_Zs") { Tmin= 13; Tmax=21;}
     else if (OBS == "RV") { Tmin=10; Tmax=23;}
     else if (OBS == "RA") { Tmin=5; Tmax=14;}
     else if (OBS == "RA_bis") {Tmin=9; Tmax=20;}
    else crash("OBS :"+OBS+" not yet implemented");
  }
   else if(Ensemble_tag=="A80.24_48") {
     if(OBS == "fp") {    Tmin=12; Tmax=23;      }
     else if (OBS=="Mpi") { Tmin=12;Tmax=22;        }
     else if (OBS == "Dm") { Tmin=14; Tmax=20;        }
     else if (OBS == "Dm_exch") {   Tmin=16; Tmax=23;               }
     else if (OBS == "Dm_hand") {   Tmin=13; Tmax=22;                }
     else if (OBS == "Dm_exch_untwist") {Tmin=13; Tmax=20;}
     else if (OBS == "Dm_exch_cons") { Tmin= 12; Tmax=20;              }
     else if (OBS == "Dm_exch_Land") {   Tmin=13; Tmax=22;               }
     else if (OBS == "Dm_hand_Land") {   Tmin=11; Tmax=22;                }
     else if (OBS == "Mpi_Z") { Tmin= 13; Tmax=21;}
     else if (OBS == "Mpi_OS_Z") { Tmin= 15; Tmax=18;}  //terrible
     else if (OBS == "Zp_ov_Zs") { Tmin= 14; Tmax=20;}
     else if (OBS == "RV") { Tmin=10; Tmax=23;}
     else if (OBS == "RA") { Tmin=8; Tmax=19;}
     else if (OBS == "RA_bis") {Tmin=15; Tmax=19;}
    else crash("OBS :"+OBS+" not yet implemented");
  }
   else if(Ensemble_tag=="A100.24_48") {
     if(OBS == "fp") {   Tmin=15; Tmax=22;       }
     else if (OBS=="Mpi") {   Tmin=15; Tmax=22;      }
     else if (OBS == "Dm") {  Tmin=12; Tmax=22;       }
     else if (OBS == "Dm_exch") {   Tmin=13; Tmax=22;               }
     else if (OBS == "Dm_hand") {   Tmin=11; Tmax=22;                }
     else if (OBS == "Dm_exch_untwist") {Tmin=10; Tmax=18;}
     else if (OBS == "Dm_exch_cons") { Tmin= 11; Tmax=22;              }
     else if (OBS == "Dm_exch_Land") {   Tmin=12; Tmax=22;               }
     else if (OBS == "Dm_hand_Land") {   Tmin=14; Tmax=22;                }
     else if (OBS == "Mpi_Z") { Tmin= 11; Tmax=19;}
     else if (OBS == "Mpi_OS_Z") { Tmin= 12; Tmax=18;}
     else if (OBS == "Zp_ov_Zs") { Tmin= 12; Tmax=20;}
     else if (OBS == "RV") { Tmin=12; Tmax=22;}
     else if (OBS == "RA") { Tmin=8; Tmax=18;}
     else if (OBS == "RA_bis") {Tmin=12; Tmax=18;}
    else crash("OBS :"+OBS+" not yet implemented");
  }
   else if(Ensemble_tag=="B25.32_64") {
     if(OBS == "fp") {   Tmin=17; Tmax=30;       }
     else if (OBS=="Mpi") {  Tmin=16; Tmax=26;       }
     else if (OBS == "Dm") {  Tmin=10; Tmax=28;       }
     else if (OBS == "Dm_exch") {   Tmin=13; Tmax=25;               }
     else if (OBS == "Dm_hand") {   Tmin=15; Tmax=29;                }
    else crash("OBS :"+OBS+" not yet implemented for Ensemble: "+Ensemble_tag);
  }
   else if(Ensemble_tag=="B35.32_64") {
     if(OBS == "fp") {    Tmin=17; Tmax=30;      }
     else if (OBS=="Mpi") {  Tmin=21; Tmax=28;       }
     else if (OBS == "Dm") {  Tmin=15; Tmax=28;       } //also Tmin=15 ?
     else if (OBS == "Dm_exch") {   Tmin=17; Tmax=28;               }
     else if (OBS == "Dm_hand") {   Tmin=14; Tmax=29;                }
     else if (OBS == "Dm_exch_untwist") {Tmin=11; Tmax=26;}
     else if (OBS == "Dm_exch_cons") { Tmin= 15; Tmax=25;              }
     else if (OBS == "Dm_exch_Land") {   Tmin=16; Tmax=29;               }
     else if (OBS == "Dm_hand_Land") {   Tmin=14; Tmax=29;                }
     else if (OBS == "Mpi_Z") { Tmin= 19; Tmax=29;}
     else if (OBS == "Mpi_OS_Z") { Tmin= 14; Tmax=24;}
     else if (OBS == "Zp_ov_Zs") { Tmin= 15; Tmax=24;}
     else if (OBS == "RV") { Tmin=10; Tmax=30;}
     else if (OBS == "RA") { Tmin=8; Tmax=26;}
     else if (OBS == "RA_bis") {Tmin=5; Tmax=26;}
    else crash("OBS :"+OBS+" not yet implemented");
  }
   else if(Ensemble_tag=="B55.32_64") {
     if(OBS == "fp") {  Tmin=17; Tmax=30;        }
     else if (OBS=="Mpi") {  Tmin=15; Tmax=30;       }
     else if (OBS == "Dm") { Tmin=15; Tmax=24;        }
     else if (OBS == "Dm_exch") {   Tmin=15; Tmax=24;               }
     else if (OBS == "Dm_hand") {   Tmin=12; Tmax=29;                }
     else if (OBS == "Dm_exch_untwist") {Tmin=15; Tmax=29;}
     else if (OBS == "Dm_exch_cons") { Tmin= 18; Tmax=28;              }
     else if (OBS == "Dm_exch_Land") {   Tmin=13; Tmax=27;               }
     else if (OBS == "Dm_hand_Land") {   Tmin=15; Tmax=30;                }
     else if (OBS == "Mpi_Z") { Tmin= 14; Tmax=27;}
     else if (OBS == "Mpi_OS_Z") { Tmin= 13; Tmax=22;}
     else if (OBS == "Zp_ov_Zs") { Tmin= 12; Tmax=23;}
     else if (OBS == "RV") { Tmin=11; Tmax=27;}
     else if (OBS == "RA") { Tmin=13; Tmax=22;}
     else if (OBS == "RA_bis") {Tmin=13; Tmax=21;}
    else crash("OBS :"+OBS+" not yet implemented");
  }
   else if(Ensemble_tag=="B75.32_64") {
     if(OBS == "fp") {  Tmin=14; Tmax=30;        }
     else if (OBS=="Mpi") { Tmin=18; Tmax=27;        }
     else if (OBS == "Dm") {  Tmin=19; Tmax=30;       }
     else if (OBS == "Dm_exch") {   Tmin=19; Tmax=30;               }
     else if (OBS == "Dm_hand") {   Tmin=12; Tmax=30;                }
     else if (OBS == "Dm_exch_untwist") {Tmin=11; Tmax=20;}
     else if (OBS == "Dm_exch_cons") { Tmin= 18; Tmax=30;              }
     else if (OBS == "Dm_exch_Land") {   Tmin=18; Tmax=30;               }
     else if (OBS == "Dm_hand_Land") {   Tmin=12; Tmax=25;                }
     else if (OBS == "Mpi_Z") { Tmin= 22; Tmax=27;}
     else if (OBS == "Mpi_OS_Z") { Tmin= 17; Tmax=26;}
     else if (OBS == "Zp_ov_Zs") { Tmin= 14; Tmax=27;}
     else if (OBS == "RV") { Tmin=10; Tmax=30;}
     else if (OBS == "RA") { Tmin=11; Tmax=30;}
     else if (OBS == "RA_bis") {Tmin=11; Tmax=25;}
    else crash("OBS :"+OBS+" not yet implemented");
  }
   else if(Ensemble_tag=="B85.24_48") {
     if(OBS == "fp") {  Tmin=14; Tmax=22;        }
     else if (OBS=="Mpi") {  Tmin=15; Tmax=22;       }
     else if (OBS == "Dm") {  Tmin=11; Tmax=21;       }
     else if (OBS == "Dm_exch") {   Tmin=12; Tmax=21;               }
     else if (OBS == "Dm_hand") {   Tmin=12; Tmax=22;                }
     else if (OBS == "Dm_exch_untwist") {Tmin=13; Tmax=17;}
     else if (OBS == "Dm_exch_cons") { Tmin= 14; Tmax=21;              }
     else if (OBS == "Dm_exch_Land") {   Tmin=15; Tmax=22;               }
     else if (OBS == "Dm_hand_Land") {   Tmin=12; Tmax=22;                }
     else if (OBS == "Mpi_Z") { Tmin= 15; Tmax=22;}
     else if (OBS == "Mpi_OS_Z") { Tmin= 16; Tmax=21;}
     else if (OBS == "Zp_ov_Zs") { Tmin= 10; Tmax=22;}
     else if (OBS == "RV") { Tmin=11; Tmax=38;}
     else if (OBS == "RA") { Tmin=11; Tmax=16;}
     else if (OBS == "RA_bis") {Tmin=11; Tmax=19;}
    else crash("OBS :"+OBS+" not yet implemented");
  }
   else if(Ensemble_tag=="D15.48_96") {
     if(OBS == "fp") { Tmin=20; Tmax=40;         }
     else if (OBS=="Mpi") { Tmin=21; Tmax=42;        }
     else if (OBS == "Dm") { Tmin=21; Tmax=40;        }
     else if (OBS == "Dm_exch") {   Tmin=26; Tmax=40;               }
     else if (OBS == "Dm_hand") {   Tmin=19; Tmax=42;                }
     else if (OBS == "Dm_exch_untwist") {Tmin=23; Tmax=44;}
     else if (OBS == "Dm_exch_cons") { Tmin= 17; Tmax=43;              }
     else if (OBS == "Dm_exch_Land") {   Tmin=26; Tmax=40;               }
     else if (OBS == "Dm_hand_Land") {   Tmin=19; Tmax=42;                }
     else if (OBS == "Mpi_Z") { Tmin= 22; Tmax=40;}
     else if (OBS == "Mpi_OS_Z") { Tmin= 22; Tmax=43;}
     else if (OBS == "Zp_ov_Zs") { Tmin= 24; Tmax=48;}
     else if (OBS == "RV") { Tmin=10; Tmax=47;}
     else if (OBS == "RA") { Tmin=10; Tmax=36;}
     else if (OBS == "RA_bis") {Tmin=10; Tmax=36;}
    else crash("OBS :"+OBS+" not yet implemented");
  }
   else if(Ensemble_tag=="D20.48_96") {
     if(OBS == "fp") {  Tmin=20; Tmax=40;        }
     else if (OBS=="Mpi") { Tmin=24; Tmax=36;        }
     else if (OBS == "Dm") {  Tmin=24; Tmax=38;       }
     else if (OBS == "Dm_exch") {   Tmin=26; Tmax=38;               }   
     else if (OBS == "Dm_hand") {   Tmin=17; Tmax=40;                }
     else if (OBS == "Dm_exch_untwist") {Tmin=20; Tmax=39;}
     else if (OBS == "Dm_exch_cons") { Tmin= 28; Tmax=39;              }
     else if (OBS == "Dm_exch_Land") {   Tmin=20; Tmax=40;               }   
     else if (OBS == "Dm_hand_Land") {   Tmin=15; Tmax=43;                }
     else if (OBS == "Mpi_Z") { Tmin= 24; Tmax=38;}
     else if (OBS == "Mpi_OS_Z") { Tmin= 25; Tmax=41;}
     else if (OBS == "Zp_ov_Zs") { Tmin= 24; Tmax=48;}
     else if (OBS == "RV") { Tmin=10; Tmax=47;}
     else if (OBS == "RA") { Tmin=7; Tmax=28;}
     else if (OBS == "RA_bis") {Tmin=19; Tmax=28;}
    else crash("OBS :"+OBS+" not yet implemented");
  }
   else if(Ensemble_tag=="D30.48_96") {
     if(OBS == "fp") { Tmin=23; Tmax=42;         }
     else if (OBS=="Mpi") { Tmin=23; Tmax=42;        }
     else if (OBS == "Dm") {  Tmin=20; Tmax=35;       }
     else if (OBS == "Dm_exch") {   Tmin=20; Tmax=35;               } 
     else if (OBS == "Dm_hand") {   Tmin=21; Tmax=44;                }
     else if (OBS == "Dm_exch_untwist") { Tmin=20; Tmax=40;}
     else if (OBS == "Dm_exch_cons") { Tmin= 20; Tmax=32;              }
     else if (OBS == "Dm_exch_Land") {   Tmin=20; Tmax=39;               } 
     else if (OBS == "Dm_hand_Land") {   Tmin=15; Tmax=44;                }
     else if (OBS == "Mpi_Z") { Tmin= 21; Tmax=45;}
     else if (OBS == "Mpi_OS_Z") { Tmin= 23; Tmax=32;}
     else if (OBS == "Zp_ov_Zs") { Tmin= 20; Tmax=34;}
     else if (OBS == "RV") { Tmin=9; Tmax=47;}
     else if (OBS == "RA") { Tmin=11; Tmax=31;}
     else if (OBS == "RA_bis") {Tmin=21; Tmax=30;}
     else crash("OBS :"+OBS+" not yet implemented");
  }
   else crash("Ensemble: "+Ensemble_tag+" not found in Get_plateaux");
  

  return;

}



void Pion_mass_analysis_twisted_ov_X(string CURRENT_TYPE, bool IncludeDisconnected) {

  omp_set_num_threads(1);
  
  data_t m_data, dm_exch_data, dm_hand_data, m_data_hand_run, m_twisted_data, m_data_conserved, dm_exch_untwist_data, dm_exch_cons_data, m_data_c, m_data_l;
  if (CURRENT_TYPE=="LOCAL") { //current is local
    // m_data.Read("../datasets", "mes_contr_00", "P5P5");
    //dm_exch_data.Read("../datasets", "mes_contr_LL", "P5P5");
    m_data.Read("../datasets_RTM_QED/data", "mes_contr_M0_R0_0_M0_R0_0", "P5P5");
    m_data_conserved.Read("../datasets_RTM_QED/data", "mes_contr_00_conserved", "P5P5");
    m_data_c.Read("../datasets/data", "mes_contr_00", "P5P5");
    m_data_l.Read("../datasetslocal/data", "mes_contr_M0_R0_0_M0_R0_0", "P5P5");
    //PAY ATTENTION TO WHAT 2pt YOU ARE LOADING
    dm_exch_data.Read("../datasets_RTM_QED/data", "mes_contr_M0_R0_F_M0_R0_F", "P5P5");
    dm_exch_untwist_data.Read( "../datasetslocal/data", "mes_contr_M0_R0_F_M0_R0_F", "P5P5"  );
    dm_exch_cons_data.Read("../datasets/data", "mes_contr_LL", "P5P5");
    if(IncludeDisconnected) {
      m_data_hand_run.Read("../datasets_RTM_QED/data", "mes_contr_M0_R0_0_M0_R0_0_handcuffs", "P5P5");        //PAY ATTENTION TO WHAT 2pt YOU ARE LOADING
      dm_hand_data.Read("../datasets_RTM_QED/data", "handcuffs", "P5P5");
    }
  }
  else crash("CURRENT_TYPE: "+CURRENT_TYPE+ " not yet implemented. Exiting...");

  
  //#######################################################################################################
  //load data from Landau gauge
  data_t hand_Land,  exch_Land_tw, pt2_exch_Land, pt2_hand_Land;



  hand_Land.Read("../Gauge_fix_m100_PionQed_bis/data/DISCO_LANDAU", "handcuffs", "P5P5");
  exch_Land_tw.Read("../Gauge_fix_m100_PionQed_bis/data/CONN_LANDAU", "mes_contr_QED_6", "P5P5");
  pt2_exch_Land.Read("../Gauge_fix_m100_PionQed_bis/data/CONN_LANDAU", "mes_contr_LO_1", "P5P5");
  pt2_hand_Land.Read("../Gauge_fix_m100_PionQed_bis/data/DISCO_LANDAU", "mes_contr_hand_LO", "P5P5");
  
  //#######################################################################################################


  //#######################################################################################################
  //load data from galileo100
  data_t hand_Land_Gal, hand_Land_Gal_new_norm;



  hand_Land_Gal.Read("../test_disc_galileo/data", "handcuffs", "P5P5");
  hand_Land_Gal_new_norm.Read("../test_disc_galileo/data_new_norm", "handcuffs", "P5P5");
  
  

  //#######################################################################################################
  
  
  
  //to compute ZV
  data_t corr_A0P5, corr_P5P5;
  //to compute ZA
  data_t corr_A0P5_OS, corr_P5P5_OS;


  corr_A0P5.Read("../ZaForPiQed/data", "mes_contr_M0_R0_0_M0_R0_0", "A0P5");
  corr_P5P5.Read("../ZaForPiQed/data", "mes_contr_M0_R0_0_M0_R0_0", "P5P5");
  corr_A0P5_OS.Read("../ZaForPiQed/data", "mes_contr_M0_R0_0_M0_R1_0", "A0P5");
  corr_P5P5_OS.Read("../ZaForPiQed/data", "mes_contr_M0_R0_0_M0_R1_0", "P5P5");

  
  
  int Nens = m_data.size; 
  cout<<"N_ens: "<<Nens<<endl;



  distr_t_list Mpi_distr_list(UseJack), Dm_ov_X_distr_list(UseJack), X_pi_distr_list(UseJack), Dm2_ren_univ_sub_distr_list(UseJack),  Dm_ov_X_ren_univ_sub_distr_list(UseJack),
    Dm2_ren_univ_sub_M1_distr_list(UseJack), diff_Dm2_ren_univ_sub_M1_M2_distr_list(UseJack);
  distr_t_list  ratio_exch_tw_untw(UseJack), ratio_exch_loc_cons(UseJack), ratio_exch_tw_cons(UseJack);
  distr_t_list Mpi_distr_dim_list(UseJack);
  distr_t_list Dm2_disc(UseJack);
  distr_t_list Dm2_disc_Landau(UseJack);
  distr_t_list fp_fit_distr_list(UseJack);
  distr_t_list fp_fit_naive_list(UseJack);
  distr_t_list fp_fit_GL_list(UseJack);
  distr_t_list fp_fit_CDH_list(UseJack);
  distr_t_list Mpi_fit_naive_list(UseJack);
  distr_t_list Mpi_fit_GL_list(UseJack);
  distr_t_list Mpi_fit_CDH_list(UseJack);
  distr_t_list Ratio_CDH_naive_fp(UseJack);
  distr_t_list Ratio_CDH_naive_Mpi(UseJack);
  vector<Eigen::MatrixXd> CovMatrixPion(Nens);
  distr_t_list ratio_LF_exch(Nens), ratio_LF_hand(Nens), ratio_LF_tot(Nens);
  distr_t_list ratio_LF_exch_Univ_sub(Nens), ratio_LF_tot_Univ_sub(Nens);

  Vfloat a_list;

  for(int i=0; i < Nens; i++) {

    if(dm_exch_data.Tag[i] != dm_exch_untwist_data.Tag[i]) crash("twisted and untwisted ens do not coincide");

    if(dm_exch_data.Tag[i] != dm_exch_cons_data.Tag[i]) crash("local twisted and conserved ens do not coincide");

    if(dm_exch_data.Tag[i] != corr_P5P5.Tag[i]) crash("local twisted and RCs ens do not coincide");
   
    CorrAnalysis Corr(UseJack, Njacks,Nboots);

 
    double ll, mm, tt;
    Read_pars_from_ensemble_tag(m_data.Tag[i], mm, ll, tt);

    if(m_data.Tag[i].substr(0,1) == "A") {
      if((int)ll==20) Corr.Tmin=13;
      else Corr.Tmin = 14;
    }
    else if(m_data.Tag[i].substr(0,1) =="B") Corr.Tmin = 15;
    else Corr.Tmin= 21;
    if(m_data.Tag[i].substr(0,1)=="D") Corr.Tmax = 36;
    else {
      if(m_data.Tag[i].substr(0,1) == "A") {
	if((int)ll == 20) Corr.Tmax= 17; 
	else if((int)ll == 24) Corr.Tmax= 20;
	else if((int)ll == 32) Corr.Tmax= 25;
	else if((int)ll == 40) Corr.Tmax= 30;
	else if((int)ll == 48) Corr.Tmax= 35;
	else crash("volume does not exist V: "+to_string_with_precision(ll, 4));
      }
      else if(m_data.Tag[i].substr(0,1) == "B") {
	if((int)ll == 24) Corr.Tmax=20;
	else Corr.Tmax = m_data.nrows[i]/2 -6;

      }
      else crash("lattice spacing not found when fixing time intervals");
    }
    Corr.Nt = m_data.nrows[i];
    boost::filesystem::create_directory("../data");
    boost::filesystem::create_directory("../data/Mpi_twisted_ov_X");
    string p= "local";
    boost::filesystem::create_directory("../data/Mpi_twisted_ov_X/"+p);
    boost::filesystem::create_directory("../data/Mpi_twisted_ov_X/"+p+LorF+"/");

   
    //read Za from RI-MOM
    LatticeInfo LL("LOCAL");
    LL.LatInfo(m_data.Tag[i].substr(0,1));
    double Za = LL.Za_M2;
    double Za_e = LL.Za_M2_err;
    double Za_M1 = LL.Za;
    double Za_M1_e= LL.Za_err;
    double Zv = LL.Zv_M2;
    double Zv_e = LL.Zv_M2_err;
    double Zv_M1 = LL.Zv;
    double Zv_M1_e = LL.Zv_err;
    double lat_spacing_inv = LL.ainv;
    a_list.push_back(1.0/lat_spacing_inv);
    distr_t ainv_distr(UseJack);
    GaussianMersenne G(444363);
    distr_t Za_distr(UseJack);
    distr_t Za_distr_M1(UseJack);
    distr_t Zv_distr(UseJack);
    distr_t Zv_distr_M1(UseJack);
    for(int nj=0;nj<Njacks;nj++) Za_distr_M1.distr.push_back( Za_M1 + G()*Za_M1_e/sqrt(Njacks-1));
    for(int nj=0;nj<Njacks;nj++) Za_distr.distr.push_back( Za + G()*Za_e/sqrt(Njacks-1));
    for(int nj=0;nj<Njacks;nj++) Zv_distr_M1.distr.push_back( Zv_M1 + 0.0*G()*Zv_M1_e/sqrt(Njacks-1));
    for(int nj=0;nj<Njacks;nj++) Zv_distr.distr.push_back( Zv + G()*Zv_e/sqrt(Njacks-1));


    //generate fake Jack-distribution of inverse lattice spacing
    for(int nj=0; nj<Njacks;nj++) ainv_distr.distr.push_back( LL.ainv  + G()*LL.ainv_err/sqrt(Njacks -1.0));
    

    
    distr_t_list Pi_iso_correlated_distr = Corr.corr_t(m_data.col(0)[i], "");
   
    distr_t_list Pi_iso_distr = Corr.corr_t(m_data_conserved.col(0)[i], "");

    distr_t_list Pi_iso_c = Corr.corr_t(m_data_c.col(0)[i], "");

    distr_t_list Pi_iso_l = Corr.corr_t(m_data_l.col(0)[i], "");

    distr_t_list Mpi_eff_distr = Corr.effective_mass_t(Pi_iso_distr,  "../data/Mpi_twisted_ov_X/"+p+"/mass."+m_data.Tag[i]);

    distr_t_list Mpi_eff_distr_corr= Corr.effective_mass_t(Pi_iso_correlated_distr, "../data/Mpi_twisted_ov_X/"+p+"/mass_corr."+m_data.Tag[i]);

    
    if(use_tailored_interval) Get_plateaux(m_data.Tag[i], "Mpi",Corr.Tmin,Corr.Tmax);

    
    distr_t_list Exch_distr =  Corr.corr_t(dm_exch_data.col(0)[i], "");

 
    distr_t_list Exch_untwist_distr = Corr.corr_t(dm_exch_untwist_data.col(0)[i], "");

    distr_t_list Exch_cons_distr = Corr.corr_t(dm_exch_cons_data.col(0)[i], "");

   
    

    distr_t_list fp_distr= Corr.decay_constant_t(pow(2.0*mm,2)*Pi_iso_distr, "../data/Mpi_twisted_ov_X/"+p+"/fp."+m_data.Tag[i]);

    distr_t_list fp_distr_corr= Corr.decay_constant_t(pow(2.0*mm,2)*Pi_iso_correlated_distr, "../data/Mpi_twisted_ov_X/"+p+"/fp_corr."+m_data.Tag[i]);

    

   



    if(use_tailored_interval) Get_plateaux(m_data.Tag[i], "fp",Corr.Tmin,Corr.Tmax);
  
    distr_t fp_fit= Corr.Fit_distr(fp_distr);
    distr_t fp_fit_naive = fp_fit;
    distr_t fp_fit_naive_corr = Corr.Fit_distr(fp_distr_corr);
    distr_t fp_fit_corr = fp_fit_naive_corr;

    if(use_tailored_interval) Get_plateaux(m_data.Tag[i], "Mpi",Corr.Tmin,Corr.Tmax);
  
    distr_t Mpi_fit = Corr.Fit_distr(Mpi_eff_distr);
    distr_t Mpi_fit_naive= Mpi_fit;
    distr_t Mpi_fit_naive_corr = Corr.Fit_distr(Mpi_eff_distr_corr);
    distr_t Mpi_fit_corr = Mpi_fit_naive_corr;

    auto F_15 = [](double x) { return pow(x,1.0/5);};
    auto LOG = [](double x) { return log(x);};

    distr_t X_pi_fit = distr_t::f_of_distr(F_15, fp_fit*Mpi_fit*Mpi_fit*Mpi_fit*Mpi_fit);


    //############# CORRECT Fpi and Mpi   Usi G&L and CDH formulae ################
    distr_t fp_fit_GL, fp_fit_CDH, Mpi_fit_GL, Mpi_fit_CDH, fp_fit_corr_GL, fp_fit_corr_CDH, Mpi_fit_corr_GL, Mpi_fit_corr_CDH;
    distr_t csi_L = Mpi_fit*Mpi_fit/(pow(4.0*M_PI,2)*fp_fit*fp_fit);
    distr_t g1 = distr_t::f_of_distr(g1_l, Mpi_fit*ll);
    distr_t g2 = distr_t::f_of_distr(g2_l, Mpi_fit*ll);
    distr_t log_l = log(csi_phys) - distr_t::f_of_distr(LOG, csi_L);

    distr_t csi_L_C = Mpi_fit_corr*Mpi_fit_corr/(pow(4.0*M_PI,2)*fp_fit_corr*fp_fit_corr);
    distr_t g1_C = distr_t::f_of_distr(g1_l, Mpi_fit_corr*ll);
    distr_t g2_C = distr_t::f_of_distr(g2_l, Mpi_fit_corr*ll);
    distr_t log_l_C = log(csi_phys) - distr_t::f_of_distr(LOG, csi_L_C);

    fp_fit_GL = fp_fit/(1.0-1.2*2.0*csi_L*g1);
    Mpi_fit_GL = Mpi_fit/(1.0 + 1.2*0.5*csi_L*g1);
    fp_fit_CDH = fp_fit/(1.0 -2.0*csi_L*g1 +2.0*csi_L*csi_L*( (Cf1(l1ph,l2ph,l3ph,l4ph) + Sf1(s0,s1,s2,s3) + Cf1_log()*log_l)*g1 + (Cf2(l1ph,l2ph,l3ph,l4ph) + Sf2(s0,s1,s2,s3) + Cf2_log()*log_l)*g2));
    Mpi_fit_CDH = Mpi_fit/(1.0 + 0.5*csi_L*g1 - csi_L*csi_L*( (Cm1(l1ph,l2ph,l3ph,l4ph) + Sm1(s0,s1,s2,s3) + Cm1_log()*log_l)*g1 + (Cm2(l1ph,l2ph,l3ph,l4ph) + Sm2(s0,s1,s2,s3) + Cm2_log()*log_l)*g2));

    fp_fit_corr_GL = fp_fit_corr/(1.0-1.2*2.0*csi_L_C*g1_C);
    Mpi_fit_corr_GL = Mpi_fit_corr/(1.0 + 1.2*0.5*csi_L_C*g1_C);
    fp_fit_corr_CDH = fp_fit_corr/(1.0 -2.0*csi_L_C*g1_C +2.0*csi_L_C*csi_L_C*( (Cf1(l1ph,l2ph,l3ph,l4ph) + Sf1(s0,s1,s2,s3) + Cf1_log()*log_l_C)*g1_C + (Cf2(l1ph,l2ph,l3ph,l4ph) + Sf2(s0,s1,s2,s3) + Cf2_log()*log_l_C)*g2_C));
    Mpi_fit_corr_CDH = Mpi_fit_corr/(1.0 + 0.5*csi_L_C*g1_C - csi_L_C*csi_L_C*( (Cm1(l1ph,l2ph,l3ph,l4ph) + Sm1(s0,s1,s2,s3) + Cm1_log()*log_l_C)*g1_C + (Cm2(l1ph,l2ph,l3ph,l4ph) + Sm2(s0,s1,s2,s3) + Cm2_log()*log_l_C)*g2_C));


    //#################### STORE NAIVE, GL and CDH corrected data######################
    Mpi_fit_naive_list.distr_list.push_back(Mpi_fit);
    Mpi_fit_GL_list.distr_list.push_back(Mpi_fit_GL);
    Mpi_fit_CDH_list.distr_list.push_back(Mpi_fit_CDH);
    fp_fit_naive_list.distr_list.push_back(fp_fit);
    fp_fit_GL_list.distr_list.push_back(fp_fit_GL);
    fp_fit_CDH_list.distr_list.push_back(fp_fit_CDH);
    Ratio_CDH_naive_Mpi.distr_list.push_back(Mpi_fit_CDH/Mpi_fit);
    Ratio_CDH_naive_fp.distr_list.push_back(fp_fit_CDH/fp_fit);

    //############ CHOOSE WHICH ONE TO USE IN THE FIT ############################

    if(GL_correct_FVE || CDH_correct_FVE) {
      if(CDH_correct_FVE) { fp_fit = fp_fit_CDH; Mpi_fit = Mpi_fit_CDH; fp_fit_corr= fp_fit_corr_CDH; Mpi_fit_corr= Mpi_fit_corr_CDH;}
      else { fp_fit = fp_fit_GL; Mpi_fit = Mpi_fit_GL; fp_fit_corr= fp_fit_GL; Mpi_fit_corr = Mpi_fit_corr_GL;} 
    }

  
    

    //############ END FVE CORRECTIONS ON Fpi and Mpi #############################



    //########################## COMPUTE ZA AND ZV USING HADRONIC DEFINITION  (method M3)   #######################################

    //define lambdas
    auto sqr= [=](double a, double b) {return sqrt(a);};
    auto SINH= [](double m) -> double  {return sinh(m);};
    auto SINH2 = [](double m, double t) -> double {return sinh(m);};
  


    distr_t_list P5P5_distr, overlap_P5P5_distr, overlap_P5P5_OS_distr;
    distr_t_list A0P5_distr, A0P5_OS_distr;
    distr_t_list P5P5_OS_distr,  RV, RA, ratio_P5P5_overlap_OS_tm, Zp_ov_Zs_distr, RA_bis, Zp_ov_Zs_distr_bis;
    distr_t_list Mpi_OS_distr_Z, Mpi_distr_Z;
    distr_t Zv_M3, Za_M3, Zp_ov_Zs, Mpi_OS_fit_Z, Mpi_fit_Z, Za_M4;

    //tm pion sector
    P5P5_distr= Corr.corr_t(corr_P5P5.col(0)[i], "");

    Get_plateaux(m_data.Tag[i], "Mpi_Z", Corr.Tmin, Corr.Tmax);
    Mpi_distr_Z= Corr.effective_mass_t(P5P5_distr, "../data/Mpi_twisted_ov_X/"+p+"/Mpi_Z."+m_data.Tag[i]);
    Mpi_fit_Z=Corr.Fit_distr(Mpi_distr_Z);
    overlap_P5P5_distr = Corr.residue_t(P5P5_distr, "");

    //OS sector
    P5P5_OS_distr = Corr.corr_t(corr_P5P5_OS.col(0)[i], "");
    Get_plateaux(m_data.Tag[i], "Mpi_OS_Z", Corr.Tmin, Corr.Tmax);
    Mpi_OS_distr_Z= Corr.effective_mass_t(P5P5_OS_distr, "../data/Mpi_twisted_ov_X/"+p+"/Mpi_OS_Z."+m_data.Tag[i]);
    Mpi_OS_fit_Z= Corr.Fit_distr(Mpi_OS_distr_Z);
    overlap_P5P5_OS_distr= Corr.residue_t(P5P5_OS_distr, "");

    //take ratio between OS and tm pion amplitude to compute Zp/Zs RC.
    ratio_P5P5_overlap_OS_tm= overlap_P5P5_OS_distr/overlap_P5P5_distr;
    Zp_ov_Zs_distr = distr_t_list::f_of_distr_list(sqr, ratio_P5P5_overlap_OS_tm);
    Zp_ov_Zs_distr_bis = distr_t_list::f_of_distr_list(sqr, Corr.residue_t_wo_m_fit(P5P5_OS_distr,"")/Corr.residue_t_wo_m_fit(P5P5_distr, ""));


    //antysymmetrize A0P5 correlators
    Corr.Reflection_sign = -1;
    A0P5_distr= Corr.corr_t(corr_A0P5.col(0)[i], "");
    A0P5_OS_distr = Corr.corr_t(corr_A0P5_OS.col(0)[i], "");
    //restore symmetrization
    Corr.Reflection_sign = 1;

  
    Get_plateaux(m_data.Tag[i], "Zp_ov_Zs", Corr.Tmin, Corr.Tmax);
    Zp_ov_Zs = Corr.Fit_distr(Zp_ov_Zs_distr);
    Get_plateaux(m_data.Tag[i], "RV", Corr.Tmin, Corr.Tmax);
    RV= 2.0*mm*P5P5_distr/distr_t_list::derivative(A0P5_distr, 0);
    Zv_M3= Corr.Fit_distr(RV);
    Get_plateaux(m_data.Tag[i], "RA", Corr.Tmin, Corr.Tmax);
    RA = 2.0*mm*(P5P5_OS_distr/distr_t_list::derivative(A0P5_OS_distr, 0))*(Mpi_OS_fit_Z/Mpi_fit_Z)*(distr_t::f_of_distr(SINH, Mpi_OS_fit_Z)/distr_t::f_of_distr(SINH, Mpi_fit_Z))*(1.0/Zp_ov_Zs);
    Za_M3 = Corr.Fit_distr(RA);
    RA_bis = 2.0*mm*(P5P5_OS_distr/distr_t_list::derivative(A0P5_OS_distr, 0))*(Mpi_OS_distr_Z/Mpi_distr_Z)*(distr_t_list::f_of_distr_list(SINH2, Mpi_OS_distr_Z)/distr_t_list::f_of_distr_list(SINH2, Mpi_distr_Z))*(1.0/Zp_ov_Zs_distr_bis);
    Get_plateaux(m_data.Tag[i], "RA_bis", Corr.Tmin, Corr.Tmax);
    Za_M4 = Corr.Fit_distr(RA_bis);
    


    //print Ra, Rv, Zp_ov_Zs

    Print_To_File({}, {RA_bis.ave(), RA_bis.err()}, "../data/Mpi_twisted_ov_X/"+p+"/RA."+m_data.Tag[i]+".t", "", "");
    Print_To_File({}, {RV.ave(), RV.err()}, "../data/Mpi_twisted_ov_X/"+p+"/RV."+m_data.Tag[i]+".t", "", "");
    Print_To_File({}, {Zp_ov_Zs_distr.ave(), Zp_ov_Zs_distr.err()}, "../data/Mpi_twisted_ov_X/"+p+"/Zp_ov_Zs."+m_data.Tag[i]+".t", "", "");

    

    //############################# COMPUTATION IN LANDAU GAUGE ################################################
    //#########################################################################################################//
    //LOAD LANDAU GAUGE
    //load 2pts
    if(pt2_exch_Land.Tag[i] != m_data.Tag[i] || hand_Land.Tag[i] != m_data.Tag[i] || exch_Land_tw.Tag[i] != m_data.Tag[i] || pt2_hand_Land.Tag[i] != m_data.Tag[i]) crash("Ordering of ensembles in Landau and Feynman Gauge is not the same");
    distr_t_list Land_2pt_exch = Corr.corr_t(pt2_exch_Land.col(0)[i], "");
    distr_t_list Land_2pt_hand = Corr.corr_t(pt2_hand_Land.col(0)[i], "");
    distr_t_list Land_Mpi_exch = Corr.effective_mass_t(Land_2pt_exch, "../data/Mpi_twisted_ov_X/"+p+"/Landau_gauge/mass_Landau."+m_data.Tag[i]);
          
    //load exchange diagrams

    distr_t_list Exch_Land_tw = Corr.corr_t(exch_Land_tw.col(0)[i], "");


    
    distr_t_list Hand_Land(UseJack); 
    //compute effective slopes
    distr_t_list eff_slope_exch_Land_tw = Corr.effective_slope_t(Exch_Land_tw, Land_2pt_exch, "../data/Mpi_twisted_ov_X/"+p+"/Landau_gauge/dm_exch_Landau_tw."+m_data.Tag[i]);
    distr_t_list eff_slope_hand_Land;
    if(IncludeDisconnected) {
      Hand_Land= Corr.corr_t(hand_Land.col(0)[i], "../data/Mpi_twisted_ov_X/"+p+"/Landau_gauge/corr_hand_Landau_tw."+m_data.Tag[i]);
      eff_slope_hand_Land= Corr.effective_slope_t(Hand_Land, Land_2pt_exch, "../data/Mpi_twisted_ov_X/"+p+"/Landau_gauge/dm_hand_Landau_tw."+m_data.Tag[i]);
      if(use_tailored_interval) Get_plateaux(m_data.Tag[i], "Dm_hand_Land",Corr.Tmin,Corr.Tmax);
      Dm2_disc_Landau.distr_list.push_back(Za_distr*Za_distr*e2*Mpi_fit*LL.ainv*LL.ainv*Corr.Fit_distr(eff_slope_hand_Land));
      distr_t_list Dm2_hand_eff_distr_Landau= eff_slope_hand_Land*2.0*Mpi_eff_distr;
      Print_To_File({}, {Dm2_hand_eff_distr_Landau.ave(),Dm2_hand_eff_distr_Landau.err()}, "../data/Mpi_twisted_ov_X/"+p+"/Landau_gauge/dm2_hand."+m_data.Tag[i], "", "");
    }

    //##############################################################################################################
    
          
     


    //############################################################################################################################

   

    
    distr_t_list Dm_exch_eff_distr = Corr.effective_slope_t(Exch_distr, Pi_iso_correlated_distr, "../data/Mpi_twisted_ov_X/"+p+"/dm_exch."+m_data.Tag[i]);
    distr_t_list Dm_exch_eff_distr_uncorr= Corr.effective_slope_t(Exch_distr, Pi_iso_distr, "../data/Mpi_twisted_ov_X/"+p+"/dm_exch_uncorr."+m_data.Tag[i]);
    distr_t_list Dm_exch_untwist_eff_distr = Corr.effective_slope_t(Exch_untwist_distr, Pi_iso_l, "../data/Mpi_twisted_ov_X/"+p+"/dm_exch_untwist."+m_data.Tag[i]);
    distr_t_list Dm_exch_cons_eff_distr = Corr.effective_slope_t(Exch_cons_distr, Pi_iso_c, "../data/Mpi_twisted_ov_X/"+p+"/dm_exch_cons."+m_data.Tag[i]);

    

    distr_t_list Dm_exch_eff_distr_ren= Corr.effective_slope_t(Za_distr*Za_distr*Exch_distr, Pi_iso_correlated_distr, "../data/Mpi_twisted_ov_X/"+p+"/dm_exch_ren."+m_data.Tag[i]);
    distr_t_list Dm_exch_untwist_eff_distr_ren = Corr.effective_slope_t(Zv_M3*Zv_M3*Exch_untwist_distr, Pi_iso_l, "../data/Mpi_twisted_ov_X/"+p+"/dm_exch_untwist_ren."+m_data.Tag[i]);




    Get_plateaux(m_data.Tag[i], "Dm_exch", Corr.Tmin, Corr.Tmax);
    distr_t Dm_exch_fit_Feyn = Corr.Fit_distr(Dm_exch_eff_distr_ren) +(kappa/(4.0*M_PI))*(1.0/ll)*( 1.0 + 2.0/(Mpi_fit*ll)) ;
    Get_plateaux(m_data.Tag[i], "Dm_exch_untwist", Corr.Tmin, Corr.Tmax);
    distr_t Dm_exch_untwist_fit_Feyn = Corr.Fit_distr(Dm_exch_untwist_eff_distr_ren)  +(kappa/(4.0*M_PI))*(1.0/ll)*( 1.0 + 2.0/(Mpi_fit*ll));
    Get_plateaux(m_data.Tag[i], "Dm_exch_cons", Corr.Tmin, Corr.Tmax);
    distr_t Dm_cons_fit = Corr.Fit_distr(Dm_exch_cons_eff_distr)  +(kappa/(4.0*M_PI))*(1.0/ll)*( 1.0 + 2.0/(Mpi_fit*ll));
   
  
    distr_t Dm_hand_fit(UseJack);
    distr_t_list Dm_hand_eff_distr(UseJack), Mpi_eff_hand_run_distr(UseJack), Hand_distr(UseJack), Pi_iso_distr_hand_run(UseJack);
    if(IncludeDisconnected) {
      Get_plateaux(m_data.Tag[i], "Dm_hand", Corr.Tmin, Corr.Tmax);
      Hand_distr = Corr.corr_t(dm_hand_data.col(0)[i], "");
      distr_t_list Mpi_hand_run = Corr.corr_t(m_data_hand_run.col(0)[i], "");
      distr_t_list ratio_Hand_Pi = Hand_distr/Pi_iso_correlated_distr;
      distr_t_list ratio_Hand_Pi_hand_run = Hand_distr/Mpi_hand_run;
      Print_To_File({}, {ratio_Hand_Pi.ave(), ratio_Hand_Pi.err()} ,"../data/Mpi_twisted_ov_X/"+p+"/ratio_Hand_pi."+m_data.Tag[i]+".dat", "", "");
      Print_To_File({}, {ratio_Hand_Pi_hand_run.ave(), ratio_Hand_Pi_hand_run.err()} ,"../data/Mpi_twisted_ov_X/"+p+"/ratio_Hand_pi_hand_run."+m_data.Tag[i]+".dat", "", "");
      Dm_hand_eff_distr= Corr.effective_slope_t(Hand_distr, Pi_iso_correlated_distr, "../data/Mpi_twisted_ov_X/"+p+"/dm_hand."+m_data.Tag[i]);   
      Dm_hand_fit = Corr.Fit_distr(Dm_hand_eff_distr);
      distr_t_list Dm2_hand_eff_distr = Dm_hand_eff_distr*2.0*Mpi_eff_distr;
      Print_To_File({}, {Dm2_hand_eff_distr.ave(),Dm2_hand_eff_distr.err()}, "../data/Mpi_twisted_ov_X/"+p+"/dm2_hand."+m_data.Tag[i], "", "");
      Dm2_disc.distr_list.push_back(Za_distr*Za_distr*e2*Mpi_fit*LL.ainv*LL.ainv*Dm_hand_fit);

      //print to file the difference between LANDAU and FEYNMAN hand eff slope
      Print_To_File({}, {(Dm_hand_eff_distr-eff_slope_hand_Land ).ave(), (Dm_hand_eff_distr-eff_slope_hand_Land).err()},"../data/Mpi_twisted_ov_X/"+p+"/Landau_gauge/dm_hand_LF_diff."+m_data.Tag[i], "", "");
    }
    

    
   
    
    Mpi_distr_list.distr_list.push_back(Mpi_fit);
    Mpi_distr_dim_list.distr_list.push_back( ainv_distr*Mpi_fit);
    fp_fit_distr_list.distr_list.push_back(fp_fit);
    
    distr_t_list Dm_tot_eff_distr=Dm_exch_eff_distr;
    distr_t_list Dm_tot_eff_distr_uncorr(UseJack);

    
  
    distr_t_list Dm_tot_eff_distr_renorm= Corr.effective_slope_t(Za_distr*Za_distr*Exch_distr , Pi_iso_correlated_distr, "../data/Mpi_twisted_ov_X/"+p+"/dm_exch_renorm_univ_sub."+m_data.Tag[i]);
    if(IncludeDisconnected) {
      Dm_tot_eff_distr = Corr.effective_slope_t(Exch_distr-Hand_distr, Pi_iso_correlated_distr, "../data/Mpi_twisted_ov_X/"+p+"/dm_tot_corr."+m_data.Tag[i]);
      Dm_tot_eff_distr_uncorr = Corr.effective_slope_t(Exch_distr- Hand_distr, Pi_iso_distr, "../data/Mpi_twisted_ov_X/"+p+"/dm_tot."+m_data.Tag[i]);
      distr_t_list Dm2_tot_distr_corr= Corr.effective_slope_t( Mpi_fit_naive_corr*(Exch_distr-Hand_distr), Pi_iso_correlated_distr, "../data/Mpi_twisted_ov_X/"+p+"/dm2_tot_corr."+m_data.Tag[i]);
      distr_t_list Dm2_tot_distr = Corr.effective_slope_t( Mpi_fit_naive*(Exch_distr-Hand_distr), Pi_iso_distr, "../data/Mpi_twisted_ov_X/"+p+"/dm2_tot."+m_data.Tag[i]);
      Dm_tot_eff_distr_renorm= Corr.effective_slope_t(Za_distr*Za_distr*(Exch_distr-Hand_distr) , Pi_iso_correlated_distr, "../data/Mpi_twisted_ov_X/"+p+"/dm_tot_renorm_univ_sub."+m_data.Tag[i]);
    }


   
    distr_t Dm_tot_fit_Feyn(UseJack), Dm_tot_fit_exch_Feyn(UseJack), Dm_tot_fit_hand_Feyn(UseJack), Dm2_ov_X_Feyn(UseJack);
    distr_t Dm_tot_fit_Land(UseJack), Dm_tot_fit_exch_Land(UseJack), Dm_tot_fit_hand_Land(UseJack), Dm2_ov_X_Land(UseJack);
    distr_t Dm_tot_fit(UseJack);
     
    //FEYNMAN FIT
    if(use_tailored_interval) Get_plateaux(m_data.Tag[i], "Dm_exch",Corr.Tmin,Corr.Tmax);
    Dm_tot_fit_exch_Feyn= Corr.Fit_distr(Dm_exch_eff_distr);
    Dm_tot_fit_Feyn=Dm_tot_fit_exch_Feyn;
    //LANDAU FIT
    if(use_tailored_interval) Get_plateaux(m_data.Tag[i], "Dm_exch_Land",Corr.Tmin,Corr.Tmax);
    Dm_tot_fit_exch_Land = Corr.Fit_distr(eff_slope_exch_Land_tw);
    Dm_tot_fit_Land = Dm_tot_fit_exch_Land;
    if(IncludeDisconnected) {
      //FEYNMAN
      if(use_tailored_interval) Get_plateaux(m_data.Tag[i], "Dm_hand",Corr.Tmin,Corr.Tmax);
      Dm_tot_fit_hand_Feyn= Corr.Fit_distr(Dm_hand_eff_distr);
      Dm_tot_fit_Feyn = Dm_tot_fit_exch_Feyn- Dm_tot_fit_hand_Feyn;    
      //LANDAU
      if(use_tailored_interval) Get_plateaux(m_data.Tag[i], "Dm_hand_Land",Corr.Tmin,Corr.Tmax);
      Dm_tot_fit_hand_Land = Corr.Fit_distr(eff_slope_hand_Land);
      Dm_tot_fit_Land = Dm_tot_fit_exch_Land- Dm_tot_fit_hand_Land;
     
    }

    ratio_LF_exch.distr_list.push_back( Dm_tot_fit_exch_Land/Dm_tot_fit_exch_Feyn);
    ratio_LF_hand.distr_list.push_back( Dm_tot_fit_hand_Land/Dm_tot_fit_hand_Feyn);
    ratio_LF_tot.distr_list.push_back( Dm_tot_fit_Land/Dm_tot_fit_Feyn);
    distr_t FSE_fact = (1.0/Mpi_fit)*(kappa/(2.0*M_PI))*(1.0/ll)*(Mpi_fit+ 2.0/ll);
    ratio_LF_exch_Univ_sub.distr_list.push_back( (Za_distr*Za_distr*Dm_tot_fit_exch_Land + FSE_fact)/(Za_distr*Za_distr*Dm_tot_fit_exch_Feyn + FSE_fact));
    ratio_LF_tot_Univ_sub.distr_list.push_back( (Za_distr*Za_distr*Dm_tot_fit_Land + FSE_fact)/(Za_distr*Za_distr*Dm_tot_fit_Feyn + FSE_fact));

    Dm2_ov_X_Feyn = 2.0*Dm_tot_fit_Feyn/(Use_fp?(fp_fit*fp_fit/Mpi_fit):Mpi_fit);
    Dm2_ov_X_Land = 2.0*Dm_tot_fit_Land/(Use_fp?(fp_fit*fp_fit/Mpi_fit):Mpi_fit);
    Dm_tot_fit = (Use_Landau_Gauge?Dm_tot_fit_Land:Dm_tot_fit_Feyn);

    //compute ratios between twisted, untwisted and conserved in Feynman Gauge
    ratio_exch_tw_untw.distr_list.push_back( (Dm_exch_fit_Feyn)/Dm_exch_untwist_fit_Feyn);
    ratio_exch_loc_cons.distr_list.push_back( Dm_exch_untwist_fit_Feyn/Dm_cons_fit);
    ratio_exch_tw_cons.distr_list.push_back( (Dm_exch_fit_Feyn)/Dm_cons_fit);

  

    //if A40.48 load disc on galileo
    if(m_data.Tag[i] == "A40.48_96") {
      distr_t_list Hand_Land_gal= Corr.corr_t(hand_Land_Gal.col(0)[0], "");
      //distr_t_list Hand_Land_gal_new_norm= Corr.corr_t(hand_Land_Gal_new_norm.col(0)[0],"");
      distr_t_list Dm_Land_gal_eff_slope= Corr.effective_slope_t(Hand_Land_gal, Pi_iso_correlated_distr , "../data/Mpi_twisted_ov_X/"+p+"/Landau_gauge/dm_gal_hand."+m_data.Tag[i]);
      //distr_t_list Dm_Land_gal_eff_slope_new_norm= Corr.effective_slope_t(Hand_Land_gal_new_norm, Pi_iso_correlated_distr , "../data/Mpi_twisted_ov_X/"+p+"/Landau_gauge/dm_gal_hand_new_norm."+m_data.Tag[i]);
    }

    
    //choose the RC to be used
    distr_t RC = Use_M3?Za_M4*Za_M4:(Za_M4/Za_M4);
     

    //push_back the results for Dm2/X  and Dm2 depending on whether Landau or Feynman gauge is used!
    Dm_ov_X_distr_list.distr_list.push_back( (2.0*RC*Dm_tot_fit*Mpi_fit  + (Use_M3?1.0:0.0)*((kappa/(2.0*M_PI))*(1.0/ll)*( Mpi_fit + 2.0/ll)))/(Use_fp?(fp_fit*fp_fit):(Mpi_fit*Mpi_fit))) ;
    Dm_ov_X_ren_univ_sub_distr_list.distr_list.push_back( (2.0*Mpi_fit*Za_distr*Za_distr*Dm_tot_fit + (kappa/(2.0*M_PI))*(1.0/ll)*( Mpi_fit + 2.0/ll))/(Use_fp?(fp_fit*fp_fit):(Mpi_fit*Mpi_fit)));
    Dm2_ren_univ_sub_distr_list.distr_list.push_back( (2.0*Mpi_fit*Za_distr*Za_distr*Dm_tot_fit + (kappa/(2.0*M_PI))*(1.0/ll)*( Mpi_fit + 2.0/ll)));
    Dm2_ren_univ_sub_M1_distr_list.distr_list.push_back( (2.0*Mpi_fit*Za_distr_M1*Za_distr_M1*Dm_tot_fit + (kappa/(2.0*M_PI))*(1.0/ll)*( Mpi_fit + 2.0/ll)));
    diff_Dm2_ren_univ_sub_M1_M2_distr_list.distr_list.push_back( 2.0*Mpi_fit*(Za_distr_M1*Za_distr_M1-Za_distr*Za_distr)*Dm_tot_fit);



    //push_back X_pi
    X_pi_distr_list.distr_list.push_back(X_pi_fit );
   

    
    if(!Use_JB_distribution) { //resample from Gaussian distribution
      Compute_covariance_matrix(UseJack, CovMatrixPion[i],3, Mpi_distr_list.distr_list[i].distr, Dm_ov_X_distr_list.distr_list[i].distr, fp_fit_distr_list.distr_list[i].distr);
      Eigen::MatrixXd Corr_Matrix_Pion;
      Compute_correlation_matrix(UseJack, Corr_Matrix_Pion,3, Mpi_distr_list.distr_list[i].distr, Dm_ov_X_distr_list.distr_list[i].distr, fp_fit_distr_list.distr_list[i].distr);
      cout<<"Printing corr matrix pion for iens: "<<i<<endl;
      cout<<Corr_Matrix_Pion<<endl;
      
    }

    cout<<"################### ANALYZING ENSEMBLE : "<<m_data.Tag[i]<<" ############################ "<<endl;
    cout<<"PRINTING INFO: "<<endl;
    cout<<"Mpi: "<<Mpi_fit_naive.ave()<<" +- "<<Mpi_fit_naive.err()<<endl;
    cout<<"Mpi (dim) : "<<(lat_spacing_inv*Mpi_fit_naive).ave()<<" +- "<<(lat_spacing_inv*Mpi_fit_naive).err()<<endl;
    cout<<"Mpi*L :"<<Mpi_fit_naive.ave()*ll<<endl;
    cout<<"L/a: "<<ll<<endl;
    cout<<"fp: "<<fp_fit_naive.ave()<<" +- "<<fp_fit_naive.err()<<endl;
    cout<<"X_pi: "<<X_pi_fit.ave()<<" +- "<<X_pi_fit.err()<<endl;
    cout<<"ainv (GeV-1): "<<lat_spacing_inv<<endl;
    cout<<"%%%%%%%%%  RCs info  %%%%%%%%%"<<endl;
    cout<<"------------    ZV  ---------------"<<endl;
    cout<<"Zv (M1): "<<Zv_M1<<" +- "<<Zv_M1_e<<endl;
    cout<<"Zv (M2): "<<Zv<<" +- "<<Zv_e<<endl;
    cout<<"Zv (M3): "<<Zv_M3.ave()<<" +- "<<Zv_M3.err()<<endl;
    cout<<"------------    ZA  ---------------"<<endl;
    cout<<"Za (M1): "<<Za_M1<<" +- "<<Za_M1_e<<endl;
    cout<<"Za (M2): "<<Za<<" +- "<<Za_e<<endl;
    cout<<"Za (M3): "<<Za_M3.ave()<<" +- "<<Za_M3.err()<<endl;
    cout<<"Za (M3-bis): "<<Za_M4.ave()<<" +- "<<Za_M4.err()<<endl;
    cout<<"-----------    Zp/Zs  -------------"<<endl;
    cout<<"Zp/Zs (M3): "<<Zp_ov_Zs.ave()<<" +- "<<Zp_ov_Zs.err()<<endl;
    cout<<"%%%%%%%%% END Rcs info %%%%%%%%"<<endl;
    cout<<"%%%%%%%%% Feynman Gauge%%%%%%%%%"<<endl;
    cout<<"Dm(exchange, bare): "<<Dm_tot_fit_exch_Feyn.ave()<<" +- "<<Dm_tot_fit_exch_Feyn.err()<<endl;
    cout<<"Dm(hand, bare): "<<Dm_tot_fit_hand_Feyn.ave()<<" +- "<<Dm_tot_fit_hand_Feyn.err()<<endl;
    cout<<"Dm(tot, bare): "<<Dm_tot_fit_Feyn.ave()<<" +- "<<Dm_tot_fit_Feyn.err()<<endl;
    cout<<"Dm2/X: "<<Dm2_ov_X_Feyn.ave()<<" +- "<<Dm2_ov_X_Feyn.err()<<endl;
    cout<<"%%%%%%%%% Landau Gauge%%%%%%%%%"<<endl;
    cout<<"Dm(exchange, bare): "<<Dm_tot_fit_exch_Land.ave()<<" +- "<<Dm_tot_fit_exch_Land.err()<<endl;
    cout<<"Dm(hand, bare): "<<Dm_tot_fit_hand_Land.ave()<<" +- "<<Dm_tot_fit_hand_Land.err()<<endl;
    cout<<"Dm(tot, bare): "<<Dm_tot_fit_Land.ave()<<" +- "<<Dm_tot_fit_Land.err()<<endl;
    cout<<"Dm2/X: "<<Dm2_ov_X_Land.ave()<<" +- "<<Dm2_ov_X_Land.err()<<endl;
    cout<<"#########################################################################################################"<<endl;
  }
 
  //everything is set up. Start fitting!
  LatticeInfo L_info(CURRENT_TYPE);
  
  bootstrap_fit<X_t_M,Y_t_M> bf(nboots);
  bf.set_warmup_lev(1);
  bf.Set_number_of_measurements(Nens);
  bf.Set_verbosity(verbose);
  
  //Add parameters
  bf.Add_par("chir", (0.5e-4)/pow(f0,2), (0.3e-6)/pow(f0,2));
  bf.Add_par("log", 3.0 , 0.2);
  bf.Add_par("A_1", -1.0, 0.02);
  bf.Add_par("D", 2e-3, 1e-5);
  bf.Add_par("lg_a", 0.1, 1e-3);
  bf.Add_par("F_m", 3.0, 0.1);
  bf.Add_par("F_4", 1.0, 0.1);
  if(!Use_Fa_prior) bf.Add_par("F_a", 1, 0.1);
  else bf.Add_prior_pars({"F_a"});
  if(!Use_A2_Dm_prior)  {bf.Add_par("Dm", 0.5, 1e-3); bf.Add_par("A_2", 1.0,0.01); }
  else bf.Add_prior_pars({"Dm", "A_2"});

  bf.Add_prior_pars({ "ainv0", "ainv1", "ainv2","Za0", "Za1", "Za2"    });

 
  //Fix some parameters to make test
  bf.Fix_par("F_m",1.0);
  bf.Fix_par("F_a",0.0);
  //bf.Fix_par("F_4",0.0);
  bf.Fix_par("Dm",0.0);
  bf.Fix_par("D",0.0);
  bf.Fix_par("lg_a", 0.0);
  bf.Fix_par("A_2", 0.0);
  //bf.Fix_par("A_1",0.0);
  if((Use_fp || !Enable_f0)) bf.Fix_par("log", 3.0);



  //Add List of parameters to be released after first minimization
  if(!Use_M3)  bf.Fix_n_release({"ainv0", "ainv1", "ainv2", "Za0", "Za1", "Za2"});

  else  {
    bf.Fix_n_release({"ainv0", "ainv1", "ainv2"});
    bf.Fix_par("Za0", 1.0);
    bf.Fix_par("Za1", 1.0);
    bf.Fix_par("Za2", 1.0);
  }

  double Npars= bf.Get_number_of_fit_pars();


   
  Vfloat L, m_l, T;
  Read_pars_from_ensemble_tag(m_data.Tag, m_l, L, T);

  


  //define fitting function

  //############################################################################################################

  bf.ansatz =  [=](const X_t_M &p, const Y_t_M &ip) -> double {
		 
		 
    double ainv = p.ainv[ip.ibeta];

    double csi_val = ip.csi*(Enable_f0?pow(ainv,2):1.0);

    double csi_resc= pow(4.0*M_PI,2)*csi_val;

    double scale = Use_fp?csi_resc:1.0;

    double scale_fve = Use_fp?pow(ip.Mpi/ip.fp,2):1.0;

    

    double a=1.0/ainv;

    double L = ip.L/ainv;

    double Mp= ip.Mpi*ainv;

    double Mp2 = Mp*Mp;

    double Mp3= Mp*Mp*Mp;

    double ML = ip.Mpi*ip.L;

         
    double SD_FVE = (L>0.0)?(e2/3.0)*scale_fve*r0*(Mp2/pow(ML,3))*(p.F_m+ p.F_4/(ML)+ (3.0/r0)*p.F_a*(pow(a,2))):0.0;  //in p.F4 ML -> Mp2

    double FVE_universal =  (Use_M3?0.0:((L>0.0)?(scale_fve*(kappa*alpha/ML)*( 1.0 + 2.0/ML)):0.0));


    double fitted_value = scale*(e2/pow(4.0*M_PI,2))*( p.chir/csi_val -p.log*log(csi_val) + p.A_1 + p.A_2*csi_val + (p.D/(csi_val))*a*a + p.Dm*a*a);


    fitted_value += -FVE_universal + SD_FVE;
     

    
    

    return X_phys_val*fitted_value;
  };
  

  bf.measurement = [=](const X_t_M& p,const Y_t_M& ip) -> double {

    double Za = p.Za[ip.ibeta];
            
    double m1 = (e2/2)*pow(Za,2)*ip.Dm_ov_X_val;
      
    return X_phys_val*m1;
  };
  

  bf.error =  [=](const X_t_M& p,const  Y_t_M &ip) -> double {

    double Za = p.Za[ip.ibeta];

    //Za=Zv;

       
    double m1 = X_phys_val*(e2/2)*pow(Za,2)*ip.Dm_ov_X_err;

    return m1;
  };

 

  
  //############################################################################################################

  //init random Number Generators

  GaussianMersenne GM(54353);
  RandomMersenne RM(23423, UseJack?(Njacks-1):(Nboots-1));

  Eigen::MatrixXd CovMatrixInput(9,9); //covariance matrix of input parameters
  Eigen::VectorXd Ave_input_parameters(9);

  vector<boot_fit_data<X_t_M>> Bt_fit;

  //we store here all resampled data. 
  vector<vector<vector<Y_t_M>>> all_ens(Nbranches);
  for(auto &all_ens_br : all_ens) {
    all_ens_br.resize(nboots);
    for(auto & all_ens_ib: all_ens_br) {
      all_ens_ib.resize(Nens);
      for(int iens=0; iens<Nens;iens++) {

	all_ens_ib[iens].L=L[iens];
	if(m_data.Tag[iens].substr(0,1) == "A") all_ens_ib[iens].ibeta=0;
	else if(m_data.Tag[iens].substr(0,1) == "B") all_ens_ib[iens].ibeta=1;
	else all_ens_ib[iens].ibeta=2;
	all_ens_ib[iens].Mpi_err = Mpi_distr_list.err(iens);
	all_ens_ib[iens].Dm_ov_X_err = Dm_ov_X_distr_list.err(iens);
	

      }
    }
  }

  //store the chi2 evaluated over mean values
  Vfloat chi2_MV;
 
  for(int ibranch=0; ibranch < Nbranches; ibranch++) {


    
    ReadBranch(ibranch, CovMatrixInput, Ave_input_parameters );


    //clear all measurements
    bf.Clear_priors();
    bf.Clear_input_pars();
    
   
    //generate bootstrap data
    for(int iboot=0; iboot<nboots; iboot++) {
      bf.ib= &iboot;
      Vfloat lat_input = Covariate(CovMatrixInput, Ave_input_parameters, GM);
      if(iboot==0) for(unsigned int i=0; i<lat_input.size();i++)  lat_input[i] = Ave_input_parameters[i];
      if(Use_JB_distribution && iboot != 0) {
	for(int icov=0;icov<(signed)Ave_input_parameters.size();icov++) {
	  lat_input[icov] = Ave_input_parameters[icov] + (lat_input[icov]-Ave_input_parameters[icov])/sqrt(Njacks-1.0);
	}
      }
    
     
      double f0_resampled= lat_input[2];
      double Za0 = L_info.Retrieve_Za("A",ibranch).first;
      double Za1 = L_info.Retrieve_Za("B",ibranch).first;
      double Za2 = L_info.Retrieve_Za("D", ibranch).first;

      if(iboot != 0) {
	if(!Use_JB_distribution) {
	Za0= gauss(L_info.Retrieve_Za("A", ibranch),GM);
	Za1= gauss(L_info.Retrieve_Za("B", ibranch),GM);
	Za2= gauss(L_info.Retrieve_Za("D", ibranch) ,GM);
	}
	else {
	  Za0= gauss(L_info.Retrieve_Za("A", ibranch).first, (1.0/sqrt(Njacks-1))*L_info.Retrieve_Za("A",ibranch).second,GM);
	  Za1= gauss(L_info.Retrieve_Za("B", ibranch).first, (1.0/sqrt(Njacks-1))*L_info.Retrieve_Za("B",ibranch).second,GM);
	  Za2= gauss(L_info.Retrieve_Za("D", ibranch).first, (1.0/sqrt(Njacks-1))*L_info.Retrieve_Za("D",ibranch).second,GM);
	}
      }
      
      if(Use_M3) { Za0 =1.0; Za1= 1.0; Za2= 1.0;}
      if(Use_Fa_prior) bf.Append_to_prior("F_a", 0.0,10);
      if(Use_A2_Dm_prior){ bf.Append_to_prior("Dm", 0.0, 10);  bf.Append_to_prior("A_2", 0.0, 10);} //prior on A_2 and Dm
      bf.Append_to_prior("ainv0", lat_input[3], sqrt(CovMatrixInput(3,3)));
      bf.Append_to_prior("ainv1", lat_input[4], sqrt(CovMatrixInput(4,4)));
      bf.Append_to_prior("ainv2", lat_input[5], sqrt(CovMatrixInput(5,5)));
      bf.Append_to_prior("Za0", Za0, L_info.Retrieve_Za("A", ibranch).second);
      bf.Append_to_prior("Za1", Za1, L_info.Retrieve_Za("B", ibranch).second);
      bf.Append_to_prior("Za2", Za2, L_info.Retrieve_Za("D", ibranch).second);


      int k= RM();
      for(int imeas=0; imeas <Nens; imeas++) {
	if(iboot==0) {
	    all_ens[ibranch][iboot][imeas].Mpi= Mpi_distr_list.ave(imeas);
	    all_ens[ibranch][iboot][imeas].Dm_ov_X_val = Dm_ov_X_distr_list.ave(imeas);
	    all_ens[ibranch][iboot][imeas].fp = fp_fit_distr_list.ave(imeas);
	    all_ens[ibranch][iboot][imeas].f0 = Ave_input_parameters[2];
	    all_ens[ibranch][iboot][imeas].csi = pow(Mpi_distr_list.ave(imeas),2)/pow(4.0*M_PI*(Enable_f0?Ave_input_parameters[2]:fp_fit_distr_list.ave(imeas)),2);

	}
	else if(!Use_JB_distribution) {
	  Eigen::VectorXd vec(CovMatrixPion[imeas].rows());
	 
	  vec<<Mpi_distr_list.ave(imeas),Dm_ov_X_distr_list.ave(imeas),fp_fit_distr_list.ave(imeas);
	  Vfloat res_meas(3,0.0);
	 
	  res_meas= Covariate(CovMatrixPion[imeas],vec, GM);
	  all_ens[ibranch][iboot][imeas].Mpi= res_meas[0];
	  all_ens[ibranch][iboot][imeas].Dm_ov_X_val = res_meas[1];
	  all_ens[ibranch][iboot][imeas].fp = res_meas[2];
	  all_ens[ibranch][iboot][imeas].f0 = f0_resampled;
	  all_ens[ibranch][iboot][imeas].csi = pow(res_meas[0],2)/pow(4.0*M_PI*(Enable_f0?f0_resampled:res_meas[2]),2);
	 	 
	}
      
	else {
	  double M,FP;
	  M= Mpi_distr_list.distr_list[imeas].distr[k];
	  FP=fp_fit_distr_list.distr_list[imeas].distr[k];
	  all_ens[ibranch][iboot][imeas].Mpi = M;
	  all_ens[ibranch][iboot][imeas].Dm_ov_X_val = Dm_ov_X_distr_list.distr_list[imeas].distr[k];
	  all_ens[ibranch][iboot][imeas].fp = FP;
	  all_ens[ibranch][iboot][imeas].f0 = f0_resampled;
	  all_ens[ibranch][iboot][imeas].csi = pow(M,2)/pow(4.0*M_PI*(Enable_f0?f0_resampled:FP),2);
	  
	}
      }
      
    }
    
    bf.Append_to_input_par(all_ens[ibranch]);
    Bt_fit.push_back(bf.Perform_bootstrap_fit());

    chi2_MV.push_back( Bt_fit[ibranch].chi2[0]);
    
  }
  
 
  //print the data
  
  //define lambda functions

  auto FVE = [=](double ML, double fp, double mp) {

	       if(Use_M3) return 0.0;
	       double resc= Use_fp?pow(mp/fp,2):1.0;
	       return resc*X_phys_val*(kappa*alpha/ML)*( 1.0 + 2.0/ML);
  };

  auto SD_FVE = [=](double ML, double M, double FP, double ainv,  double F_a, double F_m, double F_4) {

		  double resc= Use_fp?pow(M/FP,2):1.0;
		  double Mp= M*ainv;
		  double Mp2 = Mp*Mp;
		  double Mp3 = Mp*Mp*Mp;
		  double SDE =(e2/3.0)*(Mp2/pow(ML,3))*r0*(F_m + F_4/(ML) +  (3.0/r0)*F_a*(1.0/pow(ainv,2))); //in F_4 ML -> Mp2
		  return X_phys_val*SDE*resc;
  };


  auto cont_ansatz = [=](X_t_M p, double csi) -> double {


		       double scale = Use_fp?pow(4.0*M_PI,2)*csi:1.0;
		       double val = scale*X_phys_val*(e2/pow(4.0*M_PI,2))*( p.chir/csi - p.log*log(csi) + p.A_1 + p.A_2*csi);

		    
		       return val;
  };

  auto ANSATZ = [=](const X_t_M &p, const Y_t_M &ip) -> double {
		 
		 
    double ainv = p.ainv[ip.ibeta];

    double csi_val = ip.csi;

    double csi_resc= pow(4.0*M_PI,2)*csi_val;

    double scale = Use_fp?csi_resc:1.0;

    double scale_fve = Use_fp?pow(4.0*M_PI,2)*ip.csi:1.0;

    

    double a=1.0/ainv;

    double L = ip.L/ainv;

    double Mp= ip.Mpi*ainv;

    double Mp2 = Mp*Mp;

    double Mp3= Mp*Mp*Mp;

    double ML = ip.Mpi*ip.L;

         
    double SD_FVE= (L>0.0)?(e2/3.0)*scale_fve*r0*(Mp2/pow(ML,3))*(p.F_m+ p.F_4/(ML) + (3.0/r0)*p.F_a*(pow(a,2))):0.0;   //in p.F_4 ML -> Mp2

    double FVE_universal =  (L>0.0)?(scale_fve*(kappa*alpha/ML)*( 1.0 + 2.0/ML)):0.0;


    double fitted_value = scale*(e2/pow(4.0*M_PI,2))*( p.chir/csi_val -p.log*log(csi_val) + p.A_1 + p.A_2*csi_val + (p.D/(csi_val))*a*a + p.Dm*a*a);


    fitted_value += -FVE_universal + SD_FVE;
     

    
    

    return X_phys_val*fitted_value;
  };

  auto Artifacts = [=](const X_t_M &p, double csi, int ibeta) {

		    double a = 1.0/p.ainv[ibeta];

		    double csi_resc = pow(4.0*M_PI,2)*csi;
		    double scale= Use_fp?csi_resc:1.0;
		    return X_phys_val*scale*(e2/pow(4.0*M_PI,2))*( (p.D/csi)*a*a + p.Dm*a*a);

		  };



 
  
  
  Vfloat Csi(1000);
  Vfloat vols(1000);
  for(unsigned int m=0; m<Csi.size();m++) Csi[m] = 0.0001*m+ 0.0002;
  for(unsigned int vol=0;vol<vols.size();vol++) vols[vol] = 10 +vol*0.5;
  VVVfloat SD_subtracted_data, Csi_measured, raw_data, Dm2_raw_data;
  VVVfloat Univ_subtracted_data, Dm2_Univ_subtracted_data;
  VVVfloat Dm2_A40_slice;
  VVVVfloat Fitted_func;
  VVfloat Physical_point;
  VVVfloat A40_slice;
  VVfloat C,LOG,A1,A2,K,F_4, Fa,D, Dm, Ainv0, Ainv1, Ainv2,  Za0, Za1, Za2;
  cascade_resize(SD_subtracted_data, Vint{bf.Get_number_of_measurements(), Nbranches, nboots});
  cascade_resize(Univ_subtracted_data, Vint{bf.Get_number_of_measurements(), Nbranches, nboots});
  cascade_resize(Dm2_Univ_subtracted_data, Vint{bf.Get_number_of_measurements(), Nbranches, nboots});
  cascade_resize(Csi_measured, Vint{bf.Get_number_of_measurements(), Nbranches, nboots});
  cascade_resize(Fitted_func, Vint{4, (int)Csi.size(), Nbranches, nboots});
  cascade_resize(Physical_point, Vint{Nbranches, nboots});
  cascade_resize(A40_slice, Vint{ (int)vols.size(), Nbranches,nboots});
  cascade_resize(Dm2_A40_slice, Vint{ (int)vols.size(), Nbranches,nboots});
  cascade_resize(raw_data, Vint{bf.Get_number_of_measurements(), Nbranches, nboots});
  cascade_resize(Dm2_raw_data, Vint{bf.Get_number_of_measurements(), Nbranches, nboots});
  cascade_resize(C, Vint{Nbranches, nboots});
  cascade_resize(LOG, Vint{Nbranches, nboots});
  cascade_resize(A1, Vint{Nbranches, nboots});
  cascade_resize(A2, Vint{Nbranches, nboots});
  cascade_resize(K, Vint{Nbranches, nboots});
  cascade_resize(F_4, Vint{Nbranches, nboots});
  cascade_resize(Fa, Vint{Nbranches, nboots});
  cascade_resize(D, Vint{Nbranches, nboots});
  cascade_resize(Dm, Vint{Nbranches, nboots});
  cascade_resize(Ainv0, Vint{Nbranches, nboots});
  cascade_resize(Ainv1, Vint{Nbranches, nboots});
  cascade_resize(Ainv2, Vint{Nbranches, nboots});
  cascade_resize(Za0, Vint{Nbranches, nboots});
  cascade_resize(Za1, Vint{Nbranches, nboots});
  cascade_resize(Za2, Vint{Nbranches, nboots});
 
 
  
  for(int ibranch=0; ibranch<Nbranches;ibranch++) {
    for(int iboot=0; iboot<nboots;iboot++) {
      X_t_M P = Bt_fit[ibranch].par[iboot];
         //store fitted params
      C[ibranch][iboot]= P.chir;
      LOG[ibranch][iboot] =P.log;
      A1[ibranch][iboot] =P.A_1;
      A2[ibranch][iboot]= P.A_2;
      K[ibranch][iboot]= P.F_m;
      F_4[ibranch][iboot] = P.F_4;
      Fa[ibranch][iboot] =P.F_a;
      D[ibranch][iboot]= P.D;
      Dm[ibranch][iboot] =P.Dm;
      Ainv0[ibranch][iboot] = P.ainv[0];
      Ainv1[ibranch][iboot] = P.ainv[1];
      Ainv2[ibranch][iboot] = P.ainv[2];
      Za0[ibranch][iboot] = P.Za[0];
      Za1[ibranch][iboot] = P.Za[1];
      Za2[ibranch][iboot] = P.Za[2];
      for(unsigned int ivol=0; ivol<vols.size();ivol++) {
	Y_t_M new_p;
	//determine A40.48 ensemble
	int Tag_A40_48;
	bool Tag_found=false;
	for(int im=0;im<Nens;im++) {
	    if(m_data.Tag[im]=="A40.48_96")
	      {
		Tag_A40_48 = im;
		Tag_found=true;
		break;
	      }
	}
	if(!Tag_found) {
	  for(int im=0;im<Nens;im++) {
	    if(m_data.Tag[im]=="A40.40_80")
	      {
		Tag_A40_48 = im;
		Tag_found=true;
		break;
	      }
	  }
	}
	if(!Tag_found)  crash("Cannot find ensemble A40.48_96 in the ensemble list");
	double f0_boot = all_ens[ibranch][iboot][0].f0/P.ainv[0];
	new_p.ibeta=0;
	new_p.Mpi= all_ens[ibranch][iboot][Tag_A40_48].Mpi;  //pion mass of A40.XX in lattice units
	new_p.L = vols[ivol];
	new_p.fp = all_ens[ibranch][iboot][Tag_A40_48].fp;
	double csi_A40 = pow(new_p.Mpi/(4.0*M_PI*new_p.fp),2);
	new_p.csi = Enable_f0?pow(new_p.Mpi,2)/pow(4.0*M_PI*f0_boot,2):csi_A40;
	//new_p.fp = new_p.Mpi/(4.0*M_PI*sqrt(new_p.csi));
	A40_slice[ivol][ibranch][iboot] = cont_ansatz(P,csi_A40) + Artifacts( P, csi_A40, 0)  + SD_FVE(vols[ivol]*new_p.Mpi, new_p.Mpi, new_p.fp, P.ainv[0], P.F_a, P.F_m, P.F_4);
	Dm2_A40_slice[ivol][ibranch][iboot] = pow(new_p.fp/fpi_phys,2)*A40_slice[ivol][ibranch][iboot];
      }
      
      for(int imeas=0; imeas< bf.Get_number_of_measurements(); imeas ++) {
	Y_t_M E = all_ens[ibranch][iboot][imeas];
	double fval = Enable_f0?(E.f0/P.ainv[E.ibeta]):E.fp;
	SD_subtracted_data[imeas][ibranch][iboot] = bf.measurement(P, E) + FVE(E.L*E.Mpi, E.fp, E.Mpi) - SD_FVE(E.L*E.Mpi, E.Mpi, E.fp, P.ainv[E.ibeta], P.F_a, P.F_m, P.F_4);
	Univ_subtracted_data[imeas][ibranch][iboot] = bf.measurement(P,E)+ FVE(E.L*E.Mpi, E.fp, E.Mpi);
	Dm2_Univ_subtracted_data[imeas][ibranch][iboot] = pow(E.fp/fpi_phys,2)*bf.measurement(P,E)+ pow(E.fp/fpi_phys,2)*FVE(E.L*E.Mpi, E.fp, E.Mpi);
	Dm2_raw_data[imeas][ibranch][iboot] = pow(E.fp/fpi_phys,2)*bf.measurement(P,E);
	raw_data[imeas][ibranch][iboot] = bf.measurement(P,E);
	Csi_measured[imeas][ibranch][iboot] = pow(E.Mpi,2)/pow(4.0*M_PI*fval,2);
      }
     
      for(int ibeta=0; ibeta < 4;ibeta++) {
	for(unsigned int m=0; m<Csi.size(); m++){
	  if(ibeta<3) {
	  
	    Y_t_M point;
	    point.L = -1.0;
	    point.ibeta=ibeta;
	    point.Mpi = 1.0;
	    point.csi= Csi[m];
	    Fitted_func[ibeta][m][ibranch][iboot] = ANSATZ(P,point);
	    
	    
	    
	  }
	  else Fitted_func[ibeta][m][ibranch][iboot] = cont_ansatz(P, Csi[m]) ;
	}
      }
  
      //Add physical point
      Physical_point[ibranch][iboot] = cont_ansatz(P, csi_phys );
   
    }
  }




 

  Vfloat SD_subtracted_val, SD_subtracted_err, CSI, raw_data_val, raw_data_err;
  Vfloat Univ_subtracted_val, Univ_subtracted_err;
  Vfloat Dm2_Univ_subtracted_val, Dm2_Univ_subtracted_err, Dm2_raw_data_val, Dm2_raw_data_err;
  for(int imeas=0;imeas<bf.Get_number_of_measurements();imeas++) {
    SD_subtracted_val.push_back( Boot_ave(SD_subtracted_data[imeas]));
    SD_subtracted_err.push_back( Boot_err(SD_subtracted_data[imeas], Rerr,1));
    Univ_subtracted_val.push_back(Boot_ave(Univ_subtracted_data[imeas]));
    Univ_subtracted_err.push_back(Boot_err(Univ_subtracted_data[imeas], Rerr,1));
    
    Dm2_Univ_subtracted_val.push_back(Boot_ave(Dm2_Univ_subtracted_data[imeas]));
    Dm2_Univ_subtracted_err.push_back(Boot_err(Dm2_Univ_subtracted_data[imeas],Rerr,1));
    Dm2_raw_data_val.push_back(Boot_ave(Dm2_raw_data[imeas]));
    Dm2_raw_data_err.push_back(Boot_err(Dm2_raw_data[imeas], Rerr,1));
    raw_data_val.push_back(Boot_ave(raw_data[imeas]));
    raw_data_err.push_back(Boot_err(raw_data[imeas], Rerr,1));
    CSI.push_back( Boot_ave(Csi_measured[imeas]));
  }
  
  VVfloat Fitted_func_val, Fitted_func_err;
  Vfloat A40_slice_val, A40_slice_err, A40_raw_val, A40_raw_err;
  Vfloat Dm2_A40_slice_val, Dm2_A40_slice_err;
 
  cascade_resize(Fitted_func_val, {4,(int)Csi.size()});
  cascade_resize(Fitted_func_err, {4,(int)Csi.size()});


  for(int ibeta=0;ibeta<4;ibeta++) {
    for(unsigned int m=0; m<Csi.size();m++) {
	
	Fitted_func_val[ibeta][m] = Boot_ave(Fitted_func[ibeta][m]);
	Fitted_func_err[ibeta][m] = Boot_err(Fitted_func[ibeta][m], Rerr, 1);

    }
  }

  for(auto &A40_boot : A40_slice) {
    A40_slice_val.push_back( Boot_ave(A40_boot));
    A40_slice_err.push_back(Boot_err(A40_boot, Rerr,1));
  }

  for(auto &Dm2_A40_boot : Dm2_A40_slice) {
    Dm2_A40_slice_val.push_back(Boot_ave(Dm2_A40_boot));
    Dm2_A40_slice_err.push_back(Boot_err(Dm2_A40_boot, Rerr,1));
  }
  //plot the result

  Vfloat Mpi_L_A40, meas_A40, err_meas_A40, vol_A40, Dm2_A40_raw_val, Dm2_A40_raw_err, meas_Dm2_A40, err_meas_Dm2_A40;
  for(int iens=0; iens<Nens;iens++)
    if(m_data.Tag[iens].substr(0,3)=="A40") {
      vol_A40.push_back( L[iens]);
      Mpi_L_A40.push_back( L[iens]*Mpi_distr_list.ave(iens));
      meas_A40.push_back(Univ_subtracted_val[iens]);
      err_meas_A40.push_back(Univ_subtracted_err[iens]);
      A40_raw_val.push_back(raw_data_val[iens]);
      A40_raw_err.push_back(raw_data_err[iens]);
      meas_Dm2_A40.push_back(Dm2_Univ_subtracted_val[iens]);
      err_meas_Dm2_A40.push_back(Dm2_Univ_subtracted_err[iens]);
      Dm2_A40_raw_val.push_back(Dm2_raw_data_val[iens]);
      Dm2_A40_raw_err.push_back(Dm2_raw_data_err[iens]);
    }



  string print_path = (CURRENT_TYPE=="CONSERVED")?"conserved":"local";
  string exch_or_tot = IncludeDisconnected?"tot":"exch";

    
  


   
  //save data in files
  boost::filesystem::create_directory("../data");
  boost::filesystem::create_directory("../data/Mpi_twisted_ov_X/"+print_path);
  Print_To_File({}, {Csi, Fitted_func_val[0], Fitted_func_err[0],  Fitted_func_val[1], Fitted_func_err[1], Fitted_func_val[2], Fitted_func_err[2], Fitted_func_val[3], Fitted_func_err[3]} , "../data/Mpi_twisted_ov_X/"+print_path+LorF+"/Fitted_func_"+exch_or_tot+".dat", "OUT", "#Csi     beta=1.90      beta=1.95      beta=2.10           cont");
 
  
  Print_To_File(m_data.Tag, { CSI, Mpi_distr_dim_list.ave(), L, (Mpi_distr_list*L).ave(),  SD_subtracted_val, SD_subtracted_err}, "../data/Mpi_twisted_ov_X/"+print_path+LorF+"/sd_subtracted_data_"+exch_or_tot+".dat", "", "#Ens   #Csi     #Mpi #L/a  #Mpi L   #Mpi^2_+ - Mpi^2_0      #err");

  Print_To_File(m_data.Tag, { CSI, Mpi_distr_dim_list.ave(), L, (Mpi_distr_list*L).ave(), X_pi_distr_list.ave(), X_pi_distr_list.err(), Mpi_distr_list.ave(), Mpi_distr_list.err(), fp_fit_distr_list.ave(), fp_fit_distr_list.err()}, "../data/Mpi_twisted_ov_X/"+print_path+"/Xi_pi_"+exch_or_tot+".dat", "", "#Ens   #Csi     #Mpi #L/a     #Mpi*L       #Xi     #Mpi    #fp");
 
  Print_To_File(m_data.Tag, { CSI, Mpi_distr_dim_list.ave(), Mpi_distr_dim_list.err(),  L, (Mpi_distr_list*L).ave(), Dm2_Univ_subtracted_val, Dm2_Univ_subtracted_err,  Univ_subtracted_val, Univ_subtracted_err, raw_data_val, raw_data_err, ((X_phys_val*e2/2.0)*Dm_ov_X_distr_list).ave(), ( (X_phys_val*e2/2.0)*Dm_ov_X_distr_list).err()}, "../data/Mpi_twisted_ov_X/"+print_path+LorF+"/data_"+exch_or_tot+".dat", "", "#Ens      Csi   #Mpi #Mpi_err   #L/a   #Mpi*L Dm2_univ_sub(resampled)  X_univ_sub(resampled)      X(resampled)     X(bare)");

  Print_To_File(m_data.Tag, { CSI, Mpi_distr_dim_list.ave(), L, (Mpi_distr_list*L).ave(), ((e2/2.0)*(1.0/(1e-3*2.0*MPiPhys))*Dm2_ren_univ_sub_distr_list).ave(),  ((e2/2.0)*(1.0/(1e-3*2.0*MPiPhys))*Dm2_ren_univ_sub_distr_list).err(),  ((e2/2.0)*(1.0/(1e-3*2.0*MPiPhys))*Dm2_ren_univ_sub_M1_distr_list).ave(),  ((e2/2.0)*(1.0/(1e-3*2.0*MPiPhys))*Dm2_ren_univ_sub_M1_distr_list).err(),   ((e2/2.0)*(1.0/(1e-3*2.0*MPiPhys))*diff_Dm2_ren_univ_sub_M1_M2_distr_list).ave(),  ((e2/2.0)*(1.0/(1e-3*2.0*MPiPhys))*diff_Dm2_ren_univ_sub_M1_M2_distr_list).err()  , ((X_phys_val*e2/2.0)*Dm_ov_X_ren_univ_sub_distr_list).ave(), ((X_phys_val*e2/2.0)*Dm_ov_X_ren_univ_sub_distr_list).err()  }, "../data/Mpi_twisted_ov_X/"+print_path+LorF+"/data_wo_fit_"+exch_or_tot+".dat", "", "#Ens      Csi   #Mpi   #L/a   #Mpi*L  Dm2_univ_sub_M2  Dm2_univ_sub_M1 diff_Dm2_univ_sub_M1_M2      X_univ_sub  ");
 
  Print_To_File({}, {vol_A40, Mpi_L_A40,  meas_A40, err_meas_A40, A40_raw_val, A40_raw_err, meas_Dm2_A40, err_meas_Dm2_A40, Dm2_A40_raw_val, Dm2_A40_raw_err}, "../data/Mpi_twisted_ov_X/"+print_path+LorF+"/A40_"+exch_or_tot+".dat","", " #L    #Mpi*L         univ_meas            univ_err       raw_meas      raw_err         univ meas(Dm2)    univ_err(Dm2)    raw(Dm2)    err_raw(Dm2)");
  //print A40 prediction for
  Print_To_File({},{vols, A40_slice_val, A40_slice_err, Dm2_A40_slice_val, Dm2_A40_slice_err},"../data/Mpi_twisted_ov_X/"+print_path+LorF+"/A40_pred_"+exch_or_tot+".dat","", "#L   Mpi/X    Dm2   ");

  Print_To_File(m_data.Tag, {L, (Mpi_distr_list*L).ave(),  Mpi_fit_naive_list.ave(), Mpi_fit_naive_list.err(), Mpi_fit_GL_list.ave(), Mpi_fit_GL_list.err(), Mpi_fit_CDH_list.ave(), Mpi_fit_CDH_list.err(), Ratio_CDH_naive_Mpi.ave(), Ratio_CDH_naive_Mpi.err()}, "../data/Mpi_twisted_ov_X/"+print_path+"/FVE_Mpi.dat", "", "#Ens #L  #Mpi*L   #Mpi_naive #Mpi_GL #Mpi_CDH    ratio_CDH_naive");

  Print_To_File(m_data.Tag, {L, (Mpi_distr_list*L).ave(),  fp_fit_naive_list.ave(), fp_fit_naive_list.err(), fp_fit_GL_list.ave(), fp_fit_GL_list.err(), fp_fit_CDH_list.ave(), fp_fit_CDH_list.err(), Ratio_CDH_naive_fp.ave(), Ratio_CDH_naive_fp.err()}, "../data/Mpi_twisted_ov_X/"+print_path+"/FVE_fp.dat", "", "#Ens #L  #Mpi*L   #fp_naive #fp_GL #fp_CDH     ratio_CDH_naive");

  Print_To_File(m_data.Tag, { CSI, Mpi_fit_naive_list.ave(), ratio_LF_exch.ave(), ratio_LF_exch.err(), ratio_LF_hand.ave(), ratio_LF_hand.err(), ratio_LF_tot.ave(), ratio_LF_tot.err(), a_list}, "../data/Mpi_twisted_ov_X/"+print_path+"/Landau_gauge/ratio_LF.dat", "", "#Ens #Csi #Mpi  #exch  #disc  #tot #a_list");

  Print_To_File(m_data.Tag, { CSI, Mpi_fit_naive_list.ave(), ratio_LF_exch_Univ_sub.ave(), ratio_LF_exch_Univ_sub.err(), ratio_LF_tot_Univ_sub.ave(), ratio_LF_tot_Univ_sub.err(), a_list}, "../data/Mpi_twisted_ov_X/"+print_path+"/Landau_gauge/ratio_LF_Univ_sub.dat", "", "#Ens #Csi #Mpi  #exch  #tot #a_list");
  
  
  string command = "echo "+to_string_with_precision(csi_phys,8)+"\t\t"+to_string_with_precision(Boot_ave(Physical_point),8)+"\t\t"+to_string_with_precision(Boot_err(Physical_point, Rerr,1),8)+" > ../data/Mpi_twisted_ov_X/"+print_path.c_str()+"/Phys_val_"+exch_or_tot+".dat";
  system(command.c_str());

  if(IncludeDisconnected) {
    Print_To_File(m_data.Tag, {CSI,Mpi_distr_dim_list.ave(), L, (Mpi_distr_list*L).ave(),  Dm2_disc.ave(), Dm2_disc.err()}, "../data/Mpi_twisted_ov_X/"+print_path+"/dm2_disc.data","","#Ens Csi   Mp L   val err    ");
    Print_To_File(m_data.Tag, {CSI,Mpi_distr_dim_list.ave(), L, (Mpi_distr_list*L).ave(),   Dm2_disc_Landau.ave(), Dm2_disc_Landau.err()}, "../data/Mpi_twisted_ov_X/"+print_path+"/Landau_gauge/dm2_disc.data","","#Ens Csi   Mp L Mpi*L   val err    ");
  }


  //print ratio exch twisted and untwisted
  Print_To_File(m_data.Tag,{CSI, Mpi_distr_dim_list.ave(),a_list, L, (Mpi_distr_list*L).ave(), ratio_exch_tw_untw.ave(), ratio_exch_tw_untw.err(), ratio_exch_loc_cons.ave(), ratio_exch_loc_cons.err(), ratio_exch_tw_cons.ave(), ratio_exch_tw_cons.err()}, "../data/Mpi_twisted_ov_X/"+print_path+"/ratio_exch_tw_untw.data"  , "", "#ENS Csi Mpi a L Mpi*L twist/untwist   loc/cons     twist_loc/cons");



  

   

  cout<<"Physical Pion mass difference: "<< Boot_ave(Physical_point)<<"  +-  "<<Boot_err(Physical_point, Rerr,1)<<" [ MeV ]"<<endl;
  cout<<"On single branches: "<<endl;
  for(int ibr=0; ibr<Nbranches;ibr++) cout<<"Branch: "<<ibr<<": "<<Boot_ave(Physical_point[ibr])<<" +-  "<<Rerr*Boot_err(Physical_point[ibr])<<" [ MeV ]"<<endl;
  cout<<"Average chi2:"<<endl;
  double tot_ch2=0.0;
  for(int ibr=0;ibr<Nbranches;ibr++) {
    double ch2_ibr = accumulate(Bt_fit[ibr].chi2.begin(), Bt_fit[ibr].chi2.end(), 0.0)/Bt_fit[ibr].chi2.size();
    tot_ch2 += ch2_ibr;
    cout<<"Branch: "<<ibr<<" chi2 = "<<ch2_ibr<<endl;
  }
  cout<<"Average over branches: "<<tot_ch2/Nbranches<<endl;

  
  cout<<"Chi2 fitting over mean values: "<<accumulate(chi2_MV.begin(), chi2_MV.end(),0.0)/Nbranches<<endl;
  
  

  //print averaged parameters
  cout<<"Printing averaged params: "<<endl;
  cout<<"C: "<<Boot_ave(C)<<" +- "<<Boot_err(C, Rerr,1)<<endl;
  cout<<"log: "<<Boot_ave(LOG)<<" +- "<<Boot_err(LOG, Rerr,1)<<endl;
  cout<<"A1: "<<Boot_ave(A1)<<" +- "<<Boot_err(A1, Rerr,1)<<endl;
  cout<<"A2: "<<Boot_ave(A2)<<" +- "<<Boot_err(A2, Rerr,1)<<endl;
  cout<<"K: "<<Boot_ave(K)<<" +- "<<Boot_err(K, Rerr,1)<<endl;
  cout<<"F_4: "<<Boot_ave(F_4)<<" +- "<<Boot_err(F_4, Rerr,1)<<endl;
  cout<<"F_a: "<<Boot_ave(Fa)<<" +- "<<Boot_err(Fa, Rerr,1)<<endl;
  cout<<"D: "<<Boot_ave(D)<<" +- "<<Boot_err(D, Rerr,1)<<endl;
  cout<<"Dm: "<<Boot_ave(Dm)<<" +- "<<Boot_err(Dm, Rerr,1)<<endl;
 

  cout<<"############# PRINTING LATTICE SPACING INFOS #################"<<endl;
  
  Eigen::MatrixXd Cov(9,9); //covariance matrix of input parameters
  Eigen::VectorXd Ave(9);
  ReadBranch(0, Cov, Ave );
  cout<<"ainv0 branch 0: "<<Boot_ave(Ainv0[0])<<" +- "<<Rerr*Boot_err(Ainv0[0])<<" Original from branch nr. 0: "<<Ave[3] <<" +- "<<sqrt(Cov(3,3))<<endl;
  cout<<"ainv1 branch 0: "<<Boot_ave(Ainv1[0])<<" +- "<<Rerr*Boot_err(Ainv1[0])<<" Original from branch nr. 0: "<<Ave[4] <<" +- "<<sqrt(Cov(4,4))<<endl;
  cout<<"ainv2 branch 0: "<<Boot_ave(Ainv2[0])<<" +- "<<Rerr*Boot_err(Ainv2[0])<<" Original from branch nr. 0: "<<Ave[5] <<" +- "<<sqrt(Cov(5,5))<<endl;

  cout<<"##############ZA INFO############"<<endl;
  cout<<"beta = 1.90"<<endl;
  cout<<"Za0 branch 0: "<<Boot_ave(Za0[0])<<" +- "<<Rerr*Boot_err(Za0[0])<<" Original from branch nr. 0: "<<L_info.Retrieve_Za("A", 0).first<<" +- "<<L_info.Retrieve_Za("A",0).second<<endl;
  cout<<"Za0 branch 1: "<<Boot_ave(Za0[4])<<" +- "<<Rerr*Boot_err(Za0[4])<<" Original from branch nr. 0: "<<L_info.Retrieve_Za("A", 4).first<<" +- "<<L_info.Retrieve_Za("A",4).second<<endl;
  cout<<"beta = 1.95"<<endl;
  cout<<"Za0 branch 0: "<<Boot_ave(Za1[0])<<" +- "<<Rerr*Boot_err(Za1[0])<<" Original from branch nr. 0: "<<L_info.Retrieve_Za("B", 0).first<<" +- "<<L_info.Retrieve_Za("B",0).second<<endl;
  cout<<"Za0 branch 1: "<<Boot_ave(Za1[4])<<" +- "<<Rerr*Boot_err(Za1[4])<<" Original from branch nr. 0: "<<L_info.Retrieve_Za("B", 4).first<<" +- "<<L_info.Retrieve_Za("B",4).second<<endl;
  cout<<"beta = 2.10"<<endl;
  cout<<"Za0 branch 0: "<<Boot_ave(Za2[0])<<" +- "<<Rerr*Boot_err(Za2[0])<<" Original from branch nr. 0: "<<L_info.Retrieve_Za("D", 0).first<<" +- "<<L_info.Retrieve_Za("D",0).second<<endl;
  cout<<"Za0 branch 1: "<<Boot_ave(Za2[4])<<" +- "<<Rerr*Boot_err(Za2[4])<<" Original from branch nr. 0: "<<L_info.Retrieve_Za("D", 4).first<<" +- "<<L_info.Retrieve_Za("D",4).second<<endl;


  //print infos to file

  ofstream Print_Info("../data/Mpi_twisted_ov_X/"+print_path+"/fit_result.dat", ofstream::app);
  double ch2_ave= accumulate(chi2_MV.begin(),chi2_MV.end(),0.0)/Nbranches;
  double Akaike_weight = exp(-(ch2_ave/2.0 + Npars));
  Print_Info<< Akaike_weight<<"  "<<Boot_ave(Physical_point)<<" "<<Boot_err(Physical_point,Rerr,1)<<" "<<ch2_ave<<"  "<<Npars<<endl;
  Print_Info.close();
  

    
  return;
}
 
					
  






















