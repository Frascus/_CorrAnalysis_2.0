#include "../include/tau_decay.h"
#include "Corr_analysis.h"

// Use (void) to silence unused warnings.
#define assertm(exp, msg) assert(((void)msg, exp))

const double MPiPhys=0.135;
const double alpha = 1.0/137.04;
const double e2 = alpha*4.0*M_PI;
const bool UseJack=1;
const int Njacks=75;
const int Nboots=800;
const double qu= 2.0/3.0;
const double qd= -1.0/3.0;
const double qs= qd;
const double qc= qu;
const double ln2_10=3.32192809489;
const double fm_to_inv_Gev= 1.0/0.197327;
const int prec = 128;
const double Nc=3;
bool tau_verbosity_lev=1;
const double m_mu= 0.10565837; // [ GeV ]
const double GF= 1.1663787*1e-5; //[GeV^-2]
//CKM matrix elements
const double Vud = 0.97373;
const double m_tau= 1.77686;
const double m_Jpsi= 3.0969;
const double m_Jpsi_err= 0.001;
const double m_phi= 1.019461;
const double m_phi_err= 0.000016;
const double m_etac= 2.9839;
const double m_etac_err = 0.004;
const double m_etas = 0.68989;
const double m_etas_err= 0.00050;
const double m_pi_plus = 0.13957039;
const double E0_l = 0.05*m_tau;
Vfloat sigma_list;
const double C_V = 2*M_PI/(pow(m_tau,3));
const double GAMMA_FACT= 12*M_PI; //12*M_PI*pow(Vud*GF,2);
const string MODE="TANT";
bool Use_t_up_to_T_half=true;
const int sm_func_mode= 0;
const string SM_TYPE_0= "KL_"+to_string(sm_func_mode);
const string SM_TYPE_1= "KT_"+to_string(sm_func_mode);
VVfloat covariance_fake;
int Num_LUSCH=3; //17; //3;;
int Nres= 2; //15; //2;
int pts_spline=200;
const double QCD_scale= 0.3*fm_to_inv_Gev;
bool COMPUTE_SPEC_DENS_FREE=false;
bool Skip_spectral_density_analysis=true;
const bool Perform_continuum_extrapolation=true;
bool Use_Customized_plateaux=true;
using namespace std;


double Customized_plateaux_tau_spectre( double alpha, double Emax, string channel, string reg, double s, string Ens ) {

  double Ra0=-1;
  int alpha_m= (int)(alpha+1);

  if(alpha_m < 3) crash("Customized_plateaux_tau_spectre called with alpha = "+to_string_with_precision(alpha,2));
  if( (alpha_m != 3) && (alpha_m != 4) && (alpha_m != 5) ) crash("Customized_plateaux_tau_spectre called with (int)alpha: "+to_string(alpha_m));


  //SMALL    SIGMA 
  if(s < 0.15) {

    if( reg=="tm") {
      if(channel=="Aii") {
	if(Ens == "cB211b.072.64") {  Ra0 =1e6;  }
	else if(Ens == "cB211b.072.96") { Ra0= ((alpha_m==3)?2e5:2e6);   }
	else if(Ens == "cC211a.06.80") { Ra0= 3e7;   }  //THIS IS TRICKY, OLD VALUE WAS Ra0= 1e5
	else if(Ens == "cD211a.054.96") {  Ra0=1e6;  }
	else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
      }
      else if(channel=="Vii") {
	if(Ens == "cB211b.072.64") { Ra0=3e5;   }
	else if(Ens == "cB211b.072.96") { Ra0= ((alpha_m==3)?1e7:2e7);    } //THIS IS TRICKY, OLD VALUE WAS Ra0=1e5
	else if(Ens == "cC211a.06.80") {  Ra0=1e7;  }  //PREVIOUS WAS Ra0=1e5
	else if(Ens == "cD211a.054.96") {  Ra0= ((alpha_m==3)?1e5:1e6);   }
	else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
      }
      else if(channel=="A0") {
	if(Ens == "cB211b.072.64") { Ra0=1e7;   }
	else if(Ens == "cB211b.072.96") { Ra0=1e7;   }
	else if(Ens == "cC211a.06.80") {  Ra0=3e6;  }
	else if(Ens == "cD211a.054.96") { Ra0=9e7;   } //THIS IS TRICKY, OLD VALUE WAS 1e7
	else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
      }
      else crash("In Customized_plateaux_tau_spectre, channel: "+channel+" not yet implemented");

    }
    else if(reg=="OS") {

      if(channel=="Aii") {
	if(Ens == "cB211b.072.64") { Ra0=1e5;   }
	else if(Ens == "cB211b.072.96") { Ra0=1e5;   }
	else if(Ens == "cC211a.06.80") { Ra0=1e7;   }  //PREVIOUS WAS 1e5
	else if(Ens == "cD211a.054.96") {  Ra0=1e7;  } //PREVIOUS WAS 2e6
	else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
      }
      else if(channel=="Vii") {
	if(Ens == "cB211b.072.64") {   Ra0= ((alpha_m==3)?7e5:1e6);   }
	else if(Ens == "cB211b.072.96") {  Ra0= ((alpha_m==3)?1e5:1e6);    }
	else if(Ens == "cC211a.06.80") {   Ra0= ((alpha_m==3)?1e5:1e6);   }
	else if(Ens == "cD211a.054.96") {  Ra0=1e6;  }
	else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
      }
      else if(channel=="A0") {
	if(Ens == "cB211b.072.64") { Ra0=1e6;   }
	else if(Ens == "cB211b.072.96") { Ra0=2e6;   }
	else if(Ens == "cC211a.06.80") { Ra0=1e6;   }
	else if(Ens == "cD211a.054.96") { if( ( (Emax > (5.0-1e-4)) && (Emax < (5.0 +1e-4))) && (alpha_m==3)  ) Ra0= 2e6; else  Ra0=9e6;   } //THIS IS TRICKY, OLD VALUE WAS 1e6
	else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
      }
      else crash("In Customized_plateaux_tau_spectre, channel: "+channel+" not yet implemented");

    }

    else crash("In Customized_plateaux_tau_spectre, reg: "+reg+" not yet implemented");


  }

  //LARGE   SIGMA 
  else {

    if( reg=="tm") {

      if(channel=="Aii") {
	if(Ens == "cB211b.072.64") {  Ra0=1e4;  }
	else if(Ens == "cB211b.072.96") { Ra0=1e4;   }
	else if(Ens == "cC211a.06.80") {  Ra0=1e4;  }
	else if(Ens == "cD211a.054.96") {  Ra0=1e4;  }
	else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
      }
      else if(channel=="Vii") {
	if(Ens == "cB211b.072.64") { Ra0=3e5;   }
	else if(Ens == "cB211b.072.96") { Ra0=1e5;   }
	else if(Ens == "cC211a.06.80") {  Ra0=1e4;  }
	else if(Ens == "cD211a.054.96") {  Ra0=1e4;  }
	else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
      }
      else if(channel=="A0") {
	if(Ens == "cB211b.072.64") { Ra0=1e7;   }
	else if(Ens == "cB211b.072.96") { Ra0=1e7;   }
	else if(Ens == "cC211a.06.80") { Ra0=4e6;   }
	else if(Ens == "cD211a.054.96") { Ra0=1e7;   }
	else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
      }
      else crash("In Customized_plateaux_tau_spectre, channel: "+channel+" not yet implemented");
       
    }
     
    else if(reg=="OS") {
      if(channel=="Aii") {
	if(Ens == "cB211b.072.64") {  Ra0=1e5;  }
	else if(Ens == "cB211b.072.96") { Ra0=1e4;   }
	else if(Ens == "cC211a.06.80") { Ra0=1e5;   }
	else if(Ens == "cD211a.054.96") { Ra0=1e5;   }
	else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
      }
      else if(channel=="Vii") {
	if(Ens == "cB211b.072.64") {  Ra0=1e5;  }
	else if(Ens == "cB211b.072.96") { Ra0=1e5;   }
	else if(Ens == "cC211a.06.80") { Ra0=1e5;   }
	else if(Ens == "cD211a.054.96") { Ra0=1e5;   }
	else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
      }
      else if(channel=="A0") {
	if(Ens == "cB211b.072.64") { Ra0=1e6;   }
	else if(Ens == "cB211b.072.96") { Ra0=1e6;   }
	else if(Ens == "cC211a.06.80") { Ra0=1e6;   }
	else if(Ens == "cD211a.054.96") { Ra0=1e6;   }
	else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
      }
      else crash("In Customized_plateaux_tau_spectre, channel: "+channel+" not yet implemented");
      
    }

    else crash("In Customized_plateaux_tau_spectre, reg: "+reg+" not yet implemented");




  }

  assertm(Ra0 > 0, "Assertion Ra0 > 0 failed");
  return Ra0;


}



void get_sigma_list() {

  bool test=false;
  if(test) { sigma_list.push_back(0.05); return;}


  double s_max= 0.2;
  double s_min= 0.004420;

  sigma_list= {0.004, 0.02, 0.03, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.20};
  return;

  //sigma_list.push_back(s_max);
  //double s=s_max;
  //while (s>= s_min +1.0e-6) { s /= pow(2,0.25); sigma_list.push_back(s);}

  return;

}



void tau_decay_analysis() {



  get_sigma_list();

  Vfloat betas({ 1.99, 2.99, 3.99, 4.99, 2.99, 3.99, 2.99, 4.99});
  Vfloat Emax_list({4.0, 4.0, 4.0, 4.0, 5.0, 5.0, 6.0, 5.0 });
  vector<bool> Is_Emax_Finite({0,1,1,1,1,1,1,1});
  
  

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

  
   


  //Init LL_functions;
  //find first  zeros of the Lusher functions
  Vfloat Luscher_zeroes;
  Zeta_function_zeroes(Num_LUSCH, Luscher_zeroes);
  

  //############################################INTERPOLATE PHI FUNCTION AND DERIVATIVES#############################
  cout<<"Computing Luscher-zeros"<<flush;

  VVfloat phi_data, phi_der_data;
  Vfloat sx_int;
  Vfloat sx_der, dx_der;
  Vfloat Dz;
  
  for(int L_zero=0;L_zero<Nres+1;L_zero++) {
    cout<<"."<<flush;
    double sx, dx;
    //interpolating between the Luscher_zero[L_zero-1] and Luscher_zero[L_zero];
    if(L_zero==0) { sx_int.push_back(0.0); sx=0.0;}
    else {sx=Luscher_zeroes[L_zero-1];  sx_int.push_back(sx);}
    dx= Luscher_zeroes[L_zero];
    phi_data.resize(L_zero+1);
    phi_der_data.resize(L_zero+1);
    phi_data[L_zero].push_back(L_zero==0?0.0:-M_PI/2.0);
    //divide interval into thousand points;
    double dz = (dx-sx)/pts_spline;
    Dz.push_back(dz);


    for(int istep=1;istep<=pts_spline-1;istep++) { double pt= sx+dz*istep; phi_data[L_zero].push_back( phi(sqrt(pt)));}

    phi_data[L_zero].push_back(M_PI/2.0);
    double sx_der_loc =  phi_der_for_back(sqrt(sx)+1e-14, 1);
    double dx_der_loc =  phi_der_for_back(sqrt(dx)-1e-14, -1);
    sx_der.push_back(sx_der_loc);
    dx_der.push_back(dx_der_loc);

    phi_der_data[L_zero].push_back(sx_der_loc);
    for(int istep=1;istep<=pts_spline-1;istep++) { double pt= sx+dz*istep; phi_der_data[L_zero].push_back( phi_der(sqrt(pt)));}
    phi_der_data[L_zero].push_back(dx_der_loc);
    
  }
  cout<<"done!"<<endl;



  //###########################################END INTERPOLATION PHI FUNCTION AND DERIVATIVES################################
 
   

  LL_functions LL(phi_data,phi_der_data,sx_der, dx_der, sx_int, Dz, Nres, Luscher_zeroes);

  


  
  for(int i=rank*N/size; i<(rank+1)*N/size;i++) {Compute_tau_decay_width(Is_Emax_Finite[i], Emax_list[i], betas[i], LL); COMPUTE_SPEC_DENS_FREE=false;}

}


void Compute_tau_decay_width(bool Is_Emax_Finite, double Emax, double beta,LL_functions &LL) {

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
  string Is_Spec_free_computed= ((COMPUTE_SPEC_DENS_FREE==0)?"no":"yes");
  string alpha_Emax_Tag= "{"+to_string_with_precision(beta,2)+","+((Is_Emax_Finite==0)?"inf":to_string_with_precision(Emax,1))+"}";


  cout<<"STARTING COMPUTATION OF: {alpha,Emax} : "<<alpha_Emax_Tag<<endl;
  cout<<"COMPUTE FREE SPECTRAL DENSITY: "<<Is_Spec_free_computed<<endl;
  cout<<"Creating output directories...";
  
  cout.precision(5);



  if(COMPUTE_SPEC_DENS_FREE) Get_spec_dens_free({0.00072, 0.00060, 0.00054}, "tau_decay");

 
    


  //create directories

  boost::filesystem::create_directory("../data/tau_decay");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type);
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/light");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/light/continuum");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/light/continuum/fit_func");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/light/continuum/fit_data");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/light/continuum/cov");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/light/continuum/corr");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/light/continuum/AIC");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/light/continuum/fit_data/A0A0");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/light/continuum/fit_data/A0A0/tm");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/light/continuum/fit_data/A0A0/OS");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/light/continuum/fit_data/AkAk");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/light/continuum/fit_data/AkAk/tm");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/light/continuum/fit_data/AkAk/OS");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/light/continuum/fit_data/VkVk");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/light/continuum/fit_data/VkVk/tm");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/light/continuum/fit_data/VkVk/OS");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/light/continuum/fit_data/tot");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/light/continuum/fit_data/tot/tm");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/light/continuum/fit_data/tot/OS");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/light/continuum/fit_data/VA");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/light/continuum/fit_data/VA/tm");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/light/continuum/fit_data/VA/OS");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/light/continuum/fit_data/VMA");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/light/continuum/fit_data/VMA/tm");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/light/continuum/fit_data/VMA/OS");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/light/continuum/fit_data/AX");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/light/continuum/fit_data/AX/tm");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/light/continuum/fit_data/AX/OS");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/light/continuum/fit_data/T");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/light/continuum/fit_data/T/tm");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/light/continuum/fit_data/T/OS");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/light/jackknife");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/light/jackknife/tm");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/light/jackknife/tm/A0A0");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/light/jackknife/tm/AkAk");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/light/jackknife/tm/VkVk");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/light/jackknife/OS");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/light/jackknife/OS/A0A0");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/light/jackknife/OS/AkAk");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/light/jackknife/OS/VkVk");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/light/Br");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/light/corr");
  boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/light/covariance");

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

 
   data_t Vk_data_tm, V0_data_tm, Ak_data_tm, A0_data_tm, P5_data_tm;
   data_t Vk_data_OS, V0_data_OS, Ak_data_OS, A0_data_OS;
  

  
  //light 
  //tm
  Vk_data_tm.Read("../R_ratio_data/light", "mes_contr_2pts_ll_1", "VKVK", Sort_light_confs);
  V0_data_tm.Read("../R_ratio_data/light", "mes_contr_2pts_ll_1", "V0V0", Sort_light_confs);
  Ak_data_tm.Read("../R_ratio_data/light", "mes_contr_2pts_ll_1", "AKAK", Sort_light_confs);
  A0_data_tm.Read("../R_ratio_data/light", "mes_contr_2pts_ll_1", "A0A0", Sort_light_confs);
  P5_data_tm.Read("../R_ratio_data/light", "mes_contr_2pts_ll_1", "P5P5", Sort_light_confs);
  //OS
  Vk_data_OS.Read("../R_ratio_data/light", "mes_contr_2pts_ll_2", "VKVK", Sort_light_confs);
  V0_data_OS.Read("../R_ratio_data/light", "mes_contr_2pts_ll_2", "V0V0", Sort_light_confs);
  Ak_data_OS.Read("../R_ratio_data/light", "mes_contr_2pts_ll_2", "AKAK", Sort_light_confs);
  A0_data_OS.Read("../R_ratio_data/light", "mes_contr_2pts_ll_2", "A0A0", Sort_light_confs);

  
  data_t Vk_ph_data_tm, V0_ph_data_tm, Ak_ph_data_tm, A0_ph_data_tm;
  data_t Vk_ph_data_OS, V0_ph_data_OS, Ak_ph_data_OS, A0_ph_data_OS;


  //light 
  //tm
  Vk_ph_data_tm.Read("../tau_decay_data/light_mass_correction/phys", "mes_contr_2pts_ll_1", "VKVK", Sort_light_confs);
  V0_ph_data_tm.Read("../tau_decay_data/light_mass_correction/phys", "mes_contr_2pts_ll_1", "V0V0", Sort_light_confs);
  Ak_ph_data_tm.Read("../tau_decay_data/light_mass_correction/phys", "mes_contr_2pts_ll_1", "AKAK", Sort_light_confs);
  A0_ph_data_tm.Read("../tau_decay_data/light_mass_correction/phys", "mes_contr_2pts_ll_1", "A0A0", Sort_light_confs);
  //OS
  Vk_ph_data_OS.Read("../tau_decay_data/light_mass_correction/phys", "mes_contr_2pts_ll_2", "VKVK", Sort_light_confs);
  V0_ph_data_OS.Read("../tau_decay_data/light_mass_correction/phys", "mes_contr_2pts_ll_2", "V0V0", Sort_light_confs);
  Ak_ph_data_OS.Read("../tau_decay_data/light_mass_correction/phys", "mes_contr_2pts_ll_2", "AKAK", Sort_light_confs);
  A0_ph_data_OS.Read("../tau_decay_data/light_mass_correction/phys", "mes_contr_2pts_ll_2", "A0A0", Sort_light_confs);
  

  data_t Vk_uni_data_tm, V0_uni_data_tm, Ak_uni_data_tm, A0_uni_data_tm;
  data_t Vk_uni_data_OS, V0_uni_data_OS, Ak_uni_data_OS, A0_uni_data_OS;


  //light 
  //tm
  Vk_uni_data_tm.Read("../tau_decay_data/light_mass_correction/unitary", "mes_contr_2pts_ll_1", "VKVK", Sort_light_confs);
  V0_uni_data_tm.Read("../tau_decay_data/light_mass_correction/unitary", "mes_contr_2pts_ll_1", "V0V0", Sort_light_confs);
  Ak_uni_data_tm.Read("../tau_decay_data/light_mass_correction/unitary", "mes_contr_2pts_ll_1", "AKAK", Sort_light_confs);
  A0_uni_data_tm.Read("../tau_decay_data/light_mass_correction/unitary", "mes_contr_2pts_ll_1", "A0A0", Sort_light_confs);
  //OS
  Vk_uni_data_OS.Read("../tau_decay_data/light_mass_correction/unitary", "mes_contr_2pts_ll_2", "VKVK", Sort_light_confs);
  V0_uni_data_OS.Read("../tau_decay_data/light_mass_correction/unitary", "mes_contr_2pts_ll_2", "V0V0", Sort_light_confs);
  Ak_uni_data_OS.Read("../tau_decay_data/light_mass_correction/unitary", "mes_contr_2pts_ll_2", "AKAK", Sort_light_confs);
  A0_uni_data_OS.Read("../tau_decay_data/light_mass_correction/unitary", "mes_contr_2pts_ll_2", "A0A0", Sort_light_confs);



  
  



  //############################################################################################
  //generate fake jack_distr for lattice spacing a_A a_B, a_C, a_D and RENORMALIZATION CONSTANT
  GaussianMersenne GM(36551294);
  LatticeInfo a_info;
  distr_t a_A(UseJack), a_B(UseJack), a_C(UseJack), a_D(UseJack);
  distr_t ZV_A(UseJack), ZV_B(UseJack), ZV_C(UseJack), ZV_D(UseJack);
  distr_t ZA_A(UseJack), ZA_B(UseJack), ZA_C(UseJack), ZA_D(UseJack);
  double a_A_ave, a_A_err, a_B_ave, a_B_err, a_C_ave, a_C_err, a_D_ave, a_D_err;
  double ZV_A_ave, ZV_A_err, ZV_B_ave, ZV_B_err, ZV_C_ave, ZV_C_err, ZV_D_ave, ZV_D_err;
  double ZA_A_ave, ZA_A_err, ZA_B_ave, ZA_B_err, ZA_C_ave, ZA_C_err, ZA_D_ave, ZA_D_err;
  a_info.LatInfo_new_ens("cA211a.53.24");
  a_A_ave= a_info.a_from_afp;
  a_A_err= a_info.a_from_afp_err;
  ZA_A_ave = a_info.Za_WI_strange;
  ZA_A_err = a_info.Za_WI_strange_err;
  ZV_A_ave = a_info.Zv_WI_strange;
  ZV_A_err = a_info.Zv_WI_strange_err;
  a_info.LatInfo_new_ens("cB211b.072.64");
  a_B_ave= a_info.a_from_afp;
  a_B_err= a_info.a_from_afp_err;
  ZA_B_ave = a_info.Za_WI_strange;
  ZA_B_err = a_info.Za_WI_strange_err;
  ZV_B_ave = a_info.Zv_WI_strange;
  ZV_B_err = a_info.Zv_WI_strange_err;
  a_info.LatInfo_new_ens("cC211a.06.80");
  a_C_ave= a_info.a_from_afp;
  a_C_err= a_info.a_from_afp_err;
  ZA_C_ave = a_info.Za_WI_strange;
  ZA_C_err = a_info.Za_WI_strange_err;
  ZV_C_ave = a_info.Zv_WI_strange;
  ZV_C_err = a_info.Zv_WI_strange_err;
  a_info.LatInfo_new_ens("cD211a.054.96");
  a_D_ave= a_info.a_from_afp;
  a_D_err= a_info.a_from_afp_err;
  ZA_D_ave = a_info.Za_WI_strange;
  ZA_D_err = a_info.Za_WI_strange_err;
  ZV_D_ave = a_info.Zv_WI_strange;
  ZV_D_err = a_info.Zv_WI_strange_err;
  
  if(UseJack)  { for(int ijack=0;ijack<Njacks;ijack++) {
      a_A.distr.push_back( fm_to_inv_Gev*( a_A_ave + GM()*a_A_err*(1.0/sqrt(Njacks-1.0))));
      a_B.distr.push_back( fm_to_inv_Gev*( a_B_ave + GM()*a_B_err*(1.0/sqrt(Njacks-1.0))));
      a_C.distr.push_back( fm_to_inv_Gev*( a_C_ave + GM()*a_C_err*(1.0/sqrt(Njacks-1.0))));
      a_D.distr.push_back( fm_to_inv_Gev*( a_D_ave + GM()*a_D_err*(1.0/sqrt(Njacks-1.0))));
      ZA_A.distr.push_back(  ZA_A_ave + GM()*ZA_A_err*(1.0/sqrt(Njacks -1.0)));
      ZV_A.distr.push_back(  ZV_A_ave + GM()*ZV_A_err*(1.0/sqrt(Njacks -1.0)));
      ZA_B.distr.push_back(  ZA_B_ave + GM()*ZA_B_err*(1.0/sqrt(Njacks -1.0)));
      ZV_B.distr.push_back(  ZV_B_ave + GM()*ZV_B_err*(1.0/sqrt(Njacks -1.0)));
      ZA_C.distr.push_back(  ZA_C_ave + GM()*ZA_C_err*(1.0/sqrt(Njacks -1.0)));
      ZV_C.distr.push_back(  ZV_C_ave + GM()*ZV_C_err*(1.0/sqrt(Njacks -1.0)));
      ZA_D.distr.push_back(  ZA_D_ave + GM()*ZA_D_err*(1.0/sqrt(Njacks -1.0)));
      ZV_D.distr.push_back(  ZV_D_ave + GM()*ZV_D_err*(1.0/sqrt(Njacks -1.0)));
      
    }
  }
  else {
    for (int iboot=0; iboot<Nboots;iboot++) {
      a_A.distr.push_back( fm_to_inv_Gev*( a_A_ave + GM()*a_A_err));
      a_B.distr.push_back( fm_to_inv_Gev*( a_B_ave + GM()*a_B_err));
      a_C.distr.push_back( fm_to_inv_Gev*( a_C_ave + GM()*a_C_err));
      a_D.distr.push_back( fm_to_inv_Gev*( a_D_ave + GM()*a_D_err));
      ZA_A.distr.push_back(  ZA_A_ave + GM()*ZA_A_err);
      ZV_A.distr.push_back(  ZV_A_ave + GM()*ZV_A_err);
      ZA_B.distr.push_back(  ZA_B_ave + GM()*ZA_B_err);
      ZV_B.distr.push_back(  ZV_B_ave + GM()*ZV_B_err);
      ZA_C.distr.push_back(  ZA_C_ave + GM()*ZA_C_err);
      ZV_C.distr.push_back(  ZV_C_ave + GM()*ZV_C_err);
      ZA_D.distr.push_back(  ZA_D_ave + GM()*ZA_D_err);
      ZV_D.distr.push_back(  ZV_D_ave + GM()*ZV_D_err);
      
    }
  }


  //############################################################################################




  int Nens= Vk_data_tm.size;

  //resize vector with systematic errors
  vector<distr_t_list> Br_tau_tm, Br_tau_OS;
  vector<distr_t_list> Br_A0_tau_tm, Br_Aii_tau_tm, Br_Vii_tau_tm;
  vector<distr_t_list> Br_A0_tau_OS, Br_Aii_tau_OS, Br_Vii_tau_OS;

  distr_t_list RESC_FP(UseJack);
  
  VVfloat syst_per_ens_tm_A0(Nens);
  VVfloat syst_per_ens_tm_Ak(Nens);
  VVfloat syst_per_ens_tm_Vk(Nens);
  VVfloat syst_per_ens_OS_A0(Nens);
  VVfloat syst_per_ens_OS_Ak(Nens);
  VVfloat syst_per_ens_OS_Vk(Nens);
  VVfloat syst_per_ens_tm(Nens);
  VVfloat syst_per_ens_OS(Nens);


  for(int iens=0; iens<Nens; iens++) {
    Br_tau_tm.emplace_back( UseJack, sigma_list.size());
    Br_A0_tau_tm.emplace_back( UseJack, sigma_list.size());
    Br_Aii_tau_tm.emplace_back( UseJack, sigma_list.size());
    Br_Vii_tau_tm.emplace_back( UseJack, sigma_list.size());
    
    Br_tau_OS.emplace_back( UseJack, sigma_list.size());
    Br_A0_tau_OS.emplace_back( UseJack, sigma_list.size());
    Br_Aii_tau_OS.emplace_back( UseJack, sigma_list.size());
    Br_Vii_tau_OS.emplace_back( UseJack, sigma_list.size());

    syst_per_ens_tm_A0[iens].resize(sigma_list.size());
    syst_per_ens_tm_Ak[iens].resize(sigma_list.size());
    syst_per_ens_tm_Vk[iens].resize(sigma_list.size());
    syst_per_ens_tm[iens].resize(sigma_list.size());

    syst_per_ens_OS_A0[iens].resize(sigma_list.size());
    syst_per_ens_OS_Ak[iens].resize(sigma_list.size());
    syst_per_ens_OS_Vk[iens].resize(sigma_list.size());
    syst_per_ens_OS[iens].resize(sigma_list.size());
       
    

  }

 
  
  
 
  if(!Skip_spectral_density_analysis) {


  //loop over the ensembles
  for(int iens=0; iens<Nens;iens++) {


    //print number of gauge configurations
    cout<<"################### Analyzing Ensemble: "<<Vk_data_tm.Tag[iens]<<endl;
    cout<<"Nconfs : "<<Vk_data_OS.Nconfs[iens]<<endl;

    auto SQRT = [](double x) {return sqrt(x);};

    LatticeInfo L_info;
    
    L_info.LatInfo_new_ens(Vk_data_tm.Tag[iens]);

    double aml= L_info.ml;


    //Read perturbative data for OS and tm
    Vfloat Spec_tm = Read_From_File("../data/tau_decay/spec_dens_free/tm/am_"+to_string_with_precision(L_info.ml,5), 2, 4);
    Vfloat Spec_OS = Read_From_File("../data/tau_decay/spec_dens_free/OS/am_"+to_string_with_precision(L_info.ml,5), 2, 4);
    Vfloat Ergs_pert = Read_From_File("../data/tau_decay/spec_dens_free/tm/am_"+to_string_with_precision(L_info.ml,5), 1, 4);

    cout<<"perturbative data for Ensemble: "<<Vk_data_tm.Tag[iens]<<" READ! "<<endl;
    

    //interpolate perturbative data
    boost::math::interpolators::cardinal_cubic_b_spline<double> F_boost_tm(Spec_tm.begin(), Spec_tm.end(), Ergs_pert[0], 2.0*Ergs_pert[0]);
    boost::math::interpolators::cardinal_cubic_b_spline<double> F_boost_OS(Spec_OS.begin(), Spec_OS.end(), Ergs_pert[0], 2.0*Ergs_pert[0]);

    cout<<"Cubic spline for perturbative data for Ensemble: "<<Vk_data_tm.Tag[iens]<<" produced! "<<endl;

    auto F_free_tm = [&F_boost_tm](double E) { return F_boost_tm(E);};
    auto F_free_OS = [&F_boost_OS](double E) { return F_boost_OS(E);};
    
     
    CorrAnalysis Corr(UseJack, Njacks,Nboots);
    CorrAnalysis Corr_block_1(0, Vk_data_tm.Nconfs[iens],Nboots, iens);
    Corr_block_1.Nt= Vk_data_tm.nrows[iens];
    Corr.Nt = Vk_data_tm.nrows[iens];
    int T = Corr.Nt;
    

   
       
    //get lattice spacing
    distr_t a_distr(UseJack);
    distr_t Zv(UseJack), Za(UseJack);
    double Mpi=0.0;
    double Mpi_err=0.0;
    double fpi=0.0;
    double Mpi_OS=0.0;
    double fpi_OS=0.0;
    double dm=0.0;
    if(Vk_data_tm.Tag[iens].substr(1,1)=="B") {a_distr=a_B; Zv = ZV_B; Za = ZA_B; Mpi=0.05653312833; Mpi_err=1.430196186e-05; fpi=0.05278353769; Mpi_OS=0.1203989717; dm= (0.00072-0.0006675);}
    else if(Vk_data_tm.Tag[iens].substr(1,1)=="C") {a_distr=a_C; Zv = ZV_C; Za = ZA_C; Mpi=0.04722061628; Mpi_err=3.492993579e-05; fpi=0.0450246; Mpi_OS=0.08597942324; dm = (0.00060-0.000585);}
    else if(Vk_data_tm.Tag[iens].substr(1,1)=="D") {a_distr=a_D; Zv = ZV_D; Za = ZA_D; Mpi=0.04062107883; Mpi_err= 2.973916243e-05; fpi=0.03766423429; Mpi_OS=0.06064150466; dm = (0.00054 -0.0004964);}
    else crash("lattice spacing distribution for Ens: "+Vk_data_tm.Tag[iens]+" not found");
    fpi_OS= fpi*L_info.Za_WI_strange/L_info.Za_WI;

  

   
    //jack distr for Mpi
    distr_t Mpi_distr(UseJack);
    for(int ij=0;ij<Njacks;ij++) Mpi_distr.distr.push_back( Mpi + GM()*Mpi_err/sqrt(Njacks-1.0));

 

    //tm
    distr_t_list Vk_tm_distr, V0_tm_distr, Ak_tm_distr, A0_tm_distr;
    //tm block1
    distr_t_list Vk_tm_block_1_distr, V0_tm_block_1_distr, Ak_tm_block_1_distr, A0_tm_block_1_distr;
    //OS
    distr_t_list Vk_OS_distr, V0_OS_distr, Ak_OS_distr, A0_OS_distr;
    //OS block1
    distr_t_list Vk_OS_block_1_distr, V0_OS_block_1_distr, Ak_OS_block_1_distr, A0_OS_block_1_distr;

    distr_t_list P5_tm_distr;


    
   
    //light-tm sector
    Vk_tm_distr = Corr.corr_t( summ_master(Vk_data_tm.col(0)[iens], Vk_ph_data_tm.col(0)[iens], Multiply_Vvector_by_scalar(Vk_uni_data_tm.col(0)[iens],-1.0)) , "../data/tau_decay/"+Tag_reco_type+"/light/corr/Vk_tm_"+Vk_data_tm.Tag[iens]+".dat");
    V0_tm_distr = Corr.corr_t(summ_master(V0_data_tm.col(0)[iens], V0_ph_data_tm.col(0)[iens], Multiply_Vvector_by_scalar(V0_uni_data_tm.col(0)[iens],-1.0)), "../data/tau_decay/"+Tag_reco_type+"/light/corr/V0_tm_"+Vk_data_tm.Tag[iens]+".dat");
    Ak_tm_distr = Corr.corr_t(summ_master(Ak_data_tm.col(0)[iens], Ak_ph_data_tm.col(0)[iens], Multiply_Vvector_by_scalar(Ak_uni_data_tm.col(0)[iens],-1.0)), "../data/tau_decay/"+Tag_reco_type+"/light/corr/Ak_tm_"+Vk_data_tm.Tag[iens]+".dat");
    A0_tm_distr = Corr.corr_t(summ_master(A0_data_tm.col(0)[iens], A0_ph_data_tm.col(0)[iens], Multiply_Vvector_by_scalar(A0_uni_data_tm.col(0)[iens],-1.0)), "../data/tau_decay/"+Tag_reco_type+"/light/corr/A0_tm_"+Vk_data_tm.Tag[iens]+".dat");
    P5_tm_distr = Corr.corr_t(P5_data_tm.col(0)[iens], "../data/tau_decay/"+Tag_reco_type+"/light/corr/P5_tm_"+Vk_data_tm.Tag[iens]+".dat");

    //light-OS sector
    Vk_OS_distr = Corr.corr_t( summ_master(Vk_data_OS.col(0)[iens], Vk_ph_data_OS.col(0)[iens], Multiply_Vvector_by_scalar(Vk_uni_data_OS.col(0)[iens],-1.0)) , "../data/tau_decay/"+Tag_reco_type+"/light/corr/Vk_OS_"+Vk_data_OS.Tag[iens]+".dat");
    V0_OS_distr = Corr.corr_t(summ_master(V0_data_OS.col(0)[iens], V0_ph_data_OS.col(0)[iens], Multiply_Vvector_by_scalar(V0_uni_data_OS.col(0)[iens],-1.0)), "../data/tau_decay/"+Tag_reco_type+"/light/corr/V0_OS_"+Vk_data_OS.Tag[iens]+".dat");
    Ak_OS_distr = Corr.corr_t(summ_master(Ak_data_OS.col(0)[iens], Ak_ph_data_OS.col(0)[iens], Multiply_Vvector_by_scalar(Ak_uni_data_OS.col(0)[iens],-1.0)), "../data/tau_decay/"+Tag_reco_type+"/light/corr/Ak_OS_"+Vk_data_OS.Tag[iens]+".dat");
    A0_OS_distr = Corr.corr_t(summ_master(A0_data_OS.col(0)[iens], A0_ph_data_OS.col(0)[iens], Multiply_Vvector_by_scalar(A0_uni_data_OS.col(0)[iens],-1.0)), "../data/tau_decay/"+Tag_reco_type+"/light/corr/A0_OS_"+Vk_data_OS.Tag[iens]+".dat");
  
    
    //Vk_OS_distr = Corr.corr_t(Vk_data_OS.col(0)[iens], "../data/tau_decay/"+Tag_reco_type+"/light/corr/Vk_OS_"+Vk_data_tm.Tag[iens]+".dat");
    //V0_OS_distr = Corr.corr_t(V0_data_OS.col(0)[iens], "../data/tau_decay/"+Tag_reco_type+"/light/corr/V0_OS_"+Vk_data_tm.Tag[iens]+".dat");
    //Ak_OS_distr = Corr.corr_t(Ak_data_OS.col(0)[iens], "../data/tau_decay/"+Tag_reco_type+"/light/corr/Ak_OS_"+Vk_data_tm.Tag[iens]+".dat");
    //A0_OS_distr = Corr.corr_t(A0_data_OS.col(0)[iens], "../data/tau_decay/"+Tag_reco_type+"/light/corr/A0_OS_"+Vk_data_tm.Tag[iens]+".dat");


    //analyze data with Njacks=Nconfs
    //light-tm sector
    Vk_tm_block_1_distr = Corr_block_1.corr_t(summ_master(Vk_data_tm.col(0)[iens], Vk_ph_data_tm.col(0)[iens], Multiply_Vvector_by_scalar(Vk_uni_data_tm.col(0)[iens],-1.0)), "");
    V0_tm_block_1_distr = Corr_block_1.corr_t(summ_master(V0_data_tm.col(0)[iens], V0_ph_data_tm.col(0)[iens], Multiply_Vvector_by_scalar(V0_uni_data_tm.col(0)[iens],-1.0)), "");
    Ak_tm_block_1_distr = Corr_block_1.corr_t(summ_master(Ak_data_tm.col(0)[iens], Ak_ph_data_tm.col(0)[iens], Multiply_Vvector_by_scalar(Ak_uni_data_tm.col(0)[iens],-1.0)), "");
    A0_tm_block_1_distr = Corr_block_1.corr_t(summ_master(A0_data_tm.col(0)[iens], A0_ph_data_tm.col(0)[iens], Multiply_Vvector_by_scalar(A0_uni_data_tm.col(0)[iens],-1.0)), "");

    //light-OS sector
    Vk_OS_block_1_distr = Corr_block_1.corr_t(summ_master(Vk_data_OS.col(0)[iens], Vk_ph_data_OS.col(0)[iens], Multiply_Vvector_by_scalar(Vk_uni_data_OS.col(0)[iens],-1.0)), "");
    V0_OS_block_1_distr = Corr_block_1.corr_t(summ_master(V0_data_OS.col(0)[iens], V0_ph_data_OS.col(0)[iens], Multiply_Vvector_by_scalar(V0_uni_data_OS.col(0)[iens],-1.0)), "");
    Ak_OS_block_1_distr = Corr_block_1.corr_t(summ_master(Ak_data_OS.col(0)[iens], Ak_ph_data_OS.col(0)[iens], Multiply_Vvector_by_scalar(Ak_uni_data_OS.col(0)[iens],-1.0)), "");
    A0_OS_block_1_distr = Corr_block_1.corr_t(summ_master(A0_data_OS.col(0)[iens], A0_ph_data_OS.col(0)[iens], Multiply_Vvector_by_scalar(A0_uni_data_OS.col(0)[iens],-1.0)), "");



    int Tmin_P5=0;
    int Tmax_P5=0;
    if( P5_data_tm.Tag[iens] =="cB211b.072.96")     { Tmin_P5=30; Tmax_P5=70;}
    else if(P5_data_tm.Tag[iens] =="cB211b.072.64") { Tmin_P5=27; Tmax_P5=50;}
    else if(P5_data_tm.Tag[iens].substr(1,1)=="C")  { Tmin_P5=40; Tmax_P5=60;}
    else if(P5_data_tm.Tag[iens].substr(1,1)=="D")  { Tmin_P5=41; Tmax_P5=80;}
    else crash("Cannot recognize the ensemble: "+P5_data_tm.Tag[iens]+" in assigning Tmin_P5,Tmax_P5 for ensemble: ");

    Corr.Tmin = Tmin_P5; Corr.Tmax= Tmax_P5;

    distr_t MP= Corr.Fit_distr( Corr.effective_mass_t( P5_tm_distr, ""));
    distr_t FP= Corr.Fit_distr( Corr.decay_constant_t( pow(2*aml, 2)*P5_tm_distr, ""));

    double L= L_info.L;
    auto LOG = [](double x) { return log(x);};
    double fp_phys= 0.1304;
    double Mp_phys= 0.1350;
    double csi_phys= pow(Mp_phys/(4.0*M_PI*fp_phys),2);
    double l1ph= -0.4; //-0.4
    double l2ph= 4.3; //4.3
    double l3ph= 3.2; //3.2
    double l4ph= 4.4; //4.4
    double s0= 2.0-M_PI/2.0;
    double s1 = M_PI/4.0 - 0.5;
    double s2 = 0.5 - M_PI/8.0;
    double s3 = 3.0*M_PI/16.0 - 0.5;
    distr_t csi_L = MP*MP/(pow(4.0*M_PI,2)*FP*FP);
    distr_t g1 = distr_t::f_of_distr(g1_l, MP*L);
    distr_t g2 = distr_t::f_of_distr(g2_l, MP*L);
    distr_t log_l = log(csi_phys) - distr_t::f_of_distr(LOG, csi_L);
    distr_t FP_CDH = FP/(1.0 -2.0*csi_L*g1 +2.0*csi_L*csi_L*( (Cf1(l1ph,l2ph,l3ph,l4ph) + Sf1(s0,s1,s2,s3) + Cf1_log()*log_l)*g1 + (Cf2(l1ph,l2ph,l3ph,l4ph) + Sf2(s0,s1,s2,s3) + Cf2_log()*log_l)*g2));



    distr_t resc_GeV = C_V*GAMMA_FACT/(a_distr.ave()*a_distr*a_distr);

    distr_t resc_fpi = a_distr*a_distr*(0.13041*0.13041)/(FP_CDH*FP_CDH);
    
    RESC_FP.distr_list.push_back(resc_fpi);
    

  
    //print covariance matrix

    Vfloat cov_A0_tm, cov_Ak_tm, cov_Vk_tm, cov_A0_OS, cov_Ak_OS, cov_Vk_OS, TT, RR;
    Vfloat corr_m_A0_tm, corr_m_Ak_tm, corr_m_Vk_tm, corr_m_A0_OS, corr_m_Ak_OS, corr_m_Vk_OS;
    for(int tt=0;tt<Corr.Nt;tt++)
      for(int rr=0;rr<Corr.Nt;rr++) {
	TT.push_back(tt);
	RR.push_back(rr);
	cov_A0_tm.push_back( A0_tm_block_1_distr.distr_list[tt]%A0_tm_block_1_distr.distr_list[rr]);
	cov_Ak_tm.push_back( Ak_tm_block_1_distr.distr_list[tt]%Ak_tm_block_1_distr.distr_list[rr]);
	cov_Vk_tm.push_back( Vk_tm_block_1_distr.distr_list[tt]%Vk_tm_block_1_distr.distr_list[rr]);
	cov_A0_OS.push_back( A0_OS_block_1_distr.distr_list[tt]%A0_OS_block_1_distr.distr_list[rr]);
	cov_Ak_OS.push_back( Ak_OS_block_1_distr.distr_list[tt]%Ak_OS_block_1_distr.distr_list[rr]);
	cov_Vk_OS.push_back( Vk_OS_block_1_distr.distr_list[tt]%Vk_OS_block_1_distr.distr_list[rr]);


	corr_m_A0_tm.push_back( (A0_tm_block_1_distr.distr_list[tt]%A0_tm_block_1_distr.distr_list[rr])/(A0_tm_block_1_distr.err(tt)*A0_tm_block_1_distr.err(rr)));
	corr_m_Ak_tm.push_back( (Ak_tm_block_1_distr.distr_list[tt]%Ak_tm_block_1_distr.distr_list[rr])/( Ak_tm_block_1_distr.err(tt)*Ak_tm_block_1_distr.err(rr)));
	corr_m_Vk_tm.push_back( (Vk_tm_block_1_distr.distr_list[tt]%Vk_tm_block_1_distr.distr_list[rr])/( Vk_tm_block_1_distr.err(tt)*Vk_tm_block_1_distr.err(rr)));
	corr_m_A0_OS.push_back( (A0_OS_block_1_distr.distr_list[tt]%A0_OS_block_1_distr.distr_list[rr])/( A0_OS_block_1_distr.err(tt)*A0_OS_block_1_distr.err(rr)));
	corr_m_Ak_OS.push_back( (Ak_OS_block_1_distr.distr_list[tt]%Ak_OS_block_1_distr.distr_list[rr])/( Ak_OS_block_1_distr.err(tt)*Ak_OS_block_1_distr.err(rr)));
	corr_m_Vk_OS.push_back( (Vk_OS_block_1_distr.distr_list[tt]%Vk_OS_block_1_distr.distr_list[rr])/( Vk_OS_block_1_distr.err(tt)*Vk_OS_block_1_distr.err(rr)));
      }

    Print_To_File({}, {TT,RR,cov_A0_tm, corr_m_A0_tm}, "../data/tau_decay/"+Tag_reco_type+"/light/covariance/A0_tm_"+Vk_data_tm.Tag[iens]+".dat", "", "");
    Print_To_File({}, {TT,RR,cov_Ak_tm, corr_m_Ak_tm}, "../data/tau_decay/"+Tag_reco_type+"/light/covariance/Aii_tm_"+Vk_data_tm.Tag[iens]+".dat", "", "");
    Print_To_File({}, {TT,RR,cov_Vk_tm, corr_m_Vk_tm}, "../data/tau_decay/"+Tag_reco_type+"/light/covariance/Vii_tm_"+Vk_data_tm.Tag[iens]+".dat", "", "");
    Print_To_File({}, {TT,RR,cov_A0_OS, corr_m_A0_OS}, "../data/tau_decay/"+Tag_reco_type+"/light/covariance/A0_OS_"+Vk_data_tm.Tag[iens]+".dat", "", "");
    Print_To_File({}, {TT,RR,cov_Ak_OS, corr_m_Ak_OS}, "../data/tau_decay/"+Tag_reco_type+"/light/covariance/Aii_OS_"+Vk_data_tm.Tag[iens]+".dat", "", "");
    Print_To_File({}, {TT,RR,cov_Vk_OS, corr_m_Vk_OS}, "../data/tau_decay/"+Tag_reco_type+"/light/covariance/Vii_OS_"+Vk_data_tm.Tag[iens]+".dat", "", "");

       

    distr_t_list A0_tm, Aii_tm, A0_OS, Aii_OS, Vii_tm, Vii_OS;

    
    
    //######### DEFINE 0th and ii component of C^munu ###########
    //tm
    A0_tm = A0_tm_distr;
    Aii_tm =Ak_tm_distr;
    Vii_tm = Vk_tm_distr;
    //OS
    A0_OS = A0_OS_distr;
    Aii_OS = Ak_OS_distr;
    Vii_OS = Vk_OS_distr;
    //###########################################################

    //############ INTERPOLATE CORRELATORS ######################
    vector<boost::math::interpolators::cardinal_cubic_b_spline<double>> Aii_tm_interpol_func;
    vector<boost::math::interpolators::cardinal_cubic_b_spline<double>> Aii_OS_interpol_func;
    for(int t=1;t< Corr.Nt;t++) {
      for(int ijack=0;ijack<Njacks;ijack++) {
	Vfloat Aii_tm_ijack, Aii_OS_ijack;
	for(int t=1;t< Corr.Nt;t++) { Aii_tm_ijack.push_back( Aii_tm.distr_list[t].distr[ijack]*pow(t,3)); Aii_OS_ijack.push_back( Aii_OS.distr_list[t].distr[ijack]*pow(t,3));}

      Aii_tm_interpol_func.emplace_back( Aii_tm_ijack.begin(), Aii_tm_ijack.end(), 1.0, 1.0);
      Aii_OS_interpol_func.emplace_back( Aii_OS_ijack.begin(), Aii_OS_ijack.end(), 1.0, 1.0);
      }
    }

    auto Aii_tm_interpolated =  [&Aii_tm_interpol_func, &a_distr](double t) -> distr_t {
				  distr_t ret;
				  for(int i=0;i<Njacks;i++) ret.distr.push_back( Aii_tm_interpol_func[i](t/a_distr.distr[i]));
				  return ret;
			       };

    auto Aii_OS_interpolated =  [&Aii_OS_interpol_func, &a_distr](double t) -> distr_t {
				  distr_t ret;
				  for(int i=0;i<Njacks;i++) ret.distr.push_back( Aii_OS_interpol_func[i](t/a_distr.distr[i]));
				  return ret;
			       };

    cout<<"Printing Za: "<<endl;
    cout<<"ZA: (from Aii at t=0.6fm):"<< distr_t::f_of_distr(SQRT, Zv*Zv*Aii_tm_interpolated(0.6/0.197327)/Aii_OS_interpolated(0.6/0.197327)).ave()<<" +- "<<distr_t::f_of_distr(SQRT, Zv*Zv*Aii_tm_interpolated(0.6/0.197327)/Aii_OS_interpolated(0.6/0.197327)).err()<<endl;
    cout<<"ZA (from A0P5): "<<L_info.Za_WI_strange<<" +- "<<L_info.Za_WI_strange_err<<endl;
    
 


    

    bool Found_error_less_x_percent=false;
    double x=30;
    //tm
    int tmax_tm_0=1;
    while(!Found_error_less_x_percent && tmax_tm_0 < Corr.Nt/2 -1 ) {
   
      if( (A0_tm.distr_list[tmax_tm_0]).err()/fabs( (A0_tm.distr_list[tmax_tm_0]).ave()) <  0.01*x) tmax_tm_0++;
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

    //OS

    int tmax_OS_0=1;
    while(!Found_error_less_x_percent && tmax_OS_0 < Corr.Nt/2 -1 ) {
   
      if( (A0_OS.distr_list[tmax_OS_0]).err()/fabs( (A0_OS.distr_list[tmax_OS_0]).ave()) <  0.01*x) tmax_OS_0++;
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


    distr_t Edual_tm, Edual_OS, Mrho_tm, Mrho_OS, Rdual_tm, Rdual_OS, grpp_tm, grpp_OS;
    LL.MLLGS_fit_to_corr(Za*Za*Vii_tm, Mpi_distr, a_distr, L_info.L, Edual_tm, Rdual_tm, Mrho_tm, grpp_tm, 5, tmax_tm_1_Vii, Vk_data_tm.Tag[iens]+"_tm");
    LL.MLLGS_fit_to_corr(Zv*Zv*Vii_OS, Mpi_distr, a_distr, L_info.L, Edual_OS, Rdual_OS, Mrho_OS, grpp_OS, 7, tmax_OS_1_Vii, Vk_data_OS.Tag[iens]+"_OS");
    
    
    
       
    Vfloat gppis({grpp_tm.ave(), grpp_OS.ave()});
    Vfloat Mrhos({Mrho_tm.ave(), Mrho_OS.ave()});
    Vfloat Eduals({Edual_tm.ave(),Edual_OS.ave()});
    Vfloat Rduals({Rdual_tm.ave(), Rdual_OS.ave()});
    Vfloat En_tm, Ampl_tm;
    Vfloat En_OS, Ampl_OS;
    LL.Find_pipi_energy_lev( L_info.L  , Mrhos[0],  gppis[0], Mpi, 0.0, En_tm);
    LL.Find_pipi_energy_lev( L_info.L  , Mrhos[1],  gppis[1], Mpi, 0.0, En_OS);
    int N=En_tm.size();
    for(int n=0; n<N;n++) {
      Ampl_tm.push_back( 2.0*resc_GeV.ave()*LL.Amplitude( En_tm[n], L_info.L, Mrhos[0], gppis[0], Mpi, 0.0)); En_tm[n] = 2.0*sqrt( En_tm[n]*En_tm[n] + Mpi*Mpi);
      Ampl_OS.push_back( 2.0*resc_GeV.ave()*LL.Amplitude( En_OS[n], L_info.L, Mrhos[1], gppis[1], Mpi, 0.0)); En_OS[n] = 2.0*sqrt( En_OS[n]*En_OS[n] + Mpi*Mpi);
    }
    VVfloat Ergs({En_tm, En_OS});
    VVfloat Amplitudes({Ampl_tm, Ampl_OS});

    
    
    auto GS_V = [ &a_distr, &resc_GeV](double E, Vfloat &En, Vfloat &Ampl, double Mrho, double Edual, double Rdual) -> double {
		      
		  //build a spectral density with resonances up to 1.5 GeV, from 1.5 GeV use pQCD result. two-pion peaks are smeared over a few MeV interval

		  double result=0.0;
		  double DE= 0.003*a_distr.ave();
		     

		  //pi-pi states
		  for(int n=0; n < (signed)En.size();n++) {
		    if(En[n]< 1.5*a_distr.ave()) {
		      result += Ampl[n]*(1.0/sqrt( 2.0*M_PI*DE*DE))*exp( - ( En[n] - E)*(En[n]-E)/(2.0*DE*DE));
		    }
		  }
		  //pQCD part
		  double res_pQCD=0.0;
		 
		  double sth= Mrho+Edual;
		  if(E> sth) {
		    res_pQCD += resc_GeV.ave()*Rdual*(1.0/(2*M_PI*M_PI))*(0.5*pow(E-sth,2) + 0.5*pow(sth,2)+ sth*(E-sth));
		  }
		  result += res_pQCD;
		  
		  
		  return result;
		    };

    
    auto f_syst_A = [](const function<double(double)> &F) { return 0.0;};
    auto f_syst_A0_tm = [&resc_GeV, &Mpi, &fpi](const function<double(double)> &F) -> double { return resc_GeV.ave()*F(Mpi)*pow(fpi,2)*Mpi/(2.0);};
    auto f_syst_A0_OS = [&resc_GeV, &Mpi_OS, &fpi_OS](const function<double(double)> &F) -> double { return resc_GeV.ave()*F(Mpi_OS)*pow(fpi_OS,2)*Mpi_OS/(2.0);};
      
    auto f_syst_V_tm = [&Ergs, &Amplitudes, &Mrhos, &gppis, &Eduals, &Rduals, &GS_V, &a_distr, &Mpi, &resc_GeV, &LL, &F_free_tm](const function<double(double)> &F) ->double {


		     		     
		      Vfloat systs;
		     		      
		      auto FS_asympt = [ &Ergs, &Amplitudes, &Mrhos, &Eduals, &Rduals, &F, &GS_V, &F_free_tm, &resc_GeV](double E) {
					       double syst = F(E)*GS_V(E, Ergs[0], Amplitudes[0], Mrhos[0], Eduals[0], Rduals[0]);
					       if( E > 1) syst = F(E)*resc_GeV.ave()*F_free_tm(E);
					       return syst;
					     };
		      /*
		      //compute model estimate of vector contribution to R_tau
		      gsl_function_pp<decltype(FS)> SYST(FS);
		      gsl_integration_workspace * w_SYST = gsl_integration_workspace_alloc (1000);
		      gsl_function *G_SYST = static_cast<gsl_function*>(&SYST);
		      double val_mod,err_mod;
		      gsl_integration_qags(G_SYST, E0_l*a_distr.ave(), 4.0,  0.0, 1e-3, 1000, w_SYST, &val_mod, &err_mod);
		      if(err_mod/fabs(val_mod) > 1e-2) crash("Cannot reach accuracy in evaluating systematic");
		      gsl_integration_workspace_free(w_SYST);
		      systs.push_back(fabs(val_mod));
		      */
		      double val_mod, err_mod;
		      gsl_function_pp<decltype(FS_asympt)> SYST_asympt(FS_asympt);
		      gsl_integration_workspace * w_SYST_asympt = gsl_integration_workspace_alloc (1000);
		      gsl_function *G_SYST_asympt = static_cast<gsl_function*>(&SYST_asympt);
		      gsl_integration_qags(G_SYST_asympt, E0_l*a_distr.ave(), 4.0,  0.0, 5e-3, 1000, w_SYST_asympt, &val_mod, &err_mod);
		      //if(err_mod/fabs(val_mod) > 1e-2) crash("Cannot reach accuracy in evaluating systematic");
		      gsl_integration_workspace_free(w_SYST_asympt);
		      systs.push_back( fabs(val_mod));
					   
	      
		      //evaluate GS in infinite volume limit
		      auto FS_infL = [&F, &Mpi, &a_distr, &resc_GeV, &LL, &Mrhos, &gppis, &Eduals,  &Rduals, &F_free_tm](double E) {
				       double sth= Mrhos[0]+Eduals[0];
				       double syst = F(E)*resc_GeV.ave()*(
						    (1.0/(24.0*pow(M_PI,2)))*pow(E,2)*pow(1.0- pow(2.0*Mpi/E,2), 3.0/2.0)*pow(LL.F_pi_GS_mod(E, Mrhos[0], gppis[0],Mpi,0),2)
						    + (E> sth)?(Rduals[0]*(1.0/(2*M_PI*M_PI))*(0.5*pow(E-sth,2) + 0.5*pow(sth,2)+ sth*(E-sth))):0.0  );
				       if( E>1) return F(E)*resc_GeV.ave()*F_free_tm(E);
				       return syst;
				     };
		      gsl_function_pp<decltype(FS_infL)> SYST_infL(FS_infL);
		      gsl_integration_workspace * w_SYST_infL = gsl_integration_workspace_alloc (1000);
		      gsl_function *G_SYST_infL = static_cast<gsl_function*>(&SYST_infL);
		      double val_mod_infL,err_mod_infL;
		      gsl_integration_qags(G_SYST_infL, E0_l*a_distr.ave(), 4.0, 0.0, 5e-3, 1000, w_SYST_infL, &val_mod_infL, &err_mod_infL);
		      gsl_integration_workspace_free(w_SYST_infL);
		      systs.push_back(fabs(val_mod_infL));

			
		      double syst_ave=0.0;
		      double syst_dev=0.0;
		      for (int n=0;n<(signed)systs.size();n++) { syst_ave += systs[n]/systs.size(); syst_dev += systs[n]*systs[n];}

		      syst_dev =  sqrt( (1.0/(systs.size()-1.0))*(syst_dev - systs.size()*syst_ave*syst_ave) );
		      double max = *max_element(systs.begin(), systs.end());
		      return max;
		   
		    };



    auto f_syst_V_OS = [&Ergs, &Amplitudes, &Mrhos, &gppis, &Eduals, &Rduals, &GS_V, &a_distr, &Mpi, &resc_GeV, &LL, &F_free_OS](const function<double(double)> &F) ->double {


		     		     
		      Vfloat systs;
		     		      
		      auto FS_asympt = [ &Ergs, &Amplitudes, &Mrhos, &Eduals, &Rduals, &F, &GS_V, &F_free_OS, &resc_GeV](double E) {
					       double syst = F(E)*GS_V(E, Ergs[1], Amplitudes[1], Mrhos[1], Eduals[1], Rduals[1]);
					       if( E > 1) syst = F(E)*resc_GeV.ave()*F_free_OS(E);
					       return syst;
					     };
		      /*
		      //compute model estimate of vector contribution to R_tau
		      gsl_function_pp<decltype(FS)> SYST(FS);
		      gsl_integration_workspace * w_SYST = gsl_integration_workspace_alloc (1000);
		      gsl_function *G_SYST = static_cast<gsl_function*>(&SYST);
		      double val_mod,err_mod;
		      gsl_integration_qags(G_SYST, E0_l*a_distr.ave(), 4.0,  0.0, 1e-3, 1000, w_SYST, &val_mod, &err_mod);
		      if(err_mod/fabs(val_mod) > 1e-2) crash("Cannot reach accuracy in evaluating systematic");
		      gsl_integration_workspace_free(w_SYST);
		      systs.push_back(fabs(val_mod));
		      */
		      double val_mod, err_mod;
		      gsl_function_pp<decltype(FS_asympt)> SYST_asympt(FS_asympt);
		      gsl_integration_workspace * w_SYST_asympt = gsl_integration_workspace_alloc (1000);
		      gsl_function *G_SYST_asympt = static_cast<gsl_function*>(&SYST_asympt);
		      gsl_integration_qags(G_SYST_asympt, E0_l*a_distr.ave(), 4.0,  0.0, 5e-3, 1000, w_SYST_asympt, &val_mod, &err_mod);
		      //if(err_mod/fabs(val_mod) > 1e-2) crash("Cannot reach accuracy in evaluating systematic");
		      gsl_integration_workspace_free(w_SYST_asympt);
		      systs.push_back( fabs(val_mod));
					   
	      
		      //evaluate GS in infinite volume limit
		      auto FS_infL = [&F, &Mpi, &a_distr, &resc_GeV, &LL, &Mrhos, &gppis, &Eduals, &Rduals, &F_free_OS](double E) {
				       double sth= Mrhos[1]+Eduals[1];
				       double syst = F(E)*resc_GeV.ave()*(
						    (1.0/(24.0*pow(M_PI,2)))*pow(E,2)*pow(1.0- pow(2.0*Mpi/E,2), 3.0/2.0)*pow(LL.F_pi_GS_mod(E, Mrhos[1], gppis[1],Mpi,0),2)
						    + (E> sth)?(Rduals[1]*(1.0/(2*M_PI*M_PI))*(0.5*pow(E-sth,2) + 0.5*pow(sth,2)+ sth*(E-sth))):0.0  );
				       if( E>1) return F(E)*resc_GeV.ave()*F_free_OS(E);
				       return syst;
				     };
		      gsl_function_pp<decltype(FS_infL)> SYST_infL(FS_infL);
		      gsl_integration_workspace * w_SYST_infL = gsl_integration_workspace_alloc (1000);
		      gsl_function *G_SYST_infL = static_cast<gsl_function*>(&SYST_infL);
		      double val_mod_infL,err_mod_infL;
		      gsl_integration_qags(G_SYST_infL, E0_l*a_distr.ave(), 4.0, 0.0, 5e-3, 1000, w_SYST_infL, &val_mod_infL, &err_mod_infL);
		      gsl_integration_workspace_free(w_SYST_infL);
		      systs.push_back(fabs(val_mod_infL));

			
		      double syst_ave=0.0;
		      double syst_dev=0.0;
		      for (int n=0;n<(signed)systs.size();n++) { syst_ave += systs[n]/systs.size(); syst_dev += systs[n]*systs[n];}

		      syst_dev =  sqrt( (1.0/(systs.size()-1.0))*(syst_dev - systs.size()*syst_ave*syst_ave) );
		      double max = *max_element(systs.begin(), systs.end());
		      return max;
		   
		    };



    // lambda function to be used as a smearing func.


      
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



    const auto K0_shifted = [&a_distr](const PrecFloat &E, const PrecFloat &m, const PrecFloat &s, const PrecFloat &E0, int ijack) -> PrecFloat {

      
		      
		      PrecFloat X;
		      PrecFloat X_ave = E/(m_tau*(a_distr.ave() + a_distr.err()));
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

    const auto K1_shifted = [&a_distr](const PrecFloat &E, const PrecFloat &m, const PrecFloat &s, const PrecFloat &E0, int ijack) -> PrecFloat {

		      PrecFloat X;
		      PrecFloat X_ave = E/(m_tau*(a_distr.ave() + a_distr.err()));
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
    

    const auto K0_sharp= [&a_distr](const double &E, const double &E0) {

      if( E > m_tau*a_distr.ave() || E < E0) return 0.0;
      double X= E/(m_tau*a_distr.ave());
      return (1/X)*pow(( 1 -pow(X,2)),2);
    };

    const auto K1_sharp= [&a_distr](const double &E, const double &E0) {

      if( E > m_tau*a_distr.ave() || E < E0) return 0.0;



      double X= E/(m_tau*a_distr.ave());
      return  (1 + 2*pow(X,2))*(1/(X))*pow(( 1 -pow(X,2)),2); 
    };
    


    if(Use_t_up_to_T_half) {
      tmax_tm_1_Aii = tmax_tm_1_Vii= tmax_tm_0 = tmax_OS_0 = tmax_OS_1_Aii = tmax_OS_1_Vii= Corr.Nt/2 -1;
    }

  
  
    
    

    //############# MODEL ESTIMATE ###############
    const auto model_estimate_V_tm = [&K1_sharp, &GS_V,  &a_distr, &Ergs, &Amplitudes, &Mrhos, &Eduals, &Rduals ](double E) -> double {  return K1_sharp(E,E0_l*a_distr.ave())*GS_V(E, Ergs[0], Amplitudes[0], Mrhos[0], Eduals[0], Rduals[0]);};
    const auto model_V_tm =  [&GS_V, &Ergs, &Amplitudes, &a_distr, &Mrhos, &Eduals, &Rduals](double E) -> double {  return GS_V(E, Ergs[0], Amplitudes[0], Mrhos[0], Eduals[0], Rduals[0]);};
    const auto model_estimate_V_OS = [&K1_sharp, &GS_V, &a_distr, &Ergs, &Amplitudes, &Mrhos, &Eduals, &Rduals ](double E) -> double {  return K1_sharp(E,E0_l*a_distr.ave())*GS_V(E, Ergs[1], Amplitudes[1], Mrhos[1], Eduals[1], Rduals[1]);};
    const auto model_V_OS =  [&GS_V, &Ergs, &Amplitudes, &a_distr, &Mrhos, &Eduals, &Rduals](double E) -> double {  return GS_V(E, Ergs[1], Amplitudes[1], Mrhos[1], Eduals[1], Rduals[1]);};
    const auto model_A0_tm = [&Mpi, &fpi, &resc_GeV ](double E) -> double { double s=Mpi*0.001; return resc_GeV.ave()*(fpi*fpi*Mpi/(2.0))*(1.0/sqrt( 2*M_PI*s*s))*exp( -0.5*pow((E-Mpi)/s,2));};
    const auto model_A0_OS = [&Mpi_OS, &fpi_OS, &resc_GeV ](double E) -> double { double s=Mpi_OS*0.001; return resc_GeV.ave()*(fpi_OS*fpi_OS*Mpi_OS/(2.0))*(1.0/sqrt( 2*M_PI*s*s))*exp( -0.5*pow((E-Mpi_OS)/s,2));};
    //compute model estimate of vector contribution to R_tau tm
    gsl_function_pp<decltype(model_estimate_V_tm)> MOD_tm(model_estimate_V_tm);
    gsl_integration_workspace * w_MOD_tm = gsl_integration_workspace_alloc (10000);
    gsl_function *G_MOD_tm = static_cast<gsl_function*>(&MOD_tm);
    double val_mod_tm,err_mod_tm;
    gsl_integration_qagiu(G_MOD_tm, 0.0, 0.0, 1e-4, 10000, w_MOD_tm, &val_mod_tm, &err_mod_tm);
    gsl_integration_workspace_free(w_MOD_tm);
   
    //compute model estimate of vector contribution to R_tau OS
    gsl_function_pp<decltype(model_estimate_V_OS)> MOD_OS(model_estimate_V_OS);
    gsl_integration_workspace * w_MOD_OS = gsl_integration_workspace_alloc (10000);
    gsl_function *G_MOD_OS = static_cast<gsl_function*>(&MOD_OS);
    double val_mod_OS,err_mod_OS;
    gsl_integration_qagiu(G_MOD_OS, 0.0, 0.0, 1e-4, 10000, w_MOD_OS, &val_mod_OS, &err_mod_OS);
    gsl_integration_workspace_free(w_MOD_OS);
    cout<<"############### PRINTING MODEL ESTIMATE FOR V-A contributions ##############"<<endl;
    cout<<"Model estimate R_t(V,tm): "<<val_mod_tm<<" +- "<<err_mod_tm<<endl;
    cout<<"Model estimate R_t(V,OS): "<<val_mod_OS<<" +- "<<err_mod_OS<<endl;
    cout<<"Model estimate R_t(A0, tm): "<<K0_sharp(Mpi,E0_l*a_distr.ave())*resc_GeV.ave()*pow(fpi,2)*Mpi/2.0<<endl;
    cout<<"Model estimate R_t(A0, OS): "<<K0_sharp(Mpi_OS, E0_l*a_distr.ave())*resc_GeV.ave()*pow(fpi_OS,2)*Mpi_OS/2.0<<endl;
    cout<<"############################################################################"<<endl;
    //#############################################

    cout<<"sigma list : {"<<sigma_list[0];
    for(int is=1;is<(signed)sigma_list.size();is++) { cout<<","<<sigma_list[is]<<flush;}
    cout<<"}"<<endl<<flush;
    cout<<"Looping over sigma"<<flush;

    //loop over sigma
    vector<tuple<int,double, double, double>> thread_times_tm(sigma_list.size()), thread_times_OS(sigma_list.size());
    
    #pragma omp parallel for schedule(dynamic)
    for(int is=0; is < (signed)sigma_list.size(); is++) {

      double s= sigma_list[is];
      
      distr_t Br_sigma_A0_tm;
      distr_t Br_sigma_Aii_tm;
      distr_t Br_sigma_Vii_tm;
      distr_t Br_sigma_A0_OS;
      distr_t Br_sigma_Aii_OS;
      distr_t Br_sigma_Vii_OS;

      distr_t Br_s_sigma_A0_tm;
      distr_t Br_s_sigma_Aii_tm;
      distr_t Br_s_sigma_Vii_tm;
      distr_t Br_s_sigma_A0_OS;
      distr_t Br_s_sigma_Aii_OS;
      distr_t Br_s_sigma_Vii_OS;
   
      //int tmax= T/2 -4;
      double lA0_tm, lAii_tm, lVii_tm;
      double lA0_OS, lAii_OS, lVii_OS;
      double syst_A0_tm, syst_Aii_tm, syst_Vii_tm;
      double syst_A0_OS, syst_Aii_OS, syst_Vii_OS;

      double syst_s_A0_tm, syst_s_Aii_tm, syst_s_Vii_tm;
      double syst_s_A0_OS, syst_s_Aii_OS, syst_s_Vii_OS;
      
      double mult=1e4;
      if( (beta > 2) && (s < 0.15) ) mult=1e5;

    
     
  
      auto start = chrono::system_clock::now();
      if( (beta > 2) && Use_Customized_plateaux) mult=  Customized_plateaux_tau_spectre( beta, Emax,  "Aii", "tm" , s, Vk_data_tm.Tag[iens] );
      Br_sigma_Aii_tm = Get_Laplace_transfo(  0.0,  s, E0_l*a_distr.ave(),  T, tmax_tm_1_Aii, prec, SM_TYPE_1,K1, Aii_tm, syst_Aii_tm, mult, lAii_tm, MODE, "tm", "Aii_light_"+Vk_data_tm.Tag[iens], -1,0, resc_GeV*Zv*Zv, 0.0, "tau_decay", cov_Ak_tm, fake_func,0, fake_func_d ,  Is_Emax_Finite, Emax, beta);
      Br_s_sigma_Aii_tm = Get_Laplace_transfo(  0.0,  s, E0_l*a_distr.ave(),  T, tmax_tm_1_Aii, prec, SM_TYPE_1,K1_shifted, Aii_tm, syst_s_Aii_tm, mult, lAii_tm, MODE, "tm", "Aii_s_light_"+Vk_data_tm.Tag[iens], -1,0, resc_GeV*Zv*Zv, 0.0, "tau_decay", cov_Ak_tm, fake_func,0, fake_func_d ,  Is_Emax_Finite, Emax, beta);
      syst_s_Aii_tm = fabs( Br_s_sigma_Aii_tm.ave() - Br_sigma_Aii_tm.ave());
      auto end = chrono::system_clock::now();
      cout<<"node: "<<_hostname<<", rank: "<<rank<<", thread_id: "<<omp_get_thread_num()<<" core-id: "<<sched_getcpu()<<endl<<flush;
      chrono::duration<double> elapsed_seconds = end-start;
      double time_Aii_tm= elapsed_seconds.count();
      if(tau_verbosity_lev) {
	cout<<endl<<flush;
	cout<<"Elapsed time[Aii_tm, sigma: "<<s<<", Ens: "<<Vk_data_tm.Tag[iens]<<", rank: "<<rank<<", #thread="<<omp_get_thread_num()<<"] : "<<time_Aii_tm<<" s"<<endl<<flush;
      }
      else cout<<"."<<flush;

     
      start = chrono::system_clock::now();
      if( (beta > 2) && Use_Customized_plateaux) mult=  Customized_plateaux_tau_spectre( beta, Emax, "Aii", "OS" , s, Vk_data_tm.Tag[iens] );
      Br_sigma_Aii_OS = Get_Laplace_transfo(  0.0,  s, E0_l*a_distr.ave(),  T, tmax_OS_1_Aii, prec, SM_TYPE_1,K1, Aii_OS, syst_Aii_OS, mult, lAii_OS, MODE, "OS", "Aii_light_"+Vk_data_tm.Tag[iens], -1,0,resc_GeV*Za*Za, 0.0, "tau_decay", cov_Ak_OS, fake_func,0, fake_func_d ,  Is_Emax_Finite, Emax, beta);
      Br_s_sigma_Aii_OS = Get_Laplace_transfo(  0.0,  s, E0_l*a_distr.ave(),  T, tmax_OS_1_Aii, prec, SM_TYPE_1,K1_shifted, Aii_OS, syst_s_Aii_OS, mult, lAii_OS, MODE, "OS", "Aii_s_light_"+Vk_data_tm.Tag[iens], -1,0,resc_GeV*Za*Za, 0.0, "tau_decay", cov_Ak_OS, fake_func,0, fake_func_d ,  Is_Emax_Finite, Emax, beta);
      syst_s_Aii_OS = fabs( Br_s_sigma_Aii_OS.ave() - Br_sigma_Aii_OS.ave());
     
      end = chrono::system_clock::now();
      elapsed_seconds = end-start;
      double time_Aii_OS= elapsed_seconds.count();
      if(tau_verbosity_lev) {
	cout<<endl<<flush;
	cout<<"Elapsed time[Aii_OS, sigma: "<<s<<", Ens: "<<Vk_data_tm.Tag[iens]<<", rank: "<<rank<<", #thread="<<omp_get_thread_num()<<"] : "<<time_Aii_OS<<" s"<<endl<<flush;
      }
      else cout<<"."<<flush;
      

     
      start = chrono::system_clock::now();
      if( (beta > 2) && Use_Customized_plateaux) mult=  Customized_plateaux_tau_spectre( beta, Emax, "Vii", "tm" , s, Vk_data_tm.Tag[iens] );
      Br_sigma_Vii_tm = Get_Laplace_transfo(  0.0,  s, E0_l*a_distr.ave(),  T, tmax_tm_1_Vii, prec, SM_TYPE_1,K1, Vii_tm, syst_Vii_tm, mult, lVii_tm, MODE, "tm", "Vii_light_"+Vk_data_tm.Tag[iens], -1,0, resc_GeV*Za*Za, 0.0, "tau_decay", cov_Vk_tm, f_syst_V_tm,1, model_V_tm, Is_Emax_Finite, Emax, beta);
      Br_s_sigma_Vii_tm = Get_Laplace_transfo(  0.0,  s, E0_l*a_distr.ave(),  T, tmax_tm_1_Vii, prec, SM_TYPE_1,K1_shifted, Vii_tm, syst_s_Vii_tm, mult, lVii_tm, MODE, "tm", "Vii_s_light_"+Vk_data_tm.Tag[iens], -1,0, resc_GeV*Za*Za, 0.0, "tau_decay", cov_Vk_tm, f_syst_V_tm,1, model_V_tm, Is_Emax_Finite, Emax, beta);
      syst_s_Vii_tm = fabs( Br_s_sigma_Vii_tm.ave() - Br_sigma_Vii_tm.ave());
      end = chrono::system_clock::now();
      elapsed_seconds = end-start;
      double time_Vii_tm= elapsed_seconds.count();
      if(tau_verbosity_lev) {
	cout<<endl<<flush;
	cout<<"Elapsed time[Vii_tm, sigma: "<<s<<", Ens: "<<Vk_data_tm.Tag[iens]<<", rank: "<<rank<<", #thread="<<omp_get_thread_num()<<"] : "<<time_Vii_tm<<" s"<<endl<<flush;
      }
      else cout<<"."<<flush;

     
      start = chrono::system_clock::now();
      if( (beta > 2) && Use_Customized_plateaux) mult=  Customized_plateaux_tau_spectre( beta, Emax,  "Vii", "OS" , s, Vk_data_tm.Tag[iens] );
      Br_sigma_Vii_OS = Get_Laplace_transfo(  0.0,  s, E0_l*a_distr.ave(),  T, tmax_OS_1_Vii, prec, SM_TYPE_1,K1, Vii_OS, syst_Vii_OS, mult, lVii_OS, MODE, "OS", "Vii_light_"+Vk_data_tm.Tag[iens],-1,0, resc_GeV*Zv*Zv, 0.0, "tau_decay", cov_Vk_OS, f_syst_V_OS,1, model_V_OS ,  Is_Emax_Finite, Emax, beta);
      Br_s_sigma_Vii_OS = Get_Laplace_transfo(  0.0,  s, E0_l*a_distr.ave(),  T, tmax_OS_1_Vii, prec, SM_TYPE_1,K1_shifted, Vii_OS, syst_s_Vii_OS, mult, lVii_OS, MODE, "OS", "Vii_s_light_"+Vk_data_tm.Tag[iens],-1,0, resc_GeV*Zv*Zv, 0.0, "tau_decay", cov_Vk_OS, f_syst_V_OS,1, model_V_OS ,  Is_Emax_Finite, Emax, beta);
      syst_s_Vii_OS = fabs( Br_s_sigma_Vii_OS.ave() - Br_sigma_Vii_OS.ave());
      end = chrono::system_clock::now();
      elapsed_seconds = end-start;
      double time_Vii_OS= elapsed_seconds.count();
      if(tau_verbosity_lev) {
	cout<<endl<<flush;
	cout<<"Elapsed time[Vii_OS, sigma: "<<s<<", Ens: "<<Vk_data_tm.Tag[iens]<<", rank: "<<rank<<", #thread="<<omp_get_thread_num()<<"] : "<<time_Vii_OS<<" s"<<endl<<flush;
      }
      else cout<<"."<<flush;



      start = chrono::system_clock::now();
      if( (beta > 2) && Use_Customized_plateaux) mult=  Customized_plateaux_tau_spectre( beta, Emax, "A0", "tm" , s, Vk_data_tm.Tag[iens] );
      Br_sigma_A0_tm = Get_Laplace_transfo(  0.0,  s, E0_l*a_distr.ave(),  T, tmax_tm_0, prec, SM_TYPE_0,K0, -1*A0_tm, syst_A0_tm, mult, lA0_tm, MODE, "tm", "A0_light_"+Vk_data_tm.Tag[iens], -1, 0, resc_GeV*Zv*Zv, 0.0, "tau_decay", cov_A0_tm, f_syst_A0_tm,1, model_A0_tm,  Is_Emax_Finite, Emax, beta );
      Br_s_sigma_A0_tm = Get_Laplace_transfo(  0.0,  s, E0_l*a_distr.ave(),  T, tmax_tm_0, prec, SM_TYPE_0,K0_shifted, -1*A0_tm, syst_s_A0_tm, mult, lA0_tm, MODE, "tm", "A0_s_light_"+Vk_data_tm.Tag[iens], -1, 0, resc_GeV*Zv*Zv, 0.0, "tau_decay", cov_A0_tm, f_syst_A0_tm,1, model_A0_tm,  Is_Emax_Finite, Emax, beta );
      syst_s_A0_tm = fabs( Br_s_sigma_A0_tm.ave() - Br_sigma_A0_tm.ave());
      end = chrono::system_clock::now();
      elapsed_seconds = end-start;
      double time_A0_tm= elapsed_seconds.count();
      if(tau_verbosity_lev) {
	cout<<endl<<flush;
	cout<<"Elapsed time[A0_tm, sigma: "<<s<<", Ens: "<<Vk_data_tm.Tag[iens]<<", rank: "<<rank<<", #thread="<<omp_get_thread_num()<<"] : "<<time_A0_tm<<" s"<<endl<<flush;
      }
      else cout<<"."<<flush;
      

      start = chrono::system_clock::now();
      if( (beta > 2) && Use_Customized_plateaux) mult=  Customized_plateaux_tau_spectre( beta, Emax, "A0", "OS" , s, Vk_data_tm.Tag[iens] );
      Br_sigma_A0_OS = Get_Laplace_transfo(  0.0,  s, E0_l*a_distr.ave(),  T, tmax_OS_0, prec, SM_TYPE_0,K0, -1*A0_OS, syst_A0_OS, mult, lA0_OS, MODE, "OS", "A0_light_"+Vk_data_tm.Tag[iens], -1, 0, resc_GeV*Za*Za, 0.0, "tau_decay", cov_A0_OS, f_syst_A0_OS,1, model_A0_OS,  Is_Emax_Finite, Emax, beta );
      Br_s_sigma_A0_OS = Get_Laplace_transfo(  0.0,  s, E0_l*a_distr.ave(),  T, tmax_OS_0, prec, SM_TYPE_0,K0_shifted, -1*A0_OS, syst_s_A0_OS, mult, lA0_OS, MODE, "OS", "A0_s_light_"+Vk_data_tm.Tag[iens], -1, 0, resc_GeV*Za*Za, 0.0, "tau_decay", cov_A0_OS, f_syst_A0_OS,1, model_A0_OS,  Is_Emax_Finite, Emax, beta );
      syst_s_A0_OS = fabs( Br_s_sigma_A0_OS.ave() - Br_sigma_A0_OS.ave());
      end = chrono::system_clock::now();
      elapsed_seconds = end-start;
      double time_A0_OS= elapsed_seconds.count();
      if(tau_verbosity_lev) {
	cout<<endl<<flush;
	cout<<"Elapsed time[A0_OS, sigma: "<<s<<", Ens: "<<Vk_data_tm.Tag[iens]<<", rank: "<<rank<<", #thread="<<omp_get_thread_num()<<"] : "<<time_A0_OS<<" s"<<endl<<flush;
      }
      else cout<<"."<<flush;

      thread_times_tm[is] = make_tuple(omp_get_thread_num(), time_A0_tm, time_Aii_tm, time_Vii_tm);
      thread_times_OS[is] = make_tuple(omp_get_thread_num(), time_A0_OS, time_Aii_OS, time_Vii_OS);
      
           
      syst_per_ens_tm_A0[iens][is] = sqrt( pow(syst_A0_tm,2)+ pow(syst_s_A0_tm,2)); // syst_A0_tm;
      syst_per_ens_tm_Ak[iens][is]=  sqrt( pow(syst_Aii_tm,2) + pow(syst_s_Aii_tm,2)); //syst_Aii_tm;
      syst_per_ens_tm_Vk[iens][is]=  sqrt( pow(syst_Vii_tm,2) + pow(syst_s_Vii_tm,2));//syst_Vii_tm;
      syst_per_ens_tm[iens][is]= sqrt( pow(syst_A0_tm,2)+ pow(syst_Aii_tm,2)+ pow(syst_Vii_tm,2));
      syst_per_ens_OS_A0[iens][is]= sqrt( pow(syst_A0_OS,2) + pow(syst_s_A0_OS,2));  //syst_A0_OS;
      syst_per_ens_OS_Ak[iens][is]= sqrt( pow(syst_Aii_OS,2) + pow(syst_s_Aii_OS,2));  //syst_Aii_OS;
      syst_per_ens_OS_Vk[iens][is]= sqrt( pow(syst_Vii_OS,2) + pow(syst_s_Vii_OS,2));//syst_Vii_OS;
      syst_per_ens_OS[iens][is]= sqrt( pow(syst_A0_OS,2)+ pow(syst_Aii_OS,2)+ pow(syst_Vii_OS,2));

     
      
      distr_t Br_sigma_tm = Br_sigma_Aii_tm + Br_sigma_Vii_tm + Br_sigma_A0_tm;
      distr_t Br_sigma_OS = Br_sigma_Aii_OS + Br_sigma_Vii_OS + Br_sigma_A0_OS;
      Br_tau_tm[iens].distr_list[is] = Br_sigma_tm;
      Br_Aii_tau_tm[iens].distr_list[is] = Br_sigma_Aii_tm;
      Br_Vii_tau_tm[iens].distr_list[is] = Br_sigma_Vii_tm;
      Br_A0_tau_tm[iens].distr_list[is] = Br_sigma_A0_tm;
      Br_tau_OS[iens].distr_list[is] = Br_sigma_OS;
      Br_Aii_tau_OS[iens].distr_list[is] = Br_sigma_Aii_OS;
      Br_Vii_tau_OS[iens].distr_list[is] = Br_sigma_Vii_OS;
      Br_A0_tau_OS[iens].distr_list[is] = Br_sigma_A0_OS;


      cout<<"Ensemble: "<<Vk_data_tm.Tag[iens]<<", sigma: "<<s<<" completed!"<<endl<<flush;
     
   
     
    }
    
    cout<<endl;
    cout<<"Finished ensemble: "<<Vk_data_tm.Tag[iens]<<"########################"<<endl<<flush;
    cout<<"Summary of performances: "<<endl<<flush;
    cout<<"sigma #thread  A0    Aii     Vii"<<endl<<flush;
    cout<<"- - - - - - - - - - - - - - - - - - - - - "<<endl<<flush;
    for(int is=0; is < (signed)sigma_list.size();is++) {
      cout<<sigma_list[is]<<", "<<get<0>(thread_times_tm[is])<<": "<<get<1>(thread_times_tm[is])<<" s, "<<get<2>(thread_times_tm[is])<<" s, "<<get<3>(thread_times_tm[is])<<" s"<<endl<<flush;
      cout<<sigma_list[is]<<", "<<get<0>(thread_times_OS[is])<<": "<<get<1>(thread_times_OS[is])<<" s, "<<get<2>(thread_times_OS[is])<<" s, "<<get<3>(thread_times_OS[is])<<" s"<<endl<<flush;
      cout<<"- - - - - - - - - - - - - - - - - - - - "<<endl<<flush;
    }
    
    
   
  }
  
    
   
  



  //Print to File
  for(int iens=0; iens<Nens;iens++) {

    
    
    Print_To_File({}, {sigma_list, Br_tau_tm[iens].ave(), Br_tau_tm[iens].err(), syst_per_ens_tm[iens] , Br_tau_OS[iens].ave(), Br_tau_OS[iens].err(), syst_per_ens_OS[iens]}, "../data/tau_decay/"+Tag_reco_type+"/light/Br/br_"+MODE+"_sm_func_mode_"+to_string(sm_func_mode)+"_"+Vk_data_tm.Tag[iens]+".dat", "", "#sigma Br[tm] Br[OS]");
    Print_To_File({}, {sigma_list, Br_A0_tau_tm[iens].ave(), Br_A0_tau_tm[iens].err(), syst_per_ens_tm_A0[iens], Br_Aii_tau_tm[iens].ave(), Br_Aii_tau_tm[iens].err(), syst_per_ens_tm_Ak[iens], Br_Vii_tau_tm[iens].ave(), Br_Vii_tau_tm[iens].err(), syst_per_ens_tm_Vk[iens]}, "../data/tau_decay/"+Tag_reco_type+"/light/Br/br_contrib_tm_"+MODE+"_sm_func_mode_"+to_string(sm_func_mode)+"_"+Vk_data_tm.Tag[iens]+".dat", "", "#sigma Br_A0 Br_Aii Br_Vii");
    Print_To_File({}, {sigma_list, Br_A0_tau_OS[iens].ave(), Br_A0_tau_OS[iens].err(), syst_per_ens_OS_A0[iens], Br_Aii_tau_OS[iens].ave(), Br_Aii_tau_OS[iens].err(), syst_per_ens_OS_Ak[iens], Br_Vii_tau_OS[iens].ave(), Br_Vii_tau_OS[iens].err(), syst_per_ens_OS_Vk[iens]}, "../data/tau_decay/"+Tag_reco_type+"/light/Br/br_contrib_OS_"+MODE+"_sm_func_mode_"+to_string(sm_func_mode)+"_"+Vk_data_tm.Tag[iens]+".dat", "", "#sigma Br_A0 Br_Aii Br_Vii");

    
    //RESCALED SYSTEMATIC
    Vfloat syst_per_ens_tm_fpi_resc, syst_per_ens_OS_fpi_resc;
    Vfloat syst_per_ens_tm_A0_fpi_resc, syst_per_ens_tm_Ak_fpi_resc, syst_per_ens_tm_Vk_fpi_resc;
    Vfloat syst_per_ens_OS_A0_fpi_resc, syst_per_ens_OS_Ak_fpi_resc, syst_per_ens_OS_Vk_fpi_resc;

    for(int is=0;is<(signed)sigma_list.size();is++) {

      syst_per_ens_tm_fpi_resc.push_back( syst_per_ens_tm[iens][is]*RESC_FP.ave(iens));
      syst_per_ens_OS_fpi_resc.push_back( syst_per_ens_OS[iens][is]*RESC_FP.ave(iens));

      syst_per_ens_tm_A0_fpi_resc.push_back( syst_per_ens_tm_A0[iens][is]*RESC_FP.ave(iens));
      syst_per_ens_tm_Ak_fpi_resc.push_back( syst_per_ens_tm_Ak[iens][is]*RESC_FP.ave(iens));
      syst_per_ens_tm_Vk_fpi_resc.push_back( syst_per_ens_tm_Vk[iens][is]*RESC_FP.ave(iens));

      syst_per_ens_OS_A0_fpi_resc.push_back( syst_per_ens_OS_A0[iens][is]*RESC_FP.ave(iens));
      syst_per_ens_OS_Ak_fpi_resc.push_back( syst_per_ens_OS_Ak[iens][is]*RESC_FP.ave(iens));
      syst_per_ens_OS_Vk_fpi_resc.push_back( syst_per_ens_OS_Vk[iens][is]*RESC_FP.ave(iens));

    }
    

    Print_To_File({}, {sigma_list, (RESC_FP.distr_list[iens]*Br_tau_tm[iens]).ave(), (RESC_FP.distr_list[iens]*Br_tau_tm[iens]).err(), syst_per_ens_tm_fpi_resc , (RESC_FP.distr_list[iens]*Br_tau_OS[iens]).ave(), (RESC_FP.distr_list[iens]*Br_tau_OS[iens]).err(), syst_per_ens_OS_fpi_resc}, "../data/tau_decay/"+Tag_reco_type+"/light/Br/br_resc_fp_"+MODE+"_sm_func_mode_"+to_string(sm_func_mode)+"_"+Vk_data_tm.Tag[iens]+".dat", "", "#sigma Br[tm] Br[OS]");
    Print_To_File({}, {sigma_list, (RESC_FP.distr_list[iens]*Br_A0_tau_tm[iens]).ave(), (RESC_FP.distr_list[iens]*Br_A0_tau_tm[iens]).err(), syst_per_ens_tm_A0_fpi_resc, (RESC_FP.distr_list[iens]*Br_Aii_tau_tm[iens]).ave(), (RESC_FP.distr_list[iens]*Br_Aii_tau_tm[iens]).err(), syst_per_ens_tm_Ak_fpi_resc, (RESC_FP.distr_list[iens]*Br_Vii_tau_tm[iens]).ave(), (RESC_FP.distr_list[iens]*Br_Vii_tau_tm[iens]).err(), syst_per_ens_tm_Vk_fpi_resc}, "../data/tau_decay/"+Tag_reco_type+"/light/Br/br_resc_fp_contrib_tm_"+MODE+"_sm_func_mode_"+to_string(sm_func_mode)+"_"+Vk_data_tm.Tag[iens]+".dat", "", "#sigma Br_A0 Br_Aii Br_Vii");
    Print_To_File({}, {sigma_list, (RESC_FP.distr_list[iens]*Br_A0_tau_OS[iens]).ave(), (RESC_FP.distr_list[iens]*Br_A0_tau_OS[iens]).err(), syst_per_ens_OS_A0_fpi_resc, (RESC_FP.distr_list[iens]*Br_Aii_tau_OS[iens]).ave(), (RESC_FP.distr_list[iens]*Br_Aii_tau_OS[iens]).err(), syst_per_ens_OS_Ak_fpi_resc, (RESC_FP.distr_list[iens]*Br_Vii_tau_OS[iens]).ave(), (RESC_FP.distr_list[iens]*Br_Vii_tau_OS[iens]).err(), syst_per_ens_OS_Vk_fpi_resc}, "../data/tau_decay/"+Tag_reco_type+"/light/Br/br_resc_fp_contrib_OS_"+MODE+"_sm_func_mode_"+to_string(sm_func_mode)+"_"+Vk_data_tm.Tag[iens]+".dat", "", "#sigma Br_A0 Br_Aii Br_Vii");
    
  }

  cout<<"Output per ensemble printed!"<<endl;

  //Print all ens for each sigma
  for(int is=0; is<(signed)sigma_list.size();is++) {

    Vfloat syst_tm, syst_OS;
    Vfloat syst_tm_fpi_resc,  syst_OS_fpi_resc;
    distr_t_list Br_tau_tm_fixed_sigma(UseJack), Br_tau_OS_fixed_sigma(UseJack);
    for(int iens=0;iens<Nens;iens++) {
      Br_tau_tm_fixed_sigma.distr_list.push_back( Br_tau_tm[iens].distr_list[is]);
      Br_tau_OS_fixed_sigma.distr_list.push_back( Br_tau_OS[iens].distr_list[is]);
      syst_tm.push_back( syst_per_ens_tm[iens][is]);
      syst_OS.push_back( syst_per_ens_OS[iens][is]);
      syst_tm_fpi_resc.push_back( syst_per_ens_tm[iens][is]*RESC_FP.ave(iens));
      syst_OS_fpi_resc.push_back( syst_per_ens_OS[iens][is]*RESC_FP.ave(iens));
    }

    Print_To_File(Vk_data_tm.Tag,{ Br_tau_tm_fixed_sigma.ave(), Br_tau_tm_fixed_sigma.err(), syst_tm, Br_tau_OS_fixed_sigma.ave(), Br_tau_OS_fixed_sigma.err(), syst_OS},"../data/tau_decay/"+Tag_reco_type+"/light/Br/br_sigma_"+to_string_with_precision(sigma_list[is],3)+"_sm_func_mode_"+to_string(sm_func_mode)+".dat", "", "#Ens tm  OS");

    Print_To_File(Vk_data_tm.Tag,{ (RESC_FP*Br_tau_tm_fixed_sigma).ave(), (RESC_FP*Br_tau_tm_fixed_sigma).err(), syst_tm_fpi_resc, (RESC_FP*Br_tau_OS_fixed_sigma).ave(), (RESC_FP*Br_tau_OS_fixed_sigma).err(), syst_OS_fpi_resc},"../data/tau_decay/"+Tag_reco_type+"/light/Br/br_resc_fp_sigma_"+to_string_with_precision(sigma_list[is],3)+"_sm_func_mode_"+to_string(sm_func_mode)+".dat", "", "#Ens tm  OS");
    
  }

  cout<<"output per sigma printed!"<<endl;



  //STORE JACKKNIFE DISTRIBUTIONS
  //light
  for(int i_ens=0; i_ens<Nens;i_ens++) {
    boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/light/jackknife/tm/A0A0/"+Vk_data_tm.Tag[i_ens]);
    boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/light/jackknife/tm/AkAk/"+Vk_data_tm.Tag[i_ens]);
    boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/light/jackknife/tm/VkVk/"+Vk_data_tm.Tag[i_ens]);
    boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/light/jackknife/OS/A0A0/"+Vk_data_tm.Tag[i_ens]);
    boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/light/jackknife/OS/AkAk/"+Vk_data_tm.Tag[i_ens]);
    boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/light/jackknife/OS/VkVk/"+Vk_data_tm.Tag[i_ens]);

    boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/light/jackknife/tm/A0A0/"+Vk_data_tm.Tag[i_ens]+"/sm_func_mode_"+to_string(sm_func_mode));
    boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/light/jackknife/tm/AkAk/"+Vk_data_tm.Tag[i_ens]+"/sm_func_mode_"+to_string(sm_func_mode));
    boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/light/jackknife/tm/VkVk/"+Vk_data_tm.Tag[i_ens]+"/sm_func_mode_"+to_string(sm_func_mode));
    boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/light/jackknife/OS/A0A0/"+Vk_data_tm.Tag[i_ens]+"/sm_func_mode_"+to_string(sm_func_mode));
    boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/light/jackknife/OS/AkAk/"+Vk_data_tm.Tag[i_ens]+"/sm_func_mode_"+to_string(sm_func_mode));
    boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/light/jackknife/OS/VkVk/"+Vk_data_tm.Tag[i_ens]+"/sm_func_mode_"+to_string(sm_func_mode));

  
    for(int is=0; is<(signed)sigma_list.size();is++) {
      //print jackknife distribution for tm and OS
      ofstream Print_tm_A0A0("../data/tau_decay/"+Tag_reco_type+"/light/jackknife/tm/A0A0/"+Vk_data_tm.Tag[i_ens]+"/sm_func_mode_"+to_string(sm_func_mode)+"/sigma_"+to_string_with_precision(sigma_list[is],3)+".jack");
      ofstream Print_tm_AkAk("../data/tau_decay/"+Tag_reco_type+"/light/jackknife/tm/AkAk/"+Vk_data_tm.Tag[i_ens]+"/sm_func_mode_"+to_string(sm_func_mode)+"/sigma_"+to_string_with_precision(sigma_list[is],3)+".jack");
      ofstream Print_tm_VkVk("../data/tau_decay/"+Tag_reco_type+"/light/jackknife/tm/VkVk/"+Vk_data_tm.Tag[i_ens]+"/sm_func_mode_"+to_string(sm_func_mode)+"/sigma_"+to_string_with_precision(sigma_list[is],3)+".jack");

      ofstream Print_OS_A0A0("../data/tau_decay/"+Tag_reco_type+"/light/jackknife/OS/A0A0/"+Vk_data_tm.Tag[i_ens]+"/sm_func_mode_"+to_string(sm_func_mode)+"/sigma_"+to_string_with_precision(sigma_list[is],3)+".jack");
      ofstream Print_OS_AkAk("../data/tau_decay/"+Tag_reco_type+"/light/jackknife/OS/AkAk/"+Vk_data_tm.Tag[i_ens]+"/sm_func_mode_"+to_string(sm_func_mode)+"/sigma_"+to_string_with_precision(sigma_list[is],3)+".jack");
      ofstream Print_OS_VkVk("../data/tau_decay/"+Tag_reco_type+"/light/jackknife/OS/VkVk/"+Vk_data_tm.Tag[i_ens]+"/sm_func_mode_"+to_string(sm_func_mode)+"/sigma_"+to_string_with_precision(sigma_list[is],3)+".jack");
     
      for(int ijack=0; ijack<Njacks;ijack++) {
	Print_tm_A0A0<<Br_A0_tau_tm[i_ens].distr_list[is].distr[ijack]<<endl;
	Print_tm_AkAk<<Br_Aii_tau_tm[i_ens].distr_list[is].distr[ijack]<<endl;
	Print_tm_VkVk<<Br_Vii_tau_tm[i_ens].distr_list[is].distr[ijack]<<endl;
	Print_OS_A0A0<<Br_A0_tau_OS[i_ens].distr_list[is].distr[ijack]<<endl;
	Print_OS_AkAk<<Br_Aii_tau_OS[i_ens].distr_list[is].distr[ijack]<<endl;
	Print_OS_VkVk<<Br_Vii_tau_OS[i_ens].distr_list[is].distr[ijack]<<endl;
      }
      Print_tm_A0A0.close();
      Print_tm_AkAk.close();
      Print_tm_VkVk.close();
      Print_OS_A0A0.close();
      Print_OS_AkAk.close();
      Print_OS_VkVk.close();
     
    }

   
  }
  cout<<"Jackknives printed!"<<endl;
  }


  



  if(Perform_continuum_extrapolation) {

    cout<<"Performing continuum limit extrapolation"<<endl;


    
    vector<string> Contribs({"A0A0", "AkAk", "VkVk", "tot", "VA", "VMA", "AX", "T"});
    //vector<string> Fit_types({"tm", "OS", "comb"});
    vector<string> Fit_types({"comb"});
    //vector<string> poly_types({"const", "linear"});
    vector<string> poly_types({"const", "linear", "tm_linear", "OS_linear"});
    
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

    vector<distr_t_list> A0A0_tm_all_s(sigma_list.size());
    vector<distr_t_list> AkAk_tm_all_s(sigma_list.size());
    vector<distr_t_list> VkVk_tm_all_s(sigma_list.size());

    vector<distr_t_list> A0A0_OS_all_s(sigma_list.size());
    vector<distr_t_list> AkAk_OS_all_s(sigma_list.size());
    vector<distr_t_list> VkVk_OS_all_s(sigma_list.size());

    //load systematic errors

     
    VVfloat syst_A0A0_tm(Nens), syst_AkAk_tm(Nens), syst_VkVk_tm(Nens);
    VVfloat syst_A0A0_OS(Nens), syst_AkAk_OS(Nens), syst_VkVk_OS(Nens);
    for(int iens=0;iens<Nens;iens++) {
      //A0A0
      syst_A0A0_tm[iens] = Read_From_File("../data/tau_decay/"+Tag_reco_type+"/light/Br/br_contrib_tm_"+MODE+"_sm_func_mode_"+to_string(sm_func_mode)+"_"+Vk_data_tm.Tag[iens]+".dat", 4,  11,1);
      syst_A0A0_OS[iens] = Read_From_File("../data/tau_decay/"+Tag_reco_type+"/light/Br/br_contrib_OS_"+MODE+"_sm_func_mode_"+to_string(sm_func_mode)+"_"+Vk_data_tm.Tag[iens]+".dat", 4,  11,1);
      //AkAk
      syst_AkAk_tm[iens] = Read_From_File("../data/tau_decay/"+Tag_reco_type+"/light/Br/br_contrib_tm_"+MODE+"_sm_func_mode_"+to_string(sm_func_mode)+"_"+Vk_data_tm.Tag[iens]+".dat", 7,  11,1);
      syst_AkAk_OS[iens] = Read_From_File("../data/tau_decay/"+Tag_reco_type+"/light/Br/br_contrib_OS_"+MODE+"_sm_func_mode_"+to_string(sm_func_mode)+"_"+Vk_data_tm.Tag[iens]+".dat", 7,  11,1);
      //VkVk
      syst_VkVk_tm[iens] = Read_From_File("../data/tau_decay/"+Tag_reco_type+"/light/Br/br_contrib_tm_"+MODE+"_sm_func_mode_"+to_string(sm_func_mode)+"_"+Vk_data_tm.Tag[iens]+".dat", 10,  11,1);
      syst_VkVk_OS[iens] = Read_From_File("../data/tau_decay/"+Tag_reco_type+"/light/Br/br_contrib_OS_"+MODE+"_sm_func_mode_"+to_string(sm_func_mode)+"_"+Vk_data_tm.Tag[iens]+".dat", 10,  11,1);


      cout<<"Ens: "<<Vk_data_tm.Tag[iens]<<" "<<Tag_reco_type<<"  syst tm A0, Ak Vk: "<<syst_A0A0_tm[iens][0]<<" "<<syst_AkAk_tm[iens][0]<<" "<<syst_VkVk_tm[iens][0]<<endl;
    }


     

    for(int is=0; is < (signed)sigma_list.size() ; is++) { //loop over sigma

      //load data
      distr_t_list  A0A0_tm(UseJack), AkAk_tm(UseJack), VkVk_tm(UseJack);
      distr_t_list  A0A0_OS(UseJack), AkAk_OS(UseJack), VkVk_OS(UseJack);

      GaussianMersenne GM_sigma(654324);

      for(int iens=0;iens<Nens;iens++) {


	//load A0A0 jack distribution
	A0A0_tm.distr_list.emplace_back(UseJack, Read_From_File("../data/tau_decay/"+Tag_reco_type+"/light/jackknife/tm/A0A0/"+Vk_data_tm.Tag[iens]+"/sm_func_mode_"+to_string(sm_func_mode)+"/sigma_"+to_string_with_precision(sigma_list[is],3)+".jack", 0,1));
	A0A0_OS.distr_list.emplace_back(UseJack, Read_From_File("../data/tau_decay/"+Tag_reco_type+"/light/jackknife/OS/A0A0/"+Vk_data_tm.Tag[iens]+"/sm_func_mode_"+to_string(sm_func_mode)+"/sigma_"+to_string_with_precision(sigma_list[is],3)+".jack", 0,1));

	
	//load AkAk jack distribution
	AkAk_tm.distr_list.emplace_back(UseJack, Read_From_File("../data/tau_decay/"+Tag_reco_type+"/light/jackknife/tm/AkAk/"+Vk_data_tm.Tag[iens]+"/sm_func_mode_"+to_string(sm_func_mode)+"/sigma_"+to_string_with_precision(sigma_list[is],3)+".jack", 0,1));
	AkAk_OS.distr_list.emplace_back(UseJack, Read_From_File("../data/tau_decay/"+Tag_reco_type+"/light/jackknife/OS/AkAk/"+Vk_data_tm.Tag[iens]+"/sm_func_mode_"+to_string(sm_func_mode)+"/sigma_"+to_string_with_precision(sigma_list[is],3)+".jack", 0,1));

	
	//load VkVk jack distribution
	VkVk_tm.distr_list.emplace_back(UseJack, Read_From_File("../data/tau_decay/"+Tag_reco_type+"/light/jackknife/tm/VkVk/"+Vk_data_tm.Tag[iens]+"/sm_func_mode_"+to_string(sm_func_mode)+"/sigma_"+to_string_with_precision(sigma_list[is],3)+".jack", 0,1));
	VkVk_OS.distr_list.emplace_back(UseJack, Read_From_File("../data/tau_decay/"+Tag_reco_type+"/light/jackknife/OS/VkVk/"+Vk_data_tm.Tag[iens]+"/sm_func_mode_"+to_string(sm_func_mode)+"/sigma_"+to_string_with_precision(sigma_list[is],3)+".jack", 0,1));

	
	//generate fake jack-distributions corresponding to systematic errors
	distr_t distr_syst_A0A0_tm(UseJack), distr_syst_A0A0_OS(UseJack);
	distr_t distr_syst_AkAk_tm(UseJack), distr_syst_AkAk_OS(UseJack);
	distr_t distr_syst_VkVk_tm(UseJack), distr_syst_VkVk_OS(UseJack);
	
	for(int ijack=0;ijack<Njacks;ijack++) {
	  //A0A0
	  distr_syst_A0A0_tm.distr.push_back( GM_sigma()*syst_A0A0_tm[iens][is]/sqrt(Njacks-1.0));
	  distr_syst_A0A0_OS.distr.push_back( GM_sigma()*syst_A0A0_OS[iens][is]/sqrt(Njacks-1.0));
	  //AkAk
	  distr_syst_AkAk_tm.distr.push_back( GM_sigma()*syst_AkAk_tm[iens][is]/sqrt(Njacks-1.0));
	  distr_syst_AkAk_OS.distr.push_back( GM_sigma()*syst_AkAk_OS[iens][is]/sqrt(Njacks-1.0));
	  //VkVk
	  distr_syst_VkVk_tm.distr.push_back( GM_sigma()*syst_VkVk_tm[iens][is]/sqrt(Njacks-1.0));
	  distr_syst_VkVk_OS.distr.push_back( GM_sigma()*syst_VkVk_OS[iens][is]/sqrt(Njacks-1.0));
	}
	//A0A0
	A0A0_tm.distr_list[iens] = A0A0_tm.distr_list[iens] + distr_syst_A0A0_tm;
	A0A0_OS.distr_list[iens] = A0A0_OS.distr_list[iens] + distr_syst_A0A0_OS;
	//AkAk
	AkAk_tm.distr_list[iens] = AkAk_tm.distr_list[iens] + distr_syst_AkAk_tm;
	AkAk_OS.distr_list[iens] = AkAk_OS.distr_list[iens] + distr_syst_AkAk_OS;
	//VkVk
	VkVk_tm.distr_list[iens] = VkVk_tm.distr_list[iens] + distr_syst_VkVk_tm;
	VkVk_OS.distr_list[iens] = VkVk_OS.distr_list[iens] + distr_syst_VkVk_OS;
	
      }

   
      //push back to all_s
      A0A0_tm_all_s[is] = A0A0_tm;
      A0A0_OS_all_s[is] = A0A0_OS;
      AkAk_tm_all_s[is] = AkAk_tm;
      AkAk_OS_all_s[is] = AkAk_OS;
      VkVk_tm_all_s[is] = VkVk_tm;
      VkVk_OS_all_s[is] = VkVk_OS;

      distr_t_list a_distr_list(UseJack);


      //estimate FSEs at fixed sigma from difference between B96 and B64


      distr_t_list A0A0_tm_red(UseJack), A0A0_OS_red(UseJack);
      distr_t_list AkAk_tm_red(UseJack), AkAk_OS_red(UseJack);
      distr_t_list VkVk_tm_red(UseJack), VkVk_OS_red(UseJack);

      distr_t_list a_distr_list_red(UseJack);
      vector<string> Tag_ens_red;

      
       
      //###################################################
      //############## COMPUTE FSE ########################
      double corr_fact_A0A0_FSE, corr_fact_A0A0_FSE_tm, corr_fact_A0A0_FSE_OS;
      double corr_fact_AkAk_FSE, corr_fact_AkAk_FSE_tm, corr_fact_AkAk_FSE_OS;
      double corr_fact_VkVk_FSE, corr_fact_VkVk_FSE_tm, corr_fact_VkVk_FSE_OS;
      int id_B64=-1;
      int id_B96=-1;

      for(int iens=0;iens<Nens;iens++) {
	if(Vk_data_tm.Tag[iens] == "cB211b.072.64") { id_B64=iens;}
	else if(Vk_data_tm.Tag[iens] == "cB211b.072.96") { id_B96=iens;}
      }

      if( (id_B64==-1) || (id_B96==-1)) crash("Cannot find id_B64 and/or id_B96");

      if(id_B64 == id_B96) crash("Error: id_B64 == id_B96");

      //A0A0
      corr_fact_A0A0_FSE_tm = fabs(((A0A0_tm.ave(id_B96) - A0A0_tm.ave(id_B64))/(A0A0_tm.ave(id_B64)))*erf( (A0A0_tm.ave(id_B96)-A0A0_tm.ave(id_B64))/(sqrt(2.0*( pow( A0A0_tm.err(id_B96)  ,2)  + pow( A0A0_tm.err(id_B64) ,2)   )))));
      corr_fact_A0A0_FSE_OS = fabs(((A0A0_OS.ave(id_B96) - A0A0_OS.ave(id_B64))/(A0A0_OS.ave(id_B64)))*erf( (A0A0_OS.ave(id_B96)-A0A0_OS.ave(id_B64))/(sqrt(2.0*( pow( A0A0_OS.err(id_B96)  ,2)  + pow( A0A0_OS.err(id_B64) ,2)   )))));
      corr_fact_A0A0_FSE = max(corr_fact_A0A0_FSE_tm, corr_fact_A0A0_FSE_OS);
      distr_t distr_syst_FSE_A0A0(UseJack);
      

      //AkAk
      corr_fact_AkAk_FSE_tm = fabs(((AkAk_tm.ave(id_B96) - AkAk_tm.ave(id_B64))/(AkAk_tm.ave(id_B64)))*erf( (AkAk_tm.ave(id_B96)-AkAk_tm.ave(id_B64))/(sqrt(2.0*( pow( AkAk_tm.err(id_B96)  ,2)  + pow( AkAk_tm.err(id_B64) ,2)   )))));
      corr_fact_AkAk_FSE_OS = fabs(((AkAk_OS.ave(id_B96) - AkAk_OS.ave(id_B64))/(AkAk_OS.ave(id_B64)))*erf( (AkAk_OS.ave(id_B96)-AkAk_OS.ave(id_B64))/(sqrt(2.0*( pow( AkAk_OS.err(id_B96)  ,2)  + pow( AkAk_OS.err(id_B64) ,2)   )))));
      corr_fact_AkAk_FSE = max(corr_fact_AkAk_FSE_tm, corr_fact_AkAk_FSE_OS);
      distr_t distr_syst_FSE_AkAk(UseJack);

      

      //VkVk
      corr_fact_VkVk_FSE_tm = fabs(((VkVk_tm.ave(id_B96) - VkVk_tm.ave(id_B64))/(VkVk_tm.ave(id_B64)))*erf( (VkVk_tm.ave(id_B96)-VkVk_tm.ave(id_B64))/(sqrt(2.0*( pow( VkVk_tm.err(id_B96)  ,2)  + pow( VkVk_tm.err(id_B64) ,2)   )))));
      corr_fact_VkVk_FSE_OS = fabs(((VkVk_OS.ave(id_B96) - VkVk_OS.ave(id_B64))/(VkVk_OS.ave(id_B64)))*erf( (VkVk_OS.ave(id_B96)-VkVk_OS.ave(id_B64))/(sqrt(2.0*( pow( VkVk_OS.err(id_B96)  ,2)  + pow( VkVk_OS.err(id_B64) ,2)   )))));
      corr_fact_VkVk_FSE = max(corr_fact_VkVk_FSE_tm, corr_fact_VkVk_FSE_OS);
      cout<<"corr_facts (tm, OS, max): "<<corr_fact_VkVk_FSE_tm<<" "<<corr_fact_VkVk_FSE_OS<<" "<<corr_fact_VkVk_FSE<<endl;
      distr_t distr_syst_FSE_VkVk(UseJack);

      for(int ijack=0;ijack<Njacks;ijack++) {
	distr_syst_FSE_A0A0.distr.push_back( 1.0 + GM_sigma()*corr_fact_A0A0_FSE/sqrt(Njacks-1.0));
	distr_syst_FSE_AkAk.distr.push_back( 1.0 + GM_sigma()*corr_fact_AkAk_FSE/sqrt(Njacks-1.0));
	distr_syst_FSE_VkVk.distr.push_back( 1.0 + GM_sigma()*corr_fact_VkVk_FSE/sqrt(Njacks-1.0));
      }
      
      for(int iens=0;iens<Nens;iens++) {

	if(Vk_data_tm.Tag[iens] != "cB211b.072.96") {
	  //A0A0 
	  A0A0_tm_red.distr_list.push_back( A0A0_tm.distr_list[iens]*distr_syst_FSE_A0A0);
	  A0A0_OS_red.distr_list.push_back( A0A0_OS.distr_list[iens]*distr_syst_FSE_A0A0);
	  //AkAk
	  AkAk_tm_red.distr_list.push_back( AkAk_tm.distr_list[iens]*distr_syst_FSE_AkAk);
	  AkAk_OS_red.distr_list.push_back( AkAk_OS.distr_list[iens]*distr_syst_FSE_AkAk);
	  //VkVk
	  VkVk_tm_red.distr_list.push_back( VkVk_tm.distr_list[iens]*distr_syst_FSE_VkVk);
	  VkVk_OS_red.distr_list.push_back( VkVk_OS.distr_list[iens]*distr_syst_FSE_VkVk);
	  if(Vk_data_tm.Tag[iens].substr(1,1) == "B") a_distr_list_red.distr_list.push_back( a_B/fm_to_inv_Gev); //lattice spacing is in fm
	  else if(Vk_data_tm.Tag[iens].substr(1,1) == "C") a_distr_list_red.distr_list.push_back( a_C/fm_to_inv_Gev); //lattice spacing is in fm
	  else if(Vk_data_tm.Tag[iens].substr(1,1) == "D") a_distr_list_red.distr_list.push_back( a_D/fm_to_inv_Gev); //lattice spacing is in fm
	  else crash("While building a_distr_list_red cannot recognize ensemble: "+Vk_data_tm.Tag[iens]);
	
	  Tag_ens_red.push_back( Vk_data_tm.Tag[iens]);

	}

	if(Vk_data_tm.Tag[iens].substr(1,1) == "B") a_distr_list.distr_list.push_back(a_B/fm_to_inv_Gev);  //lattice spacing is in fm
	else if(Vk_data_tm.Tag[iens].substr(1,1) == "C") a_distr_list.distr_list.push_back( a_C/fm_to_inv_Gev); //lattice spacing is in fm
	else if(Vk_data_tm.Tag[iens].substr(1,1) == "D") a_distr_list.distr_list.push_back( a_D/fm_to_inv_Gev); //lattice spacing is in fm
	else crash("When building a_distr_list cannot recognize ensemble: "+Vk_data_tm.Tag[iens]);
	
      }
      //###################################################


      
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
	  if((signed)par.size() != 3) crash("In class fpar_TAU, class constructor Vfloat par has size != 3");
	  D=par[0];
	  D2_tm=par[1];
	  D2_OS=par[2];
	}

	double D,D2_tm, D2_OS;
      };

       for( auto &contr: Contribs) {
	for( auto &fit_type: Fit_types) {
	  for( auto &poly_type: poly_types) {

	    if( (fit_type != "comb") && (poly_type.substr(0,2)=="tm" || poly_type.substr(0,2) == "OS")) crash("Cannot use tm/OS_linear with fit type: "+fit_type);
	    cout<<"###########################################################"<<endl;
	    cout<<"Performing continuum limit extrapolation for sigma: "<<sigma_list[is]<<endl;
	    cout<<"Contribution: "<<contr<<endl;
	    cout<<"Fit type: "<<fit_type<<endl;
	    cout<<"polynomial: "<<poly_type<<endl;

	    //Depending on fit considered, determine Nmeas, Npars and Ndof
	    int Nmeas= ((fit_type=="comb")?(2*Nens_eff):Nens_eff);
	    int deg_pars= ((fit_type=="comb")?2:1);
	    int Npars= 1 + (poly_type=="linear")*deg_pars +(poly_type.substr(0,2)=="tm") + (poly_type.substr(0,2)=="OS");
	    int Ndof= Nmeas-Npars;

	    cout<<"Nmeas: "<<Nmeas<<endl;
	    cout<<"Npars: "<<Npars<<endl;
	    cout<<"Ndof: "<<Ndof<<endl;
	    cout<<"Nens: "<<Nens_eff<<endl;
	    
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
	    bf_TAU.Add_par("D2_OS", 2, 0.1);
	    //ch2
	    bf_TAU_ch2.Add_par("D", 3.0, 0.1);
	    bf_TAU_ch2.Add_par("D2_tm", 2, 0.1);
	    bf_TAU_ch2.Add_par("D2_OS", 2, 0.1);

	    //fix parameters depending on fit type
	    if(poly_type=="const") {
	      bf_TAU.Fix_par("D2_tm", 0);
	      bf_TAU.Fix_par("D2_OS", 0);
	      //ch2
	      bf_TAU_ch2.Fix_par("D2_tm", 0);
	      bf_TAU_ch2.Fix_par("D2_OS", 0);
	    }
	    else if(poly_type=="linear") {
	      if(fit_type=="OS") { bf_TAU.Fix_par("D2_tm", 0); bf_TAU_ch2.Fix_par("D2_tm", 0); }
	      else if(fit_type=="tm") { bf_TAU.Fix_par("D2_OS",0); bf_TAU_ch2.Fix_par("D2_OS",0); }
	    }
	    else if(poly_type=="tm_linear") {
	      bf_TAU.Fix_par("D2_OS",0); bf_TAU_ch2.Fix_par("D2_OS",0);
	    }
	    else if(poly_type=="OS_linear") {
	       bf_TAU.Fix_par("D2_tm",0); bf_TAU_ch2.Fix_par("D2_tm",0);
	    }
	    else crash("poly_type: "+poly_type+" not yet implemented");


	    //ansatz
	    bf_TAU.ansatz=  [ ](const fpar_TAU &p, const ipar_TAU &ip) {
	      double D2=0.0;
	      if( ip.Is_tm==true ) D2=p.D2_tm;
	      else D2=p.D2_OS;
	      
	      return p.D + D2*pow(ip.a*QCD_scale,2);
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
		else if(contr=="AkAk") {
		  Cov_Matrix(iens,iens) = pow(AkAk_tm_red.err(iens),2); Cov_Matrix(iens+off_OS,iens+off_OS) = pow(AkAk_OS_red.err(iens),2); 
		  Cov_Matrix(iens, off_OS+iens) = (AkAk_tm_red.distr_list[iens]%AkAk_OS_red.distr_list[iens]);
		}

		else if(contr=="VkVk") {
		  Cov_Matrix(iens,iens) = pow(VkVk_tm_red.err(iens),2); Cov_Matrix(iens+off_OS,iens+off_OS) = pow(VkVk_OS_red.err(iens),2);
		  Cov_Matrix(iens, off_OS+iens) = (VkVk_tm_red.distr_list[iens]%VkVk_OS_red.distr_list[iens]);
		}

		else if(contr=="VA") {
		  Cov_Matrix(iens,iens) = pow((VkVk_tm_red/(AkAk_tm_red+A0A0_tm_red)).err(iens),2); Cov_Matrix(iens+off_OS,iens+off_OS) = pow((VkVk_OS_red/(AkAk_OS_red+A0A0_OS_red)).err(iens),2);
		  Cov_Matrix(iens, off_OS+iens) = ((VkVk_tm_red/(AkAk_tm_red+A0A0_tm_red)).distr_list[iens]%(VkVk_OS_red/(AkAk_OS_red+A0A0_OS_red)).distr_list[iens]);
		}

		else if(contr=="VMA") {
		  Cov_Matrix(iens,iens) = pow((( VkVk_tm_red-AkAk_tm_red-A0A0_tm_red )/(VkVk_tm_red+AkAk_tm_red+A0A0_tm_red)).err(iens),2); Cov_Matrix(iens+off_OS,iens+off_OS) = pow( (( VkVk_OS_red-AkAk_OS_red-A0A0_OS_red )/(VkVk_OS_red+AkAk_OS_red+A0A0_OS_red)).err(iens),2);
		  Cov_Matrix(iens, off_OS+iens) = ((( VkVk_tm_red-AkAk_tm_red-A0A0_tm_red )/(VkVk_tm_red+AkAk_tm_red+A0A0_tm_red)).distr_list[iens]%(( VkVk_OS_red-AkAk_OS_red-A0A0_OS_red )/(VkVk_OS_red+AkAk_OS_red+A0A0_OS_red)).distr_list[iens]);
		}

		else if(contr=="AX") {
		  Cov_Matrix(iens,iens) = pow( (AkAk_tm_red+A0A0_tm_red).err(iens),2); Cov_Matrix(iens+off_OS,iens+off_OS) = pow( (AkAk_OS_red+ A0A0_OS_red).err(iens),2); 
		  Cov_Matrix(iens, off_OS+iens) = ( (AkAk_tm_red+A0A0_tm_red).distr_list[iens]% (AkAk_OS_red+A0A0_OS_red).distr_list[iens]);
		}

		else if(contr=="T") {
		  Cov_Matrix(iens,iens) = pow( (AkAk_tm_red+VkVk_tm_red).err(iens),2); Cov_Matrix(iens+off_OS,iens+off_OS) = pow( (AkAk_OS_red+ VkVk_OS_red).err(iens),2); 
		  Cov_Matrix(iens, off_OS+iens) = ( (AkAk_tm_red+VkVk_tm_red).distr_list[iens]% (AkAk_OS_red+VkVk_OS_red).distr_list[iens]);
		}
		

		else if(contr=="tot") {
		  Cov_Matrix(iens,iens) = pow((VkVk_tm_red+AkAk_tm_red+A0A0_tm_red).err(iens),2); Cov_Matrix(iens+off_OS,iens+off_OS) = pow( (VkVk_OS_red+AkAk_OS_red+A0A0_OS_red).err(iens),2);
		  Cov_Matrix(iens, off_OS+iens) = ((VkVk_tm_red+AkAk_tm_red+A0A0_tm_red).distr_list[iens]%(VkVk_OS_red+AkAk_OS_red+A0A0_OS_red).distr_list[iens]);
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
	      ofstream Print_Cov("../data/tau_decay/"+Tag_reco_type+"/light/continuum/cov/"+contr+"_sigma_"+to_string_with_precision(sigma_list[is],3)+".cov");
	      ofstream Print_Corr("../data/tau_decay/"+Tag_reco_type+"/light/continuum/corr/"+contr+"_sigma_"+to_string_with_precision(sigma_list[is],3)+".corr");

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
		    data[ijack][iens].Br = (VkVk_tm_red/(AkAk_tm_red+A0A0_tm_red)).distr_list[iens].distr[ijack];
		    data[ijack][iens].Br_err= (VkVk_tm_red/(AkAk_tm_red+A0A0_tm_red)).err(iens);
		    data[ijack][iens].Is_tm = true;
		    data[ijack][iens].a = a_distr_list_red.distr_list[iens].distr[ijack];
		  }
		  
		  if((fit_type=="OS") || (fit_type=="comb")) {
		    data[ijack][iens+off_OS].Br = (VkVk_OS_red/(AkAk_OS_red+A0A0_OS_red)).distr_list[iens].distr[ijack];
		    data[ijack][iens+off_OS].Br_err= (VkVk_OS_red/(AkAk_OS_red+A0A0_OS_red)).err(iens);
		    data[ijack][iens+off_OS].Is_tm = false;
		    data[ijack][iens+off_OS].a = a_distr_list_red.distr_list[iens].distr[ijack];
		  }
		  
		}

		else if(contr=="VMA") {
		  if((fit_type=="tm") || (fit_type=="comb")) {
		    data[ijack][iens].Br = ( ( VkVk_tm_red-AkAk_tm_red-A0A0_tm_red )/(VkVk_tm_red+AkAk_tm_red+A0A0_tm_red)).distr_list[iens].distr[ijack];
		    data[ijack][iens].Br_err= ( ( VkVk_tm_red-AkAk_tm_red-A0A0_tm_red )/(VkVk_tm_red+AkAk_tm_red+A0A0_tm_red)).err(iens);
		    data[ijack][iens].Is_tm = true;
		    data[ijack][iens].a = a_distr_list_red.distr_list[iens].distr[ijack];
		  }

		  if((fit_type=="OS") || (fit_type=="comb")) {
		    data[ijack][iens+off_OS].Br = ( ( VkVk_OS_red-AkAk_OS_red-A0A0_OS_red )/(VkVk_OS_red+AkAk_OS_red+A0A0_OS_red)).distr_list[iens].distr[ijack];
		    data[ijack][iens+off_OS].Br_err=  ( ( VkVk_OS_red-AkAk_OS_red-A0A0_OS_red )/(VkVk_OS_red+AkAk_OS_red+A0A0_OS_red)).err(iens);
		    data[ijack][iens+off_OS].Is_tm = false;
		    data[ijack][iens+off_OS].a = a_distr_list_red.distr_list[iens].distr[ijack];
		  }
		 
		}

		else if(contr=="T") {
		   if((fit_type=="tm") || (fit_type=="comb")) {
		     data[ijack][iens].Br = (AkAk_tm_red+VkVk_tm_red).distr_list[iens].distr[ijack];
		     data[ijack][iens].Br_err= (AkAk_tm_red+VkVk_tm_red).err(iens);
		     data[ijack][iens].Is_tm = true;
		     data[ijack][iens].a = a_distr_list_red.distr_list[iens].distr[ijack];
		   }

		   if((fit_type=="OS") || (fit_type=="comb")) {
		    data[ijack][iens+off_OS].Br = (AkAk_OS_red+VkVk_OS_red).distr_list[iens].distr[ijack];
		    data[ijack][iens+off_OS].Br_err= (AkAk_OS_red+VkVk_OS_red).err(iens);
		    data[ijack][iens+off_OS].Is_tm = false;
		    data[ijack][iens+off_OS].a = a_distr_list_red.distr_list[iens].distr[ijack];
		   }
		   
		}

		

				
		else if(contr=="tot") {
		  if((fit_type=="tm") || (fit_type=="comb")) {
		    data[ijack][iens].Br = (VkVk_tm_red+AkAk_tm_red+A0A0_tm_red).distr_list[iens].distr[ijack];
		    data[ijack][iens].Br_err= (VkVk_tm_red+AkAk_tm_red+A0A0_tm_red).err(iens);
		    data[ijack][iens].Is_tm = true;
		    data[ijack][iens].a = a_distr_list_red.distr_list[iens].distr[ijack];
		  }

		  if((fit_type=="OS") || (fit_type=="comb")) {
		    data[ijack][iens+off_OS].Br = (VkVk_OS_red+AkAk_OS_red+A0A0_OS_red).distr_list[iens].distr[ijack];
		    data[ijack][iens+off_OS].Br_err= (VkVk_OS_red+AkAk_OS_red+A0A0_OS_red).err(iens);
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
		      data_ch2[ijack][iens].Br =  (VkVk_tm_red/(AkAk_tm_red+A0A0_tm_red)).ave(iens);
		      data_ch2[ijack][iens].Br_err=  (VkVk_tm_red/(AkAk_tm_red+A0A0_tm_red)).err(iens);
		      data_ch2[ijack][iens].Is_tm = true;
		      data_ch2[ijack][iens].a = a_distr_list_red.ave(iens);
		    }

		    if((fit_type=="OS") || (fit_type=="comb")) {
		      data_ch2[ijack][iens+off_OS].Br =  (VkVk_OS_red/(AkAk_OS_red+A0A0_OS_red)).ave(iens);
		      data_ch2[ijack][iens+off_OS].Br_err=  (VkVk_OS_red/(AkAk_OS_red+A0A0_OS_red)).err(iens);
		      data_ch2[ijack][iens+off_OS].Is_tm = false;
		      data_ch2[ijack][iens+off_OS].a = a_distr_list_red.ave(iens);
		    }
		    
		  }

		  else if(contr=="VMA") {
		    if((fit_type=="tm") || (fit_type=="comb")) {
		      data_ch2[ijack][iens].Br = ((VkVk_tm_red-AkAk_tm_red-A0A0_tm_red)/(VkVk_tm_red+AkAk_tm_red+A0A0_tm_red)).ave(iens);
		      data_ch2[ijack][iens].Br_err=  ((VkVk_tm_red-AkAk_tm_red-A0A0_tm_red)/(VkVk_tm_red+AkAk_tm_red+A0A0_tm_red)).err(iens);
		      data_ch2[ijack][iens].Is_tm = true;
		      data_ch2[ijack][iens].a = a_distr_list_red.ave(iens);
		    }
		    
		    if((fit_type=="OS") || (fit_type=="comb")) {
		      data_ch2[ijack][iens+off_OS].Br =  ((VkVk_OS_red-AkAk_OS_red-A0A0_OS_red)/(VkVk_OS_red+AkAk_OS_red+A0A0_OS_red)).ave(iens);
		      data_ch2[ijack][iens+off_OS].Br_err=  ((VkVk_OS_red-AkAk_OS_red-A0A0_OS_red)/(VkVk_OS_red+AkAk_OS_red+A0A0_OS_red)).err(iens);
		      data_ch2[ijack][iens+off_OS].Is_tm = false;
		      data_ch2[ijack][iens+off_OS].a = a_distr_list_red.ave(iens);
		    }
		    
		  }

		  else if(contr=="T") {
		    if((fit_type=="tm") || (fit_type=="comb")) {
		      data_ch2[ijack][iens].Br = (AkAk_tm_red+VkVk_tm_red).ave(iens);
		      data_ch2[ijack][iens].Br_err= (AkAk_tm_red+VkVk_tm_red).err(iens);
		      data_ch2[ijack][iens].Is_tm = true;
		      data_ch2[ijack][iens].a = a_distr_list_red.ave(iens);
		    }

		    if((fit_type=="OS") || (fit_type=="comb")) {
		      data_ch2[ijack][iens+off_OS].Br = (AkAk_OS_red+VkVk_OS_red).ave(iens);
		      data_ch2[ijack][iens+off_OS].Br_err= (AkAk_OS_red+VkVk_OS_red).err(iens);
		      data_ch2[ijack][iens+off_OS].Is_tm = false;
		      data_ch2[ijack][iens+off_OS].a = a_distr_list_red.ave(iens);
		    }
		    
		  }
		  
		  else if(contr=="tot") {
		    if((fit_type=="tm") || (fit_type=="comb")) {
		      data_ch2[ijack][iens].Br = (VkVk_tm_red+AkAk_tm_red+A0A0_tm_red).ave(iens);
		      data_ch2[ijack][iens].Br_err= (VkVk_tm_red+AkAk_tm_red+A0A0_tm_red).err(iens);
		      data_ch2[ijack][iens].Is_tm = true;
		      data_ch2[ijack][iens].a = a_distr_list_red.ave(iens);
		    }
		    
		    if((fit_type=="OS") || (fit_type=="comb")) {
		      data_ch2[ijack][iens+off_OS].Br = (VkVk_OS_red+AkAk_OS_red+A0A0_OS_red).ave(iens);
		      data_ch2[ijack][iens+off_OS].Br_err= (VkVk_OS_red+AkAk_OS_red+A0A0_OS_red).err(iens);
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
	    distr_t D(UseJack), D2_tm(UseJack), D2_OS(UseJack);
	    for(int ijack=0;ijack<Njacks;ijack++) { D.distr.push_back( Bt_fit.par[ijack].D); D2_tm.distr.push_back( Bt_fit.par[ijack].D2_tm); D2_OS.distr.push_back( Bt_fit.par[ijack].D2_OS);}
	    //reduced ch2
	    double ch2= Bt_fit_ch2.get_ch2_ave()/Ndof;


	    //print fit function
	    distr_t_list FF_tm_to_print(UseJack), FF_OS_to_print(UseJack);

	    
	    for(auto &a: a_to_print) {
	      FF_tm_to_print.distr_list.push_back( D + D2_tm*pow(a*QCD_scale,2));
	      FF_OS_to_print.distr_list.push_back( D + D2_OS*pow(a*QCD_scale,2));
	    }
	    string Fit_tag= "../data/tau_decay/"+Tag_reco_type+"/light/continuum/fit_func/sigma_"+to_string_with_precision(sigma_list[is],3)+"_contr_"+contr+"_ftype_"+fit_type+"_"+poly_type+".dat";
	    Print_To_File({}, {a_to_print, FF_tm_to_print.ave(), FF_tm_to_print.err(), FF_OS_to_print.ave(), FF_OS_to_print.err()},Fit_tag, "", "#a[fm] tm OS, ch2/dof: "+to_string_with_precision(ch2,4));
	    

	    //push back information on ch2  and on fit result

	    res_map.find({contr, fit_type, poly_type})->second.distr_list.push_back(D);
	    ch2_map.find({contr, fit_type, poly_type})->second.push_back(ch2);
	    Nmeas_map.find({contr, fit_type, poly_type})->second.push_back(Nmeas);
	    Ndof_map.find({contr, fit_type, poly_type})->second.push_back(Ndof);
	    



	    
	  }
	}
       }

       //print data
       //all
       //tm
       Print_To_File({Vk_data_tm.Tag}, { a_distr_list.ave(), A0A0_tm.ave(), A0A0_tm.err() } ,  "../data/tau_decay/"+Tag_reco_type+"/light/continuum/fit_data/A0A0/tm/sigma_"+to_string_with_precision(sigma_list[is],3)+".dat"  , "", "");
       Print_To_File({Vk_data_tm.Tag}, { a_distr_list.ave(), AkAk_tm.ave(), AkAk_tm.err() } ,  "../data/tau_decay/"+Tag_reco_type+"/light/continuum/fit_data/AkAk/tm/sigma_"+to_string_with_precision(sigma_list[is],3)+".dat"  , "", "");
       Print_To_File({Vk_data_tm.Tag}, { a_distr_list.ave(), (AkAk_tm+A0A0_tm).ave(), (AkAk_tm+A0A0_tm).err() } ,  "../data/tau_decay/"+Tag_reco_type+"/light/continuum/fit_data/AX/tm/sigma_"+to_string_with_precision(sigma_list[is],3)+".dat"  , "", "");
       Print_To_File({Vk_data_tm.Tag}, { a_distr_list.ave(), (AkAk_tm+VkVk_tm).ave(), (AkAk_tm+VkVk_tm).err() } ,  "../data/tau_decay/"+Tag_reco_type+"/light/continuum/fit_data/T/tm/sigma_"+to_string_with_precision(sigma_list[is],3)+".dat"  , "", "");
       Print_To_File({Vk_data_tm.Tag}, { a_distr_list.ave(), VkVk_tm.ave(), VkVk_tm.err() } ,  "../data/tau_decay/"+Tag_reco_type+"/light/continuum/fit_data/VkVk/tm/sigma_"+to_string_with_precision(sigma_list[is],3)+".dat"  , "", "");
       Print_To_File({Vk_data_tm.Tag}, { a_distr_list.ave(), (VkVk_tm+AkAk_tm+A0A0_tm).ave(), (VkVk_tm+AkAk_tm+A0A0_tm).err() } ,  "../data/tau_decay/"+Tag_reco_type+"/light/continuum/fit_data/tot/tm/sigma_"+to_string_with_precision(sigma_list[is],3)+".dat"  , "", "");
       Print_To_File({Vk_data_tm.Tag}, { a_distr_list.ave(), ((VkVk_tm/(AkAk_tm+A0A0_tm))).ave(), ((VkVk_tm/(AkAk_tm+A0A0_tm))).err() } ,  "../data/tau_decay/"+Tag_reco_type+"/light/continuum/fit_data/VA/tm/sigma_"+to_string_with_precision(sigma_list[is],3)+".dat"  , "", "");
       Print_To_File({Vk_data_tm.Tag}, { a_distr_list.ave(), (((VkVk_tm-AkAk_tm-A0A0_tm)/(VkVk_tm+AkAk_tm+A0A0_tm))).ave(), (((VkVk_tm-AkAk_tm-A0A0_tm)/(VkVk_tm+AkAk_tm+A0A0_tm))).err() } ,  "../data/tau_decay/"+Tag_reco_type+"/light/continuum/fit_data/VMA/tm/sigma_"+to_string_with_precision(sigma_list[is],3)+".dat"  , "", "");
       //OS
       Print_To_File({Vk_data_OS.Tag}, { a_distr_list.ave(), A0A0_OS.ave(), A0A0_OS.err() } ,  "../data/tau_decay/"+Tag_reco_type+"/light/continuum/fit_data/A0A0/OS/sigma_"+to_string_with_precision(sigma_list[is],3)+".dat"  , "", "");
       Print_To_File({Vk_data_OS.Tag}, { a_distr_list.ave(), AkAk_OS.ave(), AkAk_OS.err() } ,  "../data/tau_decay/"+Tag_reco_type+"/light/continuum/fit_data/AkAk/OS/sigma_"+to_string_with_precision(sigma_list[is],3)+".dat"  , "", "");
       Print_To_File({Vk_data_OS.Tag}, { a_distr_list.ave(), (AkAk_OS+A0A0_OS).ave(), (AkAk_OS+A0A0_OS).err() } ,  "../data/tau_decay/"+Tag_reco_type+"/light/continuum/fit_data/AX/OS/sigma_"+to_string_with_precision(sigma_list[is],3)+".dat"  , "", "");
       Print_To_File({Vk_data_OS.Tag}, { a_distr_list.ave(), (AkAk_OS+VkVk_OS).ave(), (AkAk_OS+VkVk_OS).err() } ,  "../data/tau_decay/"+Tag_reco_type+"/light/continuum/fit_data/T/OS/sigma_"+to_string_with_precision(sigma_list[is],3)+".dat"  , "", "");
       Print_To_File({Vk_data_OS.Tag}, { a_distr_list.ave(), VkVk_OS.ave(), VkVk_OS.err() } ,  "../data/tau_decay/"+Tag_reco_type+"/light/continuum/fit_data/VkVk/OS/sigma_"+to_string_with_precision(sigma_list[is],3)+".dat"  , "", "");
       Print_To_File({Vk_data_OS.Tag}, { a_distr_list.ave(), (VkVk_OS+AkAk_OS+A0A0_OS).ave(), (VkVk_OS+AkAk_OS+A0A0_OS).err() } ,  "../data/tau_decay/"+Tag_reco_type+"/light/continuum/fit_data/tot/OS/sigma_"+to_string_with_precision(sigma_list[is],3)+".dat"  , "", "");
       Print_To_File({Vk_data_OS.Tag}, { a_distr_list.ave(), ((VkVk_OS/(AkAk_OS+A0A0_OS))).ave(), ((VkVk_OS/(AkAk_OS+A0A0_OS))).err() } ,  "../data/tau_decay/"+Tag_reco_type+"/light/continuum/fit_data/VA/OS/sigma_"+to_string_with_precision(sigma_list[is],3)+".dat"  , "", "");
       Print_To_File({Vk_data_OS.Tag}, { a_distr_list.ave(), (((VkVk_OS-AkAk_OS-A0A0_OS)/(VkVk_OS+AkAk_OS+A0A0_OS))).ave(), (((VkVk_OS-AkAk_OS-A0A0_OS)/(VkVk_OS+AkAk_OS+A0A0_OS))).err() } ,  "../data/tau_decay/"+Tag_reco_type+"/light/continuum/fit_data/VMA/OS/sigma_"+to_string_with_precision(sigma_list[is],3)+".dat"  , "", "");

       
       //red
       Print_To_File({Tag_ens_red}, { a_distr_list_red.ave(), A0A0_tm_red.ave(), A0A0_tm_red.err() } ,  "../data/tau_decay/"+Tag_reco_type+"/light/continuum/fit_data/A0A0/tm/FSE_corr_sigma_"+to_string_with_precision(sigma_list[is],3)+".dat"  , "", "");
       Print_To_File({Tag_ens_red}, { a_distr_list_red.ave(), AkAk_tm_red.ave(), AkAk_tm_red.err() } ,  "../data/tau_decay/"+Tag_reco_type+"/light/continuum/fit_data/AkAk/tm/FSE_corr_sigma_"+to_string_with_precision(sigma_list[is],3)+".dat"  , "", "");
       Print_To_File({Tag_ens_red}, { a_distr_list_red.ave(), (AkAk_tm_red+A0A0_tm_red).ave(), (AkAk_tm_red+A0A0_tm_red).err() } ,  "../data/tau_decay/"+Tag_reco_type+"/light/continuum/fit_data/AX/tm/FSE_corr_sigma_"+to_string_with_precision(sigma_list[is],3)+".dat"  , "", "");
       Print_To_File({Tag_ens_red}, { a_distr_list_red.ave(), (AkAk_tm_red+VkVk_tm_red).ave(), (AkAk_tm_red+VkVk_tm_red).err() } ,  "../data/tau_decay/"+Tag_reco_type+"/light/continuum/fit_data/T/tm/FSE_corr_sigma_"+to_string_with_precision(sigma_list[is],3)+".dat"  , "", "");
       Print_To_File({Tag_ens_red}, { a_distr_list_red.ave(), VkVk_tm_red.ave(), VkVk_tm_red.err() } ,  "../data/tau_decay/"+Tag_reco_type+"/light/continuum/fit_data/VkVk/tm/FSE_corr_sigma_"+to_string_with_precision(sigma_list[is],3)+".dat"  , "", "");
       Print_To_File({Tag_ens_red}, { a_distr_list_red.ave(), (VkVk_tm_red+AkAk_tm_red+A0A0_tm_red).ave(), (VkVk_tm_red+AkAk_tm_red+A0A0_tm_red).err() } ,  "../data/tau_decay/"+Tag_reco_type+"/light/continuum/fit_data/tot/tm/FSE_corr_sigma_"+to_string_with_precision(sigma_list[is],3)+".dat"  , "", "");
       Print_To_File({Tag_ens_red}, { a_distr_list_red.ave(), ((VkVk_tm_red/(AkAk_tm_red+A0A0_tm_red))).ave(), ((VkVk_tm_red/(AkAk_tm_red+A0A0_tm_red))).err() } ,  "../data/tau_decay/"+Tag_reco_type+"/light/continuum/fit_data/VA/tm/FSE_corr_sigma_"+to_string_with_precision(sigma_list[is],3)+".dat"  , "", "");
       Print_To_File({Tag_ens_red}, { a_distr_list_red.ave(), (((VkVk_tm_red -AkAk_tm_red - A0A0_tm_red )/(VkVk_tm_red+ AkAk_tm_red+A0A0_tm_red))).ave(), (((VkVk_tm_red -AkAk_tm_red - A0A0_tm_red )/(VkVk_tm_red+ AkAk_tm_red+A0A0_tm_red))).err() } ,  "../data/tau_decay/"+Tag_reco_type+"/light/continuum/fit_data/VMA/tm/FSE_corr_sigma_"+to_string_with_precision(sigma_list[is],3)+".dat"  , "", "");
       //OS
       Print_To_File({Tag_ens_red}, { a_distr_list_red.ave(), A0A0_OS_red.ave(), A0A0_OS_red.err() } ,  "../data/tau_decay/"+Tag_reco_type+"/light/continuum/fit_data/A0A0/OS/FSE_corr_sigma_"+to_string_with_precision(sigma_list[is],3)+".dat"  , "", "");
       Print_To_File({Tag_ens_red}, { a_distr_list_red.ave(), AkAk_OS_red.ave(), AkAk_OS_red.err() } ,  "../data/tau_decay/"+Tag_reco_type+"/light/continuum/fit_data/AkAk/OS/FSE_corr_sigma_"+to_string_with_precision(sigma_list[is],3)+".dat"  , "", "");
       Print_To_File({Tag_ens_red}, { a_distr_list_red.ave(), (AkAk_OS_red+A0A0_OS_red).ave(), (AkAk_OS_red+A0A0_OS_red).err() } ,  "../data/tau_decay/"+Tag_reco_type+"/light/continuum/fit_data/AX/OS/FSE_corr_sigma_"+to_string_with_precision(sigma_list[is],3)+".dat"  , "", "");
       Print_To_File({Tag_ens_red}, { a_distr_list_red.ave(), (AkAk_OS_red+VkVk_OS_red).ave(), (AkAk_OS_red+VkVk_OS_red).err() } ,  "../data/tau_decay/"+Tag_reco_type+"/light/continuum/fit_data/T/OS/FSE_corr_sigma_"+to_string_with_precision(sigma_list[is],3)+".dat"  , "", "");
       Print_To_File({Tag_ens_red}, { a_distr_list_red.ave(), VkVk_OS_red.ave(), VkVk_OS_red.err() } ,  "../data/tau_decay/"+Tag_reco_type+"/light/continuum/fit_data/VkVk/OS/FSE_corr_sigma_"+to_string_with_precision(sigma_list[is],3)+".dat"  , "", "");
       Print_To_File({Tag_ens_red}, { a_distr_list_red.ave(), (VkVk_OS_red+AkAk_OS_red+A0A0_OS_red).ave(), (VkVk_OS_red+AkAk_OS_red+A0A0_OS_red).err() } ,  "../data/tau_decay/"+Tag_reco_type+"/light/continuum/fit_data/tot/OS/FSE_corr_sigma_"+to_string_with_precision(sigma_list[is],3)+".dat"  , "", "");
       Print_To_File({Tag_ens_red}, { a_distr_list_red.ave(), ((VkVk_OS_red/(AkAk_OS_red+A0A0_OS_red))).ave(), ((VkVk_OS_red/(AkAk_OS_red+A0A0_OS_red))).err() } ,  "../data/tau_decay/"+Tag_reco_type+"/light/continuum/fit_data/VA/OS/FSE_corr_sigma_"+to_string_with_precision(sigma_list[is],3)+".dat"  , "", "");
       Print_To_File({Tag_ens_red}, { a_distr_list_red.ave(), (((VkVk_OS_red -AkAk_OS_red - A0A0_OS_red )/(VkVk_OS_red+ AkAk_OS_red+A0A0_OS_red))).ave(), (((VkVk_OS_red -AkAk_OS_red - A0A0_OS_red )/(VkVk_OS_red+ AkAk_OS_red+A0A0_OS_red))).err() } ,  "../data/tau_decay/"+Tag_reco_type+"/light/continuum/fit_data/VMA/OS/FSE_corr_sigma_"+to_string_with_precision(sigma_list[is],3)+".dat"  , "", "");


       
      
       
    }
    

    //here I should print the result
    for (auto const& [tag,Br] : res_map) {
      string out_tag= get<0>(tag)+"_"+get<1>(tag)+"_"+get<2>(tag);
      Vfloat ch2 = ch2_map.find(tag)->second;
      Print_To_File({}, { sigma_list, Br.ave(), Br.err(), ch2 }, "../data/tau_decay/"+Tag_reco_type+"/light/continuum/Extr_"+out_tag+".dat", "", "#sigma ave err ch2/dof ");
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

    vector<distr_t_list> Br_finals;
    VVfloat Br_final_systs(8);
    for(int c=0; c<(signed)Contribs.size(); c++)  { Br_finals.emplace_back(UseJack, sigma_list.size()); Br_final_systs[c].resize(sigma_list.size(), 0);}

 
    for(int is=0; is < (signed)sigma_list.size();is++) {

      //create directory for AIC output
      boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/light/continuum/AIC/sigma_"+to_string_with_precision(sigma_list[is],3));
       
      int c=0;
      //loop over contributions
      for(auto &contr: Contribs ) {

	map< tuple<string,string> , double> AIC_map;
	Vfloat weight_list, ch2_list, Nmeas_list, Npars_list;
	vector<string> ftpt_list;
	
	double w_tot=0;

	distr_t_list Print_Res_partial(UseJack);
	
	for(auto &ftype: Fit_types)
	  for(auto &ptype: poly_types) {

	    
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
	    Print_Res_partial.distr_list.push_back( res_map.find({contr, ftype,ptype})->second.distr_list[is]);
	    
	    AIC_map.insert( { {ftype, ptype}, w });
	    w_tot += w;
		
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
	
	  

	
	//print details on AIC
	Print_To_File({ftpt_list}, {weight_list, ch2_list, Nmeas_list, Npars_list, Print_Res_partial.ave(), Print_Res_partial.err()  } , "../data/tau_decay/"+Tag_reco_type+"/light/continuum/AIC/sigma_"+to_string_with_precision(sigma_list[is],3)+"/"+contr+".out", "", "#ftpt   w(AIC)   ch2/dof  Nmeas  Npars");
	cout<<"generating histograms for "<<Tag_reco_type<<endl<<flush;
	cout<<"Contrib: "<<Contribs[c]<<endl;

	//generate histograms
	int Nboots_histo=1000;
	int k=50;
	Vfloat x_list, s_list;
	double x_min= 0;
	double x_max= 2*Br.ave();
	if( Contribs[c] == "VMA") { x_min= -2000*fabs(Br.ave()-Br.err()); x_max = 2000*fabs(Br.ave()+Br.err()) ;}
	

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
	Print_To_File({}, {x_list, hist, s_list}, "../data/tau_decay/"+Tag_reco_type+"/light/continuum/AIC/sigma_"+to_string_with_precision(sigma_list[is],3)+"/"+contr+"_histo.out", "", "# x histo[x] ");

	c++;
      }
    }

    

    //print results
    for(int c=0;c<(signed)Contribs.size();c++) {
      string contr= Contribs[c];
      Print_To_File({}, { sigma_list, Br_finals[c].ave(), Br_finals[c].err(), Br_final_systs[c]}, "../data/tau_decay/"+Tag_reco_type+"/light/continuum/Extr_AIC_"+contr+".dat", "", "#sigma ave err_stat  err_syst ");
      
    }

    
    //###############################################
    //###############################################
    //###############################################


    //perform sigma to zero extrapolation of the various contributions
    boost::filesystem::create_directory("../data/tau_decay/"+Tag_reco_type+"/light/sigma_extr");

    vector<distr_t> Final_extr_data;

    int sigma_to_exclude=3;
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

      int Nmeas= sigma_list.size() - sigma_to_exclude;

      bootstrap_fit<fpar_SIGMA,ipar_SIGMA> bf_SIGMA(Njacks);
      bootstrap_fit<fpar_SIGMA,ipar_SIGMA> bf_SIGMA_ch2(1);
      bf_SIGMA.Set_number_of_measurements(Nmeas);
      bf_SIGMA.Set_verbosity(1);
      //ch2
      bf_SIGMA_ch2.Set_number_of_measurements(Nmeas);
      bf_SIGMA_ch2.Set_verbosity(1);

      //add fit parameters
      bf_SIGMA.Add_par("D", 2.0, 0.1);
      bf_SIGMA.Add_par("D4", 2, 0.1);
      bf_SIGMA.Add_par("D6", 2, 0.1);
      //ch2
      bf_SIGMA_ch2.Add_par("D", 2.0, 0.1);
      bf_SIGMA_ch2.Add_par("D4", 2, 0.1);
      bf_SIGMA_ch2.Add_par("D6", 2, 0.1);

      bool Fix_D6=true;
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
      
      for(auto &data_iboot: data) data_iboot.resize(Nmeas);
      for(auto &data_iboot: data_ch2) data_iboot.resize(Nmeas);


      //correlated fit
      Eigen::MatrixXd Cov_Matrix(Nmeas,Nmeas);
      Eigen::MatrixXd Corr_Matrix(Nmeas,Nmeas);
      for(int i=0;i<Nmeas;i++)
	for(int j=0;j<Nmeas;j++) {
	  Cov_Matrix(i,j) = Br_finals[c].distr_list[i]%Br_finals[c].distr_list[j];
	  Corr_Matrix(i,j)= Cov_Matrix(i,j)/(Br_finals[c].err(i)*Br_finals[c].err(j));
	}
      
      //add cov matrix to bootstrap fit
      //bf_SIGMA.Add_covariance_matrix(Cov_Matrix);
      //bf_SIGMA_ch2.Add_covariance_matrix(Cov_Matrix);

      //print covariance matrix
      ofstream Print_Cov("../data/tau_decay/"+Tag_reco_type+"/light/sigma_extr/contr_"+Contribs[c]+".cov");
      ofstream Print_Corr("../data/tau_decay/"+Tag_reco_type+"/light/sigma_extr/contr_"+Contribs[c]+".corr");
      
      Print_Cov<<Cov_Matrix<<endl;  Print_Corr<<Corr_Matrix<<endl;
      Print_Cov.close();            Print_Corr.close();
      
      

      GaussianMersenne GS(224223); //GS(15431); //;
         
      for(int ijack=0;ijack<Njacks;ijack++) {

	double rn= GS();
	for(int is=0;is<Nmeas;is++) {

	  data[ijack][is].Br = Br_finals[c].distr_list[is].distr[ijack] + rn*Br_final_systs[c][is]/sqrt(Njacks-1.0);
	  data[ijack][is].Br_err = sqrt( pow(Br_finals[c].err(is),2) + pow(Br_final_systs[c][is],2));
	  data[ijack][is].sigma= sigma_list[is];

	  if(ijack==0) {

	    data_ch2[ijack][is].Br= Br_finals[c].ave(is);
	    data_ch2[ijack][is].Br_err = sqrt( pow(Br_finals[c].err(is),2) + pow(Br_final_systs[c][is],2));
	    data_ch2[ijack][is].sigma = sigma_list[is];

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
      int Ndof= sigma_list.size() - sigma_to_exclude - ( (Fix_D6==true)?2:3 );
      double ch2= Bt_fit_ch2.get_ch2_ave()/Ndof;

      Final_extr_data.push_back(D);
      
      
      //print fit function
      distr_t_list R_at_sigma_to_print(UseJack);

      Vfloat sigma_to_print;
      for(int iss=0;iss<1000;iss++) { sigma_to_print.push_back( iss*0.0002);}
      
      for(auto &ss: sigma_to_print) {
	R_at_sigma_to_print.distr_list.push_back( D + D4*pow(ss,4) + D6*pow(ss,6));
      }
      string Fit_tag= "../data/tau_decay/"+Tag_reco_type+"/light/sigma_extr/contr_"+Contribs[c];
      Print_To_File({}, {sigma_to_print, R_at_sigma_to_print.ave(), R_at_sigma_to_print.err()},Fit_tag, "", "#sigma Br Br_err "+to_string_with_precision(ch2,4));
      
      
      
    }
    
    if( Tag_reco_type == "Beta_3.99_Emax_4.0" ) {
      double Sew= 1.0201;
      double Vud2= 1.0; // 0.97373;
      Vud2 *= Vud2;
      cout<<"L  TA    TV   T   A   tot"<<endl;
      cout<<" ~ $"<<to_string_with_precision(Vud2*Sew*Final_extr_data[0].ave(),3)<<"~("<<to_string_with_precision(1000*Vud2*Sew*Final_extr_data[0].err(),0)<<")$ ~ &";
      cout<<" ~ $"<<to_string_with_precision(Vud2*Sew*Final_extr_data[1].ave(),3)<<"~("<<to_string_with_precision(1000*Vud2*Sew*Final_extr_data[1].err(),0)<<")$ ~ &";
      cout<<" ~ $"<<to_string_with_precision(Vud2*Sew*Final_extr_data[2].ave(),3)<<"~("<<to_string_with_precision(1000*Vud2*Sew*Final_extr_data[2].err(),0)<<")$ ~ &";
      cout<<" ~ $"<<to_string_with_precision(Vud2*Sew*Final_extr_data[7].ave(),3)<<"~("<<to_string_with_precision(1000*Vud2*Sew*Final_extr_data[7].err(),0)<<")$ ~ &";
      cout<<" ~ $"<<to_string_with_precision(Vud2*Sew*Final_extr_data[6].ave(),3)<<"~("<<to_string_with_precision(1000*Vud2*Sew*Final_extr_data[6].err(),0)<<")$ ~ &";
      cout<<" ~ $"<<to_string_with_precision(Vud2*Sew*Final_extr_data[3].ave(),3)<<"~("<<to_string_with_precision(1000*Vud2*Sew*Final_extr_data[3].err(),0)<<")$ ~ \\"<<endl;


      cout<<"Ratio V/A "<<endl;
      cout<<( (Final_extr_data[2]/(Final_extr_data[0]+Final_extr_data[1]))).ave()<<" +- "<<( (Final_extr_data[2]/(Final_extr_data[0]+Final_extr_data[1]))).err()<<endl;
      cout<<"Direct calculation: "<<Final_extr_data[4].ave()<<" +- "<<Final_extr_data[4].err()<<endl;
      cout<<"Ratio (V-A)/(V+A)"<<endl;
      cout<<"Direct calculation: "<<Final_extr_data[5].ave()<<" +- "<<Final_extr_data[5].err()<<endl;
      

      
    }

  }


	
      

      



  





  
  cout<<"Bye"<<endl;


  
    
 
  


  return;

}
