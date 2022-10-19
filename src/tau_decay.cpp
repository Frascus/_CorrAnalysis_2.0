#include "../include/tau_decay.h"


const double MPiPhys=0.135;
const double alpha = 1.0/137.04;
const double e2 = alpha*4.0*M_PI;
const bool UseJack=1;
const int Njacks=50;
const int Nboots=200;
const double qu= 2.0/3.0;
const double qd= -1.0/3.0;
const double qs= qd;
const double qc= qu;
const double fm_to_inv_Gev= 1.0/0.197327;
const int prec = 256;
const double Nc=3;
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
const double E0_l = 0.05;
Vfloat sigma_list;
const double C_V = 2*M_PI/(pow(m_tau,3));
const double GAMMA_FACT= 12*M_PI; //12*M_PI*pow(Vud*GF,2);
const string MODE="TANT";
bool Use_t_up_to_T_half=true;
const int sm_func_mode= 0;
const string SM_TYPE_0= "KL_"+to_string(sm_func_mode);
const string SM_TYPE_1= "KT_"+to_string(sm_func_mode);
VVfloat covariance_fake;
int Num_LUSCH=17;
int Nres= 15;
int pts_spline=200;
bool COMPUTE_SPEC_DENS_FREE=false;
using namespace std;



void get_sigma_list() {

  bool test=true;
  if(test) { sigma_list.push_back(0.05); return;}

  double s_max= 0.2;
  double s_min= 0.004420;

  sigma_list.push_back(s_max);
  double s=s_max;
  while (s>= s_min +1.0e-6) { s /= pow(2,0.25); sigma_list.push_back(s);}

  return;

}

void Get_spec_dens_free() {

  Vfloat ams({0.00072, 0.00060, 0.00054});

  for( auto &am: ams)  {
    Compute_free_spectral_density(3, am, 1, 0.01, "tau_decay");
    Compute_free_spectral_density(3, am, -1, 0.01, "tau_decay");
    cout<<"am: "<<am<<" computed!"<<endl;
  }
  
  exit(-1);
  return;
}


void Compute_tau_decay_width() {

  PrecFloat::setDefaultPrecision(prec);

  //generate Luscher's energy level

  if(COMPUTE_SPEC_DENS_FREE) Get_spec_dens_free();

  get_sigma_list();


  //Init LL_functions;
  //find first  zeros of the Lusher functions
  Vfloat Luscher_zeroes;
  Zeta_function_zeroes(Num_LUSCH, Luscher_zeroes);

  //############################################INTERPOLATE PHI FUNCTION AND DERIVATIVES#############################

  VVfloat phi_data, phi_der_data;
  Vfloat sx_int;
  Vfloat sx_der, dx_der;
  Vfloat Dz;

  for(int L_zero=0;L_zero<Nres+1;L_zero++) {
    cout<<"Computing n(Lusch): "<<L_zero<<endl;
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



 
   

  LL_functions LL(phi_data,phi_der_data,sx_der, dx_der, sx_int, Dz, Nres, Luscher_zeroes);
    
  //###########################################END INTERPOLATION PHI FUNCTION AND DERIVATIVES################################
  cout<<"####Spline for phi(x) and phi'(x) successfully generated!"<<endl;

  //create directories

  boost::filesystem::create_directory("../data/tau_decay");
  boost::filesystem::create_directory("../data/tau_decay/light");
  boost::filesystem::create_directory("../data/tau_decay/light/Br");
  boost::filesystem::create_directory("../data/tau_decay/light/corr");
  boost::filesystem::create_directory("../data/tau_decay/light/covariance");




 




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



  data_t Vk_data_tm, V0_data_tm, Ak_data_tm, A0_data_tm;
  data_t Vk_data_OS, V0_data_OS, Ak_data_OS, A0_data_OS;


  //light 
  //tm
  Vk_data_tm.Read("../tau_decay_data/light", "mes_contr_2pts_ll_1", "VKVK", Sort_light_confs);
  V0_data_tm.Read("../tau_decay_data/light", "mes_contr_2pts_ll_1", "V0V0", Sort_light_confs);
  Ak_data_tm.Read("../tau_decay_data/light", "mes_contr_2pts_ll_1", "AKAK", Sort_light_confs);
  A0_data_tm.Read("../tau_decay_data/light", "mes_contr_2pts_ll_1", "A0A0", Sort_light_confs);
  //OS
  Vk_data_OS.Read("../tau_decay_data/light", "mes_contr_2pts_ll_2", "VKVK", Sort_light_confs);
  V0_data_OS.Read("../tau_decay_data/light", "mes_contr_2pts_ll_2", "V0V0", Sort_light_confs);
  Ak_data_OS.Read("../tau_decay_data/light", "mes_contr_2pts_ll_2", "AKAK", Sort_light_confs);
  A0_data_OS.Read("../tau_decay_data/light", "mes_contr_2pts_ll_2", "A0A0", Sort_light_confs);


 








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


  vector<distr_t_list> Br_tau_tm, Br_tau_OS;
  vector<distr_t_list> Br_A0_tau_tm, Br_Aii_tau_tm, Br_Vii_tau_tm;
  vector<distr_t_list> Br_A0_tau_OS, Br_Aii_tau_OS, Br_Vii_tau_OS;

  for(int iens=0; iens<Nens; iens++) {
    Br_tau_tm.emplace_back( UseJack, sigma_list.size());
    Br_A0_tau_tm.emplace_back( UseJack, sigma_list.size());
    Br_Aii_tau_tm.emplace_back( UseJack, sigma_list.size());
    Br_Vii_tau_tm.emplace_back( UseJack, sigma_list.size());
    
    Br_tau_OS.emplace_back( UseJack, sigma_list.size());
    Br_A0_tau_OS.emplace_back( UseJack, sigma_list.size());
    Br_Aii_tau_OS.emplace_back( UseJack, sigma_list.size());
    Br_Vii_tau_OS.emplace_back( UseJack, sigma_list.size());

  }

  //resize vector with systematic errors
  VVfloat syst_per_ens_tm_A0(Nens);
  VVfloat syst_per_ens_tm_Ak(Nens);
  VVfloat syst_per_ens_tm_Vk(Nens);
  VVfloat syst_per_ens_OS_A0(Nens);
  VVfloat syst_per_ens_OS_Ak(Nens);
  VVfloat syst_per_ens_OS_Vk(Nens);
  VVfloat syst_per_ens_tm(Nens);
  VVfloat syst_per_ens_OS(Nens);
  
 



  //loop over the ensembles
  #pragma omp parallel for
  for(int iens=0; iens<Nens;iens++) {


    //print number of gauge configurations
    cout<<"Nconfs : "<<Vk_data_OS.Nconfs[iens]<<endl;

    LatticeInfo L_info;
    L_info.LatInfo_new_ens(Vk_data_tm.Tag[iens]);



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
    CorrAnalysis Corr_block_1(UseJack, Vk_data_tm.Nconfs[iens],Nboots);
    Corr_block_1.Nt= Vk_data_tm.nrows[iens];
    Corr.Nt = Vk_data_tm.nrows[iens];
    int T = Corr.Nt;
    

    cout<<"Analyzing Ensemble: "<<Vk_data_tm.Tag[iens]<<endl;
    cout<<"NT: "<<T<<endl;
     
    //get lattice spacing
    distr_t a_distr(UseJack);
    distr_t Zv(UseJack), Za(UseJack);
    double Mpi=0.0;
    double Mpi_err=0.0;
    double fpi=0.0;
    if(Vk_data_tm.Tag[iens].substr(1,1)=="A") {a_distr=a_A; Zv = ZV_A; Za = ZA_A;}
    else if(Vk_data_tm.Tag[iens].substr(1,1)=="B") {a_distr=a_B; Zv = ZV_B; Za = ZA_B; Mpi=0.05653312833; Mpi_err=1.430196186e-05; fpi=0.05278353769;}
    else if(Vk_data_tm.Tag[iens].substr(1,1)=="C") {a_distr=a_C; Zv = ZV_C; Za = ZA_C; Mpi=0.04722061628; Mpi_err=3.492993579e-05; fpi=0.0450246;}
    else if(Vk_data_tm.Tag[iens].substr(1,1)=="D") {a_distr=a_D; Zv = ZV_D; Za = ZA_D; Mpi=0.04062107883; Mpi_err= 2.973916243e-05; fpi=0.03766423429;}
    else crash("lattice spacing distribution for Ens: "+Vk_data_tm.Tag[iens]+" not found");

    //jack distr for Mpi
    distr_t Mpi_distr(UseJack);
    for(int ij=0;ij<Njacks;ij++) Mpi_distr.distr.push_back( Mpi + GM()*Mpi_err/sqrt(Njacks-1.0));

    distr_t resc_GeV = C_V*GAMMA_FACT/(a_distr*a_distr*a_distr);

    //tm
    distr_t_list Vk_tm_distr, V0_tm_distr, Ak_tm_distr, A0_tm_distr;
    //tm block1
    distr_t_list Vk_tm_block_1_distr, V0_tm_block_1_distr, Ak_tm_block_1_distr, A0_tm_block_1_distr;
    //OS
    distr_t_list Vk_OS_distr, V0_OS_distr, Ak_OS_distr, A0_OS_distr;
    //OS block1
    distr_t_list Vk_OS_block_1_distr, V0_OS_block_1_distr, Ak_OS_block_1_distr, A0_OS_block_1_distr;
   
    //light-tm sector
    Vk_tm_distr = Corr.corr_t(Vk_data_tm.col(0)[iens], "../data/tau_decay/light/corr/Vk_tm_"+Vk_data_tm.Tag[iens]+".dat");
    V0_tm_distr = Corr.corr_t(V0_data_tm.col(0)[iens], "../data/tau_decay/light/corr/V0_tm_"+Vk_data_tm.Tag[iens]+".dat");
    Ak_tm_distr = Corr.corr_t(Ak_data_tm.col(0)[iens], "../data/tau_decay/light/corr/Ak_tm_"+Vk_data_tm.Tag[iens]+".dat");
    A0_tm_distr = Corr.corr_t(A0_data_tm.col(0)[iens], "../data/tau_decay/light/corr/A0_tm_"+Vk_data_tm.Tag[iens]+".dat");

    //light-OS sector
    Vk_OS_distr = Corr.corr_t(Vk_data_OS.col(0)[iens], "../data/tau_decay/light/corr/Vk_OS_"+Vk_data_tm.Tag[iens]+".dat");
    V0_OS_distr = Corr.corr_t(V0_data_OS.col(0)[iens], "../data/tau_decay/light/corr/V0_OS_"+Vk_data_tm.Tag[iens]+".dat");
    Ak_OS_distr = Corr.corr_t(Ak_data_OS.col(0)[iens], "../data/tau_decay/light/corr/Ak_OS_"+Vk_data_tm.Tag[iens]+".dat");
    A0_OS_distr = Corr.corr_t(A0_data_OS.col(0)[iens], "../data/tau_decay/light/corr/A0_OS_"+Vk_data_tm.Tag[iens]+".dat");


    //analyze data with Njacks=Nconfs
    //light-tm sector
    Vk_tm_block_1_distr = Corr_block_1.corr_t(Vk_data_tm.col(0)[iens], "");
    V0_tm_block_1_distr = Corr_block_1.corr_t(V0_data_tm.col(0)[iens], "");
    Ak_tm_block_1_distr = Corr_block_1.corr_t(Ak_data_tm.col(0)[iens], "");
    A0_tm_block_1_distr = Corr_block_1.corr_t(A0_data_tm.col(0)[iens], "");

    //light-OS sector
    Vk_OS_block_1_distr = Corr_block_1.corr_t(Vk_data_OS.col(0)[iens], "");
    V0_OS_block_1_distr = Corr_block_1.corr_t(V0_data_OS.col(0)[iens], "");
    Ak_OS_block_1_distr = Corr_block_1.corr_t(Ak_data_OS.col(0)[iens], "");
    A0_OS_block_1_distr = Corr_block_1.corr_t(A0_data_OS.col(0)[iens], "");

  
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

    Print_To_File({}, {TT,RR,cov_A0_tm, corr_m_A0_tm}, "../data/tau_decay/light/covariance/A0_tm_"+Vk_data_tm.Tag[iens]+".dat", "", "");
    Print_To_File({}, {TT,RR,cov_Ak_tm, corr_m_Ak_tm}, "../data/tau_decay/light/covariance/Ak_tm_"+Vk_data_tm.Tag[iens]+".dat", "", "");
    Print_To_File({}, {TT,RR,cov_Vk_tm, corr_m_Vk_tm}, "../data/tau_decay/light/covariance/Vk_tm_"+Vk_data_tm.Tag[iens]+".dat", "", "");
    Print_To_File({}, {TT,RR,cov_A0_OS, corr_m_A0_OS}, "../data/tau_decay/light/covariance/A0_OS_"+Vk_data_tm.Tag[iens]+".dat", "", "");
    Print_To_File({}, {TT,RR,cov_Ak_OS, corr_m_Ak_OS}, "../data/tau_decay/light/covariance/Ak_OS_"+Vk_data_tm.Tag[iens]+".dat", "", "");
    Print_To_File({}, {TT,RR,cov_Vk_OS, corr_m_Vk_OS}, "../data/tau_decay/light/covariance/Vk_OS_"+Vk_data_tm.Tag[iens]+".dat", "", "");

       

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

    bool Found_error_less_x_percent=false;
    double x=15;
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
    auto f_syst_A0_OS = [&resc_GeV, &Mpi, &fpi](const function<double(double)> &F) -> double { return resc_GeV.ave()*F(Mpi)*pow(fpi,2)*Mpi/(2.0);};
      
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
		      gsl_integration_qags(G_SYST_asympt, E0_l*a_distr.ave(), 4.0,  0.0, 1e-3, 1000, w_SYST_asympt, &val_mod, &err_mod);
		      if(err_mod/fabs(val_mod) > 1e-2) crash("Cannot reach accuracy in evaluating systematic");
		      gsl_integration_workspace_free(w_SYST_asympt);
		      systs.push_back( fabs(val_mod));
					   
	      
		      //evaluate GS in infinite volume limit
		      auto FS_infL = [&F, &Mpi, &a_distr, &resc_GeV, &LL, &Mrhos, &gppis, &Eduals, &Mpi, &Rduals, &F_free_tm](double E) {
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
		      gsl_integration_qags(G_SYST_infL, E0_l*a_distr.ave(), 4.0, 0.0, 1e-3, 1000, w_SYST_infL, &val_mod_infL, &err_mod_infL);
		      gsl_integration_workspace_free(w_SYST_infL);
		      systs.push_back(fabs(val_mod_infL));

			
		      double syst_ave=0.0;
		      double syst_dev=0.0;
		      for (int n=0;n<(signed)systs.size();n++) { syst_ave += systs[n]/systs.size(); syst_dev += systs[n]*systs[n];}

		      printV(systs, "SYST", 0);
		      syst_dev =  sqrt( (1.0/(systs.size()-1.0))*(syst_dev - systs.size()*syst_ave*syst_ave) );
		      double max = *max_element(systs.begin(), systs.end());
		      cout<<"systematic computed!!: average: "<<syst_ave<<" max: "<<max<<endl;
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
		      gsl_integration_qags(G_SYST_asympt, E0_l*a_distr.ave(), 4.0,  0.0, 1e-3, 1000, w_SYST_asympt, &val_mod, &err_mod);
		      if(err_mod/fabs(val_mod) > 1e-2) crash("Cannot reach accuracy in evaluating systematic");
		      gsl_integration_workspace_free(w_SYST_asympt);
		      systs.push_back( fabs(val_mod));
					   
	      
		      //evaluate GS in infinite volume limit
		      auto FS_infL = [&F, &Mpi, &a_distr, &resc_GeV, &LL, &Mrhos, &gppis, &Eduals, &Mpi, &Rduals, &F_free_OS](double E) {
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
		      gsl_integration_qags(G_SYST_infL, E0_l*a_distr.ave(), 4.0, 0.0, 1e-3, 1000, w_SYST_infL, &val_mod_infL, &err_mod_infL);
		      gsl_integration_workspace_free(w_SYST_infL);
		      systs.push_back(fabs(val_mod_infL));

			
		      double syst_ave=0.0;
		      double syst_dev=0.0;
		      for (int n=0;n<(signed)systs.size();n++) { syst_ave += systs[n]/systs.size(); syst_dev += systs[n]*systs[n];}

		      printV(systs, "SYST", 0);
		      syst_dev =  sqrt( (1.0/(systs.size()-1.0))*(syst_dev - systs.size()*syst_ave*syst_ave) );
		      double max = *max_element(systs.begin(), systs.end());
		      cout<<"systematic computed!!: average: "<<syst_ave<<" max: "<<max<<endl;
		      return max;
		   
		    };




    // lambda function to be used as a smearing func.
  
    const auto K0 = [&a_distr](const PrecFloat &E, const PrecFloat &m, const PrecFloat &s, const PrecFloat &E0, int ijack) -> PrecFloat {


		      
		      PrecFloat X;
		      PrecFloat X_ave = E/(m_tau*a_distr.ave());
		      if(X_ave < E0) return 0.0;

		      
		      
		      if(ijack==-1) {
			X=X_ave;
		      }
		      else X= E/(m_tau*a_distr.distr[ijack]);
		     
		      PrecFloat sm_theta;

		      if(sm_func_mode==0) sm_theta= 1/(1+ exp(-(1-X)/s));
		      else if(sm_func_mode==1) sm_theta= 1/(1+ exp(-sinh((1-X)/s)));
		      else if(sm_func_mode==2) sm_theta= (1+erf((1-X)/s))/2;
		      else crash("sm_func_mode: "+to_string(sm_func_mode)+" not yet implemented");
						 
		      return (1/X)*pow(( 1 -pow(X,2)),2)*sm_theta;
		   
		 };

    const auto K1 = [&a_distr](const PrecFloat &E, const PrecFloat &m, const PrecFloat &s, const PrecFloat &E0, int ijack) -> PrecFloat {

		      PrecFloat X;
		      PrecFloat X_ave = E/(m_tau*a_distr.ave());
		      if( X_ave < E0) return 0.0;


		      if(ijack==-1) {
			X=X_ave;
		      }
		      else X= E/(m_tau*a_distr.distr[ijack]);

		      
		      PrecFloat sm_theta;
		      if(sm_func_mode==0) sm_theta= 1/(1+ exp(-(1-X)/s));
		      else if(sm_func_mode==1) sm_theta= 1/(1+ exp(-sinh((1-X)/s)));
		      else if(sm_func_mode==2) sm_theta= (1+erf((1-X)/s))/2;
		      else crash("sm_func_mode: "+to_string(sm_func_mode)+" not yet implemented");
		      
		      return (1 + 2*pow(X,2))*(1/(X))*pow(( 1 -pow(X,2)),2)*sm_theta;
		   
		    };


    if(Use_t_up_to_T_half) {
      tmax_tm_1_Aii = tmax_tm_1_Vii= tmax_tm_0 = tmax_OS_0 = tmax_OS_1_Aii = tmax_OS_1_Vii= Corr.Nt/2 -1;
    }

    cout<<"tmax_tm_0: "<<tmax_tm_0<<" tmax_tm_1_Aii: "<<tmax_tm_1_Aii<<" tmax_tm_1_Vii: "<<tmax_tm_1_Vii<<" tmax_OS_0: "<<tmax_OS_0<<" tmax_OS_1_Aii: "<<tmax_OS_1_Aii<<" tmax_OS_1_Vii: "<<tmax_OS_1_Vii<<endl;
  
    cout<<"sigma : {";
    for(int is=0;is<(signed)sigma_list.size();is++) { cout<<sigma_list[is]<<", ";}
    cout<<"}"<<endl;

    //loop over sigma
    #pragma omp parallel for
    for(int is=0; is < (signed)sigma_list.size(); is++) {

      distr_t Br_sigma_A0_tm;
      distr_t Br_sigma_Aii_tm;
      distr_t Br_sigma_Vii_tm;
      distr_t Br_sigma_A0_OS;
      distr_t Br_sigma_Aii_OS;
      distr_t Br_sigma_Vii_OS;
      double s= sigma_list[is];
      //int tmax= T/2 -4;
      double lA0_tm, lAii_tm, lVii_tm;
      double lA0_OS, lAii_OS, lVii_OS;
      double syst_A0_tm, syst_Aii_tm, syst_Vii_tm;
      double syst_A0_OS, syst_Aii_OS, syst_Vii_OS;
      double mult=1e4;
      if(MODE=="SANF") mult= 1e3;
      const auto model_estimate_V_tm = [&K1, &GS_V, &s, &a_distr, &Ergs, &Amplitudes, &Mrhos, &Eduals, &Rduals ](double E) -> double {  return K1(E,0.0,s,E0_l*a_distr.ave(), -1).get()*GS_V(E, Ergs[0], Amplitudes[0], Mrhos[0], Eduals[0], Rduals[0]);};
      const auto model_V_tm =  [&GS_V, &Ergs, &Amplitudes, &a_distr, &Mrhos, &Eduals, &Rduals](double E) -> double {  return GS_V(E, Ergs[0], Amplitudes[0], Mrhos[0], Eduals[0], Rduals[0]);};
      const auto model_estimate_V_OS = [&K1, &GS_V, &s, &a_distr, &Ergs, &Amplitudes, &Mrhos, &Eduals, &Rduals ](double E) -> double {  return K1(E,0.0,s,E0_l*a_distr.ave(), -1).get()*GS_V(E, Ergs[1], Amplitudes[1], Mrhos[1], Eduals[1], Rduals[1]);};
      const auto model_V_OS =  [&GS_V, &Ergs, &Amplitudes, &a_distr, &Mrhos, &Eduals, &Rduals](double E) -> double {  return GS_V(E, Ergs[1], Amplitudes[1], Mrhos[1], Eduals[1], Rduals[1]);};
      const auto model_A0 = [&Mpi, &fpi, &resc_GeV ](double E) -> double { double s=Mpi*0.001; return resc_GeV.ave()*(fpi*fpi*Mpi/(2.0))*(1.0/sqrt( 2*M_PI*s*s))*exp( -0.5*pow((E-Mpi)/s,2));};
      //compute model estimate of vector contribution to R_tau tm
      gsl_function_pp<decltype(model_estimate_V_tm)> MOD_tm(model_estimate_V_tm);
      gsl_integration_workspace * w_MOD_tm = gsl_integration_workspace_alloc (10000);
      gsl_function *G_MOD_tm = static_cast<gsl_function*>(&MOD_tm);
      double val_mod_tm,err_mod_tm;
      gsl_integration_qagiu(G_MOD_tm, 0.0, 0.0, 1e-4, 10000, w_MOD_tm, &val_mod_tm, &err_mod_tm);
      gsl_integration_workspace_free(w_MOD_tm);
      cout<<"Model estimate R_t(V,tm): "<<val_mod_tm<<" +- "<<err_mod_tm<<endl;
      //compute model estimate of vector contribution to R_tau OS
      gsl_function_pp<decltype(model_estimate_V_OS)> MOD_OS(model_estimate_V_OS);
      gsl_integration_workspace * w_MOD_OS = gsl_integration_workspace_alloc (10000);
      gsl_function *G_MOD_OS = static_cast<gsl_function*>(&MOD_OS);
      double val_mod_OS,err_mod_OS;
      gsl_integration_qagiu(G_MOD_OS, 0.0, 0.0, 1e-4, 10000, w_MOD_OS, &val_mod_OS, &err_mod_OS);
      gsl_integration_workspace_free(w_MOD_OS);
      cout<<"Model estimate R_t(V,OS): "<<val_mod_OS<<" +- "<<err_mod_OS<<endl;
      cout<<"Model estimate R_t(A0, tm&OS): "<<K0(Mpi, 0.0, s ,E0_l*a_distr.ave(),-1).get()*resc_GeV.ave()*pow(fpi,2)*Mpi/2.0<<endl;
      
      

      Br_sigma_Vii_tm = Get_Laplace_transfo(  0.0,  s, E0_l*a_distr.ave(),  T, tmax_tm_1_Vii, prec, SM_TYPE_1,K1, Vii_tm, syst_Vii_tm, 1e4, lVii_tm, MODE, "tm", "Vii_light_"+Vk_data_tm.Tag[iens], -1,0, resc_GeV*Za*Za, "tau_decay", cov_Vk_tm, f_syst_V_tm,1, model_V_tm );
      cout<<"BrVii_tm["<<Vk_data_tm.Tag[iens]<<"] , s= "<<s<<", l= "<<lVii_tm<<" : "<<Br_sigma_Vii_tm.ave()<<" +- "<<Br_sigma_Vii_tm.err()<<endl;
      Br_sigma_Vii_OS = Get_Laplace_transfo(  0.0,  s, E0_l*a_distr.ave(),  T, tmax_OS_1_Vii, prec, SM_TYPE_1,K1, Vii_OS, syst_Vii_OS, 1e4, lVii_OS, MODE, "OS", "Vii_light_"+Vk_data_tm.Tag[iens],-1,0, resc_GeV*Zv*Zv, "tau_decay", cov_Vk_OS, f_syst_V_OS,1, model_V_OS );
      cout<<"BrVii_OS["<<Vk_data_tm.Tag[iens]<<"] , s= "<<s<<", l= "<<lVii_OS<<" : "<<Br_sigma_Vii_OS.ave()<<" +- "<<Br_sigma_Vii_OS.err()<<endl;
      Br_sigma_A0_tm = Get_Laplace_transfo(  0.0,  s, E0_l*a_distr.ave(),  T, tmax_tm_0, prec, SM_TYPE_0,K0, -1*A0_tm, syst_A0_tm, 1e4, lA0_tm, MODE, "tm", "A0_light_"+Vk_data_tm.Tag[iens], -1, 0, resc_GeV*Zv*Zv, "tau_decay", cov_A0_tm, f_syst_A0_tm,1, model_A0 );
      cout<<"BrA0_tm["<<Vk_data_tm.Tag[iens]<<"] , s= "<<s<<", l= "<<lA0_tm<<" : "<<Br_sigma_A0_tm.ave()<<" +- "<<Br_sigma_A0_tm.err()<<endl; 
      Br_sigma_A0_OS = Get_Laplace_transfo(  0.0,  s, E0_l*a_distr.ave(),  T, tmax_OS_0, prec, SM_TYPE_0,K0, -1*A0_OS, syst_A0_OS, 1e4, lA0_OS, MODE, "OS", "A0_light_"+Vk_data_tm.Tag[iens], -1, 0, resc_GeV*Za*Za, "tau_decay", cov_A0_OS, f_syst_A0_OS,1, model_A0 );
      cout<<"BrA0_OS["<<Vk_data_tm.Tag[iens]<<"] , s= "<<s<<", l= "<<lA0_OS<<" : "<<Br_sigma_A0_OS.ave()<<" +- "<<Br_sigma_A0_OS.err()<<endl; 
      Br_sigma_Aii_tm = Get_Laplace_transfo(  0.0,  s, E0_l*a_distr.ave(),  T, tmax_tm_1_Aii, prec, SM_TYPE_1,K1, Aii_tm, syst_Aii_tm, 1e4, lAii_tm, MODE, "tm", "Aii_light_"+Vk_data_tm.Tag[iens], -1,0, resc_GeV*Zv*Zv, "tau_decay", cov_Ak_tm, fake_func,0, fake_func_d );
      cout<<"BrAii_tm["<<Vk_data_tm.Tag[iens]<<"] , s= "<<s<<", l= "<<lAii_tm<<" : "<<Br_sigma_Aii_tm.ave()<<" +- "<<Br_sigma_Aii_tm.err()<<endl;
      Br_sigma_Aii_OS = Get_Laplace_transfo(  0.0,  s, E0_l*a_distr.ave(),  T, tmax_OS_1_Aii, prec, SM_TYPE_1,K1, Aii_OS, syst_Aii_OS, 1e4, lAii_OS, MODE, "OS", "Aii_light_"+Vk_data_tm.Tag[iens], -1,0,resc_GeV*Za*Za, "tau_decay", cov_Ak_OS, fake_func,0, fake_func_d );
      cout<<"BrAii_OS["<<Vk_data_tm.Tag[iens]<<"] , s= "<<s<<", l= "<<lAii_OS<<" : "<<Br_sigma_Aii_OS.ave()<<" +- "<<Br_sigma_Aii_OS.err()<<endl;
      


      //push_back systematic error
     

      
      syst_per_ens_tm_A0[iens].push_back( syst_A0_tm);
      syst_per_ens_tm_Ak[iens].push_back( syst_Aii_tm);
      syst_per_ens_tm_Vk[iens].push_back( syst_Vii_tm);
      syst_per_ens_tm[iens].push_back( sqrt( pow(syst_A0_tm,2)+ pow(syst_Aii_tm,2)+ pow(syst_Vii_tm,2)));
      syst_per_ens_OS_A0[iens].push_back( syst_A0_OS);
      syst_per_ens_OS_Ak[iens].push_back( syst_Aii_OS);
      syst_per_ens_OS_Vk[iens].push_back( syst_Vii_OS);
      syst_per_ens_OS[iens].push_back( sqrt( pow(syst_A0_OS,2)+ pow(syst_Aii_OS,2)+ pow(syst_Vii_OS,2)));

      

      
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
   

    }
    
    
  }
  
    
   




  //Print to File
  for(int iens=0; iens<Nens;iens++) {
    Print_To_File({}, {sigma_list, Br_tau_tm[iens].ave(), Br_tau_tm[iens].err(), syst_per_ens_tm[iens] , Br_tau_OS[iens].ave(), Br_tau_OS[iens].err(), syst_per_ens_OS[iens]}, "../data/tau_decay/light/Br/br_"+MODE+"_"+Vk_data_tm.Tag[iens]+".dat", "", "#sigma Br[tm] Br[OS]");
    Print_To_File({}, {sigma_list, Br_A0_tau_tm[iens].ave(), Br_A0_tau_tm[iens].err(), syst_per_ens_tm_A0[iens], Br_Aii_tau_tm[iens].ave(), Br_Aii_tau_tm[iens].err(), syst_per_ens_tm_Ak[iens], Br_Vii_tau_tm[iens].ave(), Br_Vii_tau_tm[iens].err(), syst_per_ens_tm_Vk[iens]}, "../data/tau_decay/light/Br/br_contrib_tm_"+MODE+"_"+Vk_data_tm.Tag[iens]+".dat", "", "#sigma Br_A0 Br_Aii Br_Vii");
    Print_To_File({}, {sigma_list, Br_A0_tau_OS[iens].ave(), Br_A0_tau_OS[iens].err(), syst_per_ens_OS_A0[iens], Br_Aii_tau_OS[iens].ave(), Br_Aii_tau_OS[iens].err(), syst_per_ens_OS_Ak[iens], Br_Vii_tau_OS[iens].ave(), Br_Vii_tau_OS[iens].err(), syst_per_ens_OS_Vk[iens]}, "../data/tau_decay/light/Br/br_contrib_OS_"+MODE+"_"+Vk_data_tm.Tag[iens]+".dat", "", "#sigma Br_A0 Br_Aii Br_Vii");
  }

  cout<<"Output per ensemble printed..."<<endl;

  //Print all ens for each sigma
  for(int is=0; is<(signed)sigma_list.size();is++) {

    Vfloat syst_tm, syst_OS;
    distr_t_list Br_tau_tm_fixed_sigma(UseJack), Br_tau_OS_fixed_sigma(UseJack);
    for(int iens=0;iens<Nens;iens++) {
      Br_tau_tm_fixed_sigma.distr_list.push_back( Br_tau_tm[iens].distr_list[is]);
      Br_tau_OS_fixed_sigma.distr_list.push_back( Br_tau_OS[iens].distr_list[is]);
      syst_tm.push_back( syst_per_ens_tm[iens][is]);
      syst_OS.push_back( syst_per_ens_OS[iens][is]);
    }

    Print_To_File(Vk_data_tm.Tag,{ Br_tau_tm_fixed_sigma.ave(), Br_tau_tm_fixed_sigma.err(), syst_tm, Br_tau_OS_fixed_sigma.ave(), Br_tau_OS_fixed_sigma.err(), syst_OS},"../data/tau_decay/light/Br/br_sigma_"+to_string_with_precision(sigma_list[is],3)+".dat", "", "#Ens tm  OS");
    
  }


  
  
    
 
  


  return;

}
