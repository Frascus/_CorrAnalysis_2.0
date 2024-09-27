#include "../include/virtual_07.h"
#include "Corr_analysis.h"
#include "Spectral.h"
#include "numerics.h"
#include "stat.h"
using namespace std;


bool verbose_lev_07=1;
//Vfloat sigmas_07({1.75, 2, 2.25, 2.5, 2.75, 3.0, 3.5});
Vfloat sigmas_07({1.75, 2, 2.25, 2.5, 2.75, 3.0, 3.5, 4.0, 4.5});  
//Vfloat sigmas_07({ 2, 3.5});
Vfloat sigmas_07_w0({0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29,  0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.25, 1.5, 1.75, 2.0,2.25,2.5,2.75,3.0});
int prec_07=128;
const string MODE_FF="TANT";
const bool Skip_spectral_reconstruction_07 = true;
const bool virtuality_scan = false;
const bool Use_preconditioning = true;
const string  preco_tag= (Use_preconditioning)?"prec_":"";
const double Mjpsi= 3.0969; //GeV
const double Mphi= 1.019461; //GeV
const double MDs_phys = 1.96847; // GeV
const double MBs = 5.36692;
const double E0_fact = 0.90;
const bool CONS_EM_CURRENT = true;

const double QU = -1.0/3.0;
const double QD = -1.0 / 3.0;
const double mh0 = 1.0 / 0.5074431945705498;
const double mh1 = 1.0 / 0.4007281135883006;
const double mh3 = 1.0 / 0.2845905478176752;
const Vfloat masses({mh0, mh1, mh3});


rt_07_Bs Get_virtual_tensor_FF(int n_xg, bool UseJack, int Njacks, string MESON,  string Corr_path,string path_out) {

  Njacks=25;

  Vfloat virtualities;
  for(int i=0;i<75;i++) { virtualities.push_back( 3*i/74.0) ; }
  
  rt_07_Bs return_class;

  //old 65, 60, 57, 55
  //new 67  62  55  50
  Vfloat y_eff({0.65, 0.60, 0.57, 0.55});
  return_class.y_eff= y_eff;



  PrecFloat::setDefaultPrecision(prec_07);
  cout<<"max possible exponent: "<<PrecFloat::getEmax_max()<<endl;
  cout<<"current max exponent: "<<PrecFloat::getEmax()<<endl;
  cout<<"Number of xg to analyze: "<<n_xg<<endl;

  int t_07_s, t_07_s_HLT, t_07_c;
 
  string TAG_CURR="";
  if(CONS_EM_CURRENT==false) TAG_CURR="LOC_";

 
  int size_mu_nu= 4;

  double sign_kz=-1.0; //correct one is -1

  //BK
  vector<vector<vector<data_t>>> C_B_d_data(size_mu_nu);
  //T
  vector<vector<vector<data_t>>>  C_T_d_data(size_mu_nu);


  vector<vector<vector<data_t>>> C_B_u_data_std(size_mu_nu), C_B_d_data_std(size_mu_nu);
  vector<vector<vector<data_t>>> C_T_u_data_std(size_mu_nu), C_T_d_data_std(size_mu_nu);



  data_t data_2pts_SM, data_2pts_SMSM;

 
  


  
  for(int mu=0;mu<size_mu_nu;mu++) {

    C_B_d_data[mu].resize(size_mu_nu);
    C_T_d_data[mu].resize(size_mu_nu);

    C_B_u_data_std[mu].resize(size_mu_nu);
    C_B_d_data_std[mu].resize(size_mu_nu);
    C_T_u_data_std[mu].resize(size_mu_nu);
    C_T_d_data_std[mu].resize(size_mu_nu);

   
    for(int nu=0;nu<size_mu_nu;nu++) {

      C_B_d_data[mu][nu].resize(n_xg);
      C_T_d_data[mu][nu].resize(n_xg);

      C_B_u_data_std[mu][nu].resize(n_xg);
      C_B_d_data_std[mu][nu].resize(n_xg);
      C_T_u_data_std[mu][nu].resize(n_xg);
      C_T_d_data_std[mu][nu].resize(n_xg);

    
    }
  }

  //custom sorting of gauge confs
  auto Sort_confs = [](string A, string B) {

			   

    int conf_length_A= A.length();
    int conf_length_B= B.length();

    int pos_a_slash=-1;
    int pos_b_slash=-1;
    for(int i=0;i<conf_length_A;i++) if(A.substr(i,1)=="/") pos_a_slash=i;
    for(int j=0;j<conf_length_B;j++) if(B.substr(j,1)=="/") pos_b_slash=j;

    string A_bis= A.substr(pos_a_slash+1);
    string B_bis= B.substr(pos_b_slash+1);

					     
    string conf_num_A = A_bis.substr(0,4);
    string conf_num_B = B_bis.substr(0,4);
							       
		      
    string rA = A_bis.substr(A_bis.length()-2);
    string rB = B_bis.substr(B_bis.length()-2);
    if(rA.substr(0,1) == "r") { 
      int n1 = stoi(A_bis.substr(A_bis.length()-1));
      int n2 = stoi(B_bis.substr(B_bis.length()-1));
      if(rA == rB) {
	if(rA=="r0" || rA=="r2") return conf_num_A > conf_num_B;
	else if(rA=="r1" || rA=="r3") return conf_num_A < conf_num_B;
	else crash("stream not recognized");
      }
      else return n1<n2;
    }
    return A_bis<B_bis;
  };
  
  //read data

  data_2pts_SM.Read(Corr_path+"/spectre", "mes_contr_2pts_SM_3", "P5P5", Sort_confs);
  data_2pts_SMSM.Read(Corr_path+"/spectre", "mes_contr_2pts_SMSM_3", "P5P5", Sort_confs);

  
 
  //loop over mu and nu axial
  vector<pair<int,int>> mu_nu_pair_B({make_pair(1,1),make_pair(2,2)});
  vector<pair<int,int>> mu_nu_pair_T({make_pair(1,2),make_pair(2,1)});
   
  
   
  for(int ixg=0;ixg<n_xg;ixg++) {


   
    //B
    for(auto &pair_B : mu_nu_pair_B) {
      int mu=pair_B.first;
      int nu=pair_B.second;
  
      string Tag_contr="S0P5";
      if(CONS_EM_CURRENT==false) Tag_contr="V"+to_string(mu)+"P5";
      //B
      //d
      C_B_d_data[mu][nu][ixg].Read(Corr_path+"/spectre", TAG_CURR+"C_d_B_nu_"+to_string(nu)+"_mu_"+to_string(mu)+"_ixg_"+to_string(ixg), Tag_contr, Sort_confs);

      //B
      //u
      C_B_u_data_std[mu][nu][ixg].Read(Corr_path+"/standard", TAG_CURR+"C_u_B_nu_"+to_string(nu)+"_mu_"+to_string(mu)+"_ixg_"+to_string(ixg), Tag_contr, Sort_confs);
      //d
      C_B_d_data_std[mu][nu][ixg].Read(Corr_path+"/standard", TAG_CURR+"C_d_B_nu_"+to_string(nu)+"_mu_"+to_string(mu)+"_ixg_"+to_string(ixg), Tag_contr, Sort_confs);
      
    
    }

    //T
    for(auto &pair_T : mu_nu_pair_T) {
      int mu=pair_T.first;
      int nu=pair_T.second;
    
      string Tag_contr="S0P5";
      if(CONS_EM_CURRENT==false) Tag_contr="V"+to_string(mu)+"P5";
      //T
      //d
      C_T_d_data[mu][nu][ixg].Read(Corr_path+"/spectre", TAG_CURR+"C_d_T_nu_"+to_string(nu)+"_mu_"+to_string(mu)+"_ixg_"+to_string(ixg), Tag_contr, Sort_confs);

      //u
      C_T_u_data_std[mu][nu][ixg].Read(Corr_path+"/standard", TAG_CURR+"C_u_T_nu_"+to_string(nu)+"_mu_"+to_string(mu)+"_ixg_"+to_string(ixg), Tag_contr, Sort_confs);
      //d
      C_T_d_data_std[mu][nu][ixg].Read(Corr_path+"/standard", TAG_CURR+"C_d_T_nu_"+to_string(nu)+"_mu_"+to_string(mu)+"_ixg_"+to_string(ixg), Tag_contr, Sort_confs);
      
    
    }
  }

  
   
  GaussianMersenne GM_07(4455);

  //resample RCs
  distr_t ZT_A(UseJack), ZT_B(UseJack), ZT_C(UseJack), ZT_D(UseJack);
  distr_t ZT_A_boot(0), ZT_B_boot(0), ZT_C_boot(0), ZT_D_boot(0);
  distr_t a_A(UseJack), a_B(UseJack), a_C(UseJack), a_D(UseJack);
  distr_t ZV_A(UseJack), ZV_B(UseJack), ZV_C(UseJack), ZV_D(UseJack);

  
  double fmTGeV= 1.0/0.197327;

  LatticeInfo L_info_A, L_info_B, L_info_C, L_info_D;
  L_info_A.LatInfo_new_ens("cA211a.12.48");
  L_info_B.LatInfo_new_ens("cB211b.072.96");
  L_info_C.LatInfo_new_ens("cC211a.06.80");
  L_info_D.LatInfo_new_ens("cD211a.054.96");

   

  for(int ijack=0; ijack<Njacks;ijack++) {

   

    
    a_A.distr.push_back( L_info_A.a_from_afp*fmTGeV + GM_07()*L_info_A.a_from_afp_err*fmTGeV/((UseJack==true)?sqrt(Njacks -1.0):1.0));
    a_B.distr.push_back( L_info_B.a_from_afp*fmTGeV + GM_07()*L_info_B.a_from_afp_err*fmTGeV/((UseJack==true)?sqrt(Njacks -1.0):1.0));
    a_C.distr.push_back( L_info_C.a_from_afp*fmTGeV + GM_07()*L_info_C.a_from_afp_err*fmTGeV/((UseJack==true)?sqrt(Njacks -1.0):1.0));
    a_D.distr.push_back( L_info_D.a_from_afp*fmTGeV + GM_07()*L_info_D.a_from_afp_err*fmTGeV/((UseJack==true)?sqrt(Njacks -1.0):1.0));

    
    ZT_A.distr.push_back( L_info_A.ZT_RI2 + GM_07()*L_info_A.ZT_RI2_err/((UseJack==true)?sqrt(Njacks -1.0):1.0));
    ZT_B.distr.push_back( L_info_B.ZT_RI2 + GM_07()*L_info_B.ZT_RI2_err/((UseJack==true)?sqrt(Njacks -1.0):1.0));
    ZT_C.distr.push_back( L_info_C.ZT_RI2 + GM_07()*L_info_C.ZT_RI2_err/((UseJack==true)?sqrt(Njacks -1.0):1.0));
    ZT_D.distr.push_back( L_info_D.ZT_RI2 + GM_07()*L_info_D.ZT_RI2_err/((UseJack==true)?sqrt(Njacks -1.0):1.0));

   

       
    ZV_A.distr.push_back( L_info_A.Zv_WI_strange + GM_07()*L_info_A.Zv_WI_strange_err/((UseJack==true)?sqrt(Njacks -1.0):1.0));
    ZV_B.distr.push_back( L_info_B.Zv_WI_strange + GM_07()*L_info_B.Zv_WI_strange_err/((UseJack==true)?sqrt(Njacks -1.0):1.0));
    ZV_C.distr.push_back( L_info_C.Zv_WI_strange + GM_07()*L_info_C.Zv_WI_strange_err/((UseJack==true)?sqrt(Njacks -1.0):1.0));
    ZV_D.distr.push_back( L_info_D.Zv_WI_strange + GM_07()*L_info_D.Zv_WI_strange_err/((UseJack==true)?sqrt(Njacks -1.0):1.0));
    
    
  }

  for(int iboot=0;iboot<1000;iboot++) {

    ZT_A_boot.distr.push_back( L_info_A.ZT_RI2 + GM_07()*L_info_A.ZT_RI2_err);
    ZT_B_boot.distr.push_back( L_info_B.ZT_RI2 + GM_07()*L_info_B.ZT_RI2_err);
    ZT_C_boot.distr.push_back( L_info_C.ZT_RI2 + GM_07()*L_info_C.ZT_RI2_err);
    ZT_D_boot.distr.push_back( L_info_D.ZT_RI2 + GM_07()*L_info_D.ZT_RI2_err);

  }

  cout<<"a_A: "<<a_A.ave()/fmTGeV << " +- "<<a_A.err()/fmTGeV<<endl;
  cout<<"a_B: "<<a_B.ave()/fmTGeV << " +- "<<a_B.err()/fmTGeV<<endl;
  cout<<"a_C: "<<a_C.ave()/fmTGeV << " +- "<<a_C.err()/fmTGeV<<endl;
  cout<<"a_D: "<<a_D.ave()/fmTGeV << " +- "<<a_D.err()/fmTGeV<<endl;

  cout<<"ZT_A : "<<ZT_A.ave()<< " +- "<<ZT_A.err()<<endl;
  cout<<"ZT_B : "<<ZT_B.ave()<< " +- "<<ZT_B.err()<<endl;
  cout<<"ZT_C : "<<ZT_C.ave()<< " +- "<<ZT_C.err()<<endl;
  cout<<"ZT_D : "<<ZT_D.ave()<< " +- "<<ZT_D.err()<<endl;

  cout<<"ZV_A : "<<ZV_A.ave()<< " +- "<<ZV_A.err()<<endl;
  cout<<"ZV_B : "<<ZV_B.ave()<< " +- "<<ZV_B.err()<<endl;
  cout<<"ZV_C : "<<ZV_C.ave()<< " +- "<<ZV_C.err()<<endl;
  cout<<"ZV_D : "<<ZV_D.ave()<< " +- "<<ZV_D.err()<<endl;


  
  
  
 
  int Nens= C_B_d_data[1][1][0].size;
  vector<string> Ens_tags= C_B_d_data[1][1][0].Tag;
  Vint Nts=C_B_d_data[1][1][0].nrows;

  distr_t_list xg_list(UseJack);
  for(int ixg=0;ixg<4;ixg++) xg_list.distr_list.push_back( Get_id_distr(Njacks, UseJack)*(ixg+1)*0.1);
  vector<distr_t_list> F_T_u_list;
  vector<distr_t_list> F_T_d_I_list;
  vector<distr_t_list> F_T_d_I_ori_list;
  vector<distr_t_list> F_T_d_I_ori_sp_list;
  vector<distr_t_list> F_T_d_I_list_15;
  vector<distr_t_list> F_T_d_I_list_7;
  vector<distr_t_list> F_T_d_sp_I_list;
  vector<distr_t_list> FV_T_u_real_list;
  vector<distr_t_list> FV_T_d_real_list;
  vector<distr_t_list> FA_T_u_real_list;
  vector<distr_t_list> FA_T_d_real_list;

  vector<distr_t_list> FA_T_d_I_real_list;
  vector<distr_t_list> FA_T_d_II_real_list;

  vector<distr_t_list> FV_T_d_I_real_list;
  vector<distr_t_list> FV_T_d_II_real_list;

  vector<distr_t_list> FV_T_d_sp_real_list;
  vector<distr_t_list> FA_T_d_sp_real_list;

  vector<distr_t_list> reminder_F_T_d_list;
  vector<distr_t_list> reminder_F_T_d_sp_list;
  

  vector<vector<distr_t_list>> F_T_d_MB_RE_sm_list(n_xg);
  vector<vector<distr_t_list>> F_T_d_MB_IM_sm_list(n_xg);
  vector<vector<distr_t_list>> F_T_d_MB_IM_Gauss_sm_list(n_xg);

  vector<vector<distr_t_list>> F_T_d_MB_RE_sm_list_no_sub(n_xg);
   

  vector<vector<distr_t_list>> F_T_d_MB_RE_sm_list_15(n_xg);
  vector<vector<distr_t_list>> F_T_d_MB_IM_sm_list_15(n_xg);

  vector<vector<distr_t_list>> F_T_d_MB_RE_sm_list_7(n_xg);
  vector<vector<distr_t_list>> F_T_d_MB_IM_sm_list_7(n_xg);

  vector<vector<distr_t_list>> F_T_d_MB_RE_VMD_sm_list(n_xg);
  vector<vector<distr_t_list>> F_T_d_MB_IM_VMD_sm_list(n_xg);

  vector<vector<distr_t_list>> F_T_d_MB_RE_VMD_II_state_sm_list(n_xg);
  vector<vector<distr_t_list>> F_T_d_MB_IM_VMD_II_state_sm_list(n_xg);

  vector<vector<distr_t_list>> F_T_d_MB_RE_VMD_III_state_sm_list(n_xg);
  vector<vector<distr_t_list>> F_T_d_MB_IM_VMD_III_state_sm_list(n_xg);

  vector<distr_t_list> F_T_u_VMD_list;
  vector<distr_t_list> F_T_u_VMD_spectre_list;

  vector<vector<double>> sigma_simulated(n_xg);

  for(int ixg=0; ixg<n_xg;ixg++) {

    sigma_simulated[ixg].resize(sigmas_07.size());

    F_T_u_list.emplace_back(UseJack);
    F_T_d_I_list.emplace_back(UseJack);
    F_T_d_I_ori_list.emplace_back(UseJack);
    F_T_d_I_ori_sp_list.emplace_back(UseJack);
    F_T_d_I_list_15.emplace_back(UseJack);
    F_T_d_I_list_7.emplace_back(UseJack);
    F_T_d_sp_I_list.emplace_back(UseJack);
    FV_T_u_real_list.emplace_back(UseJack);
    FV_T_d_real_list.emplace_back(UseJack);
    FA_T_u_real_list.emplace_back(UseJack);
    FA_T_d_real_list.emplace_back(UseJack);

    FA_T_d_I_real_list.emplace_back(UseJack);
    FA_T_d_II_real_list.emplace_back(UseJack);

    FV_T_d_I_real_list.emplace_back(UseJack);
    FV_T_d_II_real_list.emplace_back(UseJack);

    F_T_u_VMD_list.emplace_back(UseJack);
    F_T_u_VMD_spectre_list.emplace_back(UseJack);

    FV_T_d_sp_real_list.emplace_back(UseJack);
    FA_T_d_sp_real_list.emplace_back(UseJack);

    reminder_F_T_d_list.emplace_back(UseJack);
    reminder_F_T_d_sp_list.emplace_back(UseJack);

    for(int iss=0; iss<(signed)sigmas_07_w0.size(); iss ++) {

   
    F_T_d_MB_RE_VMD_sm_list[ixg].emplace_back(UseJack);
    F_T_d_MB_IM_VMD_sm_list[ixg].emplace_back(UseJack);

    F_T_d_MB_RE_VMD_II_state_sm_list[ixg].emplace_back(UseJack);
    F_T_d_MB_IM_VMD_II_state_sm_list[ixg].emplace_back(UseJack);

    F_T_d_MB_RE_VMD_III_state_sm_list[ixg].emplace_back(UseJack);
    F_T_d_MB_IM_VMD_III_state_sm_list[ixg].emplace_back(UseJack);


    }
    
    for(int is=0; is<(signed)sigmas_07.size(); is++) {
   
      F_T_d_MB_RE_sm_list[ixg].emplace_back(UseJack);
      F_T_d_MB_IM_sm_list[ixg].emplace_back(UseJack);
      F_T_d_MB_IM_Gauss_sm_list[ixg].emplace_back(UseJack);

      F_T_d_MB_RE_sm_list_no_sub[ixg].emplace_back(UseJack);

      F_T_d_MB_RE_sm_list_15[ixg].emplace_back(UseJack);
      F_T_d_MB_IM_sm_list_15[ixg].emplace_back(UseJack);

      F_T_d_MB_RE_sm_list_7[ixg].emplace_back(UseJack);
      F_T_d_MB_IM_sm_list_7[ixg].emplace_back(UseJack);

      
    }
  }

 
  
  
  boost::filesystem::create_directory( path_out+"/corr_2pts");
  boost::filesystem::create_directory( path_out+"/corr_3pts");
  boost::filesystem::create_directory( path_out+"/mass");
  boost::filesystem::create_directory( path_out+"/covariance");
  boost::filesystem::create_directory( path_out+"/FF_d_I");
  boost::filesystem::create_directory( path_out+"/FF_d_II");
  boost::filesystem::create_directory( path_out+"/FF_u");
  boost::filesystem::create_directory( path_out+"/FF_d");
  boost::filesystem::create_directory( path_out+"/FF");


  
  
  //loop over ensembles
  for(int iens=0; iens<Nens;iens++) {

    double Mp_cont;
    
    if(MESON=="B0s") {
      if(Ens_tags[iens] == "cB211b.072.64") {  t_07_s=25; t_07_s_HLT=14; t_07_c=25;    }
      else if(Ens_tags[iens] == "cD211a.054.96") {  t_07_s=35; t_07_s_HLT=22; t_07_c=35;    } //t_07_s_HLT_std=22
      else crash("B0s-crash-virtual");
      Mp_cont = masses[0];
    }
    else if(MESON=="B1s") {
      if(Ens_tags[iens] == "cB211b.072.64") {  t_07_s=25; t_07_s_HLT=14; t_07_c=25;    }
      else if(Ens_tags[iens] == "cD211a.054.96") {  t_07_s=35; t_07_s_HLT=22; t_07_c=35;    }
      else crash("B1s-crash-virtual");
      Mp_cont = masses[1];
    }
    else if(MESON=="B3s") {
      if(Ens_tags[iens] == "cB211b.072.64") {  t_07_s=25; t_07_s_HLT=10; t_07_c=25;    }
      else if(Ens_tags[iens] == "cD211a.054.96") {  t_07_s=35; t_07_s_HLT=21; t_07_c=35;    }
      else crash("B3s-crash-virtual");
      Mp_cont = masses[2];
    }
    else crash("Meson: "+MESON+" not yet simulated");

    boost::filesystem::create_directory( path_out+"/corr");
    boost::filesystem::create_directory( path_out+"/mass/"+data_2pts_SM.Tag[iens]);
    boost::filesystem::create_directory( path_out+"/FF_d_I/"+data_2pts_SM.Tag[iens]);
    boost::filesystem::create_directory( path_out+"/FF_d_II/"+data_2pts_SM.Tag[iens]);
    boost::filesystem::create_directory( path_out+"/FF_d_II/"+data_2pts_SM.Tag[iens]+"/VMD_virt_scan");
    boost::filesystem::create_directory( path_out+"/FF_u/"+data_2pts_SM.Tag[iens]);
    boost::filesystem::create_directory( path_out+"/FF_d/"+data_2pts_SM.Tag[iens]);
    boost::filesystem::create_directory( path_out+"/FF/"+data_2pts_SM.Tag[iens]);
    
    CorrAnalysis Corr(UseJack,Njacks,1000);
    Corr.Nt = Nts[iens];
    CorrAnalysis Corr_boot(false,Njacks, 1000);
    Corr_boot.Nt= Nts[iens];
   

    LatticeInfo L_info;
    L_info.LatInfo_new_ens(data_2pts_SM.Tag[iens]);
   

    //read theta values and loop over them
    Vfloat thetas, masses_u, masses_d;

    thetas= Read_From_File(Corr_path+"/spectre/"+Ens_tags[iens]+"/pars_list.dat", 1 , 5);
    masses_u= Read_From_File(Corr_path+"/spectre/"+Ens_tags[iens]+"/pars_list.dat", 3 , 5);
    masses_d= Read_From_File( Corr_path+"/spectre/"+Ens_tags[iens]+"/pars_list.dat", 4 , 5);

    double mu= masses_u[0];
    double md= masses_d[0];


    //RCs
    
    distr_t ZT, a_distr, ZV;
    distr_t ZT_boot;
    if(data_2pts_SM.Tag[iens].substr(1,1)=="A") { ZT=ZT_A; ZV=ZV_A; ZT_boot = ZT_A_boot; a_distr=a_A;}
    else if(data_2pts_SM.Tag[iens].substr(1,1)=="B") { ZT= ZT_B; ZV=ZV_B; ZT_boot = ZT_B_boot; a_distr=a_B;}
    else if(data_2pts_SM.Tag[iens].substr(1,1)=="C") { ZT= ZT_C; ZV=ZV_C; ZT_boot = ZT_C_boot; a_distr=a_C;}
    else if(data_2pts_SM.Tag[iens].substr(1,1)=="D") { ZT= ZT_D; ZV=ZV_D; ZT_boot = ZT_D_boot; a_distr=a_D;}
    else crash("Ensemble: "+data_2pts_SM.Tag[iens]+" not recognised");


    if(CONS_EM_CURRENT==false) ZT = ZT*ZV;


    //2pts plateaux
    if(data_2pts_SM.Tag[iens].substr(1,1) =="A") {
      if(MESON == "B1s") { Corr.Tmin=14; Corr.Tmax= 24;}
      else if(MESON == "B2s" ) { Corr.Tmin=13; Corr.Tmax= 24;}
      else if(MESON=="B3s") { Corr.Tmin=14; Corr.Tmax=20;}
      else if(MESON=="B4s") { Corr.Tmin=14; Corr.Tmax=19;}
      else { Corr.Tmin=14; Corr.Tmax=27;}
      
    }
    else if(data_2pts_SM.Tag[iens] =="cB211b.072.64") {
      if(MESON == "B1s") {Corr.Tmin=24; Corr.Tmax=38;}
      else if(MESON == "B2s") {Corr.Tmin=22; Corr.Tmax=34;}
      else if(MESON=="B3s") {Corr.Tmin=17; Corr.Tmax=23;}
      else if(MESON=="B4s") { Corr.Tmin=16; Corr.Tmax=22;}
      else {Corr.Tmin=23; Corr.Tmax=33;}
    }
    else if(data_2pts_SM.Tag[iens] =="cB211b.072.96") {Corr.Tmin=20; Corr.Tmax=36;}
    
    else if(data_2pts_SM.Tag[iens].substr(1,1) == "C")  {
      if(MESON == "B1s") {Corr.Tmin=27; Corr.Tmax=36;}
      else if(MESON == "B2s") {Corr.Tmin=27; Corr.Tmax=36;}
      else if(MESON=="B3s") {Corr.Tmin=18; Corr.Tmax=31;}
      else if(MESON=="B4s") { Corr.Tmin=18; Corr.Tmax=30;}
      else {Corr.Tmin=29; Corr.Tmax=42;}
    }
    else if(data_2pts_SM.Tag[iens].substr(1,1) == "D")  {
      if(MESON == "B1s") {Corr.Tmin=32; Corr.Tmax=53;}
      else if(MESON == "B2s") { Corr.Tmin=28; Corr.Tmax= 37;}
      else if(MESON=="B3s") {Corr.Tmin=28; Corr.Tmax=37;}
      else if(MESON=="B4s") { Corr.Tmin=27; Corr.Tmax=37;}
      else {Corr.Tmin=38; Corr.Tmax=52;}
    }
    else crash("In fixing [Tmin, Tmax] for smeared MP, Ensemble: "+data_2pts_SM.Tag[iens]+" not recognized");

    Corr_boot.Tmin= Corr.Tmin; Corr_boot.Tmax=Corr.Tmax;


    //mass and decay constants
    distr_t_list pt2_distr_SM= Corr.corr_t(data_2pts_SM.col(0)[iens], path_out+"/corr_2pts/"+data_2pts_SM.Tag[iens]+"/corr_2pt_SM.dat");
    distr_t_list eff_mass_SM = Corr.effective_mass_t(pt2_distr_SM, path_out+"/mass/"+data_2pts_SM.Tag[iens]+"/eff_mass_SM.dat");
    distr_t_list pt2_distr_SMSM= Corr.corr_t(data_2pts_SMSM.col(0)[iens], path_out+"/corr_2pts/"+data_2pts_SM.Tag[iens]+"/corr_2pt_SMSM.dat");
    distr_t_list eff_mass_SMSM = Corr.effective_mass_t(pt2_distr_SMSM, path_out+"/mass/"+data_2pts_SM.Tag[iens]+"/eff_mass_SMSM.dat");
    distr_t mel_SMSM= Corr.Fit_distr(Corr.mel_ov_mass_t(pt2_distr_SMSM, ""))/2.0;

    distr_t_list pt2_distr_SM_boot= Corr_boot.corr_t(data_2pts_SM.col(0)[iens], "");
    distr_t_list eff_mass_SM_boot = Corr_boot.effective_mass_t(pt2_distr_SM_boot, "");
    distr_t_list pt2_distr_SMSM_boot= Corr_boot.corr_t(data_2pts_SMSM.col(0)[iens], "");
    distr_t_list eff_mass_SMSM_boot = Corr_boot.effective_mass_t(pt2_distr_SMSM_boot,"");


    distr_t M_P_SM = Corr.Fit_distr(eff_mass_SM);
    distr_t M_P_SMSM = Corr.Fit_distr(eff_mass_SMSM);
    distr_t M_P= 0.5*(M_P_SM+M_P_SMSM);
    auto SINH= [](double x) { return sinh(x);};
    distr_t_list FP_SM_distr_list = (mu + md )*Corr.residue_t( pt2_distr_SM, "")/(M_P*distr_t::f_of_distr(SINH, M_P)*Corr.matrix_element_t(pt2_distr_SMSM, ""));
    Print_To_File({}, {FP_SM_distr_list.ave(), FP_SM_distr_list.err()}, path_out+"/decay_const/"+data_2pts_SM.Tag[iens]+"/decay_const_SM.dat.t", "", "");
    distr_t FP_SM= Corr.Fit_distr( FP_SM_distr_list  );

    distr_t M_P_boot = 0.5*( Corr_boot.Fit_distr( eff_mass_SMSM_boot) + Corr_boot.Fit_distr(eff_mass_SM_boot));

  
    distr_t mel_SMSM_boot= Corr_boot.Fit_distr(Corr_boot.mel_ov_mass_t(pt2_distr_SMSM_boot, ""))/2.0;
  


    cout<<"MP: "<<(M_P/a_distr).ave()<<" +- "<<(M_P/a_distr).err()<<endl;
    cout<<"FP_SM: "<<(FP_SM/a_distr).ave()<<" +- "<<(FP_SM/a_distr).err()<<endl;

    
    distr_t F_P= FP_SM;

    //loop over photon momentum
    for(int ixg=0;ixg<n_xg;ixg++) {

      //get xg, Eg, kz from thetas
      double theta=thetas[ixg];
      pt3_momenta pt3_mom_07(0.0, 0.0, thetas[ixg]/2.0, masses_u[ixg], masses_d[ixg], 0.0, L_info.L, L_info.T);
      double Eg= pt3_mom_07.Egamma();
      distr_t Eg_off = M_P - Eg;
      double kz = pt3_mom_07.k()[2];
      distr_t xg= pt3_mom_07.x_gamma(M_P);
      distr_t xg_boot= pt3_mom_07.x_gamma(M_P_boot);
      distr_t Eg_MB_off = MBs*a_distr*(1.0 -xg/2.0);
      //double y_eff= y*(1-xg.ave()/2.0)/(1-0.1/2.0);
     
      distr_t Eg_MB_off_new= (y_eff[ixg]*Mp_cont + (1-y_eff[ixg])*MBs)*(1.0 -0.5*xg);
      Eg_MB_off = y_eff[ixg]*Eg_off + (1-y_eff[ixg])*Eg_MB_off;
      //xg_list.distr_list.push_back(xg);

      
    
      cout<<"##### Considering kinematic with..."<<endl;
      cout<<"Eg (on-shell): "<<Eg<<endl;
      cout<<"Eg (virt): "<<Eg_off.ave()<<" +- "<<Eg_off.err()<<endl;
      cout<<"thz: "<<theta<<endl;
      cout<<"kz: "<<kz<<endl;
      cout<<"xg: "<<xg.ave()<<" +- "<<xg.err()<<endl;
      int Im_Re;

      //tensor- electric part
      Corr.Reflection_sign= 1;
      Im_Re=1;
      Corr.Perform_Nt_t_average = 0;
      Corr_boot.Reflection_sign=1;
      Corr_boot.Perform_Nt_t_average=0;
	 
	 
      distr_t_list T_u_std = 0.5*QU*Corr.corr_t(summ_master(C_T_u_data_std[1][2][ixg].col(Im_Re)[iens], Multiply_Vvector_by_scalar(C_T_u_data_std[2][1][ixg].col(Im_Re)[iens], -1.0)), path_out+"/corr_3pts/"+TAG_CURR+""+data_2pts_SM.Tag[iens]+"_T_u_xg_"+to_string(ixg));
      	 
      distr_t_list T_d_std = 0.5*QD*Corr.corr_t(summ_master(C_T_d_data_std[1][2][ixg].col(Im_Re)[iens], Multiply_Vvector_by_scalar(C_T_d_data_std[2][1][ixg].col(Im_Re)[iens], -1.0)), path_out+"/corr_3pts/"+TAG_CURR+""+data_2pts_SM.Tag[iens]+"_T_d_xg_"+to_string(ixg));
      
      distr_t_list T_d = 0.5*QD*Corr.corr_t(summ_master(C_T_d_data[1][2][ixg].col(Im_Re)[iens], Multiply_Vvector_by_scalar(C_T_d_data[2][1][ixg].col(Im_Re)[iens], -1.0)), path_out+"/corr_3pts/"+TAG_CURR+""+data_2pts_SM.Tag[iens]+"_T_d_sp_xg_"+to_string(ixg));
      distr_t_list T_d_boot= 0.5*QD*Corr_boot.corr_t(summ_master(C_T_d_data[1][2][ixg].col(Im_Re)[iens], Multiply_Vvector_by_scalar(C_T_d_data[2][1][ixg].col(Im_Re)[iens], -1.0)),"");

      Im_Re=0;

      distr_t_list B_u_std = 0.5*QU*Corr.corr_t(summ_master(C_B_u_data_std[1][1][ixg].col(Im_Re)[iens], C_B_u_data_std[2][2][ixg].col(Im_Re)[iens]), path_out+"/corr_3pts/"+TAG_CURR+""+data_2pts_SM.Tag[iens]+"_B_u_xg_"+to_string(ixg));
      distr_t_list B_d_std = 0.5*QD*Corr.corr_t(summ_master(C_B_d_data_std[1][1][ixg].col(Im_Re)[iens], C_B_d_data_std[2][2][ixg].col(Im_Re)[iens]), path_out+"/corr_3pts/"+TAG_CURR+""+data_2pts_SM.Tag[iens]+"_B_d_xg_"+to_string(ixg));
      distr_t_list B_d = 0.5*QD*Corr.corr_t(summ_master(C_B_d_data[1][1][ixg].col(Im_Re)[iens], C_B_d_data[2][2][ixg].col(Im_Re)[iens]), path_out+"/corr_3pts/"+TAG_CURR+""+data_2pts_SM.Tag[iens]+"_B_d_sp_xg_"+to_string(ixg));
      distr_t_list B_d_boot= 0.5*QD*Corr_boot.corr_t(summ_master(C_B_d_data[1][1][ixg].col(Im_Re)[iens], C_B_d_data[2][2][ixg].col(Im_Re)[iens]),"");


     
      //determine the form factors that do not need spectral-reconstruction techniques, i.e. up-type quark (heavy) contribution and down-type quark (s-quark) contribution in I-TO


      //#############     VMD PRED  FT_u  #####################
      //define corr for Fu in 2nd TO for VMD pred
      distr_t_list Corr_Tu_2TO(UseJack);
      for(int t=t_07_c; t <=Corr.Nt/2;t++) {
	Corr_Tu_2TO.distr_list.push_back( (xg/2.0)*(sign_kz*T_u_std.distr_list[t] + B_u_std.distr_list[t]));
      }
      CorrAnalysis Corr_VMD_anal(UseJack, Njacks,1000);
      if(MESON=="B0s") {
	
	if(Ens_tags[iens] == "cB211b.072.64") {
	  Corr_VMD_anal.Tmin= 17;
	  Corr_VMD_anal.Tmax= 30;
	}
	else if(Ens_tags[iens] == "cD211a.054.96") {
	  Corr_VMD_anal.Tmin=25;
	  Corr_VMD_anal.Tmax=46;
	}
	else crash("Ensemble: "+Ens_tags[iens]+" not yet implemented");
      }
      else if(MESON=="B1s") {
	if(Ens_tags[iens] == "cB211b.072.64") {
	Corr_VMD_anal.Tmin= 17;
	Corr_VMD_anal.Tmax= 30;
	}
	else if(Ens_tags[iens] == "cD211a.054.96") {
	  Corr_VMD_anal.Tmin=27;
	  Corr_VMD_anal.Tmax=50;
	}
	else crash("Ensemble: "+Ens_tags[iens]+" not yet implemented");
      }
      else if(MESON=="B3s") {
	if(Ens_tags[iens] == "cB211b.072.64") {
	  Corr_VMD_anal.Tmin= 17;
	  Corr_VMD_anal.Tmax= 30;
	}
	else if(Ens_tags[iens] == "cD211a.054.96") {
	  Corr_VMD_anal.Tmin=25;
	  Corr_VMD_anal.Tmax=50;
	}
	else crash("Ensemble: "+Ens_tags[iens]+" not yet implemented");
      }
      else crash("Meson: "+MESON+" not yet implemented");
      
      Corr_VMD_anal.Nt=2*(Corr_Tu_2TO.size()-1);
      distr_t_list Corr_Tu_2TO_symm= Corr_Tu_2TO;
      for(int t=Corr_VMD_anal.Nt/2 +1;t<Corr_VMD_anal.Nt;t++) Corr_Tu_2TO_symm.distr_list.push_back( Corr_Tu_2TO.distr_list[Corr_VMD_anal.Nt -t]);
      distr_t_list eff_M_Tu_distr= Corr_VMD_anal.effective_mass_t(Corr_Tu_2TO_symm, "" );
      distr_t eff_M_Tu= Corr_VMD_anal.Fit_distr(eff_M_Tu_distr)/a_distr;
      Print_To_File({}, {eff_M_Tu_distr.ave(), eff_M_Tu_distr.err(), (eff_M_Tu_distr-eff_M_Tu_distr[30]).ave(),  (eff_M_Tu_distr-eff_M_Tu_distr[30]).err() }, path_out+"/mass/"+Ens_tags[iens]+"/"+TAG_CURR+"Vb_eff_xg_"+to_string_with_precision(xg_list.ave(ixg),2), "", "");
      distr_t Ups_MT = Corr_VMD_anal.Fit_distr( Corr_Tu_2TO_symm*(EXPT_DL(eff_M_Tu_distr)));

       auto K_RE_distr_u_VMD= [](const distr_t &E, const distr_t &m, double s) -> distr_t {
	     
	     distr_t ret(1);
	   
	       for(int ijack=0;ijack<E.size();ijack++) {
		 double x= (E.distr[ijack]-m.distr[ijack]);
		 double y= E.distr[ijack];
		 ret.distr.push_back( 2*sinh(x)*cos(s)/(cosh(2*x) - cos(2*s)) - 2*sinh(y)*cos(0)/(cosh(2*y) - cos(2*0))  );
	       }
	   	     
	     return ret;
	     
       };

      
       F_T_u_VMD_spectre_list[ixg].distr_list.push_back( Ups_MT*ZT*(1.0/(mel_SMSM*Eg))*EXP_D( M_P*t_07_c)*K_RE_distr_u_VMD( eff_M_Tu*a_distr, Eg_off, 1e-5));
       distr_t_list Corr_F_T_u_VMD = Ups_MT*EXPT_D(-1.0*eff_M_Tu*a_distr, Corr_Tu_2TO.size());
       distr_t F_T_u_VMD_distr(UseJack, Njacks);
       for(int ty=1;ty<Corr_Tu_2TO.size(); ty++) {
	 F_T_u_VMD_distr = F_T_u_VMD_distr + Corr_F_T_u_VMD[ty]*EXP_D(ty*Eg_off);
       }
       F_T_u_VMD_list[ixg].distr_list.push_back( F_T_u_VMD_distr*ZT*(1.0/(mel_SMSM*Eg))*EXP_D( M_P*t_07_c));
       
       
       
       //#######################################################
       
      
      distr_t F_T_u(UseJack, Njacks);
      distr_t F_T_d_I(UseJack, Njacks);
      distr_t F_T_d_I_ori(UseJack,Njacks);
      distr_t F_T_d_I_ori_sp(UseJack,Njacks);
      distr_t F_T_d_I_15(UseJack, Njacks);
      distr_t F_T_d_I_7(UseJack, Njacks);
      
      distr_t F_T_d_I_sp(UseJack, Njacks);
      
      distr_t FV_T_u_real(UseJack, Njacks);
      distr_t FV_T_d_real(UseJack, Njacks);
      distr_t FA_T_u_real(UseJack, Njacks);
      distr_t FA_T_d_real(UseJack, Njacks);
      distr_t FV_T_d_sp_real(UseJack,Njacks);
      distr_t FA_T_d_sp_real(UseJack,Njacks);

      distr_t FA_T_d_I_real(UseJack, Njacks);
      distr_t FA_T_d_II_real(UseJack, Njacks);

      distr_t FV_T_d_I_real(UseJack, Njacks);
      distr_t FV_T_d_II_real(UseJack, Njacks);

      distr_t_list F_T_d_ty(UseJack), F_T_d_sp_ty(UseJack), F_T_u_ty(UseJack);
      distr_t_list F_T_d_ty_psum(UseJack), F_T_d_sp_ty_psum(UseJack), F_T_u_ty_psum(UseJack);

      distr_t_list FA_T_d_real_ty(UseJack), FA_T_d_sp_real_ty(UseJack), FA_T_u_real_ty(UseJack);
      distr_t_list FA_T_d_real_ty_psum(UseJack), FA_T_d_sp_real_ty_psum(UseJack), FA_T_u_real_ty_psum(UseJack);

      distr_t_list FV_T_d_real_ty(UseJack), FV_T_d_sp_real_ty(UseJack), FV_T_u_real_ty(UseJack);
      distr_t_list FV_T_d_real_ty_psum(UseJack), FV_T_d_sp_real_ty_psum(UseJack), FV_T_u_real_ty_psum(UseJack);

      //reminder of

      distr_t reminder_FV_T_d_I_real(UseJack,Njacks), reminder_FA_T_d_I_real(UseJack,Njacks) , reminder_F_T_d_I(UseJack,Njacks), reminder_F_T_d_I_sp(UseJack,Njacks);

      distr_t TO2_reminder(UseJack,Njacks), TO2_reminder_sp(UseJack,Njacks);

      distr_t reminder_F_T_d_I_15(UseJack,Njacks), reminder_F_T_d_I_7(UseJack,Njacks);

      distr_t_list C_0_corr(UseJack);

      
      auto HeavyTheta=[](const int x) {  return ((x>=0)+(x>0))/2.0;   };
      auto Exp= [&Njacks, &UseJack](const distr_t &A) -> distr_t { distr_t ret(UseJack); for(int ijack=0;ijack<Njacks;ijack++) ret.distr.push_back( exp(A.distr[ijack])); return ret;};
      double T=Corr.Nt;
      
      for(int ty=0; ty <= Corr.Nt/2; ty++) {
	
	const distr_t f1=  Exp(-(T/2-ty)*Eg_off);

	const distr_t f2= Exp(-((3*T/2)-ty)*Eg_off);

	const distr_t f1_s=  Exp(-(T/2-ty)*Eg_off);

	const distr_t f2_s= Exp(-((3*T/2)-ty)*Eg_off);

	const double f1_real= exp(-(T/2-ty)*Eg);
	const double f2_real= exp(-((3*T/2)-ty)*Eg);

	const double h1=HeavyTheta((T/2)-ty);
	
	const double h2=HeavyTheta(ty-(T/2));

	const distr_t f1_sub=  Exp(-(T/2-ty)*Eg_off) - Exp( -1.0*Eg_off*abs(T/2 - t_07_s));
	const distr_t f2_sub=  Exp(-((3*T/2)-ty)*Eg_off) - Exp( -1.0*Eg_off*abs(T/2 - t_07_s));

	const distr_t f1_sub_15=  Exp(-(T/2-ty)*1.5*a_distr) - Exp( -1.0*1.5*a_distr*abs(T/2 - t_07_s));
	const distr_t f2_sub_15=  Exp(-((3*T/2)-ty)*1.5*a_distr) - Exp( -1.0*1.5*a_distr*abs(T/2 - t_07_s));

	const distr_t f1_sub_7=  Exp(-(T/2-ty)*0.7*a_distr) - Exp( -1.0*0.7*a_distr*abs(T/2 - t_07_s));
	const distr_t f2_sub_7=  Exp(-((3*T/2)-ty)*0.7*a_distr) - Exp( -1.0*0.7*a_distr*abs(T/2 - t_07_s));

	const distr_t f1_sub_sp= Exp(-(T/2-ty)*Eg_off) - Exp( -1.0*Eg_off*abs(T/2 - t_07_s_HLT));
	const distr_t f2_sub_sp= Exp(-((3*T/2)-ty)*Eg_off) - Exp( -1.0*Eg_off*abs(T/2 - t_07_s_HLT));

	const double f1_real_sub= f1_real - exp( -1.0*Eg*abs(T/2 - t_07_s));
	const double f2_real_sub= f2_real - exp( -1.0*Eg*abs(T/2 - t_07_s));



	F_T_d_ty.distr_list.push_back( (sign_kz*T_d_std.distr_list[ty]*(xg/2.0)   + B_d_std.distr_list[ty]*xg/2.0 )*(h1*f1+ h2*f2)            );
	F_T_d_sp_ty.distr_list.push_back( (sign_kz*T_d.distr_list[ty]*(xg/2.0) + B_d.distr_list[ty]*xg/2.0)*(h1*f1+h2*f2)   );
	F_T_u_ty.distr_list.push_back(  (sign_kz*T_u_std.distr_list[ty]*(xg/2.0)   + B_u_std.distr_list[ty]*xg/2.0 )*(h1*f1+ h2*f2));


	F_T_d_ty_psum.distr_list.push_back( ((ty==0)?F_T_d_ty[ty]:(F_T_d_ty_psum[ty-1] + F_T_d_ty[ty])));
	F_T_d_sp_ty_psum.distr_list.push_back( ((ty==0)?F_T_d_sp_ty[ty]:(F_T_d_sp_ty_psum[ty-1] + F_T_d_sp_ty[ty])));
	F_T_u_ty_psum.distr_list.push_back( ((ty==0)?F_T_u_ty[ty]:(F_T_u_ty_psum[ty-1] + F_T_u_ty[ty])));



	FA_T_d_real_ty.distr_list.push_back( ( B_d_std.distr_list[ty]*(1.0- xg/2.0)   + sign_kz*T_d_std.distr_list[ty]*xg/2.0 )*(h1*f1_real+ h2*f2_real)    )      ;

        FA_T_d_sp_real_ty.distr_list.push_back( ( B_d.distr_list[ty]*(1.0- xg/2.0)   + sign_kz*T_d.distr_list[ty]*xg/2.0 )*(h1*f1_real+ h2*f2_real));
	FA_T_u_real_ty.distr_list.push_back(  ( B_u_std.distr_list[ty]*(1.0- xg/2.0)   + sign_kz*T_u_std.distr_list[ty]*xg/2.0 )*(h1*f1_real+ h2*f2_real));



        FA_T_d_real_ty_psum.distr_list.push_back( ((ty==0)?FA_T_d_real_ty[ty]:(FA_T_d_real_ty_psum[ty-1] + FA_T_d_real_ty[ty])));
	FA_T_d_sp_real_ty_psum.distr_list.push_back( ((ty==0)?FA_T_d_sp_real_ty[ty]:(FA_T_d_sp_real_ty_psum[ty-1] + FA_T_d_sp_real_ty[ty])));
	FA_T_u_real_ty_psum.distr_list.push_back( ((ty==0)?FA_T_u_real_ty[ty]:(FA_T_u_real_ty_psum[ty-1] + FA_T_u_real_ty[ty])));



	FV_T_d_real_ty.distr_list.push_back( (sign_kz*T_d_std.distr_list[ty]*(1.0- xg/2.0)   + B_d_std.distr_list[ty]*xg/2.0 )*(h1*f1_real+ h2*f2_real));
	FV_T_d_sp_real_ty.distr_list.push_back( (sign_kz*T_d.distr_list[ty]*(1.0- xg/2.0)   + B_d.distr_list[ty]*xg/2.0 )*(h1*f1_real+ h2*f2_real));
	FV_T_u_real_ty.distr_list.push_back(  (sign_kz*T_u_std.distr_list[ty]*(1.0- xg/2.0)   + B_u_std.distr_list[ty]*xg/2.0 )*(h1*f1_real+ h2*f2_real));


	FV_T_d_real_ty_psum.distr_list.push_back( ((ty==0)?FV_T_d_real_ty[ty]:(FV_T_d_real_ty_psum[ty-1] + FV_T_d_real_ty[ty])));
	FV_T_d_sp_real_ty_psum.distr_list.push_back( ((ty==0)?FV_T_d_sp_real_ty[ty]:(FV_T_d_sp_real_ty_psum[ty-1] + FV_T_d_sp_real_ty[ty])));
	FV_T_u_real_ty_psum.distr_list.push_back( ((ty==0)?FV_T_u_real_ty[ty]:(FV_T_u_real_ty_psum[ty-1] + FV_T_u_real_ty[ty])));

	//if(ty != t_07_c) {
	//if(ty== t_07_c +1 || ty== t_07_c-1)  F_T_u = F_T_u +  0.5*(sign_kz*T_u_std.distr_list[ty]*(xg/2.0)   + B_u_std.distr_list[ty]*xg/2.0 )*(h1*f1+ h2*f2);
	// else F_T_u = F_T_u +  (sign_kz*T_u_std.distr_list[ty]*(xg/2.0)   + B_u_std.distr_list[ty]*xg/2.0 )*(h1*f1+ h2*f2);
	// }
	F_T_u = F_T_u +  (sign_kz*T_u_std.distr_list[ty]*(xg/2.0)   + B_u_std.distr_list[ty]*xg/2.0 )*(h1*f1+ h2*f2);
	
	FV_T_u_real = FV_T_u_real + (sign_kz*T_u_std.distr_list[ty]*(1.0- xg/2.0)   + B_u_std.distr_list[ty]*xg/2.0 )*(h1*f1_real+ h2*f2_real);
	FV_T_d_real = FV_T_d_real + (sign_kz*T_d_std.distr_list[ty]*(1.0- xg/2.0)   + B_d_std.distr_list[ty]*xg/2.0 )*(h1*f1_real+ h2*f2_real);
	FA_T_u_real = FA_T_u_real + ( B_u_std.distr_list[ty]*(1.0- xg/2.0)   + sign_kz*T_u_std.distr_list[ty]*xg/2.0 )*(h1*f1_real+ h2*f2_real);
	FA_T_d_real = FA_T_d_real + ( B_d_std.distr_list[ty]*(1.0- xg/2.0)   + sign_kz*T_d_std.distr_list[ty]*xg/2.0 )*(h1*f1_real+ h2*f2_real);
	FA_T_d_sp_real = FA_T_d_sp_real + ( B_d.distr_list[ty]*(1.0- xg/2.0)   + sign_kz*T_d.distr_list[ty]*xg/2.0 )*(h1*f1_real+ h2*f2_real);
	FV_T_d_sp_real = FV_T_d_sp_real + (sign_kz*T_d.distr_list[ty]*(1.0- xg/2.0)   + B_d.distr_list[ty]*xg/2.0 )*(h1*f1_real+ h2*f2_real);

	if(ty<=t_07_s) {
	FA_T_d_I_real = FA_T_d_I_real + ( B_d_std.distr_list[ty]*(1.0- xg/2.0)   + sign_kz*T_d_std.distr_list[ty]*xg/2.0 )*(h1*f1_real_sub+ h2*f2_real_sub);
	FV_T_d_I_real = FV_T_d_I_real +  (sign_kz*T_d_std.distr_list[ty]*(1.0- xg/2.0)   + B_d_std.distr_list[ty]*xg/2.0 )*(h1*f1_real_sub+ h2*f2_real_sub);
	}
	else {
	  FA_T_d_II_real = FA_T_d_II_real + ( B_d_std.distr_list[ty]*(1.0- xg/2.0)   + sign_kz*T_d_std.distr_list[ty]*xg/2.0 )*(h1*f1_real_sub+ h2*f2_real_sub);
	  FV_T_d_II_real = FV_T_d_II_real +  (sign_kz*T_d_std.distr_list[ty]*(1.0- xg/2.0)   + B_d_std.distr_list[ty]*xg/2.0 )*(h1*f1_real_sub+ h2*f2_real_sub);
	}
	
	
	if(ty<= t_07_s) {  //SHOULD BE <= t_07_s

	 

	  //interpolate the correlator at t=t_07_s
	  distr_t corr_tw_minus(UseJack);
	  distr_t corr_tw_plus(UseJack);
	  distr_t corr_tw_minus_2(UseJack);
	  distr_t corr_tw_plus_2(UseJack);
	  if(ty==t_07_s) {
	    
	    int t1= t_07_s-1;
	    int t2= t_07_s-2;
	    int t3= t_07_s-3;

	    int t4= t_07_s+1;
	    int t5= t_07_s+2;
	    int t6= t_07_s+3;

	    distr_t C0= (sign_kz*T_d_std.distr_list[ty]*(xg/2.0)   + B_d_std.distr_list[ty]*xg/2.0);

	    distr_t C1= (sign_kz*T_d_std.distr_list[t1]*(xg/2.0)   + B_d_std.distr_list[t1]*xg/2.0);
	    distr_t C2= (sign_kz*T_d_std.distr_list[t2]*(xg/2.0)   + B_d_std.distr_list[t2]*xg/2.0);
	    distr_t C3= (sign_kz*T_d_std.distr_list[t3]*(xg/2.0)   + B_d_std.distr_list[t3]*xg/2.0);
	    
	    distr_t C4= (sign_kz*T_d_std.distr_list[t4]*(xg/2.0)   + B_d_std.distr_list[t4]*xg/2.0);
	    distr_t C5= (sign_kz*T_d_std.distr_list[t5]*(xg/2.0)   + B_d_std.distr_list[t5]*xg/2.0);
	    distr_t C6= (sign_kz*T_d_std.distr_list[t6]*(xg/2.0)   + B_d_std.distr_list[t6]*xg/2.0);

	   
	    //if(C1.ave()*C2.ave() < 0 || C2.ave()*C3.ave() < 0) crash("While finding C_s(tw-), correlators do not have same sign");
	    //if(C4.ave()*C5.ave() < 0 || C5.ave()*C6.ave() < 0) crash("While finding C_s(tw+), correlators do not have same sign");
	    int sign=1;
	    if(C1.ave()<0) sign=-1;
	    
	    //for(int ijack=0;ijack<Njacks;ijack++) corr_tw_minus.distr.push_back( sign*exp(  quad_interpolator(log(fabs(C1.distr[ijack])), log(fabs(C2.distr[ijack])), log(fabs(C3.distr[ijack])), -1,-2,-3,0)));
	    
	    for(int ijack=0;ijack<Njacks;ijack++) corr_tw_minus.distr.push_back(  quad_interpolator(C1.distr[ijack], C2.distr[ijack], C3.distr[ijack], -1,-2,-3,0));
	    for(int ijack=0;ijack<Njacks;ijack++) corr_tw_minus_2.distr.push_back(  lin_interpolator(C1.distr[ijack], C2.distr[ijack],-1,-2,0));
	  
	    cout<<"sign for C(tw-) ixg: "<<ixg<<" MESON: "<<MESON<<" : "<<sign<<endl;
	    cout<<"C1: "<<(C1).ave()<<" C2: "<<(C2).ave()<<" C3: "<<(C3).ave()<<endl;
	    cout<<"C_minus quadratic: "<<(corr_tw_minus).ave()<<endl;
	    cout<<"C_minus linear: "<<(corr_tw_minus_2).ave()<<endl;

	    sign=1;
	    if(C4.ave()<0) sign=-1;
	    
	    //for(int ijack=0;ijack<Njacks;ijack++) corr_tw_plus.distr.push_back( sign*exp(  quad_interpolator(log(fabs(C4.distr[ijack])), log(fabs(C5.distr[ijack])), log(fabs(C6.distr[ijack])), 1,2,3,0)));

	    for(int ijack=0;ijack<Njacks;ijack++) corr_tw_plus.distr.push_back(  quad_interpolator(C4.distr[ijack], C5.distr[ijack], C6.distr[ijack], 1,2,3,0));
	    for(int ijack=0;ijack<Njacks;ijack++) corr_tw_plus_2.distr.push_back(  lin_interpolator(C4.distr[ijack], C5.distr[ijack],1,2,0));
	   
	    
	    cout<<"sign for C(tw+): ixg"<<ixg<<" MESON: "<<MESON<<" : "<<sign<<endl;
	    cout<<"C4: "<<(C4).ave()<<" C5: "<<(C5).ave()<<" C6: "<<(C6).ave()<<endl;
	    cout<<"C_plus quadratic: "<<corr_tw_plus.ave()<<endl;
	    cout<<"C_plus linear: "<<corr_tw_plus_2.ave()<<endl;

	    double syst_min= 0.5*fabs(corr_tw_minus.ave()-corr_tw_minus_2.ave());
	    corr_tw_minus= 0.5*(corr_tw_minus+corr_tw_minus_2);
	    corr_tw_minus= corr_tw_minus.ave() + (corr_tw_minus -corr_tw_minus.ave())*sqrt( pow(corr_tw_minus.err(),2)+ pow(syst_min,2))/corr_tw_minus.err();

	    double syst_max= 0.5*fabs(corr_tw_plus.ave()-corr_tw_plus_2.ave());
	    corr_tw_plus= 0.5*(corr_tw_plus+corr_tw_plus_2);
	    corr_tw_plus= corr_tw_plus.ave() + (corr_tw_plus -corr_tw_plus.ave())*sqrt( pow(corr_tw_plus.err(),2)+ pow(syst_max,2))/corr_tw_plus.err();

	    cout<<"C(tw-)+C(tw+)/2 : "<<(corr_tw_minus+corr_tw_plus).ave()*0.5<<" +- "<<(corr_tw_minus+corr_tw_plus).err()*0.5<<endl;
	    cout<<"C(0): "<<(C0).ave()<<" +- "<<(C0).err()<<endl;

	    corr_tw_minus= corr_tw_minus*(Exp(-(T/2-ty)*Eg_off));
	    corr_tw_plus= corr_tw_plus*(Exp(-(T/2-ty)*Eg_off));


	    C_0_corr.distr_list.push_back( corr_tw_minus);
	    C_0_corr.distr_list.push_back( 0.5*(corr_tw_minus+corr_tw_plus));
	    C_0_corr.distr_list.push_back( corr_tw_plus);


	    
	    
	  }
	  else corr_tw_minus = 0.0*Get_id_distr(Njacks, UseJack);
	  
	  
	  F_T_d_I = F_T_d_I + (sign_kz*T_d_std.distr_list[ty]*(xg/2.0)   + B_d_std.distr_list[ty]*xg/2.0)*(h1*f1_sub+ h2*f2_sub);
	  F_T_d_I_ori = F_T_d_I_ori + ((ty != t_07_s)?(sign_kz*T_d_std.distr_list[ty]*(xg/2.0)   + B_d_std.distr_list[ty]*xg/2.0)*(h1*f1_s+ h2*f2_s):0.5*corr_tw_minus);
	  F_T_d_I_15 = F_T_d_I_15 + (sign_kz*T_d_std.distr_list[ty]*(xg/2.0)   + B_d_std.distr_list[ty]*xg/2.0 )*(h1*f1_sub_15+ h2*f2_sub_15);
	  F_T_d_I_7 = F_T_d_I_7 + (sign_kz*T_d_std.distr_list[ty]*(xg/2.0)   + B_d_std.distr_list[ty]*xg/2.0 )*(h1*f1_sub_7+ h2*f2_sub_7);
	}
	
	if(ty <= t_07_s_HLT) {	  
	  F_T_d_I_sp = F_T_d_I_sp + (sign_kz*T_d.distr_list[ty]*(xg/2.0) + B_d.distr_list[ty]*xg/2.0)*(h1*f1_sub_sp+h2*f2_sub_sp);
	  F_T_d_I_ori_sp = F_T_d_I_ori_sp + (sign_kz*T_d.distr_list[ty]*(xg/2.0) + B_d.distr_list[ty]*xg/2.0)*(h1*f1_s+h2*f2_s);
	}

	//collect reminders
	reminder_FA_T_d_I_real = reminder_FA_T_d_I_real + ( B_d_std.distr_list[ty]*(1.0- xg/2.0)   + sign_kz*T_d_std.distr_list[ty]*xg/2.0 )*(h1 + h2)*exp( -Eg*abs(T/2 - t_07_s));
	reminder_FV_T_d_I_real = reminder_FV_T_d_I_real + ( sign_kz*T_d_std.distr_list[ty]*(1.0- xg/2.0)   + B_d_std.distr_list[ty]*xg/2.0 )*(h1 + h2)*exp( -Eg*abs(T/2 - t_07_s));

	reminder_F_T_d_I= reminder_F_T_d_I +  (sign_kz*T_d_std.distr_list[ty]*(xg/2.0)   + B_d_std.distr_list[ty]*xg/2.0 )*(h1+h2)*Exp( -1.0*Eg_off*abs(T/2 - t_07_s));
	reminder_F_T_d_I_15= reminder_F_T_d_I_15 +  (sign_kz*T_d_std.distr_list[ty]*(xg/2.0)   + B_d_std.distr_list[ty]*xg/2.0 )*(h1+h2)*Exp( -1.0*1.5*a_distr*abs(T/2 - t_07_s));
	reminder_F_T_d_I_7= reminder_F_T_d_I_7 +  (sign_kz*T_d_std.distr_list[ty]*(xg/2.0)   + B_d_std.distr_list[ty]*xg/2.0 )*(h1+h2)*Exp( -1.0*0.7*a_distr*abs(T/2 - t_07_s));
	reminder_F_T_d_I_sp= reminder_F_T_d_I_sp +   (sign_kz*T_d.distr_list[ty]*(xg/2.0) + B_d.distr_list[ty]*xg/2.0)*(h1+h2)*Exp( -1.0*Eg_off*abs(T/2 - t_07_s_HLT));

	if(ty> t_07_s) {
	  TO2_reminder = TO2_reminder +  (sign_kz*T_d_std.distr_list[ty]*(xg/2.0)   + B_d_std.distr_list[ty]*xg/2.0 )*(h1+h2)*Exp( -1.0*Eg_off*abs(T/2 - t_07_s));
	}
	if(ty> t_07_s_HLT) {
	  TO2_reminder_sp = TO2_reminder_sp +  (sign_kz*T_d.distr_list[ty]*(xg/2.0)   + B_d.distr_list[ty]*xg/2.0 )*(h1+h2)*Exp( -1.0*Eg_off*abs(T/2 - t_07_s_HLT));
	}
	
      }

      //sum reminders 
      FA_T_d_I_real = FA_T_d_I_real + reminder_FA_T_d_I_real;
      FV_T_d_I_real = FV_T_d_I_real + reminder_FV_T_d_I_real;
      F_T_d_I = F_T_d_I + reminder_F_T_d_I;
      F_T_d_I_15 = F_T_d_I_15 + reminder_F_T_d_I_15;
      F_T_d_I_7 = F_T_d_I_7 + reminder_F_T_d_I_7;
      F_T_d_I_sp = F_T_d_I_sp + reminder_F_T_d_I_sp;


      reminder_F_T_d_I = reminder_F_T_d_I*ZT*(1.0/(mel_SMSM*Eg))*Exp( Eg_off*abs(T/2 - t_07_s))*Exp( M_P*t_07_s) ;
      reminder_F_T_d_I_sp = reminder_F_T_d_I_sp*ZT*(1.0/(mel_SMSM*Eg))*Exp( Eg_off*abs(T/2 - t_07_s_HLT))*Exp( M_P*t_07_s_HLT);

      TO2_reminder = TO2_reminder*ZT*(1.0/(mel_SMSM*Eg))*Exp( Eg_off*abs(T/2 - t_07_s))*Exp( M_P*t_07_s) ;
      TO2_reminder_sp = TO2_reminder_sp*ZT*(1.0/(mel_SMSM*Eg))*Exp( Eg_off*abs(T/2 - t_07_s_HLT))*Exp( M_P*t_07_s_HLT);

      distr_t comb_F_d = ZT*(1.0/(mel_SMSM*Eg))*Exp( Eg_off*abs(T/2 - t_07_s))*Exp( M_P*t_07_s);

      distr_t comb_F_u = ZT*(1.0/(mel_SMSM*Eg))*Exp( Eg_off*abs(T/2 - t_07_c))*Exp( M_P*t_07_c);

      
   
      Print_To_File({}, {((comb_F_d/a_distr)*(sign_kz*T_d_std*(xg/2.0)   + B_d_std*xg/2.0)*Exp(-T*Eg_off/2)*EXPT_D( Eg_off,T_d_std.size())).ave(), ( (comb_F_d/a_distr)*(sign_kz*T_d_std*(xg/2.0)   + B_d_std*xg/2.0)*Exp(-T*Eg_off/2)*EXPT_D( Eg_off,T_d_std.size())).err() }, path_out+"/corr_2pts/CT_d_xg_"+to_string(ixg)+"_"+Ens_tags[iens], "", "");

      Print_To_File({}, {((comb_F_u/a_distr)*(sign_kz*T_u_std*(xg/2.0)   + B_u_std*xg/2.0)*Exp(-T*Eg_off/2)*EXPT_D( Eg_off,T_u_std.size())).ave(), ( (comb_F_u/a_distr)*(sign_kz*T_u_std*(xg/2.0)   + B_u_std*xg/2.0)*Exp(-T*Eg_off/2)*EXPT_D( Eg_off,T_u_std.size())).err() }, path_out+"/corr_2pts/CT_u_xg_"+to_string(ixg)+"_"+Ens_tags[iens], "", "");

      Print_To_File({}, { Vfloat({-0.0001,0,0.0001}), ((comb_F_d/a_distr)*C_0_corr).ave(), ((comb_F_d/a_distr)*C_0_corr).err()}, path_out+"/corr_2pts/C0T_d_xg_"+to_string(ixg)+"_"+Ens_tags[iens], "", "");
      
         
      //normalize
      F_T_u = F_T_u*ZT*(1.0/(mel_SMSM*Eg))*Exp( Eg_off*abs(T/2 - t_07_c))*Exp( M_P*t_07_c) ;
      F_T_d_I = F_T_d_I*ZT*(1.0/(mel_SMSM*Eg))*Exp( Eg_off*abs(T/2 - t_07_s))*Exp( M_P*t_07_s) ;
      F_T_d_I_ori = F_T_d_I_ori*ZT*(1.0/(mel_SMSM*Eg))*Exp( Eg_off*abs(T/2 - t_07_s))*Exp( M_P*t_07_s) ;
      F_T_d_I_ori_sp = F_T_d_I_ori_sp*ZT*(1.0/(mel_SMSM*Eg))*Exp( Eg_off*abs(T/2 - t_07_s_HLT))*Exp( M_P*t_07_s_HLT) ;
      F_T_d_I_15 = F_T_d_I_15*ZT*(1.0/(mel_SMSM*Eg))*Exp( 1.5*a_distr*abs(T/2 - t_07_s))*Exp( M_P*t_07_s) ;
      F_T_d_I_7 = F_T_d_I_7*ZT*(1.0/(mel_SMSM*Eg))*Exp( 0.7*a_distr*abs(T/2 - t_07_s))*Exp( M_P*t_07_s) ;
	  
      FV_T_u_real= FV_T_u_real*ZT*(1.0/(mel_SMSM*Eg))*exp( Eg*abs(T/2 - t_07_c))*Exp( M_P*t_07_c) ;
      FV_T_d_real= FV_T_d_real*ZT*(1.0/(mel_SMSM*Eg))*exp( Eg*abs(T/2 - t_07_s))*Exp( M_P*t_07_s) ;
      FA_T_u_real= FA_T_u_real*ZT*(1.0/(mel_SMSM*Eg))*exp( Eg*abs(T/2 - t_07_c))*Exp( M_P*t_07_c) ;
      FA_T_d_real= FA_T_d_real*ZT*(1.0/(mel_SMSM*Eg))*exp( Eg*abs(T/2 - t_07_s))*Exp( M_P*t_07_s) ;
      FA_T_d_I_real= FA_T_d_I_real*ZT*(1.0/(mel_SMSM*Eg))*exp( Eg*abs(T/2 - t_07_s))*Exp( M_P*t_07_s) ;
      FA_T_d_II_real= FA_T_d_II_real*ZT*(1.0/(mel_SMSM*Eg))*exp( Eg*abs(T/2 - t_07_s))*Exp( M_P*t_07_s) ;

      FV_T_d_I_real = FV_T_d_I_real*ZT*(1.0/(mel_SMSM*Eg))*exp( Eg*abs(T/2 - t_07_s))*Exp( M_P*t_07_s) ;
      FV_T_d_II_real= FV_T_d_II_real*ZT*(1.0/(mel_SMSM*Eg))*exp( Eg*abs(T/2 - t_07_s))*Exp( M_P*t_07_s) ;

      FV_T_d_sp_real= FV_T_d_sp_real*ZT*(1.0/(mel_SMSM*Eg))*exp( Eg*abs(T/2 - t_07_s_HLT))*Exp( M_P*t_07_s_HLT) ;
      FA_T_d_sp_real= FA_T_d_sp_real*ZT*(1.0/(mel_SMSM*Eg))*exp( Eg*abs(T/2 - t_07_s_HLT))*Exp( M_P*t_07_s_HLT) ;
      F_T_d_I_sp = F_T_d_I_sp*ZT*(1.0/(mel_SMSM*Eg))*Exp( Eg_off*abs(T/2 - t_07_s_HLT))*Exp( M_P*t_07_s_HLT);


      F_T_d_ty = F_T_d_ty*ZT*(1.0/(mel_SMSM*Eg))*Exp( Eg_off*abs(T/2 - t_07_s))*Exp( M_P*t_07_s);
      F_T_d_sp_ty = F_T_d_sp_ty*ZT*(1.0/(mel_SMSM*Eg))*Exp( Eg_off*abs(T/2 - t_07_s_HLT))*Exp( M_P*t_07_s_HLT);
      F_T_u_ty= F_T_u_ty*ZT*(1.0/(mel_SMSM*Eg))*Exp( Eg_off*abs(T/2 - t_07_c))*Exp( M_P*t_07_c);

      F_T_d_ty_psum = F_T_d_ty_psum*ZT*(1.0/(mel_SMSM*Eg))*Exp( Eg_off*abs(T/2 - t_07_s))*Exp( M_P*t_07_s);
      F_T_d_sp_ty_psum = F_T_d_sp_ty_psum*ZT*(1.0/(mel_SMSM*Eg))*Exp( Eg_off*abs(T/2 - t_07_s_HLT))*Exp( M_P*t_07_s_HLT);
      F_T_u_ty_psum= F_T_u_ty_psum*ZT*(1.0/(mel_SMSM*Eg))*Exp( Eg_off*abs(T/2 - t_07_c))*Exp( M_P*t_07_c);


      FA_T_d_real_ty = FA_T_d_real_ty*ZT*(1.0/(mel_SMSM*Eg))*exp( Eg*abs(T/2 - t_07_s))*Exp( M_P*t_07_s);
      FA_T_d_sp_real_ty = FA_T_d_sp_real_ty*ZT*(1.0/(mel_SMSM*Eg))*exp( Eg*abs(T/2 - t_07_s_HLT))*Exp( M_P*t_07_s_HLT);
      FA_T_u_real_ty= FA_T_u_real_ty*ZT*(1.0/(mel_SMSM*Eg))*exp( Eg*abs(T/2 - t_07_c))*Exp( M_P*t_07_c);

      FA_T_d_real_ty_psum = FA_T_d_real_ty_psum*ZT*(1.0/(mel_SMSM*Eg))*exp( Eg*abs(T/2 - t_07_s))*Exp( M_P*t_07_s);
      FA_T_d_sp_real_ty_psum = FA_T_d_sp_real_ty_psum*ZT*(1.0/(mel_SMSM*Eg))*exp( Eg*abs(T/2 - t_07_s_HLT))*Exp( M_P*t_07_s_HLT);
      FA_T_u_real_ty_psum= FA_T_u_real_ty_psum*ZT*(1.0/(mel_SMSM*Eg))*exp( Eg*abs(T/2 - t_07_c))*Exp( M_P*t_07_c);

      FV_T_d_real_ty = FV_T_d_real_ty*ZT*(1.0/(mel_SMSM*Eg))*exp( Eg*abs(T/2 - t_07_s))*Exp( M_P*t_07_s);
      FV_T_d_sp_real_ty = FV_T_d_sp_real_ty*ZT*(1.0/(mel_SMSM*Eg))*exp( Eg*abs(T/2 - t_07_s_HLT))*Exp( M_P*t_07_s_HLT);
      FV_T_u_real_ty= FV_T_u_real_ty*ZT*(1.0/(mel_SMSM*Eg))*exp( Eg*abs(T/2 - t_07_c))*Exp( M_P*t_07_c);

      FV_T_d_real_ty_psum = FV_T_d_real_ty_psum*ZT*(1.0/(mel_SMSM*Eg))*exp( Eg*abs(T/2 - t_07_s))*Exp( M_P*t_07_s);
      FV_T_d_sp_real_ty_psum = FV_T_d_sp_real_ty_psum*ZT*(1.0/(mel_SMSM*Eg))*exp( Eg*abs(T/2 - t_07_s_HLT))*Exp( M_P*t_07_s_HLT);
      FV_T_u_real_ty_psum= FV_T_u_real_ty_psum*ZT*(1.0/(mel_SMSM*Eg))*exp( Eg*abs(T/2 - t_07_c))*Exp( M_P*t_07_c);

           
      //push_back
      F_T_u_list[ixg].distr_list.push_back( F_T_u);
      F_T_d_I_list[ixg].distr_list.push_back( F_T_d_I);
      F_T_d_I_ori_list[ixg].distr_list.push_back( F_T_d_I_ori);
      F_T_d_I_ori_sp_list[ixg].distr_list.push_back( F_T_d_I_ori_sp);
      F_T_d_I_list_15[ixg].distr_list.push_back( F_T_d_I_15);
      F_T_d_I_list_7[ixg].distr_list.push_back( F_T_d_I_7);
      F_T_d_sp_I_list[ixg].distr_list.push_back( F_T_d_I_sp);
     
      FV_T_u_real_list[ixg].distr_list.push_back( FV_T_u_real);
      FV_T_d_real_list[ixg].distr_list.push_back( FV_T_d_real);
      FA_T_u_real_list[ixg].distr_list.push_back( FA_T_u_real);
      FA_T_d_real_list[ixg].distr_list.push_back( FA_T_d_real);
      FA_T_d_I_real_list[ixg].distr_list.push_back( FA_T_d_I_real);
      FA_T_d_II_real_list[ixg].distr_list.push_back( FA_T_d_II_real);
      FV_T_d_I_real_list[ixg].distr_list.push_back( FV_T_d_I_real);
      FV_T_d_II_real_list[ixg].distr_list.push_back( FV_T_d_II_real);

      FV_T_d_sp_real_list[ixg].distr_list.push_back( FV_T_d_sp_real);
      FA_T_d_sp_real_list[ixg].distr_list.push_back( FA_T_d_sp_real);

      reminder_F_T_d_list[ixg].distr_list.push_back( reminder_F_T_d_I);
      reminder_F_T_d_sp_list[ixg].distr_list.push_back( reminder_F_T_d_I_sp);

      //print to File ty analysis
      Print_To_File( { }, { F_T_d_ty.ave(), F_T_d_ty.err(), F_T_d_ty_psum.ave(), F_T_d_ty_psum.err()}, path_out+"/FF_d/"+Ens_tags[iens]+"/"+TAG_CURR+"ty_analysis_xg_"+to_string_with_precision(xg_list.ave(ixg),2), "", "");
      Print_To_File( { }, { F_T_d_sp_ty.ave(), F_T_d_sp_ty.err(), F_T_d_sp_ty_psum.ave(), F_T_d_sp_ty_psum.err()}, path_out+"/FF_d/"+Ens_tags[iens]+"/"+TAG_CURR+"ty_analysis_sp_xg_"+to_string_with_precision(xg_list.ave(ixg),2), "", "");
      Print_To_File( { }, { F_T_u_ty.ave(), F_T_u_ty.err(),  F_T_u_ty_psum.ave(), F_T_u_ty_psum.err() }, path_out+"/FF_u/"+Ens_tags[iens]+"/"+TAG_CURR+"ty_analysis_xg_"+to_string_with_precision(xg_list.ave(ixg),2), "", "");

      Print_To_File( { }, { FA_T_d_real_ty.ave(), FA_T_d_real_ty.err(), FA_T_d_real_ty_psum.ave(), FA_T_d_real_ty_psum.err()}, path_out+"/FF_d/"+Ens_tags[iens]+"/"+TAG_CURR+"TA_real_ty_analysis_xg_"+to_string_with_precision(xg_list.ave(ixg),2), "", "");
      Print_To_File( { }, { FA_T_d_sp_real_ty.ave(), FA_T_d_sp_real_ty.err(), FA_T_d_sp_real_ty_psum.ave(), FA_T_d_sp_real_ty_psum.err()}, path_out+"/FF_d/"+Ens_tags[iens]+"/"+TAG_CURR+"TA_real_ty_analysis_sp_xg_"+to_string_with_precision(xg_list.ave(ixg),2), "", "");
      Print_To_File( { }, { FA_T_u_real_ty.ave(), FA_T_u_real_ty.err(),  FA_T_u_real_ty_psum.ave(), FA_T_u_real_ty_psum.err() }, path_out+"/FF_u/"+Ens_tags[iens]+"/"+TAG_CURR+"TA_real_ty_analysis_xg_"+to_string_with_precision(xg_list.ave(ixg),2), "", "");

      Print_To_File( { }, { FV_T_d_real_ty.ave(), FV_T_d_real_ty.err(), FV_T_d_real_ty_psum.ave(), FV_T_d_real_ty_psum.err()}, path_out+"/FF_d/"+Ens_tags[iens]+"/"+TAG_CURR+"TV_real_ty_analysis_xg_"+to_string_with_precision(xg_list.ave(ixg),2), "", "");
      Print_To_File( { }, { FV_T_d_sp_real_ty.ave(), FV_T_d_sp_real_ty.err(), FV_T_d_sp_real_ty_psum.ave(), FV_T_d_sp_real_ty_psum.err()}, path_out+"/FF_d/"+Ens_tags[iens]+"/"+TAG_CURR+"TV_real_ty_analysis_sp_xg_"+to_string_with_precision(xg_list.ave(ixg),2), "", "");
      Print_To_File( { }, { FV_T_u_real_ty.ave(), FV_T_u_real_ty.err(),  FV_T_u_real_ty_psum.ave(), FV_T_u_real_ty_psum.err() }, path_out+"/FF_u/"+Ens_tags[iens]+"/"+TAG_CURR+"TV_real_ty_analysis_xg_"+to_string_with_precision(xg_list.ave(ixg),2), "", "");


      Print_To_File( { }, { (FA_T_d_real_ty+FA_T_u_real_ty).ave(), (FA_T_d_real_ty+FA_T_u_real_ty).err(), (FA_T_d_real_ty_psum+FA_T_u_real_ty_psum).ave(), (FA_T_d_real_ty_psum+FA_T_u_real_ty).err()}, path_out+"/FF/"+Ens_tags[iens]+"/"+TAG_CURR+"TA_real_ty_analysis_xg_"+to_string_with_precision(xg_list.ave(ixg),2), "", "");
      Print_To_File( { }, { (FV_T_d_real_ty+FV_T_u_real_ty).ave(), (FV_T_d_real_ty+FV_T_u_real_ty).err(), (FV_T_d_real_ty_psum+FV_T_u_real_ty_psum).ave(), (FV_T_d_real_ty_psum+FV_T_u_real_ty).err()}, path_out+"/FF/"+Ens_tags[iens]+"/"+TAG_CURR+"TV_real_ty_analysis_xg_"+to_string_with_precision(xg_list.ave(ixg),2), "", "");
        
    
      
   

   
            

      //define corrs to be used in spec dens reconstruction
      distr_t_list Corr_T(UseJack);
      distr_t_list Corr_T_boot(0);

      distr_t_list Corr_TV_real(UseJack), Corr_TV_real_boot(0);

      for(int t=t_07_s_HLT; t <= Corr.Nt/2; t++) {
	
	Corr_T.distr_list.push_back( sign_kz*T_d.distr_list[t] + B_d.distr_list[t]);
	Corr_T_boot.distr_list.push_back( sign_kz*T_d_boot.distr_list[t] + B_d_boot.distr_list[t]);

	Corr_TV_real.distr_list.push_back(   sign_kz*T_d.distr_list[t]*(1 - xg/2.0) + B_d.distr_list[t]*(xg/2.0));
	Corr_TV_real_boot.distr_list.push_back(   sign_kz*T_d_boot.distr_list[t]*(1 - xg_boot/2.0) + B_d_boot.distr_list[t]*(xg_boot/2.0));
      }

  

      cout<<"ZT: "<<ZT_boot.ave()<<" +- "<<ZT_boot.err()<<" [boot], "<<ZT.ave()<<" +- "<<ZT.err()<<" [jack]"<<endl;
      cout<<"mel: "<<mel_SMSM_boot.ave()<<" +- "<<mel_SMSM_boot.err()<<" [boot], "<<mel_SMSM.ave()<<" +- "<<mel_SMSM.err()<<" [jack]"<<endl;
      cout<<"M_P: "<<M_P_boot.ave()<<" +- "<<M_P_boot.err()<<" [boot], "<<M_P.ave()<<" +- "<<M_P.err()<<" [jack]"<<endl;
      cout<<"xg: "<<xg_boot.ave()<<" +- "<<xg_boot.err()<<" [boot], "<<xg.ave()<<" +- "<<xg.err()<<" [jack]"<<endl;
      cout<<"Exp: "<<EXP_D(M_P*t_07_s_HLT).ave()<<" +- "<<EXP_D(M_P*t_07_s_HLT).err()<<" [boot], "<<Exp(M_P*t_07_s_HLT).ave()<<" +- "<<Exp(M_P*t_07_s_HLT).err()<<" [jack]"<<endl;

      
      
      cout<<"Exp(m*t)_boot: "<<EXP_D(M_P_boot*t_07_s_HLT).size()<<endl;

      //rescale error of corr_boot

    
      distr_t FACT_boot= ZT_boot*(1.0/(mel_SMSM_boot*Eg))*EXP_D( M_P_boot*t_07_s_HLT)*xg_boot/2.0;
      distr_t_list Corr_T_boot_to_print=FACT_boot*Corr_T_boot;
      distr_t FACT= ZT*(1.0/(mel_SMSM*Eg) )*Exp( M_P*t_07_s_HLT)*xg/2.0;
      distr_t FACT_real = ZT*(1.0/(mel_SMSM*Eg) )*Exp( M_P*t_07_s_HLT);

      for(int t=t_07_s_HLT; t<= Corr.Nt/2; t++) {

	Corr_T_boot_to_print.distr_list[t-t_07_s_HLT] = Corr_T_boot_to_print.ave(t-t_07_s_HLT) + (Corr_T_boot_to_print.distr_list[t-t_07_s_HLT] - Corr_T_boot_to_print.ave(t-t_07_s_HLT))*(FACT*Corr_T).err(t-t_07_s_HLT)/Corr_T_boot_to_print.err(t-t_07_s_HLT);

	Corr_T_boot.distr_list[t-t_07_s_HLT] =  Corr_T_boot.ave(t-t_07_s_HLT) + (Corr_T_boot.distr_list[t-t_07_s_HLT] - Corr_T_boot.ave(t-t_07_s_HLT))*(Corr_T).err(t-t_07_s_HLT)/Corr_T_boot.err(t-t_07_s_HLT);
	
      }
      
      
      

   

      //print correlator to file
      boost::filesystem::create_directory("../Nazario_data_C_HLT");
      ofstream Print_boot("../Nazario_data_C_HLT/C_mh_"+to_string_with_precision(M_P.ave(), 6)+"_"+to_string_with_precision(xg_list.ave(ixg),2)+"_"+Ens_tags[iens]+".boot");
      Print_boot.precision(10);
      Print_boot<<"#Nconfs=1000     T="<<Corr.Nt/2 - t_07_s_HLT<<endl;
      for(int iboot=0;iboot<1000;iboot++) {
	for(int t=t_07_s_HLT;t<=Corr.Nt/2; t++) Print_boot<<Corr_T_boot_to_print.distr_list[t-t_07_s_HLT].distr[iboot]<<endl ;
      }
      Print_boot.close();

      int tmax= Corr_T.size();

     
      Print_To_File({}, { (Corr_T*FACT).ave(),  (Corr_T*FACT).err(),  (Corr_T_boot_to_print).ave(),  (Corr_T_boot_to_print).err()},  path_out+"/corr_2pts/HLT_boot_jack_xg_"+to_string_with_precision(xg_list.ave(ixg),2), "",""); 

      
      CorrAnalysis Corr_HLT(UseJack, Njacks,1000);
      CorrAnalysis Corr_HLT_boot(0, Njacks, 1000);
      if(MESON=="B0s") {
	if(ixg==0) {
	  if(Ens_tags[iens]=="cB211b.072.64") {
	    Corr_HLT.Tmin= 14;
	    Corr_HLT.Tmax= 25;
	  }
	  else if(Ens_tags[iens]=="cD211a.054.96") {
	    Corr_HLT.Tmin= 22; //27; //22;
	    Corr_HLT.Tmax= 26; //34; //26;
	  }
	}
	else if(ixg==1) {
	  if(Ens_tags[iens]=="cB211b.072.64") {
	    Corr_HLT.Tmin= 14;
	    Corr_HLT.Tmax= 25;
	  }
	  else if(Ens_tags[iens]=="cD211a.054.96") {
	    Corr_HLT.Tmin= 22; //27; //22;
	    Corr_HLT.Tmax= 26; //34; //26;
	  }
	}
	else if(ixg==2) {
	  if(Ens_tags[iens]=="cB211b.072.64") {
	    Corr_HLT.Tmin= 13;
	    Corr_HLT.Tmax= 19;
	  }
	  else if(Ens_tags[iens]=="cD211a.054.96") {
	    Corr_HLT.Tmin= 22; //27; //22;
	    Corr_HLT.Tmax= 26; //34; //26;
	  }
	}
	else if(ixg==3) {
	  if(Ens_tags[iens]=="cB211b.072.64") {
	    Corr_HLT.Tmin= 13;
	    Corr_HLT.Tmax= 18;
	  }
	  else if(Ens_tags[iens]=="cD211a.054.96") {
	    Corr_HLT.Tmin= 23; //27; //23;
	    Corr_HLT.Tmax= 26; //34; //26;
	  }
	  
	}
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
	
      }
      else if(MESON=="B1s") {
	 if(ixg==0) {
	   if(Ens_tags[iens]=="cB211b.072.64") {
	     Corr_HLT.Tmin= 17;
	     Corr_HLT.Tmax= 21;
	   }
	   else if(Ens_tags[iens]=="cD211a.054.96") {
	    Corr_HLT.Tmin= 19;
	    Corr_HLT.Tmax= 26;
	   }
	 }
	 else if(ixg==1) {
	    if(Ens_tags[iens]=="cB211b.072.64") {
	      Corr_HLT.Tmin= 17;
	      Corr_HLT.Tmax= 22;
	    }
	    else if(Ens_tags[iens]=="cD211a.054.96") {
	      Corr_HLT.Tmin= 20;
	      Corr_HLT.Tmax= 26;
	    }
	 }
	 else if(ixg==2) {
	    if(Ens_tags[iens]=="cB211b.072.64") {
	      Corr_HLT.Tmin= 16;
	      Corr_HLT.Tmax= 21;
	    }
	    else if(Ens_tags[iens]=="cD211a.054.96") {
	      Corr_HLT.Tmin= 20;
	      Corr_HLT.Tmax= 26;
	    }
	 }
	 else if(ixg==3) {
	    if(Ens_tags[iens]=="cB211b.072.64") {
	      Corr_HLT.Tmin= 16;
	      Corr_HLT.Tmax= 22;
	    }
	    else if(Ens_tags[iens]=="cD211a.054.96") {
	      Corr_HLT.Tmin= 20;
	      Corr_HLT.Tmax= 26;
	    }
	 }
	 else crash("ixg: "+to_string(ixg)+" not yet implemented");
	 
      }
      else if(MESON=="B3s") {
	
	 if(ixg==0) {
	   if(Ens_tags[iens]=="cB211b.072.64") {
	     Corr_HLT.Tmin= 13;
	     Corr_HLT.Tmax= 18;
	   }
	   else if(Ens_tags[iens]=="cD211a.054.96") {
	     Corr_HLT.Tmin= 19;
	     Corr_HLT.Tmax= 27;
	   }
	 }
	 else if(ixg==1) {
	   if(Ens_tags[iens]=="cB211b.072.64") {
	     Corr_HLT.Tmin= 13;
	     Corr_HLT.Tmax= 18;
	   }
	   else if(Ens_tags[iens]=="cD211a.054.96") {
	     Corr_HLT.Tmin= 15;
	     Corr_HLT.Tmax= 24;
	   }
	 }
	 else if(ixg==2) {
	    if(Ens_tags[iens]=="cB211b.072.64") {
	      Corr_HLT.Tmin= 13;
	      Corr_HLT.Tmax= 18;
	    }
	    else if(Ens_tags[iens]=="cD211a.054.96") {
	      Corr_HLT.Tmin= 15;
	      Corr_HLT.Tmax= 21;
	    }
	 }
	 else if(ixg==3) {
	   if(Ens_tags[iens]=="cB211b.072.64") {
	     Corr_HLT.Tmin= 12;
	     Corr_HLT.Tmax= 18;
	   }
	   else if(Ens_tags[iens]=="cD211a.054.96") {
	     Corr_HLT.Tmin= 15;
	     Corr_HLT.Tmax= 21;
	   }
	 }
	 else crash("ixg: "+to_string(ixg)+" not yet implemented");
	 
      }
      else crash("Meson: "+MESON+" not yet implemented");

      Corr_HLT_boot.Tmin= Corr_HLT.Tmin;
      Corr_HLT_boot.Tmax= Corr_HLT.Tmax;

      
      Corr_HLT.Nt=2*(Corr_T.size()-1);
      Corr_HLT_boot.Nt= Corr_HLT.Nt;
      distr_t_list Corr_T_symm= Corr_T;
      distr_t_list Corr_T_symm_boot= Corr_T_boot;
      for(int t=Corr_HLT.Nt/2 +1;t<Corr_HLT.Nt;t++) {
	Corr_T_symm.distr_list.push_back( Corr_T.distr_list[Corr_HLT.Nt -t]);
	Corr_T_symm_boot.distr_list.push_back( Corr_T_boot.distr_list[Corr_HLT_boot.Nt -t]);
      }
      distr_t_list eff_M_Td_distr= Corr_HLT.effective_mass_t(Corr_T_symm, "" );
      distr_t_list eff_M_Td_boot_distr= Corr_HLT_boot.effective_mass_t( Corr_T_symm_boot, "");
      distr_t eff_M_Td= Corr_HLT.Fit_distr(eff_M_Td_distr)/a_distr;
      distr_t eff_M_Td_boot= Corr_HLT.Fit_distr(eff_M_Td_boot_distr);
      double Mphi_motion= sqrt( Mphi*Mphi + pow(kz/a_distr.ave(),2));
     
      Print_To_File({}, {eff_M_Td_distr.ave(), eff_M_Td_distr.err()}, path_out+"/mass/"+Ens_tags[iens]+"/"+TAG_CURR+"Vs_eff_xg_"+to_string_with_precision(xg_list.ave(ixg),2), "", "");
     
      Print_To_File({}, { (Corr_T*FACT).ave(), (Corr_T*FACT).err() }, path_out+"/corr_2pts/"+TAG_CURR+"HLT_"+Ens_tags[iens]+"_xg_"+to_string_with_precision(xg_list.ave(ixg),2), "", "");

      distr_t M_phi_eff= SQRT_D( POW_D(eff_M_Td*a_distr,2) - pow(kz,2))/a_distr;
      distr_t Mphi_at_mb_half= SQRT_D(M_phi_eff*M_phi_eff + pow(0.5*5.367,2));

      cout<<"eff_M_T: "<<eff_M_Td.ave()<<" +- "<<eff_M_Td.err()<<" expected: "<<Mphi_motion<<endl;
      cout<<"M_phi: "<<M_phi_eff.ave()<<" "<<M_phi_eff.err()<<endl;
      distr_t_list phi_MT_distr= Corr_T_symm*EXPT_DL(eff_M_Td_distr);
      distr_t phi_MT = Corr_HLT.Fit_distr( Corr_T_symm*(EXPT_DL(eff_M_Td_distr)));
      distr_t phi_MT_boot= Corr_HLT_boot.Fit_distr( Corr_T_symm_boot*(EXPT_DL(eff_M_Td_boot_distr)));


      distr_t_list fake_distr_with_hole= Get_id_distr_list(Corr_T.size(), Njacks, UseJack);
      fake_distr_with_hole.distr_list[0] = 0.0*fake_distr_with_hole.distr_list[0];

      distr_t_list fake_boot_distr_with_hole= Get_id_distr_list(Corr_T.size(), 1000, 0);
      fake_boot_distr_with_hole.distr_list[0] = 0.0*fake_boot_distr_with_hole.distr_list[0];
      
      distr_t_list Corr_T_VMD = fake_distr_with_hole*phi_MT*EXPT_D(-1.0*eff_M_Td*a_distr, Corr_T.size());
      distr_t_list Corr_T_VMD_boot = fake_boot_distr_with_hole*phi_MT_boot*EXPT_D(-1.0*eff_M_Td_boot, Corr_T_boot.size());

     
      for(int t=1; t < Corr_T_VMD_boot.size(); t++) Corr_T_VMD_boot.distr_list[t] = Corr_T_VMD_boot.ave(t) + (Corr_T_VMD_boot.distr_list[t] - Corr_T_VMD_boot.ave(t))*Corr_T_VMD.err(t)/Corr_T_VMD_boot.err(t);
     
           
      distr_t_list Corr_T_sub = Corr_T - Corr_T_VMD;
      distr_t_list Corr_T_boot_sub = Corr_T_boot - Corr_T_VMD_boot;


      cout<<"PRINTING VMD PREDICTION FOR THE COUPLING ofr meson: "<<MESON+" xg: "<<xg_list.ave(ixg)<<" Ensemble: "<<Ens_tags[iens]<<endl;
      cout<<"g+ * f_V : "<<(FACT*phi_MT/(a_distr)).ave()<<" +- " <<(FACT*phi_MT/(a_distr)).err()<<endl;
      cout<<"Form factor at q^2=q'^2 =0: "<<(FACT*phi_MT*(eff_M_Td/Mphi_at_mb_half)*1.0/(  Mphi_at_mb_half*a_distr -  0.5*5.367*a_distr)).ave()<<endl;
      cout<<"##############################################################################"<<endl;


      //#########################################################################################################



      //Determine second-state 
      
      //ari-symmetrize
      distr_t_list Corr_T_ary_symm= Corr_T_sub;
      distr_t_list Corr_T_ary_symm_boot= Corr_T_boot_sub;
      for(int t=Corr_HLT.Nt/2 +1;t<Corr_HLT.Nt;t++) {
	Corr_T_ary_symm.distr_list.push_back( Corr_T_sub.distr_list[Corr_HLT.Nt -t]);
	Corr_T_ary_symm_boot.distr_list.push_back( Corr_T_boot_sub.distr_list[Corr_HLT.Nt -t]);
      }
      distr_t_list eff_M_prime_Td_distr= Corr_HLT.effective_mass_t( Corr_T_ary_symm, path_out+"/mass/"+Ens_tags[iens]+"/"+TAG_CURR+"Vs_exc_eff_xg_"+to_string_with_precision(xg_list.ave(ixg),2));
      distr_t_list eff_M_prime_Td_boot_distr= Corr_HLT.effective_mass_t( Corr_T_ary_symm_boot, "");
      int Tmin_old= Corr_HLT.Tmin; int Tmax_old= Corr_HLT.Tmax;

      
      if(MESON=="B0s") {
	if(ixg==0) {
	  if(Ens_tags[iens]=="cB211b.072.64") {
	    Corr_HLT.Tmin= 5;
	    Corr_HLT.Tmax= 9;
	  }
	  else if(Ens_tags[iens]=="cD211a.054.96") {
	    Corr_HLT.Tmin=9;
	    Corr_HLT.Tmax=14;
	  }
	}
	else if(ixg==1) {
	  if(Ens_tags[iens]=="cB211b.072.64") {
	    Corr_HLT.Tmin= 5;
	    Corr_HLT.Tmax= 9;
	  }
	  else if(Ens_tags[iens]=="cD211a.054.96") {
	    Corr_HLT.Tmin=10;
	    Corr_HLT.Tmax=14;
	  }
	}
	else if(ixg==2) {
	  if(Ens_tags[iens]=="cB211b.072.64") {
	    Corr_HLT.Tmin= 4;
	    Corr_HLT.Tmax= 6;
	  }
	  else if(Ens_tags[iens]=="cD211a.054.96") {
	    Corr_HLT.Tmin=10;
	    Corr_HLT.Tmax=14;
	  }

	}
	else if(ixg==3) {
	  if(Ens_tags[iens]=="cB211b.072.64") {
	    Corr_HLT.Tmin= 4;
	    Corr_HLT.Tmax= 6;
	  }
	  else if(Ens_tags[iens]=="cD211a.054.96") {
	    Corr_HLT.Tmin=10;
	    Corr_HLT.Tmax=14;
	  }
	}
	else crash("ixg: "+to_string(ixg)+" not yet implemented");

      }
      else if(MESON=="B1s") {
	 
	 if(ixg==0) {
	   if(Ens_tags[iens]=="cB211b.072.64") {
	     Corr_HLT.Tmin= 6;
	     Corr_HLT.Tmax= 9;
	   }
	   else if(Ens_tags[iens]=="cD211a.054.96") {
	     Corr_HLT.Tmin=11;
	     Corr_HLT.Tmax=19;
	   }
	 }
	 else if(ixg==1) {
	   if(Ens_tags[iens]=="cB211b.072.64") {
	     Corr_HLT.Tmin= 6;
	     Corr_HLT.Tmax= 9;
	   }
	   else if(Ens_tags[iens]=="cD211a.054.96") {
	     Corr_HLT.Tmin=11;
	     Corr_HLT.Tmax=17;
	   }
	 }
	 else if(ixg==2) {
	   if(Ens_tags[iens]=="cB211b.072.64") {
	     Corr_HLT.Tmin= 8;
	     Corr_HLT.Tmax= 11;
	   }
	   else if(Ens_tags[iens]=="cD211a.054.96") {
	     Corr_HLT.Tmin=10;
	     Corr_HLT.Tmax=17;
	   }
	 }
	 else if(ixg==3) {
	   if(Ens_tags[iens]=="cB211b.072.64") {
	     Corr_HLT.Tmin= 8;
	     Corr_HLT.Tmax= 11;
	   }
	   else if(Ens_tags[iens]=="cD211a.054.96") {
	     Corr_HLT.Tmin=10;
	     Corr_HLT.Tmax=17;
	   }
	 }
	 else crash("ixg: "+to_string(ixg)+" not yet implemented");
	 
      }
      else if(MESON=="B3s") {
	
	 if(ixg==0) {
	   if(Ens_tags[iens]=="cB211b.072.64") {
	     Corr_HLT.Tmin= 5;
	     Corr_HLT.Tmax= 8;
	   }
	   else if(Ens_tags[iens]=="cD211a.054.96") {
	     Corr_HLT.Tmin=10;
	     Corr_HLT.Tmax=16;
	   }
	 }
	 else if(ixg==1) {
	   if(Ens_tags[iens]=="cB211b.072.64") {
	     Corr_HLT.Tmin= 5;
	     Corr_HLT.Tmax= 7;
	   }
	   else if(Ens_tags[iens]=="cD211a.054.96") {
	     Corr_HLT.Tmin=7;
	     Corr_HLT.Tmax=11;
	   }
	 }
	 else if(ixg==2) {
	   if(Ens_tags[iens]=="cB211b.072.64") {
	     Corr_HLT.Tmin= 5;
	     Corr_HLT.Tmax= 7;
	   }
	   else if(Ens_tags[iens]=="cD211a.054.96") {
	     Corr_HLT.Tmin=4;
	     Corr_HLT.Tmax=8;
	   }
	 }
	 else if(ixg==3) {
	   if(Ens_tags[iens]=="cB211b.072.64") {
	     Corr_HLT.Tmin= 5;
	     Corr_HLT.Tmax= 7;
	   }
	   else if(Ens_tags[iens]=="cD211a.054.96") {
	     Corr_HLT.Tmin=4;
	     Corr_HLT.Tmax=8;
	   }
	 }
	 else crash("ixg: "+to_string(ixg)+" not yet implemented");
	 
      }
      else crash("Meson: "+MESON+" not yet implemented");

      Corr_HLT_boot.Tmin= Corr_HLT.Tmin;
      Corr_HLT_boot.Tmax= Corr_HLT.Tmax;

      
      
      distr_t eff_M_prime_Td= Corr_HLT.Fit_distr(eff_M_prime_Td_distr)/a_distr;
      distr_t eff_M_prime_Td_boot= Corr_HLT_boot.Fit_distr(eff_M_prime_Td_boot_distr);
      distr_t_list phi_prime_MT_distr=  Corr_T_ary_symm*(EXPT_DL(eff_M_prime_Td_distr));
      distr_t_list phi_prime_MT_boot_distr=  Corr_T_ary_symm_boot*(EXPT_DL(eff_M_prime_Td_boot_distr));
      distr_t phi_prime_MT= Corr_HLT.Fit_distr( Corr_T_ary_symm*(EXPT_DL(eff_M_prime_Td_distr)));
      distr_t phi_prime_MT_boot= Corr_HLT_boot.Fit_distr( Corr_T_ary_symm_boot*(EXPT_DL(eff_M_prime_Td_boot_distr)));
      
      Corr_HLT.Tmin = Tmin_old; Corr_HLT.Tmax= Tmax_old;
      Corr_HLT_boot.Tmin = Tmin_old; Corr_HLT_boot.Tmax= Tmax_old;
      
      distr_t_list Corr_T_VMD_II = fake_distr_with_hole*phi_prime_MT*EXPT_D(-1.0*eff_M_prime_Td*a_distr, Corr_T.size());
      distr_t_list Corr_T_sub_II = Corr_T_sub - Corr_T_VMD_II;
      distr_t_list Corr_T_boot_VMD_II = fake_boot_distr_with_hole*phi_prime_MT_boot*EXPT_D(-1.0*eff_M_prime_Td_boot, Corr_T.size());
    
      for(int t=1; t < Corr_T_boot_VMD_II.size(); t++) Corr_T_boot_VMD_II.distr_list[t] = Corr_T_boot_VMD_II.ave(t) + (Corr_T_boot_VMD_II.distr_list[t] - Corr_T_boot_VMD_II.ave(t))*Corr_T_VMD_II.err(t)/Corr_T_boot_VMD_II.err(t);
    
      
      distr_t_list Corr_T_boot_sub_II = Corr_T_boot_sub - Corr_T_boot_VMD_II;



      //#########################################################################################################



      //Determine third-state 
      
      //ari-symmetrize
      distr_t_list Corr_T_aryary_symm= Corr_T_sub_II;
      distr_t_list Corr_T_aryary_symm_boot= Corr_T_boot_sub_II;
      for(int t=Corr_HLT.Nt/2 +1;t<Corr_HLT.Nt;t++) {
	Corr_T_aryary_symm.distr_list.push_back( Corr_T_sub_II.distr_list[Corr_HLT.Nt -t]);
	Corr_T_aryary_symm_boot.distr_list.push_back( Corr_T_boot_sub_II.distr_list[Corr_HLT.Nt -t]);
      }
      distr_t_list eff_M_second_Td_distr= Corr_HLT.effective_mass_t( Corr_T_aryary_symm, path_out+"/mass/"+Ens_tags[iens]+"/"+TAG_CURR+"Vs_exc2_eff_xg_"+to_string_with_precision(xg_list.ave(ixg),2));
      distr_t_list eff_M_second_Td_boot_distr= Corr_HLT.effective_mass_t( Corr_T_aryary_symm_boot, "");
      Tmin_old= Corr_HLT.Tmin; Tmax_old= Corr_HLT.Tmax;
      if(MESON=="B0s") {
	if(ixg==0) {
	  if(Ens_tags[iens]=="cB211b.072.64") {
	    Corr_HLT.Tmin= 1;
	    Corr_HLT.Tmax= 3;
	  }
	  else if(Ens_tags[iens]=="cD211a.054.96") {
	    Corr_HLT.Tmin=2;
	    Corr_HLT.Tmax=4;
	  }
	}
	else if(ixg==1) {
	  if(Ens_tags[iens]=="cB211b.072.64") {
	    Corr_HLT.Tmin= 1;
	    Corr_HLT.Tmax= 3;
	  }
	  else if(Ens_tags[iens]=="cD211a.054.96") {
	    Corr_HLT.Tmin=2;
	    Corr_HLT.Tmax=4;
	  }
	}
	else if(ixg==2) {
	  if(Ens_tags[iens]=="cB211b.072.64") {
	    Corr_HLT.Tmin= 1;
	    Corr_HLT.Tmax= 3;
	  }
	  else if(Ens_tags[iens]=="cD211a.054.96") {
	    Corr_HLT.Tmin=2;
	    Corr_HLT.Tmax=4;
	  }
	}
	else if(ixg==3) {
	  if(Ens_tags[iens]=="cB211b.072.64") {
	    Corr_HLT.Tmin= 1;
	    Corr_HLT.Tmax= 3;
	  }
	  else if(Ens_tags[iens]=="cD211a.054.96") {
	    Corr_HLT.Tmin=2;
	    Corr_HLT.Tmax=4;
	  }
	}
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
	
      }
      else if(MESON=="B1s") {
	 
	 if(ixg==0) {
	    if(Ens_tags[iens]=="cB211b.072.64") {
	      Corr_HLT.Tmin= 1;
	      Corr_HLT.Tmax= 3;
	    }
	    else if(Ens_tags[iens]=="cD211a.054.96") {
	      Corr_HLT.Tmin=6;
	      Corr_HLT.Tmax=9;
	    }
	 }
	 else if(ixg==1) {
	   if(Ens_tags[iens]=="cB211b.072.64") {
	     Corr_HLT.Tmin= 1;
	     Corr_HLT.Tmax= 3;
	   }
	   else if(Ens_tags[iens]=="cD211a.054.96") {
	     Corr_HLT.Tmin=4;
	     Corr_HLT.Tmax=8;
	   }
	 }
	 else if(ixg==2) {
	   if(Ens_tags[iens]=="cB211b.072.64") {
	     Corr_HLT.Tmin= 1;
	     Corr_HLT.Tmax= 3;
	   }
	   else if(Ens_tags[iens]=="cD211a.054.96") {
	     Corr_HLT.Tmin=4;
	     Corr_HLT.Tmax=8;
	   }
	 }
	 else if(ixg==3) {
	   if(Ens_tags[iens]=="cB211b.072.64") {
	     Corr_HLT.Tmin= 1;
	     Corr_HLT.Tmax= 3;
	   }
	   else if(Ens_tags[iens]=="cD211a.054.96") {
	     Corr_HLT.Tmin=4;
	     Corr_HLT.Tmax=7;
	   }
	 }
	 else crash("ixg: "+to_string(ixg)+" not yet implemented");
	 
      }
      else if(MESON=="B3s") {
	
	 if(ixg==0) {
	   if(Ens_tags[iens]=="cB211b.072.64") {
	     Corr_HLT.Tmin= 1;
	     Corr_HLT.Tmax= 3;
	   }
	   else if(Ens_tags[iens]=="cD211a.054.96") {
	     Corr_HLT.Tmin=2;
	     Corr_HLT.Tmax=6;
	   }
	 }
	 else if(ixg==1) {
	   if(Ens_tags[iens]=="cB211b.072.64") {
	     Corr_HLT.Tmin= 1;
	     Corr_HLT.Tmax= 3;
	   }
	   else if(Ens_tags[iens]=="cD211a.054.96") {
	     Corr_HLT.Tmin=1;
	     Corr_HLT.Tmax=4;
	   }
	 }
	 else if(ixg==2) {
	   if(Ens_tags[iens]=="cB211b.072.64") {
	     Corr_HLT.Tmin= 1;
	     Corr_HLT.Tmax= 3;
	   }
	   else if(Ens_tags[iens]=="cD211a.054.96") {
	     Corr_HLT.Tmin=1;
	     Corr_HLT.Tmax=3;
	   }
	 }
	 else if(ixg==3) {
	   if(Ens_tags[iens]=="cB211b.072.64") {
	     Corr_HLT.Tmin= 1;
	     Corr_HLT.Tmax= 3;
	   }
	   else if(Ens_tags[iens]=="cD211a.054.96") {
	     Corr_HLT.Tmin=1;
	     Corr_HLT.Tmax=3;
	   }
	 }
	 else crash("ixg: "+to_string(ixg)+" not yet implemented");
	 
      }
      else crash("Meson: "+MESON+" not yet implemented");

      Corr_HLT_boot.Tmin= Corr_HLT.Tmin;
      Corr_HLT_boot.Tmax= Corr_HLT.Tmax;

      
      
      distr_t eff_M_second_Td= Corr_HLT.Fit_distr(eff_M_second_Td_distr)/a_distr;
      distr_t eff_M_second_Td_boot= Corr_HLT_boot.Fit_distr(eff_M_second_Td_boot_distr);
      distr_t_list phi_second_MT_distr=  Corr_T_aryary_symm*(EXPT_DL(eff_M_second_Td_distr));
      distr_t_list phi_second_MT_boot_distr=  Corr_T_aryary_symm_boot*(EXPT_DL(eff_M_second_Td_boot_distr));
      distr_t phi_second_MT= Corr_HLT.Fit_distr( Corr_T_aryary_symm*(EXPT_DL(eff_M_second_Td_distr)));
      distr_t phi_second_MT_boot= Corr_HLT_boot.Fit_distr( Corr_T_aryary_symm_boot*(EXPT_DL(eff_M_second_Td_boot_distr)));
      
      Corr_HLT.Tmin = Tmin_old; Corr_HLT.Tmax= Tmax_old;
      Corr_HLT_boot.Tmin = Tmin_old; Corr_HLT_boot.Tmax= Tmax_old;

    
      
      distr_t_list Corr_T_VMD_III = fake_distr_with_hole*phi_second_MT*EXPT_D(-1.0*eff_M_second_Td*a_distr, Corr_T.size());
      distr_t_list Corr_T_sub_III = Corr_T_sub_II - Corr_T_VMD_III;
      distr_t_list Corr_T_boot_VMD_III = fake_boot_distr_with_hole*phi_second_MT_boot*EXPT_D(-1.0*eff_M_second_Td_boot, Corr_T.size());
  
      for(int t=1; t < Corr_T_boot_VMD_III.size(); t++) Corr_T_boot_VMD_III.distr_list[t] = Corr_T_boot_VMD_III.ave(t) + (Corr_T_boot_VMD_III.distr_list[t] - Corr_T_boot_VMD_III.ave(t))*Corr_T_VMD_III.err(t)/Corr_T_boot_VMD_III.err(t);
  
      distr_t_list Corr_T_boot_sub_III = Corr_T_boot_sub_II - Corr_T_boot_VMD_III;





      //#########################################################################################################

      
      Print_To_File({}, { (Corr_T_sub*FACT).ave(), (Corr_T_sub*FACT).err(), (Corr_T_sub_II*FACT).ave(), (Corr_T_sub_II*FACT).err(), (Corr_T_sub_III*FACT).ave(), (Corr_T_sub_II*FACT).err() }, path_out+"/corr_2pts/"+TAG_CURR+"sub_HLT_"+Ens_tags[iens]+"_xg_"+to_string_with_precision(xg_list.ave(ixg),2), "", "");

      distr_t_list Corr_T_ori= Corr_T;

      if(Use_preconditioning) {
	Corr_T= Corr_T_sub;
	Corr_T_boot= Corr_T_boot_sub;
      }

      Print_To_File({}, { (FACT*phi_MT_distr/a_distr).ave(), (FACT*phi_MT_distr/a_distr).err(), (FACT*phi_prime_MT_distr/a_distr).ave(), (FACT*phi_prime_MT_distr/a_distr).err(), (FACT*phi_second_MT_distr/a_distr).ave(), (FACT*phi_second_MT_distr/a_distr).err()}, path_out+"/FF_d_II/"+Ens_tags[iens]+"/VMD_MT_xg_"+to_string_with_precision(xg_list.ave(ixg),2), "", "");
 
      //Generate Covariance matrix
      Vfloat cov_T, corr_T;
      Vfloat TT, RR;
      for(int tt=0;tt< tmax;tt++)
	for(int rr=0;rr< tmax;rr++) {
	  TT.push_back(tt);
	  RR.push_back(rr);
	  
	  corr_T.push_back( ( Corr_T_boot.distr_list[tt]%Corr_T_boot.distr_list[rr])/( Corr_T_boot.err(tt)*Corr_T_boot.err(rr)));
	  cov_T.push_back( (Corr_T_boot.distr_list[tt]%Corr_T_boot.distr_list[rr])*(Corr_T.err(tt)*Corr_T.err(rr)/(Corr_T_boot.err(tt)*Corr_T_boot.err(rr))));
	  
	  
	}

      int tmax_new=0;
      //tmax_new = tmax-1;
      for(int t=1;t<tmax;t++) if( Corr_T_ori.err(t)/fabs(Corr_T_ori.ave(t)) < 0.4) tmax_new++;
      int tmax_new_im=0;
      for(int t=1;t<tmax;t++) if( Corr_T_ori.err(t)/fabs(Corr_T_ori.ave(t)) < 0.4) tmax_new_im++;
   
   

      //print covariance matrix
      Print_To_File({},{TT,RR, cov_T, corr_T}, path_out+"/covariance/"+Ens_tags[iens]+"/"+preco_tag+"cov_T_xg_"+to_string_with_precision(xg.ave(),2)+".cov", "" , "");
  

      cout<<"Starting spectral reconstruction:..."<<endl;

      for(int iss=0; iss <(signed)sigmas_07_w0.size(); iss++) {



	string SM_FF="FF_Exp";
	//#############################################################

	   auto K_RE_distr_sub= [&SM_FF](const distr_t &E, const distr_t &m, double s) -> distr_t {
	     
	     distr_t ret(1);
	     if(SM_FF=="FF_Exp")  {
	       for(int ijack=0;ijack<E.size();ijack++) {
		 double x= (E.distr[ijack]-m.distr[ijack]);
		 double y= E.distr[ijack];
		 ret.distr.push_back( 2*sinh(x)*cos(s)/(cosh(2*x) - cos(2*s)) - 2*sinh(y)*cos(0)/(cosh(2*y) - cos(2*0))  );
	       }
	     }
	     
	     
	     else if(SM_FF=="FF_Gauss") {
	       for(int ijack=0;ijack<E.size();ijack++) {
		 double x= (E.distr[ijack]-m.distr[ijack]);
		 double y= E.distr[ijack];
		 ret.distr.push_back( 2*sinh(x)*cos(s)/(cosh(2*x) - cos(2*s)) - 2*sinh(y)*cos(0)/(cosh(2*y) - cos(2*0))  );
	       }
	     }
	     else crash("SM_FF: "+SM_FF+" not yet implemented");
	     
	     return ret;
	     
	   };

	   auto K_IM_distr= [&SM_FF](const distr_t &E, const distr_t &m, double s) -> distr_t {
	     
	     
	     distr_t ret(1);
	     if(SM_FF=="FF_Exp") {
	       for(int ijack=0;ijack<E.size();ijack++) {
		 double x= (E.distr[ijack]-m.distr[ijack]);
		 ret.distr.push_back ( 2*cosh(x)*sin(s)/(cosh(2*x) - cos(2*s)) );
		 //ret.distr.push_back(exp(-x)*sin(s)/( 1 + exp(-2*x) -2*cos(s)*exp(-x)));
	       }
	     }
	     else if(SM_FF=="FF_Gauss") {
	       for(int ijack=0;ijack<E.size();ijack++) {
		 ret.distr.push_back( M_PI*Get_exact_gauss( E.distr[ijack], m.distr[ijack],  s, 0.0 ));
	       }
	     }
	     else crash("SM_FF: "+SM_FF+" not yet implemented");
	     
	     
	     return ret;
	     
	   };







	//#############################################################


	F_T_d_MB_RE_VMD_sm_list[ixg][iss].distr_list.push_back( FACT*phi_MT*K_RE_distr_sub( eff_M_Td*a_distr, Eg_MB_off, sigmas_07_w0[iss]*a_distr.ave()));
	F_T_d_MB_IM_VMD_sm_list[ixg][iss].distr_list.push_back( FACT*phi_MT*K_IM_distr( eff_M_Td*a_distr, Eg_MB_off,  sigmas_07_w0[iss]*a_distr.ave()));

	F_T_d_MB_RE_VMD_II_state_sm_list[ixg][iss].distr_list.push_back( FACT*phi_MT*K_RE_distr_sub( eff_M_Td*a_distr, Eg_MB_off, sigmas_07_w0[iss]*a_distr.ave())   +  FACT*phi_prime_MT*K_RE_distr_sub( eff_M_prime_Td*a_distr, Eg_MB_off, sigmas_07_w0[iss]*a_distr.ave())   );
	F_T_d_MB_IM_VMD_II_state_sm_list[ixg][iss].distr_list.push_back( FACT*phi_MT*K_IM_distr( eff_M_Td*a_distr, Eg_MB_off,  sigmas_07_w0[iss]*a_distr.ave())  + FACT*phi_prime_MT*K_IM_distr( eff_M_prime_Td*a_distr, Eg_MB_off, sigmas_07_w0[iss]*a_distr.ave())    );

	F_T_d_MB_RE_VMD_III_state_sm_list[ixg][iss].distr_list.push_back( FACT*phi_MT*K_RE_distr_sub( eff_M_Td*a_distr, Eg_MB_off, sigmas_07_w0[iss]*a_distr.ave())   + FACT*phi_prime_MT*K_RE_distr_sub( eff_M_prime_Td*a_distr, Eg_MB_off, sigmas_07_w0[iss]*a_distr.ave()) +  FACT*phi_second_MT*K_RE_distr_sub( eff_M_second_Td*a_distr, Eg_MB_off, sigmas_07_w0[iss]*a_distr.ave())  );
	F_T_d_MB_IM_VMD_III_state_sm_list[ixg][iss].distr_list.push_back( FACT*phi_MT*K_IM_distr( eff_M_Td*a_distr, Eg_MB_off,  sigmas_07_w0[iss]*a_distr.ave())  + FACT*phi_prime_MT*K_IM_distr( eff_M_prime_Td*a_distr, Eg_MB_off, sigmas_07_w0[iss]*a_distr.ave()) +  FACT*phi_second_MT*K_IM_distr( eff_M_second_Td*a_distr, Eg_MB_off, sigmas_07_w0[iss]*a_distr.ave())    );

	distr_t_list VMD_VIRT_SCAN_RE(UseJack), VMD_VIRT_SCAN_IM(UseJack);
	distr_t_list VMD_II_VIRT_SCAN_RE(UseJack), VMD_II_VIRT_SCAN_IM(UseJack);
	distr_t_list VMD_III_VIRT_SCAN_RE(UseJack), VMD_III_VIRT_SCAN_IM(UseJack);

	for(int vir=0;vir<(signed)virtualities.size(); vir++) {

	  VMD_VIRT_SCAN_RE.distr_list.push_back(  FACT*phi_MT*K_RE_distr_sub( eff_M_Td*a_distr , M_P*virtualities[vir] , sigmas_07_w0[iss]*a_distr.ave()));
	  VMD_VIRT_SCAN_IM.distr_list.push_back( FACT*phi_MT*K_IM_distr( eff_M_Td*a_distr, M_P*virtualities[vir],  sigmas_07_w0[iss]*a_distr.ave()) );

	  VMD_II_VIRT_SCAN_RE.distr_list.push_back(  FACT*phi_MT*K_RE_distr_sub( eff_M_Td*a_distr , M_P*virtualities[vir] , sigmas_07_w0[iss]*a_distr.ave())   +  FACT*phi_prime_MT*K_RE_distr_sub( eff_M_prime_Td*a_distr, M_P*virtualities[vir], sigmas_07_w0[iss]*a_distr.ave())   );
	  VMD_II_VIRT_SCAN_IM.distr_list.push_back(  FACT*phi_MT*K_IM_distr( eff_M_Td*a_distr , M_P*virtualities[vir],  sigmas_07_w0[iss]*a_distr.ave())   +  FACT*phi_prime_MT*K_IM_distr( eff_M_prime_Td*a_distr, M_P*virtualities[vir], sigmas_07_w0[iss]*a_distr.ave())    );

	  VMD_III_VIRT_SCAN_RE.distr_list.push_back(  FACT*phi_MT*K_RE_distr_sub( eff_M_Td*a_distr , M_P*virtualities[vir] , sigmas_07_w0[iss]*a_distr.ave())   +  FACT*phi_prime_MT*K_RE_distr_sub( eff_M_prime_Td*a_distr, M_P*virtualities[vir], sigmas_07_w0[iss]*a_distr.ave())   +  FACT*phi_second_MT*K_RE_distr_sub( eff_M_second_Td*a_distr, M_P*virtualities[vir], sigmas_07_w0[iss]*a_distr.ave())   );
	  VMD_III_VIRT_SCAN_IM.distr_list.push_back( FACT*phi_MT*K_IM_distr( eff_M_Td*a_distr, M_P*virtualities[vir],  sigmas_07_w0[iss]*a_distr.ave())  + FACT*phi_prime_MT*K_IM_distr( eff_M_prime_Td*a_distr, M_P*virtualities[vir], sigmas_07_w0[iss]*a_distr.ave())  + FACT*phi_second_MT*K_IM_distr( eff_M_second_Td*a_distr, M_P*virtualities[vir], sigmas_07_w0[iss]*a_distr.ave())    );

	}

	Print_To_File({}, {(virtualities*Get_id_jack_distr_list(virtualities.size(), Njacks)*M_P/a_distr).ave(), VMD_VIRT_SCAN_RE.ave(), VMD_VIRT_SCAN_RE.err(), VMD_VIRT_SCAN_IM.ave(), VMD_VIRT_SCAN_IM.err()}, path_out+"/FF_d_II/"+Ens_tags[iens]+"/VMD_virt_scan/xg_"+to_string_with_precision(xg_list.ave(ixg),2)+"_sm_"+to_string_with_precision(sigmas_07_w0[iss],3)+"_"+Ens_tags[iens], "", "");

	Print_To_File({}, {(virtualities*Get_id_jack_distr_list(virtualities.size(), Njacks)*M_P/a_distr).ave(), VMD_II_VIRT_SCAN_RE.ave(), VMD_II_VIRT_SCAN_RE.err(), VMD_II_VIRT_SCAN_IM.ave(), VMD_II_VIRT_SCAN_IM.err()}, path_out+"/FF_d_II/"+Ens_tags[iens]+"/VMD_virt_scan/II_xg_"+to_string_with_precision(xg_list.ave(ixg),2)+"_sm_"+to_string_with_precision(sigmas_07_w0[iss],3)+"_"+Ens_tags[iens], "", "");

	Print_To_File({}, {(virtualities*Get_id_jack_distr_list(virtualities.size(), Njacks)*M_P/a_distr).ave(), VMD_III_VIRT_SCAN_RE.ave(), VMD_III_VIRT_SCAN_RE.err(), VMD_III_VIRT_SCAN_IM.ave(), VMD_III_VIRT_SCAN_IM.err()}, path_out+"/FF_d_II/"+Ens_tags[iens]+"/VMD_virt_scan/III_xg_"+to_string_with_precision(xg_list.ave(ixg),2)+"_sm_"+to_string_with_precision(sigmas_07_w0[iss],3)+"_"+Ens_tags[iens], "", "");
      }

      for(int isg=0;isg<(signed)sigmas_07.size();isg++) {  if(iens==0) sigma_simulated[ixg][isg] = sigmas_07[isg]*Eg_MB_off_new.ave()/MBs; }

      string SM_FF_0 = "FF_Exp";
      
       auto K_RE_0= [&SM_FF_0](const PrecFloat &E, const PrecFloat &m, const PrecFloat &s, const PrecFloat &E0, int ijack) -> PrecFloat {
    

	     if(SM_FF_0=="FF_Gauss") {
	       PrecFloat x= (E-m);
	       
	       PrecFloat cosh_ov_sinh_half= (exp(x) + exp(-3*x))/(1-exp(-2*x)); 
      
	       return  2*cos(s)/(cosh_ov_sinh_half - cos(2*s)/sinh(x));
      
	     }


	     
	     if(SM_FF_0=="FF_Exp") {
      
	       PrecFloat x= (E-m);
	       
	       PrecFloat cosh_ov_sinh_half= (exp(x) + exp(-3*x))/(1-exp(-2*x)); 
	       
	       return  2*cos(s)/(cosh_ov_sinh_half - cos(2*s)/sinh(x));
	       
	      
	     }
         
	     exit(-1);
	     return 0;
       };


         auto K_RE_0_distr= [&SM_FF_0](const distr_t &E, const distr_t &m, double s) -> distr_t {
	     
	     distr_t ret(1);
	     if(SM_FF_0=="FF_Exp")  {
	       for(int ijack=0;ijack<E.size();ijack++) {
		 double x= (E.distr[ijack]-m.distr[ijack]);
		 double y= E.distr[ijack];
		 ret.distr.push_back(2*sinh(y)*cos(0)/(cosh(2*y) - cos(2*0))  );
	       }
	     }
	     
	     
	     else if(SM_FF_0=="FF_Gauss") {
	       for(int ijack=0;ijack<E.size();ijack++) {
		 double x= (E.distr[ijack]-m.distr[ijack]);
		 double y= E.distr[ijack];
		 ret.distr.push_back(2*sinh(y)*cos(0)/(cosh(2*y) - cos(2*0))  );
	       }
	     }
	     else crash("SM_FF: "+SM_FF_0+" not yet implemented");
	     
	     return ret;
	     
	   };

      //0 REMINDER FROM HLT
       double th_0= E0_fact*Mphi_motion;
       double l_re_T_0, syst_T_0;
       distr_t REMINDER_MB_from_HLT=  Get_Laplace_transfo_tmin( 0.0,  0.0, th_0*a_distr.ave(),  Nts[iens], 1,  tmax_new, prec_07, SM_FF_0+"_RE_ixg_"+to_string(ixg+1),K_RE_0, Corr_T, syst_T_0, 1e-5 ,  l_re_T_0, MODE_FF, "E0_"+to_string_with_precision(E0_fact,1), TAG_CURR+"MB_T_"+Ens_tags[iens], 1e-6,0, FACT, Get_id_jack_distr(Njacks)*0.0,  "zero_energy_"+preco_tag+MESON+"_07_FF_Tw_"+to_string(t_07_s_HLT), cov_T, fake_func,0, fake_func_d ,  1 , 4.0, 0.00,1);

       distr_t REMINDER_MB_PHI_from_HLT=  ((Use_preconditioning==false)?0.0:1.0)*(FACT*phi_MT*K_RE_0_distr( eff_M_Td*a_distr, 0.0*Eg_MB_off, 0.0)) ;

       distr_t REMINDER_MB(UseJack,Njacks);
       
       for(int t=1;t<Corr_T.size();t++) REMINDER_MB = REMINDER_MB + FACT*Corr_T[t];

       distr_t REMINDER_MB_TOT(UseJack,Njacks);
       
       for(int t=1;t<Corr_T_ori.size();t++) REMINDER_MB_TOT = REMINDER_MB_TOT + FACT*Corr_T_ori[t];

       cout<<"REM MB: "<<REMINDER_MB.ave()<<" +- "<<REMINDER_MB.err()<<" REM MB(HLT): "<<REMINDER_MB_from_HLT.ave()<<" +- "<<REMINDER_MB_from_HLT.err()<<endl;
       cout<<"REM MB(inc phi:): "<<(REMINDER_MB+REMINDER_MB_PHI_from_HLT).ave()<<" +- "<<(REMINDER_MB+REMINDER_MB_PHI_from_HLT).err()<<" REM MB(inc phi: HLT): "<<(REMINDER_MB_from_HLT+REMINDER_MB_PHI_from_HLT).ave()<<" +- "<<(REMINDER_MB_from_HLT+REMINDER_MB_PHI_from_HLT).err()<<endl;
       cout<<"REM MB(full sum:): "<<(REMINDER_MB_TOT).ave()<<" +- "<<(REMINDER_MB_TOT).err()<<" REM MB(inc phi: HLT): "<<(REMINDER_MB_from_HLT+REMINDER_MB_PHI_from_HLT).ave()<<" +- "<<(REMINDER_MB_from_HLT+REMINDER_MB_PHI_from_HLT).err()<<endl;

       cout<<"REM TOT: "<<reminder_F_T_d_I.ave()<<" +- "<<reminder_F_T_d_I.err()<<" SP: "<<reminder_F_T_d_I_sp.ave()<<" +- "<<reminder_F_T_d_I_sp.err()<<endl;
       cout<<"REM 2TO: "<<TO2_reminder.ave()<<" +- "<<TO2_reminder.err()<<" SP: "<<TO2_reminder_sp.ave()<<" +- "<<TO2_reminder_sp.err()<<endl;

       REMINDER_MB= REMINDER_MB_from_HLT;

      if(!Skip_spectral_reconstruction_07) {

	vector<distr_t_list> F_T_d_RE_virt_scan;
	vector<distr_t_list> F_T_d_IM_virt_scan;

	for(int isg=0;isg<(signed)sigmas_07.size();isg++) {
	  F_T_d_RE_virt_scan.emplace_back(UseJack, virtualities.size(), Njacks);
	  F_T_d_IM_virt_scan.emplace_back(UseJack, virtualities.size(), Njacks);
	}

	#pragma omp parallel for schedule(dynamic)
	//spectral reconstruction for second time ordering
	for(int isg=0;isg<(signed)sigmas_07.size();isg++) {



	  //#########################################################################



	   string SM_FF = "FF_Exp";



	   auto K_RE_distr= [&SM_FF](const distr_t &E, const distr_t &m, double s) -> distr_t {
	     
	     distr_t ret(1);
	     if(SM_FF=="FF_Exp")  {
	       for(int ijack=0;ijack<E.size();ijack++) {
		 double x= (E.distr[ijack]-m.distr[ijack]);
		 ret.distr.push_back( 2*sinh(x)*cos(s)/(cosh(2*x) - cos(2*s)));
	       }
	     }
	     
	     
	     else if(SM_FF=="FF_Gauss") {
	       for(int ijack=0;ijack<E.size();ijack++) {
		 double x= (E.distr[ijack]-m.distr[ijack]);
		 ret.distr.push_back(  2*sinh(x)*cos(s)/(cosh(2*x) - cos(2*s)));
	       }
	     }
	     else crash("SM_FF: "+SM_FF+" not yet implemented");
	     
	     return ret;
	     
	   };
	   
	   
	   auto K_RE_distr_sub= [&SM_FF](const distr_t &E, const distr_t &m, double s) -> distr_t {
	     
	     distr_t ret(1);
	     if(SM_FF=="FF_Exp")  {
	       for(int ijack=0;ijack<E.size();ijack++) {
		 double x= (E.distr[ijack]-m.distr[ijack]);
		 double y= E.distr[ijack];
		 ret.distr.push_back( 2*sinh(x)*cos(s)/(cosh(2*x) - cos(2*s)) - 2*sinh(y)*cos(0)/(cosh(2*y) - cos(2*0))  );
	       }
	     }
	     
	     
	     else if(SM_FF=="FF_Gauss") {
	       for(int ijack=0;ijack<E.size();ijack++) {
		 double x= (E.distr[ijack]-m.distr[ijack]);
		 double y= E.distr[ijack];
		 ret.distr.push_back( 2*sinh(x)*cos(s)/(cosh(2*x) - cos(2*s)) - 2*sinh(y)*cos(0)/(cosh(2*y) - cos(2*0))  );
	       }
	     }
	     else crash("SM_FF: "+SM_FF+" not yet implemented");
	     
	     return ret;
	     
	   };

	   auto K_IM_distr= [&SM_FF](const distr_t &E, const distr_t &m, double s) -> distr_t {
	     
	     
	     distr_t ret(1);
	     if(SM_FF=="FF_Exp") {
	       for(int ijack=0;ijack<E.size();ijack++) {
		 double x= (E.distr[ijack]-m.distr[ijack]);
		 ret.distr.push_back ( 2*cosh(x)*sin(s)/(cosh(2*x) - cos(2*s)) );
		 //ret.distr.push_back(exp(-x)*sin(s)/( 1 + exp(-2*x) -2*cos(s)*exp(-x)));
	       }
	     }
	     else if(SM_FF=="FF_Gauss") {
	       for(int ijack=0;ijack<E.size();ijack++) {
		 ret.distr.push_back( M_PI*Get_exact_gauss( E.distr[ijack], m.distr[ijack],  s, 0.0 ));
	       }
	     }
	     else crash("SM_FF: "+SM_FF+" not yet implemented");
	     
	     
	     return ret;
	     
	   };
	   
 
	   
	   auto K_RE= [&SM_FF](const PrecFloat &E, const PrecFloat &m, const PrecFloat &s, const PrecFloat &E0, int ijack) -> PrecFloat {
    

	     if(SM_FF=="FF_Gauss") {
	       PrecFloat x= (E-m);
	       
	       PrecFloat cosh_ov_sinh_half= (exp(x) + exp(-3*x))/(1-exp(-2*x)); 
      
	       return  2*cos(s)/(cosh_ov_sinh_half - cos(2*s)/sinh(x));
      
	     }


	     
	     if(SM_FF=="FF_Exp") {
      
	       PrecFloat x= (E-m);
	       
	       PrecFloat cosh_ov_sinh_half= (exp(x) + exp(-3*x))/(1-exp(-2*x)); 
	       
	       return  2*cos(s)/(cosh_ov_sinh_half - cos(2*s)/sinh(x));
	       
	      
	     }
         
	     exit(-1);
	     return 0;
	   };
	   

	   auto K_IM = [&SM_FF](const PrecFloat &E, const PrecFloat &m, const PrecFloat &s, const PrecFloat &E0, int ijack) -> PrecFloat {

	     if(SM_FF=="FF_Gauss") return precPi()*Get_exact_gauss(E, m, s, E0);
	     
    
	     
	     else if(SM_FF=="FF_Exp") {
	       
	       PrecFloat x= (E-m);

	       PrecFloat cosh_ov_cosh_half= (exp(x) + exp(-3*x))/(1+exp(-2*x)); 

	       return 2*sin(s)/(cosh_ov_cosh_half - cos(2*s)/cosh(x));
	       
	     }

	     return PrecFloat(0);
    
	   };

	  //###########################################################################








	  

	  double s= sigmas_07[isg]*a_distr.ave();

	  double syst_T;
	  double mult_T_IM= 0.5;

	 
	  double mult_T_RE=0.05;
	  double Ag_target=1e-2;
	  if(sigmas_07[isg] < 0.5) Ag_target=1e-2;
	  else if(sigmas_07[isg] < 1.5) Ag_target=1e-2;
	  double th= E0_fact*Mphi_motion;
	  double l_re_T;
	  mult_T_RE=1;
	  
	  
	
	  if(SM_FF=="FF_Exp") {
	  //select multiplicities
	    if(MESON=="B0s" || MESON=="B1s") mult_T_IM=0.1;
	    if(MESON=="B3s") mult_T_IM=0.3;
	    //if(MESON=="B3s" && ixg == 0) mult_T_RE=0.005;
	    //if(MESON=="B3s" && ixg == 1) mult_T_RE=0.01;
	    if(MESON=="B3s" && ixg >= 2 && Ens_tags[iens] == "cB211b.072.64") mult_T_RE=0.05;
	    if(MESON=="B0s" && ixg==3 && Ens_tags[iens] == "cD211a.054.96") mult_T_RE=10;
	    if(MESON=="B1s" && ixg >= 2 ) mult_T_RE=6;

	    if( sigmas_07[isg] < 2.1 && MESON == "B1s"  ) mult_T_IM=2.5;
	    else if( sigmas_07[isg] > 2 && MESON=="B1s") mult_T_IM=1;

	    
	    if( MESON == "B3s") mult_T_IM=4;

	    if( (MESON =="B3s") && (sigmas_07[isg] > 3.9) && (Ens_tags[iens] == "cB211b.072.64") ) mult_T_RE=10;


	    

	    if( MESON=="B3s" && sigmas_07[isg] > 3.9) {
	      if(ixg==0) mult_T_IM=1.0;
	      if(ixg==1) mult_T_IM=0.4;
	      if(ixg==2) mult_T_IM=0.1;
	      if(ixg==3) mult_T_IM=0.1;
	      
	    }

	    //Vfloat sigmas_07({1.75, 2, 2.25, 2.5, 2.75, 3.0, 3.5, 4.0, 4.5});  

	    if(MESON=="B3s" && ixg==0) {
	      if(sigmas_07[isg] < 1.8) mult_T_IM=2.5;
	      if(sigmas_07[isg] > 1.99 && sigmas_07[isg] < 2.28) mult_T_IM=1.0;
	      if(sigmas_07[isg] > 2.4 && sigmas_07[isg] < 4.6) mult_T_IM=0.4;
	      
	    }
	    
	    if(MESON=="B3s" && ixg==1) {
	      if(sigmas_07[isg] > 2.74 || sigmas_07[isg] < 3.1) mult_T_IM=1.0;
	      if(sigmas_07[isg] > 3.49) mult_T_IM=0.1;

	    }

	  }
	  
	 
    

	  distr_t RE_sm= Get_id_jack_distr(Njacks);
	  distr_t IM_sm= Get_id_jack_distr(Njacks);
	  

	  cout<<"#########################  E= MB - Egamma ######################### "<<endl<<flush;

	  //modify s in order to make it scale with the energy oE

	  s= s*Eg_MB_off.ave()/(MBs*a_distr.ave());
	 

	

	
	  distr_t OFFSET_RE=  ((Use_preconditioning==false)?0.0:1.0)*(FACT*phi_MT*K_RE_distr_sub( eff_M_Td*a_distr, Eg_MB_off, s) + 0.0*FACT*phi_prime_MT*K_RE_distr_sub( eff_M_prime_Td*a_distr, Eg_MB_off, s) ) + F_T_d_I -REMINDER_MB;


	  
	  distr_t RE_MB_sm = Get_Laplace_transfo_tmin( Eg_MB_off.ave(),  s, th*a_distr.ave(),  Nts[iens], 1, tmax_new, prec_07, SM_FF+"_RE_ixg_"+to_string(ixg+1),K_RE, Corr_T, syst_T, mult_T_RE ,  l_re_T, MODE_FF, "E0_"+to_string_with_precision(E0_fact,1), TAG_CURR+"MB_T_"+Ens_tags[iens], Ag_target,0, FACT, OFFSET_RE,  preco_tag+MESON+"_07_FF_Tw_"+to_string(t_07_s_HLT), cov_T, fake_func,0, fake_func_d ,  1 , 4.0, 0.00,1);
	  
	  RE_MB_sm = RE_MB_sm.ave() +  (RE_MB_sm - RE_MB_sm.ave())*(sqrt( pow(syst_T,2) + pow(RE_MB_sm.err(),2)))/RE_MB_sm.err() -REMINDER_MB;


	  distr_t OFFSET_RE_15=  ((Use_preconditioning==false)?0.0:1.0)*(FACT*phi_MT*K_RE_distr_sub( eff_M_Td*a_distr, 1.5*a_distr, s) + 0.0*FACT*phi_prime_MT*K_RE_distr_sub( eff_M_prime_Td*a_distr, Eg_MB_off, s) ) + F_T_d_I_15 -REMINDER_MB;

	  
	  distr_t RE_MB_sm_15 = Get_Laplace_transfo_tmin( (1.5*a_distr).ave(),  s, th*a_distr.ave(),  Nts[iens], 1, tmax_new, prec_07, SM_FF+"_RE_15_ixg_"+to_string(ixg+1),K_RE, Corr_T, syst_T, mult_T_RE ,  l_re_T, MODE_FF, "E0_"+to_string_with_precision(E0_fact,1), TAG_CURR+"MB_T_"+Ens_tags[iens], Ag_target,0, FACT, OFFSET_RE,  preco_tag+MESON+"_07_FF_Tw_"+to_string(t_07_s_HLT), cov_T, fake_func,0, fake_func_d ,  1 , 4.0, 0.00,1);
	  
	  RE_MB_sm_15 = RE_MB_sm_15.ave() +  (RE_MB_sm_15 - RE_MB_sm_15.ave())*(sqrt( pow(syst_T,2) + pow(RE_MB_sm_15.err(),2)))/RE_MB_sm_15.err() -REMINDER_MB;



	  distr_t OFFSET_RE_7=  ((Use_preconditioning==false)?0.0:1.0)*(FACT*phi_MT*K_RE_distr_sub( eff_M_Td*a_distr, 0.7*a_distr, s) + 0.0*FACT*phi_prime_MT*K_RE_distr_sub( eff_M_prime_Td*a_distr, Eg_MB_off, s) ) + F_T_d_I_7 -REMINDER_MB;

	  
	  distr_t RE_MB_sm_7 = Get_Laplace_transfo_tmin( (0.7*a_distr).ave(),  s, th*a_distr.ave(),  Nts[iens], 1, tmax_new, prec_07, SM_FF+"_RE_7_ixg_"+to_string(ixg+1),K_RE, Corr_T, syst_T, mult_T_RE ,  l_re_T, MODE_FF, "E0_"+to_string_with_precision(E0_fact,1), TAG_CURR+"MB_T_"+Ens_tags[iens], Ag_target,0, FACT, OFFSET_RE,  preco_tag+MESON+"_07_FF_Tw_"+to_string(t_07_s_HLT), cov_T, fake_func,0, fake_func_d ,  1 , 4.0, 0.00,1);
	  
	  RE_MB_sm_7 = RE_MB_sm_7.ave() +  (RE_MB_sm_7 - RE_MB_sm_7.ave())*(sqrt( pow(syst_T,2) + pow(RE_MB_sm_7.err(),2)))/RE_MB_sm_7.err() -REMINDER_MB;

	  
	  distr_t OFFSET_IM=   ((Use_preconditioning==false)?0.0:1.0)*(FACT*phi_MT*K_IM_distr( eff_M_Td*a_distr, Eg_MB_off, s)+ 0.0*FACT*phi_prime_MT*K_IM_distr( eff_M_prime_Td*a_distr, Eg_MB_off, s) );

	  
	  distr_t IM_MB_sm = Get_Laplace_transfo_tmin( Eg_MB_off.ave(),  s, th*a_distr.ave(),  Nts[iens], 1, tmax_new_im, prec_07, SM_FF+"_IM_ixg_"+to_string(ixg+1),K_IM, Corr_T, syst_T, mult_T_IM ,  l_re_T, MODE_FF, "E0_"+to_string_with_precision(E0_fact,1), TAG_CURR+"MB_T_"+Ens_tags[iens], Ag_target,0, FACT, OFFSET_IM , preco_tag+MESON+"_07_FF_Tw_"+to_string(t_07_s_HLT), cov_T, fake_func,0, fake_func_d ,  1 , 4.0, 0.0,1);
	  
	  IM_MB_sm = IM_MB_sm.ave() +  (IM_MB_sm - IM_MB_sm.ave())*(sqrt( pow(syst_T,2) + pow(IM_MB_sm.err(),2)))/IM_MB_sm.err();


	  distr_t OFFSET_IM_15=   ((Use_preconditioning==false)?0.0:1.0)*(FACT*phi_MT*K_IM_distr( eff_M_Td*a_distr, 1.5*a_distr, s)+ 0.0*FACT*phi_prime_MT*K_IM_distr( eff_M_prime_Td*a_distr, Eg_MB_off, s) );
	  
	  distr_t IM_MB_sm_15 = Get_Laplace_transfo_tmin( (1.5*a_distr).ave(),  s, th*a_distr.ave(),  Nts[iens], 1, tmax_new_im, prec_07, SM_FF+"_IM_15_ixg_"+to_string(ixg+1),K_IM, Corr_T, syst_T, mult_T_IM ,  l_re_T, MODE_FF, "E0_"+to_string_with_precision(E0_fact,1), TAG_CURR+"MB_T_"+Ens_tags[iens], Ag_target,0, FACT, OFFSET_IM , preco_tag+MESON+"_07_FF_Tw_"+to_string(t_07_s_HLT), cov_T, fake_func,0, fake_func_d ,  1 , 4.0, 0.0,1);
	  
	  IM_MB_sm_15 = IM_MB_sm_15.ave() +  (IM_MB_sm_15 - IM_MB_sm_15.ave())*(sqrt( pow(syst_T,2) + pow(IM_MB_sm_15.err(),2)))/IM_MB_sm_15.err();


	  distr_t OFFSET_IM_7=   ((Use_preconditioning==false)?0.0:1.0)*(FACT*phi_MT*K_IM_distr( eff_M_Td*a_distr, 0.7*a_distr, s)+ 0.0*FACT*phi_prime_MT*K_IM_distr( eff_M_prime_Td*a_distr, Eg_MB_off, s) );

	  distr_t IM_MB_sm_7 = Get_Laplace_transfo_tmin( (0.7*a_distr).ave(),  s, th*a_distr.ave(),  Nts[iens], 1, tmax_new_im, prec_07, SM_FF+"_IM_7_ixg_"+to_string(ixg+1),K_IM, Corr_T, syst_T, mult_T_IM ,  l_re_T, MODE_FF, "E0_"+to_string_with_precision(E0_fact,1), TAG_CURR+"MB_T_"+Ens_tags[iens], Ag_target,0, FACT, OFFSET_IM , preco_tag+MESON+"_07_FF_Tw_"+to_string(t_07_s_HLT), cov_T, fake_func,0, fake_func_d ,  1 , 4.0, 0.0,1);
	  
	  IM_MB_sm_7 = IM_MB_sm_7.ave() +  (IM_MB_sm_7 - IM_MB_sm_7.ave())*(sqrt( pow(syst_T,2) + pow(IM_MB_sm_7.err(),2)))/IM_MB_sm_7.err();
	  
	  
	  F_T_d_MB_RE_sm_list[ixg][isg].distr_list.push_back( RE_MB_sm + ((Use_preconditioning==false)?0.0:1.0)*(FACT*phi_MT*K_RE_distr_sub( eff_M_Td*a_distr, Eg_MB_off, s) + 0.0*FACT*phi_prime_MT*K_RE_distr_sub( eff_M_prime_Td*a_distr, Eg_MB_off, s) ));

	  F_T_d_MB_RE_sm_list_no_sub[ixg][isg].distr_list.push_back( RE_MB_sm + REMINDER_MB+  ((Use_preconditioning==false)?0.0:1.0)*(FACT*phi_MT*K_RE_distr( eff_M_Td*a_distr, Eg_MB_off, s) + 0.0*FACT*phi_prime_MT*K_RE_distr( eff_M_prime_Td*a_distr, Eg_MB_off, s) ));
	  
	  F_T_d_MB_IM_sm_list[ixg][isg].distr_list.push_back( IM_MB_sm +  ((Use_preconditioning==false)?0.0:1.0)*(FACT*phi_MT*K_IM_distr( eff_M_Td*a_distr, Eg_MB_off, s)+ 0.0*FACT*phi_prime_MT*K_IM_distr( eff_M_prime_Td*a_distr, Eg_MB_off, s) ));


	  F_T_d_MB_RE_sm_list_15[ixg][isg].distr_list.push_back( RE_MB_sm_15 + ((Use_preconditioning==false)?0.0:1.0)*(FACT*phi_MT*K_RE_distr_sub( eff_M_Td*a_distr, 1.5*a_distr, s) + 0.0*FACT*phi_prime_MT*K_RE_distr_sub( eff_M_prime_Td*a_distr, Eg_MB_off, s) ));
	  
	  F_T_d_MB_IM_sm_list_15[ixg][isg].distr_list.push_back( IM_MB_sm_15 +  ((Use_preconditioning==false)?0.0:1.0)*(FACT*phi_MT*K_IM_distr( eff_M_Td*a_distr, 1.5*a_distr, s)+ 0.0*FACT*phi_prime_MT*K_IM_distr( eff_M_prime_Td*a_distr, Eg_MB_off, s) ));

	  F_T_d_MB_RE_sm_list_7[ixg][isg].distr_list.push_back( RE_MB_sm_7 + ((Use_preconditioning==false)?0.0:1.0)*(FACT*phi_MT*K_RE_distr_sub( eff_M_Td*a_distr, 0.7*a_distr, s) + 0.0*FACT*phi_prime_MT*K_RE_distr_sub( eff_M_prime_Td*a_distr, Eg_MB_off, s) ));
	  
	  F_T_d_MB_IM_sm_list_7[ixg][isg].distr_list.push_back( IM_MB_sm_7 +  ((Use_preconditioning==false)?0.0:1.0)*(FACT*phi_MT*K_IM_distr( eff_M_Td*a_distr, 0.7*a_distr, s)+ 0.0*FACT*phi_prime_MT*K_IM_distr( eff_M_prime_Td*a_distr, Eg_MB_off, s) ));




	  //set multiplicity

	  //sigmas_07({  1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3.0, 3.5});
	  
	  /*
	  
	  SM_FF = "FF_Gauss";

	  mult_T_IM=0.3;
	  if(MESON=="B0s") {

	    if(ixg<3 && (isg==(signed)sigmas_07.size()-1)) { mult_T_IM=0.08;}

	    if(ixg==0 && (isg==(signed)sigmas_07.size()-1) && Ens_tags[iens] == "cB211b.072.64") { mult_T_IM=0.04;}
	    if(ixg==0 && (isg==(signed)sigmas_07.size()-1) && Ens_tags[iens] == "cD211a.054.96") { mult_T_IM=0.2;}

	    if(ixg==0 && (isg==(signed)sigmas_07.size()-2) && Ens_tags[iens] == "cB211b.072.64") { mult_T_IM=0.2;}

	    if(ixg==0 && (isg==(signed)sigmas_07.size()-3) && Ens_tags[iens] == "cB211b.072.64") { mult_T_IM=0.15;}


	    if(ixg==0 && (isg==(signed)sigmas_07.size()-7) && Ens_tags[iens] == "cB211b.072.64") { mult_T_IM=0.11;}
	    if(ixg==0 && (isg==(signed)sigmas_07.size()-7) && Ens_tags[iens] == "cD211a.054.96") { mult_T_IM=0.5;}

	    if(ixg==0 && (isg==(signed)sigmas_07.size()-8) && Ens_tags[iens] == "cB211b.072.64") { mult_T_IM=0.11;}
	    if(ixg==0 && (isg==(signed)sigmas_07.size()-8) && Ens_tags[iens] == "cD211a.054.96") { mult_T_IM=0.5;}

	    if(ixg==0 && (isg==(signed)sigmas_07.size()-9) && Ens_tags[iens] == "cB211b.072.64") { mult_T_IM=0.11;}
	    if(ixg==0 && (isg==(signed)sigmas_07.size()-9) && Ens_tags[iens] == "cD211a.054.96") { mult_T_IM=0.5;}



	    if(ixg==1 && (isg==(signed)sigmas_07.size()-2) && Ens_tags[iens] == "cD211a.054.96") { mult_T_IM=1;}

	    if(ixg==1 && (isg==(signed)sigmas_07.size()-3) && Ens_tags[iens] == "cD211a.054.96") { mult_T_IM=1.5;}

	    if(ixg==1 && (isg==(signed)sigmas_07.size()-4) && Ens_tags[iens] == "cB211b.072.64") { mult_T_IM=0.5;}
	    if(ixg==1 && (isg==(signed)sigmas_07.size()-4) && Ens_tags[iens] == "cD211a.054.96") { mult_T_IM=1;}

	    if(ixg==1 && (isg==(signed)sigmas_07.size()-5) && Ens_tags[iens] == "cB211b.072.64") { mult_T_IM=1;}
	    if(ixg==1 && (isg==(signed)sigmas_07.size()-5) && Ens_tags[iens] == "cD211a.054.96") { mult_T_IM=0.1;}

	    if(ixg==1 && (isg==(signed)sigmas_07.size()-6) && Ens_tags[iens] == "cB211b.072.64") { mult_T_IM=0.3;}
	    if(ixg==1 && (isg==(signed)sigmas_07.size()-6) && Ens_tags[iens] == "cD211a.054.96") { mult_T_IM=0.06;}

	    if(ixg==1 && (isg==(signed)sigmas_07.size()-7) && Ens_tags[iens] == "cB211b.072.64") { mult_T_IM=0.3;}
	    if(ixg==1 && (isg==(signed)sigmas_07.size()-7) && Ens_tags[iens] == "cD211a.054.96") { mult_T_IM=0.3;}

	    if(ixg==1 && (isg==(signed)sigmas_07.size()-8) && Ens_tags[iens] == "cB211b.072.64") { mult_T_IM=0.3;}
	    if(ixg==1 && (isg==(signed)sigmas_07.size()-8) && Ens_tags[iens] == "cD211a.054.96") { mult_T_IM=0.3;}

	    if(ixg==1 && (isg==(signed)sigmas_07.size()-9) && Ens_tags[iens] == "cB211b.072.64") { mult_T_IM=0.3;}
	    if(ixg==1 && (isg==(signed)sigmas_07.size()-9) && Ens_tags[iens] == "cD211a.054.96") { mult_T_IM=0.3;}


	    if(ixg==2 && (isg==(signed)sigmas_07.size()-2) && Ens_tags[iens] == "cB211b.072.64") { mult_T_IM=0.5;}
	    if(ixg==2 && (isg==(signed)sigmas_07.size()-2) && Ens_tags[iens] == "cD211a.054.96") { mult_T_IM=1;}

	    if(ixg==2 && (isg==(signed)sigmas_07.size()-3) && Ens_tags[iens] == "cB211b.072.64") { mult_T_IM=2;}
	    if(ixg==2 && (isg==(signed)sigmas_07.size()-3) && Ens_tags[iens] == "cD211a.054.96") { mult_T_IM=2;}

	    if(ixg==2 && (isg==(signed)sigmas_07.size()-4) && Ens_tags[iens] == "cB211b.072.64") { mult_T_IM=6;}
	    if(ixg==2 && (isg==(signed)sigmas_07.size()-4) && Ens_tags[iens] == "cD211a.054.96") { mult_T_IM=0.5;}

	    if(ixg==2 && (isg==(signed)sigmas_07.size()-5) && Ens_tags[iens] == "cB211b.072.64") { mult_T_IM=0.1;}
	    if(ixg==2 && (isg==(signed)sigmas_07.size()-5) && Ens_tags[iens] == "cD211a.054.96") { mult_T_IM=0.25;}

	    if(ixg==2 && (isg==(signed)sigmas_07.size()-6) && Ens_tags[iens] == "cB211b.072.64") { mult_T_IM=0.03;}
	    if(ixg==2 && (isg==(signed)sigmas_07.size()-6) && Ens_tags[iens] == "cD211a.054.96") { mult_T_IM=0.03;}

	    if(ixg==2 && (isg==(signed)sigmas_07.size()-7) && Ens_tags[iens] == "cB211b.072.64") { mult_T_IM=0.05;}
	    if(ixg==2 && (isg==(signed)sigmas_07.size()-7) && Ens_tags[iens] == "cD211a.054.96") { mult_T_IM=0.05;}

	    if(ixg==2 && (isg==(signed)sigmas_07.size()-8) && Ens_tags[iens] == "cB211b.072.64") { mult_T_IM=0.05;}
	    if(ixg==2 && (isg==(signed)sigmas_07.size()-8) && Ens_tags[iens] == "cD211a.054.96") { mult_T_IM=0.05;}
	    
	    if(ixg==2 && (isg==(signed)sigmas_07.size()-9) && Ens_tags[iens] == "cB211b.072.64") { mult_T_IM=0.05;}
	    if(ixg==2 && (isg==(signed)sigmas_07.size()-9) && Ens_tags[iens] == "cD211a.054.96") { mult_T_IM=0.05;}



	    if(ixg==3 && (isg==(signed)sigmas_07.size()-1) && Ens_tags[iens] == "cB211b.072.64") { mult_T_IM=0.3;}
	    if(ixg==3 && (isg==(signed)sigmas_07.size()-1) && Ens_tags[iens] == "cD211a.054.96") { mult_T_IM=3;}


	    if(ixg==3 && (isg==(signed)sigmas_07.size()-2) && Ens_tags[iens] == "cB211b.072.64") { mult_T_IM=3;}
	    if(ixg==3 && (isg==(signed)sigmas_07.size()-2) && Ens_tags[iens] == "cD211a.054.96") { mult_T_IM=1;}


	    if(ixg==3 && (isg==(signed)sigmas_07.size()-3) && Ens_tags[iens] == "cB211b.072.64") { mult_T_IM=3;}
	    if(ixg==3 && (isg==(signed)sigmas_07.size()-3) && Ens_tags[iens] == "cD211a.054.96") { mult_T_IM=0.013;}


	    if(ixg==3 && (isg==(signed)sigmas_07.size()-4) && Ens_tags[iens] == "cB211b.072.64") { mult_T_IM=2;}
	    if(ixg==3 && (isg==(signed)sigmas_07.size()-4) && Ens_tags[iens] == "cD211a.054.96") { mult_T_IM=0.013;}

	      if(ixg==3 && (isg==(signed)sigmas_07.size()-5) && Ens_tags[iens] == "cB211b.072.64") { mult_T_IM=0.2;}
	    if(ixg==3 && (isg==(signed)sigmas_07.size()-5) && Ens_tags[iens] == "cD211a.054.96") { mult_T_IM=0.013;}

	      if(ixg==3 && (isg==(signed)sigmas_07.size()-6) && Ens_tags[iens] == "cB211b.072.64") { mult_T_IM=0.01;}
	    if(ixg==3 && (isg==(signed)sigmas_07.size()-6) && Ens_tags[iens] == "cD211a.054.96") { mult_T_IM=0.1;}

	      if(ixg==3 && (isg==(signed)sigmas_07.size()-7) && Ens_tags[iens] == "cB211b.072.64") { mult_T_IM=0.01;}
	    if(ixg==3 && (isg==(signed)sigmas_07.size()-7) && Ens_tags[iens] == "cD211a.054.96") { mult_T_IM=0.1;}

	      if(ixg==3 && (isg==(signed)sigmas_07.size()-8) && Ens_tags[iens] == "cB211b.072.64") { mult_T_IM=0.01;}
	    if(ixg==3 && (isg==(signed)sigmas_07.size()-8) && Ens_tags[iens] == "cD211a.054.96") { mult_T_IM=0.1;}

	      if(ixg==3 && (isg==(signed)sigmas_07.size()-9) && Ens_tags[iens] == "cB211b.072.64") { mult_T_IM=0.01;}
	    if(ixg==3 && (isg==(signed)sigmas_07.size()-9) && Ens_tags[iens] == "cD211a.054.96") { mult_T_IM=0.1;}



	    
	    
	    
	   
	  }
	  if(MESON=="B1s") {
	    if(ixg==0 && (isg==(signed)sigmas_07.size()-3) && Ens_tags[iens] == "cB211b.072.64") { mult_T_IM=0.04;}
	    if(ixg==0 && (isg==(signed)sigmas_07.size()-3) && Ens_tags[iens] == "cD211a.054.96") { mult_T_IM=1;}

	    if(ixg==0 && (isg==(signed)sigmas_07.size()-4) && Ens_tags[iens] == "cB211b.072.64") { mult_T_IM=0.04;}
	    if(ixg==0 && (isg==(signed)sigmas_07.size()-4) && Ens_tags[iens] == "cD211a.054.96") { mult_T_IM=1;}
	    
	    if(ixg==0 && (isg==(signed)sigmas_07.size()-5) && Ens_tags[iens] == "cB211b.072.64") { mult_T_IM=0.05;}
	    if(ixg==0 && (isg==(signed)sigmas_07.size()-5) && Ens_tags[iens] == "cD211a.054.96") { mult_T_IM=1;}

	    if(ixg==0 && (isg==(signed)sigmas_07.size()-6) && Ens_tags[iens] == "cD211a.054.96") { mult_T_IM=1;}
	    if(ixg==0 && (isg==(signed)sigmas_07.size()-7) && Ens_tags[iens] == "cD211a.054.96") { mult_T_IM=1;}

	    if(ixg==0 && (isg==(signed)sigmas_07.size()-8) && Ens_tags[iens] == "cB211b.072.64") { mult_T_IM=0.8;}
	    if(ixg==0 && (isg==(signed)sigmas_07.size()-8) && Ens_tags[iens] == "cD211a.054.96") { mult_T_IM=0.8;}

	    if(ixg==0 && (isg==(signed)sigmas_07.size()-9) && Ens_tags[iens] == "cB211b.072.64") { mult_T_IM=1;}
	    if(ixg==0 && (isg==(signed)sigmas_07.size()-9) && Ens_tags[iens] == "cD211a.054.96") { mult_T_IM=1;}


	    if(ixg==1 && (isg==(signed)sigmas_07.size()-4) && Ens_tags[iens] == "cB211b.072.64") { mult_T_IM=0.06;}

	    if(ixg==1 && (isg==(signed)sigmas_07.size()-5) && Ens_tags[iens] == "cB211b.072.64") { mult_T_IM=0.1;}
	    if(ixg==1 && (isg==(signed)sigmas_07.size()-5) && Ens_tags[iens] == "cD211a.054.96") { mult_T_IM=1;}

	    
	    if(ixg==1 && (isg==(signed)sigmas_07.size()-6) && Ens_tags[iens] == "cB211b.072.64") { mult_T_IM=0.5;}
	    if(ixg==1 && (isg==(signed)sigmas_07.size()-6) && Ens_tags[iens] == "cD211a.054.96") { mult_T_IM=1;}

	    if(ixg==1 && (isg==(signed)sigmas_07.size()-7) && Ens_tags[iens] == "cB211b.072.64") { mult_T_IM=1;}
	    if(ixg==1 && (isg==(signed)sigmas_07.size()-7) && Ens_tags[iens] == "cD211a.054.96") { mult_T_IM=1;}

	    if(ixg==1 && (isg==(signed)sigmas_07.size()-8) && Ens_tags[iens] == "cB211b.072.64") { mult_T_IM=1;}
	    if(ixg==1 && (isg==(signed)sigmas_07.size()-8) && Ens_tags[iens] == "cD211a.054.96") { mult_T_IM=1;}

	    if(ixg==1 && (isg==(signed)sigmas_07.size()-9) && Ens_tags[iens] == "cB211b.072.64") { mult_T_IM=1;}
	    if(ixg==1 && (isg==(signed)sigmas_07.size()-9) && Ens_tags[iens] == "cD211a.054.96") { mult_T_IM=1;}


	    if(ixg==2 && (isg==(signed)sigmas_07.size()-1) && Ens_tags[iens] == "cD211a.054.96") { mult_T_IM=1.26;}

	    
	    if(ixg==2 && (isg==(signed)sigmas_07.size()-2) && Ens_tags[iens] == "cB211b.072.64") { mult_T_IM=0.08;}
	     if(ixg==2 && (isg==(signed)sigmas_07.size()-2) && Ens_tags[iens] == "cD211a.054.96") { mult_T_IM=2;}

	    if(ixg==2 && (isg==(signed)sigmas_07.size()-3) && Ens_tags[iens] == "cD211a.054.96") { mult_T_IM=2;}
	    if(ixg==2 && (isg==(signed)sigmas_07.size()-4) && Ens_tags[iens] == "cD211a.054.96") { mult_T_IM=2;}

	    if(ixg==2 && (isg==(signed)sigmas_07.size()-5) && Ens_tags[iens] == "cB211b.072.64") { mult_T_IM=0.7;}
	    if(ixg==2 && (isg==(signed)sigmas_07.size()-5) && Ens_tags[iens] == "cD211a.054.96") { mult_T_IM=0.7;}

	    if(ixg==2 && (isg==(signed)sigmas_07.size()-6) && Ens_tags[iens] == "cB211b.072.64") { mult_T_IM=1;}
	    if(ixg==2 && (isg==(signed)sigmas_07.size()-6) && Ens_tags[iens] == "cD211a.054.96") { mult_T_IM=1;}

	    if(ixg==2 && (isg==(signed)sigmas_07.size()-7) && Ens_tags[iens] == "cB211b.072.64") { mult_T_IM=1;}
	    if(ixg==2 && (isg==(signed)sigmas_07.size()-7) && Ens_tags[iens] == "cD211a.054.96") { mult_T_IM=1;}

	    if(ixg==2 && (isg==(signed)sigmas_07.size()-8) && Ens_tags[iens] == "cB211b.072.64") { mult_T_IM=1;}
	    if(ixg==2 && (isg==(signed)sigmas_07.size()-8) && Ens_tags[iens] == "cD211a.054.96") { mult_T_IM=1;}

	    if(ixg==2 && (isg==(signed)sigmas_07.size()-9) && Ens_tags[iens] == "cB211b.072.64") { mult_T_IM=1;}
	    if(ixg==2 && (isg==(signed)sigmas_07.size()-9) && Ens_tags[iens] == "cD211a.054.96") { mult_T_IM=1;}


	    if(ixg==3 && (isg==(signed)sigmas_07.size()-1) && Ens_tags[iens] == "cB211b.072.64") { mult_T_IM=1;}
	    if(ixg==3 && (isg==(signed)sigmas_07.size()-1) && Ens_tags[iens] == "cD211a.054.96") { mult_T_IM=1;}

	    if(ixg==3 && (isg==(signed)sigmas_07.size()-2) && Ens_tags[iens] == "cB211b.072.64") { mult_T_IM=0.1;}
	    if(ixg==3 && (isg==(signed)sigmas_07.size()-2) && Ens_tags[iens] == "cD211a.054.96") { mult_T_IM=2;}

	  
	    if(ixg==3 && (isg==(signed)sigmas_07.size()-3) && Ens_tags[iens] == "cD211a.054.96") { mult_T_IM=2;}

	    if(ixg==3 && (isg==(signed)sigmas_07.size()-4) && Ens_tags[iens] == "cB211b.072.64") { mult_T_IM=1;}
	    if(ixg==3 && (isg==(signed)sigmas_07.size()-4) && Ens_tags[iens] == "cD211a.054.96") { mult_T_IM=1;}

	    if(ixg==3 && (isg==(signed)sigmas_07.size()-5) && Ens_tags[iens] == "cB211b.072.64") { mult_T_IM=1;}
	    if(ixg==3 && (isg==(signed)sigmas_07.size()-5) && Ens_tags[iens] == "cD211a.054.96") { mult_T_IM=1;}

	    if(ixg==3 && (isg==(signed)sigmas_07.size()-6) && Ens_tags[iens] == "cB211b.072.64") { mult_T_IM=2;}
	    if(ixg==3 && (isg==(signed)sigmas_07.size()-6) && Ens_tags[iens] == "cD211a.054.96") { mult_T_IM=2;}

	    if(ixg==3 && (isg==(signed)sigmas_07.size()-7) && Ens_tags[iens] == "cB211b.072.64") { mult_T_IM=4;}
	    if(ixg==3 && (isg==(signed)sigmas_07.size()-7) && Ens_tags[iens] == "cD211a.054.96") { mult_T_IM=4;}

	    if(ixg==3 && (isg==(signed)sigmas_07.size()-8) && Ens_tags[iens] == "cB211b.072.64") { mult_T_IM=3;}
	    if(ixg==3 && (isg==(signed)sigmas_07.size()-8) && Ens_tags[iens] == "cD211a.054.96") { mult_T_IM=3;}

	    if(ixg==3 && (isg==(signed)sigmas_07.size()-9) && Ens_tags[iens] == "cB211b.072.64") { mult_T_IM=3;}
	    if(ixg==3 && (isg==(signed)sigmas_07.size()-9) && Ens_tags[iens] == "cD211a.054.96") { mult_T_IM=3;}

	  
	    
	  
	  }
	  if(MESON=="B3s") {
	    if(ixg==0 && (isg==(signed)sigmas_07.size()-1) && Ens_tags[iens] == "cB211b.072.64") { mult_T_IM=0.03;}
	    
	    if(ixg==0 &&  (isg == (signed)sigmas_07.size()-2) && Ens_tags[iens] == "cB211b.072.64" ) { mult_T_IM=0.04;}
	    if(ixg==0 &&  (isg == (signed)sigmas_07.size()-2) && Ens_tags[iens] == "cD211a.054.96" ) { mult_T_IM=2.35;}
	    
	    if(ixg==0 &&  (isg == (signed)sigmas_07.size()-3) && Ens_tags[iens] == "cB211b.072.64" ) { mult_T_IM=0.05;}
	    if(ixg==0 &&  (isg == (signed)sigmas_07.size()-3) && Ens_tags[iens] == "cD211a.054.96" ) { mult_T_IM=2.5;}

	    if(ixg==0 &&  (isg == (signed)sigmas_07.size()-4) && Ens_tags[iens] == "cB211b.072.64" ) { mult_T_IM=0.04;}
	    if(ixg==0 && (isg==(signed)sigmas_07.size()-4) && Ens_tags[iens] == "cD211a.054.96") { mult_T_IM=3;}

	    if(ixg==0 &&  (isg == (signed)sigmas_07.size()-5) && Ens_tags[iens] == "cB211b.072.64" ) { mult_T_IM=0.04;}
	    if(ixg==0 && (isg==(signed)sigmas_07.size()-5) && Ens_tags[iens] == "cD211a.054.96") { mult_T_IM=3;}

	    if(ixg==0 &&  (isg == (signed)sigmas_07.size()-6) && Ens_tags[iens] == "cB211b.072.64" ) { mult_T_IM=0.04;}
	    if(ixg==0 && (isg==(signed)sigmas_07.size()-6) && Ens_tags[iens] == "cD211a.054.96") { mult_T_IM=3;}

	    
	    if(ixg==0 && (isg==(signed)sigmas_07.size()-7) && Ens_tags[iens] == "cD211a.054.96") { mult_T_IM=2;}

	    if(ixg==0 &&  (isg == (signed)sigmas_07.size()-8) && Ens_tags[iens] == "cB211b.072.64" ) { mult_T_IM=1;}
	    if(ixg==0 && (isg==(signed)sigmas_07.size()-8) && Ens_tags[iens] == "cD211a.054.96") { mult_T_IM=1;}

	    if(ixg==0 &&  (isg == (signed)sigmas_07.size()-9) && Ens_tags[iens] == "cB211b.072.64" ) { mult_T_IM=1;}
	    if(ixg==0 && (isg==(signed)sigmas_07.size()-9) && Ens_tags[iens] == "cD211a.054.96") { mult_T_IM=1;}


	    
	    if(ixg==1 && (isg==(signed)sigmas_07.size()-1) && Ens_tags[iens] == "cB211b.072.64") { mult_T_IM=0.07;}
	    if(ixg==1 &&  (isg == (signed)sigmas_07.size()-2) && Ens_tags[iens] == "cB211b.072.64" ) { mult_T_IM=0.04;}
	    
	    if(ixg==1 &&  (isg == (signed)sigmas_07.size()-3) && Ens_tags[iens] == "cB211b.072.64" ) { mult_T_IM=0.05;}
	    if(ixg==1 &&  (isg == (signed)sigmas_07.size()-3) && Ens_tags[iens] == "cD211a.054.96" ) { mult_T_IM=1;}
	    
	    if(ixg==1 &&  (isg == (signed)sigmas_07.size()-4) && Ens_tags[iens] == "cB211b.072.64" ) { mult_T_IM=0.02;}
	    if(ixg==1 &&  (isg == (signed)sigmas_07.size()-4) && Ens_tags[iens] == "cD211a.054.96" ) { mult_T_IM=1;}
	    
	    if(ixg==1 &&  (isg == (signed)sigmas_07.size()-5) && Ens_tags[iens] == "cB211b.072.64" ) { mult_T_IM=0.01;}
	    if(ixg==1 &&  (isg == (signed)sigmas_07.size()-5) && Ens_tags[iens] == "cD211a.054.96" ) { mult_T_IM=1;}

	    if(ixg==1 &&  (isg == (signed)sigmas_07.size()-6) && Ens_tags[iens] == "cD211a.054.96" ) { mult_T_IM=1;}

	    if(ixg==1 &&  (isg == (signed)sigmas_07.size()-7) && Ens_tags[iens] == "cB211b.072.64" ) { mult_T_IM=1;}
	    if(ixg==1 &&  (isg == (signed)sigmas_07.size()-7) && Ens_tags[iens] == "cD211a.054.96" ) { mult_T_IM=1;}
	    
	    if(ixg==1 &&  (isg == (signed)sigmas_07.size()-8) && Ens_tags[iens] == "cB211b.072.64" ) { mult_T_IM=1;}
	    if(ixg==1 &&  (isg == (signed)sigmas_07.size()-8) && Ens_tags[iens] == "cD211a.054.96" ) { mult_T_IM=1;}
	    
	    if(ixg==1 &&  (isg == (signed)sigmas_07.size()-9) && Ens_tags[iens] == "cB211b.072.64" ) { mult_T_IM=1;}
	    if(ixg==1 &&  (isg == (signed)sigmas_07.size()-9) && Ens_tags[iens] == "cD211a.054.96" ) { mult_T_IM=1;}

	    
	    if(ixg==2 &&  (isg == (signed)sigmas_07.size()-1) && Ens_tags[iens] == "cB211b.072.64" ) { mult_T_IM=0.03;}
	    if(ixg==2 &&  (isg == (signed)sigmas_07.size()-1) && Ens_tags[iens] == "cD211a.054.96" ) { mult_T_IM=0.2;}
	    
	    if(ixg==2 &&  (isg == (signed)sigmas_07.size()-2) && Ens_tags[iens] == "cB211b.072.64") { mult_T_IM=0.02;}
	    if(ixg==2 &&  (isg == (signed)sigmas_07.size()-2) && Ens_tags[iens] == "cD211a.054.96" ) { mult_T_IM=0.2;}
	    
	    if(ixg==2 &&  (isg == (signed)sigmas_07.size()-3) && Ens_tags[iens] == "cB211b.072.64") { mult_T_IM=0.003;}
	    if(ixg==2 &&  (isg == (signed)sigmas_07.size()-3) && Ens_tags[iens] == "cD211a.054.96" ) { mult_T_IM=0.09;}
	    
	    if(ixg==2 &&  (isg == (signed)sigmas_07.size()-4) && Ens_tags[iens] == "cB211b.072.64") { mult_T_IM=0.01;}
	    if(ixg==2 &&  (isg == (signed)sigmas_07.size()-4) && Ens_tags[iens] == "cD211a.054.96" ) { mult_T_IM=0.2;}

	    
	    if(ixg==2 &&  (isg == (signed)sigmas_07.size()-5) && Ens_tags[iens] == "cB211b.072.64") { mult_T_IM=0.04;}
	    if(ixg==2 &&  (isg == (signed)sigmas_07.size()-5) && Ens_tags[iens] == "cD211a.054.96" ) { mult_T_IM=3;}

	    if(ixg==2 &&  (isg == (signed)sigmas_07.size()-6) && Ens_tags[iens] == "cD211a.054.96" ) { mult_T_IM=1;}
	    
	    if(ixg==2 &&  (isg == (signed)sigmas_07.size()-7) && Ens_tags[iens] == "cB211b.072.64") { mult_T_IM=1;}
	    if(ixg==2 &&  (isg == (signed)sigmas_07.size()-7) && Ens_tags[iens] == "cD211a.054.96" ) { mult_T_IM=1;}
	    
	    if(ixg==2 &&  (isg == (signed)sigmas_07.size()-8) && Ens_tags[iens] == "cB211b.072.64") { mult_T_IM=1;}
	    if(ixg==2 &&  (isg == (signed)sigmas_07.size()-8) && Ens_tags[iens] == "cD211a.054.96" ) { mult_T_IM=1;}
	    
	    if(ixg==2 &&  (isg == (signed)sigmas_07.size()-9) && Ens_tags[iens] == "cB211b.072.64") { mult_T_IM=1;}
	    if(ixg==2 &&  (isg == (signed)sigmas_07.size()-9) && Ens_tags[iens] == "cD211a.054.96" ) { mult_T_IM=1;}


	    if(ixg==3 &&  (isg == (signed)sigmas_07.size()-1) && Ens_tags[iens] == "cB211b.072.64" ) { mult_T_IM=0.01;}
	    if(ixg==3 &&  (isg == (signed)sigmas_07.size()-1) && Ens_tags[iens] == "cD211a.054.96" ) { mult_T_IM=0.1;}
	    
	    if(ixg==3 &&  (isg == (signed)sigmas_07.size()-2) && Ens_tags[iens] == "cB211b.072.64") { mult_T_IM=0.001;}
	    if(ixg==3 &&  (isg == (signed)sigmas_07.size()-2) && Ens_tags[iens] == "cD211a.054.96" ) { mult_T_IM=0.15;}
	    
	    if(ixg==3 &&  (isg == (signed)sigmas_07.size()-3) && Ens_tags[iens] == "cB211b.072.64") { mult_T_IM=0.003;}

	    
	    if(ixg==3 &&  (isg == (signed)sigmas_07.size()-4) && Ens_tags[iens] == "cB211b.072.64") { mult_T_IM=1;}
	    if(ixg==3 &&  (isg == (signed)sigmas_07.size()-4) && Ens_tags[iens] == "cD211a.054.96" ) { mult_T_IM=1;}

	    if(ixg==3 &&  (isg == (signed)sigmas_07.size()-5) && Ens_tags[iens] == "cB211b.072.64") { mult_T_IM=1;}
	    if(ixg==3 &&  (isg == (signed)sigmas_07.size()-5) && Ens_tags[iens] == "cD211a.054.96" ) { mult_T_IM=1;}

	    if(ixg==3 &&  (isg == (signed)sigmas_07.size()-6) && Ens_tags[iens] == "cB211b.072.64") { mult_T_IM=1;}
	    if(ixg==3 &&  (isg == (signed)sigmas_07.size()-6) && Ens_tags[iens] == "cD211a.054.96" ) { mult_T_IM=1;}
	    
	    if(ixg==3 &&  (isg == (signed)sigmas_07.size()-7) && Ens_tags[iens] == "cB211b.072.64") { mult_T_IM=1;}
	    if(ixg==3 &&  (isg == (signed)sigmas_07.size()-7) && Ens_tags[iens] == "cD211a.054.96" ) { mult_T_IM=1;}
	    
	    if(ixg==3 &&  (isg == (signed)sigmas_07.size()-8) && Ens_tags[iens] == "cB211b.072.64") { mult_T_IM=1;}
	    if(ixg==3 &&  (isg == (signed)sigmas_07.size()-8) && Ens_tags[iens] == "cD211a.054.96" ) { mult_T_IM=1;}

	    if(ixg==3 &&  (isg == (signed)sigmas_07.size()-9) && Ens_tags[iens] == "cB211b.072.64") { mult_T_IM=1;}
	    if(ixg==3 &&  (isg == (signed)sigmas_07.size()-9) && Ens_tags[iens] == "cD211a.054.96" ) { mult_T_IM=1;}

	    	   


	  }


	  distr_t OFFSET_IM_GAUSS= ((Use_preconditioning==false)?0.0:1.0)*(FACT*phi_MT*K_IM_distr( eff_M_Td*a_distr, Eg_MB_off, s)+ 0.0*FACT*phi_prime_MT*K_IM_distr( eff_M_prime_Td*a_distr, Eg_MB_off, s) );  
	  
	  distr_t IM_MB_Gauss_sm = Get_Laplace_transfo_tmin( Eg_MB_off.ave(),  s, th*a_distr.ave(),  Nts[iens], 1, tmax_new_im, prec_07, SM_FF+"_IM_ixg_"+to_string(ixg+1),K_IM, Corr_T, syst_T, mult_T_IM ,  l_re_T, MODE_FF, "E0_"+to_string_with_precision(E0_fact,1), TAG_CURR+"MB_T_"+Ens_tags[iens], Ag_target,0, FACT, OFFSET_IM_GAUSS , preco_tag+MESON+"_07_FF_Tw_"+to_string(t_07_s_HLT), cov_T, fake_func,0, fake_func_d ,  1 , 4.0, 0.0,1);

	  IM_MB_Gauss_sm = IM_MB_Gauss_sm.ave() +  (IM_MB_Gauss_sm - IM_MB_Gauss_sm.ave())*(sqrt( pow(syst_T,2) + pow(IM_MB_Gauss_sm.err(),2)))/IM_MB_Gauss_sm.err();
	  

	  F_T_d_MB_IM_Gauss_sm_list[ixg][isg].distr_list.push_back( IM_MB_Gauss_sm +  ((Use_preconditioning==false)?0.0:1.0)*(FACT*phi_MT*K_IM_distr( eff_M_Td*a_distr, Eg_MB_off, s)+ 0.0*FACT*phi_prime_MT*K_IM_distr( eff_M_prime_Td*a_distr, Eg_MB_off, s) ));
          */

	  
	  F_T_d_MB_IM_Gauss_sm_list[ixg][isg].distr_list.push_back( 0.0*a_distr);


	  SM_FF = "FF_Exp";
	  

	  if(virtuality_scan) {
	    
	    for(int vir=0;vir<(signed)virtualities.size(); vir++) {

	      F_T_d_RE_virt_scan[isg].distr_list[vir] = Get_Laplace_transfo ( M_P.ave()*virtualities[vir],  s, th*a_distr.ave(),  Nts[iens], tmax_new, prec_07, SM_FF+"_RE",K_RE, Corr_T, syst_T, mult_T_RE ,  l_re_T, MODE_FF, "E0_"+to_string_with_precision(E0_fact,1), TAG_CURR+"VIRT_SCAN_T_"+Ens_tags[iens], Ag_target,0, FACT, 0.0 , MESON+"_07_FF_Tw_"+to_string(t_07_s_HLT), cov_T, fake_func,0, fake_func_d ,  1 , 4.0, 0.0,1) +  ((Use_preconditioning==false)?0.0:1.0)*(FACT*phi_MT*K_RE_distr_sub( eff_M_Td*a_distr, M_P*virtualities[vir], s) + 0.0*FACT*phi_prime_MT*K_RE_distr_sub( eff_M_prime_Td*a_distr, M_P*virtualities[vir], s) - REMINDER_MB);
	      
	      F_T_d_IM_virt_scan[isg].distr_list[vir] = Get_Laplace_transfo ( M_P.ave()*virtualities[vir],  s, th*a_distr.ave(),  Nts[iens], tmax_new_im, prec_07, SM_FF+"_IM",K_IM, Corr_T, syst_T, mult_T_IM ,  l_re_T, MODE_FF, "E0_"+to_string_with_precision(E0_fact,1), TAG_CURR+"VIRT_SCAN_T_"+Ens_tags[iens], Ag_target,0, FACT, 0.0 , MESON+"_07_FF_Tw_"+to_string(t_07_s_HLT), cov_T, fake_func,0, fake_func_d ,  1 , 4.0, 0.0,1)  +  ((Use_preconditioning==false)?0.0:1.0)*(FACT*phi_MT*K_IM_distr( eff_M_Td*a_distr, M_P*virtualities[vir], s)+0.0*FACT*phi_prime_MT*K_IM_distr( eff_M_prime_Td*a_distr, M_P*virtualities[vir], s));

	    }

	  }
	 	  
	  
	}
	//if virtuality scan, print results
	if(virtuality_scan) {
	  for(int isg=0;isg<(signed)sigmas_07.size();isg++) {
	    Print_To_File({}, {(virtualities*Get_id_jack_distr_list(virtualities.size(),Njacks)*M_P/a_distr).ave(), F_T_d_RE_virt_scan[isg].ave(), F_T_d_RE_virt_scan[isg].err(), F_T_d_IM_virt_scan[isg].ave(), F_T_d_IM_virt_scan[isg].err() }, path_out+"/FF_d_II/"+Ens_tags[iens]+"/VIRT_SCAN_sm_"+to_string_with_precision(sigmas_07[isg],3)+"_xg_"+to_string_with_precision(0.1+0.1*ixg,2), "", "");
	  }
	}
      }  
      
    }

  }

  

  //print results
  //per kinematic
  for(int ixg=0;ixg<n_xg;ixg++) {

    Print_To_File( Ens_tags, { F_T_u_list[ixg].ave(), F_T_u_list[ixg].err()     } , path_out+"/FF_u/"+TAG_CURR+"F_T_xg_"+to_string_with_precision(xg_list.ave(ixg),2), "", "");
    Print_To_File( Ens_tags, { F_T_d_I_list[ixg].ave(), F_T_d_I_list[ixg].err() , F_T_d_sp_I_list[ixg].ave(), F_T_d_sp_I_list[ixg].err(), F_T_d_I_ori_list[ixg].ave(), F_T_d_I_ori_list[ixg].err()    } , path_out+"/FF_d_I/"+TAG_CURR+"F_T_xg_"+to_string_with_precision(xg_list.ave(ixg),2), "", "");
    Print_To_File( Ens_tags, { FV_T_u_real_list[ixg].ave(), FV_T_u_real_list[ixg].err()     } , path_out+"/FF_u/"+TAG_CURR+"FV_T_real_xg_"+to_string_with_precision(xg_list.ave(ixg),2), "", "");
    Print_To_File( Ens_tags, { FV_T_d_real_list[ixg].ave(), FV_T_d_real_list[ixg].err()     } , path_out+"/FF_d/"+TAG_CURR+"FV_T_real_xg_"+to_string_with_precision(xg_list.ave(ixg),2), "", "");
    Print_To_File( Ens_tags, { FA_T_u_real_list[ixg].ave(), FA_T_u_real_list[ixg].err()     } , path_out+"/FF_u/"+TAG_CURR+"FA_T_real_xg_"+to_string_with_precision(xg_list.ave(ixg),2), "", "");
    Print_To_File( Ens_tags, { FA_T_d_real_list[ixg].ave(), FA_T_d_real_list[ixg].err()     } , path_out+"/FF_d/"+TAG_CURR+"FA_T_real_xg_"+to_string_with_precision(xg_list.ave(ixg),2), "", "");

    Print_To_File( Ens_tags, { FV_T_d_sp_real_list[ixg].ave(), FV_T_d_sp_real_list[ixg].err()     } , path_out+"/FF_d/"+TAG_CURR+"FV_T_sp_real_xg_"+to_string_with_precision(xg_list.ave(ixg),2), "", "");
    Print_To_File( Ens_tags, { FA_T_d_sp_real_list[ixg].ave(), FA_T_d_sp_real_list[ixg].err()     } , path_out+"/FF_d/"+TAG_CURR+"FA_T_sp_real_xg_"+to_string_with_precision(xg_list.ave(ixg),2), "", "");

    Print_To_File( Ens_tags, { F_T_u_VMD_list[ixg].ave(), F_T_u_VMD_list[ixg].err(), F_T_u_VMD_spectre_list[ixg].ave(), F_T_u_VMD_spectre_list[ixg].err()     } , path_out+"/FF_u/"+TAG_CURR+"VMD_F_T_xg_"+to_string_with_precision(xg_list.ave(ixg),2), "", "");

    
    if(!Skip_spectral_reconstruction_07) {
    for(int isg=0;isg<(signed)sigmas_07.size();isg++) {

   
      Print_To_File( Ens_tags, { F_T_d_MB_RE_sm_list[ixg][isg].ave(), F_T_d_MB_RE_sm_list[ixg][isg].err(), F_T_d_MB_IM_sm_list[ixg][isg].ave(), F_T_d_MB_IM_sm_list[ixg][isg].err()     } , path_out+"/FF_d_II/"+TAG_CURR+preco_tag+"F_MB_T_sm_"+to_string_with_precision(sigmas_07[isg],3)+"_xg_"+to_string_with_precision(xg_list.ave(ixg),2), "", "");

      Print_To_File( Ens_tags, {  F_T_d_MB_IM_Gauss_sm_list[ixg][isg].ave(), F_T_d_MB_IM_Gauss_sm_list[ixg][isg].err()     } , path_out+"/FF_d_II/"+TAG_CURR+preco_tag+"F_MB_T_Gauss_sm_"+to_string_with_precision(sigmas_07[isg],3)+"_xg_"+to_string_with_precision(xg_list.ave(ixg),2), "", "");


    }
    }
  }
  //per ensemble
  for(int iens=0; iens<Nens;iens++) {
    distr_t_list F_T_u_per_ens(UseJack),  F_T_d_I_per_ens(UseJack);
    distr_t_list FV_T_u_real_per_ens(UseJack),  FV_T_d_real_per_ens(UseJack);
    distr_t_list FA_T_u_real_per_ens(UseJack),  FA_T_d_real_per_ens(UseJack);
    distr_t_list FA_T_d_I_real_per_ens(UseJack), FA_T_d_II_real_per_ens(UseJack);
    distr_t_list FV_T_d_I_real_per_ens(UseJack), FV_T_d_II_real_per_ens(UseJack);
    distr_t_list FA_T_d_sp_real_per_ens(UseJack), FV_T_d_sp_real_per_ens(UseJack);
    distr_t_list F_T_u_VMD_per_ens(UseJack), F_T_u_VMD_spectre_per_ens(UseJack);
    distr_t_list F_T_d_sp_I_per_ens(UseJack);
    distr_t_list F_T_d_I_ori_per_ens(UseJack);
    distr_t_list F_T_d_I_ori_sp_per_ens(UseJack);

    distr_t_list reminder_F_T_d_per_ens(UseJack), reminder_F_T_d_sp_per_ens(UseJack);


    distr_t_list F_T_d_I_per_ens_15(UseJack), F_T_d_I_per_ens_7(UseJack);

    for(int ixg=0;ixg<n_xg;ixg++) {
      F_T_u_per_ens.distr_list.push_back( F_T_u_list[ixg].distr_list[iens]);
      F_T_d_I_per_ens.distr_list.push_back( F_T_d_I_list[ixg].distr_list[iens]);
      F_T_d_I_ori_per_ens.distr_list.push_back( F_T_d_I_ori_list[ixg].distr_list[iens]);
      F_T_d_I_ori_sp_per_ens.distr_list.push_back( F_T_d_I_ori_sp_list[ixg].distr_list[iens]);
      F_T_d_I_per_ens_15.distr_list.push_back( F_T_d_I_list_15[ixg].distr_list[iens]);
      F_T_d_I_per_ens_7.distr_list.push_back( F_T_d_I_list_7[ixg].distr_list[iens]);
      
      F_T_d_sp_I_per_ens.distr_list.push_back( F_T_d_sp_I_list[ixg].distr_list[iens]);
      FV_T_u_real_per_ens.distr_list.push_back( FV_T_u_real_list[ixg].distr_list[iens]);
      FV_T_d_real_per_ens.distr_list.push_back( FV_T_d_real_list[ixg].distr_list[iens]);
      FA_T_u_real_per_ens.distr_list.push_back( FA_T_u_real_list[ixg].distr_list[iens]);
      FA_T_d_real_per_ens.distr_list.push_back( FA_T_d_real_list[ixg].distr_list[iens]);

      FA_T_d_I_real_per_ens.distr_list.push_back( FA_T_d_I_real_list[ixg].distr_list[iens]);
      FA_T_d_II_real_per_ens.distr_list.push_back( FA_T_d_II_real_list[ixg].distr_list[iens]);

      FV_T_d_I_real_per_ens.distr_list.push_back( FV_T_d_I_real_list[ixg].distr_list[iens]);
      FV_T_d_II_real_per_ens.distr_list.push_back( FV_T_d_II_real_list[ixg].distr_list[iens]);

      FV_T_d_sp_real_per_ens.distr_list.push_back( FV_T_d_sp_real_list[ixg].distr_list[iens]);
      FA_T_d_sp_real_per_ens.distr_list.push_back( FA_T_d_sp_real_list[ixg].distr_list[iens]);

      F_T_u_VMD_per_ens.distr_list.push_back( F_T_u_VMD_list[ixg].distr_list[iens]);
      F_T_u_VMD_spectre_per_ens.distr_list.push_back( F_T_u_VMD_spectre_list[ixg].distr_list[iens]);

      reminder_F_T_d_per_ens.distr_list.push_back( reminder_F_T_d_list[ixg].distr_list[iens]);
      reminder_F_T_d_sp_per_ens.distr_list.push_back( reminder_F_T_d_sp_list[ixg].distr_list[iens]);

      
      distr_t_list F_T_d_MB_RE_sm_per_ens_per_kin, F_T_d_MB_IM_sm_per_ens_per_kin;

      distr_t_list F_T_d_MB_RE_sm_per_ens_per_kin_15, F_T_d_MB_IM_sm_per_ens_per_kin_15;
      distr_t_list F_T_d_MB_RE_sm_per_ens_per_kin_7, F_T_d_MB_IM_sm_per_ens_per_kin_7;
      
      
      distr_t_list F_T_d_MB_RE_VMD_sm_per_ens_per_kin, F_T_d_MB_IM_VMD_sm_per_ens_per_kin;
      
      
      distr_t_list F_T_d_MB_RE_VMD_II_state_sm_per_ens_per_kin, F_T_d_MB_IM_VMD_II_state_sm_per_ens_per_kin;
      
      
      distr_t_list F_T_d_MB_RE_VMD_III_state_sm_per_ens_per_kin, F_T_d_MB_IM_VMD_III_state_sm_per_ens_per_kin;

      distr_t_list  F_T_d_MB_IM_Gauss_sm_per_ens_per_kin;
       
      
      if(!Skip_spectral_reconstruction_07) {
	for(int iss=0; iss<(signed)sigmas_07_w0.size(); iss++) {
	
	  F_T_d_MB_RE_VMD_sm_per_ens_per_kin.distr_list.push_back( F_T_d_MB_RE_VMD_sm_list[ixg][iss].distr_list[iens]);
	  F_T_d_MB_IM_VMD_sm_per_ens_per_kin.distr_list.push_back( F_T_d_MB_IM_VMD_sm_list[ixg][iss].distr_list[iens]);

	  F_T_d_MB_RE_VMD_II_state_sm_per_ens_per_kin.distr_list.push_back( F_T_d_MB_RE_VMD_II_state_sm_list[ixg][iss].distr_list[iens]);
	  F_T_d_MB_IM_VMD_II_state_sm_per_ens_per_kin.distr_list.push_back( F_T_d_MB_IM_VMD_II_state_sm_list[ixg][iss].distr_list[iens]);

	 
	  F_T_d_MB_RE_VMD_III_state_sm_per_ens_per_kin.distr_list.push_back( F_T_d_MB_RE_VMD_III_state_sm_list[ixg][iss].distr_list[iens]);
	  F_T_d_MB_IM_VMD_III_state_sm_per_ens_per_kin.distr_list.push_back( F_T_d_MB_IM_VMD_III_state_sm_list[ixg][iss].distr_list[iens]);
	}
	for(int isg=0;isg<(signed)sigmas_07.size(); isg++) {
	

	  F_T_d_MB_RE_sm_per_ens_per_kin.distr_list.push_back( F_T_d_MB_RE_sm_list[ixg][isg].distr_list[iens]);
	  F_T_d_MB_IM_sm_per_ens_per_kin.distr_list.push_back( F_T_d_MB_IM_sm_list[ixg][isg].distr_list[iens]);
	  F_T_d_MB_IM_Gauss_sm_per_ens_per_kin.distr_list.push_back( F_T_d_MB_IM_Gauss_sm_list[ixg][isg].distr_list[iens]);
	  
	  F_T_d_MB_RE_sm_per_ens_per_kin_15.distr_list.push_back( F_T_d_MB_RE_sm_list_15[ixg][isg].distr_list[iens]);
	  F_T_d_MB_IM_sm_per_ens_per_kin_15.distr_list.push_back( F_T_d_MB_IM_sm_list_15[ixg][isg].distr_list[iens]);

	  F_T_d_MB_RE_sm_per_ens_per_kin_7.distr_list.push_back( F_T_d_MB_RE_sm_list_7[ixg][isg].distr_list[iens]);
	  F_T_d_MB_IM_sm_per_ens_per_kin_7.distr_list.push_back( F_T_d_MB_IM_sm_list_7[ixg][isg].distr_list[iens]);

	  
	  
	}

	
      
	Print_To_File( {}, {sigmas_07, sigma_simulated[ixg], F_T_d_MB_RE_sm_per_ens_per_kin.ave(), F_T_d_MB_RE_sm_per_ens_per_kin.err(), F_T_d_MB_IM_sm_per_ens_per_kin.ave(), F_T_d_MB_IM_sm_per_ens_per_kin.err()}, path_out+"/FF_d_II/"+Ens_tags[iens]+"/"+TAG_CURR+preco_tag+"F_MB_T_xg_"+to_string_with_precision(xg_list.ave(ixg),2), "", "");
	Print_To_File( {}, {sigmas_07, sigma_simulated[ixg], F_T_d_MB_IM_Gauss_sm_per_ens_per_kin.ave(), F_T_d_MB_IM_Gauss_sm_per_ens_per_kin.err()}, path_out+"/FF_d_II/"+Ens_tags[iens]+"/"+TAG_CURR+preco_tag+"F_MB_Gauss_T_xg_"+to_string_with_precision(xg_list.ave(ixg),2), "", "");
	Print_To_File( {}, {sigmas_07_w0, F_T_d_MB_RE_VMD_sm_per_ens_per_kin.ave(), F_T_d_MB_RE_VMD_sm_per_ens_per_kin.err(), F_T_d_MB_IM_VMD_sm_per_ens_per_kin.ave(), F_T_d_MB_IM_VMD_sm_per_ens_per_kin.err()}, path_out+"/FF_d_II/"+Ens_tags[iens]+"/"+TAG_CURR+"VMD_F_MB_T_xg_"+to_string_with_precision(xg_list.ave(ixg),2), "", "");
	Print_To_File( {}, {sigmas_07_w0, F_T_d_MB_RE_VMD_II_state_sm_per_ens_per_kin.ave(), F_T_d_MB_RE_VMD_II_state_sm_per_ens_per_kin.err(), F_T_d_MB_IM_VMD_II_state_sm_per_ens_per_kin.ave(), F_T_d_MB_IM_VMD_II_state_sm_per_ens_per_kin.err()}, path_out+"/FF_d_II/"+Ens_tags[iens]+"/"+TAG_CURR+"VMD_MB_II_state_F_T_xg_"+to_string_with_precision(xg_list.ave(ixg),2), "", "");
	Print_To_File( {}, {sigmas_07_w0, F_T_d_MB_RE_VMD_III_state_sm_per_ens_per_kin.ave(), F_T_d_MB_RE_VMD_III_state_sm_per_ens_per_kin.err(), F_T_d_MB_IM_VMD_III_state_sm_per_ens_per_kin.ave(), F_T_d_MB_IM_VMD_III_state_sm_per_ens_per_kin.err()}, path_out+"/FF_d_II/"+Ens_tags[iens]+"/"+TAG_CURR+"VMD_MB_III_state_F_T_xg_"+to_string_with_precision(xg_list.ave(ixg),2), "", "");
	Print_To_File( {}, {sigmas_07, sigma_simulated[ixg], (F_T_d_MB_RE_sm_per_ens_per_kin+F_T_u_per_ens[ixg]+ F_T_d_I_per_ens[ixg]).ave(), (F_T_d_MB_RE_sm_per_ens_per_kin+F_T_u_per_ens[ixg]+F_T_d_I_per_ens[ixg]).err(), F_T_d_MB_IM_sm_per_ens_per_kin.ave(), F_T_d_MB_IM_sm_per_ens_per_kin.err()}, path_out+"/FF/"+Ens_tags[iens]+"/"+TAG_CURR+preco_tag+"F_MB_T_xg_"+to_string_with_precision(xg_list.ave(ixg),2), "", "");
	Print_To_File( {}, {sigmas_07, sigma_simulated[ixg], (F_T_d_MB_RE_sm_per_ens_per_kin+F_T_d_I_per_ens[ixg]).ave(), (F_T_d_MB_RE_sm_per_ens_per_kin+F_T_d_I_per_ens[ixg]).err(), F_T_d_MB_IM_sm_per_ens_per_kin.ave(), F_T_d_MB_IM_sm_per_ens_per_kin.err()}, path_out+"/FF_d/"+Ens_tags[iens]+"/"+TAG_CURR+preco_tag+"F_MB_T_xg_"+to_string_with_precision(xg_list.ave(ixg),2), "", "");

	Print_To_File( {}, {sigmas_07, sigma_simulated[ixg], (F_T_d_MB_RE_sm_per_ens_per_kin_15+F_T_d_I_per_ens_15[ixg]).ave(), (F_T_d_MB_RE_sm_per_ens_per_kin_15+F_T_d_I_per_ens_15[ixg]).err(), F_T_d_MB_IM_sm_per_ens_per_kin_15.ave(), F_T_d_MB_IM_sm_per_ens_per_kin_15.err()}, path_out+"/FF_d/"+Ens_tags[iens]+"/"+TAG_CURR+preco_tag+"F_15_MB_T_xg_"+to_string_with_precision(xg_list.ave(ixg),2), "", "");

	Print_To_File( {}, {sigmas_07, sigma_simulated[ixg], (F_T_d_MB_RE_sm_per_ens_per_kin_7+F_T_d_I_per_ens_7[ixg]).ave(), (F_T_d_MB_RE_sm_per_ens_per_kin_7+F_T_d_I_per_ens_7[ixg]).err(), F_T_d_MB_IM_sm_per_ens_per_kin_7.ave(), F_T_d_MB_IM_sm_per_ens_per_kin_7.err()}, path_out+"/FF_d/"+Ens_tags[iens]+"/"+TAG_CURR+preco_tag+"F_7_MB_T_xg_"+to_string_with_precision(xg_list.ave(ixg),2), "", "");
		
		
      }
    }

    Print_To_File( {}, {xg_list.ave(), F_T_u_per_ens.ave(), F_T_u_per_ens.err()     } , path_out+"/FF_u/"+Ens_tags[iens]+"/"+TAG_CURR+"F_T", "", "");
    Print_To_File( {}, {xg_list.ave(), F_T_d_I_per_ens.ave(), F_T_d_I_per_ens.err(), F_T_d_sp_I_per_ens.ave(), F_T_d_sp_I_per_ens.err(), F_T_d_I_ori_per_ens.ave(), F_T_d_I_ori_per_ens.err() , F_T_d_I_ori_sp_per_ens.ave(), F_T_d_I_ori_sp_per_ens.err()    } , path_out+"/FF_d_I/"+Ens_tags[iens]+"/"+TAG_CURR+"F_T", "", "");
    Print_To_File( {}, {xg_list.ave(),  FV_T_u_real_per_ens.ave(), FV_T_u_real_per_ens.err()     } , path_out+"/FF_u/"+Ens_tags[iens]+"/"+TAG_CURR+"FV_T_real", "", "");
    Print_To_File( {}, {xg_list.ave(), FV_T_d_real_per_ens.ave(), FV_T_d_real_per_ens.err()     } , path_out+"/FF_d/"+Ens_tags[iens]+"/"+TAG_CURR+"FV_T_real", "", "");

    Print_To_File( {}, {xg_list.ave(),  FA_T_u_real_per_ens.ave(), FA_T_u_real_per_ens.err()     } , path_out+"/FF_u/"+Ens_tags[iens]+"/"+TAG_CURR+"FA_T_real", "", "");
    Print_To_File( {}, {xg_list.ave(), FA_T_d_real_per_ens.ave(), FA_T_d_real_per_ens.err()     } , path_out+"/FF_d/"+Ens_tags[iens]+"/"+TAG_CURR+"FA_T_real", "", "");

    Print_To_File( {}, {xg_list.ave(), FA_T_d_I_real_per_ens.ave(), FA_T_d_I_real_per_ens.err()     } , path_out+"/FF_d_I/"+Ens_tags[iens]+"/"+TAG_CURR+"FA_T_real", "", "");
    Print_To_File( {}, {xg_list.ave(), FA_T_d_II_real_per_ens.ave(), FA_T_d_II_real_per_ens.err()     } , path_out+"/FF_d_II/"+Ens_tags[iens]+"/"+TAG_CURR+"FA_T_real", "", "");

    Print_To_File( {}, {xg_list.ave(), FV_T_d_I_real_per_ens.ave(), FV_T_d_I_real_per_ens.err()     } , path_out+"/FF_d_I/"+Ens_tags[iens]+"/"+TAG_CURR+"FV_T_real", "", "");
    Print_To_File( {}, {xg_list.ave(), FV_T_d_II_real_per_ens.ave(), FV_T_d_II_real_per_ens.err()     } , path_out+"/FF_d_II/"+Ens_tags[iens]+"/"+TAG_CURR+"FV_T_real", "", "");

    Print_To_File( {}, {xg_list.ave(),  FA_T_d_sp_real_per_ens.ave(), FA_T_d_sp_real_per_ens.err()     } , path_out+"/FF_d/"+Ens_tags[iens]+"/"+TAG_CURR+"FA_T_sp_real", "", "");
    Print_To_File( {}, {xg_list.ave(), FV_T_d_sp_real_per_ens.ave(), FV_T_d_sp_real_per_ens.err()     } , path_out+"/FF_d/"+Ens_tags[iens]+"/"+TAG_CURR+"FV_T_sp_real", "", "");
    
    Print_To_File( {}, {xg_list.ave(),  F_T_u_VMD_per_ens.ave(), F_T_u_VMD_per_ens.err(), F_T_u_VMD_spectre_per_ens.ave(), F_T_u_VMD_spectre_per_ens.err()     } , path_out+"/FF_u/"+Ens_tags[iens]+"/"+TAG_CURR+"VMD_F_T", "", "");

      Print_To_File( {}, {xg_list.ave(), reminder_F_T_d_per_ens.ave(), reminder_F_T_d_per_ens.err(), reminder_F_T_d_sp_per_ens.ave(), reminder_F_T_d_sp_per_ens.err() } , path_out+"/FF_d_I/"+Ens_tags[iens]+"/"+TAG_CURR+"reminder_F_T", "", "");

    
    
    if(!Skip_spectral_reconstruction_07) {
    for(int isg=0;isg<(signed)sigmas_07.size(); isg++) {

    
      distr_t_list F_T_d_MB_RE_sm_per_ens_per_sigma, F_T_d_MB_IM_sm_per_ens_per_sigma;
      distr_t_list F_T_d_MB_IM_Gauss_sm_per_ens_per_sigma;
      
       for(int ixg=0;ixg<n_xg;ixg++) {

	 F_T_d_MB_RE_sm_per_ens_per_sigma.distr_list.push_back( F_T_d_MB_RE_sm_list[ixg][isg].distr_list[iens]);
	 F_T_d_MB_IM_sm_per_ens_per_sigma.distr_list.push_back( F_T_d_MB_IM_sm_list[ixg][isg].distr_list[iens]);
	 F_T_d_MB_IM_Gauss_sm_per_ens_per_sigma.distr_list.push_back( F_T_d_MB_IM_Gauss_sm_list[ixg][isg].distr_list[iens]);
       }


       Print_To_File( {}, {xg_list.ave(), F_T_d_MB_RE_sm_per_ens_per_sigma.ave(), F_T_d_MB_RE_sm_per_ens_per_sigma.err(), F_T_d_MB_IM_sm_per_ens_per_sigma.ave(), F_T_d_MB_IM_sm_per_ens_per_sigma.err()}, path_out+"/FF_d_II/"+Ens_tags[iens]+"/"+TAG_CURR+preco_tag+"F_MB_T_sm_"+to_string_with_precision(sigmas_07[isg],3), "", "");
       Print_To_File( {}, {xg_list.ave(),  F_T_d_MB_IM_Gauss_sm_per_ens_per_sigma.ave(), F_T_d_MB_IM_Gauss_sm_per_ens_per_sigma.err()}, path_out+"/FF_d_II/"+Ens_tags[iens]+"/"+TAG_CURR+preco_tag+"F_MB_Gauss_T_sm_"+to_string_with_precision(sigmas_07[isg],3), "", "");
       Print_To_File( {}, {xg_list.ave(), (F_T_d_MB_RE_sm_per_ens_per_sigma + F_T_u_per_ens+ F_T_d_I_per_ens).ave(), (F_T_d_MB_RE_sm_per_ens_per_sigma + F_T_u_per_ens+ F_T_d_I_per_ens).err(),  F_T_d_MB_IM_sm_per_ens_per_sigma.ave(), F_T_d_MB_IM_sm_per_ens_per_sigma.err()}, path_out+"/FF/"+Ens_tags[iens]+"/"+TAG_CURR+preco_tag+"F_MB_T_sm_"+to_string_with_precision(sigmas_07[isg],3), "", "");
       Print_To_File( {}, {xg_list.ave(), (F_T_d_MB_RE_sm_per_ens_per_sigma + F_T_d_I_per_ens).ave(), (F_T_d_MB_RE_sm_per_ens_per_sigma + F_T_d_I_per_ens).err(),  F_T_d_MB_IM_sm_per_ens_per_sigma.ave(), F_T_d_MB_IM_sm_per_ens_per_sigma.err()}, path_out+"/FF_d/"+Ens_tags[iens]+"/"+TAG_CURR+preco_tag+"F_MB_T_sm_"+to_string_with_precision(sigmas_07[isg],3), "", "");
       
    }
    //save jackknives
    boost::filesystem::create_directory(path_out+"/FF_d/"+Ens_tags[iens]+"/jackknife");
    for(int ixg=0;ixg<n_xg;ixg++) {

      boost::filesystem::create_directory(path_out+"/FF_d/"+Ens_tags[iens]+"/jackknife/ixg_"+to_string(ixg));
      for(int isg=0;isg<(signed)sigmas_07.size();isg++) {
	Print_To_File({},{F_T_d_MB_RE_sm_list[ixg][isg].distr_list[iens].distr, F_T_d_MB_IM_sm_list[ixg][isg].distr_list[iens].distr, F_T_d_MB_IM_Gauss_sm_list[ixg][isg].distr_list[iens].distr}, path_out+"/FF_d/"+Ens_tags[iens]+"/jackknife/ixg_"+to_string(ixg)+"/FF_MB_isg_"+to_string(isg)+".dat","",  "Njacks: "+to_string(Njacks));
	Print_To_File({},{F_T_d_MB_RE_sm_list_15[ixg][isg].distr_list[iens].distr, F_T_d_MB_IM_sm_list_15[ixg][isg].distr_list[iens].distr}, path_out+"/FF_d/"+Ens_tags[iens]+"/jackknife/ixg_"+to_string(ixg)+"/FF_15_MB_isg_"+to_string(isg)+".dat","", "Njacks: "+to_string(Njacks));
	Print_To_File({},{F_T_d_MB_RE_sm_list_7[ixg][isg].distr_list[iens].distr, F_T_d_MB_IM_sm_list_7[ixg][isg].distr_list[iens].distr}, path_out+"/FF_d/"+Ens_tags[iens]+"/jackknife/ixg_"+to_string(ixg)+"/FF_7_MB_isg_"+to_string(isg)+".dat","",  "Njacks: "+to_string(Njacks));
	Print_To_File({},{F_T_d_MB_RE_sm_list_no_sub[ixg][isg].distr_list[iens].distr}, path_out+"/FF_d/"+Ens_tags[iens]+"/jackknife/ixg_"+to_string(ixg)+"/FF_no_sub_MB_isg_"+to_string(isg)+".dat","",  "Njacks: "+to_string(Njacks));
	
      }
      

    }
    
    
    }
  }


  //read jackknifes

  for(int ixg=0;ixg<n_xg;ixg++) {
    for(int isg=0;isg<(signed)sigmas_07.size();isg++) {

      
      F_T_d_MB_RE_sm_list[ixg][isg] = Get_id_distr_list(Nens, Njacks, UseJack);
      F_T_d_MB_IM_sm_list[ixg][isg] = Get_id_distr_list(Nens, Njacks, UseJack);
      F_T_d_MB_IM_Gauss_sm_list[ixg][isg] = Get_id_distr_list(Nens, Njacks, UseJack);
      F_T_d_MB_RE_sm_list_15[ixg][isg] = Get_id_distr_list(Nens, Njacks, UseJack);
      F_T_d_MB_IM_sm_list_15[ixg][isg] = Get_id_distr_list(Nens, Njacks, UseJack);
      F_T_d_MB_RE_sm_list_7[ixg][isg] = Get_id_distr_list(Nens, Njacks, UseJack);
      F_T_d_MB_IM_sm_list_7[ixg][isg] = Get_id_distr_list(Nens, Njacks, UseJack);

      F_T_d_MB_RE_sm_list_no_sub[ixg][isg] = Get_id_distr_list(Nens, Njacks, UseJack);
      
      for(int iens=0;iens<Nens;iens++) {
	 
	 
	F_T_d_MB_RE_sm_list[ixg][isg].distr_list[iens].distr = Read_From_File( path_out+"/FF_d/"+Ens_tags[iens]+"/jackknife/ixg_"+to_string(ixg)+"/FF_MB_isg_"+to_string(isg)+".dat", 1, 4,1);
	F_T_d_MB_IM_sm_list[ixg][isg].distr_list[iens].distr = Read_From_File( path_out+"/FF_d/"+Ens_tags[iens]+"/jackknife/ixg_"+to_string(ixg)+"/FF_MB_isg_"+to_string(isg)+".dat", 2, 4,1);
	F_T_d_MB_IM_Gauss_sm_list[ixg][isg].distr_list[iens].distr = Read_From_File( path_out+"/FF_d/"+Ens_tags[iens]+"/jackknife/ixg_"+to_string(ixg)+"/FF_MB_isg_"+to_string(isg)+".dat", 3, 4,1);
	
	F_T_d_MB_RE_sm_list_15[ixg][isg].distr_list[iens].distr= Read_From_File( path_out+"/FF_d/"+Ens_tags[iens]+"/jackknife/ixg_"+to_string(ixg)+"/FF_15_MB_isg_"+to_string(isg)+".dat", 1, 3,1);
	F_T_d_MB_IM_sm_list_15[ixg][isg].distr_list[iens].distr= Read_From_File( path_out+"/FF_d/"+Ens_tags[iens]+"/jackknife/ixg_"+to_string(ixg)+"/FF_15_MB_isg_"+to_string(isg)+".dat", 2, 3,1);
	
	F_T_d_MB_RE_sm_list_7[ixg][isg].distr_list[iens].distr= Read_From_File( path_out+"/FF_d/"+Ens_tags[iens]+"/jackknife/ixg_"+to_string(ixg)+"/FF_7_MB_isg_"+to_string(isg)+".dat", 1, 3,1);
	F_T_d_MB_IM_sm_list_7[ixg][isg].distr_list[iens].distr= Read_From_File( path_out+"/FF_d/"+Ens_tags[iens]+"/jackknife/ixg_"+to_string(ixg)+"/FF_7_MB_isg_"+to_string(isg)+".dat", 2, 3,1);


	F_T_d_MB_RE_sm_list_no_sub[ixg][isg].distr_list[iens].distr = Read_From_File( path_out+"/FF_d/"+Ens_tags[iens]+"/jackknife/ixg_"+to_string(ixg)+"/FF_no_sub_MB_isg_"+to_string(isg)+".dat", 1, 2,1);
      }
    }
  }
    
  cout<<"Reading RE"<<endl;
  printV(F_T_d_MB_RE_sm_list[0][0].distr_list[0].distr, "RE", 1);
  cout<<"Reading_IM"<<endl,
  printV(F_T_d_MB_IM_sm_list[0][0].distr_list[0].distr, "IM", 1);
  //extrapolate the results for each sigma and xg to the continuum limit, employing either a constant or linear fit in a^2 (to be combined with BMA with uniform weight 1/2 )

  distr_t_list F_T_d_I_ori_extr(UseJack);
  distr_t_list F_T_d_I_sub_extr(UseJack);

  for(int ixg=0; ixg<n_xg;ixg++) {

    //################################################################################
    // b-quark contribution
    distr_t b_B64(UseJack), b_D96(UseJack) ;
    for(int iens=0;iens<Ens_tags.size();iens++) {
      if(Ens_tags[iens] == "cB211b.072.64") {
	b_B64= F_T_u_list[ixg][iens];
      }
      else b_D96= F_T_u_list[ixg][iens];
    }

         
    distr_t FT_b_const= (1/pow(b_B64.err(),2))*b_B64 + (1/pow(b_D96.err(),2))*b_D96;
    FT_b_const = FT_b_const/( (1/pow(b_B64.err(),2)) +  (1/pow(b_D96.err(),2)));
    double ch2_const = pow((FT_b_const.ave() - b_B64.ave())/b_B64.err(),2) +  pow( (FT_b_const.ave() - b_D96.ave())/b_D96.err(),2) ;
       
    distr_t FT_b_a2= b_D96*a_B*a_B - b_B64*a_D*a_D;
    FT_b_a2 = FT_b_a2/( a_B*a_B - a_D*a_D);


        
    //combine using Eq. 29
    distr_t FT_b(UseJack);
    if(ch2_const < 2) {
      FT_b  = BMA_Eq_29( {FT_b_const, FT_b_a2});
    }
    else FT_b= FT_b_a2;
    //################################################################################

    
   
    // s-quark contribution

     
      distr_t FT_I_ori(UseJack);

      distr_t s_I_ori_B64(UseJack), s_I_ori_D96(UseJack) ;
      for(int iens=0;iens<Ens_tags.size();iens++) {
	if(Ens_tags[iens] == "cB211b.072.64") {
	  s_I_ori_B64= F_T_d_I_ori_list[ixg][iens];
	}
	else {
	  s_I_ori_D96= F_T_d_I_ori_list[ixg][iens];
	}
      }

      distr_t FT_I_ori_const =  (1/pow(s_I_ori_B64.err(),2))*s_I_ori_B64 + (1/pow(s_I_ori_D96.err(),2))*s_I_ori_D96;
      FT_I_ori_const = FT_I_ori_const/(  (1/pow(s_I_ori_B64.err(),2)) + (1/pow(s_I_ori_D96.err(),2))   );
      distr_t FT_I_ori_a2= s_I_ori_D96*a_B*a_B - s_I_ori_B64*a_D*a_D;
      FT_I_ori_a2 = FT_I_ori_a2/( a_B*a_B - a_D*a_D);
      ch2_const = pow( (FT_I_ori_const.ave() - s_I_ori_B64.ave())/s_I_ori_B64.err(),2) + pow( (FT_I_ori_const.ave() - s_I_ori_D96.ave())/s_I_ori_D96.err(),2);
      //combine using Eq. 29
      if(ch2_const < 2) {
	FT_I_ori  = BMA_Eq_29( {FT_I_ori_const, FT_I_ori_a2});
      }
      else FT_I_ori= FT_I_ori_a2;

      F_T_d_I_ori_extr.distr_list.push_back(FT_I_ori);


      distr_t FT_I_sub(UseJack);

      distr_t s_I_sub_B64(UseJack), s_I_sub_D96(UseJack) ;
      for(int iens=0;iens<Ens_tags.size();iens++) {
	if(Ens_tags[iens] == "cB211b.072.64") {
	  s_I_sub_B64= F_T_d_I_list[ixg][iens];
	}
	else {
	  s_I_sub_D96= F_T_d_I_list[ixg][iens];
	}
      }

      distr_t FT_I_sub_const =  (1/pow(s_I_sub_B64.err(),2))*s_I_sub_B64 + (1/pow(s_I_sub_D96.err(),2))*s_I_sub_D96;
      FT_I_sub_const = FT_I_sub_const/(  (1/pow(s_I_sub_B64.err(),2)) + (1/pow(s_I_sub_D96.err(),2))   );
      distr_t FT_I_sub_a2= s_I_sub_D96*a_B*a_B - s_I_sub_B64*a_D*a_D;
      FT_I_sub_a2 = FT_I_sub_a2/( a_B*a_B - a_D*a_D);
      ch2_const = pow( (FT_I_sub_const.ave() - s_I_sub_B64.ave())/s_I_sub_B64.err(),2) + pow( (FT_I_sub_const.ave() - s_I_sub_D96.ave())/s_I_sub_D96.err(),2);
      //combine using Eq. 29
      if(ch2_const < 2) {
	FT_I_sub  = BMA_Eq_29( {FT_I_sub_const, FT_I_sub_a2});
      }
      else FT_I_sub= FT_I_sub_a2;

      FT_I_sub= FT_I_sub_a2;

      F_T_d_I_sub_extr.distr_list.push_back(FT_I_sub);
      
	

      distr_t_list FT_s_RE_sigma(UseJack), FT_s_IM_sigma(UseJack), FT_s_IM_Gauss_sigma(UseJack);
      distr_t_list FT_s_RE_sigma_15(UseJack), FT_s_IM_sigma_15(UseJack);
      distr_t_list FT_s_RE_sigma_7(UseJack), FT_s_IM_sigma_7(UseJack);
      distr_t_list FT_s_RE_sigma_no_sub(UseJack);
      distr_t_list FT_s_RE_sigma_no_sub_B64(UseJack);
      distr_t_list FT_s_RE_sigma_no_sub_D96(UseJack);

      
    for(int is=0;is<(signed)sigma_simulated[ixg].size();is++) {

     
      
	distr_t s_I_B64(UseJack), s_I_D96(UseJack) ;
	

	distr_t s_II_B64(UseJack), s_II_D96(UseJack);
	distr_t s_B64_IM(UseJack), s_D96_IM(UseJack);

	distr_t s_II_B64_no_sub(UseJack), s_II_D96_no_sub(UseJack);

	distr_t s_I_B64_15(UseJack), s_I_D96_15(UseJack) ;
	distr_t s_II_B64_15(UseJack), s_II_D96_15(UseJack);
	distr_t s_B64_IM_15(UseJack), s_D96_IM_15(UseJack);

	distr_t s_I_B64_7(UseJack), s_I_D96_7(UseJack) ;
	distr_t s_II_B64_7(UseJack), s_II_D96_7(UseJack);
	distr_t s_B64_IM_7(UseJack), s_D96_IM_7(UseJack);

	
	distr_t s_B64_IM_Gauss(UseJack), s_D96_IM_Gauss(UseJack);
	for(int iens=0;iens<Ens_tags.size();iens++) {
	  if(Ens_tags[iens] == "cB211b.072.64") {
	    s_I_B64= F_T_d_I_list[ixg][iens];
            s_II_B64= F_T_d_MB_RE_sm_list[ixg][is][iens];
	    s_B64_IM= F_T_d_MB_IM_sm_list[ixg][is][iens];
	    s_B64_IM_Gauss= F_T_d_MB_IM_Gauss_sm_list[ixg][is][iens];

	  	   
	    s_I_B64_15= F_T_d_I_list_15[ixg][iens];
            s_II_B64_15= F_T_d_MB_RE_sm_list_15[ixg][is][iens];
	    s_B64_IM_15= F_T_d_MB_IM_sm_list_15[ixg][is][iens];

	    s_I_B64_7= F_T_d_I_list_7[ixg][iens];
            s_II_B64_7= F_T_d_MB_RE_sm_list_7[ixg][is][iens];
	    s_B64_IM_7= F_T_d_MB_IM_sm_list_7[ixg][is][iens];

	    s_II_B64_no_sub = F_T_d_MB_RE_sm_list_no_sub[ixg][is][iens];

	    FT_s_RE_sigma_no_sub_B64.distr_list.push_back( s_II_B64_no_sub);
	    
	  }
	  else {
	    s_I_D96= F_T_d_I_list[ixg][iens];
	    s_II_D96= F_T_d_MB_RE_sm_list[ixg][is][iens];
	    s_D96_IM= F_T_d_MB_IM_sm_list[ixg][is][iens];
	    s_D96_IM_Gauss= F_T_d_MB_IM_Gauss_sm_list[ixg][is][iens];


	    s_I_D96_15= F_T_d_I_list_15[ixg][iens];
	    s_II_D96_15= F_T_d_MB_RE_sm_list_15[ixg][is][iens];
	    s_D96_IM_15= F_T_d_MB_IM_sm_list_15[ixg][is][iens];

	    s_I_D96_7= F_T_d_I_list_7[ixg][iens];
	    s_II_D96_7= F_T_d_MB_RE_sm_list_7[ixg][is][iens];
	    s_D96_IM_7= F_T_d_MB_IM_sm_list_7[ixg][is][iens];

	    s_II_D96_no_sub = F_T_d_MB_RE_sm_list_no_sub[ixg][is][iens];

	    FT_s_RE_sigma_no_sub_D96.distr_list.push_back( s_II_D96_no_sub);
	    
	    
	  }
	}

	//add sigma contribution
	distr_t s_B64_RE = s_I_B64 + s_II_B64;
	distr_t s_D96_RE = s_I_D96 + s_II_D96;


	distr_t s_B64_RE_15 = s_I_B64_15 + s_II_B64_15;
	distr_t s_D96_RE_15 = s_I_D96_15 + s_II_D96_15;
	
	distr_t s_B64_RE_7 = s_I_B64_7 + s_II_B64_7;
	distr_t s_D96_RE_7 = s_I_D96_7 + s_II_D96_7;

	distr_t s_B64_RE_no_sub= s_II_B64_no_sub;
	distr_t s_D96_RE_no_sub= s_II_D96_no_sub;
	
	//extrapolate


	distr_t FT_s_RE_const =  (1/pow(s_B64_RE.err(),2))*s_B64_RE + (1/pow(s_D96_RE.err(),2))*s_D96_RE;
	FT_s_RE_const = FT_s_RE_const/(  (1/pow(s_B64_RE.err(),2)) + (1/pow(s_D96_RE.err(),2))   );
	distr_t FT_s_RE_a2= s_D96_RE*a_B*a_B - s_B64_RE*a_D*a_D;
	FT_s_RE_a2 = FT_s_RE_a2/( a_B*a_B - a_D*a_D);
	ch2_const = pow( (FT_s_RE_const.ave() - s_B64_RE.ave())/s_B64_RE.err(),2) + pow( (FT_s_RE_const.ave() - s_D96_RE.ave())/s_D96_RE.err(),2);
	//combine using Eq. 29
	distr_t FT_s_RE(UseJack);
	if(ch2_const < 2) {
	  FT_s_RE  = BMA_Eq_29( {FT_s_RE_const, FT_s_RE_a2});
	}
	else FT_s_RE= FT_s_RE_a2;






	


	distr_t FT_s_RE_const_15 =  (1/pow(s_B64_RE_15.err(),2))*s_B64_RE_15 + (1/pow(s_D96_RE_15.err(),2))*s_D96_RE_15;
	FT_s_RE_const_15 = FT_s_RE_const_15/(  (1/pow(s_B64_RE_15.err(),2)) + (1/pow(s_D96_RE_15.err(),2))   );
	distr_t FT_s_RE_a2_15= s_D96_RE_15*a_B*a_B - s_B64_RE_15*a_D*a_D;
	FT_s_RE_a2_15 = FT_s_RE_a2_15/( a_B*a_B - a_D*a_D);
	ch2_const = pow( (FT_s_RE_const_15.ave() - s_B64_RE_15.ave())/s_B64_RE_15.err(),2) + pow( (FT_s_RE_const_15.ave() - s_D96_RE_15.ave())/s_D96_RE_15.err(),2);
	//combine using Eq. 29
	distr_t FT_s_RE_15(UseJack);
	if(ch2_const < 2) {
	  FT_s_RE_15  = BMA_Eq_29( {FT_s_RE_const_15, FT_s_RE_a2_15});
	}
	else FT_s_RE_15= FT_s_RE_a2_15;


	distr_t FT_s_RE_const_7 =  (1/pow(s_B64_RE_7.err(),2))*s_B64_RE_7 + (1/pow(s_D96_RE_7.err(),2))*s_D96_RE_7;
	FT_s_RE_const_7 = FT_s_RE_const_7/(  (1/pow(s_B64_RE_7.err(),2)) + (1/pow(s_D96_RE_7.err(),2))   );
	distr_t FT_s_RE_a2_7= s_D96_RE_7*a_B*a_B - s_B64_RE_7*a_D*a_D;
	FT_s_RE_a2_7 = FT_s_RE_a2_7/( a_B*a_B - a_D*a_D);
	ch2_const = pow( (FT_s_RE_const_7.ave() - s_B64_RE_7.ave())/s_B64_RE_7.err(),2) + pow( (FT_s_RE_const_7.ave() - s_D96_RE_7.ave())/s_D96_RE_7.err(),2);
	//combine using Eq. 29
	distr_t FT_s_RE_7(UseJack);
	if(ch2_const < 2) {
	  FT_s_RE_7  = BMA_Eq_29( {FT_s_RE_const_7, FT_s_RE_a2_7});
	}
	else FT_s_RE_7= FT_s_RE_a2_7;


	//no sub
 
	distr_t FT_s_RE_const_no_sub =  (1/pow(s_B64_RE_no_sub.err(),2))*s_B64_RE_no_sub + (1/pow(s_D96_RE_no_sub.err(),2))*s_D96_RE_no_sub;
	FT_s_RE_const_no_sub = FT_s_RE_const_no_sub/(  (1/pow(s_B64_RE_no_sub.err(),2)) + (1/pow(s_D96_RE_no_sub.err(),2))   );
	distr_t FT_s_RE_a2_no_sub= s_D96_RE_no_sub*a_B*a_B - s_B64_RE_no_sub*a_D*a_D;
	FT_s_RE_a2_no_sub = FT_s_RE_a2_no_sub/( a_B*a_B - a_D*a_D);
	ch2_const = pow( (FT_s_RE_const_no_sub.ave() - s_B64_RE_no_sub.ave())/s_B64_RE_no_sub.err(),2) + pow( (FT_s_RE_const_no_sub.ave() - s_D96_RE_no_sub.ave())/s_D96_RE_no_sub.err(),2);
	//combine using Eq. 29
	distr_t FT_s_RE_no_sub(UseJack);
	if(ch2_const < 2) {
	  FT_s_RE_no_sub  = BMA_Eq_29( {FT_s_RE_const_no_sub, FT_s_RE_a2_no_sub});
	}
	else FT_s_RE_no_sub= FT_s_RE_a2_no_sub;

	


	distr_t FT_s_IM_const =  (1/pow(s_B64_IM.err(),2))*s_B64_IM + (1/pow(s_D96_IM.err(),2))*s_D96_IM;
	FT_s_IM_const= FT_s_IM_const/( (1/pow(s_B64_IM.err(),2)) +  (1/pow(s_D96_IM.err(),2)));
	distr_t FT_s_IM_a2= s_D96_IM*a_B*a_B - s_B64_IM*a_D*a_D;
	FT_s_IM_a2 = FT_s_IM_a2/( a_B*a_B - a_D*a_D);

	ch2_const = pow( (FT_s_IM_const.ave() - s_B64_IM.ave())/s_B64_IM.err(),2) + pow( (FT_s_IM_const.ave() - s_D96_IM.ave())/s_D96_IM.err(),2);
	//combine using Eq. 29
	distr_t FT_s_IM(UseJack);
	if(ch2_const < 2) {
	  FT_s_IM  = BMA_Eq_29( {FT_s_IM_const, FT_s_IM_a2});
	}
	else FT_s_IM= FT_s_IM_a2;



	distr_t FT_s_IM_const_15 =  (1/pow(s_B64_IM_15.err(),2))*s_B64_IM_15 + (1/pow(s_D96_IM_15.err(),2))*s_D96_IM_15;
	FT_s_IM_const_15= FT_s_IM_const_15/( (1/pow(s_B64_IM_15.err(),2)) +  (1/pow(s_D96_IM_15.err(),2)));
	distr_t FT_s_IM_a2_15= s_D96_IM_15*a_B*a_B - s_B64_IM_15*a_D*a_D;
	FT_s_IM_a2_15 = FT_s_IM_a2_15/( a_B*a_B - a_D*a_D);

	ch2_const = pow( (FT_s_IM_const_15.ave() - s_B64_IM_15.ave())/s_B64_IM_15.err(),2) + pow( (FT_s_IM_const_15.ave() - s_D96_IM_15.ave())/s_D96_IM_15.err(),2);
	//combine using Eq. 29
	distr_t FT_s_IM_15(UseJack);
	if(ch2_const < 2) {
	  FT_s_IM_15  = BMA_Eq_29( {FT_s_IM_const_15, FT_s_IM_a2_15});
	}
	else FT_s_IM_15= FT_s_IM_a2_15;


	
	distr_t FT_s_IM_const_7 =  (1/pow(s_B64_IM_7.err(),2))*s_B64_IM_7 + (1/pow(s_D96_IM_7.err(),2))*s_D96_IM_7;
	FT_s_IM_const_7= FT_s_IM_const_7/( (1/pow(s_B64_IM_7.err(),2)) +  (1/pow(s_D96_IM_7.err(),2)));
	distr_t FT_s_IM_a2_7= s_D96_IM_7*a_B*a_B - s_B64_IM_7*a_D*a_D;
	FT_s_IM_a2_7 = FT_s_IM_a2_7/( a_B*a_B - a_D*a_D);

	ch2_const = pow( (FT_s_IM_const_7.ave() - s_B64_IM_7.ave())/s_B64_IM_7.err(),2) + pow( (FT_s_IM_const_7.ave() - s_D96_IM_7.ave())/s_D96_IM_7.err(),2);
	//combine using Eq. 29
	distr_t FT_s_IM_7(UseJack);
	if(ch2_const < 2) {
	  FT_s_IM_7  = BMA_Eq_29( {FT_s_IM_const_7, FT_s_IM_a2_7});
	}
	else FT_s_IM_7= FT_s_IM_a2_7;



	//Gauss

	distr_t FT_s_IM_Gauss_const =  (1/pow(s_B64_IM_Gauss.err(),2))*s_B64_IM_Gauss + (1/pow(s_D96_IM_Gauss.err(),2))*s_D96_IM_Gauss;
	FT_s_IM_Gauss_const= FT_s_IM_Gauss_const/( (1/pow(s_B64_IM_Gauss.err(),2)) +  (1/pow(s_D96_IM_Gauss.err(),2)));
	distr_t FT_s_IM_Gauss_a2= s_D96_IM_Gauss*a_B*a_B - s_B64_IM_Gauss*a_D*a_D;
	FT_s_IM_Gauss_a2 = FT_s_IM_Gauss_a2/( a_B*a_B - a_D*a_D);

	ch2_const = pow( (FT_s_IM_Gauss_const.ave() - s_B64_IM_Gauss.ave())/s_B64_IM_Gauss.err(),2) + pow( (FT_s_IM_Gauss_const.ave() - s_D96_IM_Gauss.ave())/s_D96_IM_Gauss.err(),2);
	//combine using Eq. 29
	distr_t FT_s_IM_Gauss(UseJack);
	if(ch2_const < 2) {
	  FT_s_IM_Gauss  = BMA_Eq_29( {FT_s_IM_Gauss_const, FT_s_IM_Gauss_a2});
	}
	else FT_s_IM_Gauss= FT_s_IM_Gauss_a2;

	

	FT_s_RE_sigma.distr_list.push_back(FT_s_RE);
	FT_s_IM_sigma.distr_list.push_back(FT_s_IM);
	FT_s_IM_Gauss_sigma.distr_list.push_back(FT_s_IM_Gauss);


	FT_s_RE_sigma_15.distr_list.push_back(FT_s_RE_15);
	FT_s_IM_sigma_15.distr_list.push_back(FT_s_IM_15);

	FT_s_RE_sigma_7.distr_list.push_back(FT_s_RE_7);
	FT_s_IM_sigma_7.distr_list.push_back(FT_s_IM_7);


	FT_s_RE_sigma_no_sub.distr_list.push_back(FT_s_RE_no_sub);

	

	

    }


    //print result for nosub

    Print_To_File({}, {sigma_simulated[ixg], FT_s_RE_sigma_no_sub_B64.ave(), FT_s_RE_sigma_no_sub_B64.err(),  FT_s_RE_sigma_no_sub_D96.ave(), FT_s_RE_sigma_no_sub_D96.err() , FT_s_RE_sigma_no_sub.ave(), FT_s_RE_sigma_no_sub.err()}, path_out+"/FF_d_extr/FT_s_sub_ixg_"+to_string(ixg)+".dat", "", "");

  


    //extrapolate to vanishing sigma
    class ipar_07 {
    public:
      ipar_07() : FF_RE(0.0), FF_RE_err(0.0), FF_IM(0.0), FF_IM_err(0.0), Is_Gauss(false), Is_five_finest(true) {}
      double FF_RE, FF_RE_err, FF_IM, FF_IM_err, sigma;
      bool Is_Gauss;
      bool Is_five_finest;
    };
  
    class fpar_07 {
    public:
      fpar_07() {}
      fpar_07(const Vfloat &par) {
	if((signed)par.size() != 5) crash("In class fpar_07  class constructor Vfloat par has size != 5");
	R=par[0];
	A=par[1];
	B=par[2];
	C=par[3];
	D=par[4];
	
      }
      double R,A,B,C,D;
    };
    
    int sigma_to_fit=sigmas_07.size();
    //init bootstrap fit
    bootstrap_fit<fpar_07,ipar_07> bf_07(Njacks);
    bf_07.set_warmup_lev(1); //sets warmup
    bf_07.Set_number_of_measurements(2*sigma_to_fit);
    bf_07.Set_verbosity(1);

    bf_07.Add_par("R", -0.05, 0.1);
    bf_07.Add_par("A", 1.0, 0.1);
    bf_07.Add_par("B", 1.0 , 0.1);
    bf_07.Add_par("C", 1.0 , 0.1);
    bf_07.Add_par("D", 1.0 , 0.1);
   
  
    //fit on mean values to get ch2
    bootstrap_fit<fpar_07,ipar_07> bf_07_ch2(1);
    bf_07_ch2.set_warmup_lev(1); //sets warmup
    bf_07_ch2.Set_number_of_measurements(2*sigma_to_fit);
    bf_07_ch2.Set_verbosity(1);
    bf_07_ch2.Add_par("R", -0.05, 0.1);
    bf_07_ch2.Add_par("A", 1.0, 0.1);
    bf_07_ch2.Add_par("B", 1.0, 0.1);
    bf_07_ch2.Add_par("C", 1.0, 0.1);
    bf_07_ch2.Add_par("D", 1.0, 0.1);


    bool Is_RE=true;

    if(MESON=="B3s") {
      bf_07.Fix_par("B",0.0);
      bf_07_ch2.Fix_par("B",0.0);
    }

    bf_07.Fix_par("C",0.0);
    bf_07_ch2.Fix_par("C",0.0);

    bf_07.Fix_par("D",0.0);
    bf_07_ch2.Fix_par("D",0.0);


    bool Use_five_finest=false;

    //ansatz
    bf_07.ansatz=  [&Is_RE , &Use_five_finest](const fpar_07 &p, const ipar_07 &ip) {

      if( (Use_five_finest==true) && (ip.Is_five_finest == false) ) return 0.0;

      if(Is_RE) {
	if(ip.Is_Gauss==false)
	  return p.R + p.A*(ip.sigma) + p.B*pow(ip.sigma,2);
	else return 0.0;
      }
      else {

	if(ip.Is_Gauss==false)
	  return p.R + p.A*(ip.sigma) + p.B*pow(ip.sigma,2);
	
	else return 0.0; //uncomment for global Gauss-Cauchy fit  p.R + p.C*pow(ip.sigma,2) + p.D*pow(ip.sigma,4);

      }

      return 0.0;
    	 	  
    };

  
    bf_07.measurement=  [&Is_RE, &Use_five_finest](const fpar_07 &p, const ipar_07 &ip) {

      if( (Use_five_finest==true) && (ip.Is_five_finest == false) ) return 0.0;
      
      if(Is_RE) return ip.FF_RE;

      if(ip.Is_Gauss==true) return 0.0; //to uncomment for global Gauss-Cauchy fit

      return ip.FF_IM;
      
    };
    
    bf_07.error=  [&Is_RE ](const fpar_07 &p, const ipar_07 &ip) {

      if(Is_RE) return ip.FF_RE_err;

      return ip.FF_IM_err;
    
  };
	
  bf_07_ch2.ansatz= bf_07.ansatz;
  bf_07_ch2.measurement = bf_07.measurement;
  bf_07_ch2.error = bf_07.error;

 
  //start fitting
  //fill the data
  vector<vector<ipar_07>> data_07(Njacks);
  vector<vector<ipar_07>> data_07_ch2(1);
  //allocate space for output result
  boot_fit_data<fpar_07> Bt_fit_07_RE;
  boot_fit_data<fpar_07> Bt_fit_07_RE_ch2;
  boot_fit_data<fpar_07> Bt_fit_07_IM;
  boot_fit_data<fpar_07> Bt_fit_07_IM_ch2;

  boot_fit_data<fpar_07> Bt_fit_07_RE_red;
  boot_fit_data<fpar_07> Bt_fit_07_RE_red_ch2;
  boot_fit_data<fpar_07> Bt_fit_07_IM_red;
  boot_fit_data<fpar_07> Bt_fit_07_IM_red_ch2;
  
  for(auto &data_iboot: data_07) data_iboot.resize(2*sigma_to_fit);
  for(auto &data_iboot: data_07_ch2) data_iboot.resize(2*sigma_to_fit);
  for(int ijack=0;ijack<Njacks;ijack++) {
    for(int is=0;is<sigma_to_fit;is++) {
      data_07[ijack][is].FF_RE = (FT_s_RE_sigma[is]).distr[ijack];
      data_07[ijack][is].FF_RE_err= (FT_s_RE_sigma[is]).err();
      data_07[ijack][is].FF_IM = (FT_s_IM_sigma[is]).distr[ijack];
      data_07[ijack][is].FF_IM_err= (FT_s_IM_sigma[is]).err();
      data_07[ijack][is].sigma = sigma_simulated[ixg][is];
      data_07[ijack][is].Is_Gauss = false;
      if(is < 5) { data_07[ijack][is].Is_five_finest= true;}
      else data_07[ijack][is].Is_five_finest= false;


      data_07[ijack][is+sigma_to_fit].FF_RE = 0.0;
      data_07[ijack][is+sigma_to_fit].FF_RE_err= 1.0;
      data_07[ijack][is+sigma_to_fit].FF_IM = (FT_s_IM_Gauss_sigma[is]).distr[ijack];
      data_07[ijack][is+sigma_to_fit].FF_IM_err= (FT_s_IM_Gauss_sigma[is]).err() + 1e-4;
      data_07[ijack][is+sigma_to_fit].sigma = sigma_simulated[ixg][is];
      data_07[ijack][is+sigma_to_fit].Is_Gauss = true;
      if(is < 5) { data_07[ijack][is+sigma_to_fit].Is_five_finest= true;}
      else data_07[ijack][is+sigma_to_fit].Is_five_finest= false;
      
      if(ijack==0) {
	data_07_ch2[ijack][is].FF_RE = (FT_s_RE_sigma[is]).ave();
	data_07_ch2[ijack][is].FF_RE_err= (FT_s_RE_sigma[is]).err();
	data_07_ch2[ijack][is].FF_IM = (FT_s_IM_sigma[is]).ave();
	data_07_ch2[ijack][is].FF_IM_err= (FT_s_IM_sigma[is]).err();
	data_07_ch2[ijack][is].sigma = sigma_simulated[ixg][is];
	data_07_ch2[ijack][is].Is_Gauss = false;
	if(is < 5) { data_07_ch2[ijack][is].Is_five_finest= true;}
	else data_07_ch2[ijack][is].Is_five_finest= false;


	data_07_ch2[ijack][is+sigma_to_fit].FF_RE = 0.0;
	data_07_ch2[ijack][is+sigma_to_fit].FF_RE_err= 1.0;
	data_07_ch2[ijack][is+sigma_to_fit].FF_IM = (FT_s_IM_Gauss_sigma[is]).ave();
	data_07_ch2[ijack][is+sigma_to_fit].FF_IM_err= (FT_s_IM_Gauss_sigma[is]).err() + 1e-4;
	data_07_ch2[ijack][is+sigma_to_fit].sigma = sigma_simulated[ixg][is];
	data_07_ch2[ijack][is+sigma_to_fit].Is_Gauss = true;
	if(is < 5) { data_07_ch2[ijack][is+sigma_to_fit].Is_five_finest= true;}
	else data_07_ch2[ijack][is+sigma_to_fit].Is_five_finest= false;
	
      }
    }
  }
  
  //append
  bf_07.Append_to_input_par(data_07);
  bf_07_ch2.Append_to_input_par(data_07_ch2);
  //fit
  cout<<"Fitting real part of the form factor for ixg: "<<ixg<<endl;



  Bt_fit_07_RE= bf_07.Perform_bootstrap_fit();
  Bt_fit_07_RE_ch2= bf_07_ch2.Perform_bootstrap_fit();
  Use_five_finest=true;
  bf_07.Fix_par("B",0.0);
  bf_07_ch2.Fix_par("B",0.0);
  bf_07.Set_par_val("A", -1.0, 0.03);
  bf_07_ch2.Set_par_val("A", -1.0, 0.03);
  Bt_fit_07_RE_red= bf_07.Perform_bootstrap_fit();
  Bt_fit_07_RE_red_ch2= bf_07_ch2.Perform_bootstrap_fit();
  Use_five_finest=false;
  bf_07.Release_par("B");
  bf_07_ch2.Release_par("B");
  
  double ch2_red_RE= Bt_fit_07_RE_ch2.get_ch2_ave()/( sigma_to_fit - bf_07.Get_number_of_fit_pars());
  cout<<"Reduced ch2: "<<ch2_red_RE<<endl;
  Is_RE=false;
  //bf_07.Fix_par("B",0.0);
  //bf_07_ch2.Fix_par("B",0.0);
  bf_07.Release_par("C");
  bf_07_ch2.Release_par("C");
  bf_07.Set_par_val("C", 1.0);
  bf_07_ch2.Set_par_val("C",1.0);

  bf_07.Release_par("D");
  bf_07_ch2.Release_par("D");
  bf_07.Set_par_val("D", 1.0);
  bf_07_ch2.Set_par_val("D",1.0);
  
  bf_07.Set_par_val("R", 0.05, 0.01);
  bf_07_ch2.Set_par_val("R", 0.05,0.01);
  

  if(MESON=="B3s") {
    bf_07.Fix_par("B",0.0);
    bf_07_ch2.Fix_par("B",0.0);
  }
  
  if(ixg==3) { bf_07.Fix_par("B",0.0); bf_07_ch2.Fix_par("B",0.0);}
  if(ixg > 0 && MESON != "B3s") { bf_07.Fix_par("B",0.0); bf_07_ch2.Fix_par("B",0.0);}
  
  cout<<"Fitting imag part of the form factor for ixg: "<<ixg<<endl;
  Bt_fit_07_IM= bf_07.Perform_bootstrap_fit();
  Bt_fit_07_IM_ch2= bf_07_ch2.Perform_bootstrap_fit();
  double ch2_red_IM= Bt_fit_07_IM_ch2.get_ch2_ave()/( 2*sigma_to_fit - bf_07.Get_number_of_fit_pars());
  cout<<"Reduced ch2: "<<ch2_red_IM<<endl;

  Use_five_finest=true;
  bf_07.Fix_par("B",0.0);
  bf_07_ch2.Fix_par("B",0.0);
  Bt_fit_07_IM_red= bf_07.Perform_bootstrap_fit();
  Bt_fit_07_IM_red_ch2= bf_07_ch2.Perform_bootstrap_fit();
  Use_five_finest=false;
  bf_07.Release_par("B");
  bf_07_ch2.Release_par("B");
  
  



  //retrieve parameters

  distr_t R_RE(UseJack), R_IM(UseJack), A_RE(UseJack), A_IM(UseJack), B_RE(UseJack), B_IM(UseJack), C_IM(UseJack), D_IM(UseJack);

  distr_t R_RE_red(UseJack), R_IM_red(UseJack);

  for(int ijack=0;ijack<Njacks;ijack++) {
    R_RE.distr.push_back( Bt_fit_07_RE.par[ijack].R);
    A_RE.distr.push_back( Bt_fit_07_RE.par[ijack].A);
    B_RE.distr.push_back( Bt_fit_07_RE.par[ijack].B);

    R_IM.distr.push_back( Bt_fit_07_IM.par[ijack].R);
    A_IM.distr.push_back( Bt_fit_07_IM.par[ijack].A);
    B_IM.distr.push_back( Bt_fit_07_IM.par[ijack].B);
    C_IM.distr.push_back( Bt_fit_07_IM.par[ijack].C);
    D_IM.distr.push_back( Bt_fit_07_IM.par[ijack].D);

    R_RE_red.distr.push_back( Bt_fit_07_RE_red.par[ijack].R);
    R_IM_red.distr.push_back( Bt_fit_07_IM_red.par[ijack].R);
    
  }

  //get fit function

  Vfloat sigmas_to_print;
  Vfloat sigmas_to_print_0;
  int Nsigmas_to_print=400;
  if(MESON== "B0s") sigmas_to_print_0.push_back( -0.03);
  if(MESON== "B1s") sigmas_to_print_0.push_back( -0.02);
  if(MESON== "B3s") sigmas_to_print_0.push_back( -0.01);

  
  
  for(int n=0;n<Nsigmas_to_print;n++) sigmas_to_print.push_back( n*6.0/(Nsigmas_to_print-1));

  distr_t_list F_print_RE(UseJack, Nsigmas_to_print, Njacks), F_print_IM(UseJack, Nsigmas_to_print, Njacks), F_print_IM_Gauss(UseJack,Nsigmas_to_print, Njacks);

  distr_t_list F_print_RE_0(UseJack), F_print_IM_0(UseJack);

  double sRE= fabs( R_RE.ave() - R_RE_red.ave())*erf(  fabs(R_RE.ave() - R_RE_red.ave())/(sqrt(2.0)*R_RE_red.err())   );
  double sIM= fabs( R_IM.ave() - R_IM_red.ave())*erf(  fabs(R_IM.ave() - R_IM_red.ave())/(sqrt(2.0)*R_IM_red.err())   );

   
  F_print_RE_0.distr_list.push_back(  R_RE.ave() + (R_RE -R_RE.ave())*sqrt( pow(R_RE.err(),2) + pow(sRE,2))/R_RE.err());
  F_print_IM_0.distr_list.push_back(  R_IM.ave() + (R_IM -R_IM.ave())*sqrt( pow(R_IM.err(),2) + pow(sIM,2))/R_IM.err());
    
  for(int n=0;n<Nsigmas_to_print;n++) {

    ipar_07 fS;
    fS.sigma= sigmas_to_print[n];
   
    for(int ijack=0;ijack<Njacks;ijack++) {
      Is_RE=true;
      fS.Is_Gauss=false;
      F_print_RE.distr_list[n].distr[ijack] = bf_07.ansatz( Bt_fit_07_RE.par[ijack], fS);
      Is_RE=false;
      F_print_IM.distr_list[n].distr[ijack] = bf_07.ansatz( Bt_fit_07_IM.par[ijack], fS);
      fS.Is_Gauss=true;
      F_print_IM_Gauss.distr_list[n].distr[ijack] = bf_07.ansatz( Bt_fit_07_IM.par[ijack], fS);

    }
    
  }



  //########################## COMBINED   MODEL-DEPENDENT EXTRAPOLATION #############################################
  //#################################################################################################################


    //extrapolate to vanishing sigma
    class ipar_07_comb {
    public:
      ipar_07_comb() : FF(0.0), FF_err(0), Is_RE(0), G1(0.0), G2(0.0), G3(0.0) {}
      double FF, FF_err, sigma, I_ori;
      bool Is_RE;
      double G1, G2, G3;
     
    };
  
    class fpar_07_comb {
    public:
      fpar_07_comb() {}
      fpar_07_comb(const Vfloat &par) {
	if((signed)par.size() != 5) crash("In class fpar_07_comb  class constructor Vfloat par has size != 5");

	TH=par[0];
	A=par[1];
	//priori pars
	A1=par[2];
	A2=par[3];
	A3=par[4];
      }
      double TH,A, A1, A2, A3;
    };

    bool fix_TH=false;
    bool Fit_I_TO=false;
    
    //init bootstrap fit
    bootstrap_fit<fpar_07_comb,ipar_07_comb> bf_07_comb(Njacks);
    bf_07_comb.set_warmup_lev(1); //sets warmup
    bf_07_comb.Set_number_of_measurements(2*sigma_to_fit);
    bf_07_comb.Set_verbosity(1);

    bf_07_comb.Add_par("TH", 0.3, 0.05);
    bf_07_comb.Add_par("A", -0.003, 0.001);
    bf_07_comb.Add_par("A1", 0.05, 0.005);
    bf_07_comb.Add_par("A2", 0.03, 0.003);
    bf_07_comb.Add_par("A3", 0.003, 0.001);
    bf_07_comb.Add_prior_par("TH", 0.5, 0.05);
    bf_07_comb.Add_prior_par("A", -0.003, 0.001);
    bf_07_comb.Add_prior_par("A1", 0.05, 0.005);
    bf_07_comb.Add_prior_par("A2", 0.03, 0.003);
    bf_07_comb.Add_prior_par("A3", 0.003, 0.001);

   
    

   
  
    //fit on mean values to get ch2
    bootstrap_fit<fpar_07_comb,ipar_07_comb> bf_07_comb_ch2(1);
    bf_07_comb_ch2.set_warmup_lev(1); //sets warmup
    bf_07_comb_ch2.Set_number_of_measurements(2*sigma_to_fit);
    bf_07_comb_ch2.Set_verbosity(1);

    bf_07_comb_ch2.Add_par("TH", 0.3, 0.05);
    bf_07_comb_ch2.Add_par("A", -0.003, 0.001);
    bf_07_comb_ch2.Add_par("A1", 0.05, 0.005);
    bf_07_comb_ch2.Add_par("A2", 0.03, 0.003);
    bf_07_comb_ch2.Add_par("A3", 0.003, 0.001);

    bf_07_comb_ch2.Add_prior_par("TH", 0.5, 0.05);
    bf_07_comb_ch2.Add_prior_par("A", -0.003, 0.001);
    bf_07_comb_ch2.Add_prior_par("A1", 0.05, 0.005);
    bf_07_comb_ch2.Add_prior_par("A2", 0.03, 0.003);
    bf_07_comb_ch2.Add_prior_par("A3", 0.003, 0.001);


    if(fix_TH) {

      bf_07_comb.Fix_par("TH",1e3);
      bf_07_comb_ch2.Fix_par("TH",1e3);
      
    }

    
    if(!Fit_I_TO) {
      bf_07_comb.Fix_par("A",0.0);
      bf_07_comb_ch2.Fix_par("A",0.0);
    }
    
      

   
    double Mp_cont;
    if(MESON=="B0s") Mp_cont= masses[0];
    if(MESON=="B1s") Mp_cont= masses[1];
    if(MESON=="B3s") Mp_cont= masses[2];
    double xg= 0.1*(ixg+1);
    double Et= (y_eff[ixg]*Mp_cont + (1-y_eff[ixg])*MBs)*(1.0 -0.5*xg);
    double Eph1= sqrt(pow(1.019,2) + pow(0.5*xg*Mp_cont,2));
    double Eph2= sqrt(pow(1.680,2) + pow(0.5*xg*Mp_cont,2));
    double Eph3= sqrt(pow(2.163,2) + pow(0.5*xg*Mp_cont,2));
    distr_t G1_distr(UseJack), G2_distr(UseJack), G3_distr(UseJack);
    double G1_ave= 0.004249/2.0; double G1_err= 0.000013/2.0;
    double G2_ave= 0.150/2.0;    double G2_err= 0.050/2.0;
    double G3_ave= 0.103/2.0;    double G3_err=0.025/2.0;

    for(int ijack=0;ijack<Njacks;ijack++) {
      G1_distr.distr.push_back( G1_ave+ GM_07()*G1_err/sqrt(Njacks-1.0));
      G2_distr.distr.push_back( G2_ave+ GM_07()*G2_err/sqrt(Njacks-1.0));
      G3_distr.distr.push_back( G3_ave+ GM_07()*G3_err/sqrt(Njacks-1.0));
    }
      

    string MODE_SPECTR="";
    
   
    //ansatz
    bf_07_comb.ansatz =  [&Et, &Eph1, &Eph2, &Eph3,  &MODE_SPECTR, &fix_TH](const fpar_07_comb &p, const ipar_07_comb &ip) {

      double TH= Eph3 + sqrt(p.TH*p.TH);
      Vfloat As({p.A1, p.A2, p.A3});
      Vfloat Gs({ip.G1, ip.G2, ip.G3});
      Vfloat Ephs({Eph1,Eph2,Eph3});
      
      //coupling D #####################
      double D=  TH*( p.A1*(ip.G1)/(pow(TH-Eph1,2) + pow(ip.G1,2)) +  p.A2*(ip.G2)/(pow(TH-Eph2,2) + pow(ip.G2,2)) +p.A3*(ip.G3)/(pow(TH-Eph3,2) + pow(ip.G3,2)) );
      if(MODE_SPECTR == "sqrt")  D= sqrt(TH)*( p.A1*(ip.G1)/(pow(TH-Eph1,2) + pow(ip.G1,2)) +  p.A2*(ip.G2)/(pow(TH-Eph2,2) + pow(ip.G2,2)) +p.A3*(ip.G3)/(pow(TH-Eph3,2) + pow(ip.G3,2)) );
      if(MODE_SPECTR == "quad")  D= pow(TH,2)*( p.A1*(ip.G1)/(pow(TH-Eph1,2) + pow(ip.G1,2)) +  p.A2*(ip.G2)/(pow(TH-Eph2,2) + pow(ip.G2,2)) +p.A3*(ip.G3)/(pow(TH-Eph3,2) + pow(ip.G3,2)) );
      //###############################
      
      //######## real part ###############
      if(ip.Is_RE) {

	double x_res_ori= -p.A1*(Et-Eph1)/(pow(Et-Eph1,2) + pow(ip.G1+ip.sigma,2)) +  -p.A2*(Et-Eph2)/(pow(Et-Eph2,2) + pow(ip.G2+ip.sigma,2)) +  -p.A3*(Et-Eph3)/(pow(Et-Eph3,2) + pow(ip.G3+ip.sigma,2)) ;
	double x_res_0_ori= -p.A1*(0.0-Eph1)/(pow(0.0-Eph1,2) + pow(ip.G1,2)) +  -p.A2*(0.0-Eph2)/(pow(0.0-Eph2,2) + pow(ip.G2,2)) +  -p.A3*(0.0-Eph3)/(pow(0.0-Eph3,2) + pow(ip.G3,2)) ;

	double x_phi= -p.A1*(Et-Eph1)/(pow(Et-Eph1,2) + pow(ip.G1+ip.sigma,2));
	double x_phi_1680= -p.A2*(Et-Eph2)/(pow(Et-Eph2,2) + pow(ip.G2+ip.sigma,2));
	double x_phi_2200= -p.A3*(Et-Eph3)/(pow(Et-Eph3,2) + pow(ip.G3+ip.sigma,2));



	double x_res=0;
	double x_res_0=0;

	complex I(0.0,1.0);
	double e=ip.sigma;

	for(int n=0;n<(signed)As.size();n++) {

	  complex r= ((Ephs[n]-I*e -Et + ((TH > Et)?2.0*I*Gs[n]:0.0))*M_PI + 2.0*(-Ephs[n] + I*e + Et)*atan( (Ephs[n] -TH)/Gs[n]) + Gs[n]*( 2.0*I*atan( e/(Et-TH)) -log( pow(Ephs[n]-TH,2) + pow(Gs[n],2)) + log( pow(e,2) + pow(Et-TH,2))))/( (pow(Gs[n],2) + (Et -Ephs[n] +I*e)*(Et-Ephs[n]+I*e)));

	  complex r0 = ((Ephs[n] + 2.0*I*Gs[n])*M_PI -2.0*I*Gs[n]*atan(e/TH) -2.0*(Ephs[n])*atan( (Ephs[n]-TH)/Gs[n]) + Gs[n]*log( ( pow(TH,2))/( pow(Gs[n],2) + pow(TH-Ephs[n],2))))/( pow(Gs[n],2) + (Ephs[n])*(Ephs[n]));

	  x_res += real(r)*As[n]/(2.0*M_PI);

	  //if(n==0) cout<<"RE phi: original: "<<x_phi<<" , new: "<<real(r)*As[n]/(2.0*M_PI)<<endl;
	  //if(n==1) cout<<"RE phi(1680): original: "<<x_phi_1680<<" , new: "<<real(r)*As[n]/(2.0*M_PI)<<endl;
	  //if(n==2) cout<<"RE phi(2200): original: "<<x_phi_2200<<" , new: "<<real(r)*As[n]/(2.0*M_PI)<<endl;
	  
	  x_res_0 += real(r0)*As[n]/(2.0*M_PI);

	  //if(n==0) cout<<"RE0 phi: original: "<<-p.A1*(0.0-Eph1)/(pow(0.0-Eph1,2) + pow(ip.G1,2))<<" , new: "<<real(r0)*As[n]/(2.0*M_PI)<<endl;
	  //if(n==1) cout<<"RE0 phi(1680): "<<real(r0)*As[n]/(2.0*M_PI)<<endl;
	  //if(n==2) cout<<"RE0 phi(2200): "<<real(r0)*As[n]/(2.0*M_PI)<<endl;

	}

	//cout<<"x_res(RE) original: "<<x_res_ori<<" , new: "<<x_res<<endl;
	//cout<<"x_res0(RE) original: "<<x_res_0_ori<<" , new: "<<x_res_0<<endl;


        complex z(Et,ip.sigma);
	double x_cont= real( -D*log(1.0 - z/TH)/z)/M_PI;
	double x_cont_0 = D/(M_PI*TH);

	//if(fix_TH)  return p.A + x_res;

	if(MODE_SPECTR=="sqrt") {
	  x_cont= real( 2*D*atanh(sqrt(z/TH))/sqrt(z))/M_PI;
	  x_cont_0 = 2*D/(M_PI*TH);
	}
	if(MODE_SPECTR=="quad") {
	  x_cont = real( -D*log(1.0-z/TH)/(z*z)  -D/(TH*z))/M_PI;
	  x_cont_0 = D/(M_PI*2*TH*TH);
	}

	return (p.A+ ip.I_ori) +  x_res + x_cont -x_cont_0 -x_res_0;
      }
      //#################################################################################################

      
      
      double x_res_ori= p.A1*(ip.G1+ip.sigma)/(pow(Et-Eph1,2) + pow(ip.G1+ip.sigma,2)) +  p.A2*(ip.G2+ip.sigma)/(pow(Et-Eph2,2) + pow(ip.G2+ip.sigma,2)) + p.A3*(ip.G3+ip.sigma)/(pow(Et-Eph3,2) + pow(ip.G3+ip.sigma,2))  ;


      double x_phi_IM= p.A1*(ip.G1+ip.sigma)/(pow(Et-Eph1,2) + pow(ip.G1+ip.sigma,2));
      double x_phi_1680_IM= p.A2*(ip.G2+ip.sigma)/(pow(Et-Eph2,2) + pow(ip.G2+ip.sigma,2));
      double x_phi_2200_IM= p.A3*(ip.G3+ip.sigma)/(pow(Et-Eph3,2) + pow(ip.G3+ip.sigma,2));

      double x_res=0;

      complex I(0.0,1.0);
      double e=ip.sigma;

      for(int n=0;n<(signed)As.size();n++) {

	complex r= ((Ephs[n]-I*e -Et + ((TH > Et)?2.0*I*Gs[n]:0.0))*M_PI + 2.0*(-Ephs[n] + I*e + Et)*atan( (Ephs[n] -TH)/Gs[n]) + Gs[n]*( 2.0*I*atan( e/(Et-TH)) -log( pow(Ephs[n]-TH,2) + pow(Gs[n],2)) + log( pow(e,2) + pow(Et-TH,2))))/( (pow(Gs[n],2) + (Et -Ephs[n] +I*e)*(Et-Ephs[n]+I*e)));

	x_res += imag(r)*As[n]/(2.0*M_PI);

	//if(ip.sigma==0.0001) {
	//if(n==0) cout<<"IM phi : original: "<<x_phi_IM<<" , new: "<<imag(r)*As[n]/(2.0*M_PI)<<endl;
	//if(n==1) cout<<"IM phi(1680) : original: "<<x_phi_1680_IM<<" , new: "<<imag(r)*As[n]/(2.0*M_PI)<<endl;
	//if(n==2) cout<<"IM phi(2200) : original: "<<x_phi_2200_IM<<" , new: "<<imag(r)*As[n]/(2.0*M_PI)<<endl;
	//}

      }

      //if(ip.sigma==0.0001) cout<<"x_res(IM) original: "<<x_res_ori<<" , new: "<<x_res<<endl;

      complex z(Et,ip.sigma);
      double x_cont= imag( -D*log(1.0 - z/TH)/z)/M_PI;

      if(fix_TH) return x_res;
      
      if(MODE_SPECTR=="sqrt") x_cont= imag( 2*D*atanh(sqrt(z/TH))/sqrt(z))/M_PI;
      if(MODE_SPECTR=="quad") x_cont = imag( -D*log(1.0-z/TH)/(z*z)  -D/(TH*z))/M_PI;

      return  x_res + x_cont;
      
          	 	  
    };

  
    bf_07_comb.measurement=  [](const fpar_07_comb &p, const ipar_07_comb &ip) {

      return ip.FF;
      
    };
    
    bf_07_comb.error=  [](const fpar_07_comb &p, const ipar_07_comb &ip) {

      return ip.FF_err;
      
    };
    
  bf_07_comb_ch2.ansatz= bf_07_comb.ansatz;
  bf_07_comb_ch2.measurement = bf_07_comb.measurement;
  bf_07_comb_ch2.error = bf_07_comb.error;

 
  //start fitting
  //fill the data
  vector<vector<ipar_07_comb>> data_07_comb(Njacks);
  vector<vector<ipar_07_comb>> data_07_comb_ch2(1);
  //allocate space for output result
  boot_fit_data<fpar_07_comb> Bt_fit_07_comb;
  boot_fit_data<fpar_07_comb> Bt_fit_07_comb_ch2;

  boot_fit_data<fpar_07_comb> Bt_fit_07_comb_sqrt;
  boot_fit_data<fpar_07_comb> Bt_fit_07_comb_sqrt_ch2;

  
  
  for(auto &data_iboot: data_07_comb) data_iboot.resize(2*sigma_to_fit);
  for(auto &data_iboot: data_07_comb_ch2) data_iboot.resize(2*sigma_to_fit);
  for(int ijack=0;ijack<Njacks;ijack++) {
    for(int is=0;is<sigma_to_fit;is++) {
      data_07_comb[ijack][is].FF = (FT_s_RE_sigma[is]).distr[ijack];
      data_07_comb[ijack][is].FF_err= (FT_s_RE_sigma[is]).err();
      data_07_comb[ijack][is].sigma = sigma_simulated[ixg][is];
      data_07_comb[ijack][is].I_ori = FT_I_sub.distr[ijack];
      data_07_comb[ijack][is].G1= G1_distr.distr[ijack];
      data_07_comb[ijack][is].G2= G2_distr.distr[ijack];
      data_07_comb[ijack][is].G3= G3_distr.distr[ijack];
      data_07_comb[ijack][is].Is_RE = true;

      data_07_comb[ijack][is+sigma_to_fit].FF = (FT_s_IM_sigma[is]).distr[ijack];
      data_07_comb[ijack][is+sigma_to_fit].FF_err= (FT_s_IM_sigma[is]).err();
      data_07_comb[ijack][is+sigma_to_fit].sigma = sigma_simulated[ixg][is];
      data_07_comb[ijack][is+sigma_to_fit].I_ori = FT_I_sub.distr[ijack];
      data_07_comb[ijack][is+sigma_to_fit].G1= G1_distr.distr[ijack];
      data_07_comb[ijack][is+sigma_to_fit].G2= G2_distr.distr[ijack];
      data_07_comb[ijack][is+sigma_to_fit].G3= G3_distr.distr[ijack];
      data_07_comb[ijack][is+sigma_to_fit].Is_RE = false;
     
        
      if(ijack==0) {

	data_07_comb_ch2[ijack][is].FF = (FT_s_RE_sigma[is]).ave();
	data_07_comb_ch2[ijack][is].FF_err= (FT_s_RE_sigma[is]).err();
	data_07_comb_ch2[ijack][is].sigma = sigma_simulated[ixg][is];
	data_07_comb_ch2[ijack][is].I_ori = FT_I_sub.ave();
	data_07_comb_ch2[ijack][is].G1= G1_distr.ave();
	data_07_comb_ch2[ijack][is].G2= G2_distr.ave();
	data_07_comb_ch2[ijack][is].G3= G3_distr.ave();
	data_07_comb_ch2[ijack][is].Is_RE = true;
	
	data_07_comb_ch2[ijack][is+sigma_to_fit].FF = (FT_s_IM_sigma[is]).ave();
	data_07_comb_ch2[ijack][is+sigma_to_fit].FF_err= (FT_s_IM_sigma[is]).err();
	data_07_comb_ch2[ijack][is+sigma_to_fit].sigma = sigma_simulated[ixg][is];
	data_07_comb_ch2[ijack][is+sigma_to_fit].I_ori = FT_I_sub.ave();
	data_07_comb_ch2[ijack][is+sigma_to_fit].G1= G1_distr.ave();
	data_07_comb_ch2[ijack][is+sigma_to_fit].G2= G2_distr.ave();
	data_07_comb_ch2[ijack][is+sigma_to_fit].G3= G3_distr.ave();
	data_07_comb_ch2[ijack][is+sigma_to_fit].Is_RE = false;

	
      }
    }
  }


  for(int ijack=0;ijack<Njacks;ijack++) {

    double gc=0.047;
    if(MESON=="B1s") gc=0.042;
    if(MESON=="B3s") gc=0.040;

    bf_07_comb.Append_to_prior("TH", (fix_TH)?1e3:0.5, 0.5); bf_07_comb_ch2.Append_to_prior("TH", (fix_TH)?1e3:0.5 ,0.5);
    bf_07_comb.Append_to_prior("A", (Fit_I_TO)?-0.005:0.0, 0.005); bf_07_comb_ch2.Append_to_prior("A", (Fit_I_TO)?-0.005:0.0 ,0.005);
    bf_07_comb.Append_to_prior("A1", gc, 0.10*gc); bf_07_comb_ch2.Append_to_prior("A1", gc,0.10*gc);
    bf_07_comb.Append_to_prior("A2", 0.5*gc, 0.5*gc); bf_07_comb_ch2.Append_to_prior("A2", 0.5*gc,0.5*gc);
    bf_07_comb.Append_to_prior("A3", 0.5*gc, 0.5*gc); bf_07_comb_ch2.Append_to_prior("A3", 0.5*gc,0.5*gc);
    
  }
  
  //append
  bf_07_comb.Append_to_input_par(data_07_comb);
  bf_07_comb_ch2.Append_to_input_par(data_07_comb_ch2);
  //fit
  cout<<"Fitting real part of the form factor for ixg: "<<ixg<<endl;

  cout<<"I TO: "<<FT_I_ori.ave()<<" +- "<<FT_I_ori.err()<<endl;
  Bt_fit_07_comb= bf_07_comb.Perform_bootstrap_fit();
  Bt_fit_07_comb_ch2= bf_07_comb_ch2.Perform_bootstrap_fit();
  MODE_SPECTR="sqrt";
  Bt_fit_07_comb_sqrt= bf_07_comb.Perform_bootstrap_fit();
  Bt_fit_07_comb_sqrt_ch2= bf_07_comb_ch2.Perform_bootstrap_fit();
  Use_five_finest=true;

  distr_t A(UseJack),  A1(UseJack), A2(UseJack), A3(UseJack), TH(UseJack);

  distr_t A_sqrt(UseJack),  A1_sqrt(UseJack), A2_sqrt(UseJack), A3_sqrt(UseJack), TH_sqrt(UseJack);

 

  for(int ijack=0;ijack<Njacks;ijack++) {

    A.distr.push_back( Bt_fit_07_comb.par[ijack].A);
    A1.distr.push_back( Bt_fit_07_comb.par[ijack].A1);
    A2.distr.push_back( Bt_fit_07_comb.par[ijack].A2);
    A3.distr.push_back( Bt_fit_07_comb.par[ijack].A3);
    TH.distr.push_back( Bt_fit_07_comb.par[ijack].TH);


    A_sqrt.distr.push_back( Bt_fit_07_comb_sqrt.par[ijack].A);
    A1_sqrt.distr.push_back( Bt_fit_07_comb_sqrt.par[ijack].A1);
    A2_sqrt.distr.push_back( Bt_fit_07_comb_sqrt.par[ijack].A2);
    A3_sqrt.distr.push_back( Bt_fit_07_comb_sqrt.par[ijack].A3);
    TH_sqrt.distr.push_back( Bt_fit_07_comb_sqrt.par[ijack].TH);
  }

  //get fit function


  distr_t_list F_print_COMB_RE(UseJack, Nsigmas_to_print, Njacks), F_print_COMB_IM(UseJack, Nsigmas_to_print, Njacks);

  distr_t_list F_print_COMB_SQRT_RE(UseJack, Nsigmas_to_print, Njacks), F_print_COMB_SQRT_IM(UseJack, Nsigmas_to_print, Njacks);


  distr_t F_RE_lin(UseJack); distr_t F_RE_sqrt(UseJack);
  distr_t F_IM_lin(UseJack); distr_t F_IM_sqrt(UseJack);

  for(int ijack=0;ijack<Njacks;ijack++) {
    ipar_07_comb fS;
    fS.sigma= 1e-10;
    fS.I_ori=FT_I_sub.distr[ijack];
    fS.G1 = G1_distr.distr[ijack];
    fS.G2 = G2_distr.distr[ijack];
    fS.G3 = G3_distr.distr[ijack];

    MODE_SPECTR="";
    fS.Is_RE=true;
    F_RE_lin.distr.push_back( bf_07_comb.ansatz( Bt_fit_07_comb.par[ijack], fS) );
    fS.Is_RE=false;
    F_IM_lin.distr.push_back( bf_07_comb.ansatz( Bt_fit_07_comb.par[ijack], fS) );
    
    MODE_SPECTR="sqrt";
    fS.Is_RE=true;
    F_RE_sqrt.distr.push_back( bf_07_comb.ansatz( Bt_fit_07_comb_sqrt.par[ijack], fS));
    fS.Is_RE=false;
    F_IM_sqrt.distr.push_back( bf_07_comb.ansatz( Bt_fit_07_comb_sqrt.par[ijack], fS));

  }
  
  for(int n=0;n<Nsigmas_to_print;n++) {

    ipar_07_comb fS;
    fS.sigma= sigmas_to_print[n];
   
    for(int ijack=0;ijack<Njacks;ijack++) {

      fS.I_ori= FT_I_sub.distr[ijack];
      fS.G1 = G1_distr.distr[ijack];
      fS.G2 = G2_distr.distr[ijack];
      fS.G3 = G3_distr.distr[ijack];

      MODE_SPECTR="";
      fS.Is_RE=true;
      F_print_COMB_RE.distr_list[n].distr[ijack] = bf_07_comb.ansatz( Bt_fit_07_comb.par[ijack], fS);
      fS.Is_RE=false;
      F_print_COMB_IM.distr_list[n].distr[ijack] = bf_07_comb.ansatz( Bt_fit_07_comb.par[ijack], fS);
      
      MODE_SPECTR="sqrt";
      fS.Is_RE=true;
      F_print_COMB_SQRT_RE.distr_list[n].distr[ijack] = bf_07_comb.ansatz( Bt_fit_07_comb_sqrt.par[ijack], fS);
      fS.Is_RE=false;
      F_print_COMB_SQRT_IM.distr_list[n].distr[ijack] = bf_07_comb.ansatz( Bt_fit_07_comb_sqrt.par[ijack], fS);
   
    }
    
  }

  //print fit function
  


  //###############################################################################################################
  //###############################################################################################################
  //###############################################################################################################


  //######################################################
  //######################################################

  






  //print spectral density at epsilon=0
  Vfloat Ergs;
  int Nergs=5001;
  for(int ierg=0;ierg<Nergs;ierg++) Ergs.push_back(ierg*0.001);

  distr_t_list sp_dens(UseJack, Nergs, Njacks);
  distr_t_list sp_dens_s_small(UseJack, Nergs,Njacks);

  MODE_SPECTR="";

  for(int ierg=0;ierg<Nergs;ierg++) {

    double E= Ergs[ierg];

    for(int ijack=0;ijack<Njacks;ijack++) {

     
      

      double G1= G1_distr.distr[ijack];
      double G2= G2_distr.distr[ijack];
      double G3= G3_distr.distr[ijack];
      double A1_v= A1.distr[ijack];
      double A2_v= A2.distr[ijack];
      double A3_v= A3.distr[ijack];
      double TH_v= TH.distr[ijack];
      
      Et = E;
      ipar_07_comb Ip_s;
      Ip_s.Is_RE = false;
      Ip_s.G1= G1;
      Ip_s.G2= G2;
      Ip_s.G3= G3;
      Ip_s.sigma=0.0001;

      double TH= Eph3 + sqrt(TH_v*TH_v);

      fpar_07_comb Fp_s({TH_v, 0.0, A1_v, A2_v, A3_v});
      
      double res_s= bf_07_comb.ansatz(Fp_s, Ip_s);

     
      double res= 2*A1_v*G1/( pow(E-Eph1,2) + pow(G1,2)) + 2*A2_v*G2/( pow(E-Eph2,2) + pow(G2,2)) + 2*A3_v*G3/( pow(E-Eph3,2) + pow(G3,2));
      

      double D= 2*A1_v*G1/( pow(TH-Eph1,2) + pow(G1,2)) + 2*A2_v*G2/( pow(TH-Eph2,2) + pow(G2,2)) + 2*A3_v*G3/( pow(TH-Eph3,2) + pow(G3,2));
      double rho_cont=0;
      
      if(MODE_SPECTR=="sqrt") { D *= sqrt(TH); rho_cont = D/sqrt(E); }
      if(MODE_SPECTR=="quad") { D *= pow(TH,2);  rho_cont = D/pow(E,2);}
      else { D *= TH; rho_cont = D/E;}
      
      
      if(E > TH) sp_dens.distr_list[ierg].distr[ijack] = rho_cont;
      else  sp_dens.distr_list[ierg].distr[ijack] =  res;

      
      sp_dens_s_small.distr_list[ierg].distr[ijack] = res_s; 
      
    }

  }

  Print_To_File({}, {Ergs, sp_dens.ave(), sp_dens.err(), sp_dens_s_small.ave(), sp_dens_s_small.err()}, path_out+"/FF_d_extr/sp_dens_ixg_"+to_string(ixg)+".data", "", "");


  //##############################################################################################################################################




  
    


  //print data and fit parameters
  boost::filesystem::create_directory(path_out+"/FF_d_extr");
  string out_path_data = path_out+"/FF_d_extr/FT_ixg_"+to_string(ixg)+".data" ;
  string out_path_fit_func = path_out+"/FF_d_extr/FT_ixg_"+to_string(ixg)+".fit_func";
  string out_path_fit_func_comb = path_out+"/FF_d_extr/FT_ixg_"+to_string(ixg)+"_comb.fit_func";
  string out_path_fit_func_comb_sqrt = path_out+"/FF_d_extr/FT_ixg_"+to_string(ixg)+"_comb_sqrt.fit_func";
  
  Print_To_File({}, {sigma_simulated[ixg], FT_s_RE_sigma.ave(), FT_s_RE_sigma.err(), FT_s_IM_sigma.ave(), FT_s_IM_sigma.err(), FT_s_IM_Gauss_sigma.ave(), FT_s_IM_Gauss_sigma.err()}, out_path_data, "", "");
  Print_To_File({}, {sigma_simulated[ixg], FT_s_RE_sigma_15.ave(), FT_s_RE_sigma_15.err(), FT_s_IM_sigma_15.ave(), FT_s_IM_sigma_15.err()},  path_out+"/FF_d_extr/FT_15_ixg_"+to_string(ixg)+".data", "", "");
  Print_To_File({}, {sigma_simulated[ixg], FT_s_RE_sigma_7.ave(), FT_s_RE_sigma_7.err(), FT_s_IM_sigma_7.ave(), FT_s_IM_sigma_7.err()},  path_out+"/FF_d_extr/FT_7_ixg_"+to_string(ixg)+".data", "", "");
  Print_To_File({}, {sigma_simulated[ixg], FT_s_RE_sigma_no_sub.ave(), FT_s_RE_sigma_no_sub.err()},  path_out+"/FF_d_extr/FT_no_sub_ixg_"+to_string(ixg)+".data", "", "");
  Print_To_File({}, {sigmas_to_print, F_print_RE.ave(), F_print_RE.err(), F_print_IM.ave(), F_print_IM.err(), F_print_IM_Gauss.ave(), F_print_IM_Gauss.err()}, out_path_fit_func, "", "");
  Print_To_File({}, {sigmas_to_print_0, F_print_RE_0.ave(), F_print_RE_0.err(), F_print_IM_0.ave(), F_print_IM_0.err()}, out_path_fit_func+"_0", "", "");

  Print_To_File({}, {sigmas_to_print, F_print_COMB_RE.ave(), F_print_COMB_RE.err(), F_print_COMB_IM.ave(), F_print_COMB_IM.err()}, out_path_fit_func_comb, "", "");

  Print_To_File({}, {sigmas_to_print, F_print_COMB_SQRT_RE.ave(), F_print_COMB_SQRT_RE.err(), F_print_COMB_SQRT_IM.ave(), F_print_COMB_SQRT_IM.err()}, out_path_fit_func_comb_sqrt, "", "");
  
  
  double syst_RE= fabs( R_RE.ave() - R_RE_red.ave())*erf(  fabs(R_RE.ave() - R_RE_red.ave())/(sqrt(2.0)*R_RE_red.err())   );
  double syst_IM= fabs( R_IM.ave() - R_IM_red.ave())*erf(  fabs(R_IM.ave() - R_IM_red.ave())/(sqrt(2.0)*R_IM_red.err())   );

  R_RE = R_RE.ave() + (R_RE -R_RE.ave())*sqrt( pow(R_RE.err(),2) + pow(syst_RE,2))/R_RE.err();
  R_IM = R_IM.ave() + (R_IM -R_IM.ave())*sqrt( pow(R_IM.err(),2) + pow(syst_IM,2))/R_IM.err();

  distr_t ave_RE = (1.0/2.0)*( R_RE +  F_RE_sqrt);
  distr_t ave_IM = (1.0/2.0)*( R_IM +  F_IM_sqrt);

  double syst_RE_tot= 0.5*fabs( R_RE.ave() - F_RE_sqrt.ave());
  double syst_IM_tot= 0.5*fabs( R_IM.ave() - F_IM_sqrt.ave());

  ave_RE = ave_RE.ave() + (ave_RE - ave_RE.ave())*sqrt( pow(ave_RE.err(),2) + pow(syst_RE_tot,2))/ave_RE.err();

  ave_IM = ave_IM.ave() + (ave_IM - ave_IM.ave())*fabs( ave_IM.err() + syst_IM_tot)/ave_IM.err();
  
  return_class.F_T_d_RE.distr_list.push_back( ave_RE );
  return_class.F_T_d_IM.distr_list.push_back( ave_IM );
  
  
  
  
    
    
   return_class.F_T_u.distr_list.push_back( FT_b);
    
    
    
    
  }

  //print extrapolated ori
  Print_To_File({}, {xg_list.ave() ,F_T_d_I_ori_extr.ave(), F_T_d_I_ori_extr.err()}, path_out+"/FF_d_extr/F_T_d_ori.data", "", "");
  
  return_class.num_xg= n_xg;
  return_class.MESON= MESON;
  

    
 




  return return_class;

}
