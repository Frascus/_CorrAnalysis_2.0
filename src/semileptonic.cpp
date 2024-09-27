#include "../include/semileptonic.h"




using namespace std;



const double M2PiPhys=pow(0.135,2);
const double alpha = 1/137.04;
const double e2 = alpha*4.0*M_PI;
const int Nboots= 800;
const bool UseJack=1;
const int Njacks=((UseJack==true)?25:Nboots);
const string channel="D_Klnu";
//Tseps 64, 44, 36, 30


void semileptonic_FF_analysis() {

  Vint Tseps({28,30,36,40});
  for(auto &Tsep: Tseps) Get_semileptonic_FF(Tsep);
}


void Get_semileptonic_FF(int Tsep) {


  //lambda functions to be used later on
  auto LOG= [](double m) { return log(m);};
  auto SQRT_DL= [](double m, double t) { return sqrt(m);};
  auto EXP = [](double m) { return exp(m);};

  //create directories
  boost::filesystem::create_directory("../data/semileptonic");
  boost::filesystem::create_directory("../data/semileptonic/"+channel);
  boost::filesystem::create_directory("../data/semileptonic/"+channel+"/masses");
  boost::filesystem::create_directory("../data/semileptonic/"+channel+"/FF");


  data_t data_2pts_SM_P, data_2pts_SSM_P;
  data_t data_2pts_SM_K, data_2pts_SSM_K;
  data_t data_3pts_V0_F, data_3pts_V0_B;
  data_t data_3pts_V0_P, data_3pts_V0_K;
  data_t data_3pts_S0_F, data_3pts_S0_B;


  //Read data

  string tg_3pt= to_string(Tsep);

  //2pts function
    

  data_2pts_SM_P.Read("../D_Klnu_data", "mes_contr_MES_CONTR_2_PT_D", "P5P5");
  data_2pts_SSM_P.Read("../D_Klnu_data", "mes_contr_MES_CONTR_2_PT_D_SM", "P5P5");
  data_2pts_SM_K.Read("../D_Klnu_data", "mes_contr_MES_CONTR_2_PT_K", "P5P5");
  data_2pts_SSM_K.Read("../D_Klnu_data", "mes_contr_MES_CONTR_2_PT_K_SM", "P5P5");
  
  //3pts function
  data_3pts_V0_F.Read("../D_Klnu_data", "mes_contr_MES_CONTR_3_PT_F_"+tg_3pt, "V0P5");
  data_3pts_V0_B.Read("../D_Klnu_data", "mes_contr_MES_CONTR_3_PT_B_"+tg_3pt, "V0P5");

  data_3pts_S0_F.Read("../D_Klnu_data", "mes_contr_MES_CONTR_3_PT_F_"+tg_3pt, "S0P5");
  data_3pts_S0_B.Read("../D_Klnu_data", "mes_contr_MES_CONTR_3_PT_B_"+tg_3pt, "S0P5");

  data_3pts_V0_P.Read("../D_Klnu_data", "mes_contr_MES_CONTR_3_PT_DIAG_c_"+tg_3pt, "V0P5");
  data_3pts_V0_K.Read("../D_Klnu_data", "mes_contr_MES_CONTR_3_PT_DIAG_s_"+tg_3pt, "V0P5");
  

  


  int Nens = data_2pts_SM_P.size;
  GaussianMersenne GM(94538123);

  //resample RCs
  distr_t ZA_A(UseJack), ZA_B(UseJack), ZA_C(UseJack), ZA_D(UseJack);
  distr_t ZV_A(UseJack), ZV_B(UseJack), ZV_C(UseJack), ZV_D(UseJack);
  distr_t a_A(UseJack), a_B(UseJack), a_C(UseJack), a_D(UseJack);

  double fmTGeV= 1.0/0.197327;

  LatticeInfo L_info_A, L_info_B, L_info_C, L_info_D;
  L_info_A.LatInfo_new_ens("cA211a.12.48");
  L_info_B.LatInfo_new_ens("cB211b.072.96");
  L_info_C.LatInfo_new_ens("cC211a.06.80");
  L_info_D.LatInfo_new_ens("cD211a.054.96");
  

  for(int ijack=0; ijack<Njacks;ijack++) {

    ZA_A.distr.push_back( L_info_A.Za_WI_strange + GM()*L_info_A.Za_WI_strange_err/sqrt( Njacks -1.0));
    ZV_A.distr.push_back( L_info_A.Zv_WI_strange + GM()*L_info_A.Zv_WI_strange_err/sqrt( Njacks -1.0));

    ZA_B.distr.push_back( L_info_B.Za_WI_strange + GM()*L_info_B.Za_WI_strange_err/sqrt( Njacks -1.0));
    ZV_B.distr.push_back( L_info_B.Zv_WI_strange + GM()*L_info_B.Zv_WI_strange_err/sqrt( Njacks -1.0));

    ZA_C.distr.push_back( L_info_C.Za_WI_strange + GM()*L_info_C.Za_WI_strange_err/sqrt( Njacks -1.0));
    ZV_C.distr.push_back( L_info_C.Zv_WI_strange + GM()*L_info_C.Zv_WI_strange_err/sqrt( Njacks -1.0));

    ZA_D.distr.push_back( L_info_D.Za_WI_strange + GM()*L_info_D.Za_WI_strange_err/sqrt( Njacks -1.0));
    ZV_D.distr.push_back( L_info_D.Zv_WI_strange + GM()*L_info_D.Zv_WI_strange_err/sqrt( Njacks -1.0));

    a_A.distr.push_back( L_info_A.a_from_afp*fmTGeV + GM()*L_info_A.a_from_afp_err*fmTGeV/sqrt(Njacks-1.0));
    a_B.distr.push_back( L_info_B.a_from_afp*fmTGeV + GM()*L_info_B.a_from_afp_err*fmTGeV/sqrt(Njacks-1.0));
    a_C.distr.push_back( L_info_C.a_from_afp*fmTGeV + GM()*L_info_C.a_from_afp_err*fmTGeV/sqrt(Njacks-1.0));
    a_D.distr.push_back( L_info_D.a_from_afp*fmTGeV + GM()*L_info_D.a_from_afp_err*fmTGeV/sqrt(Njacks-1.0));

  }



  for(int iens=0;iens<Nens;iens++) { //loop over ensemble

    cout<<"Analyzing ensemble: "<<data_2pts_SM_P.Tag[iens]<<endl;

    boost::filesystem::create_directory("../data/semileptonic/"+channel+"/masses/"+data_2pts_SM_P.Tag[iens]);
    boost::filesystem::create_directory("../data/semileptonic/"+channel+"/FF/"+data_2pts_SM_P.Tag[iens]);


    //setup correlator analysis
    //Lattice info
    LatticeInfo L_info;
    L_info.LatInfo_new_ens(data_2pts_SM_P.Tag[iens]);
    CorrAnalysis Corr(UseJack, Njacks,Nboots, iens);
    Corr.Nt = data_2pts_SM_P.nrows[iens];
    Corr.Reflection_sign=1;
    Corr.Perform_Nt_t_average=1;


    double mu_ut=0;
    double mu_dt=0;
    //RCs
    distr_t Za, Zv, a_distr;
    if(data_2pts_SM_P.Tag[iens].substr(1,1)=="A") { Za= ZA_A; Zv=ZV_A;a_distr=a_A;}
    else if(data_2pts_SM_P.Tag[iens].substr(1,1)=="B") { Za= ZA_B; Zv=ZV_B;a_distr=a_B; mu_ut= 0.237; mu_dt=0.0184  ;   }
    else if(data_2pts_SM_P.Tag[iens].substr(1,1)=="C") { Za= ZA_C; Zv=ZV_C; a_distr=a_C;}
    else if(data_2pts_SM_P.Tag[iens].substr(1,1)=="D") { Za= ZA_D; Zv=ZV_D;a_distr=a_D;}
    else crash("Ensemble: "+data_2pts_SM_P.Tag[iens]+" not recognised");



    //set time interval for 2pt obs
    int Tmin_P, Tmax_P, Tmin_K, Tmax_K;
    if(data_2pts_SM_P.Tag[iens].substr(1,1) =="A") {Tmin_P=19; Tmax_P=35; Tmin_K = 19; Tmax_K=35;}
    else if(data_2pts_SM_P.Tag[iens] =="cB211b.072.64") {Tmin_P=14; Tmax_P=20; Tmin_K = 13; Tmax_K=30;}
    else if(data_2pts_SM_P.Tag[iens] =="cB211b.072.96") {Tmin_P=19; Tmax_P=35; Tmin_K = 19; Tmax_K=35;}
    else if(data_2pts_SM_P.Tag[iens].substr(1,1) == "C") {Tmin_P=19; Tmax_P=35; Tmin_K = 19; Tmax_K=35;}
    else if(data_2pts_SM_P.Tag[iens].substr(1,1) == "D") {Tmin_P=19; Tmax_P=35; Tmin_K = 19; Tmax_K=35;}
    else crash("In fixing [Tmin, Tmax] for P and K, Ensemble: "+data_2pts_SM_P.Tag[iens]+" not recognized");


    distr_t_list pt2_distr_SM_P = Corr.corr_t(data_2pts_SM_P.col(0)[iens], "../data/semileptonic/"+channel+"/masses/"+data_2pts_SM_P.Tag[iens]+"/corr_2pt_P_SM.dat");
    distr_t_list pt2_distr_SSM_P = Corr.corr_t(data_2pts_SSM_P.col(0)[iens], "../data/semileptonic/"+channel+"/masses/"+data_2pts_SM_P.Tag[iens]+"/corr_2pt_P_SSM.dat");
    distr_t_list pt2_distr_SSM_P_tilde= 0.5*( pt2_distr_SSM_P + distr_t_list::f_of_distr_list(SQRT_DL,  pt2_distr_SSM_P*pt2_distr_SSM_P - pt2_distr_SSM_P.distr_list[Corr.Nt/2]*pt2_distr_SSM_P.distr_list[Corr.Nt/2]));
    distr_t_list eff_mass_SM_P = Corr.effective_mass_t(pt2_distr_SM_P, "../data/semileptonic/"+channel+"/masses/"+data_2pts_SM_P.Tag[iens]+"/eff_mass_P_SM.dat");
    distr_t_list eff_mass_SSM_P = Corr.effective_mass_t(pt2_distr_SSM_P, "../data/semileptonic/"+channel+"/masses/"+data_2pts_SM_P.Tag[iens]+"/eff_mass_P_SSM.dat");
   
    distr_t_list pt2_distr_SM_K = Corr.corr_t(data_2pts_SM_K.col(0)[iens], "../data/semileptonic/"+channel+"/masses/"+data_2pts_SM_P.Tag[iens]+"/corr_2pt_K_SM.dat");
    distr_t_list pt2_distr_SSM_K = Corr.corr_t(data_2pts_SSM_K.col(0)[iens], "../data/semileptonic/"+channel+"/masses/"+data_2pts_SM_P.Tag[iens]+"/corr_2pt_K_SSM.dat");
    distr_t_list pt2_distr_SSM_K_tilde= 0.5*( pt2_distr_SSM_K + distr_t_list::f_of_distr_list(SQRT_DL,  pt2_distr_SSM_K*pt2_distr_SSM_K - pt2_distr_SSM_K.distr_list[Corr.Nt/2]*pt2_distr_SSM_K.distr_list[Corr.Nt/2]));
    distr_t_list eff_mass_SM_K = Corr.effective_mass_t(pt2_distr_SM_K, "../data/semileptonic/"+channel+"/masses/"+data_2pts_SM_P.Tag[iens]+"/eff_mass_K_SM.dat");
    distr_t_list eff_mass_SSM_K = Corr.effective_mass_t(pt2_distr_SSM_K, "../data/semileptonic/"+channel+"/masses/"+data_2pts_SM_P.Tag[iens]+"/eff_mass_K_SSM.dat");


    
    if(Tsep != Corr.Nt/2) Corr.Perform_Nt_t_average=0;
  
    distr_t_list pt3_distr_V0_F_RE = Corr.corr_t(data_3pts_V0_F.col(0)[iens], "../data/semileptonic/"+channel+"/FF/"+data_2pts_SM_P.Tag[iens]+"/corr_3pt_V0_F_RE_Tsep_"+to_string(Tsep)+".dat");
    distr_t_list pt3_distr_V0_B_RE = Corr.corr_t(data_3pts_V0_B.col(0)[iens], "../data/semileptonic/"+channel+"/FF/"+data_2pts_SM_P.Tag[iens]+"/corr_3pt_V0_B_RE_Tsep_"+to_string(Tsep)+".dat");
    distr_t_list pt3_distr_S0_F_RE = Corr.corr_t(data_3pts_S0_F.col(0)[iens], "../data/semileptonic/"+channel+"/FF/"+data_2pts_SM_P.Tag[iens]+"/corr_3pt_S0_F_RE_Tsep_"+to_string(Tsep)+".dat");
    distr_t_list pt3_distr_S0_B_RE = Corr.corr_t(data_3pts_S0_B.col(0)[iens], "../data/semileptonic/"+channel+"/FF/"+data_2pts_SM_P.Tag[iens]+"/corr_3pt_S0_B_RE_Tsep_"+to_string(Tsep)+".dat");
    
    distr_t_list pt3_distr_V0_F_IM = Corr.corr_t(data_3pts_V0_F.col(1)[iens], "../data/semileptonic/"+channel+"/FF/"+data_2pts_SM_P.Tag[iens]+"/corr_3pt_V0_F_IM_Tsep_"+to_string(Tsep)+".dat");
    distr_t_list pt3_distr_V0_B_IM = Corr.corr_t(data_3pts_V0_B.col(1)[iens], "../data/semileptonic/"+channel+"/FF/"+data_2pts_SM_P.Tag[iens]+"/corr_3pt_V0_B_IM_Tsep_"+to_string(Tsep)+".dat");
    distr_t_list pt3_distr_S0_F_IM = Corr.corr_t(data_3pts_S0_F.col(1)[iens], "../data/semileptonic/"+channel+"/FF/"+data_2pts_SM_P.Tag[iens]+"/corr_3pt_S0_F_IM_Tsep_"+to_string(Tsep)+".dat");
    distr_t_list pt3_distr_S0_B_IM = Corr.corr_t(data_3pts_S0_B.col(1)[iens], "../data/semileptonic/"+channel+"/FF/"+data_2pts_SM_P.Tag[iens]+"/corr_3pt_S0_B_IM_Tsep_"+to_string(Tsep)+".dat");

    distr_t_list pt3_distr_V0_P_RE = Corr.corr_t(data_3pts_V0_P.col(0)[iens], "../data/semileptonic/"+channel+"/FF/"+data_2pts_SM_P.Tag[iens]+"/corr_3pt_V0_P_RE_Tsep_"+to_string(Tsep)+".dat");
    distr_t_list pt3_distr_V0_K_RE = Corr.corr_t(data_3pts_V0_K.col(0)[iens], "../data/semileptonic/"+channel+"/FF/"+data_2pts_SM_K.Tag[iens]+"/corr_3pt_V0_K_RE_Tsep_"+to_string(Tsep)+".dat");


    distr_t_list pt3_distr_V0_P_IM = Corr.corr_t(data_3pts_V0_P.col(1)[iens], "../data/semileptonic/"+channel+"/FF/"+data_2pts_SM_P.Tag[iens]+"/corr_3pt_V0_P_IM_Tsep_"+to_string(Tsep)+".dat");
    distr_t_list pt3_distr_V0_K_IM = Corr.corr_t(data_3pts_V0_K.col(1)[iens], "../data/semileptonic/"+channel+"/FF/"+data_2pts_SM_K.Tag[iens]+"/corr_3pt_V0_K_IM_Tsep_"+to_string(Tsep)+".dat");
    
    Corr.Perform_Nt_t_average=1;


    distr_t_list pt3_distr_V0_F= pt3_distr_V0_F_RE;
    distr_t_list pt3_distr_V0_B= pt3_distr_V0_B_RE;
    distr_t_list pt3_distr_S0_F= pt3_distr_S0_F_RE;
    distr_t_list pt3_distr_S0_B= pt3_distr_S0_B_RE;
    distr_t_list pt3_distr_V0_P= pt3_distr_V0_P_RE;
    distr_t_list pt3_distr_V0_K= pt3_distr_V0_K_RE;
    

    distr_t MP_SM, MP_SSM, ZP2_SM,  ZP2_SSM;
    Corr.Tmin= Tmin_P; Corr.Tmax= Tmax_P;
    MP_SM = Corr.Fit_distr(eff_mass_SM_P);
    MP_SSM= Corr.Fit_distr(eff_mass_SSM_P);

    distr_t_list residue_SSM_P = Corr.residue_t(pt2_distr_SSM_P, "../data/semileptonic/"+channel+"/masses/"+data_2pts_SM_P.Tag[iens]+"/residue_P_SSM.dat");
    distr_t_list residue_SM_P = Corr.residue_t(pt2_distr_SM_P, "../data/semileptonic/"+channel+"/masses/"+data_2pts_SM_P.Tag[iens]+"/residue_P_SM.dat");
    ZP2_SSM= Corr.Fit_distr(residue_SSM_P);
    ZP2_SM=  Corr.Fit_distr(residue_SM_P);
    

    distr_t MK_SM, MK_SSM, ZK2_SSM, ZK2_SM;
    Corr.Tmin= Tmin_K; Corr.Tmax= Tmax_K;
    MK_SM = Corr.Fit_distr(eff_mass_SM_K);
    MK_SSM = Corr.Fit_distr(eff_mass_SSM_K);
    distr_t_list residue_SSM_K = Corr.residue_t(pt2_distr_SSM_K, "../data/semileptonic/"+channel+"/masses/"+data_2pts_SM_P.Tag[iens]+"/residue_K_SSM.dat");
    distr_t_list residue_SM_K = Corr.residue_t(pt2_distr_SM_K, "../data/semileptonic/"+channel+"/masses/"+data_2pts_SM_P.Tag[iens]+"/residue_K_SM.dat");
    ZK2_SSM= Corr.Fit_distr(residue_SSM_K);
    ZK2_SM= Corr.Fit_distr(residue_SM_K);

    distr_t_list pt3_distr_V0_F_ave(UseJack);
    distr_t_list pt3_distr_V0_B_ave(UseJack);
    distr_t_list pt3_distr_S0_F_ave(UseJack);
    distr_t_list pt3_distr_S0_B_ave(UseJack);
    for(int t=0;t<Corr.Nt;t++) {
      pt3_distr_V0_F_ave.distr_list.push_back(  ( pt3_distr_V0_F.distr_list[t] + pt3_distr_V0_B.distr_list[ (Tsep-t + Corr.Nt)%Corr.Nt] )*0.5);
      pt3_distr_V0_B_ave.distr_list.push_back(  ( pt3_distr_V0_B.distr_list[t] + pt3_distr_V0_F.distr_list[ (Tsep-t + Corr.Nt)%Corr.Nt] )*0.5);
      pt3_distr_S0_F_ave.distr_list.push_back(  ( pt3_distr_S0_F.distr_list[t] + pt3_distr_S0_B.distr_list[ (Tsep-t + Corr.Nt)%Corr.Nt] )*0.5);
      pt3_distr_S0_B_ave.distr_list.push_back(  ( pt3_distr_S0_B.distr_list[t] + pt3_distr_S0_F.distr_list[ (Tsep-t + Corr.Nt)%Corr.Nt] )*0.5);
    }


    distr_t_list EST_V0 = 4.0*MP_SSM*MK_SSM*Zv*Zv*pt3_distr_V0_F*pt3_distr_V0_B/(pt2_distr_SSM_P_tilde.distr_list[Tsep]*pt2_distr_SSM_K_tilde.distr_list[Tsep]);
    distr_t_list EST_IMPR_V0 = 4.0*MP_SSM*MK_SSM*pt3_distr_V0_F*pt3_distr_V0_B/(pt3_distr_V0_P*pt3_distr_V0_K);
    distr_t_list EST_S0 = 4.0*MP_SSM*MK_SSM*pt3_distr_S0_F*pt3_distr_S0_B/(pt2_distr_SSM_P_tilde.distr_list[Tsep]*pt2_distr_SSM_K_tilde.distr_list[Tsep]);

    distr_t_list EST_V0_ave = 4.0*MP_SSM*MK_SSM*Zv*Zv*pt3_distr_V0_F_ave*pt3_distr_V0_B_ave/(pt2_distr_SSM_P_tilde.distr_list[Tsep]*pt2_distr_SSM_K_tilde.distr_list[Tsep]);
    distr_t_list EST_IMPR_V0_ave = 4.0*MP_SSM*MK_SSM*pt3_distr_V0_F_ave*pt3_distr_V0_B_ave/(pt3_distr_V0_P*pt3_distr_V0_K);
    distr_t_list EST_S0_ave = 4.0*MP_SSM*MK_SSM*pt3_distr_S0_F_ave*pt3_distr_S0_B_ave/(pt2_distr_SSM_P_tilde.distr_list[Tsep]*pt2_distr_SSM_K_tilde.distr_list[Tsep]);


    distr_t_list ZV_EST_P= 2.0*MP_SSM*pt2_distr_SSM_P.distr_list[Tsep]/pt3_distr_V0_P;
    distr_t_list ZV_EST_K= 2.0*MK_SSM*pt2_distr_SSM_K.distr_list[Tsep]/pt3_distr_V0_K;
    
    distr_t_list f0 = distr_t_list::f_of_distr_list(SQRT_DL, EST_V0)/(MP_SSM+MK_SSM);
    distr_t_list f0_IMPR = distr_t_list::f_of_distr_list(SQRT_DL, EST_IMPR_V0)/(MP_SSM+MK_SSM);
    distr_t_list f0_from_S0 = distr_t_list::f_of_distr_list(SQRT_DL, EST_S0)*(mu_ut - mu_dt)/(MP_SSM*MP_SSM - MK_SSM*MK_SSM);


    distr_t_list f0_ave = distr_t_list::f_of_distr_list(SQRT_DL, EST_V0_ave)/(MP_SSM+MK_SSM);
    distr_t_list f0_IMPR_ave = distr_t_list::f_of_distr_list(SQRT_DL, EST_IMPR_V0_ave)/(MP_SSM+MK_SSM);
    distr_t_list f0_from_S0_ave = distr_t_list::f_of_distr_list(SQRT_DL, EST_S0_ave)*(mu_ut - mu_dt)/(MP_SSM*MP_SSM - MK_SSM*MK_SSM);


    //


    


    

    Print_To_File({}, { EST_V0.ave(), EST_V0.err(), EST_IMPR_V0.ave(), EST_IMPR_V0.err()}, "../data/semileptonic/"+channel+"/FF/"+data_2pts_SM_P.Tag[iens]+"/V02_Tsep_"+to_string(Tsep)+".dat", "", "");
    Print_To_File({}, {EST_S0.ave(), EST_S0.err()}, "../data/semileptonic/"+channel+"/FF/"+data_2pts_SM_P.Tag[iens]+"/S02_Tsep_"+to_string(Tsep)+".dat", "", "");
    
    Print_To_File({}, { f0.ave(), f0.err(), f0_IMPR.ave(), f0_IMPR.err(), f0_from_S0.ave(), f0_from_S0.err()}, "../data/semileptonic/"+channel+"/FF/"+data_2pts_SM_P.Tag[iens]+"/f0_Tsep_"+to_string(Tsep)+".dat", "", "");

    Print_To_File({}, { f0_ave.ave(), f0_ave.err(), f0_IMPR_ave.ave(), f0_IMPR_ave.err(), f0_from_S0_ave.ave(), f0_from_S0_ave.err()}, "../data/semileptonic/"+channel+"/FF/"+data_2pts_SM_P.Tag[iens]+"/f0_averaged_Tsep_"+to_string(Tsep)+".dat", "", "");

    Print_To_File({ }, { (ZV_EST_P/Zv).ave(), (ZV_EST_P/Zv).err(), (ZV_EST_K/Zv).ave(), (ZV_EST_K/Zv).err()}, "../data/semileptonic/"+channel+"/FF/"+data_2pts_SM_P.Tag[iens]+"/ZV_est_Tsep_"+to_string(Tsep)+".dat", "", "");


    



  }
  



  return;
}
