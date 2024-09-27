#include "../include/Axion_l7.h"


using namespace std;
namespace plt = matplotlibcpp;

const double MPiPhys=0.135;
const bool Use_JB_distribution= false;
const bool UseJack=1;
const int Njacks=100;
const int Nboots=200;


void Neutral_pi_TM() {

  data_t pi_conn, pi_disc, pi_bubble, pi_plus;

   
  pi_conn.Read("../datasets_pi_0_TM", "conn_contr", "P5P5");
  pi_plus.Read("../datasets_pi_0_TM", "conn_pi_plus", "P5P5");
  pi_disc.Read("../datasets_pi_0_TM", "disc_contr", "P5");
  pi_bubble.Read("../datasets_pi_0_TM", "bubble", "P5");

  int Nens= pi_conn.size;
  cout<<"N_ens: "<<Nens<<endl;

 

  for(int i=0; i<Nens;i++) {
    distr_t_list pi_conn_distr(UseJack), pi_disc_distr(UseJack), pi_plus_distr(UseJack);
    CorrAnalysis Corr(UseJack, Njacks, Nboots);
    Corr.Perform_Nt_t_average=1;
    Corr.Reflection_sign=1;
    if(pi_conn.Tag[i].substr(0,1)=="A") Corr.Tmin=12;
    else if(pi_conn.Tag[i].substr(0,1)=="B") Corr.Tmin=16;
    else if(pi_conn.Tag[i].substr(0,1)=="D") Corr.Tmin=21;
    else if(pi_conn.Tag[i].substr(0,1)=="a") Corr.Tmin=12;
    else Corr.Tmin=21;
    Corr.Tmax = pi_conn.nrows[1]/2 - 6;
    Corr.Nt = pi_conn.nrows[i];
    //LatticeInfo L_info("LOCAL");

    //L_info.LatInfo(pi_conn.Tag[i].substr(0,1));
    double L, NNT, m;
    Read_pars_from_ensemble_tag(pi_conn.Tag[i], m, L, NNT);
    boost::filesystem::create_directory("../data");
    boost::filesystem::create_directory("../data/pi_0_TM");
    pi_conn_distr=Corr.corr_t(pi_conn.col(0)[i], "../data/pi_0_TM/pi_conn."+pi_conn.Tag[i]);
    pi_plus_distr=Corr.corr_t(pi_plus.col(0)[i], "../data/pi_0_TM/pi_plus."+pi_conn.Tag[i]);
    
    distr_t_list Pi_conn_eff_mass = Corr.effective_mass_t(pi_conn_distr, "../data/pi_0_TM/pi_conn_eff_mass."+pi_conn.Tag[i]);

    auto F= [&](const Vfloat& par) { if((signed)par.size() != 2) crash("Lambda function expects par[2], but par["+to_string((signed)par.size())+"] provided"); return par[0] -par[1]*par[1];};
    Vfloat bubble_average(pi_bubble.col(0)[i][0].size(),0);
    for(int t=0;t<Corr.Nt;t++) {
      if(bubble_average.size() != pi_bubble.col(0)[i][t].size()) crash("Number of confs not costants over time in Neutral_pi_TM");
      for(unsigned int iconf=0; iconf < bubble_average.size();iconf++) {
	bubble_average[iconf] += pi_bubble.col(0)[i][t][iconf]/Corr.Nt;
      }
    }
    for(int t=0;t<Corr.Nt;t++) {
	   if(UseJack==1) {
	     Jackknife J(10000,Njacks);
	     pi_disc_distr.distr_list.push_back(J.DoJack(F, 2, pi_disc.col(0)[i][t], bubble_average));
	   }
	   else {
	     Bootstrap B(Nboots, 432423, Corr.Nt);
	     pi_disc_distr.distr_list.push_back(B.DoBoot(F, 2, pi_disc.col(0)[i][t], bubble_average));
	   }
    }
    Print_To_File({}, {pi_disc_distr.ave(), pi_disc_distr.err()},"../data/pi_0_TM/pi_disc."+pi_disc.Tag[i]+".t"  ,"OUT", "");

  
    distr_t_list Pi_0_eff_mass= Corr.effective_mass_t(pi_conn_distr+2.0*pi_disc_distr, "../data/pi_0_TM/pi_0_eff_mass."+pi_disc.Tag[i]);
    distr_t_list Pi_plus_eff_mass= Corr.effective_mass_t(pi_plus_distr, "../data/pi_0_TM/pi_plus_eff_mass."+pi_disc.Tag[i]);

   

 
    Pfloat res_charged_pion= Corr.Fit_(Pi_conn_eff_mass);
    Pfloat res_neutral_pion= Corr.Fit_(Pi_0_eff_mass);
    Pfloat res_pi_plus= Corr.Fit_(Pi_plus_eff_mass);
    cout<<"M_pi_plus("<<pi_disc.Tag[i]<<"): "<<res_pi_plus.first<<"("<<res_pi_plus.second<<")"<<endl;
    cout<<"M_pi_0_conn("<<pi_disc.Tag[i]<<"): "<<res_charged_pion.first<<"("<<res_charged_pion.second<<")"<<endl;
    cout<<"M_pi0("<<pi_disc.Tag[i]<<"): "<<res_neutral_pion.first<<"("<<res_neutral_pion.second<<")"<<endl;
    cout<<"##########"<<endl;
  }
  return;
}

void Eta_TM() {


  data_t eta_conn, eta_disc, eta_bubble, pi_plus;

   
  eta_conn.Read("../datasets_eta_TM", "conn_contr", "P5P5");
  eta_disc.Read("../datasets_eta_TM", "disc_contr", "P5");
  eta_bubble.Read("../datasets_eta_TM", "bubble", "P5");
  pi_plus.Read("../datasets_eta_TM", "conn_pi_plus", "P5P5");

  int Nens= eta_conn.size;
  cout<<"N_ens: "<<Nens<<endl;

 

  for(int i=0; i<Nens;i++) {
    distr_t_list eta_conn_distr(UseJack), eta_disc_distr(UseJack), pi_plus_distr(UseJack);
    CorrAnalysis Corr(UseJack, Njacks, Nboots);
    Corr.Perform_Nt_t_average=1;
    Corr.Reflection_sign=1;
    if(eta_conn.Tag[i].substr(0,1)=="A") Corr.Tmin=12;
    else if(eta_conn.Tag[i].substr(0,1)=="B") Corr.Tmin=16;
    else if(eta_conn.Tag[i].substr(0,1)=="D") Corr.Tmin=21;
    else if(eta_conn.Tag[i].substr(0,2)=="a") Corr.Tmin=12;
    else Corr.Tmin=21;
    Corr.Tmax = eta_conn.nrows[i]/2-6;
    Corr.Nt = eta_conn.nrows[i];
    //LatticeInfo L_info("LOCAL");

    //L_info.LatInfo(eta_conn.Tag[i].substr(0,1));
    double L, NNT, m;
    Read_pars_from_ensemble_tag(eta_conn.Tag[i], m, L, NNT);
    boost::filesystem::create_directory("../data");
    boost::filesystem::create_directory("../data/eta_TM");
    eta_conn_distr=Corr.corr_t(eta_conn.col(0)[i], "../data/eta_TM/eta_conn."+eta_conn.Tag[i]);
    pi_plus_distr=Corr.corr_t(pi_plus.col(0)[i], "../data/eta_TM/pi_plus."+eta_conn.Tag[i]);
    
    distr_t_list Eta_conn_eff_mass = Corr.effective_mass_t(eta_conn_distr, "../data/eta_TM/eta_conn_eff_mass."+eta_conn.Tag[i]);
    distr_t_list Pi_plus_eff_mass = Corr.effective_mass_t(pi_plus_distr, "../data/eta_TM/pi_plus_eff_mass."+eta_conn.Tag[i]);

    auto F= [&](const Vfloat& par) { if((signed)par.size() != 2) crash("Lambda function expects par[2], but par["+to_string((signed)par.size())+"] provided"); return par[0] -1.0*par[1]*par[1];};
    Vfloat bubble_average(eta_bubble.col(0)[i][0].size(),0);
    for(int t=0;t<Corr.Nt;t++) {
      if(bubble_average.size() != eta_bubble.col(0)[i][t].size()) crash("Number of confs not costants over time in eta_TM");
      for(unsigned int iconf=0; iconf < bubble_average.size();iconf++) {
	bubble_average[iconf] += eta_bubble.col(0)[i][t][iconf]/Corr.Nt;
      }
    }
    for(int t=0;t<Corr.Nt;t++) {
	   if(UseJack==1) {
	     Jackknife J(10000,Njacks);
	     eta_disc_distr.distr_list.push_back(J.DoJack(F, 2, eta_disc.col(0)[i][t], bubble_average));
	   }
	   else {
	     Bootstrap B(Nboots, 432423, Corr.Nt);
	     eta_disc_distr.distr_list.push_back(B.DoBoot(F, 2, eta_disc.col(0)[i][t], bubble_average));
	   }
    }
    Print_To_File({}, {eta_disc_distr.ave(), eta_disc_distr.err()},"../data/eta_TM/eta_disc."+eta_disc.Tag[i]+".t"  ,"OUT", "");

    
 
    distr_t_list Eta_eff_mass= Corr.effective_mass_t(eta_conn_distr-2.0*eta_disc_distr, "../data/eta_TM/eta_eff_mass."+eta_disc.Tag[i]);
    Pfloat res_eta_conn= Corr.Fit_(Eta_conn_eff_mass);
    Pfloat res_eta= Corr.Fit_(Eta_eff_mass);
    Pfloat res_pi_plus= Corr.Fit_(Pi_plus_eff_mass);
    cout<<"M_pi_plus("<<eta_disc.Tag[i]<<"): "<<res_pi_plus.first<<"("<<res_pi_plus.second<<")"<<endl;
    cout<<"M_eta_conn("<<eta_disc.Tag[i]<<"): "<<res_eta_conn.first<<"("<<res_eta_conn.second<<")"<<endl;
    cout<<"M_eta("<<eta_disc.Tag[i]<<"): "<<res_eta.first<<"("<<res_eta.second<<")"<<endl;
    cout<<"##############"<<endl;
  }
  return;
}

distr_t_list M2_disc(VVfloat& disco, VVfloat& b1, VVfloat& b2) {

   auto F= [&](const Vfloat& par) { if((signed)par.size() != 3) crash("Lambda function expects par[3], but par["+to_string((signed)par.size())+"] provided"); return par[0] -1.0*par[1]*par[2];};

   distr_t_list M2_disc(UseJack);
   distr_t b1_jack(UseJack), b2_jack(UseJack);
   distr_t_list M2_quark_line_disc(UseJack);

  

   if( (b1.size() != b2.size()) || (b1.size() != disco.size())) crash("In M2_disc size of b1,b2 and disco do not coincide");

   int Nt= b1.size();
   if(Nt ==0) crash("M2_disc called with Nt=0");
   if( (b1[0].size() != b2[0].size() ) || (b1[0].size() != disco[0].size() )) crash("In M2_disc, b1,b2 and disco do not have same number of configs");
   
   int Nconfs= b1[0].size();
   Vfloat b1_average(Nconfs,0);
   Vfloat b2_average(Nconfs,0);

   //symmetrize disco
   VVfloat disco_symm(Nt);
   for(int t=0; t<Nt;t++) {
     disco_symm[t].resize(Nconfs,0);
     for(int iconf=0;iconf<Nconfs;iconf++) disco_symm[t][iconf] = (disco[t][iconf] + disco[(Nt-t)%Nt][iconf])/2.0;
   }
      
     for(int t=0;t<Nt;t++) {
       for(int iconf=0; iconf < Nconfs;iconf++) {
	b1_average[iconf] += b1[t][iconf]/Nt;
	b2_average[iconf] += b2[t][iconf]/Nt;
      }
    }


     //call jackknife or bootstrap
      for(int t=0;t<Nt;t++) {
	   if(UseJack==1) {
	     Jackknife J(10000,Njacks);
	     M2_disc.distr_list.push_back(J.DoJack(F, 3, disco_symm[t], b1_average, b2_average));
	     b1_jack= J.DoJack(1,b1_average);
	     b2_jack= J.DoJack(1,b2_average);
	   }
	   else {
	     Bootstrap B(Nboots, 15163, Nt);
	     M2_disc.distr_list.push_back(B.DoBoot(F, 3, disco_symm[t], b1_average, b2_average));
	     b1_jack= B.DoBoot(1, b1_average);
	     b2_jack= B.DoBoot(1, b2_average);
	   }
      }

      cout<<"b1_ave: "<<b1_jack.ave()<<"("<<b1_jack.err()<<")"<<endl;
      cout<<"b2_ave: "<<b2_jack.ave()<<"("<<b2_jack.err()<<")"<<endl;

      
      
      return M2_disc;

}

distr_t_list M1_disc(VVfloat& disco, VVfloat& b1, VVfloat& b2) {

   auto F= [&](const Vfloat& par) { if((signed)par.size() != 3) crash("Lambda function expects par[3], but par["+to_string((signed)par.size())+"] provided"); return par[0] -1.0*par[1]*par[1] +0*par[2]*par[2];};

   distr_t_list M1_disc(UseJack);
   distr_t b1_jack(UseJack), b2_jack(UseJack);
   distr_t_list M1_quark_line_disc(UseJack);

  

   if( (b1.size() != b2.size()) || (b1.size() != disco.size())) crash("In M1_disc size of b1,b2 and disco do not coincide");

   int Nt= b1.size();
   if(Nt ==0) crash("M1_disc called with Nt=0");
   if( (b1[0].size() != b2[0].size() ) || (b1[0].size() != disco[0].size() )) crash("In M1_disc, b1,b2 and disco do not have same number of configs");
   
   int Nconfs= b1[0].size();
   Vfloat b1_average(Nconfs,0);
   Vfloat b2_average(Nconfs,0);

   //symmetrize disco
   VVfloat disco_symm(Nt);
   for(int t=0; t<Nt;t++) {
     disco_symm[t].resize(Nconfs,0);
     for(int iconf=0;iconf<Nconfs;iconf++) disco_symm[t][iconf] = (disco[t][iconf] + disco[(Nt-t)%Nt][iconf])/2.0;
   }
      
     for(int t=0;t<Nt;t++) {
       for(int iconf=0; iconf < Nconfs;iconf++) {
	b1_average[iconf] += b1[t][iconf]/Nt;
	b2_average[iconf] += b2[t][iconf]/Nt;
      }
    }


     //call jackknife or bootstrap
      for(int t=0;t<Nt;t++) {
	   if(UseJack==1) {
	     Jackknife J(10000,Njacks);
	     M1_disc.distr_list.push_back(J.DoJack(F, 3, disco_symm[t], b1_average, b2_average));
	     b1_jack= J.DoJack(1,b1_average);
	     b2_jack= J.DoJack(1,b2_average);
	   }
	   else {
	     Bootstrap B(Nboots, 1655163, Nt);
	     M1_disc.distr_list.push_back(B.DoBoot(F, 3, disco_symm[t], b1_average, b2_average));
	     b1_jack= B.DoBoot(1, b1_average);
	     b2_jack= B.DoBoot(1, b2_average);
	   }
      }

      cout<<"b1_ave: "<<b1_jack.ave()<<"("<<b1_jack.err()<<")"<<endl;
      cout<<"b2_ave: "<<b2_jack.ave()<<"("<<b2_jack.err()<<")"<<endl;

      
      
      return M1_disc;

}


void Axion_l7_analysis() {




  data_t pi_data_r0,pi_data_r1, IE1_data, IE2_data, IH1_data, IH2_data, M2_conn_data, M2_disc_data, M2_b1_data, M2_b2_data, IB_S0P5_data, pi_data_OS;
   pi_data_r0.Read("../datasets_axion_M2", "mes_contr_00_r0", "P5P5");

   IE1_data.Read("../datasets_axion_M2", "mes_contr_IE1", "P5P5");

   IH1_data.Read("../datasets_axion_M2", "mes_contr_IH1", "S0P5");

   IB_S0P5_data.Read("../datasets_axion_M2", "mes_contr_IB_S0P5", "S0P5");
   pi_data_OS.Read("../datasets_axion_M2", "mes_contr_01_r0", "P5P5");
 
   //IS1_data.Read("../datasets_axion", "mes_contr_IS1", "P5P5");


   //############M2##########
   M2_conn_data.Read("../datasets_axion_M2", "M2_conn", "P5P5");
   M2_disc_data.Read("../datasets_axion_M2", "M2_disc", "P5P5");
   M2_b1_data.Read("../datasets_axion_M2", "BUB_S", "S0P5");
   M2_b2_data.Read("../datasets_axion_M2", "BUB_0", "S0P5");
  

  

    
  
  int Nens = pi_data_r0.size;  //== dm_exch_data.size .....
  cout<<"N_ens: "<<Nens<<endl;

  distr_t_list Pi_r0_corr(UseJack),  Exch_r0_corr(UseJack), Hand_r0_corr(UseJack), Pi_OS_corr(UseJack);
  //distr_t_list IS1_corr(UseJack);
  distr_t_list M2_conn_corr(UseJack);
  distr_t_list Pi_plus_corr(UseJack), Pi_plus_eff_mass(UseJack);
  

 
  distr_t_list Exch_r0_eff_slope(UseJack), Hand_r0_eff_slope(UseJack), Global_r0_eff_slope(UseJack);
  
  for(int i=0; i < Nens; i++) {
   
    CorrAnalysis Corr(UseJack, Njacks,Nboots);
    if(pi_data_r0.Tag[i].substr(0,1) == "A")  Corr.Tmin = 12;
      else if(pi_data_r0.Tag[i].substr(0,1) =="B") Corr.Tmin = 13;
      else if(pi_data_r0.Tag[i].substr(0,1) =="a") Corr.Tmin=15;
      else Corr.Tmin= 12;
    Corr.Tmax=pi_data_r0.nrows[i]/2 - 4;
    Corr.Nt = pi_data_r0.nrows[i];
    LatticeInfo L_info("LOCAL");
    
    L_info.LatInfo(pi_data_r0.Tag[i].substr(0,1));
    double L, NNT, m;
    Read_pars_from_ensemble_tag(pi_data_r0.Tag[i], m, L, NNT);
    boost::filesystem::create_directory("../data");
    boost::filesystem::create_directory("../data/axion_l7");
    

    //compute correlators

    Pi_r0_corr = Corr.corr_t(pi_data_r0.col(0)[i], "../data/axion_l7/pi_corr_r0."+pi_data_r0.Tag[i]);
    Pi_OS_corr = Corr.corr_t(pi_data_OS.col(0)[i], "../data/axion_l7/pi_corr_OS."+pi_data_r0.Tag[i]);
    Exch_r0_corr = Corr.corr_t(IE1_data.col(0)[i], "../data/axion_l7/Exch_r0."+pi_data_r0.Tag[i]);
    
    //Hand_r0_corr = Corr.corr_t(IH1_data.col(0)[i], "../data/axion_l7/Hand_r0."+pi_data_r0.Tag[i]);
    Hand_r0_corr= M1_disc(IH1_data.col(0)[i], IB_S0P5_data.col(0)[i], IB_S0P5_data.col(1)[i]);

    Print_To_File({}, {Hand_r0_corr.ave(),Hand_r0_corr.err()}, "../data/axion_l7/Hand_r0."+pi_data_r0.Tag[i]+".t", "", "");

 
    //####################M2####################à   
    
    M2_conn_corr = Corr.corr_t(M2_conn_data.col(0)[i],"");

    M2_conn_corr = M2_conn_corr*-1.0;

    Print_To_File({}, {M2_conn_corr.ave(), M2_conn_corr.err()}, "../data/axion_l7/M2_conn_corr."+pi_data_r0.Tag[i]+".t","","");

    distr_t_list M2_conn_mixed= M2_conn_corr/Pi_r0_corr;

    distr_t_list effective_slope_M2 = Corr.effective_slope_t(M2_conn_corr, Pi_r0_corr, "../data/axion_l7/M2_conn_effective_slope."+pi_data_r0.Tag[i]+".t");

    Print_To_File({},{M2_conn_mixed.ave(), M2_conn_mixed.err()}, "../data/axion_l7/M2_conn_ratio_corr."+pi_data_r0.Tag[i], "", "");

    distr_t_list M2_conn_eff_mass= Corr.effective_mass_t(M2_conn_corr, "../data/axion_l7/M2_conn_effective_mass."+pi_data_r0.Tag[i]);


    
    

    //determine disconnected contribution in M2 method

    distr_t_list M2_disc_tot_corr= -1.0*M2_disc(M2_disc_data.col(0)[i],M2_b1_data.col(0)[i], M2_b2_data.col(0)[i]);

    //combine it with connected contribution and divide by isosymmetric pion correlator

    distr_t_list M2_corr = M2_conn_corr -1.0*M2_disc_tot_corr;

    distr_t_list M2_eff_mass = Corr.effective_mass_t(M2_corr, "../data/axion_l7/M2_effective_mass."+pi_data_r0.Tag[i]);

    //print to file disconnected and total contrib

    Print_To_File({}, {M2_disc_tot_corr.ave(), M2_disc_tot_corr.err()}, "../data/axion_l7/M2_disc_corr."+pi_data_r0.Tag[i]+".t", "", "");
    Print_To_File({}, {M2_corr.ave(), M2_corr.err()}, "../data/axion_l7/M2_corr."+pi_data_r0.Tag[i]+".t", "", "");
     


    //###################M2#################################
    

    

    

    //compute effective masses of neutral and charged pion

    Pi_plus_corr=Pi_r0_corr;
    Pi_plus_eff_mass = Corr.effective_mass_t(Pi_plus_corr, "../data/axion_l7/pi_plus_eff_mass."+pi_data_r0.Tag[i]);
    distr_t_list Pi_OS_eff_mass = Corr.effective_mass_t(Pi_OS_corr, "../data/axion_l7/pi_OS_eff_mass."+pi_data_r0.Tag[i]);
    
    
    auto sqr = [&](double x) -> double { return sqrt(x);};
    auto sqr_t = [&](double x, int t) -> double {return sqrt(x);};
    distr_t_list f0_distr_pi_plus = Corr.decay_constant_t(pow(2.0*m,2)*Pi_plus_corr, "../data/axion_l7/pi_plus_decay_const."+pi_data_r0.Tag[i]);
    distr_t_list f0_distr_pi_OS =  Corr.decay_constant_t(pow(2.0*m,2)*Pi_OS_corr, "../data/axion_l7/pi_OS_decay_const."+pi_data_r0.Tag[i]);
    distr_t_list overlap_pi_plus_distr = Corr.residue_t(Pi_plus_corr, "../data/axion_l7/pi_plus_overlap."+pi_data_r0.Tag[i]);
    distr_t_list sqrt_overlap_pi_plus_distr = distr_t_list::f_of_distr_list(sqr_t, overlap_pi_plus_distr);
    distr_t overlap_pi= Corr.Fit_distr(overlap_pi_plus_distr);
    distr_t sqrt_overlap_pi = distr_t::f_of_distr(sqr, overlap_pi); 
    distr_t f0 = Corr.Fit_distr(f0_distr_pi_plus);
    distr_t Pi_plus_fit = Corr.Fit_distr(Pi_plus_eff_mass);
    distr_t Pi_OS_fit = Corr.Fit_distr(Pi_OS_eff_mass);
    distr_t_list overlap_OS_distr= Corr.residue_t(Pi_OS_corr, "../data/axion_l7/pi_OS_overlap."+pi_data_r0.Tag[i]);
    distr_t_list sqrt_overlap_pi_OS_distr= distr_t_list::f_of_distr_list(sqr_t, overlap_OS_distr);
    distr_t_list Zp_ov_Zs = sqrt_overlap_pi_OS_distr/sqrt_overlap_pi_plus_distr;
    Print_To_File({}, {Zp_ov_Zs.ave(), Zp_ov_Zs.err()}, "../data/axion_l7/Zp_ov_Zs."+pi_data_r0.Tag[i]+".t", "", "");
    distr_t Zp_ov_Zs_fit= Corr.Fit_distr(Zp_ov_Zs);

      
    

    //M2
    distr_t M2_conn_mass_fit=Corr.Fit_distr(M2_conn_eff_mass);
    distr_t M2_mass_fit = Corr.Fit_distr(M2_eff_mass);
    distr_t_list overlap_eta_pi=Corr.residue_t(M2_corr, "../data/axion_l7/eta_pi_overlap."+pi_data_r0.Tag[i]);
    distr_t_list ratio_eta_pi_pi = M2_corr/Pi_r0_corr;
    //M2
       
    
     
    distr_t_list ratio_corr_exch_r0 = Exch_r0_corr/Pi_plus_corr;
    //distr_t_list ratio_corr_exch_r1 = Exch_r1_corr/Pi_0_corr;

    Print_To_File({}, {ratio_corr_exch_r0.ave(), ratio_corr_exch_r0.err()},"../data/axion_l7/ratio_corr_exch_r0."+pi_data_r0.Tag[i]  ,"OUT", "");
    //Print_To_File({}, {ratio_corr_exch_r1.ave(), ratio_corr_exch_r1.err()},"../data/axion_l7/ratio_corr_exch_r1."+pi_data_r0.Tag[i]  ,"OUT", "");

    distr_t_list ratio_corr_hand_r0 = Hand_r0_corr/Pi_plus_corr;
    //distr_t_list ratio_corr_hand_r1 = Hand_r1_corr/Pi_0_corr;

    Print_To_File({}, {ratio_corr_hand_r0.ave(), ratio_corr_hand_r0.err()},"../data/axion_l7/ratio_corr_hand_r0."+pi_data_r0.Tag[i]  ,"OUT", "");
    //Print_To_File({}, {ratio_corr_hand_r1.ave(), ratio_corr_hand_r1.err()},"../data/axion_l7/ratio_corr_hand_r1."+pi_data_r0.Tag[i]  ,"OUT", "");

    
    // distr_t_list ratio_IS1;
    //ratio_IS1 = IS1_corr/Pi_plus_corr;
   
    //distr_t_list eff_mass_IS1;

    //eff_mass_IS1 = Corr.effective_slope_t(IS1_corr, Pi_plus_corr, "../data/axion_l7/IS_r0r0_eff_slope."+pi_data_r0.Tag[i]);
   
    //Print_To_File({}, {ratio_IS1.ave(), ratio_IS1.err()},"../data/axion_l7/ratio_corr_IS_r0r0."+pi_data_r0.Tag[i]  ,"OUT", "");
        
  
    
   
    distr_t ZP_jack(UseJack), ZS_jack(UseJack);
    GaussianMersenne Gauss_ZP(43423,L_info.ZP_err/sqrt(Njacks-1));
    GaussianMersenne Gauss_ZS(65443,L_info.ZS_err/sqrt(Njacks-1));
    for(int nj=0;nj<Njacks;nj++) ZP_jack.distr.push_back( L_info.ZP +Gauss_ZP());
    for(int nj=0;nj<Njacks;nj++) ZS_jack.distr.push_back(L_info.ZS +Gauss_ZS()); 
    
  
     
    Exch_r0_eff_slope = Corr.effective_slope_t(Exch_r0_corr, Pi_plus_corr, "../data/axion_l7/Exch_r0_eff_slope."+pi_data_r0.Tag[i]);

    Hand_r0_eff_slope = Corr.effective_slope_t(Hand_r0_corr, Pi_plus_corr, "../data/axion_l7/Hand_r0_eff_slope."+pi_data_r0.Tag[i]);

    distr_t_list Global_r0 = Exch_r0_corr- Hand_r0_corr;

    distr_t_list ratio_corr_Global_r0 = Global_r0/Pi_plus_corr;

    Print_To_File({}, {Global_r0.ave(), Global_r0.err()}, "../data/axion_l7/Global_r0_corr."+pi_data_r0.Tag[i], "", "");
    Print_To_File({}, {ratio_corr_Global_r0.ave(), ratio_corr_Global_r0.err()},"../data/axion_l7/ratio_corr_Global_r0."+pi_data_r0.Tag[i]  ,"OUT", "");
    

    Global_r0_eff_slope = Corr.effective_slope_t(Exch_r0_corr-Hand_r0_corr, Pi_plus_corr, "../data/axion_l7/Global_r0_eff_slope."+pi_data_r0.Tag[i]);


    int Tmin_old= Corr.Tmin;
    int Tmax_old= Corr.Tmax;
    Corr.Tmin=2;
    Corr.Tmax=30;
    distr_t Global_r0_eff_slope_from_tanh = Corr.effective_slope_t_tanh_fit(Exch_r0_corr-Hand_r0_corr, Pi_plus_corr, "../data/axion_l7/Global_r0_tanh_fit_func."+pi_data_r0.Tag[i]);
    Corr.Tmin=Tmin_old;
    Corr.Tmax=Tmax_old;

    
    
    //fit Global_r0_eff_slope


    Corr.Tmin=5;
    Corr.Tmax=12;
    distr_t Fit_result = Corr.Fit_distr(Global_r0_eff_slope);

    distr_t_list ell_7= 1.0*Global_r0_eff_slope*(1.0/(Zp_ov_Zs_fit*Zp_ov_Zs_fit))*f0*f0*(1.0/(Pi_plus_fit*Pi_plus_fit*Pi_plus_fit))*(m*m);


    //################M2################
    //fit M2 and normalize it properly to extract l7
    distr_t M2_fit_res = -sqrt(2)*Corr.Fit_distr(overlap_eta_pi)*(m*m*f0/sqrt(2))*(1.0/(Pi_plus_fit*Pi_plus_fit*Pi_plus_fit*Pi_plus_fit))*(1.0/sqrt_overlap_pi)/(Zp_ov_Zs_fit) ; //check for possible sqrt(2) !!!!!!
    distr_t M2_fit_res_2 =  -sqrt(2)*Corr.Fit_distr(ratio_eta_pi_pi)*(m*m*f0/sqrt(2))*(1.0/(Pi_plus_fit*Pi_plus_fit*Pi_plus_fit*Pi_plus_fit))*(sqrt_overlap_pi)/(Zp_ov_Zs_fit);
    distr_t_list l7_M2 =  -sqrt(2)*ratio_eta_pi_pi*(m*m*f0/sqrt(2))*(1.0/(Pi_plus_fit*Pi_plus_fit*Pi_plus_fit*Pi_plus_fit))*(sqrt_overlap_pi)/(Zp_ov_Zs_fit);
    distr_t_list l7_M2_new_v= -sqrt(2)*ratio_eta_pi_pi*(m*m*f0/sqrt(2))*(1.0/(Pi_plus_fit*Pi_plus_fit*Pi_plus_fit*Pi_plus_fit))*(sqrt_overlap_pi)/(Zp_ov_Zs_fit);
    Print_To_File({}, {l7_M2_new_v.ave(), l7_M2_new_v.err(), l7_M2.ave(), l7_M2.err()}, "../data/axion_l7/ell_7_M2."+pi_data_r0.Tag[i]+".t", "", "");
    //################M2###############
   
    Print_To_File({}, {ell_7.ave(), ell_7.err()}, "../data/axion_l7/ell_7."+pi_data_r0.Tag[i]+".t", "OUT", "");
    
    //plots
    plt::clf();
    plt::figure_size(1000,1000);
    plt::xlim(1.0, (double)Corr.Nt/2.0 +1);
    plt::xlabel("t");
    plt::ylabel("effective slope");
    Vfloat times;
    for(int tt=0;tt<Corr.Nt;tt++) times.push_back(tt);
    plt::errorbar(times, Exch_r0_eff_slope.ave(),Exch_r0_eff_slope.err(), { {"c", "black"}, {"marker", "."} , {"ls" , "-"}, {"label", "connected"}});
    plt::errorbar(times, Global_r0_eff_slope.ave(), Global_r0_eff_slope.err(), { {"c", "red"}, {"marker", "."} , {"ls" , "-"}, {"label", "global"}});
    plt::errorbar(times, Hand_r0_eff_slope.ave(), Hand_r0_eff_slope.err(), { {"c", "green"}, {"marker", "."}, {"ls", "-"}, {"label", "disconnected"}});
    plt::legend();
    string path="../plots/axion_l7/eff_slopes_"+pi_data_r0.Tag[i]+".png";
    plt::save(path.c_str());
    plt::close();


    
    
    cout<<"#####################"<<endl;
    cout<<"Ensemble: "<<pi_data_r0.Tag[i]<<endl;
    //try to reconstruct l7
    cout<<"m: "<<m<<endl;
    cout<<"Zp/Zs (FRA): "<<L_info.ZP/L_info.ZS<<endl;
    cout<<"Zp/Zs: "<<Zp_ov_Zs_fit.ave()<<"("<<Zp_ov_Zs_fit.err()<<")"<<endl;
    cout<<"pi+ mass: "<<Pi_plus_fit.ave()<<" "<<Pi_plus_fit.err()<<endl;
    cout<<"pi OS mass: "<<Pi_OS_fit.ave()<<" "<<Pi_OS_fit.err()<<endl;
    cout<<"f0 (charged pion): "<<f0.ave()<<"("<<f0.err()<<")"<<endl;
    cout<<"Z_pipi (charged pion): "<<sqrt_overlap_pi.ave()<<"("<<sqrt_overlap_pi.err()<<")"<<endl;
    //cout<<"dm: "<<Fit_result.ave()<<"("<<Fit_result.err()<<")"<<endl;
    distr_t l7 =  1.0*(1.0/(Zp_ov_Zs_fit*Zp_ov_Zs_fit))*f0*f0*Fit_result*(1.0/(Pi_plus_fit*Pi_plus_fit*Pi_plus_fit))*(m*m);
    distr_t l7_tanh=  1.0*(1.0/(Zp_ov_Zs_fit*Zp_ov_Zs_fit))*f0*f0*Global_r0_eff_slope_from_tanh*(1.0/(Pi_plus_fit*Pi_plus_fit*Pi_plus_fit))*(m*m);
   
    cout<<"l7(M1): "<<l7.ave()<<"("<<l7.err()<<")"<<endl;
    cout<<"l7(M1) [tanh fit]: "<<l7_tanh.ave()<<"("<<l7_tanh.err()<<")"<<endl;

    //###############M2##################
    cout<<"M2_conn_mass: "<<M2_conn_mass_fit.ave()<<" "<<M2_conn_mass_fit.err()<<endl;
    cout<<"M2_mass: "<<M2_mass_fit.ave()<<" "<<M2_mass_fit.err()<<endl;
    cout<<"l7(M2): "<<M2_fit_res.ave()<<"("<<M2_fit_res.err()<<")"<<endl;
    cout<<"l7(M2) new version: "<<M2_fit_res_2.ave()<<"("<<M2_fit_res_2.err()<<")"<<endl;
    //###############M2################à#
    cout<<"#####################"<<endl;
     
  }


  
     
  return;
}
 
					
  






















