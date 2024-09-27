#include "../include/HVP.h"


const double DTT = 0.5;
const double alpha = 1.0/137.035999;	     

using namespace std;

void Bounding_HVP(distr_t &amu_HVP, int &Tcut_opt,  const distr_t_list &V, const distr_t &a, string path,distr_t lowest_mass) {


  int Njacks=158;
  bool UseJack=true;
 
  int Simps_ord=2;
  double fm_to_inv_Gev= 1.0/0.197327;


  double T= V.size(); 

  auto LOG = [](double R_G, double t) { return log(fabs(R_G));};

  auto K = [&](double Mv, double t, double size) -> double { return kernel_K(t, Mv);};

  distr_t_list Ker = distr_t_list::f_of_distr(K, a , T/2);

  int T_ext_max= 300;

  distr_t_list Ker_extended= distr_t_list::f_of_distr(K, a , T_ext_max);

  int Tdatas_opt=-1;

  
  auto amu_HVP_func = [&V, &Ker, &UseJack, &Njacks,  &Simps_ord](double tcut) -> distr_t {
    
    distr_t ret_win(UseJack,UseJack?Njacks:800);
    for(int t=1;t<tcut;t++) ret_win = ret_win + 4.0*w(t,Simps_ord)*pow(alpha,2)*V.distr_list[t]*Ker.distr_list[t];
    return ret_win;
    
  };

  auto amu_HVP_int = [ &Ker_extended, &Simps_ord](double t) -> distr_t {

    return 4.0*w(t,Simps_ord)*pow(alpha,2)*Ker_extended.distr_list[t];

  };
  
 
 
  CorrAnalysis Corr(UseJack, Njacks,800);
  Corr.Nt = V.size();
  distr_t_list ratio_corr_V(UseJack);
  for(int t=0; t<Corr.Nt;t++) ratio_corr_V.distr_list.push_back( V.distr_list[t]/V.distr_list[(t+1)%Corr.Nt]);
  //decide whether to use m_eff(t) or log( V(t)/V(t+1)) (true = use log)  
  bool update_min_Tdata_from_log_ratio=true;    
  distr_t_list eff_mass_V = (update_min_Tdata_from_log_ratio)?distr_t_list::f_of_distr_list(LOG, ratio_corr_V):Corr.effective_mass_t(V, "");

  

  distr_t_list amu_HVP_min_Tdata(UseJack), amu_HVP_max_Tdata(UseJack), amu_HVP_T_2(UseJack);

  Vfloat TCUTS;

  Vfloat Is_T_data_opt;
  int slice_to_use_for_eff_mass=1;
    
    
  //loop over tcut
  for(int tcut=1; tcut<T/2;tcut++) {


    TCUTS.push_back( (double)tcut);
    distr_t amu_HVP_up_to_tcut = amu_HVP_func(tcut+1);
    amu_HVP_min_Tdata.distr_list.push_back(amu_HVP_up_to_tcut);
    amu_HVP_max_Tdata.distr_list.push_back(amu_HVP_up_to_tcut);
    amu_HVP_T_2.distr_list.push_back(amu_HVP_up_to_tcut);
    distr_t V_tcut = V.distr_list[tcut];
    Is_T_data_opt.push_back( 0.0);
      
    bool eff_mass_is_nan= isnan( eff_mass_V.ave(tcut));
    bool update_min_Tdata=true;
    if(eff_mass_is_nan || (eff_mass_V.err(tcut)/eff_mass_V.ave(tcut) > 0.10) || (eff_mass_V.ave(tcut) < 0 )) update_min_Tdata=false;

    if(update_min_Tdata) slice_to_use_for_eff_mass=tcut;
      
      

    for(int t=tcut+1; t < T_ext_max;t++) {

      //lambda function for lower and upper limit of single exp V(t)
      auto EXP_MIN = [&tcut, &t](double E) { return exp(-E*(t-tcut));};
	
      distr_t ker_val = V_tcut*amu_HVP_int(t);
      int size_min= amu_HVP_min_Tdata.size();
      int size_max= amu_HVP_max_Tdata.size();
      if(size_min != tcut || size_max != tcut) crash("size_min or size_max is different from tcut");
      distr_t lower_exp = distr_t::f_of_distr(EXP_MIN,eff_mass_V.distr_list[slice_to_use_for_eff_mass]);
      amu_HVP_min_Tdata.distr_list[size_min-1] = amu_HVP_min_Tdata.distr_list[size_min-1] + ker_val*lower_exp;
      amu_HVP_max_Tdata.distr_list[size_max-1] = amu_HVP_max_Tdata.distr_list[size_max-1] + ker_val*distr_t::f_of_distr(EXP_MIN, lowest_mass);
    }

      
  }


   


  //find interval where difference between amu_HVP_min_Tdata and amu_HVP_max_Tdata is smaller than 0.2 sigma
  //average
  bool Found_Tdata_opt=false;
  int tdata_opt=1;

   
  while(!Found_Tdata_opt && tdata_opt < T/2) {

    distr_t diff_max_min = amu_HVP_max_Tdata.distr_list[tdata_opt-1] - amu_HVP_min_Tdata.distr_list[tdata_opt-1];
    if( fabs(diff_max_min.ave())/min( amu_HVP_max_Tdata.err(tdata_opt-1) , amu_HVP_min_Tdata.err(tdata_opt-1)) < 0.5) Found_Tdata_opt=true;
    else tdata_opt++;
  }

      

  //if tdata_opt has not been found return tdata = -1 and amu_HVP =amu_HVP(T/2)
  if(!Found_Tdata_opt) {
	
    Tdatas_opt = -1;
    amu_HVP =  amu_HVP_func(T/2);
	
  }
  else { //tdata_opt has been found


    bool fit_to_constant=true;
    //average over min max over an interval of 1 fm
    if( amu_HVP_max_Tdata.size() != amu_HVP_min_Tdata.size()) crash("Size of amu_HVP_max_Tdata and amu_HVP_min_Tdata are not equal");
    int S_Tdata= amu_HVP_min_Tdata.size();
    int DT = min( (int)(DTT*fm_to_inv_Gev/a.ave()), (S_Tdata -tdata_opt));
    distr_t amu_HVP_ave_max= amu_HVP_max_Tdata.distr_list[tdata_opt-1];
    distr_t amu_HVP_ave_min= amu_HVP_min_Tdata.distr_list[tdata_opt-1];
    if(fit_to_constant) {
      amu_HVP_ave_max = amu_HVP_ave_max*1.0/pow(amu_HVP_max_Tdata.err(tdata_opt-1),2);
      amu_HVP_ave_min = amu_HVP_ave_min*1.0/pow(amu_HVP_min_Tdata.err(tdata_opt-1),2);
    }
    double amu_HVP_weight_max = (fit_to_constant)?1.0/(pow(amu_HVP_max_Tdata.err(tdata_opt-1),2)):1.0;
    double amu_HVP_weight_min = (fit_to_constant)?1.0/(pow(amu_HVP_min_Tdata.err(tdata_opt-1),2)):1.0;
    for(int t=tdata_opt; t< tdata_opt+DT; t++) {

      double weight_max_t = (fit_to_constant)?1.0/pow(amu_HVP_max_Tdata.err(t),2):1.0;
      double weight_min_t = (fit_to_constant)?1.0/pow(amu_HVP_min_Tdata.err(t),2):1.0;
	 
      amu_HVP_ave_max = amu_HVP_ave_max + amu_HVP_max_Tdata.distr_list[t]*weight_max_t;
      amu_HVP_ave_min = amu_HVP_ave_min + amu_HVP_min_Tdata.distr_list[t]*weight_min_t;
      amu_HVP_weight_max += weight_max_t;
      amu_HVP_weight_min += weight_min_t;

    }

    amu_HVP_ave_max= amu_HVP_ave_max/amu_HVP_weight_max;
    amu_HVP_ave_min= amu_HVP_ave_min/amu_HVP_weight_min;


    double den_weight= 1.0/(pow(amu_HVP_ave_max.err(),2))  + 1.0/(pow(amu_HVP_ave_min.err(),2));
    double fin_weight_max= (fit_to_constant)?(1.0/(pow(amu_HVP_ave_max.err(),2)*den_weight)):0.5;
    double fin_weight_min= (fit_to_constant)?(1.0/(pow(amu_HVP_ave_min.err(),2)*den_weight)):0.5;

    Tdatas_opt =tdata_opt;
    amu_HVP = fin_weight_min*amu_HVP_ave_min + fin_weight_max*amu_HVP_ave_max;
    Is_T_data_opt[tdata_opt-1] = 1.0;

  }

  distr_t_list amu_HVP_VEC(UseJack);
  vector<double> Is_in_fit_int;
  for(int i=0;i<(signed)TCUTS.size(); i++) {
    amu_HVP_VEC.distr_list.push_back( amu_HVP);
    Is_in_fit_int.push_back(( (i >= Tdatas_opt) && (i <= (Tdatas_opt+  DTT*fm_to_inv_Gev/a.ave())))?1.0:0.0);
  }

  Tcut_opt= Tdatas_opt;
      
  //Print to File
  Print_To_File({}, {TCUTS, amu_HVP_min_Tdata.ave(), amu_HVP_min_Tdata.err(), amu_HVP_max_Tdata.ave(), amu_HVP_max_Tdata.err(), amu_HVP_T_2.ave(), amu_HVP_T_2.err(), amu_HVP_VEC.ave(), amu_HVP_VEC.err(), Is_in_fit_int, Is_T_data_opt}, path+".bound", "", "#tcut   lower    upper     T/2  result   Is_in_fit_int?     Is_Tcut_opt.       Tcut_opt= "+to_string(Tdatas_opt));



  return;

}








void HVP() {

  Generate_free_corr_data();
  exit(-1);


  int Njacks=158;
  bool UseJack=true;
  double qu= 2.0/3.0;
  double qd= -1.0/3.0;
  double fm_to_inv_Gev= 1.0/0.197327;


  int NB=20000;

  double GF=  1.1663787*1e-5;
  double hbar= 6.582119569*1e-25;
  double mt= 1.77686;
  double mk= 0.493677;
  double mpi=0.13957039;
  double mpi0=0.1349768;
  distr_t B_TK(false);
  distr_t B_TP(false);
  distr_t FK(false);
  distr_t FP(false);
  distr_t SEW(false);
  distr_t dR_TK(false);
  distr_t dR_TP(false);
  distr_t LT_tau(false);
  GaussianMersenne GV(645611);
  distr_t Vus(false);
  distr_t Vud(false);


  /*
  //Gounaris Sakurai model

  int npts_spline= 400;
  int Luscher_num_zeroes= 22;
  int Nresonances=20;
  //Init LL_functions;
  //find first  zeros of the Lusher functions
  Vfloat Luscher_zeroes;
  Zeta_function_zeroes(Luscher_num_zeroes, Luscher_zeroes);

  //############################################INTERPOLATE PHI FUNCTION AND DERIVATIVES#############################

  
  VVfloat phi_data, phi_der_data;
  Vfloat sx_int;
  Vfloat sx_der, dx_der;
  Vfloat Dz;

  for(int L_zero=0;L_zero<Nresonances+1;L_zero++) {
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
    double dz = (dx-sx)/npts_spline;
    Dz.push_back(dz);


    for(int istep=1;istep<=npts_spline-1;istep++) { double pt= sx+dz*istep; phi_data[L_zero].push_back( phi(sqrt(pt)));}

    phi_data[L_zero].push_back(M_PI/2.0);
    double sx_der_loc =  phi_der_for_back(sqrt(sx)+1e-14, 1);
    double dx_der_loc =  phi_der_for_back(sqrt(dx)-1e-14, -1);
    sx_der.push_back(sx_der_loc);
    dx_der.push_back(dx_der_loc);

    phi_der_data[L_zero].push_back(sx_der_loc);
    for(int istep=1;istep<=npts_spline-1;istep++) { double pt= sx+dz*istep; phi_der_data[L_zero].push_back( phi_der(sqrt(pt)));}
    phi_der_data[L_zero].push_back(dx_der_loc);
    
  }



 
   

  LL_functions LL(phi_data,phi_der_data,sx_der, dx_der, sx_int, Dz, Nresonances, Luscher_zeroes);
    
  //###########################################END INTERPOLATION PHI FUNCTION AND DERIVATIVES################################
  cout<<"####Spline for phi(x) and phi'(x) successfully generated!"<<endl;

  double vol1= 5.1*fm_to_inv_Gev;
  double vol2= 9*fm_to_inv_Gev;

  Vfloat En1;
  LL.Find_pipi_energy_lev(vol1, 0.770, 5.5, 0.140, 0.0, En1);

  Vfloat En2;
  LL.Find_pipi_energy_lev(vol2, 0.770, 5.5, 0.140, 0.0, En2);
  ofstream PGS("../data/GS_Humb.dat");
  for(int i=0; i<Nresonances;i++) PGS<<2.0*sqrt( pow(En1[i],2) + pow(0.139,2))<<" "<<LL.Amplitude(En1[i],vol1 , 0.770, 5.5, 0.140, 0.0)<<" "<<2.0*sqrt( pow(En2[i],2) + pow(0.140,2))<<" "<<LL.Amplitude(En2[i],vol2 , 0.770, 5.5, 0.140, 0.0)<<endl;

  PGS.close();

  Vfloat epsilons({0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0.05});

  int Num_ergs=400;
  Vfloat Ergs;
  for(int i=0;i<Num_ergs;i++) Ergs.push_back( i*2.0/(Num_ergs-1.0));

  VVfloat  IM_H1(epsilons.size()), IM_H2(epsilons.size());

  for(int i=0; i<(signed)epsilons.size(); i++) {
    double eps = epsilons[i];

    for(int j=0;j<Num_ergs;j++) {
      double res1=0;
      double res2=0;
      for(int ires=0;ires<Nresonances;ires++) {

	double E1= 2.0*sqrt( pow(En1[ires],2) + pow(0.139,2));
	double E2= 2.0*sqrt( pow(En2[ires],2) + pow(0.139,2)); 
	
	res1 += LL.Amplitude(En1[ires],vol1 , 0.770, 5.5, 0.140, 0.0)*(eps/2.0)*(1.0/( pow((Ergs[j] - E1),2) + eps*eps));
	res2 += LL.Amplitude(En2[ires],vol2 , 0.770, 5.5, 0.140, 0.0)*(eps/2.0)*(1.0/( pow((Ergs[j] - E2),2) + eps*eps));
      }

      IM_H1[i].push_back(res1);
      IM_H2[i].push_back(res2);
    }
    
    Print_To_File({}, {Ergs, IM_H1[i], IM_H2[i]}, "../data/eps_"+to_string(i)+".dat", "", "");

  }


  //print correlator
  Vfloat ts;
  for(int t=0;t<101;t++) ts.push_back( 2.0 + 2.0*t/100.0);

  Vfloat corr;
  for(int t=0;t<101;t++) corr.push_back( (10.0/9.0)*4.0*pow(alpha,2)*LL.V_pipi(ts[t]*fm_to_inv_Gev, vol1, 0.770, 5.5, 0.140, 0.0,  En1));
    
  for (int t=0;t<101;t++) {
    cout<<ts[t]<<" "<<corr[t]<<endl;
  }

  exit(-1);

  
  */
  

  for(int iboot=0;iboot<NB;iboot++) {
    B_TK.distr.push_back(   (6.957 + GV()*0.096)*1e-3);
    B_TP.distr.push_back(   (0.10808 + GV()*0.00053));
    FK.distr.push_back( 0.1557 + GV()*0.0003);
    FP.distr.push_back( 0.1302 + GV()*0.0008);
    SEW.distr.push_back( 1.0232 + GV()*0.0003);
    LT_tau.distr.push_back( (2.903+ GV()*0.005)*1e-13);
    dR_TK.distr.push_back(  1 + (-0.15 + GV()*0.57)*1e-2);
    dR_TP.distr.push_back(  1 + (-0.24 + GV()*0.56)*1e-2);

  }
  
  Vus = GF*GF*(1.0/(16*M_PI*hbar))*FK*FK*LT_tau*pow(mt,3)*pow( 1 - pow(mk/mt,2),2)*SEW*dR_TK/B_TK;

  Vus= 1.0/Vus;
  Vus= SQRT_D(Vus);

  Vud = GF*GF*(1.0/(16*M_PI*hbar))*FP*FP*LT_tau*pow(mt,3)*pow( 1 - pow(mpi/mt,2),2)*SEW*dR_TP/B_TP;

  Vud= 1.0/Vud;
  Vud= SQRT_D(Vud);


  distr_t FK_T_Vus_2 =  (1.0/( GF*GF*(1.0/(16*M_PI*hbar))*LT_tau*pow(mt,3)*pow( 1 - pow(mk/mt,2),2)*SEW*dR_TK/B_TK));

  cout<<"|f_k * Vus |^2_tau : "<<FK_T_Vus_2.ave()<<" +- "<<FK_T_Vus_2.err()<<endl;


  cout<<"Vus: "<<Vus.ave()<<" +- "<<Vus.err()<<endl;
  cout<<"Vud: "<<Vud.ave()<<" +- "<<Vud.err()<<endl;

  //exit(-1);
   

  

  
  

  

 
  //create directories
  boost::filesystem::create_directory("../data/HVP");
  boost::filesystem::create_directory("../data/HVP/Bounding");
  boost::filesystem::create_directory("../data/HVP/eff_mass");
  boost::filesystem::create_directory("../data/HVP/Corr");

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


  auto Sort_old_light_confs = [](string A, string B) {

			   

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
  
  
  data_t Vk_data_tm, Vk_data_OS, P5P5_data_tm;

  Vk_data_tm.Read("../R_ratio_data/light", "mes_contr_2pts_ll_1", "VKVK", Sort_light_confs);
  P5P5_data_tm.Read("../R_ratio_data/light", "mes_contr_2pts_ll_1", "P5P5", Sort_light_confs);
  Vk_data_OS.Read("../R_ratio_data/light", "mes_contr_2pts_ll_2", "VKVK", Sort_light_confs);

  data_t Vk_data_defl_tm, Vk_data_defl_OS;

  Vk_data_defl_tm.Read("../new_gm2_project/semid_new_nev400_nts1024/AkAk_r+1_r-1", "corr", "" ,  Sort_old_light_confs);
  Vk_data_defl_OS.Read("../new_gm2_project/semid_new_nev400_nts1024/VkVk_r+1_r+1", "corr", "" ,  Sort_old_light_confs);
  


  //############################################################################################
  //generate fake jack_distr for lattice spacing a_A a_B, a_C, a_D and RENORMALIZATION CONSTANT
  GaussianMersenne GM(36551294);
  LatticeInfo a_info;
  distr_t a_A(UseJack), a_B(UseJack), a_C(UseJack), a_D(UseJack), a_Z(UseJack), a_E(UseJack);
  distr_t ZV_A(UseJack), ZV_B(UseJack), ZV_C(UseJack), ZV_D(UseJack);
  distr_t ZA_A(UseJack), ZA_B(UseJack), ZA_C(UseJack), ZA_D(UseJack);
  double a_A_ave, a_A_err, a_B_ave, a_B_err, a_C_ave, a_C_err, a_D_ave, a_D_err, a_Z_ave, a_Z_err, a_E_ave, a_E_err;
  double ZV_A_ave, ZV_A_err, ZV_B_ave, ZV_B_err, ZV_C_ave, ZV_C_err, ZV_D_ave, ZV_D_err;
  double ZA_A_ave, ZA_A_err, ZA_B_ave, ZA_B_err, ZA_C_ave, ZA_C_err, ZA_D_ave, ZA_D_err;
  a_info.LatInfo_new_ens("cA211a.53.24");
  a_A_ave= a_info.a_from_afp_FLAG;
  a_A_err= a_info.a_from_afp_FLAG_err;
  ZA_A_ave = a_info.Za_WI_strange;
  ZA_A_err = a_info.Za_WI_strange_err;
  ZV_A_ave = a_info.Zv_WI_strange;
  ZV_A_err = a_info.Zv_WI_strange_err;
  a_info.LatInfo_new_ens("cB211b.072.64");
  a_B_ave= a_info.a_from_afp_FLAG;
  a_B_err= a_info.a_from_afp_FLAG_err;
  ZA_B_ave = a_info.Za_WI_strange;
  ZA_B_err = a_info.Za_WI_strange_err;
  ZV_B_ave = a_info.Zv_WI_strange;
  ZV_B_err = a_info.Zv_WI_strange_err;
  a_info.LatInfo_new_ens("cC211a.06.80");
  a_C_ave= a_info.a_from_afp_FLAG;
  a_C_err= a_info.a_from_afp_FLAG_err;
  ZA_C_ave = a_info.Za_WI_strange;
  ZA_C_err = a_info.Za_WI_strange_err;
  ZV_C_ave = a_info.Zv_WI_strange;
  ZV_C_err = a_info.Zv_WI_strange_err;
  a_info.LatInfo_new_ens("cD211a.054.96");
  a_D_ave= a_info.a_from_afp_FLAG;
  a_D_err= a_info.a_from_afp_FLAG_err;
  ZA_D_ave = a_info.Za_WI_strange;
  ZA_D_err = a_info.Za_WI_strange_err;
  ZV_D_ave = a_info.Zv_WI_strange;
  ZV_D_err = a_info.Zv_WI_strange_err;
  a_info.LatInfo_new_ens("cZ211a.077.64");
  a_Z_ave= a_info.a_from_afp_FLAG;
  a_Z_err= a_info.a_from_afp_FLAG_err;
  a_info.LatInfo_new_ens("cE211a.044.112");
  a_E_ave= a_info.a_from_afp_FLAG;
  a_E_err= a_info.a_from_afp_FLAG_err;

  for(int ijack=0;ijack<Njacks;ijack++) {

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

  }
      

  //############################################################################################




  
  VVfloat C_tm, C_OS;


  //###################################################
  //###########    READING  BINARY FILE ################

  
  //read binary
  vector<string> Corr_tags({"VkVk_OS.bin", "VkVk_tm.bin"});
  
  
  for(int id=0; id<(signed)Corr_tags.size(); id++) {

    
    string ch= Corr_tags[id];
	
    FILE *stream = fopen(("../blinding_test/cB211b.072.64/"+ch).c_str(), "rb");
    size_t Nconfs, TT, Nhits, Nsubs;
    bin_read(Nconfs, stream);
    bin_read(Nhits, stream);
    bin_read(TT, stream);
    bin_read(Nsubs, stream);

    if(id==0) {
      C_OS.resize(TT);
      for(auto & cc: C_OS) cc.resize(Nconfs,0);
    }
    else {
      C_tm.resize(TT);
      for(auto & cc: C_tm) cc.resize(Nconfs,0);
    }
   
     
    
   
     cout<<"Nconfs: "<<Nconfs<<endl;
     cout<<"TT: "<<TT<<" "<<TT/2+1<<endl;
     cout<<"Nhits: "<<Nhits<<endl;
     cout<<"Nsubs: "<<Nsubs<<endl;
    
  
    
    
    for(size_t iconf=0;iconf<Nconfs;iconf++) {
      for(size_t t=0;t<TT/2+1;t++) {
	for(size_t is=0;is<Nsubs;is++) {
	  
	  double c;
	  bin_read(c, stream);
	  if(id==0) {
	    C_OS[t][iconf] += c/Nsubs;
	  }
	  else {
	    C_tm[t][iconf] += c/Nsubs;
	  }
	}
	//symmetrize
	if( t != 0) {
	  if(id==0) {
	    C_OS[TT-t][iconf] = C_OS[t][iconf];
	  }
	  else {
	    C_tm[TT-t][iconf] = C_tm[t][iconf];
	  }
	}
      }
    }
    
   
    
    fclose(stream);
    
  }
  
  
  CorrAnalysis Corr_BL(UseJack, Njacks,100);
  Corr_BL.Nt = 128;
  
  boost::filesystem::create_directory("../data/test_blinding");
  
  distr_t_list C_tm_BLIND= Corr_BL.corr_t( C_tm, "../data/test_blinding/C_tm");
  distr_t_list C_OS_BLIND= Corr_BL.corr_t( C_OS, "../data/test_blinding/C_OS");
  
  
  

  exit(-1);




  //################################## E ENSEMBLE ANALYSIS ####################################


  bool Get_ASCII= true;

    if(Get_ASCII) {
    //read binary files
    boost::filesystem::create_directory("../E112_RCs");
    

    vector<string> Ens_T1({"E.44.112"});
    vector<string> Ens_TT1({"cE211a.044.112"});

    for( int it=0; it<(signed)Ens_T1.size(); it++) {

      vector<string> channels({"mix_l_l",  "mix_l_s1", "mix_l_s2", "mix_s1_s1", "mix_s2_s2"});

      for(auto &channel : channels) {
	boost::filesystem::create_directory("../E112_RCs/"+channel);
	boost::filesystem::create_directory("../E112_RCs/"+channel+"/"+Ens_TT1[it]);
      }
      //read binary
      vector<string> Corr_tags({"TM_VKVK", "TM_AKAK", "TM_A0A0", "TM_V0V0", "TM_P5P5", "TM_A0P5", "OS_VKVK", "OS_AKAK", "OS_A0A0", "OS_V0V0", "OS_P5P5", "OS_A0P5"});

          
      for(int id=0; id<(signed)Corr_tags.size(); id++) {

	cout<<"Corr: "<<Corr_tags[id]<<endl;
	for( auto &channel: channels) {

	FILE *stream = fopen( ("../gm2_tau_rep_bin/"+Ens_T1[it]+"/"+channel+"_"+Corr_tags[id]).c_str(), "rb");
        size_t Nconfs, T, Nhits;
	bin_read(Nconfs, stream);
	bin_read(Nhits, stream);
	bin_read(T, stream);
	cout<<"channel: "<<channel<<endl;
	cout<<"Nconfs: "<<Nconfs<<endl;
	cout<<"T: "<<T<<" "<<T/2+1<<endl;
	cout<<"Nhits: "<<Nhits<<endl;
	for(size_t iconf=0;iconf<Nconfs;iconf++) {
	  vector<double> C(T/2+1);
	  for(size_t t=0;t<T/2+1;t++) bin_read(C[t], stream);
	  boost::filesystem::create_directory("../E112_RCs/"+channel+"/"+Ens_TT1[it]+"/"+to_string(iconf));
	  ofstream PrintCorr("../E112_RCs/"+channel+"/"+Ens_TT1[it]+"/"+to_string(iconf)+"/mes_contr_"+channel+"_"+Corr_tags[id]);
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


    bool Get_ASCII_charm= false;

    if(Get_ASCII_charm) {
    //read binary files
    boost::filesystem::create_directory("../charm_E");
    

    vector<string> Ens_T1({"E.44.112"});
    vector<string> Ens_TT1({"cE211a.044.112"});

    for( int it=0; it<(signed)Ens_T1.size(); it++) {

      vector<string> channels({"mix_s_c1",  "mix_s_c2", "mix_s_c3"});

      for(auto &channel : channels) {
	boost::filesystem::create_directory("../charm_E/"+channel);
	boost::filesystem::create_directory("../charm_E/"+channel+"/"+Ens_TT1[it]);
      }
      //read binary
      vector<string> Corr_tags({"TM_P5P5", "OS_P5P5"});

          
      for(int id=0; id<(signed)Corr_tags.size(); id++) {

	cout<<"Corr: "<<Corr_tags[id]<<endl;
	for( auto &channel: channels) {

	FILE *stream = fopen( ("../charm_E_bin/"+Ens_T1[it]+"/"+channel+"_"+Corr_tags[id]).c_str(), "rb");
        size_t Nconfs, T, Nhits;
	bin_read(Nconfs, stream);
	bin_read(Nhits, stream);
	bin_read(T, stream);
	cout<<"channel: "<<channel<<endl;
	cout<<"Nconfs: "<<Nconfs<<endl;
	cout<<"T: "<<T<<" "<<T/2+1<<endl;
	cout<<"Nhits: "<<Nhits<<endl;
	for(size_t iconf=0;iconf<Nconfs;iconf++) {
	  vector<double> C(T/2+1);
	  for(size_t t=0;t<T/2+1;t++) bin_read(C[t], stream);
	  boost::filesystem::create_directory("../charm_E/"+channel+"/"+Ens_TT1[it]+"/"+to_string(iconf));
	  ofstream PrintCorr("../charm_E/"+channel+"/"+Ens_TT1[it]+"/"+to_string(iconf)+"/mes_contr_"+channel+"_"+Corr_tags[id]);
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
   

    

    data_t V_strange_L, V_strange_OS_L, V_strange_M, V_strange_OS_M;
    data_t corr_P5P5_strange, corr_P5P5_OS_strange, corr_P5P5_strange_heavy, corr_P5P5_OS_strange_heavy;
    data_t corr_A0P5_strange, corr_A0P5_OS_strange, corr_A0P5_strange_heavy, corr_A0P5_OS_strange_heavy;
    //strange
    //L
    V_strange_L.Read("../E112_RCs/mix_s1_s1", "mes_contr_mix_s1_s1_TM_VKVK", "VKVK", Sort_easy);
    V_strange_OS_L.Read("../E112_RCs/mix_s1_s1", "mes_contr_mix_s1_s1_OS_VKVK", "VKVK", Sort_easy);
    //M
    V_strange_M.Read("../E112_RCs/mix_s2_s2", "mes_contr_mix_s2_s2_TM_VKVK", "VKVK", Sort_easy); 
    V_strange_OS_M.Read("../E112_RCs/mix_s2_s2", "mes_contr_mix_s2_s2_OS_VKVK", "VKVK", Sort_easy);
    //P5P5
    corr_P5P5_strange.Read("../E112_RCs/mix_s1_s1", "mes_contr_mix_s1_s1_TM_P5P5", "P5P5", Sort_easy);
    corr_P5P5_OS_strange.Read("../E112_RCs/mix_s1_s1", "mes_contr_mix_s1_s1_OS_P5P5", "P5P5", Sort_easy);
    corr_P5P5_strange_heavy.Read("../E112_RCs/mix_s2_s2", "mes_contr_mix_s2_s2_TM_P5P5", "P5P5", Sort_easy);
    corr_P5P5_OS_strange_heavy.Read("../E112_RCs/mix_s2_s2", "mes_contr_mix_s2_s2_OS_P5P5","P5P5",Sort_easy);
    //A0P5
    corr_A0P5_strange.Read("../E112_RCs/mix_s1_s1", "mes_contr_mix_s1_s1_TM_A0P5", "A0P5", Sort_easy);
    corr_A0P5_OS_strange.Read("../E112_RCs/mix_s1_s1", "mes_contr_mix_s1_s1_OS_A0P5", "A0P5", Sort_easy);
    corr_A0P5_strange_heavy.Read("../E112_RCs/mix_s2_s2", "mes_contr_mix_s2_s2_TM_A0P5", "A0P5", Sort_easy);
    corr_A0P5_OS_strange_heavy.Read("../E112_RCs/mix_s2_s2", "mes_contr_mix_s2_s2_OS_A0P5", "A0P5", Sort_easy);
    //VKVK light
    data_t V_light_tm, V_light_OS;
    V_light_tm.Read("../E112_RCs/mix_l_l", "mes_contr_mix_l_l_TM_VKVK", "VKVK", Sort_easy);
    V_light_OS.Read("../E112_RCs/mix_l_l", "mes_contr_mix_l_l_OS_VKVK", "VKVK", Sort_easy);
    int id_E=0;
    
    CorrAnalysis Corr(UseJack, Njacks,100);
    Corr.Nt = V_strange_L.nrows[id_E];

    boost::filesystem::create_directory("../data/RC_analysis_E");

    
    data_t c1_s, c2_s, c3_s  ;


    distr_t aE(UseJack);
    for(int ijack=0;ijack<Njacks;ijack++) { aE.distr.push_back( 0.0489042 + GM()*6.01389e-05/sqrt(Njacks-1.0));}

    aE = aE*fm_to_inv_Gev;

    c1_s.Read("../charm_E/mix_s_c1", "mes_contr_mix_s_c1_TM_P5P5", "P5P5", Sort_easy);
    c2_s.Read("../charm_E/mix_s_c2", "mes_contr_mix_s_c2_TM_P5P5", "P5P5", Sort_easy);
    c3_s.Read("../charm_E/mix_s_c3", "mes_contr_mix_s_c3_TM_P5P5", "P5P5", Sort_easy);

    distr_t_list corr_s_c1 = Corr.corr_t( c1_s.col(0)[id_E], "");
    distr_t_list corr_s_c2 = Corr.corr_t( c2_s.col(0)[id_E], "");
    distr_t_list corr_s_c3 = Corr.corr_t( c3_s.col(0)[id_E], "");

    Corr.Tmin=49;
    Corr.Tmax=85;


    boost::filesystem::create_directory("../data/charmed_new_ens");
    

    distr_t eff_mass_Ds_1 = Corr.Fit_distr( Corr.effective_mass_t( corr_s_c1, "../data/charmed_new_ens/Ds_1"));
    distr_t eff_mass_Ds_2 = Corr.Fit_distr( Corr.effective_mass_t( corr_s_c2, "../data/charmed_new_ens/Ds_2"));
    distr_t eff_mass_Ds_3 = Corr.Fit_distr( Corr.effective_mass_t( corr_s_c3, "../data/charmed_new_ens/Ds_3"));

    cout<<" Ds1: "<<( eff_mass_Ds_1/aE).ave()<<" +- "<<( eff_mass_Ds_1/aE).err()<<endl;

    cout<<" Ds2: "<<( eff_mass_Ds_2/aE).ave()<<" +- "<<( eff_mass_Ds_3/aE).err()<<endl;

    cout<<" Ds3: "<<( eff_mass_Ds_3/aE).ave()<<" +- "<<( eff_mass_Ds_3/aE).err()<<endl;

    
  

    
    
  distr_t_list corr_s_V_L_tm = Corr.corr_t( V_strange_L.col(0)[id_E], "");
  distr_t_list corr_s_V_L_OS = Corr.corr_t( V_strange_OS_L.col(0)[id_E], "");

  distr_t_list corr_s_V_M_tm = Corr.corr_t( V_strange_M.col(0)[id_E], "");
  distr_t_list corr_s_V_M_OS = Corr.corr_t( V_strange_OS_M.col(0)[id_E], "");

  distr_t_list corr_s_P_L_tm = Corr.corr_t( corr_P5P5_strange.col(0)[id_E], "");
  distr_t_list corr_s_P_L_OS = Corr.corr_t( corr_P5P5_OS_strange.col(0)[id_E], "");

  distr_t_list corr_s_P_M_tm = Corr.corr_t( corr_P5P5_strange_heavy.col(0)[id_E], "");
  distr_t_list corr_s_P_M_OS = Corr.corr_t( corr_P5P5_OS_strange_heavy.col(0)[id_E], "");

  distr_t_list corr_s_AP_L_tm = Corr.corr_t( corr_A0P5_strange.col(0)[id_E], "");
  distr_t_list corr_s_AP_L_OS = Corr.corr_t( corr_A0P5_OS_strange.col(0)[id_E], "");

  distr_t_list corr_s_AP_M_tm = Corr.corr_t( corr_A0P5_strange_heavy.col(0)[id_E], "");
  distr_t_list corr_s_AP_M_OS = Corr.corr_t( corr_A0P5_OS_strange_heavy.col(0)[id_E], "");

  distr_t_list corr_V_tm = (5.0/9.0)*Corr.corr_t( V_light_tm.col(0)[id_E], "");
  distr_t_list corr_V_OS = (5.0/9.0)*Corr.corr_t( V_light_OS.col(0)[id_E], "");


  Corr.Tmin=62;
  Corr.Tmax=100;


   
  LatticeInfo L_info;
     
  L_info.LatInfo_new_ens(V_light_tm.Tag[id_E]);
     
  double aml= L_info.ml;
  double ams1= L_info.ms_L_new;
  double ams2= L_info.ms_M_new;

  distr_t Meta_tm_L = Corr.Fit_distr( Corr.effective_mass_t(corr_s_P_L_tm, "../data/RC_analysis_E/eta_L_tm"));
  distr_t Meta_tm_M = Corr.Fit_distr( Corr.effective_mass_t(corr_s_P_M_tm, "../data/RC_analysis_E/eta_M_tm"));
  distr_t Meta_OS_L = Corr.Fit_distr( Corr.effective_mass_t(corr_s_P_L_OS, "../data/RC_analysis_E/eta_L_OS"));
  distr_t Meta_OS_M = Corr.Fit_distr( Corr.effective_mass_t(corr_s_P_M_OS, "../data/RC_analysis_E/eta_M_OS"));


  
  distr_t_list RV_L = 2*ams1*corr_s_P_L_tm/(distr_t_list::derivative(corr_s_AP_L_tm, 0));  
  distr_t_list RV_M = 2*ams2*corr_s_P_M_tm/(distr_t_list::derivative(corr_s_AP_M_tm, 0));

  distr_t_list RA_L = 2*ams1*corr_s_P_L_OS/(distr_t_list::derivative(corr_s_AP_L_OS, 0)); 
  distr_t_list RA_M = 2*ams2*corr_s_P_M_OS/(distr_t_list::derivative(corr_s_AP_M_OS, 0));

  D(1);


  RA_L = RA_L*(Meta_OS_L/Meta_tm_L)*(SINH_D(Meta_OS_L)/SINH_D(Meta_tm_L))*Corr.matrix_element_t(corr_s_P_L_tm, "")/Corr.matrix_element_t(corr_s_P_L_OS, "");
  RA_M = RA_M*(Meta_OS_M/Meta_tm_M)*(SINH_D(Meta_OS_M)/SINH_D(Meta_tm_M))*Corr.matrix_element_t(corr_s_P_M_tm, "")/Corr.matrix_element_t(corr_s_P_M_OS, "");


  distr_t_list Zp_ov_Zs_L_distr = 1.0/(Corr.matrix_element_t(corr_s_P_L_tm, "")/Corr.matrix_element_t(corr_s_P_L_OS, "") );
  distr_t_list Zp_ov_Zs_M_distr = 1.0/(Corr.matrix_element_t(corr_s_P_M_tm, "")/Corr.matrix_element_t(corr_s_P_M_OS, "") );

  D(2);

  Print_To_File({}, {RV_L.ave(), RV_L.err(), RV_M.ave(), RV_M.err()} , "../data/RC_analysis_E/RV.dat", "", "");
  Print_To_File({}, {RA_L.ave(), RA_L.err(), RA_M.ave(), RA_M.err()} , "../data/RC_analysis_E/RA.dat", "", "");

  Print_To_File({}, {Zp_ov_Zs_L_distr.ave(), Zp_ov_Zs_L_distr.err(), Zp_ov_Zs_M_distr.ave(), Zp_ov_Zs_M_distr.err()} , "../data/RC_analysis_E/Zp_ov_Zs.dat", "", "");

  D(3);

  Corr.Tmin=47;
  Corr.Tmax=100;

  distr_t Zv_L = Corr.Fit_distr(RV_L);
  distr_t Zv_M = Corr.Fit_distr(RV_M);
  distr_t Za_L = Corr.Fit_distr(RA_L);
  distr_t Za_M = Corr.Fit_distr(RA_M);


  distr_t Zp_ov_Zs_L = Corr.Fit_distr(Zp_ov_Zs_L_distr);
  distr_t Zp_ov_Zs_M = Corr.Fit_distr(Zp_ov_Zs_M_distr);


   


  
  vector<distr_t> M2etas_fit({Meta_tm_L*Meta_tm_L/(a_E*a_E), Meta_tm_M*Meta_tm_M/(a_E*a_E)});

  distr_t m_etas_phys_distr;
  const double m_etas = 0.68989;
  const double m_etas_err= 0.00050;
  
 for(int ijack=0;ijack<Njacks;ijack++) m_etas_phys_distr.distr.push_back( m_etas + GM()*m_etas_err/sqrt(Njacks-1.0));
  
  //Generate fake ms_distr
  distr_t ms_light_distr;
  distr_t ms_heavy_distr;
  for(int ijack=0;ijack<Njacks;ijack++) ms_light_distr.distr.push_back( ams1 );
  for(int ijack=0;ijack<Njacks;ijack++) ms_heavy_distr.distr.push_back( ams2 );
  
  //estrapolate ms phys
  vector<distr_t> ms_list( {ms_light_distr, ms_heavy_distr});
  distr_t ms_phys_extr = Obs_extrapolation_meson_mass(ms_list, M2etas_fit, m_etas_phys_distr*m_etas_phys_distr,  "../data/RC_analysis_E", "ms_extrapolation_etas", UseJack, "SPLINE");

  vector<distr_t> Za_hadr_list, Zv_hadr_list, Zp_ov_Zs_hadr_list;
  Za_hadr_list = {Za_L, Za_M};
  Zv_hadr_list = {Zv_L, Zv_M};
  Zp_ov_Zs_hadr_list ={ Zp_ov_Zs_L, Zp_ov_Zs_M};
  

  distr_t Za = Obs_extrapolation_meson_mass( Za_hadr_list, ms_list, ms_phys_extr, "../data/RC_analysis_E", "Za_extr_quark_mass", UseJack, "SPLINE");
  distr_t Zv = Obs_extrapolation_meson_mass( Zv_hadr_list, ms_list, ms_phys_extr, "../data/RC_analysis_E", "Zv_extr_quark_mass", UseJack, "SPLINE");
  distr_t Zp_ov_Zs = Obs_extrapolation_meson_mass( Zp_ov_Zs_hadr_list, ms_list, ms_phys_extr, "../data/RC_analysis_E", "Zp_ov_Zs_extr_quark_mass", UseJack, "SPLINE");

  //evaluate correlator

  int T= Corr.Nt;

  distr_t agm2_SD_tm(UseJack, UseJack?Njacks:1000);
  distr_t agm2_SD_OS(UseJack, UseJack?Njacks:1000);

  double t0=  0.4*fm_to_inv_Gev;
  double Delta= 0.15*fm_to_inv_Gev;

  auto K = [&](double Mv, double t, double size) -> double { return kernel_K(t, Mv);};
  auto th0 = [&t0, &Delta](double ta) ->double { return 1.0/(1.0 + exp(-2.0*(ta-t0)/Delta));};

  distr_t_list Ker = distr_t_list::f_of_distr(K, a_E , T/2);


  //free correlator
  //Generate_free_corr_data();

  
  string Pt_free_oppor= "../Vkvk_cont/"+to_string(Corr.Nt/2)+"_m"+to_string_with_precision(L_info.ml,5)+"/OPPOR";
  Vfloat VV_free_oppor= Read_From_File(Pt_free_oppor, 1, 4);
  if(VV_free_oppor.size() != Corr.Nt) crash("Failed to read properly free VV correlator ml w opposite r");
  //################## READ FREE THEORY VECTOR-VECTOR CORRELATOR SAME R ####################################
  string Pt_free_samer= "../Vkvk_cont/"+to_string(Corr.Nt/2)+"_m"+to_string_with_precision(L_info.ml,5)+"/SAMER";
  Vfloat VV_free_samer= Read_From_File(Pt_free_samer, 1, 4);
  if(VV_free_samer.size() != Corr.Nt) crash("Failed to read properly free VV correlator ml  w same r");
  //Insert electric charges
  for( auto & OP:VV_free_oppor) OP *= 5.0/9.0;
  for( auto & SA:VV_free_samer) SA *= 5.0/9.0;

   corr_V_tm = (1.0/(Za*Za))*(Za*Za*corr_V_tm + VV_free_oppor); //free_corr_log_art
   corr_V_OS = (1.0/(Zv*Zv))*(Zv*Zv*corr_V_OS + VV_free_samer); //free_corr_log_art

 

  for(int t=1; t <Corr.Nt/2; t++) {

    agm2_SD_tm = agm2_SD_tm + 4.0*pow(alpha,2)*w(t,3)*Za*Za*corr_V_tm[t]*Ker[t]*( 1.0 - distr_t::f_of_distr(th0, t*a_E));
    agm2_SD_OS = agm2_SD_OS + 4.0*pow(alpha,2)*w(t,3)*Zv*Zv*corr_V_OS[t]*Ker[t]*( 1.0 - distr_t::f_of_distr(th0, t*a_E));
  }


  cout<<"Za: "<<Za.ave()<<" "<<Za.err()<<endl;
  cout<<"Zv: "<<Zv.ave()<<" "<<Zv.err()<<endl;
  cout<<"Zp/Zs: "<<Zp_ov_Zs.ave()<<" "<<Zp_ov_Zs.err()<<endl;
  
  
  cout<<"tm: "<<agm2_SD_tm.ave()<<" "<<agm2_SD_tm.err()<<endl;
  cout<<"OS: "<<agm2_SD_OS.ave()<<" "<<agm2_SD_OS.err()<<endl;

  cout<<"a: "<<(a_E/fm_to_inv_Gev).ave()<<" "<<(a_E/fm_to_inv_Gev).err()<<endl;


  

















  //#############################################################################################




  int Nens= Vk_data_tm.size;


  for(int iens=0;iens<Nens;iens++) {

    cout<<"Analyzing ensemble: "<<Vk_data_tm.Tag[iens]<<endl; 

    CorrAnalysis Corr(UseJack, Njacks,800);
    Corr.Nt = Vk_data_tm.nrows[iens];
     
    distr_t a_distr(UseJack);
    distr_t Zv(UseJack), Za(UseJack);

    if(Vk_data_tm.Tag[iens].substr(1,1)=="B") {a_distr=a_B; Zv = ZV_B; Za = ZA_B; }
    else if(Vk_data_tm.Tag[iens].substr(1,1)=="C") {a_distr=a_C; Zv = ZV_C; Za = ZA_C;}
    else if(Vk_data_tm.Tag[iens].substr(1,1)=="D") {a_distr=a_D; Zv = ZV_D; Za = ZA_D;}
    else crash("Ensemble not found");

    LatticeInfo L_info;
    L_info.LatInfo_new_ens(Vk_data_tm.Tag[iens]);

    int L= L_info.L;

    distr_t_list Vk_DD_tm_distr(UseJack), Vk_SS_tm_distr(UseJack), Vk_SD_tm_distr(UseJack);
    distr_t_list Vk_DD_OS_distr(UseJack), Vk_SS_OS_distr(UseJack), Vk_SD_OS_distr(UseJack);
    distr_t_list Vk_DEFL_tm_distr(UseJack), Vk_DEFL_OS_distr(UseJack);
      
       
    distr_t_list Vk_tm_distr = 1e10*Za*Za*(pow(qu,2)+pow(qd,2))*Corr.corr_t( Vk_data_tm.col(0)[iens] , "../data/HVP/Corr/Vk_tm_"+Vk_data_tm.Tag[iens]+".dat");
    distr_t_list Vk_OS_distr = 1e10*Zv*Zv*(pow(qu,2)+pow(qd,2))*Corr.corr_t( Vk_data_OS.col(0)[iens] , "../data/HVP/Corr/Vk_OS_"+Vk_data_OS.Tag[iens]+".dat");
    distr_t_list P5_distr= Corr.corr_t( P5P5_data_tm.col(0)[iens] , "");

    if(Vk_data_tm.Tag[iens].substr(1,12) == "B211b.072.64") {

      distr_t R_TM= 1.0/( 1e10*Za*Za*(pow(qu,2)+pow(qd,2)));
      distr_t R_OS= 1.0/( 1e10*Zv*Zv*(pow(qu,2)+pow(qd,2)));
      
      Vk_DD_tm_distr = 1e10*Za*Za*(pow(qu,2)+pow(qd,2))*Corr.corr_t( Vk_data_defl_tm.col(1)[iens] , "../data/HVP/Corr/Vk_DD_tm_"+Vk_data_tm.Tag[iens]+".dat");
      Vk_SD_tm_distr = 2*1e10*Za*Za*(pow(qu,2)+pow(qd,2))*Corr.corr_t( Vk_data_defl_tm.col(3)[iens] , "../data/HVP/Corr/Vk_SD_tm_"+Vk_data_tm.Tag[iens]+".dat");
      Vk_SS_tm_distr = 1e10*Za*Za*(pow(qu,2)+pow(qd,2))*Corr.corr_t( Vk_data_defl_tm.col(5)[iens] , "../data/HVP/Corr/Vk_SS_tm_"+Vk_data_tm.Tag[iens]+".dat");
      Vk_DEFL_tm_distr = Vk_DD_tm_distr+Vk_SD_tm_distr+Vk_SS_tm_distr;
      Print_To_File({},{(R_TM*Vk_DEFL_tm_distr).ave(),(R_TM*Vk_DEFL_tm_distr).err()}, "../data/HVP/Corr/Vk_DEF_tm_"+Vk_data_tm.Tag[iens]+".dat.t", "", ""); 

      Vk_DD_OS_distr = 1e10*Zv*Zv*(pow(qu,2)+pow(qd,2))*Corr.corr_t( Vk_data_defl_OS.col(1)[iens] , "../data/HVP/Corr/Vk_DD_OS_"+Vk_data_tm.Tag[iens]+".dat");
      Vk_SD_OS_distr = 2*1e10*Zv*Zv*(pow(qu,2)+pow(qd,2))*Corr.corr_t( Vk_data_defl_OS.col(3)[iens] , "../data/HVP/Corr/Vk_SD_OS_"+Vk_data_tm.Tag[iens]+".dat");
      Vk_SS_OS_distr = 1e10*Zv*Zv*(pow(qu,2)+pow(qd,2))*Corr.corr_t( Vk_data_defl_OS.col(5)[iens] , "../data/HVP/Corr/Vk_SS_OS_"+Vk_data_tm.Tag[iens]+".dat");
      Vk_DEFL_OS_distr = Vk_DD_OS_distr+Vk_SD_OS_distr+Vk_SS_OS_distr;
      
      Print_To_File({},{(R_OS*Vk_DEFL_OS_distr).ave(),(R_OS*Vk_DEFL_OS_distr).err()}, "../data/HVP/Corr/Vk_DEF_OS_"+Vk_data_tm.Tag[iens]+".dat.t", "", ""); 
      
    }

    if(Vk_data_tm.Tag[iens].substr(1,12)=="B211b.072.96") {Corr.Tmin=30; Corr.Tmax=70;}
    else if(Vk_data_tm.Tag[iens].substr(1,12)=="B211b.072.64") { Corr.Tmin=35; Corr.Tmax=55;}
    else if(Vk_data_tm.Tag[iens].substr(1,1)=="C") {Corr.Tmin=40; Corr.Tmax=60;}
    else if(Vk_data_tm.Tag[iens].substr(1,1)=="D") {Corr.Tmin=41; Corr.Tmax=80;}
    else crash("In scale setting analysis cannot find Tmin,Tmax for ensemble: "+Vk_data_tm.Tag[iens]);
    
    distr_t aMpi_distr= Corr.Fit_distr(Corr.effective_mass_t( P5_distr, ""));

    distr_t p2_rest= 2*aMpi_distr;

    distr_t p2_mot= 2*SQRT_D( aMpi_distr*aMpi_distr + pow( 2*M_PI/L,2)); 

    //B64 110 MeV
    distr_t C= (0.110 - 0.135)/(a_B*a_B);
    distr_t aMpi_N= (0.135 +  C*a_distr*a_distr)*a_distr;
    distr_t p2_mot_tm = 2*SQRT_D( (0.5*(aMpi_distr + aMpi_N))*(0.5*(aMpi_distr + aMpi_N)) + pow(2*M_PI/L,2));

    cout<<"pi neutral: "<< (aMpi_N/a_distr).ave()<< " +- "<< (aMpi_N/a_distr).err()<<endl;

   
    distr_t amu_HVP_tm(UseJack);
    distr_t amu_HVP_OS(UseJack);
    int Tcut_opt_tm, Tcut_opt_OS;

   
          
    Bounding_HVP(amu_HVP_tm, Tcut_opt_tm,  Vk_tm_distr, a_distr,"../data/HVP/Bounding/"+Vk_data_tm.Tag[iens]+"_tm" , p2_mot);
    Bounding_HVP(amu_HVP_OS, Tcut_opt_OS, Vk_OS_distr, a_distr, "../data/HVP/Bounding/"+Vk_data_tm.Tag[iens]+"_OS" , p2_mot);

    distr_t amu_HVP_defl_tm(UseJack), amu_HVP_defl_OS(UseJack);
    distr_t amu_HVP_defl_tm_Blind(UseJack), amu_HVP_defl_OS_Blind(UseJack);
    int Tcut_opt_defl_tm, Tcut_opt_defl_OS;

    Vfloat thetas;
    for(int i=0;i<100;i++) thetas.push_back( 2*M_PI*i/(99.0));
    distr_t_list amu_HVP_theta(UseJack);
    vector<int> Tcut_opt_theta;

    double M1 = 0.750*a_distr.ave(); double M2= 0.850*a_distr.ave(); double M3=0.925*a_distr.ave();
    double M4=1.2*a_distr.ave(), M5 = 1.5*a_distr.ave(), M6 = 1.75*a_distr.ave();

    double M1p = 0.775*a_distr.ave(); double M2p= 0.875*a_distr.ave(); double M3p=0.975*a_distr.ave();
    double M4p=1.4*a_distr.ave(), M5p = 1.7*a_distr.ave(), M6p = 2.2*a_distr.ave();

    auto C1 = [&M1, &a_distr, T=Corr.Nt](double t) { return 0.6*1e10*pow(a_distr.ave(),3)*(1e-3)*(exp(-M1*t) + exp(-M1*(T-t))); };
    auto C2 = [&M2, &a_distr, T=Corr.Nt](double t) { return 1.4e10*pow(a_distr.ave(),3)*(1e-3)*(exp(-M2*t) + exp(-M2*(T-t))); };
    auto C3 = [&M3, &a_distr, T=Corr.Nt](double t) { return 2e10*pow(a_distr.ave(),3)*(1e-3)*(exp(-M3*t) + exp(-M3*(T-t))); };

    auto C4 = [&M4, &a_distr, T=Corr.Nt](double t) { return 1.4e10*pow(a_distr.ave(),3)*(1e-3)*(exp(-M4*t) + exp(-M4*(T-t))); };
    auto C5 = [&M5, &a_distr, T=Corr.Nt](double t) { return 1.5e10*pow(a_distr.ave(),3)*(1e-3)*(exp(-M5*t) + exp(-M5*(T-t))); };
    auto C6 = [&M6, &a_distr, T=Corr.Nt](double t) { return 3e10*pow(a_distr.ave(),3)*(1e-3)*(exp(-M6*t) + exp(-M6*(T-t))); };

    auto C1p = [&M1p, &a_distr, T=Corr.Nt](double t) { return 0*1e10*pow(a_distr.ave(),3)*(1e-3)*(exp(-M1p*t) + exp(-M1p*(T-t))); };
    auto C2p = [&M2p, &a_distr, T=Corr.Nt](double t) { return 0*1e10*pow(a_distr.ave(),3)*(1e-3)*(exp(-M2p*t) + exp(-M2p*(T-t))); };
    auto C3p = [&M3p, &a_distr, T=Corr.Nt](double t) { return 0*1e10*pow(a_distr.ave(),3)*(1e-3)*(exp(-M3p*t) + exp(-M3p*(T-t))); };

    auto C4p = [&M4p, &a_distr, T=Corr.Nt](double t) { return 0*1e10*pow(a_distr.ave(),3)*(1e-3)*(exp(-M4p*t) + exp(-M4p*(T-t))); };
    auto C5p = [&M5p, &a_distr, T=Corr.Nt](double t) { return 0*1e10*pow(a_distr.ave(),3)*(1e-3)*(exp(-M5p*t) + exp(-M5p*(T-t))); };
    auto C6p = [&M6p, &a_distr, T=Corr.Nt](double t) { return 0*1e10*pow(a_distr.ave(),3)*(1e-3)*(exp(-M6p*t) + exp(-M6p*(T-t))); };


    distr_t_list Vk_DEFL_bound_tm= Vk_DEFL_tm_distr;
    distr_t_list Vk_DEFL_bound_OS= Vk_DEFL_OS_distr;

    for(int t=0; t<Vk_DEFL_tm_distr.size();t++ ) {
      Vk_DEFL_bound_tm.distr_list[t] = Vk_DEFL_bound_tm.distr_list[t] + C1(t) + C2(t) + C3(t) + C4(t) + C5(t) + C6(t) +  C1p(t) + C2p(t) + C3p(t) + C4p(t) + C5p(t) + C6p(t);
      Vk_DEFL_bound_OS.distr_list[t] = Vk_DEFL_bound_OS.distr_list[t] + C1(t) + C2(t) + C3(t) + C4(t) + C5(t) + C6(t)  + C1p(t) + C2p(t) + C3p(t) + C4p(t) + C5p(t) + C6p(t);
    }
    
    if(Vk_data_tm.Tag[iens].substr(1,12) == "B211b.072.64") {

      
      Bounding_HVP(amu_HVP_defl_tm, Tcut_opt_defl_tm, Vk_DEFL_tm_distr, a_distr,"../data/HVP/Bounding/"+Vk_data_tm.Tag[iens]+"_DEFL_tm" , p2_mot);
      Bounding_HVP(amu_HVP_defl_OS, Tcut_opt_defl_OS, Vk_DEFL_OS_distr, a_distr, "../data/HVP/Bounding/"+Vk_data_tm.Tag[iens]+"_DEFL_OS" , p2_mot);

      Bounding_HVP(amu_HVP_defl_tm_Blind, Tcut_opt_defl_tm, Vk_DEFL_bound_tm, a_distr,"../data/HVP/Bounding/"+Vk_data_tm.Tag[iens]+"_DEFL_BLIND_tm" , p2_mot);
      Bounding_HVP(amu_HVP_defl_OS_Blind, Tcut_opt_defl_OS, Vk_DEFL_bound_OS, a_distr, "../data/HVP/Bounding/"+Vk_data_tm.Tag[iens]+"_DEFL_BLIND_OS" , p2_mot);

      //plot effective mass
      distr_t_list eff_mass_tm = Corr.effective_mass_t(Vk_DEFL_tm_distr, "../data/HVP/eff_mass/"+Vk_data_tm.Tag[iens]+"_tm");
      distr_t_list eff_mass_OS = Corr.effective_mass_t(Vk_DEFL_OS_distr, "../data/HVP/eff_mass/"+Vk_data_tm.Tag[iens]+"_OS");
      distr_t_list eff_mass_BOUND_tm = Corr.effective_mass_t(Vk_DEFL_bound_tm, "../data/HVP/eff_mass/"+Vk_data_tm.Tag[iens]+"_BLIND_tm");
      distr_t_list eff_mass_BOUND_OS = Corr.effective_mass_t(Vk_DEFL_bound_OS, "../data/HVP/eff_mass/"+Vk_data_tm.Tag[iens]+"_BLIND_OS");


    
      for(int i=0;i<100;i++) {
	int tcut; distr_t amu(UseJack);
	distr_t_list VK =  pow(cos(thetas[i]),2)*Vk_DEFL_tm_distr +  pow(sin(thetas[i]),2)*Vk_DEFL_OS_distr;
	Bounding_HVP(amu, tcut,  VK, a_distr,"../data/HVP/Bounding/"+Vk_data_tm.Tag[iens]+"_TH_"+to_string(thetas[i]) , p2_mot);
	Tcut_opt_theta.push_back(tcut);
	amu_HVP_theta.distr_list.push_back( amu);
	
      }
      
    }


   

   
    vector<string> Tags({"tm", "OS"});
    distr_t_list amu_HVP;
    amu_HVP.distr_list.push_back(amu_HVP_tm);
    amu_HVP.distr_list.push_back(amu_HVP_OS);
    vector<double> Tmins({(double)Tcut_opt_tm, (double)Tcut_opt_OS});
    vector<double> Tmaxs({ Tcut_opt_tm + DTT*fm_to_inv_Gev/a_distr.ave(), Tcut_opt_OS + DTT*fm_to_inv_Gev/a_distr.ave()});

  

    Print_To_File(Tags, { amu_HVP.ave(), amu_HVP.err(), Tmins, Tmaxs}, "../data/HVP/Bounding/"+Vk_data_tm.Tag[iens]+".res", "", "");


    

    
    distr_t_list Ker = distr_t_list::f_of_distr(K, a_distr , Corr.Nt/2);

    distr_t HVP_B_1(UseJack,UseJack?Njacks:800);
    distr_t HVP_B_2(UseJack,UseJack?Njacks:800);
    distr_t HVP_B_3(UseJack,UseJack?Njacks:800);
    distr_t HVP_B_4(UseJack,UseJack?Njacks:800);
    distr_t HVP_B_5(UseJack,UseJack?Njacks:800);
    distr_t HVP_B_6(UseJack,UseJack?Njacks:800);

    distr_t HVP_B_1p(UseJack,UseJack?Njacks:800);
    distr_t HVP_B_2p(UseJack,UseJack?Njacks:800);
    distr_t HVP_B_3p(UseJack,UseJack?Njacks:800);
    distr_t HVP_B_4p(UseJack,UseJack?Njacks:800);
    distr_t HVP_B_5p(UseJack,UseJack?Njacks:800);
    distr_t HVP_B_6p(UseJack,UseJack?Njacks:800);
    
    for(int t=1;t<Corr.Nt/2;t++) HVP_B_1 = HVP_B_1 + 4.0*w(t,2)*pow(alpha,2)*C1(t)*Ker.distr_list[t];
    for(int t=1;t<Corr.Nt/2;t++) HVP_B_2 = HVP_B_2 + 4.0*w(t,2)*pow(alpha,2)*C2(t)*Ker.distr_list[t];
    for(int t=1;t<Corr.Nt/2;t++) HVP_B_3 = HVP_B_3 + 4.0*w(t,2)*pow(alpha,2)*C3(t)*Ker.distr_list[t];

    for(int t=1;t<Corr.Nt/2;t++) HVP_B_4 = HVP_B_4 + 4.0*w(t,2)*pow(alpha,2)*C4(t)*Ker.distr_list[t];
    for(int t=1;t<Corr.Nt/2;t++) HVP_B_5 = HVP_B_5 + 4.0*w(t,2)*pow(alpha,2)*C5(t)*Ker.distr_list[t];
    for(int t=1;t<Corr.Nt/2;t++) HVP_B_6 = HVP_B_6 + 4.0*w(t,2)*pow(alpha,2)*C6(t)*Ker.distr_list[t];


    for(int t=1;t<Corr.Nt/2;t++) HVP_B_1p = HVP_B_1p + 4.0*w(t,2)*pow(alpha,2)*C1p(t)*Ker.distr_list[t];
    for(int t=1;t<Corr.Nt/2;t++) HVP_B_2p = HVP_B_2p + 4.0*w(t,2)*pow(alpha,2)*C2p(t)*Ker.distr_list[t];
    for(int t=1;t<Corr.Nt/2;t++) HVP_B_3p = HVP_B_3p + 4.0*w(t,2)*pow(alpha,2)*C3p(t)*Ker.distr_list[t];

    for(int t=1;t<Corr.Nt/2;t++) HVP_B_4p = HVP_B_4p + 4.0*w(t,2)*pow(alpha,2)*C4p(t)*Ker.distr_list[t];
    for(int t=1;t<Corr.Nt/2;t++) HVP_B_5p = HVP_B_5p + 4.0*w(t,2)*pow(alpha,2)*C5p(t)*Ker.distr_list[t];
    for(int t=1;t<Corr.Nt/2;t++) HVP_B_6p = HVP_B_6p + 4.0*w(t,2)*pow(alpha,2)*C6p(t)*Ker.distr_list[t];


    
    cout<<"########BOUNDING ENS: "<<Vk_data_tm.Tag[iens]<<"########"<<endl;
    cout<<"HVP(1): "<<HVP_B_1.ave()<<" +- "<<HVP_B_1.err()<<endl;
    cout<<"HVP(2): "<<HVP_B_2.ave()<<" +- "<<HVP_B_2.err()<<endl;
    cout<<"HVP(3): "<<HVP_B_3.ave()<<" +- "<<HVP_B_3.err()<<endl;
    cout<<"HVP(4): "<<HVP_B_4.ave()<<" +- "<<HVP_B_4.err()<<endl;
    cout<<"HVP(5): "<<HVP_B_5.ave()<<" +- "<<HVP_B_5.err()<<endl;
    cout<<"HVP(6): "<<HVP_B_6.ave()<<" +- "<<HVP_B_6.err()<<endl;

    cout<<"HVP(1p): "<<HVP_B_1p.ave()<<" +- "<<HVP_B_1p.err()<<endl;
    cout<<"HVP(2p): "<<HVP_B_2p.ave()<<" +- "<<HVP_B_2p.err()<<endl;
    cout<<"HVP(3p): "<<HVP_B_3p.ave()<<" +- "<<HVP_B_3p.err()<<endl;
    cout<<"HVP(4p): "<<HVP_B_4p.ave()<<" +- "<<HVP_B_4p.err()<<endl;
    cout<<"HVP(5p): "<<HVP_B_5p.ave()<<" +- "<<HVP_B_5p.err()<<endl;
    cout<<"HVP(6p): "<<HVP_B_6p.ave()<<" +- "<<HVP_B_6p.err()<<endl;

    if(Vk_data_tm.Tag[iens].substr(1,12) == "B211b.072.64") {
      //fit difference between new and old data to infer blinding parameters

       class ipar_BL {
	
       public:
	 ipar_BL()  {}
	 
	 double dC, dC_err;
	 double t;
       };
      
      
      class fpar_BL {

      public:
	fpar_BL() {}
	fpar_BL(const Vfloat &par) {
	  if((signed)par.size() != 24) crash("In class fpar_BL, class constructor Vfloat par has size != 24");
	  M1=par[0];
	  M2=par[1];
	  M3=par[2];
	  M4=par[3];
	  M5=par[4];
	  M6=par[5];

	  M1p=par[6];
	  M2p=par[7];
	  M3p=par[8];
	  M4p=par[9];
	  M5p=par[10];
	  M6p=par[11];



	  A1=par[12];
	  A2=par[13];
	  A3=par[14];
	  A4=par[15];
	  A5=par[16];
	  A6=par[17];

	  A1p=par[18];
	  A2p=par[19];
	  A3p=par[20];
	  A4p=par[21];
	  A5p=par[22];
	  A6p=par[23];
	}

	double M1, M2, M3, M4, M5, M6 , M1p, M2p, M3p, M4p, M5p, M6p;
	double A1, A2, A3, A4, A5, A6 , A1p, A2p, A3p, A4p, A5p, A6p;
      };



      bootstrap_fit<fpar_BL,ipar_BL> bf_BL(Njacks);
      bootstrap_fit<fpar_BL,ipar_BL> bf_BL_ch2(1);
      bf_BL.Set_number_of_measurements(24);
      bf_BL.Set_verbosity(1);
      //add fit parameters
      bf_BL.Add_par("M1", 0.8, 0.1);
      bf_BL.Add_par("M2", 0.8, 0.1);
      bf_BL.Add_par("M3", 0.8, 0.1);
      bf_BL.Add_par("M4", 1.5, 0.1);
      bf_BL.Add_par("M5", 1.5, 0.1);
      bf_BL.Add_par("M6", 1.5, 0.1);

      bf_BL.Add_par("M1p", 0.8, 0.1);
      bf_BL.Add_par("M2p", 0.8, 0.1);
      bf_BL.Add_par("M3p", 0.8, 0.1);
      bf_BL.Add_par("M4p", 1.5, 0.1);
      bf_BL.Add_par("M5p", 1.5, 0.1);
      bf_BL.Add_par("M6p", 1.5, 0.1);


      bf_BL.Add_par("A1", 0.25, 0.01);
      bf_BL.Add_par("A2", 0.5, 0.01);
      bf_BL.Add_par("A3", 0.5, 0.01);
      bf_BL.Add_par("A4", 0.5, 0.01);
      bf_BL.Add_par("A5", 0.5, 0.01);
      bf_BL.Add_par("A6", 0.5, 0.01);

      bf_BL.Add_par("A1p", 0.25, 0.01);
      bf_BL.Add_par("A2p", 0.5, 0.01);
      bf_BL.Add_par("A3p", 0.5, 0.01);
      bf_BL.Add_par("A4p", 0.5, 0.01);
      bf_BL.Add_par("A5p", 0.5, 0.01);
      bf_BL.Add_par("A6p", 0.5, 0.01);
      
      //ansatz
      bf_BL.ansatz=  [&a_distr , T=Corr.Nt ](const fpar_BL &p, const ipar_BL &ip) {

	Vfloat Ms= {p.M1, p.M2, p.M3, p.M4, p.M5, p.M6, p.M1p, p.M2p, p.M3p, p.M4p, p.M5p, p.M6p};
	Vfloat As= {p.A1, p.A2, p.A3, p.A4, p.A5, p.A6, p.A1p, p.A2p, p.A3p, p.A4p, p.A5p, p.A6p};

	double an=0;
	for (int i=0;i<(signed)Ms.size();i++) an += As[i]*( 1e10*pow(a_distr.ave(),3)*(1e-3)*(exp(-Ms[i]*a_distr.ave()*ip.t) + exp(-Ms[i]*a_distr.ave()*(T-ip.t))));

	return an;
	
      };
      //meas
      bf_BL.measurement=  [ ](const fpar_BL&p, const ipar_BL &ip) {
	return ip.dC;
      };
      //err
      bf_BL.error=  [ ](const fpar_BL &p, const ipar_BL &ip) {
	return ip.dC_err;
      };

      //fill the data
      vector<vector<ipar_BL>> data(Njacks);
      //allocate space for output result
      boot_fit_data<fpar_BL> Bt_fit;
      for(auto &data_iboot: data) data_iboot.resize(24);

      for(int ijack=0;ijack<Njacks;ijack++) {

	for(int t=0; t <24;t++) {

	  data[ijack][t].t = t;
	  data[ijack][t].dC =  0.5*(Vk_DEFL_bound_tm.distr_list[t].distr[ijack] -  Vk_tm_distr.distr_list[t].distr[ijack]); // + Vk_DEFL_bound_OS.distr_list[t].distr[ijack] - Vk_OS_distr.distr_list[t].distr[ijack]);
	  data[ijack][t].dC_err= 0.5*(Vk_DEFL_bound_tm.distr_list[t] -  Vk_tm_distr.distr_list[t]).err() ; // + Vk_DEFL_bound_OS.distr_list[t] - Vk_OS_distr.distr_list[t]).err(); 
	

	}
      }


      bf_BL.Append_to_input_par(data);
      
      //fit
      cout<<"Fitting...."<<endl;
      Bt_fit= bf_BL.Perform_bootstrap_fit();


      //determine fit parameters
      distr_t m1(UseJack), m2(UseJack), m3(UseJack), m4(UseJack), m5(UseJack), m6(UseJack);
      distr_t m1p(UseJack), m2p(UseJack), m3p(UseJack), m4p(UseJack), m5p(UseJack), m6p(UseJack);
      
      distr_t a1(UseJack), a2(UseJack), a3(UseJack), a4(UseJack), a5(UseJack), a6(UseJack);
      distr_t a1p(UseJack), a2p(UseJack), a3p(UseJack), a4p(UseJack), a5p(UseJack), a6p(UseJack);

     
      
      for(int ijack=0;ijack<Njacks;ijack++) {
	
	m1.distr.push_back( Bt_fit.par[ijack].M1);
	m2.distr.push_back( Bt_fit.par[ijack].M2);
	m3.distr.push_back( Bt_fit.par[ijack].M3);
	m4.distr.push_back( Bt_fit.par[ijack].M4);
	m5.distr.push_back( Bt_fit.par[ijack].M5);
	m6.distr.push_back( Bt_fit.par[ijack].M6);

	m1p.distr.push_back( Bt_fit.par[ijack].M1p);
	m2p.distr.push_back( Bt_fit.par[ijack].M2p);
	m3p.distr.push_back( Bt_fit.par[ijack].M3p);
	m4p.distr.push_back( Bt_fit.par[ijack].M4p);
	m5p.distr.push_back( Bt_fit.par[ijack].M5p);
	m6p.distr.push_back( Bt_fit.par[ijack].M6p);


	a1.distr.push_back( Bt_fit.par[ijack].A1);
	a2.distr.push_back( Bt_fit.par[ijack].A2);
	a3.distr.push_back( Bt_fit.par[ijack].A3);
	a4.distr.push_back( Bt_fit.par[ijack].A4);
	a5.distr.push_back( Bt_fit.par[ijack].A5);
	a6.distr.push_back( Bt_fit.par[ijack].A6);

	a1p.distr.push_back( Bt_fit.par[ijack].A1p);
	a2p.distr.push_back( Bt_fit.par[ijack].A2p);
	a3p.distr.push_back( Bt_fit.par[ijack].A3p);
	a4p.distr.push_back( Bt_fit.par[ijack].A4p);
	a5p.distr.push_back( Bt_fit.par[ijack].A5p);
	a6p.distr.push_back( Bt_fit.par[ijack].A6p);

	
      }


      distr_t_list Corr_Blind_RECO(UseJack);

      for(int t=0;t<Corr.Nt/2;t++) {
	ipar_BL iBL;
	iBL.t = t;
	distr_t A(UseJack);
	for(int ijack=0;ijack<Njacks;ijack++) {
	  A.distr.push_back( bf_BL.ansatz( Bt_fit.par[ijack], iBL));
	  
	}
	Corr_Blind_RECO.distr_list.push_back( A);
      }
      
      distr_t HVP_Blind_RECO(UseJack,UseJack?Njacks:800);
      

      for(int t=1;t<Corr.Nt/2;t++) HVP_Blind_RECO = HVP_Blind_RECO + 4.0*w(t,2)*pow(alpha,2)*Corr_Blind_RECO.distr_list[t]*Ker.distr_list[t];
      

      cout<<"HVP FROM BLINDING: "<<HVP_Blind_RECO.ave()<<" +- "<<HVP_Blind_RECO.err()<<endl;
      
      
     
    }

    

    if(Vk_data_tm.Tag[iens].substr(1,12) == "B211b.072.64") {
      distr_t_list amu_defl_HVP;
      amu_defl_HVP.distr_list.push_back(amu_HVP_defl_tm);
      amu_defl_HVP.distr_list.push_back(amu_HVP_defl_OS);
      vector<double> Tmins_defl({(double)Tcut_opt_defl_tm, (double)Tcut_opt_defl_OS});
      vector<double> Tmaxs_defl({ Tcut_opt_defl_tm+ DTT*fm_to_inv_Gev/a_distr.ave(), Tcut_opt_defl_OS + DTT*fm_to_inv_Gev/a_distr.ave()});
      Print_To_File(Tags, { amu_defl_HVP.ave(), amu_defl_HVP.err(), Tmins_defl, Tmaxs_defl}, "../data/HVP/Bounding/"+Vk_data_tm.Tag[iens]+"_defl.res", "", "");
    }

    cout<<"#### "<<Vk_data_tm.Tag[iens]<<" ###"<<endl;
    cout<<"HVP tm: "<<amu_HVP_tm.ave()<<" +- "<<amu_HVP_tm.err()<<" stat. "<< (amu_HVP_tm.err()*100/amu_HVP_tm.ave())<<"%"<<endl;
    cout<<"HVP OS: "<<amu_HVP_OS.ave()<<" +- "<<amu_HVP_OS.err()<<" stat. "<< (amu_HVP_OS.err()*100/amu_HVP_OS.ave())<<"%"<<endl;
    if(Vk_data_tm.Tag[iens].substr(1,12) == "B211b.072.64") {

      distr_t HVP_Blind_tot= HVP_B_1+HVP_B_2+HVP_B_3+HVP_B_4+HVP_B_5+HVP_B_6 + HVP_B_1p+HVP_B_2p+HVP_B_3p+HVP_B_4p+HVP_B_5p+HVP_B_6p;
      cout<<"HVP(defl) tm: "<<(amu_HVP_defl_tm_Blind -HVP_Blind_tot).ave()<<" +- "<<(amu_HVP_defl_tm_Blind - HVP_Blind_tot).err()<<" stat. "<< ((amu_HVP_defl_tm_Blind -HVP_Blind_tot).err()*100/(amu_HVP_defl_tm_Blind - HVP_Blind_tot).ave())<<"%"<<endl;
      cout<<"HVP(defl) OS: "<<(amu_HVP_defl_OS_Blind - HVP_Blind_tot).ave()<<" +- "<<(amu_HVP_defl_OS_Blind - HVP_Blind_tot).err()<<" stat. "<< ((amu_HVP_defl_OS_Blind - HVP_Blind_tot).err()*100/(amu_HVP_defl_OS_Blind - HVP_Blind_tot).ave())<<"%"<<endl;
      cout<<"###### FROM BLIND ANALYSIS ##########: "<<endl;
      cout<<"HVP(defl,BLIND) tm: "<<amu_HVP_defl_tm.ave()<<" +- "<<amu_HVP_defl_tm.err()<<" stat. "<< (amu_HVP_defl_tm.err()*100/amu_HVP_defl_tm.ave())<<"%"<<endl;
      cout<<"HVP(defl,BLIND) OS: "<<amu_HVP_defl_OS.ave()<<" +- "<<amu_HVP_defl_OS.err()<<" stat. "<< (amu_HVP_defl_OS.err()*100/amu_HVP_defl_OS.ave())<<"%"<<endl;
      cout<<"#######THETA#######"<<endl;
      //for(int i=0;i<(signed)thetas.size();i++) {
       //cout<<"HVP(defl, th: "<<thetas[i]<<" ) :"<<amu_HVP_theta.ave(i)<<" +- "<<amu_HVP_theta.err(i)<<" stat. "<< (amu_HVP_theta.err(i)*100/amu_HVP_theta.ave(i))<<"%"<<endl;
       //}
    }
    cout<<"#######"<<endl;

    auto K = [&](double Mv, double t, double size) -> double { return kernel_K(t, Mv);};

   
    
    
    //distr_t_list Ker_el = 4*pow(alpha,2)*distr_t_list::f_of_distr(K, a_distr*(0.0005/0.10565837) , Corr.Nt);
    //distr_t_list Ker_mu = 4*pow(alpha,2)*distr_t_list::f_of_distr(K, a_distr*(1.77/0.1056837) , Corr.Nt);



    

   
       
  }


  return;

}
