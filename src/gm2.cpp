#include "../include/gm2.h"



using namespace std;

//set constants
namespace plt = matplotlibcpp;

constexpr double kappa=2.837297;
const double MPiPhys=0.135;
const double m_pi_charged = 0.139;
const double alpha = 1.0/137.035999;
const double e2 = alpha*4.0*M_PI;
const int nboots= 200;
const bool UseJack=1;
const int Njacks=50; //50
const int Nboots=200;
const double qu= 2.0/3.0;
const double qd= -1.0/3.0;
const double qs= qd;
const double qc= qu;
const int Upper_Limit_Time_Integral_strange= 300;
const int Upper_Limit_Time_Integral_charm= 300;
const int Upper_Limit_Time_Integral_light=300;
const double fm_to_inv_Gev= 1.0/0.197327;
const bool verbosity=1;
const double Nresonances= 2; //15; //50; //50; //25; //2; //normally you used 12
const int Luscher_num_zeroes= 4; //17; //55; // 55; //30; // 4; //normally you used  20
const int npts_spline= 10; //200; //1000; //10; //normally you used 1000 
bool Use_Mpi_OS=false;
bool Include_light_disco= true;
bool Include_strange_disco= true;
bool Include_charm_disco=true;
bool Include_off_diagonal_disco=true;
double Qfact= 10.0/9.0;
double grpp_phys = 5.5;
double grpp_phys_BMW = 5.95;
const double m_Jpsi= 3.0969;
const double m_Jpsi_err= 0.001;
const double m_phi= 1.019461;
const double m_phi_err= 0.000016;
const double m_pi = 0.135;
const double m_rho= 0.7754;
const double m_k= 0.497611;
const double m_d= 1.86484;
const double m_etac= 2.9839;
const double m_etac_err = 0.004;
const double m_etas = 0.68989;
const double m_etas_err= 0.00050;
const double fp_phys= 0.1304;
const double Gamma_rho=0.1491;
const double X_pi_phys = pow( pow(m_pi,4)*fp_phys,1.0/5.0);
const double fk_phys= 0.1557; //FLAG REV ("2021")
const double fd_phys= 0.212; //FLAG REV ("2021")
const double csi_phys= pow(m_pi/(4.0*M_PI*fp_phys),2);
const double t0 = 0.4*fm_to_inv_Gev;
const double t1 = 1.0*fm_to_inv_Gev;
const double Delta= 0.15*fm_to_inv_Gev;
const double Delta_prime= 0.30*fm_to_inv_Gev;
const double l1ph= -0.4; //-0.4
const double l2ph= 4.3; //4.3
const double l3ph= 3.2; //3.2
const double l4ph= 4.4; //4.4
const double s0= 2.0-M_PI/2.0;
const double s1 = M_PI/4.0 - 0.5;
const double s2 = 0.5 - M_PI/8.0;
const double s3 = 3.0*M_PI/16.0 - 0.5;
const double Csi_K_phys= m_k*m_k/(fk_phys*fk_phys) ;
const double Csi_D_phys= m_d*m_d/(fd_phys*fd_phys) ;
bool Skip_total_light_calc= false;
bool Fit_only_light_corr_tail=true;
bool Reco_agm2_total=false;
string Mod="FLEE_THEORY";
bool Use_Za_Zv_from_strange_run=true;
bool Use_Za_Zv_from_charm_run = false;
bool Use_Za_Zv_from_strange_run_in_light_corr=true;
bool Use_Extrapolated_Za_Zv_strange = true;
bool Use_Extrapolated_Za_Zv_charm= true;
string Extrapolation_strange_mode="etas";
string ELM_mass_strange = "phi";
string Extrapolation_charm_mode= "etac";
string ELM_mass_charm = "Jpsi";
const int pert_corr_charm_on_off= 1;
const int pert_corr_strange_on_off=1;
const int pert_corr_light_on_off=1;
const int sum_pert_corr_charm_to_bare_corr= 0;
const int sum_pert_corr_strange_to_bare_corr=0;
const bool Gen_free_corr_data=false;
double add_pert_corr_charm_up_to = 2.0;
double add_pert_corr_strange_up_to = 2.0;
double add_pert_corr_light_up_to = 3.0;
int Simps_ord= 1 ;
bool Determine_scale_setting=true; //true
bool Use_three_finest_in_scale_setting_w0X=false;
bool Use_three_finest_in_scale_setting_fp=false;
bool scale_setting_from_w0X=false;
bool scale_setting_from_fp = true; //true 
bool Use_scale_setting_from_this_analysis=true; //true
bool Print_FSEs_from_GSLL_param=false;
bool Print_Mpi_dep_from_GSLL_param=false;
//Vfloat tmins({0.08});
Vfloat tmins({ 0.08, 0.09, 0.100, 0.11, 0.125, 0.130, 0.140, 0.150}); //tmins for SD extrapolation in fermi
Vfloat Qs2;
const int eps_win_size=12;
const Vfloat eps_win({ 0.0*fm_to_inv_Gev, -0.05*fm_to_inv_Gev, -0.10*fm_to_inv_Gev, -0.20*fm_to_inv_Gev, 0.05*fm_to_inv_Gev, 0.1*fm_to_inv_Gev, 0.15*fm_to_inv_Gev, 0.2*fm_to_inv_Gev, 0.25*fm_to_inv_Gev, 0.3*fm_to_inv_Gev, 0.35*fm_to_inv_Gev, 0.4*fm_to_inv_Gev});


class ipar {

public:
  ipar() : V_light(0.0), V_light_err(0.0) {}

  
  double Mp, Mp_OS;
  double t,  L;
  double V_light, V_light_err;
  double fp;
   


};

class fit_par {

public:
  fit_par() {}
  fit_par(const Vfloat &par) {
    if((signed)par.size() != 6) crash("In class fit_par in fitting analytic representation of V(t)_light, class constructor Vfloat par has size != 6");
    Rd=par[0];
    Ed=par[1];
    Mrho=par[2];
    gpi=par[3];
    kappa= par[4];
    Mpi_corr=par[5];
    
  }

  double Rd,Ed, Mrho, gpi,kappa, Mpi_corr;



};

void Init_Qs2() {


  bool mode_skip=1;
  
  double q2_1=0.01;
  double q2_2=0.1;
  double q2_3=0.5;
  double q2_4=2.0;
  double q2_5=5.0;
  double q2_6=7.0;
  double q2_7=14.0;

  Qs2.push_back( q2_1);
  if(mode_skip) return;
  while( q2_1 <= q2_2+1e-6) { q2_1+= 0.01; Qs2.push_back( q2_1);};
  while( q2_2 <= q2_3+1e-6) { q2_2+= 0.02; Qs2.push_back( q2_2);};
  while( q2_3 <= q2_4+1e-6) { q2_3+= 0.05; Qs2.push_back( q2_3);};
  while( q2_4 <= q2_5+1e-6) { q2_4+= 0.10; Qs2.push_back( q2_4);};
  while( q2_5 <= q2_6+1e-6) { q2_5+= 0.20; Qs2.push_back( q2_5);};
  while( q2_6 <= q2_7+1e-6) { q2_6+= 0.50; Qs2.push_back( q2_6);};

 

  return;

}


distr_t PI_q2( const distr_t_list &V,const distr_t &a, double Q2, int Tmax) {

  double T= V.size();
  distr_t res= 0.0*a;

  distr_t_list Ker_PI(UseJack);

  for(int t=0; t < T; t++) {

    distr_t Ker_PI_t(UseJack);
    for(int ijack=0;ijack<Njacks;ijack++) { Ker_PI_t.distr.push_back( Kernel_Pi_q2( 1.0*t, sqrt( Q2), a.distr[ijack]));}

    Ker_PI.distr_list.push_back( Ker_PI_t);
  }
    
  
  for(int t=1;t < Tmax; t++) {

    res = res + Ker_PI.distr_list[t]*V.distr_list[t];

  }

  return res;


}

distr_t PI_q2( const distr_t_list &V,const distr_t &a, const distr_t &Q2, int Tmax) {

  double T= V.size();
  distr_t res= 0.0*a;

  distr_t_list Ker_PI(UseJack);

  for(int t=0; t < T; t++) {

    distr_t Ker_PI_t(UseJack);
    for(int ijack=0;ijack<Njacks;ijack++) { Ker_PI_t.distr.push_back( Kernel_Pi_q2( 1.0*t, sqrt( Q2.distr[ijack]), a.distr[ijack]));}

    Ker_PI.distr_list.push_back( Ker_PI_t);
  }
    
  
  for(int t=1;t < Tmax; t++) {

    res = res + Ker_PI.distr_list[t]*V.distr_list[t];

  }

  return res;


}

distr_t PI_q2_fixed_t(const distr_t &Vt, const distr_t &a, const distr_t &Q2, int t) {

 
  distr_t Ker_PI_t(UseJack);
  for(int ijack=0;ijack<Njacks;ijack++) { Ker_PI_t.distr.push_back( Kernel_Pi_q2( 1.0*t, sqrt( Q2.distr[ijack]), a.distr[ijack]));}

  return Ker_PI_t*Vt;


}

void Add_ens_val_PI_q2( vector<distr_t_list> &PI_master, const distr_t_list &V, const distr_t &a, int Tmax) {

  //sanity check
  if(PI_master.size() != Qs2.size()) crash("In add_ens_val_PI_q2, size of PI_master and Qs are not equal");

  int Q_size= Qs2.size();

  for(int q=0; q < Q_size; q++) {

    PI_master[q].distr_list.push_back( PI_q2(V,a,Qs2[q]*a*a, Tmax));

  }

  return;

}

void Add_ens_val_PI_q2( vector<distr_t_list> &PI_master, const distr_t_list &PI_per_ens) {

  //sanity check
  if(PI_master.size() != Qs2.size()) crash("In add_ens_val_PI_q2, size of PI_master and Qs are not equal");
  if(PI_per_ens.size() != (signed)Qs2.size()) crash("In add_ens_val_PI_q2, size of PI_per_ens and Qs are not equal");

  int Q_size= Qs2.size();

  for(int q=0; q < Q_size; q++) {

    PI_master[q].distr_list.push_back( PI_per_ens.distr_list[q]);

  }

  return;

}

void Get_PI_q2( distr_t_list &PI_per_ens, const distr_t_list &V, const distr_t &a, int Tmax) {

  int Q_size= Qs2.size();

  for(int q=0; q<Q_size;q++) {

    PI_per_ens.distr_list.push_back( PI_q2(V,a, Qs2[q]*a*a, Tmax));

  }

  return;

}

void Get_PI_q2(distr_t_list &PI_per_ens, const Vfloat &V, const distr_t &a, int Tmax) {

  GaussianMersenne G(435435);
  int size= (signed)V.size();
  distr_t_list V_distr(UseJack, size);
  for(int t=0; t<size;t++) {

    for(int ijack=0;ijack<Njacks;ijack++) V_distr.distr_list[t].distr.push_back( V[t] + G()*1e-16*V[t]/sqrt(Njacks-1));
    

  }

  Get_PI_q2(PI_per_ens, V_distr, a, Tmax);
  return;

}

void Bounding_amu_W(distr_t &amu_W, const distr_t_list &V, const distr_t &a, string path,const distr_t &Z, distr_t lowest_mass) {


  double T= V.size(); 

  auto LOG = [](double R_G, double t) { return log(fabs(R_G));};

  auto K = [&](double Mv, double t, double size) -> double { return kernel_K(t, Mv);};

  distr_t_list Ker = distr_t_list::f_of_distr(K, a , T/2);

  int T_ext_max= 300;

  distr_t_list Ker_extended= distr_t_list::f_of_distr(K, a , T_ext_max);

  int Tdatas_opt=-1;

 //define lambdas for the theta func
    auto th0 = [](double ta) ->double { return 1.0/(1.0 + exp(-2.0*(ta-t0)/Delta));};
    auto th1 = [](double ta) ->double { return 1.0/(1.0 + exp(-2.0*(ta-t1)/Delta));};

    auto amu_W_func = [&V, &a,&th0, &th1, &Ker, &Z](double tcut) -> distr_t {

			distr_t ret_win(UseJack,UseJack?Njacks:Nboots);
			for(int t=1;t<tcut;t++) ret_win = ret_win + 4.0*w(t,Simps_ord)*pow(alpha,2)*Z*Z*V.distr_list[t]*Ker.distr_list[t]*( distr_t::f_of_distr(th0, t*a) - distr_t::f_of_distr(th1, t*a));
			return ret_win;

		      };

    auto amu_W_int = [&a, &th0, &th1, &Ker_extended, &Z](double t) -> distr_t {



		       return 4.0*w(t,Simps_ord)*pow(alpha,2)*Z*Z*Ker_extended.distr_list[t]*( distr_t::f_of_distr(th0, t*a) - distr_t::f_of_distr(th1, t*a));


		     };
  
 
  
  CorrAnalysis Corr(UseJack, Njacks,Nboots);
  Corr.Nt = V.size();
  distr_t_list ratio_corr_V(UseJack);
  for(int t=0; t<Corr.Nt;t++) ratio_corr_V.distr_list.push_back( V.distr_list[t]/V.distr_list[(t+1)%Corr.Nt]);
  //decide whether to use m_eff(t) or log( V(t)/V(t+1)) (true = use log)  
  bool update_min_Tdata_from_log_ratio=true;    
  distr_t_list eff_mass_V = (update_min_Tdata_from_log_ratio)?distr_t_list::f_of_distr_list(LOG, ratio_corr_V):Corr.effective_mass_t(V, "");

  

  distr_t_list amu_W_min_Tdata(UseJack), amu_W_max_Tdata(UseJack), amu_W_T_2(UseJack);

    Vfloat TCUTS;

    Vfloat Is_T_data_opt;
    int slice_to_use_for_eff_mass=1;
    
    
    //loop over tcut
    for(int tcut=1; tcut<T/2;tcut++) {


      TCUTS.push_back( (double)tcut);
      distr_t amu_W_up_to_tcut = amu_W_func(tcut+1);
      amu_W_min_Tdata.distr_list.push_back(amu_W_up_to_tcut);
      amu_W_max_Tdata.distr_list.push_back(amu_W_up_to_tcut);
      amu_W_T_2.distr_list.push_back(amu_W_up_to_tcut);
      distr_t V_tcut = V.distr_list[tcut];
      Is_T_data_opt.push_back( 0.0);
      
      bool eff_mass_is_nan= isnan( eff_mass_V.ave(tcut));
      bool update_min_Tdata=true;
      if(eff_mass_is_nan || (eff_mass_V.err(tcut)/eff_mass_V.ave(tcut) > 0.05) || (eff_mass_V.ave(tcut) < 0 )) update_min_Tdata=false;

      if(update_min_Tdata) slice_to_use_for_eff_mass=tcut;
      
      

      for(int t=tcut+1; t < T_ext_max;t++) {

	//lambda function for lower and upper limit of single exp V(t)
	auto EXP_MIN = [&tcut, &t](double E) { return exp(-E*(t-tcut));};
	
	distr_t ker_val = V_tcut*amu_W_int(t);
	int size_min= amu_W_min_Tdata.size();
	int size_max= amu_W_max_Tdata.size();
	if(size_min != tcut || size_max != tcut) crash("size_min or size_max is different from tcut");
	distr_t lower_exp = distr_t::f_of_distr(EXP_MIN,eff_mass_V.distr_list[slice_to_use_for_eff_mass]);
	amu_W_min_Tdata.distr_list[size_min-1] = amu_W_min_Tdata.distr_list[size_min-1] + ker_val*lower_exp;
	amu_W_max_Tdata.distr_list[size_max-1] = amu_W_max_Tdata.distr_list[size_max-1] + ker_val*distr_t::f_of_distr(EXP_MIN, lowest_mass);
      }

      
    }


   


      //find interval where difference between amu_W_min_Tdata and amu_W_max_Tdata is smaller than 0.3 sigma
      //average
      bool Found_Tdata_opt=false;
      int tdata_opt=1;

   
      while(!Found_Tdata_opt && tdata_opt < T/2) {

	distr_t diff_max_min = amu_W_max_Tdata.distr_list[tdata_opt-1] - amu_W_min_Tdata.distr_list[tdata_opt-1];
	if( diff_max_min.ave()/min( amu_W_max_Tdata.err(tdata_opt-1) , amu_W_min_Tdata.err(tdata_opt-1)) < 0.3) Found_Tdata_opt=true;
	else tdata_opt++;
      }

      

      //if tdata_opt has not been found return tdata = -1 and amu_W =amu_W(T/2)
      if(!Found_Tdata_opt) {
	
	Tdatas_opt = -1;
	amu_W =  amu_W_func(T/2);
	
      }
      else { //tdata_opt has been found


	bool fit_to_constant=true;
	//average over min max over an interval of 0.5 fm
	if( amu_W_max_Tdata.size() != amu_W_min_Tdata.size()) crash("Size of amu_W_max_Tdata and amu_W_min_Tdata are not equal");
	int S_Tdata= amu_W_min_Tdata.size();
	int DT = min( (int)(0.5*fm_to_inv_Gev/a.ave()), (S_Tdata -tdata_opt));
	distr_t amu_W_ave_max= amu_W_max_Tdata.distr_list[tdata_opt-1];
	distr_t amu_W_ave_min= amu_W_min_Tdata.distr_list[tdata_opt-1];
	if(fit_to_constant) {
	  amu_W_ave_max = amu_W_ave_max*1.0/pow(amu_W_max_Tdata.err(tdata_opt-1),2);
	  amu_W_ave_min = amu_W_ave_min*1.0/pow(amu_W_min_Tdata.err(tdata_opt-1),2);
	}
	double amu_W_weight_max = (fit_to_constant)?1.0/(pow(amu_W_max_Tdata.err(tdata_opt-1),2)):1.0;
	double amu_W_weight_min = (fit_to_constant)?1.0/(pow(amu_W_min_Tdata.err(tdata_opt-1),2)):1.0;
	for(int t=tdata_opt; t< tdata_opt+DT; t++) {

	  double weight_max_t = (fit_to_constant)?1.0/pow(amu_W_max_Tdata.err(t),2):1.0;
	  double weight_min_t = (fit_to_constant)?1.0/pow(amu_W_min_Tdata.err(t),2):1.0;
	 
	  amu_W_ave_max = amu_W_ave_max + amu_W_max_Tdata.distr_list[t]*weight_max_t;
	  amu_W_ave_min = amu_W_ave_min + amu_W_min_Tdata.distr_list[t]*weight_min_t;
	  amu_W_weight_max += weight_max_t;
	  amu_W_weight_min += weight_min_t;

	}

	amu_W_ave_max= amu_W_ave_max/amu_W_weight_max;
	amu_W_ave_min= amu_W_ave_min/amu_W_weight_min;


	double den_weight= 1.0/(pow(amu_W_ave_max.err(),2))  + 1.0/(pow(amu_W_ave_min.err(),2));
	double fin_weight_max= (fit_to_constant)?(1.0/(pow(amu_W_ave_max.err(),2)*den_weight)):0.5;
	double fin_weight_min= (fit_to_constant)?(1.0/(pow(amu_W_ave_min.err(),2)*den_weight)):0.5;

	Tdatas_opt =tdata_opt;
	amu_W = fin_weight_min*amu_W_ave_min + fin_weight_max*amu_W_ave_max;
	Is_T_data_opt[tdata_opt-1] = 1.0;

      }


      //Print to File
      Print_To_File({}, {TCUTS, amu_W_min_Tdata.ave(), amu_W_min_Tdata.err(), amu_W_max_Tdata.ave(), amu_W_max_Tdata.err(), amu_W_T_2.ave(), amu_W_T_2.err(), Is_T_data_opt}, path+".dat", "", "#tcut   lower    upper     T/2      Is_Tcut_opt.       Tcut_opt= "+to_string(Tdatas_opt));



  return;

}




void Bounding_PI_q2(distr_t_list &PI_per_ens, const distr_t_list &V, const distr_t &a, string path, Vint &Tdatas_opt, distr_t &lowest_mass) {


  auto LOG = [](double R_G, double t) { return log(fabs(R_G));};
  

 
    
  int T_ext_max= 300;
 
  
  CorrAnalysis Corr(UseJack, Njacks,Nboots);
  Corr.Nt = V.size();
  distr_t_list ratio_corr_V(UseJack);
  for(int t=0; t<Corr.Nt;t++) ratio_corr_V.distr_list.push_back( V.distr_list[t]/V.distr_list[(t+1)%Corr.Nt]);
  //decide whether to use m_eff(t) or log( V(t)/V(t+1)) (true = use log)  
  bool update_min_Tdata_from_log_ratio=true;    
  distr_t_list eff_mass_V = (update_min_Tdata_from_log_ratio)?distr_t_list::f_of_distr_list(LOG, ratio_corr_V):Corr.effective_mass_t(V, "");
  

  //loop over Qs2
  
  for(int q=0; q < (signed)Qs2.size(); q++) {

   


    distr_t_list PI_q2_min_Tdata(UseJack), PI_q2_max_Tdata(UseJack), PI_q2_T_2(UseJack);

    Vfloat TCUTS;

    Vfloat Is_T_data_opt;
    int slice_to_use_for_eff_mass=1;
    
    
    //loop over tcut
    for(int tcut=1; tcut<V.size()/2 -2 ;tcut++) {


      TCUTS.push_back( (double)tcut);
      distr_t PI_q2_up_to_tcut = PI_q2(V, a, Qs2[q]*a*a, tcut+1);
      PI_q2_min_Tdata.distr_list.push_back(PI_q2_up_to_tcut);
      PI_q2_max_Tdata.distr_list.push_back(PI_q2_up_to_tcut);
      PI_q2_T_2.distr_list.push_back(PI_q2_up_to_tcut);
      distr_t V_tcut = V.distr_list[tcut];
      Is_T_data_opt.push_back( 0.0);
      
      bool eff_mass_is_nan= isnan( eff_mass_V.ave(tcut));
      bool update_min_Tdata=true;
      if(eff_mass_is_nan || (eff_mass_V.err(tcut)/eff_mass_V.ave(tcut) > 0.05) || (eff_mass_V.ave(tcut) < 0 )) update_min_Tdata=false;

      if(update_min_Tdata) slice_to_use_for_eff_mass=tcut;
      
      

      for(int t=tcut+1; t < T_ext_max;t++) {

	//lambda function for lower and upper limit of single exp V(t)
	auto EXP_MIN = [&tcut, &t](double E) { return exp(-E*(t-tcut));};
	
	distr_t ker_val = PI_q2_fixed_t(V_tcut, a, Qs2[q]*a*a , t);
	int size_min= PI_q2_min_Tdata.size();
	int size_max= PI_q2_max_Tdata.size();
	if(size_min != tcut || size_max != tcut) crash("size_min or size_max is different from tcut");
	distr_t lower_exp = distr_t::f_of_distr(EXP_MIN,eff_mass_V.distr_list[slice_to_use_for_eff_mass]);
	PI_q2_min_Tdata.distr_list[size_min-1] = PI_q2_min_Tdata.distr_list[size_min-1] + ker_val*lower_exp;
	PI_q2_max_Tdata.distr_list[size_max-1] = PI_q2_max_Tdata.distr_list[size_max-1] + ker_val*distr_t::f_of_distr(EXP_MIN, lowest_mass);
      }

      
    }


   


      //find interval where difference between PI_q2_min_Tdata and PI_q2_max_Tdata is smaller than 0.5 sigma
      //average
      bool Found_Tdata_opt=false;
      int tdata_opt=1;
      int NT= V.size();
      while(!Found_Tdata_opt && tdata_opt < NT/2 -2 ) {

	distr_t diff_max_min = PI_q2_max_Tdata.distr_list[tdata_opt-1] - PI_q2_min_Tdata.distr_list[tdata_opt-1];
	if( diff_max_min.ave()/min( PI_q2_max_Tdata.err(tdata_opt-1) , PI_q2_min_Tdata.err(tdata_opt-1)) < 0.5) Found_Tdata_opt=true;
	else tdata_opt++;
      }

      

      //if tdata_opt has not been found return tdata = -1 and PI(q^2) = PI(q^2, T/2)
      if(!Found_Tdata_opt) {
	
	Tdatas_opt.push_back( -1);
	PI_per_ens.distr_list.push_back(  PI_q2(V,a,Qs2[q]*a*a, V.size()/2));
	
      }
      else { //tdata_opt has been found

	bool fit_to_constant=true;
	//average over min max over an interval of 0.5 fm
	if( PI_q2_max_Tdata.size() != PI_q2_min_Tdata.size()) crash("Size of PI_q2_max_Tdata and PI_q2_min_Tdata are not equal");
	int S_Tdata= PI_q2_min_Tdata.size();
	int DT = min( (int)(0.5*fm_to_inv_Gev/a.ave()), (S_Tdata -tdata_opt));
	distr_t PI_ave_max= PI_q2_max_Tdata.distr_list[tdata_opt-1];
	distr_t PI_ave_min= PI_q2_min_Tdata.distr_list[tdata_opt-1];
	if(fit_to_constant) {
	  PI_ave_max = PI_ave_max*1.0/pow(PI_q2_max_Tdata.err(tdata_opt-1),2);
	  PI_ave_min = PI_ave_min*1.0/pow(PI_q2_min_Tdata.err(tdata_opt-1),2);
	}
	double PI_weight_max = (fit_to_constant)?1.0/(pow(PI_q2_max_Tdata.err(tdata_opt-1),2)):1.0;
	double PI_weight_min = (fit_to_constant)?1.0/(pow(PI_q2_min_Tdata.err(tdata_opt-1),2)):1.0;
	for(int t=tdata_opt; t< tdata_opt+DT; t++) {

	  double weight_max_t = (fit_to_constant)?1.0/pow(PI_q2_max_Tdata.err(t),2):1.0;
	  double weight_min_t = (fit_to_constant)?1.0/pow(PI_q2_min_Tdata.err(t),2):1.0;
	 
	  PI_ave_max = PI_ave_max + PI_q2_max_Tdata.distr_list[t]*weight_max_t;
	  PI_ave_min = PI_ave_min + PI_q2_min_Tdata.distr_list[t]*weight_min_t;
	  PI_weight_max += weight_max_t;
	  PI_weight_min += weight_min_t;

	}

	PI_ave_max= PI_ave_max/PI_weight_max;
	PI_ave_min= PI_ave_min/PI_weight_min;

	double den_weight= 1.0/(pow(PI_ave_max.err(),2))  + 1.0/(pow(PI_ave_min.err(),2));
	double fin_weight_max= (fit_to_constant)?(1.0/(pow(PI_ave_max.err(),2)*den_weight)):0.5;
	double fin_weight_min= (fit_to_constant)?(1.0/(pow(PI_ave_min.err(),2)*den_weight)):0.5;

	Tdatas_opt.push_back(tdata_opt);
	PI_per_ens.distr_list.push_back( fin_weight_min*PI_ave_min + fin_weight_max*PI_ave_max);
	Is_T_data_opt[tdata_opt-1] = 1.0;

      }


      //Print to File
      Print_To_File({}, {TCUTS, PI_q2_min_Tdata.ave(), PI_q2_min_Tdata.err(), PI_q2_max_Tdata.ave(), PI_q2_max_Tdata.err(), PI_q2_T_2.ave(), PI_q2_T_2.err(), Is_T_data_opt}, path+"_Q2_"+to_string_with_precision(Qs2[q],5)+".dat", "", "#tcut   lower    upper     T/2      Is_Tcut_opt.       Tcut_opt= "+to_string(Tdatas_opt[q]));


  }

  return;

}


void Bounding_PI_q2_disco(distr_t_list &PI_per_ens, const distr_t_list &V, const distr_t_list &Conn_guess, const distr_t &a, string path, Vint &Tdatas_opt, distr_t m_rho_GS) {

  Bounding_PI_q2( PI_per_ens, V, a, path, Tdatas_opt, m_rho_GS);

  PI_per_ens = PI_per_ens - Conn_guess;

  return;

}

void Get_amu_W_eps( vector<distr_t_list> &eps_win_list, const distr_t_list &V, const distr_t &a) {

  double T= V.size();

  vector<distr_t> eps_win_per_ens;

  
  auto K = [&](double Mv, double t, double size) -> double { return kernel_K(t, Mv);};

  distr_t_list Ker = distr_t_list::f_of_distr(K, a , T/2);


  for(int id_eps=0; id_eps<eps_win_size;id_eps++) {

    double eps= eps_win[id_eps]; 
  
    //define lambdas for the theta func
    auto th0 = [&eps](double ta) ->double { return 1.0/(1.0 + exp(-2.0*(ta-(t0-eps))/Delta));};
    auto th1 = [&eps](double ta) ->double { return 1.0/(1.0 + exp(-2.0*(ta-(t1-eps))/Delta));};

    auto amu_W_func_eps = [&V, &a,&th0, &th1, &Ker, &T]() -> distr_t {

			distr_t ret_win(UseJack,UseJack?Njacks:Nboots);
			for(int t=1;t<T/2;t++) ret_win = ret_win + 4.0*w(t,Simps_ord)*pow(alpha,2)*V.distr_list[t]*Ker.distr_list[t]*( distr_t::f_of_distr(th0, t*a) - distr_t::f_of_distr(th1, t*a));
			return ret_win;

		      };

    eps_win_list[id_eps].distr_list.push_back(amu_W_func_eps());

   
  }

  return;  

}


double Get_GS_artifacts_win(double M, double L, LL_functions &LL, bool mode_mixed=0) {

  double res, err;
  Vfloat Ergs_cont;
  Vfloat Ergs_lat;
  LL.Find_pipi_energy_lev(L, m_rho, 5.95, M, 0.0, Ergs_lat);
  LL.Find_pipi_energy_lev(L, m_rho, 5.95, 0.135, 0.0, Ergs_cont);

  auto Bt_coeff= [](double t) { return 0.5*(1.0 - exp(-pow( Get_4l_alpha_s( 1.0/t, 4, 0.340)/1.0,2)));};

  
  auto GS_art =[&LL, &M, &L, &Ergs_lat, &Ergs_cont, &Bt_coeff, &mode_mixed](double t) {
        	 
		 if(mode_mixed) return  4.0*pow(alpha,2)*(10.0/9.0)*(1.0-Bt_coeff(t))*(LL.V_pipi(t, L , m_rho, 5.95, 0.135, 0.0, Ergs_cont)  -LL.V_pipi(t, L , m_rho, 5.95, M, 0.0, Ergs_lat))*kernel_K(t, 1.0)*(  1.0/(1.0 + exp(-2.0*(t-t0)/Delta)) -  1.0/(1.0 + exp(-2.0*(t-t1)/Delta)));


		 return 4.0*pow(alpha,2)*(10.0/9.0)*( LL.V_pipi(t, L , m_rho, 5.95, 0.135, 0.0, Ergs_cont) - LL.V_pipi(t, L , m_rho, 5.95, M, 0.0, Ergs_lat))*kernel_K(t, 1.0)*(  1.0/(1.0 + exp(-2.0*(t-t0)/Delta)) -  1.0/(1.0 + exp(-2.0*(t-t1)/Delta)));
	       };


  double prec=1e-4;
  gsl_function_pp<decltype(GS_art)> gsl_GS_art(GS_art);
  gsl_integration_workspace *w_GS_art = gsl_integration_workspace_alloc(10000);
  gsl_function *G_GS_art = static_cast<gsl_function*>(&gsl_GS_art);
  gsl_integration_qagiu(G_GS_art, 0.3, 0.0, prec, 10000, w_GS_art, &res, &err);
  gsl_integration_workspace_free (w_GS_art);
  if(fabs(err/res) > 10*prec) crash("In Get_GS_artifacts_win precision reached is: "+to_string_with_precision( fabs(err/res), 5));
  

  return res;

  

}




void Gm2() {

  omp_set_num_threads(1);

  if(Mod=="FREE_THEORY") { Compute_SD_window_Free(); exit(-1);}

  if(Gen_free_corr_data) {Generate_free_corr_data(); exit(-1);}


  Plot_Energy_windows_K();

  Init_Qs2();

  Print_4l_alpha_s() ;

  


  //create directories g-2
  boost::filesystem::create_directory("../data/gm2");
  boost::filesystem::create_directory("../data/gm2/scale_setting");
  boost::filesystem::create_directory("../data/gm2/scale_setting/Mp");
  boost::filesystem::create_directory("../data/gm2/scale_setting/fp");
  boost::filesystem::create_directory("../data/gm2/light");
  boost::filesystem::create_directory("../data/gm2/light/mass_dep");
  boost::filesystem::create_directory("../data/gm2/light/mass_dep/tm");
  boost::filesystem::create_directory("../data/gm2/light/mass_dep/OS");
  boost::filesystem::create_directory("../data/gm2/light/GS_FSEs");
  boost::filesystem::create_directory("../data/gm2/light/autocorr");
  boost::filesystem::create_directory("../data/gm2/light/OS");
  boost::filesystem::create_directory("../data/gm2/light/tm");
  boost::filesystem::create_directory("../data/gm2/light/disco");
  boost::filesystem::create_directory("../data/gm2/strange");
  boost::filesystem::create_directory("../data/gm2/strange/mass_dep");
  boost::filesystem::create_directory("../data/gm2/strange/mass_dep/tm");
  boost::filesystem::create_directory("../data/gm2/strange/mass_dep/OS");
  boost::filesystem::create_directory("../data/gm2/strange/autocorr");
  boost::filesystem::create_directory("../data/gm2/strange/tm_phi");
  boost::filesystem::create_directory("../data/gm2/strange/OS_phi");
  boost::filesystem::create_directory("../data/gm2/strange/tm_etas");
  boost::filesystem::create_directory("../data/gm2/strange/OS_etas");
  boost::filesystem::create_directory("../data/gm2/strange/Z_phi");
  boost::filesystem::create_directory("../data/gm2/strange/Z_etas");
  boost::filesystem::create_directory("../data/gm2/strange/disco");
  boost::filesystem::create_directory("../data/gm2/charm");
  boost::filesystem::create_directory("../data/gm2/charm/autocorr");
  boost::filesystem::create_directory("../data/gm2/charm/tm_Jpsi");
  boost::filesystem::create_directory("../data/gm2/charm/OS_Jpsi");
  boost::filesystem::create_directory("../data/gm2/charm/tm_etac");
  boost::filesystem::create_directory("../data/gm2/charm/OS_etac");
  boost::filesystem::create_directory("../data/gm2/charm/Z_Jpsi");
  boost::filesystem::create_directory("../data/gm2/charm/Z_etac");
  boost::filesystem::create_directory("../data/gm2/charm/disco");
  boost::filesystem::create_directory("../data/gm2/light_strange");
  boost::filesystem::create_directory("../data/gm2/light_charm");
  boost::filesystem::create_directory("../data/gm2/strange_charm");
  boost::filesystem::create_directory("../data/gm2/light_strange/disco");
  boost::filesystem::create_directory("../data/gm2/light_charm/disco");
  boost::filesystem::create_directory("../data/gm2/strange_charm/disco");
  boost::filesystem::create_directory("../data/gm2/total_disco");
  boost::filesystem::create_directory("../data/gm2/total_disco/disco");

  //create directories PI(q^2)
  boost::filesystem::create_directory("../data/PI_Q2");
  boost::filesystem::create_directory("../data/PI_Q2/light");
  boost::filesystem::create_directory("../data/PI_Q2/light/tm");
  boost::filesystem::create_directory("../data/PI_Q2/light/OS");
  boost::filesystem::create_directory("../data/PI_Q2/light/disco");
  boost::filesystem::create_directory("../data/PI_Q2/strange");
  boost::filesystem::create_directory("../data/PI_Q2/strange/tm_phi");
  boost::filesystem::create_directory("../data/PI_Q2/strange/OS_phi");
  boost::filesystem::create_directory("../data/PI_Q2/strange/tm_etas");
  boost::filesystem::create_directory("../data/PI_Q2/strange/OS_etas");
  boost::filesystem::create_directory("../data/PI_Q2/strange/disco");
  boost::filesystem::create_directory("../data/PI_Q2/charm");
  boost::filesystem::create_directory("../data/PI_Q2/charm/tm_Jpsi");
  boost::filesystem::create_directory("../data/PI_Q2/charm/OS_Jpsi");
  boost::filesystem::create_directory("../data/PI_Q2/charm/tm_etac");
  boost::filesystem::create_directory("../data/PI_Q2/charm/OS_etac");
  boost::filesystem::create_directory("../data/PI_Q2/charm/disco");
  boost::filesystem::create_directory("../data/PI_Q2/light_strange");
  boost::filesystem::create_directory("../data/PI_Q2/light_charm");
  boost::filesystem::create_directory("../data/PI_Q2/strange_charm");
  boost::filesystem::create_directory("../data/PI_Q2/light_strange/disco");
  boost::filesystem::create_directory("../data/PI_Q2/light_charm/disco");
  boost::filesystem::create_directory("../data/PI_Q2/strange_charm/disco");
  boost::filesystem::create_directory("../data/PI_Q2/total_disco");
  boost::filesystem::create_directory("../data/PI_Q2/total_disco/disco");


 

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


  //print phi function and its derivative

  double Npoints=10000;
  Vfloat phi_sample, phi_prime_sample, zetas;
  for(int i=0;i<Npoints;i++) {
    zetas.push_back( i*0.001);
    phi_sample.push_back( LL.phi_spline( i*0.001));
    phi_prime_sample.push_back(LL.phi_der_spline( i*0.001));

  }

  Print_To_File({}, {zetas, phi_sample, phi_prime_sample}, "../data/gm2/phi.dat", "", "#z  phi   phi'");



  //test
  
  double vol = 6.272*fm_to_inv_Gev;

  double k= 0.065333333*fm_to_inv_Gev;
  Vfloat En_test;
  LL.Find_pipi_energy_lev(vol, m_rho, grpp_phys_BMW, 0.139, 0.0, En_test);

  for(int i=0; i<Nresonances;i++) cout<<"En: "<<2.0*sqrt( pow(En_test[i],2) + pow(0.139,2))<<" A: "<<LL.Amplitude(En_test[i],vol , m_rho, grpp_phys_BMW, 0.139, 0.0)<<endl;

  for(int t=1; t<19;t++) cout<<"V("<<t<<"): "<<LL.V_pipi(k*t, vol , m_rho, grpp_phys_BMW, 0.139, 0.0, En_test)<<" "<<LL.V_pipi_infL(k*t, m_rho, grpp_phys_BMW, 0.139, 0.0)<<endl;


  //end test


  

       

  //######################################################################################################


  //######################################################################################################
  if(Print_Mpi_dep_from_GSLL_param) {
    
    double dm= 0.0005;
    double Mp= 0.130;
    double L= 5.5*fm_to_inv_Gev;
    ofstream Print_Mpi("../data/gm2/light/GS_FSEs/Mp_dep.dat");
    Print_Mpi<<"#Mpi   L=5.5     L=infty"<<endl;
    while( Mp <= 0.149999) {
      Vfloat En_lev;
      LL.Find_pipi_energy_lev(L , 0.775, 5.95, Mp, 0.0, En_lev);
      auto W_L= [&Mp, &LL, &L, &En_lev] (double t) -> double {
		  double res_V_pipi = LL.V_pipi(t, L, 0.775, 5.95, Mp, 0.0, En_lev);
		  return 1e10*(10.0/9.0)*4.0*pow(alpha,2)*res_V_pipi*kernel_K(t, 1.0)*(  1.0/(1.0 + exp(-2.0*(t-t0)/Delta)) -  1.0/(1.0 + exp(-2.0*(t-t1)/Delta)));

		};
      auto W_inf = [&Mp, &LL](double t) -> double {
		 if(t<1e-1) return 0.0;    
		 double res_V_pipi_infL =  LL.V_pipi_infL(t, 0.775, 5.95, Mp, 0.0);
		 return 1e10*(10.0/9.0)*4.0*pow(alpha,2)*(res_V_pipi_infL)*kernel_K(t, 1.0)*(  1.0/(1.0 + exp(-2.0*(t-t0)/Delta)) -  1.0/(1.0 + exp(-2.0*(t-t1)/Delta)));
	       };

      double val_L, err_L;
      double val_inf, err_inf;
      double prec=1e-6;

      gsl_function_pp<decltype(W_L)> amu_W_L(W_L);
      gsl_integration_workspace *w_L = gsl_integration_workspace_alloc(10000);
      gsl_function *G_amu_W_L = static_cast<gsl_function*>(&amu_W_L);
      gsl_integration_qagiu(G_amu_W_L, 0.1, 0.0, prec, 10000, w_L, &val_L, &err_L);
      gsl_integration_workspace_free (w_L);

      gsl_function_pp<decltype(W_inf)> amu_W_inf(W_inf);
      gsl_integration_workspace *w_inf = gsl_integration_workspace_alloc(10000);
      gsl_function *G_amu_W_inf = static_cast<gsl_function*>(&amu_W_inf);
      gsl_integration_qagiu(G_amu_W_inf, 0.1, 0.0, prec, 10000, w_inf, &val_inf, &err_inf);
      gsl_integration_workspace_free (w_inf);

      Print_Mpi<<Mp<<"\t"<<val_L<<"\t"<<val_inf<<endl;


      Mp+=dm;
    }

    Print_Mpi.close();
    exit(-1);
  }

    



  


  
 
  //init Gaussian number generator
  GaussianMersenne GM(943832);


  GaussianMersenne GM_tests(54193);
  string channel="";

  
  data_t  V_light_1, V_light_2, V_light_3, pt2_pion, pt2_pion_B25, pt2_pion_A;
  data_t  V_light_OS_1, V_light_OS_2, V_light_OS_3;
  data_t corr_P5P5_strange, corr_P5P5_strange_heavy,  pt2_pion_charm, corr_P5P5_OS_strange, corr_P5P5_OS_strange_heavy,  pt2_pion_OS_charm;
  data_t V_light_1_m_small;
  data_t V_light_1_m_big;
  data_t V_light_OS_1_m_small;
  data_t V_light_OS_1_m_big;
  data_t P5P5_m_small;
  data_t P5P5_m_big;
  data_t correlated_V_tm_bubble, correlated_V_OS_bubble, chiral_condensate;
  data_t correlated_P5P5_tm_bubble, correlated_P5P5_OS_bubble;
  

  //L
  data_t  V_strange_1_L, V_strange_2_L, V_strange_3_L, V_strange_OS_1_L, V_strange_OS_2_L, V_strange_OS_3_L;
  data_t  V_charm_1_L, V_charm_2_L, V_charm_3_L, V_charm_OS_1_L, V_charm_OS_2_L, V_charm_OS_3_L;
  data_t  pt2_etaC_L, pt2_etaC_OS_L;
  //M
  data_t  V_strange_1_M, V_strange_2_M, V_strange_3_M, V_strange_OS_1_M, V_strange_OS_2_M, V_strange_OS_3_M;
  data_t  V_charm_1_M, V_charm_2_M, V_charm_3_M, V_charm_OS_1_M, V_charm_OS_2_M, V_charm_OS_3_M;
  data_t  pt2_etaC_M, pt2_etaC_OS_M;
  //H
  data_t  V_charm_1_H, V_charm_2_H, V_charm_3_H, V_charm_OS_1_H, V_charm_OS_2_H, V_charm_OS_3_H;
  data_t  pt2_etaC_H, pt2_etaC_OS_H;
  //to compute ZV 
  data_t corr_A0P5;
  //to compute ZA
  data_t corr_A0P5_OS, corr_P5P5_OS;
  //disco
  data_t disco_light, disco_strange, disco_charm;

  //disco_improved  (Simone) 
  data_t disco_impr_light, disco_impr_lightD, disco_impr_lightDD;
  data_t disco_impr_light1, disco_impr_light2, disco_impr_light4, disco_impr_light_no_rep;
  data_t disco_impr_charm, disco_impr_strange;
  data_t disco_impr_light_strange, disco_impr_lightD_strange;
  data_t disco_impr_light_charm, disco_impr_lightD_charm;
  data_t disco_impr_strange_charm;

  //to compute ZV (strange)
  data_t corr_A0P5_strange, corr_A0P5_strange_heavy;
  //to compute ZA (strange)
  data_t corr_A0P5_OS_strange, corr_A0P5_OS_strange_heavy;


  
  //to compute ZV and ZA (charm)
  data_t corr_A0P5_charm_L , corr_A0P5_OS_charm_L;
  data_t corr_A0P5_charm_M , corr_A0P5_OS_charm_M;
  data_t corr_A0P5_charm_H , corr_A0P5_OS_charm_H;
  data_t corr_A0P5_charm_pion, corr_A0P5_OS_charm_pion;
  
  //Read data

  //Custom sorting for V_light to account for the two replica r0 and r1, ...., rn

 
  
  auto Sort_light_confs = [](string A, string B) {

			   

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


  //#################################END CUSTOM SORTING#################
  V_light_1.Read("../gm2_data/light", "mes_contr_2pts_ll_1", "V1V1", Sort_light_confs);
  V_light_2.Read("../gm2_data/light", "mes_contr_2pts_ll_1", "V2V2", Sort_light_confs);
  V_light_3.Read("../gm2_data/light", "mes_contr_2pts_ll_1", "V3V3", Sort_light_confs);

  V_light_OS_1.Read("../gm2_data/light", "mes_contr_2pts_ll_2", "V1V1", Sort_light_confs);
  V_light_OS_2.Read("../gm2_data/light", "mes_contr_2pts_ll_2", "V2V2", Sort_light_confs);
  V_light_OS_3.Read("../gm2_data/light", "mes_contr_2pts_ll_2", "V3V3", Sort_light_confs);

  //#####################################
  //for test of mass dependence on V_light
  // standard
  V_light_1_m_small.Read("../gm2_new_run_correlated/physical", "mes_contr_2pts_ll_1", "VKVK", Sort_light_confs);
  V_light_1_m_big.Read("../gm2_new_run_correlated/simulated", "mes_contr_2pts_ll_1", "VKVK", Sort_light_confs);
  P5P5_m_small.Read("../gm2_new_run_correlated/physical", "mes_contr_2pts_ll_1", "P5P5", Sort_light_confs);
  P5P5_m_big.Read("../gm2_new_run_correlated/simulated", "mes_contr_2pts_ll_1", "P5P5", Sort_light_confs);
  
  V_light_OS_1_m_small.Read("../gm2_new_run_correlated/physical", "mes_contr_2pts_ll_2", "VKVK", Sort_light_confs);
  V_light_OS_1_m_big.Read("../gm2_new_run_correlated/simulated", "mes_contr_2pts_ll_2", "VKVK", Sort_light_confs);
  correlated_V_tm_bubble.Read("../gm2_sea_effects/conn", "mes_contr_2pts_ll_1", "V1V1", Sort_light_confs);
  correlated_V_OS_bubble.Read("../gm2_sea_effects/conn", "mes_contr_2pts_ll_2", "V1V1", Sort_light_confs);
  correlated_P5P5_tm_bubble.Read("../gm2_sea_effects/conn", "mes_contr_2pts_ll_1", "P5P5", Sort_light_confs);
  correlated_P5P5_OS_bubble.Read("../gm2_sea_effects/conn", "mes_contr_2pts_ll_2", "P5P5", Sort_light_confs);
  chiral_condensate.Read("../gm2_sea_effects/bubble", "disco", "", Sort_light_confs);
  

  //Use for Nazario
  /*
  V_light_1_m_small.Read("../mass_corrections_gm2/valence/physical", "mes_contr_2pts_ll_1", "VKVK", Sort_light_confs);
  V_light_1_m_big.Read("../mass_corrections_gm2/valence/simulated", "mes_contr_2pts_ll_1", "VKVK", Sort_light_confs);
  P5P5_m_small.Read("../mass_corrections_gm2/valence/physical", "mes_contr_2pts_ll_1", "P5P5", Sort_light_confs);
  P5P5_m_big.Read("../mass_corrections_gm2/valence/simulated", "mes_contr_2pts_ll_1", "P5P5", Sort_light_confs);
  
  V_light_OS_1_m_small.Read("../mass_corrections_gm2/valence/physical", "mes_contr_2pts_ll_2", "VKVK", Sort_light_confs);
  V_light_OS_1_m_big.Read("../mass_corrections_gm2/valence/simulated", "mes_contr_2pts_ll_2", "VKVK", Sort_light_confs);
  correlated_V_tm_bubble.Read("../mass_corrections_gm2/original", "mes_contr_2pts_ll_1", "V1V1", Sort_light_confs);
  correlated_V_OS_bubble.Read("../mass_corrections_gm2/original", "mes_contr_2pts_ll_2", "V1V1", Sort_light_confs);
  correlated_P5P5_tm_bubble.Read("../mass_corrections_gm2/original", "mes_contr_2pts_ll_1", "P5P5", Sort_light_confs);
  correlated_P5P5_OS_bubble.Read("../mass_corrections_gm2/original", "mes_contr_2pts_ll_2", "P5P5", Sort_light_confs);
  chiral_condensate.Read("../mass_corrections_gm2/bubble", "disco", "", Sort_light_confs);
  */
  //#####################################
  
  pt2_pion.Read("../gm2_data/light", "mes_contr_2pts_ll_1", "P5P5", Sort_light_confs);
  pt2_pion_B25.Read("../gm2_data/B25_light_ens", "mes_contr_2pts_ll_1", "P5P5", Sort_light_confs);
  pt2_pion_A.Read("../gm2_data/A_light_ens", "mes_contr_2pts_ll_1", "P5P5", Sort_light_confs);
  corr_A0P5.Read("../gm2_data/light", "mes_contr_2pts_ll_1", "P5A0", Sort_light_confs);
  corr_A0P5_OS.Read("../gm2_data/light", "mes_contr_2pts_ll_2", "P5A0", Sort_light_confs);
  corr_P5P5_OS.Read("../gm2_data/light", "mes_contr_2pts_ll_2", "P5P5", Sort_light_confs);

  //disco light
  if(Include_light_disco) {
    disco_light.Read("../gm2_data/disco_light/data", "disco", "", Sort_light_confs);

    //stochastic-stochastic
    disco_impr_light.Read("../gm2_data_disc_Simone/loops/data/light_light", "disco",  "", Sort_light_confs);
    //stochastic-deflated
    disco_impr_lightD.Read("../gm2_data_disc_Simone/loops/data/light_light_D", "disco",  "", Sort_light_confs);
    //deflated-deflated
    disco_impr_lightDD.Read("../gm2_data_disc_Simone/loops/data/light_D_light_D", "disco", "", Sort_light_confs);

    //light_1
    disco_impr_light1.Read("../gm2_data_disc_Simone/loops/data/light_light_1", "disco", "", Sort_light_confs);

    //light_2
    disco_impr_light2.Read("../gm2_data_disc_Simone/loops/data/light_light_2", "disco", "", Sort_light_confs);

    //light_4
    disco_impr_light4.Read("../gm2_data_disc_Simone/loops/data/light_light_4", "disco", "", Sort_light_confs);

    //light_no_replica

    disco_impr_light_no_rep.Read("../gm2_data_disc_Simone/loops/data/light_light_no_rep", "disco", "", Sort_light_confs);
    
  }

  
  
 

  //strange
  //L
  V_strange_1_L.Read("../gm2_data/strange_Nhits64/light", "mes_contr_2pts_ll_1", "V1V1", Sort_light_confs);
  V_strange_2_L.Read("../gm2_data/strange_Nhits64/light", "mes_contr_2pts_ll_1", "V2V2", Sort_light_confs);
  V_strange_3_L.Read("../gm2_data/strange_Nhits64/light", "mes_contr_2pts_ll_1", "V3V3", Sort_light_confs);
  V_strange_OS_1_L.Read("../gm2_data/strange_Nhits64/light", "mes_contr_2pts_ll_2", "V1V1", Sort_light_confs); 
  V_strange_OS_2_L.Read("../gm2_data/strange_Nhits64/light", "mes_contr_2pts_ll_2", "V2V2", Sort_light_confs);
  V_strange_OS_3_L.Read("../gm2_data/strange_Nhits64/light", "mes_contr_2pts_ll_2", "V3V3", Sort_light_confs);
  //M
  V_strange_1_M.Read("../gm2_data/strange_Nhits64/heavy", "mes_contr_2pts_ll_1", "V1V1", Sort_light_confs); 
  V_strange_2_M.Read("../gm2_data/strange_Nhits64/heavy", "mes_contr_2pts_ll_1", "V2V2", Sort_light_confs);
  V_strange_3_M.Read("../gm2_data/strange_Nhits64/heavy", "mes_contr_2pts_ll_1", "V3V3", Sort_light_confs);
  V_strange_OS_1_M.Read("../gm2_data/strange_Nhits64/heavy", "mes_contr_2pts_ll_2", "V1V1", Sort_light_confs); 
  V_strange_OS_2_M.Read("../gm2_data/strange_Nhits64/heavy", "mes_contr_2pts_ll_2", "V2V2", Sort_light_confs);
  V_strange_OS_3_M.Read("../gm2_data/strange_Nhits64/heavy", "mes_contr_2pts_ll_2", "V3V3", Sort_light_confs);
  //P5P5
  corr_P5P5_strange.Read("../gm2_data/strange_Nhits64/light", "mes_contr_2pts_ll_1", "P5P5", Sort_light_confs);
  corr_P5P5_OS_strange.Read("../gm2_data/strange_Nhits64/light", "mes_contr_2pts_ll_2", "P5P5", Sort_light_confs);
  corr_P5P5_strange_heavy.Read("../gm2_data/strange_Nhits64/heavy", "mes_contr_2pts_ll_1", "P5P5", Sort_light_confs);
  corr_P5P5_OS_strange_heavy.Read("../gm2_data/strange_Nhits64/heavy", "mes_contr_2pts_ll_2", "P5P5", Sort_light_confs); 
  //A0P5
  corr_A0P5_strange.Read("../gm2_data/strange_Nhits64/light", "mes_contr_2pts_ll_1", "P5A0", Sort_light_confs);
  corr_A0P5_OS_strange.Read("../gm2_data/strange_Nhits64/light", "mes_contr_2pts_ll_2", "P5A0", Sort_light_confs);
  corr_A0P5_strange_heavy.Read("../gm2_data/strange_Nhits64/heavy", "mes_contr_2pts_ll_1", "P5A0", Sort_light_confs);
  corr_A0P5_OS_strange_heavy.Read("../gm2_data/strange_Nhits64/heavy", "mes_contr_2pts_ll_2", "P5A0", Sort_light_confs);


  //disco strange
  if(Include_strange_disco) {
    disco_strange.Read("../gm2_data/disco_strange/data", "disco", "", Sort_light_confs);
    disco_impr_strange.Read("../gm2_data_disc_Simone/loops/data/strange_strange", "disco", "", Sort_light_confs);
  }
  

  //charm
  //L
  V_charm_1_L.Read("../gm2_data/charm_Nhits20/light", "mes_contr_2pts_ll_1", "V1V1",  Sort_light_confs); 
  V_charm_2_L.Read("../gm2_data/charm_Nhits20/light", "mes_contr_2pts_ll_1", "V2V2",  Sort_light_confs);
  V_charm_3_L.Read("../gm2_data/charm_Nhits20/light", "mes_contr_2pts_ll_1", "V3V3",  Sort_light_confs);
  V_charm_OS_1_L.Read("../gm2_data/charm_Nhits20/light", "mes_contr_2pts_ll_2", "V1V1",  Sort_light_confs); 
  V_charm_OS_2_L.Read("../gm2_data/charm_Nhits20/light", "mes_contr_2pts_ll_2", "V2V2",  Sort_light_confs);
  V_charm_OS_3_L.Read("../gm2_data/charm_Nhits20/light", "mes_contr_2pts_ll_2", "V3V3",  Sort_light_confs);
  //M
  V_charm_1_M.Read("../gm2_data/charm_Nhits20/medium", "mes_contr_2pts_ll_1", "V1V1",  Sort_light_confs); 
  V_charm_2_M.Read("../gm2_data/charm_Nhits20/medium", "mes_contr_2pts_ll_1", "V2V2",  Sort_light_confs);
  V_charm_3_M.Read("../gm2_data/charm_Nhits20/medium", "mes_contr_2pts_ll_1", "V3V3",  Sort_light_confs);
  V_charm_OS_1_M.Read("../gm2_data/charm_Nhits20/medium", "mes_contr_2pts_ll_2", "V1V1",  Sort_light_confs); 
  V_charm_OS_2_M.Read("../gm2_data/charm_Nhits20/medium", "mes_contr_2pts_ll_2", "V2V2",  Sort_light_confs);
  V_charm_OS_3_M.Read("../gm2_data/charm_Nhits20/medium", "mes_contr_2pts_ll_2", "V3V3",  Sort_light_confs);
  //H
  V_charm_1_H.Read("../gm2_data/charm_Nhits20/heavy", "mes_contr_2pts_ll_1", "V1V1",  Sort_light_confs);
  V_charm_2_H.Read("../gm2_data/charm_Nhits20/heavy", "mes_contr_2pts_ll_1", "V2V2",  Sort_light_confs);
  V_charm_3_H.Read("../gm2_data/charm_Nhits20/heavy", "mes_contr_2pts_ll_1", "V3V3",  Sort_light_confs);
  V_charm_OS_1_H.Read("../gm2_data/charm_Nhits20/heavy", "mes_contr_2pts_ll_2", "V1V1",  Sort_light_confs); 
  V_charm_OS_2_H.Read("../gm2_data/charm_Nhits20/heavy", "mes_contr_2pts_ll_2", "V2V2",  Sort_light_confs);
  V_charm_OS_3_H.Read("../gm2_data/charm_Nhits20/heavy", "mes_contr_2pts_ll_2", "V3V3",  Sort_light_confs);

  

  
  //P5P5
  pt2_pion_charm.Read("../gm2_data/light_silvano_A_ens", "mes_contr_2pts_ll_1", "P5P5",  Sort_light_confs);
  pt2_pion_OS_charm.Read("../gm2_data/light_silvano_A_ens", "mes_contr_2pts_ll_2", "P5P5",  Sort_light_confs);
 
  
  //eta_c
  //L
  pt2_etaC_L.Read("../gm2_data/charm_Nhits20/light", "mes_contr_2pts_ll_1", "P5P5",  Sort_light_confs);
  pt2_etaC_OS_L.Read("../gm2_data/charm_Nhits20/light", "mes_contr_2pts_ll_2", "P5P5",  Sort_light_confs);
  //M
  pt2_etaC_M.Read("../gm2_data/charm_Nhits20/medium", "mes_contr_2pts_ll_1", "P5P5",  Sort_light_confs);
  pt2_etaC_OS_M.Read("../gm2_data/charm_Nhits20/medium", "mes_contr_2pts_ll_2", "P5P5",  Sort_light_confs);
  //H
  pt2_etaC_H.Read("../gm2_data/charm_Nhits20/heavy", "mes_contr_2pts_ll_1", "P5P5",  Sort_light_confs);
  pt2_etaC_OS_H.Read("../gm2_data/charm_Nhits20/heavy", "mes_contr_2pts_ll_2", "P5P5",  Sort_light_confs);

  //A0P5 L
  corr_A0P5_charm_L.Read("../gm2_data/charm_Nhits20/light", "mes_contr_2pts_ll_1", "P5A0",  Sort_light_confs);
  corr_A0P5_OS_charm_L.Read("../gm2_data/charm_Nhits20/light", "mes_contr_2pts_ll_2", "P5A0",  Sort_light_confs);

  //A0P5 M
  corr_A0P5_charm_M.Read("../gm2_data/charm_Nhits20/medium", "mes_contr_2pts_ll_1", "P5A0",  Sort_light_confs);
  corr_A0P5_OS_charm_M.Read("../gm2_data/charm_Nhits20/medium", "mes_contr_2pts_ll_2", "P5A0",  Sort_light_confs);

  //A0P5 H
  corr_A0P5_charm_H.Read("../gm2_data/charm_Nhits20/heavy", "mes_contr_2pts_ll_1", "P5A0",  Sort_light_confs);
  corr_A0P5_OS_charm_H.Read("../gm2_data/charm_Nhits20/heavy", "mes_contr_2pts_ll_2", "P5A0",  Sort_light_confs);


  //A0P5 charm-pion
  corr_A0P5_charm_pion.Read("../gm2_data/light_silvano_A_ens", "mes_contr_2pts_ll_1", "P5A0", Sort_light_confs);
  corr_A0P5_OS_charm_pion.Read("../gm2_data/light_silvano_A_ens", "mes_contr_2pts_ll_2", "P5A0", Sort_light_confs);
  
  //disco charm
  if(Include_charm_disco) {
    disco_charm.Read("../gm2_data/disco_charm/data", "disco", "", Sort_light_confs);
    disco_impr_charm.Read("../gm2_data_disc_Simone/loops/data/charm_charm", "disco", "", Sort_light_confs);
  }
  



  //include off-diagonal disconnected
  if(Include_off_diagonal_disco) {

    //stochastic-stochastic
    disco_impr_light_strange.Read("../gm2_data_disc_Simone/loops/data/light_strange", "disco", "", Sort_light_confs);
    disco_impr_light_charm.Read("../gm2_data_disc_Simone/loops/data/light_charm", "disco", "", Sort_light_confs);
    disco_impr_strange_charm.Read("../gm2_data_disc_Simone/loops/data/strange_charm", "disco", "", Sort_light_confs);
    //stochastic-deflated
    disco_impr_lightD_strange.Read("../gm2_data_disc_Simone/loops/data/light_D_strange", "disco", "", Sort_light_confs);
    disco_impr_lightD_charm.Read("../gm2_data_disc_Simone/loops/data/light_D_charm", "disco",  "", Sort_light_confs);


  }
  
  
  
 



  //##################################################################################
  //###################    SCALE SETTING  USING w0, fp, Mp    ########################


 
  distr_t a_from_w0_A(UseJack), a_from_w0_B(UseJack), a_from_w0_C(UseJack), a_from_w0_D(UseJack);
  distr_t a_from_fp_A(UseJack), a_from_fp_B(UseJack), a_from_fp_C(UseJack), a_from_fp_D(UseJack);

  if(Determine_scale_setting) {

    distr_t_list Mpi_scale_setting(UseJack), fpi_scale_setting(UseJack), w0_scale_setting(UseJack);
    distr_t_list Mpi_scale_setting_phys_point_ens(UseJack), fpi_scale_setting_phys_point_ens(UseJack);
    distr_t_list Mpi_scale_setting_Bens(UseJack), fpi_scale_setting_Bens(UseJack);
    distr_t_list Mpi_scale_setting_Aens(UseJack), fpi_scale_setting_Aens(UseJack);
    vector<string> Ensemble_scale_setting_tag_list;
    vector<string> Ensemble_phys_point_tag_list;
    vector<string> Ensemble_B_tag_list;
    vector<string> Ensemble_A_tag_list;
    Vfloat L_scale_setting_list;
    Vfloat L_phys_point, L_B_ens, L_A_ens;
    //physical point ensemble
    for(int iens=0;iens< pt2_pion.size; iens++) {

      //Lattice info
      LatticeInfo L_info;
      L_info.LatInfo_new_ens(pt2_pion.Tag[iens]);
      CorrAnalysis Corr(UseJack, Njacks,Nboots);
      Corr.Nt = pt2_pion.nrows[iens];
      //Analyze correlators
      if(pt2_pion.Tag[iens].substr(1,12)=="B211b.072.96") {Corr.Tmin=30; Corr.Tmax=70;}
      else if(pt2_pion.Tag[iens].substr(1,12)=="B211b.072.64") { Corr.Tmin=27; Corr.Tmax=50;}
      else if(pt2_pion.Tag[iens].substr(1,1)=="C") {Corr.Tmin=40; Corr.Tmax=60;}
      else if(pt2_pion.Tag[iens].substr(1,1)=="D") {Corr.Tmin=41; Corr.Tmax=80;}
      else crash("In scale setting analysis cannot find Tmin,Tmax for ensemble: "+pt2_pion.Tag[iens]);
      distr_t_list pion_corr = Corr.corr_t(pt2_pion.col(0)[iens], "");
      if(pt2_pion.Tag[iens] !="cB211b.072.96") {
	auto LOG = [](double x) { return log(x);};
	distr_t_list delta_corr_sea_tm, delta_corr_sea_OS;
	distr_t_list delta_corr_sea_strange_tm, delta_corr_sea_strange_OS;


	//############  FIND SEA AND VALENCE ENSEMBLE ID #######################
	int id_ens=-1;
	for(int i_ens_chir=0; i_ens_chir< (signed)chiral_condensate.Tag.size();i_ens_chir++) {
	if( chiral_condensate.Tag[i_ens_chir] == pt2_pion.Tag[iens]) id_ens = i_ens_chir;
	}
	if(id_ens==-1) crash("In evaluating sea quark mistuning in scale setting: cannot find ensemble: "+pt2_pion.Tag[iens]);
	int id_ens_val=-1;
	for(int i_ens_val=0; i_ens_val<(signed)P5P5_m_big.Tag.size(); i_ens_val++) {
	  if( P5P5_m_big.Tag[i_ens_val] == pt2_pion.Tag[iens]) id_ens_val= i_ens_val;
	}
	if(id_ens_val==-1) crash("In evaluating valence quark mistuning in scale setting: cannot find ensemble: "+pt2_pion.Tag[iens]);
	//######################################################################
	

	//DEFINE DM ########################
	double dm=0, amphys=0;
	double alat_res_ave, alat_res_err;
	distr_t alat_res(UseJack);
        
	
	if( pt2_pion.Tag[iens] == "cB211b.072.64") { dm= (0.00072-0.0006675); amphys=0.0006675; alat_res_ave= 0.0795738*fm_to_inv_Gev  ; alat_res_err= 0.000131636*fm_to_inv_Gev ;}
	else if(pt2_pion.Tag[iens] == "cC211a.06.80") {dm = (0.00060-0.000585); amphys= 0.000585; alat_res_ave = 0.0682083*fm_to_inv_Gev  ; alat_res_err= 0.000134938*fm_to_inv_Gev; }
	else if( pt2_pion.Tag[iens] == "cD211a.054.96") {dm = (0.00054 -0.0004964); amphys= 0.0004964; alat_res_ave = 0.0569183*fm_to_inv_Gev; alat_res_err= 0.000115387*fm_to_inv_Gev;}
	else crash("DM shift in valence quark mistuning (scale setting) cannot be found ");
	for(int ijack=0;ijack<Njacks;ijack++) alat_res.distr.push_back( alat_res_ave + GM_tests()*alat_res_err/sqrt(Njacks -1.0));
	//##################################

        //sea effects
	VVfloat chiral_condensate_correlated_to_P5P5_tm(Corr.Nt);
	VVfloat chiral_condensate_correlated_to_P5P5_OS(Corr.Nt);
	//sanity check
	int Nconfs_tm_0 = correlated_P5P5_tm_bubble.col(0)[id_ens][0].size();
	int Nconfs_OS_0 = correlated_P5P5_OS_bubble.col(0)[id_ens][0].size();
	int Nconfs_chir_0 = chiral_condensate.col(1)[id_ens][0].size();
	if((Nconfs_tm_0 != Nconfs_OS_0) || (Nconfs_tm_0 != Nconfs_chir_0)) crash("In evaluating sea quark mistuning: number of confs between connected(tm,OS) and chiral condensate do not match");

	
	Vfloat chiral_condensate_averaged(Nconfs_chir_0,0.0);
	for(int t=0;t<Corr.Nt;t++) {
	  int Nconfs_chir = chiral_condensate.col(1)[id_ens][t].size();
	  if( (Nconfs_chir != Nconfs_chir_0)) crash("In evaluating sea quark mistuning: Nconfs chiral not constant over time");
	  for(int iconf=0; iconf<Nconfs_chir;iconf++) chiral_condensate_averaged[iconf] += -1.0*chiral_condensate.col(1)[id_ens][t][iconf];
	}

	Vfloat chiral_condensate_averaged_light;
	for(auto &chir: chiral_condensate_averaged) {
	  chiral_condensate_averaged_light.push_back( chir*dm);
	}
	
	auto F_disco_jack= [&](const Vfloat& par) { if((signed)par.size() != 3) crash("Lambda function F_disco_jack expects par[3], but par["+to_string((signed)par.size())+"] provided"); return par[0] -par[1]*par[2];};
      
      for(int t=0; t < Corr.Nt;t++) {

	//sanity chekcs
	int Nconfs_tm = correlated_P5P5_tm_bubble.col(0)[id_ens][t].size();
	int Nconfs_OS = correlated_P5P5_OS_bubble.col(0)[id_ens][t].size();
	if( (Nconfs_tm != Nconfs_tm_0)) crash("In evaluating sea quark mistuning: Nconfs tm P5P5 light not constant over time");
	if( (Nconfs_OS != Nconfs_OS_0)) crash("In evaluating sea quark mistuning: Nconfs OS P5P5 light not constant over time");
		
	for(int iconf=0; iconf<Nconfs_tm_0;iconf++) {
	  chiral_condensate_correlated_to_P5P5_tm[t].push_back( dm*correlated_P5P5_tm_bubble.col(0)[id_ens][t][iconf]*chiral_condensate_averaged[iconf]);
	  chiral_condensate_correlated_to_P5P5_OS[t].push_back( dm*correlated_P5P5_OS_bubble.col(0)[id_ens][t][iconf]*chiral_condensate_averaged[iconf]);
	}

	//jackknife
	Jackknife J_tm(10000,Njacks);
	delta_corr_sea_tm.distr_list.push_back(J_tm.DoJack(F_disco_jack, 3, chiral_condensate_correlated_to_P5P5_tm[t], chiral_condensate_averaged_light, correlated_P5P5_tm_bubble.col(0)[id_ens][t]));
	Jackknife J_OS(10000,Njacks);
	delta_corr_sea_OS.distr_list.push_back(J_OS.DoJack(F_disco_jack, 3, chiral_condensate_correlated_to_P5P5_OS[t], chiral_condensate_averaged_light, correlated_P5P5_OS_bubble.col(0)[id_ens][t]));
      }

      
      
      //Print to File
      Print_To_File({}, {delta_corr_sea_tm.ave(), delta_corr_sea_tm.err()}, "../data/gm2/light/mass_dep/tm/delta_corr_sea_pion_"+pt2_pion.Tag[iens]+".dat", "", "");
      Print_To_File({}, {delta_corr_sea_OS.ave(), delta_corr_sea_OS.err()}, "../data/gm2/light/mass_dep/OS/delta_corr_sea_pion_"+pt2_pion.Tag[iens]+".dat", "", "");
     
      distr_t_list P5P5_corr_phys_point_val= Corr.corr_t(P5P5_m_small.col(0)[id_ens_val], "") - Corr.corr_t(P5P5_m_big.col(0)[id_ens_val], "") + pion_corr;
      distr_t_list P5P5_corr_phys_point_sea= pion_corr + delta_corr_sea_tm;
      distr_t_list P5P5_corr_phys_point= P5P5_corr_phys_point_val + delta_corr_sea_tm;

      distr_t_list eff_mass_original= Corr.effective_mass_t( pion_corr, "../data/gm2/light/mass_dep/tm/eff_mass_original_"+pt2_pion.Tag[iens]+".dat");
      distr_t_list eff_mass_sea= Corr.effective_mass_t( P5P5_corr_phys_point_sea, "../data/gm2/light/mass_dep/tm/eff_mass_sea_"+pt2_pion.Tag[iens]+".dat");
      distr_t_list eff_mass_valence= Corr.effective_mass_t( P5P5_corr_phys_point_val, "../data/gm2/light/mass_dep/tm/eff_mass_valence_"+pt2_pion.Tag[iens]+".dat");
      distr_t_list eff_mass_corrected= Corr.effective_mass_t( P5P5_corr_phys_point, "../data/gm2/light/mass_dep/tm/eff_mass_corrected_"+pt2_pion.Tag[iens]+".dat");

      distr_t_list eff_mass_original_OS= Corr.effective_mass_t( pion_corr, "../data/gm2/light/mass_dep/OS/eff_mass_original_"+pt2_pion.Tag[iens]+".dat");
      distr_t_list eff_mass_sea_OS= Corr.effective_mass_t( P5P5_corr_phys_point_sea, "../data/gm2/light/mass_dep/tm/eff_mass_sea_"+pt2_pion.Tag[iens]+".dat");
      distr_t_list eff_mass_valence_OS= Corr.effective_mass_t( P5P5_corr_phys_point_val, "../data/gm2/light/mass_dep/tm/eff_mass_valence_"+pt2_pion.Tag[iens]+".dat");
      distr_t_list eff_mass_corrected_OS= Corr.effective_mass_t( P5P5_corr_phys_point, "../data/gm2/light/mass_dep/tm/eff_mass_corrected_"+pt2_pion.Tag[iens]+".dat");
      
	
      distr_t_list Mpi_eff_phys_point = Corr.effective_mass_t(P5P5_corr_phys_point, "../data/gm2/scale_setting/Mp/Mpi_phys_"+pt2_pion.Tag[iens]+".dat");
      distr_t_list Mpi_eff_phys_point_val= Corr.effective_mass_t(P5P5_corr_phys_point_val,"");
      distr_t_list Mpi_eff_phys_point_wo_bias_correction= Corr.effective_mass_t(P5P5_m_small.col(0)[id_ens_val], "../data/gm2/scale_setting/Mp/Mpi_phys_w0_bias_correction_"+pt2_pion.Tag[iens]+".dat");
      distr_t_list fpi_eff_phys_point = Corr.decay_constant_t(pow(2.0*amphys,2)*P5P5_corr_phys_point, "../data/gm2/scale_setting/fp/fpi_phys_"+pt2_pion.Tag[iens]+".dat");
      distr_t_list fpi_eff_phys_point_val= Corr.decay_constant_t( pow(2.0*amphys,2)*P5P5_corr_phys_point_val,"");
      cout<<"##################  AT ml=ml^phys on "<<pt2_pion.Tag[iens]<<" ensemble #########"<<endl;
      distr_t Mpi_fit_phys_point= Corr.Fit_distr(Mpi_eff_phys_point);
      distr_t Mpi_fit_phys_point_val= Corr.Fit_distr(Mpi_eff_phys_point_val);
      distr_t fpi_fit_phys_point= Corr.Fit_distr(fpi_eff_phys_point);
      distr_t fpi_fit_phys_point_val= Corr.Fit_distr(fpi_eff_phys_point_val);
      distr_t csi_fit_phys_point= Mpi_fit_phys_point*Mpi_fit_phys_point/(16*M_PI*M_PI*fpi_fit_phys_point*fpi_fit_phys_point);
      distr_t csi_fit_phys_point_val= Mpi_fit_phys_point_val*Mpi_fit_phys_point_val/(16*M_PI*M_PI*fpi_fit_phys_point_val*fpi_fit_phys_point_val);
      //GL and CDH corrections
      distr_t csi_L = Mpi_fit_phys_point*Mpi_fit_phys_point/(pow(4.0*M_PI,2)*fpi_fit_phys_point*fpi_fit_phys_point);
      distr_t csi_L_val= csi_fit_phys_point_val;
      distr_t g1 = distr_t::f_of_distr(g1_l, Mpi_fit_phys_point*L_info.L);
      distr_t g2 = distr_t::f_of_distr(g2_l, Mpi_fit_phys_point*L_info.L);
      distr_t g1_val = distr_t::f_of_distr(g1_l, Mpi_fit_phys_point_val*L_info.L);
      distr_t g2_val = distr_t::f_of_distr(g2_l, Mpi_fit_phys_point_val*L_info.L);
      double csi_true_phys= (1.0/(16.0*M_PI*M_PI))*pow(0.1350/0.1304,2);
      distr_t log_l = log(csi_true_phys) - distr_t::f_of_distr(LOG, csi_L);
      distr_t log_l_val = log(csi_true_phys) - distr_t::f_of_distr(LOG, csi_L_val);
      distr_t fpi_GL_phys_point, Mpi_GL_phys_point, csi_GL_phys_point;
      distr_t fpi_CDH_phys_point, Mpi_CDH_phys_point, csi_CDH_phys_point;
      distr_t fpi_CDH_phys_point_val, Mpi_CDH_phys_point_val, csi_CDH_phys_point_val;
      fpi_GL_phys_point= fpi_fit_phys_point/(1.0 -2.0*csi_L*g1);
      Mpi_GL_phys_point= Mpi_fit_phys_point/(1.0 + 0.5*csi_L*g1);
      distr_t fpi_isoQCD(UseJack);
      if(pt2_pion.Tag[iens] != "cC211a.06.80") {
      for(int ijack=0;ijack<Njacks;ijack++) fpi_isoQCD.distr.push_back( 0.1304 + GM()*2e-4/sqrt(Njacks-1));
      }
      else {
	for(int ijack=0;ijack<Njacks;ijack++) fpi_isoQCD.distr.push_back( 0.1304 + GM_tests()*2e-4/(Njacks-1));
      }
      csi_GL_phys_point= Mpi_GL_phys_point*Mpi_GL_phys_point/(16.0*M_PI*M_PI*fpi_GL_phys_point*fpi_GL_phys_point);
      fpi_CDH_phys_point = fpi_fit_phys_point/(1.0 -2.0*csi_L*g1 +2.0*csi_L*csi_L*( (Cf1(l1ph,l2ph,l3ph,l4ph) + Sf1(s0,s1,s2,s3) + Cf1_log()*log_l)*g1 + (Cf2(l1ph,l2ph,l3ph,l4ph) + Sf2(s0,s1,s2,s3) + Cf2_log()*log_l)*g2));
      Mpi_CDH_phys_point = (Mpi_fit_phys_point/(1.0 + 0.5*csi_L*g1 - csi_L*csi_L*( (Cm1(l1ph,l2ph,l3ph,l4ph) + Sm1(s0,s1,s2,s3) + Cm1_log()*log_l)*g1 + (Cm2(l1ph,l2ph,l3ph,l4ph) + Sm2(s0,s1,s2,s3) + Cm2_log()*log_l)*g2)));
      csi_CDH_phys_point = Mpi_CDH_phys_point*Mpi_CDH_phys_point/(16.0*M_PI*M_PI*fpi_CDH_phys_point*fpi_CDH_phys_point);

      fpi_CDH_phys_point_val = fpi_fit_phys_point_val/(1.0 -2.0*csi_L_val*g1_val +2.0*csi_L_val*csi_L_val*( (Cf1(l1ph,l2ph,l3ph,l4ph) + Sf1(s0,s1,s2,s3) + Cf1_log()*log_l_val)*g1_val + (Cf2(l1ph,l2ph,l3ph,l4ph) + Sf2(s0,s1,s2,s3) + Cf2_log()*log_l_val)*g2_val));
      Mpi_CDH_phys_point_val = (Mpi_fit_phys_point_val/(1.0 + 0.5*csi_L_val*g1_val - csi_L_val*csi_L_val*( (Cm1(l1ph,l2ph,l3ph,l4ph) + Sm1(s0,s1,s2,s3) + Cm1_log()*log_l_val)*g1_val + (Cm2(l1ph,l2ph,l3ph,l4ph) + Sm2(s0,s1,s2,s3) + Cm2_log()*log_l_val)*g2_val)));
      csi_CDH_phys_point_val = Mpi_CDH_phys_point_val*Mpi_CDH_phys_point_val/(16.0*M_PI*M_PI*fpi_CDH_phys_point_val*fpi_CDH_phys_point_val);

	
	cout<<"af: "<< fpi_fit_phys_point.ave()<< "+- "<<fpi_fit_phys_point.err()<<endl;
	cout<<"af: (only valence) "<< fpi_fit_phys_point_val.ave()<< "+- "<<fpi_fit_phys_point_val.err()<<endl;
	cout<<"af(GL-corrected): "<<fpi_GL_phys_point.ave()<<" +- "<<fpi_GL_phys_point.err()<<endl;
	cout<<"af(CDH-corrected): "<<fpi_CDH_phys_point.ave()<< "+- "<<fpi_CDH_phys_point.err()<<endl;
	cout<<"af(CDH-corrected) (only valence): "<<fpi_CDH_phys_point_val.ave()<< "+- "<<fpi_CDH_phys_point_val.err()<<endl;
	cout<<"aMpi: "<< (Mpi_fit_phys_point/alat_res).ave()<< "+- "<<(Mpi_fit_phys_point/alat_res).err()<<endl;
	cout<<"aMpi: (only valence) "<< (Mpi_fit_phys_point_val/alat_res).ave()<< "+- "<<(Mpi_fit_phys_point_val/alat_res).err()<<endl;
	cout<<"aMpi(GL-corrected): "<<(Mpi_GL_phys_point/alat_res).ave()<<" +- "<<(Mpi_GL_phys_point/alat_res).err()<<endl;
	cout<<"aMpi(CDH-corrected): "<<(Mpi_CDH_phys_point/alat_res).ave()<<" +- "<<(Mpi_CDH_phys_point/alat_res).err()<<endl;
	cout<<"aMpi(CDH-corrected) (only valence): "<<(Mpi_CDH_phys_point_val/alat_res).ave()<<" +- "<<(Mpi_CDH_phys_point_val/alat_res).err()<<endl;
	cout<<"xi: "<<csi_fit_phys_point.ave()<< "+- "<<csi_fit_phys_point.err()<<endl;
	cout<<"xi: (only valence) "<<csi_fit_phys_point_val.ave()<< "+- "<<csi_fit_phys_point_val.err()<<endl;
	cout<<"xi(GL-corrected): "<<csi_GL_phys_point.ave()<<" +- "<<csi_GL_phys_point.err()<<endl;
	cout<<"xi(CDH-corrected): "<<csi_CDH_phys_point.ave()<< "+- "<<csi_CDH_phys_point.err()<<endl;
	cout<<"xi(CDH-corrected) (only valence): "<<csi_CDH_phys_point_val.ave()<< "+- "<<csi_CDH_phys_point_val.err()<<endl;
	cout<<"xi_isoQCD: "<<csi_true_phys<<" +- 0.00002 "<<endl;
	cout<<"Prediction for lattice spacing: "<<(fpi_CDH_phys_point/(fpi_isoQCD*fm_to_inv_Gev)).ave()<<" +- "<<(fpi_CDH_phys_point/(fm_to_inv_Gev*fpi_isoQCD)).err()<<endl;
	cout<<"Prediction for lattice spacing(only valence): "<<(fpi_CDH_phys_point_val/(fpi_isoQCD*fm_to_inv_Gev)).ave()<<" +- "<<(fpi_CDH_phys_point_val/(fm_to_inv_Gev*fpi_isoQCD)).err()<<endl;
	cout<<"Global fit estimate: ";
	if(pt2_pion.Tag[iens].substr(1,1)=="B") cout<< "0.0795738 +- 0.000131636 [fm] "<<endl;
	else if(pt2_pion.Tag[iens].substr(1,1)=="C") cout<<"0.0682083 +- 0.000134938 [fm]   "<<endl;
        else if (pt2_pion.Tag[iens].substr(1,1)=="D") cout<<"0.0569183 +- .000115387 [fm] "<<endl;
	cout<<"#########################################################################"<<endl;

      }
      
      distr_t_list Mpi_eff_distr = Corr.effective_mass_t(pt2_pion.col(0)[iens], "../data/gm2/scale_setting/Mp/Mpi_"+pt2_pion.Tag[iens]+".dat");
      distr_t_list fpi_eff_distr = Corr.decay_constant_t(pow(2.0*L_info.ml,2)*pion_corr, "../data/gm2/scale_setting/fp/fpi_"+pt2_pion.Tag[iens]+".dat");
      distr_t Mpi_fit = Corr.Fit_distr(Mpi_eff_distr);
      distr_t fpi_fit = Corr.Fit_distr(fpi_eff_distr);
      Mpi_scale_setting.distr_list.push_back(Mpi_fit);
      fpi_scale_setting.distr_list.push_back(fpi_fit);
      Ensemble_scale_setting_tag_list.push_back(pt2_pion.Tag[iens]);
      L_scale_setting_list.push_back( L_info.L);
      Mpi_scale_setting_phys_point_ens.distr_list.push_back(Mpi_fit);
      fpi_scale_setting_phys_point_ens.distr_list.push_back(fpi_fit);
      Ensemble_phys_point_tag_list.push_back(pt2_pion.Tag[iens]);
      L_phys_point.push_back( L_info.L);
      if(pt2_pion.Tag[iens].substr(1,1) == "B") {
	Mpi_scale_setting_Bens.distr_list.push_back( Mpi_fit);
	fpi_scale_setting_Bens.distr_list.push_back( fpi_fit);
	Ensemble_B_tag_list.push_back(pt2_pion.Tag[iens]);
	L_B_ens.push_back(L_info.L);
      }
      //Read w0/a
      distr_t w0a_distr(UseJack);
      Vfloat confs_aw0;
      if(scale_setting_from_w0X) {
      confs_aw0 = Read_From_File("../gm2_data/data_w0/"+pt2_pion.Tag[iens], 1, 2);
      Compute_autocorrelation_time(confs_aw0, "../data/gm2/scale_setting", "w0A_"+pt2_pion.Tag[iens]);
      }
    }


    //B25 ensembles
    for(int iens=0;iens< pt2_pion_B25.size; iens++) {

      //Lattice info
      LatticeInfo L_info;
      L_info.LatInfo_new_ens(pt2_pion_B25.Tag[iens]);
      CorrAnalysis Corr(UseJack, Njacks,Nboots);
      Corr.Nt = pt2_pion_B25.nrows[iens];
      if(pt2_pion_B25.Tag[iens].substr(1,11)=="B211a.25.24") {Corr.Tmin=15; Corr.Tmax=20;}
      else if(pt2_pion_B25.Tag[iens].substr(1,11)=="B211a.25.32") {Corr.Tmin=23; Corr.Tmax=30;}
      else if(pt2_pion_B25.Tag[iens].substr(1,11)=="B211a.25.48") {Corr.Tmin=25; Corr.Tmax=44;}
      else if(pt2_pion_B25.Tag[iens].substr(1,11)=="B211a.14.64") {Corr.Tmin = 39; Corr.Tmax=52;}
      else crash("In scale setting analysis cannot find Tmin,Tmax for ensemble: "+pt2_pion_B25.Tag[iens]);
      distr_t_list pion_corr = Corr.corr_t(pt2_pion_B25.col(0)[iens], "");
      distr_t_list Mpi_eff_distr = Corr.effective_mass_t(pt2_pion_B25.col(0)[iens], "../data/gm2/scale_setting/Mp/Mpi_"+pt2_pion_B25.Tag[iens]+".dat");
      distr_t_list fpi_eff_distr = Corr.decay_constant_t(pow(2.0*L_info.ml,2)*pion_corr, "../data/gm2/scale_setting/fp/fpi_"+pt2_pion_B25.Tag[iens]+".dat");
      distr_t Mpi_fit = Corr.Fit_distr(Mpi_eff_distr);
      distr_t fpi_fit = Corr.Fit_distr(fpi_eff_distr);
      Mpi_scale_setting.distr_list.push_back(Mpi_fit);
      fpi_scale_setting.distr_list.push_back(fpi_fit);
      Ensemble_scale_setting_tag_list.push_back(pt2_pion_B25.Tag[iens]);
      Mpi_scale_setting_Bens.distr_list.push_back(Mpi_fit);
      fpi_scale_setting_Bens.distr_list.push_back(fpi_fit);
      Ensemble_B_tag_list.push_back( pt2_pion_B25.Tag[iens]);
      L_scale_setting_list.push_back( L_info.L);
      L_B_ens.push_back( L_info.L);
      //Read w0/a
      distr_t w0a_distr(UseJack);
      Vfloat confs_aw0;
      if(scale_setting_from_w0X) {
      confs_aw0 = Read_From_File("../gm2_data/data_w0/"+pt2_pion_B25.Tag[iens], 1, 2);
      Compute_autocorrelation_time(confs_aw0, "../data/gm2/scale_setting", "w0A_"+pt2_pion_B25.Tag[iens]);
      }
    }

  

    //A ensembles
    for(int iens=0;iens< pt2_pion_A.size; iens++) {

      if(pt2_pion_A.Tag[iens].substr(1,1) == "A" && !Use_three_finest_in_scale_setting_w0X) {
      //Lattice info
      LatticeInfo L_info;
      L_info.LatInfo_new_ens(pt2_pion_A.Tag[iens]);
      CorrAnalysis Corr(UseJack, Njacks,Nboots);
      Corr.Nt = pt2_pion_A.nrows[iens];
      if(pt2_pion_A.Tag[iens] == "cA211a.40.24") {Corr.Tmin=15; Corr.Tmax=21;}
      else if(pt2_pion_A.Tag[iens] == "cA211a.53.24") {Corr.Tmin=15; Corr.Tmax=22;}
      else if(pt2_pion_A.Tag[iens] == "cA211ab.30.32") {Corr.Tmin=16; Corr.Tmax=29;}
      else if(pt2_pion_A.Tag[iens] == "cA211a.12.48") {Corr.Tmin=21; Corr.Tmax=34;}
      else crash("In scale setting analysis cannot find Tmin,Tmax for ensemble: "+pt2_pion_A.Tag[iens]);
      distr_t_list pion_corr = Corr.corr_t(pt2_pion_A.col(0)[iens], "");
      distr_t_list Mpi_eff_distr = Corr.effective_mass_t(pt2_pion_A.col(0)[iens], "../data/gm2/scale_setting/Mp/Mpi_"+pt2_pion_A.Tag[iens]+".dat");
      distr_t_list fpi_eff_distr = Corr.decay_constant_t(pow(2.0*L_info.ml,2)*pion_corr, "../data/gm2/scale_setting/fp/fpi_"+pt2_pion_A.Tag[iens]+".dat");
      distr_t Mpi_fit = Corr.Fit_distr(Mpi_eff_distr);
      distr_t fpi_fit = Corr.Fit_distr(fpi_eff_distr);
      Mpi_scale_setting.distr_list.push_back(Mpi_fit);
      fpi_scale_setting.distr_list.push_back(fpi_fit);
      Mpi_scale_setting_Aens.distr_list.push_back(Mpi_fit);
      fpi_scale_setting_Aens.distr_list.push_back(fpi_fit);
      Ensemble_A_tag_list.push_back( pt2_pion_A.Tag[iens]);
      L_scale_setting_list.push_back( L_info.L);
      L_A_ens.push_back( L_info.L);
      }
      if(pt2_pion_A.Tag[iens].substr(1,1) == "A") {
	Ensemble_scale_setting_tag_list.push_back(pt2_pion_A.Tag[iens]);
	//Read w0/a
	distr_t w0a_distr(UseJack);
	Vfloat confs_aw0;
	if(scale_setting_from_w0X) {
	confs_aw0 = Read_From_File("../gm2_data/data_w0/"+pt2_pion_A.Tag[iens], 1, 2);
	Compute_autocorrelation_time(confs_aw0, "../data/gm2/scale_setting", "w0A_"+pt2_pion_A.Tag[iens]);
	}
      }
    }
  
 
  
    if(scale_setting_from_w0X) {
    Determine_scale_from_w0_fp(Mpi_scale_setting, fpi_scale_setting, w0_scale_setting, L_scale_setting_list, Ensemble_scale_setting_tag_list, a_from_w0_A, a_from_w0_B, a_from_w0_C, a_from_w0_D, UseJack, Use_three_finest_in_scale_setting_w0X);
    }
    else if(scale_setting_from_fp) {
      Determine_scale_from_fp(Mpi_scale_setting_phys_point_ens, fpi_scale_setting_phys_point_ens, Mpi_scale_setting_Bens, fpi_scale_setting_Bens, Mpi_scale_setting_Aens, fpi_scale_setting_Aens, L_phys_point, L_B_ens, L_A_ens, Ensemble_phys_point_tag_list, Ensemble_B_tag_list, Ensemble_A_tag_list, a_from_fp_A, a_from_fp_B, a_from_fp_C, a_from_fp_D, UseJack, Use_three_finest_in_scale_setting_fp);
    }
    else crash("scale_setting_from_w0X is not w0X or fp, but Determine_scale_setting = true");
  }
  //###################################################################################



  //###########################################################################################
  //#######################  generate fake jack_distr for lattice spacing a_A a_B, a_C, a_D ###########################
  LatticeInfo a_info;
  distr_t a_A(UseJack), a_B(UseJack), a_C(UseJack), a_D(UseJack);
  distr_t a_A_afp(UseJack), a_B_afp(UseJack), a_C_afp(UseJack), a_D_afp(UseJack);

 
  double a_A_ave, a_A_err, a_B_ave, a_B_err, a_C_ave, a_C_err, a_D_ave, a_D_err;
  double a_A_afp_ave, a_A_afp_err, a_B_afp_ave, a_B_afp_err, a_C_afp_ave, a_C_afp_err, a_D_afp_ave, a_D_afp_err;
  a_info.LatInfo_new_ens("cA211a.53.24");
  a_A_ave= a_info.a;
  a_A_err= a_info.a_err;
  a_A_afp_ave = a_info.a_from_afp;
  a_A_afp_err = a_info.a_from_afp_err;
  a_info.LatInfo_new_ens("cB211b.072.64");
  a_B_ave= a_info.a;
  a_B_err= a_info.a_err;
  a_B_afp_ave = a_info.a_from_afp;
  a_B_afp_err = a_info.a_from_afp_err;
  a_info.LatInfo_new_ens("cC211a.06.80");
  a_C_ave= a_info.a;
  a_C_err= a_info.a_err;
  a_C_afp_ave = a_info.a_from_afp;
  a_C_afp_err = a_info.a_from_afp_err;
  a_info.LatInfo_new_ens("cD211a.054.96");
  a_D_ave= a_info.a;
  a_D_err= a_info.a_err;
  a_D_afp_ave = a_info.a_from_afp;
  a_D_afp_err = a_info.a_from_afp_err;
  
  if(UseJack)  { for(int ijack=0;ijack<Njacks;ijack++) {
      a_A.distr.push_back( fm_to_inv_Gev*( a_A_ave + GM()*a_A_err*(1.0/sqrt(Njacks-1.0))));
      a_B.distr.push_back( fm_to_inv_Gev*( a_B_ave + GM()*a_B_err*(1.0/sqrt(Njacks-1.0))));
      a_C.distr.push_back( fm_to_inv_Gev*( a_C_ave + GM()*a_C_err*(1.0/sqrt(Njacks-1.0))));
      a_D.distr.push_back( fm_to_inv_Gev*( a_D_ave + GM()*a_D_err*(1.0/sqrt(Njacks-1.0))));
      a_A_afp.distr.push_back( fm_to_inv_Gev*( a_A_afp_ave+ 2*GM()*a_A_afp_err*(1.0/sqrt(Njacks -1.0))));
      a_B_afp.distr.push_back( fm_to_inv_Gev*( a_B_afp_ave+ 2*GM()*a_B_afp_err*(1.0/sqrt(Njacks -1.0))));
      a_C_afp.distr.push_back( fm_to_inv_Gev*( a_C_afp_ave+ 2*GM()*a_C_afp_err*(1.0/sqrt(Njacks -1.0))));
      a_D_afp.distr.push_back( fm_to_inv_Gev*( a_D_afp_ave+ 2*GM()*a_D_afp_err*(1.0/sqrt(Njacks -1.0))));

    }
    }
  else {
    for (int iboot=0; iboot<Nboots;iboot++) {
      a_A.distr.push_back( fm_to_inv_Gev*( a_A_ave + GM()*a_A_err));
      a_B.distr.push_back( fm_to_inv_Gev*( a_B_ave + GM()*a_B_err));
      a_C.distr.push_back( fm_to_inv_Gev*( a_C_ave + GM()*a_C_err));
      a_D.distr.push_back( fm_to_inv_Gev*( a_D_ave + GM()*a_D_err));
      a_A_afp.distr.push_back( fm_to_inv_Gev*( a_A_afp_ave+ 2*GM()*a_A_afp_err));
      a_B_afp.distr.push_back( fm_to_inv_Gev*( a_B_afp_ave+ 2*GM()*a_B_afp_err));
      a_C_afp.distr.push_back( fm_to_inv_Gev*( a_C_afp_ave+ 2*GM()*a_C_afp_err));
      a_D_afp.distr.push_back( fm_to_inv_Gev*( a_D_afp_ave+ 2*GM()*a_D_afp_err));
    }
  }
  
  

  //Set a_X to the one determined in this scale setting analysis, if required.
  if(Use_scale_setting_from_this_analysis) {
    if(scale_setting_from_w0X) {
      a_A = a_from_w0_A;
      a_B = a_from_w0_B;
      a_C = a_from_w0_C;
      a_D = a_from_w0_D;
    }
    else if (scale_setting_from_fp) {
      a_A = a_from_fp_A;
      a_B = a_from_fp_B;
      a_C = a_from_fp_C;
      a_D = a_from_fp_D;
    }
    else crash("Use_scale_setting_from_this_analysis = true but scale_setting_from_w0X and scale_setting_from_fp are both false");
  }


  cout<<"###################################################"<<endl;
  cout<<"Info on lattice spacing distribution used: "<<endl;
  cout<<"Ensemble A: (quark mass paper): ("<<a_A_ave*fm_to_inv_Gev<<" +- "<<a_A_err*fm_to_inv_Gev<<") [GeV-1] used: ("<<a_A.ave()<<" +- "<<a_A.err()<<") [GeV-1]"<<endl;
  cout<<"Ensemble B: (quark mass paper): ("<<a_B_ave*fm_to_inv_Gev<<" +- "<<a_B_err*fm_to_inv_Gev<<") [GeV-1] used: ("<<a_B.ave()<<" +- "<<a_B.err()<<") [GeV-1]"<<endl;
  cout<<"Ensemble C: (quark mass paper): ("<<a_C_ave*fm_to_inv_Gev<<" +- "<<a_C_err*fm_to_inv_Gev<<") [GeV-1] used: ("<<a_C.ave()<<" +- "<<a_C.err()<<") [GeV-1]"<<endl;
  cout<<"Ensemble D: (quark mass paper): ("<<a_D_ave*fm_to_inv_Gev<<" +- "<<a_D_err*fm_to_inv_Gev<<") [GeV-1] used: ("<<a_D.ave()<<" +- "<<a_D.err()<<") [GeV-1]"<<endl;
  cout<<"###################################################"<<endl;

 

  //#####################################################################################################################


 
  int Nens_light= V_light_1.size;
  int Nens_strange = V_strange_1_L.size;
  int Nens_charm = V_charm_1_L.size;
  int Nens_disco_light = disco_light.size;
  int Nens_disco_strange= disco_strange.size;
  int Nens_disco_charm = disco_charm.size;
  int Nens_light_silvano = pt2_pion_charm.size;
  vector<string> disco_light_Tags;
  vector<string> disco_impr_light_Tags;
  vector<string> disco_impr_lightD_Tags;
  vector<string> disco_impr_lightDD_Tags;
  vector<string> disco_strange_Tags;
  vector<string> disco_impr_strange_Tags;
  vector<string> disco_charm_Tags;
  vector<string> disco_impr_charm_Tags;
  //off-diagonal
  vector<string> disco_impr_strange_charm_Tags;
  vector<string> disco_impr_light_charm_Tags;
  vector<string> disco_impr_lightD_charm_Tags;
  vector<string> disco_impr_light_strange_Tags;
  vector<string> disco_impr_lightD_strange_Tags;
  cout<<"N_ens light: "<<Nens_light<<endl;
  cout<<"N_ens disco_light: "<<Nens_disco_light<<endl;
  cout<<"N_ens strange: "<<Nens_strange<<endl;
  cout<<"N_ens disco_strange: "<<Nens_disco_strange<<endl;
  cout<<"N_ens charm: "<<Nens_charm<<endl;
  cout<<"N_ens disco_charm: "<<Nens_disco_charm<<endl;
  cout<<"N_ens silvano light: "<<Nens_light_silvano<<endl;
  cout<<"Printing Number of configs for each ensemble: "<<endl;
  cout<<"Ensemble light: "<<endl;
  for(int ins=0; ins<Nens_light;ins++) cout<<V_light_1.Tag[ins]<<" : "<<V_light_1.Nconfs[ins]<<" Used: "<<Njacks*( V_light_1.Nconfs[ins]/Njacks)<<endl;
  cout<<"Ensemble strange: "<<endl;
  for(int ins=0; ins<Nens_strange;ins++) cout<<V_strange_1_L.Tag[ins]<<" : "<<V_strange_1_L.Nconfs[ins]<<" Used: "<<Njacks*( V_strange_1_L.Nconfs[ins]/Njacks)<<endl;
  cout<<"Ensemble charm: "<<endl;
  for(int ins=0; ins<Nens_charm;ins++) cout<<V_charm_1_L.Tag[ins]<<" : "<<V_charm_1_L.Nconfs[ins]<<" Used: "<<Njacks*( V_charm_1_L.Nconfs[ins]/Njacks)<<endl;



 

  //####################################### GET AUTOCORRELATION FOR VECTOR CORRELATOR #####################################
  for(int ins=0;ins<Nens_light;ins++) {
    int Tmax= V_light_1.col(0)[ins].size();
    for(int t=0;t<Tmax/2;t++) {
      Compute_autocorrelation_time(V_light_1.col(0)[ins][t],"../data/gm2/light/autocorr", "VV_t_"+to_string(t)+"_"+V_light_1.Tag[ins]);  
    }
  }
  //#######################################################################################################################
  



  //define distr_t_list to be used in chiral+continuum analysis

  //eps-light-analysis

  vector<distr_t_list> amu_W_eps_list_tm(eps_win_size);
  vector<distr_t_list> amu_W_eps_list_OS(eps_win_size);
  vector<distr_t_list> amu_W_eps_val_list_tm(eps_win_size);
  vector<distr_t_list> amu_W_eps_val_list_OS(eps_win_size);
  vector<distr_t_list> amu_W_eps_sea_list_tm(eps_win_size);
  vector<distr_t_list> amu_W_eps_sea_list_OS(eps_win_size);


  //GS-Sakurai prediction for cut-off effects
  Vfloat GS_artifacts_W_win_tm;
  Vfloat GS_artifacts_W_win_OS;

  //strange and charm
  //L
  distr_t_list agm2_strange_L(UseJack), agm2_charm_L(UseJack), agm2_strange_OS_L(UseJack), agm2_charm_OS_L(UseJack);
  //M
  distr_t_list agm2_strange_M(UseJack), agm2_charm_M(UseJack), agm2_strange_OS_M(UseJack), agm2_charm_OS_M(UseJack);
  //H
  distr_t_list  agm2_charm_H(UseJack), agm2_charm_OS_H(UseJack);
  //Extr
  distr_t_list agm2_strange_Extr(UseJack), agm2_charm_Extr(UseJack), agm2_strange_OS_Extr(UseJack), agm2_charm_OS_Extr(UseJack);

  //disco strange and charm
  //distr_t_list agm2_disco_strange(UseJack), agm2_disco_charm(UseJack);

  //strange and charm NON_ELM
  //L
  distr_t_list agm2_strange_No_ELM_L(UseJack), agm2_charm_No_ELM_L(UseJack), agm2_strange_OS_No_ELM_L(UseJack), agm2_charm_OS_No_ELM_L(UseJack);
  //M
  distr_t_list agm2_strange_No_ELM_M(UseJack), agm2_charm_No_ELM_M(UseJack), agm2_strange_OS_No_ELM_M(UseJack), agm2_charm_OS_No_ELM_M(UseJack);
  //H
  distr_t_list  agm2_charm_No_ELM_H(UseJack),  agm2_charm_OS_No_ELM_H(UseJack);
  //Extr
  distr_t_list agm2_strange_No_ELM_Extr(UseJack), agm2_charm_No_ELM_Extr(UseJack), agm2_strange_OS_No_ELM_Extr(UseJack), agm2_charm_OS_No_ELM_Extr(UseJack);

  //disco strange and charm
  distr_t_list agm2_disco_strange_No_ELM(UseJack), agm2_disco_charm_No_ELM(UseJack);
  //improved
  distr_t_list agm2_disco_impr_strange_No_ELM(UseJack), agm2_disco_impr_charm_No_ELM(UseJack);
  

  //light
  distr_t_list agm2_light(UseJack), agm2_light_fit(UseJack), agm2_light_2L_fit(UseJack), agm2_light_Lprime_fit(UseJack)  , agm2_light_infL_fit(UseJack);
  distr_t_list agm2_light_OS(UseJack), agm2_light_fit_OS(UseJack), agm2_light_2L_fit_OS(UseJack), agm2_light_Lprime_fit_OS(UseJack), agm2_light_infL_fit_OS(UseJack);


 
  //disco light
  distr_t_list agm2_disco_light(UseJack), agm2_disco_light_ELM(UseJack);
  //improved
  distr_t_list agm2_disco_impr_light(UseJack), agm2_disco_impr_light_ELM(UseJack);
  //improved deflated-stochastic
  distr_t_list agm2_disco_impr_lightD(UseJack), agm2_disco_impr_lightD_ELM(UseJack);
  //improved deflated-deflated
  distr_t_list amg2_disco_impr_lightDD(UseJack), agm2_disco_impr_lightDD_ELM(UseJack);


  //off-diagonal
  //improved
  distr_t_list agm2_disco_impr_light_charm_No_ELM(UseJack), agm2_disco_impr_light_strange_No_ELM(UseJack), agm2_disco_impr_strange_charm_No_ELM(UseJack);
  //improved D
  distr_t_list agm2_disco_impr_lightD_charm_No_ELM(UseJack), agm2_disco_impr_lightD_strange_No_ELM(UseJack);

  //#############################################################################################################
  //window contributions
  //light
  distr_t_list agm2_light_W(UseJack), agm2_light_SD(UseJack), agm2_light_W_ELM(UseJack), agm2_light_SD_ELM(UseJack), agm2_light_W_2a(UseJack), agm2_light_W_der_err_scale_setting(UseJack), agm2_light_SD_der_err_scale_setting(UseJack), agm2_light_tot_der_err_scale_setting(UseJack);
  vector<distr_t_list> agm2_light_SD_tmins_distr_list(tmins.size());
  Vfloat pert_result_SD_list;
  distr_t_list agm2_light_W_OS(UseJack), agm2_light_SD_OS(UseJack), agm2_light_W_ELM_OS(UseJack), agm2_light_SD_ELM_OS(UseJack);
  vector<distr_t_list> agm2_light_SD_OS_tmins_distr_list(tmins.size());
  //light bounding amu_W
  distr_t_list amu_W_bounding_tm_list(UseJack);
  distr_t_list amu_W_bounding_OS_list(UseJack);
  //pion mass correction
  distr_t_list amu_W_mass_corr_tm_list(UseJack);
  distr_t_list amu_W_mass_corr_OS_list(UseJack);
  distr_t_list amu_SD_mass_corr_tm_list(UseJack);
  distr_t_list amu_SD_mass_corr_OS_list(UseJack);
  vector<string> amu_W_val_Ens_tag;
  distr_t_list amu_W_mass_corr_sea_tm_list(UseJack);
  distr_t_list amu_W_mass_corr_sea_OS_list(UseJack);
  distr_t_list amu_SD_mass_corr_sea_tm_list(UseJack);
  distr_t_list amu_SD_mass_corr_sea_OS_list(UseJack);
  vector<string> amu_W_sea_Ens_tag;

  //disco light
  distr_t_list agm2_disco_light_W(UseJack), agm2_disco_light_SD(UseJack), agm2_disco_light_W_ELM(UseJack), agm2_disco_light_SD_ELM(UseJack);
  //improved
  distr_t_list agm2_disco_impr_light_W(UseJack), agm2_disco_impr_light_SD(UseJack), agm2_disco_impr_light_W_ELM(UseJack), agm2_disco_impr_light_SD_ELM(UseJack);
  //improved deflated-stochastic
  distr_t_list agm2_disco_impr_lightD_W(UseJack), agm2_disco_impr_lightD_SD(UseJack), agm2_disco_impr_lightD_W_ELM(UseJack), agm2_disco_impr_lightD_SD_ELM(UseJack);
  //improved deflated-deflated
  distr_t_list agm2_disco_impr_lightDD_W(UseJack), agm2_disco_impr_lightDD_SD(UseJack), agm2_disco_impr_lightDD_W_ELM(UseJack), agm2_disco_impr_lightDD_SD_ELM(UseJack);
  
  
  //strange
  //L
  distr_t_list agm2_strange_W_L(UseJack), agm2_strange_SD_L(UseJack), agm2_strange_W_ELM_L(UseJack), agm2_strange_SD_ELM_L(UseJack);
  distr_t_list agm2_strange_W_OS_L(UseJack), agm2_strange_SD_OS_L(UseJack), agm2_strange_W_ELM_OS_L(UseJack), agm2_strange_SD_ELM_OS_L(UseJack);
  //M
  distr_t_list agm2_strange_W_M(UseJack), agm2_strange_SD_M(UseJack), agm2_strange_W_ELM_M(UseJack), agm2_strange_SD_ELM_M(UseJack);
  distr_t_list agm2_strange_W_OS_M(UseJack), agm2_strange_SD_OS_M(UseJack), agm2_strange_W_ELM_OS_M(UseJack), agm2_strange_SD_ELM_OS_M(UseJack);
  //Extr
  distr_t_list agm2_strange_W_Extr(UseJack), agm2_strange_SD_Extr(UseJack), agm2_strange_W_ELM_Extr(UseJack), agm2_strange_SD_ELM_Extr(UseJack);
  distr_t_list agm2_strange_W_OS_Extr(UseJack), agm2_strange_SD_OS_Extr(UseJack), agm2_strange_W_ELM_OS_Extr(UseJack), agm2_strange_SD_ELM_OS_Extr(UseJack);

  //disco strange
  distr_t_list agm2_disco_strange_W(UseJack), agm2_disco_strange_SD(UseJack); // agm2_disco_strange_W_ELM(UseJack), agm2_disco_strange_SD_ELM(UseJack);
  //improved
  distr_t_list agm2_disco_impr_strange_W(UseJack), agm2_disco_impr_strange_SD(UseJack);
  
  //charm
  //L
  distr_t_list agm2_charm_W_L(UseJack), agm2_charm_SD_L(UseJack), agm2_charm_W_ELM_L(UseJack), agm2_charm_SD_ELM_L(UseJack);
  distr_t_list agm2_charm_W_OS_L(UseJack), agm2_charm_SD_OS_L(UseJack), agm2_charm_W_ELM_OS_L(UseJack), agm2_charm_SD_ELM_OS_L(UseJack);
  //M
  distr_t_list agm2_charm_W_M(UseJack), agm2_charm_SD_M(UseJack), agm2_charm_W_ELM_M(UseJack), agm2_charm_SD_ELM_M(UseJack);
  distr_t_list agm2_charm_W_OS_M(UseJack), agm2_charm_SD_OS_M(UseJack), agm2_charm_W_ELM_OS_M(UseJack), agm2_charm_SD_ELM_OS_M(UseJack);
  //H
  distr_t_list agm2_charm_W_H(UseJack), agm2_charm_SD_H(UseJack), agm2_charm_W_ELM_H(UseJack), agm2_charm_SD_ELM_H(UseJack);
  distr_t_list agm2_charm_W_OS_H(UseJack), agm2_charm_SD_OS_H(UseJack), agm2_charm_W_ELM_OS_H(UseJack), agm2_charm_SD_ELM_OS_H(UseJack);
  //Extr
  distr_t_list agm2_charm_W_Extr(UseJack), agm2_charm_SD_Extr(UseJack), agm2_charm_W_ELM_Extr(UseJack), agm2_charm_SD_ELM_Extr(UseJack);
  distr_t_list agm2_charm_W_OS_Extr(UseJack), agm2_charm_SD_OS_Extr(UseJack), agm2_charm_W_ELM_OS_Extr(UseJack), agm2_charm_SD_ELM_OS_Extr(UseJack);


  //disco charm
  distr_t_list agm2_disco_charm_W(UseJack), agm2_disco_charm_SD(UseJack); // agm2_disco_charm_W_ELM(UseJack), agm2_disco_charm_SD_ELM(UseJack);
  //improved
  distr_t_list agm2_disco_impr_charm_W(UseJack), agm2_disco_impr_charm_SD(UseJack);


  //disco off-diagonal (light-charm)
  distr_t_list agm2_disco_impr_light_charm_W(UseJack), agm2_disco_impr_light_charm_SD(UseJack);
  //disco off-diagonal (lightD-charm)
  distr_t_list agm2_disco_impr_lightD_charm_W(UseJack), agm2_disco_impr_lightD_charm_SD(UseJack);
  //disco off-diagonal (light-strange)
  distr_t_list agm2_disco_impr_light_strange_W(UseJack), agm2_disco_impr_light_strange_SD(UseJack);
  //disco off-diagonal (lightD-strange)
  distr_t_list agm2_disco_impr_lightD_strange_W(UseJack), agm2_disco_impr_lightD_strange_SD(UseJack);
  //disco off-diagonal (strange-charm)
  distr_t_list agm2_disco_impr_strange_charm_W(UseJack), agm2_disco_impr_strange_charm_SD(UseJack);
  //#############################################################################################################

 
  
  vector<vector<fit_par>> par_list_anal_repr;//only used for light quark contribution with opposite r 
  vector<vector<fit_par>> par_list_anal_repr_OS; //only used for light quark contribution with same r

  //light MV and ZV
  distr_t_list MV_fit_light(UseJack), MV_fit_light_OS(UseJack), ZV_fit_light(UseJack), ZV_fit_light_OS(UseJack);

  //strange and charm MV and ZV
  //L
  distr_t_list MV_fit_strange_L(UseJack), MV_fit_charm_L(UseJack), MV_fit_strange_OS_L(UseJack), MV_fit_charm_OS_L(UseJack);
  distr_t_list ZV_fit_strange_L(UseJack), ZV_fit_charm_L(UseJack), ZV_fit_strange_OS_L(UseJack), ZV_fit_charm_OS_L(UseJack);
  //M
  distr_t_list MV_fit_strange_M(UseJack), MV_fit_charm_M(UseJack), MV_fit_strange_OS_M(UseJack), MV_fit_charm_OS_M(UseJack);
  distr_t_list ZV_fit_strange_M(UseJack), ZV_fit_charm_M(UseJack), ZV_fit_strange_OS_M(UseJack), ZV_fit_charm_OS_M(UseJack);
  //H
  distr_t_list  MV_fit_charm_H(UseJack),  MV_fit_charm_OS_H(UseJack);
  distr_t_list  ZV_fit_charm_H(UseJack),  ZV_fit_charm_OS_H(UseJack);

    
  distr_t_list Mpi_fit_charm(UseJack), Mpi_OS_fit_charm(UseJack), fp_fit_charm(UseJack);
  distr_t_list Mpi_fit_strange(UseJack); 
  distr_t_list Mpi_fit(UseJack), Mpi_OS_fit(UseJack),  fp_fit(UseJack), Mpi_fit_disco(UseJack), Mpi_OS_fit_disco(UseJack), fp_fit_disco(UseJack);
  distr_t_list Zv_fit(UseJack), Za_fit(UseJack), Zp_ov_Zs_fit(UseJack), Zv_fit_disco(UseJack);
  distr_t_list Zv_fit_strange(UseJack), Za_fit_strange(UseJack);
  distr_t_list Zv_fit_strange_heavy(UseJack), Za_fit_strange_heavy(UseJack);
  distr_t_list Zv_fit_strange_Extr(UseJack), Za_fit_strange_Extr(UseJack);
  distr_t_list Zv_diff_strange(UseJack), Za_diff_strange(UseJack), Zv_diff_RIMOM_strange(UseJack), Za_diff_RIMOM_strange(UseJack);
  distr_t_list Zv_fit_silvano_strange_A_ens(UseJack), Za_fit_silvano_strange_A_ens(UseJack), Zp_ov_Zs_fit_silvano_strange_A_ens(UseJack);
  distr_t_list ms_extr_list(UseJack), ms_extr_etas_list(UseJack), ms_extr_phi_list(UseJack);
  distr_t_list mc_extr_etac_list(UseJack), mc_extr_Jpsi_list(UseJack);
  distr_t_list ms_extr_diff_list(UseJack);
  distr_t_list Zv_fit_charm_L(UseJack), Za_fit_charm_L(UseJack);
  distr_t_list Zv_fit_charm_M(UseJack), Za_fit_charm_M(UseJack);
  distr_t_list Zv_fit_charm_H(UseJack), Za_fit_charm_H(UseJack);
  distr_t_list Zv_fit_charm_Extr(UseJack), Za_fit_charm_Extr(UseJack);
  distr_t_list Zv_fit_charm_light(UseJack), Za_fit_charm_light(UseJack);
  distr_t_list Zv_diff_charm(UseJack), Za_diff_charm(UseJack), Zv_diff_strange_charm(UseJack), Za_diff_strange_charm(UseJack);
  distr_t_list Zv_fit_silvano_light_A_ens(UseJack), Za_fit_silvano_light_A_ens(UseJack);
  distr_t_list Zv_fit_silvano_charm_A_ens(UseJack), Za_fit_silvano_charm_A_ens(UseJack);
  distr_t_list Mpi_fit_silvano_A_ens(UseJack);
  vector<string> A_ens_charm_silvano_tags;
  vector<string> A_ens_strange_silvano_tags;
  distr_t_list mc_extr_list(UseJack);
  distr_t_list Zv_RIMOM(UseJack), Za_RIMOM(UseJack), Zv_WI(UseJack), Za_WI(UseJack);
  Vfloat L_list, a_list, ml_list, L_list_disco, a_list_disco, ml_list_disco;
  Vfloat a_list_disco_impr, a_list_disco_impr_D, a_list_disco_impr_DD;
  distr_t_list a_distr_list(UseJack), a_distr_list_strange(UseJack), a_distr_list_charm(UseJack), a_distr_list_disco_light(UseJack);
  Vfloat L_strange_list, a_strange_list, ml_strange_list;
  Vfloat L_charm_list, a_charm_list, ml_charm_list;


  //####################### PI(Q^2) analysis ################################//

  //objects where to store Pi(Q^2) data
  vector<distr_t_list> PI_Q2_light_tm(Qs2.size()), PI_Q2_light_OS(Qs2.size());
  vector<distr_t_list> PI_Q2_strange_tm(Qs2.size()), PI_Q2_strange_OS(Qs2.size());
  vector<distr_t_list> PI_Q2_charm_tm(Qs2.size()), PI_Q2_charm_OS(Qs2.size());
  vector<distr_t_list> PI_Q2_disco(Qs2.size());
  vector<distr_t_list> PI_Q2_light_tm_pert_sub(Qs2.size()), PI_Q2_light_OS_pert_sub(Qs2.size());
  vector<distr_t_list> PI_Q2_strange_tm_pert_sub(Qs2.size()), PI_Q2_strange_OS_pert_sub(Qs2.size());
  vector<distr_t_list> PI_Q2_charm_tm_pert_sub(Qs2.size()), PI_Q2_charm_OS_pert_sub(Qs2.size());


  vector<distr_t_list> CORR_DISCO_FOR_PI_Q2(3);
  vector<distr_t_list> CORR_LIGHT_FOR_PI_Q2(3);
  vector<distr_t_list> CORR_STRANGE_FOR_PI_Q2(3);
  vector<distr_t_list> CORR_CHARM_FOR_PI_Q2(3);

  vector<distr_t_list> STRANGE_PI_Q2_FOR_DISCO(3);
  vector<distr_t_list> CHARM_PI_Q2_FOR_DISCO(3);
  vector<distr_t_list> LIGHT_PI_Q2_FOR_DISCO(3);

  

  distr_t_list a_disc_PI_Q2(UseJack,3);
  vector<string> Ens_list_disc_PI_Q2(3);
  Vfloat L_list_disco_PI_Q2(3);
  distr_t_list Mpi_fit_disc_PI_Q2(UseJack,3);
  distr_t_list fpi_fit_disc_PI_Q2(UseJack,3);
  vector<distr_t> Mrho_from_GS(3);
  //####################### PI(Q^2) analysis ################################//


  //store pseudoscalar masses for physical point ensembles
  distr_t_list Mp_light_tm(UseJack,4), Mp_s1_tm(UseJack,4), Mp_s2_tm(UseJack,4), Mp_c1_tm(UseJack,4), Mp_c2_tm(UseJack,4), Mp_c3_tm(UseJack,4);
  distr_t_list Mp_light_OS(UseJack,4), Mp_s1_OS(UseJack,4), Mp_s2_OS(UseJack,4), Mp_c1_OS(UseJack,4), Mp_c2_OS(UseJack,4), Mp_c3_OS(UseJack,4); 
  vector<string> ens_masses_id(4);
  
  //define lambda for convolution with kernel
  auto K = [&](double Mv, double t, double size) -> double { return kernel_K(t, Mv);};
  auto exp_MV = [&](double Mv, double t, double size) -> double { return exp(-Mv*t);};


































  
  
  
  //strange
  channel="s";
  for(int i_ens=0;i_ens<Nens_strange;i_ens++) {
    
    CorrAnalysis Corr(UseJack, Njacks,Nboots);
    Corr.Nt = V_strange_1_L.nrows[i_ens];


    //resample lattice spacing
    distr_t a_distr(UseJack), Za(UseJack), Zv(UseJack), Za_RIMOM_distr(UseJack), Zv_RIMOM_distr(UseJack), Za_WI_distr(UseJack), Zv_WI_distr(UseJack);
    distr_t Zv_heavy(UseJack), Za_heavy(UseJack);
    LatticeInfo L_info;
    L_info.LatInfo_new_ens(V_strange_1_L.Tag[i_ens]);
    //generate jackknife sample of input parameters
  
    if(UseJack)  { for(int ijack=0;ijack<Njacks;ijack++) {
	Za_RIMOM_distr.distr.push_back(  L_info.Za + GM()*L_info.Za_err*(1.0/sqrt(Njacks-1.0)));
	Zv_RIMOM_distr.distr.push_back(  L_info.Zv + GM()*L_info.Zv_err*(1.0/sqrt(Njacks-1.0)));
	Za_WI_distr.distr.push_back(  L_info.Za_WI + GM()*L_info.Za_WI_err*(1.0/sqrt(Njacks-1.0)));
	Zv_WI_distr.distr.push_back(  L_info.Zv_WI + GM()*L_info.Zv_WI_err*(1.0/sqrt(Njacks-1.0)));
      }
    }
    else {
      for (int iboot=0; iboot<Nboots;iboot++) {
	Za_RIMOM_distr.distr.push_back(  L_info.Za + GM()*L_info.Za_err);
	Zv_RIMOM_distr.distr.push_back(  L_info.Zv + GM()*L_info.Zv_err);
	Za_WI_distr.distr.push_back(  L_info.Za_WI + GM()*L_info.Za_WI_err);
	Zv_WI_distr.distr.push_back(  L_info.Zv_WI + GM()*L_info.Zv_WI_err);
      }
    }
  
   

    //push_back Zv and Za RI-MOM
    Zv_RIMOM.distr_list.push_back(Zv_RIMOM_distr);
    Za_RIMOM.distr_list.push_back(Za_RIMOM_distr);
    Zv_WI.distr_list.push_back(Zv_WI_distr);
    Za_WI.distr_list.push_back(Za_WI_distr);
    Zv = Zv_WI_distr;
    Za = Za_WI_distr;
    Zv_heavy = Zv_WI_distr;
    Za_heavy = Za_WI_distr;

    int Tmin_P5P5;
    int Tmax_P5P5;
    int Tmin_VV;
    int Tmax_VV;
    int Tmin_VV_OS;
    int Tmax_VV_OS;


    //set time intervals for pseudoscalar obs
    if(V_strange_1_L.Tag[i_ens].substr(1,1) == "C") {
      a_distr = a_C;
      if(V_strange_1_L.Tag[i_ens]=="cC211a.06.80") { Corr.Tmin=40; Corr.Tmax=70;}
      else crash("Cannot find ensemble tag: "+V_strange_1_L.Tag[i_ens]);
    }
    else if(V_strange_1_L.Tag[i_ens].substr(1,1) == "B") {
      a_distr = a_B;
      if(V_strange_1_L.Tag[i_ens]== "cB211a.14.64") {Corr.Tmin=31; Corr.Tmax=58;}
      else if(V_strange_1_L.Tag[i_ens] == "cB211a.25.48") {Corr.Tmin=23;Corr.Tmax=44;}
      else if(V_strange_1_L.Tag[i_ens] == "cB211b.072.64") {Corr.Tmin=36; Corr.Tmax= 57;}
      else if(V_strange_1_L.Tag[i_ens] == "cB211b.072.96") {Corr.Tmin=40; Corr.Tmax= 80;}
      else crash("Cannot find ensemble tag: "+V_strange_1_L.Tag[i_ens]);
    }
    else if(V_strange_1_L.Tag[i_ens].substr(1,1) == "A") {
      a_distr = a_A;
      if(V_strange_1_L.Tag[i_ens] == "cA211a.12.48") {Corr.Tmin=19; Corr.Tmax=33;}
      else if(V_strange_1_L.Tag[i_ens] == "cA211a.40.24") {Corr.Tmin=18; Corr.Tmax=23;}
      else if(V_strange_1_L.Tag[i_ens] == "cA211a.53.24") {Corr.Tmin=16; Corr.Tmax=22;}
      else if(V_strange_1_L.Tag[i_ens] == "cA211ab.30.32") {Corr.Tmin=23; Corr.Tmax=30;}
      else crash("Cannot find ensemble tag: "+V_strange_1_L.Tag[i_ens]);
    }
    else if(V_strange_1_L.Tag[i_ens].substr(1,1) == "D") {
      a_distr = a_D;
      if(V_strange_1_L.Tag[i_ens] == "cD211a.054.96") {Corr.Tmin=55; Corr.Tmax=88;}
      else crash("Cannot find ensemble tag: "+V_strange_1_L.Tag[i_ens]);
    }
    else crash("Ensemble tag not valid");

    //push_back lattice info
    L_strange_list.push_back(L_info.L);
    a_strange_list.push_back(a_distr.ave()/fm_to_inv_Gev);
    ml_strange_list.push_back(L_info.ml);

    a_distr_list_strange.distr_list.push_back(a_distr);

    //set Tmin_P5P5 and Tmax_P5P5 to the values Corr.Tmin and Corr.Tmax
    Tmin_P5P5 = Corr.Tmin;
    Tmax_P5P5 = Corr.Tmax;
  
    //set time intervals for vector obs
    if(V_strange_1_L.Tag[i_ens].substr(1,1) == "C") {
      if(V_strange_1_L.Tag[i_ens]=="cC211a.06.80") { Tmin_VV=26; Tmax_VV=35; Tmin_VV_OS=24; Tmax_VV_OS= 33;}
      else crash("Cannot find ensemble tag: "+V_strange_1_L.Tag[i_ens]);
    }
    else if(V_strange_1_L.Tag[i_ens].substr(1,1) == "B") {
      if(V_strange_1_L.Tag[i_ens]== "cB211a.14.64") {Tmin_VV=18; Tmax_VV=28; Tmin_VV_OS=Tmin_VV; Tmax_VV_OS=Tmax_VV;}
      else if(V_strange_1_L.Tag[i_ens] == "cB211a.25.48") {Tmin_VV=19;Tmax_VV=25;  Tmin_VV_OS=Tmin_VV; Tmax_VV_OS=Tmax_VV;}
      else if(V_strange_1_L.Tag[i_ens] == "cB211b.072.64") {Tmin_VV=25; Tmax_VV= 34; Tmin_VV_OS= 26; Tmax_VV_OS=34;} //25:32, 25:32
      else if(V_strange_1_L.Tag[i_ens] == "cB211b.072.96") {Tmin_VV=23; Tmax_VV= 31; Tmin_VV_OS=25; Tmax_VV_OS=34;} //25:32 , 25:32
      else crash("Cannot find ensemble tag: "+V_strange_1_L.Tag[i_ens]);
    }
    else if(V_strange_1_L.Tag[i_ens].substr(1,1) == "A") {
      if(V_strange_1_L.Tag[i_ens] == "cA211a.12.48") {Tmin_VV=18; Tmax_VV=28;  Tmin_VV_OS=Tmin_VV; Tmax_VV_OS=Tmax_VV;}
      else if(V_strange_1_L.Tag[i_ens] == "cA211a.40.24") {Tmin_VV=13; Tmax_VV=20;  Tmin_VV_OS=Tmin_VV; Tmax_VV_OS=Tmax_VV;}
      else if(V_strange_1_L.Tag[i_ens] == "cA211a.53.24") {Tmin_VV=15; Tmax_VV=20;  Tmin_VV_OS=Tmin_VV; Tmax_VV_OS=Tmax_VV;}
      else if(V_strange_1_L.Tag[i_ens] == "cA211ab.30.32") {Tmin_VV=17;Tmax_VV=23;  Tmin_VV_OS=Tmin_VV; Tmax_VV_OS=Tmax_VV;}
      else crash("Cannot find ensemble tag: "+V_strange_1_L.Tag[i_ens]);
    }
    else if(V_strange_1_L.Tag[i_ens].substr(1,1) == "D") {
      if(V_strange_1_L.Tag[i_ens]=="cD211a.054.96") { Tmin_VV=32; Tmax_VV=40;  Tmin_VV_OS=36; Tmax_VV_OS=45;}
      else crash("Cannot find ensemble tag: "+V_strange_1_L.Tag[i_ens]);
    }
    else crash("Ensemble tag not valid");

  

  
  
    //tm
    //L
    distr_t_list  V_strange_1_distr_L, V_strange_2_distr_L, V_strange_3_distr_L;
    distr_t_list  V_strange_distr_L, MV_strange_distr_L, ZV_strange_distr_L;
    distr_t_list M_etas_distr, M_etas_distr_heavy;
    distr_t MV_strange_L, ZV_strange_L, M_etas, M_etas_heavy;
    //M
    distr_t_list  V_strange_1_distr_M, V_strange_2_distr_M, V_strange_3_distr_M;
    distr_t_list  V_strange_distr_M, MV_strange_distr_M, ZV_strange_distr_M;
    distr_t MV_strange_M, ZV_strange_M;

    //OS
    //L
    distr_t_list  V_strange_OS_1_distr_L, V_strange_OS_2_distr_L, V_strange_OS_3_distr_L;
    distr_t_list  V_strange_OS_distr_L, MV_strange_OS_distr_L, ZV_strange_OS_distr_L;
    distr_t_list M_etas_OS_distr, M_etas_OS_distr_heavy;
    distr_t MV_strange_OS_L, ZV_strange_OS_L, M_etas_OS, M_etas_OS_heavy;
    //M
    distr_t_list  V_strange_OS_1_distr_M, V_strange_OS_2_distr_M, V_strange_OS_3_distr_M;
    distr_t_list  V_strange_OS_distr_M, MV_strange_OS_distr_M, ZV_strange_OS_distr_M;
    distr_t MV_strange_OS_M, ZV_strange_OS_M;

    //To compute Za and Zv (hadronic method)
    distr_t_list overlap_P5P5_distr, fetas_distr, overlap_P5P5_distr_heavy, fetas_distr_heavy;
    distr_t_list etas_corr, P5P5_OS_distr;
    distr_t_list etas_corr_heavy, P5P5_OS_distr_heavy;
    distr_t_list overlap_P5P5_OS_distr, overlap_P5P5_OS_distr_heavy;
    distr_t_list ratio_P5P5_overlap_OS_tm,  Zp_ov_Zs_distr, ratio_P5P5_overlap_OS_tm_heavy, Zp_ov_Zs_distr_heavy;
    distr_t_list A0P5_distr, A0P5_OS_distr, A0P5_distr_heavy, A0P5_OS_distr_heavy;
    distr_t_list RV, RA, RV_heavy, RA_heavy;
    distr_t Zv_hadr, Za_hadr, Zp_ov_Zs;
    distr_t Zv_hadr_heavy, Za_hadr_heavy, Zp_ov_Zs_heavy;
    distr_t fetas, fetas_heavy;

    //disco
    distr_t_list disco_distr;

    //disco improved
    distr_t_list disco_impr_distr;

    //off-diagonal light-strange
    distr_t_list disco_impr_light_strange_distr;
    distr_t_list disco_impr_lightD_strange_distr;

    //off-diagonal strange-charm
    distr_t_list disco_impr_strange_charm_distr;


    
    //tm
  
    etas_corr= Corr.corr_t(corr_P5P5_strange.col(0)[i_ens], "");
    M_etas_distr= Corr.effective_mass_t(corr_P5P5_strange.col(0)[i_ens], "../data/gm2/strange/tm_"+Extrapolation_strange_mode+"/P5P5_ss_mass_"+V_strange_1_L.Tag[i_ens]+".dat");
    etas_corr_heavy = Corr.corr_t(corr_P5P5_strange_heavy.col(0)[i_ens], "");
    M_etas_distr_heavy = Corr.effective_mass_t(corr_P5P5_strange_heavy.col(0)[i_ens], "../data/gm2/strange/tm_"+Extrapolation_strange_mode+"/P5P5_ss_mass_heavy_"+V_strange_1_L.Tag[i_ens]+".dat");
  
  
    overlap_P5P5_distr = Corr.residue_t(corr_P5P5_strange.col(0)[i_ens], "");
    overlap_P5P5_distr_heavy = Corr.residue_t(corr_P5P5_strange_heavy.col(0)[i_ens], "");
    fetas_distr = 2.0*L_info.ms_L*Corr.decay_constant_t(corr_P5P5_strange.col(0)[i_ens], "");
    fetas_distr_heavy= 2.0*L_info.ms_M*Corr.decay_constant_t(corr_P5P5_strange_heavy.col(0)[i_ens], "");
    fetas = Corr.Fit_distr(fetas_distr);
    fetas_heavy = Corr.Fit_distr(fetas_distr_heavy);
 
    //L
    V_strange_1_distr_L = Corr.corr_t(V_strange_1_L.col(0)[i_ens], "../data/gm2/strange/tm_"+Extrapolation_strange_mode+"/corr_1_"+V_strange_1_L.Tag[i_ens]+"_L.dat");
    V_strange_2_distr_L = Corr.corr_t(V_strange_2_L.col(0)[i_ens], "../data/gm2/strange/tm_"+Extrapolation_strange_mode+"/corr_2_"+V_strange_2_L.Tag[i_ens]+"_L.dat");
    V_strange_3_distr_L = Corr.corr_t(V_strange_3_L.col(0)[i_ens], "../data/gm2/strange/tm_"+Extrapolation_strange_mode+"/corr_3_"+V_strange_3_L.Tag[i_ens]+"_L.dat");
 
  
    //M
    V_strange_1_distr_M = Corr.corr_t(V_strange_1_M.col(0)[i_ens], "../data/gm2/strange/tm_"+Extrapolation_strange_mode+"/corr_1_"+V_strange_1_M.Tag[i_ens]+"_M.dat");
    V_strange_2_distr_M = Corr.corr_t(V_strange_2_M.col(0)[i_ens], "../data/gm2/strange/tm_"+Extrapolation_strange_mode+"/corr_2_"+V_strange_2_M.Tag[i_ens]+"_M.dat");
    V_strange_3_distr_M = Corr.corr_t(V_strange_3_M.col(0)[i_ens], "../data/gm2/strange/tm_"+Extrapolation_strange_mode+"/corr_3_"+V_strange_3_M.Tag[i_ens]+"_M.dat");
 
 
    //OS
 
    P5P5_OS_distr = Corr.corr_t(corr_P5P5_OS_strange.col(0)[i_ens], "");
    P5P5_OS_distr_heavy = Corr.corr_t(corr_P5P5_OS_strange_heavy.col(0)[i_ens], "");
    M_etas_OS_distr= Corr.effective_mass_t(corr_P5P5_OS_strange.col(0)[i_ens], "../data/gm2/strange/OS_"+Extrapolation_strange_mode+"/P5P5_ss_mass_"+V_strange_1_L.Tag[i_ens]+".dat");
    M_etas_OS_distr_heavy = Corr.effective_mass_t(corr_P5P5_OS_strange_heavy.col(0)[i_ens], "../data/gm2/strange/OS_"+Extrapolation_strange_mode+"/P5P5_ss_mass_heavy_"+V_strange_1_L.Tag[i_ens]+".dat");

 
    overlap_P5P5_OS_distr= Corr.residue_t(corr_P5P5_OS_strange.col(0)[i_ens], "");
    overlap_P5P5_OS_distr_heavy = Corr.residue_t(corr_P5P5_OS_strange_heavy.col(0)[i_ens], "");

    //L
    V_strange_OS_1_distr_L = Corr.corr_t(V_strange_OS_1_L.col(0)[i_ens], "../data/gm2/strange/OS_"+Extrapolation_strange_mode+"/corr_1_"+V_strange_1_L.Tag[i_ens]+"_L.dat");
    V_strange_OS_2_distr_L = Corr.corr_t(V_strange_OS_2_L.col(0)[i_ens], "../data/gm2/strange/OS_"+Extrapolation_strange_mode+"/corr_2_"+V_strange_2_L.Tag[i_ens]+"_L.dat");
    V_strange_OS_3_distr_L = Corr.corr_t(V_strange_OS_3_L.col(0)[i_ens], "../data/gm2/strange/OS_"+Extrapolation_strange_mode+"/corr_3_"+V_strange_3_L.Tag[i_ens]+"_L.dat");
    //M
    V_strange_OS_1_distr_M = Corr.corr_t(V_strange_OS_1_M.col(0)[i_ens], "../data/gm2/strange/OS_"+Extrapolation_strange_mode+"/corr_1_"+V_strange_1_M.Tag[i_ens]+"_M.dat");
    V_strange_OS_2_distr_M = Corr.corr_t(V_strange_OS_2_M.col(0)[i_ens], "../data/gm2/strange/OS_"+Extrapolation_strange_mode+"/corr_2_"+V_strange_2_M.Tag[i_ens]+"_M.dat");
    V_strange_OS_3_distr_M = Corr.corr_t(V_strange_OS_3_M.col(0)[i_ens], "../data/gm2/strange/OS_"+Extrapolation_strange_mode+"/corr_3_"+V_strange_3_M.Tag[i_ens]+"_M.dat");



 

   
    //#######################################  COMPUTATION OF ZV AND ZA (Hadronic method) ##################################
    //define lambda functions to be used
    auto sqr= [=](double a, double b) {return sqrt(a);};
    auto SINH= [](double m) -> double  {return sinh(m);};

  
    //take ratio between OS and tm pion amplitude to compute Zp/Zs RC.
    ratio_P5P5_overlap_OS_tm= overlap_P5P5_OS_distr/overlap_P5P5_distr;
    Zp_ov_Zs_distr = distr_t_list::f_of_distr_list(sqr, ratio_P5P5_overlap_OS_tm);
    ratio_P5P5_overlap_OS_tm_heavy= overlap_P5P5_OS_distr_heavy/overlap_P5P5_distr_heavy;
    Zp_ov_Zs_distr_heavy = distr_t_list::f_of_distr_list(sqr, ratio_P5P5_overlap_OS_tm_heavy);
  

  
    //antysymmetrize w.r.t. t -> T-t for A0P5 correlators
    if(V_strange_1_L.Tag[i_ens].substr(1,1) == "A") Corr.Reflection_sign = -1;
    A0P5_distr= Corr.corr_t(corr_A0P5_strange.col(0)[i_ens], "");
    A0P5_OS_distr = Corr.corr_t(corr_A0P5_OS_strange.col(0)[i_ens], "");
    A0P5_distr_heavy = Corr.corr_t(corr_A0P5_strange_heavy.col(0)[i_ens], "");
    A0P5_OS_distr_heavy = Corr.corr_t(corr_A0P5_OS_strange_heavy.col(0)[i_ens],"");

 


  
    //restore symmetrization
    Corr.Reflection_sign = 1;

    //compute RV (estimator for Zv)


    RV= 2.0*L_info.ms_L*etas_corr/distr_t_list::derivative(A0P5_distr, 0); //central derivative
    RV_heavy = 2.0*L_info.ms_M*etas_corr_heavy/distr_t_list::derivative(A0P5_distr_heavy, 0); //central derivative

   
    //tm and OS P5P5
    M_etas = Corr.Fit_distr(M_etas_distr);
    M_etas_heavy = Corr.Fit_distr(M_etas_distr_heavy);

    cout<<"M_etas (L)  ["<<V_strange_1_L.Tag[i_ens]<<" ]  : "<<M_etas.ave()<<" "<<M_etas.err()<<" "<<Corr.Tmin<<" "<<Corr.Tmax<<endl;
    cout<<"M_etas (H)  ["<<V_strange_1_L.Tag[i_ens]<<" ] : "<<M_etas_heavy.ave()<<" "<<M_etas_heavy.err()<<" "<<Corr.Tmin<<" "<<Corr.Tmax<<endl;
    M_etas_OS= Corr.Fit_distr(M_etas_OS_distr);
    M_etas_OS_heavy = Corr.Fit_distr(M_etas_OS_distr_heavy);

    int id=0;
    if(V_strange_1_L.Tag[i_ens].substr(1,1) != "A") {
      if(V_strange_1_L.Tag[i_ens] == "cB211b.072.64") id=0;
      else if(V_strange_1_L.Tag[i_ens] == "cB211b.072.96") id=1;
      else if(V_strange_1_L.Tag[i_ens].substr(1,1) == "C") id=2;
      else if(V_strange_1_L.Tag[i_ens].substr(1,1) == "D") id=3;
      else crash("Ensemble not recognized");
      ens_masses_id[id] == V_strange_1_L.Tag[i_ens];
      Mp_s1_tm.distr_list[id] = M_etas;
      Mp_s2_tm.distr_list[id] = M_etas_heavy;
      Mp_s1_OS.distr_list[id] = M_etas_OS;
      Mp_s2_OS.distr_list[id] = M_etas_OS_heavy;
    }
    //fit obs to compute Zv and Za (hadronic method)
    Zp_ov_Zs = Corr.Fit_distr(Zp_ov_Zs_distr);
    Zp_ov_Zs_heavy = Corr.Fit_distr(Zp_ov_Zs_distr_heavy);
  
    //set plateaux for RV
    if(V_strange_1_L.Tag[i_ens].substr(1,1) == "A") {
      if( V_strange_1_L.Tag[i_ens] == "cA211a.40.24") {Corr.Tmin=20; Corr.Tmax = 23;}
      else if( V_strange_1_L.Tag[i_ens] == "cA211a.53.24") { Corr.Tmin=20; Corr.Tmax= 23;}
      else if( V_strange_1_L.Tag[i_ens] == "cA211ab.30.32") { Corr.Tmin =23; Corr.Tmax=31;}
      else crash("Ensemble A: "+V_strange_1_L.Tag[i_ens]+" not found");
    }
    else if(V_strange_1_L.Tag[i_ens].substr(1,1) == "C") { Corr.Tmin=36; Corr.Tmax=74;}
    else if(V_strange_1_L.Tag[i_ens].substr(1,1) == "B") {
      if(V_strange_1_L.Tag[i_ens] == "cB211b.072.64") {    Corr.Tmin=32; Corr.Tmax = 60 ;}
      else if(V_strange_1_L.Tag[i_ens] == "cB211b.072.96") { Corr.Tmin=35; Corr.Tmax= 80;    }
      else crash("Ensemble: "+V_strange_1_L.Tag[i_ens]+" not found");
    }
    else if(V_strange_1_L.Tag[i_ens].substr(1,1) == "D") { Corr.Tmin=45; Corr.Tmax=91;}
    else crash("Ensemble: "+V_strange_1_L.Tag[i_ens]+" not found");
    Zv_hadr= Corr.Fit_distr(RV);
    Zv_hadr_heavy = Corr.Fit_distr(RV_heavy);
    RA = 2.0*L_info.ms_L*(P5P5_OS_distr/distr_t_list::derivative(A0P5_OS_distr, 0))*(M_etas_OS/M_etas)*(distr_t::f_of_distr(SINH, M_etas_OS)/distr_t::f_of_distr(SINH, M_etas))*(1.0/Zp_ov_Zs);
    RA_heavy = 2.0*L_info.ms_M*(P5P5_OS_distr_heavy/distr_t_list::derivative(A0P5_OS_distr_heavy, 0))*(M_etas_OS_heavy/M_etas_heavy)*(distr_t::f_of_distr(SINH, M_etas_OS_heavy)/distr_t::f_of_distr(SINH, M_etas_heavy))*(1.0/Zp_ov_Zs_heavy);
    //set plateaux for RA
    int Tmin_RA=0;
    int Tmax_RA=0;
    //set time intervals for RA
    if(V_strange_1_L.Tag[i_ens].substr(1,1) == "C") {
      if(V_strange_1_L.Tag[i_ens]=="cC211a.06.80") { Tmin_RA=40; Tmax_RA=65;}
      else crash("Cannot find ensemble tag: "+V_strange_1_L.Tag[i_ens]);
    }
    else if(V_strange_1_L.Tag[i_ens].substr(1,1) == "B") {
      if(V_strange_1_L.Tag[i_ens]== "cB211a.14.64") {Tmin_RA=18; Tmax_RA=26;}
      else if(V_strange_1_L.Tag[i_ens] == "cB211a.25.48") {Tmin_RA=10;Tmax_RA=21;}
      else if(V_strange_1_L.Tag[i_ens] == "cB211b.072.64") {Tmin_RA=30; Tmax_RA= 57;}
      else if(V_strange_1_L.Tag[i_ens] == "cB211b.072.96") {Tmin_RA=30; Tmax_RA= 80;}
      else crash("Cannot find ensemble tag: "+V_strange_1_L.Tag[i_ens]);
    }
    else if(V_strange_1_L.Tag[i_ens].substr(1,1) == "A") {
      if(V_strange_1_L.Tag[i_ens] == "cA211a.12.48") {Tmin_RA=15; Tmax_RA=25;}
      else if(V_strange_1_L.Tag[i_ens] == "cA211a.40.24") {Tmin_RA=11; Tmax_RA=22;}
      else if(V_strange_1_L.Tag[i_ens] == "cA211a.53.24") {Tmin_RA=10; Tmax_RA=21;}
      else if(V_strange_1_L.Tag[i_ens] == "cA211ab.30.32") {Tmin_RA=17;Tmax_RA=30;}
      else crash("Cannot find ensemble tag: "+V_strange_1_L.Tag[i_ens]);
    }
    else if(V_strange_1_L.Tag[i_ens].substr(1,1) == "D") {
      if(V_strange_1_L.Tag[i_ens] == "cD211a.054.96") {Tmin_RA=51; Tmax_RA=88;}
      else crash("Cannot find ensemble tag: "+V_strange_1_L.Tag[i_ens]);
    }
    else crash("Ensemble tag not valid");

    Corr.Tmin=Tmin_RA;
    Corr.Tmax=Tmax_RA;
  
    Za_hadr= Corr.Fit_distr(RA);
    Za_hadr_heavy = Corr.Fit_distr(RA_heavy);

    //print Rv and RA
    //print RV
    Print_To_File({}, {RV.ave(), RV.err(), RV_heavy.ave(), RV_heavy.err()}, "../data/gm2/strange/Z_"+Extrapolation_strange_mode+"/RV_ss"+V_strange_1_L.Tag[i_ens]+".dat.t", "", "");
    //print RA
    Print_To_File({}, {RA.ave(), RA.err(), RA_heavy.ave(), RA_heavy.err()}, "../data/gm2/strange/Z_"+Extrapolation_strange_mode+"/RA_ss"+V_strange_1_L.Tag[i_ens]+".dat.t", "", "");
    //print Zp_ov_Zs
    Print_To_File({}, {Zp_ov_Zs_distr.ave(), Zp_ov_Zs_distr.err(), Zp_ov_Zs_distr_heavy.ave(), Zp_ov_Zs_distr_heavy.err()}, "../data/gm2/strange/Z_"+Extrapolation_strange_mode+"/Zp_ov_Zs_ss"+V_strange_1_L.Tag[i_ens]+".dat.t", "", "");

    //push_back

    Zv_fit_strange.distr_list.push_back(Zv_hadr);
    Za_fit_strange.distr_list.push_back(Za_hadr);
    Zv_fit_strange_heavy.distr_list.push_back(Zv_hadr_heavy);
    Za_fit_strange_heavy.distr_list.push_back(Za_hadr_heavy);

    if(Use_Za_Zv_from_strange_run) { Zv= Zv_hadr; Za=Za_hadr; Zv_heavy = Zv_hadr_heavy; Za_heavy= Za_hadr_heavy;}

    //################################################ END OF COMPUTATION OF ZV AND ZA (Hadronic method) #####################################à



 
    //sum over the Lorenz indices of the e.m. current
    //L
    //tm
    V_strange_distr_L= (pow(qs,2)/3.0)*(V_strange_1_distr_L+ V_strange_2_distr_L + V_strange_3_distr_L);
    //OS
    V_strange_OS_distr_L= (pow(qs,2)/3.0)*(V_strange_OS_1_distr_L+ V_strange_OS_2_distr_L + V_strange_OS_3_distr_L);
  
    //M
    //tm
    V_strange_distr_M= (pow(qs,2)/3.0)*(V_strange_1_distr_M+ V_strange_2_distr_M + V_strange_3_distr_M);
    //OS
    V_strange_OS_distr_M= (pow(qs,2)/3.0)*(V_strange_OS_1_distr_M+ V_strange_OS_2_distr_M + V_strange_OS_3_distr_M);


    //extract effective masses, overlap from V and fit

    //tm
    //L
    Corr.Tmin=Tmin_VV;
    Corr.Tmax=Tmax_VV;
    MV_strange_distr_L= Corr.effective_mass_t(V_strange_distr_L, "../data/gm2/strange/tm_"+Extrapolation_strange_mode+"/MV_mass_"+V_strange_1_L.Tag[i_ens]+"_L.dat");
    ZV_strange_distr_L= Corr.residue_t(V_strange_distr_L, "../data/gm2/strange/tm_"+Extrapolation_strange_mode+"/ZV_overlap_"+V_strange_1_L.Tag[i_ens]+"_L.dat");
    MV_strange_L = Corr.Fit_distr(MV_strange_distr_L);
    ZV_strange_L = Corr.Fit_distr(ZV_strange_distr_L);

    cout<<"M_phi (L) ["<<V_strange_1_L.Tag[i_ens]<<" ] : "<<MV_strange_L.ave()<<" "<<MV_strange_L.err()<<" "<<Corr.Tmin<<" "<<Corr.Tmax<<endl;
    
  
    //M
    Corr.Tmin=Tmin_VV;
    Corr.Tmax=Tmax_VV;
    MV_strange_distr_M= Corr.effective_mass_t(V_strange_distr_M, "../data/gm2/strange/tm_"+Extrapolation_strange_mode+"/MV_mass_"+V_strange_1_M.Tag[i_ens]+"_M.dat");
    ZV_strange_distr_M= Corr.residue_t(V_strange_distr_M, "../data/gm2/strange/tm_"+Extrapolation_strange_mode+"/ZV_overlap_"+V_strange_1_M.Tag[i_ens]+"_M.dat");
    MV_strange_M = Corr.Fit_distr(MV_strange_distr_M);
    ZV_strange_M = Corr.Fit_distr(ZV_strange_distr_M);

    cout<<"M_phi (H)  ["<<V_strange_1_L.Tag[i_ens]<<" ] : "<<MV_strange_M.ave()<<" "<<MV_strange_M.err()<<" "<<Corr.Tmin<<" "<<Corr.Tmax<<endl;


    //Print to File M_etas and Mphi

    distr_t_list etas_masses_list(UseJack);
    distr_t_list phi_masses_list(UseJack);

    etas_masses_list.distr_list.push_back( M_etas);
    etas_masses_list.distr_list.push_back( M_etas_heavy);
    phi_masses_list.distr_list.push_back( MV_strange_L);
    phi_masses_list.distr_list.push_back( MV_strange_M);

    Print_To_File({}, {etas_masses_list.ave(), etas_masses_list.err(), phi_masses_list.ave(), phi_masses_list.err()}, "../data/gm2/strange/masses_"+V_strange_1_L.Tag[i_ens]+".list", "", "# etas   phi");
 
   
    //OS
    //L
    Corr.Tmin= Tmin_VV_OS;
    Corr.Tmax= Tmax_VV_OS;
    MV_strange_OS_distr_L= Corr.effective_mass_t(V_strange_OS_distr_L, "../data/gm2/strange/OS_"+Extrapolation_strange_mode+"/MV_mass_"+V_strange_1_L.Tag[i_ens]+"_L.dat");
    ZV_strange_OS_distr_L= Corr.residue_t(V_strange_OS_distr_L, "../data/gm2/strange/OS_"+Extrapolation_strange_mode+"/ZV_overlap_"+V_strange_1_L.Tag[i_ens]+"_L.dat");
    MV_strange_OS_L = Corr.Fit_distr(MV_strange_OS_distr_L);
    ZV_strange_OS_L = Corr.Fit_distr(ZV_strange_OS_distr_L);
    //M
    MV_strange_OS_distr_M= Corr.effective_mass_t(V_strange_OS_distr_M, "../data/gm2/strange/OS_"+Extrapolation_strange_mode+"/MV_mass_"+V_strange_1_M.Tag[i_ens]+"_M.dat");
    ZV_strange_OS_distr_M= Corr.residue_t(V_strange_OS_distr_M, "../data/gm2/strange/OS_"+Extrapolation_strange_mode+"/ZV_overlap_"+V_strange_1_M.Tag[i_ens]+"_M.dat");
    MV_strange_OS_M = Corr.Fit_distr(MV_strange_OS_distr_M);
    ZV_strange_OS_M = Corr.Fit_distr(ZV_strange_OS_distr_M);




    bool Found_disco_ens=false;
    bool Found_disco_impr_ens=false;
    if(Include_strange_disco) {
      int i_ens_disco=0;
      
      for(int j=0;j<Nens_disco_strange;j++) if(disco_strange.Tag[j] == V_strange_1_L.Tag[i_ens]) { Found_disco_ens=true; i_ens_disco=j;disco_strange_Tags.push_back(disco_strange.Tag[j]) ;break;}
      if(Found_disco_ens) {
	disco_distr = Corr.corr_t(disco_strange.col(0)[i_ens_disco], "");
	disco_distr = disco_distr*(pow(qs,2));
      }

      //improved
      for(int j=0;j<disco_impr_strange.size;j++) if( disco_impr_strange.Tag[j] == V_strange_1_L.Tag[i_ens]) {Found_disco_impr_ens=true; i_ens_disco=j;disco_impr_strange_Tags.push_back(disco_impr_strange.Tag[j]) ;break;}
      if(Found_disco_impr_ens) {
	disco_impr_distr = Corr.corr_t(disco_impr_strange.col(0)[i_ens_disco], "");
	disco_impr_distr = disco_impr_distr*(pow(qs,2));
      }
      
	  
    }

    bool Found_disco_impr_light_strange_ens=false;
    bool Found_disco_impr_lightD_strange_ens=false;
    bool Found_disco_impr_strange_charm_ens=false;
    if(Include_off_diagonal_disco) {

      int i_ens_disco=0;

      //improved light-strange
      for(int j=0;j<disco_impr_light_strange.size;j++) if( disco_impr_light_strange.Tag[j] == V_strange_1_L.Tag[i_ens]) {Found_disco_impr_light_strange_ens=true; i_ens_disco=j;disco_impr_light_strange_Tags.push_back(disco_impr_light_strange.Tag[j]) ;break;}
      if(Found_disco_impr_light_strange_ens) {
	disco_impr_light_strange_distr = Corr.corr_t(disco_impr_light_strange.col(0)[i_ens_disco], "");
	disco_impr_light_strange_distr = disco_impr_light_strange_distr*( 2.0*qs*( qu+qd)  );
      }

      //improved lightD-strange
      for(int j=0;j<disco_impr_lightD_strange.size;j++) if( disco_impr_lightD_strange.Tag[j] == V_strange_1_L.Tag[i_ens]) {Found_disco_impr_lightD_strange_ens=true; i_ens_disco=j;disco_impr_lightD_strange_Tags.push_back(disco_impr_lightD_strange.Tag[j]) ;break;}
      if(Found_disco_impr_lightD_strange_ens) {
	disco_impr_lightD_strange_distr = Corr.corr_t(disco_impr_lightD_strange.col(0)[i_ens_disco], "");
	disco_impr_lightD_strange_distr = disco_impr_lightD_strange_distr*(2.0*qs*(qu+qd));
      }


      //improved strange-charm
      for(int j=0;j<disco_impr_strange_charm.size;j++) if( disco_impr_strange_charm.Tag[j] == V_strange_1_L.Tag[i_ens]) {Found_disco_impr_strange_charm_ens=true; i_ens_disco=j;disco_impr_strange_charm_Tags.push_back(disco_impr_strange_charm.Tag[j]) ;break;}
      if(Found_disco_impr_strange_charm_ens) {
	disco_impr_strange_charm_distr = Corr.corr_t(disco_impr_strange_charm.col(0)[i_ens_disco], "");
	disco_impr_strange_charm_distr = disco_impr_strange_charm_distr*(2.0*qs*qc);
      }
      
    }



  
    
    //###########   EXTRAPOLATE m_s phys  Za_phys , Zv_phys  #############

 
   
    distr_t m_etas_phys_distr, m_phi_phys_distr;
  
    for(int ijack=0;ijack<Njacks;ijack++) m_etas_phys_distr.distr.push_back( m_etas + GM()*m_etas_err/sqrt(Njacks-1.0));
    for(int ijack=0;ijack<Njacks;ijack++) m_phi_phys_distr.distr.push_back( m_phi+ GM()*m_phi_err/sqrt(Njacks-1.0));

    vector<distr_t> X_2_fit;
    vector<distr_t> X_4_fit;
    distr_t X_2_phys;
    distr_t X_4_phys;
    vector<distr_t> M2phi_fit({MV_strange_L/a_distr, MV_strange_M/a_distr});
    vector<distr_t> M2etas_fit({M_etas/a_distr, M_etas_heavy/a_distr});

    if(Extrapolation_strange_mode == "phi") {
      X_2_fit =  {MV_strange_L/a_distr, MV_strange_M/a_distr};
      X_4_fit =  {MV_strange_L*MV_strange_L/(a_distr*a_distr), MV_strange_M*MV_strange_M/(a_distr*a_distr)};
      X_2_phys = m_phi_phys_distr;
      X_4_phys = m_phi_phys_distr*m_phi_phys_distr;
    }
    else if(Extrapolation_strange_mode == "etas") {
      X_2_fit =  { M_etas/a_distr, M_etas_heavy/a_distr};
      X_4_fit =  { M_etas*M_etas/(a_distr*a_distr), M_etas_heavy*M_etas_heavy/(a_distr*a_distr)};
      X_2_phys = m_etas_phys_distr;
      X_4_phys = m_etas_phys_distr*m_etas_phys_distr;
    }
    else crash("Extrapolation strange mode: "+Extrapolation_strange_mode+" not yet implemented");

  
    vector<distr_t> Za_hadr_list, Zv_hadr_list;
    vector<distr_t> Zp_ov_Zs_hadr_list;
    Za_hadr_list = {Za_hadr, Za_hadr_heavy};
    Zv_hadr_list = {Zv_hadr, Zv_hadr_heavy};
    Zp_ov_Zs_hadr_list = {Zp_ov_Zs, Zp_ov_Zs_heavy};

  
    
    //Generate fake ms_distr
    distr_t ms_light_distr;
    distr_t ms_heavy_distr;
    for(int ijack=0;ijack<Njacks;ijack++) ms_light_distr.distr.push_back( L_info.ms_L );
    for(int ijack=0;ijack<Njacks;ijack++) ms_heavy_distr.distr.push_back( L_info.ms_M );
  
    //estrapolate ms phys
    vector<distr_t> ms_list( {ms_light_distr, ms_heavy_distr});
    distr_t ms_phys_etas_extr = Obs_extrapolation_meson_mass(ms_list, M2etas_fit, m_etas_phys_distr,  "../data/gm2/strange", "ms_extrapolation_etas_"+V_strange_1_L.Tag[i_ens], UseJack, "SPLINE");
    distr_t ms_phys_phi_extr = Obs_extrapolation_meson_mass(ms_list, M2phi_fit, m_phi_phys_distr,  "../data/gm2/strange", "ms_extrapolation_phi_"+V_strange_1_L.Tag[i_ens], UseJack, "SPLINE");
    //push_back
    distr_t ms_phys_extr=(Extrapolation_strange_mode=="etas")?ms_phys_etas_extr:ms_phys_phi_extr;
    ms_extr_list.distr_list.push_back(ms_phys_extr/a_distr);
    ms_extr_etas_list.distr_list.push_back(ms_phys_etas_extr/a_distr);
    ms_extr_phi_list.distr_list.push_back( ms_phys_phi_extr/a_distr);
    ms_extr_diff_list.distr_list.push_back( (ms_phys_etas_extr-ms_phys_phi_extr)/a_distr);


    //Extrapolate f_etas using ms_phys_extr
    vector<distr_t> fetas_list({ fetas, fetas_heavy});
    distr_t fetas_extr = Obs_extrapolation_meson_mass( fetas_list, ms_list, ms_phys_extr, "../data/gm2/strange", "f_etas_extr_quark_mass_"+Extrapolation_strange_mode+"_"+V_strange_1_L.Tag[i_ens], UseJack, "SPLINE");
    

    //Extrapolate Zv and Za using ms_phys_extr
    distr_t Za_hadr_extr = Obs_extrapolation_meson_mass( Za_hadr_list, ms_list, ms_phys_extr, "../data/gm2/strange", "Za_extr_quark_mass_"+Extrapolation_strange_mode+"_"+V_strange_1_L.Tag[i_ens], UseJack, "SPLINE");
    distr_t Zv_hadr_extr = Obs_extrapolation_meson_mass( Zv_hadr_list, ms_list, ms_phys_extr, "../data/gm2/strange", "Zv_extr_quark_mass_"+Extrapolation_strange_mode+"_"+V_strange_1_L.Tag[i_ens], UseJack, "SPLINE");

    distr_t Zp_ov_Zs_hadr_extr = Obs_extrapolation_meson_mass( Zp_ov_Zs_hadr_list, ms_list, ms_phys_extr, "../data/gm2/strange", "Zp_ov_Zs_extr_quark_mass_"+Extrapolation_strange_mode+"_"+V_strange_1_L.Tag[i_ens], UseJack, "SPLINE");

    //push_back the result
    Zv_fit_strange_Extr.distr_list.push_back(Zv_hadr_extr);
    Za_fit_strange_Extr.distr_list.push_back(Za_hadr_extr);
    Zv_diff_strange.distr_list.push_back(Zv_hadr_extr- Zv_WI_distr);
    Za_diff_strange.distr_list.push_back(Za_hadr_extr-Za_WI_distr);
    Zv_diff_RIMOM_strange.distr_list.push_back(Zv_hadr_extr -Zv_RIMOM_distr);
    Za_diff_RIMOM_strange.distr_list.push_back(Za_hadr_extr -Za_RIMOM_distr);


    cout<<"Zp_ov_Zs ["<<V_strange_1_L.Tag[i_ens]<<"] : "<<Zp_ov_Zs_hadr_extr.ave()<<" +- "<<Zp_ov_Zs_hadr_extr.err()<<endl;
    cout<<"Zv["<<V_strange_1_L.Tag[i_ens]<<"] : "<<Zv_hadr_extr.ave()<<" +- "<<Zv_hadr_extr.err()<<endl;
    cout<<"Za["<<V_strange_1_L.Tag[i_ens]<<"] : "<<Za_hadr_extr.ave()<<" +- "<<Za_hadr_extr.err()<<endl;
    cout<<"Za_ov_Zv["<<V_strange_1_L.Tag[i_ens]<<"] : "<<(Za_hadr_extr/Zv_hadr_extr).ave()<<" +- "<<(Za_hadr_extr/Zv_hadr_extr).err()<<endl;

 

    if(Use_Za_Zv_from_strange_run && Use_Extrapolated_Za_Zv_strange) { Za = Za_hadr_extr; Za_heavy= Za_hadr_extr; Zv= Zv_hadr_extr; Zv_heavy= Zv_hadr_extr;}

    //push_back Za and Zv hadr extr if ensemble A

    if(V_strange_1_L.Tag[i_ens].substr(1,1) == "A") {

      Zv_fit_silvano_strange_A_ens.distr_list.push_back( Zv_hadr_extr);
      Za_fit_silvano_strange_A_ens.distr_list.push_back( Za_hadr_extr);
      Zp_ov_Zs_fit_silvano_strange_A_ens.distr_list.push_back( Zp_ov_Zs_hadr_extr);
      A_ens_strange_silvano_tags.push_back( V_strange_1_L.Tag[i_ens]);

    }


    //free corr LO artifacts
    //################## READ FREE THEORY VECTOR-VECTOR CORRELATOR OPPOSITE R ####################################
    string Pt_free_oppor_L= "../Vkvk_cont/"+to_string(Corr.Nt/2)+"_m"+to_string_with_precision(L_info.ms_L,3)+"/OPPOR";
    string Pt_free_oppor_M= "../Vkvk_cont/"+to_string(Corr.Nt/2)+"_m"+to_string_with_precision(L_info.ms_M,3)+"/OPPOR";
    Vfloat VV_free_oppor_L= Read_From_File(Pt_free_oppor_L, 1, 4);
    Vfloat VV_free_oppor_M= Read_From_File(Pt_free_oppor_M, 1, 4);
    if(VV_free_oppor_L.size() != Corr.Nt) crash("Failed to read properly free VV correlator mc_L w opposite r");
    if(VV_free_oppor_M.size() != Corr.Nt) crash("Failed to read properly free VV correlator mc_M w opposite r");
    //################## READ FREE THEORY VECTOR-VECTOR CORRELATOR SAME R ####################################
    string Pt_free_samer_L= "../Vkvk_cont/"+to_string(Corr.Nt/2)+"_m"+to_string_with_precision(L_info.ms_L,3)+"/SAMER";
    string Pt_free_samer_M= "../Vkvk_cont/"+to_string(Corr.Nt/2)+"_m"+to_string_with_precision(L_info.ms_M,3)+"/SAMER";
    Vfloat VV_free_samer_L= Read_From_File(Pt_free_samer_L, 1, 4);
    Vfloat VV_free_samer_M= Read_From_File(Pt_free_samer_M, 1, 4);
    if(VV_free_samer_L.size() != Corr.Nt) crash("Failed to read properly free VV correlator mc_L  w same r");
    if(VV_free_samer_M.size() != Corr.Nt) crash("Failed to read properly free VV correlator mc_M  w same r");
    //Insert electric charges
    for( auto & OP:VV_free_oppor_L) OP *= pert_corr_strange_on_off*qs*qs;
    for( auto & SA:VV_free_samer_L) SA *= pert_corr_strange_on_off*qs*qs;
    for( auto & OP:VV_free_oppor_M) OP *= pert_corr_strange_on_off*qs*qs;
    for( auto & SA:VV_free_samer_M) SA *= pert_corr_strange_on_off*qs*qs;
 
 


    //free corr LO artifacts
    Vfloat free_corr_log_art(Corr.Nt, 0.0);
    for(int t=0;t<Corr.Nt;t++) { if( t*a_distr.ave() < add_pert_corr_strange_up_to*fm_to_inv_Gev && t != 0) {   free_corr_log_art[t] = -1.0*pert_corr_strange_on_off*(qs*qs)*(1.0/(2.0*M_PI*M_PI*pow(t,5)));} else free_corr_log_art[t] = 0.0;

       if(t==0 || t*a_distr.ave() > add_pert_corr_strange_up_to*fm_to_inv_Gev) { VV_free_samer_L[t] =0; VV_free_samer_M[t] = 0;  VV_free_oppor_L[t] = 0; VV_free_oppor_M[t]=0;}


    }


    /*
    VV_free_oppor_L = free_corr_log_art;
    VV_free_oppor_M = free_corr_log_art;
    VV_free_samer_L = free_corr_log_art;
    VV_free_samer_M = free_corr_log_art; */ 
    

    distr_t_list V_strange_distr_L_pert_sub, V_strange_OS_distr_L_pert_sub;
    distr_t_list V_strange_distr_M_pert_sub, V_strange_OS_distr_M_pert_sub;

    if(!sum_pert_corr_strange_to_bare_corr) { //sum to renormalized correlator
      //L
      V_strange_distr_L_pert_sub = (1.0/(Za*Za))*(Za*Za*V_strange_distr_L + VV_free_oppor_L);
      V_strange_OS_distr_L_pert_sub = (1.0/(Zv*Zv))*(Zv*Zv*V_strange_OS_distr_L + VV_free_samer_L);
      //M
      V_strange_distr_M_pert_sub = (1.0/(Za_heavy*Za_heavy))*(Za_heavy*Za_heavy*V_strange_distr_M + VV_free_oppor_M);
      V_strange_OS_distr_M_pert_sub =  (1.0/(Zv_heavy*Zv_heavy))*(Zv_heavy*Zv_heavy*V_strange_OS_distr_M + VV_free_samer_M);
    }
    else { //sum to bare correlator
      //L
      V_strange_distr_L_pert_sub = (V_strange_distr_L + VV_free_oppor_L);
      V_strange_OS_distr_L_pert_sub = (V_strange_OS_distr_L + VV_free_samer_L);
      //M
      V_strange_distr_M_pert_sub = (V_strange_distr_M + VV_free_oppor_M);
      V_strange_OS_distr_M_pert_sub =  (V_strange_OS_distr_M + VV_free_samer_M);
    }
  
    // print summed correlators to file
    //L
    //tm
    Print_To_File({}, {V_strange_distr_L.ave(), V_strange_distr_L.err(), (Za*Za*V_strange_distr_L).ave(), (Za*Za*V_strange_distr_L).err()}, "../data/gm2/strange/tm_"+Extrapolation_strange_mode+"/corr_sum_"+V_strange_1_L.Tag[i_ens]+"_L.dat.t", "", "# t    bare   renormalized");
    //OS
    Print_To_File({}, {V_strange_OS_distr_L.ave(), V_strange_OS_distr_L.err(), (Zv*Zv*V_strange_OS_distr_L).ave(), (Zv*Zv*V_strange_OS_distr_L).err()}, "../data/gm2/strange/OS_"+Extrapolation_strange_mode+"/corr_sum_"+V_strange_1_L.Tag[i_ens]+"_L.dat.t", "", "#t   bare   renormalized");
    //tm
    Print_To_File({}, {(Za*Za*V_strange_distr_L_pert_sub).ave(), (Za*Za*V_strange_distr_L_pert_sub).err()}, "../data/gm2/strange/tm_"+Extrapolation_strange_mode+"/pert_sub_corr_"+V_strange_1_L.Tag[i_ens]+"_L.dat.t", "", "# t    bare   renormalized");
    //OS
    Print_To_File({}, {(Zv*Zv*V_strange_OS_distr_L_pert_sub).ave(), (Zv*Zv*V_strange_OS_distr_L_pert_sub).err()}, "../data/gm2/strange/OS_"+Extrapolation_strange_mode+"/pert_sub_corr_"+V_strange_1_L.Tag[i_ens]+"_L.dat.t", "", "#t   bare   renormalized");

    //M
    //tm
    Print_To_File({}, {V_strange_distr_M.ave(), V_strange_distr_M.err(), (Za_heavy*Za_heavy*V_strange_distr_M).ave(), (Za_heavy*Za_heavy*V_strange_distr_M).err()}, "../data/gm2/strange/tm_"+Extrapolation_strange_mode+"/corr_sum_"+V_strange_1_M.Tag[i_ens]+"_M.dat.t", "", "# t    bare   renormalized");
    //OS
    Print_To_File({}, {V_strange_OS_distr_M.ave(), V_strange_OS_distr_M.err(), (Zv_heavy*Zv_heavy*V_strange_OS_distr_M).ave(), (Zv_heavy*Zv_heavy*V_strange_OS_distr_M).err()}, "../data/gm2/strange/OS_"+Extrapolation_strange_mode+"/corr_sum_"+V_strange_1_M.Tag[i_ens]+"_M.dat.t", "", "#t   bare   renormalized");
    //tm
    Print_To_File({}, {(Za_heavy*Za_heavy*V_strange_distr_M_pert_sub).ave(), (Za_heavy*Za_heavy*V_strange_distr_M_pert_sub).err()}, "../data/gm2/strange/tm_"+Extrapolation_strange_mode+"/pert_sub_corr_"+V_strange_1_L.Tag[i_ens]+"_M.dat.t", "", "# t    bare   renormalized");
    //OS
    Print_To_File({}, {(Zv_heavy*Zv_heavy*V_strange_OS_distr_M_pert_sub).ave(), (Zv_heavy*Zv_heavy*V_strange_OS_distr_M_pert_sub).err()}, "../data/gm2/strange/OS_"+Extrapolation_strange_mode+"/pert_sub_corr_"+V_strange_1_L.Tag[i_ens]+"_M.dat.t", "", "#t   bare   renormalized");

    //print diff L-M
    Print_To_File({},{ (V_strange_distr_L- V_strange_distr_M).ave(), (V_strange_distr_L-V_strange_distr_M).err(), V_strange_distr_L.ave(), V_strange_distr_L.err(), V_strange_distr_M.ave(), V_strange_distr_M.err()}, "../data/gm2/strange/tm_"+Extrapolation_strange_mode+"/corr_sum_"+V_strange_1_L.Tag[i_ens]+"_diff_LM.dat.t","", "#dC(t) err_dC(t)   C(t)_L  dC(t)_L, C(t)_M, dC(t)_M");
    Print_To_File({},{ (V_strange_OS_distr_L- V_strange_OS_distr_M).ave(), (V_strange_OS_distr_L-V_strange_OS_distr_M).err(), V_strange_OS_distr_L.ave(), V_strange_OS_distr_L.err(), V_strange_OS_distr_M.ave(), V_strange_OS_distr_M.err()}, "../data/gm2/strange/OS_"+Extrapolation_strange_mode+"/corr_sum_"+V_strange_1_L.Tag[i_ens]+"_diff_LM.dat.t","", "#dC(t) err_dC(t)   C(t)_L  dC(t)_L, C(t)_M, dC(t)_M");


    //print disco strange
    if(Include_strange_disco && Found_disco_ens) Print_To_File({}, {disco_distr.ave(), disco_distr.err(), (Zv*Zv*disco_distr).ave(), (Zv*Zv*disco_distr).err()}, "../data/gm2/strange/disco/disc_"+V_strange_1_L.Tag[i_ens]+".dat.t","","# bare renormalized");

    //print disco improved
    if(Include_strange_disco && Found_disco_impr_ens) Print_To_File({}, {disco_impr_distr.ave(), disco_impr_distr.err(), (Zv*Zv*disco_impr_distr).ave(), (Zv*Zv*disco_impr_distr).err()}, "../data/gm2/strange/disco/disc_impr_"+V_strange_1_L.Tag[i_ens]+".dat.t","","# bare renormalized");


    //print off-diagonal light-strange
    if(Include_off_diagonal_disco && Found_disco_impr_light_strange_ens) Print_To_File({}, {disco_impr_light_strange_distr.ave(), disco_impr_light_strange_distr.err(), (Zv*Zv*disco_impr_light_strange_distr).ave(), (Zv*Zv*disco_impr_light_strange_distr).err()}, "../data/gm2/light_strange/disco/disc_impr_"+V_strange_1_L.Tag[i_ens]+".dat.t","","# bare renormalized");
    //print off-diagonal lightD-strange
    if(Include_off_diagonal_disco && Found_disco_impr_lightD_strange_ens) Print_To_File({}, {disco_impr_lightD_strange_distr.ave(), disco_impr_lightD_strange_distr.err(), (Zv*Zv*disco_impr_lightD_strange_distr).ave(), (Zv*Zv*disco_impr_lightD_strange_distr).err()}, "../data/gm2/light_strange/disco/disc_impr_D_"+V_strange_1_L.Tag[i_ens]+".dat.t","","# bare renormalized");
    //print off-diagonal strange-charm
    if(Include_off_diagonal_disco && Found_disco_impr_strange_charm_ens) Print_To_File({}, {disco_impr_strange_charm_distr.ave(), disco_impr_strange_charm_distr.err(), (Zv*Zv*disco_impr_strange_charm_distr).ave(), (Zv*Zv*disco_impr_strange_charm_distr).err()}, "../data/gm2/strange_charm/disco/disc_impr_"+V_strange_1_L.Tag[i_ens]+".dat.t","","# bare renormalized");
    
    


    //compute effective mass of phi meson including disconnected diagram

    if(Found_disco_impr_ens) {

      //tm
      distr_t_list MV_L_full_eff_distr= Corr.effective_mass_t( Za*Za*V_strange_distr_L + Zv*Zv*disco_impr_distr, "../data/gm2/strange/tm_"+Extrapolation_strange_mode+"/MV_mass_w_disco_"+V_strange_1_L.Tag[i_ens]+"_L.dat");
      distr_t_list MV_M_full_eff_distr= Corr.effective_mass_t( Za*Za*V_strange_distr_M + Zv*Zv*disco_impr_distr, "../data/gm2/strange/tm_"+Extrapolation_strange_mode+"/MV_mass_w_disco_"+V_strange_1_L.Tag[i_ens]+"_M.dat");
      distr_t_list MV_L_disco_eff_distr = Corr.effective_mass_t( Zv*Zv*disco_impr_distr, "../data/gm2/strange/tm_"+Extrapolation_strange_mode+"/MV_mass_only_disco_"+V_strange_1_L.Tag[i_ens]+"_L.dat") ;
      distr_t_list MV_M_disco_eff_distr = Corr.effective_mass_t( Zv*Zv*disco_impr_distr, "../data/gm2/strange/tm_"+Extrapolation_strange_mode+"/MV_mass_only_disco_"+V_strange_1_L.Tag[i_ens]+"_M.dat");
      //OS
      distr_t_list MV_L_full_eff_OS_distr= Corr.effective_mass_t( Zv*Zv*V_strange_OS_distr_L + Zv*Zv*disco_impr_distr, "../data/gm2/strange/OS_"+Extrapolation_strange_mode+"/MV_mass_w_disco_"+V_strange_1_L.Tag[i_ens]+"_L.dat");
      distr_t_list MV_M_full_eff_OS_distr= Corr.effective_mass_t( Zv*Zv*V_strange_OS_distr_M + Zv*Zv*disco_impr_distr, "../data/gm2/strange/OS_"+Extrapolation_strange_mode+"/MV_mass_w_disco_"+V_strange_1_L.Tag[i_ens]+"_M.dat");
      distr_t_list MV_L_disco_eff_OS_distr = Corr.effective_mass_t( Zv*Zv*disco_impr_distr, "../data/gm2/strange/OS_"+Extrapolation_strange_mode+"/MV_mass_only_disco_"+V_strange_1_L.Tag[i_ens]+"_L.dat") ;
      distr_t_list MV_M_disco_eff_OS_distr = Corr.effective_mass_t( Zv*Zv*disco_impr_distr, "../data/gm2/strange/OS_"+Extrapolation_strange_mode+"/MV_mass_only_disco_"+V_strange_1_L.Tag[i_ens]+"_M.dat");
    }
  


 

    //push_back MV_strange, ZV_strange, tm
 
    //L
    MV_fit_strange_L.distr_list.push_back(MV_strange_L);
    ZV_fit_strange_L.distr_list.push_back(Za*Za*ZV_strange_L); 
 
    //M
    MV_fit_strange_M.distr_list.push_back(MV_strange_M);
    ZV_fit_strange_M.distr_list.push_back(Za_heavy*Za_heavy*ZV_strange_M);

    //push_back MV_strange and ZV_strange, OS
  
    //L
    MV_fit_strange_OS_L.distr_list.push_back(MV_strange_OS_L);
    ZV_fit_strange_OS_L.distr_list.push_back(Zv*Zv*ZV_strange_OS_L);
    //M
    MV_fit_strange_OS_M.distr_list.push_back(MV_strange_OS_M);
    ZV_fit_strange_OS_M.distr_list.push_back(Zv_heavy*Zv_heavy*ZV_strange_OS_M);

 


  
  

 


 
    int Tdata_min= 8;
    int Tdata_max = Corr.Nt/2.0 -2;
    int Tdata_fit = (Tdata_min + Tdata_max)/2;

    distr_t ELM_mass_L;
    distr_t ELM_mass_M;
    distr_t ELM_mass_OS_L;
    distr_t ELM_mass_OS_M;
  
    if(ELM_mass_strange == "phi") { ELM_mass_L = MV_strange_L/m_phi_phys_distr; ELM_mass_M = MV_strange_M/m_phi_phys_distr; ELM_mass_OS_L= MV_strange_OS_L/m_phi_phys_distr; ELM_mass_OS_M= MV_strange_OS_M/m_phi_phys_distr;}
    else if(ELM_mass_strange == "etas") {ELM_mass_L = M_etas/m_etas_phys_distr; ELM_mass_M = M_etas_heavy/m_etas_phys_distr; ELM_mass_OS_L = M_etas_OS/m_etas_phys_distr; ELM_mass_OS_M = M_etas_OS_heavy/m_etas_phys_distr;}
    else crash("ELM mass strange: "+ELM_mass_strange+" not yet implemented");
  
    //compute kernel distribution
    distr_t_list Kernel_distr_No_ELM= distr_t_list::f_of_distr(K,a_distr, Upper_Limit_Time_Integral_strange+1);
    //tm
    //L
    distr_t_list Kernel_distr_L = distr_t_list::f_of_distr(K,ELM_mass_L, Upper_Limit_Time_Integral_strange+1);
    //M
    distr_t_list Kernel_distr_M = distr_t_list::f_of_distr(K,ELM_mass_M, Upper_Limit_Time_Integral_strange+1);
 
  
    //OS
    //L
    distr_t_list Kernel_OS_distr_L = distr_t_list::f_of_distr(K, ELM_mass_OS_L, Upper_Limit_Time_Integral_strange +1);
    //M
    distr_t_list Kernel_OS_distr_M = distr_t_list::f_of_distr(K, ELM_mass_OS_M, Upper_Limit_Time_Integral_strange +1);
  
    //compute exp(-Mv*t) distribution
    //tm
    //L
    distr_t_list exp_MVs_L = distr_t_list::f_of_distr(exp_MV, MV_strange_L, Upper_Limit_Time_Integral_strange+1);
    //M
    distr_t_list exp_MVs_M = distr_t_list::f_of_distr(exp_MV, MV_strange_M, Upper_Limit_Time_Integral_strange+1);
 
  

    //OS
    //L
    distr_t_list exp_OS_MVs_L = distr_t_list::f_of_distr(exp_MV, MV_strange_OS_L, Upper_Limit_Time_Integral_strange+1);
    //M
    distr_t_list exp_OS_MVs_M = distr_t_list::f_of_distr(exp_MV, MV_strange_OS_M, Upper_Limit_Time_Integral_strange+1);
 
  
    //#######################################################################################################################################à
    //Print single-exponential prediction to file
    //tm
    //L
    Print_To_File({}, {(exp_MVs_L*(ZV_strange_L/(2.0*MV_strange_L))).ave(), (exp_MVs_L*(ZV_strange_L/(2.0*MV_strange_L))).err(), (Za*Za*exp_MVs_L*(ZV_strange_L/(2.0*MV_strange_L))).ave(), (Za*Za*exp_MVs_L*(ZV_strange_L/(2.0*MV_strange_L))).err() }, "../data/gm2/strange/tm_"+Extrapolation_strange_mode+"/corr_gsd_sum_"+V_strange_1_L.Tag[i_ens]+"_L.dat.t", "", "#t   bare   renormalized");
    //M
    Print_To_File({}, {(exp_MVs_M*(ZV_strange_M/(2.0*MV_strange_M))).ave(), (exp_MVs_M*(ZV_strange_M/(2.0*MV_strange_M))).err(), (Za_heavy*Za_heavy*exp_MVs_M*(ZV_strange_M/(2.0*MV_strange_M))).ave(), (Za_heavy*Za_heavy*exp_MVs_M*(ZV_strange_M/(2.0*MV_strange_M))).err() }, "../data/gm2/strange/tm_"+Extrapolation_strange_mode+"/corr_gsd_sum_"+V_strange_1_M.Tag[i_ens]+"_M.dat.t", "", "#t   bare   renormalized");

   
    //OS
    //L
    Print_To_File({}, {(exp_OS_MVs_L*(ZV_strange_OS_L/(2.0*MV_strange_OS_L))).ave(), (exp_OS_MVs_L*(ZV_strange_OS_L/(2.0*MV_strange_OS_L))).err(), (Zv*Zv*exp_OS_MVs_L*(ZV_strange_OS_L/(2.0*MV_strange_OS_L))).ave(), (Zv*Zv*exp_OS_MVs_L*(ZV_strange_OS_L/(2.0*MV_strange_OS_L))).err()}, "../data/gm2/strange/OS_"+Extrapolation_strange_mode+"/corr_gsd_sum_"+V_strange_1_L.Tag[i_ens]+"_L.dat.t", "", "#t   bare   renormalized");
    //M
    Print_To_File({}, {(exp_OS_MVs_M*(ZV_strange_OS_M/(2.0*MV_strange_OS_M))).ave(), (exp_OS_MVs_M*(ZV_strange_OS_M/(2.0*MV_strange_OS_M))).err(), (Zv_heavy*Zv_heavy*exp_OS_MVs_M*(ZV_strange_OS_M/(2.0*MV_strange_OS_M))).ave(), (Zv_heavy*Zv_heavy*exp_OS_MVs_M*(ZV_strange_OS_M/(2.0*MV_strange_OS_M))).err()}, "../data/gm2/strange/OS_"+Extrapolation_strange_mode+"/corr_gsd_sum_"+V_strange_1_M.Tag[i_ens]+"_M.dat.t", "", "#t   bare   renormalized");
    //########################################################################################################################################
  
    distr_t_list agm2_distr_Tdata_L(UseJack), agm2_OS_distr_Tdata_L(UseJack);
    distr_t_list agm2_distr_Tdata_No_ELM_L(UseJack), agm2_OS_distr_Tdata_No_ELM_L(UseJack);
    distr_t_list agm2_distr_Tdata_M(UseJack), agm2_OS_distr_Tdata_M(UseJack);
    distr_t_list agm2_distr_Tdata_No_ELM_M(UseJack), agm2_OS_distr_Tdata_No_ELM_M(UseJack);
    Vfloat Tdata_vec;
    bool Find_Tdata_fit=false;

   
    for(int Tdata=Tdata_min;Tdata<Tdata_max;Tdata++) {
      //compute 4\pia^2 using lattice data up to Tdata (included)
   
      distr_t agm2_L(UseJack, UseJack?Njacks:Nboots); //constructor sets it to 0
      distr_t agm2_No_ELM_L(UseJack, UseJack?Njacks:Nboots); //constructor sets it to 0
      distr_t agm2_OS_L(UseJack, UseJack?Njacks:Nboots); //constructor sets it to 0
      distr_t agm2_OS_No_ELM_L(UseJack, UseJack?Njacks:Nboots); //constructor sets it to 0
      distr_t agm2_M(UseJack, UseJack?Njacks:Nboots); //constructor sets it to 0
      distr_t agm2_No_ELM_M(UseJack, UseJack?Njacks:Nboots); //constructor sets it to 0
      distr_t agm2_OS_M(UseJack, UseJack?Njacks:Nboots); //constructor sets it to 0
      distr_t agm2_OS_No_ELM_M(UseJack, UseJack?Njacks:Nboots); //constructor sets it to 0
       
      for(int t=1;t<=Upper_Limit_Time_Integral_strange;t++) {
	if(t<=Tdata) {
	  //L
	  agm2_L = agm2_L + w(t,Simps_ord)*4.0*pow(alpha,2)*V_strange_distr_L_pert_sub.distr_list[t]*Kernel_distr_L.distr_list[t];
	  agm2_No_ELM_L = agm2_No_ELM_L +  w(t,Simps_ord)*4.0*pow(alpha,2)*V_strange_distr_L_pert_sub.distr_list[t]*Kernel_distr_No_ELM.distr_list[t];
	  agm2_OS_L = agm2_OS_L + w(t,Simps_ord)*4.0*pow(alpha,2)*V_strange_OS_distr_L_pert_sub.distr_list[t]*Kernel_OS_distr_L.distr_list[t];
	  agm2_OS_No_ELM_L = agm2_OS_No_ELM_L +  w(t,Simps_ord)*4.0*pow(alpha,2)*V_strange_OS_distr_L_pert_sub.distr_list[t]*Kernel_distr_No_ELM.distr_list[t];
	  //M
	  agm2_M = agm2_M +  w(t,Simps_ord)*4.0*pow(alpha,2)*V_strange_distr_M_pert_sub.distr_list[t]*Kernel_distr_M.distr_list[t];
	  agm2_No_ELM_M = agm2_No_ELM_M +  w(t,Simps_ord)*4.0*pow(alpha,2)*V_strange_distr_M_pert_sub.distr_list[t]*Kernel_distr_No_ELM.distr_list[t];
	  agm2_OS_M = agm2_OS_M +  w(t,Simps_ord)*4.0*pow(alpha,2)*V_strange_OS_distr_M_pert_sub.distr_list[t]*Kernel_OS_distr_M.distr_list[t];
	  agm2_OS_No_ELM_M = agm2_OS_No_ELM_M +  w(t,Simps_ord)*4.0*pow(alpha,2)*V_strange_OS_distr_M_pert_sub.distr_list[t]*Kernel_distr_No_ELM.distr_list[t];
	
	}
	else {
	  //L
	  agm2_L= agm2_L +  w(t,Simps_ord)*4.0*pow(alpha,2)*(ZV_strange_L/(2.0*MV_strange_L))*exp_MVs_L.distr_list[t]*Kernel_distr_L.distr_list[t];
	  agm2_No_ELM_L= agm2_No_ELM_L +  w(t,Simps_ord)*4.0*pow(alpha,2)*(ZV_strange_L/(2.0*MV_strange_L))*exp_MVs_L.distr_list[t]*Kernel_distr_No_ELM.distr_list[t];
	  agm2_OS_L= agm2_OS_L +  w(t,Simps_ord)*4.0*pow(alpha,2)*(ZV_strange_OS_L/(2.0*MV_strange_OS_L))*exp_OS_MVs_L.distr_list[t]*Kernel_OS_distr_L.distr_list[t];
	  agm2_OS_No_ELM_L= agm2_OS_No_ELM_L +  w(t,Simps_ord)*4.0*pow(alpha,2)*(ZV_strange_OS_L/(2.0*MV_strange_OS_L))*exp_OS_MVs_L.distr_list[t]*Kernel_distr_No_ELM.distr_list[t];
	  //M
	  agm2_M= agm2_M +  w(t,Simps_ord)*4.0*pow(alpha,2)*(ZV_strange_M/(2.0*MV_strange_M))*exp_MVs_M.distr_list[t]*Kernel_distr_M.distr_list[t];
	  agm2_No_ELM_M= agm2_No_ELM_M +  w(t,Simps_ord)*4.0*pow(alpha,2)*(ZV_strange_M/(2.0*MV_strange_M))*exp_MVs_M.distr_list[t]*Kernel_distr_No_ELM.distr_list[t];
	  agm2_OS_M= agm2_OS_M +  w(t,Simps_ord)*4.0*pow(alpha,2)*(ZV_strange_OS_M/(2.0*MV_strange_OS_M))*exp_OS_MVs_M.distr_list[t]*Kernel_OS_distr_M.distr_list[t];
	  agm2_OS_No_ELM_M= agm2_OS_No_ELM_M +  w(t,Simps_ord)*4.0*pow(alpha,2)*(ZV_strange_OS_M/(2.0*MV_strange_OS_M))*exp_OS_MVs_M.distr_list[t]*Kernel_distr_No_ELM.distr_list[t];

	}
      }
    
      Tdata_vec.push_back((double)Tdata);
      //L
      agm2_distr_Tdata_L.distr_list.push_back(Za*Za*agm2_L);
      agm2_OS_distr_Tdata_L.distr_list.push_back(Zv*Zv*agm2_OS_L);
      agm2_distr_Tdata_No_ELM_L.distr_list.push_back(Za*Za*agm2_No_ELM_L);
      agm2_OS_distr_Tdata_No_ELM_L.distr_list.push_back(Zv*Zv*agm2_OS_No_ELM_L);
      //M
      agm2_distr_Tdata_M.distr_list.push_back(Za_heavy*Za_heavy*agm2_M);
      agm2_OS_distr_Tdata_M.distr_list.push_back(Zv_heavy*Zv_heavy*agm2_OS_M);
      agm2_distr_Tdata_No_ELM_M.distr_list.push_back(Za_heavy*Za_heavy*agm2_No_ELM_M);
      agm2_OS_distr_Tdata_No_ELM_M.distr_list.push_back(Zv_heavy*Zv_heavy*agm2_OS_No_ELM_M);
 

    
    
      if(Tdata==Tdata_fit){
	//push back L
	agm2_strange_L.distr_list.push_back(Za*Za*agm2_L);
	agm2_strange_No_ELM_L.distr_list.push_back(Za*Za*agm2_No_ELM_L);
	agm2_strange_OS_L.distr_list.push_back(Zv*Zv*agm2_OS_L);
	agm2_strange_OS_No_ELM_L.distr_list.push_back(Zv*Zv*agm2_OS_No_ELM_L);
	//push back M
	agm2_strange_M.distr_list.push_back(Za_heavy*Za_heavy*agm2_M);
	agm2_strange_No_ELM_M.distr_list.push_back(Za_heavy*Za_heavy*agm2_No_ELM_M);
	agm2_strange_OS_M.distr_list.push_back(Zv_heavy*Zv_heavy*agm2_OS_M);
	agm2_strange_OS_No_ELM_M.distr_list.push_back(Zv_heavy*Zv_heavy*agm2_OS_No_ELM_M);
	//extrapolate to the physical kaon point
	vector<distr_t> agm2s_strange({Za*Za*agm2_L, Za_heavy*Za_heavy*agm2_M});
	vector<distr_t> agm2s_strange_OS({ Zv*Zv*agm2_OS_L, Zv_heavy*Zv_heavy*agm2_OS_M});
	vector<distr_t> agm2s_strange_No_ELM({ Za*Za*agm2_No_ELM_L, Za_heavy*Za_heavy*agm2_No_ELM_M});
	vector<distr_t> agm2s_strange_OS_No_ELM({ Zv*Zv*agm2_OS_No_ELM_L, Zv_heavy*Zv_heavy*agm2_OS_No_ELM_M});
	agm2_strange_Extr.distr_list.push_back( Obs_extrapolation_meson_mass(agm2s_strange, X_2_fit, X_2_phys,  "../data/gm2/strange", "agm2_ELM_Extrapolation_"+Extrapolation_strange_mode+"_"+V_strange_1_L.Tag[i_ens], UseJack, "SPLINE"));
	agm2_strange_OS_Extr.distr_list.push_back( Obs_extrapolation_meson_mass(agm2s_strange_OS, X_2_fit, X_2_phys,  "../data/gm2/strange", "agm2_OS_ELM_Extrapolation_"+Extrapolation_strange_mode+"_"+V_strange_1_L.Tag[i_ens], UseJack, "SPLINE"));
	agm2_strange_No_ELM_Extr.distr_list.push_back( Obs_extrapolation_meson_mass(agm2s_strange_No_ELM, X_2_fit, X_2_phys,  "../data/gm2/strange", "agm2_Extrapolation_"+Extrapolation_strange_mode+"_"+V_strange_1_L.Tag[i_ens], UseJack, "SPLINE"));
	agm2_strange_OS_No_ELM_Extr.distr_list.push_back( Obs_extrapolation_meson_mass(agm2s_strange_OS_No_ELM, X_2_fit, X_2_phys,  "../data/gm2/strange", "agm2_OS_Extrapolation_"+Extrapolation_strange_mode+"_"+V_strange_1_L.Tag[i_ens], UseJack, "SPLINE"));
      
	Find_Tdata_fit=true;
      }
    }

    if(!Find_Tdata_fit) crash("Cannot find Tdata fit value: "+to_string(Tdata_fit));


   
    //#######################  INTERMEDIATE AND SHORT-DISTANCE WINDOW ###################################

    
    //############################   TWISTED MASS ######################################################

    //L
    distr_t agm2_W_L(UseJack, UseJack?Njacks:Nboots);  //constructor sets agm2_W to zero by default
    distr_t agm2_SD_L(UseJack, UseJack?Njacks:Nboots);  //constructor sets agm2_SD to zero by default
    distr_t agm2_W_ELM_L(UseJack, UseJack?Njacks:Nboots); //constructor sets agm2_W_ELM to zero by default
    distr_t agm2_SD_ELM_L(UseJack, UseJack?Njacks:Nboots);  //constructor sets agm2_SD_ELM to zero by default
    //M
    distr_t agm2_W_M(UseJack, UseJack?Njacks:Nboots);  //constructor sets agm2_W to zero by default
    distr_t agm2_SD_M(UseJack, UseJack?Njacks:Nboots);  //constructor sets agm2_SD to zero by default
    distr_t agm2_W_ELM_M(UseJack, UseJack?Njacks:Nboots); //constructor sets agm2_W_ELM to zero by default
    distr_t agm2_SD_ELM_M(UseJack, UseJack?Njacks:Nboots);  //constructor sets agm2_SD_ELM to zero by default


    //#################################################################################################


  
    distr_t_list Ker = distr_t_list::f_of_distr(K, a_distr , Corr.Nt/2);
    //L
    distr_t_list Ker_ELM_tm_L = distr_t_list::f_of_distr(K, ELM_mass_L, Corr.Nt/2);
    distr_t_list Ker_ELM_OS_L = distr_t_list::f_of_distr(K, ELM_mass_OS_L, Corr.Nt/2);
    //M
    distr_t_list Ker_ELM_tm_M = distr_t_list::f_of_distr(K, ELM_mass_M, Corr.Nt/2);
    distr_t_list Ker_ELM_OS_M = distr_t_list::f_of_distr(K, ELM_mass_OS_M, Corr.Nt/2);

    
    //define lambdas for the theta func
    auto th0 = [](double ta) ->double { return 1.0/(1.0 + exp(-2.0*(ta-t0)/Delta));};
    auto th1 = [](double ta) ->double { return 1.0/(1.0 + exp(-2.0*(ta-t1)/Delta));};

  
    for(int t=1; t< Corr.Nt/2; t++) {
      //L
      agm2_W_L = agm2_W_L +  w(t,Simps_ord)*4.0*pow(alpha,2)*Za*Za*V_strange_distr_L.distr_list[t]*Ker.distr_list[t]*( distr_t::f_of_distr(th0, t*a_distr) - distr_t::f_of_distr(th1, t*a_distr));
      agm2_SD_L = agm2_SD_L +  w(t,Simps_ord)*4.0*pow(alpha,2)*Za*Za*(V_strange_distr_L_pert_sub.distr_list[t])*Ker.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*a_distr));
      agm2_W_ELM_L = agm2_W_ELM_L +  w(t,Simps_ord)*4.0*pow(alpha,2)*Za*Za*V_strange_distr_L.distr_list[t]*Ker_ELM_tm_L.distr_list[t]*( distr_t::f_of_distr(th0, t*ELM_mass_L) - distr_t::f_of_distr(th1, t*ELM_mass_L));
      agm2_SD_ELM_L = agm2_SD_ELM_L +  w(t,Simps_ord)*4.0*pow(alpha,2)*Za*Za*(V_strange_distr_L_pert_sub.distr_list[t])*Ker_ELM_tm_L.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*ELM_mass_L));
      //M
      agm2_W_M = agm2_W_M +  w(t,Simps_ord)*4.0*pow(alpha,2)*Za_heavy*Za_heavy*V_strange_distr_M.distr_list[t]*Ker.distr_list[t]*( distr_t::f_of_distr(th0, t*a_distr) - distr_t::f_of_distr(th1, t*a_distr));
      agm2_SD_M = agm2_SD_M +  w(t,Simps_ord)*4.0*pow(alpha,2)*Za_heavy*Za_heavy*(V_strange_distr_M_pert_sub.distr_list[t])*Ker.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*a_distr));
      agm2_W_ELM_M = agm2_W_ELM_M +  w(t,Simps_ord)*4.0*pow(alpha,2)*Za_heavy*Za_heavy*V_strange_distr_M.distr_list[t]*Ker_ELM_tm_M.distr_list[t]*( distr_t::f_of_distr(th0, t*ELM_mass_M) - distr_t::f_of_distr(th1, t*ELM_mass_M));
      agm2_SD_ELM_M = agm2_SD_ELM_M +  w(t,Simps_ord)*4.0*pow(alpha,2)*Za_heavy*Za_heavy*(V_strange_distr_M_pert_sub.distr_list[t])*Ker_ELM_tm_M.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*ELM_mass_M));

    }
  
    //push_back the result

    //L
    agm2_strange_W_L.distr_list.push_back(agm2_W_L);
    agm2_strange_SD_L.distr_list.push_back(agm2_SD_L);
    agm2_strange_W_ELM_L.distr_list.push_back(agm2_W_ELM_L);
    agm2_strange_SD_ELM_L.distr_list.push_back(agm2_SD_ELM_L);
    //M
    agm2_strange_W_M.distr_list.push_back(agm2_W_M);
    agm2_strange_SD_M.distr_list.push_back(agm2_SD_M);
    agm2_strange_W_ELM_M.distr_list.push_back(agm2_W_ELM_M);
    agm2_strange_SD_ELM_M.distr_list.push_back(agm2_SD_ELM_M);


    //extrapolate the result to the physical kaon point

    vector<distr_t> agm2s_strange_W({agm2_W_L, agm2_W_M});
    vector<distr_t> agm2s_strange_SD({agm2_SD_L, agm2_SD_M});
    vector<distr_t> agm2s_strange_W_ELM({agm2_W_ELM_L, agm2_W_ELM_M});
    vector<distr_t> agm2s_strange_SD_ELM({agm2_SD_ELM_L, agm2_SD_ELM_M});

  
    agm2_strange_W_Extr.distr_list.push_back( Obs_extrapolation_meson_mass(agm2s_strange_W, X_2_fit, X_2_phys, "../data/gm2/strange", "agm2_W_Extrapolation_"+Extrapolation_strange_mode+"_"+V_strange_1_L.Tag[i_ens], UseJack, "SPLINE"));
    agm2_strange_SD_Extr.distr_list.push_back( Obs_extrapolation_meson_mass(agm2s_strange_SD, X_2_fit, X_2_phys, "../data/gm2/strange", "agm2_SD_Extrapolation_"+Extrapolation_strange_mode+"_"+V_strange_1_L.Tag[i_ens], UseJack, "SPLINE"));
    agm2_strange_W_ELM_Extr.distr_list.push_back( Obs_extrapolation_meson_mass(agm2s_strange_W_ELM, X_2_fit,X_2_phys, "../data/gm2/strange", "agm2_W_ELM_Extrapolation_"+Extrapolation_strange_mode+"_"+V_strange_1_L.Tag[i_ens], UseJack, "SPLINE"));
    agm2_strange_SD_ELM_Extr.distr_list.push_back( Obs_extrapolation_meson_mass(agm2s_strange_SD_ELM, X_2_fit, X_2_phys,  "../data/gm2/strange", "agm2_SD_ELM_Extrapolation_"+Extrapolation_strange_mode+"_"+V_strange_1_L.Tag[i_ens], UseJack, "SPLINE"));
 
  
  

    //####################################################################################################



    // ######################################## OSTERWALDER-SEILER #######################################

    //L
    distr_t agm2_W_OS_L(UseJack, UseJack?Njacks:Nboots);  //constructor sets  to zero by default
    distr_t agm2_SD_OS_L(UseJack, UseJack?Njacks:Nboots);  //constructor sets  to zero by default
    distr_t agm2_W_ELM_OS_L(UseJack, UseJack?Njacks:Nboots); //constructor sets  to zero by default
    distr_t agm2_SD_ELM_OS_L(UseJack, UseJack?Njacks:Nboots);  //constructor sets  to zero by default
    //M
    distr_t agm2_W_OS_M(UseJack, UseJack?Njacks:Nboots);  //constructor sets  to zero by default
    distr_t agm2_SD_OS_M(UseJack, UseJack?Njacks:Nboots);  //constructor sets  to zero by default
    distr_t agm2_W_ELM_OS_M(UseJack, UseJack?Njacks:Nboots); //constructor sets  to zero by default
    distr_t agm2_SD_ELM_OS_M(UseJack, UseJack?Njacks:Nboots);  //constructor sets  to zero by default
 
    //#################################################################################################

    for(int t=1; t< Corr.Nt/2; t++) {
      //L
      agm2_W_OS_L = agm2_W_OS_L +  w(t,Simps_ord)*4.0*pow(alpha,2)*Zv*Zv*V_strange_OS_distr_L.distr_list[t]*Ker.distr_list[t]*( distr_t::f_of_distr(th0, t*a_distr) - distr_t::f_of_distr(th1, t*a_distr));
      agm2_SD_OS_L = agm2_SD_OS_L +  w(t,Simps_ord)*4.0*pow(alpha,2)*Zv*Zv*(V_strange_OS_distr_L_pert_sub.distr_list[t])*Ker.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*a_distr));
      agm2_W_ELM_OS_L = agm2_W_ELM_OS_L +  w(t,Simps_ord)*4.0*pow(alpha,2)*Zv*Zv*V_strange_OS_distr_L.distr_list[t]*Ker_ELM_OS_L.distr_list[t]*( distr_t::f_of_distr(th0, t*ELM_mass_OS_L) - distr_t::f_of_distr(th1, t*ELM_mass_OS_L));
      agm2_SD_ELM_OS_L = agm2_SD_ELM_OS_L +  w(t,Simps_ord)*4.0*pow(alpha,2)*Zv*Zv*V_strange_OS_distr_L_pert_sub.distr_list[t]*Ker_ELM_OS_L.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*ELM_mass_OS_L));
      //M
      agm2_W_OS_M = agm2_W_OS_M +  w(t,Simps_ord)*4.0*pow(alpha,2)*Zv_heavy*Zv_heavy*V_strange_OS_distr_M.distr_list[t]*Ker.distr_list[t]*( distr_t::f_of_distr(th0, t*a_distr) - distr_t::f_of_distr(th1, t*a_distr));
      agm2_SD_OS_M = agm2_SD_OS_M +  w(t,Simps_ord)*4.0*pow(alpha,2)*Zv_heavy*Zv_heavy*(V_strange_OS_distr_M_pert_sub.distr_list[t])*Ker.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*a_distr));
      agm2_W_ELM_OS_M = agm2_W_ELM_OS_M +  w(t,Simps_ord)*4.0*pow(alpha,2)*Zv_heavy*Zv_heavy*V_strange_OS_distr_M.distr_list[t]*Ker_ELM_OS_M.distr_list[t]*( distr_t::f_of_distr(th0, t*ELM_mass_OS_M) - distr_t::f_of_distr(th1, t*ELM_mass_OS_M));
      agm2_SD_ELM_OS_M = agm2_SD_ELM_OS_M +  w(t,Simps_ord)*4.0*pow(alpha,2)*Zv_heavy*Zv_heavy*V_strange_OS_distr_M_pert_sub.distr_list[t]*Ker_ELM_OS_M.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*ELM_mass_OS_M));
    }
  
    //push_back the result

    //L
    agm2_strange_W_OS_L.distr_list.push_back(agm2_W_OS_L);
    agm2_strange_SD_OS_L.distr_list.push_back(agm2_SD_OS_L);
    agm2_strange_W_ELM_OS_L.distr_list.push_back(agm2_W_ELM_OS_L);
    agm2_strange_SD_ELM_OS_L.distr_list.push_back(agm2_SD_ELM_OS_L);
    //M
    agm2_strange_W_OS_M.distr_list.push_back(agm2_W_OS_M);
    agm2_strange_SD_OS_M.distr_list.push_back(agm2_SD_OS_M);
    agm2_strange_W_ELM_OS_M.distr_list.push_back(agm2_W_ELM_OS_M);
    agm2_strange_SD_ELM_OS_M.distr_list.push_back(agm2_SD_ELM_OS_M);
 

    //extrapolate the result to the physical kaon point
    vector<distr_t> agm2s_strange_W_OS({agm2_W_OS_L, agm2_W_OS_M});
    vector<distr_t> agm2s_strange_SD_OS({agm2_SD_OS_L, agm2_SD_OS_M});
    vector<distr_t> agm2s_strange_W_ELM_OS({agm2_W_ELM_OS_L, agm2_W_ELM_OS_M});
    vector<distr_t> agm2s_strange_SD_ELM_OS({agm2_SD_ELM_OS_L, agm2_SD_ELM_OS_M});

  
    agm2_strange_W_OS_Extr.distr_list.push_back( Obs_extrapolation_meson_mass(agm2s_strange_W_OS, X_2_fit, X_2_phys, "../data/gm2/strange", "agm2_W_OS_Extrapolation_"+Extrapolation_strange_mode+"_"+V_strange_1_L.Tag[i_ens], UseJack, "SPLINE"));
    agm2_strange_SD_OS_Extr.distr_list.push_back( Obs_extrapolation_meson_mass(agm2s_strange_SD_OS, X_2_fit, X_2_phys, "../data/gm2/strange", "agm2_SD_OS_Extrapolation_"+Extrapolation_strange_mode+"_"+V_strange_1_L.Tag[i_ens], UseJack, "SPLINE"));
    agm2_strange_W_ELM_OS_Extr.distr_list.push_back( Obs_extrapolation_meson_mass(agm2s_strange_W_ELM_OS, X_2_fit, X_2_phys, "../data/gm2/strange", "agm2_W_ELM_OS_Extrapolation_"+Extrapolation_strange_mode+"_"+V_strange_1_L.Tag[i_ens], UseJack, "SPLINE"));
    agm2_strange_SD_ELM_OS_Extr.distr_list.push_back( Obs_extrapolation_meson_mass(agm2s_strange_SD_ELM_OS, X_2_fit, X_2_phys, "../data/gm2/strange", "agm2_SD_ELM_OS_Extrapolation_"+Extrapolation_strange_mode+"_"+V_strange_1_L.Tag[i_ens], UseJack, "SPLINE"));

    //####################################################################################################



    //####################################### DISCO STRANGE ###############################################
    distr_t agm2_disco_W(UseJack, UseJack?Njacks:Nboots); //constructor sets it to 0
    distr_t agm2_disco_SD(UseJack, UseJack?Njacks:Nboots); //constructor sets it to 0
    distr_t agm2_disco_full(UseJack, UseJack?Njacks:Nboots); //constructor sets it to 0

    //improved
    distr_t agm2_disco_impr_W(UseJack, UseJack?Njacks:Nboots); //constructor sets it to 0
    distr_t agm2_disco_impr_SD(UseJack, UseJack?Njacks:Nboots); //constructor sets it to 0
    distr_t agm2_disco_impr_full(UseJack, UseJack?Njacks:Nboots); //constructor sets it to 0


    //off diagonal impr light-strange
    distr_t agm2_disco_impr_ls_W(UseJack, UseJack?Njacks:Nboots); //constructor sets it to 0
    distr_t agm2_disco_impr_ls_SD(UseJack, UseJack?Njacks:Nboots); //constructor sets it to 0
    distr_t agm2_disco_impr_ls_full(UseJack, UseJack?Njacks:Nboots); //constructor sets it to 0

    //off diagonal impr lightD-strange
    distr_t agm2_disco_impr_lDs_W(UseJack, UseJack?Njacks:Nboots); //constructor sets it to 0
    distr_t agm2_disco_impr_lDs_SD(UseJack, UseJack?Njacks:Nboots); //constructor sets it to 0
    distr_t agm2_disco_impr_lDs_full(UseJack, UseJack?Njacks:Nboots); //constructor sets it to 0


    //off diagonal impr strange-charm
    distr_t agm2_disco_impr_sc_W(UseJack, UseJack?Njacks:Nboots); //constructor sets it to 0
    distr_t agm2_disco_impr_sc_SD(UseJack, UseJack?Njacks:Nboots); //constructor sets it to 0
    distr_t agm2_disco_impr_sc_full(UseJack, UseJack?Njacks:Nboots); //constructor sets it to 0
    
 
    if(Include_strange_disco && Found_disco_ens) {
    
      for(int t=1; t< Corr.Nt/2; t++) {
	agm2_disco_W = agm2_disco_W +  w(t,Simps_ord)*4.0*pow(alpha,2)*Zv*Zv*disco_distr.distr_list[t]*Ker.distr_list[t]*( distr_t::f_of_distr(th0, t*a_distr) - distr_t::f_of_distr(th1, t*a_distr));
	agm2_disco_SD = agm2_disco_SD +  w(t,Simps_ord)*4.0*pow(alpha,2)*Zv*Zv*disco_distr.distr_list[t]*Ker.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*a_distr));
	agm2_disco_full = agm2_disco_full +  w(t,Simps_ord)*4.0*pow(alpha,2)*Zv*Zv*disco_distr.distr_list[t]*Ker.distr_list[t];
      }
  
      //push_back the result

      agm2_disco_strange_W.distr_list.push_back(agm2_disco_W);
      agm2_disco_strange_SD.distr_list.push_back(agm2_disco_SD);
      agm2_disco_strange_No_ELM.distr_list.push_back(agm2_disco_full);
    }

    //improved
    if(Include_strange_disco && Found_disco_impr_ens) {
    
      for(int t=1; t< Corr.Nt/2; t++) {
	agm2_disco_impr_W = agm2_disco_impr_W +  w(t,Simps_ord)*4.0*pow(alpha,2)*Zv*Zv*disco_impr_distr.distr_list[t]*Ker.distr_list[t]*( distr_t::f_of_distr(th0, t*a_distr) - distr_t::f_of_distr(th1, t*a_distr));
	agm2_disco_impr_SD = agm2_disco_impr_SD +  w(t,Simps_ord)*4.0*pow(alpha,2)*Zv*Zv*disco_impr_distr.distr_list[t]*Ker.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*a_distr));
	agm2_disco_impr_full = agm2_disco_impr_full +  w(t,Simps_ord)*4.0*pow(alpha,2)*Zv*Zv*disco_impr_distr.distr_list[t]*Ker.distr_list[t];
      }
  
      //push_back the result
      agm2_disco_impr_strange_W.distr_list.push_back(agm2_disco_impr_W);
      agm2_disco_impr_strange_SD.distr_list.push_back(agm2_disco_impr_SD);
      agm2_disco_impr_strange_No_ELM.distr_list.push_back(agm2_disco_impr_full);
    }

    

    //off-diagonal impr light-strange
    if(Include_off_diagonal_disco && Found_disco_impr_light_strange_ens) {
    
      for(int t=1; t< Corr.Nt/2; t++) {
	agm2_disco_impr_ls_W = agm2_disco_impr_ls_W +  w(t,Simps_ord)*4.0*pow(alpha,2)*Zv*Zv*disco_impr_light_strange_distr.distr_list[t]*Ker.distr_list[t]*( distr_t::f_of_distr(th0, t*a_distr) - distr_t::f_of_distr(th1, t*a_distr));
	agm2_disco_impr_ls_SD = agm2_disco_impr_ls_SD +  w(t,Simps_ord)*4.0*pow(alpha,2)*Zv*Zv*disco_impr_light_strange_distr.distr_list[t]*Ker.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*a_distr));
	agm2_disco_impr_ls_full = agm2_disco_impr_ls_full +  w(t,Simps_ord)*4.0*pow(alpha,2)*Zv*Zv*disco_impr_light_strange_distr.distr_list[t]*Ker.distr_list[t];
      }
  
      //push_back the result
      agm2_disco_impr_light_strange_W.distr_list.push_back(agm2_disco_impr_ls_W);
      agm2_disco_impr_light_strange_SD.distr_list.push_back(agm2_disco_impr_ls_SD);
      agm2_disco_impr_light_strange_No_ELM.distr_list.push_back(agm2_disco_impr_ls_full);
    }


    //off-diagonal impr lightD-strange
    if(Include_off_diagonal_disco && Found_disco_impr_lightD_strange_ens) {
    
      for(int t=1; t< Corr.Nt/2; t++) {
	agm2_disco_impr_lDs_W = agm2_disco_impr_lDs_W +  w(t,Simps_ord)*4.0*pow(alpha,2)*Zv*Zv*disco_impr_lightD_strange_distr.distr_list[t]*Ker.distr_list[t]*( distr_t::f_of_distr(th0, t*a_distr) - distr_t::f_of_distr(th1, t*a_distr));
	agm2_disco_impr_lDs_SD = agm2_disco_impr_lDs_SD +  w(t,Simps_ord)*4.0*pow(alpha,2)*Zv*Zv*disco_impr_lightD_strange_distr.distr_list[t]*Ker.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*a_distr));
	agm2_disco_impr_lDs_full = agm2_disco_impr_lDs_full +  w(t,Simps_ord)*4.0*pow(alpha,2)*Zv*Zv*disco_impr_lightD_strange_distr.distr_list[t]*Ker.distr_list[t];
      }
  
      //push_back the result
      agm2_disco_impr_lightD_strange_W.distr_list.push_back(agm2_disco_impr_lDs_W);
      agm2_disco_impr_lightD_strange_SD.distr_list.push_back(agm2_disco_impr_lDs_SD);
      agm2_disco_impr_lightD_strange_No_ELM.distr_list.push_back(agm2_disco_impr_lDs_full);
    }


    //off-diagonal impr strange-charm
    if(Include_off_diagonal_disco && Found_disco_impr_strange_charm_ens) {
    
      for(int t=1; t< Corr.Nt/2; t++) {
	agm2_disco_impr_sc_W = agm2_disco_impr_sc_W +  w(t,Simps_ord)*4.0*pow(alpha,2)*Zv*Zv*disco_impr_strange_charm_distr.distr_list[t]*Ker.distr_list[t]*( distr_t::f_of_distr(th0, t*a_distr) - distr_t::f_of_distr(th1, t*a_distr));
	agm2_disco_impr_sc_SD = agm2_disco_impr_sc_SD +  w(t,Simps_ord)*4.0*pow(alpha,2)*Zv*Zv*disco_impr_strange_charm_distr.distr_list[t]*Ker.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*a_distr));
	agm2_disco_impr_sc_full = agm2_disco_impr_sc_full +  w(t,Simps_ord)*4.0*pow(alpha,2)*Zv*Zv*disco_impr_strange_charm_distr.distr_list[t]*Ker.distr_list[t];
      }
      
      //push_back the result
      agm2_disco_impr_strange_charm_W.distr_list.push_back(agm2_disco_impr_sc_W);
      agm2_disco_impr_strange_charm_SD.distr_list.push_back(agm2_disco_impr_sc_SD);
      agm2_disco_impr_strange_charm_No_ELM.distr_list.push_back(agm2_disco_impr_sc_full);
    }

    

    
    

  
  
    //print total contribution to file
    //L
    Print_To_File({}, {Tdata_vec, agm2_distr_Tdata_L.ave(), agm2_distr_Tdata_L.err()}, "../data/gm2/strange/tm_"+Extrapolation_strange_mode+"/agm2_Tdata_ELM_"+V_strange_1_L.Tag[i_ens]+"_L.dat.t", "", "#id  Tdata   ag2m agm2_err");
    Print_To_File({}, {Tdata_vec, agm2_OS_distr_Tdata_L.ave(), agm2_OS_distr_Tdata_L.err()}, "../data/gm2/strange/OS_"+Extrapolation_strange_mode+"/agm2_Tdata_ELM_"+V_strange_1_L.Tag[i_ens]+"_L.dat.t", "", "#id  Tdata   ag2m agm2_err");
    Print_To_File({}, {Tdata_vec, agm2_distr_Tdata_No_ELM_L.ave(), agm2_distr_Tdata_No_ELM_L.err()}, "../data/gm2/strange/tm_"+Extrapolation_strange_mode+"/agm2_Tdata_"+V_strange_1_L.Tag[i_ens]+"_L.dat.t", "", "#id  Tdata   ag2m agm2_err");
    Print_To_File({}, {Tdata_vec, agm2_OS_distr_Tdata_No_ELM_L.ave(), agm2_OS_distr_Tdata_No_ELM_L.err()}, "../data/gm2/strange/OS_"+Extrapolation_strange_mode+"/agm2_Tdata_"+V_strange_1_L.Tag[i_ens]+"_L.dat.t", "", "#id  Tdata   ag2m agm2_err");
    //M
    Print_To_File({}, {Tdata_vec, agm2_distr_Tdata_M.ave(), agm2_distr_Tdata_M.err()}, "../data/gm2/strange/tm_"+Extrapolation_strange_mode+"/agm2_Tdata_ELM_"+V_strange_1_M.Tag[i_ens]+"_M.dat.t", "", "#id  Tdata   ag2m agm2_err");
    Print_To_File({}, {Tdata_vec, agm2_OS_distr_Tdata_M.ave(), agm2_OS_distr_Tdata_M.err()}, "../data/gm2/strange/OS_"+Extrapolation_strange_mode+"/agm2_Tdata_ELM_"+V_strange_1_M.Tag[i_ens]+"_M.dat.t", "", "#id  Tdata   ag2m agm2_err");
    Print_To_File({}, {Tdata_vec, agm2_distr_Tdata_No_ELM_M.ave(), agm2_distr_Tdata_No_ELM_M.err()}, "../data/gm2/strange/tm_"+Extrapolation_strange_mode+"/agm2_Tdata_"+V_strange_1_M.Tag[i_ens]+"_M.dat.t", "", "#id  Tdata   ag2m agm2_err");
    Print_To_File({}, {Tdata_vec, agm2_OS_distr_Tdata_No_ELM_M.ave(), agm2_OS_distr_Tdata_No_ELM_M.err()}, "../data/gm2/strange/OS_"+Extrapolation_strange_mode+"/agm2_Tdata_"+V_strange_1_M.Tag[i_ens]+"_M.dat.t", "", "#id  Tdata   ag2m agm2_err");




    
    //######################   PI(Q^2) analysis   ###########################


    distr_t_list PI_Q2_L_tm, PI_Q2_L_OS, PI_Q2_M_tm, PI_Q2_M_OS;
    distr_t_list PI_Q2_L_tm_pert_sub, PI_Q2_L_OS_pert_sub, PI_Q2_M_tm_pert_sub, PI_Q2_M_OS_pert_sub;

    
    Vint Tdatas_opt_L_tm, Tdatas_opt_L_OS, Tdatas_opt_M_tm, Tdatas_opt_M_OS;
    Vint Tdatas_opt_L_tm_pert_sub, Tdatas_opt_L_OS_pert_sub, Tdatas_opt_M_tm_pert_sub, Tdatas_opt_M_OS_pert_sub;

   

    

    
    //Apply bounding method
    Bounding_PI_q2(PI_Q2_L_tm, Za*Za*V_strange_distr_L, a_distr, "../data/PI_Q2/strange/tm_"+Extrapolation_strange_mode+"/PI_Q2_Tdata_"+V_strange_1_L.Tag[i_ens]+"_L"  , Tdatas_opt_L_tm, MV_strange_L );
    Bounding_PI_q2(PI_Q2_M_tm, Za*Za*V_strange_distr_M, a_distr, "../data/PI_Q2/strange/tm_"+Extrapolation_strange_mode+"/PI_Q2_Tdata_"+V_strange_1_M.Tag[i_ens]+"_M"  , Tdatas_opt_M_tm, MV_strange_M);
    Bounding_PI_q2(PI_Q2_L_OS, Zv*Zv*V_strange_OS_distr_L, a_distr, "../data/PI_Q2/strange/OS_"+Extrapolation_strange_mode+"/PI_Q2_Tdata_"+V_strange_1_L.Tag[i_ens]+"_L"  , Tdatas_opt_L_OS, MV_strange_OS_L);
    Bounding_PI_q2(PI_Q2_M_OS, Zv*Zv*V_strange_OS_distr_M, a_distr, "../data/PI_Q2/strange/OS_"+Extrapolation_strange_mode+"/PI_Q2_Tdata_"+V_strange_1_M.Tag[i_ens]+"_M"  , Tdatas_opt_M_OS, MV_strange_OS_M);

    //include perturbative subtraction
    Get_PI_q2(PI_Q2_L_tm_pert_sub, VV_free_oppor_L, a_distr, Corr.Nt/2);
    Get_PI_q2(PI_Q2_M_tm_pert_sub, VV_free_oppor_M, a_distr, Corr.Nt/2);
    Get_PI_q2(PI_Q2_L_OS_pert_sub, VV_free_samer_L, a_distr, Corr.Nt/2);
    Get_PI_q2(PI_Q2_M_OS_pert_sub, VV_free_samer_M, a_distr, Corr.Nt/2);

    PI_Q2_L_tm_pert_sub = PI_Q2_L_tm_pert_sub + PI_Q2_L_tm;
    PI_Q2_M_tm_pert_sub = PI_Q2_M_tm_pert_sub + PI_Q2_M_tm;
    PI_Q2_L_OS_pert_sub = PI_Q2_L_OS_pert_sub + PI_Q2_L_OS;
    PI_Q2_M_OS_pert_sub = PI_Q2_M_OS_pert_sub + PI_Q2_M_OS;
    


    //push back OS strange correlator for bounding on disconnected, and the disconnected correlator 

    if(V_strange_1_L.Tag[i_ens] == "cD211a.054.96" || V_strange_1_L.Tag[i_ens] == "cC211a.06.80" || V_strange_1_L.Tag[i_ens] == "cB211b.072.64") {
      int id_disco_ens=0;
      if(V_strange_1_L.Tag[i_ens] == "cD211a.054.96") id_disco_ens=2;
      else if(V_strange_1_L.Tag[i_ens] == "cC211a.06.80") id_disco_ens=1;
      else if(V_strange_1_L.Tag[i_ens] == "cB211b.072.64") id_disco_ens=0;
      else crash("what");
      
      CORR_STRANGE_FOR_PI_Q2[id_disco_ens] = Zv*Zv*V_strange_OS_distr_M;
      STRANGE_PI_Q2_FOR_DISCO[id_disco_ens] =PI_Q2_M_OS;
      if(Include_strange_disco && Include_off_diagonal_disco) {
	if(!Found_disco_impr_ens || !Found_disco_impr_lightD_strange_ens || !Found_disco_impr_strange_charm_ens) crash("disconnected for PI(Q^2) not found in strange ens");
	CORR_DISCO_FOR_PI_Q2[id_disco_ens] = Zv*Zv*(disco_impr_distr + disco_impr_lightD_strange_distr+ disco_impr_strange_charm_distr);
      }
      else crash("what?");
    }



    //interpolate to the physical point

    distr_t_list PI_Q2_Extr_tm, PI_Q2_Extr_OS;
    distr_t_list PI_Q2_Extr_tm_pert_sub, PI_Q2_Extr_OS_pert_sub;

    Vfloat Tcut_f_tm, Tcut_f_OS, Tcut_f_tm_pert_sub, Tcut_f_OS_pert_sub;

    
    int Q_size= Qs2.size();
    for(int q=0; q < Q_size;q++) {

      vector<distr_t> PIs_q2_strange_tm({ PI_Q2_L_tm.distr_list[q], PI_Q2_M_tm.distr_list[q]});
      vector<distr_t> PIs_q2_strange_OS({ PI_Q2_L_OS.distr_list[q], PI_Q2_M_OS.distr_list[q]});
      vector<distr_t> PIs_q2_strange_tm_pert_sub({ PI_Q2_L_tm_pert_sub.distr_list[q], PI_Q2_M_tm_pert_sub.distr_list[q]});
      vector<distr_t> PIs_q2_strange_OS_pert_sub({ PI_Q2_L_OS_pert_sub.distr_list[q], PI_Q2_M_OS_pert_sub.distr_list[q]});
      
      PI_Q2_Extr_tm.distr_list.push_back(Obs_extrapolation_meson_mass(PIs_q2_strange_tm, X_2_fit, X_2_phys,  "../data/PI_Q2/strange", "PI_Q2_tm_extr_"+Extrapolation_strange_mode+"_"+V_strange_1_L.Tag[i_ens]+"_Q2_"+to_string_with_precision(Qs2[q],5), UseJack, "SPLINE"));

      PI_Q2_Extr_OS.distr_list.push_back(Obs_extrapolation_meson_mass(PIs_q2_strange_OS, X_2_fit, X_2_phys,  "../data/PI_Q2/strange", "PI_Q2_OS_extr_"+Extrapolation_strange_mode+"_"+V_strange_1_M.Tag[i_ens]+"_Q2_"+to_string_with_precision(Qs2[q],5), UseJack, "SPLINE"));
      PI_Q2_Extr_tm_pert_sub.distr_list.push_back(Obs_extrapolation_meson_mass(PIs_q2_strange_tm_pert_sub, X_2_fit, X_2_phys,  "../data/PI_Q2/strange", "PI_Q2_tm_pert_sub_extr_"+Extrapolation_strange_mode+"_"+V_strange_1_L.Tag[i_ens]+"_Q2_"+to_string_with_precision(Qs2[q],5), UseJack, "SPLINE"));
      PI_Q2_Extr_OS_pert_sub.distr_list.push_back(Obs_extrapolation_meson_mass(PIs_q2_strange_OS_pert_sub, X_2_fit, X_2_phys,  "../data/PI_Q2/strange", "PI_Q2_OS_pert_sub_extr_"+Extrapolation_strange_mode+"_"+V_strange_1_M.Tag[i_ens]+"_Q2_"+to_string_with_precision(Qs2[q],5), UseJack, "SPLINE"));



      //check if bounding method worked for both L and M (1=worked, 0= not worked)
      if( Tdatas_opt_L_tm[q] > 0 && Tdatas_opt_M_tm[q] > 0)  Tcut_f_tm.push_back( 1.0);
      else Tcut_f_tm.push_back( 0.0);
      if( Tdatas_opt_L_OS[q] > 0 && Tdatas_opt_M_OS[q] > 0) Tcut_f_OS.push_back( 1.0);
      else Tcut_f_OS.push_back( 0.0);
      //if( Tdatas_opt_L_tm_pert_sub[q] > 0 && Tdatas_opt_M_tm_pert_sub[q] > 0) Tcut_f_tm_pert_sub.push_back( 1.0);
      //else Tcut_f_tm_pert_sub.push_back( 0.0);
      //if( Tdatas_opt_L_OS_pert_sub[q] > 0 && Tdatas_opt_M_OS_pert_sub[q] > 0) Tcut_f_OS_pert_sub.push_back( 1.0);
      //else Tcut_f_OS_pert_sub.push_back( 0.0);

      
    }
    


    //push_back and print to file

    Add_ens_val_PI_q2( PI_Q2_strange_tm, PI_Q2_Extr_tm);
    Add_ens_val_PI_q2( PI_Q2_strange_OS, PI_Q2_Extr_OS);
    Add_ens_val_PI_q2( PI_Q2_strange_tm_pert_sub, PI_Q2_Extr_tm_pert_sub);
    Add_ens_val_PI_q2( PI_Q2_strange_OS_pert_sub, PI_Q2_Extr_OS_pert_sub);


    Print_To_File({}, {Qs2, PI_Q2_Extr_tm.ave(), PI_Q2_Extr_tm.err(), Tcut_f_tm } , "../data/PI_Q2/strange/tm_"+Extrapolation_strange_mode+"/PI_Q2_extr_"+V_strange_1_L.Tag[i_ens]+".t", "", "# Q2[GeV2]   PI(Q^2)   Tcut_f");
    Print_To_File({}, {Qs2, PI_Q2_Extr_OS.ave(), PI_Q2_Extr_OS.err(), Tcut_f_OS } , "../data/PI_Q2/strange/OS_"+Extrapolation_strange_mode+"/PI_Q2_extr_"+V_strange_1_L.Tag[i_ens]+".t", "", "# Q2[GeV2]   PI(Q^2)   Tcut_f");
    Print_To_File({}, {Qs2, PI_Q2_Extr_tm_pert_sub.ave(), PI_Q2_Extr_tm_pert_sub.err(), Tcut_f_tm } , "../data/PI_Q2/strange/tm_"+Extrapolation_strange_mode+"/PI_Q2_extr_pert_sub_"+V_strange_1_L.Tag[i_ens]+".t", "", "# Q2[GeV2]   PI(Q^2)   Tcut_f");
    Print_To_File({}, {Qs2, PI_Q2_Extr_OS_pert_sub.ave(), PI_Q2_Extr_OS_pert_sub.err(), Tcut_f_OS } , "../data/PI_Q2/strange/OS_"+Extrapolation_strange_mode+"/PI_Q2_extr_pert_sub_"+V_strange_1_L.Tag[i_ens]+".t", "", "# Q2[GeV2]   PI(Q^2)   Tcut_f");




    //#######################################################################
  
  }

  cout<<"strange quark correlator analyzed!"<<endl;



   

  //charm
  channel="c";
  for(int i_ens=0;i_ens<Nens_charm;i_ens++) { 
    
    CorrAnalysis Corr(UseJack, Njacks,Nboots);
    Corr.Nt = V_charm_1_L.nrows[i_ens];


  

    //resample lattice spacing
    distr_t a_distr(UseJack), Za_WI_distr(UseJack), Zv_WI_distr(UseJack), Za_RIMOM_distr(UseJack), Zv_RIMOM_distr(UseJack), Zv_WI_strange_distr(UseJack), Za_WI_strange_distr(UseJack);
    distr_t Zv_WI_charm_extr_distr(UseJack), Za_WI_charm_extr_distr(UseJack);
    LatticeInfo L_info;
    L_info.LatInfo_new_ens(V_charm_1_L.Tag[i_ens]);
    //generate jackknife sample of input parameters
  
    if(UseJack)  { for(int ijack=0;ijack<Njacks;ijack++) {
	Za_WI_distr.distr.push_back(  L_info.Za_WI + GM()*L_info.Za_WI_err*(1.0/sqrt(Njacks-1.0)));
	Zv_WI_distr.distr.push_back(  L_info.Zv_WI + GM()*L_info.Zv_WI_err*(1.0/sqrt(Njacks-1.0)));
	Za_RIMOM_distr.distr.push_back( L_info.Za + GM()*L_info.Za_err*(1.0/sqrt(Njacks-1.0)));
	Zv_RIMOM_distr.distr.push_back( L_info.Zv + GM()*L_info.Zv_err*(1.0/sqrt(Njacks-1.0)));
	Za_WI_strange_distr.distr.push_back( L_info.Za_WI_strange + GM()*L_info.Za_WI_strange_err*(1.0/sqrt(Njacks-1.0)));
	Zv_WI_strange_distr.distr.push_back( L_info.Zv_WI_strange + GM()*L_info.Zv_WI_strange_err*(1.0/sqrt(Njacks-1.0)));
	if(V_charm_1_L.Tag[i_ens].substr(1,1)=="A") {
	  Zv_WI_charm_extr_distr.distr.push_back( L_info.Zv_WI_charm_extr + GM()*L_info.Zv_WI_charm_extr_err*(1.0/sqrt(Njacks-1.0)));
	  Za_WI_charm_extr_distr.distr.push_back( L_info.Za_WI_charm_extr + GM()*L_info.Za_WI_charm_extr_err*(1.0/sqrt(Njacks-1.0)));
	}
      }
    }
    else {
      for (int iboot=0; iboot<Nboots;iboot++) {
	Za_WI_distr.distr.push_back(  L_info.Za_WI + GM()*L_info.Za_WI_err);
	Zv_WI_distr.distr.push_back(  L_info.Zv_WI + GM()*L_info.Zv_WI_err);
	Za_RIMOM_distr.distr.push_back(  L_info.Za + GM()*L_info.Za_err);
	Zv_RIMOM_distr.distr.push_back(  L_info.Zv + GM()*L_info.Zv_err);
	Za_WI_strange_distr.distr.push_back( L_info.Za_WI_strange + GM()*L_info.Za_WI_strange_err);
	Zv_WI_strange_distr.distr.push_back( L_info.Zv_WI_strange + GM()*L_info.Zv_WI_strange_err);
	if(V_charm_1_L.Tag[i_ens].substr(1,1)=="A") {
	  Zv_WI_charm_extr_distr.distr.push_back( L_info.Zv_WI_charm_extr + GM()*L_info.Zv_WI_charm_extr_err);
	  Za_WI_charm_extr_distr.distr.push_back( L_info.Za_WI_charm_extr + GM()*L_info.Za_WI_charm_extr_err);
	}
      }
    }

    distr_t Za_L = Za_WI_strange_distr;
    distr_t Za_M = Za_WI_strange_distr;
    distr_t Za_H = Za_WI_strange_distr;
    distr_t Zv_L = Zv_WI_strange_distr;
    distr_t Zv_M = Zv_WI_strange_distr;
    distr_t Zv_H = Zv_WI_strange_distr;

  



    int Tmin_P5P5;
    int Tmax_P5P5;
    int Tmin_VV;
    int Tmax_VV;

  
    //set time intervals for pseudoscalar obs
    if(V_charm_1_L.Tag[i_ens].substr(1,1) == "C") {
      a_distr = a_C;
      if(V_charm_1_L.Tag[i_ens]=="cC211a.06.80") { Corr.Tmin=23; Corr.Tmax=39;}
      else crash("Cannot find ensemble tag: "+V_charm_1_L.Tag[i_ens]);
    }
    else if(V_charm_1_L.Tag[i_ens].substr(1,1) == "B") {
      a_distr = a_B;
      if(V_charm_1_L.Tag[i_ens]== "cB211a.14.64") {Corr.Tmin=19; Corr.Tmax=36;}
      else if(V_charm_1_L.Tag[i_ens] == "cB211a.25.48") {Corr.Tmin=21;Corr.Tmax=34;}
      else if(V_charm_1_L.Tag[i_ens] == "cB211b.072.64") {Corr.Tmin=16; Corr.Tmax= 27;}
      else crash("Cannot find ensemble tag: "+V_charm_1_L.Tag[i_ens]);
    }
    else if(V_charm_1_L.Tag[i_ens].substr(1,1) == "A") {
      a_distr = a_A;
      if(V_charm_1_L.Tag[i_ens] == "cA211a.12.48") {Corr.Tmin=15; Corr.Tmax=28;}
      else if(V_charm_1_L.Tag[i_ens] == "cA211a.40.24") {Corr.Tmin=14; Corr.Tmax=21;}
      else if(V_charm_1_L.Tag[i_ens] == "cA211a.53.24") {Corr.Tmin=15; Corr.Tmax=22;}
      else if(V_charm_1_L.Tag[i_ens] == "cA211ab.30.32") {Corr.Tmin=15; Corr.Tmax=26;}
      else crash("Cannot find ensemble tag: "+V_charm_1_L.Tag[i_ens]);
    }
    else if(V_charm_1_L.Tag[i_ens].substr(1,1) == "D") {
      a_distr = a_D;
      if(V_charm_1_L.Tag[i_ens]== "cD211a.054.96") {Corr.Tmin=36; Corr.Tmax=80;}
      else crash("Cannot find ensemble tag: "+V_charm_1_L.Tag[i_ens]);
    }
  
    else crash("Ensemble tag not valid");

    //push_back lattice info
    L_charm_list.push_back(L_info.L);
    a_charm_list.push_back(a_distr.ave()/fm_to_inv_Gev);
    ml_charm_list.push_back(L_info.ml);

    a_distr_list_charm.distr_list.push_back(a_distr);

    //set Tmin_P5P5 and Tmax_P5P5 to the values Corr.Tmin and Corr.Tmax
    Tmin_P5P5 = Corr.Tmin;
    Tmax_P5P5 = Corr.Tmax;

  
    //set time intervals for vector obs
    if(V_charm_1_L.Tag[i_ens].substr(1,1) == "C") {
      if(V_charm_1_L.Tag[i_ens]=="cC211a.06.80") { Tmin_VV=43; Tmax_VV=58;}
      else crash("Cannot find ensemble tag: "+V_charm_1_L.Tag[i_ens]);
    }
    else if(V_charm_1_L.Tag[i_ens].substr(1,1) == "B") {
      if(V_charm_1_L.Tag[i_ens]== "cB211a.14.64") {Tmin_VV=31; Tmax_VV=43;}
      else if(V_charm_1_L.Tag[i_ens] == "cB211a.25.48") {Tmin_VV=32;Tmax_VV=42;}
      else if(V_charm_1_L.Tag[i_ens] == "cB211b.072.64") {Tmin_VV=28; Tmax_VV= 40;}
      else crash("Cannot find ensemble tag: "+V_charm_1_L.Tag[i_ens]);
    }
    else if(V_charm_1_L.Tag[i_ens].substr(1,1) == "A") {
      if(V_charm_1_L.Tag[i_ens] == "cA211a.12.48") {Tmin_VV=27; Tmax_VV=42;}
      else if(V_charm_1_L.Tag[i_ens] == "cA211a.40.24") {Tmin_VV=19; Tmax_VV=24;}
      else if(V_charm_1_L.Tag[i_ens] == "cA211a.53.24") {Tmin_VV=20; Tmax_VV=24;}
      else if(V_charm_1_L.Tag[i_ens] == "cA211ab.30.32") {Tmin_VV=24;Tmax_VV=31;}
      else crash("Cannot find ensemble tag: "+V_charm_1_L.Tag[i_ens]);
    }
    else if(V_charm_1_L.Tag[i_ens].substr(1,1) == "D") {
      if(V_charm_1_L.Tag[i_ens]== "cD211a.054.96") {Tmin_VV=49; Tmax_VV=65;}
      else crash("Cannot find ensemble tag: "+V_charm_1_L.Tag[i_ens]);
    }

    else crash("Ensemble tag not valid");

  

    //set Tmin and Tmax for the eta_C
    int Tmin_etaC, Tmax_etaC;
    if(V_charm_1_L.Tag[i_ens].substr(1,1) == "C") {
      if(V_charm_1_L.Tag[i_ens]=="cC211a.06.80") { Tmin_etaC=40; Tmax_etaC=65;}
      else crash("Cannot find ensemble tag: "+V_charm_1_L.Tag[i_ens]);
    }
    else if(V_charm_1_L.Tag[i_ens].substr(1,1) == "B") {
      if(V_charm_1_L.Tag[i_ens]== "cB211a.14.64") {Tmin_etaC=35; Tmax_etaC=50;}
      else if(V_charm_1_L.Tag[i_ens] == "cB211a.25.48") {Tmin_etaC=32;Tmax_etaC=45;}
      else if(V_charm_1_L.Tag[i_ens] == "cB211b.072.64") {Tmin_etaC=41; Tmax_etaC= 57;}
      else crash("Cannot find ensemble tag: "+V_charm_1_L.Tag[i_ens]);
    }
    else if(V_charm_1_L.Tag[i_ens].substr(1,1) == "A") {
      if(V_charm_1_L.Tag[i_ens] == "cA211a.12.48") {Tmin_etaC=28; Tmax_etaC=45;}
      else if(V_charm_1_L.Tag[i_ens] == "cA211a.40.24") {Tmin_etaC=20; Tmax_etaC=23;}
      else if(V_charm_1_L.Tag[i_ens] == "cA211a.53.24") {Tmin_etaC=20; Tmax_etaC=23;}
      else if(V_charm_1_L.Tag[i_ens] == "cA211ab.30.32") {Tmin_etaC=25;Tmax_etaC=31;}
      else crash("Cannot find ensemble tag: "+V_charm_1_L.Tag[i_ens]);
    }
    else if(V_charm_1_L.Tag[i_ens].substr(1,1) == "D") {
      if(V_charm_1_L.Tag[i_ens]== "cD211a.054.96") {Tmin_etaC=50; Tmax_etaC=71;}
      else crash("Cannot find ensemble tag: "+V_charm_1_L.Tag[i_ens]);
    }
  
    else crash("Ensemble tag not valid");
  

    //tm
    distr_t_list Mpi_distr, fpi_distr;
    distr_t Mpi, fpi;
    //L
    distr_t_list  V_charm_1_distr_L, V_charm_2_distr_L, V_charm_3_distr_L;
    distr_t_list  V_charm_distr_L, MV_charm_distr_L, ZV_charm_distr_L;
    distr_t_list  M_etaC_distr_L, M_etaC_distr_OS_L;
    distr_t MV_charm_L , ZV_charm_L,M_etaC_L;
    distr_t_list fetac_L_distr;
    distr_t fetac_L;
    //M
    distr_t_list  V_charm_1_distr_M, V_charm_2_distr_M, V_charm_3_distr_M;
    distr_t_list  V_charm_distr_M, MV_charm_distr_M, ZV_charm_distr_M;
    distr_t_list  M_etaC_distr_M, M_etaC_distr_OS_M;
    distr_t MV_charm_M, ZV_charm_M,M_etaC_M;
    distr_t_list fetac_M_distr;
    distr_t fetac_M;
    //H
    distr_t_list  V_charm_1_distr_H, V_charm_2_distr_H, V_charm_3_distr_H;
    distr_t_list  V_charm_distr_H, MV_charm_distr_H, ZV_charm_distr_H;
    distr_t_list   M_etaC_distr_H, M_etaC_distr_OS_H;
    distr_t MV_charm_H , ZV_charm_H, M_etaC_H;
    distr_t_list fetac_H_distr;
    distr_t fetac_H;


    //OS
    distr_t_list Mpi_OS_distr;
    distr_t Mpi_OS;
    //L
    distr_t_list  V_charm_OS_1_distr_L, V_charm_OS_2_distr_L, V_charm_OS_3_distr_L;
    distr_t_list  V_charm_OS_distr_L, MV_charm_OS_distr_L, ZV_charm_OS_distr_L;
    distr_t MV_charm_OS_L , ZV_charm_OS_L;
    //M
    distr_t_list  V_charm_OS_1_distr_M, V_charm_OS_2_distr_M, V_charm_OS_3_distr_M;
    distr_t_list  V_charm_OS_distr_M, MV_charm_OS_distr_M, ZV_charm_OS_distr_M;
    distr_t MV_charm_OS_M, ZV_charm_OS_M;
    //H
    distr_t_list  V_charm_OS_1_distr_H, V_charm_OS_2_distr_H, V_charm_OS_3_distr_H;
    distr_t_list  V_charm_OS_distr_H, MV_charm_OS_distr_H, ZV_charm_OS_distr_H;
    distr_t MV_charm_OS_H , ZV_charm_OS_H;


    //observables to extract Zv and Za (Hadronic method)  (L)
    distr_t_list cbar_c_distr_L, cbar_c_OS_distr_L;
    distr_t_list overlap_P5P5_distr_L, overlap_P5P5_OS_distr_L;
    distr_t_list ratio_P5P5_overlap_OS_tm_L,  Zp_ov_Zs_distr_L;
    distr_t_list  A0P5_distr_L, A0P5_OS_distr_L;
    distr_t_list RA_L, RV_L;
    distr_t cbar_c_mass_L, cbar_c_OS_mass_L;
    distr_t Zp_ov_Zs_L, Zv_hadr_L, Za_hadr_L;

    //observables to extract Zv and Za (Hadronic method)  (M)
    distr_t_list cbar_c_distr_M, cbar_c_OS_distr_M;
    distr_t_list overlap_P5P5_distr_M, overlap_P5P5_OS_distr_M;
    distr_t_list ratio_P5P5_overlap_OS_tm_M,  Zp_ov_Zs_distr_M;
    distr_t_list  A0P5_distr_M, A0P5_OS_distr_M;
    distr_t_list RA_M, RV_M;
    distr_t cbar_c_mass_M, cbar_c_OS_mass_M;
    distr_t Zp_ov_Zs_M, Zv_hadr_M, Za_hadr_M;

    //observables to extract Zv and Za (Hadronic method)  (H)
    distr_t_list cbar_c_distr_H, cbar_c_OS_distr_H;
    distr_t_list overlap_P5P5_distr_H, overlap_P5P5_OS_distr_H;
    distr_t_list ratio_P5P5_overlap_OS_tm_H,  Zp_ov_Zs_distr_H;
    distr_t_list  A0P5_distr_H, A0P5_OS_distr_H;
    distr_t_list RA_H, RV_H;
    distr_t cbar_c_mass_H, cbar_c_OS_mass_H;
    distr_t Zp_ov_Zs_H, Zv_hadr_H, Za_hadr_H;

    //observables to extract Zv and Za from light correlator (Hadronic method)
    distr_t_list P5P5_charm_pion_distr, P5P5_OS_charm_pion_distr;
    distr_t_list A0P5_charm_pion_distr, A0P5_OS_charm_pion_distr;
    distr_t_list overlap_P5P5_distr_light, overlap_P5P5_OS_distr_light;
    distr_t_list ratio_P5P5_overlap_OS_tm_light, Zp_ov_Zs_distr_light;
    distr_t_list RA_light, RV_light, RA0_light;
    distr_t Zp_ov_Zs_light, Zv_hadr_light, Za_hadr_light;


    //disco
    distr_t_list disco_distr;
    //improved
    distr_t_list disco_impr_distr;
    //off_diagonal light charm
    //improved
    distr_t_list disco_impr_light_charm_distr, disco_impr_lightD_charm_distr;


  
    //Analyze correlators
  

  
    //tm and OS pion
    if(V_charm_1_L.Tag[i_ens].substr(1,1) == "A") {
      //search run
      int j_ens=0;
      bool found_jens=false;
      while( !found_jens) {
	if(pt2_pion_charm.Tag[j_ens] != V_charm_1_L.Tag[i_ens]) j_ens++;
	else found_jens=true;
      
	if(!found_jens && j_ens == Nens_light_silvano) crash("Cannot find Silvano-light Ensemble: "+V_charm_1_L.Tag[i_ens]);
      }
      Mpi_distr = Corr.effective_mass_t(pt2_pion_charm.col(0)[j_ens], "../data/gm2/charm/tm_"+Extrapolation_charm_mode+"/Mpi_"+V_charm_1_L.Tag[i_ens]+".dat");
      Mpi_OS_distr = Corr.effective_mass_t(pt2_pion_OS_charm.col(0)[j_ens], "../data/gm2/charm/OS_"+Extrapolation_charm_mode+"/Mpi_"+V_charm_1_L.Tag[i_ens]+".dat");  
      //fpi
      fpi_distr = 2.0*L_info.ml*Corr.decay_constant_t(pt2_pion_charm.col(0)[j_ens], "../data/gm2/charm/tm_"+Extrapolation_charm_mode+"/fpi_"+V_charm_1_L.Tag[i_ens]+".dat");
      P5P5_charm_pion_distr = Corr.corr_t(pt2_pion_charm.col(0)[j_ens], "");
      P5P5_OS_charm_pion_distr = Corr.corr_t(pt2_pion_OS_charm.col(0)[j_ens], "");
      if(V_charm_1_L.Tag[i_ens] != "cA211ab.30.32") Corr.Reflection_sign = -1;
      A0P5_charm_pion_distr= Corr.corr_t(corr_A0P5_charm_pion.col(0)[j_ens], "");
      A0P5_OS_charm_pion_distr = Corr.corr_t(corr_A0P5_OS_charm_pion.col(0)[j_ens], "");
      Corr.Reflection_sign = 1;
    }
    else { //read from light run
      //search run
      int j_ens=0;
      bool found_jens=false;
      while( !found_jens) { if(V_light_1.Tag[j_ens] != V_charm_1_L.Tag[i_ens]) j_ens++; else { found_jens=true;} if(!found_jens && j_ens == Nens_light) crash("Cannot find light Ensemble: "+V_charm_1_L.Tag[i_ens]);}
      Mpi_distr= Corr.effective_mass_t(pt2_pion.col(0)[j_ens], "../data/gm2/charm/tm_"+Extrapolation_charm_mode+"/Mpi_"+V_charm_1_L.Tag[i_ens]+".dat");
      Mpi_OS_distr= Corr.effective_mass_t(corr_P5P5_OS.col(0)[j_ens], "../data/gm2/charm/OS_"+Extrapolation_charm_mode+"/Mpi_"+V_charm_1_L.Tag[i_ens]+".dat");
      //fpi
      fpi_distr = 2.0*L_info.ml*Corr.decay_constant_t(pt2_pion.col(0)[j_ens], "../data/gm2/charm/tm_"+Extrapolation_charm_mode+"/fpi_"+V_charm_1_L.Tag[i_ens]+".dat");
      P5P5_charm_pion_distr = Corr.corr_t(pt2_pion.col(0)[j_ens], "");
      P5P5_OS_charm_pion_distr = Corr.corr_t(corr_P5P5_OS.col(0)[j_ens], "");
      A0P5_charm_pion_distr = Corr.corr_t(corr_A0P5.col(0)[j_ens], "");
      A0P5_OS_charm_pion_distr = Corr.corr_t(corr_A0P5_OS.col(0)[j_ens], "");
    }



  
    //L
    V_charm_1_distr_L = Corr.corr_t(V_charm_1_L.col(0)[i_ens], "../data/gm2/charm/tm_"+Extrapolation_charm_mode+"/corr_1_"+V_charm_1_L.Tag[i_ens]+"_L.dat");
    V_charm_2_distr_L = Corr.corr_t(V_charm_2_L.col(0)[i_ens], "../data/gm2/charm/tm_"+Extrapolation_charm_mode+"/corr_2_"+V_charm_2_L.Tag[i_ens]+"_L.dat");
    V_charm_3_distr_L = Corr.corr_t(V_charm_3_L.col(0)[i_ens], "../data/gm2/charm/tm_"+Extrapolation_charm_mode+"/corr_3_"+V_charm_3_L.Tag[i_ens]+"_L.dat");

    Corr.Tmin= Tmin_etaC;
    Corr.Tmax= Tmax_etaC;
    cbar_c_distr_L = Corr.corr_t(pt2_etaC_L.col(0)[i_ens], "");
    cbar_c_distr_M = Corr.corr_t(pt2_etaC_M.col(0)[i_ens], "");
    cbar_c_distr_H = Corr.corr_t(pt2_etaC_H.col(0)[i_ens], "");
    M_etaC_distr_L = Corr.effective_mass_t(cbar_c_distr_L, "../data/gm2/charm/tm_"+Extrapolation_charm_mode+"/M_etaC_"+V_charm_1_L.Tag[i_ens]+"_L.dat");
    M_etaC_distr_M = Corr.effective_mass_t(cbar_c_distr_M, "../data/gm2/charm/tm_"+Extrapolation_charm_mode+"/M_etaC_"+V_charm_1_L.Tag[i_ens]+"_M.dat");
    M_etaC_distr_H = Corr.effective_mass_t(cbar_c_distr_H, "../data/gm2/charm/tm_"+Extrapolation_charm_mode+"/M_etaC_"+V_charm_1_L.Tag[i_ens]+"_H.dat");
    overlap_P5P5_distr_L = Corr.residue_t(cbar_c_distr_L, "");
    overlap_P5P5_distr_M = Corr.residue_t(cbar_c_distr_M, "");
    overlap_P5P5_distr_H = Corr.residue_t(cbar_c_distr_H, "");
    fetac_L_distr = 2.0*L_info.mc_L*Corr.decay_constant_t(cbar_c_distr_L, "");
    fetac_M_distr = 2.0*L_info.mc_M*Corr.decay_constant_t(cbar_c_distr_M, "");
    fetac_H_distr = 2.0*L_info.mc_H*Corr.decay_constant_t(cbar_c_distr_H, "");
    fetac_L = Corr.Fit_distr(fetac_L_distr);
    fetac_M = Corr.Fit_distr(fetac_M_distr);
    fetac_H = Corr.Fit_distr(fetac_H_distr);
    cbar_c_OS_distr_L = Corr.corr_t(pt2_etaC_OS_L.col(0)[i_ens], "");
    cbar_c_OS_distr_M = Corr.corr_t(pt2_etaC_OS_M.col(0)[i_ens], "");
    cbar_c_OS_distr_H = Corr.corr_t(pt2_etaC_OS_H.col(0)[i_ens], "");
    M_etaC_distr_OS_L = Corr.effective_mass_t(cbar_c_OS_distr_L, "../data/gm2/charm/OS_"+Extrapolation_charm_mode+"/M_etaC_"+V_charm_1_L.Tag[i_ens]+"_L.dat");
    M_etaC_distr_OS_M = Corr.effective_mass_t(cbar_c_OS_distr_M, "../data/gm2/charm/OS_"+Extrapolation_charm_mode+"/M_etaC_"+V_charm_1_M.Tag[i_ens]+"_M.dat");
    M_etaC_distr_OS_H = Corr.effective_mass_t(cbar_c_OS_distr_H, "../data/gm2/charm/OS_"+Extrapolation_charm_mode+"/M_etaC_"+V_charm_1_H.Tag[i_ens]+"_H.dat");
    overlap_P5P5_OS_distr_L = Corr.residue_t(cbar_c_OS_distr_L, "");
    overlap_P5P5_OS_distr_M = Corr.residue_t(cbar_c_OS_distr_M, "");
    overlap_P5P5_OS_distr_H = Corr.residue_t(cbar_c_OS_distr_H, "");
    Corr.Tmin = Tmin_P5P5;
    Corr.Tmax = Tmax_P5P5;
    overlap_P5P5_distr_light = Corr.residue_t(P5P5_charm_pion_distr, "");
    overlap_P5P5_OS_distr_light = Corr.residue_t(P5P5_OS_charm_pion_distr, "");
    //fit pion
    Mpi = Corr.Fit_distr(Mpi_distr);
    fpi = Corr.Fit_distr(fpi_distr);
    Mpi_OS = Corr.Fit_distr(Mpi_OS_distr);
  
    //M
    V_charm_1_distr_M = Corr.corr_t(V_charm_1_M.col(0)[i_ens], "../data/gm2/charm/tm_"+Extrapolation_charm_mode+"/corr_1_"+V_charm_1_M.Tag[i_ens]+"_M.dat");
    V_charm_2_distr_M = Corr.corr_t(V_charm_2_M.col(0)[i_ens], "../data/gm2/charm/tm_"+Extrapolation_charm_mode+"/corr_2_"+V_charm_2_M.Tag[i_ens]+"_M.dat");
    V_charm_3_distr_M = Corr.corr_t(V_charm_3_M.col(0)[i_ens], "../data/gm2/charm/tm_"+Extrapolation_charm_mode+"/corr_3_"+V_charm_3_M.Tag[i_ens]+"_M.dat");
    //H
    V_charm_1_distr_H = Corr.corr_t(V_charm_1_H.col(0)[i_ens], "../data/gm2/charm/tm_"+Extrapolation_charm_mode+"/corr_1_"+V_charm_1_H.Tag[i_ens]+"_H.dat");
    V_charm_2_distr_H = Corr.corr_t(V_charm_2_H.col(0)[i_ens], "../data/gm2/charm/tm_"+Extrapolation_charm_mode+"/corr_2_"+V_charm_2_H.Tag[i_ens]+"_H.dat");
    V_charm_3_distr_H = Corr.corr_t(V_charm_3_H.col(0)[i_ens], "../data/gm2/charm/tm_"+Extrapolation_charm_mode+"/corr_3_"+V_charm_3_H.Tag[i_ens]+"_H.dat");
  
    //OS
    //L
    V_charm_OS_1_distr_L = Corr.corr_t(V_charm_OS_1_L.col(0)[i_ens], "../data/gm2/charm/OS_"+Extrapolation_charm_mode+"/corr_1_"+V_charm_1_L.Tag[i_ens]+"_L.dat");
    V_charm_OS_2_distr_L = Corr.corr_t(V_charm_OS_2_L.col(0)[i_ens], "../data/gm2/charm/OS_"+Extrapolation_charm_mode+"/corr_2_"+V_charm_2_L.Tag[i_ens]+"_L.dat");
    V_charm_OS_3_distr_L = Corr.corr_t(V_charm_OS_3_L.col(0)[i_ens], "../data/gm2/charm/OS_"+Extrapolation_charm_mode+"/corr_3_"+V_charm_3_L.Tag[i_ens]+"_L.dat");
    //M
    V_charm_OS_1_distr_M = Corr.corr_t(V_charm_OS_1_M.col(0)[i_ens], "../data/gm2/charm/OS_"+Extrapolation_charm_mode+"/corr_1_"+V_charm_1_M.Tag[i_ens]+"_M.dat");
    V_charm_OS_2_distr_M = Corr.corr_t(V_charm_OS_2_M.col(0)[i_ens], "../data/gm2/charm/OS_"+Extrapolation_charm_mode+"/corr_2_"+V_charm_2_M.Tag[i_ens]+"_M.dat");
    V_charm_OS_3_distr_M = Corr.corr_t(V_charm_OS_3_M.col(0)[i_ens], "../data/gm2/charm/OS_"+Extrapolation_charm_mode+"/corr_3_"+V_charm_3_M.Tag[i_ens]+"_M.dat");
    //H
    V_charm_OS_1_distr_H = Corr.corr_t(V_charm_OS_1_H.col(0)[i_ens], "../data/gm2/charm/OS_"+Extrapolation_charm_mode+"/corr_1_"+V_charm_1_H.Tag[i_ens]+"_H.dat");
    V_charm_OS_2_distr_H = Corr.corr_t(V_charm_OS_2_H.col(0)[i_ens], "../data/gm2/charm/OS_"+Extrapolation_charm_mode+"/corr_2_"+V_charm_2_H.Tag[i_ens]+"_H.dat");
    V_charm_OS_3_distr_H = Corr.corr_t(V_charm_OS_3_H.col(0)[i_ens], "../data/gm2/charm/OS_"+Extrapolation_charm_mode+"/corr_3_"+V_charm_3_H.Tag[i_ens]+"_H.dat");


   
    //sum over the Lorenz indices of the e.m. current
    //tm
    //L
    V_charm_distr_L= (pow(qc,2)/3.0)*(V_charm_1_distr_L+ V_charm_2_distr_L + V_charm_3_distr_L);
    //M
    V_charm_distr_M= (pow(qc,2)/3.0)*(V_charm_1_distr_M+ V_charm_2_distr_M + V_charm_3_distr_M);
    //H
    V_charm_distr_H= (pow(qc,2)/3.0)*(V_charm_1_distr_H+ V_charm_2_distr_H + V_charm_3_distr_H);
  
    //OS
    //L
    V_charm_OS_distr_L= (pow(qc,2)/3.0)*(V_charm_OS_1_distr_L+ V_charm_OS_2_distr_L + V_charm_OS_3_distr_L);
    //M
    V_charm_OS_distr_M= (pow(qc,2)/3.0)*(V_charm_OS_1_distr_M+ V_charm_OS_2_distr_M + V_charm_OS_3_distr_M);
    //H
    V_charm_OS_distr_H= (pow(qc,2)/3.0)*(V_charm_OS_1_distr_H+ V_charm_OS_2_distr_H + V_charm_OS_3_distr_H);



    bool Found_disco_ens=false;
    bool Found_disco_impr_ens=false;
    if(Include_charm_disco) {
      int i_ens_disco=0;
      for(int j=0;j<Nens_disco_charm;j++) if(disco_charm.Tag[j] == V_charm_1_L.Tag[i_ens]) { Found_disco_ens=true; i_ens_disco=j;disco_charm_Tags.push_back(disco_charm.Tag[j]) ;break;}
      if(Found_disco_ens) {
	disco_distr = Corr.corr_t(disco_charm.col(0)[i_ens_disco], "");
	disco_distr = disco_distr*(pow(qc,2));
      }

      //improved
      for(int j=0;j<disco_impr_charm.size;j++) if(disco_impr_charm.Tag[j] == V_charm_1_L.Tag[i_ens]) { Found_disco_impr_ens=true; i_ens_disco=j;disco_impr_charm_Tags.push_back(disco_impr_charm.Tag[j]) ;break;}
      if(Found_disco_impr_ens) {
	disco_impr_distr = Corr.corr_t(disco_impr_charm.col(0)[i_ens_disco], "");
	disco_impr_distr = disco_impr_distr*(pow(qc,2));
      }
      
    }


    //off-diagonal light-charm
    //improved
    bool Found_disco_impr_light_charm_ens=false;
    if(Include_off_diagonal_disco) {

      int i_ens_disco=0;
      for(int j=0;j<disco_impr_light_charm.size;j++) if(disco_impr_light_charm.Tag[j] == V_charm_1_L.Tag[i_ens]) { Found_disco_impr_light_charm_ens=true; i_ens_disco=j;disco_impr_light_charm_Tags.push_back(disco_impr_light_charm.Tag[j]) ;break;}
      if(Found_disco_impr_light_charm_ens) {
	disco_impr_light_charm_distr = Corr.corr_t(disco_impr_light_charm.col(0)[i_ens_disco], "");
	disco_impr_light_charm_distr = disco_impr_light_charm_distr*(2.0*qc*(qu+qd));
      }
    }

    //improved D
    bool Found_disco_impr_lightD_charm_ens=false;
    if(Include_off_diagonal_disco) {

      int i_ens_disco=0;
      for(int j=0;j<disco_impr_lightD_charm.size;j++) if(disco_impr_lightD_charm.Tag[j] == V_charm_1_L.Tag[i_ens]) { Found_disco_impr_lightD_charm_ens=true; i_ens_disco=j;disco_impr_lightD_charm_Tags.push_back(disco_impr_lightD_charm.Tag[j]) ;break;}
      if(Found_disco_impr_lightD_charm_ens) {
	disco_impr_lightD_charm_distr = Corr.corr_t(disco_impr_lightD_charm.col(0)[i_ens_disco], "");
	disco_impr_lightD_charm_distr = disco_impr_lightD_charm_distr*(2.0*qc*(qu+qd));
      }
    }


   

  
  
    //#######################################  COMPUTATION OF ZV AND ZA (Hadronic method) ##################################
    //define lambda functions to be used
    auto sqr= [=](double a, double b) {return sqrt(a);};
    auto SINH= [](double m) -> double  {return sinh(m);};
    auto SINH2 = [](double m, double t) -> double {return sinh(m);};

  
    //take ratio between OS and tm pion amplitude to compute Zp/Zs RC.
    ratio_P5P5_overlap_OS_tm_L= overlap_P5P5_OS_distr_L/overlap_P5P5_distr_L;
    Zp_ov_Zs_distr_L = distr_t_list::f_of_distr_list(sqr, ratio_P5P5_overlap_OS_tm_L);
    ratio_P5P5_overlap_OS_tm_M= overlap_P5P5_OS_distr_M/overlap_P5P5_distr_M;
    Zp_ov_Zs_distr_M = distr_t_list::f_of_distr_list(sqr, ratio_P5P5_overlap_OS_tm_M);
    ratio_P5P5_overlap_OS_tm_H= overlap_P5P5_OS_distr_H/overlap_P5P5_distr_H;
    Zp_ov_Zs_distr_H = distr_t_list::f_of_distr_list(sqr, ratio_P5P5_overlap_OS_tm_H);

    //for light Zs
    ratio_P5P5_overlap_OS_tm_light = overlap_P5P5_OS_distr_light/overlap_P5P5_distr_light;
    Zp_ov_Zs_distr_light = distr_t_list::f_of_distr_list(sqr, ratio_P5P5_overlap_OS_tm_light);



    //antysymmetrize w.r.t. t -> T-t for A0P5 correlators
    if(V_charm_1_L.Tag[i_ens].substr(1,1) != "D") Corr.Reflection_sign = -1;
    A0P5_distr_L= Corr.corr_t(corr_A0P5_charm_L.col(0)[i_ens], "../data/gm2/charm/tm_"+Extrapolation_charm_mode+"/A0P5_"+V_charm_1_L.Tag[i_ens]+"_L.dat");
    A0P5_OS_distr_L = Corr.corr_t(corr_A0P5_OS_charm_L.col(0)[i_ens], "../data/gm2/charm/OS_"+Extrapolation_charm_mode+"/A0P5_"+V_charm_1_L.Tag[i_ens]+"_L.dat");
    A0P5_distr_M= Corr.corr_t(corr_A0P5_charm_M.col(0)[i_ens], "../data/gm2/charm/tm_"+Extrapolation_charm_mode+"/A0P5_"+V_charm_1_L.Tag[i_ens]+"_M.dat");
    A0P5_OS_distr_M = Corr.corr_t(corr_A0P5_OS_charm_M.col(0)[i_ens], "../data/gm2/charm/OS_"+Extrapolation_charm_mode+"/A0P5_"+V_charm_1_L.Tag[i_ens]+"_M.dat");
    A0P5_distr_H= Corr.corr_t(corr_A0P5_charm_H.col(0)[i_ens], "../data/gm2/charm/tm_"+Extrapolation_charm_mode+"/A0P5_"+V_charm_1_L.Tag[i_ens]+"_H.dat");
    A0P5_OS_distr_H = Corr.corr_t(corr_A0P5_OS_charm_H.col(0)[i_ens], "../data/gm2/charm/OS_"+Extrapolation_charm_mode+"/A0P5_"+V_charm_1_L.Tag[i_ens]+"_H.dat");

  
    //restore symmetrization
    Corr.Reflection_sign = 1;

    //compute RV (estimator for Zv)


    RV_L= 2.0*L_info.mc_L*cbar_c_distr_L/distr_t_list::derivative(A0P5_distr_L, 0); //central derivative
    RV_M= 2.0*L_info.mc_M*cbar_c_distr_M/distr_t_list::derivative(A0P5_distr_M, 0); //central derivative
    RV_H= 2.0*L_info.mc_H*cbar_c_distr_H/distr_t_list::derivative(A0P5_distr_H, 0); //central derivative

    RV_light = 2.0*L_info.ml*P5P5_charm_pion_distr/distr_t_list::derivative(A0P5_charm_pion_distr,0); //central derivative 
  
    //tm and OS P5P5
    Corr.Tmin=Tmin_etaC;
    Corr.Tmax=Tmax_etaC;
    cbar_c_mass_L = Corr.Fit_distr(M_etaC_distr_L);
    cbar_c_OS_mass_L= Corr.Fit_distr(M_etaC_distr_OS_L);
    cbar_c_mass_M = Corr.Fit_distr(M_etaC_distr_M);
    cbar_c_OS_mass_M= Corr.Fit_distr(M_etaC_distr_OS_M);
    cbar_c_mass_H = Corr.Fit_distr(M_etaC_distr_H);
    cbar_c_OS_mass_H= Corr.Fit_distr(M_etaC_distr_OS_H);

    if( V_charm_1_L.Tag[i_ens].substr(1,1) != "A") {

      int id=0;
      if( V_charm_1_L.Tag[i_ens].substr(1,1) == "D") id = 3;
      else if(V_charm_1_L.Tag[i_ens].substr(1,1) == "C") id=2;
      else if(V_charm_1_L.Tag[i_ens].substr(1,1) == "B") id=0;
      else crash("Ensemble not recognised");
      Mp_c1_tm.distr_list[id] = cbar_c_mass_L;
      Mp_c2_tm.distr_list[id] =cbar_c_mass_M;
      Mp_c3_tm.distr_list[id] = cbar_c_mass_H;
      Mp_c1_OS.distr_list[id] =cbar_c_OS_mass_L;
      Mp_c2_OS.distr_list[id] =cbar_c_OS_mass_M;
      Mp_c3_OS.distr_list[id] =cbar_c_OS_mass_H;
      if (id == 0) {
	Mp_c1_tm.distr_list[id+1] = cbar_c_mass_L;
	Mp_c2_tm.distr_list[id+1] =cbar_c_mass_M;
	Mp_c3_tm.distr_list[id+1] = cbar_c_mass_H;
	Mp_c1_OS.distr_list[id+1] =cbar_c_OS_mass_L;
	Mp_c2_OS.distr_list[id+1] =cbar_c_OS_mass_M;
	Mp_c3_OS.distr_list[id+1] =cbar_c_OS_mass_H;
      }
    }
    //fit obs to compute Zv and Za (hadronic method)
    Zp_ov_Zs_L = Corr.Fit_distr(Zp_ov_Zs_distr_L);
    Zp_ov_Zs_M = Corr.Fit_distr(Zp_ov_Zs_distr_M);
    Zp_ov_Zs_H = Corr.Fit_distr(Zp_ov_Zs_distr_H);
    Corr.Tmin= Tmin_P5P5;
    Corr.Tmax= Tmax_P5P5;
    Zp_ov_Zs_light= Corr.Fit_distr(Zp_ov_Zs_distr_light);
    //set plateaux for RV
    int Tmin_RV, Tmax_RV;
    if(V_charm_1_L.Tag[i_ens].substr(1,1) == "C") {
      if(V_charm_1_L.Tag[i_ens]=="cC211a.06.80") { Tmin_RV=36; Tmax_RV=73;}
      else crash("Cannot find ensemble tag: "+V_charm_1_L.Tag[i_ens]);
    }
    else if(V_charm_1_L.Tag[i_ens].substr(1,1) == "B") {
      if(V_charm_1_L.Tag[i_ens]== "cB211a.14.64") {Tmin_RV=28; Tmax_RV=57;}
      else if(V_charm_1_L.Tag[i_ens] == "cB211a.25.48") {Tmin_RV=34;Tmax_RV=46;}
      else if(V_charm_1_L.Tag[i_ens] == "cB211b.072.64") {Tmin_RV=35; Tmax_RV= 61;}
      else if(V_charm_1_L.Tag[i_ens] == "cB211b.072.96") {Tmin_RV=35; Tmax_RV= 61;}
      else crash("Cannot find ensemble tag: "+V_charm_1_L.Tag[i_ens]);
    }
    else if(V_charm_1_L.Tag[i_ens].substr(1,1) == "A") {
      if(V_charm_1_L.Tag[i_ens] == "cA211a.12.48") {Tmin_RV=28; Tmax_RV=46;}
      else if(V_charm_1_L.Tag[i_ens] == "cA211a.40.24") {Tmin_RV=17; Tmax_RV=23;}
      else if(V_charm_1_L.Tag[i_ens] == "cA211a.53.24") {Tmin_RV=16; Tmax_RV=23;}
      else if(V_charm_1_L.Tag[i_ens] == "cA211ab.30.32") {Tmin_RV=21;Tmax_RV=31;}
      else crash("Cannot find ensemble tag: "+V_charm_1_L.Tag[i_ens]);
    }
    else if(V_charm_1_L.Tag[i_ens].substr(1,1) == "D") {
      if(V_charm_1_L.Tag[i_ens] == "cD211a.054.96") {Tmin_RV=44; Tmax_RV=70;}
      else crash("Cannot find ensemble tag: "+V_charm_1_L.Tag[i_ens]);
    }
    else crash("Ensemble tag not valid");

    Corr.Tmin = Tmin_RV;
    Corr.Tmax = Tmax_RV;
  

    Zv_hadr_L= Corr.Fit_distr(RV_L);
    Zv_hadr_M= Corr.Fit_distr(RV_M);
    Zv_hadr_H= Corr.Fit_distr(RV_H);

    int Tmin_RV_light, Tmax_RV_light;

    //set plateaux for RV_light
    if(V_charm_1_L.Tag[i_ens].substr(1,1) == "C") {
      if(V_charm_1_L.Tag[i_ens]=="cC211a.06.80") { Tmin_RV_light=36; Tmax_RV_light=73;}
      else crash("Cannot find ensemble tag: "+V_charm_1_L.Tag[i_ens]);
    }
    else if(V_charm_1_L.Tag[i_ens].substr(1,1) == "B") {
      if(V_charm_1_L.Tag[i_ens]== "cB211a.14.64") {Tmin_RV_light=28; Tmax_RV_light=57;}
      else if(V_charm_1_L.Tag[i_ens] == "cB211a.25.48") {Tmin_RV_light=34;Tmax_RV_light=46;}
      else if(V_charm_1_L.Tag[i_ens] == "cB211b.072.64") {Tmin_RV_light=35; Tmax_RV_light= 61;}
      else if(V_charm_1_L.Tag[i_ens] == "cB211b.072.96") {Tmin_RV_light=35; Tmax_RV_light= 61;}
      else crash("Cannot find ensemble tag: "+V_charm_1_L.Tag[i_ens]);
    }
    else if(V_charm_1_L.Tag[i_ens].substr(1,1) == "A") {
      if(V_charm_1_L.Tag[i_ens] == "cA211a.12.48") {Tmin_RV_light=21; Tmax_RV_light=46;}
      else if(V_charm_1_L.Tag[i_ens] == "cA211a.40.24") {Tmin_RV_light=10; Tmax_RV_light=23;}
      else if(V_charm_1_L.Tag[i_ens] == "cA211a.53.24") {Tmin_RV_light=10; Tmax_RV_light=23;}
      else if(V_charm_1_L.Tag[i_ens] == "cA211ab.30.32") {Tmin_RV_light=12;Tmax_RV_light=30;}
      else crash("Cannot find ensemble tag: "+V_charm_1_L.Tag[i_ens]);
    }
    else if(V_charm_1_L.Tag[i_ens].substr(1,1) == "D") {
      if(V_charm_1_L.Tag[i_ens] == "cD211a.054.96") {Tmin_RV_light=44; Tmax_RV_light=70;}
      else crash("Cannot find ensemble tag: "+V_charm_1_L.Tag[i_ens]);
    }
    else crash("Ensemble tag not valid");

    Corr.Tmin= Tmin_RV_light;
    Corr.Tmax= Tmax_RV_light;
  
    Zv_hadr_light = Corr.Fit_distr(RV_light);
  
    RA_L = 2.0*L_info.mc_L*(cbar_c_OS_distr_L/distr_t_list::derivative(A0P5_OS_distr_L, 0))*(cbar_c_OS_mass_L/cbar_c_mass_L)*(distr_t::f_of_distr(SINH, cbar_c_OS_mass_L)/distr_t::f_of_distr(SINH, cbar_c_mass_L))*(1.0/Zp_ov_Zs_L);
    RA_M = 2.0*L_info.mc_M*(cbar_c_OS_distr_M/distr_t_list::derivative(A0P5_OS_distr_M, 0))*(cbar_c_OS_mass_M/cbar_c_mass_M)*(distr_t::f_of_distr(SINH, cbar_c_OS_mass_M)/distr_t::f_of_distr(SINH, cbar_c_mass_M))*(1.0/Zp_ov_Zs_M);
    RA_H = 2.0*L_info.mc_H*(cbar_c_OS_distr_H/distr_t_list::derivative(A0P5_OS_distr_H, 0))*(cbar_c_OS_mass_H/cbar_c_mass_H)*(distr_t::f_of_distr(SINH, cbar_c_OS_mass_H)/distr_t::f_of_distr(SINH, cbar_c_mass_H))*(1.0/Zp_ov_Zs_H);

    //RA_light
    RA0_light = 2.0*L_info.ml*(P5P5_OS_charm_pion_distr/distr_t_list::derivative(A0P5_OS_charm_pion_distr,0))*(Mpi_OS_distr/Mpi_distr)*(1.0/Zp_ov_Zs_distr_light)*distr_t_list::f_of_distr_list(SINH2, Mpi_OS_distr)/distr_t_list::f_of_distr_list(SINH2, Mpi_distr);
    RA_light = 2.0*L_info.ml*(P5P5_OS_charm_pion_distr/distr_t_list::derivative(A0P5_OS_charm_pion_distr,0))*(Mpi_OS/Mpi)*(distr_t::f_of_distr(SINH, Mpi_OS)/distr_t::f_of_distr(SINH, Mpi))*(1.0/Zp_ov_Zs_light);
    //set plateaux for RA
    int Tmin_RA_light=0;
    int Tmax_RA_light=0;
    //set time intervals for RA
    int Tmin_RA= Tmin_etaC;
    int Tmax_RA= Tmax_etaC;
  
    if(V_charm_1_L.Tag[i_ens].substr(1,1) == "C") {
      if(V_charm_1_L.Tag[i_ens]=="cC211a.06.80") { Tmin_RA_light=25; Tmax_RA_light=50;}
      else crash("Cannot find ensemble tag: "+V_charm_1_L.Tag[i_ens]);
    }
    else if(V_charm_1_L.Tag[i_ens].substr(1,1) == "B") {
      if(V_charm_1_L.Tag[i_ens]== "cB211a.14.64") {Tmin_RA_light=18; Tmax_RA_light=26;}
      else if(V_charm_1_L.Tag[i_ens] == "cB211a.25.48") {Tmin_RA_light=10;Tmax_RA_light=21;}
      else if(V_charm_1_L.Tag[i_ens] == "cB211b.072.64") {Tmin_RA_light=20; Tmax_RA_light= 50;}
      else if(V_charm_1_L.Tag[i_ens] == "cB211b.072.96") {Tmin_RA_light=20; Tmax_RA_light= 50;}
      else crash("Cannot find ensemble tag: "+V_charm_1_L.Tag[i_ens]);
    }
    else if(V_charm_1_L.Tag[i_ens].substr(1,1) == "A") {
      if(V_charm_1_L.Tag[i_ens] == "cA211a.12.48") {Tmin_RA_light=9; Tmax_RA_light=20;}
      else if(V_charm_1_L.Tag[i_ens] == "cA211a.40.24") {Tmin_RA_light=11; Tmax_RA_light=20; Tmin_RA = 17; Tmax_RA = 23 ;}
      else if(V_charm_1_L.Tag[i_ens] == "cA211a.53.24") {Tmin_RA_light=10; Tmax_RA_light=20; Tmin_RA = 17; Tmax_RA = 23;}
      else if(V_charm_1_L.Tag[i_ens] == "cA211ab.30.32") {Tmin_RA_light=14;Tmax_RA_light=27; Tmin_RA = 19; Tmax_RA = 30;}
      else crash("Cannot find ensemble tag: "+V_charm_1_L.Tag[i_ens]);
    }
    else if(V_charm_1_L.Tag[i_ens].substr(1,1) == "D") {
      if(V_charm_1_L.Tag[i_ens] == "cD211a.054.96") {Tmin_RA_light=32; Tmax_RA_light=70;}
      else crash("Cannot find ensemble tag: "+V_charm_1_L.Tag[i_ens]);
    }
    else crash("Ensemble tag not valid");
  


    Corr.Tmin =Tmin_RA;
    Corr.Tmax =Tmax_RA;
  
    Za_hadr_L= Corr.Fit_distr(RA_L);
    Za_hadr_M= Corr.Fit_distr(RA_M);
    Za_hadr_H= Corr.Fit_distr(RA_H);

    Corr.Tmin=Tmin_RA_light;
    Corr.Tmax=Tmax_RA_light;

    //fit light
    Za_hadr_light = Corr.Fit_distr(RA_light);

    Corr.Tmin= Tmin_etaC;
    Corr.Tmax= Tmax_etaC;

    //print Rv and RA
    //print RV
    Print_To_File({}, {RV_L.ave(), RV_L.err(), RV_M.ave(), RV_M.err(), RV_H.ave(), RV_H.err(), RV_light.ave(), RV_light.err()}, "../data/gm2/charm/Z_"+Extrapolation_charm_mode+"/RV_"+V_charm_1_L.Tag[i_ens]+".dat.t", "", "");
    //print RA
    Print_To_File({}, {RA_L.ave(), RA_L.err(), RA_M.ave(), RA_M.err(), RA_H.ave(), RA_H.err(), RA_light.ave(), RA_light.err(), RA0_light.ave(), RA0_light.err()}, "../data/gm2/charm/Z_"+Extrapolation_charm_mode+"/RA_"+V_charm_1_L.Tag[i_ens]+".dat.t", "", "");
    //print Zp_ov_Zs
    Print_To_File({}, {Zp_ov_Zs_distr_L.ave(), Zp_ov_Zs_distr_L.err(), Zp_ov_Zs_distr_M.ave(), Zp_ov_Zs_distr_M.err(), Zp_ov_Zs_distr_H.ave(), Zp_ov_Zs_distr_H.err(), Zp_ov_Zs_distr_light.ave(), Zp_ov_Zs_distr_light.err()}, "../data/gm2/charm/Z_"+Extrapolation_charm_mode+"/Zp_ov_Zs_"+V_charm_1_L.Tag[i_ens]+".dat.t", "", "");

    //push_back

    Zv_fit_charm_L.distr_list.push_back(Zv_hadr_L);
    Za_fit_charm_L.distr_list.push_back(Za_hadr_L);
    Zv_fit_charm_M.distr_list.push_back(Zv_hadr_M);
    Za_fit_charm_M.distr_list.push_back(Za_hadr_M);
    Zv_fit_charm_H.distr_list.push_back(Zv_hadr_H);
    Za_fit_charm_H.distr_list.push_back(Za_hadr_H);

    Za_fit_charm_light.distr_list.push_back(Za_hadr_light);
    Zv_fit_charm_light.distr_list.push_back(Zv_hadr_light);

 

    if(Use_Za_Zv_from_charm_run) { Zv_L= Zv_hadr_L; Za_L=Za_hadr_L; Zv_M = Zv_hadr_M; Za_M = Za_hadr_M; Zv_H = Zv_hadr_H; Za_H = Za_hadr_H;}


    //################################################ END OF COMPUTATION OF ZV AND ZA (Hadronic method) #####################################à


    //extract effective masses, overlap from V and fit



    
  


    //tm

    //L
    Corr.Tmin=Tmin_P5P5;
    Corr.Tmax=Tmax_P5P5;
 
    Corr.Tmin=Tmin_etaC;
    Corr.Tmax=Tmax_etaC;
    M_etaC_L = Corr.Fit_distr(M_etaC_distr_L);
    Corr.Tmin= Tmin_VV;
    Corr.Tmax= Tmax_VV;
    MV_charm_distr_L= Corr.effective_mass_t(V_charm_distr_L, "../data/gm2/charm/tm_"+Extrapolation_charm_mode+"/MV_mass_"+V_charm_1_L.Tag[i_ens]+"_L.dat");
    ZV_charm_distr_L= Corr.residue_t(V_charm_distr_L, "../data/gm2/charm/tm_"+Extrapolation_charm_mode+"/ZV_overlap_"+V_charm_1_L.Tag[i_ens]+"_L.dat");
    MV_charm_L = Corr.Fit_distr(MV_charm_distr_L);
    ZV_charm_L = Corr.Fit_distr(ZV_charm_distr_L);
 
    //M
    Corr.Tmin=Tmin_P5P5;
    Corr.Tmax=Tmax_P5P5;
    Corr.Tmin=Tmin_etaC;
    Corr.Tmax=Tmax_etaC;
    M_etaC_M = Corr.Fit_distr(M_etaC_distr_M);
    Corr.Tmin=Tmin_VV;
    Corr.Tmax=Tmax_VV;
    MV_charm_distr_M= Corr.effective_mass_t(V_charm_distr_M, "../data/gm2/charm/tm_"+Extrapolation_charm_mode+"/MV_mass_"+V_charm_1_M.Tag[i_ens]+"_M.dat");
    ZV_charm_distr_M= Corr.residue_t(V_charm_distr_M, "../data/gm2/charm/tm_"+Extrapolation_charm_mode+"/ZV_overlap_"+V_charm_1_M.Tag[i_ens]+"_M.dat");
    MV_charm_M = Corr.Fit_distr(MV_charm_distr_M);
    ZV_charm_M = Corr.Fit_distr(ZV_charm_distr_M);
 
    //H
    Corr.Tmin=Tmin_P5P5;
    Corr.Tmax=Tmax_P5P5;
    Corr.Tmin=Tmin_etaC;
    Corr.Tmax=Tmax_etaC;
    M_etaC_H = Corr.Fit_distr(M_etaC_distr_H);
    Corr.Tmin=Tmin_VV;
    Corr.Tmax=Tmax_VV;
    MV_charm_distr_H= Corr.effective_mass_t(V_charm_distr_H, "../data/gm2/charm/tm_"+Extrapolation_charm_mode+"/MV_mass_"+V_charm_1_H.Tag[i_ens]+"_H.dat");
    ZV_charm_distr_H= Corr.residue_t(V_charm_distr_H, "../data/gm2/charm/tm_"+Extrapolation_charm_mode+"/ZV_overlap_"+V_charm_1_H.Tag[i_ens]+"_H.dat");
    MV_charm_H = Corr.Fit_distr(MV_charm_distr_H);
    ZV_charm_H = Corr.Fit_distr(ZV_charm_distr_H);

    cout<<"Printing etac masses for Ensemble: "<<V_charm_1_L.Tag[i_ens]<<endl;
    cout.precision(10);
    cout<<" M_etac(L): "<<M_etaC_L.ave()<<" +- "<<M_etaC_L.err()<<endl;
    cout<<" M_etac(M): "<<M_etaC_M.ave()<<" +- "<<M_etaC_M.err()<<endl;
    cout<<" M_etac(H): "<<M_etaC_H.ave()<<" +- "<<M_etaC_H.err()<<endl;
    cout<<"#######################################"<<endl;
    cout.precision(10);
    cout<<"Printing Jpsi masses for Ensemble: "<<V_charm_1_L.Tag[i_ens]<<endl;
    cout<<" M_Jpsi(L): "<<MV_charm_L.ave()<<" +- "<<MV_charm_L.err()<<endl;
    cout<<" M_Jpsi(M): "<<MV_charm_M.ave()<<" +- "<<MV_charm_M.err()<<endl;
    cout<<" M_Jpsi(H): "<<MV_charm_H.ave()<<" +- "<<MV_charm_H.err()<<endl;
    cout<<"#######################################"<<endl;



    //Print to File M_etac and M_Jpsi

    distr_t_list etac_masses_list(UseJack);
    distr_t_list Jpsi_masses_list(UseJack);

    etac_masses_list.distr_list.push_back( M_etaC_L);
    etac_masses_list.distr_list.push_back( M_etaC_M);
    etac_masses_list.distr_list.push_back( M_etaC_H);
    Jpsi_masses_list.distr_list.push_back( MV_charm_L);
    Jpsi_masses_list.distr_list.push_back( MV_charm_M);
    Jpsi_masses_list.distr_list.push_back( MV_charm_H);

    Print_To_File({}, {etac_masses_list.ave(), etac_masses_list.err(), Jpsi_masses_list.ave(), Jpsi_masses_list.err()}, "../data/gm2/charm/masses_"+V_charm_1_L.Tag[i_ens]+".list", "", "# etac   Jpsi");



    //compute effective mass of J/psi including disconnected diagram
    if(Found_disco_impr_ens) {

    //tm  
    distr_t_list MV_L_full_eff_distr= Corr.effective_mass_t( Za_L*Za_L*V_charm_distr_L + Zv_L*Zv_L*disco_impr_distr, "../data/gm2/charm/tm_"+Extrapolation_charm_mode+"/MV_mass_w_disco_"+V_charm_1_L.Tag[i_ens]+"_L.dat");
    distr_t_list MV_M_full_eff_distr= Corr.effective_mass_t( Za_M*Za_M*V_charm_distr_M + Zv_M*Zv_M*disco_impr_distr, "../data/gm2/charm/tm_"+Extrapolation_charm_mode+"/MV_mass_w_disco_"+V_charm_1_L.Tag[i_ens]+"_M.dat");
    distr_t_list MV_H_full_eff_distr= Corr.effective_mass_t( Za_H*Za_H*V_charm_distr_H + Zv_H*Zv_H*disco_impr_distr, "../data/gm2/charm/tm_"+Extrapolation_charm_mode+"/MV_mass_w_disco_"+V_charm_1_L.Tag[i_ens]+"_H.dat");
    distr_t_list MV_L_disco_eff_distr = Corr.effective_mass_t( Zv_L*Zv_L*disco_impr_distr, "../data/gm2/charm/tm_"+Extrapolation_charm_mode+"/MV_mass_only_disco_"+V_charm_1_L.Tag[i_ens]+"_L.dat") ;
    distr_t_list MV_M_disco_eff_distr = Corr.effective_mass_t( Zv_M*Zv_M*disco_impr_distr, "../data/gm2/charm/tm_"+Extrapolation_charm_mode+"/MV_mass_only_disco_"+V_charm_1_L.Tag[i_ens]+"_M.dat");
    distr_t_list MV_H_disco_eff_distr = Corr.effective_mass_t( Zv_H*Zv_H*disco_impr_distr, "../data/gm2/charm/tm_"+Extrapolation_charm_mode+"/MV_mass_only_disco_"+V_charm_1_L.Tag[i_ens]+"_H.dat");
    //OS
    distr_t_list MV_L_full_eff_OS_distr= Corr.effective_mass_t( Zv_L*Zv_L*V_charm_OS_distr_L + Zv_L*Zv_L*disco_impr_distr, "../data/gm2/charm/OS_"+Extrapolation_charm_mode+"/MV_mass_w_disco_"+V_charm_1_L.Tag[i_ens]+"_L.dat");
    distr_t_list MV_M_full_eff_OS_distr= Corr.effective_mass_t( Zv_M*Zv_M*V_charm_OS_distr_M + Zv_M*Zv_M*disco_impr_distr, "../data/gm2/charm/OS_"+Extrapolation_charm_mode+"/MV_mass_w_disco_"+V_charm_1_L.Tag[i_ens]+"_M.dat");
    distr_t_list MV_H_full_eff_OS_distr= Corr.effective_mass_t( Zv_H*Zv_H*V_charm_OS_distr_H + Zv_H*Zv_H*disco_impr_distr, "../data/gm2/charm/OS_"+Extrapolation_charm_mode+"/MV_mass_w_disco_"+V_charm_1_L.Tag[i_ens]+"_H.dat");
    distr_t_list MV_L_disco_eff_OS_distr = Corr.effective_mass_t( Zv_L*Zv_L*disco_impr_distr, "../data/gm2/charm/OS_"+Extrapolation_charm_mode+"/MV_mass_only_disco_"+V_charm_1_L.Tag[i_ens]+"_L.dat") ;
    distr_t_list MV_M_disco_eff_OS_distr = Corr.effective_mass_t( Zv_M*Zv_M*disco_impr_distr, "../data/gm2/charm/OS_"+Extrapolation_charm_mode+"/MV_mass_only_disco_"+V_charm_1_L.Tag[i_ens]+"_M.dat");
    distr_t_list MV_H_disco_eff_OS_distr = Corr.effective_mass_t( Zv_H*Zv_H*disco_impr_distr, "../data/gm2/charm/OS_"+Extrapolation_charm_mode+"/MV_mass_only_disco_"+V_charm_1_L.Tag[i_ens]+"_H.dat");
    
    }
    





    



    //push back to M_etac^2 vector (physical units)
    vector<distr_t> M2_etaC, Jpsi2, M_etaC, Jpsi, mc_list, M_etaC_ov_Jpsi, M_etaC_ov_fpi, Jpsi_ov_fpi, Za_hadr_list, Zv_hadr_list;
    distr_t mc_L_distr, mc_M_distr, mc_H_distr;
    if(V_charm_1_L.Tag[i_ens].substr(1,1) != "D") {
      M2_etaC = { M_etaC_L*M_etaC_L/(a_distr*a_distr), M_etaC_M*M_etaC_M/(a_distr*a_distr), M_etaC_H*M_etaC_H/(a_distr*a_distr)};
      Jpsi2=  { (MV_charm_L*MV_charm_L)/(a_distr*a_distr),  (MV_charm_M*MV_charm_M)/(a_distr*a_distr),  (MV_charm_H*MV_charm_H)/(a_distr*a_distr)};
      M_etaC = { M_etaC_L/a_distr, M_etaC_M/a_distr, M_etaC_H/a_distr};
      Jpsi = { MV_charm_L/a_distr, MV_charm_M/a_distr, MV_charm_H/a_distr};
      M_etaC_ov_Jpsi = {M_etaC_L/MV_charm_L, M_etaC_M/MV_charm_M, M_etaC_H/MV_charm_L};
      M_etaC_ov_fpi = {M_etaC_L/fpi, M_etaC_M/fpi, M_etaC_H/fpi};
      Jpsi_ov_fpi= { MV_charm_L/fpi, MV_charm_M/fpi, MV_charm_L/fpi};
      Za_hadr_list = {Za_hadr_L, Za_hadr_M, Za_hadr_H};
      Zv_hadr_list = {Zv_hadr_L, Zv_hadr_M, Zv_hadr_H};
      for(int ijack=0;ijack<Njacks;ijack++) {
	mc_L_distr.distr.push_back( L_info.mc_L);
	mc_M_distr.distr.push_back( L_info.mc_M);
	mc_H_distr.distr.push_back( L_info.mc_H);
      }
      mc_list = {mc_L_distr, mc_M_distr, mc_H_distr};
    }
    else {
      M2_etaC = { M_etaC_L*M_etaC_L/(a_distr*a_distr), M_etaC_M*M_etaC_M/(a_distr*a_distr)};
      Jpsi2=  { (MV_charm_L*MV_charm_L)/(a_distr*a_distr),  (MV_charm_M*MV_charm_M)/(a_distr*a_distr)};
      M_etaC ={ M_etaC_L/a_distr, M_etaC_M/a_distr};
      Jpsi = {MV_charm_L/a_distr, MV_charm_M/a_distr};
      M_etaC_ov_Jpsi = { M_etaC_L/MV_charm_L, M_etaC_M/MV_charm_M};
      M_etaC_ov_fpi = {M_etaC_L/fpi, M_etaC_M/fpi};
      Jpsi_ov_fpi= { MV_charm_L/fpi, MV_charm_M/fpi};
      Za_hadr_list = {Za_hadr_L, Za_hadr_M};
      Zv_hadr_list = {Zv_hadr_L, Zv_hadr_M};
      for(int ijack=0;ijack<Njacks;ijack++) {
	mc_L_distr.distr.push_back( L_info.mc_L);
	mc_M_distr.distr.push_back( L_info.mc_M);
      }
      mc_list ={mc_L_distr, mc_M_distr};
    }


    distr_t m_Jpsi_phys_distr(UseJack), m_etac_phys_distr(UseJack);
    vector<distr_t> X_2_fit; 
    distr_t X_2_phys;
    vector<distr_t> X_4_fit;
    distr_t X_4_phys;

    //generate fake jackknife distribution for etac and Jpsi
    for(int ijack=0;ijack<Njacks;ijack++) { m_Jpsi_phys_distr.distr.push_back( m_Jpsi + GM()*m_Jpsi_err/sqrt( Njacks -1.0)); m_etac_phys_distr.distr.push_back( m_etac + GM()*m_etac_err/sqrt(Njacks -1.0));}

    if(Extrapolation_charm_mode == "Jpsi") {  X_2_fit = Jpsi; X_2_phys = m_Jpsi_phys_distr; X_4_fit = Jpsi2; X_4_phys = m_Jpsi_phys_distr*m_Jpsi_phys_distr; }
    else if(Extrapolation_charm_mode== "etac") { X_2_fit = M_etaC; X_2_phys = m_etac_phys_distr; X_4_fit = M2_etaC; X_4_phys = m_etac_phys_distr*m_etac_phys_distr;}
    //else if(Extrapolation_charm_mode=="etac_ov_Jpsi") {X_2_fit = M_etaC_ov_Jpsi; X_2_phys = m_etac/m_Jpsi;}
    //else if(Extrapolation_charm_mode=="etac_ov_fpi") { X_2_fit = M_etaC_ov_fpi; X_2_phys = m_etac/fp_phys;}
    //else if(Extrapolation_charm_mode=="Jpsi_ov_fpi") { X_2_fit = Jpsi_ov_fpi; X_2_phys = m_Jpsi/fp_phys;}
    else crash("Extrapolation charm mode: "+Extrapolation_charm_mode+" not yet implemented");

    string tag_mc_extrapolation= Extrapolation_charm_mode;
    distr_t mc_phys_Jpsi_extr = Obs_extrapolation_meson_mass(mc_list, Jpsi, m_Jpsi_phys_distr,  "../data/gm2/charm", "mc_extrapolation_Jpsi_"+V_charm_1_L.Tag[i_ens], UseJack, "SPLINE");
    distr_t mc_phys_etac_extr = Obs_extrapolation_meson_mass(mc_list, M_etaC, m_etac_phys_distr,  "../data/gm2/charm", "mc_extrapolation_etac_"+V_charm_1_L.Tag[i_ens], UseJack, "SPLINE");
    distr_t mc_phys_extr;
    if(Extrapolation_charm_mode=="Jpsi") mc_phys_extr = mc_phys_Jpsi_extr;
    else mc_phys_extr = mc_phys_etac_extr;
    
    //push_back
    mc_extr_list.distr_list.push_back(mc_phys_extr/a_distr);
    mc_extr_Jpsi_list.distr_list.push_back( mc_phys_Jpsi_extr/a_distr);
    mc_extr_etac_list.distr_list.push_back( mc_phys_etac_extr/a_distr);


    //Extrapolate f_etac using mc_phys_extr
    vector<distr_t> fetac_list;
    if(V_charm_1_L.Tag[i_ens].substr(1,1) != "D") fetac_list = {fetac_L, fetac_M, fetac_H};
    else fetac_list = {fetac_L, fetac_M};
    distr_t fetac_extr= Obs_extrapolation_meson_mass( fetac_list, mc_list, mc_phys_extr, "../data/gm2/charm", "f_etac_extr_quark_mass_"+tag_mc_extrapolation+"_"+V_charm_1_L.Tag[i_ens], UseJack, "SPLINE");


    //Extrapolate Zv and Za using mc_phys_extr
    distr_t Za_hadr_extr = Obs_extrapolation_meson_mass( Za_hadr_list, mc_list, mc_phys_extr, "../data/gm2/charm", "Za_extr_quark_mass_"+tag_mc_extrapolation+"_"+V_charm_1_L.Tag[i_ens], UseJack, "SPLINE");
    distr_t Zv_hadr_extr = Obs_extrapolation_meson_mass( Zv_hadr_list, mc_list, mc_phys_extr, "../data/gm2/charm", "Zv_extr_quark_mass_"+tag_mc_extrapolation+"_"+V_charm_1_L.Tag[i_ens], UseJack, "SPLINE");

    //push_back the result
    Zv_fit_charm_Extr.distr_list.push_back(Zv_hadr_extr);
    Za_fit_charm_Extr.distr_list.push_back(Za_hadr_extr);
    Zv_diff_charm.distr_list.push_back( Zv_hadr_extr- Zv_WI_distr);
    Za_diff_charm.distr_list.push_back( Za_hadr_extr-Za_WI_distr);
    Zv_diff_strange_charm.distr_list.push_back( Zv_hadr_extr-Zv_WI_strange_distr);
    Za_diff_strange_charm.distr_list.push_back( Za_hadr_extr-Za_WI_strange_distr);


    if(Use_Za_Zv_from_charm_run && Use_Extrapolated_Za_Zv_charm) {
      Za_L = Za_hadr_extr; Za_M = Za_hadr_extr; Za_H = Za_hadr_extr; Zv_L= Zv_hadr_extr; Zv_M = Zv_hadr_extr; Zv_H = Zv_hadr_extr;
      if(V_charm_1_L.Tag[i_ens].substr(1,1) == "A") {
	Za_L = Za_WI_charm_extr_distr ; Za_M = Za_WI_charm_extr_distr; Za_H = Za_WI_charm_extr_distr ;
	Zv_L = Zv_WI_charm_extr_distr ; Zv_M = Zv_WI_charm_extr_distr; Zv_H = Zv_WI_charm_extr_distr ;
      }
    }
  




    if(V_charm_1_L.Tag[i_ens].substr(1,1) == "A") { //push_back Zv and Za light from ensembles of type A}
      Za_fit_silvano_light_A_ens.distr_list.push_back(Za_hadr_light);
      Zv_fit_silvano_light_A_ens.distr_list.push_back(Zv_hadr_light);
      Za_fit_silvano_charm_A_ens.distr_list.push_back(Za_hadr_extr);
      Zv_fit_silvano_charm_A_ens.distr_list.push_back(Zv_hadr_extr);
      Mpi_fit_silvano_A_ens.distr_list.push_back(Mpi*Mpi/(a_distr*a_distr));
      A_ens_charm_silvano_tags.push_back( V_charm_1_L.Tag[i_ens]);
    }




  
    //push_back MV_charm, ZV_charm, Mpi
 
    Mpi_fit_charm.distr_list.push_back(Mpi);
    fp_fit_charm.distr_list.push_back(fpi);
    //L
    MV_fit_charm_L.distr_list.push_back(MV_charm_L);
    ZV_fit_charm_L.distr_list.push_back(Za_L*Za_L*ZV_charm_L);
    //M
    MV_fit_charm_M.distr_list.push_back(MV_charm_M);
    ZV_fit_charm_M.distr_list.push_back(Za_M*Za_M*ZV_charm_M);
    //H
    MV_fit_charm_H.distr_list.push_back(MV_charm_H);
    ZV_fit_charm_H.distr_list.push_back(Za_H*Za_H*ZV_charm_H);

  


    double Zv2_pert=1.0;
    double Za2_pert=1.0;
  
    

    //free corr LO artifacts
    //################## READ FREE THEORY VECTOR-VECTOR CORRELATOR OPPOSITE R ####################################
    string Pt_free_oppor_L= "../Vkvk_cont/"+to_string(Corr.Nt/2)+"_m"+to_string_with_precision(L_info.mc_L,3)+"/OPPOR";
    string Pt_free_oppor_M= "../Vkvk_cont/"+to_string(Corr.Nt/2)+"_m"+to_string_with_precision(L_info.mc_M,3)+"/OPPOR";
    string Pt_free_oppor_H= "../Vkvk_cont/"+to_string(Corr.Nt/2)+"_m"+to_string_with_precision(L_info.mc_H,3)+"/OPPOR";
    Vfloat VV_free_oppor_L= Read_From_File(Pt_free_oppor_L, 1, 4);
    Vfloat VV_free_oppor_M= Read_From_File(Pt_free_oppor_M, 1, 4);
    Vfloat VV_free_oppor_H= Read_From_File(Pt_free_oppor_H, 1, 4);
    if(VV_free_oppor_L.size() != Corr.Nt) crash("Failed to read properly free VV correlator mc_L w opposite r");
    if(VV_free_oppor_M.size() != Corr.Nt) crash("Failed to read properly free VV correlator mc_M w opposite r");
    if(VV_free_oppor_H.size() != Corr.Nt) crash("Failed to read properly free VV correlator mc_H w opposite r");
    //################## READ FREE THEORY VECTOR-VECTOR CORRELATOR SAME R ####################################
    string Pt_free_samer_L= "../Vkvk_cont/"+to_string(Corr.Nt/2)+"_m"+to_string_with_precision(L_info.mc_L,3)+"/SAMER";
    string Pt_free_samer_M= "../Vkvk_cont/"+to_string(Corr.Nt/2)+"_m"+to_string_with_precision(L_info.mc_M,3)+"/SAMER";
    string Pt_free_samer_H= "../Vkvk_cont/"+to_string(Corr.Nt/2)+"_m"+to_string_with_precision(L_info.mc_H,3)+"/SAMER";
    Vfloat VV_free_samer_L= Read_From_File(Pt_free_samer_L, 1, 4);
    Vfloat VV_free_samer_M= Read_From_File(Pt_free_samer_M, 1, 4);
    Vfloat VV_free_samer_H= Read_From_File(Pt_free_samer_H, 1, 4);
    if(VV_free_samer_L.size() != Corr.Nt) crash("Failed to read properly free VV correlator mc_L  w same r");
    if(VV_free_samer_M.size() != Corr.Nt) crash("Failed to read properly free VV correlator mc_M  w same r");
    if(VV_free_samer_H.size() != Corr.Nt) crash("Failed to read properly free VV correlator mc_H  w same r");
    //Insert electric charges
    for( auto & OP:VV_free_oppor_L) OP *= Za2_pert*pert_corr_charm_on_off*qc*qc;
    for( auto & SA:VV_free_samer_L) SA *= Zv2_pert*pert_corr_charm_on_off*qc*qc;
    for( auto & OP:VV_free_oppor_M) OP *= Za2_pert*pert_corr_charm_on_off*qc*qc;
    for( auto & SA:VV_free_samer_M) SA *= Zv2_pert*pert_corr_charm_on_off*qc*qc;
    for( auto & OP:VV_free_oppor_H) OP *= Za2_pert*pert_corr_charm_on_off*qc*qc;
    for( auto & SA:VV_free_samer_H) SA *= Zv2_pert*pert_corr_charm_on_off*qc*qc;
  

  
    Vfloat free_corr_log_art(Corr.Nt, 0.0);
    for(int t=0;t<Corr.Nt;t++) {  if( t*a_distr.ave() < 1.0*fm_to_inv_Gev && t != 0) {
	free_corr_log_art[t] = -1.0*pert_corr_charm_on_off*(qc*qc)*(1.0/(2.0*M_PI*M_PI*pow(t,5)));
      }
      if(t==0 || t*a_distr.ave() > add_pert_corr_charm_up_to*fm_to_inv_Gev) { VV_free_samer_L[t] =0; VV_free_samer_M[t] = 0; VV_free_samer_H[t] =0; VV_free_oppor_L[t] = 0; VV_free_oppor_M[t]=0; VV_free_oppor_H[t]=0;}
    }

    distr_t_list V_charm_distr_L_pert_sub, V_charm_OS_distr_L_pert_sub;
    distr_t_list V_charm_distr_M_pert_sub, V_charm_OS_distr_M_pert_sub;
    distr_t_list V_charm_distr_H_pert_sub, V_charm_OS_distr_H_pert_sub;

    if(!sum_pert_corr_charm_to_bare_corr) { //sum to renormalized correlator
      //L
      V_charm_distr_L_pert_sub = (1.0/(Za_L*Za_L))*(Za_L*Za_L*V_charm_distr_L + VV_free_oppor_L);
      V_charm_OS_distr_L_pert_sub = (1.0/(Zv_L*Zv_L))*(Zv_L*Zv_L*V_charm_OS_distr_L + VV_free_samer_L);
      //M
      V_charm_distr_M_pert_sub = (1.0/(Za_M*Za_M))*(Za_M*Za_M*V_charm_distr_M + VV_free_oppor_M);
      V_charm_OS_distr_M_pert_sub =  (1.0/(Zv_M*Zv_M))*(Zv_M*Zv_M*V_charm_OS_distr_M + VV_free_samer_M);
      //H
      V_charm_distr_H_pert_sub = (1.0/(Za_H*Za_H))*(Za_H*Za_H*V_charm_distr_H + VV_free_oppor_H);
      V_charm_OS_distr_H_pert_sub =  (1.0/(Zv_H*Zv_H))*(Zv_H*Zv_H*V_charm_OS_distr_H + VV_free_samer_H);
    }
    else { //sum to bare correlator
      //L
      V_charm_distr_L_pert_sub = (V_charm_distr_L + VV_free_oppor_L);
      V_charm_OS_distr_L_pert_sub = (V_charm_OS_distr_L + VV_free_samer_L);
      //M
      V_charm_distr_M_pert_sub = (V_charm_distr_M + VV_free_oppor_M);
      V_charm_OS_distr_M_pert_sub =  (V_charm_OS_distr_M + VV_free_samer_M);
      //H
      V_charm_distr_H_pert_sub = (V_charm_distr_H + VV_free_oppor_H);
      V_charm_OS_distr_H_pert_sub =  (V_charm_OS_distr_H + VV_free_samer_H);
    }

 

    // print summed correlators to file
    //tm
    //L
    Print_To_File({}, {V_charm_distr_L.ave(), V_charm_distr_L.err(), (Za_L*Za_L*V_charm_distr_L).ave(), (Za_L*Za_L*V_charm_distr_L).err()}, "../data/gm2/charm/tm_"+Extrapolation_charm_mode+"/corr_sum_"+V_charm_1_L.Tag[i_ens]+"_L.dat.t", "", "#t   bare  renormalized");
    //M
    Print_To_File({}, {V_charm_distr_M.ave(), V_charm_distr_M.err(), (Za_M*Za_M*V_charm_distr_M).ave(), (Za_M*Za_M*V_charm_distr_M).err()}, "../data/gm2/charm/tm_"+Extrapolation_charm_mode+"/corr_sum_"+V_charm_1_M.Tag[i_ens]+"_M.dat.t", "", "#t   bare  renormalized");
    //H
    Print_To_File({}, {V_charm_distr_H.ave(), V_charm_distr_H.err(), (Za_H*Za_H*V_charm_distr_H).ave(), (Za_H*Za_H*V_charm_distr_H).err()}, "../data/gm2/charm/tm_"+Extrapolation_charm_mode+"/corr_sum_"+V_charm_1_H.Tag[i_ens]+"_H.dat.t", "", "#t   bare  renormalized");
    //pert sub
    Print_To_File({}, { (Za_L*Za_L*V_charm_distr_L_pert_sub).ave(), (Za_L*Za_L*V_charm_distr_L_pert_sub).err()}, "../data/gm2/charm/tm_"+Extrapolation_charm_mode+"/pert_sub_corr_"+V_charm_1_L.Tag[i_ens]+"_L.dat.t", "", "#t   bare  renormalized");
    //M
    Print_To_File({}, {(Za_M*Za_M*V_charm_distr_M_pert_sub).ave(), (Za_M*Za_M*V_charm_distr_M_pert_sub).err()}, "../data/gm2/charm/tm_"+Extrapolation_charm_mode+"/pert_sub_corr_"+V_charm_1_M.Tag[i_ens]+"_M.dat.t", "", "#t   bare  renormalized");
    //H
    Print_To_File({}, {(Za_H*Za_H*V_charm_distr_H_pert_sub).ave(), (Za_H*Za_H*V_charm_distr_H_pert_sub).err()}, "../data/gm2/charm/tm_"+Extrapolation_charm_mode+"/pert_sub_corr_"+V_charm_1_H.Tag[i_ens]+"_H.dat.t", "", "#t   bare  renormalized");

  
    //OS
    //L
    Print_To_File({}, {V_charm_OS_distr_L.ave(), V_charm_OS_distr_L.err(), (Zv_L*Zv_L*V_charm_OS_distr_L).ave(), (Zv_H*Zv_H*V_charm_OS_distr_L).err()}, "../data/gm2/charm/OS_"+Extrapolation_charm_mode+"/corr_sum_"+V_charm_1_L.Tag[i_ens]+"_L.dat.t", "", "#t   bare  renormalized");
    //M
    Print_To_File({}, {V_charm_OS_distr_M.ave(), V_charm_OS_distr_M.err(), (Zv_M*Zv_M*V_charm_OS_distr_M).ave(), (Zv_H*Zv_H*V_charm_OS_distr_M).err()}, "../data/gm2/charm/OS_"+Extrapolation_charm_mode+"/corr_sum_"+V_charm_1_M.Tag[i_ens]+"_M.dat.t", "", "#t   bare  renormalized");
    //H
    Print_To_File({}, {V_charm_OS_distr_H.ave(), V_charm_OS_distr_H.err(), (Zv_H*Zv_H*V_charm_OS_distr_H).ave(), (Zv_H*Zv_H*V_charm_OS_distr_H).err()}, "../data/gm2/charm/OS_"+Extrapolation_charm_mode+"/corr_sum_"+V_charm_1_H.Tag[i_ens]+"_H.dat.t", "", "#t   bare  renormalized");
    //pert sub
    Print_To_File({}, {(Zv_L*Zv_L*V_charm_OS_distr_L_pert_sub).ave(), (Zv_H*Zv_H*V_charm_OS_distr_L_pert_sub).err()}, "../data/gm2/charm/OS_"+Extrapolation_charm_mode+"/pert_sub_corr_"+V_charm_1_L.Tag[i_ens]+"_L.dat.t", "", "#t   bare  renormalized");
    //M
    Print_To_File({}, {(Zv_M*Zv_M*V_charm_OS_distr_M_pert_sub).ave(), (Zv_H*Zv_H*V_charm_OS_distr_M_pert_sub).err()}, "../data/gm2/charm/OS_"+Extrapolation_charm_mode+"/pert_sub_corr_"+V_charm_1_M.Tag[i_ens]+"_M.dat.t", "", "#t   bare  renormalized");
    //H
    Print_To_File({}, {(Zv_H*Zv_H*V_charm_OS_distr_H_pert_sub).ave(), (Zv_H*Zv_H*V_charm_OS_distr_H_pert_sub).err()}, "../data/gm2/charm/OS_"+Extrapolation_charm_mode+"/pert_sub_corr_"+V_charm_1_H.Tag[i_ens]+"_H.dat.t", "", "#t   bare  renormalized");


    //print disco charm
    if(Include_charm_disco && Found_disco_ens) Print_To_File({}, {disco_distr.ave(), disco_distr.err(), (Zv_L*Zv_L*disco_distr).ave(), (Zv_L*Zv_L*disco_distr).err()}, "../data/gm2/charm/disco/disc_"+V_charm_1_L.Tag[i_ens]+".dat.t","","# bare renormalized");

    //improved
    if(Include_charm_disco && Found_disco_impr_ens) Print_To_File({}, {disco_impr_distr.ave(), disco_impr_distr.err(), (Zv_L*Zv_L*disco_impr_distr).ave(), (Zv_L*Zv_L*disco_impr_distr).err()}, "../data/gm2/charm/disco/disc_impr_"+V_charm_1_L.Tag[i_ens]+".dat.t","","# bare renormalized");

    //off-diagonal light-charm
    if(Include_off_diagonal_disco && Found_disco_impr_light_charm_ens) Print_To_File({}, {disco_impr_light_charm_distr.ave(), disco_impr_light_charm_distr.err(), (Zv_L*Zv_L*disco_impr_light_charm_distr).ave(), (Zv_L*Zv_L*disco_impr_light_charm_distr).err()}, "../data/gm2/light_charm/disco/disc_impr_"+V_charm_1_L.Tag[i_ens]+".dat.t","","# bare renormalized");
    //off-diagonal lightD-charm
    if(Include_off_diagonal_disco && Found_disco_impr_lightD_charm_ens) Print_To_File({}, {disco_impr_lightD_charm_distr.ave(), disco_impr_lightD_charm_distr.err(), (Zv_L*Zv_L*disco_impr_lightD_charm_distr).ave(), (Zv_L*Zv_L*disco_impr_lightD_charm_distr).err()}, "../data/gm2/light_charm/disco/disc_impr_D_"+V_charm_1_L.Tag[i_ens]+".dat.t","","# bare renormalized");
    


  

  
 
    //OS  
    Corr.Tmin=Tmin_VV;
    Corr.Tmax=Tmax_VV;
    //L
    MV_charm_OS_distr_L= Corr.effective_mass_t(V_charm_OS_distr_L, "../data/gm2/charm/OS_"+Extrapolation_charm_mode+"/MV_mass_"+V_charm_1_L.Tag[i_ens]+"_L.dat");
    ZV_charm_OS_distr_L= Corr.residue_t(V_charm_OS_distr_L, "../data/gm2/charm/OS_"+Extrapolation_charm_mode+"/ZV_overlap_"+V_charm_1_L.Tag[i_ens]+"_L.dat");
    MV_charm_OS_L = Corr.Fit_distr(MV_charm_OS_distr_L);
    ZV_charm_OS_L = Corr.Fit_distr(ZV_charm_OS_distr_L);
    //M
    MV_charm_OS_distr_M= Corr.effective_mass_t(V_charm_OS_distr_M, "../data/gm2/charm/OS_"+Extrapolation_charm_mode+"/MV_mass_"+V_charm_1_M.Tag[i_ens]+"_M.dat");
    ZV_charm_OS_distr_M= Corr.residue_t(V_charm_OS_distr_M, "../data/gm2/charm/OS_"+Extrapolation_charm_mode+"/ZV_overlap_"+V_charm_1_M.Tag[i_ens]+"_M.dat");
    MV_charm_OS_M = Corr.Fit_distr(MV_charm_OS_distr_M);
    ZV_charm_OS_M = Corr.Fit_distr(ZV_charm_OS_distr_M);
    //H
    MV_charm_OS_distr_H= Corr.effective_mass_t(V_charm_OS_distr_H, "../data/gm2/charm/OS_"+Extrapolation_charm_mode+"/MV_mass_"+V_charm_1_H.Tag[i_ens]+"_H.dat");
    ZV_charm_OS_distr_H= Corr.residue_t(V_charm_OS_distr_H, "../data/gm2/charm/OS_"+Extrapolation_charm_mode+"/ZV_overlap_"+V_charm_1_H.Tag[i_ens]+"_H.dat");
    MV_charm_OS_H = Corr.Fit_distr(MV_charm_OS_distr_H);
    ZV_charm_OS_H = Corr.Fit_distr(ZV_charm_OS_distr_H);

  

    //push_back MV_charm, ZV_charm Mpi_OS
    Mpi_OS_fit_charm.distr_list.push_back(Mpi_OS);
    //L
    MV_fit_charm_OS_L.distr_list.push_back(MV_charm_OS_L);
    ZV_fit_charm_OS_L.distr_list.push_back(Zv_L*Zv_L*ZV_charm_OS_L);
    //M
    MV_fit_charm_OS_M.distr_list.push_back(MV_charm_OS_M);
    ZV_fit_charm_OS_M.distr_list.push_back(Zv_M*Zv_M*ZV_charm_OS_M);
    //H
    MV_fit_charm_OS_H.distr_list.push_back(MV_charm_OS_H);
    ZV_fit_charm_OS_H.distr_list.push_back(Zv_H*Zv_H*ZV_charm_OS_H);
 


    int Tdata_min= 8;
    int Tdata_max = Corr.Nt/2.0 -1;
    int Tdata_fit = (Tmin_VV+Tmax_VV)/2;

    distr_t ELM_mass_L, ELM_mass_M, ELM_mass_H;
    distr_t ELM_mass_OS_L, ELM_mass_OS_M, ELM_mass_OS_H;

    if(ELM_mass_charm=="Jpsi") { ELM_mass_L = MV_charm_L/m_Jpsi_phys_distr; ELM_mass_M = MV_charm_M/m_Jpsi_phys_distr; ELM_mass_H = MV_charm_H/m_Jpsi_phys_distr; ELM_mass_OS_L = MV_charm_OS_L/m_Jpsi_phys_distr; ELM_mass_OS_M = MV_charm_OS_M/m_Jpsi_phys_distr; ELM_mass_OS_H = MV_charm_OS_H/m_Jpsi_phys_distr;}
    else if(ELM_mass_charm=="etac") { ELM_mass_L = M_etaC_L/m_etac_phys_distr; ELM_mass_M = M_etaC_M/m_etac_phys_distr; ELM_mass_H = M_etaC_H/m_etac_phys_distr; ELM_mass_OS_L = cbar_c_OS_mass_L/m_etac_phys_distr; ELM_mass_OS_M = cbar_c_OS_mass_M/m_etac_phys_distr; ELM_mass_OS_H = cbar_c_OS_mass_H/m_etac_phys_distr;}
    else crash("ELM mass charm is invalid: "+ELM_mass_charm);
  
    //compute kernel distribution
    //NON-ELM
    distr_t_list Kernel_distr = distr_t_list::f_of_distr(K, a_distr, Upper_Limit_Time_Integral_charm+1);

  
    //tm
    //L
    distr_t_list Kernel_distr_L = distr_t_list::f_of_distr(K,ELM_mass_L, Upper_Limit_Time_Integral_charm+1);
    //M
    distr_t_list Kernel_distr_M = distr_t_list::f_of_distr(K,ELM_mass_M, Upper_Limit_Time_Integral_charm+1);
    //H
    distr_t_list Kernel_distr_H = distr_t_list::f_of_distr(K,ELM_mass_H, Upper_Limit_Time_Integral_charm+1);

  
    //OS
    //L
    distr_t_list Kernel_OS_distr_L = distr_t_list::f_of_distr(K,ELM_mass_OS_L, Upper_Limit_Time_Integral_charm+1);
    //M
    distr_t_list Kernel_OS_distr_M = distr_t_list::f_of_distr(K,ELM_mass_OS_M, Upper_Limit_Time_Integral_charm+1);
    //H
    distr_t_list Kernel_OS_distr_H = distr_t_list::f_of_distr(K,ELM_mass_OS_H, Upper_Limit_Time_Integral_charm+1);


  
    //compute exp(-Mv*t) distribution
    //tm
    //L
    distr_t_list exp_MVc_L = distr_t_list::f_of_distr(exp_MV, MV_charm_L, Upper_Limit_Time_Integral_charm+1);
    //M
    distr_t_list exp_MVc_M = distr_t_list::f_of_distr(exp_MV, MV_charm_M, Upper_Limit_Time_Integral_charm+1);
    //H
    distr_t_list exp_MVc_H = distr_t_list::f_of_distr(exp_MV, MV_charm_H, Upper_Limit_Time_Integral_charm+1);
 
    //OS
    //L
    distr_t_list exp_OS_MVc_L = distr_t_list::f_of_distr(exp_MV, MV_charm_OS_L, Upper_Limit_Time_Integral_charm+1);
    //M
    distr_t_list exp_OS_MVc_M = distr_t_list::f_of_distr(exp_MV, MV_charm_OS_M, Upper_Limit_Time_Integral_charm+1);
    //H
    distr_t_list exp_OS_MVc_H = distr_t_list::f_of_distr(exp_MV, MV_charm_OS_H, Upper_Limit_Time_Integral_charm+1);

  
    //Print single-exponential prediction to file
    //tm
    //L
    Print_To_File({}, {(exp_MVc_L*(ZV_charm_L/(2.0*MV_charm_L))).ave(), (exp_MVc_L*(ZV_charm_L/(2.0*MV_charm_L))).err(), (Za_L*Za_L*exp_MVc_L*(ZV_charm_L/(2.0*MV_charm_L))).ave(), (Za_L*Za_L*exp_MVc_L*(ZV_charm_L/(2.0*MV_charm_L))).err() }, "../data/gm2/charm/tm_"+Extrapolation_charm_mode+"/corr_gsd_sum_"+V_charm_1_L.Tag[i_ens]+"_L.dat.t", "", "#t  bare  renormalized");
    //M
    Print_To_File({}, {(exp_MVc_M*(ZV_charm_M/(2.0*MV_charm_M))).ave(), (exp_MVc_M*(ZV_charm_M/(2.0*MV_charm_M))).err(), (Za_M*Za_M*exp_MVc_M*(ZV_charm_M/(2.0*MV_charm_M))).ave(), (Za_M*Za_M*exp_MVc_M*(ZV_charm_M/(2.0*MV_charm_M))).err() }, "../data/gm2/charm/tm_"+Extrapolation_charm_mode+"/corr_gsd_sum_"+V_charm_1_M.Tag[i_ens]+"_M.dat.t", "", "#t  bare  renormalized");
    //H
    Print_To_File({}, {(exp_MVc_H*(ZV_charm_H/(2.0*MV_charm_H))).ave(), (exp_MVc_H*(ZV_charm_H/(2.0*MV_charm_H))).err(), (Za_H*Za_H*exp_MVc_H*(ZV_charm_H/(2.0*MV_charm_H))).ave(), (Za_H*Za_H*exp_MVc_H*(ZV_charm_H/(2.0*MV_charm_H))).err() }, "../data/gm2/charm/tm_"+Extrapolation_charm_mode+"/corr_gsd_sum_"+V_charm_1_H.Tag[i_ens]+"_H.dat.t", "", "#t  bare  renormalized");

  
    //OS
    //L
    Print_To_File({}, {(exp_OS_MVc_L*(ZV_charm_OS_L/(2.0*MV_charm_OS_L))).ave(), (exp_OS_MVc_L*(ZV_charm_OS_L/(2.0*MV_charm_OS_L))).err(), (Zv_L*Zv_L*exp_OS_MVc_L*(ZV_charm_OS_L/(2.0*MV_charm_OS_L))).ave(), (Zv_L*Zv_L*exp_OS_MVc_L*(ZV_charm_OS_L/(2.0*MV_charm_OS_L))).err() }, "../data/gm2/charm/OS_"+Extrapolation_charm_mode+"/corr_gsd_sum_"+V_charm_1_L.Tag[i_ens]+"_L.dat.t", "", "#t  bare  renormalized");
    //M
    Print_To_File({}, {(exp_OS_MVc_M*(ZV_charm_OS_M/(2.0*MV_charm_OS_M))).ave(), (exp_OS_MVc_M*(ZV_charm_OS_M/(2.0*MV_charm_OS_M))).err(), (Zv_M*Zv_M*exp_OS_MVc_M*(ZV_charm_OS_M/(2.0*MV_charm_OS_M))).ave(), (Zv_M*Zv_M*exp_OS_MVc_M*(ZV_charm_OS_M/(2.0*MV_charm_OS_M))).err() }, "../data/gm2/charm/OS_"+Extrapolation_charm_mode+"/corr_gsd_sum_"+V_charm_1_M.Tag[i_ens]+"_M.dat.t", "", "#t  bare  renormalized");
    //H
    Print_To_File({}, {(exp_OS_MVc_H*(ZV_charm_OS_H/(2.0*MV_charm_OS_H))).ave(), (exp_OS_MVc_H*(ZV_charm_OS_H/(2.0*MV_charm_OS_H))).err(), (Zv_H*Zv_H*exp_OS_MVc_H*(ZV_charm_OS_H/(2.0*MV_charm_OS_H))).ave(), (Zv_H*Zv_H*exp_OS_MVc_H*(ZV_charm_OS_H/(2.0*MV_charm_OS_H))).err() }, "../data/gm2/charm/OS_"+Extrapolation_charm_mode+"/corr_gsd_sum_"+V_charm_1_H.Tag[i_ens]+"_H.dat.t", "", "#t  bare  renormalized");
  

    distr_t_list agm2_distr_Tdata_L(UseJack), agm2_OS_distr_Tdata_L(UseJack);
    distr_t_list agm2_distr_Tdata_M(UseJack), agm2_OS_distr_Tdata_M(UseJack);
    distr_t_list agm2_distr_Tdata_H(UseJack), agm2_OS_distr_Tdata_H(UseJack);
    distr_t_list agm2_distr_Tdata_No_ELM_L(UseJack), agm2_OS_distr_Tdata_No_ELM_L(UseJack);
    distr_t_list agm2_distr_Tdata_No_ELM_M(UseJack), agm2_OS_distr_Tdata_No_ELM_M(UseJack);
    distr_t_list agm2_distr_Tdata_No_ELM_H(UseJack), agm2_OS_distr_Tdata_No_ELM_H(UseJack);
    Vfloat Tdata_vec;
    bool Find_Tdata_fit=false;
  
    for(int Tdata=Tdata_min;Tdata<Tdata_max;Tdata++) {
      //compute 4\pia^2 using lattice data up to Tdata (included)

      //L
      distr_t agm2_L(UseJack, UseJack?Njacks:Nboots);
      distr_t agm2_OS_L(UseJack, UseJack?Njacks:Nboots);
      //M
      distr_t agm2_M(UseJack, UseJack?Njacks:Nboots);
      distr_t agm2_OS_M(UseJack, UseJack?Njacks:Nboots);
      //H
      distr_t agm2_H(UseJack, UseJack?Njacks:Nboots);
      distr_t agm2_OS_H(UseJack, UseJack?Njacks:Nboots);

      //NON-ELM
      //L
      distr_t agm2_No_ELM_L(UseJack, UseJack?Njacks:Nboots);
      distr_t agm2_OS_No_ELM_L(UseJack, UseJack?Njacks:Nboots);
      //M
      distr_t agm2_No_ELM_M(UseJack, UseJack?Njacks:Nboots);
      distr_t agm2_OS_No_ELM_M(UseJack, UseJack?Njacks:Nboots);
      //H
      distr_t agm2_No_ELM_H(UseJack, UseJack?Njacks:Nboots);
      distr_t agm2_OS_No_ELM_H(UseJack, UseJack?Njacks:Nboots);

    
      for(int t=1;t<=Upper_Limit_Time_Integral_charm;t++) {
	if(t<=Tdata) {
	  //L
	  agm2_L = agm2_L + w(t,Simps_ord)*4.0*pow(alpha,2)*V_charm_distr_L_pert_sub.distr_list[t]*Kernel_distr_L.distr_list[t];
	  agm2_OS_L = agm2_OS_L + w(t,Simps_ord)*4.0*pow(alpha,2)*V_charm_OS_distr_L_pert_sub.distr_list[t]*Kernel_OS_distr_L.distr_list[t];
	  //M
	  agm2_M = agm2_M +w(t,Simps_ord)*4.0*pow(alpha,2)*V_charm_distr_M_pert_sub.distr_list[t]*Kernel_distr_M.distr_list[t];
	  agm2_OS_M = agm2_OS_M + w(t,Simps_ord)*4.0*pow(alpha,2)*V_charm_OS_distr_M_pert_sub.distr_list[t]*Kernel_OS_distr_M.distr_list[t];
	  //H
	  agm2_H = agm2_H + w(t,Simps_ord)*4.0*pow(alpha,2)*V_charm_distr_H_pert_sub.distr_list[t]*Kernel_distr_H.distr_list[t];
	  agm2_OS_H = agm2_OS_H + w(t,Simps_ord)*4.0*pow(alpha,2)*V_charm_OS_distr_H_pert_sub.distr_list[t]*Kernel_OS_distr_H.distr_list[t];

	  //NON_ELM
	  //L
	  agm2_No_ELM_L = agm2_No_ELM_L + w(t,Simps_ord)*4.0*pow(alpha,2)*V_charm_distr_L_pert_sub.distr_list[t]*Kernel_distr.distr_list[t];
	  agm2_OS_No_ELM_L = agm2_OS_No_ELM_L + w(t,Simps_ord)*4.0*pow(alpha,2)*V_charm_OS_distr_L_pert_sub.distr_list[t]*Kernel_distr.distr_list[t];
	  //M
	  agm2_No_ELM_M = agm2_No_ELM_M + w(t,Simps_ord)*4.0*pow(alpha,2)*V_charm_distr_M_pert_sub.distr_list[t]*Kernel_distr.distr_list[t];
	  agm2_OS_No_ELM_M = agm2_OS_No_ELM_M + w(t,Simps_ord)*4.0*pow(alpha,2)*V_charm_OS_distr_M_pert_sub.distr_list[t]*Kernel_distr.distr_list[t];
	  //H
	  agm2_No_ELM_H = agm2_No_ELM_H + w(t,Simps_ord)*4.0*pow(alpha,2)*V_charm_distr_H_pert_sub.distr_list[t]*Kernel_distr.distr_list[t];
	  agm2_OS_No_ELM_H = agm2_OS_No_ELM_H + w(t,Simps_ord)*4.0*pow(alpha,2)*V_charm_OS_distr_H_pert_sub.distr_list[t]*Kernel_distr.distr_list[t];
	
	}
	else {
	  //L
	  agm2_L= agm2_L + w(t,Simps_ord)*4.0*pow(alpha,2)*(ZV_charm_L/(2.0*MV_charm_L))*exp_MVc_L.distr_list[t]*Kernel_distr_L.distr_list[t];
	  agm2_OS_L= agm2_OS_L + w(t,Simps_ord)*4.0*pow(alpha,2)*(ZV_charm_OS_L/(2.0*MV_charm_OS_L))*exp_OS_MVc_L.distr_list[t]*Kernel_OS_distr_L.distr_list[t];
	  //M
	  agm2_M= agm2_M + w(t,Simps_ord)*4.0*pow(alpha,2)*(ZV_charm_M/(2.0*MV_charm_M))*exp_MVc_M.distr_list[t]*Kernel_distr_M.distr_list[t];
	  agm2_OS_M= agm2_OS_M + w(t,Simps_ord)*4.0*pow(alpha,2)*(ZV_charm_OS_M/(2.0*MV_charm_OS_M))*exp_OS_MVc_M.distr_list[t]*Kernel_OS_distr_M.distr_list[t];
	  //H
	  agm2_H= agm2_H + w(t,Simps_ord)*4.0*pow(alpha,2)*(ZV_charm_H/(2.0*MV_charm_H))*exp_MVc_H.distr_list[t]*Kernel_distr_H.distr_list[t];
	  agm2_OS_H= agm2_OS_H + w(t,Simps_ord)*4.0*pow(alpha,2)*(ZV_charm_OS_H/(2.0*MV_charm_OS_H))*exp_OS_MVc_H.distr_list[t]*Kernel_OS_distr_H.distr_list[t];

	  //NON_ELM
	  //L
	  agm2_No_ELM_L= agm2_No_ELM_L + w(t,Simps_ord)*4.0*pow(alpha,2)*(ZV_charm_L/(2.0*MV_charm_L))*exp_MVc_L.distr_list[t]*Kernel_distr.distr_list[t];
	  agm2_OS_No_ELM_L= agm2_OS_No_ELM_L + w(t,Simps_ord)*4.0*pow(alpha,2)*(ZV_charm_OS_L/(2.0*MV_charm_OS_L))*exp_OS_MVc_L.distr_list[t]*Kernel_distr.distr_list[t];
	  //M
	  agm2_No_ELM_M= agm2_No_ELM_M + w(t,Simps_ord)*4.0*pow(alpha,2)*(ZV_charm_M/(2.0*MV_charm_M))*exp_MVc_M.distr_list[t]*Kernel_distr.distr_list[t];
	  agm2_OS_No_ELM_M= agm2_OS_No_ELM_M + w(t,Simps_ord)*4.0*pow(alpha,2)*(ZV_charm_OS_M/(2.0*MV_charm_OS_M))*exp_OS_MVc_M.distr_list[t]*Kernel_distr.distr_list[t];
	  //H
	  agm2_No_ELM_H= agm2_No_ELM_H + w(t,Simps_ord)*4.0*pow(alpha,2)*(ZV_charm_H/(2.0*MV_charm_H))*exp_MVc_H.distr_list[t]*Kernel_distr.distr_list[t];
	  agm2_OS_No_ELM_H= agm2_OS_No_ELM_H + w(t,Simps_ord)*4.0*pow(alpha,2)*(ZV_charm_OS_H/(2.0*MV_charm_OS_H))*exp_OS_MVc_H.distr_list[t]*Kernel_distr.distr_list[t];
	
	}
      }
    
      Tdata_vec.push_back((double)Tdata);
      //L
      agm2_distr_Tdata_L.distr_list.push_back(Za_L*Za_L*agm2_L);
      agm2_OS_distr_Tdata_L.distr_list.push_back(Zv_L*Zv_L*agm2_OS_L);
      //M
      agm2_distr_Tdata_M.distr_list.push_back(Za_M*Za_M*agm2_M);
      agm2_OS_distr_Tdata_M.distr_list.push_back(Zv_M*Zv_M*agm2_OS_M);
      //H
      agm2_distr_Tdata_H.distr_list.push_back(Za_H*Za_H*agm2_H);
      agm2_OS_distr_Tdata_H.distr_list.push_back(Zv_H*Zv_H*agm2_OS_H);

      //NON_ELM
      //L
      agm2_distr_Tdata_No_ELM_L.distr_list.push_back(Za_L*Za_L*agm2_No_ELM_L);
      agm2_OS_distr_Tdata_No_ELM_L.distr_list.push_back(Zv_L*Zv_L*agm2_OS_No_ELM_L);
      //M
      agm2_distr_Tdata_No_ELM_M.distr_list.push_back(Za_M*Za_M*agm2_No_ELM_M);
      agm2_OS_distr_Tdata_No_ELM_M.distr_list.push_back(Zv_M*Zv_M*agm2_OS_No_ELM_M);
      //H
      agm2_distr_Tdata_No_ELM_H.distr_list.push_back(Za_H*Za_H*agm2_No_ELM_H);
      agm2_OS_distr_Tdata_No_ELM_H.distr_list.push_back(Zv_H*Zv_H*agm2_OS_No_ELM_H);

    
      if(Tdata==Tdata_fit) {
	//L
	agm2_charm_L.distr_list.push_back(Za_L*Za_L*agm2_L);
	agm2_charm_OS_L.distr_list.push_back(Zv_L*Zv_L*agm2_OS_L);
	//M
	agm2_charm_M.distr_list.push_back(Za_M*Za_M*agm2_M);
	agm2_charm_OS_M.distr_list.push_back(Zv_M*Zv_M*agm2_OS_M);
	//H
	agm2_charm_H.distr_list.push_back(Za_H*Za_H*agm2_H);
	agm2_charm_OS_H.distr_list.push_back(Zv_H*Zv_H*agm2_OS_H);
      
	//NON_ELM
	//L
	agm2_charm_No_ELM_L.distr_list.push_back(Za_L*Za_L*agm2_No_ELM_L);
	agm2_charm_OS_No_ELM_L.distr_list.push_back(Zv_L*Zv_L*agm2_OS_No_ELM_L);
	//M
	agm2_charm_No_ELM_M.distr_list.push_back(Za_M*Za_M*agm2_No_ELM_M);
	agm2_charm_OS_No_ELM_M.distr_list.push_back(Zv_M*Zv_M*agm2_OS_No_ELM_M);
	//H
	agm2_charm_No_ELM_H.distr_list.push_back(Za_H*Za_H*agm2_No_ELM_H);
	agm2_charm_OS_No_ELM_H.distr_list.push_back(Zv_H*Zv_H*agm2_OS_No_ELM_H);

   
      
	//extrapolate to the physical point
	vector<distr_t> agm2s_charm, agm2s_charm_OS, agm2s_charm_No_ELM, agm2s_charm_OS_No_ELM;
	if(V_charm_1_L.Tag[i_ens].substr(1,1) != "D") {
	  agm2s_charm = {Za_L*Za_L*agm2_L, Za_M*Za_M*agm2_M, Za_H*Za_H*agm2_H};
	  agm2s_charm_OS = { Zv_L*Zv_L*agm2_OS_L, Zv_M*Zv_M*agm2_OS_M, Zv_H*Zv_H*agm2_OS_H};
	  agm2s_charm_No_ELM = {Za_L*Za_L*agm2_No_ELM_L, Za_M*Za_M*agm2_No_ELM_M, Za_H*Za_H*agm2_No_ELM_H};
	  agm2s_charm_OS_No_ELM = { Zv_L*Zv_L*agm2_OS_No_ELM_L, Zv_M*Zv_M*agm2_OS_No_ELM_M, Zv_H*Zv_H*agm2_OS_No_ELM_H};
	}
	else {
	  agm2s_charm = {Za_L*Za_L*agm2_L, Za_M*Za_M*agm2_M};
	  agm2s_charm_OS = { Zv_L*Zv_L*agm2_OS_L, Zv_M*Zv_M*agm2_OS_M};
	  agm2s_charm_No_ELM = {Za_L*Za_L*agm2_No_ELM_L, Za_M*Za_M*agm2_No_ELM_M};
	  agm2s_charm_OS_No_ELM = { Zv_L*Zv_L*agm2_OS_No_ELM_L, Zv_M*Zv_M*agm2_OS_No_ELM_M};
	}
	agm2_charm_Extr.distr_list.push_back( Obs_extrapolation_meson_mass(agm2s_charm, X_2_fit, X_2_phys, "../data/gm2/charm", "agm2_ELM_Extrapolation_"+Extrapolation_charm_mode+"_"+V_charm_1_L.Tag[i_ens], UseJack, "SPLINE"));
	agm2_charm_OS_Extr.distr_list.push_back( Obs_extrapolation_meson_mass(agm2s_charm_OS, X_2_fit, X_2_phys, "../data/gm2/charm", "agm2_OS_ELM_Extrapolation_"+Extrapolation_charm_mode+"_"+V_charm_1_L.Tag[i_ens], UseJack, "SPLINE"));
	agm2_charm_No_ELM_Extr.distr_list.push_back( Obs_extrapolation_meson_mass(agm2s_charm_No_ELM, X_2_fit, X_2_phys, "../data/gm2/charm", "agm2_Extrapolation_"+Extrapolation_charm_mode+"_"+V_charm_1_L.Tag[i_ens], UseJack, "SPLINE"));
	agm2_charm_OS_No_ELM_Extr.distr_list.push_back( Obs_extrapolation_meson_mass(agm2s_charm_OS_No_ELM, X_2_fit, X_2_phys, "../data/gm2/charm", "agm2_OS_Extrapolation_"+Extrapolation_charm_mode+"_"+V_charm_1_L.Tag[i_ens], UseJack, "SPLINE"));
	Find_Tdata_fit=true;
      }
    }

    if(!Find_Tdata_fit) crash("Cannot find the Tdata value: "+to_string(Tdata_fit));



  
    //#######################  INTERMEDIATE AND SHORT-DISTANCE WINDOW ###################################

    
    //############################   TWISTED MASS ######################################################

    //L
    distr_t agm2_W_L(UseJack, UseJack?Njacks:Nboots);  //constructor sets agm2_W to zero by default
    distr_t agm2_SD_L(UseJack, UseJack?Njacks:Nboots);  //constructor sets agm2_SD to zero by default
    distr_t agm2_W_ELM_L(UseJack, UseJack?Njacks:Nboots); //constructor sets agm2_W_ELM to zero by default
    distr_t agm2_SD_ELM_L(UseJack, UseJack?Njacks:Nboots);  //constructor sets agm2_SD_ELM to zero by default
    //M
    distr_t agm2_W_M(UseJack, UseJack?Njacks:Nboots);  //constructor sets agm2_W to zero by default
    distr_t agm2_SD_M(UseJack, UseJack?Njacks:Nboots);  //constructor sets agm2_SD to zero by default
    distr_t agm2_W_ELM_M(UseJack, UseJack?Njacks:Nboots); //constructor sets agm2_W_ELM to zero by default
    distr_t agm2_SD_ELM_M(UseJack, UseJack?Njacks:Nboots);  //constructor sets agm2_SD_ELM to zero by default
    //H
    distr_t agm2_W_H(UseJack, UseJack?Njacks:Nboots);  //constructor sets agm2_W to zero by default
    distr_t agm2_SD_H(UseJack, UseJack?Njacks:Nboots);  //constructor sets agm2_SD to zero by default
    distr_t agm2_W_ELM_H(UseJack, UseJack?Njacks:Nboots); //constructor sets agm2_W_ELM to zero by default
    distr_t agm2_SD_ELM_H(UseJack, UseJack?Njacks:Nboots);  //constructor sets agm2_SD_ELM to zero by default


    //#################################################################################################


  
    distr_t_list Ker = distr_t_list::f_of_distr(K, a_distr , Corr.Nt/2);
    //L
    distr_t_list Ker_ELM_tm_L = distr_t_list::f_of_distr(K, ELM_mass_L, Corr.Nt/2);
    distr_t_list Ker_ELM_OS_L = distr_t_list::f_of_distr(K, ELM_mass_OS_L, Corr.Nt/2);
    //M
    distr_t_list Ker_ELM_tm_M = distr_t_list::f_of_distr(K, ELM_mass_M, Corr.Nt/2);
    distr_t_list Ker_ELM_OS_M = distr_t_list::f_of_distr(K, ELM_mass_OS_M, Corr.Nt/2);
    //H
    distr_t_list Ker_ELM_tm_H = distr_t_list::f_of_distr(K, ELM_mass_H, Corr.Nt/2);
    distr_t_list Ker_ELM_OS_H = distr_t_list::f_of_distr(K, ELM_mass_OS_H, Corr.Nt/2);
    
    //define lambdas for the theta func
    auto th0 = [](double ta) ->double { return 1.0/(1.0 + exp(-2.0*(ta-t0)/Delta));};
    auto th1 = [](double ta) ->double { return 1.0/(1.0 + exp(-2.0*(ta-t1)/Delta));};

  
    for(int t=1; t< Corr.Nt/2; t++) {
      //L
      agm2_W_L = agm2_W_L + w(t,Simps_ord)*4.0*pow(alpha,2)*Za_L*Za_L*V_charm_distr_L.distr_list[t]*Ker.distr_list[t]*( distr_t::f_of_distr(th0, t*a_distr) - distr_t::f_of_distr(th1, t*a_distr));
      agm2_SD_L = agm2_SD_L + w(t,Simps_ord)*4.0*pow(alpha,2)*Za_L*Za_L*(V_charm_distr_L_pert_sub.distr_list[t])*Ker.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*a_distr));
      agm2_W_ELM_L = agm2_W_ELM_L + w(t,Simps_ord)*4.0*pow(alpha,2)*Za_L*Za_L*V_charm_distr_L.distr_list[t]*Ker_ELM_tm_L.distr_list[t]*( distr_t::f_of_distr(th0, t*ELM_mass_L) - distr_t::f_of_distr(th1, t*ELM_mass_L));
      agm2_SD_ELM_L = agm2_SD_ELM_L + w(t,Simps_ord)*4.0*pow(alpha,2)*(Za_L*Za_L*V_charm_distr_L_pert_sub.distr_list[t])*Ker_ELM_tm_L.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*ELM_mass_L));
      //M
      agm2_W_M = agm2_W_M + w(t,Simps_ord)*4.0*pow(alpha,2)*Za_M*Za_M*V_charm_distr_M.distr_list[t]*Ker.distr_list[t]*( distr_t::f_of_distr(th0, t*a_distr) - distr_t::f_of_distr(th1, t*a_distr));
      agm2_SD_M = agm2_SD_M + w(t,Simps_ord)*4.0*pow(alpha,2)*Za_M*Za_M*(V_charm_distr_M_pert_sub.distr_list[t])*Ker.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*a_distr));
      agm2_W_ELM_M = agm2_W_ELM_M + w(t,Simps_ord)*4.0*pow(alpha,2)*Za_M*Za_M*V_charm_distr_M.distr_list[t]*Ker_ELM_tm_M.distr_list[t]*( distr_t::f_of_distr(th0, t*ELM_mass_M) - distr_t::f_of_distr(th1, t*ELM_mass_M));
      agm2_SD_ELM_M = agm2_SD_ELM_M + w(t,Simps_ord)*4.0*pow(alpha,2)*(Za_M*Za_M*V_charm_distr_M_pert_sub.distr_list[t])*Ker_ELM_tm_M.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*ELM_mass_M));
      //H
      agm2_W_H = agm2_W_H + w(t,Simps_ord)*4.0*pow(alpha,2)*Za_H*Za_H*V_charm_distr_H.distr_list[t]*Ker.distr_list[t]*( distr_t::f_of_distr(th0, t*a_distr) - distr_t::f_of_distr(th1, t*a_distr));
      agm2_SD_H = agm2_SD_H + w(t,Simps_ord)*4.0*pow(alpha,2)*Za_H*Za_H*(V_charm_distr_H_pert_sub.distr_list[t])*Ker.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*a_distr));
      agm2_W_ELM_H = agm2_W_ELM_H + w(t,Simps_ord)*4.0*pow(alpha,2)*Za_H*Za_H*V_charm_distr_H.distr_list[t]*Ker_ELM_tm_H.distr_list[t]*( distr_t::f_of_distr(th0, t*ELM_mass_H) - distr_t::f_of_distr(th1, t*ELM_mass_H));
      agm2_SD_ELM_H = agm2_SD_ELM_H + w(t,Simps_ord)*4.0*pow(alpha,2)*(Za_H*Za_H*V_charm_distr_H_pert_sub.distr_list[t])*Ker_ELM_tm_H.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*ELM_mass_H));
    }
  
    //push_back the result

    //L
    agm2_charm_W_L.distr_list.push_back(agm2_W_L);
    agm2_charm_SD_L.distr_list.push_back(agm2_SD_L);
    agm2_charm_W_ELM_L.distr_list.push_back(agm2_W_ELM_L);
    agm2_charm_SD_ELM_L.distr_list.push_back(agm2_SD_ELM_L);
    //M
    agm2_charm_W_M.distr_list.push_back(agm2_W_M);
    agm2_charm_SD_M.distr_list.push_back(agm2_SD_M);
    agm2_charm_W_ELM_M.distr_list.push_back(agm2_W_ELM_M);
    agm2_charm_SD_ELM_M.distr_list.push_back(agm2_SD_ELM_M);
    //H
    agm2_charm_W_H.distr_list.push_back(agm2_W_H);
    agm2_charm_SD_H.distr_list.push_back(agm2_SD_H);
    agm2_charm_W_ELM_H.distr_list.push_back(agm2_W_ELM_H);
    agm2_charm_SD_ELM_H.distr_list.push_back(agm2_SD_ELM_H);

  
    //extrapolate the result to the physical point

    vector<distr_t> agm2s_charm_W, agm2s_charm_SD, agm2s_charm_W_ELM, agm2s_charm_SD_ELM;
    if(V_charm_1_L.Tag[i_ens].substr(1,1) != "D") {
      agm2s_charm_W = {agm2_W_L, agm2_W_M, agm2_W_H};
      agm2s_charm_SD = {agm2_SD_L, agm2_SD_M, agm2_SD_H};
      agm2s_charm_W_ELM = {agm2_W_ELM_L, agm2_W_ELM_M, agm2_W_ELM_H};
      agm2s_charm_SD_ELM = {agm2_SD_ELM_L, agm2_SD_ELM_M, agm2_SD_ELM_H};
    }
    else {
      agm2s_charm_W = {agm2_W_L, agm2_W_M};
      agm2s_charm_SD = {agm2_SD_L, agm2_SD_M};
      agm2s_charm_W_ELM = {agm2_W_ELM_L, agm2_W_ELM_M};
      agm2s_charm_SD_ELM = {agm2_SD_ELM_L, agm2_SD_ELM_M};
    }

  
    agm2_charm_W_Extr.distr_list.push_back( Obs_extrapolation_meson_mass(agm2s_charm_W, X_2_fit, X_2_phys, "../data/gm2/charm", "agm2_W_Extrapolation_"+Extrapolation_charm_mode+"_"+V_charm_1_L.Tag[i_ens], UseJack, "SPLINE"));
    agm2_charm_SD_Extr.distr_list.push_back( Obs_extrapolation_meson_mass(agm2s_charm_SD, X_2_fit, X_2_phys, "../data/gm2/charm", "agm2_SD_Extrapolation_"+Extrapolation_charm_mode+"_"+V_charm_1_L.Tag[i_ens], UseJack, "SPLINE"));
    agm2_charm_W_ELM_Extr.distr_list.push_back( Obs_extrapolation_meson_mass(agm2s_charm_W_ELM, X_2_fit, X_2_phys, "../data/gm2/charm", "agm2_W_ELM_Extrapolation_"+Extrapolation_charm_mode+"_"+V_charm_1_L.Tag[i_ens], UseJack, "SPLINE"));
    agm2_charm_SD_ELM_Extr.distr_list.push_back( Obs_extrapolation_meson_mass(agm2s_charm_SD_ELM, X_2_fit, X_2_phys, "../data/gm2/charm", "agm2_SD_ELM_Extrapolation_"+Extrapolation_charm_mode+"_"+V_charm_1_L.Tag[i_ens], UseJack, "SPLINE"));
    //####################################################################################################



    // ######################################## OSTERWALDER-SEILER #######################################

    //L
    distr_t agm2_W_OS_L(UseJack, UseJack?Njacks:Nboots);  //constructor sets agm2_W_OS to zero by default
    distr_t agm2_SD_OS_L(UseJack, UseJack?Njacks:Nboots);  //constructor sets agm2_SD_OS to zero by default
    distr_t agm2_W_ELM_OS_L(UseJack, UseJack?Njacks:Nboots); //constructor sets agm2_W_ELM_OS to zero by default
    distr_t agm2_SD_ELM_OS_L(UseJack, UseJack?Njacks:Nboots);  //constructor sets agm2_SD_ELM_OS to zero by default
    //M
    distr_t agm2_W_OS_M(UseJack, UseJack?Njacks:Nboots);  //constructor sets agm2_W_OS to zero by default
    distr_t agm2_SD_OS_M(UseJack, UseJack?Njacks:Nboots);  //constructor sets agm2_SD_OS to zero by default
    distr_t agm2_W_ELM_OS_M(UseJack, UseJack?Njacks:Nboots); //constructor sets agm2_W_ELM_OS to zero by default
    distr_t agm2_SD_ELM_OS_M(UseJack, UseJack?Njacks:Nboots);  //constructor sets agm2_SD_ELM_OS to zero by default
    //H
    distr_t agm2_W_OS_H(UseJack, UseJack?Njacks:Nboots);  //constructor sets agm2_W_OS to zero by default
    distr_t agm2_SD_OS_H(UseJack, UseJack?Njacks:Nboots);  //constructor sets agm2_SD_OS to zero by default
    distr_t agm2_W_ELM_OS_H(UseJack, UseJack?Njacks:Nboots); //constructor sets agm2_W_ELM_OS to zero by default
    distr_t agm2_SD_ELM_OS_H(UseJack, UseJack?Njacks:Nboots);  //constructor sets agm2_SD_ELM_OS to zero by default


    //#################################################################################################

    for(int t=1; t< Corr.Nt/2; t++) {
      //L
      agm2_W_OS_L = agm2_W_OS_L + w(t,Simps_ord)*4.0*pow(alpha,2)*Zv_L*Zv_L*V_charm_OS_distr_L.distr_list[t]*Ker.distr_list[t]*( distr_t::f_of_distr(th0, t*a_distr) - distr_t::f_of_distr(th1, t*a_distr));
      agm2_SD_OS_L = agm2_SD_OS_L + w(t,Simps_ord)*4.0*pow(alpha,2)*Zv_L*Zv_L*V_charm_OS_distr_L_pert_sub.distr_list[t]*Ker.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*a_distr));
      agm2_W_ELM_OS_L = agm2_W_ELM_OS_L + w(t,Simps_ord)*4.0*pow(alpha,2)*Zv_L*Zv_L*V_charm_OS_distr_L.distr_list[t]*Ker_ELM_OS_L.distr_list[t]*( distr_t::f_of_distr(th0, t*ELM_mass_OS_L) - distr_t::f_of_distr(th1, t*ELM_mass_OS_L));
      agm2_SD_ELM_OS_L = agm2_SD_ELM_OS_L + w(t,Simps_ord)*4.0*pow(alpha,2)*Zv_L*Zv_L*V_charm_OS_distr_L_pert_sub.distr_list[t]*Ker_ELM_OS_L.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*ELM_mass_OS_L));
      //M
      agm2_W_OS_M = agm2_W_OS_M + w(t,Simps_ord)*4.0*pow(alpha,2)*Zv_M*Zv_M*V_charm_OS_distr_M.distr_list[t]*Ker.distr_list[t]*( distr_t::f_of_distr(th0, t*a_distr) - distr_t::f_of_distr(th1, t*a_distr));
      agm2_SD_OS_M = agm2_SD_OS_M + w(t,Simps_ord)*4.0*pow(alpha,2)*Zv_M*Zv_M*V_charm_OS_distr_M_pert_sub.distr_list[t]*Ker.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*a_distr));
      agm2_W_ELM_OS_M = agm2_W_ELM_OS_M + w(t,Simps_ord)*4.0*pow(alpha,2)*Zv_M*Zv_M*V_charm_OS_distr_M.distr_list[t]*Ker_ELM_OS_M.distr_list[t]*( distr_t::f_of_distr(th0, t*ELM_mass_OS_M) - distr_t::f_of_distr(th1, t*ELM_mass_OS_M));
      agm2_SD_ELM_OS_M = agm2_SD_ELM_OS_M + w(t,Simps_ord)*4.0*pow(alpha,2)*Zv_M*Zv_M*V_charm_OS_distr_M_pert_sub.distr_list[t]*Ker_ELM_OS_M.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*ELM_mass_OS_M));
      //H
      agm2_W_OS_H = agm2_W_OS_H + w(t,Simps_ord)*4.0*pow(alpha,2)*Zv_H*Zv_H*V_charm_OS_distr_H.distr_list[t]*Ker.distr_list[t]*( distr_t::f_of_distr(th0, t*a_distr) - distr_t::f_of_distr(th1, t*a_distr));
      agm2_SD_OS_H = agm2_SD_OS_H + w(t,Simps_ord)*4.0*pow(alpha,2)*Zv_H*Zv_H*V_charm_OS_distr_H_pert_sub.distr_list[t]*Ker.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*a_distr));
      agm2_W_ELM_OS_H = agm2_W_ELM_OS_H + w(t,Simps_ord)*4.0*pow(alpha,2)*Zv_H*Zv_H*V_charm_OS_distr_H.distr_list[t]*Ker_ELM_OS_H.distr_list[t]*( distr_t::f_of_distr(th0, t*ELM_mass_OS_H) - distr_t::f_of_distr(th1, t*ELM_mass_OS_H));
      agm2_SD_ELM_OS_H = agm2_SD_ELM_OS_H + w(t,Simps_ord)*4.0*pow(alpha,2)*Zv_H*Zv_H*V_charm_OS_distr_H_pert_sub.distr_list[t]*Ker_ELM_OS_H.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*ELM_mass_OS_H));
    }
  
    //push_back the result
    //L
    agm2_charm_W_OS_L.distr_list.push_back(agm2_W_OS_L);
    agm2_charm_SD_OS_L.distr_list.push_back(agm2_SD_OS_L);
    agm2_charm_W_ELM_OS_L.distr_list.push_back(agm2_W_ELM_OS_L);
    agm2_charm_SD_ELM_OS_L.distr_list.push_back(agm2_SD_ELM_OS_L);
    //M
    agm2_charm_W_OS_M.distr_list.push_back(agm2_W_OS_M);
    agm2_charm_SD_OS_M.distr_list.push_back(agm2_SD_OS_M);
    agm2_charm_W_ELM_OS_M.distr_list.push_back(agm2_W_ELM_OS_M);
    agm2_charm_SD_ELM_OS_M.distr_list.push_back(agm2_SD_ELM_OS_M);
    //H
    agm2_charm_W_OS_H.distr_list.push_back(agm2_W_OS_H);
    agm2_charm_SD_OS_H.distr_list.push_back(agm2_SD_OS_H);
    agm2_charm_W_ELM_OS_H.distr_list.push_back(agm2_W_ELM_OS_H);
    agm2_charm_SD_ELM_OS_H.distr_list.push_back(agm2_SD_ELM_OS_H);


  
    //extrapolate the result to the physical point
    vector<distr_t> agm2s_charm_W_OS, agm2s_charm_SD_OS, agm2s_charm_W_ELM_OS, agm2s_charm_SD_ELM_OS;

    if(V_charm_1_L.Tag[i_ens].substr(1,1) != "D") {
      agm2s_charm_W_OS = {agm2_W_OS_L, agm2_W_OS_M, agm2_W_OS_H};
      agm2s_charm_SD_OS = {agm2_SD_OS_L, agm2_SD_OS_M, agm2_SD_OS_H};
      agm2s_charm_W_ELM_OS = {agm2_W_ELM_OS_L, agm2_W_ELM_OS_M, agm2_W_ELM_OS_H};
      agm2s_charm_SD_ELM_OS = {agm2_SD_ELM_OS_L, agm2_SD_ELM_OS_M, agm2_SD_ELM_OS_H};
    }
    else {
      agm2s_charm_W_OS = {agm2_W_OS_L, agm2_W_OS_M};
      agm2s_charm_SD_OS = {agm2_SD_OS_L, agm2_SD_OS_M};
      agm2s_charm_W_ELM_OS = {agm2_W_ELM_OS_L, agm2_W_ELM_OS_M};
      agm2s_charm_SD_ELM_OS = {agm2_SD_ELM_OS_L, agm2_SD_ELM_OS_M};
    }

  
    agm2_charm_W_OS_Extr.distr_list.push_back( Obs_extrapolation_meson_mass(agm2s_charm_W_OS, X_2_fit, X_2_phys, "../data/gm2/charm", "agm2_W_Extrapolation_OS_"+Extrapolation_charm_mode+"_"+V_charm_1_L.Tag[i_ens], UseJack, "SPLINE"));
    agm2_charm_SD_OS_Extr.distr_list.push_back( Obs_extrapolation_meson_mass(agm2s_charm_SD_OS, X_2_fit, X_2_phys, "../data/gm2/charm", "agm2_SD_Extrapolation_OS_"+Extrapolation_charm_mode+"_"+V_charm_1_L.Tag[i_ens], UseJack, "SPLINE"));
    agm2_charm_W_ELM_OS_Extr.distr_list.push_back( Obs_extrapolation_meson_mass(agm2s_charm_W_ELM_OS,X_2_fit, X_2_phys, "../data/gm2/charm", "agm2_W_ELM_Extrapolation_OS_"+Extrapolation_charm_mode+"_"+V_charm_1_L.Tag[i_ens], UseJack, "SPLINE"));
    agm2_charm_SD_ELM_OS_Extr.distr_list.push_back( Obs_extrapolation_meson_mass(agm2s_charm_SD_ELM_OS, X_2_fit, X_2_phys, "../data/gm2/charm", "agm2_SD_ELM_Extrapolation_OS_"+Extrapolation_charm_mode+"_"+V_charm_1_L.Tag[i_ens], UseJack, "SPLINE"));

    //####################################################################################################

    //####################################### DISCO CHARM ###############################################
    distr_t agm2_disco_W(UseJack, UseJack?Njacks:Nboots); //constructor sets it to 0
    distr_t agm2_disco_SD(UseJack, UseJack?Njacks:Nboots); //constructor sets it to 0
    distr_t agm2_disco_full(UseJack, UseJack?Njacks:Nboots); //constructor sets it to 0
    //improved
    distr_t agm2_disco_impr_W(UseJack, UseJack?Njacks:Nboots); //constructor sets it to 0
    distr_t agm2_disco_impr_SD(UseJack, UseJack?Njacks:Nboots); //constructor sets it to 0
    distr_t agm2_disco_impr_full(UseJack, UseJack?Njacks:Nboots); //constructor sets it to 0

    //off-diagonal light-charm
    distr_t agm2_disco_impr_lc_W(UseJack, UseJack?Njacks:Nboots); //constructor sets it to 0
    distr_t agm2_disco_impr_lc_SD(UseJack, UseJack?Njacks:Nboots); //constructor sets it to 0
    distr_t agm2_disco_impr_lc_full(UseJack, UseJack?Njacks:Nboots); //constructor sets it to 0

    //off-diagonal lightD-charm
    distr_t agm2_disco_impr_lDc_W(UseJack, UseJack?Njacks:Nboots); //constructor sets it to 0
    distr_t agm2_disco_impr_lDc_SD(UseJack, UseJack?Njacks:Nboots); //constructor sets it to 0
    distr_t agm2_disco_impr_lDc_full(UseJack, UseJack?Njacks:Nboots); //constructor sets it to 0
    
 
    if(Include_charm_disco && Found_disco_ens) {
    
      for(int t=1; t< Corr.Nt/2; t++) {
	agm2_disco_W = agm2_disco_W + w(t,Simps_ord)*4.0*pow(alpha,2)*Zv_L*Zv_L*disco_distr.distr_list[t]*Ker.distr_list[t]*( distr_t::f_of_distr(th0, t*a_distr) - distr_t::f_of_distr(th1, t*a_distr));
	agm2_disco_SD = agm2_disco_SD + w(t,Simps_ord)*4.0*pow(alpha,2)*Zv_L*Zv_L*disco_distr.distr_list[t]*Ker.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*a_distr));
	agm2_disco_full = agm2_disco_full + w(t,Simps_ord)*4.0*pow(alpha,2)*Zv_L*Zv_L*disco_distr.distr_list[t]*Ker.distr_list[t];
      }
  
      //push_back the result

      agm2_disco_charm_W.distr_list.push_back(agm2_disco_W);
      agm2_disco_charm_SD.distr_list.push_back(agm2_disco_SD);
      agm2_disco_charm_No_ELM.distr_list.push_back(agm2_disco_full);
    }

    //improved
    if(Include_charm_disco && Found_disco_impr_ens) {
    
      for(int t=1; t< Corr.Nt/2; t++) {
	agm2_disco_impr_W = agm2_disco_impr_W + w(t,Simps_ord)*4.0*pow(alpha,2)*Zv_L*Zv_L*disco_impr_distr.distr_list[t]*Ker.distr_list[t]*( distr_t::f_of_distr(th0, t*a_distr) - distr_t::f_of_distr(th1, t*a_distr));
	agm2_disco_impr_SD = agm2_disco_impr_SD + w(t,Simps_ord)*4.0*pow(alpha,2)*Zv_L*Zv_L*disco_impr_distr.distr_list[t]*Ker.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*a_distr));
	agm2_disco_impr_full = agm2_disco_impr_full + w(t,Simps_ord)*4.0*pow(alpha,2)*Zv_L*Zv_L*disco_impr_distr.distr_list[t]*Ker.distr_list[t];
      }
  
      //push_back the result

      agm2_disco_impr_charm_W.distr_list.push_back(agm2_disco_impr_W);
      agm2_disco_impr_charm_SD.distr_list.push_back(agm2_disco_impr_SD);
      agm2_disco_impr_charm_No_ELM.distr_list.push_back(agm2_disco_impr_full);
    }


  

    //off-diagonal light-charm improved
    if(Include_charm_disco && Found_disco_impr_light_charm_ens) {
    
      for(int t=1; t< Corr.Nt/2; t++) {
	agm2_disco_impr_lc_W = agm2_disco_impr_lc_W + w(t,Simps_ord)*4.0*pow(alpha,2)*Zv_L*Zv_L*disco_impr_light_charm_distr.distr_list[t]*Ker.distr_list[t]*( distr_t::f_of_distr(th0, t*a_distr) - distr_t::f_of_distr(th1, t*a_distr));
	agm2_disco_impr_lc_SD = agm2_disco_impr_lc_SD + w(t,Simps_ord)*4.0*pow(alpha,2)*Zv_L*Zv_L*disco_impr_light_charm_distr.distr_list[t]*Ker.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*a_distr));
	agm2_disco_impr_lc_full = agm2_disco_impr_lc_full + w(t,Simps_ord)*4.0*pow(alpha,2)*Zv_L*Zv_L*disco_impr_light_charm_distr.distr_list[t]*Ker.distr_list[t];
      }
  
      //push_back the result
      
      agm2_disco_impr_light_charm_W.distr_list.push_back(agm2_disco_impr_lc_W);
      agm2_disco_impr_light_charm_SD.distr_list.push_back(agm2_disco_impr_lc_SD);
      agm2_disco_impr_light_charm_No_ELM.distr_list.push_back(agm2_disco_impr_lc_full);
    }


    //off-diagonal lightD-charm improved
    if(Include_charm_disco && Found_disco_impr_lightD_charm_ens) {
    
      for(int t=1; t< Corr.Nt/2; t++) {
	agm2_disco_impr_lDc_W = agm2_disco_impr_lDc_W + w(t,Simps_ord)*4.0*pow(alpha,2)*Zv_L*Zv_L*disco_impr_lightD_charm_distr.distr_list[t]*Ker.distr_list[t]*( distr_t::f_of_distr(th0, t*a_distr) - distr_t::f_of_distr(th1, t*a_distr));
	agm2_disco_impr_lDc_SD = agm2_disco_impr_lDc_SD + w(t,Simps_ord)*4.0*pow(alpha,2)*Zv_L*Zv_L*disco_impr_lightD_charm_distr.distr_list[t]*Ker.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*a_distr));
	agm2_disco_impr_lDc_full = agm2_disco_impr_lDc_full + w(t,Simps_ord)*4.0*pow(alpha,2)*Zv_L*Zv_L*disco_impr_lightD_charm_distr.distr_list[t]*Ker.distr_list[t];
      }
  
      //push_back the result
      
      agm2_disco_impr_lightD_charm_W.distr_list.push_back(agm2_disco_impr_lDc_W);
      agm2_disco_impr_lightD_charm_SD.distr_list.push_back(agm2_disco_impr_lDc_SD);
      agm2_disco_impr_lightD_charm_No_ELM.distr_list.push_back(agm2_disco_impr_lDc_full);
    }

    
    


  

  
    //print full contribution to file
    //L
    Print_To_File({}, {Tdata_vec, agm2_distr_Tdata_L.ave(), agm2_distr_Tdata_L.err()}, "../data/gm2/charm/tm_"+Extrapolation_charm_mode+"/agm2_Tdata_"+V_charm_1_L.Tag[i_ens]+"_ELM_L.dat.t", "", "#id  Tdata   ag2m agm2_err");
    Print_To_File({}, {Tdata_vec, agm2_OS_distr_Tdata_L.ave(), agm2_OS_distr_Tdata_L.err()}, "../data/gm2/charm/OS_"+Extrapolation_charm_mode+"/agm2_Tdata_"+V_charm_1_L.Tag[i_ens]+"_ELM_L.dat.t", "", "#id  Tdata   ag2m agm2_err");
    //M
    Print_To_File({}, {Tdata_vec, agm2_distr_Tdata_M.ave(), agm2_distr_Tdata_M.err()}, "../data/gm2/charm/tm_"+Extrapolation_charm_mode+"/agm2_Tdata_"+V_charm_1_M.Tag[i_ens]+"_ELM_M.dat.t", "", "#id  Tdata   ag2m agm2_err");
    Print_To_File({}, {Tdata_vec, agm2_OS_distr_Tdata_M.ave(), agm2_OS_distr_Tdata_M.err()}, "../data/gm2/charm/OS_"+Extrapolation_charm_mode+"/agm2_Tdata_"+V_charm_1_M.Tag[i_ens]+"_ELM_M.dat.t", "", "#id  Tdata   ag2m agm2_err");
    //H
    Print_To_File({}, {Tdata_vec, agm2_distr_Tdata_H.ave(), agm2_distr_Tdata_H.err()}, "../data/gm2/charm/tm_"+Extrapolation_charm_mode+"/agm2_Tdata_"+V_charm_1_H.Tag[i_ens]+"_ELM_H.dat.t", "", "#id  Tdata   ag2m agm2_err");
    Print_To_File({}, {Tdata_vec, agm2_OS_distr_Tdata_H.ave(), agm2_OS_distr_Tdata_H.err()}, "../data/gm2/charm/OS_"+Extrapolation_charm_mode+"/agm2_Tdata_"+V_charm_1_H.Tag[i_ens]+"_ELM_H.dat.t", "", "#id  Tdata   ag2m agm2_err");

    //NON_ELM
    //L
    Print_To_File({}, {Tdata_vec, agm2_distr_Tdata_No_ELM_L.ave(), agm2_distr_Tdata_No_ELM_L.err()}, "../data/gm2/charm/tm_"+Extrapolation_charm_mode+"/agm2_Tdata_"+V_charm_1_L.Tag[i_ens]+"_L.dat.t", "", "#id  Tdata   ag2m agm2_err");
    Print_To_File({}, {Tdata_vec, agm2_OS_distr_Tdata_No_ELM_L.ave(), agm2_OS_distr_Tdata_No_ELM_L.err()}, "../data/gm2/charm/OS_"+Extrapolation_charm_mode+"/agm2_Tdata_"+V_charm_1_L.Tag[i_ens]+"_L.dat.t", "", "#id  Tdata   ag2m agm2_err");
    //M
    Print_To_File({}, {Tdata_vec, agm2_distr_Tdata_No_ELM_M.ave(), agm2_distr_Tdata_No_ELM_M.err()}, "../data/gm2/charm/tm_"+Extrapolation_charm_mode+"/agm2_Tdata_"+V_charm_1_M.Tag[i_ens]+"_M.dat.t", "", "#id  Tdata   ag2m agm2_err");
    Print_To_File({}, {Tdata_vec, agm2_OS_distr_Tdata_No_ELM_M.ave(), agm2_OS_distr_Tdata_No_ELM_M.err()}, "../data/gm2/charm/OS_"+Extrapolation_charm_mode+"/agm2_Tdata_"+V_charm_1_M.Tag[i_ens]+"_M.dat.t", "", "#id  Tdata   ag2m agm2_err");
    //H
    Print_To_File({}, {Tdata_vec, agm2_distr_Tdata_No_ELM_H.ave(), agm2_distr_Tdata_No_ELM_H.err()}, "../data/gm2/charm/tm_"+Extrapolation_charm_mode+"/agm2_Tdata_"+V_charm_1_H.Tag[i_ens]+"_H.dat.t", "", "#id  Tdata   ag2m agm2_err");
    Print_To_File({}, {Tdata_vec, agm2_OS_distr_Tdata_No_ELM_H.ave(), agm2_OS_distr_Tdata_No_ELM_H.err()}, "../data/gm2/charm/OS_"+Extrapolation_charm_mode+"/agm2_Tdata_"+V_charm_1_H.Tag[i_ens]+"_H.dat.t", "", "#id  Tdata   ag2m agm2_err");






    
    //######################   PI(Q^2) analysis   ###########################


    distr_t_list PI_Q2_L_tm, PI_Q2_L_OS, PI_Q2_M_tm, PI_Q2_M_OS, PI_Q2_H_tm, PI_Q2_H_OS;
    distr_t_list PI_Q2_L_tm_pert_sub, PI_Q2_L_OS_pert_sub, PI_Q2_M_tm_pert_sub, PI_Q2_M_OS_pert_sub, PI_Q2_H_tm_pert_sub, PI_Q2_H_OS_pert_sub;

    
    Vint Tdatas_opt_L_tm, Tdatas_opt_L_OS, Tdatas_opt_M_tm, Tdatas_opt_M_OS, Tdatas_opt_H_tm, Tdatas_opt_H_OS;
    Vint Tdatas_opt_L_tm_pert_sub, Tdatas_opt_L_OS_pert_sub, Tdatas_opt_M_tm_pert_sub, Tdatas_opt_M_OS_pert_sub, Tdatas_opt_H_tm_pert_sub, Tdatas_opt_H_OS_pert_sub;

    

    // Apply bounding method
    Bounding_PI_q2(PI_Q2_L_tm, Za_L*Za_L*V_charm_distr_L, a_distr, "../data/PI_Q2/charm/tm_"+Extrapolation_charm_mode+"/PI_Q2_Tdata_"+V_charm_1_L.Tag[i_ens]+"_L"  , Tdatas_opt_L_tm, MV_charm_L);
    Bounding_PI_q2(PI_Q2_M_tm, Za_M*Za_M*V_charm_distr_M, a_distr, "../data/PI_Q2/charm/tm_"+Extrapolation_charm_mode+"/PI_Q2_Tdata_"+V_charm_1_M.Tag[i_ens]+"_M"  , Tdatas_opt_M_tm, MV_charm_M);
    Bounding_PI_q2(PI_Q2_H_tm, Za_H*Za_H*V_charm_distr_H, a_distr, "../data/PI_Q2/charm/tm_"+Extrapolation_charm_mode+"/PI_Q2_Tdata_"+V_charm_1_H.Tag[i_ens]+"_H"  , Tdatas_opt_H_tm, MV_charm_H);
    Bounding_PI_q2(PI_Q2_L_OS, Zv_L*Zv_L*V_charm_OS_distr_L, a_distr, "../data/PI_Q2/charm/OS_"+Extrapolation_charm_mode+"/PI_Q2_Tdata_"+V_charm_1_L.Tag[i_ens]+"_L"  , Tdatas_opt_L_OS, MV_charm_OS_L);
    Bounding_PI_q2(PI_Q2_M_OS, Zv_M*Zv_M*V_charm_OS_distr_M, a_distr, "../data/PI_Q2/charm/OS_"+Extrapolation_charm_mode+"/PI_Q2_Tdata_"+V_charm_1_M.Tag[i_ens]+"_M"  , Tdatas_opt_M_OS, MV_charm_OS_M);
    Bounding_PI_q2(PI_Q2_H_OS, Zv_H*Zv_H*V_charm_OS_distr_H, a_distr, "../data/PI_Q2/charm/OS_"+Extrapolation_charm_mode+"/PI_Q2_Tdata_"+V_charm_1_H.Tag[i_ens]+"_H"  , Tdatas_opt_H_OS, MV_charm_OS_H);

  
    //include perturbative subtraction
    Get_PI_q2(PI_Q2_L_tm_pert_sub, VV_free_oppor_L, a_distr, Corr.Nt/2);
    Get_PI_q2(PI_Q2_M_tm_pert_sub, VV_free_oppor_M, a_distr, Corr.Nt/2);
    Get_PI_q2(PI_Q2_H_tm_pert_sub, VV_free_oppor_H, a_distr, Corr.Nt/2);
    Get_PI_q2(PI_Q2_L_OS_pert_sub, VV_free_samer_L, a_distr, Corr.Nt/2);
    Get_PI_q2(PI_Q2_M_OS_pert_sub, VV_free_samer_M, a_distr, Corr.Nt/2);
    Get_PI_q2(PI_Q2_H_OS_pert_sub, VV_free_samer_H, a_distr, Corr.Nt/2);

    PI_Q2_L_tm_pert_sub = PI_Q2_L_tm_pert_sub + PI_Q2_L_tm;
    PI_Q2_M_tm_pert_sub = PI_Q2_M_tm_pert_sub + PI_Q2_M_tm;
    PI_Q2_H_tm_pert_sub = PI_Q2_H_tm_pert_sub + PI_Q2_H_tm;
    PI_Q2_L_OS_pert_sub = PI_Q2_L_OS_pert_sub + PI_Q2_L_OS;
    PI_Q2_M_OS_pert_sub = PI_Q2_M_OS_pert_sub + PI_Q2_M_OS;
    PI_Q2_H_OS_pert_sub = PI_Q2_H_OS_pert_sub + PI_Q2_H_OS;



     if(V_charm_1_M.Tag[i_ens] == "cD211a.054.96" || V_charm_1_M.Tag[i_ens] == "cC211a.06.80" || V_charm_1_M.Tag[i_ens] == "cB211b.072.64") {
       int id_disco_ens=0;
       if(V_charm_1_M.Tag[i_ens]== "cD211a.054.96") id_disco_ens=2;
       else if(V_charm_1_M.Tag[i_ens]== "cC211a.06.80") id_disco_ens=1;
       else if(V_charm_1_M.Tag[i_ens]== "cB211b.072.64") id_disco_ens=0;
       else crash("what?");
       
       CORR_CHARM_FOR_PI_Q2[id_disco_ens]= Zv_H*Zv_H*V_charm_OS_distr_H;
       CHARM_PI_Q2_FOR_DISCO[id_disco_ens] = PI_Q2_H_OS;
       if(Include_charm_disco && Include_off_diagonal_disco) {
	 if(!Found_disco_impr_ens || !Found_disco_impr_lightD_charm_ens) crash("Disconnected for PI(Q^2) not found in charm ens");
	 CORR_DISCO_FOR_PI_Q2[id_disco_ens] = CORR_DISCO_FOR_PI_Q2[id_disco_ens] + Zv_M*Zv_M*( disco_impr_distr + disco_impr_lightD_charm_distr);

       }
       else crash("what?");

     }

    //interpolate to the physical point

    distr_t_list PI_Q2_Extr_tm, PI_Q2_Extr_OS;
    distr_t_list PI_Q2_Extr_tm_pert_sub, PI_Q2_Extr_OS_pert_sub;

    Vfloat Tcut_f_tm, Tcut_f_OS, Tcut_f_tm_pert_sub, Tcut_f_OS_pert_sub;

    
    int Q_size= Qs2.size();
    for(int q=0; q < Q_size;q++) {

      vector<distr_t> PIs_q2_charm_tm, PIs_q2_charm_OS, PIs_q2_charm_tm_pert_sub, PIs_q2_charm_OS_pert_sub;

      if(V_charm_1_L.Tag[i_ens].substr(1,1)=="D") {
	PIs_q2_charm_tm = { PI_Q2_L_tm.distr_list[q], PI_Q2_M_tm.distr_list[q]};
	PIs_q2_charm_OS = { PI_Q2_L_OS.distr_list[q], PI_Q2_M_OS.distr_list[q]};
	PIs_q2_charm_tm_pert_sub = { PI_Q2_L_tm_pert_sub.distr_list[q], PI_Q2_M_tm_pert_sub.distr_list[q]};
	PIs_q2_charm_OS_pert_sub = { PI_Q2_L_OS_pert_sub.distr_list[q], PI_Q2_M_OS_pert_sub.distr_list[q]};
      }
      else {

	PIs_q2_charm_tm = { PI_Q2_L_tm.distr_list[q], PI_Q2_M_tm.distr_list[q], PI_Q2_H_tm.distr_list[q] };
	PIs_q2_charm_OS = { PI_Q2_L_OS.distr_list[q], PI_Q2_M_OS.distr_list[q], PI_Q2_H_OS.distr_list[q]};
	PIs_q2_charm_tm_pert_sub = { PI_Q2_L_tm_pert_sub.distr_list[q], PI_Q2_M_tm_pert_sub.distr_list[q], PI_Q2_H_tm_pert_sub.distr_list[q]};
	PIs_q2_charm_OS_pert_sub = { PI_Q2_L_OS_pert_sub.distr_list[q], PI_Q2_M_OS_pert_sub.distr_list[q], PI_Q2_H_OS_pert_sub.distr_list[q] };

      }
      
      PI_Q2_Extr_tm.distr_list.push_back(Obs_extrapolation_meson_mass(PIs_q2_charm_tm, X_2_fit, X_2_phys,  "../data/PI_Q2/charm", "PI_Q2_tm_extr_"+Extrapolation_charm_mode+"_"+V_charm_1_L.Tag[i_ens]+"_Q2_"+to_string_with_precision(Qs2[q],5), UseJack, "SPLINE"));

      PI_Q2_Extr_OS.distr_list.push_back(Obs_extrapolation_meson_mass(PIs_q2_charm_OS, X_2_fit, X_2_phys,  "../data/PI_Q2/charm", "PI_Q2_OS_extr_"+Extrapolation_charm_mode+"_"+V_charm_1_L.Tag[i_ens]+"_Q2_"+to_string_with_precision(Qs2[q],5), UseJack, "SPLINE"));
      PI_Q2_Extr_tm_pert_sub.distr_list.push_back(Obs_extrapolation_meson_mass(PIs_q2_charm_tm_pert_sub, X_2_fit, X_2_phys,  "../data/PI_Q2/charm", "PI_Q2_tm_pert_sub_extr_"+Extrapolation_charm_mode+"_"+V_charm_1_L.Tag[i_ens]+"_Q2_"+to_string_with_precision(Qs2[q],5), UseJack, "SPLINE"));
      PI_Q2_Extr_OS_pert_sub.distr_list.push_back(Obs_extrapolation_meson_mass(PIs_q2_charm_OS_pert_sub, X_2_fit, X_2_phys,  "../data/PI_Q2/charm", "PI_Q2_OS_pert_sub_extr_"+Extrapolation_charm_mode+"_"+V_charm_1_L.Tag[i_ens]+"_Q2_"+to_string_with_precision(Qs2[q],5), UseJack, "SPLINE"));



      //check if bounding method worked for both L and M (1=worked, 0= not worked)
      if( Tdatas_opt_L_tm[q] > 0 && Tdatas_opt_M_tm[q] > 0 && Tdatas_opt_H_tm[q] > 0) Tcut_f_tm.push_back( 1.0);
      else Tcut_f_tm.push_back( 0.0);
      if( Tdatas_opt_L_OS[q] > 0 && Tdatas_opt_M_OS[q] > 0 && Tdatas_opt_H_OS[q] > 0) Tcut_f_OS.push_back( 1.0);
      else Tcut_f_OS.push_back( 0.0);
      //if( Tdatas_opt_L_tm_pert_sub[q] > 0 && Tdatas_opt_M_tm_pert_sub[q] > 0 && Tdatas_opt_H_tm_pert_sub[q] > 0) Tcut_f_tm_pert_sub.push_back( 1.0);
      //else Tcut_f_tm_pert_sub.push_back( 0.0);
      //if( Tdatas_opt_L_OS_pert_sub[q] > 0 && Tdatas_opt_M_OS_pert_sub[q] > 0 && Tdatas_opt_H_OS_pert_sub[q] > 0) Tcut_f_OS_pert_sub.push_back( 1.0);
      //else Tcut_f_OS_pert_sub.push_back( 0.0);

      
    }
    


    //push_back and print to file

    Add_ens_val_PI_q2( PI_Q2_charm_tm, PI_Q2_Extr_tm);
    Add_ens_val_PI_q2( PI_Q2_charm_OS, PI_Q2_Extr_OS);
    Add_ens_val_PI_q2( PI_Q2_charm_tm_pert_sub, PI_Q2_Extr_tm_pert_sub);
    Add_ens_val_PI_q2( PI_Q2_charm_OS_pert_sub, PI_Q2_Extr_OS_pert_sub);


    Print_To_File({}, {Qs2, PI_Q2_Extr_tm.ave(), PI_Q2_Extr_tm.err(), Tcut_f_tm } , "../data/PI_Q2/charm/tm_"+Extrapolation_charm_mode+"/PI_Q2_extr_"+V_charm_1_L.Tag[i_ens]+".t", "", "# Q2[GeV2]   PI(Q^2)   Tcut_f");
    Print_To_File({}, {Qs2, PI_Q2_Extr_OS.ave(), PI_Q2_Extr_OS.err(), Tcut_f_OS } , "../data/PI_Q2/charm/OS_"+Extrapolation_charm_mode+"/PI_Q2_extr_"+V_charm_1_L.Tag[i_ens]+".t", "", "# Q2[GeV2]   PI(Q^2)   Tcut_f");
    Print_To_File({}, {Qs2, PI_Q2_Extr_tm_pert_sub.ave(), PI_Q2_Extr_tm_pert_sub.err(), Tcut_f_tm } , "../data/PI_Q2/charm/tm_"+Extrapolation_charm_mode+"/PI_Q2_extr_pert_sub_"+V_charm_1_L.Tag[i_ens]+".t", "", "# Q2[GeV2]   PI(Q^2)   Tcut_f");
    Print_To_File({}, {Qs2, PI_Q2_Extr_OS_pert_sub.ave(), PI_Q2_Extr_OS_pert_sub.err(), Tcut_f_OS } , "../data/PI_Q2/charm/OS_"+Extrapolation_charm_mode+"/PI_Q2_extr_pert_sub_"+V_charm_1_L.Tag[i_ens]+".t", "", "# Q2[GeV2]   PI(Q^2)   Tcut_f");




    //#######################################################################

















    
  
  }


  cout<<"charm quark correlator analyzed!"<<endl;

 
  //light
  channel="l";
  

 
  
  for(int i_ens=0;i_ens<Nens_light;i_ens++) {
    
    CorrAnalysis Corr(UseJack, Njacks,Nboots);
    Corr.Nt = V_light_1.nrows[i_ens];

    //resample lattice spacing
    distr_t a_distr(UseJack), a_distr_2a(UseJack);
    LatticeInfo L_info;
    L_info.LatInfo_new_ens(V_light_1.Tag[i_ens]);
 
   
    distr_t_list  V_light_1_distr, V_light_2_distr, V_light_3_distr;
    
    //disconnected
    distr_t_list disco_distr;
    //improved
    distr_t_list disco_impr_distr, disco_impr_D_distr, disco_impr_DD_distr;

    //light1 light2 light4
    distr_t_list disco_impr_1_distr, disco_impr_2_distr, disco_impr_4_distr;

    //light no replica
    distr_t_list disco_impr_no_rep_distr;

   
    
									  
    distr_t_list  V_light_OS_1_distr, V_light_OS_2_distr, V_light_OS_3_distr, V_light_OS_distr;  
    distr_t_list  V_light_distr, MV_light_distr, ZV_light_distr;
    distr_t_list  MV_light_OS_distr, ZV_light_OS_distr;
    distr_t_list Mpi_distr, Mpi_OS_distr, fp_distr, overlap_P5P5_distr, overlap_P5P5_OS_distr, pion_corr;
    distr_t_list A0P5_distr, A0P5_OS_distr;
    distr_t MV_light, ZV_light, Mpi, Mpi_OS, fp;
    distr_t MV_light_OS, ZV_light_OS;
    distr_t_list P5P5_OS_distr,  RV, RA,RA0, ratio_P5P5_overlap_OS_tm, Zp_ov_Zs_distr;
    distr_t Zv, Za, Zp_ov_Zs, fp_ov_Mpi;
  

    //Analyze correlators
    if(V_light_1.Tag[i_ens].substr(1,1)=="A") {Corr.Tmin=14;Corr.Tmax=19; a_distr=a_A; a_distr_2a= a_A_afp;}
    else if(V_light_1.Tag[i_ens].substr(1,12)=="B211b.072.96") {Corr.Tmin=26; Corr.Tmax=60; a_distr=a_B; a_distr_2a= a_B_afp;}
    else if(V_light_1.Tag[i_ens].substr(1,12)=="B211b.072.64") { Corr.Tmin=24; Corr.Tmax=36; a_distr = a_B; a_distr_2a= a_B_afp;}
    else if(V_light_1.Tag[i_ens].substr(1,11)=="B211a.25.24") {Corr.Tmin=15; Corr.Tmax=20; a_distr=a_B; a_distr_2a= a_B_afp;}
    else if(V_light_1.Tag[i_ens].substr(1,11)=="B211a.25.32") {Corr.Tmin=21; Corr.Tmax=30; a_distr = a_B; a_distr_2a= a_B_afp;}
    else if(V_light_1.Tag[i_ens].substr(1,11)=="B211a.25.48") {Corr.Tmin=21; Corr.Tmax=40; a_distr = a_B; a_distr_2a= a_B_afp;}
    else if(V_light_1.Tag[i_ens].substr(1,1)=="C") {Corr.Tmin=40; Corr.Tmax=60; a_distr=a_C; a_distr_2a= a_C_afp;}
    else if(V_light_1.Tag[i_ens].substr(1,1)=="D") {Corr.Tmin=41; Corr.Tmax=80; a_distr=a_D; a_distr_2a= a_D_afp;} 
    else crash("Cannot find ensemble: "+V_light_1.Tag[i_ens]+" while modifying fit interval for MV_light");


    distr_t Zv_WI_strange_distr, Za_WI_strange_distr;

    if(UseJack) {
      for(int ijack=0;ijack<Njacks;ijack++) {
	Zv_WI_strange_distr.distr.push_back( L_info.Zv_WI_strange + GM()*L_info.Zv_WI_strange_err/(sqrt(Njacks -1.0)));
	Za_WI_strange_distr.distr.push_back( L_info.Za_WI_strange + GM()*L_info.Za_WI_strange_err/(sqrt(Njacks -1.0)));
      }
    }
    else {
      for(int iboot=0;iboot<Nboots;iboot++) {
	Zv_WI_strange_distr.distr.push_back( L_info.Zv_WI_strange + GM()*L_info.Zv_WI_strange_err);
	Za_WI_strange_distr.distr.push_back( L_info.Za_WI_strange + GM()*L_info.Za_WI_strange_err);
      }
    }
  



    

    //define lambda function to compute function of distributions
    auto sqr= [=](double a, double b) {return sqrt(a);};
    auto SINH= [](double m) -> double  {return sinh(m);};
    auto SINH2 = [](double m, double t) -> double {return sinh(m);};

  
    //vector light sector
    V_light_1_distr = Corr.corr_t(V_light_1.col(0)[i_ens], "../data/gm2/light/tm/corr_1_"+V_light_1.Tag[i_ens]+".dat");
    V_light_2_distr = Corr.corr_t(V_light_2.col(0)[i_ens], "../data/gm2/light/tm/corr_2_"+V_light_2.Tag[i_ens]+".dat");
    V_light_3_distr = Corr.corr_t(V_light_3.col(0)[i_ens], "../data/gm2/light/tm/corr_3_"+V_light_3.Tag[i_ens]+".dat");
    V_light_OS_1_distr = Corr.corr_t(V_light_OS_1.col(0)[i_ens], "../data/gm2/light/OS/corr_1_"+V_light_1.Tag[i_ens]+".dat");
    V_light_OS_2_distr = Corr.corr_t(V_light_OS_2.col(0)[i_ens], "../data/gm2/light/OS/corr_2_"+V_light_2.Tag[i_ens]+".dat");
    V_light_OS_3_distr = Corr.corr_t(V_light_OS_3.col(0)[i_ens], "../data/gm2/light/OS/corr_3_"+V_light_3.Tag[i_ens]+".dat");
    //sum over the Lorentz indices of the e.m. current
    V_light_distr= ((pow(qu,2)+pow(qd,2))/3.0)*(V_light_1_distr+ V_light_2_distr + V_light_3_distr);
    V_light_OS_distr =  ((pow(qu,2)+pow(qd,2))/3.0)*(V_light_OS_1_distr+ V_light_OS_2_distr + V_light_OS_3_distr);

   
    //JACK_TEST
    if(V_light_1.Tag[i_ens].substr(1,1)=="D") printV(V_light_1_distr.distr_list[20].distr, "JACK_TEST", 1);
    
  

    //Read disconnected 
    bool Found_disco_ens=false;
    bool Found_disco_impr_ens=false;
    bool Found_disco_impr_D_ens=false;
    bool Found_disco_impr_DD_ens=false;
    bool Found_disco_impr_1_ens=false;
    bool Found_disco_impr_2_ens=false;
    bool Found_disco_impr_4_ens=false;
    bool Found_disco_impr_no_rep_ens=false;
    
    if(Include_light_disco) {

      int i_ens_disco=0;
      
      for(int j=0;j<Nens_disco_light;j++) if(disco_light.Tag[j] == V_light_1.Tag[i_ens]) { Found_disco_ens=true; i_ens_disco=j;disco_light_Tags.push_back(disco_light.Tag[j]) ;break;}
      if(Found_disco_ens) {
	disco_distr = Corr.corr_t(disco_light.col(0)[i_ens_disco], "");
	disco_distr = disco_distr*(pow(qu+qd,2));
      }

      
      //improved
      for(int j=0;j<disco_impr_light.size;j++) if(disco_impr_light.Tag[j] == V_light_1.Tag[i_ens]) { Found_disco_impr_ens=true; i_ens_disco=j;disco_impr_light_Tags.push_back(disco_impr_light.Tag[j]) ;break;}
      if(Found_disco_impr_ens) {
	disco_impr_distr = Corr.corr_t(disco_impr_light.col(0)[i_ens_disco], "");
	disco_impr_distr = disco_impr_distr*(pow(qu+qd,2));
      }

      
      //improved D
      for(int j=0;j<disco_impr_lightD.size;j++) if(disco_impr_lightD.Tag[j] == V_light_1.Tag[i_ens]) { Found_disco_impr_D_ens=true; i_ens_disco=j;disco_impr_lightD_Tags.push_back(disco_impr_lightD.Tag[j]) ;break;}
      if(Found_disco_impr_D_ens) {
	disco_impr_D_distr = Corr.corr_t(disco_impr_lightD.col(0)[i_ens_disco], "");
	disco_impr_D_distr = disco_impr_D_distr*(pow(qu+qd,2));
      }

      
      //improved DD
      for(int j=0;j<disco_impr_lightDD.size;j++) if(disco_impr_lightDD.Tag[j] == V_light_1.Tag[i_ens]) { Found_disco_impr_DD_ens=true; i_ens_disco=j;disco_impr_lightDD_Tags.push_back(disco_impr_lightDD.Tag[j]) ;break;}
      if(Found_disco_impr_DD_ens) {
	disco_impr_DD_distr = Corr.corr_t(disco_impr_lightDD.col(0)[i_ens_disco], "");
	disco_impr_DD_distr = disco_impr_DD_distr*(pow(qu+qd,2));
      }

      //improved 1
      for(int j=0;j<disco_impr_light1.size ;j++) if(disco_impr_light1.Tag[j] == V_light_1.Tag[i_ens]) { Found_disco_impr_1_ens=true; i_ens_disco=j; ;break;}
      if(Found_disco_impr_1_ens) {
	disco_impr_1_distr = Corr.corr_t(disco_impr_light1.col(0)[i_ens_disco], "");
	disco_impr_1_distr = disco_impr_1_distr*(pow(qu+qd,2));
      }

      //improved 2
      for(int j=0;j<disco_impr_light2.size ;j++) if(disco_impr_light2.Tag[j] == V_light_1.Tag[i_ens]) { Found_disco_impr_2_ens=true; i_ens_disco=j; ;break;}
      if(Found_disco_impr_2_ens) {
	disco_impr_2_distr = Corr.corr_t(disco_impr_light2.col(0)[i_ens_disco], "");
	disco_impr_2_distr = disco_impr_2_distr*(pow(qu+qd,2));
      }

      //improved 4
      for(int j=0;j<disco_impr_light4.size ;j++) if(disco_impr_light4.Tag[j] == V_light_1.Tag[i_ens]) { Found_disco_impr_4_ens=true; i_ens_disco=j; ;break;}
      if(Found_disco_impr_4_ens) {
	disco_impr_4_distr = Corr.corr_t(disco_impr_light4.col(0)[i_ens_disco], "");
	disco_impr_4_distr = disco_impr_4_distr*(pow(qu+qd,2));
      }

      //improved no_replica
      for(int j=0;j<disco_impr_light_no_rep.size ;j++) if(disco_impr_light_no_rep.Tag[j] == V_light_1.Tag[i_ens]) { Found_disco_impr_no_rep_ens=true; i_ens_disco=j; ;break;}
      if(Found_disco_impr_no_rep_ens) {
	disco_impr_no_rep_distr = Corr.corr_t(disco_impr_light_no_rep.col(0)[i_ens_disco], "");
	disco_impr_no_rep_distr = disco_impr_no_rep_distr*(pow(qu+qd,2));
      }
      
    }
    





    
    MV_light_distr= Corr.effective_mass_t(V_light_distr, "../data/gm2/light/tm/MV_"+V_light_1.Tag[i_ens]+".dat");
    ZV_light_distr= Corr.residue_t(V_light_distr, "../data/gm2/light/tm/ZV_overlap_"+V_light_1.Tag[i_ens]+"dat");
    MV_light_OS_distr= Corr.effective_mass_t(V_light_OS_distr, "../data/gm2/light/OS/MV_"+V_light_1.Tag[i_ens]+".dat");
    ZV_light_OS_distr= Corr.residue_t(V_light_OS_distr, "../data/gm2/light/OS/ZV_overlap_"+V_light_1.Tag[i_ens]+"dat");
  


    //tm pion sector
    pion_corr= Corr.corr_t(pt2_pion.col(0)[i_ens], "");
    Mpi_distr= Corr.effective_mass_t(pion_corr, "../data/gm2/light/Mpi_"+V_light_1.Tag[i_ens]+".dat");
    overlap_P5P5_distr = Corr.residue_t(pion_corr, "");
    fp_distr = Corr.decay_constant_t(pow(2.0*L_info.ml,2)*pion_corr, "../data/gm2/light/fp_"+V_light_1.Tag[i_ens]+".dat");


    //OS sector
    P5P5_OS_distr = Corr.corr_t(corr_P5P5_OS.col(0)[i_ens], "");
    Mpi_OS_distr= Corr.effective_mass_t(P5P5_OS_distr, "../data/gm2/light/Mpi_OS_"+V_light_1.Tag[i_ens]+".dat");
    overlap_P5P5_OS_distr= Corr.residue_t(P5P5_OS_distr, "");



    //take ratio between OS and tm pion amplitude to compute Zp/Zs RC.
    ratio_P5P5_overlap_OS_tm= overlap_P5P5_OS_distr/overlap_P5P5_distr;
    Zp_ov_Zs_distr = distr_t_list::f_of_distr_list(sqr, ratio_P5P5_overlap_OS_tm);



    //antysymmetrize w.r.t. t -> T-t for A0P5 correlators
    //Corr.Reflection_sign = -1;
    A0P5_distr= Corr.corr_t(corr_A0P5.col(0)[i_ens], "");
    A0P5_OS_distr = Corr.corr_t(corr_A0P5_OS.col(0)[i_ens], "");

  
    //restore symmetrization
    // Corr.Reflection_sign = 1;

    //compute RV (estimator for Zv)


    RV= 2.0*L_info.ml*pion_corr/distr_t_list::derivative(A0P5_distr, 0); //central derivative

 
  

    //extract effective masses, overlap from V and fit
    int Tmin_old=Corr.Tmin;
    int Tmax_old=Corr.Tmax;
    //define restricted time intervals since Rho plateaux are tricky
    if(V_light_1.Tag[i_ens].substr(1,1)=="C" || V_light_1.Tag[i_ens].substr(1,1) == "D") {Corr.Tmin=14;Corr.Tmax=19;}
    else if(V_light_1.Tag[i_ens].substr(1,1)=="B" || V_light_1.Tag[i_ens].substr(1,1) =="A") {Corr.Tmin=16; Corr.Tmax=22;}
    else crash("Cannot find ensemble: "+V_light_1.Tag[i_ens]+" while modifying fit interval for MV_light");
    //################################################

  
    MV_light = Corr.Fit_distr(MV_light_distr);
    ZV_light = Corr.Fit_distr(ZV_light_distr);
    MV_light_OS = Corr.Fit_distr(MV_light_OS_distr);
    ZV_light_OS = Corr.Fit_distr(ZV_light_OS_distr);
    //push_back MV_light and ZV_light
    MV_fit_light.distr_list.push_back(MV_light);
    ZV_fit_light.distr_list.push_back(ZV_light);
    MV_fit_light_OS.distr_list.push_back(MV_light_OS);
    ZV_fit_light_OS.distr_list.push_back(ZV_light_OS);


    //restore old time intervals
    Corr.Tmin=Tmin_old;
    Corr.Tmax=Tmax_old;
    //########################

  
    Mpi=Corr.Fit_distr(Mpi_distr);
    Mpi_OS= Corr.Fit_distr(Mpi_OS_distr);
    int id=0;
    if(V_light_1.Tag[i_ens] == "cB211b.072.64") id=0;
    else if(V_light_1.Tag[i_ens]== "cB211b.072.96") id=1;
    else if(V_light_1.Tag[i_ens].substr(1,1) == "C") id=2;
    else if(V_light_1.Tag[i_ens].substr(1,1) == "D") id=3;
    else crash("Ensemble not recognised");
    Mp_light_tm.distr_list[id] = Mpi;
    Mp_light_OS.distr_list[id] = Mpi_OS;
    fp= Corr.Fit_distr(fp_distr);
    fp_ov_Mpi= fp/Mpi;
    Zp_ov_Zs = Corr.Fit_distr(Zp_ov_Zs_distr);
    Zv= Corr.Fit_distr(RV);
    RA0 = (Mpi_OS_distr/Mpi_distr)*(distr_t_list::f_of_distr_list(SINH2, Mpi_OS_distr)/distr_t_list::f_of_distr_list(SINH2, Mpi_distr));
    RA = 2.0*L_info.ml*(P5P5_OS_distr/distr_t_list::derivative(A0P5_OS_distr, 0))*(Mpi_OS/Mpi)*(distr_t::f_of_distr(SINH, Mpi_OS)/distr_t::f_of_distr(SINH, Mpi))*(1.0/Zp_ov_Zs);

    //set plateaux for RA
  
    //Analyze correlators
    if(V_light_1.Tag[i_ens].substr(1,1)=="A") {Corr.Tmin=14;Corr.Tmax=19; a_distr=a_A;}
    else if(V_light_1.Tag[i_ens].substr(1,12)=="B211b.072.96") {Corr.Tmin=12; Corr.Tmax=35; a_distr=a_B;}
    else if(V_light_1.Tag[i_ens].substr(1,12)=="B211b.072.64") { Corr.Tmin=15; Corr.Tmax=38; a_distr = a_B;}
    else if(V_light_1.Tag[i_ens].substr(1,11)=="B211a.25.24") {Corr.Tmin=12; Corr.Tmax=22; a_distr=a_B;}
    else if(V_light_1.Tag[i_ens].substr(1,11)=="B211a.25.32") {Corr.Tmin=19; Corr.Tmax=29; a_distr = a_B;}
    else if(V_light_1.Tag[i_ens].substr(1,11)=="B211a.25.48") {Corr.Tmin=16; Corr.Tmax=32; a_distr = a_B;}
    else if(V_light_1.Tag[i_ens].substr(1,1)=="C") {Corr.Tmin=21; Corr.Tmax=43; a_distr=a_C;}
    else if(V_light_1.Tag[i_ens].substr(1,1)=="D") {Corr.Tmin=28; Corr.Tmax=51; a_distr=a_D;}
    else crash("Cannot find ensemble: "+V_light_1.Tag[i_ens]+" while modifying fit interval for MV_light");
    Za = Corr.Fit_distr(RA);


    //if Use_Zv_Za from strange correlator modify Za and Zv
    if(Use_Za_Zv_from_strange_run_in_light_corr) { Zv = Zv_WI_strange_distr; Za=Za_WI_strange_distr;}

    a_distr_list.distr_list.push_back(a_distr);
    //#################################################################



    //################### GS PREDICTION FOR FSE #######################

    distr_t Mpi_0(UseJack);
    distr_t Mpi_ave_OS_neutral(UseJack);
    double c_neutral= 1.0/pow(fm_to_inv_Gev,2);
    for(int ijack=0; ijack<Njacks;ijack++) {
      Mpi_0.distr.push_back( sqrt( pow(0.135,2) - pow(a_distr.distr[ijack],2)*c_neutral));
      Mpi_ave_OS_neutral.distr.push_back( 0.5*(Mpi_OS.distr[ijack]/a_distr.distr[ijack] + Mpi_0.distr[ijack]));
    }
    
    if(V_light_1.Tag[i_ens] != "cB211b.072.96") {
      GS_artifacts_W_win_tm.push_back(Get_GS_artifacts_win( (0.135+Mpi_0.ave())/2.0, 5.46*fm_to_inv_Gev, LL));
      GS_artifacts_W_win_OS.push_back(Get_GS_artifacts_win( Mpi_ave_OS_neutral.ave(), 5.46*fm_to_inv_Gev, LL, 1));
    }




    //#############################VALENCE AND SEA PION MASS MISTUNING CORRECTIONS#############################################################
     //valence mistuning effects
    distr_t_list V_light_m_small_distr, V_light_m_big_distr, V_light_OS_m_small_distr, V_light_OS_m_big_distr;
    distr_t_list V_light_m_small_distr_bootstrap(0), V_light_m_big_distr_bootstrap(0), V_light_OS_m_small_distr_bootstrap(0), V_light_OS_m_big_distr_bootstrap(0);


    //find id for ensemble B64
    int i_ens_B64=-1;
    for(int i=0;i<(signed)V_light_1.Tag.size();i++) if( V_light_1.Tag[i] == "cB211b.072.64") { i_ens_B64=i;}
    if(i_ens_B64==-1) crash("Cannot find ensemble B64");
    
    CorrAnalysis Corr_boot(0, Njacks, 800, i_ens);
    CorrAnalysis Corr_boot_b64(0, Njacks, 800, i_ens_B64);
    Corr_boot_b64.Nt= 128;
    Corr_boot.Nt= V_light_1.nrows[i_ens];
    int seed_boot= (V_light_1.Tag[i_ens] != "cB211b.072.96")?Corr_boot.seed:Corr_boot_b64.seed;
    int Nt_chir= (V_light_1.Tag[i_ens] == "cB211b.072.96")?128:Corr.Nt;
    
    
    
    int id_ens_val=-1;
    

    //find valence ensemble id
    for(int i_ens_val=0; i_ens_val<(signed)V_light_1_m_small.Tag.size(); i_ens_val++) {
      if( V_light_1_m_small.Tag[i_ens_val] == V_light_1.Tag[i_ens]) id_ens_val=i_ens_val;
    }
    if(id_ens_val==-1) crash("In evaluating valence quark mass correction to vector correlator, cannot find ensemble: "+V_light_1.Tag[i_ens]);
      
    //valence
    V_light_m_small_distr = (pow(qu,2)+ pow(qd,2))*Corr.corr_t(V_light_1_m_small.col(0)[id_ens_val], "");
    V_light_m_big_distr = (pow(qu,2)+ pow(qd,2))*Corr.corr_t(V_light_1_m_big.col(0)[id_ens_val], "");
    V_light_OS_m_small_distr = (pow(qu,2)+ pow(qd,2))*Corr.corr_t(V_light_OS_1_m_small.col(0)[id_ens_val], "");
    V_light_OS_m_big_distr = (pow(qu,2)+ pow(qd,2))*Corr.corr_t(V_light_OS_1_m_big.col(0)[id_ens_val], "");
    
    //valence bootstrap
    V_light_m_small_distr_bootstrap = (pow(qu,2)+ pow(qd,2))*Corr_boot.corr_t(V_light_1_m_small.col(0)[id_ens_val], "");
    V_light_m_big_distr_bootstrap = (pow(qu,2)+ pow(qd,2))*Corr_boot.corr_t(V_light_1_m_big.col(0)[id_ens_val], "");
    V_light_OS_m_small_distr_bootstrap = (pow(qu,2)+ pow(qd,2))*Corr_boot.corr_t(V_light_OS_1_m_small.col(0)[id_ens_val], "");
    V_light_OS_m_big_distr_bootstrap = (pow(qu,2)+ pow(qd,2))*Corr_boot.corr_t(V_light_OS_1_m_big.col(0)[id_ens_val], "");
    
  


    //sea quark mistuning effect
    distr_t_list delta_corr_sea_tm, delta_corr_sea_OS;
    distr_t_list delta_corr_sea_tm_bootstrap(0), delta_corr_sea_OS_bootstrap(0);
   
    distr_t_list chiral_condensate_distr;

    int id_ens=-1;
    int id_ens_corr_to_bubble=-1;
    for(int i_ens_chir=0; i_ens_chir< (signed)chiral_condensate.Tag.size();i_ens_chir++) {
      if(V_light_1.Tag[i_ens]=="cB211b.072.96") {
	if(chiral_condensate.Tag[i_ens_chir] == "cB211b.072.64") id_ens= i_ens_chir;
      }
      else {
	if( chiral_condensate.Tag[i_ens_chir] == V_light_1.Tag[i_ens]) id_ens = i_ens_chir;
      }
    }
    if(id_ens==-1) crash("In evaluating sea quark mistuning: cannot find ensemble: "+V_light_1.Tag[i_ens]);
    
    for(int i_ens_chir=0; i_ens_chir< (signed)correlated_V_tm_bubble.Tag.size();i_ens_chir++) {
      if( correlated_V_tm_bubble.Tag[i_ens_chir] == V_light_1.Tag[i_ens]) id_ens_corr_to_bubble = i_ens_chir;
    }
    if(id_ens_corr_to_bubble==-1) crash("In evaluating sea quark mistuning (2): cannot find ensemble: "+V_light_1.Tag[i_ens]);
    
    int id_ens_corr_to_bubble_bis= -1;
    if(V_light_1.Tag[i_ens] == "cB211b.072.96") {
      for(int i=0;i<(signed)correlated_V_tm_bubble.Tag.size();i++) {
	if( correlated_V_tm_bubble.Tag[i] == "cB211b.072.64") id_ens_corr_to_bubble_bis= i;
      }
    }
    else id_ens_corr_to_bubble_bis=id_ens_corr_to_bubble;
    
    if(id_ens_corr_to_bubble_bis==-1) crash("In evaluating sea quark mistuning (3): cannot find ensemble: "+V_light_1.Tag[i_ens]);
    
    double correction_fact=1.0;
    if(V_light_1.Tag[i_ens] == "cB211b.072.64" || V_light_1.Tag[i_ens] == "cB211b.072.96") correction_fact= (0.00072-0.0006675);
    else if(V_light_1.Tag[i_ens] == "cC211a.06.80") correction_fact= (0.00060 - 0.000585  );
    else if(V_light_1.Tag[i_ens] == "cD211a.054.96") correction_fact= (0.00054 - 0.0004964);
    else crash("correction fact in sea quark mistuning, cannot be determined");
    


    //sea effects
    //chiral_condensate_distr= -1.0*Corr.corr_t( chiral_condensate.col(1)[id_ens], "../data/gm2/light/mass_dep/chiral_cond_"+V_light_1.Tag[i_ens]+".t");
    
    VVfloat chiral_condensate_correlated_to_V_tm(Nt_chir);
    VVfloat chiral_condensate_correlated_to_V_OS(Nt_chir);
    
    int Nconfs_tm_0 = correlated_V_tm_bubble.col(0)[id_ens_corr_to_bubble_bis][0].size();
    int Nconfs_OS_0 = correlated_V_OS_bubble.col(0)[id_ens_corr_to_bubble_bis][0].size();
    int Nconfs_chir_0 = chiral_condensate.col(1)[id_ens][0].size();
    if((Nconfs_tm_0 != Nconfs_OS_0) || (Nconfs_tm_0 != Nconfs_chir_0)) crash("In evaluating sea quark mistuning: number of confs between connected(tm,OS) and chiral condensate do not match"); 
    Vfloat chiral_condensate_averaged(Nconfs_chir_0,0.0);
    for(int t=0;t<Nt_chir;t++) {
      int Nconfs_chir = chiral_condensate.col(1)[id_ens][t].size();
      if( (Nconfs_chir != Nconfs_chir_0)) crash("In evaluating sea quark mistuning: Nconfs chiral not constant over time");
      for(int iconf=0; iconf<Nconfs_chir;iconf++) chiral_condensate_averaged[iconf] += -1.0*correction_fact*chiral_condensate.col(1)[id_ens][t][iconf];
    }
    
    auto F_disco_jack= [&](const Vfloat& par) { if((signed)par.size() != 3) crash("Lambda function F_disco_jack expects par[3], but par["+to_string((signed)par.size())+"] provided"); return par[0] -par[1]*par[2];};
    
    for(int t=0; t < Nt_chir;t++) {
      
      int Nconfs_tm = correlated_V_tm_bubble.col(0)[id_ens_corr_to_bubble_bis][t].size();
      int Nconfs_OS = correlated_V_OS_bubble.col(0)[id_ens_corr_to_bubble_bis][t].size();
      
      if( (Nconfs_tm != Nconfs_tm_0)) crash("In evaluating sea quark mistuning: Nconfs tm not constant over time");
      if( (Nconfs_OS != Nconfs_OS_0)) crash("In evaluating sea quark mistuning: Nconfs OS not constant over time");
      
      for(int iconf=0; iconf<Nconfs_tm_0;iconf++) {
	chiral_condensate_correlated_to_V_tm[t].push_back( correlated_V_tm_bubble.col(0)[id_ens_corr_to_bubble_bis][t][iconf]*chiral_condensate_averaged[iconf]);
	chiral_condensate_correlated_to_V_OS[t].push_back( correlated_V_OS_bubble.col(0)[id_ens_corr_to_bubble_bis][t][iconf]*chiral_condensate_averaged[iconf]);
      }
      
      
      
      //jackknife
      Jackknife J_tm(10000,Njacks);
      delta_corr_sea_tm.distr_list.push_back(J_tm.DoJack(F_disco_jack, 3, chiral_condensate_correlated_to_V_tm[t], chiral_condensate_averaged, correlated_V_tm_bubble.col(0)[id_ens_corr_to_bubble_bis][t]));
      Jackknife J_OS(10000,Njacks);
      delta_corr_sea_OS.distr_list.push_back(J_OS.DoJack(F_disco_jack, 3, chiral_condensate_correlated_to_V_OS[t], chiral_condensate_averaged, correlated_V_OS_bubble.col(0)[id_ens_corr_to_bubble_bis][t]));

      //bootstrap
      Bootstrap B_tm(800, seed_boot, chiral_condensate_averaged.size());
      delta_corr_sea_tm_bootstrap.distr_list.push_back(B_tm.DoBoot(F_disco_jack, 3, chiral_condensate_correlated_to_V_tm[t], chiral_condensate_averaged, correlated_V_tm_bubble.col(0)[id_ens_corr_to_bubble_bis][t]));
      Bootstrap B_OS(800, seed_boot, chiral_condensate_averaged.size() );
      delta_corr_sea_OS_bootstrap.distr_list.push_back(B_OS.DoBoot(F_disco_jack, 3, chiral_condensate_correlated_to_V_OS[t], chiral_condensate_averaged, correlated_V_OS_bubble.col(0)[id_ens_corr_to_bubble_bis][t]));
            	
    }


    if(V_light_1.Tag[i_ens]=="cB211b.072.96" ) {

      cout<<"I am here"<<endl;
      cout<<"size of delta corr before resizing: "<<delta_corr_sea_tm.size()<<endl;
      cout<<"size of delta corr before resizing: "<<delta_corr_sea_OS.size()<<endl;
      cout<<"size of delta corr before resizing: "<<delta_corr_sea_tm_bootstrap.size()<<endl;
      cout<<"size of delta corr before resizing: "<<delta_corr_sea_OS_bootstrap.size()<<endl;

      //resize delta_corr_sea_tm and delta_corr_sea_OS for jackknife and bootstrap
      distr_t fake_jack(1);
      distr_t fake_boot(0);
      for(int ijack=0;ijack<Njacks;ijack++) fake_jack.distr.push_back(0.0);
      for(int ib=0;ib<800;ib++) fake_boot.distr.push_back(0.0);


      for(int t=0;t<Corr.Nt;t++) {


	if(t<64) { //Do Nothing
	}
	else if(t<128) {
	  delta_corr_sea_tm.distr_list[t] = fake_jack;
	  delta_corr_sea_OS.distr_list[t] = fake_jack;
	  delta_corr_sea_tm_bootstrap.distr_list[t]= fake_boot;
	  delta_corr_sea_OS_bootstrap.distr_list[t]= fake_boot;
	}
	else {
	  delta_corr_sea_tm.distr_list.push_back(fake_jack);
	  delta_corr_sea_OS.distr_list.push_back(fake_jack);
	  delta_corr_sea_tm_bootstrap.distr_list.push_back(fake_boot);
	  delta_corr_sea_OS_bootstrap.distr_list.push_back(fake_boot);
	}
      }

      //now symmetrize
      for(int t=0;t<Corr.Nt;t++) {
	delta_corr_sea_tm.distr_list[t] = delta_corr_sea_tm.distr_list[t] + delta_corr_sea_tm.distr_list[ (Corr.Nt-t)%Corr.Nt];
	delta_corr_sea_OS.distr_list[t] = delta_corr_sea_OS.distr_list[t] + delta_corr_sea_OS.distr_list[ (Corr.Nt-t)%Corr.Nt];
	delta_corr_sea_tm_bootstrap.distr_list[t] = delta_corr_sea_tm_bootstrap.distr_list[t] + delta_corr_sea_tm_bootstrap.distr_list[ (Corr.Nt-t)%Corr.Nt];
	delta_corr_sea_OS_bootstrap.distr_list[t] = delta_corr_sea_OS_bootstrap.distr_list[t] + delta_corr_sea_OS_bootstrap.distr_list[ (Corr.Nt-t)%Corr.Nt];
      }
      
      

      cout<<"size of delta corr after resizing: "<<delta_corr_sea_tm.size()<<endl;
    }
  


   
    
    Jackknife JJ(10000, Njacks);
    //cout<<"chiral cond ave for ens: "<<V_light_1.Tag[i_ens]<<" : "<<JJ.DoJack(1, chiral_condensate_averaged).ave()<<" +- "<<JJ.DoJack(1,chiral_condensate_averaged).err()<<endl;

    delta_corr_sea_tm= (pow(qu,2)+pow(qd,2))*delta_corr_sea_tm;
    delta_corr_sea_OS= (pow(qu,2)+pow(qd,2))*delta_corr_sea_OS;
    delta_corr_sea_tm_bootstrap= (pow(qu,2)+ pow(qd,2))*delta_corr_sea_tm_bootstrap;
    delta_corr_sea_OS_bootstrap= (pow(qu,2) + pow(qd,2))*delta_corr_sea_OS_bootstrap;


    distr_t_list V_uncorrected_corr_to_bubble_tm= (pow(qu,2)+pow(qd,2))*Corr.corr_t(correlated_V_tm_bubble.col(0)[id_ens_corr_to_bubble],"");
    distr_t_list V_uncorrected_corr_to_bubble_OS= (pow(qu,2)+pow(qd,2))*Corr.corr_t(correlated_V_OS_bubble.col(0)[id_ens_corr_to_bubble],"");
    distr_t_list V_uncorrected_corr_to_bubble_tm_bootstrap=  (pow(qu,2)+pow(qd,2))*Corr_boot.corr_t(correlated_V_tm_bubble.col(0)[id_ens_corr_to_bubble],"");
    distr_t_list V_uncorrected_corr_to_bubble_OS_bootstrap= (pow(qu,2)+pow(qd,2))*Corr_boot.corr_t(correlated_V_OS_bubble.col(0)[id_ens_corr_to_bubble],"");
    
    /*
      if(V_light_1.Tag[i_ens] == "cB211b.072.64") {
      //load B96
      int n_id_B96=-1;
      for(int idB=0; idB<(signed)V_light_1.Tag.size();idB++) if( V_light_1.Tag[idB]== "cB211b.072.96") n_id_B96= idB;
      if(n_id_B96==-1) crash("Cannot find ensemble B96 when interpolating B64 correlator at L ~ 5.46 fm");
      CorrAnalysis Corr_B96(0, Njacks, 800, n_id_B96);
      Corr_B96.Nt= V_light_1.nrows[n_id_B96];
      distr_t_list corr_B96_uncorr_tm= (pow(qu,2)+pow(qd,2))*Corr_B96.corr_t(V_light_1.col(0)[n_id_B96], "");
      distr_t_list corr_B96_uncorr_OS= (pow(qu,2)+pow(qd,2))*Corr_B96.corr_t(V_light_OS_1.col(0)[n_id_B96], "");
      distr_t_list corr_B64_interpolated_to_546_tm(0);
      distr_t_list corr_B64_interpolated_to_546_OS(0);
      
      
	for(int t=0;t<Corr.Nt;t++) {
	distr_t c96_t_tm= corr_B96_uncorr_tm.distr_list[t];
	distr_t c64_t_tm= V_uncorrected_corr_to_bubble_tm_bootstrap.distr_list[t];
	distr_t c96_t_OS= corr_B96_uncorr_OS.distr_list[t];
	distr_t c64_t_OS= V_uncorrected_corr_to_bubble_OS_bootstrap.distr_list[t];
	
	distr_t tm_interpol;
	distr_t OS_interpol;

	//interpolate
	for(int ijack=0;ijack<800;ijack++) {
	double exp_2_ML_B64= exp( -1.0*Mpi.ave()*64);
	double exp_2_ML_B96= exp(-1.0*Mpi.ave()*96);
	double exp_2_ML_extr= exp(-1.0*Mpi.ave()*5.46*fm_to_inv_Gev/a_distr.ave());
	double c96_t_tm_jk= c96_t_tm.distr[ijack];
	double c64_t_tm_jk= c64_t_tm.distr[ijack];
	double c96_t_OS_jk= c96_t_OS.distr[ijack];
	double c64_t_OS_jk= c64_t_OS.distr[ijack];
	double amu_Linf_tm = (c96_t_tm_jk*exp_2_ML_B64 - c64_t_tm_jk*exp_2_ML_B96)/(exp_2_ML_B64 - exp_2_ML_B96);
	double bmu_Linf_tm = (c64_t_tm_jk - amu_Linf_tm)/exp_2_ML_B64;
	double amu_Linf_OS = (c96_t_OS_jk*exp_2_ML_B64 - c64_t_OS_jk*exp_2_ML_B96)/(exp_2_ML_B64 - exp_2_ML_B96);
	double bmu_Linf_OS = (c64_t_OS_jk - amu_Linf_OS)/exp_2_ML_B64;
	tm_interpol.distr.push_back( amu_Linf_tm+ bmu_Linf_tm*exp_2_ML_extr);
	OS_interpol.distr.push_back( amu_Linf_OS+ bmu_Linf_OS*exp_2_ML_extr);
	}
	
	corr_B64_interpolated_to_546_tm.distr_list.push_back( tm_interpol);
	corr_B64_interpolated_to_546_OS.distr_list.push_back( OS_interpol);
	}
	
	
	//V_uncorrected_corr_to_bubble_tm_bootstrap = corr_B64_interpolated_to_546_tm;
	//V_uncorrected_corr_to_bubble_OS_bootstrap = corr_B64_interpolated_to_546_OS;
	}
      */

      distr_t_list correlator_bare_corrected_tm_bootstrap = delta_corr_sea_tm_bootstrap + V_uncorrected_corr_to_bubble_tm_bootstrap ;
      distr_t_list correlator_bare_corrected_OS_bootstrap = delta_corr_sea_OS_bootstrap + V_uncorrected_corr_to_bubble_OS_bootstrap ;
      correlator_bare_corrected_tm_bootstrap = correlator_bare_corrected_tm_bootstrap + V_light_m_small_distr_bootstrap - V_light_m_big_distr_bootstrap;
      correlator_bare_corrected_OS_bootstrap = correlator_bare_corrected_OS_bootstrap + V_light_OS_m_small_distr_bootstrap - V_light_OS_m_big_distr_bootstrap;

      

      //print jackknives for correlator_bare_correctred_tm

      boost::filesystem::create_directory("../data/gm2/light/mass_dep/tm/bootstrap");
      boost::filesystem::create_directory("../data/gm2/light/mass_dep/OS/bootstrap");
      ofstream Print_boot_tm("../data/gm2/light/mass_dep/tm/bootstrap/boot_"+V_light_1.Tag[i_ens]+".dat");
      ofstream Print_boot_OS("../data/gm2/light/mass_dep/OS/bootstrap/boot_"+V_light_1.Tag[i_ens]+".dat");
      Print_boot_tm<<"#Nboots=800     T="<<Corr.Nt<<endl;
      Print_boot_OS<<"#Nboots=800     T="<<Corr.Nt<<endl;
      Print_boot_tm.precision(10);
      Print_boot_OS.precision(10);
      for(int ib=0;ib<800;ib++) {
	
	for(int t=0;t<Corr.Nt;t++) {
	  Print_boot_tm<<correlator_bare_corrected_tm_bootstrap.distr_list[t].distr[ib]<<endl;
	  Print_boot_OS<<correlator_bare_corrected_OS_bootstrap.distr_list[t].distr[ib]<<endl;
	}
      }
      Print_boot_tm.close();
      Print_boot_OS.close();
      
      //Print to File
      Print_To_File({}, {delta_corr_sea_tm.ave(), delta_corr_sea_tm.err(), (delta_corr_sea_tm*a_distr/delta_corr_sea_tm).ave() }, "../data/gm2/light/mass_dep/tm/delta_corr_sea_"+V_light_1.Tag[i_ens]+".dat", "", "#t deltaC   a[GeV-1]");
      Print_To_File({}, {delta_corr_sea_OS.ave(), delta_corr_sea_OS.err(), (delta_corr_sea_OS*a_distr/delta_corr_sea_OS).ave()}, "../data/gm2/light/mass_dep/OS/delta_corr_sea_"+V_light_1.Tag[i_ens]+".dat", "", "t deltaC  a[GeV-1]");



    
    //####################################################################################################




    //free corr LO artifacts
    //################## READ FREE THEORY VECTOR-VECTOR CORRELATOR OPPOSITE R ####################################
    string Pt_free_oppor= "../Vkvk_cont/"+to_string(Corr.Nt/2)+"_m"+to_string_with_precision(L_info.ml,5)+"/OPPOR";
    Vfloat VV_free_oppor= Read_From_File(Pt_free_oppor, 1, 4);
    if(VV_free_oppor.size() != Corr.Nt) crash("Failed to read properly free VV correlator ml w opposite r");
    //################## READ FREE THEORY VECTOR-VECTOR CORRELATOR SAME R ####################################
    string Pt_free_samer= "../Vkvk_cont/"+to_string(Corr.Nt/2)+"_m"+to_string_with_precision(L_info.ml,5)+"/SAMER";
    Vfloat VV_free_samer= Read_From_File(Pt_free_samer, 1, 4);
    if(VV_free_samer.size() != Corr.Nt) crash("Failed to read properly free VV correlator ml  w same r");
    //Insert electric charges
    for( auto & OP:VV_free_oppor) OP *= pert_corr_light_on_off*(qu*qu + qd*qd);
    for( auto & SA:VV_free_samer) SA *= pert_corr_light_on_off*(qu*qu + qd*qd);
    

    

    Vfloat free_corr_log_art(Corr.Nt,0);
    for(int t=1;t<Corr.Nt;t++) {
      if(t*a_distr.ave() < add_pert_corr_light_up_to*fm_to_inv_Gev)  free_corr_log_art[t] = -1.0*(qu*qu +qd*qd)*(1.0/(2.0*M_PI*M_PI*pow(t,5)));

      //if(t==0 || t*a_distr.ave() > add_pert_corr_light_up_to*fm_to_inv_Gev) { VV_free_samer[t] =0; VV_free_oppor[t] = 0;}


    }

    /*VV_free_samer= free_corr_log_art;
      VV_free_oppor= free_corr_log_art; */

    auto F_exp_sigma= [](double a, double t, double NT) {return exp(-t*t*pow(a/(0.4*fm_to_inv_Gev),2));};
    distr_t_list Exp_sigma_pert = distr_t_list::f_of_distr(F_exp_sigma, a_distr, Corr.Nt);

    distr_t_list V_light_distr_tm_corr = (1.0/(Za*Za))*(Za*Za*V_light_distr + VV_free_oppor); //free_corr_log_art
    distr_t_list V_light_distr_OS_corr = (1.0/(Zv*Zv))*(Zv*Zv*V_light_OS_distr + VV_free_samer); //free_corr_log_art


    distr_t_list V_light_distr_tm_Mpi_corr= V_light_distr_tm_corr+ delta_corr_sea_tm + V_light_m_small_distr - V_light_m_big_distr;
    distr_t_list V_light_distr_OS_Mpi_corr= V_light_distr_OS_corr+ delta_corr_sea_OS + V_light_OS_m_small_distr - V_light_OS_m_big_distr;



    //interpolate V_light_distr_tm_corr and V_light_distr_OS_corr + pion mass mistuning effects
    //################################################################
    vector<boost::math::interpolators::cardinal_cubic_b_spline<double>> Vt3_light_tm_interpol_func;
    vector<boost::math::interpolators::cardinal_cubic_b_spline<double>> Vt3_light_OS_interpol_func;
    for(int ijack=0;ijack<Njacks;ijack++) {
      Vfloat Vt3_tm_ijack, Vt3_OS_ijack;
      for(int t=1;t< Corr.Nt;t++) { Vt3_tm_ijack.push_back( V_light_distr_tm_Mpi_corr.distr_list[t].distr[ijack]*pow(t,3)); Vt3_OS_ijack.push_back( V_light_distr_OS_Mpi_corr.distr_list[t].distr[ijack]*pow(t,3));}

      Vt3_light_tm_interpol_func.emplace_back( Vt3_tm_ijack.begin(), Vt3_tm_ijack.end(), 1.0, 1.0);
      Vt3_light_OS_interpol_func.emplace_back( Vt3_OS_ijack.begin(), Vt3_OS_ijack.end(), 1.0, 1.0);
    }

    //define lambda functions that return the jack-distribution of the interpolated correlators
    auto Vt3_light_tm_interpol = [&Vt3_light_tm_interpol_func](double t) -> distr_t {

				 distr_t ret;

				 for(int i=0;i<Njacks;i++) ret.distr.push_back( Vt3_light_tm_interpol_func[i](t));

				 return ret;

			       };


    auto Vt3_light_OS_interpol = [&Vt3_light_OS_interpol_func](double t) -> distr_t {

				 distr_t ret;

				 for(int i=0;i<Njacks;i++) ret.distr.push_back( Vt3_light_OS_interpol_func[i](t));

				 return ret;

			       };
    //#################################################################
    




    //Print correlation matrix
    Eigen::MatrixXd Corr_matrix_OS(Corr.Nt/2, Corr.Nt/2);
    Eigen::MatrixXd Corr_matrix_tm(Corr.Nt/2, Corr.Nt/2);
    Eigen::MatrixXd Cov_matrix_OS(Corr.Nt/2, Corr.Nt/2);
    Eigen::MatrixXd Cov_matrix_tm(Corr.Nt/2, Corr.Nt/2);
    Eigen::MatrixXd Cov_matrix_OS_pert_sub(Corr.Nt/2, Corr.Nt/2);
    Eigen::MatrixXd Cov_matrix_tm_pert_sub(Corr.Nt/2, Corr.Nt/2);

    for(int i=0;i<Corr.Nt/2;i++) {
      for(int j=0;j<Corr.Nt/2;j++) {
	Corr_matrix_tm(i,j) = (Za*Za*V_light_distr).distr_list[i]%(Za*Za*V_light_distr).distr_list[j]/((Za*Za*V_light_distr).err(i)*(Za*Za*V_light_distr).err(j));
	Corr_matrix_OS(i,j) = (Zv*Zv*V_light_OS_distr).distr_list[i]%(Zv*Zv*V_light_OS_distr).distr_list[j]/((Zv*Zv*V_light_OS_distr).err(i)*(Zv*Zv*V_light_OS_distr).err(j));
	Cov_matrix_tm(i,j) = (Za*Za*V_light_distr).distr_list[i]%(Za*Za*V_light_distr).distr_list[j];
	Cov_matrix_OS(i,j) = (Zv*Zv*V_light_OS_distr).distr_list[i]%(Zv*Zv*V_light_OS_distr).distr_list[j];
	Cov_matrix_tm_pert_sub(i,j) = (Za*Za*V_light_distr_tm_corr).distr_list[i]%(Za*Za*V_light_distr_tm_corr).distr_list[j];
	Cov_matrix_OS_pert_sub(i,j) = (Zv*Zv*V_light_distr_OS_corr).distr_list[i]%(Zv*Zv*V_light_distr_OS_corr).distr_list[j];
      }
    }

    
    ofstream Print_Corr_tm("../data/gm2/light/tm/corr_matrix_vec_corr_"+V_light_1.Tag[i_ens]);
    ofstream Print_Corr_OS("../data/gm2/light/OS/corr_matrix_vec_corr_"+V_light_1.Tag[i_ens]);
    ofstream Print_Cov_tm("../data/gm2/light/tm/cov_matrix_vec_corr_"+V_light_1.Tag[i_ens]);
    ofstream Print_Cov_OS("../data/gm2/light/OS/cov_matrix_vec_corr_"+V_light_1.Tag[i_ens]);
    ofstream Print_Cov_tm_pert_sub("../data/gm2/light/tm/cov_matrix_vec_corr_pert_sub_"+V_light_1.Tag[i_ens]);
    ofstream Print_Cov_OS_pert_sub("../data/gm2/light/OS/cov_matrix_vec_corr_pert_sub_"+V_light_1.Tag[i_ens]);


    Print_Corr_tm<<Corr_matrix_tm<<endl;
    Print_Corr_OS<<Corr_matrix_OS<<endl;
    Print_Cov_tm<<Cov_matrix_tm<<endl;
    Print_Cov_OS<<Cov_matrix_OS<<endl;
    Print_Cov_tm_pert_sub<<Cov_matrix_tm_pert_sub<<endl;
    Print_Cov_OS_pert_sub<<Cov_matrix_OS_pert_sub<<endl;


    Print_Corr_tm.close();
    Print_Corr_OS.close();
    Print_Cov_tm.close();
    Print_Cov_OS.close();
    Print_Cov_tm_pert_sub.close();
    Print_Cov_OS_pert_sub.close();
    


    

    


    
    //Print To File additional observables
    // print summed connected correlators to file
    Print_To_File({}, {V_light_distr.ave(), V_light_distr.err(), (Za*Za*V_light_distr).ave(), (Za*Za*V_light_distr).err(),  ( (Za*Za*V_light_distr- Zv*Zv*V_light_OS_distr)/(Za*Za*V_light_distr)).ave(),  ( (Za*Za*V_light_distr- Zv*Zv*V_light_OS_distr)/(Za*Za*V_light_distr)).err(), (Za*Za*V_light_distr_tm_corr).ave(), (Za*Za*V_light_distr_tm_corr).err()}, "../data/gm2/light/tm/corr_sum_"+V_light_1.Tag[i_ens]+".dat.t", "", "#time   V(t)^tm    V(t)^tm(renormalized)    DV(t)(tm-os/tm)  V(t)^tm(VV_free corrected) ");
    Print_To_File({}, {V_light_OS_distr.ave(), V_light_OS_distr.err(), (Zv*Zv*V_light_OS_distr).ave(), (Zv*Zv*V_light_OS_distr).err(), (Zv*Zv*V_light_distr_OS_corr).ave() , (Zv*Zv*V_light_distr_OS_corr).err()}, "../data/gm2/light/OS/corr_sum_"+V_light_1.Tag[i_ens]+".dat.t", "", "#time V(t)^OS   V(t)^OS(renormalized) V(t)^OS(VV_free_corrected)");
   
    Print_To_File({}, { (Za*Za*(V_light_m_small_distr-V_light_m_big_distr)).ave(), (Za*Za*(V_light_m_small_distr-V_light_m_big_distr)).err(),(Za*Za*V_light_m_small_distr).ave(), (Za*Za*V_light_m_small_distr).err(), (Za*Za*V_light_m_big_distr).ave(), (Za*Za*V_light_m_big_distr).err()}, "../data/gm2/light/mass_dep/tm/corrs_"+V_light_1.Tag[i_ens]+".dat.t", "", "#t diff  0.0006675 0.00072");
      Print_To_File({}, { (Zv*Zv*(V_light_OS_m_small_distr-V_light_OS_m_big_distr)).ave(), (Zv*Zv*(V_light_OS_m_small_distr-V_light_OS_m_big_distr)).err(),(Zv*Zv*V_light_OS_m_small_distr).ave(), (Zv*Zv*V_light_OS_m_small_distr).err(), (Zv*Zv*V_light_OS_m_big_distr).ave(), (Zv*Zv*V_light_OS_m_big_distr).err()}, "../data/gm2/light/mass_dep/OS/corrs_"+V_light_1.Tag[i_ens]+".dat.t", "", "#t diff  0.0006675 0.00072");
   
    

    //print disco light
    if(Include_light_disco && Found_disco_ens) Print_To_File({}, {disco_distr.ave(), disco_distr.err(), (Zv*Zv*disco_distr).ave(), (Zv*Zv*disco_distr).err()}, "../data/gm2/light/disco/disc_"+V_light_1.Tag[i_ens]+".dat.t","","# bare renormalized");
    //improved
    if(Include_light_disco && Found_disco_impr_ens) Print_To_File({}, {disco_impr_distr.ave(), disco_impr_distr.err(), (Zv*Zv*disco_impr_distr).ave(), (Zv*Zv*disco_impr_distr).err()}, "../data/gm2/light/disco/disc_impr_"+V_light_1.Tag[i_ens]+".dat.t","","# bare renormalized");
    //improved D
    if(Include_light_disco && Found_disco_impr_D_ens) Print_To_File({}, {disco_impr_D_distr.ave(), disco_impr_D_distr.err(), (Zv*Zv*disco_impr_D_distr).ave(), (Zv*Zv*disco_impr_D_distr).err()}, "../data/gm2/light/disco/disc_impr_D_"+V_light_1.Tag[i_ens]+".dat.t","","# bare renormalized");
    //improved DD
    if(Include_light_disco && Found_disco_impr_DD_ens) Print_To_File({}, {disco_impr_DD_distr.ave(), disco_impr_DD_distr.err(), (Zv*Zv*disco_impr_DD_distr).ave(), (Zv*Zv*disco_impr_DD_distr).err()}, "../data/gm2/light/disco/disc_impr_DD_"+V_light_1.Tag[i_ens]+".dat.t","","# bare renormalized");
    //improved 1
    if(Include_light_disco && Found_disco_impr_1_ens) Print_To_File({}, {disco_impr_1_distr.ave(), disco_impr_1_distr.err(), (Zv*Zv*disco_impr_1_distr).ave(), (Zv*Zv*disco_impr_1_distr).err()}, "../data/gm2/light/disco/disc_impr_1_"+V_light_1.Tag[i_ens]+".dat.t","","# bare renormalized");
     //improved 2
    if(Include_light_disco && Found_disco_impr_2_ens) Print_To_File({}, {disco_impr_2_distr.ave(), disco_impr_2_distr.err(), (Zv*Zv*disco_impr_2_distr).ave(), (Zv*Zv*disco_impr_2_distr).err()}, "../data/gm2/light/disco/disc_impr_2_"+V_light_1.Tag[i_ens]+".dat.t","","# bare renormalized");
     //improved 4
    if(Include_light_disco && Found_disco_impr_4_ens) Print_To_File({}, {disco_impr_4_distr.ave(), disco_impr_4_distr.err(), (Zv*Zv*disco_impr_4_distr).ave(), (Zv*Zv*disco_impr_4_distr).err()}, "../data/gm2/light/disco/disc_impr_4_"+V_light_1.Tag[i_ens]+".dat.t","","# bare renormalized");
    //improved no_replica
    if(Include_light_disco && Found_disco_impr_no_rep_ens) Print_To_File({}, {disco_impr_no_rep_distr.ave(), disco_impr_no_rep_distr.err(), (Zv*Zv*disco_impr_no_rep_distr).ave(), (Zv*Zv*disco_impr_no_rep_distr).err()}, "../data/gm2/light/disco/disc_impr_no_rep_"+V_light_1.Tag[i_ens]+".dat.t","","# bare renormalized");
    
    
    
    



    
    //print RV
    Print_To_File({}, {RV.ave(), RV.err()}, "../data/gm2/light/RV_"+V_light_1.Tag[i_ens]+".dat.t", "", "");
    //print RA
    Print_To_File({}, {RA.ave(), RA.err(), RA0.ave(), RA0.err()}, "../data/gm2/light/RA_"+V_light_1.Tag[i_ens]+".dat.t", "", "");
    //print Zp_ov_Zs
    Print_To_File({}, {Zp_ov_Zs_distr.ave(), Zp_ov_Zs_distr.err()}, "../data/gm2/light/Zp_ov_Zs_"+V_light_1.Tag[i_ens]+".dat.t", "", "");


    cout<<"//#########################ENSEMBLE INFO: "<<endl;
    cout<<"Ensemble Tag: "<<V_light_1.Tag[i_ens]<<endl; 
    cout<<"L: "<<L_info.L<<" T: "<<L_info.T<<endl;
    cout<<"lat spacing: "<<a_distr.ave()/fm_to_inv_Gev<<"+- "<<a_distr.err()/fm_to_inv_Gev<<" (fm)"<<endl;
    cout<<"lat spacing from Mpi: "<<Mpi.ave()/(m_pi*fm_to_inv_Gev)<<" +- "<<Mpi.err()/(m_pi*fm_to_inv_Gev)<<" (fm)"<<endl;
    cout<<"lat spacing from fp: "<<fp.ave()/(fp_phys*fm_to_inv_Gev)<<" +- "<<fp.err()/(fp_phys*fm_to_inv_Gev)<<" (fm)"<<endl;
    cout<<"ml: "<<L_info.ml<<endl;
    cout<<"Mpi: "<<Mpi.ave()<<" +- "<<Mpi.err()<<endl;
    cout<<"Mpi*L: "<<Mpi.ave()*L_info.L<<" +- "<<Mpi.err()*L_info.L<<endl;
    cout<<"Mpi OS: "<<Mpi_OS.ave()<<" +- "<<Mpi_OS.err()<<endl;
    cout<<"fp: "<<fp.ave()<<" +- "<<fp.err()<<endl;
    cout<<"fp/Mpi: "<<fp_ov_Mpi.ave()<<" +- "<<fp_ov_Mpi.err()<<" Phys point: "<<0.966<<endl;
    cout<<"Zp/Zs: "<<Zp_ov_Zs.ave()<<" +- "<<Zp_ov_Zs.err()<<endl;
    cout<<"Za: "<<Za.ave()<<" +- "<<Za.err()<<endl;
    cout<<"Zv: "<<Zv.ave()<<" +- "<<Zv.err()<<endl;
    cout<<"MV: "<<MV_light.ave()<<" +- "<<MV_light.err()<<endl;
 

  
    //push_back Mpi Mpi_OS and fp, Rcs a and L
    Mpi_fit.distr_list.push_back(Mpi);
    Mpi_OS_fit.distr_list.push_back(Mpi_OS);
    fp_fit.distr_list.push_back(fp);
    Za_fit.distr_list.push_back(Za);
    Zv_fit.distr_list.push_back(Zv);
    Zp_ov_Zs_fit.distr_list.push_back(Zp_ov_Zs);
    ml_list.push_back(L_info.ml);
    a_list.push_back(a_distr.ave()/fm_to_inv_Gev);
    L_list.push_back((double)L_info.L);

    //push_back Mpi, Mpi_OS, fp, Zv a and L if Include_light_disco && Found_disco
    if(Include_light_disco && Found_disco_ens) {
      Mpi_fit_disco.distr_list.push_back(Mpi);
      Mpi_OS_fit_disco.distr_list.push_back(Mpi_OS);
      fp_fit_disco.distr_list.push_back(fp);
      Zv_fit_disco.distr_list.push_back(Zv);
      ml_list_disco.push_back(L_info.ml);
      a_list_disco.push_back(a_distr.ave()/fm_to_inv_Gev);
      a_distr_list_disco_light.distr_list.push_back(a_distr);
      L_list_disco.push_back((double)L_info.L);
    }

    //push-back info on lattice spacing for improved disconnected
    if(Include_light_disco && Found_disco_impr_ens) {
      a_list_disco_impr.push_back( a_distr.ave()/fm_to_inv_Gev);
    }
    //push-back info on lattice spacing for improved D disconnected
    if(Include_light_disco && Found_disco_impr_D_ens) {
      a_list_disco_impr_D.push_back( a_distr.ave()/fm_to_inv_Gev);
    }
    //push-back info on lattice spacing for improved D disconnected
    if(Include_light_disco && Found_disco_impr_DD_ens) {
      a_list_disco_impr_DD.push_back( a_distr.ave()/fm_to_inv_Gev);
    }
  

  


    int Tdata_min= 8;
    int Tdata_max = Corr.Nt/2.0 -4;
    int Tdata_fit = (Corr.Tmax+Corr.Tmin)/2;
  
    //compute kernel distribution
    distr_t_list Kernel_distr = distr_t_list::f_of_distr(K,Mpi/m_pi, Upper_Limit_Time_Integral_light+1);
    //compute exp(-Mv*t) distribution
    distr_t_list exp_MVl = distr_t_list::f_of_distr(exp_MV, MV_light, Upper_Limit_Time_Integral_light+1);
    distr_t_list exp_MVl_OS = distr_t_list::f_of_distr(exp_MV, MV_light_OS, Upper_Limit_Time_Integral_light+1);


    //Print single-exponential prediction to file
    Print_To_File({}, {(exp_MVl*(ZV_light/(2.0*MV_light))).ave(), (exp_MVl*(ZV_light/(2.0*MV_light))).err()}, "../data/gm2/light/tm/corr_gsd_sum_"+V_light_1.Tag[i_ens]+".dat.t", "", "");
    Print_To_File({}, {(exp_MVl_OS*(ZV_light_OS/(2.0*MV_light_OS))).ave(), (exp_MVl_OS*(ZV_light_OS/(2.0*MV_light_OS))).err()}, "../data/gm2/light/OS/corr_gsd_sum_"+V_light_1.Tag[i_ens]+".dat.t", "", "");

    distr_t_list agm2_distr_Tdata(UseJack);
    distr_t_list agm2_distr_OS_Tdata(UseJack);
    distr_t_list agm2_distr_Tdata_lower(UseJack);
    distr_t_list agm2_distr_OS_Tdata_lower(UseJack);
    Vfloat Tdata_vec;

  
    for(int Tdata=Tdata_min;Tdata<Tdata_max;Tdata++) {
      //compute 4\pia^2 using lattice data up to Tdata (included)
   
      distr_t agm2(UseJack, UseJack?Njacks:Nboots); //constructor sets agm2 to zero by default
      distr_t agm2_OS(UseJack, UseJack?Njacks:Nboots); //constructor sets agm2_OS to zero by default
      distr_t agm2_low(UseJack, UseJack?Njacks:Nboots); //constructor sets agm2 to zero by default
      distr_t agm2_OS_low(UseJack, UseJack?Njacks:Nboots); //constructor sets agm2_OS to zero by default 
    
    
      for(int t=1;t<=Upper_Limit_Time_Integral_light;t++) {
	if(t<=Tdata) {
	  agm2 = agm2 + 4.0*w(t,Simps_ord)*pow(alpha,2)*V_light_distr.distr_list[t]*Kernel_distr.distr_list[t];
	  agm2_OS = agm2_OS + 4.0*w(t,Simps_ord)*pow(alpha,2)*V_light_OS_distr.distr_list[t]*Kernel_distr.distr_list[t];
	  agm2_low = agm2_low + 4.0*w(t,Simps_ord)*pow(alpha,2)*V_light_distr.distr_list[t]*Kernel_distr.distr_list[t];
	  agm2_OS_low = agm2_OS_low + 4.0*w(t,Simps_ord)*pow(alpha,2)*V_light_OS_distr.distr_list[t]*Kernel_distr.distr_list[t];
	}
	else {
	  agm2= agm2 + 4.0*w(t,Simps_ord)*pow(alpha,2)*(ZV_light/(2.0*MV_light))*exp_MVl.distr_list[t]*Kernel_distr.distr_list[t];
	  agm2_OS= agm2_OS + 4.0*w(t,Simps_ord)*pow(alpha,2)*(ZV_light_OS/(2.0*MV_light_OS))*exp_MVl_OS.distr_list[t]*Kernel_distr.distr_list[t];
	}
      }
      Tdata_vec.push_back((double)Tdata);
      agm2 = agm2*Za*Za;
      agm2_OS = agm2_OS*Zv*Zv;
      agm2_low = agm2_low*Za*Za;
      agm2_OS_low = agm2_OS_low*Zv*Zv;
      agm2_distr_Tdata.distr_list.push_back(agm2);
      agm2_distr_OS_Tdata.distr_list.push_back(agm2_OS);
      agm2_distr_Tdata_lower.distr_list.push_back(agm2_low);
      agm2_distr_OS_Tdata_lower.distr_list.push_back(agm2_OS_low);
      if(Tdata==Tdata_fit) {
	agm2_light.distr_list.push_back(agm2);
	agm2_light_OS.distr_list.push_back(agm2_OS);
      }
    }
    //print to file
    Print_To_File({}, {Tdata_vec, agm2_distr_Tdata.ave(), agm2_distr_Tdata.err(), agm2_distr_Tdata_lower.ave(), agm2_distr_Tdata_lower.err()}, "../data/gm2/light/tm/agm2_Tdata_"+V_light_1.Tag[i_ens]+".dat.t", "", "#id  Tdata   ag2m agm2_err  agm2_lower agm2_lower_err");
    Print_To_File({}, {Tdata_vec, agm2_distr_OS_Tdata.ave(), agm2_distr_OS_Tdata.err(), agm2_distr_OS_Tdata_lower.ave(), agm2_distr_OS_Tdata_lower.err()}, "../data/gm2/light/OS/agm2_Tdata_"+V_light_1.Tag[i_ens]+".dat.t", "", "#id  Tdata   ag2m agm2_err agm2_lower agm2_lower_err");
  
  
 

    //#######################  INTERMEDIATE AND SHORT-DISTANCE WINDOW ###################################

    
    //############################   TWISTED MASS ######################################################

    distr_t agm2_W(UseJack, UseJack?Njacks:Nboots);  //constructor sets agm2_W to zero by default
    distr_t agm2_SD(UseJack, UseJack?Njacks:Nboots);  //constructor sets agm2_SD to zero by default
    distr_t agm2_W_ELM(UseJack, UseJack?Njacks:Nboots); //constructor sets agm2_W_ELM to zero by default
    distr_t agm2_SD_ELM(UseJack, UseJack?Njacks:Nboots);  //constructor sets agm2_SD_ELM to zero by default

    distr_t agm2_W_m_small(UseJack, UseJack?Njacks:Nboots);
    distr_t agm2_W_m_big(UseJack, UseJack?Njacks:Nboots);
    distr_t agm2_SD_m_small(UseJack, UseJack?Njacks:Nboots);
    distr_t agm2_SD_m_big(UseJack, UseJack?Njacks:Nboots);
    distr_t agm2_SD_sea_dep(UseJack, UseJack?Njacks:Nboots);
    distr_t agm2_W_sea_dep(UseJack, UseJack?Njacks:Nboots);
   
    //test behavior with lattice spacing
    distr_t agm2_W_2a(UseJack, UseJack?Njacks:Nboots);
    distr_t_list agm2_W_func_a(UseJack);
    double a_min= 0.03*fm_to_inv_Gev;
    double a_max= 0.2*fm_to_inv_Gev;
    int Nsteps_a= (int)((a_max-a_min)/(0.005*fm_to_inv_Gev));
    distr_t agm2_W_rel_err(UseJack, UseJack?Njacks:Nboots);
    distr_t agm2_SD_rel_err(UseJack, UseJack?Njacks:Nboots);
    distr_t agm2_tot(UseJack, UseJack?Njacks:Nboots);
    distr_t agm2_tot_rel_err(UseJack, UseJack?Njacks:Nboots);


    //#################################################################################################

    auto pow_1_5 = [](double x) { return pow(x,1.0/5.0);};
    distr_t X_pi = distr_t::f_of_distr(pow_1_5,Mpi*Mpi*Mpi*Mpi*fp);


  
    distr_t X = X_pi/X_pi_phys;  //Mpi/( (Mpi/a_distr).ave()); //variable to be used in ELT

    distr_t_list Ker = distr_t_list::f_of_distr(K, a_distr , Corr.Nt);
    distr_t_list Ker_2a= distr_t_list::f_of_distr(K, a_distr_2a, Corr.Nt);
    distr_t_list Ker_ELM = distr_t_list::f_of_distr(K, X, Corr.Nt);
    auto K_der_W_func_times_a = [&](double Mv, double t, double size) -> double { return Mv*der_kernel_K_W_win(t, Mv);};
    auto K_der_SD_func_times_a = [&](double Mv, double t, double size) -> double { return Mv*der_kernel_K_SD_win(t, Mv);};
    auto K_der_func_times_a = [&](double Mv, double t, double size) -> double { return Mv*der_kernel_K(t, Mv);};
    //auto K = [&](double Mv, double t, double size) -> double { return kernel_K(t, Mv);};
    distr_t_list Ker_der_W= distr_t_list::f_of_distr(K_der_W_func_times_a, a_distr, Corr.Nt);
    distr_t_list Ker_der_SD= distr_t_list::f_of_distr(K_der_SD_func_times_a, a_distr, Corr.Nt);
    distr_t_list Ker_der= distr_t_list::f_of_distr(K_der_func_times_a, a_distr, Corr.Nt);

    
    //define lambdas for the theta func
    auto th0 = [](double ta) ->double { return 1.0/(1.0 + exp(-2.0*(ta-t0)/Delta));};
    auto th1 = [](double ta) ->double { return 1.0/(1.0 + exp(-2.0*(ta-t1)/Delta));};

    //get epsilon-window
    Get_amu_W_eps(amu_W_eps_list_tm, Za*Za*V_light_distr, a_distr);
    Get_amu_W_eps(amu_W_eps_val_list_tm, Za*Za*(V_light_m_small_distr-V_light_m_big_distr), a_distr);
    if(V_light_1.Tag[i_ens] != "cB211b.072.96") 	Get_amu_W_eps(amu_W_eps_sea_list_tm, Za*Za*delta_corr_sea_tm, a_distr);
    
    for(int t=1; t< Corr.Nt/2; t++) {
      agm2_W = agm2_W + 4.0*w(t,Simps_ord)*pow(alpha,2)*Za*Za*V_light_distr.distr_list[t]*Ker.distr_list[t]*( distr_t::f_of_distr(th0, t*a_distr) - distr_t::f_of_distr(th1, t*a_distr));
      agm2_tot= agm2_tot + 4.0*w(t,Simps_ord)*pow(alpha,2)*Za*Za*V_light_distr.distr_list[t]*Ker.distr_list[t];
      agm2_W_rel_err = agm2_W_rel_err + 4.0*w(t,Simps_ord)*pow(alpha,2)*Za*Za*V_light_distr.distr_list[t]*Ker_der_W.distr_list[t];
      agm2_SD_rel_err = agm2_SD_rel_err + 4.0*w(t,Simps_ord)*pow(alpha,2)*Za*Za*V_light_distr.distr_list[t]*Ker_der_SD.distr_list[t];
      agm2_tot_rel_err = agm2_tot_rel_err + 4.0*w(t,Simps_ord)*pow(alpha,2)*Za*Za*V_light_distr.distr_list[t]*Ker_der.distr_list[t];
      agm2_W_2a = agm2_W_2a + 4.0*w(t,Simps_ord)*pow(alpha,2)*Za*Za*V_light_distr.distr_list[t]*Ker_2a.distr_list[t]*( distr_t::f_of_distr(th0, t*a_distr_2a) - distr_t::f_of_distr(th1, t*a_distr_2a));
      agm2_SD = agm2_SD + 4.0*w(t,Simps_ord)*pow(alpha,2)*Za*Za*(V_light_distr_tm_corr.distr_list[t])*Ker.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*a_distr));
      agm2_W_ELM = agm2_W_ELM + 4.0*w(t,Simps_ord)*pow(alpha,2)*Za*Za*V_light_distr.distr_list[t]*Ker_ELM.distr_list[t]*( distr_t::f_of_distr(th0, t*X) - distr_t::f_of_distr(th1, t*X));
      agm2_SD_ELM = agm2_SD_ELM + 4.0*w(t,Simps_ord)*pow(alpha,2)*(Za*Za*V_light_distr_tm_corr.distr_list[t])*Ker.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*a_distr));
      
      
      agm2_W_m_small = agm2_W_m_small + 4.0*w(t,Simps_ord)*pow(alpha,2)*Za*Za*V_light_m_small_distr.distr_list[t]*Ker.distr_list[t]*( distr_t::f_of_distr(th0, t*a_distr) - distr_t::f_of_distr(th1, t*a_distr));
      agm2_SD_m_small = agm2_SD_m_small + 4.0*w(t,Simps_ord)*pow(alpha,2)*Za*Za*(V_light_m_small_distr.distr_list[t])*Ker.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*a_distr));
      agm2_W_m_big = agm2_W_m_big + 4.0*w(t,Simps_ord)*pow(alpha,2)*Za*Za*V_light_m_big_distr.distr_list[t]*Ker.distr_list[t]*( distr_t::f_of_distr(th0, t*a_distr) - distr_t::f_of_distr(th1, t*a_distr));
      agm2_SD_m_big = agm2_SD_m_big + 4.0*w(t,Simps_ord)*pow(alpha,2)*Za*Za*(V_light_m_big_distr.distr_list[t])*Ker.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*a_distr));
	
      
      if(V_light_1.Tag[i_ens] != "cB211b.072.96") {
	agm2_W_sea_dep = agm2_W_sea_dep + 4.0*w(t,Simps_ord)*pow(alpha,2)*Za*Za*delta_corr_sea_tm.distr_list[t]*Ker.distr_list[t]*( distr_t::f_of_distr(th0, t*a_distr) - distr_t::f_of_distr(th1, t*a_distr));
	agm2_SD_sea_dep = agm2_SD_sea_dep + 4.0*w(t,Simps_ord)*pow(alpha,2)*Za*Za*(delta_corr_sea_tm.distr_list[t])*Ker.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*a_distr));
      }
    }

    agm2_W_rel_err= agm2_W_rel_err/agm2_W;
    agm2_SD_rel_err= agm2_SD_rel_err/agm2_SD;
    agm2_tot_rel_err = agm2_tot_rel_err/agm2_tot;

  
    //test, compute a_mu^W as a function of lattice spacing
    Vfloat a_test_list;
    for(int is=0;is<Nsteps_a;is++) {
      agm2_W_func_a.distr_list.push_back(0.0*agm2_W);
      double a_to_use= a_min + is*0.005*fm_to_inv_Gev;
      a_test_list.push_back(a_to_use);
      for(int t=1;t < Corr.Nt/2;t++) {
	agm2_W_func_a.distr_list[is] = agm2_W_func_a.distr_list[is] + 4.0*w(t,Simps_ord)*pow(alpha,2)*Za*Za*V_light_distr.distr_list[t]*kernel_K(t,a_to_use)*( th0(t*a_to_use) - th1(t*a_to_use));

      }
    }

  
    //print a_mu^W as a function of lattice spacing
    Print_To_File({},{a_test_list, agm2_W_func_a.ave(), agm2_W_func_a.err() }, "../data/gm2/light/tm/agm2_W_func_a_"+V_light_1.Tag[i_ens]+".dat.t", "", "a[GeV] a_mu^W a_mu^W(err)");

    vector<distr_t> agm2_SD_tmins( tmins.size());

    cout<<"Starting computation of a_mu^SD-tm for t0s: ";
    printV(tmins, "", 0);
    cout<<endl;

    int t_counter=0;
    double tolerance = 1e-14;
    double err;

    //load R_ratio at order alpha_s^4(1/tmin)

    vector<Vfloat> R_ratio_alpha4(tmins.size());
    vector<Vfloat> R_ratio_ergs_alpha4(tmins.size());
    Vfloat SD_at_alpha_4;
    for(int tm=0; tm<(signed)tmins.size();tm++) {

      R_ratio_alpha4[tm]= Read_From_File("../rhad_pert_results/table_scan_"+to_string(tm)+".dat", 1, 2);
      R_ratio_ergs_alpha4[tm] = Read_From_File("../rhad_pert_results/table_scan_"+to_string(tm)+".dat", 0, 2);

    }

    double pert_result_SD, pert_err_SD;

    auto integrand_SD_pt_cont = [&L_info, &a_distr, &th0](double t) {
				  return (t<1e-10)?0.0:(4.0*pow(alpha,2)*( pow(qu,2)+ pow(qd,2))*free_vector_corr_cont(3,0.0,t)*(1.0 - th0(t))*kernel_K(t,1.0));

				};
    
     gsl_function_pp<decltype(integrand_SD_pt_cont)> F_SD_pert(integrand_SD_pt_cont);
     gsl_integration_workspace * w_SD_pert = gsl_integration_workspace_alloc (10000);
     gsl_function *G_SD_pert = static_cast<gsl_function*>(&F_SD_pert);
     gsl_integration_qags(G_SD_pert, 0.0, 100.0, 0.0, 1e-8, 10000, w_SD_pert, &pert_result_SD, &pert_err_SD);
     gsl_integration_workspace_free (w_SD_pert);

     pert_result_SD_list.push_back(pert_result_SD);

     
			   
    for(auto & tm0: tmins) {

      if(tm0*fm_to_inv_Gev < a_distr.ave()) crash("tm0 < a");

      auto integrand_pt_cont = [&L_info, &a_distr, &th0, &R_ratio_alpha4, &R_ratio_ergs_alpha4, &t_counter](double t) {
				 //return (t<1e-10)?0.0:(4.0*pow(alpha,2)*( pow(qu,2)+ pow(qd,2))*free_vector_corr_cont(3, 0.0,t)*(1.0 - th0(t))*kernel_K(t,1.0));
				 double fact= 0.0005*(5.0/4.0)*(1.0/(12*M_PI*M_PI));
				 double corr=0.0;
				 if(t<1e-10) return 0.0;
				 for(int j=0;j<(signed)R_ratio_ergs_alpha4[t_counter].size();j++) {
				   double Erg= R_ratio_ergs_alpha4[t_counter][j];
				   double R_rat= R_ratio_alpha4[t_counter][j];
				   if( Erg > 0.2) corr+= fact*pow(Erg,2)*exp(-Erg*t)*R_rat;
				 }
				 return 4.0*pow(alpha,2)*corr*(1.0-th0(t))*kernel_K(t,1.0);
				 

			       };


     
      double agm2_pert_up_to_t0, agm2_pert_up_to_t0_err;


      //set numerical integration with gsl

      gsl_function_pp<decltype(integrand_pt_cont)> Fp(integrand_pt_cont);
      gsl_integration_workspace * w_t0 = gsl_integration_workspace_alloc (10000);
      gsl_function *G = static_cast<gsl_function*>(&Fp);
      gsl_integration_qags(G, 0.0, tm0*fm_to_inv_Gev, 0.0, 1e-8, 10000, w_t0, &agm2_pert_up_to_t0, &agm2_pert_up_to_t0_err);
      gsl_integration_workspace_free (w_t0);
      if(agm2_pert_up_to_t0_err/agm2_pert_up_to_t0 > 1e-5) crash("In determining continuum integral between 0 and tmin, accuracy achieved is only: "+to_string_with_precision(agm2_pert_up_to_t0_err/agm2_pert_up_to_t0,6));

      //double agm2_pert_up_to_t0 = boost::math::quadrature::gauss_kronrod<double, 15>::integrate( integrand_pt_cont, 0.0, tm0*fm_to_inv_Gev, 5,tolerance, &err);

      for(int ijack=0;ijack<Njacks;ijack++) {

	auto integrand_SD_light_interpolated = [&ijack, &Za, &Vt3_light_tm_interpol, &a_distr, &th0, &tm0, &L_info, &Corr](double t) {


						 double corr=0.0;
						 if( t < tm0*fm_to_inv_Gev/a_distr.distr[ijack]) {
						   corr= 0.0;
						 }
						 else {
						   corr = Vt3_light_tm_interpol(t).distr[ijack]/pow(t,3) ;
						 }

						 return 4.0*pow(alpha,2)*pow(Za.distr[ijack],2)*corr*kernel_K(t, a_distr.distr[ijack])*( 1.0 - distr_t::f_of_distr(th0, t*a_distr).distr[ijack]);

					       };

	gsl_function_pp<decltype(integrand_SD_light_interpolated)> F_SD(integrand_SD_light_interpolated);
	gsl_integration_workspace * w_SD = gsl_integration_workspace_alloc (10000);
	gsl_function *G_SD = static_cast<gsl_function*>(&F_SD);
	double agm2_SD_t0, agm2_SD_t0_err;
	gsl_integration_qags(G_SD, 0.0, Corr.Nt/2 -1.0, 0.0, 1e-9, 10000, w_SD, &agm2_SD_t0, &agm2_SD_t0_err);
	gsl_integration_workspace_free (w_SD);

	//double agm2_SD_t0 = boost::math::quadrature::gauss_kronrod<double, 15>::integrate( integrand_SD_light_interpolated, 0.0, Corr.Nt/2 -1.0, 5,tolerance, &err);

	

	agm2_SD_tmins[t_counter].distr.push_back(agm2_SD_t0+ agm2_pert_up_to_t0);
	
      }

      SD_at_alpha_4.push_back( agm2_pert_up_to_t0);
      t_counter++;
      
    }

    //print interpolated tm-corr
    //###############################################################################
    double dx =0.1;
    int Nsteps= (1.0/dx)*(Corr.Nt/2 -1);
    Vfloat tt_interpolation;
    distr_t_list interpolated_corr;
    for(int istep=0; istep< Nsteps;istep++) {
      double t= 1.0+istep*dx;
      tt_interpolation.push_back(t);
      interpolated_corr.distr_list.push_back( Vt3_light_tm_interpol(t)/pow(t,3));
    }
    
    Print_To_File({}, {tt_interpolation, interpolated_corr.ave(), interpolated_corr.err(), (Za*Za*interpolated_corr).ave(), (Za*Za*interpolated_corr).err()}, "../data/gm2/light/tm/interpolated_corr_"+V_light_1.Tag[i_ens]+".dat" , "", "#t corr_bare corr_renorm");
    Print_To_File({}, {V_light_distr_tm_corr.ave(), V_light_distr_tm_corr.err(), (Za*Za*V_light_distr_tm_corr).ave(), (Za*Za*V_light_distr_tm_corr).err()},  "../data/gm2/light/tm/corr_sum_pert_sub_"+V_light_1.Tag[i_ens]+".dat" , "", "#t corr_bare corr_renorm");
    //###############################################################################


    
  
    //push_back the result

    agm2_light_W.distr_list.push_back(agm2_W);
    agm2_light_SD.distr_list.push_back(agm2_SD);
    agm2_light_W_ELM.distr_list.push_back(agm2_W_ELM);
    agm2_light_SD_ELM.distr_list.push_back(agm2_SD_ELM);

    agm2_light_W_2a.distr_list.push_back(agm2_W_2a);
    agm2_light_W_der_err_scale_setting.distr_list.push_back(agm2_W_rel_err);
    agm2_light_SD_der_err_scale_setting.distr_list.push_back(agm2_SD_rel_err);
    agm2_light_tot_der_err_scale_setting.distr_list.push_back(agm2_tot_rel_err);

  
    //push_back the result for tmins
    for(int tmin_id=0; tmin_id< (signed)tmins.size();tmin_id++) agm2_light_SD_tmins_distr_list[tmin_id].distr_list.push_back( agm2_SD_tmins[tmin_id]);

    //####################################################################################################



    // ######################################## OSTERWALDER-SEILER #######################################


    distr_t agm2_W_OS(UseJack, UseJack?Njacks:Nboots);  //constructor sets agm2_W_OS to zero by default
    distr_t agm2_SD_OS(UseJack, UseJack?Njacks:Nboots);  //constructor sets agm2_SD_OS to zero by default
    distr_t agm2_W_ELM_OS(UseJack, UseJack?Njacks:Nboots); //constructor sets agm2_W_ELM_OS to zero by default
    distr_t agm2_SD_ELM_OS(UseJack, UseJack?Njacks:Nboots);  //constructor sets agm2_SD_ELM_OS to zero by default

    distr_t agm2_W_OS_m_small(UseJack, UseJack?Njacks:Nboots);
    distr_t agm2_W_OS_m_big(UseJack, UseJack?Njacks:Nboots);
    distr_t agm2_W_OS_sea_dep(UseJack, UseJack?Njacks:Nboots);
    distr_t agm2_SD_OS_m_small(UseJack, UseJack?Njacks:Nboots);
    distr_t agm2_SD_OS_m_big(UseJack, UseJack?Njacks:Nboots);
    distr_t agm2_SD_OS_sea_dep(UseJack, UseJack?Njacks:Nboots);


    //get epsilon-window
    Get_amu_W_eps(amu_W_eps_list_OS, Zv*Zv*V_light_OS_distr, a_distr);
    Get_amu_W_eps(amu_W_eps_val_list_OS, Zv*Zv*(V_light_OS_m_small_distr-V_light_OS_m_big_distr), a_distr);
    if(V_light_1.Tag[i_ens] != "cB211b.072.96")  	Get_amu_W_eps(amu_W_eps_sea_list_OS, Zv*Zv*delta_corr_sea_OS, a_distr);
    //#################################################################################################

    for(int t=1; t< Corr.Nt/2; t++) {
      agm2_W_OS = agm2_W_OS + 4.0*w(t,Simps_ord)*pow(alpha,2)*Zv*Zv*V_light_OS_distr.distr_list[t]*Ker.distr_list[t]*( distr_t::f_of_distr(th0, t*a_distr) - distr_t::f_of_distr(th1, t*a_distr));
      agm2_SD_OS = agm2_SD_OS + 4.0*w(t,Simps_ord)*pow(alpha,2)*Zv*Zv*(V_light_distr_OS_corr.distr_list[t])*Ker.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*a_distr));
      agm2_W_ELM_OS = agm2_W_ELM_OS + 4.0*w(t,Simps_ord)*pow(alpha,2)*Zv*Zv*V_light_OS_distr.distr_list[t]*Ker_ELM.distr_list[t]*( distr_t::f_of_distr(th0, t*X) - distr_t::f_of_distr(th1, t*X));
      agm2_SD_ELM_OS = agm2_SD_ELM_OS + 4.0*w(t,Simps_ord)*pow(alpha,2)*(Zv*Zv*V_light_distr_OS_corr.distr_list[t])*Ker.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*a_distr));
      

       
      agm2_W_OS_m_small = agm2_W_OS_m_small + 4.0*w(t,Simps_ord)*pow(alpha,2)*Zv*Zv*V_light_OS_m_small_distr.distr_list[t]*Ker.distr_list[t]*( distr_t::f_of_distr(th0, t*a_distr) - distr_t::f_of_distr(th1, t*a_distr));
      agm2_SD_OS_m_small = agm2_SD_OS_m_small + 4.0*w(t,Simps_ord)*pow(alpha,2)*Zv*Zv*(V_light_OS_m_small_distr.distr_list[t])*Ker.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*a_distr));
      agm2_W_OS_m_big = agm2_W_OS_m_big + 4.0*w(t,Simps_ord)*pow(alpha,2)*Zv*Zv*V_light_OS_m_big_distr.distr_list[t]*Ker.distr_list[t]*( distr_t::f_of_distr(th0, t*a_distr) - distr_t::f_of_distr(th1, t*a_distr));
      agm2_SD_OS_m_big = agm2_SD_OS_m_big + 4.0*w(t,Simps_ord)*pow(alpha,2)*Zv*Zv*(V_light_OS_m_big_distr.distr_list[t])*Ker.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*a_distr));
      
      if(V_light_1.Tag[i_ens] != "cB211b.072.96") {
	agm2_SD_OS_sea_dep = agm2_SD_OS_sea_dep + 4.0*w(t,Simps_ord)*pow(alpha,2)*Zv*Zv*(delta_corr_sea_OS.distr_list[t])*Ker.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*a_distr));
	agm2_W_OS_sea_dep = agm2_W_OS_sea_dep +  4.0*w(t,Simps_ord)*pow(alpha,2)*Zv*Zv*delta_corr_sea_OS.distr_list[t]*Ker.distr_list[t]*( distr_t::f_of_distr(th0, t*a_distr) - distr_t::f_of_distr(th1, t*a_distr));
       }
    }


    //#############if ensemble cB64 print valence mass correction to the window
    if(i_ens==0) {ofstream Print_mass_dep_SD_val_tm("../data/gm2/light/mass_dep/tm/win_SD_val.dat"); ofstream Print_mass_dep_W_val_tm("../data/gm2/light/mass_dep/tm/win_W_val.dat"); Print_mass_dep_SD_val_tm.close(); Print_mass_dep_W_val_tm.close(); ofstream Print_mass_dep_SD_val_OS("../data/gm2/light/mass_dep/OS/win_SD_val.dat"); ofstream Print_mass_dep_W_val_OS("../data/gm2/light/mass_dep/OS/win_W_val.dat"); Print_mass_dep_SD_val_OS.close(); Print_mass_dep_W_val_OS.close();}
 
   
    //valence tm
    ofstream Print_mass_dep_SD_val_tm, Print_mass_dep_W_val_tm;
      
    Print_mass_dep_SD_val_tm.open("../data/gm2/light/mass_dep/tm/win_SD_val.dat", ofstream::app);
    Print_mass_dep_W_val_tm.open("../data/gm2/light/mass_dep/tm/win_W_val.dat", ofstream::app);
    Print_mass_dep_SD_val_tm<<V_light_1.Tag[i_ens]<<" "<<(a_distr/fm_to_inv_Gev).ave()<<" "<<1e10*(agm2_SD_m_small-agm2_SD_m_big).ave()<<"  "<<1e10*(agm2_SD_m_small-agm2_SD_m_big).err()<<" "<<1e10*(agm2_SD_m_small-agm2_SD_m_big+agm2_SD).ave()<<" "<<1e10*(agm2_SD_m_small-agm2_SD_m_big+agm2_SD).err()<<" "<<1e10*agm2_SD_m_small.ave()<<" "<<1e10*agm2_SD_m_small.err()<<" "<<1e10*agm2_SD_m_big.ave()<<" "<<1e10*agm2_SD_m_big.err()<<endl;
    Print_mass_dep_W_val_tm<<V_light_1.Tag[i_ens]<<" "<<(a_distr/fm_to_inv_Gev).ave()<<" "<<1e10*(agm2_W_m_small-agm2_W_m_big).ave()<<"  "<<1e10*(agm2_W_m_small-agm2_W_m_big).err()<<" "<<1e10*(agm2_W_m_small-agm2_W_m_big+ agm2_W).ave()<<" "<<1e10*(agm2_W_m_small-agm2_W_m_big+agm2_W).err()<<" "<<1e10*agm2_W_m_small.ave()<<" "<<1e10*agm2_W_m_small.err()<<" "<<1e10*agm2_W_m_big.ave()<<" "<<1e10*agm2_W_m_big.err()<<endl;
    Print_mass_dep_SD_val_tm.close();
    Print_mass_dep_W_val_tm.close();
    //push_back
    amu_W_mass_corr_tm_list.distr_list.push_back( agm2_W_m_small-agm2_W_m_big);
    amu_SD_mass_corr_tm_list.distr_list.push_back( agm2_SD_m_small- agm2_SD_m_big);

    //valence OS
    ofstream Print_mass_dep_SD_val_OS, Print_mass_dep_W_val_OS;
    Print_mass_dep_SD_val_OS.open("../data/gm2/light/mass_dep/OS/win_SD_val.dat", ofstream::app);
    Print_mass_dep_W_val_OS.open("../data/gm2/light/mass_dep/OS/win_W_val.dat", ofstream::app); 

    Print_mass_dep_SD_val_OS<<V_light_1.Tag[i_ens]<<" "<<(a_distr/fm_to_inv_Gev).ave()<<" "<<1e10*(agm2_SD_OS_m_small-agm2_SD_OS_m_big).ave()<<"  "<<1e10*(agm2_SD_OS_m_small-agm2_SD_OS_m_big).err()<<" "<<1e10*(agm2_SD_OS_m_small-agm2_SD_OS_m_big+agm2_SD_OS).ave()<<" "<<1e10*(agm2_SD_OS_m_small-agm2_SD_OS_m_big+agm2_SD_OS).err()<<" "<<1e10*agm2_SD_OS_m_small.ave()<<" "<<1e10*agm2_SD_OS_m_small.err()<<" "<<1e10*agm2_SD_OS_m_big.ave()<<" "<<1e10*agm2_SD_OS_m_big.err()<<endl;
    Print_mass_dep_W_val_OS<<V_light_1.Tag[i_ens]<<" "<<(a_distr/fm_to_inv_Gev).ave()<<" "<<1e10*(agm2_W_OS_m_small-agm2_W_OS_m_big).ave()<<"  "<<1e10*(agm2_W_OS_m_small-agm2_W_OS_m_big).err()<<" "<<1e10*(agm2_W_OS_m_small-agm2_W_OS_m_big+ agm2_W_OS).ave()<<" "<<1e10*(agm2_W_OS_m_small-agm2_W_OS_m_big+agm2_W_OS).err()<<" "<<1e10*agm2_W_OS_m_small.ave()<<" "<<1e10*agm2_W_OS_m_small.err()<<" "<<1e10*agm2_W_OS_m_big.ave()<<" "<<1e10*agm2_W_OS_m_big.err()<<endl;
    Print_mass_dep_SD_val_OS.close();
    Print_mass_dep_W_val_OS.close();
    //push_back
    amu_W_mass_corr_OS_list.distr_list.push_back( agm2_W_OS_m_small-agm2_W_OS_m_big);
    amu_SD_mass_corr_OS_list.distr_list.push_back( agm2_SD_OS_m_small - agm2_SD_OS_m_big);

    
    amu_W_val_Ens_tag.push_back( V_light_1.Tag[i_ens]);

  
    

    //print sea quark mistuning contributions to the window
    if(i_ens==0) {
       ofstream Print_mass_dep_SD_sea_tm("../data/gm2/light/mass_dep/tm/win_SD_sea.dat");
       ofstream Print_mass_dep_W_sea_tm("../data/gm2/light/mass_dep/tm/win_W_sea.dat");
       ofstream Print_mass_dep_SD_sea_OS("../data/gm2/light/mass_dep/OS/win_SD_sea.dat");
       ofstream Print_mass_dep_W_sea_OS("../data/gm2/light/mass_dep/OS/win_W_sea.dat");
       Print_mass_dep_SD_sea_tm.close();
       Print_mass_dep_W_sea_tm.close();
       Print_mass_dep_SD_sea_OS.close();
       Print_mass_dep_W_sea_OS.close();

    }

    
    if(V_light_1.Tag[i_ens] != "cB211b.072.96") {


      //sea tm
      ofstream Print_mass_dep_SD_sea_tm, Print_mass_dep_W_sea_tm;
      Print_mass_dep_SD_sea_tm.open("../data/gm2/light/mass_dep/tm/win_SD_sea.dat", ofstream::app);
      Print_mass_dep_W_sea_tm.open("../data/gm2/light/mass_dep/tm/win_W_sea.dat", ofstream::app);
      Print_mass_dep_SD_sea_tm<<V_light_1.Tag[i_ens]<<" "<<(a_distr/fm_to_inv_Gev).ave()<<" "<<1e10*(agm2_SD_sea_dep).ave()<<" "<<1e10*(agm2_SD_sea_dep).err()<<" "<<1e10*(agm2_SD_sea_dep+agm2_SD).ave()<<" "<<1e10*(agm2_SD_sea_dep+agm2_SD).err()<<endl;
      Print_mass_dep_W_sea_tm<<V_light_1.Tag[i_ens]<<" "<<(a_distr/fm_to_inv_Gev).ave()<<" "<<1e10*(agm2_W_sea_dep).ave()<<" "<<1e10*(agm2_W_sea_dep).err()<<" "<<1e10*(agm2_W_sea_dep+agm2_W).ave()<<" "<<1e10*(agm2_W_sea_dep+agm2_W).err()<<endl;
      Print_mass_dep_SD_sea_tm.close();
      Print_mass_dep_W_sea_tm.close();

      amu_W_mass_corr_sea_tm_list.distr_list.push_back(agm2_W_sea_dep);
      amu_SD_mass_corr_sea_tm_list.distr_list.push_back(agm2_SD_sea_dep);
   

      //sea OS
      ofstream Print_mass_dep_SD_sea_OS, Print_mass_dep_W_sea_OS;
      Print_mass_dep_SD_sea_OS.open("../data/gm2/light/mass_dep/OS/win_SD_sea.dat", ofstream::app);
      Print_mass_dep_W_sea_OS.open("../data/gm2/light/mass_dep/OS/win_W_sea.dat", ofstream::app); 
      Print_mass_dep_SD_sea_OS<<V_light_1.Tag[i_ens]<<" "<<(a_distr/fm_to_inv_Gev).ave()<<" "<<1e10*(agm2_SD_OS_sea_dep).ave()<<" "<<1e10*(agm2_SD_OS_sea_dep).err()<<" "<<1e10*(agm2_SD_OS_sea_dep+ agm2_SD_OS).ave()<<" "<<1e10*(agm2_SD_OS_sea_dep+agm2_SD_OS).err()<<endl;
      Print_mass_dep_W_sea_OS<<V_light_1.Tag[i_ens]<<" "<<(a_distr/fm_to_inv_Gev).ave()<<" "<<1e10*(agm2_W_OS_sea_dep).ave()<<" "<<1e10*(agm2_W_OS_sea_dep).err()<<" "<<1e10*(agm2_W_OS_sea_dep+ agm2_W_OS).ave()<<" "<<1e10*(agm2_W_OS_sea_dep+agm2_W_OS).err()<<endl;
      Print_mass_dep_SD_sea_OS.close();
      Print_mass_dep_W_sea_OS.close();

      amu_W_mass_corr_sea_OS_list.distr_list.push_back(agm2_W_OS_sea_dep);
      amu_SD_mass_corr_sea_OS_list.distr_list.push_back(agm2_SD_OS_sea_dep);

      amu_W_sea_Ens_tag.push_back( V_light_1.Tag[i_ens]);
      
    }


    //##################### COMPUTE COVARIANCE MATRIX FOR TM OS CORRELATORS ############################


    double cov_W = agm2_W%agm2_W_OS/(agm2_W_OS.err()*agm2_W.err());
    double cov_SD = agm2_SD%agm2_SD_OS/(agm2_SD.err()*agm2_SD_OS.err());


    Print_To_File({}, {Vfloat({cov_SD}), Vfloat({cov_W})}, "../data/gm2/light/corr_matrix_"+V_light_1.Tag[i_ens]+".dat", "", "#SD W");

    //##################################################################################################


    vector<distr_t> agm2_SD_OS_tmins( tmins.size());

    cout<<"Starting computation of a_mu^SD-OS for t0s: ";
    printV(tmins, "", 0);
    cout<<endl;

    t_counter=0;
    tolerance = 1e-14;
 
			   
    for(auto & tm0: tmins) {

      if(tm0*fm_to_inv_Gev < a_distr.ave()) crash("tm0 < a");

      /*

      auto integrand_pt_cont = [&L_info, &a_distr, &th0](double t) { return (t<1e-10)?0.0:4.0*pow(alpha,2)*( pow(qu,2)+ pow(qd,2)  )*free_vector_corr_cont(3, L_info.ml/a_distr.ave(),t)*(1.0 - th0(t))*kernel_K(t,1.0);   };

      double agm2_pert_up_to_t0, agm2_pert_up_to_t0_err;

      //set numerical integration with gsl

      gsl_function_pp<decltype(integrand_pt_cont)> Fp(integrand_pt_cont);
      gsl_integration_workspace * w_t0 = gsl_integration_workspace_alloc (10000);
      gsl_function *G = static_cast<gsl_function*>(&Fp);
      gsl_integration_qags(G, 0.0, tm0*fm_to_inv_Gev, 0.0, 1e-9, 10000, w_t0, &agm2_pert_up_to_t0, &agm2_pert_up_to_t0_err);
      gsl_integration_workspace_free (w_t0);

      //double agm2_pert_up_to_t0 = boost::math::quadrature::gauss_kronrod<double, 15>::integrate( integrand_pt_cont, 0.0, tm0*fm_to_inv_Gev, 5,tolerance, &err);

      */

      double agm2_pert_up_to_t0 = SD_at_alpha_4[t_counter];

      for(int ijack=0;ijack<Njacks;ijack++) {

	auto integrand_SD_light_interpolated = [&ijack, &Zv, &Vt3_light_OS_interpol, &a_distr, &tm0, &th0, &L_info, &Corr](double t) {
																    
						 double corr=0.0;

					
						 if( t < tm0*fm_to_inv_Gev/a_distr.distr[ijack]) {
						   corr= 0.0;
						 }
						 else {
						   corr = Vt3_light_OS_interpol(t).distr[ijack]/pow(t,3);
						 }
						 return 4.0*pow(alpha,2)*pow(Zv.distr[ijack],2)*corr*kernel_K(t, a_distr.distr[ijack])*( 1.0 - distr_t::f_of_distr(th0, t*a_distr).distr[ijack]);
					       };

	gsl_function_pp<decltype(integrand_SD_light_interpolated)> F_SD(integrand_SD_light_interpolated);
	gsl_integration_workspace * w_SD = gsl_integration_workspace_alloc (10000);
	gsl_function *G_SD = static_cast<gsl_function*>(&F_SD);
	double agm2_SD_t0, agm2_SD_t0_err;
	gsl_integration_qags(G_SD, 0.0, Corr.Nt/2 -1.0, 0.0, 1e-8, 10000, w_SD, &agm2_SD_t0, &agm2_SD_t0_err);
	gsl_integration_workspace_free (w_SD);


	//agm2_SD_t0 = boost::math::quadrature::gauss_kronrod<double, 15>::integrate( integrand_SD_light_interpolated, 0.0, Corr.Nt/2 -1.0, 5,tolerance, &err);

	agm2_SD_OS_tmins[t_counter].distr.push_back( agm2_SD_t0+ agm2_pert_up_to_t0);
      }

      t_counter++;
      
    }


    //print interpolated OS-corr
    //###############################################################################
    dx =0.1;
    Nsteps= (1.0/dx)*(Corr.Nt/2 -1);
    Vfloat tt_interpolation_OS;
    distr_t_list interpolated_corr_OS;
    for(int istep=0; istep< Nsteps;istep++) {
      double t= 1.0+istep*dx;
      tt_interpolation_OS.push_back(t);
      interpolated_corr_OS.distr_list.push_back(Vt3_light_OS_interpol(t)/pow(t,3));
    }
    Print_To_File({}, {tt_interpolation_OS, interpolated_corr_OS.ave(), interpolated_corr_OS.err(), (Zv*Zv*interpolated_corr_OS).ave(), (Zv*Zv*interpolated_corr_OS).err()}, "../data/gm2/light/OS/interpolated_corr_"+V_light_1.Tag[i_ens]+".dat" , "", "#t corr_bare corr_renorm");
    Print_To_File({}, {V_light_distr_OS_corr.ave(), V_light_distr_OS_corr.err(), (Zv*Zv*V_light_distr_OS_corr).ave(), (Zv*Zv*V_light_distr_OS_corr).err()},  "../data/gm2/light/OS/corr_sum_pert_sub_"+V_light_1.Tag[i_ens]+".dat" , "", "#t corr_bare corr_renorm");
    //###############################################################################
  
    //push_back the result

    agm2_light_W_OS.distr_list.push_back(agm2_W_OS);
    agm2_light_SD_OS.distr_list.push_back(agm2_SD_OS);
    agm2_light_W_ELM_OS.distr_list.push_back(agm2_W_ELM_OS);
    agm2_light_SD_ELM_OS.distr_list.push_back(agm2_SD_ELM_OS);

    //push_back the result for tmins
    for(int tmin_id=0; tmin_id< (signed)tmins.size();tmin_id++) agm2_light_SD_OS_tmins_distr_list[tmin_id].distr_list.push_back( agm2_SD_OS_tmins[tmin_id]);

    //####################################################################################################




    //####################################### DISCO LIGHT ###############################################

    //INTERMEDIATE AND SHORT-DISTANCE
    //#################################################################################
    distr_t agm2_disco_W(UseJack, UseJack?Njacks:Nboots); //constructor sets it to 0
    distr_t agm2_disco_W_ELM(UseJack, UseJack?Njacks:Nboots); //constructor sets it to 0
    distr_t agm2_disco_SD(UseJack, UseJack?Njacks:Nboots); //constructor sets it to 0
    distr_t agm2_disco_SD_ELM(UseJack, UseJack?Njacks:Nboots); //constructor sets it to 0

    //improved
    distr_t agm2_disco_impr_W(UseJack, UseJack?Njacks:Nboots); //constructor sets it to 0
    distr_t agm2_disco_impr_W_ELM(UseJack, UseJack?Njacks:Nboots); //constructor sets it to 0
    distr_t agm2_disco_impr_SD(UseJack, UseJack?Njacks:Nboots); //constructor sets it to 0
    distr_t agm2_disco_impr_SD_ELM(UseJack, UseJack?Njacks:Nboots); //constructor sets it to 0

    //improved D
    distr_t agm2_disco_impr_D_W(UseJack, UseJack?Njacks:Nboots); //constructor sets it to 0
    distr_t agm2_disco_impr_D_W_ELM(UseJack, UseJack?Njacks:Nboots); //constructor sets it to 0
    distr_t agm2_disco_impr_D_SD(UseJack, UseJack?Njacks:Nboots); //constructor sets it to 0
    distr_t agm2_disco_impr_D_SD_ELM(UseJack, UseJack?Njacks:Nboots); //constructor sets it to 0

    //improved DD
    distr_t agm2_disco_impr_DD_W(UseJack, UseJack?Njacks:Nboots); //constructor sets it to 0
    distr_t agm2_disco_impr_DD_W_ELM(UseJack, UseJack?Njacks:Nboots); //constructor sets it to 0
    distr_t agm2_disco_impr_DD_SD(UseJack, UseJack?Njacks:Nboots); //constructor sets it to 0
    distr_t agm2_disco_impr_DD_SD_ELM(UseJack, UseJack?Njacks:Nboots); //constructor sets it to 0
    //#####################################################################################

    
    //TOTAL
    distr_t_list agm2_disco_full_Tdata_low; 
    distr_t_list agm2_disco_full_ELM_Tdata_low;
    distr_t_list agm2_disco_full_Tdata_upp;
    distr_t_list agm2_disco_full_ELM_Tdata_upp;
    //improved
    distr_t_list agm2_disco_impr_full_Tdata;
    distr_t_list agm2_disco_impr_full_ELM_Tdata;
    //improved D
    distr_t_list agm2_disco_impr_D_full_Tdata;
    distr_t_list agm2_disco_impr_D_full_ELM_Tdata;
    //improved DD
    distr_t_list agm2_disco_impr_DD_full_Tdata;
    distr_t_list agm2_disco_impr_DD_full_ELM_Tdata;
    

    
    int Tdata_min_disco=1;
    int Tdata_fit_disco;
    Vfloat Tdata_vec_disco;


    //define two-pion-energy state with lower k
    //define phi_c
    //####################################################################################################
    distr_t E_2pi;
    for(int ijack=0;ijack<Njacks;ijack++) E_2pi.distr.push_back( 2.0*sqrt( Mpi.distr[ijack]*Mpi.distr[ijack] + pow(2.0*M_PI/L_info.L,2)));
    auto phi_c = [&E_2pi](double t1, double t2, double hT) -> distr_t {
		   distr_t ret(E_2pi.UseJack);
		   for(int i=0;i<Njacks;i++) ret.distr.push_back( (cosh( E_2pi.distr[i]*(t1-hT)) +1.0)/(cosh( E_2pi.distr[i]*(t2-hT)) +1.0));
		   return ret;
		 };

    //####################################################################################################

    if(Include_light_disco && Found_disco_ens) {


      if(V_light_1.Tag[i_ens].substr(1,1)=="B") Tdata_fit_disco = 28;
      else if(V_light_1.Tag[i_ens].substr(1,1)=="C") Tdata_fit_disco = 32;
      else if(V_light_1.Tag[i_ens].substr(1,1)=="D") Tdata_fit_disco = 41;
      else crash("Tdata_fit_disco for disco_light cannot be set, Ensemble: "+V_light_1.Tag[i_ens]+" not recognized");
      bool Find_Tdata_fit_disco=false;
    
    
      for(int t=1; t< Corr.Nt/2; t++) {
	agm2_disco_W = agm2_disco_W + 4.0*w(t,Simps_ord)*pow(alpha,2)*Zv*Zv*disco_distr.distr_list[t]*Ker.distr_list[t]*( distr_t::f_of_distr(th0, t*a_distr) - distr_t::f_of_distr(th1, t*a_distr));
	agm2_disco_SD = agm2_disco_SD + 4.0*w(t,Simps_ord)*pow(alpha,2)*Zv*Zv*disco_distr.distr_list[t]*Ker.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*a_distr));
	agm2_disco_W_ELM= agm2_disco_W_ELM + 4.0*w(t,Simps_ord)*pow(alpha,2)*Zv*Zv*disco_distr.distr_list[t]*Ker_ELM.distr_list[t]*( distr_t::f_of_distr(th0, t*X) - distr_t::f_of_distr(th1, t*X));
	agm2_disco_SD_ELM= agm2_disco_SD_ELM + 4.0*w(t,Simps_ord)*pow(alpha,2)*Zv*Zv*disco_distr.distr_list[t]*Ker_ELM.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*X));
      }
    

      for(int Tdata=Tdata_min_disco; Tdata< Corr.Nt/2; Tdata++) {

	distr_t agm2_disco_full_low(UseJack, UseJack?Njacks:Nboots); //constructor sets it to 0
	distr_t agm2_disco_full_ELM_low(UseJack, UseJack?Njacks:Nboots); //constructor sets it to 0
	distr_t agm2_disco_full_upp(UseJack, UseJack?Njacks:Nboots); //constructor sets it to 0
	distr_t agm2_disco_full_ELM_upp(UseJack, UseJack?Njacks:Nboots); //constructor sets it to 0
	Tdata_vec_disco.push_back(Tdata);
      
	for(int t=1; t<= Tdata; t++) {
	  agm2_disco_full_low = agm2_disco_full_low + 4.0*w(t,Simps_ord)*pow(alpha,2)*Zv*Zv*disco_distr.distr_list[t]*Ker.distr_list[t];
	  agm2_disco_full_ELM_low = agm2_disco_full_ELM_low + 4.0*w(t,Simps_ord)*pow(alpha,2)*Zv*Zv*disco_distr.distr_list[t]*Ker_ELM.distr_list[t];
	  agm2_disco_full_upp = agm2_disco_full_upp + 4.0*w(t,Simps_ord)*pow(alpha,2)*Zv*Zv*disco_distr.distr_list[t]*Ker.distr_list[t];
	  agm2_disco_full_ELM_upp = agm2_disco_full_ELM_upp + 4.0*w(t,Simps_ord)*pow(alpha,2)*Zv*Zv*disco_distr.distr_list[t]*Ker_ELM.distr_list[t];
	}

	for(int t=Tdata+1;t<Corr.Nt; t++) {

	  agm2_disco_full_upp = agm2_disco_full_upp -4.0*w(t,Simps_ord)*pow(alpha,2)*Zv*Zv*V_light_OS_distr.distr_list[Tdata]*(1.0/10.0)*Ker.distr_list[t];
	  agm2_disco_full_ELM_upp = agm2_disco_full_ELM_upp -4.0*w(t,Simps_ord)*pow(alpha,2)*Zv*Zv*V_light_OS_distr.distr_list[Tdata]*(1.0/10.0)*Ker_ELM.distr_list[t];


	}

	agm2_disco_full_Tdata_low.distr_list.push_back(agm2_disco_full_low);
	agm2_disco_full_ELM_Tdata_low.distr_list.push_back(agm2_disco_full_ELM_low);
	agm2_disco_full_Tdata_upp.distr_list.push_back( agm2_disco_full_upp);
	agm2_disco_full_ELM_Tdata_upp.distr_list.push_back(agm2_disco_full_ELM_upp);

	if(Tdata==Tdata_fit_disco) { // push_back the result

	  agm2_disco_light.distr_list.push_back(agm2_disco_full_low);
	  agm2_disco_light_ELM.distr_list.push_back(agm2_disco_full_ELM_low);
	  Find_Tdata_fit_disco = true;
	}

      }

      if(!Find_Tdata_fit_disco) crash("Tdata_fit_disco for Ensemble: "+V_light_1.Tag[i_ens]+" is not in ["+to_string(Tdata_min_disco)+","+to_string(Corr.Nt/2)+"]");

    
      //push_back the result

      agm2_disco_light_W.distr_list.push_back(agm2_disco_W);
      agm2_disco_light_SD.distr_list.push_back(agm2_disco_SD);
      agm2_disco_light_W_ELM.distr_list.push_back(agm2_disco_W_ELM);
      agm2_disco_light_SD_ELM.distr_list.push_back(agm2_disco_SD_ELM);
			      
			      
      //Print Tdata
      Print_To_File({}, {Tdata_vec_disco, agm2_disco_full_Tdata_low.ave(), agm2_disco_full_Tdata_low.err(), agm2_disco_full_Tdata_upp.ave(), agm2_disco_full_Tdata_upp.err()}, "../data/gm2/light/disco/agm2_Tdata_"+V_light_1.Tag[i_ens]+".dat.t", "", "#id  Tdata   ag2m_lower   agm2_upper");
      Print_To_File({}, {Tdata_vec_disco, agm2_disco_full_ELM_Tdata_low.ave(), agm2_disco_full_ELM_Tdata_low.err(), agm2_disco_full_ELM_Tdata_upp.ave(), agm2_disco_full_ELM_Tdata_upp.err()}, "../data/gm2/light/disco/agm2_Tdata_ELM_"+V_light_1.Tag[i_ens]+".dat.t", "", "#id  Tdata   ag2m_lower  agm2_upper");
        
    }


    //improved
    if(Include_light_disco && Found_disco_impr_ens) {
      for(int t=1; t< Corr.Nt/2; t++) {
	
	agm2_disco_impr_W = agm2_disco_impr_W + 4.0*w(t,Simps_ord)*pow(alpha,2)*Zv*Zv*disco_impr_distr.distr_list[t]*Ker.distr_list[t]*( distr_t::f_of_distr(th0, t*a_distr) - distr_t::f_of_distr(th1, t*a_distr));
	agm2_disco_impr_SD = agm2_disco_impr_SD + 4.0*w(t,Simps_ord)*pow(alpha,2)*Zv*Zv*disco_impr_distr.distr_list[t]*Ker.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*a_distr));
	agm2_disco_impr_W_ELM= agm2_disco_impr_W_ELM + 4.0*w(t,Simps_ord)*pow(alpha,2)*Zv*Zv*disco_impr_distr.distr_list[t]*Ker_ELM.distr_list[t]*( distr_t::f_of_distr(th0, t*X) - distr_t::f_of_distr(th1, t*X));
	agm2_disco_impr_SD_ELM= agm2_disco_impr_SD_ELM + 4.0*w(t,Simps_ord)*pow(alpha,2)*Zv*Zv*disco_impr_distr.distr_list[t]*Ker_ELM.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*X));
      }


       for(int Tdata=Tdata_min_disco; Tdata< Corr.Nt/2; Tdata++) {

	 distr_t agm2_disco_impr_full(UseJack, UseJack?Njacks:Nboots); //constructor sets it to 0
	 distr_t agm2_disco_impr_full_ELM(UseJack, UseJack?Njacks:Nboots); //constructor sets it to 0
		 
	 for(int t=1; t<= Tdata; t++) {
	   agm2_disco_impr_full= agm2_disco_impr_full + 4.0*w(t,Simps_ord)*pow(alpha,2)*Zv*Zv*disco_impr_distr.distr_list[t]*Ker.distr_list[t];
	   agm2_disco_impr_full_ELM = agm2_disco_impr_full_ELM + 4.0*w(t,Simps_ord)*pow(alpha,2)*Zv*Zv*disco_impr_distr.distr_list[t]*Ker_ELM.distr_list[t];
	 }


	 agm2_disco_impr_full_Tdata.distr_list.push_back(agm2_disco_impr_full);
	 agm2_disco_impr_full_ELM_Tdata.distr_list.push_back(agm2_disco_impr_full_ELM);
       }

       //push_back the result
       

       agm2_disco_impr_light_W.distr_list.push_back(agm2_disco_impr_W);
       agm2_disco_impr_light_SD.distr_list.push_back(agm2_disco_impr_SD);
       agm2_disco_impr_light_W_ELM.distr_list.push_back(agm2_disco_impr_W_ELM);
       agm2_disco_impr_light_SD_ELM.distr_list.push_back(agm2_disco_impr_SD_ELM);


       //Print Tdata
       Print_To_File({}, {Tdata_vec_disco, agm2_disco_impr_full_Tdata.ave(), agm2_disco_impr_full_Tdata.err()}, "../data/gm2/light/disco/agm2_impr_Tdata_"+V_light_1.Tag[i_ens]+".dat.t", "", "#id  Tdata   ag2m ");
       Print_To_File({}, {Tdata_vec_disco, agm2_disco_impr_full_ELM_Tdata.ave(), agm2_disco_impr_full_ELM_Tdata.err()}, "../data/gm2/light/disco/agm2_impr_Tdata_ELM_"+V_light_1.Tag[i_ens]+".dat.t", "", "#id  Tdata   ag2m");
    }

    

    //improved D
    if(Include_light_disco && Found_disco_impr_D_ens) {

      
      for(int t=1; t< Corr.Nt/2; t++) {
	agm2_disco_impr_D_W = agm2_disco_impr_D_W + 4.0*w(t,Simps_ord)*pow(alpha,2)*Zv*Zv*disco_impr_D_distr.distr_list[t]*Ker.distr_list[t]*( distr_t::f_of_distr(th0, t*a_distr) - distr_t::f_of_distr(th1, t*a_distr));
	agm2_disco_impr_D_SD = agm2_disco_impr_D_SD + 4.0*w(t,Simps_ord)*pow(alpha,2)*Zv*Zv*disco_impr_D_distr.distr_list[t]*Ker.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*a_distr));
	agm2_disco_impr_D_W_ELM= agm2_disco_impr_D_W_ELM + 4.0*w(t,Simps_ord)*pow(alpha,2)*Zv*Zv*disco_impr_D_distr.distr_list[t]*Ker_ELM.distr_list[t]*( distr_t::f_of_distr(th0, t*X) - distr_t::f_of_distr(th1, t*X));
	agm2_disco_impr_D_SD_ELM= agm2_disco_impr_D_SD_ELM + 4.0*w(t,Simps_ord)*pow(alpha,2)*Zv*Zv*disco_impr_D_distr.distr_list[t]*Ker_ELM.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*X));
      }



      for(int Tdata=Tdata_min_disco; Tdata< Corr.Nt/2; Tdata++) {

	 distr_t agm2_disco_impr_D_full(UseJack, UseJack?Njacks:Nboots); //constructor sets it to 0
	 distr_t agm2_disco_impr_D_full_ELM(UseJack, UseJack?Njacks:Nboots); //constructor sets it to 0
		 
	 for(int t=1; t<= Tdata; t++) {
	   agm2_disco_impr_D_full= agm2_disco_impr_D_full + 4.0*w(t,Simps_ord)*pow(alpha,2)*Zv*Zv*disco_impr_D_distr.distr_list[t]*Ker.distr_list[t];
	   agm2_disco_impr_D_full_ELM = agm2_disco_impr_D_full_ELM + 4.0*w(t,Simps_ord)*pow(alpha,2)*Zv*Zv*disco_impr_D_distr.distr_list[t]*Ker_ELM.distr_list[t];
	 }


	 agm2_disco_impr_D_full_Tdata.distr_list.push_back(agm2_disco_impr_D_full);
	 agm2_disco_impr_D_full_ELM_Tdata.distr_list.push_back(agm2_disco_impr_D_full_ELM);
      }

      //push_back the result

      agm2_disco_impr_lightD_W.distr_list.push_back(agm2_disco_impr_D_W);
      agm2_disco_impr_lightD_SD.distr_list.push_back(agm2_disco_impr_D_SD);
      agm2_disco_impr_lightD_W_ELM.distr_list.push_back(agm2_disco_impr_D_W_ELM);
      agm2_disco_impr_lightD_SD_ELM.distr_list.push_back(agm2_disco_impr_D_SD_ELM);


      //Print Tdata
      Print_To_File({}, {Tdata_vec_disco, agm2_disco_impr_D_full_Tdata.ave(), agm2_disco_impr_D_full_Tdata.err()}, "../data/gm2/light/disco/agm2_impr_D_Tdata_"+V_light_1.Tag[i_ens]+".dat.t", "", "#id  Tdata   ag2m ");
      Print_To_File({}, {Tdata_vec_disco, agm2_disco_impr_D_full_ELM_Tdata.ave(), agm2_disco_impr_D_full_ELM_Tdata.err()}, "../data/gm2/light/disco/agm2_impr_D_Tdata_ELM_"+V_light_1.Tag[i_ens]+".dat.t", "", "#id  Tdata   ag2m");
    }

    

    //improved DD
    if(Include_light_disco && Found_disco_impr_DD_ens) {

      
      for(int t=1; t< Corr.Nt/2; t++) {
	agm2_disco_impr_DD_W = agm2_disco_impr_DD_W + 4.0*w(t,Simps_ord)*pow(alpha,2)*Zv*Zv*disco_impr_DD_distr.distr_list[t]*Ker.distr_list[t]*( distr_t::f_of_distr(th0, t*a_distr) - distr_t::f_of_distr(th1, t*a_distr));
	agm2_disco_impr_DD_SD = agm2_disco_impr_DD_SD + 4.0*w(t,Simps_ord)*pow(alpha,2)*Zv*Zv*disco_impr_DD_distr.distr_list[t]*Ker.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*a_distr));
	agm2_disco_impr_DD_W_ELM= agm2_disco_impr_DD_W_ELM + 4.0*w(t,Simps_ord)*pow(alpha,2)*Zv*Zv*disco_impr_DD_distr.distr_list[t]*Ker_ELM.distr_list[t]*( distr_t::f_of_distr(th0, t*X) - distr_t::f_of_distr(th1, t*X));
	agm2_disco_impr_DD_SD_ELM= agm2_disco_impr_DD_SD_ELM + 4.0*w(t,Simps_ord)*pow(alpha,2)*Zv*Zv*disco_impr_DD_distr.distr_list[t]*Ker_ELM.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*X));
      }

      for(int Tdata=Tdata_min_disco; Tdata< Corr.Nt/2; Tdata++) {

	distr_t agm2_disco_impr_DD_full(UseJack, UseJack?Njacks:Nboots); //constructor sets it to 0
	distr_t agm2_disco_impr_DD_full_ELM(UseJack, UseJack?Njacks:Nboots); //constructor sets it to 0
		 
	for(int t=1; t<= Tdata; t++) {
	  agm2_disco_impr_DD_full= agm2_disco_impr_DD_full + 4.0*w(t,Simps_ord)*pow(alpha,2)*Zv*Zv*disco_impr_DD_distr.distr_list[t]*Ker.distr_list[t];
	  agm2_disco_impr_DD_full_ELM = agm2_disco_impr_DD_full_ELM + 4.0*w(t,Simps_ord)*pow(alpha,2)*Zv*Zv*disco_impr_DD_distr.distr_list[t]*Ker_ELM.distr_list[t];
	}
	

	 agm2_disco_impr_DD_full_Tdata.distr_list.push_back(agm2_disco_impr_DD_full);
	 agm2_disco_impr_DD_full_ELM_Tdata.distr_list.push_back(agm2_disco_impr_DD_full_ELM);
      }

      //push_back the result

      agm2_disco_impr_lightDD_W.distr_list.push_back(agm2_disco_impr_DD_W);
      agm2_disco_impr_lightDD_SD.distr_list.push_back(agm2_disco_impr_DD_SD);
      agm2_disco_impr_lightDD_W_ELM.distr_list.push_back(agm2_disco_impr_DD_W_ELM);
      agm2_disco_impr_lightDD_SD_ELM.distr_list.push_back(agm2_disco_impr_DD_SD_ELM);

      //Print Tdata
      Print_To_File({}, {Tdata_vec_disco, agm2_disco_impr_DD_full_Tdata.ave(), agm2_disco_impr_DD_full_Tdata.err()}, "../data/gm2/light/disco/agm2_impr_DD_Tdata_"+V_light_1.Tag[i_ens]+".dat.t", "", "#id  Tdata   ag2m ");
      Print_To_File({}, {Tdata_vec_disco, agm2_disco_impr_DD_full_ELM_Tdata.ave(), agm2_disco_impr_DD_full_ELM_Tdata.err()}, "../data/gm2/light/disco/agm2_impr_DD_Tdata_ELM_"+V_light_1.Tag[i_ens]+".dat.t", "", "#id  Tdata   ag2m");
    }
    


  


  
    //#################################################################################
    //#################################################################################
    //#################################################################################
    //#################################################################################
    //fit lattice data using analytic representation for V(t)_light


    //#########################################  TWISTED MASS CORRELATOR ###################################################################################################


    if(!Skip_total_light_calc) {

      bootstrap_fit<fit_par,ipar> bf(Njacks);
      bf.set_warmup_lev(4); //sets warmup

      //######################FIT INTERVALS####################
      //##############INCLUDE only times t for which dC(t)/C(t) < 0.1#####################

      
      int Tfit_min;
      if(!Fit_only_light_corr_tail) {
      Tfit_min= 5;
      if(V_light_1.Tag[i_ens].substr(1,1) == "D") Tfit_min++;
      }
      else { //start fitting from 1 fermi
	Tfit_min = (int)(1.2*fm_to_inv_Gev/a_distr.ave());
      }
      int Tfit_max= Tfit_min;  
      bool Found_error_less_10_percent=false;
      while(!Found_error_less_10_percent && Tfit_max < Corr.Nt/2) {
   
	if( (Za*Za*V_light_distr_tm_corr.distr_list[Tfit_max]).err()/fabs( (Za*Za*V_light_distr_tm_corr.distr_list[Tfit_max]).ave()) <  0.10) Tfit_max++;
	else Found_error_less_10_percent=true;
      }
      int Tfit_points = Tfit_max +1 - Tfit_min;
      //###################################################################################

  
      bf.Set_number_of_measurements(Tfit_points);
      bf.Set_verbosity(verbosity);

      cout<<"FITTING WITH PIPI+DUAL REPRESENTATION ENSEMBLE: "<<V_light_1.Tag[i_ens]<<endl;
      cout<<"Tmin: "<<Tfit_min<<", Tmax: "<<Tfit_max<<endl;

      //initial guesses are based on the results on the old ETMC confs
      bf.Add_par("Rd", 1.1, 0.1);
      bf.Set_limits("Rd", 0.7, 1.8);
      if(!Use_Mpi_OS) {
	bf.Add_par("Ed", 2.5, 0.2);
	bf.Set_limits("Ed", 0.7, 5.0);
	bf.Add_par("Mrho",5.5, 0.1);
	bf.Set_limits("Mrho", 4.5 , 7.0);
      }
      else {
	bf.Add_par("Ed", 1.3, 0.2);
	bf.Set_limits("Ed", 0.5, 3.0);
	bf.Add_par("Mrho",3.0, 0.05);
	bf.Set_limits("Mrho", 2.4, 5.0);
      }
      bf.Add_par("gpi", 1.0, 0.01);
      bf.Add_par("kappa", -3.0, 0.1);
      bf.Add_par("Mpi_corr", 1.0, 0.05);
      bf.Set_limits("gpi",0.6, 1.5);
      bf.Set_limits("Mpi_corr", 0.5, 1.5); 
      //bf.Set_limits("kappa", -5.0, 1.0);


      //bf.Fix_par("gpi",1.0);
      bf.Fix_par("kappa", 0.0);
      bf.Fix_par("Mpi_corr", 1.0);
      if(Fit_only_light_corr_tail) { //set to zero parameters of dual representation

	bf.Fix_par("Ed", 0.0);
	bf.Fix_par("Rd", 0.0);

      }

      map<pair<pair<double,double>, pair<double,double>>,Vfloat> Energy_lev_list;
  

      bf.ansatz =  [&](const fit_par &p, const ipar &ip) -> double {


		     double Pi_M = Use_Mpi_OS?ip.Mp_OS:ip.Mp;

		     double GPI = p.gpi*p.Mrho*Pi_M/ip.fp;
		     Vfloat Knpp;
		     pair<double,double> Mass_par= make_pair(p.Mrho*Pi_M, p.Mpi_corr*Pi_M);
		     pair<double,double> Couplings = make_pair(GPI, p.kappa);
		     pair< pair<double,double>,pair<double, double>> input_pars = make_pair( Mass_par, Couplings);
		     map<pair<pair<double,double>, pair<double, double>>,Vfloat>::iterator it;
		     it= Energy_lev_list.find(input_pars);
		     if(it != Energy_lev_list.end()) Knpp= it->second;
		     else {
		       LL.Find_pipi_energy_lev(ip.L,p.Mrho*Pi_M, GPI, p.Mpi_corr*Pi_M, p.kappa, Knpp);
		       //add to Energy_lev_list
		       Energy_lev_list.insert( make_pair(input_pars, Knpp));
		     }
		
		     return Qfact*LL.V_pipi(ip.t, ip.L, p.Mrho*Pi_M, GPI, p.Mpi_corr*Pi_M, p.kappa, Knpp) + LL.Vdual(ip.t, p.Mrho*Pi_M, p.Ed*Pi_M, p.Rd);
		   };
      bf.measurement = [&](const fit_par& p,const ipar& ip) -> double {
			 return ip.V_light;
		       };
      bf.error =  [&](const fit_par& p,const ipar &ip) -> double {
		    return ip.V_light_err;
		  };

      //fill the data
      vector<vector<ipar>> data(Njacks);
      //allocate space for output result
      boot_fit_data<fit_par> Bt_fit;

  

      for(auto &data_iboot: data) data_iboot.resize(Tfit_points);
  
      for(int ijack=0;ijack<Njacks;ijack++) {
	for(int t=Tfit_min;t<Tfit_points+Tfit_min;t++) {
	  int tt=t-Tfit_min;
	  data[ijack][tt].V_light = pow(Za.distr[ijack],2)*V_light_distr_tm_corr.distr_list[t].distr[ijack];
	  data[ijack][tt].V_light_err = (Za*Za*V_light_distr_tm_corr.distr_list[t]).err();
	  data[ijack][tt].L = L_info.L;
	  data[ijack][tt].Mp = Mpi.distr[ijack];
	  data[ijack][tt].Mp_OS = Mpi_OS.distr[ijack];
	  data[ijack][tt].t = t;
	  data[ijack][tt].fp = fp.distr[ijack];
	} 
      }

    
      //add prior values
      for(int ijack=0;ijack<Njacks;ijack++) {
	bf.ib= &ijack;
	//bf.Append_to_prior("gpi", 1.0, 0.2);
	//bf.Append_to_prior("kappa", 0.0, 3.0);
      }

    
      //append
      bf.Append_to_input_par(data);

    
    
      //fit
      Bt_fit= bf.Perform_bootstrap_fit();
    
      Bt_fit.ch2_ave();


      //get energy levels and store the distribution of the ground state (2-pi)
      distr_t gs_tm(UseJack);
      for(int ijack=0;ijack<Njacks;ijack++) {
	double M= Mpi.distr[ijack];
	double MR= Bt_fit.par[ijack].Mrho;
	double FP= fp.distr[ijack];
	double GG= Bt_fit.par[ijack].gpi;
	Vfloat Enn;
	LL.Find_pipi_energy_lev((double)L_info.L,M*MR, GG*MR*M/FP, M, 0.0, Enn);
	gs_tm.distr.push_back(2.0*sqrt( pow(Enn[0],2) + M*M));
      }

      //###################################################
      //print fitted func
    
      double step_size_time= 0.1;
      int nsteps= (int)((1.0/step_size_time)*Corr.Nt); //500 points in time
      cout<<"Printing fit function with: istep: "<<step_size_time<<" and nstep: "<<nsteps<<endl;
      distr_t_list func(UseJack, nsteps);
      Vfloat times;
    

      //generate data
      for(int istep=0; istep<nsteps;istep++) {
	ipar my_ipar;
	my_ipar.t = step_size_time*(1.0 + istep);
	my_ipar.L = (double)L_info.L;
	times.push_back(my_ipar.t);
	for(int ijack=0;ijack<Njacks;ijack++) {
	  my_ipar.Mp=Mpi.distr[ijack];
	  my_ipar.Mp_OS=Mpi_OS.distr[ijack];
	  my_ipar.fp = fp.distr[ijack];
	  func.distr_list[istep].distr.push_back( bf.ansatz( Bt_fit.par[ijack], my_ipar ));
	  
	  
	  
	}
      }

      //print to file
      Print_To_File({}, {times, func.ave(), func.err()}, "../data/gm2/light/tm/Corr_anal_rep_"+V_light_1.Tag[i_ens]+".dat.t", "", "#id_col  time  val   err  (fit interval: ["+to_string(Tfit_min)+","+to_string(Tfit_points+Tfit_min-1)+"])  a = "+to_string_with_precision(a_distr.ave()/fm_to_inv_Gev,6)+" +- "+to_string_with_precision(a_distr.err()/fm_to_inv_Gev, 6)+" fm");
      //####################################################
    
    
      //push_back jack distribution of fitted params
      par_list_anal_repr.push_back(Bt_fit.par);
      cout<<"SIZE OF ENERGY LEVEL MAP FOR TM-FIT IS: "<<Energy_lev_list.size()<<endl;
      cout<<"#########################END ENSEMBLE FIT"<<endl;

      //##################################################################################################################################################################################################


  
      //#########################################################  OSTERWALDER-SEILER CORRELATOR  ########################################################################################################




      bootstrap_fit<fit_par,ipar> bf_OS(Njacks);
      bf_OS.set_warmup_lev(2); //sets warmup

      //######################FIT INTERVALS####################
      //##############INCLUDE only times t for which dC(t)/C(t) < 0.1#####################

      if(!Fit_only_light_corr_tail) {
      Tfit_min= 6;
      if(V_light_1.Tag[i_ens].substr(1,1) == "D") Tfit_min++;
      }
      else { //start fitting from 1 fm
	Tfit_min = (int)(1.2*fm_to_inv_Gev/a_distr.ave());
      }
      Tfit_max= Tfit_min;  
      Found_error_less_10_percent=false;
      while(!Found_error_less_10_percent && Tfit_max < Corr.Nt/2) {
   
	if( (Zv*Zv*V_light_distr_OS_corr.distr_list[Tfit_max]).err()/fabs( (Zv*Zv*V_light_distr_OS_corr.distr_list[Tfit_max]).ave()) <  0.10) Tfit_max++;
	else Found_error_less_10_percent=true;
      }
      Tfit_points = Tfit_max +1 - Tfit_min;
      //###################################################################################

  
      bf_OS.Set_number_of_measurements(Tfit_points);
      bf_OS.Set_verbosity(verbosity);

      cout<<"FITTING WITH PIPI+DUAL REPRESENTATION ENSEMBLE  (OS correlator): "<<V_light_1.Tag[i_ens]<<endl;
      cout<<"Tmin: "<<Tfit_min<<", Tmax: "<<Tfit_max<<endl;

      //initial guesses are based on the results on the old ETMC confs
      bf_OS.Add_par("Rd", 1.1, 0.1);
      bf_OS.Set_limits("Rd", 0.6, 2.0);
      if(!Use_Mpi_OS) {
	bf_OS.Add_par("Ed", 2.1, 0.2);
	bf_OS.Set_limits("Ed", 0.5, 5.0);
	bf_OS.Add_par("Mrho",5.5, 0.1);
	bf_OS.Set_limits("Mrho", 4.5 , 7.0);
      }
      else {
	bf_OS.Add_par("Ed", 1.3, 0.2);
	bf_OS.Set_limits("Ed", 0.5, 3.0);
	bf_OS.Add_par("Mrho",3.0, 0.05);
	bf_OS.Set_limits("Mrho", 2.4, 5.0);
      }
      bf_OS.Add_par("gpi", 1.0, 0.01);
      bf_OS.Add_par("kappa", -3.0, 0.1);
      bf_OS.Add_par("Mpi_corr", 1.0, 0.05);
      bf_OS.Set_limits("gpi",0.6, 1.5);
      bf_OS.Set_limits("Mpi_corr", 0.5, 1.5);
      //bf_OS.Set_limits("kappa", -5.0, 1.0);


      //bf_OS.Fix_par("gpi",1.0);
      bf_OS.Fix_par("Mpi_corr", 1.0);
      bf_OS.Fix_par("kappa", 0.0);
      if(Fit_only_light_corr_tail) { //fix to zero params of dual representation
	bf_OS.Fix_par("Ed",0.0);
	bf_OS.Fix_par("Rd", 0.0);
      }

      map<pair<pair<double,double>, pair<double,double>>,Vfloat> Energy_lev_list_OS;
  

      bf_OS.ansatz =  [&](const fit_par &p, const ipar &ip) -> double {


			double Pi_M = Use_Mpi_OS?ip.Mp_OS:ip.Mp;

			double GPI = p.gpi*p.Mrho*Pi_M/ip.fp;
			Vfloat Knpp;
			pair<double,double> Mass_par= make_pair(p.Mrho*Pi_M, p.Mpi_corr*Pi_M);
			pair<double,double> Couplings = make_pair(GPI, p.kappa);
			pair< pair<double,double>,pair<double, double>> input_pars = make_pair( Mass_par, Couplings);
			map<pair<pair<double,double>, pair<double, double>>,Vfloat>::iterator it;
			it= Energy_lev_list_OS.find(input_pars);
			if(it != Energy_lev_list_OS.end()) Knpp= it->second;
			else {
			  LL.Find_pipi_energy_lev(ip.L,p.Mrho*Pi_M, GPI, p.Mpi_corr*Pi_M, p.kappa, Knpp);
			  //add to Energy_lev_list
			  Energy_lev_list_OS.insert( make_pair(input_pars, Knpp));
			}
		
			return Qfact*LL.V_pipi(ip.t, ip.L, p.Mrho*Pi_M, GPI, p.Mpi_corr*Pi_M, p.kappa, Knpp) + LL.Vdual(ip.t, p.Mrho*Pi_M, p.Ed*Pi_M, p.Rd);
		      };
      bf_OS.measurement = [&](const fit_par& p,const ipar& ip) -> double {
			    return ip.V_light;
			  };
      bf_OS.error =  [&](const fit_par& p,const ipar &ip) -> double {
		       return ip.V_light_err;
		     };

      //fill the data
      vector<vector<ipar>> data_OS(Njacks);
      //allocate space for output result
      boot_fit_data<fit_par> Bt_fit_OS;

  

      for(auto &data_iboot: data_OS) data_iboot.resize(Tfit_points);
  
      for(int ijack=0;ijack<Njacks;ijack++) {
	for(int t=Tfit_min;t<Tfit_points+Tfit_min;t++) {
	  int tt=t-Tfit_min;
	  data_OS[ijack][tt].V_light = pow(Zv.distr[ijack],2)*V_light_distr_OS_corr.distr_list[t].distr[ijack];
	  data_OS[ijack][tt].V_light_err = (Zv*Zv*V_light_distr_OS_corr.distr_list[t]).err();
	  data_OS[ijack][tt].L = L_info.L;
	  data_OS[ijack][tt].Mp = Mpi.distr[ijack];
	  data_OS[ijack][tt].Mp_OS = Mpi_OS.distr[ijack];
	  data_OS[ijack][tt].t = t;
	  data_OS[ijack][tt].fp = fp.distr[ijack];
	} 
      }

    
      //add prior values
      for(int ijack=0;ijack<Njacks;ijack++) {
	bf_OS.ib= &ijack;
	//bf.Append_to_prior("gpi", 1.0, 0.2);
	//bf_OS.Append_to_prior("kappa", 0.0, 3.0);
      }

    
      //append
      bf_OS.Append_to_input_par(data_OS);

    
    
      //fit
      Bt_fit_OS= bf_OS.Perform_bootstrap_fit();
    
      Bt_fit_OS.ch2_ave();


      //get energy levels and store the distribution of the ground state (2-pi)
      distr_t gs_OS(UseJack);
      distr_t mrho_GS(UseJack);
      for(int ijack=0;ijack<Njacks;ijack++) {
	double M= Mpi.distr[ijack];
	double MR= Bt_fit_OS.par[ijack].Mrho;
	double FP= fp.distr[ijack];
	mrho_GS.distr.push_back( M*MR);
	double GG= Bt_fit_OS.par[ijack].gpi;
	Vfloat Enn;
	LL.Find_pipi_energy_lev((double)L_info.L,M*MR, GG*MR*M/FP, M, 0.0, Enn);
	gs_OS.distr.push_back(2.0*sqrt( pow(Enn[0],2) + M*M));
      }

      if(V_light_1.Tag[i_ens] == "cD211a.054.96") Mrho_from_GS[2] = mrho_GS;
      else if(V_light_1.Tag[i_ens] == "cC211a.06.80") Mrho_from_GS[1] = mrho_GS;
      else if(V_light_1.Tag[i_ens] == "cB211b.072.64") Mrho_from_GS[0]= mrho_GS;

      //###################################################
      //print fitted func

      cout<<"Printing fit function (OS) with: istep: "<<step_size_time<<" and nstep: "<<nsteps<<endl;
      distr_t_list func_OS(UseJack, nsteps);
      Vfloat times_OS;

      //generate data
      for(int istep=0; istep<nsteps;istep++) {
	ipar my_ipar;
	my_ipar.t = step_size_time*(1.0 + istep);
	my_ipar.L = (double)L_info.L;
	times_OS.push_back(my_ipar.t);
	for(int ijack=0;ijack<Njacks;ijack++) {
	  my_ipar.Mp=Mpi.distr[ijack];
	  my_ipar.Mp_OS=Mpi_OS.distr[ijack];
	  my_ipar.fp = fp.distr[ijack];
	  func_OS.distr_list[istep].distr.push_back( bf.ansatz( Bt_fit_OS.par[ijack], my_ipar ));
	}
      }

      //print to file
      Print_To_File({}, {times_OS, func_OS.ave(), func_OS.err()}, "../data/gm2/light/OS/Corr_anal_rep_"+V_light_1.Tag[i_ens]+".dat.t", "", "#id_col  time  val   err  (fit interval: ["+to_string(Tfit_min)+","+to_string(Tfit_points+Tfit_min-1)+"])  a = "+to_string_with_precision(a_distr.ave()/fm_to_inv_Gev,6)+" +- "+to_string_with_precision(a_distr.err()/fm_to_inv_Gev, 6)+" fm");
      //####################################################
    
    
      //push_back jack distribution of fitted params
      par_list_anal_repr_OS.push_back(Bt_fit_OS.par);
      cout<<"SIZE OF ENERGY LEVEL MAP FOR OS-FIT IS: "<<Energy_lev_list_OS.size()<<endl;
      cout<<"#########################END ENSEMBLE FIT"<<endl;

      //#################################################################################################################################################################################################

    

      //BOUNDING ON INTERMEDIATE WINDOW #################
      distr_t amu_W_bound_tm(UseJack), amu_W_bound_OS(UseJack);

      Bounding_amu_W(amu_W_bound_tm, V_light_distr, a_distr, "../data/gm2/light/tm/window_W_Tdata_"+V_light_1.Tag[i_ens], Za, 0.8*gs_tm);
      Bounding_amu_W(amu_W_bound_OS, V_light_OS_distr, a_distr, "../data/gm2/light/OS/window_W_Tdata_"+V_light_1.Tag[i_ens], Zv, 0.8*gs_OS);
      //push_back the result
      amu_W_bounding_tm_list.distr_list.push_back(amu_W_bound_tm);
      amu_W_bounding_OS_list.distr_list.push_back(amu_W_bound_OS);
      //################################################





      //######################   PI(Q^2) analysis   ###########################


      distr_t_list PI_Q2_tm, PI_Q2_OS;
      distr_t_list PI_Q2_tm_pert_sub, PI_Q2_OS_pert_sub;

    
      Vint Tdatas_opt_tm, Tdatas_opt_OS;
      Vint Tdatas_opt_tm_pert_sub, Tdatas_opt_OS_pert_sub;

    
      // Apply bounding method
      Bounding_PI_q2(PI_Q2_tm, Za*Za*V_light_distr, a_distr, "../data/PI_Q2/light/tm/PI_Q2_Tdata_"+V_light_1.Tag[i_ens]  , Tdatas_opt_tm, gs_tm);
  
      Bounding_PI_q2(PI_Q2_OS, Zv*Zv*V_light_OS_distr, a_distr, "../data/PI_Q2/light/OS/PI_Q2_Tdata_"+V_light_1.Tag[i_ens], Tdatas_opt_OS, gs_OS);
  
      //include perturbative subtraction
      Bounding_PI_q2(PI_Q2_tm_pert_sub, Za*Za*V_light_distr_tm_corr, a_distr, "../data/PI_Q2/light/tm/PI_Q2_Tdata_"+V_light_1.Tag[i_ens]+"_pert_sub", Tdatas_opt_tm_pert_sub, gs_tm);
  
      Bounding_PI_q2(PI_Q2_OS_pert_sub, Zv*Zv*V_light_distr_OS_corr, a_distr, "../data/PI_Q2/light/OS/PI_Q2_Tdata_"+V_light_1.Tag[i_ens]+"_pert_sub", Tdatas_opt_OS_pert_sub, gs_OS);


      //push_back light corr for bounding on disconnected
      if(V_light_1.Tag[i_ens] == "cD211a.054.96" || V_light_1.Tag[i_ens] == "cC211a.06.80" || V_light_1.Tag[i_ens] == "cB211b.072.64") {
	int id_disco_ens=0;
	if(V_light_1.Tag[i_ens]== "cD211a.054.96") id_disco_ens=2;
	else if(V_light_1.Tag[i_ens]== "cC211a.06.80") id_disco_ens=1;
	else if(V_light_1.Tag[i_ens]== "cB211b.072.64") id_disco_ens=0;
	else crash("what?");
       
	CORR_LIGHT_FOR_PI_Q2[id_disco_ens]= Zv*Zv*V_light_OS_distr;
	LIGHT_PI_Q2_FOR_DISCO[id_disco_ens] = PI_Q2_OS;

	if(Include_light_disco && Include_off_diagonal_disco) {
	  if(!Found_disco_impr_D_ens) crash("Disconnected for PI(Q^2) not found in light ens");
	  CORR_DISCO_FOR_PI_Q2[id_disco_ens] = CORR_DISCO_FOR_PI_Q2[id_disco_ens] + Zv*Zv*disco_impr_D_distr;
	  L_list_disco_PI_Q2[id_disco_ens] = L_info.L;
	  Mpi_fit_disc_PI_Q2.distr_list[id_disco_ens] = Mpi ;
	  fpi_fit_disc_PI_Q2.distr_list[id_disco_ens] =fp;
	  a_disc_PI_Q2.distr_list[id_disco_ens] = a_distr;
	  Ens_list_disc_PI_Q2[id_disco_ens] = V_light_1.Tag[i_ens];

	}
	else crash("what?");

      }

  


      
      Vfloat Tcut_f_tm, Tcut_f_OS, Tcut_f_tm_pert_sub, Tcut_f_OS_pert_sub;

    
      int Q_size= Qs2.size();
      for(int q=0; q < Q_size;q++) {
	Tcut_f_tm.push_back( (double)Tdatas_opt_tm[q]);
	Tcut_f_OS.push_back( (double)Tdatas_opt_OS[q]);
	Tcut_f_tm_pert_sub.push_back( (double)Tdatas_opt_tm_pert_sub[q]);
	Tcut_f_OS_pert_sub.push_back( (double)Tdatas_opt_OS_pert_sub[q]);
      }
    
      //push_back and print to file

      Add_ens_val_PI_q2( PI_Q2_light_tm, PI_Q2_tm);
      Add_ens_val_PI_q2( PI_Q2_light_OS, PI_Q2_OS);
      Add_ens_val_PI_q2( PI_Q2_light_tm_pert_sub, PI_Q2_tm_pert_sub);
      Add_ens_val_PI_q2( PI_Q2_light_OS_pert_sub, PI_Q2_OS_pert_sub);


      Print_To_File({}, {Qs2, PI_Q2_tm.ave(), PI_Q2_tm.err(), Tcut_f_tm } , "../data/PI_Q2/light/tm/PI_Q2_extr_"+V_light_1.Tag[i_ens]+".t", "", "# Q2[GeV2]   PI(Q^2)   Tcut_f. Mpi: "+to_string_with_precision(Mpi.ave(),6)+" +- "+to_string_with_precision(Mpi.err(), 6)+" lightest 2pi-state: "+to_string_with_precision(gs_tm.ave(),6)+" +- "+to_string_with_precision(gs_tm.err(),6));
    Print_To_File({}, {Qs2, PI_Q2_OS.ave(), PI_Q2_OS.err(), Tcut_f_OS } , "../data/PI_Q2/light/OS/PI_Q2_extr_"+V_light_1.Tag[i_ens]+".t", "", "# Q2[GeV2]   PI(Q^2)   Tcut_f. Mpi: "+to_string_with_precision(Mpi.ave(),6)+" +- "+to_string_with_precision(Mpi.err(), 6)+" lightest 2pi-state: "+to_string_with_precision(gs_OS.ave(),6)+" +- "+to_string_with_precision(gs_OS.err(),6));
    Print_To_File({}, {Qs2, PI_Q2_tm_pert_sub.ave(), PI_Q2_tm_pert_sub.err(), Tcut_f_tm_pert_sub } , "../data/PI_Q2/light/tm/PI_Q2_extr_pert_sub_"+V_light_1.Tag[i_ens]+".t", "", "# Q2[GeV2]   PI(Q^2)   Tcut_f. Mpi: "+to_string_with_precision(Mpi.ave(),6)+" +- "+to_string_with_precision(Mpi.err(), 6)+" lightest 2pi-state: "+to_string_with_precision(gs_tm.ave(),6)+" +- "+to_string_with_precision(gs_tm.err(),6));
    Print_To_File({}, {Qs2, PI_Q2_OS_pert_sub.ave(), PI_Q2_OS_pert_sub.err(), Tcut_f_OS_pert_sub } , "../data/PI_Q2/light/OS/PI_Q2_extr_pert_sub_"+V_light_1.Tag[i_ens]+".t", "", "# Q2[GeV2]   PI(Q^2)   Tcut_f. Mpi: "+to_string_with_precision(Mpi.ave(),6)+" +- "+to_string_with_precision(Mpi.err(), 6)+" lightest 2pi-state: "+to_string_with_precision(gs_OS.ave(),6)+" +- "+to_string_with_precision(gs_OS.err(),6));




    //#######################################################################

    }
  }
  
  cout<<"light quark correlator analyzed!"<<endl;



 

  

  distr_t_list Edual(UseJack), Rdual(UseJack), Mrho(UseJack), gpi(UseJack), Kappa(UseJack), M_corr(UseJack);
  distr_t_list Edual_OS(UseJack), Rdual_OS(UseJack), Mrho_OS(UseJack), gpi_OS(UseJack), Kappa_OS(UseJack), M_corr_OS(UseJack);

  if(!Skip_total_light_calc && Reco_agm2_total) {





    cout<<"####### Reconstructing agm2_light using fit paramters Ed, Rd, Mrho, gpi....."<<endl;


    //#############################################   TWISTED MASS CORRELATOR ##########################################################################

  
    cout<<"#########################  TWISTED MASS CORRELATOR ###########################"<<endl;



  
 
    for(int i_ens=0; i_ens<Nens_light;i_ens++) {
      cout<<"######### ENSEMBLE: "<<V_light_1.Tag[i_ens]<<endl;
      LatticeInfo L_info;
      L_info.LatInfo_new_ens(V_light_1.Tag[i_ens]);
  

      distr_t a_distr = a_distr_list[i_ens];
    
    
      distr_t Ed(UseJack), Rd(UseJack), Mr(UseJack), g(UseJack), kap(UseJack), mpi_corr(UseJack), agm2_dual(UseJack), mp(UseJack), agm2_2L_dual(UseJack), agm2_infL_dual(UseJack), agm2_Lprime_dual(UseJack);
      distr_t agm2_W_FSEs_GS_distr(UseJack);
      for(int ijack=0;ijack<Njacks;ijack++) {
	//retrieve fit_parameters
	fit_par my_fit_par = par_list_anal_repr[i_ens][ijack];
	Ed.distr.push_back(my_fit_par.Ed);
	Rd.distr.push_back(my_fit_par.Rd);
	Mr.distr.push_back(my_fit_par.Mrho);
	g.distr.push_back(my_fit_par.gpi);
	kap.distr.push_back(my_fit_par.kappa);
	mpi_corr.distr.push_back(my_fit_par.Mpi_corr);
       

	double Mp,fp, csi;
	double L= 1.0*L_list[i_ens];
	fp= fp_fit[i_ens].distr[ijack];
	//double csi = pow(Mp,2)/pow(4.0*M_PI*fp,2);
	if(!Use_Mpi_OS) Mp=Mpi_fit.distr_list[i_ens].distr[ijack];
	else Mp=Mpi_OS_fit.distr_list[i_ens].distr[ijack];
	mp.distr.push_back(Mp);



    
	Vfloat Elev;
    
	LL.Find_pipi_energy_lev(L,my_fit_par.Mrho*Mp, my_fit_par.gpi*(my_fit_par.Mrho*Mp/fp), my_fit_par.Mpi_corr*Mp, my_fit_par.kappa, Elev);

	//Find energy levs corresponding to 2L;

	Vfloat Elev_2L;

	LL.Find_pipi_energy_lev(1.5*L, my_fit_par.Mrho*Mp, my_fit_par.gpi*(my_fit_par.Mrho*Mp/fp), my_fit_par.Mpi_corr*Mp, my_fit_par.kappa, Elev_2L);


	Vfloat Elev_MpiL_4dot2;

	double Lprime= 4.2/(!Use_Mpi_OS?Mpi_fit.ave(i_ens):Mpi_OS_fit.ave(i_ens));

	LL.Find_pipi_energy_lev(Lprime, my_fit_par.Mrho*Mp, my_fit_par.gpi*(my_fit_par.Mrho*Mp/fp), my_fit_par.Mpi_corr*Mp, my_fit_par.kappa, Elev_MpiL_4dot2);



	auto FSEs_win_tm = [&](double t) {

			     double kern_val = 4.0*pow(alpha,2)*w(t,Simps_ord)*kernel_K(t, a_distr.distr[ijack] );
			    double func_val =  Qfact*( LL.V_pipi_infL(t, my_fit_par.Mrho*Mp, my_fit_par.gpi*(my_fit_par.Mrho*Mp/fp), my_fit_par.Mpi_corr*Mp, my_fit_par.kappa) - LL.V_pipi(t, L, my_fit_par.Mrho*Mp, my_fit_par.gpi*(my_fit_par.Mrho*Mp/fp), my_fit_par.Mpi_corr*Mp, my_fit_par.kappa, Elev));
			    double win_fact= (  1.0/(1.0 + exp(-2.0*(t*a_distr.distr[ijack]-t0)/Delta)) -  1.0/(1.0 + exp(-2.0*(t*a_distr.distr[ijack]-t1)/Delta)));

			    return kern_val*func_val*win_fact;


			  };



	auto F_int= [&](double t) {


		  
		      //double kern_val = 4.0*pow(alpha,2)*sqrt(m_rho*a_distr.distr[ijack]/(my_fit_par.Mrho*Mp))*kernel_K(t, my_fit_par.Mrho*Mp/m_rho);
		  
		      double kern_val = 4.0*pow(alpha,2)*w(t,Simps_ord)*kernel_K(t, a_distr.distr[ijack] );
		      //double kern_val= 4.0*pow(alpha,2)*kernel_K(t, my_fit_par.Mrho*Mp/m_rho);
		      double func_val =  Qfact*LL.V_pipi(t, L, my_fit_par.Mrho*Mp, my_fit_par.gpi*(my_fit_par.Mrho*Mp/fp),my_fit_par.Mpi_corr*Mp, my_fit_par.kappa, Elev)+ LL.Vdual(t, my_fit_par.Mrho*Mp, my_fit_par.Ed*Mp, my_fit_par.Rd);

		      double F_int_val = kern_val*func_val;
        

		      return F_int_val;
		    };


	auto F_int_2L = [&](double t) {

			  //double kern_val = 4.0*pow(alpha,2)*sqrt(m_rho*a_distr.distr[ijack]/(my_fit_par.Mrho*Mp))*kernel_K(t, my_fit_par.Mrho*Mp/m_rho);
		  
			  double kern_val = 4.0*pow(alpha,2)*w(t,Simps_ord)*kernel_K(t, a_distr.distr[ijack] );
			  //double kern_val= 4.0*pow(alpha,2)*kernel_K(t, my_fit_par.Mrho*Mp/m_rho);
			  double func_val =  Qfact*LL.V_pipi(t, 1.5*L, my_fit_par.Mrho*Mp, my_fit_par.gpi*(my_fit_par.Mrho*Mp/fp), my_fit_par.Mpi_corr*Mp, my_fit_par.kappa, Elev_2L)+ LL.Vdual(t, my_fit_par.Mrho*Mp, my_fit_par.Ed*Mp, my_fit_par.Rd);

			  double F_int_val = kern_val*func_val;
        

			  return F_int_val;
			};

	auto F_int_Lprime = [&](double t) {

			      //  double kern_val = 4.0*pow(alpha,2)*sqrt(m_rho*a_distr.distr[ijack]/(my_fit_par.Mrho*Mp))*kernel_K(t, my_fit_par.Mrho*Mp/m_rho);
		  
			      double kern_val = 4.0*pow(alpha,2)*w(t,Simps_ord)*kernel_K(t, a_distr.distr[ijack] );
			      //double kern_val= 4.0*pow(alpha,2)*kernel_K(t, my_fit_par.Mrho*Mp/m_rho);
			      double func_val =  Qfact*LL.V_pipi(t, Lprime, my_fit_par.Mrho*Mp, my_fit_par.gpi*(my_fit_par.Mrho*Mp/fp), my_fit_par.Mpi_corr*Mp, my_fit_par.kappa,  Elev_MpiL_4dot2)+ LL.Vdual(t, my_fit_par.Mrho*Mp, my_fit_par.Ed*Mp, my_fit_par.Rd);

			      double F_int_val = kern_val*func_val;
        

			      return F_int_val;
			    };

	auto F_infL = [&](double t) {

			//double kern_val = 4.0*pow(alpha,2)*sqrt(m_rho*a_distr.distr[ijack]/(my_fit_par.Mrho*Mp))*kernel_K(t, my_fit_par.Mrho*Mp/m_rho);
		  
			//double kern_val = 4.0*pow(alpha,2)*kernel_K(t, a_distr.distr[ijack] );
			double kern_val= 4.0*pow(alpha,2)*w(t,Simps_ord)*kernel_K(t, my_fit_par.Mrho*Mp/m_rho);

			double func_val =  Qfact*LL.V_pipi_infL(t, my_fit_par.Mrho*Mp, my_fit_par.gpi*(my_fit_par.Mrho*Mp/fp), my_fit_par.Mpi_corr*Mp, my_fit_par.kappa)+ LL.Vdual(t, my_fit_par.Mrho*Mp, my_fit_par.Ed*Mp, my_fit_par.Rd);

			double F_int_val = kern_val*func_val;
        

			return F_int_val;

		      };

	

   

    
    
	//Vfloat Sum_Series_Terms;
	double agm2_summ=0.0;
	double agm2_2L_summ=0.0;
	double agm2_infL_summ=0.0;
	double agm2_Lprime_summ=0.0;
	double agm2_W_FSEs_GS=0.0;
	int Nterms = 2000;
	int Nterms_2L= 2000;
	int Nterms_win= 100;

	for(int iterm=1;iterm<Nterms_win;iterm++) {

	  double F_n = FSEs_win_tm((double)iterm);
	  agm2_W_FSEs_GS += F_n;
	}
	agm2_W_FSEs_GS_distr.distr.push_back(agm2_W_FSEs_GS);


	//prediction for FSE in intermediate window
	
    
	//actual volume
	for(int iterm=1;iterm<Nterms;iterm++) {
	  double F_n = F_int((double)iterm);
	  //Sum_Series 
	  agm2_summ += F_n;}
	agm2_dual.distr.push_back(agm2_summ);
   

	//1.5L volume
	for(int iterm=1;iterm<Nterms_2L;iterm++) {
	  double F_n_2L = F_int_2L((double)iterm);
	  //Sum Series 
	  agm2_2L_summ += F_n_2L;}
	agm2_2L_dual.distr.push_back(agm2_2L_summ);


	//Mpi*L=4.2 volume
	for(int iterm=1;iterm<Nterms_2L;iterm++) {
	  double F_n_Lprime = F_int_Lprime((double)iterm);
	  //Sum Series 
	  agm2_Lprime_summ += F_n_Lprime;}
	agm2_Lprime_dual.distr.push_back(agm2_Lprime_summ);
    

	//inf volume
	for(int iterm=1; iterm<Nterms_2L;iterm++) {
	  double F_n_infL = F_infL((double)iterm);
	  //Sum series
	  agm2_infL_summ += F_n_infL;
	}
	agm2_infL_dual.distr.push_back(agm2_infL_summ);

      
    

      }

     

      //push_back
      Edual.distr_list.push_back(Ed);
      Rdual.distr_list.push_back(Rd);
      Mrho.distr_list.push_back(Mr);
      gpi.distr_list.push_back(g);
      Kappa.distr_list.push_back(kap);
      M_corr.distr_list.push_back(mpi_corr);
      agm2_light_fit.distr_list.push_back(agm2_dual);
      agm2_light_2L_fit.distr_list.push_back(agm2_2L_dual);
      agm2_light_infL_fit.distr_list.push_back(agm2_infL_dual);
      agm2_light_Lprime_fit.distr_list.push_back(agm2_Lprime_dual);



      //print ensemble infos

      //intermediate and SD window info
      cout<<"%%%%%%% Print info from W and SD windows TM-correlator %%%%%%"<<endl;
      cout<<"SD: "<<agm2_light_SD.ave(i_ens)<<" +- "<<agm2_light_SD.err(i_ens)<<endl;
      cout<<"SD (ELM): "<<agm2_light_SD_ELM.ave(i_ens)<<" +- "<<agm2_light_SD_ELM.err(i_ens)<<endl;
      cout<<"W: "<<agm2_light_W.ave(i_ens)<<" +- "<<agm2_light_W.err(i_ens)<<endl;
      cout<<"W (ELM): "<<agm2_light_W_ELM.ave(i_ens)<<" +- "<<agm2_light_W_ELM.err(i_ens)<<endl;
      cout<<"%%%%%%% End windows info %%%%%%%"<<endl;
      cout<<"Edual/Mpi: "<<(Ed).ave()<<" +- "<<(Ed).err()<<endl;
      cout<<"Edual: "<<(Ed*mp).ave()<<" +- "<<(Ed*mp).err()<<endl;
      cout<<"Rdual: "<<Rd.ave()<<" +- "<<Rd.err()<<endl;
      cout<<"Mrho/Mpi: "<<(Mr).ave()<<" +- "<<(Mr).err()<<endl;
      cout<<"Mrho (lattice units): "<<(Mr*mp).ave()<<" +- "<<(Mr*mp).err()<<endl;
      cout<<"Mrho (GeV): "<<(Mr*mp/a_distr).ave()<<" +- "<<(Mr*mp/a_distr).err()<<endl;
      cout<<"Mrho/fp: "<<(Mr*mp/fp_fit[i_ens]).ave()<<" +- "<<(Mr*mp/fp_fit[i_ens]).err()<<endl;
      cout<<"gpi: "<<g.ave()<<" +- "<<g.err()<<endl;
      cout<<"kappa: "<<kap.ave()<<" +- "<<kap.err()<<endl;
      cout<<"Mpi corr: "<<mpi_corr.ave()<<" +- "<<mpi_corr.err()<<endl;
      cout<<"Mpi: "<<mp.ave()<<" +- "<<mp.err()<<endl;
      cout<<"Mpi*L : "<<mp.ave()*L_list[i_ens]<<" +- "<<mp.err()*L_list[i_ens]<<endl;
      cout<<"L: "<<L_list[i_ens]<<endl;
      cout<<"agm2 pp+dual (L): "<<agm2_dual.ave()<<" +- "<<agm2_dual.err()<<endl;
      cout<<"agm2 pp+dual (1.5*L): "<<agm2_2L_dual.ave()<<" +- "<<agm2_2L_dual.err()<<endl;
      cout<<"agm2 pp+dual (Mpi*L=4.2): "<<agm2_Lprime_dual.ave()<<" +- "<<agm2_Lprime_dual.err()<<endl;
      cout<<"agm2 pp+dual (L -> infty): "<<agm2_infL_dual.ave()<<" +- "<<agm2_infL_dual.err()<<endl;
      cout<<"%%%%% FSEs GS-prediction for light-connected contribution to the intermediate window %%%%%%"<<endl;
      cout<<"a_mu^W (infty) - a_mu^W(L): "<<agm2_W_FSEs_GS_distr.ave()<<" +- "<<agm2_W_FSEs_GS_distr.err()<<endl;
      cout<<"#######################################"<<endl;
    
    



    }


    //########################################      OSTERWALDER-SEILER CORRELATOR        ########################################################



    cout<<"#########################  OSTERWALDER-SEILER CORRELATOR ###########################"<<endl;



  
 
    for(int i_ens=0; i_ens<Nens_light;i_ens++) {
      cout<<"######### ENSEMBLE: "<<V_light_1.Tag[i_ens]<<endl;
      LatticeInfo L_info;
      L_info.LatInfo_new_ens(V_light_1.Tag[i_ens]);
  

      distr_t a_distr = a_distr_list[i_ens];
    
      distr_t Ed(UseJack), Rd(UseJack), Mr(UseJack), g(UseJack), kap(UseJack), mpi_corr(UseJack), agm2_dual(UseJack), mp(UseJack), agm2_2L_dual(UseJack), agm2_infL_dual(UseJack), agm2_Lprime_dual(UseJack);
      distr_t agm2_W_FSEs_GS_distr(UseJack);
      for(int ijack=0;ijack<Njacks;ijack++) {
	//retrieve fit_parameters
	fit_par my_fit_par = par_list_anal_repr_OS[i_ens][ijack];
	Ed.distr.push_back(my_fit_par.Ed);
	Rd.distr.push_back(my_fit_par.Rd);
	Mr.distr.push_back(my_fit_par.Mrho);
	g.distr.push_back(my_fit_par.gpi);
	kap.distr.push_back(my_fit_par.kappa);
	mpi_corr.distr.push_back(my_fit_par.Mpi_corr);

   

	double Mp,fp, csi;
	double L= 1.0*L_list[i_ens];
	fp= fp_fit[i_ens].distr[ijack];
	//double csi = pow(Mp,2)/pow(4.0*M_PI*fp,2);
	if(!Use_Mpi_OS) Mp=Mpi_fit.distr_list[i_ens].distr[ijack];
	else Mp=Mpi_OS_fit.distr_list[i_ens].distr[ijack];
	mp.distr.push_back(Mp);

    
	Vfloat Elev;
    
	LL.Find_pipi_energy_lev(L,my_fit_par.Mrho*Mp, my_fit_par.gpi*(my_fit_par.Mrho*Mp/fp), my_fit_par.Mpi_corr*Mp, my_fit_par.kappa, Elev);

	//Find energy levs corresponding to 2L;

	Vfloat Elev_2L;

	LL.Find_pipi_energy_lev(1.5*L, my_fit_par.Mrho*Mp, my_fit_par.gpi*(my_fit_par.Mrho*Mp/fp),  my_fit_par.Mpi_corr*Mp, my_fit_par.kappa, Elev_2L);


	Vfloat Elev_MpiL_4dot2;

	double Lprime= 4.2/(!Use_Mpi_OS?Mpi_fit.ave(i_ens):Mpi_OS_fit.ave(i_ens));

	LL.Find_pipi_energy_lev(Lprime, my_fit_par.Mrho*Mp, my_fit_par.gpi*(my_fit_par.Mrho*Mp/fp),  my_fit_par.Mpi_corr*Mp, my_fit_par.kappa, Elev_MpiL_4dot2);

    


	auto FSEs_win_OS = [&](double t) {

			    double kern_val = 4.0*pow(alpha,2)*w(t,Simps_ord)*kernel_K(t, a_distr.distr[ijack] );
			    double func_val =  Qfact*( LL.V_pipi_infL(t, my_fit_par.Mrho*Mp, my_fit_par.gpi*(my_fit_par.Mrho*Mp/fp),  my_fit_par.Mpi_corr*Mp, my_fit_par.kappa) - LL.V_pipi(t, L, my_fit_par.Mrho*Mp, my_fit_par.gpi*(my_fit_par.Mrho*Mp/fp),  my_fit_par.Mpi_corr*Mp, my_fit_par.kappa, Elev));
			    double win_term = (  1.0/(1.0 + exp(-2.0*(t*a_distr.distr[ijack]-t0)/Delta)) -  1.0/(1.0 + exp(-2.0*(t*a_distr.distr[ijack]-t1)/Delta)));
			    

			    return kern_val*func_val*win_term;


			  };


	auto F_int= [&](double t) {


		  
		      //double kern_val = 4.0*pow(alpha,2)*sqrt(m_rho*a_distr.distr[ijack]/(my_fit_par.Mrho*Mp))*kernel_K(t, my_fit_par.Mrho*Mp/m_rho);
		  
		      double kern_val = 4.0*pow(alpha,2)*w(t,Simps_ord)*kernel_K(t, a_distr.distr[ijack] );
		      //double kern_val= 4.0*pow(alpha,2)*kernel_K(t, my_fit_par.Mrho*Mp/m_rho);
		      double func_val =  Qfact*LL.V_pipi(t, L, my_fit_par.Mrho*Mp, my_fit_par.gpi*(my_fit_par.Mrho*Mp/fp),  my_fit_par.Mpi_corr*Mp, my_fit_par.kappa, Elev)+ LL.Vdual(t, my_fit_par.Mrho*Mp, my_fit_par.Ed*Mp, my_fit_par.Rd);

		      double F_int_val = kern_val*func_val;
        

		      return F_int_val;
		    };


	auto F_int_2L = [&](double t) {

			  //double kern_val = 4.0*pow(alpha,2)*sqrt(m_rho*a_distr.distr[ijack]/(my_fit_par.Mrho*Mp))*kernel_K(t, my_fit_par.Mrho*Mp/m_rho);
		  
			  double kern_val = 4.0*pow(alpha,2)*w(t,Simps_ord)*kernel_K(t, a_distr.distr[ijack] );
			  //double kern_val= 4.0*pow(alpha,2)*kernel_K(t, my_fit_par.Mrho*Mp/m_rho);
			  double func_val =  Qfact*LL.V_pipi(t, 1.5*L, my_fit_par.Mrho*Mp, my_fit_par.gpi*(my_fit_par.Mrho*Mp/fp),  my_fit_par.Mpi_corr*Mp, my_fit_par.kappa, Elev_2L)+ LL.Vdual(t, my_fit_par.Mrho*Mp, my_fit_par.Ed*Mp, my_fit_par.Rd);

			  double F_int_val = kern_val*func_val;
        

			  return F_int_val;
			};

	auto F_int_Lprime = [&](double t) {

			      //  double kern_val = 4.0*pow(alpha,2)*sqrt(m_rho*a_distr.distr[ijack]/(my_fit_par.Mrho*Mp))*kernel_K(t, my_fit_par.Mrho*Mp/m_rho);
		  
			      double kern_val = 4.0*pow(alpha,2)*w(t,Simps_ord)*kernel_K(t, a_distr.distr[ijack] );
			      //double kern_val= 4.0*pow(alpha,2)*kernel_K(t, my_fit_par.Mrho*Mp/m_rho);
			      double func_val =  Qfact*LL.V_pipi(t, Lprime, my_fit_par.Mrho*Mp, my_fit_par.gpi*(my_fit_par.Mrho*Mp/fp),  my_fit_par.Mpi_corr*Mp, my_fit_par.kappa,  Elev_MpiL_4dot2)+ LL.Vdual(t, my_fit_par.Mrho*Mp, my_fit_par.Ed*Mp, my_fit_par.Rd);

			      double F_int_val = kern_val*func_val;
        

			      return F_int_val;
			    };

	auto F_infL = [&](double t) {

			//double kern_val = 4.0*pow(alpha,2)*sqrt(m_rho*a_distr.distr[ijack]/(my_fit_par.Mrho*Mp))*kernel_K(t, my_fit_par.Mrho*Mp/m_rho);
		  
			double kern_val = 4.0*pow(alpha,2)*w(t,Simps_ord)*kernel_K(t, a_distr.distr[ijack] );
			//double kern_val= 4.0*pow(alpha,2)*kernel_K(t, my_fit_par.Mrho*Mp/m_rho);

			double func_val =  Qfact*LL.V_pipi_infL(t, my_fit_par.Mrho*Mp, my_fit_par.gpi*(my_fit_par.Mrho*Mp/fp),  my_fit_par.Mpi_corr*Mp, my_fit_par.kappa)+ LL.Vdual(t, my_fit_par.Mrho*Mp, my_fit_par.Ed*Mp, my_fit_par.Rd);

			double F_int_val = kern_val*func_val;
        

			return F_int_val;

		      };

   

    
    
	//Vfloat Sum_Series_Terms;
	double agm2_W_FSEs_GS = 0.0;
	double agm2_summ=0.0;
	double agm2_2L_summ=0.0;
	double agm2_infL_summ=0.0;
	double agm2_Lprime_summ=0.0;
	int Nterms = 2000;
	int Nterms_2L= 2000;
	int Nterms_win= 100;

	for(int iterm=1;iterm<Nterms_win;iterm++) {

	  double F_n = FSEs_win_OS((double)iterm);
	  agm2_W_FSEs_GS += F_n;
	}
	agm2_W_FSEs_GS_distr.distr.push_back(agm2_W_FSEs_GS);
    
	//actual volume
	for(int iterm=1;iterm<Nterms;iterm++) {
	  double F_n = F_int((double)iterm);
	  //Sum_Series 
	  agm2_summ += F_n;}
	agm2_dual.distr.push_back(agm2_summ);
   

	//1.5L volume
	for(int iterm=1;iterm<Nterms_2L;iterm++) {
	  double F_n_2L = F_int_2L((double)iterm);
	  //Sum Series 
	  agm2_2L_summ += F_n_2L;}
	agm2_2L_dual.distr.push_back(agm2_2L_summ);


	//Mpi*L=4.2 volume
	for(int iterm=1;iterm<Nterms_2L;iterm++) {
	  double F_n_Lprime = F_int_Lprime((double)iterm);
	  //Sum Series 
	  agm2_Lprime_summ += F_n_Lprime;}
	agm2_Lprime_dual.distr.push_back(agm2_Lprime_summ);
    

	//inf volume
	for(int iterm=1; iterm<Nterms_2L;iterm++) {
	  double F_n_infL = F_infL((double)iterm);
	  //Sum series
	  agm2_infL_summ += F_n_infL;
	}
	agm2_infL_dual.distr.push_back(agm2_infL_summ);

      
    

      }

  
    

      //push_back
      Edual_OS.distr_list.push_back(Ed);
      Rdual_OS.distr_list.push_back(Rd);
      Mrho_OS.distr_list.push_back(Mr);
      gpi_OS.distr_list.push_back(g);
      Kappa_OS.distr_list.push_back(kap);
      M_corr_OS.distr_list.push_back(mpi_corr);
      agm2_light_fit_OS.distr_list.push_back(agm2_dual);
      agm2_light_2L_fit_OS.distr_list.push_back(agm2_2L_dual);
      agm2_light_infL_fit_OS.distr_list.push_back(agm2_infL_dual);
      agm2_light_Lprime_fit_OS.distr_list.push_back(agm2_Lprime_dual);



      //print ensemble infos

      //find whether the disconnected diagram corresponding to i_ens has been analyzed
      bool Found_disco_tag=false;
      int disco_ens=0;
      if(Include_light_disco) {
	for(int j=0; j<(signed)disco_light_Tags.size();j++) { if(disco_light_Tags[j] == V_light_1.Tag[i_ens]) {Found_disco_tag=true; disco_ens=j; break;}}
      }

      //intermediate and SD window info
      cout<<"%%%%%%% Print info from W and SD windows OS-correlator %%%%%%"<<endl;
      cout<<"SD: "<<agm2_light_SD_OS.ave(i_ens)<<" +- "<<agm2_light_SD_OS.err(i_ens)<<endl;
      cout<<"SD (ELM): "<<agm2_light_SD_ELM_OS.ave(i_ens)<<" +- "<<agm2_light_SD_ELM_OS.err(i_ens)<<endl;
      cout<<"W: "<<agm2_light_W_OS.ave(i_ens)<<" +- "<<agm2_light_W_OS.err(i_ens)<<endl;
      cout<<"W (ELM): "<<agm2_light_W_ELM_OS.ave(i_ens)<<" +- "<<agm2_light_W_ELM_OS.err(i_ens)<<endl;
      if(Found_disco_tag) {
	cout<<"%%%%%% disconnected contribution to the windows  %%%%%%"<<endl;
	cout<<"SD: "<<agm2_disco_light_SD.ave(disco_ens)<<" +- "<<agm2_disco_light_SD.err(disco_ens)<<endl;
	cout<<"SD (ELM): "<<agm2_disco_light_SD_ELM.ave(disco_ens)<<" +- "<<agm2_disco_light_SD_ELM.err(disco_ens)<<endl;
	cout<<"W: "<<agm2_disco_light_W.ave(disco_ens)<<" +- "<<agm2_disco_light_W.err(disco_ens)<<endl;
	cout<<"W (ELM): "<<agm2_disco_light_W_ELM.ave(disco_ens)<<" +- "<<agm2_disco_light_W_ELM.err(disco_ens)<<endl;
      }
      cout<<"%%%%%%% End windows info %%%%%%%"<<endl;
      cout<<"Edual/Mpi: "<<(Ed).ave()<<" +- "<<(Ed).err()<<endl;
      cout<<"Edual: "<<(Ed*mp).ave()<<" +- "<<(Ed*mp).err()<<endl;
      cout<<"Rdual: "<<Rd.ave()<<" +- "<<Rd.err()<<endl;
      cout<<"Mrho/Mpi: "<<(Mr).ave()<<" +- "<<(Mr).err()<<endl;
      cout<<"Mrho (lattice units): "<<(Mr*mp).ave()<<" +- "<<(Mr*mp).err()<<endl;
      cout<<"Mrho (GeV): "<<(Mr*mp/a_distr).ave()<<" +- "<<(Mr*mp/a_distr).err()<<endl;
      cout<<"Mrho/fp: "<<(Mr*mp/fp_fit[i_ens]).ave()<<" +- "<<(Mr*mp/fp_fit[i_ens]).err()<<endl;
      cout<<"gpi: "<<g.ave()<<" +- "<<g.err()<<endl;
      cout<<"kappa: "<<kap.ave()<<" +- "<<kap.err()<<endl;
      cout<<"Mpi corr: "<<mpi_corr.ave()<<" +- "<<mpi_corr.err()<<endl;
      cout<<"Mpi: "<<mp.ave()<<" +- "<<mp.err()<<endl;
      cout<<"Mpi*L : "<<mp.ave()*L_list[i_ens]<<" +- "<<mp.err()*L_list[i_ens]<<endl;
      cout<<"L: "<<L_list[i_ens]<<endl;
      cout<<"agm2 pp+dual (L): "<<agm2_dual.ave()<<" +- "<<agm2_dual.err()<<endl;
      cout<<"agm2 pp+dual (1.5*L): "<<agm2_2L_dual.ave()<<" +- "<<agm2_2L_dual.err()<<endl;
      cout<<"agm2 pp+dual (Mpi*L=4.2): "<<agm2_Lprime_dual.ave()<<" +- "<<agm2_Lprime_dual.err()<<endl;
      cout<<"agm2 pp+dual (L -> infty): "<<agm2_infL_dual.ave()<<" +- "<<agm2_infL_dual.err()<<endl;
      cout<<"%%%%% FSEs GS-prediction for light-connected contribution to the intermediate window %%%%%%"<<endl;
      cout<<"a_mu^W ( infty) - a_mu^W(L) ) : "<<agm2_W_FSEs_GS_distr.ave()<<" +- "<<agm2_W_FSEs_GS_distr.err()<<endl;
      cout<<"#######################################"<<endl;

    }

  }

  
  //############################################################################################################################################


  
  cout<<"########### DONE #############"<<endl;

  //print pseudoscalar masses for physical point ensembles
  Print_To_File(ens_masses_id, {Mp_light_tm.ave(), Mp_light_tm.err(), Mp_s1_tm.ave(), Mp_s1_tm.err(), Mp_s2_tm.ave(), Mp_s2_tm.err(), Mp_c1_tm.ave(), Mp_c1_tm.err(), Mp_c2_tm.ave(), Mp_c2_tm.err(), Mp_c3_tm.ave(), Mp_c3_tm.err()}, "../data/gm2/MP_masses_tm.dat", "", "#Ensemble   Mp_l   Mp_s1   Mp_s2   Mp_c1    Mp_c2   Mp_c3");
  
  Print_To_File(ens_masses_id, {Mp_light_OS.ave(), Mp_light_OS.err(), Mp_s1_OS.ave(), Mp_s1_OS.err(), Mp_s2_OS.ave(), Mp_s2_OS.err(), Mp_c1_OS.ave(), Mp_c1_OS.err(), Mp_c2_OS.ave(), Mp_c2_OS.err(), Mp_c3_OS.ave(), Mp_c3_OS.err()}, "../data/gm2/MP_masses_OS.dat", "", "#Ensemble   Mp_l   Mp_s1   Mp_s2   Mp_c1    Mp_c2   Mp_c3");
  

 
 
  //get Mpi corresponding to the strange run
  for(int tag_i=0; tag_i < (signed)V_strange_1_L.Tag.size(); tag_i ++) {
    bool find_light_tag=false;
    for(int tag_j=0; tag_j < (signed)V_light_1.Tag.size(); tag_j++) {
      if( V_strange_1_L.Tag[tag_i] == V_light_1.Tag[tag_j] ) {  find_light_tag=true; Mpi_fit_strange.distr_list.push_back( Mpi_fit.distr_list[tag_j]); break;}
    }
    if( !find_light_tag) {
      bool find_charm_tag=false;
      for(int tag_j=0; tag_j < (signed)V_charm_1_L.Tag.size(); tag_j++) {
	if( V_strange_1_L.Tag[tag_i] == V_charm_1_L.Tag[tag_j] ) { find_charm_tag = true; Mpi_fit_strange.distr_list.push_back( Mpi_fit_charm.distr_list[tag_j]); break;}
      }

      if(!find_charm_tag) crash("Cannot find strange ensemble: "+V_strange_1_L.Tag[tag_i]+" in either light or charm ensemble list, in order to get Mpi");
    }
  }

   

  //strange
  //print RCs
  Print_To_File(V_strange_1_L.Tag, {L_strange_list, a_strange_list, ms_extr_list.ave(), ms_extr_list.err(), Za_fit_strange_Extr.ave(), Za_fit_strange_Extr.err()} , "../data/gm2/strange/Z_"+Extrapolation_strange_mode+"/Za.list", "", "#Ens L a ms Za" );
  Print_To_File(V_strange_1_L.Tag, {L_strange_list, a_strange_list, ms_extr_list.ave(), ms_extr_list.err(), Zv_fit_strange_Extr.ave(), Zv_fit_strange_Extr.err()} , "../data/gm2/strange/Z_"+Extrapolation_strange_mode+"/Zv.list", "", "#Ens L a ms Zv" );

  //tm
  //L
  Print_To_File(V_strange_1_L.Tag, {L_strange_list, a_strange_list, ml_strange_list, Mpi_fit_strange.ave(), Mpi_fit_strange.err(),  ZV_fit_strange_L.ave(), ZV_fit_strange_L.err(), MV_fit_strange_L.ave(), MV_fit_strange_L.err(), agm2_strange_L.ave(), agm2_strange_L.err()}, "../data/gm2/strange/tm_"+Extrapolation_strange_mode+"/agm2_fit_ELM_L.list", "", "# Ens L a ml Mpi  ZV   MV    agm2(ELM)");
  //M
  Print_To_File(V_strange_1_M.Tag, {L_strange_list, a_strange_list, ml_strange_list, Mpi_fit_strange.ave(), Mpi_fit_strange.err(),   ZV_fit_strange_M.ave(), ZV_fit_strange_M.err(), MV_fit_strange_M.ave(), MV_fit_strange_M.err(), agm2_strange_M.ave(), agm2_strange_M.err()}, "../data/gm2/strange/tm_"+Extrapolation_strange_mode+"/agm2_fit_ELM_M.list", "", "# Ens L a ml Mpi  ZV   MV    agm2(ELM)");
  //Extr
  Print_To_File(V_strange_1_L.Tag, {L_strange_list, a_strange_list, ml_strange_list, Mpi_fit_strange.ave(), Mpi_fit_strange.err(),  agm2_strange_Extr.ave(), agm2_strange_Extr.err()}, "../data/gm2/strange/tm_"+Extrapolation_strange_mode+"/agm2_fit_ELM_Extr.list", "", "# Ens L a ml Mpi Mpi_OS agm2(ELM)");

  
  //OS
  //L
  Print_To_File(V_strange_OS_1_L.Tag, {L_strange_list, a_strange_list, ml_strange_list, Mpi_fit_strange.ave(), Mpi_fit_strange.err(),  ZV_fit_strange_OS_L.ave(), ZV_fit_strange_OS_L.err(), MV_fit_strange_OS_L.ave(), MV_fit_strange_OS_L.err(), agm2_strange_OS_L.ave(), agm2_strange_OS_L.err()}, "../data/gm2/strange/OS_"+Extrapolation_strange_mode+"/agm2_fit_ELM_L.list", "", "# Ens L a ml Mpi  ZV   MV    agm2(ELM)");
  //M
  Print_To_File(V_strange_OS_1_M.Tag, {L_strange_list, a_strange_list, ml_strange_list, Mpi_fit_strange.ave(), Mpi_fit_strange.err(),  ZV_fit_strange_OS_M.ave(), ZV_fit_strange_OS_M.err(), MV_fit_strange_OS_M.ave(), MV_fit_strange_OS_M.err(), agm2_strange_OS_M.ave(), agm2_strange_OS_M.err()}, "../data/gm2/strange/OS_"+Extrapolation_strange_mode+"/agm2_fit_ELM_M.list", "", "# Ens L a ml Mpi  ZV   MV    agm2(ELM)");
  //Extr
  Print_To_File(V_strange_OS_1_L.Tag, {L_strange_list, a_strange_list, ml_strange_list, Mpi_fit_strange.ave(), Mpi_fit_strange.err(), agm2_strange_OS_Extr.ave(), agm2_strange_OS_Extr.err()}, "../data/gm2/strange/OS_"+Extrapolation_strange_mode+"/agm2_fit_ELM_Extr.list", "", "# Ens L a ml Mpi agm2(ELM)");


  //strange non-ELM
  //tm
  //L
  Print_To_File(V_strange_1_L.Tag, {L_strange_list, a_strange_list, ml_strange_list, Mpi_fit_strange.ave(), Mpi_fit_strange.err(),  ZV_fit_strange_L.ave(), ZV_fit_strange_L.err(), MV_fit_strange_L.ave(), MV_fit_strange_L.err(), agm2_strange_No_ELM_L.ave(), agm2_strange_No_ELM_L.err()}, "../data/gm2/strange/tm_"+Extrapolation_strange_mode+"/agm2_fit_L.list", "", "# Ens L a ml Mpi ZV   MV  agm2");
  //M
  Print_To_File(V_strange_1_M.Tag, {L_strange_list, a_strange_list, ml_strange_list, Mpi_fit_strange.ave(), Mpi_fit_strange.err(),  ZV_fit_strange_M.ave(), ZV_fit_strange_M.err(), MV_fit_strange_M.ave(), MV_fit_strange_M.err(), agm2_strange_No_ELM_M.ave(), agm2_strange_No_ELM_M.err()}, "../data/gm2/strange/tm_"+Extrapolation_strange_mode+"/agm2_fit_M.list", "", "# Ens L a ml Mpi  ZV   MV agm2");
  //Extr
  Print_To_File(V_strange_1_L.Tag, {L_strange_list, a_strange_list, ml_strange_list, Mpi_fit_strange.ave(), Mpi_fit_strange.err(), agm2_strange_No_ELM_Extr.ave(), agm2_strange_No_ELM_Extr.err()}, "../data/gm2/strange/tm_"+Extrapolation_strange_mode+"/agm2_fit_Extr.list", "", "# Ens L a ml Mpi  agm2");
  
  //OS
  //L
  Print_To_File(V_strange_OS_1_L.Tag, {L_strange_list, a_strange_list, ml_strange_list, Mpi_fit_strange.ave(), Mpi_fit_strange.err(),  ZV_fit_strange_OS_L.ave(), ZV_fit_strange_OS_L.err(), MV_fit_strange_OS_L.ave(), MV_fit_strange_OS_L.err(), agm2_strange_OS_No_ELM_L.ave(), agm2_strange_OS_No_ELM_L.err()}, "../data/gm2/strange/OS_"+Extrapolation_strange_mode+"/agm2_fit_L.list", "", "# Ens L a ml Mpi ZV   MV    agm2");
  //M
  Print_To_File(V_strange_OS_1_M.Tag, {L_strange_list, a_strange_list, ml_strange_list, Mpi_fit_strange.ave(), Mpi_fit_strange.err(), ZV_fit_strange_OS_M.ave(), ZV_fit_strange_OS_M.err(), MV_fit_strange_OS_M.ave(), MV_fit_strange_OS_M.err(), agm2_strange_OS_No_ELM_M.ave(), agm2_strange_OS_No_ELM_M.err()}, "../data/gm2/strange/OS_"+Extrapolation_strange_mode+"/agm2_fit_M.list", "", "# Ens L a ml Mpi  ZV   MV    agm2");
  //Extr
  Print_To_File(V_strange_OS_1_L.Tag, {L_strange_list, a_strange_list, ml_strange_list, Mpi_fit_strange.ave(), Mpi_fit_strange.err(), agm2_strange_OS_No_ELM_Extr.ave(), agm2_strange_OS_No_ELM_Extr.err()}, "../data/gm2/strange/OS_"+Extrapolation_strange_mode+"/agm2_fit_Extr.list", "", "# Ens L a ml Mpi  agm2");


  
  //disco strange
  if(Include_strange_disco) {
    Print_To_File(disco_strange_Tags,{agm2_disco_strange_No_ELM.ave(), agm2_disco_strange_No_ELM.err()} , "../data/gm2/strange/disco/agm2_fit.list", "", "#ENS agm2");
    //improved
    Print_To_File(disco_impr_strange_Tags,{agm2_disco_impr_strange_No_ELM.ave(), agm2_disco_impr_strange_No_ELM.err()} , "../data/gm2/strange/disco/agm2_impr_fit.list", "", "#ENS agm2");
  }

     

  //charm
  //print RCs
  Print_To_File(V_charm_1_L.Tag, {L_charm_list, a_charm_list, mc_extr_list.ave(), mc_extr_list.err(), Za_fit_charm_Extr.ave(), Za_fit_charm_Extr.err(), Za_diff_charm.ave(), Za_diff_charm.err(), Za_diff_strange_charm.ave(), Za_diff_strange_charm.err()} , "../data/gm2/charm/Z_"+Extrapolation_charm_mode+"/Za.list", "", "#Ens L a mc Za Za_diff  Za_diff( w strange)" );
  Print_To_File(V_charm_1_L.Tag, {L_charm_list, a_charm_list, mc_extr_list.ave(), mc_extr_list.err(), Zv_fit_charm_Extr.ave(), Zv_fit_charm_Extr.err(), Zv_diff_charm.ave(), Zv_diff_charm.err(), Zv_diff_strange_charm.ave(), Zv_diff_strange_charm.err()} , "../data/gm2/charm/Z_"+Extrapolation_charm_mode+"/Zv.list", "", "#Ens L a mc Zv Zv_diff Zv_diff( w strange)" );

  //tm
  //L
  Print_To_File(V_charm_1_L.Tag, {L_charm_list, a_charm_list, ml_charm_list, Mpi_fit_charm.ave(), Mpi_fit_charm.err(), Mpi_OS_fit_charm.ave(), Mpi_OS_fit_charm.err(),  ZV_fit_charm_L.ave(), ZV_fit_charm_L.err(), MV_fit_charm_L.ave(), MV_fit_charm_L.err(), agm2_charm_L.ave(), agm2_charm_L.err()}, "../data/gm2/charm/tm_"+Extrapolation_charm_mode+"/agm2_fit_ELM_L.list", "", "# Ens  L a ml Mpi Mpi_OS ZV MV  agm2");
  //M
  Print_To_File(V_charm_1_M.Tag, {L_charm_list, a_charm_list, ml_charm_list, Mpi_fit_charm.ave(), Mpi_fit_charm.err(), Mpi_OS_fit_charm.ave(), Mpi_OS_fit_charm.err(),  ZV_fit_charm_M.ave(), ZV_fit_charm_M.err(), MV_fit_charm_M.ave(), MV_fit_charm_M.err(), agm2_charm_M.ave(), agm2_charm_M.err()}, "../data/gm2/charm/tm_"+Extrapolation_charm_mode+"/agm2_fit_ELM_M.list", "", "# Ens  L a ml Mpi Mpi_OS  ZV MV  agm2");
  //H
  Print_To_File(V_charm_1_H.Tag, {L_charm_list, a_charm_list, ml_charm_list, Mpi_fit_charm.ave(), Mpi_fit_charm.err(), Mpi_OS_fit_charm.ave(), Mpi_OS_fit_charm.err(),  ZV_fit_charm_H.ave(), ZV_fit_charm_H.err(), MV_fit_charm_H.ave(), MV_fit_charm_H.err(), agm2_charm_H.ave(), agm2_charm_H.err()}, "../data/gm2/charm/tm_"+Extrapolation_charm_mode+"/agm2_fit_ELM_H.list", "", "# Ens  L a ml Mpi Mpi_OS ZV MV  agm2");
  //Extr
  Print_To_File(V_charm_1_L.Tag, {L_charm_list, a_charm_list, ml_charm_list, Mpi_fit_charm.ave(), Mpi_fit_charm.err(), Mpi_OS_fit_charm.ave(), Mpi_OS_fit_charm.err(), agm2_charm_Extr.ave(), agm2_charm_Extr.err()}, "../data/gm2/charm/tm_"+Extrapolation_charm_mode+"/agm2_fit_ELM_Extr.list", "", "# Ens L a ml Mpi Mpi_OS agm2");

  //OS
  //L
  Print_To_File(V_charm_OS_1_L.Tag, {L_charm_list, a_charm_list, ml_charm_list, Mpi_fit_charm.ave(), Mpi_fit_charm.err(), Mpi_OS_fit_charm.ave(), Mpi_OS_fit_charm.err(), ZV_fit_charm_OS_L.ave(), ZV_fit_charm_OS_L.err(), MV_fit_charm_OS_L.ave(), MV_fit_charm_OS_L.err(), agm2_charm_OS_L.ave(), agm2_charm_OS_L.err()}, "../data/gm2/charm/OS_"+Extrapolation_charm_mode+"/agm2_fit_ELM_L.list", "", "# Ens  L a ml Mpi Mpi_OS  ZV MV  agm2");
  //M
  Print_To_File(V_charm_OS_1_M.Tag, {L_charm_list, a_charm_list, ml_charm_list, Mpi_fit_charm.ave(), Mpi_fit_charm.err(), Mpi_OS_fit_charm.ave(), Mpi_OS_fit_charm.err(), ZV_fit_charm_OS_M.ave(), ZV_fit_charm_OS_M.err(), MV_fit_charm_OS_M.ave(), MV_fit_charm_OS_M.err(), agm2_charm_OS_M.ave(), agm2_charm_OS_M.err()}, "../data/gm2/charm/OS_"+Extrapolation_charm_mode+"/agm2_fit_ELM_M.list", "", "# Ens  L a ml Mpi Mpi_OS  ZV MV  agm2");
  //H
  Print_To_File(V_charm_OS_1_H.Tag, {L_charm_list, a_charm_list, ml_charm_list, Mpi_fit_charm.ave(), Mpi_fit_charm.err(), Mpi_OS_fit_charm.ave(), Mpi_OS_fit_charm.err(), ZV_fit_charm_OS_H.ave(), ZV_fit_charm_OS_H.err(), MV_fit_charm_OS_H.ave(), MV_fit_charm_OS_H.err(), agm2_charm_OS_H.ave(), agm2_charm_OS_H.err()}, "../data/gm2/charm/OS_"+Extrapolation_charm_mode+"/agm2_fit_ELM_H.list", "", "# Ens  L a ml Mpi Mpi_OS  ZV MV  agm2");
  //Extr
  Print_To_File(V_charm_OS_1_L.Tag, {L_charm_list, a_charm_list, ml_charm_list, Mpi_fit_charm.ave(), Mpi_fit_charm.err(), Mpi_OS_fit_charm.ave(), Mpi_OS_fit_charm.err(), agm2_charm_OS_Extr.ave(), agm2_charm_OS_Extr.err()}, "../data/gm2/charm/OS_"+Extrapolation_charm_mode+"/agm2_fit_ELM_Extr.list", "", "# Ens L a ml Mpi Mpi_OS agm2");


  //charm NON_ELM
  Print_To_File(V_charm_1_L.Tag, {L_charm_list, a_charm_list, ml_charm_list, Mpi_fit_charm.ave(), Mpi_fit_charm.err(), Mpi_OS_fit_charm.ave(), Mpi_OS_fit_charm.err(),  ZV_fit_charm_L.ave(), ZV_fit_charm_L.err(), MV_fit_charm_L.ave(), MV_fit_charm_L.err(), agm2_charm_No_ELM_L.ave(), agm2_charm_No_ELM_L.err()}, "../data/gm2/charm/tm_"+Extrapolation_charm_mode+"/agm2_fit_L.list", "", "# Ens  L a ml Mpi Mpi_OS  ZV MV  agm2");
  //M
  Print_To_File(V_charm_1_M.Tag, {L_charm_list, a_charm_list, ml_charm_list, Mpi_fit_charm.ave(), Mpi_fit_charm.err(), Mpi_OS_fit_charm.ave(), Mpi_OS_fit_charm.err(),  ZV_fit_charm_M.ave(), ZV_fit_charm_M.err(), MV_fit_charm_M.ave(), MV_fit_charm_M.err(), agm2_charm_No_ELM_M.ave(), agm2_charm_No_ELM_M.err()}, "../data/gm2/charm/tm_"+Extrapolation_charm_mode+"/agm2_fit_M.list", "", "# Ens  L a ml Mpi Mpi_OS ZV MV  agm2");
  //H
  Print_To_File(V_charm_1_H.Tag, {L_charm_list, a_charm_list, ml_charm_list, Mpi_fit_charm.ave(), Mpi_fit_charm.err(), Mpi_OS_fit_charm.ave(), Mpi_OS_fit_charm.err(),  ZV_fit_charm_H.ave(), ZV_fit_charm_H.err(), MV_fit_charm_H.ave(), MV_fit_charm_H.err(), agm2_charm_No_ELM_H.ave(), agm2_charm_No_ELM_H.err()}, "../data/gm2/charm/tm_"+Extrapolation_charm_mode+"/agm2_fit_H.list", "", "# Ens  L a ml Mpi Mpi_OS  ZV MV  agm2");
  //Extr
  Print_To_File(V_charm_1_L.Tag, {L_charm_list, a_charm_list, ml_charm_list, Mpi_fit_charm.ave(), Mpi_fit_charm.err(), Mpi_OS_fit_charm.ave(), Mpi_OS_fit_charm.err(), agm2_charm_No_ELM_Extr.ave(), agm2_charm_No_ELM_Extr.err()}, "../data/gm2/charm/tm_"+Extrapolation_charm_mode+"/agm2_fit_Extr.list", "", "# Ens L a ml Mpi Mpi_OS agm2");

  //OS
  //L
  Print_To_File(V_charm_OS_1_L.Tag, {L_charm_list, a_charm_list, ml_charm_list, Mpi_fit_charm.ave(), Mpi_fit_charm.err(), Mpi_OS_fit_charm.ave(), Mpi_OS_fit_charm.err(), ZV_fit_charm_OS_L.ave(), ZV_fit_charm_OS_L.err(), MV_fit_charm_OS_L.ave(), MV_fit_charm_OS_L.err(), agm2_charm_OS_No_ELM_L.ave(), agm2_charm_OS_No_ELM_L.err()}, "../data/gm2/charm/OS_"+Extrapolation_charm_mode+"/agm2_fit_L.list", "", "# Ens  L a ml Mpi Mpi_OS ZV MV  agm2");
  //M
  Print_To_File(V_charm_OS_1_M.Tag, {L_charm_list, a_charm_list, ml_charm_list, Mpi_fit_charm.ave(), Mpi_fit_charm.err(), Mpi_OS_fit_charm.ave(), Mpi_OS_fit_charm.err(), ZV_fit_charm_OS_M.ave(), ZV_fit_charm_OS_M.err(), MV_fit_charm_OS_M.ave(), MV_fit_charm_OS_M.err(), agm2_charm_OS_No_ELM_M.ave(), agm2_charm_OS_No_ELM_M.err()}, "../data/gm2/charm/OS_"+Extrapolation_charm_mode+"/agm2_fit_M.list", "", "# Ens  L a ml Mpi Mpi_OS ZV MV  agm2");
  //H
  Print_To_File(V_charm_OS_1_H.Tag, {L_charm_list, a_charm_list, ml_charm_list, Mpi_fit_charm.ave(), Mpi_fit_charm.err(), Mpi_OS_fit_charm.ave(), Mpi_OS_fit_charm.err(), ZV_fit_charm_OS_H.ave(), ZV_fit_charm_OS_H.err(), MV_fit_charm_OS_H.ave(), MV_fit_charm_OS_H.err(), agm2_charm_OS_No_ELM_H.ave(), agm2_charm_OS_No_ELM_H.err()}, "../data/gm2/charm/OS_"+Extrapolation_charm_mode+"/agm2_fit_H.list", "", "# Ens  L a ml Mpi Mpi_OS ZV MV  agm2");
  //Extr
  Print_To_File(V_charm_OS_1_L.Tag, {L_charm_list, a_charm_list, ml_charm_list, Mpi_fit_charm.ave(), Mpi_fit_charm.err(), Mpi_OS_fit_charm.ave(), Mpi_OS_fit_charm.err(), agm2_charm_OS_No_ELM_Extr.ave(), agm2_charm_OS_No_ELM_Extr.err()}, "../data/gm2/charm/OS_"+Extrapolation_charm_mode+"/agm2_fit_Extr.list", "", "# Ens L a ml Mpi Mpi_OS agm2");

    
  //disco charm
  if(Include_charm_disco) {
    Print_To_File(disco_charm_Tags,{agm2_disco_charm_No_ELM.ave(), agm2_disco_charm_No_ELM.err()} , "../data/gm2/charm/disco/agm2_fit.list", "", "#ENS agm2");
    //improved
    Print_To_File(disco_impr_charm_Tags,{agm2_disco_impr_charm_No_ELM.ave(), agm2_disco_impr_charm_No_ELM.err()} , "../data/gm2/charm/disco/agm2_impr_fit.list", "", "#ENS agm2");
  }

  
  //light
  //tm
  Print_To_File(V_light_1.Tag, {ZV_fit_light.ave(), ZV_fit_light.err(), MV_fit_light.ave(), MV_fit_light.err()}, "../data/gm2/light/tm/ZV_MV_fitted.list", "", "#Ens ZV MV");
  //OS
  Print_To_File(V_light_1.Tag, {ZV_fit_light_OS.ave(), ZV_fit_light_OS.err(), MV_fit_light_OS.ave(), MV_fit_light_OS.err()}, "../data/gm2/light/OS/ZV_MV_fitted.list", "", "#Ens ZV MV");

  //disco light
  if(Include_light_disco) {
    Print_To_File(disco_light_Tags,{L_list_disco, a_list_disco, ml_list_disco, Mpi_fit_disco.ave(), Mpi_fit_disco.err(), Mpi_OS_fit_disco.ave(), Mpi_OS_fit_disco.err(), fp_fit_disco.ave(), fp_fit_disco.err(), Zv_fit_disco.ave(), Zv_fit_disco.err(), agm2_disco_light.ave(), agm2_disco_light.err(), agm2_disco_light_ELM.ave(), agm2_disco_light_ELM.err()} , "../data/gm2/light/disco/agm2_fit.list", "", "#ENS L a ml Mpi_tm Mpi_OS fp Zv  agm2");
  }

  //off-diagonal disconnected
  if(Include_off_diagonal_disco) {
    
    //light -strange
    Print_To_File(disco_impr_light_strange_Tags, {agm2_disco_impr_light_strange_No_ELM.ave(), agm2_disco_impr_light_strange_No_ELM.err()} , "../data/gm2/light_strange/disco/agm2_impr_fit.list", "", "");
    Print_To_File(disco_impr_lightD_strange_Tags, {agm2_disco_impr_lightD_strange_No_ELM.ave(), agm2_disco_impr_lightD_strange_No_ELM.err()} , "../data/gm2/light_strange/disco/agm2_impr_D_fit.list", "", "");

    //light-charm
    Print_To_File(disco_impr_light_charm_Tags, {agm2_disco_impr_light_charm_No_ELM.ave(), agm2_disco_impr_light_charm_No_ELM.err()} , "../data/gm2/light_charm/disco/agm2_impr_fit.list", "", "");
    Print_To_File(disco_impr_lightD_charm_Tags, {agm2_disco_impr_lightD_charm_No_ELM.ave(), agm2_disco_impr_lightD_charm_No_ELM.err()} , "../data/gm2/light_charm/disco/agm2_impr_D_fit.list", "", "");

    //strange-charm
    Print_To_File(disco_impr_strange_charm_Tags, {agm2_disco_impr_strange_charm_No_ELM.ave(), agm2_disco_impr_strange_charm_No_ELM.err()} , "../data/gm2/strange_charm/disco/agm2_impr_fit.list", "", "");
    
  }



  //print informations on the windows
  //light
  //tm
  Print_To_File(V_light_1.Tag, {L_list, a_list, ml_list, Mpi_fit.ave(), Mpi_fit.err(), Mpi_OS_fit.ave(), Mpi_OS_fit.err(),  fp_fit.ave(), fp_fit.err(), Zv_fit.ave(), Zv_fit.err(), Za_fit.ave(), Za_fit.err(), Zp_ov_Zs_fit.ave(), Zp_ov_Zs_fit.err(), agm2_light_W.ave(), agm2_light_W.err(), agm2_light_W_ELM.ave(), agm2_light_W_ELM.err(), agm2_light_SD.ave(), agm2_light_SD.err(), agm2_light_SD_ELM.ave(), agm2_light_SD_ELM.err()}, "../data/gm2/light/tm/windows.list", "", "#ENS L a ml  Mpi_tm  Mpi_OS fp  Zv   Za   Zp/Zs   W    W(ELM)     SD     SD(ELM) ");
  //tm test 2a
  Print_To_File(V_light_1.Tag, {agm2_light_W_2a.ave(), agm2_light_W_2a.err()}, "../data/gm2/light/tm/windows_test_2a.list", "", "#ENS  W");
  Print_To_File(V_light_1.Tag, {agm2_light_W_der_err_scale_setting.ave(), agm2_light_W_der_err_scale_setting.err(), agm2_light_SD_der_err_scale_setting.ave(), agm2_light_SD_der_err_scale_setting.err(), agm2_light_tot_der_err_scale_setting.ave(), agm2_light_tot_der_err_scale_setting.err() }, "../data/gm2/light/tm/win_scale_setting_err.list", "", "#ENS  W   SD  LD");
  //OS
  Print_To_File(V_light_1.Tag, {L_list, a_list, ml_list, Mpi_fit.ave(), Mpi_fit.err(), Mpi_OS_fit.ave(), Mpi_OS_fit.err(),  fp_fit.ave(), fp_fit.err(), Zv_fit.ave(), Zv_fit.err(), Za_fit.ave(), Za_fit.err(), Zp_ov_Zs_fit.ave(), Zp_ov_Zs_fit.err(), agm2_light_W_OS.ave(), agm2_light_W_OS.err(), agm2_light_W_ELM_OS.ave(), agm2_light_W_ELM_OS.err(), agm2_light_SD_OS.ave(), agm2_light_SD_OS.err(), agm2_light_SD_ELM_OS.ave(), agm2_light_SD_ELM_OS.err()}, "../data/gm2/light/OS/windows.list", "", "#ENS L a ml  Mpi_tm  Mpi_OS fp  Zv   Za   Zp/Zs   W    W(ELM)     SD     SD(ELM) ");
  //result from bounding
  Print_To_File(V_light_1.Tag, {amu_W_bounding_tm_list.ave(), amu_W_bounding_tm_list.err()}, "../data/gm2/light/tm/windows_bounding.list", "", "#ENS W ");
  Print_To_File(V_light_1.Tag, {amu_W_bounding_OS_list.ave(), amu_W_bounding_OS_list.err()}, "../data/gm2/light/OS/windows_bounding.list", "", "#ENS W ");










  //##########################################################################################################################################################
  //compute physical point estimate for light window (tm and OS)
  distr_t_list amu_W_physical_point_tm(UseJack), amu_W_physical_point_OS(UseJack);
  distr_t_list amu_SD_physical_point_tm(UseJack), amu_SD_physical_point_OS(UseJack);
  distr_t_list corr_W_physical_point_tm(UseJack), corr_W_physical_point_OS(UseJack);
  distr_t_list amu_W_physical_point_tm_red(UseJack), amu_W_physical_point_OS_red(UseJack);
  distr_t_list amu_SD_physical_point_tm_red(UseJack), amu_SD_physical_point_OS_red(UseJack);
  vector<distr_t_list> amu_W_eps_physical_point_tm(eps_win_size);
  vector<distr_t_list> amu_W_eps_physical_point_OS(eps_win_size);
  vector<distr_t_list> amu_W_eps_physical_point_tm_red(eps_win_size);
  vector<distr_t_list> amu_W_eps_physical_point_OS_red(eps_win_size);
  vector<distr_t_list> amu_SD_tmins_physical_point_tm_red(tmins.size());
  vector<distr_t_list> amu_SD_tmins_physical_point_OS_red(tmins.size());
  distr_t_list a_distr_list_fixed_L(UseJack);
  distr_t_list Mpi_fit_fixed_L(UseJack);
  distr_t_list fp_fit_fixed_L(UseJack);
  Vfloat L_list_fixed_L;
  vector<string> V_light_Tag_fixed_L;
  for(int iens=0;iens<(signed)V_light_1.Tag.size();iens++) {

    cout<<"Computing physical point val of amu^W(light) and amu^SD(light) for ensemble: "<<V_light_1.Tag[iens]<<endl;

    distr_t amu_tm_ph, amu_OS_ph, amu_tm_ph_red, amu_OS_ph_red;
    distr_t amu_SD_tm_ph, amu_SD_OS_ph, amu_SD_tm_ph_red, amu_SD_OS_ph_red;
    vector<distr_t> amu_eps_tm_ph;
    vector<distr_t> amu_eps_OS_ph;
    vector<distr_t> amu_eps_tm_ph_red;
    vector<distr_t> amu_eps_OS_ph_red;
    for(int id_eps=0; id_eps<eps_win_size;id_eps++) { amu_eps_tm_ph_red.emplace_back(UseJack); amu_eps_OS_ph_red.emplace_back(UseJack);}
    int id_sea=-1;
    int id_val=-1;

    if(V_light_1.Tag[iens] != "cB211b.072.96") {

      //find sea_id
      for(int isea=0; isea<(signed)amu_W_sea_Ens_tag.size();isea++) {
	if(amu_W_sea_Ens_tag[isea] == V_light_1.Tag[iens]) id_sea=isea;
      }
      if(id_sea==-1) crash("Cannot find sea ensemble: "+V_light_1.Tag[iens]+" when computing physical point value of amu_W");

      //find val_id
      for(int ival=0; ival<(signed)amu_W_val_Ens_tag.size();ival++) {
	if(amu_W_val_Ens_tag[ival] == V_light_1.Tag[iens]) id_val=ival;
      }
      if(id_val==-1) crash("Cannot find val ensemble: "+V_light_1.Tag[iens]+" when computing physical point value of amu_W");

      //W
      amu_tm_ph= agm2_light_W.distr_list[iens] + amu_W_mass_corr_sea_tm_list.distr_list[id_sea] + amu_W_mass_corr_tm_list.distr_list[id_val];
      amu_OS_ph= agm2_light_W_OS.distr_list[iens] + amu_W_mass_corr_sea_OS_list.distr_list[id_sea] + amu_W_mass_corr_OS_list.distr_list[id_val];
      //SD
      amu_SD_tm_ph= agm2_light_SD.distr_list[iens] + amu_SD_mass_corr_sea_tm_list.distr_list[id_sea] + amu_SD_mass_corr_tm_list.distr_list[id_val];
      amu_SD_OS_ph= agm2_light_SD_OS.distr_list[iens] + amu_SD_mass_corr_sea_OS_list.distr_list[id_sea] + amu_SD_mass_corr_OS_list.distr_list[id_val];

      //correction
      corr_W_physical_point_tm.distr_list.push_back(amu_W_mass_corr_sea_tm_list.distr_list[id_sea] + amu_W_mass_corr_tm_list.distr_list[id_val]);
      corr_W_physical_point_OS.distr_list.push_back(amu_W_mass_corr_sea_OS_list.distr_list[id_sea] + amu_W_mass_corr_OS_list.distr_list[id_val]);

      
      for(int id_eps=0;id_eps<eps_win_size;id_eps++) {
	amu_eps_tm_ph.push_back( amu_W_eps_list_tm[id_eps].distr_list[iens] + amu_W_eps_sea_list_tm[id_eps].distr_list[id_sea] + amu_W_eps_val_list_tm[id_eps].distr_list[id_val]);
	amu_eps_OS_ph.push_back( amu_W_eps_list_OS[id_eps].distr_list[iens] + amu_W_eps_sea_list_OS[id_eps].distr_list[id_sea] + amu_W_eps_val_list_OS[id_eps].distr_list[id_val]);
      }

 
    }
    else if(V_light_1.Tag[iens] == "cB211b.072.96") {

      //find sea_id
      for(int isea=0; isea<(signed)amu_W_sea_Ens_tag.size();isea++) {
	if(amu_W_sea_Ens_tag[isea] == "cB211b.072.64") id_sea=isea;
      }
      if(id_sea==-1) crash("Cannot find sea ensemble: cB211b.072.64  when computing physical point value of amu_W");

      //find val_id
      for(int ival=0; ival<(signed)amu_W_val_Ens_tag.size();ival++) {
	if(amu_W_val_Ens_tag[ival] == V_light_1.Tag[iens]) id_val=ival;
      }
      if(id_val==-1) crash("Cannot find val ensemble: cB211b.072.64 when computing physical point value of amu_W");

      //W
      amu_tm_ph= agm2_light_W.distr_list[iens] + amu_W_mass_corr_sea_tm_list.distr_list[id_sea] + amu_W_mass_corr_tm_list.distr_list[id_val];
      amu_OS_ph= agm2_light_W_OS.distr_list[iens] + amu_W_mass_corr_sea_OS_list.distr_list[id_sea] + amu_W_mass_corr_OS_list.distr_list[id_val];
      //SD
      amu_SD_tm_ph= agm2_light_SD.distr_list[iens] + amu_SD_mass_corr_sea_tm_list.distr_list[id_sea] + amu_SD_mass_corr_tm_list.distr_list[id_val];
      amu_SD_OS_ph= agm2_light_SD_OS.distr_list[iens] + amu_SD_mass_corr_sea_OS_list.distr_list[id_sea] + amu_SD_mass_corr_OS_list.distr_list[id_val];

      //corrections
      corr_W_physical_point_tm.distr_list.push_back(amu_W_mass_corr_sea_tm_list.distr_list[id_sea] + amu_W_mass_corr_tm_list.distr_list[id_val]);
      corr_W_physical_point_OS.distr_list.push_back(amu_W_mass_corr_sea_OS_list.distr_list[id_sea] + amu_W_mass_corr_OS_list.distr_list[id_val]);
      
      
      for(int id_eps=0;id_eps<eps_win_size;id_eps++) {
	amu_eps_tm_ph.push_back( amu_W_eps_list_tm[id_eps].distr_list[iens] + amu_W_eps_sea_list_tm[id_eps].distr_list[id_sea] + amu_W_eps_val_list_tm[id_eps].distr_list[id_val]);
	amu_eps_OS_ph.push_back( amu_W_eps_list_OS[id_eps].distr_list[iens] + amu_W_eps_sea_list_OS[id_eps].distr_list[id_sea] + amu_W_eps_val_list_OS[id_eps].distr_list[id_val]);
      }


    }
    /*
    else if(V_light_1.Tag[iens].substr(1,1)=="C") {
      
      //find sea_id
      for(int isea=0; isea<(signed)amu_W_sea_Ens_tag.size();isea++) {
	if(amu_W_sea_Ens_tag[isea] == V_light_1.Tag[iens]) id_sea=isea;
      }
      if(id_sea==-1) crash("Cannot find sea ensemble: "+V_light_1.Tag[iens]+" when computing physical point value of amu_W");


      int id_val_ens_D96=-1;
      int id_val_ens_B64= -1;
      int id_mass_ens_D96=-1;
      int id_mass_ens_B64=-1;

      for(int i_mass=0; i_mass<(signed)V_light_1.Tag.size();i_mass++) {
	if( V_light_1.Tag[i_mass] == "cB211b.072.64") id_mass_ens_B64= i_mass;
	if( V_light_1.Tag[i_mass] == "cD211a.054.96") id_mass_ens_D96= i_mass;
	
      }
      if(id_mass_ens_D96==-1 || id_mass_ens_B64== -1) crash("Cannot find ensemble D96 or B64 when evaluating physical point value of amu_W for ensemble: "+V_light_1.Tag[iens]);

      for(int ival_ens=0; ival_ens<(signed)amu_W_val_Ens_tag.size(); ival_ens++) {
	if(amu_W_val_Ens_tag[ival_ens] == "cB211b.072.64") id_val_ens_B64=ival_ens;
	if(amu_W_val_Ens_tag[ival_ens] == "cD211a.054.96") id_val_ens_D96=ival_ens;
      }
      if(id_val_ens_D96 == -1 || id_val_ens_B64 == -1) crash("Cannot find ensemble D96 or B64 when evaluating physical point value of amu_W for ensemble: "+V_light_1.Tag[iens]);

      distr_t Mpi_B= Mpi_fit.distr_list[id_mass_ens_B64];
      distr_t Mpi_D= Mpi_fit.distr_list[id_mass_ens_D96];
      distr_t Mpi_C= Mpi_fit.distr_list[iens];

      distr_t corr_fact_B64=  (( Mpi_C*Mpi_C - 0.135*0.135*a_C*a_C)/(Mpi_B*Mpi_B - 0.135*0.135*a_B*a_B))*(a_C*a_C - a_D*a_D)/(a_B*a_B - a_D*a_D);
      distr_t corr_fact_D96=  (( Mpi_C*Mpi_C - 0.135*0.135*a_C*a_C)/(Mpi_D*Mpi_D - 0.135*0.135*a_D*a_D))*(a_B*a_B - a_C*a_C)/(a_B*a_B - a_D*a_D);


      amu_tm_ph= agm2_light_W.distr_list[iens] + amu_W_mass_corr_sea_tm_list.distr_list[id_sea] + corr_fact_B64*amu_W_mass_corr_tm_list.distr_list[id_val_ens_B64] + corr_fact_D96*amu_W_mass_corr_tm_list.distr_list[id_val_ens_D96];
      amu_OS_ph= agm2_light_W_OS.distr_list[iens] + amu_W_mass_corr_sea_OS_list.distr_list[id_sea] + corr_fact_B64*amu_W_mass_corr_OS_list.distr_list[id_val_ens_B64] + corr_fact_D96*amu_W_mass_corr_OS_list.distr_list[id_val_ens_D96];


      amu_SD_tm_ph= agm2_light_SD.distr_list[iens] + amu_SD_mass_corr_sea_tm_list.distr_list[id_sea] + corr_fact_B64*amu_SD_mass_corr_tm_list.distr_list[id_val_ens_B64] + corr_fact_D96*amu_SD_mass_corr_tm_list.distr_list[id_val_ens_D96];
      amu_SD_OS_ph= agm2_light_SD_OS.distr_list[iens] + amu_SD_mass_corr_sea_OS_list.distr_list[id_sea] + corr_fact_B64*amu_SD_mass_corr_OS_list.distr_list[id_val_ens_B64] + corr_fact_D96*amu_SD_mass_corr_OS_list.distr_list[id_val_ens_D96];

      
      
      for(int id_eps=0;id_eps<eps_win_size;id_eps++) {
	amu_eps_tm_ph.push_back( amu_W_eps_list_tm[id_eps].distr_list[iens] + amu_W_eps_sea_list_tm[id_eps].distr_list[id_sea] + corr_fact_B64*amu_W_eps_val_list_tm[id_eps].distr_list[id_val_ens_B64]+ corr_fact_D96*amu_W_eps_val_list_tm[id_eps].distr_list[id_val_ens_D96] );
	amu_eps_OS_ph.push_back( amu_W_eps_list_OS[id_eps].distr_list[iens] + amu_W_eps_sea_list_OS[id_eps].distr_list[id_sea] + corr_fact_B64*amu_W_eps_val_list_OS[id_eps].distr_list[id_val_ens_B64]+ corr_fact_D96*amu_W_eps_val_list_OS[id_eps].distr_list[id_val_ens_D96] );
      }

      
    } */
    else crash("When computing physical point value of amu_W, ensemble: "+V_light_1.Tag[iens]+" not recognized");

    //W
    amu_W_physical_point_tm.distr_list.push_back(amu_tm_ph);
    amu_W_physical_point_OS.distr_list.push_back(amu_OS_ph);
    //SD
    amu_SD_physical_point_tm.distr_list.push_back(amu_SD_tm_ph);
    amu_SD_physical_point_OS.distr_list.push_back(amu_SD_OS_ph);

    
    for(int id_eps=0; id_eps<eps_win_size;id_eps++) {
      amu_W_eps_physical_point_tm[id_eps].distr_list.push_back( amu_eps_tm_ph[id_eps]);
      amu_W_eps_physical_point_OS[id_eps].distr_list.push_back( amu_eps_OS_ph[id_eps]);
    }



    //################ COMPUTE VALUES TO PERFORM EXTRAPOLATION AT A FIXED VOLUME OF L^3 WITH L=5.46fm

    if(V_light_1.Tag[iens] != "cB211b.072.96") {

      //Interpolate B64 to L~ 5.46 fm using B96
      if(V_light_1.Tag[iens] == "cB211b.072.64") {
	int b96_id=0;
	int b96_id_val=0;
	bool Find_b96=false;
	bool Find_b96_val=false;
	for(int jens=0; jens<V_light_1.Tag.size(); jens++) {
	  if( V_light_1.Tag[jens] == "cB211b.072.96") { b96_id=jens; Find_b96=true; break;}
	}
	if(!Find_b96) crash("Cannot find ensemble cB211b.072.96 in W_red construction");

	for(int jens=0; jens<V_light_1.Tag.size(); jens++) {
	  if( V_light_1.Tag[jens] == "cB211b.072.96") { b96_id_val=jens; Find_b96_val=true; break;}
	}
	if(!Find_b96_val) crash("Cannot find ensemble cB211b.072.96 in W_red val construction");
	

	//compute B96 ph values
	//W
	distr_t amu_B96_tm_ph= agm2_light_W.distr_list[b96_id] + amu_W_mass_corr_sea_tm_list.distr_list[id_sea] + amu_W_mass_corr_tm_list.distr_list[b96_id_val];  
	distr_t amu_B96_OS_ph= agm2_light_W_OS.distr_list[b96_id] + amu_W_mass_corr_sea_OS_list.distr_list[id_sea] + amu_W_mass_corr_OS_list.distr_list[b96_id_val];
	//SD
	distr_t amu_SD_B96_tm_ph= agm2_light_SD.distr_list[b96_id] + amu_SD_mass_corr_sea_tm_list.distr_list[id_sea] + amu_SD_mass_corr_tm_list.distr_list[b96_id_val];  
	distr_t amu_SD_B96_OS_ph= agm2_light_SD_OS.distr_list[b96_id] + amu_SD_mass_corr_sea_OS_list.distr_list[id_sea] + amu_SD_mass_corr_OS_list.distr_list[b96_id_val];

	
	vector<distr_t> amu_B96_eps_tm_ph(eps_win_size);
	vector<distr_t> amu_B96_eps_OS_ph(eps_win_size);
	for(int id_eps=0;id_eps<eps_win_size;id_eps++) {
	  amu_B96_eps_tm_ph[id_eps] = amu_W_eps_list_tm[id_eps].distr_list[b96_id] + amu_W_eps_sea_list_tm[id_eps].distr_list[id_sea] + amu_W_eps_val_list_tm[id_eps].distr_list[b96_id_val] ;
	  amu_B96_eps_OS_ph[id_eps] = amu_W_eps_list_OS[id_eps].distr_list[b96_id] + amu_W_eps_sea_list_OS[id_eps].distr_list[id_sea] + amu_W_eps_val_list_OS[id_eps].distr_list[b96_id_val]  ; 
	}

	vector<distr_t> amu_SD_B96_tmins_tm_ph(tmins.size());
	vector<distr_t> amu_SD_B96_tmins_OS_ph(tmins.size());
	for(int tmins_id=0; tmins_id<(signed)tmins.size();tmins_id++) {
	  amu_SD_B96_tmins_tm_ph[tmins_id] = agm2_light_SD_tmins_distr_list[tmins_id].distr_list[b96_id];
	  amu_SD_B96_tmins_OS_ph[tmins_id] = agm2_light_SD_OS_tmins_distr_list[tmins_id].distr_list[b96_id];
	}

	vector<distr_t> amu_SD_Bextr_tmins_tm_ph(tmins.size());
	vector<distr_t> amu_SD_Bextr_tmins_OS_ph(tmins.size());


	distr_t Mpi_phys_fake_distr(UseJack);
	for(int ijack=0;ijack<Njacks;ijack++) Mpi_phys_fake_distr.distr.push_back( (0.135+GM()*0.0002/sqrt(Njacks-1.0))*a_distr_list.distr_list[iens].distr[ijack]);

	for(int ijck=0; ijck< Njacks;ijck++) {
	  //double exp_ML_B64= exp( -1.0*Mpi_fit.distr_list[iens].distr[ijck]*L_list[iens]);
	  //double exp_ML_B96= exp(-1.0*Mpi_fit.distr_list[b96_id].distr[ijck]*L_list[b96_id]);
	  //double exp_ML_extr= exp(-1.0*Mpi_fit.distr_list[iens].distr[ijck]*5.46*fm_to_inv_Gev/a_distr_list.distr_list[iens].distr[ijck]);

	  double exp_ML_B64= exp( -1.0*Mpi_phys_fake_distr.distr[ijck]*L_list[iens]);
	  double exp_ML_B96= exp(-1.0*Mpi_phys_fake_distr.distr[ijck]*L_list[b96_id]);
	  double exp_ML_extr= exp(-1.0*Mpi_phys_fake_distr.distr[ijck]*5.46*fm_to_inv_Gev/a_distr_list.distr_list[iens].distr[ijck]);

	  //tm W
	  double amu_96_tm= amu_B96_tm_ph.distr[ijck];
	  double amu_64_tm = amu_tm_ph.distr[ijck];
	  double amu_Linf_tm = (amu_96_tm*exp_ML_B64 - amu_64_tm*exp_ML_B96)/(exp_ML_B64 - exp_ML_B96);
	  double bmu_Linf_tm = (amu_64_tm - amu_Linf_tm)/exp_ML_B64;
	  amu_tm_ph_red.distr.push_back( amu_Linf_tm+ bmu_Linf_tm*exp_ML_extr);
	  //tm SD
	  double amu_SD_96_tm= amu_SD_B96_tm_ph.distr[ijck];
	  double amu_SD_64_tm = amu_SD_tm_ph.distr[ijck];
	  double amu_SD_Linf_tm = (amu_SD_96_tm*exp_ML_B64 - amu_SD_64_tm*exp_ML_B96)/(exp_ML_B64 - exp_ML_B96);
	  double bmu_SD_Linf_tm = (amu_SD_64_tm - amu_SD_Linf_tm)/exp_ML_B64;
	  amu_SD_tm_ph_red.distr.push_back( amu_SD_Linf_tm+ bmu_SD_Linf_tm*exp_ML_extr);
	  
	  //OS W
	  double amu_96_OS= amu_B96_OS_ph.distr[ijck];
	  double amu_64_OS = amu_OS_ph.distr[ijck];
	  double amu_Linf_OS = (amu_96_OS*exp_ML_B64 - amu_64_OS*exp_ML_B96)/(exp_ML_B64 - exp_ML_B96);
	  double bmu_Linf_OS = (amu_64_OS - amu_Linf_OS)/exp_ML_B64;
	  amu_OS_ph_red.distr.push_back( amu_Linf_OS+ bmu_Linf_OS*exp_ML_extr);

	  //OS SD
	  double amu_SD_96_OS= amu_SD_B96_OS_ph.distr[ijck];
	  double amu_SD_64_OS = amu_SD_OS_ph.distr[ijck];
	  double amu_SD_Linf_OS = (amu_SD_96_OS*exp_ML_B64 - amu_SD_64_OS*exp_ML_B96)/(exp_ML_B64 - exp_ML_B96);
	  double bmu_SD_Linf_OS = (amu_SD_64_OS - amu_SD_Linf_OS)/exp_ML_B64;
	  amu_SD_OS_ph_red.distr.push_back( amu_SD_Linf_OS+ bmu_SD_Linf_OS*exp_ML_extr);



	  //SD tmins
	  for(int tmin_id=0; tmin_id<(signed)tmins.size();tmin_id++) {
	    //tm
	    double amu_SD_96_tm_tmins= amu_SD_B96_tmins_tm_ph[tmin_id].distr[ijck];
	    double amu_SD_64_tm_tmins= agm2_light_SD_tmins_distr_list[tmin_id].distr_list[iens].distr[ijck];
	    double amu_SD_Linf_tm_tmins= (amu_SD_96_tm_tmins*exp_ML_B64 - amu_SD_64_tm_tmins*exp_ML_B96)/(exp_ML_B64 - exp_ML_B96);
	    double bmu_SD_Linf_tm_tmins = (amu_SD_64_tm_tmins - amu_SD_Linf_tm_tmins)/exp_ML_B64;
	    amu_SD_Bextr_tmins_tm_ph[tmin_id].distr.push_back( amu_SD_Linf_tm_tmins+ bmu_SD_Linf_tm_tmins*exp_ML_extr);

	    //OS
	    double amu_SD_96_OS_tmins= amu_SD_B96_tmins_OS_ph[tmin_id].distr[ijck];
	    double amu_SD_64_OS_tmins= agm2_light_SD_OS_tmins_distr_list[tmin_id].distr_list[iens].distr[ijck];
	    double amu_SD_Linf_OS_tmins = (amu_SD_96_OS_tmins*exp_ML_B64 - amu_SD_64_OS_tmins*exp_ML_B96)/(exp_ML_B64 - exp_ML_B96);
	    double bmu_SD_Linf_OS_tmins = (amu_SD_64_OS_tmins - amu_SD_Linf_OS_tmins)/exp_ML_B64;
	    amu_SD_Bextr_tmins_OS_ph[tmin_id].distr.push_back(  amu_SD_Linf_OS_tmins+ bmu_SD_Linf_OS_tmins*exp_ML_extr);
	    
	  }
	  

	  for(int id_eps=0;id_eps<eps_win_size;id_eps++) {

	    //double exp_2ML_B64= exp( -1.0*Mpi_fit.distr_list[iens].distr[ijck]*L_list[iens]);
	    //double exp_2ML_B96= exp(-1.0*Mpi_fit.distr_list[b96_id].distr[ijck]*L_list[b96_id]);
	    //double exp_2ML_extr= exp(-1.0*Mpi_fit.distr_list[iens].distr[ijck]*5.46*fm_to_inv_Gev/a_distr_list.distr_list[iens].distr[ijck]);

	    double exp_2ML_B64= exp( -1.0*Mpi_phys_fake_distr.distr[ijck]*L_list[iens]);
	    double exp_2ML_B96= exp(-1.0*Mpi_phys_fake_distr.distr[ijck]*L_list[b96_id]);
	    double exp_2ML_extr= exp(-1.0*Mpi_phys_fake_distr.distr[ijck]*5.46*fm_to_inv_Gev/a_distr_list.distr_list[iens].distr[ijck]);

	    //tm
	    double amu_96_tm_eps= amu_B96_eps_tm_ph[id_eps].distr[ijck];
	    double amu_64_tm_eps = amu_eps_tm_ph[id_eps].distr[ijck];
	    double amu_Linf_tm_eps = (amu_96_tm_eps*exp_2ML_B64 - amu_64_tm_eps*exp_2ML_B96)/(exp_2ML_B64 - exp_2ML_B96);
	    double bmu_Linf_tm_eps = (amu_64_tm_eps - amu_Linf_tm_eps)/exp_2ML_B64;
	    amu_eps_tm_ph_red[id_eps].distr.push_back( amu_Linf_tm_eps+ bmu_Linf_tm_eps*exp_2ML_extr);
	    //OS
	    double amu_96_OS_eps= amu_B96_eps_OS_ph[id_eps].distr[ijck];
	    double amu_64_OS_eps = amu_eps_OS_ph[id_eps].distr[ijck];
	    double amu_Linf_OS_eps = (amu_96_OS_eps*exp_2ML_B64 - amu_64_OS_eps*exp_2ML_B96)/(exp_2ML_B64 - exp_2ML_B96);
	    double bmu_Linf_OS_eps = (amu_64_OS_eps - amu_Linf_OS_eps)/exp_2ML_B64;
	    amu_eps_OS_ph_red[id_eps].distr.push_back( amu_Linf_OS_eps+ bmu_Linf_OS_eps*exp_2ML_extr);
	  }
	  
	}

	//push back tmins
	for(int tmin_id=0; tmin_id<(signed)tmins.size();tmin_id++) {
	  amu_SD_tmins_physical_point_tm_red[tmin_id].distr_list.push_back( amu_SD_Bextr_tmins_tm_ph[tmin_id]);
	  amu_SD_tmins_physical_point_OS_red[tmin_id].distr_list.push_back( amu_SD_Bextr_tmins_OS_ph[tmin_id]);
	}
	
      }
      else {
	amu_tm_ph_red=amu_tm_ph; amu_OS_ph_red=amu_OS_ph; amu_eps_tm_ph_red= amu_eps_tm_ph; amu_eps_OS_ph_red= amu_eps_OS_ph;
	amu_SD_tm_ph_red=amu_SD_tm_ph; amu_SD_OS_ph_red=amu_SD_OS_ph;
	for(int tmin_id=0; tmin_id<(signed)tmins.size();tmin_id++) {
	  amu_SD_tmins_physical_point_tm_red[tmin_id].distr_list.push_back( agm2_light_SD_tmins_distr_list[tmin_id].distr_list[iens]);
	  amu_SD_tmins_physical_point_OS_red[tmin_id].distr_list.push_back( agm2_light_SD_OS_tmins_distr_list[tmin_id].distr_list[iens]);
	}
      }


        
      //push_back
      //W
      amu_W_physical_point_tm_red.distr_list.push_back(amu_tm_ph_red);
      amu_W_physical_point_OS_red.distr_list.push_back(amu_OS_ph_red);
      //SD
      amu_SD_physical_point_tm_red.distr_list.push_back(amu_SD_tm_ph_red);
      amu_SD_physical_point_OS_red.distr_list.push_back(amu_SD_OS_ph_red);
      
      for(int id_eps=0;id_eps<eps_win_size;id_eps++) {
	amu_W_eps_physical_point_tm_red[id_eps].distr_list.push_back(amu_eps_tm_ph_red[id_eps]);
	amu_W_eps_physical_point_OS_red[id_eps].distr_list.push_back(amu_eps_OS_ph_red[id_eps]);
      }
      
   
      

      
      a_distr_list_fixed_L.distr_list.push_back( a_distr_list.distr_list[iens]);
      Mpi_fit_fixed_L.distr_list.push_back( Mpi_fit.distr_list[iens]);
      fp_fit_fixed_L.distr_list.push_back( fp_fit.distr_list[iens]);
      L_list_fixed_L.push_back( L_list[iens]);
      V_light_Tag_fixed_L.push_back( V_light_1.Tag[iens]);

    }


    
    

  }


  //print raw
  //W
  Print_To_File(V_light_1.Tag , {L_list, a_list, amu_W_physical_point_tm.ave(), amu_W_physical_point_tm.err(), corr_W_physical_point_tm.ave(), corr_W_physical_point_tm.err()}   , "../data/gm2/light/tm/windows_tm_phys_point", "", "#ENS L a val err");
  Print_To_File(V_light_1.Tag , {L_list, a_list, amu_W_physical_point_OS.ave(), amu_W_physical_point_OS.err(), corr_W_physical_point_OS.ave(), corr_W_physical_point_OS.err()}   , "../data/gm2/light/OS/windows_OS_phys_point", "", "#ENS L a val err");
  //SD
  Print_To_File(V_light_1.Tag , {L_list, a_list, amu_SD_physical_point_tm.ave(), amu_SD_physical_point_tm.err()}   , "../data/gm2/light/tm/SD_tm_phys_point", "", "#ENS L a val err");
  Print_To_File(V_light_1.Tag , {L_list, a_list, amu_SD_physical_point_OS.ave(), amu_SD_physical_point_OS.err()}   , "../data/gm2/light/OS/SD_OS_phys_point", "", "#ENS L a val err");
  for(int id_eps=0;id_eps<eps_win_size;id_eps++) {
    Print_To_File(V_light_1.Tag, {L_list, a_list, amu_W_eps_list_tm[id_eps].ave(), amu_W_eps_list_tm[id_eps].err()}, "../data/gm2/light/tm/windows_eps_"+to_string_with_precision(eps_win[id_eps]/fm_to_inv_Gev,2)+".list", " ", "#ENS L a val err");
    Print_To_File(V_light_1.Tag, {L_list, a_list, amu_W_eps_list_OS[id_eps].ave(), amu_W_eps_list_OS[id_eps].err()}, "../data/gm2/light/OS/windows_eps_"+to_string_with_precision(eps_win[id_eps]/fm_to_inv_Gev,2)+".list", "", "#ENS L a val err");
    Print_To_File(V_light_1.Tag , {L_list, a_list, amu_W_eps_physical_point_tm[id_eps].ave(), amu_W_eps_physical_point_tm[id_eps].err(), (amu_W_eps_physical_point_tm[id_eps]-amu_W_eps_physical_point_OS[id_eps]).ave(), (amu_W_eps_physical_point_tm[id_eps]-amu_W_eps_physical_point_OS[id_eps]).err()}   , "../data/gm2/light/tm/windows_tm_phys_point_eps_"+to_string_with_precision(eps_win[id_eps]/fm_to_inv_Gev,2), "", "#ENS L a val err Delta(tm-OS)");
    Print_To_File(V_light_1.Tag , {L_list, a_list, amu_W_eps_physical_point_OS[id_eps].ave(), amu_W_eps_physical_point_OS[id_eps].err()}   , "../data/gm2/light/OS/windows_OS_phys_point_eps_"+to_string_with_precision(eps_win[id_eps]/fm_to_inv_Gev,2), "", "#ENS L a val err");
  }
  //print at fixed L~5.46fm
  //W
  Print_To_File(V_light_Tag_fixed_L,{L_list_fixed_L,(a_distr_list_fixed_L/fm_to_inv_Gev).ave(),  amu_W_physical_point_tm_red.ave(), amu_W_physical_point_tm_red.err(), GS_artifacts_W_win_tm}, "../data/gm2/light/tm/windows_tm_fixed_L.dat", "", "#Ens L a  agm2 GS_art");
  Print_To_File(V_light_Tag_fixed_L,{L_list_fixed_L,(a_distr_list_fixed_L/fm_to_inv_Gev).ave(),  amu_W_physical_point_OS_red.ave(), amu_W_physical_point_OS_red.err(), GS_artifacts_W_win_OS}, "../data/gm2/light/OS/windows_OS_fixed_L.dat", "", "#Ens L a  agm2");
  //SD
  Print_To_File(V_light_Tag_fixed_L,{L_list_fixed_L,(a_distr_list_fixed_L/fm_to_inv_Gev).ave(),  amu_SD_physical_point_tm_red.ave(), amu_SD_physical_point_tm_red.err()}, "../data/gm2/light/tm/SD_tm_fixed_L.dat", "", "#Ens L a  agm2 GS_art");
  Print_To_File(V_light_Tag_fixed_L,{L_list_fixed_L,(a_distr_list_fixed_L/fm_to_inv_Gev).ave(),  amu_SD_physical_point_OS_red.ave(), amu_SD_physical_point_OS_red.err()}, "../data/gm2/light/OS/SD_OS_fixed_L.dat", "", "#Ens L a  agm2");
  
  for(int id_eps=0; id_eps<eps_win_size;id_eps++) {

    distr_t_list Diff= amu_W_eps_physical_point_tm_red[id_eps] - amu_W_eps_physical_point_OS_red[id_eps];
    distr_t_list Ratio_minus_1 = amu_W_eps_physical_point_tm_red[id_eps]/amu_W_eps_physical_point_OS_red[id_eps] -1;
    distr_t_list Diff_over_Ratio_minus_1 = Diff/Ratio_minus_1;
    Print_To_File(V_light_Tag_fixed_L,{L_list_fixed_L,(a_distr_list_fixed_L/fm_to_inv_Gev).ave(), amu_W_eps_physical_point_tm_red[id_eps].ave(), amu_W_eps_physical_point_tm_red[id_eps].err(), (amu_W_eps_physical_point_tm_red[id_eps]-amu_W_eps_physical_point_OS_red[id_eps]).ave(),  (amu_W_eps_physical_point_tm_red[id_eps]-amu_W_eps_physical_point_OS_red[id_eps]).err(), (amu_W_eps_physical_point_tm_red[id_eps]/amu_W_eps_physical_point_OS_red[id_eps]).ave(),  (amu_W_eps_physical_point_tm_red[id_eps]/amu_W_eps_physical_point_OS_red[id_eps]).err(), Diff_over_Ratio_minus_1.ave()     , Diff_over_Ratio_minus_1.err()    }, "../data/gm2/light/tm/windows_tm_fixed_L_eps_"+to_string_with_precision(eps_win[id_eps]/fm_to_inv_Gev,2), "", "#Ens L a  agm2 Delta(tm-OS)   R(tm/OS)     D/(R-1)");
    Print_To_File(V_light_Tag_fixed_L,{L_list_fixed_L,(a_distr_list_fixed_L/fm_to_inv_Gev).ave(),  amu_W_eps_physical_point_OS_red[id_eps].ave(), amu_W_eps_physical_point_OS_red[id_eps].err(),  (amu_W_eps_physical_point_OS_red[id_eps]/amu_W_eps_physical_point_tm_red[id_eps]).ave(),  (amu_W_eps_physical_point_OS_red[id_eps]/amu_W_eps_physical_point_tm_red[id_eps]).err() }, "../data/gm2/light/OS/windows_OS_fixed_L_eps_"+to_string_with_precision(eps_win[id_eps]/fm_to_inv_Gev,2), "", "#Ens L a  agm2  R(OS/tm)");
  
  }

  for(int tmin_id=0;tmin_id<(signed)tmins.size();tmin_id++) {

    Print_To_File(V_light_Tag_fixed_L,{L_list_fixed_L,(a_distr_list_fixed_L/fm_to_inv_Gev).ave(),  amu_SD_tmins_physical_point_tm_red[tmin_id].ave(), amu_SD_tmins_physical_point_tm_red[tmin_id].err()}, "../data/gm2/light/tm/SD_tmins_"+to_string_with_precision(tmins[tmin_id],5)+"_tm_fixed_L.dat", "", "#Ens L a  agm2");
    Print_To_File(V_light_Tag_fixed_L,{L_list_fixed_L,(a_distr_list_fixed_L/fm_to_inv_Gev).ave(),  amu_SD_tmins_physical_point_OS_red[tmin_id].ave(), amu_SD_tmins_physical_point_OS_red[tmin_id].err()}, "../data/gm2/light/OS/SD_tmins_"+to_string_with_precision(tmins[tmin_id],5)+"_OS_fixed_L.dat", "", "#Ens L a  agm2");


  }
  
  //###########################################################################################################################################



  


  exit(-1);

















  

 
  //print informations on the SD windows cut at tmins
  for(int tmins_id=0; tmins_id<(signed)tmins.size(); tmins_id++) {

    //tm
    Print_To_File(V_light_1.Tag, {L_list, a_list, ml_list, Mpi_fit.ave(), Mpi_fit.err(), Mpi_OS_fit.ave(), Mpi_OS_fit.err(),  fp_fit.ave(), fp_fit.err(), Zv_fit.ave(), Zv_fit.err(), Za_fit.ave(), Za_fit.err(), Zp_ov_Zs_fit.ave(), Zp_ov_Zs_fit.err(), agm2_light_SD_tmins_distr_list[tmins_id].ave(), agm2_light_SD_tmins_distr_list[tmins_id].err(), pert_result_SD_list}, "../data/gm2/light/tm/windows_tmin_"+to_string_with_precision(tmins[tmins_id],4)+".list", "", "#ENS L a ml  Mpi_tm  Mpi_OS fp  Zv   Za   Zp/Zs   SD   SD_pert");
    //OS
    Print_To_File(V_light_1.Tag, {L_list, a_list, ml_list, Mpi_fit.ave(), Mpi_fit.err(), Mpi_OS_fit.ave(), Mpi_OS_fit.err(),  fp_fit.ave(), fp_fit.err(), Zv_fit.ave(), Zv_fit.err(), Za_fit.ave(), Za_fit.err(), Zp_ov_Zs_fit.ave(), Zp_ov_Zs_fit.err(), agm2_light_SD_OS_tmins_distr_list[tmins_id].ave(), agm2_light_SD_OS_tmins_distr_list[tmins_id].err(), pert_result_SD_list}, "../data/gm2/light/OS/windows_tmin_"+to_string_with_precision(tmins[tmins_id],4)+".list", "", "#ENS L a ml  Mpi_tm  Mpi_OS fp  Zv   Za   Zp/Zs  SD   SD_pert");


  }

  
  //disco light
  if(Include_light_disco) {
    
    Print_To_File(disco_light_Tags,{L_list_disco, a_list_disco, ml_list_disco, Mpi_fit_disco.ave(), Mpi_fit_disco.err(), Mpi_OS_fit_disco.ave(), Mpi_OS_fit_disco.err(), fp_fit_disco.ave(), fp_fit_disco.err(), Zv_fit_disco.ave(), Zv_fit_disco.err(), agm2_disco_light_W.ave(), agm2_disco_light_W.err(), agm2_disco_light_W_ELM.ave(), agm2_disco_light_W_ELM.err(), agm2_disco_light_SD.ave(), agm2_disco_light_SD.err(), agm2_disco_light_SD_ELM.ave(), agm2_disco_light_SD_ELM.err()} , "../data/gm2/light/disco/windows.list", "", "#ENS L a ml Mpi_tm Mpi_OS fp Zv W W(ELM) SD SD(ELM)");
    
    //improved
    Print_To_File(disco_impr_light_Tags,{a_list_disco_impr, agm2_disco_impr_light_W.ave(), agm2_disco_impr_light_W.err(), agm2_disco_impr_light_W_ELM.ave(), agm2_disco_impr_light_W_ELM.err(), agm2_disco_impr_light_SD.ave(), agm2_disco_impr_light_SD.err(), agm2_disco_impr_light_SD_ELM.ave(), agm2_disco_impr_light_SD_ELM.err()} , "../data/gm2/light/disco/impr_windows.list", "", "#ENS a  W W(ELM) SD SD(ELM)");

    //improved D
    Print_To_File(disco_impr_lightD_Tags,{a_list_disco_impr_D, agm2_disco_impr_lightD_W.ave(), agm2_disco_impr_lightD_W.err(), agm2_disco_impr_lightD_W_ELM.ave(), agm2_disco_impr_lightD_W_ELM.err(), agm2_disco_impr_lightD_SD.ave(), agm2_disco_impr_lightD_SD.err(), agm2_disco_impr_lightD_SD_ELM.ave(), agm2_disco_impr_lightD_SD_ELM.err()} , "../data/gm2/light/disco/impr_D_windows.list", "", "#ENS a  W W(ELM) SD SD(ELM)");
    
    //improved DD
    Print_To_File(disco_impr_lightDD_Tags,{a_list_disco_impr_DD, agm2_disco_impr_lightDD_W.ave(), agm2_disco_impr_lightDD_W.err(), agm2_disco_impr_lightDD_W_ELM.ave(), agm2_disco_impr_lightDD_W_ELM.err(), agm2_disco_impr_lightDD_SD.ave(), agm2_disco_impr_lightDD_SD.err(), agm2_disco_impr_lightDD_SD_ELM.ave(), agm2_disco_impr_lightDD_SD_ELM.err()} , "../data/gm2/light/disco/impr_DD_windows.list", "", "#ENS a  W W(ELM) SD SD(ELM)");
  }

  
  //strange
  //tm
  //L
  Print_To_File(V_strange_1_L.Tag, {L_strange_list, a_strange_list, ml_strange_list, Mpi_fit_strange.ave(), Mpi_fit_strange.err(),  MV_fit_strange_L.ave(), MV_fit_strange_L.err(),  agm2_strange_W_L.ave(), agm2_strange_W_L.err(), agm2_strange_W_ELM_L.ave(), agm2_strange_W_ELM_L.err(), agm2_strange_SD_L.ave(), agm2_strange_SD_L.err(), agm2_strange_SD_ELM_L.ave(), agm2_strange_SD_ELM_L.err()}, "../data/gm2/strange/tm_"+Extrapolation_strange_mode+"/windows_L.list", "", "#ENS L a ml Mpi_tm  MV W   W(ELM)   SD    SD(ELM)");
  //M
  Print_To_File(V_strange_1_M.Tag, {L_strange_list, a_strange_list, ml_strange_list, Mpi_fit_strange.ave(), Mpi_fit_strange.err(),  MV_fit_strange_M.ave(), MV_fit_strange_M.err(),  agm2_strange_W_M.ave(), agm2_strange_W_M.err(), agm2_strange_W_ELM_M.ave(), agm2_strange_W_ELM_M.err(), agm2_strange_SD_M.ave(), agm2_strange_SD_M.err(), agm2_strange_SD_ELM_M.ave(), agm2_strange_SD_ELM_M.err()}, "../data/gm2/strange/tm_"+Extrapolation_strange_mode+"/windows_M.list", "", "#ENS L a ml Mpi_tm  MV W   W(ELM)   SD    SD(ELM)");
  //Extr
  Print_To_File(V_strange_1_L.Tag, {L_strange_list, a_strange_list, ml_strange_list, Mpi_fit_strange.ave(), Mpi_fit_strange.err(),   agm2_strange_W_Extr.ave(), agm2_strange_W_Extr.err(), agm2_strange_W_ELM_Extr.ave(), agm2_strange_W_ELM_Extr.err(), agm2_strange_SD_Extr.ave(), agm2_strange_SD_Extr.err(), agm2_strange_SD_ELM_Extr.ave(), agm2_strange_SD_ELM_Extr.err()}, "../data/gm2/strange/tm_"+Extrapolation_strange_mode+"/windows_Extr.list", "", "#ENS L a ml Mpi_tm   W   W(ELM)   SD    SD(ELM)");

  
  //OS
  //L
  Print_To_File(V_strange_OS_1_L.Tag, {L_strange_list, a_strange_list, ml_strange_list, Mpi_fit_strange.ave(), Mpi_fit_strange.err(),  MV_fit_strange_OS_L.ave(), MV_fit_strange_OS_L.err(), agm2_strange_W_OS_L.ave(), agm2_strange_W_OS_L.err(), agm2_strange_W_ELM_OS_L.ave(), agm2_strange_W_ELM_OS_L.err(), agm2_strange_SD_OS_L.ave(), agm2_strange_SD_OS_L.err(), agm2_strange_SD_ELM_OS_L.ave(), agm2_strange_SD_ELM_OS_L.err()}, "../data/gm2/strange/OS_"+Extrapolation_strange_mode+"/windows_L.list", "", "#ENS L a ml Mpi_tm  MV W   W(ELM)   SD    SD(ELM)");
  //M
  Print_To_File(V_strange_OS_1_M.Tag, {L_strange_list, a_strange_list, ml_strange_list, Mpi_fit_strange.ave(), Mpi_fit_strange.err(), MV_fit_strange_OS_M.ave(), MV_fit_strange_OS_M.err(), agm2_strange_W_OS_M.ave(), agm2_strange_W_OS_M.err(), agm2_strange_W_ELM_OS_M.ave(), agm2_strange_W_ELM_OS_M.err(), agm2_strange_SD_OS_M.ave(), agm2_strange_SD_OS_M.err(), agm2_strange_SD_ELM_OS_M.ave(), agm2_strange_SD_ELM_OS_M.err()}, "../data/gm2/strange/OS_"+Extrapolation_strange_mode+"/windows_M.list", "", "#ENS L a ml Mpi_tm  MV W   W(ELM)   SD    SD(ELM)");
  //Extr
  Print_To_File(V_strange_OS_1_L.Tag, {L_strange_list, a_strange_list, ml_strange_list, Mpi_fit_strange.ave(), Mpi_fit_strange.err(),  agm2_strange_W_OS_Extr.ave(), agm2_strange_W_OS_Extr.err(), agm2_strange_W_ELM_OS_Extr.ave(), agm2_strange_W_ELM_OS_Extr.err(), agm2_strange_SD_OS_Extr.ave(), agm2_strange_SD_OS_Extr.err(), agm2_strange_SD_ELM_OS_Extr.ave(), agm2_strange_SD_ELM_OS_Extr.err()}, "../data/gm2/strange/OS_"+Extrapolation_strange_mode+"/windows_Extr.list", "", "#ENS L a ml Mpi_tm   W   W(ELM)   SD    SD(ELM)");

  //disco strange
  if(Include_strange_disco) {
    
    Print_To_File(disco_strange_Tags,{agm2_disco_strange_W.ave(), agm2_disco_strange_W.err(), agm2_disco_strange_SD.ave(), agm2_disco_strange_SD.err()} , "../data/gm2/strange/disco/windows.list", "", "#ENS W SD");
    
    //improved
    Print_To_File(disco_impr_strange_Tags,{agm2_disco_impr_strange_W.ave(), agm2_disco_impr_strange_W.err(), agm2_disco_impr_strange_SD.ave(), agm2_disco_impr_strange_SD.err()} , "../data/gm2/strange/disco/impr_windows.list", "", "#ENS W SD");
  }



  
  //charm  
  //tm
  //L
  Print_To_File(V_charm_1_L.Tag, {L_charm_list, a_charm_list, ml_charm_list, Mpi_fit_charm.ave(), Mpi_fit_charm.err(), Mpi_OS_fit_charm.ave(), Mpi_OS_fit_charm.err(),  MV_fit_charm_L.ave(), MV_fit_charm_L.err(),  agm2_charm_W_L.ave(), agm2_charm_W_L.err(), agm2_charm_W_ELM_L.ave(), agm2_charm_W_ELM_L.err(), agm2_charm_SD_L.ave(), agm2_charm_SD_L.err(), agm2_charm_SD_ELM_L.ave(), agm2_charm_SD_ELM_L.err()}, "../data/gm2/charm/tm_"+Extrapolation_charm_mode+"/windows_L.list", "", "#ENS L a ml Mpi_tm Mpi_OS  MV W   W(ELM)   SD    SD(ELM)");
  //M
  Print_To_File(V_charm_1_M.Tag, {L_charm_list, a_charm_list, ml_charm_list, Mpi_fit_charm.ave(), Mpi_fit_charm.err(), Mpi_OS_fit_charm.ave(), Mpi_OS_fit_charm.err(),  MV_fit_charm_M.ave(), MV_fit_charm_M.err(),  agm2_charm_W_M.ave(), agm2_charm_W_M.err(), agm2_charm_W_ELM_M.ave(), agm2_charm_W_ELM_M.err(), agm2_charm_SD_M.ave(), agm2_charm_SD_M.err(), agm2_charm_SD_ELM_M.ave(), agm2_charm_SD_ELM_M.err()}, "../data/gm2/charm/tm_"+Extrapolation_charm_mode+"/windows_M.list", "", "#ENS L a ml Mpi_tm Mpi_OS  MV W   W(ELM)   SD    SD(ELM)");
  //H
  Print_To_File(V_charm_1_H.Tag, {L_charm_list, a_charm_list, ml_charm_list, Mpi_fit_charm.ave(), Mpi_fit_charm.err(), Mpi_OS_fit_charm.ave(), Mpi_OS_fit_charm.err(),  MV_fit_charm_H.ave(), MV_fit_charm_H.err(),  agm2_charm_W_H.ave(), agm2_charm_W_H.err(), agm2_charm_W_ELM_H.ave(), agm2_charm_W_ELM_H.err(), agm2_charm_SD_H.ave(), agm2_charm_SD_H.err(), agm2_charm_SD_ELM_H.ave(), agm2_charm_SD_ELM_H.err()}, "../data/gm2/charm/tm_"+Extrapolation_charm_mode+"/windows_H.list", "", "#ENS L a ml Mpi_tm Mpi_OS MV W   W(ELM)   SD    SD(ELM)");
  //Extr
  Print_To_File(V_charm_1_L.Tag, {L_charm_list, a_charm_list, ml_charm_list, Mpi_fit_charm.ave(), Mpi_fit_charm.err(), Mpi_OS_fit_charm.ave(), Mpi_OS_fit_charm.err(),  agm2_charm_W_Extr.ave(), agm2_charm_W_Extr.err(), agm2_charm_W_ELM_Extr.ave(), agm2_charm_W_ELM_Extr.err(), agm2_charm_SD_Extr.ave(), agm2_charm_SD_Extr.err(), agm2_charm_SD_ELM_Extr.ave(), agm2_charm_SD_ELM_Extr.err()}, "../data/gm2/charm/tm_"+Extrapolation_charm_mode+"/windows_Extr.list", "", "#ENS L a ml Mpi_tm Mpi_OS  W   W(ELM)   SD    SD(ELM)");
  
  //OS
  //L
  Print_To_File(V_charm_OS_1_L.Tag, {L_charm_list, a_charm_list, ml_charm_list, Mpi_fit_charm.ave(), Mpi_fit_charm.err(), Mpi_OS_fit_charm.ave(), Mpi_OS_fit_charm.err(),  MV_fit_charm_OS_L.ave(), MV_fit_charm_OS_L.err(),  agm2_charm_W_OS_L.ave(), agm2_charm_W_OS_L.err(), agm2_charm_W_ELM_OS_L.ave(), agm2_charm_W_ELM_OS_L.err(), agm2_charm_SD_OS_L.ave(), agm2_charm_SD_OS_L.err(), agm2_charm_SD_ELM_OS_L.ave(), agm2_charm_SD_ELM_OS_L.err()}, "../data/gm2/charm/OS_"+Extrapolation_charm_mode+"/windows_L.list", "", "#ENS L a ml Mpi_tm Mpi_OS   MV   W   W(ELM)   SD    SD(ELM)");
  //M
  Print_To_File(V_charm_OS_1_M.Tag, {L_charm_list, a_charm_list, ml_charm_list, Mpi_fit_charm.ave(), Mpi_fit_charm.err(), Mpi_OS_fit_charm.ave(), Mpi_OS_fit_charm.err(), MV_fit_charm_OS_M.ave(), MV_fit_charm_OS_M.err(),  agm2_charm_W_OS_M.ave(), agm2_charm_W_OS_M.err(), agm2_charm_W_ELM_OS_M.ave(), agm2_charm_W_ELM_OS_M.err(), agm2_charm_SD_OS_M.ave(), agm2_charm_SD_OS_M.err(), agm2_charm_SD_ELM_OS_M.ave(), agm2_charm_SD_ELM_OS_M.err()}, "../data/gm2/charm/OS_"+Extrapolation_charm_mode+"/windows_M.list", "", "#ENS L a ml Mpi_tm Mpi_OS  MV   W   W(ELM)   SD    SD(ELM)");
  //H
  Print_To_File(V_charm_OS_1_H.Tag, {L_charm_list, a_charm_list, ml_charm_list, Mpi_fit_charm.ave(), Mpi_fit_charm.err(), Mpi_OS_fit_charm.ave(), Mpi_OS_fit_charm.err(),  MV_fit_charm_OS_H.ave(), MV_fit_charm_OS_H.err(),  agm2_charm_W_OS_H.ave(), agm2_charm_W_OS_H.err(), agm2_charm_W_ELM_OS_H.ave(), agm2_charm_W_ELM_OS_H.err(), agm2_charm_SD_OS_H.ave(), agm2_charm_SD_OS_H.err(), agm2_charm_SD_ELM_OS_H.ave(), agm2_charm_SD_ELM_OS_H.err()}, "../data/gm2/charm/OS_"+Extrapolation_charm_mode+"/windows_H.list", "", "#ENS L a ml Mpi_tm Mpi_OS  MV   W   W(ELM)   SD    SD(ELM)");
  //Extr
  Print_To_File(V_charm_OS_1_L.Tag, {L_charm_list, a_charm_list, ml_charm_list, Mpi_fit_charm.ave(), Mpi_fit_charm.err(), Mpi_OS_fit_charm.ave(), Mpi_OS_fit_charm.err(),  agm2_charm_W_OS_Extr.ave(), agm2_charm_W_OS_Extr.err(), agm2_charm_W_ELM_OS_Extr.ave(), agm2_charm_W_ELM_OS_Extr.err(), agm2_charm_SD_OS_Extr.ave(), agm2_charm_SD_OS_Extr.err(), agm2_charm_SD_ELM_OS_Extr.ave(), agm2_charm_SD_ELM_OS_Extr.err()}, "../data/gm2/charm/OS_"+Extrapolation_charm_mode+"/windows_Extr.list", "", "#ENS L a ml Mpi_tm Mpi_OS  W   W(ELM)   SD    SD(ELM)");

  //disco charm
  if(Include_charm_disco) {
    Print_To_File(disco_charm_Tags,{agm2_disco_charm_W.ave(), agm2_disco_charm_W.err(), agm2_disco_charm_SD.ave(), agm2_disco_charm_SD.err()} , "../data/gm2/charm/disco/windows.list", "", "#ENS W SD");
    //improved
    Print_To_File(disco_impr_charm_Tags,{agm2_disco_impr_charm_W.ave(), agm2_disco_impr_charm_W.err(), agm2_disco_impr_charm_SD.ave(), agm2_disco_impr_charm_SD.err()} , "../data/gm2/charm/disco/impr_windows.list", "", "#ENS W SD");
  }



  //off-diagonal disconnected
  if(Include_off_diagonal_disco) {
    
    //improved light-strange
    Print_To_File(disco_impr_light_strange_Tags,{agm2_disco_impr_light_strange_W.ave(), agm2_disco_impr_light_strange_W.err(), agm2_disco_impr_light_strange_SD.ave(), agm2_disco_impr_light_strange_SD.err()} , "../data/gm2/light_strange/disco/impr_windows.list", "", "#ENS W SD");
    //improved lightD-strange
    Print_To_File(disco_impr_lightD_strange_Tags,{agm2_disco_impr_lightD_strange_W.ave(), agm2_disco_impr_lightD_strange_W.err(), agm2_disco_impr_lightD_strange_SD.ave(), agm2_disco_impr_lightD_strange_SD.err()} , "../data/gm2/light_strange/disco/impr_D_windows.list", "", "#ENS W SD");

    //improved light-charm
    Print_To_File(disco_impr_light_charm_Tags,{agm2_disco_impr_light_charm_W.ave(), agm2_disco_impr_light_charm_W.err(), agm2_disco_impr_light_charm_SD.ave(), agm2_disco_impr_light_charm_SD.err()} , "../data/gm2/light_charm/disco/impr_windows.list", "", "#ENS W SD");
    //improved lightD-charm
    Print_To_File(disco_impr_lightD_charm_Tags,{agm2_disco_impr_lightD_charm_W.ave(), agm2_disco_impr_lightD_charm_W.err(), agm2_disco_impr_lightD_charm_SD.ave(), agm2_disco_impr_lightD_charm_SD.err()} , "../data/gm2/light_charm/disco/impr_D_windows.list", "", "#ENS W SD");

    //improved strange-charm
    Print_To_File(disco_impr_strange_charm_Tags,{agm2_disco_impr_strange_charm_W.ave(), agm2_disco_impr_strange_charm_W.err(), agm2_disco_impr_strange_charm_SD.ave(), agm2_disco_impr_strange_charm_SD.err()} , "../data/gm2/strange_charm/disco/impr_windows.list", "", "#ENS W SD");
    
  }


  distr_t_list comb_disco_W, comb_disco_SD;
  vector<string> comb_Tags;
  distr_t_list a_lat;
  Vfloat vols_disco;
  if(Include_charm_disco && Include_light_disco && Include_strange_disco && Include_off_diagonal_disco) {

    int counter=0;
    //combine all disconnected when possible
  
    for(int i_ens=0; i_ens<(signed)disco_impr_lightD_strange_Tags.size(); i_ens++) {

      //get started
      distr_t accumulate_disco_W = agm2_disco_impr_lightD_strange_W.distr_list[i_ens];
      distr_t accumulate_disco_SD = agm2_disco_impr_lightD_strange_SD.distr_list[i_ens];

      bool Add_ensemble=true;

      //get lattice spacing
      distr_t lat;
      if( disco_impr_lightD_strange_Tags[i_ens].substr(1,1) == "A") lat =a_A;
      else if( disco_impr_lightD_strange_Tags[i_ens].substr(1,1) == "B") lat =a_B;
      else if( disco_impr_lightD_strange_Tags[i_ens].substr(1,1) == "C") lat =a_C;
      else if( disco_impr_lightD_strange_Tags[i_ens].substr(1,1) == "D") lat =a_D;
      else crash("Ensemble: "+disco_impr_lightD_strange_Tags[i_ens]+" not recognised");

    
      //find lightD_charm
      int id_lightD_charm;
      bool Find_lightD_charm=false;
      for(int j_ens=0; j_ens < (signed)disco_impr_lightD_charm_Tags.size(); j_ens++) {
	if( (disco_impr_lightD_charm_Tags[j_ens] == disco_impr_lightD_strange_Tags[i_ens]) ) { Find_lightD_charm=true; id_lightD_charm= j_ens; break;}}
      if(!Find_lightD_charm) Add_ensemble=false;
      else { //add lightD_charm
	accumulate_disco_W = accumulate_disco_W + agm2_disco_impr_lightD_charm_W.distr_list[id_lightD_charm];
	accumulate_disco_SD = accumulate_disco_SD + agm2_disco_impr_lightD_charm_SD.distr_list[id_lightD_charm];
      }
    
      //find strange_charm
      int id_strange_charm;
      bool Find_strange_charm=false;
      for(int j_ens=0; j_ens < (signed)disco_impr_strange_charm_Tags.size(); j_ens++) {
	if( (disco_impr_strange_charm_Tags[j_ens] == disco_impr_lightD_strange_Tags[i_ens] ) ) { Find_strange_charm=true; id_strange_charm= j_ens; break;}}
      if(!Find_strange_charm) Add_ensemble=false;
      else { //add strange_charm
	accumulate_disco_W = accumulate_disco_W + agm2_disco_impr_strange_charm_W.distr_list[id_strange_charm];
	accumulate_disco_SD = accumulate_disco_SD + agm2_disco_impr_strange_charm_SD.distr_list[id_strange_charm];
      }    

      //find light
      int id_light;
      bool Find_light=false;
      for(int j_ens=0; j_ens < (signed)disco_impr_lightD.size; j_ens++) {
	if( disco_impr_lightD_Tags[j_ens] == disco_impr_lightD_strange_Tags[i_ens]) { Find_light=true; id_light= j_ens; break;}}
      if(!Find_light) Add_ensemble=false;
      else { //add light
	accumulate_disco_W = accumulate_disco_W + agm2_disco_impr_lightD_W.distr_list[id_light];
	accumulate_disco_SD = accumulate_disco_SD + agm2_disco_impr_lightD_SD.distr_list[id_light];
      }

      //find strange
      int id_strange;
      bool Find_strange=false;
      for(int j_ens=0; j_ens < (signed)disco_impr_strange.size; j_ens++) {
	if( disco_impr_strange_Tags[j_ens] == disco_impr_lightD_strange_Tags[i_ens]) { Find_strange=true; id_strange= j_ens; break;}}
      if(!Find_strange) Add_ensemble=false;
      else { //add strange
	accumulate_disco_W = accumulate_disco_W + agm2_disco_impr_strange_W.distr_list[id_strange];
	accumulate_disco_SD = accumulate_disco_SD + agm2_disco_impr_strange_SD.distr_list[id_strange];
      }

      //find charm
      int id_charm;
      bool Find_charm=false;
      for(int j_ens=0; j_ens < (signed)disco_impr_charm.size;j_ens++) {
	if( (disco_impr_charm_Tags[j_ens] == disco_impr_lightD_strange_Tags[i_ens])) {
	  Find_charm=true; id_charm= j_ens; break;}}
    
      if(!Find_charm) Add_ensemble=false;
      else { //add charm
	accumulate_disco_W = accumulate_disco_W + agm2_disco_impr_charm_W.distr_list[id_charm];
	accumulate_disco_SD = accumulate_disco_SD + agm2_disco_impr_charm_SD.distr_list[id_charm];
      }

      if(Add_ensemble) {
	comb_disco_W.distr_list.push_back( accumulate_disco_W);
	comb_disco_SD.distr_list.push_back( accumulate_disco_SD);
	comb_Tags.push_back( disco_impr_lightD_strange_Tags[i_ens]);
	a_lat.distr_list.push_back( lat/fm_to_inv_Gev);
	LatticeInfo LI;
	LI.LatInfo_new_ens(disco_impr_lightD_strange_Tags[i_ens]);
	vols_disco.push_back( LI.L);
	counter++;
	cout<<"Adding disco comb, Ensemble: "<<disco_impr_lightD_strange_Tags[i_ens]<<endl;
      }
    
    }

    //Print the result
    Print_To_File(comb_Tags, { a_lat.ave(), comb_disco_W.ave(), comb_disco_W.err(), comb_disco_SD.ave(), comb_disco_SD.err()} , "../data/gm2/total_disco/disco/windows.list", "", "#Ens  a   W    SD");

  }
  
  
  

  if(!Skip_total_light_calc && Reco_agm2_total) {
 
    //print fitted pars for all ensembles
    //tm
    Print_To_File(V_light_1.Tag, {L_list, a_list, ml_list, Mpi_fit.ave(), Mpi_fit.err(), Mpi_OS_fit.ave(), Mpi_OS_fit.err(),  fp_fit.ave(), fp_fit.err(), Zv_fit.ave(), Zv_fit.err(), Za_fit.ave(), Za_fit.err(), Zp_ov_Zs_fit.ave(), Zp_ov_Zs_fit.err(), Edual.ave(), Edual.err(), Rdual.ave(), Rdual.err(), Mrho.ave(), Mrho.err(), (Mrho*Mpi_fit/a_distr_list).ave(), (Mrho*Mpi_fit/a_distr_list).err(), gpi.ave(), gpi.err(), Kappa.ave(), Kappa.err(), M_corr.ave(), M_corr.err()}, "../data/gm2/light/tm/fit_pars.list", "", "#ENS L a ml  Mpi_tm  Mpi_OS fp  Zv   Za   Zp/Zs  Edual Rdual Mrho g kappa Mpi_corr");
    //OS
    Print_To_File(V_light_1.Tag, {L_list, a_list, ml_list, Mpi_fit.ave(), Mpi_fit.err(), Mpi_OS_fit.ave(), Mpi_OS_fit.err(),  fp_fit.ave(), fp_fit.err(), Zv_fit.ave(), Zv_fit.err(), Za_fit.ave(), Za_fit.err(), Zp_ov_Zs_fit.ave(), Zp_ov_Zs_fit.err(), Edual_OS.ave(), Edual_OS.err(), Rdual_OS.ave(), Rdual_OS.err(), Mrho_OS.ave(), Mrho_OS.err(), (Mrho_OS*Mpi_fit/a_distr_list).ave(), (Mrho_OS*Mpi_fit/a_distr_list).err(), gpi_OS.ave(), gpi_OS.err(), Kappa_OS.ave(), Kappa_OS.err(), M_corr_OS.ave(), M_corr_OS.err()}, "../data/gm2/light/OS/fit_pars.list", "", "#ENS L a ml  Mpi_tm  Mpi_OS fp  Zv   Za   Zp/Zs  Edual Rdual Mrho g kappa Mpi_corr");

  
  
    //tm and OS pion mass, decay constant, RCs and agm2
    //tm
    Print_To_File(V_light_1.Tag, {L_list, a_list, ml_list, Mpi_fit.ave(), Mpi_fit.err(), Mpi_OS_fit.ave(), Mpi_OS_fit.err(),  fp_fit.ave(), fp_fit.err(), Zv_fit.ave(), Zv_fit.err(), Za_fit.ave(), Za_fit.err(), Zp_ov_Zs_fit.ave(), Zp_ov_Zs_fit.err(), agm2_light_fit.ave(), agm2_light_fit.err(), agm2_light_2L_fit.ave(), agm2_light_2L_fit.err(), agm2_light_Lprime_fit.ave(), agm2_light_Lprime_fit.err()   , agm2_light_infL_fit.ave(), agm2_light_infL_fit.err()}, "../data/gm2/light/tm/agm2_fit.list", "", "#ENS L a ml  Mpi_tm  Mpi_OS fp  Zv   Za   Zp/Zs    agm2(L)    agm2(1.5L) agm2(Mpi*L=4.2)   agm2(infL)");
    //OS
    Print_To_File(V_light_1.Tag, {L_list, a_list, ml_list, Mpi_fit.ave(), Mpi_fit.err(), Mpi_OS_fit.ave(), Mpi_OS_fit.err(),  fp_fit.ave(), fp_fit.err(), Zv_fit.ave(), Zv_fit.err(), Za_fit.ave(), Za_fit.err(), Zp_ov_Zs_fit.ave(), Zp_ov_Zs_fit.err(), agm2_light_fit_OS.ave(), agm2_light_fit_OS.err(), agm2_light_2L_fit_OS.ave(), agm2_light_2L_fit_OS.err(), agm2_light_Lprime_fit_OS.ave(), agm2_light_Lprime_fit_OS.err()   , agm2_light_infL_fit_OS.ave(), agm2_light_infL_fit_OS.err()}, "../data/gm2/light/OS/agm2_fit.list", "", "#ENS L a ml  Mpi_tm  Mpi_OS fp  Zv   Za   Zp/Zs    agm2(L)    agm2(1.5L) agm2(Mpi*L=4.2)   agm2(infL)");


  }




  //########################### PRINT TO FILE PI(Q^2) #######################################

  int Qs_size= Qs2.size();

  for(int q=0; q < Qs_size; q++) {

    //tm light
    Print_To_File(V_light_1.Tag,{L_list, a_list, PI_Q2_light_tm[q].ave(), PI_Q2_light_tm[q].err(), PI_Q2_light_tm_pert_sub[q].ave(), PI_Q2_light_tm_pert_sub[q].err()}, "../data/PI_Q2/light/tm/PI_Q2_"+to_string_with_precision(Qs2[q],5)+".t", "", "#Ens L a PI(Q^2)   PI(Q^2,sub)"  );
    //OS light
    Print_To_File(V_light_1.Tag,{L_list, a_list, PI_Q2_light_OS[q].ave(), PI_Q2_light_OS[q].err(), PI_Q2_light_OS_pert_sub[q].ave(), PI_Q2_light_OS_pert_sub[q].err()}, "../data/PI_Q2/light/OS/PI_Q2_"+to_string_with_precision(Qs2[q],5)+".t", "", "#Ens L a PI(Q^2)   PI(Q^2,sub)"  );

    //tm strange
    Print_To_File(V_strange_1_L.Tag,{L_strange_list, a_strange_list, PI_Q2_strange_tm[q].ave(), PI_Q2_strange_tm[q].err(), PI_Q2_strange_tm_pert_sub[q].ave(), PI_Q2_strange_tm_pert_sub[q].err()}, "../data/PI_Q2/strange/tm_"+Extrapolation_strange_mode+"/PI_Q2_"+to_string_with_precision(Qs2[q],5)+".t", "", "#Ens L a PI(Q^2)   PI(Q^2,sub)"  );
    //OS strange
    Print_To_File(V_strange_1_L.Tag,{L_strange_list, a_strange_list, PI_Q2_strange_OS[q].ave(), PI_Q2_strange_OS[q].err(), PI_Q2_strange_OS_pert_sub[q].ave(), PI_Q2_strange_OS_pert_sub[q].err()}, "../data/PI_Q2/strange/OS_"+Extrapolation_strange_mode+"/PI_Q2_"+to_string_with_precision(Qs2[q],5)+".t", "", "#Ens L a PI(Q^2)   PI(Q^2,sub)"  );

    //tm charm
    Print_To_File(V_charm_1_L.Tag,{L_charm_list, a_charm_list, PI_Q2_charm_tm[q].ave(), PI_Q2_charm_tm[q].err(), PI_Q2_charm_tm_pert_sub[q].ave(), PI_Q2_charm_tm_pert_sub[q].err()}, "../data/PI_Q2/charm/tm_"+Extrapolation_charm_mode+"/PI_Q2_"+to_string_with_precision(Qs2[q],5)+".t", "", "#Ens L a PI(Q^2)   PI(Q^2,sub)"  );
    //OS charm
    Print_To_File(V_charm_1_L.Tag,{L_charm_list, a_charm_list, PI_Q2_charm_OS[q].ave(), PI_Q2_charm_OS[q].err(), PI_Q2_charm_OS_pert_sub[q].ave(), PI_Q2_charm_OS_pert_sub[q].err()}, "../data/PI_Q2/charm/OS_"+Extrapolation_charm_mode+"/PI_Q2_"+to_string_with_precision(Qs2[q],5)+".t", "", "#Ens L a PI(Q^2)   PI(Q^2,sub)"  );



  }


  //#########################################################################################





  
  //PRINT DISCONNECTED BOUNDING FOR PI(Q^2) analysis
  //###########################################################
  if(Include_light_disco && Include_charm_disco && Include_strange_disco && Include_off_diagonal_disco) {
  for(int iens_disc=0; iens_disc<(signed)CORR_DISCO_FOR_PI_Q2.size();iens_disc++) {
    distr_t_list DISC_PI_PER_ENS(UseJack);
    Vint Tdatas_opt_disc;
    Bounding_PI_q2_disco(DISC_PI_PER_ENS, CORR_DISCO_FOR_PI_Q2[iens_disc] +0.1*CORR_LIGHT_FOR_PI_Q2[iens_disc]+CORR_STRANGE_FOR_PI_Q2[iens_disc]+CORR_CHARM_FOR_PI_Q2[iens_disc], 0.1*LIGHT_PI_Q2_FOR_DISCO[iens_disc] + STRANGE_PI_Q2_FOR_DISCO[iens_disc] + CHARM_PI_Q2_FOR_DISCO[iens_disc] ,  a_disc_PI_Q2.distr_list[iens_disc], "../data/PI_Q2/disco/PI_Q2_Tdata_"+Ens_list_disc_PI_Q2[iens_disc], Tdatas_opt_disc, m_rho*a_disc_PI_Q2.distr_list[iens_disc]);
    Vfloat Tcut_f_disco;
    for(int q=0;q<Qs_size;q++) {
      if( Tdatas_opt_disc[q] > 0) Tcut_f_disco.push_back( 1.0);
    }
    Add_ens_val_PI_q2(PI_Q2_disco, DISC_PI_PER_ENS);

    //Print
    Print_To_File({}, {Qs2, DISC_PI_PER_ENS.ave(), DISC_PI_PER_ENS.err(), Tcut_f_disco } , "../data/PI_Q2/disco/PI_Q2_extr_"+Ens_list_disc_PI_Q2[iens_disc]+".t", "", "# Q2[GeV2]   PI(Q^2)   Tcut_f");
   
  }
  }
  //###########################################################


  //print Zv and Za (Hadronic and RI-MOM) from strange correlators
  for(int is=0;is<Nens_strange;is++) {
    cout<<"####################"<<endl;
    cout<<"Ens: "<<V_strange_1_L.Tag[is]<<endl;
    cout<<"Za hadr-RIMOM (s (l), s (h), Extr,  u) "<<Za_fit_strange.ave(is)<<" +- "<<Za_fit_strange.err(is)<<" , "<<Za_fit_strange_heavy.ave(is)<<" +- "<<Za_fit_strange_heavy.err(is)<<" , "<<Za_fit_strange_Extr.ave(is)<<" +- "<<Za_fit_strange_Extr.err(is)<<" , "<<Za_WI.ave(is)<<" +- "<<Za_WI.err(is)<<endl;
    cout<<"Zv hadr-RIMOM (s (l), s (h), Extr,  u) "<<Zv_fit_strange.ave(is)<<" +- "<<Zv_fit_strange.err(is)<<" , "<<Zv_fit_strange_heavy.ave(is)<<" +- "<<Zv_fit_strange_heavy.err(is)<<" , "<<Zv_fit_strange_Extr.ave(is)<<" +- "<<Zv_fit_strange_Extr.err(is)<<" , "<<Zv_WI.ave(is)<<" +- "<<Zv_WI.err(is)<<endl;
    cout<<"####################"<<endl;
  }
  //print Zv and Za (Hadronic) from charm correlators
  for(int is=0;is<Nens_charm;is++) {
    cout<<"####################"<<endl;
    cout<<"Ens: "<<V_charm_1_L.Tag[is]<<endl;
    cout<<"#######   ZA      ######### "<<endl;
    cout<<"Za hadr (c) (L) "<<Za_fit_charm_L.ave(is)<<" +- "<<Za_fit_charm_L.err(is)<<endl;
    cout<<"Za hadr (c) (M) "<<Za_fit_charm_M.ave(is)<<" +- "<<Za_fit_charm_M.err(is)<<endl;
    cout<<"Za hadr (c) (H) "<<Za_fit_charm_H.ave(is)<<" +- "<<Za_fit_charm_H.err(is)<<endl;
    cout<<"Za hadr (c) Extrapolated: "<<Za_fit_charm_Extr.ave(is)<<" +- "<<Za_fit_charm_Extr.err(is)<<endl;
    cout<<"#######   Zv      ######### "<<endl;
    cout<<"Zv hadr (c) (L) "<<Zv_fit_charm_L.ave(is)<<" +- "<<Zv_fit_charm_L.err(is)<<endl;
    cout<<"Zv hadr (c) (M) "<<Zv_fit_charm_M.ave(is)<<" +- "<<Zv_fit_charm_M.err(is)<<endl;
    cout<<"Zv hadr (c) (H) "<<Zv_fit_charm_H.ave(is)<<" +- "<<Zv_fit_charm_H.err(is)<<endl;
    cout<<"Zv hadr (c) Extrapolated: "<<Zv_fit_charm_Extr.ave(is)<<" +- "<<Zv_fit_charm_Extr.err(is)<<endl;
    cout<<"####################"<<endl;
  }

  cout<<"#### Zv and Za (Hadronic) from light correlator (charm run)"<<endl;
  //print Zv and Za (Hadronic) from charm run (light)
  for(int is=0;is<Nens_charm;is++) {
    cout<<"####################"<<endl;
    cout<<"Ens: "<<V_charm_1_L.Tag[is]<<endl;
    cout<<"#######   ZA      ######### "<<endl;
    cout<<Za_fit_charm_light.ave(is)<<" +- "<<Za_fit_charm_light.err(is)<<endl;
    cout<<"#######   ZV      ######### "<<endl;
    cout<<Zv_fit_charm_light.ave(is)<<" +- "<<Zv_fit_charm_light.err(is)<<endl;
  }


  //extrapolate Silvano's Za and Zv light and strange on the A ensembles to the physical pion point
  vector<distr_t> ZAs_silva_light({Za_fit_silvano_light_A_ens.distr_list[0], Za_fit_silvano_light_A_ens.distr_list[1], Za_fit_silvano_light_A_ens.distr_list[2]});
  vector<distr_t> ZVs_silva_light({Zv_fit_silvano_light_A_ens.distr_list[0], Zv_fit_silvano_light_A_ens.distr_list[1], Zv_fit_silvano_light_A_ens.distr_list[2]});
  vector<distr_t> ZAs_silva_charm({Za_fit_silvano_charm_A_ens.distr_list[0], Za_fit_silvano_charm_A_ens.distr_list[1], Za_fit_silvano_charm_A_ens.distr_list[2]});
  vector<distr_t> ZVs_silva_charm({Zv_fit_silvano_charm_A_ens.distr_list[0], Zv_fit_silvano_charm_A_ens.distr_list[1], Zv_fit_silvano_charm_A_ens.distr_list[2]});
  vector<distr_t> Mpis_silva({ Mpi_fit_silvano_A_ens.distr_list[0], Mpi_fit_silvano_A_ens.distr_list[1], Mpi_fit_silvano_A_ens.distr_list[2]});
   
  vector<distr_t> Mpis_silva_strange, ZAs_silva_strange, ZVs_silva_strange, Zp_ov_Zss_silva_strange;
  if(A_ens_strange_silvano_tags.size() != 0) {
    ZAs_silva_strange = {Za_fit_silvano_strange_A_ens.distr_list[0], Za_fit_silvano_strange_A_ens.distr_list[1], Za_fit_silvano_strange_A_ens.distr_list[2]};
    ZVs_silva_strange = {Zv_fit_silvano_strange_A_ens.distr_list[0], Zv_fit_silvano_strange_A_ens.distr_list[1], Zv_fit_silvano_strange_A_ens.distr_list[2]};
    Zp_ov_Zss_silva_strange =  {Zp_ov_Zs_fit_silvano_strange_A_ens.distr_list[0], Zp_ov_Zs_fit_silvano_strange_A_ens.distr_list[1], Zp_ov_Zs_fit_silvano_strange_A_ens.distr_list[2]};
    for(int tag_i = 0; tag_i < (signed)A_ens_strange_silvano_tags.size(); tag_i++) {
      bool find_tag_A_ens_charm=false;
      int tag_id=0;
      for(int tag_j =0 ;tag_j< (signed)A_ens_charm_silvano_tags.size(); tag_j++) {
	if(A_ens_strange_silvano_tags[tag_i] == A_ens_charm_silvano_tags[tag_j]) { tag_id = tag_j; find_tag_A_ens_charm=true; break;}
      }
      if(!find_tag_A_ens_charm) crash("Cannot find Ensemble tag : "+A_ens_strange_silvano_tags[tag_i]+ " in A_ens_charm_silvano_tags");
      Mpis_silva_strange.push_back( Mpis_silva[tag_id]);
    }
  }
    
    
    
    

  distr_t Zv_extr_silvano_light = Obs_extrapolation_meson_mass( ZVs_silva_light, Mpis_silva, m_pi*m_pi, "../data/gm2/charm", "Zv_light_A_ens_extr", UseJack, "FIT" );
  distr_t Za_extr_silvano_light = Obs_extrapolation_meson_mass( ZAs_silva_light, Mpis_silva, m_pi*m_pi, "../data/gm2/charm", "Za_light_A_ens_extr", UseJack, "FIT" );
  distr_t Zv_extr_silvano_charm = Obs_extrapolation_meson_mass( ZVs_silva_charm, Mpis_silva, m_pi*m_pi, "../data/gm2/charm", "Zv_charm_A_ens_extr_"+Extrapolation_charm_mode, UseJack, "FIT" );
  distr_t Za_extr_silvano_charm = Obs_extrapolation_meson_mass( ZAs_silva_charm, Mpis_silva, m_pi*m_pi, "../data/gm2/charm", "Za_charm_A_ens_extr_"+Extrapolation_charm_mode, UseJack, "FIT" );
  distr_t Zv_extr_silvano_strange, Za_extr_silvano_strange, Zp_ov_Zs_extr_silvano_strange;
  if(A_ens_strange_silvano_tags.size() != 0) {
    Zv_extr_silvano_strange = Obs_extrapolation_meson_mass( ZVs_silva_strange, Mpis_silva_strange, m_pi*m_pi, "../data/gm2/strange", "Zv_strange_A_ens_extr_"+Extrapolation_strange_mode, UseJack, "FIT" );
    Za_extr_silvano_strange = Obs_extrapolation_meson_mass( ZAs_silva_strange, Mpis_silva_strange, m_pi*m_pi, "../data/gm2/strange", "Za_strange_A_ens_extr_"+Extrapolation_strange_mode, UseJack, "FIT" );
    Zp_ov_Zs_extr_silvano_strange = Obs_extrapolation_meson_mass( Zp_ov_Zss_silva_strange, Mpis_silva_strange, m_pi*m_pi, "../data/gm2/strange", "Zp_ov_Zs_strange_A_ens_extr_"+Extrapolation_strange_mode, UseJack, "FIT");
  }

  cout<<"Extrapolated Zv and Za light from Silvano's runs: "<<endl;
  cout<<"Za (A) : "<<Za_extr_silvano_light.ave()<<" +- "<<Za_extr_silvano_light.err()<<endl;
  cout<<"Zv (A) : "<<Zv_extr_silvano_light.ave()<<" +- "<<Zv_extr_silvano_light.err()<<endl;
  cout<<"#########################################"<<endl;


  if(A_ens_strange_silvano_tags.size() != 0 ) {
    cout<<"Extrapolated Zv, Za, Zp_ov_Zs strange from Silvano's runs: "<<endl;
    cout<<"Za (A) : "<<Za_extr_silvano_strange.ave()<<" +- "<<Za_extr_silvano_strange.err()<<endl;
    cout<<"Zv (A) : "<<Zv_extr_silvano_strange.ave()<<" +- "<<Zv_extr_silvano_strange.err()<<endl;
    cout<<"Zp_ov_Zs (A) : "<<Zp_ov_Zs_extr_silvano_strange.ave()<<" +- "<<Zp_ov_Zs_extr_silvano_strange.err()<<endl;
    cout<<"#########################################"<<endl;
  }

  cout<<"Extrapolated Zv and Za charm from Silvano's runs: "<<endl;
  cout<<"Za (A) : "<<Za_extr_silvano_charm.ave()<<" +- "<<Za_extr_silvano_charm.err()<<endl;
  cout<<"Zv (A) : "<<Zv_extr_silvano_charm.ave()<<" +- "<<Zv_extr_silvano_charm.err()<<endl;
  cout<<"#########################################"<<endl;



  // Extrapolate mc diff from J/psi and D meson on A,B,C ensembles
  distr_t mc_A_Jpsi, mc_B_Jpsi, mc_C_Jpsi;
  distr_t mc_A_D, mc_B_D, mc_C_D;

  for(int ijack=0; ijack<Njacks;ijack++) {
    mc_A_Jpsi.distr.push_back( 0.2914580943 + GM()*0.002032600526/sqrt(Njacks -1.0));
    mc_B_Jpsi.distr.push_back( 0.2467282673 + GM()*0.001476345086/sqrt(Njacks -1.0));
    mc_C_Jpsi.distr.push_back( 0.2075874971 + GM()*0.001211427752/sqrt(Njacks -1.0));

    mc_A_D.distr.push_back( 0.281766 + GM()*0.00378329/sqrt(Njacks -1.0));
    mc_B_D.distr.push_back( 0.241232 + GM()*0.00296929/sqrt(Njacks -1.0));
    mc_C_D.distr.push_back( 0.204887  + GM()*0.00239961/sqrt(Njacks -1.0));

  }

  vector<distr_t> mc_diff_list({ (mc_A_Jpsi-mc_A_D)/a_A, (mc_B_Jpsi-mc_B_D)/a_B, (mc_C_Jpsi-mc_C_D)/a_C });
  vector<distr_t> a2_ABC_list( { a_A*a_A, a_B*a_B, a_C*a_C});
  distr_t res_mc_diff = Obs_extrapolation_meson_mass( mc_diff_list, a2_ABC_list,  0.0, "../data/gm2/charm", "mc_diff_D_Jpsi_fit", UseJack, "FIT");


  //########## QUARK MASS ANALYSIS FOR STRANGE AND CHARM ########################


  //COMPUTE ZP
  distr_t ZP_A(UseJack), ZP_B(UseJack), ZP_C(UseJack), ZP_D(UseJack);
  distr_t Conv_to_3_GeV(UseJack), Conv_to_2_GeV(UseJack);

  for(int ijack=0; ijack<Njacks;ijack++) {
    ZP_A.distr.push_back( 0.569 + GM()*0.00583/sqrt(Njacks -1.0));
    ZP_B.distr.push_back( 0.574 + GM()*0.00640/sqrt(Njacks -1.0));
    ZP_C.distr.push_back( 0.584 + GM()*0.00583/sqrt(Njacks -1.0));
    ZP_D.distr.push_back( 0.584*1.0046948 + GM()*0.0052/sqrt(Njacks -1.0));
    Conv_to_3_GeV.distr.push_back( 0.92570 + GM()*0.00034/sqrt(Njacks -1.0));
    Conv_to_2_GeV.distr.push_back(  0.83416 + GM()*0.00086/sqrt(Njacks -1.0));

  }


  distr_t_list ams_phi_list_ren(UseJack), ams_etas_list_ren(UseJack);
  vector<distr_t> ms_phi_list, ms_etas_list;
  vector<distr_t> ms_diff_list;
  vector<distr_t> a2_ms_diff_list;
  vector<distr_t> ms_phi_etas_ratio_list;
  for(int iens_strange=0;iens_strange<Nens_strange;iens_strange++) {

    distr_t Z_M(UseJack);
    if(V_strange_1_L.Tag[iens_strange].substr(1,1) == "A") Z_M = 1.0/(ZP_A*Conv_to_2_GeV);
    else if(V_strange_1_L.Tag[iens_strange].substr(1,1) == "B") Z_M = 1.0/(ZP_B*Conv_to_2_GeV);
    else if(V_strange_1_L.Tag[iens_strange].substr(1,1) == "C") Z_M = 1.0/(ZP_C*Conv_to_2_GeV);
    else if(V_strange_1_L.Tag[iens_strange].substr(1,1) == "D") Z_M = 1.0/(ZP_D*Conv_to_2_GeV);
    else crash("In Determining Z_M (strange), Ensemble : "+V_strange_1_L.Tag[iens_strange]+" not found");

    ams_etas_list_ren.distr_list.push_back( Z_M*ms_extr_etas_list.distr_list[iens_strange]*a_distr_list_strange.distr_list[iens_strange]);
    ams_phi_list_ren.distr_list.push_back( Z_M*ms_extr_phi_list.distr_list[iens_strange]*a_distr_list_strange.distr_list[iens_strange]);
      
    if( V_strange_1_L.Tag[iens_strange].substr(1,1) != "A" && V_strange_1_L.Tag[iens_strange] != "cB211b.072.96" ) {

      
      
      ms_phi_list.push_back( Z_M*ms_extr_phi_list.distr_list[iens_strange]);
      ms_etas_list.push_back( Z_M*ms_extr_etas_list.distr_list[iens_strange]);
      ms_diff_list.push_back( Z_M*ms_extr_diff_list.distr_list[iens_strange]);
      ms_phi_etas_ratio_list.push_back( ms_extr_phi_list.distr_list[iens_strange]/ms_extr_etas_list.distr_list[iens_strange]);
      a2_ms_diff_list.push_back( a_distr_list_strange.distr_list[iens_strange]*a_distr_list_strange.distr_list[iens_strange]);
    }
  }
  distr_t res_ms_diff = Obs_extrapolation_meson_mass( ms_diff_list, a2_ms_diff_list,  0.0, "../data/gm2/strange", "ms_diff_etas_phi_fit", UseJack, "FIT");
  distr_t res_ms_phi = Obs_extrapolation_meson_mass( ms_phi_list, a2_ms_diff_list,  0.0, "../data/gm2/strange", "ms_cont_limit_phi", UseJack, "FIT");
  distr_t res_ms_etas = Obs_extrapolation_meson_mass( ms_etas_list, a2_ms_diff_list,  0.0, "../data/gm2/strange", "ms_cont_limit_etas", UseJack, "FIT");
  distr_t res_ms_phi_etas_ratio = Obs_extrapolation_meson_mass( ms_phi_etas_ratio_list, a2_ms_diff_list,  0.0, "../data/gm2/strange", "ms_phi_etas_ratio_cont_limit", UseJack, "FIT");

  Print_To_File({V_strange_1_L.Tag}, { ams_etas_list_ren.ave(), ams_etas_list_ren.err(), ams_phi_list_ren.ave(), ams_phi_list_ren.err()}, "../data/gm2/strange/ms_ren_2_GeV.list", "", "#Ens ms(etas)  ms(phi)");


  distr_t_list amc_Jpsi_list_ren(UseJack), amc_etac_list_ren(UseJack);
  vector<distr_t> mc_etac_list, mc_Jpsi_list;
  vector<distr_t> mc_Jpsi_etac_ratio_list;
  vector<distr_t> mc_etac_finest_list, mc_Jpsi_finest_list;
  vector<distr_t> mc_Jpsi_etac_ratio_finest_list;
  vector<distr_t> a2_mc_diff_list, a2_mc_diff_finest_list;
  for(int iens_charm=0;iens_charm<Nens_charm;iens_charm++) {

    distr_t Z_M(UseJack);
    if(V_charm_1_L.Tag[iens_charm].substr(1,1) == "A") Z_M = 1.0/(ZP_A*Conv_to_3_GeV);
    else if(V_charm_1_L.Tag[iens_charm].substr(1,1) == "B") Z_M = 1.0/(ZP_B*Conv_to_3_GeV);
    else if(V_charm_1_L.Tag[iens_charm].substr(1,1) == "C") Z_M = 1.0/(ZP_C*Conv_to_3_GeV);
    else if(V_charm_1_L.Tag[iens_charm].substr(1,1) == "D") Z_M = 1.0/(ZP_D*Conv_to_3_GeV);
    else crash("In Determining Z_M (charm), Ensemble : "+V_charm_1_L.Tag[iens_charm]+" not found");


    amc_etac_list_ren.distr_list.push_back( Z_M*mc_extr_etac_list.distr_list[iens_charm]*a_distr_list_charm.distr_list[iens_charm]);
    amc_Jpsi_list_ren.distr_list.push_back( Z_M*mc_extr_Jpsi_list.distr_list[iens_charm]*a_distr_list_charm.distr_list[iens_charm]);
    
    if( V_charm_1_L.Tag[iens_charm].substr(1,1) != "A" || V_charm_1_L.Tag[iens_charm] == "cA211ab.30.32") {


      
      
      mc_Jpsi_list.push_back( Z_M*mc_extr_Jpsi_list.distr_list[iens_charm]);
      mc_etac_list.push_back( Z_M*mc_extr_etac_list.distr_list[iens_charm]);
      mc_Jpsi_etac_ratio_list.push_back( mc_extr_Jpsi_list.distr_list[iens_charm]/mc_extr_etac_list.distr_list[iens_charm]);
      a2_mc_diff_list.push_back( a_distr_list_charm.distr_list[iens_charm]*a_distr_list_charm.distr_list[iens_charm]);
      if( V_charm_1_L.Tag[iens_charm].substr(1,1) != "A") {
	mc_Jpsi_finest_list.push_back( Z_M*mc_extr_Jpsi_list.distr_list[iens_charm]);
	mc_etac_finest_list.push_back( Z_M*mc_extr_etac_list.distr_list[iens_charm]);
	mc_Jpsi_etac_ratio_finest_list.push_back( mc_extr_Jpsi_list.distr_list[iens_charm]/mc_extr_etac_list.distr_list[iens_charm]);
	a2_mc_diff_finest_list.push_back( a_distr_list_charm.distr_list[iens_charm]*a_distr_list_charm.distr_list[iens_charm]);

      }
   
    }
    
  }
  distr_t res_mc_Jpsi = Obs_extrapolation_meson_mass( mc_Jpsi_list, a2_mc_diff_list,  0.0, "../data/gm2/charm", "mc_cont_limit_Jpsi", UseJack, "FIT");
  distr_t res_mc_etac = Obs_extrapolation_meson_mass( mc_etac_list, a2_mc_diff_list,  0.0, "../data/gm2/charm", "mc_cont_limit_etac", UseJack, "FIT");
  distr_t res_mc_Jpsi_etac_ratio = Obs_extrapolation_meson_mass( mc_Jpsi_etac_ratio_list, a2_mc_diff_list,  0.0, "../data/gm2/charm", "mc_Jpsi_etac_ratio_cont_limit", UseJack, "FIT");
  distr_t res_mc_Jpsi_finest = Obs_extrapolation_meson_mass( mc_Jpsi_finest_list, a2_mc_diff_finest_list,  0.0, "../data/gm2/charm", "mc_finest_cont_limit_Jpsi", UseJack, "FIT");
  distr_t res_mc_etac_finest = Obs_extrapolation_meson_mass( mc_etac_finest_list, a2_mc_diff_finest_list,  0.0, "../data/gm2/charm", "mc_finest_cont_limit_etac", UseJack, "FIT");
  distr_t res_mc_Jpsi_etac_ratio_finest = Obs_extrapolation_meson_mass( mc_Jpsi_etac_ratio_finest_list, a2_mc_diff_finest_list,  0.0, "../data/gm2/charm", "mc_finest_Jpsi_etac_ratio_cont_limit", UseJack, "FIT");


  Print_To_File({V_charm_1_L.Tag}, { amc_etac_list_ren.ave(), amc_etac_list_ren.err(), amc_Jpsi_list_ren.ave(), amc_Jpsi_list_ren.err()}, "../data/gm2/charm/mc_ren_3_GeV.list", "", "#Ens mc(etac)  mc(Jpsi)");



  //#############################################        CONTINUUM/THERMODYNAMIC/PHYSICAL-POINT EXTRAPOLATION       ###############################################
  vector<string> a2_list({"OS"});
  vector<string> FSEs_list({"off"});
  vector<string> a4_list({"off"});
  vector<string> mass_extr_list({"off"});
  vector<string> single_fit_list({"OS"});
  VPfloat n_m_pair_list({make_pair(0,0), make_pair(3,0), make_pair(0,3), make_pair(3,3), make_pair(1,0), make_pair(0,1), make_pair(1,1), make_pair(2,0), make_pair(0,2), make_pair(2,2)});
  bool allow_a4_and_log= false;
  bool allow_only_finest= false;
  int w0s_mult=2;
  distr_t return4_val;


  



  
  //light
  //#####################################################################################################################
  //#####################################################################################################################
  //#####################################################################################################################

  //epsilon-window computation

 
  for(int id_eps=0;id_eps<eps_win_size;id_eps++) {

     allow_a4_and_log=false;
     mass_extr_list={"off"};
     FSEs_list = {"off"};
     single_fit_list ={"off"};
     a4_list={"on"};
     a2_list={"on"};
     n_m_pair_list= {make_pair(0,0)};

    distr_t return_val;
    Perform_Akaike_fits(amu_W_eps_physical_point_tm_red[id_eps], amu_W_eps_physical_point_OS_red[id_eps], a_A, a_B, a_C, a_D, L_list_fixed_L, a_distr_list_fixed_L, Mpi_fit_fixed_L,fp_fit_fixed_L, V_light_Tag_fixed_L, UseJack, Njacks, Nboots, "W_win_fixed_L_eps_"+to_string_with_precision(eps_win[id_eps]/fm_to_inv_Gev,2), "light",a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list,allow_a4_and_log,allow_only_finest, 1, false, 0, 0, 200.0, LL, 0.0, return_val, return4_val);

    single_fit_list={"tm"};
    a4_list={"tm"};
    a2_list={"tm"};
    distr_t D1_diff, D4_diff;
    Perform_Akaike_fits(amu_W_eps_physical_point_tm_red[id_eps]-amu_W_eps_physical_point_OS_red[id_eps], amu_W_eps_physical_point_OS_red[id_eps], a_A, a_B, a_C, a_D, L_list_fixed_L, a_distr_list_fixed_L, Mpi_fit_fixed_L,fp_fit_fixed_L, V_light_Tag_fixed_L, UseJack, Njacks, Nboots, "diff_tm_OS_fixed_L_eps_"+to_string_with_precision(eps_win[id_eps]/fm_to_inv_Gev,2), "light",a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list,allow_a4_and_log,allow_only_finest, 1, false, 0, 0, 200.0, LL, 0.0, D1_diff, D4_diff);
    distr_t R1_ratio, R4_ratio;
    Perform_Akaike_fits(1.0e-10*amu_W_eps_physical_point_tm_red[id_eps]/amu_W_eps_physical_point_OS_red[id_eps], amu_W_eps_physical_point_OS_red[id_eps], a_A, a_B, a_C, a_D, L_list_fixed_L, a_distr_list_fixed_L, Mpi_fit_fixed_L,fp_fit_fixed_L, V_light_Tag_fixed_L, UseJack, Njacks, Nboots, "ratio_tm_OS_fixed_L_eps_"+to_string_with_precision(eps_win[id_eps]/fm_to_inv_Gev,2), "light",a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list,allow_a4_and_log,allow_only_finest, 1, false, 0, 0, 200.0, LL, 0.0, R1_ratio, R4_ratio);
    distr_t R2_ratio, R4_bis_ratio;
    Perform_Akaike_fits(1.0e-10*amu_W_eps_physical_point_OS_red[id_eps]/amu_W_eps_physical_point_tm_red[id_eps], amu_W_eps_physical_point_OS_red[id_eps], a_A, a_B, a_C, a_D, L_list_fixed_L, a_distr_list_fixed_L, Mpi_fit_fixed_L,fp_fit_fixed_L, V_light_Tag_fixed_L, UseJack, Njacks, Nboots, "ratio_tm_OS_rev_fixed_L_eps_"+to_string_with_precision(eps_win[id_eps]/fm_to_inv_Gev,2), "light",a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list,allow_a4_and_log,allow_only_finest, 1, false, 0, 0, 200.0, LL, 0.0, R2_ratio, R4_bis_ratio);


    auto f_func_tm= [&D1_diff, &D4_diff, &R1_ratio, &R4_ratio, &R2_ratio, &R4_bis_ratio](double a) -> distr_t {

		   distr_t R1= 1.0+ R1_ratio*a*a + R4_ratio*a*a*a*a;
		   distr_t R2= 1.0 +R2_ratio*a*a + R4_bis_ratio*a*a*a*a;
		   distr_t D = D1_diff*a*a + D4_diff*a*a*a*a;

		   return 0.5*D*( R1/(R1-1.0) -1.0/(R2-1.0));

		    };

    auto f_func_OS= [&D1_diff, &D4_diff, &R1_ratio, &R4_ratio, &R2_ratio, &R4_bis_ratio](double a) -> distr_t {

		   distr_t R1= 1.0+ R1_ratio*a*a + R4_ratio*a*a*a*a;
		   distr_t R2= 1.0 +R2_ratio*a*a + R4_bis_ratio*a*a*a*a;
		   distr_t D = D1_diff*a*a + D4_diff*a*a*a*a;

		   return 0.5*D*( 1.0/(R1-1.0) -R2/(R2-1.0));

		    };


    int alat_points=40;
    Vfloat apoints;
    distr_t_list sampling_tm(UseJack);
    distr_t_list sampling_OS(UseJack);
    for(int ilat=0;ilat<alat_points;ilat++) {
      double al = ilat*2.0*a_B.ave()/alat_points;
      apoints.push_back(al);
      sampling_tm.distr_list.push_back( f_func_tm(al));
      sampling_OS.distr_list.push_back( f_func_OS(al));
    }

    //print
    Print_To_File({}, {apoints,sampling_tm.ave(),sampling_tm.err(), sampling_OS.ave(), sampling_OS.err()}, "../data/gm2/light/windows_fit_func/ratio_tm_OS_fixed_L_eps_"+to_string_with_precision(eps_win[id_eps]/fm_to_inv_Gev,2)+"/windows_fitting_func.dat",   "", "#a  tm   OS");

      

    cout<<"EPSILON: "<<eps_win[id_eps]/fm_to_inv_Gev<<" d1/r1: "<< (D1_diff/R1_ratio).ave()<<" "<<(D1_diff/R1_ratio).err()<<" "<<(-1*D1_diff/R2_ratio).ave()<<" "<<(-1*D1_diff/R2_ratio).err()<<endl;

    single_fit_list={"off"};
    a4_list={"off"};
    a2_list={"on"};
    Perform_Akaike_fits(amu_W_eps_physical_point_tm_red[id_eps], amu_W_eps_physical_point_OS_red[id_eps], a_A, a_B, a_C, a_D, L_list_fixed_L, a_distr_list_fixed_L, Mpi_fit_fixed_L,fp_fit_fixed_L, V_light_Tag_fixed_L, UseJack, Njacks, Nboots, "W_win_fixed_L_a2_eps_"+to_string_with_precision(eps_win[id_eps]/fm_to_inv_Gev,2), "light",a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list,allow_a4_and_log,allow_only_finest, 1, false, 0, 0, 200.0, LL, 0.0, return_val, return4_val);


    single_fit_list={"tm", "OS", "off"};
    a2_list={"tm", "OS", "on"};
    a4_list={"off", "OS"};
    Perform_Akaike_fits(amu_W_eps_physical_point_tm_red[id_eps]+ ((id_eps==0)?1.0:0.0)*Get_id_jack_distr_list(amu_W_eps_physical_point_tm_red[id_eps].size(),Njacks)*GS_artifacts_W_win_tm, amu_W_eps_physical_point_OS_red[id_eps] + ((id_eps==0)?1.0:0.0)*Get_id_jack_distr_list(amu_W_eps_physical_point_OS_red[id_eps].size(),Njacks)*GS_artifacts_W_win_OS, a_A, a_B, a_C, a_D, L_list_fixed_L, a_distr_list_fixed_L, Mpi_fit_fixed_L,fp_fit_fixed_L, V_light_Tag_fixed_L, UseJack, Njacks, Nboots, "W_win_fixed_L_single_a2_eps_"+to_string_with_precision(eps_win[id_eps]/fm_to_inv_Gev,2), "light",a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list,allow_a4_and_log,allow_only_finest, 1, false, 0, 0, 200.0, LL, 0.0, return_val, return4_val);
    

    single_fit_list={"off"};
    a4_list={"off"};
    a2_list={"on"};

    distr_t_list amu_W_eps_tm_two_finest(UseJack);
    distr_t_list amu_W_eps_OS_two_finest(UseJack);
    Vfloat L_list_fixed_L_two_finest;
    distr_t_list a_distr_list_fixed_L_two_finest(UseJack);
    distr_t_list Mpi_fit_fixed_L_two_finest(UseJack);
    distr_t_list fp_fit_fixed_L_two_finest(UseJack);
    vector<string> V_light_Tag_fixed_L_two_finest;

    for(int iens_tf=0; iens_tf<V_light_Tag_fixed_L.size();iens_tf++) {
      if(V_light_Tag_fixed_L[iens_tf].substr(1,1) != "B") {
	amu_W_eps_tm_two_finest.distr_list.push_back( amu_W_eps_physical_point_tm_red[id_eps].distr_list[iens_tf]);
	amu_W_eps_OS_two_finest.distr_list.push_back( amu_W_eps_physical_point_OS_red[id_eps].distr_list[iens_tf]);
	L_list_fixed_L_two_finest.push_back(L_list_fixed_L[iens_tf]);
	a_distr_list_fixed_L_two_finest.distr_list.push_back( a_distr_list_fixed_L.distr_list[iens_tf]);
	Mpi_fit_fixed_L_two_finest.distr_list.push_back( Mpi_fit_fixed_L.distr_list[iens_tf]);
	fp_fit_fixed_L_two_finest.distr_list.push_back( fp_fit_fixed_L.distr_list[iens_tf]);
	V_light_Tag_fixed_L_two_finest.push_back( V_light_Tag_fixed_L[iens_tf]);
      }
    }

    Perform_Akaike_fits(amu_W_eps_tm_two_finest, amu_W_eps_OS_two_finest, a_A, a_B, a_C, a_D, L_list_fixed_L_two_finest, a_distr_list_fixed_L_two_finest, Mpi_fit_fixed_L_two_finest,fp_fit_fixed_L_two_finest, V_light_Tag_fixed_L_two_finest, UseJack, Njacks, Nboots, "W_win_fixed_L_a2_two_finest_eps_"+to_string_with_precision(eps_win[id_eps]/fm_to_inv_Gev,2), "light",a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list,allow_a4_and_log,allow_only_finest, 1, false, 0, 0, 200.0, LL, 0.0, return_val, return4_val);

    /*
    single_fit_list={"tm"};
    a4_list={"tm"};
    a2_list={"tm"};
    Perform_Akaike_fits(amu_W_eps_physical_point_tm_red[id_eps], amu_W_eps_physical_point_OS_red[id_eps], a_A, a_B, a_C, a_D, L_list_fixed_L, a_distr_list_fixed_L, Mpi_fit_fixed_L,fp_fit_fixed_L, V_light_Tag_fixed_L, UseJack, Njacks, Nboots, "W_win_fixed_L_tm_eps_"+to_string_with_precision(eps_win[id_eps]/fm_to_inv_Gev,2), "light",a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list,allow_a4_and_log,allow_only_finest, 1, false, 0, 0, 200.0, LL, 0.0, return_val, return4_val);

    single_fit_list={"OS"};
    a4_list={"OS"};
    a2_list={"OS"};
    Perform_Akaike_fits(amu_W_eps_physical_point_tm_red[id_eps], amu_W_eps_physical_point_OS_red[id_eps], a_A, a_B, a_C, a_D, L_list_fixed_L, a_distr_list_fixed_L, Mpi_fit_fixed_L,fp_fit_fixed_L, V_light_Tag_fixed_L, UseJack, Njacks, Nboots, "W_win_fixed_L_tm_eps_"+to_string_with_precision(eps_win[id_eps]/fm_to_inv_Gev,2), "light",a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list,allow_a4_and_log,allow_only_finest, 1, false, 0, 0, 200.0, LL, 0.0, return_val, return4_val);
    */

    
    

  }


  //##########################    FINAL EXTRAPOLATION OF LIGHT INTERMEDIATE-WINDOW CONTRIBUTION ################################
  cout<<"Continuum extrapolation of standard Intermediate window"<<endl;
  //continuum extrapolation of the standard intermediate window
  distr_t return_val;
  
  //n_m_pair_list= { make_pair(0,0), make_pair(1,1), make_pair(2,2), make_pair(3,3)};
  mass_extr_list={"off"};
  FSEs_list = {"off"};
  allow_only_finest=false;
  allow_a4_and_log=true;
  single_fit_list = {"tm"};
  a2_list = {"tm"};
  a4_list = {"off"};
  n_m_pair_list = { make_pair(0,0), make_pair(1,0), make_pair(2,0), make_pair(3,0) };
  //single fit tm
  Perform_Akaike_fits(amu_W_physical_point_tm_red, amu_W_physical_point_OS_red, a_A, a_B, a_C, a_D, L_list_fixed_L, a_distr_list_fixed_L, Mpi_fit_fixed_L,fp_fit_fixed_L, V_light_Tag_fixed_L, UseJack, Njacks, Nboots, "W_win_final_single_tm", "light",a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list,allow_a4_and_log,allow_only_finest, 2, false, 0, 0, 200.0, LL, 0.0, return_val, return4_val);

  single_fit_list={"OS"};
  a2_list= {"OS"};
  n_m_pair_list = { make_pair(0,0), make_pair(0,1), make_pair(0,2), make_pair(0,3) };

  //single fit OS
  Perform_Akaike_fits(amu_W_physical_point_tm_red, amu_W_physical_point_OS_red, a_A, a_B, a_C, a_D, L_list_fixed_L, a_distr_list_fixed_L, Mpi_fit_fixed_L,fp_fit_fixed_L, V_light_Tag_fixed_L, UseJack, Njacks, Nboots, "W_win_final_single_OS", "light",a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list,allow_a4_and_log,allow_only_finest, 2, false, 0, 0, 200.0, LL, 0.0, return_val, return4_val);

  n_m_pair_list = {make_pair(0,0), make_pair(3,0), make_pair(0,3), make_pair(3,3), make_pair(1,0), make_pair(0,1), make_pair(1,1), make_pair(2,0), make_pair(0,2), make_pair(2,2)};
  //combined fit standard
  single_fit_list={"off"};
  a2_list={"on"};
  a4_list={"on", "off", "tm", "OS"};
  Perform_Akaike_fits(amu_W_physical_point_tm_red, amu_W_physical_point_OS_red, a_A, a_B, a_C, a_D, L_list_fixed_L, a_distr_list_fixed_L, Mpi_fit_fixed_L,fp_fit_fixed_L, V_light_Tag_fixed_L, UseJack, Njacks, Nboots, "W_win_final", "light",a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list,allow_a4_and_log,allow_only_finest, 2, true, 0, 0, 200.0, LL, 0.0, return_val, return4_val);


  //linear a^2 leaving 1 or two measurements 
  a4_list={"off"};
  Perform_Akaike_fits(amu_W_physical_point_tm_red, amu_W_physical_point_OS_red, a_A, a_B, a_C, a_D, L_list_fixed_L, a_distr_list_fixed_L, Mpi_fit_fixed_L,fp_fit_fixed_L, V_light_Tag_fixed_L, UseJack, Njacks, Nboots, "W_win_final_leave_tm_B", "light",a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list,allow_a4_and_log,allow_only_finest, 2, false, 0, 0, 200.0, LL, 0.0, return_val, return4_val);
  Perform_Akaike_fits(amu_W_physical_point_tm_red, amu_W_physical_point_OS_red, a_A, a_B, a_C, a_D, L_list_fixed_L, a_distr_list_fixed_L, Mpi_fit_fixed_L,fp_fit_fixed_L, V_light_Tag_fixed_L, UseJack, Njacks, Nboots, "W_win_final_leave_OS_B", "light",a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list,allow_a4_and_log,allow_only_finest, 2, false, 0, 0, 200.0, LL, 0.0, return_val, return4_val);
  Perform_Akaike_fits(amu_W_physical_point_tm_red, amu_W_physical_point_OS_red, a_A, a_B, a_C, a_D, L_list_fixed_L, a_distr_list_fixed_L, Mpi_fit_fixed_L,fp_fit_fixed_L, V_light_Tag_fixed_L, UseJack, Njacks, Nboots, "W_win_final_leave_B", "light",a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list,allow_a4_and_log,allow_only_finest, 2, false, 0, 0, 200.0, LL, 0.0, return_val, return4_val);
  


  //Pade'
  a4_list = {"tm", "OS"};
  //Perform_Akaike_fits(amu_W_physical_point_tm_red, amu_W_physical_point_OS_red, a_A, a_B, a_C, a_D, L_list_fixed_L, a_distr_list_fixed_L, Mpi_fit_fixed_L,fp_fit_fixed_L, V_light_Tag_fixed_L, UseJack, Njacks, Nboots, "W_win_final_Pade", "light",a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list,allow_a4_and_log,allow_only_finest, 2, false, 0, 0, 200.0, LL, 0.0, return_val, return4_val);


  //Ratio/Difference method
  a4_list = {"on", "off"};
  Perform_Akaike_fits(amu_W_physical_point_tm_red-amu_W_physical_point_OS_red, 1.0e-10*amu_W_physical_point_tm_red/amu_W_physical_point_OS_red, a_A, a_B, a_C, a_D, L_list_fixed_L, a_distr_list_fixed_L, Mpi_fit_fixed_L,fp_fit_fixed_L, V_light_Tag_fixed_L, UseJack, Njacks, Nboots, "W_win_final_ratio_diff_tm", "light",a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list,allow_a4_and_log,allow_only_finest, 2, true, 0, 0, 200.0, LL, 0.0, return_val, return4_val);
  Perform_Akaike_fits(amu_W_physical_point_tm_red-amu_W_physical_point_OS_red, 1.0e-10*amu_W_physical_point_OS_red/amu_W_physical_point_tm_red, a_A, a_B, a_C, a_D, L_list_fixed_L, a_distr_list_fixed_L, Mpi_fit_fixed_L,fp_fit_fixed_L, V_light_Tag_fixed_L, UseJack, Njacks, Nboots, "W_win_final_ratio_diff_OS", "light",a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list,allow_a4_and_log,allow_only_finest, 2, true, 0, 0, 200.0, LL, 0.0, return_val, return4_val);


  //##########################################################################################################################################################################################




  //####################### FINAL CONTINUUM EXTRAPOLATION OF SHORT-DISTANCE LIGHT #####################################################à
  

  mass_extr_list={"off"};
  FSEs_list = {"off"};
  allow_only_finest=false;
  allow_a4_and_log=true;
  single_fit_list = {"tm"};
  a2_list = {"tm"};
  a4_list = {"off"};
  n_m_pair_list = { make_pair(0,0), make_pair(1,0), make_pair(2,0), make_pair(3,0) };
  //single fit tm
  Perform_Akaike_fits(amu_SD_physical_point_tm_red, amu_SD_physical_point_OS_red, a_A, a_B, a_C, a_D, L_list_fixed_L, a_distr_list_fixed_L, Mpi_fit_fixed_L,fp_fit_fixed_L, V_light_Tag_fixed_L, UseJack, Njacks, Nboots, "SD_win_final_single_tm", "light",a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list,allow_a4_and_log,allow_only_finest, 2, false, 0, 0, 48, LL, 0.0, return_val, return4_val);

  single_fit_list={"OS"};
  a2_list= {"OS"};
  n_m_pair_list = { make_pair(0,0), make_pair(0,1), make_pair(0,2), make_pair(0,3) };

  //single fit OS
  Perform_Akaike_fits(amu_SD_physical_point_tm_red, amu_SD_physical_point_OS_red, a_A, a_B, a_C, a_D, L_list_fixed_L, a_distr_list_fixed_L, Mpi_fit_fixed_L,fp_fit_fixed_L, V_light_Tag_fixed_L, UseJack, Njacks, Nboots, "SD_win_final_single_OS", "light",a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list,allow_a4_and_log,allow_only_finest, 2, false, 0, 0, 48, LL, 0.0, return_val, return4_val);


  //combined fit standard
  single_fit_list={"off"};
  a2_list={"on"};
  a4_list={"on", "off", "tm", "OS"};
  n_m_pair_list = {make_pair(0,0), make_pair(3,0), make_pair(0,3), make_pair(3,3), make_pair(1,0), make_pair(0,1), make_pair(1,1), make_pair(2,0), make_pair(0,2), make_pair(2,2)};
  Perform_Akaike_fits(amu_SD_physical_point_tm_red, amu_SD_physical_point_OS_red, a_A, a_B, a_C, a_D, L_list_fixed_L, a_distr_list_fixed_L, Mpi_fit_fixed_L,fp_fit_fixed_L, V_light_Tag_fixed_L, UseJack, Njacks, Nboots, "SD_win_final", "light",a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list,allow_a4_and_log,allow_only_finest, 2, true, 0, 0, 48, LL, 0.0, return_val, return4_val);


  //linear a^2 leaving 1 or two measurements 
  a4_list={"off"};
  Perform_Akaike_fits(amu_SD_physical_point_tm_red, amu_SD_physical_point_OS_red, a_A, a_B, a_C, a_D, L_list_fixed_L, a_distr_list_fixed_L, Mpi_fit_fixed_L,fp_fit_fixed_L, V_light_Tag_fixed_L, UseJack, Njacks, Nboots, "SD_win_final_leave_tm_B", "light",a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list,allow_a4_and_log,allow_only_finest, 2, false, 0, 0, 48, LL, 0.0, return_val, return4_val);
  Perform_Akaike_fits(amu_SD_physical_point_tm_red, amu_SD_physical_point_OS_red, a_A, a_B, a_C, a_D, L_list_fixed_L, a_distr_list_fixed_L, Mpi_fit_fixed_L,fp_fit_fixed_L, V_light_Tag_fixed_L, UseJack, Njacks, Nboots, "SD_win_final_leave_OS_B", "light",a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list,allow_a4_and_log,allow_only_finest, 2, false, 0, 0, 48, LL, 0.0, return_val, return4_val);
  Perform_Akaike_fits(amu_SD_physical_point_tm_red, amu_SD_physical_point_OS_red, a_A, a_B, a_C, a_D, L_list_fixed_L, a_distr_list_fixed_L, Mpi_fit_fixed_L,fp_fit_fixed_L, V_light_Tag_fixed_L, UseJack, Njacks, Nboots, "SD_win_final_leave_B", "light",a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list,allow_a4_and_log,allow_only_finest, 2, false, 0, 0, 48, LL, 0.0, return_val, return4_val);
  


  //Pade'
  a4_list = {"tm", "OS"};
  //Perform_Akaike_fits(amu_SD_physical_point_tm_red, amu_SD_physical_point_OS_red, a_A, a_B, a_C, a_D, L_list_fixed_L, a_distr_list_fixed_L, Mpi_fit_fixed_L,fp_fit_fixed_L, V_light_Tag_fixed_L, UseJack, Njacks, Nboots, "SD_win_final_Pade", "light",a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list,allow_a4_and_log,allow_only_finest, 2, false, 0, 0, 48, LL, 0.0, return_val, return4_val);


  //Ratio/Difference method
  a4_list = {"on", "off"};
  Perform_Akaike_fits(amu_SD_physical_point_tm_red-amu_SD_physical_point_OS_red, 1.0e-10*amu_SD_physical_point_tm_red/amu_SD_physical_point_OS_red, a_A, a_B, a_C, a_D, L_list_fixed_L, a_distr_list_fixed_L, Mpi_fit_fixed_L,fp_fit_fixed_L, V_light_Tag_fixed_L, UseJack, Njacks, Nboots, "SD_win_final_ratio_diff_tm", "light",a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list,allow_a4_and_log,allow_only_finest, 2, true, 0, 0, 48, LL, 0.0, return_val, return4_val);
  Perform_Akaike_fits(amu_SD_physical_point_tm_red-amu_SD_physical_point_OS_red, 1.0e-10*amu_SD_physical_point_OS_red/amu_SD_physical_point_tm_red, a_A, a_B, a_C, a_D, L_list_fixed_L, a_distr_list_fixed_L, Mpi_fit_fixed_L,fp_fit_fixed_L, V_light_Tag_fixed_L, UseJack, Njacks, Nboots, "SD_win_final_ratio_diff_OS", "light",a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list,allow_a4_and_log,allow_only_finest, 2, true, 0, 0, 48, LL, 0.0, return_val, return4_val);



  //####################################################################################################################################




  //CONTINUUM EXTRAPOLATION OF SD-TMINS ################################################################################################

  vector<distr_t> ret_distr_SD_tmins;

  for(int tmins_id=0;tmins_id<(signed)tmins.size(); tmins_id++) {

    ret_distr_SD_tmins.emplace_back(UseJack);
    
    mass_extr_list={"off"};
    FSEs_list = {"off"};
    allow_only_finest=false;
    allow_a4_and_log=true;
    single_fit_list = {"tm"};
    a2_list = {"tm"};
    a4_list = {"off"};
    n_m_pair_list = { make_pair(0,0), make_pair(1,0), make_pair(2,0), make_pair(3,0) };
    //single fit tm
    Perform_Akaike_fits(amu_SD_tmins_physical_point_tm_red[tmins_id], amu_SD_tmins_physical_point_OS_red[tmins_id], a_A, a_B, a_C, a_D, L_list_fixed_L, a_distr_list_fixed_L, Mpi_fit_fixed_L,fp_fit_fixed_L, V_light_Tag_fixed_L, UseJack, Njacks, Nboots, "SD_win_final_single_tm_tmins_"+to_string_with_precision(tmins[tmins_id],4), "light",a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list,allow_a4_and_log,allow_only_finest, 2, false, 0, 0, 48, LL, 0.0, ret_distr_SD_tmins[tmins_id], return4_val);

    single_fit_list={"OS"};
    a2_list= {"OS"};
    n_m_pair_list = { make_pair(0,0), make_pair(0,1), make_pair(0,2), make_pair(0,3) };

    //single fit OS
    Perform_Akaike_fits(amu_SD_tmins_physical_point_tm_red[tmins_id], amu_SD_tmins_physical_point_OS_red[tmins_id], a_A, a_B, a_C, a_D, L_list_fixed_L, a_distr_list_fixed_L, Mpi_fit_fixed_L,fp_fit_fixed_L, V_light_Tag_fixed_L, UseJack, Njacks, Nboots, "SD_win_final_single_OS_tmins_"+to_string_with_precision(tmins[tmins_id],4), "light",a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list,allow_a4_and_log,allow_only_finest, 2, false, 0, 0, 48, LL, 0.0, ret_distr_SD_tmins[tmins_id], return4_val);


    //combined fit standard
    single_fit_list={"off"};
    a2_list={"on"};
    a4_list={"on", "off", "tm", "OS"};
    n_m_pair_list = {make_pair(0,0), make_pair(3,0), make_pair(0,3), make_pair(3,3), make_pair(1,0), make_pair(0,1), make_pair(1,1), make_pair(2,0), make_pair(0,2), make_pair(2,2)};
    Perform_Akaike_fits(amu_SD_tmins_physical_point_tm_red[tmins_id], amu_SD_tmins_physical_point_OS_red[tmins_id], a_A, a_B, a_C, a_D, L_list_fixed_L, a_distr_list_fixed_L, Mpi_fit_fixed_L,fp_fit_fixed_L, V_light_Tag_fixed_L, UseJack, Njacks, Nboots, "SD_win_final_tmins_"+to_string_with_precision(tmins[tmins_id],4), "light",a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list,allow_a4_and_log,allow_only_finest, 2, true, 0, 0, 48, LL, 0.0, ret_distr_SD_tmins[tmins_id], return4_val);


    //linear a^2 leaving 1 or two measurements 
    a4_list={"off"};
    Perform_Akaike_fits(amu_SD_tmins_physical_point_tm_red[tmins_id], amu_SD_tmins_physical_point_OS_red[tmins_id], a_A, a_B, a_C, a_D, L_list_fixed_L, a_distr_list_fixed_L, Mpi_fit_fixed_L,fp_fit_fixed_L, V_light_Tag_fixed_L, UseJack, Njacks, Nboots, "SD_win_final_leave_tm_B_tmins_"+to_string_with_precision(tmins[tmins_id],4), "light",a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list,allow_a4_and_log,allow_only_finest, 2, false, 0, 0, 48, LL, 0.0, ret_distr_SD_tmins[tmins_id], return4_val);
    Perform_Akaike_fits(amu_SD_tmins_physical_point_tm_red[tmins_id], amu_SD_tmins_physical_point_OS_red[tmins_id], a_A, a_B, a_C, a_D, L_list_fixed_L, a_distr_list_fixed_L, Mpi_fit_fixed_L,fp_fit_fixed_L, V_light_Tag_fixed_L, UseJack, Njacks, Nboots, "SD_win_final_leave_OS_B_tmins_"+to_string_with_precision(tmins[tmins_id],4), "light",a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list,allow_a4_and_log,allow_only_finest, 2, false, 0, 0, 48, LL, 0.0, ret_distr_SD_tmins[tmins_id], return4_val);
    Perform_Akaike_fits(amu_SD_tmins_physical_point_tm_red[tmins_id], amu_SD_tmins_physical_point_OS_red[tmins_id], a_A, a_B, a_C, a_D, L_list_fixed_L, a_distr_list_fixed_L, Mpi_fit_fixed_L,fp_fit_fixed_L, V_light_Tag_fixed_L, UseJack, Njacks, Nboots, "SD_win_final_leave_B_tmins_"+to_string_with_precision(tmins[tmins_id],4), "light",a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list,allow_a4_and_log,allow_only_finest, 2, false, 0, 0, 48, LL, 0.0, ret_distr_SD_tmins[tmins_id], return4_val);
  



    //Ratio/Difference method
    a4_list = {"on", "off"};
    Perform_Akaike_fits(amu_SD_tmins_physical_point_tm_red[tmins_id]-amu_SD_tmins_physical_point_OS_red[tmins_id], 1.0e-10*amu_SD_tmins_physical_point_tm_red[tmins_id]/amu_SD_tmins_physical_point_OS_red[tmins_id], a_A, a_B, a_C, a_D, L_list_fixed_L, a_distr_list_fixed_L, Mpi_fit_fixed_L,fp_fit_fixed_L, V_light_Tag_fixed_L, UseJack, Njacks, Nboots, "SD_win_final_ratio_diff_tm_tmins_"+to_string_with_precision(tmins[tmins_id],4), "light",a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list,allow_a4_and_log,allow_only_finest, 2, true, 0, 0, 48, LL, 0.0, ret_distr_SD_tmins[tmins_id], return4_val);
    Perform_Akaike_fits(amu_SD_tmins_physical_point_tm_red[tmins_id]-amu_SD_tmins_physical_point_OS_red[tmins_id], 1.0e-10*amu_SD_tmins_physical_point_OS_red[tmins_id]/amu_SD_tmins_physical_point_tm_red[tmins_id], a_A, a_B, a_C, a_D, L_list_fixed_L, a_distr_list_fixed_L, Mpi_fit_fixed_L,fp_fit_fixed_L, V_light_Tag_fixed_L, UseJack, Njacks, Nboots, "SD_win_final_ratio_diff_OS_tmins_"+to_string_with_precision(tmins[tmins_id],4), "light",a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list,allow_a4_and_log,allow_only_finest, 2, true, 0, 0, 48, LL, 0.0, ret_distr_SD_tmins[tmins_id], return4_val);


  }

  distr_t_list ret_distr_list_SD_tmins(UseJack);
  for(int ttt=0; ttt<(signed)tmins.size(); ttt++) ret_distr_list_SD_tmins.distr_list.push_back( ret_distr_SD_tmins[ttt]);

  //print tmins data
  Print_To_File({}, {tmins, ret_distr_list_SD_tmins.ave(), ret_distr_list_SD_tmins.err()  }, "../data/gm2/light/windows_fit_func/tmins_SD_light.dat" , "", "#tmins[fm] val err");

  //perform extrapolation at tmin=0

  class ipar_tmin {

  public:
    ipar_tmin() {}


    double val, err, tmin;

  };

  class fpar_tmin {

  public:
    fpar_tmin() {}
    fpar_tmin(const Vfloat &par)  {
    if((signed)par.size() != 3) crash("In class fpar_tmin constructor fpar_tmin(vector<double>) called with vector.size != 3 ");
    w0=par[0];
    D2=par[1];
    D4=par[2];
    }

    double w0, D2, D4;

  };

  bootstrap_fit<fpar_tmin, ipar_tmin> bf_tmin(Njacks);
  bf_tmin.set_warmup_lev(0);
  bf_tmin.Set_number_of_measurements(tmins.size());
  bf_tmin.Set_verbosity(1);
  bf_tmin.Add_par("w0", 48.0, 0.5);
  bf_tmin.Add_par("D2", 1.0, 0.1);
  bf_tmin.Add_par("D4", 1.0, 0.1);
  bf_tmin.Fix_par("D4", 0.0);
  //fit on mean values to get ch2
  bootstrap_fit<fpar_tmin, ipar_tmin> bf_ch2_tmin(1);
  bf_ch2_tmin.set_warmup_lev(0);
  bf_ch2_tmin.Set_number_of_measurements(tmins.size());
  bf_ch2_tmin.Set_verbosity(1);
  bf_ch2_tmin.Add_par("w0", 48.0, 0.5);
  bf_ch2_tmin.Add_par("D2", 1.0, 0.1);
  bf_ch2_tmin.Add_par("D4", 1.0, 0.1);
  bf_ch2_tmin.Fix_par("D4", 0.0);

  double w0_g=48.0;

  //ansatz
  bf_tmin.ansatz=  [&w0_g ](const fpar_tmin &p, const ipar_tmin &ip) {

		  
		     return p.w0 + p.D2*w0_g*pow(ip.tmin*0.3,2)/pow(log(ip.tmin*0.3),1) + p.D4*w0_g*pow(ip.tmin*0.3,4);
		   }; 
  bf_tmin.measurement=  [ ](const fpar_tmin &p, const ipar_tmin &ip) {

		 return ip.val;
		 };
  bf_tmin.error=  [ ](const fpar_tmin &p, const ipar_tmin &ip) {

		 return ip.err;
		 };
  

  
  bf_ch2_tmin.ansatz= bf_tmin.ansatz;
  bf_ch2_tmin.measurement = bf_tmin.measurement;
  bf_ch2_tmin.error = bf_tmin.error;


  vector<vector<ipar_tmin>> data_tmin(Njacks);
  vector<vector<ipar_tmin>> data_ch2_tmin(1);
  //allocate space for output result
  boot_fit_data<fpar_tmin> Bt_fit_tmin;
  boot_fit_data<fpar_tmin> Bt_fit_ch2_tmin;
  for(auto &data_iboot: data_tmin) data_iboot.resize(tmins.size());
  for(auto &data_iboot: data_ch2_tmin) data_iboot.resize(tmins.size());


  //covariance matrix
  Eigen::MatrixXd Cov_Matrix_tmin(tmins.size(),tmins.size());
  Eigen::MatrixXd Corr_Matrix_tmin(tmins.size(), tmins.size());
  for(int i=0;i<(signed)tmins.size();i++)
    for(int j=0;j<(signed)tmins.size();j++) {
      Cov_Matrix_tmin(i,j) = ret_distr_SD_tmins[i]%ret_distr_SD_tmins[j];
      Corr_Matrix_tmin(i,j) = Cov_Matrix_tmin(i,j)/(ret_distr_SD_tmins[i].err()*ret_distr_SD_tmins[j].err());
    }
  //bf_tmin.Add_covariance_matrix(Cov_Matrix_tmin);
  //bf_ch2_tmin.Add_covariance_matrix(Cov_Matrix_tmin);


  Eigen::MatrixXd Corr_Matrix_inverse_tmin = Corr_Matrix_tmin.inverse();

  ofstream Print_cov_tmin("../data/gm2/Covariance_matrix/SD_light_tmins.cov");
  ofstream Print_corr_tmin("../data/gm2/Correlation_matrix/SD_light_tmins.cov");
  ofstream Print_corr_inv_tmin("../data/gm2/Correlation_matrix_inverse/SD_light_tmins.cov");
  Print_cov_tmin<<Cov_Matrix_tmin<<endl;
  Print_corr_tmin<<Corr_Matrix_tmin<<endl;
  Print_corr_inv_tmin<<Corr_Matrix_inverse_tmin<<endl;
  Print_cov_tmin.close();
  Print_corr_tmin.close();
  Print_corr_inv_tmin.close();

  for(int ijack=0;ijack<Njacks;ijack++) {
    for(int tmin_id=0;tmin_id<(signed)tmins.size();tmin_id++) {

      data_tmin[ijack][tmin_id].val = ret_distr_SD_tmins[tmin_id].distr[ijack];
      data_tmin[ijack][tmin_id].err = ret_distr_SD_tmins[tmin_id].err();
      data_tmin[ijack][tmin_id].tmin = tmins[tmin_id]*fm_to_inv_Gev;

      if(ijack==0) {

	data_ch2_tmin[ijack][tmin_id].val = ret_distr_SD_tmins[tmin_id].ave();
	data_ch2_tmin[ijack][tmin_id].err = ret_distr_SD_tmins[tmin_id].err();
	data_ch2_tmin[ijack][tmin_id].tmin = tmins[tmin_id]*fm_to_inv_Gev;
      }
    }
  }

  //append
  bf_tmin.Append_to_input_par(data_tmin);
  bf_ch2_tmin.Append_to_input_par(data_ch2_tmin);
  //fit
  cout<<"Fitting tmins... "<<endl;
  Bt_fit_tmin= bf_tmin.Perform_bootstrap_fit();
  Bt_fit_ch2_tmin= bf_ch2_tmin.Perform_bootstrap_fit();

  //retrieve parameters
  distr_t w0_tmin(UseJack), D2_tmin(UseJack), D4_tmin(UseJack);
  for(int ijack=0;ijack<Njacks;ijack++) { w0_tmin.distr.push_back( Bt_fit_tmin.par[ijack].w0); D2_tmin.distr.push_back( Bt_fit_tmin.par[ijack].D2); D4_tmin.distr.push_back( Bt_fit_tmin.par[ijack].D4); }
  //push_back ch2
  double ch2_tmin= Bt_fit_ch2_tmin.get_ch2_ave();

 
  //print fit func
  distr_t_list SD_tmins_to_print(UseJack);
  Vfloat tmin_list_to_print(400);
  for(int tt=0;tt<400;tt++) tmin_list_to_print[tt] = tt*0.2/(400-1.0);
  for(auto &tmin: tmin_list_to_print) {
    
    SD_tmins_to_print.distr_list.push_back( w0_tmin + D2_tmin*w0_g*pow(tmin*fm_to_inv_Gev*0.3,2)/pow(log(tmin*fm_to_inv_Gev*0.3),1) + D4_tmin*w0_g*pow(tmin*fm_to_inv_Gev*0.3,4) );

      }
  Print_To_File({}, {tmin_list_to_print, SD_tmins_to_print.ave(), SD_tmins_to_print.err()}, "../data/gm2/light/windows_fit_func/tmins.fit_func", "", "#tmin[fm] SD SD_err, ch2= "+to_string_with_precision(ch2_tmin,5));
  
  
  //######################################################################################################################################









  //####################################################################################################################################



  distr_t_list amu_strange_W_tm_red(UseJack), amu_strange_W_OS_red(UseJack);
  distr_t_list amu_strange_SD_tm_red(UseJack), amu_strange_SD_OS_red(UseJack);
  Vfloat L_strange_list_red;
  distr_t_list a_distr_list_strange_red(UseJack);
  distr_t_list Mpi_fit_strange_red(UseJack);
  distr_t_list fp_fit_strange_red(UseJack);
  vector<string> V_strange_tag_red;
  for(int ist=0;ist<(signed)V_strange_1_L.Tag.size();ist++) {
    if(V_strange_1_L.Tag[ist] != "cB211b.072.64") {
      if(V_strange_1_L.Tag[ist] != "cB211b.072.96") {
      amu_strange_W_tm_red.distr_list.push_back( agm2_strange_W_Extr.distr_list[ist]);
      amu_strange_W_OS_red.distr_list.push_back( agm2_strange_W_OS_Extr.distr_list[ist]);
      amu_strange_SD_tm_red.distr_list.push_back( agm2_strange_SD_Extr.distr_list[ist]);
      amu_strange_SD_OS_red.distr_list.push_back( agm2_strange_SD_OS_Extr.distr_list[ist]);
      L_strange_list_red.push_back( L_strange_list[ist]);
      a_distr_list_strange_red.distr_list.push_back( a_distr_list_strange.distr_list[ist]);
      Mpi_fit_strange_red.distr_list.push_back(Mpi_fit.distr_list[ist]);
      fp_fit_strange_red.distr_list.push_back(fp_fit.distr_list[ist]);
      V_strange_tag_red.push_back(V_strange_1_L.Tag[ist]);
      }
      else {
	//find B64 tag
	bool find_B64=false;
	int B64_id=0;
	for(int ib=0; ib<(signed)V_strange_1_L.Tag.size();ib++) {
	  if(V_strange_1_L.Tag[ib] == "cB211b.072.64") { find_B64=true; B64_id= ib; break;}
	}
	if(!find_B64) crash("Cannot find ensemble B64 when computing strange window contribution");

	double weight_W_tm_B64= (1.0/pow( agm2_strange_W_Extr.err(B64_id),2))/( (1.0/pow(agm2_strange_W_Extr.err(B64_id),2)) + (1.0/pow(agm2_strange_W_Extr.err(ist),2)));
	double weight_W_tm_B96= (1.0/pow( agm2_strange_W_Extr.err(ist),2))/( (1.0/pow(agm2_strange_W_Extr.err(B64_id),2)) + (1.0/pow(agm2_strange_W_Extr.err(ist),2)));

	double weight_W_OS_B64= (1.0/pow( agm2_strange_W_OS_Extr.err(B64_id),2))/( (1.0/pow(agm2_strange_W_OS_Extr.err(B64_id),2)) + (1.0/pow(agm2_strange_W_OS_Extr.err(ist),2)));
	double weight_W_OS_B96= (1.0/pow( agm2_strange_W_OS_Extr.err(ist),2))/( (1.0/pow(agm2_strange_W_OS_Extr.err(B64_id),2)) + (1.0/pow(agm2_strange_W_OS_Extr.err(ist),2)));

	double weight_SD_tm_B64= (1.0/pow( agm2_strange_SD_Extr.err(B64_id),2))/( (1.0/pow(agm2_strange_SD_Extr.err(B64_id),2)) + (1.0/pow(agm2_strange_SD_Extr.err(ist),2)));
	double weight_SD_tm_B96= (1.0/pow( agm2_strange_SD_Extr.err(ist),2))/( (1.0/pow(agm2_strange_SD_Extr.err(B64_id),2)) + (1.0/pow(agm2_strange_SD_Extr.err(ist),2)));

	double weight_SD_OS_B64= (1.0/pow( agm2_strange_SD_OS_Extr.err(B64_id),2))/( (1.0/pow(agm2_strange_SD_OS_Extr.err(B64_id),2)) + (1.0/pow(agm2_strange_SD_OS_Extr.err(ist),2)));
	double weight_SD_OS_B96= (1.0/pow( agm2_strange_SD_OS_Extr.err(ist),2))/( (1.0/pow(agm2_strange_SD_OS_Extr.err(B64_id),2)) + (1.0/pow(agm2_strange_SD_OS_Extr.err(ist),2)));

	//push_back
	//W
	amu_strange_W_tm_red.distr_list.push_back( weight_W_tm_B64*agm2_strange_W_Extr.distr_list[B64_id] + weight_W_tm_B96*agm2_strange_W_Extr.distr_list[ist]);
	amu_strange_W_OS_red.distr_list.push_back( weight_W_OS_B64*agm2_strange_W_OS_Extr.distr_list[B64_id] + weight_W_OS_B96*agm2_strange_W_OS_Extr.distr_list[ist]);

	//SD
	amu_strange_SD_tm_red.distr_list.push_back( weight_SD_tm_B64*agm2_strange_SD_Extr.distr_list[B64_id] + weight_SD_tm_B96*agm2_strange_SD_Extr.distr_list[ist]);
	amu_strange_SD_OS_red.distr_list.push_back( weight_SD_OS_B64*agm2_strange_SD_OS_Extr.distr_list[B64_id] + weight_SD_OS_B96*agm2_strange_SD_OS_Extr.distr_list[ist]);
	
	L_strange_list_red.push_back( L_strange_list[ist]);
	a_distr_list_strange_red.distr_list.push_back( a_distr_list_strange.distr_list[ist]);
	Mpi_fit_strange_red.distr_list.push_back(Mpi_fit.distr_list[ist]);
	fp_fit_strange_red.distr_list.push_back(fp_fit.distr_list[ist]);
	V_strange_tag_red.push_back(V_strange_1_L.Tag[ist]);
      }
    }
  }


  //tm
  Print_To_File(V_strange_tag_red, { L_strange_list_red, (a_distr_list_strange_red/fm_to_inv_Gev).ave() , amu_strange_W_tm_red.ave(), amu_strange_W_tm_red.err()}  , "../data/gm2/strange/tm_"+Extrapolation_strange_mode+"/windows_tm_red", "", "#Ens L a amuW");
  Print_To_File(V_strange_tag_red, { L_strange_list_red, (a_distr_list_strange_red/fm_to_inv_Gev).ave() , amu_strange_SD_tm_red.ave(), amu_strange_SD_tm_red.err()}  , "../data/gm2/strange/tm_"+Extrapolation_strange_mode+"/SD_tm_red", "", "#Ens L a amuSD");

  //OS
  Print_To_File(V_strange_tag_red, { L_strange_list_red, (a_distr_list_strange_red/fm_to_inv_Gev).ave() , amu_strange_W_OS_red.ave(), amu_strange_W_OS_red.err()}  , "../data/gm2/strange/OS_"+Extrapolation_strange_mode+"/windows_OS_red", "", "#Ens L a amuW");
  Print_To_File(V_strange_tag_red, { L_strange_list_red, (a_distr_list_strange_red/fm_to_inv_Gev).ave() , amu_strange_SD_OS_red.ave(), amu_strange_SD_OS_red.err()}  , "../data/gm2/strange/OS_"+Extrapolation_strange_mode+"/SD_OS_red", "", "#Ens L a amuSD");

  //FINAL CONTINUUM EXTRAPOLATION OF INTERMEDIATE DISTANCE STRANGE #####################################################################
 

  
  //n_m_pair_list= { make_pair(0,0), make_pair(1,1), make_pair(2,2), make_pair(3,3)};
  mass_extr_list={"off"};
  FSEs_list = {"off"};
  allow_only_finest=false;
  allow_a4_and_log=true;
  single_fit_list = {"tm"};
  a2_list = {"tm"};
  a4_list = {"off"};
  n_m_pair_list = { make_pair(0,0), make_pair(1,0), make_pair(2,0), make_pair(3,0) };
  //single fit tm
  Perform_Akaike_fits(amu_strange_W_tm_red, amu_strange_W_OS_red, a_A, a_B, a_C, a_D, L_strange_list_red, a_distr_list_strange_red, Mpi_fit_strange_red,fp_fit_strange_red, V_strange_tag_red, UseJack, Njacks, Nboots, Extrapolation_strange_mode+"_W_win_final_single_tm", "strange",a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list,allow_a4_and_log,allow_only_finest, 2, false, 0, 0, 27.2, LL, 0.0, return_val, return4_val);

  single_fit_list={"OS"};
  a2_list= {"OS"};
  n_m_pair_list = { make_pair(0,0), make_pair(0,1), make_pair(0,2), make_pair(0,3) };

  //single fit OS
  Perform_Akaike_fits(amu_strange_W_tm_red, amu_strange_W_OS_red, a_A, a_B, a_C, a_D, L_strange_list_red, a_distr_list_strange_red, Mpi_fit_strange_red,fp_fit_strange_red, V_strange_tag_red, UseJack, Njacks, Nboots, Extrapolation_strange_mode+"_W_win_final_single_OS", "strange",a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list,allow_a4_and_log,allow_only_finest, 2, false, 0, 0, 27.2, LL, 0.0, return_val, return4_val);


  //combined fit standard
  single_fit_list={"off"};
  a2_list={"on"};
  a4_list={"on", "off", "tm", "OS"};
  n_m_pair_list = {make_pair(0,0), make_pair(3,0), make_pair(0,3), make_pair(3,3), make_pair(1,0), make_pair(0,1), make_pair(1,1), make_pair(2,0), make_pair(0,2), make_pair(2,2)};
  Perform_Akaike_fits(amu_strange_W_tm_red, amu_strange_W_OS_red, a_A, a_B, a_C, a_D, L_strange_list_red, a_distr_list_strange_red, Mpi_fit_strange_red,fp_fit_strange_red, V_strange_tag_red, UseJack, Njacks, Nboots, Extrapolation_strange_mode+"_W_win_final", "strange",a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list,allow_a4_and_log,allow_only_finest, 2, true, 0, 0, 27.2, LL, 0.0, return_val, return4_val);


  //linear a^2 leaving 1 or two measurements 
  a4_list={"off"};
  Perform_Akaike_fits(amu_strange_W_tm_red, amu_strange_W_OS_red, a_A, a_B, a_C, a_D, L_strange_list_red, a_distr_list_strange_red, Mpi_fit_strange_red,fp_fit_strange_red, V_strange_tag_red, UseJack, Njacks, Nboots, Extrapolation_strange_mode+"_W_win_final_leave_tm_B", "strange",a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list,allow_a4_and_log,allow_only_finest, 2, false, 0, 0, 27.2, LL, 0.0, return_val, return4_val);
  Perform_Akaike_fits(amu_strange_W_tm_red, amu_strange_W_OS_red, a_A, a_B, a_C, a_D, L_strange_list_red, a_distr_list_strange_red, Mpi_fit_strange_red,fp_fit_strange_red, V_strange_tag_red, UseJack, Njacks, Nboots, Extrapolation_strange_mode+"_W_win_final_leave_OS_B", "strange",a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list,allow_a4_and_log,allow_only_finest, 2, false, 0, 0, 27.2, LL, 0.0, return_val, return4_val);
  Perform_Akaike_fits(amu_strange_W_tm_red, amu_strange_W_OS_red, a_A, a_B, a_C, a_D, L_strange_list_red, a_distr_list_strange_red, Mpi_fit_strange_red,fp_fit_strange_red, V_strange_tag_red, UseJack, Njacks, Nboots, Extrapolation_strange_mode+"_W_win_final_leave_B", "strange",a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list,allow_a4_and_log,allow_only_finest, 2, false, 0, 0, 27.2, LL, 0.0, return_val, return4_val);
  


  //Pade'
  a4_list = {"tm", "OS"};
  //Perform_Akaike_fits(amu_strange_W_tm_red, amu_strange_W_OS_red, a_A, a_B, a_C, a_D, L_strange_list_red, a_distr_list_strange_red, Mpi_fit_strange_red,fp_fit_strange_red, V_strange_tag_red, UseJack, Njacks, Nboots, Extrapolation_strange_mode+"_W_win_final_Pade", "strange",a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list,allow_a4_and_log,allow_only_finest, 2, false, 0, 0, 27.2, LL, 0.0, return_val, return4_val);


  //Ratio/Difference method
  a4_list = {"on", "off"};
  Perform_Akaike_fits(amu_strange_W_tm_red-amu_strange_W_OS_red, 1.0e-10*amu_strange_W_tm_red/amu_strange_W_OS_red, a_A, a_B, a_C, a_D, L_strange_list_red, a_distr_list_strange_red, Mpi_fit_strange_red,fp_fit_strange_red, V_strange_tag_red, UseJack, Njacks, Nboots, Extrapolation_strange_mode+"_W_win_final_ratio_diff_tm", "strange",a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list,allow_a4_and_log,allow_only_finest, 2, true, 0, 0, 27.2, LL, 0.0, return_val, return4_val);
  Perform_Akaike_fits(amu_strange_W_tm_red-amu_strange_W_OS_red, 1.0e-10*amu_strange_W_OS_red/amu_strange_W_tm_red, a_A, a_B, a_C, a_D, L_strange_list_red, a_distr_list_strange_red, Mpi_fit_strange_red,fp_fit_strange_red, V_strange_tag_red, UseJack, Njacks, Nboots, Extrapolation_strange_mode+"_W_win_final_ratio_diff_OS", "strange",a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list,allow_a4_and_log,allow_only_finest, 2, true, 0, 0, 27.2, LL, 0.0, return_val, return4_val);
  

  //####################################################################################################################################





   

  //FINAL CONTINUUM EXTRAPOLATION OF SHORT-DISTANCE STRANGE #####################################################################
  //n_m_pair_list= { make_pair(0,0), make_pair(1,1), make_pair(2,2), make_pair(3,3)};
  mass_extr_list={"off"};
  FSEs_list = {"off"};
  allow_only_finest=false;
  allow_a4_and_log=true;
  single_fit_list = {"tm"};
  a2_list = {"tm"};
  a4_list = {"off"};
  n_m_pair_list = { make_pair(0,0), make_pair(1,0), make_pair(2,0), make_pair(3,0) };
  //single fit tm
  Perform_Akaike_fits(amu_strange_SD_tm_red, amu_strange_SD_OS_red, a_A, a_B, a_C, a_D, L_strange_list_red, a_distr_list_strange_red, Mpi_fit_strange_red,fp_fit_strange_red, V_strange_tag_red, UseJack, Njacks, Nboots, Extrapolation_strange_mode+"_SD_win_final_single_tm", "strange",a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list,allow_a4_and_log,allow_only_finest, 2, false, 0, 0, 8, LL, 0.0, return_val, return4_val);

  single_fit_list={"OS"};
  a2_list= {"OS"};
  n_m_pair_list = { make_pair(0,0), make_pair(0,1), make_pair(0,2), make_pair(0,3) };

  //single fit OS
  Perform_Akaike_fits(amu_strange_SD_tm_red, amu_strange_SD_OS_red, a_A, a_B, a_C, a_D, L_strange_list_red, a_distr_list_strange_red, Mpi_fit_strange_red,fp_fit_strange_red, V_strange_tag_red, UseJack, Njacks, Nboots, Extrapolation_strange_mode+"_SD_win_final_single_OS", "strange",a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list,allow_a4_and_log,allow_only_finest, 2, false, 0, 0, 8, LL, 0.0, return_val, return4_val);


  //combined fit standard
  single_fit_list={"off"};
  a2_list={"on"};
  a4_list={"on", "off", "tm", "OS"};
  n_m_pair_list = {make_pair(0,0), make_pair(3,0), make_pair(0,3), make_pair(3,3), make_pair(1,0), make_pair(0,1), make_pair(1,1), make_pair(2,0), make_pair(0,2), make_pair(2,2)};
  Perform_Akaike_fits(amu_strange_SD_tm_red, amu_strange_SD_OS_red, a_A, a_B, a_C, a_D, L_strange_list_red, a_distr_list_strange_red, Mpi_fit_strange_red,fp_fit_strange_red, V_strange_tag_red, UseJack, Njacks, Nboots, Extrapolation_strange_mode+"_SD_win_final", "strange",a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list,allow_a4_and_log,allow_only_finest, 2, true, 0, 0, 8, LL, 0.0, return_val, return4_val);


  //linear a^2 leaving 1 or two measurements 
  a4_list={"off"};
  Perform_Akaike_fits(amu_strange_SD_tm_red, amu_strange_SD_OS_red, a_A, a_B, a_C, a_D, L_strange_list_red, a_distr_list_strange_red, Mpi_fit_strange_red,fp_fit_strange_red, V_strange_tag_red, UseJack, Njacks, Nboots, Extrapolation_strange_mode+"_SD_win_final_leave_tm_B", "strange",a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list,allow_a4_and_log,allow_only_finest, 2, false, 0, 0, 8, LL, 0.0, return_val, return4_val);
  Perform_Akaike_fits(amu_strange_SD_tm_red, amu_strange_SD_OS_red, a_A, a_B, a_C, a_D, L_strange_list_red, a_distr_list_strange_red, Mpi_fit_strange_red,fp_fit_strange_red, V_strange_tag_red, UseJack, Njacks, Nboots, Extrapolation_strange_mode+"_SD_win_final_leave_OS_B", "strange",a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list,allow_a4_and_log,allow_only_finest, 2, false, 0, 0, 8, LL, 0.0, return_val, return4_val);
  Perform_Akaike_fits(amu_strange_SD_tm_red, amu_strange_SD_OS_red, a_A, a_B, a_C, a_D, L_strange_list_red, a_distr_list_strange_red, Mpi_fit_strange_red,fp_fit_strange_red, V_strange_tag_red, UseJack, Njacks, Nboots, Extrapolation_strange_mode+"_SD_win_final_leave_B", "strange",a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list,allow_a4_and_log,allow_only_finest, 2, false, 0, 0, 8, LL, 0.0, return_val, return4_val);
  


  //Pade'
  a4_list = {"tm", "OS"};
  //Perform_Akaike_fits(amu_strange_SD_tm_red, amu_strange_SD_OS_red, a_A, a_B, a_C, a_D, L_strange_list_red, a_distr_list_strange_red, Mpi_fit_strange_red,fp_fit_strange_red, V_strange_tag_red, UseJack, Njacks, Nboots, Extrapolation_strange_mode+"_SD_win_final_Pade", "strange",a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list,allow_a4_and_log,allow_only_finest, 2, false, 0, 0, 8, LL, 0.0, return_val, return4_val);


  //Ratio/Difference method
  a4_list = {"on", "off"};
  Perform_Akaike_fits(amu_strange_SD_tm_red-amu_strange_SD_OS_red, 1.0e-10*amu_strange_SD_tm_red/amu_strange_SD_OS_red, a_A, a_B, a_C, a_D, L_strange_list_red, a_distr_list_strange_red, Mpi_fit_strange_red,fp_fit_strange_red, V_strange_tag_red, UseJack, Njacks, Nboots, Extrapolation_strange_mode+"_SD_win_final_ratio_diff_tm", "strange",a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list,allow_a4_and_log,allow_only_finest, 2, true, 0, 0, 8, LL, 0.0, return_val, return4_val);
  Perform_Akaike_fits(amu_strange_SD_tm_red-amu_strange_SD_OS_red, 1.0e-10*amu_strange_SD_OS_red/amu_strange_SD_tm_red, a_A, a_B, a_C, a_D, L_strange_list_red, a_distr_list_strange_red, Mpi_fit_strange_red,fp_fit_strange_red, V_strange_tag_red, UseJack, Njacks, Nboots, Extrapolation_strange_mode+"_SD_win_final_ratio_diff_OS", "strange",a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list,allow_a4_and_log,allow_only_finest, 2, true, 0, 0, 8, LL, 0.0, return_val, return4_val);
  

  //####################################################################################################################################



  //##############################################   CHARM #############################################################################



  bool first_occurrence_A_ens=true;
  distr_t_list amu_charm_W_tm_red(UseJack), amu_charm_W_OS_red(UseJack);
  distr_t_list amu_charm_SD_tm_red(UseJack), amu_charm_SD_OS_red(UseJack);
  Vfloat L_charm_list_red;
  distr_t_list a_distr_list_charm_red(UseJack);
  distr_t_list Mpi_fit_charm_red(UseJack);
  distr_t_list fp_fit_charm_red(UseJack);
  vector<string> V_charm_tag_red;
  for(int ist=0;ist<(signed)V_charm_1_L.Tag.size();ist++) {
    if(V_charm_1_L.Tag[ist].substr(1,1) != "A") {
      amu_charm_W_tm_red.distr_list.push_back( agm2_charm_W_Extr.distr_list[ist]);
      amu_charm_W_OS_red.distr_list.push_back( agm2_charm_W_OS_Extr.distr_list[ist]);
      amu_charm_SD_tm_red.distr_list.push_back( agm2_charm_SD_Extr.distr_list[ist]);
      amu_charm_SD_OS_red.distr_list.push_back( agm2_charm_SD_OS_Extr.distr_list[ist]);
      L_charm_list_red.push_back( L_charm_list[ist]);
      a_distr_list_charm_red.distr_list.push_back( a_distr_list_charm.distr_list[ist]);
      Mpi_fit_charm_red.distr_list.push_back(Mpi_fit_charm.distr_list[ist]);
      fp_fit_charm_red.distr_list.push_back(fp_fit_charm.distr_list[ist]);
      V_charm_tag_red.push_back(V_charm_1_L.Tag[ist]);
    }
    else  {
      if(first_occurrence_A_ens) {
	//loop over A ensemble
	double mult_A_W_tm_ens=0;
	double mult_A_W_OS_ens=0;
	double mult_A_SD_tm_ens=0;
	double mult_A_SD_OS_ens=0;
	distr_t W_tm(UseJack,UseJack?Njacks:Nboots), W_OS(UseJack,UseJack?Njacks:Nboots), SD_OS(UseJack,UseJack?Njacks:Nboots), SD_tm(UseJack,UseJack?Njacks:Nboots); //constructor sets to 0 
	for(int isA=0; isA<(signed)V_charm_1_L.Tag.size();isA++) {
	  if(V_charm_1_L.Tag[isA].substr(1,1) == "A") {
	    W_tm = W_tm + (1.0/pow(agm2_charm_W_Extr.err(isA),2))*agm2_charm_W_Extr.distr_list[isA];
	    W_OS = W_OS + (1.0/pow(agm2_charm_W_OS_Extr.err(isA),2))*agm2_charm_W_OS_Extr.distr_list[isA];
	    SD_tm = SD_tm + (1.0/pow(agm2_charm_SD_Extr.err(isA),2))*agm2_charm_SD_Extr.distr_list[isA];
	    SD_OS = SD_OS + (1.0/pow(agm2_charm_SD_OS_Extr.err(isA),2))*agm2_charm_SD_OS_Extr.distr_list[isA];
	    mult_A_W_tm_ens += 1.0/pow(agm2_charm_W_Extr.err(isA),2);
	    mult_A_W_OS_ens += 1.0/pow(agm2_charm_W_OS_Extr.err(isA),2);
	    mult_A_SD_tm_ens += 1.0/pow(agm2_charm_SD_Extr.err(isA),2);
	    mult_A_SD_OS_ens += 1.0/pow(agm2_charm_SD_OS_Extr.err(isA),2);
	  }
	}

	//push_back
	amu_charm_W_tm_red.distr_list.push_back(W_tm/mult_A_W_tm_ens );
	amu_charm_W_OS_red.distr_list.push_back( W_OS/mult_A_W_OS_ens);
	amu_charm_SD_tm_red.distr_list.push_back( SD_tm/mult_A_SD_tm_ens);
	amu_charm_SD_OS_red.distr_list.push_back( SD_OS/mult_A_SD_OS_ens);
	L_charm_list_red.push_back( L_charm_list[ist]);
	a_distr_list_charm_red.distr_list.push_back( a_distr_list_charm.distr_list[ist]);
	Mpi_fit_charm_red.distr_list.push_back(Mpi_fit_charm.distr_list[ist]);
	fp_fit_charm_red.distr_list.push_back(fp_fit_charm.distr_list[ist]);
	V_charm_tag_red.push_back(V_charm_1_L.Tag[ist]);
	
	
	first_occurrence_A_ens=false;
      }
    }
  }


  cout<<"Charm fit-data produced"<<endl;




  //tm
  Print_To_File(V_charm_tag_red, { L_charm_list_red, (a_distr_list_charm_red/fm_to_inv_Gev).ave() , amu_charm_W_tm_red.ave(), amu_charm_W_tm_red.err()}  , "../data/gm2/charm/tm_"+Extrapolation_charm_mode+"/windows_tm_red", "", "#Ens L a amuW");
  Print_To_File(V_charm_tag_red, { L_charm_list_red, (a_distr_list_charm_red/fm_to_inv_Gev).ave() , amu_charm_SD_tm_red.ave(), amu_charm_SD_tm_red.err()}  , "../data/gm2/charm/tm_"+Extrapolation_charm_mode+"/SD_tm_red", "", "#Ens L a amuSD");

  //OS
  Print_To_File(V_charm_tag_red, { L_charm_list_red, (a_distr_list_charm_red/fm_to_inv_Gev).ave() , amu_charm_W_OS_red.ave(), amu_charm_W_OS_red.err()}  , "../data/gm2/charm/OS_"+Extrapolation_charm_mode+"/windows_OS_red", "", "#Ens L a amuW");
  Print_To_File(V_charm_tag_red, { L_charm_list_red, (a_distr_list_charm_red/fm_to_inv_Gev).ave() , amu_charm_SD_OS_red.ave(), amu_charm_SD_OS_red.err()}  , "../data/gm2/charm/OS_"+Extrapolation_charm_mode+"/SD_OS_red", "", "#Ens L a amuSD");
		

  
  //FINAL CONTINUUM EXTRAPOLATION OF INTERMEDIATE WINDOW CHARM

  //n_m_pair_list= { make_pair(0,0), make_pair(1,1), make_pair(2,2), make_pair(3,3)};
  mass_extr_list={"off"};
  FSEs_list = {"off"};
  allow_only_finest=false;
  allow_a4_and_log=true;
  single_fit_list = {"tm"};
  a2_list = { "tm"};
  a4_list = {"off", "tm"};
  n_m_pair_list = { make_pair(0,0), make_pair(1,0), make_pair(2,0), make_pair(3,0) };
  //single fit tm
  Perform_Akaike_fits(amu_charm_W_tm_red, amu_charm_W_OS_red, a_A, a_B, a_C, a_D, L_charm_list_red, a_distr_list_charm_red, Mpi_fit_charm_red,fp_fit_charm_red, V_charm_tag_red, UseJack, Njacks, Nboots, Extrapolation_charm_mode+"_W_win_final_single_tm", "charm",a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list,allow_a4_and_log,allow_only_finest, 2, true, 0, 0, 2.9, LL, 0.0, return_val, return4_val);

  single_fit_list={"OS"};
  a2_list= { "OS"};
  a4_list= {"off", "OS"};
  n_m_pair_list = { make_pair(0,0), make_pair(0,1), make_pair(0,2), make_pair(0,3) };

  //single fit OS
  Perform_Akaike_fits(amu_charm_W_tm_red, amu_charm_W_OS_red, a_A, a_B, a_C, a_D, L_charm_list_red, a_distr_list_charm_red, Mpi_fit_charm_red,fp_fit_charm_red, V_charm_tag_red, UseJack, Njacks, Nboots, Extrapolation_charm_mode+"_W_win_final_single_OS", "charm",a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list,allow_a4_and_log,allow_only_finest, 2, true, 0, 0, 2.9, LL, 0.0, return_val, return4_val);


  //combined fit standard
  single_fit_list={"off"};
  a2_list={"on"};
  a4_list={"on", "off", "tm", "OS"};
  n_m_pair_list = {make_pair(0,0), make_pair(3,0), make_pair(0,3), make_pair(3,3), make_pair(1,0), make_pair(0,1), make_pair(1,1), make_pair(2,0), make_pair(0,2), make_pair(2,2)};
  Perform_Akaike_fits(amu_charm_W_tm_red, amu_charm_W_OS_red, a_A, a_B, a_C, a_D, L_charm_list_red, a_distr_list_charm_red, Mpi_fit_charm_red,fp_fit_charm_red, V_charm_tag_red, UseJack, Njacks, Nboots, Extrapolation_charm_mode+"_W_win_final", "charm",a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list,allow_a4_and_log,allow_only_finest, 2, true, 0, 0, 2.9, LL, 0.0, return_val, return4_val);


  //linear a^2 leaving 1 or two measurements 
  a4_list={"off", "on", "tm", "OS"};
  Perform_Akaike_fits(amu_charm_W_tm_red, amu_charm_W_OS_red, a_A, a_B, a_C, a_D, L_charm_list_red, a_distr_list_charm_red, Mpi_fit_charm_red,fp_fit_charm_red, V_charm_tag_red, UseJack, Njacks, Nboots, Extrapolation_charm_mode+"_W_win_final_leave_tm_A", "charm",a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list,allow_a4_and_log,allow_only_finest, 2, true, 0, 0, 2.9, LL, 0.0, return_val, return4_val);
  Perform_Akaike_fits(amu_charm_W_tm_red, amu_charm_W_OS_red, a_A, a_B, a_C, a_D, L_charm_list_red, a_distr_list_charm_red, Mpi_fit_charm_red,fp_fit_charm_red, V_charm_tag_red, UseJack, Njacks, Nboots, Extrapolation_charm_mode+"_W_win_final_leave_OS_A", "charm",a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list,allow_a4_and_log,allow_only_finest, 2, true, 0, 0, 2.9, LL, 0.0, return_val, return4_val);
  Perform_Akaike_fits(amu_charm_W_tm_red, amu_charm_W_OS_red, a_A, a_B, a_C, a_D, L_charm_list_red, a_distr_list_charm_red, Mpi_fit_charm_red,fp_fit_charm_red, V_charm_tag_red, UseJack, Njacks, Nboots, Extrapolation_charm_mode+"_W_win_final_leave_A", "charm",a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list,allow_a4_and_log,allow_only_finest, 2, true, 0, 0, 2.9, LL, 0.0, return_val, return4_val);
  


  //Pade'
  //a4_list = {"tm", "OS"};
  //Perform_Akaike_fits(amu_charm_W_tm_red, amu_charm_W_OS_red, a_A, a_B, a_C, a_D, L_charm_list_red, a_distr_list_charm_red, Mpi_fit_charm_red,fp_fit_charm_red, V_charm_tag_red, UseJack, Njacks, Nboots, Extrapolation_charm_mode+"_W_win_final_Pade", "charm",a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list,allow_a4_and_log,allow_only_finest, 2, false, 0, 0, 2.9, LL, 0.0, return_val, return4_val);


  //Ratio/Difference method
  a4_list = {"on", "off"};
  Perform_Akaike_fits(amu_charm_W_tm_red-amu_charm_W_OS_red, 1.0e-10*amu_charm_W_tm_red/amu_charm_W_OS_red, a_A, a_B, a_C, a_D, L_charm_list_red, a_distr_list_charm_red, Mpi_fit_charm_red,fp_fit_charm_red, V_charm_tag_red, UseJack, Njacks, Nboots, Extrapolation_charm_mode+"_W_win_final_ratio_diff_tm", "charm",a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list,allow_a4_and_log,allow_only_finest, 2, true, 0, 0, 2.9, LL, 0.0, return_val, return4_val);
  Perform_Akaike_fits(amu_charm_W_tm_red-amu_charm_W_OS_red, 1.0e-10*amu_charm_W_OS_red/amu_charm_W_tm_red, a_A, a_B, a_C, a_D, L_charm_list_red, a_distr_list_charm_red, Mpi_fit_charm_red,fp_fit_charm_red, V_charm_tag_red, UseJack, Njacks, Nboots, Extrapolation_charm_mode+"_W_win_final_ratio_diff_OS", "charm",a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list,allow_a4_and_log,allow_only_finest, 2, true, 0, 0, 2.9, LL, 0.0, return_val, return4_val);



  //FINAL CONTINUUM EXTRAPOLATION OF SHORT-DISTANCE WINDOW CHARM
  

 
  //n_m_pair_list= { make_pair(0,0), make_pair(1,1), make_pair(2,2), make_pair(3,3)};
  mass_extr_list={"off"};
  FSEs_list = {"off"};
  allow_only_finest=false;
  allow_a4_and_log=true;
  single_fit_list = {"tm"};
  a2_list = {"tm"};
  a4_list = {"off", "tm"};
  n_m_pair_list = { make_pair(0,0), make_pair(1,0), make_pair(2,0), make_pair(3,0) };
  //single fit tm
  Perform_Akaike_fits(amu_charm_SD_tm_red, amu_charm_SD_OS_red, a_A, a_B, a_C, a_D, L_charm_list_red, a_distr_list_charm_red, Mpi_fit_charm_red,fp_fit_charm_red, V_charm_tag_red, UseJack, Njacks, Nboots, Extrapolation_charm_mode+"_SD_win_final_single_tm", "charm",a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list,allow_a4_and_log,allow_only_finest, 2, true, 0, 0, 11.0, LL, 0.0, return_val, return4_val);

  single_fit_list={"OS"};
  a2_list= {"OS"};
  a4_list= {"off", "OS"};
  n_m_pair_list = { make_pair(0,0), make_pair(0,1), make_pair(0,2), make_pair(0,3) };

  //single fit OS
  Perform_Akaike_fits(amu_charm_SD_tm_red, amu_charm_SD_OS_red, a_A, a_B, a_C, a_D, L_charm_list_red, a_distr_list_charm_red, Mpi_fit_charm_red,fp_fit_charm_red, V_charm_tag_red, UseJack, Njacks, Nboots, Extrapolation_charm_mode+"_SD_win_final_single_OS", "charm",a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list,allow_a4_and_log,allow_only_finest, 2, true, 0, 0, 11.0, LL, 0.0, return_val, return4_val);


  //combined fit standard
  single_fit_list={"off"};
  a2_list={"on"};
  a4_list={"on", "off", "tm", "OS"};
  n_m_pair_list = {make_pair(0,0), make_pair(3,0), make_pair(0,3), make_pair(3,3), make_pair(1,0), make_pair(0,1), make_pair(1,1), make_pair(2,0), make_pair(0,2), make_pair(2,2)};
  Perform_Akaike_fits(amu_charm_SD_tm_red, amu_charm_SD_OS_red, a_A, a_B, a_C, a_D, L_charm_list_red, a_distr_list_charm_red, Mpi_fit_charm_red,fp_fit_charm_red, V_charm_tag_red, UseJack, Njacks, Nboots, Extrapolation_charm_mode+"_SD_win_final", "charm",a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list,allow_a4_and_log,allow_only_finest, 2, true, 0, 0, 11.0, LL, 0.0, return_val, return4_val);


  //linear a^2 leaving 1 or two measurements 
  a4_list={"off", "on", "tm", "OS"};
  Perform_Akaike_fits(amu_charm_SD_tm_red, amu_charm_SD_OS_red, a_A, a_B, a_C, a_D, L_charm_list_red, a_distr_list_charm_red, Mpi_fit_charm_red,fp_fit_charm_red, V_charm_tag_red, UseJack, Njacks, Nboots, Extrapolation_charm_mode+"_SD_win_final_leave_tm_A", "charm",a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list,allow_a4_and_log,allow_only_finest, 2, true, 0, 0, 11.0, LL, 0.0, return_val, return4_val);
  Perform_Akaike_fits(amu_charm_SD_tm_red, amu_charm_SD_OS_red, a_A, a_B, a_C, a_D, L_charm_list_red, a_distr_list_charm_red, Mpi_fit_charm_red,fp_fit_charm_red, V_charm_tag_red, UseJack, Njacks, Nboots, Extrapolation_charm_mode+"_SD_win_final_leave_OS_A", "charm",a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list,allow_a4_and_log,allow_only_finest, 2, true, 0, 0, 11.0, LL, 0.0, return_val, return4_val);
  Perform_Akaike_fits(amu_charm_SD_tm_red, amu_charm_SD_OS_red, a_A, a_B, a_C, a_D, L_charm_list_red, a_distr_list_charm_red, Mpi_fit_charm_red,fp_fit_charm_red, V_charm_tag_red, UseJack, Njacks, Nboots, Extrapolation_charm_mode+"_SD_win_final_leave_A", "charm",a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list,allow_a4_and_log,allow_only_finest, 2, true, 0, 0, 11.0, LL, 0.0, return_val, return4_val);
  


  //Pade'
  //a4_list = {"tm", "OS"};
  //Perform_Akaike_fits(amu_charm_SD_tm_red, amu_charm_SD_OS_red, a_A, a_B, a_C, a_D, L_charm_list_red, a_distr_list_charm_red, Mpi_fit_charm_red,fp_fit_charm_red, V_charm_tag_red, UseJack, Njacks, Nboots, Extrapolation_charm_mode+"_SD_win_final_Pade", "charm",a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list,allow_a4_and_log,allow_only_finest, 2, false, 0, 0, 11.0, LL, 0.0, return_val, return4_val);


  //Ratio/Difference method
  a4_list = {"on", "off"};
  Perform_Akaike_fits(amu_charm_SD_tm_red-amu_charm_SD_OS_red, 1.0e-10*amu_charm_SD_tm_red/amu_charm_SD_OS_red, a_A, a_B, a_C, a_D, L_charm_list_red, a_distr_list_charm_red, Mpi_fit_charm_red,fp_fit_charm_red, V_charm_tag_red, UseJack, Njacks, Nboots, Extrapolation_charm_mode+"_SD_win_final_ratio_diff_tm", "charm",a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list,allow_a4_and_log,allow_only_finest, 2, true, 0, 0, 11.0, LL, 0.0, return_val, return4_val);
  Perform_Akaike_fits(amu_charm_SD_tm_red-amu_charm_SD_OS_red, 1.0e-10*amu_charm_SD_OS_red/amu_charm_SD_tm_red, a_A, a_B, a_C, a_D, L_charm_list_red, a_distr_list_charm_red, Mpi_fit_charm_red,fp_fit_charm_red, V_charm_tag_red, UseJack, Njacks, Nboots, Extrapolation_charm_mode+"_SD_win_final_ratio_diff_OS", "charm",a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list,allow_a4_and_log,allow_only_finest, 2, true, 0, 0, 11.0, LL, 0.0, return_val, return4_val);
  





  //################################################################################################



  //disco


  distr_t ret_distr_SD_disco(UseJack), ret_distr_W_disco(UseJack);


  if(Include_light_disco) {


    single_fit_list={"OS"};
    FSEs_list={"off"};
    a4_list={"off"};
    n_m_pair_list ={make_pair(0,0)};
    a2_list={"off", "OS"};
    mass_extr_list={"off"};
    allow_a4_and_log = false;


    //W disco
    Perform_Akaike_fits(comb_disco_W, comb_disco_W, a_A, a_B, a_C, a_D, vols_disco, a_lat, Mpi_fit_disco, fp_fit_disco, comb_Tags, UseJack, Njacks, Nboots, "W_win_disco", "total_disco", a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list, allow_a4_and_log, allow_only_finest, 1, false, 0,0, -0.9, LL, 0.0, ret_distr_W_disco, return4_val);


    //SD disco
    Perform_Akaike_fits(comb_disco_SD, comb_disco_SD, a_A, a_B, a_C, a_D, vols_disco, a_lat, Mpi_fit_disco, fp_fit_disco, comb_Tags, UseJack, Njacks, Nboots, "SD_win_disco", "total_disco", a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list, allow_a4_and_log, allow_only_finest, 1, false,  0,0, -0.001, LL, 0.0, ret_distr_SD_disco, return4_val);



  }


  
  exit(-1);


  

  distr_t ret_distr_SD_light(UseJack), ret_distr_W_light(UseJack); 
  //full disco light
  //Perform_Akaike_fits(agm2_disco_light, agm2_disco_light, a_A, a_B, a_C, a_D, L_list_disco, a_distr_list_disco_light, Mpi_fit_disco, fp_fit_disco, disco_light_Tags, UseJack, Njacks, Nboots, "total_disco", "light", a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list, allow_a4_and_log, allow_only_finest,0,0, -10.0, LL, 0.0);


  allow_a4_and_log=true;
  mass_extr_list = {"on", "off"};
  single_fit_list = {"off"};
  a4_list = {"off", "on", "tm", "OS"};
  FSEs_list = {"off", "comb_GS"}; //add comb_GS
  a2_list = {"on"};

  //light connected
    

  //NO ELM
  //Perform_Akaike_fits(agm2_light_SD, agm2_light_SD_OS, a_A, a_B, a_C, a_D, L_list, a_distr_list, Mpi_fit,fp_fit, V_light_1.Tag, UseJack, Njacks, Nboots, "SD_win", "light",a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list,allow_a4_and_log,allow_only_finest, w0s_mult, true, 0, 0, 200.0, LL, 0.0, ret_distr_SD_light);

  /*
 
  vector<distr_t> ret_distr_SD_tmins;

  FSEs_list = {"off"};

  

  for(int tmins_id=0;tmins_id<(signed)tmins.size(); tmins_id++) {

    ret_distr_SD_tmins.emplace_back(UseJack);
    Perform_Akaike_fits(agm2_light_SD_tmins_distr_list[tmins_id], agm2_light_SD_OS_tmins_distr_list[tmins_id], a_A, a_B, a_C, a_D, L_list, a_distr_list, Mpi_fit,fp_fit, V_light_1.Tag, UseJack, Njacks, Nboots, "SD_win_tmins_"+to_string_with_precision(tmins[tmins_id],4), "light",a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list,allow_a4_and_log,allow_only_finest, w0s_mult, true, 0, 0, 200.0, LL, tmins[tmins_id]*fm_to_inv_Gev, ret_distr_SD_tmins[tmins_id]);

  }


  //Extrapolate distr_SD_tmins to zero
  vector<distr_t> tmins2_distr;
  for(int tt=0;tt<(signed)tmins.size();tt++) {
    tmins2_distr.emplace_back(UseJack);
    for(int ijack=0; ijack<Njacks;ijack++) {
      tmins2_distr[tt].distr.push_back( pow(tmins[tt],2)*(1.0 + GM()*1e-14));
    }
  }
  Obs_extrapolation_meson_mass(ret_distr_SD_tmins, tmins2_distr, 0.0, "../data/gm2/light", "tmins_cont_lim", UseJack, "FIT");
  
  */


 
  FSEs_list = {"comb_GS"};
  mass_extr_list= {"off"};
  a4_list = {"off", "tm", "OS", "on"};
  a2_list = {"on"};
  n_m_pair_list = {make_pair(0,0), make_pair(3,0), make_pair(0,3), make_pair(3,3), make_pair(1,0), make_pair(0,1), make_pair(1,1), make_pair(2,0), make_pair(0,2), make_pair(2,2)};
  //ELM
  //Perform_Akaike_fits(agm2_light_W_ELM, agm2_light_W_ELM_OS, a_A, a_B, a_C, a_D, L_list, a_distr_list, Mpi_fit,fp_fit, V_light_1.Tag, UseJack, Njacks, Nboots, "W_win_ELM", "light",a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list,allow_a4_and_log,allow_only_finest, 0, 0, 200.0, LL, 0.0);
  //NO ELM

  int id_64_ens=-1;
  /*
  //shift the data because of pion mass mistuning
  //determine cB64ens
  
  for(int i=0;i<(signed)V_light_1.Tag.size();i++) if(V_light_1.Tag[i]=="cB211b.072.64") id_64_ens=i;

  for(int i=0;i<(signed)V_light_1.Tag.size();i++) {

    if(V_light_1.Tag[i]=="cB211b.072.64") {
      agm2_light_W.distr_list[i] = agm2_light_W.distr_list[i]+ amu_W_mass_corr_tm_list.distr_list[0];
      agm2_light_W_OS.distr_list[i] = agm2_light_W_OS.distr_list[i]+ amu_W_mass_corr_OS_list.distr_list[0];
    }
    else {

      distr_t Mp_corr= (Mpi_fit.distr_list[i]/a_distr_list[i]-m_pi)/(Mpi_fit.distr_list[id_64_ens]/a_distr_list[id_64_ens] -m_pi);
      agm2_light_W.distr_list[i] = agm2_light_W.distr_list[i] + amu_W_mass_corr_tm_list.distr_list[0]*Mp_corr;
      agm2_light_W_OS.distr_list[i] = agm2_light_W_OS.distr_list[i] + amu_W_mass_corr_OS_list.distr_list[0]*Mp_corr;


    }
  }
  */
  
  //Perform_Akaike_fits(agm2_light_W, agm2_light_W_OS, a_A, a_B, a_C, a_D, L_list, a_distr_list, Mpi_fit,fp_fit, V_light_1.Tag, UseJack, Njacks, Nboots, "W_win", "light",a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list,allow_a4_and_log,allow_only_finest, w0s_mult, true, 0, 0, 200.0, LL, 0.0, ret_distr_W_light);

  //without B96

  distr_t_list agm2_light_W_red(UseJack), agm2_light_W_OS_red(UseJack), a_distr_list_red(UseJack), Mpi_fit_red(UseJack), fp_fit_red(UseJack);
  Vfloat L_list_red;
  vector<string> V_light_Tag_red;
  distr_t ret_distr_W_light_red(UseJack);

  for(int iens=0; iens<(signed)V_light_1.Tag.size(); iens++) {

    if(V_light_1.Tag[iens] != "cB211b.072.96") {

      //Interpolate B64 to L~ 5.46 fm using B96
      distr_t agm2_W_tm_interpol(UseJack);
      if(V_light_1.Tag[iens] == "cB211b.072.64") {
	int b96_id=0;
	bool Find_b96=false;
	for(int jens=0; jens<V_light_1.Tag.size(); jens++) {
	  if( V_light_1.Tag[jens] == "cB211b.072.96") { b96_id=jens; Find_b96=true; break;}
	}
	if(!Find_b96) crash("Cannot find ensemble cB211b.072.96 in W_red construction");

	for(int ijck=0; ijck< Njacks;ijck++) {
	  double exp_2_ML_B64= exp( -1.0*Mpi_fit.distr_list[iens].distr[ijck]*L_list[iens]);
	  double exp_2_ML_B96= exp(-1.0*Mpi_fit.distr_list[b96_id].distr[ijck]*L_list[b96_id]);
	  double exp_2_ML_extr= exp(-1.0*Mpi_fit.distr_list[iens].distr[ijck]*5.46*fm_to_inv_Gev/a_distr_list.distr_list[iens].distr[ijck]);
	  double amu_96= agm2_light_W.distr_list[b96_id].distr[ijck];
	  double amu_64 = agm2_light_W.distr_list[iens].distr[ijck];
	  double amu_Linf = (amu_96*exp_2_ML_B64 - amu_64*exp_2_ML_B96)/(exp_2_ML_B64 - exp_2_ML_B96);
	  double bmu_Linf = (amu_64 - amu_Linf)/exp_2_ML_B64;
	  agm2_W_tm_interpol.distr.push_back( amu_Linf+ bmu_Linf*exp_2_ML_extr);
	  
	}
      }

      
      
      if( V_light_1.Tag[iens] != "cB211b.072.64") agm2_light_W_red.distr_list.push_back( agm2_light_W.distr_list[iens]);
      else agm2_light_W_red.distr_list.push_back( agm2_W_tm_interpol);
      agm2_light_W_OS_red.distr_list.push_back( agm2_light_W_OS.distr_list[iens]);
      a_distr_list_red.distr_list.push_back( a_distr_list.distr_list[iens]);
      Mpi_fit_red.distr_list.push_back( Mpi_fit.distr_list[iens]);
      fp_fit_red.distr_list.push_back( fp_fit.distr_list[iens]);
      L_list_red.push_back( L_list[iens]);
      V_light_Tag_red.push_back( V_light_1.Tag[iens]);

    }
    

  }

  //Print To File
  Print_To_File(V_light_Tag_red,{a_distr_list_red.ave(), L_list_red, agm2_light_W_red.ave(), agm2_light_W_red.err()}, "../data/gm2/light/tm/windows_red.dat", "", "#Ens a  L    agm2");
  Print_To_File(V_light_Tag_red,{a_distr_list_red.ave(), L_list_red, agm2_light_W_OS_red.ave(), agm2_light_W_OS_red.err()}, "../data/gm2/light/OS/windows_red.dat", "", "#Ens a  L    agm2");


  //perform extrapolation at fixed volume
  FSEs_list = {"off"};



  //Perform_Akaike_fits(agm2_light_W_red, agm2_light_W_OS_red, a_A, a_B, a_C, a_D, L_list_red, a_distr_list_red, Mpi_fit_red,fp_fit_red, V_light_Tag_red, UseJack, Njacks, Nboots, "W_win_red", "light",a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list,allow_a4_and_log,allow_only_finest, w0s_mult, true, 0, 0, 200.0, LL, 0.0, ret_distr_W_light_red);

 
    
   
    


  //#####################################################################################################################
  //#####################################################################################################################
  //#####################################################################################################################






    

    


  //strange
  //#####################################################################################################################
  //#####################################################################################################################
  //#####################################################################################################################

  FSEs_list = {"off"};
  a4_list={"off", "on", "tm", "OS"};
  mass_extr_list={"off","on"};
  single_fit_list = {"off"};
  a2_list = {"on", "off", "tm", "OS"};


  distr_t ret_distr_SD_strange(UseJack), ret_distr_W_strange(UseJack);

   
  //W strange
  //ELM
  //Perform_Akaike_fits(agm2_strange_W_ELM_Extr, agm2_strange_W_ELM_OS_Extr, a_A, a_B, a_C, a_D, L_strange_list, a_distr_list_strange, Mpi_fit, fp_fit, V_strange_1_L.Tag, UseJack, Njacks, Nboots, "W_win_ELM_"+Extrapolation_strange_mode, "strange",a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list,allow_a4_and_log,allow_only_finest,0, 0, 27.0, LL, 0.0);
  //NO ELM
  Perform_Akaike_fits(agm2_strange_W_Extr, agm2_strange_W_OS_Extr, a_A, a_B, a_C, a_D, L_strange_list, a_distr_list_strange, Mpi_fit, fp_fit, V_strange_1_L.Tag, UseJack, Njacks, Nboots, "W_win_"+Extrapolation_strange_mode, "strange",a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list,allow_a4_and_log,allow_only_finest, w0s_mult, true, 0, 0, 27.0, LL, 0.0, ret_distr_W_strange, return4_val);


  n_m_pair_list = {make_pair(0,0), make_pair(3,0), make_pair(0,3), make_pair(3,3), make_pair(1,0), make_pair(0,1), make_pair(1,1), make_pair(2,0), make_pair(0,2), make_pair(2,2)};
    
  //SD strange
  //ELM
  //Perform_Akaike_fits(agm2_strange_SD_ELM_Extr, agm2_strange_SD_ELM_OS_Extr, a_A, a_B, a_C, a_D, L_strange_list, a_distr_list_strange, Mpi_fit, fp_fit, V_strange_1_L.Tag, UseJack, Njacks, Nboots, "SD_win_ELM_"+Extrapolation_strange_mode, "strange", a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list,allow_a4_and_log,allow_only_finest,0,0, 9.0, LL, 0.0);
  //NO ELM
  Perform_Akaike_fits(agm2_strange_SD_Extr, agm2_strange_SD_OS_Extr, a_A, a_B, a_C, a_D, L_strange_list, a_distr_list_strange, Mpi_fit, fp_fit, V_strange_1_L.Tag, UseJack, Njacks, Nboots, "SD_win_"+Extrapolation_strange_mode, "strange", a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list,allow_a4_and_log,allow_only_finest, w0s_mult, true, 0,0, 9.0, LL, 0.0, ret_distr_SD_strange, return4_val);


  n_m_pair_list = {make_pair(0,0), make_pair(3,0), make_pair(0,3), make_pair(3,3), make_pair(1,0), make_pair(0,1), make_pair(1,1), make_pair(2,0), make_pair(0,2), make_pair(2,2)};

  //total strange
  //ELM
  //Perform_Akaike_fits(agm2_strange_No_ELM_Extr, agm2_strange_OS_No_ELM_Extr, a_A, a_B, a_C, a_D, L_strange_list, a_distr_list_strange, Mpi_fit, fp_fit, V_strange_1_L.Tag, UseJack, Njacks, Nboots, "total_ELM_"+Extrapolation_strange_mode, "strange", a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list, allow_a4_and_log,allow_only_finest, 0,0, 50.0, LL, 0.0);
  //NO ELM
  //Perform_Akaike_fits(agm2_strange_Extr, agm2_strange_OS_Extr, a_A, a_B, a_C, a_D, L_strange_list, a_distr_list_strange, Mpi_fit, fp_fit, V_strange_1_L.Tag, UseJack, Njacks, Nboots, "total_"+Extrapolation_strange_mode, "strange", a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list, allow_a4_and_log,allow_only_finest, 0,0, 50.0, LL, 0.0);

  //#####################################################################################################################
  //#####################################################################################################################
  //#####################################################################################################################









    

    

  //charm
  //#####################################################################################################################
  //#####################################################################################################################
  //#####################################################################################################################


  a2_list ={"on", "tm", "OS"};
  FSEs_list ={"off"};
  mass_extr_list={"off"};
  single_fit_list={"off"};
  a4_list ={"on", "off", "tm", "OS"};
  allow_a4_and_log=true;
  allow_only_finest=true;

  n_m_pair_list = {make_pair(0,0), make_pair(3,0), make_pair(0,3), make_pair(3,3), make_pair(1,0), make_pair(0,1), make_pair(1,1), make_pair(2,0), make_pair(0,2), make_pair(2,2)};

  distr_t ret_distr_SD_charm(UseJack), ret_distr_W_charm(UseJack);

 
    
  //W charm
  //ELM
  //Perform_Akaike_fits(agm2_charm_W_ELM_Extr, agm2_charm_W_ELM_OS_Extr, a_A, a_B, a_C, a_D, L_charm_list, a_distr_list_charm, Mpi_fit_charm, fp_fit_charm, V_charm_1_L.Tag, UseJack, Njacks, Nboots, "W_win_ELM_"+Extrapolation_charm_mode, "charm", a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list, allow_a4_and_log, allow_only_finest, 0,0, 2.5, LL, 0.0);
  //NO ELM
  Perform_Akaike_fits(agm2_charm_W_Extr, agm2_charm_W_OS_Extr, a_A, a_B, a_C, a_D, L_charm_list, a_distr_list_charm, Mpi_fit_charm, fp_fit_charm, V_charm_1_L.Tag, UseJack, Njacks, Nboots, "W_win_"+Extrapolation_charm_mode, "charm", a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list, allow_a4_and_log,allow_only_finest, w0s_mult, true, 0,0, 2.5, LL, 0.0, ret_distr_W_charm, return4_val);


  n_m_pair_list = {make_pair(0,0), make_pair(3,0), make_pair(0,3), make_pair(3,3), make_pair(1,0), make_pair(0,1), make_pair(1,1), make_pair(2,0), make_pair(0,2), make_pair(2,2)};
  //SD charm
  //ELM
  //Perform_Akaike_fits(agm2_charm_SD_ELM_Extr, agm2_charm_SD_ELM_OS_Extr, a_A, a_B, a_C, a_D, L_charm_list, a_distr_list_charm, Mpi_fit_charm, fp_fit_charm, V_charm_1_L.Tag, UseJack, Njacks, Nboots, "SD_win_ELM_"+Extrapolation_charm_mode, "charm", a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list, allow_a4_and_log,allow_only_finest, 0,0, 11.5, LL, 0.0);
  //NO ELM

  a2_list={"on"};
  
  Perform_Akaike_fits(agm2_charm_SD_Extr, agm2_charm_SD_OS_Extr, a_A, a_B, a_C, a_D, L_charm_list, a_distr_list_charm, Mpi_fit_charm, fp_fit_charm, V_charm_1_L.Tag, UseJack, Njacks, Nboots, "SD_win_"+Extrapolation_charm_mode, "charm", a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list, allow_a4_and_log, allow_only_finest, w0s_mult, true, 0,0, 11.5, LL, 0.0, ret_distr_SD_charm, return4_val);

  

  //total charm
  //ELM
  //Perform_Akaike_fits(agm2_charm_Extr, agm2_charm_OS_Extr, a_A, a_B, a_C, a_D, L_charm_list, a_distr_list_charm, Mpi_fit_charm, fp_fit_charm, V_charm_1_L.Tag, UseJack, Njacks, Nboots, "total_ELM_"+Extrapolation_charm_mode, "charm", a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list, allow_a4_and_log,allow_only_finest, 0,0, 14.0, LL, 0.0);
  //NO ELM
  //Perform_Akaike_fits(agm2_charm_No_ELM_Extr, agm2_charm_OS_No_ELM_Extr, a_A, a_B, a_C, a_D, L_charm_list, a_distr_list_charm, Mpi_fit_charm, fp_fit_charm, V_charm_1_L.Tag, UseJack, Njacks, Nboots, "total_"+Extrapolation_charm_mode, "charm", a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list, allow_a4_and_log,allow_only_finest, 0,0, 14.0, LL, 0.0);
    

  allow_only_finest=false;
   

  
  //#####################################################################################################################
  //#####################################################################################################################
  //#####################################################################################################################


  
 







  //disco
    
  //#####################################################################################################################
  //#####################################################################################################################
  //#####################################################################################################################






  Print_To_File({}, { ret_distr_SD_light.distr, ret_distr_SD_strange.distr, ret_distr_SD_charm.distr, ret_distr_SD_disco.distr}  , "../data/gm2/fit_corr_SD_"+Extrapolation_strange_mode+"_"+Extrapolation_charm_mode+".dat" ,"", "#light strange charm disco");
  Print_To_File({}, { ret_distr_W_light.distr, ret_distr_W_strange.distr, ret_distr_W_charm.distr, ret_distr_W_disco.distr}  , "../data/gm2/fit_corr_W_"+Extrapolation_strange_mode+"_"+Extrapolation_charm_mode+".dat" ,"", "#light strange charm disco");



  exit(-1);


  //#####################################################################################################
  //continuum extrapolation of PI(Q^2)

  
  for(int j=0; j< (signed)Qs2.size();j++) {

    distr_t ret_distr_PI_Q2_light(UseJack), ret_distr_PI_Q2_light_ps(UseJack);
    allow_a4_and_log=true;
    allow_only_finest=false;
    mass_extr_list = {"on", "off"};
    single_fit_list = {"off"};
    a4_list = {"off", "on", "tm", "OS"};
    FSEs_list = {"off", "comb_GS"}; //add comb_GS
    a2_list = {"on", "off", "tm", "OS"};
    n_m_pair_list = {make_pair(0,0), make_pair(3,0), make_pair(0,3), make_pair(3,3), make_pair(1,0), make_pair(0,1), make_pair(1,1), make_pair(2,0), make_pair(0,2), make_pair(2,2)};
    //light 
    Perform_Akaike_fits_PI_Q2(Qs2[j], PI_Q2_light_tm[j], PI_Q2_light_OS[j], a_A, a_B, a_C, a_D, L_list, a_distr_list, Mpi_fit,fp_fit, V_light_1.Tag, UseJack, Njacks, Nboots, "Q2_"+to_string_with_precision(Qs2[j],5), "light", a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list,allow_a4_and_log,allow_only_finest, w0s_mult, true, 0, 0, 0.1, LL, ret_distr_PI_Q2_light);
    //light pert sub
    Perform_Akaike_fits_PI_Q2(Qs2[j],PI_Q2_light_tm_pert_sub[j], PI_Q2_light_OS_pert_sub[j], a_A, a_B, a_C, a_D, L_list, a_distr_list, Mpi_fit,fp_fit, V_light_1.Tag, UseJack, Njacks, Nboots, "Q2_"+to_string_with_precision(Qs2[j],5)+"_ps", "light", a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list,allow_a4_and_log,allow_only_finest, w0s_mult, true, 0, 0, 0.1, LL, ret_distr_PI_Q2_light_ps);

    distr_t ret_distr_PI_Q2_strange(UseJack), ret_distr_PI_Q2_strange_ps(UseJack);
    FSEs_list = {"off"};
    mass_extr_list = {"off"};
    //strange
    Perform_Akaike_fits_PI_Q2(Qs2[j], PI_Q2_strange_tm[j], PI_Q2_strange_OS[j], a_A, a_B, a_C, a_D, L_strange_list, a_distr_list_strange, Mpi_fit, fp_fit, V_strange_1_L.Tag, UseJack, Njacks, Nboots, "Q2_"+to_string_with_precision(Qs2[j],5)+"_"+Extrapolation_strange_mode, "strange",a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list,allow_a4_and_log,allow_only_finest, w0s_mult, true,0, 0, 0.05, LL,  ret_distr_PI_Q2_strange);
    //strange pert sub
    Perform_Akaike_fits_PI_Q2(Qs2[j],PI_Q2_strange_tm_pert_sub[j], PI_Q2_strange_OS_pert_sub[j], a_A, a_B, a_C, a_D, L_strange_list, a_distr_list_strange, Mpi_fit, fp_fit, V_strange_1_L.Tag, UseJack, Njacks, Nboots, "Q2_"+to_string_with_precision(Qs2[j],5)+"_"+Extrapolation_strange_mode+"_ps", "strange",a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list,allow_a4_and_log,allow_only_finest, w0s_mult,true, 0, 0, 0.05, LL,  ret_distr_PI_Q2_strange_ps);

    allow_only_finest=true;

    distr_t ret_distr_PI_Q2_charm(UseJack), ret_distr_PI_Q2_charm_ps(UseJack);
    //charm
    Perform_Akaike_fits_PI_Q2(Qs2[j],PI_Q2_charm_tm[j], PI_Q2_charm_OS[j], a_A, a_B, a_C, a_D, L_charm_list, a_distr_list_charm, Mpi_fit_charm, fp_fit_charm, V_charm_1_L.Tag, UseJack, Njacks, Nboots, "Q2_"+to_string_with_precision(Qs2[j],5)+"_"+Extrapolation_charm_mode, "charm", a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list, allow_a4_and_log,allow_only_finest, w0s_mult, true, 0,0, 0.01, LL,  ret_distr_PI_Q2_charm);
    //charm pert sub
    Perform_Akaike_fits_PI_Q2(Qs2[j],PI_Q2_charm_tm[j], PI_Q2_charm_OS[j], a_A, a_B, a_C, a_D, L_charm_list, a_distr_list_charm, Mpi_fit_charm, fp_fit_charm, V_charm_1_L.Tag, UseJack, Njacks, Nboots, "Q2_"+to_string_with_precision(Qs2[j],5)+"_"+Extrapolation_charm_mode+"_ps", "charm", a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list, allow_a4_and_log,allow_only_finest, w0s_mult, true, 0,0, 0.01, LL,  ret_distr_PI_Q2_charm_ps);


    single_fit_list={"OS"};
    FSEs_list={"off"};
    a4_list={"off"};
    n_m_pair_list ={make_pair(0,0)};
    a2_list={"off", "OS"};
    mass_extr_list={"off"};
    allow_a4_and_log = false;
    allow_only_finest= false;

    distr_t ret_distr_PI_Q2_disc(UseJack);

    //disconnected
    if(Include_charm_disco && Include_light_disco && Include_off_diagonal_disco && Include_strange_disco) {
      Perform_Akaike_fits_PI_Q2(Qs2[j],PI_Q2_disco[j], PI_Q2_disco[j], a_A, a_B, a_C, a_D, L_list_disco_PI_Q2, a_disc_PI_Q2, Mpi_fit_disc_PI_Q2, fpi_fit_disc_PI_Q2, Ens_list_disc_PI_Q2, UseJack, Njacks, Nboots, "Q2_"+to_string_with_precision(Qs2[j],5), "disco", a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list, allow_a4_and_log,allow_only_finest, 1, false, 0,0, -0.01, LL,  ret_distr_PI_Q2_disc);
    }
    
    
  }


  //######################################################################################################


  //print correlation matrix on screen

  vector<distr_t> SD_corr({ ret_distr_SD_light, ret_distr_SD_strange, ret_distr_SD_charm, ret_distr_SD_disco  });
  vector<distr_t> W_corr({ ret_distr_W_light, ret_distr_W_strange, ret_distr_W_charm, ret_distr_W_disco  });


  //print SD corr
  cout<<"Printing correlation matrix of single contributions for SD"<<endl;
  for(int i=0; i< (signed)SD_corr.size();i++) {
    
    for(int j=0; j<(signed)SD_corr.size(); j++) {

      double val= (SD_corr[i]%SD_corr[j])/( SD_corr[i].err()*SD_corr[j].err());
      cout<<val<<" ";

    }
    cout<<endl;
  }

  //print W corr
  cout<<"Printing correlation matrix of single contributions for W"<<endl;
  for(int i=0; i< (signed)W_corr.size();i++) {
    for(int j=0; j<(signed)W_corr.size(); j++) {

      double val= (W_corr[i]%W_corr[j])/( W_corr[i].err()*W_corr[j].err());
      cout<<val<<" ";

    }
    cout<<endl;
    
  }





  

  cout<<"Bye"<<endl;

    
  //#####################################################################################################################
  //#####################################################################################################################
  //#####################################################################################################################

   


  

 
  

  //print kernel
  //Plot_kernel_K(2000);
  

  
   



  return; 

}
