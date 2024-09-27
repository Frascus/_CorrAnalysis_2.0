#include "../include/VMD_Regge_charmed.h"



using namespace std;



class P_mod_val{
public:
  P_mod_val() : g1(0), m1(0), mp(0), xg(0), meas(0), err(0) {}
  double g1,m1,mp,xg, meas, err;
};

class P_mod_fit{

public:
  P_mod_fit() : A(0), geff(0), m2(0) {}
  P_mod_fit(const Vfloat & par) { if((signed)par.size() != 3) crash("par size != 3");
    A=par[0];
    geff=par[1];
    m2=par[2];
  }
  
  double A, geff,m2;
};



void Regge() {


  int Nj= 500;

  bool UseJack=0;
  int verbosity=1;
  int Nboots=100;
  int Njacks=15;

  //print prediction of form factors from VMD and Reggeized poles. s0 ~ 1GeV, untwisted trajectory.

  distr_t a0_Ds_V(UseJack), a1_Ds_V(UseJack), f_Ds(UseJack), M_Ds(UseJack), f_Ds_star(UseJack), M_Ds_star(UseJack), gV_lp(UseJack), f_Ds_prime(UseJack), gA_lp(UseJack), M_Ds_prime(UseJack), gA_lp_hq(UseJack);

  distr_t a0_D_V(UseJack), a1_D_V(UseJack), f_D(UseJack), M_D(UseJack), f_D_star(UseJack), M_D_star(UseJack), gV_D_lp(UseJack), f_D_prime(UseJack), gA_D_lp(UseJack), M_D_prime(UseJack), gA_D_lp_hq(UseJack);

  int Ngammas=100;

  Vfloat xg;
  for(int g=0;g<Ngammas;g++) xg.push_back( g*1.0/(double)Ngammas);
  
  distr_t_list Regge_propagator(UseJack,Ngammas), FV_Ds_VMD(UseJack,Ngammas), FV_Ds_R_VMD(UseJack,Ngammas), FA_Ds_VMD(UseJack, Ngammas), FA_Ds_VMD_hq(UseJack,Ngammas);
  distr_t_list Regge_propagator_D(UseJack,Ngammas), FV_D_VMD(UseJack,Ngammas), FV_D_R_VMD(UseJack,Ngammas), FA_D_VMD(UseJack,Ngammas), FA_D_VMD_hq(UseJack,Ngammas);

  GaussianMersenne Gauss_a0_Ds_V(45433,0.1);
  GaussianMersenne Gauss_a1_Ds_V(652157423,0.02);

  //GaussianMersenne Gauss_a0_Ds_V(45433,0.003/sqrt(Njacks-1));
  //GaussianMersenne Gauss_a1_Ds_V(652157423,0.05/sqrt(Njacks-1));
  
  GaussianMersenne Gauss_f_Ds(7654433,0.012);
  GaussianMersenne Gauss_M_Ds(76533,0.0001);
  GaussianMersenne Gauss_M_Ds_star(765212133,0.0001);
  GaussianMersenne Gauss_f_Ds_star(79914131,0.016);
  GaussianMersenne Gauss_gV_lp(15175931, 0.02);
  GaussianMersenne Gauss_M_Ds_prime(1165212133,0.0001);
  GaussianMersenne Gauss_f_Ds_prime(65414131,0.030);
  GaussianMersenne Gauss_gA_lp(7771931, 0.025); //assuming ~20-30% errors on qm estimate
  GaussianMersenne Gauss_gA_lp_hq(5341121,0.04); //assuming ~ 20-30% errors on qm estimate
  
  //GaussianMersenne Gauss_a0_D_V(45433,0.08);
  //GaussianMersenne Gauss_a1_D_V(652157423,0.03);
  GaussianMersenne Gauss_a0_D_V(45433,0.08);
  GaussianMersenne Gauss_a1_D_V(652157423,0.02);
  GaussianMersenne Gauss_f_D(7654433,0.014);
  GaussianMersenne Gauss_M_D(76533,0.0001);
  GaussianMersenne Gauss_M_D_star(765212133,0.0001);
  GaussianMersenne Gauss_f_D_star(79914131,0.020);
  GaussianMersenne Gauss_gV_D_lp(15175931, 0.06);
  GaussianMersenne Gauss_M_D_prime(11652123,0.0001);
  GaussianMersenne Gauss_f_D_prime(65414131,0.030);
  GaussianMersenne Gauss_gA_D_lp(7771931, 0.02);  //assuming ~20-30% errors on qm estimate
  GaussianMersenne Gauss_gA_D_lp_hq(5341121,0.06); //assuming ~20-30% errors on qm estimate


  
  auto R = [&](double t, double a0, double a1, double m) -> double { return 1.0*a1*pow(m*m, a0+a1*t-1)*tgamma(1-a0-a1*t);};  //Regge propagator




  //read form factors from file
  //##################################################################################
  file_t FV_D_cont_data, FA_D_cont_data, FV_Ds_cont_data, FA_Ds_cont_data;


  FV_D_cont_data.Read("../data/form_factors/VMD/cont_FV_D.dat");
  FA_D_cont_data.Read("../data/form_factors/VMD/cont_FA_D.dat");
  FV_Ds_cont_data.Read("../data/form_factors/VMD/cont_FV_Ds.dat");
  FA_Ds_cont_data.Read("../data/form_factors/VMD/cont_FA_Ds.dat");

  Vfloat FV_D_cont, FA_D_cont, FV_Ds_cont, FA_Ds_cont;
  Vfloat xg_V_D, xg_A_D, xg_V_Ds, xg_A_Ds;
  Vfloat FV_D_cont_err, FA_D_cont_err, FV_Ds_cont_err, FA_Ds_cont_err;

  FV_D_cont = FV_D_cont_data.col(1);
  xg_V_D= FV_D_cont_data.col(0);
  FV_D_cont_err = FV_D_cont_data.col(2);
  
  FA_D_cont = FA_D_cont_data.col(1);
  xg_A_D= FA_D_cont_data.col(0);
  FA_D_cont_err = FA_D_cont_data.col(2);

  FV_Ds_cont = FV_Ds_cont_data.col(1);
  xg_V_Ds= FV_Ds_cont_data.col(0);
  FV_Ds_cont_err = FV_Ds_cont_data.col(2);
  
  FA_Ds_cont = FA_Ds_cont_data.col(1);
  xg_A_Ds= FA_Ds_cont_data.col(0);
  FA_Ds_cont_err = FA_Ds_cont_data.col(2);
  //#################################################################################
  

  
  //BOOTSTRAP POLE MODEL FIT
  //#################################################################################
  bootstrap_fit<P_mod_fit, P_mod_val> fit_Ds_V(Nboots);
  bootstrap_fit<P_mod_fit, P_mod_val> fit_Ds_A(Nboots);
  bootstrap_fit<P_mod_fit, P_mod_val> fit_D_V(Nboots);
  bootstrap_fit<P_mod_fit, P_mod_val> fit_D_A(Nboots);
  auto MEAS = [](const P_mod_fit& par, const P_mod_val& val) -> double { return val.meas;};
  auto ERR = [](const P_mod_fit& par, const P_mod_val& val) -> double {return val.err;};
  auto ANSATZ =  [](const P_mod_fit& par, const P_mod_val& val) -> double {

    double k= val.xg*val.mp/2.0;

    double ret_val =  par.A+ val.g1/(2.0*sqrt( pow(val.m1,2) + k*k)*(sqrt(pow(val.m1,2) +k*k) + k - val.mp));

    ret_val += par.geff/(2.0*sqrt( pow(par.m2,2) + k*k)*(sqrt( pow(par.m2,2) +k*k) +k- val.mp));

    return ret_val;
    
  };

  fit_Ds_V.measurement= MEAS;
  fit_Ds_V.error= ERR;
  fit_Ds_V.ansatz= ANSATZ;
  fit_Ds_A.measurement= MEAS;
  fit_Ds_A.error= ERR;
  fit_Ds_A.ansatz= ANSATZ;
  fit_D_V.measurement= MEAS;
  fit_D_V.error= ERR;
  fit_D_V.ansatz= ANSATZ;
  fit_D_A.measurement= MEAS;
  fit_D_A.error= ERR;
  fit_D_A.ansatz= ANSATZ;


  //couplings
  double gV_Ds = -1.968*2.112*0.272;
  double gA_Ds = 1.968*2.460*0.231*0.075;
  double gV_D = -1.869*2.010*0.245*0.47;
  double gA_D = 1.869*2.421*0.211*0.20;
  
  

  fit_Ds_V.Set_number_of_measurements(7);
  fit_Ds_V.Set_verbosity(verbosity);
  fit_Ds_V.Add_par("A",0.1, 0.01);
  fit_Ds_V.Add_par("geff", -gV_Ds, 0.01);
  fit_Ds_V.Add_par("meff", 5.0, 0.1);
  //fit_Ds_V.Fix_par("geff",-gV_Ds);
  //fit_Ds_V.Fix_par("meff",2.57);
  fit_Ds_V.Fix_par("A",0.0);
  
  fit_Ds_A.Set_number_of_measurements(7);
  fit_Ds_A.Set_verbosity(verbosity);
  fit_Ds_A.Add_par("A",0.05, 0.002);
  fit_Ds_A.Add_par("geff", -gA_Ds, 0.01);
  fit_Ds_A.Add_par("meff", 2.0, 0.1);
  //fit_Ds_A.Fix_par("A",0.0);
  fit_Ds_A.Fix_par("geff",0.0);
  fit_Ds_A.Fix_par("meff",2.6);

  fit_D_V.Set_number_of_measurements(5);
  fit_D_V.Set_verbosity(verbosity);
  fit_D_V.Add_par("A",0.2, 0.01);
  fit_D_V.Add_par("geff",-gV_D, 0.1);
  fit_D_V.Add_par("meff", 4, 0.01);
  fit_D_V.Fix_par("A",0.00);
  //fit_D_V.Fix_par("meff", 2.4 );
  //fit_D_V.Fix_par("geff", -gV_D);
  
  fit_D_A.Set_number_of_measurements(7);
  fit_D_A.Set_verbosity(verbosity);
  fit_D_A.Add_par("A",0.2, 0.01);
  fit_D_A.Add_par("geff",-gA_D, 0.01);
  fit_D_A.Add_par("meff", 2.7, 0.1);
  fit_D_A.Fix_par("geff",0.0);
  fit_D_A.Fix_par("meff",1.0);
  //fit_D_A.Fix_par("A",0.0);
	      
  

  //prepare gaussian number generator
  GaussianMersenne GM(432);
 
 
  //###########################################################################################
  //###########################################################################################
  //generate bootstrap sample for Ds_V
  vector<P_mod_val>  boot_par_Ds_V;
 
  for(int iboot=0;iboot<Nboots;iboot++) {
    
    fit_Ds_V.ib =&iboot;
    double mp= 1.968+GM()*0.0005;
    double m1= 2.112+GM()*0.0005;
    double f= 0.272+GM()*0.016;
    double g=-0.10+GM()*0.02;
    for(int imeas=0;imeas<fit_Ds_V.Get_number_of_measurements();imeas++) {
      double meas= FV_Ds_cont[imeas];
      double err = FV_Ds_cont_err[imeas];
      P_mod_val X;
      X.xg= xg_V_Ds[imeas];
      X.mp= mp;
      X.m1= m1;
      X.g1 = mp*m1*f*g;	
      X.meas=meas+ GM()*err;
      X.err=err;
      fit_Ds_V.Append_to_input_par(X);
      boot_par_Ds_V.push_back(X);
    }
  }
  //Perform bootstrap fit for Ds_V
  boot_fit_data<P_mod_fit> Fit_output_Ds_V= fit_Ds_V.Perform_bootstrap_fit();
  Vfloat Fit_func_Ds_V, Fit_func_Ds_V_err;
  VVfloat pars_Ds_V(4);
  for(int iboot=10;iboot<Nboots;iboot++) { pars_Ds_V[3].push_back(Fit_output_Ds_V.chi2[iboot]);  pars_Ds_V[0].push_back(Fit_output_Ds_V.par[iboot].A); pars_Ds_V[1].push_back(Fit_output_Ds_V.par[iboot].geff); pars_Ds_V[2].push_back(Fit_output_Ds_V.par[iboot].m2);}
  for(auto &x: xg) {
    Vfloat boot_data;
    for(int iboot=0;iboot<Nboots;iboot++) {
      P_mod_val v_b;
      v_b.xg = x;
      v_b.g1 = boot_par_Ds_V[iboot].g1;
      v_b.mp = boot_par_Ds_V[iboot].mp;
      v_b.m1 = boot_par_Ds_V[iboot].m1;
      boot_data.push_back(ANSATZ(Fit_output_Ds_V.par[iboot], v_b));
    }
    Fit_func_Ds_V.push_back(Boot_ave(boot_data));
    Fit_func_Ds_V_err.push_back(Boot_err(boot_data));
  }
 
  //###########################################################################################
  //###########################################################################################








  
  

  //###########################################################################################
  //###########################################################################################
  //generate bootstrap sample for Ds_A
  vector<P_mod_val>  boot_par_Ds_A;
  for(int iboot=0;iboot<Nboots;iboot++) {
    fit_Ds_A.ib =&iboot;
    double mp= 1.968+GM()*0.0005;
    double m1= 2.460+GM()*0.0005;
    double f= 0.231+GM()*0.030;
    double g=0.075+GM()*0.025;
    for(int imeas=0;imeas<fit_Ds_A.Get_number_of_measurements();imeas++) {
      double meas= FA_Ds_cont[imeas];
      double err = FA_Ds_cont_err[imeas];
      P_mod_val X;
      X.xg= xg_A_Ds[imeas];
      X.mp= mp;
      X.m1= m1;
      X.g1 = mp*m1*f*g;	
      X.meas=meas+ GM()*err;
      X.err=err;
      fit_Ds_A.Append_to_input_par(X);
      boot_par_Ds_A.push_back(X);
    }
  }
  //Perform bootstrap fit for Ds_A
  boot_fit_data<P_mod_fit> Fit_output_Ds_A= fit_Ds_A.Perform_bootstrap_fit();
  Vfloat Fit_func_Ds_A, Fit_func_Ds_A_err;
  VVfloat pars_Ds_A(4);
  for(int iboot=10;iboot<Nboots;iboot++) { pars_Ds_A[3].push_back(Fit_output_Ds_A.chi2[iboot]); pars_Ds_A[0].push_back(Fit_output_Ds_A.par[iboot].A); pars_Ds_A[1].push_back(Fit_output_Ds_A.par[iboot].geff); pars_Ds_A[2].push_back(Fit_output_Ds_A.par[iboot].m2);}
  for(auto &x: xg) {
    Vfloat boot_data;
    for(int iboot=0;iboot<Nboots;iboot++) {
      P_mod_val v_b;
      v_b.xg = x;
      v_b.g1 = boot_par_Ds_A[iboot].g1;
      v_b.mp = boot_par_Ds_A[iboot].mp;
      v_b.m1 = boot_par_Ds_A[iboot].m1;
      boot_data.push_back(ANSATZ(Fit_output_Ds_A.par[iboot], v_b));
    }
    Fit_func_Ds_A.push_back(Boot_ave(boot_data));
    Fit_func_Ds_A_err.push_back(Boot_err(boot_data));
  }
  //###########################################################################################
  //###########################################################################################

  

  
  
  //###########################################################################################
  //###########################################################################################
  //generate bootstrap sample for D_V
  vector<P_mod_val> boot_par_D_V;
  for(int iboot=0;iboot<Nboots;iboot++) {
    fit_D_V.ib =&iboot;
    double mp= 1.869+GM()*0.0005;
    double m1= 2.010+GM()*0.0005;
    double f= 0.245+GM()*0.020;
    double g=-0.47+GM()*0.06;
    for(int imeas=0;imeas<fit_D_V.Get_number_of_measurements();imeas++) {
      double meas= FV_D_cont[imeas];
      double err = FV_D_cont_err[imeas];
      P_mod_val X;
      X.xg= xg_V_D[imeas];
      X.mp= mp;
      X.m1= m1;
      X.g1 = mp*m1*f*g;	
      X.meas=meas+ GM()*err;
      X.err=err;
      fit_D_V.Append_to_input_par(X);
      boot_par_D_V.push_back(X);
    }
  }

  
  //Perform bootstrap fit for D_V
  boot_fit_data<P_mod_fit> Fit_output_D_V= fit_D_V.Perform_bootstrap_fit();
  Vfloat Fit_func_D_V, Fit_func_D_V_err;
  VVfloat pars_D_V(4);
  for(int iboot=10;iboot<Nboots;iboot++) {  pars_D_V[3].push_back(Fit_output_D_V.chi2[iboot]); pars_D_V[0].push_back(Fit_output_D_V.par[iboot].A); pars_D_V[1].push_back(Fit_output_D_V.par[iboot].geff); pars_D_V[2].push_back(Fit_output_D_V.par[iboot].m2);}
  for(auto &x: xg) {
    Vfloat boot_data;
    for(int iboot=0;iboot<Nboots;iboot++) {
      P_mod_val v_b;
      v_b.xg = x;
      v_b.g1 = boot_par_D_V[iboot].g1;
      v_b.mp = boot_par_D_V[iboot].mp;
      v_b.m1 = boot_par_D_V[iboot].m1;
      boot_data.push_back(ANSATZ(Fit_output_D_V.par[iboot], v_b));
    }
    Fit_func_D_V.push_back(Boot_ave(boot_data));
    Fit_func_D_V_err.push_back(Boot_err(boot_data));
  }
  
  //###########################################################################################
  //###########################################################################################





  
  
  
  //###########################################################################################
  //###########################################################################################
 
  //generate bootstrap sample for D_A
  vector<P_mod_val> boot_par_D_A;

  for(int iboot=0;iboot<Nboots;iboot++) {
    fit_D_A.ib =&iboot;
    double mp= 1.869+GM()*0.0005;
    double m1= 2.421+GM()*0.0005;
    double f= 0.211+GM()*0.030;
    double g= 0.20+GM()*0.02;
    for(int imeas=0;imeas<fit_D_A.Get_number_of_measurements();imeas++) {
      double meas= FA_D_cont[imeas];
      double err = FA_D_cont_err[imeas];
      P_mod_val X;
      X.xg= xg_A_D[imeas];
      X.mp= mp;
      X.m1= m1;
      X.g1 = mp*m1*f*g;	
      X.meas=meas+ GM()*err;
      X.err=err;
      fit_D_A.Append_to_input_par(X);
      boot_par_D_A.push_back(X);
    }
  }
  //Perform bootstrap fit for D_A
  boot_fit_data<P_mod_fit> Fit_output_D_A= fit_D_A.Perform_bootstrap_fit();
  Vfloat Fit_func_D_A, Fit_func_D_A_err;
  VVfloat pars_D_A(4);
  for(int iboot=10;iboot<Nboots;iboot++) {  pars_D_A[3].push_back(Fit_output_D_A.chi2[iboot]); pars_D_A[0].push_back(Fit_output_D_A.par[iboot].A); pars_D_A[1].push_back(Fit_output_D_A.par[iboot].geff); pars_D_A[2].push_back(Fit_output_D_A.par[iboot].m2);}
  for(auto &x: xg) {
    Vfloat boot_data;
    for(int iboot=0;iboot<Nboots;iboot++) {
      P_mod_val v_b;
      v_b.xg = x;
      v_b.g1 = boot_par_D_A[iboot].g1;
      v_b.mp = boot_par_D_A[iboot].mp;
      v_b.m1 = boot_par_D_A[iboot].m1;
      boot_data.push_back(ANSATZ(Fit_output_D_A.par[iboot], v_b));
    }
    Fit_func_D_A.push_back(Boot_ave(boot_data));
    Fit_func_D_A_err.push_back(Boot_err(boot_data));
  }

  //###########################################################################################
  //###########################################################################################
 

  //Print meff, geff and A for the 4 form factors
  cout<<"#Printing fit parameters#"<<endl;
  cout<<"########################"<<endl;
  cout<<"A(F_V_Ds): "<<Boot_ave(pars_Ds_V[0])<<"   "<<Boot_err(pars_Ds_V[0])<<endl;
  cout<<"g(F_V_Ds): "<<gV_Ds<<endl;
  cout<<"geff(F_V_Ds): "<<Boot_ave(pars_Ds_V[1])<<"   "<<Boot_err(pars_Ds_V[1])<<endl;
  cout<<"meff(F_V_Ds): "<<Boot_ave(pars_Ds_V[2])<<"   "<<Boot_err(pars_Ds_V[2])<<endl;
  cout<<"average ch2: "<<Boot_ave(pars_Ds_V[3])<<"    "<<Boot_err(pars_Ds_V[3])<<endl;
  cout<<"########################"<<endl;
  cout<<"A(F_A_Ds): "<<Boot_ave(pars_Ds_A[0])<<"   "<<Boot_err(pars_Ds_A[0])<<endl;
  cout<<"g(F_A_Ds): "<<gA_Ds<<endl;
  cout<<"geff(F_A_Ds): "<<Boot_ave(pars_Ds_A[1])<<"   "<<Boot_err(pars_Ds_A[1])<<endl;
  cout<<"meff(F_A_Ds): "<<Boot_ave(pars_Ds_A[2])<<"   "<<Boot_err(pars_Ds_A[2])<<endl;
  cout<<"average ch2: "<<Boot_ave(pars_Ds_A[3])<<"    "<<Boot_err(pars_Ds_A[3])<<endl;
  cout<<"########################"<<endl;
  cout<<"A(F_V_D): "<<Boot_ave(pars_D_V[0])<<"   "<<Boot_err(pars_D_V[0])<<endl;
  cout<<"g(F_V_D): "<<gV_D<<endl;
  cout<<"geff(F_V_D): "<<Boot_ave(pars_D_V[1])<<"   "<<Boot_err(pars_D_V[1])<<endl;
  cout<<"meff(F_V_D): "<<Boot_ave(pars_D_V[2])<<"   "<<Boot_err(pars_D_V[2])<<endl;
  cout<<"average ch2: "<<Boot_ave(pars_D_V[3])<<"    "<<Boot_err(pars_D_V[3])<<endl;
  cout<<"########################"<<endl;
  cout<<"A(F_A_D): "<<Boot_ave(pars_D_A[0])<<"   "<<Boot_err(pars_D_A[0])<<endl;
  cout<<"g(F_A_D): "<<gA_D<<endl;
  cout<<"geff(F_A_D): "<<Boot_ave(pars_D_A[1])<<"   "<<Boot_err(pars_D_A[1])<<endl;
  cout<<"meff(F_A_D): "<<Boot_ave(pars_D_A[2])<<"   "<<Boot_err(pars_D_A[2])<<endl;
  cout<<"average ch2: "<<Boot_ave(pars_D_A[3])<<"    "<<Boot_err(pars_D_A[3])<<endl;
  cout<<"########################"<<endl;









				

  

  //generate jackknife data for VMD and Reggeized_VMD_fit
  
  for(int nj=0;nj<Nj;nj++) {

    //double a0 = -1.09 + Gauss_a0_Ds_V();
    //double a0=-1.42 + Gauss_a0_Ds_V();
    double a0 = -1.32 + Gauss_a0_Ds_V();
    //double a1 = 0.47 + Gauss_a1_Ds_V();
  
    double fds = 0.231 +Gauss_f_Ds();
    double Mds= 1.968 + Gauss_M_Ds();
    double fds_star= 0.272 + Gauss_f_Ds_star();
    double Mds_star =2.112 + Gauss_M_Ds_star();
    double a1= (1.0- a0)/(Mds_star*Mds_star);
    double gVlp= -0.10 + Gauss_gV_lp();
    double gAlp = 0.075 + Gauss_gA_lp();
    double fds_prime = 0.231 +Gauss_f_Ds_prime();
    double Mds_prime = 2.460 + Gauss_M_Ds_prime();
    double gAlp_hq = 0.23 + Gauss_gA_lp_hq();
    a0_Ds_V.distr.push_back(a0);
    a1_Ds_V.distr.push_back(a1);
    f_Ds.distr.push_back( fds);
    M_Ds.distr.push_back(Mds);
    f_Ds_star.distr.push_back(fds_star);
    M_Ds_star.distr.push_back(Mds_star);
    gV_lp.distr.push_back(gVlp);
    f_Ds_prime.distr.push_back(fds_prime);
    M_Ds_prime.distr.push_back(Mds_prime);
    gA_lp.distr.push_back(gAlp);
    gA_lp_hq.distr.push_back(gAlp_hq);

    //double a0_D = -1.05 + Gauss_a0_D_V();
    double a0_D = -1.18 + Gauss_a0_D_V();
   
    // double a0_D = -1.27  + Gauss_a0_D_V();
    //double a1_D = 0.50 + Gauss_a1_D_V();
    //double a1_D = 0.55 + Gauss_a1_D_V();
    double fd = 0.211 +Gauss_f_D();
    double Md= 1.869 + Gauss_M_D();
    
    double fd_star= 0.245 + Gauss_f_D_star();
    double Md_star = 2.010 + Gauss_M_D_star();
    double a1_D = (1.0-a0_D)/(Md_star*Md_star);
    //double gVlpD= -0.47 + Gauss_gV_D_lp();
    double gVlpD= -0.10+ Gauss_gV_lp();
    double gAlpD = 0.20 + Gauss_gA_D_lp();
    double fd_prime = 0.211 +Gauss_f_D_prime();
    double Md_prime = 2.421 + Gauss_M_D_prime();
    double gAlpD_hq = 0.35 + Gauss_gA_D_lp_hq();
    a0_D_V.distr.push_back(a0_D);
    a1_D_V.distr.push_back(a1_D);
    f_D.distr.push_back( fd);
    M_D.distr.push_back(Md);
    f_D_star.distr.push_back(fd_star);
    M_D_star.distr.push_back(Md_star);
    gV_D_lp.distr.push_back(gVlpD);
    f_D_prime.distr.push_back(fd_prime);
    M_D_prime.distr.push_back(Md_prime);
    gA_D_lp.distr.push_back(gAlpD);
    gA_D_lp_hq.distr.push_back(gAlpD_hq);
    
    int xg_count=0;
    for(auto &x: xg) {

      double k = x*Mds/2.0;

    
      //Ds
      //cout<<"regge Ds"<<endl<<flush;
      Regge_propagator.distr_list[xg_count].distr.push_back( R( pow(Mds,2)*(1-x), a0, a1, Mds));
      FV_Ds_VMD.distr_list[xg_count].distr.push_back( fds_star*Mds*Mds_star*gVlp/(2.0*sqrt( pow(Mds_star,2)+ k*k)*(sqrt( pow(Mds_star,2)+ k*k) + k - Mds))) ;
      
      FV_Ds_R_VMD.distr_list[xg_count].distr.push_back( fds_star*Mds*Mds_star*gVlp*R(pow(Mds,2)*(1-x), a0,a1,Mds));
      //cout<<"#####  Ds"<<endl<<flush;
      FA_Ds_VMD.distr_list[xg_count].distr.push_back(fds_prime*Mds*Mds_prime*gAlp/(2.0*sqrt( pow(Mds_prime,2)+ k*k)*(sqrt( pow(Mds_prime,2)+ k*k) + k - Mds)));
      FA_Ds_VMD_hq.distr_list[xg_count].distr.push_back(fds_prime*Mds*Mds_prime*gAlp_hq/(2.0*sqrt( pow(Mds_prime,2)+ k*k)*(sqrt( pow(Mds_prime,2)+ k*k) + k - Mds)));

      k= x*Md/2.0;

     
      //D
      //cout<<"Regge D"<<endl<<flush;
      Regge_propagator_D.distr_list[xg_count].distr.push_back( R( pow(Md,2)*(1-x), a0_D, a1_D, Md));
      FV_D_VMD.distr_list[xg_count].distr.push_back( fd_star*Md*Md_star*gVlpD/(2.0*sqrt( pow(Md_star,2)+ k*k)*(sqrt( pow(Md_star,2)+ k*k) + k - Md))) ;
      //cout<<"#####   D xg: "<<x<<endl<<flush;
      FV_D_R_VMD.distr_list[xg_count].distr.push_back( fd_star*Md*Md_star*gVlpD*R(pow(Md,2)*(1-x), a0_D,a1_D,Md));
     
      FA_D_VMD.distr_list[xg_count].distr.push_back(fd_prime*Md*Md_prime*gAlpD/(2.0*sqrt( pow(Md_prime,2)+ k*k)*(sqrt( pow(Md_prime,2)+ k*k) + k - Md)));
      FA_D_VMD_hq.distr_list[xg_count].distr.push_back(fd_prime*Md*Md_prime*gAlpD_hq/(2.0*sqrt( pow(Md_prime,2)+ k*k)*(sqrt( pow(Md_prime,2)+ k*k) + k - Md)));
     
           

      xg_count++;
	

    }
  }

  //print predictions to FILE
  
  boost::filesystem::create_directory("../data");
  boost::filesystem::create_directory("../data/form_factors");
  boost::filesystem::create_directory("../data/form_factors/VMD");

  //Ds
  Print_To_File({}, {xg, Regge_propagator.ave(), Regge_propagator.err(), FV_Ds_VMD.ave(), FV_Ds_VMD.err(), FV_Ds_R_VMD.ave(), FV_Ds_R_VMD.err(), Fit_func_Ds_V, Fit_func_Ds_V_err}, "../data/form_factors/VMD/FV_Ds.dat", "", "# int   xg    regge_prop    VMD    VMD_err     Reggeized_VMD       Reggeized_VMD_err     2-pole          2-pole_err" );
  Print_To_File({}, {xg, Regge_propagator.ave(), Regge_propagator.err(), FA_Ds_VMD.ave(), FA_Ds_VMD.err(),FA_Ds_VMD_hq.ave(), FA_Ds_VMD_hq.err(), Fit_func_Ds_A, Fit_func_Ds_A_err }, "../data/form_factors/VMD/FA_Ds.dat", "", "# int   xg    regge_prop    VMD    VMD_err    VMD_hq      VMD_hq_err      2-pole      2-pole_err" );


  //D
  Print_To_File({}, {xg, Regge_propagator_D.ave(), Regge_propagator_D.err(), FV_D_VMD.ave(), FV_D_VMD.err(), FV_D_R_VMD.ave(), FV_D_R_VMD.err(), Fit_func_D_V, Fit_func_D_V_err}, "../data/form_factors/VMD/FV_D.dat", "", "# int   xg    regge_prop    VMD    VMD_err     Reggeized_VMD       Reggeized_VMD_err     2-pole      2-pole_err" );
  Print_To_File({}, {xg, Regge_propagator_D.ave(), Regge_propagator_D.err(), FA_D_VMD.ave(), FA_D_VMD.err(),FA_D_VMD_hq.ave(), FA_D_VMD_hq.err(), Fit_func_D_A, Fit_func_D_A_err}, "../data/form_factors/VMD/FA_D.dat", "", "# int   xg    regge_prop    VMD    VMD_err    VMD_hq      VMD_hq_err       2-pole          2-pole_err" );
  
  
 

 
  
    


  return ;

}
