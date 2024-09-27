#include "../include/3pts_mes_gamma_W.h"

using namespace std;
namespace plt = matplotlibcpp;


const double M2PiPhys=pow(0.135,2);
const double alpha_em = 1/137.04;
const double e2 = alpha_em*4.0*M_PI;
const int Nbranches = 8;
const int Nboots= 100;
bool UseJack=1;
const int nboots=150;
const int Njacks=30;
int sm_lev=1;
const double e_f1 = 2.0/3.0; //electric charge of u-type quark
const double e_f2 = -1.0/3.0; //electric charge of d-type quark
bool Include_k0_noise=1; //1
bool VIRTUAL_RUN=1;
bool Determine_contaminations=0;
int SUB_ZERO_MOMENTUM_VIRTUAL_INT= 0; //1
double SUB_ZERO_MOMENTUM_VIRTUAL=(double)SUB_ZERO_MOMENTUM_VIRTUAL_INT;
int verbosity=1; //used in Fit Contaminations to print infos about minimization
bool FIT_VIRTUAL_FF=0; //1
bool USE_FITTED_FF=0; //1
bool COMPUTE_l4_DECAY_RATE=0; //1
int VIRTUAL_ESTIMATOR_SET=3;
bool FIND_k0_KIN=false;
bool FIND_k0_kz0_KIN=false;
const string Meson="K_def";
const int Include_FSES=0;

class ExpFit {

 
public:
  ExpFit(const Vfloat &par)  { if((signed)par.size() != 5) crash("par size != 5");
    a1=par[0];
    a2=par[1];
    a3=par[2];
    a4=par[3];
    a5=par[4];

  }

  double a1;
  double a2;
  double a3;
  double a4;
  double a5;
};

class ExpFit_der {

 
public:
  ExpFit_der(const Vfloat &par)  { if((signed)par.size() != 4) crash("In ExpFit_der par size != 5");
    a2=par[1];
    a3=par[2];
    a4=par[3];
    a5=par[4];
  }

  double a2;
  double a3;
  double a4;
  double a5;
};



class Val{
public:
  Val() : t(0), meas(0), err(0) {}
  double t;
  double meas;
  double err;
};



void Plot_form_factors(string W, distr_t_list& F, distr_t& F_fit, int Tmin, int Tmax, int Nt, string Ens_tag, double xg, double offsh, int smearing) {

  boost::filesystem::create_directory("../plots");
  boost::filesystem::create_directory("../plots/form_factors");
  boost::filesystem::create_directory("../plots/form_factors/"+Meson);
  plt::clf();
  plt::figure_size(1500,1000);
  plt::xlabel("$t/a$");
  if(W=="V") plt::ylabel("$F_{V}$");
  else plt::ylabel("$F_{A}$");
  plt::xlim(3, Nt/2 -3);
  double y_min=1e10, y_max= -1e10;
  for(int t=3;t<=Nt/2 -3;t++) { y_min = min( y_min, F.ave()[t] - 1.2*F.err()[t]); y_max= max( y_max, F.ave()[t] + 1.2*F.err()[t]);}  
  plt::ylim(y_min, y_max);
  plt::grid(4);
  Vfloat TT;
  Vfloat TT_tot;
  for(int t=0;t<Nt;t++) TT_tot.push_back(t);
  //generate data_point for fit;
  Vfloat F_val(Tmax - Tmin +1, F_fit.ave());
  Vfloat F_err(Tmax - Tmin +1, F_fit.err());
  for(int t=Tmin;t<=Tmax;t++)  TT.push_back(t);
  plt::errorbar(TT_tot, F.ave(), F.err(), { {"c", "black"}, {"marker", "."} , {"ls" , "-"}, {"label", "xg="+to_string_with_precision(xg,4)}});
  plt::errorbar(TT, F_val, F_err, { {"c", "red"}, {"marker", "."}, {"ls", "-"}, {"label", "fitted value"}});
  // Enable legend.
  plt::legend();
  string figure_path;
  if(W=="V") figure_path= "../plots/form_factors/"+Meson+"/FV_"+Ens_tag+"_xg_"+to_string_with_precision(xg,4)+"_offsh_"+to_string_with_precision(offsh,4)+"_smlev_"+to_string(smearing)+"_k0_noise_"+to_string(Include_k0_noise)+".png";
  else figure_path= "../plots/form_factors/"+Meson+"/FA_"+Ens_tag+"_xg_"+to_string_with_precision(xg,4)+"_offsh_"+to_string_with_precision(offsh,4)+"_smlev_"+to_string(smearing)+"_k0_noise_"+to_string(Include_k0_noise)+".png";
  plt::save(figure_path.c_str());
  plt::close();
  return;
}

void Fit_contaminations(string W, distr_t_list& F, distr_t &F_fit, string Ens_tag, int Nt, double xg, double offsh, int smearing) {

  

  double T_min, T_max_sx, T_min_dx, T_max;
  T_min=2;
  int t_shift=1;
  if(W=="RV" || W=="FV") T_max= Nt/2 - 3;
  else T_max=Nt/2 -2;
  T_max_sx= (double)min(14, Nt/4);
  T_min_dx= T_max_sx+1;
  int ndof= (T_max -T_min+1 - 5);
  //choose the t-interval
  auto ansatz_sx= [&](const Vfloat& par, double t, int imeas) -> double { return par[0] +par[1]*exp(-1*par[2]*t);};
  auto ansatz_dx= [&](const Vfloat& par, double t, int imeas) -> double { return par[0] +par[1]*exp(-1*par[2]*(Nt/2-t));};
  auto ansatz_tot=[&](const Vfloat& par, double t, int imeas) -> double { return par[0] +par[1]*exp(-1*par[2]*t) + par[3]*exp(-1*par[4]*(Nt/2-t));};
  auto ansatz = [=](const ExpFit& par, const Val& val) -> double {return par.a1 + par.a2*exp(-1*par.a3*val.t) +  par.a4*exp(-1*par.a5*(Nt/2-val.t));};
  Vfloat T_dx, T_sx, Y_dx, Y_sx, Y_err_dx, Y_err_sx, Y_dx_sx;
  for(int t=T_min;t<=T_max_sx;t++) {T_sx.push_back(t); Y_sx.push_back( F.ave()[t]); Y_err_sx.push_back( F.err()[t]);}
  for(int t=T_min_dx;t<=T_max;t++) {T_dx.push_back(t); Y_dx.push_back( F.ave()[t]); Y_err_dx.push_back( F.err()[t]);}


  double a1_guess, a1_guess_dx, a2_guess, a3_guess, a4_guess,a5_guess;
  pair<double,double> fit_interval_a3;
  pair<double, double> fit_interval_a5;
  int sign_a2, sign_a4;


  //choose guess intervals depending on the observable
  //###############################################################
  fit_interval_a3= make_pair(0.01, 1.0);
  fit_interval_a5= make_pair(0.01, 1.0);
  if( (W=="RV" || W=="FV")) {
    if(smearing) {
      if(W=="FV") {
	sign_a2=-1;
	sign_a4=1;
      }
      else {
	sign_a2=-1;
	sign_a4=1;
      }
    }
    else {
      if(W=="FV") {
	sign_a2=1;
	sign_a4=1;
      }
      else {
	sign_a2=1;
	sign_a4=-1;
      }

    }
  }
  else {
    if(W=="FA") {
      sign_a2=1;
      sign_a4=-1;
    }
    else {
      if(smearing) {
	sign_a2=-1;
	sign_a4=1;
      }
      else {
	sign_a2=-1;
	sign_a4=1;
      }
    }
  }
  //################################################################


  bool Check_min_validity=false;


  while(!Check_min_validity && t_shift < Nt/18) {

    Check_min_validity=true;
    //SX_FIT
    UniformMersenne  A2(6545534,log(0.001), log(10.0));
    UniformMersenne A4(87756756,log(0.001), log(10.0));
    UniformMersenne A3(659878,fit_interval_a3.first, fit_interval_a3.second);
    UniformMersenne A5(546353,fit_interval_a5.first, fit_interval_a5.second);
    int its_sx=0;
    bool sx_min_found=false;
    Vfloat Ch2_sx_list;
    bool sign_sx_flipped=false;
    while(!sx_min_found) {
      its_sx++;
      T_fit fit_sx(T_sx, Y_sx, Y_err_sx);
      fit_sx.ansatz = ansatz_sx;
      int ndof = (T_max_sx-T_min +1 - 3);
      double sx_1 = sign_a2*exp( A2());
      double sx_2 = A3();
      fit_sx.add_pars({a1_guess, sx_1, sx_2});
      fit_sx.add_par_errs({F_fit.err(), sx_1/10, sx_2/10});
      fit_t_res sx_fit_res = fit_sx.fit();
      Ch2_sx_list.push_back(sx_fit_res.chi2/ndof);
      if(sx_fit_res.chi2/ndof < 1) { cout<<"ch2_sx for OBS: "<<W<<" xg: "<<xg<<" ch2: "<<sx_fit_res.chi2/ndof<<"  ITs: "<<its_sx<<endl; a1_guess= sx_fit_res.pars[0];a2_guess= sx_fit_res.pars[1]; a3_guess=sx_fit_res.pars[2]; sx_min_found=true;}
      if (its_sx % 100 == 0) cout<<"IT(sx): "<<its_sx<<" OBS: "<<W<<" ,chi2/ndof: "<<sx_fit_res.chi2/ndof<<endl;
      if( (its_sx > 2000 && !sign_sx_flipped) || (its_sx > 7000 && sign_sx_flipped)) {
	//cout<<"ch2_sx for OBS: "<<W<<" xg: "<<xg<<" ch2: "<<sx_fit_res.chi2/ndof<<" ITs: "<<its_sx<<endl;
	//cout<<"do you want to continue?"<<endl;
	//cin >> sx_min_found;
	double actual_ch2 = Ch2_sx_list[Ch2_sx_list.size()-1];
	double threshold = Ch2_sx_list.size()*0.50;
	int count_ch2=0;
	if(actual_ch2 < 4) {for( auto &c : Ch2_sx_list) { if( fabs(c-actual_ch2) < 1.0e-3) count_ch2++;}}
	if((double)count_ch2 > threshold) sx_min_found=true; //ch2 is bad but seems that Minuit found the right minimum
	if(sx_min_found) { a1_guess= sx_fit_res.pars[0]; a2_guess=sx_fit_res.pars[1]; a3_guess=sx_fit_res.pars[2];}
      }
      if(its_sx > 5000 && !sign_sx_flipped) {sign_a2 *= -1; Ch2_sx_list.clear(); sign_sx_flipped=true;  } //change sign
      if(its_sx > 10000) return;  //not able to perform the fit
    }

 

    //DX_FIT
    bool dx_min_found=false;
    int its_dx=0;
    Vfloat Ch2_dx_list;
    bool sign_dx_flipped=false;
    while(!dx_min_found) {
      its_dx++;
      T_fit fit_dx(T_dx, Y_dx, Y_err_dx);
      fit_dx.ansatz= ansatz_dx;
      int ndof = (T_max-T_min_dx +1 - 3);
      double dx_1= sign_a4*exp(A4());
      double dx_2= A5();
      fit_dx.add_pars({F_fit.ave(), dx_1, dx_2});
      fit_dx.add_par_errs({F_fit.err(), dx_1/10, dx_2/10});
      fit_t_res dx_fit_res = fit_dx.fit();
      Ch2_dx_list.push_back(dx_fit_res.chi2/ndof);
      if(dx_fit_res.chi2/ndof < 1) { cout<<"ch2_dx for OBS: "<<W<<" xg: "<<xg<<" ch2: "<<dx_fit_res.chi2/ndof<<" ITs: "<<its_dx<<endl; a4_guess= dx_fit_res.pars[1]; a5_guess=dx_fit_res.pars[2]; dx_min_found=true;}
      if (its_dx % 100 == 0) cout<<"IT(dx): "<<its_dx<<" OBS: "<<W<<" ,chi2/ndof: "<<dx_fit_res.chi2/ndof<<endl;
      if((its_dx > 2000 && !sign_dx_flipped) || (its_dx > 7000 && sign_dx_flipped)) {
	//cout<<"ch2_dx for OBS: "<<W<<" xg: "<<xg<<" ch2: "<<dx_fit_res.chi2/ndof<<" ITs: "<<its_dx<<endl;
	//cout<<"do you want to continue?"<<endl;
	//cin >> dx_min_found;
	double actual_ch2 = Ch2_dx_list[Ch2_dx_list.size()-1];
	double threshold = Ch2_dx_list.size()*0.50;
	int count_ch2=0;
	if(actual_ch2 < 4) {for( auto &c : Ch2_dx_list) { if( fabs(c-actual_ch2) < 1.0e-3) count_ch2++;}}
	if((double)count_ch2 > threshold) dx_min_found=true;
	if(dx_min_found) { a1_guess_dx=dx_fit_res.pars[0] ;a4_guess=dx_fit_res.pars[1]; a5_guess=dx_fit_res.pars[2];}
    
      }
      if(its_dx > 5000 && !sign_dx_flipped) {sign_a4*= -1; Ch2_dx_list.clear(); sign_dx_flipped=true;}
      if(its_dx > 10000) return; //not able to perform the fit 
    }

  
    //BOOTSTRAP FIT
    //##############################################################################
    bootstrap_fit<ExpFit, Val> fit(Nboots);
    ndof = T_max -T_min + 1- 5 -t_shift;
    fit.measurement = [=](const ExpFit& par, const Val& val) -> double { return val.meas;};
    fit.error = [=](const ExpFit& par, const Val& val) -> double {return val.err;};
    fit.ansatz =  ansatz;
    fit.Set_number_of_measurements(T_max - T_min +1-t_shift);
    fit.Set_verbosity(verbosity);
  
    cout<<"a1_g: "<<a1_guess<<endl;
    cout<<"a2_g: "<<a2_guess<<endl;
    cout<<"a3_g: "<<a3_guess<<endl;
    cout<<"a4_g: "<<a4_guess<<endl;
    cout<<"a5_g: "<<a5_guess<<endl;
    fit.Add_par("a1", a1_guess, fabs(a1_guess)/30);
    fit.Add_par("a2", a2_guess, fabs(a2_guess)/30);
    fit.Add_par("a3", a3_guess, fabs(a3_guess)/30);
    fit.Add_par("a4", a4_guess, fabs(a4_guess)/30);
    fit.Add_par("a5", a5_guess, fabs(a5_guess)/30);
  
  
    //prepare gaussian number generator
    GaussianMersenne GM(4355);
    double F_resc=1.0;
  
    //generate bootstrap sample
    for(int iboot=0;iboot<Nboots;iboot++) {
      fit.ib =&iboot;
      for(int i=T_min+t_shift;i<=T_max;i++) {
	double meas= F.ave()[i];
	double err = F.err()[i];
	Val X;
	X.t= (double)i;
	X.meas=meas+ GM()*err/F_resc;
	X.err=err;
	fit.Append_to_input_par(X);
      }
    }
  
  
  
    //Perform bootstrap fit
    boot_fit_data<ExpFit> Fit_output= fit.Perform_bootstrap_fit();
  
  
  
    //retrieve output params
    //##############################################
  
    double ave_ch2=0;
    Vfloat a1_b, a2_b, a3_b, a4_b, a5_b;
    VVfloat Check_Min(Nt/2);
  
    int iboot_min=5;
    for(int iboot=iboot_min;iboot<Nboots;iboot++) {
      a1_b.push_back(Fit_output.par[iboot].a1);
      a2_b.push_back(Fit_output.par[iboot].a2);
      a3_b.push_back(Fit_output.par[iboot].a3);
      a4_b.push_back(Fit_output.par[iboot].a4);
      a5_b.push_back(Fit_output.par[iboot].a5);
      ave_ch2+= Fit_output.chi2[iboot]/(Nboots-iboot_min);
      for(int tt=0; tt<Nt/2;tt++) {
	Val X;
	X.t = tt;
	Check_Min[tt].push_back( (ansatz(Fit_output.par[iboot],X)/Fit_output.par[iboot].a1) -1); 
      }
    }
  
    //check reliability of the estimate
    //###############################################################à
    //determine aposteriori fit interval
    Vint No_cont_T(Nt/2,0);
    double threshold_contamination = 1;
    for(int tt=0; tt<Nt/2;tt++) {
      if (fabs(Boot_ave(Check_Min[tt])) < threshold_contamination) No_cont_T[tt]=1;
      cout<<"tt: "<<tt<<" "<<Boot_ave(Check_Min[tt])<<" "<<Boot_err(Check_Min[tt])<<endl;
    }
  
    //printV(No_cont_T, "No_cont_T", 1);  
    //find optimized interval for fit
    int T_min_apost, T_max_apost;
  
  
    bool T_min_apost_found=false;
    bool T_max_apost_found=false;

    for(int tt=0; tt<Nt/2;tt++) {
      if (No_cont_T[tt] == 1) {
	if(!T_min_apost_found) {
	  T_min_apost = tt;
	  T_min_apost_found= true;
	}
	if(!T_max_apost_found) T_max_apost = tt;
      }
      if(No_cont_T[tt]==0 && T_min_apost_found) {T_max_apost_found=true; break;}
    
    }


    double F_fit_opt_ave, F_fit_opt_err;
    if(T_min_apost_found && T_max_apost_found) {
      //Refit in the new_interval
      CorrAnalysis corr(UseJack, Njacks, nboots);
      corr.Nt=Nt;
      corr.Tmin= T_min_apost;
      corr.Tmax= T_max_apost;
      distr_t F_fit_opt = corr.Fit_distr(F);
      double error_threshold=0.2;
      //check validity of exp fit
      double a1_b_ave= Boot_ave(a1_b);
      double a1_b_err= Boot_err(a1_b);
      F_fit_opt_ave= F_fit_opt.ave();
      F_fit_opt_err= F_fit_opt.err();
    
      if( a1_b_err > error_threshold*F_fit_opt_err) {
	if( a1_b_ave < F_fit_opt_ave) {
	  if(a1_b_ave + 1.5*a1_b_err > F_fit_opt_ave - 1.5*F_fit_opt_err) Check_min_validity=true;
	}
	if (a1_b_ave >= F_fit_opt_ave) {
	  if(F_fit_opt_ave + 1.5*F_fit_opt_err > a1_b_ave - 1.5*a1_b_err) Check_min_validity=true;
	}
      }
    }
    



	//###############################################################


	if(Check_min_validity) {
	  //print output parameters
	  ofstream Print_out("../data/form_factors/"+Meson+"/Contaminations/contaminations_"+W+"_"+Ens_tag+".dat", ofstream::app);
	  Print_out<<xg<<setw(20)<<offsh<<setw(20)<<smearing<<setw(20)<<Boot_ave(a1_b)<<setw(20)<<F_resc*Boot_err(a1_b)<<setw(20)<<Boot_ave(a2_b)<<setw(20)<<F_resc*Boot_err(a2_b)<<setw(20)<<Boot_ave(a3_b)<<setw(20)<<F_resc*Boot_err(a3_b)<<setw(20)<<Boot_ave(a5_b)<<setw(20)<<F_resc*Boot_err(a5_b)<<setw(20)<<ave_ch2/ndof<<endl;
	  Print_out.close();

	  //Print form factor from costant fit in optimized interval
	  ofstream Print_opt_FF("../data/form_factors/"+Meson+"/FORM_FACTOR_LIST/FF_optimized_"+W+"_"+Ens_tag+".dat", ofstream::app);
	  Print_opt_FF<<xg<<setw(20)<<offsh<<setw(20)<<smearing<<setw(20)<<F_fit_opt_ave<<setw(20)<<F_fit_opt_err<<setw(20)<<T_min_apost<<setw(20)<<T_max_apost<<endl;
	  Print_opt_FF.close();


	  //print fitting function
	  boost::filesystem::create_directory("../data/form_factors/"+Meson+"/Contaminations/fit_func");
	  ofstream Print_func("../data/form_factors/"+Meson+"/Contaminations/fit_func/"+W+"_"+Ens_tag+"_xg_"+to_string_with_precision(xg,4)+"_offsh_"+to_string_with_precision(offsh,4)+"_smlev_"+to_string(smearing)+"_k0_noise_"+to_string(Include_k0_noise)+".dat");

	  int tsteps=500;

	  for(int is=0;is<tsteps;is++) {

	    Vfloat fitted_t;
	    double t= ((double)is)*(Nt/(2.0*(double)tsteps));
	    Val X;
	    X.t = t;
	    for(int iboot=iboot_min;iboot<Nboots;iboot++) fitted_t.push_back( ansatz(Fit_output.par[iboot], X  ));
	    Print_func<<t<<setw(20)<<Boot_ave(fitted_t)<<setw(20)<<F_resc*Boot_err(fitted_t)<<endl;
	  }


	  Print_func.close();
	}

      t_shift++;
  }
    //###############################################################################

    return;
}

void Fit_contaminations_from_derivative(string W, distr_t_list& F, distr_t &F_fit, string Ens_tag, int Nt, double xg, double offsh, int smearing) {



  //compute numerical derivative of F
  distr_t_list F_der(UseJack);
  for(int t=0;t<F.size();t++) {
    //forward derivative
    F_der.distr_list.push_back( F.distr_list[(t+1)%Nt] -F.distr_list[t]);
  }

  //Print F_der
  boost::filesystem::create_directory("../data/form_factors/"+Meson+"/corr_der");
  Print_To_File({}, {F_der.ave(), F_der.err()}, "../data/form_factors/corr_der/"+W+"_"+Ens_tag+"_xg_"+to_string_with_precision(xg,4)+"_offsh_"+to_string_with_precision(offsh,4)+"_smlev_"+to_string(smearing)+"_k0_noise_"+to_string(Include_k0_noise)+".dat", "", "");

  double T_min, T_max_sx, T_min_dx, T_max;
  T_min=2;
  int t_shift=0;
  if(W=="RV" || W=="FV") T_max= Nt/2 - 3;
  else T_max=Nt/2 -2;
  T_max_sx= (double)min(14, Nt/4);
  T_min_dx= T_max_sx+1;
  int ndof= (T_max -T_min+1 - 5);
  //choose the t-interval
  auto ansatz_sx= [&](const Vfloat& par, double t, int imeas) -> double { return par[0]*exp(-1*par[1]*t);};
  auto ansatz_dx= [&](const Vfloat& par, double t, int imeas) -> double { return par[0]*exp(-1*par[1]*(Nt/2-t));};
  auto ansatz_tot=[&](const Vfloat& par, double t, int imeas) -> double { return par[0]*exp(-1*par[1]*t) + par[2]*exp(-1*par[3]*(Nt/2-t));};
  auto ansatz = [=](const ExpFit_der& par, const Val& val) -> double {return par.a2*exp(-1*par.a3*val.t) +  par.a4*exp(-1*par.a5*(Nt/2-val.t));};
  Vfloat T_dx, T_sx, Y_dx, Y_sx, Y_err_dx, Y_err_sx, Y_dx_sx;
  for(int t=T_min;t<=T_max_sx;t++) {T_sx.push_back(t); Y_sx.push_back( F_der.ave()[t]); Y_err_sx.push_back( F_der.err()[t]);}
  for(int t=T_min_dx;t<=T_max;t++) {T_dx.push_back(t); Y_dx.push_back( F_der.ave()[t]); Y_err_dx.push_back( F_der.err()[t]);}
 
  double  a2_guess, a3_guess, a4_guess,a5_guess;
  pair<double,double> fit_interval_a3;
  pair<double, double> fit_interval_a5;
  int sign_a2, sign_a4;


  //choose guess intervals depending on the observable
  //###############################################################
  fit_interval_a3= make_pair(0.01, 1.0);
  fit_interval_a5= make_pair(0.01, 1.0);
  if( (W=="RV" || W=="FV")) {
    if(smearing) {
      if(W=="FV") {
	sign_a2=1;
	sign_a4=1;
      }
      else {
	sign_a2=1;
	sign_a4=1;
      }
    }
    else {
      if(W=="FV") {
	sign_a2=-1;
	sign_a4=1;
      }
      else {
	sign_a2=-1;
	sign_a4=-1;
      }

    }
  }
  else {
    if(W=="FA") {
      sign_a2=-1;
      sign_a4=-1;
    }
    else {
      if(smearing) {
	sign_a2=+1;
	sign_a4=1;
      }
      else {
	sign_a2=+1;
	sign_a4=1;
      }
    }
  }
  //################################################################
  
  //SX_FIT
    UniformMersenne  A2(6545534,log(0.001), log(10.0));
    UniformMersenne A4(87756756,log(0.001), log(10.0));
    UniformMersenne A3(659878,fit_interval_a3.first, fit_interval_a3.second);
    UniformMersenne A5(546353,fit_interval_a5.first, fit_interval_a5.second);
    int its_sx=0;
    bool sx_min_found=false;
    Vfloat Ch2_sx_list;
    bool sign_sx_flipped=false;
    while(!sx_min_found) {
      its_sx++;
      T_fit fit_sx(T_sx, Y_sx, Y_err_sx);
      fit_sx.ansatz = ansatz_sx;
      int ndof = (T_max_sx-T_min +1 - 3);
      double sx_1 = sign_a2*exp( A2());
      double sx_2 = A3();
      fit_sx.add_pars({sx_1, sx_2});
      fit_sx.add_par_errs({sx_1/10, sx_2/10});
      fit_t_res sx_fit_res = fit_sx.fit();
      Ch2_sx_list.push_back(sx_fit_res.chi2/ndof);
      if(sx_fit_res.chi2/ndof < 1) { cout<<"ch2_sx der for OBS: "<<W<<" xg: "<<xg<<" ch2: "<<sx_fit_res.chi2/ndof<<"  ITs: "<<its_sx<<endl; a2_guess= sx_fit_res.pars[0]; a3_guess=sx_fit_res.pars[1]; sx_min_found=true;}
      if (its_sx % 100 == 0) cout<<"IT(sx) der: "<<its_sx<<" OBS: "<<W<<" ,chi2/ndof: "<<sx_fit_res.chi2/ndof<<endl;
      if( (its_sx > 2000 && !sign_sx_flipped) || (its_sx > 7000 && sign_sx_flipped)) {
	//cout<<"ch2_sx for OBS: "<<W<<" xg: "<<xg<<" ch2: "<<sx_fit_res.chi2/ndof<<" ITs: "<<its_sx<<endl;
	//cout<<"do you want to continue?"<<endl;
	//cin >> sx_min_found;
	double actual_ch2 = Ch2_sx_list[Ch2_sx_list.size()-1];
	double threshold = Ch2_sx_list.size()*0.75;
	int count_ch2=0;
	if(actual_ch2 < 4) {for( auto &c : Ch2_sx_list) { if( fabs(c-actual_ch2) < 1.0e-3) count_ch2++;}}
	if((double)count_ch2 > threshold) sx_min_found=true; //ch2 is bad but seems that Minuit found the right minimum
	if(sx_min_found) { a2_guess=sx_fit_res.pars[0]; a3_guess=sx_fit_res.pars[1];}
      }
      if(its_sx > 5000 && !sign_sx_flipped) {sign_a2 *= -1; Ch2_sx_list.clear(); sign_sx_flipped=true;  } //change sign
      if(its_sx > 10000) return;  //not able to perform the fit
    }

 

    //DX_FIT
    bool dx_min_found=false;
    int its_dx=0;
    Vfloat Ch2_dx_list;
    bool sign_dx_flipped=false;
    while(!dx_min_found) {
      its_dx++;
      T_fit fit_dx(T_dx, Y_dx, Y_err_dx);
      fit_dx.ansatz= ansatz_dx;
      int ndof = (T_max-T_min_dx +1 - 3);
      double dx_1= sign_a4*exp(A4());
      double dx_2= A5();
      fit_dx.add_pars({dx_1, dx_2});
      fit_dx.add_par_errs({ dx_1/10, dx_2/10});
      fit_t_res dx_fit_res = fit_dx.fit();
      Ch2_dx_list.push_back(dx_fit_res.chi2/ndof);
      if(dx_fit_res.chi2/ndof < 1) { cout<<"ch2_dx der for OBS: "<<W<<" xg: "<<xg<<" ch2: "<<dx_fit_res.chi2/ndof<<" ITs: "<<its_dx<<endl; a4_guess= dx_fit_res.pars[0]; a5_guess=dx_fit_res.pars[1]; dx_min_found=true;}
      if (its_dx % 100 == 0) cout<<"IT(dx) der: "<<its_dx<<" OBS: "<<W<<" ,chi2/ndof: "<<dx_fit_res.chi2/ndof<<endl;
      if((its_dx > 2000 && !sign_dx_flipped) || (its_dx > 7000 && sign_dx_flipped)) {
	//cout<<"ch2_dx for OBS: "<<W<<" xg: "<<xg<<" ch2: "<<dx_fit_res.chi2/ndof<<" ITs: "<<its_dx<<endl;
	//cout<<"do you want to continue?"<<endl;
	//cin >> dx_min_found;
	double actual_ch2 = Ch2_dx_list[Ch2_dx_list.size()-1];
	double threshold = Ch2_dx_list.size()*0.75;
	int count_ch2=0;
	if(actual_ch2 < 4) {for( auto &c : Ch2_dx_list) { if( fabs(c-actual_ch2) < 1.0e-3) count_ch2++;}}
	if((double)count_ch2 > threshold) dx_min_found=true;
	if(dx_min_found) { a4_guess=dx_fit_res.pars[0]; a5_guess=dx_fit_res.pars[1];}
    
      }
      if(its_dx > 5000 && !sign_dx_flipped) {sign_a4*= -1; Ch2_dx_list.clear(); sign_dx_flipped=true;}
      if(its_dx > 10000) return; //not able to perform the fit 
    }

  
    //BOOTSTRAP FIT BASED ON JACKKNIFE RESAMPLING
    
    //##############################################################################
    bootstrap_fit<ExpFit_der, Val> fit(Njacks);
    ndof = T_max -T_min + 1- 5 -t_shift;
    fit.measurement = [=](const ExpFit_der& par, const Val& val) -> double { return val.meas;};
    fit.error = [=](const ExpFit_der& par, const Val& val) -> double {return val.err;};
    fit.ansatz =  ansatz;
    fit.Set_number_of_measurements(T_max - T_min +1-t_shift);
    fit.Set_verbosity(verbosity);
  
    cout<<"a2_g_der: "<<a2_guess<<endl;
    cout<<"a3_g_der: "<<a3_guess<<endl;
    cout<<"a4_g_der: "<<a4_guess<<endl;
    cout<<"a5_g_der: "<<a5_guess<<endl;
    fit.Add_par("a2", a2_guess, fabs(a2_guess)/30);
    fit.Add_par("a3", a3_guess, fabs(a3_guess)/30);
    fit.Add_par("a4", a4_guess, fabs(a4_guess)/30);
    fit.Add_par("a5", a5_guess, fabs(a5_guess)/30);
  
  
    //generate bootstrap sample
    for(int iboot=0;iboot<Njacks;iboot++) {
      fit.ib =&iboot;
      for(int i=T_min+t_shift;i<=T_max;i++) {
	double meas= F_der.distr_list[i].distr[iboot];
	double err = F_der.err()[i];
	Val X;
	X.t= (double)i;
	X.meas=meas;
	X.err=err;
	fit.Append_to_input_par(X);
      }
    }
    
  
  
    //Perform bootstrap fit
    boot_fit_data<ExpFit_der> Fit_output= fit.Perform_bootstrap_fit();
  
  
  
    //retrieve output params
    //##############################################
  
    double ave_ch2=0;
    Vfloat a2_der, a3_der, a4_der, a5_der;
 
  
    for(int iboot=0;iboot<Njacks;iboot++) {
      a2_der.push_back(Fit_output.par[iboot].a2);
      a3_der.push_back(Fit_output.par[iboot].a3);
      a4_der.push_back(Fit_output.par[iboot].a4);
      a5_der.push_back(Fit_output.par[iboot].a5);
      ave_ch2+= Fit_output.chi2[iboot]/(Njacks);
    }

    //determine C(t) - a1 from fit and substract it to F
    distr_t_list F_sub(UseJack, Nt);

    for(int t=0;t<Nt;t++) {
      for(int ijack=0;ijack<Njacks;ijack++) {

	double a2 = -1.0*Fit_output.par[ijack].a2/Fit_output.par[ijack].a3;
	double a3 = Fit_output.par[ijack].a3;
	double a4 = Fit_output.par[ijack].a4/Fit_output.par[ijack].a5;
	double a5 = Fit_output.par[ijack].a5;

	F_sub.distr_list[t].distr.push_back( ansatz_tot({a2,a3,a4,a5}, t, 0));
      }
    }


    cout<<"F size: "<<F.distr_list[0].size()<<endl;
    cout<<"F_sub size: "<<F_sub.distr_list[0].size()<<endl;
    distr_t_list improved_constant = F- F_sub;

    boost::filesystem::create_directory("../data/form_factors/"+Meson+"/Contaminations/improved_constant");

    Print_To_File({}, {improved_constant.ave(), improved_constant.err()}, "../data/form_factors/"+Meson+"/Contaminations/improved_constant/"+W+"_"+Ens_tag+"_xg_"+to_string_with_precision(xg,4)+"_offsh_"+to_string_with_precision(offsh,4)+"_smlev_"+to_string(smearing)+"_k0_noise_"+to_string(Include_k0_noise)+".dat", "", "");
    
    return;
}


















void Get_Tmin_Tmax(string corr_type, string Ens_tag, CorrAnalysis &corr, double xg, string W) {

    if(corr_type=="2pt") {
      corr.Tmax = 20;
      if(Ens_tag.substr(0,1) =="A" || Ens_tag.substr(0,1) == "C" || Ens_tag.substr(0,1) == "B" || Ens_tag.substr(0,1)=="D" || Ens_tag.substr(0,1) =="E" ) corr.Tmin = 13;

    }

    if(corr_type=="3pt") {
      if(Ens_tag.substr(0,6)=="A40.32" || Ens_tag.substr(0,6) == "C40.32") {
	if(W=="A") { 
	  if( xg > -0.05 && xg < 0.05)  {corr.Tmax = 23; corr.Tmin= 18;}
	  else if( xg > 0.09 && xg < 0.11)  {corr.Tmax = 20; corr.Tmin= 15;}
	  else if( xg > 0.19 && xg < 0.21)  {corr.Tmax = 20; corr.Tmin=14;}
	  else if (xg > 0.39 && xg < 0.41)  {corr.Tmax = 19; corr.Tmin= 14;}
	  else if( xg > 0.59 && xg < 0.61)  {corr.Tmax = 16; corr.Tmin= 13;}
	  else if( xg > 0.79 && xg < 0.81)  {corr.Tmax = 15; corr.Tmin= 11;}
	  else if( xg > 0.9 && xg < 1.1)  {corr.Tmax = 15; corr.Tmin= 11;}
	  else {corr.Tmin= 9; corr.Tmax=16;}
	}
	else if (W=="V") {
	  if( xg > -0.05 && xg < 0.05)  {corr.Tmax = 23; corr.Tmin= 18;}
	  else if (xg >0.09 && xg < 0.11) { corr.Tmax = 20; corr.Tmin= 15;}
	  else if( xg > 0.19 && xg < 0.21) { corr.Tmax = 23; corr.Tmin= 15;}
	  else if (xg > 0.39 && xg < 0.41)  {corr.Tmax = 20; corr.Tmin= 15;}
	  else if( xg > 0.59 && xg < 0.61)  {corr.Tmax = 25; corr.Tmin= 12;}
	  else if( xg > 0.79 && xg < 0.81)  {corr.Tmax = 23; corr.Tmin= 12;}
	  else if( xg > 0.9 && xg < 1.1)  {corr.Tmax = 17; corr.Tmin= 13;}
	  else {corr.Tmin= 9; corr.Tmax=16;}
	}
      
      }
      else  if(Ens_tag.substr(0,6)=="B40.32") {
	if(W=="A") { 
	  if( xg > -0.05 && xg < 0.05)  {corr.Tmax = 23; corr.Tmin= 18;}
	  else if( xg > 0.09 && xg < 0.11)  {corr.Tmax = 17; corr.Tmin= 13;}
	  else if( xg > 0.19 && xg < 0.21)  {corr.Tmax = 16; corr.Tmin=12;}
	  else if (xg > 0.39 && xg < 0.41)  {corr.Tmax = 17; corr.Tmin= 13;}
	  else if( xg > 0.59 && xg < 0.61)  {corr.Tmax = 16; corr.Tmin= 13;}
	  else if( xg > 0.79 && xg < 0.81)  {corr.Tmax = 15; corr.Tmin= 11;}
	  else if( xg > 0.9 && xg < 1.1)  {corr.Tmax = 15; corr.Tmin= 11;}
	  else {corr.Tmin= 9; corr.Tmax=16;}
	}
	else if (W=="V") {
	  if( xg > -0.05 && xg < 0.05)  {corr.Tmax = 23; corr.Tmin= 18;}
	  else if (xg >0.09 && xg < 0.11) { corr.Tmax = 20; corr.Tmin= 15;}
	  else if( xg > 0.19 && xg < 0.21) { corr.Tmax = 20; corr.Tmin= 15;}
	  else if (xg > 0.39 && xg < 0.41)  {corr.Tmax = 20; corr.Tmin= 12;}
	  else if( xg > 0.59 && xg < 0.61)  {corr.Tmax = 16; corr.Tmin= 12;}
	  else if( xg > 0.79 && xg < 0.81)  {corr.Tmax = 16; corr.Tmin= 11;}
	  else if( xg > 0.9 && xg < 1.1)  {corr.Tmax = 16; corr.Tmin= 11;}
	  else {corr.Tmin= 9; corr.Tmax=16;}
	}
      
      }
      else if(Ens_tag.substr(0,6)=="A40.24" || Ens_tag.substr(0,6)=="B40.24" || Ens_tag.substr(0,6)=="C40.24") {
	if(W=="A") {
	  if( xg > -0.05 && xg < 0.05)  {corr.Tmax = 23; corr.Tmin= 18;}
	  else if( xg > 0.19 && xg < 0.21)  {corr.Tmax = 16; corr.Tmin= 9;}
	  else if (xg > 0.39 && xg < 0.41)  {corr.Tmax = 17; corr.Tmin= 9;}
	  else if( xg > 0.59 && xg < 0.61)  {corr.Tmax = 14; corr.Tmin= 8;}
	  else if( xg > 0.79 && xg < 0.81)  {corr.Tmax = 13; corr.Tmin= 8;}
	  else if( xg > 0.9 && xg < 1.1)  {corr.Tmax = 14; corr.Tmin= 8;}
	  else {corr.Tmin= 9; corr.Tmax=16;}
	}
	else if (W=="V") {
	  if( xg > -0.05 && xg < 0.05)  {corr.Tmax = 23; corr.Tmin= 18;}
	  else if( xg > 0.09 && xg < 0.11) { corr.Tmax = 18; corr.Tmin= 14;}
	  else if( xg > 0.19 && xg < 0.21) { corr.Tmax = 18; corr.Tmin= 14;}
	  else if (xg > 0.39 && xg < 0.41)  {corr.Tmax = 17; corr.Tmin= 12;}
	  else if( xg > 0.59 && xg < 0.61)  {corr.Tmax = 13; corr.Tmin= 12;}
	  else if( xg > 0.79 && xg < 0.81)  {corr.Tmax = 13; corr.Tmin= 10;}
	  else if( xg > 0.9 && xg < 1.1)  {corr.Tmax = 13; corr.Tmin= 10;}
	  else {corr.Tmin= 9; corr.Tmax=16;}
	}
	else crash("cannot find the ensemble");
      }
      else if(Ens_tag.substr(0,6)=="A40.40" || Ens_tag.substr(0,6)=="B40.40") {
	if(W=="A") {
	  if( xg > -0.05 && xg < 0.05)  {corr.Tmax = 25; corr.Tmin= 18;}
	  else if( xg > 0.09 && xg < 0.11)  {corr.Tmax = 20; corr.Tmin= 15;}
	  else if( xg > 0.19 && xg < 0.21)  {corr.Tmax = 20; corr.Tmin= 15;}
	  else if (xg > 0.39 && xg < 0.41)  {corr.Tmax = 20; corr.Tmin= 15;}
	  else if( xg > 0.59 && xg < 0.61)  {corr.Tmax = 20; corr.Tmin= 15;}
	  else if( xg > 0.79 && xg < 0.81)  {corr.Tmax = 20; corr.Tmin= 13;}
	  else if( xg > 0.9 && xg < 1.1)  {corr.Tmax = 22; corr.Tmin= 13;}
	  else {corr.Tmin= 9; corr.Tmax=16;}
	}
	else if (W=="V") {
	  if( xg > -0.05 && xg < 0.05)  {corr.Tmax = 24; corr.Tmin= 18;}
	  else if( xg > 0.09 && xg < 0.11)  {corr.Tmax = 24; corr.Tmin= 18;}
	  else if( xg > 0.19 && xg < 0.21) { corr.Tmax = 24; corr.Tmin= 18;}
	  else if (xg > 0.39 && xg < 0.41)  {corr.Tmax = 24; corr.Tmin= 18;}
	  else if( xg > 0.59 && xg < 0.61)  {corr.Tmax = 25; corr.Tmin= 18;}
	  else if( xg > 0.79 && xg < 0.81)  {corr.Tmax = 24; corr.Tmin= 15;}
	  else if( xg > 0.9 && xg < 1.1)  {corr.Tmax = 23; corr.Tmin= 17;}
	  else {corr.Tmin= 9; corr.Tmax=16;}
	}
	else crash("cannot find the ensemble");
      }
      else if(Ens_tag.substr(0,6)=="A40.48") {
	if(W=="A") { corr.Tmax=20; corr.Tmin=14;}
	else if(W=="V") {corr.Tmax= 24; corr.Tmin=18;}
      }

    
      else { corr.Tmin= 9; corr.Tmax=16;}
    }

    return;

}




distr_t_list V_ave_unpolarized(vector<vector<distr_t_list>>& distr_mom_k, vector<vector<distr_t_list>>& distr_mom_0, distr_t& Meson_mass, pt3_momenta& Mom,int twall) {

  VVfloat e(4);
  e[0] = {0.0, 0.0, 0.0, 0.0};
  e[3] = {0.0, 0.0, 0.0, 0.0};
  e[1] = {0.0, -1.0/sqrt(2), -1.0/sqrt(2), 0.0};
  e[2] = {0.0, 1.0/sqrt(2), -1.0/sqrt(2), 0.0};

  



  if(distr_mom_k.size() != 4 || distr_mom_0.size() != 4) crash("V_ave_unpolarized called with a v<v<distr_t_list>> with size != alphas=4");
  int size= distr_mom_k[0][0].size();
  if(size==0) crash("In V_ave_unpolarized, distr_mom_k[0][0] has size zero, exiting....");
  int distr_size = distr_mom_k[0][0].distr_list[0].size();
  if(distr_size==0) crash("In V_ave_unpolarized, distr_mom_k[0][0] has zero distr_size, exiting....");


  distr_t_list V_ave(UseJack, size, distr_size);

  Vfloat exp_Nt_t;
  double Egamma=Mom.Egamma();
  double Egamma_T=sinh(Egamma)*(1-exp(-Egamma*Mom.Nt()));
  for(int t=0; t<Mom.Nt();t++) exp_Nt_t.push_back( exp( fabs((twall-t))*Egamma));
  //for(int t=0;t<Mom.Nt(); t++) exp_Nt_t.push_back( exp( ((Mom.Nt()/2.0)-t)*Mom.Egamma()) - exp((t- Mom.Nt()/2.0)*Mom.Egamma()));

  bool is_mom_zero=false;
  Vfloat theta_t_old(Mom.Theta(2));
  if(Mom.k_mod() < eps(5)) { is_mom_zero=true; Vfloat theta_t = {0.0, 0.0, 1.0}; Mom.Set_theta(theta_t, 2);}

  
  
  for(int alpha=1; alpha <= 2; alpha ++) {
   
    if(distr_mom_k[alpha].size() != 4 || distr_mom_0[alpha].size() != 4) crash("V_ave_unpolarized called with a v<v<distr_t_list>> having an element with size != 4");
    for(int mu=0; mu<4;mu++) {

      if(mu != alpha) {
    
	distr_t_list den1(3, Meson_mass);
	distr_t_list den2(3, Meson_mass);
	
	
	den1 = den1*Multiply_vector_by_scalar(external_prod(slicing(e[1],1,3),Mom.k()),-1.0) + Multiply_vector_by_scalar(external_prod(slicing(e[1],1,3), Mom.p()), Mom.Egamma());
    
	den2 = den2*Multiply_vector_by_scalar(external_prod(slicing(e[2],1,3),Mom.k()), -1.0) + Multiply_vector_by_scalar(external_prod(slicing(e[2],1,3), Mom.p()), Mom.Egamma());
    
	V_ave = V_ave+ e[1][mu]*(distr_mom_k[alpha][mu]*exp_Nt_t - ((double)Include_k0_noise)*distr_mom_0[alpha][mu])/den1[alpha-1];

	V_ave = V_ave+ e[2][mu]*(distr_mom_k[alpha][mu]*exp_Nt_t - ((double)Include_k0_noise)*distr_mom_0[alpha][mu])/den2[alpha-1];


      }

     
    }
  }


 
  
  
  distr_t_list DIFF2 = 2.0*( distr_mom_k[2][1]*exp_Nt_t -((double)Include_k0_noise)*distr_mom_0[2][1] - distr_mom_k[1][2]*exp_Nt_t + ((double)Include_k0_noise)*distr_mom_0[1][2]);
  distr_t_list DIFF= DIFF2/(Mom.k()[2]*Meson_mass);
  //for(int t=0; t < V_ave.size(); t++) cout<<V_ave.distr_list[t].ave()<<" "<<V_ave.distr_list[t].err()<<" "<<DIFF.distr_list[t].ave()<<" "<<DIFF.distr_list[t].err()<<endl;
  
  if(is_mom_zero) Mom.Set_theta(theta_t_old, 2);
  
  return V_ave;
  //return DIFF2;


}

distr_t_list A_ave_unpolarized(vector<vector<distr_t_list>>& distr_mom_k) {

  VVfloat e(4);
  e[0] = {0.0, 0.0, 0.0, 0.0};
  e[3] = {0.0, 0.0, 0.0, 0.0};
  e[1] = {0.0, -1.0/sqrt(2), -1.0/sqrt(2), 0.0};
  e[2] = {0.0, 1.0/sqrt(2), -1.0/sqrt(2), 0.0};



  if(distr_mom_k.size() != 4) crash("A_ave_unpolarized called with a v<v<distr_t_list>> with size != alphas= 4");
  for(auto & distr_mom_k_alpha : distr_mom_k) if(distr_mom_k_alpha.size() != 4) crash("A_ave unpolarized called with a with a vv<distr_t_list>> having an element of size != 4");
  int size = distr_mom_k[0][0].size();
  if(size==0) crash(" A_ave_unpolarized called with Nt=0");
  int distr_size = distr_mom_k[0][0].distr_list[0].size();
  if(distr_size ==0) crash(" A_ave unpolarized called with a vv<distr_t_list>> having zero statistics");
  distr_t_list A_ave(UseJack, size, distr_size);

 
  for(int alpha=1; alpha <= 2; alpha ++) {
 
    if(distr_mom_k[alpha].size() != 4) crash(" A_ave_unpolarized called with a v<v<distr_t_list>> having an element with size != 4");
 
    for(int mu=0; mu< 4; mu ++) {

 
      A_ave = A_ave+ distr_mom_k[alpha][mu]*(e[1][mu]/e[1][alpha]);

      A_ave = A_ave+ distr_mom_k[alpha][mu]*(e[2][mu]/e[2][alpha]);
    }
  }


  /*
    cout<<"A##########"<<endl;
    distr_t_list DIFF = 2.0*( distr_mom_k[1][1] +distr_mom_k[2][2]);
    for(int t=0; t < A_ave.size(); t++) cout<<A_ave.distr_list[t].ave()<<" "<<A_ave.distr_list[t].err()<<" "<<DIFF.distr_list[t].ave()<<" "<<DIFF.distr_list[t].err()<<endl;
    cout<<"############"<<endl;
  */

  //return distr_mom_k[1][1];
  return A_ave;
  
}


distr_t_list H_2(const distr_t_list& H30,const distr_t_list& H11,const distr_t_list& H33,const distr_t_list& H30_0,const distr_t_list& H11_0,const distr_t_list& H33_0,pt3_momenta& Mom,const distr_t& m) {//valid for p=0, k= kz

  auto Power = [&](const distr_t& A, int n) -> distr_t{
    distr_t ret_val = A;
    for(int i=1;i<n;i++) ret_val= ret_val*A;
    return ret_val;
  };

  double Eg= Mom.Egamma();
  double kz= Mom.k()[2];
  double ksq = pow(Mom.virt(),2);

  
  distr_t_list H30_sub = H30 - H11_0*(m-Eg)/( (2.0*m*Eg/kz) -(ksq/kz));
  distr_t_list H11_sub = H11 - H11_0 ;
  distr_t_list H33_sub = H33 - H33_0*(2.0*m*Eg -ksq-pow(kz,2))/(2.0*m*Eg -ksq);

  distr_t_list H_2 = +1.0*H30_sub/(Eg*kz*Power(m,2));
  H_2 = H_2 - H33_sub/(pow(kz,2)*Power(m,2));
  H_2 = H_2+ H11_sub/(pow(kz,2)*Power(m,2));
  return -1.0*m*H_2*(2.0*m*Eg-ksq  );

}


distr_t_list H_1(const distr_t_list& H30,const distr_t_list& H11,const distr_t_list& H33,const distr_t_list& H30_0,const distr_t_list& H11_0,const distr_t_list& H33_0,pt3_momenta& Mom,const distr_t& m) //valid for p=0, k=kz
{

  auto Power = [&](const distr_t& A, int n) -> distr_t{
    distr_t ret_val = A;
    for(int i=1;i<n;i++) ret_val= ret_val*A;
    return ret_val;
  };
  
  double Eg= Mom.Egamma();
  double kz= Mom.k()[2];
  double ksq = pow(Mom.virt(),2);

  
  distr_t_list H30_sub = H30 - H11_0*(m-Eg)/( (2.0*m*Eg/kz) -(ksq/kz));
  distr_t_list H11_sub = H11 - H11_0 ;
  distr_t_list H33_sub = H33 - H33_0*(2.0*m*Eg -ksq-pow(kz,2))/(2.0*m*Eg -ksq);


  distr_t_list H_1 = H30_sub*(ksq-Eg*m)/(Eg*kz*Power(m,2));
  H_1 = H_1+ H33_sub*( (-Eg+m)*(ksq-Eg*m)/(Eg*pow(kz,2)*Power(m,2)));
  H_1 = H_1+ H11_sub*( -pow(Eg,2)*m - (ksq+pow(kz,2))*m +Eg*(ksq +Power(m,2)))/(Eg*pow(kz,2)*Power(m,2));

  
  return H_1*m;



}


distr_t_list FA_off(const distr_t_list& H30,const distr_t_list& H11,const distr_t_list& H33,const distr_t_list& H30_0,const distr_t_list& H11_0,const distr_t_list& H33_0,pt3_momenta& Mom,const distr_t &m) //valid for p=0, k=kz
{

  double Eg= Mom.Egamma();
  double kz= Mom.k()[2]; 
  double ksq = pow(Mom.virt(),2);
  
  distr_t_list H30_sub = H30 - H11_0*(m-Eg)/( (2.0*m*Eg/kz) -(ksq/kz));
  distr_t_list H11_sub = H11 - H11_0 ;
  distr_t_list H33_sub = H33 - H33_0*(2.0*m*Eg -ksq-pow(kz,2))/(2.0*m*Eg -ksq);

  distr_t_list FA_off = H30_sub*ksq/(Eg*kz*m);
  FA_off = FA_off+ H33_sub*ksq*(m-Eg)/(Eg*pow(kz,2)*m);
  FA_off = FA_off+ H11_sub*( Eg*ksq - (ksq+pow(kz,2))*m)/(Eg*pow(kz,2)*m);
 
 
  return FA_off;
}

distr_t_list H_2_impr(const distr_t_list& H30,const distr_t_list& H03, const distr_t_list& H11,const distr_t_list& H33, const distr_t_list& H11_0,const distr_t_list& H33_0,pt3_momenta& Mom,const distr_t& m) {//valid for p=0, k= kz

  auto Power = [&](const distr_t& A, int n) -> distr_t{
    distr_t ret_val = A;
    for(int i=1;i<n;i++) ret_val= ret_val*A;
    return ret_val;
  };

  double Eg= Mom.Egamma();
  double kz= Mom.k()[2];
  double ksq = pow(Mom.virt(),2);

  CorrAnalysis C;
  C.Nt = 64;
  
  
  distr_t pole_corr;
  
  for(int ijack=0;ijack<Njacks;ijack++) pole_corr.distr.push_back( 0.5*( 1 + (m.distr[ijack] - Eg)/( 2*asinh( 0.5*sqrt(  kz*kz + pow( 2*sinh(m.distr[ijack]/2),2))))));


  distr_t corr_3003;
  for(int ijack=0;ijack < Njacks;ijack++) corr_3003.distr.push_back( -0.5*(kz*m.distr[ijack]/(2*m.distr[ijack] - Eg))*(1.0/( 2*asinh( 0.5*sqrt(  kz*kz + pow( 2*sinh(m.distr[ijack]/2),2))))));

  distr_t_list resc_FSES;
  distr_t mk_k;
  for(int ijack=0;ijack<Njacks;ijack++) mk_k.distr.push_back(  2*asinh( 0.5*sqrt(  kz*kz + pow( 2*sinh(m.distr[ijack]/2),2))) + Eg - m.distr[ijack]);
  auto exp_M= [](double mass, double t, double N) { return exp( -mass*t);};
  distr_t_list exp_mk = distr_t_list::f_of_distr(exp_M, mk_k, C.Nt);
  resc_FSES= 1- pole_corr*(1.0*Include_FSES)*exp_mk;

  



 distr_t_list H30_n = H30/H11_0;
 distr_t_list H03_n = H03/H11_0;
 distr_t_list H11_n = H11/H11_0;
 distr_t_list H33_n = H33/H11_0;
 distr_t_list H11_0_n = H11_0/H11_0;
 distr_t_list H33_0_n = H33_0/H11_0;

 distr_t den = pow(Eg*kz,2)*Power(m,2)*(Eg-2.0*m);
  
 distr_t_list H30_sub = H30_n - H03_n*(m-Eg)/(2.0*m-Eg) - Include_FSES*corr_3003*exp_mk ;
 distr_t_list H11_sub = H11_n - H11_0_n ;
 distr_t_list H33_sub = H33_n - H33_0_n*(2.0*m*Eg -ksq-resc_FSES*pow(kz,2))/(2.0*m*Eg -ksq);
 
 distr_t_list H_2 = +1.0*H30_sub/(Eg*kz*Power(m,2));
 H_2 = H_2 + H33_sub*( (pow(Eg,2)+ksq)*m - Eg*ksq)/den;
 H_2 = H_2+ H11_sub*( Eg*(ksq+pow(kz,2))-1.0*(pow(Eg,2)+ ksq+pow(kz,2))*m)/den;


  
 
  Get_virtual_ff_fit_interval_v1("H2", Mom.virt(), Mom.Theta(2)[2], C.Tmin, C.Tmax);
  distr_t_list contr_1 = (H30_sub/(Eg*kz*Power(m,2)));
  distr_t_list contr_2 = (H33_sub*( (pow(Eg,2)+ksq)*m - Eg*ksq)/den);
  distr_t_list contr_3 = (H11_sub*( Eg*(ksq+pow(kz,2))-1.0*(pow(Eg,2)+ ksq+pow(kz,2))*m)/den);
  distr_t_list contr_11_33 = ( ((pow(Eg,2)+ksq)*m - Eg*ksq)/den   +   ( Eg*(ksq+pow(kz,2))-1.0*(pow(Eg,2)+ ksq+pow(kz,2))*m)/den      )*(H33_sub + H11_sub)/2.0;
  distr_t_list contr_33_m_11 =    ( ((pow(Eg,2)+ksq)*m - Eg*ksq)/den   -   ( Eg*(ksq+pow(kz,2))-1.0*(pow(Eg,2)+ ksq+pow(kz,2))*m)/den      )*(H33_sub - H11_sub)/2.0;
  distr_t_list total= H_2;
  cout<<"########## H2 impr infos: "<<endl;
  cout<<"contr 30: "<<C.Fit_distr(contr_1).ave()<<" +- "<<C.Fit_distr(contr_1).err()<<endl;
  cout<<"contr 33: "<<C.Fit_distr(contr_2).ave()<<" +- "<<C.Fit_distr(contr_2).err()<<endl;
  cout<<"contr 11: "<<C.Fit_distr(contr_3).ave()<<" +- "<<C.Fit_distr(contr_3).err()<<endl;
  cout<<"contr 33+11: "<<C.Fit_distr(contr_11_33).ave()<<" +- "<<C.Fit_distr(contr_11_33).err()<<endl;
  cout<<"contr 33-11: "<<C.Fit_distr(contr_33_m_11).ave()<<" +- "<<C.Fit_distr(contr_33_m_11).err()<<endl;
  cout<<"total: "<<C.Fit_distr(total).ave()<<" +- "<<C.Fit_distr(total).err()<<endl;
  
  cout<<"########## END H2 impr infos "<<endl;



  
  return -1.0*m*H_2*(2.0*m*Eg-ksq  );

}


distr_t_list H_1_impr(const distr_t_list& H30,const distr_t_list& H03, const distr_t_list& H11,const distr_t_list& H33, const distr_t_list& H11_0,const distr_t_list& H33_0,pt3_momenta& Mom,const distr_t& m) //valid for p=0, k=kz
{
  
  double Eg= Mom.Egamma();
  double kz= Mom.k()[2]; 
  double ksq = pow(Mom.virt(),2);

  auto Power = [&](const distr_t& A, int n) -> distr_t{
    distr_t ret_val = A;
    for(int i=1;i<n;i++) ret_val= ret_val*A;
    return ret_val;
  };

  
  CorrAnalysis C;
  C.Nt = 64;

  distr_t pole_corr;
  for(int ijack=0;ijack<Njacks;ijack++) pole_corr.distr.push_back( 0.5*( 1 + (m.distr[ijack] - Eg)/( 2*asinh( 0.5*sqrt(  kz*kz + pow( 2*sinh(m.distr[ijack]/2),2))))));

  distr_t corr_3003;
  for(int ijack=0;ijack < Njacks;ijack++) corr_3003.distr.push_back( -0.5*(kz*m.distr[ijack]/(2*m.distr[ijack] - Eg))*(1.0/( 2*asinh( 0.5*sqrt(  kz*kz + pow( 2*sinh(m.distr[ijack]/2),2))))));

  distr_t_list resc_FSES;
  distr_t mk_k;
  for(int ijack=0;ijack<Njacks;ijack++) mk_k.distr.push_back(  2*asinh( 0.5*sqrt(  kz*kz + pow( 2*sinh(m.distr[ijack]/2),2))) + Eg -m.distr[ijack]);
  auto exp_M= [](double mass, double t, double N) { return exp( -mass*t);};
  distr_t_list exp_mk = distr_t_list::f_of_distr(exp_M, mk_k, C.Nt);
  resc_FSES= 1-pole_corr*(1.0*Include_FSES)*exp_mk;

  distr_t_list H30_n = H30/H11_0;
  distr_t_list H03_n = H03/H11_0;
  distr_t_list H11_n = H11/H11_0;
  distr_t_list H33_n = H33/H11_0;
  distr_t_list H11_0_n = H11_0/H11_0;
  distr_t_list H33_0_n = H33_0/H11_0;
  
  distr_t den = pow(Eg*kz,2)*Power(m,2)*(Eg-2.0*m);
  
  
  distr_t_list H30_sub = H30_n - H03_n*(m-Eg)/(2.0*m-Eg) - Include_FSES*corr_3003*exp_mk;
  distr_t_list H11_sub = H11_n - H11_0_n ;
  distr_t_list H33_sub = H33_n - H33_0_n*(2.0*m*Eg -ksq-resc_FSES*pow(kz,2))/(2.0*m*Eg -ksq);


  distr_t_list H_1 = H30_sub*(ksq - Eg*m)/(Eg*kz*Power(m,2));
  H_1 = H_1 -H33_sub*(Eg-m)*(Eg*m-ksq)*(-ksq+2.0*Eg*m)/den;
  H_1 = H_1+ H11_sub*(-1.0*ksq*m*(ksq +pow(kz,2)) + 2.0*pow(Eg,3)*Power(m,2) + Eg*(ksq+pow(kz,2))*(ksq+3.0*Power(m,2)) -pow(Eg,2)*m*(3.0*ksq + 2.0*(pow(kz,2)+Power(m,2))))/den;


  Get_virtual_ff_fit_interval_v1("H1", Mom.virt(), Mom.Theta(2)[2], C.Tmin, C.Tmax);
  distr_t_list contr_1 = H30_sub*((ksq - Eg*m)/(Eg*kz*Power(m,2)));
  distr_t_list contr_2 = (-1.0*H33_sub*(Eg-m)*(Eg*m-ksq)*(-ksq+2.0*Eg*m)/den);
  distr_t_list contr_3 = (H11_sub*(-1.0*ksq*m*(ksq +pow(kz,2)) + 2.0*pow(Eg,3)*Power(m,2) + Eg*(ksq+pow(kz,2))*(ksq+3.0*Power(m,2)) -pow(Eg,2)*m*(3.0*ksq + 2.0*(pow(kz,2)+Power(m,2))))/den);
  distr_t_list contr_11_33 = ( -1.0*(Eg-m)*(Eg*m-ksq)*(-ksq+2.0*Eg*m)/den  +  (-1.0*ksq*m*(ksq +pow(kz,2)) + 2.0*pow(Eg,3)*Power(m,2) + Eg*(ksq+pow(kz,2))*(ksq+3.0*Power(m,2)) -pow(Eg,2)*m*(3.0*ksq + 2.0*(pow(kz,2)+Power(m,2))))/den)*(H33_sub + H11_sub)/2.0;
  distr_t_list contr_33_m_11 =   ( -1.0*(Eg-m)*(Eg*m-ksq)*(-ksq+2.0*Eg*m)/den -  (-1.0*ksq*m*(ksq +pow(kz,2)) + 2.0*pow(Eg,3)*Power(m,2) + Eg*(ksq+pow(kz,2))*(ksq+3.0*Power(m,2)) -pow(Eg,2)*m*(3.0*ksq + 2.0*(pow(kz,2)+Power(m,2))))/den)*(H33_sub - H11_sub)/2.0;
  distr_t_list total= H_1;
  cout<<"########## H1 impr infos: "<<endl;
  cout<<"contr 30: "<<C.Fit_distr(contr_1).ave()<<" +- "<<C.Fit_distr(contr_1).err()<<endl;
  cout<<"contr 33: "<<C.Fit_distr(contr_2).ave()<<" +- "<<C.Fit_distr(contr_2).err()<<endl;
  cout<<"contr 11: "<<C.Fit_distr(contr_3).ave()<<" +- "<<C.Fit_distr(contr_3).err()<<endl;
  cout<<"contr 33+11: "<<C.Fit_distr(contr_11_33).ave()<<" +- "<<C.Fit_distr(contr_11_33).err()<<endl;
  cout<<"contr 33-11: "<<C.Fit_distr(contr_33_m_11).ave()<<" +- "<<C.Fit_distr(contr_33_m_11).err()<<endl;
  cout<<"total: "<<C.Fit_distr(total).ave()<<" +- "<<C.Fit_distr(total).err()<<endl;
  cout<<"########## END H1 impr infos "<<endl;

  
  return H_1*m;



}


distr_t_list FA_off_impr(const distr_t_list& H30,const distr_t_list& H03, const distr_t_list& H11,const distr_t_list& H33,const distr_t_list& H11_0,const distr_t_list& H33_0,pt3_momenta& Mom,const distr_t &m) //valid for p=0, k=kz
{

  double Eg= Mom.Egamma();
  double kz= Mom.k()[2]; 
  double ksq = pow(Mom.virt(),2);

  auto Power = [&](const distr_t& A, int n) -> distr_t{
    distr_t ret_val = A;
    for(int i=1;i<n;i++) ret_val= ret_val*A;
    return ret_val;
  };

  
  CorrAnalysis C;
  C.Nt = 64;

  distr_t pole_corr;
  for(int ijack=0;ijack<Njacks;ijack++) pole_corr.distr.push_back( 0.5*( 1 + (m.distr[ijack] - Eg)/( 2*asinh( 0.5*sqrt(  kz*kz + pow( 2*sinh(m.distr[ijack]/2),2))))));

  distr_t corr_3003;
  for(int ijack=0;ijack < Njacks;ijack++) corr_3003.distr.push_back( -0.5*(kz*m.distr[ijack]/(2*m.distr[ijack] - Eg))*(1.0/( 2*asinh( 0.5*sqrt(  kz*kz + pow( 2*sinh(m.distr[ijack]/2),2))))));

  distr_t_list resc_FSES;
  distr_t mk_k;
  for(int ijack=0;ijack<Njacks;ijack++) mk_k.distr.push_back(  2*asinh( 0.5*sqrt(  kz*kz + pow( 2*sinh(m.distr[ijack]/2),2))) + Eg - m.distr[ijack]);
  auto exp_M= [](double mass, double t, double N) { return exp( -mass*t);};
  distr_t_list exp_mk = distr_t_list::f_of_distr(exp_M, mk_k, C.Nt);
  resc_FSES= 1-pole_corr*(1.0*Include_FSES)*exp_mk;

  distr_t_list H0_ave = (2.0*H11_0 + 1.0*H33_0)/3.0;

  double xk = (Mom.virt()/m).ave();

  distr_t_list H_to_use = (xk > 0.46)?H11_0:H33_0;

  distr_t_list H30_n = H30/H11_0;
  distr_t_list H03_n = H03/H11_0;
  distr_t_list H11_n = H11/H11_0;
  distr_t_list H33_n = H33/H11_0;
  distr_t_list H11_0_n = H11_0/H11_0;
  distr_t_list H33_0_n = H_to_use/H11_0;

  distr_t den = pow(Eg*kz,2)*Power(m,2)*(Eg-2.0*m);
   
  distr_t_list H30_sub = H30_n - H03_n*(m-Eg)/(2.0*m-Eg) - Include_FSES*corr_3003*exp_mk;
  distr_t_list H11_sub = H11_n - H11_0_n ;
  distr_t_list H33_sub = H33_n - H33_0_n*(2.0*m*Eg -ksq-resc_FSES*pow(kz,2))/(2.0*m*Eg -ksq);  


  distr_t_list FA_off = H30_sub*ksq/(Eg*kz*m);
  FA_off = FA_off+ H33_sub*m*ksq*(m-Eg)*(ksq-2.0*Eg*m)/den;
  FA_off = FA_off+ H11_sub*m*( -ksq*m*(ksq+pow(kz,2)) -pow(Eg,2)*m*(2.0*ksq+pow(kz,2)) +Eg*(ksq + pow(kz,2))*(ksq+2.0*Power(m,2)))/den;


 
  Get_virtual_ff_fit_interval_v1("FA", Mom.virt(), Mom.Theta(2)[2], C.Tmin, C.Tmax);
  distr_t_list contr_1 = (H30_sub*ksq/(Eg*kz*m));
  distr_t_list contr_2 = (H33_sub*m*ksq*(m-Eg)*(ksq-2.0*Eg*m)/den);
  distr_t_list contr_3 = (H11_sub*m*( -ksq*m*(ksq+pow(kz,2)) -pow(Eg,2)*m*(2.0*ksq+pow(kz,2)) +Eg*(ksq + pow(kz,2))*(ksq+2.0*Power(m,2)))/den);
  distr_t_list contr_11_33 = ((m*ksq*(m-Eg)*(ksq-2.0*Eg*m) + m*( -ksq*m*(ksq+pow(kz,2)) -pow(Eg,2)*m*(2.0*ksq+pow(kz,2)) +Eg*(ksq + pow(kz,2))*(ksq+2.0*Power(m,2))))/den)*(H33_sub + H11_sub)/2.0;
  distr_t_list contr_33_m_11 =  ((m*ksq*(m-Eg)*(ksq-2.0*Eg*m) - 1.0*m*( -ksq*m*(ksq+pow(kz,2)) -pow(Eg,2)*m*(2.0*ksq+pow(kz,2)) +Eg*(ksq + pow(kz,2))*(ksq+2.0*Power(m,2))))/den)*(H33_sub - H11_sub)/2.0;
  distr_t_list total= FA_off;
  cout<<"########## FA impr infos: "<<endl;
  cout<<"contr 30: "<<C.Fit_distr(contr_1).ave()<<" +- "<<C.Fit_distr(contr_1).err()<<endl;
  cout<<"contr 33: "<<C.Fit_distr(contr_2).ave()<<" +- "<<C.Fit_distr(contr_2).err()<<endl;
  cout<<"contr 11: "<<C.Fit_distr(contr_3).ave()<<" +- "<<C.Fit_distr(contr_3).err()<<endl;
  cout<<"contr 33+11: "<<C.Fit_distr(contr_11_33).ave()<<" +- "<<C.Fit_distr(contr_11_33).err()<<endl;
  cout<<"contr 33-11: "<<C.Fit_distr(contr_33_m_11).ave()<<" +- "<<C.Fit_distr(contr_33_m_11).err()<<endl;
  cout<<"total: "<<C.Fit_distr(total).ave()<<" +- "<<C.Fit_distr(total).err()<<endl;
  cout<<"########## END FA impr infos "<<endl;
 
  return FA_off;
}


distr_t_list H_2_mixed_diag(const distr_t_list& H30,const distr_t_list& H03, const distr_t_list& H11,const distr_t_list& H33, const distr_t_list& H11_0,const distr_t_list& H33_0, const distr_t_list& H33_kz0, const distr_t_list& H11_kz0, pt3_momenta& Mom,const distr_t& m) {//valid for p=0, k= kz

  auto Power = [&](const distr_t& A, int n) -> distr_t{
    distr_t ret_val = A;
    for(int i=1;i<n;i++) ret_val= ret_val*A;
    return ret_val;
  };

  double Eg= Mom.Egamma();
  double kz= Mom.k()[2];
  double ksq = pow(Mom.virt(),2);

  CorrAnalysis C;
  C.Nt = 64;

  distr_t pole_corr;
  for(int ijack=0;ijack<Njacks;ijack++) pole_corr.distr.push_back( 0.5*( 1 + (m.distr[ijack] - Eg)/( 2*asinh( 0.5*sqrt(  kz*kz + pow( 2*sinh(m.distr[ijack]/2),2))))));

  distr_t corr_3003;
  for(int ijack=0;ijack < Njacks;ijack++) corr_3003.distr.push_back( -0.5*(kz*m.distr[ijack]/(2*m.distr[ijack] - Eg))*(1.0/( 2*asinh( 0.5*sqrt(  kz*kz + pow( 2*sinh(m.distr[ijack]/2),2))))));


  distr_t_list resc_FSES;
  distr_t mk_k;
  for(int ijack=0;ijack<Njacks;ijack++) mk_k.distr.push_back( 2*asinh( 0.5*sqrt(  kz*kz + pow( 2*sinh(m.distr[ijack]/2),2))) + Eg - m.distr[ijack]);
  auto exp_M= [](double mass, double t, double N) { return exp( -mass*t);};
  distr_t_list exp_mk = distr_t_list::f_of_distr(exp_M, mk_k, C.Nt);
  resc_FSES= 1-pole_corr*(1.0*Include_FSES)*exp_mk;

  

  distr_t_list H30_n = H30/H11_0;
  distr_t_list H03_n = H03/H11_0;
  distr_t_list H11_n = H11/H11_0;
  distr_t_list H33_n = H33/H11_0;
  distr_t_list H11_0_n = H11_0/H11_0;
  distr_t_list H33_0_n = H33_0/H11_0;
  distr_t_list H11_kz0_n = (H11_kz0 -H11_0)/H11_0;
  distr_t_list H33_kz0_n = (H33_kz0 - H33_0)/H11_0;

  distr_t den = pow(Eg*kz,2)*Power(m,2)*(Eg-2.0*m);
  
  distr_t_list H30_sub = H30_n - H03_n*(m-Eg)/(2.0*m-Eg) - Include_FSES*corr_3003*exp_mk;
  distr_t_list H11_sub_temp = H11_n - H11_0_n ;
  distr_t_list H33_sub_temp = H33_n - H33_0_n*(2.0*m*Eg -ksq-resc_FSES*pow(kz,2))/(2.0*m*Eg -ksq);
  distr_t_list H33_sub = (H11_sub_temp+H33_sub_temp)/2.0;
  distr_t_list H11_sub = (H33_sub_temp-H11_sub_temp -H33_kz0_n + H11_kz0_n)/2.0;

  distr_t_list H_2 = +1.0*H30_sub/(Eg*kz*Power(m,2));
  H_2 = H_2 + H33_sub*kz*kz*(Eg-m)/den;
  H_2 = H_2+ H11_sub*(-Eg*(2.0*ksq+pow(kz,2)) + 2.0*pow(Eg,2)*m +m*(2*ksq+pow(kz,2)))/den;

 
  Get_virtual_ff_fit_interval_v1("H2", Mom.virt(), Mom.Theta(2)[2], C.Tmin, C.Tmax);
  distr_t_list contr_1 = (1.0*H30_sub/(Eg*kz*Power(m,2)));
  distr_t_list contr_2 = (H33_sub*kz*kz*(Eg-m)/den);
  distr_t_list contr_3 = (H11_sub*(-Eg*(2.0*ksq+pow(kz,2)) + 2.0*pow(Eg,2)*m +m*(2*ksq+pow(kz,2)))/den);
  distr_t_list total= H_2;
  cout<<"########## H2 mixed diag infos: "<<endl;
  cout<<"contr 30: "<<C.Fit_distr(contr_1).ave()<<" +- "<<C.Fit_distr(contr_1).err()<<endl;
  cout<<"contr 33+11: "<<C.Fit_distr(contr_2).ave()<<" +- "<<C.Fit_distr(contr_2).err()<<endl;
  cout<<"contr 33-11: "<<C.Fit_distr(contr_3).ave()<<" +- "<<C.Fit_distr(contr_3).err()<<endl;
  cout<<"total: "<<C.Fit_distr(total).ave()<<" +- "<<C.Fit_distr(total).err()<<endl;
  cout<<"########## END H2 mixed diag infos "<<endl;
  return -1.0*m*H_2*(2.0*m*Eg-ksq  );

}


distr_t_list H_1_mixed_diag(const distr_t_list& H30,const distr_t_list& H03, const distr_t_list& H11,const distr_t_list& H33, const distr_t_list& H11_0,const distr_t_list& H33_0, const distr_t_list& H33_kz0, const distr_t_list& H11_kz0, pt3_momenta& Mom,const distr_t& m) //valid for p=0, k=kz
{
  
  double Eg= Mom.Egamma();
  double kz= Mom.k()[2]; 
  double ksq = pow(Mom.virt(),2);

  auto Power = [&](const distr_t& A, int n) -> distr_t{
    distr_t ret_val = A;
    for(int i=1;i<n;i++) ret_val= ret_val*A;
    return ret_val;
  };

  
  CorrAnalysis C;
  C.Nt = 64;

  distr_t pole_corr;
  for(int ijack=0;ijack<Njacks;ijack++) pole_corr.distr.push_back( 0.5*( 1 + (m.distr[ijack] - Eg)/( 2*asinh( 0.5*sqrt(  kz*kz + pow( 2*sinh(m.distr[ijack]/2),2))))));

  distr_t corr_3003;
  for(int ijack=0;ijack < Njacks;ijack++) corr_3003.distr.push_back( -0.5*(kz*m.distr[ijack]/(2*m.distr[ijack] - Eg))*(1.0/( 2*asinh( 0.5*sqrt(  kz*kz + pow( 2*sinh(m.distr[ijack]/2),2))))));

  distr_t_list resc_FSES;
  distr_t mk_k, mk_k_no_lat1;
  for(int ijack=0;ijack<Njacks;ijack++) {
    mk_k.distr.push_back( 2*asinh( 0.5*sqrt(  kz*kz + pow( 2*sinh(m.distr[ijack]/2),2))) + Eg - m.distr[ijack]);
    mk_k_no_lat1.distr.push_back( sqrt( kz*kz + m.distr[ijack]*m.distr[ijack]));
  }
  auto exp_M= [](double mass, double t, double N) { return exp( -mass*t);};
  distr_t_list exp_mk = distr_t_list::f_of_distr(exp_M, mk_k, C.Nt);
  resc_FSES= 1-pole_corr*(1.0*Include_FSES)*exp_mk;

  distr_t_list H30_n = H30/H11_0;
  distr_t_list H03_n = H03/H11_0;
  distr_t_list H11_n = H11/H11_0;
  distr_t_list H33_n = H33/H11_0;
  distr_t_list H11_0_n = H11_0/H11_0;
  distr_t_list H33_0_n = H33_0/H11_0;
  distr_t_list H11_kz0_n = (H11_kz0-H11_0)/H11_0;
  distr_t_list H33_kz0_n = (H33_kz0-H33_0)/H11_0;
  
  distr_t den = pow(Eg*kz,2)*Power(m,2)*(Eg-2.0*m);
  
  
  distr_t_list H30_sub = H30_n - H03_n*(m-Eg)/(2.0*m-Eg) - Include_FSES*corr_3003*exp_mk;
  distr_t_list H11_sub_temp = H11_n - H11_0_n ;
  distr_t_list H33_sub_temp = H33_n - H33_0_n*(2.0*m*Eg -ksq-resc_FSES*pow(kz,2))/(2.0*m*Eg -ksq);
  distr_t_list H33_sub = (H11_sub_temp+H33_sub_temp)/2.0;
  distr_t_list H11_sub = (H33_sub_temp-H11_sub_temp-H33_kz0_n + H11_kz0_n)/2.0;


  distr_t_list H_1 = H30_sub*(ksq - Eg*m)/(Eg*kz*Power(m,2));
  H_1 = H_1 +H33_sub*kz*kz*(-2.0*pow(Eg,2)*m -ksq*m+Eg*(ksq+3.0*Power(m,2)))/den;
  H_1 = H_1+ H11_sub*(ksq*(2.0*ksq+pow(kz,2))*m -4*pow(Eg,3)*Power(m,2)+2.0*pow(Eg,2)*m*(3.0*ksq+pow(kz,2)+2.0*Power(m,2))-Eg*(2.0*pow(ksq,2)+3.0*pow(kz,2)*Power(m,2)+ksq*(pow(kz,2)+6.0*Power(m,2))))/den;


  Get_virtual_ff_fit_interval_v1("H1", Mom.virt(), Mom.Theta(2)[2], C.Tmin, C.Tmax);
  distr_t_list contr_1 = (H30_sub*(ksq - Eg*m)/(Eg*kz*Power(m,2)));
  distr_t_list contr_2 = (H33_sub*kz*kz*(-2.0*pow(Eg,2)*m -ksq*m+Eg*(ksq+3.0*Power(m,2)))/den);
  distr_t_list contr_3 = (H11_sub*(ksq*(2.0*ksq+pow(kz,2))*m -4*pow(Eg,3)*Power(m,2)+2.0*pow(Eg,2)*m*(3.0*ksq+pow(kz,2)+2.0*Power(m,2))-Eg*(2.0*pow(ksq,2)+3.0*pow(kz,2)*Power(m,2)+ksq*(pow(kz,2)+6.0*Power(m,2))))/den);
  distr_t_list total= H_1;
  cout<<"########## H1 mixed_diag infos: "<<endl;
  cout<<"contr 30: "<<C.Fit_distr(contr_1).ave()<<" +- "<<C.Fit_distr(contr_1).err()<<endl;
  cout<<"contr 33+11: "<<C.Fit_distr(contr_2).ave()<<" +- "<<C.Fit_distr(contr_2).err()<<endl;
  cout<<"contr 33-11: "<<C.Fit_distr(contr_3).ave()<<" +- "<<C.Fit_distr(contr_3).err()<<endl;
  cout<<"total: "<<C.Fit_distr(total).ave()<<" +- "<<C.Fit_distr(total).err()<<endl;
  cout<<"########## END H1 mixed_diag infos "<<endl;

  
  return H_1*m;



}


distr_t_list FA_off_mixed_diag(const distr_t_list& H30,const distr_t_list& H03, const distr_t_list& H11,const distr_t_list& H33,const distr_t_list& H11_0,const distr_t_list& H33_0, const distr_t_list& H33_kz0, const distr_t_list& H11_kz0, pt3_momenta& Mom,const distr_t &m) //valid for p=0, k=kz
{

  double Eg= Mom.Egamma();
  double kz= Mom.k()[2]; 
  double ksq = pow(Mom.virt(),2);

  auto Power = [&](const distr_t& A, int n) -> distr_t{
    distr_t ret_val = A;
    for(int i=1;i<n;i++) ret_val= ret_val*A;
    return ret_val;
  };

  CorrAnalysis C;
  C.Nt = 64;

  distr_t pole_corr;
  for(int ijack=0;ijack<Njacks;ijack++) pole_corr.distr.push_back( 0.5*( 1 + (m.distr[ijack] - Eg)/( 2*asinh( 0.5*sqrt(  kz*kz + pow( 2*sinh(m.distr[ijack]/2),2))))));

  distr_t corr_3003;
  for(int ijack=0;ijack < Njacks;ijack++) corr_3003.distr.push_back( -0.5*(kz*m.distr[ijack]/(2*m.distr[ijack] - Eg))*(1.0/( 2*asinh( 0.5*sqrt(  kz*kz + pow( 2*sinh(m.distr[ijack]/2),2))))));

  distr_t_list resc_FSES;
  distr_t mk_k;
  for(int ijack=0;ijack<Njacks;ijack++) mk_k.distr.push_back(  2*asinh( 0.5*sqrt(  kz*kz + pow( 2*sinh(m.distr[ijack]/2),2))) + Eg - m.distr[ijack]);
  auto exp_M= [](double mass, double t, double N) { return exp( -mass*t);};
  distr_t_list exp_mk = distr_t_list::f_of_distr(exp_M, mk_k, C.Nt);
  resc_FSES= 1-pole_corr*(1.0*Include_FSES)*exp_mk;

  distr_t den = pow(Eg*kz,2)*Power(m,2)*(Eg-2.0*m);

  distr_t_list H0_ave = (2.0*H11_0 + 1.0*H33_0)/3.0;

  distr_t_list H30_n = H30/H11_0;
  distr_t_list H03_n = H03/H11_0;
  distr_t_list H11_n = H11/H11_0;
  distr_t_list H33_n = H33/H11_0;
  distr_t_list H11_0_n = H11_0/H11_0;
  distr_t_list H33_0_n = H33_0/H11_0;
  distr_t_list H11_kz0_n = (H11_kz0-H11_0)/H11_0;
  distr_t_list H33_kz0_n = (H33_kz0-H33_0)/H11_0;
   
  distr_t_list H30_sub = H30_n - H03_n*(m-Eg)/(2.0*m-Eg) - Include_FSES*corr_3003*exp_mk;
  distr_t_list H11_sub_temp = H11_n - H11_0_n ;
  distr_t_list H33_sub_temp = H33_n - H33_0_n*(2.0*m*Eg -ksq-resc_FSES*pow(kz,2))/(2.0*m*Eg -ksq);
  distr_t_list H33_sub = (H11_sub_temp+H33_sub_temp)/2.0;
  distr_t_list H11_sub = (H33_sub_temp-H11_sub_temp-H33_kz0_n+H11_kz0_n)/2.0;


  distr_t_list FA_off = H30_sub*ksq/(Eg*kz*m);
  FA_off = FA_off+ H33_sub*m*kz*kz*(-pow(Eg,2)*m -ksq*m +Eg*(ksq+2.0*Power(m,2)))/den;
  FA_off = FA_off+ H11_sub*m*(ksq*(2.0*ksq+pow(kz,2))*m+ pow(Eg,2)*(4*ksq+pow(kz,2))*m -Eg*(2.0*pow(ksq,2)+2.0*pow(kz,2)*Power(m,2)+ksq*(pow(kz,2)+4.0*Power(m,2))) )/den;



  Get_virtual_ff_fit_interval_v1("FA", Mom.virt(), Mom.Theta(2)[2], C.Tmin, C.Tmax);
  distr_t_list contr_1 = (H30_sub*ksq/(Eg*kz*m));
  distr_t_list contr_2 = (H33_sub*m*kz*kz*(-pow(Eg,2)*m -ksq*m +Eg*(ksq+2.0*Power(m,2)))/den);
  distr_t_list contr_3 = (H11_sub*m*(ksq*(2.0*ksq+pow(kz,2))*m+ pow(Eg,2)*(4*ksq+pow(kz,2))*m -Eg*(2.0*pow(ksq,2)+2.0*pow(kz,2)*Power(m,2)+ksq*(pow(kz,2)+4.0*Power(m,2))) )/den);
  distr_t_list total= FA_off;
  cout<<"########## FA mixed_diag infos: "<<endl;
  cout<<"contr 30: "<<C.Fit_distr(contr_1).ave()<<" +- "<<C.Fit_distr(contr_1).err()<<endl;
  cout<<"contr 33+11: "<<C.Fit_distr(contr_2).ave()<<" +- "<<C.Fit_distr(contr_2).err()<<endl;
  cout<<"contr 33-11: "<<C.Fit_distr(contr_3).ave()<<" +- "<<C.Fit_distr(contr_3).err()<<endl;
  cout<<"total: "<<C.Fit_distr(total).ave()<<" +- "<<C.Fit_distr(total).err()<<endl;
  cout<<"########## END FA mixed_diag infos "<<endl;
 
 
  return FA_off;
}














void Compute_form_factors() {




  
  string dir= "../3pt_data/"+Meson;
  string dir_test = "../3pt_data/K_no_sm_new_seed";
  int Nens=0;
  vector<string> Ens_tag;
  

  //read number of ensembles and store the paths

  boost::filesystem::directory_iterator end_itr;

  for (boost::filesystem::directory_iterator itr(dir); itr != end_itr; ++itr) {
  

    if(!boost::filesystem::is_regular_file(itr->path())) {

      Nens++;
      Ens_tag.push_back(itr->path().string().substr(dir.length()+1));
    }
  }



  
  pt3_momenta_list mom3_l;
  pt2_momenta_list mom2_l;
  vector<vector<pt3_momenta>> mom3_l_analyzed_only;

  vector<vector<distr_t_list>> H1_list;
  vector<vector<distr_t_list>> H2_list;
  vector<vector<distr_t_list>> FA_off_list;
  vector<vector<distr_t_list>> FV_off_list;
  vector<distr_t> f_p, m_p, Za_ov_Zv, Zv_distr;
  vector<vector<distr_t>> H1_const_fit_list;
  vector<vector<distr_t>> H2_const_fit_list;
  vector<vector<distr_t>> FA_off_const_fit_list;
  vector<vector<distr_t>> FV_off_const_fit_list;
  vector<distr_t> ainv;

  
  //init RNG
  GaussianMersenne Gauss_RC(91068231);
  
  
  
  

  for (int i_ens=0; i_ens<Nens;i_ens++) {


    FILE *stream_3pt;
    FILE *stream_2pt;

    FILE *stream_2pt_test;

    struct header_virph header_3pt;
    struct header_virph header_2pt;

    struct header_virph header_2pt_test;
   
    

   

    string Path3pt = dir+"/"+Ens_tag[i_ens]+"/conf.virtualph.dat";
    string Path2pt = dir+"/"+Ens_tag[i_ens]+"/conf.virtualph.dat2";


    string Path2pt_test = dir_test+"/"+Ens_tag[i_ens]+"/conf.virtualph.dat2";



    stream_3pt= fopen(Path3pt.c_str(), "rb");
    stream_2pt= fopen(Path2pt.c_str(), "rb");

    stream_2pt_test= fopen(Path2pt_test.c_str(), "rb");



    read_header_bin(stream_3pt, header_3pt);
    read_header_bin(stream_2pt, header_2pt);

    read_header_bin(stream_2pt_test, header_2pt_test);


  
    //Forward lattice  informations to 3pt_momenta_list
    double ml, l, Nt;
    Read_pars_from_ensemble_tag(Ens_tag[i_ens], ml, l , Nt);
    Add_to_mom_list(mom3_l, header_3pt, l);
    Add_to_mom_list(mom2_l, header_2pt, l);
    mom2_l.Nt.push_back(Nt);
    mom3_l.Nt.push_back(Nt);


    //print size of header file
    cout<<"Ens: "<<Ens_tag[i_ens]<<", header size 2pt (Bytes): "<<header_2pt.header_size<<endl;
    cout<<"Number of configs (from 2pt): "<<Get_number_of_configs_2pt(stream_2pt, header_2pt)<<endl;
    cout<<"Ens: "<<Ens_tag[i_ens]<<", header size 3pt (Bytes): "<<header_3pt.header_size<<endl;
    cout<<"Number of configs (from 3pt): "<<Get_number_of_configs_3pt(stream_3pt, header_3pt)<<endl;
    cout<<"Number of smearing levels present: "<<header_3pt.nqsml<<endl;


    boost::filesystem::create_directory("../data");
    boost::filesystem::create_directory("../data/form_factors");
    boost::filesystem::create_directory("../data/form_factors/"+Meson);
    boost::filesystem::create_directory("../data/form_factors/"+Meson+"/H_"+Ens_tag[i_ens]);
    boost::filesystem::create_directory("../data/form_factors/"+Meson+"/C_"+Ens_tag[i_ens]);
    boost::filesystem::create_directory("../data/form_factors/"+Meson+"/FORM_FACTOR_LIST");
    boost::filesystem::create_directory("../data/form_factors/"+Meson+"/VIRTUAL_FORM_FACTOR_LIST");
    boost::filesystem::create_directory("../data/form_factors/"+Meson+"/Contaminations");
    
    //PRINT HEADERS
    //#########################################################
    ofstream P_head;
    if(Determine_contaminations) {
      P_head.open("../data/form_factors/"+Meson+"/Contaminations/contaminations_RV_"+Ens_tag[i_ens]+".dat");
      P_head<<"#xg"<<setw(20)<<"offsh"<<setw(20)<<"sm_lev"<<setw(20)<<"a1"<<setw(20)<<"a1_err"<<setw(20)<<"a2"<<setw(20)<<"a2_err"<<setw(20)<<"a3"<<setw(20)<<"a3_err"<<setw(20)<<"a5"<<setw(20)<<"a5_err"<<setw(20)<<"ch2/ndof"<<endl;
      P_head.close();
      P_head.open("../data/form_factors/"+Meson+"/Contaminations/contaminations_RA_"+Ens_tag[i_ens]+".dat");
      P_head<<"#xg"<<setw(20)<<"offsh"<<setw(20)<<"sm_lev"<<setw(20)<<"a1"<<setw(20)<<"a1_err"<<setw(20)<<"a2"<<setw(20)<<"a2_err"<<setw(20)<<"a3"<<setw(20)<<"a3_err"<<setw(20)<<"a5"<<setw(20)<<"a5_err"<<setw(20)<<"ch2/ndof"<<endl;
      P_head.close();
      P_head.open("../data/form_factors/"+Meson+"/Contaminations/contaminations_FA_"+Ens_tag[i_ens]+".dat");
      P_head<<"#xg"<<setw(20)<<"offsh"<<setw(20)<<"sm_lev"<<setw(20)<<"a1"<<setw(20)<<"a1_err"<<setw(20)<<"a2"<<setw(20)<<"a2_err"<<setw(20)<<"a3"<<setw(20)<<"a3_err"<<setw(20)<<"a5"<<setw(20)<<"a5_err"<<setw(20)<<"ch2/ndof"<<endl;
      P_head.close();
      P_head.open("../data/form_factors/"+Meson+"/Contaminations/contaminations_FV_"+Ens_tag[i_ens]+".dat");
      P_head<<"#xg"<<setw(20)<<"offsh"<<setw(20)<<"sm_lev"<<setw(20)<<"a1"<<setw(20)<<"a1_err"<<setw(20)<<"a2"<<setw(20)<<"a2_err"<<setw(20)<<"a3"<<setw(20)<<"a3_err"<<setw(20)<<"a5"<<setw(20)<<"a5_err"<<setw(20)<<"ch2/ndof"<<endl;
      P_head.close();
    }
    P_head.open("../data/form_factors/"+Meson+"/FORM_FACTOR_LIST/vector_form_factors_list_"+Ens_tag[i_ens]+".dat");
    P_head<<"#xg"<<setw(20)<<"offsh"<<setw(20)<<"sm_lev"<<setw(20)<<"F"<<setw(20)<<"F_err"<<endl;
    P_head.close();
    P_head.open("../data/form_factors/"+Meson+"/FORM_FACTOR_LIST/axial_form_factors_list_"+Ens_tag[i_ens]+".dat");
    P_head<<"#xg"<<setw(20)<<"offsh"<<setw(20)<<"sm_lev"<<setw(20)<<"F"<<setw(20)<<"F_err"<<endl;
    P_head.close();
    P_head.open("../data/form_factors/"+Meson+"/FORM_FACTOR_LIST/FF_optimized_FV_"+Ens_tag[i_ens]+".dat");
    P_head<<"#xg"<<setw(20)<<"offsh"<<setw(20)<<"sm_lev"<<setw(20)<<"F"<<setw(20)<<"F_err"<<setw(20)<<"T_min"<<setw(20)<<"T_max"<<endl;
    P_head.close();
    P_head.open("../data/form_factors/"+Meson+"/FORM_FACTOR_LIST/FF_optimized_FA_"+Ens_tag[i_ens]+".dat");
    P_head<<"#xg"<<setw(20)<<"offsh"<<setw(20)<<"sm_lev"<<setw(20)<<"F"<<setw(20)<<"F_err"<<setw(20)<<"T_min"<<setw(20)<<"T_max"<<endl;
    P_head.close();
    P_head.open("../data/form_factors/"+Meson+"/FORM_FACTOR_LIST/FF_optimized_RV_"+Ens_tag[i_ens]+".dat");
    P_head<<"#xg"<<setw(20)<<"offsh"<<setw(20)<<"sm_lev"<<setw(20)<<"F"<<setw(20)<<"F_err"<<setw(20)<<"T_min"<<setw(20)<<"T_max"<<endl;
    P_head.close();
    P_head.open("../data/form_factors/"+Meson+"/FORM_FACTOR_LIST/FF_optimized_RA_"+Ens_tag[i_ens]+".dat");
    P_head<<"#xg"<<setw(20)<<"offsh"<<setw(20)<<"sm_lev"<<setw(20)<<"F"<<setw(20)<<"F_err"<<setw(20)<<"T_min"<<setw(20)<<"T_max"<<endl;
    P_head.close();
    P_head.open("../data/form_factors/"+Meson+"/VIRTUAL_FORM_FACTOR_LIST/"+Ens_tag[i_ens]+".dat");
    P_head<<"#xk"<<setw(20)<<"xq"<<setw(20)<<"H1"<<setw(20)<<"H1_err"<<setw(20)<<"H2"<<setw(20)<<"H2_err"<<setw(20)<<"FA"<<setw(20)<<"FA_err"<<setw(20)<<"FV"<<setw(20)<<"FV_err"<<endl;
    P_head.close();
    P_head.open("../data/form_factors/"+Meson+"/VIRTUAL_FORM_FACTOR_LIST/"+Ens_tag[i_ens]+"_impr.dat");
    P_head<<"#xk"<<setw(20)<<"xq"<<setw(20)<<"H1"<<setw(20)<<"H1_err"<<setw(20)<<"H2"<<setw(20)<<"H2_err"<<setw(20)<<"FA"<<setw(20)<<"FA_err"<<setw(20)<<"FV"<<setw(20)<<"FV_err"<<endl;
    P_head.close();
    P_head.open("../data/form_factors/"+Meson+"/VIRTUAL_FORM_FACTOR_LIST/"+Ens_tag[i_ens]+"_mixed_diag.dat");
    P_head<<"#xk"<<setw(20)<<"xq"<<setw(20)<<"H1"<<setw(20)<<"H1_err"<<setw(20)<<"H2"<<setw(20)<<"H2_err"<<setw(20)<<"FA"<<setw(20)<<"FA_err"<<setw(20)<<"FV"<<setw(20)<<"FV_err"<<endl;
    P_head.close();
        

    //##########################################################


    
    for(int nqsm=0; nqsm<header_3pt.nqsml;nqsm++) {

      int Analyze_sm_lev=1;
      if(header_3pt.nqsml == 1) sm_lev=0;
      else {
	cout<<"Should I analyze smearing level: "<<nqsm<<" ?"<<endl;
	cout<<"1->yes, 0->no"<<endl;
	cin>>Analyze_sm_lev;
      }

      
      if(Analyze_sm_lev) {
	sm_lev=nqsm;

	//resize FF_LSIT
	H1_list.resize(H1_list.size()+1);
	H2_list.resize(H2_list.size()+1);
	FA_off_list.resize(FA_off_list.size()+1);
	FV_off_list.resize(FV_off_list.size()+1);
	H1_const_fit_list.resize(H1_const_fit_list.size()+1);
	H2_const_fit_list.resize(H2_const_fit_list.size()+1);
	FA_off_const_fit_list.resize(FA_off_const_fit_list.size()+1);
	FV_off_const_fit_list.resize(FV_off_const_fit_list.size()+1);
	
	mom3_l_analyzed_only.resize(mom3_l_analyzed_only.size()+1);
	
	//STORAGE FOR 2-pt FUNC OBS
	vector<distr_t_list> fp_distr_list, m_distr_list;
	vector<distr_t> fp_fit_distr, m_fit_distr, sqrt_overlap_fit_distr;
	vector<distr_t_list> sqrt_overlap_distr, m_distr;


	 
	
	
	
	



	//init CorrAnalysis class
	//CorrAnalysis corr(UseJack, Get_number_of_configs_3pt(stream_3pt, header_3pt), nboots);
	CorrAnalysis corr(UseJack, Njacks, nboots); //in case you want to use a different number of NJacks
	corr.Nt=Nt;


	//get lattice spacing
	//###############################################################################
	LatticeInfo L_info("LOCAL"); //CURRENT_TYPE == LOCAL
	L_info.LatInfo("A");
	distr_t a_temp(UseJack);
	int jackmax= (UseJack)?corr.Njacks:corr.Nboots;
	for(int ijack=0;ijack<jackmax;ijack++) {
	  if(UseJack) a_temp.distr.push_back( L_info.ainv + L_info.ainv_err*(1.0/(sqrt(corr.Njacks-1)))*Gauss_RC());
	  else a_temp.distr.push_back( L_info.ainv + L_info.ainv_err*Gauss_RC());
	}
	ainv.push_back(a_temp);
	//###############################################################################


	
	//DEFINE LAMBDA TO REMOVE EXPONENTIAL DEPENDENCE ON MESON ENERGY FROM THE CORRELATORS
	//###############################################################################
	auto Exp_lat = [&] (double a, double b, double c) { return (b<c/2)?1.0/(exp(-a*b)):1.0/(exp(-a*(c-b)));};
	auto Exp_lat2 = [&] (double a, int t) { return (t<Nt/2)?exp(a*t):exp(a*(Nt-t));};
	//#############################################################################

    

	//2pts
	Get_Tmin_Tmax("2pt", Ens_tag[i_ens], corr, 0, "V");
	//#########################################################
	

	int ncomb2pt= header_2pt.ncomb;
	int ncomb3pt= header_3pt.ncomb;

	cout<<"#########################"<<endl;
	cout<<"number of combs in 2pt: "<<ncomb2pt<<endl;
	cout<<"number of combs in 3pt: "<<ncomb3pt<<endl;
	cout<<"#########################"<<endl;



	//THIS PART OF THE CODE HANDLES THE TWO-PT FUNCTION. IT COMPUTES EFFECTIVE MASSES AND MESON'S DECAY CONSTANT, FOR ALL KINEMATICS.
	//#############################################################################
	//#############################################################################
	//#############################################################################

	for(int icomb2pt=0; icomb2pt<ncomb2pt;icomb2pt++) {

	  corr.Perform_Nt_t_average=1;
	  corr.Reflection_sign=1;
	  
	  auto c = header_2pt.comb[icomb2pt];	  
	  string tag=to_string(icomb2pt);

	  distr_t_list corr_2pt_distr = corr.corr_t(Get_obs_2pt(stream_2pt,header_2pt,0,icomb2pt,0,sm_lev), "../data/form_factors/"+Meson+"/corr_2pt_"+Ens_tag[i_ens]+"_sm_lev"+to_string(sm_lev)+"_"+tag);

	  distr_t_list corr_2pt_distr_test = corr.corr_t(Get_obs_2pt(stream_2pt_test,header_2pt_test,0,icomb2pt,0,sm_lev), "");

	  //Print difference
	  Print_To_File({}, { (corr_2pt_distr-corr_2pt_distr_test).ave(), (corr_2pt_distr-corr_2pt_distr_test).err()}, "../data/form_factors/"+Meson+"/corr_diff_2pt_"+Ens_tag[i_ens]+"_sm_lev"+to_string(sm_lev)+"_"+tag, "","");
	  
	  distr_t_list overlap_distr_list = corr.residue_t(Get_obs_2pt(stream_2pt,header_2pt,0,icomb2pt,0,0), "");

	  distr_t_list overlap_smeared_distr_list = corr.residue_t(Get_obs_2pt(stream_2pt,header_2pt,0,icomb2pt,0,sm_lev), "../data/form_factors/"+Meson+"/overlap_"+Ens_tag[i_ens]+"_sm_lev"+to_string(sm_lev)+"_"+tag);

	  distr_t overlap_fit_distr =corr.Fit_distr(overlap_distr_list);
	  distr_t overlap_smeared_fit_distr=corr.Fit_distr( overlap_smeared_distr_list);

	  fp_distr_list.push_back( (c.mu1+c.mu2)*corr.decay_constant_t( Get_obs_2pt(stream_2pt, header_2pt, 0, icomb2pt, 0, 0), "../data/form_factors/"+Meson+"/fp_"+Ens_tag[i_ens]+"_sm_lev"+to_string(sm_lev)+"_"+tag));

	  m_distr_list.push_back( corr.effective_mass_t(Get_obs_2pt(stream_2pt, header_2pt, 0, icomb2pt, 0, sm_lev), "../data/form_factors/"+Meson+"/meson_mass_"+Ens_tag[i_ens]+"_sm_lev"+to_string(sm_lev)+"_"+tag));

	  fp_fit_distr.push_back(corr.Fit_distr(fp_distr_list[icomb2pt]));

	  m_fit_distr.push_back(corr.Fit_distr(m_distr_list[icomb2pt]));

	  
	  //print 2pt info
	  cout<<"#####################"<<endl;
	  cout<<"icomb_2pt: "<<icomb2pt<<endl;
	  cout<<"m1: "<<c.mu1<<endl;
	  cout<<"m2: "<<c.mu2<<endl;
	  cout<<"off: "<<c.off<<endl;
	  cout<<"th_0: "<<c.th0[0]<<", "<<c.th0[1]<<", "<<c.th0[2]<<endl;
	  cout<<"th_s: "<<c.ths[0]<<", "<<c.ths[1]<<", "<<c.ths[2]<<endl;
	  cout<<"th_t: "<<c.tht[0]<<", "<<c.tht[1]<<", "<<c.tht[2]<<endl;
	  
	  cout<<"#####################"<<endl;
	  
	  if(icomb2pt==0) {
	    f_p.push_back(fp_fit_distr[0]); m_p.push_back(m_fit_distr[0]);
	    //read RC
	   
	    double Zv, Zv_err, Za, Za_err;
	    Zv=L_info.Retrieve_Zv("A",0).first;
	    Zv_err= L_info.Retrieve_Zv("A",0).second;
	    Za=L_info.Retrieve_Za("A",0).first;
	    Za_err= L_info.Retrieve_Za("A",0).second;
	    //define Za_ov_Zv distr
	    distr_t Za_ov_Zv_distr(UseJack);
	    distr_t Zv_fake_distr(UseJack);
	    
	    if(UseJack) {
	      for(int ijack=0;ijack<m_fit_distr[0].size();ijack++) {
		Za_ov_Zv_distr.distr.push_back( (Za + (1.0/sqrt( m_fit_distr[0].size() -1))*Za_err*Gauss_RC())/(Zv+(1.0/sqrt(m_fit_distr[0].size()-1))*Zv_err*Gauss_RC()));
		//Za_ov_Zv_distr.distr.push_back( 1 + 0.0001*Gauss_RC());
		Zv_fake_distr.distr.push_back( (Zv + Zv_err*Gauss_RC()/(sqrt(m_fit_distr[0].size() -1))));
	      }
	    }
	    else { //using bootstrap
	      for(int iboot=0;iboot<m_fit_distr[0].size();iboot++) {
			Za_ov_Zv_distr.distr.push_back( (Za + Za_err*Gauss_RC())/(Zv+Zv_err*Gauss_RC()));
			Zv_fake_distr.distr.push_back( (Zv + Zv_err*Gauss_RC()));
	      }
	    }
	    
	    Za_ov_Zv.push_back(Za_ov_Zv_distr);
	    Zv_distr.push_back(Zv_fake_distr);
	  }
	  
	  auto sq= [](double x)->double {return sqrt(x);};
	  auto sq_t= [](double x, int t) -> double { return sqrt(x);};
	  sqrt_overlap_distr.push_back( overlap_smeared_distr_list/(distr_t_list::f_of_distr_list(sq_t, overlap_distr_list)));
	  sqrt_overlap_fit_distr.push_back( overlap_smeared_fit_distr/(distr_t::f_of_distr(sq,overlap_fit_distr)));
	}

	//#############################################################################
	//#############################################################################
	//#############################################################################
	//end anaysis 2pts

    
	//BEGINNING ANALYSIS OF 3PT-FUNCTIONS
	//#############################################################################
	//#############################################################################
	//#############################################################################
      
	//loop over kinematics
	for(int icomb3pt=0; icomb3pt<ncomb3pt;icomb3pt++) {

	  auto c = header_3pt.comb[icomb3pt];

	  //if(icomb3pt==0) corr.Perform_Nt_t_average=1;
	  //else corr.Perform_Nt_t_average=0;
	  
	  vector<vector<distr_t_list>> distr_V(4), distr_V_k0(4), distr_A(4), distr_A_k0(4), distr_A_exp(4), distr_A_exp_k0(4), distr_V_exp(4);
	  vector<vector<distr_t_list>> distr_no_symm_V(4), distr_no_symm_V_k0(4), distr_no_symm_A(4), distr_no_symm_A_k0(4);
	  vector<vector<distr_t_list>> distr_no_symm_A_kz0_k2(4), distr_A_kz0_k2(4);
	  for(int i=0;i<4;i++) { distr_V[i].resize(4); distr_V_k0[i].resize(4); distr_A[i].resize(4); distr_A_k0[i].resize(4);distr_A_exp[i].resize(4); distr_V_exp[i].resize(4); distr_A_exp_k0[i].resize(4);}
	  for(int i=0;i<4;i++) { distr_no_symm_V[i].resize(4); distr_no_symm_V_k0[i].resize(4); distr_no_symm_A[i].resize(4); distr_no_symm_A_k0[i].resize(4); distr_no_symm_A_kz0_k2[i].resize(4); distr_A_kz0_k2[i].resize(4); }
	  
	 

	  //infos about the kinematic
	  //####################################################

      
	  int icomb_k0= 0;
	  if(FIND_k0_KIN) icomb_k0= Get_comb_k0(header_3pt, icomb3pt);
	  int symmetric_comb=Get_symmetric_comb(header_3pt, icomb3pt);
	  int symmetric_comb_k0 = 0;
	  if(FIND_k0_KIN) symmetric_comb_k0 = Get_symmetric_comb(header_3pt, icomb_k0);
	  int icomb_kz0_k2=0;
	  int symmetric_comb_kz0_k2=0;
	  if(VIRTUAL_RUN && FIND_k0_kz0_KIN) {
	  icomb_kz0_k2 = Get_comb_k0_same_off(header_3pt, icomb3pt);
	  symmetric_comb_kz0_k2 = Get_symmetric_comb(header_3pt, icomb_kz0_k2);
	  }
	  else { icomb_kz0_k2 = -1; symmetric_comb_kz0_k2=-1;}
	  
	  int pt2_k0p0 = 0;
	  if(FIND_k0_KIN) Get_2pt_k0p0(header_2pt, c.mu1, c.mu2);
	  int pt2_p = Get_2pt_p(header_2pt, c.i0, c.is);
	  double Egamma = mom3_l.mom[i_ens][icomb3pt].Egamma();
	  double Egamma_ash = mom3_l.mom[i_ens][icomb3pt].Egamma_asinh();
	  double xg = mom3_l.mom[i_ens][icomb3pt].x_gamma(m_fit_distr[pt2_k0p0]).ave();
	  double xg_off =  mom3_l.mom[i_ens][icomb3pt].x_gamma_off(m_fit_distr[pt2_k0p0]).ave();
	  double offsh= mom3_l.mom[i_ens][icomb3pt].virt();
	  double Egamma_T= sinh(Egamma)*(1.0-exp(-Egamma*Nt));
	  double meson_mass = m_fit_distr[pt2_k0p0].ave();
	  double xk= sqrt(offsh*offsh/pow(meson_mass,2));
	  double xq= sqrt((1+ offsh*offsh/pow(meson_mass,2) -xg));
	  
	 
	  	 
     
     
	  //PRINT INFOS
	  cout<<"###BEG###"<<endl;
	  cout<<"icomb: "<<icomb3pt<<" th_t: "<<c.tht[2]<<endl;
	  cout<<"nconfs: "<<Get_number_of_configs_3pt(stream_3pt, header_3pt)<<endl;
	  cout<<"c.mu1: "<<c.mu1<<endl;
	  cout<<"c.mu2: "<<c.mu2<<endl;
	  if(FIND_k0_KIN) {
	  cout<<"icomb_k0: "<<icomb_k0<<" th_t: "<<header_3pt.comb[icomb_k0].tht[2]<<endl;
	  }
	  if(VIRTUAL_RUN && FIND_k0_kz0_KIN) {
	  cout<<"icomb_kz0_k2: "<<icomb_kz0_k2<<" th_t: "<<header_3pt.comb[icomb_kz0_k2].tht[2]<<", off: "<<header_3pt.comb[icomb_kz0_k2].off<<endl;
	  }
	  cout<<"icomb_symm: "<<symmetric_comb<<" th_t: "<<header_3pt.comb[symmetric_comb].tht[2]<<endl;
	  if(FIND_k0_KIN) {
	  cout<<"icomb_symm_k0: "<<symmetric_comb_k0<<" th_t: "<<header_3pt.comb[symmetric_comb_k0].tht[2]<<endl;
	  }
	  if(VIRTUAL_RUN) {
	  cout<<"icomb_symm_kz0_k2: "<<symmetric_comb_kz0_k2<<" th_t: "<<header_3pt.comb[symmetric_comb_kz0_k2].tht[2]<<", off= "<<header_3pt.comb[symmetric_comb_kz0_k2].off<<endl;
	  }
	  cout<<"fp: "<<fp_fit_distr[pt2_k0p0].ave()<<"("<<fp_fit_distr[pt2_k0p0].err()<<")"<<endl;
	  cout<<"k_z: "<<mom3_l.mom[i_ens][icomb3pt].k()[2]<<endl;
	  cout<<"Egamma["<<icomb3pt<<"] :"<<Egamma<<endl;
	  cout<<"Egamma asinh["<<icomb3pt<<"] :"<<Egamma_ash<<endl;
	  cout<<"EgammaT["<<icomb3pt<<"] :"<<Egamma_T<<endl;
	  cout<<"Mass: "<<meson_mass<<"("<<m_fit_distr[pt2_k0p0].err()<<")"<<endl;
	  cout<<"fp^3pt*Mp: "<<(fp_fit_distr[pt2_k0p0]*(1.0/Zv_distr[0])*m_fit_distr[pt2_k0p0]).ave()<<" +- "<<(fp_fit_distr[pt2_k0p0]*(1.0/Zv_distr[0])*m_fit_distr[pt2_k0p0]).err()<<endl;
	  cout<<"x_gamma: "<<xg<<endl;
	  cout<<"x_gamma_off: "<<xg_off<<endl;
	  cout<<"xk: "<<xk<<endl;
	  cout<<"xq: "<<xq<<endl;
	  cout<<"kz/mk: "<<mom3_l.mom[i_ens][icomb3pt].k()[2]/meson_mass<<endl;
	  cout<<"Symmetrize 3pt: "<<corr.Perform_Nt_t_average<<endl;
	  cout<<"Sub kz=0 with same off: "<<SUB_ZERO_MOMENTUM_VIRTUAL_INT<<endl;
	  
	  cout<<"##############"<<endl;
	


	  //####################################################

	  //TO REMOVE EXPONENTIAL IN THE PHOTON ENERGY
	  //####################################################
	  Vfloat exp_Nt_t;
	  for(int t=0; t<Nt;t++) exp_Nt_t.push_back( exp( abs((Nt/2-t))*Egamma));

	  Vfloat exp_Nt_t_k0;
	  for(int t=0; t<Nt;t++) exp_Nt_t_k0.push_back( exp(   abs((Nt/2-t))*offsh));
	  
	  
	  //####################################################


	  //loop over alpha and mu
	  for(int alpha=0; alpha<4;alpha++) {
	    for(int mu=0; mu<4;mu++) {

	  
	      // compute_vector_correlators C_V_alpha^mu for each alpha and mu
	      //######################################################################

	      corr.Reflection_sign = -1;
   
	      distr_no_symm_V[alpha][mu] =e_f1*corr.corr_t(Get_obs_3pt(stream_3pt, header_3pt, 1, icomb3pt, alpha, mu, "V", sm_lev), "");
	      if(FIND_k0_KIN) {
	      distr_no_symm_V_k0[alpha][mu] =e_f1*corr.corr_t(Get_obs_3pt(stream_3pt, header_3pt, 1,icomb_k0, alpha, mu, "V", sm_lev), "");
	      }
	      distr_V[alpha][mu] = distr_no_symm_V[alpha][mu] + e_f2*corr.corr_t(Get_obs_3pt(stream_3pt, header_3pt, 1, symmetric_comb, alpha, mu, "V", sm_lev), "");
	      if(FIND_k0_KIN) {
	      distr_V_k0[alpha][mu] = distr_no_symm_V_k0[alpha][mu]  + e_f2*corr.corr_t(Get_obs_3pt(stream_3pt, header_3pt, 1, symmetric_comb_k0, alpha, mu, "V", sm_lev), "");
	      }

	      //######################################################################


	      //compute_axial_correlators C_A_alpha^\mu for each alpha and mu
	      //######################################################################
	      
	      corr.Reflection_sign = 1;
	  
	      int Im_Re= 0; double parity=1.0; double sign=1;
	      if( (alpha==0 || mu==0) && (alpha != 0 || mu != 0)  ) {Im_Re=1; corr.Reflection_sign=-1;sign=1; parity=-1;}
	      distr_no_symm_A[alpha][mu] = parity*e_f1*corr.corr_t(Get_obs_3pt(stream_3pt, header_3pt, Im_Re, icomb3pt, alpha, mu, "A", sm_lev), "");
	      if(FIND_k0_KIN) {
	      distr_no_symm_A_k0[alpha][mu] =parity*e_f1*corr.corr_t(Get_obs_3pt(stream_3pt, header_3pt,Im_Re, icomb_k0, alpha, mu, "A", sm_lev), "");
	      }
	      distr_A[alpha][mu] = distr_no_symm_A[alpha][mu]  - parity*sign*e_f2*corr.corr_t(Get_obs_3pt(stream_3pt, header_3pt, Im_Re, symmetric_comb, alpha, mu, "A", sm_lev), "");
	      if(FIND_k0_KIN) {
	      distr_A_k0[alpha][mu] = distr_no_symm_A_k0[alpha][mu]   - parity*sign*e_f2*corr.corr_t(Get_obs_3pt(stream_3pt, header_3pt, Im_Re, symmetric_comb_k0, alpha, mu, "A", sm_lev), "");
	      }


	      
	      //compute correlators for the kinematics with same k2 but kz=0. Do it only if correctly found and SUBTRACT_ZERO_MOMENTUM_VIRTUAL is nonzero.
	      if(SUB_ZERO_MOMENTUM_VIRTUAL_INT && (icomb_kz0_k2 >= 0) && (symmetric_comb_kz0_k2 >=0)) { 
	      distr_no_symm_A_kz0_k2[alpha][mu] = parity*e_f1*corr.corr_t(Get_obs_3pt(stream_3pt, header_3pt,Im_Re, icomb_kz0_k2, alpha, mu, "A", sm_lev), "");
	      distr_A_kz0_k2[alpha][mu] = distr_no_symm_A_kz0_k2[alpha][mu]   - parity*sign*e_f2*corr.corr_t(Get_obs_3pt(stream_3pt, header_3pt, Im_Re, symmetric_comb_kz0_k2, alpha, mu, "A", sm_lev), "");
	      }
	      else {
		distr_no_symm_A_kz0_k2[alpha][mu] = 0.0*distr_A_k0[alpha][mu];
		distr_A_kz0_k2[alpha][mu] = 0.0*distr_no_symm_A_kz0_k2[alpha][mu];
	      }
	      //#########################################################################################################################################
	    
	      
	      //######################################################################
	      corr.Reflection_sign= 1;
	      sign=1;
	      parity=1.0;
	      

	      //PRINT TO FILE THE TENSOR C_W_alpha^\mu
	      Print_To_File({}, {distr_no_symm_A[alpha][mu].ave(), distr_no_symm_A[alpha][mu].err(), distr_no_symm_V[alpha][mu].ave(), distr_no_symm_V[alpha][mu].err()}, "../data/form_factors/"+Meson+"/C_"+Ens_tag[i_ens]+"/no_symm_icomb"+to_string(icomb3pt)+"_alpha"+to_string(alpha)+"_mu"+to_string(mu)+"_sm_lev_"+to_string(sm_lev)+".dat","", "#");

	      Print_To_File({}, {distr_A[alpha][mu].ave(), distr_A[alpha][mu].err(), distr_V[alpha][mu].ave(), distr_V[alpha][mu].err()}, "../data/form_factors/"+Meson+"/C_"+Ens_tag[i_ens]+"/icomb_"+to_string(icomb3pt)+"_alpha"+to_string(alpha)+"_mu"+to_string(mu)+"_sm_lev_"+to_string(sm_lev)+".dat","", "#");


	      
	      //COMPUTE THE TENSOR H_W_alpha^\mu
	      //########################################################################


	      distr_t_list distr_no_symm_A_exp= (m_fit_distr[pt2_p]/sqrt_overlap_fit_distr[pt2_p])*distr_t_list::f_of_distr(Exp_lat, m_fit_distr[pt2_p],Nt)*(distr_no_symm_A[alpha][mu]*exp_Nt_t); 
	      distr_t_list distr_no_symm_V_exp=  (m_fit_distr[pt2_p]/sqrt_overlap_fit_distr[pt2_p])*distr_t_list::f_of_distr(Exp_lat, m_fit_distr[pt2_p],Nt)*(distr_no_symm_V[alpha][mu]*exp_Nt_t);
	      distr_A_exp[alpha][mu] = (2.0*m_fit_distr[pt2_p]/sqrt_overlap_fit_distr[pt2_p])*distr_t_list::f_of_distr(Exp_lat, m_fit_distr[pt2_p],Nt)*(distr_A[alpha][mu]*exp_Nt_t);
	      if(FIND_k0_KIN) {
	      distr_A_exp_k0[alpha][mu] = (1.0*m_fit_distr[pt2_k0p0]/sqrt_overlap_fit_distr[pt2_k0p0])*distr_t_list::f_of_distr(Exp_lat, m_fit_distr[pt2_k0p0],Nt)*(distr_A_k0[alpha][mu]);
	      }
	      distr_V_exp[alpha][mu] =  (2.0*m_fit_distr[pt2_p]/sqrt_overlap_fit_distr[pt2_p])*distr_t_list::f_of_distr(Exp_lat, m_fit_distr[pt2_p],Nt)*(distr_V[alpha][mu]*exp_Nt_t);

	      //PRINT TO FILE THE TENSOR H_alpha^mu
	      Print_To_File({}, {distr_no_symm_A_exp.ave(), distr_no_symm_A_exp.err(), distr_no_symm_V_exp.ave(), distr_no_symm_V_exp.err()}, "../data/form_factors/"+Meson+"/H_"+Ens_tag[i_ens]+"/no_symm_icomb"+to_string(icomb3pt)+"_alpha"+to_string(alpha)+"_mu"+to_string(mu)+"_sm_lev_"+to_string(sm_lev)+".dat","", "#");
	      Print_To_File({}, {distr_A_exp[alpha][mu].ave(), distr_A_exp[alpha][mu].err(), distr_V_exp[alpha][mu].ave(), distr_V_exp[alpha][mu].err()}, "../data/form_factors/"+Meson+"/H_"+Ens_tag[i_ens]+"/icomb_"+to_string(icomb3pt)+"_alpha"+to_string(alpha)+"_mu"+to_string(mu)+"_sm_lev_"+to_string(sm_lev)+".dat","", "#");


	      //########################################################################

	 
	    }
	  }

	  exit(-1);
	 
	 	      
	  //COMPUTE R and bar_R and form factors
	  distr_t_list bar_R_A_distr, bar_R_V_distr, R_A_distr, R_V_distr, bar_R_A_distr_exp;

	  if( (Meson.substr(0,1) == "D" && c.mu1 > c.mu2  ) || (Meson.substr(0,1) =="K" && c.mu2 > c.mu1   )) {
      
	    //####################################################################################################################################
	    if(offsh == 0) {
	      distr_t_list V_ave = V_ave_unpolarized(distr_V,distr_V_k0, m_fit_distr[pt2_p], mom3_l.mom[i_ens][icomb3pt],Nt/2);
	      distr_t_list A_ave = A_ave_unpolarized(distr_A);
	      distr_t_list A_ave_exp = A_ave_unpolarized(distr_A_exp);
	      distr_t_list A_ave_zero_mom = A_ave_unpolarized(distr_A_k0);
	      distr_t_list A_ave_exp_zero_mom= A_ave_unpolarized(distr_A_exp_k0);
	      distr_t F_A_distr_factor= (2.0*fp_fit_distr[pt2_k0p0]/(m_fit_distr[pt2_k0p0]*mom3_l.mom[i_ens][icomb3pt].x_gamma(m_fit_distr[pt2_k0p0])));
	      auto Exp_photon= [&](double a, int t) { return exp(fabs(corr.Nt/2-t)*a);};

	      // distr_t_list eff_photon_energy= corr.effective_mass_t( A_ave_zero_mom/A_ave,"");
	      
	      //bar_R_A_distr=(distr_t_list::f_of_distr_list(Exp_photon, eff_photon_energy)*A_ave/A_ave_zero_mom);
	      distr_t fp_from_3pt = corr.Fit_distr(A_ave_exp_zero_mom);
	      bar_R_A_distr=(exp_Nt_t*A_ave/A_ave_zero_mom -1.0);
	      bar_R_A_distr_exp= ((A_ave_exp -fp_from_3pt)/fp_from_3pt);
	    	      
	      bar_R_V_distr= (V_ave/A_ave_zero_mom)*m_fit_distr[pt2_k0p0]*fp_fit_distr[pt2_k0p0];
	      R_A_distr= (1.0*m_fit_distr[pt2_p]/(2.0*m_fit_distr[pt2_k0p0]))*(A_ave/sqrt_overlap_fit_distr[icomb3pt])*distr_t_list::f_of_distr(Exp_lat, m_fit_distr[pt2_p], Nt)*exp_Nt_t;
	      R_V_distr= (m_fit_distr[pt2_p]*m_fit_distr[pt2_k0p0]/4.0)*(V_ave/sqrt_overlap_fit_distr[icomb3pt])*2.0*distr_t_list::f_of_distr(Exp_lat, m_fit_distr[pt2_p], Nt);

	      distr_t_list F_A_distr = bar_R_A_distr*F_A_distr_factor;
	      distr_t_list F_A_distr_exp = bar_R_A_distr_exp*F_A_distr_factor;
	      
	      //printV(A_ave_exp_zero_mom.ave(), "fp_3pt",1);
	      
	      //#####################################################################################################################################
	      
	      //FIT form factors and R estimators
	      // #########################################
	      Get_Tmin_Tmax("3pt", Ens_tag[i_ens], corr,xg , "V");
	      int T_min_V= corr.Tmin;
	      int T_max_V= corr.Tmax;
	      distr_t fit_result_V = corr.Fit_distr(bar_R_V_distr);
	      distr_t fit_result_R_V = corr.Fit_distr(R_V_distr);
	      Get_Tmin_Tmax("3pt", Ens_tag[i_ens], corr,xg , "A");
	      int T_min_A= corr.Tmin;
	      int T_max_A=corr.Tmax;
	      distr_t fit_result_A = corr.Fit_distr(F_A_distr);
	      distr_t fit_result_R_A = corr.Fit_distr(R_A_distr);
	      // #####################################
      

      
	      //PRINT R-estimators and form factors to file
	      //###############################################################################
	      VVfloat To_print_bar_R({bar_R_A_distr.ave(), bar_R_A_distr.err(), bar_R_A_distr_exp.ave(), bar_R_A_distr_exp.err(), bar_R_V_distr.ave(), bar_R_V_distr.err(), R_A_distr.ave(), R_A_distr.err(), R_V_distr.ave(), R_V_distr.err()});
	       
	     
	      Print_To_File({}, To_print_bar_R, "../data/form_factors/"+Meson+"/R_"+Ens_tag[i_ens]+"_"+mom3_l.mom[i_ens][icomb3pt].name()+"_smlev_"+to_string(sm_lev)+"_k0_noise_"+to_string(Include_k0_noise)+".dat", "", "#t       bar_RA    bar_RA_err   bar_R_A_exp     bar_R_A_exp_err    bar_RV      bar_RV_err    RA      RA_err        RV      RV_err");
	      Print_To_File({}, {F_A_distr.ave(), F_A_distr.err(), F_A_distr_exp.ave(), F_A_distr_exp.err(), bar_R_V_distr.ave(), bar_R_V_distr.err() } , "../data/form_factors/"+Meson+"/F_"+Ens_tag[i_ens]+"_"+mom3_l.mom[i_ens][icomb3pt].name()+"_smlev_"+to_string(sm_lev)+"_k0_noise_"+to_string(Include_k0_noise)+".dat", "", "#t       FA    FA_err   FA_exp    FA_exp_err     FV      FV_err");
	      if(xg > 0) { //print form factors only for xg > 0
		ofstream Print_Fitted_form_factors("../data/form_factors/"+Meson+"/FORM_FACTOR_LIST/axial_form_factors_list_"+Ens_tag[i_ens]+".dat",ofstream::app);
		Print_Fitted_form_factors<<xg<<setw(20)<<offsh<<setw(20)<<sm_lev<<setw(20)<<fit_result_A.ave()<<setw(20)<<fit_result_A.err()<<setw(20)<<"["<<T_min_A<<","<<T_max_A<<"]"<<endl;
		Print_Fitted_form_factors.close();
		Print_Fitted_form_factors.open("../data/form_factors/"+Meson+"/FORM_FACTOR_LIST/vector_form_factors_list_"+Ens_tag[i_ens]+".dat",ofstream::app);
		Print_Fitted_form_factors<<xg<<setw(20)<<offsh<<setw(20)<<sm_lev<<setw(20)<<fit_result_V.ave()<<setw(20)<<fit_result_V.err()<<setw(20)<<"["<<T_min_V<<","<<T_max_V<<"]"<<endl;
		Print_Fitted_form_factors.close();
	      }
	      //###############################################################################
      
	      
	      //plot form factors
	      //##################################################################################################
	      if(xg > 1.0e-7) { //plot only if xg > 0
		Plot_form_factors("V", bar_R_V_distr, fit_result_V, T_min_V, T_max_V, Nt, Ens_tag[i_ens], xg, offsh, sm_lev);
		Plot_form_factors("A", F_A_distr, fit_result_A, T_min_A, T_max_A, Nt, Ens_tag[i_ens], xg, offsh, sm_lev);
	      }
	      //##################################################################################################
      

	      //FIT EXPONENTIAL CONTAMINATIONS FROM R and bar_R 
	      if(Determine_contaminations && xg > 1.0e-7) { //fit contaminations only if xg > 0
	   
		Fit_contaminations("FA", F_A_distr, fit_result_A, Ens_tag[i_ens], Nt, xg, offsh, sm_lev);
		Fit_contaminations("FV", bar_R_V_distr, fit_result_V, Ens_tag[i_ens], Nt, xg, offsh, sm_lev);
		Fit_contaminations("RA", R_A_distr, fit_result_R_A, Ens_tag[i_ens], Nt, xg, offsh, sm_lev);
		Fit_contaminations("RV", R_V_distr, fit_result_R_V, Ens_tag[i_ens],  Nt, xg, offsh, sm_lev);
		//Fit_contaminations_from_derivative("FA", F_A_distr, fit_result_A, Ens_tag[i_ens], Nt, xg, offsh, sm_lev);
		//Fit_contaminations_from_derivative("FV", bar_R_V_distr, fit_result_V, Ens_tag[i_ens], Nt, xg, offsh, sm_lev);
		//Fit_contaminations_from_derivative("RA", R_A_distr, fit_result_R_A, Ens_tag[i_ens], Nt, xg, offsh, sm_lev);
		//Fit_contaminations_from_derivative("RV", R_V_distr, fit_result_R_V, Ens_tag[i_ens],  Nt, xg, offsh, sm_lev); 
	      }
	    }

	    //determine form_factors H1 H2 e FA_offsh (only valid if meson at rest and k=kx)
	    if( xg > 1.0e-7 && mom3_l.mom[i_ens][icomb3pt].k()[2] !=  0 && VIRTUAL_RUN) { //only if xg, kz != 0



	      
	      //distr_t_list H_resc= (m_fit_distr[pt2_p]/sqrt_overlap_fit_distr[pt2_p])*distr_t_list::f_of_distr(Exp_lat, m_fit_distr[pt2_p], Nt)*exp_Nt_t;
	      //distr_t_list H_resc_0 = (m_fit_distr[pt2_k0p0]/sqrt_overlap_fit_distr[pt2_k0p0])*distr_t_list::f_of_distr(Exp_lat, m_fit_distr[pt2_k0p0], Nt);


	      //#####################UNIMPROVED ESTIMATORS###########################################
	      
	       distr_t_list H1 = H_1(distr_A[0][3]*exp_Nt_t -SUB_ZERO_MOMENTUM_VIRTUAL*distr_A_kz0_k2[0][3]*exp_Nt_t_k0 , 0.5*(distr_A[1][1]+distr_A[2][2])*exp_Nt_t, distr_A[3][3]*exp_Nt_t, distr_A_k0[0][3], 0.5*(distr_A_k0[1][1]+distr_A_k0[2][2]), distr_A_k0[3][3], mom3_l.mom[i_ens][icomb3pt], m_fit_distr[pt2_k0p0])*(-2.0*fp_fit_distr[pt2_k0p0]/(distr_A_k0[1][1]+distr_A_k0[2][2]));
	      
	      distr_t_list H2 = H_2(distr_A[0][3]*exp_Nt_t -SUB_ZERO_MOMENTUM_VIRTUAL*distr_A_kz0_k2[0][3]*exp_Nt_t_k0, 0.5*(distr_A[1][1]+distr_A[2][2])*exp_Nt_t, distr_A[3][3]*exp_Nt_t, distr_A_k0[0][3], 0.5*(distr_A_k0[1][1]+distr_A_k0[2][2]), distr_A_k0[3][3], mom3_l.mom[i_ens][icomb3pt], m_fit_distr[pt2_k0p0])*(-2.0*fp_fit_distr[pt2_k0p0]/(distr_A_k0[1][1]+distr_A_k0[2][2]));
	      distr_t_list FAoff = FA_off(distr_A[0][3]*exp_Nt_t-SUB_ZERO_MOMENTUM_VIRTUAL*distr_A_kz0_k2[0][3]*exp_Nt_t_k0, 0.5*(distr_A[1][1]+distr_A[2][2])*exp_Nt_t, distr_A[3][3]*exp_Nt_t, distr_A_k0[0][3], 0.5*(distr_A_k0[1][1]+distr_A_k0[2][2]), distr_A_k0[3][3], mom3_l.mom[i_ens][icomb3pt], m_fit_distr[pt2_k0p0])*(-2.0*fp_fit_distr[pt2_k0p0]/(distr_A_k0[1][1]+distr_A_k0[2][2]));
	      
	      //#####################################################################################


	      

	      //###############IMPROVED ESTIMATORS###################################################à
	      
	      distr_t_list H1_impr = H_1_impr(distr_A[0][3]*exp_Nt_t-SUB_ZERO_MOMENTUM_VIRTUAL*distr_A_kz0_k2[0][3]*exp_Nt_t_k0, distr_A[3][0]*exp_Nt_t-SUB_ZERO_MOMENTUM_VIRTUAL*distr_A_kz0_k2[3][0]*exp_Nt_t_k0, 0.5*(distr_A[1][1]+distr_A[2][2])*exp_Nt_t, distr_A[3][3]*exp_Nt_t,  0.5*(distr_A_k0[1][1]+distr_A_k0[2][2]), distr_A_k0[3][3], mom3_l.mom[i_ens][icomb3pt], m_fit_distr[pt2_k0p0])*(-1.0*fp_fit_distr[pt2_k0p0]);

      
	      distr_t_list H2_impr = H_2_impr(distr_A[0][3]*exp_Nt_t-SUB_ZERO_MOMENTUM_VIRTUAL*distr_A_kz0_k2[0][3]*exp_Nt_t_k0, distr_A[3][0]*exp_Nt_t-SUB_ZERO_MOMENTUM_VIRTUAL*distr_A_kz0_k2[3][0]*exp_Nt_t_k0, 0.5*(distr_A[1][1]+distr_A[2][2])*exp_Nt_t, distr_A[3][3]*exp_Nt_t,  0.5*(distr_A_k0[1][1]+distr_A_k0[2][2]), distr_A_k0[3][3], mom3_l.mom[i_ens][icomb3pt], m_fit_distr[pt2_k0p0])*(-1.0*fp_fit_distr[pt2_k0p0]);
	      distr_t_list FAoff_impr = FA_off_impr(distr_A[0][3]*exp_Nt_t-SUB_ZERO_MOMENTUM_VIRTUAL*((distr_A_kz0_k2[0][3])*exp_Nt_t_k0), distr_A[3][0]*exp_Nt_t-SUB_ZERO_MOMENTUM_VIRTUAL*((distr_A_kz0_k2[3][0])*exp_Nt_t_k0), 0.5*(distr_A[1][1]+distr_A[2][2])*exp_Nt_t, distr_A[3][3]*exp_Nt_t,  0.5*(distr_A_k0[1][1]+distr_A_k0[2][2]), distr_A_k0[3][3], mom3_l.mom[i_ens][icomb3pt], m_fit_distr[pt2_k0p0])*(-1.0*fp_fit_distr[pt2_k0p0]);
	      
	      //######################################################################################


	       //###############MIXED 11 AND 33 ESTIMATORS###################################################
	      
	      distr_t_list H1_mixed_diag = H_1_mixed_diag(distr_A[0][3]*exp_Nt_t-SUB_ZERO_MOMENTUM_VIRTUAL*distr_A_kz0_k2[0][3]*exp_Nt_t_k0, distr_A[3][0]*exp_Nt_t-SUB_ZERO_MOMENTUM_VIRTUAL*distr_A_kz0_k2[3][0]*exp_Nt_t_k0, 0.5*(distr_A[1][1]+distr_A[2][2])*exp_Nt_t, distr_A[3][3]*exp_Nt_t,  0.5*(distr_A_k0[1][1]+distr_A_k0[2][2]), distr_A_k0[3][3], 1.0*SUB_ZERO_MOMENTUM_VIRTUAL*( distr_A_kz0_k2[3][3]*exp_Nt_t_k0),  1.0*SUB_ZERO_MOMENTUM_VIRTUAL*( 0.5*(distr_A_kz0_k2[1][1]+ distr_A_kz0_k2[2][2] )*exp_Nt_t_k0), mom3_l.mom[i_ens][icomb3pt], m_fit_distr[pt2_k0p0])*(-1.0*fp_fit_distr[pt2_k0p0]);
	      distr_t_list H2_mixed_diag = H_2_mixed_diag(distr_A[0][3]*exp_Nt_t-SUB_ZERO_MOMENTUM_VIRTUAL*distr_A_kz0_k2[0][3]*exp_Nt_t_k0, distr_A[3][0]*exp_Nt_t-SUB_ZERO_MOMENTUM_VIRTUAL*distr_A_kz0_k2[3][0]*exp_Nt_t_k0, 0.5*(distr_A[1][1]+distr_A[2][2])*exp_Nt_t, distr_A[3][3]*exp_Nt_t,  0.5*(distr_A_k0[1][1]+distr_A_k0[2][2]),distr_A_k0[3][3], 1.0*SUB_ZERO_MOMENTUM_VIRTUAL*( distr_A_kz0_k2[3][3]*exp_Nt_t_k0),  1.0*SUB_ZERO_MOMENTUM_VIRTUAL*( 0.5*(distr_A_kz0_k2[1][1]+ distr_A_kz0_k2[2][2] )*exp_Nt_t_k0) , mom3_l.mom[i_ens][icomb3pt], m_fit_distr[pt2_k0p0])*(-1.0*fp_fit_distr[pt2_k0p0]);
	      distr_t_list FAoff_mixed_diag = FA_off_mixed_diag(distr_A[0][3]*exp_Nt_t-SUB_ZERO_MOMENTUM_VIRTUAL*distr_A_kz0_k2[0][3]*exp_Nt_t_k0, distr_A[3][0]*exp_Nt_t-SUB_ZERO_MOMENTUM_VIRTUAL*distr_A_kz0_k2[3][0]*exp_Nt_t_k0, 0.5*(distr_A[1][1]+distr_A[2][2])*exp_Nt_t, distr_A[3][3]*exp_Nt_t,  0.5*(distr_A_k0[1][1]+distr_A_k0[2][2]), distr_A_k0[3][3], 1.0*SUB_ZERO_MOMENTUM_VIRTUAL*( distr_A_kz0_k2[3][3]*exp_Nt_t_k0),  1.0*SUB_ZERO_MOMENTUM_VIRTUAL*( 0.5*(distr_A_kz0_k2[1][1]+ distr_A_kz0_k2[2][2] )*exp_Nt_t_k0) , mom3_l.mom[i_ens][icomb3pt], m_fit_distr[pt2_k0p0])*(-1.0*fp_fit_distr[pt2_k0p0]);
	      
	     
	      distr_t_list FVoff = (V_ave_unpolarized(distr_V,distr_V_k0, m_fit_distr[pt2_p], mom3_l.mom[i_ens][icomb3pt],Nt/2)/A_ave_unpolarized(distr_A_k0))*m_fit_distr[pt2_k0p0]*fp_fit_distr[pt2_k0p0];

	      //###########################################

	      //fit virtual FF
	      Get_virtual_ff_fit_interval_v1("H1", offsh, c.tht[2], corr.Tmin, corr.Tmax);
	      distr_t H1_fit = corr.Fit_distr(H1);
	      distr_t H1_impr_fit = corr.Fit_distr(H1_impr);
	      distr_t H1_mixed_diag_fit= corr.Fit_distr(H1_mixed_diag);
	      Get_virtual_ff_fit_interval_v1("H2", offsh, c.tht[2], corr.Tmin, corr.Tmax);
	      distr_t H2_fit = corr.Fit_distr(H2);
	      distr_t H2_impr_fit = corr.Fit_distr(H2_impr);
	      distr_t H2_mixed_diag_fit= corr.Fit_distr(H2_mixed_diag);
	      Get_virtual_ff_fit_interval_v1("FA", offsh, c.tht[2], corr.Tmin, corr.Tmax);
	      distr_t FAoff_fit = corr.Fit_distr(FAoff);
	      distr_t FAoff_impr_fit = corr.Fit_distr(FAoff_impr);
	      distr_t FAoff_mixed_diag_fit= corr.Fit_distr(FAoff_mixed_diag);
	      Get_virtual_ff_fit_interval_v1("FV", offsh, c.tht[2], corr.Tmin, corr.Tmax);
	      distr_t FVoff_fit = corr.Fit_distr(FVoff);

	   
	      	   

	      //IF FIT_FORM_FACTORS ADD THEM TO FF_LIST
	     
	      if(FIT_VIRTUAL_FF) {
		mom3_l_analyzed_only[i_ens].push_back(mom3_l.mom[i_ens][icomb3pt]);
		  FV_off_list[i_ens].push_back(FVoff);
		  FV_off_const_fit_list[i_ens].push_back(FVoff_fit);
		  if(VIRTUAL_ESTIMATOR_SET==0) {
		  H1_list[i_ens].push_back(H1);
		  H2_list[i_ens].push_back(H2);
		  FA_off_list[i_ens].push_back(FAoff);
		  H1_const_fit_list[i_ens].push_back(H1_fit);
		  H2_const_fit_list[i_ens].push_back(H2_fit);
		  FA_off_const_fit_list[i_ens].push_back(FAoff_fit);
		  }
		else if(VIRTUAL_ESTIMATOR_SET==1) {
		  H1_list[i_ens].push_back(H1_impr);
		  H2_list[i_ens].push_back(H2_impr);
		  FA_off_list[i_ens].push_back(FAoff_impr);
		  H1_const_fit_list[i_ens].push_back(H1_impr_fit);
		  H2_const_fit_list[i_ens].push_back(H2_impr_fit);
		  FA_off_const_fit_list[i_ens].push_back(FAoff_impr_fit);
		}
		else if(VIRTUAL_ESTIMATOR_SET==2) {
		  H1_list[i_ens].push_back(H1_impr);
		  H2_list[i_ens].push_back(H2_mixed_diag);
		  FA_off_list[i_ens].push_back(FAoff_mixed_diag);
		  H1_const_fit_list[i_ens].push_back(H1_mixed_diag_fit);
		  H2_const_fit_list[i_ens].push_back(H2_mixed_diag_fit);
		  FA_off_const_fit_list[i_ens].push_back(FAoff_mixed_diag_fit);
		}

		else if(VIRTUAL_ESTIMATOR_SET==3) {
		   H2_list[i_ens].push_back(H2_mixed_diag);
		   H1_list[i_ens].push_back(H1_mixed_diag);
		   FA_off_list[i_ens].push_back(FAoff_impr);
		   H1_const_fit_list[i_ens].push_back(H1_mixed_diag_fit);
		   H2_const_fit_list[i_ens].push_back(H2_mixed_diag_fit);
		   FA_off_const_fit_list[i_ens].push_back(FAoff_impr_fit);
		}
		else crash("VIRTUAL_ESTIMATOR_SET is neither 0,1,2,3");
	      }
	      


	      //rescale FV
	      FVoff= FVoff*Za_ov_Zv[Za_ov_Zv.size() -1];
	      FVoff_fit= FVoff_fit*Za_ov_Zv[Za_ov_Zv.size()-1];

	      //print to file

	      ofstream Print_Fitted_virtual_form_factors("../data/form_factors/"+Meson+"/VIRTUAL_FORM_FACTOR_LIST/"+Ens_tag[i_ens]+".dat",ofstream::app);
	      Print_Fitted_virtual_form_factors<<xk<<setw(20)<<xq<<setw(20)<<H1_fit.ave()<<setw(20)<<H1_fit.err()<<setw(20)<<H2_fit.ave()<<setw(20)<<H2_fit.err()<<setw(20)<<FAoff_fit.ave()<<setw(20)<<FAoff_fit.err()<<setw(20)<<FVoff_fit.ave()<<setw(20)<<FVoff_fit.err()<<endl;
	      Print_Fitted_virtual_form_factors.close();
	      Print_Fitted_virtual_form_factors.open("../data/form_factors/"+Meson+"/VIRTUAL_FORM_FACTOR_LIST/"+Ens_tag[i_ens]+"_impr.dat",ofstream::app);
	      Print_Fitted_virtual_form_factors<<xk<<setw(20)<<xq<<setw(20)<<H1_impr_fit.ave()<<setw(20)<<H1_impr_fit.err()<<setw(20)<<H2_impr_fit.ave()<<setw(20)<<H2_impr_fit.err()<<setw(20)<<FAoff_impr_fit.ave()<<setw(20)<<FAoff_impr_fit.err()<<setw(20)<<FVoff_fit.ave()<<setw(20)<<FVoff_fit.err()<<endl;
	      Print_Fitted_virtual_form_factors.close();
	      Print_Fitted_virtual_form_factors.open("../data/form_factors/"+Meson+"/VIRTUAL_FORM_FACTOR_LIST/"+Ens_tag[i_ens]+"_mixed_diag.dat",ofstream::app);
	      Print_Fitted_virtual_form_factors<<xk<<setw(20)<<xq<<setw(20)<<H1_mixed_diag_fit.ave()<<setw(20)<<H1_mixed_diag_fit.err()<<setw(20)<<H2_mixed_diag_fit.ave()<<setw(20)<<H2_mixed_diag_fit.err()<<setw(20)<<FAoff_mixed_diag_fit.ave()<<setw(20)<<FAoff_mixed_diag_fit.err()<<setw(20)<<FVoff_fit.ave()<<setw(20)<<FVoff_fit.err()<<endl;
	      Print_Fitted_virtual_form_factors.close();

          
	     
	      Print_To_File({}, {H1.ave(), H1.err(), H2.ave(), H2.err(), FAoff.ave(), FAoff.err(), FVoff.ave(), FVoff.err()}, "../data/form_factors/"+Meson+"/FH_off_"+Ens_tag[i_ens]+"_"+mom3_l.mom[i_ens][icomb3pt].name()+"_smlev_"+to_string(sm_lev)+"_k0_noise_"+to_string(Include_k0_noise)+".dat", "", "#t       H1    H1_err     H2      H2_err    FA_off      FA_off_err     FV_off     FV_off_err ");
	      Print_To_File({}, {H1_impr.ave(), H1_impr.err(), H2_impr.ave(), H2_impr.err(), FAoff_impr.ave(), FAoff_impr.err(), FVoff.ave(), FVoff.err()}, "../data/form_factors/"+Meson+"/FH_off_impr_"+Ens_tag[i_ens]+"_"+mom3_l.mom[i_ens][icomb3pt].name()+"_smlev_"+to_string(sm_lev)+"_k0_noise_"+to_string(Include_k0_noise)+".dat", "", "#t       H1    H1_err     H2      H2_err    FA_off      FA_off_err     FV_off     FV_off_err ");
	       Print_To_File({}, {H1_mixed_diag.ave(), H1_mixed_diag.err(), H2_mixed_diag.ave(), H2_mixed_diag.err(), FAoff_mixed_diag.ave(), FAoff_mixed_diag.err(), FVoff.ave(), FVoff.err()}, "../data/form_factors/"+Meson+"/FH_off_mixed_diag_"+Ens_tag[i_ens]+"_"+mom3_l.mom[i_ens][icomb3pt].name()+"_smlev_"+to_string(sm_lev)+"_k0_noise_"+to_string(Include_k0_noise)+".dat", "", "#t       H1    H1_err     H2      H2_err    FA_off      FA_off_err     FV_off     FV_off_err ");
	      //fit in time interval
	     
	      //save in file
	      if(xg > 1.0e-7) { //print form factors only for xg > 0
		ofstream Print_virtual_form_factors("../data/form_factors/"+Meson+"/FORM_FACTOR_LIST/virtual_form_factors_list_"+Ens_tag[i_ens]+".dat",ofstream::app);
		Print_virtual_form_factors<<xg<<setw(20)<<offsh<<setw(20)<<sm_lev<<setw(20)<<H1_fit.ave()<<setw(20)<<H1_fit.err()<<setw(20)<<H2_fit.ave()<<setw(20)<<H2_fit.err()<<setw(20)<<FAoff_fit.ave()<<setw(20)<<FAoff_fit.err()<<setw(20)<<FVoff_fit.ave()<<setw(20)<<FVoff_fit.err()<<endl;
		Print_virtual_form_factors.close();


		//print subtracted hadronic tensor
	
		
		double Eg= mom3_l.mom[i_ens][icomb3pt].Egamma();
		double kz= mom3_l.mom[i_ens][icomb3pt].k()[2];
		double ksq = pow(mom3_l.mom[i_ens][icomb3pt].virt(),2);
		distr_t m=m_fit_distr[pt2_k0p0];
		distr_t_list H_resc_0 = (m_fit_distr[pt2_k0p0]/sqrt_overlap_fit_distr[pt2_k0p0])*distr_t_list::f_of_distr(Exp_lat, m_fit_distr[pt2_k0p0], Nt);
		distr_t_list H_resc= (m_fit_distr[pt2_p]/sqrt_overlap_fit_distr[pt2_p])*distr_t_list::f_of_distr(Exp_lat, m_fit_distr[pt2_p], Nt)*exp_Nt_t;

		distr_t_list H30_sub = distr_A[0][3]*H_resc - distr_A_k0[1][1]*H_resc_0*(m-Eg)/( (2.0*m*Eg/kz) -(ksq/kz));
		distr_t_list H11_sub = distr_A[1][1]*H_resc - distr_A_k0[1][1]*H_resc_0;
		distr_t_list H33_sub = distr_A[3][3]*H_resc - distr_A_k0[3][3]*H_resc_0*(2.0*m*Eg -ksq -pow(kz,2))/(2.0*m*Eg - ksq);
		distr_t_list H03_sub = distr_A[3][0]*H_resc - distr_A_k0[1][1]*H_resc_0*(2.0*m -Eg)/(2.0*m*(Eg/kz) -(ksq/kz));
		distr_t_list H00_sub = distr_A[0][0]*H_resc - distr_A_k0[0][0]*H_resc_0*(2.0*m -Eg + pow(Eg,2)/m -ksq/m)/(2.0*m -ksq/Eg);

		Print_To_File({}, {H30_sub.ave(), H30_sub.err(), H11_sub.ave(), H11_sub.err(), H33_sub.ave(), H33_sub.err(), H03_sub.ave(), H03_sub.err(), H00_sub.ave(), H00_sub.err()},  "../data/form_factors/"+Meson+"/Hadr_tens_subtracted_"+Ens_tag[i_ens]+"_"+mom3_l.mom[i_ens][icomb3pt].name()+"_smlev_"+to_string(sm_lev)+"_k0_noise_"+to_string(Include_k0_noise)+".dat", "", "#t     H30         H30_err        H11        H11_err     H33       H33_err       H03         H03_err            H00               H00_err");
        
		
	      }

	      
	    
	    }

	  
	   
	  }
	  //corr.Perform_Nt_t_average=1;
	}
	
	
      }



      
    }
    fclose(stream_3pt);
    fclose(stream_2pt);
  }

 

  if(FIT_VIRTUAL_FF) {

    Vfloat MC_ee_ChPT, MC_mumu_ChPT, QUAD_ee_ChPT, QUAD_mumu_ChPT,MC_ee_VMD, MC_mumu_VMD, QUAD_ee_VMD, QUAD_mumu_VMD;
    Vfloat MC_ee_ChPT_err, MC_mumu_ChPT_err, QUAD_ee_ChPT_err, QUAD_mumu_ChPT_err,MC_ee_VMD_err, MC_mumu_VMD_err, QUAD_ee_VMD_err, QUAD_mumu_VMD_err;
    Vfloat MC_eee_ChPT, MC_mumumu_ChPT, MC_eee_VMD, MC_mumumu_VMD;
    Vfloat MC_eee_ChPT_err, MC_mumumu_ChPT_err, MC_eee_VMD_err, MC_mumumu_VMD_err;
    Vfloat times;

    if(H1_list.size() == 0 || H2_list.size() == 0 || FA_off_list.size() == 0 || FV_off_list.size() == 0 || f_p.size() == 0 || m_p.size() == 0) crash("At least one between ff,m_p,f_p list has size zero. Exiting...");
    if(H1_const_fit_list.size() == 0 || H2_const_fit_list.size() == 0 || FA_off_const_fit_list.size() == 0 || FV_off_const_fit_list.size() == 0) crash("At least one FF const fit list has size zero. Exiting...");
    int t_min_ff_extr, t_max_ff_extr;
    string Tag=Ens_tag[0];
    if(!USE_FITTED_FF) {
    t_min_ff_extr=5;
    t_max_ff_extr=22;
    }
    else {t_min_ff_extr=0;t_max_ff_extr=0;}
    for(int t=t_min_ff_extr;t<=t_max_ff_extr;t++) {
      vector<distr_t> H1_temp_list, H2_temp_list, FA_off_temp_list, FV_off_temp_list;
     

      distr_t f_p_temp= f_p[0];
      distr_t m_p_temp= m_p[0];
      distr_t Za_ov_Zv_temp= Za_ov_Zv[0];

      //I shall assume now that we analyze a single ensemble
      if(!USE_FITTED_FF) {
	for(int icomb=0;icomb<(signed)H1_list[0].size();icomb++) {
	H1_temp_list.push_back(H1_list[0][icomb].distr_list[t]);
	H2_temp_list.push_back(H2_list[0][icomb].distr_list[t]);
	FA_off_temp_list.push_back(FA_off_list[0][icomb].distr_list[t]);
	FV_off_temp_list.push_back(FV_off_list[0][icomb].distr_list[t]);
	}
      }
      else { //forward result of constant fit
	H1_temp_list = H1_const_fit_list[0];
	H2_temp_list = H2_const_fit_list[0];
	FA_off_temp_list = FA_off_const_fit_list[0];
	FV_off_temp_list = FV_off_const_fit_list[0];
      }
     
	if(mom3_l_analyzed_only.size() == 0) crash("3pt_momenta list to use in FIT_VIRTUAL_FF has size 0");
	vector<pt3_momenta> Tmom= mom3_l_analyzed_only[0];
	
      
	vector<function<double(double xk, double xq)>> H1_VMD, H2_VMD, FA_off_VMD, FV_off_VMD, H1_ChPT, H2_ChPT, FA_off_ChPT, FV_off_ChPT;
	//define jack_distr for fit parameters
	//polynomial fit
	distr_t a0_H1, ak_H1, aq_H1; //H1
	distr_t a0_H2, ak_H2, aq_H2; //H2
	distr_t a0_FA, ak_FA, aq_FA; //FA
	distr_t a0_FV, ak_FV, aq_FV; //FV
	//pole-like fit
	distr_t A_H1, rq_H1, rk_H1; //H1
	distr_t A_H2, rq_H2, rk_H2; //H2
	distr_t A_FA, rq_FA, rk_FA; //FA
	distr_t A_FV, rq_FV, rk_FV; //FV
	
	cout<<"####Fitting form factor H1  VMD Ansatz.....####"<<endl;
	cout<<"USE_FITTED_FF: "<<USE_FITTED_FF;
	if(!USE_FITTED_FF) cout<<"time: "<<t<<endl;
	else cout<<endl;
	//
	Fit_virtual_FF_VMD(H1_VMD, A_H1, rk_H1, rq_H1,  H1_temp_list, f_p_temp, m_p_temp, Za_ov_Zv_temp,  Tmom, "H1", "A", Tag, Meson, UseJack, USE_FITTED_FF, t);// VECTOR MESON DOMINANCE FIT
	cout<<"####Fitting form factor H2  VMD Ansatz.....####"<<endl;
	cout<<"USE_FITTED_FF: "<<USE_FITTED_FF;
	if(!USE_FITTED_FF) cout<<"time: "<<t<<endl;
	else cout<<endl;
	//
	Fit_virtual_FF_VMD(H2_VMD, A_H2, rk_H2, rq_H2 , H2_temp_list, f_p_temp, m_p_temp, Za_ov_Zv_temp,Tmom, "H2", "A", Tag, Meson, UseJack,  USE_FITTED_FF, t);// VECTOR MESON DOMINANCE FIT
	cout<<"####Fitting form factor FA  VMD Ansatz.....####"<<endl;
	cout<<"USE_FITTED_FF: "<<USE_FITTED_FF;
	if(!USE_FITTED_FF) cout<<"time: "<<t<<endl;
	else cout<<endl;
	//
	Fit_virtual_FF_VMD(FA_off_VMD, A_FA, rk_FA, rq_FA, FA_off_temp_list, f_p_temp, m_p_temp, Za_ov_Zv_temp, Tmom, "FA", "A", Tag, Meson, UseJack,  USE_FITTED_FF, t);// VECTOR MESON DOMINANCE FIT
	cout<<"####Fitting form factor FV  VMD Ansatz.....####"<<endl;
	cout<<"USE_FITTED_FF: "<<USE_FITTED_FF;
	if(!USE_FITTED_FF) cout<<"time: "<<t<<endl;
	else cout<<endl;
	//
	Fit_virtual_FF_VMD(FV_off_VMD, A_FV, rk_FV, rq_FV, FV_off_temp_list, f_p_temp, m_p_temp, Za_ov_Zv_temp,  Tmom, "FV", "V", Tag, Meson, UseJack,  USE_FITTED_FF, t);// VECTOR MESON DOMINANCE FIT
	cout<<"####Fitting form factor H1  ChPT Ansatz.....####"<<endl;
	cout<<"USE_FITTED_FF: "<<USE_FITTED_FF;
	if(!USE_FITTED_FF) cout<<"time: "<<t<<endl;
	else cout<<endl;
	//
	Fit_virtual_FF_ChPT(H1_ChPT, a0_H1, ak_H1, aq_H1, H1_temp_list, f_p_temp, m_p_temp, Za_ov_Zv_temp,  Tmom, "H1", "A", Tag, Meson, UseJack,  USE_FITTED_FF, t);// CHPT FIT
	cout<<"####Fitting form factor H2  ChPT Ansatz.....####"<<endl;
	cout<<"USE_FITTED_FF: "<<USE_FITTED_FF;
	if(!USE_FITTED_FF) cout<<"time: "<<t<<endl;
	else cout<<endl;
	//
	Fit_virtual_FF_ChPT(H2_ChPT, a0_H2, ak_H2, aq_H2, H2_temp_list, f_p_temp, m_p_temp, Za_ov_Zv_temp,  Tmom, "H2", "A", Tag, Meson, UseJack,  USE_FITTED_FF, t);// CHPT FIT
	cout<<"####Fitting form factor FA  ChPT Ansatz.....####"<<endl;
	cout<<"USE_FITTED_FF: "<<USE_FITTED_FF;
	if(!USE_FITTED_FF) cout<<"time: "<<t<<endl;
	else cout<<endl;
	//
	Fit_virtual_FF_ChPT(FA_off_ChPT, a0_FA, ak_FA, aq_FA, FA_off_temp_list, f_p_temp, m_p_temp, Za_ov_Zv_temp,  Tmom, "FA", "A", Tag, Meson, UseJack,  USE_FITTED_FF, t);// CHPT FIT
	cout<<"####Fitting form factor FV  ChPT Ansatz.....####"<<endl;
	cout<<"USE_FITTED_FF: "<<USE_FITTED_FF;
	if(!USE_FITTED_FF) cout<<"time: "<<t<<endl;
	else cout<<endl;
	//
	Fit_virtual_FF_ChPT(FV_off_ChPT, a0_FV, ak_FV, aq_FV, FV_off_temp_list, f_p_temp, m_p_temp, Za_ov_Zv_temp,  Tmom, "FV", "V", Tag, Meson, UseJack,  USE_FITTED_FF, t);// CHPT FIT


	//print correlation matrix of fit parameters

	//pole-like fit
	VVfloat Corr_mat_pole_like;
	vector<distr_t> fit_pars_pole_like({ A_H1, rk_H1, rq_H1, A_H2, rk_H2, rq_H2, A_FA, rk_FA, rq_FA, A_FV, rk_FV, rq_FV});
	int N= (signed)fit_pars_pole_like.size();

	//set dimensions of Corr_mat_pole_like   ( NxN matrix )
	Corr_mat_pole_like.resize(N);
	for( auto &corr_mat_row: Corr_mat_pole_like) corr_mat_row.resize(N);

	cout<<"Printing correlation matrix for pole-like fit"<<endl;
	for(int i=0; i < N; i++) {
	  cout<<"( "<<endl;
	  for(int j=0; j < N; j++) {
	    Corr_mat_pole_like[i][j] = (fit_pars_pole_like[i]%fit_pars_pole_like[j])/( sqrt( ( fit_pars_pole_like[i]%fit_pars_pole_like[i])*( fit_pars_pole_like[j]%fit_pars_pole_like[j]         )));
	    cout<<Corr_mat_pole_like[i][j]<<", ";
	  }
	  cout<<" )"<<endl;
	}
	cout<<"done!"<<endl;
	//polynomial fit
	VVfloat Corr_mat_polynomial;
	vector<distr_t> fit_pars_polynomial({ a0_H1, ak_H1, aq_H1, a0_H2, ak_H2, aq_H2, a0_FA, ak_FA, aq_FA, a0_FV, ak_FV, aq_FV});
	N= (signed)fit_pars_polynomial.size();

	//set dimensions of Corr_mat_pole_like   ( NxN matrix )
	Corr_mat_polynomial.resize(N);
	for( auto &corr_mat_row: Corr_mat_polynomial) corr_mat_row.resize(N);
	cout<<"Printing correlation matrix for polynomial fit"<<endl;
	for(int i=0; i < N; i++) {
	  cout<<"( "<<endl;
	  for(int j=0; j < N; j++) {
	    Corr_mat_polynomial[i][j] = (fit_pars_polynomial[i]%fit_pars_polynomial[j])/( sqrt( ( fit_pars_polynomial[i]%fit_pars_polynomial[i])*( fit_pars_polynomial[j]%fit_pars_polynomial[j]         )));
	    cout<<Corr_mat_polynomial[i][j]<<", ";
	  }
	  cout<<" )"<<endl;
	}
	cout<<"done!"<<endl;


	


	

	


	
	if(COMPUTE_l4_DECAY_RATE) {
	  cout<<"#####STARTING RATE COMPUTATION#####"<<endl;
	  cout<<"USE_FITTED_FF: "<<USE_FITTED_FF;
	  if(!USE_FITTED_FF) cout<<"time: "<<t<<endl;
	  else cout<<endl;
	  cout<<"#####COMPUTATION OF THE RATE USING VMD ANSATZ#####"<<endl;
	  //set ff to zero to compute tree level contribution

	  
	  auto zero_func= [](double xk, double xq) -> double {return 0.0;};
	  auto id_func = [](double xk, double xq) -> double {return 1.0;};
	  auto id_func2 = [](double xk, double xq) -> double {return 20.0;};


	  
	  /* only to make tests
	  for(unsigned int i=0;i<H1_VMD.size();i++) {
	    H1_VMD[i]= zero_func;
	    H2_VMD[i]= zero_func;
	    FA_off_VMD[i]= zero_func;
	    FV_off_VMD[i]= id_func;
	   
	  }
	  only to make tests */ 
	  
	 
	  distr_t f_p_temp= f_p[0];
	  distr_t m_p_temp= m_p[0];
	  cout<<"####COMPUTATION OF THE RATE USING ChPT ANSATZ#####"<<endl;
	  Decay_Rate_Integration_Result RATE_ChPT= Num_Integrate_Decay_Rate(H1_ChPT, H2_ChPT, FA_off_ChPT, FV_off_ChPT, m_p_temp, f_p_temp, UseJack, Meson, false);
	  cout<<"####COMPUTATION OF THE RATE USING ChPT ANSATZ#####"<<endl;
	  Decay_Rate_Integration_Result RATE_VMD = Num_Integrate_Decay_Rate(H1_VMD, H2_VMD, FA_off_VMD, FV_off_VMD, m_p_temp, f_p_temp, UseJack, Meson, false);
	  
	  cout<<"####END RATE COMPUTATION####"<<endl;
	  cout<<"####FINAL RATE ESTIMATES:####"<<endl;
	  cout<<"e+e- (ChPT) : (Quad) "<<RATE_ChPT.Int_Quad_val_ee<<"("<<RATE_ChPT.Int_Quad_err_ee<<"),  (MC) "<<RATE_ChPT.Int_MonteCarlo_val_ee<<"("<<RATE_ChPT.Int_MonteCarlo_err_ee<<")"<<endl;
          cout<<"mu+mu- (ChPT) : (Quad) "<<RATE_ChPT.Int_Quad_val_mumu<<"("<<RATE_ChPT.Int_Quad_err_mumu<<"),  (MC) "<<RATE_ChPT.Int_MonteCarlo_val_mumu<<"("<<RATE_ChPT.Int_MonteCarlo_err_mumu<<")"<<endl;
	  cout<<"e+e-e+ (ChPT) : (MC) "<<RATE_ChPT.Int_MonteCarlo_val_eee<<"("<<RATE_ChPT.Int_MonteCarlo_err_eee<<")"<<endl;
          cout<<"mu+mu-mu+ (ChPT) : (MC) "<<RATE_ChPT.Int_MonteCarlo_val_mumumu<<"("<<RATE_ChPT.Int_MonteCarlo_err_mumumu<<")"<<endl;

	  
	  cout<<"e+e- (VMD) : (Quad) "<<RATE_VMD.Int_Quad_val_ee<<"("<<RATE_VMD.Int_Quad_err_ee<<"),  (MC) "<<RATE_VMD.Int_MonteCarlo_val_ee<<"("<<RATE_VMD.Int_MonteCarlo_err_ee<<")"<<endl;
          cout<<"mu+mu- (VMD) : (Quad) "<<RATE_VMD.Int_Quad_val_mumu<<"("<<RATE_VMD.Int_Quad_err_mumu<<"),  (MC) "<<RATE_VMD.Int_MonteCarlo_val_mumu<<"("<<RATE_VMD.Int_MonteCarlo_err_mumu<<")"<<endl;

	  
	  cout<<"e+e-e+e- (VMD) : (MC) "<<RATE_VMD.Int_MonteCarlo_val_eee<<"("<<RATE_VMD.Int_MonteCarlo_err_eee<<")"<<endl;
          cout<<"mu+mu-mu- (VMD) : (MC) "<<RATE_VMD.Int_MonteCarlo_val_mumumu<<"("<<RATE_VMD.Int_MonteCarlo_err_mumumu<<")"<<endl;

	  
	  if(!USE_FITTED_FF) {
	    //val
	    QUAD_ee_ChPT.push_back( RATE_ChPT.Int_Quad_val_ee);
	    QUAD_mumu_ChPT.push_back( RATE_ChPT.Int_Quad_val_mumu);
	    MC_ee_ChPT.push_back( RATE_ChPT.Int_MonteCarlo_val_ee);
	    MC_mumu_ChPT.push_back( RATE_ChPT.Int_MonteCarlo_val_mumu);
	    QUAD_ee_VMD.push_back(RATE_VMD.Int_Quad_val_ee);
	    QUAD_mumu_VMD.push_back( RATE_VMD.Int_Quad_val_mumu);
	    MC_ee_VMD.push_back( RATE_VMD.Int_MonteCarlo_val_ee);
	    MC_mumu_VMD.push_back( RATE_VMD.Int_MonteCarlo_val_mumu);
	    MC_eee_ChPT.push_back(RATE_ChPT.Int_MonteCarlo_val_eee);
	    MC_mumumu_ChPT.push_back(RATE_ChPT.Int_MonteCarlo_val_mumumu);
	    MC_eee_VMD.push_back(RATE_VMD.Int_MonteCarlo_val_eee);
	    MC_mumumu_VMD.push_back(RATE_VMD.Int_MonteCarlo_val_mumumu);

	    //errors
	    QUAD_ee_ChPT_err.push_back( RATE_ChPT.Int_Quad_err_ee);
	    QUAD_mumu_ChPT_err.push_back( RATE_ChPT.Int_Quad_err_mumu);
	    MC_ee_ChPT_err.push_back( RATE_ChPT.Int_MonteCarlo_err_ee);
	    MC_mumu_ChPT_err.push_back( RATE_ChPT.Int_MonteCarlo_err_mumu);
	    QUAD_ee_VMD_err.push_back(RATE_VMD.Int_Quad_err_ee);
	    QUAD_mumu_VMD_err.push_back( RATE_VMD.Int_Quad_err_mumu);
	    MC_ee_VMD_err.push_back( RATE_VMD.Int_MonteCarlo_err_ee);
	    MC_mumu_VMD_err.push_back( RATE_VMD.Int_MonteCarlo_err_mumu);
	    MC_eee_ChPT_err.push_back(RATE_ChPT.Int_MonteCarlo_err_eee);
	    MC_mumumu_ChPT_err.push_back(RATE_ChPT.Int_MonteCarlo_err_mumumu);
	    MC_eee_VMD_err.push_back(RATE_VMD.Int_MonteCarlo_err_eee);
	    MC_mumumu_VMD_err.push_back(RATE_VMD.Int_MonteCarlo_err_mumumu);
	    times.push_back((double)t);
	  }

	}
      
      
    

      
    }
    if(!USE_FITTED_FF) {   //print rates
      Print_To_File({}, {times, QUAD_ee_ChPT, QUAD_ee_ChPT_err, QUAD_ee_VMD, QUAD_ee_VMD_err, MC_ee_ChPT, MC_ee_ChPT_err, MC_ee_VMD, MC_ee_VMD_err, QUAD_mumu_ChPT, QUAD_mumu_ChPT_err, QUAD_mumu_VMD, QUAD_mumu_VMD_err, MC_mumu_ChPT, MC_mumu_ChPT_err, MC_mumu_VMD, MC_mumu_VMD_err, MC_eee_ChPT, MC_eee_ChPT_err, MC_eee_VMD, MC_eee_VMD_err, MC_mumumu_ChPT, MC_mumumu_ChPT_err, MC_mumumu_VMD, MC_mumu_VMD_err} , "../data/form_factors/"+Meson+"/decay_rate_plateaux_"+Tag+".dat", "", "#times, QUAD_ee_ChPT, QUAD_ee_ChPT_err, QUAD_ee_VMD, QUAD_ee_VMD_err, MC_ee_ChPT, MC_ee_ChPT_err, MC_ee_VMD, MC_ee_VMD_err, QUAD_mumu_ChPT, QUAD_mumu_ChPT_err, QUAD_mumu_VMD, QUAD_mumu_VMD_err, MC_mumu_ChPT, MC_mumu_ChPT_err, MC_mumu_VMD, MC_mumu_VMD_err, MC_eee_ChPT, MC_eee_ChPT_err, MC_eee_VMD, MC_eee_VMD_err, MC_mumumu_ChPT, MC_mumumu_ChPT_err, MC_mumumu_VMD, MC_mumumu_VMD_err");
    }
  }
  
  return;
}

