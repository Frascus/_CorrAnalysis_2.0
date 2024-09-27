#include "../include/Spectral.h"

double step_size = 0.01; //in units of sigma
const bool INCLUDE_ERRORS= true;
const double lambda= INCLUDE_ERRORS?0.9:0.0;
const bool FIND_OPTIMAL_LAMBDA= true;
const string COV_MATRIX_MODE = "";
const int Nmoms=0;
const int alpha=0;
const int verbosity_lev=1;
double Beta= 0.0; //1.99;
const bool mult_estimated_from_norm0=false;
const bool print_reco_in_stability_analysis=false;
const bool extended_analysis_of_syst=false;
bool Integrate_up_to_max_energy=true;
double Emax_int=4.0;
bool ONLY_FORWARD=false;
bool USE_GENERALIZED_NORM=false;


using namespace std;


double Get_exact_gauss(const double &E, const double &m , const double &s, const double &E0) {

  double e= exp(-0.5*(E-m)*(E-m)/(s*s));
  double norm= s*sqrt(2.0*M_PI) ;

  return e/norm;

}

PrecFloat Get_exact_gauss(const PrecFloat &E,const PrecFloat &m,const PrecFloat &s,const PrecFloat &E0) {
  PrecFloat e = exp( -0.5*(E-m)*(E-m)/(s*s));
  PrecFloat norm= s*sqrt(precPi()*PrecFloat(2.0)) ;
  return e/norm;
}



PrecFloat BaseFunc(const PrecFloat& E, int t, int T) { return exp(-E*t) + ((ONLY_FORWARD)?PrecFloat(0.0):exp( -E*T+ E*t));  }




PrecFloat aE0(const PrecFloat &E0,const PrecFloat &t, int n) {
  if(n!=0) crash("In AE0 n!=0");

  
  if(!Integrate_up_to_max_energy)  return exp(-E0*t)/t;



  return exp(-E0*t)/t - exp(-PrecFloat(Emax_int)*t)/t;
   
 
}

PrecFloat aE0_std(const PrecFloat &E0,const PrecFloat &t, int n) {
  if(n!=0) crash("In AE0_std n!=0");

  return exp(-E0*t)/t;


}

PrecFloat aE0_std_Emax(const PrecFloat &E0,const PrecFloat &t, int n) {
  if(n!=0) crash("In AE0_std_Emax n!=0");

  return exp(-E0*t)/t - exp(-PrecFloat(Emax_int)*t)/t;
 
 
}


void Get_Atr(PrecMatr& Atr, const PrecFloat &E0, int T, int tmin, int tmax, const function<PrecFloat( PrecFloat )> &Atr_gen_NORM )  {

  Atr.resize(tmax-tmin+1, tmax-tmin+1);

  for(int t=tmin;t<= tmax; t++)
    for(int r=tmin; r<= tmax; r++)
      {
	if(USE_GENERALIZED_NORM) {
	  if(ONLY_FORWARD) Atr(t-tmin,r-tmin)= Atr_gen_NORM( -PrecFloat(Beta)+ t+r);
	  else  Atr(t-tmin,r-tmin) = Atr_gen_NORM(-PrecFloat(Beta)+t+r) + Atr_gen_NORM(-PrecFloat(Beta)+ T-t+r) + Atr_gen_NORM(-PrecFloat(Beta)+T-r + t) + Atr_gen_NORM(-PrecFloat(Beta)+ 2*T -t -r);
	}
	else {
	  if(ONLY_FORWARD) Atr(t-tmin,r-tmin) = aE0(E0,-PrecFloat(Beta) + t+r,alpha) ;
	  else Atr(t-tmin,r-tmin) = aE0(E0,-PrecFloat(Beta) + t+r,alpha) + aE0(E0, -PrecFloat(Beta) + T - t +r, alpha) + aE0(E0,-PrecFloat(Beta)+ T+t -r, alpha) + aE0(E0,-PrecFloat(Beta)+ 2*T -t -r, alpha);
	}
      }
   
  
  return;
}

void Get_Atr_std(PrecMatr& Atr, const PrecFloat &E0, int T, int tmin, int tmax)  {

  Atr.resize(tmax-tmin+1, tmax-tmin+1);

  for(int t=tmin;t<= tmax; t++)
    for(int r=tmin; r<= tmax; r++) {
      if(ONLY_FORWARD) Atr(t-tmin,r-tmin) = aE0_std(E0, PrecFloat(t+r),0);
      else  Atr(t-tmin,r-tmin) = aE0_std(E0, PrecFloat(t+r),0) + aE0_std(E0, PrecFloat(T - t +r), 0) + aE0_std(E0, PrecFloat(T+t -r), 0) + aE0_std(E0,PrecFloat(2*T -t -r), 0);
    }
  
  
  return;
}

void Get_Atr_std_Emax(PrecMatr& Atr, const PrecFloat &E0, int T, int tmin, int tmax)  {

  Atr.resize(tmax-tmin+1, tmax-tmin+1);

  for(int t=tmin;t<= tmax; t++)
    for(int r=tmin; r<= tmax; r++) {
      if(ONLY_FORWARD) Atr(t-tmin,r-tmin) = aE0_std_Emax(E0, PrecFloat(t+r),0) ;
      else Atr(t-tmin,r-tmin) = aE0_std_Emax(E0, PrecFloat(t+r),0) + aE0_std_Emax(E0, PrecFloat(T - t +r), 0) + aE0_std_Emax(E0, PrecFloat(T+t -r), 0) + aE0_std_Emax(E0,PrecFloat(2*T -t -r), 0);

    }
    
  
  return;
}



void Get_Rt(PrecVect& Rt, const PrecFloat &E0,  int T, int tmin, int tmax) {

  Rt.resize(tmax-tmin+1);

  

  for(int t=tmin;t<=tmax;t++) {
    PrecFloat x= exp(-E0*t)/t;
    PrecFloat y = exp(-E0*(T-t))/(T-t);
    Rt(t-tmin) =  x+y;
  }

  return;

}


void Get_ft(PrecVect& ft, const PrecFloat &E0, const PrecFloat &m, const PrecFloat &s, int jack_id, int T, int tmin, int tmax, string SMEARING_FUNC, const function<PrecFloat(const PrecFloat&, const PrecFloat&,const PrecFloat&,const PrecFloat&, int)> &f, const function<PrecFloat(const PrecFloat &, const PrecFloat &, const PrecFloat &, const PrecFloat &, int)> &F_NORM) {

  ft.resize(tmax-tmin+1);

  if(verbosity_lev>=2) cout<<"Entering get_ft!"<<endl<<flush;
 
  for(int t=tmin;t<=tmax;t++) {

    if(verbosity_lev>=2) cout<<"Computing time: "<<t<<"..."<<flush;
    
    PrecFloat t_n= t - PrecFloat(Beta);
    PrecFloat T_n= T-t -PrecFloat(Beta);
    
    const auto ftT=
      [&f, &m,&s,&E0, &t, &T, &jack_id, &F_NORM](const PrecFloat& x) -> PrecFloat
      {
	if(USE_GENERALIZED_NORM) return F_NORM(x,m,s,E0,jack_id)*f(x,m,s,E0, jack_id)*(exp(-x*(-PrecFloat(Beta)+t)) + ((ONLY_FORWARD)?PrecFloat(0):exp(-x*(-PrecFloat(Beta)+T-t)))) ;


	return exp(PrecFloat(Beta)*x)*f(x,m,s,E0,jack_id)*( exp(-x*t) + ((ONLY_FORWARD)?PrecFloat(0.0):exp(-x*(T-t))));
	//return f(x,m,s,E0,jack_id)*( exp(-x*(-PrecFloat(Beta)+PrecFloat(t))) + ((ONLY_FORWARD)?PrecFloat(0.0):exp(-x*(-PrecFloat(Beta)+PrecFloat(T-t)))));
      };

    if(Integrate_up_to_max_energy) {

      /*if(SMEARING_FUNC=="FF_Gauss_IM") {
	ft(t-tmin) = precPi()*(exp(t_n*(t_n*s*s -2*m)/2)*(erf((PrecFloat(Emax_int) -m+t_n*s*s)/(sqrt(PrecFloat(2))*s))  - erf((E0-m+t_n*s*s)/(sqrt(PrecFloat(2))*s)))/2
	+ ((ONLY_FORWARD)?PrecFloat(0):(exp(T_n*(T_n*s*s -2*m)/2)*(erf((PrecFloat(Emax_int) -m+T_n*s*s)/(sqrt(PrecFloat(2))*s))  - erf((E0-m+T_n*s*s)/(sqrt(PrecFloat(2))*s)))/2 ) )); } */
      //else
	ft(t-tmin) = integrateUpToXmax( ftT, E0.get(), Emax_int, (verbosity_lev > 1));
      
    }
    else {
      /*
      if(SMEARING_FUNC=="FF_Gauss_IM") { ft(t-tmin) = precPi()*( exp(t_n*(t_n*s*s -2*m)/2)*( 1 - erf((E0-m+t_n*s*s)/(sqrt(PrecFloat(2))*s)))/2    + ((ONLY_FORWARD)?PrecFloat(0):( exp(T_n*(T_n*s*s -2*m)/2)*( 1 - erf((E0-m+T_n*s*s)/(sqrt(PrecFloat(2))*s)))/2  )))  ;
	PrecFloat test= integrateUpToInfinite(ftT, E0.get(), (verbosity_lev > 1));
	cout.precision(100);
	cout<<"f("<<t<<"): EXACT: "<<ft(t-tmin)<<", NUM: "<<test<<endl<<flush;
	} */
      //else
      ft(t-tmin) =   integrateUpToInfinite(ftT, E0.get(), (verbosity_lev > 1));

    }
    
    if(verbosity_lev>=2) cout<<"done!"<<endl<<flush;
  }
 
  if(verbosity_lev>=2) cout<<"Exiting get_ft!"<<endl<<flush;
  
  return;
}

void Get_ft_std(PrecVect& ft, const PrecFloat &E0, const PrecFloat &m, const PrecFloat &s, int jack_id, int T, int tmin, int tmax, string SMEARING_FUNC, const function<PrecFloat(const PrecFloat&, const PrecFloat&,const PrecFloat&,const PrecFloat&, int)> &f) {

  ft.resize(tmax-tmin+1);

 
  for(int t=tmin;t<=tmax;t++) {

    PrecFloat tr = T-t;
   
    const auto ftT=
      [&f, &m,&s,&E0, &t, &T, &jack_id](const PrecFloat& x) -> PrecFloat
      {
    
	return f(x,m,s,E0, jack_id)*(exp(-x*t) + ((ONLY_FORWARD)?PrecFloat(0.0):exp(-x*(T-t))) ) ;
      };

    //if(SMEARING_FUNC=="FF_Gauss_IM")  ft(t-tmin) = precPi()*( exp(t*(t*s*s -2*m)/2)*(1 - erf((E0-m+t*s*s)/(sqrt(PrecFloat(2))*s)))/2  + ((ONLY_FORWARD)?PrecFloat(0):(exp(tr*(tr*s*s -2*m)/2)*(1 - erf((E0-m+tr*s*s)/(sqrt(PrecFloat(2))*s)))/2   )))  ;
    //else
    ft(t-tmin) =   integrateUpToInfinite(ftT, E0.get(), (verbosity_lev > 1));
  }


  return;
}

void Get_ft_std_Emax(PrecVect& ft, const PrecFloat &E0, const PrecFloat &m, const PrecFloat &s, int jack_id, int T, int tmin, int tmax, string SMEARING_FUNC, const function<PrecFloat(const PrecFloat&, const PrecFloat&,const PrecFloat&,const PrecFloat&, int)> &f) {

  ft.resize(tmax-tmin+1);

 
 
  for(int t=tmin;t<=tmax;t++) {

    PrecFloat tr= T-t;
   
    const auto ftT=
      [&f, &m,&s,&E0, &t, &T, &jack_id](const PrecFloat& x) -> PrecFloat
      {
    
	return f(x,m,s,E0, jack_id)*(exp(-x*t) + ((ONLY_FORWARD)?PrecFloat(0.0):exp(-x*(T-t)))  ) ;
      };

    /*if(SMEARING_FUNC=="FF_Gauss_IM") {
      ft(t-tmin) = precPi()*( exp(t*(t*s*s -2*m)/2)*(erf((PrecFloat(Emax_int) -m+t*s*s)/(sqrt(PrecFloat(2))*s)) - erf((E0-m+t*s*s)/(sqrt(PrecFloat(2))*s)))/2
      + ((ONLY_FORWARD)?PrecFloat(0):( exp(tr*(tr*s*s -2*m)/2)*(erf((PrecFloat(Emax_int) -m+tr*s*s)/(sqrt(PrecFloat(2))*s)) - erf((E0-m+tr*s*s)/(sqrt(PrecFloat(2))*s)))/2  )));   } */
    //else
    ft(t-tmin) = integrateUpToXmax( ftT, E0.get(), Emax_int, (verbosity_lev > 1));
  }
  

  return;
}




void Get_bt(PrecVect& bt,const PrecFloat &E,  int T, int tmin, int tmax) {

  bt.resize(tmax-tmin+1);

  for(int t=tmin;t<=tmax;t++) bt(t-tmin) = BaseFunc(E, t, T);

  return;


}

PrecFloat Get_norm_constraint(PrecFloat &m, PrecFloat &s, PrecFloat &E0, int jack_id, string SMEARING_FUNC, const function<PrecFloat(const PrecFloat&, const PrecFloat&,const PrecFloat&,const PrecFloat&,int)> &f) {

  PrecFloat res=1;

  
  const auto f1=
    [&f, &m,&s,&E0, &jack_id](const PrecFloat& x) -> PrecFloat
    {
      return f(x, m, s, E0, jack_id) ;
    };

  if(Integrate_up_to_max_energy) return integrateUpToXmax(f1, E0.get(), Emax_int, (verbosity_lev > 1));
  else return  integrateUpToInfinite(f1, E0.get(), (verbosity_lev > 1));
  
}


PrecFloat Get_M2(PrecFloat &m, PrecFloat &s, PrecFloat &E0, int jack_id,  const function<PrecFloat(const PrecFloat&, const PrecFloat&,const PrecFloat&,const PrecFloat&, int)> &f, const function<PrecFloat(const PrecFloat &, const PrecFloat &, const PrecFloat &, const PrecFloat &, int)> &F_NORM) {


  const auto f2=
    [&f, &m,&s,&E0, &jack_id, &F_NORM](const PrecFloat& x) -> PrecFloat
    {
      if(USE_GENERALIZED_NORM) return F_NORM(x,m,s,E0,jack_id)*exp(PrecFloat(Beta)*x)*pow(f(x, m, s, E0, jack_id),2);
      return exp(PrecFloat(Beta)*x)*pow(f(x, m, s, E0, jack_id),2);
    };

  if(USE_GENERALIZED_NORM) {
    if(Integrate_up_to_max_energy) return integrateUpToXmax(f2, 0.0, Emax_int, (verbosity_lev > 1));
    else return integrateUpToInfinite(f2, 0.0, (verbosity_lev > 1));
  }
  else{
    if(Integrate_up_to_max_energy) return integrateUpToXmax(f2, E0.get(), Emax_int, (verbosity_lev > 1));
    else return  integrateUpToInfinite(f2, E0.get(), (verbosity_lev > 1));
  }

  return PrecFloat(0.0);

}

PrecFloat Get_M2_std_norm(PrecFloat &m, PrecFloat &s, PrecFloat &E0, int jack_id,  const function<PrecFloat(const PrecFloat&, const PrecFloat&,const PrecFloat&,const PrecFloat&, int)> &f) {


  const auto f2=
    [&f, &m,&s,&E0, &jack_id](const PrecFloat& x) -> PrecFloat
    {
      return pow(f(x, m, s, E0, jack_id),2);
    };

  return  integrateUpToInfinite(f2, E0.get(), (verbosity_lev > 1));

}

PrecFloat Get_M2_std_norm_Emax(PrecFloat &m, PrecFloat &s, PrecFloat &E0, int jack_id,  const function<PrecFloat(const PrecFloat&, const PrecFloat&,const PrecFloat&,const PrecFloat&, int)> &f) {


  const auto f2=
    [&f, &m,&s,&E0, &jack_id](const PrecFloat& x) -> PrecFloat
    {
      return pow(f(x, m, s, E0, jack_id),2);
    };

  return integrateUpToXmax(f2, E0.get(), Emax_int, (verbosity_lev > 1));
  integrateUpToInfinite(f2, E0.get(), (verbosity_lev > 1));

}

void Get_M_N(PrecFloat &m, PrecFloat &s, PrecFloat &E0, int jack_id,  const function<PrecFloat(const PrecFloat&, const PrecFloat&,const PrecFloat&,const PrecFloat&, int)> &f,  PrecVect &M_n) {

  if(Nmoms==0) return;
  
  M_n.resize(Nmoms);

  for(int n=0;n<Nmoms;n++) {

    const auto fn=
      [&f, &m,&s,&E0, &n, &jack_id](const PrecFloat& x) -> PrecFloat
      {
	return f(x, m, s, E0, jack_id)*pow(x,PrecFloat(n)) ;
      };

    if(Integrate_up_to_max_energy) M_n(n) = integrateUpToXmax(fn, E0.get(), Emax_int, (verbosity_lev > 1));
    else M_n(n) = integrateUpToInfinite(fn,E0.get(), (verbosity_lev > 1));   
     
  }

  return;

}

void Get_M_tilde_N(const PrecVect &ft, const PrecMatr &Atr_inv, vector<PrecVect> &Rt_n,  PrecVect &M_tilde_n) {

  if(Nmoms==0) return;

  M_tilde_n.resize(Nmoms);

  for(int n=0;n<Nmoms;n++) {

    M_tilde_n(n) = Rt_n[n].transpose()*Atr_inv*ft;   

  }

  return;

}


void Get_Rt_up_to_N(PrecFloat &E0, int T, int tmin, int tmax, vector<PrecVect> &Rt_n) {

  if(Nmoms==0) return;

  Rt_n.resize(Nmoms);


  for(int n=0;n<Nmoms;n++) {

    Rt_n[n].resize(tmax-tmin+1);

    for(int t=tmin;t<=tmax;t++) Rt_n[n](t-tmin) = aE0(E0,PrecFloat(t),n) + aE0(E0,PrecFloat(T-t),n);

  }


  return;
 
}

void Get_G_matrix(PrecMatr &G,const PrecMatr &Atr_inv, vector<PrecVect> &Rt_n) {

  G.resize(Nmoms,Nmoms);

  for(int i=0;i<Nmoms;i++) {
    for(int j=0;j<Nmoms;j++) {

      G(i,j) = Rt_n[i].transpose()*Atr_inv*Rt_n[j];

    }
  }

  return;
  
}

void Compute_covariance_matrix(PrecMatr &B,const PrecMatr &Atr, const distr_t_list &corr, int tmin, int tmax, PrecFloat m, PrecFloat s,  string analysis_name,  string MODE, Vfloat &covariance) {

  
  PrecFloat norm= sqr(PrecFloat(corr.ave(0)));

  
  if(analysis_name.substr(0,7)=="R_ratio") norm = sqr(PrecFloat(corr.ave(1)))*pow(m,-6);
  if(analysis_name.substr(0,10)=="virtual_FF") norm= sqr(PrecFloat(corr.ave(1)))*pow(m+PrecFloat(0.08),-2);

  B.resize(tmax-tmin+1, tmax-tmin+1);


  for(int t=tmin;t<=tmax;t++) {
    for(int r=t;r<=tmax;r++) {
      if(MODE=="TANT") {
	B(t-tmin,r-tmin) = covariance[t*corr.size()+ r];
	B(r-tmin, t-tmin) = B(t-tmin,r-tmin);
      }
      else {
	B(t-tmin, r-tmin) = covariance[t*corr.size()+r]/(corr.ave(t)*corr.ave(r));
	B(r-tmin,t-tmin) = B(t-tmin,r-tmin);
      }
    }
  }

  /*
  for(int t=tmin;t<=tmax;t++)
    for(int r=t;r<=tmax;r++) if ( fabs( ((B(t-tmin,r-tmin) - B(r-tmin,t-tmin))/(0.5*( B(t-tmin,r-tmin)+B(r-tmin,t-tmin)))).get()) > 1e-14 ) crash("covariance matrix is not symmetric, Cov(i,j) = "+to_string_with_precision(B(t-tmin,r-tmin).get(),50)+" != Cov(j,i)= "+to_string_with_precision(B(r-tmin,t-tmin).get(),50));
  */

   
  if(MODE != "TANT") B= Atr*B;
  else B= (1.0/norm)*B;

 
      
 
  return;
}


void Get_optimal_lambda(const PrecMatr &Atr, const PrecMatr &Atr_std_norm, const PrecMatr &Atr_std_norm_Emax, const PrecMatr &B,const PrecVect &ft, vector<PrecVect>& ft_jack, const PrecVect &ft_std_norm, const PrecVect &ft_std_norm_Emax, const PrecFloat & M2, const PrecFloat &M2_std_norm, const PrecFloat &M2_std_norm_Emax, const double &mean, const double &sigma, const double &Estart,  double& lambda_opt, double& lambda_opt_10, vector<PrecVect> Rt_n, const PrecVect &M_n , vector<PrecVect>& M_n_jack, const distr_t_list & corr,int T, int tmin, int tmax,const double mult,  string MODE, string curr_type, string SMEARING_FUNC, string CORR_NAME, double Ag_ov_A0_tg, bool JackOnKer,const distr_t& Prefact, const distr_t &offset,  string analysis_name,  const function<PrecFloat(const PrecFloat&, const PrecFloat&,const PrecFloat&,const PrecFloat&, int)> &f,  const function<double(const function<double(double)>&)> &syst_func, bool Use_guess_density, const function<double(double)> &guess_density ) {


 

  int MAX_Iters = 1000;

  int Global_id=0;

  string Emax_val= (Integrate_up_to_max_energy==0)?"inf":to_string_with_precision(Emax_int,1);

 
  string out_path = MODE+"_"+CORR_NAME+"_"+curr_type+"_"+SMEARING_FUNC+"_E*_"+to_string_with_precision(mean,3)+"_sigma_"+to_string_with_precision(sigma,3)+"_E0_"+to_string_with_precision(Estart,3)+"_T_"+to_string(T)+"_tmin_"+to_string(tmin)+"_tmax_"+to_string(tmax)+"_alpha_"+to_string(alpha)+".dat";
  
  ofstream Print_R_at_lambda("../data/spectral_reconstruction/"+analysis_name+"/lambda_stability/beta_"+to_string_with_precision(Beta,2)+"_Emax_"+Emax_val+"/"+out_path);

  //Print header
  if(Use_guess_density) Print_R_at_lambda<<"# $1=lambda,  $2= A[g]/A[0], $3= A^E[g]/A^E[0] , $4=B[g], $5=B^E[g],  $6=Aa[g]/Aa[0], $7= val, $8=err, $9=syst, $10=S_FLAG $11=mult"<<endl; 
  else Print_R_at_lambda<<"# $1=lambda,  $2= A[g]/A[0], $3= A^E[g]/A^E[0] , $4=B[g], $5=B^E[g], $6=Aa[g]/Aa[0], $7= val, $8=err, $9=syst,  $10=S_FLAG $11=mult"<<endl;

  Print_R_at_lambda.precision(10);
    
  if(verbosity_lev) cout<<"Finding optimal lambda* ..."<<endl;

  int Njacks= corr.distr_list[0].distr.size();

 
   
  const auto A1=
    [&M2, &Atr, &ft](const PrecVect& gmin) -> PrecFloat
    {
      PrecVect Atr_g = Atr*gmin;
      PrecFloat g_Atr_g = gmin.transpose()*Atr_g;
      PrecFloat ft_g = ft.transpose()*gmin;
      return (M2 + g_Atr_g -2*ft_g)/M2  ;
    };

  const auto B1=
    [&B, &M2, &MODE](const PrecVect& gmin) -> PrecFloat
    {
      PrecVect B_g = B*gmin;
      PrecFloat g_B_g = gmin.transpose()*B_g;
      if(MODE=="SANF") g_B_g = g_B_g/M2; 
      return g_B_g ;
    };

  const auto A1_std=
    [&M2_std_norm, &Atr_std_norm, &ft_std_norm](const PrecVect& gmin) -> PrecFloat
    {
      PrecVect Atr_g = Atr_std_norm*gmin;
      PrecFloat g_Atr_g = gmin.transpose()*Atr_g;
      PrecFloat ft_g = ft_std_norm.transpose()*gmin;
      return (M2_std_norm + g_Atr_g -2*ft_g)/M2_std_norm  ;
    };

  const auto B1_std=
    [&B, &M2_std_norm, &MODE](const PrecVect& gmin) -> PrecFloat
    {
      PrecVect B_g = B*gmin;
      PrecFloat g_B_g = gmin.transpose()*B_g;
      if(MODE=="SANF") g_B_g = g_B_g/M2_std_norm; 
      return g_B_g ;
    };

  const auto A1_std_Emax=
    [&M2_std_norm_Emax, &Atr_std_norm_Emax, &ft_std_norm_Emax](const PrecVect& gmin) -> PrecFloat
    {
      PrecVect Atr_g = Atr_std_norm_Emax*gmin;
      PrecFloat g_Atr_g = gmin.transpose()*Atr_g;
      PrecFloat ft_g = ft_std_norm_Emax.transpose()*gmin;
      return (M2_std_norm_Emax + g_Atr_g -2*ft_g)/M2_std_norm_Emax  ;
    };

  const auto B1_std_Emax=
    [&B, &M2_std_norm_Emax, &MODE](const PrecVect& gmin) -> PrecFloat
    {
      PrecVect B_g = B*gmin;
      PrecFloat g_B_g = gmin.transpose()*B_g;
      if(MODE=="SANF") g_B_g = g_B_g/M2_std_norm_Emax; 
      return g_B_g ;
    };


  //bisection search
  PrecFloat l_low =0;
  PrecFloat l_up = 1;
  PrecFloat l_start=1.0;
  int Nit=0;
  int Nit_Ag0=0;
  int Nit_10=0;
  int Nit_100=0;
  double lambda_balance;
  double lambda_balance_10;
  double lambda_balance_100;
  bool lambda_balance_found=false;
  bool lambda_balance_found_10=false;
  PrecFloat Ag_ov_A0_target=1e-4;
  if(Ag_ov_A0_tg > 0) Ag_ov_A0_target= PrecFloat(Ag_ov_A0_tg);
  Vfloat Ags_mult ={10.0, 1.0}; 
 
  if(verbosity_lev) cout<<"Finding lambdas corresponding to A[g]/A[0] in  "<<Ag_ov_A0_target<<" * { 1, 10}"<<endl;

  //###############################################################################################
  

  
  






  //################################################################################################
  
  int counter_Ag_m=0;
 
  for(auto & Ag_m: Ags_mult) {
  
    bool lambda_found_Ag_A0=false;
    Nit_Ag0=0;
    l_low=0;
    
    //bisection search for given A[g]/A[0]
    while( !lambda_found_Ag_A0 ) {

     
      PrecFloat lambda_mid  =  (Nit_Ag0==0 && counter_Ag_m==0)?l_start:(l_up+l_low)/2;
      PrecMatr C = Atr*(1-lambda_mid)/M2 + B*lambda_mid/(MODE=="SANF"?M2:1.0);
      PrecMatr C_inv = C.inverse();
      PrecVect ft_l = ft*(1-lambda_mid)/M2;
      PrecVect M_tilde_n;
      PrecMatr G_n;
      vector<PrecVect> ft_l_jack;
      if(JackOnKer) for(int ijack=0;ijack<Njacks;ijack++) ft_l_jack.push_back( ft_jack[ijack]*(1-lambda_mid)/M2);


    


       
      //Lagrangian multipliers
    
      Get_M_tilde_N(ft_l, C_inv, Rt_n, M_tilde_n);
      Get_G_matrix(G_n, C_inv, Rt_n);
    
         

      PrecVect gm = C_inv*ft_l;
      vector<PrecVect> gm_jack;
      PrecVect pn;
      PrecMatr G_n_inv;
      PrecVect M_n_diff;
      if(Nmoms > 0) {
	G_n_inv= G_n.inverse();
	M_n_diff = M_n - M_tilde_n;
	pn = G_n_inv*M_n_diff;
      }

      for(int n=0;n<Nmoms;n++) {
	PrecVect lmult_n = C_inv*Rt_n[n]; 
	gm = gm + lmult_n*pn(n);
      }

      if(JackOnKer) {

	for(int ijack=0; ijack<Njacks;ijack++) {
	  PrecVect gm_ij = C_inv*ft_l_jack[ijack];
	  PrecVect M_tilde_n_ij;
	  Get_M_tilde_N(ft_l_jack[ijack], C_inv, Rt_n, M_tilde_n_ij);
	  PrecVect M_n_diff_ij, pn_ij;
	  if(Nmoms > 0) {
	    M_n_diff_ij= M_n_jack[ijack] - M_tilde_n_ij;
	    pn_ij = G_n_inv*M_n_diff_ij;
	  }
	  for(int n=0; n<Nmoms;n++) {
	    PrecVect lmult_n= C_inv*Rt_n[n];
	    gm_ij = gm_ij + lmult_n*pn_ij(n);
	  }
	  gm_jack.push_back(gm_ij);
	}

      }

   
      PrecFloat A1_val = A1(gm);
      PrecFloat B1_val = B1(gm);
      PrecFloat W_val = (1-lambda_mid)*A1_val + lambda_mid*B1_val;

      PrecFloat A1_val_std = A1_std(gm);
      PrecFloat B1_val_std = B1_std(gm);
      PrecFloat W_val_std = (1-lambda_mid)*A1_val_std + lambda_mid*B1_val_std;

      PrecFloat A1_val_std_Emax= A1_std_Emax(gm);
      PrecFloat B1_val_std_Emax= B1_std_Emax(gm);

      double syst=0.0;

   
      if(Use_guess_density) {

	auto integrand_syst = [&tmin, &tmax, &T, &gm, &f, &mean, &sigma, &Estart](double E) ->double {
	    PrecVect bt;
	    Get_bt(bt, E, T, tmin,tmax);
	    return (f(E,mean, sigma, Estart, -1) - gm.transpose()*bt).get();
	  };

	if(extended_analysis_of_syst)   syst= syst_func(integrand_syst);
	

	if(print_reco_in_stability_analysis) {
	  int Npoints;
	  if(analysis_name == "tau_decay") {Npoints= 2000; step_size=0.005;}
	  else Npoints= (int)(((mean+20*sigma -Estart)/(sigma*step_size)));
	  Vfloat Erg, Error, Spec_dens_guess, Agvs;

	  for(int ip=0; ip<Npoints;ip++) {
	    
	    double E;
	    Agvs.push_back(A1_val_std_Emax.get());
	    if(analysis_name == "tau_decay") E= (Estart+ ip*step_size);
	    else E = (Estart + (ip*sigma)*step_size);
	    Error.push_back(integrand_syst(E));
	    Erg.push_back(E);
	    Spec_dens_guess.push_back(guess_density(E));
	  }

	  Print_To_File({}, {Agvs, Erg, Error, Spec_dens_guess}, "../data/spectral_reconstruction/"+analysis_name+"/error_funcs/beta_"+to_string_with_precision(Beta,2)+"/"+to_string(Global_id)+"_"+out_path,"", "#E   |reco-exact|   guess ");


	}
	
	
      }

      cout.precision(20);

         
     


      //##########################################################################################
      //compute anti-Laplace transform corresponding to lambda_mid:
      distr_t R_E_lambda;
      
      for(int ijack=0; ijack< corr.distr_list[1].size(); ijack++) {

	PrecFloat spec_lambda_d_jack=0;

	for(int t=tmin;t<=tmax;t++) spec_lambda_d_jack = spec_lambda_d_jack + ((JackOnKer)?gm_jack[ijack](t-tmin):gm(t-tmin))*corr.distr_list[t].distr[ijack]; 

	R_E_lambda.distr.push_back( spec_lambda_d_jack.get());
      }

      PrecFloat mult_est = (mult_estimated_from_norm0==false)?A1_val/B1_val:A1_val_std/B1_val_std;
     
      Print_R_at_lambda<<lambda_mid<<" "<<A1_val_std<<" "<<A1_val_std_Emax<<" "<<B1_val_std<<" "<<B1_val_std_Emax<<" "<<A1_val<<" "<<(Prefact*(R_E_lambda)+offset).ave()<<" "<<(Prefact*(R_E_lambda)+offset).err();
      Print_R_at_lambda<<" "<<syst;
      Print_R_at_lambda<<" 0 "<<mult_est<<endl;

      //##########################################################################################

     
   

      if(A1_val_std > Ag_ov_A0_target*Ag_m) { // lambda_mid is new l_low
	l_up =lambda_mid;
      }
      else { //lambda_mid is new l_up
	l_low = lambda_mid;
      }

      
      
      Nit_Ag0++;
      if( (A1_val_std < 1.5*Ag_ov_A0_target*Ag_m && A1_val_std > 0.7*Ag_ov_A0_target*Ag_m)) lambda_found_Ag_A0=true;

      if(Nit_Ag0 >= 60) { //R_ratio used Nit_Ag0=10
	//cout<<"WARNING: A[g]/A[0]: "<<Ag_ov_A0_target*Ag_m<<" cannot be obtained after "<<Nit_Ag0<<" iterations...Skipping!"<<endl;
	lambda_found_Ag_A0=true;
      }
      Global_id++;
    }
   
    counter_Ag_m++;
  }

  if(verbosity_lev) cout<<"A[g]/A[0] scan completed!"<<endl;
  cout.precision(10);
 


  //#############################################################################################################################################
  //#############################################################################################################################################
  //#############################################################################################################################################
  //#############################################################################################################################################
  
  if(verbosity_lev) cout<<"Finding lambda from balance condition A=mult*B, mult= "<<mult<<endl;
  l_up=1.0;
  l_low =0.0;
  
 
  //bisection search for condition A = mult*B
  while( !lambda_balance_found ) {

   
    //evaluate the minimum at midpoint
    PrecFloat lambda_mid = (Nit==0)?l_start:(l_up+l_low)/2;
    PrecMatr C = Atr*(1-lambda_mid)/M2 + B*lambda_mid/(MODE=="SANF"?M2:1);
    PrecMatr C_inv = C.inverse();
    PrecVect ft_l = ft*(1-lambda_mid)/M2;
    vector<PrecVect> ft_l_jack;
    if(JackOnKer) for(int ijack=0;ijack<Njacks;ijack++) ft_l_jack.push_back( ft_jack[ijack]*(1-lambda_mid)/M2);
    PrecVect M_tilde_n;
    PrecMatr G_n;

    //Lagrangian multipliers
    
    Get_M_tilde_N(ft_l, C_inv, Rt_n, M_tilde_n);
    Get_G_matrix(G_n, C_inv, Rt_n);
    
      
    
    PrecVect gm = C_inv*ft_l;

    vector<PrecVect> gm_jack;

    PrecVect pn, M_n_diff;

    PrecMatr G_n_inv;


    if(Nmoms > 0) {
      G_n_inv = G_n.inverse();
      M_n_diff = M_n - M_tilde_n;
      pn = G_n_inv*M_n_diff;
    }

    for(int n=0;n<Nmoms;n++) {
      PrecVect lmult_n = C_inv*Rt_n[n]; 
      gm = gm + lmult_n*pn(n);
    }


    if(JackOnKer) {

      for(int ijack=0; ijack<Njacks;ijack++) {
	PrecVect gm_ij = C_inv*ft_l_jack[ijack];
	PrecVect M_tilde_n_ij, M_n_diff_ij, pn_ij;
	Get_M_tilde_N(ft_l_jack[ijack], C_inv, Rt_n, M_tilde_n_ij);
	if(Nmoms>0) {
	  PrecVect M_n_diff_ij = M_n_jack[ijack] - M_tilde_n_ij;
	  PrecVect pn_ij = G_n_inv*M_n_diff_ij;
	}
	for(int n=0; n<Nmoms;n++) {
	  PrecVect lmult_n= C_inv*Rt_n[n];
	  gm_ij = gm_ij + lmult_n*pn_ij(n);
	}
	gm_jack.push_back(gm_ij);
      }
      
    }

    PrecFloat A1_val = A1(gm);
    PrecFloat B1_val = B1(gm);
    PrecFloat W_val = (1-lambda_mid)*A1_val + lambda_mid*B1_val;
    
    PrecFloat A1_val_std = A1_std(gm);
    PrecFloat B1_val_std = B1_std(gm);
    PrecFloat W_val_std = (1-lambda_mid)*A1_val_std + lambda_mid*B1_val_std;


    PrecFloat A1_val_std_Emax= A1_std_Emax(gm);
    PrecFloat B1_val_std_Emax= B1_std_Emax(gm);


  


    
      
    double mult_est = (mult_estimated_from_norm0==false)?(A1_val/B1_val).get():(A1_val_std/B1_val_std).get();
         
    if(mult> mult_est) { // lambda_mid is new l_low
      l_low =lambda_mid;
    }
    else { //lambda_mid is new l_up
      l_up = lambda_mid;
    }
    
    lambda_balance= lambda_mid.get();
    
    Nit++;
    if(fabs(mult_est - mult)/(mult) < 0.01) lambda_balance_found=true;



    double syst=0.0;

    if(Use_guess_density) {

      auto integrand_syst = [&tmin, &tmax, &T, &gm, &f, &mean, &sigma, &Estart](double E) ->double {
	  PrecVect bt;
	  Get_bt(bt, E, T, tmin,tmax);
	  return (f(E,mean, sigma, Estart, -1) - gm.transpose()*bt).get();
	};

      if(extended_analysis_of_syst || lambda_balance_found)  syst= syst_func(integrand_syst);


      if(print_reco_in_stability_analysis) {
	int Npoints;
	if(analysis_name == "tau_decay") {Npoints= 2000; step_size=0.005;}
	else Npoints= (int)(((mean+20*sigma -Estart)/(sigma*step_size)));
	Vfloat Erg, Error, Spec_dens_guess, Agvs;
	
	for(int ip=0; ip<Npoints;ip++) {
	  
	  double E;
	  Agvs.push_back(A1_val_std.get());
	  if(analysis_name == "tau_decay") E= (Estart+ ip*step_size);
	  else E = (Estart + (ip*sigma)*step_size);
	  Error.push_back(integrand_syst(E));
	  Erg.push_back(E);
	  Spec_dens_guess.push_back(guess_density(E));
	}
	
	Print_To_File({}, {Agvs,Erg, Error, Spec_dens_guess}, "../data/spectral_reconstruction/"+analysis_name+"/error_funcs/beta_"+to_string_with_precision(Beta,2)+"/"+to_string(Global_id)+"_"+out_path,"", "#E   |reco-exact|   guess ");

	
      }
    }

    
    if(Nit > MAX_Iters) {
      cout<<"###### FAILED CONVERGENCE #########"<<endl;
      cout<<"###### INFO #######################"<<endl;
      cout<<"-----------------------------"<<endl;
      cout<<"Current values of g[t] & f[t]:"<<endl;
      for(int t=tmin;t<=tmax;t++) cout<<"t: "<<t<<"   "<<gm(t-tmin).get()<<"      "<<ft(t-tmin)<<endl;
      cout<<"-----------------------------"<<endl;
      cout<<"M2: "<<M2.get()<<endl;
      cout<<"M2(std. norm): "<<M2_std_norm.get()<<endl;
      cout<<"M2(std. norm Emax): "<<M2_std_norm_Emax.get()<<endl;
      cout<<"A[g]: "<<A1_val<<", B[g]: "<<B1_val<<endl;
      cout<<"Printing inverse W_tr = A_tr*(1-lambda) + B_tr*lambda matrix: "<<endl;
      cout<<"-------------------------------------------------------------"<<endl;
      cout<<C_inv<<endl;
      cout<<"-------------------------------------------------------------"<<endl;
      cout<<"Printing inverse A_tr matrix: "<<endl;
      cout<<"-------------------------------------------------------------"<<endl;
      cout<<Atr.inverse()<<endl;
      cout<<"-------------------------------------------------------------"<<endl;
      cout<<"lambda_low: "<<l_low<<" lambda_up: "<<l_up<<endl;
      cout<<"####################################"<<endl;
      crash("After "+to_string(Nit)+" iterations, balance condition A = mult*B cannot be obtained for CORR: "+CORR_NAME+" , MODE: "+MODE+", CURR_TYPE: "+curr_type+", SM_TYPE: "+SMEARING_FUNC+", Beta: "+to_string_with_precision(Beta,3)+", mult(target) = "+to_string_with_precision( mult, 8)+", a*sigma: "+to_string_with_precision(sigma,5)+", aE*: "+to_string_with_precision(mean, 5)+" lambda: "+to_string_with_precision(lambda_mid, 5)+" , mult: "+to_string_with_precision(mult_est,5));
    
      
    }


    //##########################################################################################
    //compute anti-Laplace transform corresponding to lambda_mid:
    distr_t R_E_lambda;
      
    for(int ijack=0; ijack< corr.distr_list[1].size(); ijack++) {

      PrecFloat spec_lambda_d_jack=0;

      for(int t=tmin;t<=tmax;t++) spec_lambda_d_jack = spec_lambda_d_jack + ((JackOnKer==true)?gm_jack[ijack](t-tmin):gm(t-tmin))*corr.distr_list[t].distr[ijack]; 

      R_E_lambda.distr.push_back( spec_lambda_d_jack.get());
    }

      

    Print_R_at_lambda<<lambda_mid<<" "<<A1_val_std<<" "<<A1_val_std_Emax<<" "<<B1_val_std<<" "<<B1_val_std_Emax<<" "<<A1_val<<" "<<(Prefact*(R_E_lambda)+offset).ave()<<" "<<(Prefact*(R_E_lambda)+offset).err();
    Print_R_at_lambda<<" "<<syst;
    Print_R_at_lambda<<" "<<lambda_balance_found<<" "<<mult_est<<endl;

    //##########################################################################################

     
    
    Global_id++;
  }


  if(verbosity_lev) cout<<"lambda_opt = "<<lambda_balance<<endl;

  //#############################################################################################################################################
  //#############################################################################################################################################
  //#############################################################################################################################################
  //#############################################################################################################################################







  //#############################################################################################################################################
  //#############################################################################################################################################
  //#############################################################################################################################################
  //#############################################################################################################################################

  double k=0.1;
  l_low =0.0;
  if(verbosity_lev) cout<<"Finding lambda from balance condition A=(k="<<to_string_with_precision(k,2)<<")*mult*B, mult= "<<mult<<endl;
 
 
  //bisection search for condition A = mult*B
  while( !lambda_balance_found_10 ) {

   

    //evaluate the minimum at midpoint
    PrecFloat lambda_mid = (Nit_10==0)?l_up:(l_up+l_low)/2;
    PrecMatr C = Atr*(1-lambda_mid)/M2 + B*lambda_mid/(MODE=="SANF"?M2:1);
    PrecMatr C_inv = C.inverse();
    PrecVect ft_l = ft*(1-lambda_mid)/M2;
    vector<PrecVect> ft_l_jack;
    if(JackOnKer) for(int ijack=0;ijack<Njacks;ijack++) ft_l_jack.push_back( ft_jack[ijack]*(1-lambda_mid)/M2);
    PrecVect M_tilde_n;
    PrecMatr G_n;

    //Lagrangian multipliers
    
    Get_M_tilde_N(ft_l, C_inv, Rt_n, M_tilde_n);
    Get_G_matrix(G_n, C_inv, Rt_n);
    
      
    
    PrecVect gm = C_inv*ft_l;

    vector<PrecVect> gm_jack;

    PrecVect pn, M_n_diff;

    PrecMatr G_n_inv;


    if(Nmoms > 0) {
      G_n_inv = G_n.inverse();
      M_n_diff = M_n - M_tilde_n;
      pn = G_n_inv*M_n_diff;
    }

    for(int n=0;n<Nmoms;n++) {
      PrecVect lmult_n = C_inv*Rt_n[n]; 
      gm = gm + lmult_n*pn(n);
    }


    if(JackOnKer) {

      for(int ijack=0; ijack<Njacks;ijack++) {
	PrecVect gm_ij = C_inv*ft_l_jack[ijack];
	PrecVect M_tilde_n_ij, M_n_diff_ij, pn_ij;
	Get_M_tilde_N(ft_l_jack[ijack], C_inv, Rt_n, M_tilde_n_ij);
	if(Nmoms>0) {
	  PrecVect M_n_diff_ij = M_n_jack[ijack] - M_tilde_n_ij;
	  PrecVect pn_ij = G_n_inv*M_n_diff_ij;
	}
	for(int n=0; n<Nmoms;n++) {
	  PrecVect lmult_n= C_inv*Rt_n[n];
	  gm_ij = gm_ij + lmult_n*pn_ij(n);
	}
	gm_jack.push_back(gm_ij);
      }
      
    }

    PrecFloat A1_val = A1(gm);
    PrecFloat B1_val = B1(gm);
    PrecFloat W_val = (1-lambda_mid)*A1_val + lambda_mid*B1_val;
    
    PrecFloat A1_val_std = A1_std(gm);
    PrecFloat B1_val_std = B1_std(gm);
    PrecFloat W_val_std = (1-lambda_mid)*A1_val_std + lambda_mid*B1_val_std;

    PrecFloat A1_val_std_Emax = A1_std_Emax(gm);
    PrecFloat B1_val_std_Emax = B1_std_Emax(gm);


  

      
    double mult_est = (mult_estimated_from_norm0==false)?(A1_val/B1_val).get():(A1_val_std/B1_val_std).get();
         
    if(k*mult> mult_est) { // lambda_mid is new l_low
      l_low =lambda_mid;
    }
    else { //lambda_mid is new l_up
      l_up = lambda_mid;
    }

    lambda_balance_10= lambda_mid.get();
      
    Nit_10++;
    if(fabs(mult_est - k*mult)/(k*mult) < 0.01) lambda_balance_found_10=true;


    double syst=0.0;

    if(Use_guess_density) {

      auto integrand_syst = [&tmin, &tmax, &T, &gm, &f, &mean, &sigma, &Estart](double E) ->double {
	  PrecVect bt;
	  Get_bt(bt, E, T, tmin,tmax);
	  return (f(E,mean, sigma, Estart, -1) - gm.transpose()*bt).get();
	};

      if(extended_analysis_of_syst || lambda_balance_found_10) 	syst= syst_func(integrand_syst);
      

	
      if(print_reco_in_stability_analysis) {
	int Npoints;
	if(analysis_name == "tau_decay") {Npoints= 2000; step_size=0.005;}
	else Npoints= (int)(((mean+20*sigma -Estart)/(sigma*step_size)));
	Vfloat Erg, Error, Spec_dens_guess, Agvs;
	  
	for(int ip=0; ip<Npoints;ip++) {
	    
	  double E;
	  Agvs.push_back(A1_val_std.get());
	  if(analysis_name == "tau_decay") E= (Estart+ ip*step_size);
	  else E = (Estart + (ip*sigma)*step_size);
	  Error.push_back(integrand_syst(E));
	  Erg.push_back(E);
	  Spec_dens_guess.push_back(guess_density(E));
	}
	  
	Print_To_File({}, {Agvs, Erg, Error, Spec_dens_guess}, "../data/spectral_reconstruction/"+analysis_name+"/error_funcs/beta_"+to_string_with_precision(Beta,2)+"/"+to_string(Global_id)+"_"+out_path,"", "#E   |reco-exact|   guess ");
	  
	  
      }
    }


      

    if(Nit_10 > MAX_Iters) {
      cout<<"###### FAILED CONVERGENCE #########"<<endl;
      cout<<"###### INFO #######################"<<endl;
      cout<<"-----------------------------"<<endl;
      cout<<"Current values of g[t] & f[t]:"<<endl;
      for(int t=tmin;t<=tmax;t++) cout<<"t: "<<t<<"   "<<gm(t-tmin).get()<<"      "<<ft(t-tmin)<<endl;
      cout<<"-----------------------------"<<endl;
      cout<<"M2: "<<M2.get()<<endl;
      cout<<"M2(std. norm): "<<M2_std_norm.get()<<endl;
      cout<<"M2(std. norm Emax): "<<M2_std_norm_Emax.get()<<endl;
      cout<<"A[g]: "<<A1_val<<", B[g]: "<<B1_val<<endl;
      cout<<"Printing inverse W_tr = A_tr*(1-lambda) + B_tr*lambda matrix: "<<endl;
      cout<<"-------------------------------------------------------------"<<endl;
      cout<<C_inv<<endl;
      cout<<"-------------------------------------------------------------"<<endl;
      cout<<"Printing inverse A_tr matrix: "<<endl;
      cout<<"-------------------------------------------------------------"<<endl;
      cout<<Atr.inverse()<<endl;
      cout<<"-------------------------------------------------------------"<<endl;
      cout<<"lambda_low: "<<l_low<<" lambda_up: "<<l_up<<endl;
      cout<<"####################################"<<endl;
      crash("After "+to_string(Nit_10)+" iterations, balance condition A = k*mult*B cannot be obtained for CORR: "+CORR_NAME+" , MODE: "+MODE+", CURR_TYPE: "+curr_type+", SM_TYPE: "+SMEARING_FUNC+", Beta: "+to_string_with_precision(Beta,3)+", mult(target) = "+to_string_with_precision( k*mult, 8)+" a*sigma: "+to_string_with_precision(sigma,5)+", aE*: "+to_string_with_precision(mean, 5)+" lambda: "+to_string_with_precision(lambda_mid, 5)+" , mult: "+to_string_with_precision(mult_est,5));
    }


    //##########################################################################################
    //compute anti-Laplace transform corresponding to lambda_mid:
    distr_t R_E_lambda;
      
    for(int ijack=0; ijack< corr.distr_list[1].size(); ijack++) {

      PrecFloat spec_lambda_d_jack=0;

      for(int t=tmin;t<=tmax;t++) spec_lambda_d_jack = spec_lambda_d_jack + ((JackOnKer==true)?gm_jack[ijack](t-tmin):gm(t-tmin))*corr.distr_list[t].distr[ijack]; 

      R_E_lambda.distr.push_back( spec_lambda_d_jack.get());
    }

    Print_R_at_lambda<<lambda_mid<<" "<<A1_val_std<<" "<<A1_val_std_Emax<<" "<<B1_val_std<<" "<<B1_val_std_Emax<<" "<<A1_val<<" "<<(Prefact*(R_E_lambda)+offset).ave()<<" "<<(Prefact*(R_E_lambda)+offset).err();
    Print_R_at_lambda<<" "<<syst;
    Print_R_at_lambda<<" "<<2*(lambda_balance_found_10==1)<<" "<<mult_est<<endl;

    //##########################################################################################

    Global_id++;
    
   
  }


  if(verbosity_lev) cout<<"lambda_opt_10 = "<<lambda_balance_10<<endl;

  //#############################################################################################################################################
  //#############################################################################################################################################
  //#############################################################################################################################################
  //#############################################################################################################################################



  k=0.01;
  l_low =0.0;
  if(verbosity_lev) cout<<"Finding lambda from balance condition A=(k="<<to_string_with_precision(k,2)<<")*mult*B, mult: "<<mult<<endl;
  bool lambda_balance_found_100=false;
 
 
  //bisection search for condition A = mult*B
  while( !lambda_balance_found_100 ) {

   

    //evaluate the minimum at midpoint
    PrecFloat lambda_mid = (Nit_100==0)?l_up:(l_up+l_low)/2;
    PrecMatr C = Atr*(1-lambda_mid)/M2 + B*lambda_mid/(MODE=="SANF"?M2:1);
    PrecMatr C_inv = C.inverse();
    PrecVect ft_l = ft*(1-lambda_mid)/M2;
    vector<PrecVect> ft_l_jack;
    if(JackOnKer) for(int ijack=0;ijack<Njacks;ijack++) ft_l_jack.push_back( ft_jack[ijack]*(1-lambda_mid)/M2);
    PrecVect M_tilde_n;
    PrecMatr G_n;

    //Lagrangian multipliers
    
    Get_M_tilde_N(ft_l, C_inv, Rt_n, M_tilde_n);
    Get_G_matrix(G_n, C_inv, Rt_n);
    
      
    
    PrecVect gm = C_inv*ft_l;

    vector<PrecVect> gm_jack;

    PrecVect pn, M_n_diff;

    PrecMatr G_n_inv;


    if(Nmoms > 0) {
      G_n_inv = G_n.inverse();
      M_n_diff = M_n - M_tilde_n;
      pn = G_n_inv*M_n_diff;
    }

    for(int n=0;n<Nmoms;n++) {
      PrecVect lmult_n = C_inv*Rt_n[n]; 
      gm = gm + lmult_n*pn(n);
    }


    if(JackOnKer) {

      for(int ijack=0; ijack<Njacks;ijack++) {
	PrecVect gm_ij = C_inv*ft_l_jack[ijack];
	PrecVect M_tilde_n_ij, M_n_diff_ij, pn_ij;
	Get_M_tilde_N(ft_l_jack[ijack], C_inv, Rt_n, M_tilde_n_ij);
	if(Nmoms>0) {
	  PrecVect M_n_diff_ij = M_n_jack[ijack] - M_tilde_n_ij;
	  PrecVect pn_ij = G_n_inv*M_n_diff_ij;
	}
	for(int n=0; n<Nmoms;n++) {
	  PrecVect lmult_n= C_inv*Rt_n[n];
	  gm_ij = gm_ij + lmult_n*pn_ij(n);
	}
	gm_jack.push_back(gm_ij);
      }
      
    }

    PrecFloat A1_val = A1(gm);
    PrecFloat B1_val = B1(gm);
    PrecFloat W_val = (1-lambda_mid)*A1_val + lambda_mid*B1_val;
    
    PrecFloat A1_val_std = A1_std(gm);
    PrecFloat B1_val_std = B1_std(gm);
    PrecFloat W_val_std = (1-lambda_mid)*A1_val_std + lambda_mid*B1_val_std;

    PrecFloat A1_val_std_Emax = A1_std_Emax(gm);
    PrecFloat B1_val_std_Emax = B1_std_Emax(gm);


  
      
    double mult_est = (mult_estimated_from_norm0==false)?(A1_val/B1_val).get():(A1_val_std/B1_val_std).get();
         
    if(k*mult> mult_est) { // lambda_mid is new l_low
      l_low =lambda_mid;
    }
    else { //lambda_mid is new l_up
      l_up = lambda_mid;
    }

    lambda_balance_100= lambda_mid.get();
      
    Nit_100++;
    if(fabs(mult_est - k*mult)/(k*mult) < 0.01) lambda_balance_found_100=true;


    double syst=0.0;

    if(Use_guess_density) {

      auto integrand_syst = [&tmin, &tmax, &T, &gm, &f, &mean, &sigma, &Estart](double E) ->double {
	PrecVect bt;
	Get_bt(bt, E, T, tmin,tmax);
	return (f(E,mean, sigma, Estart, -1) - gm.transpose()*bt).get();
      };

      if(extended_analysis_of_syst)  syst= syst_func(integrand_syst);


      if(print_reco_in_stability_analysis) {
	int Npoints;
	if(analysis_name == "tau_decay") {Npoints= 2000; step_size=0.005;}
	else Npoints= (int)(((mean+20*sigma -Estart)/(sigma*step_size)));
	Vfloat Erg, Error, Spec_dens_guess, Agvs;

	for(int ip=0; ip<Npoints;ip++) {
	    
	  double E;
	  Agvs.push_back(A1_val_std.get());
	  if(analysis_name == "tau_decay") E= (Estart+ ip*step_size);
	  else E = (Estart + (ip*sigma)*step_size);
	  Error.push_back(integrand_syst(E));
	  Erg.push_back(E);
	  Spec_dens_guess.push_back(guess_density(E));
	}

	Print_To_File({}, {Agvs, Erg, Error, Spec_dens_guess}, "../data/spectral_reconstruction/"+analysis_name+"/error_funcs/beta_"+to_string_with_precision(Beta,2)+"/"+to_string(Global_id)+"_"+out_path,"", "#E   |reco-exact|   guess ");


      }
    }

    bool skipping_N_100=false;

    if(Nit_100 > 20) {
     
      //string msg = "Warning: After "+to_string(Nit_100)+" iterations, balance condition A = k*mult*B cannot be obtained for CORR: "+CORR_NAME+" , MODE: "+MODE+", CURR_TYPE: "+curr_type+", mult(target) = "+to_string_with_precision( k*mult, 8)+" a*sigma: "+to_string_with_precision(sigma,5)+", aE*: "+to_string_with_precision(mean, 5)+" lambda: "+to_string_with_precision(lambda_mid, 5)+" , mult: "+to_string_with_precision(mult_est,5);

      //cout<<msg<<endl;
      lambda_balance_found_100=true;
      skipping_N_100=true;
    }


    //##########################################################################################
    //compute anti-Laplace transform corresponding to lambda_mid:
    distr_t R_E_lambda;
      
    for(int ijack=0; ijack< corr.distr_list[1].size(); ijack++) {

      PrecFloat spec_lambda_d_jack=0;

      for(int t=tmin;t<=tmax;t++) spec_lambda_d_jack = spec_lambda_d_jack + ((JackOnKer==true)?gm_jack[ijack](t-tmin):gm(t-tmin))*corr.distr_list[t].distr[ijack]; 

      R_E_lambda.distr.push_back( spec_lambda_d_jack.get());
    }

    Print_R_at_lambda<<lambda_mid<<" "<<A1_val_std<<" "<<A1_val_std_Emax<<" "<<B1_val_std<<" "<<B1_val_std_Emax<<" "<<A1_val<<" "<<(Prefact*(R_E_lambda)+offset).ave()<<" "<<(Prefact*(R_E_lambda)+offset).err();
    Print_R_at_lambda<<" "<<syst;
    Print_R_at_lambda<<" "<<3*(lambda_balance_found_100==1 && skipping_N_100==0)<<" "<<mult_est<<endl;

    //##########################################################################################

     
    Global_id++;
   
  }

  if(verbosity_lev) cout<<"lambda_opt_100 = "<<lambda_balance_100<<endl;


 
  //#############################################################################################################################################
  //#############################################################################################################################################
  //#############################################################################################################################################
  //#############################################################################################################################################

  

  lambda_opt_10 = lambda_balance_10;
  lambda_opt =  lambda_balance;
 
  Print_R_at_lambda.close();

  
  return;

}

distr_t Get_Laplace_transfo( double mean, double sigma, double Estart, int T, int tmax, int prec, string SMEARING_FUNC, const function<PrecFloat(const PrecFloat&, const PrecFloat&,const PrecFloat&,const PrecFloat&, int)> &f, const distr_t_list &corr, double &syst,const double mult, double& pull, string MODE, string reg_type, string CORR_NAME, double Ag_ov_A0_target, bool JackOnKer, const distr_t &Prefact,const double &offset,  string analysis_name, Vfloat &covariance, const function<double(const function<double(double)>&)> &syst_func, bool Use_guess_density, const function<double(double)> &guess_density, bool Int_up_to_Max, double Max_Erg, double b, bool ONLY_FW,  bool GENERALIZED_NORM,  const function<PrecFloat(const PrecFloat &, const PrecFloat &, const PrecFloat &, const PrecFloat &, int)> F_NORM, const function<PrecFloat( PrecFloat )> Atr_gen_NORM) {


  cout<<"precision used: "<<prec<<endl;

  Integrate_up_to_max_energy= Int_up_to_Max;
  Emax_int=Max_Erg;
  Beta=b;
  ONLY_FORWARD=ONLY_FW;
  USE_GENERALIZED_NORM= GENERALIZED_NORM;

 
  if(MODE != "TANT" && MODE != "SANF") crash("MODE: "+MODE+" not recognized");

  if( (Integrate_up_to_max_energy==false) && (Beta >= 2.0)) crash("Cannot use Beta >=2.0 without a finite Emax. Beta: "+to_string_with_precision(Beta,3));

  int Njacks= corr.distr_list[0].distr.size();

  string Emax_val= (Integrate_up_to_max_energy==0)?"inf":to_string_with_precision(Emax_int,1);
 

  //create output directory
  boost::filesystem::create_directory("../data/spectral_reconstruction");
  boost::filesystem::create_directory("../data/spectral_reconstruction/"+analysis_name);
  boost::filesystem::create_directory("../data/spectral_reconstruction/"+analysis_name+"/smearing_func");
  boost::filesystem::create_directory("../data/spectral_reconstruction/"+analysis_name+"/lambda_stability");
  boost::filesystem::create_directory("../data/spectral_reconstruction/"+analysis_name+"/error_funcs");
  boost::filesystem::create_directory("../data/spectral_reconstruction/"+analysis_name+"/smearing_func/beta_"+to_string_with_precision(Beta,2)+"_Emax_"+Emax_val);
  boost::filesystem::create_directory("../data/spectral_reconstruction/"+analysis_name+"/lambda_stability/beta_"+to_string_with_precision(Beta,2)+"_Emax_"+Emax_val);
  boost::filesystem::create_directory("../data/spectral_reconstruction/"+analysis_name+"/error_funcs/beta_"+to_string_with_precision(Beta,2)+"_Emax_"+Emax_val);

  
  PrecFloat::setDefaultPrecision(prec);

  PrecFloat s = sigma;
  PrecFloat m = mean;
  PrecFloat E0 = Estart;

 

  PrecMatr Atr, Atr_std, Atr_std_Emax, B;
  PrecMatr Atr_10;
  PrecVect ft, ft_std, ft_std_Emax, ft_10;
  PrecFloat M2, M2_std, M2_std_Emax;


  //vectors for Lagrangian multiplier
  vector<PrecVect> Rt_n;
  PrecVect M_n;
  PrecVect M_tilde_n, M_tilde_n_10;
  PrecMatr G_n, G_n_10;

   
  //get vectors ft, and matrix Atr

  if(verbosity_lev) cout<<"Njacks: "<<Njacks<<" T: "<<corr.size()<<endl<<flush;

  if(verbosity_lev) cout<<"computing f(t)..."<<flush;

  
  Get_ft(ft, E0, m, s, -1, T, 1, tmax, SMEARING_FUNC, f, F_NORM);
  M2=Get_M2(m,s,E0,-1,f, F_NORM);

  if(verbosity_lev) { cout<<"done!"<<endl<<flush;}

   

  
  if( (Beta != 0) || Integrate_up_to_max_energy) {
    if(verbosity_lev) cout<<"computing f(t)_std..."<<flush;
    Get_ft_std(ft_std, E0, m, s, -1, T, 1, tmax, SMEARING_FUNC, f);
    M2_std= Get_M2_std_norm(m,s,E0,-1,f);
    if(verbosity_lev) cout<<"done!"<<endl<<flush;
  }
  else { ft_std= ft; M2_std=M2; }
  

  if( (Beta != 0) && Integrate_up_to_max_energy) {
    if(verbosity_lev) cout<<"computing f(t)_std_Emax..."<<flush;
    Get_ft_std_Emax(ft_std_Emax, E0, m, s, -1, T, 1, tmax, SMEARING_FUNC, f);
    M2_std_Emax=Get_M2_std_norm_Emax(m,s,E0,-1,f);
    if(verbosity_lev) cout<<"done!"<<endl<<flush;
  }
  else { ft_std_Emax= ft; M2_std_Emax= M2; }

  
  if(verbosity_lev>=2) {
    cout.precision(PrecFloat::getNDigits());
    cout<<"Printing ft: "<<endl;
    for(int t=1;t<=tmax;t++) cout<<"f("<<t<<") : "<<ft(t-1)<<" "<<ft_std(t-1)<<" "<<ft_std_Emax(t-1)<<endl<<flush;
  }
    
  Get_Atr(Atr, E0, T, 1, tmax, Atr_gen_NORM);

  Get_Atr_std(Atr_std, E0, T, 1, tmax);

  Get_Atr_std_Emax(Atr_std_Emax, E0, T, 1, tmax);

  Get_Rt_up_to_N(E0, T, 1, tmax, Rt_n);

  Get_M_N(m,s,E0,-1,f, M_n);

  if(verbosity_lev>=2) {
    cout.precision(PrecFloat::getNDigits());
    cout<<"M2 : "<<M2<<endl<<flush;
    cout<<"M2 std: "<<M2_std<<endl<<flush;
    cout<<"M2_std_Emax: "<<M2_std_Emax<<endl<<flush;
  }
 

  cout.precision(10);
  //get ft, M_N in case a statistical fluctuations in kernel function must be computed
  vector<PrecVect> ft_jack, ft_jack_10;
  vector<PrecVect> M_n_jack, M_n_diff_jack;
  vector<PrecVect> M_n_jack_10, M_n_diff_jack_10;
 

  if(JackOnKer) {
    //compute jackknife distribution of ft,  M_n
    for(int ijack=0;ijack<Njacks;ijack++) {
      PrecVect ft_ij;
      PrecVect M_n_ij;

      Get_ft(ft_ij, E0, m , s, ijack, T,1,tmax, SMEARING_FUNC, f, F_NORM);
      Get_M_N(m,s,E0,ijack, f, M_n_ij);

      ft_jack.push_back(ft_ij);
      ft_jack_10.push_back(ft_ij);
      M_n_jack.push_back(M_n_ij);
    }

  }

  if(INCLUDE_ERRORS) Compute_covariance_matrix(B,Atr, corr,1,tmax, m,s, analysis_name, MODE, covariance);

  double lambda_opt= lambda;
  double lambda_opt_10=lambda;

  distr_t offset_distr= Get_id_jack_distr(Njacks)*offset;

  if(INCLUDE_ERRORS && FIND_OPTIMAL_LAMBDA) Get_optimal_lambda(Atr, Atr_std, Atr_std_Emax,  B, ft, ft_jack, ft_std, ft_std_Emax, M2, M2_std, M2_std_Emax, mean, sigma, Estart, lambda_opt , lambda_opt_10, Rt_n, M_n, M_n_jack, corr, T , 1 , tmax, mult,  MODE, reg_type, SMEARING_FUNC,  CORR_NAME, Ag_ov_A0_target, JackOnKer, Prefact, offset_distr, analysis_name, f, syst_func, Use_guess_density, guess_density);


  if(verbosity_lev>=2) cout<<"Stability analysis completed!"<<endl;
    							         
  if(INCLUDE_ERRORS) {
    Atr_10 = Atr*(1-lambda_opt_10)/M2 + B*lambda_opt_10/((MODE=="SANF")?M2:1);
    Atr = Atr*(1-lambda_opt)/M2 + B*lambda_opt/((MODE=="SANF")?M2:1);
    ft_10 = ft*(1-lambda_opt_10)/M2;
    ft= ft*(1-lambda_opt)/M2;
    if(JackOnKer) {
      for(int ijack=0;ijack<Njacks;ijack++) {
	ft_jack_10[ijack] = ft_jack_10[ijack]*(1-lambda_opt_10)/M2;
	ft_jack[ijack] = ft_jack[ijack]*(1-lambda_opt)/M2;
      }
    }
  }


  //invert Atr

  const PrecMatr Atr_inv = Atr.inverse();
  const PrecMatr Atr_inv_10= Atr_10.inverse();


  
  //to compute Lagrangian multiplier

  
  Get_M_tilde_N(ft, Atr_inv, Rt_n, M_tilde_n);
  Get_M_tilde_N(ft_10, Atr_inv_10, Rt_n, M_tilde_n_10);
  Get_G_matrix(G_n, Atr_inv, Rt_n);
  Get_G_matrix(G_n_10, Atr_inv_10, Rt_n);

  if(JackOnKer && Nmoms > 0) {
    for(int ijack=0; ijack<Njacks;ijack++) {
      PrecVect M_tilde_n_ij;
      PrecVect M_tilde_n_ij_10;
      Get_M_tilde_N(ft_jack[ijack], Atr_inv, Rt_n, M_tilde_n_ij);
      Get_M_tilde_N(ft_jack_10[ijack], Atr_inv_10, Rt_n, M_tilde_n_ij_10);
      M_n_diff_jack.push_back(M_n_jack[ijack] - M_tilde_n_ij);
      M_n_diff_jack_10.push_back( M_n_jack[ijack] - M_tilde_n_ij_10);
    }
  }
  
  

  //get g(t) coefficient vector
  

  PrecVect g = Atr_inv*ft;
  PrecVect g_10= Atr_inv_10*ft_10;

  vector<PrecVect> g_jack;
  vector<PrecVect> g_jack_10;

  //add Lagrangian multipliers


  PrecMatr G_n_inv, G_n_inv_10;
  PrecVect pn, pn_10;

  if(Nmoms>0) { G_n_inv= G_n.inverse(); G_n_inv_10=G_n_10.inverse();}

  PrecVect M_n_diff, M_n_diff_10;

  if(Nmoms > 0) {

    M_n_diff= M_n - M_tilde_n;
    M_n_diff_10= M_n - M_tilde_n_10;
    pn= G_n_inv*M_n_diff;
    pn_10= G_n_inv_10*M_n_diff_10;
  }


  for(int n=0;n <Nmoms;n++) {
    PrecVect lmult_n = Atr_inv*Rt_n[n];
    PrecVect lmult_n_10= Atr_inv_10*Rt_n[n];
    g = g + lmult_n*pn(n);
    g_10= g_10 + lmult_n_10*pn_10(n);
  }

  if(JackOnKer) {
    for(int ijack=0;ijack<Njacks;ijack++) {
      PrecVect g_ij = Atr_inv*ft_jack[ijack];
      PrecVect g_ij_10= Atr_inv_10*ft_jack_10[ijack];
      PrecVect pn_ij, pn_ij_10;
      if(Nmoms>0) { pn_ij= G_n_inv*M_n_diff_jack[ijack]; pn_ij_10= G_n_inv_10*M_n_diff_jack_10[ijack];}
      for(int n=0; n < Nmoms; n++) {
	PrecVect lmult_n= Atr_inv*Rt_n[n];
	PrecVect lmult_n_10= Atr_inv_10*Rt_n[10];
	g_ij = g_ij + lmult_n*pn_ij(n);
	g_ij_10= g_ij_10 + lmult_n_10*pn_ij_10(n);
      }
      g_jack.push_back(g_ij);
      g_jack_10.push_back(g_ij_10);
    }
  }


  if(verbosity_lev>=2) cout<<"printing target & reco smearing func"<<endl<<flush; 
  //print reconstructed smearing func, exact smearing func, diff

  //compute it from E0 up to 10 E* with a step size of sigma * step_size

  vector<double> Reco;
  vector<double> Exact;
  vector<double> Err;
  vector<double> Erg;
  vector<double> Spec_dens_guess;


  int Npoints;
  if(analysis_name == "tau_decay") {Npoints= 20000; step_size=0.001;}
  else Npoints= (int)(((m+20*s)/(s*step_size)).get());
  for(int ip=0; ip<Npoints;ip++) {
    PrecVect bt;
    PrecFloat E;
    if(analysis_name == "tau_decay") E= 0.25*E0+ ip*step_size;
    else E = (ip*s)*step_size;
    //if(verbosity_lev==2) cout<<"ip: "<<ip<<"/"<<Npoints<<"  (E-m): "<<((E-m)/s).get()<<flush;
    Get_bt(bt, E, T, 1, tmax);  
    const PrecFloat reco_result = g.transpose()*bt;
    const PrecFloat exact_result = f(E,m,s,E0,-1);
    //if(verbosity_lev==2) cout<<"Exact and reco for ip: "<<ip<<" computed"<<endl<<flush;
    Reco.push_back( reco_result.get());
    Exact.push_back( exact_result.get());
    Err.push_back(  (exact_result -reco_result).get());
    Erg.push_back(E.get());
    if(Use_guess_density) Spec_dens_guess.push_back( guess_density(E.get()));
  }

  if(verbosity_lev==2) cout<<"data created.."<<endl<<flush;
  
  //print to file
  if(!Use_guess_density) {
    Print_To_File({}, { Erg, Reco, Exact, Err}, "../data/spectral_reconstruction/"+analysis_name+"/smearing_func/beta_"+to_string_with_precision(Beta,2)+"_Emax_"+Emax_val+"/"+MODE+"_"+CORR_NAME+"_"+reg_type+"_"+SMEARING_FUNC+"_func_E*_"+to_string_with_precision(mean,3)+"_sigma_"+to_string_with_precision(sigma,3)+"_E0_"+to_string_with_precision(E0.get(),3)+"_T_"+to_string(T)+"_tmax_"+to_string(tmax)+"_alpha_"+to_string(alpha)+"_beta_"+to_string_with_precision(Beta,2)+".dat", "", "#id  E   reco    exact    error");
  }
  else {
    Print_To_File({}, { Erg, Reco, Exact, Err, Spec_dens_guess}, "../data/spectral_reconstruction/"+analysis_name+"/smearing_func/beta_"+to_string_with_precision(Beta,2)+"_Emax_"+Emax_val+"/"+MODE+"_"+CORR_NAME+"_"+reg_type+"_"+SMEARING_FUNC+"_func_E*_"+to_string_with_precision(mean,3)+"_sigma_"+to_string_with_precision(sigma,3)+"_E0_"+to_string_with_precision(E0.get(),3)+"_T_"+to_string(T)+"_tmax_"+to_string(tmax)+"_alpha_"+to_string(alpha)+"_beta_"+to_string_with_precision(Beta,2)+".dat", "", "#id  E   reco    exact    error guess_density");

  }

  if(verbosity_lev>=2) cout<<"output written to file!"<<endl<<flush;

  distr_t Spec_dens_at_E_star; //Uses Jackknife distr by default
  distr_t Spec_dens_at_E_star_10;

  if(verbosity_lev>=2) cout<<"Computing ave+stat+syst on spec density"<<endl<<flush;


  if(INCLUDE_ERRORS) {

    //build the spectral density at E*
    

    for(int ijack=0; ijack< Njacks; ijack++) {

      PrecFloat spec_d_jack=0;
      PrecFloat spec_d_jack_10=0;
      for(int t=1;t<=tmax;t++)  {
	spec_d_jack = spec_d_jack + ((JackOnKer)?g_jack[ijack](t-1):g(t-1))*corr.distr_list[t].distr[ijack];
	spec_d_jack_10 = spec_d_jack_10 + ((JackOnKer)?g_jack_10[ijack](t-1):g_10(t-1))*corr.distr_list[t].distr[ijack];
      }

      Spec_dens_at_E_star.distr.push_back( spec_d_jack.get());
      Spec_dens_at_E_star_10.distr.push_back( spec_d_jack_10.get());
    }

    Spec_dens_at_E_star = (Spec_dens_at_E_star)*Prefact;
    Spec_dens_at_E_star_10 = (Spec_dens_at_E_star_10)*Prefact;

    
    cout<<"done"<<endl<<flush;
    
    //print coefficient
    ofstream PrintCoeff("../data/spectral_reconstruction/"+analysis_name+"/smearing_func/beta_"+to_string_with_precision(Beta,2)+"_Emax_"+Emax_val+"/"+MODE+"_"+CORR_NAME+"_"+reg_type+"_"+SMEARING_FUNC+"_coeff_E*"+to_string_with_precision(mean,3)+"_sigma_"+to_string_with_precision(sigma,3)+"_E0_"+to_string_with_precision(E0.get(),3)+"_T_"+to_string(T)+"_tmax_"+to_string(tmax)+"_alpha_"+to_string(alpha)+"_beta_"+to_string_with_precision(Beta,3)+".dat");
    if(!INCLUDE_ERRORS) {
      PrintCoeff<<g<<endl;
    }
    else {
      if(!JackOnKer) {
	for(int t=1;t<=tmax;t++) PrintCoeff<<g(t-1)<<"\t"<<g(t-1)*(Prefact*corr).ave(t)<<"\t"<<Spec_dens_at_E_star.ave()<<"\t"<<Spec_dens_at_E_star.err()<<endl;
      }
      else {
	for(int t=1; t<=tmax;t++) {
	  distr_t g_t;
	  for(int ijack=0;ijack<Njacks;ijack++) g_t.distr.push_back( g_jack[ijack](t-1).get());
	  PrintCoeff<<t<<"\t"<<g_t.ave()<<"  "<<g_t.err()<<"\t"<<(g_t*Prefact*corr).ave(t)<<"   "<<(g_t*Prefact*corr).err(t)<<"\t"<<Spec_dens_at_E_star.ave()<<"\t"<<Spec_dens_at_E_star.err()<<endl;
	}
      }
    }
    PrintCoeff.close();





    
        
    syst = erf(fabs( (Spec_dens_at_E_star - Spec_dens_at_E_star_10).ave()/(sqrt(2)*Spec_dens_at_E_star_10.err())))*fabs( (Spec_dens_at_E_star - Spec_dens_at_E_star_10).ave());
    pull= (Spec_dens_at_E_star - Spec_dens_at_E_star_10).ave()/(Spec_dens_at_E_star_10.err());
  }
    

  //set lambda to lambda_optimal


  
  return Spec_dens_at_E_star;  
  
}





distr_t Get_Laplace_transfo_tmin( double mean, double sigma, double Estart, int T,int tmin,  int tmax, int prec, string SMEARING_FUNC, const function<PrecFloat(const PrecFloat&, const PrecFloat&,const PrecFloat&,const PrecFloat&, int)> &f, const distr_t_list &corr, double &syst,const double mult, double& pull, string MODE, string reg_type, string CORR_NAME, double Ag_ov_A0_target, bool JackOnKer, const distr_t &Prefact,const distr_t &offset,  string analysis_name, Vfloat &covariance, const function<double(const function<double(double)>&)> &syst_func, bool Use_guess_density, const function<double(double)> &guess_density, bool Int_up_to_Max, double Max_Erg, double b, bool ONLY_FW,  bool GENERALIZED_NORM,  const function<PrecFloat(const PrecFloat &, const PrecFloat &, const PrecFloat &, const PrecFloat &, int)> F_NORM, const function<PrecFloat( PrecFloat )> Atr_gen_NORM) {


  cout<<"precision used: "<<prec<<endl;

  Integrate_up_to_max_energy= Int_up_to_Max;
  Emax_int=Max_Erg;
  Beta=b;
  ONLY_FORWARD=ONLY_FW;
  USE_GENERALIZED_NORM= GENERALIZED_NORM;

 
  if(MODE != "TANT" && MODE != "SANF") crash("MODE: "+MODE+" not recognized");

  if( (Integrate_up_to_max_energy==false) && (Beta >= 2.0)) crash("Cannot use Beta >=2.0 without a finite Emax. Beta: "+to_string_with_precision(Beta,3));

  int Njacks= corr.distr_list[0].distr.size();

  string Emax_val= (Integrate_up_to_max_energy==0)?"inf":to_string_with_precision(Emax_int,1);
 

  //create output directory
  boost::filesystem::create_directory("../data/spectral_reconstruction");
  boost::filesystem::create_directory("../data/spectral_reconstruction/"+analysis_name);
  boost::filesystem::create_directory("../data/spectral_reconstruction/"+analysis_name+"/smearing_func");
  boost::filesystem::create_directory("../data/spectral_reconstruction/"+analysis_name+"/lambda_stability");
  boost::filesystem::create_directory("../data/spectral_reconstruction/"+analysis_name+"/error_funcs");
  boost::filesystem::create_directory("../data/spectral_reconstruction/"+analysis_name+"/smearing_func/beta_"+to_string_with_precision(Beta,2)+"_Emax_"+Emax_val);
  boost::filesystem::create_directory("../data/spectral_reconstruction/"+analysis_name+"/lambda_stability/beta_"+to_string_with_precision(Beta,2)+"_Emax_"+Emax_val);
  boost::filesystem::create_directory("../data/spectral_reconstruction/"+analysis_name+"/error_funcs/beta_"+to_string_with_precision(Beta,2)+"_Emax_"+Emax_val);

  
  PrecFloat::setDefaultPrecision(prec);

  PrecFloat s = sigma;
  PrecFloat m = mean;
  PrecFloat E0 = Estart;

 

  PrecMatr Atr, Atr_std, Atr_std_Emax, B;
  PrecMatr Atr_10;
  PrecVect ft, ft_std, ft_std_Emax, ft_10;
  PrecFloat M2, M2_std, M2_std_Emax;


  //vectors for Lagrangian multiplier
  vector<PrecVect> Rt_n;
  PrecVect M_n;
  PrecVect M_tilde_n, M_tilde_n_10;
  PrecMatr G_n, G_n_10;

   
  //get vectors ft, and matrix Atr

  if(verbosity_lev) cout<<"Njacks: "<<Njacks<<" T: "<<corr.size()<<endl<<flush;

  if(verbosity_lev) cout<<"computing f(t)..."<<flush;

  
  Get_ft(ft, E0, m, s, -1, T, tmin, tmax, SMEARING_FUNC, f, F_NORM);
  M2=Get_M2(m,s,E0,-1,f, F_NORM);

  if(verbosity_lev) { cout<<"done!"<<endl<<flush;}

   

  
  if( (Beta != 0) || Integrate_up_to_max_energy) {
    if(verbosity_lev) cout<<"computing f(t)_std..."<<flush;
    Get_ft_std(ft_std, E0, m, s, -1, T, tmin, tmax, SMEARING_FUNC, f);
    M2_std= Get_M2_std_norm(m,s,E0,-1,f);
    if(verbosity_lev) cout<<"done!"<<endl<<flush;
  }
  else { ft_std= ft; M2_std=M2; }
  

  if( (Beta != 0) && Integrate_up_to_max_energy) {
    if(verbosity_lev) cout<<"computing f(t)_std_Emax..."<<flush;
    Get_ft_std_Emax(ft_std_Emax, E0, m, s, -1, T, tmin, tmax, SMEARING_FUNC, f);
    M2_std_Emax=Get_M2_std_norm_Emax(m,s,E0,-1,f);
    if(verbosity_lev) cout<<"done!"<<endl<<flush;
  }
  else { ft_std_Emax= ft; M2_std_Emax= M2; }

  
  if(verbosity_lev>=2) {
    cout.precision(PrecFloat::getNDigits());
    cout<<"Printing ft: "<<endl;
    for(int t=tmin;t<=tmax;t++) cout<<"f("<<t<<") : "<<ft(t-tmin)<<" "<<ft_std(t-tmin)<<" "<<ft_std_Emax(t-tmin)<<endl<<flush;
  }
    
  Get_Atr(Atr, E0, T, tmin, tmax, Atr_gen_NORM);

  Get_Atr_std(Atr_std, E0, T, tmin, tmax);

  Get_Atr_std_Emax(Atr_std_Emax, E0, T, tmin, tmax);

  Get_Rt_up_to_N(E0, T, tmin, tmax, Rt_n);

  Get_M_N(m,s,E0,-1,f, M_n);

  if(verbosity_lev>=2) {
    cout.precision(PrecFloat::getNDigits());
    cout<<"M2 : "<<M2<<endl<<flush;
    cout<<"M2 std: "<<M2_std<<endl<<flush;
    cout<<"M2_std_Emax: "<<M2_std_Emax<<endl<<flush;
  }
 

  cout.precision(10);
  //get ft, M_N in case a statistical fluctuations in kernel function must be computed
  vector<PrecVect> ft_jack, ft_jack_10;
  vector<PrecVect> M_n_jack, M_n_diff_jack;
  vector<PrecVect> M_n_jack_10, M_n_diff_jack_10;
 

  if(JackOnKer) {
    //compute jackknife distribution of ft,  M_n
    for(int ijack=0;ijack<Njacks;ijack++) {
      PrecVect ft_ij;
      PrecVect M_n_ij;

      Get_ft(ft_ij, E0, m , s, ijack, T,tmin,tmax, SMEARING_FUNC, f, F_NORM);
      Get_M_N(m,s,E0,ijack, f, M_n_ij);

      ft_jack.push_back(ft_ij);
      ft_jack_10.push_back(ft_ij);
      M_n_jack.push_back(M_n_ij);
    }

  }

  if(INCLUDE_ERRORS) Compute_covariance_matrix(B,Atr, corr,tmin,tmax, m,s, analysis_name, MODE, covariance);

  double lambda_opt= lambda;
  double lambda_opt_10=lambda;

  if(INCLUDE_ERRORS && FIND_OPTIMAL_LAMBDA) Get_optimal_lambda(Atr, Atr_std, Atr_std_Emax,  B, ft, ft_jack, ft_std, ft_std_Emax, M2, M2_std, M2_std_Emax, mean, sigma, Estart, lambda_opt , lambda_opt_10, Rt_n, M_n, M_n_jack, corr, T , tmin , tmax, mult,  MODE, reg_type, SMEARING_FUNC,  CORR_NAME, Ag_ov_A0_target, JackOnKer, Prefact, offset, analysis_name, f, syst_func, Use_guess_density, guess_density);


  if(verbosity_lev>=2) cout<<"Stability analysis completed!"<<endl;
    							         
  if(INCLUDE_ERRORS) {
    Atr_10 = Atr*(1-lambda_opt_10)/M2 + B*lambda_opt_10/((MODE=="SANF")?M2:1);
    Atr = Atr*(1-lambda_opt)/M2 + B*lambda_opt/((MODE=="SANF")?M2:1);
    ft_10 = ft*(1-lambda_opt_10)/M2;
    ft= ft*(1-lambda_opt)/M2;
    if(JackOnKer) {
      for(int ijack=0;ijack<Njacks;ijack++) {
	ft_jack_10[ijack] = ft_jack_10[ijack]*(1-lambda_opt_10)/M2;
	ft_jack[ijack] = ft_jack[ijack]*(1-lambda_opt)/M2;
      }
    }
  }


  //invert Atr

  const PrecMatr Atr_inv = Atr.inverse();
  const PrecMatr Atr_inv_10= Atr_10.inverse();


  
  //to compute Lagrangian multiplier

  
  Get_M_tilde_N(ft, Atr_inv, Rt_n, M_tilde_n);
  Get_M_tilde_N(ft_10, Atr_inv_10, Rt_n, M_tilde_n_10);
  Get_G_matrix(G_n, Atr_inv, Rt_n);
  Get_G_matrix(G_n_10, Atr_inv_10, Rt_n);

  if(JackOnKer && Nmoms > 0) {
    for(int ijack=0; ijack<Njacks;ijack++) {
      PrecVect M_tilde_n_ij;
      PrecVect M_tilde_n_ij_10;
      Get_M_tilde_N(ft_jack[ijack], Atr_inv, Rt_n, M_tilde_n_ij);
      Get_M_tilde_N(ft_jack_10[ijack], Atr_inv_10, Rt_n, M_tilde_n_ij_10);
      M_n_diff_jack.push_back(M_n_jack[ijack] - M_tilde_n_ij);
      M_n_diff_jack_10.push_back( M_n_jack[ijack] - M_tilde_n_ij_10);
    }
  }
  
  

  //get g(t) coefficient vector
  

  PrecVect g = Atr_inv*ft;
  PrecVect g_10= Atr_inv_10*ft_10;

  vector<PrecVect> g_jack;
  vector<PrecVect> g_jack_10;

  //add Lagrangian multipliers


  PrecMatr G_n_inv, G_n_inv_10;
  PrecVect pn, pn_10;

  if(Nmoms>0) { G_n_inv= G_n.inverse(); G_n_inv_10=G_n_10.inverse();}

  PrecVect M_n_diff, M_n_diff_10;

  if(Nmoms > 0) {

    M_n_diff= M_n - M_tilde_n;
    M_n_diff_10= M_n - M_tilde_n_10;
    pn= G_n_inv*M_n_diff;
    pn_10= G_n_inv_10*M_n_diff_10;
  }


  for(int n=0;n <Nmoms;n++) {
    PrecVect lmult_n = Atr_inv*Rt_n[n];
    PrecVect lmult_n_10= Atr_inv_10*Rt_n[n];
    g = g + lmult_n*pn(n);
    g_10= g_10 + lmult_n_10*pn_10(n);
  }

  if(JackOnKer) {
    for(int ijack=0;ijack<Njacks;ijack++) {
      PrecVect g_ij = Atr_inv*ft_jack[ijack];
      PrecVect g_ij_10= Atr_inv_10*ft_jack_10[ijack];
      PrecVect pn_ij, pn_ij_10;
      if(Nmoms>0) { pn_ij= G_n_inv*M_n_diff_jack[ijack]; pn_ij_10= G_n_inv_10*M_n_diff_jack_10[ijack];}
      for(int n=0; n < Nmoms; n++) {
	PrecVect lmult_n= Atr_inv*Rt_n[n];
	PrecVect lmult_n_10= Atr_inv_10*Rt_n[10];
	g_ij = g_ij + lmult_n*pn_ij(n);
	g_ij_10= g_ij_10 + lmult_n_10*pn_ij_10(n);
      }
      g_jack.push_back(g_ij);
      g_jack_10.push_back(g_ij_10);
    }
  }


  if(verbosity_lev>=2) cout<<"printing target & reco smearing func"<<endl<<flush; 
  //print reconstructed smearing func, exact smearing func, diff

  //compute it from E0 up to 10 E* with a step size of sigma * step_size

  vector<double> Reco;
  vector<double> Exact;
  vector<double> Err;
  vector<double> Erg;
  vector<double> Spec_dens_guess;


  int Npoints;
  if(analysis_name == "tau_decay") {Npoints= 20000; step_size=0.001;}
  else Npoints= (int)(((m+20*s)/(s*step_size)).get());
  for(int ip=0; ip<Npoints;ip++) {
    PrecVect bt;
    PrecFloat E;
    if(analysis_name == "tau_decay") E= 0.25*E0+ ip*step_size;
    else E = (ip*s)*step_size;
    //if(verbosity_lev==2) cout<<"ip: "<<ip<<"/"<<Npoints<<"  (E-m): "<<((E-m)/s).get()<<flush;
    Get_bt(bt, E, T, tmin, tmax);  
    const PrecFloat reco_result = g.transpose()*bt;
    const PrecFloat exact_result = f(E,m,s,E0,-1);
    //if(verbosity_lev==2) cout<<"Exact and reco for ip: "<<ip<<" computed"<<endl<<flush;
    Reco.push_back( reco_result.get());
    Exact.push_back( exact_result.get());
    Err.push_back(  (exact_result -reco_result).get());
    Erg.push_back(E.get());
    if(Use_guess_density) Spec_dens_guess.push_back( guess_density(E.get()));
  }

  if(verbosity_lev==2) cout<<"data created.."<<endl<<flush;
  
  //print to file
  if(!Use_guess_density) {
    Print_To_File({}, { Erg, Reco, Exact, Err}, "../data/spectral_reconstruction/"+analysis_name+"/smearing_func/beta_"+to_string_with_precision(Beta,2)+"_Emax_"+Emax_val+"/"+MODE+"_"+CORR_NAME+"_"+reg_type+"_"+SMEARING_FUNC+"_func_E*_"+to_string_with_precision(mean,3)+"_sigma_"+to_string_with_precision(sigma,3)+"_E0_"+to_string_with_precision(E0.get(),3)+"_T_"+to_string(T)+"_tmin_"+to_string(tmin)+"_tmax_"+to_string(tmax)+"_alpha_"+to_string(alpha)+"_beta_"+to_string_with_precision(Beta,2)+".dat", "", "#id  E   reco    exact    error");
  }
  else {
    Print_To_File({}, { Erg, Reco, Exact, Err, Spec_dens_guess}, "../data/spectral_reconstruction/"+analysis_name+"/smearing_func/beta_"+to_string_with_precision(Beta,2)+"_Emax_"+Emax_val+"/"+MODE+"_"+CORR_NAME+"_"+reg_type+"_"+SMEARING_FUNC+"_func_E*_"+to_string_with_precision(mean,3)+"_sigma_"+to_string_with_precision(sigma,3)+"_E0_"+to_string_with_precision(E0.get(),3)+"_T_"+to_string(T)+"_tmin_"+to_string(tmin)+"_tmax_"+to_string(tmax)+"_alpha_"+to_string(alpha)+"_beta_"+to_string_with_precision(Beta,2)+".dat", "", "#id  E   reco    exact    error guess_density");

  }

  if(verbosity_lev>=2) cout<<"output written to file!"<<endl<<flush;

  distr_t Spec_dens_at_E_star; //Uses Jackknife distr by default
  distr_t Spec_dens_at_E_star_10;

  if(verbosity_lev>=2) cout<<"Computing ave+stat+syst on spec density"<<endl<<flush;


  if(INCLUDE_ERRORS) {

    //build the spectral density at E*
    

    for(int ijack=0; ijack< Njacks; ijack++) {

      PrecFloat spec_d_jack=0;
      PrecFloat spec_d_jack_10=0;
      for(int t=tmin;t<=tmax;t++)  {
	spec_d_jack = spec_d_jack + ((JackOnKer)?g_jack[ijack](t-tmin):g(t-tmin))*corr.distr_list[t].distr[ijack];
	spec_d_jack_10 = spec_d_jack_10 + ((JackOnKer)?g_jack_10[ijack](t-tmin):g_10(t-tmin))*corr.distr_list[t].distr[ijack];
      }

      Spec_dens_at_E_star.distr.push_back( spec_d_jack.get());
      Spec_dens_at_E_star_10.distr.push_back( spec_d_jack_10.get());
    }

    Spec_dens_at_E_star = (Spec_dens_at_E_star)*Prefact;
    Spec_dens_at_E_star_10 = (Spec_dens_at_E_star_10)*Prefact;

    
    cout<<"done"<<endl<<flush;
    
    //print coefficient
    ofstream PrintCoeff("../data/spectral_reconstruction/"+analysis_name+"/smearing_func/beta_"+to_string_with_precision(Beta,2)+"_Emax_"+Emax_val+"/"+MODE+"_"+CORR_NAME+"_"+reg_type+"_"+SMEARING_FUNC+"_coeff_E*"+to_string_with_precision(mean,3)+"_sigma_"+to_string_with_precision(sigma,3)+"_E0_"+to_string_with_precision(E0.get(),3)+"_T_"+to_string(T)+"_tmin_"+to_string(tmin)+"_tmax_"+to_string(tmax)+"_alpha_"+to_string(alpha)+"_beta_"+to_string_with_precision(Beta,3)+".dat");
    if(!INCLUDE_ERRORS) {
      PrintCoeff<<g<<endl;
    }
    else {
      if(!JackOnKer) {
	for(int t=tmin;t<=tmax;t++) PrintCoeff<<g(t-tmin)<<"\t"<<g(t-tmin)*(Prefact*corr).ave(t)<<"\t"<<Spec_dens_at_E_star.ave()<<"\t"<<Spec_dens_at_E_star.err()<<endl;
      }
      else {
	for(int t=tmin; t<=tmax;t++) {
	  distr_t g_t;
	  for(int ijack=0;ijack<Njacks;ijack++) g_t.distr.push_back( g_jack[ijack](t-tmin).get());
	  PrintCoeff<<t<<"\t"<<g_t.ave()<<"  "<<g_t.err()<<"\t"<<(g_t*Prefact*corr).ave(t)<<"   "<<(g_t*Prefact*corr).err(t)<<"\t"<<Spec_dens_at_E_star.ave()<<"\t"<<Spec_dens_at_E_star.err()<<endl;
	}
      }
    }
    PrintCoeff.close();





    
        
    syst = erf(fabs( (Spec_dens_at_E_star - Spec_dens_at_E_star_10).ave()/(sqrt(2)*Spec_dens_at_E_star_10.err())))*fabs( (Spec_dens_at_E_star - Spec_dens_at_E_star_10).ave());

    pull =  (Spec_dens_at_E_star - Spec_dens_at_E_star_10).ave()/(Spec_dens_at_E_star_10.err());
  }
    

  //set lambda to lambda_optimal

  //lambda= lambda_opt;
  
  return Spec_dens_at_E_star;  
  
}
