#include "../include/g_minus_2_utilities.h"


using namespace std;

const double m_muon= 0.10565837;  //1.7768 ;
const double relerr= 1e-12;



#pragma omp declare reduction(+: Vfloat : omp_out= Sum_vectors(omp_out,omp_in)) initializer(omp_priv=Vfloat(omp_orig.size(),0))





double kernel_K(double t, double MV) {

  //PrecFloat::setDefaultPrecision(80);

  double m_mu= m_muon;

  if(t < 1e-18) return 0.0; 

  auto F = [&](double x) -> double {

    if (x<1e-10) return 0;
    return (4.0/pow(m_mu*MV,2))*(1.0/sqrt(4.0 + pow(x,2)))*pow(  (sqrt(4.0+pow(x,2))-x)/(sqrt(4.0+pow(x,2))+x),2)*( (cos(m_mu*MV*t*x)-1)/pow(x,2) + (1.0/2.0)*pow(t*m_mu*MV,2));
  };

  double prec= 1e-8; //old_prec is 1e-8
  if ( t*MV*m_mu < 5e-4) prec= 1e-6; 

  if( t*MV*m_mu < 1e-4) prec= 5e-5;

  auto F2 = [&](double x) -> double {


	      double z= t*MV*m_mu;

	      //PrecFloat y_MP= 0.5*z*x/sqrt( 1.0- PrecFloat(x)); 

	      double y= 0.5*z*x/sqrt(1.0-x);

	      if( z < 8e-6) { cout<<"Warning: kernel function evaluated at too small value of the argument and has been set to zero"<<endl; return 0.0;}

	      return (1.0-x)*( 1.0- pow(sin(y)/y,2) );

	    };

  double err;
  //double tol_kernel= 1e-16;
  //double result= boost::math::quadrature::gauss_kronrod<double,61>::integrate(F, 0.0, numeric_limits<double>::infinity() , 5, tol_kernel, &err);


  gsl_function_pp<decltype(F2)> Fp(F2);

 
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (10000);


  gsl_function *G = static_cast<gsl_function*>(&Fp);
  double res_GSL, err_GSL;
 

  gsl_integration_qags(G, 0.0, 1.0, 0.0, prec, 10000, w, &res_GSL, &err_GSL);
  gsl_integration_workspace_free (w);
  if(err_GSL/fabs(res_GSL) > prec*5 && err_GSL/fabs(res_GSL) > 1e-5) crash("GSL integrator did not achieve target precision: "+to_string_with_precision(5*prec,10)+". prec achieved: "+to_string_with_precision( err_GSL/fabs(res_GSL), 7)+"at time t: "+to_string_with_precision(t,8));
  if(res_GSL < 0) crash("kernel_K < 0, for t: "+to_string_with_precision(t,5)+" K(t): "+to_string_with_precision(res_GSL, 13));

  return res_GSL*(1.0/pow(m_mu*MV,2))*pow(t*MV*m_mu,2);

  //if(result < 0) crash("kernel_K < 0, for t: "+to_string_with_precision(t,5)+" K(t): "+to_string_with_precision(result, 13));

  //return result;


}


double kernel_K_old(double t, double MV) {

  double m_mu= m_muon;

  auto F = [&](double x) -> double {

    if (x<1e-40) return 0;
    return (4.0/pow(m_mu*MV,2))*(1.0/sqrt(4.0 + pow(x,2)))*pow(  (sqrt(4.0+pow(x,2))-x)/(sqrt(4.0+pow(x,2))+x),2)*( (cos(m_mu*MV*t*x)-1)/pow(x,2) + (1.0/2.0)*pow(t*m_mu*MV,2));
  };


  double err;
  double tol_kernel= 1e-14;
  return boost::math::quadrature::gauss_kronrod<double,61>::integrate(F, 0.0, numeric_limits<double>::infinity() , 5, tol_kernel, &err);


}


double der_kernel_K_W_win(double t, double MV) {

  auto der_K = [&t](double x) {
		 double tt0=0.4*1.0/0.197327;
		 double tt1=1.0*1.0/0.197327;
		 double dd= 0.15*1.0/0.197327;
		 auto tth0 = [&tt0, &dd](double ta) ->double { return 1.0/(1.0 + exp(-2.0*(ta-tt0)/dd));};
		 auto tth1 = [&tt1,&dd](double ta) ->double { return 1.0/(1.0 + exp(-2.0*(ta-tt1)/dd));};
		 return kernel_K(t,x)*(tth0(t*x) - tth1(t*x));};

  gsl_function_pp<decltype(der_K)> dK(der_K);

  gsl_function *G = static_cast<gsl_function*>(&dK);

  double result, err;
  gsl_deriv_central(G, MV, MV*0.01, &result, &err);

  if(err/result > 1e-2) crash("Failed to evaluate kernel  derivative");

  cout<<"(t, MV): ("<<t<<","<<MV<<")"<<endl;
  cout<<"ker_W(t,MV):"<<der_K(MV)<<endl;
  cout<<"ker_W_derivative: "<<result<<" +- "<<err<<endl;

  return result;
  
  

}

double der_kernel_K_SD_win(double t, double MV) {

  auto der_K = [&t](double x) {
		 double tt0=0.4*1.0/0.197327;
		 double tt1=1.0*1.0/0.197327;
		 double dd= 0.15*1.0/0.197327;
		 auto tth0 = [&tt0, &dd](double ta) ->double { return 1.0/(1.0 + exp(-2.0*(ta-tt0)/dd));};
		 return kernel_K(t,x)*(1.0 -tth0(t*x));};

  gsl_function_pp<decltype(der_K)> dK(der_K);

  gsl_function *G = static_cast<gsl_function*>(&dK);

  double result, err;
  gsl_deriv_central(G, MV, MV*0.01, &result, &err);

  if(err/result > 1e-2) crash("Failed to evaluate kernel  derivative");

  cout<<"(t, MV): ("<<t<<","<<MV<<")"<<endl;
  cout<<"ker_SD(t,MV):"<<der_K(MV)<<endl;
  cout<<"ker_SD_derivative: "<<result<<" +- "<<err<<endl;

  return result;
  
  

}

double der_kernel_K(double t, double MV) {

  auto der_K = [&t](double x) {
		 
		 return kernel_K(t,x);};

  gsl_function_pp<decltype(der_K)> dK(der_K);

  gsl_function *G = static_cast<gsl_function*>(&dK);

  double result, err;
  gsl_deriv_central(G, MV, MV*0.01, &result, &err);

  if(err/result > 1e-2) crash("Failed to evaluate kernel  derivative");

  cout<<"(t, MV): ("<<t<<","<<MV<<")"<<endl;
  cout<<"ker_tot(t,MV):"<<der_K(MV)<<endl;
  cout<<"ker_tot_derivative: "<<result<<" +- "<<err<<endl;

  return result;
  
  

}


double Kernel_Pi_q2(double t, double Q, double a) {

  
  return t*t - 4.0*pow(sin(Q*t/2.0),2)/pow(Q,2);

}

void Plot_kernel_K(int npoints) {

  double tmin=0.0;
  double tmax= 50*5.0; //time is in fermi

  Vfloat pts;
  Vfloat val;

  for(int i=0; i<=npoints;i++) {
    double pt = tmin + i*(tmax-tmin)/npoints;
    pts.push_back(pt);
    val.push_back( kernel_K( pt/0.197327, 1.0));
  }

  Print_To_File({}, {pts, val}, "../data/gm2/light/kernel.dat", "", "");

  
  return;

}


void Plot_Energy_windows_K() {

  //open Silvano's file

  //tau=m_\mu t    t(fm)       M(SD)       M(W)        M(LD)      x=E/m_\mu   E(GeV)     Mtilde(SD)  Mtilde(W)  Mtilde(LD)

  Vfloat t_adim = Read_From_File("../energy_window_ker/modulating.dat", 0, 10);
  Vfloat t_fm = Read_From_File("../energy_window_ker/modulating.dat", 1, 10);
  Vfloat M_SD_t = Read_From_File("../energy_window_ker/modulating.dat", 2, 10);
  Vfloat M_W_t = Read_From_File("../energy_window_ker/modulating.dat", 3, 10);
  Vfloat M_LD_t = Read_From_File("../energy_window_ker/modulating.dat", 4, 10);
  Vfloat E_adim= Read_From_File("../energy_window_ker/modulating.dat", 5, 10);
  Vfloat E_GeV= Read_From_File("../energy_window_ker/modulating.dat", 6, 10);
  Vfloat M_SD_E = Read_From_File("../energy_window_ker/modulating.dat", 7, 10);
  Vfloat M_W_E = Read_From_File("../energy_window_ker/modulating.dat", 8, 10);
  Vfloat M_LD_E = Read_From_File("../energy_window_ker/modulating.dat", 9, 10);


  //multiply M by kernel functions

  auto F_t = [&](double t) -> double {


	    
	       auto F = [&t](double x) -> double {
			  
			  double z= t*m_muon;
			  
			  double y= 0.5*z*x/sqrt(1.0-x);
			  
			  if( z < 8e-6) { cout<<"Warning: kernel function evaluated at too small value of the argument and has been set to zero"<<endl; return 0.0;}

			  return 2*(1.0-x)*(1.0- pow(sin(y)/y,2));

			};

	       gsl_function_pp<decltype(F)> Fp(F);
 	       gsl_integration_workspace * w = gsl_integration_workspace_alloc (10000);
	       gsl_function *G = static_cast<gsl_function*>(&Fp);
	       double res_GSL, err_GSL;
	       gsl_integration_qags(G, 0.0, 1.0, 0.0, 1e-5, 10000, w, &res_GSL, &err_GSL);
	       gsl_integration_workspace_free (w);

	       return res_GSL;
  
	     };


   auto F_E = [&](double E) -> double { 

	       auto F = [&E](double x) -> double {
			  
			  double x2= pow(E/m_muon,2);

			  double resc= pow(m_muon/E,3);
			  
			  return 3*resc*x2*(1.0-x)*x*x/(x*x + (1-x)*x2);

			};

	       gsl_function_pp<decltype(F)> Fp(F);
 	       gsl_integration_workspace * w = gsl_integration_workspace_alloc (10000);
	       gsl_function *G = static_cast<gsl_function*>(&Fp);
	       double res_GSL, err_GSL;
	       gsl_integration_qags(G, 0.0, 1.0, 0.0, 1e-5, 10000, w, &res_GSL, &err_GSL);
	       gsl_integration_workspace_free (w);

	       return res_GSL;
  
	     };


   Vfloat M_SD_t_N, M_W_t_N, M_LD_t_N;

   for(unsigned int tt=0; tt< t_fm.size();tt++) {

     double integrand= F_t( t_fm[tt]*1.0/0.197327);
     M_SD_t_N.push_back( M_SD_t[tt]*integrand);
     M_W_t_N.push_back( M_W_t[tt]*integrand);
     M_LD_t_N.push_back( M_LD_t[tt]*integrand);

   }

   Vfloat M_SD_E_N, M_W_E_N, M_LD_E_N;

   for(unsigned int ee=0; ee< E_GeV.size();ee++) {

     double integrand= F_E( E_GeV[ee]);
     M_SD_E_N.push_back( M_SD_E[ee]*integrand); 
     M_W_E_N.push_back( M_W_E[ee]*integrand);
     M_LD_E_N.push_back( M_LD_E[ee]*integrand);

   }
   

   Print_To_File({}, {t_adim, t_fm, M_SD_t_N, M_W_t_N, M_LD_t_N, E_adim, E_GeV, M_SD_E_N, M_W_E_N, M_LD_E_N} , "../energy_window_ker/modulating_bis.dat","", "//tau=m_\mu t    t(fm)       M(SD)       M(W)        M(LD)      x^2=(E/m_\mu)^2   E(GeV)     Mtilde(SD)  Mtilde(W)  Mtilde(LD)" );


  return;  
}


double Zeta_function_laplacian_Luscher(double z) {


  
  double inf = 1e-6;
  double tol_zeta= 1e-10;

  assert(z>= 0) ;

  //compute Zeta_function using Luscher parametrization

  double z2= pow(z,2);

  //l2 is cutoff
  int l2= z2 + 2;
  // if(degeneracy(l2) == 0) l2++;

  
  int thresh = l2+20;  //l2+ 10;

  double t0= z2>2.0?M_PI/z2:0.5*M_PI;

  //t0=1.0;
  

  //sum first terms
  double val=0;
  double sum0=0.0;
  for(int m=0; m < l2; m++) {
    if(degeneracy(m) > 0 ) {
      sum0 += (1.0/sqrt(4*M_PI))*degeneracy(m)*(1.0/(m-z2));
    }
  }


  auto F2 = [&](double t) -> double {

	      double sum_t=0.0;
	      double exp_t= exp(-t);
	      double exp_qt= exp(-pow(M_PI,2)/t);
	      for(int m=0; m<l2;m++) sum_t += (1.0/sqrt(4.0*M_PI))*degeneracy(m)*(1.0/pow(2.0*M_PI,3))*pow(exp_t,m);
	      
	      double heat_kernel_red= (exp(t*z2)-1.0)/(pow(4.0*M_PI,2)*pow(t,1.5)); //m=0 term
	      for(int m=1; m<= thresh;m++) {
		heat_kernel_red += ( exp(t*z2)*(1.0/sqrt(4*M_PI)))*degeneracy(m)*( 1.0/pow(4.0*M_PI*t,3.0/2.0))*pow(exp_qt,m);
	      }
	      
	      double res= pow(2.0*M_PI,3)*(heat_kernel_red -exp(t*z2 + log(sum_t)));

	     
	      if( isnan(res)) crash("In integral repr. Zeta function_1 t: "+to_string_with_precision(t,3)+" m: "+to_string_with_precision(z2, 3)+" , integrand is nan");
	    
	      return res;

	      
	    };
    
  auto F = [&](double t) -> double {
	     
	     double heat_kernel_red=0;
	     double exp_t= exp(-t);
	     for(int m=l2;m<=thresh; m++) {
	       heat_kernel_red += (1.0/pow(2.0*M_PI,3))*degeneracy(m)*pow(exp_t,m);
	     }
	     
	     
	     double res= pow(2.0*M_PI,3)*(1.0/sqrt(4*M_PI))*exp(t*z2 + log(heat_kernel_red));

	     
	     if( isnan(res)) crash("In integral repr. Zeta function_2 t: "+to_string_with_precision(t,3)+" m: "+to_string_with_precision(z2, 3)+" , integrand is nan");
	    
	     return res;
	   };


  double int_val_2, int_err_2, int_val_1, int_err_1;
  double int_val_3= 0;
  double int_err_3=0;

 
  gsl_function_pp<decltype(F2)> F_1(F2);
  gsl_integration_workspace * w_1 = gsl_integration_workspace_alloc (10000);
  gsl_function *G_1 = static_cast<gsl_function*>(&F_1);
  gsl_integration_qags(G_1, inf, t0, 0.0, 1e-6, 10000, w_1, &int_val_1, &int_err_1);
  gsl_integration_workspace_free (w_1);


  gsl_function_pp<decltype(F2)> F_3(F2);
  gsl_integration_workspace * w_3 = gsl_integration_workspace_alloc (100000);
  gsl_function *G_3 = static_cast<gsl_function*>(&F_3);
  gsl_integration_qags(G_3, 0.0, inf, 0.0, 3e-4, 100000, w_3, &int_val_3, &int_err_3);
  gsl_integration_workspace_free (w_3);


  

  gsl_function_pp<decltype(F)> F_2(F);
  gsl_integration_workspace * w_2 = gsl_integration_workspace_alloc (10000);
  gsl_function *G_2 = static_cast<gsl_function*>(&F_2);
  gsl_integration_qags(G_2, t0, 40.0, 0.0, 1e-6, 10000, w_2, &int_val_2, &int_err_2);
  gsl_integration_workspace_free (w_2);

  if( (int_err_1+ int_err_3)/(fabs(int_val_1 + int_val_3)) > 1e-5) {

    gsl_function_pp<decltype(F2)> F_4(F2);
    gsl_integration_workspace * w_4 = gsl_integration_workspace_alloc (100000);
    gsl_function *G_4 = static_cast<gsl_function*>(&F_4);
    gsl_integration_qags(G_4, 0.0, t0, 0.0, 1e-6, 100000, w_4, &int_val_3, &int_err_3);
    gsl_integration_workspace_free (w_4);
    int_val_1 =0.0;
    int_err_1 =0.0;
    if( (int_err_1+ int_err_3)/(fabs(int_val_1 + int_val_3)) > 5e-6) {
    cout<<"val1 : "<<int_val_1<<" +- "<<int_err_1<<endl;
    cout<<"val3 : "<<int_val_3<<" +- "<<int_err_3<<endl;
    crash("In Z_function_laplacian_luscher unable to get target accuracy of 1e-5. Current precision: "+to_string_with_precision( (int_err_1+ int_err_3)/(fabs(int_val_1)+fabs(int_val_3)),10));
    }
  }
  if( int_err_2/fabs(int_val_2) > 5e-6) crash("In Z_function_laplacian_luscher unable to get target accuracy of 5e-6. Current precision: "+to_string_with_precision( int_err_2/fabs(int_val_2),10));
  //double int_val_2= boost::math::quadrature::gauss_kronrod<double, 15>::integrate(F, t0, 30.0 , 5, tol_zeta);
  //double int_val_1= boost::math::quadrature::gauss_kronrod<double, 15>::integrate(F2, inf, t0 , 5, tol_zeta);


  //cout<<"z: "<<z<<" [0,inf] : "<<int_val_3<<" [inf,t0] : "<<int_val_1<<endl;

  double int_val = int_val_3+ int_val_2 + int_val_1 + sum0 -M_PI/sqrt(t0);

  if( isnan( int_val)) crash("integral in Generalized zeta is nan");
  return int_val;
}


double tan_phi(double z) {

  int n2= (int)pow(z,2); //check whether is integer.
  int n2_p= n2+1;
  if(  (fabs( pow(z,2) - n2) < 1e-18 && degeneracy(n2*n2) > 0) ||  (fabs(pow(z,2) - n2_p) < 1e-18 && degeneracy(n2_p*n2_p) > 0) ) return 0.0;

  return -1.0*pow(M_PI, 1.5)*z/Zeta_function_laplacian_Luscher(z);

}

double tan_phi_der(double z) {

  auto tf = [&](double t) { return tan_phi(t);};

  //double x=z;

  double result, err;

  gsl_function_pp<decltype(tf)> Fp(tf);

  gsl_function *G = static_cast<gsl_function*>(&Fp);
  if(z < 2e-3) {
     gsl_deriv_forward(G, z, 1e-3, &result, &err);
  }
  else {
  gsl_deriv_central(G, z, 1e-3, &result, &err);
  }

  //double result= boost::math::differentiation::finite_difference_derivative(tf,x, &err);

  if(err/result > 1.0e-2) crash("numerical derivative in tan_phi'(z) has reached a low level of accuracy in phi_der. z: "+to_string_with_precision(z,10)+" rel_err: "+to_string_with_precision(err/result,10));


  if(isnan(result)) crash("Function tan_phi_der(z) returned nan");

  return result;

}

double phi_der(double z) {

  auto tf = [&](double t) { return tan_phi(t);};

  //double x=z;

  double result, err;

  gsl_function_pp<decltype(tf)> Fp(tf);

  gsl_function *G = static_cast<gsl_function*>(&Fp);
  if(z < 2e-3) {
     gsl_deriv_forward(G, z, 1e-4, &result, &err);
  }
  else {
  gsl_deriv_central(G, z, 1e-4, &result, &err);
  }

  //double result= boost::math::differentiation::finite_difference_derivative(tf,x, &err);

  if(err/result > 2.0e-2) crash("numerical derivative in phi'(z) has reached a low level of accuracy in phi_der. z: "+to_string_with_precision(z,10)+" rel_err: "+to_string_with_precision(err/result,10));

  double res= result/(1.0+ pow(tan_phi(z),2));

  if(isnan(res)) crash("Function phi_der(z) returned nan");

  return res;

}

double phi_der_for_back(double z, int mode) {

  auto tf = [&](double t) { return phi(t);};

  //double x=z;

  double result, err;

  gsl_function_pp<decltype(tf)> Fp(tf);

  gsl_function *G = static_cast<gsl_function*>(&Fp);
  if(mode==1)  gsl_deriv_forward(G, z, 1e-3, &result, &err);
  else if(mode==-1)  gsl_deriv_backward(G, z, 1e-3, &result, &err);
  else crash("mode: "+to_string( mode)+" not yet implemented");

  //double result= boost::math::differentiation::finite_difference_derivative(tf,x, &err);

  if(err/result > 1.0e-2) crash("numerical for_back derivative in phi'(z) has reached a low level of accuracy in phi_der. z: "+to_string_with_precision(z,10)+" rel_err: "+to_string_with_precision(err/result,10));


  if(isnan(result)) crash("Function phi_der_for_back(z) returned nan");

  return result;

}

void Zeta_function_zeroes(int Nzeros, Vfloat &res) { //using Newton method

  res.clear();

  auto F = [&](double z) -> double {

	     return Zeta_function_laplacian_Luscher(sqrt(z));

	   };

  auto dF = [&](double z) -> double {

	      auto tf = [&](double t) {return Zeta_function_laplacian_Luscher(sqrt(t));};
	      double der, err_der;
	      gsl_function_pp<decltype(tf)> Fp2(tf);

	      gsl_function *G = static_cast<gsl_function*>(&Fp2);
	      gsl_deriv_forward(G, z, 1e-5, &der, &err_der);
	      return der;
	     

	   };

  auto fdF = [&](double z, double* f, double *df) -> void {
	       *f = F(z);
	       *df = dF(z);
	       return;
	     };


  double Precision= 1e-9;
  double offset = 1e-10;

  int N2old=0;
  
  for(int izero=0; izero<Nzeros;izero++) {

    int status;   
    int iter = 0, max_iter = 1000;
    const gsl_root_fdfsolver_type *TT;
    gsl_root_fdfsolver *ss;
    double x0,x;

    if(izero==0) x= offset;
    else {
      int temp=N2old+1;
      while(degeneracy(temp) == 0) temp++;
      x= temp + offset;
      N2old=temp;
    }

    
    gsl_function_fdf_pp<decltype(F), decltype(dF), decltype(fdF)> Fp(F, dF, fdF);
    
    gsl_function_fdf *FDF = static_cast<gsl_function_fdf*>(&Fp);

    

    
    TT = gsl_root_fdfsolver_newton;
    ss = gsl_root_fdfsolver_alloc (TT);
    gsl_root_fdfsolver_set(ss, FDF, x);

    
    do
    {
      iter++;
      status = gsl_root_fdfsolver_iterate(ss);
      x0 = x;
      x = gsl_root_fdfsolver_root (ss);
      status = gsl_root_test_delta (x, x0, 0, Precision);

           
    }
    while (status == GSL_CONTINUE && iter < max_iter);

    if(status != GSL_SUCCESS) {

      cout<<"gsl_root_fdf_solver_newton was unable to find zero of the Luscher zeta function for  izero: "<<izero<<endl;
      crash("Aborting...");
    }

    gsl_root_fdfsolver_free(ss);

    res.push_back(x);
  }

  

  return;
}

double phi(double z) {

  int n2= (int)pow(z,2); //check whether is integer.
  int n2_p= n2+1;
  if( (fabs( pow(z,2) - n2) < 1e-16 && degeneracy(n2*n2) > 0) || ( fabs(pow(z,2) - n2_p) < 1e-16 && degeneracy(n2_p*n2_p) > 0) ) return 0.0;
  double sum= Zeta_function_laplacian_Luscher(z);

  double res= atan( -1.0*pow(M_PI,1.5)*z/sum);

  if(isnan(res)) crash("Function phi(z) returned nan");
  
  return res;

}

void Generate_free_corr_data() {

  //generate pairs

  vector<pair<double, int>> amu_Tmax;

  //amu_Tmax.push_back( make_pair( 0.00060, 112));

  amu_Tmax.push_back( make_pair( 0.011, 112));
  amu_Tmax.push_back( make_pair( 0.012, 112));
  amu_Tmax.push_back( make_pair( 0.016, 112));
  amu_Tmax.push_back( make_pair( 0.018, 112));
  
  
  //amu_Tmax.push_back( make_pair( 0.00044, 112));

  //amu_Tmax.push_back( make_pair( 0.00077, 64));

  //amu_Tmax.push_back( make_pair(0.018,64));
  //amu_Tmax.push_back( make_pair(0.020,64));
  //amu_Tmax.push_back( make_pair(0.018,96));
  //amu_Tmax.push_back( make_pair(0.020,96));
  //amu_Tmax.push_back( make_pair(0.013,96));
  

  /*
  //charm

  // A ensembles
  //amu_Tmax.push_back( make_pair( 0.240, 24));
  //amu_Tmax.push_back( make_pair( 0.240, 32));
  //amu_Tmax.push_back( make_pair( 0.240, 48));
  amu_Tmax.push_back( make_pair( 0.265, 24));
  amu_Tmax.push_back( make_pair( 0.265, 32));
  amu_Tmax.push_back( make_pair( 0.265, 48));
  amu_Tmax.push_back( make_pair( 0.290, 24));
  amu_Tmax.push_back( make_pair( 0.290, 32));
  amu_Tmax.push_back( make_pair( 0.290, 48));
  amu_Tmax.push_back( make_pair( 0.300, 24));
  amu_Tmax.push_back( make_pair( 0.300, 32));
  amu_Tmax.push_back( make_pair( 0.300, 48));

  //B ensembles
  amu_Tmax.push_back( make_pair( 0.21, 48));
  amu_Tmax.push_back( make_pair( 0.21, 64));
  amu_Tmax.push_back( make_pair( 0.23, 48));
  amu_Tmax.push_back( make_pair( 0.23, 64));
  amu_Tmax.push_back( make_pair( 0.25, 48));
  amu_Tmax.push_back( make_pair( 0.25, 64));

  //C ensembles
  amu_Tmax.push_back( make_pair( 0.175, 80));
  amu_Tmax.push_back( make_pair( 0.195, 80));
  amu_Tmax.push_back( make_pair( 0.215, 80));


  //D ensembles
  amu_Tmax.push_back( make_pair( 0.165, 96));
  amu_Tmax.push_back( make_pair( 0.175, 96));


  //strange

  //B ensembles
  
  amu_Tmax.push_back( make_pair(0.019, 64));
  amu_Tmax.push_back( make_pair(0.019, 96));
  amu_Tmax.push_back( make_pair(0.021, 64));
  amu_Tmax.push_back( make_pair(0.021, 96));
  //C ensembles
  amu_Tmax.push_back( make_pair(0.016, 80));
  amu_Tmax.push_back( make_pair(0.018, 80));
  //D ensembles
  amu_Tmax.push_back( make_pair( 0.014, 96));
  amu_Tmax.push_back( make_pair( 0.015, 96));


  //A ensembles
  /*amu_Tmax.push_back( make_pair( 2*0.0205, 24));
  amu_Tmax.push_back( make_pair( 2*0.023, 24));
  amu_Tmax.push_back( make_pair( 2*0.0205, 32));
  amu_Tmax.push_back( make_pair( 2*0.023, 32));
  amu_Tmax.push_back( make_pair( 0.0205, 24));
  amu_Tmax.push_back( make_pair( 0.023, 24));
  amu_Tmax.push_back( make_pair( 0.0205, 32));
  amu_Tmax.push_back( make_pair( 0.023, 32));
  amu_Tmax.push_back(make_pair( 2*0.300, 24));
  amu_Tmax.push_back(make_pair( 0.300, 24));
  amu_Tmax.push_back(make_pair( 2*0.300, 32));
  amu_Tmax.push_back(make_pair( 0.300, 32));
  

  //light
  //amu_Tmax.push_back( make_pair( 2*0.00054, 96));
  //amu_Tmax.push_back( make_pair( 0.00054, 96));
  //amu_Tmax.push_back( make_pair( 2*0.00060, 80));
  //amu_Tmax.push_back( make_pair( 0.00060, 80));
  //amu_Tmax.push_back( make_pair( 2*0.00072, 64));
  //amu_Tmax.push_back( make_pair( 0.00072, 64));
  //amu_Tmax.push_back( make_pair( 2*0.00072, 96));
  //amu_Tmax.push_back( make_pair( 0.00072, 96));

  //amu_Tmax.push_back( make_pair(0.0, 200));


  
  */

  
  
  

  for(auto &p: amu_Tmax) { cout<<"Executing amu: "<<p.first<<" Tmax: "<<p.second<<endl; Compute_free_corr(p.first, p.second);}

  
  return;
}


void Get_spec_dens_free(const Vfloat &ams, string out_path) {

 
  for( auto &am: ams)  {
    Compute_free_spectral_density(3, am, 1, 0.01, out_path);
    Compute_free_spectral_density(3, am, -1, 0.01, out_path);
    cout<<"am: "<<am<<" computed!"<<endl;
  }
  
  return;
}


void Compute_free_corr(double am, int Tmax) {

 
  double Nc=3; //three colors

  auto C_cont = [&Nc, &am](int t) -> double {

		  //double tolerance=1e-16;
		  //double err;

		  auto f = [&am, &t, &Nc](double x) {  return (Nc*2.0/pow(M_PI,2))*exp(-2.0*t*sqrt( pow(x,2) + pow(am,2)))*pow(x,2)*( 1.0/3 + pow(am,2)/( 6.0*( pow(am,2) + pow(x,2))));};

		  double val;
		  double tolerance=1e-9;
		  double err;
		  


		  gsl_function_pp<decltype(f)> F_corr(f);
		  gsl_integration_workspace * w = gsl_integration_workspace_alloc (10000);
		  gsl_function *G = static_cast<gsl_function*>(&F_corr);
		  gsl_integration_qagiu(G, 0.0, 0.0, tolerance, 10000, w, &val, &err);
		  gsl_integration_workspace_free (w);

		  if( err/fabs(val) > 5*tolerance) crash("In free_vector_corr_cont gls integration not able to achieve target precision");
		  
		  return val;

		  //return boost::math::quadrature::gauss_kronrod<double, 15>::integrate( f, 0, numeric_limits<double>::infinity(), 5,tolerance, &err);
		  
		};

  auto ptm2 = [](double p1,double p2, double p3) { return pow(sin(p1),2)+ pow(sin(p2),2)+ pow(sin(p3),2);};
  auto phm2 = [](double p1,double p2, double p3) { return pow(2*sin(p1/2),2) + pow(2*sin(p2/2),2) + pow(2*sin(p3/2),2);};
  auto fa = [&phm2](double p1, double p2, double p3) { return 1.0 + 0.5*phm2(p1,p2,p3);};
  auto fb = [&phm2](double p1, double p2, double p3) { return phm2(p1,p2,p3) + 0.5*( pow(2*sin(p1/2),2)*(pow(2*sin(p2/2),2)+pow(2*sin(p3/2),2)) + pow(2*sin(p2/2),2)*pow(2*sin(p3/2),2));};
  auto D2 = [&am, &fa, &fb](double p1,double p2, double p3) { return ( fb(p1,p2,p3) + pow(am,2))*( 4*fa(p1,p2,p3) + fb(p1,p2,p3) + pow(am,2));};
  auto W = [&am, &phm2, &fa, &fb](double p1, double p2, double p3) { return 0.25*pow( phm2(p1,p2,p3) - (fb(p1,p2,p3) + pow(am,2))/fa(p1,p2,p3)  ,2);};
  auto shaEp2 = [&am, &fa, &fb](double p1,double p2, double p3) { return pow( (fb(p1,p2,p3) + pow(am,2))/(2*fa(p1,p2,p3)),2) + (fb(p1,p2,p3)+pow(am,2))/fa(p1,p2,p3);};
  auto shaEp =[&shaEp2](double p1,double p2, double p3) { return sqrt(shaEp2(p1,p2,p3));};
  auto Ep =[&shaEp](double p1,double p2, double p3) {return asinh(shaEp(p1,p2,p3));};

  auto corr = [&am, &Nc, &ptm2, &D2, &W, &shaEp2, &Ep](double p1,double p2,double p3, double t, int r) -> double {

		return (32.0*Nc/(pow(2.0*M_PI,3)))*exp(-2*Ep(p1,p2,p3)*t)*( shaEp2(p1,p2,p3) +(1.0/3.0)*ptm2(p1,p2,p3) +pow(am,2)+ r*W(p1,p2,p3))/D2(p1,p2,p3);
	      };

  vector<double> corr_pert_res_tm(2*Tmax,0.0);
  vector<double> corr_pert_res_OS(2*Tmax, 0.0);

  for(int i=0;i<2;i++) { //loop over r

    int r= 2*i-1; //set Wilson parameter

   

    int time;
    double tol = 1e-9;
    
    for(int t=1;t<=Tmax;t++) { //loop over time

      cout<<"r: "<<r<<" t: "<<t<<endl;
      time =t;

      double err;
    
      auto corrp1= [&corr, &time, &r, &tol](double p1) -> double {  //boost only performs 1d integrals. 
     
		     auto corrp2 = [&corr, &p1, &time, &r, &tol](double p2) {

				     auto corrp3 = [&corr, &p1, &p2, &time, &r](double p3) {
						     return corr(p1, p2, p3, time,r);
						   };

				     
				     double err_3;
				     double tol3= 1e-9;
				     double val_3;
				     //double val_3 = boost::math::quadrature::gauss_kronrod<double, 15>::integrate( corrp3, 0, M_PI, 5,tol, &err_3);
				     gsl_function_pp<decltype(corrp3)> F_corr3(corrp3);
				     gsl_integration_workspace * w3 = gsl_integration_workspace_alloc(10000);
				     gsl_function *G3 = static_cast<gsl_function*>(&F_corr3);
				     gsl_integration_qags(G3, 0.0, M_PI, 0.0, tol3, 10000, w3, &val_3, &err_3);
				     gsl_integration_workspace_free (w3);

				     if( err_3/fabs(val_3) > 2*tol3) crash("corr_p3 did not achieve the target accuracy of "+to_string_with_precision(2*tol3, 10)+". Current precision: "+to_string_with_precision( err_3/fabs(val_3), 10));
				     return val_3;
				   };
		     
		     double err_2;
		     double tol2=1e-9;
		     double val_2;
		     //double val_2 = boost::math::quadrature::gauss_kronrod<double, 15>::integrate( corrp2, 0, M_PI, 5,tol, &err_2);
		     gsl_function_pp<decltype(corrp2)> F_corr2(corrp2);
		     gsl_integration_workspace * w2 = gsl_integration_workspace_alloc(10000);
		     gsl_function *G2 = static_cast<gsl_function*>(&F_corr2);
		     gsl_integration_qags(G2, 0.0,M_PI, 0.0, tol2, 10000, w2, &val_2, &err_2);
		     gsl_integration_workspace_free (w2);
		     
		     if( err_2/fabs(val_2) > 2*tol2) crash("corr_p2 did not achieve the target accuracy of "+to_string_with_precision(2*tol2, 10)+". Current precision: "+to_string_with_precision( err_2/fabs(val_2), 10));
		     return val_2;
		   };

      double val;
      gsl_function_pp<decltype(corrp1)> F_corr1(corrp1);
      gsl_integration_workspace * w1 = gsl_integration_workspace_alloc(10000);
      gsl_function *G1 = static_cast<gsl_function*>(&F_corr1);
      gsl_integration_qags(G1, 0.0,M_PI, 0.0, tol, 10000, w1, &val, &err);
      gsl_integration_workspace_free (w1);

      if( err/fabs(val) > 2*tol) crash("corr_p1 did not achieve the target accuracy of "+to_string_with_precision(2*tol, 10)+". Current precision: "+to_string_with_precision( err/fabs(val), 10));
      //double val = boost::math::quadrature::gauss_kronrod<double, 15>::integrate(corrp1, 0, M_PI, 5, tol, &err);
           
  
      if(r==1)  corr_pert_res_OS[t] = val;
      else if(r==-1) corr_pert_res_tm[t] = val;
      else crash("Wilson parameter r: "+to_string(r)+" not recognized");

    }

  }


  //compute correlator in the continuum limit
  vector<double> corr_pert_cont(2*Tmax,0.0);
  for(int t=1;t<=Tmax;t++) corr_pert_cont[t] = C_cont(t);


  //take difference between corr_pert_cont and corr_pert_res_tm(OS)
  vector<double> a2corr_pert_res_tm(2*Tmax,0.0);
  vector<double> a2corr_pert_res_OS(2*Tmax,0.0);
  
  for(int t=0;t<2*Tmax;t++) {
    a2corr_pert_res_tm[t] = corr_pert_cont[t] - corr_pert_res_tm[t];
    a2corr_pert_res_OS[t] = corr_pert_cont[t] - corr_pert_res_OS[t];
  }
 
  //Print the result

  //create directory
  boost::filesystem::create_directory("../Vkvk_cont");
  boost::filesystem::create_directory("../Vkvk_cont/a_mu_SD_scaling");
  boost::filesystem::create_directory("../Vkvk_cont/"+to_string(Tmax)+"_m"+to_string_with_precision(am,5));


  Print_To_File({}, {a2corr_pert_res_tm, corr_pert_res_tm, corr_pert_cont}, "../Vkvk_cont/"+to_string(Tmax)+"_m"+to_string_with_precision(am,5)+"/OPPOR", "" , "");
  Print_To_File({}, {a2corr_pert_res_OS, corr_pert_res_OS, corr_pert_cont}, "../Vkvk_cont/"+to_string(Tmax)+"_m"+to_string_with_precision(am,5)+"/SAMER", "" , "");
  

  Vfloat amu_SD_OPPOR_list, amu_SD_OPPOR_diff_list, amu_SD_SAMER_list, amu_SD_SAMER_diff_list, a_lat_list;

  double a_coarse= 0.0908026;
  double a_finest= 0.0100000;
  int Niter=100;

  //get amu_SD_scaling

    double tt0 = 0.4/0.197327;
    double tDelta= 0.15/0.197327;
    double aem= 1.0/137.035999;
    

    for(int it=0; it < Niter;it++) {

      double ag= (a_finest+ it*(a_coarse-a_finest)/(Niter-1.0));
      a_lat_list.push_back(ag);
      ag /= 0.197327; //conversion in Gev^-1

      double amu_SD_OPPOR_diff=0.0;
      double amu_SD_SAMER_diff=0.0;
      double amu_SD_OPPOR=0.0;
      double amu_SD_SAMER=0.0;

      for(int tt=1; tt< Tmax; tt++) {

	amu_SD_OPPOR_diff += 4.0*pow(aem,2)*( pow(2.0/3.0,2) + pow(1.0/3.0, 2))*a2corr_pert_res_tm[tt]*(1.0 -  1.0/(1.0 + exp(-2.0*(tt*ag-tt0)/tDelta)))*kernel_K(tt, ag);
	amu_SD_SAMER_diff += 4.0*pow(aem,2)*( pow(2.0/3.0,2) + pow(1.0/3.0, 2))*a2corr_pert_res_OS[tt]*(1.0 -  1.0/(1.0 + exp(-2.0*(tt*ag-tt0)/tDelta)))*kernel_K(tt, ag);
	amu_SD_OPPOR +=  4.0*pow(aem,2)*( pow(2.0/3.0,2) + pow(1.0/3.0, 2))*corr_pert_res_tm[tt]*(1.0 -  1.0/(1.0 + exp(-2.0*(tt*ag-tt0)/tDelta)))*kernel_K(tt, ag);
	amu_SD_SAMER +=  4.0*pow(aem,2)*( pow(2.0/3.0,2) + pow(1.0/3.0, 2))*corr_pert_res_OS[tt]*(1.0 -  1.0/(1.0 + exp(-2.0*(tt*ag-tt0)/tDelta)))*kernel_K(tt, ag);
      }

      amu_SD_OPPOR_diff_list.push_back(amu_SD_OPPOR_diff);
      amu_SD_SAMER_diff_list.push_back(amu_SD_SAMER_diff);
      amu_SD_OPPOR_list.push_back( amu_SD_OPPOR);
      amu_SD_SAMER_list.push_back( amu_SD_SAMER);
    }




    Print_To_File( {}, { a_lat_list, amu_SD_OPPOR_list, amu_SD_OPPOR_diff_list,  amu_SD_SAMER_list, amu_SD_SAMER_diff_list}, "../Vkvk_cont/a_mu_SD_scaling/scaling_am_"+to_string_with_precision(am,5)+".dat", "", "#a[fm] a_mu[tm]  Da_mu[tm] a_mu[OS]   Da_mu[OS]");
  


  return;
}

void Compute_free_spectral_density(int Nc, double am, int reg, double step_size_erg, string dir_out) {


  boost::filesystem::create_directory("../data/"+dir_out);
  boost::filesystem::create_directory("../data/"+dir_out+"/spec_dens_free");
  boost::filesystem::create_directory("../data/"+dir_out+"/spec_dens_free/tm");
  boost::filesystem::create_directory("../data/"+dir_out+"/spec_dens_free/OS");
 


  if( (reg != -1) && (reg != 1)) crash("Compute_free_spectral_density called with an unknown regularization reg: "+to_string(reg));
  string tag_reg;
  if(reg==-1) tag_reg= "tm";
  else tag_reg="OS";

  //assume that energy max is smaller than 2*M_PI in lattice units;
  //divide energies in step_size of 0.01
  int Erg_size = 2*M_PI/step_size_erg;
  Vfloat spec_dens_hist(Erg_size, 0.0);
  Vfloat energy_list, spec_dens_cont;
  for(int i=0; i<(signed)spec_dens_hist.size(); i++) energy_list.push_back( (i+0.5)*step_size_erg);
			
  auto ptm2 = [](double p1,double p2, double p3) { return pow(sin(p1),2)+ pow(sin(p2),2)+ pow(sin(p3),2);};
  auto phm2 = [](double p1,double p2, double p3) { return pow(2*sin(p1/2),2) + pow(2*sin(p2/2),2) + pow(2*sin(p3/2),2);};
  auto fa = [&phm2](double p1, double p2, double p3) { return 1.0 + 0.5*phm2(p1,p2,p3);};
  auto fb = [&phm2](double p1, double p2, double p3) { return phm2(p1,p2,p3) + 0.5*( pow(2*sin(p1/2),2)*(pow(2*sin(p2/2),2)+pow(2*sin(p3/2),2)) + pow(2*sin(p2/2),2)*pow(2*sin(p3/2),2));};
  auto D2 = [&am, &fa, &fb](double p1,double p2, double p3) { return ( fb(p1,p2,p3) + pow(am,2))*( 4*fa(p1,p2,p3) + fb(p1,p2,p3) + pow(am,2));};
  auto W = [&am, &phm2, &fa, &fb](double p1, double p2, double p3) { return 0.25*pow( phm2(p1,p2,p3) - (fb(p1,p2,p3) + pow(am,2))/fa(p1,p2,p3)  ,2);};
  auto shaEp2 = [&am, &fa, &fb](double p1,double p2, double p3) { return pow( (fb(p1,p2,p3) + pow(am,2))/(2*fa(p1,p2,p3)),2) + (fb(p1,p2,p3)+pow(am,2))/fa(p1,p2,p3);};
  auto shaEp =[&shaEp2](double p1,double p2, double p3) { return sqrt(shaEp2(p1,p2,p3));};
  auto Ep =[&shaEp](double p1,double p2, double p3) {return asinh(shaEp(p1,p2,p3));};


 //lattice spectral density_int
 auto spec_int_lat = [&am, &Nc, &ptm2, &D2, &W, &shaEp2, &Ep, &reg](double p1,double p2,double p3) -> double {

		return (32.0*Nc/(pow(2.0*M_PI,3)))*(shaEp2(p1,p2,p3) +(1.0/3.0)*ptm2(p1,p2,p3) +pow(am,2)+ reg*W(p1,p2,p3))/D2(p1,p2,p3);
	      };


 auto spec_cont = [&am, &Nc](double E) {
		    double x=0;
		    if( E < 2.0*am) return 0.0;
		    else x= sqrt( pow(E/2.0,2) - pow(am,2));

		    double jaco= 2.0*x/sqrt( pow(x,2) + pow(am,2));
		    return (1.0/jaco)*(Nc*2.0/pow(M_PI,2))*pow(x,2)*( 1.0/3 + pow(am,2)/( 6.0*( pow(am,2) + pow(x,2))));
		  };

 int Npoints_dir=1000;
 double step_size_px = M_PI/(Npoints_dir);
 Vfloat pxs, pys, pzs;
 for(int i=0; i<Npoints_dir;i++) { pxs.push_back( i*M_PI/(Npoints_dir-1.0)); pys.push_back( i*M_PI/(Npoints_dir -1.0)); pzs.push_back( i*M_PI/(Npoints_dir-1.0));}

 #pragma omp parallel for reduction (+:spec_dens_hist)
 for( auto &px: pxs)
   for( auto &py: pys)
     for( auto &pz: pzs) {
       double Erg = 2*Ep(px,py,pz);
       double Ampl= pow(step_size_px,3)*spec_int_lat(px,py,pz)/step_size_erg;
       spec_dens_hist[ (int)(Erg/step_size_erg) ] += Ampl;
     }
 

 //boost::math::interpolators::cardinal_cubic_b_spline<double> f_interp_boost(spec_dens_hist.begin(), spec_dens_hist.end(), 0.5*step_size_erg, step_size_erg);

 //evaluate spectral density in the continuum
 for(auto &e: energy_list) spec_dens_cont.push_back( spec_cont(e));

 //print_to_file

 Print_To_File({}, {energy_list, spec_dens_hist, spec_dens_cont}, "../data/"+dir_out+"/spec_dens_free/"+tag_reg+"/am_"+to_string_with_precision(am,5), "", "");


}

double free_vector_corr_cont(int Nc, double am, double t) {

  
  double val;
  double tolerance=1e-9;
  double err;

  auto f = [&am, &t, &Nc](double x) {  return (Nc*2.0/pow(M_PI,2))*exp(-2.0*t*sqrt( pow(x,2) + pow(am,2)))*pow(x,2)*( 1.0/3 + pow(am,2)/( 6.0*( pow(am,2) + pow(x,2))));};


   gsl_function_pp<decltype(f)> F_corr(f);
   gsl_integration_workspace * w = gsl_integration_workspace_alloc (10000);
   gsl_function *G = static_cast<gsl_function*>(&F_corr);
   gsl_integration_qagiu(G, 0.0, 0.0, tolerance, 10000, w, &val, &err);
   gsl_integration_workspace_free (w);

   if( err/fabs(val) > 5*tolerance) crash("In free_vector_corr_cont gls integration not able to achieve target precision");

   return val;

   //return boost::math::quadrature::gauss_kronrod<double, 15>::integrate( f, 0, numeric_limits<double>::infinity(), 5,tolerance, &err);
		  


}



void Compute_SD_window_Free() {

  VVfloat C_samer;
  VVfloat C_opper;
  VVfloat C_cons;

   double t0= 0.4; //fm
   double a0=1.0;
   double Q= pow(2.0/3.0,2)+ pow(1.0/3.0,2);
   double Delta = 0.15; //fm
   double Alpha=1.0/137.04;
   auto SD_th = [&t0, &Delta](double t) ->double { return 1.0- 1.0/(1.0 + exp(-2.0*(t-t0)/Delta));};
   auto int_SD_cont = [&](double t) -> double { return 4.0*pow(Alpha,2)*(1.0/(2.0*M_PI*M_PI*t*t*t))*SD_th(t)*kernel_K(t,1.0);};

   double error_amu_cont;
   double eps_cont = 1.0e-10;
   double a_mu_SD_cont= boost::math::quadrature::gauss_kronrod<double,15>::integrate(int_SD_cont, 0.0, numeric_limits<double>::infinity() , 5, eps_cont, &error_amu_cont);

  for(int k=1;k<=15;k++) {
    Vfloat rr= Read_From_File("../theo_RM123/VkVk_free_data/corr_r1_1_r2_1_L_48_T_96_mu_0.00000_scale_"+to_string(k), 1, 2);
    Vfloat rmr= Read_From_File("../theo_RM123/VkVk_free_data/corr_r1_1_r2_0_L_48_T_96_mu_0.00000_scale_"+to_string(k), 1, 2);
    Vfloat rcons= Read_From_File("../theo_RM123/VkVk_free_data/corr_r1_1_r2_2_L_48_T_96_mu_0.00000_scale_"+to_string(k), 1, 2);
    C_samer.push_back(rr);
    C_opper.push_back(rmr);
    C_cons.push_back(rcons);
  }

  Vfloat a2_Corr_48_96_samer;
  Vfloat a2_Corr_48_96_opper;

  a2_Corr_48_96_samer= Read_From_File("../theo_RM123/VkVk_free_data/48_SAMER", 1, 2);
  a2_Corr_48_96_opper= Read_From_File("../theo_RM123/VkVk_free_data/48_OPPOR", 1, 2);

  Vfloat cont_lim_samer = Sum_vectors<double>(a2_Corr_48_96_samer, C_opper[0]);
  Vfloat cont_lim_opper = Sum_vectors<double>(a2_Corr_48_96_opper, C_samer[0]);


  for(int t=1; t< (signed)cont_lim_samer.size();t++) cont_lim_samer[t]*=2.0*M_PI*M_PI*t*t*t;
  for(int t=1; t< (signed)cont_lim_opper.size();t++) cont_lim_opper[t]*=2.0*M_PI*M_PI*t*t*t;



  Vfloat a2_coeff_samer, a2_coeff_opper, a2_coeff_cons;
  Vfloat cont_samer, cont_opper, cont_cons;

  //fit correlator for each t
  for(int r=0;r<3;r++) {
  for(int t=1;t<10;t++) {
    

    Vfloat X,Y,Z;
    for(int i=5;i<15;i++) {X.push_back( a0/(1.0+i)); Y.push_back( (1.0/a0/a0/a0)*C_samer[i][t]); Z.push_back( (1.0/a0/a0/a0)*C_opper[i][t]);}
    //glue Y and Z
    Vfloat V=Y;
    Vfloat Xd = X;
    Vfloat Y_cons;
    for(int i=5;i<15;i++) {Y_cons.push_back( (1.0/a0/a0/a0)*C_cons[i][t]);}
    V.insert(V.end(), Z.begin(), Z.end());
    Xd.insert(Xd.end(), Xd.begin(), Xd.end());
    T_fit Fit_f((r==2)?X:Xd, (r==2)?Y_cons:V);
    if(r==2) Fit_f.add_pars(3.0); //cont. lim
    Fit_f.add_pars(2.5); //a^2 samer
    if(r != 2) Fit_f.add_pars(2.5); //a^2 opper
    //Fit_f.add_pars(1.0); //a^4 samer
    //Fit_f.add_pars(1.0); //a^4 opper
  

    Fit_f.ansatz = [&](const Vfloat &ip, double x, int imeas) -> double {

		     

		     double tph = t*a0;
		     double K= 1.0/(2.0*M_PI*M_PI*pow(tph,3));
		     if(r==2) return (1.0/a0/a0/a0)*K*ip[0] + 1.0*K*(1.0/a0/a0/a0)*(x/tph)*(x/tph)*ip[1] + 0.0*K*(1.0/a0*a0*a0)*pow((x/tph),4);
		     else return (1.0/a0/a0/a0)*K*(imeas<Y.size()?cont_lim_samer[t]:cont_lim_opper[t])  + 1.0*(imeas<Y.size()?ip[0]:ip[1])*K*(x/tph)*(x/tph) + 0.0*K*pow((x/tph),4);
		   };

    fit_t_res fit_result = Fit_f.fit();
    if(r==0) {a2_coeff_samer.push_back(fit_result.pars[0]); cont_samer.push_back(fit_result.pars[0]);}
    else if(r==1) {a2_coeff_opper.push_back(fit_result.pars[1]); cont_opper.push_back(fit_result.pars[0]);}
    else {a2_coeff_cons.push_back( fit_result.pars[1]); cont_cons.push_back(fit_result.pars[0]);}

   

  }
  }

  cout<<"Printing a^2 fit coefficients (same r)"<<endl;
  printV(a2_coeff_samer, "", 1);
  cout<<"Printing a^2 fit coefficients (opposite r)"<<endl;
  printV(a2_coeff_opper, "", 1);
  cout<<"Printing a^2 fit coefficients conserved"<<endl;
  printV(a2_coeff_cons, "", 1);
  cout<<"Printing cont fit coefficients (same r)"<<endl;
  printV(cont_samer, "", 1);
  cout<<"Printing cont fit coefficients (opposite r)"<<endl;
  printV(cont_opper, "", 1);
  cout<<"Printing cont fit coefficients conserved"<<endl;
  printV(cont_cons, "", 1);
  

  

  
  Vfloat alist, aSD_samer_list, aSD_opper_list;
  Vfloat Clog_SD_samer_list, Clog_SD_opper_list;
  //Compute_SD_window
  for(int k=0;k<15;k++) {
    double aSD_samer=0.0;
    double aSD_opper=0.0;
    double a= a0/(1.0+k);
    for(int t=1; t< 30;t++) { aSD_samer += a0*(1.0/pow(a0,3))*4.0*pow(Alpha,2)*SD_th(t*a0)*C_samer[k][t]*kernel_K(t*a0,1.0);}
    for(int t=1; t< 30;t++) { aSD_opper += a0*(1.0/pow(a0,3))*4.0*pow(Alpha,2)*SD_th(t*a0)*C_opper[k][t]*kernel_K(t*a0,1.0);}
    alist.push_back(a);
    aSD_samer_list.push_back(aSD_samer);
    aSD_opper_list.push_back(aSD_opper);
    Clog_SD_samer_list.push_back(aSD_samer/(a*a*log(1.0/a)));
    Clog_SD_opper_list.push_back(aSD_opper/(a*a*log(1.0/a)));
      

  }

  
  
  Vfloat cont_val(15,a_mu_SD_cont);
 


   //fit SD window
  Vfloat a2log_samer, a2log_opper, SD_cont_fit_samer, SD_cont_fit_opper;
  for(int r=0;r<2;r++) {
    

    Vfloat X,Y,Z;
    for(int i=3;i<15;i++) {X.push_back( a0/(1.0+i)); Y.push_back( aSD_samer_list[i]); Z.push_back( aSD_opper_list[i]);}
    T_fit Fit_f_SD(X, r==0?Y:Z);
    Fit_f_SD.add_pars(1.3); //cont. lim
    Fit_f_SD.add_pars(0.5); //a^2*log(1/a)
    Fit_f_SD.add_pars(1.0); //a^2

    Fit_f_SD.ansatz = [&](const Vfloat &ip, double x, int imeas) -> double {

			double K= cont_val[0];
			double xd = x/a0;
			return 1.0*ip[0]*K + 1.0*ip[1]*xd*xd*log(1/xd) + ip[2]*xd*xd;
		   };

    fit_t_res fit_result_SD = Fit_f_SD.fit();
    if(r==0) {a2log_samer.push_back(fit_result_SD.pars[1]); SD_cont_fit_samer.push_back(fit_result_SD.pars[0]);}
    else {a2log_opper.push_back(fit_result_SD.pars[1]); SD_cont_fit_opper.push_back(fit_result_SD.pars[0]);}
  }

  //print
  cout<<"Printing a^2*log(1/a) fit coefficients (same r)"<<endl;
  printV(a2log_samer, "", 1);
  cout<<"Printing a^2*log(1/a) fit coefficients (opposite r)"<<endl;
  printV(a2log_opper, "", 1);
  cout<<"Printing cont fit SD (same r)"<<endl;
  printV(SD_cont_fit_samer, "", 1);
  cout<<"Printing cont fit SD (opposite r)"<<endl;
  printV(SD_cont_fit_opper, "", 1);

  VVfloat C_samer_swap= C_samer;
  VVfloat C_opper_swap= C_opper;
  Transpose_VV<double>(C_samer_swap);
  Transpose_VV<double>(C_opper_swap);

  for(int t=1;t<C_samer[0].size()/2;t++) { Print_To_File( {}, {alist, C_samer_swap[t], C_opper_swap[t]}, "../theo_RM123/VkVk_free_data/Corr_"+to_string(t)+".dat", "", "");}


  Print_To_File({}, {alist, aSD_samer_list, aSD_opper_list, cont_val }, "../theo_RM123/VkVk_free_data/SD.out", "", "# a   SD(samer)   SD(opper)   Clog(samer)   Clog(opper)");

  return;

};



double LL_functions::Vdual(double t, double m_rho, double Edual, double Rdual) {


  return (5.0/(18.0*pow(M_PI,2)))*(Rdual/pow(t,3))*exp(-(m_rho+Edual)*t)*( 1.0 + (m_rho+Edual)*t + 0.5*pow( (m_rho+Edual)*t,2));

}


double LL_functions::Gamma_rpp(double omega, double g_rho_pipi, double Mpi) {

  double k= sqrt( pow(omega,2)/4.0 - pow(Mpi,2));
  return pow(g_rho_pipi,2)*pow(k,3)/(6.0*M_PI*pow(omega,2));  
}

double LL_functions::Gamma_rpp_der(double omega, double g_rho_pipi, double Mpi) {

  double k= sqrt( pow(omega,2)/4.0 - pow(Mpi,2));

  double dk_domega = omega/(4.0*k);

  return (pow(g_rho_pipi,2)/(6.0*M_PI))*( -2.0*pow(k,3)/pow(omega,3) + 3.0*pow(k,2)*dk_domega/pow(omega,2));
  //return 2*k*( (8*pow(Mpi,2) + pow(omega,2))/(4*pow(omega,3)))*pow(g_rho_pipi,2)/(6.0*M_PI);
  
}

double LL_functions::h(double omega, double g_rho_pipi, double Mpi) {
  double k= sqrt( pow(omega,2)/4.0 - pow(Mpi,2));
  return pow(g_rho_pipi,2)*pow(k,3)*(2.0/(omega*6.0*pow(M_PI,2)))*log( (omega + 2*k)/(2.0*Mpi));

}

double LL_functions::h_prime(double omega, double g_rho_pipi, double Mpi) {

  double k= sqrt( pow(omega,2)/4.0 - pow(Mpi,2));
  return pow(g_rho_pipi,2)*pow(k,2)*(1.0/(6.0*pow(M_PI,2)*omega))*( 1.0 + (1.0+ 2.0*pow(Mpi/omega,2))*(omega/k)*log( (omega+2*k)/(2.0*Mpi)));

}

double LL_functions::h_second(double omega, double g_rho_pipi, double Mpi) {

   double k= sqrt( pow(omega,2)/4.0 - pow(Mpi,2));

   return (1.0/6.0)*pow(g_rho_pipi/M_PI,2)*(3.0/4.0 + 2.0*pow(Mpi/omega,2) + (1.0-2.0*pow(Mpi/omega,2))*(k/omega)*log( (omega + 2.0*k)/(2.0*Mpi)));

}



double LL_functions::A_pipi_0(double m_rho, double g_rho_pipi, double Mpi, double kappa) {

  return h(m_rho,g_rho_pipi, Mpi) - (m_rho/2.0)*h_prime(m_rho, g_rho_pipi, Mpi) + pow(g_rho_pipi,2)*pow(Mpi,2)/(6.0*pow(M_PI,2)) + kappa*(pow(m_rho,2)/8.0)*(h_second(m_rho, g_rho_pipi, Mpi) - h_prime(m_rho, g_rho_pipi, Mpi)/m_rho) ;  

}


double LL_functions::Re_A_pipi(double omega, double m_rho, double g_rho_pipi, double Mpi, double kappa) {

  return h(m_rho, g_rho_pipi, Mpi) + (pow(omega,2)- pow(m_rho,2))*(h_prime(m_rho, g_rho_pipi,Mpi)/(2.0*m_rho)) - h(omega, g_rho_pipi,Mpi) + (kappa/2.0)*pow((pow(omega,2)-pow(m_rho,2)),2)*(1.0/(pow(2.0*m_rho,2)))*( h_second(m_rho,g_rho_pipi,Mpi) - (1.0/m_rho)*h_prime(m_rho, g_rho_pipi,Mpi));

  
}

double LL_functions::Im_A_pipi(double omega, double g_rho_pipi, double Mpi) {

  return omega*Gamma_rpp(omega,g_rho_pipi,Mpi);

}

Pfloat LL_functions::A_pipi(double omega, double m_rho, double g_rho_pipi, double Mpi, double kappa) {

  return make_pair(Re_A_pipi(omega, m_rho, g_rho_pipi,Mpi, kappa), Im_A_pipi(omega, g_rho_pipi,Mpi));

}

double LL_functions::F_pi_GS_mod(double omega, double m_rho, double g_rho_pipi, double Mpi, double kappa) {

  double Re = pow(m_rho,2) -pow(omega,2) - Re_A_pipi(omega,m_rho, g_rho_pipi,Mpi, kappa);
  double Im = -1.0*Im_A_pipi(omega, g_rho_pipi, Mpi);
  return fabs((pow(m_rho,2) - A_pipi_0(m_rho, g_rho_pipi,Mpi, kappa)))/( sqrt( pow(Re,2)+ pow(Im,2))); 

}

double LL_functions::cot_d_11(double k, double m_rho, double g_rho_pipi, double Mpi, double kappa) {


  double omega= 2.0*sqrt( pow(Mpi,2) + pow(k,2));

  double Im_a_pipi= Im_A_pipi(omega, g_rho_pipi, Mpi);

  double res = (pow(m_rho,2) - pow(omega,2) - Re_A_pipi(omega, m_rho, g_rho_pipi,Mpi, kappa))/Im_a_pipi;

  return res;

}

double LL_functions::d_11(double k, double m_rho, double g_rho_pipi, double Mpi, double kappa) {


  double d_11=  cot_d_11(k, m_rho, g_rho_pipi, Mpi, kappa);
  
  //double sgn= 1.0;

  //if(d_11 < 0.0) sgn = -1.0;

  double ret_val = atan(1.0/d_11);
  
  //ret_val=  asin( sgn/sqrt(1.0 + pow(d_11,2)));

  return ret_val;

}


double LL_functions::cot_d_11_der(double k, double m_rho, double g_rho_pipi, double Mpi, double kappa) {


  double omega = 2.0*sqrt( pow(Mpi,2) + pow(k,2));

  double d_omega_d_k = 4.0*k/omega;

  double num= ( pow(m_rho,2) -pow(omega,2) - Re_A_pipi(omega, m_rho, g_rho_pipi, Mpi, kappa));


   double den = Im_A_pipi(omega, g_rho_pipi, Mpi);

   double num_der = (-2.0*omega -2.0*omega*(h_prime(m_rho, g_rho_pipi,Mpi)/(2.0*m_rho)) + h_prime(omega,g_rho_pipi,Mpi)    -(2.0*kappa*omega*( pow(omega,2) -pow(m_rho,2)   ))*(1.0/(pow(2.0*m_rho,2)))*( h_second(m_rho,g_rho_pipi,Mpi) - (1.0/m_rho)*h_prime(m_rho, g_rho_pipi,Mpi)))*d_omega_d_k;

  
  double den_der = (Gamma_rpp(omega, g_rho_pipi, Mpi) + omega*Gamma_rpp_der(omega, g_rho_pipi, Mpi))*d_omega_d_k;

 
  double res_2 = cot_d_11(k, m_rho, g_rho_pipi,Mpi, kappa)*( num_der/num - den_der/den);

 

  return res_2;


}


double LL_functions::d_11_der(double k, double m_rho, double g_rho_pipi, double Mpi, double kappa) {



  return -1.0*(1.0/(1.0 + pow(cot_d_11(k,m_rho,g_rho_pipi,Mpi, kappa),2)))*cot_d_11_der(k, m_rho, g_rho_pipi, Mpi, kappa);

}

double LL_functions::d_11_der_num(double k, double m_rho, double g_rho_pipi, double Mpi, double kappa) {

  auto F = [&](double x) { return d_11(x, m_rho, g_rho_pipi,Mpi, kappa);};

  //double x=z;

  double result, err;

  gsl_function_pp<decltype(F)> Fp(F);

  gsl_function *G = static_cast<gsl_function*>(&Fp);

  gsl_deriv_central(G, k, 1e-3, &result, &err);


  if(err/result > 1.0e-5) crash("numerical derivative in d11'(z) has reached a low level of accuracy in d_11_der_num. rel_err: "+to_string_with_precision(err/result,4));

  return result;


}





double LL_functions::Amplitude(double k, double L, double m_rho, double g_rho_pipi, double Mpi, double kappa) {

  
  double omega= 2.0*sqrt( pow(k,2)+pow(Mpi,2));


  return ((2.0*pow(k,5)/(3.0*M_PI*pow(omega,2))))*pow(F_pi_GS_mod(omega, m_rho,g_rho_pipi,Mpi, kappa),2)/(k*d_11_der(k, m_rho, g_rho_pipi,Mpi, kappa) + (k*L/(2.0*M_PI))*phi_der_spline(k*L/(2.0*M_PI)));

}



void LL_functions::Find_pipi_energy_lev(double L, double m_rho, double g_rho_pipi, double Mpi, double kappa,  Vfloat &res) {

  //find energy level n of pi-pi bound state

  Vfloat divergent_levels= Luscher_zeroes;

  double z_crit = (L/(2.0*M_PI))*sqrt( pow(m_rho,2)/4 - pow(Mpi,2));
  double z2_crit= z_crit*z_crit;


  //Mpi=0.135;
  //add z2_crit to divergent_levels
  divergent_levels.push_back(z2_crit);
  sort(divergent_levels.begin(), divergent_levels.end());


  res.clear();

  auto F = [&](double k) -> double {
	     return d_11(k, m_rho, g_rho_pipi, Mpi, kappa) + phi_spline(k*L/(2.0*M_PI));
	   };

  auto F_exact = [&](double k) -> double {
		   return d_11(k, m_rho, g_rho_pipi, Mpi,kappa) + phi(k*L/(2.0*M_PI));
	   };

  auto dF = [&](double k) -> double {

	      return -cot_d_11_der(k, m_rho, g_rho_pipi, Mpi,kappa)/pow(cot_d_11(k, m_rho, g_rho_pipi, Mpi,kappa),2) + (L/(2.0*M_PI))*tan_phi_der(k*L/(2.0*M_PI));
	   };

  auto fdF = [&](double k, double* f, double *df) -> void {
	       *f = F(k);
	       *df = dF(k);
	       return;
	     };


  double Precision = 1e-4;

  
  for(int ilev=1; ilev<=Nres;ilev++) {
  
    int status;
    int iter = 0, max_iter = 100;
    // const gsl_root_fdfsolver_type *T;
    const gsl_root_fsolver_type *T_Br;
    //gsl_root_fdfsolver *s;
    gsl_root_fsolver *s_Br;
    double x0,x;
  
    
    // gsl_function_fdf_pp<decltype(F), decltype(dF), decltype(fdF)> Fp(F, dF, fdF);
    gsl_function_pp<decltype(F)> Fp_Brent(F);
    gsl_function_pp<decltype(F_exact)> Fp_Brent_exact(F_exact);
    
    // gsl_function_fdf *FDF = static_cast<gsl_function_fdf*>(&Fp);
    gsl_function *F_Brent = static_cast<gsl_function*>(&Fp_Brent);
    gsl_function *F_Brent_exact = static_cast<gsl_function*>(&Fp_Brent_exact);
    
 
    
    //   T = gsl_root_fdfsolver_newton;
    // s = gsl_root_fdfsolver_alloc (T);
    T_Br= gsl_root_fsolver_brent;
    s_Br= gsl_root_fsolver_alloc(T_Br);
   

    double XL= (2.0*M_PI/L)*sqrt(divergent_levels[ilev-1]);
    double XH= (2.0*M_PI/L)*sqrt(divergent_levels[ilev]);
    double XL_start= XL;
    double XH_start = XH;
    double DX= (XH -XL)*1e-3;
    XL+=(XH-XL)*1e-6;
    XH-=(XH-XL)*1e-6;

  
    gsl_set_error_handler_off();
    int Brent_status=gsl_root_fsolver_set(s_Br,F_Brent, XL, XH);

    //gsl_root_fdfsolver_set(s, FDF, x);
    
    while(Brent_status != GSL_SUCCESS) {
      XL+=DX;
      XH-=DX;
      Brent_status= gsl_root_fsolver_set(s_Br,F_Brent, XL, XH);
     

      
      if(XL > XH){
	cout.precision(15);
	cout<<"Brent algo is not able to bracket the root"<<endl;
	cout<<"status: "<<Brent_status<<endl;
	cout<<"PRINTING INFO: "<<endl;
	cout<<"Mrho: "<<m_rho<<endl;
	cout<<"Mpi: "<<Mpi<<endl;
	cout<<"g: "<<g_rho_pipi<<endl;
	cout<<"kappa: "<<kappa<<endl;
	cout<<"L: "<<L<<endl;
	cout<<"ilev: "<<ilev<<endl;
	cout<<"Bracketing values: "<<XL<<" "<<XH<<endl;
	cout<<"Starting bracketing values: "<<XL_start<<" "<<XH_start<<endl;
	cout<<"Values of F at bracketing values (spline) : "<<F(XL)<<" "<<F(XH)<<endl;
	cout<<"Values of F at bracketing values (exact) : "<<F_exact(XL)<<" "<<F_exact(XH)<<endl;
	cout<<"Values of F at starting bracketing values (spline) : "<<F(XL_start)<<" "<<F(XH_start)<<endl;
	cout<<"Values of F at starting bracketing values (exact) : "<<F_exact(XL_start)<<" "<<F_exact(XH_start)<<endl;
	printV(divergent_levels, "Printing diverging levels....", 1);
	cout<<"rho resonance divergence is at z^2: "<<z2_crit<<endl;
	crash("Exiting....");
      }
    }

    double x_lo,x_hi;
    
    do
    {
      iter++;
      status = gsl_root_fsolver_iterate(s_Br);
      //x0 = x;
      x = gsl_root_fsolver_root (s_Br);
      //status = gsl_root_test_delta (x, x0, 0, Precision);

      x_lo = gsl_root_fsolver_x_lower(s_Br);
      x_hi = gsl_root_fsolver_x_upper(s_Br);
      status = gsl_root_test_interval (x_lo, x_hi, 0.0, Precision);
      
    }
    while (status == GSL_CONTINUE && iter < max_iter);

    if(status != GSL_SUCCESS) {

      cout<<"gsl_root_solver was unable to find a solution to the Luscher quantization equation for ilev: "<<ilev<<endl;
      cout<<"Printing info: "<<endl;
      cout<<"L: "<<L<<endl;
      cout<<"Mrho: "<<m_rho<<endl;
      cout<<"g: "<<g_rho_pipi<<endl;
      cout<<"Mpi: "<<Mpi<<endl;
      cout<<"kappa: "<<kappa<<endl;
      cout<<"Last root position: "<<x<<endl;
      cout<<"divergent level position 1: "<<(2.0*M_PI/L)*sqrt(divergent_levels[ilev-1]+2e-2)<<endl;
      cout<<"divergent level position 2: "<<(2.0*M_PI/L)*sqrt(divergent_levels[ilev]+2e-2)<<endl;
      crash("Aborting...");
    }


    //if(switch_to_exact_function) {
    //cout.precision(12);
    // cout<<"ilev: "<<ilev<<" switched to exact function evaluation"<<endl;
    // cout<<"root(approximate): "<<z_failed<<" root exact: "<<x<<endl;
    //}

    //if(!switch_to_exact_function) { root_found= F_exact(x_lo)*F_exact(x_hi) <= 0; z_failed=x;}
    //if(!root_found) switch_to_exact_function=true;
    
    // gsl_root_fdfsolver_free (s);
    gsl_root_fsolver_free(s_Br);

    res.push_back(x);

  }

  

  return;
  
}


double LL_functions::V_pipi(double t, double L, double m_rho, double g_rho_pipi, double Mpi, double kappa,  Vfloat &Knpp) {

  double ret_val=0;

 
 
  for(int i_lev=0;i_lev<Nres;i_lev++) {
    double k_n= Knpp[i_lev];  
    double omega_n = 2.0*sqrt( pow(Mpi,2) + pow(k_n,2));
    double Ampl= Amplitude(k_n, L, m_rho, g_rho_pipi, Mpi,kappa);
    ret_val += Ampl*exp(-omega_n*t);
   

  }
  

  return ret_val;

}


double LL_functions::V_pipi_infL(double t, double m_rho_infL, double g_rho_pipi_infL, double Mpi_infL, double kappa_infL) {

  double tol_infL= 1e-14;

  auto Integrand = [&](double omega) -> double {
		     return (1.0/(48.0*pow(M_PI,2)))*pow(omega,2)*pow(1.0- pow(2.0*Mpi_infL/omega,2), 3.0/2.0)*exp(-omega*t)*pow(F_pi_GS_mod(omega, m_rho_infL, g_rho_pipi_infL,Mpi_infL,kappa_infL),2);};


  double val,err;
  double tolerance = 1e-9;
  gsl_function_pp<decltype(Integrand)> F_corr(Integrand);
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (10000);
  gsl_function *G = static_cast<gsl_function*>(&F_corr);
  gsl_integration_qagiu(G, 2*Mpi_infL, 0.0, tolerance, 10000, w, &val, &err);
  gsl_integration_workspace_free (w);
  
  if( err/fabs(val) > 5*tolerance) crash("In free_vector_corr_cont gls integration not able to achieve target precision");

  return val;

  

  
}



void LL_functions::MLLGS_fit_to_corr(const distr_t_list &Corr,const distr_t &Mpi,const distr_t &a_distr, double L, distr_t &Edual, distr_t &Rdual, distr_t &Mrho, distr_t &grpp, int tmin, int tmax, string Tag) {

  int NUMM_THREADS= omp_get_max_threads();
  omp_set_num_threads(1);
   
  class ipar_MLLGS {

  public:
    ipar_MLLGS() : V_light(0.0), V_light_err(0.0) {}
    double Mp;
    double t,  L;
    double V_light, V_light_err;
  };

  class fit_par_MLLGS {

  public:
    fit_par_MLLGS() {}
    fit_par_MLLGS(const Vfloat &par) {
      if((signed)par.size() != 5) crash("In class fit_par_MLLGS in fitting analytic representation of V(t)_light, class constructor Vfloat par has size != 5");
      Rd=par[0];
      Ed=par[1];
      Mrho=par[2];
      gpi=par[3];
      kappa= par[4];
    }

    double Rd,Ed, Mrho, gpi,kappa;
  };
  
  int Njacks= a_distr.size();
  
  bootstrap_fit<fit_par_MLLGS,ipar_MLLGS> bf_MLLGS(Njacks);
  //bf_MLLGS.set_warmup_lev(2); //sets warmup
  int Tfit_points= tmax +1 - tmin;
  
  

  
  bf_MLLGS.Set_number_of_measurements(Tfit_points);
  bf_MLLGS.Set_verbosity(1);

  cout<<"FITTING WITH PIPI+DUAL REPRESENTATION: "<<Tag<<endl;
  

  //initial guesses are based on the results on the old ETMC confs
  bf_MLLGS.Add_par("Rd", 1.1, 0.1);
  bf_MLLGS.Set_limits("Rd", 0.7, 1.8);
  bf_MLLGS.Add_par("Ed", 2.5, 0.2);
  bf_MLLGS.Set_limits("Ed", 0.7, 5.0);
  bf_MLLGS.Add_par("Mrho",5.5, 0.1);
  bf_MLLGS.Set_limits("Mrho", 4.5 , 6.7);
  bf_MLLGS.Add_par("gpi", 1.0, 0.01);
  bf_MLLGS.Add_par("kappa", -3.0, 0.1);
  bf_MLLGS.Fix_par("kappa", 0.0);
  bf_MLLGS.Set_limits("gpi",0.6, 1.4);
  
  map<pair<pair<double,double>, pair<double,double>>,Vfloat> Energy_lev_list;
  

  bf_MLLGS.ansatz =  [&Energy_lev_list, this](const fit_par_MLLGS &p, const ipar_MLLGS &ip) -> double {

		     double Pi_M = ip.Mp;

		     double GPI = p.gpi*5.95;
		     Vfloat Knpp;
		     pair<double,double> Mass_par= make_pair(p.Mrho*Pi_M, Pi_M);
		     pair<double,double> Couplings = make_pair(GPI, p.kappa);
		     pair< pair<double,double>,pair<double, double>> input_pars = make_pair( Mass_par, Couplings);
		     map<pair<pair<double,double>, pair<double, double>>,Vfloat>::iterator it;
		     it= Energy_lev_list.find(input_pars);
		     if(it != Energy_lev_list.end()) Knpp= it->second;
		     else {
		       Find_pipi_energy_lev(ip.L,p.Mrho*Pi_M, GPI, Pi_M, p.kappa, Knpp);
		       //add to Energy_lev_list
		       Energy_lev_list.insert( make_pair(input_pars, Knpp));
		     }
		
		     return 2.0*V_pipi(ip.t, ip.L, p.Mrho*Pi_M, GPI, Pi_M, p.kappa, Knpp) + (9.0/5.0)*Vdual(ip.t, p.Mrho*Pi_M, p.Ed*Pi_M, p.Rd);
		   };
      bf_MLLGS.measurement = [](const fit_par_MLLGS& p,const ipar_MLLGS& ip) -> double {
			 return ip.V_light;
		       };
      bf_MLLGS.error =  [](const fit_par_MLLGS& p,const ipar_MLLGS &ip) -> double {
		    return ip.V_light_err;
		  };

      //fill the data
      vector<vector<ipar_MLLGS>> data(Njacks);
      //allocate space for output result
      boot_fit_data<fit_par_MLLGS> Bt_fit;

  

      for(auto &data_iboot: data) data_iboot.resize(Tfit_points);
  
      for(int ijack=0;ijack<Njacks;ijack++) {
	for(int t=tmin;t<Tfit_points+tmin;t++) {
	  int tt=t-tmin;
	  data[ijack][tt].V_light = Corr.distr_list[t].distr[ijack];
	  data[ijack][tt].V_light_err = Corr.err(t);
	  data[ijack][tt].L = L;
	  data[ijack][tt].Mp = Mpi.distr[ijack];
	  data[ijack][tt].t = t;
	} 
      }
    
      //append
      bf_MLLGS.Append_to_input_par(data);
      //fit
      Bt_fit= bf_MLLGS.Perform_bootstrap_fit();

      for(int ijack=0;ijack<Njacks;ijack++) {
        Mrho.distr.push_back(Bt_fit.par[ijack].Mrho*Mpi.distr[ijack]);
	grpp.distr.push_back( Bt_fit.par[ijack].gpi*5.95);
	Edual.distr.push_back( Bt_fit.par[ijack].Ed*Mpi.distr[ijack]);
	Rdual.distr.push_back( Bt_fit.par[ijack].Rd);
      }


      Energy_lev_list.clear();

      omp_set_num_threads(NUMM_THREADS);

      return; 
}
