#include "../include/Num_integrate_l4_decay_rate.h"

using namespace std;


const double r_mu =  pow(0.10565837/0.493677,2); //0.046 //squared ratio m_mu / m_K
const double r_el = pow( 0.000510998950/0.493677,2); //1.1e-6 //squared ratio m_el / mK
const double Gf = 1.1663787*1e-5; //Fermi constant in GeV^-2
const double Vus=  0.2245; //0.220; //  0.2245; //Vus CKM matrix element
const double Gamma = 5.317*1e-17; // Total decay width of the K^+ in GeV
const double alpha = 1/137.036; //alpha constant
const double eps_mumu= 1e-15;
const double eps_ee = 1e-15;
const double MkPh=  0.493677; //GeV
const double fkPh =0.155; // 0.155;// 0.1136*sqrt(2); //0.155 ;//GeV
string MODE= "TOTAL";
double cut = pow(0.140/MkPh,2); //pow(0.140/MkPh,2); //pow(0.050/MkPh,2); //  //4.0*r_el; // //0.001; // // pow(0.020/MkPh,2); 0.001; //
double cut_mumu= 4.0*r_mu;
double error_goal_MC_ee = 100*1.0e-2*(Gamma/pow(Gf*Vus*alpha,2))*8.0*1.0e-8;
double error_goal_MC_mumu= 100*1.0e-2*(Gamma/pow(Gf*Vus*alpha,2))*1.0e-8;
double F_same_lepton= pow(Gf*Vus*alpha,2)*(1.0/Gamma)*(1.0/(pow(M_PI,4)*pow(2.0,12)));
double error_goal_MC_eee = 0.5e-6*(3.0/F_same_lepton)*1e-8;
double error_goal_MC_mumumu= 0.5e-6*(1.2/F_same_lepton)*1e-8;
double error_goal_MC_ee_extended= 0.1e-8*(8.5/F_same_lepton)*1e-8; 
double error_goal_MC_mumu_extended = 0.5e-8*(0.8/F_same_lepton)*1e-8;
bool USE_CHPT_FORM_FACTORS=false;
const int SEED_VAL= 366432;


void display_results(char *title, double result, double error)
{
  printf ("%s ==================\n", title);
  printf ("result = % .16f\n", result);
  printf ("relative error  = % .16f\n ", 100*error/result);
}


double MonteCarlo_integration_extended_phase_space(const function<double(double,double)> &H1, const function<double(double,double)> &H2, const function<double(double,double)> &FA, const function<double(double,double)> &FV,double mk, double fk, string channel, double xk_inf,  bool same_lepton, int MySeed, string mode) {

  VPfloat bounds;
  double bound_l[5];
  double bound_u[5];
  double rll,rl, error_goal;
  if(channel=="e+e-") {
    if(same_lepton) {
    rll= r_el;
    rl=  r_el;
    }
    else {
      rl=r_mu;
      rll=r_el;
    }
    error_goal= error_goal_MC_eee;
    bounds = VPfloat{ { xk_inf,1-sqrt(rl)}, {sqrt(rl), 1.0}, {-1.0, 1.0}, {-1.0,1.0}, {0.0, 2.0*M_PI}};

    bound_l[0] = xk_inf;
    bound_l[1]= sqrt(rl);
    bound_l[2] = -1.0;
    bound_l[3] = -1.0;
    bound_l[4] = 0.0;
    bound_u[0] = 1.0-sqrt(rl);
    bound_u[1] = 1.0;
    bound_u[2] = 1.0;
    bound_u[3] = 1.0;
    bound_u[4] = 2.0*M_PI;
    

    //bounds = VPfloat{ { cut, pow(1.0-sqrt(rl),2)}, {rl, 1.0}, {-1.0,1.0}, {-1.0,1.0}, {0.0, 2.0*M_PI}};
  }
  else if(channel=="mu+mu-") {
    if(same_lepton) {
    rll=r_mu;
    rl= r_mu;
    }
    else {
      rl= r_el;
      rll=r_mu;
    }
    error_goal= error_goal_MC_mumumu;
    bound_l[0] = xk_inf; //2.0*sqrt(rll);
    bound_l[1]= sqrt(rl);
    bound_l[2] = -1.0;
    bound_l[3] = -1.0;
    bound_l[4] = 0.0;
    bound_u[0] = 1.0-sqrt(rl);
    bound_u[1] = 1.0;
    bound_u[2] = 1.0;
    bound_u[3] = 1.0;
    bound_u[4] = 2.0*M_PI;
    bounds = VPfloat{ { xk_inf,1-sqrt(rl)}, {sqrt(rl), 1.0}, {-1.0, 1.0}, {-1.0,1.0}, {0.0, 2.0*M_PI} };

    
    //bounds = VPfloat{ { 4.0*rll, pow(1.0-sqrt(rl),2)}, {rl, 1.0}, {-1.0,1.0}, {-1.0,1.0}, {0.0, 2.0*M_PI}};
  }
  else crash("In Montecarlo_integration_extended_phase_space channel: "+channel+" not yet implemented");

  //define lambda function to define energies and momenta in terms of the 5 integration variables, xk, xq, y12, y34 and phi. All dimensionful quantities are normalized over the Kaon Mass

  auto lambda = [&](double xk, double xq) -> double { return sqrt( pow((1.0 -pow(xk,2)-pow(xq,2)),2) -4.0*pow(xk*xq,2));};

  auto l_12 = [&](double xk) -> double { return sqrt( 1-4.0*rll/pow(xk,2));};

  auto l_34 = [&](double xq) -> double { return 1-(rl/pow(xq,2));};

  auto delta= [&](double xk, double xq) -> double { return pow(xk,2) -pow(xq,2);};

  auto delta34 = [&](double xq) -> double {return - rl/pow(xq,2);};

  auto y= [&](double xk, double xq, double y12, double y34, double phi) -> double { return ((1.0 - delta(xk,xq))*(1.0 -delta34(xq)) - lambda(xk,xq)*y34)/2.0 -rl;};
  
  auto A= [&](double xk, double xq, double y12, double y34, double phi) -> double {
	    
	    double E1 = ((1.0 + delta(xk,xq)) + lambda(xk,xq)*y12)/4.0;
	    double E4 =  ((1.0 - delta(xk,xq))*(1.0 - delta34(xq)) - lambda(xk,xq)*y34)/4.0;
	   
	    double p1x = -(xk/2.0)*sqrt( pow(l_12(xk),2) - pow(y12,2));
	    double p1y = (lambda(xk,xq) + (1.0 + delta(xk,xq))*y12)/4.0;

	    double p4x = (xq/2.0)*sqrt( pow(l_34(xq),2) - pow(y34,2))*cos(phi);
	    double p4y = -(lambda(xk,xq)*(1.0 - delta34(xq)) - (1.0 - delta(xk,xq))*y34)/4.0;

	    return 2.0*( E1*E4 - p1x*p4x - p1y*p4y);
	  };
  
  auto B= [&](double xk, double xq, double y12, double y34, double phi) -> double {


	    double E2 =  ((1.0 + delta(xk,xq)) - lambda(xk,xq)*y12)/4.0; 
	    double E3 = ((1.0 - delta(xk,xq))*(1.0 + delta34(xq)) + lambda(xk,xq)*y34)/4.0;
	    
	    double p2x = (xk/2.0)*sqrt( pow(l_12(xk),2) - pow(y12,2));
	    double p2y = (lambda(xk,xq) -(1.0 + delta(xk,xq))*y12)/4.0;
	    
	    double p3x = -(xq/2.0)*sqrt( pow(l_34(xq),2) - pow(y34,2))*cos(phi);
	    double p3y = -(lambda(xk,xq)*(1.0 + delta34(xq)) + (1.0 - delta(xk,xq))*y34)/4.0;
	   

	    return 2.0*(E2*E3 - p2x*p3x - p2y*p3y);
	      
	  };

  
  auto C= [&](double xk, double xq, double y12, double y34, double phi) -> double { return (lambda(xk,xq)*xk*xq/8.0)*sin(phi)*sqrt(  (pow(l_12(xk),2) -pow(y12,2))*(pow(l_34(xq),2) -pow(y34,2)));};


  //perform Monte Carlo integration


  auto square_amplitude = [&](vector<double> const &par) -> double {
			      
			     if((signed)par.size() != 5) crash("squared_amplitude in Num_Integrate_Decay_Rate_extendend_phase_space only accepts vector<double> of size 5. "+to_string( (signed)par.size())+" provided.");

			    
			     double xk= par[0];
			     double xq= par[1];
			     if(xq > 1.0-xk) return 0.0;
			     double y12= par[2]*l_12(xk);
			     double y34= par[3]*l_34(xq);
			     double phi=par[4];

			
			     
			     double a= A(xk,xq, y12, y34, phi);
			     double b= B(xk,xq, y12, y34, phi);
			     double Y= y(xk,xq, y12, y34, phi);

			     double xk_prime= sqrt(2*rll + a);
			     double xq_prime= sqrt(rl + b);

			     if(xk_prime < 2.0*sqrt(rll) && same_lepton) crash("invalid xk_prime generated");
			     if(xk < 2.0*sqrt(rll)) crash("invalid xk generated");

			     if(xk_prime < xk_inf && same_lepton) return 0.0;
			  			  
			 
			     double jacobian= 4.0*xk*xq;

			     double square_ampl;

			     if(same_lepton) {

			       /*
			      long double amplitude_v1 = Compute_square_amplitude_extended(H1(xk,xq), H2(xk,xq), FA(xk,xq), FV(xk,xq), H1(xk_prime,xq_prime), H2(xk_prime, xq_prime), FA(xk_prime, xq_prime), FV(xk_prime,xq_prime), mk, fk,rl,xk, xq, a,b,Y, MODE );
			      long double amplitude_v2 = Compute_square_amplitude_extended_v2(H1(xk,xq), H2(xk,xq), FA(xk,xq), FV(xk,xq), H1(xk_prime,xq_prime), H2(xk_prime, xq_prime), FA(xk_prime, xq_prime), FV(xk_prime,xq_prime), mk, fk,rl,xk, xq, a,b,Y, MODE );
			      
			      if( amplitude_v1/amplitude_v2 < 1 -1e-4 || amplitude_v1/amplitude_v2 > 1+ 1e-4) {
				cout.precision(10);
				cout<<"v1: "<<amplitude_v1<<endl;
				cout<<"v2: "<<amplitude_v2<<endl;
				crash("v1 and v2 do not match");
			      }
			       */
			      square_ampl= Compute_square_amplitude_extended_v2(H1(xk,xq), H2(xk,xq), FA(xk,xq), FV(xk,xq), H1(xk_prime,xq_prime), H2(xk_prime, xq_prime), FA(xk_prime, xq_prime), FV(xk_prime,xq_prime), mk, fk,rl,xk, xq, a,b,Y, mode )*(jacobian*l_12(xk)*l_34(xq)*lambda(xk,xq)/mk)*pow(MkPh/mk,5);

			     }

			     else {
			       
			       long double amplitude_v1 = Compute_square_amplitude_different_lepton(H1(xk,xq), H2(xk,xq), FA(xk,xq), FV(xk,xq), mk, fk,rl,rll,xk, xq, a,b,Y, MODE );
			       long double amplitude_v2 = Compute_square_amplitude_different_lepton_v2(H1(xk,xq), H2(xk,xq), FA(xk,xq), FV(xk,xq), mk, fk,rl,rll,xk, xq, a,b,Y, MODE );
			       if( amplitude_v1/amplitude_v2 < 1 -1e-6 || amplitude_v1/amplitude_v2 > 1+ 1e-6) {
				cout.precision(10);
				cout<<"v1: "<<amplitude_v1<<endl;
				cout<<"v2: "<<amplitude_v2<<endl;
				crash("v1 and v2 do not match");
			       }
			       
			       square_ampl = 2.0*(jacobian*l_12(xk)*l_34(xq)*lambda(xk,xq)/mk)*pow(MkPh/mk,5)*Compute_square_amplitude_different_lepton(H1(xk,xq), H2(xk,xq), FA(xk,xq), FV(xk,xq), mk, fk,rl,rll,xk, xq, a,b,Y, mode );

			     }
			     
			     return square_ampl;
			   };
   


    //VEGAS GLS INTEGRATION

    double res_vegas, err_vegas;
    size_t calls = 500000; //old was 1000000; 

    const gsl_rng_type *T;
     gsl_rng *r;

    gsl_rng_env_setup();

    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

   
    //to enter a generic seed mySeed
    gsl_rng_set(r, MySeed);

    gsl_monte_function_pp<decltype(square_amplitude)> Fp(square_amplitude);

    gsl_monte_function *G = static_cast<gsl_monte_function*>(&Fp);


    // G.f = #c function



    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(5);

    int warmup_calls= 50000; //old was 100000;

    gsl_monte_vegas_integrate (G, bound_l, bound_u, 5, warmup_calls, r, s,
                               &res_vegas, &err_vegas);

    //display_results ("vegas warm-up", res_vegas, err_vegas);


    printf ("converging...\n");

    do
      {
        gsl_monte_vegas_integrate (G,bound_l, bound_u, 5, calls/5, r, s,
                                   &res_vegas, &err_vegas);
        //printf ("result = % .6f sigma = % .6f "
	//      "chisq/dof = %.1f\n", res_vegas, err_vegas, gsl_monte_vegas_chisq (s));
      }
    while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);

    display_results("vegas final", res_vegas, err_vegas);

    gsl_monte_vegas_free(s);

    gsl_rng_free(r);


    return res_vegas*F_same_lepton;

    
}


Decay_Rate_Integration_Result Num_Integrate_Decay_Rate(vector<function<double(double, double)>> &H1, vector<function<double(double, double)>> &H2  , vector<function<double(double, double)>> &FA,  vector<function<double(double, double)>> &FV, distr_t &m_distr, distr_t & fp_distr, bool UseJack, string Meson, bool Print_rate) {


  

  
  if( (H1.size() != H2.size()) || (H1.size() != FA.size()) || (H1.size() != FV.size()) || (H1.size() != (unsigned)m_distr.size()) || (H1.size() != (unsigned)fp_distr.size())) crash("Sizes of form factors in Decay_Rate_Integration_Result are different.");

  int Njacks= H1.size();
  Decay_Rate_Integration_Result Res(UseJack);

  double eps;

  //init RandomMersenne to be used in Extended phase space integration
  RandomMersenne   RM(SEED_VAL, ipow(2,13));

  

  //relative accuracy in integration
  Res.eps_rel_mumu= eps_mumu;
  Res.eps_rel_ee = eps_ee;


  //define distr_t_list to store the differential rate (unequal leptons decay)
  int Npoints_d=100;
  double off_d_e = 2.0*sqrt(r_el);
  double step_p_d_e= 1e-2;
  double step_p_d_e_2 = 3e-6;
  double off_d_mu = 2.0*sqrt(r_mu);
  double step_p_d_mu = (1.0-sqrt(r_mu) - 2.0*sqrt(r_mu))/Npoints_d;
  double step_p_d_mu_2 = step_p_d_mu/1000;
  Vfloat p_list_d_e;
  Vfloat p_list_d_mu;
  for(int ip=0;ip<Npoints_d;ip++) p_list_d_e.push_back( off_d_e + ip*step_p_d_e_2);
  for(int ip=1;ip<Npoints_d;ip++) p_list_d_e.push_back( off_d_e + ip*step_p_d_e);
  for(int ip=0;ip<Npoints_d;ip++) p_list_d_mu.push_back( off_d_mu + ip*step_p_d_mu_2);
  for(int ip=1;ip<Npoints_d;ip++) p_list_d_mu.push_back( off_d_mu + ip*step_p_d_mu);
  distr_t_list Plot_diff_rate_ee_tot(UseJack, (signed)p_list_d_e.size() );
  distr_t_list Plot_diff_rate_ee_SD(UseJack, (signed)p_list_d_e.size() );
  distr_t_list Plot_diff_rate_ee_pt(UseJack, (signed)p_list_d_e.size() );
  distr_t_list Plot_diff_rate_ee_Q(UseJack, (signed)p_list_d_e.size());
  distr_t_list Plot_diff_rate_ee_int(UseJack, (signed)p_list_d_e.size());
  distr_t_list Plot_diff_rate_mumu_tot(UseJack, (signed)p_list_d_mu.size());
  distr_t_list Plot_diff_rate_mumu_SD(UseJack, (signed)p_list_d_mu.size());
  distr_t_list Plot_diff_rate_mumu_pt(UseJack, (signed)p_list_d_mu.size());
  distr_t_list Plot_diff_rate_mumu_Q(UseJack, (signed)p_list_d_mu.size());
  distr_t_list Plot_diff_rate_mumu_int(UseJack, (signed)p_list_d_mu.size());

  //define distr_t_list to store the total rate as a function of the lower cut
  int Npoints_c= 20;
  double off_c_e = 3.5e-2;
  double step_p_c_e = 3.5e-2;
  double off_c_mu = 2.0*sqrt(r_mu);
  double step_p_c_mu = (0.67 -2.0*sqrt(r_mu))/Npoints_c;
  Vfloat p_list_c_e;
  Vfloat p_list_c_mu;
  for(int ip=0; ip<Npoints_c;ip++) p_list_c_e.push_back( off_c_e + ip*step_p_c_e);
  for(int ip=0; ip<Npoints_c;ip++) p_list_c_mu.push_back( off_c_mu + ip*step_p_c_mu);
  distr_t_list Plot_cut_rate_eee_tot(UseJack, Npoints_c );
  distr_t_list Plot_cut_rate_eee_SD(UseJack, Npoints_c );
  distr_t_list Plot_cut_rate_eee_pt(UseJack, Npoints_c );
  distr_t_list Plot_cut_rate_mumumu_tot(UseJack, Npoints_c);
  distr_t_list Plot_cut_rate_mumumu_SD(UseJack, Npoints_c);
  distr_t_list Plot_cut_rate_mumumu_pt(UseJack, Npoints_c);

  //loop over jackknives
  for(int ijack=0; ijack< Njacks; ijack++) {
    //Compute decay rates for a given jackknife

    double m= m_distr.distr[ijack];
    double fp= fp_distr.distr[ijack];

    
    if(USE_CHPT_FORM_FACTORS) {
      m=MkPh;
      fp=fkPh;
      if(MODE !="PT") { 
      Compute_ChPT_form_factors(H1[ijack],H2[ijack], FA[ijack], FV[ijack]);
      }
    }


    
    
    //m=MkPh;
    //fp=fkPh;

    
 

    

    //define double differential decay rate
    auto Rt_diff = [&](double xk, double xq, double l, double ll) -> double {

		     double Int= ptrate(xk,xq, l, ll, m, fp);
		     double interference= H1[ijack](xk,xq)*kern1(xk, xq, l, ll, m,fp) + H2[ijack](xk,xq)*kern2(xk, xq, l, ll, m,fp) + FA[ijack](xk,xq)*kernA(xk,xq,l,ll,m,fp) +FV[ijack](xk,xq)*kernV(xk,xq,l,ll,m,fp);
		     double quadratic= pow(H1[ijack](xk,xq),2)*kern11(xk,xq,l,ll,m) + pow(H2[ijack](xk,xq),2)*kern22(xk,xq,l,ll,m) + pow(FA[ijack](xk,xq),2)*kernAA(xk,xq,l,ll,m) + pow(FV[ijack](xk,xq),2)*kernVV(xk,xq,l,ll,m) + H1[ijack](xk,xq)*H2[ijack](xk,xq)*kern12(xk,xq,l,ll,m) + H1[ijack](xk,xq)*FA[ijack](xk,xq)*kernA1(xk,xq,l,ll,m);
		     double jacobian= 4.0*xk*xq;
		     
		     if(MODE=="PT") return Int*jacobian*pow(MkPh/m,5);
		     else if(MODE=="INTERFERENCE") return interference*jacobian*pow(MkPh/m,5);
		     else if(MODE=="QUADRATIC") return quadratic*jacobian*pow(MkPh/m,5);
		     else if(MODE=="SD") return (quadratic+interference)*jacobian*pow(MkPh/m,5);
		     else if(MODE=="TOTAL") return (Int+interference+quadratic)*jacobian*pow(MkPh/m,5);
		     else crash(" In Rt_diff MODE: "+MODE+" not yet implemented");
    };

    cout<<"#####Quad integration#####"<<endl;
    //e+e-

    double rl=r_mu;
    double rll=r_el;

 
    eps=eps_ee;

    auto Fxk= [&](double xk) -> double {
      double inf = sqrt(rl);
      double sup = 1-xk;
      auto g = [&](double xq) {
	return Rt_diff(xk, xq, rl, rll);
        };
      return boost::math::quadrature::gauss_kronrod<double, 15>::integrate(g, inf, sup, 5, eps);
    };

    double error_ee;
    double res_ee= (pow(alpha*Gf*Vus,2)/Gamma)*boost::math::quadrature::gauss_kronrod<double,15>::integrate(Fxk, sqrt(cut), 1-sqrt(rl), 5, eps, &error_ee);
    Res.Jack_Distr_Int_Quad_ee.distr.push_back(res_ee );
    cout<<"Relative error achieved in Quad e+e- for jack: "<<ijack<<" is: "<<(pow(alpha*Gf*Vus,2)/Gamma)*error_ee<<endl;
    cout<<"Branching ratio Quad e+e-: "<<res_ee<<endl;
    cout<<"################"<<endl;

    //compute diff rate e+e- for plotting

    
    if(Print_rate) {
      string OLD_MODE= MODE;
      for(int ip=0; ip<(signed)p_list_d_e.size(); ip++) {
	double xk_p= p_list_d_e[ip];
	MODE= "TOTAL";
	Plot_diff_rate_ee_tot.distr_list[ip].distr.push_back( Fxk(xk_p)*pow(alpha*Gf*Vus,2)/Gamma);
	MODE= "SD";
	Plot_diff_rate_ee_SD.distr_list[ip].distr.push_back( Fxk(xk_p)*pow(alpha*Gf*Vus,2)/Gamma);
	MODE= "QUADRATIC";
	Plot_diff_rate_ee_Q.distr_list[ip].distr.push_back( Fxk(xk_p)*pow(alpha*Gf*Vus,2)/Gamma);
	MODE= "PT";
	Plot_diff_rate_ee_pt.distr_list[ip].distr.push_back( Fxk(xk_p)*pow(alpha*Gf*Vus,2)/Gamma);
	MODE="INTERFERENCE";
	Plot_diff_rate_ee_int.distr_list[ip].distr.push_back( Fxk(xk_p)*pow(alpha*Gf*Vus,2)/Gamma);
      }
      //restore old mode
      MODE=OLD_MODE;
    }

    

      
    
   
    //mu+mu-

    rl=r_el;
    rll=r_mu;
    eps=eps_mumu;
    

    double error_mm;
    double res_mm = (pow(alpha*Gf*Vus,2)/Gamma)*boost::math::quadrature::gauss_kronrod<double,15>::integrate(Fxk, 2*sqrt(rll), 1-sqrt(rl), 5, eps, &error_mm);
    Res.Jack_Distr_Int_Quad_mumu.distr.push_back(res_mm);
    cout<<"Relative error achieved in Quad mu+mu- for jack: "<<ijack<<" is: "<<(pow(alpha*Gf*Vus,2)/Gamma)*error_mm<<endl;
    cout<<"Branching ratio Quad  mu+mu-: "<<res_mm<<endl;
    cout<<"################"<<endl;
    //perform MonteCarlo Integration
    cout<<"#####Monte Carlo integration#####"<<endl;

    //compute diff rate mu+mu- for plotting

    if(Print_rate) {
    string OLD_MODE= MODE;
    for(int ip=0; ip<(signed)p_list_d_mu.size(); ip++) {
	double xk_p= p_list_d_mu[ip];
	MODE= "TOTAL";
	Plot_diff_rate_mumu_tot.distr_list[ip].distr.push_back( Fxk(xk_p)*pow(alpha*Gf*Vus,2)/Gamma);
	MODE= "SD";
	Plot_diff_rate_mumu_SD.distr_list[ip].distr.push_back( Fxk(xk_p)*pow(alpha*Gf*Vus,2)/Gamma);
	MODE= "QUADRATIC";
	Plot_diff_rate_mumu_Q.distr_list[ip].distr.push_back( Fxk(xk_p)*pow(alpha*Gf*Vus,2)/Gamma);
	MODE= "PT";
	Plot_diff_rate_mumu_pt.distr_list[ip].distr.push_back( Fxk(xk_p)*pow(alpha*Gf*Vus,2)/Gamma);
	MODE="INTERFERENCE";
	Plot_diff_rate_mumu_int.distr_list[ip].distr.push_back( Fxk(xk_p)*pow(alpha*Gf*Vus,2)/Gamma);
      }
      //restore old mode
      MODE=OLD_MODE;

    }
  




    
    //define lambda function for Monte Carlo Integration


    auto Rt_diff_MC = [&](Vfloat const &par) -> double {
      if(par.size() != 2) crash("Rt_diff_MC in Num_Integrate_Decay_Rate only accepts vector<double> of size 2. "+to_string( (signed)par.size())+" provided.");
      double xk= par[0];
      double xq= par[1];
      if(xq +xk > 1) return 0;
      
      return Rt_diff(xk,xq,rl,rll);
    };



    //e+e-

    rl=r_mu;
    rll=r_el;
    VPfloat bounds_ee{ {sqrt(cut),1.0-sqrt(rl)}, { sqrt(rl),1.0-sqrt(cut)}};

    /*

    
    boost::math::quadrature::naive_monte_carlo<double, decltype(Rt_diff_MC)> mc_ee(Rt_diff_MC, bounds_ee, error_goal_MC_ee,  true, thread::hardware_concurrency() - 1);

    
    future<double> task_ee = mc_ee.integrate();
    while (task_ee.wait_for(chrono::seconds(1)) != future_status::ready)
      {
	// do nothing
      }

    double res_ee_MC = (pow(alpha*Gf*Vus,2)/Gamma)*task_ee.get();
    Res.Jack_Distr_Int_MonteCarlo_ee.distr.push_back(res_ee_MC);
    Res.Stat_err_MonteCarlo_ee.push_back( (pow(alpha*Gf*Vus,2)/Gamma)*error_goal_MC_ee/res_ee_MC);
    cout<<"Relative error achieved in Monte Carlo e+e- for jack: "<<ijack<<" is: "<<(pow(alpha*Gf*Vus,2)/Gamma)*error_goal_MC_ee/res_ee_MC<<endl;
    cout<<"Branching ratio Monte Carlo e+e-: "<<res_ee_MC<<endl;
    cout<<"################"<<endl;


    //mu+mu-

    rl=r_el;
    rll=r_mu;

    VPfloat bounds_mumu{ {2.0*sqrt(rll), 1-sqrt(rl)},  {sqrt(rl) , 1.0- 2.0*sqrt(rll) }};
    
    boost::math::quadrature::naive_monte_carlo<double, decltype(Rt_diff_MC)> mc_mumu(Rt_diff_MC, bounds_mumu, error_goal_MC_mumu,  true, thread::hardware_concurrency() - 1);

    
    future<double> task_mumu = mc_mumu.integrate();
    while (task_mumu.wait_for(chrono::seconds(1)) != future_status::ready)
      {
	// do nothing
      }

    double res_mumu_MC = (pow(alpha*Gf*Vus,2)/Gamma)*task_mumu.get();
    Res.Jack_Distr_Int_MonteCarlo_mumu.distr.push_back(res_mumu_MC);
    Res.Stat_err_MonteCarlo_mumu.push_back(error_goal_MC_mumu/res_mumu_MC);
    cout<<"Relative error achieved in Monte Carlo mu+mu- for jack: "<<ijack<<" is: "<<(pow(alpha*Gf*Vus,2)/Gamma)*error_goal_MC_mumu/res_mumu_MC<<endl;
    cout<<"Branching ratio Monte Carlo mu+mu-: "<<res_mumu_MC<<endl;
    cout<<"################"<<endl;

    */


    

    //e+e-e+

    
    

    double res_eee_MC = MonteCarlo_integration_extended_phase_space(H1[ijack], H2[ijack], FA[ijack], FV[ijack], m, fp, "e+e-", sqrt(cut),1,  RM(), MODE);

    Res.Jack_Distr_Int_MonteCarlo_eee.distr.push_back(res_eee_MC);
    Res.Stat_err_MonteCarlo_eee.push_back(error_goal_MC_eee/res_eee_MC);
    cout<<"Relative error achieved in Monte Carlo e+e-e+ for jack: "<<ijack<<" is: "<<F_same_lepton*error_goal_MC_eee/res_eee_MC<<endl;
    cout<<"Branching ratio Monte Carlo e+e-e+: "<<res_eee_MC<<endl;
    cout<<"################"<<endl;


    //store total rate eee as a function of the inf cut for eee
    if(Print_rate) {
      for(int ip=0; ip<Npoints_c;ip++) {
	double xk_cut = off_c_e+ ip*step_p_c_e;
	string MODE_NEW="TOTAL";
	cout<<"computing  e+e-e+ for xk_cut: "<<xk_cut<<" MODE: "<<MODE_NEW<<endl;
	double res_eee_tot_MC = MonteCarlo_integration_extended_phase_space(H1[ijack], H2[ijack], FA[ijack], FV[ijack], m, fp, "e+e-", xk_cut, 1, RM(), MODE_NEW);
	Plot_cut_rate_eee_tot.distr_list[ip].distr.push_back(res_eee_tot_MC);
	MODE_NEW="PT";
	cout<<"computing e+e-e+ for xk_cut: "<<xk_cut<<" MODE: "<<MODE_NEW<<endl;
	double res_eee_pt_MC = MonteCarlo_integration_extended_phase_space(H1[ijack], H2[ijack], FA[ijack], FV[ijack], m, fp, "e+e-", xk_cut, 1, RM(), MODE_NEW);
	Plot_cut_rate_eee_pt.distr_list[ip].distr.push_back(res_eee_pt_MC);
	//MODE_NEW="SD";
	//cout<<"computing  e+e-e+ for xk_cut: "<<xk_cut<<" MODE: "<<MODE_NEW<<endl;
	//double res_eee_SD_MC = MonteCarlo_integration_extended_phase_space(H1[ijack], H2[ijack], FA[ijack], FV[ijack], m, fp, "e+e-", xk_cut, 1, RM(), MODE_NEW);
	Plot_cut_rate_eee_SD.distr_list[ip].distr.push_back(res_eee_tot_MC-res_eee_pt_MC);
      }
    }



    
    //redo e+e-

    double res_ee_MC_extended= MonteCarlo_integration_extended_phase_space(H1[ijack], H2[ijack], FA[ijack], FV[ijack], m, fp, "e+e-", sqrt(cut), 0, RM(), MODE);
    cout<<"Relative error achieved in Monte Carlo e+e- (Extended) for jack: "<<ijack<<" is: "<<F_same_lepton*error_goal_MC_ee_extended/res_ee_MC_extended<<endl;
    cout<<"Branching ratio Monte Carlo e+e- (Extended): "<<res_ee_MC_extended<<endl;
    cout<<"################"<<endl;

    

   

   


    //mu+mu-mu+

    double res_mumumu_MC = MonteCarlo_integration_extended_phase_space(H1[ijack], H2[ijack], FA[ijack], FV[ijack], m, fp, "mu+mu-", sqrt(cut_mumu), 1, RM(), MODE);

    Res.Jack_Distr_Int_MonteCarlo_mumumu.distr.push_back(res_mumumu_MC);
    Res.Stat_err_MonteCarlo_mumumu.push_back(error_goal_MC_mumumu/res_mumumu_MC);
    cout<<"Relative error achieved in Monte Carlo mu+mu-mu+ for jack: "<<ijack<<" is: "<<F_same_lepton*error_goal_MC_mumumu/res_mumumu_MC<<endl;
    cout<<"Branching ratio Monte Carlo mu+mu-mu+: "<<res_mumumu_MC<<endl;
    cout<<"################"<<endl;


    
    //store total rate mumumu as a function of the inf cut for eee
    if(Print_rate) {
      for(int ip=0; ip<Npoints_c;ip++) {
	double xk_cut = off_c_mu+ ip*step_p_c_mu;
	string MODE_NEW="TOTAL";
	cout<<"computing mu+mu-mu+ for xk_cut: "<<xk_cut<<" MODE: "<<MODE_NEW<<endl;
	double res_mumumu_tot_MC = MonteCarlo_integration_extended_phase_space(H1[ijack], H2[ijack], FA[ijack], FV[ijack], m, fp, "mu+mu-", xk_cut, 1, RM(), MODE_NEW);
	Plot_cut_rate_mumumu_tot.distr_list[ip].distr.push_back(res_mumumu_tot_MC);
	MODE_NEW="PT";
	cout<<"computing mu+mu-mu+ for xk_cut: "<<xk_cut<<" MODE: "<<MODE_NEW<<endl;
	double res_mumumu_pt_MC = MonteCarlo_integration_extended_phase_space(H1[ijack], H2[ijack], FA[ijack], FV[ijack], m, fp, "mu+mu-", xk_cut, 1, RM(), MODE_NEW);
	Plot_cut_rate_mumumu_pt.distr_list[ip].distr.push_back(res_mumumu_pt_MC);
	//MODE_NEW="SD";
	//cout<<"computing mu+mu-mu+ for xk_cut: "<<xk_cut<<" MODE: "<<MODE_NEW<<endl;
	//double res_mumumu_SD_MC = MonteCarlo_integration_extended_phase_space(H1[ijack], H2[ijack], FA[ijack], FV[ijack], m, fp, "mu+mu-", xk_cut, 1, RM(), MODE_NEW);
	Plot_cut_rate_mumumu_SD.distr_list[ip].distr.push_back(res_mumumu_tot_MC- res_mumumu_pt_MC);
      }
    }


    

    //redo mu+mu-

    double res_mumu_MC_extended= MonteCarlo_integration_extended_phase_space(H1[ijack], H2[ijack], FA[ijack], FV[ijack], m, fp, "mu+mu-", sqrt(cut_mumu), 0, RM(), MODE);
    cout<<"Relative error achieved in Monte Carlo mu+mu- (Extended) for jack: "<<ijack<<" is: "<<F_same_lepton*error_goal_MC_mumu_extended/res_mumu_MC_extended<<endl;
    cout<<"Branching ratio Monte Carlo mu+mu- (Extended): "<<res_mumu_MC_extended<<endl;
    cout<<"################"<<endl;


    

    



  }


  //compute expectation value and std err
  Res.Int_Quad_val_ee= Res.Jack_Distr_Int_Quad_ee.ave();
  Res.Int_Quad_val_mumu= Res.Jack_Distr_Int_Quad_mumu.ave();
  Res.Int_MonteCarlo_val_ee = Res.Jack_Distr_Int_MonteCarlo_ee.ave();
  Res.Int_MonteCarlo_val_mumu= Res.Jack_Distr_Int_MonteCarlo_mumu.ave();
  Res.Int_Quad_err_ee= Res.Jack_Distr_Int_Quad_ee.err();
  Res.Int_Quad_err_mumu= Res.Jack_Distr_Int_Quad_mumu.err();
  Res.Int_MonteCarlo_err_ee = Res.Jack_Distr_Int_MonteCarlo_ee.err();
  Res.Int_MonteCarlo_err_mumu= Res.Jack_Distr_Int_MonteCarlo_mumu.err();
  Res.Int_MonteCarlo_val_eee = Res.Jack_Distr_Int_MonteCarlo_eee.ave();
  Res.Int_MonteCarlo_val_mumumu = Res.Jack_Distr_Int_MonteCarlo_mumumu.ave();
  Res.Int_MonteCarlo_err_eee = Res.Jack_Distr_Int_MonteCarlo_eee.err();
  Res.Int_MonteCarlo_err_mumumu = Res.Jack_Distr_Int_MonteCarlo_mumumu.err();


 
  //plot diff rate for e+e- and mu+mu-
  if(Print_rate) {
  boost::filesystem::create_directory("../data/form_factors/"+Meson+"/diff_rate");
  Print_To_File({}, {p_list_d_e, Plot_diff_rate_ee_tot.ave(), Plot_diff_rate_ee_tot.err(), Plot_diff_rate_ee_SD.ave(), Plot_diff_rate_ee_SD.err(), Plot_diff_rate_ee_Q.ave(), Plot_diff_rate_ee_Q.err(), Plot_diff_rate_ee_pt.ave(), Plot_diff_rate_ee_pt.err(), Plot_diff_rate_ee_int.ave(), Plot_diff_rate_ee_int.err()}, "../data/form_factors/"+Meson+"/diff_rate/diff_rate_ee.dat", "", "#xk    tot    sd  quadratic   pt     int");
  Print_To_File({}, {p_list_d_mu, Plot_diff_rate_mumu_tot.ave(), Plot_diff_rate_mumu_tot.err(), Plot_diff_rate_mumu_SD.ave(), Plot_diff_rate_mumu_SD.err(), Plot_diff_rate_mumu_Q.ave(), Plot_diff_rate_mumu_Q.err(), Plot_diff_rate_mumu_pt.ave(), Plot_diff_rate_mumu_pt.err(), Plot_diff_rate_mumu_int.ave(), Plot_diff_rate_mumu_int.err()}, "../data/form_factors/"+Meson+"/diff_rate/diff_rate_mumu.dat", "", "#xk   tot    sd  quadratic   pt    int");
  }


   if(Print_rate) {
  //plot total rate for e+e-e+ and mu+mu-mu+ as a function of the lower cut
   boost::filesystem::create_directory("../data/form_factors/"+Meson+"/tot_rate_lcut");
   Print_To_File({}, {p_list_c_e, Plot_cut_rate_eee_tot.ave(), Plot_cut_rate_eee_tot.err(), Plot_cut_rate_eee_SD.ave(), Plot_cut_rate_eee_SD.err(), Plot_cut_rate_eee_pt.ave(), Plot_cut_rate_eee_pt.err()}, "../data/form_factors/"+Meson+"/tot_rate_lcut/rate_eee.dat", "", "#xk    tot    sd     pt");
   Print_To_File({}, {p_list_c_mu, Plot_cut_rate_mumumu_tot.ave(), Plot_cut_rate_mumumu_tot.err(), Plot_cut_rate_mumumu_SD.ave(), Plot_cut_rate_mumumu_SD.err(), Plot_cut_rate_mumumu_pt.ave(), Plot_cut_rate_mumumu_pt.err()}, "../data/form_factors/"+Meson+"/tot_rate_lcut/rate_mumumu.dat", "", "#xk   tot    sd     pt");
   }
  


  return Res;
}
