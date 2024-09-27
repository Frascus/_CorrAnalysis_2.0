#include "../include/chi_mag.h"
#include "g_minus_2_utilities.h"
#include "numerics.h"
using namespace std;

const int Njacks = 50;
const bool UseJack = 1;
bool Get_ASCII = false;
const double fmTGeV= 1.0/0.197327;
const double QCD_scale = 0.5 * fmTGeV;
const bool Gen_free_corr_data = false;
const bool Use_tree_level_sub = false;
const double Qp = 2.0 / 3.0;
const double Qn = -1.0 / 3.0;
const bool Include_disco = true;
//const bool Include_sea_quark_mass_derivative=true;


double Ker_sub(double x) {

  return (sinh(x) - x);
  
  double sum=0;
  double sum_p;
  bool converged= false;
  double  prec= 1e-10;
  double n=3;
  double fact=-x;
  while(!converged) {
    sum_p=sum;
    fact *= -1*pow(x,2)/(n*(n-1));
    sum += fact;
    if( fabs(sum- sum_p) < prec) converged=true;
    n+= 2;
  }

  return sum;
}

distr_t Ker_sub_distr(const distr_t &x) {

  double Nj= x.size();
  distr_t ret(UseJack, Nj);
  for(int ij=0;ij<Nj;ij++) ret.distr[ij] = Ker_sub(x.distr[ij]);

  
  return ret;
}

double susc_pert( double p, double m) {

  m= 2*m;

  double Nc=3;
  double val, err;
  auto Kv = [&Nc, &p, &m](double x) { return (Nc/(2*M_PI*M_PI))*( cosh(p*x/(2*m)) - sinh(p*x/(2*m))/(p/(2*m)))*boost::math::cyl_bessel_k(1,x);};

  gsl_function_pp<decltype(Kv)> F_corr(Kv);
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (10000);
  gsl_function *G = static_cast<gsl_function*>(&F_corr);
  gsl_integration_qagiu(G, 1e-6*(2*m/p), 0.0, 1e-8, 10000, w, &val, &err);
  gsl_integration_workspace_free (w);
  

  return val*1000*m;

}


double Get_tree_lev_der(double x) {


  auto Fb = [](double x) {

    auto Kv = [](double n) { return boost::math::cyl_bessel_k(1,n);};


    double val, err;
    gsl_function_pp<decltype(Kv)> F_corr(Kv);
		  gsl_integration_workspace * w = gsl_integration_workspace_alloc (10000);
		  gsl_function *G = static_cast<gsl_function*>(&F_corr);
		  gsl_integration_qagiu(G, x, 0.0, 1e-8, 10000, w, &val, &err);
		  gsl_integration_workspace_free (w);

    
		  return val;
      };
  
  double Nc=3;
  double F1=  (Nc/(2*M_PI*M_PI))*Fb(x);
  double F2=  (Nc/(2*M_PI*M_PI))*x*boost::math::cyl_bessel_k(1,x);

  return F1-F2;

}



void Generate_free_corr_data_VT() {

  vector<tuple<double, int, int>> amu_Tmax;

 amu_Tmax.push_back( tuple<double,int,int>( 0.01, 96,48));
 amu_Tmax.push_back( tuple<double,int,int>( 0.01, 200,100));
 amu_Tmax.push_back( tuple<double,int,int>( 0.01, 400,200));
 amu_Tmax.push_back( tuple<double,int,int>( 0.01, 100,-1));

 
  //strang

  //B ensembles

 /*

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




  //light

  amu_Tmax.push_back( make_pair( 0.00054, 96));
  amu_Tmax.push_back( make_pair( 0.00060, 80));
  amu_Tmax.push_back( make_pair( 0.00072, 64));
  amu_Tmax.push_back( make_pair( 0.00072, 96));

 */
    

 for(auto &p: amu_Tmax) { cout<<"Executing amu: "<<get<0>(p)<<" Tmax: "<<get<1>(p)<<" L: "<<get<2>(p)<<endl; Compute_free_VT_corr(get<0>(p), get<1>(p), get<2>(p));}

}


void Compute_free_VT_corr(double am, int Tmax, int L) {


  int Tmax_n = (L==-1)?Tmax:2*L;

  int fwbw = 0;

  //create directory
  boost::filesystem::create_directory("../data/magnetic_susc/free_theory");

 
  double Nc=3; //three colors

  auto C_cont = [&Nc, &am](int t) -> double {

    //double tolerance=1e-16;
    //double err;

    auto f = [&am, &t, &Nc](double x) {  return (Nc*4.0/pow(M_PI,2))*exp(-2.0*t*sqrt( pow(x,2) + pow(am,2)))*(am)*pow(x,2)/(4*sqrt((pow(x,2) + pow(am,2))));};

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

  auto corr = [&am, &Nc, &ptm2, &D2, &W, &shaEp2, &Ep, &fwbw, &Tmax_n](double p1,double p2,double p3, double t, int r) -> double {

    return (8*8*am*Nc/(pow(2.0*M_PI,3)))*(exp(-2*Ep(p1,p2,p3)*t) + fwbw*exp(-2*Ep(p1,p2,p3)*(Tmax_n-t)))*( sqrt(shaEp2(p1,p2,p3)))/D2(p1,p2,p3);
	      };

  vector<double> corr_pert_res(2*Tmax_n,0.0);

  int r=1;
   

  int time;
  double tol = 1e-9;
    
  for(int t=1;t<=Tmax_n;t++) { //loop over time


    if(L==-1) cout<<" t: "<<t<<endl;
    if(L != -1) cout<<"t: "<<t<<" L: "<<L<<endl;
    time =t;

    if( L == -1) {
    double err;
    
    auto corrp1= [&corr, &time, &r, &tol](double p1) -> double {  //gsl only performs 1d integrals
     
      auto corrp2 = [&corr, &p1, &time, &r, &tol](double p2) {
	
	auto corrp3 = [&corr, &p1, &p2, &time, &r](double p3) {
	  return corr(p1, p2, p3, time,r);
	};
	
	
	double err_3;
	double tol3= 1e-9;
	double val_3;
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
      
      
      corr_pert_res[t] = val;

    }
    else {

      for(int iP1=0;iP1<=L/2;iP1++)
	for(int iP2=0;iP2<=iP1;iP2++)
	  for(int iP3=0;iP3<=iP2;iP3++)
	    {
	      /// Space parity
	      const int parMult=
		((iP1 and iP1!=L/2)+1)*
		((iP2 and iP2!=L/2)+1)*
		((iP3 and iP3!=L/2)+1);
	      
	      /// Multiplicity of the permutation, given the number of equal
	      /// components (1 is actually never used)
	      constexpr double permMultTable[]=
		{1,1,3,6};
	      
	      /// Space components multiplicity
	      const int permMult=
		permMultTable
		[(iP1!=iP2)+
		 (iP2!=iP3)+
		 (iP3!=iP1)];

	      corr_pert_res[t] += parMult*permMult*pow((2*M_PI/L),3)*corr(2*M_PI*iP1/L, 2*M_PI*iP2/L, 2*M_PI*iP3/L, t, 1)/8;

	    }
    }
      
  }
  



  //compute correlator in the continuum limit
  vector<double> corr_pert_cont(2*Tmax_n,0.0);
  for(int t=1;t<=Tmax_n;t++) corr_pert_cont[t] = C_cont(t);


  //take difference between corr_pert_cont and corr_pert_res_tm(OS)
  vector<double> a2corr_pert_res(2*Tmax_n,0.0);
  vector<double> susc_lat_integrand(2*Tmax_n,0.0);
  vector<double> susc_cont_integrand(2*Tmax_n,0.0);
  vector<double> susc_lat(2*Tmax_n,0.0);
  vector<double> susc_cont(2*Tmax_n,0.0);
   
  for(int t=0;t<2*Tmax_n;t++) {
    a2corr_pert_res[t] = corr_pert_cont[t] - corr_pert_res[t];
    susc_lat_integrand[t] = corr_pert_res[t]*t;
    susc_cont_integrand[t] = corr_pert_cont[t]*t;
    if( t > 0) {
    susc_lat[t] = susc_lat[t-1] + corr_pert_res[t]*t;
    susc_cont[t] = susc_cont[t-1] + corr_pert_cont[t]*t;
    }
  }
 
  //Print the result

 

  if(L==-1) {
  Print_To_File({}, {a2corr_pert_res, corr_pert_res, corr_pert_cont}, "../data/magnetic_susc/free_theory/corr_T"+to_string(Tmax_n)+"_m_"+to_string_with_precision(am,3), "" , "");
  Print_To_File({}, {susc_lat_integrand,susc_cont_integrand, susc_lat, susc_cont}, "../data/magnetic_susc/free_theory/susc_T_"+to_string(Tmax_n)+"_m_"+to_string_with_precision(am,3), "" , "");
  }
  else {

    Print_To_File({}, {a2corr_pert_res, corr_pert_res, corr_pert_cont}, "../data/magnetic_susc/free_theory/corr_T"+to_string(Tmax_n)+"_L_"+to_string(L)+"_m_"+to_string_with_precision(am,3), "" , "");
    Print_To_File({}, {susc_lat_integrand,susc_cont_integrand, susc_lat, susc_cont}, "../data/magnetic_susc/free_theory/susc_T_"+to_string(Tmax_n)+"_L_"+to_string(L)+"_m_"+to_string_with_precision(am,3), "" , "");

  }



  return;
}


void Compute_magnetic_susc() {


  //read binary files
    boost::filesystem::create_directory("../new_masses_strange_B64");
    

    vector<string> Ens_T1({"B.72.64"});
    vector<string> Ens_TT1({ "cB211b.072.64"});

    for( int it=0; it<(signed)Ens_T1.size(); it++) {

      vector<string> channels({"mix_s1_s1", "mix_s2_s2"});

      for(auto &channel : channels) {
	boost::filesystem::create_directory("../new_masses_strange_B64/"+channel);
	boost::filesystem::create_directory("../new_masses_strange_B64/"+channel+"/"+Ens_TT1[it]);
      }
      //read binary
      vector<string> Corr_tags({"TM_VKVK", "TM_AKAK", "TM_A0A0", "TM_V0V0", "TM_P5P5", "TM_S0S0", "TM_TKVK", "TM_VKTK", "TM_TKTK", "OS_VKVK", "OS_AKAK", "OS_A0A0", "OS_V0V0", "OS_P5P5", "OS_S0S0", "OS_TKVK", "OS_VKTK", "OS_TKTK"});

          
      for(int id=0; id<(signed)Corr_tags.size(); id++) {



	for( auto &channel: channels) {

	FILE *stream = fopen( ("../gm2_tau_rep_bin/"+Ens_T1[it]+"/"+channel+"_"+Corr_tags[id]).c_str(), "rb");
        size_t Nconfs, T, Nhits;
	bin_read(Nconfs, stream);
	bin_read(Nhits, stream);
	bin_read(T, stream);

	cout<<"Nconfs: "<<Nconfs<<endl;
	cout<<"T: "<<T<<" "<<T/2+1<<endl;
	cout<<"Nhits: "<<Nhits<<endl;
	for(size_t iconf=0;iconf<Nconfs;iconf++) {
	  vector<double> C(T/2+1);
	  for(size_t t=0;t<T/2+1;t++) bin_read(C[t], stream);
	  boost::filesystem::create_directory("../new_masses_strange_B64/"+channel+"/"+Ens_TT1[it]+"/"+to_string(iconf));
	  ofstream PrintCorr("../new_masses_strange_B64/"+channel+"/"+Ens_TT1[it]+"/"+to_string(iconf)+"/mes_contr_"+channel+"_"+Corr_tags[id]);
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

  Get_magnetic_susc(false);
  Get_magnetic_susc(true);


}

void Get_magnetic_susc(bool Include_sea_quark_mass_derivative) {


  string Tag_val=  (Include_sea_quark_mass_derivative?"":"_only_val");

  if(Gen_free_corr_data) Generate_free_corr_data_VT();

  GaussianMersenne GM(356923432);

  //define lambda function to combine FF and BB

  //resample RCs
  distr_t ZA_A(UseJack), ZA_B(UseJack), ZA_C(UseJack), ZA_D(UseJack);
  distr_t ZV_A(UseJack), ZV_B(UseJack), ZV_C(UseJack), ZV_D(UseJack);
  distr_t ZT_A(UseJack), ZT_B(UseJack), ZT_C(UseJack), ZT_D(UseJack);
  distr_t a_A(UseJack), a_B(UseJack), a_C(UseJack), a_D(UseJack);

 
  double fmTGeV= 1.0/0.197327;

  LatticeInfo L_info_A, L_info_B, L_info_C, L_info_D;
  L_info_A.LatInfo_new_ens("cA211a.12.48");
  L_info_B.LatInfo_new_ens("cB211b.072.96");
  L_info_C.LatInfo_new_ens("cC211a.06.80");
  L_info_D.LatInfo_new_ens("cD211a.054.96");
  

  for(int ijack=0; ijack<Njacks;ijack++) {

    ZA_A.distr.push_back( L_info_A.Za_WI_strange + GM()*L_info_A.Za_WI_strange_err/((UseJack==true)?sqrt(Njacks -1.0):1.0));
    ZV_A.distr.push_back( L_info_A.Zv_WI_strange + GM()*L_info_A.Zv_WI_strange_err/((UseJack==true)?sqrt(Njacks -1.0):1.0));

    ZA_B.distr.push_back( L_info_B.Za_WI_strange + GM()*L_info_B.Za_WI_strange_err/((UseJack==true)?sqrt(Njacks -1.0):1.0));
    ZV_B.distr.push_back( L_info_B.Zv_WI_strange + GM()*L_info_B.Zv_WI_strange_err/((UseJack==true)?sqrt(Njacks -1.0):1.0));

    ZA_C.distr.push_back( L_info_C.Za_WI_strange + GM()*L_info_C.Za_WI_strange_err/((UseJack==true)?sqrt(Njacks -1.0):1.0));
    ZV_C.distr.push_back( L_info_C.Zv_WI_strange + GM()*L_info_C.Zv_WI_strange_err/((UseJack==true)?sqrt(Njacks -1.0):1.0));

    ZA_D.distr.push_back( L_info_D.Za_WI_strange + GM()*L_info_D.Za_WI_strange_err/((UseJack==true)?sqrt(Njacks -1.0):1.0));
    ZV_D.distr.push_back( L_info_D.Zv_WI_strange + GM()*L_info_D.Zv_WI_strange_err/((UseJack==true)?sqrt(Njacks -1.0):1.0));

    a_A.distr.push_back( L_info_A.a_from_afp*fmTGeV + GM()*L_info_A.a_from_afp_err*fmTGeV/((UseJack==true)?sqrt(Njacks -1.0):1.0));
    a_B.distr.push_back( L_info_B.a_from_afp*fmTGeV + GM()*L_info_B.a_from_afp_err*fmTGeV/((UseJack==true)?sqrt(Njacks -1.0):1.0));
    a_C.distr.push_back( L_info_C.a_from_afp*fmTGeV + GM()*L_info_C.a_from_afp_err*fmTGeV/((UseJack==true)?sqrt(Njacks -1.0):1.0));
    a_D.distr.push_back( L_info_D.a_from_afp*fmTGeV + GM()*L_info_D.a_from_afp_err*fmTGeV/((UseJack==true)?sqrt(Njacks -1.0):1.0));

       
    ZT_A.distr.push_back( L_info_A.ZT_RI2 + GM()*L_info_A.ZT_RI2_err/((UseJack==true)?sqrt(Njacks -1.0):1.0));
    ZT_B.distr.push_back( L_info_B.ZT_RI2 + GM()*L_info_B.ZT_RI2_err/((UseJack==true)?sqrt(Njacks -1.0):1.0));
    ZT_C.distr.push_back( L_info_C.ZT_RI2 + GM()*L_info_C.ZT_RI2_err/((UseJack==true)?sqrt(Njacks -1.0):1.0));
    ZT_D.distr.push_back( L_info_D.ZT_RI2 + GM()*L_info_D.ZT_RI2_err/((UseJack==true)?sqrt(Njacks -1.0):1.0));
    

  }


  cout<<"Lattice spacings used: "<<endl;
  cout<<"A: "<< (a_A/fmTGeV).ave()<<" +- "<<(a_A/fmTGeV).err()<<endl;
  cout<<"B: "<< (a_B/fmTGeV).ave()<<" +- "<<(a_B/fmTGeV).err()<<endl;
  cout<<"C: "<< (a_C/fmTGeV).ave()<<" +- "<<(a_C/fmTGeV).err()<<endl;
  cout<<"D: "<< (a_D/fmTGeV).ave()<<" +- "<<(a_D/fmTGeV).err()<<endl;
  cout<<"Renormalization constants used: "<<endl;
  cout<<"A: "<< (ZT_A).ave()<<" +- "<<(ZT_A).err()<<endl;
  cout<<"B: "<< (ZT_B).ave()<<" +- "<<(ZT_B).err()<<endl;
  cout<<"C: "<< (ZT_C).ave()<<" +- "<<(ZT_C).err()<<endl;
  cout<<"D: "<< (ZT_D).ave()<<" +- "<<(ZT_D).err()<<endl;
  

  
  if(Get_ASCII) {
    //read binary files
    boost::filesystem::create_directory("../magnetic_susc_data");
    

    vector<string> Ens_T1({"B.72.64", "B.72.96", "C.06.80", "D.54.96"});
    vector<string> Ens_TT1({"cB211b.072.64", "cB211b.072.96", "cC211a.06.80", "cD211a.054.96"});

    for( int it=0; it<(signed)Ens_T1.size(); it++) {

      vector<string> channels({"ll", "ss1", "ss2"});

      for(auto &channel : channels) {
	boost::filesystem::create_directory("../magnetic_susc_data/"+channel);
	boost::filesystem::create_directory("../magnetic_susc_data/"+channel+"/"+Ens_TT1[it]);
      }
      //read binary
      vector<string> Corr_tags({"TM_TKTK", "TM_TKVK", "TM_VKTK", "TM_P5P5", "TM_VKVK", "OS_TKTK", "OS_TKVK", "OS_VKTK", "OS_P5P5", "OS_VKVK"});

 
      
      for(int id=0; id<(signed)Corr_tags.size(); id++) {



	for( auto &channel: channels) {

	FILE *stream = fopen( ("../magnetic_susc_bin/magn_susc_data/"+Ens_T1[it]+"/data/"+channel+"_"+Corr_tags[id]).c_str(), "rb");
        size_t Nconfs, T, Nhits;
	bin_read(Nconfs, stream);
	bin_read(Nhits, stream);
	bin_read(T, stream);

	cout<<"Nconfs: "<<Nconfs<<endl;
	cout<<"T: "<<T<<" "<<T/2+1<<endl;
	cout<<"Nhits: "<<Nhits<<endl;
	ifstream Read_confs_tag( "../magnetic_susc_bin/magn_susc_data/"+Ens_T1[it]+"/data/confsList"+channel+".txt");
	for(size_t iconf=0;iconf<Nconfs;iconf++) {
	  vector<double> C(T/2+1);
	  for(size_t t=0;t<T/2+1;t++) bin_read(C[t], stream);
	  string Confs_tag;
	  Read_confs_tag>>Confs_tag;
	  boost::filesystem::create_directory("../magnetic_susc_data/"+channel+"/"+Ens_TT1[it]+"/"+Confs_tag);
	  ofstream PrintCorr("../magnetic_susc_data/"+channel+"/"+Ens_TT1[it]+"/"+Confs_tag+"/mes_contr_"+channel+"_"+Corr_tags[id]);
	  PrintCorr.precision(16);
	  PrintCorr<<"# "<<Corr_tags[id].substr(3,4)<<endl;
	  for(size_t t=0;t<(T/2+1);t++) PrintCorr<<C[t]<<endl;
	  if(Corr_tags[id].substr(3,4) == "VKTK" || Corr_tags[id].substr(3,4) == "TKVK") { for(size_t t=T/2+1; t<T;t++) PrintCorr<<-1*C[T-t]<<endl;   }
	  else  {for(size_t t=T/2+1; t<T;t++) PrintCorr<<C[T-t]<<endl;}
	  PrintCorr.close();

	}
	Read_confs_tag.close();

	fclose(stream);

	}
	
      }
      //eliminate unmatching confs between ss1 and ss2
      ifstream Read_confs_ss1( "../magnetic_susc_bin/magn_susc_data/"+Ens_T1[it]+"/data/confsListss1.txt");
      ifstream Read_confs_ss2( "../magnetic_susc_bin/magn_susc_data/"+Ens_T1[it]+"/data/confsListss2.txt");

      while(!Read_confs_ss1.eof()) {
	string conf_ss1_tag;
	Read_confs_ss1>>conf_ss1_tag;
	if(!Read_confs_ss1.eof()) {
	  if( ! boost::filesystem::exists( "../magnetic_susc_data/ss2/"+Ens_TT1[it]+"/"+conf_ss1_tag) ) boost::filesystem::remove_all( "../magnetic_susc_data/ss1/"+Ens_TT1[it]+"/"+conf_ss1_tag);
	}
      }

      while(!Read_confs_ss2.eof()) {
	string conf_ss2_tag;
	Read_confs_ss2>>conf_ss2_tag;
	if(!Read_confs_ss2.eof()) {
	  if( ! boost::filesystem::exists( "../magnetic_susc_data/ss1/"+Ens_TT1[it]+"/"+conf_ss2_tag) ) boost::filesystem::remove_all( "../magnetic_susc_data/ss2/"+Ens_TT1[it]+"/"+conf_ss2_tag);
	}
      }
      
      Read_confs_ss1.close();
      Read_confs_ss2.close();
    }
  }

  
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





  //light
  data_t data_TKTK_tm, data_TKVK_tm, data_VKTK_tm, data_P5P5_tm;
  data_t data_TKTK_OS, data_TKVK_OS, data_VKTK_OS, data_P5P5_OS;

  //light to estimate der
  data_t data_TKVK_tm_m1, data_VKTK_tm_m1, data_TKVK_OS_m1, data_VKTK_OS_m1;
  data_t data_TKVK_tm_m2, data_VKTK_tm_m2, data_TKVK_OS_m2, data_VKTK_OS_m2;

  //strange 1 and 2
  data_t data_s1_TKTK_tm, data_s1_TKVK_tm, data_s1_VKTK_tm, data_s1_P5P5_tm;
  data_t data_s1_TKTK_OS, data_s1_TKVK_OS, data_s1_VKTK_OS, data_s1_P5P5_OS;

  data_t data_s2_TKTK_tm, data_s2_TKVK_tm, data_s2_VKTK_tm, data_s2_P5P5_tm;
  data_t data_s2_TKTK_OS, data_s2_TKVK_OS, data_s2_VKTK_OS, data_s2_P5P5_OS;

  //strange 3 and 4 for B64
  data_t data_s3_TKTK_tm, data_s3_TKVK_tm, data_s3_VKTK_tm, data_s3_P5P5_tm;
  data_t data_s3_TKTK_OS, data_s3_TKVK_OS, data_s3_VKTK_OS, data_s3_P5P5_OS;

  data_t data_s4_TKTK_tm, data_s4_TKVK_tm, data_s4_VKTK_tm, data_s4_P5P5_tm;
  data_t data_s4_TKTK_OS, data_s4_TKVK_OS, data_s4_VKTK_OS, data_s4_P5P5_OS;
  

  //disconnected
  data_t data_disc_VKTK_light, data_disc_VKTK_strange, data_disc_JJ_light;
  data_t data_disc_V1T1_light, data_disc_V2T2_light, data_disc_V3T3_light;


  data_t data_barXX_light, data_barXX_strange;
  data_t data_VKTK_light_tm_correlated, data_TKVK_light_tm_correlated, data_VKTK_light_OS_correlated, data_TKVK_light_OS_correlated;
  data_t data_VKTK_s1_tm_correlated, data_TKVK_s1_tm_correlated, data_VKTK_s1_OS_correlated, data_TKVK_s1_OS_correlated;
  data_t data_VKTK_s2_tm_correlated, data_TKVK_s2_tm_correlated, data_VKTK_s2_OS_correlated, data_TKVK_s2_OS_correlated;

  data_t  VKTK_barXX_tm_light, TKVK_barXX_tm_light, VKTK_barXX_OS_light, TKVK_barXX_OS_light;
  data_t  VKTK_barXX_tm_s1, TKVK_barXX_tm_s1, VKTK_barXX_OS_s1, TKVK_barXX_OS_s1;
  data_t  VKTK_barXX_tm_s2, TKVK_barXX_tm_s2, VKTK_barXX_OS_s2, TKVK_barXX_OS_s2;

  data_TKTK_tm.Read("../magnetic_susc_data/ll", "mes_contr_ll_TM_TKTK", "TKTK", Sort_light_confs);
  data_TKVK_tm.Read("../magnetic_susc_data/ll", "mes_contr_ll_TM_TKVK", "TKVK", Sort_light_confs);
  data_VKTK_tm.Read("../magnetic_susc_data/ll", "mes_contr_ll_TM_VKTK", "VKTK", Sort_light_confs);
  data_P5P5_tm.Read("../magnetic_susc_data/ll", "mes_contr_ll_TM_P5P5", "P5P5", Sort_light_confs);

  data_TKTK_OS.Read("../magnetic_susc_data/ll", "mes_contr_ll_OS_TKTK", "TKTK", Sort_light_confs);
  data_TKVK_OS.Read("../magnetic_susc_data/ll", "mes_contr_ll_OS_TKVK", "TKVK", Sort_light_confs);
  data_VKTK_OS.Read("../magnetic_susc_data/ll", "mes_contr_ll_OS_VKTK", "VKTK", Sort_light_confs);
  data_P5P5_OS.Read("../magnetic_susc_data/ll", "mes_contr_ll_OS_P5P5", "P5P5", Sort_light_confs);

  //to estimate der
  data_TKVK_tm_m1.Read("../magnetic_susc_data/valence/physical", "mes_contr_2pts_ll_1", "TKVK", Sort_light_confs);
  data_VKTK_tm_m1.Read("../magnetic_susc_data/valence/physical", "mes_contr_2pts_ll_1", "VKTK", Sort_light_confs);
  data_TKVK_OS_m1.Read("../magnetic_susc_data/valence/physical", "mes_contr_2pts_ll_2", "TKVK", Sort_light_confs);
  data_VKTK_OS_m1.Read("../magnetic_susc_data/valence/physical", "mes_contr_2pts_ll_2", "VKTK", Sort_light_confs);

  data_TKVK_tm_m2.Read("../magnetic_susc_data/valence/simulated", "mes_contr_2pts_ll_1", "TKVK", Sort_light_confs);
  data_VKTK_tm_m2.Read("../magnetic_susc_data/valence/simulated", "mes_contr_2pts_ll_1", "VKTK", Sort_light_confs);
  data_TKVK_OS_m2.Read("../magnetic_susc_data/valence/simulated", "mes_contr_2pts_ll_2", "TKVK", Sort_light_confs);
  data_VKTK_OS_m2.Read("../magnetic_susc_data/valence/simulated", "mes_contr_2pts_ll_2", "VKTK", Sort_light_confs);


  data_s1_TKTK_tm.Read("../magnetic_susc_data/ss1", "mes_contr_ss1_TM_TKTK", "TKTK", Sort_light_confs);
  data_s1_TKVK_tm.Read("../magnetic_susc_data/ss1", "mes_contr_ss1_TM_TKVK", "TKVK", Sort_light_confs);
  data_s1_VKTK_tm.Read("../magnetic_susc_data/ss1", "mes_contr_ss1_TM_VKTK", "VKTK", Sort_light_confs);
  data_s1_P5P5_tm.Read("../magnetic_susc_data/ss1", "mes_contr_ss1_TM_P5P5", "P5P5", Sort_light_confs);

  data_s1_TKTK_OS.Read("../magnetic_susc_data/ss1", "mes_contr_ss1_OS_TKTK", "TKTK", Sort_light_confs);
  data_s1_TKVK_OS.Read("../magnetic_susc_data/ss1", "mes_contr_ss1_OS_TKVK", "TKVK", Sort_light_confs);
  data_s1_VKTK_OS.Read("../magnetic_susc_data/ss1", "mes_contr_ss1_OS_VKTK", "VKTK", Sort_light_confs);
  data_s1_P5P5_OS.Read("../magnetic_susc_data/ss1", "mes_contr_ss1_OS_P5P5", "P5P5", Sort_light_confs);


  data_s2_TKTK_tm.Read("../magnetic_susc_data/ss2", "mes_contr_ss2_TM_TKTK", "TKTK", Sort_light_confs);
  data_s2_TKVK_tm.Read("../magnetic_susc_data/ss2", "mes_contr_ss2_TM_TKVK", "TKVK", Sort_light_confs);
  data_s2_VKTK_tm.Read("../magnetic_susc_data/ss2", "mes_contr_ss2_TM_VKTK", "VKTK", Sort_light_confs);
  data_s2_P5P5_tm.Read("../magnetic_susc_data/ss2", "mes_contr_ss2_TM_P5P5", "P5P5", Sort_light_confs);

  data_s2_TKTK_OS.Read("../magnetic_susc_data/ss2", "mes_contr_ss2_OS_TKTK", "TKTK", Sort_light_confs);
  data_s2_TKVK_OS.Read("../magnetic_susc_data/ss2", "mes_contr_ss2_OS_TKVK", "TKVK", Sort_light_confs);
  data_s2_VKTK_OS.Read("../magnetic_susc_data/ss2", "mes_contr_ss2_OS_VKTK", "VKTK", Sort_light_confs);
  data_s2_P5P5_OS.Read("../magnetic_susc_data/ss2", "mes_contr_ss2_OS_P5P5", "P5P5", Sort_light_confs);

  //
  
  data_s3_TKTK_tm.Read("../new_masses_strange_B64/mix_s1_s1", "mes_contr_mix_s1_s1_TM_TKTK", "TKTK", Sort_easy);
  data_s3_TKVK_tm.Read("../new_masses_strange_B64/mix_s1_s1", "mes_contr_mix_s1_s1_TM_TKVK", "TKVK", Sort_easy);
  data_s3_VKTK_tm.Read("../new_masses_strange_B64/mix_s1_s1", "mes_contr_mix_s1_s1_TM_VKTK", "VKTK", Sort_easy);
  data_s3_P5P5_tm.Read("../new_masses_strange_B64/mix_s1_s1", "mes_contr_mix_s1_s1_TM_P5P5", "P5P5", Sort_easy);

  data_s3_TKTK_OS.Read("../new_masses_strange_B64/mix_s1_s1", "mes_contr_mix_s1_s1_OS_TKTK", "TKTK", Sort_easy);
  data_s3_TKVK_OS.Read("../new_masses_strange_B64/mix_s1_s1", "mes_contr_mix_s1_s1_OS_TKVK", "TKVK", Sort_easy);
  data_s3_VKTK_OS.Read("../new_masses_strange_B64/mix_s1_s1", "mes_contr_mix_s1_s1_OS_VKTK", "VKTK", Sort_easy);
  data_s3_P5P5_OS.Read("../new_masses_strange_B64/mix_s1_s1", "mes_contr_mix_s1_s1_OS_P5P5", "P5P5", Sort_easy);


  data_s4_TKTK_tm.Read("../new_masses_strange_B64/mix_s2_s2", "mes_contr_mix_s2_s2_TM_TKTK", "TKTK", Sort_easy);
  data_s4_TKVK_tm.Read("../new_masses_strange_B64/mix_s2_s2", "mes_contr_mix_s2_s2_TM_TKVK", "TKVK", Sort_easy);
  data_s4_VKTK_tm.Read("../new_masses_strange_B64/mix_s2_s2", "mes_contr_mix_s2_s2_TM_VKTK", "VKTK", Sort_easy);
  data_s4_P5P5_tm.Read("../new_masses_strange_B64/mix_s2_s2", "mes_contr_mix_s2_s2_TM_P5P5", "P5P5", Sort_easy);

  data_s4_TKTK_OS.Read("../new_masses_strange_B64/mix_s2_s2", "mes_contr_mix_s2_s2_OS_TKTK", "TKTK", Sort_easy);
  data_s4_TKVK_OS.Read("../new_masses_strange_B64/mix_s2_s2", "mes_contr_mix_s2_s2_OS_TKVK", "TKVK", Sort_easy);
  data_s4_VKTK_OS.Read("../new_masses_strange_B64/mix_s2_s2", "mes_contr_mix_s2_s2_OS_VKTK", "VKTK", Sort_easy);
  data_s4_P5P5_OS.Read("../new_masses_strange_B64/mix_s2_s2", "mes_contr_mix_s2_s2_OS_P5P5", "P5P5", Sort_easy);

  

  

  //disconnected
  data_disc_VKTK_light.Read("../magnetic_susc_disco/disco", "light_TKVK.txt", "TKVK", Sort_light_confs);
  data_disc_V1T1_light.Read("../magnetic_susc_disco/disco", "light_T1V1.txt", "T1V1", Sort_light_confs);
  data_disc_V2T2_light.Read("../magnetic_susc_disco/disco", "light_T2V2.txt", "T2V2", Sort_light_confs);
  data_disc_V3T3_light.Read("../magnetic_susc_disco/disco", "light_T3V3.txt", "T3V3", Sort_light_confs);
  data_disc_JJ_light.Read("../magnetic_susc_disco/disco", "JJ_light.txt", "VKVK", Sort_light_confs);
  data_disc_VKTK_strange.Read("../magnetic_susc_disco/disco", "strange_TKVK.txt", "TKVK", Sort_light_confs);
 


  //chiral condensate
  data_barXX_light.Read("../magnetic_susc_disco/conn_XX_correlated/ll", "scalar_light.txt", "S", Sort_light_confs);
  data_barXX_strange.Read("../magnetic_susc_disco/conn_XX_correlated/ss1", "scalar_strange.txt", "S", Sort_light_confs);


  //connected TKVK correlated to chiral condensate
  data_VKTK_light_tm_correlated.Read("../magnetic_susc_disco/conn_XX_correlated/ll", "mes_contr_ll_TM_VKTK", "VKTK", Sort_light_confs);
  data_TKVK_light_tm_correlated.Read("../magnetic_susc_disco/conn_XX_correlated/ll", "mes_contr_ll_TM_TKVK", "TKVK", Sort_light_confs);
  data_VKTK_light_OS_correlated.Read("../magnetic_susc_disco/conn_XX_correlated/ll", "mes_contr_ll_OS_VKTK", "VKTK", Sort_light_confs);
  data_TKVK_light_OS_correlated.Read("../magnetic_susc_disco/conn_XX_correlated/ll", "mes_contr_ll_OS_TKVK", "TKVK", Sort_light_confs);


  data_VKTK_s1_tm_correlated.Read("../magnetic_susc_disco/conn_XX_correlated/ss1", "mes_contr_ss1_TM_VKTK", "VKTK", Sort_light_confs);
  data_TKVK_s1_tm_correlated.Read("../magnetic_susc_disco/conn_XX_correlated/ss1", "mes_contr_ss1_TM_TKVK", "TKVK", Sort_light_confs);
  data_VKTK_s1_OS_correlated.Read("../magnetic_susc_disco/conn_XX_correlated/ss1", "mes_contr_ss1_OS_VKTK", "VKTK", Sort_light_confs);
  data_TKVK_s1_OS_correlated.Read("../magnetic_susc_disco/conn_XX_correlated/ss1", "mes_contr_ss1_OS_TKVK", "TKVK", Sort_light_confs);

  data_VKTK_s2_tm_correlated.Read("../magnetic_susc_disco/conn_XX_correlated/ss2", "mes_contr_ss2_TM_VKTK", "VKTK", Sort_light_confs);
  data_TKVK_s2_tm_correlated.Read("../magnetic_susc_disco/conn_XX_correlated/ss2", "mes_contr_ss2_TM_TKVK", "TKVK", Sort_light_confs);
  data_VKTK_s2_OS_correlated.Read("../magnetic_susc_disco/conn_XX_correlated/ss2", "mes_contr_ss2_OS_VKTK", "VKTK", Sort_light_confs);
  data_TKVK_s2_OS_correlated.Read("../magnetic_susc_disco/conn_XX_correlated/ss2", "mes_contr_ss2_OS_TKVK", "TKVK", Sort_light_confs);


  VKTK_barXX_tm_light.Read("../magnetic_susc_disco/conn_XX_correlated/ll", "ll_TM_VKTK_barXX_corr", "VKTK", Sort_light_confs);
  TKVK_barXX_tm_light.Read("../magnetic_susc_disco/conn_XX_correlated/ll", "ll_TM_TKVK_barXX_corr", "TKVK", Sort_light_confs);
  VKTK_barXX_OS_light.Read("../magnetic_susc_disco/conn_XX_correlated/ll", "ll_OS_VKTK_barXX_corr", "VKTK", Sort_light_confs);
  TKVK_barXX_OS_light.Read("../magnetic_susc_disco/conn_XX_correlated/ll", "ll_OS_TKVK_barXX_corr", "TKVK", Sort_light_confs);

  VKTK_barXX_tm_s1.Read("../magnetic_susc_disco/conn_XX_correlated/ss1", "ss1_TM_VKTK_barXX_corr", "VKTK", Sort_light_confs);
  TKVK_barXX_tm_s1.Read("../magnetic_susc_disco/conn_XX_correlated/ss1", "ss1_TM_TKVK_barXX_corr", "TKVK", Sort_light_confs);
  VKTK_barXX_OS_s1.Read("../magnetic_susc_disco/conn_XX_correlated/ss1", "ss1_OS_VKTK_barXX_corr", "VKTK", Sort_light_confs);
  TKVK_barXX_OS_s1.Read("../magnetic_susc_disco/conn_XX_correlated/ss1", "ss1_OS_TKVK_barXX_corr", "TKVK", Sort_light_confs);

  VKTK_barXX_tm_s2.Read("../magnetic_susc_disco/conn_XX_correlated/ss2", "ss2_TM_VKTK_barXX_corr", "VKTK", Sort_light_confs);
  TKVK_barXX_tm_s2.Read("../magnetic_susc_disco/conn_XX_correlated/ss2", "ss2_TM_TKVK_barXX_corr", "TKVK", Sort_light_confs);
  VKTK_barXX_OS_s2.Read("../magnetic_susc_disco/conn_XX_correlated/ss2", "ss2_OS_VKTK_barXX_corr", "VKTK", Sort_light_confs);
  TKVK_barXX_OS_s2.Read("../magnetic_susc_disco/conn_XX_correlated/ss2", "ss2_OS_TKVK_barXX_corr", "TKVK", Sort_light_confs);


  int Nens = data_TKTK_tm.size;

  //light
  //data where to store ensembles info
  distr_t_list a_distr_list(UseJack);
  distr_t_list susc_TV_tm_list(UseJack);
  distr_t_list susc_VT_tm_list(UseJack);
  distr_t_list susc_VTV_tm_list(UseJack);
  distr_t_list susc_TV_OS_list(UseJack);
  distr_t_list susc_VT_OS_list(UseJack);
  distr_t_list susc_VTV_OS_list(UseJack);
  //data where to store ensembles info red
  distr_t_list a_distr_list_red(UseJack);
  distr_t_list susc_TV_tm_list_red(UseJack);
  distr_t_list susc_VT_tm_list_red(UseJack);
  distr_t_list susc_VTV_tm_list_red(UseJack);
  distr_t_list susc_TV_OS_list_red(UseJack);
  distr_t_list susc_VT_OS_list_red(UseJack);
  distr_t_list susc_VTV_OS_list_red(UseJack);
       

  //strange
  //data where to store ensembles info
  distr_t_list susc_s_TV_tm_list(UseJack);
  distr_t_list susc_s_VT_tm_list(UseJack);
  distr_t_list susc_s_VTV_tm_list(UseJack);
  distr_t_list susc_s_TV_OS_list(UseJack);
  distr_t_list susc_s_VT_OS_list(UseJack);
  distr_t_list susc_s_VTV_OS_list(UseJack);
  //data where to store ensembles info red
  distr_t_list susc_s_TV_tm_list_red(UseJack);
  distr_t_list susc_s_VT_tm_list_red(UseJack);
  distr_t_list susc_s_VTV_tm_list_red(UseJack);
  distr_t_list susc_s_TV_OS_list_red(UseJack);
  distr_t_list susc_s_VT_OS_list_red(UseJack);
  distr_t_list susc_s_VTV_OS_list_red(UseJack);
      

  
  

  vector<string> Ensemble_list;
  vector<string> Ensemble_list_red;

  distr_t chiral_cond(UseJack);
  double val_cond= 1;
  double err_cond= 0;
  for(int ijack=0;ijack<Njacks;ijack++) chiral_cond.distr.push_back( val_cond + err_cond/sqrt(Njacks-1.0)); 
  //final values
  distr_t chi_light_TV(UseJack), chi_light_VT(UseJack), chi_light_VTV(UseJack);
  distr_t chi_strange_TV(UseJack), chi_strange_VT(UseJack), chi_strange_VTV(UseJack);
   

  boost::filesystem::create_directory("../data/magnetic_susc");

  for(int iens=0; iens<Nens;iens++) {
    
    boost::filesystem::create_directory("../data/magnetic_susc/"+data_TKTK_tm.Tag[iens]);

    LatticeInfo L_info;
    L_info.LatInfo_new_ens(data_TKTK_tm.Tag[iens]);
    double ml= L_info.ml;
    double ms1= L_info.ms_L;
    double ms2= L_info.ms_M;
    //define dm
    double dm=0;


    cout<<"Analyzing ensemble: "<<data_TKTK_tm.Tag[iens]<<endl;
    cout<<"aml: "<<ml<<" ams1: "<<ms1<<" ams2: "<<ms2<<endl;
    distr_t Za, Zv, Z_T, a_distr;
    if(data_TKTK_tm.Tag[iens].substr(1,1)=="A") {      Za= ZA_A; Zv=ZV_A; Z_T=ZT_A; a_distr=a_A; dm =0;}
    else if(data_TKTK_tm.Tag[iens].substr(1,1)=="B") { Za= ZA_B; Zv=ZV_B; Z_T=ZT_B; a_distr=a_B; dm= 0.00072-0.0006675;}
    else if(data_TKTK_tm.Tag[iens].substr(1,1)=="C") { Za= ZA_C; Zv=ZV_C; Z_T=ZT_C; a_distr=a_C; dm= 0.00060-0.000585;}
    else if(data_TKTK_tm.Tag[iens].substr(1,1)=="D") { Za= ZA_D; Zv=ZV_D; Z_T=ZT_D; a_distr=a_D; dm= 0.00054-0.0004964;}
    else crash("Ensemble: "+data_TKTK_tm.Tag[iens]+" not recognised");

    double dms= ms2 -ms1;
    double ms_ave= (ms1+ms2)/2;

    //evaluate two-pt function
    CorrAnalysis Corr(UseJack, Njacks, 800, iens);
    Corr.Nt = data_TKTK_tm.nrows[iens];
    Corr.Reflection_sign=1;
    Corr.Perform_Nt_t_average=1;

    //load chiral condensate and correlate it to TKVK
    auto F_disco_jack= [&](const Vfloat& par) { if((signed)par.size() != 3) crash("Lambda function F_disco_jack expects par[3], but par["+to_string((signed)par.size())+"] provided"); return par[0] -par[1]*par[2];};

    

    //sanity check
    if( (data_VKTK_tm.Tag[iens] != data_disc_VKTK_light.Tag[iens] )   || (data_VKTK_tm.Tag[iens] != data_barXX_light.Tag[iens] )) crash("Connected and disconnected ensemble list do not match");
    if( (data_VKTK_tm.Tag[iens] != data_disc_VKTK_strange.Tag[iens] )   || (data_VKTK_tm.Tag[iens] != data_barXX_strange.Tag[iens] )) crash("Connected and disconnected ensemble list do not match");


    //-------------------------   FOR SEA QUARK MASS DERIVATIVE ---------------------------//
    //-------   -------//
    int iens_XX=iens;
    int T_XX= Corr.Nt;
    if(data_VKTK_tm.Tag[iens] == "cB211b.072.96") {
      T_XX= 128;
      for(int j=0; j < (signed)data_VKTK_tm.Tag.size(); j++)
	if ( data_VKTK_tm.Tag[j] == "cB211b.072.64") iens_XX = j;
    }
 
    
    int Nconfs_barXX_light= data_barXX_light.col(1)[iens_XX][0].size();
    int Nconfs_barXX_strange= data_barXX_strange.col(1)[iens_XX][0].size();
    Vfloat barXX_light(Nconfs_barXX_light,0.0);
    Vfloat barXX_strange(Nconfs_barXX_strange,0.0);
  
    //final observerbable to evaluate sea-quark mass derivative
    distr_t_list VKTK_sea_der_tm_light(UseJack), TKVK_sea_der_tm_light(UseJack), VKTK_sea_der_OS_light(UseJack), TKVK_sea_der_OS_light(UseJack);
    distr_t_list VKTK_sea_der_tm_s1(UseJack), TKVK_sea_der_tm_s1(UseJack), VKTK_sea_der_OS_s1(UseJack), TKVK_sea_der_OS_s1(UseJack);
    distr_t_list VKTK_sea_der_tm_s2(UseJack), TKVK_sea_der_tm_s2(UseJack), VKTK_sea_der_OS_s2(UseJack), TKVK_sea_der_OS_s2(UseJack);

    for(int t=0;t<T_XX;t++) {
      for(int iconf=0; iconf<Nconfs_barXX_light;iconf++) barXX_light[iconf] += -1.0*data_barXX_light.col(1)[iens_XX][t][iconf];
      for(int iconf=0; iconf<Nconfs_barXX_strange;iconf++) barXX_strange[iconf] += -1.0*data_barXX_strange.col(1)[iens_XX][t][iconf];
    }

    for(int t=0; t < T_XX;t++) {
      

	//cout<<"SIZES ENS: "<<data_VKTK_tm.Tag[iens_XX]<<" : "<<VKTK_barXX_tm_light[t].size()<<" "<<barXX_light.size()<<" "<<data_VKTK_light_tm_correlated.col(0)[iens_XX][t].size()<<endl;
	//jackknife
       	Jackknife J_light(10000,Njacks);
	VKTK_sea_der_tm_light.distr_list.push_back(-dm*J_light.DoJack(F_disco_jack, 3, VKTK_barXX_tm_light.col(0)[iens_XX][t], barXX_light, data_VKTK_light_tm_correlated.col(0)[iens_XX][t]));
	TKVK_sea_der_tm_light.distr_list.push_back(-dm*J_light.DoJack(F_disco_jack, 3, TKVK_barXX_tm_light.col(0)[iens_XX][t], barXX_light, data_TKVK_light_tm_correlated.col(0)[iens_XX][t]));
	VKTK_sea_der_OS_light.distr_list.push_back(-dm*J_light.DoJack(F_disco_jack, 3, VKTK_barXX_OS_light.col(0)[iens_XX][t], barXX_light, data_VKTK_light_OS_correlated.col(0)[iens_XX][t]));
	TKVK_sea_der_OS_light.distr_list.push_back(-dm*J_light.DoJack(F_disco_jack, 3, TKVK_barXX_OS_light.col(0)[iens_XX][t], barXX_light, data_TKVK_light_OS_correlated.col(0)[iens_XX][t]));
	Jackknife J_s1(10000,Njacks);
	VKTK_sea_der_tm_s1.distr_list.push_back(-dms*J_s1.DoJack(F_disco_jack, 3, VKTK_barXX_tm_s1.col(0)[iens_XX][t], barXX_strange, data_VKTK_s1_tm_correlated.col(0)[iens_XX][t]));
	TKVK_sea_der_tm_s1.distr_list.push_back(-dms*J_s1.DoJack(F_disco_jack, 3, TKVK_barXX_tm_s1.col(0)[iens_XX][t], barXX_strange, data_TKVK_s1_tm_correlated.col(0)[iens_XX][t]));
	VKTK_sea_der_OS_s1.distr_list.push_back(-dms*J_s1.DoJack(F_disco_jack, 3, VKTK_barXX_OS_s1.col(0)[iens_XX][t], barXX_strange, data_VKTK_s1_OS_correlated.col(0)[iens_XX][t]));
	TKVK_sea_der_OS_s1.distr_list.push_back(-dms*J_s1.DoJack(F_disco_jack, 3, TKVK_barXX_OS_s1.col(0)[iens_XX][t], barXX_strange, data_TKVK_s1_OS_correlated.col(0)[iens_XX][t]));
	Jackknife J_s2(10000,Njacks);
	VKTK_sea_der_tm_s2.distr_list.push_back(-dms*J_s2.DoJack(F_disco_jack, 3, VKTK_barXX_tm_s2.col(0)[iens_XX][t], barXX_strange, data_VKTK_s2_tm_correlated.col(0)[iens_XX][t]));
	TKVK_sea_der_tm_s2.distr_list.push_back(-dms*J_s2.DoJack(F_disco_jack, 3, TKVK_barXX_tm_s2.col(0)[iens_XX][t], barXX_strange, data_TKVK_s2_tm_correlated.col(0)[iens_XX][t]));
	VKTK_sea_der_OS_s2.distr_list.push_back(-dms*J_s2.DoJack(F_disco_jack, 3, VKTK_barXX_OS_s2.col(0)[iens_XX][t], barXX_strange, data_VKTK_s2_OS_correlated.col(0)[iens_XX][t]));
	TKVK_sea_der_OS_s2.distr_list.push_back(-dms*J_s2.DoJack(F_disco_jack, 3, TKVK_barXX_OS_s2.col(0)[iens_XX][t], barXX_strange, data_TKVK_s2_OS_correlated.col(0)[iens_XX][t]));
    }

    
       
    //antysymmetrize
    for(int t=0; t < T_XX; t++) {
      VKTK_sea_der_tm_light.distr_list[t] = 0.5*(VKTK_sea_der_tm_light.distr_list[t] -  VKTK_sea_der_tm_light.distr_list[(T_XX -t+T_XX)%T_XX]);
      TKVK_sea_der_tm_light.distr_list[t] = 0.5*(TKVK_sea_der_tm_light.distr_list[t] -  TKVK_sea_der_tm_light.distr_list[(T_XX -t+T_XX)%T_XX]);
      VKTK_sea_der_OS_light.distr_list[t] = 0.5*(VKTK_sea_der_OS_light.distr_list[t] -  VKTK_sea_der_OS_light.distr_list[(T_XX -t+T_XX)%T_XX]);
      TKVK_sea_der_OS_light.distr_list[t] = 0.5*(TKVK_sea_der_OS_light.distr_list[t] -  TKVK_sea_der_OS_light.distr_list[(T_XX -t+T_XX)%T_XX]);

      VKTK_sea_der_tm_s1.distr_list[t] = 0.5*(VKTK_sea_der_tm_s1.distr_list[t] -  VKTK_sea_der_tm_s1.distr_list[(T_XX -t+T_XX)%T_XX]);
      TKVK_sea_der_tm_s1.distr_list[t] = 0.5*(TKVK_sea_der_tm_s1.distr_list[t] -  TKVK_sea_der_tm_s1.distr_list[(T_XX -t+T_XX)%T_XX]);
      VKTK_sea_der_OS_s1.distr_list[t] = 0.5*(VKTK_sea_der_OS_s1.distr_list[t] -  VKTK_sea_der_OS_s1.distr_list[(T_XX -t+T_XX)%T_XX]);
      TKVK_sea_der_OS_s1.distr_list[t] = 0.5*(TKVK_sea_der_OS_s1.distr_list[t] -  TKVK_sea_der_OS_s1.distr_list[(T_XX -t+T_XX)%T_XX]);

      VKTK_sea_der_tm_s2.distr_list[t] = 0.5*(VKTK_sea_der_tm_s2.distr_list[t] -  VKTK_sea_der_tm_s2.distr_list[(T_XX -t+T_XX)%T_XX]);
      TKVK_sea_der_tm_s2.distr_list[t] = 0.5*(TKVK_sea_der_tm_s2.distr_list[t] -  TKVK_sea_der_tm_s2.distr_list[(T_XX -t+T_XX)%T_XX]);
      VKTK_sea_der_OS_s2.distr_list[t] = 0.5*(VKTK_sea_der_OS_s2.distr_list[t] -  VKTK_sea_der_OS_s2.distr_list[(T_XX -t+T_XX)%T_XX]);
      TKVK_sea_der_OS_s2.distr_list[t] = 0.5*(TKVK_sea_der_OS_s2.distr_list[t] -  TKVK_sea_der_OS_s2.distr_list[(T_XX -t+T_XX)%T_XX]);
    }
    
    if(data_VKTK_tm.Tag[iens] == "cB211b.072.96") {
      distr_t_list FAKE_ZERO= 0.0*Get_id_jack_distr_list( T_XX/2 , Njacks) ;
      distr_t_list VKTK_ll_tm_tmp(UseJack), TKVK_ll_tm_tmp(UseJack), VKTK_ll_OS_tmp(UseJack), TKVK_ll_OS_tmp(UseJack); 
      distr_t_list VKTK_s1_tm_tmp(UseJack), TKVK_s1_tm_tmp(UseJack), VKTK_s1_OS_tmp(UseJack), TKVK_s1_OS_tmp(UseJack);
      distr_t_list VKTK_s2_tm_tmp(UseJack), TKVK_s2_tm_tmp(UseJack), VKTK_s2_OS_tmp(UseJack), TKVK_s2_OS_tmp(UseJack);

      for(int t=0; t < T_XX/2; t++) {
	VKTK_ll_tm_tmp.distr_list.push_back( VKTK_sea_der_tm_light.distr_list[t]);
	TKVK_ll_tm_tmp.distr_list.push_back( TKVK_sea_der_tm_light.distr_list[t]);
	VKTK_ll_OS_tmp.distr_list.push_back( VKTK_sea_der_OS_light.distr_list[t]);
	TKVK_ll_OS_tmp.distr_list.push_back( TKVK_sea_der_OS_light.distr_list[t]);
	VKTK_s1_tm_tmp.distr_list.push_back( VKTK_sea_der_tm_s1.distr_list[t]);
	TKVK_s1_tm_tmp.distr_list.push_back( TKVK_sea_der_tm_s1.distr_list[t]);
	VKTK_s1_OS_tmp.distr_list.push_back( VKTK_sea_der_OS_s1.distr_list[t]);
	TKVK_s1_OS_tmp.distr_list.push_back( TKVK_sea_der_OS_s1.distr_list[t]);
	VKTK_s2_tm_tmp.distr_list.push_back( VKTK_sea_der_tm_s2.distr_list[t]);
	TKVK_s2_tm_tmp.distr_list.push_back( TKVK_sea_der_tm_s2.distr_list[t]);
	VKTK_s2_OS_tmp.distr_list.push_back( VKTK_sea_der_OS_s2.distr_list[t]);
	TKVK_s2_OS_tmp.distr_list.push_back( TKVK_sea_der_OS_s2.distr_list[t]);
      }
      for(int t=0; t <= T_XX/4; t++) {
	VKTK_ll_tm_tmp.distr_list.push_back(  FAKE_ZERO.distr_list[t]);
	TKVK_ll_tm_tmp.distr_list.push_back(  FAKE_ZERO.distr_list[t]);
	VKTK_ll_OS_tmp.distr_list.push_back(  FAKE_ZERO.distr_list[t]);
	TKVK_ll_OS_tmp.distr_list.push_back(  FAKE_ZERO.distr_list[t]); 
	VKTK_s1_tm_tmp.distr_list.push_back(  FAKE_ZERO.distr_list[t]); 
	TKVK_s1_tm_tmp.distr_list.push_back(  FAKE_ZERO.distr_list[t]); 
	VKTK_s1_OS_tmp.distr_list.push_back(  FAKE_ZERO.distr_list[t]); 
	TKVK_s1_OS_tmp.distr_list.push_back(  FAKE_ZERO.distr_list[t]);
	VKTK_s2_tm_tmp.distr_list.push_back(  FAKE_ZERO.distr_list[t]); 
	TKVK_s2_tm_tmp.distr_list.push_back(  FAKE_ZERO.distr_list[t]); 
	VKTK_s2_OS_tmp.distr_list.push_back(  FAKE_ZERO.distr_list[t]); 
	TKVK_s2_OS_tmp.distr_list.push_back(  FAKE_ZERO.distr_list[t]); 
      }

      for(int t= Corr.Nt/2 + 1; t < Corr.Nt; t++) {
	VKTK_ll_tm_tmp.distr_list.push_back(  -1.0*VKTK_ll_tm_tmp.distr_list[ (Corr.Nt -t) ] );
	TKVK_ll_tm_tmp.distr_list.push_back(  -1.0*TKVK_ll_tm_tmp.distr_list[ (Corr.Nt -t) ] );
	VKTK_ll_OS_tmp.distr_list.push_back(  -1.0*VKTK_ll_OS_tmp.distr_list[ (Corr.Nt -t) ] );
	TKVK_ll_OS_tmp.distr_list.push_back(  -1.0*TKVK_ll_OS_tmp.distr_list[ (Corr.Nt -t) ] );
	VKTK_s1_tm_tmp.distr_list.push_back(  -1.0*VKTK_s1_tm_tmp.distr_list[ (Corr.Nt -t) ] );
	TKVK_s1_tm_tmp.distr_list.push_back(  -1.0*TKVK_s1_tm_tmp.distr_list[ (Corr.Nt -t) ] );
	VKTK_s1_OS_tmp.distr_list.push_back(  -1.0*VKTK_s1_OS_tmp.distr_list[ (Corr.Nt -t) ] );
	TKVK_s1_OS_tmp.distr_list.push_back(  -1.0*TKVK_s1_OS_tmp.distr_list[ (Corr.Nt -t) ] );
	VKTK_s2_tm_tmp.distr_list.push_back(  -1.0*VKTK_s2_tm_tmp.distr_list[ (Corr.Nt -t) ] );
	TKVK_s2_tm_tmp.distr_list.push_back(  -1.0*TKVK_s2_tm_tmp.distr_list[ (Corr.Nt -t) ] );
	VKTK_s2_OS_tmp.distr_list.push_back(  -1.0*VKTK_s2_OS_tmp.distr_list[ (Corr.Nt -t) ] );
	TKVK_s2_OS_tmp.distr_list.push_back(  -1.0*TKVK_s2_OS_tmp.distr_list[ (Corr.Nt -t) ] );
      }

      VKTK_sea_der_tm_light = VKTK_ll_tm_tmp;
      TKVK_sea_der_tm_light = TKVK_ll_tm_tmp;
      VKTK_sea_der_OS_light = VKTK_ll_OS_tmp;
      TKVK_sea_der_OS_light = TKVK_ll_OS_tmp;
      VKTK_sea_der_tm_s1 = VKTK_s1_tm_tmp;
      TKVK_sea_der_tm_s1 = TKVK_s1_tm_tmp;
      VKTK_sea_der_OS_s1 = VKTK_s1_OS_tmp;
      TKVK_sea_der_OS_s1 = TKVK_s1_OS_tmp;
      VKTK_sea_der_tm_s2 = VKTK_s2_tm_tmp;
      TKVK_sea_der_tm_s2 = TKVK_s2_tm_tmp;
      VKTK_sea_der_OS_s2 = VKTK_s2_OS_tmp;
      TKVK_sea_der_OS_s2 = TKVK_s2_OS_tmp;
    }
    //-------   -------//
    //-------------------------   FOR SEA QUARK MASS DERIVATIVE ---------------------------//
    cout<<"sea quark mass derivative computed..."<<endl;

    Print_To_File( {}, {VKTK_sea_der_tm_light.ave(), VKTK_sea_der_tm_light.err()}, "../data/magnetic_susc/"+data_TKTK_tm.Tag[iens]+"/ll_sea_der_VKTK_tm.dat.t", "", "");
    Print_To_File( {}, {TKVK_sea_der_tm_light.ave(), TKVK_sea_der_tm_light.err()}, "../data/magnetic_susc/"+data_TKTK_tm.Tag[iens]+"/ll_sea_der_TKVK_tm.dat.t", "", "");
    Print_To_File( {}, {VKTK_sea_der_OS_light.ave(), VKTK_sea_der_OS_light.err()}, "../data/magnetic_susc/"+data_TKTK_OS.Tag[iens]+"/ll_sea_der_VKTK_OS.dat.t", "", "");
    Print_To_File( {}, {TKVK_sea_der_OS_light.ave(), TKVK_sea_der_OS_light.err()}, "../data/magnetic_susc/"+data_TKTK_OS.Tag[iens]+"/ll_sea_der_TKVK_OS.dat.t", "", "");

    Print_To_File( {}, {VKTK_sea_der_tm_s1.ave(), VKTK_sea_der_tm_s1.err()}, "../data/magnetic_susc/"+data_TKTK_tm.Tag[iens]+"/s1_sea_der_VKTK_tm.dat.t", "", "");
    Print_To_File( {}, {TKVK_sea_der_tm_s1.ave(), TKVK_sea_der_tm_s1.err()}, "../data/magnetic_susc/"+data_TKTK_tm.Tag[iens]+"/s1_sea_der_TKVK_tm.dat.t", "", "");
    Print_To_File( {}, {VKTK_sea_der_OS_s1.ave(), VKTK_sea_der_OS_s1.err()}, "../data/magnetic_susc/"+data_TKTK_OS.Tag[iens]+"/s1_sea_der_VKTK_OS.dat.t", "", "");
    Print_To_File( {}, {TKVK_sea_der_OS_s1.ave(), TKVK_sea_der_OS_s1.err()}, "../data/magnetic_susc/"+data_TKTK_OS.Tag[iens]+"/s1_sea_der_TKVK_OS.dat.t", "", "");

    Print_To_File( {}, {VKTK_sea_der_tm_s2.ave(), VKTK_sea_der_tm_s2.err()}, "../data/magnetic_susc/"+data_TKTK_tm.Tag[iens]+"/s2_sea_der_VKTK_tm.dat.t", "", "");
    Print_To_File( {}, {TKVK_sea_der_tm_s2.ave(), TKVK_sea_der_tm_s2.err()}, "../data/magnetic_susc/"+data_TKTK_tm.Tag[iens]+"/s2_sea_der_TKVK_tm.dat.t", "", "");
    Print_To_File( {}, {VKTK_sea_der_OS_s2.ave(), VKTK_sea_der_OS_s2.err()}, "../data/magnetic_susc/"+data_TKTK_OS.Tag[iens]+"/s2_sea_der_VKTK_OS.dat.t", "", "");
    Print_To_File( {}, {TKVK_sea_der_OS_s2.ave(), TKVK_sea_der_OS_s2.err()}, "../data/magnetic_susc/"+data_TKTK_OS.Tag[iens]+"/s2_sea_der_TKVK_OS.dat.t", "", "");
      
      
    //tm
    distr_t_list Corr_TKTK_tm= Corr.corr_t(data_TKTK_tm.col(0)[iens], "../data/magnetic_susc/"+data_TKTK_tm.Tag[iens]+"/ll_TKTK_tm.dat");
    Corr.Reflection_sign=-1;
    
    distr_t_list Corr_VKTK_tm= Corr.corr_t(data_VKTK_tm.col(0)[iens], "../data/magnetic_susc/"+data_VKTK_tm.Tag[iens]+"/ll_VKTK_tm.dat") + Include_disco*(Zv/Za)*Corr.corr_t(data_disc_VKTK_light.col(0)[iens], "")/Qp;
    distr_t_list Corr_TKVK_tm= Corr.corr_t(data_TKVK_tm.col(0)[iens], "../data/magnetic_susc/"+data_TKVK_tm.Tag[iens]+"/ll_TKVK_tm.dat") + Include_disco*(Zv/Za)*Corr.corr_t(data_disc_VKTK_light.col(0)[iens], "")/Qp;
    distr_t_list Corr_VKTKVK_tm(UseJack);
    for(int t=0;t<Corr.Nt;t++) { Corr_VKTKVK_tm.distr_list.push_back( (t==0 || t==Corr.Nt/2)?Corr_VKTK_tm.distr_list[0]:((1.0/pow(Corr_VKTK_tm.err(t),2))*Corr_VKTK_tm.distr_list[t] +  (1.0/pow(Corr_TKVK_tm.err(t),2))*Corr_TKVK_tm.distr_list[t])/(       (1.0/pow(Corr_VKTK_tm.err(t),2)) + (1.0/pow(Corr_TKVK_tm.err(t),2)) ));}


    distr_t_list Corr_VKTK_tm_PLOT = Corr.corr_t( data_VKTK_tm.col(0)[iens], "");
    distr_t_list Corr_TKVK_tm_PLOT = Corr.corr_t( data_TKVK_tm.col(0)[iens], "");
    distr_t_list Corr_VKTKVK_tm_PLOT(UseJack);
    for(int t=0;t<Corr.Nt;t++) { Corr_VKTKVK_tm_PLOT.distr_list.push_back( (t==0 || t==Corr.Nt/2)?Corr_VKTK_tm_PLOT.distr_list[0]:((1.0/pow(Corr_VKTK_tm_PLOT.err(t),2))*Corr_VKTK_tm_PLOT.distr_list[t] +  (1.0/pow(Corr_TKVK_tm_PLOT.err(t),2))*Corr_TKVK_tm_PLOT.distr_list[t])/(       (1.0/pow(Corr_VKTK_tm_PLOT.err(t),2)) + (1.0/pow(Corr_TKVK_tm_PLOT.err(t),2)) ));}

    distr_t_list Corr_DISCO_VKTK_tm_PLOT = (Zv/Za)*Corr.corr_t(data_disc_VKTK_light.col(0)[iens], "")/Qp;
    

    Print_To_File({}, { (Corr_DISCO_VKTK_tm_PLOT/Corr_VKTKVK_tm_PLOT).ave(), (Corr_DISCO_VKTK_tm_PLOT/Corr_VKTKVK_tm_PLOT).err() , Corr_DISCO_VKTK_tm_PLOT.ave(), Corr_DISCO_VKTK_tm_PLOT.err()}, "../data/magnetic_susc/"+data_TKTK_OS.Tag[iens]+"/dC_DC_tm.dat", "", "");



    cout<<"TEST DISCO: "<<data_TKVK_tm.Tag[iens]<<endl;
    distr_t test_l_CONN(UseJack, UseJack?Njacks:800), test_l_DISC(UseJack,UseJack?Njacks:800);

    for(int t=1; t < Corr.Nt/2 -10; t++) {

      test_l_CONN = test_l_CONN - 1000*(2/a_distr)*Za*Z_T*Corr_VKTKVK_tm_PLOT.distr_list[t]*t;
      test_l_DISC = test_l_DISC - 1000*(2/a_distr)*Za*Z_T*Corr_DISCO_VKTK_tm_PLOT.distr_list[t]*t;

    }


    cout<<"CONN: "<<test_l_CONN.ave()<<" "<<test_l_CONN.err()<<endl;
    cout<<"DISC: "<<test_l_DISC.ave()<<" "<<test_l_DISC.err()<<endl;
    cout<<"################"<<endl;

  

   
     
    
    Corr.Reflection_sign=1;
    distr_t_list Corr_P5P5_tm= Corr.corr_t(data_P5P5_tm.col(0)[iens], "../data/magnetic_susc/"+data_P5P5_tm.Tag[iens]+"/ll_P5P5_tm.dat");
    //OS
    distr_t_list Corr_TKTK_OS= Corr.corr_t(data_TKTK_OS.col(0)[iens], "../data/magnetic_susc/"+data_TKTK_OS.Tag[iens]+"/ll_TKTK_OS.dat");
    
    Corr.Reflection_sign=-1;
    distr_t_list Corr_VKTK_OS= Corr.corr_t(data_VKTK_OS.col(0)[iens], "../data/magnetic_susc/"+data_VKTK_OS.Tag[iens]+"/ll_VKTK_OS.dat") + Include_disco*Corr.corr_t(data_disc_VKTK_light.col(0)[iens], "../data/magnetic_susc/"+data_TKTK_OS.Tag[iens]+"/disc_light_VKTK.dat")/Qp; 
    distr_t_list Corr_TKVK_OS= Corr.corr_t(data_TKVK_OS.col(0)[iens], "../data/magnetic_susc/"+data_TKVK_OS.Tag[iens]+"/ll_TKVK_OS.dat")  +  Include_disco*Corr.corr_t(data_disc_VKTK_light.col(0)[iens], "")/Qp; 
    distr_t_list Corr_VKTKVK_OS(UseJack);
    for(int t=0;t<Corr.Nt;t++) { Corr_VKTKVK_OS.distr_list.push_back( (t==0 || t==Corr.Nt/2)?Corr_VKTK_OS.distr_list[0]:((1.0/pow(Corr_VKTK_OS.err(t),2))*Corr_VKTK_OS.distr_list[t] +  (1.0/pow(Corr_TKVK_OS.err(t),2))*Corr_TKVK_OS.distr_list[t] )/(       (1.0/pow(Corr_VKTK_OS.err(t),2)) + (1.0/pow(Corr_TKVK_OS.err(t),2)) ));}
    //disco T1V1, T2V2, T3V3
    distr_t_list Corr_disco_V1T1= Corr.corr_t( data_disc_V1T1_light.col(0)[iens], "../data/magnetic_susc/"+data_TKTK_OS.Tag[iens]+"/disc_light_V1T1.dat");
    distr_t_list Corr_disco_V2T2= Corr.corr_t( data_disc_V2T2_light.col(0)[iens], "../data/magnetic_susc/"+data_TKTK_OS.Tag[iens]+"/disc_light_V2T2.dat");
    distr_t_list Corr_disco_V3T3= Corr.corr_t( data_disc_V3T3_light.col(0)[iens], "../data/magnetic_susc/"+data_TKTK_OS.Tag[iens]+"/disc_light_V3T3.dat");
    
    Corr.Reflection_sign=1;
    distr_t_list Corr_P5P5_OS= Corr.corr_t(data_P5P5_OS.col(0)[iens], "../data/magnetic_susc/"+data_P5P5_OS.Tag[iens]+"/ll_P5P5_OS.dat");
    //disco JJ
    distr_t_list Corr_disco_JJ= Corr.corr_t( data_disc_JJ_light.col(0)[iens], "../data/magnetic_susc/"+data_TKTK_OS.Tag[iens]+"/disc_JJ_light.dat");
   

    //###########################################################################
    //to estimate derivative
    Corr.Reflection_sign=1;
    distr_t_list Corr_der_TKVK_tm= -1*(Corr.corr_t(data_TKVK_tm_m2.col(0)[iens], "") - Corr.corr_t(data_TKVK_tm_m1.col(0)[iens], "")) + Include_sea_quark_mass_derivative*TKVK_sea_der_tm_light ;
    distr_t_list Corr_der_VKTK_tm= Corr.corr_t(data_VKTK_tm_m2.col(0)[iens], "") - Corr.corr_t(data_VKTK_tm_m1.col(0)[iens], "") + Include_sea_quark_mass_derivative*VKTK_sea_der_tm_light;
    distr_t_list Corr_der_VKTKVK_tm(UseJack);
    for(int t=0;t<Corr.Nt;t++) { Corr_der_VKTKVK_tm.distr_list.push_back( (t==0 || t==Corr.Nt/2)?Corr_der_VKTK_tm.distr_list[0]:((1.0/pow(Corr_der_VKTK_tm.err(t),2))*Corr_der_VKTK_tm.distr_list[t] +  (1.0/pow(Corr_der_TKVK_tm.err(t),2))*Corr_der_TKVK_tm.distr_list[t])/(       (1.0/pow(Corr_der_VKTK_tm.err(t),2)) + (1.0/pow(Corr_der_TKVK_tm.err(t),2)) ));}
    distr_t_list Corr_der_TKVK_OS= -1*(Corr.corr_t(data_TKVK_OS_m2.col(0)[iens], "") - Corr.corr_t(data_TKVK_OS_m1.col(0)[iens], "")) + Include_sea_quark_mass_derivative*TKVK_sea_der_OS_light;
    distr_t_list Corr_der_VKTK_OS= Corr.corr_t(data_VKTK_OS_m2.col(0)[iens], "") - Corr.corr_t(data_VKTK_OS_m1.col(0)[iens], "") + Include_sea_quark_mass_derivative*VKTK_sea_der_OS_light;
    distr_t_list Corr_der_VKTKVK_OS(UseJack);
    for(int t=0;t<Corr.Nt;t++) { Corr_der_VKTKVK_OS.distr_list.push_back( (t==0 || t==Corr.Nt/2)?Corr_der_VKTK_OS.distr_list[0]:((1.0/pow(Corr_der_VKTK_OS.err(t),2))*Corr_der_VKTK_OS.distr_list[t] +  (1.0/pow(Corr_der_TKVK_OS.err(t),2))*Corr_der_TKVK_OS.distr_list[t])/(       (1.0/pow(Corr_der_VKTK_OS.err(t),2)) + (1.0/pow(Corr_der_TKVK_OS.err(t),2)) ));}
    Corr.Reflection_sign=1;

    

    
    if(data_TKVK_tm.Tag[iens] != data_TKVK_tm_m1.Tag[iens]) crash("ensemble used to compute C and dC do not match");
   
    

    //strange 1
    //tm
    distr_t_list Corr_s1_TKTK_tm= Corr.corr_t(data_s1_TKTK_tm.col(0)[iens], "../data/magnetic_susc/"+data_TKTK_tm.Tag[iens]+"/s1_TKTK_tm.dat");
    
    Corr.Reflection_sign=-1;
    distr_t_list Corr_s1_VKTK_tm= Corr.corr_t(data_s1_VKTK_tm.col(0)[iens], "../data/magnetic_susc/"+data_VKTK_tm.Tag[iens]+"/s1_VKTK_tm.dat") + Include_disco*(Zv/Za)*Corr.corr_t(data_disc_VKTK_strange.col(0)[iens],  "")/Qn;
    distr_t_list Corr_s1_TKVK_tm= Corr.corr_t(data_s1_TKVK_tm.col(0)[iens], "../data/magnetic_susc/"+data_TKVK_tm.Tag[iens]+"/s1_TKVK_tm.dat") + Include_disco*(Zv/Za)*Corr.corr_t(data_disc_VKTK_strange.col(0)[iens],  "")/Qn;
    distr_t_list Corr_s1_VKTKVK_tm(UseJack);
    for(int t=0;t<Corr.Nt;t++) { Corr_s1_VKTKVK_tm.distr_list.push_back( (t==0 || t==Corr.Nt/2)?Corr_s1_VKTK_tm.distr_list[0]:((1.0/pow(Corr_s1_VKTK_tm.err(t),2))*Corr_s1_VKTK_tm.distr_list[t] +  (1.0/pow(Corr_s1_TKVK_tm.err(t),2))*Corr_s1_TKVK_tm.distr_list[t])/(       (1.0/pow(Corr_s1_VKTK_tm.err(t),2)) + (1.0/pow(Corr_s1_TKVK_tm.err(t),2)) ));}


    distr_t_list Corr_s1_VKTK_tm_PLOT = Corr.corr_t(data_s1_VKTK_tm.col(0)[iens], "");
    distr_t_list Corr_s1_TKVK_tm_PLOT = Corr.corr_t(data_s1_TKVK_tm.col(0)[iens], "");
    distr_t_list Corr_s1_VKTKVK_tm_PLOT(UseJack);
    for(int t=0;t<Corr.Nt;t++) { Corr_s1_VKTKVK_tm_PLOT.distr_list.push_back( (t==0 || t==Corr.Nt/2)?Corr_s1_VKTK_tm_PLOT.distr_list[0]:((1.0/pow(Corr_s1_VKTK_tm_PLOT.err(t),2))*Corr_s1_VKTK_tm_PLOT.distr_list[t] +  (1.0/pow(Corr_s1_TKVK_tm_PLOT.err(t),2))*Corr_s1_TKVK_tm_PLOT.distr_list[t])/(       (1.0/pow(Corr_s1_VKTK_tm_PLOT.err(t),2)) + (1.0/pow(Corr_s1_TKVK_tm_PLOT.err(t),2)) ));}

    distr_t_list Corr_DISCO_s1_VKTK_tm_PLOT= (Zv/Za)*Corr.corr_t(data_disc_VKTK_strange.col(0)[iens],"")/Qn;

    Print_To_File({}, { (Corr_DISCO_s1_VKTK_tm_PLOT/Corr_s1_VKTKVK_tm_PLOT).ave(), (Corr_DISCO_s1_VKTK_tm_PLOT/Corr_s1_VKTKVK_tm_PLOT).err() , Corr_DISCO_s1_VKTK_tm_PLOT.ave(), Corr_DISCO_s1_VKTK_tm_PLOT.err()}, "../data/magnetic_susc/"+data_TKTK_OS.Tag[iens]+"/dC_DC_strange_tm.dat", "", "");

    cout<<"TEST DISCO STRANGE: "<<data_TKVK_tm.Tag[iens]<<endl;
    distr_t test_s_CONN(UseJack, UseJack?Njacks:800), test_s_DISC(UseJack,UseJack?Njacks:800);

    for(int t=1; t < Corr.Nt/2 -10; t++) {

      test_s_CONN = test_s_CONN - 1000*(2/a_distr)*Za*Z_T*Corr_s1_VKTKVK_tm_PLOT.distr_list[t]*t;
      test_s_DISC = test_s_DISC - 1000*(2/a_distr)*Za*Z_T*Corr_DISCO_s1_VKTK_tm_PLOT.distr_list[t]*t;

    }


    cout<<"CONN: "<<test_s_CONN.ave()<<" "<<test_s_CONN.err()<<endl;
    cout<<"DISC: "<<test_s_DISC.ave()<<" "<<test_s_DISC.err()<<endl;
    cout<<"################"<<endl;
    
    
    
    Corr.Reflection_sign=1;
    distr_t_list Corr_s1_P5P5_tm= Corr.corr_t(data_s1_P5P5_tm.col(0)[iens], "../data/magnetic_susc/"+data_P5P5_tm.Tag[iens]+"/s1_P5P5_tm.dat");
    //OS
    distr_t_list Corr_s1_TKTK_OS= Corr.corr_t(data_s1_TKTK_OS.col(0)[iens], "../data/magnetic_susc/"+data_TKTK_OS.Tag[iens]+"/s1_TKTK_OS.dat");
    
    Corr.Reflection_sign=-1;
    distr_t_list Corr_s1_VKTK_OS= Corr.corr_t(data_s1_VKTK_OS.col(0)[iens], "../data/magnetic_susc/"+data_VKTK_OS.Tag[iens]+"/s1_VKTK_OS.dat")  + Include_disco*Corr.corr_t(data_disc_VKTK_strange.col(0)[iens], "../data/magnetic_susc/"+data_TKTK_OS.Tag[iens]+"/disc_stange_VKTK.dat")/Qn;
    distr_t_list Corr_s1_TKVK_OS= Corr.corr_t(data_s1_TKVK_OS.col(0)[iens], "../data/magnetic_susc/"+data_TKVK_OS.Tag[iens]+"/s1_TKVK_OS.dat") + Include_disco*Corr.corr_t(data_disc_VKTK_strange.col(0)[iens], "")/Qn;
    distr_t_list Corr_s1_VKTKVK_OS(UseJack);
    for(int t=0;t<Corr.Nt;t++) { Corr_s1_VKTKVK_OS.distr_list.push_back( (t==0 || t==Corr.Nt/2)?Corr_s1_VKTK_OS.distr_list[0]:((1.0/pow(Corr_s1_VKTK_OS.err(t),2))*Corr_s1_VKTK_OS.distr_list[t] +  (1.0/pow(Corr_s1_TKVK_OS.err(t),2))*Corr_s1_TKVK_OS.distr_list[t])/(       (1.0/pow(Corr_s1_VKTK_OS.err(t),2)) + (1.0/pow(Corr_s1_TKVK_OS.err(t),2))));}



    
    
    
    Corr.Reflection_sign=1;
    distr_t_list Corr_s1_P5P5_OS= Corr.corr_t(data_s1_P5P5_OS.col(0)[iens], "../data/magnetic_susc/"+data_P5P5_OS.Tag[iens]+"/s1_P5P5_OS.dat");


    //strange 2
    //tm
    distr_t_list Corr_s2_TKTK_tm= Corr.corr_t(data_s2_TKTK_tm.col(0)[iens], "../data/magnetic_susc/"+data_TKTK_tm.Tag[iens]+"/s2_TKTK_tm.dat");
    
    Corr.Reflection_sign=-1;
    distr_t_list Corr_s2_VKTK_tm= Corr.corr_t(data_s2_VKTK_tm.col(0)[iens], "../data/magnetic_susc/"+data_VKTK_tm.Tag[iens]+"/s2_VKTK_tm.dat") + Include_disco*(Zv/Za)*Corr.corr_t(data_disc_VKTK_strange.col(0)[iens], "")/Qn;
    distr_t_list Corr_s2_TKVK_tm= Corr.corr_t(data_s2_TKVK_tm.col(0)[iens], "../data/magnetic_susc/"+data_TKVK_tm.Tag[iens]+"/s2_TKVK_tm.dat") + Include_disco*(Zv/Za)*Corr.corr_t(data_disc_VKTK_strange.col(0)[iens], "")/Qn;
    distr_t_list Corr_s2_VKTKVK_tm(UseJack);
    for(int t=0;t<Corr.Nt;t++) { Corr_s2_VKTKVK_tm.distr_list.push_back( (t==0 || t==Corr.Nt/2)?Corr_s2_VKTK_tm.distr_list[0]:((1.0/pow(Corr_s2_VKTK_tm.err(t),2))*Corr_s2_VKTK_tm.distr_list[t] +  (1.0/pow(Corr_s2_TKVK_tm.err(t),2))*Corr_s2_TKVK_tm.distr_list[t])/(       (1.0/pow(Corr_s2_VKTK_tm.err(t),2)) + (1.0/pow(Corr_s2_TKVK_tm.err(t),2)) ));}

    
    Corr.Reflection_sign=1;
    distr_t_list Corr_s2_P5P5_tm= Corr.corr_t(data_s2_P5P5_tm.col(0)[iens], "../data/magnetic_susc/"+data_P5P5_tm.Tag[iens]+"/s2_P5P5_tm.dat");
    //OS
    distr_t_list Corr_s2_TKTK_OS= Corr.corr_t(data_s2_TKTK_OS.col(0)[iens], "../data/magnetic_susc/"+data_TKTK_OS.Tag[iens]+"/s2_TKTK_OS.dat");

    
    Corr.Reflection_sign=-1;
    distr_t_list Corr_s2_VKTK_OS= Corr.corr_t(data_s2_VKTK_OS.col(0)[iens], "../data/magnetic_susc/"+data_VKTK_OS.Tag[iens]+"/s2_VKTK_OS.dat") + Include_disco*Corr.corr_t(data_disc_VKTK_strange.col(0)[iens], "")/Qn;
    distr_t_list Corr_s2_TKVK_OS= Corr.corr_t(data_s2_TKVK_OS.col(0)[iens], "../data/magnetic_susc/"+data_TKVK_OS.Tag[iens]+"/s2_TKVK_OS.dat") + Include_disco*Corr.corr_t(data_disc_VKTK_strange.col(0)[iens], "")/Qn;
    distr_t_list Corr_s2_VKTKVK_OS(UseJack);
    for(int t=0;t<Corr.Nt;t++) { Corr_s2_VKTKVK_OS.distr_list.push_back( (t==0 || t==Corr.Nt/2)?Corr_s2_VKTK_OS.distr_list[0]:((1.0/pow(Corr_s2_VKTK_OS.err(t),2))*Corr_s2_VKTK_OS.distr_list[t] +  (1.0/pow(Corr_s2_TKVK_OS.err(t),2))*Corr_s2_TKVK_OS.distr_list[t])/(       (1.0/pow(Corr_s2_VKTK_OS.err(t),2)) + (1.0/pow(Corr_s2_TKVK_OS.err(t),2)) ));}

    
    Corr.Reflection_sign=1;
    distr_t_list Corr_s2_P5P5_OS= Corr.corr_t(data_s2_P5P5_OS.col(0)[iens], "../data/magnetic_susc/"+data_P5P5_OS.Tag[iens]+"/s2_P5P5_OS.dat");


    distr_t_list Corr_s1_der_TKVK_tm = Corr_s2_VKTK_tm -Corr_s1_VKTK_tm + Include_sea_quark_mass_derivative*TKVK_sea_der_tm_s1;
    distr_t_list Corr_s1_der_VKTK_tm = Corr_s2_TKVK_tm -Corr_s1_TKVK_tm + Include_sea_quark_mass_derivative*VKTK_sea_der_tm_s1;
    distr_t_list Corr_s1_der_VKTKVK_tm = Corr_s2_VKTKVK_tm - Corr_s1_VKTKVK_tm;
    distr_t_list Corr_s1_der_TKVK_OS = Corr_s2_VKTK_OS -Corr_s1_VKTK_OS + Include_sea_quark_mass_derivative*TKVK_sea_der_OS_s1;
    distr_t_list Corr_s1_der_VKTK_OS = Corr_s2_TKVK_OS -Corr_s1_TKVK_OS + Include_sea_quark_mass_derivative*VKTK_sea_der_OS_s1;
    distr_t_list Corr_s1_der_VKTKVK_OS = Corr_s2_VKTKVK_OS - Corr_s1_VKTKVK_OS;
    distr_t_list sea_quark_mass_der_s1_tm(UseJack, Corr.Nt,  UseJack?Njacks:800);
    distr_t_list sea_quark_mass_der_s1_OS(UseJack, Corr.Nt,  UseJack?Njacks:800);
    if(Include_sea_quark_mass_derivative ) {
      for(int t=0;t<Corr.Nt;t++) {
	sea_quark_mass_der_s1_tm.distr_list[t] = (( (t==0) || (t==Corr.Nt/2)  || ( VKTK_sea_der_tm_s1.ave(t)==0 && VKTK_sea_der_tm_s1.err(t) == 0 ) )?VKTK_sea_der_tm_s1.distr_list[0]:((1.0/pow(VKTK_sea_der_tm_s1.err(t),2))*VKTK_sea_der_tm_s1.distr_list[t] +  (1.0/pow(TKVK_sea_der_tm_s1.err(t),2))*TKVK_sea_der_tm_s1.distr_list[t])/(       (1.0/pow(VKTK_sea_der_tm_s1.err(t),2)) + (1.0/pow(TKVK_sea_der_tm_s1.err(t),2))));
	sea_quark_mass_der_s1_OS.distr_list[t] =  (( (t==0) || (t==Corr.Nt/2)  || ( VKTK_sea_der_OS_s1.ave(t)==0 && VKTK_sea_der_OS_s1.err(t) == 0 ) )?VKTK_sea_der_OS_s1.distr_list[0]:((1.0/pow(VKTK_sea_der_OS_s1.err(t),2))*VKTK_sea_der_OS_s1.distr_list[t] +  (1.0/pow(TKVK_sea_der_OS_s1.err(t),2))*TKVK_sea_der_OS_s1.distr_list[t])/(       (1.0/pow(VKTK_sea_der_OS_s1.err(t),2)) + (1.0/pow(TKVK_sea_der_OS_s1.err(t),2))));
	
	Corr_s1_der_VKTKVK_tm.distr_list[t] = Corr_s1_der_VKTKVK_tm.distr_list[t] + (( (t==0) || (t==Corr.Nt/2)  || ( VKTK_sea_der_tm_s1.ave(t)==0 && VKTK_sea_der_tm_s1.err(t) == 0 ) )?VKTK_sea_der_tm_s1.distr_list[0]:((1.0/pow(VKTK_sea_der_tm_s1.err(t),2))*VKTK_sea_der_tm_s1.distr_list[t] +  (1.0/pow(TKVK_sea_der_tm_s1.err(t),2))*TKVK_sea_der_tm_s1.distr_list[t])/(       (1.0/pow(VKTK_sea_der_tm_s1.err(t),2)) + (1.0/pow(TKVK_sea_der_tm_s1.err(t),2))));
	Corr_s1_der_VKTKVK_OS.distr_list[t] = Corr_s1_der_VKTKVK_OS.distr_list[t] + (( (t==0) || (t==Corr.Nt/2)  || ( VKTK_sea_der_OS_s1.ave(t)==0 && VKTK_sea_der_OS_s1.err(t) == 0 ) )?VKTK_sea_der_OS_s1.distr_list[0]:((1.0/pow(VKTK_sea_der_OS_s1.err(t),2))*VKTK_sea_der_OS_s1.distr_list[t] +  (1.0/pow(TKVK_sea_der_OS_s1.err(t),2))*TKVK_sea_der_OS_s1.distr_list[t])/(       (1.0/pow(VKTK_sea_der_OS_s1.err(t),2)) + (1.0/pow(TKVK_sea_der_OS_s1.err(t),2))));
      }
    }


    distr_t_list Corr_s2_der_TKVK_tm = Corr_s2_VKTK_tm -Corr_s1_VKTK_tm + Include_sea_quark_mass_derivative*TKVK_sea_der_tm_s2;
    distr_t_list Corr_s2_der_VKTK_tm = Corr_s2_TKVK_tm -Corr_s1_TKVK_tm + Include_sea_quark_mass_derivative*VKTK_sea_der_tm_s2;
    distr_t_list Corr_s2_der_VKTKVK_tm = Corr_s2_VKTKVK_tm - Corr_s1_VKTKVK_tm;
    distr_t_list Corr_s2_der_TKVK_OS = Corr_s2_VKTK_OS -Corr_s1_VKTK_OS + Include_sea_quark_mass_derivative*TKVK_sea_der_OS_s2;
    distr_t_list Corr_s2_der_VKTK_OS = Corr_s2_TKVK_OS -Corr_s1_TKVK_OS + Include_sea_quark_mass_derivative*VKTK_sea_der_OS_s2;
    distr_t_list Corr_s2_der_VKTKVK_OS = Corr_s2_VKTKVK_OS - Corr_s1_VKTKVK_OS;
    if(Include_sea_quark_mass_derivative ) {
      for(int t=0;t<Corr.Nt;t++) {

	
	Corr_s2_der_VKTKVK_tm.distr_list[t] = Corr_s2_der_VKTKVK_tm.distr_list[t] + (( (t==0) || (t==Corr.Nt/2)  || ( VKTK_sea_der_tm_s2.ave(t)==0 && VKTK_sea_der_tm_s2.err(t) == 0 ) )?VKTK_sea_der_tm_s2.distr_list[0]:((1.0/pow(VKTK_sea_der_tm_s2.err(t),2))*VKTK_sea_der_tm_s2.distr_list[t] +  (1.0/pow(TKVK_sea_der_tm_s2.err(t),2))*TKVK_sea_der_tm_s2.distr_list[t])/(       (1.0/pow(VKTK_sea_der_tm_s2.err(t),2)) + (1.0/pow(TKVK_sea_der_tm_s2.err(t),2))));
	Corr_s2_der_VKTKVK_OS.distr_list[t] = Corr_s2_der_VKTKVK_OS.distr_list[t] + (( (t==0) || (t==Corr.Nt/2)  || ( VKTK_sea_der_OS_s2.ave(t)==0 && VKTK_sea_der_OS_s2.err(t) == 0 ) )?VKTK_sea_der_OS_s2.distr_list[0]:((1.0/pow(VKTK_sea_der_OS_s2.err(t),2))*VKTK_sea_der_OS_s2.distr_list[t] +  (1.0/pow(TKVK_sea_der_OS_s2.err(t),2))*TKVK_sea_der_OS_s2.distr_list[t])/(       (1.0/pow(VKTK_sea_der_OS_s2.err(t),2)) + (1.0/pow(TKVK_sea_der_OS_s2.err(t),2))));
      }
    }

    distr_t syst_der(UseJack);
    for(int ijack=0;ijack<Njacks;ijack++) syst_der.distr.push_back( GM()*0.02/sqrt(Njacks-1.0));
    //add 2% systematic error or Corr_s1_der_VKTKVK_tm/OS and Corr_s2_der_VKTKVK_tm/OS
    Corr_s1_der_VKTKVK_tm = Corr_s1_der_VKTKVK_tm*( 1 + syst_der);
    Corr_s2_der_VKTKVK_tm = Corr_s2_der_VKTKVK_tm*( 1 + syst_der);
    Corr_s1_der_VKTKVK_OS = Corr_s1_der_VKTKVK_OS*( 1 + syst_der);
    Corr_s2_der_VKTKVK_OS = Corr_s2_der_VKTKVK_OS*( 1 + syst_der);
    

    //##########################################################
    //to estimate derivative approx.
    distr_t_list Corr_s_der_A_tm(UseJack), Corr_s_der_B_tm(UseJack), Corr_s_der_C_tm(UseJack), Corr_s_der_D_tm(UseJack), Corr_s_der_E_tm(UseJack), Corr_s_der_F_tm(UseJack);
    distr_t_list Corr_s_der_A_OS(UseJack), Corr_s_der_B_OS(UseJack), Corr_s_der_C_OS(UseJack), Corr_s_der_D_OS(UseJack), Corr_s_der_E_OS(UseJack), Corr_s_der_F_OS(UseJack);    
    
    if(data_P5P5_OS.Tag[iens] == "cB211b.072.64") {
      //strange 3
      //tm
      distr_t_list Corr_s3_TKTK_tm= Corr.corr_t(data_s3_TKTK_tm.col(0)[0], "../data/magnetic_susc/"+data_TKTK_tm.Tag[iens]+"/s3_TKTK_tm.dat");
    
      Corr.Reflection_sign=-1;
      distr_t_list Corr_s3_VKTK_tm= Corr.corr_t(data_s3_VKTK_tm.col(0)[0], "../data/magnetic_susc/"+data_VKTK_tm.Tag[iens]+"/s3_VKTK_tm.dat") + Include_disco*(Zv/Za)*Corr.corr_t(data_disc_VKTK_strange.col(0)[iens],  "")/Qn;
      distr_t_list Corr_s3_TKVK_tm= Corr.corr_t(data_s3_TKVK_tm.col(0)[0], "../data/magnetic_susc/"+data_TKVK_tm.Tag[iens]+"/s3_TKVK_tm.dat") + Include_disco*(Zv/Za)*Corr.corr_t(data_disc_VKTK_strange.col(0)[iens],  "")/Qn;
      distr_t_list Corr_s3_VKTKVK_tm(UseJack);
      for(int t=0;t<Corr.Nt;t++) { Corr_s3_VKTKVK_tm.distr_list.push_back( (t==0 || t==Corr.Nt/2)?Corr_s3_VKTK_tm.distr_list[0]:((1.0/pow(Corr_s3_VKTK_tm.err(t),2))*Corr_s3_VKTK_tm.distr_list[t] +  (1.0/pow(Corr_s3_TKVK_tm.err(t),2))*Corr_s3_TKVK_tm.distr_list[t])/(       (1.0/pow(Corr_s3_VKTK_tm.err(t),2)) + (1.0/pow(Corr_s3_TKVK_tm.err(t),2)) ));}
    
      Corr.Reflection_sign=1;
      distr_t_list Corr_s3_P5P5_tm= Corr.corr_t(data_s3_P5P5_tm.col(0)[0], "../data/magnetic_susc/"+data_P5P5_tm.Tag[iens]+"/s3_P5P5_tm.dat");
      //OS
      distr_t_list Corr_s3_TKTK_OS= Corr.corr_t(data_s3_TKTK_OS.col(0)[0], "../data/magnetic_susc/"+data_TKTK_OS.Tag[iens]+"/s3_TKTK_OS.dat");
    
      Corr.Reflection_sign=-1;
      distr_t_list Corr_s3_VKTK_OS= Corr.corr_t(data_s3_VKTK_OS.col(0)[0], "../data/magnetic_susc/"+data_VKTK_OS.Tag[iens]+"/s3_VKTK_OS.dat")  + Include_disco*Corr.corr_t(data_disc_VKTK_strange.col(0)[iens], "../data/magnetic_susc/"+data_TKTK_OS.Tag[iens]+"/disc_stange_VKTK.dat")/Qn;
      distr_t_list Corr_s3_TKVK_OS= Corr.corr_t(data_s3_TKVK_OS.col(0)[0], "../data/magnetic_susc/"+data_TKVK_OS.Tag[iens]+"/s3_TKVK_OS.dat") + Include_disco*Corr.corr_t(data_disc_VKTK_strange.col(0)[iens], "")/Qn;
      distr_t_list Corr_s3_VKTKVK_OS(UseJack);
      for(int t=0;t<Corr.Nt;t++) { Corr_s3_VKTKVK_OS.distr_list.push_back( (t==0 || t==Corr.Nt/2)?Corr_s3_VKTK_OS.distr_list[0]:((1.0/pow(Corr_s3_VKTK_OS.err(t),2))*Corr_s3_VKTK_OS.distr_list[t] +  (1.0/pow(Corr_s3_TKVK_OS.err(t),2))*Corr_s3_TKVK_OS.distr_list[t])/(       (1.0/pow(Corr_s3_VKTK_OS.err(t),2)) + (1.0/pow(Corr_s3_TKVK_OS.err(t),2))));}
    
      Corr.Reflection_sign=1;
      distr_t_list Corr_s3_P5P5_OS= Corr.corr_t(data_s3_P5P5_OS.col(0)[0], "../data/magnetic_susc/"+data_P5P5_OS.Tag[iens]+"/s3_P5P5_OS.dat");


      //strange 4
      //tm
      distr_t_list Corr_s4_TKTK_tm= Corr.corr_t(data_s4_TKTK_tm.col(0)[0], "../data/magnetic_susc/"+data_TKTK_tm.Tag[iens]+"/s4_TKTK_tm.dat");
    
      Corr.Reflection_sign=-1;
      distr_t_list Corr_s4_VKTK_tm= Corr.corr_t(data_s4_VKTK_tm.col(0)[0], "../data/magnetic_susc/"+data_VKTK_tm.Tag[iens]+"/s4_VKTK_tm.dat") + Include_disco*(Zv/Za)*Corr.corr_t(data_disc_VKTK_strange.col(0)[iens], "")/Qn;
      distr_t_list Corr_s4_TKVK_tm= Corr.corr_t(data_s4_TKVK_tm.col(0)[0], "../data/magnetic_susc/"+data_TKVK_tm.Tag[iens]+"/s4_TKVK_tm.dat") + Include_disco*(Zv/Za)*Corr.corr_t(data_disc_VKTK_strange.col(0)[iens], "")/Qn;
      distr_t_list Corr_s4_VKTKVK_tm(UseJack);
      for(int t=0;t<Corr.Nt;t++) { Corr_s4_VKTKVK_tm.distr_list.push_back( (t==0 || t==Corr.Nt/2)?Corr_s4_VKTK_tm.distr_list[0]:((1.0/pow(Corr_s4_VKTK_tm.err(t),2))*Corr_s4_VKTK_tm.distr_list[t] +  (1.0/pow(Corr_s4_TKVK_tm.err(t),2))*Corr_s4_TKVK_tm.distr_list[t])/(       (1.0/pow(Corr_s4_VKTK_tm.err(t),2)) + (1.0/pow(Corr_s4_TKVK_tm.err(t),2)) ));}

    
      Corr.Reflection_sign=1;
      distr_t_list Corr_s4_P5P5_tm= Corr.corr_t(data_s4_P5P5_tm.col(0)[0], "../data/magnetic_susc/"+data_P5P5_tm.Tag[iens]+"/s4_P5P5_tm.dat");
      //OS
      distr_t_list Corr_s4_TKTK_OS= Corr.corr_t(data_s4_TKTK_OS.col(0)[0], "../data/magnetic_susc/"+data_TKTK_OS.Tag[iens]+"/s4_TKTK_OS.dat");

    
      Corr.Reflection_sign=-1;
      distr_t_list Corr_s4_VKTK_OS= Corr.corr_t(data_s4_VKTK_OS.col(0)[0], "../data/magnetic_susc/"+data_VKTK_OS.Tag[iens]+"/s4_VKTK_OS.dat") + Include_disco*Corr.corr_t(data_disc_VKTK_strange.col(0)[iens], "")/Qn;
      distr_t_list Corr_s4_TKVK_OS= Corr.corr_t(data_s4_TKVK_OS.col(0)[0], "../data/magnetic_susc/"+data_TKVK_OS.Tag[iens]+"/s4_TKVK_OS.dat") + Include_disco*Corr.corr_t(data_disc_VKTK_strange.col(0)[iens], "")/Qn;
      distr_t_list Corr_s4_VKTKVK_OS(UseJack);
      for(int t=0;t<Corr.Nt;t++) { Corr_s4_VKTKVK_OS.distr_list.push_back( (t==0 || t==Corr.Nt/2)?Corr_s4_VKTK_OS.distr_list[0]:((1.0/pow(Corr_s4_VKTK_OS.err(t),2))*Corr_s4_VKTK_OS.distr_list[t] +  (1.0/pow(Corr_s4_TKVK_OS.err(t),2))*Corr_s4_TKVK_OS.distr_list[t])/(       (1.0/pow(Corr_s4_VKTK_OS.err(t),2)) + (1.0/pow(Corr_s4_TKVK_OS.err(t),2)) ));}

    
      Corr.Reflection_sign=1;
      distr_t_list Corr_s4_P5P5_OS= Corr.corr_t(data_s4_P5P5_OS.col(0)[0], "../data/magnetic_susc/"+data_P5P5_OS.Tag[iens]+"/s4_P5P5_OS.dat");


    


  

      Corr_s_der_A_tm = Za*Z_T*(Corr_s2_VKTKVK_tm - Corr_s1_VKTKVK_tm + sea_quark_mass_der_s1_tm);
      Corr_s_der_B_tm = Za*Z_T*( 2*(Corr_s1_VKTKVK_tm - Corr_s3_VKTKVK_tm) + sea_quark_mass_der_s1_tm);
      Corr_s_der_C_tm = Za*Z_T*( 2*(Corr_s4_VKTKVK_tm - Corr_s1_VKTKVK_tm) + sea_quark_mass_der_s1_tm);
      Corr_s_der_D_tm = Za*Z_T*( 2*(Corr_s2_VKTKVK_tm - Corr_s4_VKTKVK_tm) + sea_quark_mass_der_s1_tm);
      Corr_s_der_E_tm = Za*Z_T*( Corr_s4_VKTKVK_tm - Corr_s3_VKTKVK_tm + sea_quark_mass_der_s1_tm);
  
      Corr_s_der_A_OS = Zv*Z_T*(Corr_s2_VKTKVK_OS - Corr_s1_VKTKVK_OS + sea_quark_mass_der_s1_OS);
      Corr_s_der_B_OS = Zv*Z_T*( 2*(Corr_s1_VKTKVK_OS - Corr_s3_VKTKVK_OS) + sea_quark_mass_der_s1_OS);
      Corr_s_der_C_OS = Zv*Z_T*( 2*(Corr_s4_VKTKVK_OS - Corr_s1_VKTKVK_OS) + sea_quark_mass_der_s1_OS);
      Corr_s_der_D_OS = Zv*Z_T*( 2*(Corr_s2_VKTKVK_OS - Corr_s4_VKTKVK_OS) + sea_quark_mass_der_s1_OS);
      Corr_s_der_E_OS = Zv*Z_T*(Corr_s4_VKTKVK_OS - Corr_s3_VKTKVK_OS + sea_quark_mass_der_s1_OS);
      

    }

   
    
    //##########################################################
    

    //renormalize correlators
    //light
    Corr_TKTK_tm = Z_T*Z_T*Corr_TKTK_tm;
    Corr_VKTK_tm = Za*Z_T*Corr_VKTK_tm;
    Corr_TKVK_tm = Za*Z_T*Corr_TKVK_tm;
    Corr_VKTKVK_tm = Za*Z_T*Corr_VKTKVK_tm;
    Corr_TKTK_OS = Z_T*Z_T*Corr_TKTK_OS;
    Corr_VKTK_OS = Zv*Z_T*Corr_VKTK_OS;
    Corr_TKVK_OS = Zv*Z_T*Corr_TKVK_OS;
    Corr_VKTKVK_OS = Zv*Z_T*Corr_VKTKVK_OS;
    //to compute derivative
    Corr_der_VKTK_tm = Za*Z_T*Corr_der_VKTK_tm;
    Corr_der_TKVK_tm = Za*Z_T*Corr_der_TKVK_tm;
    Corr_der_VKTKVK_tm = Za*Z_T*Corr_der_VKTKVK_tm;
    Corr_der_VKTK_OS = Zv*Z_T*Corr_der_VKTK_OS;
    Corr_der_TKVK_OS = Zv*Z_T*Corr_der_TKVK_OS;
    Corr_der_VKTKVK_OS = Zv*Z_T*Corr_der_VKTKVK_OS;
    
    //strange 1
    Corr_s1_TKTK_tm = Z_T*Z_T*Corr_s1_TKTK_tm;
    Corr_s1_VKTK_tm = Za*Z_T*Corr_s1_VKTK_tm;
    Corr_s1_TKVK_tm = Za*Z_T*Corr_s1_TKVK_tm;
    Corr_s1_VKTKVK_tm = Za*Z_T*Corr_s1_VKTKVK_tm;
    Corr_s1_TKTK_OS = Z_T*Z_T*Corr_s1_TKTK_OS;
    Corr_s1_VKTK_OS = Zv*Z_T*Corr_s1_VKTK_OS;
    Corr_s1_TKVK_OS = Zv*Z_T*Corr_s1_TKVK_OS;
    Corr_s1_VKTKVK_OS = Zv*Z_T*Corr_s1_VKTKVK_OS;
    //strange 2
    Corr_s2_TKTK_tm = Z_T*Z_T*Corr_s2_TKTK_tm;
    Corr_s2_VKTK_tm = Za*Z_T*Corr_s2_VKTK_tm;
    Corr_s2_TKVK_tm = Za*Z_T*Corr_s2_TKVK_tm;
    Corr_s2_VKTKVK_tm = Za*Z_T*Corr_s2_VKTKVK_tm;
    Corr_s2_TKTK_OS = Z_T*Z_T*Corr_s2_TKTK_OS;
    Corr_s2_VKTK_OS = Zv*Z_T*Corr_s2_VKTK_OS;
    Corr_s2_TKVK_OS = Zv*Z_T*Corr_s2_TKVK_OS;
    Corr_s2_VKTKVK_OS = Zv*Z_T*Corr_s2_VKTKVK_OS;
    //to compute derivative
    //strange 1
    Corr_s1_der_VKTK_tm = Za*Z_T*Corr_s1_der_VKTK_tm;
    Corr_s1_der_TKVK_tm = Za*Z_T*Corr_s1_der_TKVK_tm;
    Corr_s1_der_VKTKVK_tm = Za*Z_T*Corr_s1_der_VKTKVK_tm;
    Corr_s1_der_VKTK_OS = Zv*Z_T*Corr_s1_der_VKTK_OS;
    Corr_s1_der_TKVK_OS = Zv*Z_T*Corr_s1_der_TKVK_OS;
    Corr_s1_der_VKTKVK_OS = Zv*Z_T*Corr_s1_der_VKTKVK_OS;
    //strange 2
    Corr_s2_der_VKTK_tm = Za*Z_T*Corr_s2_der_VKTK_tm;
    Corr_s2_der_TKVK_tm = Za*Z_T*Corr_s2_der_TKVK_tm;
    Corr_s2_der_VKTKVK_tm = Za*Z_T*Corr_s2_der_VKTKVK_tm;
    Corr_s2_der_VKTK_OS = Zv*Z_T*Corr_s2_der_VKTK_OS;
    Corr_s2_der_TKVK_OS = Zv*Z_T*Corr_s2_der_TKVK_OS;
    Corr_s2_der_VKTKVK_OS = Zv*Z_T*Corr_s2_der_VKTKVK_OS;


   


    //subtracted corr light
    distr_t_list Corr_sub_TV_tm  = Corr_TKVK_tm - (ml/dm)*Corr_der_TKVK_tm;
    distr_t_list Corr_sub_VT_tm  = Corr_VKTK_tm - (ml/dm)*Corr_der_VKTK_tm;
    distr_t_list Corr_sub_VTV_tm  = Corr_VKTKVK_tm - (ml/dm)*Corr_der_VKTKVK_tm;
    distr_t_list Corr_sub_TV_OS  = Corr_TKVK_OS - (ml/dm)*Corr_der_TKVK_OS;
    distr_t_list Corr_sub_VT_OS  = Corr_VKTK_OS - (ml/dm)*Corr_der_VKTK_OS;
    distr_t_list Corr_sub_VTV_OS  = Corr_VKTKVK_OS - (ml/dm)*Corr_der_VKTKVK_OS;

    //subtracted corr strange 1
    distr_t_list Corr_s1_sub_TV_tm  = Corr_s1_TKVK_tm - (ms1/dms)*Corr_s1_der_TKVK_tm;
    distr_t_list Corr_s1_sub_VT_tm  = Corr_s1_VKTK_tm - (ms1/dms)*Corr_s1_der_VKTK_tm;
    distr_t_list Corr_s1_sub_VTV_tm  = Corr_s1_VKTKVK_tm - (ms1/dms)*Corr_s1_der_VKTKVK_tm;
    distr_t_list Corr_s1_sub_TV_OS  = Corr_s1_TKVK_OS - (ms1/dms)*Corr_s1_der_TKVK_OS;
    distr_t_list Corr_s1_sub_VT_OS  = Corr_s1_VKTK_OS - (ms1/dms)*Corr_s1_der_VKTK_OS;
    distr_t_list Corr_s1_sub_VTV_OS  = Corr_s1_VKTKVK_OS - (ms1/dms)*Corr_s1_der_VKTKVK_OS;

    //subtracted corr strange 2
    distr_t_list Corr_s2_sub_TV_tm  = Corr_s2_TKVK_tm - (ms2/dms)*Corr_s2_der_TKVK_tm;
    distr_t_list Corr_s2_sub_VT_tm  = Corr_s2_VKTK_tm - (ms2/dms)*Corr_s2_der_VKTK_tm;
    distr_t_list Corr_s2_sub_VTV_tm  = Corr_s2_VKTKVK_tm - (ms2/dms)*Corr_s2_der_VKTKVK_tm;
    distr_t_list Corr_s2_sub_TV_OS  = Corr_s2_TKVK_OS - (ms2/dms)*Corr_s2_der_TKVK_OS;
    distr_t_list Corr_s2_sub_VT_OS  = Corr_s2_VKTK_OS - (ms2/dms)*Corr_s2_der_VKTK_OS;
    distr_t_list Corr_s2_sub_VTV_OS  = Corr_s2_VKTKVK_OS - (ms2/dms)*Corr_s2_der_VKTKVK_OS;

 
    
    //##### light
    distr_t susc_TV_tm(UseJack, UseJack?Njacks:800), susc_VT_tm(UseJack, UseJack?Njacks:800), susc_VTV_tm(UseJack, UseJack?Njacks:800);

    distr_t susc_TV_OS(UseJack, UseJack?Njacks:800), susc_VT_OS(UseJack, UseJack?Njacks:800), susc_VTV_OS(UseJack, UseJack?Njacks:800);
    //to evaluate derivative
    distr_t susc_der_TV_tm(UseJack, UseJack?Njacks:800), susc_der_VT_tm(UseJack, UseJack?Njacks:800),  susc_der_VTV_tm(UseJack, UseJack?Njacks:800);
    distr_t susc_der_TV_OS(UseJack, UseJack?Njacks:800), susc_der_VT_OS(UseJack, UseJack?Njacks:800),  susc_der_VTV_OS(UseJack, UseJack?Njacks:800);
    //integrated susc
    distr_t_list susc_TV_tm_int(UseJack, Corr.Nt/2 +1,  UseJack?Njacks:800), susc_VT_tm_int(UseJack, Corr.Nt/2+1,  UseJack?Njacks:800), susc_VTV_tm_int(UseJack, Corr.Nt/2+1,  UseJack?Njacks:800);
    distr_t_list susc_TV_OS_int(UseJack, Corr.Nt/2 +1,  UseJack?Njacks:800), susc_VT_OS_int(UseJack, Corr.Nt/2+1,  UseJack?Njacks:800), susc_VTV_OS_int(UseJack, Corr.Nt/2+1,  UseJack?Njacks:800);
 
    //integrated sub_susc
    distr_t_list susc_sub_TV_tm_int(UseJack, Corr.Nt/2 +1,  UseJack?Njacks:800), susc_sub_VT_tm_int(UseJack, Corr.Nt/2+1,  UseJack?Njacks:800), susc_sub_VTV_tm_int(UseJack, Corr.Nt/2+1,  UseJack?Njacks:800) ;
    distr_t_list susc_sub_TV_OS_int(UseJack, Corr.Nt/2 +1,  UseJack?Njacks:800), susc_sub_VT_OS_int(UseJack, Corr.Nt/2+1,  UseJack?Njacks:800), susc_sub_VTV_OS_int(UseJack, Corr.Nt/2+1,  UseJack?Njacks:800) ;
    //integrate derivative
    distr_t_list susc_der_TV_tm_int(UseJack, Corr.Nt/2 +1,  UseJack?Njacks:800), susc_der_VT_tm_int(UseJack, Corr.Nt/2+1,  UseJack?Njacks:800), susc_der_VTV_tm_int(UseJack, Corr.Nt/2+1,  UseJack?Njacks:800);
    distr_t_list susc_der_TV_OS_int(UseJack, Corr.Nt/2 +1,  UseJack?Njacks:800), susc_der_VT_OS_int(UseJack, Corr.Nt/2+1,  UseJack?Njacks:800), susc_der_VTV_OS_int(UseJack, Corr.Nt/2+1,  UseJack?Njacks:800);

    
    //###### strange 
    distr_t susc_s1_TV_tm(UseJack, UseJack?Njacks:800), susc_s1_VT_tm(UseJack, UseJack?Njacks:800), susc_s1_VTV_tm(UseJack, UseJack?Njacks:800);
    distr_t susc_s1_TV_OS(UseJack, UseJack?Njacks:800), susc_s1_VT_OS(UseJack, UseJack?Njacks:800), susc_s1_VTV_OS(UseJack, UseJack?Njacks:800);
    distr_t susc_s2_TV_tm(UseJack, UseJack?Njacks:800), susc_s2_VT_tm(UseJack, UseJack?Njacks:800), susc_s2_VTV_tm(UseJack, UseJack?Njacks:800);
    distr_t susc_s2_TV_OS(UseJack, UseJack?Njacks:800), susc_s2_VT_OS(UseJack, UseJack?Njacks:800), susc_s2_VTV_OS(UseJack, UseJack?Njacks:800);
    //to evaluate derivative 
    distr_t susc_s1_der_TV_tm(UseJack, UseJack?Njacks:800), susc_s1_der_VT_tm(UseJack, UseJack?Njacks:800), susc_s1_der_VTV_tm(UseJack, UseJack?Njacks:800);
    distr_t susc_s1_der_TV_OS(UseJack, UseJack?Njacks:800), susc_s1_der_VT_OS(UseJack, UseJack?Njacks:800), susc_s1_der_VTV_OS(UseJack, UseJack?Njacks:800);
    distr_t susc_s2_der_TV_tm(UseJack, UseJack?Njacks:800), susc_s2_der_VT_tm(UseJack, UseJack?Njacks:800), susc_s2_der_VTV_tm(UseJack, UseJack?Njacks:800);
    distr_t susc_s2_der_TV_OS(UseJack, UseJack?Njacks:800), susc_s2_der_VT_OS(UseJack, UseJack?Njacks:800), susc_s2_der_VTV_OS(UseJack, UseJack?Njacks:800);
    //integrated susc
    distr_t_list susc_s1_TV_tm_int(UseJack, Corr.Nt/2 +1,  UseJack?Njacks:800), susc_s1_VT_tm_int(UseJack, Corr.Nt/2+1,  UseJack?Njacks:800),  susc_s1_VTV_tm_int(UseJack, Corr.Nt/2+1,  UseJack?Njacks:800);
    distr_t_list susc_s1_TV_OS_int(UseJack, Corr.Nt/2 +1,  UseJack?Njacks:800), susc_s1_VT_OS_int(UseJack, Corr.Nt/2+1,  UseJack?Njacks:800),  susc_s1_VTV_OS_int(UseJack, Corr.Nt/2+1,  UseJack?Njacks:800);
    //integrated sub_susc
    distr_t_list susc_s1_sub_TV_tm_int(UseJack, Corr.Nt/2 +1,  UseJack?Njacks:800), susc_s1_sub_VT_tm_int(UseJack, Corr.Nt/2+1,  UseJack?Njacks:800),  susc_s1_sub_VTV_tm_int(UseJack, Corr.Nt/2+1,  UseJack?Njacks:800);
    distr_t_list susc_s1_sub_TV_OS_int(UseJack, Corr.Nt/2 +1,  UseJack?Njacks:800), susc_s1_sub_VT_OS_int(UseJack, Corr.Nt/2+1,  UseJack?Njacks:800),  susc_s1_sub_VTV_OS_int(UseJack, Corr.Nt/2+1,  UseJack?Njacks:800);
    //integrate derivative
    distr_t_list susc_s1_der_TV_tm_int(UseJack, Corr.Nt/2 + 1, UseJack?Njacks:800), susc_s1_der_VT_tm_int(UseJack, Corr.Nt/2 + 1, UseJack?Njacks:800), susc_s1_der_VTV_tm_int(UseJack, Corr.Nt/2 + 1, UseJack?Njacks:800);
    distr_t_list susc_s1_der_TV_OS_int(UseJack, Corr.Nt/2 + 1, UseJack?Njacks:800), susc_s1_der_VT_OS_int(UseJack, Corr.Nt/2 + 1, UseJack?Njacks:800), susc_s1_der_VTV_OS_int(UseJack, Corr.Nt/2 + 1, UseJack?Njacks:800);

    distr_t_list susc_s2_der_TV_tm_int(UseJack, Corr.Nt/2 + 1, UseJack?Njacks:800), susc_s2_der_VT_tm_int(UseJack, Corr.Nt/2 + 1, UseJack?Njacks:800), susc_s2_der_VTV_tm_int(UseJack, Corr.Nt/2 + 1, UseJack?Njacks:800);
    distr_t_list susc_s2_der_TV_OS_int(UseJack, Corr.Nt/2 + 1, UseJack?Njacks:800), susc_s2_der_VT_OS_int(UseJack, Corr.Nt/2 + 1, UseJack?Njacks:800), susc_s2_der_VTV_OS_int(UseJack, Corr.Nt/2 + 1, UseJack?Njacks:800);
      

    
 

    distr_t susc_s_der_A_tm(UseJack, UseJack?Njacks:800), susc_s_der_B_tm(UseJack, UseJack?Njacks:800), susc_s_der_C_tm(UseJack, UseJack?Njacks:800), susc_s_der_D_tm(UseJack, UseJack?Njacks:800), susc_s_der_E_tm(UseJack, UseJack?Njacks:800);
    distr_t susc_s_der_A_OS(UseJack, UseJack?Njacks:800), susc_s_der_B_OS(UseJack, UseJack?Njacks:800), susc_s_der_C_OS(UseJack, UseJack?Njacks:800), susc_s_der_D_OS(UseJack, UseJack?Njacks:800),  susc_s_der_E_OS(UseJack, UseJack?Njacks:800);

    //set time limit for integral ov derivative
    int Tmax_der_tm=Corr.Nt/2-10;
    int Tmax_der_OS=Corr.Nt/2-10;
    if(data_VKTK_tm.Tag[iens] == "cB211b.072.64") {    Tmax_der_tm=26;   Tmax_der_OS= 26;          }
    else if(data_VKTK_tm.Tag[iens] == "cB211b.072.96") {    Tmax_der_tm=26;   Tmax_der_OS= 44;          }
    else if( data_VKTK_tm.Tag[iens].substr(1,1) == "C") { Tmax_der_tm= 23; Tmax_der_OS=23;} //23;         
    else if( data_VKTK_tm.Tag[iens].substr(1,1) == "D") { Tmax_der_tm=24;  Tmax_der_OS=24;        }
    else crash("Ensemble: "+data_VKTK_tm.Tag[iens]+" not found");

   
    

    for(int t=1;t < Corr.Nt/2 - 10; t++) {

      

      //light
      susc_TV_tm = susc_TV_tm - 1000*(2/a_distr)*Corr_TKVK_tm.distr_list[t]*t;
      susc_VT_tm = susc_VT_tm - 1000*(2/a_distr)*Corr_VKTK_tm.distr_list[t]*t;
      susc_VTV_tm = susc_VTV_tm - 1000*(2/a_distr)*Corr_VKTKVK_tm.distr_list[t]*t;
     
      susc_TV_OS = susc_TV_OS - 1000*(2/a_distr)*Corr_TKVK_OS.distr_list[t]*t;
      susc_VT_OS = susc_VT_OS - 1000*(2/a_distr)*Corr_VKTK_OS.distr_list[t]*t;
      susc_VTV_OS = susc_VTV_OS - 1000*(2/a_distr)*Corr_VKTKVK_OS.distr_list[t]*t;
    

      //to evaluate derivative
      if(t<= Tmax_der_tm) {
	susc_der_TV_tm = susc_der_TV_tm - 1000*(2/a_distr)*Corr_der_TKVK_tm.distr_list[t]*t;
	susc_der_VT_tm = susc_der_VT_tm - 1000*(2/a_distr)*Corr_der_VKTK_tm.distr_list[t]*t;
	susc_der_VTV_tm = susc_der_VTV_tm - 1000*(2/a_distr)*Corr_der_VKTKVK_tm.distr_list[t]*t;
	
      }
      if(t<= Tmax_der_OS) {
	susc_der_TV_OS = susc_der_TV_OS - 1000*(2/a_distr)*Corr_der_TKVK_OS.distr_list[t]*t;
	susc_der_VT_OS = susc_der_VT_OS - 1000*(2/a_distr)*Corr_der_VKTK_OS.distr_list[t]*t;
	susc_der_VTV_OS = susc_der_VTV_OS - 1000*(2/a_distr)*Corr_der_VKTKVK_OS.distr_list[t]*t;

      }
      
      //integrated susc
      susc_TV_tm_int.distr_list[t] = susc_TV_tm_int.distr_list[t-1]  -1000*(2/a_distr)*Corr_TKVK_tm.distr_list[t]*t;
      susc_VT_tm_int.distr_list[t] = susc_VT_tm_int.distr_list[t-1]  -1000*(2/a_distr)*Corr_VKTK_tm.distr_list[t]*t;
      susc_VTV_tm_int.distr_list[t] = susc_VTV_tm_int.distr_list[t-1]  -1000*(2/a_distr)*Corr_VKTKVK_tm.distr_list[t]*t;
      susc_TV_OS_int.distr_list[t] = susc_TV_OS_int.distr_list[t-1]  -1000*(2/a_distr)*Corr_TKVK_OS.distr_list[t]*t;
      susc_VT_OS_int.distr_list[t] = susc_VT_OS_int.distr_list[t-1]  -1000*(2/a_distr)*Corr_VKTK_OS.distr_list[t]*t;
      susc_VTV_OS_int.distr_list[t] = susc_VTV_OS_int.distr_list[t-1]  -1000*(2/a_distr)*Corr_VKTKVK_OS.distr_list[t]*t;

      //integrated sub susc
      susc_sub_TV_tm_int.distr_list[t] = susc_sub_TV_tm_int.distr_list[t-1]  -1000*(2/a_distr)*((t<=Tmax_der_tm)?Corr_sub_TV_tm.distr_list[t]:Corr_TKVK_tm.distr_list[t])*t;
      susc_sub_VT_tm_int.distr_list[t] = susc_sub_VT_tm_int.distr_list[t-1]  -1000*(2/a_distr)*((t<=Tmax_der_tm)?Corr_sub_VT_tm.distr_list[t]:Corr_VKTK_tm.distr_list[t])*t;
      susc_sub_VTV_tm_int.distr_list[t] = susc_sub_VTV_tm_int.distr_list[t-1]  -1000*(2/a_distr)*((t<=Tmax_der_tm)?Corr_sub_VTV_tm.distr_list[t]:Corr_VKTKVK_tm.distr_list[t])*t;
      susc_sub_TV_OS_int.distr_list[t] = susc_sub_TV_OS_int.distr_list[t-1]  -1000*(2/a_distr)*((t<=Tmax_der_tm)?Corr_sub_TV_OS.distr_list[t]:Corr_TKVK_OS.distr_list[t])*t;
      susc_sub_VT_OS_int.distr_list[t] = susc_sub_VT_OS_int.distr_list[t-1]  -1000*(2/a_distr)*((t<=Tmax_der_tm)?Corr_sub_VT_OS.distr_list[t]:Corr_VKTK_OS.distr_list[t])*t;
      susc_sub_VTV_OS_int.distr_list[t] = susc_sub_VTV_OS_int.distr_list[t-1]  -1000*(2/a_distr)*((t<=Tmax_der_tm)?Corr_sub_VTV_OS.distr_list[t]:Corr_VKTKVK_OS.distr_list[t])*t;
    

      //integrated derivative
      susc_der_TV_tm_int.distr_list[t] = susc_der_TV_tm_int.distr_list[t-1] - 1000*(ml/dm)*(2/a_distr)*Corr_der_TKVK_tm.distr_list[t]*t;
      susc_der_VT_tm_int.distr_list[t] = susc_der_VT_tm_int.distr_list[t-1] - 1000*(ml/dm)*(2/a_distr)*Corr_der_VKTK_tm.distr_list[t]*t;
      susc_der_VTV_tm_int.distr_list[t] = susc_der_VTV_tm_int.distr_list[t-1] - 1000*(ml/dm)*(2/a_distr)*Corr_der_VKTKVK_tm.distr_list[t]*t;
      susc_der_TV_OS_int.distr_list[t] = susc_der_TV_OS_int.distr_list[t-1] - 1000*(ml/dm)*(2/a_distr)*Corr_der_TKVK_OS.distr_list[t]*t;
      susc_der_VT_OS_int.distr_list[t] = susc_der_VT_OS_int.distr_list[t-1] - 1000*(ml/dm)*(2/a_distr)*Corr_der_VKTK_OS.distr_list[t]*t;
      susc_der_VTV_OS_int.distr_list[t] = susc_der_VTV_OS_int.distr_list[t-1] - 1000*(ml/dm)*(2/a_distr)*Corr_der_VKTKVK_OS.distr_list[t]*t;
      

      

      //strange 1
      susc_s1_TV_tm = susc_s1_TV_tm - 1000*(2/a_distr)*Corr_s1_TKVK_tm.distr_list[t]*t;
      susc_s1_VT_tm = susc_s1_VT_tm - 1000*(2/a_distr)*Corr_s1_VKTK_tm.distr_list[t]*t;
      susc_s1_VTV_tm = susc_s1_VTV_tm - 1000*(2/a_distr)*Corr_s1_VKTKVK_tm.distr_list[t]*t;
   
      susc_s1_TV_OS = susc_s1_TV_OS - 1000*(2/a_distr)*Corr_s1_TKVK_OS.distr_list[t]*t;
      susc_s1_VT_OS = susc_s1_VT_OS - 1000*(2/a_distr)*Corr_s1_VKTK_OS.distr_list[t]*t;
      susc_s1_VTV_OS = susc_s1_VTV_OS - 1000*(2/a_distr)*Corr_s1_VKTKVK_OS.distr_list[t]*t;
     
      //integrated susc
      susc_s1_TV_tm_int.distr_list[t] = susc_s1_TV_tm_int.distr_list[t-1]  -1000*(2/a_distr)*Corr_s1_TKVK_tm.distr_list[t]*t;
      susc_s1_VT_tm_int.distr_list[t] = susc_s1_VT_tm_int.distr_list[t-1]  -1000*(2/a_distr)*Corr_s1_VKTK_tm.distr_list[t]*t;
      susc_s1_VTV_tm_int.distr_list[t] = susc_s1_VTV_tm_int.distr_list[t-1]  -1000*(2/a_distr)*Corr_s1_VKTKVK_tm.distr_list[t]*t;
      susc_s1_TV_OS_int.distr_list[t] = susc_s1_TV_OS_int.distr_list[t-1]  -1000*(2/a_distr)*Corr_s1_TKVK_OS.distr_list[t]*t;
      susc_s1_VT_OS_int.distr_list[t] = susc_s1_VT_OS_int.distr_list[t-1]  -1000*(2/a_distr)*Corr_s1_VKTK_OS.distr_list[t]*t;
      susc_s1_VTV_OS_int.distr_list[t] = susc_s1_VTV_OS_int.distr_list[t-1]  -1000*(2/a_distr)*Corr_s1_VKTKVK_OS.distr_list[t]*t;
      //integrated sub susc
      susc_s1_sub_TV_tm_int.distr_list[t] = susc_s1_sub_TV_tm_int.distr_list[t-1]  -1000*(2/a_distr)*Corr_s1_sub_TV_tm.distr_list[t]*t;
      susc_s1_sub_VT_tm_int.distr_list[t] = susc_s1_sub_VT_tm_int.distr_list[t-1]  -1000*(2/a_distr)*Corr_s1_sub_VT_tm.distr_list[t]*t;
      susc_s1_sub_VTV_tm_int.distr_list[t] = susc_s1_sub_VTV_tm_int.distr_list[t-1]  -1000*(2/a_distr)*Corr_s1_sub_VTV_tm.distr_list[t]*t;
      susc_s1_sub_TV_OS_int.distr_list[t] = susc_s1_sub_TV_OS_int.distr_list[t-1]  -1000*(2/a_distr)*Corr_s1_sub_TV_OS.distr_list[t]*t;
      susc_s1_sub_VT_OS_int.distr_list[t] = susc_s1_sub_VT_OS_int.distr_list[t-1]  -1000*(2/a_distr)*Corr_s1_sub_VT_OS.distr_list[t]*t;
      susc_s1_sub_VTV_OS_int.distr_list[t] = susc_s1_sub_VTV_OS_int.distr_list[t-1]  -1000*(2/a_distr)*Corr_s1_sub_VTV_OS.distr_list[t]*t;
      //integrated derivative
      susc_s1_der_TV_tm_int.distr_list[t] = susc_s1_der_TV_tm_int.distr_list[t-1] - 1000*(ms1/dms)*(2/a_distr)*Corr_s1_der_TKVK_tm.distr_list[t]*t;
      susc_s1_der_VT_tm_int.distr_list[t] = susc_s1_der_VT_tm_int.distr_list[t-1] - 1000*(ms1/dms)*(2/a_distr)*Corr_s1_der_VKTK_tm.distr_list[t]*t;
      susc_s1_der_VTV_tm_int.distr_list[t] = susc_s1_der_VTV_tm_int.distr_list[t-1] - 1000*(ms1/dms)*(2/a_distr)*Corr_s1_der_VKTKVK_tm.distr_list[t]*t;
      susc_s1_der_TV_OS_int.distr_list[t] = susc_s1_der_TV_OS_int.distr_list[t-1] - 1000*(ms1/dms)*(2/a_distr)*Corr_s1_der_TKVK_OS.distr_list[t]*t;
      susc_s1_der_VT_OS_int.distr_list[t] = susc_s1_der_VT_OS_int.distr_list[t-1] - 1000*(ms1/dms)*(2/a_distr)*Corr_s1_der_VKTK_OS.distr_list[t]*t;
      susc_s1_der_VTV_OS_int.distr_list[t] = susc_s1_der_VTV_OS_int.distr_list[t-1] - 1000*(ms1/dms)*(2/a_distr)*Corr_s1_der_VKTKVK_OS.distr_list[t]*t;
      

      //strange 2
      susc_s2_TV_tm = susc_s2_TV_tm - 1000*(2/a_distr)*Corr_s2_TKVK_tm.distr_list[t]*t;
      susc_s2_VT_tm = susc_s2_VT_tm - 1000*(2/a_distr)*Corr_s2_VKTK_tm.distr_list[t]*t;
      susc_s2_VTV_tm = susc_s2_VTV_tm - 1000*(2/a_distr)*Corr_s2_VKTKVK_tm.distr_list[t]*t;
      
      susc_s2_TV_OS = susc_s2_TV_OS - 1000*(2/a_distr)*Corr_s2_TKVK_OS.distr_list[t]*t;
      susc_s2_VT_OS = susc_s2_VT_OS - 1000*(2/a_distr)*Corr_s2_VKTK_OS.distr_list[t]*t;
      susc_s2_VTV_OS = susc_s2_VTV_OS - 1000*(2/a_distr)*Corr_s2_VKTKVK_OS.distr_list[t]*t;


      //integrated derivative
      susc_s2_der_TV_tm_int.distr_list[t] = susc_s2_der_TV_tm_int.distr_list[t-1] - 1000*(ms2/dms)*(2/a_distr)*Corr_s2_der_TKVK_tm.distr_list[t]*t;
      susc_s2_der_VT_tm_int.distr_list[t] = susc_s2_der_VT_tm_int.distr_list[t-1] - 1000*(ms2/dms)*(2/a_distr)*Corr_s2_der_VKTK_tm.distr_list[t]*t;
      susc_s2_der_VTV_tm_int.distr_list[t] = susc_s2_der_VTV_tm_int.distr_list[t-1] - 1000*(ms2/dms)*(2/a_distr)*Corr_s2_der_VKTKVK_tm.distr_list[t]*t;
      susc_s2_der_TV_OS_int.distr_list[t] = susc_s2_der_TV_OS_int.distr_list[t-1] - 1000*(ms2/dms)*(2/a_distr)*Corr_s2_der_TKVK_OS.distr_list[t]*t;
      susc_s2_der_VT_OS_int.distr_list[t] = susc_s2_der_VT_OS_int.distr_list[t-1] - 1000*(ms2/dms)*(2/a_distr)*Corr_s2_der_VKTK_OS.distr_list[t]*t;
      susc_s2_der_VTV_OS_int.distr_list[t] = susc_s2_der_VTV_OS_int.distr_list[t-1] - 1000*(ms2/dms)*(2/a_distr)*Corr_s2_der_VKTKVK_OS.distr_list[t]*t;



      
    
      //to evaluate derivative
      //strange 1
      susc_s1_der_TV_tm = susc_s1_der_TV_tm - 1000*(2/a_distr)*Corr_s1_der_TKVK_tm.distr_list[t]*t;
      susc_s1_der_VT_tm = susc_s1_der_VT_tm - 1000*(2/a_distr)*Corr_s1_der_VKTK_tm.distr_list[t]*t;
      susc_s1_der_VTV_tm = susc_s1_der_VTV_tm - 1000*(2/a_distr)*Corr_s1_der_VKTKVK_tm.distr_list[t]*t;
     
      susc_s1_der_TV_OS = susc_s1_der_TV_OS - 1000*(2/a_distr)*Corr_s1_der_TKVK_OS.distr_list[t]*t;
      susc_s1_der_VT_OS = susc_s1_der_VT_OS - 1000*(2/a_distr)*Corr_s1_der_VKTK_OS.distr_list[t]*t;
      susc_s1_der_VTV_OS = susc_s1_der_VTV_OS - 1000*(2/a_distr)*Corr_s1_der_VKTKVK_OS.distr_list[t]*t;

      //strange 2
      susc_s2_der_TV_tm = susc_s2_der_TV_tm - 1000*(2/a_distr)*Corr_s2_der_TKVK_tm.distr_list[t]*t;
      susc_s2_der_VT_tm = susc_s2_der_VT_tm - 1000*(2/a_distr)*Corr_s2_der_VKTK_tm.distr_list[t]*t;
      susc_s2_der_VTV_tm = susc_s2_der_VTV_tm - 1000*(2/a_distr)*Corr_s2_der_VKTKVK_tm.distr_list[t]*t;
     
      susc_s2_der_TV_OS = susc_s2_der_TV_OS - 1000*(2/a_distr)*Corr_s2_der_TKVK_OS.distr_list[t]*t;
      susc_s2_der_VT_OS = susc_s2_der_VT_OS - 1000*(2/a_distr)*Corr_s2_der_VKTK_OS.distr_list[t]*t;
      susc_s2_der_VTV_OS = susc_s2_der_VTV_OS - 1000*(2/a_distr)*Corr_s2_der_VKTKVK_OS.distr_list[t]*t;


      
      if(data_VKTK_OS.Tag[iens] == "cB211b.072.64") {

	susc_s_der_A_tm = susc_s_der_A_tm - (ms1/dms)*1000*(2/a_distr)*Corr_s_der_A_tm.distr_list[t]*t;
	susc_s_der_B_tm = susc_s_der_B_tm - (ms1/dms)*1000*(2/a_distr)*Corr_s_der_B_tm.distr_list[t]*t;
	susc_s_der_C_tm = susc_s_der_C_tm - (ms1/dms)*1000*(2/a_distr)*Corr_s_der_C_tm.distr_list[t]*t;
	susc_s_der_D_tm = susc_s_der_D_tm - (ms1/dms)*1000*(2/a_distr)*Corr_s_der_D_tm.distr_list[t]*t;
	susc_s_der_E_tm = susc_s_der_E_tm - (ms1/dms)*1000*(2/a_distr)*Corr_s_der_E_tm.distr_list[t]*t;

	susc_s_der_A_OS = susc_s_der_A_OS - (ms1/dms)*1000*(2/a_distr)*Corr_s_der_A_OS.distr_list[t]*t;
	susc_s_der_B_OS = susc_s_der_B_OS - (ms1/dms)*1000*(2/a_distr)*Corr_s_der_B_OS.distr_list[t]*t;
	susc_s_der_C_OS = susc_s_der_C_OS - (ms1/dms)*1000*(2/a_distr)*Corr_s_der_C_OS.distr_list[t]*t;
	susc_s_der_D_OS = susc_s_der_D_OS - (ms1/dms)*1000*(2/a_distr)*Corr_s_der_D_OS.distr_list[t]*t;
	susc_s_der_E_OS = susc_s_der_E_OS - (ms1/dms)*1000*(2/a_distr)*Corr_s_der_E_OS.distr_list[t]*t;

      
      }
      
    }

    
        

    //perform substraction of three-level log divergence
    
    
    distr_t tree_level_log_naive_ml= 1000*(3/(2*M_PI*M_PI*a_distr))*ml*log(ml);
    distr_t tree_level_log_naive_ms1= 1000*(3/(2*M_PI*M_PI*a_distr))*ms1*log(ms1);
    distr_t tree_level_log_naive_ms2= 1000*(3/(2*M_PI*M_PI*a_distr))*ms2*log(ms2);

    distr_t tree_level_log_ml = -1000*(ml/a_distr)*Get_tree_lev_der(2*ml);
    distr_t tree_level_log_ms1 = -1000*(ms1/a_distr)*Get_tree_lev_der(2*ms1);
    distr_t tree_level_log_ms2 = -1000*(ms2/a_distr)*Get_tree_lev_der(2*ms2);

    if(Use_tree_level_sub) {
      susc_TV_tm = susc_TV_tm + tree_level_log_ml;
      susc_VT_tm = susc_VT_tm + tree_level_log_ml;
      susc_VTV_tm = susc_VTV_tm + tree_level_log_ml;
      susc_TV_OS = susc_TV_OS + tree_level_log_ml;
      susc_VT_OS = susc_VT_OS + tree_level_log_ml;
      susc_VTV_OS = susc_VTV_OS + tree_level_log_ml;

      susc_s1_TV_tm = susc_s1_TV_tm + tree_level_log_ms1;
      susc_s1_VT_tm = susc_s1_VT_tm + tree_level_log_ms1;
      susc_s1_VTV_tm = susc_s1_VTV_tm + tree_level_log_ms1;
      susc_s1_TV_OS = susc_s1_TV_OS + tree_level_log_ms1;
      susc_s1_VT_OS = susc_s1_VT_OS + tree_level_log_ms1;
      susc_s1_VTV_OS = susc_s1_VTV_OS + tree_level_log_ms1;

      susc_s2_TV_tm = susc_s2_TV_tm + tree_level_log_ms2;
      susc_s2_VT_tm = susc_s2_VT_tm + tree_level_log_ms2;
      susc_s2_VTV_tm = susc_s2_VTV_tm + tree_level_log_ms2;
      susc_s2_TV_OS = susc_s2_TV_OS + tree_level_log_ms2;
      susc_s2_VT_OS = susc_s2_VT_OS + tree_level_log_ms2;
      susc_s2_VTV_OS = susc_s2_VTV_OS + tree_level_log_ms2;
    }
    else {

      double dlog_ml= 1 + 1/(log(ml));
      double dlog_ms1= 1 + 1/(log(ms1));
      double dlog_ms2= 1 + 1/(log(ms2));

      

      distr_t coeff_light_TV_tm = (a_distr/1000)*susc_der_TV_tm/(log(ml)*dm);
      distr_t coeff_light_VT_tm = (a_distr/1000)*susc_der_VT_tm/(log(ml)*dm);
      distr_t coeff_light_VTV_tm = (a_distr/1000)*susc_der_VTV_tm/(log(ml)*dm);
      distr_t coeff_light_TV_OS = (a_distr/1000)*susc_der_TV_OS/(log(ml)*dm);
      distr_t coeff_light_VT_OS = (a_distr/1000)*susc_der_VT_OS/(log(ml)*dm);
      distr_t coeff_light_VTV_OS = (a_distr/1000)*susc_der_VTV_OS/(log(ml)*dm);

     
      

      distr_t coeff_strange_TV_tm = (a_distr/1000)*susc_s1_der_TV_tm/(log(ms1)*dms);
      distr_t coeff_strange_VT_tm = (a_distr/1000)*susc_s1_der_VT_tm/(log(ms1)*dms);
      distr_t coeff_strange_VTV_tm = (a_distr/1000)*susc_s1_der_VTV_tm/(log(ms1)*dms);
      distr_t coeff_strange_TV_OS = (a_distr/1000)*susc_s1_der_TV_OS/(log(ms1)*dms);
      distr_t coeff_strange_VT_OS = (a_distr/1000)*susc_s1_der_VT_OS/(log(ms1)*dms);
      distr_t coeff_strange_VTV_OS = (a_distr/1000)*susc_s1_der_VTV_OS/(log(ms1)*dms);

      susc_TV_tm = susc_TV_tm + (-1000/a_distr)*coeff_light_TV_tm*ml*log(ml);
      susc_VT_tm = susc_VT_tm + (-1000/a_distr)*coeff_light_VT_tm*ml*log(ml);
      susc_VTV_tm = susc_VTV_tm + (-1000/a_distr)*coeff_light_VTV_tm*ml*log(ml);   
      susc_TV_OS = susc_TV_OS + (-1000/a_distr)*coeff_light_TV_OS*ml*log(ml);
      susc_VT_OS = susc_VT_OS + (-1000/a_distr)*coeff_light_VT_OS*ml*log(ml);
      susc_VTV_OS = susc_VTV_OS + (-1000/a_distr)*coeff_light_VTV_OS*ml*log(ml);
      
      susc_s1_TV_tm = susc_s1_TV_tm + (-ms1/dms)*susc_s1_der_TV_tm;
      susc_s1_VT_tm = susc_s1_VT_tm + (-ms1/dms)*susc_s1_der_VT_tm;
      susc_s1_VTV_tm = susc_s1_VTV_tm + (-ms1/dms)*susc_s1_der_VTV_tm;
      susc_s1_TV_OS = susc_s1_TV_OS + (-ms1/dms)*susc_s1_der_TV_OS;
      susc_s1_VT_OS = susc_s1_VT_OS + (-ms1/dms)*susc_s1_der_VT_OS;
      susc_s1_VTV_OS = susc_s1_VTV_OS + (-ms1/dms)*susc_s1_der_VTV_OS;

      susc_s2_TV_tm = susc_s2_TV_tm + (-ms2/dms)*susc_s2_der_TV_tm;
      susc_s2_VT_tm = susc_s2_VT_tm + (-ms2/dms)*susc_s2_der_VT_tm;
      susc_s2_VTV_tm = susc_s2_VTV_tm + (-ms2/dms)*susc_s2_der_VTV_tm;
      susc_s2_TV_OS = susc_s2_TV_OS + (-ms2/dms)*susc_s2_der_TV_OS;
      susc_s2_VT_OS = susc_s2_VT_OS + (-ms2/dms)*susc_s2_der_VT_OS;
      susc_s2_VTV_OS = susc_s2_VTV_OS + (-ms2/dms)*susc_s2_der_VTV_OS;
      

      cout<<"Ensemble: "<<data_TKVK_tm.Tag[iens]<<endl;
      cout<<"##### light #####"<<endl;
      cout<<"SUSC(l,tm): "<< susc_VTV_tm.ave()<<" +- "<<susc_VTV_tm.err()<<endl;
      cout<<"SUSC(l,OS): "<< susc_VTV_OS.ave()<<" +- "<<susc_VTV_OS.err()<<endl; 
      cout<<"der_TV(l,tm): "<<(-ml/dm)*susc_der_TV_tm.ave()<<" +- "<<(-ml/dm)*susc_der_TV_tm.err()<<" tree: "<<tree_level_log_ml.ave()<<" tree naive: "<<tree_level_log_naive_ml.ave()<<endl;
      cout<<"der_VT(l,tm): "<<(-ml/dm)*susc_der_VT_tm.ave()<<" +- "<<(-ml/dm)*susc_der_VT_tm.err()<<" tree: "<<tree_level_log_ml.ave()<<" tree naive: "<<tree_level_log_naive_ml.ave()<<endl;
      cout<<"der_VTV(l,tm): "<<(-ml/dm)*susc_der_VTV_tm.ave()<<" +- "<<(-ml/dm)*susc_der_VTV_tm.err()<<" tree: "<<tree_level_log_ml.ave()<<" tree naive: "<<tree_level_log_naive_ml.ave()<<endl;
      cout<<"der_TV(l,OS): "<<(-ml/dm)*susc_der_TV_OS.ave()<<" +- "<<(-ml/dm)*susc_der_TV_OS.err()<<" tree: "<<tree_level_log_ml.ave()<<" tree naive: "<<tree_level_log_naive_ml.ave()<<endl;
      cout<<"der_VT(l,OS): "<<(-ml/dm)*susc_der_VT_OS.ave()<<" +- "<<(-ml/dm)*susc_der_VT_OS.err()<<" tree: "<<tree_level_log_ml.ave()<<" tree naive: "<<tree_level_log_naive_ml.ave()<<endl;
      cout<<"der_VTV(l,OS): "<<(-ml/dm)*susc_der_VTV_OS.ave()<<" +- "<<(-ml/dm)*susc_der_VTV_OS.err()<<" tree: "<<tree_level_log_ml.ave()<<" tree naive: "<<tree_level_log_naive_ml.ave()<<endl;

      
      cout<<"###### strange 1 #####"<<endl;
      cout<<"SUSC(s,tm): "<< susc_s1_VTV_tm.ave()<<" +- "<<susc_s1_VTV_tm.err()<<endl;
      cout<<"SUSC(s,OS): "<< susc_s1_VTV_OS.ave()<<" +- "<<susc_s1_VTV_OS.err()<<endl;      
      cout<<"der_TV(s1,tm): "<<(-ms1/dms)*susc_s1_der_TV_tm.ave()<<" +- "<<(-ms1/dms)*susc_s1_der_TV_tm.err()<<" tree: "<<tree_level_log_ms1.ave()<<" tree naive: "<<tree_level_log_naive_ms1.ave()<<endl;
      cout<<"der_VT(s1,tm): "<<(-ms1/dms)*susc_s1_der_VT_tm.ave()<<" +- "<<(-ms1/dms)*susc_s1_der_VT_tm.err()<<" tree: "<<tree_level_log_ms1.ave()<<" tree naive: "<<tree_level_log_naive_ms1.ave()<<endl;
      cout<<"der_VTV(s1,tm): "<<(-ms1/dms)*susc_s1_der_VTV_tm.ave()<<" +- "<<(-ms1/dms)*susc_s1_der_VTV_tm.err()<<" tree: "<<tree_level_log_ms1.ave()<<" tree naive: "<<tree_level_log_naive_ms1.ave()<<endl;
      cout<<"der_TV(s1,OS): "<<(-ms1/dms)*susc_s1_der_TV_OS.ave()<<" +- "<<(-ms1/dms)*susc_s1_der_TV_OS.err()<<" tree: "<<tree_level_log_ms1.ave()<<" tree naive: "<<tree_level_log_naive_ms1.ave()<<endl;
      cout<<"der_VT(s1,OS): "<<(-ms1/dms)*susc_s1_der_VT_OS.ave()<<" +- "<<(-ms1/dms)*susc_s1_der_VT_OS.err()<<" tree: "<<tree_level_log_ms1.ave()<<" tree naive: "<<tree_level_log_naive_ms1.ave()<<endl;
      cout<<"der_VTV(s1,OS): "<<(-ms1/dms)*susc_s1_der_VTV_OS.ave()<<" +- "<<(-ms1/dms)*susc_s1_der_VTV_OS.err()<<" tree: "<<tree_level_log_ms1.ave()<<" tree naive: "<<tree_level_log_naive_ms1.ave()<<endl;
      
      cout<<"###### strange 2 #####"<<endl;      
      cout<<"SUSC(s,tm): "<< susc_s2_VTV_tm.ave()<<" +- "<<susc_s2_VTV_tm.err()<<endl;
      cout<<"SUSC(s,OS): "<< susc_s2_VTV_OS.ave()<<" +- "<<susc_s2_VTV_OS.err()<<endl;       
      cout<<"der_TV(s2,tm): "<<(-ms2/dms)*susc_s2_der_TV_tm.ave()<<" +- "<<(-ms2/dms)*susc_s2_der_TV_tm.err()<<" tree: "<<tree_level_log_ms2.ave()<<" tree naive: "<<tree_level_log_naive_ms2.ave()<<endl;
      cout<<"der_VT(s2,tm): "<<(-ms2/dms)*susc_s2_der_VT_tm.ave()<<" +- "<<(-ms2/dms)*susc_s2_der_VT_tm.err()<<" tree: "<<tree_level_log_ms2.ave()<<" tree naive: "<<tree_level_log_naive_ms2.ave()<<endl;
      cout<<"der_VTV(s2,tm): "<<(-ms2/dms)*susc_s2_der_VTV_tm.ave()<<" +- "<<(-ms2/dms)*susc_s2_der_VTV_tm.err()<<" tree: "<<tree_level_log_ms2.ave()<<" tree naive: "<<tree_level_log_naive_ms2.ave()<<endl;
      cout<<"der_TV(s2,OS): "<<(-ms2/dms)*susc_s2_der_TV_OS.ave()<<" +- "<<(-ms2/dms)*susc_s2_der_TV_OS.err()<<" tree: "<<tree_level_log_ms2.ave()<<" tree naive: "<<tree_level_log_naive_ms2.ave()<<endl;
      cout<<"der_VT(s2,OS): "<<(-ms2/dms)*susc_s2_der_VT_OS.ave()<<" +- "<<(-ms2/dms)*susc_s2_der_VT_OS.err()<<" tree: "<<tree_level_log_ms2.ave()<<" tree naive: "<<tree_level_log_naive_ms2.ave()<<endl;
      cout<<"der_VTV(s2,OS): "<<(-ms2/dms)*susc_s2_der_VTV_OS.ave()<<" +- "<<(-ms2/dms)*susc_s2_der_VTV_OS.err()<<" tree: "<<tree_level_log_ms2.ave()<<" tree naive: "<<tree_level_log_naive_ms2.ave()<<endl;
      cout<<"##### coeffs #####"<<endl;
      cout<<"TV(tm):  (l) "<<coeff_light_TV_tm.ave()<<" +- "<<coeff_light_TV_tm.err()<<" (s) "<<coeff_strange_TV_tm.ave()<<" +- "<<coeff_strange_TV_tm.err()<<endl;
      cout<<"VT(tm):  (l) "<<coeff_light_VT_tm.ave()<<" +- "<<coeff_light_VT_tm.err()<<" (s) "<<coeff_strange_VT_tm.ave()<<" +- "<<coeff_strange_VT_tm.err()<<endl;
      cout<<"VTV(tm):  (l) "<<coeff_light_VTV_tm.ave()<<" +- "<<coeff_light_VTV_tm.err()<<" (s) "<<coeff_strange_VTV_tm.ave()<<" +- "<<coeff_strange_VTV_tm.err()<<endl;
      cout<<"TV(OS):  (l) "<<coeff_light_TV_OS.ave()<<" +- "<<coeff_light_TV_OS.err()<<" (s) "<<coeff_strange_TV_OS.ave()<<" +- "<<coeff_strange_TV_OS.err()<<endl;
      cout<<"VT(OS):  (l) "<<coeff_light_VT_OS.ave()<<" +- "<<coeff_light_VT_OS.err()<<" (s) "<<coeff_strange_VT_OS.ave()<<" +- "<<coeff_strange_VT_OS.err()<<endl;
      cout<<"VTV(OS):  (l) "<<coeff_light_VTV_OS.ave()<<" +- "<<coeff_light_VTV_OS.err()<<" (s) "<<coeff_strange_VTV_OS.ave()<<" +- "<<coeff_strange_VTV_OS.err()<<endl;


      if(data_VKTK_OS.Tag[iens] == "cB211b.072.64") {

	cout<<"#### TESTING LINEAR ASSUMPTION IN ESTIMATING DERIVATIVE ##########"<<endl;
	cout<<"A(tm, OS) "<<susc_s_der_A_tm.ave()<<" +- "<<susc_s_der_A_tm.err()<<" , "<<susc_s_der_A_OS.ave()<<" +- "<<susc_s_der_A_OS.err()<<endl;
	cout<<"B(tm, OS) "<<susc_s_der_B_tm.ave()<<" +- "<<susc_s_der_B_tm.err()<<" , "<<susc_s_der_B_OS.ave()<<" +- "<<susc_s_der_B_OS.err()<<endl;
	cout<<"C(tm, OS) "<<susc_s_der_C_tm.ave()<<" +- "<<susc_s_der_C_tm.err()<<" , "<<susc_s_der_C_OS.ave()<<" +- "<<susc_s_der_C_OS.err()<<endl;
	cout<<"D(tm, OS) "<<susc_s_der_D_tm.ave()<<" +- "<<susc_s_der_D_tm.err()<<" , "<<susc_s_der_D_OS.ave()<<" +- "<<susc_s_der_D_OS.err()<<endl;
	cout<<"D(tm, OS) "<<susc_s_der_E_tm.ave()<<" +- "<<susc_s_der_E_tm.err()<<" , "<<susc_s_der_E_OS.ave()<<" +- "<<susc_s_der_E_OS.err()<<endl;
	cout<<"##################################################################"<<endl;
      }
      
    
    }

   
       
    //print susc info
    //light
    Print_To_File({}, {susc_VTV_tm_int.ave(), susc_VTV_tm_int.err(), susc_TV_tm_int.ave(), susc_TV_tm_int.err(), susc_VT_tm_int.ave(), susc_VT_tm_int.err()}, "../data/magnetic_susc/"+data_P5P5_OS.Tag[iens]+"/susc_data_tm.dat", "", "#VTV TV  VT");
    Print_To_File({}, {susc_VTV_OS_int.ave(), susc_VTV_OS_int.err(), susc_TV_OS_int.ave(), susc_TV_OS_int.err(), susc_VT_OS_int.ave(), susc_VT_OS_int.err()}, "../data/magnetic_susc/"+data_P5P5_OS.Tag[iens]+"/susc_data_OS.dat", "", "#VTV TV  VT");
    Print_To_File({}, {susc_sub_VTV_tm_int.ave(), susc_sub_VTV_tm_int.err(), susc_sub_TV_tm_int.ave(), susc_sub_TV_tm_int.err(), susc_sub_VT_tm_int.ave(), susc_sub_VT_tm_int.err()}, "../data/magnetic_susc/"+data_P5P5_OS.Tag[iens]+"/susc_sub_data_tm"+Tag_val+".dat", "", "#VTV TV  VT");
    Print_To_File({}, {susc_sub_VTV_OS_int.ave(), susc_sub_VTV_OS_int.err(), susc_sub_TV_OS_int.ave(), susc_sub_TV_OS_int.err(), susc_sub_VT_OS_int.ave(), susc_sub_VT_OS_int.err()}, "../data/magnetic_susc/"+data_P5P5_OS.Tag[iens]+"/susc_sub_data_OS"+Tag_val+".dat", "", "#VTV TV  VT");
    Print_To_File({}, {susc_der_VTV_tm_int.ave(), susc_der_VTV_tm_int.err(), susc_der_TV_tm_int.ave(), susc_der_TV_tm_int.err(), susc_der_VT_tm_int.ave(), susc_der_VT_tm_int.err()}, "../data/magnetic_susc/"+data_P5P5_OS.Tag[iens]+"/susc_der_data_tm"+Tag_val+".dat", "", "#VTV TV  VT");
    Print_To_File({}, {susc_der_VTV_OS_int.ave(), susc_der_VTV_OS_int.err(), susc_der_TV_OS_int.ave(), susc_der_TV_OS_int.err(), susc_der_VT_OS_int.ave(), susc_der_VT_OS_int.err()}, "../data/magnetic_susc/"+data_P5P5_OS.Tag[iens]+"/susc_der_data_OS"+Tag_val+".dat", "", "#VTV TV  VT");
    //strange 1
    Print_To_File({}, {susc_s1_VTV_tm_int.ave(), susc_s1_VTV_tm_int.err(), susc_s1_TV_tm_int.ave(), susc_s1_TV_tm_int.err(), susc_s1_VT_tm_int.ave(), susc_s1_VT_tm_int.err()}, "../data/magnetic_susc/"+data_P5P5_OS.Tag[iens]+"/susc_s1_data_tm.dat", "", "#VTV TV  VT");
    Print_To_File({}, {susc_s1_VTV_OS_int.ave(), susc_s1_VTV_OS_int.err(),  susc_s1_TV_OS_int.ave(), susc_s1_TV_OS_int.err(), susc_s1_VT_OS_int.ave(), susc_s1_VT_OS_int.err()}, "../data/magnetic_susc/"+data_P5P5_OS.Tag[iens]+"/susc_s1_data_OS.dat", "", "#VTV TV  VT");
    Print_To_File({}, {susc_s1_sub_VTV_tm_int.ave(), susc_s1_sub_VTV_tm_int.err(), susc_s1_sub_TV_tm_int.ave(), susc_s1_sub_TV_tm_int.err(), susc_s1_sub_VT_tm_int.ave(), susc_s1_sub_VT_tm_int.err()}, "../data/magnetic_susc/"+data_P5P5_OS.Tag[iens]+"/susc_s1_sub_data_tm"+Tag_val+".dat", "", "#VTV TV  VT");
    Print_To_File({}, {susc_s1_sub_VTV_OS_int.ave(), susc_s1_sub_VTV_OS_int.err(), susc_s1_sub_TV_OS_int.ave(), susc_s1_sub_TV_OS_int.err(), susc_s1_sub_VT_OS_int.ave(), susc_s1_sub_VT_OS_int.err()}, "../data/magnetic_susc/"+data_P5P5_OS.Tag[iens]+"/susc_s1_sub_data_OS"+Tag_val+".dat", "", "#VTV TV  VT");
    Print_To_File({}, {susc_s1_der_VTV_tm_int.ave(), susc_s1_der_VTV_tm_int.err(), susc_s1_der_TV_tm_int.ave(), susc_s1_der_TV_tm_int.err(), susc_s1_der_VT_tm_int.ave(), susc_s1_der_VT_tm_int.err()}, "../data/magnetic_susc/"+data_P5P5_OS.Tag[iens]+"/susc_s1_der_data_tm"+Tag_val+".dat", "", "#VTV TV  VT");
    Print_To_File({}, {susc_s1_der_VTV_OS_int.ave(), susc_s1_der_VTV_OS_int.err(), susc_s1_der_TV_OS_int.ave(), susc_s1_der_TV_OS_int.err(), susc_s1_der_VT_OS_int.ave(), susc_s1_der_VT_OS_int.err()}, "../data/magnetic_susc/"+data_P5P5_OS.Tag[iens]+"/susc_s1_der_data_OS"+Tag_val+".dat", "", "#VTV TV  VT");
    Print_To_File({}, {susc_s2_der_VTV_tm_int.ave(), susc_s2_der_VTV_tm_int.err(), susc_s2_der_TV_tm_int.ave(), susc_s2_der_TV_tm_int.err(), susc_s2_der_VT_tm_int.ave(), susc_s2_der_VT_tm_int.err()}, "../data/magnetic_susc/"+data_P5P5_OS.Tag[iens]+"/susc_s2_der_data_tm"+Tag_val+".dat", "", "#VTV TV  VT");
    Print_To_File({}, {susc_s2_der_VTV_OS_int.ave(), susc_s2_der_VTV_OS_int.err(), susc_s2_der_TV_OS_int.ave(), susc_s2_der_TV_OS_int.err(), susc_s2_der_VT_OS_int.ave(), susc_s2_der_VT_OS_int.err()}, "../data/magnetic_susc/"+data_P5P5_OS.Tag[iens]+"/susc_s2_der_data_OS"+Tag_val+".dat", "", "#VTV TV  VT");

    //print corr_info
    //light
    Print_To_File({}, {((-2*1000/a_distr)*Corr_VKTKVK_tm).ave(), ((-2*1000/a_distr)*Corr_VKTKVK_tm).err(), ((-2*1000/a_distr)*Corr_TKVK_tm).ave(), ((-2*1000/a_distr)*Corr_TKVK_tm).err(), ( (-2*1000/a_distr)*Corr_VKTK_tm).ave(), ( (-2*1000/(a_distr))*Corr_VKTK_tm).err()}, "../data/magnetic_susc/"+data_P5P5_OS.Tag[iens]+"/corr_data_tm.dat", "", "#VTV TV  VT");
    Print_To_File({}, {((-2*1000/a_distr)*Corr_VKTKVK_OS).ave(), ((-2*1000/a_distr)*Corr_VKTKVK_OS).err(), ((-2*1000/a_distr)*Corr_TKVK_OS).ave(), ((-2*1000/a_distr)*Corr_TKVK_OS).err(), ((-2*1000/a_distr)*Corr_VKTK_OS).ave(), ((-2*1000/a_distr)*Corr_VKTK_OS).err()}, "../data/magnetic_susc/"+data_P5P5_OS.Tag[iens]+"/corr_data_OS.dat", "", "#VTV TV  VT");
    
    Print_To_File({}, {((-2*1000/a_distr)*Corr_sub_VTV_tm).ave(), ((-2*1000/a_distr)*Corr_sub_VTV_tm).err(), ((-2*1000/a_distr)*Corr_sub_TV_tm).ave(), ((-2*1000/a_distr)*Corr_sub_TV_tm).err(), ( (-2*1000/a_distr)*Corr_sub_VT_tm).ave(), ( (-2*1000/(a_distr))*Corr_sub_VT_tm).err()}, "../data/magnetic_susc/"+data_P5P5_OS.Tag[iens]+"/corr_sub_data_tm"+Tag_val+".dat", "", "#VTV TV  VT");
    Print_To_File({}, {((-2*1000/a_distr)*Corr_sub_VTV_OS).ave(), ((-2*1000/a_distr)*Corr_sub_VTV_OS).err(), ((-2*1000/a_distr)*Corr_sub_TV_OS).ave(), ((-2*1000/a_distr)*Corr_sub_TV_OS).err(), ((-2*1000/a_distr)*Corr_sub_VT_OS).ave(), ((-2*1000/a_distr)*Corr_sub_VT_OS).err()}, "../data/magnetic_susc/"+data_P5P5_OS.Tag[iens]+"/corr_sub_data_OS"+Tag_val+".dat", "", "#VTV TV  VT");
    //strange
    Print_To_File({}, {((-2*1000/a_distr)*Corr_s1_VKTKVK_tm).ave(), ((-2*1000/a_distr)*Corr_s1_VKTKVK_tm).err(), ((-2*1000/a_distr)*Corr_s1_TKVK_tm).ave(), ((-2*1000/a_distr)*Corr_s1_TKVK_tm).err(), ((-2*1000/a_distr)*Corr_s1_VKTK_tm).ave(), ((-2*1000/a_distr)*Corr_s1_VKTK_tm).err()}, "../data/magnetic_susc/"+data_P5P5_OS.Tag[iens]+"/corr_s1_data_tm.dat", "", "#VTV TV  VT");
    Print_To_File({}, {((-2*1000/a_distr)*Corr_s1_VKTKVK_OS).ave(), ((-2*1000/a_distr)*Corr_s1_VKTKVK_OS).err(), ((-2*1000/a_distr)*Corr_s1_TKVK_OS).ave(), ((-2*1000/a_distr)*Corr_s1_TKVK_OS).err(), ((-2*1000/a_distr)*Corr_s1_VKTK_OS).ave(), ((-2*1000/a_distr)*Corr_s1_VKTK_OS).err()}, "../data/magnetic_susc/"+data_P5P5_OS.Tag[iens]+"/corr_s1_data_OS.dat", "", "#VTV TV  VT");
    Print_To_File({}, {((-2*1000/a_distr)*Corr_s1_sub_VTV_tm).ave(), ((-2*1000/a_distr)*Corr_s1_sub_VTV_tm).err(),  ((-2*1000/a_distr)*Corr_s1_sub_TV_tm).ave(), ((-2*1000/a_distr)*Corr_s1_sub_TV_tm).err(), ((-2*1000/a_distr)*Corr_s1_sub_VT_tm).ave(), ((-2*1000/a_distr)*Corr_s1_sub_VT_tm).err()}, "../data/magnetic_susc/"+data_P5P5_OS.Tag[iens]+"/corr_s1_sub_data_tm"+Tag_val+".dat", "", "#VTV TV  VT");
    Print_To_File({}, {((-2*1000/a_distr)*Corr_s1_sub_VTV_OS).ave(), ((-2*1000/a_distr)*Corr_s1_sub_VTV_OS).err(), ((-2*1000/a_distr)*Corr_s1_sub_TV_OS).ave(), ((-2*1000/a_distr)*Corr_s1_sub_TV_OS).err(), ((-2*1000/a_distr)*Corr_s1_sub_VT_OS).ave(), ((-2*1000/a_distr)*Corr_s1_sub_VT_OS).err()}, "../data/magnetic_susc/"+data_P5P5_OS.Tag[iens]+"/corr_s1_sub_data_OS"+Tag_val+".dat", "", "#VTV TV  VT");

    //push_back light info
    susc_TV_tm_list.distr_list.push_back( susc_TV_tm);
    susc_VT_tm_list.distr_list.push_back( susc_VT_tm);
    susc_VTV_tm_list.distr_list.push_back( susc_VTV_tm);
   
    susc_TV_OS_list.distr_list.push_back(susc_TV_OS);
    susc_VT_OS_list.distr_list.push_back(susc_VT_OS);
    susc_VTV_OS_list.distr_list.push_back(susc_VTV_OS);

     
    
    
  
    //interpolate at physical strange quark mass
    //set time intervals for pseudoscalar obs
    if(data_TKTK_tm.Tag[iens].substr(1,1) == "C") {
      if(data_TKTK_tm.Tag[iens]=="cC211a.06.80") { Corr.Tmin=40; Corr.Tmax=70;}
      else crash("Cannot find ensemble tag: "+data_TKTK_tm.Tag[iens]);
    }
    else if(data_TKTK_tm.Tag[iens].substr(1,1) == "B") {
      if(data_TKTK_tm.Tag[iens]== "cB211a.14.64") {Corr.Tmin=31; Corr.Tmax=58;}
      else if(data_TKTK_tm.Tag[iens] == "cB211a.25.48") {Corr.Tmin=23;Corr.Tmax=44;}
      else if(data_TKTK_tm.Tag[iens] == "cB211b.072.64") {Corr.Tmin=36; Corr.Tmax= 57;}
      else if(data_TKTK_tm.Tag[iens] == "cB211b.072.96") {Corr.Tmin=40; Corr.Tmax= 80;}
      else crash("Cannot find ensemble tag: "+data_TKTK_tm.Tag[iens]);
    }
    else if(data_TKTK_tm.Tag[iens].substr(1,1) == "A") {
      if(data_TKTK_tm.Tag[iens] == "cA211a.12.48") {Corr.Tmin=19; Corr.Tmax=33;}
      else if(data_TKTK_tm.Tag[iens] == "cA211a.40.24") {Corr.Tmin=18; Corr.Tmax=23;}
      else if(data_TKTK_tm.Tag[iens] == "cA211a.53.24") {Corr.Tmin=16; Corr.Tmax=22;}
      else if(data_TKTK_tm.Tag[iens] == "cA211ab.30.32") {Corr.Tmin=23; Corr.Tmax=30;}
      else crash("Cannot find ensemble tag: "+data_TKTK_tm.Tag[iens]);
    }
    else if(data_TKTK_tm.Tag[iens].substr(1,1) == "D") {
      if(data_TKTK_tm.Tag[iens] == "cD211a.054.96") {Corr.Tmin=55; Corr.Tmax=88;}
      else crash("Cannot find ensemble tag: "+data_TKTK_tm.Tag[iens]);
    }
    else crash("Ensemble tag not valid");

    distr_t etas1_M = Corr.Fit_distr( Corr.effective_mass_t(Corr_s1_P5P5_tm, ""))/a_distr;
    distr_t etas2_M = Corr.Fit_distr( Corr.effective_mass_t(Corr_s2_P5P5_tm, ""))/a_distr;
    distr_t etas_phys(UseJack);
    for(int ijack=0;ijack<Njacks;ijack++) etas_phys.distr.push_back( 0.68989 + GM()*0.00050/sqrt(Njacks-1.0) );
    distr_t etas_phys2= etas_phys*etas_phys;
    vector<distr_t> X_2_fit({etas1_M*etas1_M, etas2_M*etas2_M});
    vector<distr_t> susc_TV_tm_2_fit({susc_s1_TV_tm, susc_s2_TV_tm});
    vector<distr_t> susc_VT_tm_2_fit({susc_s1_VT_tm, susc_s2_VT_tm});
    vector<distr_t> susc_VTV_tm_2_fit({susc_s1_VTV_tm, susc_s2_VTV_tm});
 
    vector<distr_t> susc_TV_OS_2_fit({susc_s1_TV_OS, susc_s2_TV_OS});
    vector<distr_t> susc_VT_OS_2_fit({susc_s1_VT_OS, susc_s2_VT_OS});
    vector<distr_t> susc_VTV_OS_2_fit({susc_s1_VTV_OS, susc_s2_VTV_OS});
    
    
    
    susc_s_TV_tm_list.distr_list.push_back( Obs_extrapolation_meson_mass(susc_TV_tm_2_fit, X_2_fit, etas_phys2 ,  "../data/magnetic_susc", "susc_s_TV_tm_"+data_TKTK_tm.Tag[iens], UseJack, "SPLINE"));
    susc_s_VT_tm_list.distr_list.push_back( Obs_extrapolation_meson_mass(susc_VT_tm_2_fit, X_2_fit, etas_phys2,  "../data/magnetic_susc", "susc_s_VT_tm_"+data_TKTK_tm.Tag[iens], UseJack, "SPLINE"));
    susc_s_VTV_tm_list.distr_list.push_back( Obs_extrapolation_meson_mass(susc_VTV_tm_2_fit, X_2_fit, etas_phys2,  "../data/magnetic_susc", "susc_s_VTV_tm_"+data_TKTK_tm.Tag[iens], UseJack, "SPLINE"));

   
       
    susc_s_TV_OS_list.distr_list.push_back( Obs_extrapolation_meson_mass(susc_TV_OS_2_fit, X_2_fit, etas_phys2,  "../data/magnetic_susc", "susc_s_TV_OS_"+data_TKTK_OS.Tag[iens], UseJack, "SPLINE"));
    susc_s_VT_OS_list.distr_list.push_back( Obs_extrapolation_meson_mass(susc_VT_OS_2_fit, X_2_fit, etas_phys2,  "../data/magnetic_susc", "susc_s_VT_OS_"+data_TKTK_OS.Tag[iens], UseJack, "SPLINE"));
    susc_s_VTV_OS_list.distr_list.push_back( Obs_extrapolation_meson_mass(susc_VTV_OS_2_fit, X_2_fit, etas_phys2,  "../data/magnetic_susc", "susc_s_VTV_OS_"+data_TKTK_OS.Tag[iens], UseJack, "SPLINE"));

    //interpolate unsubtracted strange correlator (tm-only for figure paper)

    distr_t_list Corr_s_VKTKVK_tm_interpol(UseJack);
    for(int t=0;t<Corr.Nt;t++) {
      vector<distr_t> Corr_fit({ Corr_s1_VKTKVK_tm.distr_list[t], Corr_s2_VKTKVK_tm.distr_list[t] });
      Corr_s_VKTKVK_tm_interpol.distr_list.push_back( Obs_extrapolation_meson_mass(Corr_fit, X_2_fit, etas_phys2,  "../data/magnetic_susc", "corr_s_VTV_tm_t_"+to_string(t)+"_"+data_TKTK_OS.Tag[iens], UseJack, "SPLINE"));
    }

        
    distr_t_list susc_s_int_VTV_tm_interpol(UseJack, Corr.Nt/2 +1,  UseJack?Njacks:800);
    for(int t=1; t < Corr.Nt/2 -10; t++) {
      susc_s_int_VTV_tm_interpol.distr_list[t] = susc_s_int_VTV_tm_interpol.distr_list[t-1]  -1000*(2/a_distr)*Corr_s_VKTKVK_tm_interpol.distr_list[t]*t;
    }

    Print_To_File({}, {((-2*1000/a_distr)*Corr_s_VKTKVK_tm_interpol).ave(), ((-2*1000/a_distr)*Corr_s_VKTKVK_tm_interpol).err()}, "../data/magnetic_susc/"+data_P5P5_OS.Tag[iens]+"/corr_s_interpol_data_tm.dat", "", "#VTV");
    Print_To_File({}, {susc_s_int_VTV_tm_interpol.ave(), susc_s_int_VTV_tm_interpol.err()}, "../data/magnetic_susc/"+data_P5P5_OS.Tag[iens]+"/susc_s_interpol_data_tm.dat", "", "#VTV");
    

    
      
    
    
    a_distr_list.distr_list.push_back( a_distr/fmTGeV);
    Ensemble_list.push_back( data_TKTK_tm.Tag[iens]);
    
    
    

    
    

  }

  //Print to File
  boost::filesystem::create_directory("../data/magnetic_susc/cont");
  Print_To_File(Ensemble_list, { a_distr_list.ave(), susc_TV_tm_list.ave(), susc_TV_tm_list.err()}, "../data/magnetic_susc/cont/susc_ll_TV_tm.dat.t", "", "");
  Print_To_File(Ensemble_list, { a_distr_list.ave(), susc_VT_tm_list.ave(), susc_VT_tm_list.err()}, "../data/magnetic_susc/cont/susc_ll_VT_tm.dat.t", "", "");
  Print_To_File(Ensemble_list, { a_distr_list.ave(), susc_VTV_tm_list.ave(), susc_VTV_tm_list.err()}, "../data/magnetic_susc/cont/susc_ll_VTV_tm.dat.t", "", "");
 
  Print_To_File(Ensemble_list, { a_distr_list.ave(), susc_TV_OS_list.ave(), susc_TV_OS_list.err()}, "../data/magnetic_susc/cont/susc_ll_TV_OS.dat.t", "", "");
  Print_To_File(Ensemble_list, { a_distr_list.ave(), susc_VT_OS_list.ave(), susc_VT_OS_list.err()}, "../data/magnetic_susc/cont/susc_ll_VT_OS.dat.t", "", "");
  Print_To_File(Ensemble_list, { a_distr_list.ave(), susc_VTV_OS_list.ave(), susc_VTV_OS_list.err()}, "../data/magnetic_susc/cont/susc_ll_VTV_OS.dat.t", "", "");
 

  Print_To_File(Ensemble_list, { a_distr_list.ave(), susc_s_TV_tm_list.ave(), susc_s_TV_tm_list.err()}, "../data/magnetic_susc/cont/susc_s_TV_tm.dat.t", "", "");
  Print_To_File(Ensemble_list, { a_distr_list.ave(), susc_s_VT_tm_list.ave(), susc_s_VT_tm_list.err()}, "../data/magnetic_susc/cont/susc_s_VT_tm.dat.t", "", "");
  Print_To_File(Ensemble_list, { a_distr_list.ave(), susc_s_VTV_tm_list.ave(), susc_s_VTV_tm_list.err()}, "../data/magnetic_susc/cont/susc_s_VTV_tm.dat.t", "", "");
 
  Print_To_File(Ensemble_list, { a_distr_list.ave(), susc_s_TV_OS_list.ave(), susc_s_TV_OS_list.err()}, "../data/magnetic_susc/cont/susc_s_TV_OS.dat.t", "", "");
  Print_To_File(Ensemble_list, { a_distr_list.ave(), susc_s_VT_OS_list.ave(), susc_s_VT_OS_list.err()}, "../data/magnetic_susc/cont/susc_s_VT_OS.dat.t", "", "");
  Print_To_File(Ensemble_list, { a_distr_list.ave(), susc_s_VTV_OS_list.ave(), susc_s_VTV_OS_list.err()}, "../data/magnetic_susc/cont/susc_s_VTV_OS.dat.t", "", "");


  


 


 
  
  //perform infinite volume extrapolation
  

  //find B64 and B96 ensembles 
  int B64=-1;
  int B96=-1;

  for(int iens=0; iens<Nens; iens++) {
    if( data_VKTK_tm.Tag[iens] == "cB211b.072.64") B64=iens;
    if( data_VKTK_tm.Tag[iens] == "cB211b.072.96") B96=iens;
  }

  double sigma_comb_TV_tm= sqrt( pow(susc_TV_tm_list.err(B64),2) + pow(susc_TV_tm_list.err(B96),2));
  double sigma_comb_VT_tm= sqrt( pow(susc_VT_tm_list.err(B64),2) + pow(susc_VT_tm_list.err(B96),2));
  double sigma_comb_VTV_tm= sqrt( pow(susc_VTV_tm_list.err(B64),2) + pow(susc_VTV_tm_list.err(B96),2));

  
  double sigma_comb_TV_OS= sqrt( pow(susc_TV_OS_list.err(B64),2) + pow(susc_TV_OS_list.err(B96),2));
  double sigma_comb_VT_OS= sqrt( pow(susc_VT_OS_list.err(B64),2) + pow(susc_VT_OS_list.err(B96),2));
  double sigma_comb_VTV_OS= sqrt( pow(susc_VTV_OS_list.err(B64),2) + pow(susc_VTV_OS_list.err(B96),2));
  
    
  double syst_TV_tm = (fabs(susc_TV_tm_list.ave(B64) - susc_TV_tm_list.ave(B96))/fabs(susc_TV_tm_list.ave(B64)))*erf( fabs(susc_TV_tm_list.ave(B64) - susc_TV_tm_list.ave(B96))/(sqrt(2)*sigma_comb_TV_tm));
  double syst_VT_tm = (fabs(susc_VT_tm_list.ave(B64) - susc_VT_tm_list.ave(B96))/fabs(susc_VT_tm_list.ave(B64)))*erf( fabs(susc_VT_tm_list.ave(B64) - susc_VT_tm_list.ave(B96))/(sqrt(2)*sigma_comb_VT_tm));
  double syst_VTV_tm = (fabs(susc_VTV_tm_list.ave(B64) - susc_VTV_tm_list.ave(B96))/fabs(susc_VTV_tm_list.ave(B64)))*erf( fabs(susc_VTV_tm_list.ave(B64) - susc_VTV_tm_list.ave(B96))/(sqrt(2)*sigma_comb_VTV_tm));

  

  double syst_TV_OS = (fabs(susc_TV_OS_list.ave(B64) - susc_TV_OS_list.ave(B96))/fabs(susc_TV_OS_list.ave(B64)))*erf( fabs(susc_TV_OS_list.ave(B64) - susc_TV_OS_list.ave(B96))/(sqrt(2)*sigma_comb_TV_OS));
  double syst_VT_OS = (fabs(susc_VT_OS_list.ave(B64) - susc_VT_OS_list.ave(B96))/fabs(susc_VT_OS_list.ave(B64)))*erf( fabs(susc_VT_OS_list.ave(B64) - susc_VT_OS_list.ave(B96))/(sqrt(2)*sigma_comb_VT_OS));
  double syst_VTV_OS = (fabs(susc_VTV_OS_list.ave(B64) - susc_VTV_OS_list.ave(B96))/fabs(susc_VTV_OS_list.ave(B64)))*erf( fabs(susc_VTV_OS_list.ave(B64) - susc_VTV_OS_list.ave(B96))/(sqrt(2)*sigma_comb_VTV_OS));

   
  double syst_TV = max(syst_TV_tm, syst_TV_OS);
  double syst_VT = max(syst_VT_tm, syst_VT_OS);
  double syst_VTV = max(syst_VTV_tm, syst_VTV_OS);
  

  distr_t distr_syst_FSE_TV(UseJack), distr_syst_FSE_VT(UseJack), distr_syst_FSE_VTV(UseJack);
  
 
  
  for(int ijack=0;ijack<Njacks;ijack++) {
    distr_syst_FSE_TV.distr.push_back( 1.0 + 0.0*GM()*syst_TV/sqrt(Njacks-1.0));
    distr_syst_FSE_VT.distr.push_back( 1.0 + 0.0*GM()*syst_VT/sqrt(Njacks-1.0));
    distr_syst_FSE_VTV.distr.push_back( 1.0 + 0.0*GM()*syst_VTV/sqrt(Njacks-1.0));
  
  }


  for(int iens=0;iens<Nens;iens++) {

    //to interpolate between B64 and B96
    distr_t corr_FSE_light_TV_tm = Get_id_jack_distr( Njacks);
    distr_t corr_FSE_light_VT_tm = Get_id_jack_distr(Njacks);
    distr_t corr_FSE_light_VTV_tm = Get_id_jack_distr(Njacks);
    distr_t corr_FSE_light_TV_OS = Get_id_jack_distr(Njacks);
    distr_t corr_FSE_light_VT_OS = Get_id_jack_distr(Njacks);
    distr_t corr_FSE_light_VTV_OS = Get_id_jack_distr(Njacks);
      
    if(data_TKTK_tm.Tag[iens] != "cB211b.072.96") {


      if(data_TKTK_tm.Tag[iens] == "cB211b.072.64") {

	double DC= exp(-5.43)-exp(-3.62);
	double DC_ref= exp(-3.62*5.46/5.09 ) -exp(-3.62);
	double DC_rel= DC_ref/DC;
	//find B96 tag
	int iB=-1;
	for(int jens=0;jens<Nens;jens++) if( data_TKTK_tm.Tag[jens] == "cB211b.072.96") iB=jens;
	if(iB == -1) crash("iB = -1");
	corr_FSE_light_TV_tm = corr_FSE_light_TV_tm + DC_rel*(susc_TV_tm_list.distr_list[iB]/susc_TV_tm_list.distr_list[iens] -1);
	corr_FSE_light_VT_tm = corr_FSE_light_VT_tm + DC_rel*(susc_VT_tm_list.distr_list[iB]/susc_VT_tm_list.distr_list[iens] -1);
	corr_FSE_light_VTV_tm = corr_FSE_light_VTV_tm + DC_rel*(susc_VTV_tm_list.distr_list[iB]/susc_VTV_tm_list.distr_list[iens] -1);

	corr_FSE_light_TV_OS = corr_FSE_light_TV_OS +  DC_rel*(susc_TV_OS_list.distr_list[iB]/susc_TV_OS_list.distr_list[iens] -1);
	corr_FSE_light_VT_OS = corr_FSE_light_VT_OS + DC_rel*(susc_VT_OS_list.distr_list[iB]/susc_VT_OS_list.distr_list[iens] -1);
	corr_FSE_light_VTV_OS = corr_FSE_light_VTV_OS+ DC_rel*(susc_VTV_OS_list.distr_list[iB]/susc_VTV_OS_list.distr_list[iens] -1);



      }

      susc_TV_tm_list_red.distr_list.push_back( susc_TV_tm_list.distr_list[iens]*corr_FSE_light_TV_tm);
      susc_VT_tm_list_red.distr_list.push_back( susc_VT_tm_list.distr_list[iens]*corr_FSE_light_VT_tm);
      susc_VTV_tm_list_red.distr_list.push_back( susc_VTV_tm_list.distr_list[iens]*corr_FSE_light_VTV_tm);
    

      susc_TV_OS_list_red.distr_list.push_back( susc_TV_OS_list.distr_list[iens]*corr_FSE_light_TV_OS);
      susc_VT_OS_list_red.distr_list.push_back( susc_VT_OS_list.distr_list[iens]*corr_FSE_light_VT_OS);
      susc_VTV_OS_list_red.distr_list.push_back( susc_VTV_OS_list.distr_list[iens]*corr_FSE_light_VTV_OS);
   

      susc_s_TV_tm_list_red.distr_list.push_back( susc_s_TV_tm_list.distr_list[iens]);
      susc_s_VT_tm_list_red.distr_list.push_back( susc_s_VT_tm_list.distr_list[iens]);
      susc_s_VTV_tm_list_red.distr_list.push_back( susc_s_VTV_tm_list.distr_list[iens]);
     
      
      susc_s_TV_OS_list_red.distr_list.push_back( susc_s_TV_OS_list.distr_list[iens]);
      susc_s_VT_OS_list_red.distr_list.push_back( susc_s_VT_OS_list.distr_list[iens]);
      susc_s_VTV_OS_list_red.distr_list.push_back( susc_s_VTV_OS_list.distr_list[iens]);
    

      
      
      a_distr_list_red.distr_list.push_back( a_distr_list.distr_list[iens]);
      Ensemble_list_red.push_back(  data_TKTK_tm.Tag[iens]);
    }
  }


  //print reduced points for light
  Print_To_File(Ensemble_list_red, { a_distr_list_red.ave(), susc_TV_tm_list_red.ave(), susc_TV_tm_list_red.err()}, "../data/magnetic_susc/cont/susc_ll_red_TV_tm.dat.t", "", "");
  Print_To_File(Ensemble_list_red, { a_distr_list_red.ave(), susc_VT_tm_list_red.ave(), susc_VT_tm_list_red.err()}, "../data/magnetic_susc/cont/susc_ll_red_VT_tm.dat.t", "", "");
  Print_To_File(Ensemble_list_red, { a_distr_list_red.ave(), susc_VTV_tm_list_red.ave(), susc_VTV_tm_list_red.err()}, "../data/magnetic_susc/cont/susc_ll_red_VTV_tm.dat.t", "", "");

  Print_To_File(Ensemble_list_red, { a_distr_list_red.ave(), susc_TV_OS_list_red.ave(), susc_TV_OS_list_red.err()}, "../data/magnetic_susc/cont/susc_ll_red_TV_OS.dat.t", "", "");
  Print_To_File(Ensemble_list_red, { a_distr_list_red.ave(), susc_VT_OS_list_red.ave(), susc_VT_OS_list_red.err()}, "../data/magnetic_susc/cont/susc_ll_red_VT_OS.dat.t", "", "");
  Print_To_File(Ensemble_list_red, { a_distr_list_red.ave(), susc_VTV_OS_list_red.ave(), susc_VTV_OS_list_red.err()}, "../data/magnetic_susc/cont/susc_ll_red_VTV_OS.dat.t", "", "");
 

  

   
 

  //perform continuum limit extrapolation

  
  int Nens_eff=Nens-1;

     

  class ipar_chi {
	
  public:
    ipar_chi()  {}
	
    double chi, chi_err, a; //lattice spacing a is in fm
    bool Is_tm;
  };
  
      
  class fpar_chi {
    
  public:
    fpar_chi() {}
    fpar_chi(const Vfloat &par) {
      if((signed)par.size() != 3) crash("In class fpar_chi, class constructor Vfloat par has size != 3");
      D=par[0];
      D2_tm=par[1];
      D2_OS=par[2];
    }
    
    double D,D2_tm, D2_OS;
  };

  vector<string> flavs({"light", "strange"});
  vector<distr_t_list> susc_flav_TV_tm({susc_TV_tm_list_red, susc_s_TV_tm_list_red});
  vector<distr_t_list> susc_flav_VT_tm({susc_VT_tm_list_red, susc_s_VT_tm_list_red});
  vector<distr_t_list> susc_flav_VTV_tm({susc_VTV_tm_list_red, susc_s_VTV_tm_list_red});
  vector<distr_t_list> susc_flav_TV_OS({susc_TV_OS_list_red, susc_s_TV_OS_list_red});
  vector<distr_t_list> susc_flav_VT_OS({susc_VT_OS_list_red, susc_s_VT_OS_list_red});
  vector<distr_t_list> susc_flav_VTV_OS({susc_VTV_OS_list_red, susc_s_VTV_OS_list_red});
  
  


  vector<string> contribs({"VTV"});

  VVfloat ch2_ndof(flavs.size());
  for(auto &ch : ch2_ndof) ch.resize(contribs.size());


  for(int iflav=0; iflav<(signed)flavs.size(); iflav++) {

  
    int icontr=-1;

 

    for(auto &contr: contribs) {
      icontr++;
      cout<<"###########################################################"<<endl;
      cout<<"Performing continuum limit extrapolation for type: "<<flavs[iflav]<<" Dirac structure: "<<contr<<endl;
      cout<<"Fit_type: combined fit"<<endl;

       
     int Nmeas= 2*Nens_eff;
     int Npars= 3;
     int Ndof= Nmeas-Npars;

     cout<<"Nmeas: "<<Nmeas<<endl;
     cout<<"Npars: "<<Npars<<endl;
     cout<<"Ndof: "<<Ndof<<endl;
     cout<<"Nens: "<<Nens_eff<<endl;
	    
     bootstrap_fit<fpar_chi,ipar_chi> bf_chi(Njacks);
     bootstrap_fit<fpar_chi,ipar_chi> bf_chi_ch2(1);
     //bf_TAU.Disable_correlated_fit();
     //bf_TAU_ch2.Disable_correlated_fit();
     bf_chi.Set_number_of_measurements(Nmeas);
     bf_chi.Set_verbosity(1);
     //ch2
     bf_chi_ch2.Set_number_of_measurements(Nmeas);
     bf_chi_ch2.Set_verbosity(1);

     //add fit parameters
     bf_chi.Add_par("D", 3.0, 0.1);
     bf_chi.Add_par("D2_tm", 2, 0.1);
     bf_chi.Add_par("D2_OS", 2, 0.1);
     //ch2
     bf_chi_ch2.Add_par("D", 3.0, 0.1);
     bf_chi_ch2.Add_par("D2_tm", 2, 0.1);
     bf_chi_ch2.Add_par("D2_OS", 2, 0.1);


     //ansatz
     bf_chi.ansatz=  [ ](const fpar_chi &p, const ipar_chi &ip) {
	      double D2=0.0;
	      if( ip.Is_tm==true ) D2=p.D2_tm;
	      else D2=p.D2_OS;
	      
	      return p.D + D2*pow(ip.a*QCD_scale,2);
	    };
     //meas
     bf_chi.measurement=  [ ](const fpar_chi &p, const ipar_chi &ip) {
	      return ip.chi;
	    };
     //err
     bf_chi.error=  [ ](const fpar_chi &p, const ipar_chi &ip) {
	      return ip.chi_err;
     };
     //ch2
     bf_chi_ch2.ansatz= bf_chi.ansatz;
     bf_chi_ch2.measurement= bf_chi.measurement;
     bf_chi_ch2.error= bf_chi.error;

     //fill the data
     int off_OS = Nens_eff;
     vector<vector<ipar_chi>> data(Njacks);
     vector<vector<ipar_chi>> data_ch2(1);
     //allocate space for output result
     boot_fit_data<fpar_chi> Bt_fit;
     boot_fit_data<fpar_chi> Bt_fit_ch2;
     for(auto &data_iboot: data) data_iboot.resize(Nmeas);
     for(auto &data_iboot: data_ch2) data_iboot.resize(Nmeas);



     //add covariance matrix
     Eigen::MatrixXd Cov_Matrix(Nmeas,Nmeas);
     Eigen::MatrixXd Corr_Matrix(Nmeas,Nmeas);
     for(int i=0;i<Nmeas;i++) for(int j=0;j<Nmeas;j++) {Cov_Matrix(i,j)=0; Corr_Matrix(i,j)=0;}

     
     //compute cov matrix between tm and OS
     for(int iens=0; iens<Nens_eff;iens++) {
       
       Corr_Matrix(iens,iens) = 1;
       Corr_Matrix(iens+off_OS,iens+off_OS) = 1;
       
       if(contr=="TV") {
	 Cov_Matrix(iens,iens) = pow(susc_flav_TV_tm[iflav].err(iens),2); Cov_Matrix(iens+off_OS,iens+off_OS) = pow(susc_flav_TV_OS[iflav].err(iens),2);
	 Cov_Matrix(iens, off_OS+iens) = (susc_flav_TV_tm[iflav].distr_list[iens]%susc_flav_TV_OS[iflav].distr_list[iens]);
       }
       else if(contr=="VT") {
	 Cov_Matrix(iens,iens) = pow(susc_flav_VT_tm[iflav].err(iens),2); Cov_Matrix(iens+off_OS,iens+off_OS) = pow(susc_flav_VT_OS[iflav].err(iens),2); 
	 Cov_Matrix(iens, off_OS+iens) = (susc_flav_VT_tm[iflav].distr_list[iens]%susc_flav_VT_OS[iflav].distr_list[iens]);
       }

       else if(contr=="VTV") {
	 Cov_Matrix(iens,iens) = pow(susc_flav_VTV_tm[iflav].err(iens),2); Cov_Matrix(iens+off_OS,iens+off_OS) = pow(susc_flav_VTV_OS[iflav].err(iens),2); 
	 Cov_Matrix(iens, off_OS+iens) = (susc_flav_VTV_tm[iflav].distr_list[iens]%susc_flav_VTV_OS[iflav].distr_list[iens]);
       }

     
       else crash("contr: "+contr+" not recognized");

       Corr_Matrix(iens, off_OS+iens) = Cov_Matrix(iens,off_OS+iens)/sqrt( Cov_Matrix(iens,iens)*Cov_Matrix(off_OS+iens,off_OS+iens));
       
       
       //symmetrize
       Cov_Matrix(off_OS+iens, iens) = Cov_Matrix(iens, off_OS+iens);
       Corr_Matrix(off_OS+iens, iens) = Corr_Matrix(iens, off_OS+iens);
       


     }

     //add cov matrix to bootstrap fit
     bf_chi.Add_covariance_matrix(Cov_Matrix);
     bf_chi_ch2.Add_covariance_matrix(Cov_Matrix);

     //print covariance matrix
     boost::filesystem::create_directory("../data/magnetic_susc/cont/cov");
     boost::filesystem::create_directory("../data/magnetic_susc/cont/corr");
     ofstream Print_Cov("../data/magnetic_susc/cont/cov/"+flavs[iflav]+"_"+contr+".cov");
     ofstream Print_Corr("../data/magnetic_susc/cont/corr/"+flavs[iflav]+"_"+contr+".corr");

     Print_Cov<<Cov_Matrix<<endl;  Print_Corr<<Corr_Matrix<<endl;
     Print_Cov.close();            Print_Corr.close();



     //fill with data
      for(int ijack=0;ijack<Njacks;ijack++) {
	      for(int iens=0;iens<Nens_eff;iens++) {

		if(contr=="TV") {
		   data[ijack][iens].chi = susc_flav_TV_tm[iflav].distr_list[iens].distr[ijack];
		   data[ijack][iens].chi_err= susc_flav_TV_tm[iflav].err(iens);
		   data[ijack][iens].Is_tm = true;
		   data[ijack][iens].a = a_distr_list_red.distr_list[iens].distr[ijack];
		   data[ijack][iens+off_OS].chi = susc_flav_TV_OS[iflav].distr_list[iens].distr[ijack];
		   data[ijack][iens+off_OS].chi_err= susc_flav_TV_OS[iflav].err(iens);
		   data[ijack][iens+off_OS].Is_tm = false;
		   data[ijack][iens+off_OS].a = a_distr_list_red.distr_list[iens].distr[ijack];
		}

		else if (contr =="VT") {

		  data[ijack][iens].chi = susc_flav_VT_tm[iflav].distr_list[iens].distr[ijack];
		  data[ijack][iens].chi_err= susc_flav_VT_tm[iflav].err(iens);
		  data[ijack][iens].Is_tm = true;
		  data[ijack][iens].a = a_distr_list_red.distr_list[iens].distr[ijack];
		  data[ijack][iens+off_OS].chi = susc_flav_VT_OS[iflav].distr_list[iens].distr[ijack];
		  data[ijack][iens+off_OS].chi_err= susc_flav_VT_OS[iflav].err(iens);
		  data[ijack][iens+off_OS].Is_tm = false;
		  data[ijack][iens+off_OS].a = a_distr_list_red.distr_list[iens].distr[ijack];
		  
		}

		else if (contr =="VTV") {

		  data[ijack][iens].chi = susc_flav_VTV_tm[iflav].distr_list[iens].distr[ijack];
		  data[ijack][iens].chi_err= susc_flav_VTV_tm[iflav].err(iens);
		  data[ijack][iens].Is_tm = true;
		  data[ijack][iens].a = a_distr_list_red.distr_list[iens].distr[ijack];
		  data[ijack][iens+off_OS].chi = susc_flav_VTV_OS[iflav].distr_list[iens].distr[ijack];
		  data[ijack][iens+off_OS].chi_err= susc_flav_VTV_OS[iflav].err(iens);
		  data[ijack][iens+off_OS].Is_tm = false;
		  data[ijack][iens+off_OS].a = a_distr_list_red.distr_list[iens].distr[ijack];
		  
		}

		
        
		else crash("contr : "+contr+" not yet implemented");


		//mean values
		if(ijack==0) {

		  if(contr=="TV") {
		    data_ch2[ijack][iens].chi = susc_flav_TV_tm[iflav].ave(iens);
		    data_ch2[ijack][iens].chi_err= susc_flav_TV_tm[iflav].err(iens);
		    data_ch2[ijack][iens].Is_tm = true;
		    data_ch2[ijack][iens].a = a_distr_list_red.ave(iens);
		    data_ch2[ijack][iens+off_OS].chi = susc_flav_TV_OS[iflav].ave(iens);
		    data_ch2[ijack][iens+off_OS].chi_err= susc_flav_TV_OS[iflav].err(iens);
		    data_ch2[ijack][iens+off_OS].Is_tm = false;
		    data_ch2[ijack][iens+off_OS].a = a_distr_list_red.ave(iens);
		  }

		  else if (contr =="VT") {

		    data_ch2[ijack][iens].chi = susc_flav_VT_tm[iflav].ave(iens);
		    data_ch2[ijack][iens].chi_err= susc_flav_VT_tm[iflav].err(iens);
		    data_ch2[ijack][iens].Is_tm = true;
		    data_ch2[ijack][iens].a = a_distr_list_red.ave(iens);
		    data_ch2[ijack][iens+off_OS].chi = susc_flav_VT_OS[iflav].ave(iens);
		    data_ch2[ijack][iens+off_OS].chi_err= susc_flav_VT_OS[iflav].err(iens);
		    data_ch2[ijack][iens+off_OS].Is_tm = false;
		    data_ch2[ijack][iens+off_OS].a = a_distr_list_red.ave(iens);
		    
		  }

		  else if (contr =="VTV") {

		    data_ch2[ijack][iens].chi = susc_flav_VTV_tm[iflav].ave(iens);
		    data_ch2[ijack][iens].chi_err= susc_flav_VTV_tm[iflav].err(iens);
		    data_ch2[ijack][iens].Is_tm = true;
		    data_ch2[ijack][iens].a = a_distr_list_red.ave(iens);
		    data_ch2[ijack][iens+off_OS].chi = susc_flav_VTV_OS[iflav].ave(iens);
		    data_ch2[ijack][iens+off_OS].chi_err= susc_flav_VTV_OS[iflav].err(iens);
		    data_ch2[ijack][iens+off_OS].Is_tm = false;
		    data_ch2[ijack][iens+off_OS].a = a_distr_list_red.ave(iens);
		    
		  }

	        

		  else crash("contr : "+contr+" not yet implemented");

		  
		}
		
	      }
      }


      //append
      bf_chi.Append_to_input_par(data);
      bf_chi_ch2.Append_to_input_par(data_ch2);
      //fit
      cout<<"Fitting...."<<endl;
      Bt_fit= bf_chi.Perform_bootstrap_fit();
      Bt_fit_ch2= bf_chi_ch2.Perform_bootstrap_fit();

	   
      //retrieve parameters
      distr_t D(UseJack), D2_tm(UseJack), D2_OS(UseJack);
      for(int ijack=0;ijack<Njacks;ijack++) { D.distr.push_back( Bt_fit.par[ijack].D); D2_tm.distr.push_back( Bt_fit.par[ijack].D2_tm); D2_OS.distr.push_back( Bt_fit.par[ijack].D2_OS);}
      //reduced ch2
      double ch2= Bt_fit_ch2.get_ch2_ave()/Ndof;

      ch2_ndof[iflav][icontr] = ch2;

      //print fit function
      distr_t_list chi_tm_to_print(UseJack), chi_OS_to_print(UseJack);

      int Nlat=500;
      Vfloat a_to_print_chi;
      double sx= 0.08*1.5/(Nlat-1.0); //fm
      for(int pp=0;pp<Nlat;pp++) a_to_print_chi.push_back( sx*pp);
      
      for(auto &a: a_to_print_chi) {
	chi_tm_to_print.distr_list.push_back( D + D2_tm*pow(a*QCD_scale,2));
	chi_OS_to_print.distr_list.push_back( D + D2_OS*pow(a*QCD_scale,2));
      }
      boost::filesystem::create_directory("../data/magnetic_susc/cont/fit_func");
      string Fit_tag= "../data/magnetic_susc/cont/fit_func/"+flavs[iflav]+"_"+contr+".dat";
      Print_To_File({}, {a_to_print_chi, chi_tm_to_print.ave(), chi_tm_to_print.err(), chi_OS_to_print.ave(), chi_OS_to_print.err()},Fit_tag, "", "#a[fm] tm OS, ch2/dof: "+to_string_with_precision(ch2,4));
      
      if(flavs[iflav]== "light") {
	if(contr == "TV") chi_light_TV = D/chiral_cond;
	if(contr == "VT") chi_light_VT = D/chiral_cond;
	if(contr == "VTV") chi_light_VTV = D/chiral_cond;

      }
      else if(flavs[iflav]=="strange") {
	if(contr == "TV") chi_strange_TV = D/chiral_cond;
	if(contr == "VT") chi_strange_VT = D/chiral_cond;
	if(contr == "VTV") chi_strange_VTV = D/chiral_cond;

      }

      

  }
  }

  //print final values

  cout<<"#######"<<endl;
  cout<<"##### light ######"<<endl;
  //cout<<" chiral susc from TV: "<<chi_light_TV.ave()<<" +- "<<chi_light_TV.err()<<"("<<syst_TV*chi_light_TV.ave()<<" ) [GeV^-2]"<<endl;
  //cout<<" chiral susc from VT: "<<chi_light_VT.ave()<<" +- "<<chi_light_VT.err()<<"("<<syst_VT*chi_light_VT.ave()<<" ) [GeV^-2]"<<endl;
  cout<<" chiral susc from VTV: "<<chi_light_VTV.ave()<<" +- "<<chi_light_VTV.err()<<"("<<syst_VTV*chi_light_VTV.ave()<<" ) [GeV^-2]"<<endl;
  cout<<" chiral susc Bali: "<<"u: 2.08(8) GeV^-2 , d: -2.02(9) GeV^-2"<<endl;
  cout<<"##### strange ######"<<endl;
  //cout<<" chiral susc from TV: "<<chi_strange_TV.ave()<<" +- "<<chi_strange_TV.err()<<" [GeV^-2]"<<endl;
  //cout<<" chiral susc from VT: "<<chi_strange_VT.ave()<<" +- "<<chi_strange_VT.err()<<" [GeV^-2]"<<endl;
  cout<<" chiral susc from VTV: "<<chi_strange_VTV.ave()<<" +- "<<chi_strange_VTV.err()<<" [GeV^-2]"<<endl;
  

  

  cout<<"###### PRINTING CH2/NDOF: "<<endl;

  for(int iflav=0; iflav<(signed)flavs.size(); iflav++) {
    for(int icontr=0; icontr<(signed)contribs.size(); icontr++) {
      cout<<flavs[iflav]<<" , "<<contribs[icontr]<<" ch2/dof: "<<ch2_ndof[iflav][icontr]<<endl;
    }
  }
  

  return;

  
}
