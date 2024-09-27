#include "../include/R_ratio.h"


const double MPiPhys=0.135;
const double alpha = 1.0/137.04;
const double e2 = alpha*4.0*M_PI;
const bool UseJack=1;
const int Njacks=50;
const int Nboots=800;
const double qu= 2.0/3.0;
const double qd= -1.0/3.0;
const double qs= qd;
const double qc= qu;
const double fm_to_inv_Gev= 1.0/0.197327;
const int ln2_10=3.32192809489;
const int prec = 200;
const int prec_charm=512;
const double rho_R= 12*M_PI*M_PI;
const double Nc=3;
const double Rpert= Nc*( qu*qu + qd*qd);
const double Rpert_strange = Nc*(qs*qs);
const double Rpert_charm = Nc*(qc*qc);
const double m_mu= 0.10565837; // [ GeV ]
//const Vfloat sigmas({4.2*m_mu});
const Vfloat sigmas({4.2*m_mu, 5.0*m_mu, 6.0*m_mu});
const string SM_TYPE= "GAUSSIAN";
const bool add_pert_corr_light_up_to = 1.0;
const bool add_pert_corr_strange_up_to = 1.0;
const bool add_pert_corr_charm_up_to = 1.0;
const double max_energy = 33.5*m_mu; //33.5*m_mu; //33.5*m_mu // [ GeV ]
const double min_energy = 2.5*m_mu;//2.5*m_mu;  //2.5*m_mu;
const double Eth = 2.0*m_mu;
const double step_size= 0.5*m_mu; //step size
const double m_Jpsi= 3.0969;
const double m_Jpsi_err= 0.001;
const double m_phi= 1.019461;
const double m_phi_err= 0.000016;
const double m_etac= 2.9839;
const double m_etac_err = 0.004;
const double m_etas = 0.68989;
const double m_etas_err= 0.00050;
const string Extrapolation_strange_mode="etas";
const string Extrapolation_charm_mode="etac";
const int pert_corr_strange_on_off = 0;
const int pert_corr_charm_on_off = 0;
const int pert_corr_light_on_off = 0;
const double Lambda_QCD=0.3;
bool Compute_experimental_smeared_R_ratio=false;
bool Compute_free_spec_dens=false;
const int Sim_ord=4;
const bool Use_t_up_to_T_half=true;
const bool R_ratio_verbosity_lev=1;
Vfloat Ergs_GeV_list;
bool SANF_MODE_OFF=true;
bool skip_light=true;
bool skip_strange=true;
bool skip_charm=false;
bool skip_disconnected=true;
Vfloat cov_fake;
int Num_LUSCH_R_ratio=17;
int Nres_R_ratio= 15;
int pts_spline_R_ratio=200;
bool test_mode=false;
bool only_continuum_extrapolation=true;
Vfloat lat_to_print; //fm
using namespace std;


void Get_Ergs_list() {

  Ergs_GeV_list.clear();
  
  double Erg=min_energy;

  while(Erg<=max_energy+1e-10) { Ergs_GeV_list.push_back(Erg); Erg+= step_size;}

  return;


}

void Get_lat_to_print() {

  int Nlat=300;
  double sx= 0.08*1.5/(Nlat-1.0);
  for(int a=0; a < Nlat;a++) lat_to_print.push_back(sx*a);

}


void Get_exp_smeared_R_ratio(const Vfloat& Ergs_GeV_list_exp, double sigma) {


  const auto f_R = [](const double &E, const double &m, const double &s, const double &E0) -> double {
		     
		     if(SM_TYPE == "WIN")  return ((1.0/(1 + exp(-2*(E-(m-s))/s)) - 1.0/(1 + exp(-2*(E-(m+s))/s)))/(2*s));
		     
		     else if (SM_TYPE == "GAUSSIAN") return Get_exact_gauss(E, m, s, E0);
		   
		     else crash("SM_TYPE: "+SM_TYPE+" not yet implemented");
		   
		     return 0.0;
		   
		   };  

  //load File
  Vfloat E = Read_From_File("../dispersive_data/ratio.dati", 0,5);
  Vfloat RE_exp = Read_From_File("../dispersive_data/ratio.dati", 1,5);
  Vfloat RE_exp_stat = Read_From_File("../dispersive_data/ratio.dati", 2,5);
  Vfloat RE_exp_syst = Read_From_File("../dispersive_data/ratio.dati", 3,5);

  //define intervals around resonances where to sum absolute systematic error
  Vfloat Res_MIN({ 0.40, 3.09, 3.74, 9.45, 10.021, 10.352, 10.54}); 
  Vfloat Res_MAX({ 1.06, 3.10, 3.81, 9.47, 10.026, 10.358, 10.62});

  //assign E values to interval  [ Res_MIN[i]: Res_MAX[i] ]
  int Nsteps= E.size();
  vector<int> E_intervals(Nsteps, -1);
  for(int istep=0;istep<Nsteps;istep++) {

    for(int j=0;j<(signed)Res_MIN.size();j++) {
      if( E[istep] >= Res_MIN[j] && E[istep] <= Res_MAX[j]) {
	E_intervals[istep] = j;
	break;
      }
      if( E[istep] < Res_MIN[j]) break;
      if( E[istep] > Res_MAX[Res_MAX.size() -1]) break;
    }

  }


  //start integrating using Simpson Rule 8/3
  
  double dE = 0.00025; //GeV
  Vfloat RE_sm_exp_ave(Ergs_GeV_list_exp.size(), 0.0);
  Vfloat RE_sm_exp_stat(Ergs_GeV_list_exp.size(),0.0);
  Vfloat RE_sm_exp_syst(Ergs_GeV_list_exp.size(),0.0);
  Vfloat RE_sm_exp_err(Ergs_GeV_list_exp.size(), 0.0);

  for(int ierg=0;ierg<(signed)Ergs_GeV_list_exp.size();ierg++) {

    Vfloat DE_corr(Res_MIN.size(),0.0);

    //Simpson integration
    for(int istep=0;istep<Nsteps;istep++) {

      double sm_coeff = f_R(E[istep], Ergs_GeV_list_exp[ierg],sigma, Eth);
      RE_sm_exp_ave[ierg] += w(istep,Sim_ord)*dE*RE_exp[istep]*sm_coeff;
      RE_sm_exp_stat[ierg] += pow(w(istep,Sim_ord)*dE*RE_exp_stat[istep]*sm_coeff,2);
      if( E_intervals[istep] != -1) DE_corr[ E_intervals[istep]] += fabs(w(istep,Sim_ord)*dE*RE_exp_syst[istep]*sm_coeff);
      else RE_sm_exp_syst[ierg] += pow( w(istep,Sim_ord)*dE*RE_exp_syst[istep]*sm_coeff,2);
    }

    //add sum of absolute systematic error around resonances
    for(int ires=0;ires< (signed)Res_MIN.size();ires++) RE_sm_exp_syst[ierg] += pow( DE_corr[ires],2);

       
    RE_sm_exp_err[ierg] = RE_sm_exp_stat[ierg] + RE_sm_exp_syst[ierg];

    //take sqrt of stat^2 syst^2 and tot^2
    RE_sm_exp_err[ierg] = sqrt( RE_sm_exp_err[ierg]);
    RE_sm_exp_stat[ierg] = sqrt( RE_sm_exp_stat[ierg]);
    RE_sm_exp_syst[ierg] = sqrt( RE_sm_exp_syst[ierg]);

  }

     
     
  cout<<"Smeared experimental R(e+e- -> hadrons) with sigma: "<<sigma<<" computed!"<<endl;
  //Print To File
  Print_To_File({}, {Ergs_GeV_list_exp ,RE_sm_exp_ave, RE_sm_exp_stat, RE_sm_exp_syst, RE_sm_exp_err}, "../dispersive_data/smeared_data/smeared_R_sigma_"+to_string_with_precision(sigma,3)+".dat", "", "#E[GeV]  R[E]    stat. syst.  total ");
}
   


void R_ratio_analysis() {


  Get_Ergs_list();


  if(test_mode && ( !skip_charm || !skip_strange || ! skip_disconnected)) crash("test mode is safe only if skip_light=1 and skip_charm=skip_strange=skip_disconnected=0");


  if(only_continuum_extrapolation) { R_ratio_cont_extrapolation(); exit(-1);}



  if(skip_light==true) { //do not perform GS-analysis of systematics
    Num_LUSCH_R_ratio=2;
    Nres_R_ratio= 1;
    pts_spline_R_ratio=10;
  }

  Vfloat Ergs_GeV_list_exp;
  double Dx_e= 0.01; //10 MeV                                                                                                                                                                             
  int Nps_exp= (int)((5.0 - Eth)/Dx_e); //up to 5 GeV                                                                                                                                                     
  for(int ip=0;ip<Nps_exp;ip++) Ergs_GeV_list_exp.push_back(  Eth+ ip*((5.0- Eth)/((double)Nps_exp)));
  for( auto &sigma: sigmas)
    if(Compute_experimental_smeared_R_ratio) Get_exp_smeared_R_ratio(Ergs_GeV_list_exp, sigma);



  
  Vfloat betas({ 0.0, 1.0, 1.99, 2.99, 3.99, 0.0, 1.0, 1.99});
  Vfloat Emax_list({ 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0});
  vector<bool> Is_Emax_Finite({1,1,1,1,1, 0,0,0});
  int N= betas.size();

  int rank,size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);


  if(N%size != 0) crash("MPI called with -np= "+to_string(size)+". np does not divide vector size N="+to_string(N));

  cout<<"################# DETERMINATION OF THE SMEARED R-ratio #################"<<endl;
  cout<<"Rank: "<<rank<<" pid: "<<getpid()<<" core id: "<<"("<<sched_getcpu()<<")"<<endl;
  cout<<"SMEARING_FUNCTION: "<<SM_TYPE<<endl;
  cout<<"INVERSE LAPLACE RECONSTRUCTION CALLED FOR:"<<endl;
  for(int i=rank*N/size;i<(rank+1)*N/size;i++) {
    string alpha_Emax_Tag= "{"+to_string_with_precision(betas[i],2)+","+((Is_Emax_Finite[i]==0)?"inf":to_string_with_precision(Emax_list[i],1))+"}";
    cout<<"{alpha,Emax} = "<<alpha_Emax_Tag<<endl;
  }
  cout<<"##########################################"<<endl;




  
  for(int i=rank*N/size;i<(rank+1)*N/size;i++) {Compute_R_ratio(Is_Emax_Finite[i], Emax_list[i], betas[i]); Compute_free_spec_dens=false;}

  return;

}



void Compute_R_ratio(bool Is_Emax_Finite, double Emax, double beta) {

  string Tag_reco_type="Beta_"+to_string_with_precision(beta,2);
  Tag_reco_type+="_Emax_"+((Is_Emax_Finite==0)?"inf":to_string_with_precision(Emax,1));

  string alpha_Emax_tag= "{"+to_string_with_precision(beta,2)+","+((Is_Emax_Finite==0)?"inf":to_string_with_precision(Emax,1))+"}";

  cout<<"STARTING COMPUTATION OF: {alpha,Emax} : "<<alpha_Emax_tag<<endl;


  //COMPUTE FREE SPECTRAL DENSITY IS NEEDED
  if(Compute_free_spec_dens) {
    cout<<"Computing free-theory spectral density"<<endl<<flush;
    Vfloat ams({0.00072,0.00060,0.00054,0.019,0.021,0.01600,0.01800,0.014,0.015,0.21000,0.23000,0.25000,0.17500,0.19500,0.21500,0.165,0.175});
    Get_spec_dens_free(ams,"R_ratio");
    cout<<"done!"<<endl<<flush;
  }

  
  cout<<"Creating output directories...";


  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int _hostname_len;
  char _hostname[MPI_MAX_PROCESSOR_NAME];
  MPI_Get_processor_name(_hostname, &_hostname_len);
  
  
  
  //create output directories
  boost::filesystem::create_directory("../dispersive_data");
  boost::filesystem::create_directory("../dispersive_data/smeared_data");
  boost::filesystem::create_directory("../data/R_ratio");
  boost::filesystem::create_directory("../data/R_ratio/"+Tag_reco_type);
  boost::filesystem::create_directory("../data/R_ratio/"+Tag_reco_type+"/light");
  boost::filesystem::create_directory("../data/R_ratio/"+Tag_reco_type+"/light/jackknife");
  boost::filesystem::create_directory("../data/R_ratio/"+Tag_reco_type+"/light/jackknife/tm");
  boost::filesystem::create_directory("../data/R_ratio/"+Tag_reco_type+"/light/jackknife/OS");
  boost::filesystem::create_directory("../data/R_ratio/"+Tag_reco_type+"/disco");
  boost::filesystem::create_directory("../data/R_ratio/"+Tag_reco_type+"/disco/jackknife");
  boost::filesystem::create_directory("../data/R_ratio/"+Tag_reco_type+"/strange");
  boost::filesystem::create_directory("../data/R_ratio/"+Tag_reco_type+"/strange/jackknife");
  boost::filesystem::create_directory("../data/R_ratio/"+Tag_reco_type+"/strange/jackknife/tm");
  boost::filesystem::create_directory("../data/R_ratio/"+Tag_reco_type+"/strange/jackknife/OS");
  boost::filesystem::create_directory("../data/R_ratio/"+Tag_reco_type+"/charm");
  boost::filesystem::create_directory("../data/R_ratio/"+Tag_reco_type+"/charm/jackknife");
  boost::filesystem::create_directory("../data/R_ratio/"+Tag_reco_type+"/charm/jackknife/tm");
  boost::filesystem::create_directory("../data/R_ratio/"+Tag_reco_type+"/charm/jackknife/OS");
  boost::filesystem::create_directory("../data/R_ratio/"+Tag_reco_type+"/corr");
  boost::filesystem::create_directory("../data/R_ratio/"+Tag_reco_type+"/corr/light");
  boost::filesystem::create_directory("../data/R_ratio/"+Tag_reco_type+"/corr/disco");
  boost::filesystem::create_directory("../data/R_ratio/"+Tag_reco_type+"/corr/strange");
  boost::filesystem::create_directory("../data/R_ratio/"+Tag_reco_type+"/corr/charm");
  boost::filesystem::create_directory("../data/R_ratio/"+Tag_reco_type+"/total");
  boost::filesystem::create_directory("../data/R_ratio/"+Tag_reco_type+"/covariance");
  boost::filesystem::create_directory("../data/R_ratio/"+Tag_reco_type+"/covariance/light");
  boost::filesystem::create_directory("../data/R_ratio/"+Tag_reco_type+"/covariance/disco");
  boost::filesystem::create_directory("../data/R_ratio/"+Tag_reco_type+"/covariance/strange");
  boost::filesystem::create_directory("../data/R_ratio/"+Tag_reco_type+"/covariance/charm");

  cout<<"done!"<<endl;


  

   //Init LL_functions;
  //find first  zeros of the Lusher functions
  Vfloat Lusch_zeroes;
  Zeta_function_zeroes(Num_LUSCH_R_ratio, Lusch_zeroes);
  

  //############################################INTERPOLATE PHI FUNCTION AND DERIVATIVES#############################
  cout<<"Computing Luscher-zeros"<<flush;

  VVfloat phi_data, phi_der_data;
  Vfloat sx_int;
  Vfloat sx_der, dx_der;
  Vfloat Dz;
  
  for(int L_zero=0;L_zero<Nres_R_ratio+1;L_zero++) {
    cout<<"."<<flush;
    double sx, dx;
    //interpolating between the Luscher_zero[L_zero-1] and Luscher_zero[L_zero];
    if(L_zero==0) { sx_int.push_back(0.0); sx=0.0;}
    else {sx=Lusch_zeroes[L_zero-1];  sx_int.push_back(sx);}
    dx= Lusch_zeroes[L_zero];
    phi_data.resize(L_zero+1);
    phi_der_data.resize(L_zero+1);
    phi_data[L_zero].push_back(L_zero==0?0.0:-M_PI/2.0);
    //divide interval into thousand points;
    double dz = (dx-sx)/pts_spline_R_ratio;
    Dz.push_back(dz);


    for(int istep=1;istep<=pts_spline_R_ratio-1;istep++) { double pt= sx+dz*istep; phi_data[L_zero].push_back( phi(sqrt(pt)));}

    phi_data[L_zero].push_back(M_PI/2.0);
    double sx_der_loc =  phi_der_for_back(sqrt(sx)+1e-14, 1);
    double dx_der_loc =  phi_der_for_back(sqrt(dx)-1e-14, -1);
    sx_der.push_back(sx_der_loc);
    dx_der.push_back(dx_der_loc);

    phi_der_data[L_zero].push_back(sx_der_loc);
    for(int istep=1;istep<=pts_spline_R_ratio-1;istep++) { double pt= sx+dz*istep; phi_der_data[L_zero].push_back( phi_der(sqrt(pt)));}
    phi_der_data[L_zero].push_back(dx_der_loc);
    
  }
  cout<<"done!"<<endl;



  //###########################################END INTERPOLATION PHI FUNCTION AND DERIVATIVES################################
 
   

  LL_functions LL(phi_data,phi_der_data,sx_der, dx_der, sx_int, Dz, Nres_R_ratio, Lusch_zeroes);

  //LOAD DATA
  GaussianMersenne GM(981832);

  //light
  data_t  V_light_1, V_light_OS_1;
  data_t disco_light, disco_strange, disco_charm, disco_light_strange, disco_light_charm, disco_strange_charm;

  //strange and charm
  
  //L
  data_t  V_strange_1_L, V_strange_OS_1_L;
  data_t  V_charm_1_L, V_charm_2_L, V_charm_3_L, V_charm_OS_1_L, V_charm_OS_2_L, V_charm_OS_3_L;
  data_t  pt2_etaC_L, pt2_etaC_OS_L, pt2_etaS_L, pt2_etaS_OS_L;
  //M
  data_t  V_strange_1_M,  V_strange_OS_1_M;
  data_t  V_charm_1_M, V_charm_2_M, V_charm_3_M, V_charm_OS_1_M, V_charm_OS_2_M, V_charm_OS_3_M;
  data_t  pt2_etaC_M, pt2_etaC_OS_M, pt2_etaS_M, pt2_etaS_OS_M;
  //H
  data_t  V_charm_1_H, V_charm_2_H, V_charm_3_H, V_charm_OS_1_H, V_charm_OS_2_H, V_charm_OS_3_H;
  data_t  pt2_etaC_H, pt2_etaC_OS_H;




  //Read data

  //Custom sorting for V_light to account for the two replica r0 and r1
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


  auto Sort_light_confs_h5 = [](string A, string B) {


			     //return A<B;
			     
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


  // lambda function to be used as a smearing func.

  const auto f = [](const PrecFloat &E, const PrecFloat &m, const PrecFloat &s, const PrecFloat &E0, int Nj) -> PrecFloat {
		   
		   if(SM_TYPE == "WIN")  return (1.0/(1 + exp(-2*(E-(m-s))/s)) - 1.0/(1 + exp(-2*(E-(m+s))/s)))/(2*s*sqr(E));
		   
		   else if (SM_TYPE == "GAUSSIAN") return Get_exact_gauss(E, m, s, E0)/sqr(E);
		   
		   else crash("SM_TYPE: "+SM_TYPE+" not yet implemented");
		   
		   return PrecFloat(0);
		   
		 };

  

  cout<<"Loading data...";



  //light 
  //#################################END CUSTOM SORTING#################
  V_light_1.Read("../R_ratio_data/light", "mes_contr_2pts_ll_1", "VKVK", Sort_light_confs_h5);
  V_light_OS_1.Read("../R_ratio_data/light", "mes_contr_2pts_ll_2", "VKVK", Sort_light_confs_h5);
 

  //disco_light
  disco_light.Read("../R_ratio_data/disconnected/light_light", "disco", "", Sort_light_confs);
  //disco_strange
  disco_strange.Read("../R_ratio_data/disconnected/strange_strange", "disco", "", Sort_light_confs);
  //disco_charm
  disco_charm.Read("../R_ratio_data/disconnected/charm_charm", "disco", "", Sort_light_confs);
  //disco light-strange
  disco_light_strange.Read("../R_ratio_data/disconnected/light_strange", "disco", "", Sort_light_confs);
  //disco light-charm
  disco_light_charm.Read("../R_ratio_data/disconnected/light_charm", "disco", "", Sort_light_confs);
  //disco strange-charm
  disco_strange_charm.Read("../R_ratio_data/disconnected/strange_charm", "disco", "", Sort_light_confs);
  //#################################END READING LIGHT########################

  //strange
  //L
  V_strange_1_L.Read("../R_ratio_data/strange/light", "mes_contr_2pts_ll_1", "V1V1", Sort_light_confs);
  pt2_etaS_L.Read("../R_ratio_data/strange/light", "mes_contr_2pts_ll_1", "P5P5", Sort_light_confs);
  V_strange_OS_1_L.Read("../R_ratio_data/strange/light", "mes_contr_2pts_ll_2", "V1V1", Sort_light_confs); 
  pt2_etaS_OS_L.Read("../R_ratio_data/strange/light", "mes_contr_2pts_ll_2", "P5P5", Sort_light_confs);
  //M
  V_strange_1_M.Read("../R_ratio_data/strange/heavy", "mes_contr_2pts_ll_1", "V1V1", Sort_light_confs); 
  pt2_etaS_M.Read("../R_ratio_data/strange/heavy", "mes_contr_2pts_ll_1", "P5P5", Sort_light_confs);
  V_strange_OS_1_M.Read("../R_ratio_data/strange/heavy", "mes_contr_2pts_ll_2", "V1V1", Sort_light_confs); 
  pt2_etaS_OS_M.Read("../R_ratio_data/strange/heavy", "mes_contr_2pts_ll_2", "P5P5", Sort_light_confs);
  //#################################END READING STRANGE########################



  //charm
  //L
  V_charm_1_L.Read("../R_ratio_data/charm/light", "mes_contr_2pts_ll_1", "V1V1",  Sort_light_confs); 
  V_charm_2_L.Read("../R_ratio_data/charm/light", "mes_contr_2pts_ll_1", "V2V2",  Sort_light_confs);
  V_charm_3_L.Read("../R_ratio_data/charm/light", "mes_contr_2pts_ll_1", "V3V3",  Sort_light_confs);
  pt2_etaC_L.Read( "../R_ratio_data/charm/light", "mes_contr_2pts_ll_1", "P5P5", Sort_light_confs); 
  V_charm_OS_1_L.Read("../R_ratio_data/charm/light", "mes_contr_2pts_ll_2", "V1V1",  Sort_light_confs); 
  V_charm_OS_2_L.Read("../R_ratio_data/charm/light", "mes_contr_2pts_ll_2", "V2V2",  Sort_light_confs);
  V_charm_OS_3_L.Read("../R_ratio_data/charm/light", "mes_contr_2pts_ll_2", "V3V3",  Sort_light_confs);
  pt2_etaC_OS_L.Read( "../R_ratio_data/charm/light", "mes_contr_2pts_ll_2", "P5P5", Sort_light_confs); 
  //M
  V_charm_1_M.Read("../R_ratio_data/charm/medium", "mes_contr_2pts_ll_1", "V1V1",  Sort_light_confs); 
  V_charm_2_M.Read("../R_ratio_data/charm/medium", "mes_contr_2pts_ll_1", "V2V2",  Sort_light_confs);
  V_charm_3_M.Read("../R_ratio_data/charm/medium", "mes_contr_2pts_ll_1", "V3V3",  Sort_light_confs);
  pt2_etaC_M.Read( "../R_ratio_data/charm/medium", "mes_contr_2pts_ll_1", "P5P5", Sort_light_confs); 
  V_charm_OS_1_M.Read("../R_ratio_data/charm/medium", "mes_contr_2pts_ll_2", "V1V1",  Sort_light_confs); 
  V_charm_OS_2_M.Read("../R_ratio_data/charm/medium", "mes_contr_2pts_ll_2", "V2V2",  Sort_light_confs);
  V_charm_OS_3_M.Read("../R_ratio_data/charm/medium", "mes_contr_2pts_ll_2", "V3V3",  Sort_light_confs);
  pt2_etaC_OS_M.Read( "../R_ratio_data/charm/medium", "mes_contr_2pts_ll_2", "P5P5", Sort_light_confs); 
  //H
  V_charm_1_H.Read("../R_ratio_data/charm/heavy", "mes_contr_2pts_ll_1", "V1V1",  Sort_light_confs);
  V_charm_2_H.Read("../R_ratio_data/charm/heavy", "mes_contr_2pts_ll_1", "V2V2",  Sort_light_confs);
  V_charm_3_H.Read("../R_ratio_data/charm/heavy", "mes_contr_2pts_ll_1", "V3V3",  Sort_light_confs);
  pt2_etaC_H.Read( "../R_ratio_data/charm/heavy", "mes_contr_2pts_ll_1", "P5P5", Sort_light_confs); 
  V_charm_OS_1_H.Read("../R_ratio_data/charm/heavy", "mes_contr_2pts_ll_2", "V1V1",  Sort_light_confs); 
  V_charm_OS_2_H.Read("../R_ratio_data/charm/heavy", "mes_contr_2pts_ll_2", "V2V2",  Sort_light_confs);
  V_charm_OS_3_H.Read("../R_ratio_data/charm/heavy", "mes_contr_2pts_ll_2", "V3V3",  Sort_light_confs);
  pt2_etaC_OS_H.Read( "../R_ratio_data/charm/heavy", "mes_contr_2pts_ll_2", "P5P5", Sort_light_confs); 

  cout<<"done!"<<endl;


  //###########################################################################################


  //Ensemble Tags

  //light Tags
  vector<string> Light_Ens_Tags;
  vector<string> Strange_Ens_Tags;
  vector<string> Charm_Ens_Tags;
  vector<string> Disco_Ens_Tags;


  //############################################################################################
  //generate fake jack_distr for lattice spacing a_A a_B, a_C, a_D and RENORMALIZATION CONSTANT
  LatticeInfo a_info;
  distr_t a_A(UseJack), a_B(UseJack), a_C(UseJack), a_D(UseJack);
  distr_t ZV_A(UseJack), ZV_B(UseJack), ZV_C(UseJack), ZV_D(UseJack);
  distr_t ZA_A(UseJack), ZA_B(UseJack), ZA_C(UseJack), ZA_D(UseJack);
  double a_A_ave, a_A_err, a_B_ave, a_B_err, a_C_ave, a_C_err, a_D_ave, a_D_err;
  double ZV_A_ave, ZV_A_err, ZV_B_ave, ZV_B_err, ZV_C_ave, ZV_C_err, ZV_D_ave, ZV_D_err;
  double ZA_A_ave, ZA_A_err, ZA_B_ave, ZA_B_err, ZA_C_ave, ZA_C_err, ZA_D_ave, ZA_D_err;
  a_info.LatInfo_new_ens("cA211a.53.24");
  a_A_ave= a_info.a_from_afp;
  a_A_err= a_info.a_from_afp_err;
  ZA_A_ave = a_info.Za_WI_strange;
  ZA_A_err = a_info.Za_WI_strange_err;
  ZV_A_ave = a_info.Zv_WI_strange;
  ZV_A_err = a_info.Zv_WI_strange_err;
  a_info.LatInfo_new_ens("cB211b.072.64");
  a_B_ave= a_info.a_from_afp;
  a_B_err= a_info.a_from_afp_err;
  ZA_B_ave = a_info.Za_WI_strange;
  ZA_B_err = a_info.Za_WI_strange_err;
  ZV_B_ave = a_info.Zv_WI_strange;
  ZV_B_err = a_info.Zv_WI_strange_err;
  a_info.LatInfo_new_ens("cC211a.06.80");
  a_C_ave= a_info.a_from_afp;
  a_C_err= a_info.a_from_afp_err;
  ZA_C_ave = a_info.Za_WI_strange;
  ZA_C_err = a_info.Za_WI_strange_err;
  ZV_C_ave = a_info.Zv_WI_strange;
  ZV_C_err = a_info.Zv_WI_strange_err;
  a_info.LatInfo_new_ens("cD211a.054.96");
  a_D_ave= a_info.a_from_afp;
  a_D_err= a_info.a_from_afp_err;
  ZA_D_ave = a_info.Za_WI_strange;
  ZA_D_err = a_info.Za_WI_strange_err;
  ZV_D_ave = a_info.Zv_WI_strange;
  ZV_D_err = a_info.Zv_WI_strange_err;
  
  if(UseJack)  { for(int ijack=0;ijack<Njacks;ijack++) {
      a_A.distr.push_back( fm_to_inv_Gev*( a_A_ave + GM()*a_A_err*(1.0/sqrt(Njacks-1.0))));
      a_B.distr.push_back( fm_to_inv_Gev*( a_B_ave + GM()*a_B_err*(1.0/sqrt(Njacks-1.0))));
      a_C.distr.push_back( fm_to_inv_Gev*( a_C_ave + GM()*a_C_err*(1.0/sqrt(Njacks-1.0))));
      a_D.distr.push_back( fm_to_inv_Gev*( a_D_ave + GM()*a_D_err*(1.0/sqrt(Njacks-1.0))));
      ZA_A.distr.push_back(  ZA_A_ave + GM()*ZA_A_err*(1.0/sqrt(Njacks -1.0)));
      ZV_A.distr.push_back(  ZV_A_ave + GM()*ZV_A_err*(1.0/sqrt(Njacks -1.0)));
      ZA_B.distr.push_back(  ZA_B_ave + GM()*ZA_B_err*(1.0/sqrt(Njacks -1.0)));
      ZV_B.distr.push_back(  ZV_B_ave + GM()*ZV_B_err*(1.0/sqrt(Njacks -1.0)));
      ZA_C.distr.push_back(  ZA_C_ave + GM()*ZA_C_err*(1.0/sqrt(Njacks -1.0)));
      ZV_C.distr.push_back(  ZV_C_ave + GM()*ZV_C_err*(1.0/sqrt(Njacks -1.0)));
      ZA_D.distr.push_back(  ZA_D_ave + GM()*ZA_D_err*(1.0/sqrt(Njacks -1.0)));
      ZV_D.distr.push_back(  ZV_D_ave + GM()*ZV_D_err*(1.0/sqrt(Njacks -1.0)));
      
    }
  }
  else {
    for (int iboot=0; iboot<Nboots;iboot++) {
      a_A.distr.push_back( fm_to_inv_Gev*( a_A_ave + GM()*a_A_err));
      a_B.distr.push_back( fm_to_inv_Gev*( a_B_ave + GM()*a_B_err));
      a_C.distr.push_back( fm_to_inv_Gev*( a_C_ave + GM()*a_C_err));
      a_D.distr.push_back( fm_to_inv_Gev*( a_D_ave + GM()*a_D_err));
      ZA_A.distr.push_back(  ZA_A_ave + GM()*ZA_A_err);
      ZV_A.distr.push_back(  ZV_A_ave + GM()*ZV_A_err);
      ZA_B.distr.push_back(  ZA_B_ave + GM()*ZA_B_err);
      ZV_B.distr.push_back(  ZV_B_ave + GM()*ZV_B_err);
      ZA_C.distr.push_back(  ZA_C_ave + GM()*ZA_C_err);
      ZV_C.distr.push_back(  ZV_C_ave + GM()*ZV_C_err);
      ZA_D.distr.push_back(  ZA_D_ave + GM()*ZA_D_err);
      ZV_D.distr.push_back(  ZV_D_ave + GM()*ZV_D_err);
      
    }
  }


  //############################################################################################


  //############################################################################################
  //generate fake jack_distr for M_etas and M_phi phys
  distr_t Mphi_phys_distr(UseJack), Metas_phys_distr(UseJack), MJpsi_phys_distr(UseJack), Metac_phys_distr(UseJack);
  for(int ijack=0; ijack<Njacks;ijack++) {
    Mphi_phys_distr.distr.push_back( m_phi + GM()*m_phi_err*(1.0/sqrt(Njacks -1.0)));
    Metas_phys_distr.distr.push_back( m_etas + GM()*m_etas_err*(1.0/sqrt(Njacks -1.0)));
    MJpsi_phys_distr.distr.push_back( m_Jpsi + GM()*m_Jpsi_err*(1.0/sqrt(Njacks -1.0)));
    Metac_phys_distr.distr.push_back( m_etac + GM()*m_etac_err*(1.0/sqrt(Njacks -1.0)));
  }



 

  int Nens_light =  V_light_1.size;
  int Nens_disco = disco_light.size;
  int Nens_strange = V_strange_1_L.size;
  int Nens_charm = V_charm_1_L.size;


  //vector of distr_t_list to store the R(E) results for all ensemble (different vectors for different contributions)

  //light
  vector<vector<distr_t_list>> RE_light_TANT_OS(sigmas.size());
  vector<vector<distr_t_list>> RE_light_TANT_tm(sigmas.size());
  vector<vector<distr_t_list>> RE_light_SANF_OS(sigmas.size());
  vector<vector<distr_t_list>> RE_light_SANF_tm(sigmas.size());

 

  //strange
  vector<vector<distr_t_list>> RE_strange_TANT_OS(sigmas.size());
  vector<vector<distr_t_list>> RE_strange_TANT_tm(sigmas.size());
  vector<vector<distr_t_list>> RE_strange_SANF_OS(sigmas.size());
  vector<vector<distr_t_list>> RE_strange_SANF_tm(sigmas.size());

  //charm
  vector<vector<distr_t_list>> RE_charm_TANT_OS(sigmas.size());
  vector<vector<distr_t_list>> RE_charm_TANT_tm(sigmas.size());
  vector<vector<distr_t_list>> RE_charm_SANF_OS(sigmas.size());
  vector<vector<distr_t_list>> RE_charm_SANF_tm(sigmas.size());

  //disconnected
  vector<vector<distr_t_list>> RE_disco_TANT(sigmas.size());
  vector<vector<distr_t_list>> RE_disco_SANF(sigmas.size());

  for(int isg=0;isg<(signed)sigmas.size();isg++) {

    //light 
    RE_light_TANT_OS[isg].resize(Nens_light);
    RE_light_TANT_tm[isg].resize(Nens_light);
    RE_light_SANF_OS[isg].resize(Nens_light);
    RE_light_SANF_tm[isg].resize(Nens_light);

    //strange
    RE_strange_TANT_OS[isg].resize(Nens_strange);
    RE_strange_TANT_tm[isg].resize(Nens_strange);
    RE_strange_SANF_OS[isg].resize(Nens_strange);
    RE_strange_SANF_tm[isg].resize(Nens_strange);

    //charm
    RE_charm_TANT_OS[isg].resize(Nens_charm);
    RE_charm_TANT_tm[isg].resize(Nens_charm);
    RE_charm_SANF_OS[isg].resize(Nens_charm);
    RE_charm_SANF_tm[isg].resize(Nens_charm);

    //disconnected
    RE_disco_TANT[isg].resize(Nens_light);
    RE_disco_SANF[isg].resize(Nens_light);     
  }

  Light_Ens_Tags.resize(Nens_light);
  Strange_Ens_Tags.resize(Nens_strange);
  Charm_Ens_Tags.resize(Nens_charm);










  //#############################################################################################
  //#############################################################################################
  //#############################################################################################
  //#############################################################################################
  
  //################## CALCULATION OF THE CHARM CONTRIBUTION ##################################

  
  //#############################################################################################
  //#############################################################################################
  //#############################################################################################
  //#############################################################################################

  //get charm correlators only on B,C,D physical point ensemble

  if(!skip_charm) {

    cout<<"STARTING COMPUTATION OF CHARM CONTRIBUTION:"<<endl;
    cout<<"PRECISION: "<<prec_charm/ln2_10<<" digits."<<endl;
  
  for(int i_ens=0;i_ens<Nens_charm;i_ens++) {

      
    Charm_Ens_Tags[i_ens] = V_charm_1_L.Tag[i_ens];


    string FLAV= "charm";
    string charm_tag = "charm";
    
    CorrAnalysis Corr(UseJack, Njacks,Nboots);
    CorrAnalysis Corr_block_1(0,V_charm_1_L.Nconfs[i_ens], Nboots, i_ens);
    Corr_block_1.Nt= V_charm_1_L.nrows[i_ens];
    Corr.Nt = V_charm_1_L.nrows[i_ens];
    int T = Corr.Nt;

    cout<<"Analyzing Ensemble: "<<V_charm_1_L.Tag[i_ens]<<endl;

    //get lattice spacing
    distr_t a_distr(UseJack);
    distr_t Zv(UseJack), Za(UseJack);
    if(V_charm_1_L.Tag[i_ens].substr(1,1)=="A") {a_distr=a_A; Zv = ZV_A; Za= ZA_A;}
    else if(V_charm_1_L.Tag[i_ens].substr(1,1)=="B") {a_distr=a_B; Zv = ZV_B; Za= ZA_B;}
    else if(V_charm_1_L.Tag[i_ens].substr(1,1)=="C") {a_distr=a_C; Zv = ZV_C; Za= ZA_C;}
    else if(V_charm_1_L.Tag[i_ens].substr(1,1)=="D") {a_distr=a_D; Zv = ZV_D; Za= ZA_D;}
    else crash("lattice spacing distribution for Ens: "+V_charm_1_L.Tag[i_ens]+" not found");



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

    Corr.Tmin = Tmin_etaC; Corr.Tmax= Tmax_etaC;

    int Tmin_VV, Tmax_VV;
 
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

    

      
  
    LatticeInfo L_info;
    L_info.LatInfo_new_ens(V_charm_1_L.Tag[i_ens]);
 
    //L
    distr_t_list   V_charm_L_distr, V_charm_L_bin_distr;
    distr_t_list   V_charm_OS_L_distr, V_charm_OS_L_bin_distr;
    distr_t  Metac_L_distr, Metac_OS_L_distr, MV_L_distr, MV_OS_L_distr;
    distr_t overlap_V_L_distr, overlap_OS_V_L_distr;
    //M
    distr_t_list   V_charm_M_distr, V_charm_M_bin_distr;
    distr_t_list   V_charm_OS_M_distr, V_charm_OS_M_bin_distr;
    distr_t Metac_M_distr, Metac_OS_M_distr, MV_M_distr, MV_OS_M_distr;
    distr_t overlap_V_M_distr, overlap_OS_V_M_distr;
    //H
    distr_t_list   V_charm_H_distr, V_charm_H_bin_distr;
    distr_t_list   V_charm_OS_H_distr, V_charm_OS_H_bin_distr;
    distr_t Metac_H_distr, Metac_OS_H_distr, MV_H_distr, MV_OS_H_distr;
    distr_t overlap_V_H_distr, overlap_OS_V_H_distr;


     
  

    //vector charm
    //L
    V_charm_L_distr = (1.0/3.0)*Corr.corr_t(summ_master(V_charm_1_L.col(0)[i_ens], V_charm_2_L.col(0)[i_ens], V_charm_3_L.col(0)[i_ens]),"");
    V_charm_OS_L_distr = (1.0/3.0)*Corr.corr_t(summ_master(V_charm_OS_1_L.col(0)[i_ens], V_charm_OS_2_L.col(0)[i_ens], V_charm_OS_3_L.col(0)[i_ens]),"");
    V_charm_L_bin_distr = (1.0/3.0)*Corr_block_1.corr_t(summ_master(V_charm_1_L.col(0)[i_ens], V_charm_2_L.col(0)[i_ens], V_charm_3_L.col(0)[i_ens]),"");
    V_charm_OS_L_bin_distr = (1.0/3.0)*Corr_block_1.corr_t(summ_master(V_charm_OS_1_L.col(0)[i_ens], V_charm_OS_2_L.col(0)[i_ens], V_charm_OS_3_L.col(0)[i_ens]),"");
    //M
    V_charm_M_distr = (1.0/3.0)*Corr.corr_t(summ_master(V_charm_1_M.col(0)[i_ens], V_charm_2_M.col(0)[i_ens], V_charm_3_M.col(0)[i_ens]),"");
    V_charm_OS_M_distr = (1.0/3.0)*Corr.corr_t(summ_master(V_charm_OS_1_M.col(0)[i_ens], V_charm_OS_2_M.col(0)[i_ens], V_charm_OS_3_M.col(0)[i_ens]),"");
    V_charm_M_bin_distr = (1.0/3.0)*Corr_block_1.corr_t(summ_master(V_charm_1_M.col(0)[i_ens], V_charm_2_M.col(0)[i_ens], V_charm_3_M.col(0)[i_ens]),"");
    V_charm_OS_M_bin_distr = (1.0/3.0)*Corr_block_1.corr_t(summ_master(V_charm_OS_1_M.col(0)[i_ens], V_charm_OS_2_M.col(0)[i_ens], V_charm_OS_3_M.col(0)[i_ens]),"");
    //H
    V_charm_H_distr = (1.0/3.0)*Corr.corr_t(summ_master(V_charm_1_H.col(0)[i_ens], V_charm_2_H.col(0)[i_ens], V_charm_3_H.col(0)[i_ens]),"");
    V_charm_OS_H_distr = (1.0/3.0)*Corr.corr_t(summ_master(V_charm_OS_1_H.col(0)[i_ens], V_charm_OS_2_H.col(0)[i_ens], V_charm_OS_3_H.col(0)[i_ens]),"");
    V_charm_H_bin_distr = (1.0/3.0)*Corr_block_1.corr_t(summ_master(V_charm_1_H.col(0)[i_ens], V_charm_2_H.col(0)[i_ens], V_charm_3_H.col(0)[i_ens]),"");
    V_charm_OS_H_bin_distr = (1.0/3.0)*Corr_block_1.corr_t(summ_master(V_charm_OS_1_H.col(0)[i_ens], V_charm_OS_2_H.col(0)[i_ens], V_charm_OS_3_H.col(0)[i_ens]),"");


    //print covariance matrix

    Vfloat cov_tm_L, cov_OS_L, cov_tm_M, cov_OS_M, cov_tm_H, cov_OS_H, TT, RR;
    Vfloat corr_tm_L, corr_OS_L, corr_tm_M, corr_OS_M, corr_tm_H, corr_OS_H;
    
    for(int tt=0;tt<Corr.Nt;tt++)
      for(int rr=0;rr<Corr.Nt;rr++) {
	TT.push_back(tt);
	RR.push_back(rr);
	//L
	cov_tm_L.push_back( (V_charm_L_bin_distr.distr_list[tt]%V_charm_L_bin_distr.distr_list[rr]));
	cov_OS_L.push_back( (V_charm_OS_L_bin_distr.distr_list[tt]%V_charm_OS_L_bin_distr.distr_list[rr]));
	corr_tm_L.push_back( (V_charm_L_bin_distr.distr_list[tt]%V_charm_L_bin_distr.distr_list[rr])/(V_charm_L_bin_distr.err(tt)*V_charm_L_bin_distr.err(rr)));
	corr_OS_L.push_back( (V_charm_OS_L_bin_distr.distr_list[tt]%V_charm_OS_L_bin_distr.distr_list[rr])/( V_charm_OS_L_bin_distr.err(tt)*V_charm_OS_L_bin_distr.err(rr)));

	//M
	cov_tm_M.push_back( (V_charm_M_bin_distr.distr_list[tt]%V_charm_M_bin_distr.distr_list[rr]));
	cov_OS_M.push_back( (V_charm_OS_M_bin_distr.distr_list[tt]%V_charm_OS_M_bin_distr.distr_list[rr]));
	corr_tm_M.push_back( (V_charm_M_bin_distr.distr_list[tt]%V_charm_M_bin_distr.distr_list[rr])/(V_charm_M_bin_distr.err(tt)*V_charm_M_bin_distr.err(rr)));
	corr_OS_M.push_back( (V_charm_OS_M_bin_distr.distr_list[tt]%V_charm_OS_M_bin_distr.distr_list[rr])/( V_charm_OS_M_bin_distr.err(tt)*V_charm_OS_M_bin_distr.err(rr)));


	//H
	cov_tm_H.push_back( (V_charm_H_bin_distr.distr_list[tt]%V_charm_H_bin_distr.distr_list[rr]));
	cov_OS_H.push_back( (V_charm_OS_H_bin_distr.distr_list[tt]%V_charm_OS_H_bin_distr.distr_list[rr]));
	corr_tm_H.push_back( (V_charm_H_bin_distr.distr_list[tt]%V_charm_H_bin_distr.distr_list[rr])/(V_charm_H_bin_distr.err(tt)*V_charm_H_bin_distr.err(rr)));
	corr_OS_H.push_back( (V_charm_OS_H_bin_distr.distr_list[tt]%V_charm_OS_H_bin_distr.distr_list[rr])/( V_charm_OS_H_bin_distr.err(tt)*V_charm_OS_H_bin_distr.err(rr)));
      }

    Print_To_File({}, {TT,RR,cov_tm_L, corr_tm_L}, "../data/R_ratio/"+Tag_reco_type+"/covariance/"+charm_tag+"/cov_L"+V_charm_1_L.Tag[i_ens]+"_tm.dat", "", "");
    Print_To_File({}, {TT,RR,cov_OS_L, corr_OS_L}, "../data/R_ratio/"+Tag_reco_type+"/covariance/"+charm_tag+"/cov_L"+V_charm_1_L.Tag[i_ens]+"_OS.dat", "", "");
    Print_To_File({}, {TT,RR,cov_tm_M, corr_tm_M}, "../data/R_ratio/"+Tag_reco_type+"/covariance/"+charm_tag+"/cov_M"+V_charm_1_L.Tag[i_ens]+"_tm.dat", "", "");
    Print_To_File({}, {TT,RR,cov_OS_M, corr_OS_M}, "../data/R_ratio/"+Tag_reco_type+"/covariance/"+charm_tag+"/cov_M"+V_charm_1_L.Tag[i_ens]+"_OS.dat", "", "");
    Print_To_File({}, {TT,RR,cov_tm_H, corr_tm_H}, "../data/R_ratio/"+Tag_reco_type+"/covariance/"+charm_tag+"/cov_H"+V_charm_1_L.Tag[i_ens]+"_tm.dat", "", "");
    Print_To_File({}, {TT,RR,cov_OS_H, corr_OS_H}, "../data/R_ratio/"+Tag_reco_type+"/covariance/"+charm_tag+"/cov_H"+V_charm_1_L.Tag[i_ens]+"_OS.dat", "", "");

    
  
    Corr.Tmin = Tmin_etaC; Corr.Tmax = Tmax_etaC;
    //P5P5
    //L
    Metac_L_distr = Corr.Fit_distr(Corr.effective_mass_t(pt2_etaC_L.col(0)[i_ens], ""));
    Metac_OS_L_distr = Corr.Fit_distr(Corr.effective_mass_t(pt2_etaC_OS_L.col(0)[i_ens], ""));
    //M
    Metac_M_distr = Corr.Fit_distr(Corr.effective_mass_t(pt2_etaC_M.col(0)[i_ens], ""));
    Metac_OS_M_distr = Corr.Fit_distr(Corr.effective_mass_t(pt2_etaC_OS_M.col(0)[i_ens], ""));
    //H
    Metac_H_distr = Corr.Fit_distr(Corr.effective_mass_t(pt2_etaC_H.col(0)[i_ens], ""));
    Metac_OS_H_distr = Corr.Fit_distr(Corr.effective_mass_t(pt2_etaC_OS_H.col(0)[i_ens], ""));




    //get MV and overlap from vector correlator
    Corr.Tmin= Tmin_VV; Corr.Tmax = Tmax_VV;
    //L
    MV_L_distr = Corr.Fit_distr( Corr.effective_mass_t(V_charm_L_distr, ""));
    overlap_V_L_distr = Za*Za*Corr.Fit_distr(  Corr.residue_t(V_charm_L_distr, ""))/(2.0*MV_L_distr);
    MV_OS_L_distr = Corr.Fit_distr( Corr.effective_mass_t(V_charm_OS_L_distr, ""));
    overlap_OS_V_L_distr = Zv*Zv*Corr.Fit_distr(  Corr.residue_t(V_charm_OS_L_distr, ""))/(2.0*MV_OS_L_distr);
    //M
    MV_M_distr = Corr.Fit_distr( Corr.effective_mass_t(V_charm_M_distr, ""));
    overlap_V_M_distr = Za*Za*Corr.Fit_distr(  Corr.residue_t(V_charm_M_distr, ""))/(2.0*MV_M_distr);
    MV_OS_M_distr = Corr.Fit_distr( Corr.effective_mass_t(V_charm_OS_M_distr, ""));
    overlap_OS_V_M_distr = Zv*Zv*Corr.Fit_distr(  Corr.residue_t(V_charm_OS_M_distr, ""))/(2.0*MV_OS_M_distr);
    //H
    MV_H_distr = Corr.Fit_distr( Corr.effective_mass_t(V_charm_H_distr, ""));
    overlap_V_H_distr = Za*Za*Corr.Fit_distr(  Corr.residue_t(V_charm_H_distr, ""))/(2.0*MV_H_distr);
    MV_OS_H_distr = Corr.Fit_distr( Corr.effective_mass_t(V_charm_OS_H_distr, ""));
    overlap_OS_V_H_distr = Zv*Zv*Corr.Fit_distr(  Corr.residue_t(V_charm_OS_H_distr, ""))/(2.0*MV_OS_H_distr);


    vector<distr_t> Metac_vec, MJpsi_vec;
    if(V_charm_1_L.Tag[i_ens].substr(1,1) != "D") {
      Metac_vec = { Metac_L_distr/a_distr, Metac_M_distr/a_distr, Metac_H_distr/a_distr};
      MJpsi_vec = { MV_L_distr/a_distr, MV_M_distr/a_distr, MV_H_distr/a_distr};
    }
    else {
      Metac_vec = { Metac_L_distr/a_distr, Metac_M_distr/a_distr};
      MJpsi_vec = { MV_L_distr/a_distr, MV_M_distr/a_distr};
    }
    vector<distr_t> X_fit;
    if(Extrapolation_charm_mode == "Jpsi") X_fit = MJpsi_vec;
    else if(Extrapolation_charm_mode == "etac") X_fit = Metac_vec;
    else crash("Extrapolation charm mode: "+Extrapolation_charm_mode+" is invalid");
    distr_t X_phys = (Extrapolation_charm_mode=="Jpsi")?MJpsi_phys_distr:Metac_phys_distr;


   


    //############################################################################################


    //perturbative substraction of lattice artifacts##############################################

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
   
  
    Vfloat free_corr_log_art(Corr.Nt, 0.0);
    for(int t=0;t<Corr.Nt;t++) {  if( t*a_distr.ave() < 1.0*fm_to_inv_Gev && t != 0) {
	free_corr_log_art[t] = -1.0*pert_corr_charm_on_off*(qc*qc)*(1.0/(2.0*M_PI*M_PI*pow(t,5)));
      }
      if(t==0 || t*a_distr.ave() > add_pert_corr_charm_up_to*fm_to_inv_Gev) { VV_free_samer_L[t] =0; VV_free_samer_M[t] = 0; VV_free_samer_H[t] =0; VV_free_oppor_L[t] = 0; VV_free_oppor_M[t]=0; VV_free_oppor_H[t]=0;}
    }


    //####################################################################
    //Read perturbative spectral_density for OS and tm
    //####################################################################
    
    Vfloat Spec_tm_L = Read_From_File("../data/R_ratio/spec_dens_free/tm/am_"+to_string_with_precision(L_info.mc_L,5), 2, 4);
    Vfloat Spec_OS_L = Read_From_File("../data/R_ratio/spec_dens_free/OS/am_"+to_string_with_precision(L_info.mc_L,5), 2, 4);
    Vfloat Spec_tm_M = Read_From_File("../data/R_ratio/spec_dens_free/tm/am_"+to_string_with_precision(L_info.mc_M,5), 2, 4);
    Vfloat Spec_OS_M = Read_From_File("../data/R_ratio/spec_dens_free/OS/am_"+to_string_with_precision(L_info.mc_M,5), 2, 4);
    Vfloat Spec_tm_H = Read_From_File("../data/R_ratio/spec_dens_free/tm/am_"+to_string_with_precision(L_info.mc_H,5), 2, 4);
    Vfloat Spec_OS_H = Read_From_File("../data/R_ratio/spec_dens_free/OS/am_"+to_string_with_precision(L_info.mc_H,5), 2, 4);
    Vfloat Ergs_pert = Read_From_File("../data/R_ratio/spec_dens_free/tm/am_"+to_string_with_precision(L_info.mc_L,5), 1, 4);

    cout<<"perturbative spectral density for Ensemble: "<<V_charm_1_L.Tag[i_ens]<<" READ! "<<endl;
    

    //interpolate perturbative data
    boost::math::interpolators::cardinal_cubic_b_spline<double> F_boost_tm_L(Spec_tm_L.begin(), Spec_tm_L.end(), Ergs_pert[0], 2.0*Ergs_pert[0]);
    boost::math::interpolators::cardinal_cubic_b_spline<double> F_boost_OS_L(Spec_OS_L.begin(), Spec_OS_L.end(), Ergs_pert[0], 2.0*Ergs_pert[0]);
    boost::math::interpolators::cardinal_cubic_b_spline<double> F_boost_tm_M(Spec_tm_M.begin(), Spec_tm_M.end(), Ergs_pert[0], 2.0*Ergs_pert[0]);
    boost::math::interpolators::cardinal_cubic_b_spline<double> F_boost_OS_M(Spec_OS_M.begin(), Spec_OS_M.end(), Ergs_pert[0], 2.0*Ergs_pert[0]);
    boost::math::interpolators::cardinal_cubic_b_spline<double> F_boost_tm_H(Spec_tm_H.begin(), Spec_tm_H.end(), Ergs_pert[0], 2.0*Ergs_pert[0]);
    boost::math::interpolators::cardinal_cubic_b_spline<double> F_boost_OS_H(Spec_OS_H.begin(), Spec_OS_H.end(), Ergs_pert[0], 2.0*Ergs_pert[0]);

    cout<<"Cubic spline for perturbative spectral density for Ensemble: "<<V_charm_1_L.Tag[i_ens]<<" produced! "<<endl;

    auto F_free_tm_L = [&F_boost_tm_L](double E) { return F_boost_tm_L(E);};
    auto F_free_OS_L = [&F_boost_OS_L](double E) { return F_boost_OS_L(E);};
    auto F_free_tm_M = [&F_boost_tm_M](double E) { return F_boost_tm_M(E);};
    auto F_free_OS_M = [&F_boost_OS_M](double E) { return F_boost_OS_M(E);};
    auto F_free_tm_H = [&F_boost_tm_H](double E) { return F_boost_tm_H(E);};
    auto F_free_OS_H = [&F_boost_OS_H](double E) { return F_boost_OS_H(E);};


    //multiply corr using Zv and Za
    //L
    V_charm_L_distr = (V_charm_L_distr)*(1.0 + pert_corr_charm_on_off*(1.0/(Za*Za*V_charm_L_distr))*VV_free_oppor_L) ;
    V_charm_OS_L_distr = (V_charm_OS_L_distr)*(1.0 + pert_corr_charm_on_off*(1.0/(Zv*Zv*V_charm_OS_L_distr))*VV_free_samer_L);
    //M
    V_charm_M_distr = (V_charm_M_distr)*(1.0 + pert_corr_charm_on_off*(1.0/(Za*Za*V_charm_M_distr))*VV_free_oppor_M) ;
    V_charm_OS_M_distr = (V_charm_OS_M_distr)*(1.0 + pert_corr_charm_on_off*(1.0/(Zv*Zv*V_charm_OS_M_distr))*VV_free_samer_M);
    //H
    V_charm_H_distr = (V_charm_H_distr)*(1.0 + pert_corr_charm_on_off*(1.0/(Za*Za*V_charm_H_distr))*VV_free_oppor_H) ;
    V_charm_OS_H_distr = (V_charm_OS_H_distr)*(1.0 + pert_corr_charm_on_off*(1.0/(Zv*Zv*V_charm_OS_H_distr))*VV_free_samer_H);


      


    //print correlator

    //L
    Print_To_File({}, {V_charm_L_distr.ave(), V_charm_L_distr.err(), V_charm_OS_L_distr.ave(), V_charm_OS_L_distr.err()}, "../data/R_ratio/"+Tag_reco_type+"/corr/"+charm_tag+"/corr_L_"+V_charm_1_L.Tag[i_ens]+".dat", "", "# t  tm  OS");
    //M
    Print_To_File({}, {V_charm_M_distr.ave(), V_charm_M_distr.err(), V_charm_OS_M_distr.ave(), V_charm_OS_M_distr.err()}, "../data/R_ratio/"+Tag_reco_type+"/corr/"+charm_tag+"/corr_M_"+V_charm_1_L.Tag[i_ens]+".dat", "", "# t  tm  OS");
    //H
    Print_To_File({}, {V_charm_H_distr.ave(), V_charm_H_distr.err(), V_charm_OS_H_distr.ave(), V_charm_OS_H_distr.err()}, "../data/R_ratio/"+Tag_reco_type+"/corr/"+charm_tag+"/corr_H_"+V_charm_1_H.Tag[i_ens]+".dat", "", "# t  tm  OS");

   
    //####################################################################################################################





    //#########################   RECONSTRUCT THE SMEARED R-RATIO ################################


    //set tmax to the value where the error on V(t) is larger than x%
    //#############################################################################################################################

    //L
    bool Found_error_less_x_percent=false;
    double x=5;
    double tmax_L=1;
    while(!Found_error_less_x_percent && tmax_L < Corr.Nt/2) {
   
      if( (V_charm_L_distr.distr_list[tmax_L]).err()/fabs( (V_charm_L_distr.distr_list[tmax_L]).ave()) <  0.01*x) tmax_L++;
      else Found_error_less_x_percent=true;
    }

    Found_error_less_x_percent=false;
    double tmax_OS_L=1;
    while(!Found_error_less_x_percent && tmax_OS_L < Corr.Nt/2) {
   
      if( (V_charm_OS_L_distr.distr_list[tmax_OS_L]).err()/fabs( (V_charm_OS_L_distr.distr_list[tmax_OS_L]).ave()) <  0.01*x) tmax_OS_L++;
      else Found_error_less_x_percent=true;
    }

    //M
    Found_error_less_x_percent=false;
    double tmax_M=1;
    while(!Found_error_less_x_percent && tmax_M < Corr.Nt/2) {
   
      if( (V_charm_M_distr.distr_list[tmax_M]).err()/fabs( (V_charm_M_distr.distr_list[tmax_M]).ave()) <  0.01*x) tmax_M++;
      else Found_error_less_x_percent=true;
    }

    Found_error_less_x_percent=false;
    double tmax_OS_M=1;
    while(!Found_error_less_x_percent && tmax_OS_M < Corr.Nt/2) {
   
      if( (V_charm_OS_M_distr.distr_list[tmax_OS_M]).err()/fabs( (V_charm_OS_M_distr.distr_list[tmax_OS_M]).ave()) <  0.01*x) tmax_OS_M++;
      else Found_error_less_x_percent=true;
    }

    //H
    Found_error_less_x_percent=false;
    double tmax_H=1;
    while(!Found_error_less_x_percent && tmax_H < Corr.Nt/2) {
   
      if( (V_charm_H_distr.distr_list[tmax_H]).err()/fabs( (V_charm_H_distr.distr_list[tmax_H]).ave()) <  0.01*x) tmax_H++;
      else Found_error_less_x_percent=true;
    }

    Found_error_less_x_percent=false;
    double tmax_OS_H=1;
    while(!Found_error_less_x_percent && tmax_OS_H < Corr.Nt/2) {
   
      if( (V_charm_OS_H_distr.distr_list[tmax_OS_H]).err()/fabs( (V_charm_OS_H_distr.distr_list[tmax_OS_H]).ave()) <  0.01*x) tmax_OS_H++;
      else Found_error_less_x_percent=true;
    }

    if(Use_t_up_to_T_half) {
      tmax_L= tmax_OS_L = tmax_M = tmax_OS_M=tmax_H=tmax_OS_H= T/2 -1;
    }
    
    //#############################################################################################################################



    //#####################################    MODEL ESTIMATE FOR SYSTEMATIC ERRORS ##############################################

    double resc_charm= rho_R*pow(qc,2);

    //tm
    auto model_charm_tm_L = [&a_distr,&MV_L_distr, &F_free_tm_L, &overlap_V_L_distr, &resc_charm](double E) {
      double DE= 0.003*a_distr.ave();
      bool Is_pert= (E > 3.1*a_distr.ave());
      return resc_charm*(overlap_V_L_distr.ave()*((1.0/sqrt(2*M_PI*DE*DE))*exp(- pow(E-MV_L_distr.ave(),2)/(2*DE*DE))) + Is_pert*F_free_tm_L(E));
    };
    auto model_charm_tm_M = [&a_distr,&MV_M_distr, &F_free_tm_M, &overlap_V_M_distr, &resc_charm](double E) {
      double DE= 0.003*a_distr.ave();
      bool Is_pert= (E > 3.1*a_distr.ave());
      return resc_charm*(overlap_V_M_distr.ave()*((1.0/sqrt(2*M_PI*DE*DE))*exp(- pow(E-MV_M_distr.ave(),2)/(2*DE*DE))) + Is_pert*F_free_tm_M(E));
    };
    auto model_charm_tm_H = [&a_distr,&MV_H_distr, &F_free_tm_H, &overlap_V_H_distr, &resc_charm](double E) {
      double DE= 0.003*a_distr.ave();
      bool Is_pert= (E > 3.0*a_distr.ave());
      return resc_charm*(overlap_V_H_distr.ave()*((1.0/sqrt(2*M_PI*DE*DE))*exp(- pow(E-MV_H_distr.ave(),2)/(2*DE*DE))) + Is_pert*F_free_tm_H(E));
    };
    //OS
    auto model_charm_OS_L = [&a_distr,&MV_OS_L_distr, &F_free_OS_L, &overlap_OS_V_L_distr, &resc_charm](double E) {
      double DE= 0.003*a_distr.ave();
      bool Is_pert= (E > 3.1*a_distr.ave());
      return resc_charm*(overlap_OS_V_L_distr.ave()*((1.0/sqrt(2*M_PI*DE*DE))*exp(- pow(E-MV_OS_L_distr.ave(),2)/(2*DE*DE))) + Is_pert*F_free_OS_L(E));
    };
    auto model_charm_OS_M = [&a_distr,&MV_OS_M_distr, &F_free_OS_M, &overlap_OS_V_M_distr, &resc_charm](double E) {
      double DE= 0.003*a_distr.ave();
      bool Is_pert= (E > 3.1*a_distr.ave());
      return resc_charm*(overlap_OS_V_M_distr.ave()*((1.0/sqrt(2*M_PI*DE*DE))*exp(- pow(E-MV_OS_M_distr.ave(),2)/(2*DE*DE))) + Is_pert*F_free_OS_M(E));
    };
    auto model_charm_OS_H = [&a_distr,&MV_OS_H_distr, &F_free_OS_H, &overlap_OS_V_H_distr, &resc_charm](double E) {
      double DE= 0.003*a_distr.ave();
      bool Is_pert= (E > 3.1*a_distr.ave());
      return resc_charm*(overlap_OS_V_H_distr.ave()*((1.0/sqrt(2*M_PI*DE*DE))*exp(- pow(E-MV_OS_H_distr.ave(),2)/(2*DE*DE))) + Is_pert*F_free_OS_H(E));
    };
      
      


    //tm
    auto f_syst_tm_L = [&overlap_V_L_distr, &MV_L_distr, &F_free_tm_L, &resc_charm, &a_distr](const function<double(double)> &F) ->double {
      double val_mod, err_mod;
      auto FS= [&F, &a_distr, &F_free_tm_L](double E) { if (E> 3.1*a_distr.ave()) return F_free_tm_L(E)*F(E); return 0.0;}; 
      gsl_function_pp<decltype(FS)> SYST(FS);
      gsl_integration_workspace * w_SYST = gsl_integration_workspace_alloc (1000);
      gsl_function *G_SYST = static_cast<gsl_function*>(&SYST);
      gsl_integration_qags(G_SYST, Eth*a_distr.ave(), 4.0,  0.0, 7e-3, 1000, w_SYST, &val_mod, &err_mod);
      return resc_charm*fabs((val_mod + F(MV_L_distr.ave())*overlap_V_L_distr.ave()));
      
    };

    auto f_syst_tm_M = [&overlap_V_M_distr, &MV_M_distr, &F_free_tm_M, &resc_charm, &a_distr](const function<double(double)> &F) ->double {
      double val_mod, err_mod;
      auto FS= [&F, &a_distr, &F_free_tm_M](double E) { if (E> 3.1*a_distr.ave()) return F_free_tm_M(E)*F(E); return 0.0;}; 
      gsl_function_pp<decltype(FS)> SYST(FS);
      gsl_integration_workspace * w_SYST = gsl_integration_workspace_alloc (1000);
      gsl_function *G_SYST = static_cast<gsl_function*>(&SYST);
      gsl_integration_qags(G_SYST, Eth*a_distr.ave(), 4.0,  0.0, 7e-3, 1000, w_SYST, &val_mod, &err_mod);
      return resc_charm*fabs((val_mod + F(MV_M_distr.ave())*overlap_V_M_distr.ave()));
      
    };

    auto f_syst_tm_H = [&overlap_V_H_distr, &MV_H_distr, &F_free_tm_H, &resc_charm, &a_distr](const function<double(double)> &F) ->double {
      double val_mod, err_mod;
      auto FS= [&F, &a_distr, &F_free_tm_H](double E) { if (E> 3.1*a_distr.ave()) return F_free_tm_H(E)*F(E); return 0.0;}; 
      gsl_function_pp<decltype(FS)> SYST(FS);
      gsl_integration_workspace * w_SYST = gsl_integration_workspace_alloc (1000);
      gsl_function *G_SYST = static_cast<gsl_function*>(&SYST);
      gsl_integration_qags(G_SYST, Eth*a_distr.ave(), 4.0,  0.0, 7e-3, 1000, w_SYST, &val_mod, &err_mod);
      return resc_charm*fabs((val_mod + F(MV_H_distr.ave())*overlap_V_H_distr.ave()));
      
    };


    //OS
    auto f_syst_OS_L = [&overlap_OS_V_L_distr, &MV_OS_L_distr, &F_free_OS_L, &resc_charm, &a_distr](const function<double(double)> &F) ->double {
      double val_mod, err_mod;
      auto FS= [&F, &a_distr, &F_free_OS_L](double E) { if (E> 3.1*a_distr.ave()) return F_free_OS_L(E)*F(E); return 0.0;}; 
      gsl_function_pp<decltype(FS)> SYST(FS);
      gsl_integration_workspace * w_SYST = gsl_integration_workspace_alloc (1000);
      gsl_function *G_SYST = static_cast<gsl_function*>(&SYST);
      gsl_integration_qags(G_SYST, Eth*a_distr.ave(), 4.0,  0.0, 7e-3, 1000, w_SYST, &val_mod, &err_mod);
      return resc_charm*fabs((val_mod + F(MV_OS_L_distr.ave())*overlap_OS_V_L_distr.ave()));
      
    };

    auto f_syst_OS_M = [&overlap_OS_V_M_distr, &MV_OS_M_distr, &F_free_OS_M, &resc_charm, &a_distr](const function<double(double)> &F) ->double {
      double val_mod, err_mod;
      auto FS= [&F, &a_distr, &F_free_OS_M](double E) { if (E> 3.1*a_distr.ave()) return F_free_OS_M(E)*F(E); return 0.0;}; 
      gsl_function_pp<decltype(FS)> SYST(FS);
      gsl_integration_workspace * w_SYST = gsl_integration_workspace_alloc (1000);
      gsl_function *G_SYST = static_cast<gsl_function*>(&SYST);
      gsl_integration_qags(G_SYST, Eth*a_distr.ave(), 4.0,  0.0, 7e-3, 1000, w_SYST, &val_mod, &err_mod);
      return resc_charm*fabs((val_mod + F(MV_OS_M_distr.ave())*overlap_OS_V_M_distr.ave()));
      
    };

    auto f_syst_OS_H = [&overlap_OS_V_H_distr, &MV_OS_H_distr, &F_free_OS_H, &resc_charm, &a_distr](const function<double(double)> &F) ->double {
      double val_mod, err_mod;
      auto FS= [&F, &a_distr, &F_free_OS_H](double E) { if (E> 3.1*a_distr.ave()) return F_free_OS_H(E)*F(E); return 0.0;}; 
      gsl_function_pp<decltype(FS)> SYST(FS);
      gsl_integration_workspace * w_SYST = gsl_integration_workspace_alloc (1000);
      gsl_function *G_SYST = static_cast<gsl_function*>(&SYST);
      gsl_integration_qags(G_SYST, Eth*a_distr.ave(), 4.0,  0.0, 7e-3, 1000, w_SYST, &val_mod, &err_mod);
      return resc_charm*fabs((val_mod + F(MV_OS_H_distr.ave())*overlap_OS_V_H_distr.ave()));
      
    };

    //##############################################################################################################################





    

  

    //##################################################
    //#############################SET GENERAL PARAMETERS######################################
    double E0 =  Eth*a_distr.ave(); // two pion treshold ~ 280 MeV
    double mult_SANF = 1e-4;
    double mult_TANT=10.0;
    //#########################################################################################

   

    //loop over sigma

    for(int isg=0; isg<(signed)sigmas.size();isg++) {

      double sigma= sigmas[isg]*a_distr.ave();
      
      cout<<"Reconstructing R(E) for:"<<endl;;
      cout<<"a*E0: "<<E0<<" -> E0: "<<Eth<<" [GeV] "<<endl;
      cout<<"SM_TYPE: "<<SM_TYPE<<endl;
      cout<<"a*sigma(E*): "<<sigma<<" -> sigma(E*): "<<sigma/a_distr.ave()<<" [GeV] "<<endl;
    
        
    
    

      distr_t_list Spectral_dens_L(UseJack, Ergs_GeV_list.size());
      distr_t_list Spectral_dens_OS_L(UseJack, Ergs_GeV_list.size());
      distr_t_list Spectral_dens_SANF_L(UseJack, Ergs_GeV_list.size());
      distr_t_list Spectral_dens_OS_SANF_L(UseJack, Ergs_GeV_list.size());
      distr_t_list Spectral_dens_M(UseJack, Ergs_GeV_list.size());
      distr_t_list Spectral_dens_OS_M(UseJack, Ergs_GeV_list.size());
      distr_t_list Spectral_dens_SANF_M(UseJack, Ergs_GeV_list.size());
      distr_t_list Spectral_dens_OS_SANF_M(UseJack, Ergs_GeV_list.size());
      distr_t_list Spectral_dens_H(UseJack, Ergs_GeV_list.size());
      distr_t_list Spectral_dens_OS_H(UseJack, Ergs_GeV_list.size());
      distr_t_list Spectral_dens_SANF_H(UseJack, Ergs_GeV_list.size());
      distr_t_list Spectral_dens_OS_SANF_H(UseJack, Ergs_GeV_list.size());
      distr_t_list Spectral_dens_Extr(UseJack, Ergs_GeV_list.size());
      distr_t_list Spectral_dens_OS_Extr(UseJack, Ergs_GeV_list.size());
      distr_t_list Spectral_dens_SANF_Extr(UseJack, Ergs_GeV_list.size());
      distr_t_list Spectral_dens_OS_SANF_Extr(UseJack, Ergs_GeV_list.size());

      //######### vector where to store systematic errors
      Vfloat syst_tm_L(Ergs_GeV_list.size());
      Vfloat syst_OS_L(Ergs_GeV_list.size());
      Vfloat syst_tm_SANF_L(Ergs_GeV_list.size());
      Vfloat syst_OS_SANF_L(Ergs_GeV_list.size());
      Vfloat syst_tm_M(Ergs_GeV_list.size());
      Vfloat syst_OS_M(Ergs_GeV_list.size());
      Vfloat syst_tm_SANF_M(Ergs_GeV_list.size());
      Vfloat syst_OS_SANF_M(Ergs_GeV_list.size());
      Vfloat syst_tm_H(Ergs_GeV_list.size());
      Vfloat syst_OS_H(Ergs_GeV_list.size());
      Vfloat syst_tm_SANF_H(Ergs_GeV_list.size());
      Vfloat syst_OS_SANF_H(Ergs_GeV_list.size());

      vector<tuple<int, double, double, double>> thread_times_tm(Ergs_GeV_list.size()), thread_times_OS(Ergs_GeV_list.size());


      //COMPUTE THE SMEARED R-RATIO

      cout<<"Looping over energies"<<flush;
      
      #pragma omp parallel for schedule(dynamic)
      for(int ip=0; ip<(signed)Ergs_GeV_list.size();ip++) {

	double mean = Ergs_GeV_list[ip]*a_distr.ave();
	double lambda_Estar;
	double lambda_Estar_SANF;
	double Ag_ov_A0_target= 1e-6;
 
	


	//L (T)
	//define jackknife distribution to account for systematic error:
	distr_t syst_L_T_tm(UseJack), syst_L_T_OS(UseJack);
	for(int ijack=0; ijack<Njacks;ijack++) {syst_L_T_tm.distr.push_back( GM()/sqrt(Njacks-1.0)); syst_L_T_OS.distr.push_back( GM()/sqrt(Njacks-1.0));}

	auto start = chrono::system_clock::now();
	Spectral_dens_L.distr_list[ip] = Get_Laplace_transfo(  mean,  sigma, E0,  T, tmax_L, prec_charm, SM_TYPE+"_ov_E2",f, V_charm_L_distr, syst_tm_L[ip], mult_TANT, lambda_Estar, "TANT", "tm", "L_"+V_charm_1_L.Tag[i_ens], Ag_ov_A0_target, 0, rho_R*Za*Za*pow(qc,2), 0.0, "R_ratio_charm", cov_tm_L, f_syst_tm_L, 0, model_charm_tm_L, Is_Emax_Finite, Emax,beta ) + syst_L_T_tm*syst_tm_L[ip] ;
	auto end = chrono::system_clock::now();
	chrono::duration<double> elapsed_seconds = end-start;
	double time_L_tm= elapsed_seconds.count();
	if(R_ratio_verbosity_lev) {
	  cout<<endl<<flush;
	  cout<<"Elapsed time[L_tm, #thread="<<omp_get_thread_num()<<"] : "<<time_L_tm<<" s"<<endl<<flush;
	}
	cout<<"."<<flush;

	start = chrono::system_clock::now();
	Spectral_dens_OS_L.distr_list[ip] = Get_Laplace_transfo(  mean,  sigma, E0,  T, tmax_OS_L, prec_charm, SM_TYPE+"_ov_E2",f, V_charm_OS_L_distr, syst_OS_L[ip], mult_TANT, lambda_Estar, "TANT", "OS", "L_"+V_charm_1_L.Tag[i_ens], Ag_ov_A0_target, 0, rho_R*Zv*Zv*pow(qc,2), 0.0, "R_ratio_charm", cov_OS_L, f_syst_OS_L, 0, model_charm_OS_L, Is_Emax_Finite, Emax,beta  ) + syst_L_T_OS*syst_OS_L[ip] ;
	end = chrono::system_clock::now();
	elapsed_seconds = end-start;
	double time_L_OS= elapsed_seconds.count();
	if(R_ratio_verbosity_lev) {
	  cout<<endl<<flush;
	  cout<<"Elapsed time[L_OS, #thread="<<omp_get_thread_num()<<"] : "<<time_L_OS<<" s"<<endl<<flush;
	}
	cout<<"."<<flush;


	//L (S)
	if(!SANF_MODE_OFF) {
	  //define jackknife distribution to account for systematic error:
	  distr_t syst_L_S_tm(UseJack), syst_L_S_OS(UseJack);
	  for(int ijack=0; ijack<Njacks;ijack++) {syst_L_S_tm.distr.push_back( GM()/sqrt(Njacks-1.0)); syst_L_S_OS.distr.push_back( GM()/sqrt(Njacks-1.0));}
	  Spectral_dens_SANF_L.distr_list[ip] = Get_Laplace_transfo(  mean,  sigma, E0,  T, tmax_L, prec_charm, SM_TYPE+"_ov_E2",f, V_charm_L_distr, syst_tm_SANF_L[ip], mult_SANF, lambda_Estar_SANF, "SANF", "tm", "L_"+V_charm_1_L.Tag[i_ens], Ag_ov_A0_target, 0, rho_R*Za*Za*pow(qc,2), 0.0, "R_ratio_charm", cov_tm_L, fake_func, 0, fake_func_d, Is_Emax_Finite, Emax,beta  )+ syst_tm_SANF_L[ip]*syst_L_S_tm;
	  Spectral_dens_OS_SANF_L.distr_list[ip] = Get_Laplace_transfo(  mean,  sigma, E0,  T, tmax_OS_L, prec_charm, SM_TYPE+"_ov_E2",f, V_charm_OS_L_distr, syst_OS_SANF_L[ip], mult_SANF, lambda_Estar_SANF, "SANF", "OS", "L_"+V_charm_1_L.Tag[i_ens], Ag_ov_A0_target, 0, rho_R*Zv*Zv*pow(qc,2), 0.0, "R_ratio_charm", cov_OS_L, fake_func, 0, fake_func_d, Is_Emax_Finite, Emax,beta  )+ syst_OS_SANF_L[ip]*syst_L_S_OS ;
	}



	
	//M (T)
	//define jackknife distribution to account for systematic error:
	distr_t syst_M_T_tm(UseJack), syst_M_T_OS(UseJack);
	for(int ijack=0; ijack<Njacks;ijack++) {syst_M_T_tm.distr.push_back( GM()/sqrt(Njacks-1.0)); syst_M_T_OS.distr.push_back( GM()/sqrt(Njacks-1.0));}
	
	start = chrono::system_clock::now();
	Spectral_dens_M.distr_list[ip] = Get_Laplace_transfo(  mean,  sigma, E0,  T, tmax_M, prec_charm, SM_TYPE+"_ov_E2",f, V_charm_M_distr, syst_tm_M[ip], mult_TANT, lambda_Estar, "TANT", "tm", "M_"+V_charm_1_M.Tag[i_ens], Ag_ov_A0_target , 0, rho_R*Za*Za*pow(qc,2), 0.0, "R_ratio_charm", cov_tm_M, f_syst_tm_M, 0, model_charm_tm_M, Is_Emax_Finite, Emax,beta ) + syst_tm_M[ip]*syst_M_T_tm ;
	end = chrono::system_clock::now();
	elapsed_seconds = end-start;
	double time_M_tm= elapsed_seconds.count();
	if(R_ratio_verbosity_lev) {
	  cout<<endl<<flush;
	  cout<<"Elapsed time[M_tm, #thread="<<omp_get_thread_num()<<"] : "<<time_M_tm<<" s"<<endl<<flush;
	}
	cout<<"."<<flush;

	start = chrono::system_clock::now();
	Spectral_dens_OS_M.distr_list[ip] = Get_Laplace_transfo(  mean,  sigma, E0,  T, tmax_OS_M, prec_charm, SM_TYPE+"_ov_E2",f, V_charm_OS_M_distr, syst_OS_M[ip], mult_TANT, lambda_Estar, "TANT", "OS", "M_"+V_charm_1_M.Tag[i_ens], Ag_ov_A0_target , 0, rho_R*Zv*Zv*pow(qc,2), 0.0, "R_ratio_charm", cov_OS_M, f_syst_OS_M, 0, model_charm_OS_M, Is_Emax_Finite, Emax,beta  )+ syst_OS_M[ip]*syst_M_T_OS ;
	end = chrono::system_clock::now();
	elapsed_seconds = end-start;
	double time_M_OS= elapsed_seconds.count();
	if(R_ratio_verbosity_lev) {
	  cout<<endl<<flush;
	  cout<<"Elapsed time[M_OS, #thread="<<omp_get_thread_num()<<"] : "<<time_M_OS<<" s"<<endl<<flush;
	}
	cout<<"."<<flush;



	//M (S)
	if(!SANF_MODE_OFF) {
	  //define jackknife distribution to account for systematic error:
	  distr_t syst_M_S_tm(UseJack), syst_M_S_OS(UseJack);
	  for(int ijack=0; ijack<Njacks;ijack++) {syst_M_S_tm.distr.push_back( GM()/sqrt(Njacks-1.0)); syst_M_S_OS.distr.push_back( GM()/sqrt(Njacks-1.0));}
	  Spectral_dens_SANF_M.distr_list[ip] = Get_Laplace_transfo(  mean,  sigma, E0,  T, tmax_M, prec_charm, SM_TYPE+"_ov_E2",f, V_charm_M_distr, syst_tm_SANF_M[ip], mult_SANF, lambda_Estar_SANF, "SANF", "tm", "M_"+V_charm_1_M.Tag[i_ens], Ag_ov_A0_target , 0, rho_R*Za*Za*pow(qc,2), 0.0, "R_ratio_charm", cov_tm_M, fake_func, 0, fake_func_d, Is_Emax_Finite, Emax,beta  )+ syst_tm_SANF_M[ip]*syst_M_S_tm;
	  Spectral_dens_OS_SANF_M.distr_list[ip] = Get_Laplace_transfo(  mean,  sigma, E0,  T, tmax_OS_M, prec_charm, SM_TYPE+"_ov_E2",f, V_charm_OS_M_distr, syst_OS_SANF_M[ip], mult_SANF, lambda_Estar_SANF, "SANF", "OS", "M_"+V_charm_1_M.Tag[i_ens], Ag_ov_A0_target , 0, rho_R*Zv*Zv*pow(qc,2), 0.0, "R_ratio_charm", cov_OS_M, fake_func, 0, fake_func_d, Is_Emax_Finite, Emax,beta )+ syst_OS_SANF_M[ip]*syst_M_S_OS;
	}

	
	//H (T)
	//define jackknife distribution to account for systematic error:
	distr_t syst_H_T_tm(UseJack), syst_H_T_OS(UseJack);
	for(int ijack=0; ijack<Njacks;ijack++) {syst_H_T_tm.distr.push_back( GM()/sqrt(Njacks-1.0)); syst_H_T_OS.distr.push_back( GM()/sqrt(Njacks-1.0));}

	start = chrono::system_clock::now();
	Spectral_dens_H.distr_list[ip] = Get_Laplace_transfo(  mean,  sigma, E0,  T, tmax_H, prec_charm, SM_TYPE+"_ov_E2",f, V_charm_H_distr, syst_tm_H[ip], mult_TANT, lambda_Estar, "TANT", "tm", "H_"+V_charm_1_H.Tag[i_ens], Ag_ov_A0_target, 0, rho_R*Za*Za*pow(qc,2), 0.0, "R_ratio_charm", cov_tm_H, f_syst_tm_H, 0, model_charm_tm_H, Is_Emax_Finite, Emax,beta  )+ syst_tm_H[ip]*syst_H_T_tm ;
	end = chrono::system_clock::now();
	elapsed_seconds = end-start;
	double time_H_tm= elapsed_seconds.count();
	if(R_ratio_verbosity_lev) {
	  cout<<endl<<flush;
	  cout<<"Elapsed time[H_tm, #thread="<<omp_get_thread_num()<<"] : "<<time_H_tm<<" s"<<endl<<flush;
	}
	cout<<"."<<flush;

	start = chrono::system_clock::now();
	Spectral_dens_OS_H.distr_list[ip] = Get_Laplace_transfo(  mean,  sigma, E0,  T, tmax_OS_H, prec_charm, SM_TYPE+"_ov_E2",f, V_charm_OS_H_distr, syst_OS_H[ip], mult_TANT, lambda_Estar, "TANT", "OS","H_"+V_charm_1_H.Tag[i_ens], Ag_ov_A0_target, 0, rho_R*Zv*Zv*pow(qc,2), 0.0, "R_ratio_charm", cov_OS_H, f_syst_OS_H, 0, model_charm_OS_H, Is_Emax_Finite, Emax,beta  )+ syst_OS_H[ip]*syst_H_T_OS ;
	end = chrono::system_clock::now();
	elapsed_seconds = end-start;
	double time_H_OS= elapsed_seconds.count();
	if(R_ratio_verbosity_lev) {
	  cout<<endl<<flush;
	  cout<<"Elapsed time[H_OS, #thread="<<omp_get_thread_num()<<"] : "<<time_H_OS<<" s"<<endl<<flush;
	}
	cout<<"."<<flush;

	//H (S)
	if(!SANF_MODE_OFF) {
	  distr_t syst_H_S_tm(UseJack), syst_H_S_OS(UseJack);
	  for(int ijack=0; ijack<Njacks;ijack++) {syst_H_S_tm.distr.push_back( GM()/sqrt(Njacks-1.0)); syst_H_S_OS.distr.push_back( GM()/sqrt(Njacks-1.0));}
	  Spectral_dens_SANF_H.distr_list[ip] = Get_Laplace_transfo(  mean,  sigma, E0,  T, tmax_H, prec_charm, SM_TYPE+"_ov_E2",f, V_charm_H_distr, syst_tm_SANF_H[ip], mult_SANF, lambda_Estar_SANF, "SANF", "tm", "H_"+V_charm_1_H.Tag[i_ens], Ag_ov_A0_target, 0, rho_R*Za*Za*pow(qc,2), 0.0, "R_ratio_charm", cov_tm_H, fake_func, 0, fake_func_d, Is_Emax_Finite, Emax,beta  )+ syst_tm_SANF_H[ip]*syst_H_S_tm;
	  Spectral_dens_OS_SANF_H.distr_list[ip] = Get_Laplace_transfo(  mean,  sigma, E0,  T, tmax_OS_H, prec_charm, SM_TYPE+"_ov_E2",f, V_charm_OS_H_distr, syst_OS_SANF_H[ip], mult_SANF, lambda_Estar_SANF, "SANF", "OS", "H_"+V_charm_1_H.Tag[i_ens], Ag_ov_A0_target , 0, rho_R*Zv*Zv*pow(qc,2), 0.0, "R_ratio_charm", cov_OS_H, fake_func, 0, fake_func_d, Is_Emax_Finite, Emax,beta )+ syst_OS_SANF_H[ip]*syst_H_S_OS ;
	}
      

	thread_times_tm[ip] = make_tuple( omp_get_thread_num(), time_L_tm, time_M_tm, time_H_tm);
	thread_times_tm[ip] = make_tuple( omp_get_thread_num(), time_L_OS, time_M_OS, time_H_OS);

	//Extrapolate to the physical charm point
	vector<distr_t> Y_fit_tm, Y_fit_OS, Y_fit_SANF_tm, Y_fit_SANF_OS;
	if(V_charm_1_L.Tag[i_ens].substr(1,1) != "D") {
	  Y_fit_tm = {Spectral_dens_L.distr_list[ip], Spectral_dens_M.distr_list[ip], Spectral_dens_H.distr_list[ip]};
	  Y_fit_OS = {Spectral_dens_OS_L.distr_list[ip], Spectral_dens_OS_M.distr_list[ip], Spectral_dens_OS_H.distr_list[ip]};
	  if(!SANF_MODE_OFF) {
	  Y_fit_SANF_tm = {Spectral_dens_SANF_L.distr_list[ip], Spectral_dens_SANF_M.distr_list[ip], Spectral_dens_SANF_H.distr_list[ip]};
	  Y_fit_SANF_OS = {Spectral_dens_OS_SANF_L.distr_list[ip], Spectral_dens_OS_SANF_M.distr_list[ip], Spectral_dens_OS_SANF_H.distr_list[ip]};
	  }
	}
	else {
	  Y_fit_tm = {Spectral_dens_L.distr_list[ip], Spectral_dens_M.distr_list[ip]};
	  Y_fit_OS = {Spectral_dens_OS_L.distr_list[ip], Spectral_dens_OS_M.distr_list[ip]};
	  if(!SANF_MODE_OFF) {
	  Y_fit_SANF_tm = {Spectral_dens_SANF_L.distr_list[ip], Spectral_dens_SANF_M.distr_list[ip]};
	  Y_fit_SANF_OS = {Spectral_dens_OS_SANF_L.distr_list[ip], Spectral_dens_OS_SANF_M.distr_list[ip]};
	  }
	}
	
	//Extrapolation at physical charm quark mass (T)
	Spectral_dens_Extr.distr_list[ip] =  Obs_extrapolation_meson_mass( Y_fit_tm, X_fit, X_phys ,  "../data/R_ratio/"+Tag_reco_type+"/"+charm_tag  , "R_extr_tm_TANT_"+V_charm_1_L.Tag[i_ens]+"_sigma_"+to_string_with_precision(sigmas[isg],3)+"_Estar_"+to_string_with_precision(Ergs_GeV_list[ip],3)+".dat",  UseJack, "SPLINE" );
	Spectral_dens_OS_Extr.distr_list[ip] =   Obs_extrapolation_meson_mass( Y_fit_OS, X_fit, X_phys ,  "../data/R_ratio/"+Tag_reco_type+"/"+charm_tag  , "R_extr_OS_TANT_"+V_charm_1_L.Tag[i_ens]+"_sigma_"+to_string_with_precision(sigmas[isg],3)+"_Estar_"+to_string_with_precision(Ergs_GeV_list[ip],3)+".dat",  UseJack, "SPLINE" );

	//Extrapolation at physical charm quark mass (S)
	if(!SANF_MODE_OFF) {
	Spectral_dens_SANF_Extr.distr_list[ip] =  Obs_extrapolation_meson_mass( Y_fit_SANF_tm, X_fit, X_phys ,  "../data/R_ratio/"+Tag_reco_type+"/"+charm_tag  , "R_extr_tm_SANF_"+V_charm_1_L.Tag[i_ens]+"_sigma_"+to_string_with_precision(sigmas[isg],3)+"_Estar_"+to_string_with_precision(Ergs_GeV_list[ip],3)+".dat",  UseJack, "SPLINE" );
	Spectral_dens_OS_SANF_Extr.distr_list[ip]=  Obs_extrapolation_meson_mass( Y_fit_SANF_OS, X_fit, X_phys ,  "../data/R_ratio/"+Tag_reco_type+"/"+charm_tag  , "R_extr_OS_SANF_"+V_charm_1_L.Tag[i_ens]+"_sigma_"+to_string_with_precision(sigmas[isg],3)+"_Estar_"+to_string_with_precision(Ergs_GeV_list[ip],3)+".dat",  UseJack, "SPLINE" );
	}
	
	

     
	//############################################################################################


      }

      cout<<"done"<<endl<<flush;
      cout<<"Summary of performances: "<<endl<<flush;
      cout<<"Erg #thread  L    M    H"<<endl<<flush;
      cout<<"- - - - - - - - - - - - - - - - - - - - - "<<endl<<flush;
      for(int ip=0; ip < (signed)Ergs_GeV_list.size();ip++) {
	cout<<Ergs_GeV_list[ip]<<","<<get<0>(thread_times_tm[ip])<<": "<<get<1>(thread_times_tm[ip])<<" s, "<<get<2>(thread_times_tm[ip])<<" s, "<<get<3>(thread_times_tm[ip])<<" s"<<endl<<flush;
	cout<<Ergs_GeV_list[ip]<<","<<get<0>(thread_times_OS[ip])<<": "<<get<1>(thread_times_OS[ip])<<" s, "<<get<2>(thread_times_OS[ip])<<" s, "<<get<3>(thread_times_OS[ip])<<" s"<<endl<<flush;
	cout<<"- - - - - - - - - - - - - - - - - - - - "<<endl<<flush;
      }

      if(SANF_MODE_OFF) {
	Spectral_dens_SANF_Extr= Spectral_dens_Extr; Spectral_dens_OS_SANF_Extr = Spectral_dens_OS_Extr;
	Spectral_dens_SANF_L= Spectral_dens_L; Spectral_dens_OS_SANF_L = Spectral_dens_OS_L;
	Spectral_dens_SANF_M= Spectral_dens_M; Spectral_dens_OS_SANF_M = Spectral_dens_OS_M;
	Spectral_dens_SANF_H= Spectral_dens_H; Spectral_dens_OS_SANF_H = Spectral_dens_OS_H;
      }

      RE_charm_TANT_tm[isg][i_ens] = Spectral_dens_Extr;
      RE_charm_TANT_OS[isg][i_ens] = Spectral_dens_OS_Extr;
      RE_charm_SANF_tm[isg][i_ens] = Spectral_dens_SANF_Extr;
      RE_charm_SANF_OS[isg][i_ens] = Spectral_dens_OS_SANF_Extr;

      cout<<endl;
      cout<<"printing output charm for sigma "<<sigmas[isg]<<" ..."<<flush;

      //print to file
      //L
      Print_To_File({}, {Ergs_GeV_list,  Spectral_dens_L.ave(), Spectral_dens_L.err(),  Spectral_dens_SANF_L.ave(), Spectral_dens_SANF_L.err(), Spectral_dens_OS_L.ave(), Spectral_dens_OS_L.err(),  Spectral_dens_OS_SANF_L.ave(), Spectral_dens_OS_SANF_L.err()}, "../data/R_ratio/"+Tag_reco_type+"/"+charm_tag+"/L_"+SM_TYPE+"_"+V_charm_1_L.Tag[i_ens]+"_sigma_"+to_string_with_precision(sigmas[isg],3)+".dat", "", "#E*(GeV)   R(E)_tm[T]  R(E)_tm (<> stat. ) [S]   R(E)_OS (<> stat. ) [T]  R(E)_OS (<> stat. ) [S]");
      //M
      Print_To_File({}, {Ergs_GeV_list,  Spectral_dens_M.ave(), Spectral_dens_M.err(), Spectral_dens_SANF_M.ave(), Spectral_dens_SANF_M.err(), Spectral_dens_OS_M.ave(), Spectral_dens_OS_M.err(),  Spectral_dens_OS_SANF_M.ave(), Spectral_dens_OS_SANF_M.err()}, "../data/R_ratio/"+Tag_reco_type+"/"+charm_tag+"/M_"+SM_TYPE+"_"+V_charm_1_M.Tag[i_ens]+"_sigma_"+to_string_with_precision(sigmas[isg],3)+".dat", "", "#E*(GeV)   R(E)_tm (<> stat ) [T]  R(E)_tm[S]   R(E)_OS[T]  R(E)_OS[S]");
      //H
      Print_To_File({}, {Ergs_GeV_list,  Spectral_dens_H.ave(), Spectral_dens_H.err(), Spectral_dens_SANF_H.ave(), Spectral_dens_SANF_H.err(), Spectral_dens_OS_H.ave(), Spectral_dens_OS_H.err(), Spectral_dens_OS_SANF_H.ave(), Spectral_dens_OS_SANF_H.err()}, "../data/R_ratio/"+Tag_reco_type+"/"+charm_tag+"/H_"+SM_TYPE+"_"+V_charm_1_H.Tag[i_ens]+"_sigma_"+to_string_with_precision(sigmas[isg],3)+".dat", "", "#E*(GeV)   R(E)_tm[T]  R(E)_tm[S]   R(E)_OS[T]  R(E)_OS[S]");
      //Extr
      Print_To_File({}, {Ergs_GeV_list,  Spectral_dens_Extr.ave(), Spectral_dens_Extr.err(), Spectral_dens_SANF_Extr.ave(), Spectral_dens_SANF_Extr.err(), Spectral_dens_OS_Extr.ave(), Spectral_dens_OS_Extr.err(), Spectral_dens_OS_SANF_Extr.ave(), Spectral_dens_OS_SANF_Extr.err() }, "../data/R_ratio/"+Tag_reco_type+"/"+charm_tag+"/Extr_"+SM_TYPE+"_"+V_charm_1_M.Tag[i_ens]+"_sigma_"+to_string_with_precision(sigmas[isg],3)+".dat", "", "#E*(GeV)   R(E)_tm[T]  R(E)_tm[S]   R(E)_OS[T]  R(E)_OS[S]");

      cout<<"done!"<<endl<<flush;


    }

  }
  

  cout<<"Calculation of charm contribution COMPLETED"<<endl;

  }


  //#############################################################################################
  //#############################################################################################
  //#############################################################################################
  //#############################################################################################
  
  //################## CALCULATION OF THE LIGHT CONTRIBUTION ##################################

  
  //#############################################################################################
  //#############################################################################################
  //#############################################################################################
  //#############################################################################################

if(!skip_light) {

  cout<<"STARTING COMPUTATION OF LIGHT CONTRIBUTION:"<<endl;
  cout<<"PRECISION: "<<prec/ln2_10<<" digits."<<endl;
  
  for(int i_ens=0;i_ens<Nens_light;i_ens++) {

 
    
    Light_Ens_Tags[i_ens] = V_light_1.Tag[i_ens];

    string FLAV= "light";
    string light_tag = "light";
    
    CorrAnalysis Corr(UseJack, Njacks,Nboots);
    CorrAnalysis Corr_block_1(0,V_light_1.Nconfs[i_ens], Nboots, i_ens);
    Corr_block_1.Nt = V_light_1.nrows[i_ens];
    Corr.Nt = V_light_1.nrows[i_ens];
    int T = Corr.Nt;

    cout<<"Analyzing Ensemble: "<<V_light_1.Tag[i_ens]<<endl;

    //get lattice spacing
    distr_t a_distr(UseJack);
    distr_t Zv(UseJack), Za(UseJack);
    double Mpi=0;
    double Mpi_err=0;
    if(V_light_1.Tag[i_ens].substr(1,1)=="B") {a_distr=a_B; Zv = ZV_B; Za = ZA_B; Mpi=0.05653312833; Mpi_err=1.430196186e-05;}
    else if(V_light_1.Tag[i_ens].substr(1,1)=="C") {a_distr=a_C; Zv = ZV_C; Za = ZA_C; Mpi=0.04722061628; Mpi_err=3.492993579e-05;}
    else if(V_light_1.Tag[i_ens].substr(1,1)=="D") {a_distr=a_D; Zv = ZV_D; Za = ZA_D; Mpi=0.04062107883; Mpi_err= 2.973916243e-05;}
    else crash("lattice spacing distribution for Ens: "+V_light_1.Tag[i_ens]+" not found");


    //jack distr for Mpi
    distr_t Mpi_distr(UseJack);
    for(int ij=0;ij<Njacks;ij++) Mpi_distr.distr.push_back( Mpi + GM()*Mpi_err/sqrt(Njacks-1.0));
  
    LatticeInfo L_info;
    L_info.LatInfo_new_ens(V_light_1.Tag[i_ens]);
 
       
    distr_t_list  V_light_distr, V_light_OS_distr;
    distr_t_list  V_light_bin_OS_distr, V_light_bin_tm_distr;
    
  

    //vector light sector
    V_light_distr = Corr.corr_t(V_light_1.col(0)[i_ens], ""); //ViVi already averaged
    V_light_OS_distr = Corr.corr_t(V_light_OS_1.col(0)[i_ens], ""); //ViVi already averaged
    V_light_bin_tm_distr=  Corr_block_1.corr_t(V_light_OS_1.col(0)[i_ens], ""); //ViVi already averaged
    V_light_bin_OS_distr = Corr_block_1.corr_t(V_light_OS_1.col(0)[i_ens], ""); //ViVi already averaged



    //print covariance matrix

    Vfloat cov_tm, cov_OS, TT, RR;
    Vfloat corr_tm, corr_OS, cov_red_tm, cov_red_OS;
    //double k_fact= pow(1.0/(pow(qu,2)+ pow(qd,2)),2);
 
    for(int tt=0;tt<Corr.Nt;tt++)
      for(int rr=0;rr<Corr.Nt;rr++) {
	TT.push_back(tt);
	RR.push_back(rr);
	cov_tm.push_back( (V_light_bin_tm_distr.distr_list[tt]%V_light_bin_tm_distr.distr_list[rr]));
	cov_OS.push_back( (V_light_bin_OS_distr.distr_list[tt]%V_light_bin_OS_distr.distr_list[rr]));
	cov_red_tm.push_back( cov_tm[tt*Corr.Nt + rr]/pow(V_light_bin_tm_distr.ave(1),2));
	cov_red_OS.push_back( cov_OS[tt*Corr.Nt + rr]/pow(V_light_bin_OS_distr.ave(1),2));
	corr_tm.push_back( (V_light_bin_tm_distr.distr_list[tt]%V_light_bin_tm_distr.distr_list[rr])/(V_light_bin_tm_distr.err(tt)*V_light_bin_tm_distr.err(rr)));
	corr_OS.push_back( (V_light_bin_OS_distr.distr_list[tt]%V_light_bin_OS_distr.distr_list[rr])/( V_light_bin_OS_distr.err(tt)*V_light_bin_OS_distr.err(rr)));
      }

    Print_To_File({}, {TT,RR,cov_tm, cov_red_tm, corr_tm}, "../data/R_ratio/"+Tag_reco_type+"/covariance/"+light_tag+"/cov_"+V_light_1.Tag[i_ens]+"_tm.dat", "", "");
    Print_To_File({}, {TT,RR,cov_OS, cov_red_OS, corr_OS}, "../data/R_ratio/"+Tag_reco_type+"/covariance/"+light_tag+"/cov_"+V_light_1.Tag[i_ens]+"_OS.dat", "", "");
     

      
    


    //############################################################################################


    //perturbative substraction of lattice artifacts##############################################


    
    //free corr LO artifacts
    //################## READ FREE THEORY VECTOR-VECTOR CORRELATOR OPPOSITE R ####################################
    string Pt_free_oppor= "../Vkvk_cont/"+to_string(Corr.Nt/2)+"_m"+to_string_with_precision(L_info.ml,5)+"/OPPOR";
    Vfloat VV_free_oppor= Read_From_File(Pt_free_oppor, 1, 4);
    if(VV_free_oppor.size() != Corr.Nt) crash("Failed to read properly free VV correlator ml w opposite r");
    //################## READ FREE THEORY VECTOR-VECTOR CORRELATOR SAME R ####################################
    string Pt_free_samer= "../Vkvk_cont/"+to_string(Corr.Nt/2)+"_m"+to_string_with_precision(L_info.ml,5)+"/SAMER";
    Vfloat VV_free_samer= Read_From_File(Pt_free_samer, 1, 4);
    if(VV_free_samer.size() != Corr.Nt) crash("Failed to read properly free VV correlator ml  w same r");
    

    Vfloat free_corr_log_art(Corr.Nt,0);
    for(int t=1;t<Corr.Nt;t++) {
      if(t*a_distr.ave() < add_pert_corr_light_up_to*fm_to_inv_Gev)  free_corr_log_art[t] = -1.0*(qu*qu +qd*qd)*(1.0/(2.0*M_PI*M_PI*pow(t,5)));

      if(t==0 || t*a_distr.ave() > add_pert_corr_light_up_to*fm_to_inv_Gev) { VV_free_samer[t] =0; VV_free_oppor[t] = 0;}


    }


    //####################################################################
    //Read perturbative spectral_density for OS and tm
    //####################################################################
    
    Vfloat Spec_tm = Read_From_File("../data/R_ratio/spec_dens_free/tm/am_"+to_string_with_precision(L_info.ml,5), 2, 4);
    Vfloat Spec_OS = Read_From_File("../data/R_ratio/spec_dens_free/OS/am_"+to_string_with_precision(L_info.ml,5), 2, 4);
    Vfloat Ergs_pert = Read_From_File("../data/R_ratio/spec_dens_free/tm/am_"+to_string_with_precision(L_info.ml,5), 1, 4);

    cout<<"perturbative spectral density for Ensemble: "<<V_light_1.Tag[i_ens]<<" READ! "<<endl;
    

    //interpolate perturbative data
    boost::math::interpolators::cardinal_cubic_b_spline<double> F_boost_tm(Spec_tm.begin(), Spec_tm.end(), Ergs_pert[0], 2.0*Ergs_pert[0]);
    boost::math::interpolators::cardinal_cubic_b_spline<double> F_boost_OS(Spec_OS.begin(), Spec_OS.end(), Ergs_pert[0], 2.0*Ergs_pert[0]);
   
    cout<<"Cubic spline for perturbative spectral density for Ensemble: "<<V_light_1.Tag[i_ens]<<" produced! "<<endl;

    auto F_free_tm = [&F_boost_tm](double E) { return F_boost_tm(E);};
    auto F_free_OS = [&F_boost_OS](double E) { return F_boost_OS(E);};
   
     
    //sum perturbative corrections
    V_light_distr = (V_light_distr)*(1.0 + pert_corr_light_on_off*(1.0/(Za*Za))*(1.0/V_light_distr)*VV_free_oppor);
    V_light_OS_distr = (V_light_OS_distr)*( 1.0+  pert_corr_light_on_off*(1.0/(Zv*Zv))*(1.0/V_light_OS_distr)*VV_free_samer);

  

    //print correlator

    Print_To_File({}, {V_light_distr.ave(), V_light_distr.err(), V_light_OS_distr.ave(), V_light_OS_distr.err()}, "../data/R_ratio/"+Tag_reco_type+"/corr/"+light_tag+"/corr_"+V_light_1.Tag[i_ens]+".dat", "", "# t  tm  OS");


    //############################################################################################





    //#########################   RECONSTRUCT THE SMEARED R-RATIO ################################


    //set tmax to the value where the error on V(t) is larger than x%
    //#############################################################################################################################

    bool Found_error_less_x_percent=false;
    double x=10;
    double tmax=1;
    while(!Found_error_less_x_percent && tmax < Corr.Nt/2) {
   
      if( (V_light_distr.distr_list[tmax]).err()/fabs( (V_light_distr.distr_list[tmax]).ave()) <  0.01*x) tmax++;
      else Found_error_less_x_percent=true;
    }

    Found_error_less_x_percent=false;
    double tmax_OS=1;
    while(!Found_error_less_x_percent && tmax_OS < Corr.Nt/2) {
   
      if( (V_light_OS_distr.distr_list[tmax_OS]).err()/fabs( (V_light_OS_distr.distr_list[tmax_OS]).ave()) <  0.01*x) tmax_OS++;
      else Found_error_less_x_percent=true;
    }

    //###########################################################################################################################





    //#####################################    MODEL ESTIMATE FOR SYSTEMATIC ERRORS ##############################################

    distr_t Edual_tm, Edual_OS, Mrho_tm, Mrho_OS, Rdual_tm, Rdual_OS, grpp_tm, grpp_OS;
    LL.MLLGS_fit_to_corr(Za*Za*V_light_distr, Mpi_distr, a_distr, L_info.L, Edual_tm, Rdual_tm, Mrho_tm, grpp_tm, 5, tmax, V_light_1.Tag[i_ens]+"_tm");
    LL.MLLGS_fit_to_corr(Zv*Zv*V_light_OS_distr, Mpi_distr, a_distr, L_info.L, Edual_OS, Rdual_OS, Mrho_OS, grpp_OS, 7, tmax_OS, V_light_1.Tag[i_ens]+"_OS");

    Vfloat gppis({grpp_tm.ave(), grpp_OS.ave()});
    Vfloat Mrhos({Mrho_tm.ave(), Mrho_OS.ave()});
    Vfloat Eduals({Edual_tm.ave(),Edual_OS.ave()});
    Vfloat Rduals({Rdual_tm.ave(), Rdual_OS.ave()});
    Vfloat En_tm, Ampl_tm;
    Vfloat En_OS, Ampl_OS;
    LL.Find_pipi_energy_lev( L_info.L  , Mrhos[0],  gppis[0], Mpi, 0.0, En_tm);
    LL.Find_pipi_energy_lev( L_info.L  , Mrhos[1],  gppis[1], Mpi, 0.0, En_OS);
    int N=En_tm.size();
    for(int n=0; n<N;n++) {
      Ampl_tm.push_back( 2.0*LL.Amplitude( En_tm[n], L_info.L, Mrhos[0], gppis[0], Mpi, 0.0)); En_tm[n] = 2.0*sqrt( En_tm[n]*En_tm[n] + Mpi*Mpi);
      Ampl_OS.push_back( 2.0*LL.Amplitude( En_OS[n], L_info.L, Mrhos[1], gppis[1], Mpi, 0.0)); En_OS[n] = 2.0*sqrt( En_OS[n]*En_OS[n] + Mpi*Mpi);
    }
    VVfloat Ergs({En_tm, En_OS});
    VVfloat Amplitudes({Ampl_tm, Ampl_OS});

    double resc_light= rho_R*(pow(qu,2)+pow(qd,2));
    
    
    auto GS_V = [ &a_distr, &resc_light](double E, Vfloat &En, Vfloat &Ampl, double Mrho, double Edual, double Rdual) -> double {
		      
		  //build a spectral density with resonances up to 1.5 GeV, from 1.5 GeV use pQCD result. two-pion peaks are smeared over a few MeV interval

		  double result=0.0;
		  double DE= 0.003*a_distr.ave();
		     

		  //pi-pi states
		  for(int n=0; n < (signed)En.size();n++) {
		    if(En[n]< 1.5*a_distr.ave()) {
		      result += resc_light*Ampl[n]*(1.0/sqrt( 2.0*M_PI*DE*DE))*exp( - ( En[n] - E)*(En[n]-E)/(2.0*DE*DE));
		    }
		  }
		  //pQCD part
		  double res_pQCD=0.0;
		 
		  double sth= Mrho+Edual;
		  if(E> sth) {
		    res_pQCD += resc_light*Rdual*(1.0/(2*M_PI*M_PI))*(0.5*pow(E-sth,2) + 0.5*pow(sth,2)+ sth*(E-sth));
		  }
		  result += res_pQCD;
		  
		  
		  return result;
		    };


    
    const auto model_light_tm =  [&GS_V, &Ergs, &Amplitudes, &a_distr, &Mrhos, &Eduals, &Rduals](double E) -> double {  return GS_V(E, Ergs[0], Amplitudes[0], Mrhos[0], Eduals[0], Rduals[0]);};
    const auto model_light_OS =  [&GS_V, &Ergs, &Amplitudes, &a_distr, &Mrhos, &Eduals, &Rduals](double E) -> double {  return GS_V(E, Ergs[1], Amplitudes[1], Mrhos[1], Eduals[1], Rduals[1]);};

    auto f_syst_tm = [&Ergs, &Amplitudes, &Mrhos, &gppis, &Eduals, &Rduals, &GS_V, &a_distr, &Mpi, &resc_light, &F_free_tm](const function<double(double)> &F) ->double {


		     		     
			     		      
		      auto FS = [ &Ergs, &Amplitudes, &Mrhos, &Eduals, &Rduals, &F, &GS_V, &F_free_tm, &resc_light](double E) {
					       double syst = F(E)*GS_V(E, Ergs[0], Amplitudes[0], Mrhos[0], Eduals[0], Rduals[0]);
					       if( E > 1) syst = F(E)*resc_light*F_free_tm(E);
					       return syst;
					     };
		     
		      double val_mod, err_mod;
		      gsl_function_pp<decltype(FS)> SYST(FS);
		      gsl_integration_workspace * w_SYST = gsl_integration_workspace_alloc (1000);
		      gsl_function *G_SYST = static_cast<gsl_function*>(&SYST);
		      gsl_integration_qags(G_SYST, Eth*a_distr.ave(), 4.0,  0.0, 5e-3, 1000, w_SYST, &val_mod, &err_mod);
		      gsl_integration_workspace_free(w_SYST);
		      return fabs(val_mod);
					   
	     
		   
		    };



    auto f_syst_OS = [&Ergs, &Amplitudes, &Mrhos, &gppis, &Eduals, &Rduals, &GS_V, &a_distr, &Mpi, &resc_light, &F_free_OS](const function<double(double)> &F) ->double {

		     		      
		      auto FS = [ &Ergs, &Amplitudes, &Mrhos, &Eduals, &Rduals, &F, &GS_V, &F_free_OS, &resc_light](double E) {
					       double syst = F(E)*GS_V(E, Ergs[1], Amplitudes[1], Mrhos[1], Eduals[1], Rduals[1]);
					       if( E > 1) syst = F(E)*resc_light*F_free_OS(E);
					       return syst;
					     };

		      double val_mod, err_mod;
		      gsl_function_pp<decltype(FS)> SYST(FS);
		      gsl_integration_workspace * w_SYST = gsl_integration_workspace_alloc (1000);
		      gsl_function *G_SYST = static_cast<gsl_function*>(&SYST);
		      gsl_integration_qags(G_SYST, Eth*a_distr.ave(), 4.0,  0.0, 5e-3, 1000, w_SYST, &val_mod, &err_mod);
		      gsl_integration_workspace_free(w_SYST);
		      return fabs(val_mod);
					   
			   
		    };


    //###########################################################################################################################






    
    if(Use_t_up_to_T_half) tmax=tmax_OS= T/2 -1;
    
 
    //##################################################
    //#############################SET GENERAL PARAMETERS######################################
    double E0= Eth*a_distr.ave(); 
    double mult_TANT= 10.0;
    double mult_SANF = 1e-4;
    //#########################################################################################


  
    

    //loop over sigma

    for(int isg=0;isg<(signed)sigmas.size();isg++) {
  

      double sigma=sigmas[isg]*a_distr.ave();

      cout<<"Reconstructing R(E)....."<<endl;
      cout<<"aE0: "<<E0<<" -> E0: "<<E0/a_distr.ave()<<" [GeV] "<<endl;
      cout<<"SM_TYPE: "<<SM_TYPE<<endl;
      cout<<"a*sigma(E*): "<<sigma<<" -> sigma(E*): "<<sigma/a_distr.ave()<<" [GeV] "<<endl;

      
      distr_t_list Spectral_dens(UseJack,Ergs_GeV_list.size());
      distr_t_list Spectral_dens_OS(UseJack,Ergs_GeV_list.size());
      distr_t_list Spectral_dens_SANF(UseJack,Ergs_GeV_list.size());
      distr_t_list Spectral_dens_OS_SANF(UseJack,Ergs_GeV_list.size());
      Vfloat syst_tm(Ergs_GeV_list.size());
      Vfloat syst_OS(Ergs_GeV_list.size());
      Vfloat syst_tm_SANF(Ergs_GeV_list.size());
      Vfloat syst_OS_SANF(Ergs_GeV_list.size());

      vector<tuple<int,  double>> thread_times_tm(Ergs_GeV_list.size()), thread_times_OS(Ergs_GeV_list.size());


      //COMPUTE THE SMEARED R-RATIO

      cout<<"Looping over energies"<<flush;

      #pragma omp parallel for schedule(dynamic)

      for(int ip=0; ip<(signed)Ergs_GeV_list.size();ip++) {

	double mean = Ergs_GeV_list[ip]*a_distr.ave();
	double lambda_Estar;
	double lambda_Estar_SANF;
   
	
	// (T)
        //define jackknife distribution to account for systematic error:
	distr_t syst_T_tm(UseJack), syst_T_OS(UseJack);
	for(int ijack=0; ijack<Njacks;ijack++) {syst_T_tm.distr.push_back( GM()/sqrt(Njacks-1.0)); syst_T_OS.distr.push_back( GM()/sqrt(Njacks-1.0));}

	auto start= chrono::system_clock::now();
	Spectral_dens.distr_list[ip] = Get_Laplace_transfo(  mean,  sigma, E0,  T, tmax, prec, SM_TYPE+"_ov_E2",f, V_light_distr, syst_tm[ip], mult_TANT, lambda_Estar, "TANT", "tm", V_light_1.Tag[i_ens], -1 , 0, rho_R*Za*Za*(pow(qu,2)+pow(qd,2)), 0.0, "R_ratio_light" , cov_tm, f_syst_tm, 1, model_light_tm, Is_Emax_Finite, Emax,beta ) + syst_tm[ip]*syst_T_tm ;
	auto end = chrono::system_clock::now();
	cout<<"node: "<<_hostname<<", rank: "<<rank<<", thread_id: "<<omp_get_thread_num()<<" core-id: "<<sched_getcpu()<<endl<<flush; 
	chrono::duration<double> elapsed_seconds = end-start;
	double time_tm= elapsed_seconds.count();
	if(R_ratio_verbosity_lev) {
	  cout<<endl<<flush;
	  cout<<"Elapsed time[tm, #thread="<<omp_get_thread_num()<<"] : "<<time_tm<<" s"<<endl<<flush;
	}
	cout<<"."<<flush;

	start= chrono::system_clock::now();
	Spectral_dens_OS.distr_list[ip] = Get_Laplace_transfo(  mean,  sigma, E0,  T, tmax_OS, prec, SM_TYPE+"_ov_E2",f, V_light_OS_distr, syst_OS[ip], mult_TANT, lambda_Estar, "TANT", "OS", V_light_1.Tag[i_ens], -1 , 0,  rho_R*Zv*Zv*(pow(qu,2)+pow(qd,2)), 0.0, "R_ratio_light", cov_OS, f_syst_OS, 1, model_light_OS, Is_Emax_Finite, Emax,beta  )+ syst_OS[ip]*syst_T_OS ;
	end = chrono::system_clock::now();
	elapsed_seconds = end-start;
	double time_OS= elapsed_seconds.count();
	if(R_ratio_verbosity_lev) {
	  cout<<endl<<flush;
	  cout<<"Elapsed time[OS, #thread="<<omp_get_thread_num()<<"] : "<<time_OS<<" s"<<endl<<flush;
	}
	cout<<"."<<flush;

	// (S)
	if(!SANF_MODE_OFF) {
	  //define jackknife distribution to account for systematic error:
	  distr_t syst_S_tm(UseJack), syst_S_OS(UseJack);
	  for(int ijack=0; ijack<Njacks;ijack++) {syst_S_tm.distr.push_back( GM()/sqrt(Njacks-1.0)); syst_S_OS.distr.push_back( GM()/sqrt(Njacks-1.0));}
	  Spectral_dens_OS_SANF.distr_list[ip] = Get_Laplace_transfo(  mean,  sigma, E0,  T, tmax_OS, prec, SM_TYPE+"_ov_E2",f, V_light_OS_distr, syst_OS_SANF[ip], mult_SANF, lambda_Estar_SANF, "SANF", "OS", V_light_1.Tag[i_ens], -1 , 0,  rho_R*Zv*Zv*(pow(qu,2)+pow(qd,2)), 0.0, "R_ratio_light", cov_OS, fake_func, 0, fake_func_d, Is_Emax_Finite, Emax,beta )+ syst_OS_SANF[ip]*syst_S_tm ;
	  Spectral_dens_SANF.distr_list[ip] = Get_Laplace_transfo(  mean,  sigma, E0,  T, tmax, prec, SM_TYPE+"_ov_E2",f, V_light_distr, syst_tm_SANF[ip], mult_SANF, lambda_Estar_SANF, "SANF", "tm", V_light_1.Tag[i_ens], -1 , 0,  rho_R*Za*Za*(pow(qu,2)+pow(qd,2)), 0.0, "R_ratio_light", cov_tm, fake_func, 0, fake_func_d, Is_Emax_Finite, Emax,beta  ) + syst_tm_SANF[ip]*syst_S_OS ;
	}
	
     
	//############################################################################################

	thread_times_tm[ip] = make_tuple( omp_get_thread_num(), time_tm);
	thread_times_OS[ip] = make_tuple( omp_get_thread_num(), time_OS);
	

      }

      cout<<"done!"<<endl;
      cout<<"Summary of performances: "<<endl<<flush;
      cout<<"Erg #thread  L"<<endl<<flush;
      cout<<"- - - - - - - - - - - - - - - - - - - - - "<<endl<<flush;
      for(int ip=0; ip < (signed)Ergs_GeV_list.size();ip++) {
	cout<<Ergs_GeV_list[ip]<<","<<get<0>(thread_times_tm[ip])<<": "<<get<1>(thread_times_tm[ip])<<" s"<<endl<<flush;
	cout<<Ergs_GeV_list[ip]<<","<<get<0>(thread_times_OS[ip])<<": "<<get<1>(thread_times_OS[ip])<<" s"<<endl<<flush;
	cout<<"- - - - - - - - - - - - - - - - - - - - "<<endl<<flush;
      }


      if(!test_mode) {

      if(SANF_MODE_OFF) {
	Spectral_dens_SANF= Spectral_dens; Spectral_dens_OS_SANF = Spectral_dens_OS;
      }

    


      RE_light_TANT_tm[isg][i_ens] = Spectral_dens;
      RE_light_TANT_OS[isg][i_ens] = Spectral_dens_OS;
      RE_light_SANF_tm[isg][i_ens] = Spectral_dens_SANF;
      RE_light_SANF_OS[isg][i_ens] = Spectral_dens_OS_SANF;

      cout<<"printing output light for sigma: "<<sigmas[isg]<<" ..."<<flush;

      //light
      //print to file
      Print_To_File({}, {Ergs_GeV_list, Spectral_dens.ave(), Spectral_dens.err(), Spectral_dens_SANF.ave(), Spectral_dens_SANF.err(), Spectral_dens_OS.ave(), Spectral_dens_OS.err(), Spectral_dens_OS_SANF.ave(), Spectral_dens_OS_SANF.err()}, "../data/R_ratio/"+Tag_reco_type+"/"+light_tag+"/"+SM_TYPE+"_"+V_light_1.Tag[i_ens]+"_sigma_"+to_string_with_precision(sigmas[isg],3)+".dat", "", "#E*(GeV)   R(E)_tm [T]  R(E)_tm[S]   R(E)_OS[T]  R(E)_OS[S]");
      }

      cout<<"done!"<<endl;
    }
    


  }


  cout<<"Calculation of light contribution COMPLETED"<<endl;

}

  
  //#############################################################################################
  //#############################################################################################
  //#############################################################################################
  //#############################################################################################
  
  //################## CALCULATION OF THE STRANGE CONTRIBUTION ##################################

  
  //#############################################################################################
  //#############################################################################################
  //#############################################################################################
  //#############################################################################################

 if(!skip_strange) { 

   cout<<"STARTING COMPUTATION OF STRANGE CONTRIBUTION:"<<endl;
   cout<<"PRECISION: "<<prec/ln2_10<<" digits."<<endl;

  for(int i_ens=0;i_ens<Nens_strange;i_ens++) {


    
    
    Strange_Ens_Tags[i_ens] = V_strange_1_L.Tag[i_ens];

    string FLAV= "strange";
    string strange_tag = "strange";
    
    CorrAnalysis Corr(UseJack, Njacks,Nboots);
    CorrAnalysis Corr_block_1(0, V_strange_1_L.Nconfs[i_ens], Nboots, i_ens);
    Corr_block_1.Nt = V_strange_1_L.nrows[i_ens];
    Corr.Nt = V_strange_1_L.nrows[i_ens];
    int T = Corr.Nt;

    cout<<"Analyzing Ensemble: "<<V_strange_1_L.Tag[i_ens]<<endl;
   
    //get lattice spacing
    distr_t a_distr(UseJack);
    distr_t Zv(UseJack), Za(UseJack);
    if(V_strange_1_L.Tag[i_ens].substr(1,1)=="A") {a_distr=a_A; Zv= ZV_A; Za = ZA_A;}
    else if(V_strange_1_L.Tag[i_ens].substr(1,1)=="B") {a_distr=a_B; Zv= ZV_B; Za = ZA_B;}
    else if(V_strange_1_L.Tag[i_ens].substr(1,1)=="C") {a_distr=a_C; Zv= ZV_C; Za = ZA_C;}
    else if(V_strange_1_L.Tag[i_ens].substr(1,1)=="D") {a_distr=a_D; Zv= ZV_D; Za = ZA_D;}
    else crash("lattice spacing distribution for Ens: "+V_strange_1_L.Tag[i_ens]+" not found");



    //set time intervals for pseudoscalar obs
    if(V_strange_1_L.Tag[i_ens].substr(1,1) == "C") {
      if(V_strange_1_L.Tag[i_ens]=="cC211a.06.80") { Corr.Tmin=40; Corr.Tmax=70;}
      else crash("Cannot find ensemble tag: "+V_strange_1_L.Tag[i_ens]);
    }
    else if(V_strange_1_L.Tag[i_ens].substr(1,1) == "B") {
      if(V_strange_1_L.Tag[i_ens]== "cB211a.14.64") {Corr.Tmin=31; Corr.Tmax=58;}
      else if(V_strange_1_L.Tag[i_ens] == "cB211a.25.48") {Corr.Tmin=23;Corr.Tmax=44;}
      else if(V_strange_1_L.Tag[i_ens] == "cB211b.072.64") {Corr.Tmin=36; Corr.Tmax= 57;}
      else if(V_strange_1_L.Tag[i_ens] == "cB211b.072.96") {Corr.Tmin=40; Corr.Tmax= 80;}
      else crash("Cannot find ensemble tag: "+V_strange_1_L.Tag[i_ens]);
    }
    else if(V_strange_1_L.Tag[i_ens].substr(1,1) == "A") {
      if(V_strange_1_L.Tag[i_ens] == "cA211a.12.48") {Corr.Tmin=19; Corr.Tmax=33;}
      else if(V_strange_1_L.Tag[i_ens] == "cA211a.40.24") {Corr.Tmin=18; Corr.Tmax=23;}
      else if(V_strange_1_L.Tag[i_ens] == "cA211a.53.24") {Corr.Tmin=16; Corr.Tmax=22;}
      else if(V_strange_1_L.Tag[i_ens] == "cA211ab.30.32") {Corr.Tmin=23; Corr.Tmax=30;}
      else crash("Cannot find ensemble tag: "+V_strange_1_L.Tag[i_ens]);
    }
    else if(V_strange_1_L.Tag[i_ens].substr(1,1) == "D") {
      if(V_strange_1_L.Tag[i_ens] == "cD211a.054.96") {Corr.Tmin=55; Corr.Tmax=88;}
      else crash("Cannot find ensemble tag: "+V_strange_1_L.Tag[i_ens]);
    }
    else crash("Ensemble tag not valid");

    int Tmin_P5P5 = Corr.Tmin;
    int Tmax_P5P5 = Corr.Tmax;

    int Tmin_VV, Tmax_VV;
    int Tmin_VV_OS, Tmax_VV_OS;

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

    

      
  
    LatticeInfo L_info;
    L_info.LatInfo_new_ens(V_strange_1_L.Tag[i_ens]);
 
    //L
    distr_t_list   V_strange_L_distr, V_strange_L_bin_distr;
    distr_t_list   V_strange_OS_L_distr, V_strange_OS_L_bin_distr;
    distr_t  Metas_L_distr, Metas_OS_L_distr, MV_L_distr, MV_OS_L_distr;
    distr_t overlap_V_L_distr, overlap_OS_V_L_distr;
    //M
    distr_t_list   V_strange_M_distr, V_strange_M_bin_distr;
    distr_t_list   V_strange_OS_M_distr, V_strange_OS_M_bin_distr;
    distr_t Metas_M_distr, Metas_OS_M_distr, MV_M_distr, MV_OS_M_distr;
    distr_t overlap_V_M_distr, overlap_OS_V_M_distr;

     
    

    //vector light sector
    //L
    V_strange_L_distr = Corr.corr_t(V_strange_1_L.col(0)[i_ens], "");
    V_strange_OS_L_distr = Corr.corr_t(V_strange_OS_1_L.col(0)[i_ens], "");
    V_strange_L_bin_distr = Corr_block_1.corr_t(V_strange_1_L.col(0)[i_ens], "");
    V_strange_OS_L_bin_distr = Corr_block_1.corr_t(V_strange_OS_1_L.col(0)[i_ens], "");
    //M
    V_strange_M_distr = Corr.corr_t(V_strange_1_M.col(0)[i_ens], "");
    V_strange_OS_M_distr = Corr.corr_t(V_strange_OS_1_M.col(0)[i_ens], "");
    V_strange_M_bin_distr = Corr_block_1.corr_t(V_strange_1_M.col(0)[i_ens], "");
    V_strange_OS_M_bin_distr = Corr_block_1.corr_t(V_strange_OS_1_M.col(0)[i_ens], "");




    //print covariance matrix

    Vfloat cov_tm_L, cov_OS_L, cov_tm_M, cov_OS_M,  TT, RR;
    Vfloat corr_tm_L, corr_OS_L, corr_tm_M, corr_OS_M;
    
    for(int tt=0;tt<Corr.Nt;tt++)
      for(int rr=0;rr<Corr.Nt;rr++) {
	TT.push_back(tt);
	RR.push_back(rr);
	//L
	cov_tm_L.push_back( (V_strange_L_bin_distr.distr_list[tt]%V_strange_L_bin_distr.distr_list[rr]));
	cov_OS_L.push_back( (V_strange_OS_L_bin_distr.distr_list[tt]%V_strange_OS_L_bin_distr.distr_list[rr]));
	corr_tm_L.push_back( (V_strange_L_bin_distr.distr_list[tt]%V_strange_L_bin_distr.distr_list[rr])/(V_strange_L_bin_distr.err(tt)*V_strange_L_bin_distr.err(rr)));
	corr_OS_L.push_back( (V_strange_OS_L_bin_distr.distr_list[tt]%V_strange_OS_L_bin_distr.distr_list[rr])/( V_strange_OS_L_bin_distr.err(tt)*V_strange_OS_L_bin_distr.err(rr)));

	//M
	cov_tm_M.push_back( (V_strange_M_bin_distr.distr_list[tt]%V_strange_M_bin_distr.distr_list[rr]));
	cov_OS_M.push_back( (V_strange_OS_M_bin_distr.distr_list[tt]%V_strange_OS_M_bin_distr.distr_list[rr]));
	corr_tm_M.push_back( (V_strange_M_bin_distr.distr_list[tt]%V_strange_M_bin_distr.distr_list[rr])/(V_strange_M_bin_distr.err(tt)*V_strange_M_bin_distr.err(rr)));
	corr_OS_M.push_back( (V_strange_OS_M_bin_distr.distr_list[tt]%V_strange_OS_M_bin_distr.distr_list[rr])/( V_strange_OS_M_bin_distr.err(tt)*V_strange_OS_M_bin_distr.err(rr)));


      }

    Print_To_File({}, {TT,RR,cov_tm_L, corr_tm_L}, "../data/R_ratio/"+Tag_reco_type+"/covariance/"+strange_tag+"/cov_L"+V_strange_1_L.Tag[i_ens]+"_tm.dat", "", "");
    Print_To_File({}, {TT,RR,cov_OS_L, corr_OS_L}, "../data/R_ratio/"+Tag_reco_type+"/covariance/"+strange_tag+"/cov_L"+V_strange_1_L.Tag[i_ens]+"_OS.dat", "", "");
    Print_To_File({}, {TT,RR,cov_tm_M, corr_tm_M}, "../data/R_ratio/"+Tag_reco_type+"/covariance/"+strange_tag+"/cov_M"+V_strange_1_L.Tag[i_ens]+"_tm.dat", "", "");
    Print_To_File({}, {TT,RR,cov_OS_M, corr_OS_M}, "../data/R_ratio/"+Tag_reco_type+"/covariance/"+strange_tag+"/cov_M"+V_strange_1_L.Tag[i_ens]+"_OS.dat", "", "");
   

   

    //P5P5
    //L
    Metas_L_distr = Corr.Fit_distr(Corr.effective_mass_t(pt2_etaS_L.col(0)[i_ens], ""));
    Metas_OS_L_distr = Corr.Fit_distr(Corr.effective_mass_t(pt2_etaS_OS_L.col(0)[i_ens], ""));
    //M
    Metas_M_distr = Corr.Fit_distr(Corr.effective_mass_t(pt2_etaS_M.col(0)[i_ens], ""));
    Metas_OS_M_distr = Corr.Fit_distr(Corr.effective_mass_t(pt2_etaS_OS_M.col(0)[i_ens], ""));


 



    //get MV from vector correlator
    //L
    Corr.Tmin = Tmin_VV; Corr.Tmax= Tmax_VV;
    MV_L_distr = Corr.Fit_distr( Corr.effective_mass_t(V_strange_L_distr, ""));
    overlap_V_L_distr = Za*Za*Corr.Fit_distr(  Corr.residue_t(V_strange_L_distr, ""))/(2.0*MV_L_distr);
    Corr.Tmin = Tmin_VV_OS; Corr.Tmax = Tmax_VV_OS;
    MV_OS_L_distr = Corr.Fit_distr( Corr.effective_mass_t(V_strange_OS_L_distr, ""));
    overlap_OS_V_L_distr = Zv*Zv*Corr.Fit_distr(  Corr.residue_t(V_strange_OS_L_distr, ""))/(2.0*MV_OS_L_distr);
    //M
    Corr.Tmin = Tmin_VV; Corr.Tmax = Tmax_VV;
    MV_M_distr = Corr.Fit_distr( Corr.effective_mass_t(V_strange_M_distr, ""));
    overlap_V_M_distr = Za*Za*Corr.Fit_distr(  Corr.residue_t(V_strange_M_distr, ""))/(2.0*MV_M_distr);
    Corr.Tmin = Tmin_VV_OS; Corr.Tmax = Tmax_VV_OS;
    MV_OS_M_distr = Corr.Fit_distr( Corr.effective_mass_t(V_strange_OS_M_distr, ""));
    overlap_OS_V_M_distr = Zv*Zv*Corr.Fit_distr(  Corr.residue_t(V_strange_OS_M_distr, ""))/(2.0*MV_OS_M_distr);


    vector<distr_t> Metas_vec({ Metas_L_distr/a_distr, Metas_M_distr/a_distr});
    vector<distr_t> Mphi_vec({ MV_L_distr/a_distr, MV_M_distr/a_distr});
    vector<distr_t> X_fit;
    if(Extrapolation_strange_mode == "phi") X_fit = Mphi_vec;
    else if(Extrapolation_strange_mode == "etas") X_fit = Metas_vec;
    else crash("Extrapolation strange mode: "+Extrapolation_strange_mode+" is invalid");
    distr_t X_phys = (Extrapolation_strange_mode=="phi")?Mphi_phys_distr:Metas_phys_distr;



    //############################################################################################
    //perturbative substraction of lattice artifacts##############################################

    //free corr LO artifacts
    //################## READ FREE THEORY VECTOR-VECTOR CORRELATOR OPPOSITE R ####################################
    string Pt_free_oppor_L= "../Vkvk_cont/"+to_string(Corr.Nt/2)+"_m"+to_string_with_precision(2*L_info.ms_L,3)+"/OPPOR";
    string Pt_free_oppor_M= "../Vkvk_cont/"+to_string(Corr.Nt/2)+"_m"+to_string_with_precision(2*L_info.ms_M,3)+"/OPPOR";
    Vfloat VV_free_oppor_L= Read_From_File(Pt_free_oppor_L, 1, 4);
    Vfloat VV_free_oppor_M= Read_From_File(Pt_free_oppor_M, 1, 4);
    if(VV_free_oppor_L.size() != Corr.Nt) crash("Failed to read properly free VV correlator mc_L w opposite r");
    if(VV_free_oppor_M.size() != Corr.Nt) crash("Failed to read properly free VV correlator mc_M w opposite r");
    //################## READ FREE THEORY VECTOR-VECTOR CORRELATOR SAME R ####################################
    string Pt_free_samer_L= "../Vkvk_cont/"+to_string(Corr.Nt/2)+"_m"+to_string_with_precision(2*L_info.ms_L,3)+"/SAMER";
    string Pt_free_samer_M= "../Vkvk_cont/"+to_string(Corr.Nt/2)+"_m"+to_string_with_precision(2*L_info.ms_M,3)+"/SAMER";
    Vfloat VV_free_samer_L= Read_From_File(Pt_free_samer_L, 1, 4);
    Vfloat VV_free_samer_M= Read_From_File(Pt_free_samer_M, 1, 4);
    if(VV_free_samer_L.size() != Corr.Nt) crash("Failed to read properly free VV correlator mc_L  w same r");
    if(VV_free_samer_M.size() != Corr.Nt) crash("Failed to read properly free VV correlator mc_M  w same r");
     
    //free corr LO artifacts
    Vfloat free_corr_log_art(Corr.Nt, 0.0);
    for(int t=0;t<Corr.Nt;t++) { if( t*a_distr.ave() < add_pert_corr_strange_up_to*fm_to_inv_Gev && t != 0) {   free_corr_log_art[t] = -1.0*pert_corr_strange_on_off*(qs*qs)*(1.0/(2.0*M_PI*M_PI*pow(t,5)));} else free_corr_log_art[t] = 0.0;

      if(t==0 || t*a_distr.ave() > add_pert_corr_strange_up_to*fm_to_inv_Gev) { VV_free_samer_L[t] =0; VV_free_samer_M[t] = 0;  VV_free_oppor_L[t] = 0; VV_free_oppor_M[t]=0;}


    }


    //####################################################################
    //Read perturbative spectral_density for OS and tm
    //####################################################################
    
    Vfloat Spec_tm_L = Read_From_File("../data/R_ratio/spec_dens_free/tm/am_"+to_string_with_precision(L_info.ms_L,5), 2, 4);
    Vfloat Spec_OS_L = Read_From_File("../data/R_ratio/spec_dens_free/OS/am_"+to_string_with_precision(L_info.ms_L,5), 2, 4);
    Vfloat Spec_tm_M = Read_From_File("../data/R_ratio/spec_dens_free/tm/am_"+to_string_with_precision(L_info.ms_M,5), 2, 4);
    Vfloat Spec_OS_M = Read_From_File("../data/R_ratio/spec_dens_free/OS/am_"+to_string_with_precision(L_info.ms_M,5), 2, 4);
    Vfloat Ergs_pert = Read_From_File("../data/R_ratio/spec_dens_free/tm/am_"+to_string_with_precision(L_info.ms_L,5), 1, 4);

    cout<<"perturbative spectral density for Ensemble: "<<V_strange_1_L.Tag[i_ens]<<" READ! "<<endl;
    

    //interpolate perturbative data
    boost::math::interpolators::cardinal_cubic_b_spline<double> F_boost_tm_L(Spec_tm_L.begin(), Spec_tm_L.end(), Ergs_pert[0], 2.0*Ergs_pert[0]);
    boost::math::interpolators::cardinal_cubic_b_spline<double> F_boost_OS_L(Spec_OS_L.begin(), Spec_OS_L.end(), Ergs_pert[0], 2.0*Ergs_pert[0]);
    boost::math::interpolators::cardinal_cubic_b_spline<double> F_boost_tm_M(Spec_tm_M.begin(), Spec_tm_M.end(), Ergs_pert[0], 2.0*Ergs_pert[0]);
    boost::math::interpolators::cardinal_cubic_b_spline<double> F_boost_OS_M(Spec_OS_M.begin(), Spec_OS_M.end(), Ergs_pert[0], 2.0*Ergs_pert[0]);
   
    cout<<"Cubic spline for perturbative spectral density for Ensemble: "<<V_strange_1_L.Tag[i_ens]<<" produced! "<<endl;

    auto F_free_tm_L = [&F_boost_tm_L](double E) { return F_boost_tm_L(E);};
    auto F_free_OS_L = [&F_boost_OS_L](double E) { return F_boost_OS_L(E);};
    auto F_free_tm_M = [&F_boost_tm_M](double E) { return F_boost_tm_M(E);};
    auto F_free_OS_M = [&F_boost_OS_M](double E) { return F_boost_OS_M(E);};
  

    


    //multiply corr using Zv and Za
    //L
    V_strange_L_distr = (V_strange_L_distr)*( 1.0 + pert_corr_strange_on_off*(1.0/(Za*Za*V_strange_L_distr))*VV_free_oppor_L) ;
    V_strange_OS_L_distr = (V_strange_OS_L_distr)*(1.0 + pert_corr_strange_on_off*(1.0/(Zv*Zv*V_strange_OS_L_distr))*VV_free_samer_L);
    //M
    V_strange_M_distr = (V_strange_M_distr)*(1.0 + pert_corr_strange_on_off*(1.0/(Za*Za*V_strange_M_distr))*VV_free_oppor_M) ;
    V_strange_OS_M_distr = (V_strange_OS_M_distr)*( 1.0 + pert_corr_strange_on_off*(1.0/(Zv*Zv*V_strange_OS_M_distr))*VV_free_samer_M);


    //print correlator

    //L
    Print_To_File({}, {V_strange_L_distr.ave(), V_strange_L_distr.err(), V_strange_OS_L_distr.ave(), V_strange_OS_L_distr.err()}, "../data/R_ratio/"+Tag_reco_type+"/corr/"+strange_tag+"/corr_L_"+V_strange_1_L.Tag[i_ens]+".dat", "", "# t  tm  OS");
    //M
    Print_To_File({}, {V_strange_M_distr.ave(), V_strange_M_distr.err(), V_strange_OS_M_distr.ave(), V_strange_OS_M_distr.err()}, "../data/R_ratio/"+Tag_reco_type+"/corr/"+strange_tag+"/corr_M_"+V_strange_1_L.Tag[i_ens]+".dat", "", "# t  tm  OS");

    //############################################################################################################

  


    //#########################   RECONSTRUCT THE SMEARED R-RATIO ################################


    //set tmax to the value where the error on V(t) is larger than x%
    //#############################################################################################################################

    //L
    bool Found_error_less_x_percent=false;
    double x=5;
    double tmax_L=1;
    while(!Found_error_less_x_percent && tmax_L < Corr.Nt/2) {
   
      if( (V_strange_L_distr.distr_list[tmax_L]).err()/fabs( (V_strange_L_distr.distr_list[tmax_L]).ave()) <  0.01*x) tmax_L++;
      else Found_error_less_x_percent=true;
    }

    Found_error_less_x_percent=false;
    double tmax_OS_L=1;
    while(!Found_error_less_x_percent && tmax_OS_L < Corr.Nt/2) {
   
      if( (V_strange_OS_L_distr.distr_list[tmax_OS_L]).err()/fabs( (V_strange_OS_L_distr.distr_list[tmax_OS_L]).ave()) <  0.01*x) tmax_OS_L++;
      else Found_error_less_x_percent=true;
    }

    //M
    Found_error_less_x_percent=false;
    double tmax_M=1;
    while(!Found_error_less_x_percent && tmax_M < Corr.Nt/2) {
   
      if( (V_strange_M_distr.distr_list[tmax_M]).err()/fabs( (V_strange_M_distr.distr_list[tmax_M]).ave()) <  0.01*x) tmax_M++;
      else Found_error_less_x_percent=true;
    }

    Found_error_less_x_percent=false;
    double tmax_OS_M=1;
    while(!Found_error_less_x_percent && tmax_OS_M < Corr.Nt/2) {
   
      if( (V_strange_OS_M_distr.distr_list[tmax_OS_M]).err()/fabs( (V_strange_OS_M_distr.distr_list[tmax_OS_M]).ave()) <  0.01*x) tmax_OS_M++;
      else Found_error_less_x_percent=true;
    }

    if(Use_t_up_to_T_half) tmax_L=tmax_OS_L= tmax_M=tmax_OS_M = T/2-1;
    
    //#############################################################################################################################


    //#####################################    MODEL ESTIMATE FOR SYSTEMATIC ERRORS ##############################################

    double resc_strange= rho_R*pow(qs,2);

    //tm
    auto model_strange_tm_L = [&a_distr,&MV_L_distr, &F_free_tm_L, &overlap_V_L_distr, &resc_strange](double E) {
      double DE= 0.003*a_distr.ave();
      bool Is_pert= (E >= 1.3*a_distr.ave());
      return resc_strange*(overlap_V_L_distr.ave()*((1.0/sqrt(2*M_PI*DE*DE))*exp(- pow(E-MV_L_distr.ave(),2)/(2*DE*DE))) + Is_pert*F_free_tm_L(E));
    };
    auto model_strange_tm_M = [&a_distr,&MV_M_distr, &F_free_tm_M, &overlap_V_M_distr, &resc_strange](double E) {
      double DE= 0.003*a_distr.ave();
      bool Is_pert= (E >= 1.3*a_distr.ave());
      return resc_strange*(overlap_V_M_distr.ave()*((1.0/sqrt(2*M_PI*DE*DE))*exp(- pow(E-MV_M_distr.ave(),2)/(2*DE*DE))) + Is_pert*F_free_tm_M(E));
    };
    //OS
    auto model_strange_OS_L = [&a_distr,&MV_OS_L_distr, &F_free_OS_L, &overlap_OS_V_L_distr, &resc_strange](double E) {
      double DE= 0.003*a_distr.ave();
      bool Is_pert= (E >= 1.3*a_distr.ave());
      return resc_strange*(overlap_OS_V_L_distr.ave()*((1.0/sqrt(2*M_PI*DE*DE))*exp(- pow(E-MV_OS_L_distr.ave(),2)/(2*DE*DE))) + Is_pert*F_free_OS_L(E));
    };
    auto model_strange_OS_M = [&a_distr,&MV_OS_M_distr, &F_free_OS_M, &overlap_OS_V_M_distr, &resc_strange](double E) {
      double DE= 0.003*a_distr.ave();
      bool Is_pert= (E >= 1.3*a_distr.ave());
      return resc_strange*(overlap_OS_V_M_distr.ave()*((1.0/sqrt(2*M_PI*DE*DE))*exp(- pow(E-MV_OS_M_distr.ave(),2)/(2*DE*DE))) + Is_pert*F_free_OS_M(E));
    };
  

    //tm
      auto f_syst_tm_L = [&overlap_V_L_distr, &MV_L_distr, &F_free_tm_L, &resc_strange, &a_distr](const function<double(double)> &F) ->double {
      double val_mod, err_mod;
      auto FS= [&F, &a_distr, &F_free_tm_L](double E) { if (E>= 1.3*a_distr.ave()) return F_free_tm_L(E)*F(E); return 0.0;}; 
      gsl_function_pp<decltype(FS)> SYST(FS);
      gsl_integration_workspace * w_SYST = gsl_integration_workspace_alloc (1000);
      gsl_function *G_SYST = static_cast<gsl_function*>(&SYST);
      gsl_integration_qags(G_SYST, Eth*a_distr.ave(), 4.0,  0.0, 3e-2, 1000, w_SYST, &val_mod, &err_mod);
      if(err_mod/fabs(val_mod) > 5e-2) crash("Cannot reach accuracy in evaluating systematic");
      return resc_strange*fabs((val_mod + F(MV_L_distr.ave())*overlap_V_L_distr.ave()));
      
    };

    auto f_syst_tm_M = [&overlap_V_M_distr, &MV_M_distr, &F_free_tm_M, &resc_strange, &a_distr](const function<double(double)> &F) ->double {
      double val_mod, err_mod;
      auto FS= [&F, &a_distr, &F_free_tm_M](double E) { if (E>= 1.3*a_distr.ave()) return F_free_tm_M(E)*F(E); return 0.0;}; 
      gsl_function_pp<decltype(FS)> SYST(FS);
      gsl_integration_workspace * w_SYST = gsl_integration_workspace_alloc (1000);
      gsl_function *G_SYST = static_cast<gsl_function*>(&SYST);
      gsl_integration_qags(G_SYST, Eth*a_distr.ave(), 4.0,  0.0, 3e-2, 1000, w_SYST, &val_mod, &err_mod);
      if(err_mod/fabs(val_mod) > 5e-2) crash("Cannot reach accuracy in evaluating systematic");
      return resc_strange*fabs((val_mod + F(MV_M_distr.ave())*overlap_V_M_distr.ave()));
      
    };

  
    //OS
    auto f_syst_OS_L = [&overlap_OS_V_L_distr, &MV_OS_L_distr, &F_free_OS_L, &resc_strange, &a_distr](const function<double(double)> &F) ->double {
      double val_mod, err_mod;
      auto FS= [&F, &a_distr, &F_free_OS_L](double E) { if (E>= 1.3*a_distr.ave()) return F_free_OS_L(E)*F(E); return 0.0;}; 
      gsl_function_pp<decltype(FS)> SYST(FS);
      gsl_integration_workspace * w_SYST = gsl_integration_workspace_alloc (1000);
      gsl_function *G_SYST = static_cast<gsl_function*>(&SYST);
      gsl_integration_qags(G_SYST, Eth*a_distr.ave(), 4.0,  0.0, 3e-2, 1000, w_SYST, &val_mod, &err_mod);
      if(err_mod/fabs(val_mod) > 5e-2) crash("Cannot reach accuracy in evaluating systematic");
      return resc_strange*fabs((val_mod + F(MV_OS_L_distr.ave())*overlap_OS_V_L_distr.ave()));
      
    };

    auto f_syst_OS_M = [&overlap_OS_V_M_distr, &MV_OS_M_distr, &F_free_OS_M, &resc_strange, &a_distr](const function<double(double)> &F) ->double {
      double val_mod, err_mod;
      auto FS= [&F, &a_distr, &F_free_OS_M](double E) { if (E> 1.3*a_distr.ave()) return F_free_OS_M(E)*F(E); return 0.0;}; 
      gsl_function_pp<decltype(FS)> SYST(FS);
      gsl_integration_workspace * w_SYST = gsl_integration_workspace_alloc (1000);
      gsl_function *G_SYST = static_cast<gsl_function*>(&SYST);
      gsl_integration_qags(G_SYST, Eth*a_distr.ave(), 4.0,  0.0, 3e-2, 1000, w_SYST, &val_mod, &err_mod);
      if(err_mod/fabs(val_mod) > 5e-2) crash("Cannot reach accuracy in evaluating systematic");
      return resc_strange*fabs((val_mod + F(MV_OS_M_distr.ave())*overlap_OS_V_M_distr.ave()));
      
    };

 
    //##############################################################################################################################




    

  

    //##################################################
    //#############################SET GENERAL PARAMETERS######################################
    double E0 =  Eth*a_distr.ave(); 
    double mult_SANF=1e-4;
    double mult_TANT=10.0;
    //#########################################################################################


   


    for(int isg=0;isg<(signed)sigmas.size();isg++) {


      double sigma=sigmas[isg]*a_distr.ave();
      cout<<"Reconstructing R(E)....."<<endl;
      cout<<"aE0: "<<E0<<" -> E0: "<<E0/a_distr.ave()<<" [GeV] "<<endl;
      cout<<"SM_TYPE: "<<SM_TYPE<<endl;
      cout<<"a*sigma(E*): "<<sigma<<" -> sigma(E*): "<<sigma/a_distr.ave()<<" [GeV] "<<endl;

  
      distr_t_list Spectral_dens_L(UseJack, Ergs_GeV_list.size());
      distr_t_list Spectral_dens_OS_L(UseJack, Ergs_GeV_list.size());
      distr_t_list Spectral_dens_SANF_L(UseJack, Ergs_GeV_list.size());
      distr_t_list Spectral_dens_OS_SANF_L(UseJack, Ergs_GeV_list.size());
      distr_t_list Spectral_dens_M(UseJack, Ergs_GeV_list.size());
      distr_t_list Spectral_dens_OS_M(UseJack, Ergs_GeV_list.size());
      distr_t_list Spectral_dens_SANF_M(UseJack, Ergs_GeV_list.size());
      distr_t_list Spectral_dens_OS_SANF_M(UseJack, Ergs_GeV_list.size());
      distr_t_list Spectral_dens_Extr(UseJack, Ergs_GeV_list.size());
      distr_t_list Spectral_dens_OS_Extr(UseJack, Ergs_GeV_list.size());
      distr_t_list Spectral_dens_SANF_Extr(UseJack, Ergs_GeV_list.size());
      distr_t_list Spectral_dens_OS_SANF_Extr(UseJack, Ergs_GeV_list.size());
      Vfloat syst_tm_L(Ergs_GeV_list.size());
      Vfloat syst_OS_L(Ergs_GeV_list.size());
      Vfloat syst_tm_SANF_L(Ergs_GeV_list.size());
      Vfloat syst_OS_SANF_L(Ergs_GeV_list.size());
      Vfloat syst_tm_M(Ergs_GeV_list.size());
      Vfloat syst_OS_M(Ergs_GeV_list.size());
      Vfloat syst_tm_SANF_M(Ergs_GeV_list.size());
      Vfloat syst_OS_SANF_M(Ergs_GeV_list.size());

      vector<tuple<int,double, double>> thread_times_tm(Ergs_GeV_list.size()), thread_times_OS(Ergs_GeV_list.size());


      //COMPUTE THE SMEARED R-RATIO

      cout<<"Looping over energies"<<flush;

      #pragma omp parallel for schedule(dynamic) 
      for(int ip=0; ip<(signed)Ergs_GeV_list.size();ip++) {

	double mean = Ergs_GeV_list[ip]*a_distr.ave();
	
	double lambda_Estar;
	double lambda_Estar_SANF;
 

	//L (T)
	//define jackknife distribution to account for systematic error:
	distr_t syst_L_T_tm(UseJack), syst_L_T_OS(UseJack);
	for(int ijack=0; ijack<Njacks;ijack++) {syst_L_T_tm.distr.push_back( GM()/sqrt(Njacks-1.0)); syst_L_T_OS.distr.push_back( GM()/sqrt(Njacks-1.0));}

	auto start = chrono::system_clock::now();
	Spectral_dens_L.distr_list[ip] = Get_Laplace_transfo(  mean,  sigma, E0,  T, tmax_L, prec, SM_TYPE+"_ov_E2",f, V_strange_L_distr, syst_tm_L[ip], mult_TANT, lambda_Estar, "TANT", "tm", "L_"+V_strange_1_L.Tag[i_ens], -1 , 0, rho_R*Za*Za*pow(qs,2), 0.0, "R_ratio_strange", cov_tm_L, f_syst_tm_L, 0, model_strange_tm_L, Is_Emax_Finite, Emax,beta  )+ syst_tm_L[ip]*syst_L_T_tm ;
	auto end = chrono::system_clock::now();
	chrono::duration<double> elapsed_seconds = end-start;
	double time_L_tm= elapsed_seconds.count();
	if(R_ratio_verbosity_lev) {
	  cout<<endl<<flush;
	  cout<<"Elapsed time[L_tm, #thread="<<omp_get_thread_num()<<"] : "<<time_L_tm<<" s"<<endl<<flush;
	}
	cout<<"."<<flush;

	start = chrono::system_clock::now();
	Spectral_dens_OS_L.distr_list[ip]=Get_Laplace_transfo(  mean,  sigma, E0,  T, tmax_OS_L, prec, SM_TYPE+"_ov_E2",f, V_strange_OS_L_distr, syst_OS_L[ip], mult_TANT, lambda_Estar, "TANT", "OS", "L_"+V_strange_1_L.Tag[i_ens], -1 , 0, rho_R*Zv*Zv*pow(qs,2), 0.0, "R_ratio_strange", cov_OS_L, f_syst_OS_L, 0, model_strange_OS_L, Is_Emax_Finite, Emax,beta  )+ syst_OS_L[ip]*syst_L_T_OS ;
	end = chrono::system_clock::now();
	elapsed_seconds = end-start;
	double time_L_OS= elapsed_seconds.count();
	if(R_ratio_verbosity_lev) {
	  cout<<endl<<flush;
	  cout<<"Elapsed time[L_OS, #thread="<<omp_get_thread_num()<<"] : "<<time_L_OS<<" s"<<endl<<flush;
	}
	cout<<"."<<flush;
	
	//L (S)
	if(!SANF_MODE_OFF) {
	  //define jackknife distribution to account for systematic error:
	  distr_t syst_L_S_tm(UseJack), syst_L_S_OS(UseJack);
	  for(int ijack=0; ijack<Njacks;ijack++) {syst_L_S_tm.distr.push_back( GM()/sqrt(Njacks-1.0)); syst_L_S_OS.distr.push_back( GM()/sqrt(Njacks-1.0));}
	  Spectral_dens_SANF_L.distr_list[ip] = Get_Laplace_transfo(  mean,  sigma, E0,  T, tmax_L, prec, SM_TYPE+"_ov_E2",f, V_strange_L_distr, syst_tm_SANF_L[ip], mult_SANF, lambda_Estar_SANF, "SANF", "tm", "L_"+V_strange_1_L.Tag[i_ens], -1 , 0, rho_R*Za*Za*pow(qs,2), 0.0, "R_ratio_strange", cov_tm_L, fake_func, 0, fake_func_d, Is_Emax_Finite, Emax,beta  )+ syst_tm_SANF_L[ip]*syst_L_S_tm ;
	  Spectral_dens_OS_SANF_L[ip] = Get_Laplace_transfo(  mean,  sigma, E0,  T, tmax_OS_L, prec, SM_TYPE+"_ov_E2",f, V_strange_OS_L_distr, syst_OS_SANF_L[ip], mult_SANF, lambda_Estar_SANF, "SANF", "OS", "L_"+V_strange_1_L.Tag[i_ens], -1 , 0, rho_R*Zv*Zv*pow(qs,2), 0.0, "R_ratio_strange", cov_OS_L, fake_func, 0, fake_func_d, Is_Emax_Finite, Emax,beta  )+ syst_OS_SANF_L[ip]*syst_L_S_OS ;
	}
	
	//M (T)
	//define jackknife distribution to account for systematic error:
	distr_t syst_M_T_tm(UseJack), syst_M_T_OS(UseJack);
	for(int ijack=0; ijack<Njacks;ijack++) {syst_M_T_tm.distr.push_back( GM()/sqrt(Njacks-1.0)); syst_M_T_OS.distr.push_back( GM()/sqrt(Njacks-1.0));}

	start= chrono::system_clock::now();
	Spectral_dens_M.distr_list[ip] = Get_Laplace_transfo(  mean,  sigma, E0,  T, tmax_M, prec, SM_TYPE+"_ov_E2",f, V_strange_M_distr, syst_tm_M[ip], mult_TANT, lambda_Estar, "TANT", "tm", "M_"+V_strange_1_M.Tag[i_ens], -1 , 0, rho_R*Za*Za*pow(qs,2), 0.0, "R_ratio_strange", cov_tm_M, f_syst_tm_M, 0, model_strange_tm_M, Is_Emax_Finite, Emax,beta  )+ syst_tm_M[ip]*syst_M_T_tm ;
	end = chrono::system_clock::now();
	elapsed_seconds = end-start;
	double time_M_tm= elapsed_seconds.count();
	if(R_ratio_verbosity_lev) {
	  cout<<endl<<flush;
	  cout<<"Elapsed time[M_tm, #thread="<<omp_get_thread_num()<<"] : "<<time_M_tm<<" s"<<endl<<flush;
	}
	cout<<"."<<flush;

	start = chrono::system_clock::now();
	Spectral_dens_OS_M.distr_list[ip] = Get_Laplace_transfo(  mean,  sigma, E0,  T, tmax_OS_M, prec, SM_TYPE+"_ov_E2",f, V_strange_OS_M_distr, syst_OS_M[ip], mult_TANT, lambda_Estar, "TANT", "OS", "M_"+V_strange_1_M.Tag[i_ens], -1 , 0, rho_R*Zv*Zv*pow(qs,2) , 0.0, "R_ratio_strange", cov_OS_M, f_syst_OS_M, 0, model_strange_OS_M, Is_Emax_Finite, Emax,beta )+ syst_OS_M[ip]*syst_M_T_OS ;
	end = chrono::system_clock::now();
	elapsed_seconds = end-start;
	double time_M_OS= elapsed_seconds.count();
	if(R_ratio_verbosity_lev) {
	  cout<<endl<<flush;
	  cout<<"Elapsed time[M_OS, #thread="<<omp_get_thread_num()<<"] : "<<time_M_OS<<" s"<<endl<<flush;
	}
	cout<<"."<<flush;

	//M (S)
	if(!SANF_MODE_OFF) {
	  //define jackknife distribution to account for systematic error:
	  distr_t syst_M_S_tm(UseJack), syst_M_S_OS(UseJack);
	  for(int ijack=0; ijack<Njacks;ijack++) {syst_M_S_tm.distr.push_back( GM()/sqrt(Njacks-1.0)); syst_M_S_OS.distr.push_back( GM()/sqrt(Njacks-1.0));}
	  Spectral_dens_SANF_M.distr_list[ip] = Get_Laplace_transfo(  mean,  sigma, E0,  T, tmax_M, prec, SM_TYPE+"_ov_E2",f, V_strange_M_distr, syst_tm_SANF_M[ip], mult_SANF, lambda_Estar_SANF, "SANF", "tm", "M_"+V_strange_1_M.Tag[i_ens], -1 , 0, rho_R*Za*Za*pow(qs,2), 0.0, "R_ratio_strange", cov_tm_M, fake_func, 0, fake_func_d, Is_Emax_Finite, Emax,beta  )+ syst_tm_SANF_M[ip]*syst_M_S_tm ;
	  Spectral_dens_OS_SANF_M.distr_list[ip] = Get_Laplace_transfo(  mean,  sigma, E0,  T, tmax_OS_M, prec, SM_TYPE+"_ov_E2",f, V_strange_OS_M_distr, syst_OS_SANF_M[ip], mult_SANF, lambda_Estar_SANF, "SANF", "OS", "M_"+V_strange_1_M.Tag[i_ens], -1 , 0, rho_R*Zv*Zv*pow(qs,2), 0.0, "R_ratio_strange", cov_OS_M, fake_func, 0, fake_func_d, Is_Emax_Finite, Emax,beta  )+ syst_OS_SANF_M[ip]*syst_M_S_OS ;
	}


	thread_times_tm[ip] = make_tuple(omp_get_thread_num(), time_L_tm, time_M_tm);
	thread_times_OS[ip] = make_tuple(omp_get_thread_num(), time_L_OS, time_M_OS);


	//Extrapolate to the physical strange mass
	vector<distr_t> Y_fit_tm = {Spectral_dens_L.distr_list[ip], Spectral_dens_M.distr_list[ip]};
	vector<distr_t> Y_fit_OS = {Spectral_dens_OS_L.distr_list[ip], Spectral_dens_OS_M.distr_list[ip]};
	vector<distr_t> Y_fit_SANF_tm, Y_fit_SANF_OS;
	if(!SANF_MODE_OFF) {
	  Y_fit_SANF_tm = {Spectral_dens_SANF_L.distr_list[ip], Spectral_dens_SANF_M.distr_list[ip]};
	  Y_fit_SANF_OS = {Spectral_dens_OS_SANF_L.distr_list[ip], Spectral_dens_OS_SANF_M.distr_list[ip]};
	}

	//Extrapolation to physical strange quark mass (T)
	Spectral_dens_Extr.distr_list[ip] =  Obs_extrapolation_meson_mass( Y_fit_tm, X_fit, X_phys ,  "../data/R_ratio/"+Tag_reco_type+"/"+strange_tag  , "R_extr_tm_TANT_"+V_strange_1_L.Tag[i_ens]+"_sigma_"+to_string_with_precision(sigmas[isg],3)+"_Estar_"+to_string_with_precision(Ergs_GeV_list[ip],3)+".dat",  UseJack, "SPLINE" );
	Spectral_dens_OS_Extr.distr_list[ip] = Obs_extrapolation_meson_mass( Y_fit_OS, X_fit, X_phys ,  "../data/R_ratio/"+Tag_reco_type+"/"+strange_tag  , "R_extr_OS_TANT_"+V_strange_1_L.Tag[i_ens]+"_sigma_"+to_string_with_precision(sigmas[isg],3)+"_Estar_"+to_string_with_precision(Ergs_GeV_list[ip],3)+".dat",  UseJack, "SPLINE" );

	//Extrapolation to physical strange quark mass (S)
	if(!SANF_MODE_OFF) {
	Spectral_dens_SANF_Extr.distr_list[ip] =   Obs_extrapolation_meson_mass( Y_fit_SANF_tm, X_fit, X_phys ,  "../data/R_ratio/"+Tag_reco_type+"/"+strange_tag  , "R_extr_tm_SANF_"+V_strange_1_L.Tag[i_ens]+"_sigma_"+to_string_with_precision(sigmas[isg],3)+"_Estar_"+to_string_with_precision(Ergs_GeV_list[ip],3)+".dat",  UseJack, "SPLINE" );
	Spectral_dens_OS_SANF_Extr.distr_list[ip] =  Obs_extrapolation_meson_mass( Y_fit_SANF_OS, X_fit, X_phys ,  "../data/R_ratio/"+Tag_reco_type+"/"+strange_tag  , "R_extr_OS_SANF_"+V_strange_1_L.Tag[i_ens]+"_sigma_"+to_string_with_precision(sigmas[isg],3)+"_Estar_"+to_string_with_precision(Ergs_GeV_list[ip],3)+".dat",  UseJack, "SPLINE" );
	}
	
	

        
	//############################################################################################



      }

      cout<<"done"<<endl<<flush;
      cout<<"Summary of performances: "<<endl<<flush;
      cout<<"Erg #thread  L    M"<<endl<<flush;
      cout<<"- - - - - - - - - - - - - - - - - - - - - "<<endl<<flush;
      for(int ip=0; ip < (signed)Ergs_GeV_list.size();ip++) {
	cout<<Ergs_GeV_list[ip]<<","<<get<0>(thread_times_tm[ip])<<": "<<get<1>(thread_times_tm[ip])<<" s, "<<get<2>(thread_times_tm[ip])<<" s"<<endl<<flush;
	cout<<Ergs_GeV_list[ip]<<","<<get<0>(thread_times_OS[ip])<<": "<<get<1>(thread_times_OS[ip])<<" s, "<<get<2>(thread_times_OS[ip])<<" s"<<endl<<flush;
	cout<<"- - - - - - - - - - - - - - - - - - - - "<<endl<<flush;
      }
     


      if(SANF_MODE_OFF) {
	Spectral_dens_SANF_Extr = Spectral_dens_Extr;	Spectral_dens_OS_SANF_Extr = Spectral_dens_OS_Extr;
	Spectral_dens_SANF_L = Spectral_dens_L; Spectral_dens_OS_SANF_L = Spectral_dens_OS_L;
	Spectral_dens_SANF_M = Spectral_dens_M; Spectral_dens_OS_SANF_M = Spectral_dens_OS_M;
      }


      RE_strange_TANT_tm[isg][i_ens] = Spectral_dens_Extr;
      RE_strange_TANT_OS[isg][i_ens] = Spectral_dens_OS_Extr;
      RE_strange_SANF_tm[isg][i_ens] = Spectral_dens_SANF_Extr;
      RE_strange_SANF_OS[isg][i_ens] = Spectral_dens_OS_SANF_Extr;
   

      cout<<"printing output strange for sigma: "<<sigmas[isg]<<" ..."<<flush;

      //print to file
      //L
      Print_To_File({}, {Ergs_GeV_list,  Spectral_dens_L.ave(), Spectral_dens_L.err(),  Spectral_dens_SANF_L.ave(), Spectral_dens_SANF_L.err(), Spectral_dens_OS_L.ave(), Spectral_dens_OS_L.err(), Spectral_dens_OS_SANF_L.ave(), Spectral_dens_OS_SANF_L.err()}, "../data/R_ratio/"+Tag_reco_type+"/"+strange_tag+"/L_"+SM_TYPE+"_"+V_strange_1_L.Tag[i_ens]+"_sigma_"+to_string_with_precision(sigmas[isg],3)+".dat", "", "#E*(GeV)   R(E)_tm[T]  R(E)_tm[S]   R(E)_OS[T]  R(E)_OS[S]");
      //M
      Print_To_File({}, {Ergs_GeV_list,  Spectral_dens_M.ave(), Spectral_dens_M.err(), Spectral_dens_SANF_M.ave(), Spectral_dens_SANF_M.err(), Spectral_dens_OS_M.ave(), Spectral_dens_OS_M.err(), Spectral_dens_OS_SANF_M.ave(), Spectral_dens_OS_SANF_M.err()}, "../data/R_ratio/"+Tag_reco_type+"/"+strange_tag+"/M_"+SM_TYPE+"_"+V_strange_1_M.Tag[i_ens]+"_sigma_"+to_string_with_precision(sigmas[isg],3)+".dat", "", "#E*(GeV)   R(E)_tm[T]  R(E)_tm[S]   R(E)_OS[T]  R(E)_OS[S]");
      //Extr
      Print_To_File({}, {Ergs_GeV_list,  Spectral_dens_Extr.ave(), Spectral_dens_Extr.err(), Spectral_dens_SANF_Extr.ave(), Spectral_dens_SANF_Extr.err(), Spectral_dens_OS_Extr.ave(), Spectral_dens_OS_Extr.err(), Spectral_dens_OS_SANF_Extr.ave(), Spectral_dens_OS_SANF_Extr.err() }, "../data/R_ratio/"+Tag_reco_type+"/"+strange_tag+"/Extr_"+SM_TYPE+"_"+V_strange_1_M.Tag[i_ens]+"_sigma_"+to_string_with_precision(sigmas[isg],3)+".dat", "", "#E*(GeV)   R(E)_tm[T]  R(E)_tm[S]   R(E)_OS[T]  R(E)_OS[S]");

      cout<<"done!"<<endl;
    }


  }

  cout<<"Calculation of strange contribution COMPLETED"<<endl;

 }


  //#############################################################################################
  //#############################################################################################
  //#############################################################################################
  //#############################################################################################
  
  //################## CALCULATION OF THE DISCONNECTED CONTRIBUTIONS ############################

  
  //#############################################################################################
  //#############################################################################################
  //#############################################################################################
  //#############################################################################################

if(!skip_disconnected) { 

  cout<<"STARTING COMPUTATION OF DISCONNECTED CONTRIBUTION:"<<endl;
  cout<<"PRECISION: "<<prec/ln2_10<<" digits."<<endl;

  for(int i_ens=0;i_ens<Nens_disco;i_ens++) {

    Disco_Ens_Tags.push_back(disco_light.Tag[i_ens]);

    string FLAV= "disco";

    string disco_tag= "disco";
    
    CorrAnalysis Corr(UseJack, Njacks,Nboots);
    CorrAnalysis Corr_block_1(0, disco_light.Nconfs[i_ens],Nboots, i_ens);
    Corr_block_1.Nt= disco_light.nrows[i_ens];
    Corr.Nt = disco_light.nrows[i_ens];
    int T = Corr.Nt;

    cout<<"Analyzing Ensemble: "<<disco_light.Tag[i_ens]<<endl;
  
    //get lattice spacing
    distr_t a_distr(UseJack);
    distr_t Zv(UseJack), Za(UseJack);
    if(disco_light.Tag[i_ens].substr(1,1)=="A") {a_distr=a_A; Zv= ZV_A; Za = ZA_A;}
    else if(disco_light.Tag[i_ens].substr(1,1)=="B") {a_distr=a_B; Zv= ZV_B; Za = ZA_B;}
    else if(disco_light.Tag[i_ens].substr(1,1)=="C") {a_distr=a_C; Zv= ZV_C; Za = ZA_C;}
    else if(disco_light.Tag[i_ens].substr(1,1)=="D") {a_distr=a_D; Zv= ZV_D; Za = ZA_D;}
    else crash("lattice spacing distribution for Ens: "+disco_light.Tag[i_ens]+" not found");


    LatticeInfo L_info;
    L_info.LatInfo_new_ens(disco_light.Tag[i_ens]);
 
    distr_t_list  light_light_distr, strange_strange_distr, charm_charm_distr;
    distr_t_list  light_strange_distr, light_charm_distr, strange_charm_distr;

    distr_t_list  light_light_bin_distr, strange_strange_bin_distr, charm_charm_bin_distr;
    distr_t_list  light_strange_bin_distr, light_charm_bin_distr, strange_charm_bin_distr;

    

    if( (disco_light.Tag[i_ens] != disco_strange.Tag[i_ens]) || (disco_light.Tag[i_ens] != disco_charm.Tag[i_ens]) ||  (disco_light.Tag[i_ens] != disco_light_strange.Tag[i_ens]) || (disco_light.Tag[i_ens] != disco_light_charm.Tag[i_ens]) || (disco_light.Tag[i_ens] != disco_strange_charm.Tag[i_ens])) crash("disconnected ensembles do not match"); 
    

    //diagonal disco
    light_light_distr = Corr.corr_t(disco_light.col(0)[i_ens], "");
    strange_strange_distr = Corr.corr_t(disco_strange.col(0)[i_ens], "");
    charm_charm_distr = Corr.corr_t(disco_charm.col(0)[i_ens], "");

    light_light_bin_distr = Corr_block_1.corr_t(disco_light.col(0)[i_ens], "");
    strange_strange_bin_distr = Corr_block_1.corr_t(disco_strange.col(0)[i_ens], "");
    charm_charm_bin_distr = Corr_block_1.corr_t(disco_charm.col(0)[i_ens], "");
    //off-diagonal disco
    light_strange_distr = Corr.corr_t(disco_light_strange.col(0)[i_ens], "");
    light_charm_distr = Corr.corr_t(disco_light_charm.col(0)[i_ens], "");
    strange_charm_distr = Corr.corr_t(disco_strange_charm.col(0)[i_ens], "");

    light_strange_bin_distr = Corr_block_1.corr_t(disco_light_strange.col(0)[i_ens], "");
    light_charm_bin_distr = Corr_block_1.corr_t(disco_light_charm.col(0)[i_ens], "");
    strange_charm_bin_distr = Corr_block_1.corr_t(disco_strange_charm.col(0)[i_ens], "");


    distr_t_list disco_distr =  pow(qu+qd,2)*light_light_distr + pow(qs,2)*strange_strange_distr + pow(qc,2)*charm_charm_distr + 2.0*(qu+qd)*qs*light_strange_distr + 2.0*(qu+qd)*qc*light_charm_distr + 2.0*qs*qc*strange_charm_distr;
    distr_t_list disco_bin_distr =  pow(qu+qd,2)*light_light_bin_distr + pow(qs,2)*strange_strange_bin_distr + pow(qc,2)*charm_charm_bin_distr + 2.0*(qu+qd)*qs*light_strange_bin_distr + 2.0*(qu+qd)*qc*light_charm_bin_distr + 2.0*qs*qc*strange_charm_bin_distr;


    
    //print covariance matrix

    Vfloat cov_OS, TT, RR;
    Vfloat corr_OS;
  
    for(int tt=0;tt<Corr.Nt;tt++)
      for(int rr=0;rr<Corr.Nt;rr++) {
	TT.push_back(tt);
	RR.push_back(rr);
	cov_OS.push_back( pow(disco_bin_distr.ave(1)/disco_bin_distr.ave(0),2)*(disco_bin_distr.distr_list[tt]%disco_bin_distr.distr_list[rr]));
	corr_OS.push_back( (disco_bin_distr.distr_list[tt]%disco_bin_distr.distr_list[rr])/( disco_bin_distr.err(tt)*disco_bin_distr.err(rr)));
      }

    Print_To_File({}, {TT,RR,cov_OS, corr_OS}, "../data/R_ratio/"+Tag_reco_type+"/covariance/"+disco_tag+"/cov_"+disco_light.Tag[i_ens]+"_OS.dat", "", "");

    



    //#########################   RECONSTRUCT THE SMEARED R-RATIO ################################


    //set tmax to the value where the error on V(t) is larger than x%
    //#############################################################################################################################

    bool Found_error_less_x_percent=false;
    double x=5;
    double tmax=1;
    while(!Found_error_less_x_percent && tmax < Corr.Nt/2) {
   
      if( (disco_distr.distr_list[tmax]).err()/fabs( (disco_distr.distr_list[tmax]).ave()) <  0.01*x) tmax++;
      else Found_error_less_x_percent=true;
    }

    if(Use_t_up_to_T_half) tmax=T/2-1;
    
    //#############################################################################################################################

  

    //##################################################
    //#############################SET GENERAL PARAMETERS######################################
    double E0 = Eth*a_distr.ave();
    double mult_TANT= 10.0;
    double mult_SANF = 1e-4;
    //#########################################################################################
    
  
    
    
    for(int isg=0;isg<(signed)sigmas.size();isg++) {

      double sigma=sigmas[isg]*a_distr.ave();
        
      cout<<"Reconstructing R(E)....."<<endl;
      cout<<"aE0: "<<E0<<" -> E0: "<<E0/a_distr.ave()<<" [GeV] "<<endl;
      cout<<"SM_TYPE: "<<SM_TYPE<<endl;
      cout<<"a*sigma(E*): "<<sigma<<" -> sigma(E*): "<<sigma/a_distr.ave()<<" [GeV] "<<endl;

      distr_t_list Spectral_dens(UseJack, Ergs_GeV_list.size());
      distr_t_list Spectral_dens_SANF(UseJack, Ergs_GeV_list.size());
      Vfloat syst(Ergs_GeV_list.size());
      Vfloat syst_SANF(Ergs_GeV_list.size());

      vector<tuple<int,double>> thread_times_OS(Ergs_GeV_list.size());
    

      cout<<"Looping over energies"<<flush;
      //COMPUTE THE SMEARED R-RATIO
      #pragma omp parallel for schedule(dynamic)
      for(int ip=0; ip<(signed)Ergs_GeV_list.size();ip++) {
      
	double mean = Ergs_GeV_list[ip]*a_distr.ave();
      
	double lambda_Estar;
	double lambda_Estar_SANF;
      

	//(T)

	auto start= chrono::system_clock::now();
	Spectral_dens.distr_list[ip] = Get_Laplace_transfo(  mean,  sigma, E0,  T, tmax, prec, SM_TYPE+"_ov_E2",f, disco_distr, syst[ip], mult_TANT, lambda_Estar, "TANT", "OS", disco_light.Tag[i_ens], -1 , 0, rho_R*Zv*Zv, 0.0, "R_ratio_disco", cov_OS, fake_func, 0, fake_func_d, Is_Emax_Finite, Emax,beta  ) ;
	auto end = chrono::system_clock::now();
	chrono::duration<double> elapsed_seconds = end-start;
	double time_OS= elapsed_seconds.count();
	if(R_ratio_verbosity_lev) {
	  cout<<endl<<flush;
	  cout<<"Elapsed time[OS, #thread="<<omp_get_thread_num()<<"] : "<<time_OS<<" s"<<endl<<flush;
	}
	cout<<"."<<flush;

	//(S)
	if(!SANF_MODE_OFF) {
	  Spectral_dens_SANF.distr_list[ip] = Get_Laplace_transfo(  mean,  sigma, E0,  T, tmax, prec, SM_TYPE+"_ov_E2",f, disco_distr, syst_SANF[ip], mult_SANF, lambda_Estar_SANF, "SANF", "OS", disco_light.Tag[i_ens], -1  , 0, rho_R*Zv*Zv, 0.0, "R_ratio_disco", cov_OS, fake_func, 0, fake_func_d, Is_Emax_Finite, Emax,beta ) ;
	}
	         
	//############################################################################################


	thread_times_OS[ip] = make_tuple(omp_get_thread_num(), time_OS);
      }
      cout<<"done!"<<endl<<flush;
      cout<<"Summary of performances: "<<endl<<flush;
      cout<<"Erg #thread  L    "<<endl<<flush;
      cout<<"- - - - - - - - - - - - - - - - - - - - - "<<endl<<flush;
      for(int ip=0; ip < (signed)Ergs_GeV_list.size();ip++) {
	cout<<Ergs_GeV_list[ip]<<","<<get<0>(thread_times_OS[ip])<<": "<<get<1>(thread_times_OS[ip])<<" s"<<endl<<flush;
	cout<<"- - - - - - - - - - - - - - - - - - - - "<<endl<<flush;
      }
   
      if(SANF_MODE_OFF) Spectral_dens_SANF = Spectral_dens;


      RE_disco_TANT[isg][i_ens] = Spectral_dens;
      RE_disco_SANF[isg][i_ens] = Spectral_dens_SANF;
    

      cout<<"printing output light for sigma: "<<sigmas[isg]<<" ..."<<flush;

      //light
      //print to file
      Print_To_File({}, {Ergs_GeV_list, Spectral_dens.ave(), Spectral_dens.err(), Spectral_dens_SANF.ave(), Spectral_dens_SANF.err()}, "../data/R_ratio/"+Tag_reco_type+"/"+disco_tag+"/"+SM_TYPE+"_"+disco_light.Tag[i_ens]+"_sigma_"+to_string_with_precision(sigmas[isg],3)+".dat", "", "# aE* E*(GeV)   R(E)[T]  R(E)[S]");

      cout<<"done!"<<endl<<flush;

    }
   }
 }


 

 //######################################################################################################################################################################
 //STORE JACKKNIFE DISTRIBUTIONS
 if(!test_mode) {   
 //light
 if(!skip_light) {
 for(int i_ens=0; i_ens<Nens_light;i_ens++) {
   boost::filesystem::create_directory("../data/R_ratio/"+Tag_reco_type+"/light/jackknife/tm/"+V_light_1.Tag[i_ens]);
   boost::filesystem::create_directory("../data/R_ratio/"+Tag_reco_type+"/light/jackknife/OS/"+V_light_1.Tag[i_ens]);
   for(int is=0; is<(signed)sigmas.size();is++) {
     boost::filesystem::create_directory("../data/R_ratio/"+Tag_reco_type+"/light/jackknife/tm/"+V_light_1.Tag[i_ens]+"/"+to_string_with_precision(sigmas[is],3));
     boost::filesystem::create_directory("../data/R_ratio/"+Tag_reco_type+"/light/jackknife/OS/"+V_light_1.Tag[i_ens]+"/"+to_string_with_precision(sigmas[is],3));
     for(int id_erg=0; id_erg<(signed)Ergs_GeV_list.size();id_erg++) {
       //print jackknife distribution for tm and OS
       ofstream Print_tm("../data/R_ratio/"+Tag_reco_type+"/light/jackknife/tm/"+V_light_1.Tag[i_ens]+"/"+to_string_with_precision(sigmas[is],3)+"/Erg_"+to_string_with_precision(Ergs_GeV_list[id_erg], 3)+".jack");
       ofstream Print_OS("../data/R_ratio/"+Tag_reco_type+"/light/jackknife/OS/"+V_light_1.Tag[i_ens]+"/"+to_string_with_precision(sigmas[is],3)+"/Erg_"+to_string_with_precision(Ergs_GeV_list[id_erg], 3)+".jack");
       for(int ijack=0; ijack<Njacks;ijack++) {
	 Print_tm<<RE_light_TANT_tm[is][i_ens].distr_list[id_erg].distr[ijack]<<endl;
	 Print_OS<<RE_light_TANT_OS[is][i_ens].distr_list[id_erg].distr[ijack]<<endl;     
       }
       Print_tm.close();
       Print_OS.close();
     }
   }
 }
 }

 //strange
 if(!skip_strange) {
 for(int i_ens=0; i_ens<Nens_strange;i_ens++) {
   boost::filesystem::create_directory("../data/R_ratio/"+Tag_reco_type+"/strange/jackknife/tm/"+V_strange_1_L.Tag[i_ens]);
   boost::filesystem::create_directory("../data/R_ratio/"+Tag_reco_type+"/strange/jackknife/OS/"+V_strange_1_L.Tag[i_ens]);
   for(int is=0; is<(signed)sigmas.size();is++) {
     boost::filesystem::create_directory("../data/R_ratio/"+Tag_reco_type+"/strange/jackknife/tm/"+V_strange_1_L.Tag[i_ens]+"/"+to_string_with_precision(sigmas[is],3));
     boost::filesystem::create_directory("../data/R_ratio/"+Tag_reco_type+"/strange/jackknife/OS/"+V_strange_1_L.Tag[i_ens]+"/"+to_string_with_precision(sigmas[is],3));
     for(int id_erg=0; id_erg<(signed)Ergs_GeV_list.size();id_erg++) {
       //print jackknife distribution for tm and OS
       ofstream Print_tm("../data/R_ratio/"+Tag_reco_type+"/strange/jackknife/tm/"+V_strange_1_L.Tag[i_ens]+"/"+to_string_with_precision(sigmas[is],3)+"/Erg_"+to_string_with_precision(Ergs_GeV_list[id_erg], 3)+".jack");
       ofstream Print_OS("../data/R_ratio/"+Tag_reco_type+"/strange/jackknife/OS/"+V_strange_1_L.Tag[i_ens]+"/"+to_string_with_precision(sigmas[is],3)+"/Erg_"+to_string_with_precision(Ergs_GeV_list[id_erg], 3)+".jack");
       for(int ijack=0; ijack<Njacks;ijack++) {
	 Print_tm<<RE_strange_TANT_tm[is][i_ens].distr_list[id_erg].distr[ijack]<<endl;
	 Print_OS<<RE_strange_TANT_OS[is][i_ens].distr_list[id_erg].distr[ijack]<<endl;     
       }
       Print_tm.close();
       Print_OS.close();
     }
   }
 }
 }



 //charm
 if(!skip_charm) {
 for(int i_ens=0; i_ens<Nens_charm;i_ens++) {
   boost::filesystem::create_directory("../data/R_ratio/"+Tag_reco_type+"/charm/jackknife/tm/"+V_charm_1_L.Tag[i_ens]);
   boost::filesystem::create_directory("../data/R_ratio/"+Tag_reco_type+"/charm/jackknife/OS/"+V_charm_1_L.Tag[i_ens]);
   for(int is=0; is<(signed)sigmas.size();is++) {
     boost::filesystem::create_directory("../data/R_ratio/"+Tag_reco_type+"/charm/jackknife/tm/"+V_charm_1_L.Tag[i_ens]+"/"+to_string_with_precision(sigmas[is],3));
     boost::filesystem::create_directory("../data/R_ratio/"+Tag_reco_type+"/charm/jackknife/OS/"+V_charm_1_L.Tag[i_ens]+"/"+to_string_with_precision(sigmas[is],3));
     for(int id_erg=0; id_erg<(signed)Ergs_GeV_list.size();id_erg++) {
       //print jackknife distribution for tm and OS
       ofstream Print_tm("../data/R_ratio/"+Tag_reco_type+"/charm/jackknife/tm/"+V_charm_1_L.Tag[i_ens]+"/"+to_string_with_precision(sigmas[is],3)+"/Erg_"+to_string_with_precision(Ergs_GeV_list[id_erg], 3)+".jack");
       ofstream Print_OS("../data/R_ratio/"+Tag_reco_type+"/charm/jackknife/OS/"+V_charm_1_L.Tag[i_ens]+"/"+to_string_with_precision(sigmas[is],3)+"/Erg_"+to_string_with_precision(Ergs_GeV_list[id_erg], 3)+".jack");
       for(int ijack=0; ijack<Njacks;ijack++) {
	 Print_tm<<RE_charm_TANT_tm[is][i_ens].distr_list[id_erg].distr[ijack]<<endl;
	 Print_OS<<RE_charm_TANT_OS[is][i_ens].distr_list[id_erg].distr[ijack]<<endl;     
       }
       Print_tm.close();
       Print_OS.close();
     }
   }
 }
 }


 //disconected
 if(!skip_disconnected) {
 for(int i_ens=0; i_ens<Nens_disco;i_ens++) {
   boost::filesystem::create_directory("../data/R_ratio/"+Tag_reco_type+"/disco/jackknife/"+disco_light.Tag[i_ens]);
   for(int is=0; is<(signed)sigmas.size();is++) {
     boost::filesystem::create_directory("../data/R_ratio/"+Tag_reco_type+"/disco/jackknife/"+disco_light.Tag[i_ens]+"/"+to_string_with_precision(sigmas[is],3));
     for(int id_erg=0; id_erg<(signed)Ergs_GeV_list.size();id_erg++) {
       //print jackknife distribution for tm and OS
       ofstream Print_OS("../data/R_ratio/"+Tag_reco_type+"/disco/jackknife/"+disco_light.Tag[i_ens]+"/"+to_string_with_precision(sigmas[is],3)+"/Erg_"+to_string_with_precision(Ergs_GeV_list[id_erg], 3)+".jack");
       for(int ijack=0; ijack<Njacks;ijack++) {
	 Print_OS<<RE_disco_TANT[is][i_ens].distr_list[id_erg].distr[ijack]<<endl;     
       }
       Print_OS.close();
     }
   }
 }
 }
 
 
 

 cout<<"Jackknife distributions stored!!"<<endl;
 }
 cout<<"Bye!"<<endl;
   
    
  
 return;
}




void R_ratio_cont_extrapolation() {

  int NUMM_THREADS= omp_get_max_threads();
  omp_set_num_threads(1);


  Get_lat_to_print();

  cout<<"##############CONTINUUM EXTRAPOLATION OF VARIOUS FLAVOUR CONTRIBUTION TO SMEARED R-RATIO#############"<<endl;


  
  

  
  //############################################################################################
  //generate fake jack_distr for lattice spacing a_A a_B, a_C, a_D and RENORMALIZATION CONSTANT
  GaussianMersenne GM(981832);
  LatticeInfo a_info;
  distr_t a_A(UseJack), a_B(UseJack), a_C(UseJack), a_D(UseJack);
  distr_t ZV_A(UseJack), ZV_B(UseJack), ZV_C(UseJack), ZV_D(UseJack);
  distr_t ZA_A(UseJack), ZA_B(UseJack), ZA_C(UseJack), ZA_D(UseJack);
  double a_A_ave, a_A_err, a_B_ave, a_B_err, a_C_ave, a_C_err, a_D_ave, a_D_err;
  double ZV_A_ave, ZV_A_err, ZV_B_ave, ZV_B_err, ZV_C_ave, ZV_C_err, ZV_D_ave, ZV_D_err;
  double ZA_A_ave, ZA_A_err, ZA_B_ave, ZA_B_err, ZA_C_ave, ZA_C_err, ZA_D_ave, ZA_D_err;
  a_info.LatInfo_new_ens("cA211a.53.24");
  a_A_ave= a_info.a_from_afp;
  a_A_err= a_info.a_from_afp_err;
  ZA_A_ave = a_info.Za_WI_strange;
  ZA_A_err = a_info.Za_WI_strange_err;
  ZV_A_ave = a_info.Zv_WI_strange;
  ZV_A_err = a_info.Zv_WI_strange_err;
  a_info.LatInfo_new_ens("cB211b.072.64");
  a_B_ave= a_info.a_from_afp;
  a_B_err= a_info.a_from_afp_err;
  ZA_B_ave = a_info.Za_WI_strange;
  ZA_B_err = a_info.Za_WI_strange_err;
  ZV_B_ave = a_info.Zv_WI_strange;
  ZV_B_err = a_info.Zv_WI_strange_err;
  a_info.LatInfo_new_ens("cC211a.06.80");
  a_C_ave= a_info.a_from_afp;
  a_C_err= a_info.a_from_afp_err;
  ZA_C_ave = a_info.Za_WI_strange;
  ZA_C_err = a_info.Za_WI_strange_err;
  ZV_C_ave = a_info.Zv_WI_strange;
  ZV_C_err = a_info.Zv_WI_strange_err;
  a_info.LatInfo_new_ens("cD211a.054.96");
  a_D_ave= a_info.a_from_afp;
  a_D_err= a_info.a_from_afp_err;
  ZA_D_ave = a_info.Za_WI_strange;
  ZA_D_err = a_info.Za_WI_strange_err;
  ZV_D_ave = a_info.Zv_WI_strange;
  ZV_D_err = a_info.Zv_WI_strange_err;
  
  if(UseJack)  { for(int ijack=0;ijack<Njacks;ijack++) {
      a_A.distr.push_back( fm_to_inv_Gev*( a_A_ave + GM()*a_A_err*(1.0/sqrt(Njacks-1.0))));
      a_B.distr.push_back( fm_to_inv_Gev*( a_B_ave + GM()*a_B_err*(1.0/sqrt(Njacks-1.0))));
      a_C.distr.push_back( fm_to_inv_Gev*( a_C_ave + GM()*a_C_err*(1.0/sqrt(Njacks-1.0))));
      a_D.distr.push_back( fm_to_inv_Gev*( a_D_ave + GM()*a_D_err*(1.0/sqrt(Njacks-1.0))));
      ZA_A.distr.push_back(  ZA_A_ave + GM()*ZA_A_err*(1.0/sqrt(Njacks -1.0)));
      ZV_A.distr.push_back(  ZV_A_ave + GM()*ZV_A_err*(1.0/sqrt(Njacks -1.0)));
      ZA_B.distr.push_back(  ZA_B_ave + GM()*ZA_B_err*(1.0/sqrt(Njacks -1.0)));
      ZV_B.distr.push_back(  ZV_B_ave + GM()*ZV_B_err*(1.0/sqrt(Njacks -1.0)));
      ZA_C.distr.push_back(  ZA_C_ave + GM()*ZA_C_err*(1.0/sqrt(Njacks -1.0)));
      ZV_C.distr.push_back(  ZV_C_ave + GM()*ZV_C_err*(1.0/sqrt(Njacks -1.0)));
      ZA_D.distr.push_back(  ZA_D_ave + GM()*ZA_D_err*(1.0/sqrt(Njacks -1.0)));
      ZV_D.distr.push_back(  ZV_D_ave + GM()*ZV_D_err*(1.0/sqrt(Njacks -1.0)));
      
    }
  }
  else {
    for (int iboot=0; iboot<Nboots;iboot++) {
      a_A.distr.push_back( fm_to_inv_Gev*( a_A_ave + GM()*a_A_err));
      a_B.distr.push_back( fm_to_inv_Gev*( a_B_ave + GM()*a_B_err));
      a_C.distr.push_back( fm_to_inv_Gev*( a_C_ave + GM()*a_C_err));
      a_D.distr.push_back( fm_to_inv_Gev*( a_D_ave + GM()*a_D_err));
      ZA_A.distr.push_back(  ZA_A_ave + GM()*ZA_A_err);
      ZV_A.distr.push_back(  ZV_A_ave + GM()*ZV_A_err);
      ZA_B.distr.push_back(  ZA_B_ave + GM()*ZA_B_err);
      ZV_B.distr.push_back(  ZV_B_ave + GM()*ZV_B_err);
      ZA_C.distr.push_back(  ZA_C_ave + GM()*ZA_C_err);
      ZV_C.distr.push_back(  ZV_C_ave + GM()*ZV_C_err);
      ZA_D.distr.push_back(  ZA_D_ave + GM()*ZA_D_err);
      ZV_D.distr.push_back(  ZV_D_ave + GM()*ZV_D_err);
      
    }
  }



  //############################################################################################



  //##### TMinuit2 classes for bootstrap fit ###################

  class ipar_R {

  public:
    
    ipar_R() {}


    double meas, err, a;
    bool Is_tm;
    
  };


  class fpar_R {

  public:

    fpar_R(Vfloat par) {
      if(par.size() != 3) crash("fpar_R constructor called with Vfloat of size != 3");
      D=par[0];
      D2_tm=par[1];
      D2_OS=par[2];
    }

    double D, D2_tm, D2_OS;

  };
  
  //###########################################################





  //#############    Init bootstrap fit     ###################

  bootstrap_fit<fpar_R, ipar_R> bf_R(Njacks);
  bf_R.set_warmup_lev(0);
  
  bf_R.Set_verbosity(1);
  bf_R.Add_par("D", 1.0, 0.01);
  bf_R.Add_par("D2_tm", 1.0, 0.01);
  bf_R.Add_par("D2_OS", 1.0, 0.01);
  //fit on mean values to get ch2
  bootstrap_fit<fpar_R, ipar_R> bf_R_ch2(1);
  bf_R_ch2.set_warmup_lev(0);
  bf_R_ch2.Set_verbosity(1);
  bf_R_ch2.Add_par("D", 1.0, 0.01);
  bf_R_ch2.Add_par("D2_tm", 1.0, 0.01);
  bf_R_ch2.Add_par("D2_OS", 1.0, 0.01);



  //ansatz
  bf_R.ansatz= [](const fpar_R &p, const ipar_R &ip) {
    double D2= (ip.Is_tm==1)?p.D2_tm:p.D2_OS;
    return p.D + D2*pow(ip.a*Lambda_QCD,2);
  };
  bf_R.measurement= [](const fpar_R &p, const ipar_R &ip) {

    return ip.meas;
  };
  bf_R.error= [](const fpar_R &p, const ipar_R &ip) {

    return ip.err;
  };
  
  bf_R_ch2.ansatz= bf_R.ansatz;
  bf_R_ch2.measurement= bf_R.measurement;
  bf_R_ch2.error= bf_R.error;

  
  //#########################################################
    

  
  


  //loop over different betas analyzed
  Vfloat betas({ 0.0, 1.0, 1.99, 2.99, 3.99, 0.0, 1.0, 1.99});
  Vfloat Emax_list({ 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0});
  vector<bool> Are_Emax_Finite({1,1,1,1,1, 0,0,0});
  int N= betas.size();

 

  
  for(int i=0; i < N;i++) {


    double beta= betas[i];
    double Is_Emax_Finite= Are_Emax_Finite[i];
    double Emax= Emax_list[i];
    string Tag_reco_type="Beta_"+to_string_with_precision(beta,2);
    Tag_reco_type+="_Emax_"+((Is_Emax_Finite==0)?"inf":to_string_with_precision(Emax,1));


    //create output directories
    boost::filesystem::create_directory("../data/R_ratio/"+Tag_reco_type+"/continuum");
    boost::filesystem::create_directory("../data/R_ratio/"+Tag_reco_type+"/continuum/total");
  
    vector<string> flavors({"light", "strange", "charm", "disco"});
    for(auto &f: flavors) {
    boost::filesystem::create_directory("../data/R_ratio/"+Tag_reco_type+"/continuum/"+f);
    for(auto &sigma:sigmas) {
      boost::filesystem::create_directory("../data/R_ratio/"+Tag_reco_type+"/continuum/"+f+"/"+to_string_with_precision(sigma,3));
    }
    }
  for(auto &sigma:sigmas) boost::filesystem::create_directory("../data/R_ratio/"+Tag_reco_type+"/continuum/total/"+to_string_with_precision(sigma,3));

  

  vector<vector<distr_t_list>> R_ratio_flav;
  R_ratio_flav.resize(flavors.size());
  for(auto &R_ratio_flav_sigma: R_ratio_flav) R_ratio_flav_sigma.resize(sigmas.size());

    int count_flav=0;

    //loop over contributions
    for( auto &flav: flavors) {


      distr_t_list a_distr_list(UseJack);
      vector<string> Ensemble_list;
      if(flav == "disco" || flav=="charm") {
	a_distr_list.distr_list.push_back(a_B);
	a_distr_list.distr_list.push_back(a_C);
	a_distr_list.distr_list.push_back(a_D);
	Ensemble_list = {"cB211b.072.64", "cC211a.06.80", "cD211a.054.96"};

      }
      else {
	a_distr_list.distr_list.push_back(a_B);
	a_distr_list.distr_list.push_back(a_B);
	a_distr_list.distr_list.push_back(a_C);
	a_distr_list.distr_list.push_back(a_D);
	Ensemble_list = {"cB211b.072.64", "cB211b.072.96", "cC211a.06.80", "cD211a.054.96"};
      }


      int combined_mult=(flav != "disco")?2:1;

      //fix D2_tm if disco
      if(combined_mult==1) {bf_R.Fix_par("D2_tm", 0.0); bf_R_ch2.Fix_par("D2_tm", 0.0);}
      else { bf_R.Release_par("D2_tm"); bf_R_ch2.Release_par("D2_tm");}

      
      //loop over sigma and Energies
      for(int is=0;is<(signed)sigmas.size();is++) {


	double sigma=sigmas[is];

	distr_t_list R_ratio_fixed_sigma(UseJack);
	Vfloat ch2_list;
	
	for(int id_erg=0;id_erg<(signed)Ergs_GeV_list.size();id_erg++) {
	  
	  double Erg=Ergs_GeV_list[id_erg];


	  //###################################################################
	  //load tm and OS jack-data corresponding to given sigma and Energy for all ensembles
	  distr_t_list data_tm(UseJack), data_OS(UseJack);
	  for(int iens=0; iens<(signed)Ensemble_list.size();iens++) {
	    string TAG_TM= "../data/R_ratio/"+Tag_reco_type+"/"+flav+"/jackknife/tm/"+Ensemble_list[iens]+"/"+to_string_with_precision(sigmas[is],3)+"/Erg_"+to_string_with_precision(Ergs_GeV_list[id_erg],3)+".jack";
	    string TAG_OS= "../data/R_ratio/"+Tag_reco_type+"/"+flav+"/jackknife/OS/"+Ensemble_list[iens]+"/"+to_string_with_precision(sigmas[is],3)+"/Erg_"+to_string_with_precision(Ergs_GeV_list[id_erg],3)+".jack";
	    if(flav=="disco") {
	      TAG_OS=  "../data/R_ratio/"+Tag_reco_type+"/"+flav+"/jackknife/"+Ensemble_list[iens]+"/"+to_string_with_precision(sigmas[is],3)+"/Erg_"+to_string_with_precision(Ergs_GeV_list[id_erg],3)+".jack";
	      TAG_TM=TAG_OS;
	    }
	    distr_t tm_distr(UseJack, Read_From_File(TAG_TM, 0,1));
	    distr_t OS_distr(UseJack, Read_From_File(TAG_OS, 0,1));

	    if(tm_distr.size() != Njacks) crash("tm distr size is different from Njacks="+to_string(Njacks));
	    if(OS_distr.size() != Njacks) crash("tm distr size is different from Njacks="+to_string(Njacks));
	    //push_back to data_tm and data_OS
	    data_tm.distr_list.push_back(tm_distr);
	    data_OS.distr_list.push_back(OS_distr);
	  }
	  //##################################################################

	  	 
	  bf_R.Set_number_of_measurements(combined_mult*Ensemble_list.size());
	  bf_R_ch2.Set_number_of_measurements(combined_mult*Ensemble_list.size());

	  //#################################################################
	  //Add covariance matrix if combined fit
	  if(flav != "disco") {
	    Eigen::MatrixXd Cov_Matrix(combined_mult*Ensemble_list.size(), combined_mult*Ensemble_list.size());
	    Eigen::MatrixXd Corr_Matrix(combined_mult*Ensemble_list.size(), combined_mult*Ensemble_list.size());
	    Cov_Matrix.setZero();
	    Corr_Matrix.setZero();
	    
	    for(int iens=0; iens<(signed)Ensemble_list.size(); iens++) {
	      Cov_Matrix(iens,iens) = pow(data_OS.err(iens),2);
	      Cov_Matrix(iens+Ensemble_list.size(), iens+Ensemble_list.size()) = pow(data_tm.err(iens),2);
	      Corr_Matrix(iens,iens) = 1;
	      Corr_Matrix(iens+Ensemble_list.size(), iens+Ensemble_list.size()) = 1;

	      
	      Cov_Matrix(iens, iens + Ensemble_list.size()) = data_OS.distr_list[iens]%data_tm.distr_list[iens];
	      Cov_Matrix(iens+Ensemble_list.size(), iens) = Cov_Matrix(iens, iens+Ensemble_list.size());
	      
	      Corr_Matrix(iens, iens + Ensemble_list.size()) = data_OS.distr_list[iens]%data_tm.distr_list[iens]/(data_tm.err(iens)*data_OS.err(iens));
	      Corr_Matrix(iens+Ensemble_list.size(),iens) = Corr_Matrix(iens, iens+Ensemble_list.size());
	      
	      bf_R.Add_covariance_matrix(Cov_Matrix);
	      bf_R_ch2.Add_covariance_matrix(Cov_Matrix);
	    }
	  }
	  else {
	    bf_R.Disable_correlated_fit();
	    bf_R_ch2.Disable_correlated_fit();
	  }
	  //#################################################################
	  
	  vector<vector<ipar_R>> data_boot(Njacks);
	  vector<vector<ipar_R>> data_boot_ch2(1);
	  for(auto &dt: data_boot) dt.resize(combined_mult*Ensemble_list.size());
	  for(auto &dt: data_boot_ch2) dt.resize(combined_mult*Ensemble_list.size());
	  boot_fit_data<fpar_R>  Bt_fit;
	  boot_fit_data<fpar_R>  Bt_fit_ch2;

	  for(int ijack=0;ijack<Njacks;ijack++) {
	    for(int iens=0; iens<(signed)Ensemble_list.size();iens++) {
	      data_boot[ijack][iens].meas= data_OS.distr_list[iens].distr[ijack];
	      data_boot[ijack][iens].err= data_OS.err(iens);
	      data_boot[ijack][iens].a= a_distr_list.distr_list[iens].distr[ijack];
	      data_boot[ijack][iens].Is_tm=false;
	      if(flav != "disco") {
	      data_boot[ijack][iens+Ensemble_list.size()].meas= data_tm.distr_list[iens].distr[ijack];
	      data_boot[ijack][iens+Ensemble_list.size()].err= data_tm.err(iens);
	      data_boot[ijack][iens+Ensemble_list.size()].a= a_distr_list.distr_list[iens].distr[ijack];
	      data_boot[ijack][iens+Ensemble_list.size()].Is_tm=true;
	      }

	      if(ijack==0) { //mean values
		data_boot_ch2[ijack][iens].meas= data_OS.ave(iens);
		data_boot_ch2[ijack][iens].err= data_OS.err(iens);
		data_boot_ch2[ijack][iens].a= a_distr_list.ave(iens);
		data_boot_ch2[ijack][iens].Is_tm=false;
		if(flav != "disco") {
		data_boot_ch2[ijack][iens+Ensemble_list.size()].meas= data_tm.ave(iens);
		data_boot_ch2[ijack][iens+Ensemble_list.size()].err= data_tm.err(iens);
		data_boot_ch2[ijack][iens+Ensemble_list.size()].a= a_distr_list.ave(iens);
		data_boot_ch2[ijack][iens+Ensemble_list.size()].Is_tm=true;
		}
	      }
	    }
	  }



	  bf_R.Append_to_input_par(data_boot);
	  bf_R_ch2.Append_to_input_par(data_boot_ch2);
	  //fit
	  cout<<"####FIT TO:######"<<endl;
	  cout<<"flavor: "<<flav<<endl;
	  cout<<"(alpha,Emax): ("<<beta<<","<<Emax<<")"<<endl;
	  cout<<"sigma, E: "<<sigma<<", "<<Erg<<endl;
	  Bt_fit= bf_R.Perform_bootstrap_fit();
	  Bt_fit_ch2= bf_R_ch2.Perform_bootstrap_fit();

	

	  //Print result and store the ch2
	  distr_t D(UseJack), D2_tm(UseJack), D2_OS(UseJack);
	  distr_t_list f_func_tm(UseJack), f_func_OS(UseJack);

	  for(int ijack=0;ijack<Njacks;ijack++) {
	    D.distr.push_back( Bt_fit.par[ijack].D);
	    D2_tm.distr.push_back( Bt_fit.par[ijack].D2_tm);
	    D2_OS.distr.push_back( Bt_fit.par[ijack].D2_OS);
	  }

	  R_ratio_fixed_sigma.distr_list.push_back(D);
	  ch2_list.push_back( Bt_fit_ch2.get_ch2_ave());
	  for(auto &a: lat_to_print) { f_func_tm.distr_list.push_back( (D+ pow(a*Lambda_QCD/0.197327,2)*D2_tm)); f_func_OS.distr_list.push_back( (D+ pow(a*Lambda_QCD/0.197327,2)*D2_OS)); }
	  //print
	  Print_To_File({}, {lat_to_print, f_func_tm.ave(), f_func_tm.err(), f_func_OS.ave(), f_func_OS.err()}, "../data/R_ratio/"+Tag_reco_type+"/continuum/"+flav+"/"+to_string_with_precision(sigma,3)+"/E_"+to_string_with_precision(Erg,3)+".interpol", "", "# a[fm]   tm    OS");
	  Print_To_File({Ensemble_list}, {data_tm.ave(), data_tm.err(), data_OS.ave(), data_OS.err()}, "../data/R_ratio/"+Tag_reco_type+"/continuum/"+flav+"/"+to_string_with_precision(sigma,3)+"/E_"+to_string_with_precision(Erg,3)+".dat", "", "#Ens  tm  OS");
	  
	}

	R_ratio_flav[count_flav][is] = R_ratio_fixed_sigma;

	//print R_ratio at this sigma and ch2 info
	Print_To_File({}, {Ergs_GeV_list, R_ratio_fixed_sigma.ave(), R_ratio_fixed_sigma.err(), ch2_list}, "../data/R_ratio/"+Tag_reco_type+"/continuum/"+flav+"/"+to_string_with_precision(sigma,3)+"/continuum.dat", "", "#Erg[GeV]  R_ratio   Ch2");
      }
      count_flav++;
    }
    //print total contribution for each sigma
    //accumulate contributions
    vector<distr_t_list> R_ratio_total_per_sigma(sigmas.size());
    for(int iflav=0;iflav<(signed)flavors.size();iflav++) {
      for(int is=0;is<(signed)sigmas.size();is++) {
	if(iflav==0) R_ratio_total_per_sigma[is] = R_ratio_flav[iflav][is];
	else R_ratio_total_per_sigma[is] = R_ratio_total_per_sigma[is]+R_ratio_flav[iflav][is];
      }
    }
    //print
    for(int is=0;is<(signed)sigmas.size();is++)   Print_To_File({}, {Ergs_GeV_list, R_ratio_total_per_sigma[is].ave(), R_ratio_total_per_sigma[is].err()}, "../data/R_ratio/"+Tag_reco_type+"/continuum/total/"+to_string_with_precision(sigmas[is],3)+"/continuum.dat", "", "#Erg[GeV] R_ratio"); 
    
  }




  omp_set_num_threads(NUMM_THREADS);
  cout<<"Bye!"<<endl;
  return;
}
