#include "../include/PionMassAnalysis_ov_X.h"


using namespace std;
namespace plt = matplotlibcpp;

constexpr double kappa=2.837297;
const double MPiPhys=0.134977;
const double alpha = 1.0/137.04;
const double e2 = alpha*4.0*M_PI;
const double fpi_phys = 0.130;
const double r0 = pow(0.672/0.197,2);
const double f0 = 0.121;
const bool Enable_f0 =0;
const double csi_phys = pow(MPiPhys,2)/pow(4.0*M_PI*(Enable_f0?f0:fpi_phys),2);
const int Nbranches = 8;
const int nboots= 100;
const bool Use_JB_distribution= false;
const bool UseJack=1;
const int Njacks=15;
const int Nboots=200;
const bool verbose=1;
const bool CDH_correct_FVE=1;
const bool GL_correct_FVE =CDH_correct_FVE?0:0;
const bool Use_fp=1;
const bool USE_RATIONAL_FUNCTION=0;
const double X_phys_val= 1.0e+3*(Use_fp?pow(fpi_phys,2)/(2*MPiPhys):pow(MPiPhys,2)/(2*MPiPhys));
const double l1ph= -0.4;
const double l2ph= 4.3;
const double l3ph= 3.20;
const double l4ph= 4.4;
const double s0= 2.0-M_PI/2.0;
const double s1 = M_PI/4.0 - 0.5;
const double s2 = 0.5 - M_PI/8.0;
const double s3 = 3.0*M_PI/16.0 - 0.5;



class Y_M {

public:
  Y_M() {}
  double L, Mpi, Mpi_err, fp, csi, Dm_ov_X_err, Dm_ov_X_val, f0;
  int ibeta;
};



class X_M {

public:
  X_M(const Vfloat &par) : ainv(3) {
    if(par.size() != 13) {
      cout<<"In class X_M, invalid call to constructor"<<endl;
      crash("par_size: +"+to_string(par.size()));
    }
    this->chir=par[0];
    this->log=par[1];
    this->A_1=par[2];
    this->D=par[3];
    this->Dm=par[4];
    this->F_a=par[5];
    this->log_a=par[6];
    this->F_m = par[7];
    this->A_2=par[8];
    this->Af0 = par[9];
    for(int ibeta=0;ibeta<3;ibeta++) {
      this->ainv[ibeta] = par[10+ibeta];
    }
 
  }
  X_M() : ainv(3) {}
  vector<double> ainv;
  double D, Dm, chir, log, log_a, A_1, F_a, F_m, A_2, Af0;
};





void Pion_mass_analysis_ov_X(string CURRENT_TYPE, bool IncludeDisconnected) {

  omp_set_num_threads(1);
  
  data_t m_data, dm_exch_data, dm_hand_data, m_data_hand_run, m_twisted_data;
   if(CURRENT_TYPE=="CONSERVED")  {
    m_data.Read("../datasets/data", "mes_contr_00", "P5P5");
    dm_exch_data.Read("../datasets/data", "mes_contr_LL", "P5P5");
  }
  else crash("CURRENT_TYPE: "+CURRENT_TYPE+ " not yet implemented. Exiting...");
  
  
  //compute masses and effective masses

  
  
  int Nens = m_data.size; 
  cout<<"N_ens: "<<Nens<<endl;



  distr_t_list Mpi_distr_list(UseJack),  Dm_ov_X_distr_list(UseJack), X_pi_distr_list(UseJack), Dm2_univ_sub(UseJack);
  distr_t_list Mpi_distr_dim_list(UseJack);
  distr_t_list fp_fit_distr_list(UseJack);
   distr_t_list fp_fit_naive_list(UseJack);
  distr_t_list fp_fit_GL_list(UseJack);
  distr_t_list fp_fit_CDH_list(UseJack);
  distr_t_list Mpi_fit_naive_list(UseJack);
  distr_t_list Mpi_fit_GL_list(UseJack);
  distr_t_list Mpi_fit_CDH_list(UseJack);
  vector<Eigen::MatrixXd> CovMatrixPion(Nens);

  for(int i=0; i < Nens; i++) {

    double ll, mm, tt;
    Read_pars_from_ensemble_tag(m_data.Tag[i], mm, ll, tt);
    CorrAnalysis Corr(UseJack, Njacks,Nboots);
    if(m_data.Tag[i].substr(0,1) == "A") {
      if(int(ll) == 20) Corr.Tmin=13;
      else Corr.Tmin = 14;
    }
    else if(m_data.Tag[i].substr(0,1) =="B") Corr.Tmin = 15;
    else Corr.Tmin= 21;
    if(m_data.Tag[i].substr(0,1)=="D") Corr.Tmax = 36; //m_data.nrows[i]/2 -5;
    else {
      if(m_data.Tag[i].substr(0,1) == "A") {
	if((int)ll==20) Corr.Tmax=17;
	else if((int)ll == 24) Corr.Tmax= 20;
	else if((int)ll == 32) Corr.Tmax= 24;
	else if((int)ll == 40) Corr.Tmax=30;
	else if((int)ll == 48) Corr.Tmax= 36;
	else crash("volume does not exist V: "+to_string_with_precision(ll, 4));
      }
      else if(m_data.Tag[i].substr(0,1) == "B") {
	Corr.Tmax = m_data.nrows[i]/2 -5;

      }
      else crash("lattice spacing not found when fixing time intervals");
    }
    Corr.Nt = m_data.nrows[i];
    boost::filesystem::create_directory("../data");
    boost::filesystem::create_directory("../data/Mpi_ov_X");
    string p= "conserved";
    boost::filesystem::create_directory("../data/Mpi_ov_X/"+p);

    LatticeInfo LL("CONSERVED");
    LL.LatInfo(m_data.Tag[i].substr(0,1));
    double lat_spacing_inv = LL.ainv;
    GaussianMersenne G(444363);
       
    auto F_15 = [](double x) { return pow(x,1.0/5);};
    auto LOG = [](double x) { return log(x);};
   
   
    distr_t_list Pi_iso_distr = Corr.corr_t(m_data.col(0)[i], "");

    distr_t_list Mpi_eff_distr = Corr.effective_mass_t(Pi_iso_distr,  "../data/Mpi_ov_X/"+p+"/mass."+m_data.Tag[i]);

    
    distr_t_list Exch_distr =  Corr.corr_t(dm_exch_data.col(0)[i], "");

    distr_t_list fp_distr= Corr.decay_constant_t(pow(2.0*mm,2)*Pi_iso_distr, "../data/Mpi_ov_X/"+p+"/fp."+m_data.Tag[i]);

    distr_t fp_fit= Corr.Fit_distr(fp_distr);

    distr_t Mpi_fit = Corr.Fit_distr(Mpi_eff_distr);

    
    //############# CORRECT Fpi and Mpi   Usi G&L and CDH formulae ################
    distr_t fp_fit_GL, fp_fit_CDH, Mpi_fit_GL, Mpi_fit_CDH;
    distr_t csi_L = Mpi_fit*Mpi_fit/(pow(4.0*M_PI,2)*fp_fit*fp_fit);
    distr_t g1 = distr_t::f_of_distr(g1_l, Mpi_fit*ll);
    distr_t g2 = distr_t::f_of_distr(g2_l, Mpi_fit*ll);
    distr_t log_l = log(csi_phys) - distr_t::f_of_distr(LOG, csi_L);

    fp_fit_GL = fp_fit/(1.0-1.2*2.0*csi_L*g1);
    Mpi_fit_GL = Mpi_fit/(1.0 + 1.2*0.5*csi_L*g1);
    fp_fit_CDH = fp_fit/(1.0 -2.0*csi_L*g1 +2.0*csi_L*csi_L*( (Cf1(l1ph,l2ph,l3ph,l4ph) + Sf1(s0,s1,s2,s3) + Cf1_log()*log_l)*g1 + (Cf2(l1ph,l2ph,l3ph,l4ph) + Sf2(s0,s1,s2,s3) + Cf2_log()*log_l)*g2));
    Mpi_fit_CDH = Mpi_fit/(1.0 + 0.5*csi_L*g1 - csi_L*csi_L*( (Cm1(l1ph,l2ph,l3ph,l4ph) + Sm1(s0,s1,s2,s3) + Cm1_log()*log_l)*g1 + (Cm2(l1ph,l2ph,l3ph,l4ph) + Sm2(s0,s1,s2,s3) + Cm2_log()*log_l)*g2));


    //#################### STORE NAIVE, GL and CDH corrected data######################
    Mpi_fit_naive_list.distr_list.push_back(Mpi_fit);
    Mpi_fit_GL_list.distr_list.push_back(Mpi_fit_GL);
    Mpi_fit_CDH_list.distr_list.push_back(Mpi_fit_CDH);
    fp_fit_naive_list.distr_list.push_back(fp_fit);
    fp_fit_GL_list.distr_list.push_back(fp_fit_GL);
    fp_fit_CDH_list.distr_list.push_back(fp_fit_CDH);

    //############ CHOOSE WHICH ONE TO USE IN THE FIT ############################

    if(GL_correct_FVE || CDH_correct_FVE) {
      if(CDH_correct_FVE) { fp_fit = fp_fit_CDH; Mpi_fit = Mpi_fit_CDH;}
      else { fp_fit = fp_fit_GL; Mpi_fit = Mpi_fit_GL;} 
    }


    
    X_pi_distr_list.distr_list.push_back( distr_t::f_of_distr(F_15, fp_fit*Mpi_fit*Mpi_fit*Mpi_fit*Mpi_fit));

    distr_t_list Dm_exch_eff_distr = Corr.effective_slope_t(Exch_distr, Pi_iso_distr, "../data/Mpi_ov_X/"+p+"/dm."+m_data.Tag[i]);
   
    
    Mpi_distr_list.distr_list.push_back(Corr.Fit_distr( Mpi_eff_distr));
    Mpi_distr_dim_list.distr_list.push_back( lat_spacing_inv*Corr.Fit_distr(Mpi_eff_distr));
    distr_t_list Dm_tot_eff_distr=Dm_exch_eff_distr;

    distr_t_list Dm_Ov_X_eff_distr= (2.0*Mpi_eff_distr*Dm_tot_eff_distr + (kappa/(2.0*M_PI))*(1.0/ll)*(Mpi_eff_distr + 2.0/ll))/(Use_fp?(fp_distr*fp_distr):(Mpi_eff_distr*Mpi_eff_distr));

    Print_To_File({}, {Dm_Ov_X_eff_distr.ave(), Dm_Ov_X_eff_distr.err()},  "../data/Mpi_ov_X/"+p+"/Dm_ov_X_eff."+m_data.Tag[i]+".dat.t", "", "");


    distr_t Dm_Ov_X_fit = Corr.Fit_distr(Dm_Ov_X_eff_distr);


       
    //subtract universal FVE
    distr_t mL = Mpi_fit*ll;
    Dm2_univ_sub.distr_list.push_back( (2.0*Mpi_fit*Corr.Fit_distr(Dm_tot_eff_distr) + (kappa/(2.0*M_PI))*(1.0/ll)*( Mpi_fit + 2.0/ll)));
    Dm_ov_X_distr_list.distr_list.push_back ( ((2.0/Mpi_fit)*Corr.Fit_distr(Dm_tot_eff_distr) + (kappa/(2.0*M_PI*mL))*( 1.0 + 2.0/mL))*(Use_fp?(Mpi_fit*Mpi_fit/(fp_fit*fp_fit)):(Mpi_fit/Mpi_fit)));

    fp_fit_distr_list.distr_list.push_back(fp_fit);

    
    if(!Use_JB_distribution) { //resample from Gaussian distribution
      Compute_covariance_matrix(UseJack, CovMatrixPion[i],3, Mpi_distr_list.distr_list[i].distr, Dm_ov_X_distr_list.distr_list[i].distr, fp_fit_distr_list.distr_list[i].distr);
      
    }

    cout<<"################### ANALYZED ENSEMBLE : "<<m_data.Tag[i]<<" ############################ "<<endl;
    cout<<"PRINTING INFO: "<<endl;
    cout<<"Mpi: "<<Mpi_fit.ave()<<" +- "<<Mpi_fit.err()<<endl;
    cout<<"Mpi (dim) : "<<(lat_spacing_inv*Mpi_fit).ave()<<" +- "<<(lat_spacing_inv*Mpi_fit).err()<<endl;
    cout<<"Mpi*L :"<<Mpi_fit.ave()*ll<<endl;
    cout<<"fp: "<<fp_fit.ave()<<" +- "<<fp_fit.err()<<endl;
    cout<<"Dm: "<<Corr.Fit_distr(Dm_tot_eff_distr).ave()<<" +- "<<Corr.Fit_distr(Dm_tot_eff_distr).err()<<endl;
    cout<<"Dm2: "<<(2.0*Mpi_fit*Corr.Fit_distr(Dm_tot_eff_distr)).ave()<<" +- "<<(2.0*Mpi_fit*Corr.Fit_distr(Dm_tot_eff_distr)).err()<<endl;
    cout<<"Dm2 (univ sub): "<<Dm2_univ_sub.ave(i)<<" +- "<<Dm2_univ_sub.err(i)<<endl;
    cout<<"Dm_ov_X: "<<Dm_ov_X_distr_list.ave(i)<<" +- "<<Dm_ov_X_distr_list.err(i)<<endl;
    cout<<"Dm_ov_X (fit eff): "<<Dm_Ov_X_fit.ave()<<" +- "<<Dm_Ov_X_fit.err()<<endl;
    cout<<"ainv (GeV-1): "<<lat_spacing_inv<<endl;
    cout<<"#########################################################################################################"<<endl;
  }


  
  //everything is set up. Start fitting!
  LatticeInfo L_info(CURRENT_TYPE);
  

  bootstrap_fit<X_M,Y_M> bf(nboots);
  bf.set_warmup_lev(3);

  bf.Set_number_of_measurements(Nens);
  bf.Set_verbosity(verbose);

  //Add parameters
  bf.Add_par("chir", (1.6e-4)/pow(f0,2), (4e-6)/pow(f0,2));
  bf.Add_par("log", 3.7 , 0.7);
  bf.Add_par("A_1", -5.7, 0.02);
  bf.Add_par("D", 2e-3, 1e-5);
  bf.Add_par("Dm", 0.5, 1e-3);
  bf.Add_par("F_a", 1, 0.1);
  bf.Add_par("lg_a", 0.1, 1e-3);
  bf.Add_par("F_m", 1.0, 0.1);
  bf.Add_par("A_2", 1.0, 0.01);
  bf.Add_par("Af0", 1.0, 0.1);
  bf.Add_prior_pars({"ainv0", "ainv1", "ainv2"});
  
  //Fix some parameters to make test
  //bf.Fix_par("F_m",1.0);
  //bf.Set_limits("F_m", 0.0, 5);
  bf.Fix_par("F_a",0.0);
  bf.Fix_par("Dm",0.0);
  //bf.Fix_par("D",0.0);
  bf.Fix_par("lg_a", 0.0);
  bf.Fix_par("A_2", 0.0);
  bf.Fix_par("Af0", 0.0);
  if(Use_fp || !Enable_f0) {
    bf.Fix_par("log",3.0);
  }
  //bf.Fix_par("A_1", 0.0);
  

  //Add List of parameters to be released after first minimization
  bf.Fix_n_release({"ainv0", "ainv1", "ainv2"});


  
  
  Vfloat L, m_l, T;
  Read_pars_from_ensemble_tag(m_data.Tag, m_l, L, T);

  


  //define fitting function

  //############################################################################################################

  bf.ansatz =  [=](const X_M &p, const Y_M &ip) -> double {
		 

		 
    double ainv = p.ainv[ip.ibeta];

    double csi_val = ip.csi*(Enable_f0?pow(ainv,2):1.0);

    double scale = Use_fp?pow(4.0*M_PI,2)*csi_val:1.0;

    double scale_fve = Use_fp?pow(ip.Mpi/ip.fp,2):1.0;

    

    double a=1.0/ainv;

    double L = ip.L/ainv;

    double Mp= ip.Mpi*ainv;

    double Mp2 = Mp*Mp;

    double ML = ip.Mpi*ip.L;

    
     
    double SD_FVE = (L>0.0)?(e2/3.0)*scale_fve*r0*(Mp2/pow(ML,3))*(p.F_m + p.F_a*a*a):0.0;

     
    double den = USE_RATIONAL_FUNCTION?pow(1.0  -2.0*csi_val*log(csi_val) + p.Af0*csi_val,2):1.0;
    
    double fitted_value = scale*(e2/pow(4.0*M_PI,2))*( p.chir/csi_val -p.log*log(csi_val) + p.A_1 + p.A_2*csi_val + (p.D/(csi_val))*a*a + p.Dm*a*a)/den;


    fitted_value += SD_FVE;
     

    
    

    return X_phys_val*fitted_value;
  };
  

  bf.measurement = [=](const X_M& p,const Y_M& ip) -> double {

            
    double m1 = (e2/2)*ip.Dm_ov_X_val;
      
    return X_phys_val*m1;
  };
  

  bf.error =  [=](const X_M& p,const  Y_M &ip) -> double {

		double m1 = (e2/2)*ip.Dm_ov_X_err;
		//	double m1 = (e2/2)*max(1e-2*ip.Dm_ov_X_val, ip.Dm_ov_X_err);

    return X_phys_val*m1;
  };

 

  
  //############################################################################################################

  //init random Number Generators

  GaussianMersenne GM(54353);
  RandomMersenne RM(23423, UseJack?(Njacks-1):(Nboots-1));

  Eigen::MatrixXd CovMatrixInput(9,9); //covariance matrix of input parameters
  Eigen::VectorXd Ave_input_parameters(9);

  vector<boot_fit_data<X_M>> Bt_fit;

  //we store here all resampled data. 
  vector<vector<vector<Y_M>>> all_ens(Nbranches);
  for(auto &all_ens_br : all_ens) {
    all_ens_br.resize(nboots);
    for(auto & all_ens_ib: all_ens_br) {
      all_ens_ib.resize(Nens);
      for(int iens=0; iens<Nens;iens++) {

	all_ens_ib[iens].L=L[iens];
	if(m_data.Tag[iens].substr(0,1) == "A") all_ens_ib[iens].ibeta=0;
	else if(m_data.Tag[iens].substr(0,1) == "B") all_ens_ib[iens].ibeta=1;
	else all_ens_ib[iens].ibeta=2;
	all_ens_ib[iens].Mpi_err = Mpi_distr_list.err(iens);
	all_ens_ib[iens].Dm_ov_X_err = Dm_ov_X_distr_list.err(iens);
	

      }
    }
  }

 
  for(int ibranch=0; ibranch < Nbranches; ibranch++) {


    
    ReadBranch(ibranch, CovMatrixInput, Ave_input_parameters );


    //clear all measurements
    bf.Clear_priors();
    bf.Clear_input_pars();
    
    

    //generate bootstrap data
    for(int iboot=0; iboot<nboots; iboot++) {
      bf.ib= &iboot;
      Vfloat lat_input = Covariate(CovMatrixInput, Ave_input_parameters, GM);
      double f0_resampled = lat_input[2];
      bf.Append_to_prior("ainv0", lat_input[3], sqrt(CovMatrixInput(3,3)));
      bf.Append_to_prior("ainv1", lat_input[4], sqrt(CovMatrixInput(4,4)));
      bf.Append_to_prior("ainv2", lat_input[5], sqrt(CovMatrixInput(5,5)));


      int k= RM();
      for(int imeas=0; imeas <Nens; imeas++) {
	if(!Use_JB_distribution) {
	  Eigen::VectorXd vec(CovMatrixPion[imeas].rows());
	 
	  vec<<Mpi_distr_list.ave(imeas),Dm_ov_X_distr_list.ave(imeas),fp_fit_distr_list.ave(imeas);
	  Vfloat res_meas(3,0.0);
	 
	  res_meas= Covariate(CovMatrixPion[imeas],vec, GM);
	  all_ens[ibranch][iboot][imeas].Mpi= res_meas[0];
	  all_ens[ibranch][iboot][imeas].Dm_ov_X_val = res_meas[1];
	  all_ens[ibranch][iboot][imeas].fp = res_meas[2];
	  all_ens[ibranch][iboot][imeas].f0 = f0_resampled;
	  all_ens[ibranch][iboot][imeas].csi = pow(res_meas[0],2)/pow(4.0*M_PI*(Enable_f0?f0_resampled:res_meas[2]),2);
	 	 
	}
      
	else {
	  double M,FP;
	  M= Mpi_distr_list.distr_list[imeas].distr[k];
	  FP=fp_fit_distr_list.distr_list[imeas].distr[k];
	  all_ens[ibranch][iboot][imeas].Mpi = M;
	  all_ens[ibranch][iboot][imeas].Dm_ov_X_val = Dm_ov_X_distr_list.distr_list[imeas].distr[k];
	  all_ens[ibranch][iboot][imeas].fp = FP;
	  all_ens[ibranch][iboot][imeas].f0 = f0_resampled;
	  all_ens[ibranch][iboot][imeas].csi = pow(M,2)/pow(4.0*M_PI*(Enable_f0?f0_resampled:FP),2);
	}
      }

    }
    
    bf.Append_to_input_par(all_ens[ibranch]);
    Bt_fit.push_back(bf.Perform_bootstrap_fit());
    
  }
  
 
  //print the data
  
  //define lambda functions


  auto SD_FVE = [=](double ML, double M, double FP, double ainv,  double F_a, double F_m) {

		  double resc= Use_fp?pow(M/FP,2):1.0;
		  double Mp2= pow(M*ainv,2);
		  double SDE =(e2/3.0)*(Mp2/pow(ML,3))*r0*(F_m + F_a/pow(ainv,2));
		  return X_phys_val*SDE*resc;
  };


  auto ANSATZ =  [=](const X_M &p, const Y_M &ip) -> double {
		 

		 
    double ainv = p.ainv[ip.ibeta];

    double csi_val = ip.csi;

    double scale = Use_fp?pow(4.0*M_PI,2)*csi_val:1.0;

    double scale_fve = Use_fp?pow(4.0*M_PI,2)*ip.csi:1.0;

    

    double a=1.0/ainv;

    double L = ip.L/ainv;

    double Mp= ip.Mpi*ainv;

    double Mp2 = Mp*Mp;

    double ML = ip.Mpi*ip.L;

    
     
    double SD_FVE = (L>0.0)?(e2/3.0)*scale_fve*r0*(Mp2/pow(ML,3))*(p.F_m + p.F_a*a*a):0.0;

     
    double den = USE_RATIONAL_FUNCTION?pow(1.0  -2.0*csi_val*log(csi_val) + p.Af0*csi_val,2):1.0;
    
    double fitted_value = scale*(e2/pow(4.0*M_PI,2))*( p.chir/csi_val -p.log*log(csi_val) + p.A_1 + p.A_2*csi_val + (p.D/(csi_val))*a*a + p.Dm*a*a)/den;


    fitted_value += SD_FVE;
     

    
    

    return X_phys_val*fitted_value;
  };


  auto cont_ansatz = [=](X_M p, double csi) -> double {


		       double scale = Use_fp?pow(4.0*M_PI,2)*csi:1.0;
		       double val = scale*X_phys_val*(e2/pow(4.0*M_PI,2))*( p.chir/csi - p.log*log(csi) + p.A_1 + p.A_2*csi);

		    
		       return val;
  };



 
  
  
  Vfloat Csi(1000);
  Vfloat vols(600);
  for(unsigned int m=0; m<Csi.size();m++) Csi[m] = 0.0001*m+ 0.0002;
  for(unsigned int vol=0;vol<vols.size();vol++) vols[vol] = 2.0 + vol*0.03;
  VVVfloat SD_subtracted_data, Csi_measured, raw_data;
  VVfloat C,LOG,A1,A2,K,Fa,D, Dm, Ainv0, Ainv1, Ainv2, Af0;
  VVVVfloat Fitted_func;
  VVfloat Physical_point;
  VVVfloat A40_slice;
  cascade_resize(SD_subtracted_data, Vint{bf.Get_number_of_measurements(), Nbranches, nboots});
  cascade_resize(Csi_measured, Vint{bf.Get_number_of_measurements(), Nbranches, nboots});
  cascade_resize(Fitted_func, Vint{4, (int)Csi.size(), Nbranches, nboots});
  cascade_resize(Physical_point, Vint{Nbranches, nboots});
  cascade_resize(C, Vint{Nbranches, nboots});
  cascade_resize(LOG, Vint{Nbranches, nboots});
  cascade_resize(A1, Vint{Nbranches, nboots});
  cascade_resize(A2, Vint{Nbranches, nboots});
  cascade_resize(K, Vint{Nbranches, nboots});
  cascade_resize(Fa, Vint{Nbranches, nboots});
  cascade_resize(D, Vint{Nbranches, nboots});
  cascade_resize(Dm, Vint{Nbranches, nboots});
  cascade_resize(Ainv0, Vint{Nbranches, nboots});
  cascade_resize(Ainv1, Vint{Nbranches, nboots});
  cascade_resize(Ainv2, Vint{Nbranches, nboots});
  cascade_resize(Af0, Vint{Nbranches, nboots});
  cascade_resize(A40_slice, Vint{ (int)vols.size(), Nbranches,nboots});
  cascade_resize(raw_data, Vint{bf.Get_number_of_measurements(), Nbranches, nboots});


 
  
  for(int ibranch=0; ibranch<Nbranches;ibranch++) {
    for(int iboot=0; iboot<nboots;iboot++) {
      X_M P = Bt_fit[ibranch].par[iboot];

      //store fitted params
      C[ibranch][iboot]= P.chir;
      LOG[ibranch][iboot] =P.log;
      A1[ibranch][iboot] =P.A_1;
      A2[ibranch][iboot]= P.A_2;
      K[ibranch][iboot]= P.F_m;
      Fa[ibranch][iboot] =P.F_a;
      D[ibranch][iboot]= P.D;
      Dm[ibranch][iboot] =P.Dm;
      Af0[ibranch][iboot] = P.Af0;
      Ainv0[ibranch][iboot] = P.ainv[0];
      Ainv1[ibranch][iboot] = P.ainv[1];
      Ainv2[ibranch][iboot] = P.ainv[2];
      for(unsigned int ivol=0; ivol<vols.size();ivol++) {
	Y_M new_p;
	double f0_boot = all_ens[ibranch][iboot][0].f0/P.ainv[0];
        double csi_A40 = 0.0275;
	new_p.ibeta=0;
	new_p.Mpi= 0.142;  //pion mass of A40.XX in lattice units
	new_p.L = vols[ivol]/new_p.Mpi;
	new_p.csi = Enable_f0?pow(new_p.Mpi,2)/pow(4.0*M_PI*f0_boot,2):csi_A40;
	new_p.fp = new_p.Mpi/(4.0*M_PI*sqrt(new_p.csi));
	A40_slice[ivol][ibranch][iboot] = ANSATZ(P,new_p);
      }
    
      for(int imeas=0; imeas< bf.Get_number_of_measurements(); imeas ++) {
	Y_M E = all_ens[ibranch][iboot][imeas];
	double fval = Enable_f0?(E.f0/P.ainv[E.ibeta]):E.fp;
	SD_subtracted_data[imeas][ibranch][iboot] = bf.measurement(P, E) - SD_FVE(E.L*E.Mpi, E.Mpi, E.fp, P.ainv[E.ibeta], P.F_a, P.F_m);
	raw_data[imeas][ibranch][iboot] = bf.measurement(P,E);
	Csi_measured[imeas][ibranch][iboot] = pow(E.Mpi,2)/pow(4.0*M_PI*fval,2);
      }
     
      for(int ibeta=0; ibeta < 4;ibeta++) {
	for(unsigned int m=0; m<Csi.size(); m++){
	  if(ibeta<3) {
	  
	    Y_M point;
	    point.L = -1.0;
	    point.ibeta=ibeta;
	    point.Mpi = 1.0;
	    point.csi= Csi[m];
	    Fitted_func[ibeta][m][ibranch][iboot] = ANSATZ(P,point);
	    
	  }
	  else Fitted_func[ibeta][m][ibranch][iboot] = cont_ansatz(P, Csi[m]) ;
	}
      }
  
      //Add physical point
      Y_M E2 = all_ens[ibranch][iboot][0];
      Physical_point[ibranch][iboot] = cont_ansatz(P, csi_phys );
   
    }
  }




 

  Vfloat SD_subtracted_val, SD_subtracted_err, CSI, CSI_err,  raw_data_val, raw_data_err;
  for(int imeas=0;imeas<bf.Get_number_of_measurements();imeas++) {
    SD_subtracted_val.push_back( Boot_ave(SD_subtracted_data[imeas]));
    SD_subtracted_err.push_back( Boot_err(SD_subtracted_data[imeas]));
    raw_data_val.push_back(Boot_ave(raw_data[imeas]));
    raw_data_err.push_back(Boot_err(raw_data[imeas]));
    CSI.push_back( Boot_ave(Csi_measured[imeas]));
    CSI_err.push_back( Boot_err(Csi_measured[imeas]));
  }
  
  VVfloat Fitted_func_val, Fitted_func_err;
  Vfloat A40_slice_val, A40_slice_err, A40_raw_val, A40_raw_err;
 
  cascade_resize(Fitted_func_val, {4,(int)Csi.size()});
  cascade_resize(Fitted_func_err, {4,(int)Csi.size()});


  for(int ibeta=0;ibeta<4;ibeta++) {
    for(unsigned int m=0; m<Csi.size();m++) {
	
	Fitted_func_val[ibeta][m] = Boot_ave(Fitted_func[ibeta][m]);
	Fitted_func_err[ibeta][m] = Boot_err(Fitted_func[ibeta][m]);

    }
  }

  for(auto &A40_boot : A40_slice) {
    A40_slice_val.push_back( Boot_ave(A40_boot));
    A40_slice_err.push_back( Boot_err(A40_boot));
  }
  //plot the result

  Vfloat vol_A40, meas_A40, err_meas_A40;
  for(int iens=0; iens<Nens;iens++)
    if(m_data.Tag[iens].substr(0,3)=="A40") {
      vol_A40.push_back( L[iens]*Mpi_distr_list.ave(iens));
      meas_A40.push_back(raw_data_val[iens]);
      err_meas_A40.push_back(raw_data_err[iens]);
      A40_raw_val.push_back(raw_data_val[iens]);
      A40_raw_err.push_back(raw_data_err[iens]);
    }


  string print_path = (CURRENT_TYPE=="CONSERVED")?"conserved":"local";
  string exch_or_tot = IncludeDisconnected?"tot":"exch";

    
  


   
  //save data in files
  boost::filesystem::create_directory("../data");
  boost::filesystem::create_directory("../data/Mpi_ov_X/"+print_path);
  Print_To_File({}, {Csi, Fitted_func_val[0], Fitted_func_err[0],  Fitted_func_val[1], Fitted_func_err[1], Fitted_func_val[2], Fitted_func_err[2], Fitted_func_val[3], Fitted_func_err[3]} , "../data/Mpi_ov_X/"+print_path+"/Fitted_func_"+exch_or_tot+".dat", "OUT", "#Csi     beta=1.90      beta=1.95      beta=2.10           cont");
 
  
  Print_To_File(m_data.Tag, { CSI, CSI_err, Mpi_distr_dim_list.ave(), L, (Mpi_distr_list*L).ave(),  SD_subtracted_val, SD_subtracted_err}, "../data/Mpi_ov_X/"+print_path+"/sd_subtracted_data_"+exch_or_tot+".dat", "", "#Ens   #Csi #Csi_err     #Mpi[GeV] #L/a  #Mpi*L   #Mpi^2_+ - Mpi^2_0      #err");

  Print_To_File(m_data.Tag, { CSI, CSI_err, Mpi_distr_dim_list.ave(), L, (Mpi_distr_list*L).ave(), X_pi_distr_list.ave(), X_pi_distr_list.err(), Mpi_distr_list.ave(), Mpi_distr_list.err(), fp_fit_distr_list.ave(), fp_fit_distr_list.err()}, "../data/Mpi_ov_X/"+print_path+"/Xi_pi_"+exch_or_tot+".dat", "", "#Ens   #Csi #Csi_err    #Mpi L/a #Mpi*L       #Xi    #Mpi    #fp");
 
  Print_To_File(m_data.Tag, { CSI, CSI_err, Mpi_distr_dim_list.ave(), L, (Mpi_distr_list*L).ave(), Dm2_univ_sub.ave(), Dm2_univ_sub.err(), raw_data_val, raw_data_err, ((e2*X_phys_val/2.0)*Dm_ov_X_distr_list).ave(), ((e2*X_phys_val/2.0)*Dm_ov_X_distr_list).err()}, "../data/Mpi_ov_X/"+print_path+"/data_"+exch_or_tot+".dat", "", "#Ens      #Csi #Csi_err  #Mpi   #L/a   #Mpi*L #Dm2  Raw_data(Resampled)   Raw_data  ");
 
  Print_To_File({}, {vol_A40, meas_A40, err_meas_A40, A40_raw_val, A40_raw_err}, "../data/Mpi_ov_X/"+print_path+"/A40_"+exch_or_tot+".dat","", "#Mpi*L         univ_meas            univ_err       raw_meas      raw_err");
   //print A40 prediction
  Print_To_File({},{vols, A40_slice_val, A40_slice_err},"../data/Mpi_ov_X/"+print_path+"/A40_pred_"+exch_or_tot+".dat","", "#Mpi*L    val  err");


  Print_To_File(m_data.Tag, {L, (Mpi_distr_list*L).ave(),  Mpi_fit_naive_list.ave(), Mpi_fit_naive_list.err(), Mpi_fit_GL_list.ave(), Mpi_fit_GL_list.err(), Mpi_fit_CDH_list.ave(), Mpi_fit_CDH_list.err()}, "../data/Mpi_ov_X/"+print_path+"/FVE_Mpi.dat", "", "#Ens #L  #Mpi*L   #Mpi_naive #Mpi_GL #Mpi_CDH");

  Print_To_File(m_data.Tag, {L, (Mpi_distr_list*L).ave(),  fp_fit_naive_list.ave(), fp_fit_naive_list.err(), fp_fit_GL_list.ave(), fp_fit_GL_list.err(), fp_fit_CDH_list.ave(), fp_fit_CDH_list.err()}, "../data/Mpi_ov_X/"+print_path+"/FVE_fp.dat", "", "#Ens #L  #Mpi*L   #fp_naive #fp_GL #fp_CDH");


		
  string command = "echo "+to_string_with_precision(csi_phys,8)+"\t\t"+to_string_with_precision(Boot_ave(Physical_point),8)+"\t\t"+to_string_with_precision(Boot_err(Physical_point),8)+" > ../data/Mpi_ov_X/"+print_path.c_str()+"/Phys_val_"+exch_or_tot+".dat";
  system(command.c_str());

 

   

  cout<<"Physical Pion mass difference: "<< Boot_ave(Physical_point)<<"  +-  "<<Boot_err(Physical_point)<<" [ MeV ]"<<endl;
  cout<<"On single branches: "<<endl;
  for(int ibr=0; ibr<Nbranches;ibr++) cout<<"Branch: "<<ibr<<": "<<Boot_ave(Physical_point[ibr])<<" +-  "<<Boot_err(Physical_point[ibr])<<" [ MeV ]"<<endl;
  cout<<"Average chi2:"<<endl;
  for(int ibr=0;ibr<Nbranches;ibr++) cout<<"Branch: "<<ibr<<" chi2 = "<<accumulate(Bt_fit[ibr].chi2.begin(), Bt_fit[ibr].chi2.end(), 0.0)/Bt_fit[ibr].chi2.size()<<endl;




  //print averaged parameters
  cout<<"Printing averaged params: "<<endl;
  cout<<"C: "<<Boot_ave(C)<<" +- "<<Boot_err(C)<<endl;
  cout<<"log: "<<Boot_ave(LOG)<<" +- "<<Boot_err(LOG)<<endl;
  cout<<"A1: "<<Boot_ave(A1)<<" +- "<<Boot_err(A1)<<endl;
  cout<<"A2: "<<Boot_ave(A2)<<" +- "<<Boot_err(A2)<<endl;
  cout<<"K: "<<Boot_ave(K)<<" +- "<<Boot_err(K)<<endl;
  cout<<"F_a: "<<Boot_ave(Fa)<<" +- "<<Boot_err(Fa)<<endl;
  cout<<"D: "<<Boot_ave(D)<<" +- "<<Boot_err(D)<<endl;
  cout<<"Dm: "<<Boot_ave(Dm)<<" +- "<<Boot_err(Dm)<<endl;
  cout<<"Af0: "<<Boot_ave(Af0)<<" +- "<<Boot_err(Af0)<<endl;

  cout<<"############# PRINTING LATTICE SPACING INFOS #################"<<endl;
  
  Eigen::MatrixXd Cov(9,9); //covariance matrix of input parameters
  Eigen::VectorXd Ave(9);
  ReadBranch(0, Cov, Ave );
  cout<<"ainv0 branch 0: "<<Boot_ave(Ainv0[0])<<" +- "<<Boot_err(Ainv0[0])<<" Original from branch nr. 0: "<<Ave[3] <<" +- "<<sqrt(Cov(3,3))<<endl;
  cout<<"ainv1 branch 0: "<<Boot_ave(Ainv1[0])<<" +- "<<Boot_err(Ainv1[0])<<" Original from branch nr. 0: "<<Ave[4] <<" +- "<<sqrt(Cov(4,4))<<endl;
  cout<<"ainv2 branch 0: "<<Boot_ave(Ainv2[0])<<" +- "<<Boot_err(Ainv2[0])<<" Original from branch nr. 0: "<<Ave[5] <<" +- "<<sqrt(Cov(5,5))<<endl;

  

  
  

    
  return;
}
 
					
  






















