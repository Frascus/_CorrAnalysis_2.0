#include "../include/PionMassAnalysis_twisted.h"


using namespace std;
namespace plt = matplotlibcpp;

constexpr double kappa=2.837297;
const double MPiPhys=0.135;
const double alpha = 1.0/137.04;
const double e2 = alpha*4.0*M_PI;
const double r0 = pow(0.672/0.197,2);
const double fpi_phys = 0.1304;
const double csi_phys = pow(MPiPhys,2)/pow(4.0*M_PI*fpi_phys,2);
const int Nbranches = 8;
const bool Use_JB_distribution= false;
const int Njacks=15;
const bool UseJack=1;
const int Nboots=200;
const int nboots= Use_JB_distribution?Njacks:100;
const bool verbose=1;



class Y_tt {

public:
  Y_tt() {}
  double L, Mpi, Dm2_tot,  Mpi_err, fpi, dm2_ov_fpi2_err, fpi_err, csi;
  int ibeta;
};

class Y_tt_fp {

public:
  Y_tt_fp() {}
  double  L, Mpi, Mpi_err, fpi, fpi_err, csi;
  int ibeta;
};



class X_tt {

public:
  X_tt(const Vfloat &par) : ainv(3), Za(3) {
    if(par.size() != 17) {
      cout<<"In class X_tt, invalid call to constructor"<<endl;
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
    this->l4 = par[9];
    this->A2_fpi = par[10];
    for(int ibeta=0;ibeta<3;ibeta++) {
      this->ainv[ibeta] = par[11+ibeta];
    }
    for(int ibeta=0; ibeta<3;ibeta++) {
      this->Za[ibeta] = par[14+ibeta];
   }
  }
  
  X_tt() : ainv(3), Za(3) {}
  vector<double> ainv;
  vector<double> Za;
  double D, Dm, chir, log, log_a, A_1, F_a, F_m, A_2, l4, A2_fpi;
};

class X_tt_fp {

public:
  X_tt_fp(const Vfloat &par) : ainv(3) {
    if(par.size() != 8) {
      cout<<"In class X_tt_fp, invalid call to constructor"<<endl;
      crash("par_size: +"+to_string(par.size()));
    }
    this->fpi_chir=par[0];
    this->D=par[1];
    this->Dm=par[2];
    this->l4 = par[3];
    this->A2_fpi = par[4];  
    for(int ibeta=0;ibeta<3;ibeta++) {
      this->ainv[ibeta] = par[5+ibeta];
    }
  }
  
  X_tt_fp() : ainv(3) {}
  vector<double> ainv;
  double fpi_chir;
  double D, Dm, l4, A2_fpi;
};





void Pion_mass_analysis_twisted_adim(string CURRENT_TYPE, bool IncludeDisconnected) {

   omp_set_num_threads(1);

  
  data_t m_data, dm_exch_data, dm_hand_data, m_data_hand_run, m_twisted_data;
  if (CURRENT_TYPE=="LOCAL") { //current is local
    // m_data.Read("../datasets", "mes_contr_00", "P5P5");
    //dm_exch_data.Read("../datasets", "mes_contr_LL", "P5P5");
    m_data.Read("../datasetslocal_twisted/data", "mes_contr_M0_R0_0_M0_R0_0", "P5P5");
    dm_exch_data.Read("../datasetslocal_twisted/data", "mes_contr_M0_R0_F_M0_R0_F", "P5P5");
    if(IncludeDisconnected) {
      m_data_hand_run.Read("../datasetslocal_twisted/data", "mes_contr_M0_R0_0_M0_R0_0_handcuffs", "P5P5");
      dm_hand_data.Read("../datasetslocal_twisted/data", "handcuffs", "P5P5");
    }
  }
  else crash("CURRENT_TYPE: "+CURRENT_TYPE+ " not yet implemented. Exiting...");
  
  
  //compute masses and effective masses

  
  
  int Nens = m_data.size;  //== dm_exch_data.size .....
  cout<<"N_ens: "<<Nens<<endl;



  distr_t_list Mpi_distr_list(UseJack), Dm2_tot_distr_list(UseJack), ratio_disc_exch_mass_distr_list(UseJack), ratio_disc_exch_mass_vol_sub_distr_list(UseJack);
  distr_t_list Mpi_distr_dim_list(UseJack);
  distr_t_list Dm2_disc_ov_Mpi2(UseJack), Dm2_disc_ov_Mpi2log(UseJack), Dm2_disc(UseJack);
  distr_t_list Dm2_tot_ov_fp2_distr_list(UseJack), fp_fit_distr_list(UseJack),  Dm2_tot_ov_fp2_ren_FSEs_distr_list(UseJack);
  distr_t_list Dm2_tot_ov_fp2_FSEs_distr_list(UseJack), fp_FSEs_fit_distr_list(UseJack);
  vector<Eigen::MatrixXd> CovMatrixPion;

  for(int i=0; i < Nens; i++) {
   
    CorrAnalysis Corr(UseJack, Njacks,Nboots);
    if(m_data.Tag[i].substr(0,1) == "A") {
      Corr.Tmin = 13;
    }
    else if(m_data.Tag[i].substr(0,1) =="B") Corr.Tmin = 14;
    else Corr.Tmin= 21;
    if(m_data.Tag[i].substr(0,1)=="D") Corr.Tmax = m_data.nrows[i]/2 -4;
    else Corr.Tmax= m_data.nrows[i]/2 -4;
    Corr.Nt = m_data.nrows[i];
    boost::filesystem::create_directory("../data");
    boost::filesystem::create_directory("../data/Mpi_twisted_ov_fp2");
    string p= "local";
    boost::filesystem::create_directory("../data/Mpi_twisted_ov_fp2/"+p);

    double ll, mm, tt;

    Read_pars_from_ensemble_tag(m_data.Tag[i], mm, ll, tt);


     
    //read Za
    LatticeInfo LL("LOCAL");
    LL.LatInfo(m_data.Tag[i].substr(0,1));
    double Za = LL.Za;
    double Za_e = LL.Za_err;
    double lat_spacing_inv = LL.ainv;
    GaussianMersenne G(444363);
    distr_t Za_distr(UseJack);
    for(int nj=0;nj<Njacks;nj++) Za_distr.distr.push_back( Za + G()*Za_e/sqrt(Njacks-1));
    
   
   
    distr_t_list Pi_iso_distr = Corr.corr_t(m_data.col(0)[i], "");

    distr_t_list Mpi_eff_distr = Corr.effective_mass_t(Pi_iso_distr,  "../data/Mpi_twisted_ov_fp2/"+p+"/mass."+m_data.Tag[i]);

    
    distr_t_list Exch_distr =  Corr.corr_t(dm_exch_data.col(0)[i], "");

    distr_t_list fp_distr= Corr.decay_constant_t(4.0*mm*mm*Pi_iso_distr, "../data/Mpi_twisted_ov_fp2/"+p+"/fp."+m_data.Tag[i]);

    distr_t fp_fit= Corr.Fit_distr(fp_distr);

    distr_t Mpi_fit = Corr.Fit_distr(Mpi_eff_distr);

    distr_t csi_L = Mpi_fit*Mpi_fit/(pow(4.0*M_PI,2)*fp_fit*fp_fit);

    distr_t fp_FSEs_fit = fp_fit/(1.0- 2.0*csi_L*distr_t::f_of_distr(g1_l, Mpi_fit*ll));

 
    distr_t_list Dm_exch_eff_distr = Corr.effective_slope_t(Exch_distr, Pi_iso_distr, "");

    distr_t_list Dm_exch_eff_distr_ren= Corr.effective_slope_t(Za_distr*Za_distr*Exch_distr, Pi_iso_distr, "../data/Mpi_twisted_ov_fp2/"+p+"/dm_exch."+m_data.Tag[i]);

   
   
    distr_t_list Dm_hand_eff_distr, Mpi_eff_hand_run_distr, Hand_distr, Pi_iso_distr_hand_run;
    if(IncludeDisconnected) {
      Mpi_eff_hand_run_distr= Corr.effective_mass_t(m_data_hand_run.col(0)[i], "");
      Pi_iso_distr_hand_run = Corr.corr_t(m_data_hand_run.col(0)[i],"");
      Hand_distr = Corr.corr_t(dm_hand_data.col(0)[i], "");
      Dm_hand_eff_distr= Corr.effective_slope_t(Za_distr*Za_distr*Hand_distr, Pi_iso_distr, "../data/Mpi_twisted_ov_fp2/"+p+"/dm_hand."+m_data.Tag[i]);
      distr_t_list Dm2_hand_eff_distr = Dm_hand_eff_distr*2.0*Mpi_eff_distr;
      Print_To_File({}, {Dm2_hand_eff_distr.ave(),Dm2_hand_eff_distr.err()}, "../data/Mpi_twisted_ov_fp2/"+p+"/dm2_hand."+m_data.Tag[i], "", "");
      //cout<<"disc: "<<m_data.Tag[i]<<" "<<Corr.Fit_distr(Dm_hand_eff_distr).ave()<<" "<<Corr.Fit_distr(Dm_hand_eff_distr).err()<<endl;


      ratio_disc_exch_mass_distr_list.distr_list.push_back(Corr.Fit_distr(Dm_hand_eff_distr/Dm_exch_eff_distr));
      distr_t_list univ_sub_ratio = Za_distr*Za_distr*e2*Mpi_eff_distr*Dm_hand_eff_distr/(e2*Za_distr*Za_distr*Mpi_eff_distr*Dm_exch_eff_distr + (kappa*alpha/ll)*( Mpi_eff_distr + 2.0/ll));
      ratio_disc_exch_mass_vol_sub_distr_list.distr_list.push_back(Corr.Fit_distr(univ_sub_ratio));
      Dm2_disc_ov_Mpi2.distr_list.push_back(Corr.Fit_distr(Za_distr*Za_distr*e2*Dm_hand_eff_distr/Mpi_eff_distr));
      auto mlog= [](double x, double y) -> double { return log(x);};
      Dm2_disc_ov_Mpi2log.distr_list.push_back(Corr.Fit_distr(Za_distr*Za_distr*e2*Dm_hand_eff_distr/(Mpi_eff_distr*distr_t_list::f_of_distr_list(mlog, (Mpi_eff_distr*Mpi_eff_distr/(16.0*M_PI*M_PI*fp_distr*fp_distr))))));
      Dm2_disc.distr_list.push_back(Corr.Fit_distr(Za_distr*Za_distr*e2*Mpi_eff_distr*Dm_hand_eff_distr));
    }
    
    
    
    Mpi_distr_list.distr_list.push_back(Corr.Fit_distr( Mpi_eff_distr));
    Mpi_distr_dim_list.distr_list.push_back( lat_spacing_inv*Corr.Fit_distr(Mpi_eff_distr));
    distr_t_list Dm_tot_eff_distr=Dm_exch_eff_distr;
    
    //if(IncludeDisconnected) Dm_tot_eff_distr= Dm_exch_eff_distr-Dm_hand_eff_distr;
    if(IncludeDisconnected) {
      Dm_tot_eff_distr = Corr.effective_slope_t(Exch_distr-Hand_distr, Pi_iso_distr, "../data/Mpi_twisted_ov_fp2/"+p+"/dm_tot."+m_data.Tag[i]);
      distr_t_list Dm_tot_eff_distr_renorm= Corr.effective_slope_t(Za_distr*Za_distr*(Exch_distr-Hand_distr), Pi_iso_distr, "../data/Mpi_twisted_ov_fp2/"+p+"/dm_tot_renorm."+m_data.Tag[i]);
    }
    Dm2_tot_distr_list.distr_list.push_back ( 2.0*Corr.Fit_distr(Dm_tot_eff_distr*Mpi_eff_distr));
    Dm2_tot_ov_fp2_distr_list.distr_list.push_back(fpi_phys*fpi_phys*e2*Corr.Fit_distr(Dm_tot_eff_distr*Mpi_eff_distr)/(fp_fit*fp_fit));
    Dm2_tot_ov_fp2_ren_FSEs_distr_list.distr_list.push_back( fpi_phys*fpi_phys*(e2*Corr.Fit_distr(Za_distr*Za_distr*Dm_tot_eff_distr*Mpi_eff_distr)  + (kappa*alpha/ll)*( Corr.Fit_distr(Mpi_eff_distr) + 2.0/ll)   )/(fp_FSEs_fit*fp_FSEs_fit));
    fp_fit_distr_list.distr_list.push_back(fp_fit);
    fp_FSEs_fit_distr_list.distr_list.push_back(fp_FSEs_fit);
    Dm2_tot_ov_fp2_FSEs_distr_list.distr_list.push_back(Corr.Fit_distr(Dm_tot_eff_distr*Mpi_eff_distr)/(fp_FSEs_fit*fp_FSEs_fit));
   
    if(!Use_JB_distribution) { //resample from Gaussian distribution and push_back to CovMatrixPion
      int new_size =CovMatrixPion.size()+1;
      CovMatrixPion.resize(new_size);
      Compute_covariance_matrix(UseJack, CovMatrixPion[new_size-1], 3, Mpi_distr_list.distr_list[i].distr, Dm2_tot_distr_list.distr_list[i].distr, fp_FSEs_fit_distr_list.distr_list[i].distr);
              
    }


    //###################################### PRINT ENSEMBLE INFO    ###############################################


    cout<<"################### ANALYZED ENSEMBLE : "<<m_data.Tag[i]<<" ############################ "<<endl;
    cout<<"PRINTING INFO: "<<endl;
    cout<<"Mpi: "<<Mpi_fit.ave()<<" +- "<<Mpi_fit.err()<<endl;
    cout<<"Mpi (dim) : "<<(lat_spacing_inv*Mpi_fit).ave()<<" +- "<<(lat_spacing_inv*Mpi_fit).err()<<endl;
    cout<<"Mpi*L :"<<Mpi_fit.ave()*ll<<endl;
    cout<<"fp: "<<fp_fit.ave()<<" +- "<<fp_fit.err()<<endl;
    cout<<"fp (GL): "<<fp_FSEs_fit.ave()<<" +- "<<fp_FSEs_fit.err()<<endl;
    cout<<"ainv (GeV-1): "<<lat_spacing_inv<<endl;
    cout<<"#########################################################################################################"<<endl;

    
  }


  
  //everything is set up. Start fitting!
  LatticeInfo L_info(CURRENT_TYPE);


  //#################### FIT ANSATZ FOR PION DECAY CONSANT #####################
  auto ansatz_fp = [=](const X_tt_fp &p, const Y_tt_fp &ip) -> double {

		     double csi = ip.csi;

		     return p.fpi_chir*( 1.0 + 2.0*csi*p.l4 -2.0*csi*log(csi) +  pow(csi,2)*p.A2_fpi + p.D/pow(p.ainv[ip.ibeta],2)  + p.Dm*csi/pow(p.ainv[ip.ibeta],2));
		     
		   };

  auto meas_fp = [=](const X_tt_fp &p, const Y_tt_fp &ip) -> double {

		   return ip.fpi*p.ainv[ip.ibeta];

		 };

  auto err_fp = [=](const X_tt_fp &p, const Y_tt_fp &ip) -> double {

		  return ip.fpi_err*p.ainv[ip.ibeta];

		};

  //################### END DEFINITION FIT ANSATZ #######################



  //#####################   INITIALIZATION FIT OF PION DECAY CONSTANT ######################
  bootstrap_fit<X_tt_fp, Y_tt_fp> bf_fp(nboots);
  bf_fp.set_warmup_lev(2);

  bf_fp.Set_number_of_measurements(Nens);
  bf_fp.Set_verbosity(verbose);
  bf_fp.ansatz= ansatz_fp;
  bf_fp.measurement = meas_fp;
  bf_fp.error= err_fp;

  //add parameter for bootstrap fit of fp
  bf_fp.Add_par("fp_chir", 0.124, 0.002);
  bf_fp.Add_par("D", 1, 1e-2);
  bf_fp.Add_par("Dm", 1, 1e-2);
  bf_fp.Add_par("l4", -1.0, 0.2); //LEC determined from previous chiral fits
  bf_fp.Add_par("A2_fpi", 1.0, 1e-2);
  bf_fp.Add_prior_pars({"ainv0", "ainv1", "ainv2"});

  //Add List of parameters to be released after first minimization
  bf_fp.Fix_n_release({"ainv0", "ainv1", "ainv2"});

  //fix pars
  //bf_fp.Fix_par("A2_fpi", 0.0);
  

  //#######################   INITIALIZATION FIT OF Dm^2/fp^2  ############################
  bootstrap_fit<X_tt,Y_tt> bf(nboots);
  //bf.set_warmup_lev(2);
  bf.Set_number_of_measurements(Nens);
  bf.Set_verbosity(verbose);


    //Add parameters
    bf.Add_par("chir", 0.8, 1.0e-2);
    bf.Add_par("log", 5.7 , 0.7);
    bf.Add_par("A_1", -5.7, 0.02);
    bf.Add_par("D", 1, 1e-2);
    bf.Add_par("Dm", 1, 1e-2);
    bf.Add_par("F_a", 1, 0.1);
    bf.Add_par("lg_a", 0.1, 1e-3);
    bf.Add_par("F_m", 3.0, 0.1);
    bf.Add_par("A_2", 1.0, 0.01);
    bf.Add_par("l4", 3, 0.3);
    bf.Add_par("A2_fpi", 1.0, 1e-2);
    bf.Add_prior_pars({"ainv0", "ainv1", "ainv2","Za0", "Za1", "Za2"});
    
    //Fix some parameters to make test
    //bf.Fix_par("F_m",0.0);
    bf.Fix_par("F_a",0.0);
    //bf.Fix_par("Dm",0.0);
    //bf.Fix_par("D",0.0);
    bf.Fix_par("lg_a", 0.0);
    bf.Fix_par("A_2", 0.0);
    //bf.Fix_par("log", 0.0);
    bf.Fix_par("A2_fpi",0.0);

    //Add List of parameters to be released after first minimization
    bf.Fix_n_release({"ainv0", "ainv1", "ainv2", "Za0", "Za1", "Za2"});




  
  
  Vfloat L, m_l, T;
  Read_pars_from_ensemble_tag(m_data.Tag, m_l, L, T);


  
  


  //define fitting function for Dm^2/fp^2

  //############################################################################################################

  bf.ansatz =  [=](const X_tt &p, const Y_tt &ip) -> double {


    double ainv = p.ainv[ip.ibeta];
    double Lfr = ip.L;
    double Mp = ainv*ip.Mpi;
    double Mp2 = pow(Mp,2);
    double aMpi = ip.Mpi;
    double a=1.0/ainv;
    double L = ip.L/ainv;
    double fp_dim= ainv*ip.fpi;
    double SD_FVE = (L>0.0)?(e2/3.0)*(Mp/pow(fp_dim,2))*r0/pow(L,3)*(p.F_m+ p.F_a*(pow(a,2))):0.0;
    double FVE_universal = (L>0.0)?(kappa*alpha/Lfr)*( ip.Mpi + 2.0/Lfr):0.0;
    double log_par= p.log;
    double csi = ip.csi;
     

    double fitted_value = e2*( p.chir - (log_par + p.log_a/pow(ainv,2))*csi*log(csi) + p.A_1*csi + p.A_2*pow(csi,2)); 
     
    
    double FP= 1.0 + 2*csi*p.l4  -2.0*csi*log(csi) + pow(csi,2)*p.A2_fpi;

    return pow(fpi_phys,2)*(fitted_value/pow(FP,2) + (p.Dm*csi+ p.D)*a*a  + SD_FVE - FVE_universal/pow(ip.fpi,2));
    
  };
  

  bf.measurement = [=](const X_tt& p,const Y_tt& ip) -> double {

    double Za = p.Za[ip.ibeta];

            
    double m1 = pow(fpi_phys,2)*(e2/2)*pow(Za,2)*ip.Dm2_tot/pow(ip.fpi,2);

    return m1;
  };
  

  bf.error =  [=](const X_tt& p,const  Y_tt &ip) -> double {

    double Za = p.Za[ip.ibeta];
       
    double m1 = pow(fpi_phys,2)*(e2/2)*pow(Za,2)*ip.dm2_ov_fpi2_err;

    return m1;
  };

 
  
  
  //############################################################################################################

  //init random Number Generators

  GaussianMersenne GM(54353);
  RandomMersenne RM(23423, UseJack?(Njacks-1):(Nboots-1));

  Eigen::MatrixXd CovMatrixInput(9,9); //covariance matrix of input parameters
  Eigen::VectorXd Ave_input_parameters(9);


  vector<boot_fit_data<X_tt_fp>> Bt_fit_fp;
  vector<boot_fit_data<X_tt>> Bt_fit;


  //we store here all resampled fp data

  
  //we store here all resampled fp data. 
  vector<vector<vector<Y_tt_fp>>> all_ens_fp(Nbranches);
  for(auto &all_ens_br : all_ens_fp) {
    all_ens_br.resize(nboots);
    for(auto & all_ens_ib: all_ens_br) {
      all_ens_ib.resize(Nens);
      for(int iens=0; iens<Nens;iens++) {
	all_ens_ib[iens].L=L[iens];
	if(m_data.Tag[iens].substr(0,1) == "A") all_ens_ib[iens].ibeta=0;
	else if(m_data.Tag[iens].substr(0,1) == "B") all_ens_ib[iens].ibeta=1;
	else all_ens_ib[iens].ibeta=2;
	all_ens_ib[iens].Mpi_err = Mpi_distr_list.err(iens);
	all_ens_ib[iens].fpi_err= fp_FSEs_fit_distr_list.err(iens);
      }
    }
  }


  


  

  //we store here all resampled DM^2/fp^2 data. 
  vector<vector<vector<Y_tt>>> all_ens(Nbranches);
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
	all_ens_ib[iens].dm2_ov_fpi2_err= Dm2_tot_ov_fp2_FSEs_distr_list.err(iens);
	all_ens_ib[iens].fpi_err = fp_FSEs_fit_distr_list.err(iens);

      }
    }
  }


  for(int ibranch=0; ibranch < Nbranches; ibranch++) {


    
    ReadBranch(ibranch, CovMatrixInput, Ave_input_parameters );


    //###############################   FIT fp and Dm^2/fp^2 #####################################

    //clear all measurements
    //fp
    bf_fp.Clear_priors();
    bf_fp.Clear_input_pars();
    //Dm^2/fp^2
    bf.Clear_priors();
    bf.Clear_input_pars();




    //generate bootstrap data
    for(int iboot=0; iboot<nboots; iboot++) {
      bf.ib= &iboot;
      bf_fp.ib = &iboot;
      Vfloat lat_input = Covariate(CovMatrixInput, Ave_input_parameters, GM);
      double Za0= gauss(L_info.Retrieve_Za("A", ibranch),GM);
      double Za1= gauss(L_info.Retrieve_Za("B", ibranch),GM);
      double Za2= gauss(L_info.Retrieve_Za("D", ibranch) ,GM);
     
      bf_fp.Append_to_prior("ainv0", lat_input[3], sqrt(CovMatrixInput(3,3)));
      bf_fp.Append_to_prior("ainv1", lat_input[4], sqrt(CovMatrixInput(4,4)));
      bf_fp.Append_to_prior("ainv2", lat_input[5], sqrt(CovMatrixInput(5,5)));
      bf.Append_to_prior("ainv0", lat_input[3], sqrt(CovMatrixInput(3,3)));
      bf.Append_to_prior("ainv1", lat_input[4], sqrt(CovMatrixInput(4,4)));
      bf.Append_to_prior("ainv2", lat_input[5], sqrt(CovMatrixInput(5,5)));
      bf.Append_to_prior("Za0", Za0, L_info.Retrieve_Za("A", ibranch).second);
      bf.Append_to_prior("Za1", Za1, L_info.Retrieve_Za("B", ibranch).second);
      bf.Append_to_prior("Za2", Za2, L_info.Retrieve_Za("D", ibranch).second);
      // bf.Append_to_prior("l4", l4_ib, l4_err);



      int k= iboot;
      for(int imeas=0; imeas <Nens; imeas++) {
	if(!Use_JB_distribution) {
	  Eigen::VectorXd vec(CovMatrixPion[imeas].rows());
	 
	  vec<<Mpi_distr_list.ave(imeas),Dm2_tot_distr_list.ave(imeas),fp_FSEs_fit_distr_list.ave(imeas);
	  Vfloat res_meas(3,0.0);
	 
	  res_meas= Covariate(CovMatrixPion[imeas],vec, GM);
	  all_ens_fp[ibranch][iboot][imeas].Mpi= res_meas[0];
	  all_ens_fp[ibranch][iboot][imeas].fpi = res_meas[2];
	  all_ens_fp[ibranch][iboot][imeas].csi = pow(res_meas[0],2)/(pow(4.0*M_PI*res_meas[2],2));
	  all_ens[ibranch][iboot][imeas].Mpi= res_meas[0];
	  all_ens[ibranch][iboot][imeas].Dm2_tot = res_meas[1];
	  all_ens[ibranch][iboot][imeas].fpi = res_meas[2];
	  all_ens[ibranch][iboot][imeas].csi = pow(res_meas[0],2)/(pow(4.0*M_PI*res_meas[2],2));
	 	 
	 	 
	}
      
	else {
	  double Mp_boot, fpi_boot, Dm2_tot_boot;
	  Mp_boot=  Mpi_distr_list.distr_list[imeas].distr[k];
	  fpi_boot=  fp_FSEs_fit_distr_list.distr_list[imeas].distr[k];
	  Dm2_tot_boot = Dm2_tot_distr_list.distr_list[imeas].distr[k];
	  all_ens_fp[ibranch][iboot][imeas].Mpi = Mp_boot;
	  all_ens_fp[ibranch][iboot][imeas].fpi = fpi_boot;
	  all_ens_fp[ibranch][iboot][imeas].csi = pow(Mp_boot,2)/pow(4.0*M_PI*fpi_boot,2);
	  all_ens[ibranch][iboot][imeas].Mpi = Mp_boot;
	  all_ens[ibranch][iboot][imeas].Dm2_tot = Dm2_tot_boot;
	  all_ens[ibranch][iboot][imeas].fpi = fpi_boot;
	}
      }

    }

    //fit  fp
    cout<<"#################FITTING FP####################"<<endl;
    bf_fp.Append_to_input_par(all_ens_fp[ibranch]);
    Bt_fit_fp.push_back(bf_fp.Perform_bootstrap_fit());

    


    //###########  SET INITIAL VALUE FOR l4 and log_fp USING VALUES FROM FP FIT ###############################

    Vfloat A2fpi_guess_distr, l4_guess_distr;
    double A2fpi_guess, A2fpi_guess_err, l4_guess, l4_guess_err;
    for(int iboot=0;iboot<nboots;iboot++) {
      A2fpi_guess_distr.push_back( Bt_fit_fp[ibranch].par[iboot].A2_fpi);
      l4_guess_distr.push_back( Bt_fit_fp[ibranch].par[iboot].l4);
    }

    A2fpi_guess = Boot_ave(A2fpi_guess_distr);
    A2fpi_guess_err= Boot_err(A2fpi_guess_distr, Use_JB_distribution);
    l4_guess= Boot_ave(l4_guess_distr);
    l4_guess_err = Boot_err(l4_guess_distr, Use_JB_distribution);


    //Set parameters value  log_fpi and l4 for Dm^2/fp^2 fit 
    bf.Set_par_val("l4", l4_guess, l4_guess_err);
    //bf.Set_par_val("A2_fpi", A2fpi_guess, A2fpi_guess_err);


    //##############################  ######################################################################à



    //###########################   FIT Dm^2/fp^2 ##########################################################
  

    cout<<"################ FITTING Dm^2/fp^2  #################################"<<endl;
    bf.Append_to_input_par(all_ens[ibranch]);
    Bt_fit.push_back(bf.Perform_bootstrap_fit());


    //########################### END FIT Dm^2/fp^2 #########################################################
    
  }
  
 
  //print the data
  
  //define lambda functions

  auto FVE = [](double L, double Mp, double fpi) {
    double dm= (kappa*alpha/L)*( Mp + 2.0/L);
    return pow(fpi_phys,2)*dm/pow(fpi,2);
    
  };

  auto SD_FVE = [](double L, double Mp, double fpi,  double ainv,  double F_a, double F_m) {

		  double Ldim = L/ainv;
		  double Mpdim= Mp*ainv;
		  double fpi_dim = fpi*ainv;
		  double SDE =(e2/3.0)*(Mpdim/pow(fpi_dim,2))*(1.0/pow(Ldim,3))*r0*(F_m +  F_a*(1.0/pow(ainv,2)));
		  return SDE*pow(fpi_phys,2);
  };


  auto cont_ansatz = [](X_tt p, double csi) {


    double fitted_value = e2*( p.chir - (p.log)*csi*log(csi) + p.A_1*csi + p.A_2*pow(csi,2)) ;

     
    
    double FP= 1.0 + 2.0*csi*p.l4  -2.0*csi*log(csi) + pow(csi,2)*p.A2_fpi;

 

    return pow(fpi_phys,2)*fitted_value/pow(FP,2);

     
  };


  auto fp_cont = [](X_tt_fp p, double csi) {

		   return p.fpi_chir*( 1.0 + 2.0*csi*p.l4 -2.0*csi*log(csi) + pow(csi,2)*p.A2_fpi); 
		 };

  
  Vfloat Csi_fit_points(600);
  Vfloat vols(600);
  for(unsigned int m=0; m<Csi_fit_points.size();m++) Csi_fit_points[m] = 0.0005*m+ 0.001;
  for(unsigned int vol=0;vol<vols.size();vol++) vols[vol] = 1.0/(vol+18);
  VVVfloat SD_subtracted_data, Cs_measured, raw_data;
  VVVfloat Univ_subtracted_data;
  VVVVfloat Fitted_func, Fitted_func_fp;
  VVfloat Physical_point, Violation_SU2, Physical_point_fp, Chiral_fp;
  VVVfloat A40_slice;
  cascade_resize(SD_subtracted_data, Vint{bf.Get_number_of_measurements(), Nbranches, nboots});
  cascade_resize(Univ_subtracted_data, Vint{bf.Get_number_of_measurements(), Nbranches, nboots});
  cascade_resize(Cs_measured, Vint{bf.Get_number_of_measurements(), Nbranches, nboots});
  cascade_resize(Fitted_func, Vint{4, (int)Csi_fit_points.size(), Nbranches, nboots});
  cascade_resize(Fitted_func_fp, Vint{4, (int)Csi_fit_points.size(), Nbranches, nboots});
  cascade_resize(Physical_point, Vint{Nbranches, nboots});
  cascade_resize(Physical_point_fp, Vint{Nbranches, nboots});
  cascade_resize(Chiral_fp, Vint{Nbranches, nboots});
  Violation_SU2=Physical_point;
  cascade_resize(A40_slice, Vint{ (int)vols.size(), Nbranches,nboots});
  cascade_resize(raw_data, Vint{bf.Get_number_of_measurements(), Nbranches, nboots});


  
  for(int ibranch=0; ibranch<Nbranches;ibranch++) {
    for(int iboot=0; iboot<nboots;iboot++) {
      X_tt P = Bt_fit[ibranch].par[iboot];
      X_tt_fp P_fp= Bt_fit_fp[ibranch].par[iboot];
      for(unsigned int ivol=0; ivol<vols.size();ivol++) {
	Y_tt new_p;
	new_p.L = 1.0/vols[ivol];
	new_p.ibeta=0;
	new_p.Mpi= 0.318/P.ainv[0];
	new_p.fpi = fpi_phys/P.ainv[0];
	A40_slice[ivol][ibranch][iboot] = bf.ansatz(P,new_p) + FVE(new_p.L, new_p.Mpi, new_p.fpi);
      }
      for(int imeas=0; imeas< bf.Get_number_of_measurements(); imeas ++) {
	Y_tt E = all_ens[ibranch][iboot][imeas];
	SD_subtracted_data[imeas][ibranch][iboot] = bf.measurement(P, E) + FVE(E.L, E.Mpi, E.fpi) - SD_FVE(E.L, E.Mpi, E.fpi, P.ainv[E.ibeta], P.F_a, P.F_m);
	Univ_subtracted_data[imeas][ibranch][iboot] = bf.measurement(P,E)+ FVE(E.L, E.Mpi, E.fpi);
	raw_data[imeas][ibranch][iboot] = bf.measurement(P,E);
	Cs_measured[imeas][ibranch][iboot] = pow(E.Mpi,2)/pow(4.0*M_PI*E.fpi,2);
      }
      for(int ibeta=0; ibeta < 4;ibeta++) {
	for(unsigned int m=0; m<Csi_fit_points.size(); m++){
	  if(ibeta<3) {
	    Y_tt point;
	    Y_tt_fp point_fp;
	    point.L = -1.0;
	    point.ibeta=ibeta;
	    point.Mpi= 1.0;
	    point.fpi= 1.0;
	    point.csi = Csi_fit_points[m];
	    point_fp.L= -1.0;
	    point_fp.ibeta=ibeta;
	    point_fp.csi= Csi_fit_points[m];
	    Fitted_func[ibeta][m][ibranch][iboot] = bf.ansatz(P,point);
	    Fitted_func_fp[ibeta][m][ibranch][iboot] = bf_fp.ansatz(P_fp, point_fp)/P_fp.ainv[ibeta];
	  }
	  else {;
	    Fitted_func[ibeta][m][ibranch][iboot] = cont_ansatz(P, Csi_fit_points[m]) ;
	    Fitted_func_fp[ibeta][m][ibranch][iboot] = fp_cont(P_fp, Csi_fit_points[m]);
	}
	}
      }
      //Add physical point
	Physical_point[ibranch][iboot] = cont_ansatz(P, csi_phys);
	Physical_point_fp[ibranch][iboot] = fp_cont(P_fp,csi_phys);
	Chiral_fp[ibranch][iboot] = P_fp.fpi_chir;
	Violation_SU2[ibranch][iboot] = P.log - (3.0+(4.0*P.chir)); 
      }
    }


   

  Vfloat SD_subtracted_val, SD_subtracted_err, CSI, SD_subtracted_dim_val, SD_subtracted_dim_err, raw_data_val, raw_data_err;
  Vfloat Univ_subtracted_val, Univ_subtracted_err;
  for(int imeas=0;imeas<bf.Get_number_of_measurements();imeas++) {
    SD_subtracted_val.push_back( Boot_ave(SD_subtracted_data[imeas]));
    SD_subtracted_err.push_back( Boot_err(SD_subtracted_data[imeas], Use_JB_distribution));
    Univ_subtracted_val.push_back(Boot_ave(Univ_subtracted_data[imeas]));
    Univ_subtracted_err.push_back(Boot_err(Univ_subtracted_data[imeas], Use_JB_distribution));
    raw_data_val.push_back(Boot_ave(raw_data[imeas]));
    raw_data_err.push_back(Boot_err(raw_data[imeas], Use_JB_distribution));
    CSI.push_back( Boot_ave(Cs_measured[imeas]));
  }
  
  VVfloat Fitted_func_val, Fitted_func_err, Fitted_func_fp_val, Fitted_func_fp_err;
  Vfloat A40_slice_val, A40_slice_err, A40_raw_val, A40_raw_err;
 
  cascade_resize(Fitted_func_val, {4, (int)Csi_fit_points.size()});
  cascade_resize(Fitted_func_err, {4, (int)Csi_fit_points.size()});
  cascade_resize(Fitted_func_fp_val, {4, (int)Csi_fit_points.size()});
  cascade_resize(Fitted_func_fp_err, {4, (int)Csi_fit_points.size()});

  for(int ibeta=0;ibeta<4;ibeta++) {
    for(unsigned int m=0; m<Csi_fit_points.size();m++) {
      if(ibeta<3) {
	
	Fitted_func_val[ibeta][m] = Boot_ave(Fitted_func[ibeta][m]);
	Fitted_func_err[ibeta][m] = Boot_err(Fitted_func[ibeta][m], Use_JB_distribution);
		
	Fitted_func_fp_val[ibeta][m] = Boot_ave(Fitted_func_fp[ibeta][m]);
	Fitted_func_fp_err[ibeta][m] = Boot_err(Fitted_func_fp[ibeta][m], Use_JB_distribution);
      }
      else {
	for(int ibranch=0; ibranch<Nbranches;ibranch++)
	  {
	    ofstream print_boot_res("../data/Mpi_twisted_ov_fp2/local/correlated_disc/m_"+to_string(m)+"_ibr_"+to_string(ibranch)+"_1.dat");
            for(int iboot=0; iboot<nboots;iboot++) print_boot_res<<Csi_fit_points[m]<<setw(20)<<Fitted_func[3][m][ibranch][iboot]<<endl;
	    print_boot_res.close();
	  }
	Fitted_func_val[ibeta][m] = Boot_ave(Fitted_func[ibeta][m]);
	Fitted_func_err[ibeta][m] = Boot_err(Fitted_func[ibeta][m], Use_JB_distribution);	
	Fitted_func_fp_val[ibeta][m] = Boot_ave(Fitted_func_fp[ibeta][m]);
	Fitted_func_fp_err[ibeta][m] = Boot_err(Fitted_func_fp[ibeta][m], Use_JB_distribution);
	
      }
    }
  }

  for(auto &A40_boot : A40_slice) {
    A40_slice_val.push_back( Boot_ave(A40_boot));
    A40_slice_err.push_back( Boot_err(A40_boot, Use_JB_distribution));
  }
  //plot the result

  Vfloat vol_A40, meas_A40, err_meas_A40;
  for(int iens=0; iens<Nens;iens++)
    if(m_data.Tag[iens].substr(0,3)=="A40") {
      vol_A40.push_back( 1.0/L[iens]);
      meas_A40.push_back(Univ_subtracted_val[iens]);
      err_meas_A40.push_back(Univ_subtracted_err[iens]);
      A40_raw_val.push_back(raw_data_val[iens]);
      A40_raw_err.push_back(raw_data_err[iens]);
    }


  //

      

  string print_path = (CURRENT_TYPE=="CONSERVED")?"conserved":"local";
  string exch_or_tot = IncludeDisconnected?"tot":"exch";
  boost::filesystem::create_directory("../plots");
  boost::filesystem::create_directory("../plots/Mpi_twisted_ov_fp2");



   
  //save data in files
  boost::filesystem::create_directory("../data");
  boost::filesystem::create_directory("../data/Mpi_twisted_ov_fp2/"+print_path);
  Print_To_File({}, {Csi_fit_points, Fitted_func_val[0], Fitted_func_err[0],  Fitted_func_val[1], Fitted_func_err[1], Fitted_func_val[2], Fitted_func_err[2], Fitted_func_val[3], Fitted_func_err[3]} , "../data/Mpi_twisted_ov_fp2/"+print_path+"/Fitted_func_"+exch_or_tot+".dat", "OUT", "#Csi     beta=1.90      beta=1.95      beta=2.10           cont");

    Print_To_File({}, {Csi_fit_points, Fitted_func_fp_val[0], Fitted_func_fp_err[0],  Fitted_func_fp_val[1], Fitted_func_fp_err[1], Fitted_func_fp_val[2], Fitted_func_fp_err[2], Fitted_func_fp_val[3], Fitted_func_err[3]} , "../data/Mpi_twisted_ov_fp2/"+print_path+"/Fitted_func_fp.dat", "OUT", "#Csi     beta=1.90      beta=1.95      beta=2.10           cont_fp");
 
 
    Print_To_File(m_data.Tag, { CSI, Mpi_distr_dim_list.ave(), SD_subtracted_val, SD_subtracted_err}, "../data/Mpi_twisted_ov_fp2/"+print_path+"/sd_subtracted_data_"+exch_or_tot+".dat", "", "#Ens       #csi  #Mpi        #Mpi^2_+ - Mpi^2_0      #err");
 
    Print_To_File(m_data.Tag, { CSI, Mpi_distr_dim_list.ave(), Univ_subtracted_val, Univ_subtracted_err, raw_data_val, raw_data_err }, "../data/Mpi_twisted_ov_fp2/"+print_path+"/data_"+exch_or_tot+".dat", "", "#Ens  #csi      MP     Univ_sub       Raw_data");
 
  Print_To_File({}, {vol_A40, meas_A40, err_meas_A40, A40_raw_val, A40_raw_err}, "../data/Mpi_twisted_ov_fp2/"+print_path+"/A40_"+exch_or_tot+".dat","", "#1/L         univ_meas            univ_err       raw_meas      raw_err");
  string command = "echo "+to_string_with_precision(csi_phys,8)+"\t\t"+to_string_with_precision(Boot_ave(Physical_point),8)+"\t\t"+to_string_with_precision(Boot_err(Physical_point, Use_JB_distribution),8)+" > ../data/Mpi_twisted_ov_fp2/"+print_path.c_str()+"/Phys_val_"+exch_or_tot+".dat";
  system(command.c_str());

  if(IncludeDisconnected) {
    Print_To_File(m_data.Tag, {CSI, ratio_disc_exch_mass_distr_list.ave(), ratio_disc_exch_mass_distr_list.err(), ratio_disc_exch_mass_vol_sub_distr_list.ave(), ratio_disc_exch_mass_vol_sub_distr_list.err()}, "../data/Mpi_twisted_ov_fp2/"+print_path+"/ratio_disc_exch.data","","#Ens   Mp   ratio    ratio_err   ratio_UNIV_SUB   ratio_UNIV_SUB_err        ");

    Print_To_File(m_data.Tag, {CSI, Mpi_distr_dim_list.ave(),  Dm2_disc.ave(), Dm2_disc.err()}, "../data/Mpi_twisted_ov_fp2/"+print_path+"/dm2_disc.data","","#Ens Csi   Mp   val err    ");
    Print_To_File(m_data.Tag, {CSI, Mpi_distr_dim_list.ave(),  fp_fit_distr_list.ave(), fp_fit_distr_list.err(), fp_FSEs_fit_distr_list.ave(), fp_FSEs_fit_distr_list.err(),  Dm2_tot_ov_fp2_distr_list.ave(), Dm2_tot_ov_fp2_distr_list.err(),  Dm2_tot_ov_fp2_ren_FSEs_distr_list.ave(),  Dm2_tot_ov_fp2_ren_FSEs_distr_list.err() }, "../data/Mpi_twisted_ov_fp2/"+print_path+"/dm2_tot_ov_fp2.data", "", "#Ens Csi  Mp fp fp_err fp_FSEs  fp_FSEs_err val err  val_ren_UFSES_sub   err");

  }
  

  cout<<"Printing chiral fit pars"<<endl;
  cout<<"fpi (Physical point): "<<Boot_ave(Physical_point_fp)<<" +- "<<Boot_err(Physical_point_fp, Use_JB_distribution)<<endl;
  cout<<"f0: "<<Boot_ave(Chiral_fp)<<" +- "<<Boot_err(Chiral_fp,Use_JB_distribution)<<endl;
  cout<<"Physical Pion mass difference squared [Mev2]: "<< Boot_ave(Physical_point)*1.0e+6<<" +-   "<<Boot_err(Physical_point, Use_JB_distribution)*1.0e+6<<endl;
  cout<<"SU2 ChPT violation: "<<Boot_ave(Violation_SU2)<<"    "<<Boot_err(Violation_SU2, Use_JB_distribution)<<endl;
  cout<<"On single branches: "<<endl;
  for(int ibr=0; ibr<Nbranches;ibr++) cout<<"Branch: "<<ibr<<": "<<Boot_ave(Physical_point[ibr])*1e+6<<" +-  "<<Boot_err(Physical_point[ibr], Use_JB_distribution)*1e+6<<endl;
  cout<<"Average chi2:"<<endl;
  for(int ibr=0;ibr<Nbranches;ibr++) cout<<"Branch: "<<ibr<<" chi2 = "<<accumulate(Bt_fit[ibr].chi2.begin(), Bt_fit[ibr].chi2.end(), 0.0)/Bt_fit[ibr].chi2.size()<<endl;

  
  

    
  return;
}
 
					
  






















