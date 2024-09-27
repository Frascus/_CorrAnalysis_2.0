#include "../include/scale_setting.h"

using namespace std;
const double Mp_phys=0.1350;
const double Mp_phys_err =  0.0002e-16;
const double fp_phys= 0.1305;
const double fp_phys_err = 0.0002e-16;
const double csi_phys= pow(Mp_phys/(4.0*M_PI*fp_phys),2);
//const double csi_phys_err =  0.000029;
const double Xpi_phys= pow(fp_phys*pow(Mp_phys,4),1.0/5.0);
const double l1ph= -0.4; //-0.4
const double l2ph= 4.3; //4.3
const double l3ph= 3.2; //3.2
const double l4ph= 4.4; //4.4
const double s0= 2.0-M_PI/2.0;
const double s1 = M_PI/4.0 - 0.5;
const double s2 = 0.5 - M_PI/8.0;
const double s3 = 3.0*M_PI/16.0 - 0.5;
bool scale_setting_fit_verbosity=1;
bool Use_w0_Xpi= true;
const double fm_to_iGev=  1.0/0.197327;
const double guess_w0 = 0.17383*fm_to_iGev;
const double guess_w0_fp_phys = guess_w0*fp_phys;
const double guess_w0_Xp_phys = guess_w0*Xpi_phys;
const double guess_w0_X = (Use_w0_Xpi)?guess_w0_Xp_phys:guess_w0_fp_phys;
const double corr_fact = 0.75*0.21;
const double corr_fact_err = 0.75*0.05;
const bool Enable_A2m_cont=false;
const bool Use_fpi_or_Xpi= true;



class sc_ipar {

public:
  sc_ipar() : Mp(0.0), w0_t_fp(0.0), w0_t_Xpi(0.0), fp(0.0), Xpi(0.0), csi(0.0), w0(0.0), L(0.0), csi_ph(0.0), Mp_err(1.0), w0_t_fp_err(1.0), w0_t_Xpi_err(1.0), fp_err(1.0), Xpi_err(1.0), csi_err(1.0), w0_err(1.0) {}

  
  double Mp;
  double w0_t_fp;
  double w0_t_Xpi;
  double fp;
  double Xpi;
  double csi;
  double w0;
  double L;
  double csi_ph;
 

  double Mp_err;
  double w0_t_fp_err;
  double w0_t_Xpi_err;
  double fp_err;
  double Xpi_err;
  double csi_err;
  double w0_err;
   
};

class sc_fpar {

public:
  sc_fpar() {}
  sc_fpar(const Vfloat &par)  {
    if((signed)par.size() != 7) crash("In class sc_fpar constructor sc_fpar(vector<double>) called with vector.size != 7 ");
    w0X=par[0];
    Am=par[1];
    A2m = par[2];
    D2=par[3];
    D2m = par[4];
    D4 = par[5];
    Al = par[6];
  }

  double w0X, Am, A2m, D2, D2m, D4, Al;
 
};

class sc_fp_ipar {

public:
  sc_fp_ipar() : Mp(0.0), fp(0.0), Xpi(0.0), csi(0.0),  L(0.0), csi_ph(0.0), Mp_err(1.0),  fp_err(1.0), Xpi_err(1.0), csi_err(1.0) {}

  
  double Mp;
  double fp;
  double Xpi;
  double csi;
  double L;
  double csi_ph;
  bool Is_B;

  double Mp_err;
  double fp_err;
  double Xpi_err;
  double csi_err;
     
};

class sc_fp_fpar {

public:
  sc_fp_fpar() {}
  sc_fp_fpar(const Vfloat &par)  {
    if((signed)par.size() != 9) crash("In class sc_fp_fpar constructor sc_fp_fpar(vector<double>) called with vector.size != 9 ");
    afp_B=par[0];
    Am_B=par[1];
    A2m_B = par[2];
    Al_B = par[3];
    afp_A = par[4];
    Am_A = par[5];
    A2m_A = par[6];
    Al_A = par[7];
    Al = par[8];
  }

  double afp_B, Am_B, A2m_B, Al_B, afp_A, Am_A, A2m_A, Al_A, Al;
 
};


class rel_lat_cl {

public:
  rel_lat_cl()  {}
  
  double Y,Y_err, X;
};

class fpar_a2 {

public:
  fpar_a2() {}
  fpar_a2(const Vfloat &par) {
    if((signed)par.size() != 2) crash("In Obs_extrapolation_meson_mass, class fpar_a2 has size different from two.");
    a0=par[0];
    D=par[1];
  }

  double a0,D;

};



void Determine_scale_from_w0_fp(const distr_t_list &Mpi, const distr_t_list &fpi, const distr_t_list &w0, const vector<double> &L_list, const vector<string>  &Ensemble_tags,distr_t &a_distr_A, distr_t &a_distr_B, distr_t &a_distr_C, distr_t &a_distr_D, bool UseJack, bool Use_three_finest) {


  //init Gaussian number generator
  GaussianMersenne GM(76795699);

  
  int Nens = L_list.size();
  int Nens_total = (signed)Ensemble_tags.size();
  if(Nens==0) crash("in scale_setting.cpp, zero ensemble provided");
  int Njacks= Mpi.distr_list[0].size();

  //sanity checks
  if( Mpi.size() != fpi.size() )  crash("Number of ensembles in Mpi and fpi list, do not coincide");
 
  
  //List of lambda functions needed

  auto FVE = [](double x) -> double {return (1.0/pow(x,1.5))*exp(-x);};
  auto LOG = [](double x) { return log(x);};
  auto pow_1_5 = [](double x) -> double { return pow(x, 1.0/5.0);};
  auto pow_4 = [](double x) -> double { return pow(x, 4.0);};
  auto pow_2 = [](double x) -> double { return pow(x, 2.0);};
  auto pow_3_2 = [](double x) -> double { return pow(x, 3.0/2.0);};
  auto EXP = [](double x) -> double {return exp(x);};
  auto SQRT = [](double x) -> double { return sqrt(x);};



  //###########################

 



  //Output directories

  boost::filesystem::create_directory("../data/scale_setting");
  
  //##########################


  //Generate fake distribution for csi_phys, Mpi , fp, Xpi
  
  distr_t Mpi_phys_distr(UseJack);
  distr_t fpi_phys_distr(UseJack);
  distr_t csi_phys_distr(UseJack);
  distr_t Xpi_phys_distr(UseJack);
  
  for(int ij=0;ij<Njacks;ij++) {
    double Mfak = Mp_phys + GM()*Mp_phys_err/sqrt(Njacks -1.0);
    double Ffak = fp_phys + GM()*fp_phys_err/sqrt(Njacks -1.0);

    Mpi_phys_distr.distr.push_back( Mfak);
    fpi_phys_distr.distr.push_back(Ffak);
    csi_phys_distr.distr.push_back( pow(Mfak/(4.0*M_PI*Ffak),2.0));
    Xpi_phys_distr.distr.push_back( pow( Ffak*pow(Mfak,4.0), 1.0/5.0));
  }

 
  distr_t_list fp_CDH(UseJack), Mpi_CDH(UseJack), Xpi(UseJack), Csi_CDH(UseJack);
  
  //compute CDH corrections to Mpi and fpi
  for(int iens=0; iens<Nens;iens++) {
    //GL and CDH correction for Mpi and fp
    distr_t csi_L = Mpi.distr_list[iens]*Mpi.distr_list[iens]/(pow(4.0*M_PI,2)*fpi.distr_list[iens]*fpi.distr_list[iens]);
    distr_t g1 = distr_t::f_of_distr(g1_l, Mpi.distr_list[iens]*L_list[iens]);
    distr_t g2 = distr_t::f_of_distr(g2_l, Mpi.distr_list[iens]*L_list[iens]);
    distr_t log_l = log(csi_phys) - distr_t::f_of_distr(LOG, csi_L);
    fp_CDH.distr_list.push_back(fpi.distr_list[iens]/(1.0 -2.0*csi_L*g1 +2.0*csi_L*csi_L*( (Cf1(l1ph,l2ph,l3ph,l4ph) + Sf1(s0,s1,s2,s3) + Cf1_log()*log_l)*g1 + (Cf2(l1ph,l2ph,l3ph,l4ph) + Sf2(s0,s1,s2,s3) + Cf2_log()*log_l)*g2)));
    Mpi_CDH.distr_list.push_back(Mpi.distr_list[iens]/(1.0 + 0.5*csi_L*g1 - csi_L*csi_L*( (Cm1(l1ph,l2ph,l3ph,l4ph) + Sm1(s0,s1,s2,s3) + Cm1_log()*log_l)*g1 + (Cm2(l1ph,l2ph,l3ph,l4ph) + Sm2(s0,s1,s2,s3) + Cm2_log()*log_l)*g2)));
    Csi_CDH.distr_list.push_back( Mpi_CDH.distr_list[iens]*Mpi_CDH.distr_list[iens]/(16.0*M_PI*M_PI*fp_CDH.distr_list[iens]*fp_CDH.distr_list[iens]));
    Xpi.distr_list.push_back( distr_t::f_of_distr(pow_1_5, fpi.distr_list[iens]*distr_t::f_of_distr(pow_4, Mpi.distr_list[iens])));
  }

  //Print_To_File
  Print_To_File({Ensemble_tags}, {L_list, Mpi_CDH.ave(), Mpi_CDH.err(), fp_CDH.ave(), fp_CDH.err(), Xpi.ave(), Xpi.err(), Csi_CDH.ave(), Csi_CDH.err()}, "../data/scale_setting/pion_obs.dat", "", "#Ens L Mpi fp Xpi  Csi");



  //Perform extrapolation of w0*Xpi or w0*fp

   bootstrap_fit<sc_fpar,sc_ipar> bf(Njacks);
   bf.set_warmup_lev(1);
   bf.Set_number_of_measurements(Nens);
   bf.Set_verbosity(scale_setting_fit_verbosity);
   bf.Set_print_path("chi2_scale_setting");


   bf.Add_par("w0X", 1.0, 0.001);
   bf.Add_par("Am", -3.0, 0.1);
   bf.Add_par("A2m", 2.0, 0.1);
   bf.Add_par("D2", 1.0, 0.04);
   bf.Add_par("D2m", 1.0, 0.05);
   bf.Add_par("D4", 1.0, 0.05);
   bf.Add_par("Al", 10.0, 0.1);


   //Fix parameters depending on the fit type
   bf.Fix_par("A2m", 0.0);
   bf.Fix_par("D2m", 0.0);
   bf.Fix_par("D4", 0.0);







   //set ansatz
   bf.ansatz= [=](const sc_fpar& X, const sc_ipar& Y) {

		double ret_val=0;
		double gs_ph = guess_w0_X; 
		double MpL= Y.Mp*Y.L;

		if(Use_w0_Xpi) {

		  ret_val = X.w0X*gs_ph*pow(Y.csi/Y.csi_ph, 2.0/5.0)*pow(( 1.0  -2.0*Y.csi*log(Y.csi/Y.csi_ph)  +  X.Am*(Y.csi-Y.csi_ph) + X.A2m*( pow(Y.csi,2) - pow(Y.csi_ph,2)) + X.D2*pow(Y.w0, -2.0) + X.D2m*Y.csi*pow(Y.w0, -2.0) + X.D4*pow(Y.w0, -4.0) ), 1.0/5.0)*(1.0 + X.Al*pow(Y.csi,2)*exp(-MpL)/pow(MpL, 3.0/2.0)); 

		}

		else {

		  ret_val = X.w0X*gs_ph*( 1.0  -2.0*Y.csi*log(Y.csi/Y.csi_ph) +  X.Am*(Y.csi-Y.csi_ph) + X.A2m*( pow(Y.csi,2) - pow(Y.csi_ph,2)) + X.D2*pow(Y.w0, -2.0) + X.D2m*Y.csi*pow(Y.w0, -2.0) + X.D4*pow(Y.w0, -4.0) )*(1.0 + X.Al*Y.csi*exp(-MpL)/pow(MpL, 3.0/2.0)); 

		}

		return ret_val;
		
			     
	      };

   
   bf.measurement = [=](const sc_fpar& X, const sc_ipar& Y) {
		      double ret_val;
		      if(Use_w0_Xpi) ret_val= Y.w0_t_Xpi;
		      else ret_val = Y.w0_t_fp;
		      return ret_val;
		    };
   bf.error = [=](const sc_fpar& X, const sc_ipar& Y) {
		double ret_val;
		if(Use_w0_Xpi) ret_val= Y.w0_t_Xpi_err;
		else ret_val = Y.w0_t_fp_err;
		return ret_val;
	      };



   //populate bf
   
    vector<vector<sc_ipar>> ipar_all_ens(Njacks);
    for(auto &ipar_jack: ipar_all_ens) ipar_jack.resize(Nens);
    
    for(int iens=0; iens<Nens;iens++) {  
      for(int ijack=0;ijack<Njacks;ijack++) {
	

	
	ipar_all_ens[ijack][iens].Mp = Mpi_CDH.distr_list[iens].distr[ijack];
	ipar_all_ens[ijack][iens].Mp_err = Mpi_CDH.err(iens);
	
	
	ipar_all_ens[ijack][iens].fp = fp_CDH.distr_list[iens].distr[ijack]; 
	ipar_all_ens[ijack][iens].fp_err = fp_CDH.err(iens);
	
	
	ipar_all_ens[ijack][iens].Xpi = Xpi.distr_list[iens].distr[ijack]; 
	ipar_all_ens[ijack][iens].Xpi_err = Xpi.err(iens);
	
	
	ipar_all_ens[ijack][iens].csi= Csi_CDH.distr_list[iens].distr[ijack]; 
	ipar_all_ens[ijack][iens].csi_err = Csi_CDH.err(iens);
	
			
	ipar_all_ens[ijack][iens].w0 = w0.distr_list[iens].distr[ijack]; 
	ipar_all_ens[ijack][iens].w0_err = w0.err(iens);
	
	
	ipar_all_ens[ijack][iens].w0_t_fp = fp_CDH.distr_list[iens].distr[ijack]*w0.distr_list[iens].distr[ijack]; 
	ipar_all_ens[ijack][iens].w0_t_fp_err = (fp_CDH*w0).err(iens);
	
	
	ipar_all_ens[ijack][iens].w0_t_Xpi = Xpi.distr_list[iens].distr[ijack]*w0.distr_list[iens].distr[ijack]; 
	ipar_all_ens[ijack][iens].w0_t_Xpi_err = (Xpi*w0).err(iens);
	
	
	ipar_all_ens[ijack][iens].L = L_list[iens]; 
	ipar_all_ens[ijack][iens].csi_ph = csi_phys_distr.distr[ijack];
	
      }
      
    }
    
    
    bf.Append_to_input_par(ipar_all_ens);
    //fit
    boot_fit_data<sc_fpar> Bt_sc_fit = bf.Perform_bootstrap_fit();


    //retrieve parameters
    distr_t w0X, Am, A2m, D2, D2m, D4,  Al;  //constructor sets UseJack to 1

  
 
    for(int ijack=0;ijack<Njacks;ijack++) {

      sc_fpar my_w0_fit_pars = Bt_sc_fit.par[ijack];
      w0X.distr.push_back(my_w0_fit_pars.w0X);
      Am.distr.push_back(my_w0_fit_pars.Am);
      A2m.distr.push_back(my_w0_fit_pars.A2m);
      D2.distr.push_back(my_w0_fit_pars.D2);
      D2m.distr.push_back(my_w0_fit_pars.D2m);
      D4.distr.push_back(my_w0_fit_pars.D4);
      Al.distr.push_back(my_w0_fit_pars.Al);
    }

    //get average ch2
    double ch2  = Bt_sc_fit.get_ch2_ave();



    //cout infos

    cout<<"########  scale setting fit parameter : "<<" ##########"<<endl;
    cout<<"Use w0_Xp: "<<Use_w0_Xpi<<endl;
    cout<<"ch2 total: "<<ch2<<endl;
    cout<<"w0X: "<<w0X.ave()*guess_w0_X<<" +- "<<w0X.err()*guess_w0_X<<endl;
    cout<<"w0X: [relative to arXiv 2104:06747] "<<w0X.ave()<<" +- "<<w0X.err()<<endl;
    cout<<"Am: "<<Am.ave()<<" +- "<<Am.err()<<endl;
    cout<<"A2m: "<<A2m.ave()<<" +- "<<A2m.err()<<endl;
    cout<<"D2: "<<D2.ave()<<" +- "<<D2.err()<<endl;
    cout<<"D2m: "<<D2m.ave()<<" +- "<<D2m.err()<<endl;
    cout<<"D4: "<<D4.ave()<<" +- "<<D4.err()<<endl;
    cout<<"Al: "<<Al.ave()<<" +- "<<Al.err()<<endl;





    //print fitting function


    //start with continuum limit as a function of csi
    int num_csi_points= 100;
    Vfloat csi_p_list;
    int num_lat = (Use_three_finest)?4:5;
    VVfloat extr_val(num_lat), extr_err(num_lat);

    for(int ilat=0;ilat<num_lat;ilat++) {

      
      distr_t w0_p(UseJack);
      if(ilat==1) { // Ensemble D
	bool Ens_found=false;
	for(int iens=0;Nens_total;iens++) {
	  if(Ensemble_tags[iens] == "cD211a.054.96") {Ens_found=true; w0_p = w0.distr_list[iens];}}
	
	if(!Ens_found) crash("Cannot find Ensemble cD211a.054.96 in Ensemble_tags");
      }
      else if(ilat==2) { // Ensemble C
	bool Ens_found=false;
	for(int iens=0;Nens_total;iens++) {
	  if(Ensemble_tags[iens] == "cC211a.06.80") {Ens_found=true; w0_p = w0.distr_list[iens];}}
	
	if(!Ens_found) crash("Cannot find Ensemble cC211a.06.80 in Ensemble_tags");
      }
      else if(ilat==3) { //Ensemble B
	bool Ens_found=false;
	for(int iens=0;Nens_total;iens++) {
	  if(Ensemble_tags[iens] == "cB211b.072.96") {Ens_found=true; w0_p = w0.distr_list[iens];}}
	
	if(!Ens_found) crash("Cannot find Ensemble cB211b.072.96 in Ensemble_tags");

      }
      else if(ilat==4) { //Ensemble A
	bool Ens_found=false;
	for(int iens=0;Nens_total;iens++) {
	  if(Ensemble_tags[iens] == "cA211ab.30.32") {Ens_found=true; w0_p = w0.distr_list[iens];}}
	
	if(!Ens_found) crash("Cannot find Ensemble cA211ab.30.32 in Ensemble_tags");

      }
      else crash("ilat has an invalid value. ilat: "+to_string(ilat));
      
      
      for(int icsi=0; icsi<num_csi_points;icsi++) {

	double csi_p = icsi*10.0*(csi_phys)/num_csi_points;
	csi_p_list.push_back(csi_p);

	distr_t extr(UseJack);
	
	for(int ijack=0; ijack<Njacks;ijack++) {

	  sc_fpar Xi;
	  sc_ipar Yi;

	  //populate fit parameters
	  Xi.w0X = w0X.distr[ijack];
	  Xi.Am = Am.distr[ijack];
	  Xi.A2m = A2m.distr[ijack];
	  Xi.D2 = (ilat==0)?0.0:D2.distr[ijack];
	  Xi.D2m = (ilat==0)?0.0:D2m.distr[ijack];
	  Xi.D4 = (ilat==0)?0.0:D4.distr[ijack];
	  Xi.Al = 0.0; //we look at thermodynamic limit here

	  //populate input parameters
	  Yi.L = 10.0; //will not be used
	  Yi.csi = csi_p;
	  Yi.csi_ph = csi_phys_distr.distr[ijack];
	  Yi.w0 = w0_p.distr[ijack];
	  



	  
	  extr.distr.push_back( bf.ansatz(Xi,Yi));
	


	}
	
	extr_val[ilat].push_back( extr.ave());
	extr_err[ilat].push_back( extr.err());
      }

    }


    //print extr_val and error for each lattice spacing
    string tag_w0X = (Use_w0_Xpi)?"on":"off";

    Print_To_File({}, {csi_p_list, extr_val[0], extr_err[0]}, "../data/scale_setting/fit_cont_w0Xpi_"+tag_w0X+".dat", "", "# csi val err");
    Print_To_File({}, {csi_p_list, extr_val[1], extr_err[1]}, "../data/scale_setting/fit_D_w0Xpi_"+tag_w0X+".dat", "", "# csi val err");
    Print_To_File({}, {csi_p_list, extr_val[2], extr_err[2]}, "../data/scale_setting/fit_C_w0Xpi_"+tag_w0X+".dat", "", "# csi val err");
    Print_To_File({}, {csi_p_list, extr_val[3], extr_err[3]}, "../data/scale_setting/fit_B_w0Xpi_"+tag_w0X+".dat", "", "# csi val err");
    if(!Use_three_finest) Print_To_File({}, {csi_p_list, extr_val[4], extr_err[4]}, "../data/scale_setting/fit_A_w0Xpi_"+tag_w0X+".dat", "", "# csi val err");

    




    //correct lattice points by FSEs
    distr_t_list w0X_corr_by_FSEs(UseJack);
    distr_t_list w0X_uncorr(UseJack);
    vector<string> Corrected_ensemble_list;
    for(int iens=0;iens<Nens;iens++) {
      Corrected_ensemble_list.push_back(Ensemble_tags[iens]);
      distr_t MpL_distr = L_list[iens]*Mpi_CDH.distr_list[iens];
      if(Use_w0_Xpi) {
	w0X_corr_by_FSEs.distr_list.push_back( w0.distr_list[iens]*Xpi.distr_list[iens]/(1.0 + Al*distr_t::f_of_distr(pow_2,Csi_CDH.distr_list[iens])*distr_t::f_of_distr(FVE, MpL_distr)));
	w0X_uncorr.distr_list.push_back( w0.distr_list[iens]*Xpi.distr_list[iens]); 
      }
      else {
	w0X_corr_by_FSEs.distr_list.push_back( w0.distr_list[iens]*fp_CDH.distr_list[iens]/(1.0 + Al*Csi_CDH.distr_list[iens]*distr_t::f_of_distr(FVE, MpL_distr)));
	w0X_uncorr.distr_list.push_back( w0.distr_list[iens]*fp_CDH.distr_list[iens]); 
      }
    }

    //Print corrected lattice points to file
    Print_To_File({Corrected_ensemble_list},{w0X_corr_by_FSEs.ave(), w0X_corr_by_FSEs.err(), w0X_uncorr.ave(), w0X_uncorr.err()},  "../data/scale_setting/data_FSEs_corrected_w0Xpi_"+tag_w0X+".dat", "", "# Ens   corrected   uncorrected");


    


    //get lattice spacings

    distr_t w0_phys_units= w0X*guess_w0_X/( (Use_w0_Xpi)?Xpi_phys_distr:fpi_phys_distr);



    //get w0/a from cD211a.054.96
    distr_t w0a_D;
    bool Ens_found_D=false;
    for(int iens=0;Nens_total;iens++) {
      if(Ensemble_tags[iens] == "cD211a.054.96") {Ens_found_D=true; w0a_D = w0.distr_list[iens];}}	
    if(!Ens_found_D) crash("Cannot find Ensemble cD211a.054.96 in Ensemble_tags");


    //get w0/a from cC211a.06.80
    distr_t w0a_C;
    bool Ens_found_C=false;
    for(int iens=0;Nens_total;iens++) {
      if(Ensemble_tags[iens] == "cC211a.06.80") {Ens_found_C=true; w0a_C = w0.distr_list[iens];}}	
    if(!Ens_found_C) crash("Cannot find Ensemble cC211a.06.80 in Ensemble_tags");


    //get w0/a from cB211b.072.96
    distr_t w0a_B;
    bool Ens_found_B=false;
    for(int iens=0;Nens_total;iens++) {
      if(Ensemble_tags[iens] == "cB211b.072.96") {Ens_found_B=true; w0a_B = w0.distr_list[iens];}}	
    if(!Ens_found_B) crash("Cannot find Ensemble cB211b.072.96 in Ensemble_tags");

    //get w0/a from cA211ab.30.32
    distr_t w0a_A;
    bool Ens_found_A=false;
    for(int iens=0;Nens_total;iens++) {
      if(Ensemble_tags[iens] == "cA211ab.30.32") {Ens_found_A=true; w0a_A = w0.distr_list[iens];}}	
    if(!Ens_found_A) crash("Cannot find Ensemble cA211ab.30.32 in Ensemble_tags");
    
    

    a_distr_A = w0_phys_units/w0a_A  ;
    
    a_distr_B = w0_phys_units/w0a_B  ;

    a_distr_C = w0_phys_units/w0a_C  ;

    a_distr_D = w0_phys_units/w0a_D  ;



    //print on screen
    cout<<"#######  Printing lattice spacing determined from w0X analysis ########"<<endl;
    cout<<"A ensemble: "<< a_distr_A.ave()<<" +- "<<a_distr_A.err()<<" [GeV-1]    "<<(a_distr_A/fm_to_iGev).ave()<<" +- "<<(a_distr_A/fm_to_iGev).err()<<" [fm] "<<endl;
    cout<<"B ensemble: "<< a_distr_B.ave()<<" +- "<<a_distr_B.err()<<" [GeV-1]    "<<(a_distr_B/fm_to_iGev).ave()<<" +- "<<(a_distr_B/fm_to_iGev).err()<<" [fm] "<<endl;
    cout<<"C ensemble: "<< a_distr_C.ave()<<" +- "<<a_distr_C.err()<<" [GeV-1]    "<<(a_distr_C/fm_to_iGev).ave()<<" +- "<<(a_distr_C/fm_to_iGev).err()<<" [fm] "<<endl;
    cout<<"D ensemble: "<< a_distr_D.ave()<<" +- "<<a_distr_D.err()<<" [GeV-1]    "<<(a_distr_D/fm_to_iGev).ave()<<" +- "<<(a_distr_D/fm_to_iGev).err()<<" [fm] "<<endl;


  


    //return 


  return;
}







void Determine_scale_from_fp(const distr_t_list &Mpi_phys_point, const distr_t_list &fpi_phys_point, const distr_t_list &Mpi_Bens, const distr_t_list &fpi_Bens, const distr_t_list &Mpi_Aens, const distr_t_list &fpi_Aens,  const vector<double> &L_phys_point, const vector<double> &L_B_ens, const vector<double> &L_A_ens, const vector<string>  &Ensemble_phys_point_tag_list, const vector<string>  &Ensemble_B_tag_list, const vector<string>  &Ensemble_A_tag_list, distr_t &a_distr_A, distr_t &a_distr_B, distr_t &a_distr_C, distr_t &a_distr_D,bool UseJack,bool Use_three_finest_in_scale_setting_fp) {


  cout<<"Starting scale setting analysis using afp"<<endl;

  
  //init Gaussian number generator
  GaussianMersenne GM(76795699);


  int Nens_phys = L_phys_point.size();
  int Nens_B = L_B_ens.size();
  int Nens_A = L_A_ens.size();
  if(Nens_phys==0 || Nens_B==0 || Nens_A==0) crash("in Determine_scale_setting_from_fp Nens_phys=0 || Nens_B = 0 || Nens_A = 0");
  int Njacks= Mpi_phys_point.distr_list[0].size();


  //generate K_\ell correcting factor for cA211a.12.48 ensemble according to 2104.06747
  distr_t Kl(UseJack);
  for(int ijack=0;ijack<Njacks;ijack++) { Kl.distr.push_back( sqrt( 1.0 + pow( corr_fact + GM()*corr_fact_err/sqrt(Njacks -1.0)    ,2)));} 
  

  
  //List of lambda functions needed

  auto FVE = [](double x) -> double {return (1.0/pow(x,1.5))*exp(-x);};
  auto LOG = [](double x) { return log(x);};
  auto pow_1_5 = [](double x) -> double { return pow(x, 1.0/5.0);};
  auto pow_4 = [](double x) -> double { return pow(x, 4.0);};
  auto pow_2 = [](double x) -> double { return pow(x, 2.0);};
  auto pow_3_2 = [](double x) -> double { return pow(x, 3.0/2.0);};
  auto EXP = [](double x) -> double {return exp(x);};



  //###########################



  //Output directories

  boost::filesystem::create_directory("../data/scale_setting");
  
  //##########################


  //Generate fake distribution for csi_phys, Mpi , fp, Xpi
 
  distr_t Mpi_phys_distr(UseJack);
  distr_t fpi_phys_distr(UseJack);
  distr_t csi_phys_distr(UseJack);
  distr_t Xpi_phys_distr(UseJack);
  
  for(int ij=0;ij<Njacks;ij++) {
    double Mfak = Mp_phys + GM()*Mp_phys_err/sqrt(Njacks -1.0);
    double Ffak = fp_phys + GM()*fp_phys_err/sqrt(Njacks -1.0);

    Mpi_phys_distr.distr.push_back( Mfak);
    fpi_phys_distr.distr.push_back(Ffak);
    csi_phys_distr.distr.push_back( pow(Mfak/(4.0*M_PI*Ffak),2.0));
    Xpi_phys_distr.distr.push_back( pow( Ffak*pow(Mfak,4.0), 1.0/5.0));
  }

 
  distr_t_list fp_CDH_phys_point(UseJack), Mpi_CDH_phys_point(UseJack), Xpi_phys_point(UseJack), Csi_CDH_phys_point(UseJack);
  distr_t_list fp_CDH_B(UseJack), Mpi_CDH_B(UseJack), Xpi_B(UseJack), Csi_CDH_B(UseJack);
  distr_t_list fp_CDH_A(UseJack), Mpi_CDH_A(UseJack), Xpi_A(UseJack), Csi_CDH_A(UseJack);

   
  //compute CDH corrections to Mpi and fpi
  for(int iens=0; iens<Nens_phys;iens++) {
    //GL and CDH correction for Mpi and fp
    distr_t csi_L = Mpi_phys_point.distr_list[iens]*Mpi_phys_point.distr_list[iens]/(pow(4.0*M_PI,2)*fpi_phys_point.distr_list[iens]*fpi_phys_point.distr_list[iens]);
    distr_t g1 = distr_t::f_of_distr(g1_l, Mpi_phys_point.distr_list[iens]*L_phys_point[iens]);
    distr_t g2 = distr_t::f_of_distr(g2_l, Mpi_phys_point.distr_list[iens]*L_phys_point[iens]);
    distr_t log_l = log(csi_phys) - distr_t::f_of_distr(LOG, csi_L);
    fp_CDH_phys_point.distr_list.push_back(fpi_phys_point.distr_list[iens]/(1.0 -2.0*csi_L*g1 +2.0*csi_L*csi_L*( (Cf1(l1ph,l2ph,l3ph,l4ph) + Sf1(s0,s1,s2,s3) + Cf1_log()*log_l)*g1 + (Cf2(l1ph,l2ph,l3ph,l4ph) + Sf2(s0,s1,s2,s3) + Cf2_log()*log_l)*g2)));
    Mpi_CDH_phys_point.distr_list.push_back(Mpi_phys_point.distr_list[iens]/(1.0 + 0.5*csi_L*g1 - csi_L*csi_L*( (Cm1(l1ph,l2ph,l3ph,l4ph) + Sm1(s0,s1,s2,s3) + Cm1_log()*log_l)*g1 + (Cm2(l1ph,l2ph,l3ph,l4ph) + Sm2(s0,s1,s2,s3) + Cm2_log()*log_l)*g2)));
    Csi_CDH_phys_point.distr_list.push_back( Mpi_CDH_phys_point.distr_list[iens]*Mpi_CDH_phys_point.distr_list[iens]/(16.0*M_PI*M_PI*fp_CDH_phys_point.distr_list[iens]*fp_CDH_phys_point.distr_list[iens]));
    Xpi_phys_point.distr_list.push_back( distr_t::f_of_distr(pow_1_5, fp_CDH_phys_point.distr_list[iens]*distr_t::f_of_distr(pow_4, Mpi_CDH_phys_point.distr_list[iens])));
  }


  int cB211a2548_tag=0;
   
  for(int iens=0; iens<Nens_B;iens++) {

    //GL and CDH correction for Mpi and fp
    distr_t csi_L = Mpi_Bens.distr_list[iens]*Mpi_Bens.distr_list[iens]/(pow(4.0*M_PI,2)*fpi_Bens.distr_list[iens]*fpi_Bens.distr_list[iens]);
    distr_t g1 = distr_t::f_of_distr(g1_l, Mpi_Bens.distr_list[iens]*L_B_ens[iens]);
    distr_t g2 = distr_t::f_of_distr(g2_l, Mpi_Bens.distr_list[iens]*L_B_ens[iens]);
    distr_t log_l = log(csi_phys) - distr_t::f_of_distr(LOG, csi_L);
    fp_CDH_B.distr_list.push_back(fpi_Bens.distr_list[iens]/(1.0 -2.0*csi_L*g1 +2.0*csi_L*csi_L*( (Cf1(l1ph,l2ph,l3ph,l4ph) + Sf1(s0,s1,s2,s3) + Cf1_log()*log_l)*g1 + (Cf2(l1ph,l2ph,l3ph,l4ph) + Sf2(s0,s1,s2,s3) + Cf2_log()*log_l)*g2)));
    Mpi_CDH_B.distr_list.push_back(Mpi_Bens.distr_list[iens]/(1.0 + 0.5*csi_L*g1 - csi_L*csi_L*( (Cm1(l1ph,l2ph,l3ph,l4ph) + Sm1(s0,s1,s2,s3) + Cm1_log()*log_l)*g1 + (Cm2(l1ph,l2ph,l3ph,l4ph) + Sm2(s0,s1,s2,s3) + Cm2_log()*log_l)*g2)));
    Csi_CDH_B.distr_list.push_back( Mpi_CDH_B.distr_list[iens]*Mpi_CDH_B.distr_list[iens]/(16.0*M_PI*M_PI*fp_CDH_B.distr_list[iens]*fp_CDH_B.distr_list[iens]));
    Xpi_B.distr_list.push_back( distr_t::f_of_distr(pow_1_5, fp_CDH_B.distr_list[iens]*distr_t::f_of_distr(pow_4, Mpi_CDH_B.distr_list[iens])));
    if(Ensemble_B_tag_list[iens] == "cB211a.25.48") cB211a2548_tag = iens;

  }

  for(int iens=0; iens< Nens_B;iens++) {
    if( Ensemble_B_tag_list[iens] == "cB211a.25.24" || Ensemble_B_tag_list[iens] == "cB211a.25.32") Csi_CDH_B.distr_list[iens] = Csi_CDH_B.distr_list[cB211a2548_tag];
  }

 
  for(int iens=0; iens<Nens_A;iens++) {

    //GL and CDH correction for Mpi and fp
    distr_t csi_L = Mpi_Aens.distr_list[iens]*Mpi_Aens.distr_list[iens]/(pow(4.0*M_PI,2)*fpi_Aens.distr_list[iens]*fpi_Aens.distr_list[iens]);
    distr_t g1 = distr_t::f_of_distr(g1_l, Mpi_Aens.distr_list[iens]*L_A_ens[iens]);
    distr_t g2 = distr_t::f_of_distr(g2_l, Mpi_Aens.distr_list[iens]*L_A_ens[iens]);
    distr_t log_l = log(csi_phys) - distr_t::f_of_distr(LOG, csi_L);
    fp_CDH_A.distr_list.push_back(fpi_Aens.distr_list[iens]/(1.0 -2.0*csi_L*g1 +2.0*csi_L*csi_L*( (Cf1(l1ph,l2ph,l3ph,l4ph) + Sf1(s0,s1,s2,s3) + Cf1_log()*log_l)*g1 + (Cf2(l1ph,l2ph,l3ph,l4ph) + Sf2(s0,s1,s2,s3) + Cf2_log()*log_l)*g2)));
    //if CA211a.12.48 correct fpi for violation of maximal twist condition
    if(Ensemble_A_tag_list[iens] == "cA211a.12.48") {
      fp_CDH_A.distr_list[iens] = fp_CDH_A.distr_list[iens]*Kl;
    }
    Mpi_CDH_A.distr_list.push_back(Mpi_Aens.distr_list[iens]/(1.0 + 0.5*csi_L*g1 - csi_L*csi_L*( (Cm1(l1ph,l2ph,l3ph,l4ph) + Sm1(s0,s1,s2,s3) + Cm1_log()*log_l)*g1 + (Cm2(l1ph,l2ph,l3ph,l4ph) + Sm2(s0,s1,s2,s3) + Cm2_log()*log_l)*g2)));
    Csi_CDH_A.distr_list.push_back( Mpi_CDH_A.distr_list[iens]*Mpi_CDH_A.distr_list[iens]/(16.0*M_PI*M_PI*fp_CDH_A.distr_list[iens]*fp_CDH_A.distr_list[iens]));
    Xpi_A.distr_list.push_back( distr_t::f_of_distr(pow_1_5, fp_CDH_A.distr_list[iens]*distr_t::f_of_distr(pow_4, Mpi_CDH_A.distr_list[iens])));

  }


  //Print_To_File
  Print_To_File({Ensemble_phys_point_tag_list}, {L_phys_point, Mpi_CDH_phys_point.ave(), Mpi_CDH_phys_point.err(), fp_CDH_phys_point.ave(), fp_CDH_phys_point.err(), Xpi_phys_point.ave(), Xpi_phys_point.err(), Csi_CDH_phys_point.ave(), Csi_CDH_phys_point.err()}, "../data/scale_setting/pion_obs_phys_point.dat", "", "#Ens L Mpi fp Xpi  Csi");

  Print_To_File({Ensemble_B_tag_list}, {L_B_ens, Mpi_CDH_B.ave(), Mpi_CDH_B.err(), fp_CDH_B.ave(), fp_CDH_B.err(), Xpi_B.ave(), Xpi_B.err(), Csi_CDH_B.ave(), Csi_CDH_B.err()}, "../data/scale_setting/pion_obs_B_ens.dat", "", "#Ens L Mpi fp Xpi  Csi");

  Print_To_File({Ensemble_A_tag_list}, {L_A_ens, Mpi_CDH_A.ave(), Mpi_CDH_A.err(), fp_CDH_A.ave(), fp_CDH_A.err(), Xpi_A.ave(), Xpi_A.err(), Csi_CDH_A.ave(), Csi_CDH_A.err()}, "../data/scale_setting/pion_obs_A_ens.dat", "", "#Ens L Mpi fp Xpi  Csi");

  


  //Perform extrapolation of afp

   bootstrap_fit<sc_fp_fpar,sc_fp_ipar> bf(Njacks);
   bf.set_warmup_lev(1);
   bf.Set_number_of_measurements(Nens_B+Nens_A);
   bf.Set_verbosity(scale_setting_fit_verbosity);
   bf.Set_print_path("chi2_scale_setting_fp");

   //B ens pars
   bf.Add_par("afp_B", 0.053, 0.001);
   bf.Add_par("Am_B", 10.0, 0.1);
   bf.Add_par("A2m_B", 1.0, 0.1);
   bf.Add_par("Al_B", 10.0, 0.1);

   //A ens pars
   bf.Add_par("afp_A", 0.060, 0.001);
   bf.Add_par("Am_A", 11.0, 0.1);
   bf.Add_par("A2m_A", 1.0, 0.1);
   bf.Add_par("Al_A", 10.0, 0.1);

   bf.Add_par("Al", 10.0, 0.1);

  

   //Fix parameters depending on the fit type
   if(!Enable_A2m_cont) bf.Fix_par("A2m_B", 0.0);
   bf.Fix_par("A2m_A", 0.0);


   if(Use_three_finest_in_scale_setting_fp) {
     //fix parameters for A ensembles
     bf.Fix_par("Al",0.0);
   }
   else {
     bf.Fix_par("Al_B",0.0);
     bf.Fix_par("Al_A",0.0);

   }



   
  
   //set ansatz
   bf.ansatz= [=](const sc_fp_fpar& X, const sc_fp_ipar& Y) {

		double ret_val=0;
		double MpL= Y.Mp*Y.L;

		if(Y.Is_B) {

		  if(Use_fpi_or_Xpi) ret_val = X.afp_B*(1.0  -2.0*Y.csi*log(Y.csi/Y.csi_ph) + X.Am_B*(Y.csi - Y.csi_ph) + X.A2m_B*(pow(Y.csi,2) - pow(Y.csi_ph,2)))*(1.0 + (X.Al+ X.Al_B)*Y.csi*exp(-MpL)/pow(MpL, 3.0/2.0));
		  else crash("Xpi not yet implemented"); 

		}

		else {

		  if(Use_fpi_or_Xpi) ret_val =  X.afp_A*(1.0  -2.0*Y.csi*log(Y.csi/Y.csi_ph) + X.Am_A*(Y.csi - Y.csi_ph) + X.A2m_B*(pow(Y.csi,2) - pow(Y.csi_ph,2)))*(1.0 + (X.Al+ X.Al_A)*Y.csi*exp(-MpL)/pow(MpL, 3.0/2.0));
		  else crash("Xpi not yet implemented") ;

		}

		return ret_val;
		
			     
	      };

   
   bf.measurement = [=](const sc_fp_fpar& X, const sc_fp_ipar& Y) {
		      double ret_val= Use_fpi_or_Xpi?Y.fp:Y.Xpi;
		      return ret_val;
		    };
   bf.error = [=](const sc_fp_fpar& X, const sc_fp_ipar& Y) {
		double ret_val;
		ret_val = Use_fpi_or_Xpi?Y.fp_err:Y.Xpi_err;
		return ret_val;
	      };



   //populate bf

   int Nens_tot = (Nens_B + Nens_A);
    vector<vector<sc_fp_ipar>> ipar_all_ens(Njacks);
    for(auto &ipar_jack: ipar_all_ens) ipar_jack.resize(Nens_tot);
    
    for(int iens=0; iens<Nens_tot;iens++) {  
      for(int ijack=0;ijack<Njacks;ijack++) {
	
	if(iens < Nens_B) {

	  ipar_all_ens[ijack][iens].Mp = Mpi_CDH_B.distr_list[iens].distr[ijack]; 
	  ipar_all_ens[ijack][iens].Mp_err = Mpi_CDH_B.err(iens);
	  
	
	  ipar_all_ens[ijack][iens].fp = fp_CDH_B.distr_list[iens].distr[ijack]; 
	  ipar_all_ens[ijack][iens].fp_err = fp_CDH_B.err(iens);
	
	
	  ipar_all_ens[ijack][iens].Xpi = Xpi_B.distr_list[iens].distr[ijack]; 
	  ipar_all_ens[ijack][iens].Xpi_err = Xpi_B.err(iens);
	
	
	  ipar_all_ens[ijack][iens].csi= Csi_CDH_B.distr_list[iens].distr[ijack]; 
	  ipar_all_ens[ijack][iens].csi_err = Csi_CDH_B.err(iens);
	
			
	  ipar_all_ens[ijack][iens].L = L_B_ens[iens]; 
	  ipar_all_ens[ijack][iens].csi_ph = csi_phys_distr.distr[ijack];

	  ipar_all_ens[ijack][iens].Is_B = true;

	}

	else {

	  ipar_all_ens[ijack][iens].Mp = Mpi_CDH_A.distr_list[iens-Nens_B].distr[ijack]; 
	  ipar_all_ens[ijack][iens].Mp_err = Mpi_CDH_A.err(iens-Nens_B);

	  ipar_all_ens[ijack][iens].fp = fp_CDH_A.distr_list[iens-Nens_B].distr[ijack]; 
	  ipar_all_ens[ijack][iens].fp_err = fp_CDH_A.err(iens-Nens_B);
	
	
	  ipar_all_ens[ijack][iens].Xpi = Xpi_A.distr_list[iens-Nens_B].distr[ijack]; 
	  ipar_all_ens[ijack][iens].Xpi_err = Xpi_A.err(iens-Nens_B);
	
	
	  ipar_all_ens[ijack][iens].csi= Csi_CDH_A.distr_list[iens-Nens_B].distr[ijack]; 
	  ipar_all_ens[ijack][iens].csi_err = Csi_CDH_A.err(iens-Nens_B);
	
			
	  ipar_all_ens[ijack][iens].L = L_A_ens[iens-Nens_B]; 
	  ipar_all_ens[ijack][iens].csi_ph = csi_phys_distr.distr[ijack];

	  ipar_all_ens[ijack][iens].Is_B = false;

	}
	
      }
      
    }

      
    bf.Append_to_input_par(ipar_all_ens);
    //fit
    boot_fit_data<sc_fp_fpar> Bt_sc_fit = bf.Perform_bootstrap_fit();


    //retrieve parameters
    distr_t afp_B, Am_B, A2m_B, Al_B;  //constructor sets UseJack to 1
    distr_t afp_A, Am_A, A2m_A, Al_A;
    distr_t Al;

  
 
    for(int ijack=0;ijack<Njacks;ijack++) {

      sc_fp_fpar my_w0_fit_pars = Bt_sc_fit.par[ijack];
      afp_B.distr.push_back(my_w0_fit_pars.afp_B);
      Am_B.distr.push_back(my_w0_fit_pars.Am_B);
      A2m_B.distr.push_back(my_w0_fit_pars.A2m_B);
      Al_B.distr.push_back(my_w0_fit_pars.Al_B);
      afp_A.distr.push_back(my_w0_fit_pars.afp_A);
      Am_A.distr.push_back(my_w0_fit_pars.Am_A);
      A2m_A.distr.push_back(my_w0_fit_pars.A2m_A);
      Al_A.distr.push_back(my_w0_fit_pars.Al_A);
      Al.distr.push_back(my_w0_fit_pars.Al);
      
    }

    //get average ch2
    double ch2  = Bt_sc_fit.get_ch2_ave();



    //cout infos

    cout<<"########  scale setting fit parameter : "<<" ##########"<<endl;
    cout<<"Use 3 finest?: "<<Use_three_finest_in_scale_setting_fp<<endl;
    cout<<"ch2 total: "<<ch2<<endl;
    cout<<"afp_B: "<<afp_B.ave()<<" +- "<<afp_B.err()<<endl;
    cout<<"Am_B: "<<Am_B.ave()<<" +- "<<Am_B.err()<<endl;
    cout<<"A2m_B: "<<A2m_B.ave()<<" +- "<<A2m_B.err()<<endl;
    cout<<"Al_B: "<<Al_B.ave()<<" +- "<<Al_B.err()<<endl;
    cout<<"afp_A: "<<afp_A.ave()<<" +- "<<afp_A.err()<<endl;
    cout<<"Am_A: "<<Am_A.ave()<<" +- "<<Am_A.err()<<endl;
    cout<<"A2m_A: "<<A2m_A.ave()<<" +- "<<A2m_A.err()<<endl;
    cout<<"Al_A: "<<Al_A.ave()<<" +- "<<Al_A.err()<<endl;
    cout<<"Al: "<<Al.ave()<<" +- "<<Al.err()<<endl;





    //print fitting function


    //start with continuum limit as a function of csi
    int num_csi_points= 100;
    Vfloat csi_p_list;
    Vfloat extr_A_val, extr_A_err;
    Vfloat extr_B_val, extr_B_err;

    
      
      
    for(int icsi=0; icsi<num_csi_points;icsi++) {

      double csi_p = icsi*10.0*(csi_phys)/num_csi_points;
      csi_p_list.push_back(csi_p);
      
      distr_t extr_A(UseJack);
      distr_t extr_B(UseJack);
      
      for(int ijack=0; ijack<Njacks;ijack++) {
	
	sc_fp_fpar Xi;
	sc_fp_ipar Yi_A;
	sc_fp_ipar Yi_B;

	//populate fit parameters
	Xi.afp_B = afp_B.distr[ijack];
	Xi.Am_B = Am_B.distr[ijack];
	Xi.A2m_B = A2m_B.distr[ijack];
	Xi.afp_A = afp_A.distr[ijack];
	Xi.Am_A = Am_A.distr[ijack];
	Xi.A2m_A = A2m_A.distr[ijack];
	Xi.Al = 0.0; //we look at thermodynamic limit here
	Xi.Al_B = 0.0;
	Xi.Al_A = 0.0;
	
	//populate input parameters
	Yi_B.L = 10.0; //will not be used
	Yi_B.csi = csi_p;
	Yi_B.csi_ph = csi_phys_distr.distr[ijack];
	Yi_B.Is_B = true;
	Yi_B.Mp = 1.0;

	Yi_A.L = 10.0; //will not be used
	Yi_A.csi = csi_p;
	Yi_A.csi_ph = csi_phys_distr.distr[ijack];
	Yi_A.Mp = 1.0;
	Yi_A.Is_B = false;
	  
	extr_B.distr.push_back( bf.ansatz(Xi,Yi_B));
        extr_A.distr.push_back( bf.ansatz(Xi, Yi_A));
	


      }
	
      extr_B_val.push_back( extr_B.ave());
      extr_B_err.push_back( extr_B.err());
      extr_A_val.push_back( extr_A.ave());
      extr_A_err.push_back( extr_A.err());
      
    }

   

    //print extr_val and error for each lattice spacing
    

    Print_To_File({}, {csi_p_list, extr_B_val, extr_B_err}, "../data/scale_setting/fit_B_afp.dat", "", "# csi val err");
    Print_To_File({}, {csi_p_list, extr_A_val, extr_A_err}, "../data/scale_setting/fit_A_afp.dat", "", "# csi val err");
    
    


    //correct lattice points by FSEs
    distr_t_list afp_corr_by_FSEs_B(UseJack);
    distr_t_list afp_uncorr_B(UseJack);
    distr_t_list afp_corr_by_FSEs_A(UseJack);
    distr_t_list afp_uncorr_A(UseJack);
    
    for(int iens=0;iens<Nens_B;iens++) {
      distr_t MpL_distr = L_B_ens[iens]*Mpi_CDH_B.distr_list[iens];
      afp_corr_by_FSEs_B.distr_list.push_back( fp_CDH_B.distr_list[iens]/(1.0 + (Al+ Al_B)*Csi_CDH_B.distr_list[iens]*distr_t::f_of_distr(FVE, MpL_distr)));
      afp_uncorr_B.distr_list.push_back( fp_CDH_B.distr_list[iens]); 
    }

    for(int iens=0;iens<Nens_A;iens++) {
      distr_t MpL_distr = L_A_ens[iens]*Mpi_CDH_A.distr_list[iens];
      afp_corr_by_FSEs_A.distr_list.push_back( fp_CDH_A.distr_list[iens]/(1.0 + (Al+ Al_A)*Csi_CDH_A.distr_list[iens]*distr_t::f_of_distr(FVE, MpL_distr)));
      afp_uncorr_A.distr_list.push_back( fp_CDH_A.distr_list[iens]); 
    }
    
    //Print corrected lattice points to file

    //B ensemble
    Print_To_File({Ensemble_B_tag_list},{afp_corr_by_FSEs_B.ave(), afp_corr_by_FSEs_B.err(), afp_uncorr_B.ave(), afp_uncorr_B.err(), Csi_CDH_B.ave(), Csi_CDH_B.err()},  "../data/scale_setting/data_FSEs_corrected_afp_B.dat", "", "# Ens   corrected   uncorrected    xi ");

    //A ensemble
    Print_To_File({Ensemble_A_tag_list},{afp_corr_by_FSEs_A.ave(), afp_corr_by_FSEs_A.err(), afp_uncorr_A.ave(), afp_uncorr_A.err(), Csi_CDH_A.ave(), Csi_CDH_A.err()},  "../data/scale_setting/data_FSEs_corrected_afp_A.dat", "", "# Ens   corrected   uncorrected     xi");
      


    //get lattice spacings

    //get id of C and D ensemble
    int C_id, D_id;
    bool find_C_id=false;
    bool find_D_id=false;
    for(int iens=0;iens<Nens_phys;iens++) {
      if(Ensemble_phys_point_tag_list[iens].substr(1,1) == "C") { find_C_id = true; C_id =iens;}
      if(Ensemble_phys_point_tag_list[iens].substr(1,1) == "D") { find_D_id = true; D_id =iens;}
    }

    if(!find_C_id || !find_D_id) crash("Cannot find phyisical point ensemble C and D in Ensemble_phys_point_tag_list");
      

 
    //compute Am and Am_art coefficients
    distr_t Am_cont, Am_art, A2m_cont, A2m_art;
    if(Use_three_finest_in_scale_setting_fp) {Am_cont = Am_B; Am_art = 0.0*Am_B; A2m_cont = A2m_B; A2m_art = 0.0*A2m_B;}
    else {
      Am_cont = ( Am_B*(afp_A*afp_A) - Am_A*(afp_B*afp_B) )/( afp_A*afp_A - afp_B*afp_B);
      Am_art = (Am_B - Am_A)/( afp_B*afp_B  - afp_A*afp_A);
      if( fabs(A2m_A.ave()) < 1e-16) { A2m_cont = A2m_B; A2m_art = 0.0*A2m_B;}
      else {
      A2m_cont = ( A2m_B*(afp_A*afp_A) - A2m_A*(afp_B*afp_B) )/( afp_A*afp_A - afp_B*afp_B);
      A2m_art = (A2m_B - A2m_A)/( afp_B*afp_B  - afp_A*afp_A);
      }
    }

    //correct for xi_p mistuning on C and D ensembles

    distr_t corr_C_mass, corr_D_mass;

    corr_C_mass = 1.0  -2.0*Csi_CDH_phys_point[C_id]*distr_t::f_of_distr(LOG, Csi_CDH_phys_point[C_id]/csi_phys_distr)  + (Am_cont+ Am_art*fp_CDH_phys_point[C_id]*fp_CDH_phys_point[C_id])*(Csi_CDH_phys_point[C_id] - csi_phys_distr) +  (A2m_cont + A2m_art*fp_CDH_phys_point[C_id]*fp_CDH_phys_point[C_id])*(Csi_CDH_phys_point[C_id]*Csi_CDH_phys_point[C_id] - csi_phys_distr*csi_phys_distr);

    corr_D_mass = 1.0  -2.0*Csi_CDH_phys_point[D_id]*distr_t::f_of_distr(LOG, Csi_CDH_phys_point[D_id]/csi_phys_distr)  + (Am_cont  + Am_art*fp_CDH_phys_point[D_id]*fp_CDH_phys_point[D_id])*(Csi_CDH_phys_point[D_id] - csi_phys_distr) +  (A2m_cont + A2m_art*fp_CDH_phys_point[D_id]*fp_CDH_phys_point[D_id])*(Csi_CDH_phys_point[D_id]*Csi_CDH_phys_point[D_id] - csi_phys_distr*csi_phys_distr);


    //correct for FSEs on C and D ensembles
    
    //get MpiL for C and D ensembles

    distr_t MpL_distr_C = L_phys_point[C_id]*Mpi_CDH_phys_point[C_id];
    distr_t MpL_distr_D = L_phys_point[D_id]*Mpi_CDH_phys_point[D_id];

    
    distr_t corr_C_FSEs, corr_D_FSEs;

    corr_C_FSEs = (1.0 + (Al+ Al_B)*distr_t::f_of_distr(pow_2,Csi_CDH_phys_point.distr_list[C_id])*distr_t::f_of_distr(FVE, MpL_distr_C));
    corr_D_FSEs = (1.0 + (Al+ Al_B)*distr_t::f_of_distr(pow_2,Csi_CDH_phys_point.distr_list[D_id])*distr_t::f_of_distr(FVE, MpL_distr_D));
    

    
    

    a_distr_A = afp_A/fpi_phys_distr;
    
    a_distr_B = afp_B/fpi_phys_distr;

    a_distr_C = fp_CDH_phys_point[C_id]/(corr_C_FSEs*corr_C_mass)/fpi_phys_distr;

    a_distr_D = fp_CDH_phys_point[D_id]/(corr_D_FSEs*corr_D_mass)/fpi_phys_distr;



    //print on screen
    //print afC and afD before and after corrections:
    cout<<"C ensemble"<<endl;
    cout<<"af_C: "<<fpi_phys_point.distr_list[C_id].ave()<<" +- "<<fpi_phys_point.distr_list[C_id].err()<<endl;
    cout<<"af_C (CDH-corrected): "<<fp_CDH_phys_point[C_id].ave()<<" +- "<<fp_CDH_phys_point[C_id].err()<<endl;
    cout<<"af_C (infL): "<<(fp_CDH_phys_point[C_id]/(corr_C_FSEs)).ave()<<" +- "<<(fp_CDH_phys_point[C_id]/corr_C_FSEs).err()<<endl;
    cout<<"af_C (infL-mass corrected): "<<(fp_CDH_phys_point[C_id]/(corr_C_FSEs*corr_C_mass)).ave()<<" +- "<<(fp_CDH_phys_point[C_id]/(corr_C_FSEs*corr_C_mass)).err()<<endl;
    cout<<"D ensemble"<<endl;
    cout<<"af_D: "<<fpi_phys_point.distr_list[D_id].ave()<<" +- "<<fpi_phys_point.distr_list[D_id].err()<<endl;
    cout<<"af_D (CDH-corrected): "<<fp_CDH_phys_point[D_id].ave()<<" +- "<<fp_CDH_phys_point[D_id].err()<<endl;
    cout<<"af_D (infL): "<<(fp_CDH_phys_point[D_id]/(corr_D_FSEs)).ave()<<" +- "<<(fp_CDH_phys_point[D_id]/corr_D_FSEs).err()<<endl;
    cout<<"af_D (infL-mass corrected): "<<(fp_CDH_phys_point[D_id]/(corr_D_FSEs*corr_D_mass)).ave()<<" +- "<<(fp_CDH_phys_point[D_id]/(corr_D_FSEs*corr_D_mass)).err()<<endl;   
    cout<<"#######  Printing lattice spacing determined from afp analysis ########"<<endl;
    cout<<"A ensemble: "<< a_distr_A.ave()<<" +- "<<a_distr_A.err()<<" [GeV-1]    "<<(a_distr_A/fm_to_iGev).ave()<<" +- "<<(a_distr_A/fm_to_iGev).err()<<" [fm] "<<endl;
    cout<<"B ensemble: "<< a_distr_B.ave()<<" +- "<<a_distr_B.err()<<" [GeV-1]    "<<(a_distr_B/fm_to_iGev).ave()<<" +- "<<(a_distr_B/fm_to_iGev).err()<<" [fm] "<<endl;
    cout<<"C ensemble: "<< a_distr_C.ave()<<" +- "<<a_distr_C.err()<<" [GeV-1]    "<<(a_distr_C/fm_to_iGev).ave()<<" +- "<<(a_distr_C/fm_to_iGev).err()<<" [fm] "<<endl;
    cout<<"D ensemble: "<< a_distr_D.ave()<<" +- "<<a_distr_D.err()<<" [GeV-1]    "<<(a_distr_D/fm_to_iGev).ave()<<" +- "<<(a_distr_D/fm_to_iGev).err()<<" [fm] "<<endl;


    //get scaling w.r.t. previous analysis

    LatticeInfo a_info;
    double a_A_ave, a_A_err, a_B_ave, a_B_err, a_C_ave, a_C_err, a_D_ave, a_D_err;
    a_info.LatInfo_new_ens("cA211a.53.24");
    a_A_ave= a_info.a;
    a_A_err= a_info.a_err;
    a_info.LatInfo_new_ens("cB211b.072.64");
    a_B_ave= a_info.a;
    a_B_err= a_info.a_err;
    a_info.LatInfo_new_ens("cC211a.06.80");
    a_C_ave= a_info.a;
    a_C_err= a_info.a_err;
    a_info.LatInfo_new_ens("cD211a.054.96");
    a_D_ave= a_info.a;
    a_D_err= a_info.a_err;

    distr_t a_A_old(UseJack),  a_B_old(UseJack), a_C_old(UseJack), a_D_old(UseJack);
    for(int ijack=0;ijack<Njacks;ijack++) {
      a_A_old.distr.push_back( fm_to_iGev*( a_A_ave + GM()*a_A_err*(1.0/sqrt(Njacks-1.0))));
      a_B_old.distr.push_back( fm_to_iGev*( a_B_ave + GM()*a_B_err*(1.0/sqrt(Njacks-1.0))));
      a_C_old.distr.push_back( fm_to_iGev*( a_C_ave + GM()*a_C_err*(1.0/sqrt(Njacks-1.0))));
      a_D_old.distr.push_back( fm_to_iGev*( a_D_ave + GM()*a_D_err*(1.0/sqrt(Njacks-1.0))));
    }


    //get relative difference
    
    distr_t_list rel_lattice(UseJack);
    
    distr_t_list new_lat(UseJack);
  
    rel_lattice.distr_list.push_back(a_A_old/a_distr_A -1.0);
    rel_lattice.distr_list.push_back(a_B_old/a_distr_B -1.0);
    rel_lattice.distr_list.push_back(a_C_old/a_distr_C -1.0);
    rel_lattice.distr_list.push_back(a_D_old/a_distr_D -1.0);
    new_lat.distr_list.push_back(a_distr_A/fm_to_iGev);
    new_lat.distr_list.push_back(a_distr_B/fm_to_iGev);
    new_lat.distr_list.push_back(a_distr_C/fm_to_iGev);
    new_lat.distr_list.push_back(a_distr_D/fm_to_iGev);
  

    //print the result
    Print_To_File({}, {rel_lattice.ave(), rel_lattice.err(), new_lat.ave() , new_lat.err()}, "../data/scale_setting/a2_scaling.dat", "", "# rel  a(new)");
    //Print_To_File({}, {rel_lattice.ave(), rel_lattice.err(), new_lat.ave() , new_lat.err()}, "../data/scale_setting/a2_scaling_three_finest.dat", "", "# rel  a(new)");

    //fit rel_lattice with a linear ansatz

    bootstrap_fit<fpar_a2, rel_lat_cl> bf_rel(Njacks);
    bf_rel.set_warmup_lev(2); //sets warmup

    int Nm = rel_lattice.size();

    bf_rel.Set_number_of_measurements(Nm);
    bf_rel.Set_verbosity(false);
    bf_rel.Add_par("a0", 0.0, 1e-2);
    bf_rel.Add_par("D", 0.1, 1e-2);

    bf_rel.ansatz =  [&](const fpar_a2 &p, const rel_lat_cl &ip) -> double {
		   return p.a0 + p.D*(ip.X*ip.X);
	       };

    bf_rel.measurement = [&](const fpar_a2 &p, const rel_lat_cl &ip) ->double {
		     return ip.Y;
		   };
    bf_rel.error = [&](const fpar_a2 &p, const rel_lat_cl &ip) -> double {
	       return ip.Y_err;
	     };


    //meas and errs

  
    //fill the data
    vector<vector<rel_lat_cl>> data(Njacks);
    for(auto &data_ij: data) data_ij.resize(Nm);
    //allocate space for output result
    boot_fit_data<fpar_a2> Bt_fit;
  
    for(auto &data_iboot: data) data_iboot.resize(4);
    
    for(int ijack=0;ijack<Njacks;ijack++) {
      for(int i=0;i<Nm;i++) {
	data[ijack][i].Y = rel_lattice.distr_list[i].distr[ijack];
	data[ijack][i].Y_err = rel_lattice.err(i);
	data[ijack][i].X = new_lat.distr_list[i].distr[ijack];
      } 
    }


    
    //append
    bf_rel.Append_to_input_par(data);
    
    //fit
    Bt_fit= bf_rel.Perform_bootstrap_fit();
    Bt_fit.ch2_ave();
  

  //retrieve fit parameter a0 and D
  distr_t a0_fit(UseJack), D_fit(UseJack);

  for(int ijack=0;ijack<Njacks;ijack++) {a0_fit.distr.push_back( Bt_fit.par[ijack].a0); D_fit.distr.push_back(Bt_fit.par[ijack].D);}

 
  int Npoints = 100;
  double xmin= 0.0;
  double xmax= 1.5*a_distr_A.ave()/fm_to_iGev;
  Vfloat Xpoints;
  distr_t_list Fit_func;
  for(int ipoint=0;ipoint<Npoints;ipoint++) {
    double xp= xmin+ ipoint*(xmax-xmin)/((double)Npoints);
    Fit_func.distr_list.push_back(a0_fit + D_fit*(xp*xp));
    Xpoints.push_back(xp);
    }

  Print_To_File({}, {Xpoints, Fit_func.ave(), Fit_func.err()},"../data/scale_setting/fit_a2_scaling.dat", "", "");
  //Print_To_File({}, {Xpoints, Fit_func.ave(), Fit_func.err()},"../data/scale_setting/fit_a2_scaling_three_finest.dat", "", "");
    
    
		  
  

  return;
}







































void Determine_scale_from_fp_FLAG(const distr_t_list &Mpi_phys_point, const distr_t_list &fpi_phys_point, const distr_t_list &Mpi_Bens, const distr_t_list &fpi_Bens, const distr_t_list &Mpi_Aens, const distr_t_list &fpi_Aens,  const vector<double> &L_phys_point, const vector<double> &L_B_ens, const vector<double> &L_A_ens, const vector<string>  &Ensemble_phys_point_tag_list, const vector<string>  &Ensemble_B_tag_list, const vector<string>  &Ensemble_A_tag_list, distr_t &a_distr_A, distr_t &a_distr_B, distr_t &a_distr_C, distr_t &a_distr_D, distr_t &a_distr_E, bool UseJack,bool Use_three_finest_in_scale_setting_fp) {


  cout<<"Starting scale setting analysis using afp"<<endl;

  
  //init Gaussian number generator
  GaussianMersenne GM(76795699);


  int Nens_phys = L_phys_point.size();
  int Nens_B = L_B_ens.size();
  int Nens_A = L_A_ens.size();
  if(Nens_phys==0 || Nens_B==0 || Nens_A==0) crash("in Determine_scale_setting_from_fp Nens_phys=0 || Nens_B = 0 || Nens_A = 0");
  int Njacks= Mpi_phys_point.distr_list[0].size();


  //generate K_\ell correcting factor for cA211a.12.48 ensemble according to 2104.06747
  distr_t Kl(UseJack);
  for(int ijack=0;ijack<Njacks;ijack++) { Kl.distr.push_back( sqrt( 1.0 + pow( corr_fact + GM()*corr_fact_err/sqrt(Njacks -1.0)    ,2)));} 
  

  
  //List of lambda functions needed

  auto FVE = [](double x) -> double {return (1.0/pow(x,1.5))*exp(-x);};
  auto LOG = [](double x) { return log(x);};
  auto pow_1_5 = [](double x) -> double { return pow(x, 1.0/5.0);};
  auto pow_4 = [](double x) -> double { return pow(x, 4.0);};
  auto pow_2 = [](double x) -> double { return pow(x, 2.0);};
  auto pow_3_2 = [](double x) -> double { return pow(x, 3.0/2.0);};
  auto EXP = [](double x) -> double {return exp(x);};



  //###########################



  //Output directories

  boost::filesystem::create_directory("../data/scale_setting");
  
  //##########################


  //Generate fake distribution for csi_phys, Mpi , fp, Xpi
 
  distr_t Mpi_phys_distr(UseJack);
  distr_t fpi_phys_distr(UseJack);
  distr_t csi_phys_distr(UseJack);
  distr_t Xpi_phys_distr(UseJack);
  
  for(int ij=0;ij<Njacks;ij++) {
    double Mfak = Mp_phys + GM()*Mp_phys_err/sqrt(Njacks -1.0);
    double Ffak = fp_phys + GM()*fp_phys_err/sqrt(Njacks -1.0);

    Mpi_phys_distr.distr.push_back( Mfak);
    fpi_phys_distr.distr.push_back(Ffak);
    csi_phys_distr.distr.push_back( pow(Mfak/(4.0*M_PI*Ffak),2.0));
    Xpi_phys_distr.distr.push_back( pow( Ffak*pow(Mfak,4.0), 1.0/5.0));
  }

 
  distr_t_list fp_CDH_phys_point(UseJack), Mpi_CDH_phys_point(UseJack), Xpi_phys_point(UseJack), Csi_CDH_phys_point(UseJack);
  distr_t_list fp_CDH_B(UseJack), Mpi_CDH_B(UseJack), Xpi_B(UseJack), Csi_CDH_B(UseJack);
  distr_t_list fp_CDH_A(UseJack), Mpi_CDH_A(UseJack), Xpi_A(UseJack), Csi_CDH_A(UseJack);

   
  //compute CDH corrections to Mpi and fpi
  for(int iens=0; iens<Nens_phys;iens++) {
    //GL and CDH correction for Mpi and fp
    distr_t csi_L = Mpi_phys_point.distr_list[iens]*Mpi_phys_point.distr_list[iens]/(pow(4.0*M_PI,2)*fpi_phys_point.distr_list[iens]*fpi_phys_point.distr_list[iens]);
    distr_t g1 = distr_t::f_of_distr(g1_l, Mpi_phys_point.distr_list[iens]*L_phys_point[iens]);
    distr_t g2 = distr_t::f_of_distr(g2_l, Mpi_phys_point.distr_list[iens]*L_phys_point[iens]);
    distr_t log_l = log(csi_phys) - distr_t::f_of_distr(LOG, csi_L);
    fp_CDH_phys_point.distr_list.push_back(fpi_phys_point.distr_list[iens]/(1.0 -2.0*csi_L*g1 +2.0*csi_L*csi_L*( (Cf1(l1ph,l2ph,l3ph,l4ph) + Sf1(s0,s1,s2,s3) + Cf1_log()*log_l)*g1 + (Cf2(l1ph,l2ph,l3ph,l4ph) + Sf2(s0,s1,s2,s3) + Cf2_log()*log_l)*g2)));
    Mpi_CDH_phys_point.distr_list.push_back(Mpi_phys_point.distr_list[iens]/(1.0 + 0.5*csi_L*g1 - csi_L*csi_L*( (Cm1(l1ph,l2ph,l3ph,l4ph) + Sm1(s0,s1,s2,s3) + Cm1_log()*log_l)*g1 + (Cm2(l1ph,l2ph,l3ph,l4ph) + Sm2(s0,s1,s2,s3) + Cm2_log()*log_l)*g2)));
    Csi_CDH_phys_point.distr_list.push_back( Mpi_CDH_phys_point.distr_list[iens]*Mpi_CDH_phys_point.distr_list[iens]/(16.0*M_PI*M_PI*fp_CDH_phys_point.distr_list[iens]*fp_CDH_phys_point.distr_list[iens]));
    Xpi_phys_point.distr_list.push_back( distr_t::f_of_distr(pow_1_5, fp_CDH_phys_point.distr_list[iens]*distr_t::f_of_distr(pow_4, Mpi_CDH_phys_point.distr_list[iens])));
  }


  int cB211a2548_tag=0;
   
  for(int iens=0; iens<Nens_B;iens++) {

    //GL and CDH correction for Mpi and fp
    distr_t csi_L = Mpi_Bens.distr_list[iens]*Mpi_Bens.distr_list[iens]/(pow(4.0*M_PI,2)*fpi_Bens.distr_list[iens]*fpi_Bens.distr_list[iens]);
    distr_t g1 = distr_t::f_of_distr(g1_l, Mpi_Bens.distr_list[iens]*L_B_ens[iens]);
    distr_t g2 = distr_t::f_of_distr(g2_l, Mpi_Bens.distr_list[iens]*L_B_ens[iens]);
    distr_t log_l = log(csi_phys) - distr_t::f_of_distr(LOG, csi_L);
    fp_CDH_B.distr_list.push_back(fpi_Bens.distr_list[iens]/(1.0 -2.0*csi_L*g1 +2.0*csi_L*csi_L*( (Cf1(l1ph,l2ph,l3ph,l4ph) + Sf1(s0,s1,s2,s3) + Cf1_log()*log_l)*g1 + (Cf2(l1ph,l2ph,l3ph,l4ph) + Sf2(s0,s1,s2,s3) + Cf2_log()*log_l)*g2)));
    Mpi_CDH_B.distr_list.push_back(Mpi_Bens.distr_list[iens]/(1.0 + 0.5*csi_L*g1 - csi_L*csi_L*( (Cm1(l1ph,l2ph,l3ph,l4ph) + Sm1(s0,s1,s2,s3) + Cm1_log()*log_l)*g1 + (Cm2(l1ph,l2ph,l3ph,l4ph) + Sm2(s0,s1,s2,s3) + Cm2_log()*log_l)*g2)));
    Csi_CDH_B.distr_list.push_back( Mpi_CDH_B.distr_list[iens]*Mpi_CDH_B.distr_list[iens]/(16.0*M_PI*M_PI*fp_CDH_B.distr_list[iens]*fp_CDH_B.distr_list[iens]));
    Xpi_B.distr_list.push_back( distr_t::f_of_distr(pow_1_5, fp_CDH_B.distr_list[iens]*distr_t::f_of_distr(pow_4, Mpi_CDH_B.distr_list[iens])));
    if(Ensemble_B_tag_list[iens] == "cB211a.25.48") cB211a2548_tag = iens;

  }

  for(int iens=0; iens< Nens_B;iens++) {
    if( Ensemble_B_tag_list[iens] == "cB211a.25.24" || Ensemble_B_tag_list[iens] == "cB211a.25.32") Csi_CDH_B.distr_list[iens] = Csi_CDH_B.distr_list[cB211a2548_tag];
  }

 
  for(int iens=0; iens<Nens_A;iens++) {

    //GL and CDH correction for Mpi and fp
    distr_t csi_L = Mpi_Aens.distr_list[iens]*Mpi_Aens.distr_list[iens]/(pow(4.0*M_PI,2)*fpi_Aens.distr_list[iens]*fpi_Aens.distr_list[iens]);
    distr_t g1 = distr_t::f_of_distr(g1_l, Mpi_Aens.distr_list[iens]*L_A_ens[iens]);
    distr_t g2 = distr_t::f_of_distr(g2_l, Mpi_Aens.distr_list[iens]*L_A_ens[iens]);
    distr_t log_l = log(csi_phys) - distr_t::f_of_distr(LOG, csi_L);
    fp_CDH_A.distr_list.push_back(fpi_Aens.distr_list[iens]/(1.0 -2.0*csi_L*g1 +2.0*csi_L*csi_L*( (Cf1(l1ph,l2ph,l3ph,l4ph) + Sf1(s0,s1,s2,s3) + Cf1_log()*log_l)*g1 + (Cf2(l1ph,l2ph,l3ph,l4ph) + Sf2(s0,s1,s2,s3) + Cf2_log()*log_l)*g2)));
    //if CA211a.12.48 correct fpi for violation of maximal twist condition
    if(Ensemble_A_tag_list[iens] == "cA211a.12.48") {
      fp_CDH_A.distr_list[iens] = fp_CDH_A.distr_list[iens]*Kl;
    }
    Mpi_CDH_A.distr_list.push_back(Mpi_Aens.distr_list[iens]/(1.0 + 0.5*csi_L*g1 - csi_L*csi_L*( (Cm1(l1ph,l2ph,l3ph,l4ph) + Sm1(s0,s1,s2,s3) + Cm1_log()*log_l)*g1 + (Cm2(l1ph,l2ph,l3ph,l4ph) + Sm2(s0,s1,s2,s3) + Cm2_log()*log_l)*g2)));
    Csi_CDH_A.distr_list.push_back( Mpi_CDH_A.distr_list[iens]*Mpi_CDH_A.distr_list[iens]/(16.0*M_PI*M_PI*fp_CDH_A.distr_list[iens]*fp_CDH_A.distr_list[iens]));
    Xpi_A.distr_list.push_back( distr_t::f_of_distr(pow_1_5, fp_CDH_A.distr_list[iens]*distr_t::f_of_distr(pow_4, Mpi_CDH_A.distr_list[iens])));

  }


  //Print_To_File
  Print_To_File({Ensemble_phys_point_tag_list}, {L_phys_point, Mpi_CDH_phys_point.ave(), Mpi_CDH_phys_point.err(), fp_CDH_phys_point.ave(), fp_CDH_phys_point.err(), Xpi_phys_point.ave(), Xpi_phys_point.err(), Csi_CDH_phys_point.ave(), Csi_CDH_phys_point.err()}, "../data/scale_setting/pion_obs_phys_point.dat", "", "#Ens L Mpi fp Xpi  Csi");

  Print_To_File({Ensemble_B_tag_list}, {L_B_ens, Mpi_CDH_B.ave(), Mpi_CDH_B.err(), fp_CDH_B.ave(), fp_CDH_B.err(), Xpi_B.ave(), Xpi_B.err(), Csi_CDH_B.ave(), Csi_CDH_B.err()}, "../data/scale_setting/pion_obs_B_ens.dat", "", "#Ens L Mpi fp Xpi  Csi");

  Print_To_File({Ensemble_A_tag_list}, {L_A_ens, Mpi_CDH_A.ave(), Mpi_CDH_A.err(), fp_CDH_A.ave(), fp_CDH_A.err(), Xpi_A.ave(), Xpi_A.err(), Csi_CDH_A.ave(), Csi_CDH_A.err()}, "../data/scale_setting/pion_obs_A_ens.dat", "", "#Ens L Mpi fp Xpi  Csi");

  


  //Perform extrapolation of afp

   bootstrap_fit<sc_fp_fpar,sc_fp_ipar> bf(Njacks);
   bf.set_warmup_lev(1);
   bf.Set_number_of_measurements(Nens_B+Nens_A);
   bf.Set_verbosity(scale_setting_fit_verbosity);
   bf.Set_print_path("chi2_scale_setting_fp");


   bootstrap_fit<sc_fp_fpar,sc_fp_ipar> bf_ch2(1);
   bf_ch2.set_warmup_lev(1);
   bf_ch2.Set_number_of_measurements(Nens_B+Nens_A);
   bf_ch2.Set_verbosity(scale_setting_fit_verbosity);
  
   //B ens pars
   bf_ch2.Add_par("afp_B", 0.053, 0.001);
   bf_ch2.Add_par("Am_B", 10.0, 0.1);
   bf_ch2.Add_par("A2m_B", 1.0, 0.1);
   bf_ch2.Add_par("Al_B", 10.0, 0.1);

   //A ens pars
   bf_ch2.Add_par("afp_A", 0.060, 0.001);
   bf_ch2.Add_par("Am_A", 11.0, 0.1);
   bf_ch2.Add_par("A2m_A", 1.0, 0.1);
   bf_ch2.Add_par("Al_A", 10.0, 0.1);

   bf_ch2.Add_par("Al", 10.0, 0.1);

  

   //Fix parameters depending on the fit type
   if(!Enable_A2m_cont) bf_ch2.Fix_par("A2m_B", 0.0);
   bf_ch2.Fix_par("A2m_A", 0.0);


   if(Use_three_finest_in_scale_setting_fp) {
     //fix parameters for A ensembles
     bf_ch2.Fix_par("Al",0.0);
   }
   else {
     bf_ch2.Fix_par("Al_B",0.0);
     bf_ch2.Fix_par("Al_A",0.0);

   }



   //B ens pars
   bf.Add_par("afp_B", 0.053, 0.001);
   bf.Add_par("Am_B", 10.0, 0.1);
   bf.Add_par("A2m_B", 1.0, 0.1);
   bf.Add_par("Al_B", 10.0, 0.1);

   //A ens pars
   bf.Add_par("afp_A", 0.060, 0.001);
   bf.Add_par("Am_A", 11.0, 0.1);
   bf.Add_par("A2m_A", 1.0, 0.1);
   bf.Add_par("Al_A", 10.0, 0.1);

   bf.Add_par("Al", 10.0, 0.1);

  

   //Fix parameters depending on the fit type
   if(!Enable_A2m_cont) bf.Fix_par("A2m_B", 0.0);
   bf.Fix_par("A2m_A", 0.0);


   if(Use_three_finest_in_scale_setting_fp) {
     //fix parameters for A ensembles
     bf.Fix_par("Al",0.0);
   }
   else {
     bf.Fix_par("Al_B",0.0);
     bf.Fix_par("Al_A",0.0);

   }



   
  
   //set ansatz
   bf.ansatz= [=](const sc_fp_fpar& X, const sc_fp_ipar& Y) {

		double ret_val=0;
		double MpL= Y.Mp*Y.L;

		if(Y.Is_B) {

		  if(Use_fpi_or_Xpi) ret_val = X.afp_B*(1.0  -2.0*Y.csi*log(Y.csi/Y.csi_ph) + X.Am_B*(Y.csi - Y.csi_ph) + X.A2m_B*(pow(Y.csi,2) - pow(Y.csi_ph,2)))*(1.0 + (X.Al+ X.Al_B)*Y.csi*exp(-MpL)/pow(MpL, 3.0/2.0));
		  else crash("Xpi not yet implemented"); 

		}

		else {

		  if(Use_fpi_or_Xpi) ret_val =  X.afp_A*(1.0  -2.0*Y.csi*log(Y.csi/Y.csi_ph) + X.Am_A*(Y.csi - Y.csi_ph) + X.A2m_B*(pow(Y.csi,2) - pow(Y.csi_ph,2)))*(1.0 + (X.Al+ X.Al_A)*Y.csi*exp(-MpL)/pow(MpL, 3.0/2.0));
		  else crash("Xpi not yet implemented") ;

		}

		return ret_val;
		
			     
	      };

   
   bf.measurement = [=](const sc_fp_fpar& X, const sc_fp_ipar& Y) {
		      double ret_val= Use_fpi_or_Xpi?Y.fp:Y.Xpi;
		      return ret_val;
		    };
   bf.error = [=](const sc_fp_fpar& X, const sc_fp_ipar& Y) {
		double ret_val;
		ret_val = Use_fpi_or_Xpi?Y.fp_err:Y.Xpi_err;
		return ret_val;
	      };



   bf_ch2.ansatz= bf.ansatz;
   bf_ch2.measurement = bf.measurement;
   bf_ch2.error= bf.error;



   //populate bf

   int Nens_tot = (Nens_B + Nens_A);
    vector<vector<sc_fp_ipar>> ipar_all_ens(Njacks);
    for(auto &ipar_jack: ipar_all_ens) ipar_jack.resize(Nens_tot);


    vector<vector<sc_fp_ipar>> ipar_all_ens_ch2(1);
    for(auto &ipar_jack: ipar_all_ens_ch2) ipar_jack.resize(Nens_tot);
    
    for(int iens=0; iens<Nens_tot;iens++) {  
      for(int ijack=0;ijack<Njacks;ijack++) {
	
	if(iens < Nens_B) {

	  ipar_all_ens[ijack][iens].Mp = Mpi_CDH_B.distr_list[iens].distr[ijack]; 
	  ipar_all_ens[ijack][iens].Mp_err = Mpi_CDH_B.err(iens);
	  
	
	  ipar_all_ens[ijack][iens].fp = fp_CDH_B.distr_list[iens].distr[ijack]; 
	  ipar_all_ens[ijack][iens].fp_err = fp_CDH_B.err(iens);
	
	
	  ipar_all_ens[ijack][iens].Xpi = Xpi_B.distr_list[iens].distr[ijack]; 
	  ipar_all_ens[ijack][iens].Xpi_err = Xpi_B.err(iens);
	
	
	  ipar_all_ens[ijack][iens].csi= Csi_CDH_B.distr_list[iens].distr[ijack]; 
	  ipar_all_ens[ijack][iens].csi_err = Csi_CDH_B.err(iens);
	
			
	  ipar_all_ens[ijack][iens].L = L_B_ens[iens]; 
	  ipar_all_ens[ijack][iens].csi_ph = csi_phys_distr.distr[ijack];

	  ipar_all_ens[ijack][iens].Is_B = true;

	  if(ijack==0) {

	    ipar_all_ens_ch2[ijack][iens].Mp = Mpi_CDH_B.distr_list[iens].ave(); 
	    ipar_all_ens_ch2[ijack][iens].Mp_err = Mpi_CDH_B.err(iens);
	  
	
	    ipar_all_ens_ch2[ijack][iens].fp = fp_CDH_B.distr_list[iens].ave(); 
	    ipar_all_ens_ch2[ijack][iens].fp_err = fp_CDH_B.err(iens);
	  
	
	    ipar_all_ens_ch2[ijack][iens].Xpi = Xpi_B.distr_list[iens].ave(); 
	    ipar_all_ens_ch2[ijack][iens].Xpi_err = Xpi_B.err(iens);
	  
	
	    ipar_all_ens_ch2[ijack][iens].csi= Csi_CDH_B.distr_list[iens].ave(); 
	    ipar_all_ens_ch2[ijack][iens].csi_err = Csi_CDH_B.err(iens);
	  
			
	    ipar_all_ens_ch2[ijack][iens].L = L_B_ens[iens]; 
	    ipar_all_ens_ch2[ijack][iens].csi_ph = csi_phys_distr.ave();

	    ipar_all_ens_ch2[ijack][iens].Is_B = true;


	  }

	}

	else {

	  ipar_all_ens[ijack][iens].Mp = Mpi_CDH_A.distr_list[iens-Nens_B].distr[ijack]; 
	  ipar_all_ens[ijack][iens].Mp_err = Mpi_CDH_A.err(iens-Nens_B);

	  ipar_all_ens[ijack][iens].fp = fp_CDH_A.distr_list[iens-Nens_B].distr[ijack]; 
	  ipar_all_ens[ijack][iens].fp_err = fp_CDH_A.err(iens-Nens_B);
	
	
	  ipar_all_ens[ijack][iens].Xpi = Xpi_A.distr_list[iens-Nens_B].distr[ijack]; 
	  ipar_all_ens[ijack][iens].Xpi_err = Xpi_A.err(iens-Nens_B);
	
	
	  ipar_all_ens[ijack][iens].csi= Csi_CDH_A.distr_list[iens-Nens_B].distr[ijack]; 
	  ipar_all_ens[ijack][iens].csi_err = Csi_CDH_A.err(iens-Nens_B);
	
			
	  ipar_all_ens[ijack][iens].L = L_A_ens[iens-Nens_B]; 
	  ipar_all_ens[ijack][iens].csi_ph = csi_phys_distr.distr[ijack];

	  ipar_all_ens[ijack][iens].Is_B = false;


	  if(ijack==0) {

	    ipar_all_ens_ch2[ijack][iens].Mp = Mpi_CDH_A.distr_list[iens-Nens_B].ave(); 
	    ipar_all_ens_ch2[ijack][iens].Mp_err = Mpi_CDH_A.err(iens-Nens_B);
	    
	    ipar_all_ens_ch2[ijack][iens].fp = fp_CDH_A.distr_list[iens-Nens_B].ave(); 
	    ipar_all_ens_ch2[ijack][iens].fp_err = fp_CDH_A.err(iens-Nens_B);
	
	
	    ipar_all_ens_ch2[ijack][iens].Xpi = Xpi_A.distr_list[iens-Nens_B].ave(); 
	    ipar_all_ens_ch2[ijack][iens].Xpi_err = Xpi_A.err(iens-Nens_B);
	  
	
	    ipar_all_ens_ch2[ijack][iens].csi= Csi_CDH_A.distr_list[iens-Nens_B].ave(); 
	    ipar_all_ens_ch2[ijack][iens].csi_err = Csi_CDH_A.err(iens-Nens_B);
	      
	    
	    ipar_all_ens_ch2[ijack][iens].L = L_A_ens[iens-Nens_B]; 
	    ipar_all_ens_ch2[ijack][iens].csi_ph = csi_phys_distr.ave();
	      
	    ipar_all_ens_ch2[ijack][iens].Is_B = false;


	  }

	}
	
      }
      
    }

      
    bf.Append_to_input_par(ipar_all_ens);
    //fit
    boot_fit_data<sc_fp_fpar> Bt_sc_fit = bf.Perform_bootstrap_fit();


    bf_ch2.Append_to_input_par(ipar_all_ens_ch2);
    //fit
    boot_fit_data<sc_fp_fpar> Bt_sc_fit_ch2 = bf_ch2.Perform_bootstrap_fit();


    //retrieve parameters
    distr_t afp_B, Am_B, A2m_B, Al_B;  //constructor sets UseJack to 1
    distr_t afp_A, Am_A, A2m_A, Al_A;
    distr_t Al;

  
 
    for(int ijack=0;ijack<Njacks;ijack++) {

      sc_fp_fpar my_w0_fit_pars = Bt_sc_fit.par[ijack];
      afp_B.distr.push_back(my_w0_fit_pars.afp_B);
      Am_B.distr.push_back(my_w0_fit_pars.Am_B);
      A2m_B.distr.push_back(my_w0_fit_pars.A2m_B);
      Al_B.distr.push_back(my_w0_fit_pars.Al_B);
      afp_A.distr.push_back(my_w0_fit_pars.afp_A);
      Am_A.distr.push_back(my_w0_fit_pars.Am_A);
      A2m_A.distr.push_back(my_w0_fit_pars.A2m_A);
      Al_A.distr.push_back(my_w0_fit_pars.Al_A);
      Al.distr.push_back(my_w0_fit_pars.Al);
      
    }

    //get average ch2
    double ch2  = Bt_sc_fit_ch2.get_ch2_ave();
    



    //cout infos

    cout<<"########  scale setting fit parameter : "<<" ##########"<<endl;
    cout<<"Use 3 finest?: "<<Use_three_finest_in_scale_setting_fp<<endl;
    cout<<"ch2 total: "<<ch2<<endl;
    cout<<"afp_B: "<<afp_B.ave()<<" +- "<<afp_B.err()<<endl;
    cout<<"Am_B: "<<Am_B.ave()<<" +- "<<Am_B.err()<<endl;
    cout<<"A2m_B: "<<A2m_B.ave()<<" +- "<<A2m_B.err()<<endl;
    cout<<"Al_B: "<<Al_B.ave()<<" +- "<<Al_B.err()<<endl;
    cout<<"afp_A: "<<afp_A.ave()<<" +- "<<afp_A.err()<<endl;
    cout<<"Am_A: "<<Am_A.ave()<<" +- "<<Am_A.err()<<endl;
    cout<<"A2m_A: "<<A2m_A.ave()<<" +- "<<A2m_A.err()<<endl;
    cout<<"Al_A: "<<Al_A.ave()<<" +- "<<Al_A.err()<<endl;
    cout<<"Al: "<<Al.ave()<<" +- "<<Al.err()<<endl;

    cout<<"reduced ch2: "<<ch2/(Nens_tot - bf.Get_number_of_fit_pars())<<endl;





    //print fitting function


    //start with continuum limit as a function of csi
    int num_csi_points= 100;
    Vfloat csi_p_list;
    Vfloat extr_A_val, extr_A_err;
    Vfloat extr_B_val, extr_B_err;

    
      
      
    for(int icsi=0; icsi<num_csi_points;icsi++) {

      double csi_p = icsi*10.0*(csi_phys)/num_csi_points;
      csi_p_list.push_back(csi_p);
      
      distr_t extr_A(UseJack);
      distr_t extr_B(UseJack);
      
      for(int ijack=0; ijack<Njacks;ijack++) {
	
	sc_fp_fpar Xi;
	sc_fp_ipar Yi_A;
	sc_fp_ipar Yi_B;

	//populate fit parameters
	Xi.afp_B = afp_B.distr[ijack];
	Xi.Am_B = Am_B.distr[ijack];
	Xi.A2m_B = A2m_B.distr[ijack];
	Xi.afp_A = afp_A.distr[ijack];
	Xi.Am_A = Am_A.distr[ijack];
	Xi.A2m_A = A2m_A.distr[ijack];
	Xi.Al = 0.0; //we look at thermodynamic limit here
	Xi.Al_B = 0.0;
	Xi.Al_A = 0.0;
	
	//populate input parameters
	Yi_B.L = 10.0; //will not be used
	Yi_B.csi = csi_p;
	Yi_B.csi_ph = csi_phys_distr.distr[ijack];
	Yi_B.Is_B = true;
	Yi_B.Mp = 1.0;

	Yi_A.L = 10.0; //will not be used
	Yi_A.csi = csi_p;
	Yi_A.csi_ph = csi_phys_distr.distr[ijack];
	Yi_A.Mp = 1.0;
	Yi_A.Is_B = false;
	  
	extr_B.distr.push_back( bf.ansatz(Xi,Yi_B));
        extr_A.distr.push_back( bf.ansatz(Xi, Yi_A));
	


      }
	
      extr_B_val.push_back( extr_B.ave());
      extr_B_err.push_back( extr_B.err());
      extr_A_val.push_back( extr_A.ave());
      extr_A_err.push_back( extr_A.err());
      
    }

   

    //print extr_val and error for each lattice spacing
    

    Print_To_File({}, {csi_p_list, extr_B_val, extr_B_err}, "../data/scale_setting/fit_B_afp.dat", "", "# csi val err");
    Print_To_File({}, {csi_p_list, extr_A_val, extr_A_err}, "../data/scale_setting/fit_A_afp.dat", "", "# csi val err");
    
    


    //correct lattice points by FSEs
    distr_t_list afp_corr_by_FSEs_B(UseJack);
    distr_t_list afp_uncorr_B(UseJack);
    distr_t_list afp_corr_by_FSEs_A(UseJack);
    distr_t_list afp_uncorr_A(UseJack);
    
    for(int iens=0;iens<Nens_B;iens++) {
      distr_t MpL_distr = L_B_ens[iens]*Mpi_CDH_B.distr_list[iens];
      afp_corr_by_FSEs_B.distr_list.push_back( fp_CDH_B.distr_list[iens]/(1.0 + (Al+ Al_B)*Csi_CDH_B.distr_list[iens]*distr_t::f_of_distr(FVE, MpL_distr)));
      afp_uncorr_B.distr_list.push_back( fp_CDH_B.distr_list[iens]); 
    }

    for(int iens=0;iens<Nens_A;iens++) {
      distr_t MpL_distr = L_A_ens[iens]*Mpi_CDH_A.distr_list[iens];
      afp_corr_by_FSEs_A.distr_list.push_back( fp_CDH_A.distr_list[iens]/(1.0 + (Al+ Al_A)*Csi_CDH_A.distr_list[iens]*distr_t::f_of_distr(FVE, MpL_distr)));
      afp_uncorr_A.distr_list.push_back( fp_CDH_A.distr_list[iens]); 
    }
    
    //Print corrected lattice points to file

    //B ensemble
    Print_To_File({Ensemble_B_tag_list},{afp_corr_by_FSEs_B.ave(), afp_corr_by_FSEs_B.err(), afp_uncorr_B.ave(), afp_uncorr_B.err(), Csi_CDH_B.ave(), Csi_CDH_B.err()},  "../data/scale_setting/data_FSEs_corrected_afp_B.dat", "", "# Ens   corrected   uncorrected    xi ");

    //A ensemble
    Print_To_File({Ensemble_A_tag_list},{afp_corr_by_FSEs_A.ave(), afp_corr_by_FSEs_A.err(), afp_uncorr_A.ave(), afp_uncorr_A.err(), Csi_CDH_A.ave(), Csi_CDH_A.err()},  "../data/scale_setting/data_FSEs_corrected_afp_A.dat", "", "# Ens   corrected   uncorrected     xi");
      


    //get lattice spacings

    //get id of C and D ensemble
    int C_id, D_id, E_id;
    bool find_C_id=false;
    bool find_D_id=false;
    bool find_E_id=false;
    for(int iens=0;iens<Nens_phys;iens++) {
      if(Ensemble_phys_point_tag_list[iens].substr(1,1) == "C") { find_C_id = true; C_id =iens;}
      if(Ensemble_phys_point_tag_list[iens].substr(1,1) == "D") { find_D_id = true; D_id =iens;}
      if(Ensemble_phys_point_tag_list[iens].substr(1,1) == "E") { find_E_id = true; E_id =iens;}
    }

    if(!find_C_id || !find_D_id || !find_E_id) crash("Cannot find phyisical point ensemble C,D,E in Ensemble_phys_point_tag_list");
      

 
    //compute Am and Am_art coefficients
    distr_t Am_cont, Am_art, A2m_cont, A2m_art;
    if(Use_three_finest_in_scale_setting_fp) {Am_cont = Am_B; Am_art = 0.0*Am_B; A2m_cont = A2m_B; A2m_art = 0.0*A2m_B;}
    else {
      Am_cont = ( Am_B*(afp_A*afp_A) - Am_A*(afp_B*afp_B) )/( afp_A*afp_A - afp_B*afp_B);
      Am_art = (Am_B - Am_A)/( afp_B*afp_B  - afp_A*afp_A);
      if( fabs(A2m_A.ave()) < 1e-16) { A2m_cont = A2m_B; A2m_art = 0.0*A2m_B;}
      else {
      A2m_cont = ( A2m_B*(afp_A*afp_A) - A2m_A*(afp_B*afp_B) )/( afp_A*afp_A - afp_B*afp_B);
      A2m_art = (A2m_B - A2m_A)/( afp_B*afp_B  - afp_A*afp_A);
      }
    }

    //correct for xi_p mistuning on C and D ensembles

    distr_t corr_C_mass, corr_D_mass, corr_E_mass;

    corr_C_mass = 1.0  -2.0*Csi_CDH_phys_point[C_id]*distr_t::f_of_distr(LOG, Csi_CDH_phys_point[C_id]/csi_phys_distr)  + (Am_cont+ Am_art*fp_CDH_phys_point[C_id]*fp_CDH_phys_point[C_id])*(Csi_CDH_phys_point[C_id] - csi_phys_distr) +  (A2m_cont + A2m_art*fp_CDH_phys_point[C_id]*fp_CDH_phys_point[C_id])*(Csi_CDH_phys_point[C_id]*Csi_CDH_phys_point[C_id] - csi_phys_distr*csi_phys_distr);

    corr_D_mass = 1.0  -2.0*Csi_CDH_phys_point[D_id]*distr_t::f_of_distr(LOG, Csi_CDH_phys_point[D_id]/csi_phys_distr)  + (Am_cont  + Am_art*fp_CDH_phys_point[D_id]*fp_CDH_phys_point[D_id])*(Csi_CDH_phys_point[D_id] - csi_phys_distr) +  (A2m_cont + A2m_art*fp_CDH_phys_point[D_id]*fp_CDH_phys_point[D_id])*(Csi_CDH_phys_point[D_id]*Csi_CDH_phys_point[D_id] - csi_phys_distr*csi_phys_distr);

    corr_E_mass = 1.0  -2.0*Csi_CDH_phys_point[E_id]*distr_t::f_of_distr(LOG, Csi_CDH_phys_point[E_id]/csi_phys_distr)  + (Am_cont  + Am_art*fp_CDH_phys_point[E_id]*fp_CDH_phys_point[E_id])*(Csi_CDH_phys_point[E_id] - csi_phys_distr) +  (A2m_cont + A2m_art*fp_CDH_phys_point[E_id]*fp_CDH_phys_point[E_id])*(Csi_CDH_phys_point[E_id]*Csi_CDH_phys_point[E_id] - csi_phys_distr*csi_phys_distr);


    //correct for FSEs on C and D ensembles
    
    //get MpiL for C and D ensembles

    distr_t MpL_distr_C = L_phys_point[C_id]*Mpi_CDH_phys_point[C_id];
    distr_t MpL_distr_D = L_phys_point[D_id]*Mpi_CDH_phys_point[D_id];
    distr_t MpL_distr_E = L_phys_point[E_id]*Mpi_CDH_phys_point[E_id];

    
    distr_t corr_C_FSEs, corr_D_FSEs, corr_E_FSEs;

    corr_C_FSEs = (1.0 + (Al+ Al_B)*distr_t::f_of_distr(pow_2,Csi_CDH_phys_point.distr_list[C_id])*distr_t::f_of_distr(FVE, MpL_distr_C));
    corr_D_FSEs = (1.0 + (Al+ Al_B)*distr_t::f_of_distr(pow_2,Csi_CDH_phys_point.distr_list[D_id])*distr_t::f_of_distr(FVE, MpL_distr_D));
    corr_E_FSEs = (1.0 + (Al+ Al_B)*distr_t::f_of_distr(pow_2,Csi_CDH_phys_point.distr_list[E_id])*distr_t::f_of_distr(FVE, MpL_distr_E));
    

    
    

    a_distr_A = afp_A/fpi_phys_distr;
    
    a_distr_B = afp_B/fpi_phys_distr;

    a_distr_C = fp_CDH_phys_point[C_id]/(corr_C_FSEs*corr_C_mass)/fpi_phys_distr;

    a_distr_D = fp_CDH_phys_point[D_id]/(corr_D_FSEs*corr_D_mass)/fpi_phys_distr;

    a_distr_E = fp_CDH_phys_point[E_id]/(corr_E_FSEs*corr_E_mass)/fpi_phys_distr;



    //print on screen
    //print afC and afD before and after corrections:
    cout<<"C ensemble"<<endl;
    cout<<"af_C: "<<fpi_phys_point.distr_list[C_id].ave()<<" +- "<<fpi_phys_point.distr_list[C_id].err()<<endl;
    cout<<"af_C (CDH-corrected): "<<fp_CDH_phys_point[C_id].ave()<<" +- "<<fp_CDH_phys_point[C_id].err()<<endl;
    cout<<"af_C (infL): "<<(fp_CDH_phys_point[C_id]/(corr_C_FSEs)).ave()<<" +- "<<(fp_CDH_phys_point[C_id]/corr_C_FSEs).err()<<endl;
    cout<<"af_C (infL-mass corrected): "<<(fp_CDH_phys_point[C_id]/(corr_C_FSEs*corr_C_mass)).ave()<<" +- "<<(fp_CDH_phys_point[C_id]/(corr_C_FSEs*corr_C_mass)).err()<<endl;
    cout<<"D ensemble"<<endl;
    cout<<"af_D: "<<fpi_phys_point.distr_list[D_id].ave()<<" +- "<<fpi_phys_point.distr_list[D_id].err()<<endl;
    cout<<"af_D (CDH-corrected): "<<fp_CDH_phys_point[D_id].ave()<<" +- "<<fp_CDH_phys_point[D_id].err()<<endl;
    cout<<"af_D (infL): "<<(fp_CDH_phys_point[D_id]/(corr_D_FSEs)).ave()<<" +- "<<(fp_CDH_phys_point[D_id]/corr_D_FSEs).err()<<endl;
    cout<<"af_D (infL-mass corrected): "<<(fp_CDH_phys_point[D_id]/(corr_D_FSEs*corr_D_mass)).ave()<<" +- "<<(fp_CDH_phys_point[D_id]/(corr_D_FSEs*corr_D_mass)).err()<<endl;

    cout<<"E ensemble"<<endl;
    cout<<"af_E: "<<fpi_phys_point.distr_list[E_id].ave()<<" +- "<<fpi_phys_point.distr_list[E_id].err()<<endl;
    cout<<"af_E (CDH-corrected): "<<fp_CDH_phys_point[E_id].ave()<<" +- "<<fp_CDH_phys_point[E_id].err()<<endl;
    cout<<"af_E (infL): "<<(fp_CDH_phys_point[E_id]/(corr_E_FSEs)).ave()<<" +- "<<(fp_CDH_phys_point[E_id]/corr_E_FSEs).err()<<endl;
    cout<<"af_E (infL-mass corrected): "<<(fp_CDH_phys_point[E_id]/(corr_E_FSEs*corr_E_mass)).ave()<<" +- "<<(fp_CDH_phys_point[E_id]/(corr_E_FSEs*corr_E_mass)).err()<<endl;   
    
    cout<<"#######  Printing lattice spacing determined from afp analysis ########"<<endl;
    cout<<"A ensemble: "<< a_distr_A.ave()<<" +- "<<a_distr_A.err()<<" [GeV-1]    "<<(a_distr_A/fm_to_iGev).ave()<<" +- "<<(a_distr_A/fm_to_iGev).err()<<" [fm] "<<endl;
    cout<<"B ensemble: "<< a_distr_B.ave()<<" +- "<<a_distr_B.err()<<" [GeV-1]    "<<(a_distr_B/fm_to_iGev).ave()<<" +- "<<(a_distr_B/fm_to_iGev).err()<<" [fm] "<<endl;
    cout<<"C ensemble: "<< a_distr_C.ave()<<" +- "<<a_distr_C.err()<<" [GeV-1]    "<<(a_distr_C/fm_to_iGev).ave()<<" +- "<<(a_distr_C/fm_to_iGev).err()<<" [fm] "<<endl;
    cout<<"D ensemble: "<< a_distr_D.ave()<<" +- "<<a_distr_D.err()<<" [GeV-1]    "<<(a_distr_D/fm_to_iGev).ave()<<" +- "<<(a_distr_D/fm_to_iGev).err()<<" [fm] "<<endl;
    cout<<"E ensemble: "<< a_distr_E.ave()<<" +- "<<a_distr_E.err()<<" [GeV-1]    "<<(a_distr_E/fm_to_iGev).ave()<<" +- "<<(a_distr_E/fm_to_iGev).err()<<" [fm] "<<endl;


 

    

  return;
}
