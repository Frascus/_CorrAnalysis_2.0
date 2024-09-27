#include "../include/RC_analysis.h"
#include "input.h"
#include "numerics.h"


const double alpha = 1.0/137.04;
const bool UseJack=1;
const int Njacks=30;
const int Nboots=800;
const double fm_to_inv_Gev = 1.0 / 0.197327;
const int Nconfs=30;
using namespace std;

void Perform_RC_analysis() {

  
  vector<tuple<int,int,int,int>> p_list_A({ {1,2,2,2}, {2,2,2,2}, {3,2,2,2}, {4,2,2,2}, {5,2,2,2}, {6,2,2,2}, {1,3,3,3}, {2,3,3,3}, {3,3,3,3}, {4,3,3,3}, {5,3,3,3}, {6,3,3,3}, {7,3,3,3}, {8,3,3,3}, {2,4,4,4}, {9,3,3,3}, {3,4,4,4}, {4,4,4,4}, {10,3,3,3}, {5,4,4,4}, {6,4,4,4}, {7,4,4,4}, {8,4,4,4}, {9,4,4,4}, {10,4,4,4}, {2,5,5,5}, {3,5,5,5}, {4,5,5,5}, {11,4,4,4}, {5,5,5,5}, {6,5,5,5}, {12,4,4,4}, {7,5,5,5}, {8,5,5,5}, {13,4,4,4}, {9,5,5,5}, {14,4,4,4}, {10,5,5,5}, {11,5,5,5}, {3,6,6,6}, {4,6,6,6}, {12,5,5,5}, {5,6,6,6} });

  vector<tuple<int,int,int,int>> p_list_B({ {1,2,2,2},{2,2,2,2},{3,2,2,2},{4,2,2,2},{5,2,2,2},{6,2,2,2},{1,3,3,3},{2,3,3,3},{3,3,3,3},{4,3,3,3},{5,3,3,3},{6,3,3,3},{7,3,3,3},{8,3,3,3},{2,4,4,4},{9,3,3,3},{3,4,4,4},{4,4,4,4},{10,3,3,3},{5,4,4,4},{6,4,4,4},{7,4,4,4},{8,4,4,4},{9,4,4,4},{10,4,4,4},{2,5,5,5},{3,5,5,5},{4,5,5,5},{11,4,4,4},{5,5,5,5},{6,5,5,5},{12,4,4,4},{7,5,5,5},{8,5,5,5},{13,4,4,4},{9,5,5,5},{14,4,4,4},{10,5,5,5}});
     
  vector<tuple<int,int,int,int>> p_list_C({ {2,2,2,2},{3,2,2,2},{4,2,2,2},{5,2,2,2},{3,3,3,3},{4,3,3,3},{5,3,3,3},{6,3,3,3},{7,3,3,3},{8,3,3,3},{3,4,4,4},{4,4,4,4},{5,4,4,4},{6,4,4,4},{7,4,4,4},{8,4,4,4},{9,4,4,4},{10,4,4,4},{4,5,5,5},{11,4,4,4},{5,5,5,5},{6,5,5,5},{12,4,4,4},{7,5,5,5},{8,5,5,5},{9,5,5,5},{10,5,5,5},{11,5,5,5},{12,5,5,5},{5,6,6,6},{6,6,6,6},{13,5,5,5},{7,6,6,6},{8,6,6,6},{14,5,5,5},{9,6,6,6},{15,5,5,5},{10,6,6,6},{11,6,6,6},{16,5,5,5},{12,6,6,6} });

  vector<tuple<int,int,int,int>> p_list_D({ {6,3,3,3},{8,3,3,3},{9,3,3,3},{10,3,3,3},{5,4,4,4},{6,4,4,4},{7,4,4,4},{8,4,4,4},{9,4,4,4},{10,4,4,4},{11,4,4,4},{5,5,5,5},{6,5,5,5},{7,5,5,5},{8,5,5,5},{9,5,5,5},{10,5,5,5},{11,5,5,5},{6,6,6,6},{7,6,6,6},{8,6,6,6},{9,6,6,6},{10,6,6,6},{11,6,6,6},{12,6,6,6},{13,6,6,6},{6,7,7,7},{14,6,6,6},{7,7,7,7},{8,7,7,7},{15,6,6,6},{9,7,7,7},{10,7,7,7},{16,6,6,6},{7,8,8,8},{9,8,8,8},{10,8,8,8},{11,8,8,8},{13,8,8,8}});

  

  Vfloat musea_A_list({   0.006, 0.008, 0.01, 0.0115      });
  Vfloat muval_A_list({   0.006, 0.008, 0.01, 0.0115, 0.013, 0.015, 0.017, 0.019, 0.021  });

  Vfloat musea_B_list({0.006, 0.0075, 0.0088, 0.01});
  Vfloat muval_B_list({0.005, 0.006, 0.0075, 0.009, 0.01, 0.011, 0.013, 0.015, 0.017});

  Vfloat musea_C_list({0.005, 0.0065, 0.008, 0.0095});
  Vfloat muval_C_list({0.004, 0.005, 0.0065, 0.008, 0.0095, 0.011, 0.0125, 0.014, 0.0155});

  Vfloat musea_D_list({0.004, 0.005, 0.0065, 0.008});
  Vfloat muval_D_list({0.0035, 0.004, 0.005, 0.0065, 0.008, 0.009, 0.01, 0.011, 0.012});

  VVfloat musea_list({ musea_A_list, musea_B_list, musea_C_list, musea_D_list});
  VVfloat muval_list({ muval_A_list, muval_B_list, muval_C_list, muval_D_list});

  vector<vector<tuple<int,int,int,int>>> p_list({ p_list_A, p_list_B, p_list_C, p_list_D});


  boost::filesystem::create_directory("../data/RC_analysis");
  


  vector<string> Ensembles({"A24", "B24", "C32", "D48"});
  vector<string> RCs({"Va", "Vp", "Vs", "Vt", "Vv", "Zq"});

  for(auto &RC: RCs) {
    boost::filesystem::create_directory("../data/RC_analysis/"+RC);
    boost::filesystem::create_directory("../data/RC_analysis/"+RC+"/mu_val_extr");
    boost::filesystem::create_directory("../data/RC_analysis/"+RC+"/mu_sea_extr");
    for(auto &Ens :Ensembles) {
      boost::filesystem::create_directory("../data/RC_analysis/"+RC+"/mu_sea_extr/"+Ens);
      boost::filesystem::create_directory("../data/RC_analysis/"+RC+"/mu_val_extr/"+Ens);
      boost::filesystem::create_directory("../data/RC_analysis/"+RC+"/mu_sea_extr/"+Ens+"/fit_func");
      boost::filesystem::create_directory("../data/RC_analysis/"+RC+"/mu_val_extr/"+Ens+"/fit_func");
    }
  }

  
    

  for (int iens=0; iens<(signed)Ensembles.size();iens++) {

    string Ens=Ensembles[iens];
    vector<distr_t_list> Z_extr(RCs.size()-1);
    vector<vector<distr_t_list>> Z_list(RCs.size()-1);
    vector<vector<distr_t_list>> Z_list_rev(RCs.size()-1);
    for(auto & Z: Z_list) {
      for(int isea=0;isea<(signed)musea_list[iens].size(); isea++) Z.emplace_back(UseJack);
     
    }
    for(auto & Z: Z_list_rev) {
      for(int imom=0;imom<(signed)p_list[iens].size(); imom++)  Z.emplace_back(UseJack);
    }
    
    for(int isea=0;isea<(signed)musea_list[iens].size();isea++) {

      cout<<"Analyzing musea: "<< musea_list[iens][isea]<<endl;

      cout<<"looping over momenta"<<endl;

      vector<distr_t_list> Vert_musea_extr;
      
      for(int irc=0;irc<(signed)RCs.size();irc++) Vert_musea_extr.emplace_back(UseJack);

      for(int imom=0;imom<(signed)p_list[iens].size();imom++) {

	  
	cout<<"Looping over RCs"<<endl;
	vector<distr_t_list> Vert_musea;
	for(int irc=0;irc<(signed)RCs.size();irc++) Vert_musea.emplace_back(UseJack);
	//loop over Rcs
	for(int irc=0; irc<(signed)RCs.size();irc++) {
	  
	  cout<<"Analyzing RC: "<<RCs[irc]<<endl;

	  cout<<"looping over muval: "<<endl;
	  
	  for(int ival=0;ival<(signed)muval_list[iens].size();ival++) {
	    
	    cout<<"Analyzing muval: "<<muval_list[iens][ival]<<endl;
	    //Read files
	    Vert_musea[irc].distr_list.emplace_back( UseJack, Read_From_File( "../RC_data/data/"+Ens+"/"+RCs[irc]+"_p_"+to_string(imom)+"_musea_"+to_string(isea)+"_muval_"+to_string(ival)+".dat", 0,1));
	  }

	
	  
	  //performing chiral-valence extrapolation
	  
	  //include in the fit all data satisfying mu_val >= 0.75musea
	  int Nmeas=0;
	  distr_t_list Vert_musea_fitted(UseJack);
	  Vfloat muval_list_fitted;
	  
	  for(int ival=0;ival<(signed)muval_list[iens].size();ival++) {
	    if( (muval_list[iens][ival] >= 0.75*musea_list[iens][isea]) || RCs[irc] == "Vp") {
	      Nmeas++;
	      muval_list_fitted.push_back( muval_list[iens][ival]);
	      Vert_musea_fitted.distr_list.push_back( Vert_musea[irc][ival]);
	    }
	  }
	    
	  //print data to file
	  Print_To_File({}, {muval_list[iens], Vert_musea[irc].ave(), Vert_musea[irc].err()},  "../data/RC_analysis/"+RCs[irc]+"/mu_val_extr/"+Ens+"/p_"+to_string(imom)+"_musea_"+to_string(isea)+".list", "", "");
	  Print_To_File({}, {muval_list_fitted, Vert_musea_fitted.ave(), Vert_musea_fitted.err()},  "../data/RC_analysis/"+RCs[irc]+"/mu_val_extr/"+Ens+"/p_"+to_string(imom)+"_musea_"+to_string(isea)+".fit_list", "", "");
	  
	  //define par and ipar classes for the chiral-valence fits
	  class ipar_mval {
	  public:
	    ipar_mval() : FF(0.0), FF_err(0.0) {}
	    double FF, FF_err, mu;
	  };
	  
	  class fpar_mval {
	  public:
	    fpar_mval() {}
	    fpar_mval(const Vfloat &par) {
	      if((signed)par.size() != 4) crash("In class fpar_mval  class constructor Vfloat par has size != 4");
	      R=par[0];
	      A=par[1];
	      B=par[2];
	      C=par[3];
	    }
	    double R,A,B,C;
	  };
	  
	  //init bootstrap fit
	  bootstrap_fit<fpar_mval,ipar_mval> bf_mval(UseJack?Njacks:800);
	  bf_mval.set_warmup_lev(1); //sets warmup
	  bf_mval.Set_number_of_measurements(Nmeas);
	  bf_mval.Set_verbosity(1);
	  
	  
	  bf_mval.Add_par("R", 1.0, 0.1);
	  bf_mval.Add_par("A", 1.0, 0.1);
	  bf_mval.Add_par("B", 1.0 , 0.1);
	  bf_mval.Add_par("C", 1.0 , 0.1);
	  
	  //fit on mean values to get ch2
	  bootstrap_fit<fpar_mval,ipar_mval> bf_mval_ch2(1);
	  bf_mval_ch2.set_warmup_lev(1); //sets warmup
	  bf_mval_ch2.Set_number_of_measurements(Nmeas);
	  bf_mval_ch2.Set_verbosity(1);
	  bf_mval_ch2.Add_par("R", 1.0, 0.1);
	  bf_mval_ch2.Add_par("A", 1.0, 0.1);
	  bf_mval_ch2.Add_par("B", 1.0, 0.1);
	  bf_mval_ch2.Add_par("C", 1.0, 0.1);
	  if(RCs[irc] != "Vp") { bf_mval.Fix_par("C",0.0); bf_mval_ch2.Fix_par("C",0.0);}
	  
	  //ansatz
	  bf_mval.ansatz=  [](const fpar_mval &p, const ipar_mval &ip) {
	    return  p.R  +  p.A*ip.mu + p.B*pow(ip.mu,2) + p.C/ip.mu;
	  };
	  bf_mval.measurement=  [ ](const fpar_mval &p, const ipar_mval &ip) {
	    return ip.FF;
	     };
	  bf_mval.error=  [ ](const fpar_mval &p, const ipar_mval &ip) {
	    return ip.FF_err;
	  };

	  bf_mval_ch2.ansatz= bf_mval.ansatz;
	  bf_mval_ch2.measurement = bf_mval.measurement;
	  bf_mval_ch2.error = bf_mval.error;

	  //start fitting
	  //fill the data
	  vector<vector<ipar_mval>> data_mval(UseJack?Njacks:800);
	  vector<vector<ipar_mval>> data_mval_ch2(1);
	  //allocate space for output result
	  boot_fit_data<fpar_mval> Bt_fit_mval;
	  boot_fit_data<fpar_mval> Bt_fit_mval_ch2;
	  boot_fit_data<fpar_mval> Bt_fit_mval2;
	  boot_fit_data<fpar_mval> Bt_fit_mval2_ch2;
	  for(auto &data_iboot: data_mval) data_iboot.resize(Nmeas);
	  for(auto &data_iboot: data_mval_ch2) data_iboot.resize(Nmeas);
	  for(int ijack=0;ijack<(UseJack?Njacks:800);ijack++) {
	    int im=0;
	    for(int ival=0;ival<(signed)muval_list[iens].size();ival++) {
	      if( (muval_list[iens][ival] >= 0.75*musea_list[iens][isea]) || RCs[irc] == "Vp") {
		data_mval[ijack][im].FF = (Vert_musea[irc].distr_list[ival]).distr[ijack];
		data_mval[ijack][im].FF_err= (Vert_musea[irc].distr_list[ival]).err();
		data_mval[ijack][im].mu= muval_list[iens][ival];
		
		
		if(ijack==0) {
		  data_mval_ch2[ijack][im].FF = (Vert_musea[irc].distr_list[ival]).ave();
		  data_mval_ch2[ijack][im].FF_err= (Vert_musea[irc].distr_list[ival]).err();
		  data_mval_ch2[ijack][im].mu= muval_list[iens][ival];
		}
		im++;
	      }
	    }
	  }

	  //append
	  bf_mval.Append_to_input_par(data_mval);
	  bf_mval_ch2.Append_to_input_par(data_mval_ch2);
	  //fit
	  cout<<"Fitting RC: "<<RCs[irc]<<endl;
	  Bt_fit_mval2= bf_mval.Perform_bootstrap_fit();
	  Bt_fit_mval2_ch2= bf_mval_ch2.Perform_bootstrap_fit();
	  //fix muval^2 parameter and refit
	  bf_mval.Fix_par("B", 0.0);
	  bf_mval_ch2.Fix_par("B",0.0);
	  Bt_fit_mval= bf_mval.Perform_bootstrap_fit();
	  Bt_fit_mval_ch2= bf_mval_ch2.Perform_bootstrap_fit();
	  
	  double ch2_red_mval= Bt_fit_mval_ch2.get_ch2_ave()/( Nmeas -2.0);
	  double ch2_red_mval2=  Bt_fit_mval2_ch2.get_ch2_ave()/( Nmeas -3.0);
	  
	  //retrive params
	  distr_t R(UseJack), A(UseJack), B(UseJack), C(UseJack);
	  distr_t R2(UseJack), A2(UseJack), B2(UseJack), C2(UseJack);
	  for(int ijack=0; ijack< (UseJack?Njacks:800) ;ijack++) {
	    R.distr.push_back(   Bt_fit_mval.par[ijack].R );
	    A.distr.push_back(   Bt_fit_mval.par[ijack].A );
	    B.distr.push_back(   Bt_fit_mval.par[ijack].B );
	    C.distr.push_back(   Bt_fit_mval.par[ijack].C );
	    
	    R2.distr.push_back(   Bt_fit_mval2.par[ijack].R );
	    A2.distr.push_back(   Bt_fit_mval2.par[ijack].A );
	    B2.distr.push_back(   Bt_fit_mval2.par[ijack].B );
	    C2.distr.push_back(   Bt_fit_mval2.par[ijack].C );
	    
	  }
	  
	  //propagate systematic error
	  distr_t RES= R.ave() + (R-R.ave())*(sqrt( pow(R.err(),2) + pow( R.ave() - R2.ave(),2)))/R.err();
	  
	  Vert_musea_extr[irc].distr_list.push_back(RES);
	  
	  //print fit function
	  int Nprint=300;
	  distr_t_list V_fit(UseJack), V2_fit(UseJack);
	  Vfloat muval_to_print(Nprint);
	  for(int i=0;i<Nprint;i++) muval_to_print[i] = i*1.5*muval_list[iens][muval_list[iens].size()-1]/(Nprint -1.0);
	  for(int mm=0;mm<Nprint;mm++) {
	    double mu= muval_to_print[mm];
	    distr_t V_distr(UseJack), V2_distr(UseJack);
	    for(int ijack=0;ijack<(UseJack?Njacks:800);ijack++) {
	      ipar_mval IP;
	      IP.mu = mu;
	      V_distr.distr.push_back( bf_mval.ansatz(Bt_fit_mval.par[ijack], IP));
	      V2_distr.distr.push_back( bf_mval.ansatz(Bt_fit_mval2.par[ijack], IP));
	    }
	    V_fit.distr_list.push_back( V_distr);
	    V2_fit.distr_list.push_back( V2_distr);
	  }
	  
	  //print fit function and results
	  Print_To_File({}, {muval_to_print, V_fit.ave(), V_fit.err(), V2_fit.ave(), V2_fit.err()}, "../data/RC_analysis/"+RCs[irc]+"/mu_val_extr/"+Ens+"/fit_func/p_"+to_string(imom)+"_musea_"+to_string(isea)+".func","" , "#ch2/dof : "+to_string_with_precision(ch2_red_mval,5)+" [ A*mu ], "+to_string_with_precision(ch2_red_mval2,5)+" [A*mu + B*mu^2], Nmeas: "+to_string(Nmeas)) ;
	  
	}
	
      }
    
      
      for(int irc=0;irc<(signed)RCs.size();irc++) Print_To_File({}, {Vert_musea_extr[irc].ave(), Vert_musea_extr[irc].err()}, "../data/RC_analysis/"+RCs[irc]+"/mu_val_extr/"+Ens+"/musea_extr_"+to_string(isea)+".plist", "", "");

      //push_back result for Z
      for(int irc=0;irc<(signed)RCs.size()-1; irc++) {
	for(int imom=0;imom<(signed)p_list[iens].size(); imom++) {
	  Z_list[irc][isea].distr_list.push_back(  Vert_musea_extr[RCs.size()-1][imom]/Vert_musea_extr[irc][imom] );
	  Z_list_rev[irc][imom].distr_list.push_back( Vert_musea_extr[RCs.size()-1][imom]/Vert_musea_extr[irc][imom] );
	}
      }
    }

      


    //perform extrapolation in musea for each momentum

    //define par and ipar classes for the chiral-valence fits
    class ipar_msea {
    public:
      ipar_msea() : FF(0.0), FF_err(0.0) {}
      double FF, FF_err, mu;
    };
    
    class fpar_msea {
    public:
      fpar_msea() {}
      fpar_msea(const Vfloat &par) {
	if((signed)par.size() != 2) crash("In class fpar_msea  class constructor Vfloat par has size != 2");
	R=par[0];
	A=par[1];
      }
      double R,A;
    };
    
    
    
    //loop over RCs
    for(int irc=0; irc < (signed)RCs.size() -1 ; irc++) {
      for(int imom=0; imom < (signed)p_list[iens].size() ; imom ++) {
      
	int Nmeas= musea_list[iens].size();
	
	//init bootstrap fit
	bootstrap_fit<fpar_msea,ipar_msea> bf_msea(UseJack?Njacks:800);
	bf_msea.set_warmup_lev(1); //sets warmup
	bf_msea.Set_number_of_measurements(Nmeas);
	bf_msea.Set_verbosity(1);
	
      
	bf_msea.Add_par("R", 1.0, 0.1);
	bf_msea.Add_par("A", 1.0, 0.1);
      
	//fit on mean values to get ch2
	bootstrap_fit<fpar_msea,ipar_msea> bf_msea_ch2(1);
	bf_msea_ch2.set_warmup_lev(1); //sets warmup
	bf_msea_ch2.Set_number_of_measurements(Nmeas);
	bf_msea_ch2.Set_verbosity(1);
	bf_msea_ch2.Add_par("R", 1.0, 0.1);
	bf_msea_ch2.Add_par("A", 1.0, 0.1);
	
	//ansatz
	bf_msea.ansatz=  [](const fpar_msea &p, const ipar_msea &ip) {
	  return  p.R  +  p.A*ip.mu ;
	};
	bf_msea.measurement=  [ ](const fpar_msea &p, const ipar_msea &ip) {
	  return ip.FF;
	};
	bf_msea.error=  [ ](const fpar_msea &p, const ipar_msea &ip) {
	  return ip.FF_err;
	};
	
	bf_msea_ch2.ansatz= bf_msea.ansatz;
	bf_msea_ch2.measurement = bf_msea.measurement;
	bf_msea_ch2.error = bf_msea.error;


	//start fitting
	//fill the data
	vector<vector<ipar_msea>> data_msea(UseJack?Njacks:800);
	vector<vector<ipar_msea>> data_msea_ch2(1);
	//allocate space for output result
	boot_fit_data<fpar_msea> Bt_fit_msea;
	boot_fit_data<fpar_msea> Bt_fit_msea_ch2;
	boot_fit_data<fpar_msea> Bt_fit_msea_const;
	boot_fit_data<fpar_msea> Bt_fit_msea_const_ch2;
	  for(auto &data_iboot: data_msea) data_iboot.resize(Nmeas);
	  for(auto &data_iboot: data_msea_ch2) data_iboot.resize(Nmeas);
	  for(int ijack=0;ijack<(UseJack?Njacks:800);ijack++) {
	    for(int isea=0;isea<(signed)musea_list[iens].size();isea++) {
	      data_msea[ijack][isea].FF = (Z_list[irc][isea][imom]).distr[ijack];
	      data_msea[ijack][isea].FF_err= (Z_list[irc][isea][imom]).err();
	      data_msea[ijack][isea].mu= musea_list[iens][isea];
	      
	      if(ijack==0) {
		data_msea_ch2[ijack][isea].FF = (Z_list[irc][isea][imom]).ave();
		data_msea_ch2[ijack][isea].FF_err= (Z_list[irc][isea][imom]).err();
		data_msea_ch2[ijack][isea].mu= musea_list[iens][isea];
	      }
	    }
	  }
	  
	  //append
	  bf_msea.Append_to_input_par(data_msea);
	  bf_msea_ch2.Append_to_input_par(data_msea_ch2);
	  //fit
	  cout<<"Fitting RC: "<<RCs[irc]<<endl;
	  Bt_fit_msea= bf_msea.Perform_bootstrap_fit();
	  Bt_fit_msea_ch2= bf_msea_ch2.Perform_bootstrap_fit();
	  //fix muval^2 parameter and refit
	  bf_msea.Fix_par("A", 0.0);
	  bf_msea_ch2.Fix_par("A",0.0);
	  Bt_fit_msea_const= bf_msea.Perform_bootstrap_fit();
	  Bt_fit_msea_const_ch2= bf_msea_ch2.Perform_bootstrap_fit();

	  	  
	  double ch2_red_msea= Bt_fit_msea_ch2.get_ch2_ave()/( Nmeas -2.0);
	  double ch2_red_msea_const=  Bt_fit_msea_const_ch2.get_ch2_ave()/( Nmeas -1.0);
	  int Ndof = Nmeas-2;
	  int Ndof_const= Nmeas-1;
	  
	  //retrive params
	  distr_t R(UseJack), A(UseJack);
	  distr_t R_const(UseJack), A_const(UseJack);
	  for(int ijack=0; ijack< (UseJack?Njacks:800) ;ijack++) {
	    R.distr.push_back(   Bt_fit_msea.par[ijack].R );
	    A.distr.push_back(   Bt_fit_msea.par[ijack].A );
	   	    
	    R_const.distr.push_back(   Bt_fit_msea_const.par[ijack].R );
	    A_const.distr.push_back(   Bt_fit_msea_const.par[ijack].A );
	  }



	  //propagate systematic error
	  distr_t RES= AIC( {R, R_const}, { ch2_red_msea*Ndof, ch2_red_msea_const*Ndof_const}, {Ndof, Ndof_const}, {Nmeas, Nmeas}, 0 );
	  
	  Z_extr[irc].distr_list.push_back(RES);
	  
	  //print fit function
	  int Nprint=300;
	  distr_t_list Z_fit(UseJack), Z_const_fit(UseJack);
	  Vfloat musea_to_print(Nprint);
	  for(int i=0;i<Nprint;i++) musea_to_print[i] = i*1.5*musea_list[iens][musea_list[iens].size()-1]/(Nprint -1.0);
	  for(int mm=0;mm<Nprint;mm++) {
	    double mu= musea_to_print[mm];
	    distr_t Z_distr(UseJack), Z_const_distr(UseJack);
	    for(int ijack=0;ijack<(UseJack?Njacks:800);ijack++) {
	      ipar_msea IP;
	      IP.mu = mu;
	      Z_distr.distr.push_back( bf_msea.ansatz(Bt_fit_msea.par[ijack], IP));
	      Z_const_distr.distr.push_back( bf_msea.ansatz(Bt_fit_msea_const.par[ijack], IP));
	    }
	    Z_fit.distr_list.push_back( Z_distr);
	    Z_const_fit.distr_list.push_back( Z_const_distr);
	  }
	  
	  //print fit function and results
	  Print_To_File({}, {musea_to_print, Z_fit.ave(), Z_fit.err(), Z_const_fit.ave(), Z_const_fit.err()}, "../data/RC_analysis/"+RCs[irc]+"/mu_sea_extr/"+Ens+"/fit_func/p_"+to_string(imom)+".func","" , "#ch2/dof : "+to_string_with_precision(ch2_red_msea,5)+" [ R+A*mu ], "+to_string_with_precision(ch2_red_msea_const,5)+" [R], Nmeas: "+ to_string(Nmeas)) ;
	  Print_To_File({}, {musea_list[iens], Z_list_rev[iens][imom].ave(), Z_list_rev[iens][imom].err()}, "../data/RC_analysis/"+RCs[irc]+"/mu_sea_extr/"+Ens+"/p_"+to_string(imom)+".dat", "", "");
	  
	  
      }
      //print final results
      Print_To_File({ }, {Z_extr[irc].ave(), Z_extr[irc].err()}, "../data/RC_analysis/"+RCs[irc]+"/mu_sea_extr/"+Ens+"/extr.plist", "", "");
      
    }



    
    



    
  }




  return;
}


