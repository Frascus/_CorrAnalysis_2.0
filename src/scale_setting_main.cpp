#include "../include/scale_setting_main.h"



using namespace std;

//set constants

const bool UseJack=1;
const int Njacks=50; //50
const int Nboots=200;
const double fm_to_inv_Gev = 1.0 / 0.197327;
const bool Use_three_finest_in_scale_setting_fp=false;



void Get_scale_setting() {

  omp_set_num_threads(1);




   bool Get_ASCII= true;

    if(Get_ASCII) {
    //read binary files
    //boost::filesystem::create_directory("../E112_bin");
    

      vector<string> Ens_T1({"B.72.64", "B.72.96", "C.06.80","D.54.96","E.44.112"});
      vector<string> Ens_TT1({"cB211b.072.64", "cB211b.072.96", "cC211a.06.80",  "cD211a.054.96",  "cE211a.044.112"});


    boost::filesystem::create_directory("../corr_scale_setting/light");

    for( int it=0; it<(signed)Ens_T1.size(); it++) {


      
      
      boost::filesystem::create_directory("../corr_scale_setting/light/"+Ens_TT1[it]);
      //read binary
          
      cout<<"Analyzing Ens: "<<Ens_T1[it]<<endl<<flush;

      
      FILE *stream = fopen( ("../gm2_tau_rep_bin/"+Ens_T1[it]+"/"+((Ens_TT1[it]=="cE211a.044.112")?"mix_l_l_TM_P5P5":"ll_TM_P5P5")).c_str(), "rb");
      size_t Nconfs, T, Nhits;
      bin_read(Nconfs, stream);
      bin_read(Nhits, stream);
      bin_read(T, stream);

      cout<<"Nconfs: "<<Nconfs<<endl<<flush;
      cout<<"T: "<<T<<" "<<T/2+1<<endl<<flush;
      cout<<"Nhits: "<<Nhits<<endl<<flush;
      for(size_t iconf=0;iconf<Nconfs;iconf++) {
	
	string conf_id = to_string(iconf);
	if(conf_id.length() == 1) conf_id = "000"+conf_id;
	if(conf_id.length() == 2) conf_id = "00"+conf_id;
	if(conf_id.length() == 3) conf_id = "0"+conf_id;
	if(conf_id.length() > 4 ) crash("conf_id.size(): "+to_string(conf_id.size())+" conf: "+conf_id);

	conf_id = conf_id+"_r0";
	vector<double> C(T/2+1);
	for(size_t t=0;t<T/2+1;t++) bin_read(C[t], stream);
	boost::filesystem::create_directory("../corr_scale_setting/light/"+Ens_TT1[it]+"/"+conf_id);
	ofstream PrintCorr("../corr_scale_setting/light/"+Ens_TT1[it]+"/"+conf_id+"/mes_contr_2pts_ll_1");
	PrintCorr.precision(16);
	PrintCorr<<"# "<<"P5P5"<<endl;
	for(size_t t=0;t<(T/2+1);t++) PrintCorr<<C[t]<<"   "<<"0.000000"<<endl;
	for(size_t t=T/2+1; t<T;t++) PrintCorr<<C[T-t]<<"   "<<"0.000000"<<endl;
	PrintCorr.close();
	
      }
      
      fclose(stream);
      
	
	
	cout<<"Ens: "<<Ens_T1[it]<<" finished!"<<endl;
    }
    }

    D(1);



  


  //create directories g-2
  boost::filesystem::create_directory("../data/scale_setting");
  boost::filesystem::create_directory("../data/scale_setting/Mp");
  boost::filesystem::create_directory("../data/scale_setting/fp");
 

  
  //############################################INTERPOLATE PHI FUNCTION AND DERIVATIVES#############################

    


  
 
  //init Gaussian number generator
  GaussianMersenne GM(943832);


 
  
  data_t  pt2_pion, pt2_pion_B25, pt2_pion_A;
    


  //Custom sorting for V_light to account for the two replica r0 and r1, ...., rn

 
  
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


 
   
  pt2_pion.Read("../corr_scale_setting/light", "mes_contr_2pts_ll_1", "P5P5", Sort_light_confs);
  pt2_pion_B25.Read("../corr_scale_setting/B25_light_ens", "mes_contr_2pts_ll_1", "P5P5", Sort_light_confs);
  pt2_pion_A.Read("../corr_scale_setting/A_light_ens", "mes_contr_2pts_ll_1", "P5P5", Sort_light_confs);


 

  //##################################################################################
  //###################    SCALE SETTING  USING fp, Mp    ########################


 
  distr_t a_from_fp_A(UseJack), a_from_fp_B(UseJack), a_from_fp_C(UseJack), a_from_fp_D(UseJack), a_from_fp_E(UseJack);


  distr_t_list Mpi_scale_setting(UseJack), fpi_scale_setting(UseJack);
  distr_t_list Mpi_scale_setting_phys_point_ens(UseJack), fpi_scale_setting_phys_point_ens(UseJack);
  distr_t_list Mpi_scale_setting_Bens(UseJack), fpi_scale_setting_Bens(UseJack);
  distr_t_list Mpi_scale_setting_Aens(UseJack), fpi_scale_setting_Aens(UseJack);
  vector<string> Ensemble_scale_setting_tag_list;
  vector<string> Ensemble_phys_point_tag_list;
  vector<string> Ensemble_B_tag_list;
  vector<string> Ensemble_A_tag_list;
  Vfloat L_scale_setting_list;
  Vfloat L_phys_point, L_B_ens, L_A_ens;
  //physical point ensemble
  for(int iens=0;iens< pt2_pion.size; iens++) {

      //Lattice info
      LatticeInfo L_info;
      L_info.LatInfo_new_ens(pt2_pion.Tag[iens]);
      CorrAnalysis Corr(UseJack, Njacks,Nboots);
      Corr.Nt = pt2_pion.nrows[iens];
      //Analyze correlators
      if(pt2_pion.Tag[iens].substr(1,12)=="B211b.072.96") {Corr.Tmin=40; Corr.Tmax=80;}
      else if(pt2_pion.Tag[iens].substr(1,12)=="B211b.072.64") { Corr.Tmin=35; Corr.Tmax=55;}
      else if(pt2_pion.Tag[iens].substr(1,1)=="C") {Corr.Tmin=43; Corr.Tmax=65;}
      else if(pt2_pion.Tag[iens].substr(1,1)=="D") {Corr.Tmin=45; Corr.Tmax=80;}
      else if(pt2_pion.Tag[iens].substr(1,1)=="E") {Corr.Tmin=60; Corr.Tmax=90;}
      else crash("In scale setting analysis cannot find Tmin,Tmax for ensemble: "+pt2_pion.Tag[iens]);
      distr_t_list pion_corr = Corr.corr_t(pt2_pion.col(0)[iens], "");
      distr_t_list Mpi_eff_distr = Corr.effective_mass_t(pt2_pion.col(0)[iens], "../data/scale_setting/Mp/Mpi_"+pt2_pion.Tag[iens]+".dat");


      distr_t_list fpi_eff_distr = Corr.decay_constant_t(pow(2.0*L_info.ml,2)*pion_corr, "../data/scale_setting/fp/fpi_"+pt2_pion.Tag[iens]+".dat");
      distr_t Mpi_fit = Corr.Fit_distr(Mpi_eff_distr);
      distr_t fpi_fit = Corr.Fit_distr(fpi_eff_distr);
      Mpi_scale_setting.distr_list.push_back(Mpi_fit);
      fpi_scale_setting.distr_list.push_back(fpi_fit);
      Ensemble_scale_setting_tag_list.push_back(pt2_pion.Tag[iens]);
      L_scale_setting_list.push_back( L_info.L);
      Mpi_scale_setting_phys_point_ens.distr_list.push_back(Mpi_fit);
      fpi_scale_setting_phys_point_ens.distr_list.push_back(fpi_fit);
      Ensemble_phys_point_tag_list.push_back(pt2_pion.Tag[iens]);
      L_phys_point.push_back( L_info.L);
      if(pt2_pion.Tag[iens].substr(1,1) == "B") {
	Mpi_scale_setting_Bens.distr_list.push_back( Mpi_fit);
	fpi_scale_setting_Bens.distr_list.push_back( fpi_fit);
	Ensemble_B_tag_list.push_back(pt2_pion.Tag[iens]);
	L_B_ens.push_back(L_info.L);
      }
      
    }


    //B25 ensembles
    for(int iens=0;iens< pt2_pion_B25.size; iens++) {

      //Lattice info
      LatticeInfo L_info;
      L_info.LatInfo_new_ens(pt2_pion_B25.Tag[iens]);
      CorrAnalysis Corr(UseJack, Njacks,Nboots);
      Corr.Nt = pt2_pion_B25.nrows[iens];
      if(pt2_pion_B25.Tag[iens].substr(1,11)=="B211a.25.24") {Corr.Tmin=15; Corr.Tmax=20;}
      else if(pt2_pion_B25.Tag[iens].substr(1,11)=="B211a.25.32") {Corr.Tmin=23; Corr.Tmax=30;}
      else if(pt2_pion_B25.Tag[iens].substr(1,11)=="B211a.25.48") {Corr.Tmin=25; Corr.Tmax=44;}
      else if(pt2_pion_B25.Tag[iens].substr(1,11)=="B211a.14.64") {Corr.Tmin = 39; Corr.Tmax=52;}
      else crash("In scale setting analysis cannot find Tmin,Tmax for ensemble: "+pt2_pion_B25.Tag[iens]);
      distr_t_list pion_corr = Corr.corr_t(pt2_pion_B25.col(0)[iens], "");
      distr_t_list Mpi_eff_distr = Corr.effective_mass_t(pt2_pion_B25.col(0)[iens], "../data/scale_setting/Mp/Mpi_"+pt2_pion_B25.Tag[iens]+".dat");
      distr_t_list fpi_eff_distr = Corr.decay_constant_t(pow(2.0*L_info.ml,2)*pion_corr, "../data/scale_setting/fp/fpi_"+pt2_pion_B25.Tag[iens]+".dat");
      distr_t Mpi_fit = Corr.Fit_distr(Mpi_eff_distr);
      distr_t fpi_fit = Corr.Fit_distr(fpi_eff_distr);
      Mpi_scale_setting.distr_list.push_back(Mpi_fit);
      fpi_scale_setting.distr_list.push_back(fpi_fit);
      Ensemble_scale_setting_tag_list.push_back(pt2_pion_B25.Tag[iens]);
      Mpi_scale_setting_Bens.distr_list.push_back(Mpi_fit);
      fpi_scale_setting_Bens.distr_list.push_back(fpi_fit);
      Ensemble_B_tag_list.push_back( pt2_pion_B25.Tag[iens]);
      L_scale_setting_list.push_back( L_info.L);
      L_B_ens.push_back( L_info.L);
     
    }

  

    //A ensembles
    for(int iens=0;iens< pt2_pion_A.size; iens++) {

      if(pt2_pion_A.Tag[iens].substr(1,1) == "A") {
      //Lattice info
      LatticeInfo L_info;
      L_info.LatInfo_new_ens(pt2_pion_A.Tag[iens]);
      double aml_A= L_info.ml;
      CorrAnalysis Corr(UseJack, Njacks,Nboots);
      Corr.Nt = pt2_pion_A.nrows[iens];
      if(pt2_pion_A.Tag[iens] == "cA211a.40.24") {Corr.Tmin=15; Corr.Tmax=21;}
      else if(pt2_pion_A.Tag[iens] == "cA211a.53.24") {Corr.Tmin=15; Corr.Tmax=22;}
      else if(pt2_pion_A.Tag[iens] == "cA211ab.30.32") {Corr.Tmin=16; Corr.Tmax=29;}
      else if(pt2_pion_A.Tag[iens] == "cA211a.12.48") {Corr.Tmin=21; Corr.Tmax=34;}
      else crash("In scale setting analysis cannot find Tmin,Tmax for ensemble: "+pt2_pion_A.Tag[iens]);
      distr_t_list pion_corr = Corr.corr_t(pt2_pion_A.col(0)[iens], "");
      distr_t_list Mpi_eff_distr = Corr.effective_mass_t(pt2_pion_A.col(0)[iens], "../data/scale_setting/Mp/Mpi_"+pt2_pion_A.Tag[iens]+".dat");
      distr_t_list fpi_eff_distr = Corr.decay_constant_t(pow(2.0*L_info.ml,2)*pion_corr, "../data/scale_setting/fp/fpi_"+pt2_pion_A.Tag[iens]+".dat");
      distr_t Mpi_fit = Corr.Fit_distr(Mpi_eff_distr);
      distr_t fpi_fit = Corr.Fit_distr(fpi_eff_distr);
      Mpi_scale_setting.distr_list.push_back(Mpi_fit);
      fpi_scale_setting.distr_list.push_back(fpi_fit);
      Mpi_scale_setting_Aens.distr_list.push_back(Mpi_fit);
      fpi_scale_setting_Aens.distr_list.push_back(fpi_fit);
      Ensemble_A_tag_list.push_back( pt2_pion_A.Tag[iens]);
      L_scale_setting_list.push_back( L_info.L);
      L_A_ens.push_back( L_info.L);
      cout<<pt2_pion_A.Tag[iens]<<" am_ud: "<<aml_A<<" aMpi: "<<Mpi_fit.ave()<<" +- "<<Mpi_fit.err()<<endl;
      }
      if(pt2_pion_A.Tag[iens].substr(1,1) == "A") {
	Ensemble_scale_setting_tag_list.push_back(pt2_pion_A.Tag[iens]);
      }

     
    }
  
 
  
  
    Determine_scale_from_fp_FLAG(Mpi_scale_setting_phys_point_ens, fpi_scale_setting_phys_point_ens, Mpi_scale_setting_Bens, fpi_scale_setting_Bens, Mpi_scale_setting_Aens, fpi_scale_setting_Aens, L_phys_point, L_B_ens, L_A_ens, Ensemble_phys_point_tag_list, Ensemble_B_tag_list, Ensemble_A_tag_list, a_from_fp_A, a_from_fp_B, a_from_fp_C, a_from_fp_D, a_from_fp_E,  UseJack, Use_three_finest_in_scale_setting_fp);


    //print pion masses
    cout<<"###############  PRINTING PION MASSES: "<<endl;
    for(int i=0;i<(signed)Ensemble_phys_point_tag_list.size();i++) {

      string T= Ensemble_phys_point_tag_list[i];
      distr_t M= Mpi_scale_setting_phys_point_ens.distr_list[i];
      distr_t a(UseJack);
      if(T.substr(1,1) == "B") { a= a_from_fp_B;}
      else if(T.substr(1,1) == "C" ) { a= a_from_fp_C;}
      else if(T.substr(1,1) == "D" ) { a= a_from_fp_D;}
      else if(T.substr(1,1) == "E" ) { a= a_from_fp_E;}
      cout<<"Mpi("<<T<<") : "<<(M/a).ave()<<" +- "<<(M/a).err()<<" [lu]: "<<M.ave()<<" +- "<<M.err()<<endl;
      
    }



    
  
  
  //###################################################################################


  return ;

}
